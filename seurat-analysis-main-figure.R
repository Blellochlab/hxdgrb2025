#Load libraries
library(Seurat)
library(dplyr)
library(gridExtra)
library(Azimuth)
library(ggplot2)
library(scales)
library(patchwork)
library(ComplexHeatmap)

#Set working dir
setwd("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/")

#Read qc'd experiments
hxrb03 <- readRDS(file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxrb03-qc-cells.rds"))
hxrb04 <- readRDS(file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxrb04-qc-cells.rds"))

#Filter conditons and merge experiments
hxrb03[["batch"]] <- "hxrb03"
hxrb04[["batch"]] <- "hxrb04"
dat <- merge(hxrb03,hxrb04)
exclude <- grepl("dKO",dat[[]]$sample)
dat <- dat[,!exclude]
#ind <- grepl("Flt3|G-DC|-LT",dat[[]]$condition)
ind <- grepl("Flt3|-LT|G-DC-In",dat[[]]$condition)
dat <- dat[,!ind]
# dat[[]]$condition <- factor(dat[[]]$condition,
#                             levels = c("iPS","YS-Mye","YS-G-IL4-DC","YS-G-IL4-DC-In"))
unique(dat[[]]$condition)
dat[[]]$condition <- factor(dat[[]]$condition,
                            levels = c("iPS","YS-Mye",
                                       "YS-G-DC","YS-G-IL4-DC","YS-G-IL4-DC-In"))
rm(hxrb03,hxrb04)
gc()

#Non-integrated preprocessing and UMAP
all.genes <- rownames(dat)
dat <- NormalizeData(dat,verbose=F) %>% FindVariableFeatures(verbose=F) %>% ScaleData(verbose=F,features=all.genes) %>% RunPCA(verbose=F)

dat <- RunUMAP(dat,dims=1:30,reduction='pca',reduction.name = "umap.unintegrated")

#Integrate data
dat <- IntegrateLayers(object = dat, method = CCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.cca",verbose = FALSE) #takes ~1min
dat[["RNA"]] <- JoinLayers(dat[["RNA"]])
dat <- FindNeighbors(dat, reduction = "integrated.cca", dims = 1:30)
dat <- FindClusters(dat, resolution = 1)
dat <- RunUMAP(dat, dims = 1:30, reduction = "integrated.cca")

#Save integrated data
saveRDS(dat,file.path("main-integrated.rds"))
rm(dat)
gc()

#Plot integration results as UMAP
dat <- readRDS("main-integrated.rds")
set.seed(345)
cols <- c("darkgrey","black")
p1 <- DimPlot(dat, reduction = "umap.unintegrated", shuffle=T, group.by = c("batch"),cols=cols,pt.size=1)
p2 <- DimPlot(dat, reduction = "umap", shuffle=T, group.by = c("batch"),cols=cols,pt.size=1)
cols <- DiscretePalette(18, shuffle = TRUE)
p3 <- DimPlot(dat, reduction = "umap", shuffle=F, group.by = c("seurat_clusters"),cols=cols,pt.size=1)
cols <- colorRampPalette(c("orange", "dodgerblue"))(4)
cols <- c("grey","orange","gold","dodgerblue","darkblue")
p4 <- DimPlot(dat, reduction = "umap", shuffle=T, group.by = c("condition"),cols=cols,pt.size = 1)

p1 + p2
ggsave("plots/main_umap_batchremoval.pdf",width=200,height=75, units="mm") #factor 4 vs what is shown in save image
p3
ggsave("plots/main_umap_seuratclusters.pdf",width=90,height=75, units="mm") #factor 4 vs what is shown in save image
p4
ggsave("plots/main_umap_conditions.pdf",width=120,height=75, units="mm") #factor 4 vs what is shown in save image

#Plot all clusters
plot_list <- list()
for (i in factor(0:17)) {
  CellsToPlot <- colnames(dat)[dat$seurat_clusters == i]
  p <- DimPlot(dat, cells.highlight = CellsToPlot, reduction="umap", cols.highlight = "dodgerblue",
               cols = "gray", sizes.highlight = 0.2, order = TRUE, shuffle = TRUE, label=F) + ggtitle(i) + NoLegend()  
  plot_list[[i]] <- p
}
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 5))
rm(dat)
gc()

#Find cluster markers
main <- readRDS(file.path("main-integrated.rds"))
allmarkers <- FindAllMarkers(main, only.pos = TRUE)
allmarkers <- allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01)
write.csv(allmarkers,file.path("main-all-sign-sr-cluster-markers.csv"))

#Find condition markers
main <- readRDS(file.path("main-integrated.rds"))
main <- SetIdent(main,value="condition")
allmarkers <- FindAllMarkers(main, only.pos = TRUE)
allmarkers <- allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01)
write.csv(allmarkers,file.path("main-all-sign-sr-cluster-markers.csv"))

#Cell type inference (Azimuth)
# main <- readRDS(file.path("main-integrated.rds"))
# main.ann <- RunAzimuth(main,reference="pbmcref")
#Save integrated data
# saveRDS(main.ann,file.path("main-integrated-azimuth-mapped.rds"))
rm(main)
rm(main.ann)
gc()

main.ann <- readRDS(file.path("main-integrated-azimuth-mapped.rds"))

cols <- DiscretePalette(5, shuffle = TRUE)
p1 <- DimPlot(main.ann, reduction = "umap", group.by = "predicted.celltype.l1", label = F, shuffle=F, order=T, cols=cols, label.size = 1)
cols <- DiscretePalette(10, shuffle = TRUE)
#p2 <- DimPlot(main.ann, reduction = "umap", group.by = "predicted.celltype.l2", label = F, shuffle=F, order=T, cols=cols, label.size = 1)
#cols <- colorRampPalette(c("orange", "dodgerblue"))(4)
cols <- c("grey","orange","gold","dodgerblue","darkblue")
p2 <- DimPlot(main.ann, reduction = "umap", group.by = "condition", label = F, cols=cols, label.size = 1)
p3 <- FeaturePlot(main.ann, reduction = "umap", features = "predictionscorecelltypel1_DC")
p4 <- FeaturePlot(main.ann, reduction = "umap", features = "predictionscorecelltypel1_Mono")

p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
ggsave("plots/main_umap_azimuth_mapped.pdf",width=225,height=150, units="mm") #factor 4 vs what is shown in save image

#Violin plot representation of DC score
VlnPlot(main.ann, features = "predictionscorecelltypel1_DC", group.by = "condition",pt.size=F, cols=cols) + ylab("DC prediction score")
ggsave("plots/main_violin_azimuth_DCscore.pdf",width=100,height=75, units="mm") #factor 4 vs what is shown in save image

rm(main.ann)
gc()

#Differentiation markers
main <- readRDS(file.path("main-integrated.rds"))

markers <- list(
  "iPS" = c(
    "POU5F1",
    "SOX2",
    "LIN28A",
    "NANOG",
    "TDGF1",
    "TERF1",
    "DNMT3B",
    "EPCAM",
    "PODXL",
    "TRIM71"
  ),
  "YS-MYE" = c(
    "CD14",
    "FCGR3A",
    "ITGAM",
    "CD33",
    "FCGR1A"
  ),
  "YS-G-IL4" = c(
    "CD68",
    "CSF1R",
    "CD209",
    "CLEC4A",
    "ANPEP",
    "TREM2",
    "MRC1",
    "FCER1G",
    "S100A9"
  ),
  "YS-G-IL4-In" = c(
    "FLT3",
    "CD86",
    "CCR7",
    "LAMP3",
    "CD83",
    "IRF4",
    "BATF3",
    "CD80"
  )
)

allgenes <- c(markers$iPS,markers$`YS-MYE`,markers$`YS-G-IL4`,markers$`YS-G-IL4-In`)

#Heatmap for these genes
scaled_data <- main@assays$RNA$scale.data
ind <- rownames(scaled_data) %in% allgenes
scaled_data <- scaled_data[ind,]

# Calculate Z-scores
zscore_data <- t(scale(t(scaled_data)))  # Transpose, scale, then transpose back

# Create a heatmap with hierarchical clustering
condition_data <- main@meta.data$condition 

col_fun <- circlize::colorRamp2(c(min(zscore_data),max(zscore_data)), c("white","blue"))

condition_annotation <- HeatmapAnnotation(
  condition = condition_data,
  col = list(condition = c("iPS" = "grey", "YS-Mye" = "orange", "YS-G-DC" = "gold", "YS-G-IL4-DC" = "dodgerblue", "YS-G-IL4-DC-In" = "darkblue"))
)

Heatmap(zscore_data,
        cluster_rows = TRUE,   # Cluster genes
        cluster_columns = TRUE,  # Don't cluster conditions
        show_row_names = TRUE,  # Show gene names
        show_column_names = FALSE,  # Show condition names
        column_split= condition_data,
        cluster_column_slices = FALSE, 
        name = "Z-score",
        col=col_fun,
        top_annotation = condition_annotation,
        row_names_gp = gpar(fontsize = 5),
        heatmap_legend_param = list(title = "Z-score")
)

##Dotplot by marker genes by culture condition
Idents(main) <- main$condition
DotPlot(object = main, features=markers, cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_diffmarkers.pdf",width=300,height=100, units="mm") #factor 4 vs what is shown in save image
#Avg Expr.
DotPlot(object = main, features=markers, cluster.idents=T,scale=F) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_diffmarkers_avgexpr.pdf",width=300,height=100, units="mm") #factor 4 vs what is shown in save image

rm(list = ls()) 
gc()

#############################################################################################################################
######################### YS-G-IL4-DC-In Subclustering ######################################################################
#############################################################################################################################
#Read qc'd experiments
hxrb03 <- readRDS(file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxrb03-qc-cells.rds"))
hxrb04 <- readRDS(file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxrb04-qc-cells.rds"))

#Filter conditons and merge experiments
hxrb03[["batch"]] <- "hxrb03"
hxrb04[["batch"]] <- "hxrb04"
dat <- merge(hxrb03,hxrb04)
ind <- grepl("YS-G-IL4-DC-In",dat[[]]$condition)
dat <- dat[,ind]

#Non-integrated preprocessing and UMAP
all.genes <- rownames(dat)
dat <- NormalizeData(dat,verbose=F) %>% FindVariableFeatures(verbose=F) %>% ScaleData(verbose=F,features=all.genes) %>% RunPCA(verbose=F)
dat <- RunUMAP(dat,dims=1:30,reduction='pca',reduction.name = "umap.unintegrated")

#Integrate data
dat <- IntegrateLayers(object = dat, method = CCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.cca",verbose = FALSE) #takes ~1min
dat[["RNA"]] <- JoinLayers(dat[["RNA"]])
dat <- FindNeighbors(dat, reduction = "integrated.cca", dims = 1:30)
dat <- FindClusters(dat, resolution = 1)
dat <- RunUMAP(dat, dims = 1:30, reduction = "integrated.cca")

#Save integrated data
saveRDS(dat,file.path("main-IL4-DC-In-integrated.rds"))
rm(dat)
gc()

#Plot integration results as UMAP
dat <- readRDS("main-IL4-DC-In-integrated.rds")
set.seed(345)
cols <- c("darkgrey","black")
p1 <- DimPlot(dat, reduction = "umap.unintegrated", shuffle=T, group.by = c("batch"),cols=cols,pt.size=1)
p2 <- DimPlot(dat, reduction = "umap", shuffle=T, group.by = c("batch"),cols=cols,pt.size=1)
cols <- DiscretePalette(10, shuffle = TRUE)
p3 <- DimPlot(dat, reduction = "umap", shuffle=F, group.by = c("seurat_clusters"),cols=cols,pt.size=1)

p1 + p2
ggsave("plots/main_umap_IL4-DC-In_batchremoval.pdf",width=200,height=75, units="mm") #factor 4 vs what is shown in save image
p3
ggsave("plots/main_umap_IL4-DC-In_seuratclusters.pdf",width=90,height=75, units="mm") #factor 4 vs what is shown in save image

rm(dat)
gc()

#Find cluster markers
main <- readRDS("main-IL4-DC-In-integrated.rds")
allmarkers <- FindAllMarkers(main, only.pos = TRUE)
allmarkers <- allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01)
write.csv(allmarkers,file.path("main-IL4-DC-In-sign-sr-cluster-markers.csv"))
rm(main)
gc()

#Differentiation markers
main <- readRDS("main-IL4-DC-In-integrated.rds")

markers <- list(
  "YS-MYE" = c(
    "CD14",
    "FCGR3A",
    "ITGAM",
    "CD33",
    "FCGR1A"
  ),
  "YS-G-IL4" = c(
    "CD68",
    "CSF1R",
    "CD209",
    "CLEC4A",
    "ANPEP",
    "TREM2",
    "MRC1",
    "FCER1G",
    "S100A9"
  ),
  "YS-G-IL4-In" = c(
    "FLT3",
    "CD86",
    "CCR7",
    "LAMP3",
    "CD83",
    "IRF4",
    "BATF3",
    "CD80"
  )
)

allgenes <- c(markers$iPS,markers$`YS-MYE`,markers$`YS-G-IL4`,markers$`YS-G-IL4-In`)

#Feature, mt, counts per cluster
VlnPlot(main, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=F)
ggsave("plots/main_violin_IL4-DC-In_srcluster_qcmetrics.pdf",width=150,height=75, units="mm") #factor 4 vs what is shown in save image

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers, cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_IL4-DC-In_srcluster.pdf",width=200,height=100, units="mm") #factor 4 vs what is shown in save image

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers, cluster.idents=T,scale=F) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_IL4-DC-In_srcluster_avgexpr.pdf",width=200,height=100, units="mm") #factor 4 vs what is shown in save image

#Identify cells that are part of low information clusters
ind <- main[[]]$seurat_clusters %in% c(5,7,8,9)
saveRDS(colnames(main)[ind],file.path("main-IL4-DC-In-integrated-cluster5789-cellbarcodes.rds"))
rm(main)
gc()

#############################################################################################################################
######################### YS-G-IL4-DC-In Subclustering (with low information cluster removal) ######################################################################
#############################################################################################################################
#Read qc'd experiments
hxrb03 <- readRDS(file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxrb03-qc-cells.rds"))
hxrb04 <- readRDS(file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxrb04-qc-cells.rds"))

#Filter conditons and merge experiments
hxrb03[["batch"]] <- "hxrb03"
hxrb04[["batch"]] <- "hxrb04"
dat <- merge(hxrb03,hxrb04)
ind <- grepl("YS-G-IL4-DC-In",dat[[]]$condition)
dat <- dat[,ind]

#Filter out low quality cells from clusters 4 and 5
ind <- readRDS(file.path("main-IL4-DC-In-integrated-cluster5789-cellbarcodes.rds"))
ind <- colnames(dat) %in% ind
dat <- dat[,!ind]

rm(hxrb03,hxrb04)
gc()

#Non-integrated preprocessing and UMAP
all.genes <- rownames(dat)
dat <- NormalizeData(dat,verbose=F) %>% FindVariableFeatures(verbose=F) %>% ScaleData(verbose=F,features=all.genes) %>% RunPCA(verbose=F)
dat <- RunUMAP(dat,dims=1:30,reduction='pca',reduction.name = "umap.unintegrated")

#Integrate data
dat <- IntegrateLayers(object = dat, method = CCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.cca",verbose = FALSE) #takes ~1min
dat[["RNA"]] <- JoinLayers(dat[["RNA"]])
dat <- FindNeighbors(dat, reduction = "integrated.cca", dims = 1:30)
dat <- FindClusters(dat, resolution = 0.5)

dat <- RunUMAP(dat, dims = 1:30, reduction = "integrated.cca")

#Save integrated data
saveRDS(dat,file.path("main-IL4-DC-In-integrated-woclusters.rds"))
rm(dat)
gc()

#Plot integration results as UMAP
dat <- readRDS("main-IL4-DC-In-integrated-woclusters.rds")
set.seed(345)
cols <- c("darkgrey","black")
p1 <- DimPlot(dat, reduction = "umap.unintegrated", shuffle=T, group.by = c("batch"),cols=cols,pt.size=1)
p2 <- DimPlot(dat, reduction = "umap", shuffle=T, group.by = c("batch"),cols=cols,pt.size=1)
cols <- DiscretePalette(4, shuffle = TRUE)
p3 <- DimPlot(dat, reduction = "umap", shuffle=F, group.by = c("seurat_clusters"),cols=cols,pt.size=1)
cols <- DiscretePalette(2, shuffle = TRUE)
p4 <- DimPlot(dat, reduction = "umap", shuffle=T, group.by = c("sample"),cols=cols,pt.size=1)

p1 + p2
ggsave("plots/main_umap_IL4-DC-In_batchremoval_woclusters.pdf",width=200,height=75, units="mm") #factor 4 vs what is shown in save image
p3
ggsave("plots/main_umap_IL4-DC-In_seuratclusters_woclusters.pdf",width=90,height=75, units="mm") #factor 4 vs what is shown in save image
p4
ggsave("plots/main_umap_IL4-DC-In_genotype_woclusters.pdf",width=140,height=75, units="mm") #factor 4 vs what is shown in save image
rm(dat)
gc()

#Find cluster markers
main <- readRDS("main-IL4-DC-In-integrated-woclusters.rds")
allmarkers <- FindAllMarkers(main, only.pos = TRUE)
allmarkers <- allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01)
write.csv(allmarkers,file.path("main-IL4-DC-In-sign-sr-cluster-markers-woclusters.csv"))
rm(main)
gc()

# #Cell type inference (Azimuth)
# main <- readRDS("main-IL4-DC-In-integrated-woclusters.rds")
# main.ann <- RunAzimuth(main,reference="pbmcref")
#Save integrated data
# saveRDS(main.ann,file.path("main-IL4-DC-In-integrated-woclusters-azimuth-mapped.rds"))

#Define new clusters based on predicted cell type
main.ann <- readRDS(file.path("main-IL4-DC-In-integrated-woclusters-azimuth-mapped.rds"))

# Access the prediction scores for cell type level 1
pred_scores_l1 <- GetAssayData(object = main.ann, assay = "prediction.score.celltype.l1", layer = "data")

# Extract DC and Mono scores
if ("DC" %in% rownames(pred_scores_l1)) {
  dc_scores <- pred_scores_l1["DC", ]
} else {
  stop("Failed to find scores.")
}

if ("Mono" %in% rownames(pred_scores_l1)) {
  mono_scores <- pred_scores_l1["Mono", ]
} else {
  stop("Failed to find scores.")
}

plot(dc_scores,mono_scores)
dens <- MASS::kde2d(dc_scores,mono_scores, n = 100)
# Plot as image
image(dens, col = topo.colors(20), xlab = "DC score",ylab="Mono score")
# Add contours
contour(dens, add = TRUE)

main.ann$DC_scores <- dc_scores
main.ann$Mono_scores <- mono_scores
main.ann[[]]$annotated_clusters <- ifelse(main.ann[[]]$Mono_scores > 0.6 & main.ann[[]]$DC_scores < 0.15, "Mono-like",
                                          ifelse(main.ann[[]]$Mono_scores < 0.3 & main.ann[[]]$DC_scores > 0,"DC-like","Other"))

#main[[]]$annotated_clusters <- ifelse(main[[]]$seurat_clusters %in% c(2,5,6), "mono-like","dc-like")
saveRDS(main.ann,file.path("main-mat3-wgeno-woclusters-reannotated-azimuth-mapped.rds"))
rm(main.ann)
gc()


#Calculate fractions of cells in reannotated clusters
main <- readRDS("main-mat3-wgeno-woclusters-reannotated-azimuth-mapped.rds")
group_col <- "genotype"
stack_col <- "annotated_clusters"
metadata_df <- main@meta.data
p <- ggplot(metadata_df, aes_string(x = group_col, fill = stack_col)) +
  geom_bar(position = "fill") +                   
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c("black","darkgrey","lightgray")) +
  labs(
    x = group_col,
    y = "% cells",
    fill = stack_col
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )
p
rm(main)
gc()

p1 <- FeaturePlot(main.ann, reduction = "umap", order=T, features = "predictionscorecelltypel1_Mono")
p2 <- FeaturePlot(main.ann, reduction = "umap", order=T, features = "predictionscorecelltypel1_DC")
p1 + p2 + plot_layout(ncol = 2)
ggsave("plots/main_umap_IL4-DC-In_woclusters_azimuth_mapped.pdf",width=300,height=100, units="mm") #factor 4 vs what is shown in save image

#Violin plot representation of DC score
p1 <- VlnPlot(main.ann, features = "predictionscorecelltypel1_DC", group.by = "seurat_clusters",pt.size=F, cols=cols) + ylab("DC prediction score")
p2 <- VlnPlot(main.ann, features = "predictionscorecelltypel1_Mono", group.by = "seurat_clusters",pt.size=F, cols=cols) + ylab("Mono prediction score")
p1 + p2 + plot_layout(ncol=2)
ggsave("plots/main_violin_IL4-DC-In_woclusters_azimuth_mapped.pdf",width=100,height=75, units="mm") #factor 4 vs what is shown in save image

#Differentiation markers
main <- readRDS("main-IL4-DC-In-integrated-woclusters.rds")

markers <- list(
  "YS-MYE" = c(
    "CD14",
    "FCGR3A",
    "ITGAM",
    "CD33",
    "FCGR1A"
  ),
  "YS-G-IL4" = c(
    "CD68",
    "CSF1R",
    "CD209",
    "CLEC4A",
    "ANPEP",
    "TREM2",
    "MRC1",
    "FCER1G",
    "S100A9"
  ),
  "YS-G-IL4-In" = c(
    "FLT3",
    "CD86",
    "CCR7",
    "LAMP3",
    "CD83",
    "IRF4",
    "BATF3",
    "CD80"
  )
)

allgenes <- c(markers$iPS,markers$`YS-MYE`,markers$`YS-G-IL4`,markers$`YS-G-IL4-In`)

#Feature, mt, counts per cluster
VlnPlot(main, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=F)
ggsave("plots/main_violin_IL4-DC-In_woclusters_srcluster_qcmetrics.pdf",width=150,height=75, units="mm") #factor 4 vs what is shown in save image

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers, cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_IL4-DC-In_woclusters_srcluster.pdf",width=200,height=100, units="mm") #factor 4 vs what is shown in save image

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers, cluster.idents=T,scale=F) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_IL4-DC-In_woclusters_srcluster_avgexpr.pdf",width=200,height=100, units="mm") #factor 4 vs what is shown in save image

#Heatmap for these genes
scaled_data <- main@assays$RNA$scale.data
ind <- rownames(scaled_data) %in% allgenes
scaled_data <- scaled_data[ind,]

# Calculate Z-scores
zscore_data <- t(scale(t(scaled_data)))  # Transpose, scale, then transpose back

# Create a heatmap with hierarchical clustering
set.seed(214)
condition_data <- Idents(main) 
condition_annotation <- HeatmapAnnotation(
  condition = condition_data
)

col_fun <- circlize::colorRamp2(c(min(zscore_data),max(zscore_data)), c("white","blue"))

Heatmap(zscore_data,
        cluster_rows = TRUE,   # Cluster genes
        cluster_columns = TRUE,  # Don't cluster conditions
        show_row_names = TRUE,  # Show gene names
        show_column_names = FALSE,  # Show condition names
        column_split=condition_data,
        name = "Z-score",
        col=col_fun,
        top_annotation = condition_annotation,
        row_names_gp = gpar(fontsize = 5),
        heatmap_legend_param = list(title = "Z-score")
)

#Specific gene list
markers <- list(
  DCs= c("FLT3"),
  Monocytes = c("CSF1R","CD14", "FCGR3A", "ITGAM", "CD33", "CD68", "FCGR1A"),
  cDC1s = c("THBD", "CLEC9A", "XCR1", "CADM1", "IRF8", "BATF3", "TLR3"),
  cDC2s = c("CD1C", "ITGAX","SIRPA", "FCER1A"),
  Maturation = c("RELB", "CD83"),
  Migration = c("CCR7", "MYO1G", "FSCN1", "MARCKS", "MARCKSL1"),
  Co_stim = c("CD80", "CD86", "CD40", "ICOSLG", "TNFSF4", "TNFSF9"),
  Co_inhib = c("CD274", "PDCD1LG2", "CD276", "VTCN1", "VSIR", "LGALS9"),
  Cytokines = c("IL12A", "IL12B", "IL6", "IFNA1", "IFNB1", "TNF", "IL23B", "IL27", "EBI3", "IL15","IFNG","IL10","TFGB1"),
  Chemokines = c("CCL3", "CCL4", "CCL5", "CCL17", "CCL22", "CXCL9", "CXCL10"),
  mRegDC = c("IL4RA", "IL4I1", "BCL2", "BIRC3", "LAMP3", "IDO1")
)

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers, cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_wocluster_immune_markers_zscore.pdf",width=400,height=100, units="mm") #factor 4 vs what is shown in save image
DotPlot(object = main, features=markers, cluster.idents=T,scale=F) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_woclusters_immune_markers_avgexpr.pdf",width=400,height=100, units="mm") #factor 4 vs what is shown in save image

##Feature plots
p1 <- FeaturePlot(object=main,reduction = "umap",features="FLT3")
p2 <- FeaturePlot(object=main,reduction = "umap",features="CSF1R")
p1 + p2 + plot_layout(ncol = 2)
ggsave("plots/main_umap_IL4-DC-In_woclusters_FLT3_CSF1R_featureplots.pdf",width=300,height=100, units="mm") #factor 4 vs what is shown in save image

#Violin plot representation of DC score
p1 <- VlnPlot(main, features = "FLT3", group.by = "seurat_clusters",pt.size=F, cols=cols) + ylab("FLT3 expression")
p2 <- VlnPlot(main, features = "CSF1R", group.by = "seurat_clusters",pt.size=F, cols=cols) + ylab("CSF1R expression")
p1 + p2 + plot_layout(ncol=2)
ggsave("plots/main_violin_IL4-DC-In_woclusters_FLT3_CSF1R_expr.pdf",width=100,height=75, units="mm") #factor 4 vs what is shown in save image

#############################################################################################################################
######################### Main figure mat2 mat3 full data ######################################################################
#############################################################################################################################
#Read qc'd experiments
hxdgrb05 <- readRDS(file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxdgrb05-qc-cells.rds"))

#Filter conditons and merge experiments
hxdgrb05[["batch"]] <- "hxdgrb05"
dat <- hxdgrb05
# exclude <- grepl("qKO",dat[[]]$genotype)
# dat <- dat[,!exclude]

rm(hxdgrb05)
gc()

#Non-integrated preprocessing and UMAP
all.genes <- rownames(dat)
dat <- NormalizeData(dat,verbose=F) %>% FindVariableFeatures(verbose=F) %>% ScaleData(verbose=F,features=all.genes) %>% RunPCA(verbose=F)
dat <- RunUMAP(dat,dims=1:30,reduction='pca')
dat <- FindNeighbors(dat, reduction = "pca", dims = 1:30)
dat <- FindClusters(dat)

#Save integrated data
saveRDS(dat,file.path("main-mat2vmat3.rds"))
rm(dat)
gc()

#Plot integration results as UMAP
dat <- readRDS("main-mat2vmat3.rds")
set.seed(1)
cols <- DiscretePalette(13, shuffle = F)
p1 <- DimPlot(dat, reduction = "umap", shuffle=F, group.by = c("seurat_clusters"),cols=cols,pt.size=1)
cols <- DiscretePalette(2, shuffle = F)
p2 <- DimPlot(dat, reduction = "umap", shuffle=T, group.by = c("condition"),cols=cols,pt.size=1)

p1
ggsave("plots/main_umap_mat2vmat3_seuratclusters.pdf",width=90,height=75, units="mm") #factor 4 vs what is shown in save image
p2
ggsave("plots/main_umap_mat2vmat3_conditions.pdf",width=90,height=75, units="mm") #factor 4 vs what is shown in save image
rm(dat)
gc()

#Differentiation markers
main <- readRDS("main-mat2vmat3.rds")

markers <- list(
  "YS-MYE" = c(
    "CD14",
    "FCGR3A",
    "ITGAM",
    "CD33",
    "FCGR1A"
  ),
  "YS-G-IL4" = c(
    "CD68",
    "CSF1R",
    "CD209",
    "CLEC4A",
    "ANPEP",
    "TREM2",
    "MRC1",
    "FCER1G",
    "S100A9"
  ),
  "YS-G-IL4-In" = c(
    "FLT3",
    "CD86",
    "CCR7",
    "LAMP3",
    "CD83",
    "IRF4",
    "BATF3",
    "CD80"
  )
)

allgenes <- c(markers$iPS,markers$`YS-MYE`,markers$`YS-G-IL4`,markers$`YS-G-IL4-In`)

#Feature, mt, counts per cluster
set.seed(1)
cols <- DiscretePalette(13, shuffle = F)
VlnPlot(main, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=F,cols=cols)
ggsave("plots/main_violin_mat2vsmat3_qcmetrics.pdf",width=150,height=75, units="mm") #factor 4 vs what is shown in save image

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers, cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_mat2vsmat3_srcluster.pdf",width=200,height=100, units="mm") #factor 4 vs what is shown in save image

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers, cluster.idents=T,scale=F) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_mat2vsmat3_srcluster_avgexpr.pdf",width=200,height=100, units="mm") #factor 4 vs what is shown in save image

##Define low quality cell clusters
ind <- main[[]]$seurat_clusters %in% c(4,5,8:12)
saveRDS(colnames(main)[ind],file.path("main-mat2vsmat3-lowQcluster-cellbarcodes.rds"))
rm(main)
gc()

#############################################################################################################################
######################### Mat2 vs Mat3 (with low information cluster removal) ######################################################################
#############################################################################################################################
#Read qc'd experiments
hxdgrb05 <- readRDS(file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxdgrb05-qc-cells.rds"))

#Filter conditons and merge experiments
hxdgrb05[["batch"]] <- "hxdgrb05"
dat <- hxdgrb05
exclude <- grepl("qKO",dat[[]]$genotype)
dat <- dat[,!exclude]

#Filter out low quality cells from clusters 4 and 5
ind <- readRDS(file.path("main-mat2vsmat3-lowQcluster-cellbarcodes.rds"))
ind <- colnames(dat) %in% ind
dat <- dat[,!ind]

rm(hxdgrb05)
gc()

#Preprocessing and UMAP
all.genes <- rownames(dat)
dat <- NormalizeData(dat,verbose=F) %>% FindVariableFeatures(verbose=F) %>% ScaleData(verbose=F,features=all.genes) %>% RunPCA(verbose=F)
dat <- RunUMAP(dat,dims=1:30,reduction='pca')
dat <- FindNeighbors(dat, reduction = "pca", dims = 1:30)
dat <- FindClusters(dat)

#Save integrated data
saveRDS(dat,file.path("main-mat2vmat3-woclusters.rds"))
rm(dat)
gc()

#Plot integration results as UMAP
dat <- readRDS("main-mat2vmat3-woclusters.rds")
set.seed(1)
cols <- c("darkgrey","black")
p1 <- DimPlot(dat,reduction = "umap", shuffle=T,cols=cols,group.by=c("genotype"),split.by = "condition")
cols <- DiscretePalette(7, shuffle = FALSE)
p2 <- DimPlot(dat, reduction = "umap", shuffle=F, group.by = c("seurat_clusters"),cols=cols,pt.size=1)
cols <- DiscretePalette(2, shuffle = F)
p3 <- DimPlot(dat, reduction = "umap", shuffle=T, group.by = c("condition"),cols=cols,pt.size=1)

p1
ggsave("plots/main_umap_mat2vmat3_wtvsdko_woclusters.pdf",width=200,height=75, units="mm") #factor 4 vs what is shown in save image
p2
ggsave("plots/main_umap_mat2vmat3_seuratclusters_woclusters.pdf",width=90,height=75, units="mm") #factor 4 vs what is shown in save image
p3
ggsave("plots/main_umap_mat2vmat3_conditions_woclusters.pdf",width=140,height=75, units="mm") #factor 4 vs what is shown in save image
rm(dat)
gc()

#Find cluster markers
main <- readRDS("main-mat2vmat3-woclusters.rds")
allmarkers <- FindAllMarkers(main, only.pos = TRUE)
allmarkers <- allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01)
write.csv(allmarkers,file.path("main-mat2vmat3-woclusters-cluster-markers.csv"))
rm(main)
gc()

# #Cell type inference (Azimuth)
# main <- readRDS("main-mat2vmat3-woclusters.rds")
# main.ann <- RunAzimuth(main,reference="pbmcref")
# #Save integrated data
# saveRDS(main.ann,file.path("main-mat2vmat3-woclusters-azimuth-mapped.rds"))
main.ann <- readRDS(file.path("main-mat2vmat3-woclusters-azimuth-mapped.rds"))

p1 <- FeaturePlot(main.ann, reduction = "umap", order=T, features = "predictionscorecelltypel1_Mono")
p2 <- FeaturePlot(main.ann, reduction = "umap", order=T, features = "predictionscorecelltypel1_DC")
p1 + p2 + plot_layout(ncol = 2)
ggsave("plots/main_umap_mat2vmat3-woclusters_azimuth_mapped.pdf",width=300,height=100, units="mm") #factor 4 vs what is shown in save image

#Violin plot representation of DC score
main.ann <- readRDS(file.path("main-mat2vmat3-woclusters-azimuth-mapped.rds"))
set.seed(1)
cols <- DiscretePalette(7, shuffle = F)
p1 <- VlnPlot(main.ann, features = "predictionscorecelltypel1_DC", pt.size=F, cols=cols) + ylab("DC prediction score")
p2 <- VlnPlot(main.ann, features = "predictionscorecelltypel1_Mono", pt.size=F, cols=cols) + ylab("Mono prediction score")
p1 + p2 + plot_layout(ncol=2)
ggsave("plots/main_violin_umap_mat2vmat3-woclusters_azimuth_mapped.pdf",width=100,height=75, units="mm") #factor 4 vs what is shown in save image
rm(main.ann)
gc()

#Define new clusters based on predicted cell type
main.ann <- readRDS("main-mat2vmat3-woclusters-azimuth-mapped.rds")

# Access the prediction scores for cell type level 1
pred_scores_l1 <- GetAssayData(object = main.ann, assay = "prediction.score.celltype.l1", layer = "data")

# Extract DC and Mono scores
if ("DC" %in% rownames(pred_scores_l1)) {
  dc_scores <- pred_scores_l1["DC", ]
} else {
  stop("Failed to find scores.")
}

if ("Mono" %in% rownames(pred_scores_l1)) {
  mono_scores <- pred_scores_l1["Mono", ]
} else {
  stop("Failed to find scores.")
}

plot(dc_scores,mono_scores)
dens <- MASS::kde2d(dc_scores,mono_scores, n = 100)
# Plot as image
image(dens, col = topo.colors(20), xlab = "DC score",ylab="Mono score")
# Add contours
contour(dens, add = TRUE)

main.ann$DC_scores <- dc_scores
main.ann$Mono_scores <- mono_scores
main.ann[[]]$annotated_clusters <- ifelse(main.ann[[]]$Mono_scores > 0.5 & main.ann[[]]$DC_scores < 0.2, "Mono-like",
                                          ifelse(main.ann[[]]$Mono_scores < 0.3 & main.ann[[]]$DC_scores > 0,"DC-like","Other"))

#main[[]]$annotated_clusters <- ifelse(main[[]]$seurat_clusters %in% c(2,5,6), "mono-like","dc-like")
# saveRDS(main.ann,file.path("main-mat2vmat3-woclusters-reannotated-azimuth-mapped.rds"))
rm(main.ann)
gc()

#Plot DC/Mono scores by condition
main.ann <- readRDS("main-mat2vmat3-woclusters-azimuth-mapped.rds")

# Access the prediction scores for cell type level 1
pred_scores_l1_mat2 <- GetAssayData(object = main.ann[,main.ann[[]]$condition=="YS-G-IL4-DC-In"], assay = "prediction.score.celltype.l1", layer = "data")
pred_scores_l1_mat3 <- GetAssayData(object = main.ann[,main.ann[[]]$condition=="YS-G-IL4-DC-Mat3"], assay = "prediction.score.celltype.l1", layer = "data")
dc_scores_mat2 <- pred_scores_l1_mat2["DC", ]
mono_scores_mat2 <- pred_scores_l1_mat2["Mono", ]
dc_scores_mat3 <- pred_scores_l1_mat3["DC", ]
mono_scores_mat3 <- pred_scores_l1_mat3["Mono", ]

dens_mat2 <- MASS::kde2d(dc_scores_mat2,mono_scores_mat2, n = 100)
dens_mat3 <- MASS::kde2d(dc_scores_mat3,mono_scores_mat3, n = 100)
# Plot as image
par(mfrow = c(2, 1))
image(dens_mat2, col = topo.colors(20),main="Mat2", xlab = "DC score",ylab="Mono score")
image(dens_mat3, col = topo.colors(20),main="Mat3", xlab = "DC score",ylab="Mono score")
rm(main.ann)
gc()

#Calculate fractions of cells in reannotated clusters
main <- readRDS("main-mat2vmat3-woclusters-reannotated-azimuth-mapped.rds")
group_col <- "condition"
stack_col <- "annotated_clusters"
metadata_df <- main@meta.data
p <- ggplot(metadata_df, aes_string(x = group_col, fill = stack_col)) +
  geom_bar(position = "fill") +                   
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c("black","darkgrey","lightgray")) +
  labs(
    x = group_col,
    y = "% cells",
    fill = stack_col
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )
p
rm(main)
gc()

#Differential expression analysis for mat3 vs mat2 in mono-like and dc-like clusters
#all cells
main <- readRDS("main-mat2vmat3-woclusters-reannotated-azimuth-mapped.rds")
Idents(main) <- "condition"
de.mat3vsmat2 <- FindMarkers(main, ident.1 = "YS-G-IL4-DC-Mat3", ident.2 = "YS-G-IL4-DC-In",logfc.threshold = 0.5)
write.csv(de.mat3vsmat2,file.path("main-de.mat3vsmat2.csv"))

#dc-like
selected <- subset(main, subset = annotated_clusters == "DC-like")
de.mat3vsmat2.dclike <- FindMarkers(selected, ident.1 = "YS-G-IL4-DC-Mat3", ident.2 = "YS-G-IL4-DC-In",logfc.threshold = 0.5)
write.csv(de.mat3vsmat2.dclike,file.path("main-de.mat3vsmat2.dclike.csv"))

#mono-like
selected <- subset(main, subset = annotated_clusters == "Mono-like")
de.mat3vsmat2.monolike <- FindMarkers(selected, ident.1 = "YS-G-IL4-DC-Mat3", ident.2 = "YS-G-IL4-DC-In",logfc.threshold = 0.5)
write.csv(de.mat3vsmat2.monolike,file.path("main-de.mat3vsmat2.monolike.csv"))
rm(main)
gc()

#Interferon response scores (source of genes is GSEA Hallmark www.gsea-msigdb.org)
main <- readRDS("main-mat2vmat3-woclusters-reannotated-azimuth-mapped.rds")
ifng_response <- c(
  "ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2",
  "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5",
  "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH",
  "CIITA", "CMKLR1", "CMPK2", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "RIGI", "DDX60",
  "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "CMTR1",
  "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HERC6", "HIF1A", "HLA-A", "HLA-B",
  "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30",
  "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3",
  "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1",
  "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7",
  "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1",
  "TMT1B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1",
  "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14",
  "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML",
  "PNP", "PNPT1", "HELZ2", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9",
  "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1",
  "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP",
  "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI",
  "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP",
  "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14",
  "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8",
  "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1"
)
ifna_response <- c(
  "ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74",
  "CMPK2", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1",
  "EPSTI1", "MVB12A", "TENT5A", "CMTR1", "GBP2", "GBP4", "GMPR", "HERC6", "HLA-C",
  "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1",
  "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15",
  "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10", "MX1", "NCOA7",
  "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1",
  "PNPT1", "HELZ2", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", "RIPK2",
  "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2",
  "TAP1", "TDRD7", "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26",
  "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS1"
)
main <- AddModuleScore(object = main, features = list(ifng_response), name = "ifng_response")
main <- AddModuleScore(object = main, features = list(ifna_response), name = "ifna_response")
p1 <- FeaturePlot(object = main, features = "ifng_response1",order=T)
p2 <- FeaturePlot(object = main, features = "ifna_response1",order=T)
p1 + p2 + plot_layout(ncol=2)
ggsave("plots/main_umap_mat2vmat3-woclusters_ifn-ag_response.pdf",width=300,height=100, units="mm") #factor 4 vs what is shown in save image

cols =  c("black","darkgrey","lightgray")
p1 <- VlnPlot(main, features = "ifng_response1", split.by = "annotated_clusters", group.by = "condition",pt.size=F, cols=cols) + ylab("Module score")
p2 <- VlnPlot(main, features = "ifna_response1", split.by = "annotated_clusters", group.by = "condition",pt.size=F, cols=cols) + ylab("Module score")
p1 + p2 + plot_layout(ncol=2)
ggsave("plots/main_violin_umap_mat2vmat3-woclusters_ifn-ag_response.pdf",width=150,height=75, units="mm") #factor 4 vs what is shown in save image

#Differentiation markers
main <- readRDS("main-mat2vmat3-woclusters-reannotated-azimuth-mapped.rds")

#Specific gene list
markers <- list(
  DCs= c("FLT3"),
  Monocytes = c("CSF1R","CD14", "FCGR3A", "ITGAM", "CD33", "CD68", "FCGR1A"),
  cDC1s = c("THBD", "CLEC9A", "XCR1", "CADM1", "IRF8", "BATF3", "TLR3"),
  cDC2s = c("CD1C", "ITGAX","SIRPA", "FCER1A"),
  Maturation = c("RELB", "CD83"),
  Migration = c("CCR7", "MYO1G", "FSCN1", "MARCKS", "MARCKSL1"),
  Co_stim = c("CD80", "CD86", "CD40", "ICOSLG", "TNFSF4", "TNFSF9"),
  Co_inhib = c("CD274", "PDCD1LG2", "CD276", "VTCN1", "VSIR", "LGALS9"),
  Cytokines = c("IL12A", "IL12B", "IL6", "IFNA1", "IFNB1", "TNF", "IL23B", "IL27", "EBI3", "IL15","IFNG","IL10","TFGB1"),
  Chemokines = c("CCL3", "CCL4", "CCL5", "CCL17", "CCL22", "CXCL9", "CXCL10"),
  mRegDC = c("IL4RA", "IL4I1", "BCL2", "BIRC3", "LAMP3", "IDO1")
)

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_mat2vmat3-wocluster_immune_markers_zscore.pdf",width=400,height=100, units="mm") #factor 4 vs what is shown in save image
DotPlot(object = main, features=markers,scale=F) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_mat2vmat3-wocluster_immune_markers_avgexpr.pdf",width=400,height=100, units="mm") #factor 4 vs what is shown in save image

#############################################################################################################################
######################### Genotype analysis for mat3 (with low information cluster removal) ######################################################################
#############################################################################################################################
#Read qc'd experiments
hxdgrb05 <- readRDS(file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxdgrb05-qc-cells.rds"))

#Filter conditons and merge experiments
hxdgrb05[["batch"]] <- "hxdgrb05"
dat <- hxdgrb05
exclude <- grepl("YS-G-IL4-DC-In",dat[[]]$condition)
dat <- dat[,!exclude]

#Filter out low quality cells from clusters 4 and 5
ind <- readRDS(file.path("main-mat2vsmat3-lowQcluster-cellbarcodes.rds"))
ind <- colnames(dat) %in% ind
dat <- dat[,!ind]

#Rename genotypes
unique(dat$sample)
dat[[]]$genotype <- ifelse(grepl("8086",dat[[]]$sample),"CD80-CD86-qKO",
                           ifelse(grepl("L1L2",dat[[]]$sample),"PDL1-PDL2-qKO",
                                  ifelse(dat[[]]$genotype == "dKO","dKO","WT")))

rm(hxdgrb05)
gc()

#Preprocessing and UMAP
all.genes <- rownames(dat)
dat <- NormalizeData(dat,verbose=F) %>% FindVariableFeatures(verbose=F) %>% ScaleData(verbose=F,features=all.genes) %>% RunPCA(verbose=F)
dat <- RunUMAP(dat,dims=1:30,reduction='pca')
dat <- FindNeighbors(dat, reduction = "pca", dims = 1:30)
dat <- FindClusters(dat)

#Save integrated data
saveRDS(dat,file.path("main-mat3-wgeno-woclusters.rds"))
rm(dat)
gc()

#Plot integration results as UMAP
dat <- readRDS("main-mat3-wgeno-woclusters.rds")
set.seed(1)
cols <- DiscretePalette(4, shuffle = FALSE)
p1 <- DimPlot(dat,reduction = "umap", shuffle=T,cols=cols,group.by=c("genotype"),split.by = "condition")
cols <- DiscretePalette(9, shuffle = FALSE)
p2 <- DimPlot(dat, reduction = "umap", shuffle=F, group.by = c("seurat_clusters"),cols=cols,pt.size=1)

p1
ggsave("plots/main_umap_mat3-wgeno_genotypes_woclusters.pdf",width=120,height=75, units="mm") #factor 4 vs what is shown in save image
p2
ggsave("plots/main_umap_mat3-wgeno_seuratclusters_woclusters.pdf",width=90,height=75, units="mm") #factor 4 vs what is shown in save image
rm(dat)
gc()

#Find cluster markers
main <- readRDS("main-mat3-wgeno-woclusters.rds")
allmarkers <- FindAllMarkers(main, only.pos = TRUE)
allmarkers <- allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01)
write.csv(allmarkers,file.path("main-mat3-wgeno-woclusters-cluster-markers.csv"))
rm(main)
gc()

#Calculate fractions of cells in seurat clusters
main <- readRDS("main-mat3-wgeno-woclusters.rds")
group_col <- "genotype"
stack_col <- "seurat_clusters"
metadata_df <- main@meta.data
set.seed(1)
cols <- DiscretePalette(9, shuffle = FALSE)
p <- ggplot(metadata_df, aes_string(x = group_col, fill = stack_col)) +
  geom_bar(position = "fill") +                   
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = cols) +
  labs(
    x = group_col,
    y = "% cells",
    fill = stack_col
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )
p
rm(main)
gc()

#Immune markers per seurat cluster
main <- readRDS("main-mat3-wgeno-woclusters.rds")
markers <- list(
  DCs= c("FLT3"),
  Monocytes = c("CSF1R","CD14", "FCGR3A", "ITGAM", "CD33", "CD68", "FCGR1A"),
  cDC1s = c("THBD", "CLEC9A", "XCR1", "CADM1", "IRF8", "BATF3", "TLR3"),
  cDC2s = c("CD1C", "ITGAX","SIRPA", "FCER1A"),
  Maturation = c("RELB", "CD83"),
  Migration = c("CCR7", "MYO1G", "FSCN1", "MARCKS", "MARCKSL1"),
  Co_stim = c("CD80", "CD86", "CD40", "ICOSLG", "TNFSF4", "TNFSF9"),
  Co_inhib = c("CD274", "PDCD1LG2", "CD276", "VTCN1", "VSIR", "LGALS9"),
  Cytokines = c("IL12A", "IL12B", "IL6", "IFNA1", "IFNB1", "TNF", "IL23B", "IL27", "EBI3", "IL15","IFNG","IL10","TFGB1"),
  Chemokines = c("CCL3", "CCL4", "CCL5", "CCL17", "CCL22", "CXCL9", "CXCL10"),
  mRegDC = c("IL4RA", "IL4I1", "BCL2", "BIRC3", "LAMP3", "IDO1")
)

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers,scale=F,cluster.idents = T) + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/main_dotplot_mat3-wgeno-wocluster_immune_markers_avgexpr.pdf",width=400,height=200, units="mm") #factor 4 vs what is shown in save image

#Cell type inference (Azimuth)
# main <- readRDS("main-mat3-wgeno-woclusters.rds")
# main.ann <- RunAzimuth(main,reference="pbmcref")
# #Save integrated data
# saveRDS(main.ann,file.path("main-mat3-wgeno-woclusters-azimuth-mapped.rds"))
main.ann <- readRDS(file.path("main-mat3-wgeno-woclusters-azimuth-mapped.rds"))

p1 <- FeaturePlot(main.ann, reduction = "umap", order=T, features = "predictionscorecelltypel1_Mono")
p2 <- FeaturePlot(main.ann, reduction = "umap", order=T, features = "predictionscorecelltypel1_DC")
p1 + p2 + plot_layout(ncol = 2)
ggsave("plots/main_umap_mat3-wgeno-woclusters_azimuth_mapped.pdf",width=300,height=100, units="mm") #factor 4 vs what is shown in save image
rm(main.ann)
gc()

#Violin plot representation of DC score
main.ann <- readRDS(file.path("main-mat3-wgeno-woclusters-azimuth-mapped.rds"))
set.seed(1)
cols <- DiscretePalette(9, shuffle = F)
p1 <- VlnPlot(main.ann, features = "predictionscorecelltypel1_DC", pt.size=F, cols=cols) + ylab("DC prediction score")
p2 <- VlnPlot(main.ann, features = "predictionscorecelltypel1_Mono", pt.size=F, cols=cols) + ylab("Mono prediction score")
p1 + p2 + plot_layout(ncol=2)
ggsave("plots/main_violin_umap_mat3-wgeno-woclusters_azimuth_mapped.pdf",width=100,height=75, units="mm") #factor 4 vs what is shown in save image
rm(main.ann)
gc()

#Define new clusters based on predicted cell type
main.ann <- readRDS(file.path("main-mat3-wgeno-woclusters-azimuth-mapped.rds"))

# Access the prediction scores for cell type level 1
pred_scores_l1 <- GetAssayData(object = main.ann, assay = "prediction.score.celltype.l1", layer = "data")

# Extract DC and Mono scores
if ("DC" %in% rownames(pred_scores_l1)) {
  dc_scores <- pred_scores_l1["DC", ]
} else {
  stop("Failed to find scores.")
}

if ("Mono" %in% rownames(pred_scores_l1)) {
  mono_scores <- pred_scores_l1["Mono", ]
} else {
  stop("Failed to find scores.")
}

par(mfrow = c(1, 1))
plot(dc_scores,mono_scores)
dens <- MASS::kde2d(dc_scores,mono_scores, n = 100)
# Plot as image
image(dens, col = topo.colors(20), xlab = "DC score",ylab="Mono score")
# Add contours
contour(dens, add = TRUE)

main.ann$DC_scores <- dc_scores
main.ann$Mono_scores <- mono_scores
main.ann[[]]$annotated_clusters <- ifelse(main.ann[[]]$Mono_scores > 0.6 & main.ann[[]]$DC_scores < 0.15, "Mono-like",
                                          ifelse(main.ann[[]]$Mono_scores < 0.3 & main.ann[[]]$DC_scores > 0.2,"DC-like","Other"))

#main[[]]$annotated_clusters <- ifelse(main[[]]$seurat_clusters %in% c(2,5,6), "mono-like","dc-like")
saveRDS(main.ann,file.path("main-mat3-wgeno-woclusters-reannotated-azimuth-mapped.rds"))
rm(main.ann)
gc()


#Calculate fractions of cells in reannotated clusters
main <- readRDS("main-mat3-wgeno-woclusters-reannotated-azimuth-mapped.rds")
group_col <- "genotype"
stack_col <- "annotated_clusters"
metadata_df <- main@meta.data
p <- ggplot(metadata_df, aes_string(x = group_col, fill = stack_col)) +
  geom_bar(position = "fill") +                   
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c("black","darkgrey","lightgray")) +
  labs(
    x = group_col,
    y = "% cells",
    fill = stack_col
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )
p
rm(main)
gc()


