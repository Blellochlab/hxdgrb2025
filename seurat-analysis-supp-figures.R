#Load libraries
library(Seurat)
library(dplyr)
library(gridExtra)
library(Azimuth)
library(ggplot2)
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
ind <- grepl("YS-Mye|iPS",dat[[]]$condition)
dat <- dat[,!ind]
unique(dat[[]]$condition)

dat[[]]$condition <- factor(dat[[]]$condition,
                            levels = c("YS-G-DC","YS-G-DC-In",
                                       "YS-G-IL4-DC","YS-G-IL4-DC-LT","YS-G-IL4-DC-In",
                                       "YS-G-Flt3L-DC","YS-G-Flt3L-DC-LT","YS-G-Flt3L-DC-In"))
rm(hxrb03,hxrb04)
gc()

#Non-integrated preprocessing and UMAP
dat <- NormalizeData(dat,verbose=F) %>% FindVariableFeatures(verbose=F) %>% ScaleData(verbose=F) %>% RunPCA(verbose=F)

dat <- RunUMAP(dat,dims=1:30,reduction='pca',reduction.name = "umap.unintegrated")

#Integrate data
dat <- IntegrateLayers(object = dat, method = CCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.cca",verbose = FALSE) #takes ~1min
dat[["RNA"]] <- JoinLayers(dat[["RNA"]])
dat <- FindNeighbors(dat, reduction = "integrated.cca", dims = 1:30)
dat <- FindClusters(dat, resolution = 1)
dat <- RunUMAP(dat, dims = 1:30, reduction = "integrated.cca")

#Save integrated data
saveRDS(dat,file.path("supp-integrated.rds"))
rm(dat)
gc()

#Plot integration results as UMAP
dat <- readRDS("supp-integrated.rds")
set.seed(345)
cols <- c("darkgrey","black")
p1 <- DimPlot(dat, reduction = "umap.unintegrated", shuffle=T, group.by = c("batch"),cols=cols,pt.size=1)
p2 <- DimPlot(dat, reduction = "umap", shuffle=T, group.by = c("batch"),cols=cols,pt.size=1)
cols <- DiscretePalette(15, shuffle = TRUE)
p3 <- DimPlot(dat, reduction = "umap", shuffle=F, group.by = c("seurat_clusters"),cols=cols,pt.size=1)
cols <- colorRampPalette(c("orange", "dodgerblue"))(8)
cols <- c("gray","black","lightblue","dodgerblue","darkblue","gold","orange","red4")
p4 <- DimPlot(dat, reduction = "umap", shuffle=T, group.by = c("condition"),cols=cols,pt.size = 1)

p1 + p2
ggsave("plots/supp_umap_batchremoval.pdf",width=200,height=75, units="mm") #factor 4 vs what is shown in save image
p3
ggsave("plots/supp_umap_seuratclusters.pdf",width=90,height=75, units="mm") #factor 4 vs what is shown in save image
p4
ggsave("plots/supp_umap_conditions.pdf",width=120,height=75, units="mm") #factor 4 vs what is shown in save image

#Plot all clusters
plot_list <- list()
for (i in factor(0:14)) {
  CellsToPlot <- colnames(dat)[dat$seurat_clusters == i]
  p <- DimPlot(dat, cells.highlight = CellsToPlot, reduction="umap", cols.highlight = "dodgerblue",
               cols = "gray", sizes.highlight = 0.2, order = TRUE, shuffle = TRUE, label=F) + ggtitle(i) + NoLegend()  
  plot_list[[i]] <- p
}
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 4))

#Plot all clusters
plot_list <- list()
for (i in unique(dat[[]]$condition)) {
  CellsToPlot <- colnames(dat)[dat$condition == i]
  p <- DimPlot(dat, cells.highlight = CellsToPlot, reduction="umap", cols.highlight = "dodgerblue",
               cols = "gray", sizes.highlight = 0.2, order = TRUE, shuffle = TRUE, label=F) + ggtitle(i) + NoLegend()  
  plot_list[[i]] <- p
}
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 3))

#Find cluster markers
main <- readRDS(file.path("supp-integrated.rds"))
allmarkers <- FindAllMarkers(main, only.pos = TRUE)
allmarkers <- allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01)
write.csv(allmarkers,file.path("supp-all-sign-sr-cluster-markers.csv"))

#Find condition markers
main <- readRDS(file.path("supp-integrated.rds"))
main <- SetIdent(main,value="condition")
allmarkers <- FindAllMarkers(main, only.pos = TRUE)
allmarkers <- allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01)
write.csv(allmarkers,file.path("supp-all-sign-condition-markers.csv"))

# #Cell type inference (Azimuth)
# main <- readRDS(file.path("supp-integrated.rds"))
# main.ann <- RunAzimuth(main,reference="pbmcref")
# #Save integrated data
# saveRDS(main.ann,file.path("supp-integrated-azimuth-mapped.rds"))
set.seed(345)
main.ann <- readRDS(file.path("supp-integrated-azimuth-mapped.rds"))

cols <- DiscretePalette(8, shuffle = TRUE)
p1 <- DimPlot(main.ann, reduction = "umap", group.by = "predicted.celltype.l1", label = F, shuffle=F, order=T, cols=cols, label.size = 1)
# cols <- DiscretePalette(13, shuffle = TRUE)
# p2 <- DimPlot(main.ann, reduction = "umap", group.by = "predicted.celltype.l2", label = F, shuffle=F, order=T, cols=cols, label.size = 1)
cols <- c("gray","black","lightblue","dodgerblue","darkblue","gold","orange","red4")
p2 <- DimPlot(main.ann, reduction = "umap", group.by = "condition", label = F, cols=cols, label.size = 1)
p3 <- FeaturePlot(main.ann, reduction = "umap", features = "predictionscorecelltypel1_DC")
p4 <- FeaturePlot(main.ann, reduction = "umap", features = "predictionscorecelltypel1_Mono")

p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
ggsave("plots/supp_umap_azimuth_mapped.pdf",width=225,height=150, units="mm") #factor 4 vs what is shown in save image

#Violin plot representation of DC score
p1 <- VlnPlot(main.ann, features = "predictionscorecelltypel1_DC", group.by = "condition",pt.size=F,cols = cols) + ylab("DC prediction score")
p2 <- VlnPlot(main.ann, features = "predictionscorecelltypel1_Mono", group.by = "condition",pt.size=F,cols = cols) + ylab("Mono prediction score")
p1 + p2 + plot_layout(ncol = 1)
ggsave("plots/supp_violin_azimuth_DC_Monoscore.pdf",width=150,height=150, units="mm") #factor 4 vs what is shown in save image


#############################################################################################################################
#############################################################################################################################
#################################################Old scripts#################################################################
#############################################################################################################################
#############################################################################################################################


#Differentiation markers
main <- readRDS(file.path("supp-integrated.rds"))

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
    "CCR2",
    "FCGR1A",
    "HLA-DRA"
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
    "CD86",
    "CD274",
    "CCR7",
    "LAMP3",
    "CD83",
    "HLA-DRB1",
    "IDO1",
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
condition_colors <- condition_colors[condition_data]
condition_annotation <- HeatmapAnnotation(
  condition = condition_data,
  col = list(condition = c("iPS" = "red", "YS-Mye" = "green", "YS-G-IL4-DC" = "blue", "YS-G-IL4-DC-In" = "purple"))
)

col_fun <- circlize::colorRamp2(c(min(zscore_data),max(zscore_data)), c("white","blue"))

Heatmap(zscore_data,
        cluster_rows = TRUE,   # Cluster genes
        cluster_columns = TRUE,  # Don't cluster conditions
        show_row_names = TRUE,  # Show gene names
        show_column_names = FALSE,  # Show condition names
        name = "Z-score",
        col=col_fun,
        top_annotation = condition_annotation,
        row_names_gp = gpar(fontsize = 5),
        heatmap_legend_param = list(title = "Z-score")
)

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers, cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))

##Dotplot by marker genes by culture condition
Idents(main) <- main$condition
DotPlot(object = main, features=markers, cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))

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

#Remove cluster 6 and 4 cells based on low information and doublet bias
cellBCs <- readRDS(file.path("supp-IL4-DC-In-integrated-cluster4and6-cellbarcodes.rds"))
ind <- !colnames(dat) %in% cellBCs
dat <- dat[,ind]

#Non-integrated preprocessing and UMAP
dat <- NormalizeData(dat,verbose=F) %>% FindVariableFeatures(verbose=F) %>% ScaleData(verbose=F) %>% RunPCA(verbose=F)
ElbowPlot(dat) #choose 15 dims
ggsave("plots/elbow_supp_IL4-DC-In_unintegrated.pdf",width=50,height=50, units="mm") #factor 4 vs what is shown in save image

dat <- RunUMAP(dat,dims=1:30,reduction='pca',reduction.name = "umap.unintegrated")

#Non-int. plot
set.seed(1)
cols <- DiscretePalette(5, shuffle = TRUE)

#Integrate data
dat <- IntegrateLayers(object = dat, method = CCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.cca",verbose = FALSE) #takes ~1min
dat[["RNA"]] <- JoinLayers(dat[["RNA"]])
dat <- FindNeighbors(dat, reduction = "integrated.cca", dims = 1:30)
dat <- FindClusters(dat, resolution = 1)
dat <- RunUMAP(dat, dims = 1:30, reduction = "integrated.cca")

#Save integrated data
saveRDS(dat,file.path("supp-IL4-DC-In-integrated.rds"))

#Plot integration results as UMAP
dat <- readRDS("supp-IL4-DC-In-integrated.rds")
set.seed(345)
cols <- c("darkgrey","black")
p1 <- DimPlot(dat, reduction = "umap.unintegrated", shuffle=T, group.by = c("batch"),cols=cols,pt.size=1)
p2 <- DimPlot(dat, reduction = "umap", shuffle=T, group.by = c("batch"),cols=cols,pt.size=1)
cols <- DiscretePalette(5, shuffle = TRUE)
p3 <- DimPlot(dat, reduction = "umap", shuffle=F, group.by = c("seurat_clusters"),cols=cols,pt.size=1)

p1 + p2
ggsave("plots/supp_umap_IL4-DC-In_batchremoval.pdf",width=200,height=75, units="mm") #factor 4 vs what is shown in save image
p3
ggsave("plots/supp_umap_IL4-DC-In_seuratclusters.pdf",width=90,height=75, units="mm") #factor 4 vs what is shown in save image

#Find cluster markers
main <- readRDS("supp-IL4-DC-In-integrated.rds")
allmarkers <- FindAllMarkers(main, only.pos = TRUE)
allmarkers <- allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01)
write.csv(allmarkers,file.path("supp-IL4-DC-In-sign-sr-cluster-markers.csv"))

# #Cell type inference (Azimuth)
# main <- readRDS("supp-IL4-DC-In-integrated.rds")
# main.ann <- RunAzimuth(main,reference="pbmcref")
##Save integrated data
saveRDS(main.ann,file.path("supp-IL4-DC-In-integrated-azimuth-mapped.rds"))
main.ann <- readRDS(file.path("supp-IL4-DC-In-integrated-azimuth-mapped.rds"))

cols <- DiscretePalette(3, shuffle = TRUE)
p1 <- DimPlot(main.ann, reduction = "umap", group.by = "predicted.celltype.l1", label = F, shuffle=T,order=F, cols=cols, label.size = 1)
p2 <- FeaturePlot(main.ann, reduction = "umap", order=T, features = "predictionscorecelltypel1_DC")
p1 + p2 + plot_layout(ncol = 2)
ggsave("plots/supp_umap_IL4-DC-In_azimuth_mapped.pdf",width=200,height=150, units="mm") #factor 4 vs what is shown in save image

# #Identify cells that are part of clusters 6 and 4
# dat <- readRDS("supp-IL4-DC-In-integrated.rds")
# ind <- dat$seurat_clusters %in% c(4,6)
# saveRDS(colnames(dat)[ind],file.path("supp-IL4-DC-In-integrated-cluster4and6-cellbarcodes.rds"))

#Differentiation markers
main <- readRDS("supp-IL4-DC-In-integrated.rds")

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
    "CCR2",
    "FCGR1A",
    "HLA-DRA"
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
    "CD86",
    "CD274",
    "CCR7",
    "LAMP3",
    "CD83",
    "HLA-DRB1",
    "IDO1",
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
set.seed(214)
condition_data <- main@meta.data$seurat_clusters 
condition_annotation <- HeatmapAnnotation(
  condition = condition_data
)

col_fun <- circlize::colorRamp2(c(min(zscore_data),max(zscore_data)), c("white","blue"))

Heatmap(zscore_data,
        cluster_rows = TRUE,   # Cluster genes
        cluster_columns = TRUE,  # Don't cluster conditions
        show_row_names = TRUE,  # Show gene names
        show_column_names = FALSE,  # Show condition names
        name = "Z-score",
        col=col_fun,
        top_annotation = condition_annotation,
        row_names_gp = gpar(fontsize = 5),
        heatmap_legend_param = list(title = "Z-score")
)

#Feature, mt, counts per cluster
VlnPlot(main, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "seurat_clusters",pt.size=F)

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers, cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))

#Specific gene list
markers <- list(
  Cytokines = c("IL12A", "IL12B", "IL6", "IFNA1", "IFNB1", "TNF", "IL23A", "IL23B", "IL1B", "IL27", "EBI3", "IL15"),
  Co_stim = c("CD80", "CD86", "CD40", "ICOSLG", "TNFSF4", "TNFSF9"),
  Co_inhib = c("CD274", "PDCD1LG2", "CD276", "VTCN1", "VSIR", "LGALS9"),
  Monocytes = c("CD14", "FCGR3A", "ITGAM", "CD33", "CD68", "CSF1R", "CCR2", "FCGR1A"),
  cDC1s = c("THBD", "CLEC9A", "XCR1", "CADM1", "IRF8", "BATF3", "TLR3"),
  cDC2s = c("CD1C", "ITGAX","SIRPA", "FCER1A"),
  Maturation = c("RELB", "CD83"),
  Regulatory = c("FAS", "ALDH1A2", "SOCS1", "SOCS2"),
  Migration = c("CCR7", "MYO1G", "FSCN1", "MARCKS", "MARCKSL1"),
  Chemokines = c("CCL3", "CCL4", "CCL5", "CCL17", "CCL19", "CCL22", "CCL25", "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CXCL13"),
  mRegDC = c("IL4RA", "IL4I1", "BCL2", "BIRC3", "LAMP3", "IDO1"),
  M2_polarization = c("ARG1", "CD163", "IL10", "TGFB1", "TGFB2", "TGFB3", "CHI3L1", "PPARG", "VEGFA", "IL4", "IL13")
)

##Dotplot by marker genes by cluster
DotPlot(object = main, features=markers, cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))
