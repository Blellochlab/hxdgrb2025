# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales) # Needed for percent_format()

# Set working dir
setwd("/home/goekbugd/scratch/intern/ucsf-blelloch-collab/hx-analyses/analysis-for-manuscript/")

# Define Reference Path and Label Column
ref_path <- "fg-lab-rds/Reference_mDC_Verse.rds" 
ref_label_col <- "Mega.Clusters" 
ref <- readRDS(ref_path)

# Define the specific DC-VERSE scores of interest
target_scores <- c(
  "prediction.score.IFN.DC", "prediction.score.DC3", "prediction.score.LTB.DC", 
  "prediction.score.DC2", "prediction.score.Prolif.DC", "prediction.score.CCR7..mDC", 
  "prediction.score.DC1", "prediction.score.IL1B.DC", "prediction.score.preDC", 
  "prediction.score.CD207.DC"
)

# Create short names for clean Violin Plot titles (removes "prediction.score.")
short_target_scores <- gsub("prediction.score.", "", target_scores)

# Clean IDs that correspond to the scores above for subsetting 
# (Note: Check 'unique(ref$Mega.Clusters)' to ensure exact string matching, e.g., "CCR7+ mDC")
target_ids <- c("IFN DC", "DC3", "LTB DC", "DC2", "Prolif DC", "CCR7+ mDC", "DC1", "IL1B DC", "preDC", "CD207 DC")

# ==============================================================================
# ================================= MAT 2 ======================================
# ==============================================================================
print("Starting Mat2 Pipeline...")
mat2_path <- "main-IL4-DC-In-integrated-woclusters.rds" 
dat_mat2 <- readRDS(mat2_path)

# Preprocessing & Mapping
DefaultAssay(dat_mat2) <- "RNA"
dat_mat2 <- NormalizeData(dat_mat2, verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE)

query_var_genes <- VariableFeatures(dat_mat2)
anchoring_features <- intersect(query_var_genes, rownames(ref))

anchors_mat2 <- FindTransferAnchors(
  reference = ref, query = dat_mat2, normalization.method = "LogNormalize",
  reference.reduction = "pca", dims = 1:30, features = anchoring_features 
)

predictions_mat2 <- TransferData(anchorset = anchors_mat2, refdata = ref[[ref_label_col]][,1], dims = 1:30)
dat_mat2 <- AddMetaData(dat_mat2, metadata = predictions_mat2)

# Copy scores to short-named columns for clean plotting titles
dat_mat2@meta.data[, short_target_scores] <- dat_mat2@meta.data[, target_scores]

# 1. DimPlot: DC-Verse annotated UMAP
set.seed(123)
n_mat2 <- length(unique(dat_mat2$predicted.id))
cols_mat2 <- DiscretePalette(n_mat2, shuffle = TRUE)

p_dim_mat2 <- DimPlot(dat_mat2, reduction = "umap", shuffle = TRUE, group.by = "predicted.id", cols = cols_mat2, pt.size = 1, label = FALSE, repel = TRUE) + 
  ggtitle("Mat2: DC-VERSE Annotations")
p_dim_mat2
ggsave("plots/Mat2_DCVerse_DimPlot.pdf", plot = p_dim_mat2, width = 100, height = 75, units = "mm")

# 2. Output Percent Cells and Mean/Median Scores
mat2_summary <- dat_mat2@meta.data %>%
  group_by(predicted.id) %>%
  summarise(
    n_cells = n(),
    across(all_of(target_scores), list(mean = mean, median = median), .names = "{.col}_{.fn}")
  ) %>%
  mutate(percent_cells = (n_cells / sum(n_cells)) * 100)
str(mat2_summary)
write.csv(mat2_summary, "Mat2_DCVerse_SummaryStats.csv", row.names = FALSE)

# 3. Violin plot of target DC-VERSE scores across all cells by seurat_clusters
p_vln_mat2 <- VlnPlot(dat_mat2, features = short_target_scores, group.by = "seurat_clusters", pt.size = 0, ncol = 5, cols = c("#993F00","#191919")) & ylim(0, 1)
p_vln_mat2
ggsave("plots/Mat2_DCVerse_VlnPlots_bySeuratCluster.pdf", plot = p_vln_mat2, width = 150, height = 150, units = "mm")

print("Mat2 Processing Complete.")
rm(dat_mat2, anchors_mat2, predictions_mat2, mat2_summary)
gc()

# ==============================================================================
# ================================= MAT 3 ======================================
# ==============================================================================
print("Starting Mat3 Pipeline...")
mat3_path <- "main-mat3-wgeno-woclusters.rds"
dat_mat3 <- readRDS(mat3_path)

# Preprocessing & Mapping
DefaultAssay(dat_mat3) <- "RNA"
dat_mat3 <- NormalizeData(dat_mat3, verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE)

query_var_genes <- VariableFeatures(dat_mat3)
anchoring_features <- intersect(query_var_genes, rownames(ref))

anchors_mat3 <- FindTransferAnchors(
  reference = ref, query = dat_mat3, normalization.method = "LogNormalize",
  reference.reduction = "pca", dims = 1:30, features = anchoring_features 
)

predictions_mat3 <- TransferData(anchorset = anchors_mat3, refdata = ref[[ref_label_col]][,1], dims = 1:30)
dat_mat3 <- AddMetaData(dat_mat3, metadata = predictions_mat3)

# Copy scores to short-named columns for clean plotting titles
dat_mat3@meta.data[, short_target_scores] <- dat_mat3@meta.data[, target_scores]

# 1. DimPlot: DC-Verse annotated UMAP
n_mat3 <- length(unique(dat_mat3$predicted.id))
cols_mat3 <- DiscretePalette(n_mat3, shuffle = TRUE)

p_dim_mat3 <- DimPlot(dat_mat3, reduction = "umap", shuffle = TRUE, group.by = "predicted.id", cols = cols_mat3, pt.size = 1, label = FALSE, repel = TRUE) + 
  ggtitle("Mat3: DC-VERSE Annotations")
p_dim_mat3
ggsave("plots/Mat3_DCVerse_DimPlot.pdf", plot = p_dim_mat3, width = 120, height = 75, units = "mm")

# 2. Stacked Bar Plot: Percent cells by genotype in each DC verse cluster
df_plot <- dat_mat3@meta.data %>% filter(predicted.id %in% target_ids)

p_bar <- ggplot(df_plot, aes(x = genotype, fill = predicted.id)) +
  geom_bar(position = "fill") +                    
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c(cols_mat3)) +
  labs(x = "DC-VERSE Annotations", y = "% cells", fill = "Genotype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_bar
ggsave("plots/Mat3_DCVerse_Genotype_Fractions.pdf", plot = p_bar, width = 100, height = 75, units = "mm")

# 3. Output Mean/Median prediction score in each DC verse cluster
mat3_summary <- dat_mat3@meta.data %>%
  group_by(predicted.id) %>%
  summarise(
    n_cells = n(),
    across(all_of(target_scores), list(mean = mean, median = median), .names = "{.col}_{.fn}")
  ) %>%
  mutate(percent_cells = (n_cells / sum(n_cells)) * 100)

write.csv(mat3_summary, "Mat3_DCVerse_SummaryStats.csv", row.names = FALSE)

# 4. Violin plot of target DC-VERSE scores across all cells by seurat_clusters
p_vln_mat3 <- VlnPlot(dat_mat3, features = short_target_scores, group.by = "seurat_clusters", pt.size = 0, ncol = 5, cols=c("#F0A0FF", "#0075DC", "#993F00" ,"#4C005C", "#191919", "#005C31", "#2BCE48", "#FFCC99", "#808080")) & ylim(0, 1)
p_vln_mat3
ggsave("plots/Mat3_DCVerse_VlnPlots_bySeuratCluster.pdf", plot = p_vln_mat3, width = 400, height = 150, units = "mm")
