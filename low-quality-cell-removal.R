#Load packages and data
library(Seurat)
library(ggplot2)
#library(SeuratDisk)
#library(dplyr)
#library(patchwork)
#library(harmony)

#Working directly
setwd("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/")

#Read in preprocessed seurat objects
infiles <- list.files(path = "/blellochlab/data1/deniz/analysis/hx-analyses", pattern = "-nodups-seurat.rds", recursive = TRUE, full.names = TRUE)
infiles <- infiles[grep("hxrb03|hxrb04|hxdgrb05",infiles)]

##################################################################################################################################
# HXRB03 (HX DC differentiation Experiment #1)
dat <- readRDS(infiles[grep("hxrb03",infiles)])

#Remove unrelated conditions
exclude <- grepl("inhib",dat[[]]$condition)
dat <- dat[,!exclude]
dat[[]]$condition <- factor(dat[[]]$condition,
                            levels = c("iPS","YS-Mye","YS-G-IL4-DC","YS-G-IL4-DC-LT",
                                       "YS-G-IL4-DC-In","YS-G-Flt3L-DC","YS-G-Flt3L-DC-LT","YS-G-Flt3L-DC-In"))

ind <- grepl("Flt3|LT",dat[[]]$condition)
main <- dat[,!ind]

#Main figure data
##Quality metrics
cols <- DiscretePalette(5, shuffle = FALSE)
VlnPlot(main, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "condition",pt.size=F,cols = cols)
ggsave("plots/hxrb03_main_nFeat_nRNA_PercMito.pdf",width=150,height=75, units="mm") #factor 4 vs what is shown in save image

high_mito <- quantile(main[["percent.mt"]]$percent.mt, probs=0.98)
low_ft <- quantile(main[["nFeature_RNA"]]$nFeature_RNA, probs=0.02)
high_ft <- quantile(main[["nFeature_RNA"]]$nFeature_RNA, probs=0.98)

pdf("plots/hxrb03_main_cutoffs.pdf", width = 100 / 25.4 / 2, height = 200 / 25.4 / 2)  # Convert mm to inches
par(mar = c(7, 4, 4, 2),las = 2)
barplot(c(high_mito,low_ft,high_ft),names=c("High Mt DNA","Low Features","High Features"),
        log="y", ylim=c(1,100000), ylab="Cutoff")
dev.off()  # Close the graphics device

#All data
##Quality metrics
cols <- DiscretePalette(8, shuffle = FALSE)
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "condition",pt.size=F,cols = cols)
ggsave("plots/hxrb03_all_nFeat_nRNA_PercMito.pdf",width=200,height=75, units="mm") #factor 4 vs what is shown in save image

## Remove unwanted cells
dat <- subset(dat, subset = percent.mt < high_mito & nFeature_RNA > low_ft & nFeature_RNA < high_ft) # ~4129 cells remaining

#Export
saveRDS(dat,file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxrb03-qc-cells.rds"))

##################################################################################################################################
# HXRB04 (HX DC differentiation Experiment #2)
dat <- readRDS(infiles[grep("hxrb04",infiles)])

#Remove unrelated conditions
exclude <- grepl("PBMC",dat[[]]$condition) 
dat <- dat[,!exclude]
dat[[]]$condition <- factor(dat[[]]$condition,
                            levels = c("iPS","YS-Mye","YS-G-IL4-DC","YS-G-IL4-DC-In","YS-G-Flt3L-DC","YS-G-Flt3L-DC-In",
                                       "YS-G-DC","YS-G-DC-In"))

ind <- grepl("Flt3|G-DC|LT",dat[[]]$condition)
main <- dat[,!ind]

#Main figure data
##Quality metrics
cols <- DiscretePalette(4, shuffle = FALSE)
VlnPlot(main, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "condition",pt.size=F,cols = cols)
ggsave("plots/hxrb04_main_nFeat_nRNA_PercMito.pdf",width=120,height=75, units="mm") #factor 4 vs what is shown in save image

high_mito <- quantile(main[["percent.mt"]]$percent.mt, probs=0.98)
low_ft <- quantile(main[["nFeature_RNA"]]$nFeature_RNA, probs=0.02)
high_ft <- quantile(main[["nFeature_RNA"]]$nFeature_RNA, probs=0.98)

pdf("plots/hxrb04_main_cutoffs.pdf", width = 100 / 25.4 / 2, height = 200 / 25.4 / 2)  # Convert mm to inches
par(mar = c(7, 4, 4, 2),las = 2)
barplot(c(high_mito,low_ft,high_ft),names=c("High Mt DNA","Low Features","High Features"),
        log="y", ylim=c(1,100000), ylab="Cutoff")
dev.off()  # Close the graphics device

#All figure data
##Quality metrics
cols <- DiscretePalette(8, shuffle = FALSE)
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "condition",pt.size=F,cols = cols)
ggsave("plots/hxrb04_all_nFeat_nRNA_PercMito.pdf",width=200,height=75, units="mm") #factor 4 vs what is shown in save image

## Remove unwanted cells
dat <- subset(dat, subset = percent.mt < high_mito & nFeature_RNA > low_ft & nFeature_RNA < high_ft) # ~6520 cells remaining

#Export
saveRDS(dat,file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxrb04-qc-cells.rds"))

##################################################################################################################################
# HXDGRB05 (HX DC differentiation Experiment mat3 vs mat2 #3)
dat <- readRDS(infiles[grep("hxdgrb05",infiles)])

dat[[]]$genotype <- "WT"
dat[[]]$genotype <- ifelse(grepl("dKO",dat[[]]$sample),"dKO",
                           ifelse(grepl("8086",dat[[]]$sample),"qKO","WT"))
dat[[]]$genotype <- factor(dat[[]]$genotype,
                            levels = c("WT","dKO","qKO"))

main <- dat
rm(dat)
#Main fidat#Main figure data
##Quality metrics
cols <- DiscretePalette(4, shuffle = FALSE)
VlnPlot(main, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "condition",pt.size=F,cols = cols)
ggsave("plots/hxdgrb05_main_nFeat_nRNA_PercMito.pdf",width=120,height=75, units="mm") #factor 4 vs what is shown in save image

VlnPlot(main, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "genotype",pt.size=F,cols = cols)
ggsave("plots/hxdgrb05_main_nFeat_nRNA_PercMito.pdf",width=120,height=75, units="mm") #factor 4 vs what is shown in save image

high_mito <- quantile(main[["percent.mt"]]$percent.mt, probs=0.98)
low_ft <- quantile(main[["nFeature_RNA"]]$nFeature_RNA, probs=0.02)
high_ft <- quantile(main[["nFeature_RNA"]]$nFeature_RNA, probs=0.98)

pdf("plots/hxdgrb05_main_cutoffs.pdf", width = 100 / 25.4 / 2, height = 200 / 25.4 / 2)  # Convert mm to inches
par(mar = c(7, 4, 4, 2),las = 2)
barplot(c(high_mito,low_ft,high_ft),names=c("High Mt DNA","Low Features","High Features"),
        log="y", ylim=c(1,100000), ylab="Cutoff")
dev.off()  # Close the graphics device

#No. cells per condition
pdf("plots/hxdgrb05_cell_numbers.pdf", width = 100 / 25.4 / 2, height = 200 / 25.4 / 2)  # Convert mm to inches
par(mar = c(7, 4, 4, 2),las = 2)
barplot(table(main[[]]$sample))
dev.off()  # Close the graphics device

## Remove unwanted cells
dat <- subset(main, subset = percent.mt < high_mito & nFeature_RNA > low_ft & nFeature_RNA < high_ft) # ~4260 cells remaining

#Export
saveRDS(dat,file.path("/blellochlab/data1/deniz/analysis/hx-analyses/analysis-for-manuscript/hxdgrb05-qc-cells.rds"))






