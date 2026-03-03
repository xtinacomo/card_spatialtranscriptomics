#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
platform <- args[2]
input_path <- args[3]
output_rds <- args[4]
output_violin <- args[5]
output_umap <- args[6]
output_spatial <- args[7]

library(Seurat)
library(ggplot2)

# Load data
if (platform == "visium") {
  seurat_obj <- Load10X_Spatial(data.dir = input_path)
} else if (platform == "xenium") {
  data <- read.delim(input_path, sep="\t", row.names = 1)
  seurat_obj <- CreateSeuratObject(counts = data)
} else {
  stop("Unsupported platform")
}

# Mito % for QC
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Violin plot of QC
pdf(file = output_violin, width = 8, height = 6)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Basic filtering
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 20)

# Standard pipeline
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# UMAP plot
pdf(file = output_umap, width = 6, height = 5)
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle(sample)
dev.off()

# SpatialDimPlot (only for Visium)
if (platform == "visium") {
  pdf(file = output_spatial, width = 6, height = 6)
  SpatialDimPlot(seurat_obj, label = TRUE) + ggtitle(paste0(sample, " Spatial"))
  dev.off()
} else {
  file.create(output_spatial)  
}

# Save Seurat object
saveRDS(seurat_obj, file = output_rds)
