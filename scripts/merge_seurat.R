#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- args[1]
output_merged <- args[2]
output_plot <- args[3]

library(Seurat)
library(dplyr)
library(ggplot2)

# Load sample metadata
samples <- read.csv(samples_csv)

# Load each sample's Seurat object and tag metadata
objs <- list()
for (i in 1:nrow(samples)) {
  s <- samples$sample[i]
  cond <- samples$condition[i]
  obj <- readRDS(paste0("results/", s, "/seurat.rds"))
  obj$sample <- s
  obj$condition <- cond
  objs[[s]] <- obj
}

# Merge all samples
merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))

# Integration (optional, based on size)
merged <- SCTransform(merged, verbose = FALSE)
merged <- RunPCA(merged, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:10)
merged <- FindClusters(merged, resolution = 0.5)
merged <- RunUMAP(merged, dims = 1:10)

# Save merged object
saveRDS(merged, file = output_merged)

# Plot
pdf(file = output_plot, width = 7, height = 6)
DimPlot(merged, group.by = "sample", label = TRUE) + ggtitle("Merged UMAP by Sample")
DimPlot(merged, group.by = "condition", label = TRUE) + ggtitle("Merged UMAP by Condition")
dev.off()

# Perform DE analysis: condition (global)
de_condition <- FindMarkers(merged, ident.1 = "treatment", ident.2 = "control", group.by = "condition", verbose = FALSE)

# Save DE results
write.csv(de_condition, file = "results/de_condition.csv")

# Perform DE analysis: within each cluster
Idents(merged) <- "seurat_clusters"
cluster_ids <- levels(merged)

# Store results in a list
de_by_cluster <- list()

for (cluster in cluster_ids) {
  # Subset cluster
  cells_in_cluster <- WhichCells(merged, idents = cluster)
  sub <- subset(merged, cells = cells_in_cluster)

  # Set condition as identity
  Idents(sub) <- sub$condition

  # Run DE
  de <- tryCatch({
    FindMarkers(sub, ident.1 = "treatment", ident.2 = "control", verbose = FALSE)
  }, error = function(e) {
    return(data.frame())
  })

  de$gene <- rownames(de)
  de$cluster <- cluster
  de_by_cluster[[cluster]] <- de
}

# Combine and save
de_combined <- do.call(rbind, de_by_cluster)
write.csv(de_combined, file = "results/de_by_cluster.csv", row.names = FALSE)

