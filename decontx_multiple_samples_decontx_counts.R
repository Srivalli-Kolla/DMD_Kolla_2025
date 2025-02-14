# Install required packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("celda", "singleCellTK", "enrichR"))

# Load necessary libraries
library(celda)
library(scater)
library(singleCellTK)
library(ggplot2)

# Define sample directories
sampleDirs <- list.files("/Volumes/One Touch/Wusi/dmd_raw_data/raw_data", full.names = TRUE)

# Output directory for plots
output_dir <- "/Volumes/One Touch/Wusi/DMD_Kolla_2025/decontx/plots_decontx_counts"
dir.create(output_dir, showWarnings = FALSE)

# Define marker genes
markers <- list(
  "Cardiomyocytes" = c("Ttn", "Myh6", "Myh7", "Actc1", "Tnnt2", "Ryr2", "Cacna1c", "Gja1"),
  "Fibroblasts" = c("Dcn", "Col1a1", "Col3a1", "Pdgfra", "Postn", "Fn1", "Vim", "Thbs1"),
  "Smooth Muscle Cells" = c("Acta2", "Tagln", "Myh11", "Cnn1", "Lmod1", "Smtn", "Des", "Myl9"),
  "Lymphoid Cells" = c(
    "Cd3e", "Cd4", "Cd8a", "Trac", "Tcrb",  # T cells
    "Cd19", "Cd79a", "Ms4a1", "Ighm", "Cd22" # B cells
  ),
  "Myeloid Cells" = c(
    "Adgre1", "Cd68", "Csf1r", "Mrc1", "Itgam", # Macrophages
    "S100a8", "S100a9", "Lcn2", "Mpo", "Elane", # Neutrophils
    "Itgax", "Cd86", "H2-Ab1", "H2-Aa", "H2-Eb1" # Dendritic cells
  )
)

# Loop through each sample
for (sample in sampleDirs) {
  sample_name <- basename(sample)
  cat("\nProcessing:", sample_name, "\n")
  
  # Import data
  sce.raw <- importCellRanger(sampleDirs = sample, dataType = "raw")
  
  # Remove low-quality cells
  sce.raw <- sce.raw[, colSums(assay(sce.raw, "counts")) > 100]
  
  # Run DecontX
  decontx_sce <- decontX(sce.raw)
  
  # Replace raw counts with decontaminated counts for further analysis
  assay(decontx_sce, "counts") <- assay(decontx_sce, "decontXcounts")
  
  # Generate and save contamination plot
  p1 <- plotDecontXContamination(decontx_sce)
  ggsave(filename = file.path(output_dir, paste0(sample_name, "_contamination.png")), plot = p1)
  
  # Get UMAP and save cluster plot
  umap <- reducedDim(decontx_sce, "decontX_UMAP")
  p2 <- plotDimReduceCluster(x = decontx_sce$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2])
  ggsave(filename = file.path(output_dir, paste0(sample_name, "_clusters.png")), plot = p2)
  
  # Normalize and replace row names with feature names
  norm_sce <- logNormCounts(decontx_sce)
  rownames(norm_sce) <- rowData(norm_sce)$feature_name
  
  # Save marker expression plot
  p3 <- plotDimReduceFeature(logcounts(norm_sce), dim1 = umap[, 1], dim2 = umap[, 2], features = unlist(markers))
  ggsave(filename = file.path(output_dir, paste0(sample_name, "_marker_expression.png")), plot = p3)
  
  # Save marker percentage bar plots using only decontXcounts
  p4 <- plotDecontXMarkerPercentage(norm_sce, markers = markers, assayName = "counts") # Now "counts" is decontXcounts
  ggsave(filename = file.path(output_dir, paste0(sample_name, "_marker_percentage.png")), plot = p4)
  
  # Save marker expression violin plots
  p5 <- plotDecontXMarkerExpression(norm_sce, markers = markers, ncol = 3)
  ggsave(filename = file.path(output_dir, paste0(sample_name, "_marker_expression_violin.png")), plot = p5)
  
  cat("Finished processing:", sample_name, "\n")
}