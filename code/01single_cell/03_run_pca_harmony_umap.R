# ==========================================================
# 03_run_pca_harmony_umap.R
# GSE220243: PCA, Harmony, UMAP
# ==========================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(1234)

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(here)

# ==========================================================
# 1. Paths
# ==========================================================
processed_dir <- here("data", "processed")
figures_dir <- here("results", "figures")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(processed_dir, "GSE220243_QC_normalized_hvg_scaled.rds")
output_file <- file.path(processed_dir, "GSE220243_QC_normalized_hvg_scaled_harmony_umap.rds")

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

# ==========================================================
# 2. Read object
# ==========================================================
sce <- readRDS(input_file)
DefaultAssay(sce) <- "RNA"

# ==========================================================
# 3. PCA
# ==========================================================
if (length(VariableFeatures(sce)) == 0) {
  stop("VariableFeatures not found. Please run script 02 first.")
}

sce <- RunPCA(
  object = sce,
  features = VariableFeatures(sce),
  npcs = 30,
  verbose = TRUE
)

pdf(file.path(figures_dir, "12_ElbowPlot.pdf"), width = 6, height = 4.5)
print(ElbowPlot(sce, ndims = 30))
dev.off()

# ==========================================================
# 4. Harmony
# batch correction by sample_id
# ==========================================================
if (!"sample_id" %in% colnames(sce@meta.data)) {
  stop("sample_id not found in metadata.")
}

sce <- RunHarmony(
  object = sce,
  group.by.vars = "sample_id",
  reduction = "pca",
  dims.use = 1:20,
  assay.use = "RNA",
  verbose = TRUE
)

# ==========================================================
# 5. UMAP
# ==========================================================
sce <- RunUMAP(
  object = sce,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  seed.use = 1234,
  verbose = TRUE
)

# ==========================================================
# 6. Diagnostic UMAP plots
# ==========================================================
p_tissue <- DimPlot(
  object = sce,
  reduction = "umap",
  group.by = "tissue",
  pt.size = 0.2
) +
  ggtitle("UMAP by tissue") +
  theme_classic(base_size = 12)

p_group <- DimPlot(
  object = sce,
  reduction = "umap",
  group.by = "group",
  pt.size = 0.2
) +
  ggtitle("UMAP by group") +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(figures_dir, "13_UMAP_by_tissue.pdf"),
  plot = p_tissue,
  width = 6,
  height = 5
)

ggsave(
  filename = file.path(figures_dir, "14_UMAP_by_group.pdf"),
  plot = p_group,
  width = 6,
  height = 5
)

# ==========================================================
# 7. Save object
# ==========================================================
saveRDS(
  sce,
  file = output_file
)

message("Finished: PCA, Harmony and UMAP completed.")