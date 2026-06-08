# ==========================================================
# 02_quality_control_and_normalization.R
# GSE220243: QC, merge, normalization, HVG, scaling
# ==========================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(here)

# ==========================================================
# 1. Paths
# ==========================================================
processed_dir <- here("data", "processed")
figures_dir <- here("results", "figures")
tables_dir <- here("results", "tables")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(processed_dir, "GSE220243_26samples_datalist_raw.rds")
if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

# ==========================================================
# 2. Read objects
# ==========================================================
datalist <- readRDS(input_file)

if (length(datalist) != 26) {
  warning("datalist contains ", length(datalist), " objects, not 26.")
}

required_metadata <- c("sample_id", "sample_name", "tissue", "group", "tissue_group")
metadata_check <- vapply(
  datalist,
  function(obj) all(required_metadata %in% colnames(obj@meta.data)),
  logical(1)
)
if (!all(metadata_check)) {
  stop("Some objects are missing required metadata.")
}

# ==========================================================
# 3. Add QC metrics
# ==========================================================
for (i in seq_along(datalist)) {
  obj <- datalist[[i]]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.Ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  datalist[[i]] <- obj
  rm(obj)
}

# ==========================================================
# 4. Cell counts before QC
# ==========================================================
raw_count_df <- bind_rows(
  lapply(datalist, function(obj) {
    data.frame(
      sample_id = unique(as.character(obj$sample_id)),
      sample_name = unique(as.character(obj$sample_name)),
      tissue = unique(as.character(obj$tissue)),
      group = unique(as.character(obj$group)),
      tissue_group = unique(as.character(obj$tissue_group)),
      raw_count = ncol(obj),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  raw_count_df,
  file = file.path(tables_dir, "01_GSE220243_cell_count_before_QC.csv"),
  row.names = FALSE
)

# ==========================================================
# 5. Merge before QC for visualization
# ==========================================================
sce_before <- merge(
  x = datalist[[1]],
  y = datalist[-1],
  add.cell.ids = names(datalist),
  project = "GSE220243_before_QC"
)

sce_before$tissue_group <- factor(
  sce_before$tissue_group,
  levels = c("Cartilage_Normal", "Cartilage_OA", "Meniscus_Normal", "Meniscus_OA")
)

sce_before$QC_group <- factor(
  sce_before$tissue_group,
  levels = c("Cartilage_Normal", "Cartilage_OA", "Meniscus_Normal", "Meniscus_OA"),
  labels = c("Normal cartilage", "OA cartilage", "Normal meniscus", "OA meniscus")
)

qc_before_4groups <- VlnPlot(
  object = sce_before,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.Ribo"),
  group.by = "QC_group",
  pt.size = 0.01,
  raster = TRUE,
  ncol = 2
) &
  theme_classic(base_size = 12) &
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(figures_dir, "02_QC_before_4groups.pdf"),
  plot = qc_before_4groups,
  width = 14,
  height = 10
)

qc_summary_before <- sce_before@meta.data %>%
  group_by(tissue_group) %>%
  summarise(
    cell_number = n(),
    nFeature_median = median(nFeature_RNA, na.rm = TRUE),
    nFeature_Q1 = quantile(nFeature_RNA, 0.25, na.rm = TRUE),
    nFeature_Q3 = quantile(nFeature_RNA, 0.75, na.rm = TRUE),
    nFeature_P99 = quantile(nFeature_RNA, 0.99, na.rm = TRUE),
    nCount_median = median(nCount_RNA, na.rm = TRUE),
    nCount_Q1 = quantile(nCount_RNA, 0.25, na.rm = TRUE),
    nCount_Q3 = quantile(nCount_RNA, 0.75, na.rm = TRUE),
    nCount_P99 = quantile(nCount_RNA, 0.99, na.rm = TRUE),
    percent_mt_median = median(percent.mt, na.rm = TRUE),
    percent_mt_Q1 = quantile(percent.mt, 0.25, na.rm = TRUE),
    percent_mt_Q3 = quantile(percent.mt, 0.75, na.rm = TRUE),
    percent_mt_P98 = quantile(percent.mt, 0.98, na.rm = TRUE),
    percent_Ribo_median = median(percent.Ribo, na.rm = TRUE),
    percent_Ribo_Q1 = quantile(percent.Ribo, 0.25, na.rm = TRUE),
    percent_Ribo_Q3 = quantile(percent.Ribo, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  qc_summary_before,
  file = file.path(tables_dir, "03_GSE220243_QC_summary_before.csv"),
  row.names = FALSE
)

saveRDS(
  sce_before,
  file = file.path(processed_dir, "GSE220243_merged_before_QC.rds")
)

# ==========================================================
# 6. Sample-wise QC filtering
# ==========================================================
datalist_clean <- vector("list", length(datalist))
names(datalist_clean) <- names(datalist)

qc_threshold_list <- vector("list", length(datalist))
names(qc_threshold_list) <- names(datalist)

qc_cell_count_list <- vector("list", length(datalist))
names(qc_cell_count_list) <- names(datalist)

for (i in seq_along(datalist)) {
  current_sample <- names(datalist)[i]
  obj <- datalist[[i]]
  meta <- obj@meta.data
  
  current_tissue <- unique(as.character(obj$tissue))
  if (length(current_tissue) != 1) {
    stop("Sample ", current_sample, " has multiple tissue values.")
  }
  
  nCount_upper <- as.numeric(quantile(meta$nCount_RNA, probs = 0.99, na.rm = TRUE))
  mt_q98 <- as.numeric(quantile(meta$percent.mt, probs = 0.98, na.rm = TRUE))
  mt_absolute_cap <- ifelse(current_tissue == "Cartilage", 15, 20)
  mt_upper <- min(mt_q98, mt_absolute_cap)
  
  Ribo_lower_observed <- as.numeric(quantile(meta$percent.Ribo, probs = 0.01, na.rm = TRUE))
  Ribo_upper_observed <- as.numeric(quantile(meta$percent.Ribo, probs = 0.99, na.rm = TRUE))
  
  qc_threshold_list[[i]] <- data.frame(
    sample_id = unique(as.character(obj$sample_id)),
    sample_name = unique(as.character(obj$sample_name)),
    tissue = current_tissue,
    group = unique(as.character(obj$group)),
    tissue_group = unique(as.character(obj$tissue_group)),
    nFeature_lower = 500,
    nFeature_upper = 6000,
    nCount_lower = 1000,
    nCount_upper_P99 = nCount_upper,
    percent_mt_P98 = mt_q98,
    percent_mt_absolute_cap = mt_absolute_cap,
    percent_mt_final_upper = mt_upper,
    percent_Ribo_P01_observed = Ribo_lower_observed,
    percent_Ribo_P99_observed = Ribo_upper_observed,
    percent_Ribo_used_for_filtering = FALSE,
    stringsAsFactors = FALSE
  )
  
  keep_cells <- rownames(meta)[
    meta$nFeature_RNA > 500 &
      meta$nFeature_RNA < 6000 &
      meta$nCount_RNA > 1000 &
      meta$nCount_RNA < nCount_upper &
      meta$percent.mt < mt_upper
  ]
  
  if (length(keep_cells) == 0) {
    stop("No cells retained after QC for sample: ", current_sample)
  }
  
  obj_clean <- subset(obj, cells = keep_cells)
  
  qc_cell_count_list[[i]] <- data.frame(
    sample_id = unique(as.character(obj$sample_id)),
    sample_name = unique(as.character(obj$sample_name)),
    tissue = current_tissue,
    group = unique(as.character(obj$group)),
    tissue_group = unique(as.character(obj$tissue_group)),
    raw_count = ncol(obj),
    clean_count = ncol(obj_clean),
    removed_count = ncol(obj) - ncol(obj_clean),
    retention_percent = round(ncol(obj_clean) / ncol(obj) * 100, 2),
    stringsAsFactors = FALSE
  )
  
  datalist_clean[[i]] <- obj_clean
  rm(obj, obj_clean, meta, keep_cells)
  gc(verbose = FALSE)
}

qc_threshold_df <- bind_rows(qc_threshold_list)
write.csv(
  qc_threshold_df,
  file = file.path(tables_dir, "04_GSE220243_QC_thresholds_by_sample.csv"),
  row.names = FALSE
)

qc_cell_count_df <- bind_rows(qc_cell_count_list)
write.csv(
  qc_cell_count_df,
  file = file.path(tables_dir, "05_GSE220243_cell_count_before_after_QC.csv"),
  row.names = FALSE
)

group_retention_summary <- qc_cell_count_df %>%
  group_by(tissue_group) %>%
  summarise(
    sample_number = n(),
    total_raw_cells = sum(raw_count),
    total_clean_cells = sum(clean_count),
    overall_retention_percent = round(total_clean_cells / total_raw_cells * 100, 2),
    mean_sample_retention = round(mean(retention_percent), 2),
    minimum_sample_retention = min(retention_percent),
    maximum_sample_retention = max(retention_percent),
    .groups = "drop"
  )

write.csv(
  group_retention_summary,
  file = file.path(tables_dir, "06_GSE220243_group_retention_summary.csv"),
  row.names = FALSE
)

# ==========================================================
# 7. Plot cell counts before/after QC
# ==========================================================
count_plot_data <- qc_cell_count_df %>%
  select(sample_name, tissue_group, raw_count, clean_count) %>%
  pivot_longer(
    cols = c(raw_count, clean_count),
    names_to = "QC_status",
    values_to = "cell_count"
  ) %>%
  mutate(
    QC_status = factor(
      QC_status,
      levels = c("raw_count", "clean_count"),
      labels = c("Before QC", "After QC")
    )
  )

cell_count_plot <- ggplot(
  count_plot_data,
  aes(x = sample_name, y = cell_count, fill = QC_status)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Sample", y = "Number of cells") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(figures_dir, "07_cell_count_before_after_QC.pdf"),
  plot = cell_count_plot,
  width = 15,
  height = 7
)

# ==========================================================
# 8. Merge after QC
# ==========================================================
datalist <- datalist_clean
rm(datalist_clean, qc_threshold_list, qc_cell_count_list)
gc()

sce <- merge(
  x = datalist[[1]],
  y = datalist[-1],
  add.cell.ids = names(datalist),
  project = "GSE220243_after_QC"
)

if (exists("JoinLayers")) {
  sce <- JoinLayers(sce)
}

sce$tissue_group <- factor(
  sce$tissue_group,
  levels = c("Cartilage_Normal", "Cartilage_OA", "Meniscus_Normal", "Meniscus_OA")
)

sce$QC_group <- factor(
  sce$tissue_group,
  levels = c("Cartilage_Normal", "Cartilage_OA", "Meniscus_Normal", "Meniscus_OA"),
  labels = c("Normal cartilage", "OA cartilage", "Normal meniscus", "OA meniscus")
)

qc_after_4groups <- VlnPlot(
  object = sce,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.Ribo"),
  group.by = "QC_group",
  pt.size = 0.01,
  raster = TRUE,
  ncol = 2
) &
  theme_classic(base_size = 12) &
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(figures_dir, "08_QC_after_4groups.pdf"),
  plot = qc_after_4groups,
  width = 14,
  height = 10
)

qc_summary_after <- sce@meta.data %>%
  group_by(tissue_group) %>%
  summarise(
    cell_number = n(),
    nFeature_median = median(nFeature_RNA, na.rm = TRUE),
    nFeature_Q1 = quantile(nFeature_RNA, 0.25, na.rm = TRUE),
    nFeature_Q3 = quantile(nFeature_RNA, 0.75, na.rm = TRUE),
    nFeature_P99 = quantile(nFeature_RNA, 0.99, na.rm = TRUE),
    nCount_median = median(nCount_RNA, na.rm = TRUE),
    nCount_Q1 = quantile(nCount_RNA, 0.25, na.rm = TRUE),
    nCount_Q3 = quantile(nCount_RNA, 0.75, na.rm = TRUE),
    nCount_P99 = quantile(nCount_RNA, 0.99, na.rm = TRUE),
    percent_mt_median = median(percent.mt, na.rm = TRUE),
    percent_mt_Q1 = quantile(percent.mt, 0.25, na.rm = TRUE),
    percent_mt_Q3 = quantile(percent.mt, 0.75, na.rm = TRUE),
    percent_mt_P98 = quantile(percent.mt, 0.98, na.rm = TRUE),
    percent_Ribo_median = median(percent.Ribo, na.rm = TRUE),
    percent_Ribo_Q1 = quantile(percent.Ribo, 0.25, na.rm = TRUE),
    percent_Ribo_Q3 = quantile(percent.Ribo, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  qc_summary_after,
  file = file.path(tables_dir, "09_GSE220243_QC_summary_after.csv"),
  row.names = FALSE
)

qc_before_after_4groups <- qc_before_4groups / qc_after_4groups +
  plot_annotation(
    title = "GSE220243 quality control",
    subtitle = "Before and after quality control",
    tag_levels = "A"
  )

ggsave(
  filename = file.path(figures_dir, "10_QC_before_after_4groups.pdf"),
  plot = qc_before_after_4groups,
  width = 14,
  height = 20,
  limitsize = FALSE
)

# ==========================================================
# 9. Save QC outputs
# ==========================================================
saveRDS(
  datalist,
  file = file.path(processed_dir, "GSE220243_26samples_datalist_after_QC.rds")
)

saveRDS(
  sce,
  file = file.path(processed_dir, "GSE220243_merged_after_QC.rds")
)

# ==========================================================
# 10. Normalize, HVG, ScaleData
# ==========================================================
sce <- NormalizeData(
  object = sce,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

saveRDS(
  sce,
  file = file.path(processed_dir, "GSE220243_QC_normalized.rds")
)

sce <- FindVariableFeatures(
  object = sce,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)

top20 <- head(VariableFeatures(sce), 20)

variable_plot <- VariableFeaturePlot(sce)
variable_feature_plot <- LabelPoints(
  plot = variable_plot,
  points = top20,
  repel = TRUE,
  size = 3
)

ggsave(
  filename = file.path(figures_dir, "11_top20_variable_features.pdf"),
  plot = variable_feature_plot,
  width = 9,
  height = 7
)

sce <- ScaleData(
  object = sce,
  features = VariableFeatures(sce),
  verbose = TRUE
)

sce$group <- factor(sce$group, levels = c("Normal", "OA"))
sce$tissue <- factor(sce$tissue, levels = c("Cartilage", "Meniscus"))
sce$tissue_group <- factor(
  sce$tissue_group,
  levels = c("Cartilage_Normal", "Cartilage_OA", "Meniscus_Normal", "Meniscus_OA")
)

saveRDS(
  sce,
  file = file.path(processed_dir, "GSE220243_QC_normalized_hvg_scaled.rds")
)

writeLines(
  capture.output(sessionInfo()),
  con = here("docs", "sessionInfo.txt")
)

message("Finished: QC and normalization completed.")