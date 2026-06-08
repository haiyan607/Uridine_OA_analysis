# ==========================================================
# 01_create_single_sample_seurat_objects.R
# GSE220243: create single-sample Seurat objects from raw 10X files
# ==========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(dplyr)
library(stringr)
library(Matrix)
library(here)

# ==========================================================
# 1. Paths
# ==========================================================
raw_dir <- here("data", "raw", "GSE220243_RAW")
processed_dir <- here("data", "processed")
tables_dir <- here("results", "tables")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(raw_dir)) {
  stop("Cannot find raw data directory: ", raw_dir)
}

# ==========================================================
# 2. Sample metadata
# ==========================================================
sample_info <- data.frame(
  GSM_ID = paste0("GSM", 6797148:6797173),
  sample_name = c(
    paste0("CartNorm", 1:6),
    paste0("CartOA", 1:6),
    "MenNorm1AVAS",
    "MenNorm1VAS",
    paste0("MenNorm", 2:7),
    paste0("MenOA", 1:6)
  ),
  stringsAsFactors = FALSE
)

sample_info$sample_id <- paste(sample_info$GSM_ID, sample_info$sample_name, sep = "_")

sample_info$tissue <- ifelse(
  str_detect(sample_info$sample_name, "^Cart"),
  "Cartilage",
  "Meniscus"
)

sample_info$group <- ifelse(
  str_detect(sample_info$sample_name, "Norm"),
  "Normal",
  "OA"
)

sample_info$tissue_group <- paste(sample_info$tissue, sample_info$group, sep = "_")

sample_info$donor_id <- case_when(
  str_detect(sample_info$sample_name, "^CartNorm") ~ str_extract(sample_info$sample_name, "\\d+"),
  str_detect(sample_info$sample_name, "^CartOA")   ~ str_extract(sample_info$sample_name, "\\d+"),
  str_detect(sample_info$sample_name, "^MenNorm")  ~ str_extract(sample_info$sample_name, "\\d+"),
  str_detect(sample_info$sample_name, "^MenOA")    ~ str_extract(sample_info$sample_name, "\\d+"),
  TRUE ~ NA_character_
)

sample_info$donor_label <- paste(
  sample_info$tissue,
  sample_info$group,
  paste0("Donor", sample_info$donor_id),
  sep = "_"
)

sample_info$region <- case_when(
  str_detect(sample_info$sample_name, "AVAS$") ~ "Avascular",
  str_detect(sample_info$sample_name, "VAS$")  ~ "Vascular",
  sample_info$tissue == "Meniscus"             ~ "Unspecified",
  TRUE                                         ~ "Not_applicable"
)

sample_info$group <- factor(sample_info$group, levels = c("Normal", "OA"))
sample_info$tissue <- factor(sample_info$tissue, levels = c("Cartilage", "Meniscus"))
sample_info$tissue_group <- factor(
  sample_info$tissue_group,
  levels = c("Cartilage_Normal", "Cartilage_OA", "Meniscus_Normal", "Meniscus_OA")
)

if (nrow(sample_info) != 26) {
  stop("sample_info does not contain 26 samples.")
}

# ==========================================================
# 3. Organize standard 10X folders if needed
# ==========================================================
for (i in seq_len(nrow(sample_info))) {
  sid <- sample_info$sample_id[i]
  sample_dir <- file.path(raw_dir, sid)
  
  standard_files <- file.path(
    sample_dir,
    c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  )
  
  if (dir.exists(sample_dir) && all(file.exists(standard_files))) {
    next
  }
  
  source_files <- file.path(
    raw_dir,
    c(
      paste0(sid, "_barcodes.tsv.gz"),
      paste0(sid, "_features.tsv.gz"),
      paste0(sid, "_matrix.mtx.gz")
    )
  )
  
  if (!all(file.exists(source_files))) {
    missing_files <- basename(source_files[!file.exists(source_files)])
    stop(
      paste0(
        "Missing raw files for sample ", sid, ":\n",
        paste(missing_files, collapse = "\n")
      )
    )
  }
  
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  ok <- file.copy(from = source_files, to = standard_files, overwrite = FALSE)
  if (!all(ok)) {
    stop("Failed to copy raw files for sample: ", sid)
  }
}

# ==========================================================
# 4. Read 10X and create single-sample Seurat objects
# ==========================================================
datalist <- vector("list", nrow(sample_info))
names(datalist) <- sample_info$sample_id

for (i in seq_len(nrow(sample_info))) {
  sid <- sample_info$sample_id[i]
  message("[", i, "/", nrow(sample_info), "] Reading: ", sid)
  
  dir_10x <- file.path(raw_dir, sid)
  
  my_data <- Read10X(
    data.dir = dir_10x,
    gene.column = 2,
    unique.features = TRUE
  )
  
  if (is.list(my_data)) {
    if ("Gene Expression" %in% names(my_data)) {
      my_data <- my_data[["Gene Expression"]]
    } else {
      my_data <- my_data[[1]]
    }
  }
  
  seu <- CreateSeuratObject(
    counts = my_data,
    project = sid,
    min.cells = 3,
    min.features = 250
  )
  
  seu$GSM_ID <- sample_info$GSM_ID[i]
  seu$sample_id <- sid
  seu$sample_name <- sample_info$sample_name[i]
  seu$tissue <- as.character(sample_info$tissue[i])
  seu$group <- as.character(sample_info$group[i])
  seu$tissue_group <- as.character(sample_info$tissue_group[i])
  seu$donor_id <- sample_info$donor_id[i]
  seu$donor_label <- sample_info$donor_label[i]
  seu$region <- sample_info$region[i]
  
  datalist[[i]] <- seu
  rm(my_data, seu)
  gc(verbose = FALSE)
}

# ==========================================================
# 5. Summaries
# ==========================================================
sample_dimensions <- do.call(
  rbind,
  lapply(names(datalist), function(sid) {
    obj <- datalist[[sid]]
    data.frame(
      sample_id = sid,
      genes = nrow(obj),
      cells = ncol(obj),
      tissue = unique(obj$tissue),
      group = unique(obj$group),
      tissue_group = unique(obj$tissue_group),
      donor_label = unique(obj$donor_label),
      region = unique(obj$region),
      stringsAsFactors = FALSE
    )
  })
)

rownames(sample_dimensions) <- NULL

# ==========================================================
# 6. Save outputs
# ==========================================================
saveRDS(
  datalist,
  file = file.path(processed_dir, "GSE220243_26samples_datalist_raw.rds")
)

write.csv(
  sample_info,
  file = file.path(tables_dir, "GSE220243_sample_metadata.csv"),
  row.names = FALSE
)

write.csv(
  sample_dimensions,
  file = file.path(tables_dir, "GSE220243_raw_sample_dimensions.csv"),
  row.names = FALSE
)

message("Finished: single-sample Seurat objects created.")