# ==========================================================
# 01_differential_protein_venn_plots.R
#
# Purpose:
#   Generate up-regulated and down-regulated four-set Venn
#   diagrams from the differential protein abundance dataset.
#
# Input:
#   data/raw/proteomics/20240205-DAP.xlsx
#
# Output:
#   results/figures/
#   results/tables/
# ==========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(readxl)
library(dplyr)
library(ggvenn)
library(ggplot2)
library(patchwork)
library(here)

# ==========================================================
# 1. Analysis parameters
# ==========================================================

P_CUTOFF <- 0.05

UP_R100_CUTOFF <- 2
UP_RATIO_CUTOFF <- 1
UP_R30_CUTOFF <- 1

DOWN_R100_CUTOFF <- 1
DOWN_RATIO_CUTOFF <- 1
DOWN_R30_CUTOFF <- 1

# ==========================================================
# 2. Project paths
# ==========================================================

input_file <- here(
  "data",
  "raw",
  "proteomics",
  "20240205-DAP.xlsx"
)

figure_dir <- here(
  "results",
  "figures"
)

table_dir <- here(
  "results",
  "tables"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Input file was not found:\n",
    input_file
  )
}

# ==========================================================
# 3. Import and validate data
# ==========================================================

raw_data <- read_excel(input_file)

required_columns <- c(
  "Genes",
  "30uM-RATIO",
  "30PVALUE",
  "100uM-RATIO",
  "100PVALUE"
)

missing_columns <- setdiff(
  required_columns,
  colnames(raw_data)
)

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing:\n",
    paste(missing_columns, collapse = "\n")
  )
}

# ==========================================================
# 4. Clean data and resolve duplicated proteins
#
# For duplicated protein symbols, retain the row with the
# largest absolute fold change across the 30 and 100 uM
# treatment conditions.
# ==========================================================

protein_data <- raw_data %>%
  transmute(
    Genes = trimws(as.character(Genes)),
    R30 = suppressWarnings(
      as.numeric(`30uM-RATIO`)
    ),
    P30 = suppressWarnings(
      as.numeric(`30PVALUE`)
    ),
    R100 = suppressWarnings(
      as.numeric(`100uM-RATIO`)
    ),
    P100 = suppressWarnings(
      as.numeric(`100PVALUE`)
    )
  ) %>%
  filter(
    !is.na(Genes),
    Genes != "",
    !is.na(R30),
    !is.na(R100),
    R30 > 0,
    R100 > 0
  ) %>%
  mutate(
    change_score = pmax(
      abs(log2(R30)),
      abs(log2(R100)),
      na.rm = TRUE
    )
  ) %>%
  group_by(Genes) %>%
  slice_max(
    order_by = change_score,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    R100_div_R30 = R100 / R30
  ) %>%
  filter(
    is.finite(R100_div_R30)
  )

write.csv(
  protein_data,
  file = file.path(
    table_dir,
    "proteomics_cleaned_unique_proteins.csv"
  ),
  row.names = FALSE
)

# ==========================================================
# 5. Construct four sets for the up-regulated Venn diagram
#
# Confirmed criteria:
#   P30 < 0.05
#   R100 > 2
#   R100 / R30 > 1
#   R30 > 1
# ==========================================================

up_sets <- list(
  "P30 < 0.05" = unique(
    protein_data$Genes[
      !is.na(protein_data$P30) &
        protein_data$P30 < P_CUTOFF
    ]
  ),
  "Ratio100 > 2" = unique(
    protein_data$Genes[
      protein_data$R100 > UP_R100_CUTOFF
    ]
  ),
  "Ratio100/30 > 1" = unique(
    protein_data$Genes[
      protein_data$R100_div_R30 > UP_RATIO_CUTOFF
    ]
  ),
  "Ratio30 > 1" = unique(
    protein_data$Genes[
      protein_data$R30 > UP_R30_CUTOFF
    ]
  )
)

# ==========================================================
# 6. Construct four sets for the down-regulated Venn diagram
#
# Confirmed criteria:
#   P100 < 0.05
#   R100 < 1
#   R100 / R30 < 1
#   R30 < 1
# ==========================================================

down_sets <- list(
  "P100 < 0.05" = unique(
    protein_data$Genes[
      !is.na(protein_data$P100) &
        protein_data$P100 < P_CUTOFF
    ]
  ),
  "Ratio100 < 1" = unique(
    protein_data$Genes[
      protein_data$R100 < DOWN_R100_CUTOFF
    ]
  ),
  "Ratio100/30 < 1" = unique(
    protein_data$Genes[
      protein_data$R100_div_R30 < DOWN_RATIO_CUTOFF
    ]
  ),
  "Ratio30 < 1" = unique(
    protein_data$Genes[
      protein_data$R30 < DOWN_R30_CUTOFF
    ]
  )
)

# ==========================================================
# 7. Export set members and four-way intersections
# ==========================================================

export_sets <- function(set_list, output_file) {
  set_table <- bind_rows(
    lapply(
      names(set_list),
      function(set_name) {
        data.frame(
          Set = set_name,
          Genes = set_list[[set_name]],
          stringsAsFactors = FALSE
        )
      }
    )
  )

  write.csv(
    set_table,
    file = output_file,
    row.names = FALSE
  )
}

export_sets(
  up_sets,
  file.path(
    table_dir,
    "venn_upregulated_set_members.csv"
  )
)

export_sets(
  down_sets,
  file.path(
    table_dir,
    "venn_downregulated_set_members.csv"
  )
)

up_core <- Reduce(
  intersect,
  up_sets
)

down_core <- Reduce(
  intersect,
  down_sets
)

write.csv(
  data.frame(
    Genes = sort(up_core)
  ),
  file = file.path(
    table_dir,
    "venn_upregulated_four_way_intersection.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Genes = sort(down_core)
  ),
  file = file.path(
    table_dir,
    "venn_downregulated_four_way_intersection.csv"
  ),
  row.names = FALSE
)

set_size_summary <- bind_rows(
  data.frame(
    Direction = "Up",
    Set = names(up_sets),
    Protein_number = lengths(up_sets)
  ),
  data.frame(
    Direction = "Down",
    Set = names(down_sets),
    Protein_number = lengths(down_sets)
  )
)

write.csv(
  set_size_summary,
  file = file.path(
    table_dir,
    "venn_set_size_summary.csv"
  ),
  row.names = FALSE
)

# ==========================================================
# 8. Draw Venn diagrams
# ==========================================================

up_venn <- ggvenn(
  up_sets,
  fill_color = c(
    "#F3A6A6",
    "#A8DDB5",
    "#9ECAE1",
    "#C5B0D5"
  ),
  fill_alpha = 0.55,
  stroke_size = 0.6,
  set_name_size = 4.5,
  text_size = 4
) +
  ggtitle("Up-regulated proteins") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    )
  )

down_venn <- ggvenn(
  down_sets,
  fill_color = c(
    "#F3A6A6",
    "#A8DDB5",
    "#9ECAE1",
    "#C5B0D5"
  ),
  fill_alpha = 0.55,
  stroke_size = 0.6,
  set_name_size = 4.5,
  text_size = 4
) +
  ggtitle("Down-regulated proteins") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    )
  )

combined_venn <- up_venn + down_venn +
  plot_layout(
    ncol = 2
  )

# ==========================================================
# 9. Save figures
# ==========================================================

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_venn_upregulated.pdf"
  ),
  plot = up_venn,
  width = 7,
  height = 6,
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_venn_downregulated.pdf"
  ),
  plot = down_venn,
  width = 7,
  height = 6,
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_venn_up_down_combined.pdf"
  ),
  plot = combined_venn,
  width = 14,
  height = 6,
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_venn_up_down_combined.png"
  ),
  plot = combined_venn,
  width = 14,
  height = 6,
  units = "in",
  dpi = 600
)

message("Venn analysis completed successfully.")
