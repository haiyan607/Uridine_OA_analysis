# ==========================================================
# 02_differential_protein_scatter_plot.R
#
# Purpose:
#   Generate the differential protein scatter plot using:
#     X = Ratio100 / Ratio30
#     Y = Ratio100 / Vehicle
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
library(ggplot2)
library(ggrepel)
library(here)

# ==========================================================
# 1. Analysis parameters
# ==========================================================

P_CUTOFF <- 0.05

X_UP_CUTOFF <- 1
Y_UP_CUTOFF <- 2

X_DOWN_CUTOFF <- 1
Y_DOWN_CUTOFF <- 0.5

REMOVE_DISPLAY_OUTLIERS <- TRUE
DISPLAY_X_MAX <- 1.5
DISPLAY_Y_MAX <- 4.5

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
# 4. Data cleaning and coordinate calculation
#
# Only proteins significant at both concentrations are used:
#   P30 < 0.05 and P100 < 0.05
#
# For duplicated protein symbols, retain the row showing the
# largest deviation from 1 across X and Y.
# ==========================================================

scatter_data <- raw_data %>%
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
    !is.na(P30),
    !is.na(P100),
    R30 > 0,
    R100 > 0,
    P30 < P_CUTOFF,
    P100 < P_CUTOFF
  ) %>%
  mutate(
    X = R100 / R30,
    Y = R100,
    Regulation = case_when(
      Y > Y_UP_CUTOFF &
        X > X_UP_CUTOFF ~ "Up",
      Y < Y_DOWN_CUTOFF &
        X < X_DOWN_CUTOFF ~ "Down",
      TRUE ~ NA_character_
    ),
    change_score = pmax(
      abs(log2(X)),
      abs(log2(Y)),
      na.rm = TRUE
    )
  ) %>%
  filter(
    is.finite(X),
    is.finite(Y),
    !is.na(Regulation)
  ) %>%
  group_by(Genes) %>%
  slice_max(
    order_by = change_score,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

scatter_data$Regulation <- factor(
  scatter_data$Regulation,
  levels = c("Down", "Up")
)

write.csv(
  scatter_data,
  file = file.path(
    table_dir,
    "scatter_filtered_best_per_protein.csv"
  ),
  row.names = FALSE
)

# ==========================================================
# 5. Prepare plotting data
#
# The optional display filter reproduces the manuscript plot
# by excluding only points located beyond both display limits.
# The complete filtered table remains saved above.
# ==========================================================

if (REMOVE_DISPLAY_OUTLIERS) {
  scatter_plot_data <- scatter_data %>%
    filter(
      !(X > DISPLAY_X_MAX & Y > DISPLAY_Y_MAX)
    )
} else {
  scatter_plot_data <- scatter_data
}

write.csv(
  scatter_plot_data,
  file = file.path(
    table_dir,
    "scatter_data_used_for_plotting.csv"
  ),
  row.names = FALSE
)

# ==========================================================
# 6. Define proteins to label
# ==========================================================

label_map <- c(
  "SRSF1" = "SRSF1",
  "DDX5" = "DDX5",
  "FUS" = "FUS",
  "RBM39" = "RBM39",
  "SLC16A1" = "SLC16a1",
  "IGF2R" = "IGF2R"
)

scatter_plot_data <- scatter_plot_data %>%
  mutate(
    Gene_upper = toupper(Genes)
  )

label_data <- scatter_plot_data %>%
  filter(
    Gene_upper %in% names(label_map)
  ) %>%
  mutate(
    Label = unname(
      label_map[Gene_upper]
    ),
    label_color = ifelse(
      Regulation == "Up",
      "#B22222",
      "#1F78B4"
    )
  )

write.csv(
  label_data %>%
    dplyr::select(
      Genes,
      Gene_upper,
      Label,
      Regulation,
      X,
      Y
    ) %>%
    arrange(
      Regulation,
      desc(Y),
      desc(X)
    ),
  file = file.path(
    table_dir,
    "scatter_labeled_protein_coordinates.csv"
  ),
  row.names = FALSE
)

# ==========================================================
# 7. Draw scatter plot
# ==========================================================

scatter_plot <- ggplot(
  scatter_plot_data,
  aes(
    x = X,
    y = Y,
    colour = Regulation
  )
) +
  geom_point(
    size = 2,
    alpha = 0.35,
    shape = 16,
    stroke = 0
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_text_repel(
    data = label_data,
    aes(
      x = X,
      y = Y,
      label = Label
    ),
    inherit.aes = FALSE,
    colour = label_data$label_color,
    size = 5,
    fontface = "bold",
    min.segment.length = 0,
    segment.colour = "black",
    segment.size = 0.6,
    box.padding = 0.9,
    point.padding = 0.35,
    force = 12,
    max.overlaps = Inf,
    seed = 1234
  ) +
  scale_colour_manual(
    values = c(
      "Up" = "#C95A5A",
      "Down" = "#6BAED6"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Ratio (100 \u00b5M Uridine / 30 \u00b5M Uridine)",
    y = "Ratio (100 \u00b5M Uridine / Veh)",
    colour = NULL
  ) +
  theme_classic(
    base_size = 16
  ) +
  theme(
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    legend.position = "right"
  )

print(scatter_plot)

# ==========================================================
# 8. Save figure
# ==========================================================

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_differential_protein_scatter.pdf"
  ),
  plot = scatter_plot,
  width = 6.5,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_differential_protein_scatter.png"
  ),
  plot = scatter_plot,
  width = 6.5,
  height = 5.2,
  units = "in",
  dpi = 600
)

message("Scatter plot analysis completed successfully.")
