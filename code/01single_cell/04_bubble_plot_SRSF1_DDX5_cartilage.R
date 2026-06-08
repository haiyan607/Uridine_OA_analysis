# ==========================================================
# 04_bubble_plot_SRSF1_DDX5_cartilage.R
# Bubble plot: SRSF1 and DDX5 in Cartilage (Normal vs OA)
# ==========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(plot1cell)
library(ggplot2)
library(grid)
library(scales)
library(here)

# ==========================================================
# 1. Paths
# ==========================================================
processed_dir <- here("data", "processed")
figures_dir <- here("results", "figures")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(processed_dir, "GSE220243_QC_normalized_hvg_scaled_harmony_umap.rds")
if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

# ==========================================================
# 2. Read object
# ==========================================================
sce <- readRDS(input_file)

if (!"umap" %in% Reductions(sce)) {
  stop("Current object does not contain UMAP. Please run script 03 first.")
}

meta_cols <- colnames(sce@meta.data)

if (!"tissue" %in% meta_cols) {
  if (!"orig.ident" %in% meta_cols) {
    stop("Neither tissue nor orig.ident found in metadata.")
  }
  sce$tissue <- ifelse(
    grepl("Cart", sce$orig.ident, ignore.case = TRUE), "Cartilage",
    ifelse(grepl("Men", sce$orig.ident, ignore.case = TRUE), "Meniscus", NA)
  )
}

if (!"group" %in% meta_cols) {
  if (!"orig.ident" %in% meta_cols) {
    stop("Neither group nor orig.ident found in metadata.")
  }
  sce$group <- ifelse(
    grepl("Norm", sce$orig.ident, ignore.case = TRUE), "Normal",
    ifelse(grepl("OA", sce$orig.ident, ignore.case = TRUE), "OA", NA)
  )
}

# ==========================================================
# 3. Keep Cartilage only
# ==========================================================
sce_cart <- subset(
  sce,
  subset = tissue == "Cartilage" & !is.na(group)
)

sce_cart$Group <- factor(sce_cart$group, levels = c("Normal", "OA"))
Idents(sce_cart) <- "tissue"

# ==========================================================
# 4. Check genes
# ==========================================================
genes_target <- c("SRSF1", "DDX5")

genes_found <- rownames(sce_cart)[
  match(toupper(genes_target), toupper(rownames(sce_cart)))
]

if (any(is.na(genes_found))) {
  missing_genes <- genes_target[is.na(genes_found)]
  stop("Missing genes: ", paste(missing_genes, collapse = ", "))
}

names(genes_found) <- genes_target

# ==========================================================
# 5. plot1cell dotplot
# ==========================================================
complex_dotplot_multiple(
  seu_obj = sce_cart,
  features = unname(genes_found),
  group = "Group",
  celltypes = "Cartilage"
)

p <- ggplot2::last_plot()

# ==========================================================
# 6. Helper function
# ==========================================================
fix_avgexp_scale_keep_palette <- function(
    p,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1)
) {
  idx_colour <- which(vapply(p$scales$scales, function(s) {
    any(s$aesthetics %in% c("colour", "color"))
  }, logical(1)))
  
  idx_fill <- which(vapply(p$scales$scales, function(s) {
    any(s$aesthetics %in% c("fill"))
  }, logical(1)))
  
  idx <- if (length(idx_colour) > 0) {
    idx_colour[1]
  } else if (length(idx_fill) > 0) {
    idx_fill[1]
  } else {
    NA_integer_
  }
  
  if (is.na(idx)) stop("Cannot find colour/fill scale.")
  
  sc <- p$scales$scales[[idx]]
  sc$limits <- limits
  sc$breaks <- breaks
  sc$oob <- scales::squish
  p$scales$scales[[idx]] <- sc
  p
}

p2 <- fix_avgexp_scale_keep_palette(
  p,
  limits = c(-1, 1),
  breaks = c(-1, -0.5, 0, 0.5, 1)
)

# ==========================================================
# 7. Beautify
# ==========================================================
p3 <- p2 +
  theme(
    strip.background = element_rect(fill = "#F6B0AA", colour = NA),
    strip.text = element_text(colour = "black", face = "bold", size = 18),
    axis.text.x = element_text(
      angle = 0, hjust = 0.5, vjust = 0.5,
      colour = "black", face = "bold", size = 16
    ),
    axis.text.y = element_text(
      colour = "black", face = "bold", size = 16
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  guides(
    size = guide_legend(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(shape = 16, colour = "black")
    ),
    colour = guide_colorbar(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.75, "cm"),
      barheight = unit(2.4, "cm")
    ),
    fill = guide_colorbar(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.75, "cm"),
      barheight = unit(2.4, "cm")
    )
  )

# ==========================================================
# 8. Save
# ==========================================================
ggsave(
  filename = file.path(figures_dir, "Figure_SRSF1_DDX5_Cartilage_Normal_OA.pdf"),
  plot = p3,
  width = 8,
  height = 6,
  units = "in"
)

ggsave(
  filename = file.path(figures_dir, "Figure_SRSF1_DDX5_Cartilage_Normal_OA.png"),
  plot = p3,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

message("Finished: SRSF1/DDX5 bubble plot.")