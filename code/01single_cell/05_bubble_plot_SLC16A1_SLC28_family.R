# ==========================================================
# 05_bubble_plot_SLC16A1_SLC28_family.R
# Bubble plot: SLC16A1 / SLC28A1 / SLC28A2 / SLC28A3
# in Cartilage and Meniscus, split by Normal and OA
# ==========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(plot1cell)
library(ggplot2)
library(grid)
library(scales)
library(gtable)
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

if (!"tissue" %in% colnames(sce@meta.data)) {
  if (!"orig.ident" %in% colnames(sce@meta.data)) {
    stop("Cannot identify Cartilage/Meniscus.")
  }
  sce$tissue <- ifelse(
    grepl("Cart", sce$orig.ident, ignore.case = TRUE),
    "Cartilage",
    ifelse(grepl("Men", sce$orig.ident, ignore.case = TRUE), "Meniscus", NA)
  )
}

if (!"group" %in% colnames(sce@meta.data)) {
  if (!"orig.ident" %in% colnames(sce@meta.data)) {
    stop("Cannot identify Normal/OA.")
  }
  sce$group <- ifelse(
    grepl("Norm", sce$orig.ident, ignore.case = TRUE),
    "Normal",
    ifelse(grepl("OA", sce$orig.ident, ignore.case = TRUE), "OA", NA)
  )
}

# ==========================================================
# 3. Subset data
# ==========================================================
sce_use <- subset(
  sce,
  subset = tissue %in% c("Cartilage", "Meniscus") &
    group %in% c("Normal", "OA")
)

sce_use$tissue <- factor(sce_use$tissue, levels = c("Cartilage", "Meniscus"))
sce_use$Group <- factor(sce_use$group, levels = c("Normal", "OA"))
Idents(sce_use) <- "tissue"

# ==========================================================
# 4. Check genes
# ==========================================================
genes_target <- c("SLC16A1", "SLC28A1", "SLC28A2", "SLC28A3")

genes_found <- rownames(sce_use)[
  match(toupper(genes_target), toupper(rownames(sce_use)))
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
  seu_obj = sce_use,
  features = unname(genes_found),
  group = "Group",
  celltypes = c("Cartilage", "Meniscus")
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
  idx_colour <- which(vapply(
    p$scales$scales,
    function(s) any(s$aesthetics %in% c("colour", "color")),
    logical(1)
  ))
  
  idx_fill <- which(vapply(
    p$scales$scales,
    function(s) any(s$aesthetics %in% c("fill")),
    logical(1)
  ))
  
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
# 7. Base theme
# ==========================================================
p3 <- p2 +
  theme(
    strip.background = element_rect(fill = "white", colour = NA, linewidth = 0),
    strip.text = element_text(colour = "black", face = "bold", size = 18),
    axis.text.x = element_text(
      angle = 0, hjust = 0.5, vjust = 0.5,
      colour = "black", face = "bold", size = 16
    ),
    axis.text.y = element_text(
      colour = "black", face = "italic", size = 16
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    legend.box = "horizontal",
    legend.direction = "vertical",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.spacing.x = unit(0.8, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  guides(
    colour = guide_colorbar(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.8, "cm"),
      barheight = unit(2.5, "cm")
    ),
    fill = guide_colorbar(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.8, "cm"),
      barheight = unit(2.5, "cm")
    ),
    size = guide_legend(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(shape = 16, colour = "black")
    )
  )

# ==========================================================
# 8. Modify strip colors
# ==========================================================
g <- ggplotGrob(p3)
strip_ids <- grep("^strip-t", g$layout$name)

cartilage_fill <- "#F97D47"
meniscus_fill  <- "#D6C1DA"

for (i in seq_along(strip_ids)) {
  strip_grob <- g$grobs[[strip_ids[i]]]
  
  rect_index <- which(vapply(
    strip_grob$grobs[[1]]$children,
    function(x) inherits(x, "rect"),
    logical(1)
  ))
  
  if (length(rect_index) > 0) {
    if (i == 1) {
      strip_grob$grobs[[1]]$children[[rect_index]]$gp$fill <- cartilage_fill
      strip_grob$grobs[[1]]$children[[rect_index]]$gp$col  <- NA
      strip_grob$grobs[[1]]$children[[rect_index]]$gp$lwd  <- 0
    } else if (i == 2) {
      strip_grob$grobs[[1]]$children[[rect_index]]$gp$fill <- meniscus_fill
      strip_grob$grobs[[1]]$children[[rect_index]]$gp$col  <- NA
      strip_grob$grobs[[1]]$children[[rect_index]]$gp$lwd  <- 0
    }
  }
  
  g$grobs[[strip_ids[i]]] <- strip_grob
}

# ==========================================================
# 9. Save
# ==========================================================
pdf(
  file = file.path(figures_dir, "Figure_SLC16A1_SLC28A1_SLC28A2_SLC28A3.pdf"),
  width = 9.5,
  height = 4.8
)
grid.newpage()
grid.draw(g)
dev.off()

png(
  filename = file.path(figures_dir, "Figure_SLC16A1_SLC28A1_SLC28A2_SLC28A3.png"),
  width = 9.5,
  height = 4.8,
  units = "in",
  res = 600
)
grid.newpage()
grid.draw(g)
dev.off()

message("Finished: transporter-family bubble plot.")