###############################################################
## GO MF pie plot based on dose-dependent DIA targets
##
## Input:
##   20240205-DIA-peptides.xlsx
##   sheet = "符合剂量依赖性靶标"
##
## Output:
##   GO_MF_pieplot.pdf
###############################################################

###############################################################
## 0. Clear environment
###############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

###############################################################
## 1. Load packages
###############################################################

library(readxl)

library(dplyr)
library(tidyr)
library(stringr)

library(clusterProfiler)
library(org.Mm.eg.db)

library(ggplot2)
library(ggsci)

library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")

###############################################################
## 2. Set input and output information
###############################################################

input_file <- "20240205-DIA-peptides.xlsx"
input_sheet <- "符合剂量依赖性靶标"

output_pdf <- "GO_MF_pieplot.pdf"

###############################################################
## 3. Read data
###############################################################

data_raw <- read_excel(
  path = input_file,
  sheet = input_sheet
)

data_raw <- data_raw %>%
  as.data.frame() %>%
  as_tibble()

cat("Raw data rows:", nrow(data_raw), "\n")

###############################################################
## 4. Filter dose-dependent significant proteins
###############################################################

data_filt <- data_raw %>%
  dplyr::filter(
    !is.na(pvalue_30),
    !is.na(pvalue_100),
    pvalue_30 < 0.05,
    pvalue_100 < 0.05
  )

cat("After pvalue filter:", nrow(data_filt), "\n")

###############################################################
## 5. Deduplicate Protein.Group
##
## Rule:
##   For duplicated Protein.Group entries, keep the row with
##   the smallest pvalue_100.
###############################################################

data_pg <- data_filt %>%
  dplyr::group_by(Protein.Group) %>%
  dplyr::slice_min(
    order_by = pvalue_100,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

cat("After Protein.Group dedup:", nrow(data_pg), "\n")

###############################################################
## 6. Extract gene symbols
##
## The Genes column may contain multiple gene symbols separated
## by "," or ";". Here, all symbols are split, trimmed, and
## deduplicated.
###############################################################

gene_list <- data_pg %>%
  dplyr::select(Genes) %>%
  dplyr::filter(!is.na(Genes)) %>%
  dplyr::mutate(
    Genes = stringr::str_replace_all(Genes, ",", ";")
  ) %>%
  tidyr::separate_rows(
    Genes,
    sep = ";"
  ) %>%
  dplyr::mutate(
    Genes = stringr::str_trim(Genes)
  ) %>%
  dplyr::filter(Genes != "") %>%
  dplyr::distinct(Genes) %>%
  dplyr::pull(Genes)

cat("Final gene number:", length(gene_list), "\n")

###############################################################
## 7. Convert SYMBOL to ENTREZID
###############################################################

gene_convert <- bitr(
  geneID = gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

if (nrow(gene_convert) == 0) {
  stop("No gene was converted. Please check gene symbol format.")
}

cat("Converted gene number:", nrow(gene_convert), "\n")

###############################################################
## 8. GO enrichment analysis
##
## Ontology:
##   MF: Molecular Function
###############################################################

ego_mf <- enrichGO(
  gene = gene_convert$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_res <- as.data.frame(ego_mf)

cat("GO terms found:", nrow(go_res), "\n")

###############################################################
## 9. Select target GO terms
###############################################################

target_terms <- c(
  "ubiquitin-like protein ligase binding",
  "ubiquitin protein ligase binding",
  "ribonucleoprotein complex binding",
  "ATP hydrolysis activity",
  "unfolded protein binding",
  "protein folding chaperone",
  "mRNA binding",
  "structural constituent of ribosome"
)

ego_data <- go_res %>%
  dplyr::filter(Description %in% target_terms)

cat("Matched GO terms:", nrow(ego_data), "\n")

if (nrow(ego_data) == 0) {
  stop("No target GO terms were matched.")
}

###############################################################
## 10. Calculate GeneRatio-based pie ratio
##
## GeneRatio format:
##   "Count / input gene number"
##
## For pie plot:
##   The selected target terms are normalized internally.
###############################################################

ego_data$Description <- factor(
  ego_data$Description,
  levels = target_terms
)

ego_data$percent_GeneRatio <- sapply(
  ego_data$GeneRatio,
  function(x) {
    ratio_split <- strsplit(as.character(x), "/")[[1]]
    
    if (length(ratio_split) == 2) {
      as.numeric(ratio_split[1]) / as.numeric(ratio_split[2])
    } else {
      NA
    }
  }
)

ego_data <- ego_data %>%
  dplyr::filter(!is.na(percent_GeneRatio))

ratio_sum <- sum(ego_data$percent_GeneRatio)

ego_data <- ego_data %>%
  dplyr::mutate(
    ratio = percent_GeneRatio / ratio_sum
  ) %>%
  dplyr::filter(
    !is.na(ratio),
    ratio > 0
  ) %>%
  dplyr::arrange(Description)

cat("Sum ratio check:", sum(ego_data$ratio), "\n")

###############################################################
## 11. Define pie plot colors
###############################################################

go_colors <- c(
  "ubiquitin-like protein ligase binding" = "#E0DBF0",
  "ubiquitin protein ligase binding" = "#C3BEDB",
  "ribonucleoprotein complex binding" = "#BFD1F2",
  "ATP hydrolysis activity" = "#D1D9E6",
  "unfolded protein binding" = "#F6F4D2",
  "protein folding chaperone" = "#FBE4CB",
  "mRNA binding" = "#F2CDD0",
  "structural constituent of ribosome" = "#D9E9D0"
)

go_colors <- go_colors[target_terms]

###############################################################
## 12. Plot GO MF pie chart
###############################################################

p_go_pie <- ggplot(
  ego_data,
  aes(
    x = "",
    y = ratio,
    fill = Description
  )
) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = "white"
  ) +
  coord_polar("y") +
  scale_fill_manual(
    values = go_colors,
    breaks = target_terms,
    drop = FALSE
  ) +
  geom_text(
    aes(
      label = round(ratio, 2)
    ),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  theme_void() +
  theme(
    legend.title = element_text(
      size = 12,
      color = "black"
    ),
    legend.text = element_text(
      size = 10,
      color = "black"
    ),
    legend.position = "right"
  )

print(p_go_pie)

###############################################################
## 13. Save plot
###############################################################

ggsave(
  filename = output_pdf,
  plot = p_go_pie,
  width = 6,
  height = 4.5,
  device = "pdf"
)

cat("Pie plot saved to:", output_pdf, "\n")
