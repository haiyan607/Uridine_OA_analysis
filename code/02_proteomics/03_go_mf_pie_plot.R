###############################################################
## 0. Load required packages
###############################################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(stringr)
library(dplyr)
library(ggsci)
library(readxl)


###############################################################
## 1. Read the proteomics data
###############################################################

m <- read_xlsx(
  "20240205-DIA-peptides.xlsx",
  sheet = "符合剂量依赖性靶标"
)


###############################################################
## 2. Filter target proteins
###############################################################

# Retain rows in which columns 29 and 31 are both below 0.05
a <- m[
  m[, 29] < 0.05 &
    m[, 31] < 0.05,
]

# Retain rows in which column 30 is above 2
# and column 30 is greater than column 28
b <- m[
  m[, 30] > 2 &
    m[, 30] > m[, 28],
]

# Merge the two filtered datasets
h <- merge(a, b)


###############################################################
## 3. Extract and split gene symbols
###############################################################

gene <- h$Genes

gene <- data.frame(
  str_split(
    gene,
    ";",
    simplify = TRUE
  )
)

gene <- unlist(gene)[
  unlist(gene) != ""
]


###############################################################
## 4. Convert gene symbols to Entrez IDs
###############################################################

gene_EN <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = gene,
  keytype = "SYMBOL",
  columns = c("ENTREZID")
)


###############################################################
## 5. Perform GO enrichment analysis
###############################################################

eGO <- enrichGO(
  gene = gene_EN$ENTREZID,
  keyType = "ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


###############################################################
## 6. Select predefined GO terms
###############################################################

path <- c(
  "ubiquitin-like protein ligase binding",
  "ubiquitin protein ligase binding",
  "ribonucleoprotein complex binding",
  "ATP hydrolysis activity",
  "unfolded protein binding",
  "protein folding chaperone",
  "mRNA binding",
  "structural constituent of ribosome"
)

ego_data <- eGO@result[
  eGO@result$Description %in% path,
]

ego_data$Description <- factor(
  ego_data$Description,
  levels = rev(ego_data$Description)
)


###############################################################
## 7. Generate the GO enrichment pie plot
###############################################################

ggplot(
  ego_data,
  aes(
    x = "Description",
    y = Count,
    fill = Description
  )
) +
  geom_bar(
    width = 1,
    stat = "identity",
    color = "white",
    position = "stack"
  ) +
  coord_polar(
    "y",
    start = 0
  ) +
  scale_fill_npg() +
  theme_void() +
  geom_text(
    aes(
      y = ego_data$Count / 2 +
        c(
          0,
          cumsum(ego_data$Count)[-length(ego_data$Count)]
        ),
      label = round(
        Count / sum(ego_data$Count),
        2
      )
    ),
    size = 3,
    color = "white"
  ) +
  ggtitle("MF_0.05")


###############################################################
## 8. Save the figure and enrichment results
###############################################################

ggsave(
  "pie_GO.pdf",
  width = 922,
  height = 554,
  units = "px",
  dpi = 110
)

write.csv(
  eGO@result,
  "GO_result.csv"
)
