# Uridine-Related Osteoarthritis Analysis

This repository contains the R scripts used to reproduce selected bioinformatics analyses and figures in the associated manuscript. The current version includes single-cell RNA-sequencing analysis of GSE220243 and proteomics-based Venn and scatter plots.

Additional scripts, including the GO functional-classification pie chart, may be added in later updates.

## Repository structure

```text
Uridine_OA_analysis/
├── README.md
├── LICENSE
├── .gitignore
├── code/
│   ├── 01_single_cell/
│   │   ├── 01_create_single_sample_seurat_objects.R
│   │   ├── 02_quality_control_and_normalization.R
│   │   ├── 03_run_pca_harmony_umap.R
│   │   ├── 04_bubble_plot_SRSF1_DDX5_cartilage.R
│   │   └── 05_bubble_plot_SLC16A1_SLC28_family.R
│   └── 02_proteomics/
│       ├── 01_differential_protein_venn_plots.R
│       └── 02_differential_protein_scatter_plot.R
├── data/
│   ├── raw/
│   │   └── README.md
│   ├── metadata/
│   │   └── GSE220243_sample_metadata.csv
│   └── processed/
│       └── README.md
├── results/
│   ├── figures/
│   └── tables/
└── docs/
    └── sessionInfo.txt
```

## Analysis overview

### 1. Single-cell RNA sequencing

The single-cell RNA-seq analysis is based on the public dataset **GSE220243**. The dataset contains 26 sequencing samples:

- 6 normal cartilage samples
- 6 osteoarthritis cartilage samples
- 8 normal meniscus sequencing samples from 7 donors
- 6 osteoarthritis meniscus samples

The workflow includes:

1. organization and import of raw 10X files;
2. creation of single-sample Seurat objects;
3. calculation of mitochondrial and ribosomal proportions;
4. sample-wise quality-control filtering;
5. merging and normalization;
6. highly variable feature selection and scaling;
7. PCA, Harmony batch correction and UMAP;
8. generation of manuscript bubble plots.

The target-gene bubble plots include:

- **SRSF1 and DDX5** in cartilage, comparing Normal and OA groups;
- **SLC16A1, SLC28A1, SLC28A2 and SLC28A3** in cartilage and meniscus, comparing Normal and OA groups.

### 2. Proteomics

The proteomics scripts require a table containing the following columns:

```text
Genes
30uM-RATIO
30PVALUE
100uM-RATIO
100PVALUE
```

The scripts generate:

- up-regulated and down-regulated four-set Venn diagrams;
- four-way intersection protein lists;
- a differential-protein scatter plot;
- coordinates and labels for selected proteins.

For duplicated protein symbols, the row with the largest absolute fold change across the 30 μM and 100 μM conditions is retained.

## R scripts

### Single-cell RNA-seq

#### `01_create_single_sample_seurat_objects.R`

Reads the raw 10X matrices for all 26 samples, creates individual Seurat objects, adds sample metadata and saves the raw object list.

#### `02_quality_control_and_normalization.R`

Calculates quality-control metrics, performs sample-wise filtering, merges all samples, normalizes expression values, identifies highly variable features and scales the data.

Main quality-control criteria:

```text
500 < nFeature_RNA < 6000
1000 < nCount_RNA < sample-specific 99th percentile
Cartilage percent.mt < min(sample-specific 98th percentile, 15%)
Meniscus percent.mt < min(sample-specific 98th percentile, 20%)
```

`percent.Ribo` is calculated for inspection but is not used for filtering.

#### `03_run_pca_harmony_umap.R`

Runs PCA, Harmony batch correction and UMAP. Batch correction is performed using `sample_id`.

#### `04_bubble_plot_SRSF1_DDX5_cartilage.R`

Generates the SRSF1/DDX5 bubble plot for cartilage, comparing Normal and OA groups.

#### `05_bubble_plot_SLC16A1_SLC28_family.R`

Generates the transporter-family bubble plot for cartilage and meniscus, comparing Normal and OA groups.

### Proteomics

#### `01_differential_protein_venn_plots.R`

Generates separate up-regulated and down-regulated four-set Venn diagrams.

Up-regulated sets:

```text
30PVALUE < 0.05
100uM-RATIO > 2
100uM-RATIO / 30uM-RATIO > 1
30uM-RATIO > 1
```

Down-regulated sets:

```text
100PVALUE < 0.05
100uM-RATIO < 1
100uM-RATIO / 30uM-RATIO < 1
30uM-RATIO < 1
```

#### `02_differential_protein_scatter_plot.R`

Generates a scatter plot using:

```text
X = 100uM-RATIO / 30uM-RATIO
Y = 100uM-RATIO
```

Proteins included in the scatter analysis satisfy:

```text
30PVALUE < 0.05
100PVALUE < 0.05
```

Classification criteria:

```text
Up:   Y > 2 and X > 1
Down: Y < 0.5 and X < 1
```

Selected proteins labelled in the figure include:

```text
SRSF1
DDX5
FUS
RBM39
SLC16A1
IGF2R
```

## Data availability and repository policy

Raw single-cell sequencing files and large processed Seurat objects are not included in this repository because of their size. The single-cell dataset can be obtained from GEO under accession **GSE220243**.

Data tables that are already provided as manuscript supplementary files do not need to be uploaded again to GitHub. The repository mainly provides analysis code, metadata, software information and small derived result tables required to understand or verify the workflow.

A supplementary data file may also be included in GitHub when all of the following apply:

- public redistribution is permitted;
- the file does not contain confidential or unpublished information;
- the file is small enough for normal GitHub storage;
- including it materially improves code reproducibility.

When a required input table is not included, its expected columns and directory location are documented in the repository.

## Data preparation

### Single-cell data

Download the GSE220243 raw 10X files and place them in:

```text
data/raw/GSE220243_RAW/
```

The first script supports files named in the following format:

```text
GSM6797148_CartNorm1_barcodes.tsv.gz
GSM6797148_CartNorm1_features.tsv.gz
GSM6797148_CartNorm1_matrix.mtx.gz
```

### Proteomics data

When redistribution is permitted, place the proteomics input table at:

```text
data/raw/proteomics/20240205-DAP.xlsx
```

If the table is supplied only as a manuscript supplementary file, download or copy it to the location above before running the proteomics scripts.

## Running the analysis

Run the scripts from the repository root in the following order:

```text
code/01_single_cell/01_create_single_sample_seurat_objects.R
code/01_single_cell/02_quality_control_and_normalization.R
code/01_single_cell/03_run_pca_harmony_umap.R
code/01_single_cell/04_bubble_plot_SRSF1_DDX5_cartilage.R
code/01_single_cell/05_bubble_plot_SLC16A1_SLC28_family.R

code/02_proteomics/01_differential_protein_venn_plots.R
code/02_proteomics/02_differential_protein_scatter_plot.R
```

The single-cell and proteomics workflows are independent and can be run separately.

## Software requirements

The scripts use the following R packages:

```r
Seurat
harmony
plot1cell
readxl
dplyr
tidyr
stringr
Matrix
ggplot2
ggrepel
ggvenn
patchwork
scales
gtable
grid
here
```

The exact R and package versions used for the analysis should be recorded in:

```text
docs/sessionInfo.txt
```

## Outputs

Figures are written to:

```text
results/figures/
```

Summary tables and protein lists are written to:

```text
results/tables/
```

Large intermediate Seurat objects are written to:

```text
data/processed/
```

Large RDS files and raw sequencing data should be excluded from GitHub using `.gitignore`.

## Reproducibility notes

- Run all scripts from the repository root.
- The scripts use `here()` for relative paths.
- Keep the original input column names unchanged.
- UMAP and other stochastic procedures use a fixed random seed where applicable.
- Differences in Seurat, SeuratObject, Harmony or plot1cell versions may affect object structure or visual output.

## Citation

The manuscript citation will be added after publication.

The public single-cell RNA-seq dataset is available from GEO under accession **GSE220243**.

## Contact

GitHub: [haiyan607](https://github.com/haiyan607)

Email: zhy19840108919@163.com
