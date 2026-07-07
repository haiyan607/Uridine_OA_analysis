# Raw data

This directory is reserved for raw input data required by the analysis scripts.

Raw sequencing files, large data files and restricted datasets are not stored in this GitHub repository. Their expected locations, filenames and formats are documented below so that the analysis workflow can be reproduced.

## Directory structure

```text
data/raw/
├── README.md
├── GSE220243_RAW/
└── proteomics/
    ├── 20240205-DAP.xlsx
    └── 20240205-DIA-peptides.xlsx
```

The data files and subdirectories shown above may not appear in the public repository because they are excluded by `.gitignore`.

---

## Single-cell RNA-seq data

The single-cell RNA-seq analysis uses the public dataset:

```text
GSE220243
```

The raw sequencing data can be downloaded from the Gene Expression Omnibus.

Because the raw 10X files are large, they are not included in this repository.

After downloading the files, place them in:

```text
data/raw/GSE220243_RAW/
```

The expected filename format is:

```text
GSM6797148_CartNorm1_barcodes.tsv.gz
GSM6797148_CartNorm1_features.tsv.gz
GSM6797148_CartNorm1_matrix.mtx.gz
```

Each sample should contain the following three files:

```text
barcodes.tsv.gz
features.tsv.gz
matrix.mtx.gz
```

The script:

```text
code/01_single_cell/01_create_single_sample_seurat_objects.R
```

reads these raw files and creates the initial single-sample Seurat objects.

---

## Proteomics data

Two proteomics input tables are used for different parts of the analysis.

### Differential-protein Venn and scatter plots

The Venn diagram and differential-protein scatter plot scripts use:

```text
data/raw/proteomics/20240205-DAP.xlsx
```

The corresponding scripts are:

```text
code/02_proteomics/01_differential_protein_venn_plots.R
code/02_proteomics/02_differential_protein_scatter_plot.R
```

The required columns include:

```text
Genes
30uM-RATIO
30PVALUE
100uM-RATIO
100PVALUE
```

The original column names should not be changed before running the scripts.

### GO Molecular Function pie plot

The GO Molecular Function enrichment and pie plot script uses:

```text
data/raw/proteomics/20240205-DIA-peptides.xlsx
```

The corresponding script is:

```text
code/02_proteomics/03_go_mf_pie_plot.R
```

The script reads the following Excel worksheet:

```text
符合剂量依赖性靶标
```

The required columns include:

```text
Protein.Group
Protein.Ids
Genes
30uM-RATIO
pvalue_30
100uM-RATIO
pvalue_100
```

The original worksheet name and column names should be retained because they are referenced directly by the R script.

---

## Data availability

When a proteomics table is provided as a supplementary file associated with the manuscript, it does not need to be uploaded separately to this repository.

Before running the scripts, copy or download the required supplementary data files to:

```text
data/raw/proteomics/
```

The filenames must match those specified above unless the corresponding paths in the R scripts are also updated.

---

## Data-sharing policy

Do not upload raw, confidential, restricted or unpublished data unless public redistribution has been approved.

The following files are normally excluded from the repository:

```text
raw 10X sequencing files
FASTQ files
BAM and SAM files
large Excel input tables
RDS and other processed R objects
restricted or confidential datasets
```

Large input files should remain on the local computer. This README records their expected directory locations and formats so that the analysis can still be reproduced.
