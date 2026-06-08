# Raw data

This directory is reserved for raw input data required by the analysis scripts.

## Single-cell RNA-seq data

The single-cell RNA-seq dataset used in this project is available from the Gene Expression Omnibus under accession:

```text
GSE220243
```

Because the raw 10X files are large, they are not included in this GitHub repository.

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

The first single-cell script can automatically organize these files into standard sample-specific 10X directories.

## Proteomics data

The proteomics scripts expect the input file at:

```text
data/raw/proteomics/20240205-DAP.xlsx
```

The required columns are:

```text
Genes
30uM-RATIO
30PVALUE
100uM-RATIO
100PVALUE
```

If this dataset is already provided as a manuscript supplementary file, it does not need to be uploaded again to GitHub. Copy or download the supplementary file to the path above before running the proteomics scripts.

## Data-sharing policy

Do not upload raw or restricted data unless public redistribution is permitted.

Large files, raw sequencing data and confidential datasets should remain outside the repository. Their expected paths and formats are documented here so that the analysis can still be reproduced.
