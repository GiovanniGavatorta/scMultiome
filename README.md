# scMultiome

**scMultiome** is an R package designed for the analysis and visualization of single-cell multi-omics data. It provides tools to integrate and interpret multiple layers of genomic information from single-cell experiments.

## Features

-   Integration of transcriptomic and epigenomic single-cell data
-   Visualization functions for multi-modal datasets
-   Support for large-scale datasets with efficient data handling
-   Built-in example datasets for quick start

## Installation

Currently, **scMultiome** is available on GitHub. You can install it using `devtools`:

``` r
# Install devtools if you don't have it
install.packages("devtools")

# Install scMultiome from GitHub
devtools::install_github("GiovanniGavatorta/scMultiome")
```

## Example Workflow

This example demonstrates a complete analysis pipeline using the core functions of the `scMoltiOme` package.

``` r
# Load the package and other libraries required
library(scMultiome)
library(Matrix)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tidyr)
library(IRanges)
library(rtracklayer)
```

### 1. Data Loading and Matrix Conversion

``` r
expr_dt <- ConvertMatrix(
  mtx_path = (file = "path/to/a/file/named/filename.mtx.gz"),
  features_path = (file = "path/to/a/file/of/features/named/filenamefeatures.tsv.gz"),
  barcodes_path = (file = "path/to/a/file/of/barcodes/named/filenamebarcodes.tsv.gz")
)
```

### 2. Data Extraction and Annotation

``` r
# Step 2: Split Gene Expression and ATAC-seq Data
split_data <- SplitGeneAndAtac(expr_dt)
gene_expr_dt <- split_data$gene_expr
atac_peaks_dt <- split_data$atac_peaks

# Step3 3: Summarize Data
summary_results <- SummarizeData(gene_expr_dt, atac_peaks_dt)
gene_totals <- summary_results$gene_totals
atac_totals <- summary_results$atac_totals

# Step 4: Create Genomic Ranges
features <- fread(file = "path/to/a/file/of/features/named/filenamefeatures.tsv.gz", header = FALSE)
gene_expr_dt$sum <- rowSums(gene_expr_dt[, -1, with = FALSE])
atac_peaks_dt$sum <- rowSums(atac_peaks_dt[, -1, with = FALSE])
gr_expr <- PrepareGRanges(gene_expr_dt, features_dt = features, feature_id = "Expr")
gr_atac <- PrepareGRanges(atac_peaks_dt, features_dt = features, feature_id = "ATAC")

# Step 5: Gene Annotation for ATACseq data
gtf_path <- import("path/to/a/file/named/Homo_sapiens.GRCh38.114.gtf.gz")
annotated_atac <- AnnotateATACwithGenes(gtf_path, gr_atac)

# Step 6: Finalize Expression Data
gtf_path <- import("path/to/a/file/named/Homo_sapiens.GRCh38.114.gtf.gz")
protein_coding_gr <- CreateProteinCodingGR(gtf_path)
expr_filtered <- FinalizeExpressionData(gr_expr, protein_coding_gr)
```

### 3. Data Normalization, Integration, and Final Visualization

``` r
results <- NormalizeAndIntegrate(expr_filtered, annotated_atac)
merged_data <- results$merged_table
plot_result <- PlotExprVsATAC(merged_data, gr_expr)

# Display plots
# if (is.list(plot_result)) {
# Show first plot in the list, regardless of chromosome name
#  first_plot <- plot_result[[1]]
#  first_plot
#} else {
#  plot_result
#}
```

## Detailed Tutorial

If you need a more detailed step-by-step guide on how to use the package, please see the package vignette.

``` r
# This will open the detailed tutorial in your browser
vignette("HowTo_scMultiome", package = "scMultiome")
```

## Function Overview

| Function | Description |
|---------------------|---------------------------------------------------|
| `ConvertMatrix` | converts the sparse matrix into a full matrix and save the result as a `data.table` object, which is the input for the rest of the pipeline. |
| `SplitGeneAndAtac` | splits the full `data.table` created in the step before into two separate tables: one for gene expression (rows starting with "ENSG") and one for ATAC-seq peaks (rows starting with "chr"). |
| `SummarizeData` | compute the total counts per cell by summing across rows (genes or peaks) in both gene expression and ATAC-seq peak `data.table`s. |
| `PrepareGRanges` | processes summarized expression or ATAC-seq data by computing row sums and merging them with genomic coordinates from a features annotation file, before converting them into a `GRanges` object. |
| `CreateGRangeObj` | converts a data.frame of genomic features and metadata into a `GRanges` object. It is typically used to convert ATAC-seq peaks or gene expression features into genomic ranges. |
| `AnnotateATACwithGenes` | reads a GTF annotation file, filters for protein-coding genes, and then finds overlaps between ATAC-seq peaks and these genes. It returns a `GRanges` object of ATAC peaks annotated with gene information. |
| `CreateProteinCodingGR` | imports a GTF annotation file and filters it to keep only features of type "gene" that are protein-coding. The output is a GRanges object ready for downstream overlap analysis. |
| `FinalizeExpressionData` | take a GRanges object of expression data and filters it to retain only those entries that correspond to protein-coding genes (based on GTF annotation). |
| `NormalizeAndIntegrate` | performs normalization of gene expression and ATAC-seq counts using CPM (Counts Per Million), merges them based on shared gene IDs, and generates summary boxplots for quality control by chromosome. |
| `PlotExprVsATAC` | generates scatter plots to visualize the relationship between log2-transformed CPM-normalized gene expression and ATAC-seq signals. If the dataset is too large, it creates separate plots for each chromosome. |


# **Docker Analysis Workflow**

To execute the `scMoltiOme` analysis pipeline in a reproducible environment, you can take advantage of a Docker image hosted on Docker Hub. Using this containerized approach guarantees that all necessary software dependencies are pre-installed, making it easy to run the pipeline on any machine with Docker support—ideal for ensuring consistency in genomics workflows. 

### Requirements:

1. **Docker:** Make sure Docker is installed and operational on your machine.
2. **Project setup:** Prepare a working directory on your local system, structured as follows:

- Create a `data` folder where you will place:

- Your `filtered_feature_barcodes_matrix` output from Cell Ranger.
- Your reference GTF file (`Homo_sapiens.GRCh38.114.gtf.gz`).

- Create an empty `outputs` folder where the analysis results will be saved.

Your folder layout should resemble this:

```
my_project/
├── data/
│   ├── filtered_feature_barcodes_matrix/
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── Homo_sapiens.GRCh38.114.gtf.gz
└── outputs/
```

### How to run the pipeline:
From your terminal, navigate to the main directory (`my_project`) and run the command below. It will pull the `ggavatorta/examproject` Docker image, execute the R-based analysis inside the container, and save results to your local outputs directory.

```
docker run --rm \
  -v "path/to/the/data/folder":/data \
  -v "path/to/the/outputs/folder":/results \
  ggavatorta/examproject:gavatorta
```

Once the container completes, you'll find all `output` plots and result files inside the outputs folder.
