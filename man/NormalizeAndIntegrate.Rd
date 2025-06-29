% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Step7_DataNormalizationAndIntegration.R
\name{NormalizeAndIntegrate}
\alias{NormalizeAndIntegrate}
\title{Normalize and Integrate Expression & ATAC-seq Data}
\usage{
NormalizeAndIntegrate(expr_gr, atac_annotated_gr)
}
\arguments{
\item{expr_gr}{A \code{GRanges} object containing gene expression metadata, including \code{Expr_id} and \code{Expr_sum}.}

\item{atac_annotated_gr}{A \code{GRanges} object of ATAC-seq peaks annotated with overlapping gene IDs,
containing \code{ATAC_sum} and \code{overlapping_gene_id} metadata columns.}
}
\value{
A list containing:
\describe{
  \item{merged_table}{A \code{data.table} with CPM-normalized expression and ATAC counts for shared gene IDs}
  \item{atac_not_merged}{\code{data.table} of ATAC peaks with no matched expression gene}
  \item{expr_not_merged}{\code{data.table} of expression genes with no overlapping ATAC peak}
  \item{plot_atac_by_chr}{A \code{ggplot2} boxplot of ATAC-seq intensity by chromosome}
  \item{plot_expr_by_chr}{A \code{ggplot2} boxplot of expression intensity (genes without ATAC) by chromosome}
}
}
\description{
This function performs Step 7 of the exam:
performs normalization of gene expression and ATAC-seq counts using CPM (Counts Per Million),
merges them based on shared gene IDs, and generates summary boxplots for quality control by chromosome.
}
\examples{
\dontrun{
results <- NormalizeAndIntegrate(expr_gr, atac_annotated_gr)
print(results$plot_atac_by_chr)
print(results$plot_expr_by_chr)
head(results$merged_table)
}
}
