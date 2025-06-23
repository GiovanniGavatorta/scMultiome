#' Summarize Gene Expression and ATAC-seq Data
#'
#' This function performs Step 3 of the exam:
#' It computes the total counts per cell by summing across rows (genes or peaks)
#' in both gene expression and ATAC-seq peak `data.table`s.
#'
#' @param gene_expr_dt A `data.table` containing gene expression data (genes as rows, cells as columns).
#'                     The first column should be gene IDs or row names and is excluded from the sum.
#' @param atac_peaks_dt A `data.table` containing ATAC-seq data (peaks as rows, cells as columns).
#'                      The first column should be peak coordinates and is excluded from the sum.
#'
#' @return A list with two named numeric vectors:
#' \describe{
#'   \item{gene_totals}{Total gene expression counts per cell (numeric vector).}
#'   \item{atac_totals}{Total chromatin accessibility counts per cell (numeric vector).}
#' }
#'
#'
#' @examples
#' \dontrun{
#' # Assuming you've split your data using SplitGeneAndAtac()
#' result <- SummarizeData(gene_expr_dt, atac_peaks_dt)
#' head(result$gene_totals)
#' head(result$atac_totals)
#' }
#' @export
SummarizeData <- function(gene_expr_dt, atac_peaks_dt) {

  gene_counts <- gene_expr_dt[, -1, with = FALSE]
  atac_counts <- atac_peaks_dt[, -1, with = FALSE]

  gene_totals <- colSums(gene_counts)
  atac_totals <- colSums(atac_counts)

  return(list(gene_totals = gene_totals, atac_totals = atac_totals))
}

