#' Split Gene Expression and ATAC-seq Data
#'
#' This function performs Step 2 of the exam:
#' It splits the full `data.table` from Step 1 into two separate tables:
#' one for gene expression (rows starting with "ENSG") and one for ATAC-seq peaks
#' (rows starting with "chr").
#'
#' @param dt A `data.table` created by `ConvertMatrix()` with a column called `"row.names"`
#'          containing either Ensembl gene IDs or genomic peak coordinates.
#'
#' @return A list with two named elements:
#' \describe{
#'   \item{gene_expr}{A `data.table` of gene expression rows (Ensembl gene IDs).}
#'   \item{atac_peaks}{A `data.table` of ATAC-seq peak rows (chromosome coordinates).}
#' }
#'
#'
#' @examples
#' \dontrun{
#' # Assuming expr_dt is the output from ConvertMatrix()
#' split_data <- SplitGeneAndAtac(expr_dt)
#' gene_expr_dt <- split_data$gene_expr
#' atac_peaks_dt <- split_data$atac_peaks
#'
#' head(gene_expr_dt)
#' head(atac_peaks_dt)
#' }
#' @export
SplitGeneAndAtac <- function(dt) {

  row_ids <- dt[[1]]

  gene <- grepl("^ENSG", row_ids)
  peak <- grepl("^chr", row_ids)

  gene_expr <- dt[gene, ]
  atac_peaks <- dt[peak, ]

  return(list(gene_expr = gene_expr, atac_peaks = atac_peaks))
}

