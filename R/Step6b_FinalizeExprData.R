#' Finalize Expression Data by Filtering for Protein-Coding Genes
#'
#' This function performs Step 6 of the exam:
#' takes a GRanges object of expression data and filters it to retain
#' only those entries that correspond to protein-coding genes (based on GTF annotation).
#' It also adds gene symbol names as a metadata column for easier interpretation.
#'
#' @param expr_gr A \code{GRanges} object of expression data, with gene IDs in metadata column \code{Expr_id}.
#' @param protein_coding_gr A \code{GRanges} object of protein-coding genes, with \code{gene_id} and \code{gene_name} in metadata.
#'
#' @return A filtered \code{GRanges} object containing only protein-coding genes,
#' with an added metadata column \code{Gene_Symbol}.
#'
#' @import GenomicRanges
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' expr_filtered <- FinalizeExpressionData(expr_gr, protein_coding_gr)
#' }
#'
#' @export
FinalizeExpressionData <- function(expr_gr, protein_coding_gr) {

  protein_gene_ids <- mcols(protein_coding_gr)$gene_id

  expr_gene_ids <- mcols(expr_gr)$Expr_id

  is_protein_coding <- expr_gene_ids %in% protein_gene_ids

  expr_gr_filtered <- expr_gr[is_protein_coding]

  id_to_symbol <- setNames(mcols(protein_coding_gr)$gene_name, protein_gene_ids)

  mcols(expr_gr_filtered)$Gene_Symbol <- id_to_symbol[mcols(expr_gr_filtered)$Expr_id]

  return(expr_gr_filtered)
}

