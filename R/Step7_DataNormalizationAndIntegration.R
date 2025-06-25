#' Normalize and Integrate Expression & ATAC-seq Data
#'
#' This function performs Step 7 of the exam:
#' performs normalization of gene expression and ATAC-seq counts using CPM (Counts Per Million),
#' merges them based on shared gene IDs, and generates summary boxplots for quality control by chromosome.
#'
#' @param expr_gr A \code{GRanges} object containing gene expression metadata, including \code{Expr_id} and \code{Expr_sum}.
#' @param atac_annotated_gr A \code{GRanges} object of ATAC-seq peaks annotated with overlapping gene IDs,
#'                          containing \code{ATAC_sum} and \code{overlapping_gene_id} metadata columns.
#'
#' @return A list containing:
#' \describe{
#'   \item{merged_table}{A \code{data.table} with CPM-normalized expression and ATAC counts for shared gene IDs}
#'   \item{atac_not_merged}{\code{data.table} of ATAC peaks with no matched expression gene}
#'   \item{expr_not_merged}{\code{data.table} of expression genes with no overlapping ATAC peak}
#'   \item{plot_atac_by_chr}{A \code{ggplot2} boxplot of ATAC-seq intensity by chromosome}
#'   \item{plot_expr_by_chr}{A \code{ggplot2} boxplot of expression intensity (genes without ATAC) by chromosome}
#' }
#'
#' @import data.table
#' @import GenomicRanges
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' results <- NormalizeAndIntegrate(expr_gr, atac_annotated_gr)
#' print(results$plot_atac_by_chr)
#' print(results$plot_expr_by_chr)
#' head(results$merged_table)
#' }
#' @export
NormalizeAndIntegrate <- function(expr_gr, atac_annotated_gr) {

  expr_counts <- mcols(expr_gr)$Expr_sum
  atac_counts <- mcols(atac_annotated_gr)$ATAC_sum

  names(expr_counts) <- mcols(expr_gr)$Expr_id
  names(atac_counts) <- mcols(atac_annotated_gr)$overlapping_gene_id

  NormalizeCPM <- function(counts) {
    cpm <- (counts / sum(counts)) * 1e6
    log2(cpm + 1)
  }

  expr_cpm <- NormalizeCPM(expr_counts)
  atac_cpm <- NormalizeCPM(atac_counts)

  common_genes <- intersect(names(expr_cpm), names(atac_cpm))

  merged_dt <- data.table(
    gene_id = common_genes,
    expr_cpm = expr_cpm[common_genes],
    atac_cpm = atac_cpm[common_genes]
  )

  atac_only_genes <- setdiff(names(atac_cpm), names(expr_cpm))
  atac_only_dt <- data.table(gene_id = atac_only_genes, atac_cpm = atac_cpm[atac_only_genes])

  expr_only_genes <- setdiff(names(expr_cpm), names(atac_cpm))
  expr_only_dt <- data.table(gene_id = expr_only_genes, expr_cpm = expr_cpm[expr_only_genes])

  atac_chr <- data.table(
    gene_id = mcols(atac_annotated_gr)$overlapping_gene_id,
    chr = as.character(seqnames(atac_annotated_gr)),
    atac_cpm = atac_cpm[mcols(atac_annotated_gr)$overlapping_gene_id]
  )

  expr_chr <- data.table(
    gene_id = mcols(expr_gr)$Expr_id,
    chr = as.character(seqnames(expr_gr)),
    expr_cpm = expr_cpm[mcols(expr_gr)$Expr_id]
  )

  expr_chr_filtered <- expr_chr[expr_chr[["gene_id"]] %in% expr_only_genes]

  p1 <- ggplot(atac_chr, aes(x = chr, y = atac_cpm)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "ATAC-seq Intensity by Chromosome", y = "log2(CPM+1)", x = "Chromosome")

  p2 <- ggplot(expr_chr_filtered, aes(x = chr, y = expr_cpm)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Gene Expression (no ATAC peaks) by Chromosome", y = "log2(CPM+1)", x = "Chromosome")

  return(list(
    merged_table = merged_dt,
    atac_not_merged = atac_only_dt,
    expr_not_merged = expr_only_dt,
    plot_atac_by_chr = p1,
    plot_expr_by_chr = p2
  ))
}

