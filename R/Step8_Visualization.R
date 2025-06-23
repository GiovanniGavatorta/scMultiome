#' Plot Expression vs ATAC-seq CPM
#'
#' This function performs Step 7 of the exam:
#' generates scatter plots to visualize the relationship between
#' log2-transformed CPM-normalized gene expression and ATAC-seq signals.
#' If the dataset is too large, it creates separate plots for each chromosome.
#'
#' @param merged_data A \code{data.table} or \code{data.frame} containing at least:
#'   \describe{
#'     \item{gene_id}{Gene identifier}
#'     \item{expr_cpm}{Log2 CPM-normalized gene expression values}
#'     \item{atac_cpm}{Log2 CPM-normalized ATAC-seq values}
#'   }
#' @param expr_gr A \code{GRanges} object for expression data, used to extract chromosome information.
#' @param max_points Integer; maximum number of points to plot in a single plot before splitting by chromosome (default 5000).
#'
#' @return Either a \code{ggplot2} object (single scatter plot) or a list of \code{ggplot2} plots (one per chromosome)
#'
#' @import ggplot2
#' @import data.table
#'
#' @examples
#' \dontrun{
#' plot_all <- PlotExprVsATAC(merged_data, expr_gr)
#' print(plot_all)
#' }
#' @export
PlotExprVsATAC <- function(merged_data, expr_gr, max_points = 5000) {

  chr_info <- data.table(
    gene_id = mcols(expr_gr)$Expr_id,
    chr = as.character(seqnames(expr_gr))
  )

  plot_data <- merge(as.data.table(merged_data), chr_info, by = "gene_id", all.x = TRUE)

  if (nrow(plot_data) > max_points) {
    plots_list <- lapply(unique(plot_data$chr), function(chrom) {
      p <- ggplot(plot_data[chr == chrom], aes(x = expr_cpm, y = atac_cpm)) +
        geom_point(alpha = 0.4, size = 1) +
        theme_minimal() +
        labs(
          title = paste("Expression vs ATAC CPM - Chromosome", chrom),
          x = "log2(CPM Expression + 1)",
          y = "log2(CPM ATAC + 1)"
        ) +
        xlim(min(plot_data$expr_cpm, na.rm = TRUE), max(plot_data$expr_cpm, na.rm = TRUE)) +
        ylim(min(plot_data$atac_cpm, na.rm = TRUE), max(plot_data$atac_cpm, na.rm = TRUE))
      return(p)
    })

    return(plots_list)

  } else {
    p <- ggplot(plot_data, aes(x = expr_cpm, y = atac_cpm)) +
      geom_point(alpha = 0.4, size = 1) +
      theme_minimal() +
      labs(
        title = "Expression vs ATAC CPM - All Genes",
        x = "log2(CPM Expression + 1)",
        y = "log2(CPM ATAC + 1)"
      )
    return(p)
  }
}

