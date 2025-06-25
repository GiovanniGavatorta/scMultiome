#' Prepare a GRanges Object from Summarized Data
#'
#' This function performs Step 4b of the exam:
#' processes summarized expression or ATAC-seq data by computing row sums
#' and merging them with genomic coordinates from a features annotation file, before
#' converting them into a `GRanges` object.
#'
#' @param summarized_dt A `data.table` containing summarized data (e.g., gene_expr or atac_peaks).
#'                      The first column should be row names (feature IDs).
#' @param features_dt A `data.table` from the features.tsv.gz file, used to extract
#'                    genomic coordinates (e.g., V1 = ID, V4 = chr, V5 = start, V6 = end).
#' @param feature_id A string used to label metadata fields in the resulting GRanges
#'                   (e.g. "Expr", "ATAC").
#'
#' @return A `GRanges` object with genomic ranges and feature-level metadata
#'
#' @import data.table
#' @examples
#' \dontrun{
#' gr <- PrepareGRanges(summarized_dt = gene_expr, features_dt = features, feature_id = "Expr")
#' }
#' @export
PrepareGRanges <- function(summarized_dt, features_dt, feature_id) {
  summarized_dt$sum <- rowSums(summarized_dt[, -1, with = FALSE], na.rm = TRUE)

  merged <- merge(
    x = features_dt,
    y = summarized_dt[, list(row.names, sum)],
    by.x = "V1",
    by.y = "row.names",
    all.x = FALSE,
    all.y = TRUE
  )

  gr <- CreateGRangeObj(merged, feature_id)
  return(gr)
}

