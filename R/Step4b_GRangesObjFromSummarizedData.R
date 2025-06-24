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
  # Ensure summarized_dt is a data.table
  if (!data.table::is.data.table(summarized_dt)) {
    summarized_dt <- data.table::as.data.table(summarized_dt)
  }

  # Add 'row.names' column if missing (needed for merge)
  if (!"row.names" %in% colnames(summarized_dt)) {
    summarized_dt[, row.names := rownames(summarized_dt)]
  }

  # Calculate sum of all columns except 'row.names'
  summarized_dt[, sum := rowSums(.SD, na.rm = TRUE), .SDcols = setdiff(names(summarized_dt), "row.names")]

  # Merge with features using 'row.names' as key
  merged <- merge(
    x = features_dt,
    y = summarized_dt[, .(row.names, sum)],
    by.x = "V1",
    by.y = "row.names",
    all.x = FALSE,
    all.y = TRUE
  )

  # Create GRanges object with merged data
  gr <- CreateGRangeObj(merged, feature_id)
  return(gr)
}

