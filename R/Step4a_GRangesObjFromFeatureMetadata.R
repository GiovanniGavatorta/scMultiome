#' Create a GRanges Object from Feature Metadata
#'
#' This function performs Step 4a of the exam:
#' converts a data.frame of genomic features and metadata into a `GRanges` object.
#' It is typically used to convert ATAC-seq peaks or gene expression features into genomic ranges.
#'
#' @param specific_feature A `data.frame` that includes at least these columns:
#'   \itemize{
#'     \item \code{V1}: feature ID (e.g., gene ID)
#'     \item \code{V2}: gene symbol
#'     \item \code{V4}: chromosome (e.g., "chr1")
#'     \item \code{V5}: start coordinate
#'     \item \code{V6}: end coordinate
#'     \item Optional: \code{sum} or \code{mean} column with summarized values
#'   }
#' @param feature_id A string used to label the metadata (e.g. "Expr" or "ATAC")
#'
#' @return A `GRanges` object with genomic coordinates and associated metadata
#'
#' @examples
#' \dontrun{
#' gr <- CreateGRangeObj(specific_feature = df, feature_id = "Expr")
#' }
#' @export
CreateGRangeObj <- function (specific_feature, feature_id) {
  feature_sub <- specific_feature[startsWith(specific_feature$V4, "chr"), ]

  gr <- GRanges(
    seqnames = gsub("chr", "", feature_sub$V4),
    ranges = IRanges(
      start = feature_sub$V5,
      end = feature_sub$V6
    )
  )

  mcols(gr)[[gsub(" ", "_", paste0(feature_id, "_Symbol"))]] <- feature_sub$V2
  mcols(gr)[[gsub(" ", "_", paste0(feature_id, "_id"))]] <- feature_sub$V1

  if (!is.null(feature_sub$sum)) {
    mcols(gr)[[gsub(" ", "_", paste0(feature_id, "_sum"))]] <- feature_sub$sum
  }
  if (!is.null(feature_sub$mean)) {
    mcols(gr)[[gsub(" ", "_", paste0(feature_id, "_mean"))]] <- feature_sub$mean
  }

  return(gr)
}



