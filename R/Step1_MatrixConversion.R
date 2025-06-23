#' Convert a sparse Matrix Market file into a dense data.table
#'
#' This function performs Step 1 of the exam:
#' It loads a sparse expression matrix from a `.mtx.gz` file,
#' along with corresponding features and barcodes, and converts
#' everything into a dense data.table suitable for downstream analysis.
#'
#' @param mtx_path Character. Path to the Matrix Market file (.mtx.gz).
#' @param features_path Character. Path to the features file (e.g., features.tsv.gz).
#' @param barcodes_path Character. Path to the barcodes file (e.g., barcodes.tsv.gz).
#'
#' @return A `data.table` with rows as features (genes or peaks),
#' and columns as barcodes (cells), including a column `"row.names"` for feature IDs.
#'
#' @importFrom Matrix readMM
#' @importFrom data.table fread as.data.table
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' expr_dt <- ConvertMatrix(
#'   mtx_path = "data/expression/matrix.mtx.gz",
#'   features_path = "data/expression/features.tsv.gz",
#'   barcodes_path = "data/expression/barcodes.tsv.gz"
#' )
#' head(expr_dt)
#' }
#' @export
ConvertMatrix <- function(mtx_path, features_path, barcodes_path) {

  sparse_mat <- readMM(mtx_path)

  features <- fread(features_path, header = FALSE)
  barcodes <- fread(barcodes_path, header = FALSE)

  rownames(sparse_mat) <- features$V1
  colnames(sparse_mat) <- barcodes$V1

  full_mat <- as.matrix(sparse_mat)

  dt <- as.data.table(full_mat, keep.rownames = "row.names")

  return(dt)
}
