#' Annotate ATAC-seq Peaks with Overlapping Genes
#'
#' This function performs Step 5 of the exam:
#' reads a GTF annotation file, filters for protein-coding genes,
#' and then finds overlaps between ATAC-seq peaks and these genes. It returns a
#' `GRanges` object of ATAC peaks annotated with gene information.
#'
#' @param gtf_path A string. The path to the `Homo_sapiens.GRCh38.114.gtf.gz` annotation file.
#' @param atac_gr A `GRanges` object containing ATAC-seq peak ranges and metadata,
#'                typically generated from `PrepareGRanges()`.
#'
#' @return A `GRanges` object of overlapping ATAC peaks, with metadata columns:
#'   \itemize{
#'     \item \code{overlapping_gene_id}: Ensembl gene ID of the overlapping gene
#'     \item \code{overlapping_gene_name}: Gene name (symbol)
#'     \item \code{ATAC_sum}: Total ATAC signal per peak (inherited from input)
#'   }
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges findOverlaps
#'
#' @examples
#' \dontrun{
#' gtf_path <- "data/Homo_sapiens.GRCh38.114.gtf.gz"
#' annotated_atac <- AnnotateATACwithGenes(gtf_path, atac_gr)
#' }
#' @export
AnnotateATACwithGenes <- function(gtf_path, atac_gr) {

  gtf <- rtracklayer::import(gtf_path)

  gene_features <- gtf[gtf$type == "gene"]
  protein_coding_genes <- gene_features[gene_features$gene_biotype == "protein_coding"]

  gene_gr <- protein_coding_genes

  overlaps <- findOverlaps(atac_gr, gene_gr)

  atac_hits <- atac_gr[queryHits(overlaps)]

  mcols(atac_hits)$overlapping_gene_id <- gene_gr[subjectHits(overlaps)]$gene_id
  mcols(atac_hits)$overlapping_gene_name <- gene_gr[subjectHits(overlaps)]$gene_name

  mcols(atac_hits)$ATAC_sum <- mcols(atac_hits)$ATAC_sum

  return(atac_hits)
}

