% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Step6a_CreateProteinCodingGR.R
\name{CreateProteinCodingGR}
\alias{CreateProteinCodingGR}
\title{Create GRanges Object for Protein-Coding Genes}
\usage{
CreateProteinCodingGR(gtf_file)
}
\arguments{
\item{gtf_file}{Character string. Path to the Homo_sapiens.GRCh38.114.gtf.gz file.}
}
\value{
A \code{GRanges} object containing only protein-coding genes,
with two metadata columns: \code{gene_id} and \code{gene_name}.
}
\description{
This helper function is part of Step 6 of the analysis.
It imports a GTF annotation file and filters it to keep only
features of type "gene" that are protein-coding. The output is
a GRanges object ready for downstream overlap analysis.
}
\examples{
\dontrun{
gtf_path <- "data/Homo_sapiens.GRCh38.114.gtf.gz"
gene_gr <- CreateProteinCodingGR(gtf_path)
}

}
