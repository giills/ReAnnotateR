#' Example BED file path
#'
#' Path to an example BED file included with the package for demonstration purposes.
#' This file contains genomic intervals that can be used to test annotation functions.
#'
#' @format A character string containing the file path
#' @examples
#' \dontrun{
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' }
"example_bed"

#' Example TxDb annotation object
#'
#' A TxDb object containing gene annotations for testing purposes.
#' This is a subset of gene annotations used in examples and tests.
#'
#' @format A TxDb object
#' @examples
#' \dontrun{
#' data(example_annotation)
#' }
"example_annotation"

#' Example enrichment database
#'
#' An OrgDb-like object containing gene-to-function mappings for testing purposes.
#' Used in functional enrichment analysis examples.
#'
#' @format An OrgDb-like object
#' @examples
#' \dontrun{
#' data(example_enrichment)
#' }
"example_enrichment"