#' Get a GRanges object from a BED file
#'
#' Checks for invalid chromosome start and end positions
#' @param file path to a BED file
#' @return `GRanges` object containing the intervals from the BED file
#' @examples
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' @import rtracklayer
#' @import GenomicRanges
#' @export
read_bed <- function(file) {
  if (!file.exists(file)) {
    print("This file does not exist")
    return(NULL)
  }
  gr <- rtracklayer::import(file, format = "BED")
  
  start_positions <- start(gr)
  for (i in seq_along(start(gr))) {
    if (start(gr)[i] < 1) {
      print("All start positions must be >= 1")
      return(NULL)
    }
  }
  
  end_positions <- end(gr)
  for (i in seq_along(start(gr))) {
    if (end(gr)[i] < start(gr)[i]) {
      print("All end positions must be >= start positions")
      return(NULL)
    }
  }
  
  return(gr)
}

#' Convert genomic coordinates between assemblies
#'
#' @param gr a `GRanges` object
#' @param chain_file A path to the chain file
#' @return A `GRanges` object with coordinates lifted over to the new assembly
#' @examples
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' chain_file <- system.file("extdata", "hg19ToHg38.over.chain", package = "ReAnnotateR")
#' gr_new <- convert_coordinates(gr, chain_file)
#' @import rtracklayer
#' @import GenomicRanges
#' @export
convert_coordinates <- function(gr, chain_file) {
  if (!inherits(gr, "GRanges")) {
    print("gr must be a GRanges object")
    return(NULL)
  }
  
  if (!file.exists(chain_file)) {
    print("Chain file does not exist.")
    return(NULL)
  }
  
  chain <- rtracklayer::import.chain(chain_file)
  gr_new_list <- rtracklayer::liftOver(gr, chain)
  gr_new <- unlist(gr_new_list)
  
  return(gr_new)
}
