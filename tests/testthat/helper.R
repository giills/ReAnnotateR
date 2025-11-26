# Helper file for testthat
# This sets up test fixtures that are available to all test files

# Only load test data if the required packages are available
skip_if_no_txdb <- function() {
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    testthat::skip("TxDb.Hsapiens.UCSC.hg38.knownGene not available")
  }
}

skip_if_no_orgdb <- function() {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    testthat::skip("org.Hs.eg.db not available")
  }
}

# Create test fixture function that loads txdb only when needed
get_test_txdb <- function() {
  skip_if_no_txdb()
  TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
}

get_test_orgdb <- function() {
  skip_if_no_orgdb()
  org.Hs.eg.db::org.Hs.eg.db
}

# Create a simple test GRanges object
get_test_gr <- function() {
  GenomicRanges::GRanges(
    seqnames = rep("chr1", 5),
    ranges = IRanges::IRanges(
      start = c(10000, 50000, 100000, 150000, 200000),
      end = c(11000, 51000, 101000, 151000, 201000)
    )
  )
}