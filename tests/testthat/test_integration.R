# tests/testthat/test-basic.R
# Integration and package-level tests for ReAnnotateR

# ============================================================================
# Integration tests
# ============================================================================

test_that("full annotation pipeline works end-to-end", {
  skip_if_no_orgdb()
  
  txdb <- get_test_txdb()
  orgdb <- get_test_orgdb()
  
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  expect_s4_class(gr, "GRanges")
  
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  expect_true("feature_type" %in% colnames(S4Vectors::mcols(gr_annotated)))
  
  gr_nearest <- nearest_gene(gr_annotated, txdb_info = txdb)
  expect_true("nearest_gene" %in% colnames(S4Vectors::mcols(gr_nearest)))
  
  enrichment <- suppressMessages(
    functional_terms(gr_nearest, orgdb = orgdb, go = TRUE)
  )
  expect_type(enrichment, "list")
})

test_that("pipeline preserves region count through annotation steps", {
  txdb <- get_test_txdb()
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  n_regions <- length(gr)
  
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  expect_equal(length(gr_annotated), n_regions)
  
  gr_nearest <- nearest_gene(gr_annotated, txdb_info = txdb)
  expect_equal(length(gr_nearest), n_regions)
})

#=======================================================
# Missing shiny package test
#=======================================================

test_that("runReAnnotateR handles missing shiny package", {
  skip_if(requireNamespace("shiny", quietly = TRUE), "shiny is installed")
  
  expect_error(runReAnnotateR())
})