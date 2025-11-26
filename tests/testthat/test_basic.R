# tests/testthat/test-basic.R
# Basic functionality tests for ReAnnotateR

# ============================================================================
# read_bed tests
# ============================================================================

test_that("read_bed reads valid BED files", {
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  expect_s4_class(gr, "GRanges")
  expect_true(length(gr) > 0)
  expect_true(all(!is.na(GenomicRanges::seqnames(gr))))
  expect_true(all(GenomicRanges::start(gr) > 0))
  expect_true(all(GenomicRanges::end(gr) >= GenomicRanges::start(gr)))
})

test_that("read_bed handles missing files", {
  suppressMessages(
    expect_null(read_bed("nonexistent_file.bed"))
  )
})

test_that("read_bed handles invalid input types", {
  expect_null(read_bed(NULL))
  expect_null(read_bed(123))
  expect_null(read_bed(NA))
})

test_that("read_bed parses BED columns correctly", {
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  
  # Check basic structure
  expect_true(all(c("seqnames", "ranges", "strand") %in% slotNames(gr)))
})

# ============================================================================
# annotate_regions tests
# ============================================================================

test_that("annotate_regions adds feature annotations", {
  txdb <- get_test_txdb()
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  
  expect_s4_class(gr_annotated, "GRanges")
  expect_true("feature_type" %in% colnames(S4Vectors::mcols(gr_annotated)))
  expect_true(all(S4Vectors::mcols(gr_annotated)$feature_type %in% 
                    c("promoter", "exon", "intron", "intergenic")))
  expect_equal(length(gr_annotated), length(gr))
})

test_that("annotate_regions handles empty GRanges", {
  txdb <- get_test_txdb()
  
  empty_gr <- GenomicRanges::GRanges()
  result <- suppressWarnings(annotate_regions(empty_gr, txdb))
  
  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0)
})

test_that("annotate_regions handles custom promoter regions", {
  txdb <- get_test_txdb()
  gr <- get_test_gr()
  
  result <- suppressWarnings(
    annotate_regions(gr, txdb, promoter_up = 5000, promoter_down = 1000)
  )
  
  expect_s4_class(result, "GRanges")
  expect_true("feature_type" %in% names(S4Vectors::mcols(result)))
})

# ============================================================================
# nearest_gene tests
# ============================================================================

test_that("nearest_gene identifies closest genes", {
  txdb <- get_test_txdb()
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  gr_nearest <- nearest_gene(gr_annotated, txdb_info = txdb)
  
  expect_s4_class(gr_nearest, "GRanges")
  expect_true("nearest_gene" %in% colnames(S4Vectors::mcols(gr_nearest)))
  expect_true("distance_to_gene" %in% colnames(S4Vectors::mcols(gr_nearest)))
  expect_true("feature_type" %in% colnames(S4Vectors::mcols(gr_nearest)))
  expect_true(is.numeric(S4Vectors::mcols(gr_nearest)$distance_to_gene))
  expect_true(all(S4Vectors::mcols(gr_nearest)$distance_to_gene >= 0))
  expect_equal(length(gr_nearest), length(gr_annotated))
})

test_that("nearest_gene handles empty input", {
  txdb <- get_test_txdb()
  
  empty_gr <- GenomicRanges::GRanges()
  empty_annotated <- suppressWarnings(annotate_regions(empty_gr, txdb))
  
  # Ensure feature_type column exists even for empty GRanges
  if (length(empty_annotated) == 0 && 
      !("feature_type" %in% names(S4Vectors::mcols(empty_annotated)))) {
    S4Vectors::mcols(empty_annotated)$feature_type <- character(0)
  }
  
  result <- nearest_gene(empty_annotated, txdb_info = txdb)
  
  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0)
})

test_that("nearest_gene requires annotated regions", {
  txdb <- get_test_txdb()
  gr <- get_test_gr()
  
  expect_error(
    nearest_gene(gr, txdb_info = txdb),
    regexp = "must contain columns added by annotate_regions"
  )
})

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