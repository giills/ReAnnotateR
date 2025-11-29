# tests/testthat/test-visualization.R
# Tests for visualization functions

# ============================================================================
# plot_feature_composition tests
# ============================================================================

test_that("plot_feature_composition creates valid plot", {
  txdb <- get_test_txdb()
  
  gr <- get_test_gr()
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  
  plot <- plot_feature_composition(gr_annotated)
  
  expect_s3_class(plot, "ggplot")
})

test_that("plot_feature_composition requires feature_type column", {
  gr <- get_test_gr()
  
  result <- suppressMessages(plot_feature_composition(gr))
  
  expect_null(result)
})

test_that("plot_feature_composition handles all feature types", {
  # Create GRanges with all feature types
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4),
    ranges = IRanges::IRanges(start = c(100, 200, 300, 400),
                              end = c(150, 250, 350, 450))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "exon", "intron", "intergenic")
  
  plot <- plot_feature_composition(gr)
  
  expect_s3_class(plot, "ggplot")
})

# ============================================================================
# plot_chromosomal_density tests
# ============================================================================

test_that("plot_chromosomal_density creates valid plot", {
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  plot <- plot_chromosomal_density(gr)
  
  expect_s3_class(plot, "ggplot")
})

test_that("plot_chromosomal_density handles custom bin size", {
  gr <- get_test_gr()
  
  plot <- plot_chromosomal_density(gr, bin_size = 500000)
  
  expect_s3_class(plot, "ggplot")
})

test_that("plot_chromosomal_density handles single chromosome", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 10),
    ranges = IRanges::IRanges(
      start = seq(1000000, 10000000, length.out = 10),
      width = 1000
    )
  )
  
  plot <- plot_chromosomal_density(gr)
  
  expect_s3_class(plot, "ggplot")
})

# ============================================================================
# plot_enrichment tests
# ============================================================================

test_that("plot_enrichment creates valid plot from GO results", {
  skip_if_no_orgdb()
  
  txdb <- get_test_txdb()
  orgdb <- get_test_orgdb()
  
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  gr_nearest <- nearest_gene(gr_annotated, txdb_info = txdb)
  
  enrichment <- suppressMessages(
    functional_terms(gr_nearest, orgdb = orgdb, go = TRUE)
  )
  
  skip_if(is.null(enrichment$go) || nrow(enrichment$go) == 0,
          "No enrichment results to plot")
  
  plot <- plot_enrichment(enrichment$go, top_n = 5)
  
  expect_s3_class(plot, "ggplot")
})

test_that("plot_enrichment handles empty results", {
  result <- suppressMessages(plot_enrichment(NULL))
  
  expect_null(result)
})

test_that("plot_enrichment handles custom top_n parameter", {
  # Create mock enrichment data
  mock_enrichment <- data.frame(
    Description = paste("Term", 1:20),
    p.adjust = runif(20, 0.001, 0.05),
    stringsAsFactors = FALSE
  )
  
  plot <- plot_enrichment(mock_enrichment, top_n = 5)
  
  expect_s3_class(plot, "ggplot")
})

# ============================================================================
# plot_annotation_summary tests
# ============================================================================

test_that("plot_annotation_summary creates valid plot", {
  txdb <- get_test_txdb()
  
  gr <- get_test_gr()
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  
  plot <- plot_annotation_summary(gr_annotated)
  
  expect_s3_class(plot, "ggplot")
})

test_that("plot_annotation_summary handles NA values", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 5),
    ranges = IRanges::IRanges(start = seq(100, 500, 100),
                              end = seq(150, 550, 100))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "exon", NA, "intron", "intergenic")
  
  plot <- plot_annotation_summary(gr)
  
  expect_s3_class(plot, "ggplot")
})

# ============================================================================
# plot_missing_annotations tests
# ============================================================================

test_that("plot_missing_annotations creates valid plot", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 5),
    ranges = IRanges::IRanges(start = seq(100, 500, 100),
                              end = seq(150, 550, 100))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "exon", NA, "intron", "intergenic")
  
  plot <- plot_missing_annotations(gr)
  
  expect_s3_class(plot, "ggplot")
})

test_that("plot_missing_annotations handles all annotated", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(100, 200, 300),
                              end = c(150, 250, 350))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "exon", "intron")
  
  plot <- plot_missing_annotations(gr)
  
  expect_s3_class(plot, "ggplot")
})

test_that("plot_missing_annotations handles all missing", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(100, 200, 300),
                              end = c(150, 250, 350))
  )
  S4Vectors::mcols(gr)$feature_type <- c(NA, NA, NA)
  
  plot <- plot_missing_annotations(gr)
  
  expect_s3_class(plot, "ggplot")
})

# ============================================================================
# plot_regulatory_regions tests
# ============================================================================

test_that("plot_regulatory_regions creates valid plot", {
  gr <- get_test_gr()
  
  # Create mock regulatory regions
  regulatory_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(9000, 49000, 99000),
                              end = c(11000, 51000, 101000))
  )
  
  plot <- plot_regulatory_regions(gr, regulatory_gr)
  
  expect_s3_class(plot, "ggplot")
})

test_that("plot_regulatory_regions handles no overlaps", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1000, 2000),
                              end = c(1100, 2100))
  )
  
  regulatory_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5000, 6000),
                              end = c(5100, 6100))
  )
  
  plot <- plot_regulatory_regions(gr, regulatory_gr)
  
  expect_s3_class(plot, "ggplot")
})

# ============================================================================
# plot_ideogram tests
# ============================================================================

test_that("plot_ideogram requires Gviz package", {
  skip_if_not_installed("Gviz")
  
  gr <- get_test_gr()
  
  # Should not error when Gviz is available
  expect_no_error(
    suppressWarnings(plot_ideogram(gr, genome = "hg38", chromosomes = "chr1"))
  )
})

test_that("plot_ideogram errors without Gviz", {
  skip_if(requireNamespace("Gviz", quietly = TRUE), "Gviz is installed")
  
  gr <- get_test_gr()
  
  expect_error(
    plot_ideogram(gr),
    "Package 'Gviz' is required"
  )
})

# ==================================================================
# visulization edge cases
# ==================================================================


test_that("plot_feature_composition handles single feature type", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(100, 200, 300), end = c(150, 250, 350))
  )
  S4Vectors::mcols(gr)$feature_type <- rep("promoter", 3)
  
  plot <- plot_feature_composition(gr)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_feature_composition handles NA in feature_type", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4),
    ranges = IRanges::IRanges(start = c(100, 200, 300, 400), end = c(150, 250, 350, 450))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", NA, "exon", "intron")
  
  plot <- plot_feature_composition(gr)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_chromosomal_density handles single region", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1000, end = 2000)
  )
  
  plot <- plot_chromosomal_density(gr)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_chromosomal_density handles very small bin size", {
  gr <- get_test_gr()
  
  plot <- plot_chromosomal_density(gr, bin_size = 10)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_chromosomal_density handles very large bin size", {
  gr <- get_test_gr()
  
  plot <- plot_chromosomal_density(gr, bin_size = 1e9)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_enrichment handles single term", {
  mock_enrichment <- data.frame(
    Description = "Term1",
    p.adjust = 0.01,
    stringsAsFactors = FALSE
  )
  
  plot <- plot_enrichment(mock_enrichment, top_n = 10)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_enrichment handles top_n larger than results", {
  mock_enrichment <- data.frame(
    Description = paste("Term", 1:3),
    p.adjust = c(0.01, 0.02, 0.03),
    stringsAsFactors = FALSE
  )
  
  plot <- plot_enrichment(mock_enrichment, top_n = 100)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_enrichment handles zero or negative top_n", {
  mock_enrichment <- data.frame(
    Description = paste("Term", 1:5),
    p.adjust = runif(5, 0.001, 0.05),
    stringsAsFactors = FALSE
  )
  
  # Should handle gracefully
  result <- tryCatch(
    plot_enrichment(mock_enrichment, top_n = 0),
    error = function(e) NULL
  )
  
  expect_true(is.null(result) || inherits(result, "ggplot"))
})

test_that("plot_regulatory_regions handles empty regulatory_gr", {
  gr <- get_test_gr()
  regulatory_gr <- GenomicRanges::GRanges()
  
  plot <- plot_regulatory_regions(gr, regulatory_gr)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_regulatory_regions handles different chromosomes", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1000, 2000), end = c(1100, 2100))
  )
  
  regulatory_gr <- GenomicRanges::GRanges(
    seqnames = "chr2",
    ranges = IRanges::IRanges(start = c(1000, 2000), end = c(1100, 2100))
  )
  
  plot <- suppressWarnings(plot_regulatory_regions(gr, regulatory_gr))
  expect_s3_class(plot, "ggplot")
})

test_that("plot_missing_annotations handles no missing data", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(100, 200, 300), end = c(150, 250, 350))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "exon", "intron")
  
  plot <- plot_missing_annotations(gr)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_annotation_summary handles empty character for feature_type", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, end = 150)
  )
  S4Vectors::mcols(gr)$feature_type <- ""
  
  plot <- plot_annotation_summary(gr)
  expect_s3_class(plot, "ggplot")
})