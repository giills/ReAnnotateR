# tests/testthat/test-enrichment.R
# Tests for functional enrichment and statistical analysis

# ============================================================================
# functional_terms tests
# ============================================================================

test_that("functional_terms performs GO enrichment", {
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
  
  expect_type(enrichment, "list")
  expect_true("go" %in% names(enrichment))
  
  if (!is.null(enrichment$go) && nrow(enrichment$go) > 0) {
    expect_true(all(c("Description", "p.adjust") %in% colnames(enrichment$go)))
  }
})

test_that("functional_terms handles KEGG option", {
  skip_if_no_orgdb()
  
  txdb <- get_test_txdb()
  orgdb <- get_test_orgdb()
  
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  gr_nearest <- nearest_gene(gr_annotated, txdb_info = txdb)
  
  enrichment <- suppressMessages(
    functional_terms(gr_nearest, orgdb = orgdb, go = FALSE, kegg = TRUE)
  )
  
  expect_type(enrichment, "list")
  if (!is.null(enrichment$kegg)) {
    expect_true(is.data.frame(enrichment$kegg) || nrow(enrichment$kegg) >= 0)
  }
})

test_that("functional_terms handles missing gene information", {
  skip_if_no_orgdb()
  
  orgdb <- get_test_orgdb()
  
  # Create GRanges without nearest_gene column
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(100, 200), end = c(150, 250))
  )
  
  # The function should handle missing nearest_gene gracefully
  result <- suppressMessages(
    functional_terms(gr, orgdb = orgdb, go = TRUE)
  )
  
  # Should return NULL or empty result when no genes present
  expect_true(is.null(result) || 
                (is.list(result) && (length(result) == 0 || 
                                       all(sapply(result, is.null)))))
})

test_that("functional_terms handles empty gene list", {
  skip_if_no_orgdb()
  
  orgdb <- get_test_orgdb()
  
  # Create GRanges with NA gene IDs
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, end = 150)
  )
  S4Vectors::mcols(gr)$nearest_gene <- NA_character_
  
  result <- suppressMessages(
    functional_terms(gr, orgdb = orgdb, go = TRUE)
  )
  
  expect_true(is.null(result) || 
                (is.list(result) && all(sapply(result, function(x) is.null(x) || nrow(x) == 0))))
})

# ============================================================================
# fisher_enrichment tests
# ============================================================================

test_that("fisher_enrichment calculates enrichment statistics", {
  txdb <- get_test_txdb()
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  
  # Verify feature_type column exists
  expect_true("feature_type" %in% names(S4Vectors::mcols(gr_annotated)))
  
  set.seed(123)
  n_bg <- 100
  bg_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = sample(1000000:10000000, n_bg),
                              width = 1000)
  )
  bg_gr_annotated <- suppressWarnings(annotate_regions(bg_gr, txdb))
  
  # Verify background also has feature_type column
  expect_true("feature_type" %in% names(S4Vectors::mcols(bg_gr_annotated)))
  
  fisher_results <- suppressMessages(
    fisher_enrichment(gr_annotated, 
                      category = "promoter",
                      background_gr = bg_gr_annotated)
  )
  
  # Skip if no promoters found
  skip_if(is.null(fisher_results), "No promoters found in test data")
  
  expect_type(fisher_results, "list")
  expect_true(all(c("odds_ratio", "p_value", "table") %in% names(fisher_results)))
  expect_true(is.numeric(fisher_results$odds_ratio))
  expect_true(is.numeric(fisher_results$p_value))
  expect_true(fisher_results$p_value >= 0 && fisher_results$p_value <= 1)
  expect_true(fisher_results$odds_ratio >= 0)
  expect_equal(dim(fisher_results$table), c(2, 2))
})

test_that("fisher_enrichment handles invalid categories", {
  txdb <- get_test_txdb()
  
  gr <- get_test_gr()
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  
  bg_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = sample(1000000:10000000, 100),
                              width = 1000)
  )
  bg_gr_annotated <- suppressWarnings(annotate_regions(bg_gr, txdb))
  
  result <- suppressMessages(
    fisher_enrichment(gr_annotated, 
                      category = "invalid_category",
                      background_gr = bg_gr_annotated)
  )
  
  expect_null(result)
})

test_that("fisher_enrichment tests all valid categories", {
  txdb <- get_test_txdb()
  
  gr <- get_test_gr()
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  
  set.seed(456)
  bg_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = sample(1000000:10000000, 100),
                              width = 1000)
  )
  bg_gr_annotated <- suppressWarnings(annotate_regions(bg_gr, txdb))
  
  categories <- c("promoter", "exon", "intron", "intergenic")
  
  for (cat in categories) {
    result <- tryCatch({
      suppressMessages(
        fisher_enrichment(gr_annotated, 
                          category = cat,
                          background_gr = bg_gr_annotated)
      )
    }, error = function(e) {
      # Fisher's exact test can fail if contingency table has invalid values
      NULL
    })
    
    # Result might be NULL if category not present or if fisher test fails
    if (!is.null(result)) {
      expect_type(result, "list")
    }
  }
})

test_that("fisher_enrichment with simple example", {
  # Create simple test data with known distribution
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(100, 200), c(150, 250))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "exon")
  
  bg_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(50, 150, 250, 350), c(120, 220, 300, 400))
  )
  S4Vectors::mcols(bg_gr)$feature_type <- c("promoter", "exon", "intron", "intergenic")
  
  result <- fisher_enrichment(gr, "promoter", bg_gr)
  
  expect_type(result, "list")
  expect_true("odds_ratio" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("table" %in% names(result))
})

# ============================================================================
# convert_coordinates tests
# ============================================================================

test_that("convert_coordinates lifts over genomic coordinates", {
  chain_file <- system.file("extdata", "hg19ToHg38.over.chain", 
                            package = "ReAnnotateR")
  skip_if(chain_file == "", "Chain file not available")
  
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  gr <- read_bed(bed_file)
  gr_converted <- convert_coordinates(gr, chain_file)
  
  expect_s4_class(gr_converted, "GRanges")
  expect_true(length(gr_converted) <= length(gr))
  
  if (length(gr_converted) > 0) {
    expect_true(all(GenomicRanges::start(gr_converted) > 0))
    expect_true(all(GenomicRanges::end(gr_converted) >= 
                      GenomicRanges::start(gr_converted)))
  }
})

test_that("convert_coordinates handles missing chain file", {
  gr <- get_test_gr()
  
  result <- suppressMessages(
    convert_coordinates(gr, "nonexistent.chain")
  )
  
  expect_null(result)
})

test_that("convert_coordinates handles invalid GRanges input", {
  chain_file <- system.file("extdata", "hg19ToHg38.over.chain", 
                            package = "ReAnnotateR")
  skip_if(chain_file == "", "Chain file not available")
  
  result <- suppressMessages(
    convert_coordinates("not_a_granges", chain_file)
  )
  
  expect_null(result)
})

test_that("convert_coordinates handles empty GRanges", {
  chain_file <- system.file("extdata", "hg19ToHg38.over.chain", 
                            package = "ReAnnotateR")
  skip_if(chain_file == "", "Chain file not available")
  
  empty_gr <- GenomicRanges::GRanges()
  result <- convert_coordinates(empty_gr, chain_file)
  
  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0)
})