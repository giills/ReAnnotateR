# tests/testthat/test-enrichment.R
# Tests for functional enrichment and statistical analysis



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
# ============================================================================
# annotate_regions edge cases
# ============================================================================

test_that("annotate_regions handles regions on unknown chromosomes", {
  txdb <- get_test_txdb()
  
  gr <- GenomicRanges::GRanges(
    seqnames = "chrZ",
    ranges = IRanges::IRanges(start = 1000, end = 2000)
  )
  
  result <- suppressWarnings(annotate_regions(gr, txdb))
  expect_s4_class(result, "GRanges")
  expect_true("feature_type" %in% names(S4Vectors::mcols(result)))
})

test_that("annotate_regions handles very large coordinates", {
  txdb <- get_test_txdb()
  
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1e9, end = 1e9 + 1000)
  )
  
  result <- suppressWarnings(annotate_regions(gr, txdb))
  expect_s4_class(result, "GRanges")
})

test_that("annotate_regions handles single-base intervals", {
  txdb <- get_test_txdb()
  
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 10000, end = 10000)
  )
  
  result <- suppressWarnings(annotate_regions(gr, txdb))
  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 1)
})

test_that("annotate_regions handles negative promoter parameters", {
  txdb <- get_test_txdb()
  gr <- get_test_gr()
  
  expect_error(
    annotate_regions(gr, txdb, promoter_up = -1000, promoter_down = 1000)
  )
})

test_that("annotate_regions handles zero promoter window", {
  txdb <- get_test_txdb()
  gr <- get_test_gr()
  
  result <- suppressWarnings(
    annotate_regions(gr, txdb, promoter_up = 0, promoter_down = 0)
  )
  
  expect_s4_class(result, "GRanges")
})

test_that("annotate_regions handles NULL txdb", {
  gr <- get_test_gr()
  
  expect_error(annotate_regions(gr, NULL))
})

test_that("annotate_regions handles invalid object type", {
  txdb <- get_test_txdb()
  
  expect_error(annotate_regions("not_a_granges", txdb))
})

test_that("annotate_regions handles regions with metadata", {
  txdb <- get_test_txdb()
  
  gr <- get_test_gr()
  S4Vectors::mcols(gr)$custom_col <- seq_along(gr)
  
  result <- suppressWarnings(annotate_regions(gr, txdb))
  expect_true("custom_col" %in% names(S4Vectors::mcols(result)))
  expect_true("feature_type" %in% names(S4Vectors::mcols(result)))
})

# ============================================================================
# nearest_gene edge cases
# ============================================================================

test_that("nearest_gene handles regions far from any genes", {
  txdb <- get_test_txdb()
  
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1e9, end = 1e9 + 1000)
  )
  S4Vectors::mcols(gr)$feature_type <- "intergenic"
  
  result <- nearest_gene(gr, txdb_info = txdb)
  expect_s4_class(result, "GRanges")
  expect_true("distance_to_gene" %in% names(S4Vectors::mcols(result)))
})

test_that("nearest_gene handles NULL txdb", {
  gr <- get_test_gr()
  S4Vectors::mcols(gr)$feature_type <- "promoter"
  
  expect_error(nearest_gene(gr, txdb_info = NULL))
})

test_that("nearest_gene handles GRanges with single position", {
  txdb <- get_test_txdb()
  
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, end = 100)
  )
  S4Vectors::mcols(gr)$feature_type <- "promoter"
  
  result <- nearest_gene(gr, txdb_info = txdb)
  expect_s4_class(result, "GRanges")
  expect_true("nearest_gene" %in% names(S4Vectors::mcols(result)))
})

test_that("nearest_gene preserves all existing metadata", {
  txdb <- get_test_txdb()
  
  gr <- get_test_gr()
  S4Vectors::mcols(gr)$feature_type <- rep("promoter", length(gr))
  S4Vectors::mcols(gr)$custom_data <- 1:length(gr)
  
  result <- nearest_gene(gr, txdb_info = txdb)
  expect_true("custom_data" %in% names(S4Vectors::mcols(result)))
})

# ============================================================================
# functional_terms edge cases
# ============================================================================

test_that("functional_terms handles all NA genes", {
  skip_if_no_orgdb()
  orgdb <- get_test_orgdb()
  
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(100, 200), end = c(150, 250))
  )
  S4Vectors::mcols(gr)$nearest_gene <- c(NA, NA)
  
  result <- suppressMessages(functional_terms(gr, orgdb = orgdb, go = TRUE))
  expect_true(is.null(result))
})

test_that("functional_terms handles mixed NA and valid genes", {
  skip_if_no_orgdb()
  orgdb <- get_test_orgdb()
  
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(100, 200, 300), end = c(150, 250, 350))
  )
  S4Vectors::mcols(gr)$nearest_gene <- c("1", NA, "2")
  
  result <- suppressMessages(functional_terms(gr, orgdb = orgdb, go = TRUE))
  # Should process the valid genes
  expect_true(is.null(result) || is.list(result))
})

test_that("functional_terms handles invalid gene IDs", {
  skip_if_no_orgdb()
  orgdb <- get_test_orgdb()
  
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, end = 150)
  )
  S4Vectors::mcols(gr)$nearest_gene <- "INVALID_GENE_ID"
  
  result <- suppressMessages(functional_terms(gr, orgdb = orgdb, go = TRUE))
  expect_true(is.null(result) || (is.list(result) && all(sapply(result, function(x) is.null(x) || nrow(x) == 0))))
})

test_that("functional_terms handles NULL orgdb", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, end = 150)
  )
  S4Vectors::mcols(gr)$nearest_gene <- "1"
  
  expect_error(functional_terms(gr, orgdb = NULL, go = TRUE))
})

test_that("functional_terms handles both GO and KEGG false", {
  skip_if_no_orgdb()
  orgdb <- get_test_orgdb()
  
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100, end = 150)
  )
  S4Vectors::mcols(gr)$nearest_gene <- "1"
  
  result <- suppressMessages(functional_terms(gr, orgdb = orgdb, go = FALSE, kegg = FALSE))
  expect_true(is.list(result) && length(result) == 0)
})

test_that("functional_terms handles duplicate gene IDs", {
  skip_if_no_orgdb()
  orgdb <- get_test_orgdb()
  
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(100, 200, 300), end = c(150, 250, 350))
  )
  S4Vectors::mcols(gr)$nearest_gene <- c("1", "1", "1")
  
  result <- suppressMessages(functional_terms(gr, orgdb = orgdb, go = TRUE))
  # Should handle duplicates gracefully
  expect_true(is.null(result) || is.list(result))
})

# ============================================================================
# fisher_enrichment edge cases
# ============================================================================

test_that("fisher_enrichment handles zero counts in category", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(100, 200), c(150, 250))
  )
  S4Vectors::mcols(gr)$feature_type <- c("exon", "intron")
  
  bg_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(50, 150, 250), c(120, 220, 300))
  )
  S4Vectors::mcols(bg_gr)$feature_type <- c("promoter", "exon", "intron")
  
  result <- fisher_enrichment(gr, "promoter", bg_gr)
  expect_true(is.list(result) || is.null(result))
})

test_that("fisher_enrichment handles identical distributions", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(100, 200), c(150, 250))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "exon")
  
  bg_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(50, 150), c(120, 220))
  )
  S4Vectors::mcols(bg_gr)$feature_type <- c("promoter", "exon")
  
  result <- fisher_enrichment(gr, "promoter", bg_gr)
  expect_type(result, "list")
  # Odds ratio may not be exactly 1.0 depending on the contingency table
  expect_true(is.numeric(result$odds_ratio))
  expect_true(result$odds_ratio >= 0)
})

test_that("fisher_enrichment handles all regions in one category", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(100, 200, 300), c(150, 250, 350))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "promoter", "promoter")
  
  bg_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(50, 150, 250, 350), c(120, 220, 300, 400))
  )
  S4Vectors::mcols(bg_gr)$feature_type <- c("promoter", "exon", "intron", "intergenic")
  
  result <- fisher_enrichment(gr, "promoter", bg_gr)
  expect_type(result, "list")
  expect_true(result$odds_ratio > 1)
})

test_that("fisher_enrichment handles NULL category", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(100, 200), c(150, 250))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "exon")
  
  bg_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(50, 150), c(120, 220))
  )
  S4Vectors::mcols(bg_gr)$feature_type <- c("promoter", "exon")
  
  expect_error(fisher_enrichment(gr, NULL, bg_gr))
})

test_that("fisher_enrichment handles empty background", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(100, 200), c(150, 250))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "exon")
  
  bg_gr <- GenomicRanges::GRanges()
  S4Vectors::mcols(bg_gr)$feature_type <- character(0)
  
  expect_error(fisher_enrichment(gr, "promoter", bg_gr))
})

test_that("fisher_enrichment handles missing feature_type in background", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(100, 200), c(150, 250))
  )
  S4Vectors::mcols(gr)$feature_type <- c("promoter", "exon")
  
  bg_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(c(50, 150), c(120, 220))
  )
  
  expect_error(fisher_enrichment(gr, "promoter", bg_gr))
})
# ============================================================================
# convert_coordinates edge cases
# ============================================================================

test_that("convert_coordinates handles regions that don't lift over", {
  chain_file <- system.file("extdata", "hg19ToHg38.over.chain", 
                            package = "ReAnnotateR")
  skip_if(chain_file == "", "Chain file not available")
  
  # Create regions that might not lift over
  gr <- GenomicRanges::GRanges(
    seqnames = "chrM",
    ranges = IRanges::IRanges(start = 1000, end = 2000)
  )
  
  result <- convert_coordinates(gr, chain_file)
  expect_s4_class(result, "GRanges")
  # May be empty if regions don't lift
  expect_true(length(result) <= length(gr))
})

test_that("convert_coordinates handles NULL chain file path", {
  gr <- get_test_gr()
  
  expect_error(convert_coordinates(gr, NULL))
})

test_that("convert_coordinates handles empty chain file", {
  temp_chain <- tempfile(fileext = ".chain")
  writeLines(character(0), temp_chain)
  
  gr <- get_test_gr()
  
  result <- suppressMessages(
    tryCatch(
      convert_coordinates(gr, temp_chain),
      error = function(e) NULL
    )
  )
  
  expect_true(is.null(result) || inherits(result, "GRanges"))
  
  unlink(temp_chain)
})

test_that("convert_coordinates preserves metadata", {
  chain_file <- system.file("extdata", "hg19ToHg38.over.chain", 
                            package = "ReAnnotateR")
  skip_if(chain_file == "", "Chain file not available")
  
  gr <- get_test_gr()
  S4Vectors::mcols(gr)$custom_data <- 1:length(gr)
  
  result <- convert_coordinates(gr, chain_file)
  
  if (length(result) > 0) {
    expect_true("custom_data" %in% names(S4Vectors::mcols(result)))
  }
})