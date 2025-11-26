# tests/testthat/test-io.R
# Tests for input/output functions

# ============================================================================
# export_results tests
# ============================================================================

test_that("export_results exports GRanges to TSV", {
  gr <- get_test_gr()
  
  temp_file <- tempfile(fileext = ".tsv")
  
  suppressMessages(
    export_results(gr, temp_file, format = "tsv")
  )
  
  expect_true(file.exists(temp_file))
  
  # Read back and verify
  df <- read.table(temp_file, header = TRUE, sep = "\t")
  expect_equal(nrow(df), length(gr))
  expect_true("chr" %in% colnames(df))
  expect_true("start" %in% colnames(df))
  expect_true("end" %in% colnames(df))
  
  unlink(temp_file)
})

test_that("export_results exports GRanges to CSV", {
  gr <- get_test_gr()
  
  temp_file <- tempfile(fileext = ".csv")
  
  suppressMessages(
    export_results(gr, temp_file, format = "csv")
  )
  
  expect_true(file.exists(temp_file))
  
  # Read back and verify
  df <- read.csv(temp_file)
  expect_equal(nrow(df), length(gr))
  
  unlink(temp_file)
})

test_that("export_results exports GRanges to BED", {
  gr <- get_test_gr()
  
  temp_file <- tempfile(fileext = ".bed")
  
  suppressMessages(
    export_results(gr, temp_file, format = "bed")
  )
  
  expect_true(file.exists(temp_file))
  
  unlink(temp_file)
})

test_that("export_results exports data frames", {
  df <- data.frame(
    chr = c("chr1", "chr2"),
    start = c(100, 200),
    end = c(150, 250),
    feature = c("promoter", "exon")
  )
  
  temp_file <- tempfile(fileext = ".tsv")
  
  suppressMessages(
    export_results(df, temp_file, format = "tsv")
  )
  
  expect_true(file.exists(temp_file))
  
  df_read <- read.table(temp_file, header = TRUE, sep = "\t")
  expect_equal(nrow(df_read), nrow(df))
  
  unlink(temp_file)
})

test_that("export_results exports enrichment results list", {
  # Create mock enrichment results
  enrichment <- list(
    go = data.frame(
      Description = c("Term1", "Term2"),
      pvalue = c(0.01, 0.05),
      p.adjust = c(0.02, 0.06)
    ),
    kegg = data.frame(
      Description = c("Pathway1", "Pathway2"),
      pvalue = c(0.03, 0.04)
    )
  )
  
  temp_base <- tempfile()
  
  suppressMessages(
    export_results(enrichment, paste0(temp_base, ".tsv"), format = "tsv")
  )
  
  # Check that both GO and KEGG files were created
  expect_true(file.exists(paste0(temp_base, "_go.tsv")))
  expect_true(file.exists(paste0(temp_base, "_kegg.tsv")))
  
  # Clean up
  unlink(paste0(temp_base, "_go.tsv"))
  unlink(paste0(temp_base, "_kegg.tsv"))
})

test_that("export_results handles invalid format", {
  gr <- get_test_gr()
  temp_file <- tempfile()
  
  result <- suppressMessages(
    export_results(gr, temp_file, format = "invalid")
  )
  
  expect_null(result)
})

test_that("export_results handles invalid input type", {
  temp_file <- tempfile()
  
  result <- suppressMessages(
    export_results("not_a_valid_object", temp_file, format = "tsv")
  )
  
  expect_null(result)
})

test_that("export_results preserves metadata columns", {
  txdb <- get_test_txdb()
  
  gr <- get_test_gr()
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  
  temp_file <- tempfile(fileext = ".tsv")
  
  suppressMessages(
    export_results(gr_annotated, temp_file, format = "tsv")
  )
  
  df <- read.table(temp_file, header = TRUE, sep = "\t")
  expect_true("feature_type" %in% colnames(df))
  
  unlink(temp_file)
})

test_that("export_results handles empty GRanges", {
  empty_gr <- GenomicRanges::GRanges()
  temp_file <- tempfile(fileext = ".tsv")
  
  suppressMessages(
    export_results(empty_gr, temp_file, format = "tsv")
  )
  
  expect_true(file.exists(temp_file))
  
  df <- read.table(temp_file, header = TRUE, sep = "\t")
  expect_equal(nrow(df), 0)
  
  unlink(temp_file)
})

test_that("export_results handles enrichment list with only GO", {
  enrichment <- list(
    go = data.frame(
      Description = "Term1",
      pvalue = 0.01
    )
  )
  
  temp_base <- tempfile()
  
  suppressMessages(
    export_results(enrichment, paste0(temp_base, ".tsv"), format = "tsv")
  )
  
  expect_true(file.exists(paste0(temp_base, "_go.tsv")))
  expect_false(file.exists(paste0(temp_base, "_kegg.tsv")))
  
  unlink(paste0(temp_base, "_go.tsv"))
})

test_that("export_results handles empty enrichment list", {
  enrichment <- list()
  
  temp_file <- tempfile(fileext = ".tsv")
  
  result <- suppressMessages(
    export_results(enrichment, temp_file, format = "tsv")
  )
  
  expect_null(result)
})

# ============================================================================
# Integration test: full workflow with export
# ============================================================================

test_that("full workflow with export works end-to-end", {
  skip_if_no_orgdb()
  
  txdb <- get_test_txdb()
  orgdb <- get_test_orgdb()
  
  bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
  skip_if(bed_file == "", "Example BED file not available")
  
  # Read and annotate
  gr <- read_bed(bed_file)
  gr_annotated <- suppressWarnings(annotate_regions(gr, txdb))
  gr_nearest <- nearest_gene(gr_annotated, txdb_info = txdb)
  
  # Export annotated regions
  temp_file1 <- tempfile(fileext = ".tsv")
  suppressMessages(export_results(gr_nearest, temp_file1, format = "tsv"))
  expect_true(file.exists(temp_file1))
  
  # Functional enrichment
  enrichment <- suppressMessages(
    functional_terms(gr_nearest, orgdb = orgdb, go = TRUE)
  )
  
  # Export enrichment (if results exist)
  if (!is.null(enrichment) && length(enrichment) > 0) {
    temp_file2 <- tempfile()
    suppressMessages(
      export_results(enrichment, paste0(temp_file2, ".tsv"), format = "tsv")
    )
    
    if (!is.null(enrichment$go) && nrow(enrichment$go) > 0) {
      expect_true(file.exists(paste0(temp_file2, "_go.tsv")))
      unlink(paste0(temp_file2, "_go.tsv"))
    }
  }
  
  unlink(temp_file1)
})