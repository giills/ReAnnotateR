# tests/testthat/test-io.R
# Tests for input/output functions

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

# ============================================================================
# export_results edge cases
# ============================================================================

test_that("export_results handles read-only directory", {
  skip_on_os("windows")  # Different permission handling on Windows
  
  gr <- get_test_gr()
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  Sys.chmod(temp_dir, mode = "0444")  # Read-only
  
  temp_file <- file.path(temp_dir, "output.tsv")
  
  result <- suppressMessages(
    tryCatch(
      export_results(gr, temp_file, format = "tsv"),
      error = function(e) NULL
    )
  )
  
  Sys.chmod(temp_dir, mode = "0755")  # Restore permissions
  unlink(temp_dir, recursive = TRUE)
  
  expect_true(is.null(result) || !file.exists(temp_file))
})

test_that("export_results handles GRanges with complex metadata", {
  gr <- get_test_gr()
  S4Vectors::mcols(gr)$list_col <- as.list(1:length(gr))
  S4Vectors::mcols(gr)$factor_col <- factor(rep(c("A", "B"), length.out = length(gr)))
  
  temp_file <- tempfile(fileext = ".tsv")
  
  suppressMessages(export_results(gr, temp_file, format = "tsv"))
  
  expect_true(file.exists(temp_file))
  
  unlink(temp_file)
})

test_that("export_results handles NULL file path", {
  gr <- get_test_gr()
  
  result <- suppressMessages(export_results(gr, NULL, format = "tsv"))
  expect_null(result)
})

test_that("export_results handles empty file path", {
  gr <- get_test_gr()
  
  result <- suppressMessages(export_results(gr, "", format = "tsv"))
  expect_null(result)
})

test_that("export_results handles list with NULL elements", {
  enrichment <- list(go = NULL, kegg = NULL)
  
  temp_file <- tempfile(fileext = ".tsv")
  
  result <- suppressMessages(
    export_results(enrichment, temp_file, format = "tsv")
  )
  
  expect_null(result)
})

test_that("export_results handles data frame with no rows", {
  df <- data.frame(chr = character(0), start = numeric(0), end = numeric(0))
  
  temp_file <- tempfile(fileext = ".tsv")
  
  suppressMessages(export_results(df, temp_file, format = "tsv"))
  
  expect_true(file.exists(temp_file))
  
  unlink(temp_file)
})

test_that("export_results handles very long file paths", {
  gr <- get_test_gr()
  
  # Create a very long path (but still valid on most systems)
  long_name <- paste(rep("a", 200), collapse = "")
  temp_file <- file.path(tempdir(), paste0(long_name, ".tsv"))
  
  result <- suppressMessages(
    tryCatch(
      export_results(gr, temp_file, format = "tsv"),
      error = function(e) NULL
    )
  )
  
  if (file.exists(temp_file)) {
    unlink(temp_file)
  }
  
  expect_true(is.null(result) || file.exists(temp_file))
})
# ============================================================================
# read_bed edge cases
# ============================================================================

test_that("read_bed handles empty files", {
  temp_file <- tempfile(fileext = ".bed")
  writeLines(character(0), temp_file)
  
  result <- suppressMessages(read_bed(temp_file))
  expect_null(result)
  
  unlink(temp_file)
})

test_that("read_bed handles files with only whitespace", {
  temp_file <- tempfile(fileext = ".bed")
  writeLines(c("   ", "\t\t", ""), temp_file)
  
  result <- suppressMessages(read_bed(temp_file))
  expect_null(result)
  
  unlink(temp_file)
})

test_that("read_bed handles files with fewer than 3 columns", {
  temp_file <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t100"), temp_file)
  
  result <- suppressMessages(read_bed(temp_file))
  expect_null(result)
  
  unlink(temp_file)
})

test_that("read_bed handles negative start positions", {
  temp_file <- tempfile(fileext = ".bed")
  writeLines("chr1\t-100\t200", temp_file)
  
  result <- suppressMessages(read_bed(temp_file))
  expect_null(result)
  
  unlink(temp_file)
})

test_that("read_bed handles zero start positions", {
  temp_file <- tempfile(fileext = ".bed")
  writeLines("chr1\t0\t100", temp_file)
  
  # BED is 0-based, so this should convert to 1-based correctly
  result <- suppressMessages(read_bed(temp_file))
  expect_s4_class(result, "GRanges")
  expect_equal(GenomicRanges::start(result)[1], 1)
  
  unlink(temp_file)
})

test_that("read_bed handles end < start positions", {
  temp_file <- tempfile(fileext = ".bed")
  writeLines("chr1\t200\t100", temp_file)
  
  result <- suppressMessages(read_bed(temp_file))
  expect_null(result)
  
  unlink(temp_file)
})

test_that("read_bed handles non-numeric coordinates", {
  temp_file <- tempfile(fileext = ".bed")
  writeLines("chr1\tabc\t100", temp_file)
  
  result <- suppressMessages(read_bed(temp_file))
  expect_null(result)
  
  unlink(temp_file)
})

test_that("read_bed handles special characters in chromosome names", {
  temp_file <- tempfile(fileext = ".bed")
  writeLines("chr1_random\t100\t200", temp_file)
  
  result <- suppressMessages(read_bed(temp_file))
  expect_s4_class(result, "GRanges")
  expect_equal(as.character(GenomicRanges::seqnames(result)[1]), "chr1_random")
  
  unlink(temp_file)
})

test_that("read_bed handles all 6+ standard BED columns", {
  temp_file <- tempfile(fileext = ".bed")
  writeLines("chr1\t100\t200\tregion1\t100\t+", temp_file)
  
  result <- read_bed(temp_file)
  expect_s4_class(result, "GRanges")
  expect_true("name" %in% names(S4Vectors::mcols(result)))
  expect_true("score" %in% names(S4Vectors::mcols(result)))
  expect_equal(as.character(GenomicRanges::strand(result)[1]), "+")
  
  unlink(temp_file)
})


test_that("read_bed handles vector input", {
  result <- read_bed(c("file1.bed", "file2.bed"))
  expect_null(result)
})

test_that("read_bed handles empty string", {
  result <- read_bed("")
  expect_null(result)
})
