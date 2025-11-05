library(ReAnnotateR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load test data
bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read and annotate
gr <- ReAnnotateR::read_bed(bed_file)
gr_annotated <- ReAnnotateR::annotate_regions(gr, txdb_info = txdb)

cat("=== Checking annotated regions ===\n")
cat("Length of gr_annotated:", length(gr_annotated), "\n")
cat("Columns:", paste(colnames(GenomicRanges::mcols(gr_annotated)), collapse = ", "), "\n")
cat("Has feature_type?", "feature_type" %in% names(GenomicRanges::mcols(gr_annotated)), "\n")

if ("feature_type" %in% names(GenomicRanges::mcols(gr_annotated))) {
  cat("\nFeature type counts:\n")
  print(table(GenomicRanges::mcols(gr_annotated)$feature_type))
}

# Create background
cat("\n=== Creating background ===\n")
set.seed(123)
n_bg <- 100
bg_gr <- GenomicRanges::GRanges(
  seqnames = "chr1",
  ranges = IRanges::IRanges(start = sample(1000000:10000000, n_bg),
                            width = 1000)
)
cat("Background regions created:", length(bg_gr), "\n")

# Annotate background
cat("\n=== Annotating background ===\n")
bg_gr_annotated <- ReAnnotateR::annotate_regions(bg_gr, txdb_info = txdb)
cat("Background annotated length:", length(bg_gr_annotated), "\n")
cat("Background columns:", paste(colnames(GenomicRanges::mcols(bg_gr_annotated)), collapse = ", "), "\n")

if ("feature_type" %in% names(GenomicRanges::mcols(bg_gr_annotated))) {
  cat("\nBackground feature type counts:\n")
  print(table(GenomicRanges::mcols(bg_gr_annotated)$feature_type))
}

# Try to see the fisher_enrichment function
cat("\n=== Inspecting fisher_enrichment function ===\n")
tryCatch({
  print(ReAnnotateR::fisher_enrichment)
}, error = function(e) {
  cat("Cannot print function\n")
})

# Try calling fisher_enrichment
cat("\n=== Calling fisher_enrichment ===\n")
tryCatch({
  fisher_results <- ReAnnotateR::fisher_enrichment(
    gr_annotated, 
    category = "promoter",
    background_gr = bg_gr_annotated
  )
  
  cat("Result type:", typeof(fisher_results), "\n")
  cat("Result class:", class(fisher_results), "\n")
  cat("Is NULL?", is.null(fisher_results), "\n")
  
  if (!is.null(fisher_results)) {
    cat("Names:", paste(names(fisher_results), collapse = ", "), "\n")
    print(fisher_results)
  }
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  cat("Traceback:\n")
  print(traceback())
})

# Test different categories
cat("\n=== Testing all categories ===\n")
categories <- c("promoter", "exon", "intron", "intergenic")
for (cat_name in categories) {
  result <- ReAnnotateR::fisher_enrichment(
    gr_annotated, 
    category = cat_name,
    background_gr = bg_gr_annotated
  )
  cat(sprintf("Category '%s': %s\n", cat_name, 
              if(is.null(result)) "NULL" else "SUCCESS"))
}