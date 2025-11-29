#' Annotate each genomic region with its feature type
#'
#' Assigns the genomic feature of promoter, exon, intron, or intergenic to each interval
#'
#' @param gr GRanges object of genomic regions
#' @param txdb_info A TxDb object with gene coordinates on a certain reference genome
#' @param promoter_up Optionally, the number of bases upstream of a transcription start site to include in the promoter region (default 3000)
#' @param promoter_down Optionally, the number of bases downstream of a transcription start site to include in the promoter region (default 3000)
#' @return A GRanges object that includes a new column called 'feature_type' showing which type of genomic feature each region overlaps
#' @examples
#' library(GenomicRanges)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges(100000, 101000))
#' annotate_regions(gr, txdb)
#' @importFrom S4Vectors queryHits subjectHits mcols mcols<-
#' @importFrom GenomicFeatures promoters exons intronsByTranscript genes
#' @importFrom GenomicRanges findOverlaps distanceToNearest
#' @export
annotate_regions <- function(gr, txdb_info, promoter_up = 3000, promoter_down = 3000) {
  
  # Validate inputs
  if (!inherits(gr, "GRanges")) {
    stop("gr must be a GRanges object")
  }
  
  if (is.null(txdb_info)) {
    stop("txdb_info cannot be NULL")
  }
  
  if (promoter_up < 0 || promoter_down < 0) {
    stop("promoter_up and promoter_down must be non-negative")
  }
  
  promoters <- GenomicFeatures::promoters(txdb_info, upstream = promoter_up, downstream = promoter_down)
  exons <- GenomicFeatures::exons(txdb_info)
  introns_list <- GenomicFeatures::intronsByTranscript(txdb_info)
  introns <- unlist(introns_list, use.names = FALSE)
  
  feature_type <- rep("intergenic", length(gr))
  
  overlap_promoter <- GenomicRanges::findOverlaps(gr, promoters)
  if (length(overlap_promoter) > 0) {
    hits <- S4Vectors::queryHits(overlap_promoter)
    for (i in hits) {
      feature_type[i] <- "promoter"
    }
  }
  
  overlap_exon <- GenomicRanges::findOverlaps(gr, exons)
  if (length(overlap_exon) > 0) {
    hits <- S4Vectors::queryHits(overlap_exon)
    for (i in hits) {
      feature_type[i] <- "exon"
    }
  }
  
  overlap_intron <- GenomicRanges::findOverlaps(gr, introns)
  if (length(overlap_intron) > 0) {
    hits <- S4Vectors::queryHits(overlap_intron)
    for (i in hits) {
      feature_type[i] <- "intron"
    }
  }
  
  S4Vectors::mcols(gr)$feature_type <- feature_type
  return(gr)
}

#' Identify Nearest Genes and Calculate Distances
#'
#' For each genomic interval, identifies the closest gene (by distance from
#' interval to gene body) and calculates the genomic distance in base pairs.
#' A distance of 0 indicates the interval overlaps the gene body. Positive
#' distances indicate how far the interval is from the nearest gene. This
#' information is crucial for linking regulatory elements to their potential
#' target genes or for understanding the genomic context of your intervals.
#' Must be run after annotate_regions().
#' 
#' @param gr GRanges object of genomic regions
#' @param txdb_info A TxDb object with gene coordinates on a certain reference genome
#' @return A GRanges object with extra metadata columns 'nearest_gene' and 'distance_to_gene'
#' @examples
#' \dontrun{
#' # Requires TxDb.Hsapiens.UCSC.hg19.knownGene (not available on all systems)
#' library(GenomicRanges)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges(100000, 101000))
#' gr_annotated <- annotate_regions(gr, txdb)
#' nearest_gene(gr_annotated, txdb)
#' }
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges distanceToNearest
#' @importFrom S4Vectors subjectHits mcols mcols<-
#' @export
nearest_gene <- function(gr, txdb_info) {
  if (!("feature_type" %in% names(S4Vectors::mcols(gr)))) {
    stop("Input GRanges object must contain columns added by annotate_regions().")
  }
  
  if (is.null(txdb_info)) {
    stop("txdb_info cannot be NULL")
  }
  
  all_genes_list <- GenomicFeatures::genes(txdb_info, single.strand.genes.only = FALSE)
  all_genes <- unlist(all_genes_list)
  
  hits <- GenomicRanges::distanceToNearest(gr, all_genes)
  
  nearest_names <- names(all_genes)[S4Vectors::subjectHits(hits)]
  distances <- S4Vectors::mcols(hits)$distance
  
  S4Vectors::mcols(gr)$nearest_gene <- nearest_names
  S4Vectors::mcols(gr)$distance_to_gene <- distances
  return(gr)
}

#' Retrieve functional annotations for genes
#'
#' For each gene associated with genomic intervals, retrieve functional
#' information
#'
#' @param gr GRanges object with a metadata column 'nearest_gene' of gene IDs
#' @param orgdb An OrgDb object for mapping gene IDs to GO terms
#' @param go TRUE/FALSE, if TRUE, perform GO enrichment (default TRUE)
#' @param kegg TRUE/FALSE, if TRUE, perform KEGG pathway enrichment (default FALSE)
#' @return A list with elements 'go' and/or 'kegg', each containing an enrichment result data frame
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(100000, 101000))
#' GenomicRanges::mcols(gr)$nearest_gene <- "1"
#' functional_terms(gr, org.Hs.eg.db)
#' }
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom stats na.omit
#' @importFrom S4Vectors mcols
#' @export
functional_terms <- function(gr, orgdb, go = TRUE, kegg = FALSE) {
  if (is.null(orgdb)) {
    stop("orgdb cannot be NULL")
  }
  
  genes <- as.character(S4Vectors::mcols(gr)$nearest_gene)
  genes <- stats::na.omit(genes)
  
  if (length(genes) == 0) {
    message("No gene IDs found for enrichment analysis.")
    return(NULL)
  }
  
  results <- list()
  
  if (go) {
    go_res <- clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = orgdb,
      keyType = "ENTREZID",
      ont = "BP",
      readable = TRUE
    )
    results$go <- as.data.frame(go_res)
  }
  
  if (kegg) {
    kegg_res <- clusterProfiler::enrichKEGG(
      gene = genes,
      organism = "hsa"
    )
    results$kegg <- as.data.frame(kegg_res)
  }
  
  return(results)
}

#' Test for Statistical Enrichment of Genomic Features
#'
#' Performs Fisher's exact test to determine whether a specific genomic feature
#' type (promoter, exon, intron, or intergenic) is statistically over-represented
#' or under-represented in your dataset compared to a background set. Returns the
#' odds ratio (effect size), p-value (statistical significance), and contingency
#' table. An odds ratio > 1 indicates enrichment, while < 1 indicates depletion.
#' This is useful for determining if your genomic intervals are biased toward
#' certain functional elements.
#' 
#' @param gr GRanges object with a metadata column 'feature_type' of categories
#' @param category A character string specifying which category to test for enrichment. Valid options can be any value present in the metadata column. So, if using the 'feature_type' column, valid options are: "promoter", "exon", "intron", "intergenic".
#' @param background_gr A GRanges object representing the background set of regions.
#' @return A list with elements 'odds_ratio', 'p_value', and 'table' containing the contingency table, or NULL if the category is invalid
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges(c(100,200), c(150,250)))
#' mcols(gr)$feature_type <- c("promoter", "exon")
#' bg_gr <- GRanges(seqnames = "chr1", ranges = IRanges(c(50,150,250), c(120,220,300)))
#' mcols(bg_gr)$feature_type <- c("promoter", "exon", "intron")
#' fisher_enrichment(gr, "promoter", bg_gr)
#' @importFrom stats fisher.test
#' @importFrom S4Vectors mcols
#' @export
fisher_enrichment <- function(gr, category, background_gr) {
  
  # Validate inputs
  if (is.null(category) || !is.character(category)) {
    stop("category must be a character string")
  }
  
  if (length(background_gr) == 0) {
    stop("background_gr cannot be empty")
  }
  
  if (!("feature_type" %in% names(S4Vectors::mcols(gr)))) {
    stop("gr must have a 'feature_type' metadata column")
  }
  
  if (!("feature_type" %in% names(S4Vectors::mcols(background_gr)))) {
    stop("background_gr must have a 'feature_type' metadata column")
  }
  
  gr_feature <- as.character(S4Vectors::mcols(gr)$feature_type)
  bg_feature <- as.character(S4Vectors::mcols(background_gr)$feature_type)
  
  valid_categories <- unique(c(gr_feature, bg_feature))
  if (!(category %in% valid_categories)) {
    message("Invalid category: '", category, "'.")
    message("Valid categories are: ", paste(valid_categories, collapse = ", "))
    return(NULL)
  }
  
  # Calculate counts - background should NOT include the target regions
  in_category <- sum(gr_feature == category)
  not_in_category <- sum(gr_feature != category)
  
  # Background counts should be independent of target
  bg_in_category <- sum(bg_feature == category)
  bg_not_in_category <- sum(bg_feature != category)
  
  # Create contingency table
  contingency_table <- matrix(
    c(in_category, not_in_category, bg_in_category, bg_not_in_category),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Target", "Background"),
      c("InCategory", "NotInCategory")
    )
  )
  
  # Validate contingency table
  if (any(contingency_table < 0)) {
    message("Contingency table contains negative values")
    return(NULL)
  }
  
  test <- tryCatch({
    stats::fisher.test(contingency_table)
  }, error = function(e) {
    message("Fisher's exact test failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(test)) {
    return(NULL)
  }
  
  result <- list(
    odds_ratio = test$estimate,
    p_value = test$p.value,
    table = contingency_table
  )
  
  return(result)
}