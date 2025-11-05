#' Annotate each genomic region with it's feature type
#'
#' Assigns the genomic feature of promoter, exon, intron, or intergenic to each interval
#'
#' @param gr GRanges object of genomic regions
#' @param txdb_info A TxDb object with gene coordinates on a certain reference genome
#' @param promoter_up Optionally, the number of bases upstream of a transcription start site to include in the promoter region (default 3000)
#' @param promoter_down Optionally, the number of bases downstream of a transcription start site to include in the promoter region (default 3000)
#' @return A GRanges object that includes a new column called 'feature' showing which type of genomic feature each region overlaps
#' @examples
#' library(GenomicRanges)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges(100000, 101000))
#' annotate_regions(gr, txdb)
#' @importFrom GenomicFeatures promoters exons intronsByTranscript genes
#' @export
annotate_regions <- function(gr, txdb_info, promoter_up = 3000, promoter_down = 3000) {
  
  promoters <- GenomicFeatures::promoters(txdb_info, upstream = promoter_up, downstream = promoter_down)
  exons <- GenomicFeatures::exons(txdb_info)
  introns_list <- GenomicFeatures::intronsByTranscript(txdb_info)
  introns <- unlist(introns_list, use.names = FALSE)
  
  feature_type <- rep("intergenic", length(gr))
  
  overlap_promoter <- GenomicRanges::findOverlaps(gr, promoters)
  if (length(overlap_promoter) > 0) {
    hits <- queryHits(overlap_promoter)
    for (i in hits) {
      feature_type[i] <- "promoter"
    }
  }
  
  overlap_exon <- GenomicRanges::findOverlaps(gr, exons)
  if (length(overlap_exon) > 0) {
    hits <- queryHits(overlap_exon)
    for (i in hits) {
      feature_type[i] <- "exon"
    }
  }
  
  overlap_intron <- GenomicRanges::findOverlaps(gr, introns)
  if (length(overlap_intron) > 0) {
    hits <- queryHits(overlap_intron)
    for (i in hits) {
      feature_type[i] <- "intron"
    }
  }
  
  mcols(gr)$feature_type <- feature_type
  return(gr)
}


#' Find the closest gene to each genomic region
#'
#' For each genomic interval, this function identifies the nearest gene and calculates the distance
#' from the interval to the gene
#'
#' @param gr GRanges object of genomic regions
#' @param txdb_info A TxDb object with gene coordinates on a certain reference genome
#' @return A GRanges object with extra metadata columns 'nearest_gene' and 'distance_to_gene'
#' @examples
#' library(GenomicRanges)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges(100000, 101000))
#' nearest_gene(gr, txdb)
#' @importFrom GenomicFeatures genes
#' @export
nearest_gene <- function(gr, txdb_info) {
  if (!("feature_type" %in% names(GenomicRanges::mcols(gr)))) {
    stop("Input GRanges object must contain columns added by annotate_regions().")
  }
  
  all_genes_list <- GenomicFeatures::genes(txdb_info, single.strand.genes.only = FALSE)
  all_genes <- unlist(all_genes_list)
  
  hits <- GenomicRanges::distanceToNearest(gr, all_genes)
  
  nearest_names <- names(all_genes)[S4Vectors::subjectHits(hits)]
  distances <- S4Vectors::mcols(hits)$distance
  
  GenomicRanges::mcols(gr)$nearest_gene <- nearest_names
  GenomicRanges::mcols(gr)$distance_to_gene <- distances
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
#' library(GenomicRanges)
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(100000, 101000))
#' GenomicRanges::mcols(gr)$nearest_gene <- "1"
#' # Note: Running this example requires a valid OrgDb object like org.Hs.eg.db
#' # functional_terms(gr, org.Hs.eg.db)
#' 
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom stats na.omit
#' @importFrom GenomicRanges mcols
#' @export
functional_terms <- function(gr, orgdb, go = TRUE, kegg = FALSE) {
  genes <- as.character(GenomicRanges::mcols(gr)$nearest_gene)
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



#' Perform Fisher's exact test for enrichment!
#'
#' Performs a Fisher's exact test to see if a category (like a promoter) is overrepresented among a set of genomic regions compared to a background set or not
#'
#' @param gr GRanges object with a metadata column of categories (like feature)
#' @param category A character string specifying which category to test for enrichment. Valid options can be any value present in the metadata column. So, if using the 'feature' column, valid options are: "promoter", "exon", "intron", "intergenic".
#' @param background_gr A GRanges object representing the background set of regions.
#' @return A list with elements 'odds_ratio', 'p_value', and 'table' containing the contingency table, or NULL if the category is invalid
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges(c(100,200), c(150,250)))
#' mcols(gr)$feature <- c("promoter", "exon")
#' bg_gr <- GRanges(seqnames = "chr1", ranges = IRanges(c(50,150,250), c(120,220,300)))
#' mcols(bg_gr)$feature <- c("promoter", "exon", "intron")
#' fisher_enrichment(gr, "promoter", bg_gr)
#' @import stats
#' @export
fisher_enrichment <- function(gr, category, background_gr) {
  
  gr_feature <- as.character(GenomicRanges::mcols(gr)$feature_type)
  bg_feature <- as.character(GenomicRanges::mcols(background_gr)$feature_type)
  
  valid_categories <- unique(c(gr_feature, bg_feature))
  if (!(category %in% valid_categories)) {
    message("Invalid category: '", category, "'.")
    message("Valid categories are: ", paste(valid_categories, collapse = ", "))
    return(NULL)
  }
  
  in_category <- sum(gr_feature == category)
  not_in_category <- sum(gr_feature != category)
  bg_in_category <- sum(bg_feature == category) - in_category
  bg_not_in_category <- sum(bg_feature != category) - not_in_category
  
  contingency_table <- matrix(
    c(in_category, not_in_category, bg_in_category, bg_not_in_category),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Target", "Background"),
      c("InCategory", "NotInCategory")
    )
  )
  
  test <- stats::fisher.test(contingency_table)
  
  result <- list(
    odds_ratio = test$estimate,
    p_value = test$p.value,
    table = contingency_table
  )
  
  return(result)
}
