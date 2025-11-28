#' Visualize Genomic Feature Type Distribution
#'
#' Creates a bar plot showing the count of genomic intervals classified into
#' each feature type: promoters, exons, introns, and intergenic regions. This
#' provides a quick overview of which types of genomic elements your intervals
#' predominantly overlap. For example, if most intervals fall in promoters, your
#' dataset may be enriched for regulatory regions near transcription start sites.
#' 
#' @param gr GRanges object with a "feature_type" metadata column
#' @return ggplot2 object showing feature composition as a bar plot
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' gr_annotated <- annotate_regions(gr, txdb)
#' plot_feature_composition(gr_annotated)
#' }
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal theme
#' @importFrom S4Vectors mcols
#' @importFrom rlang .data
#' @export
plot_feature_composition <- function(gr) {
  if (!"feature_type" %in% colnames(S4Vectors::mcols(gr))) {
    message("GRanges object must have a 'feature_type' metadata column")
    return(NULL)
  }
  
  feature_counts <- table(S4Vectors::mcols(gr)$feature_type)
  feature_df <- data.frame(
    feature = names(feature_counts),
    count = as.numeric(feature_counts)
  )
  
  feature_df$feature <- factor(feature_df$feature, levels = c("promoter", "exon", "intron", "intergenic"))
  
  ggplot2::ggplot(feature_df, ggplot2::aes(x = .data$feature, y = .data$count, fill = .data$feature)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(
      x = "Genomic Feature",
      y = "Number of Intervals",
      title = "Feature Composition of Genomic Intervals"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

#' Visualize Distribution of Genomic Intervals Across Chromosomes
#'
#' Creates a multi-panel line plot showing how genomic intervals are distributed
#' along each chromosome. The genome is divided into bins (default 1Mb), and the
#' number of intervals in each bin is plotted. This visualization helps identify
#' genomic regions with high clustering of intervals and can reveal patterns like
#' enrichment near telomeres or centromeres. Each chromosome is shown in a 
#' separate panel for easy comparison.
#' 
#' @param gr GRanges object containing genomic intervals
#' @param bin_size Integer specifying bin size in base pairs (default 1e6)
#' @return ggplot2 object showing number of intervals per bin along each chromosome
#' @examples
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' if (file.exists(bed_file)) {
#'   gr <- read_bed(bed_file)
#'   plot_chromosomal_density(gr)
#' }
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal facet_wrap
#' @importFrom GenomicRanges seqnames start
#' @importFrom stats aggregate
#' @importFrom rlang .data
#' @export
plot_chromosomal_density <- function(gr, bin_size = 1e6) {
  chr_list <- as.character(GenomicRanges::seqnames(gr))
  starts <- GenomicRanges::start(gr)
  
  df <- data.frame(chr = chr_list, start = starts)
  
  df$bin <- floor(df$start / bin_size) + 1
  
  density_df <- aggregate(start ~ chr + bin, data = df, FUN = length)
  colnames(density_df)[3] <- "count"
  
  ggplot2::ggplot(density_df, ggplot2::aes(x = .data$bin * bin_size, y = .data$count, color = .data$chr)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ chr, scales = "free_x") +
    ggplot2::labs(
      x = "Genomic Position (bp)",
      y = "Number of Intervals",
      title = "Chromosomal Density of Genomic Intervals"
    ) +
    ggplot2::theme_minimal()
}

#' Visualize Functional Enrichment Results
#'
#' Creates a horizontal bar plot displaying the top enriched biological terms
#' from Gene Ontology (GO) or KEGG pathway analysis. The plot shows the 
#' negative log10 of adjusted p-values, making it easy to identify the most
#' significantly enriched terms. Longer bars indicate more significant enrichment.
#'
#' @param enrich_result Data frame or enrichment result object returned from functional_terms()
#' @param top_n Integer, number of top terms to plot (default 10)
#' @return ggplot2 object showing enrichment scores for top terms
#' @examples
#' \dontrun{
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' enrich_res <- functional_terms(gr, org.Hs.eg.db)
#' plot_enrichment(enrich_res$go)
#' }
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip labs theme_minimal
#' @importFrom rlang .data
#' @export
plot_enrichment <- function(enrich_result, top_n = 10) {
  if (is.null(enrich_result) || length(enrich_result) == 0) {
    message("No enrichment results to plot.")
    return(NULL)
  }
  
  df <- as.data.frame(enrich_result)
  
  if (nrow(df) == 0) {
    message("No enrichment results to plot.")
    return(NULL)
  }
  
  df <- df[order(df$p.adjust), ]
  df <- df[1:min(top_n, nrow(df)), ]
  
  df$Description <- factor(df$Description, levels = rev(df$Description))
  
  ggplot2::ggplot(df, ggplot2::aes(x = .data$Description, y = -log10(.data$p.adjust))) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "",
      y = "-log10 Adjusted P-value",
      title = "Top Enriched Biological Functions"
    ) +
    ggplot2::theme_minimal()
}

#' Plot Chromosome Ideogram with Intervals
#'
#' Creates an ideogram visualization showing intervals on chromosomes
#'
#' @param gr GRanges object with genomic intervals
#' @param genome Reference genome build (default "hg38")
#' @param chromosomes Character vector of chromosomes to plot (default all)
#' @return NULL (plots are created as side effect)
#' @examples
#' \dontrun{
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' plot_ideogram(gr)
#' }
#' @importFrom GenomicRanges seqnames start end
#' @export
plot_ideogram <- function(gr, genome = "hg38", chromosomes = NULL) {
  if (!requireNamespace("Gviz", quietly = TRUE)) {
    stop("Package 'Gviz' is required for this function. Please install it.")
  }
  
  if (is.null(chromosomes)) {
    chromosomes <- unique(as.character(GenomicRanges::seqnames(gr)))
  }
  
  for (chr in chromosomes) {
    chr_gr <- gr[GenomicRanges::seqnames(gr) == chr]
    if (length(chr_gr) == 0) next
    
    ideogram <- Gviz::IdeogramTrack(genome = genome, chromosome = chr)
    data_track <- Gviz::AnnotationTrack(chr_gr, name = "Intervals", genome = genome, chromosome = chr)
    
    Gviz::plotTracks(list(ideogram, data_track), 
                     from = min(GenomicRanges::start(chr_gr)), 
                     to = max(GenomicRanges::end(chr_gr)), 
                     chromosome = chr)
  }
  
  return(invisible(NULL))
}

#' Visualize Overlap with Regulatory Regions
#'
#' Creates a bar plot showing how many of your genomic intervals overlap with
#' a provided set of regulatory regions (e.g., enhancers, silencers, CTCF binding
#' sites). Requires a second GRanges object containing the regulatory regions of
#' interest. This is useful for determining if your dataset is enriched for
#' specific types of regulatory elements beyond just promoters.
#' 
#' @param gr GRanges object with your intervals
#' @param regulatory_gr GRanges object with regulatory regions
#' @return ggplot2 object showing which intervals overlap regulatory regions
#' @examples
#' \dontrun{
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' # regulatory_gr would be loaded from a data source
#' # plot_regulatory_regions(gr, regulatory_gr)
#' }
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom rlang .data
#' @export
plot_regulatory_regions <- function(gr, regulatory_gr) {
  overlaps <- GenomicRanges::findOverlaps(gr, regulatory_gr)
  overlap_flag <- rep("No overlap", length(gr))
  if (length(overlaps) > 0) {
    overlap_flag[S4Vectors::queryHits(overlaps)] <- "Overlap"
  }
  
  df <- data.frame(
    Interval = seq_along(gr),
    Overlap = overlap_flag
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = .data$Interval, fill = .data$Overlap)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Interval", y = "Count", title = "Regulatory Region Overlaps") +
    ggplot2::theme_minimal()
}

#' Assess Annotation Completeness
#'
#' Creates a simple two-bar plot comparing the number of successfully annotated
#' genomic intervals versus those that failed annotation (have NA feature types).
#' A high proportion of missing annotations may indicate problems with genome
#' version mismatch, incorrect coordinate systems, or intervals in regions not
#' covered by the annotation database (e.g., unplaced contigs).
#' 
#' @param gr GRanges object with a "feature_type" metadata column
#' @return ggplot2 barplot showing counts of annotated vs. unannotated intervals
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(c(100, 200, 300, 400), 
#'                                                             c(150, 250, 350, 450)))
#' mcols(gr)$feature_type <- c("promoter", "intergenic", NA, "exon")
#' plot_missing_annotations(gr)
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal
#' @importFrom S4Vectors mcols
#' @importFrom rlang .data
#' @export
plot_missing_annotations <- function(gr) {
  feature_status <- ifelse(is.na(S4Vectors::mcols(gr)$feature_type), "Missing", "Annotated")
  df <- data.frame(Status = feature_status)
  
  ggplot2::ggplot(df, ggplot2::aes(x = .data$Status, fill = .data$Status)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "", y = "Count", title = "Annotated vs Missing Intervals") +
    ggplot2::theme_minimal()
}

#' Summarize Annotation Completeness and Feature Distribution
#'
#' Creates a bar plot showing the count of genomic intervals in each feature
#' category, including any regions that failed to annotate (marked as "Missing").
#' This differs from plot_feature_composition() by highlighting annotation
#' completeness. Use this to quality-check your annotation results and identify
#' if any intervals could not be classified.
#' 
#' @param gr GRanges object with a "feature_type" metadata column
#' @return ggplot2 barplot showing counts of intervals per feature type
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(c(100, 200, 300, 400, 500), 
#'                                                             c(150, 250, 350, 450, 550)))
#' mcols(gr)$feature_type <- c("promoter", "exon", "intron", "intergenic", NA)
#' plot_annotation_summary(gr)
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal
#' @importFrom S4Vectors mcols
#' @importFrom rlang .data
#' @export
plot_annotation_summary <- function(gr) {
  feature_vec <- S4Vectors::mcols(gr)$feature_type
  feature_vec[is.na(feature_vec)] <- "Missing"
  
  df <- data.frame(Feature = feature_vec)
  
  ggplot2::ggplot(df, ggplot2::aes(x = .data$Feature, fill = .data$Feature)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Feature Type", y = "Count", title = "Annotation Summary") +
    ggplot2::theme_minimal()
}