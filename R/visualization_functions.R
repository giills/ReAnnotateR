#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal theme geom_line facet_wrap coord_flip
#' Plot Feature Composition of Genomic Intervals
#' @importFrom rlang .data
#'
#' @param gr GRanges object with a "feature" metadata column
#' @return ggplot2 object showing feature composition as a bar plot
#' @examples
#' gr <- read_bed(system.file("extdata", "example.bed", package = "ReAnnotateR"))
#' gr <- annotate_regions(gr, txdb_info)
#' plot_feature_composition(gr)
#' @export
plot_feature_composition <- function(gr) {
  if (!"feature" %in% colnames(GenomicRanges::mcols(gr))) {
    message("GRanges object must have a 'feature' metadata column")
    return(NULL)
  }
  
  feature_counts <- table(GenomicRanges::mcols(gr)$feature)
  feature_df <- data.frame(
    feature = names(feature_counts),
    count = as.numeric(feature_counts)
  )
  
  feature_df$feature <- factor(feature_df$feature, levels = c("promoter", "exon", "intron", "intergenic"))
  
  ggplot2::ggplot(feature_df, ggplot2::aes(x = feature, y = count, fill = feature)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(
      x = "Genomic Feature",
      y = "Number of Intervals",
      title = "Feature Composition of Genomic Intervals"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

#' Plot Chromosomal Density of Genomic Intervals
#'
#' @param gr GRanges object containing genomic intervals
#' @param bin_size Integer specifying bin size in base pairs (default 1e6)
#' @return ggplot2 object showing number of intervals per bin along each chromosome
#' @examples
#' gr <- read_bed(system.file("extdata", "example.bed", package = "ReAnnotateR"))
#' plot_chromosomal_density(gr)
#' @export
plot_chromosomal_density <- function(gr, bin_size = 1e6) {
  chr_list <- as.character(GenomicRanges::seqnames(gr))
  starts <- GenomicRanges::start(gr)
  
  df <- data.frame(chr = chr_list, start = starts)
  
  df$bin <- floor(df$start / bin_size) + 1
  
  density_df <- aggregate(start ~ chr + bin, data = df, FUN = length)
  colnames(density_df)[3] <- "count"
  
  ggplot2::ggplot(density_df, ggplot2::aes(x = bin * bin_size, y = count, color = chr)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ chr, scales = "free_x") +
    ggplot2::labs(
      x = "Genomic Position (bp)",
      y = "Number of Intervals",
      title = "Chromosomal Density of Genomic Intervals"
    ) +
    ggplot2::theme_minimal()
}

#' Plot Functional Enrichment Results
#'
#' @param enrich_result Object returned from functional_terms()
#' @param top_n Integer, number of top terms to plot (default 10)
#' @return ggplot2 object showing enrichment scores for top terms
#' @examples
#' gr <- read_bed(system.file("extdata", "example.bed", package = "ReAnnotateR"))
#' enrich_res <- functional_terms(gr)
#' plot_enrichment(enrich_res)
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
  
  ggplot2::ggplot(df, ggplot2::aes(x = Description, y = -log10(p.adjust))) +
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
#' @param gr GRanges object with genomic intervals
#' @param genome Reference genome build (default "hg38")
#' @param chromosomes Character vector of chromosomes to plot (default all)
#' @return Gviz plot displaying intervals on chromosome ideograms
#' @examples
#' gr <- read_bed(system.file("extdata", "example.bed", package = "ReAnnotateR"))
#' plot_ideogram(gr)
#' @export
plot_ideogram <- function(gr, genome = "hg38", chromosomes = NULL) {
  if (is.null(chromosomes)) {
    chromosomes <- unique(as.character(GenomicRanges::seqnames(gr)))
  }
  
  tracks <- list()
  
  for (chr in chromosomes) {
    chr_gr <- gr[GenomicRanges::seqnames(gr) == chr]
    if (length(chr_gr) == 0) next
    
    ideogram <- Gviz::IdeogramTrack(genome = genome, chromosome = chr)
    data_track <- Gviz::AnnotationTrack(chr_gr, name = "Intervals", genome = genome, chromosome = chr)
    
    Gviz::plotTracks(list(ideogram, data_track), from = min(start(chr_gr)), to = max(end(chr_gr)), chromosome = chr)
  }
}

#' Plot Regulatory Regions Overlaps
#'
#' @param gr GRanges object with your intervals
#' @param regulatory_gr GRanges object with regulatory regions
#' @return ggplot2 object showing which intervals overlap regulatory regions
#' @examples
#' gr <- read_bed(system.file("extdata", "example.bed", package = "ReAnnotateR"))
#' plot_regulatory_regions(gr, regulatory_gr)
#' @export
plot_regulatory_regions <- function(gr, regulatory_gr) {
  overlaps <- GenomicRanges::findOverlaps(gr, regulatory_gr)
  overlap_flag <- rep("No overlap", length(gr))
  if (length(overlaps) > 0) {
    overlap_flag[queryHits(overlaps)] <- "Overlap"
  }
  
  df <- data.frame(
    Interval = seq_along(gr),
    Overlap = overlap_flag
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = Interval, fill = Overlap)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Interval", y = "Count", title = "Regulatory Region Overlaps") +
    ggplot2::theme_minimal()
}

#' Plot Missing Annotations
#'
#' @param gr GRanges object with a "feature" metadata column
#' @return ggplot2 barplot showing counts of annotated vs. unannotated intervals
#' @examples
#' gr <- read_bed(system.file("extdata", "example.bed", package = "ReAnnotateR"))
#' gr$feature <- c("promoter", "intergenic", NA, "exon")
#' plot_missing_annotations(gr)
#' @export
plot_missing_annotations <- function(gr) {
  feature_status <- ifelse(is.na(GenomicRanges::mcols(gr)$feature), "Missing", "Annotated")
  df <- data.frame(Status = feature_status)
  
  ggplot2::ggplot(df, ggplot2::aes(x = Status, fill = Status)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "", y = "Count", title = "Annotated vs Missing Intervals") +
    ggplot2::theme_minimal()
}

#' Plot Annotation Summary
#'
#' @param gr GRanges object with a "feature" metadata column
#' @return ggplot2 barplot showing counts of intervals per feature type
#' @examples
#' gr <- read_bed(system.file("extdata", "example.bed", package = "ReAnnotateR"))
#' gr$feature <- c("promoter", "exon", "intron", "intergenic", NA)
#' plot_annotation_summary(gr)
#' @export
plot_annotation_summary <- function(gr) {
  feature_vec <- GenomicRanges::mcols(gr)$feature
  feature_vec[is.na(feature_vec)] <- "Missing"
  
  df <- data.frame(Feature = feature_vec)
  
  ggplot2::ggplot(df, ggplot2::aes(x = Feature, fill = Feature)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Feature Type", y = "Count", title = "Annotation Summary") +
    ggplot2::theme_minimal()
}