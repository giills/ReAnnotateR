#' Visualize Genomic Feature Type Distribution
#'
#' Creates a bar plot showing the count of genomic intervals classified into
#' each feature type: promoters, exons, introns, and intergenic regions. This
#' provides a quick overview of which types of genomic elements your intervals
#' predominantly overlap. For example, if most intervals fall in promoters, your
#' dataset may be enriched for regulatory regions near transcription start sites.
#' Returns a ggplot2 object that can be further customized or saved.
#' 
#' @param gr GRanges object with a "feature_type" metadata column
#' @return ggplot2 object showing feature composition as a bar plot, or NULL if column missing
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' gr_annotated <- annotate_regions(gr, txdb)
#' 
#' # Display plot
#' p <- plot_feature_composition(gr_annotated)
#' print(p)
#' 
#' # Customize plot
#' p + ggplot2::labs(title = "Custom Title") +
#'     ggplot2::theme_bw()
#' 
#' # Save plot
#' ggplot2::ggsave("features.png", plot = p, width = 7, height = 5)
#' }
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal theme
#' @importFrom S4Vectors mcols
#' @importFrom rlang .data
#' @export
plot_feature_composition <- function(gr) {
  if (!"feature_type" %in% colnames(S4Vectors::mcols(gr))) {
    message("GRanges object must have a 'feature_type' metadata column")
    return(invisible(NULL))
  }
  
  feature_counts <- table(S4Vectors::mcols(gr)$feature_type)
  feature_df <- data.frame(
    feature = names(feature_counts),
    count = as.numeric(feature_counts)
  )
  
  feature_df$feature <- factor(feature_df$feature, levels = c("promoter", "exon", "intron", "intergenic"))
  
  p <- ggplot2::ggplot(feature_df, ggplot2::aes(x = .data$feature, y = .data$count, fill = .data$feature)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(
      x = "Genomic Feature",
      y = "Number of Intervals",
      title = "Feature Composition of Genomic Intervals"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  
  return(p)
}

#' Visualize Distribution of Genomic Intervals Across Chromosomes
#'
#' Creates a multi-panel plot showing how genomic intervals are distributed
#' along each chromosome. The genome is divided into bins (default 1Mb), and the
#' number of intervals in each bin is plotted. This visualization helps identify
#' genomic regions with high clustering of intervals and can reveal patterns like
#' enrichment near telomeres or centromeres. Each chromosome is shown in a 
#' separate panel for easy comparison. Returns a ggplot2 object that can be 
#' further customized or saved.
#' 
#' @param gr GRanges object containing genomic intervals
#' @param bin_size Integer specifying bin size in base pairs (default 1e6)
#' @return ggplot2 object showing number of intervals per bin along each chromosome
#' @examples
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' if (file.exists(bed_file)) {
#'   gr <- read_bed(bed_file)
#'   
#'   # Display plot
#'   p <- plot_chromosomal_density(gr)
#'   print(p)
#'   
#'   # Use smaller bins
#'   p2 <- plot_chromosomal_density(gr, bin_size = 5e6)
#'   print(p2)
#'   
#'   # Save plot
#'   ggplot2::ggsave("density.png", plot = p, width = 10, height = 6)
#' }
#' @importFrom ggplot2 ggplot aes geom_col labs theme_minimal facet_wrap theme element_text scale_x_continuous
#' @importFrom GenomicRanges seqnames start
#' @importFrom stats aggregate
#' @importFrom rlang .data
#' @export
plot_chromosomal_density <- function(gr, bin_size = 1e6) {
  chr_list <- as.character(GenomicRanges::seqnames(gr))
  starts <- GenomicRanges::start(gr)
  
  df <- data.frame(chr = chr_list, start = starts, stringsAsFactors = FALSE)
  
  # Assign each interval to a bin (use bin midpoint for better visualization)
  df$bin <- floor(df$start / bin_size) * bin_size + (bin_size / 2)
  
  # Count intervals per bin per chromosome
  density_df <- aggregate(start ~ chr + bin, data = df, FUN = length)
  colnames(density_df)[3] <- "count"
  
  # Ensure bin is numeric (not scientific notation issues)
  density_df$bin <- as.numeric(density_df$bin)
  
  # Create the plot
  p <- ggplot2::ggplot(density_df, ggplot2::aes(x = .data$bin, y = .data$count)) +
    ggplot2::geom_col(fill = "steelblue", color = "steelblue", width = bin_size * 0.8) +
    ggplot2::facet_wrap(~ chr, scales = "free_x", ncol = 2) +
    ggplot2::labs(
      x = "Genomic Position (Mb)",
      y = "Number of Intervals",
      title = "Chromosomal Density of Genomic Intervals"
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) round(x / 1e6, 0),
      expand = c(0.02, 0)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      panel.grid.major.x = ggplot2::element_line(color = "grey90"),
      panel.grid.minor.x = ggplot2::element_blank()
    )
  
  return(p)
}

#' Visualize Functional Enrichment Results
#'
#' Creates a horizontal bar plot displaying the top enriched biological terms
#' from Gene Ontology (GO) or KEGG pathway analysis. The plot shows the 
#' negative log10 of adjusted p-values, making it easy to identify the most
#' significantly enriched terms. Longer bars indicate more significant enrichment.
#' Returns a ggplot2 object that can be further customized or saved.
#'
#' @param enrich_result Data frame or enrichment result object returned from functional_terms()
#' @param top_n Integer, number of top terms to plot (default 10)
#' @return ggplot2 object showing enrichment scores for top terms, or NULL if no results
#' @examples
#' \dontrun{
#' # Example 1: Basic usage with actual enrichment results
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#' 
#' # Load BED file and annotate
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' gr_annotated <- annotate_regions(gr, txdb)
#' gr_with_genes <- nearest_gene(gr_annotated, txdb)
#' 
#' # Perform enrichment analysis
#' enrich_res <- functional_terms(gr_with_genes, org.Hs.eg.db, go = TRUE)
#' 
#' # Plot GO enrichment results
#' if (!is.null(enrich_res$go) && nrow(enrich_res$go) > 0) {
#'   p <- plot_enrichment(enrich_res$go, top_n = 10)
#'   print(p)
#'   
#'   # Save plot
#'   ggplot2::ggsave("enrichment.png", plot = p, width = 8, height = 6)
#' } else {
#'   message("No significant enrichment found")
#' }
#' }
#' 
#' # Example 2: Simulated enrichment results for demonstration
#' # This shows what the plot looks like when enrichment is detected
#' simulated_results <- data.frame(
#'   Description = c("DNA repair", "cell cycle", "apoptosis", 
#'                   "signal transduction", "transcription"),
#'   p.adjust = c(0.001, 0.005, 0.01, 0.02, 0.04),
#'   Count = c(25, 20, 15, 12, 10),
#'   GeneRatio = c("25/100", "20/100", "15/100", "12/100", "10/100")
#' )
#' p <- plot_enrichment(simulated_results, top_n = 5)
#' print(p)
#' 
#' # Example 3: Customize the plot
#' p + ggplot2::labs(title = "My Custom Title") + 
#'     ggplot2::theme_bw() +
#'     ggplot2::theme(text = ggplot2::element_text(size = 12))
#' 
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip labs theme_minimal
#' @importFrom rlang .data
#' @export
plot_enrichment <- function(enrich_result, top_n = 10) {
  if (is.null(enrich_result) || length(enrich_result) == 0) {
    message("No enrichment results to plot.")
    return(invisible(NULL))
  }
  
  df <- as.data.frame(enrich_result)
  
  if (nrow(df) == 0) {
    message("No enrichment results to plot.")
    return(invisible(NULL))
  }
  
  df <- df[order(df$p.adjust), ]
  df <- df[1:min(top_n, nrow(df)), ]
  
  df$Description <- factor(df$Description, levels = rev(df$Description))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Description, y = -log10(.data$p.adjust))) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "",
      y = "-log10 Adjusted P-value",
      title = "Top Enriched Biological Functions"
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}

#' Plot Chromosome Ideogram with Intervals
#'
#' Creates an ideogram visualization showing intervals on chromosomes using Gviz.
#' Note: This function creates plots as a side effect and returns NULL invisibly.
#' To save plots, use a graphics device (e.g., png(), pdf()) before calling.
#'
#' @param gr GRanges object with genomic intervals
#' @param genome Reference genome build (default "hg38")
#' @param chromosomes Character vector of chromosomes to plot (default all)
#' @return NULL invisibly (plots are created as side effect)
#' @examples
#' \dontrun{
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' 
#' # Display ideograms for all chromosomes
#' plot_ideogram(gr)
#' 
#' # Plot specific chromosomes
#' plot_ideogram(gr, chromosomes = c("chr1", "chr2"))
#' 
#' # Save to file
#' png("ideogram.png", width = 800, height = 600)
#' plot_ideogram(gr, chromosomes = "chr1")
#' dev.off()
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
#' Returns a ggplot2 object that can be further customized or saved.
#' 
#' @param gr GRanges object with your intervals
#' @param regulatory_gr GRanges object with regulatory regions
#' @return ggplot2 object showing which intervals overlap regulatory regions
#' @examples
#' \dontrun{
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' gr <- read_bed(bed_file)
#' 
#' # Create example regulatory regions
#' regulatory_gr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = c(10500, 50500), end = c(10700, 50700))
#' )
#' 
#' # Display plot
#' p <- plot_regulatory_regions(gr, regulatory_gr)
#' print(p)
#' 
#' # Save plot
#' ggplot2::ggsave("regulatory_overlap.png", plot = p, width = 7, height = 5)
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
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Interval, fill = .data$Overlap)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Interval", y = "Count", title = "Regulatory Region Overlaps") +
    ggplot2::theme_minimal()
  
  return(p)
}

#' Assess Annotation Completeness
#'
#' Creates a simple two-bar plot comparing the number of successfully annotated
#' genomic intervals versus those that failed annotation (have NA feature types).
#' A high proportion of missing annotations may indicate problems with genome
#' version mismatch, incorrect coordinate systems, or intervals in regions not
#' covered by the annotation database (e.g., unplaced contigs).
#' Returns a ggplot2 object that can be further customized or saved.
#' 
#' @param gr GRanges object with a "feature_type" metadata column
#' @return ggplot2 barplot showing counts of annotated vs. unannotated intervals
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(c(100, 200, 300, 400), 
#'                                                             c(150, 250, 350, 450)))
#' mcols(gr)$feature_type <- c("promoter", "intergenic", NA, "exon")
#' 
#' # Display plot
#' p <- plot_missing_annotations(gr)
#' print(p)
#' 
#' # Customize and save
#' p + ggplot2::labs(title = "Quality Check: Annotation Coverage") +
#'     ggplot2::theme_bw()
#' ggplot2::ggsave("annotation_qc.png", plot = p, width = 6, height = 4)
#' 
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal
#' @importFrom S4Vectors mcols
#' @importFrom rlang .data
#' @export
plot_missing_annotations <- function(gr) {
  feature_status <- ifelse(is.na(S4Vectors::mcols(gr)$feature_type), "Missing", "Annotated")
  df <- data.frame(Status = feature_status)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Status, fill = .data$Status)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "", y = "Count", title = "Annotated vs Missing Intervals") +
    ggplot2::theme_minimal()
  
  return(p)
}

#' Summarize Annotation Completeness and Feature Distribution
#'
#' Creates a bar plot showing the count of genomic intervals in each feature
#' category, including any regions that failed to annotate (marked as "Missing").
#' This differs from plot_feature_composition() by highlighting annotation
#' completeness. Use this to quality-check your annotation results and identify
#' if any intervals could not be classified.
#' Returns a ggplot2 object that can be further customized or saved.
#' 
#' @param gr GRanges object with a "feature_type" metadata column
#' @return ggplot2 barplot showing counts of intervals per feature type
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(c(100, 200, 300, 400, 500), 
#'                                                             c(150, 250, 350, 450, 550)))
#' mcols(gr)$feature_type <- c("promoter", "exon", "intron", "intergenic", NA)
#' 
#' # Display plot
#' p <- plot_annotation_summary(gr)
#' print(p)
#' 
#' # Customize colors and theme
#' p + ggplot2::scale_fill_brewer(palette = "Set2") +
#'     ggplot2::theme_classic()
#' 
#' # Save plot
#' ggplot2::ggsave("annotation_summary.png", plot = p, width = 7, height = 5)
#' 
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal
#' @importFrom S4Vectors mcols
#' @importFrom rlang .data
#' @export
plot_annotation_summary <- function(gr) {
  feature_vec <- S4Vectors::mcols(gr)$feature_type
  feature_vec[is.na(feature_vec)] <- "Missing"
  
  df <- data.frame(Feature = feature_vec)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Feature, fill = .data$Feature)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Feature Type", y = "Count", title = "Annotation Summary") +
    ggplot2::theme_minimal()
  
  return(p)
}