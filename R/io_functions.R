#' Read BED file into GRanges object
#'
#' Loads a BED file and converts it to a GRanges object.
#' Checks for invalid chromosome start and end positions.
#'
#' @param file path to a BED file
#' @return `GRanges` object containing the intervals from the BED file
#' @examples
#' bed_path <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' if (file.exists(bed_path)) {
#'   gr <- read_bed(bed_path)
#' }
#' @importFrom utils read.table
#' @importFrom GenomicRanges GRanges strand seqnames start end
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<-
#' @export
read_bed <- function(file) {
  if (!is.character(file) || length(file) != 1 || is.na(file)) {
    return(NULL)
  }
  
  if (!file.exists(file)) {
    message("This file does not exist")
    return(NULL)
  }
  
  tryCatch({
    bed_data <- utils::read.table(file, header = FALSE, stringsAsFactors = FALSE)
  }, error = function(e) {
    message("Error reading BED file: ", e$message)
    return(NULL)
  })
  
  n_cols <- ncol(bed_data)
  
  if (n_cols < 3) {
    message("BED file must have at least 3 columns! (chr, start, end)")
    return(NULL)
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = bed_data[, 1],
    ranges = IRanges::IRanges(start = as.numeric(bed_data[, 2]) + 1,
                              end = as.numeric(bed_data[, 3]))
  )
  
  if (n_cols >= 4) {
    S4Vectors::mcols(gr)$name <- bed_data[, 4]
  }
  if (n_cols >= 5) {
    S4Vectors::mcols(gr)$score <- as.numeric(bed_data[, 5])
  }
  if (n_cols >= 6) {
    GenomicRanges::strand(gr) <- bed_data[, 6]
  }
  
  for (i in seq_along(GenomicRanges::start(gr))) {
    if (GenomicRanges::start(gr)[i] < 1) {
      message("All start positions must be >= 1")
      return(NULL)
    }
  }
  
  for (i in seq_along(GenomicRanges::start(gr))) {
    if (GenomicRanges::end(gr)[i] < GenomicRanges::start(gr)[i]) {
      message("All end positions must be >= start positions")
      return(NULL)
    }
  }
  
  return(gr)
}

#' Convert genomic coordinates between assemblies
#'
#' Uses a chain file to lift genomic coordinates from one reference assembly to another
#'
#' @param gr a `GRanges` object
#' @param chain_file A path to the chain file
#' @return A `GRanges` object with coordinates lifted over to the new assembly
#' @examples
#' bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' chain_file <- system.file("extdata", "hg19ToHg38.over.chain", package = "ReAnnotateR")
#' if (file.exists(bed_file) && file.exists(chain_file)) {
#'   gr <- read_bed(bed_file)
#'   gr_new <- convert_coordinates(gr, chain_file)
#' }
#' @importFrom rtracklayer import.chain liftOver
#' @export
convert_coordinates <- function(gr, chain_file) {
  if (!inherits(gr, "GRanges")) {
    message("gr must be a GRanges object")
    return(NULL)
  }
  
  if (!file.exists(chain_file)) {
    message("Chain file does not exist.")
    return(NULL)
  }
  
  chain <- rtracklayer::import.chain(chain_file)
  gr_new_list <- rtracklayer::liftOver(gr, chain)
  gr_new <- unlist(gr_new_list)
  
  return(gr_new)
}

#' Export annotated or enriched results
#'
#' Saves annotated genomic intervals or enrichment results to a file.
#' Supports writing GRanges, data frames, or lists in TSV, CSV, or BED format
#'
#' @param object GRanges object, data frame, or enrichment results list
#' @param file Path to output file
#' @param format Output format: "tsv", "csv", or "bed" (default "tsv")
#' @return Invisible NULL
#' @examples
#' bed_path <- system.file("extdata", "example.bed", package = "ReAnnotateR")
#' if (file.exists(bed_path)) {
#'   gr <- read_bed(bed_path)
#'   temp_file <- tempfile(fileext = ".tsv")
#'   export_results(gr, temp_file, format = "tsv")
#' }
#' @importFrom methods is
#' @importFrom utils write.table write.csv
#' @importFrom rtracklayer export
#' @importFrom GenomicRanges seqnames start end strand
#' @importFrom S4Vectors mcols
#' @export
export_results <- function(object, file, format = "tsv") {
  
  if (!format %in% c("tsv", "csv", "bed")) {
    message("Format must be 'tsv', 'csv', or 'bed'")
    return(invisible(NULL))
  }
  
  if (!is.character(file) || length(file) != 1) {
    message("File path must be a single character string")
    return(invisible(NULL))
  }
  
  if (methods::is(object, "GRanges")) {
    if (format == "bed") {
      tryCatch({
        rtracklayer::export(object, file, format = "BED")
        message("Successfully exported GRanges to BED format: ", file)
      }, error = function(e) {
        message("Error exporting to BED format: ", e$message)
        return(invisible(NULL))
      })
    } else {
      df <- data.frame(
        chr = as.character(GenomicRanges::seqnames(object)),
        start = GenomicRanges::start(object),
        end = GenomicRanges::end(object),
        strand = as.character(GenomicRanges::strand(object))
      )
      
      if (ncol(S4Vectors::mcols(object)) > 0) {
        mcols_df <- as.data.frame(S4Vectors::mcols(object))
        df <- cbind(df, mcols_df)
      }
      
      if (format == "tsv") {
        utils::write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
        message("Successfully exported GRanges to TSV: ", file)
      } else {
        utils::write.csv(df, file, row.names = FALSE)
        message("Successfully exported GRanges to CSV: ", file)
      }
    }
    
  } else if (is.data.frame(object)) {
    if (format == "bed") {
      message("BED format is only supported for GRanges objects. Using TSV instead.")
      format <- "tsv"
    }
    
    if (format == "tsv") {
      utils::write.table(object, file, sep = "\t", row.names = FALSE, quote = FALSE)
      message("Successfully exported data frame to TSV: ", file)
    } else {
      utils::write.csv(object, file, row.names = FALSE)
      message("Successfully exported data frame to CSV: ", file)
    }
    
  } else if (is.list(object)) {
    if (format == "bed") {
      message("BED format is not supported for enrichment results. Using TSV instead.")
      format <- "tsv"
    }
    
    file_base <- sub("\\.[^.]*$", "", file)
    file_ext <- if (format == "tsv") ".tsv" else ".csv"
    
    exported_any <- FALSE
    
    if ("go" %in% names(object) && !is.null(object$go)) {
      go_file <- paste0(file_base, "_go", file_ext)
      if (format == "tsv") {
        utils::write.table(object$go, go_file, sep = "\t", row.names = FALSE, quote = FALSE)
      } else {
        utils::write.csv(object$go, go_file, row.names = FALSE)
      }
      message("Successfully exported GO results to: ", go_file)
      exported_any <- TRUE
    }
    
    if ("kegg" %in% names(object) && !is.null(object$kegg)) {
      kegg_file <- paste0(file_base, "_kegg", file_ext)
      if (format == "tsv") {
        utils::write.table(object$kegg, kegg_file, sep = "\t", row.names = FALSE, quote = FALSE)
      } else {
        utils::write.csv(object$kegg, kegg_file, row.names = FALSE)
      }
      message("Successfully exported KEGG results to: ", kegg_file)
      exported_any <- TRUE
    }
    
    if (!exported_any) {
      message("No GO or KEGG results found in the list to export")
      return(invisible(NULL))
    }
    
  } else {
    message("Input must be a GRanges object, data frame, or list of enrichment results")
    return(invisible(NULL))
  }
  
  return(invisible(NULL))
}