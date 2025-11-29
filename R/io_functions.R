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
  
  bed_data <- tryCatch({
    utils::read.table(file, header = FALSE, stringsAsFactors = FALSE)
  }, error = function(e) {
    message("Error reading BED file: ", e$message)
    return(NULL)
  })
  
  # Check if read failed
  if (is.null(bed_data)) {
    return(NULL)
  }
  
  # Check if file is empty
  if (nrow(bed_data) == 0) {
    message("BED file is empty")
    return(NULL)
  }
  
  n_cols <- ncol(bed_data)
  
  if (n_cols < 3) {
    message("BED file must have at least 3 columns! (chr, start, end)")
    return(NULL)
  }
  
  # Check for non-numeric coordinates
  start_coords <- suppressWarnings(as.numeric(bed_data[, 2]))
  end_coords <- suppressWarnings(as.numeric(bed_data[, 3]))
  
  if (any(is.na(start_coords)) || any(is.na(end_coords))) {
    message("Start and end coordinates must be numeric")
    return(NULL)
  }
  
  # Convert from 0-based to 1-based coordinates
  start_coords <- start_coords + 1
  
  # Check for invalid coordinates before creating GRanges
  if (any(start_coords < 1)) {
    message("All start positions must be >= 1")
    return(NULL)
  }
  
  if (any(end_coords < start_coords)) {
    message("All end positions must be >= start positions")
    return(NULL)
  }
  
  # Create GRanges object with validated coordinates
  gr <- tryCatch({
    GenomicRanges::GRanges(
      seqnames = bed_data[, 1],
      ranges = IRanges::IRanges(start = start_coords, end = end_coords)
    )
  }, error = function(e) {
    message("Error creating GRanges object: ", e$message)
    return(NULL)
  })
  
  if (is.null(gr)) {
    return(NULL)
  }
  
  # Add optional columns
  if (n_cols >= 4) {
    S4Vectors::mcols(gr)$name <- bed_data[, 4]
  }
  if (n_cols >= 5) {
    S4Vectors::mcols(gr)$score <- as.numeric(bed_data[, 5])
  }
  if (n_cols >= 6) {
    # Validate strand values
    strand_values <- bed_data[, 6]
    if (!all(strand_values %in% c("+", "-", "*"))) {
      message("Invalid strand values detected. Strand must be '+', '-', or '*'")
      return(NULL)
    }
    GenomicRanges::strand(gr) <- strand_values
  }
  
  return(gr)
}

#' Convert Genomic Coordinates Between Genome Assemblies
#'
#' Performs liftOver to convert genomic coordinates from one reference genome 
#' assembly to another (e.g., from hg19 to hg38). This is essential when working 
#' with data from different genome builds or when updating analyses to newer 
#' genome versions. Uses UCSC chain files for accurate coordinate mapping. Some 
#' regions may not lift over if they are in areas that differ significantly 
#' between genome builds.
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
    stop("gr must be a GRanges object")
  }
  
  if (is.null(chain_file)) {
    stop("chain_file cannot be NULL")
  }
  
  if (!is.character(chain_file)) {
    stop("chain_file must be a character string path")
  }
  
  if (!file.exists(chain_file)) {
    message("Chain file does not exist.")
    return(NULL)
  }
  
  chain <- tryCatch({
    rtracklayer::import.chain(chain_file)
  }, error = function(e) {
    message("Error reading chain file: ", e$message)
    return(NULL)
  })
  
  if (is.null(chain)) {
    return(NULL)
  }
  
  gr_new_list <- rtracklayer::liftOver(gr, chain)
  gr_new <- unlist(gr_new_list)
  
  return(gr_new)
}

#' Export Annotated Results to File
#'
#' Saves annotated genomic intervals or functional enrichment results to disk
#' in multiple formats (TSV, CSV, or BED). For GRanges objects, exports include
#' all metadata columns (feature types, nearest genes, distances). For enrichment
#' results lists, creates separate files for GO and KEGG results. This allows you
#' to share results with collaborators, import into other tools, or create
#' publication-ready supplementary tables.
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
      
      # Handle metadata columns - convert complex types
      if (ncol(S4Vectors::mcols(object)) > 0) {
        mcols_df <- as.data.frame(S4Vectors::mcols(object))
        
        # Convert list columns to character representation
        for (col in names(mcols_df)) {
          if (is.list(mcols_df[[col]])) {
            mcols_df[[col]] <- sapply(mcols_df[[col]], function(x) {
              if (is.null(x)) return(NA_character_)
              paste(as.character(x), collapse = ",")
            })
          }
        }
        
        df <- cbind(df, mcols_df)
      }
      
      tryCatch({
        if (format == "tsv") {
          utils::write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
          message("Successfully exported GRanges to TSV: ", file)
        } else {
          utils::write.csv(df, file, row.names = FALSE)
          message("Successfully exported GRanges to CSV: ", file)
        }
      }, error = function(e) {
        message("Error writing file: ", e$message)
        return(invisible(NULL))
      })
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