#' Launch ReAnnotateR Shiny App
#'
#' Launches the interactive Shiny application for genomic region annotation
#'
#' @return Launches the Shiny application
#' @examples
#' \dontrun{
#' runReAnnotateR()
#' }
#' @export
runReAnnotateR <- function() {
  appDir <- system.file("shiny-scripts", package = "ReAnnotateR")
  if (appDir == "") {
    stop("Could not find Shiny app directory. Try re-installing ReAnnotateR.", 
         call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}