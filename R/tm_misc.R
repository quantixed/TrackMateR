#' Setup Output Path
#'
#' In the working directory, we need to check if a path to a directory exists, and if not create it.
#' Ensures smooth running of the automated functions in TestTrackMateR.
#'
#' @param fpath string containing a path that may or may not exist
#' @return NULL
#' @keywords internal
setupOutputPath <- function(fpath) {
  # example "Output/Data/foo/bar"
  path_components <- unlist(strsplit(fpath, split = "/"))

  for(i in 1:length(path_components)){
    test_path <- paste(path_components[1:i], collapse = "/")
    ifelse(!dir.exists(test_path), dir.create(test_path),FALSE)
  }
}

#' Merge Data Frames for export
#'
#' Properties for individual traces (tracks) from different experiments need to be merged for export.
#' Utility function, not for export.
#'
#' @param x data frame to be coerced
#' @param y data frame to be coerced
#' @return NULL
#' @keywords internal
mergeDataFramesForExport <- function(x, y) {
  merge(x, y, by = c("trace", "dataid"), all = TRUE)
}

#' Process additional parameters
#'
#' Ensure ellipsis parameters have the default values
#' This function is used in two separate workflows, so a single function to edit makes sense.
#'
#' @param input list of ellipsis arguments
#' @return list of arguments adjusted for defaults
#' @keywords internal
processEllipsis <- function(input) {
  if (is.null(input$N)) input$N <- 3
  if (is.null(input$short)) input$short <- 8
  if (is.null(input$deltaT)) input$deltaT <- 1
  if (is.null(input$mode)) input$mode <- "ECDF"
  if (is.null(input$nPop)) input$nPop <- 2
  #if (is.null(input$init)) input$init <- NULL
  if (is.null(input$timeRes)) input$timeRes <- 1
  if (is.null(input$breaks)) input$breaks <- 100
  if (is.null(input$radius)) input$radius <- 1.5
  if (is.null(input$msdplot)) input$msdplot <- "linlin"

  return(input)
}

#' Find log2 y limits for symmetrical axis
#'
#' @param input vector of log2 values
#' @return list two limits for y-axis
#' @keywords internal
findLog2YAxisLimits <- function(x) {
  lo <- min(0.5, x, na.rm = TRUE)
  hi <- max(2, x, na.rm = TRUE)
  limit <- max(abs(floor(log2(lo))), ceiling(log2(hi)))
  lims <- c(2^(-limit), 2^limit)

  return(lims)
}
