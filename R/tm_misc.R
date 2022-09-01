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
