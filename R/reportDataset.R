#' Report Dataset
#'
#' This function automates the process of making a report.
#' If more control is needed over the default parameters:
#' a user can generate MSD, jump distance and track density objects from their imported TrackMate data, and make a report manually.
#'
#' @param tmList list of TrackMate data and calibration
#' @param ... pass additional parameters to modify the defaults (N, short, deltaT, mode, nPop, init, timeRes, breaks, radius)
#' @return patchwork ggplot
#' @export
reportDataset <- function(tmList, ...) {

  if(inherits(tmList, "list")) {
    tmDF <- tmList[[1]]
    calibrationDF <- tmList[[2]]
  } else {
    cat("Function requires a list of TrackMate data and calibration data\n")
    return(NULL)
  }

  # ellipsis processing
  l <- NULL
  l <- list(...)
  l <- processEllipsis(l)

  # take the units
  units <- calibrationDF$unit[1:2]
  # calculate MSD
  msdObj <- calculateMSD(tmDF, N = l$N, short = l$short)
  msdDF <- msdObj[[1]]
  alphaDF <- msdObj[[2]]
  # jump distance calc
  jdObj <- calculateJD(dataList = tmList, deltaT = l$deltaT, nPop = l$nPop)
  # track density with a radius of 1.5 units
  tdDF <- calculateTrackDensity(dataList = tmList, radius = l$radius)
  # fractal dimension
  fdDF <- calculateFD(dataList = tmList)
  # create the report for this dataset
  p <- makeSummaryReport(tmList = tmList, msdList = msdObj, jumpList = jdObj, tddf = tdDF, fddf = fdDF, titleStr = "Report", subStr = "", auto = FALSE, summary = FALSE)

  return(p)
}
