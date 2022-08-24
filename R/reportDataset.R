#' Report Dataset
#'
#' This function automates the process of making a report.
#' If more control is needed over the default parameters:
#' a user can generate MSD, jump distance and track density objects from their imported TrackMate data, and make a report manually.
#'
#' @param tmList list of TrackMate data and calibration
#' @return patchwork ggplot
#' @export
reportDataset <- function(tmList) {

  tmDF <- tmList[[1]]
  calibrationDF <- tmList[[2]]
  # take the units
  units <- calibrationDF$unit[1:2]
  # calculate MSD
  msdObj <- calculateMSD(tmDF, N = 3, short = 8)
  msdDF <- msdObj[[1]]
  alphaDF <- msdObj[[2]]
  # jump distance calc with deltaT of 1
  deltaT <- 1
  jdObj <- calculateJD(dataList = tmList, deltaT = deltaT)
  # track density with a radius of 1.5 units
  tdDF <- calculateTrackDensity(dataList = tmList, radius = 1.5)

  # create the report for this dataset
  p <- makeSummaryReport(tmList = tmList, msdList = msdObj, jumpList = jdObj, tddf = tdDF, titleStr = "Report", subStr = "", auto = FALSE, summary = FALSE)

  return(p)
}
