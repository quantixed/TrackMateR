#' Make Summary Report
#'
#' Generate several plots to visualise TrackMate data and generate a report.
#' The use of um is because ggsave does not currently save unicode to PDF reliably.
#' The code is switchable to accommodate making a "report" (one dataset) or a "summary" (several related datasets combined)
#'
#' @param tmList list of trackmate data and calibration
#' @param msdList MSD summary and alpha list = output from calculateMSD()
#' @param jumpList list of a data frame of jump data and a variable to be passed to timeRes
#' @param tddf data frame of track density data
#' @param titleStr string used as the title for the report
#' @param subStr string used as the subtitle for the report
#' @param auto boolean which selects for returning the patchwork report or a list of the patchwork report and a data frame of summary
#' @param summary boolean which selects for summary
#' @examples
#' xmlPath <- "~/Desktop/FakeTracks.xml"
#' tmObj <- readTrackMateXML(XMLpath = xmlPath)
#' tmObj <- correctTrackMateData(tmObj, xyscalar = 0.04)
#' tmDF <- tmObj[[1]]
#' calibrationDF <- tmObj[[2]]
#' msdObj <- calculateMSD(df = tmDF, method = "ensemble", N = 3, short = 8)
#' jdObj <- calculateJD(dataList = tmObj, deltaT = 1)
#' tdDF <- calculateTrackDensity(dataList = tmObj, radius = 1.5)
#' fileName <- tools::file_path_sans_ext(basename(xmlPath))
#' reportObj <- makeSummaryReport(tmList = tmObj, msdList = msdObj, jumpList = jdObj, tddf = tdDF,
#' titleStr = "Report", subStr = fileName, auto = TRUE)
#' @return list of patchwork ggplot and data frame of summary data
#' @export

makeSummaryReport <- function(tmList, msdList, jumpList, tddf, titleStr = "", subStr = "", auto = FALSE, summary = FALSE) {
  oldw <- getOption("warn")
  options(warn = -1)

  x <- y <- displacement <- track_duration <- cumulative_distance <- speed <- dataid <- density <- NULL
  xstr <- ystr <- alphaValue <- NULL

  # get ready for plotting
  df <- tmList[[1]]
  calibration <- tmList[[2]]
  units <- calibration$unit[1:2]
  if(summary) {
    ndata <- length(unique(df$dataid))
    alphaLevel <- ifelse(ndata < 4, 0.5, ifelse(ndata < 8, 0.25,0.1))
  }

  # make msd plot
  msddf <- msdList[[1]]
  if(summary){
    p_msd <- plotNMSD(msddf)
  } else {
    msdreturn <- plotMSD(msddf, units, auto = TRUE)
    p_msd <- msdreturn[[1]]
    dee <- msdreturn[[2]]
  }

  # plot all tracks colour coded by trace number (or by dataset)
  p_allTracks <- plot_tm_allTracks(df = df, summary = summary, xstr = xstr, ystr = ystr, alphaLevel = alphaLevel)

  # plot displacement over time
  xstr <- paste0("Time (",units[2],")")
  ystr <- paste0("Displacement (",units[1],")")
  p_displacementOverTime <- plot_tm_displacementOverTime(df = df, summary = summary, xstr = xstr, ystr = ystr)

  # plot cumulative distance over time
  xstr <- paste0("Duration (",units[2],")")
  ystr <- paste0("Cumulative Distance (",units[1],")")
  p_cumdistOverTime <- plot_tm_cumdistOverTime(df = df, summary = summary, xstr = xstr, ystr = ystr, alphaLevel = alphaLevel)

  # ggplot histogram of displacements
  xstr <- paste0("Displacement (",units[1],")")
  ystr <- "Frequency"
  p_displacementHist <- plot_tm_displacementHist(df = df, xstr = xstr, ystr = ystr)

  # alpha distribution of traces
  alphas <- msdList[[2]]
  alphas <- na.omit(alphas)
  # within a sensible range (log 2)
  identify <- alphas$alpha <= 4 & alphas$alpha >= -4
  # we will take log2
  # alphas <- data.frame(alpha = alphas[identify,])
  alphas <- subset(alphas, identify)
  # convert back to real numbers
  median_alpha <- 2^(median(alphas$alpha, na.rm = TRUE))

  xstr <- "alpha (log2)"
  ystr <- "Frequency"
  p_alpha <- plot_tm_alpha(df = alphas, median_alpha = median_alpha, xstr = xstr, ystr = ystr)

  # make a plot of average speed per track
  xstr <- paste0("Speed (",units[1],"/",units[2],")")
  ystr <- "Frequency"
  if(auto == TRUE & summary == FALSE) {
    speedObj <- plot_tm_speed(df = df, summary = summary, xstr = xstr, ystr = ystr, auto = auto)
    p_speed <- speedObj[[1]]
    median_speed <- speedObj[[2]]
  } else {
    p_speed <- plot_tm_speed(df = df, summary = summary, xstr = xstr, ystr = ystr, auto = auto)
  }

  # make a plot of jump distance distribution
  jdDF <- jumpList[[1]]
  jumptime <- jumpList[[2]]
  p_jump <- fittingJD(df = jdDF, mode = "ECDF", nPop = 2, units = units, breaks = 100, timeRes = jumptime)

  # calculate neighbours within 1.5 units
  p_neighbours <- plot_tm_neighbours(df = tddf)

  # make the report (patchwork of ggplots)
  # r_report <- (p_allTracks | ((p_displacementOverTime + p_displacementHist) / (p_cumdistOverTime + p_speed))) / (p_msd + p_alpha + p_jump)
  design <- "
    1123
    1145
  "
  toprow <- p_allTracks + p_displacementOverTime + p_displacementHist + p_cumdistOverTime + p_speed + plot_layout(design = design)
  bottomrow <- p_msd + p_alpha + p_jump + p_neighbours + plot_layout(ncol = 4)
  r_report <- toprow / bottomrow + plot_layout(nrow = 2, heights = c(2,1))
  r_report <- r_report + plot_annotation(title = titleStr, subtitle = subStr)

  if(auto) {
    # make a multi-column object and pass directly back
    df_report <- data.frame(alpha = median_alpha,
                            speed = median_speed,
                            dee = dee)
    returnList <- list(r_report,df_report)
    return(returnList)
  } else {
    return(r_report)
  }

  options(warn = oldw)

}
