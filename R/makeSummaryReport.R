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
#' @param fddf data frame of fractal dimension data
#' @param titleStr string used as the title for the report
#' @param subStr string used as the subtitle for the report
#' @param auto boolean which selects for returning the patchwork report (FALSE) or a list of the patchwork report and a data frame of summary (TRUE)
#' @param summary boolean which selects for report (FALSE) or summary (TRUE)
#' @param ... additional arguments passed to methods
#' @examples
#' xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
#' tmObj <- readTrackMateXML(XMLpath = xmlPath)
#' tmObj <- correctTrackMateData(tmObj, xyscalar = 0.04)
#' tmDF <- tmObj[[1]]
#' calibrationDF <- tmObj[[2]]
#' msdObj <- calculateMSD(df = tmDF, method = "ensemble", N = 3, short = 8)
#' jdObj <- calculateJD(dataList = tmObj, deltaT = 1)
#' tdDF <- calculateTrackDensity(dataList = tmObj, radius = 1.5)
#' fdDF <- calculateFD(dataList = tmObj)
#' fileName <- tools::file_path_sans_ext(basename(xmlPath))
#' reportObj <- makeSummaryReport(tmList = tmObj, msdList = msdObj, jumpList = jdObj,
#' tddf = tdDF, fddf = fdDF, titleStr = "Report", subStr = fileName, auto = TRUE)
#' @return patchwork ggplot or a list of patchwork ggplot and data frame of summary data
#' @export

makeSummaryReport <- function(tmList, msdList, jumpList, tddf, fddf, titleStr = "", subStr = "", auto = FALSE, summary = FALSE, ...) {
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
    msdreturn <- plot_tm_NMSD(msddf, auto = TRUE)
    p_msd <- msdreturn[[1]]
    msdSummary <- msdreturn[[2]]
  } else {
    msdreturn <- plot_tm_MSD(msddf, units, auto = TRUE)
    p_msd <- msdreturn[[1]]
    dee <- msdreturn[[2]]
  }

  # plot all tracks colour coded by trace number (or by dataset)
  p_allTracks <- plot_tm_allTracks(input = df, summary = summary, alphaLevel = alphaLevel)

  # plot displacement over time
  p_displacementOverTime <- plot_tm_displacementOverTime(input = tmList, summary = summary)

  # plot cumulative distance over time
  p_cumdistOverTime <- plot_tm_cumdistOverTime(input = tmList, summary = summary, alphaLevel = alphaLevel)

  # ggplot histogram of displacements
  if(auto == TRUE & summary == FALSE) {
    dispObj <- plot_tm_displacementHist(input = tmList, auto = auto)
    p_displacementHist <- dispObj[[1]]
    median_disp <- dispObj[[2]]
  } else {
    p_displacementHist <- plot_tm_displacementHist(input = tmList, auto = FALSE)
  }

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

  p_alpha <- plot_tm_alpha(df = alphas, median_alpha = median_alpha)

  # make a plot of average speed per track
  if(auto == TRUE & summary == FALSE) {
    speedObj <- plot_tm_speed(input = tmList, summary = summary, auto = auto)
    p_speed <- speedObj[[1]]
    median_speed <- speedObj[[2]]
  } else {
    p_speed <- plot_tm_speed(input = tmList, summary = summary, auto = FALSE)
  }

  # make a plot of jump distance distribution
  jdDF <- jumpList[[1]]
  jumptime <- jumpList[[2]]
  p_jump <- fittingJD(df = jdDF, units = units, timeRes = jumptime, ...)

  # calculate neighbours within 1.5 units
  if(auto == TRUE & summary == FALSE) {
    neighbourObj <- plot_tm_neighbours(df = tddf, auto = auto)
    p_neighbours <- neighbourObj[[1]]
    median_density <- neighbourObj[[2]]
  } else {
    p_neighbours <- plot_tm_neighbours(df = tddf, auto = FALSE)
  }

  # fractal dimension and width plots
  if(auto == TRUE & summary == FALSE) {
    fdObj <- plot_tm_fd(df = fddf, auto = auto)
    p_fd <- fdObj[[1]]
    median_fd <- fdObj[[2]]
    widthObj <- plot_tm_width(df = fddf, units = units, auto = auto)
    p_width <- widthObj[[1]]
    median_width <- widthObj[[2]]
  } else {
    p_fd <- plot_tm_fd(df = fddf, auto = FALSE)
    p_width <- plot_tm_width(df = fddf, units = units, auto = FALSE)
  }

  # make the report (patchwork of ggplots)
  design <- "
    1123
    1145
  "
  toprow <- p_allTracks + p_displacementOverTime + p_displacementHist + p_cumdistOverTime + p_speed + plot_layout(design = design)
  bottomrow <- p_msd + p_alpha + p_jump + p_neighbours + p_fd + p_width + plot_layout(ncol = 3, nrow = 2)
  r_report <- toprow / bottomrow
  r_report <- r_report + plot_annotation(title = titleStr, subtitle = subStr)

  if(auto == TRUE & summary == FALSE) {
    # make a multi-column object and pass directly back
    df_report <- data.frame(alpha = median_alpha,
                            speed = median_speed,
                            dee = dee,
                            displacement = median_disp,
                            neighbours = median_density,
                            fd = median_fd,
                            width = median_width)
    returnList <- list(r_report,df_report)
    return(returnList)
  } else if(auto == TRUE & summary == TRUE) {
    returnList <- list(r_report,msdSummary)
    return(returnList)
  } else {
    return(r_report)
  }

  options(warn = oldw)

}
