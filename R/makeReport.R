#' Make Report
#'
#' Generate several plots to visualise TrackMate data and generate a report.
#'
#' @param df imported TrackMate data with correct units
#' @param msddf MSD summary = output from calculateMSD()
#' @param titleStr string used as the title for the report
#' @examples
#' xmlPath <- "~/Desktop/FakeTracks.xml"
#' data <- readTrackMateXML(XMLpath = xmlPath)
#' data <- correctTrackMateData(data, xy = 0.04)
#' msdDF <- calculateMSD(data, method = "ensemble", N = 3, short = 8)
#' fileName <- tools::file_path_sans_ext(xmlPath)
#' makeReport(data, msdDF, fileName)
#' @return patchwork ggplot
#' @export

makeReport <- function(df, msddf, titleStr) {
  oldw <- getOption("warn")
  options(warn = -1)

  x <- y <- displacement <- track_duration <- cumulative_distance <- NULL

  # make msd plot
  p_msd <- plotMSD(msddf)

  # plot all tracks colour coded by trace number
  p_allTracks <-   ggplot(data = df, aes(x = x, y = y)) +
    geom_path(aes(colour = trace, alpha = 0.5)) +
    lims(x = c(0,max(df$x,df$y)),y = c(max(df$x,df$y),0)) +
    coord_fixed() +
    theme_bw() +
    theme(legend.position = "none")

  # plot displacement over time
  p_displacementOverTime <- ggplot(data = df, aes(x = t, y = displacement)) +
    geom_path(aes(y = rollmean(displacement, 20, na.pad = TRUE), group = trace, alpha = 0.1)) +
    geom_smooth(aes(x = t, y = displacement)) +
    ylim(0,NA) +
    labs(x = "Time (s)", y = "Distance (\u03bcm)") +
    theme_classic() +
    theme(legend.position = "none")

  # plot cumulative distance over time
  p_cumdistOverTime <- ggplot(data = df, aes(x = track_duration, y = cumulative_distance, group = trace)) +
    geom_path(aes(alpha = 0.1)) +
    labs(x = "Duration (s)", y = "Cumulative Distance (\u03bcm)") +
    theme_classic() +
    theme(legend.position = "none")

  # ggplot histogram of displacements
  p_displacementHist <- ggplot(data = df, aes(x = displacement)) +
    geom_histogram() +
    labs(x = "Displacement (\u03bcm)", y = "Frequency") +
    theme_classic() +
    theme(legend.position = "none")

  r_report <- p_allTracks | (p_displacementOverTime + p_displacementHist) / (p_cumdistOverTime + p_msd)
  r_report <- r_report + plot_annotation(title = titleStr)

  options(warn = oldw)

  return(r_report)
}
