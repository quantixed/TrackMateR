#' Make Report
#'
#' Generate several plots to visualise TrackMate data and generate a report.
#' Note that the units are hard-coded as um and s.
#' The use of um is because ggsave does not currently save unicode to PDF reliably.
#'
#' @param df imported TrackMate data with correct units
#' @param msdlist MSD summary and alpha list = output from calculateMSD()
#' @param titleStr string used as the title for the report
#' @param subStr string used as the subtitle for the report
#' @examples
#' xmlPath <- "~/Desktop/FakeTracks.xml"
#' data <- readTrackMateXML(XMLpath = xmlPath)
#' data <- correctTrackMateData(data, xy = 0.04)
#' result <- calculateMSD(data, method = "ensemble", N = 3, short = 8)
#' fileName <- tools::file_path_sans_ext(basename(xmlPath))
#' makeReport(data, result, "", fileName)
#' @return patchwork ggplot
#' @export

makeReport <- function(df, msdlist, titleStr, subStr) {
  oldw <- getOption("warn")
  options(warn = -1)

  x <- y <- displacement <- track_duration <- cumulative_distance <- speed <- NULL

  msddf <- msdlist[[1]]

  # make msd plot
  p_msd <- plotMSD(msddf)

  # plot all tracks colour coded by trace number
  p_allTracks <- ggplot(data = df, aes(x = x, y = y)) +
    geom_path(aes(colour = trace, alpha = 0.5)) +
    lims(x = c(0,max(df$x,df$y)),y = c(max(df$x,df$y),0)) +
    coord_fixed() +
    theme_bw() +
    theme(legend.position = "none")

  # plot displacement over time
  p_displacementOverTime <- ggplot(data = df, aes(x = t, y = displacement)) +
    geom_path(aes(y = rollmean(displacement, 20, na.pad = TRUE), group = trace, alpha = 0.1)) +
    geom_smooth(method = "gam", formula = (y ~ s(x, bs = 'cs'))) +
    ylim(0,NA) +
    labs(x = "Time (s)", y = "Displacement (um)") +
    theme_classic() +
    theme(legend.position = "none")

  # plot cumulative distance over time
  p_cumdistOverTime <- ggplot(data = df, aes(x = track_duration, y = cumulative_distance, group = trace)) +
    geom_path(aes(alpha = 0.1)) +
    labs(x = "Duration (s)", y = "Cumulative Distance (um)") +
    theme_classic() +
    theme(legend.position = "none")

  # ggplot histogram of displacements
  nBin <- floor(1 + log2(nrow(df)))
  p_displacementHist <- ggplot(data = df, aes(x = displacement)) +
    geom_histogram(binwidth = nBin) +
    labs(x = "Displacement (um)", y = "Frequency") +
    theme_classic() +
    theme(legend.position = "none")

  alphas <- data.frame(alpha = msdlist[[2]])
  alphas <- na.omit(alphas)
  # within a sensible range
  identify <- alphas$alpha < 20 & alphas$alpha > 1/20
  # we will take log2
  alphas <- data.frame(alpha = log2(alphas[identify,]))
  median_alpha <- 2^(median(alphas$alpha, na.rm = TRUE))

  p_alpha <- ggplot(data = alphas, aes(x = alpha)) +
    geom_histogram(binwidth = 0.1) +
    geom_text(aes(label = paste0("median = ",format(round(median_alpha,3), nsmall = 3)), x = -Inf, y = Inf), hjust = 0, vjust = 1) +
    labs(x = "alpha (log2)", y = "Frequency") +
    theme_classic() +
    theme(legend.position = "none")

  # make a plot of average speed per track
  speedDF <- df %>%
    group_by(trace) %>%
    summarise(cumdist = max(cumulative_distance), cumtime = max(track_duration))
  speedDF$speed <- speedDF$cumdist / speedDF$cumtime
  nBin <- max(floor(1 + log2(nrow(speedDF))),30)
  p_speed <- ggplot(data = speedDF, aes(x = speed)) +
    geom_histogram(bins = nBin) +
    labs(x = "Average speed (um/s)", y = "Frequency") +
    theme_classic() +
    theme(legend.position = "none")

  r_report <- (p_allTracks | (p_displacementOverTime + p_displacementHist)) / (p_cumdistOverTime + p_speed + p_msd + p_alpha)
  r_report <- r_report + plot_annotation(title = titleStr, subtitle = subStr)

  options(warn = oldw)

  return(r_report)
}
