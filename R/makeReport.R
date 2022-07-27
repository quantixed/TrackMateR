#' Make Report
#'
#' Generate several plots to visualise TrackMate data and generate a report.
#' Note that the units are hard-coded as um and s.
#' The use of um is because ggsave does not currently save unicode to PDF reliably.
#'
#' @param df imported TrackMate data with correct units
#' @param msdlist MSD summary and alpha list = output from calculateMSD()
#' @param units character vector to describe units (defaults are um, micrometres and  s, seconds)
#' @param titleStr string used as the title for the report
#' @param subStr string used as the subtitle for the report
#' @examples
#' xmlPath <- "~/Desktop/FakeTracks.xml"
#' data <- readTrackMateXML(XMLpath = xmlPath)
#' datalist <- readTrackMateXML(XMLpath = xmlPath)
#' data <-  datalist[[1]]
#' msdobj <- calculateMSD(df = data, method = "ensemble", N = 3, short = 8)
#' fileName <- tools::file_path_sans_ext(basename(xmlPath))
#' makeReport(df = data, msdlist = msdobj, titleStr = "Report", subStr = fileName)
#' @return patchwork ggplot
#' @export

makeReport <- function(df, msdlist, units = c("um","s"), titleStr = "", subStr = "") {
  oldw <- getOption("warn")
  options(warn = -1)

  x <- y <- displacement <- track_duration <- cumulative_distance <- speed <- NULL
  xstr <- ystr <- NULL

  msddf <- msdlist[[1]]

  # make msd plot
  p_msd <- plotMSD(msddf, units)

  # plot all tracks colour coded by trace number
  p_allTracks <- ggplot(data = df, aes(x = x, y = y)) +
    geom_path(aes(colour = trace), alpha = 0.5) +
    lims(x = c(0,max(df$x,df$y)),y = c(max(df$x,df$y),0)) +
    labs(x = xstr, y = ystr) +
    coord_fixed() +
    theme_bw() +
    theme(legend.position = "none")

  # plot displacement over time
  xstr <- paste0("Time (",units[2],")")
  ystr <- paste0("Displacement (",units[1],")")
  p_displacementOverTime <- ggplot(data = df, aes(x = t, y = displacement)) +
    geom_path(aes(y = rollmean(displacement, 20, na.pad = TRUE), group = trace, alpha = 0.1)) +
    geom_smooth(method = "gam", formula = (y ~ s(x, bs = 'cs'))) +
    ylim(0,NA) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  # plot cumulative distance over time
  xstr <- paste0("Duration (",units[2],")")
  ystr <- paste0("Cumulative Distance (",units[1],")")
  p_cumdistOverTime <- ggplot(data = df, aes(x = track_duration, y = cumulative_distance, group = trace)) +
    geom_path(aes(alpha = 0.1)) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  # ggplot histogram of displacements
  xstr <- paste0("Displacement (",units[2],")")
  ystr <- "Frequency"
  nBin <- floor(1 + log2(nrow(df)))
  p_displacementHist <- ggplot(data = df, aes(x = displacement)) +
    geom_histogram(bins = nBin) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  # alpha distribution of traces
  alphas <- data.frame(alpha = msdlist[[2]])
  alphas <- na.omit(alphas)
  # within a sensible range (log 2)
  identify <- alphas$alpha <= 4 & alphas$alpha >= -4
  # we will take log2
  alphas <- data.frame(alpha = alphas[identify,])
  # convert back to real numbers
  median_alpha <- 2^(median(alphas$alpha, na.rm = TRUE))

  xstr <- "alpha (log2)"
  ystr <- "Frequency"
  p_alpha <- ggplot(data = alphas, aes(x = alpha)) +
    geom_histogram(binwidth = 0.1) +
    geom_text(aes(label = paste0("median = ",format(round(median_alpha,3), nsmall = 3)), x = min(alpha, na.rm = TRUE), y = Inf), hjust = 0, vjust = 1) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  # make a plot of average speed per track
  speedDF <- df %>%
    group_by(trace) %>%
    summarise(cumdist = max(cumulative_distance), cumtime = max(track_duration))
  speedDF$speed <- speedDF$cumdist / speedDF$cumtime
  nBin <- max(floor(1 + log2(nrow(speedDF))),30)

  xstr <- paste0("Speed (",units[1],"/",units[2],")")
  ystr <- "Frequency"
  p_speed <- ggplot(data = speedDF, aes(x = speed)) +
    geom_histogram(bins = nBin) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  r_report <- (p_allTracks | ((p_displacementOverTime + p_displacementHist) / (p_cumdistOverTime + p_speed))) / (p_msd + p_alpha)
  r_report <- r_report + plot_annotation(title = titleStr, subtitle = subStr)

  options(warn = oldw)

  return(r_report)
}
