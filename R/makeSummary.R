#' Make Summary
#'
#' Generate several plots to visualise TrackMate data and generate a report.
#' Note that the units are hard-coded as um and s.
#' The use of um is because ggsave does not currently save unicode to PDF reliably.
#'
#' @param df imported TrackMate data with correct units
#' @param msddf MSD summary = output from calculateMSD()
#' @param units character vector to describe units (defaults are um, micrometres and  s, seconds)
#' @param titleStr string used as the title for the summary
#' @param subStr string used as the subtitle for the summary
#' @return patchwork ggplot
#' @export

makeSummary <- function(df, msddf, units = c("um","s"), titleStr = "Summary", subStr = NULL) {
  oldw <- getOption("warn")
  options(warn = -1)

  x <- y <- displacement <- track_duration <- cumulative_distance <- dataid <- NULL
  xstr <- ystr <- NULL

  ndata <- length(unique(df$dataid))
  alphaLevel <- ifelse(ndata < 4, 0.5, ifelse(ndata < 8, 0.25,0.1))

  # make msd plot
  p_msd <- plotNMSD(msddf)

  # plot all tracks colour coded by trace number
  p_allTracks <- ggplot(data = df, aes(x = x, y = y, group = interaction(dataid, trace))) +
    geom_path(aes(colour = dataid, alpha = alphaLevel)) +
    lims(x = c(0,max(df$x,df$y)),y = c(max(df$x,df$y),0)) +
    labs(x = xstr, y = ystr) +
    coord_fixed() +
    theme_bw() +
    theme(legend.position = "none")

  # plot displacement over time
  xstr <- paste0("Time (",units[2],")")
  ystr <- paste0("Displacement (",units[1],")")
  p_displacementOverTime <- ggplot(data = df, aes(x = t, y = displacement)) +
    geom_path(aes(y = rollmean(displacement, 20, na.pad = TRUE), group = interaction(dataid, trace), alpha = 0.01)) +
    geom_smooth(method = "gam", formula = (y ~ s(x, bs = 'cs'))) +
    ylim(0,NA) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  # plot cumulative distance over time
  xstr <- paste0("Duration (",units[2],")")
  ystr <- paste0("Cumulative Distance (",units[1],")")
  p_cumdistOverTime <- ggplot(data = df, aes(x = track_duration, y = cumulative_distance, group = interaction(dataid, trace))) +
    geom_path(aes(colour = dataid, alpha = alphaLevel)) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  # ggplot histogram of displacements
  # breaks <- pretty(range(df$displacement), n = nclass.FD(df$displacement), min.n = 1)
  # bwidth <- breaks[2]-breaks[1]
  nBin <- floor(1 + log2(nrow(df)))
  xstr <- paste0("Displacement (",units[2],")")
  ystr <- "Frequency"
  p_displacementHist <- ggplot(data = df, aes(x = displacement)) +
    geom_histogram(bins = nBin) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  r_report <- p_allTracks | (p_displacementOverTime + p_displacementHist) / (p_cumdistOverTime + p_msd)
  r_report <- r_report + plot_annotation(title = titleStr, subtitle = subStr)

  options(warn = oldw)

  return(r_report)
}
