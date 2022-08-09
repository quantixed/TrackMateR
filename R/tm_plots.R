#' Make a plot of all tracks
#'
#' @param df data frame of TrackMate data
#' @param summary boolean to specify if plot is of one dataset or several related datasets
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @param alphaLevel numeric variable to set alpha for the plot
#' @return ggplot
#' @export
plot_tm_allTracks <- function(df, summary = FALSE, xstr = "x", ystr = "y", alphaLevel = 0.5) {
  x <- y <- dataid <- NULL
  if(!summary) {
    alphaLevel <- 0.5
  }
  if(summary) {
    p <- ggplot(data = df, aes(x = x, y = y, group = interaction(dataid, trace))) +
      geom_path(aes(colour = dataid), alpha = alphaLevel)
  } else {
    p <- ggplot(data = df, aes(x = x, y = y)) +
      geom_path(aes(colour = trace), alpha = alphaLevel)
  }
  p <- p + lims(x = c(0,max(df$x,df$y)),y = c(max(df$x,df$y),0)) +
    labs(x = xstr, y = ystr) +
    coord_fixed() +
    theme_bw() +
    theme(legend.position = "none")

  return(p)
}


#' Make a plot of displacement over time
#'
#' @param df data frame of TrackMate data
#' @param summary boolean to specify if plot is of one dataset or several related datasets
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @return ggplot
#' @export
plot_tm_displacementOverTime <- function(df, summary = FALSE, xstr = NULL, ystr = NULL) {
  t <- displacement <- dataid <- NULL

  p <- ggplot(data = df, aes(x = t, y = displacement))
  if(summary) {
    p <- p + geom_path(aes(y = rollmean(displacement, 20, na.pad = TRUE), group = interaction(dataid, trace), alpha = 0.01))
  } else {
    p <- p + geom_path(aes(y = rollmean(displacement, 20, na.pad = TRUE), group = trace, alpha = 0.1))
  }

  p <- p + geom_smooth(method = "gam", formula = (y ~ s(x, bs = 'cs'))) +
    ylim(0,NA) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  return(p)
}

#' Make a plot of cumulative distance over time
#'
#' @param df data frame of TrackMate data
#' @param summary boolean to specify if plot is of one dataset or several related datasets
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @param alphaLevel numeric variable to set alpha for the plot
#' @return ggplot
#' @export
plot_tm_cumdistOverTime <- function(df, summary = FALSE, xstr = NULL, ystr = NULL, alphaLevel = 0.1) {
  track_duration <- cumulative_distance <- dataid <- NULL
  if(!summary) {
    alphaLevel <- 0.1
  }
  if(summary) {
    p <- ggplot(data = df, aes(x = track_duration, y = cumulative_distance, group = interaction(dataid, trace))) +
      geom_path(aes(colour = dataid), alpha = alphaLevel)
  } else {
    p <- ggplot(data = df, aes(x = track_duration, y = cumulative_distance, group = trace)) +
      geom_path(aes(alpha = alphaLevel))
  }
  p <- p + labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  return(p)
}

#' Make a histogram of displacements
#'
#' @param df data frame of TrackMate data
#' @param summary boolean to specify if plot is of one dataset or several related datasets
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @return ggplot
#' @export
plot_tm_displacementHist <- function(df, summary = FALSE, xstr = NULL, ystr = NULL) {
  displacement <- NULL
  nBin <- floor(1 + log2(nrow(df)))

  p <- ggplot(data = df, aes(x = displacement)) +
    geom_histogram(bins = nBin) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  return(p)
}


#' Make a histogram of alpha values
#'
#' @param df data frame of TrackMate data
#' @param median_alpha variable for adding label to plot
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @return ggplot
#' @export
plot_tm_alpha <- function(df, median_alpha = NULL, xstr = NULL, ystr = NULL) {
  p <- ggplot(data = df, aes(x = alpha)) +
    geom_histogram(binwidth = 0.1) +
    geom_text(aes(label = paste0("median = ",format(round(median_alpha,3), nsmall = 3)), x = min(alpha, na.rm = TRUE), y = Inf), size = 3, hjust = 0, vjust = 1, check_overlap = TRUE) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  return(p)
}


#' Make a histogram of average speed
#'
#' @param df data frame of TrackMate data
#' @param summary boolean to specify if plot is of one dataset or several related datasets
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @param auto boolean to switch between returning a ggplot and a list of ggplot and variable
#' @return ggplot or list of ggplot and variable
#' @export
plot_tm_speed <- function(df, summary = FALSE, xstr = NULL, ystr = NULL, auto = FALSE) {
  dataid <- cumulative_distance <- track_duration <- speed <- NULL
  if(summary) {
    speedDF <- df %>%
      group_by(dataid, trace) %>%
      summarise(cumdist = max(cumulative_distance), cumtime = max(track_duration))
  } else {
    speedDF <- df %>%
      group_by(trace) %>%
      summarise(cumdist = max(cumulative_distance), cumtime = max(track_duration))
  }

  speedDF$speed <- speedDF$cumdist / speedDF$cumtime
  nBin <- max(floor(1 + log2(nrow(speedDF))),30)

  p <- ggplot(data = speedDF, aes(x = speed)) +
    geom_histogram(bins = nBin) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  if(auto) {
    median_speed = median(speedDF$speed, na.rm = TRUE)
    returnList <- list(p,median_speed)
    return(returnList)
  } else {
    return(p)
  }
}


#' Make a histogram of track density (number of neighbours)
#'
#' @param df data frame of TrackMate data
#' @return ggplot
#' @export
plot_tm_neighbours <- function(df) {
  density <- NULL
  nBin <- max(floor(1 + log2(nrow(df))),30)

  p <- ggplot(data = df, aes(x = density)) +
    geom_histogram(bins = nBin) +
    labs(x = "Density", y = "Frequency") +
    theme_classic() +
    theme(legend.position = "none")

  return(p)
}
