#' Make a plot of all tracks
#'
#' @param input either a data frame of TrackMate data or list of TrackMate data frame and calibration data frame
#' @param summary boolean to specify if plot is of one dataset or several related datasets
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @param alphaLevel numeric variable to set alpha for the plot
#' @return ggplot
#' @export
plot_tm_allTracks <- function(input, summary = FALSE, xstr = "", ystr = "", alphaLevel = 0.5) {
  x <- y <- dataid <- NULL
  if(inherits(input, "list")) {
    df <- input[[1]]
  } else {
    df <- input
  }
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
#' @param input either a data frame of TrackMate data or list of TrackMate data frame and calibration data frame
#' @param summary boolean to specify if plot is of one dataset or several related datasets
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @return ggplot
#' @export
plot_tm_displacementOverTime <- function(input, summary = FALSE, xstr = NULL, ystr = NULL) {
  t <- displacement <- dataid <- NULL
  if(inherits(input, "list")) {
    df <- input[[1]]
    calibration <- input[[2]]
    units <- calibration$unit[1:2]
    xstr <- paste0("Time (",units[2],")")
    ystr <- paste0("Displacement (",units[1],")")
  } else {
    df <- input
  }

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
#' @param input either a data frame of TrackMate data or list of TrackMate data frame and calibration data frame
#' @param summary boolean to specify if plot is of one dataset or several related datasets
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @param alphaLevel numeric variable to set alpha for the plot
#' @return ggplot
#' @export
plot_tm_cumdistOverTime <- function(input, summary = FALSE, xstr = NULL, ystr = NULL, alphaLevel = 0.1) {
  track_duration <- cumulative_distance <- dataid <- NULL

  if(inherits(input, "list")) {
    df <- input[[1]]
    calibration <- input[[2]]
    units <- calibration$unit[1:2]
    xstr <- paste0("Time (",units[2],")")
    ystr <- paste0("Displacement (",units[1],")")
  } else {
    df <- input
  }

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
#' @param input either a data frame of TrackMate data or list of TrackMate data frame and calibration data frame
#' @param summary boolean to specify if plot is of one dataset or several related datasets
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @param auto boolean to switch between returning a ggplot and a list of ggplot and variable
#' @return ggplot or list of ggplot and variable
#' @export
plot_tm_displacementHist <- function(input, summary = FALSE, xstr = NULL, ystr = NULL, auto = FALSE) {
  displacement <- NULL

  if(inherits(input, "list")) {
    df <- input[[1]]
    calibration <- input[[2]]
    units <- calibration$unit[1:2]
    xstr <- paste0("Displacement (",units[1],")")
    ystr <- "Frequency"
  } else {
    df <- input
  }

  median_disp <- median(df$displacement)
  nBin <- floor(1 + log2(nrow(df)))

  p <- ggplot(data = df, aes(x = displacement)) +
    geom_histogram(bins = nBin) +
    geom_text(aes(label = paste0("median = ",format(round(median_disp,3), nsmall = 3)), x = max(displacement, na.rm = TRUE), y = Inf), size = 3, hjust = 1, vjust = 1, check_overlap = TRUE) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  if(auto) {
    returnList <- list(p, median_disp)
    return(returnList)
  } else {
    return(p)
  }
}

#' Make a histogram of alpha values
#'
#' @param df data frame of alpha values
#' @param median_alpha variable for adding label to plot
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @return ggplot
#' @export
plot_tm_alpha <- function(df, median_alpha = NULL, xstr = "alpha (log2)", ystr = "Frequency") {

  p <- ggplot(data = df, aes(x = alpha)) +
    geom_histogram(binwidth = 0.1)

  if(!missing(median_alpha)) {
    p <- p + geom_text(aes(label = paste0("median = ",format(round(median_alpha,3), nsmall = 3)), x = min(alpha, na.rm = TRUE), y = Inf), size = 3, hjust = 0, vjust = 1, check_overlap = TRUE)
  }
  p <- p + labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  return(p)
}


#' Make a histogram of average speed
#'
#' @param input either a data frame of TrackMate data or list of TrackMate data frame and calibration data frame
#' @param summary boolean to specify if plot is of one dataset or several related datasets
#' @param xstr string to label x-axis
#' @param ystr string to label y-axis
#' @param auto boolean to switch between returning a ggplot and a list of ggplot and variable
#' @return ggplot or list of ggplot and variable
#' @export
plot_tm_speed <- function(input, summary = FALSE, xstr = NULL, ystr = NULL, auto = FALSE) {
  dataid <- cumulative_distance <- track_duration <- speed <- NULL

  if(inherits(input, "list")) {
    df <- input[[1]]
    calibration <- input[[2]]
    units <- calibration$unit[1:2]
    xstr <- paste0("Speed (",units[1],"/",units[2],")")
    ystr <- "Frequency"
  } else {
    df <- input
  }

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
  median_speed <- median(speedDF$speed, na.rm = TRUE)
  nBin <- max(floor(1 + log2(nrow(speedDF))),30)

  p <- ggplot(data = speedDF, aes(x = speed)) +
    geom_histogram(bins = nBin) +
    geom_text(aes(label = paste0("median = ",format(round(median_speed,3), nsmall = 3)), x = max(speed, na.rm = TRUE), y = Inf), size = 3, hjust = 1, vjust = 1, check_overlap = TRUE) +
    labs(x = xstr, y = ystr) +
    theme_classic() +
    theme(legend.position = "none")

  if(auto) {
    returnList <- list(p, median_speed)
    return(returnList)
  } else {
    return(p)
  }
}


#' Make a histogram of track density (number of neighbours)
#'
#' @param df data frame of TrackMate data
#' @param auto boolean to switch between returning a ggplot and a list of ggplot and variable
#' @return ggplot or list of ggplot and variable
#' @export
plot_tm_neighbours <- function(df, auto = FALSE) {
  density <- NULL
  median_density <- median(df$density, na.rm = TRUE)
  nBin <- max(floor(1 + log2(nrow(df))),30)

  p <- ggplot(data = df, aes(x = density)) +
    geom_histogram(bins = nBin) +
    geom_text(aes(label = paste0("median = ",format(round(median_density,3), nsmall = 3)), x = max(density, na.rm = TRUE), y = Inf), size = 3, hjust = 1, vjust = 1, check_overlap = TRUE) +
    labs(x = "Density", y = "Frequency") +
    theme_classic() +
    theme(legend.position = "none")

  if(auto) {
    returnList <- list(p, median_density)
    return(returnList)
  } else {
    return(p)
  }
}
