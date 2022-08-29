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
    ystr <- paste0("Cumulative distance (",units[1],")")
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
    labs(x = "Track density", y = "Frequency") +
    theme_classic() +
    theme(legend.position = "none")

  if(auto) {
    returnList <- list(p, median_density)
    return(returnList)
  } else {
    return(p)
  }
}


#' Make a plot of MSD data
#'
#' Generate a plot of MSD over a series of increasing time lags.
#' Input is the output from CalculateMSD(), so the plot will display the ensemble or time-averaged MSD (whatever was requested)
#' A fit to the first four points is displayed to evaluate alpha. Diffusion coefficient from this fit is displayed top-left.
#'
#' @param df MSD summary = output from calculateMSD()
#' @param units character vector to describe units (defaults are um, micrometres and  s, seconds)
#' @param bars boolean to request error bars (1 x SD)
#' @param xlog boolean to request log10 x axis
#' @param ylog boolean to request log10 y axis
#' @param auto boolean to request plot only, TRUE gives plot and D as a list
#' @examples
#' xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
#' datalist <- readTrackMateXML(XMLpath = xmlPath)
#' data <-  datalist[[1]]
#' # use the ensemble method and only look at tracks with more than 8 points
#' msdobj <- calculateMSD(df = data, method = "ensemble", N = 3, short = 8)
#' msddf <- msdobj[[1]]
#' plot_tm_MSD(msddf, bars = FALSE)
#' @return ggplot or ggplot and variable
#' @export
plot_tm_MSD <- function(df, units = c("um","s"), bars = FALSE, xlog = FALSE, ylog = FALSE, auto = FALSE) {
  xlab <- paste0("Time (",units[2],")")
  pred <- NULL

  # fit to first four data points
  mod <- lm(mean ~ t, weights = c(n), data = df[1:4,])
  # make a column containing model y points for each t
  df$pred <- (mod$coefficients[2] * df$t) + mod$coefficients[1]
  # calculate diffusion constant (D)
  dee <- df$pred[1] / (4 * df$t[1])
  # make ggplot
  if(bars) {
    p <- ggplot(df, aes(x = t, y = mean)) +
      geom_line() +
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0) +
      geom_point()
  } else {
    p <- ggplot(df, aes(x = t, y = mean)) +
      geom_line(linesize = 1)
  }

  p <- p + geom_line(aes(x = t, y = pred), colour = "red", linetype = 2) +
    geom_text(aes(label = paste0("D = ",format(round(dee,3), nsmall = 3)), x = min(t), y = Inf), size = 3, hjust = 0, vjust = 1, check_overlap = TRUE) +
    labs(x = xlab, y = "MSD") +
    theme_classic() +
    theme(legend.position = "none")
  if(xlog) {
    p <- p + scale_x_log10()
  }
  if(ylog) {
    p <- p + scale_y_log10()
  }

  listReturn <- list(p, dee)

  if(auto) {
    return(listReturn)
  } else {
    return(p)
  }
}



#' Plot several (n) MSD curves
#'
#' Generate a plot of several MSD curves together with a summary curve.
#' This function is used to compile multiple datasets from the same condition.
#'
#' @param df dataframe of MSD summary data from multiple datasets (labelled by dataid)
#' @param auto boolean to request plot only, TRUE gives plot and summary dataframe as a list
#' @return ggplot or list of ggplot and dataframe of summary data
#' @export
plot_tm_NMSD <- function(df, auto = FALSE) {
  dataid <- pred <- value <- size <- NULL
  # generate a mean of the MSD curve over time (lag)
  # rename mean to value so as not to upset ddplyr
  colnames(df)[1] <- "value"

  # find the parameters for interpolation
  t1 <- df %>%
    subset(size == 1)
  minT <- min(t1$t)
  maxT <- max(df$t)
  steps <- ceiling(maxT / minT)

  # need to work per dataset
  datasets <- unique(df$dataid)

  for (i in datasets) {
    temp <- df %>%
      subset(dataid == i)
    newdf <- data.frame(approx(x = temp$t, y = temp$value, xout = seq(from = minT, to = steps * minT, by = minT)))
    newdf$dataid <- i
    if(i == datasets[1]) {
      alldf <- newdf
    } else {
      alldf <- rbind(alldf,newdf)
    }
  }
  # rename header
  names(alldf) <- c("t", "value", "dataid")

  msdmean <- alldf %>%
    group_by(t) %>%
    summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), n = sum(!is.na(value)))

  minN <- ceiling(max(msdmean$n) / 3)
  tmp <- msdmean[-1,]
  tmp <- tmp %>%
    filter(n < minN)
  if(nrow(tmp) > 0) {
    maxX <- min(tmp$t)
  } else {
    maxX <- NA
  }

  p <- ggplot(data = df, aes(x = t, y = value)) +
    geom_line(aes(group = dataid), colour = "blue", alpha = 0.5) +
    geom_ribbon(data = msdmean, aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2) +
    geom_line(data = msdmean, aes(x = t, y = mean), size = 1) +
    ylim(0,NA) +
    xlim(0,maxX) +
    labs(x = "Time (s)", y = "MSD") +
    theme_classic() +
    theme(legend.position = "none")

  # fit to first four data points
  mod <- lm(mean ~ t, data = msdmean[1:4,])
  # make a column containing model y points for each t
  msdmean$pred <- (mod$coefficients[2] * msdmean$t) + mod$coefficients[1]
  # calculate diffusion constant (D)
  dee <- msdmean$pred[1] / (4 * msdmean$t[1])
  # add line to show MSD with alpha = 1
  p <-  p + geom_line(data = msdmean, aes(x = t, y = pred), colour = "red", linetype = 2) +
    geom_text(aes(label = paste0("D = ",format(round(dee,3), nsmall = 3)), x = min(msdmean, na.rm = TRUE), y = Inf), size = 3, hjust = 0, vjust = 1, check_overlap = TRUE)

  # for export we only want the summary data that is plotted i.e. cut off appropriately
  if(!is.na(maxX)) {
    cutrow <- which(msdmean$t == maxX)
    if(!is.null(cutrow) & cutrow > 3) {
      msdmean <- msdmean[1:cutrow - 1,]
    }
  }

  if(auto) {
    listReturn <- list(p, msdmean)
    return(listReturn)
  } else {
    return(p)
  }
}

