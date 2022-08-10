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
#' plotMSD(msddf, bars = FALSE)
#' @return ggplot or ggplot and variable
#' @export


plotMSD <- function(df, units = c("um","s"), bars = FALSE, xlog = FALSE, ylog = FALSE, auto = FALSE) {
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

