#' Make a plot of MSD data
#'
#' Generate a plot of MSD over a series of increasing time lags.
#' Input is the output from CalculateMSD(), so the plot will display the ensemble or time-averaged MSD (whatever was requested)
#' A fit to the first four points is displayed to evaluate alpha. Diffusion coefficient from this fit is displayed top-left.
#'
#' @param df MSD summary = output from calculateMSD()
#' @param units string to describe time units (default is s, seconds)
#' @param bars boolean to request error bars (1 x SD)
#' @param xlog boolean to request log10 x axis
#' @param ylog boolean to request log10 y axis
#' @examples
#' xmlPath <- "~/Desktop/FakeTracks.xml"
#' data <- readTrackMateXML(XMLpath = xmlPath)
#' data <- correctTrackMateData(data, xy = 0.04)
#' msdDF <- calculateMSD(data, method = "ensemble", N = 3, short = 8)
#' plotMSD(msdDF, bars = FALSE)
#' @return S3 ggplot
#' @export


plotMSD <- function(df, units = "s", bars = TRUE, xlog = FALSE, ylog = FALSE) {
  xlab <- paste0("Time (",units,")")
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
      geom_line() +
      geom_point()
  }

  p <- p + geom_line(aes(x = t, y = pred, col = "red")) +
    geom_text(aes(label = paste0("D = ",format(round(dee,3), nsmall = 3)), x = min(df$t), y = Inf), hjust = 0, vjust = 1) +
    labs(x = xlab, y = "MSD") +
    theme_classic() +
    theme(legend.position = "none")
  if(xlog) {
    p <- p + scale_x_log10()
  }
  if(ylog) {
    p <- p + scale_y_log10()
  }

  return(p)
}

