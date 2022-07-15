#' Make a plot of MSD data
#'
#' Input is the output from CalculateMSD()
#'
#' @param df MSD summary = output from calculateMSD()
#' @param units string to describe time units (default is s, seconds)
#' @param xlog boolean to request log10 x axis
#' @param ylog boolean to request log10 y axis
#' @examples
#' @return S3 ggplot
#' @export


plotMSD <- function(df, units = "s", xlog = FALSE, ylog = FALSE) {
  pred <- NULL
  xlab <- paste0("Time (",units,")")

  df$pred <- 4 * (df$mean[1] / (4 * df$t[1])) * df$t
  p <- ggplot(df, aes(x = t, y = mean)) +
    geom_line() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0) +
    geom_point() +
    geom_line(aes(x = t, y = pred, col = "red")) +
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

