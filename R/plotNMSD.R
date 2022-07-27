#' Plot several (n) MSD curves
#'
#' Generate a plot of several MSD curves together with a summary curve.
#'
#' @param df dataframe of MSD summary data from multiple datasets (labelled by dataid)
#' @return ggplot
#' @export

plotNMSD <- function(df) {
  dataid <- pred <- value <- NULL
  # generate a mean of the MSD curve over time (lag)
  # rename mean to value so as not to upset ddplyr
  colnames(df)[1] <- "value"
  msdmean <- df %>%
    group_by(t) %>%
    summarise(mean = mean(value), sd = sd(value))

  p <- ggplot(data = df, aes(x = t, y = value)) +
    geom_line(aes(group = dataid), colour = "blue", alpha = 0.5) +
    geom_ribbon(data = msdmean, aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2) +
    geom_line(data = msdmean, aes(x = t, y = mean), size = 1) +
    ylim(0,NA) +
    xlim(0,NA) +
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
    geom_text(aes(label = paste0("D = ",format(round(dee,3), nsmall = 3)), x = min(msdmean, na.rm = TRUE), y = Inf), hjust = 0, vjust = 1)

  return(p)
}

