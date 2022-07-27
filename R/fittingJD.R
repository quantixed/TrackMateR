#' Fitting jump distance (JD) data
#'
#' Jump Distances have been calculated for a given time lag.
#' They can be described by fitting curves to the data, either using a histogram or cumulative probability density function.
#' Firtting to a histogram is sensitive to binning parameters and ECDF performs better for general use.
#' The idea behind this analysis is given in:
#' - Weimann et al. (2013) A quantitative comparison of single-dye tracking analysis tools using Monte Carlo simulations. PloS One 8, e64287.
#' - Menssen & Mani (2019) A Jump-Distance-Based Parameter Inference Scheme for Particulate Trajectories, Biophysical Journal, 117: 1, 143-156.
#' The bulk of this code is taken from trackR by JuG

#'
#' @param df data frame with a column named jump of jump distances
#' @param mode string indicated ECDF (default) or hist (histogram)
#' @param nPop number of populations of diffusing species (1 or 2)
#' @param init initialisation parameters for the nls fit for example list(D2 = 200,  D1 = 0.1) or list(D2 = 0.01,  D1=0.1, D3=10, D4=100)
#' @param units character vector to describe units (defaults are um, micrometres and  s, seconds)
#' @param timeRes time resolution per unit of jump. Frame interval is 0.5 s and jump interval is two steps, timeRes = 1.
#' @param breaks number of bins for histogram. With ECDF breaks can be high e.g. 100, for mode = "hist" they should be low, perhaps 30.
#' @return ggplot
#' @examples
#' xmlPath <- "~/Desktop/FakeTracks.xml"
#' datalist <- readTrackMateXML(XMLpath = xmlPath)
#' data <-  datalist[[1]]
#' data <- correctTrackMateData(data, xy = 0.04)
#' jdvec <- calculateJD(data, deltaT = 2)
#' jdDF <- data.frame(jump = jdvec)
#' fittingJD(df = jdDF, mode = "ECDF", nPop = 2, breaks = 100, timeRes = 0.06)
#' @export

fittingJD <- function(df, mode = "ECDF", nPop = 1, init, units = c("um","s"), timeRes = 1, breaks = 150, ...) {

  if(nPop < 1 | nPop > 2) {
    return()
  }
  if(missing(init)) {
    if(mode == "ECDF") {
      if(nPop == 1) {init = list(D1 = 0.05)}
      if(nPop == 2) {init = list(D1 = 0.001, D2 = 0.4, D3 = 0.1)}
    } else {
      if(nPop == 1) {init = list(D2 = 100,  D1 = 0.1)}
      if(nPop == 2) {init = list(D2 = 0.01,  D1 = 0.1, D3 = 10, D4 = 100)}
    }
  }
  # make histogram
  hd <- hist(df$jump, breaks = breaks, plot = FALSE)
  hdata <- data.frame(counts = hd$counts,
                      mid = hd$mids,
                      countsCum = cumsum(hd$counts) / sum(hd$counts) )
  xStr <- paste0("Displacement (",units[1],")")
  if(mode == "ECDF") {
    plot(hd$mid, hdata$countsCum, ylim = c(0, 1.1),
         xlab = xStr, ylab = "Frequency", cex = 0.75)
  } else {
    hist(df$jump, breaks = breaks, xlim = c(0, 1.1 * max(hdata$mid)), ylim = c(0, 1.5 * max(hdata$counts)),
         xlab = xStr, panel.first=grid())
  }

  # fitting
  x <- seq(0, max(hdata$mid), length = 10 * length(hdata$mid))

  if(nPop == 1) {

    if(mode == "ECDF") {
      fitc <- nls(hdata$countsCum ~ 1 - exp(-hdata$mid^2 / (4 * D1 * timeRes)),
                  data = hdata,
                  start = init)
      print(fitc)
      print(confint(fitc))
      y <- 1 - exp(-x^2 / (4 * coef(fitc)[1] * timeRes))
      lines(x, y, col = "red", lwd = 2)
    } else {
      fit <- nls(hdata$counts ~ D2 * hdata$mid/(2 * D1 * timeRes) * exp(-hdata$mid^2 / (4 * D1 * timeRes)),
                 data = hdata,
                 start = init)
      print(fit)
      print(suppressMessages(confint(fit)))

      y <-(coef(fit)[[1]]) * x / (2 * (coef(fit)[[2]]) * timeRes) * exp(-x^2 / (4 * (coef(fit)[[2]]) * timeRes))
      lines(x, y, col = "red",lwd = 2)
      box()
    }
  }

  if(nPop == 2) {

    if(mode == "ECDF") {
      fitc2 <- nls(hdata$countsCum ~ (1 - D2 * exp(-hdata$mid^2 / (4 * D1 * timeRes)) - (1-D2) * exp(-hdata$mid^2 / (4 * D3 * timeRes))),
                  data = hdata,
                  start = init)
      print(fitc2)
      # y <- (1 - coef(fitc2)[2] * exp(-x^2 / (4 * coef(fitc2)[1] * timeRes)) - (1 - coef(fitc2)[2]) * exp(-x^2 / (4 * coef(fitc2)[3] * timeRes)))
      y1 <- coef(fitc2)[2]-(coef(fitc2)[2] * exp(-x^2 / (4 * coef(fitc2)[1] * timeRes)))
      y2 <- (1 - coef(fitc2)[2]) - (1 - coef(fitc2)[2]) * exp(-x^2 / (4 * coef(fitc2)[3] * timeRes))
    } else {
      fit2<- nls(hdata$counts ~ D3 * hdata$mid / (2 * D1 * timeRes) * exp(-hdata$mid^2 / (4 * D1 * timeRes)) + D4 * hdata$mid / (2 * D2 * timeRes) * exp(-hdata$mid^2 / (4 * D2 * timeRes)),
                 data = hdata,
                 start = init)
      print(fit2)
      print(suppressMessages(confint(fit2)))
      # y <- (coef(fit2)[[4]]) * x / (2 * (coef(fit2)[[1]]) * timeRes) * exp(-x^2 / (4 * (coef(fit2)[[1]]) * timeRes)) +
      #   (coef(fit2)[[3]]) * x / (2 * (coef(fit2)[[2]]) * timeRes) * exp(-x^2 / (4 * (coef(fit2)[[2]]) * timeRes))
      y1 <- (coef(fit2)[[4]]) * x / (2 * (coef(fit2)[[1]]) * timeRes) * exp(-x^2 / (4 * (coef(fit2)[[1]]) * timeRes))
      y2 <- (coef(fit2)[[3]]) * x / (2 * (coef(fit2)[[2]]) * timeRes) * exp(-x^2 / (4 * (coef(fit2)[[2]]) * timeRes))
    }
    y <- y1 + y2
    lines(x , y, col = "red", lwd = 2)
    lines(x , y1, col = "blue", lwd = 2)
    lines(x , y2, col = "green", lwd = 2)
  }


  return()
}
