#' Fitting jump distance (JD) data
#'
#' Jump Distances have been calculated for a given time lag.
#' They can be described by fitting curves to the data, either using a histogram or cumulative probability density function.
#' Firtting to a histogram is sensitive to binning parameters and ECDF performs better for general use.
#' The idea behind this analysis is given in:
#' - Weimann et al. (2013) A quantitative comparison of single-dye tracking analysis tools using Monte Carlo simulations. PloS One 8, e64287.
#' - Menssen & Mani (2019) A Jump-Distance-Based Parameter Inference Scheme for Particulate Trajectories, Biophysical Journal, 117: 1, 143-156.
#' This function is called from inside `makeSummaryReport()`, using the defaults.
#' However, you can pass additional arguments for `nPop`, `init` and `breaks` via the ellipsis.
#' Fitting is tricky and two populations (default) is written to catch errors and retry.
#' In the case of failure, try passing better guesses via `init`.
#' To fit 3 populations, you must pass `nPop = 3` as an additional argument, you are advised to also pass guesses via `init`.
#'
#' @param jumpList list of a data frame (with a column named jump of jump distances), and a list of params
#' @return ggplot
#' @examples
#' xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
#' tmObj <- readTrackMateXML(XMLpath = xmlPath)
#' tmObj <- correctTrackMateData(tmObj, xyscalar = 0.04)
#' jdObj <- calculateJD(dataList = tmObj, deltaT = 2)
#' fittingJD(jumpList = jdObj)
#' @export

fittingJD <- function(jumpList) {
  coef <- counts <- countsCum <- hist <- mid <- nls <- value <- variable <- NULL

  if(inherits(jumpList, "list")) {
    df <- jumpList[[1]]
    params <- jumpList[[2]]
  } else {
    cat("Function requires a list of jump distance data and some parameters\n")
    return(NULL)
  }

  nPop <- mode <- init <- units <- timeRes <- breaks <- NULL
  nPop <- params$nPop
  mode <- params$mode
  init <- params$init
  units <- params$units
  timeRes <- params$timeRes
  breaks <- params$breaks

  if(nPop < 1 | nPop > 3) {
    return(NULL)
  }

  if(is.null(init)) {
    if(mode == "ECDF") {
      guess <- ((mean(df$jump, na.rm = T) / timeRes) / 4) / breaks
      # if(nPop == 1) {init <- list(D1 = 0.05)}
      if(nPop == 1) {init <- list(D1 = guess)}
      # if(nPop == 2) {init <- list(D1 = 0.001, D2 = 0.4, D3 = 0.1)}
      if(nPop == 2) {init <- list(D1 = guess, D2 = 0.4, D3 = guess*10)}
      # if(nPop == 3) {init <- list(D1 = 0.1, D2 = 0.8, D3 = 0.01, D4 = 0.2, D6 = 0.2)}
      if(nPop == 3) {init <- list(D1 = guess, D2 = guess*10, D3 = guess*100, D4 = 0.2, D6 = 0.2)}
    } else {
      if(nPop == 1) {init <- list(D2 = 100,  D1 = 0.1)}
      if(nPop == 2) {init <- list(D2 = 0.01,  D1 = 0.1, D3 = 10, D4 = 100)}
      if(nPop == 3) {init <- list(D2 = 1,  D1 = 0.1, D3 = 0.01, D4 = 10, D5 = 100, D6 = 30)}
    }
  }
  # make histogram
  hd <- hist(df$jump, breaks = breaks, plot = FALSE)
  hdata <- data.frame(counts = hd$counts,
                      mid = hd$mids,
                      countsCum = cumsum(hd$counts) / sum(hd$counts) )
  xStr <- paste0("Displacement (",units[1],")")
  if(mode == "ECDF") {
    p <- ggplot(data = hdata, aes(x = mid, y = countsCum)) +
      geom_point() +
      lims(x = c(0, max(hdata$mid)), y = c(0, 1.4)) +
      labs(x = xStr, y = "Frequency") +
      theme_classic() +
      theme(legend.position = "none")
  } else {
    p <- ggplot(data = hdata, aes(x = mid, y = counts)) +
      geom_col(fill = "light grey", colour = "dark grey") +
      lims(x = c(0, 1.1 * max(hdata$mid)), y = c(0, 1.4 * max(hdata$counts))) +
      labs(x = xStr, y = "Counts") +
      theme_classic() +
      theme(legend.position = "none")
  }

  # fitting
  x <- seq(0, max(hdata$mid), length = 10 * length(hdata$mid))

  if(nPop == 1) {

    if(mode == "ECDF") {
      fitc <- NULL
      try(fitc <- nls(hdata$countsCum ~
                    1 - exp(-hdata$mid^2 / (4 * D1 * timeRes)),
                  data = hdata,
                  start = init),
          silent = TRUE);
      if(is.null(fitc)) {
        cat("Failed to fit jump distances with 1 population. Try using different parameters for `init` and/or `nPop`.\n")
        return(p)
      }
      y <- 1 - exp(-x^2 / (4 * coef(fitc)[1] * timeRes))
      fitStr <- paste0("D = ",format(round(coef(fitc)[1],4), nsmall = 4))
    } else {
      fit <- NULL
      try(fit <- nls(hdata$counts ~
                   D2 * hdata$mid/(2 * D1 * timeRes) * exp(-hdata$mid^2 / (4 * D1 * timeRes)),
                 data = hdata,
                 start = init),
          silent = TRUE);
      if(is.null(fit)) {
        cat("Failed to fit jump distances with 3 populations. Try using different parameters for `init` and/or `nPop`.\n")
        return(p)
      }
      y <-(coef(fit)[[1]]) * x / (2 * (coef(fit)[[2]]) * timeRes) * exp(-x^2 / (4 * (coef(fit)[[2]]) * timeRes))
      fitStr <- paste0("D = ",format(round(coef(fit)[2],4), nsmall = 4),"\n","A = ", format(round(coef(fit)[1],4), nsmall = 4))
    }
    fitdf <- data.frame(x = x, y = y)
  }

  if(nPop == 2) {

    if(mode == "ECDF") {
      fitc2 <- NULL
      try(fitc2 <- nls(hdata$countsCum ~
                         (1 - D2 * exp(-hdata$mid^2 / (4 * D1 * timeRes)) -
                            (1 - D2) * exp(-hdata$mid^2 / (4 * D3 * timeRes))),
                       data = hdata,
                       start = init),
          silent = TRUE);
      if(is.null(fitc2)) {
        cat("Failed to fit jump distances with 2 populations. Try using different parameters for `init` and/or `nPop`.\n")
        return(p)
      }
      y1 <- coef(fitc2)[2] - (coef(fitc2)[2] * exp(-x^2 / (4 * coef(fitc2)[1] * timeRes)))
      y2 <- (1 - coef(fitc2)[2]) - (1 - coef(fitc2)[2]) * exp(-x^2 / (4 * coef(fitc2)[3] * timeRes))
      fitStr <- paste0("D1 = ",format(round(coef(fitc2)[1],4), nsmall = 4),"\n",
                       "A1 = ", format(round(coef(fitc2)[2],4), nsmall = 4),"\n",
                       "D2 = ",format(round(coef(fitc2)[3],4), nsmall = 4),"\n",
                       "A2 = ", format(round(1 - coef(fitc2)[2],4), nsmall = 4))
    } else {
      fit2 <- NULL
      try(fit2 <- nls(hdata$counts ~
                   D3 * hdata$mid / (2 * D1 * timeRes) * exp(-hdata$mid^2 / (4 * D1 * timeRes)) +
                   D4 * hdata$mid / (2 * D2 * timeRes) * exp(-hdata$mid^2 / (4 * D2 * timeRes)),
                 data = hdata,
                 start = init),
          silent = TRUE);
      if(is.null(fit2)) {
        cat("Failed to fit jump distances with 2 populations. Try using different parameters for `init` and/or `nPop`.\n")
        return(p)
      }
      y1 <- (coef(fit2)[[4]]) * x / (2 * (coef(fit2)[[1]]) * timeRes) * exp(-x^2 / (4 * (coef(fit2)[[1]]) * timeRes))
      y2 <- (coef(fit2)[[3]]) * x / (2 * (coef(fit2)[[2]]) * timeRes) * exp(-x^2 / (4 * (coef(fit2)[[2]]) * timeRes))
      fitStr <- paste0("D1 = ",format(round(coef(fit2)[1],4), nsmall = 4),"\n",
                       "A1 = ", format(round(coef(fit2)[4],4), nsmall = 4),"\n",
                       "D2 = ",format(round(coef(fit2)[2],4), nsmall = 4),"\n",
                       "A2 = ", format(round(coef(fit2)[3],4), nsmall = 4))
    }
    y <- y1 + y2
    fitdf <- data.frame(x = x, y = y, y1 = y1, y2 = y2)
  }

  if(nPop == 3) {

    if(mode == "ECDF") {
      fitc3 <- NULL
      try(fitc3 <- nls(hdata$countsCum ~
                     (1 - D4 * exp(-hdata$mid^2 / (4 * D1 * timeRes)) -
                        (1-D4-D6) * exp(-hdata$mid^2 / (4 * D2 * timeRes)) -
                        D6 * exp(-hdata$mid^2 / (4 * D3 * timeRes))),
                   data = hdata,
                   start = init),
          silent = TRUE);
      if(is.null(fitc3)) {
        cat("Failed to fit jump distances with 3 populations. Try using different parameters for `init` and/or `nPop`.\n")
        return(p)
      }
      d5 <- 1 - (coef(fitc3)[4] + coef(fitc3)[5])
      y1 <- coef(fitc3)[4] - coef(fitc3)[4] * exp(-x^2 / (4 * coef(fitc3)[1] * timeRes))
      y2 <- d5 - d5 * exp(-x^2 / (4 * coef(fitc3)[2] * timeRes))
      y3 <- coef(fitc3)[5] - coef(fitc3)[5] * exp(-x^2 / (4 * coef(fitc3)[3] * timeRes))
      fitStr <- paste0("D1 = ",format(round(coef(fitc3)[1],4), nsmall = 4),"\n",
                       "A1 = ", format(round(coef(fitc3)[4],4), nsmall = 4),"\n",
                       "D2 = ",format(round(coef(fitc3)[2],4), nsmall = 4),"\n",
                       "A2 = ", format(round(d5,4), nsmall = 4),"\n",
                       "D3 = ",format(round(coef(fitc3)[3],4), nsmall = 4),"\n",
                       "A3 = ", format(round(coef(fitc3)[5],4), nsmall = 4),"\n")
    } else {
      fit3 <- NULL
      try(fit3 <- nls(hdata$counts ~
                    D4 * hdata$mid / (2 * D1 * timeRes) * exp(-hdata$mid^2 / (4 * D1 * timeRes)) +
                    D5 * hdata$mid / (2 * D2 * timeRes) * exp(-hdata$mid^2 / (4 * D2 * timeRes)) +
                    D6 * hdata$mid / (2 * D3 * timeRes) * exp(-hdata$mid^2 / (4 * D3 * timeRes)),
                  data = hdata,
                  start = init),
          silent = TRUE);
      if(is.null(fit3)) {
        cat("Failed to fit jump distances with 3 populations. Try using different parameters for `init` and/or `nPop`.\n")
        return(p)
      }
      y1 <- (coef(fit3)[[4]]) * x / (2 * (coef(fit3)[[1]]) * timeRes) * exp(-x^2 / (4 * (coef(fit3)[[1]]) * timeRes))
      y2 <- (coef(fit3)[[5]]) * x / (2 * (coef(fit3)[[2]]) * timeRes) * exp(-x^2 / (4 * (coef(fit3)[[2]]) * timeRes))
      y3 <- (coef(fit3)[[6]]) * x / (2 * (coef(fit3)[[3]]) * timeRes) * exp(-x^2 / (4 * (coef(fit3)[[3]]) * timeRes))
      fitStr <- paste0("D1 = ",format(round(coef(fit3)[1],4), nsmall = 4),"\n",
                       "A1 = ", format(round(coef(fit3)[4],4), nsmall = 4),"\n",
                       "D2 = ",format(round(coef(fit3)[2],4), nsmall = 4),"\n",
                       "A2 = ", format(round(coef(fit3)[5],4), nsmall = 4),"\n",
                       "D3 = ",format(round(coef(fit3)[3],4), nsmall = 4),"\n",
                       "A3 = ", format(round(coef(fit3)[6],4), nsmall = 4))
    }
    y <- y1 + y2 + y3
    fitdf <- data.frame(x = x, y = y, y1 = y1, y2 = y2, y3 = y3)
  }
  # first melt so that we can plot all fit curves regardless of how many there are
  fitdf <- melt(fitdf, id.vars="x")
  # add fit lines to ggplot
  p <- p + geom_line(data = fitdf, aes(x = x, y = value, col = variable))
  # add fitting coefficients
  if(mode == "ECDF") {
    p <- p + geom_text(aes(label = fitStr, x = 0, y = Inf), size = 2, hjust = 0, vjust = 1, check_overlap = TRUE)
    p <- p + geom_text(aes(label = paste(timeRes,units[2]), x = 0, y = Inf), size = 2, hjust = 1, vjust = 0, check_overlap = TRUE)
  } else {
    p <- p + geom_text(aes(label = fitStr, x = 0, y = Inf), size = 2, hjust = 0, vjust = 1, check_overlap = TRUE)
    p <- p + geom_text(aes(label = paste(timeRes,units[2]), x = Inf, y = 0), size = 2, hjust = 1, vjust = 1, check_overlap = TRUE)
  }

  return(p)
}
