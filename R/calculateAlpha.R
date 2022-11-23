#' Calculate alpha (relationship between MSD and normal diffusion)
#'
#' Normal diffusion is alpha = 1. Subdiffusion is alpha < 1 and superdiffusion is alpha > 1.
#' Input is a data matrix of msd curves.
#' Output is mean of log2(alpha), one value for each trace.
#'
#' @param alphaMat matrix of msd curves, each col is a track, each row is time lag (will contain NAs)
#' @param tstep variable. Time step in seconds
#' @return data frame
#' @export


calculateAlpha <- function(alphaMat,tstep) {
  # check that alphaMat is at least four rows by two columns
  if(nrow(alphaMat) < 4 | ncol(alphaMat) < 2) {
    alphaDF <- data.frame(trace = character("1"),
                          alpha = numeric(1))
    return(alphaDF)
  }
  # make time vector
  tee <- (1 : nrow(alphaMat)) * tstep
  # make a vector for the results
  alphaVec <- numeric(ncol(alphaMat))
  alphaVec[] <- NA
  # check that we have four contiguous points for each col
  check <- colSums(alphaMat[1:4,])
  for(i in 1 : ncol(alphaMat)) {
    if(is.na(check[i])) {
      next
    }
    tempdf <- data.frame(mean = alphaMat[,i],
                         t = tee)
    # fit to first four data points
    if(all(is.na(tempdf$mean)) | all(is.na(tempdf$t))) {
      next
    }
    mod <- lm(mean ~ t, data = tempdf[1:4,])
    # make a column containing model y points for each t
    tempdf$pred <- (mod$coefficients[2] * tempdf$t) + mod$coefficients[1]
    tempdf$alpha <- tempdf$mean / tempdf$pred
    tempdf$alpha <- suppressWarnings(log2(tempdf$alpha))
    alphaVec[i] <- mean(tempdf$alpha, na.rm = TRUE)
  }
  # turn into data frame
  alphaDF <- data.frame(trace = colnames(alphaMat),
                        alpha = alphaVec)

  return(alphaDF)
}

