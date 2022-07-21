#' Calculate alpha (relationship between MSD and normal diffusion)
#'
#' Normal diffusion is alpha = 1. Subdiffusion is alpha < 1 and superdiffusion is alpha > 1.
#'
#' @param alphaMat matrix of msd curves, each col is a track, each row is time lag (will contain NAs)
#' @param tstep variable. Time step in seconds
#' @return numeric vector
#' @export


calculateAlpha <- function(alphaMat,tstep) {
  # make time vector
  tee <- (1 : nrow(alphaMat)) * tstep
  # make vector for the results
  alphaVec <- numeric(ncol(alphaMat))
  alphaVec[] <- NA
  # check that we have four contiguous points for each col
  check <-colSums(alphaMat[1:4,])
  for(i in 1 : ncol(alphaMat)) {
    if(is.na(check[i])) {
      next
    }
    tempdf <- data.frame(mean = alphaMat[,i],
                         t = tee)
    # fit to first four data points
    mod <- lm(mean ~ t, data = tempdf[1:4,])
    # make a column containing model y points for each t
    tempdf$pred <- (mod$coefficients[2] * tempdf$t) + mod$coefficients[1]
    tempdf$alpha <- tempdf$pred / tempdf$mean
    alphaVec[i] <- mean(tempdf$alpha, na.rm = TRUE)
  }

  return(alphaVec)
}

