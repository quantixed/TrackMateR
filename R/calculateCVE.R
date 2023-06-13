#' Calculate CVE (covariance-based estimator)
#'
#' From Vestergaard et al., Physical Review E (2019) 89, 022726
#' Input is two data matrices of X and Y displacements for each trace.
#' Output is estimator of D and sigma for each trace.
#'
#' @param xMat matrix of x displacements, each col is a track, each row is time (will contain NAs)
#' @param yMat matrix of y displacements, each col is a track, each row is time (will contain NAs)
#' @param tstep variable. Time step in seconds
#' @return data frame
#' @export


calculateCVE <- function(xMat, yMat, tstep) {
  # check that data is at least four rows by two columns
  if(nrow(xMat) < 4 | ncol(xMat) < 2) {
    cveDF <- data.frame(trace = character("1"),
                          cve = numeric(1))
    return(cveDF)
  }
  # calculation of the estimators of D and sigma^2
  estDx <- colMeans(xMat^2, na.rm = TRUE) / 2 / tstep + colMeans(xMat[-1,] * xMat[-(nrow(xMat)),], na.rm = TRUE) / tstep
  estDy <- colMeans(yMat^2, na.rm = TRUE) / 2 / tstep + colMeans(yMat[-1,] * xMat[-(nrow(yMat)),], na.rm = TRUE) / tstep
  estD <- rowMeans(cbind(estDx,estDy))
  R <- 1/6 # R: motion blur coefficient, representing the effect of motion blur. It varies [0, 1/4], R=0 in case of no blur. R=1/6  if the camera shutter is kept open for the full duration of the time-lapse.
  estSigma2x <- R * colMeans(xMat^2, na.rm = TRUE) + (2*R-1) * colMeans(xMat[-1,] * xMat[-(nrow(xMat)),], na.rm = TRUE)
  estSigma2y <- R * colMeans(yMat^2, na.rm = TRUE) + (2*R-1) * colMeans(yMat[-1,] * xMat[-(nrow(yMat)),], na.rm = TRUE)
  estSigma2 <- rowMeans(cbind(estSigma2x, estSigma2y))

  # turn into data frame
  cveDF <- data.frame(dee = estD,
                      estsigma = sqrt(estSigma2))

  return(cveDF)
}

