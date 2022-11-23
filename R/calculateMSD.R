#' Calculate Mean Squared Displacement (MSD)
#'
#' Calculation of the MSD of multiple tracks.
#' There are two methods for everaging MSD data from multiple tracks:
#' ensemble = for each time lag average all squared displacements from all tracks
#' time-averaged = find MSD for each track and then generate the average MSD from these curves
#' The MSD curves will be identical if all tracks are the same length, and diverge if not.
#' Standard deviation will be large for ensemble and smaller for time-averaged data.
#' Input is a data frame of tracks imported using readTrackMateXML()
#'
#' @param df data frame must include at a minimum - trace (track ID), x, y and t (in real coords)
#' @param method string. Either "ensemble" or "timeaveraged" (default)
#' @param N numeric variable for MSD. dt should be up to 1/N of number of data points (4 recommended)
#' @param short numeric variable for the shortest number of points we will analyse. Note, this uses the number of frames from start, not number of points in track, i.e. a track with <short points and many gaps will remain
#' @return list of two data frames
#' @examples
#' xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
#' tmObj <- readTrackMateXML(XMLpath = xmlPath)
#' tmObj <- correctTrackMateData(tmObj, xyscalar = 0.04)
#' tmDF <-  tmObj[[1]]
#' msdObj <- calculateMSD(df = tmDF, method = "ensemble", N = 3, short = 8)
#' @export

calculateMSD <- function(df, method = "timeaveraged", N = 4, short = 0) {
  x <- y <- NULL

  if(!inherits(df, "data.frame")) {
    cat("Function requires a data frame.\n")
    return(NULL)
  }

  # make a list of all traces (this is the filtered list of traces from TrackMate XML)
  traceList <- unique(df$trace)
  tList <- list()
  # for each trace find the frame offset
  for (i in traceList) {
    a <- df %>%
      filter(trace == i) %>%
      select(frame, x)
    frame0 <- a$frame[1]
    a <- a[-1,]
    tList[[i]] <-  a$frame - frame0
  }
  # filter out short traces
  if(short > 0) {
    traceList <- subset(traceList,unlist(lapply(tList, max)) >= short)
    tList <- subset(tList,unlist(lapply(tList, max)) >= short)
  }
  # find the maximum frame offset of all traces
  tListmax <- max(unlist(lapply(tList, max)))
  # create matrices for x and y coords one row for each frame, each trace in a column
  displacementXMat <- matrix(data = NA, nrow = tListmax, ncol = length(traceList))
  colnames(displacementXMat) <- traceList
  rownames(displacementXMat) <- 1:tListmax
  displacementYMat <- displacementXMat
  # now we repeat the process to extract the x y coords for each frame offset, store in matrices.
  for (i in traceList){
    a <- df %>%
      filter(trace == i) %>%
      select(x, y, frame)
    frame0 <- a$frame[1]
    a <- a[-1,]
    a$frame <- a$frame - frame0
    displacementXMat[a$frame,i] <- a$x
    displacementYMat[a$frame,i] <- a$y
  }
  # find the 90th centile longest trace (in terms of frame offset)
  tListmax2 <- quantile(unlist(lapply(tList, max)), probs=.9)
  # dt should be up to 1/4 of number of data points (default)
  numberOfdeltaT = floor(tListmax2/N)
  # make matrix to store the averaged msd summary
  msd <- matrix(data = NA, nrow = numberOfdeltaT, ncol = 5 )
  colnames(msd) <- c("mean", "sd", "n", "size", "t")
  tstep <- df$t[match(1,df$frame)]
  # make matrix to store the msd curve for each track
  trackmsd <- matrix(data = NA, nrow = numberOfdeltaT, ncol = length(traceList) )
  colnames(trackmsd) <- traceList

  for(deltaT in 1 : numberOfdeltaT){
    # calculate displacement between points for x and y
    deltaXCoords <- displacementXMat[(1 + deltaT) : tListmax,] - displacementXMat[1 : (tListmax-deltaT),]
    deltaYCoords <- displacementYMat[(1 + deltaT) : tListmax,] - displacementYMat[1 : (tListmax-deltaT),]
    # calculate squared displacement
    squaredDisplacement <- deltaXCoords**2 + deltaYCoords**2
    # we collect the MSD for each track (column)
    if(length(dim(squaredDisplacement)) < 2) {
      next
    }
    eachSquaredDisplacement <- colMeans(squaredDisplacement, na.rm = TRUE)
    # summary statistics for each deltaT
    if(method == "ensemble") {
      # store the ensemble data
      msd[deltaT,1] <- mean(squaredDisplacement, na.rm = TRUE) # average
      msd[deltaT,2] <- sd(squaredDisplacement, na.rm = TRUE) # standard deviation
      msd[deltaT,3] <- sum(!is.na(squaredDisplacement)) # n
    } else {
      # store the time-averaged data
      msd[deltaT,1] <- mean(eachSquaredDisplacement, na.rm = TRUE) # average
      msd[deltaT,2] <- sd(eachSquaredDisplacement, na.rm = TRUE) # standard deviation
      msd[deltaT,3] <- sum(!is.na(eachSquaredDisplacement)) # n
    }
    # place the MSD for this deltaT and for every track into the matrix
    trackmsd[deltaT,] <- eachSquaredDisplacement
  }
  # send msd curves for each track to compute the diffusive behaviour
  alphas <- calculateAlpha(trackmsd, tstep)

  # format msd matrix into data frame
  msd <- as.data.frame(msd)
  msd$size <- c(1 : numberOfdeltaT)
  msd$t <- msd$size * tstep

  msdList <- list(msd,alphas)

  return(msdList)
}

