#' Calculate Jump Distance (JD)
#'
#' Calculation of the JD of multiple tracks.
#' Calculation is equivalent to a single time lag point on the ensemble MSD curve, typically represented as a histogram
#' Input is a data frame of tracks imported using readTrackMateXML()
#' The default time step is one frame - which is the equivalent to the plot generated to show displacement versus time.
#'
#' @param dataList list of data frame (must include at a minimum - trace (track ID), x, y and t (in real coords)) and calibration
#' @param deltaT integer to represent the multiple of frames that are to be analysed
#' @return a list of data frame of jump distances, NAs removed; and a numeric variable (jumptime)
#' @examples
#' xmlPath <- "~/Desktop/FakeTracks.xml"
#' tmObj <- readTrackMateXML(XMLpath = xmlPath)
#' tmObj <- correctTrackMateData(tmObj, xyscalar = 0.04)
#' jdObj <- calculateJD(dataList = tmObj, deltaT = 2)
#' @export

calculateJD <- function(dataList, deltaT = 1) {
  x <- y <- NULL

  df <- dataList[[1]]
  calibration <- dataList[[2]]

  if(deltaT < 1 | deltaT %% 1 != 0) {
    # check is an integer >= 1
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
  # find the maximum frame offset of all traces
  tListmax <- max(unlist(lapply(tList, max)))
  if(deltaT > tListmax) {
    # break because requested interval would result in no data
    return(NULL)
  }
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

  # calculate displacement between points for x and y
  deltaXCoords <- displacementXMat[(1 + deltaT) : tListmax,] - displacementXMat[1 : (tListmax-deltaT),]
  deltaYCoords <- displacementYMat[(1 + deltaT) : tListmax,] - displacementYMat[1 : (tListmax-deltaT),]
  # calculate distance
  jdMat <- sqrt(deltaXCoords**2 + deltaYCoords**2)
  # convert to vector
  jdVec <- c(jdMat)
  jdVec <- jdVec[!is.na(jdVec)]
  # convert to dataframe with column called jump
  jdDF <- data.frame(jump = jdVec)
  # record jumptime which is deltaT * time resoltuion
  jumptime <- deltaT * calibration[2,1]
  # make list of these two and return
  jumpList <- list(jdDF, jumptime)

  return(jumpList)
}

