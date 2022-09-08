#' Calculate fractal dimension (FD)
#'
#' Calculate for each track, its fractal dimension
#' Katz & George (1985) define fractal dimension (D) as log(n) / (log(n) + log(d/L))
#' where n is the number of steps, d is the longest of all possible point-to-point distances and
#' L is the cumulative length of the track.
#' D is ~1 for directed trajectories, ~2 for confined and ~3 for subdiffusion
#' Here we calculate this and store D (called fd) and d (called wide) for return.
#' Note that this method does not take into account gaps in the track.
#' For a track with many gaps, n will be lowered.
#'
#' @param dataList list of a data frame (must include at a minimum - trace (track ID), x, y and frame (in real coords)) and a calibration data frame
#' @return data frame
#' @examples
#' xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
#' tmObj <- readTrackMateXML(XMLpath = xmlPath)
#' tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = 0.04)
#' fdDF <- calculateFD(dataList = tmObj)
#' @export

calculateFD <- function(dataList) {
  x <- y <- NULL

  if(inherits(dataList, "list")) {
    df <- dataList[[1]]
  } else {
    cat("Function requires a list of TrackMate data and calibration data\n")
    return(NULL)
  }

  # make a list of all traces (this is the filtered list of traces from TrackMate XML)
  traceList <- unique(df$trace)
  # data frame to hold results
  result <- data.frame(trace = rep("", length(traceList)),
                       wide = rep(NA, length(traceList)),
                       fd = rep(NA, length(traceList)))
  ii <- 1
  for (i in traceList){
    a <- df %>%
      filter(trace == i) %>%
      select(x, y)
    b <- dist(a)
    n <- nrow(a)
    d <- max(b, na.rm = TRUE)
    l <- max(df$cumulative_distance[df$trace == i])
    # fractal dimension
    fd <- log(n) / log( n * d * l^-1)
    result$trace[ii] <- i
    result$wide[ii] <- d
    result$fd[ii] <- fd
    ii <- ii + 1
  }

  return(result)
}
