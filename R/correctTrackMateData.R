#' Correct distance and time of imported TrackMate data.
#'
#' If the TrackMate data is in pixels and/or frames, the data frame can be converted with this function.
#'
#' @param df data frame of imported track mate data
#' @param xyscalar numeric multiplier to correct pixel size of original movie. Assumes isotropic scaling, i.e. pixel height = pixel width
#' @param tscalar numeric multiplier to correct frame interval of original movie. Frame interval of tracked data.
#' @return data frame
#' @examples
#' xmlPath <- "~/Desktop/FakeTracks.xml"
#' datalist <- readTrackMateXML(XMLpath = xmlPath)
#' data <-  datalist[[1]]
#' # in the case where pixel size is 0.03 um and original data is 1 pixel, xyscalar = 0.03
#' data <- correctTrackMateData(df = data, xyscalar = 0.03)
#' @export

# correcting dataframe - assumes that scaling is currently 1 pixel and/or 1 frame
correctTrackMateData <- function(df, xyscalar = 1, tscalar = 1) {
  if(xyscalar == 1 & tscalar == 1) {
    cat("No correction applied.\n")
    return(df)
  } else {
    msg <- ""
  }
  if(xyscalar != 1) {
    # change pixel size
    df$x <- df$x * xyscalar
    df$y <- df$y * xyscalar
    df$displacement <- df$displacement * xyscalar
    df$cumulative_distance <- df$cumulative_distance * xyscalar
    msg <- paste0("Correcting XY scale. ",msg)
  }
  if(tscalar != 1) {
    df$t <- df$t * tscalar
    df$track_duration <- df$track_duration * tscalar
    msg <- paste0("Correcting timescale. ",msg)
  }
  msg <- paste0(msg,"\n")
  cat(msg)
  tstep <- df$t[match(1,df$frame)]
  df$speed <- df$displacement / tstep

  return(df)
}
