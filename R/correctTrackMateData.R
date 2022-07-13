#' Correct distance and time of imported TrackMate data.
#'
#' If the TrackMate data is in pixels and/or frames, the data frame can be converted with this function.
#'
#' @param df data frame of imported track mate data
#' @param xysize pixel size of original movie. Assumes isotropic scaling, i.e. pixel height = pixel width
#' @param tsize time. Frame interval of tracked data.
#' @return data frame
#' @examples
#' xmlPath <- "~/Desktop/FakeTracks.xml"
#' data <- readTrackMateXML(XMLpath = xmlPath)
#' data <- correctTrackMateData(data, xy = 0.03)
#' @export

# correcting dataframe - assumes that scaling is currently 1 pixel and/or 1 frame
correctTrackMateData <- function(df, xysize = 1, tsize = 1) {
  if(xysize == 1 & tsize == 1) {
    return(df)
  }
  if(xysize != 1) {
    # change pixel size
    df$x <- df$x * xysize
    df$y <- df$y * xysize
    df$displacement <- df$displacement * xysize
    df$cumulative_distance <- df$cumulative_distance * xysize
  }
  if(tsize != 1) {
    df$t <- df$t * tsize
    df$track_duration <- df$track_duration * tsize
  }
  tstep <- df$t[match(1,df$frame)]
  df$speed <- df$displacement / tstep

  return(df)
}
