#' Correct distance and time of imported TrackMate data.
#'
#' If the TrackMate data is in pixels and/or frames, the data frame can be converted with this function.
#'
#' @param dataList a list of a data frame (of track data) and a calibration data frame (from TrackMate XML)
#' @param xyscalar numeric multiplier to correct pixel size of original movie. Assumes isotropic scaling, i.e. pixel height = pixel width
#' @param tscalar numeric multiplier to correct frame interval of original movie. Frame interval of tracked data.
#' @param xyunit string to describe spatial units
#' @param tunit string to describe temporal unit
#' @return list of two data frames
#' @examples
#' xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
#' tmObj <- readTrackMateXML(XMLpath = xmlPath)
#' # in the case where pixel size is 0.03 um and original data is 1 pixel, xyscalar = 0.03
#' tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = 0.03)
#' @export

# correcting dataframe - assumes that scaling is currently 1 pixel and/or 1 frame
correctTrackMateData <- function(dataList, xyscalar = 1, tscalar = 1, xyunit = NULL, tunit = NULL) {
  if(xyscalar == 1 & tscalar == 1) {
    cat("No correction applied.\n")
    return(dataList)
  } else {
    msg <- ""
  }
  if(!inherits(dataList, "list")) {
    cat("Requires TrackMate data frame and calibration data frame, as a list.\n")
    return(dataList)
  }
  df <- dataList[[1]]
  calib <- dataList[[2]]
  if(xyscalar != 1) {
    msg <- paste0("Correcting XY scale. ",msg)
    # change pixel size
    df$x <- df$x * xyscalar
    df$y <- df$y * xyscalar
    df$displacement <- df$displacement * xyscalar
    df$cumulative_distance <- df$cumulative_distance * xyscalar
    df$radius <- df$radius * xyscalar
    # there are potentially other columns present that may need adjusting
    if("area" %in% colnames(df))
    {
      df$area <- df$area * xyscalar^2
    }
    if("perimeter" %in% colnames(df))
    {
      df$perimeter <- df$perimeter * xyscalar
    }
    if("ellipse_x0" %in% colnames(df))
    {
      df$ellipse_x0 <- df$ellipse_x0 * xyscalar
      df$ellipse_y0 <- df$ellipse_y0 * xyscalar
      df$ellipse_major <- df$ellipse_major * xyscalar
      df$ellipse_minor <- df$ellipse_minor * xyscalar
    }
    # correct calibration data
    calib[1,1] <- calib[1,1] * xyscalar
    calib[3:4,1] <- calib[3:4,1] * xyscalar # border of image
  }
  if(tscalar != 1) {
    msg <- paste0("Correcting timescale. ",msg)
    # change time scale
    df$t <- df$t * tscalar
    df$track_duration <- df$track_duration * tscalar
    # correct calibration data
    calib[2,1] <- calib[2,1] * tscalar
  }
  # will change units even if no scaling is done
  if(!is.null(xyunit)) {
    calib[1,2] <- xyunit
  }
  if(!is.null(tunit)) {
    calib[2,2] <- tunit
  }
  msg <- paste0(msg,"\n")
  cat(msg)
  tstep <- df$t[match(1,df$frame)]
  df$speed <- df$displacement / tstep

  dataList <- list(df,calib)

  return(dataList)
}
