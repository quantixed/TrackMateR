#' Compare ground truth datasets
#'
#' Requires Ground Truth csv files to be organised into subfolders named according to the condition.
#' Outputs are saved to `Output/Plots/GT/` in the working directory.
#'
#' @param path string, filepath to directory equivalent to Data in compareDatasets()
#' @param xyscale numeric, pixel size (ground truth data is in pixel units)
#' @param xyunit string, spatial units
#' @param tscale numeric, frame interval (ground truth data is in frames)
#' @param tunit string, temporal units
#' @param ... pass additional parameters to modify the defaults (N, short, deltaT, mode, nPop, init, timeRes, breaks, radius)
#'
#' @return multiple pdf reports
#' @export
compareGTDatasets <- function(path = NULL, xyscale = 0.04, xyunit = "um", tscale = 0.06, tunit = "s", ...) {

  condition <- value <- dataid <- cumulative_distance <- track_duration <- NULL

  if(is.null(path)) {
    # there is no cross-platform way to safely choose directory
    cat("Please organise your ground truth csv files in a folder and use the path as an argument\r")
    return(-1)
  } else {
    datadir <- path
  }

  # ellipsis processing
  l <- NULL
  l <- list(...)
  l <- processEllipsis(l)

  # loop through condition folders within data folder
  condFolderNames <- list.dirs(path = datadir, recursive = FALSE)
  # break if there were no folders in Data directory
  if(identical(condFolderNames, character(0)) == TRUE) {
    return(-1)
  }

  for(i in 1:length(condFolderNames)) {
    condFolderPath <- condFolderNames[i]
    condFolderName <- basename(condFolderPath)
    allTrackMateFiles <- list.files(condFolderPath, pattern = "*.csv")
    # skip if there were no csv files in this folder
    if(identical(allTrackMateFiles, character(0)) == TRUE) {
      next
    }
    calibrate <- TRUE

    cat(paste0("\n","Processing ",condFolderName,"\n"))
    for(j in 1:length(allTrackMateFiles)) {
      fileName <- allTrackMateFiles[j]
      thisFilePath <- paste0(condFolderPath, "/", fileName)
      # read dataset
      tmObj <- readGTFile(thisFilePath)
      # scale dataset if required
      if(calibrate) {
        tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = xyscale, tscalar = tscale, xyunit = xyunit, tunit = tunit)
      }
      tmDF <- tmObj[[1]]
      calibrationDF <- tmObj[[2]]
      # take the units
      units <- calibrationDF$unit[1:2]
      # we need to combine data frames
      # first add a column to id the data
      thisdataid <- paste0(condFolderName,"_",as.character(j))
      tmDF$dataid <- thisdataid
      # calculate MSD
      msdObj <- calculateMSD(tmDF, N = 3, short = 8)
      msdDF <- msdObj[[1]]
      alphaDF <- msdObj[[2]]
      # same here, combine msd summary and alpha summary
      msdDF$dataid <- thisdataid
      alphaDF$dataid <- thisdataid
      # jump distance calc with deltaT of 1
      deltaT <- 1
      jdObj <- calculateJD(dataList = tmObj, deltaT = l$deltaT, nPop = l$nPop, mode = l$mode, init = l$init, timeRes = l$timeRes, breaks = l$breaks)
      jdDF <- jdObj[[1]]
      jdDF$dataid <- thisdataid
      timeRes <- jdObj[[2]]
      jdObj <- list(jdDF,timeRes)
      # track density with a radius of 1.5 units
      tdDF <- calculateTrackDensity(dataList = tmObj, radius = l$radius)
      tdDF$dataid <- thisdataid
      # fractal dimension
      fdDF <- calculateFD(dataList = tmObj)
      fdDF$dataid <- thisdataid
      # now if it's the first one make the big dataframe, otherwise add df to bigdf
      if(j == 1) {
        bigtm <- tmDF
        bigmsd <- msdDF
        bigalpha <- alphaDF
        bigjd <- jdDF
        bigtd <- tdDF
        bigfd <- fdDF
      } else {
        bigtm <- rbind(bigtm,tmDF)
        bigmsd <- rbind(bigmsd,msdDF)
        bigalpha <- rbind(bigalpha,alphaDF)
        bigjd <- rbind(bigjd,jdDF)
        bigtd <- rbind(bigtd,tdDF)
        bigfd <- rbind(bigfd,fdDF)
      }

      # create the report for this dataset
      fileName <- tools::file_path_sans_ext(basename(thisFilePath))
      both <- makeSummaryReport(tmList = tmObj, msdList = msdObj, jumpList = jdObj, tddf = tdDF, fddf = fdDF,
                                titleStr = condFolderName, subStr = fileName, auto = TRUE, summary = FALSE,
                                msdplot = l$msdplot)
      p <- both[[1]]
      destinationDir <- paste0("Output/Plots/GT/", condFolderName)
      setupOutputPath(destinationDir)
      filePath <- paste0(destinationDir, "/report_",as.character(j),".pdf")
      ggsave(filePath, plot = p, width = 25, height = 19, units = "cm")
      # retrieve other data
      df_report <- both[[2]]
      df_report$condition <- condFolderName
      df_report$dataid <- thisdataid
      if(i == 1 & j == 1) {
        megareport <- df_report
      } else if(!exists("megareport")) {
        megareport <- df_report
      } else {
        megareport <- rbind(megareport,df_report)
      }
    }
    bigtmObj <- list(bigtm,calibrationDF)
    bigmsdObj <- list(bigmsd,bigalpha)
    bigjdObj <- list(bigjd,timeRes)
    # now we have our combined dataset we can make a summary
    # note we use the timeRes of the final dataset; so it is suitable for only when all files have the same calibration
    summaryObj <- makeSummaryReport(tmList = bigtmObj, msdList = bigmsdObj, jumpList = bigjdObj, tddf = bigtd, fddf = bigfd,
                                    titleStr = condFolderName, subStr = "Summary", auto = TRUE, summary = TRUE,
                                    msdplot = l$msdplot)
    p <- summaryObj[[1]]
    destinationDir <- paste0("Output/Plots/GT/", condFolderName)
    filePath <- paste0(destinationDir, "/combined.pdf")
    ggsave(filePath, plot = p, width = 25, height = 19, units = "cm")
    # save data as csv
    destinationDir <- paste0("Output/Data/GT/", condFolderName)
    setupOutputPath(destinationDir)
    # save each dataset-level data
    write.csv(bigtm, paste0(destinationDir, "/allTM.csv"), row.names = FALSE)
    write.csv(bigmsd, paste0(destinationDir, "/allMSD.csv"), row.names = FALSE)
    write.csv(bigjd, paste0(destinationDir, "/allJD.csv"), row.names = FALSE)
    write.csv(bigfd, paste0(destinationDir, "/allFD.csv"), row.names = FALSE)
    # mega data frame of msd averages per dataset; alpha values, track density, speed by trace/dataid/condition
    msdSummary <- summaryObj[[2]]
    msdSummary$condition <- condFolderName
    # bigalpha$condition <- condFolderName
    # bigtd$condition <- condFolderName
    bigspeed <- bigtm %>%
      group_by(dataid, trace) %>%
      summarise(cumdist = max(cumulative_distance), cumtime = max(track_duration))
    bigspeed$speed <- bigspeed$cumdist / bigspeed$cumtime
    bigspeed$condition <- condFolderName
    if(i == 1 | !exists("megamsd")) {
      megamsd <- msdSummary
      megaalpha <- bigalpha
      megatd <- bigtd
      megaspeed <- bigspeed
      megafd <- bigfd
    } else {
      megamsd <- rbind(megamsd,msdSummary)
      megaalpha <- rbind(megaalpha,bigalpha)
      megatd <- rbind(megatd,bigtd)
      megaspeed <- rbind(megaspeed,bigspeed)
      megafd <- rbind(megafd,bigfd)
    }
  }

  # save summary data as csv
  destinationDir <- "Output/Data/GT"
  write.csv(megamsd, paste0(destinationDir, "/allMSDCurves.csv"), row.names = FALSE)
  write.csv(megareport, paste0(destinationDir, "/allComparison.csv"), row.names = FALSE)

  # for alpha values, track density, speed by trace/dataid/condition we must combine into one
  megatrace <- Reduce(mergeDataFramesForExport, list(megaalpha, megatd, megaspeed, megafd))
  write.csv(megatrace, paste0(destinationDir, "/allTraceData.csv"), row.names = FALSE)

  # generate the comparison plots and save
  p <- makeComparison(df = megareport, msddf = megamsd, units = units, msdplot = l$msdplot)
  destinationDir <- "Output/Plots/GT"
  filePath <- paste0(destinationDir, "/comparison.pdf")
  ggsave(filePath, plot = p, width = 19, height = 19, units = "cm")
}


#' Read a ground truth csv file
#'
#' Ground truth csv file has four columns TrackID,x,y,frame
#' this function reads the data and formats in a way that is equivalent to the way that TrackMateR reads a TrackMate XML file.
#'
#' @param path string, filepath to ground truth csv file
#'
#' @return list of two data frames
#' @export
readGTFile <- function(path) {

  daten <- read.csv(path, header = TRUE)
  # gt data is TrackID,x,y,frame - trace is required and must be character
  daten$trace <- as.character(daten$TrackID)
  daten$t <- daten$frame

  # construct calibration file
  calibrationDF <- data.frame(value = c(1,1,199,199),
                              unit = c("pixel","frame","widthpixels","heightpixels"))

  # displacement, cumulative distance and duration
  displacement <- sqrt(c(0,diff(daten$x))^2 + c(0,diff(daten$y))^2)
  cumdist <- numeric()
  cumdist[1] <- 0
  dur <- numeric()
  dur[1] <- 0
  startdur <- daten$t[1]

  for (i in 2:nrow(daten)){
    if(daten$trace[i] == daten$trace[i-1]) {
      cumdist[i] <- cumdist[i-1] + displacement[i]
    } else {
      displacement[i] <- 0
      cumdist[i] <- 0
      startdur <- daten$t[i]
    }
    dur[i] <- daten$t[i] - startdur
  }
  daten$displacement <- displacement
  daten$cumulative_distance <- cumdist
  daten$track_duration <- dur

  # it is possible that xy coords lie outside the image(!)
  # we can detect xy coords that are less than 0,0 and then use this information to offset *all* coords by this
  # this is necessary because later code relies on the origin
  minx <- min(daten$x)
  miny <- min(daten$y)
  if(minx < 0) {
    daten$x <- daten$x - minx
  }
  if(miny < 0) {
    daten$y <- daten$y - miny
  }
  # now we need to redefine the size of the "image" because xy coords may lie outside, or they may now lie outside after offsetting
  maxx <- max(daten$x)
  maxy <- max(daten$y)
  if(maxx > calibrationDF[3,1]) {
    calibrationDF[3,1] <- ceiling(maxx)
  }
  if(maxy > calibrationDF[4,1]) {
    calibrationDF[4,1] <- ceiling(maxy)
  }

  # format as tmObj
  dfList <- list(daten,calibrationDF)

  return(dfList)
}
