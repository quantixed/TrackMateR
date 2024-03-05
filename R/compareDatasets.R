#' Compare datasets
#'
#' Requires TrackMate XML files to be organised into subfolders named according to the condition.
#' If these condition folders are in `Data/` within the working directory, the code will run automatically.
#' Since there is no easy cross-platform way for the user to interactively pick a directory, the organisation of files in `Data/` is a requirement.
#' The `Data/` folder can be elsewhere on your computer, just change the `wd` prior to running `compareDatasets()` and the routine will run.
#' The code will process all the datasets individually, compile them according to condition and compare across conditions.
#' Outputs are saved to `Output/Plots/` in the working directory.
#'
#' If TrackMate XML files require recalibration, this is possible by placing a csv file into each subfolder.
#' All xml files in that folder whose calibration does not match the calibration csv file will be altered.
#' Ideally, all conditions should have the same scaling, and within a condition they should be similar.
#' The code will run if this is not the case, but beware that these discrepancies are not detected.
#' For example, comparing two datasets in um/s with one in mm/min.
#'
#' @param ... pass additional parameters to modify the defaults (N, short, deltaT, mode, nPop, init, timeRes, breaks, radius)
#'
#' @return multiple pdf reports
#' @export
compareDatasets <- function(...) {

  condition <- value <- dataid <- cumulative_distance <- track_duration <- mean_intensity <- NULL

  if(!dir.exists("Data")) {
    # there is no cross-platform way to safely choose directory
    cat("Please organise your XML files in a folder called Data in the working directory\r")
    return(-1)
  } else {
    datadir <- "Data"
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
    allTrackMateFiles <- list.files(condFolderPath, pattern = "*.xml")
    # skip if there were no XML files in this folder
    if(identical(allTrackMateFiles, character(0)) == TRUE) {
      next
    }
    # check to see if a calibration file is present
    calibrationFiles <- list.files(condFolderPath, pattern = "*.csv")
    calibrationFile <- paste0(condFolderPath,"/",calibrationFiles[1])
    if(identical(calibrationFiles, character(0)) == TRUE) {
      calibrationXY <- 1
      calibrationT <- 1
      calibrate <- FALSE
    } else {
      # load calibration file and store calibrations
      calibDF <- read.csv(calibrationFile)
      calibrate <- TRUE
    }
    cat(paste0("\n","Processing ",condFolderName,"\n"))
    for(j in 1:length(allTrackMateFiles)) {
      fileName <- allTrackMateFiles[j]
      thisFilePath <- paste0(condFolderPath, "/", fileName)
      # read dataset
      tmObj <- readTrackMateXML(XMLpath = thisFilePath)
      # scale dataset if required
      if(calibrate) {
        calibrationDF <- tmObj[[2]]
        # scalar for conversion is new / old (units not relevant)
        calibrationXY <- calibDF[1,1] / calibrationDF[1,1]
        calibrationT <- calibDF[2,1] / calibrationDF[2,1]
        # 0 in calibDF indicates no scaling is to be done
        calibrationXY <- ifelse(calibrationXY == 0, 1, calibrationXY)
        calibrationT <- ifelse(calibrationT == 0, 1, calibrationT)
        # ignore an error of 2.5%
        calibrationXY <- ifelse(calibrationXY < 1.025 & calibrationXY > 0.975, 1, calibrationXY)
        calibrationT <- ifelse(calibrationT < 1.025 & calibrationT > 0.975, 1, calibrationT)
        if(calibrationXY != 1 & calibrationT != 1) {
          tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = calibrationXY, tscalar = calibrationT, xyunit = calibDF[1,2], tunit = calibDF[2,2])
        } else if(calibrationXY != 1 & calibrationT == 1) {
          tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = calibrationXY, xyunit = calibDF[1,2])
        } else if(calibrationXY == 1 & calibrationT != 1) {
          tmObj <- correctTrackMateData(dataList = tmObj, tscalar = calibrationT, tunit = calibDF[2,2])
        } else {
          # the final possibility is nothing needs scaling but units need to change.
          # do not test if they are the same just flush the units with the values in the csv file
          tmObj <- correctTrackMateData(dataList = tmObj, xyunit = calibDF[1,2], tunit = calibDF[2,2])
        }
      }
      # we can filter here if required - for example only analyse tracks of certain length
      tmDF <- tmObj[[1]]
      calibrationDF <- tmObj[[2]]
      # if the data is not rich enough for a summary we will skip it
      if((calibrationDF[5,1] < 3 & calibrationDF[6,1] < 10)) {
        cat(paste0("Skipping ",fileName," as it has less than 3 tracks and longest track has less than 10 frames\n"))
        next
      }
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
      deeDF <- msdObj[[3]]
      # same here, combine msd summary and alpha summary
      msdDF$dataid <- thisdataid
      alphaDF$dataid <- thisdataid
      deeDF$dataid <- thisdataid
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
        bigdee <- deeDF
        bigjd <- jdDF
        bigtd <- tdDF
        bigfd <- fdDF
      } else {
        bigtm <- rbind(bigtm,tmDF)
        bigmsd <- rbind(bigmsd,msdDF)
        bigalpha <- rbind(bigalpha,alphaDF)
        bigdee <- rbind(bigdee,deeDF)
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
      destinationDir <- paste0("Output/Plots/", condFolderName)
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
    bigmsdObj <- list(bigmsd,bigalpha,bigdee)
    bigjdObj <- list(bigjd,timeRes)
    # now we have our combined dataset we can make a summary
    # note we use the timeRes of the final dataset; so it is suitable for only when all files have the same calibration
    summaryObj <- makeSummaryReport(tmList = bigtmObj, msdList = bigmsdObj, jumpList = bigjdObj, tddf = bigtd, fddf = bigfd,
                                    titleStr = condFolderName, subStr = "Summary", auto = TRUE, summary = TRUE,
                                    msdplot = l$msdplot)
    p <- summaryObj[[1]]
    destinationDir <- paste0("Output/Plots/", condFolderName)
    filePath <- paste0(destinationDir, "/combined.pdf")
    ggsave(filePath, plot = p, width = 25, height = 19, units = "cm")
    # save data as csv
    destinationDir <- paste0("Output/Data/", condFolderName)
    setupOutputPath(destinationDir)
    # save each dataset-level data
    write.csv(bigtm, paste0(destinationDir, "/allTM.csv"), row.names = FALSE)
    write.csv(bigmsd, paste0(destinationDir, "/allMSD.csv"), row.names = FALSE)
    write.csv(bigjd, paste0(destinationDir, "/allJD.csv"), row.names = FALSE)
    write.csv(bigfd, paste0(destinationDir, "/allFD.csv"), row.names = FALSE)
    # mega data frame of msd averages per dataset; alpha values, track density, speed/duration/distance, intensity by trace/dataid/condition
    msdSummary <- summaryObj[[2]]
    msdSummary$condition <- condFolderName
    bigspeed <- bigtm %>%
      group_by(dataid, trace) %>%
      summarise(cumdist = max(cumulative_distance), cumtime = max(track_duration), intensity = max(mean_intensity))
    bigspeed$speed <- bigspeed$cumdist / bigspeed$cumtime
    bigspeed$condition <- condFolderName
    if(i == 1 | !exists("megamsd")) {
      megamsd <- msdSummary
      megaalpha <- bigalpha
      megadee <- bigdee
      megatd <- bigtd
      megaspeed <- bigspeed
      megafd <- bigfd
    } else {
      megamsd <- rbind(megamsd,msdSummary)
      megaalpha <- rbind(megaalpha,bigalpha)
      megadee <- rbind(megadee,bigdee)
      megatd <- rbind(megatd,bigtd)
      megaspeed <- rbind(megaspeed,bigspeed)
      megafd <- rbind(megafd,bigfd)
    }
  }

  # save summary data as csv
  destinationDir <- "Output/Data"
  write.csv(megamsd, paste0(destinationDir, "/allMSDCurves.csv"), row.names = FALSE)
  write.csv(megareport, paste0(destinationDir, "/allComparison.csv"), row.names = FALSE)

  # for alpha values, estimator of D, track density, peed/duration/distance, intensity by trace/dataid/condition we must combine into one
  # set the name of dee to estdee
  names(megadee)[names(megadee) == "dee"] <- "estdee"
  megatrace <- Reduce(mergeDataFramesForExport, list(megaalpha, megadee, megatd, megaspeed, megafd))
  write.csv(megatrace, paste0(destinationDir, "/allTraceData.csv"), row.names = FALSE)

  # generate the comparison plots and save
  p <- makeComparison(df = megareport, msddf = megamsd, units = units, msdplot = l$msdplot)
  destinationDir <- "Output/Plots/"
  filePath <- paste0(destinationDir, "/comparison.pdf")
  ggsave(filePath, plot = p, width = 25, height = 19, units = "cm")
}
