#' Read TrackMate XML output files.
#'
#' Produces a data frame of all spots from filtered tracks, ordered by track number.
#' A warning is generated if the scaling is in pixels rather than real units.
#'
#' @param XMLpath path to the xml file
#' @return list of two data frames
#' @examples
#' xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
#' tmObj <- readTrackMateXML(XMLpath = xmlPath)
#' # get the track data in a data frame
#' tmDF <-  tmObj[[1]]
#' # get the calibration data in a data frame
#' calibrationDF <- tmObj[[2]]
#' @export

readTrackMateXML<- function(XMLpath){

  # get necessary XMLNodeSet
  e <- xmlParse(XMLpath)
  track <- getNodeSet(e, "//Track")
  if(length(track) == 0) {
    stop("File is not a TrackMate XML file")
  }
  filtered <- getNodeSet(e, "//TrackID")
  subdoc <- getNodeSet(e,"//AllSpots//SpotsInFrame//Spot")
  sublist <- getNodeSet(e,"//FeatureDeclarations//SpotFeatures//Feature/@feature")
  attrName <- c("name",unlist(sublist))
  # what are the units?
  attrList <- c("spatialunits","timeunits")
  unitVec <- sapply(attrList, function(x) xpathSApply(e, "//Model", xmlGetAttr, x))
  unitVec <- c(unitVec,"widthpixels","heightpixels","ntraces","maxframes")
  attrList <- c("pixelwidth","timeinterval","width","height")
  valueVec <- sapply(attrList, function(x) xpathSApply(e, "//ImageData", xmlGetAttr, x))
  calibrationDF <- data.frame(value = c(as.numeric(valueVec),0,0),
                              unit = unitVec)
  # convert the width and height of the image from pixels (0-based so minus 1 from wdth/height) to whatever units the TrackMate file uses
  calibrationDF[3:4,1] <- (calibrationDF[3:4,1] - 1) * calibrationDF[1,1]
  # use s for seconds
  calibrationDF[2,2] <- ifelse(calibrationDF[2,2] == "sec", "s", calibrationDF[2,2])
  # readout what the units were and warn if spatial units are pixels
  cat("Units are: ",calibrationDF[1,1],calibrationDF[1,2],"and",calibrationDF[2,1],calibrationDF[2,2],"\n")
  if(unitVec[1] == "pixel") {
    cat("Spatial units are in pixels - consider transforming to real units\n")
  }

  # multicore processing
  numCores <- parallelly::availableCores()

  if (.Platform[["OS.type"]] == "windows") {
    ## PSOCK-based parallel processing
    cl <- parallel::makeCluster(numCores)
    on.exit(parallel::stopCluster(cl))
    registerDoParallel(cl = cl)
  } else {
    ## Forked parallel processing
    registerDoParallel(cores = numCores)
  }

  # perform parallel read
  # test if we are running on windows
  if (.Platform[["OS.type"]] == "windows") {
    cat("Collecting spot data...\n")
    dtf <- as.data.frame(foreach(i = 1:length(attrName), .packages = c("foreach","XML"), .combine = cbind) %do% {
      sapply(subdoc, xmlGetAttr, attrName[i])
    })
  } else {
    cat(paste0("Collecting spot data. Using ",numCores," cores\n"))
    dtf <- as.data.frame(foreach(i = 1:length(attrName), .combine = cbind) %dopar% {
      sapply(subdoc, xmlGetAttr, attrName[i])
    })
  }

  for (i in 2:length(attrName)){
    suppressWarnings(dtf[,i] <- as.numeric(as.character(dtf[,i])))
  }
  # more R-like headers
  headerNames <- tolower(attrName)
  # remove _CH1 or whatever from headers
  headerNames <- gsub("\\wCH\\d$","",headerNames,ignore.case = T)
  # change x y z t
  headerNames <- gsub("^position\\w","",headerNames,ignore.case = T)
  names(dtf) <- headerNames

  cat("Matching track data...\n")

  # trace is an alternative name for track
  IDtrace <- data.frame(name = NA, trace = NA, displacement = NA, speed = NA)

  for (i in seq(along = track)){
    subDoc = xmlDoc(track[[i]])
    IDvec <- unique(c(unlist(xpathApply(subDoc, "//Edge", xmlGetAttr, "SPOT_SOURCE_ID")),
                      unlist(xpathApply(subDoc, "//Edge", xmlGetAttr, "SPOT_TARGET_ID"))))
    traceVec <- rep(sapply(track[i], function(el){xmlGetAttr(el, "TRACK_ID")}), length(IDvec))
    traceDF <- data.frame(name = paste0("ID",IDvec), trace = traceVec)
    # now retrieve displacement and speed for target spots in track
    targetVec <- unlist(xpathApply(subDoc, "//Edge", xmlGetAttr, "SPOT_TARGET_ID"))
    dispVec <- unlist(xpathApply(subDoc, "//Edge", xmlGetAttr, "DISPLACEMENT"))
    speedVec <- unlist(xpathApply(subDoc, "//Edge", xmlGetAttr, "SPEED"))
    dataDF <- data.frame(name = paste0("ID",targetVec), displacement = as.numeric(dispVec), speed = as.numeric(speedVec))
    # left join (will give NA for the first spot)
    allDF <- merge(x = traceDF, y = dataDF, by = "name", all.x = TRUE)
    allDF$displacement <- ifelse(is.na(allDF$displacement), 0, allDF$displacement)
    allDF$speed <- ifelse(is.na(allDF$speed), 0, allDF$speed)
    # grow the dataframe
    IDtrace <- rbind(IDtrace, allDF)
  }

  # merge track information with spot data
  daten <- merge(IDtrace, dtf, by="name")
  # now we subset for filtered tracks
  FTvec <- sapply(filtered, xmlGetAttr, "TRACK_ID")
  daten <- subset(daten, trace %in% FTvec)
  # sort final data frame
  daten <- daten[order(daten$trace, daten$t),]

  cat("Calculating distances...\n")

  # cumulative distance and duration
  cumdist <- numeric()
  cumdist[1] <- 0
  dur <- numeric()
  dur[1] <- 0
  startdur <- daten$t[1]

  for (i in 2:nrow(daten)){
    if(daten$trace[i] == daten$trace[i-1]) {
      cumdist[i] <- cumdist[i-1] + daten$displacement[i]
    }else{
      cumdist[i] <- 0
      startdur <- daten$t[i]
    }
    dur[i] <- daten$t[i] - startdur
  }
  daten$cumulative_distance <- cumdist
  daten$track_duration <- dur

  # users have reported tracks/traces with multiple spots per frame which cannot be processed
  # possibly caused by track splitting or merging. Produce warning about this.
  b <- daten %>% group_by(trace,frame) %>%
    reframe(n = n())
  # because it is not possible to have n = 0 in this column (as it is reframe) we can do
  if(sum(b$n) > nrow(b)) {
    cat("Warning: Detected multiple spots per frame for one or more tracks.\nTrackMateR will only process single tracks. Subsequent analysis will likely fail!\n")
  }

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
  # we need to know how many traces we have and how many frames is the longest one for later
  calibrationDF[5,1] <- length(unique(daten$trace))
  calibrationDF[6,1] <- max(daten$track_duration) / calibrationDF[2,1]

  # daten is our dataframe of all data, calibrationDF is the calibration data
  dfList <- list(daten,calibrationDF)

  return(dfList)
}
