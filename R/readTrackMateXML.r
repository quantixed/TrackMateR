#' Read TrackMate XML output files.
#'
#' Produces a data frame of all spots from filtered tracks, ordered by track number.
#' A warning is generated if the scaling is in pixels rather than real units.
#'
#' @param XMLpath path to the xml file
#' @return data frame
#' @examples
#' xmlPath <- "~/Desktop/FakeTracks.xml"
#' data <- readTrackMateXML(XMLpath = xmlPath)
#' @export

readTrackMateXML<- function(XMLpath){
  # load required libraries
  # if(!require("xml2")) {
  #   install.packages("xml2")
  #   library("xml2")
  # }
  # if(!require("XML")) {
  #   install.packages("XML")
  #   library("XML")
  # }
  # if(!require("foreach")) {
  #   install.packages("foreach")
  #   library("foreach")
  # }
  # if(!require("doParallel")) {
  #   install.packages("doParallel")
  #   library("doParallel")
  # }

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
  spaceunit <- getNodeSet(e,"//Model/@spatialunits")
  timeunit <- getNodeSet(e,"//Model/@timeunits")
  unitvec <- c(unlist(spaceunit),unlist(timeunit))
  cat("Units are: ",unitvec,"\n")
  if(unitvec[1] == "pixel") {
    cat("WARNING: Consider transforming to real units\n")
  }

  # multicore processing
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    numCores <- 2L
  } else {
    # use all cores in devtools::test()
    numCores <- parallel::detectCores()
  }
  registerDoParallel(numCores)

  cat("Collecting spot data...\n")

  dtf <- as.data.frame(foreach(i = 1:length(attrName), .combine = cbind) %dopar% {
    sapply(subdoc, xmlGetAttr, attrName[i])
  })
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

  return(daten)
}
