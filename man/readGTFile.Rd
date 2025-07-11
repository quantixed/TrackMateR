% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gt_functions.R
\name{readGTFile}
\alias{readGTFile}
\title{Read a ground truth csv file}
\usage{
readGTFile(path)
}
\arguments{
\item{path}{string, filepath to ground truth csv file}
}
\value{
list of two data frames
}
\description{
Ground truth csv file has four columns TrackID,x,y,frame
this function reads the data and formats in a way that is equivalent to the way that TrackMateR reads a TrackMate XML file.
}
