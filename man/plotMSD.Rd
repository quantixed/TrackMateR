% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMSD.R
\name{plotMSD}
\alias{plotMSD}
\title{Make a plot of MSD data}
\usage{
plotMSD(df, units = "s", bars = TRUE, xlog = FALSE, ylog = FALSE)
}
\arguments{
\item{df}{MSD summary = output from calculateMSD()}

\item{units}{string to describe time units (default is s, seconds)}

\item{bars}{boolean to request error bars (1 x SD)}

\item{xlog}{boolean to request log10 x axis}

\item{ylog}{boolean to request log10 y axis}
}
\value{
S3 ggplot
}
\description{
Generate a plot of MSD over a series of increasing time lags.
Input is the output from CalculateMSD(), so the plot will display the ensemble or time-averaged MSD (whatever was requested)
A fit to the first four points is displayed to evaluate alpha. Diffusion coefficient from this fit is displayed top-left.
}
\examples{
xmlPath <- "~/Desktop/FakeTracks.xml"
data <- readTrackMateXML(XMLpath = xmlPath)
data <- correctTrackMateData(data, xy = 0.04)
msdDF <- calculateMSD(data, method = "ensemble", N = 3, short = 8)
plotMSD(msdDF, bars = FALSE)
}
