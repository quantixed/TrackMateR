% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tm_plots.R
\name{plot_tm_displacementOverTime}
\alias{plot_tm_displacementOverTime}
\title{Make a plot of displacement over time}
\usage{
plot_tm_displacementOverTime(input, summary = FALSE, xstr = NULL, ystr = NULL)
}
\arguments{
\item{input}{either a data frame of TrackMate data or list of TrackMate data frame and calibration data frame}

\item{summary}{boolean to specify if plot is of one dataset or several related datasets}

\item{xstr}{string to label x-axis}

\item{ystr}{string to label y-axis}
}
\value{
ggplot
}
\description{
Make a plot of displacement over time
}
