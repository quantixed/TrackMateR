
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TrackMateR <a href='https://quantixed.github.io/TrackMateR/'><img src='man/figures/logo.png' align="right" height="131.5" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/quantixed/TrackMateR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/quantixed/TrackMateR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Analysis of TrackMate XML outputs in R.

[TrackMate](https://imagej.net/plugins/trackmate/) is a single-particle
tracking plugin for ImageJ/Fiji. The standard output from a tracking
session is in TrackMate XML format.

The goal of this R package is to import all of the data associated with
the final filtered tracks in TrackMate for further analysis and
visualization in R.

**This package is under development and should not be considered stable
until its first release.**

## Installation

Once you have installed [R](https://cran.rstudio.com) and [RStudio
Desktop](https://www.rstudio.com/products/rstudio/download/), you can
install TrackMateR using devtools

``` r
# install.packages("devtools")
devtools::install_github("quantixed/TrackMateR")
```

## An Example

A basic example is to load one TrackMate XML file and analyse it.

``` r
library(ggplot2)
library(TrackMateR)
# an example file is provided, otherwise use file.choose()
xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
# read the TrackMate XML file into R using
tmObj <- readTrackMateXML(XMLpath = xmlPath)
#> Units are:  1 pixel and 0.07002736 s 
#> Spatial units are in pixels - consider transforming to real units
#> Collecting spot data...
#> Matching track data...
#> Calculating distances...
# have a look at the tracks
plot_tm_allTracks(tmObj)
```

![](man/figures/README-example-1.png)<!-- -->

Let’s go a bit further.

TrackMateR can generate several different types of plot individually
using commands or it can make them all automatically and create a report
for you.

``` r
# perhaps the data we loaded in was not scaled properly
# Pixel size is 0.04 um and original data was 1 pixel, xyscalar = 0.04
tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = 0.04, xyunit = "um")
#> Correcting XY scale.
# we can get a data frame of the correct TrackMate data 
tmDF <- tmObj[[1]]
# and a data frame of calibration information
calibrationDF <- tmObj[[2]]
# we can calculate mean squared displacement
msdObj <- calculateMSD(df = tmDF, method = "ensemble", N = 3, short = 8)
# and jump distance
jdObj <- calculateJD(dataList = tmObj, deltaT = 1)
# and look at the density of tracks
tdDF <- calculateTrackDensity(dataList = tmObj, radius = 1.5)
# and calculate fractal dimension
fdDF <- calculateFD(dataList = tmObj)
# if we extract the name of the file
fileName <- tools::file_path_sans_ext(basename(xmlPath))
# we can send all these things to makeSummaryReport() to get a nice report of our dataset
makeSummaryReport(tmList = tmObj, msdList = msdObj, jumpList = jdObj, tddf = tdDF, fddf = fdDF,
titleStr = "Report", subStr = fileName, auto = FALSE)
```

![](man/figures/README-unnamed-chunk-2-1.png)<!-- -->

## Credits

TrackMateR builds on initial work by Julien Godet on
[trackR](https://github.com/jgodet/trackR). Méghane Sittewelle provided
example TrackMate data and helped with testing TrackMateR.