# TrackMateR

Analysis of TrackMate XML outputs in R.

[TrackMate](https://imagej.net/plugins/trackmate/) is a single-particle
tracking plugin for ImageJ/Fiji. The standard output from a tracking
session is in TrackMate XML format.

The goal of this R package is to import all of the data associated with
the final filtered tracks in TrackMate for further analysis and
visualization in R.

## Installation

Once you have installed [R](https://cran.rstudio.com) and [RStudio
Desktop](https://www.rstudio.com/products/rstudio/download/), you can
install TrackMateR using devtools

``` r
# install.packages("devtools")
devtools::install_github("quantixed/TrackMateR")
```

## An Example

A basic example is to load one TrackMate XML file, calibrate it (if
needed) and analyse it.

``` r
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.5.2
library(TrackMateR)
# an example file is provided, otherwise use file.choose()
xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
# read the TrackMate XML file into R using
tmObj <- readTrackMateXML(XMLpath = xmlPath)
#> Units are:  1 pixel and 0.07002736 s 
#> Spatial units are in pixels - consider transforming to real units
#> Matching track data...
#> Calculating distances...
# Pixel size is actually 0.04 um and original data was 1 pixel, xyscalar = 0.04
tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = 0.04, xyunit = "um")
#> Correcting XY scale.
# generate a report
reportDataset(tmObj)
```

![](reference/figures/README-example-1.png)

TrackMateR can generate several different types of plot individually
using commands or it can make them all automatically and create a report
for you.

- For details of how to make individual plots and/or tweak the default
  parameters, see
  [`vignette("TrackMateR")`](https://quantixed.github.io/TrackMateR/articles/TrackMateR.md)
- To see how to compare different datasets, see
  [`vignette("comparison")`](https://quantixed.github.io/TrackMateR/articles/comparison.md)
- In order to rescale or recalibrate TrackMate data, see
  [`vignette("recalibration")`](https://quantixed.github.io/TrackMateR/articles/recalibration.md)

## Credits

- TrackMateR builds on initial work by Julien Godet on
  [trackR](https://github.com/jgodet/trackR).
- Méghane Sittewelle provided example TrackMate data and helped with
  testing TrackMateR.
- The Fiji plug-in
  [TrackMate](https://github.com/trackmate-sc/TrackMate) is developed by
  Jean-Yves Tinevez.

## Limitations, future development

- TrackMateR is currently written for 2D data. 3D data is read but
  analysis is currently on the first two dimensions.
- Addition of a dry run option to quickly report what data would be
  analysed by
  [`compareDatasets()`](https://quantixed.github.io/TrackMateR/reference/compareDatasets.md).
- Addition of number of tracks per frame and more advanced calculation
  of the density of tracks.
