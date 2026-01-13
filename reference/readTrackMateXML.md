# Read TrackMate XML output files.

Produces a data frame of all spots from filtered tracks, ordered by
track number. A warning is generated if the scaling is in pixels rather
than real units.

## Usage

``` r
readTrackMateXML(XMLpath, slim = FALSE)
```

## Arguments

- XMLpath:

  path to the xml file

- slim:

  if TRUE, only the minimum data required to process the tracks is
  returned

## Value

list of two data frames

## Examples

``` r
xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
tmObj <- readTrackMateXML(XMLpath = xmlPath)
#> Units are:  1 pixel and 0.07002736 s 
#> Spatial units are in pixels - consider transforming to real units
#> Extracting spot data...
#> Matching track data...
#> Calculating distances...
# get the track data in a data frame
tmDF <-  tmObj[[1]]
# get the calibration data in a data frame
calibrationDF <- tmObj[[2]]
```
