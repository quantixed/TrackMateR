# Correct distance and time of imported TrackMate data.

If the TrackMate data is in pixels and/or frames, the data frame can be
converted with this function.

## Usage

``` r
correctTrackMateData(
  dataList,
  xyscalar = 1,
  tscalar = 1,
  xyunit = NULL,
  tunit = NULL
)
```

## Arguments

- dataList:

  a list of a data frame (of track data) and a calibration data frame
  (from TrackMate XML)

- xyscalar:

  numeric multiplier to correct pixel size of original movie. Assumes
  isotropic scaling, i.e. pixel height = pixel width

- tscalar:

  numeric multiplier to correct frame interval of original movie. Frame
  interval of tracked data.

- xyunit:

  string to describe spatial units

- tunit:

  string to describe temporal unit

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
# in the case where pixel size is 0.03 um and original data is 1 pixel, xyscalar = 0.03
tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = 0.03)
#> Correcting XY scale. 
```
