# Calculate density of tracks

Calculate for each track, using its starting frame, what is the relative
density of tracks. We use a search radius to find how many tracks in the
starting frame are neighbours of the track. The number of neighbours is
normalised to the search circle, so that a track in the corner of the
image with 2 neighbours has a density of 8. Code for calculating search
area (intersection between search circle and the image border) is taken
from
https://petrelharp.github.io/circle_rectangle_intersection/circle_rectangle_intersection.html

## Usage

``` r
calculateTrackDensity(dataList, radius = 1)
```

## Arguments

- dataList:

  list of a data frame (must include at a minimum - trace (track ID), x,
  y and frame (in real coords)) and a calibration data frame

- radius:

  numeric variable for search radius (in spatial units of the data)

## Value

data frame

## Examples

``` r
xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
tmObj <- readTrackMateXML(XMLpath = xmlPath)
#> Units are:  1 pixel and 0.07002736 s 
#> Spatial units are in pixels - consider transforming to real units
#> Matching track data...
#> Calculating distances...
tmObj <- correctTrackMateData(dataList = tmObj, xyscalar = 0.04)
#> Correcting XY scale. 
tdDF <- calculateTrackDensity(dataList = tmObj, radius = 2)
```
