# Calculate fractal dimension (FD)

Calculate for each track, its fractal dimension Katz & George (1985)
define fractal dimension (D) as log(n) / (log(n) + log(d/L)) where n is
the number of steps, d is the longest of all possible point-to-point
distances and L is the cumulative length of the track. D is ~1 for
directed trajectories, ~2 for confined and ~3 for subdiffusion Here we
calculate this and store D (called fd) and d (called wide) for return.
Note that this method does not take into account gaps in the track. For
a track with many gaps, n will be lowered.

## Usage

``` r
calculateFD(dataList)
```

## Arguments

- dataList:

  list of a data frame (must include at a minimum - trace (track ID), x,
  y and frame (in real coords)) and a calibration data frame

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
fdDF <- calculateFD(dataList = tmObj)
```
