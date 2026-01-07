# Calculate Jump Distance (JD)

Calculation of the JD of multiple tracks. Calculation is equivalent to a
single time lag point on the ensemble MSD curve, typically represented
as a histogram Input is a data frame of tracks imported using
readTrackMateXML() The default time step is one frame - which is the
equivalent to the plot generated to show displacement versus time.

## Usage

``` r
calculateJD(
  dataList,
  deltaT = 1,
  nPop = 2,
  mode = "ECDF",
  init = NULL,
  timeRes = 1,
  breaks = 100
)
```

## Arguments

- dataList:

  list of data frame (must include at a minimum - trace (track ID), x, y
  and t (in real coords)) and calibration

- deltaT:

  integer to represent the multiple of frames that are to be analysed

- nPop:

  integer (either 1,2 or 3) number of populations for the jump distance
  fitting

- mode:

  string indicated ECDF (default) or hist (histogram)

- init:

  initialisation parameters for the nls fit for example list(D2 = 200,
  D1 = 0.1) or list(D2 = 0.01, D1=0.1, D3=10, D4=100)

- timeRes:

  time resolution per unit of jump. Frame interval is 0.5 s and jump
  interval is two steps, timeRes = 1.

- breaks:

  number of bins for histogram. With ECDF breaks can be high e.g. 100,
  for mode = "hist" they should be low, perhaps 30.

## Value

a list of data frame of jump distances, NAs removed; and a list of
parameters called jumpParam (jumptime, deltaT, nPop)

## Examples

``` r
xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
tmObj <- readTrackMateXML(XMLpath = xmlPath)
#> Units are:  1 pixel and 0.07002736 s 
#> Spatial units are in pixels - consider transforming to real units
#> Matching track data...
#> Calculating distances...
tmObj <- correctTrackMateData(tmObj, xyscalar = 0.04)
#> Correcting XY scale. 
jdObj <- calculateJD(dataList = tmObj, deltaT = 2)
```
