# Calculate Mean Squared Displacement (MSD)

Calculation of the MSD of multiple tracks. There are two methods for
everaging MSD data from multiple tracks: ensemble = for each time lag
average all squared displacements from all tracks time-averaged = find
MSD for each track and then generate the average MSD from these curves
The MSD curves will be identical if all tracks are the same length, and
diverge if not. Standard deviation will be large for ensemble and
smaller for time-averaged data. Input is a data frame of tracks imported
using readTrackMateXML()

## Usage

``` r
calculateMSD(df, method = "timeaveraged", N = 4, short = 0)
```

## Arguments

- df:

  data frame must include at a minimum - trace (track ID), x, y and t
  (in real coords)

- method:

  string. Either "ensemble" or "timeaveraged" (default)

- N:

  numeric variable for MSD. dt should be up to 1/N of number of data
  points (4 recommended)

- short:

  numeric variable for the shortest number of points we will analyse.
  Note, this uses the number of frames from start, not number of points
  in track, i.e. a track with \<short points and many gaps will remain

## Value

list of three data frames

## Examples

``` r
xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
tmObj <- readTrackMateXML(XMLpath = xmlPath)
#> Units are:  1 pixel and 0.07002736 s 
#> Spatial units are in pixels - consider transforming to real units
#> Extracting spot data...
#> Matching track data...
#> Calculating distances...
tmObj <- correctTrackMateData(tmObj, xyscalar = 0.04)
#> Correcting XY scale. 
tmDF <-  tmObj[[1]]
msdObj <- calculateMSD(df = tmDF, method = "ensemble", N = 3, short = 8)
```
