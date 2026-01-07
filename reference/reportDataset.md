# Report Dataset

This function automates the process of making a report. If more control
is needed over the default parameters: a user can generate MSD, jump
distance and track density objects from their imported TrackMate data,
and make a report manually.

## Usage

``` r
reportDataset(tmList, ...)
```

## Arguments

- tmList:

  list of TrackMate data and calibration

- ...:

  pass additional parameters to modify the defaults (N, short, deltaT,
  mode, nPop, init, timeRes, breaks, radius)

## Value

patchwork ggplot
