# Make a plot of all tracks

Make a plot of all tracks

## Usage

``` r
plot_tm_allTracks(
  input,
  summary = FALSE,
  xstr = "",
  ystr = "",
  alphaLevel = 0.5
)
```

## Arguments

- input:

  either a data frame of TrackMate data or list of TrackMate data frame
  and calibration data frame

- summary:

  boolean to specify if plot is of one dataset or several related
  datasets

- xstr:

  string to label x-axis

- ystr:

  string to label y-axis

- alphaLevel:

  numeric variable to set alpha for the plot

## Value

ggplot
