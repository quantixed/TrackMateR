# Make a plot of cumulative distance over time

Make a plot of cumulative distance over time

## Usage

``` r
plot_tm_cumdistOverTime(
  input,
  summary = FALSE,
  xstr = NULL,
  ystr = NULL,
  alphaLevel = 0.1
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
