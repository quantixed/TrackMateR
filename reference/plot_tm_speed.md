# Make a histogram of average speed

Make a histogram of average speed

## Usage

``` r
plot_tm_speed(input, summary = FALSE, xstr = NULL, ystr = NULL, auto = FALSE)
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

- auto:

  boolean to switch between returning a ggplot and a list of ggplot and
  variable

## Value

ggplot or list of ggplot and variable
