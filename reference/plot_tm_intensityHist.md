# Make a histogram of intensities

Function uses the maximum value of the mean_intensity column for all
traces in supplied TrackMate data

## Usage

``` r
plot_tm_intensityHist(
  input,
  summary = FALSE,
  xstr = "Intensity (AU)",
  ystr = "Frequency",
  auto = FALSE
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

- auto:

  boolean to switch between returning a ggplot and a list of ggplot and
  variable

## Value

ggplot or list of ggplot and variable
