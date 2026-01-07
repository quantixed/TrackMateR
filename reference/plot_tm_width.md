# Make a histogram of the largest point-to-point distance in each track

Make a histogram of the largest point-to-point distance in each track

## Usage

``` r
plot_tm_width(df, units = c("um", "s"), auto = FALSE)
```

## Arguments

- df:

  data frame of fractal dimension data

- units:

  character vector of space and time units (default is um and s)

- auto:

  boolean to switch between returning a ggplot and a list of ggplot and
  variable

## Value

ggplot or list of ggplot and variable
