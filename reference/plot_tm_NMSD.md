# Plot several (n) MSD curves

Generate a plot of several MSD curves together with a summary curve.
This function is used to compile multiple datasets from the same
condition.

## Usage

``` r
plot_tm_NMSD(df, xlog = FALSE, ylog = FALSE, auto = FALSE)
```

## Arguments

- df:

  dataframe of MSD summary data from multiple datasets (labelled by
  dataid)

- xlog:

  boolean to request log10 x axis

- ylog:

  boolean to request log10 y axis

- auto:

  boolean to request plot only, TRUE gives plot and summary dataframe as
  a list

## Value

ggplot or list of ggplot and dataframe of summary data
