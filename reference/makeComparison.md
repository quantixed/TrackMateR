# Make Comparison Plots

A series of ggplots to compare between conditions. Called from
\`compareDatasets()\`, this function generates plots using summary data
of datasets, per condition.

## Usage

``` r
makeComparison(
  df,
  msddf,
  units = c("um", "s"),
  msdplot = "linlin",
  titleStr = "Comparison",
  subStr = NULL
)
```

## Arguments

- df:

  data frame called megareport

- msddf:

  data frame of msd averages per dataset

- units:

  character vector of space and time units for axes

- msdplot:

  keyword used to specify msdplot axes, either "loglog", "linlin"
  (default), "loglin" or "linlog"

- titleStr:

  string used as the title for the report

- subStr:

  string used as the subtitle for the report

## Value

patchwork ggplot

## Details

Note that if one dataset has a long name the names will be wrapped.
Wrapping not graceful. If you have very long names consider renaming
them before running.
