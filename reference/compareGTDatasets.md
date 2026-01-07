# Compare ground truth datasets

Requires Ground Truth csv files to be organised into subfolders named
according to the condition. Outputs are saved to \`Output/Plots/GT/\` in
the working directory.

## Usage

``` r
compareGTDatasets(
  path = NULL,
  xyscale = 0.04,
  xyunit = "um",
  tscale = 0.06,
  tunit = "s",
  ...
)
```

## Arguments

- path:

  string, filepath to directory equivalent to Data in compareDatasets()

- xyscale:

  numeric, pixel size (ground truth data is in pixel units)

- xyunit:

  string, spatial units

- tscale:

  numeric, frame interval (ground truth data is in frames)

- tunit:

  string, temporal units

- ...:

  pass additional parameters to modify the defaults (N, short, deltaT,
  mode, nPop, init, timeRes, breaks, radius)

## Value

multiple pdf reports
