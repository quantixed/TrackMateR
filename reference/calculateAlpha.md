# Calculate alpha (relationship between MSD and normal diffusion)

Alpha is the MSD exponent. Normal diffusion is alpha = 1. Subdiffusion
is alpha \< 1 and superdiffusion is alpha \> 1. Input is a data matrix
of msd curves. Output is mean of log2(alpha), one value for each trace.
D, calculated from a fit of the first four data points, is also
outputted. The method excludes the first four data points from the alpha
calculation - may result in NAs if tracks are too short.

## Usage

``` r
calculateAlpha(alphaMat, tstep)
```

## Arguments

- alphaMat:

  matrix of msd curves, each col is a track, each row is time lag (will
  contain NAs)

- tstep:

  variable. Time step in seconds

## Value

data frame
