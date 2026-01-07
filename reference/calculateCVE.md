# Calculate CVE (covariance-based estimator)

From Vestergaard et al., Physical Review E (2019) 89, 022726 Input is
two data matrices of X and Y displacements for each trace. Output is
estimator of D and sigma for each trace.

## Usage

``` r
calculateCVE(xMat, yMat, tlist, tstep)
```

## Arguments

- xMat:

  matrix of x displacements, each col is a track, each row is time (will
  contain NAs)

- yMat:

  matrix of y displacements, each col is a track, each row is time (will
  contain NAs)

- tlist:

  list of trace names

- tstep:

  variable. Time step in seconds

## Value

data frame
