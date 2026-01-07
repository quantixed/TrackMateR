# Track Density - Find search area

Find the area of the intersection of the circle centered at xy with
radius r and the radius with vertical sides at a and horizontal sides at
b. xy, a, and b must be vectors of length 2, and xy must lie within the
rectangle.

## Usage

``` r
find_td_area(r, xy, a, b)
```

## Arguments

- r:

  radius of search circle

- xy:

  numeric vector (length 2)

- a:

  numeric vector (length 2)

- b:

  numeric vector (length 2)

## Value

numeric variable

## Examples

``` r
find_td_area(r=2, xy=c(4, 4), a=c(0, 8), b=c(0, 5))
#> [1] 10.10963
```
