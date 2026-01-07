# Read a ground truth csv file

Ground truth csv file has four columns TrackID,x,y,frame this function
reads the data and formats in a way that is equivalent to the way that
TrackMateR reads a TrackMate XML file.

## Usage

``` r
readGTFile(path)
```

## Arguments

- path:

  string, filepath to ground truth csv file

## Value

list of two data frames
