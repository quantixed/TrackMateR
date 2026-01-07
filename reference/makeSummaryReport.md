# Make Summary Report

Generate several plots to visualise TrackMate data and generate a
report. The use of um is because ggsave does not currently save unicode
to PDF reliably. The code is switchable to accommodate making a "report"
(one dataset) or a "summary" (several related datasets combined)

## Usage

``` r
makeSummaryReport(
  tmList,
  msdList,
  jumpList,
  tddf,
  fddf,
  titleStr = "",
  subStr = "",
  auto = FALSE,
  summary = FALSE,
  ...
)
```

## Arguments

- tmList:

  list of trackmate data and calibration

- msdList:

  MSD summary and alpha list = output from calculateMSD()

- jumpList:

  list of a data frame of jump data and a variable to be passed to
  timeRes

- tddf:

  data frame of track density data

- fddf:

  data frame of fractal dimension data

- titleStr:

  string used as the title for the report

- subStr:

  string used as the subtitle for the report

- auto:

  boolean which selects for returning the patchwork report (FALSE) or a
  list of the patchwork report and a data frame of summary (TRUE)

- summary:

  boolean which selects for report (FALSE) or summary (TRUE)

- ...:

  additional arguments passed to methods

## Value

patchwork ggplot or a list of patchwork ggplot and data frame of summary
data

## Examples

``` r
xmlPath <- system.file("extdata", "ExampleTrackMateData.xml", package="TrackMateR")
tmObj <- readTrackMateXML(XMLpath = xmlPath)
#> Units are:  1 pixel and 0.07002736 s 
#> Spatial units are in pixels - consider transforming to real units
#> Matching track data...
#> Calculating distances...
tmObj <- correctTrackMateData(tmObj, xyscalar = 0.04)
#> Correcting XY scale. 
tmDF <- tmObj[[1]]
calibrationDF <- tmObj[[2]]
msdObj <- calculateMSD(df = tmDF, method = "ensemble", N = 3, short = 8)
jdObj <- calculateJD(dataList = tmObj, deltaT = 1)
tdDF <- calculateTrackDensity(dataList = tmObj, radius = 1.5)
fdDF <- calculateFD(dataList = tmObj)
fileName <- tools::file_path_sans_ext(basename(xmlPath))
reportObj <- makeSummaryReport(tmList = tmObj, msdList = msdObj, jumpList = jdObj,
tddf = tdDF, fddf = fdDF, titleStr = "Report", subStr = fileName, auto = TRUE)
```
