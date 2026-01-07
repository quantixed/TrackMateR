# Package index

## Loading and correcting TrackMate data

To start, TrackMate XML data needs to be read into R. It may need
recalibrating if the data were captured with incorrect parameters.

- [`readTrackMateXML()`](https://quantixed.github.io/TrackMateR/reference/readTrackMateXML.md)
  : Read TrackMate XML output files.
- [`correctTrackMateData()`](https://quantixed.github.io/TrackMateR/reference/correctTrackMateData.md)
  : Correct distance and time of imported TrackMate data.

## Automated reports

TrackMateR can be used manually to have fine control over your analysis,
however if you just want to analyse your data with default parameters,
these two functions will process one dataset, or many

- [`reportDataset()`](https://quantixed.github.io/TrackMateR/reference/reportDataset.md)
  : Report Dataset
- [`compareDatasets()`](https://quantixed.github.io/TrackMateR/reference/compareDatasets.md)
  : Compare datasets
- [`compareDatasetValues()`](https://quantixed.github.io/TrackMateR/reference/compareDatasetValues.md)
  : Compare dataset values

## Analysing track dynamics

To analyse track dynamics, there are several functions to calculate
them. These functions are called (with set parameters) in the automated
workflows and all of them must be called to make a report.

- [`calculateAlpha()`](https://quantixed.github.io/TrackMateR/reference/calculateAlpha.md)
  : Calculate alpha (relationship between MSD and normal diffusion)
- [`calculateCVE()`](https://quantixed.github.io/TrackMateR/reference/calculateCVE.md)
  : Calculate CVE (covariance-based estimator)
- [`calculateFD()`](https://quantixed.github.io/TrackMateR/reference/calculateFD.md)
  : Calculate fractal dimension (FD)
- [`calculateJD()`](https://quantixed.github.io/TrackMateR/reference/calculateJD.md)
  : Calculate Jump Distance (JD)
- [`calculateMSD()`](https://quantixed.github.io/TrackMateR/reference/calculateMSD.md)
  : Calculate Mean Squared Displacement (MSD)
- [`calculateTrackDensity()`](https://quantixed.github.io/TrackMateR/reference/calculateTrackDensity.md)
  : Calculate density of tracks
- [`fittingJD()`](https://quantixed.github.io/TrackMateR/reference/fittingJD.md)
  : Fitting jump distance (JD) data

## Generating reports and summaries

Reports and summaries are collections of plots. They can be generated
automatically using the automated workflows, or they can be made
manually after fine-tuning the analysis of track dynamics.

- [`makeSummaryReport()`](https://quantixed.github.io/TrackMateR/reference/makeSummaryReport.md)
  : Make Summary Report
- [`makeComparison()`](https://quantixed.github.io/TrackMateR/reference/makeComparison.md)
  : Make Comparison Plots

## Plotting functions

The autoimated workflows and the reports that can be generated with
TrackMateR use several different plotting functions.

- [`plot_tm_MSD()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_MSD.md)
  : Make a plot of MSD data
- [`plot_tm_NMSD()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_NMSD.md)
  : Plot several (n) MSD curves
- [`plot_tm_allTracks()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_allTracks.md)
  : Make a plot of all tracks
- [`plot_tm_alpha()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_alpha.md)
  : Make a histogram of alpha values
- [`plot_tm_cumdistOverTime()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_cumdistOverTime.md)
  : Make a plot of cumulative distance over time
- [`plot_tm_dee()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_dee.md)
  : Make a histogram of D values
- [`plot_tm_displacementHist()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_displacementHist.md)
  : Make a histogram of displacements
- [`plot_tm_displacementOverTime()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_displacementOverTime.md)
  : Make a plot of displacement over time
- [`plot_tm_durationHist()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_durationHist.md)
  : Make a histogram of track durations
- [`plot_tm_fd()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_fd.md)
  : Make a histogram of fractal dimension
- [`plot_tm_intensityHist()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_intensityHist.md)
  : Make a histogram of intensities
- [`plot_tm_neighbours()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_neighbours.md)
  : Make a histogram of track density (number of neighbours)
- [`plot_tm_speed()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_speed.md)
  : Make a histogram of average speed
- [`plot_tm_width()`](https://quantixed.github.io/TrackMateR/reference/plot_tm_width.md)
  : Make a histogram of the largest point-to-point distance in each
  track

## Analysis of non-TrackMate data

TrackMateR can read and process csv files in a similar way to TrackMate
XML outputs. These are typically “ground truth” data sets that have been
generated synthetically, but could come from some other package.

- [`compareGTDatasets()`](https://quantixed.github.io/TrackMateR/reference/compareGTDatasets.md)
  : Compare ground truth datasets
- [`readGTFile()`](https://quantixed.github.io/TrackMateR/reference/readGTFile.md)
  : Read a ground truth csv file
