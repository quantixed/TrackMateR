# Compare dataset values

This function is related to \`compareDatasets()\`, but instead of
analysing the motions and outputting graphical results, the aim is to
convert the TrackMate XML files into data frames, compile them according
to condition and save them in CSV format. Additionally, summary
statistics are calculated and saved.

Requires TrackMate XML files to be organised into subfolders named
according to the condition. The default location for these condition
folders is in \`Data/\` within the working directory. If they are
somewhere else, use the \`datadir\` argument to supply the path. The
code will process all the datasets individually, compile them according
to condition and compare across conditions. Outputs are saved to
\`Output/Plots/\` in the working directory.

If TrackMate XML files require recalibration, this is possible by
placing a csv file into each subfolder. All xml files in that folder
whose calibration does not match the calibration csv file will be
altered. Ideally, all conditions should have the same scaling, and
within a condition they should be similar. The code will run if this is
not the case, but beware that these discrepancies are not detected. For
example, comparing two datasets in um/s with one in mm/min.

## Usage

``` r
compareDatasetValues(datadir = "Data")
```

## Arguments

- datadir:

  string path to the folder containing condition subfolders with
  TrackMate XML files. Default is "Data"

## Value

csv files saved to \`Output/Data/\`
