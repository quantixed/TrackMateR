# Convert TrackMate XML file to CSV format

This function reads a TrackMate XML file and save it as a csv file for
use in other contexts.

There are options to: - filter out traces with less than n points - save
in a location other than Output/Data - select columns for export - x and
y can be scaled (default) or converted to pixels

The typical scenario is that the XML file is in \`input_dir\` and the
converted csv will be saved to \`output_dir\`. Any subfolders below
\`input_dir\` will be recreated in \`output_dir\`. If a full path is
supplied as the \`path\`, the \`input_dir\` should be adjusted to be the
directory containing the XML file or containing the subdirectory
containing the XML file.

Note: the calibration of the input TrackMate XML file must be correct.

## Usage

``` r
convert_trackmate_xml_to_csv(
  path = NULL,
  min_points = 10,
  input_dir = "Data/",
  output_dir = "Output/Data/",
  slim = TRUE,
  columns = c("trace", "frame", "y", "x"),
  pixels = FALSE,
  no_save = FALSE
)
```

## Arguments

- path:

  path to the TrackMate XML file

- min_points:

  minimum number of points in a trace to be included in the output csv
  file (default: 10). Set to 0 to include all traces.

- input_dir:

  directory where the input XML files are located (default: "Data/")

- output_dir:

  directory where the output csv file will be saved (default:
  "Output/Data/").

- slim:

  logical, if TRUE, only the minimum data required to process the tracks
  is read (default: TRUE)

- columns:

  columns to be included in the output csv file (default: c(" trace",
  "frame", "y", "x")). Set to NULL to include all columns from the
  TrackMate XML file.

- pixels:

  logical, if TRUE, x and y will be converted to pixels using the
  calibration information from the TrackMate XML file (default: FALSE)

- no_save:

  logical, if TRUE, the output csv file will not be saved and a data
  frame returned to the user (default: FALSE)

## Value

saves a csv file

## Examples

``` r
if (FALSE) { # \dontrun{
files <- list.files(path = "Data", pattern = "\\.xml$", full.names = TRUE)
lapply(files, convert_trackmate_xml_to_csv)
} # }
```
