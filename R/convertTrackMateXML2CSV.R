#' Convert TrackMate XML file to CSV format
#'
#' @description This function reads a TrackMate XML file and save it as a csv
#'   file for use in other contexts.
#'
#'   There are options to:
#' - filter out traces with less than n points
#' - save in a location other than Output/Data
#' - select columns for export
#' - x and y can be scaled (default) or converted to pixels
#'
#'   The typical scenario is that the XML file is in `input_dir` and the
#'   converted csv will be saved to `output_dir`. Any subfolders below
#'   `input_dir` will be recreated in `output_dir`. If a full path is supplied
#'   as the `path`, the `input_dir` should be adjusted to be the directory
#'   containing the XML file or containing the subdirectory containing the XML
#'   file.
#'
#'   Note: the calibration of the input TrackMate XML file must be correct.
#'
#' @param path path to the TrackMate XML file
#' @param min_points minimum number of points in a trace to be included in the
#'   output csv file (default: 10). Set to 0 to include all traces.
#' @param input_dir directory where the input XML files are located (default:
#'   "Data/")
#' @param output_dir directory where the output csv file will be saved (default:
#'   "Output/Data/").
#' @param slim logical, if TRUE, only the minimum data required to process the
#'   tracks is read (default: TRUE)
#' @param columns columns to be included in the output csv file (default: c("
#'   trace", "frame", "y", "x")). Set to NULL to include all columns from the
#'   TrackMate XML file.
#' @param pixels logical, if TRUE, x and y will be converted to pixels using the
#'   calibration information from the TrackMate XML file (default: FALSE)
#' @param no_save logical, if TRUE, the output csv file will not be saved and a
#'   data frame returned to the user (default: FALSE)
#'
#' @returns saves a csv file
#' @export
#'
#' @examples
#' \dontrun{
#' files <- list.files(path = "Data", pattern = "\\.xml$", full.names = TRUE)
#' lapply(files, convert_trackmate_xml_to_csv)
#' }
convert_trackmate_xml_to_csv <- function(path = NULL,
                                         min_points = 10,
                                         input_dir = "Data/",
                                         output_dir = "Output/Data/",
                                         slim = TRUE,
                                         columns = c("trace", "frame", "y", "x"),
                                         pixels = FALSE,
                                         no_save = FALSE) {
  XMLpath <- path
  if (is.null(XMLpath)) {
    stop("Please provide a path to the TrackMate XML file.")
  }
  # determine the subpath of path minus the input_dir
  subpath <- gsub(paste0("^", input_dir), "", XMLpath)
  # if subpath starts with a slash, remove it
  subpath <- gsub("^/", "", subpath)
  # create the output path by combining output_dir with the subpath, but
  # replacing the .xml extension with .csv
  output_csv <- file.path(output_dir, subpath)
  output_csv <- sub("\\.xml$", ".csv", output_csv)

  # create output directory if it doesn't exist
  output_dir_full <- dirname(output_csv)
  if (!dir.exists(output_dir_full)) {
    dir.create(output_dir_full, recursive = TRUE)
  }

  # now read the trackmate XML file
  df <- readTrackMateXML(XMLpath, slim = TRUE)
  data <- df[[1]]
  calib <- df[[2]]

  # filter traces if they have less than min_points
  # and min_points is greater than 0
  if (min_points > 0) {
    ntrace <- data %>%
      group_by(trace) %>%
      summarise(n = n()) %>%
      filter(n >= min_points)
    tr <- ntrace$trace
    data <- data[data$trace %in% c(tr), ]
  }

  # select columns for export
  if (!is.null(columns)) {
    # if user has supplied a list of columns to export, check that they include
    # trace, frame, x and y
    columns_required <- c("trace", "frame", "y", "x")
    columns <- unique(c(columns, columns_required))
    if (!all(columns %in% colnames(data))) {
      stop("Some columns specified in 'columns' are not present in the
           TrackMate XML data. Please check the column names and try again.")
    }
    data <- data[, columns, drop = FALSE]
  }

  # trace names are strings but we will re-embed them as integers and ensure
  # that they are consecutive
  data$trace <- as.integer(as.factor(data$trace))
  # micron to pixel conversion if required
  if (pixels) {
    cat("Converting x and y to pixels using calibration information from
        TrackMate XML file\n")
    data$y <- round(data$y / calib[1, 1])
    data$x <- round(data$x / calib[1, 1])
  }
  if (no_save) {
    return(data)
  } else {
    write.csv(data, file = output_csv, row.names = FALSE)
  }
}
