# Read outputs from spotiflow_to_r.py. Requires install.packages("jsonlite").
# source("/path/to/read_spotiflow.R")
# result <- read_spotiflow("/path/to/results")
# points <- result$points
# details <- result$details

read_spotiflow <- function(result_dir) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop('Run install.packages("jsonlite") first.')
  }
  points <- read.delim(file.path(result_dir, "points.tsv"),
                       na.strings = c("nan", "NA"), colClasses = "numeric")
  con <- gzfile(file.path(result_dir, "details.json.gz"), open = "rb")
  on.exit(close(con))
  details <- jsonlite::fromJSON(con, simplifyMatrix = TRUE)
  # Empty point-wise outputs should remain usable numeric vectors in R.
  for (name in intersect(c("prob", "intens"), names(details))) {
    if (is.list(details[[name]]) && !length(details[[name]])) details[[name]] <- numeric()
  }
  list(points = points, details = details)
}
