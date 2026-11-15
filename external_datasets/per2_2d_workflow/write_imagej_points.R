# Export a data.frame/matrix with x_px and y_px columns as an ImageJ POINT ROI.
# Base R only. Input coordinates are 0-based pixels; do not add 1 or scale to um.
# Subpixel positions are stored as 32-bit floats, per the ImageJ ROI format.
# Format reference: https://github.com/imagej/ImageJ/blob/master/ij/io/RoiEncoder.java
#
# source("/path/to/write_imagej_points.R")
# points <- read.delim("/path/to/points.tsv")
# write_imagej_points(points, "~/Downloads/49B_Per2_spots.roi")
# write_imagej_points(subset(points, score >= 0.6), "~/Downloads/49B_Per2_spots_06.roi")

write_imagej_points <- function(points, file, overwrite = FALSE) {
  if (!(is.data.frame(points) || is.matrix(points)) ||
      !all(c("x_px", "y_px") %in% colnames(points))) {
    stop("points must be a data.frame or matrix with x_px and y_px columns.")
  }
  x <- points[, "x_px"]
  y <- points[, "y_px"]
  if (!is.numeric(x) || !is.numeric(y) || any(!is.finite(c(x, y)))) {
    stop("x_px and y_px must contain finite numeric coordinates.")
  }
  n <- length(x)
  if (!n) stop("No spots to export; no ROI file was written.")
  # Keep legacy integer bounds within ImageJ's unambiguous coordinate range.
  if (any(c(x, y) < 0) || any(c(x, y) >= 60000)) {
    stop("This writer supports 0 <= x_px, y_px < 60000 pixels.")
  }
  file <- path.expand(file)
  if (tolower(tools::file_ext(file)) != "roi") stop("Output filename must end in .roi")
  if (!dir.exists(dirname(file))) stop("Output directory does not exist: ", dirname(file))
  if (dir.exists(file) || (file.exists(file) && !overwrite)) {
    stop("Output already exists; choose another name or set overwrite = TRUE.")
  }

  # All offsets below are 0-based byte offsets in the ImageJ format.
  header2 <- 64 + 12*n
  bytes <- raw(header2 + 64)
  put <- function(offset, value) {
    bytes[seq.int(offset + 1, length.out = length(value))] <<- value
  }
  u16 <- function(value) as.raw(as.vector(rbind(value %/% 256, value %% 256)))
  int32 <- function(value) writeBin(as.integer(value), raw(), size = 4, endian = "big")
  float32 <- function(value) writeBin(as.double(value), raw(), size = 4, endian = "big")
  left <- floor(min(x)); top <- floor(min(y))

  put(0, charToRaw("Iout"))
  put(4, u16(228))                           # Supports subpixel POINT ROIs.
  put(6, as.raw(10))                         # POINT, not a polygon.
  put(8, u16(top))
  put(10, u16(left))
  put(12, u16(floor(max(y)) + 1))
  put(14, u16(floor(max(x)) + 1))
  if (n <= 65535) put(16, u16(n)) else put(18, int32(n))
  put(34, u16(2))                            # Point marker size.
  put(40, as.raw(c(255, 0, 255, 0)))         # Opaque green (ARGB).
  put(50, u16(128))                          # SUB_PIXEL_RESOLUTION; no labels.
  put(60, int32(header2))
  put(64, u16(floor(x) - left))
  put(64 + 2*n, u16(floor(y) - top))
  put(64 + 4*n, float32(x))                  # Absolute x coordinates.
  put(64 + 8*n, float32(y))                  # Absolute y coordinates.
  # Empty header2: ImageJ uses the filename as the ROI name.
  writeBin(bytes, file)
  message("Wrote ", n, " points to ", normalizePath(file))
  invisible(file)
}
