# Count segmented nuclei whose geometric centroids are inside/on an SCN ROI.
# install.packages(c("png", "RImageJROI", "sp"))
# The mask and ROIs must use the same original image size, crop and orientation.
# Each positive mask label represents one nucleus; 0 is background.
# All selected ROIs are treated as SCN. Overlaps count once in n_nuclei_SCN.
# Supports simple ImageJ Polygon/Freehand/Traced contours without fitted splines
# or subpixel extensions. These restrictions match classify_spots_scn.R.
# Coordinates follow ImageJ AREA ROI geometry: the top-left pixel's center is
# (0.5, 0.5), so centroid_x_px = mean(R column index) - 0.5, likewise for y.
# This differs by +0.5 from the array-index coordinates used for Spotiflow spots.
# Reference: https://imagej.net/ij/developer/api/ij/ij/gui/Roi.html

count_nuclei_in_scn <- function(mask_file, roi_file,
                                image_id = sub("_cp_masks$", "",
                                  tools::file_path_sans_ext(basename(mask_file))),
                                roi_names = NULL) {
  for (pkg in c("png", "RImageJROI", "sp")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf('Run install.packages("%s") first.', pkg))
    }
  }
  for (path in list(mask_file, roi_file)) {
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !file.exists(path) || dir.exists(path)) stop("Provide one existing input file per argument.")
  }
  if (!is.character(image_id) || length(image_id) != 1L ||
      is.na(image_id) || !nzchar(image_id)) stop("Provide a nonempty image_id.")

  # readPNG normalizes BOTH 8-bit and 16-bit PNGs to 0..1.
  im <- png::readPNG(mask_file, native = FALSE, info = TRUE)
  info <- attr(im, "info")
  if (length(dim(im)) != 2L || !identical(info$color.type, "gray") ||
      length(info$bit.depth) != 1L || !info$bit.depth %in% c(8L, 16L)) {
    stop("Expected an 8- or 16-bit grayscale label PNG, not an RGB overlay.")
  }
  mask <- round(im * (2^info$bit.depth - 1))
  storage.mode(mask) <- "integer"

  ext <- tolower(tools::file_ext(roi_file))
  if (ext == "zip") {
    rois <- RImageJROI::read.ijzip(roi_file)
  } else if (ext == "roi") {
    roi <- RImageJROI::read.ijroi(roi_file)
    rois <- setNames(list(roi), roi$name)
  } else stop("Expected an ImageJ .zip or .roi file.")
  if (!length(rois) || is.null(names(rois)) || anyNA(names(rois)) ||
      any(!nzchar(names(rois))) || anyDuplicated(names(rois))) {
    stop("ROIs must be nonempty and uniquely named.")
  }
  if (!is.null(roi_names)) {
    if (!is.character(roi_names) || !length(roi_names) || anyNA(roi_names) ||
        anyDuplicated(roi_names) || !all(roi_names %in% names(rois))) {
      stop("roi_names must contain unique names present in the ROI file.")
    }
    rois <- rois[roi_names]
  }
  summary_names <- c("SCN_union", "outside_SCN", "overlapping_ROIs", "total")
  if (any(names(rois) %in% summary_names)) stop("An ROI name conflicts with a summary row name.")
  polygons <- lapply(seq_along(rois), function(j) {
    roi <- rois[[j]]
    if (!roi$strType %in% c("polygon", "freehand", "traced") ||
        roi$subtype != 0L || roi$shapeRoiSize != 0L ||
        bitwAnd(as.integer(roi$options), 129L) != 0L) {
      stop(sprintf("ROI '%s': use a simple Polygon/Freehand/Traced contour without spline fitting or subpixel coordinates.", names(rois)[j]))
    }
    xy <- as.matrix(roi$coords)
    if (ncol(xy) != 2L || nrow(unique(xy)) < 3L || any(!is.finite(xy))) {
      stop(sprintf("Invalid polygon coordinates in ROI '%s'.", names(rois)[j]))
    }
    # RImageJROI already adds the left/top offsets; do not add them again.
    xy
  })

  # Group all foreground pixels by label in one pass, retaining sparse label IDs.
  pixels <- which(mask > 0L, arr.ind = TRUE)
  nuclei <- data.frame(mask_label = integer(), cell_id = character(),
                       area_px = integer(), centroid_x_px = double(),
                       centroid_y_px = double(), stringsAsFactors = FALSE)
  if (nrow(pixels)) {
    sums <- rowsum(cbind(area_px = rep.int(1L, nrow(pixels)),
                         x = pixels[, 2] - 0.5, y = pixels[, 1] - 0.5),
                   group = mask[pixels], reorder = TRUE)
    ids <- as.integer(rownames(sums))
    nuclei <- data.frame(
      mask_label = ids, cell_id = paste0(image_id, "_n", ids),
      area_px = as.integer(sums[, "area_px"]),
      centroid_x_px = sums[, "x"] / sums[, "area_px"],
      centroid_y_px = sums[, "y"] / sums[, "area_px"],
      row.names = NULL, stringsAsFactors = FALSE)
  }
  codes <- matrix(0L, nrow = nrow(nuclei), ncol = length(rois),
                   dimnames = list(nuclei$cell_id, names(rois)))
  for (j in seq_along(polygons)) {
    xy <- polygons[[j]]
    if (nrow(nuclei)) codes[, j] <- sp::point.in.polygon(
      nuclei$centroid_x_px, nuclei$centroid_y_px, xy[, 1], xy[, 2])
  }
  # 0=outside, 1=inside, 2=edge, 3=vertex. Include edges and vertices.
  hits <- codes > 0L
  nuclei$in_SCN <- rowSums(hits) > 0L
  nuclei$SCN_overlap <- rowSums(hits) > 1L
  nuclei$SCN_on_boundary <- rowSums(codes >= 2L) > 0L
  nuclei$SCN_roi <- vapply(seq_len(nrow(nuclei)), function(i) {
    if (!any(hits[i, ])) return(NA_character_)
    paste(names(rois)[hits[i, ]], collapse = ";")
  }, character(1))
  counts <- data.frame(
    region = c(names(rois), summary_names),
    n_nuclei = c(colSums(hits), sum(nuclei$in_SCN), sum(!nuclei$in_SCN),
                 sum(nuclei$SCN_overlap), nrow(nuclei)),
    row.names = NULL)
  list(n_nuclei_SCN = sum(nuclei$in_SCN), counts = counts, nuclei = nuclei,
       roi_membership = hits, rois = rois, image_id = image_id,
       method = "centroid; boundary included")
}
