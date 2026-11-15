# Standalone v3: source this file, then call assign_spots_to_mask_v3(...).
# install.packages(c("png", "RImageJROI", "sp"))
# Label lookup only: no nucleus expansion or nearest-cell extrapolation.
# Spots, mask and ROIs must share the same image size, crop and orientation.
#
# Keeps the v2 spots/cells/mask output and adds cells$in_SCN (logical).
# A nucleus is in_SCN when its geometric centroid is inside/on ANY selected ROI.
# Nuclei with no spots are classified too; overlapping ROIs count only once.
# roi_file: ImageJ .zip or .roi; roi_names: optional subset (NULL uses all ROIs).
# Supports simple Polygon/Freehand/Traced contours without spline/subpixel data.
#
# centroid_x_px/y_px retain v2's 0-based pixel-center coordinates: top-left=(0,0).
# For ImageJ AREA ROI membership only, add 0.5 to each centroid coordinate.
# This matches count_nuclei_in_scn.R without changing the exported v2 centroids.
# cells$in_SCN is determined by geometry, independently of input spots$in_SCN.
# If spots$in_SCN is supplied, n_spots_in_SCN still counts those flagged spots.

read_cellpose_labels_v3 <- function(mask_file) {
  if (!requireNamespace("png", quietly = TRUE)) stop('Run install.packages("png") first.')
  im <- png::readPNG(mask_file, native = FALSE, info = TRUE)
  info <- attr(im, "info")
  if (length(dim(im)) != 2L || !identical(info$color.type, "gray")) {
    stop("Expected a grayscale integer label PNG, not an RGB overlay or screenshot.")
  }
  depth <- info$bit.depth
  if (length(depth) != 1L || !depth %in% c(8L, 16L)) stop("Expected an 8- or 16-bit label PNG.")
  # readPNG returns 0..1, even for 16-bit labels. Preserve all 16 bits.
  labels <- round(im * (2^depth - 1))
  storage.mode(labels) <- "integer"
  attributes(labels) <- list(dim = dim(im))
  ids <- sort(unique(as.vector(labels)))
  if (length(ids[ids > 0L]) <= 1L) {
    warning("Mask has zero or one object ID; check that this is not a binary mask.")
  }
  labels
}

read_scn_polygons_v3 <- function(roi_file, roi_names = NULL) {
  for (pkg in c("RImageJROI", "sp")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf('Run install.packages("%s") first.', pkg))
    }
  }
  if (!is.character(roi_file) || length(roi_file) != 1L || is.na(roi_file) ||
      !file.exists(roi_file) || dir.exists(roi_file)) {
    stop("Provide one existing ImageJ ROI .zip or .roi file as roi_file.")
  }
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
  setNames(polygons, names(rois))
}

assign_spots_to_mask_v3 <- function(spots, mask_file, image_id,
                                    mask_kind = c("nuclei", "cells"),
                                    roi_file, roi_names = NULL) {
  mask_kind <- match.arg(mask_kind)
  if (length(image_id) != 1L || is.na(image_id) || !nzchar(image_id)) stop("Provide a unique image_id.")
  if (!all(c("x_px", "y_px") %in% names(spots)) ||
      !is.numeric(spots$x_px) || !is.numeric(spots$y_px)) stop("Need numeric x_px and y_px columns.")
  reserved <- c("mask_label", "cell_id", "mask_assignment", "mask_image_id")
  if (any(reserved %in% names(spots))) stop("Assignment columns already exist; use the unassigned table.")
  if (missing(roi_file)) stop("Provide roi_file to classify nucleus/cell centroids in SCN.")
  polygons <- read_scn_polygons_v3(roi_file, roi_names)
  mask <- read_cellpose_labels_v3(mask_file)
  # Coordinates: x=column, y=row. Pixel coordinates 0-based; R indices 1-based.
  col_idx <- floor(spots$x_px + 0.5) + 1
  row_idx <- floor(spots$y_px + 0.5) + 1
  finite <- is.finite(spots$x_px) & is.finite(spots$y_px)
  valid <- finite & spots$x_px >= 0 & spots$y_px >= 0 &
    col_idx >= 1 & col_idx <= ncol(mask) & row_idx >= 1 & row_idx <= nrow(mask)
  valid[is.na(valid)] <- FALSE
  label <- rep(NA_integer_, nrow(spots))
  label[valid] <- mask[cbind(row_idx[valid], col_idx[valid])]
  assigned <- !is.na(label) & label > 0L
  spots$mask_label <- label
  spots$cell_id <- rep(NA_character_, nrow(spots))
  spots$cell_id[assigned] <- paste0(image_id, "_", if (mask_kind == "nuclei") "n" else "c", label[assigned])
  spots$mask_assignment <- ifelse(finite, "out_of_bounds", "invalid_coordinate")
  spots$mask_assignment[valid] <- "background"
  spots$mask_assignment[assigned] <- if (mask_kind == "nuclei") "inside_nucleus" else "inside_cell"
  spots$mask_image_id <- rep(image_id, nrow(spots))
  ids <- sort(unique(as.vector(mask)))
  ids <- ids[ids > 0L]
  cell_ids <- character(length(ids))
  if (length(ids)) cell_ids <- paste0(image_id, "_", if (mask_kind == "nuclei") "n" else "c", ids)
  counts <- data.frame(mask_label = ids, cell_id = cell_ids,
                       n_spots = as.integer(table(factor(label[assigned], levels = ids))))
  if ("in_SCN" %in% names(spots)) {
    scn <- assigned & !is.na(spots$in_SCN) & spots$in_SCN
    counts$n_spots_in_SCN <- as.integer(table(factor(label[scn], levels = ids)))
  }

  counts$centroid_x_px <- rep(NA_real_, length(ids))
  counts$centroid_y_px <- rep(NA_real_, length(ids))
  if (length(ids)) {
    pixels <- which(mask > 0L, arr.ind = TRUE)
    sums <- rowsum(cbind(n_pixels = rep.int(1L, nrow(pixels)),
                         x = pixels[, 2] - 1, y = pixels[, 1] - 1),
                   group = mask[pixels], reorder = TRUE)
    # Match by label, including masks with nonconsecutive or spatially unordered IDs.
    sums <- sums[match(as.character(ids), rownames(sums)), , drop = FALSE]
    counts$centroid_x_px <- sums[, "x"] / sums[, "n_pixels"]
    counts$centroid_y_px <- sums[, "y"] / sums[, "n_pixels"]
  }

  counts$in_SCN <- rep(FALSE, length(ids))
  if (length(ids)) {
    roi_x <- counts$centroid_x_px + 0.5
    roi_y <- counts$centroid_y_px + 0.5
    for (xy in polygons) {
      # 0=outside, 1=inside, 2=edge, 3=vertex. Include edges and vertices.
      inside <- sp::point.in.polygon(roi_x, roi_y, xy[, 1], xy[, 2]) > 0L
      counts$in_SCN <- counts$in_SCN | inside
    }
  }
  list(spots = spots, cells = counts, mask = mask)
}
