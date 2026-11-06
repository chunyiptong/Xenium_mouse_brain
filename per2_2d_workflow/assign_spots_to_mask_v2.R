# Label lookup only: no nucleus expansion or nearest-cell extrapolation.
# Coordinates must refer to the same unscaled/unrotated image as the mask.
# Standalone v2: source this file, then call assign_spots_to_mask_v2(...).
# Returns the original spots/cells/mask list, with two new columns in cells:
#   centroid_x_px = mean(column index - 1) across ALL pixels of each label
#   centroid_y_px = mean(row index - 1) across ALL pixels of each label
# These are geometric centroids, not spot averages or intensity-weighted centers.
# Coordinates match spots$x_px/y_px: 0-based, x rightwards, y downwards.
# Nuclei with zero assigned spots also receive centroids. Values are not rounded.
# For ImageJ AREA ROI coordinates (as used by count_nuclei_in_scn.R), add 0.5
# to both centroid coordinates; the top-left pixel center there is (0.5, 0.5).

read_cellpose_labels_v2 <- function(mask_file) {
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

assign_spots_to_mask_v2 <- function(spots, mask_file, image_id,
                                    mask_kind = c("nuclei", "cells")) {
  mask_kind <- match.arg(mask_kind)
  if (length(image_id) != 1L || is.na(image_id) || !nzchar(image_id)) stop("Provide a unique image_id.")
  if (!all(c("x_px", "y_px") %in% names(spots)) ||
      !is.numeric(spots$x_px) || !is.numeric(spots$y_px)) stop("Need numeric x_px and y_px columns.")
  reserved <- c("mask_label", "cell_id", "mask_assignment", "mask_image_id")
  if (any(reserved %in% names(spots))) stop("Assignment columns already exist; use the unassigned table.")
  mask <- read_cellpose_labels_v2(mask_file)
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
  list(spots = spots, cells = counts, mask = mask)
}
