# ImageJ Polygon/Freehand ROI -> membership of 2D Spotiflow spots.
# install.packages(c("RImageJROI", "sp"))
# source("/Users/chunyiptong/Downloads/20261010/per2_2d_workflow/classify_spots_scn.R")
# result <- classify_spots_scn("points.tsv", "SCN_ROIs.zip", "new_output_folder")
# Also accepts one .roi. All coordinates must be in the SAME original pixel grid.
# Borders count as inside. Overlapping ROIs are flagged and never duplicate rows.

classify_spots_scn <- function(spots_file, roi_file, output_dir,
                               x_col = "x_px", y_col = "y_px", sep = "\t", gene) {
  for (pkg in c("RImageJROI", "sp")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf('Run install.packages("%s") first.', pkg))
    }
  }
  spots <- read.table(spots_file, header = TRUE, sep = sep,
                      check.names = FALSE, stringsAsFactors = FALSE,
                      comment.char = "", na.strings = c("NA", "nan", ""))
  if (!all(c(x_col, y_col) %in% names(spots))) stop("Missing X/Y column names.")
  if (!is.numeric(spots[[x_col]]) || !is.numeric(spots[[y_col]]) ||
      any(!is.finite(spots[[x_col]])) || any(!is.finite(spots[[y_col]]))) {
    stop("X/Y must be finite numeric pixel coordinates. Check missing values and separator.")
  }
  ext <- tolower(tools::file_ext(roi_file))
  if (ext == "zip") {
    rois <- RImageJROI::read.ijzip(roi_file)
  } else if (ext == "roi") {
    roi <- RImageJROI::read.ijroi(roi_file)
    rois <- setNames(list(roi), roi$name)
  } else stop("Expected an ImageJ .roi or ROI .zip file.")
  if (!length(rois) || is.null(names(rois)) || anyDuplicated(names(rois))) {
    stop("ROIs must be nonempty and uniquely named.")
  }
  polygons <- lapply(rois, function(roi) {
    if (!roi$strType %in% c("polygon", "freehand", "traced")) {
      stop("Use simple Polygon/Freehand/Traced ROI contours, not lines, points or composite ROIs.")
    }
    # RImageJROI does not reconstruct fitted splines/subpixel ROI extensions.
    if (bitwAnd(as.integer(roi$options), 129L) != 0L) {
      stop("Spline-fit or subpixel ROI detected; use a reader that preserves this geometry.")
    }
    xy <- as.matrix(roi$coords)
    if (ncol(xy) != 2L || nrow(unique(xy)) < 3L || any(!is.finite(xy))) {
      stop("Invalid polygon coordinates.")
    }
    # coords already include the ROI left/top offsets: do not add them again.
    xy
  })
  reserved <- c(paste0("in_", names(rois)), "in_SCN", "SCN_roi", "SCN_overlap", "SCN_on_boundary")
  if (anyDuplicated(reserved) || any(reserved %in% names(spots))) {
    stop("Output column names conflict with input columns or ROI names.")
  }
  codes <- matrix(0L, nrow = nrow(spots), ncol = length(polygons),
                  dimnames = list(NULL, names(rois)))
  for (j in seq_along(polygons)) {
    xy <- polygons[[j]]
    if (nrow(spots)) codes[, j] <- sp::point.in.polygon(
      point.x = spots[[x_col]], point.y = spots[[y_col]],
      pol.x = xy[, 1], pol.y = xy[, 2])
  }
  # 0=outside, 1=inside, 2=edge, 3=vertex.
  hits <- codes > 0L
  for (j in seq_along(polygons)) spots[[paste0("in_", names(rois)[j])]] <- hits[, j]
  spots$in_SCN <- rowSums(hits) > 0L
  spots$SCN_overlap <- rowSums(hits) > 1L
  spots$SCN_on_boundary <- rowSums(codes >= 2L) > 0L
  spots$SCN_roi <- vapply(seq_len(nrow(spots)), function(i) {
    if (!any(hits[i, ])) return(NA_character_)
    paste(names(rois)[hits[i, ]], collapse = ";")
  }, character(1))
  counts <- data.frame(
    region = c(names(rois), "SCN_union", "outside_SCN", "overlapping_ROIs", "total"),
    n_spots = c(colSums(hits), sum(spots$in_SCN), sum(!spots$in_SCN),
                sum(spots$SCN_overlap), nrow(spots)))
#  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
#    stop("Output folder is not empty. Choose a new directory to preserve prior results.")
#  }
#  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(spots, file.path(output_dir, sprintf("%s_spots_with_SCN.txt", gene)), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
#  write.table(spots[spots$in_SCN, , drop = FALSE], file.path(output_dir, sprintf("%s_spots_with_SCN.txt", gene)),
#              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
#  write.table(counts, file.path(output_dir, "SCN_counts.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
#  vertices <- do.call(rbind, lapply(seq_along(polygons), function(j) {
#    data.frame(roi_id = names(rois)[j], vertex_order = seq_len(nrow(polygons[[j]])),
#               x_px = polygons[[j]][, 1], y_px = polygons[[j]][, 2])
#  }))
#  write.table(vertices, file.path(output_dir, "SCN_vertices_pixels.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
#  allxy <- do.call(rbind, polygons)
#  png(file.path(output_dir, "SCN_spots_check.png"), width = 1800, height = 1200, res = 150)
#  tryCatch({
#    plot(NA, xlim = range(c(allxy[, 1], spots[[x_col]])),
#         ylim = rev(range(c(allxy[, 2], spots[[y_col]]))), asp = 1,
#         xlab = "x (original pixels)", ylab = "y (original pixels; downwards)",
#         main = "SCN contours and Per2 spots")
#    graphics::points(spots[[x_col]], spots[[y_col]], pch = 16, cex = .25,
#                     col = ifelse(spots$in_SCN, "#d62728", "#bdbdbd"))
#    for (j in seq_along(polygons)) {
#      xy <- polygons[[j]]
#      lines(rbind(xy, xy[1, ]), col = "#1565c0", lwd = 2)
#      text(mean(xy[, 1]), min(xy[, 2]), names(rois)[j], col = "#1565c0", pos = 3)
#    }
#    legend("bottomleft", c("Inside/on SCN", "Outside SCN"),
#           col = c("#d62728", "#bdbdbd"), pch = 16, bty = "n")
#  }, finally = dev.off())
#  writeLines(c(paste("spots:", normalizePath(spots_file)),
#               paste("ROI:", normalizePath(roi_file)),
#               "Coordinates: original 0-based pixels; no rescaling, swapping or offset added.",
#               "Boundary included. Multiple ROI matches retained as semicolon-separated names.",
#               "ROI names are used verbatim; anatomical left/right is not inferred.",
#               "This classification does not assign cell IDs or remeasure intensity.",
#               capture.output(sessionInfo())), file.path(output_dir, "run_info.txt"))
#  print(counts, row.names = FALSE)
  invisible(list(spots = spots, counts = counts, rois = rois))
}
