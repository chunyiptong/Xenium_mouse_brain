#!/usr/bin/env Rscript
# 2D Per2 spot photometry and mask assignment. See README_R.md.
# All exported coordinates are 0-based (x rightward, y downward).
# Only Spotiflow inference uses Python, through reticulate.

require_package <- function(name) {
  if (!requireNamespace(name, quietly = TRUE)) {
    stop(sprintf("Missing R package '%s'; run install.packages('%s')", name, name),
         call. = FALSE)
  }
}

read_image <- function(path, mask = FALSE) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tif", "tiff")) {
    require_package("ijtiff")
    img <- ijtiff::read_tif(path, msg = FALSE)
    if (is.list(img) || length(dim(img)) != 4L || any(dim(img)[3:4] != 1L)) {
      stop(path, ": expected a single grayscale 2D plane")
    }
    if (mask && !identical(attr(img, "SampleFormat"), "unsigned integer data") &&
        !identical(attr(img, "SampleFormat"), "signed integer data")) {
      stop("Mask TIFF must use integer samples")
    }
    img <- matrix(as.numeric(img), nrow = dim(img)[1], ncol = dim(img)[2])
  } else if (ext == "png") {
    require_package("png")
    img <- png::readPNG(path, native = FALSE, info = TRUE)
    info <- attr(img, "info")
    if (!identical(info$color.type, "gray") || length(dim(img)) != 2L) {
      stop(path, ": expected grayscale PNG, without a palette or alpha channel")
    }
    # readPNG scales both 8-bit and 16-bit samples to [0,1]. Restore raw values.
    img <- round(img * (2^info$bit.depth - 1))
  } else if (ext %in% c("jpg", "jpeg") && !mask) {
    require_package("jpeg")
    img <- round(jpeg::readJPEG(path, native = FALSE) * 255)
  } else {
    stop("Supported inputs: grayscale TIFF/PNG (also JPEG for raw images)")
  }
  if (length(dim(img)) != 2L || any(dim(img) == 0L)) {
    stop(path, ": expected a single grayscale 2D plane")
  }
  img
}

# Decode the simple area-ROI subset of ImageJ's big-endian binary format.
# Includes embedded UTF-16 names and subpixel coordinates; no ZIP extraction.
# Format: https://github.com/imagej/ImageJ/blob/master/ij/io/RoiDecoder.java
decode_roi <- function(bytes, fallback) {
  take <- function(offset, n) {
    if (!is.finite(offset) || !is.finite(n) || offset < 0 || n < 0 ||
        offset + n > length(bytes)) stop("Truncated or invalid ImageJ ROI")
    if (!n) return(raw())
    bytes[seq.int(offset + 1, length.out = n)]
  }
  uint <- function(offset, size = 2L) {
    sum(as.integer(take(offset, size)) * 256^rev(seq_len(size) - 1L))
  }
  coord_origin <- function(offset) {
    n <- uint(offset)
    if (n >= 65536 - 5000) n <- n - 65536
    n
  }
  if (rawToChar(take(0, 4)) != "Iout") stop("Not an ImageJ ROI")
  version <- uint(4)
  type <- uint(6, 1)
  if (!type %in% c(0, 7, 8) || uint(36, 4) != 0) {
    stop("Use simple Polygon/Freehand/Traced ROIs; composite ROIs are unsupported")
  }
  n <- uint(16)
  if (n == 0) n <- uint(18, 4)
  if (n < 3 || n > (length(bytes) - 64) / 4) stop("Invalid polygon ROI")
  if (version >= 222 && bitwAnd(as.integer(uint(50)), 128L) != 0L) {
    x <- readBin(take(64 + 4*n, 4*n), numeric(), n, size = 4, endian = "big")
    y <- readBin(take(64 + 8*n, 4*n), numeric(), n, size = 4, endian = "big")
  } else {
    x <- readBin(take(64, 2*n), integer(), n, size = 2, signed = FALSE,
                 endian = "big") + coord_origin(10)
    y <- readBin(take(64 + 2*n, 2*n), integer(), n, size = 2, signed = FALSE,
                 endian = "big") + coord_origin(8)
  }
  if (any(!is.finite(c(x, y)))) stop("Invalid polygon coordinates")
  name <- fallback
  hdr2 <- uint(60, 4)
  if (hdr2 > 0) {
    offset <- uint(hdr2 + 16, 4)
    len <- uint(hdr2 + 20, 4)
    if (offset > 0 && len > 0) {
      # Decode UTF-16BE, including surrogate pairs, without embedded NUL strings.
      units <- readBin(take(offset, 2*len), integer(), len, size = 2,
                       signed = FALSE, endian = "big")
      codepoints <- integer()
      j <- 1L
      while (j <= length(units)) {
        u <- units[j]
        if (u >= 55296L && u <= 56319L && j < length(units) &&
            units[j+1L] >= 56320L && units[j+1L] <= 57343L) {
          u <- 65536L + (u - 55296L)*1024L + units[j+1L] - 56320L
          j <- j + 1L
        }
        codepoints <- c(codepoints, u)
        j <- j + 1L
      }
      name <- intToUtf8(codepoints)
    }
  }
  if (is.na(name) || !nzchar(name)) name <- fallback
  list(name = name, vertices = cbind(x, y))
}

read_rois <- function(path) {
  if (tolower(tools::file_ext(path)) == "zip") {
    entries <- utils::unzip(path, list = TRUE)
    entries <- entries[grepl("\\.roi$", entries$Name, ignore.case = TRUE), , drop = FALSE]
    rois <- lapply(seq_len(nrow(entries)), function(i) {
      con <- unz(path, entries$Name[i], open = "rb")
      on.exit(close(con))
      decode_roi(readBin(con, raw(), n = entries$Length[i]), paste0("SCN_", i))
    })
  } else {
    rois <- list(decode_roi(readBin(path, raw(), n = file.info(path)$size), "SCN_1"))
  }
  if (!length(rois) || anyDuplicated(vapply(rois, `[[`, "", "name"))) {
    stop("At least one ROI is required, with unique names")
  }
  rois
}

roi_membership <- function(xy, rois) {
  names <- rep("", nrow(xy))
  if (!nrow(xy)) return(names)
  for (roi in rois) {
    v <- roi$vertices
    inside <- on_edge <- rep(FALSE, nrow(xy))
    j <- nrow(v)
    for (i in seq_len(nrow(v))) {
      dx <- v[i, 1] - v[j, 1]; dy <- v[i, 2] - v[j, 2]
      px <- xy[, 1] - v[j, 1]; py <- xy[, 2] - v[j, 2]
      len2 <- dx*dx + dy*dy
      if (len2 > 0) {
        t <- pmax(0, pmin(1, (px*dx + py*dy)/len2))
        on_edge <- on_edge | ((px - t*dx)^2 + (py - t*dy)^2 <= (5e-8)^2)
      }
      if (dy != 0) {
        crossing <- (v[i, 2] > xy[, 2]) != (v[j, 2] > xy[, 2])
        inside <- xor(inside, crossing & xy[, 1] < v[j, 1] + py*dx/dy)
      }
      j <- i
    }
    hit <- inside | on_edge
    names[hit] <- ifelse(nzchar(names[hit]), paste(names[hit], roi$name, sep = ";"), roi$name)
  }
  names
}

# Separable exact Euclidean distance transform in physical units.
# Equal-distance ties select the lowest x, then the lowest y source pixel.
nearest_parabola <- function(f, spacing) {
  candidates <- which(is.finite(f))
  if (!length(candidates)) return(rep(NA_integer_, length(f)))
  v <- integer(length(candidates))
  z <- rep(Inf, length(candidates) + 1L)
  k <- 1L; v[k] <- candidates[1]; z[1] <- -Inf
  a <- spacing^2
  for (q in candidates[-1L]) {
    repeat {
      p <- v[k]
      s <- ((f[q] - f[p])/a + (q - 1)^2 - (p - 1)^2) / (2*(q - p))
      if (k == 1L || s > z[k]) break
      k <- k - 1L
    }
    k <- k + 1L; v[k] <- q; z[k] <- s; z[k+1L] <- Inf
  }
  v[findInterval(seq_along(f) - 1L, z[seq_len(k+1L)], left.open = TRUE)]
}

expand_labels <- function(labels, sx, sy, radius) {
  if (radius == 0 || !any(labels == 0)) return(labels)
  h <- nrow(labels); w <- ncol(labels)
  ys <- seq_len(h)
  nearest_y <- matrix(NA_integer_, h, w)
  distance_y2 <- matrix(Inf, h, w)
  for (x in seq_len(w)) {
    seeds <- which(labels[, x] > 0)
    if (!length(seeds)) next
    pos <- findInterval(ys, seeds)
    lo <- seeds[pmax(1L, pos)]
    hi <- seeds[pmin(length(seeds), pos + 1L)]
    nearest_y[, x] <- ifelse(abs(ys - lo) <= abs(ys - hi), lo, hi)
    distance_y2[, x] <- ((ys - nearest_y[, x])*sy)^2
  }
  assigned <- labels
  xs <- seq_len(w)
  for (y in ys) {
    source_x <- nearest_parabola(distance_y2[y, ], sx)
    distance <- sqrt(distance_y2[y, source_x] + ((xs - source_x)*sx)^2)
    hit <- which(labels[y, ] == 0 & distance <= radius)
    if (length(hit)) {
      xx <- source_x[hit]
      yy <- nearest_y[cbind(rep(y, length(xx)), xx)]
      assigned[y, hit] <- labels[cbind(yy, xx)]
    }
  }
  assigned
}

mask_boundary <- function(labels) {
  h <- nrow(labels); w <- ncol(labels)
  low <- high <- labels
  for (dy in -1:1) for (dx in -1:1) {
    tile <- labels[pmax(1L, pmin(h, seq_len(h) + dy)),
                   pmax(1L, pmin(w, seq_len(w) + dx)), drop = FALSE]
    low <- pmin(low, tile); high <- pmax(high, tile)
  }
  low != high
}

empty_spots <- function() {
  data.frame(image_id = character(), spot_id = character(), roi_id = character(),
    roi_overlap = logical(), x_px = numeric(), y_px = numeric(), x_um = numeric(),
    y_um = numeric(), score = numeric(), peak_raw = numeric(), background_local = numeric(),
    integrated_bgsub = numeric(), aperture_n_pixels = integer(), background_n_pixels = integer(),
    cell_id = character(), mask_label = integer(), assignment = character(),
    near_mask_boundary = logical(), aperture_overlaps_other_spot = logical(), edge_window = logical())
}

analyze <- function(raw, labels, points, scores, rois, args) {
  sx <- args$pixel_size_x; sy <- args$pixel_size_y
  assigned <- expand_labels(labels, sx, sy, args$expand_um)
  boundary <- mask_boundary(assigned)
  roi_ids <- roi_membership(points[, 2:1, drop = FALSE], rois)
  h <- nrow(raw); w <- ncol(raw)
  mode <- if (args$mask_kind == "cells") "cell_mask" else if (args$expand_um > 0) {
    "estimated_cell_from_nucleus"
  } else "nucleus_only"
  # Spatial bins avoid an all-pairs distance matrix for dense spot detections.
  bin_width <- args$bg_outer + args$spot_radius
  bins <- floor(points/bin_width)
  keys <- paste(bins[, 1], bins[, 2], sep = ":")
  buckets <- list2env(split(seq_len(nrow(points)), keys), parent = emptyenv())
  rows <- vector("list", nrow(points))
  for (i in seq_len(nrow(points))) {
    y <- points[i, 1]; x <- points[i, 2]
    iy <- max(1L, min(h, floor(y + .5) + 1L))
    ix <- max(1L, min(w, floor(x + .5) + 1L))
    radius <- args$bg_outer
    yy <- seq.int(max(0, floor(y-radius)), min(h-1, ceiling(y+radius)))
    xx <- seq.int(max(0, floor(x-radius)), min(w-1, ceiling(x+radius)))
    d2 <- outer((yy-y)^2, (xx-x)^2, "+")
    aperture <- d2 <= args$spot_radius^2
    annulus <- d2 >= args$bg_inner^2 & d2 <= radius^2
    nearby_keys <- as.vector(outer(bins[i, 1] + (-1:1), bins[i, 2] + (-1:1),
                                  function(a, b) paste(a, b, sep = ":")))
    neighbours <- unlist(mget(nearby_keys, envir = buckets, ifnotfound = list(integer())),
                         use.names = FALSE)
    neighbours <- neighbours[neighbours != i]
    nd2 <- (points[neighbours, 1]-y)^2 + (points[neighbours, 2]-x)^2
    neighbours <- neighbours[nd2 <= bin_width^2]
    for (j in neighbours) {
      annulus <- annulus & outer((yy-points[j, 1])^2, (xx-points[j, 2])^2, "+") > args$spot_radius^2
    }
    tile <- raw[yy+1L, xx+1L, drop = FALSE]
    n_bg <- sum(annulus); n_signal <- sum(aperture)
    background <- if (n_bg >= 10) stats::median(tile[annulus]) else NaN
    integrated <- if (n_signal) sum(tile[aperture]) - background*n_signal else NaN
    overlap <- any((points[neighbours, 1]-y)^2 + (points[neighbours, 2]-x)^2 < (2*args$spot_radius)^2)
    label <- assigned[iy, ix]
    rows[[i]] <- data.frame(image_id = args$image_id, spot_id = sprintf("p%06d", i),
      roi_id = roi_ids[i], roi_overlap = grepl(";", roi_ids[i], fixed = TRUE),
      x_px = x, y_px = y, x_um = x*sx, y_um = y*sy, score = scores[i],
      peak_raw = if (n_signal) max(tile[aperture]) else NaN,
      background_local = background, integrated_bgsub = integrated,
      aperture_n_pixels = n_signal, background_n_pixels = n_bg,
      cell_id = if (label) sprintf("%s_c%05d", args$image_id, as.integer(label)) else "",
      mask_label = label, assignment = if (label) mode else "unassigned",
      near_mask_boundary = boundary[iy, ix], aperture_overlaps_other_spot = overlap,
      edge_window = x-radius < 0 || y-radius < 0 || x+radius > w-1 || y+radius > h-1)
  }
  spots <- if (length(rows)) do.call(rbind, rows) else empty_spots()
  pos <- which(assigned > 0)
  object_labels <- assigned[pos]
  sums <- rowsum(cbind(x = (pos-1L) %/% h, y = (pos-1L) %% h, n = 1), object_labels)
  ids <- as.integer(rownames(sums))
  centers <- sums[, c("x", "y"), drop = FALSE]/sums[, "n"]
  cell_roi <- roi_membership(centers, rois)
  total_counts <- tabulate(match(spots$mask_label, ids), nbins = length(ids))
  scn_counts <- tabulate(match(spots$mask_label[nzchar(spots$roi_id)], ids), nbins = length(ids))
  cells <- data.frame(image_id = args$image_id, cell_id = sprintf("%s_c%05d", args$image_id, ids),
    mask_label = ids, assignment_method = mode, roi_at_centroid = cell_roi,
    centroid_x_px = centers[, 1], centroid_y_px = centers[, 2], area_um2 = sums[, "n"]*sx*sy,
    n_spots_all = total_counts, n_spots_in_scn = scn_counts, row.names = NULL)
  list(spots = spots, cells = cells, assigned = assigned)
}

write_table <- function(path, data) {
  # Match Python's True/False and nan spellings for downstream compatibility.
  data[] <- lapply(data, function(x) if (is.logical(x)) ifelse(x, "True", "False") else x)
  utils::write.table(data, path, sep = "\t", quote = TRUE, qmethod = "double",
                     row.names = FALSE, na = "nan", fileEncoding = "UTF-8")
}

write_overlay <- function(path, raw, assigned, points, spots, rois) {
  h <- nrow(raw); w <- ncol(raw)
  limits <- stats::quantile(raw, c(.01, .998), names = FALSE, type = 7)
  gray <- if (limits[2] > limits[1]) pmax(0, pmin(1, (raw-limits[1])/diff(limits))) else raw*0
  edge <- mask_boundary(assigned)
  red <- green <- blue <- gray
  red[edge] <- .4*gray[edge]
  green[edge] <- blue[edge] <- .4*gray[edge] + .6
  raster <- as.raster(matrix(grDevices::rgb(red, green, blue), h, w))
  grDevices::png(path, width = 2400, height = 1600, res = 200)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(4, 4, 3, 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(-.5, w-.5), ylim = c(h-.5, -.5), asp = 1, xaxs = "i", yaxs = "i")
  graphics::rasterImage(raster, -.5, h-.5, w-.5, -.5, interpolate = FALSE)
  for (roi in rois) {
    v <- rbind(roi$vertices, roi$vertices[1, ])
    graphics::lines(v[, 1], v[, 2], col = "yellow", lwd = 1)
    graphics::text(v[1, 1], v[1, 2], roi$name, col = "yellow", cex = .65, adj = c(0, 0))
  }
  if (nrow(points)) graphics::points(points[, 2], points[, 1], pch = 1, cex = .5, lwd = .5,
                                    col = ifelse(nzchar(spots$cell_id), "green", "magenta"))
  graphics::axis(1); graphics::axis(2); graphics::box()
  graphics::title(main = "2D Per2: assigned=green, unassigned=magenta, labels=cyan, SCN=yellow",
                  xlab = "x (0-based pixels)", ylab = "y (0-based pixels)", cex.main = .85)
}

usage <- function() {
  cat(paste0("2D Per2 analysis (R)\n\n",
    "Rscript analyze_per2_2d.R --per2 IMAGE --mask MASK --mask-kind cells|nuclei\n",
    "  --rois ROI.zip --image-id ID --image-mode single_plane|max_projection\n",
    "  --source-description TEXT --pixel-size-x UM --pixel-size-y UM --out DIR\n\n",
    "Options: --expand-um 0 --spot-radius 2 --bg-inner 3 --bg-outer 5\n",
    "         --prob-thresh NUMBER --points-tsv FILE --help\n\n",
    "With --points-tsv, all analysis runs in R. Columns: y_px,x_px,score (tab-separated).\n",
    "Otherwise reticulate calls the Python Spotiflow general model.\n"))
}

parse_args <- function(argv) {
  required <- c("per2", "mask", "mask_kind", "rois", "image_id", "image_mode",
                "source_description", "pixel_size_x", "pixel_size_y", "out")
  args <- c(setNames(rep(list(NULL), length(required)), required),
            list(expand_um = 0, spot_radius = 2, bg_inner = 3, bg_outer = 5,
                 prob_thresh = NULL, points_tsv = NULL))
  if (any(argv %in% c("--help", "-h"))) { usage(); return(NULL) }
  i <- 1L
  while (i <= length(argv)) {
    token <- argv[i]
    if (!startsWith(token, "--")) stop("Expected --option, got: ", token)
    pair <- strsplit(substring(token, 3), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", pair[1], fixed = TRUE)
    if (!key %in% names(args)) stop("Unknown option: ", token)
    if (grepl("=", token, fixed = TRUE)) {
      value <- substring(token, regexpr("=", token, fixed = TRUE)[1] + 1L)
    } else {
      i <- i + 1L
      if (i > length(argv) || startsWith(argv[i], "--")) stop("Missing value for ", token)
      value <- argv[i]
    }
    args[[key]] <- value
    i <- i + 1L
  }
  missing <- required[vapply(args[required], function(x) is.null(x) || !nzchar(x), logical(1))]
  if (length(missing)) stop("Missing required options: ", paste(gsub("_", "-", missing), collapse = ", "))
  numeric_args <- c("pixel_size_x", "pixel_size_y", "expand_um", "spot_radius", "bg_inner", "bg_outer", "prob_thresh")
  for (key in numeric_args) {
    if (is.null(args[[key]])) next
    value <- suppressWarnings(as.numeric(args[[key]]))
    if (length(value) != 1L || !is.finite(value)) stop(key, " must be a finite number")
    args[[key]] <- value
  }
  if (!args$mask_kind %in% c("cells", "nuclei")) stop("mask-kind must be cells or nuclei")
  if (!args$image_mode %in% c("single_plane", "max_projection")) stop("Invalid image-mode")
  if (min(args$pixel_size_x, args$pixel_size_y) <= 0 || args$expand_um < 0) {
    stop("Pixel sizes must be positive and expansion nonnegative")
  }
  if (!(1 <= args$spot_radius && args$spot_radius < args$bg_inner && args$bg_inner < args$bg_outer)) {
    stop("Require 1 <= spot-radius < bg-inner < bg-outer")
  }
  if (args$mask_kind == "cells" && args$expand_um > 0) stop("Expansion is only supported for nuclei masks")
  if (!is.null(args$prob_thresh) && (args$prob_thresh < 0 || args$prob_thresh > 1)) {
    stop("prob-thresh must be between 0 and 1")
  }
  args
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  args <- parse_args(argv)
  if (is.null(args)) return(invisible(NULL))
  if (file.exists(args$out) && (!dir.exists(args$out) ||
      length(list.files(args$out, all.files = TRUE, no.. = TRUE)))) {
    stop("Output directory is nonempty; use a new directory to preserve previous results")
  }
  require_package("ijtiff"); require_package("jsonlite")
  raw <- read_image(args$per2)
  labels <- read_image(args$mask, mask = TRUE)
  if (!identical(dim(raw), dim(labels)) || any(!is.finite(raw))) {
    stop("Image/mask dimensions must match; image must contain finite values")
  }
  if (any(!is.finite(labels)) || any(labels < 0) || any(labels != floor(labels)) || max(labels) == 0) {
    stop("Mask must have positive integer object IDs and background 0")
  }
  if (max(labels) > 1000000) stop("Unexpectedly large label IDs; export a compact integer label image")
  storage.mode(labels) <- "integer"
  if (length(unique(labels[labels > 0])) == 1L) {
    message("WARNING: mask contains only one object. Verify it is a label image, not a binary mask.")
  }
  rois <- read_rois(args$rois)
  python_versions <- NULL
  if (!is.null(args$points_tsv)) {
    supplied <- utils::read.delim(args$points_tsv, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("y_px", "x_px", "score") %in% names(supplied))) {
      stop("points-tsv requires y_px, x_px and score columns")
    }
    points <- cbind(as.numeric(supplied$y_px), as.numeric(supplied$x_px))
    scores <- as.numeric(supplied$score)
  } else {
    require_package("reticulate")
    model_module <- tryCatch(reticulate::import("spotiflow.model", convert = FALSE),
      error = function(e) stop("Cannot load Python Spotiflow. Set RETICULATE_PYTHON to a compatible ",
        "Python with spotiflow installed, or supply --points-tsv. Details: ", conditionMessage(e)))
    np <- reticulate::import("numpy", convert = FALSE)
    model <- model_module$Spotiflow$from_pretrained("general")
    prediction <- model$predict(np$array(raw, dtype = "float32"), prob_thresh = args$prob_thresh)
    points <- reticulate::py_to_r(prediction[[0L]])
    points <- matrix(points, ncol = 2L)
    scores <- as.numeric(reticulate::py_to_r(prediction[[1L]]$prob))
    py_metadata <- reticulate::import("importlib.metadata")
    pkgs <- c("spotiflow", "numpy", "torch")
    python_versions <- setNames(lapply(pkgs, function(p) py_metadata$version(p)), pkgs)
  }
  if (length(scores) != nrow(points) || any(!is.finite(points)) || any(!is.finite(scores))) {
    stop("Invalid detection output: coordinates and scores must be finite and have matching lengths")
  }
  if (nrow(points) && (any(points < 0) || any(points[, 1] >= nrow(raw)) || any(points[, 2] >= ncol(raw)))) {
    stop("Spot coordinates outside image")
  }
  result <- analyze(raw, labels, points, scores, rois, args)
  dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(args$out)) stop("Could not create output directory")
  write_table(file.path(args$out, "per2_spots.txt"), result$spots)
  write_table(file.path(args$out, "cells.txt"), result$cells)
  vertices <- do.call(rbind, lapply(rois, function(roi) {
    data.frame(roi_id = roi$name, vertex_order = seq_len(nrow(roi$vertices))-1L,
               x_px = roi$vertices[, 1], y_px = roi$vertices[, 2])
  }))
  write_table(file.path(args$out, "contours.txt"), vertices)
  ijtiff::write_tif(result$assigned, file.path(args$out, "assignment_labels.tif"),
                    bits_per_sample = 32L, msg = FALSE)
  write_overlay(file.path(args$out, "qc_overlay.png"), raw, result$assigned, points, result$spots, rois)
  pkgs <- c("ijtiff", "png", "jpeg", "jsonlite", "reticulate")
  versions <- setNames(lapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NULL
  }), pkgs)
  versions$R <- R.version.string
  metadata <- list(parameters = args, shape_yx = unname(dim(raw)), versions = versions,
    python_versions = python_versions, implementation = "R",
    model = if (is.null(args$points_tsv)) "general" else "external_coordinates",
    coordinate_convention = "0-based; x rightward, y downward; um relative to image origin",
    intensity_definition = paste("Raw 2D input; aperture sum minus median annulus background times",
                                 "aperture pixel count. No deblending."),
    notes = paste("Default aperture sizes are initial values requiring validation.",
                  "Projection measurements are 2D only. Expanded nuclear territories are estimated.",
                  "Inspect boundary/overlap flags."))
  jsonlite::write_json(metadata, file.path(args$out, "run_info.json"),
                       pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
  cat(sprintf("Wrote %d spots and %d objects to %s\n", nrow(result$spots), nrow(result$cells),
              normalizePath(args$out)))
  invisible(result)
}

if (sys.nframe() == 0L) {
  tryCatch(main(), error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(status = 1L)
  })
}
