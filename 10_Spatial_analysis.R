# SCN spatial and spatiotemporal analysis: Figures 3, 4, S10, and S11 ----
# Combined sources: 10_Spatially_Variable_genes_v2.R, 11_Spatiotemporal_genes.R,
# and check_spatioltemporal_gene_only.R. Plotting helpers previously loaded from
# functions.R are included here with explicit data arguments.
# Retain figure dependencies, final gene-set definitions, selected example genes,
# and the SCN gene-expression master workbook. Omit exploratory queries, diagnostic
# plots, superseded figure versions, and whole-session saves.
# Run: Rscript 10_Spatial_analysis.R
# Comments, progress messages, and error messages are in English.

# 1. Paths and settings ----
data_root <- path.expand("~/Dropbox/Spatial_transcriptomics_new")
output_dir <- path.expand("~/Downloads")
# Save figures and SCN_gene_expr_master_table.xlsx to the same output directory.
cache_dir <- output_dir
input_files <- list(
  seurat = file.path(data_root, "sc/segmented",
                     "merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda"),
  coordinates = file.path(data_root, "sc/segmented/merged-sc-xy-coordinates.rda"),
  polygons = file.path(data_root, "rda/Brain_region_annotation/SCN_SON_PVN"),
  temporal_de = file.path(cache_dir, "DEG_dfs.rda"),
  core_shell_de = file.path(cache_dir, "cauchy_scn_core_vs_shell.Rdata"),
  scrna = file.path(cache_dir, "scRNA_SCN_line_plot_dat.rda")
)
cache_files <- list(
  expression = file.path(cache_dir, "SCN_spatial_analysis_data.rds"),
  sparkx = file.path(cache_dir, "sparkx_res_list.rda"),
  combined_sparkx = file.path(cache_dir, "df_combinedPval_all_samples_sparkx.Rdata"),
  interaction = file.path(data_root, "spatiotemporal/spatiotemp.results_2.rda"),
  master = file.path(cache_dir, "df_spatiotemporal_all_2.rda")
)
# Reuse available results by default; missing caches are rebuilt from their inputs.
# Set the appropriate flag after changing data, annotations, or statistical settings.
# Rebuilding the master table requires the existing temporal-DE and core/shell-DE
# result files. scRNA validation requires its mean, SE, and per-cell expression
# matrices in scRNA_SCN_line_plot_dat.rda; set include_scrna_figures to FALSE when
# intentionally generating only Xenium panels. None of these upstream analyses is
# rerun implicitly by this script.
rebuild_expression <- FALSE
recompute_sparkx <- FALSE
recompute_interaction <- FALSE
rebuild_master <- FALSE
include_scrna_figures <- TRUE
sparkx_cores <- 4L
alpha <- 0.05
min_expression_percent <- 5
rbf_centers <- c(5L, 10L, 15L, 20L, 25L)
rbf_bandwidths <- c(0.5, 1, 1.5, 2, 2.5)

sample_ids <- c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3")
time_samples <- list(`CT-8` = sample_ids[1:3], `CT-20` = sample_ids[4:6])
sparkx_examples <- list(CT8 = c("Prkaa2", "Nrarp"), CT20 = c("Chst8", "Rab3b"),
                        Both = c("Vip", "Gpr101"))
feature_examples <- c("Csnk1e", "Calb1", "Slc39a6", "Gad1", "Per2", "Lhx1", "Rgs16", "Grin2c")
showcase_genes <- c("Slc39a6", "Lhx1", "Rgs16", "Rpa1", "Grin2c")
# Name the representative slices explicitly instead of selecting patchwork indices.
stvg_examples <- list(
  Map1b = c("CT-8-3", "CT-20-1"), Myh9 = c("CT-8-2", "CT-20-1"),
  Cacna1g = c("CT-8-1", "CT-20-3"), Grm8 = c("CT-8-3", "CT-20-3"),
  Rbms1 = c("CT-8-2", "CT-20-2"), Gigyf2 = c("CT-8-1", "CT-20-2")
)
scrna_windows <- list(Csnk1e = c(14, 34), Calb1 = c(14, 34), Slc39a6 = c(38, 58),
                      Lhx1 = c(14, 34), Rgs16 = c(38, 58), Grin2c = c(14, 34),
                      Gad1 = c(26, 46), Per2 = c(38, 58))

required_packages <- c("tidyverse", "Seurat", "Matrix", "sp", "patchwork",
                        "ggrepel", "ggvenn", "UpSetR", "openxlsx")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                             logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Please install: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
})

# 2. Input validation, normalization, and SCN extraction ----
load_rda <- function(path, object_names) {
  if (!file.exists(path)) stop("Required input not found: ", path)
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  missing <- setdiff(object_names, loaded)
  if (length(missing)) stop(path, " is missing objects: ", paste(missing, collapse = ", "))
  mget(object_names, envir = env, inherits = FALSE)
}

require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, colnames(data))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "))
}

check_gene_table <- function(data, label) {
  require_columns(data, "gene", label)
  if (anyNA(data$gene) || anyDuplicated(data$gene)) stop(label, " contains missing or duplicate gene IDs.")
  data
}

# Match the original library-size normalization: each column sums to the median
# library size. This is neither CPM nor log normalization. Preserve sparsity.
normalize_depth <- function(counts) {
  totals <- Matrix::colSums(counts)
  if (any(!is.finite(totals) | totals <= 0)) stop("Cannot normalize a column with zero or invalid total counts.")
  normalized <- counts %*% Matrix::Diagonal(x = median(totals) / totals)
  dimnames(normalized) <- dimnames(counts)
  normalized
}

read_polygons <- function(directory) {
  files <- list.files(directory, pattern = "_regions_v5\\.csv$", full.names = TRUE)
  names(files) <- gsub("_", "-", sub("_regions.*$", "", basename(files)))
  if (anyDuplicated(names(files)) || !all(sample_ids %in% names(files))) {
    stop("Exactly one v5 polygon file is required for each of the six samples: ", directory)
  }
  lapply(files[sample_ids], readr::read_csv, skip = 2, show_col_types = FALSE)
}

# Anchor the region pattern so that SCN does not also select SCN_core or SCN_shell.
# Include polygon boundaries (> 0), and deduplicate cells in overlapping polygons.
select_region_cells <- function(coordinates, polygons, region) {
  setNames(lapply(sample_ids, function(id) {
    polygon <- polygons[[id]]
    require_columns(polygon, c("Selection", "X", "Y"), paste(id, "polygons"))
    selections <- unique(polygon$Selection[grepl(paste0("^", region, "_[0-9]+$"), polygon$Selection)])
    if (!length(selections)) stop(id, " has no polygons for ", region)
    xy <- coordinates[[id]]
    unique(unlist(lapply(selections, function(selection) {
      vertices <- polygon[polygon$Selection == selection, c("X", "Y")]
      inside <- sp::point.in.polygon(xy$x, xy$y, vertices$X, vertices$Y) > 0
      xy$cell[inside]
    }), use.names = FALSE))
  }), sample_ids)
}

aggregate_counts <- function(counts, cells_by_sample) {
  result <- vapply(sample_ids, function(id) {
    cells <- intersect(cells_by_sample[[id]], colnames(counts))
    if (!length(cells)) stop("No cells are available for pseudobulk aggregation in ", id)
    Matrix::rowSums(counts[, cells, drop = FALSE])
  }, numeric(nrow(counts)))
  dimnames(result) <- list(rownames(counts), sample_ids)
  result
}

prepare_expression_data <- function() {
  object <- load_rda(input_files$seurat, "xenium_obj_filtered")$xenium_obj_filtered
  coordinates <- load_rda(input_files$coordinates, "xy_coordinates")$xy_coordinates
  if (!all(sample_ids %in% names(coordinates))) stop("Coordinates are missing required samples.")
  coordinates <- lapply(coordinates[sample_ids], function(xy) {
    require_columns(xy, c("cell", "x", "y"), "Coordinates")
    xy[xy$cell %in% colnames(object), c("cell", "x", "y"), drop = FALSE]
  })
  polygons <- read_polygons(input_files$polygons)
  regions <- setNames(lapply(c("SCN", "SCN_core", "SCN_shell"), function(region) {
    select_region_cells(coordinates, polygons, region)
  }), c("SCN", "SCN_core", "SCN_shell"))
  cells <- unique(unlist(regions$SCN, use.names = FALSE))
  if (!length(cells)) stop("No SCN cells were selected.")
  object <- subset(object, cells = cells)
  counts <- GetAssayData(object, assay = "Xenium", layer = "counts")
  metadata <- object[[]][colnames(counts), , drop = FALSE]
  require_columns(metadata, c("orig.ident", "cluster", "predicted.id"), "Seurat metadata")
  metadata$celltype <- with(metadata, case_when(
    cluster == "5" ~ "Oligodendrocytes", cluster == "3" ~ "Astrocytes",
    cluster == "8" ~ "Ependymal", TRUE ~ as.character(predicted.id)
  ))
  coordinates <- lapply(coordinates, function(xy) xy[xy$cell %in% colnames(counts), , drop = FALSE])
  xy <- bind_rows(coordinates) %>% tibble::column_to_rownames("cell")
  if (anyNA(xy[colnames(counts), c("x", "y")])) stop("Selected SCN cells have missing coordinates.")
  pseudobulk <- aggregate_counts(counts, regions$SCN)
  core <- aggregate_counts(counts, regions$SCN_core)
  shell <- aggregate_counts(counts, regions$SCN_shell)
  colnames(core) <- paste0("SCN_core_", sample_ids)
  colnames(shell) <- paste0("SCN_shell_", sample_ids)
  # Use one median library size across all 12 core/shell pseudobulk columns.
  regional <- cbind(core, shell)
  result <- list(version = 1L, counts = counts, normalized = normalize_depth(counts),
                 metadata = metadata, coordinates = coordinates, polygons = polygons,
                 pseudobulk = normalize_depth(pseudobulk),
                 regional = normalize_depth(regional))
  saveRDS(result, cache_files$expression)
  result
}

# 3. SPARK-X and the adjusted spatial-by-time interaction model ----
# Preserve the equal-weight Cauchy combination used by both source scripts,
# including their treatment of p = 1. Missing tests remain missing until table assembly.
cauchy_combine <- function(p) {
  if (!length(p) || anyNA(p)) return(NA_real_)
  if (any(!is.finite(p) | p < 0 | p > 1)) stop("P-values must be finite and lie between 0 and 1.")
  if (any(p == 1)) p[p == 1] <- max(0.9999, 1 - 1 / length(p))
  min(pcauchy(sum(tan((0.5 - p) * pi) / length(p)), lower.tail = FALSE), 1)
}

run_sparkx <- function(data) {
  if (!requireNamespace("SPARK", quietly = TRUE)) stop("Install SPARK to recompute SPARK-X results.")
  setNames(lapply(sample_ids, function(id) {
    message("Running SPARK-X: ", id)
    xy <- data$coordinates[[id]]
    locus <- as.matrix(xy[, c("x", "y"), drop = FALSE])
    rownames(locus) <- xy$cell
    SPARK::sparkx(count_in = data$counts[, xy$cell, drop = FALSE],
                  locus_in = locus, numCores = sparkx_cores)
  }), sample_ids)
}

combine_sparkx <- function(results) {
  if (!all(sample_ids %in% names(results))) stop("SPARK-X results are missing required samples.")
  tables <- setNames(lapply(sample_ids, function(id) {
    table <- as.data.frame(results[[id]]$res_mtest)
    require_columns(table, c("combinedPval", "adjustedPval"), paste(id, "SPARK-X"))
    table$gene <- rownames(table)
    check_gene_table(table, paste(id, "SPARK-X"))
  }), sample_ids)
  per_sample <- purrr::imap(tables, function(table, id) {
    rename_with(table, ~ paste0(id, "_", .x), -gene)
  }) %>% purrr::reduce(inner_join, by = "gene")
  per_time <- purrr::imap(time_samples, function(ids, timepoint) {
    genes <- Reduce(intersect, lapply(tables[ids], function(x) x$gene))
    p <- vapply(ids, function(id) tables[[id]]$combinedPval[match(genes, tables[[id]]$gene)], numeric(length(genes)))
    combined <- apply(p, 1, cauchy_combine)
    table <- data.frame(gene = genes, cauchy_pval = combined,
                         cauchy_padj = p.adjust(combined, method = "BH"))
    rename_with(table, ~ paste0(timepoint, "_", .x), -gene)
  }) %>% purrr::reduce(inner_join, by = "gene")
  inner_join(per_sample, per_time, by = "gene")
}

# Standardize each coordinate axis within each slice to [-0.5, 0.5]. A constant
# axis is placed at 0, matching the original centered-coordinate convention.
standardize_coordinates <- function(coordinates, samples) {
  coordinates <- as.matrix(coordinates)
  if (nrow(coordinates) != length(samples) || any(!is.finite(coordinates))) {
    stop("Coordinates and sample labels must align and contain finite values.")
  }
  standardized <- coordinates
  for (id in unique(samples)) {
    rows <- which(samples == id)
    for (axis in seq_len(ncol(coordinates))) {
      values <- coordinates[rows, axis]
      span <- diff(range(values))
      standardized[rows, axis] <- if (span == 0) 0 else (values - min(values)) / span - 0.5
    }
  }
  standardized
}

gaussian_basis <- function(coordinates, centers, length_scale) {
  if (!is.finite(length_scale) || length_scale <= 0) stop("The RBF length scale must be positive.")
  vapply(seq_len(nrow(centers)), function(k) {
    delta <- sweep(coordinates, 2, centers[k, ], "-")
    exp(-rowSums(delta^2) / (2 * length_scale^2))
  }, numeric(nrow(coordinates)))
}

# Compare y ~ time + spatial basis + covariates with the same model plus
# spatial-basis-by-time interactions. Slice and cell-type indicator matrices are
# retained exactly as in the source analysis. For compatibility, use its nominal
# column-count degrees of freedom, including redundant indicator columns, rather
# than silently changing inference to rank-based degrees of freedom during cleanup.
# Cache QR decompositions to avoid refactoring the same design for every gene.
interaction_pvalues <- function(expression, basis, time, covariates, min_var = 1e-8) {
  reduced <- cbind(1, time, basis, covariates)
  full <- cbind(1, time, basis, basis * time, covariates)
  qr_reduced <- qr(reduced, tol = 1e-7)
  qr_full <- qr(full, tol = 1e-7)
  p <- vapply(seq_len(nrow(expression)), function(g) {
    y <- as.numeric(expression[g, ])
    keep <- is.finite(y)
    if (sum(keep) < 2 || var(y[keep]) <= min_var) return(NA_real_)
    if (all(keep)) {
      residual_r <- qr.resid(qr_reduced, y)
      residual_f <- qr.resid(qr_full, y)
    } else {
      residual_r <- lm.fit(reduced[keep, , drop = FALSE], y[keep])$residuals
      residual_f <- lm.fit(full[keep, , drop = FALSE], y[keep])$residuals
    }
    df_full <- sum(keep) - ncol(full)
    df_difference <- ncol(full) - ncol(reduced)
    if (df_full <= 0 || df_difference <= 0) return(NA_real_)
    rss_r <- sum(residual_r^2)
    rss_f <- sum(residual_f^2)
    statistic <- ((rss_r - rss_f) / df_difference) / (rss_f / df_full)
    pf(statistic, df_difference, df_full, lower.tail = FALSE)
  }, numeric(1))
  setNames(p, rownames(expression))
}

run_interaction_model <- function(data) {
  cells <- colnames(data$normalized)
  xy <- bind_rows(data$coordinates) %>% tibble::column_to_rownames("cell")
  coordinates <- standardize_coordinates(xy[cells, c("x", "y")], data$metadata$orig.ident)
  time <- as.numeric(data$metadata$orig.ident %in% time_samples[["CT-20"]])
  slice_id <- factor(data$metadata$orig.ident)
  cell_type <- factor(data$metadata$celltype)
  covariates <- cbind(model.matrix(~ slice_id - 1), model.matrix(~ cell_type - 1))
  expression <- data$normalized[Matrix::rowSums(data$counts) > 0, , drop = FALSE]
  grid_p <- matrix(NA_real_, nrow(expression), length(rbf_centers) * length(rbf_bandwidths),
                   dimnames = list(rownames(expression), NULL))
  i <- 0L
  for (n_centers in rbf_centers) {
    set.seed(1)
    centers <- kmeans(coordinates, centers = n_centers, nstart = 25, iter.max = 200)$centers
    # Only distances from cells to centers are needed; avoid an n_cells^2 matrix.
    distances <- vapply(seq_len(nrow(centers)), function(k) {
      sqrt(rowSums(sweep(coordinates, 2, centers[k, ], "-")^2))
    }, numeric(nrow(coordinates)))
    base_scale <- median(apply(distances, 1, min))
    for (bandwidth in rbf_bandwidths) {
      i <- i + 1L
      message("Interaction model: ", n_centers, " centers, bandwidth ", bandwidth)
      basis <- gaussian_basis(coordinates, centers, bandwidth * base_scale)
      grid_p[, i] <- interaction_pvalues(expression, basis, time, covariates)
    }
  }
  p <- apply(grid_p, 1, cauchy_combine)
  data.frame(spatiotemp.pval = p, spatiotemp.padj = p.adjust(p, method = "BH"),
             row.names = rownames(expression))
}

# Build the same gene universe used by the final analysis: interaction-tested genes
# with temporal-DE results. Recalculate BH adjustments after assembling this universe.
# Fill missing p-values with 1; leave missing effect sizes as NA instead of inventing
# log fold changes. Read named DE columns rather than relying on column positions.
build_master_table <- function(sparkx, core_shell, interaction, temporal_de) {
  sparkx <- check_gene_table(sparkx, "Combined SPARK-X results")
  core_shell <- check_gene_table(core_shell, "Core/shell DE results")
  temporal_de <- check_gene_table(temporal_de, "Temporal DE results")
  if (!"gene" %in% names(interaction)) interaction <- tibble::rownames_to_column(interaction, "gene")
  interaction <- check_gene_table(interaction, "Interaction results")
  require_columns(sparkx, c("CT-8_cauchy_pval", "CT-20_cauchy_pval"), "SPARK-X results")
  require_columns(core_shell, c("CT8_cauchy_pval", "CT20_cauchy_pval"), "Core/shell results")
  require_columns(interaction, c("spatiotemp.pval", "spatiotemp.padj"), "Interaction results")
  require_columns(temporal_de, c("p_val", "logFC"), "Temporal DE results")
  for (timepoint in names(time_samples)) {
    target <- paste0(timepoint, "_avg_log2FC")
    sample_columns <- paste0(time_samples[[timepoint]], "_avg_log2FC")
    if (all(sample_columns %in% names(core_shell))) {
      core_shell[[target]] <- rowMeans(core_shell[, sample_columns, drop = FALSE])
    } else {
      require_columns(core_shell, target, "Core/shell results")
    }
  }
  add_suffix <- function(table, suffix) rename_with(table, ~ paste0(.x, suffix), -gene)
  tables <- list(add_suffix(sparkx, "_sparkx"), add_suffix(core_shell, "_core_vs_shell"),
                 add_suffix(interaction, "_spatiotemporal"))
  master <- purrr::reduce(tables, full_join, by = "gene") %>%
    dplyr::filter(gene %in% interaction$gene) %>%
    inner_join(temporal_de %>% dplyr::select(gene, p_val, logFC) %>%
                 rename(p_val_TVG = p_val, logFC_TVG = logFC), by = "gene")
  p_columns <- grep("pval|padj|Pval|p_val|p_adj", names(master), value = TRUE)
  for (column in p_columns) master[[column]][is.na(master[[column]])] <- 1
  adjustments <- c(`CT-8_cauchy_padj_sparkx` = "CT-8_cauchy_pval_sparkx",
                    `CT-20_cauchy_padj_sparkx` = "CT-20_cauchy_pval_sparkx",
                    CT8_cauchy_padj_core_vs_shell = "CT8_cauchy_pval_core_vs_shell",
                    CT20_cauchy_padj_core_vs_shell = "CT20_cauchy_pval_core_vs_shell",
                    p_adj_TVG = "p_val_TVG")
  for (column in names(adjustments)) master[[column]] <- p.adjust(master[[adjustments[[column]]]], method = "BH")
  master
}

# These are the final definitions from check_spatioltemporal_gene_only.R:
# SVG and STVG require expression in > 5% of all SCN cells; TVG uses only FDR < 0.05.
# Do not reintroduce the older log2(1.1) TVG cutoff or a prevalence filter for TVGs.
classify_genes <- function(master, data) {
  required <- c("gene", "CT-8_cauchy_padj_sparkx", "CT-20_cauchy_padj_sparkx", "p_adj_TVG",
                 "spatiotemp.padj_spatiotemporal", "CT8_cauchy_padj_core_vs_shell",
                 "CT20_cauchy_padj_core_vs_shell", "CT-8_avg_log2FC_core_vs_shell",
                 "CT-20_avg_log2FC_core_vs_shell", "logFC_TVG")
  require_columns(master, required, "Gene master table")
  check_gene_table(master, "Gene master table")
  percent <- Matrix::rowMeans(data$normalized != 0) * 100
  expressed <- names(percent)[percent > min_expression_percent]
  significant <- function(column) master$gene[which(master[[column]] < alpha)]
  sets <- list(CT8_SVG = intersect(significant("CT-8_cauchy_padj_sparkx"), expressed),
                CT20_SVG = intersect(significant("CT-20_cauchy_padj_sparkx"), expressed),
                TVG = significant("p_adj_TVG"),
                STVG = intersect(significant("spatiotemp.padj_spatiotemporal"), expressed))
  split <- list(CT8_only = setdiff(sets$CT8_SVG, sets$CT20_SVG),
                 CT20_only = setdiff(sets$CT20_SVG, sets$CT8_SVG),
                 Both = intersect(sets$CT8_SVG, sets$CT20_SVG))
  groups <- master
  for (timepoint in c("CT8", "CT20")) {
    p <- groups[[paste0(timepoint, "_cauchy_padj_core_vs_shell")]]
    fc <- groups[[paste0(ifelse(timepoint == "CT8", "CT-8", "CT-20"), "_avg_log2FC_core_vs_shell")]]
    groups[[paste0(timepoint, "_group")]] <- case_when(
      p < alpha & fc > 0 ~ "Core", p < alpha & fc < 0 ~ "Shell", TRUE ~ "None"
    )
  }
  list(sets = sets, split = split, groups = groups, percent = percent)
}

# Export the original seven-sheet master workbook using the same gene order on
# every sheet. Include all master-table genes, not only significant or > 5%-expressed
# genes. Statistical columns retain their original names and numeric values.
build_master_workbook_tables <- function(data, master) {
  check_gene_table(master, "Gene master table")
  column_map <- c(
    gene = "gene",
    SVG_CT8_pval = "CT-8_cauchy_pval_sparkx", SVG_CT8_padj = "CT-8_cauchy_padj_sparkx",
    SVG_CT20_pval = "CT-20_cauchy_pval_sparkx", SVG_CT20_padj = "CT-20_cauchy_padj_sparkx",
    TVG_pval = "p_val_TVG", TVG_padj = "p_adj_TVG", TVG_logFC = "logFC_TVG",
    core_vs_shell_CT8_pval = "CT8_cauchy_pval_core_vs_shell",
    core_vs_shell_CT8_padj = "CT8_cauchy_padj_core_vs_shell",
    core_vs_shell_CT8_logFC = "CT-8_avg_log2FC_core_vs_shell",
    core_vs_shell_CT20_pval = "CT20_cauchy_pval_core_vs_shell",
    core_vs_shell_CT20_padj = "CT20_cauchy_padj_core_vs_shell",
    core_vs_shell_CT20_logFC = "CT-20_avg_log2FC_core_vs_shell",
    spatiotemp_pval = "spatiotemp.pval_spatiotemporal",
    spatiotemp_padj = "spatiotemp.padj_spatiotemporal"
  )
  require_columns(master, unname(column_map), "Gene master table for Excel export")
  require_columns(data$metadata, c("orig.ident", "celltype"), "SCN metadata for Excel export")
  if (!identical(colnames(data$counts), rownames(data$metadata))) {
    stop("Raw counts and metadata cell orders differ for Excel export.")
  }
  genes <- as.character(master$gene)
  if (!length(genes)) stop("The gene master table is empty.")
  if (!all(genes %in% rownames(data$counts)) || !all(genes %in% rownames(data$regional))) {
    stop("Master-table genes are missing from the expression cache; rebuild the cache or check the input files.")
  }
  statistics <- as.data.frame(master[, unname(column_map), drop = FALSE])
  names(statistics) <- names(column_map)

  # Match the original workbook: subset to master-table genes first, then normalize
  # each pseudobulk sample to 1e6. This CPM export is separate from the median-depth
  # normalization used in figures. Samples with no counts remain entirely zero.
  cpm_table <- function(matrix) {
    matrix <- as.matrix(matrix[genes, , drop = FALSE])
    totals <- colSums(matrix)
    if (any(!is.finite(totals) | totals < 0)) stop("Invalid expression totals in Excel export.")
    nonzero <- totals > 0
    if (any(nonzero)) {
      matrix[, nonzero] <- sweep(matrix[, nonzero, drop = FALSE], 2, totals[nonzero], "/") * 1e6
    }
    data.frame(gene = genes, matrix, check.names = FALSE, row.names = NULL)
  }
  tables <- list(SVG_TVG_STG_table = statistics)
  celltypes <- c(SCN_GABA_expr = "GABAergic", SCN_Astrocytes_expr = "Astrocytes",
                 SCN_Vascular_expr = "Vascular", SCN_Oligo_expr = "Oligodendrocytes")
  for (sheet in names(celltypes)) {
    # Aggregate raw counts directly instead of loading a second muscat result cache.
    # Empty cell-type/sample combinations produce zero columns rather than NaN.
    counts <- vapply(sample_ids, function(id) {
      cells <- rownames(data$metadata)[which(data$metadata$orig.ident == id &
                                             data$metadata$celltype == celltypes[[sheet]])]
      Matrix::rowSums(data$counts[genes, cells, drop = FALSE])
    }, numeric(length(genes)))
    counts <- matrix(counts, nrow = length(genes), dimnames = list(genes, sample_ids))
    tables[[sheet]] <- cpm_table(counts)
  }
  for (region in c("core", "shell")) {
    columns <- paste0("SCN_", region, "_", sample_ids)
    require_columns(data$regional, columns, "Core/shell expression cache")
    regional <- data$regional[, columns, drop = FALSE]
    colnames(regional) <- sample_ids
    # Re-normalizing each cached regional column to CPM cancels its stored
    # median-depth scaling and gives the same result as using raw pseudobulk counts.
    tables[[paste0("SCN_", region, "_expr")]] <- cpm_table(regional)
  }
  tables
}

save_master_workbook <- function(data, master,
                                 path = file.path(output_dir, "SCN_gene_expr_master_table.xlsx")) {
  tables <- build_master_workbook_tables(data, master)
  workbook <- openxlsx::createWorkbook()
  scientific <- openxlsx::createStyle(numFmt = "0.000E+00")
  decimal <- openxlsx::createStyle(numFmt = "0.0000")
  for (sheet in names(tables)) {
    table <- tables[[sheet]]
    openxlsx::addWorksheet(workbook, sheet, gridLines = FALSE)
    openxlsx::writeDataTable(workbook, sheet, table, tableStyle = "TableStyleMedium2",
                             rowNames = FALSE, keepNA = FALSE)
    openxlsx::freezePane(workbook, sheet, firstRow = TRUE, firstCol = TRUE)
    widths <- pmax(14, pmin(34, nchar(names(table)) + 2))
    widths[1] <- 24
    openxlsx::setColWidths(workbook, sheet, cols = seq_along(table), widths = widths)
    rows <- seq_len(nrow(table)) + 1L
    numeric_columns <- which(vapply(table, is.numeric, logical(1)))
    p_columns <- grep("_(pval|padj)$", names(table))
    if (length(numeric_columns)) {
      openxlsx::addStyle(workbook, sheet, decimal, rows = rows, cols = numeric_columns, gridExpand = TRUE)
    }
    if (length(p_columns)) {
      openxlsx::addStyle(workbook, sheet, scientific, rows = rows, cols = p_columns, gridExpand = TRUE)
    }
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(workbook, path, overwrite = TRUE)
  message("SCN gene-expression master workbook has been saved to: ", path)
  invisible(path)
}

# 4. Shared figure functions ----
save_figure <- function(plot, filename, width, height) {
  ggsave(file.path(output_dir, filename), plot = plot, width = width, height = height, units = "in")
}

save_pages <- function(plots, filename, width, height) {
  if (!length(plots)) return(invisible(NULL))
  grDevices::pdf(file.path(output_dir, filename), width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (plot in plots) print(plot)
  invisible(NULL)
}

save_upset <- function(sets) {
  # Construct and print within the PDF device so Rscript does not create Rplots.pdf.
  grDevices::pdf(file.path(output_dir, "upset_plot.pdf"), width = 7.65, height = 4.86)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (!length(unique(unlist(sets, use.names = FALSE)))) {
    print(ggplot() + annotate("text", x = 0, y = 0, label = "No significant genes") + theme_void())
  } else {
    print(UpSetR::upset(UpSetR::fromList(sets), sets = rev(c("TVG", "CT8_SVG", "CT20_SVG", "STVG")),
                         keep.order = TRUE, order.by = "freq"))
  }
}

feature_panels <- function(data, gene, samples = sample_ids, sparkx_labels = NULL) {
  if (!gene %in% rownames(data$normalized)) stop("Expression data are missing gene: ", gene)
  expression <- as.numeric(data$normalized[gene, ])
  names(expression) <- colnames(data$normalized)
  upper <- as.numeric(quantile(expression, 0.99, na.rm = TRUE))
  setNames(lapply(samples, function(id) {
    xy <- data$coordinates[[id]]
    df <- mutate(xy, normalized_expr = unname(expression[cell])) %>% arrange(normalized_expr)
    # Feature plots share the gene's 99th-percentile scale across all six slices.
    # SPARK-X example panels retain the original per-slice percentile scaling.
    limit <- if (is.null(sparkx_labels)) upper else as.numeric(quantile(df$normalized_expr, 0.99))
    if (limit <= 0) limit <- 1
    df$point_expr <- if (is.null(sparkx_labels)) pmin(df$normalized_expr, limit) else df$normalized_expr
    title <- id
    if (!is.null(sparkx_labels)) {
      column <- paste0(id, "_adjustedPval_sparkx")
      require_columns(sparkx_labels, column, "SPARK-X master table for example labels")
      p <- sparkx_labels[[column]][match(gene, sparkx_labels$gene)]
      title <- sprintf("%s, %s: p.adj=%.2e", id, gene, p)
    }
    plot <- ggplot(df, aes(x, y, color = normalized_expr, size = point_expr)) +
      geom_point(alpha = if (is.null(sparkx_labels)) 1 else 0.8) +
      scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, limit), oob = scales::squish) +
      scale_size_continuous(range = if (is.null(sparkx_labels)) c(1, 4) else c(2, 5)) +
      theme_void() + labs(title = title, color = "Normalized\nExpression")
    if (is.null(sparkx_labels)) {
      polygon <- data$polygons[[id]] %>% dplyr::filter(grepl("^SCN_(core|shell)_[0-9]+$", Selection))
      plot <- plot + geom_polygon(data = polygon, aes(X, Y, group = Selection),
                                  fill = NA, color = "black", linewidth = 0.3, inherit.aes = FALSE)
    } else {
      plot <- plot + theme(legend.position = "none")
    }
    plot
  }), samples)
}

regional_expression_data <- function(data, gene) {
  if (!gene %in% rownames(data$regional)) stop("Core/shell data are missing gene: ", gene)
  data.frame(group = colnames(data$regional), norm_expr = as.numeric(data$regional[gene, ])) %>%
    mutate(region = sub("_CT.*$", "", group), region = sub("^SCN_", "", region),
           CT = if_else(grepl("_CT-8-", group), "CT8", "CT20"),
           CT = factor(CT, levels = c("CT8", "CT20")),
           group_1 = factor(paste(CT, region, sep = "_"),
                            levels = c("CT8_core", "CT8_shell", "CT20_core", "CT20_shell")))
}

plot_regional_expression <- function(data, gene) {
  ggplot(regional_expression_data(data, gene), aes(group_1, norm_expr, shape = region, color = CT)) +
    geom_point(size = 4, alpha = 0.7) +
    stat_summary(aes(group = group_1), fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.4) +
    theme_classic() + labs(title = gene, x = NULL, y = "Normalized expression")
}

plot_showcase <- function(data, master, gene) {
  maps <- patchwork::wrap_plots(feature_panels(data, gene), ncol = 3, guides = "collect")
  regional <- ggplot(regional_expression_data(data, gene), aes(region, norm_expr, color = CT, shape = CT)) +
    geom_point(size = 3, alpha = 0.7) + theme_classic() + labs(title = "DEG (Core vs Shell)", y = "Normalized expression")
  temporal <- data.frame(sample = colnames(data$pseudobulk), norm_expr = as.numeric(data$pseudobulk[gene, ])) %>%
    mutate(timepoint = factor(if_else(grepl("^CT-8-", sample), "CT-8", "CT-20"), levels = c("CT-8", "CT-20"))) %>%
    ggplot(aes(timepoint, norm_expr, color = timepoint, shape = timepoint)) +
    geom_point(size = 3, alpha = 0.7) + theme_classic() + labs(title = "DEG (CT-8 vs CT-20)", y = "Normalized expression")
  stats <- master[match(gene, master$gene), , drop = FALSE]
  title <- sprintf(paste0("%s: SVG adj. P (CT8)=%.2e, SVG adj. P (CT20)=%.2e, TVG adj. P=%.2e, STVG adj. P=%.2e\n",
                          "Core vs shell adj. P: CT8=%.2e, CT20=%.2e"),
                    gene, stats[["CT-8_cauchy_padj_sparkx"]], stats[["CT-20_cauchy_padj_sparkx"]],
                    stats$p_adj_TVG, stats$spatiotemp.padj_spatiotemporal,
                    stats$CT8_cauchy_padj_core_vs_shell, stats$CT20_cauchy_padj_core_vs_shell)
  side <- patchwork::wrap_plots(regional, temporal, ncol = 1, guides = "collect")
  patchwork::wrap_plots(maps, side, ncol = 2, widths = c(3, 0.8)) + patchwork::plot_annotation(title = title)
}

# Use the final Figure 4B definition: SVG + TVG + STVG, with core/shell significance
# required at both times. The earlier absolute-log2FC > 0.2 selection is superseded.
plot_spatiotemporal_scatter <- function(classification) {
  sets <- classification$sets
  genes <- intersect(intersect(union(sets$CT8_SVG, sets$CT20_SVG), sets$TVG), sets$STVG)
  df <- classification$groups %>% dplyr::filter(gene %in% genes) %>%
    mutate(x_fc = .data[["CT-8_avg_log2FC_core_vs_shell"]], y_fc = .data[["CT-20_avg_log2FC_core_vs_shell"]],
           Sig = if_else(CT8_cauchy_padj_core_vs_shell < alpha & CT20_cauchy_padj_core_vs_shell < alpha,
                          "Sig", "Not Sig"))
  # Label up to ten significant genes from each sign quadrant, as in the final plot.
  labels <- df %>% dplyr::filter(Sig == "Sig", !is.na(x_fc), !is.na(y_fc), x_fc != 0, y_fc != 0) %>%
    group_by(x_positive = x_fc > 0, y_positive = y_fc > 0) %>%
    arrange(desc(abs(x_fc)), desc(abs(y_fc)), .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup()
  ggplot(df, aes(x_fc, y_fc, color = Sig)) + geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    ggrepel::geom_text_repel(data = labels, aes(label = gene), size = 3, box.padding = 0.4,
                            point.padding = 0.3, min.segment.length = 0, max.overlaps = Inf,
                            seed = 123, show.legend = FALSE) + theme_classic() +
    labs(x = "Core vs shell log2FC (CT8)", y = "log2FC (CT20)",
         title = "Spatiotemporal genes (SVG + TVG + STVG) enriched in core/shell")
}

plot_svg_without_tvg <- function(classification) {
  groups <- classification$groups
  split <- lapply(classification$split, setdiff, y = classification$sets$TVG)
  both <- groups %>% dplyr::filter(gene %in% split$Both, CT8_group != "None", CT20_group != "None")
  shared_example <- both %>% dplyr::filter(gene == "Slc39a6")
  single_plot <- function(timepoint) {
    other <- ifelse(timepoint == "CT8", "CT20", "CT8")
    fc <- paste0(ifelse(timepoint == "CT8", "CT-8", "CT-20"), "_avg_log2FC_core_vs_shell")
    group_column <- paste0(timepoint, "_group")
    df <- groups %>% dplyr::filter(gene %in% split[[paste0(timepoint, "_only")]],
                                   .data[[group_column]] != "None", .data[[paste0(other, "_group")]] == "None") %>%
      arrange(desc(abs(.data[[fc]]))) %>% slice_head(n = 29) %>%
      bind_rows(shared_example) %>% distinct(gene, .keep_all = TRUE) %>% arrange(.data[[fc]]) %>%
      mutate(gene = factor(gene, levels = rev(gene)))
    ggplot(df, aes(.data[[fc]], gene, fill = .data[[group_column]])) +
      geom_col() + theme_classic() + labs(title = timepoint, x = "Core vs shell log2FC", y = NULL)
  }
  both <- both %>% mutate(log2fc = (abs(.data[["CT-8_avg_log2FC_core_vs_shell"]]) +
                                    abs(.data[["CT-20_avg_log2FC_core_vs_shell"]])) / 2) %>% arrange(desc(log2fc))
  # Deduplicate switching genes that are already among the top 27.
  selected <- bind_rows(slice_head(both, n = 27), dplyr::filter(both, CT8_group != CT20_group)) %>%
    distinct(gene, .keep_all = TRUE) %>% mutate(group = paste(CT8_group, CT20_group, sep = " -> "),
                                               gene = factor(gene, levels = rev(gene)))
  p_both <- ggplot(selected, aes(log2fc, gene, fill = group)) + geom_col() + theme_classic() +
    labs(title = "Both", x = "Mean absolute log2FC", y = NULL)
  patchwork::wrap_plots(single_plot("CT8"), single_plot("CT20"), p_both, ncol = 3)
}

plot_tvg_pie <- function(svg_genes, tvg_genes, title) {
  df <- data.frame(group = factor(c("TVG", "Not_TVG"), levels = c("Not_TVG", "TVG")),
                    freq = c(sum(svg_genes %in% tvg_genes), sum(!svg_genes %in% tvg_genes))) %>%
    mutate(percent = if (sum(freq) > 0) round(100 * freq / sum(freq), 2) else 0,
           label = sprintf("%s\n(%s)", freq, percent))
  ggplot(df, aes("", freq, fill = group)) + geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") + geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    theme_void() + ggtitle(title)
}

# Independent scRNA-seq validation: retain mean +/- SE curves, original CT windows,
# and the original Welch t-test on individual cells at the selected CT6/CT18 samples.
# This comparison does not introduce an animal-level replicate model.
plot_scrna_validation <- function(scrna, gene, window) {
  extract_long <- function(table, value) {
    require_columns(table, "gene", "scRNA summary")
    table %>% dplyr::filter(.data$gene == .env$gene) %>%
      pivot_longer(-gene, names_to = "CT", values_to = value)
  }
  df <- inner_join(extract_long(scrna$scRNA_mat_norm_mean, "norm_expr"),
                    extract_long(scrna$scRNA_mat_norm_se, "se"), by = c("gene", "CT")) %>%
    mutate(CT = as.numeric(sub("^CT", "", CT))) %>%
    dplyr::filter(CT >= window[1], CT <= window[2]) %>%
    mutate(CT_cycle = ((CT - 1) %% 24) + 1)
  if (!nrow(df)) stop("No scRNA summary data for ", gene, " in its selected CT window.")
  day <- df$CT[df$CT_cycle == 6]
  night <- df$CT[df$CT_cycle == 18]
  if (length(day) != 1 || length(night) != 1) stop("Expected one CT6 and one CT18 sample for ", gene)
  sample_names <- paste0("CT", c(day, night))
  if (!all(sample_names %in% names(scrna$mat_norm_list))) stop("Missing scRNA sample matrices for ", gene)
  expression <- lapply(sample_names, function(id) {
    mat <- scrna$mat_norm_list[[id]]
    if (!gene %in% rownames(mat)) stop("Missing scRNA gene ", gene, " in ", id)
    as.numeric(mat[gene, ])
  })
  p <- t.test(expression[[1]], expression[[2]])$p.value
  night_interval <- data.frame(xmin = 12, xmax = 22)
  ggplot(df, aes(CT_cycle, norm_expr)) +
    geom_rect(data = night_interval, aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
               fill = "grey40", alpha = 0.15, inherit.aes = FALSE) +
    geom_line() + geom_ribbon(aes(ymin = norm_expr - se, ymax = norm_expr + se), alpha = 0.2) +
    geom_point() + theme_classic() + scale_x_continuous(breaks = seq(2, 22, 4)) +
    labs(title = sprintf("%s: CT6 vs CT18 pval=%.2e", gene, p), x = "CT", y = "Normalized expression")
}

# 5. Assemble inputs and results without sourcing scripts with workspace side effects ----
load_analysis_results <- function(data) {
  use_master <- file.exists(cache_files$master) &&
    !any(c(rebuild_master, recompute_sparkx, recompute_interaction, rebuild_expression))
  if (use_master) return(load_rda(cache_files$master, "df_spatiotemporal_all_2")$df_spatiotemporal_all_2)
  if (recompute_sparkx || !file.exists(cache_files$combined_sparkx)) {
    if (recompute_sparkx || !file.exists(cache_files$sparkx)) {
      sparkx_res_list <- run_sparkx(data)
      save(sparkx_res_list, file = cache_files$sparkx)
    } else {
      sparkx_res_list <- load_rda(cache_files$sparkx, "sparkx_res_list")$sparkx_res_list
    }
    df_combinedPval_all_samples_sparkx <- combine_sparkx(sparkx_res_list)
    save(df_combinedPval_all_samples_sparkx, file = cache_files$combined_sparkx)
  } else {
    df_combinedPval_all_samples_sparkx <- load_rda(
      cache_files$combined_sparkx, "df_combinedPval_all_samples_sparkx"
    )$df_combinedPval_all_samples_sparkx
  }
  if (recompute_interaction || !file.exists(cache_files$interaction)) {
    spatiotemp.results <- run_interaction_model(data)
    dir.create(dirname(cache_files$interaction), recursive = TRUE, showWarnings = FALSE)
    save(spatiotemp.results, file = cache_files$interaction)
  } else {
    spatiotemp.results <- load_rda(cache_files$interaction, "spatiotemp.results")$spatiotemp.results
  }
  temporal_de <- load_rda(input_files$temporal_de, "DEG_df_SCN")$DEG_df_SCN
  core_shell <- load_rda(input_files$core_shell_de, "cauchy_scn_core_vs_shell")$cauchy_scn_core_vs_shell
  df_spatiotemporal_all_2 <- build_master_table(df_combinedPval_all_samples_sparkx,
                                               core_shell, spatiotemp.results, temporal_de)
  save(df_spatiotemporal_all_2, file = cache_files$master)
  df_spatiotemporal_all_2
}

main <- function() {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  # These upstream DE results are inputs to the source scripts, not recomputed here.
  use_master <- file.exists(cache_files$master) &&
    !any(c(rebuild_master, recompute_sparkx, recompute_interaction, rebuild_expression))
  required_inputs <- if (use_master) character() else unlist(input_files[c("temporal_de", "core_shell_de")])
  if (include_scrna_figures) required_inputs <- c(required_inputs, input_files$scrna)
  missing <- required_inputs[!file.exists(required_inputs)]
  if (length(missing)) stop("Restore the required upstream analysis files or update input_files:\n",
                             paste(missing, collapse = "\n"))
  if (rebuild_expression || !file.exists(cache_files$expression)) {
    message("Preparing SCN expression and pseudobulk data.")
    data <- prepare_expression_data()
  } else {
    data <- readRDS(cache_files$expression)
  }
  if (!identical(colnames(data$normalized), rownames(data$metadata))) stop("Expression and metadata cell orders differ.")
  master <- load_analysis_results(data)
  classification <- classify_genes(master, data)
  sets <- classification$sets
  save_master_workbook(data, master)

  # 6. Figure 3: SVG overlap, examples, TVG overlap, and regional expression ----
  p_venn <- ggvenn::ggvenn(sets[c("CT8_SVG", "CT20_SVG")])
  save_figure(p_venn, "SVG_CT8_CT20_venn.pdf", 5.79, 3.93)
  sparkx_plots <- lapply(unlist(sparkx_examples, use.names = FALSE), function(gene) {
    patchwork::wrap_plots(feature_panels(data, gene, sparkx_labels = master), ncol = 3)
  })
  # Figure 3B/C: Vip/Gpr101 and the original CT-specific example gene lists.
  save_pages(sparkx_plots, "SPARK_X_CT8_CT20_Both_genes_v2.pdf", 9.8, 5.57)
  p_pies <- patchwork::wrap_plots(plot_tvg_pie(sets$CT8_SVG, sets$TVG, "CT8"),
                                  plot_tvg_pie(sets$CT20_SVG, sets$TVG, "CT20"),
                                  ncol = 2, guides = "collect")
  save_figure(p_pies, "SVG_TVG_pie.pdf", 8.24, 4.10) # Figure 3D
  save_figure(plot_svg_without_tvg(classification), "SVG_not_TVG.pdf", 16.53, 5.60) # Figure 3E
  save_figure(ggvenn::ggvenn(list(SVG = union(sets$CT8_SVG, sets$CT20_SVG), TVG = sets$TVG)),
               "TVG_vs_SVG.pdf", 5.79, 3.93)

  # Reuse the same feature and pseudobulk plotting functions for all named panels.
  # Gad1: Figure 3F; Per2: Figure 4C/S10A; Lhx1: Figure 4E/S10B;
  # Rgs16: Figure 4F/S10C. Other genes retain the original example exports.
  for (gene in feature_examples) {
    maps <- patchwork::wrap_plots(feature_panels(data, gene), ncol = 3, guides = "collect")
    save_figure(maps, paste0(gene, "_spatial_expression.pdf"), 14.64, 6.78)
    save_figure(plot_regional_expression(data, gene), paste0(gene, "_core_shell_expression.pdf"), 3.9, 3)
  }
  map_plots <- lapply(c("Csnk1e", "Calb1", "Slc39a6"), function(gene) {
    patchwork::wrap_plots(feature_panels(data, gene), ncol = 3, guides = "collect") +
      patchwork::plot_annotation(title = gene)
  })
  save_pages(map_plots, "SVG_map_plot_v2.pdf", 14.64, 6.78)

  # 7. Figure 4 and supplementary panels ----
  save_upset(sets) # Figure 4A: final prevalence-filtered gene sets.
  save_figure(plot_spatiotemporal_scatter(classification),
               "Spatiotemporal_TVG_and_SVG_geom_point.pdf", 11.41, 6.36) # Figure 4B
  for (gene in showcase_genes) save_figure(plot_showcase(data, master, gene), paste0(gene, ".pdf"), 21.22, 8.13)
  stvg_plots <- purrr::imap(stvg_examples, function(samples, gene) {
    maps <- patchwork::wrap_plots(feature_panels(data, gene, samples), ncol = 2, guides = "collect")
    patchwork::wrap_plots(maps, plot_regional_expression(data, gene), ncol = 2, widths = c(2, 1))
  })
  save_pages(stvg_plots, "STVG_only.pdf", 15.6, 4.80) # Includes Cacna1g, Figure S11A.
  if (include_scrna_figures) {
    scrna <- load_rda(input_files$scrna, c("scRNA_mat_norm_mean", "scRNA_mat_norm_se", "mat_norm_list"))
    for (gene in names(scrna_windows)) {
      plot <- plot_scrna_validation(scrna, gene, scrna_windows[[gene]])
      save_figure(plot, paste0("p_", gene, "_4.pdf"), 8.21, 4.27)
    }
  }
  message("Spatial-analysis figures have been saved to: ", output_dir)
  invisible(NULL)
}

main()
