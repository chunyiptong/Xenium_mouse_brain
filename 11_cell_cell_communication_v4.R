# SCN cell-cell communication: Figure 5 and Figure S12 ----
# Combined sources: 11_cell_cell_communication_v3.R and merged_cellchat_results_v4.R.
# Retain SCN filtering, GABA subtyping, CellChat inference, and figure generation.
# Remove exploratory marker/pathway/LR analyses, distance debugging, and duplicate plots.
# Use the CellChat parameters from merged_cellchat_results_v4.R throughout.
# The v3 SpatialCellChat probability and filtering workflow produces different results
# and must not be mixed with this CellChat workflow.
# Run: Rscript 11_cell_cell_communication_v4.R
# Save all PDFs and required analysis caches to output_dir.

# 1. Paths, analysis parameters, and packages ----
data_root <- path.expand("~/Dropbox/Spatial_transcriptomics_new")
output_dir <- path.expand("~/Downloads")
input_files <- list(
  seurat = file.path(data_root, "sc/segmented",
                     "merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda"),
  coordinates = file.path(data_root, "sc/segmented/merged-sc-xy-coordinates.rda"),
  polygons = file.path(data_root, "rda/Brain_region_annotation/SCN_SON_PVN")
)
# Use original coordinates for both SCN polygon filtering and CellChat.
# Aligned coordinates were not used to calculate the original figures.
cache_files <- list(
  gaba = file.path(output_dir, "GABA_subclass.rda"),
  integrated = file.path(output_dir, "xenium_obj_filter_GABA_integrated.rda"),
  cellchat = file.path(output_dir, "merged_cellchat_results_v4.RData")
)
rebuild_gaba_subtypes <- FALSE # Also rebuild automatically if either GABA cache is missing.
recompute_cellchat <- TRUE    # Set to FALSE to load the existing v4 CellChat cache.
random_seed <- 1234L
future_workers <- 8L
pval_cutoff <- 0.01
test_alternative <- "greater"
# Preserve the original stat_compare_means group order: the one-sided test is CT8 > CT20.
# Pair observations by source-target cell-type combination; these are not independent
# replicates from the three animals or sections.

sample_ids <- c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3")
celltype_levels <- c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA",
                     "Astrocyte", "Oligodendrocyte", "Vascular")
gaba_cluster_labels <- c("0" = "Avp_GABA", "1" = "Vip_Nms_C1ql3_GABA",
                         "3" = "Grp_GABA")
# Specify the original Figure 5F/I genes by name rather than using indices such as
# p[[6]], which depend on database ordering.
figure_genes <- list(
  markers = c("Avp", "Vip", "Nms", "C1ql3", "Grp"),
  gaba_synthesis_transport = c("Gad1", "Gad2", "Slc32a1"),
  gaba_a_receptors = c("Gabra1", "Gabrb3", "Gabrb2", "Gabrg2"),
  gaba_b_receptors = c("Gabbr1", "Gabbr2"),
  vip = c("Vip", "Vipr2")
)

required_packages <- c("tidyverse", "Seurat", "CellChat", "Matrix", "sp",
                        "future", "patchwork", "ggvenn", "ggpubr",
                        "ComplexHeatmap", "circlize")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Please install the following packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(CellChat)
})

# 2. Shared functions for data loading and GABA subtyping ----
# Load .rda files into a separate environment to avoid overwriting workspace objects.
load_rda_object <- function(path, object_name) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (!object_name %in% loaded) stop(path, " is missing object: ", object_name)
  env[[object_name]]
}

save_figure <- function(plot, filename, width, height) {
  ggplot2::ggsave(file.path(output_dir, filename), plot = plot,
                 width = width, height = height, units = "in")
}

# Select SCN cells using each sample's v5 polygons, including boundary cells (> 0).
get_scn_cells <- function(coordinates, polygon_dir) {
  files <- list.files(polygon_dir, pattern = "_regions_v5\\.csv$", full.names = TRUE)
  names(files) <- gsub("_", "-", sub("_regions.*$", "", basename(files)))
  if (anyDuplicated(names(files))) stop("Multiple v5 polygon files found for the same sample.")
  if (!all(sample_ids %in% names(files)) || !all(sample_ids %in% names(coordinates))) {
    stop("SCN polygons or coordinates are missing one or more of the six required CT8/CT20 samples.")
  }
  unlist(lapply(sample_ids, function(sample_id) {
    polygon <- readr::read_csv(files[[sample_id]], skip = 2, show_col_types = FALSE) %>%
      dplyr::filter(grepl("^SCN_\\d+$", Selection))
    if (!nrow(polygon)) stop(sample_id, " has no SCN polygons.")
    xy <- coordinates[[sample_id]]
    cells <- unique(unlist(lapply(unique(polygon$Selection), function(region) {
      vertices <- polygon[polygon$Selection == region, c("X", "Y")]
      inside <- sp::point.in.polygon(xy$x, xy$y, vertices$X, vertices$Y) > 0
      xy$cell[inside]
    }), use.names = FALSE))
    message(sample_id, ": ", length(cells), " SCN cells")
    cells
  }), use.names = FALSE) %>% unique()
}

# Apply the same PCA, UMAP, and clustering parameters to the initial integrated
# dataset and to reclustering after removing non-target clusters.
cluster_gaba <- function(object, integration_features) {
  DefaultAssay(object) <- "XeniumIntegrated"
  VariableFeatures(object) <- integration_features
  object <- ScaleData(object, assay = "XeniumIntegrated",
                      features = rownames(object[["XeniumIntegrated"]]), verbose = FALSE)
  object <- RunPCA(object, assay = "XeniumIntegrated", npcs = 50,
                   reduction.name = "xeniumintegrated_pca",
                   reduction.key = "xeniumintegratedpc_", verbose = FALSE)
  object <- RunUMAP(object, reduction = "xeniumintegrated_pca", dims = 1:50,
                    reduction.name = "xeniumintegrated_umap",
                    reduction.key = "xeniumintegratedumap_", seed.use = random_seed)
  object <- FindNeighbors(object, reduction = "xeniumintegrated_pca", dims = 1:50,
                          graph.name = c("xeniumintegrated_nn", "xeniumintegrated_snn"))
  # v3 incorrectly read res.0.5 after clustering at resolution 0.3; use current Idents.
  object <- FindClusters(object, resolution = 0.3, graph.name = "xeniumintegrated_snn",
                         random.seed = random_seed)
  object$seurat_clusters <- Idents(object)
  object
}

# Run the v3 GABA integration and two-stage clustering only when caches are missing
# or rebuilding is requested. Preserve the original manual marker-based mapping
# from clusters to subtypes. Recheck this mapping after changing data or clustering
# parameters, because cluster numbers do not identify fixed biological cell types.
build_gaba_subtypes <- function(scn_object) {
  cells <- rownames(scn_object[[]])[scn_object$predicted.id == "GABAergic"]
  gaba <- subset(scn_object, cells = cells)
  split_objects <- SplitObject(gaba, split.by = "orig.ident")
  split_objects <- lapply(split_objects, function(object) {
    SCTransform(object, vars.to.regress = c("nCount_Xenium", "nFeature_Xenium"),
                assay = "Xenium", return.only.var.genes = FALSE, verbose = FALSE)
  })
  common_genes <- Reduce(intersect, lapply(split_objects, function(object) {
    rownames(GetAssayData(object, assay = "SCT", layer = "data"))
  }))
  # Preserve the original table() ordering of variable genes shared by all samples.
  feature_counts <- table(unlist(lapply(split_objects, function(object) {
    VariableFeatures(object, assay = "SCT")
  })))
  integration_features <- names(feature_counts)[feature_counts == length(split_objects)]
  # Preserve the original SCT assay input and the default integration methods used
  # by FindIntegrationAnchors and IntegrateData.
  anchors <- FindIntegrationAnchors(
    object.list = split_objects, dims = 1:50,
    assay = rep("SCT", length(split_objects)), anchor.features = integration_features,
    k.anchor = 5, k.filter = 199, k.score = 30
  )
  integrated <- IntegrateData(anchorset = anchors, dims = 1:50,
                              features.to.integrate = common_genes, k.weight = 80,
                              new.assay.name = "XeniumIntegrated")
  integrated <- cluster_gaba(integrated, integration_features)
  cells <- colnames(integrated)[!as.character(integrated$seurat_clusters) %in% c("2", "3")]
  integrated <- cluster_gaba(subset(integrated, cells = cells), integration_features)
  labels <- unname(gaba_cluster_labels[as.character(integrated$seurat_clusters)])
  annotations <- data.frame(cell = colnames(integrated), annot = labels) %>%
    tidyr::drop_na(annot)
  list(annotations = annotations, integrated = integrated)
}

# Merge GABA subtype labels, remove unclassified GABA cells, and retain cell types
# with more than 100 cells in total.
prepare_scn_metadata <- function(object, annotations) {
  if (!all(c("cell", "annot") %in% names(annotations)) || anyDuplicated(annotations$cell)) {
    stop("GABA_subclass must contain cell and annot columns with unique cell IDs.")
  }
  metadata <- object[[]] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(-any_of("annot")) %>%
    left_join(annotations[, c("cell", "annot")], by = "cell") %>%
    dplyr::filter(!(predicted.id == "GABAergic" & is.na(annot))) %>%
    mutate(annot = coalesce(annot, as.character(predicted.id)),
           predicted.id = annot,
           celltype = case_when(cluster == "5" ~ "Oligodendrocyte",
                                cluster == "3" ~ "Astrocyte",
                                cluster == "8" ~ "Ependymal",
                                TRUE ~ predicted.id)) %>%
    # Correct the original Nmns/Clql3 label spelling consistently across time points.
    mutate(celltype = if_else(grepl("Vip", celltype), "Vip_Nms_C1ql3_GABA", celltype))
  object <- subset(object, cells = metadata$cell)
  metadata <- tibble::column_to_rownames(metadata, "cell")
  object@meta.data <- metadata[colnames(object), , drop = FALSE]
  DefaultAssay(object) <- "Xenium"
  # Preserve the original order: remove genes unexpressed throughout the SCN before
  # applying the cell-type abundance threshold.
  counts <- GetAssayData(object, assay = "Xenium", layer = "counts")
  object <- subset(object, features = rownames(counts)[Matrix::rowSums(counts) > 0])
  celltype_counts <- table(object$celltype)
  cells <- colnames(object)[object$celltype %in% names(celltype_counts)[celltype_counts > 100]]
  object <- subset(object, cells = cells)
  if (any(!object$celltype %in% celltype_levels)) {
    stop("Cell types outside the six original figure categories were found; check celltype_levels.")
  }
  Idents(object) <- object$celltype
  object
}

# 3. CellChat inference and shared figure data ----
infer_cellchat <- function(object, coordinates, expressed_genes) {
  DefaultAssay(object) <- "Xenium"
  object <- NormalizeData(object, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = FALSE)
  metadata <- data.frame(
    labels = factor(object$celltype, levels = celltype_levels),
    samples = factor(object$orig.ident, levels = sample_ids),
    row.names = colnames(object)
  )
  metadata$samples <- droplevels(metadata$samples)
  if (anyNA(metadata)) stop("CellChat metadata contains missing cell-type or sample labels.")
  # ratio = 1 preserves the original coordinate units; tol = 5 is half the original
  # spot.size = 10.
  spatial_factors <- data.frame(ratio = rep(1, nlevels(metadata$samples)),
                                tol = rep(5, nlevels(metadata$samples)),
                                row.names = levels(metadata$samples))
  xy <- coordinates[colnames(object), c("x", "y"), drop = FALSE]
  if (anyNA(xy)) stop("Some CellChat cells are missing spatial coordinates.")
  cellchat <- CellChat::createCellChat(
    object = object, group.by = "labels", assay = "Xenium", datatype = "spatial",
    coordinates = xy, spatial.factors = spatial_factors, meta = metadata
  )
  cellchat@DB <- CellChat::CellChatDB.mouse
  cellchat <- CellChat::subsetData(cellchat, features = expressed_genes)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat, do.DE = FALSE, min.cells = 10)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat, variable.both = FALSE)
  cellchat <- CellChat::computeCommunProb(cellchat, distance.use = TRUE,
                                         contact.range = 10, interaction.range = 250,
                                         scale.distance = 4)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 5)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  CellChat::aggregateNet(cellchat)
}

# Use pathway summaries from subsetCommunication(netP), retaining its default LR
# filtering threshold of 0.05. The reported pval is the mean of included LR p-values,
# not a new pathway-level test. Extract the full table before filtering by pathway,
# so undetected pathways return an empty table.
pathway_table <- function(cellchat) {
  CellChat::subsetCommunication(cellchat, slot.name = "netP") %>%
    dplyr::select(source, target, pathway_name, prob, pval)
}

heatmap_data <- function(communication, pathway) {
  tidyr::expand_grid(source = celltype_levels, target = celltype_levels) %>%
    left_join(dplyr::filter(communication, pathway_name == pathway), by = c("source", "target")) %>%
    mutate(prob = replace_na(prob, 0), pval = replace_na(pval, 1),
           source = factor(source, levels = rev(celltype_levels)),
           target = factor(target, levels = celltype_levels))
}

# Share a probability color scale across comparison panels. Preserve the v4
# probability values for tile colors; asterisks indicate p < 0.01.
plot_pathway_heatmaps <- function(communications, pathways, hide_top_x = FALSE) {
  panels <- tidyr::expand_grid(pathway = pathways, timepoint = names(communications))
  tables <- purrr::map2(panels$pathway, panels$timepoint, function(pathway, timepoint) {
    heatmap_data(communications[[timepoint]], pathway)
  })
  maximum <- max(vapply(tables, function(df) max(df$prob), numeric(1)))
  if (maximum == 0) maximum <- 1
  plots <- lapply(seq_len(nrow(panels)), function(i) {
    df <- tables[[i]] %>% mutate(significance = if_else(pval < pval_cutoff, "*", ""))
    title <- sprintf("%s (%s)", panels$pathway[i],
                     ifelse(panels$timepoint[i] == "CT8", "CT-8", "CT-20"))
    plot <- ggplot(df, aes(target, source, fill = prob)) +
      geom_tile() + geom_text(aes(label = significance), size = 3) +
      scale_fill_gradient(low = "white", high = "red", limits = c(0, maximum)) +
      labs(title = title, fill = "communication\nprobability") + theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    if (hide_top_x && i <= nrow(panels) - length(communications)) {
      plot <- plot + theme(axis.text.x = element_blank()) + xlab(NULL)
    }
    plot
  })
  patchwork::wrap_plots(plots, ncol = length(communications), guides = "collect")
}

# Use log-normalized expression stored in CellChat for violin plots. Expand only
# selected genes to avoid converting the entire sparse expression matrix to dense.
# Bold gene labels indicate enriched LR genes at either time point.
plot_expression_violins <- function(objects, genes, signaling, title, show_x = TRUE) {
  enriched_genes <- unique(unlist(lapply(objects, function(object) {
    if (!signaling %in% object@LR$LRsig$pathway_name) return(character())
    result <- CellChat::extractEnrichedLR(object, signaling = signaling,
                                          geneLR.return = TRUE, enriched.only = TRUE,
                                          thresh = pval_cutoff)
    if (nrow(result$pairLR) == 0) character() else result$geneLR
  }), use.names = FALSE))
  data <- purrr::imap_dfr(objects, function(object, timepoint) {
    missing <- setdiff(genes, rownames(object@data))
    if (length(missing)) stop(timepoint, " is missing genes required for plotting: ", paste(missing, collapse = ", "))
    as.data.frame(as.matrix(object@data[genes, , drop = FALSE])) %>%
      tibble::rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "cell", values_to = "expr") %>%
      mutate(celltype = object@meta[cell, "labels"], timepoint = timepoint)
  }) %>% mutate(celltype = factor(celltype, levels = celltype_levels),
                 timepoint = factor(timepoint, levels = c("CT8", "CT20")))
  plots <- lapply(seq_along(genes), function(i) {
    gene <- genes[[i]]
    plot <- ggplot(dplyr::filter(data, .data$gene == .env$gene),
                   aes(celltype, expr, fill = timepoint)) +
      geom_violin() + theme_classic() + labs(y = gene, x = NULL) +
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1,
                                        face = ifelse(gene %in% enriched_genes, "bold", "plain")),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    if (!show_x || i < length(genes)) plot <- plot + theme(axis.text.x = element_blank())
    plot
  })
  patchwork::wrap_plots(plots, ncol = 1, guides = "collect") +
    patchwork::plot_annotation(title = title)
}

# Retain source-target combinations with p <= 0.01 at either time point; fill in 0
# for the other time point if it fails the threshold. Align and sort by source/target
# to ensure that observations are correctly matched for the paired t-test.
paired_probability_data <- function(communications, pathway) {
  significant <- lapply(communications, function(df) {
    df %>% dplyr::filter(pathway_name == pathway, pval <= pval_cutoff) %>%
      dplyr::select(source, target, prob)
  })
  paired <- full_join(rename(significant$CT8, CT8 = prob),
                      rename(significant$CT20, CT20 = prob), by = c("source", "target"))
  if (nrow(paired) < 2) stop(pathway, " has too few paired combinations for a paired t-test plot.")
  paired %>% mutate(CT8 = replace_na(CT8, 0), CT20 = replace_na(CT20, 0)) %>%
    arrange(source, target) %>%
    pivot_longer(c(CT8, CT20), names_to = "timepoint", values_to = "prob") %>%
    mutate(pathway_name = pathway, timepoint = factor(timepoint, levels = c("CT8", "CT20")))
}

plot_probability_comparison <- function(communications, pathway) {
  ggplot(paired_probability_data(communications, pathway),
         aes(pathway_name, prob, fill = timepoint)) +
    geom_boxplot() + theme_classic() +
    ggpubr::stat_compare_means(aes(group = timepoint), method = "t.test",
                              method.args = list(alternative = test_alternative),
                              paired = TRUE, label = "p.format")
}

# 4. Main workflow: load SCN cells, obtain subtypes, and create CellChat objects ----
main <- function() {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(random_seed)
  options(future.globals.maxSize = 8 * 1024^3) # 8 GiB
  object <- load_rda_object(input_files$seurat, "xenium_obj_filtered")
  coordinates <- load_rda_object(input_files$coordinates, "xy_coordinates")
  scn_cells <- intersect(get_scn_cells(coordinates, input_files$polygons), colnames(object))
  if (!length(scn_cells)) stop("No cells are shared between the SCN polygons and the Seurat object.")
  object <- subset(object, cells = scn_cells)

  if (rebuild_gaba_subtypes || !all(file.exists(unlist(cache_files[c("gaba", "integrated")])))) {
    message("Rebuilding GABA subtype and UMAP caches.")
    gaba <- build_gaba_subtypes(object)
    GABA_subclass <- gaba$annotations
    xenium_obj_filter_GABA_integrated <- gaba$integrated
    save(GABA_subclass, file = cache_files$gaba)
    save(xenium_obj_filter_GABA_integrated, file = cache_files$integrated)
    rm(gaba)
  } else {
    GABA_subclass <- load_rda_object(cache_files$gaba, "GABA_subclass")
    xenium_obj_filter_GABA_integrated <- load_rda_object(
      cache_files$integrated, "xenium_obj_filter_GABA_integrated"
    )
  }
  object <- prepare_scn_metadata(object, GABA_subclass)

  # Figure 5B: show the three original GABA subtypes using integrated UMAP coordinates.
  integrated <- xenium_obj_filter_GABA_integrated
  cells <- colnames(integrated)[as.character(integrated$seurat_clusters) %in% names(gaba_cluster_labels)]
  p_umap <- DimPlot(integrated, group.by = "seurat_clusters", label = TRUE,
                    reduction = "xeniumintegrated_umap", repel = TRUE, cells = cells) +
    scale_color_manual(name = "Cell type",
                        values = c("0" = "#F8766D", "1" = "#00BA38", "3" = "#619CFF"),
                        labels = gaba_cluster_labels)
  save_figure(p_umap, "GABA_subtype_umap.pdf", 6.58, 4.16)
  rm(integrated, xenium_obj_filter_GABA_integrated)

  # Reuse the feature set after prepare_scn_metadata removes zero-expression genes
  # in the original processing order.
  counts <- GetAssayData(object, assay = "Xenium", layer = "counts")
  expressed_genes <- rownames(object)

  # GABA marker dotplot (the original subtype validation plot): convert counts to CPM,
  # then calculate subtype means and the fraction of cells expressing each gene.
  # Rescale each marker's mean expression across the three subtypes to 0-5 for color;
  # this is not a z-score.
  if (any(Matrix::colSums(counts) == 0)) stop("Some cells have zero total counts.")
  cpm <- counts %*% Matrix::Diagonal(x = 1e6 / Matrix::colSums(counts))
  dimnames(cpm) <- dimnames(counts)
  if (!all(figure_genes$markers %in% rownames(cpm))) stop("Required GABA marker genes are missing.")
  dot_data <- purrr::map_dfr(celltype_levels[grepl("GABA", celltype_levels)], function(subtype) {
    cells <- colnames(object)[object$celltype == subtype]
    if (!length(cells)) stop("Missing GABA subtype: ", subtype)
    mat <- cpm[figure_genes$markers, cells, drop = FALSE]
    data.frame(gene = figure_genes$markers, norm_expr = Matrix::rowMeans(mat),
               percent = Matrix::rowMeans(mat != 0), celltype = subtype)
  }) %>% group_by(gene) %>%
    mutate(scaled_expr = scales::rescale(norm_expr, to = c(0, 5))) %>% ungroup() %>%
    mutate(gene = factor(gene, levels = figure_genes$markers))
  p_dot <- ggplot(dot_data, aes(gene, celltype, size = percent, color = scaled_expr)) +
    geom_point() + scale_color_gradient(low = "blue", high = "red") + theme_classic()
  save_figure(p_dot, "GABA_neuron_subtype_dotplot.pdf", 4.64, 4.16)
  rm(counts, cpm)

  coordinates <- bind_rows(coordinates[sample_ids]) %>%
    dplyr::select(cell, x, y) %>% tibble::column_to_rownames("cell")
  if (recompute_cellchat || !file.exists(cache_files$cellchat)) {
    previous_plan <- future::plan(future::multisession, workers = future_workers)
    on.exit(future::plan(previous_plan), add = TRUE)
    objects <- lapply(list(CT8 = sample_ids[1:3], CT20 = sample_ids[4:6]), function(samples) {
      cells <- colnames(object)[object$orig.ident %in% samples]
      if (!length(cells)) stop("No cells are available for the time point with samples: ", paste(samples, collapse = ", "))
      infer_cellchat(subset(object, cells = cells), coordinates, expressed_genes)
    })
    cellchat_CT8 <- objects$CT8
    cellchat_CT20 <- objects$CT20
    save(cellchat_CT8, cellchat_CT20, file = cache_files$cellchat)
    future::plan(previous_plan)
  } else {
    objects <- list(CT8 = load_rda_object(cache_files$cellchat, "cellchat_CT8"),
                     CT20 = load_rda_object(cache_files$cellchat, "cellchat_CT20"))
  }
  # Subsequent figures do not use centrality, so omit netAnalysis_computeCentrality.
  rm(object, coordinates)
  communications <- lapply(objects, pathway_table)

  # 5. Overall communication: count heatmaps, Figure 5C/D, and Figure S12A ----
  # Use native ComplexHeatmap instead of the missing ~/Downloads/netVisual_heatmap_custom.R.
  # Use the same cell-type order and 0-60 color scale for both time points, one per page.
  # Values above 60 use the darkest red.
  write_count_heatmaps <- function() {
    grDevices::pdf(file.path(output_dir, "communication_count_heatmap.pdf"), 5.75, 4.94)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (timepoint in names(objects)) {
      mat <- objects[[timepoint]]@net$count
      order <- intersect(celltype_levels, rownames(mat))
      mat <- mat[order, order, drop = FALSE]
      heatmap <- ComplexHeatmap::Heatmap(
        mat, name = "Interactions", column_title = timepoint,
        col = circlize::colorRamp2(c(0, 60), c("white", "red")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_title = "Source", column_names_rot = 45
      )
      ComplexHeatmap::draw(heatmap)
    }
  }
  write_count_heatmaps()

  # Figure 5C: combine the three GABA source/target subtypes before summing LR counts.
  count_data <- purrr::imap_dfr(objects, function(object, timepoint) {
    as.data.frame(object@net$count) %>% tibble::rownames_to_column("source") %>%
      pivot_longer(-source, names_to = "target", values_to = "interaction") %>%
      mutate(source = if_else(grepl("GABA", source), "GABAergic", source),
             target = if_else(grepl("GABA", target), "GABAergic", target)) %>%
      group_by(source, target) %>% summarise(interaction = sum(interaction), .groups = "drop") %>%
      mutate(communication = sprintf("%s -> %s", source, target), timepoint = timepoint)
  })
  communication_order <- count_data %>% group_by(communication) %>%
    summarise(total = sum(interaction), .groups = "drop") %>% arrange(desc(total)) %>%
    pull(communication)
  count_data <- count_data %>%
    mutate(communication = factor(communication, levels = rev(communication_order)),
           timepoint = factor(timepoint, levels = c("CT8", "CT20")))
  p_count <- ggplot(count_data, aes(interaction, communication, fill = timepoint)) +
    geom_col() + theme_classic() + xlab("number of communications")
  save_figure(p_count, "communication_count_barplot.pdf", 5.75, 4.94)

  # Figure 5D: sets of significant LR pairs; Figure S12A: sets of significant pathways.
  p_lr_venn <- ggvenn::ggvenn(
    list(`CT8 (LRs)` = objects$CT8@net$LRs, `CT20 (LRs)` = objects$CT20@net$LRs),
    fill_color = c("#e37365", "#4cb3bb"), text_size = 6
  )
  save_figure(p_lr_venn, "communication_count_venn.pdf", 5.75, 4.94)
  p_pathway_venn <- ggvenn::ggvenn(
    list(`CT8 (LR signaling pathways)` = objects$CT8@netP$pathways,
         `CT20 (LR signaling pathways)` = objects$CT20@netP$pathways),
    fill_color = c("#e37365", "#4cb3bb"), text_size = 8
  )
  save_figure(p_pathway_venn, "communication_pathway_venn.pdf", 5.75, 4.94)

  # 6. GABA: Figure 5E/F/G ----
  p_gaba_heatmap <- plot_pathway_heatmaps(communications, c("GABA-A", "GABA-B"), hide_top_x = TRUE)
  save_figure(p_gaba_heatmap, "GABA_A_B.heatmap_V2.pdf", 8.96, 5.45)
  gaba_expression <- list(
    plot_expression_violins(objects, figure_genes$gaba_synthesis_transport, "GABA-A",
                            "GABA LR pathway (enzymes and transporters)", show_x = FALSE),
    plot_expression_violins(objects, figure_genes$gaba_a_receptors, "GABA-A",
                            "GABA-A LR pathway (receptors)", show_x = FALSE),
    plot_expression_violins(objects, figure_genes$gaba_b_receptors, "GABA-B",
                            "GABA-B LR pathway (receptors)")
  )
  p_gaba_expression <- patchwork::wrap_plots(gaba_expression, ncol = 1, heights = c(3, 4, 2))
  save_figure(p_gaba_expression, "GABA_LR_pathway_gene_expression.pdf", 6.90, 7.68)
  p_gaba_comparison <- patchwork::wrap_plots(
    plot_probability_comparison(communications, "GABA-A"),
    plot_probability_comparison(communications, "GABA-B"), ncol = 2, guides = "collect"
  )
  save_figure(p_gaba_comparison, "GABA_communication_probability.pdf", 6.90, 4.16)

  # 7. VIP: Figure 5H/I/J ----
  # Use distinct names to avoid overwriting p_VIP with a contribution plot as in the original.
  p_vip_heatmap <- plot_pathway_heatmaps(communications, "VIP")
  save_figure(p_vip_heatmap, "VIP_communication_heatmap.pdf", 8.96, 4.18)
  p_vip_expression <- plot_expression_violins(objects, figure_genes$vip, "VIP",
                                              "VIP LR pathway (gene expression)")
  save_figure(p_vip_expression, "VIP_LR_pathway_gene_expression.pdf", 7.34, 4.14)
  p_vip_comparison <- plot_probability_comparison(communications, "VIP")
  save_figure(p_vip_comparison, "VIP_communication_probability.pdf", 5.75, 4.16)

  # 8. Additional original figure exports and Figure S12B ----
  # Retain the original shared-pathway PDF. Explicitly print each plot inside the
  # loop so that Rscript writes every page.
  write_shared_heatmaps <- function() {
    shared <- intersect(objects$CT8@netP$pathways, objects$CT20@netP$pathways)
    if (!length(shared)) return(invisible(NULL))
    grDevices::pdf(file.path(output_dir, "Cellchat_pathmay_both.pdf"), 9.94, 4.18)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (pathway in shared) print(plot_pathway_heatmaps(communications, pathway))
  }
  write_shared_heatmaps()

  # Figure S12B: retain the stacked rankNet comparison selected in the original script.
  merged <- CellChat::mergeCellChat(objects, add.names = names(objects))
  p_ranknet <- CellChat::rankNet(merged, mode = "comparison", stacked = TRUE)
  save_figure(p_ranknet, "cellchat_ranknet.pdf", 8.33, 6.59)
  message("Figure 5 and Figure S12 have been saved to: ", output_dir)
  invisible(NULL)
}

main()
