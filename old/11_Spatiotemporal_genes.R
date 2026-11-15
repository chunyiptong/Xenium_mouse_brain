library(Seurat)
library(tidyverse)
setwd("~/Dropbox/Spatial_transcriptomics_new/spatiotemporal")

#1. load data ----
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

xy_coordinates %>% 
  purrr::map(function(df_){
    df_
    df_ %>% dplyr::filter(cell %in% colnames(xenium_obj_filtered)) -> df_
  }) -> xy_coordinates

meta_data = xenium_obj_filtered@meta.data

# import contours
list.files(path = "~/Dropbox/Spatial_transcriptomics_new/rda/Brain_region_annotation/SCN_SON_PVN", pattern = "_regions_v5\\.csv", full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)_regions.+", "\\1", .) %>% gsub("_", "-", .)) %>% 
  purrr::map(function(file_){
    read_csv(file_, skip = 2) -> df_
  })-> polygon_

select_cells = function(region__){
  c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
    "names<-"(.,.) %>% 
    #  .[3] %>% 
    purrr::map(function(group_){
      xy_coordinates[[group_]] -> xy_coordinates
      unique(polygon_[[group_]]$Selection) %>% 
        .[grepl(sprintf("^%s_\\d+", region__), .)] %>% 
        "names<-"(.,.) %>% 
        purrr::map(function(region_){
          polygon_[[group_]] %>% dplyr::filter(Selection == region_) -> df_ 
          print(head(df_)$Selection %>% unique())
          df_ %>% dplyr::select(X,Y) -> df_
          cells_se = sp::point.in.polygon(
            point.x = xy_coordinates$x,
            point.y = xy_coordinates$y,
            pol.x = df_$X,
            pol.y = df_$Y) > 0
          xy_coordinates[cells_se, ] -> df_1
          df_1$cell
        }) %>% unlist() %>% unname() -> cells_
      print(sprintf("%s, n_cells = %s", group_, length(cells_)))
      return(cells_)
    }) %>% unlist() %>% unname() -> cells
}

scn_cells = select_cells(region__ = "SCN")
scn_core_cells = select_cells(region__ = "SCN_core")
scn_shell_cells = select_cells(region__ = "SCN_shell")

xenium_obj_filtered = subset(xenium_obj_filtered, cells = scn_cells)
GetAssayData(xenium_obj_filtered, assay = "Xenium", layer = "counts") -> counts_matrix
counts_matrix %>% colSums() %>% median() -> median_
counts_matrix %>% as.data.frame() %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*median_})) -> norm_counts_matrix

list(
  SCN = scn_cells
) %>% 
  purrr::map(function(cells_){
    c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
      #      .[1] %>% 
      purrr::map(function(CT_){
        cells_ %>% {.[grepl(CT_, .)]} -> cells_
        counts_matrix[,cells_] %>% rowSums() %>% as.data.frame() %>% "colnames<-"(., CT_)
      }) %>% do.call(cbind, .) -> df_count_mat
  }) -> df_count_mat_list

df_count_mat_list %>% 
  purrr::map(function(df_){
    df_ %>% colSums() %>% median() -> median_
    df_ %>% 
      dplyr::mutate(across(everything(), function(x){(x/sum(x))*median_})) -> df_count_mat_norm
  }) -> df_count_mat_norm_list

save(counts_matrix, norm_counts_matrix, df_count_mat_list, df_count_mat_norm_list, file = "~/Downloads/SCN_count_matrix.rda")

list(
  SCN_core = scn_core_cells,
  SCN_shell = scn_shell_cells
) %>% 
  #  .[1] %>% 
  purrr::map(function(cells_){
    c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
      "names<-"(.,.) %>% 
      #      .[1] %>% 
      purrr::map(function(sample_){
        cells_ %>% {.[grepl(sample_, cells_)]} -> cells_
        counts_matrix[, cells_] %>% rowSums() %>% as.data.frame() %>% "colnames<-"(., sample_)
      }) %>% do.call(cbind, .) -> df_
  }) -> SCN_core_shell_matrix

SCN_core_shell_matrix %>% 
  purrr::map(function(df_){
    df_ %>% colSums()
  }) %>% unlist() %>% median() -> median_

SCN_core_shell_matrix %>% 
  purrr::map(function(df_){
    df_ %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*median_}))
  }) -> SCN_core_shell_matrix_norm

save(SCN_core_shell_matrix, SCN_core_shell_matrix_norm, file = "~/Downloads/SCN_core_shell_matrix_norm.Rdata")

#2. Generate data for spatiotemporal function ----
xenium_obj_filtered@meta.data %>% dplyr::select(orig.ident) %>% dplyr::mutate(CT = gsub("(CT-\\d+?)-.+", "\\1", orig.ident)) -> Time_df
xy_coordinates %>% do.call(rbind, .) %>% "rownames<-"(., NULL) %>% column_to_rownames("cell") %>% .[rownames(Time_df), ] -> coordinate_df
raw_counts = GetAssayData(xenium_obj_filtered, assay = "Xenium", layer = "counts")
dim(raw_counts)
raw_counts %>% colSums() %>% median() -> median_
raw_counts %>% as.data.frame() %>% 
  dplyr::mutate(across(everything(), function(x){(x/sum(x))*median_})) -> norm_conuts
dim(norm_conuts)
norm_conuts %>% t() %>% as.data.frame() %>% .[rownames(Time_df), ] -> norm_conuts
dim(norm_conuts)
norm_conuts=norm_conuts[,colSums(norm_conuts)>0]
dim(norm_conuts)

xenium_obj_filtered@meta.data[rownames(Time_df), ] %>% 
  dplyr::mutate(celltype = case_when(
    cluster == "5" ~ "Oligodendrocytes",
    cluster == "3" ~ "Astrocytes",
    cluster == "8" ~ "Ependymal",
    TRUE ~ predicted.id
  )) %>% dplyr::select(celltype) -> celltype_df

dim(Time_df)
dim(norm_conuts)
dim(coordinate_df)
dim(celltype_df)

boxplot(coordinate_df$x~Time_df$orig.ident)
boxplot(coordinate_df$y~Time_df$orig.ident)

standardize_coords_by_slice <- function(x, y, slice_id, center = FALSE) {
  # x, y: numeric vectors
  # slice_id: factor or character vector (length n)
  # center: if TRUE, also center to [-0.5, 0.5]
  
  stopifnot(length(x) == length(y))
  stopifnot(length(x) == length(slice_id))
  
  slice_id <- as.factor(slice_id)
  n <- length(x)
  
  x_std <- numeric(n)
  y_std <- numeric(n)
  
  for (s in levels(slice_id)) {
    idx <- which(slice_id == s)
    
    x_s <- x[idx]
    y_s <- y[idx]
    
    # handle constant edge case
    range_x <- max(x_s) - min(x_s)
    range_y <- max(y_s) - min(y_s)
    
    if (range_x == 0) {
      x_std[idx] <- 0.5
    } else {
      x_std[idx] <- (x_s - min(x_s)) / range_x
    }
    
    if (range_y == 0) {
      y_std[idx] <- 0.5
    } else {
      y_std[idx] <- (y_s - min(y_s)) / range_y
    }
    
    if (center) {
      x_std[idx] <- x_std[idx] - 0.5
      y_std[idx] <- y_std[idx] - 0.5
    }
  }
  
  data.frame(
    x_std = x_std,
    y_std = y_std
  )
}

coordinate_aligned=standardize_coords_by_slice(coordinate_df$x, 
                                               coordinate_df$y,
                                               Time_df$orig.ident,
                                               center=TRUE)

boxplot(coordinate_aligned$x_std~Time_df$orig.ident)
boxplot(coordinate_aligned$y_std~Time_df$orig.ident)

##2.1 plot aligned SCN ----
df_plot <- data.frame(
  x = coordinate_aligned$x_std,
  y = coordinate_aligned$y_std,
  slice = as.factor(Time_df$orig.ident)
)
ggplot(df_plot, aes(x = x, y = y, color = slice)) +
  geom_point(size = 0.4, alpha = 0.7) +
  coord_equal() +
  theme_classic() +
  labs(
    title = "Standardized Spatial Coordinates by Slice",
    x = "x (standardized)",
    y = "y (standardized)",
    color = "Slice"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

colnames(norm_conuts)
time <- ifelse(Time_df$CT == "CT-20", 1, 0)
#y <- norm_conuts[, "Vip"]

#3. Run spatiotemporal function ----
gaussian_rbf_basis <- function(coords, centers, length_scale) {
  # coords: n x 2
  # centers: K x 2
  # returns: n x K basis matrix
  
  n <- nrow(coords)
  K <- nrow(centers)
  Phi <- matrix(0, n, K)
  
  for (k in seq_len(K)) {
    dx <- coords[, 1] - centers[k, 1]
    dy <- coords[, 2] - centers[k, 2]
    d2 <- dx^2 + dy^2
    Phi[, k] <- exp(-d2 / (2 * length_scale^2))
  }
  Phi
}

####################################
# Step 1: Fast gene-wise testing
####################################
run_all_genes_fast <- function(expr_mat,
                               Phi,
                               Phi_time,
                               time,
                               add_covariates = NULL,
                               min_var = 1e-8) {
  
  expr_mat <- as.matrix(expr_mat)
  G <- ncol(expr_mat)
  n <- nrow(expr_mat)
  
  gene_names <- colnames(expr_mat)
  if (is.null(gene_names)) {
    gene_names <- paste0("gene", seq_len(G))
  }
  
  # base design matrices
  X_reduced <- cbind(1, time, Phi)
  X_full <- cbind(1, time, Phi, Phi_time)
  
  colnames(X_reduced) <- c("Intercept", "time", colnames(Phi))
  colnames(X_full) <- c("Intercept", "time", colnames(Phi), colnames(Phi_time))
  
  if (!is.null(add_covariates)) {
    add_covariates <- as.matrix(add_covariates)
    X_reduced <- cbind(X_reduced, add_covariates)
    X_full <- cbind(X_full, add_covariates)
  }
  
  pvals <- rep(NA_real_, G)
  dfs <- rep(NA_real_, G)
  
  for (g in seq_len(G)) {
    y <- expr_mat[, g]
    
    if (!all(is.na(y)) && var(y, na.rm = TRUE) > min_var) {
      
      # remove missing
      keep <- is.finite(y)
      y_sub <- y[keep]
      
      Xr <- X_reduced[keep, , drop = FALSE]
      Xf <- X_full[keep, , drop = FALSE]
      
      # fit models via lm.fit (fast)
      fit_r <- lm.fit(Xr, y_sub)
      fit_f <- lm.fit(Xf, y_sub)
      
      # residual sums of squares
      rss_r <- sum(fit_r$residuals^2)
      rss_f <- sum(fit_f$residuals^2)
      
      df_r <- length(y_sub) - ncol(Xr)
      df_f <- length(y_sub) - ncol(Xf)
      
      df_diff <- df_r - df_f
      
      if (df_diff > 0 && df_f > 0) {
        F_stat <- ((rss_r - rss_f) / df_diff) / (rss_f / df_f)
        pvals[g] <- pf(F_stat, df_diff, df_f, lower.tail = FALSE)
        dfs[g] <- df_diff
      }
    }
  }
  
  data.frame(
    gene = gene_names,
    p_value_interaction = pvals,
    df_interaction = dfs,
    fdr_interaction = p.adjust(pvals, method = "BH")
  )
}



cell_type <- as.factor(celltype_df$celltype)
celltype_mat <- model.matrix(~ cell_type - 1)
colnames(celltype_mat) <- sub("^cell_type", "", colnames(celltype_mat))
head(celltype_mat)

slice_id <- as.factor(Time_df$orig.ident)
slice_mat <- model.matrix(~ slice_id - 1)
colnames(slice_mat) <- sub("^slice_id", "slice_", colnames(slice_mat))
head(slice_mat)

add_covars <- cbind(slice_mat, celltype_mat)

results.all=results.cov.adj.all=matrix(nrow=ncol(norm_conuts), ncol=25)
rownames(results.all)=rownames(results.cov.adj.all)=colnames(norm_conuts)

i=1
for(n.centers in c(5,10,15,20,25)){
  for(bandwidth in c(0.5,1,1.5,2,2.5)){
    
    cat(i,'\t', 'n.centers', n.centers, 'bandwidth', bandwidth)
    
    # Precompute centers + basis (DO THIS ONCE)
    
    # coords_std: n x 2 matrix/data.frame
    coords_mat <- as.matrix(coordinate_aligned)
    
    # 1. choose centers once
    set.seed(1)
    centers <- kmeans(coords_mat,
                      centers = n.centers,
                      nstart = 25,
                      iter.max = 200)$centers
    
    # 2. choose length scale once
    dmat <- as.matrix(dist(rbind(coords_mat, centers)))
    n <- nrow(coords_mat)
    d_cc <- dmat[seq_len(n), n + seq_len(nrow(centers)), drop = FALSE]
    length_scale <- bandwidth*median(apply(d_cc, 1, min))
    
    # 3. compute basis once
    Phi <- gaussian_rbf_basis(coords_mat, centers, length_scale)
    colnames(Phi) <- paste0("rbf", seq_len(ncol(Phi)))
    
    # 4. interaction basis
    Phi_time <- Phi * time
    colnames(Phi_time) <- paste0(colnames(Phi), "_time")
    
    # Run it
    
    results <- run_all_genes_fast(
      expr_mat = norm_conuts,
      Phi = Phi,
      Phi_time = Phi_time,
      time = time
    )
    results.all[,i]=results$p_value_interaction
    
    # Add additional covariates including cell type and slice
    results.cov.adj <- run_all_genes_fast(
      expr_mat = norm_conuts,
      Phi = Phi,
      Phi_time = Phi_time,
      time = time,
      add_covariates = add_covars
    )
    results.cov.adj.all[,i]=results.cov.adj$p_value_interaction
    
    cat('\tresults sig. genes', sum(results$fdr_interaction<0.05),'\n')
    cat('\tresults.cov.adj sig. genes', sum(results.cov.adj$fdr_interaction<0.05),'\n\n')
    
    i=i+1
  }
}

save.image(file='~/Dropbox/Spatial_transcriptomics_new/spatiotemporal/spatiotemporal_cauchy_2.rda')

load('~/Dropbox/Spatial_transcriptomics_new/spatiotemporal/spatiotemporal_cauchy_2.rda')
cauchy.test=function(p.vector){
  d <- length(p.vector)
  
  if(1 %in% p.vector){
    p.vector[p.vector==1]=max(0.9999, 1-1/length(p.vector))
  }
  
  Sd <- sum(tan((0.5-p.vector)*pi)/d)
  p.global <- pcauchy(Sd,lower.tail = F)
  p.global <- min(p.global,1)
  
  return(p.global)
}

cauchy.test(results.all['Vip',])
cauchy.test(results.all['Grpr',])


results.pvals=apply(results.all, 1, cauchy.test)
results.padj=p.adjust(results.pvals, method = "BH")
results.cov.adj.pvals=apply(results.cov.adj.all,1,cauchy.test)
results.cov.adj.padj=p.adjust(results.cov.adj.pvals, method = "BH")


plot(results.pvals, results.cov.adj.pvals)

spatiotemp.results=data.frame(spatiotemp.pval=results.cov.adj.pvals,
                              spatiotemp.padj=results.cov.adj.padj)

save(spatiotemp.results, file='~/Dropbox/Spatial_transcriptomics_new/spatiotemporal/spatiotemp.results_2.rda')

#4. Downstream analysis ----
load('~/Dropbox/Spatial_transcriptomics_new/spatiotemporal/spatiotemp.results_2.rda')
#load("~/Downloads/brain_spatiotemporal_genes")
load("~/Downloads/df_combinedPval_all_samples_sparkx.Rdata")
load("~/Downloads/cauchy_scn_core_vs_shell.Rdata")
load("~/Downloads/DEG_dfs.rda")
#load("~/Downloads/new_spade_res.Rdata")

cauchy_scn_core_vs_shell %>%
  rowwise() %>% 
  dplyr::mutate(`CT-8_avg_log2FC` = mean(c(`CT-8-1_avg_log2FC`, `CT-8-2_avg_log2FC`, `CT-8-3_avg_log2FC`))) %>% 
  dplyr::mutate(`CT-20_avg_log2FC` = mean(c(`CT-20-1_avg_log2FC`, `CT-20-2_avg_log2FC`, `CT-20-3_avg_log2FC`))) %>% 
  as.data.frame() -> cauchy_scn_core_vs_shell

spatiotemp.results %>% rownames_to_column("gene") -> spatiotemp.results
list(
  sparkx = df_combinedPval_all_samples_sparkx,
  core_vs_shell = cauchy_scn_core_vs_shell,
  spatiotemporal = spatiotemp.results
) %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, name_){
    df_ %>% as_tibble() %>% 
      column_to_rownames("gene") -> df_
    colnames(df_) = sprintf("%s_%s", colnames(df_), name_)
    df_ %>% 
      rownames_to_column("gene") -> df_
  }) %>% purrr::reduce(., full_join, by = "gene") -> df_spatiotemporal_all
df_spatiotemporal_all[is.na(df_spatiotemporal_all)] <- 1

df_spatiotemporal_all %>% dplyr::filter(gene %in% spatiotemp.results$gene) -> df_spatiotemporal_all

DEG_df_SCN[,c(1,6,9,10)] %>% 
  column_to_rownames("gene") %>% 
  {
    df_ = .
    colnames(df_) = sprintf("%s_TVG", colnames(df_))
    df_ %>% 
      rownames_to_column("gene")
  } %>% 
  inner_join(x = df_spatiotemporal_all, y = ., by = "gene") -> df_spatiotemporal_all_1

df_spatiotemporal_all_1 %>% 
  dplyr::mutate(`CT-8_cauchy_padj_sparkx` = p.adjust(`CT-8_cauchy_pval_sparkx`, method = "BH")) %>% 
  dplyr::mutate(`CT-20_cauchy_padj_sparkx` = p.adjust(`CT-20_cauchy_pval_sparkx`, method = "BH")) %>% 
  dplyr::mutate(CT8_cauchy_padj_core_vs_shell = p.adjust(CT8_cauchy_pval_core_vs_shell, method = "BH")) %>% 
  dplyr::mutate(CT20_cauchy_padj_core_vs_shell = p.adjust(CT20_cauchy_pval_core_vs_shell, method = "BH")) %>% 
  dplyr::mutate(p_adj_TVG = p.adjust(p_val_TVG, method = "BH")) -> df_spatiotemporal_all_2

list(
  CT8_SVG = df_spatiotemporal_all_2 %>% dplyr::filter(`CT-8_cauchy_padj_sparkx` < 0.05) %>% .$gene,
  CT20_SVG = df_spatiotemporal_all_2 %>% dplyr::filter(`CT-20_cauchy_padj_sparkx` < 0.05) %>% .$gene
) %>% 
  ggvenn::ggvenn() -> SVG_CT20_CT8_venn

list(
  CT8_SVG = df_spatiotemporal_all_2 %>% dplyr::filter(`CT-8_cauchy_padj_sparkx` < 0.05) %>% .$gene,
  CT20_SVG = df_spatiotemporal_all_2 %>% dplyr::filter(`CT-20_cauchy_padj_sparkx` < 0.05) %>% .$gene,
  TVG = df_spatiotemporal_all_2 %>% dplyr::filter(p_adj_TVG < 0.05, abs(logFC_TVG) > log2(1.1)) %>% .$gene, 
  Spatiotemporal = df_spatiotemporal_all_2 %>% dplyr::filter(spatiotemp.padj_spatiotemporal < 0.05) %>% .$gene
) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::mutate(across(!matches("key"), function(x){x %>% as.integer()})) %>%
  as.data.frame() %>% 
  UpSetR::upset(sets = rev(c("CT20_SVG", "CT8_SVG", "TVG", "Spatiotemporal")), keep.order = T, order.by = "freq") -> p_upset ##4.1 upset plot (SVG vs TVG vs Spatiotemporal)----
pdf(file = "~/Downloads/sptiotemporal_upset_plot.pdf", width = 6.04, height = 3.26)
p_upset
dev.off()

list(
  SVG = df_spatiotemporal_all_2 %>% dplyr::filter(((`CT-8_cauchy_padj_sparkx` < 0.05)|(`CT-20_cauchy_padj_sparkx` < 0.05))) %>% .$gene,
  TVG = df_spatiotemporal_all_2 %>% dplyr::filter(p_adj_TVG < 0.05, abs(logFC_TVG) > log2(1.1)) %>% .$gene,
  Spatiotemporal = df_spatiotemporal_all_2 %>% dplyr::filter(spatiotemp.padj_spatiotemporal < 0.05) %>% .$gene
) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::filter(SVG == TRUE, TVG == TRUE, Spatiotemporal == TRUE) %>% .$key -> spatiotemporal_genes

df_spatiotemporal_all_2 %>% dplyr::filter(gene %in% spatiotemporal_genes) %>% 
  dplyr::mutate(TVG_group = case_when(
    (p_adj_TVG < 0.05)&(logFC_TVG > 0) ~ "CT-20", 
    (p_adj_TVG < 0.05)&(logFC_TVG < 0) ~ "CT-8", 
    TRUE ~ "None"
  )) %>% 
  dplyr::mutate(SVG_CT8_group = case_when(
    (CT8_cauchy_padj_core_vs_shell < 0.05)&(`CT-8_avg_log2FC_core_vs_shell` > 0) ~ "SCN core",
    (CT8_cauchy_padj_core_vs_shell < 0.05)&(`CT-8_avg_log2FC_core_vs_shell` < 0) ~ "SCN shell",
    TRUE ~ "None"
  )) %>% 
  dplyr::mutate(SVG_CT20_group = case_when(
    (CT20_cauchy_padj_core_vs_shell < 0.05)&(`CT-20_avg_log2FC_core_vs_shell` > 0) ~ "SCN core",
    (CT20_cauchy_padj_core_vs_shell < 0.05)&(`CT-20_avg_log2FC_core_vs_shell` < 0) ~ "SCN shell",
    TRUE ~ "None"
  )) %>% 
  dplyr::filter(((SVG_CT8_group != "None")|(SVG_CT20_group != "None"))) -> df_1

df_1 %>% dplyr::select(gene, TVG_group, SVG_CT8_group, SVG_CT20_group) %>% 
  dplyr::arrange(TVG_group, SVG_CT8_group, SVG_CT20_group, gene) %>% 
  write_csv("~/Downloads/new_spatiotemporal_genes_df.csv") ##4.2 58 spatiotemporal genes ----

load('~/Downloads/SCN_count_matrix.rda')
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")
load("~/Downloads/SCN_core_shell_matrix_norm.Rdata")

SCN_core_shell_matrix %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, name_){
    df_ %>% colSums()
  }) %>% unlist() %>% median() -> median_

SCN_core_shell_matrix %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, name_){
    colnames(df_) = sprintf("%s_%s", name_, colnames(df_))
    df_ %>% rownames_to_column("gene")
  }) %>% purrr::reduce(., inner_join, by = "gene") %>% 
  column_to_rownames("gene") %>% 
  dplyr::mutate(across(everything(), function(x){(x/sum(x))*median_})) -> SCN_core_shell_matrix_norm

list.files(path = "~/Dropbox/Spatial_transcriptomics_new/rda/Brain_region_annotation/SCN_SON_PVN", pattern = "_regions_v5\\.csv", full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)_regions.+", "\\1", .) %>% gsub("_", "-", .)) %>% 
  purrr::map(function(file_){
    read_csv(file_, skip = 2) -> df_
  })-> polygon_

c("Prkaa2", "Lhx1", "Rgs16", "Per1", "Per2", "Vip", "Pcsk1", "Slc39a6", "Vip", "Cd24a", "Rgs16", "Rpa1", "Grin2c") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    norm_counts_matrix[gene_, ] %>% as.numeric() -> gene_expr
    upper = quantile(gene_expr, 0.99, na.rm = TRUE)
    c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(sample_){
        norm_counts_matrix[gene_, ] %>% 
          t() %>% 
          as.data.frame() %>% "colnames<-"(., "normalized_expr") %>% 
          rownames_to_column("cell") -> df_
        xy_coordinates[[sample_]] -> df_xy
        inner_join(x = df_, y = df_xy, by = "cell") -> df_plot
        
        print(upper)
        df_plot %>% 
          dplyr::arrange(normalized_expr) -> df_plot
        
        df_plot %>% dplyr::mutate(normalized_expr_1 = case_when(
          normalized_expr > upper ~ upper,
          TRUE ~ normalized_expr
        )) -> df_plot
        #        upper = quantile(df_plot$normalized_expr, 0.99, na.rm = TRUE)
        df_plot %>% 
          ggplot(aes(x = x, y = y, color = normalized_expr, size = normalized_expr_1)) + 
#          geom_point(size = 1) + 
          geom_point() + 
          scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, upper), oob = scales::squish) + 
          scale_size_continuous(range = c(1,4)) +
          theme_void() + 
          labs(color = "Normalized\nExpression") +
          ggtitle(sprintf("%s", sample_)) -> p
        
        polygon_[[sample_]]$Selection %>% unique() %>% {.[grepl("shell|core", .)]} -> regions_
        for (region_ in regions_){
          polygon_[[sample_]] %>% dplyr::filter(Selection == region_) -> polygon_sub
          p = p + geom_polygon(data = polygon_sub, aes(x = X, y = Y), fill = NA, color = "black", linewidth = 0.3, inherit.aes = F)
        }
        
        return(p)
        
      }) %>% patchwork::wrap_plots(., ncol = 3, guides = "collect") -> p 
    
    SCN_core_shell_matrix_norm[gene_, ] %>% 
      t() %>% as.data.frame() %>% "colnames<-"(., "norm_expr") %>% 
      rownames_to_column("group") %>% 
      dplyr::mutate(region = gsub("(.+)_CT.+", "\\1", group)) %>% 
      dplyr::mutate(timepoint = gsub("(.+)_(CT-\\d+)-.+", "\\2", group)) -> df_plot
    df_plot$timepoint = factor(df_plot$timepoint, levels = c("CT-8", "CT-20"))
    df_plot %>% 
      ggplot(aes(x = region, y = norm_expr, color = timepoint, shape = timepoint)) + 
      geom_point(size = 3, alpha = 0.7) + 
      theme_classic() + 
      ylab("Normalized expression") -> p1
    df_spatiotemporal_all_1 %>% dplyr::filter(gene == gene_) -> df_
    df_[["CT8_cauchy_padj_core_vs_shell"]] %>% sprintf("%.2e", .) -> pval_core_shell_CT8
    df_[["CT20_cauchy_padj_core_vs_shell"]] %>% sprintf("%.2e", .) -> pval_core_shell_CT20
    p1 = p1 + 
      #      ggtitle(sprintf("DEG (SCN core vs SCN shell)\nCT-8 p.adj=%s\nCT-20 p.adj=%s", pval_core_shell_CT8, pval_core_shell_CT20))
      ggtitle("DEG (Core vs Shell)")
    
    df_count_mat_norm_list$SCN[gene_, ] %>% 
      t() %>% as.data.frame() %>% "colnames<-"(., "norm_expr") %>% 
      rownames_to_column("smaple") %>% 
      dplyr::mutate(timepoint = gsub("(CT-\\d+)-.+", "\\1", smaple)) -> df_plot
    df_plot$timepoint = factor(df_plot$timepoint, levels = c("CT-8", "CT-20"))
    
    df_plot %>% 
      ggplot(aes(x = timepoint, y = norm_expr, color = timepoint, shape = timepoint)) +
      geom_point(size = 3, alpha = 0.7) + 
      theme_classic() + 
      ylab("Normalized expression") -> p2
    p2 = p2 + ggtitle("DEG (CT-8 vs CT-20)")
    patchwork::wrap_plots(p1, p2, ncol = 1, guides = "collect") -> p3
    patchwork::wrap_plots(p, p3, widths = c(3,0.8), ncol = 2) -> p4
   
    df_spatiotemporal_all_2 %>% dplyr::filter(gene == gene_) -> df_
    df_[["CT-8_cauchy_padj_sparkx"]] %>% sprintf("%.2e", .) -> CT8_SVG_adj_pval
    df_[["CT-20_cauchy_padj_sparkx"]] %>% sprintf("%.2e", .) -> CT20_SVG_adj_pval
    df_[["p_adj.loc_TVG"]] %>% sprintf("%.2e", .) -> TVG_adj_pval
    df_[["spatiotemp.padj_spatiotemporal"]] %>% sprintf("%.2e", .) -> STVG_adj_pval
    title_ = sprintf("%s: SVG p.adj (CT-8)=%s, SVG p.adj (CT-20)=%s, TVG p.adj = %s, Spatiotemporal_p.adj=%s,\ncore_vs_shell p.adj = %s (CT-8), %s (CT-20)", gene_, CT8_SVG_adj_pval, CT20_SVG_adj_pval, TVG_adj_pval, STVG_adj_pval, pval_core_shell_CT8, pval_core_shell_CT20)
    p4 + patchwork::plot_annotation(title = title_)
  }) -> p_list ##4.3 spatiotemporal gene (show case) ----

ggsave(plot = p_list$Slc39a6, filename = "~/Downloads/Slc39a6.pdf", width = 21.22, height = 8.13)
ggsave(plot = p_list$Lhx1, filename = "~/Downloads/Lhx1.pdf", width = 21.22, height = 8.13)
ggsave(plot = p_list$Rgs16, filename = "~/Downloads/Rgs16.pdf", width = 21.22, height = 8.13)
ggsave(plot = p_list$Rpa1, filename = "~/Downloads/Rpa1.pdf", width = 21.22, height = 8.13)
ggsave(plot = p_list$Grin2c, filename = "~/Downloads/Grin2c.pdf", width = 21.22, height = 8.13)



c("Csnk1e", "Calb1", "Slc39a6") %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  purrr::map(function(gene_){
    norm_counts_matrix[gene_, ] %>% as.numeric() -> gene_expr
    upper = quantile(gene_expr, 0.99, na.rm = TRUE)
    c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(sample_){
        norm_counts_matrix[gene_, ] %>% 
          t() %>% 
          as.data.frame() %>% "colnames<-"(., "normalized_expr") %>% 
          rownames_to_column("cell") -> df_
        xy_coordinates[[sample_]] -> df_xy
        inner_join(x = df_, y = df_xy, by = "cell") -> df_plot
        
        print(upper)
        df_plot %>% 
          dplyr::arrange(normalized_expr) -> df_plot
        
#        df_plot %>% dplyr::mutate(normalized_expr = case_when(
#          cell %in% unique(c(scn_core_cells, scn_shell_cells)) ~ normalized_expr,
#          TRUE ~ 0
#        )) -> df_plot
        
        df_plot %>% dplyr::mutate(normalized_expr_1 = case_when(
          normalized_expr > upper ~ upper,
          TRUE ~ normalized_expr
        )) -> df_plot
        #        upper = quantile(df_plot$normalized_expr, 0.99, na.rm = TRUE)
        df_plot %>% 
          ggplot(aes(x = x, y = y, color = normalized_expr, size = normalized_expr_1)) + 
#          geom_point(size = 1) + 
          geom_point() + 
          scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, upper), oob = scales::squish) + 
          scale_size_continuous(range = c(1,4)) +
#          scale_size_continuous(limits = c(0,upper)) +
          theme_void() + 
          labs(color = "Normalized\nExpression") +
          ggtitle(sprintf("%s", sample_)) -> p
        
        polygon_[[sample_]]$Selection %>% unique() %>% {.[grepl("shell|core", .)]} -> regions_
        for (region_ in regions_){
          polygon_[[sample_]] %>% dplyr::filter(Selection == region_) -> polygon_sub
          p = p + geom_polygon(data = polygon_sub, aes(x = X, y = Y), fill = NA, color = "black", linewidth = 0.3, inherit.aes = F)
        }
       
        return(p)
        
      }) %>% patchwork::wrap_plots(., ncol = 3, guides = "collect") -> p 
    
#    SCN_core_shell_matrix_norm[gene_, ] %>% 
#      t() %>% as.data.frame() %>% "colnames<-"(., "norm_expr") %>% 
#      rownames_to_column("group") %>% 
#      dplyr::mutate(region = gsub("(.+)_CT.+", "\\1", group)) %>% 
#      dplyr::mutate(timepoint = gsub("(.+)_(CT-\\d+)-.+", "\\2", group)) -> df_plot
#    df_plot$timepoint = factor(df_plot$timepoint, levels = c("CT-8", "CT-20"))
#    df_plot %>% 
#      ggplot(aes(x = region, y = norm_expr, color = timepoint, shape = timepoint)) + 
#      geom_point(size = 3, alpha = 0.7) + 
#      theme_classic() + 
#      ylab("Normalized expression") -> p1
#    df_spatiotemporal_all_1 %>% dplyr::filter(gene == gene_) -> df_
#    df_[["CT8_cauchy_padj_core_vs_shell"]] %>% sprintf("%.2e", .) -> pval_core_shell_CT8
#    df_[["CT20_cauchy_padj_core_vs_shell"]] %>% sprintf("%.2e", .) -> pval_core_shell_CT20
#    p1 = p1 + 
#      #      ggtitle(sprintf("DEG (SCN core vs SCN shell)\nCT-8 p.adj=%s\nCT-20 p.adj=%s", pval_core_shell_CT8, pval_core_shell_CT20))
#      ggtitle("DEG (Core vs Shell)")
#    
#    df_count_mat_norm_list$SCN[gene_, ] %>% 
#      t() %>% as.data.frame() %>% "colnames<-"(., "norm_expr") %>% 
#      rownames_to_column("smaple") %>% 
#      dplyr::mutate(timepoint = gsub("(CT-\\d+)-.+", "\\1", smaple)) -> df_plot
#    df_plot$timepoint = factor(df_plot$timepoint, levels = c("CT-8", "CT-20"))
#    
#    df_plot %>% 
#      ggplot(aes(x = timepoint, y = norm_expr, color = timepoint, shape = timepoint)) +
#      geom_point(size = 3, alpha = 0.7) + 
#      theme_classic() + 
#      ylab("Normalized expression") -> p2
#    p2 = p2 + ggtitle("DEG (CT-8 vs CT-20)")
#    patchwork::wrap_plots(p1, p2, ncol = 1, guides = "collect") -> p3
#    patchwork::wrap_plots(p, p3, widths = c(3,0.8), ncol = 2) -> p4
#    
#    df_spatiotemporal_all_2 %>% dplyr::filter(gene == gene_) -> df_
#    df_[["CT-8_cauchy_padj_sparkx"]] %>% sprintf("%.2e", .) -> CT8_SVG_adj_pval
#    df_[["CT-20_cauchy_padj_sparkx"]] %>% sprintf("%.2e", .) -> CT20_SVG_adj_pval
#    df_[["p_adj.loc_TVG"]] %>% sprintf("%.2e", .) -> TVG_adj_pval
#    df_[["spatiotemp.padj_spatiotemporal"]] %>% sprintf("%.2e", .) -> STVG_adj_pval
#    title_ = sprintf("%s: SVG p.adj (CT-8)=%s, SVG p.adj (CT-20)=%s, TVG p.adj = %s, Spatiotemporal_p.adj=%s,\ncore_vs_shell p.adj = %s (CT-8), %s (CT-20)", gene_, CT8_SVG_adj_pval, CT20_SVG_adj_pval, TVG_adj_pval, STVG_adj_pval, pval_core_shell_CT8, pval_core_shell_CT20)
#    p4 + patchwork::plot_annotation(title = title_)
  }) -> p_list_1

p_list_1$Calb1

pdf("~/Downloads/SVG_map_plot.pdf", width = 12.54, height = 4.29)
p_list_1
dev.off()

pdf("~/Downloads/SVG_map_plot_v2.pdf", width = 14.64, height = 6.78)
p_list_1
dev.off()

c("Csnk1e", "Calb1", "Per2", "Lhx1", "Vip", "Pcsk1", "Slc39a6", "Rgs16", "Rpa1", "Grin2c") %>% 
  "names<-"(.,.) %>% 
#  .[3] %>% 
  purrr::map(function(gene_){
    SCN_core_shell_matrix_norm[gene_, ] %>% 
      pivot_longer(cols = everything(), names_to = "group", values_to = "norm_expr") %>% 
      dplyr::mutate(region = gsub("(.+)_CT.+", "\\1", group)) %>% 
      dplyr::mutate(CT = gsub("(.+)_(CT-\\d+?)-.+", "\\2", group)) %>% 
      dplyr::mutate(region = gsub("SCN_", "", region)) %>% 
      dplyr::mutate(CT = gsub("-", "", CT)) %>% 
      dplyr::mutate(group_1 = sprintf("%s_%s", CT, region)) -> df_plot
    
    df_plot$group_1 = factor(df_plot$group_1, levels = c("CT8_core", "CT8_shell", "CT20_core", "CT20_shell"))
    df_plot$CT = factor(df_plot$CT, levels = c("CT8", "CT20"))
    
    df_plot %>% 
      ggplot(aes(x = group_1, y = norm_expr, shape = region, color = CT)) + 
      geom_point(size = 4, alpha = 0.7) + 
      stat_summary(
        aes(group = group_1, color = CT),
        fun = mean,
        geom = "crossbar",
        width = 0.4,
#        color = "black",
        linewidth = 0.4
      ) +
      theme_classic() + 
      ggtitle(gene_)
  }) -> p_list_2

pdf("~/Downloads/tmp_plot.pdf", height = 3, width = 3.9)
p_list_2
dev.off()


#5. Geneerate SVG, TGV, STG master table ----
df_spatiotemporal_all_2 %>% 
  dplyr::select(gene, `CT-8_cauchy_pval_sparkx`, `CT-8_cauchy_padj_sparkx`, `CT-20_cauchy_pval_sparkx`, `CT-20_cauchy_padj_sparkx`) %>% 
  "colnames<-"(., c("gene", "SVG_CT8_pval", "SVG_CT8_padj", "SVG_CT20_pval", "SVG_CT20_padj")) -> df_SVG
df_spatiotemporal_all_2 %>% 
  dplyr::select(gene, p_val_TVG, p_adj_TVG, logFC_TVG) %>% 
  "colnames<-"(., c("gene", "TVG_pval", "TVG_padj", "TVG_logFC")) -> df_TVG

df_spatiotemporal_all_2 %>% 
  dplyr::select(gene, CT8_cauchy_pval_core_vs_shell, CT8_cauchy_padj_core_vs_shell, `CT-8_avg_log2FC_core_vs_shell`, CT20_cauchy_pval_core_vs_shell, CT20_cauchy_padj_core_vs_shell, `CT-20_avg_log2FC_core_vs_shell`) %>% 
  "colnames<-"(., c("gene", "core_vs_shell_CT8_pval", "core_vs_shell_CT8_padj", "core_vs_shell_CT8_logFC", "core_vs_shell_CT20_pval", "core_vs_shell_CT20_padj", "core_vs_shell_CT20_logFC")) -> df_core_vs_shell

colnames(df_spatiotemporal_all_2)
  
df_spatiotemporal_all_2 %>% 
  dplyr::select(gene, spatiotemp.pval_spatiotemporal, spatiotemp.padj_spatiotemporal) %>% 
  "colnames<-"(., c("gene", "spatiotemp_pval", "spatiotemp_padj")) -> df_spatiotemporal

list(
  df_SVG, df_TVG, df_core_vs_shell, df_spatiotemporal
) %>% purrr::reduce(., full_join, by = "gene") -> df_final

#muscat_res_list$cell$SCN$Astrocytes$ds_pb$data$Astrocytes$counts
load("~/Downloads/muscat_res_list.rda")
muscat_res_list$cell$SCN %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(list_, celltype_){
    list_$ds_pb$data[[celltype_]]$counts -> mat_
    (rownames(mat_) %in% df_spatiotemporal_all_2$gene) %>% all()
    (df_spatiotemporal_all_2$gene %in% rownames(mat_)) %>% all()
    
    mat_ %>% as.data.frame() %>% 
      rownames_to_column("gene") %>% 
      dplyr::filter(gene %in% df_spatiotemporal_all_2$gene) -> df_
    
    dim(df_)
    setdiff(df_spatiotemporal_all_2$gene, df_$gene) -> genes_1
    
    tibble(gene = genes_1, `CT-20-1` = 0, `CT-20-2` = 0, `CT-20-3` = 0, `CT-8-1` = 0, `CT-8-2` = 0, `CT-8-3` = 0) -> df_1
    
    rbind(
      df_, 
      df_1
    ) -> df_2
    
    df_2 %>% 
      as.data.frame() %>% 
      column_to_rownames("gene") %>% 
      dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6})) %>% 
      .[df_spatiotemporal_all_2$gene, ] %>% 
      rownames_to_column("gene") -> df_3
    
  }) -> df_SCN_celltype_count_mat_norm

load("~/Downloads/DEG_muscat_res_list.rda")
muscat_SCN_res$ds_pb$data$All_cells$counts %>% 
  {
    mat_ = .
    mat_ %>% as.data.frame() %>% 
      rownames_to_column("gene") %>% 
      dplyr::filter(gene %in% df_spatiotemporal_all_2$gene) -> df_
    
    dim(df_)
    df_ %>% column_to_rownames("gene") %>% .[df_spatiotemporal_all_2$gene, ] %>% 
      dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6})) %>% 
      rownames_to_column("gene")
  } -> df_SCN_count_mat_norm

load("~/Downloads/SCN_core_shell_matrix_norm.Rdata")
SCN_core_shell_matrix %>% 
  purrr::map(function(mat_){
    mat_ %>% as.data.frame() %>% 
      rownames_to_column("gene") %>% 
      dplyr::filter(gene %in% df_spatiotemporal_all_2$gene) -> df_
    
    dim(df_)
    df_ %>% column_to_rownames("gene") %>% .[df_spatiotemporal_all_2$gene, ] %>% 
      dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6})) %>% 
      rownames_to_column("gene")
  }) -> df_SCN_core_shell_count_mat_norm

library(openxlsx)
df_list = list(
  SVG_TVG_STG_table = df_final,
  SCN_GABA_expr = df_SCN_celltype_count_mat_norm$GABAergic,
  SCN_Astrocytes_expr = df_SCN_celltype_count_mat_norm$Astrocytes,
  SCN_Vascular_expr = df_SCN_celltype_count_mat_norm$Vascular,
  SCN_Oligo_expr = df_SCN_celltype_count_mat_norm$Oligodendrocytes,
  SCN_core_expr = df_SCN_core_shell_count_mat_norm$SCN_core,
  SCN_shell_expr = df_SCN_core_shell_count_mat_norm$SCN_shell
)
write.xlsx(df_list, "~/Dropbox/Spatial_transcriptomics_new/rda/SCN_gene_expr_master_table.xlsx", overwrite = TRUE)
save(df_spatiotemporal_all_2, file = "~/Downloads/df_spatiotemporal_all_2.rda")