library(tidyverse)
library(Seurat)
library(Matrix)
library(dplyr)
library(SingleCellExperiment)
library(muscat)

#load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
#load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

#1. import Xenium data and coordinates----
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

## 1.1 Exprot celltype for xenium viewer ----
#xenium_obj_filtered@meta.data %>% 
#  dplyr::filter(orig.ident == "CT-8-1") %>% 
#  rownames_to_column("cell_id") %>% dplyr::select(cell_id, predicted.id) %>%
#  dplyr::mutate(cell_id = gsub("CT.+_(.+)", "\\1", cell_id)) %>% 
#  "colnames<-"(., c("cell_id", "group")) %>% 
#  write_csv("~/Downloads/CT-8-1-cell-type-annotations.csv")
#
#xenium_obj_filtered@meta.data %>% 
#  dplyr::filter(orig.ident == "CT-8-2") %>% 
#  rownames_to_column("cell_id") %>% dplyr::select(cell_id, predicted.id) %>%
#  dplyr::mutate(cell_id = gsub("CT.+_(.+)", "\\1", cell_id)) %>% 
#  "colnames<-"(., c("cell_id", "group")) %>% 
#  write_csv("~/Downloads/CT-8-2-cell-type-annotations.csv")
#
#xenium_obj_filtered@meta.data %>% 
#  dplyr::filter(orig.ident == "CT-8-3") %>% 
#  rownames_to_column("cell_id") %>% dplyr::select(cell_id, predicted.id) %>%
#  dplyr::mutate(cell_id = gsub("CT.+_(.+)", "\\1", cell_id)) %>% 
#  "colnames<-"(., c("cell_id", "group")) %>% 
#  write_csv("~/Downloads/CT-8-3-cell-type-annotations.csv")
#
#xenium_obj_filtered@meta.data %>% 
#  dplyr::filter(orig.ident == "CT-20-1") %>% 
#  rownames_to_column("cell_id") %>% dplyr::select(cell_id, predicted.id) %>%
#  dplyr::mutate(cell_id = gsub("CT.+_(.+)", "\\1", cell_id)) %>% 
#  "colnames<-"(., c("cell_id", "group")) %>% 
#  write_csv("~/Downloads/CT-20-1-cell-type-annotations.csv")
#
#xenium_obj_filtered@meta.data %>% 
#  dplyr::filter(orig.ident == "CT-20-2") %>% 
#  rownames_to_column("cell_id") %>% dplyr::select(cell_id, predicted.id) %>%
#  dplyr::mutate(cell_id = gsub("CT.+_(.+)", "\\1", cell_id)) %>% 
#  "colnames<-"(., c("cell_id", "group")) %>% 
#  write_csv("~/Downloads/CT-20-2-cell-type-annotations.csv")
#
#xenium_obj_filtered@meta.data %>% 
#  dplyr::filter(orig.ident == "CT-20-3") %>% 
#  rownames_to_column("cell_id") %>% dplyr::select(cell_id, predicted.id) %>%
#  dplyr::mutate(cell_id = gsub("CT.+_(.+)", "\\1", cell_id)) %>% 
#  "colnames<-"(., c("cell_id", "group")) %>% 
#  write_csv("~/Downloads/CT-20-3-cell-type-annotations.csv")

xy_coordinates %>% 
  purrr::map(function(df_){
    df_
    df_ %>% dplyr::filter(cell %in% colnames(xenium_obj_filtered)) -> df_
  }) -> xy_coordinates

meta_data = xenium_obj_filtered@meta.data

## 1.2 Import SCN, PVN, SON contours ----
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
son_cells = select_cells(region__ = "SON")
pvn_cells = select_cells(region__ = "PVN_refined")

meta_data[scn_cells, ]$predicted.id -> names(scn_cells)
save(scn_cells, file = "~/Downloads/SCN_cells.rda")
unname(scn_cells) -> scn_cells

## 1.3 generate SCN, PVN, SON seurat object ----
c(scn_cells, scn_core_cells, scn_shell_cells, pvn_cells, son_cells) %>% unique() -> cell_selected
xenium_obj_filtered = subset(xenium_obj_filtered, cells = cell_selected)
gc()

do.call(rbind, xy_coordinates) -> xy_coordinates_1
xy_coordinates_1 %>% as_tibble() %>% as.data.frame() %>% column_to_rownames("cell") -> xy_coordinates_1

list(
  scn = scn_cells,
  scn_core = scn_core_cells,
  scn_shell = scn_shell_cells,
  pvn = pvn_cells,
  son = son_cells
) %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(cells_, regions_){
    data.frame(cell = cells_) %>% dplyr::mutate(region = regions_)
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) %>% 
  {
    df_ = .
    meta_data %>% rownames_to_column("cell") -> meta_data
    left_join(x = meta_data, y = df_, by = "cell") %>% 
      dplyr::mutate(region = ifelse(is.na(region), "other", region))
  } -> meta_data_1 #annotate cell region

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  #    .[1] %>% 
  purrr::map(function(sample_){
    meta_data_1 %>% dplyr::filter(orig.ident == sample_) %>% 
      dplyr::select(cell, region) -> df_
    left_join(x = df_, y = xy_coordinates_1 %>% rownames_to_column("cell"), by = "cell") -> df_
    unique(meta_data_1$region) %>% 
      purrr::map(function(region_){
        df_ %>% dplyr::filter(region == region_)
      }) %>% do.call(rbind, .) -> df_
    
    unique(meta_data_1$region) %>% 
      .[c(3,4,5,2,6,1)] %>% rev() -> levels_
    
    df_$region = factor(df_$region, levels = rev(levels_))
    df_ %>% 
      ggplot(aes(x = x, y = y, color = region)) + 
      ggrastr::geom_point_rast(size = 0.3, alpha = 0.5) + 
      scale_color_manual(values = c(
        "scn" = "yellow",
        "scn_shell" = "#FF5733",
        "scn_core" = "#33C1FF",
        "pvn" = "#75FF33",
        "son" = "#FF33D4",
        "other" = "gray80"
      )) + theme_classic() + 
      ggtitle(sample_)
  })  -> p_list #coronal slide highlighting SCN, SCN_core, SCN_shell, PVN and SON
dev.off()

pdf("~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_SON_PVN_cells_xy_coordinates_scatter_plots_v2_1.pdf", width = 8.86, height = 6)
p_list
dev.off()

list(
  SCN = scn_cells %>% unique(),
  SCN_core = scn_core_cells %>% unique(),
  SCN_shell = scn_shell_cells %>% unique(),
  PVN = pvn_cells %>% unique(),
  SON = son_cells %>% unique()
)  %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(cells_, cluster_){
    xenium_obj_filtered@meta.data[cells_, ] -> meta_data
    data.frame(cell_type = cluster_, nCount = meta_data$nCount_Xenium)
  }) %>% do.call(rbind, .) %>%
  "rownames<-"(., NULL) %>% 
  dplyr::mutate(cell_type = factor(cell_type, levels = rev(c("SCN", "SCN_core", "SCN_shell", "SON", "PVN")))) %>% 
  ggplot(aes(y = cell_type, x = nCount, fill = cell_type)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) + 
  theme_classic() + 
  theme(legend.position = "none") -> p_nCount
p_nCount + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_nCount_1

list(
  SCN = scn_cells %>% unique(),
  SCN_core = scn_core_cells %>% unique(),
  SCN_shell = scn_shell_cells %>% unique(),
  PVN = pvn_cells %>% unique(),
  SON = son_cells %>% unique()
)  %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(cells_, cluster_){
    xenium_obj_filtered@meta.data[cells_, ] -> meta_data
    data.frame(cell_type = cluster_, nFeature = meta_data$nFeature_Xenium)
  }) %>% do.call(rbind, .) %>%
  "rownames<-"(., NULL) %>% 
  dplyr::mutate(cell_type = factor(cell_type, levels = rev(c("SCN", "SCN_core", "SCN_shell", "SON", "PVN")))) %>% 
  ggplot(aes(y = cell_type, x = nFeature, fill = cell_type)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) + 
  theme_classic() + 
  theme(legend.position = "none") -> p_nFeature
p_nFeature + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_nFeature_1

#genes_selected = c("Vip", "Avp", "Vipr2", "Prokr2", "Nms", "Grp", "Slc32a1", "Oxt", "Trh", "Sim1", "Crh")
genes_selected = c("Vip", "Vipr2", "Prokr2", "Nms", "Grp", "Slc32a1", "Oxt", "Trh", "Sim1", "Crh", "Avp")
mat_selected = GetAssayData(object = xenium_obj_filtered, assay = "Xenium", layer = "counts")[genes_selected, ]
dim(mat_selected)

list(
  SCN = scn_cells %>% unique(),
  SCN_core = scn_core_cells %>% unique(),
  SCN_shell = scn_shell_cells %>% unique(),
  PVN = pvn_cells %>% unique(),
  SON = son_cells %>% unique()
)  %>% 
  purrr::map(function(cell_){
    xenium_obj_filtered@meta.data[cell_, ] -> df_
    df_$nCount_Xenium %>% sum()
  }) %>% unlist() %>% median() -> median_

list(
  SCN = scn_cells %>% unique(),
  SCN_core = scn_core_cells %>% unique(),
  SCN_shell = scn_shell_cells %>% unique(),
  PVN = pvn_cells %>% unique(),
  SON = son_cells %>% unique()
)  %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(cell_, cluster_){
    xenium_obj_filtered@meta.data[cell_, ] -> df_
    df_$nCount_Xenium %>% sum() -> lib_size
    cells_selected = rownames(df_)
    mat_selected[, cells_selected] %>% rowSums() %>% {(./lib_size)*median_} %>% 
      as.data.frame() %>% "colnames<-"(., "expr") %>% rownames_to_column("gene") %>% 
      dplyr::mutate(cell_type = cluster_, .before = 1) -> df_1
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_expr

list(
  SCN = c(scn_core_cells, scn_shell_cells) %>% unique(),
  SCN_core = scn_core_cells %>% unique(),
  SCN_shell = scn_shell_cells %>% unique(),
  PVN = pvn_cells %>% unique(),
  SON = son_cells %>% unique()
)  %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(cell_, cluster_){
    xenium_obj_filtered@meta.data[cell_, ] -> df_
    cells_selected = rownames(df_)
    n_cells = nrow(df_)
    mat_ = mat_selected[, cells_selected] 
    (mat_ != 0) %>% rowSums() %>% {./n_cells} %>% 
      as.data.frame() %>% "colnames<-"(., "non_zero_prop") %>% 
      rownames_to_column("gene") %>% 
      dplyr::mutate(cell_type = cluster_, .before = 1) -> df_1
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_non_zero_prop

full_join(x = df_expr, y = df_non_zero_prop, by = c("cell_type", "gene")) -> df_final

#unique(df_final$cell_type) %>% 
#  "names<-"(.,.) %>% 
#  purrr::map(function(cell_type_){
#    df_final %>% dplyr::filter(cell_type == cell_type_) -> df_
#    #    df_ %>% dplyr::mutate(scaled_expr = log2(expr+1)) -> df_
#    #    df_ %>% dplyr::mutate(scaled_expr = log1p(expr)) -> df_
#    df_ %>% 
#      dplyr::mutate(expr = case_when(
#        expr > quantile(expr, 0.99) ~ max(expr[expr < quantile(expr, 0.99)]),
#        #        expr > quantile(expr, 0.95) ~ quantile(expr, 0.95),
#        TRUE ~ expr
#      )) -> df_
#    df_ %>% 
#      dplyr::mutate(scaled_expr = scales::rescale(expr, to = c(0,1))) -> df_
#    df_ %>% 
#      dplyr::mutate(non_zero_prop = non_zero_prop*100)
#  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_plot

unique(df_final$gene) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    df_final %>% dplyr::filter(gene == gene_) -> df_
    #        df_ %>% dplyr::mutate(scaled_expr = log2(expr+1)) -> df_
    #         df_ %>% dplyr::mutate(scaled_expr = log1p(expr)) -> df_
    df_ %>% 
      dplyr::mutate(expr = case_when(
        expr > quantile(expr, 0.99) ~ max(expr[expr < quantile(expr, 0.99)]),
        #        expr > quantile(expr, 0.95) ~ quantile(expr, 0.95),
        TRUE ~ expr
      )) -> df_
    df_ %>% 
      dplyr::mutate(scaled_expr = scales::rescale(expr, to = c(0,1))) -> df_
    df_ %>% 
      dplyr::mutate(non_zero_prop = non_zero_prop*100)
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_plot

df_plot$cell_type = factor(df_plot$cell_type, levels = rev(c("SCN", "SCN_core", "SCN_shell", "SON", "PVN")))
df_plot$gene = factor(df_plot$gene, levels = genes_selected)
df_plot %>% 
  ggplot(aes(x = gene, y = cell_type)) +
  geom_point(aes(size = non_zero_prop, color = scaled_expr)) + 
  scale_size(range = c(0,15), limits = c(0,100)) + 
  scale_color_gradient(low = "grey", high = "blue") + 
  theme_classic() -> p_dotplot
p_dotplot + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_dotplot_1 #SCN, SCN_core, SCN_shell, SON, PVN marker genes annotation
gc()

patchwork::wrap_plots(p_nCount, p_nFeature_1, p_dotplot_1, nrow = 1) + 
  patchwork::plot_layout(widths = c(1,1,3)) -> p_dotplot_2
ggsave(plot = p_dotplot_2, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/regional_marker_genes_dotplot_v2.pdf", 
       width = 12.18, height = 4.99) #Figure 1E

list(
  SCN = scn_cells %>% unique(),
  SCN_core = scn_core_cells %>% unique(),
  SCN_shell = scn_shell_cells %>% unique(),
  PVN = pvn_cells %>% unique(),
  SON = son_cells %>% unique()
)  %>% 
  #    .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(cells_, region_){
    xenium_obj_filtered@meta.data[cells_, ] %>% {table(.$orig.ident, .$cluster)} %>% 
      apply(2, sd) %>% as.data.frame() %>% 
      "colnames<-"(., "sd") %>% 
      rownames_to_column("celltype") -> df_sd
    xenium_obj_filtered@meta.data[cells_, ] %>% .$cluster %>% table() %>% 
      as.data.frame() %>% "colnames<-"(., c("celltype", "n_cells")) %>% 
      dplyr::mutate(percent = (n_cells/sum(n_cells))*100) %>% 
      dplyr::arrange(desc(n_cells)) -> df_
    
    df_ %>% 
      dplyr::mutate(celltype = factor(celltype, levels = rev(c("1", "4", "10", "2", "7", "9", "13", "3", "5", "8", "6", "11", "12")))) -> df_
    
    setdiff(x = c(1:13), y = df_$celltype) %>% as.data.frame() %>% 
      "colnames<-"(., "celltype") %>% dplyr::mutate(n_cells = 0, percent = 0) -> df_1
    rbind(df_, df_1) -> df_
    
    df_ = left_join(x = df_, y = df_sd, by = "celltype")
    df_ %>% dplyr::mutate(sd = case_when(
      is.na(sd) ~ 0,
      TRUE ~ sd
    )) -> df_
    
    df_ %>% 
      dplyr::mutate(celltype = factor(celltype, levels = rev(c("1", "4", "10", "2", "7", "9", "13", "3", "5", "8", "6", "11", "12")))) -> df_
    
    df_ %>% 
      ggplot(aes(y = celltype, x = n_cells, fill = celltype)) + 
      geom_bar(stat="identity") + 
      geom_errorbar(aes(xmin = n_cells - sd, xmax = n_cells + sd), width = 0.2) + 
      scale_fill_manual(values = c("1" = "#00bf7d", "4" = "#00bf7d", "10" = "#00bf7d", "2" = "#f8766d", "7" = "#f8766d", "9" = "#f8766d"
                                   , "13" = "#f8766d", "3" = "#a3a500", "5" = "#a3a500", "8" = "#a3a500", "6" = "#e76bf3", 
                                   "11" = "#00b0f6", "12" = "#00b0f6")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ggrepel::geom_text_repel(aes(label = paste0(round(percent,2), "%")), 
                               position = position_stack(vjust = 1), color = "black") -> p
    
    p + ggtitle(region_)
  }) -> p_list

patchwork::wrap_plots(p_list, nrow = 1, guides = "collect") -> p_composition #illustrate cell-type composition in SCN, SCN_cell, SCN_core, PVN and SON
ggsave(plot = p_composition, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/regional_marker_genes_cell_composition_v2.pdf", 
       width = 15.42, height = 4.99)

#2. Identify TVGs in SCN, PVN and SON ----
library(muscat)
run_muscat = function(cell_ids, sc_obj, celltype, assay){
  sce <- as.SingleCellExperiment(sc_obj, assay = assay)
  sce = sce[,colnames(sce) %in% cell_ids]
  gc()
  
  sample_id = colData(sce)$orig.ident %>% as.character()
  group_id = gsub("(CT-\\d+?)-.+", "\\1", sample_id) %>% gsub("-", "", .)
  
  colData(sce)[["sample_id"]] <- sample_id
  colData(sce)[["group_id"]] <- group_id
  colData(sce)[["cluster_id"]] <- celltype
  
  sce <- prepSCE(
    sce,
    kid = "cluster_id",   # your column with cluster/cell-type labels
    sid = "sample_id",    # sample IDs: CT-8-1 ... CT-20-3
    gid = "group_id"      # group labels: CT8 / CT20
  )
  
  pb <- aggregateData(sce, 
                      assay = "counts", fun = "sum",
                      by = c("cluster_id", "sample_id"))
  
  metadata(pb)$experiment_info
  
  design_ <- model.matrix(~ 0 + group_id, data = metadata(pb)$experiment_info)
  colnames(design_) = c("CT20", "CT8")
  metadata(pb)$experiment_info$sample_id %>% as.character() -> rownames(design_)
  design_[colnames(pb), ] -> design_
  contrast_ <- limma::makeContrasts(CT20 - CT8, levels = design_)
  
  ds_pb <- pbDS(pb, method = "edgeR", design = design_, contrast = contrast_, filter = "none", min_cells = 3)
  gc()
  
  pb_binarized <- aggregateData(sce, 
                                assay = "counts", 
                                fun = "num.detected", # !!!
                                by = c("cluster_id", "sample_id"))
  dd_pb <- pbDD(pb_binarized, design = design_, contrast = contrast_, filter = "none", min_cells = 3)
  res = list(ds_pb = ds_pb, dd_pb = dd_pb)
  
  gc()
  return(res)
}

xenium_obj_filtered@meta.data %>% 
  dplyr::mutate(custom_class = case_when(
    cluster == "3" ~ "Astrocytes",
    cluster == "5" ~ "Oligodendrocytes",
    cluster == "8" ~ "Ependymal",
    TRUE ~ predicted.id
  )) -> meta_data_2

cell_types_ = c("Astrocytes", "GABAergic", "Glutamatergic", "Immune", "Oligodendrocytes", "Vascular", "Ependymal")
list(
  SCN = scn_cells %>% unique(),
  SCN_core = scn_core_cells %>% unique(),
  SCN_shell = scn_shell_cells %>% unique(),
  PVN = pvn_cells %>% unique(),
  SON = son_cells %>% unique()
) %>% 
  purrr::map(function(cells_){
    meta_data_2[cells_, ] -> df_
    df_ %>% dplyr::filter(custom_class %in% cell_types_) -> df_
    table(df_$orig.ident, df_$custom_class) %>% 
      apply(., 1, function(x){x>=3}) %>% 
      rowSums() %>% .[. == 6] %>% names() -> cell_types_filtered_
    cell_types_filtered_ %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(custom_class_){
        df_ %>% dplyr::filter(custom_class == custom_class_) %>% rownames() -> cell_ids_
      })
  }) -> cell_list

cell_list %>% 
  purrr::map2(.x=., .y=names(.), .f=function(list_, region_){
    list_ %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(cells_, celltype_){
        muscat_res = run_muscat(cell_ids = cells_, sc_obj = xenium_obj_filtered, celltype = celltype_, assay = "Xenium")
      })
  }) -> muscat_res_list ##2.1 muscat results (Whole cell) ----

cell_list %>% 
  purrr::map2(.x=., .y=names(.), .f=function(list_, region_){
    list_ %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(cells_, celltype_){
        muscat_res = run_muscat(cell_ids = cells_, sc_obj = xenium_obj_filtered, celltype = celltype_, assay = "Xenium_nuc")
      })
  }) -> muscat_res_nuc_list ##2.2 muscat results (Nucleus) ----

cell_raw_counts = GetAssayData(xenium_obj_filtered, assay = "Xenium", layer = "counts")
nuc_raw_counts = GetAssayData(xenium_obj_filtered, assay = "Xenium_nuc", layer = "counts")

cyto_raw_counts = cell_raw_counts - nuc_raw_counts
rm(cell_raw_counts, nuc_raw_counts)

cyto_raw_counts = CreateAssayObject(counts = cyto_raw_counts)
xenium_obj_filtered[["Xenium_cyto"]] <- cyto_raw_counts
rm(cyto_raw_counts)
gc()

cell_list %>% 
  purrr::map2(.x=., .y=names(.), .f=function(list_, region_){
    list_ %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(cells_, celltype_){
        muscat_res = run_muscat(cell_ids = cells_, sc_obj = xenium_obj_filtered, celltype = celltype_, assay = "Xenium_cyto")
      })
  }) -> muscat_res_cyto_list ##2.3 muscat results (Cytoplasm) ----

list(cell = muscat_res_list,
     cyto = muscat_res_cyto_list, 
     nuc = muscat_res_nuc_list) %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(muscat_res_list_, compartment_){
    muscat_res_list_ %>% 
      purrr::map2(.x=., .y=names(.), .f=function(list_, region_){
        list_ %>% 
          purrr::map2(.x=.,.y=names(.),.f=function(list__, celltype_){
            list__[["ds_pb"]][["table"]][[1]][[celltype_]] -> DEG_df
            DEG_df %>% dplyr::mutate(region = region_, compartment = compartment_, .after = 1) -> DEG_df
          }) %>% do.call(rbind, .) -> df_
        df_ %>% dplyr::mutate(p_adj.glb = p.adjust(p_val, method = "BH")) -> DEG_df
      }) %>% do.call(rbind, .) -> DEG_df 
  }) %>% do.call(rbind, .) -> DEG_df ##2.4 TVGs result(region (SCN, SON, PVN), compartment (cell, nuc, cyto), celltype)----
rownames(DEG_df) <- NULL
save(DEG_df, file = "~/Downloads/DEG_df.Rdata")

list(cell = muscat_res_list,
     cyto = muscat_res_cyto_list, 
     nuc = muscat_res_nuc_list) -> muscat_res_list

save(muscat_res_list, file="~/Downloads/muscat_res_list.rda")

Clock_genes_1 = c("Arntl", "Bhlhe40", "Bhlhe41", "Clock", "Npas2", "Dbp", "Nfil3", "Nr1d1", "Nr1d2", "Rora", "Rorb", "Rorc", "Cry1", "Cry2", "Ciart", "Per1", "Per2", "Per3", "Hlf", "Tef", "Csnk1d", "Csnk1e", "Rgs16", "Sirt2")
Clock_genes_1 %in% rownames(xenium_obj_filtered)
DEG_df %>% dplyr::filter(gene %in% Clock_genes_1) %>% 
  dplyr::mutate(group = sprintf("%s_%s", region, cluster_id), .after = 1) %>% 
  {
    df_ = .
    Clock_genes_1 %>% 
      "names<-"(.,.) %>% 
      #      .[1] %>% 
      purrr::map(function(gene_){
        df_ %>% dplyr::filter(gene == gene_) -> df_
        df_ %>% dplyr::arrange(region) -> df_
        df_ %>% dplyr::mutate(Sig = case_when(
          (p_adj.loc < 0.05)&(logFC > log2(1.1)) ~ "CT20_enriched",
          (p_adj.loc < 0.05)&(logFC < -log2(1.1)) ~ "CT8_enriched",
          TRUE ~ "Not_significant"
        )) -> df_plot
        df_plot$group %>% unique() -> levels_
        c("^SCN_[A-Za-z]+$", "SCN_core_", "SCN_shell_", "SON_", "PVN_") %>% 
          #        c("^SCN_[A-Za-z]+$") %>% 
          purrr::map(function(order_){
            levels_[grepl(order_, levels_)]
          }) %>% purrr::reduce(c) -> levels_
        df_plot$group = factor(df_plot$group, levels = rev(levels_))
        df_plot$Sig = factor(df_plot$Sig, levels = c("CT8_enriched", "CT20_enriched", "Not_significant"))
        df_plot %>% 
          ggplot(aes(y = group, x = compartment, fill = Sig)) + 
          geom_tile(color = "black") + 
          scale_fill_manual(values = c(CT8_enriched = "#f8766d", CT20_enriched = "#00bfc4", Not_significant = "lightgrey")) +
          ylab(NULL) +
          xlab(NULL) +
          theme_minimal() + 
          ggtitle(gene_)
      })
  } -> p_list_p_val_adj_local

patchwork::wrap_plots(p_list_p_val_adj_local, ncol = 4, guides = "collect") -> p_p_val_adj_local

#Clock_genes = c("Arntl", "Bhlhe40", "Bhlhe41", "Clock", "Npas2", "Dbp", "Nfil3", "Nr1d1", "Nr1d2", "Rora", "Rorb", "Rorc", "Cry1", "Cry2", "Ciart", "Per1", "Per2", "Per3", "Hlf", "Tef", "Csnk1d", "Csnk1e", "Rgs16", "Sirt2")
#3. Identify cell-type specific TVGs of the entire brain section ----
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")

xenium_obj_filtered@meta.data %>% 
  dplyr::mutate(predicted.id_1 = case_when(
    cluster == "3" ~ "Astrocytes", 
    cluster == "5" ~ "Oligodendrocytes",
    cluster == "8" ~ "Ependymal",
    TRUE ~ predicted.id
  )) -> meta_data_3

unique(meta_data_3$predicted.id_1) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(predicted.id_1_){
    meta_data_3 %>% dplyr::filter(predicted.id_1 == predicted.id_1_) -> df_
    cells = rownames(df_)
  }) -> cell_list_2

cell_list_2 %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(cells_, celltype_){
    muscat_res = run_muscat(cell_ids = cells_, sc_obj = xenium_obj_filtered, celltype = celltype_, assay = "Xenium")
  }) -> muscat_res_celltype_list_1

muscat_res_celltype_list_1 %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(list__, celltype_){
    list__[["ds_pb"]][["table"]][[1]][[celltype_]] -> DEG_df
  }) %>% do.call(rbind, .) %>% 
  dplyr::mutate(p_adj.glb = p.adjust(p_val, method = "BH")) %>% 
  "rownames<-"(., NULL) -> DEG_df_celltype_1 ##3.1 muscat results (Whole brain cellt-ype)----

#4. Identify TVGs in SCN (All cells) ----    
scn_cells %>% length() #4915 SCN cells
muscat_SCN_res = run_muscat(cell_ids = scn_cells, sc_obj = xenium_obj_filtered, celltype = "All_cells", assay = "Xenium")
muscat_SCN_res$ds_pb$table$`CT20 - CT8`$All_cells %>% 
  dplyr::mutate(group = "SCN", region = "SCN", compartment = "cell", .after = 1) %>%
  dplyr::mutate(Sig = case_when(
    (p_adj.loc < 0.05)&(logFC > log2(1.1)) ~ "CT20_enriched",
    (p_adj.loc < 0.05)&(logFC < -log2(1.1)) ~ "CT8_enriched",
    TRUE ~ "Not_significant"
  )) -> DEG_df_SCN

p_list_p_val_adj_local %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(p_,gene_){
    p_$data %>% dplyr::filter(region == "SCN", compartment == "cell") %>% 
#      dplyr::filter(grepl("GABA|Astro", group)) -> df_
      dplyr::filter(grepl("GABA", group)) -> df_
#    DEG_df_celltype_1 %>% dplyr::filter(gene == gene_, cluster_id %in% c("GABAergic", "Astrocytes")) %>% 
    DEG_df_celltype_1 %>% dplyr::filter(gene == gene_, cluster_id %in% c("GABAergic")) %>% 
      dplyr::mutate(group = sprintf("All_%s", cluster_id), region = "All", compartment = "cell", .after = 1) %>% 
      dplyr::mutate(Sig = case_when(
        (p_adj.loc < 0.05)&(logFC > log2(1.1)) ~ "CT20_enriched",
        (p_adj.loc < 0.05)&(logFC < -log2(1.1)) ~ "CT8_enriched",
        TRUE ~ "Not_significant"
      )) -> df_1
    DEG_df_SCN %>% dplyr::filter(gene == gene_) -> df_2
    rbind(df_, df_1, df_2) -> df_final
  }) %>% do.call(rbind, .) %>% 
  dplyr::mutate(log_p_adj = -log10(p_adj.loc)) %>% 
  dplyr::mutate(log_p_adj = case_when(
    Sig == "CT20_enriched" ~ log_p_adj,
    Sig == "CT8_enriched" ~ -log_p_adj,
    #    TRUE ~ NA
    TRUE ~ log_p_adj
  )) -> df_plot

df_plot %>% 
  dplyr::mutate(gene = factor(gene, levels = Clock_genes_1)) %>% 
  dplyr::mutate(Sig_1 = case_when(
    p_adj.loc < 0.001 ~ "***",
    p_adj.loc < 0.01 ~ "**",
    p_adj.loc < 0.05 ~ "*",
    TRUE ~ ""
  )) %>% 
#  dplyr::mutate(group = factor(group, levels = rev(c("SCN", "SCN_GABAergic", "All_GABAergic", "SCN_Astrocytes", "All_Astrocytes")))) %>% 
  dplyr::mutate(group = factor(group, levels = rev(c("All_GABAergic", "SCN", "SCN_GABAergic")))) %>% 
  ggplot(aes(x = gene , y = group, fill = log_p_adj)) + 
  geom_tile(color = "black") + 
  scale_fill_gradient2(low = "red", high = "blue", mid = "grey", midpoint = 0) + 
  geom_text(aes(label = Sig_1), size = 4, color = "black") + 
  theme_classic() -> p1
ggsave(plot = p1, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_clock_genes_heatmap_v2.pdf", width = 14.81, height = 2.01) #Figure 2C

DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic", gene %in% Clock_genes_1) %>% 
  dplyr::mutate(Sig = case_when(
    (p_adj.loc < 0.05)&(logFC > log2(1.1)) ~ "CT20_enriched",
    (p_adj.loc < 0.05)&(logFC < -log2(1.1)) ~ "CT8_enriched",
    TRUE ~ "Not_significant"
  )) %>% 
  dplyr::mutate(log_p_adj = -log10(p_adj.loc)) %>% 
  dplyr::mutate(log_p_adj = case_when(
    Sig == "CT20_enriched" ~ log_p_adj,
    Sig == "CT8_enriched" ~ -log_p_adj,
    #    TRUE ~ NA
    TRUE ~ log_p_adj
  )) -> df_plot

df_plot %>% 
  dplyr::mutate(gene = factor(gene, levels = Clock_genes_1)) %>% 
  dplyr::mutate(Sig_1 = case_when(
    p_adj.loc < 0.001 ~ "***",
    p_adj.loc < 0.01 ~ "**",
    p_adj.loc < 0.05 ~ "*",
    TRUE ~ ""
  )) %>% 
  #  dplyr::mutate(group = factor(group, levels = rev(c("SCN", "SCN_GABAergic", "All_GABAergic", "SCN_Astrocytes", "All_Astrocytes")))) %>% 
  dplyr::mutate(compartment = factor(compartment, levels = rev(c("cell", "nuc", "cyto")))) %>% 
  ggplot(aes(x = gene , y = compartment, fill = log_p_adj)) + 
  geom_tile(color = "black") + 
  scale_fill_gradient2(low = "red", high = "blue", mid = "grey", midpoint = 0) + 
  geom_text(aes(label = Sig_1), size = 4, color = "black") + 
  theme_classic() -> p1_1

#5. Identify TVGs in SON (All cells) ----    
son_cells
muscat_SON_res = run_muscat(cell_ids = son_cells, sc_obj = xenium_obj_filtered, celltype = "All_cells", assay = "Xenium")
muscat_SON_res$ds_pb$table$`CT20 - CT8`$All_cells %>% 
  dplyr::mutate(group = "SON", region = "SON", compartment = "cell", .after = 1) %>%
  dplyr::mutate(Sig = case_when(
    (p_adj.loc < 0.05)&(logFC > log2(1.1)) ~ "CT20_enriched",
    (p_adj.loc < 0.05)&(logFC < -log2(1.1)) ~ "CT8_enriched",
    TRUE ~ "Not_significant"
  )) -> DEG_df_SON

p_list_p_val_adj_local %>% 
  #    .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(p_,gene_){
    p_$data %>% dplyr::filter(region == "SON", compartment == "cell") %>% 
      dplyr::filter(grepl("Gluta|Astro", group)) -> df_
    DEG_df_celltype_1 %>% dplyr::filter(gene == gene_, cluster_id %in% c("Glutamatergic", "Astrocytes")) %>% 
      dplyr::mutate(group = sprintf("All_%s", cluster_id), region = "All", compartment = "cell", .after = 1) %>% 
      dplyr::mutate(Sig = case_when(
        (p_adj.loc < 0.05)&(logFC > log2(1.1)) ~ "CT20_enriched",
        (p_adj.loc < 0.05)&(logFC < -log2(1.1)) ~ "CT8_enriched",
        TRUE ~ "Not_significant"
      )) -> df_1
    DEG_df_SON %>% dplyr::filter(gene == gene_) -> df_2
    rbind(df_, df_1, df_2) -> df_final
  }) %>% do.call(rbind, .) %>% 
  dplyr::mutate(log_p_adj = -log10(p_adj.loc)) %>% 
  dplyr::mutate(log_p_adj = case_when(
    Sig == "CT20_enriched" ~ log_p_adj,
    Sig == "CT8_enriched" ~ -log_p_adj,
    #    TRUE ~ NA
    TRUE ~ log_p_adj
  )) -> df_plot

df_plot %>% 
  dplyr::mutate(gene = factor(gene, levels = Clock_genes_1)) %>% 
  dplyr::mutate(Sig_1 = case_when(
    p_adj.loc < 0.001 ~ "***",
    p_adj.loc < 0.01 ~ "**",
    p_adj.loc < 0.05 ~ "*",
    TRUE ~ ""
  )) %>% 
  {
    df_ = .
    df_$group = factor(df_$group, levels = rev(c("SON", "SON_Glutamatergic", "All_Glutamatergic", "SON_Astrocytes", "All_Astrocytes", 
                                                 "SON_Oligodendrocytes", "All_Oligodendrocytes", "SON_Vascular", "All_Vascular")))
    df_
  } %>% 
  ggplot(aes(x = gene , y = group, fill = log_p_adj)) + 
  geom_tile(color = "black") + 
  scale_fill_gradient2(low = "red", high = "blue", mid = "grey", midpoint = 0) + 
  geom_text(aes(label = Sig_1), size = 4, color = "black") + 
  theme_classic() -> p2
ggsave(plot = p2, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SON_clock_genes_heatmap_v2.pdf", width = 11.77, height = 2.48)

#5. Identify TVGs in PVN (All cells) ----    
muscat_PVN_res = run_muscat(cell_ids = pvn_cells, sc_obj = xenium_obj_filtered, celltype = "All_cells", assay = "Xenium")
muscat_PVN_res$ds_pb$table$`CT20 - CT8`$All_cells %>% 
  dplyr::mutate(group = "PVN", region = "PVN", compartment = "cell", .after = 1) %>%
  dplyr::mutate(Sig = case_when(
    (p_adj.loc < 0.05)&(logFC > log2(1.1)) ~ "CT20_enriched",
    (p_adj.loc < 0.05)&(logFC < -log2(1.1)) ~ "CT8_enriched",
    TRUE ~ "Not_significant"
  )) -> DEG_df_PVN

p_list_p_val_adj_local %>% 
  #    .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(p_,gene_){
    p_$data %>% dplyr::filter(region == "PVN", compartment == "cell") %>% 
      dplyr::filter(grepl("GABA|Astro|Oligo", group)) -> df_
    DEG_df_celltype_1 %>% dplyr::filter(gene == gene_, cluster_id %in% c("GABAergic", "Astrocytes", "Oligodendrocytes")) %>% 
      dplyr::mutate(group = sprintf("All_%s", cluster_id), region = "All", compartment = "cell", .after = 1) %>% 
      dplyr::mutate(Sig = case_when(
        (p_adj.loc < 0.05)&(logFC > log2(1.1)) ~ "CT20_enriched",
        (p_adj.loc < 0.05)&(logFC < -log2(1.1)) ~ "CT8_enriched",
        TRUE ~ "Not_significant"
      )) -> df_1
    DEG_df_PVN %>% dplyr::filter(gene == gene_) -> df_2
    rbind(df_, df_1, df_2) -> df_final
  }) %>% do.call(rbind, .) %>% 
  dplyr::mutate(log_p_adj = -log10(p_adj.loc)) %>% 
  dplyr::mutate(log_p_adj = case_when(
    Sig == "CT20_enriched" ~ log_p_adj,
    Sig == "CT8_enriched" ~ -log_p_adj,
    #    TRUE ~ NA
    TRUE ~ log_p_adj
  )) -> df_plot

df_plot %>% 
  dplyr::mutate(gene = factor(gene, levels = Clock_genes_1)) %>% 
  dplyr::mutate(Sig_1 = case_when(
    p_adj.loc < 0.001 ~ "***",
    p_adj.loc < 0.01 ~ "**",
    p_adj.loc < 0.05 ~ "*",
    TRUE ~ ""
  )) %>% 
  {
    df_ = .
    df_$group = factor(df_$group, levels = rev(c("PVN", "PVN_GABAergic", "All_GABAergic", "PVN_Glutamatergic", "All_Glutamatergic", "PVN_Astrocytes", "All_Astrocytes", 
                                                 "PVN_Oligodendrocytes", "All_Oligodendrocytes", "PVN_Ependymal", "All_Ependymal", "PVN_Vascular", "All_Vascular")))
    df_
  } %>% 
  ggplot(aes(x = gene , y = group, fill = log_p_adj)) + 
  geom_tile(color = "black") + 
  scale_fill_gradient2(low = "red", high = "blue", mid = "grey", midpoint = 0) + 
  geom_text(aes(label = Sig_1), size = 4, color = "black") + 
  theme_classic() -> p3
ggsave(plot = p3, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/PVN_clock_genes_heatmap_v2.pdf", width = 13.48, height = 3.41)

save(DEG_df_celltype_1, DEG_df_SCN, DEG_df_PVN, DEG_df_SON, file = "~/Downloads/DEG_dfs.rda")
save(muscat_res_celltype_list_1, muscat_SCN_res, muscat_PVN_res, muscat_SON_res, file = "~/Downloads/DEG_muscat_res_list.rda")
rm(list=ls())
gc()

# Plot number of celltype specific TVGs in SCN, SON and PVN ----
unique(DEG_df$region) %>%
  "names<-"(.,.) %>% 
  purrr::map(function(region_){
    DEG_df %>% dplyr::filter(region == region_, compartment == "cell") -> df_
    unique(df_$cluster_id) %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(celltype_){
        df_ %>% dplyr::filter(cluster_id == celltype_) -> df_
        df_ %>% dplyr::filter((p_adj.loc < 0.05)&(abs(logFC) > log2(1.1))) -> df_
        data.frame(region = region_, celltype = celltype_, n_DEGs = nrow(df_))
      }) %>% do.call(rbind, .) -> df_
    df_ %>% dplyr::filter(n_DEGs > 0) -> df_
    df_ %>% dplyr::arrange(desc(n_DEGs)) -> df_
    df_$celltype = factor(df_$celltype, levels = rev(df_$celltype))
    df_ %>% 
      ggplot(aes(y = celltype, x = n_DEGs, fill = celltype)) + 
      geom_bar(stat = "identity") + 
      geom_text(aes(label = n_DEGs), 
                hjust = -0.1, size = 5) + 
      scale_fill_manual(values = c("GABAergic" = "#f8766d", "Glutamatergic" = "#00bf7d", "Astrocytes" = "#a3a500", "Oligodendrocytes" = "#a3e500", "Ependymal" = "#c3d846", "Vascular" = "#e76bf3", "Immune" = "#00b0f6")) +
      ggtitle(region_) + 
      theme_classic()
  }) %>% patchwork::wrap_plots(ncol = 3, guides = "collect") -> p

ggsave("~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_SON_PVN_DEG_barplots.pdf", plot = p, width = 12.52, height = 6.88)

gene_list_for_GO = list()
unique(DEG_df$region) %>%
  "names<-"(.,.) %>% 
  purrr::map(function(region_){
    DEG_df %>% dplyr::filter(region == region_, compartment == "cell") -> df_
    unique(df_$cluster_id) %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(celltype_){
        df_ %>% dplyr::filter(cluster_id == celltype_) -> df_
        df_ %>% dplyr::filter((p_adj.loc < 0.05)&(abs(logFC) > log2(1.1))) -> df_
        nrow(df_) -> order_
        sum(df_$logFC > 0) -> CT20_DEGs
        sum(df_$logFC < 0) -> CT8_DEGs
        gene_list_for_GO[[region_]][[celltype_]][["CT20_enriched"]] <<- df_ %>% dplyr::filter(logFC > 0) %>% .$gene
        gene_list_for_GO[[region_]][[celltype_]][["CT8_enriched"]] <<- df_ %>% dplyr::filter(logFC < 0) %>% .$gene
        rbind(
          data.frame(region = region_, celltype = celltype_, n_DEGs = CT20_DEGs, direction = "CT-20", order = order_),
          data.frame(region = region_, celltype = celltype_, n_DEGs = -CT8_DEGs, direction = "CT-8", order = order_)
        )
      }) %>% do.call(rbind, .) -> df_
    df_ %>% dplyr::filter(n_DEGs != 0) -> df_
    df_ %>% dplyr::arrange(desc(order)) -> df_
    df_[, c("celltype", "order")] %>% distinct() %>% .$celltype %>% rev() -> levels_
    df_$celltype = factor(df_$celltype, levels = levels_)
    df_ %>% dplyr::filter(celltype %in% tail(levels_, 3)) -> df_
    df_ %>% 
      ggplot(aes(y = celltype, x = n_DEGs, fill = direction)) + 
      geom_col() + 
      scale_fill_manual(values = c("CT-8" = "#f8766d", "CT-20" = "#00bfc4")) +
      geom_text(aes(
        x = n_DEGs / 2,
        label = abs(n_DEGs)
      ),
      color = "black", size = 4
      ) + 
      theme_classic() + 
      xlab("Number of DEGs") +
      ggtitle(region_) -> p_bar
  }) %>% patchwork::wrap_plots(ncol = 5, guides = "collect") -> p_1

ggsave("~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_SON_PVN_DEG_barplots_1.pdf", plot = p_1, width = 12.52, height = 6.88)
ggsave("~/Downloads/TVGs.pdf", plot = p_1, width = 18.83, height = 3.06)

DEG_df %>% dplyr::filter(region == "SCN", compartment == "cell") %>% 
  {
    df_ = .
    celltypes_ = c("GABAergic", "Astrocytes", "Vascular")
    celltypes_ %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(celltype_){
        df_ %>% dplyr::filter(cluster_id == celltype_) %>% 
          dplyr::select(gene, logFC, p_val, p_adj.loc) -> df_plot
        df_plot %>% dplyr::mutate(sig = case_when(
          p_adj.loc < 0.05 & logFC > 0 ~ "CT-20 enriched",
          p_adj.loc < 0.05 & logFC < 0 ~ "CT-8 enriched",
          TRUE ~ "Not significant"
        )) -> df_plot
        
        df_plot %>% 
          ggplot(aes(x = logFC, y = -log10(p_adj.loc), color = sig)) + 
          geom_point() + 
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") + 
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + 
          scale_color_manual(values = c("CT-20 enriched" = "#00bfc4", "CT-8 enriched" = "#f8766d", "Not significant" = "grey")) +
          theme_classic() + 
          ylab("−log10(p.adj)") + 
          ggtitle(celltype_)
      }) -> p_list
    patchwork::wrap_plots(p_list, ncol = 3, guides = 'collect')
  } #Figure 2B

#6. GO Analysis ----
library(clusterProfiler)
library(org.Mm.eg.db)

read.table("~/Dropbox/Spatial_transcriptomics_new/rda/GRCm39_mart_export.txt", header = T, sep = "\t", stringsAsFactors = F) -> gene.ref
gene.ref %>% dplyr::filter(Gene.name != "") %>% .[,c(1,3)] %>% 
  distinct() -> gene.ref
readRDS("~/Dropbox/singulomics/github_rda/gene.ref.rds") -> gene.ref.old


gene_list_for_GO %>% 
  #  .["SCN"] %>% 
  purrr::map(function(list_){
    list_ %>% 
      #      .["GABAergic"] %>% 
      purrr::map(function(list__){
        list__ %>% 
          #          .["CT20_enriched"] %>% 
          purrr::map(function(gene_){
            gene_ = gene_
            if (length(gene_) > 0){
              gene.ref[gene.ref$Gene.name %in% gene_,] %>% 
                .$Gene.stable.ID -> ENSEMBL_id
              
              gene_[!(gene_ %in% gene.ref$Gene.name)] -> gene_1
              if (length(gene_1) > 0){
                gene.ref.old[gene.ref.old$gene_name %in% gene_1,] %>% 
                  .$gene_id -> ENSEMBL_id_1              
                ENSEMBL_id = c(ENSEMBL_id, ENSEMBL_id_1) %>% unique()
              }              
            }else{
              ENSEMBL_id = character(0)
            }
            sprintf("gene name: %s, ensembl id: %s", length(gene_), length(ENSEMBL_id)) %>% print()
            ego_CC = enrichGO(gene = ENSEMBL_id, 
                              OrgDb = org.Mm.eg.db, 
                              ont = "CC", 
                              pAdjustMethod = "BH", 
                              minGSSize = 1, 
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.05, 
                              readable = TRUE, 
                              keyType = "ENSEMBL")
            
            ego_BP = enrichGO(gene = ENSEMBL_id, 
                              OrgDb = org.Mm.eg.db, 
                              ont = "BP", 
                              pAdjustMethod = "BH", 
                              minGSSize = 1, 
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.05, 
                              readable = TRUE, 
                              keyType = "ENSEMBL")
            
            ego_MF = enrichGO(gene = ENSEMBL_id, 
                              OrgDb = org.Mm.eg.db, 
                              ont = "MF", 
                              pAdjustMethod = "BH", 
                              minGSSize = 1, 
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.05, 
                              readable = TRUE, 
                              keyType = "ENSEMBL")
            
            AnnotationDbi::select(org.Mm.eg.db, keys = ENSEMBL_id, 
                                  columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                                  keytype = "ENSEMBL") %>% drop_na() -> entrez_id
            kk <- enrichKEGG(gene = entrez_id$ENTREZID, organism ="mmu", pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05, minGSSize = 1,
                             #readable = TRUE, 
                             use_internal_data =FALSE)
            
            list(ego_CC = ego_CC, ego_BP = ego_BP, ego_MF = ego_MF, kk = kk)
          })
      })
  }) -> res_go_analysis

list(
  SCN = c(
    "GABAergic",
    "Astrocytes"
  ), 
  SCN_core = c(
    "GABAergic"
  ), 
  SCN_shell = c(
    "GABAergic",
    "Astrocytes"
  ), 
  PVN = c(
    "Astrocytes",
    "GABAergic"
  ),
  SON = c(
    "Astrocytes",
    "Glutamatergic"
  )
) %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.), .f=function(celltype_, region_){
    res_go_analysis[[region_]] -> list_
    celltype_ %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(celltype__){
        list_[[celltype__]] -> list__
        list__ %>% 
          purrr::map2(.x=.,.y=names(.),.f=function(list___, ZT_){
            list___ %>% 
              .[c("ego_BP", "kk")] %>% 
              purrr::map2(.x=.,.y=names(.),.f=function(go_, go__){
                go_@result -> df_
                df_ %>% dplyr::filter(p.adjust < 0.05) -> df_
                df_ = df_[1:10, ] %>% drop_na()
                df_$Description = factor(df_$Description, levels = rev(df_$Description))
                df_ = df_[,c("Description", "p.adjust", "geneID")]
                if (go__ == "kk"){
                  df_$geneID %>% str_split("/") -> entrez_id_list
                  entrez_id_list %>% 
                    purrr::map(function(entrez_id_){
                      AnnotationDbi::select(org.Mm.eg.db, keys = entrez_id_, 
                                            columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                                            keytype = "ENTREZID") %>% drop_na() -> gene_info
                      gene_info$SYMBOL %>% paste(collapse = "/")
                    }) %>% unlist() -> df_$geneID
                }
                df_$Description = sprintf("%s (%s)", go__, df_$Description)
                df_
              }) %>% do.call(rbind, .) %>% 
              dplyr::arrange(p.adjust) -> df_
            df_$Description = factor(df_$Description, levels = rev(df_$Description))
            df_ %>% 
              ggplot(aes(y = Description, x = -log10(p.adjust))) + 
              geom_bar(stat = "identity") + 
              ggtitle(sprintf("%s - %s - %s", region_, celltype__, ZT_)) + 
              theme_classic() -> p
          })  -> p_list_1
        patchwork::wrap_plots(p_list_1, ncol = 2)
      })
  }) -> p_list_2

pdf("~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_PVN_SON_GO_analysis.pdf", width = 15.5, height = 6.2)
p_list_2
dev.off()

p_list_2$SCN$GABAergic[[1]]$data %>% dplyr::filter(grepl("Circadian", Description)) %>% .$geneID %>% str_split("/") %>% unlist() %>% unique() -> CT8_genes
p_list_2$SCN$GABAergic[[2]]$data %>% dplyr::filter(grepl("Circadian|rhythmic", Description, ignore.case = T)) %>% .$geneID %>% str_split("/") %>% unlist() %>% unique() -> CT20_genes

DEG_df %>% dplyr::filter(region == "SCN", compartment == "cell", cluster_id == "GABAergic") %>% 
  dplyr::mutate(group = case_when(
    (abs(logFC) > log2(1.1))&(p_adj.loc < 0.05)&(logFC > 0) ~ "CT20_enriched",
    (abs(logFC) > log2(1.1))&(p_adj.loc < 0.05)&(logFC < 0) ~ "CT8_enriched",
    TRUE ~ "Not_significant"
  )) %>% 
  #  .$group %>% table()
  ggplot(aes(x = logFC, y = -log10(p_adj.loc), colour = group)) + 
  geom_point() + 
  scale_color_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d", "Not_significant" = "grey")) + 
  #add labels for circadian rhythm related genes
  ggrepel::geom_text_repel(data = . %>% dplyr::filter(gene %in% c(CT8_genes, CT20_genes, "Arntl")), aes(label = gene), size = 3) +
  theme_classic()

DEG_df %>% dplyr::filter(region == "SCN", compartment == "cell", cluster_id == "GABAergic") %>% 
  dplyr::filter(gene %in% c(CT8_genes, CT20_genes, "Arntl")) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("Circadian rhythm - SCN GABAergic") -> p_circadian_SCN_GABAergic

p_list_2$SCN$GABAergic[[2]]$data %>% 
  dplyr::filter(grepl("Neurotrophin", Description)) %>% 
  .$geneID %>% 
  str_split("/") %>% unlist() %>% unique() -> neurotrophin_genes
DEG_df %>% dplyr::filter(region == "SCN", compartment == "cell", cluster_id == "GABAergic")

DEG_df %>% dplyr::filter(region == "SCN", compartment == "cell", cluster_id == "GABAergic") %>% 
  dplyr::filter(gene %in% neurotrophin_genes) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("Neurotrophin - SCN GABAergic") -> p_neurotrophin_SCN_GABAergic

p_list_2$SCN$GABAergic[[1]]$data %>% 
  dplyr::filter(grepl("axon|neurogenesis|forebrain|nervous", Description, ignore.case = T)) %>% 
  .$geneID %>% 
  str_split("/") %>% unlist() %>% unique() -> neuro_connectivity_genes

DEG_df %>% dplyr::filter(region == "SCN", compartment == "cell", cluster_id == "GABAergic") %>% 
  dplyr::filter(gene %in% neuro_connectivity_genes) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("Neuro plasticity - SCN GABAergic") -> p_neuro_connectivity_SCN_GABAergic

patchwork::wrap_plots(p_circadian_SCN_GABAergic, p_neurotrophin_SCN_GABAergic, p_neuro_connectivity_SCN_GABAergic, ncol = 3, guides = "collect") -> p
ggsave("~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_GABAergic_GO.pdf", p, width = 15.4, height = 7.01)

p_list_2$SCN$Astrocytes[[2]]$data %>% 
  dplyr::filter(grepl("circadian|rhythm", Description, ignore.case = T)) %>% 
  .$geneID %>% 
  str_split("/") %>% unlist() %>% unique() -> circadian_genes

DEG_df %>% dplyr::filter(region == "SCN", compartment == "cell", cluster_id == "Astrocytes") %>% 
  dplyr::filter(gene %in% circadian_genes) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("Circadian rhythm - SCN Astrocytes") -> p_circadian_SCN_Astrocytes

ggsave(p_circadian_SCN_Astrocytes, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_Astrocytes_GO.pdf", width = 8.33, height = 7.01)

p_list_2$SCN_core$GABAergic[[2]]$data %>% 
  dplyr::filter(grepl("glutamate|plasticity|memory|neuroactive|MAPK", Description, ignore.case = T)) %>% 
  .$geneID %>% 
  str_split("/") %>% unlist() %>% unique() -> neuronal_signaling_genes

DEG_df %>% dplyr::filter(region == "SCN_core", compartment == "cell", cluster_id == "GABAergic") %>% 
  dplyr::filter(gene %in% neuronal_signaling_genes) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("neuronal signaling and activity-dependent plasticity - SCN core GABAergic") -> p_neuro_signaling_SCN_core_GABAergic
ggsave(p_neuro_signaling_SCN_core_GABAergic, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_core_GABA_GO.pdf", width = 8.33, height = 7.01)

p_list_2$SCN_shell$GABAergic[[2]]$data %>% 
  dplyr::filter(grepl("filopodium|projection|synapse", Description, ignore.case = T)) %>% 
  .$geneID %>% 
  str_split("/") %>% unlist() %>% unique() -> Neurite_and_synaptic_remodeling_genes

DEG_df %>% dplyr::filter(region == "SCN_shell", compartment == "cell", cluster_id == "GABAergic") %>% 
  dplyr::filter(gene %in% Neurite_and_synaptic_remodeling_genes) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("Neurite and synaptic remodeling - SCN shell GABAergic") -> p_Neurite_and_synaptic_remodeling_SCN_shell_GABAergic
ggsave(p_Neurite_and_synaptic_remodeling_SCN_shell_GABAergic, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_shell_GABA_GO.pdf", width = 8.33, height = 7.01)

p_list_2$PVN$Astrocytes[[1]]$data %>% 
  dplyr::filter(grepl("starvation|nutrient|lipid", Description, ignore.case = T)) %>% 
  .$geneID %>% 
  str_split("/") %>% unlist() %>% unique() -> nutrient_genes

DEG_df %>% dplyr::filter(region == "PVN", compartment == "cell", cluster_id == "Astrocytes") %>% 
  dplyr::filter(gene %in% nutrient_genes) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("Nutrient - PVN Astrocytes") -> p_nutrient_PVN_Astrocytes

p_list_2$PVN$Astrocytes[[2]]$data %>% 
  dplyr::filter(grepl("circadian|rhythm", Description, ignore.case = T)) %>% 
  .$geneID %>% 
  str_split("/") %>% unlist() %>% unique() -> circadian_genes

DEG_df %>% dplyr::filter(region == "PVN", compartment == "cell", cluster_id == "Astrocytes") %>% 
  dplyr::filter(gene %in% circadian_genes) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("circadian - PVN Astrocytes") -> p_circadian_PVN_Astrocytes

p_list_2$PVN$Astrocytes[[2]]$data %>% 
  dplyr::filter(grepl("synapse", Description, ignore.case = T)) %>% 
  .$geneID %>% 
  str_split("/") %>% unlist() %>% unique() -> synapse_genes

DEG_df %>% dplyr::filter(region == "PVN", compartment == "cell", cluster_id == "Astrocytes") %>% 
  dplyr::filter(gene %in% synapse_genes) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("synpase - PVN Astrocytes") -> p_synpase_PVN_Astrocytes

patchwork::wrap_plots(p_nutrient_PVN_Astrocytes, p_circadian_PVN_Astrocytes, p_synpase_PVN_Astrocytes, ncol = 3, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/PVN_Astrocytes_GO.pdf", width = 14.44, height = 4.94)

p_list_2$SON$Astrocytes[[1]]$data %>% 
  dplyr::filter(grepl("glial|glio|Wnt", Description, ignore.case = T)) %>% 
  .$geneID %>% 
  str_split("/") %>% unlist() %>% unique() -> gilo_wnt_genes

DEG_df %>% dplyr::filter(region == "SON", compartment == "cell", cluster_id == "Astrocytes") %>% 
  dplyr::filter(gene %in% gilo_wnt_genes) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("Glio/Wnt - SON Astrocytes") -> p_Glio_Wnt_SON_Astrocytes

p_list_2$SON$Astrocytes[[2]]$data %>% 
  dplyr::filter(grepl("circadian|peptide", Description, ignore.case = T)) %>% 
  .$geneID %>% 
  str_split("/") %>% unlist() %>% unique() -> circadian_piptide_genes


DEG_df %>% dplyr::filter(region == "SON", compartment == "cell", cluster_id == "Astrocytes") %>% 
  dplyr::filter(gene %in% circadian_piptide_genes) %>% dplyr::arrange(desc(logFC)) %>% 
  dplyr::mutate(group = case_when(
    logFC > 0 ~ "CT20_enriched",
    logFC < 0 ~ "CT8_enriched",
  )) %>% 
  dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
  ggplot(aes(x = logFC, y = gene, fill = group)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("CT20_enriched" = "#00bfc4", "CT8_enriched" = "#f8766d")) + 
  theme_classic() + 
  ggtitle("circadian/peptide - SON Astrocytes") -> p_circadian_peptide_SON_Astrocytes

patchwork::wrap_plots(p_Glio_Wnt_SON_Astrocytes, p_circadian_peptide_SON_Astrocytes, ncol = 2, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SON_Astrocytes_GO.pdf", width = 14.18, height = 5.53)

#8. Compartment specific TVGs in SCN GABA ----
load("~/Downloads/DEG_df.Rdata")
load("~/Downloads/muscat_res_list.rda")
load("~/Downloads/DEG_dfs.rda")
load("~/Dropbox/Spatial_transcriptomics_new/external_data/Wen_mouse_SCN_circadian/results.rda")
load("~/Downloads/DEG_muscat_res_list.rda")

DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic") %>% 
  dplyr::mutate(Sig = case_when(
    (p_adj.loc < 0.05)&(logFC > log2(1.1)) ~ "CT20_enriched",
    (p_adj.loc < 0.05)&(logFC < -log2(1.1)) ~ "CT8_enriched",
    TRUE ~ "Not_significant"
  )) %>% 
  {
    df_ = .
    df_$compartment %>% 
      unique() %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(compartment_){
        df_ %>% dplyr::filter(compartment == compartment_) %>% 
          dplyr::filter(Sig != "Not_significant") %>% 
          .$gene
      }) %>% 
      ggvenn::list_to_data_frame() -> df_1
    
    df_1 %>% dplyr::filter(cell == TRUE, cyto == TRUE, nuc == TRUE) %>% .$key -> all_sig
    df_1 %>% dplyr::filter(cell == TRUE, cyto == FALSE, nuc == FALSE) %>% .$key -> cell_sig
    df_1 %>% dplyr::filter(cell == FALSE, cyto == TRUE, nuc == FALSE) %>% .$key -> cyto_sig
    df_1 %>% dplyr::filter(cell == FALSE, cyto == FALSE, nuc == TRUE) %>% .$key -> nuc_sig
    
    list(
      all_sig = all_sig,
      cell_sig = cell_sig,
      cyto_sig = cyto_sig,
      nuc_sig = nuc_sig
    )
  } -> compartment_specific_TVGs

DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic") %>% 
  dplyr::mutate(Sig = case_when(
    (p_adj.loc < 0.05)&(logFC > log2(1.1)) ~ "CT20_enriched",
    (p_adj.loc < 0.05)&(logFC < -log2(1.1)) ~ "CT8_enriched",
    TRUE ~ "Not_significant"
  )) %>% 
  {
    df_ = .
    df_$compartment %>% 
      unique() %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(compartment_){
        df_ %>% dplyr::filter(compartment == compartment_) %>% 
          dplyr::filter(Sig != "Not_significant") %>% 
          .$gene
      }) %>% 
      ggvenn::list_to_data_frame() -> df_ggvenn
    
    df_ggvenn ->> df_ggvenn
    df_ggvenn
    
  } %>% 
  dplyr::mutate(across(!matches("key"), function(x){x %>% as.integer()})) %>%
  as.data.frame() %>% 
  UpSetR::upset(sets = rev(c("cell", "nuc", "cyto")), keep.order = T, order.by = "freq"
  ) -> p_upset

muscat_res_list %>% 
#  .[1] %>% 
  purrr::map(function(list_){
    list_$SCN$GABAergic$ds_pb$data$GABAergic$counts %>% 
      as.data.frame() %>% 
      dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6}))
  }) -> muscat_res_norm_count_matrix

c("Hdac4", "Ogt") %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  purrr::map(function(gene_){
    muscat_res_norm_count_matrix %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(df_, compartment_){
        df_[gene_, ] %>% 
          pivot_longer(cols = everything(), names_to = "sample", values_to = "TPM") %>% 
          dplyr::mutate(compartment = compartment_, timepoint = gsub("(CT-\\d+?)-.+", "\\1", sample))
      }) %>% do.call(rbind, .) -> df_plot
    df_plot %>% 
      ggplot(aes(x = compartment, y = TPM, color = timepoint, shape = timepoint)) + 
      geom_point(alpha = 0.7, size = 2) + 
      theme_classic() -> p
    DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic", gene == gene_, compartment == "cell") %>% .$p_adj.loc %>% sprintf("%.2e", .) -> cell_p.adj
    DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic", gene == gene_, compartment == "nuc") %>% .$p_adj.loc %>% sprintf("%.2e", .) -> nuc_p.adj
    DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic", gene == gene_, compartment == "cyto") %>% .$p_adj.loc %>% sprintf("%.2e", .) -> cyto_p.adj
    p + 
      ggtitle(sprintf("%s\np.adj = %s (cell), %s (cyto), %s (nuc)", gene_, cell_p.adj, cyto_p.adj, nuc_p.adj))
  }) %>% patchwork::wrap_plots(., ncol = 2, guides = "collect")

DEG_df %>% dplyr::filter(region == "SCN", gene %in% c("Hdac4", "Ogt")) %>% 
  dplyr::filter(p_adj.loc < 0.05)

DEG_df_SCN %>% dplyr::filter(gene %in% c("Hdac4", "Ogt"))

df_ggvenn %>% 
  {
    df_ = .
    list(
    nuc_only = df_ %>% dplyr::filter(cell == FALSE, cyto == FALSE, nuc == TRUE) %>% .$key,
    cyto_only = df_ %>% dplyr::filter(cell == FALSE, cyto == TRUE, nuc == FALSE) %>% .$key
    )
  } -> gene_list

gene_list$nuc_only %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(gene_){
    muscat_res_norm_count_matrix %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(df_, compartment_){
        df_[gene_, ] %>% 
          pivot_longer(cols = everything(), names_to = "sample", values_to = "TPM") %>% 
          dplyr::mutate(compartment = compartment_, timepoint = gsub("(CT-\\d+?)-.+", "\\1", sample))
      }) %>% do.call(rbind, .) -> df_plot
    df_plot %>% 
      ggplot(aes(x = compartment, y = TPM, color = timepoint, shape = timepoint)) + 
      geom_point(alpha = 0.7, size = 2) + 
      theme_classic() -> p
    DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic", gene == gene_, compartment == "cell") %>% .$p_adj.loc %>% sprintf("%.2e", .) -> cell_p.adj
    DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic", gene == gene_, compartment == "nuc") %>% .$p_adj.loc %>% sprintf("%.2e", .) -> nuc_p.adj
    DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic", gene == gene_, compartment == "cyto") %>% .$p_adj.loc %>% sprintf("%.2e", .) -> cyto_p.adj
    p + 
      ggtitle(sprintf("%s\np.adj = %s (cell), %s (cyto), %s (nuc)", gene_, cell_p.adj, cyto_p.adj, nuc_p.adj))
  }) -> p_list_nuc_only

gene_list$cyto_only %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(gene_){
    muscat_res_norm_count_matrix %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(df_, compartment_){
        df_[gene_, ] %>% 
          pivot_longer(cols = everything(), names_to = "sample", values_to = "TPM") %>% 
          dplyr::mutate(compartment = compartment_, timepoint = gsub("(CT-\\d+?)-.+", "\\1", sample))
      }) %>% do.call(rbind, .) -> df_plot
    df_plot %>% 
      ggplot(aes(x = compartment, y = TPM, color = timepoint, shape = timepoint)) + 
      geom_point(alpha = 0.7, size = 2) + 
      theme_classic() -> p
    DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic", gene == gene_, compartment == "cell") %>% .$p_adj.loc %>% sprintf("%.2e", .) -> cell_p.adj
    DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic", gene == gene_, compartment == "nuc") %>% .$p_adj.loc %>% sprintf("%.2e", .) -> nuc_p.adj
    DEG_df %>% dplyr::filter(region == "SCN", cluster_id == "GABAergic", gene == gene_, compartment == "cyto") %>% .$p_adj.loc %>% sprintf("%.2e", .) -> cyto_p.adj
    p + 
      ggtitle(sprintf("%s\np.adj = %s (cell), %s (cyto), %s (nuc)", gene_, cell_p.adj, cyto_p.adj, nuc_p.adj))
  }) -> p_list_cyto_only

gene_list %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(genes_, group_){
    genes_ %>% 
      "names<-"(.,.) %>%  
#      .[1] %>% 
      purrr::map(function(gene_){
        DEG_df %>% dplyr::filter(region == "SCN", gene == gene_) -> df_
        df_$cluster_id = factor(df_$cluster_id, levels = c("GABAergic", "Astrocytes", "Oligodendrocytes", "Vascular", "Immune"))
        df_$compartment = factor(df_$compartment, levels = rev(c("cell", "nuc", "cyto")))
        df_ %>% dplyr::mutate(log_p_adj = -log10(p_adj.loc)) -> df_
        df_ %>% dplyr::mutate(log_p_adj = case_when(
          logFC > 0 ~ log_p_adj,
          TRUE ~ -log_p_adj
        )) -> df_
        
        df_ %>% 
        dplyr::mutate(Sig_1 = case_when(
          p_adj.loc < 0.001 ~ "***",
          p_adj.loc < 0.01 ~ "**",
          p_adj.loc < 0.05 ~ "*",
          TRUE ~ ""
        )) -> df_
        
        df_ %>% 
          ggplot(aes(x = cluster_id, y = compartment, fill = log_p_adj)) + 
          geom_tile() + 
          scale_fill_gradient2(low = "red", high = "blue", mid = "grey", midpoint = 0) + 
          geom_text(aes(label = Sig_1), size = 8, color = "black") + 
          theme_classic() + 
          ylab("Compartment") + xlab("Celltype") + 
          ggtitle(gene_)
      })
  }) -> p_list_heat_map

gene_list %>% 
#  .[1] %>% 
  purrr::map(function(genes_){
    genes_ %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      purrr::map(function(gene_){
        GSE117295_bulk_norm_counts[gene_, ] %>% t() %>% 
          as.data.frame() %>% 
          "colnames<-"(., "norm_expr") %>% 
          rownames_to_column("REP") %>% 
          dplyr::mutate(ZT = gsub("ZT(\\d+).+", "\\1", REP) %>% as.numeric()) %>% 
          ggplot(aes(x = ZT, y = norm_expr)) + 
          geom_line() + 
          theme_classic() + 
          ylab("Normalized expression") + 
          xlab("CT") -> p
        JTK.i %>% dplyr::filter(CycID == gene_) %>% .$JTK_pvalue -> JTK_pvalue
        p + ggtitle(sprintf("SCN scRNA-seq\n%s: JTK_pval=%s", gene_, JTK_pvalue))
      })
  }) -> p_list_scRNA_seq

muscat_SCN_res$ds_pb$data$All_cells$counts %>% 
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6})) -> SCN_count_mat_norm
  
gene_list %>% 
#  .[1] %>%
  purrr::map(function(genes_){
    genes_ %>% 
      "names<-"(.,.) %>% 
#      .[1] %>%
      purrr::map(function(gene_){
        SCN_count_mat_norm[gene_, ] %>% 
          pivot_longer(cols = everything(), names_to = "Sample", values_to = "TPM") %>% 
          dplyr::mutate(Timepoint = gsub("(CT-\\d+?)-.+", "\\1", Sample)) -> df_
        df_$Timepoint = factor(df_$Timepoint, levels = c("CT-8", "CT-20"))
        df_ %>% 
          ggplot(aes(x = Timepoint, y = TPM, color = Timepoint, shape = Timepoint)) + 
          geom_point() + 
          theme_classic() -> p
        DEG_df_SCN %>% dplyr::filter(gene == gene_) %>% .$p_adj.loc %>% sprintf("%.2e",.) -> p_adj
        p + ggtitle(sprintf("Our study (SCN)\n%s padj=%s", gene_, p_adj))
      })
  }) -> p_list_TVGs_SCN

#nuc only candidate genes #Cyfip1, Gsk3b, Mfn2, Nr2f1, Ntsr2, P3h3, Pink1, Ppt1, Prokr2, Tppp 

DEG_df %>% dplyr::filter(region == "SCN", gene %in% c("Gsk3b")) %>% 
  dplyr::filter(p_adj.loc < 0.05)
DEG_df_SCN %>% dplyr::filter(gene == "Gsk3b")
p_list_heat_map$nuc_only$Gsk3b
patchwork::wrap_plots(p_list_scRNA_seq$nuc_only$Gsk3b, p_list_TVGs_SCN$nuc_only$Gsk3b, ncol = 2, widths = c(3,1))

#cyto only candidate genes #Fzd8, Mef2c, Ogt, Wif1
DEG_df %>% dplyr::filter(region == "SCN", gene %in% c("Ogt")) %>% 
  dplyr::filter(p_adj.loc < 0.05)
DEG_df_SCN %>% dplyr::filter(gene == "Ogt")
p_list_cyto_only$Ogt
p_list_heat_map$cyto_only$Ogt
patchwork::wrap_plots(p_list_scRNA_seq$cyto_only$Ogt, p_list_TVGs_SCN$cyto_only$Ogt, ncol = 2, widths = c(3,1))

DEG_df %>% dplyr::filter(region == "SCN", gene %in% c("Wif1", "Fzd8")) %>% 
  dplyr::filter(p_adj.loc < 0.05)
DEG_df_SCN %>% dplyr::filter(gene %in% c("Wif1", "Fzd8"))
p_list_cyto_only[c("Wif1", "Fzd8")] %>% patchwork::wrap_plots(ncol = 2, guides = "collect")




pdf("~/Downloads/nuc_only_TVG_heatmap.pdf", width = 5.74, height = 2.03)
p_list_heat_map$nuc_only
dev.off()

pdf("~/Downloads/cyto_only_TVG_heatmap.pdf", width = 5.74, height = 2.03)
p_list_heat_map$cyto_only
dev.off()

pdf("~/Downloads/p_list_nuc_only.pdf", width = 6.49, height = 4.56)
p_list_nuc_only
dev.off()

pdf("~/Downloads/p_list_cyto_only.pdf", width = 6.49, height = 4.56)
p_list_cyto_only
dev.off()

# 9. Compare p-value of core clock genes in GABA, SCN and SCN GABA----
load("~/Downloads/DEG_df.Rdata") # region x compartment x celltype
load("~/Downloads/DEG_dfs.rda") # celltypes (all cells)
load("~/Downloads/DEG_muscat_res_list.rda") # regions (all cells)
load("~/Downloads/muscat_res_list.rda")

##SCN----
c("Arntl", "Npas2", "Dbp", "Nr1d1", "Nr1d2", "Rora", "Cry1", "Per1", "Per2", "Per3", "Hlf", "Tef", "Rgs16") %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(gene_){
    DEG_df_SCN %>% dplyr::filter(gene == gene_) %>% .$p_adj.loc -> SCN_padj
    DEG_df %>% dplyr::filter(region == "SCN", compartment == "cell", cluster_id == "GABAergic", gene == gene_)  %>% .$p_adj.loc -> SCN_GABA_padj
    DEG_df_celltype_1 %>% dplyr::filter(cluster_id == "GABAergic", gene == gene_) %>% .$p_adj.loc -> GABA_padj
    data.frame(group = factor(c("GABAergic", "SCN", "SCN_GABA"), levels = c("GABAergic", "SCN", "SCN_GABA"))) %>% 
      dplyr::mutate(p_adj = c(GABA_padj, SCN_padj, SCN_GABA_padj)) -> df_plot
    df_plot %>% 
      ggplot(aes(x = group, y = -log10(p_adj))) + 
      geom_col(alpha = 0.7) + 
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + 
      theme_classic() + 
      xlab(NULL) + 
      ylab("-log10(p.adj)") + 
      ggtitle(gene_) + 
      theme(plot.title = element_text(face = "italic"))
  }) -> p_list_1

patchwork::wrap_plots(p_list_1, ncol = 4)
p_list_1 = p_list_1[c(2,4,6,7,9)]

names(p_list_1) %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(gene_){
    muscat_SCN_res$ds_pb$data$All_cells$counts %>% 
      as.data.frame() %>% 
      dplyr::mutate(across(everything(), function(x){(x/sum(x)*10^6)})) -> SCN_matrix
    muscat_res_celltype_list_1$GABAergic$ds_pb$data$GABAergic$counts %>% 
      as.data.frame() %>% 
      dplyr::mutate(across(everything(), function(x){(x/sum(x)*10^6)})) -> GABA_matrix
    muscat_res_list$cell$SCN$GABAergic$ds_pb$data$GABAergic$counts %>% 
      as.data.frame() %>% 
      dplyr::mutate(across(everything(), function(x){(x/sum(x)*10^6)})) -> SCN_GABA_matrix 
    
    rbind(
      SCN_matrix[gene_, ] %>% 
        pivot_longer(cols = everything(), names_to = "sample", values_to = "TPM") %>% 
        dplyr::mutate(group = "SCN"),
      GABA_matrix[gene_, ] %>% 
        pivot_longer(cols = everything(), names_to = "sample", values_to = "TPM") %>% 
        dplyr::mutate(group = "GABAergic"), 
      SCN_GABA_matrix[gene_, ] %>% 
        pivot_longer(cols = everything(), names_to = "sample", values_to = "TPM") %>% 
        dplyr::mutate(group = "SCN_GABA")
    ) %>% dplyr::mutate(group = factor(group, levels = c("GABAergic", "SCN", "SCN_GABA"))) %>% 
      dplyr::mutate(Timepoint = gsub("(CT-\\d+?)-.+", "\\1", sample) %>% factor(., levels = c("CT-8", "CT-20"))) -> df_plot
    df_plot %>% 
      ggplot(aes(x = group, y = TPM, color = Timepoint, shape = Timepoint)) + 
      geom_point(size = 3, alpha = 0.7) + 
      xlab(NULL) + 
      theme_classic() + 
      ggtitle(gene_) + 
      theme(plot.title = element_text(face = "italic"))
  }) -> p_list_2

patchwork::wrap_plots(p_list_1, ncol = 5) -> p1
patchwork::wrap_plots(p_list_2, ncol = 5) -> p2
patchwork::wrap_plots(p1,p2, ncol = 1, guides = "collect") -> p
ggsave(plot = p, filename = "~/Downloads/SCN_DEG_plot.pdf", height = 3.82, width = 15.72) #Figure 2D

muscat_res_list$cell$SCN$GABAergic$ds_pb$data$GABAergic$counts %>% 
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), function(x){(x/sum(x)*10^6)})) -> SCN_GABA_matrix 

c("Arntl", "Dbp", "Per1", "Per3") %>% 
#  .[1] %>% 
  purrr::map(function(gene_){
    SCN_GABA_matrix[gene_, ] %>% 
      pivot_longer(cols = everything(), names_to = "sample", values_to = "TPM") %>% 
      dplyr::mutate(group = "GABAergic") %>% 
      dplyr::mutate(Timepoint = gsub("(CT-\\d+?)-.+", "\\1", sample) %>% factor(., levels = c("CT-8", "CT-20"))) -> df_plot
    df_plot %>% 
      ggplot(aes(x = Timepoint, y = TPM, color = Timepoint, shape = Timepoint)) + 
      geom_point(size = 3, alpha = 0.7) + 
      xlab(NULL) + 
      theme_classic() + 
      ggtitle(gene_) + 
      theme(plot.title = element_text(face = "italic"))
  }) %>% patchwork::wrap_plots(., ncol = 4, guides = "collect") -> p

DEG_df %>% dplyr::filter(region == "SCN", compartment == "cell", cluster_id == "GABAergic", gene %in% c("Arntl", "Dbp", "Per1", "Per3"))

#10. SCN dot plot ----
rm(list=ls())
load("~/Downloads/DEG_df.Rdata")
load("~/Downloads/DEG_muscat_res_list.rda") # regions (all cells)

muscat_SCN_res$ds_pb$data$All_cells$counts %>% 
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), function(x){((x/sum(x))*10^6)}) ) -> SCN_matrix

SCN_matrix[c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e"), ] %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  as.data.frame() %>% rownames_to_column("sample") %>% 
  pivot_longer(cols = -sample, names_to = "gene", values_to = "expr") -> df_expr

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

xy_coordinates %>% 
  purrr::map(function(df_){
    df_
    df_ %>% dplyr::filter(cell %in% colnames(xenium_obj_filtered)) -> df_
  }) -> xy_coordinates

meta_data = xenium_obj_filtered@meta.data

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
xenium_obj_filtered = subset(xenium_obj_filtered, cells = scn_cells)
GetAssayData(object = xenium_obj_filtered, slot = "counts", assay = "Xenium") -> raw_data_SCN
class(raw_data_SCN)
raw_data_SCN[c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e"), ] %>% as.data.frame() -> raw_data_SCN

non_zero_prop_func = function(x){
  n_cells = length(x)
  n_non_zero = sum(x > 0)
  return((n_non_zero/n_cells)*100)
}
raw_data_SCN %>% rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "cell", values_to = "counts") %>% 
  dplyr::mutate(sample = gsub("(.+?)_.+", "\\1", cell)) %>% 
  dplyr::select(-cell) %>% 
  group_by(gene, sample) %>% 
  summarise(non_zero_prop = non_zero_prop_func(counts)) %>% ungroup() -> df_non_zero_prop

dim(df_expr)
dim(df_non_zero_prop)

left_join(x = df_expr, y = df_non_zero_prop, by = c("gene", "sample")) %>% 
  dplyr::mutate(sample = factor(sample, levels = rev(c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3")))) %>% 
  dplyr::mutate(gene = factor(gene, levels = c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e"))) -> df_plot

df_plot %>% 
  ggplot(aes(y = sample, x = gene)) +
  geom_point(aes(size = non_zero_prop, color = expr)) +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  scale_size_continuous(range = c(0, 15), limits = c(0,100)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Gene", color = "Expression (scaled)", size = "Non-zero proportion (%)") -> p

ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_clock_genes_expression_non_zero_prop.pdf", width = 11.93, height = 5.33) #Figure 2A

#11. PVN dot plot ----
rm(list=ls())
load("~/Downloads/DEG_df.Rdata")
load("~/Downloads/DEG_muscat_res_list.rda") # regions (all cells)

muscat_PVN_res$ds_pb$data$All_cells$counts %>% 
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), function(x){((x/sum(x))*10^6)}) ) -> PVN_matrix

PVN_matrix[c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e"), ] %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  as.data.frame() %>% rownames_to_column("sample") %>% 
  pivot_longer(cols = -sample, names_to = "gene", values_to = "expr") -> df_expr

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

xy_coordinates %>% 
  purrr::map(function(df_){
    df_
    df_ %>% dplyr::filter(cell %in% colnames(xenium_obj_filtered)) -> df_
  }) -> xy_coordinates

meta_data = xenium_obj_filtered@meta.data

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

pvn_cells = select_cells(region__ = "PVN_refined")
xenium_obj_filtered = subset(xenium_obj_filtered, cells = pvn_cells)
GetAssayData(object = xenium_obj_filtered, slot = "counts", assay = "Xenium") -> raw_data_PVN
class(raw_data_PVN)
raw_data_PVN[c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e"), ] %>% as.data.frame() -> raw_data_PVN

non_zero_prop_func = function(x){
  n_cells = length(x)
  n_non_zero = sum(x > 0)
  return((n_non_zero/n_cells)*100)
}
raw_data_PVN %>% rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "cell", values_to = "counts") %>% 
  dplyr::mutate(sample = gsub("(.+?)_.+", "\\1", cell)) %>% 
  dplyr::select(-cell) %>% 
  group_by(gene, sample) %>% 
  summarise(non_zero_prop = non_zero_prop_func(counts)) %>% ungroup() -> df_non_zero_prop

dim(df_expr)
dim(df_non_zero_prop)

left_join(x = df_expr, y = df_non_zero_prop, by = c("gene", "sample")) %>% 
  dplyr::mutate(sample = factor(sample, levels = rev(c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3")))) %>% 
  dplyr::mutate(gene = factor(gene, levels = c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e"))) -> df_plot

df_plot %>% 
  ggplot(aes(y = sample, x = gene)) +
  geom_point(aes(size = non_zero_prop, color = expr)) +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  scale_size_continuous(range = c(0, 15), limits = c(0,100)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Gene", color = "Expression (scaled)", size = "Non-zero proportion (%)") -> p

ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/PVN_clock_genes_expression_non_zero_prop.pdf", width = 11.93, height = 5.33)

#12. SON dot plot ----
rm(list=ls())
load("~/Downloads/DEG_df.Rdata")
load("~/Downloads/DEG_muscat_res_list.rda") # regions (all cells)

muscat_SON_res$ds_pb$data$All_cells$counts %>% 
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), function(x){((x/sum(x))*10^6)}) ) -> SON_matrix

SON_matrix[c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e"), ] %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  as.data.frame() %>% rownames_to_column("sample") %>% 
  pivot_longer(cols = -sample, names_to = "gene", values_to = "expr") -> df_expr

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

xy_coordinates %>% 
  purrr::map(function(df_){
    df_
    df_ %>% dplyr::filter(cell %in% colnames(xenium_obj_filtered)) -> df_
  }) -> xy_coordinates

meta_data = xenium_obj_filtered@meta.data

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

son_cells = select_cells(region__ = "SON")
xenium_obj_filtered = subset(xenium_obj_filtered, cells = son_cells)
GetAssayData(object = xenium_obj_filtered, slot = "counts", assay = "Xenium") -> raw_data_SON
class(raw_data_SON)
dim(raw_data_SON)
raw_data_SON[c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e"), ] %>% as.data.frame() -> raw_data_SON

non_zero_prop_func = function(x){
  n_cells = length(x)
  n_non_zero = sum(x > 0)
  return((n_non_zero/n_cells)*100)
}
raw_data_SON %>% rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "cell", values_to = "counts") %>% 
  dplyr::mutate(sample = gsub("(.+?)_.+", "\\1", cell)) %>% 
  dplyr::select(-cell) %>% 
  group_by(gene, sample) %>% 
  summarise(non_zero_prop = non_zero_prop_func(counts)) %>% ungroup() -> df_non_zero_prop

dim(df_expr)
dim(df_non_zero_prop)

left_join(x = df_expr, y = df_non_zero_prop, by = c("gene", "sample")) %>% 
  dplyr::mutate(sample = factor(sample, levels = rev(c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3")))) %>% 
  dplyr::mutate(gene = factor(gene, levels = c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e"))) -> df_plot

df_plot %>% 
  ggplot(aes(y = sample, x = gene)) +
  geom_point(aes(size = non_zero_prop, color = expr)) +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  scale_size_continuous(range = c(0, 15), limits = c(0,100)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Gene", color = "Expression (scaled)", size = "Non-zero proportion (%)") -> p

ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SON_clock_genes_expression_non_zero_prop.pdf", width = 11.93, height = 5.33)


# 13. Compare Xenium and scRNA-seq (SCN) ----
library(tidyverse)
load("~/Downloads/SCN_cells.rda")
#load("~/Downloads/DEG_muscat_res_list.rda")
load("~/Downloads/muscat_res_list.rda")

names(scn_cells) %in% c("GABAergic") -> se
scn_cells = scn_cells[se]
scn_cells = unname(scn_cells)
gsub("(CT-\\d+?)-.+", "\\1", scn_cells) %>% table() %>% as.numeric() -> n_cells

muscat_res_list$cell$SCN$GABAergic$dd_pb$data$GABAergic$counts %>% 
{
  df_ = .
  CT20_counts = df_[,1:3]
  CT20_counts %>% rowSums() -> CT20_counts
  CT8_counts = df_[,4:6]
  CT8_counts %>% rowSums() -> CT8_counts
  data.frame(CT20 = CT20_counts, CT8 = CT8_counts) %>% 
    apply(., 1, function(x){x/n_cells}) %>% t() %>% 
    as.data.frame()
} -> df_percent_xenium

#muscat_SCN_res$dd_pb$data$All_cells$counts %>% 
#  {
#    df_ = .
#    CT20_counts = df_[,1:3]
#    CT20_counts %>% rowSums() -> CT20_counts
#    CT8_counts = df_[,4:6]
#    CT8_counts %>% rowSums() -> CT8_counts
#    data.frame(CT20 = CT20_counts, CT8 = CT8_counts) %>% 
#      apply(., 1, function(x){x/n_cells}) %>% t() %>% 
#      as.data.frame()
#  } -> df_percent_xenium

muscat_res_list$cell$SCN$GABAergic$ds_pb$data$GABAergic$counts %>% 
  {
    df_ = .
    CT20_counts = df_[,1:3]
    CT20_counts %>% rowSums() -> CT20_counts
    CT8_counts = df_[,4:6]
    CT8_counts %>% rowSums() -> CT8_counts 
    data.frame(CT20 = CT20_counts, CT8 = CT8_counts) %>% 
      dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6}))
  } -> df_expr_xenium

#muscat_SCN_res$ds_pb$data$All_cells$counts %>% 
#  {
#    df_ = .
#    CT20_counts = df_[,1:3]
#    CT20_counts %>% rowSums() -> CT20_counts
#    CT8_counts = df_[,4:6]
#    CT8_counts %>% rowSums() -> CT8_counts 
#    data.frame(CT20 = CT20_counts, CT8 = CT8_counts) %>% 
#      dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6}))
#  } -> df_expr_xenium

#load("~/Downloads/scRNA_dat.rda")
list.files("~/Dropbox/Spatial_transcriptomics_new/external_data/Wen_mouse_SCN_circadian/GSE117295_RAW", full.names = T) %>% 
  .[grep("CT", .)] %>% 
  "names<-"(., gsub("/.+/.+(CT\\d+)\\.csv\\.gz", "\\1", .)) %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(file_, CT_){
    print(file_)
    read.csv(file = file_, row.names = 1, header = T) -> df_
  }) -> mat_list

barcode_list = list()
list.files("~/Dropbox/singulomics/github_rda/GSE117295_RAW_Mouse_SCN/cluster_vs_barcode", pattern = "\\.tsv", full.names = T) %>% .[2:3] %>% 
  purrr::map(function(files_){
    read.table(file = files_, header = T, sep = "\t", stringsAsFactors = F) -> df_
    colnames(df_) <- c("barcode", "celltypes")
    df_$celltypes %>% unique() %>% sort() -> celltypes_
    celltypes_ %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(celltype_){
        df_ %>% dplyr::filter(celltypes == celltype_) %>% .$barcode -> barcodes_
        barcode_list[[celltype_]] <<- barcodes_
      })
  })

barcode_list_1 = list()
barcode_list[c("Neurons1", "Neurons2", "Neurons3", "Neurons4", "Neurons5", "Neurons6")] %>% unlist() %>% unname() %>% unique() -> barcode_list_1[["Neurons"]]

mat_list %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_,CT_){
    colnames(df_) %in% barcode_list_1[["Neurons"]] -> se
    df_ = df_[,se]
    rowSums(df_) %>% as.data.frame() %>% "colnames<-"(.,CT_) %>% rownames_to_column("gene")
  }) %>% purrr::reduce(., inner_join, by = "gene") -> mat

mat %>% column_to_rownames("gene") -> mat

cbind(
  mat[,c("CT18", "CT22", "CT42", "CT46")] %>% rowSums() %>% as.data.frame() %>% "colnames<-"(., "night"),
  mat[,c("CT30", "CT34", "CT54", "CT58")] %>% rowSums() %>% as.data.frame() %>% "colnames<-"(., "day")
) -> mat

mat %>% 
  dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6})) %>% 
  rownames_to_column("gene") -> df_expr_scRNA

mat_list %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, CT_){
#    ncol(df_)
    colnames(df_) %in% barcode_list_1[["Neurons"]] -> se
    sum(se)
  }) %>% unlist() -> n_cells
c(
n_cells[c("CT18", "CT22", "CT42", "CT46")] %>% sum() %>% "names<-"(., "night"),
n_cells[c("CT30", "CT34", "CT54", "CT58")] %>% sum() %>% "names<-"(., "day")
) -> n_cells

mat_list %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, CT_){
    colnames(df_) %in% barcode_list_1[["Neurons"]] -> se
    df_ = df_[,se]
    df_ %>% apply(.,1,function(x){sum(x!=0)}) %>% as.data.frame() %>% "colnames<-"(.,CT_) %>% rownames_to_column("gene")
  }) %>% purrr::reduce(., inner_join, by = "gene") -> mat_1
mat_1 = mat_1 %>% column_to_rownames("gene")

mat_1 %>% 
  {
    df_ = .
    
    cbind(
    df_[,c("CT18", "CT22", "CT42", "CT46")] %>% rowSums() %>% as.data.frame() %>% "colnames<-"(., "night"),
    df_[,c("CT30", "CT34", "CT54", "CT58")] %>% rowSums() %>% as.data.frame() %>% "colnames<-"(., "day")
    )
  } -> mat_1

apply(mat_1, 1, function(x){x/n_cells}) %>% t() %>% as.data.frame() -> df_non_zero_prop_scRNA
df_expr_scRNA %>% column_to_rownames("gene") -> df_expr_scRNA

df_expr_scRNA %>% dplyr::filter((night == 0), (day == 0)) %>% dim()
df_non_zero_prop_scRNA %>% dplyr::filter((night == 0), (day == 0)) %>% dim()

df_expr_scRNA %>% dplyr::filter((night != 0)|(day != 0)) -> df_expr_scRNA
df_non_zero_prop_scRNA %>% dplyr::filter((night != 0)|(day != 0)) -> df_non_zero_prop_scRNA

df_expr_xenium %>% dplyr::filter((CT20 == 0), (CT8 == 0)) %>% dim()
df_percent_xenium %>% dplyr::filter((CT20 == 0), (CT8 == 0)) %>% dim()

common_genes = intersect(rownames(df_expr_xenium), rownames(df_expr_scRNA))

library(ggrepel)
c("night", "day") %>% 
  purrr::map(function(time_){
    df_expr_xenium[common_genes, ] %>% "colnames<-"(., c("night", "day")) -> df_xenium
    df_expr_scRNA[common_genes, ] -> df_scRNA    
    xenium = df_xenium[,time_] %>% as.numeric() #%>% log1p()
    scRNA = df_scRNA[,time_] %>% as.numeric() #%>% log1p()
    cor.test(xenium, scRNA, method = "spearman") %>% .$estimate %>% as.numeric() -> cor_
    data.frame(gene = common_genes, xenium = xenium, scRNA = scRNA) -> df_plot
    
    genes_show <- c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e")
    df_label <- subset(df_plot, gene %in% genes_show)
    
    df_plot %>% 
      ggplot(aes(x = xenium, y = scRNA)) +
#      geom_point(alpha = 0.5) + 
      geom_hex(bins = 35) + 
      scale_fill_gradient(
        low = "#DEEBF7",   # 低 count：淺藍
        high = "#08306B",  # 高 count：深藍
        name = "Count"
      ) + 
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "red",
        linewidth = 1
      ) +
      geom_point(
        data = df_label,
        shape = 21, fill = "yellow", color = "black", size = 2
      ) + 
      geom_label_repel(
        data = df_label,
        aes(label = gene),
        fill = "white",
        color = "black",
        size = 3,
        box.padding = 0.5,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 123
      ) + 
      theme_classic() + 
      ggtitle(sprintf("%s, r=%s", time_, round(cor_, 3))) + 
      ylab("scRNA TPM") +
      xlab("Xenium TPM")
#      ylab("scRNA log1p(TPM)") +
#      xlab("Xenium log1p(TPM)")
  }) %>% patchwork::wrap_plots(ncol = 2) -> p1

common_genes_1 = intersect(rownames(df_percent_xenium), rownames(df_non_zero_prop_scRNA))

c("night", "day") %>% 
  purrr::map(function(time_){
    df_percent_xenium[common_genes_1, ] %>% "colnames<-"(., c("night", "day")) -> df_xenium
    df_non_zero_prop_scRNA[common_genes_1, ] -> df_scRNA    
    xenium = df_xenium[,time_] %>% as.numeric()
    scRNA = df_scRNA[,time_] %>% as.numeric()
    cor.test(xenium, scRNA, method = "spearman") %>% .$estimate %>% as.numeric() -> cor_
    data.frame(gene = common_genes_1, xenium = xenium, scRNA = scRNA) -> df_plot
#    df_plot %>% dplyr::filter(xenium != 0) %>% dplyr::filter(scRNA != 0)  -> df_plot
    
    genes_show <- c("Arntl", "Clock", "Npas2", "Dbp", "Nr1d1", "Rora", "Rorb", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Csnk1e")
    df_label <- subset(df_plot, gene %in% genes_show)
    
    df_plot %>% 
      ggplot(aes(x = xenium, y = scRNA)) +
#      geom_point() + 
      geom_hex(aes(fill = after_stat(log10(count))), bins = 35) + 
      scale_fill_gradient(
        low = "#DEEBF7",   # 低 count：淺藍
        high = "#08306B",  # 高 count：深藍
        name = "log(Count)"
      ) + 
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "red",
        linewidth = 1
      ) +
      geom_point(
        data = df_label,
        shape = 21, fill = "yellow", color = "black", size = 2
      ) + 
      geom_label_repel(
        data = df_label,
        aes(label = gene),
        fill = "white",
        color = "black",
        size = 3,
        box.padding = 0.5,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 123
      ) + 
      theme_classic() + 
      ggtitle(sprintf("%s, r=%s", time_, round(cor_, 3))) + 
      ylab("scRNA non zero prop.") +
      xlab("Xenium non zero prop.")
  }) %>% patchwork::wrap_plots(ncol = 2) -> p2

patchwork::wrap_plots(p1, p2, ncol = 2) #Figure 2E-F


