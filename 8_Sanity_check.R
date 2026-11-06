library(Seurat)
library(tidyverse)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

table(xenium_obj_filtered@meta.data$predicted.id, xenium_obj_filtered@meta.data$orig.ident) %>% 
  t() %>% .[,c(3,1,2,4,5)] %>% 
  as.data.frame() %>% 
  "colnames<-"(., c("Samples", "Cell_type", "Freq")) %>% 
  dplyr::mutate(Samples = factor(Samples, levels = c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3"))) %>% 
  dplyr::mutate(Cell_type = factor(Cell_type, levels = c("Glutamatergic", "GABAergic", "Glial", "Vascular", "Immune"))) %>% 
  ggplot(aes(x = Samples, y = Freq, fill = Cell_type)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("Glutamatergic" = "#00bf7d", 
                               "GABAergic" = "#f8766d", 
                               "Glial" = "#a3a500", 
                               "Vascular" = "#e76bf3", 
                               "Immune" = "#00b0f6")) + 
  ylab("Number of cells") + 
  theme_classic() -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/cell_type_counts_per_sample.pdf", width = 8.86, height = 6.26)

#Coronal slide
unique(xenium_obj_filtered@meta.data$orig.ident) %>% 
  as.character() %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(CT_){
    xenium_obj_filtered@meta.data %>% dplyr::filter(orig.ident == CT_) -> df_
    df_ %>% rownames_to_column("cell") %>% 
      dplyr::select(cell, predicted.id) -> df_
    xy_coordinates[[CT_]] -> xy
    left_join(x = df_, y = xy, by = "cell") -> df_
    
    df_ %>% 
      ggplot(aes(x = x, y = y, color = predicted.id)) + 
      ggrastr::geom_point_rast(size = 0.1) + 
      theme_classic() + 
      ggtitle(CT_) -> p
    
  }) -> p_list
patchwork::wrap_plots(p_list, ncol = 3, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/xenium_coronal_celltypes.pdf", height = 6.82, width = 14.4)
rm(p, p_list);gc()

#Compare label transfer with mapmycells result
read_csv("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/mapmycells/output_annotation/xenium_wmb_annotation_hierarchical.csv", skip = 4) -> wmb_cell_annotation
table(wmb_cell_annotation$class_name) %>% names() -> class_
class_list = list(
  Glutamatergic = class_[grepl("Glut", class_)], 
  GABAergic = class_[grepl("GABA", class_)], 
  Other_neuronal = class_[grepl("21|22", class_)], 
  Glial = class_[grepl("30|31|32", class_)],
  Vascular = class_[grepl("33", class_)],
  Immune = class_[grepl("34", class_)]
)

wmb_cell_annotation %>% 
  dplyr::select(cell_id, class_name, class_bootstrapping_probability) %>% 
  dplyr::mutate(custom_class = case_when(
    class_name %in% class_list$Glutamatergic ~ "Glutamatergic",
    class_name %in% class_list$GABAergic ~ "GABAergic",
    class_name %in% class_list$Other_neuronal ~ "Other_neuronal",
    class_name %in% class_list$Glial ~ "Glial",
    class_name %in% class_list$Vascular ~ "Vascular",
    class_name %in% class_list$Immune ~ "Immune",
    TRUE ~ class_name
  )) -> wmb_cell_annotation_1
dim(wmb_cell_annotation_1) #683716

wmb_cell_annotation_1 %>% dplyr::filter(class_bootstrapping_probability >= 0.8) -> wmb_cell_annotation_1
dim(wmb_cell_annotation_1) #413607
rm(wmb_cell_annotation, class_, class_list);gc()

intersect(colnames(xenium_obj_filtered), wmb_cell_annotation_1$cell_id) -> cells_to_keep
wmb_cell_annotation_1 %>% column_to_rownames("cell_id") %>% .[cells_to_keep, ] -> wmb_cell_annotation_1
xenium_obj_filtered@meta.data[cells_to_keep, c("predicted.id", "predicted.id.score")] -> rpca_cell_annotation_1
cbind(rpca_cell_annotation_1, wmb_cell_annotation_1) %>% 
  {
    df_ = .
    table(df_$predicted.id, df_$custom_class)
  }
rm(cells_to_keep, wmb_cell_annotation_1, rpca_cell_annotation_1);gc()

#Find celltype specific marker genes ----
Idents(xenium_obj_filtered) <- xenium_obj_filtered$predicted.id
markers <- FindAllMarkers(
  xenium_obj_filtered,
  only.pos        = TRUE,           # markers enriched in each cell type
  test.use        = "wilcox",       # fast & standard; or "MAST" for covariates
  min.pct         = 0.01,            # gene expressed in ≥10% cells of a group
  logfc.threshold = 0.1,           # ~1.2x fold-change
  return.thresh   = 0.05            # FDR filter (p_val_adj)
)

save(markers, file="~/Dropbox/Spatial_transcriptomics_new/rda/xenium_celltype_specific_markers.rda")
gc()

levels(markers$cluster) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(cluster_){
    markers %>% dplyr::filter(cluster == cluster_) %>% .$gene %>% .[1:20]
  })

levels(markers$cluster) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(cluster_){
    markers %>% dplyr::filter(cluster == cluster_) %>% dplyr::arrange(p_val_adj) -> df_
#    df_[df_$gene %in% c("Slc17a7","Gad1", "Olig1", "Kcnj10", "Flt1", "Cx3cr1"), ] %>% drop_na() -> df_
    df_[df_$gene %in% c("Slc17a7", "Nrn1","Gad1", "Slc32a1", "Olig1", "Kcnj10", "Flt1", "Cldn5", "Cx3cr1", "Tmem119"), ] %>% drop_na() -> df_
    df_$gene
  }) %>% unlist() %>% unname() -> genes_selected

unique(xenium_obj_filtered@meta.data$predicted.id) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(predicted.id_){
    xenium_obj_filtered@meta.data %>% dplyr::filter(predicted.id == predicted.id_) -> meta_data
    data.frame(cell_type = predicted.id_, nCount = meta_data$nCount_Xenium)
  }) %>% do.call(rbind, .) %>% 
  "rownames<-"(., NULL) %>% 
  dplyr::mutate(cell_type = factor(cell_type, levels = rev(c("Glutamatergic", "GABAergic", "Glial", "Vascular", "Immune")))) %>% 
  ggplot(aes(y = cell_type, x = nCount, fill = cell_type)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) + 
  theme_classic() + 
  theme(legend.position = "none") -> p_nCount
p_nCount + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_nCount_1

unique(xenium_obj_filtered@meta.data$predicted.id) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(predicted.id_){
    xenium_obj_filtered@meta.data %>% dplyr::filter(predicted.id == predicted.id_) -> meta_data
    data.frame(cell_type = predicted.id_, nFeature = meta_data$nFeature_Xenium)
  }) %>% do.call(rbind, .) %>% 
  "rownames<-"(., NULL) %>% 
  dplyr::mutate(cell_type = factor(cell_type, levels = rev(c("Glutamatergic", "GABAergic", "Glial", "Vascular", "Immune")))) %>% 
  ggplot(aes(y = cell_type, x = nFeature, fill = cell_type)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) + 
  theme_classic() + 
  theme(legend.position = "none") -> p_nFeature
p_nFeature + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_nFeature_1

xenium_obj_filtered@meta.data -> meta_data
mat_selected = GetAssayData(object = xenium_obj_filtered, assay = "Xenium", layer = "counts")[genes_selected, ]
dim(mat_selected)

unique(meta_data$predicted.id) %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(predicted.id_){
    meta_data %>% dplyr::filter(predicted.id == predicted.id_) -> df_
    df_$nCount_Xenium %>% sum()
  }) %>% unlist() %>% median() -> median_

unique(meta_data$predicted.id) %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(predicted.id_){
    meta_data %>% dplyr::filter(predicted.id == predicted.id_) -> df_
    df_$nCount_Xenium %>% sum() -> lib_size
    cells_selected = rownames(df_)
    mat_selected[, cells_selected] %>% rowSums() %>% {(./lib_size)*median_} %>% 
      as.data.frame() %>% "colnames<-"(., "expr") %>% rownames_to_column("gene") %>% 
      dplyr::mutate(cell_type = predicted.id_, .before = 1) -> df_1
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_expr

unique(meta_data$predicted.id) %>% 
  "names<-"(.,.) %>% 
  #    .[1] %>% 
  purrr::map(function(predicted.id_){
    meta_data %>% dplyr::filter(predicted.id == predicted.id_) -> df_
    cells_selected = rownames(df_)
    n_cells = nrow(df_)
    mat_ = mat_selected[, cells_selected] 
    (mat_ != 0) %>% rowSums() %>% {./n_cells} %>% 
      as.data.frame() %>% "colnames<-"(., "non_zero_prop") %>% 
      rownames_to_column("gene") %>% 
      dplyr::mutate(cell_type = predicted.id_, .before = 1) -> df_1
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_non_zero_prop

full_join(x = df_expr, y = df_non_zero_prop, by = c("cell_type", "gene")) -> df_final
unique(df_final$gene) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    df_final %>% dplyr::filter(gene == gene_) %>% 
      dplyr::mutate(scaled_expr = scales::rescale(expr, to = c(0,1))) %>% 
      dplyr::mutate(non_zero_prop = non_zero_prop*100)
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_plot

df_plot$cell_type = factor(df_plot$cell_type, levels = rev(c("Glutamatergic", "GABAergic", "Glial", "Vascular", "Immune")))
df_plot$gene = factor(df_plot$gene, levels = genes_selected)
df_plot %>% 
  ggplot(aes(x = gene, y = cell_type)) +
  geom_point(aes(size = non_zero_prop, color = scaled_expr)) + 
  scale_size(range = c(0,15), limits = c(0,100)) + 
  scale_color_gradient(low = "grey", high = "blue") + 
  theme_classic() -> p_dotplot
p_dotplot + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_dotplot_1
save(p_dotplot_1, file = "~/Downloads/celltype_marker_genes.rda")
gc()

patchwork::wrap_plots(p_nCount, p_nFeature_1, p_dotplot_1, nrow = 1) + 
  patchwork::plot_layout(widths = c(1,1,3)) -> p_dotplot_2
ggsave(plot = p_dotplot_2, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/celltype_specfic_marker_genes_dotplot.pdf", width = 12.1, height = 4.68)

rm(df_expr, df_final, df_non_zero_prop, df_plot, markers, mat_selected, meta_data, p_dotplot, p_dotplot_1, p_dotplot_2, p_nCount, p_nFeature, p_nCount_1, p_nFeature_1, genes_selected, median_)
gc()

load("~/Downloads/cluster_markers_list.rda")
xenium_obj_filtered@meta.data$cluster

levels(xenium_obj_filtered@meta.data$cluster) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(cluster_){
    xenium_obj_filtered@meta.data %>% dplyr::filter(cluster == cluster_) -> meta_data
    data.frame(cell_type = cluster_, nCount = meta_data$nCount_Xenium)
  }) %>% do.call(rbind, .) %>% 
  "rownames<-"(., NULL) %>% 
  dplyr::mutate(cell_type = factor(cell_type, levels = rev(c("1", "4", "10", "2", "7", "9", "13", "3", "5", "8", "6", "11", "12")))) %>% 
  ggplot(aes(y = cell_type, x = nCount, fill = cell_type)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) + 
  scale_fill_manual(values = c("1" = "#00bf7d", "4" = "#00bf7d", "10" = "#00bf7d", "2" = "#f8766d", "7" = "#f8766d", "9" = "#f8766d"
                               , "13" = "#f8766d", "3" = "#a3a500", "5" = "#a3a500", "8" = "#a3a500", "6" = "#e76bf3", 
                               "11" = "#00b0f6", "12" = "#00b0f6")) +
  theme_classic() + 
  theme(legend.position = "none") -> p_nCount
p_nCount + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_nCount_1

levels(xenium_obj_filtered@meta.data$cluster) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(cluster_){
    xenium_obj_filtered@meta.data %>% dplyr::filter(cluster == cluster_) -> meta_data
    data.frame(cell_type = cluster_, nFeature = meta_data$nFeature_Xenium)
  }) %>% do.call(rbind, .) %>% 
  "rownames<-"(., NULL) %>% 
  dplyr::mutate(cell_type = factor(cell_type, levels = rev(c("1", "4", "10", "2", "7", "9", "13", "3", "5", "8", "6", "11", "12")))) %>% 
  ggplot(aes(y = cell_type, x = nFeature, fill = cell_type)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) + 
  scale_fill_manual(values = c("1" = "#00bf7d", "4" = "#00bf7d", "10" = "#00bf7d", "2" = "#f8766d", "7" = "#f8766d", "9" = "#f8766d"
                               , "13" = "#f8766d", "3" = "#a3a500", "5" = "#a3a500", "8" = "#a3a500", "6" = "#e76bf3", 
                               "11" = "#00b0f6", "12" = "#00b0f6")) +
  theme_classic() + 
  theme(legend.position = "none") -> p_nFeature
p_nFeature + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_nFeature_1

xenium_obj_filtered@meta.data -> meta_data
genes_selected = c("Slc17a7", "Oxt", "Slc17a6", "Gad1", "Ppp1r1b", "Vip", "Dlx1", 
                   "Fgfr3", "Mag", "Tmem212", "Flt1", "Mrc1", "Siglech")
mat_selected = GetAssayData(object = xenium_obj_filtered, assay = "Xenium", layer = "counts")[genes_selected, ]
dim(mat_selected)

levels(meta_data$cluster) %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(cluster_){
    meta_data %>% dplyr::filter(cluster == cluster_) -> df_
    df_$nCount_Xenium %>% sum()
  }) %>% unlist() %>% median() -> median_

levels(meta_data$cluster) %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(cluster_){
    meta_data %>% dplyr::filter(cluster == cluster_) -> df_
    df_$nCount_Xenium %>% sum() -> lib_size
    cells_selected = rownames(df_)
    mat_selected[, cells_selected] %>% rowSums() %>% {(./lib_size)*median_} %>% 
      as.data.frame() %>% "colnames<-"(., "expr") %>% rownames_to_column("gene") %>% 
      dplyr::mutate(cell_type = cluster_, .before = 1) -> df_1
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_expr

levels(meta_data$cluster) %>% 
  "names<-"(.,.) %>% 
  #    .[1] %>% 
  purrr::map(function(cluster_){
    meta_data %>% dplyr::filter(cluster == cluster_) -> df_
    cells_selected = rownames(df_)
    n_cells = nrow(df_)
    mat_ = mat_selected[, cells_selected] 
    (mat_ != 0) %>% rowSums() %>% {./n_cells} %>% 
      as.data.frame() %>% "colnames<-"(., "non_zero_prop") %>% 
      rownames_to_column("gene") %>% 
      dplyr::mutate(cell_type = cluster_, .before = 1) -> df_1
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_non_zero_prop

full_join(x = df_expr, y = df_non_zero_prop, by = c("cell_type", "gene")) -> df_final
unique(df_final$gene) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    df_final %>% dplyr::filter(gene == gene_) %>% 
      dplyr::mutate(scaled_expr = scales::rescale(expr, to = c(0,1))) %>% 
      dplyr::mutate(non_zero_prop = non_zero_prop*100)
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_plot

df_plot$cell_type = factor(df_plot$cell_type, levels = rev(c("1", "4", "10", "2", "7", "9", "13", "3", "5", "8", "6", "11", "12")))
df_plot$gene = factor(df_plot$gene, levels = genes_selected)
df_plot %>% 
  ggplot(aes(x = gene, y = cell_type)) +
  geom_point(aes(size = non_zero_prop, color = scaled_expr)) + 
  scale_size(range = c(0,15), limits = c(0,100)) + 
  scale_color_gradient(low = "grey", high = "blue") + 
  theme_classic() -> p_dotplot
p_dotplot + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_dotplot_1
gc()

patchwork::wrap_plots(p_nCount, p_nFeature_1, p_dotplot_1, nrow = 1) + 
  patchwork::plot_layout(widths = c(1,1,3)) -> p_dotplot_2
ggsave(plot = p_dotplot_2, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/celltype_subclass_specfic_marker_genes_dotplot.pdf", width = 12.38, height = 6.05) #Figure 1E

#Show SCN markers
SCN_markers = c("Vip", "Avp", "Grp", "Nms")
mat_selected = GetAssayData(object = xenium_obj_filtered, assay = "Xenium", layer = "counts")[SCN_markers, ]
apply(mat_selected, 1, function(x){x[x!=0]} %>% names()) -> cells_with_markers
setdiff(
  x = colnames(xenium_obj_filtered),
  y = (unlist(cells_with_markers) %>% unname() %>% unique())
  ) -> cells_wo_markers
cells_with_markers[["Others"]] <- cells_wo_markers
rm(cells_wo_markers)


unique(xenium_obj_filtered@meta.data$orig.ident) %>% 
  as.character() %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(CT_){
    xenium_obj_filtered@meta.data %>% dplyr::filter(orig.ident == CT_) -> df_
    cells_to_keep = rownames(df_)
    xy_coordinates[[CT_]] -> xy
    xy = xy %>% dplyr::filter(cell %in% cells_to_keep)
    cells_with_markers %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(cells_, gene_){
        xy %>% dplyr::filter(cell %in% cells_) %>% 
          dplyr::mutate(marker = gene_)
      }) %>% do.call(rbind, .)
  }) -> df_plot

df_plot %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, CT_){
    df_ %>% 
      dplyr::mutate(marker = factor(marker, levels = c("Vip", "Avp", "Others"))) %>% 
      ggplot(aes(x = x, y = y, color = marker)) + 
      ggrastr::geom_jitter_rast(size = 0.1, alpha = 0.2) + 
      theme_classic() + 
      ggtitle(CT_) -> p
    p + scale_color_manual(values = c(Vip = "blue", Avp = "red", Others = "grey"))
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 3, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_marker_dimplot.pdf", 
       height = 8.49, width = 14.96)

p_list %>% 
  purrr::map(function(p_){
    p_ + 
      coord_cartesian(xlim = c(2500,7500), ylim = c(300, 1500)) -> p
  }) -> p_list_1

patchwork::wrap_plots(p_list_1, ncol = 3, guides = "collect") -> p_1
ggsave(plot = p_1, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_marker_dimplot_zoomed.pdf", 
       height = 8.49, width = 14.96)

rm(cells_with_markers, df_plot, mat_selected, p, p_1, p_list, p_list_1, SCN_markers);gc()
rm(list=ls())
gc()

library(tidyverse)
library(Seurat)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

#list.files(path = "~/Dropbox/Spatial_transcriptomics_new/rda", pattern = "VIP_enriched_regions\\.csv", full.names = T) %>% 
list.files(path = "~/Dropbox/Spatial_transcriptomics_new/rda", pattern = "SCN_regions\\.csv", full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)_SCN.+", "\\1", .) %>% gsub("_", "-", .)) %>% 
  purrr::map(function(file_){
    read_csv(file_, skip = 2) -> df_
  })-> SCN_polygon

make_spatial_dimplot = function(xenium_obj_filtered, xy_coordinates, markers){
  
  mat_selected = GetAssayData(object = xenium_obj_filtered, assay = "Xenium", layer = "counts")[markers, ]
  apply(mat_selected, 1, function(x){x[x!=0]} %>% names()) -> cells_with_markers
  setdiff(
    x = colnames(xenium_obj_filtered),
    y = (unlist(cells_with_markers) %>% unname() %>% unique())
  ) -> cells_wo_markers
  cells_with_markers[["Others"]] <- cells_wo_markers
  rm(cells_wo_markers)
  
  
  unique(xenium_obj_filtered@meta.data$orig.ident) %>% 
    as.character() %>% 
    "names<-"(.,.) %>% 
    purrr::map(function(CT_){
      xenium_obj_filtered@meta.data %>% dplyr::filter(orig.ident == CT_) -> df_
      cells_to_keep = rownames(df_)
      xy_coordinates[[CT_]] -> xy
      xy = xy %>% dplyr::filter(cell %in% cells_to_keep)
      cells_with_markers %>% 
        purrr::map2(.x=.,.y=names(.),.f=function(cells_, gene_){
          xy %>% dplyr::filter(cell %in% cells_) %>% 
            dplyr::mutate(marker = gene_)
        }) %>% do.call(rbind, .)
    }) -> df_plot
  
  df_plot %>% 
    purrr::map2(.x=.,.y=names(.),.f=function(df_, CT_){
      df_ %>% 
        dplyr::mutate(group = case_when(
          marker == "Others" ~ 1, 
          TRUE ~ 2
        )) %>% dplyr::arrange(group) -> df_
      df_ %>% 
        dplyr::mutate(marker = factor(marker, levels = c("Others", "Vip", "Avp", "Grp", "Nms"))) %>% 
        ggplot(aes(x = x, y = y, color = marker)) + 
#        geom_jitter(size = 0.2, alpha = 0.5) + 
        ggrastr::geom_jitter_rast(size = 0.2, alpha = 0.5) + 
#        ggrastr::geom_point_rast(size = 0.5, alpha = 1) + 
        theme_classic() + 
        ggtitle(CT_) -> p
      p + scale_color_manual(values = c(Others = "grey", Vip = "blue", Avp = "red", Grp = "green", Nms = "orange"))
    }) -> p_list 
  return(p_list)
}

SCN_markers = c("Vip", "Avp", "Grp", "Nms")
SCN_Dimplot_p_list = make_spatial_dimplot(xenium_obj_filtered = xenium_obj_filtered, xy_coordinates = xy_coordinates, markers = SCN_markers)

patchwork::wrap_plots(SCN_Dimplot_p_list, ncol = 3, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_marker_dimplot.pdf", 
       height = 8.49, width = 14.96)

SCN_Dimplot_p_list %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(p_, sample_){
    p_ + 
#      geom_polygon(data = SCN_polygon[[sample_]], aes(x = X, y = Y), fill = NA, color = "black", alpha = 0.5) + 
      coord_cartesian(xlim = c(2500,7500), ylim = c(200, 1500)) -> p
  }) -> SCN_Dimplot_p_list_1

patchwork::wrap_plots(SCN_Dimplot_p_list_1, ncol = 3, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_polygon_selection.pdf", height = 8.49, width = 14.96)

rm(SCN_Dimplot_p_list_1, p);gc()

SCN_Dimplot_p_list_1 %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(p_, sample_){
    SCN_polygon[[sample_]]$Selection %>% unique() -> selection_
    for (x in selection_){
      SCN_polygon[[sample_]] %>% 
        dplyr::filter(Selection == x) -> polygon_df_
      p_ = p_ + geom_polygon(data = polygon_df_, aes(x = X, y = Y), fill = NA, color = "black", alpha = 0.5)
    }
    return(p_)
  }) -> SCN_Dimplot_p_list_2

patchwork::wrap_plots(SCN_Dimplot_p_list_2, ncol = 3, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_polygon_selection_1.pdf", height = 8.49, width = 14.96)

rm(SCN_Dimplot_p_list_1, p);gc()

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  purrr::map(function(group_){
    xy_coordinates[[group_]] -> xy_coordinates
    unique(SCN_polygon[[group_]]$Selection) %>% 
      "names<-"(.,.) %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(region_, group__){
        SCN_polygon[[group_]] %>% dplyr::filter(Selection == region_) %>% dplyr::select(X,Y) -> df_
        cells_se = sp::point.in.polygon(
          point.x = xy_coordinates$x,
          point.y = xy_coordinates$y,
          pol.x = df_$X,
          pol.y = df_$Y) > 0
        xy_coordinates[cells_se, ] -> df_1
        df_1 %>% dplyr::mutate(group = group__) %>% 
          dplyr::mutate(group = case_when(
            grepl("AVP", group) ~ "Shell",
            TRUE ~ "Core"
          )) -> df_1
        df_1 %>% dplyr::select(-x, -y) -> df_1
      }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_final
    print(sprintf("%s, n_cells = %s", group_, length(df_final$cell)))
    return(df_final)
  }) %>% do.call(rbind, .) -> SCN_cells
SCN_cells %>% dplyr::filter(cell %in% colnames(xenium_obj_filtered)) -> SCN_cells

c("SCN", "Core", "Shell") %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  purrr::map(function(region_){
    if (region_ == "SCN"){
      df_ = SCN_cells
    }else{
      df_ = SCN_cells %>% dplyr::filter(group == region_)
    }
    SCN_cells = df_$cell
    xenium_obj_filtered@meta.data[SCN_cells,]$predicted.id %>% 
    table() %>% 
      as.data.frame() %>% "colnames<-"(., c("celltype", "n_cells")) %>% 
      dplyr::mutate(Percent = n_cells/sum(n_cells)*100) %>% 
      dplyr::arrange(desc(n_cells)) %>% 
      dplyr::mutate(celltype = factor(celltype, levels = c("GABAergic", "Glutamatergic", "Glial", "Vascular", "Immune"))) %>% 
      ggplot(aes(y = n_cells, x = celltype, fill = celltype)) + 
      geom_bar(stat="identity") + 
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ggrepel::geom_text_repel(aes(label = paste0(round(Percent,2), "%")), 
                               position = position_stack(vjust = 1), color = "black") -> p
    p + ggtitle(region_) -> p
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 3, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_region_celltype_composition.pdf")

c("SCN", "Core", "Shell") %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(region_){
    if (region_ == "SCN"){
      df_ = SCN_cells
    }else{
      df_ = SCN_cells %>% dplyr::filter(group == region_)
    }
    SCN_cells = df_$cell
    xenium_obj_filtered@meta.data[SCN_cells,]$cluster %>% 
      table() %>% 
      as.data.frame() %>% "colnames<-"(., c("celltype", "n_cells")) %>% 
      dplyr::mutate(Percent = n_cells/sum(n_cells)*100) %>% 
      dplyr::arrange(desc(n_cells)) %>% 
      dplyr::mutate(celltype = factor(celltype, levels = c("1", "4", "10", "2", "7", "9", "13", "3", "5", "8", "6", "11", "12"))) %>% 
      ggplot(aes(y = n_cells, x = celltype, fill = celltype)) + 
      geom_bar(stat="identity") + 
      scale_fill_manual(values = c("1" = "#00bf7d", "4" = "#00bf7d", "10" = "#00bf7d", "2" = "#f8766d", "7" = "#f8766d", "9" = "#f8766d"
                                   , "13" = "#f8766d", "3" = "#a3a500", "5" = "#a3a500", "8" = "#a3a500", "6" = "#e76bf3", 
                                   "11" = "#00b0f6", "12" = "#00b0f6")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ggrepel::geom_text_repel(aes(label = paste0(round(Percent,2), "%")), 
                               position = position_stack(vjust = 1), color = "black") -> p
    p + ggtitle(region_) -> p
  }) -> p_list
patchwork::wrap_plots(p_list, nrow = 3, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_region_cluster_composition.pdf", height = 10.84, width = 7.94)


xenium_obj_filtered@meta.data %>% dplyr::mutate(group = case_when(
  rownames(.) %in% SCN_cells ~ "SCN_region",
  TRUE ~ "Other_regions"
)) -> xenium_obj_filtered@meta.data

levels(xenium_obj_filtered$orig.ident)
xenium_obj_filtered@meta.data %>% dplyr::filter(group == "SCN_region") %>% 
  {
    df_ = .
    c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      purrr::map(function(CT_){
        df_ %>% dplyr::filter(orig.ident == CT_) %>% dplyr::select(orig.ident, predicted.id) -> meta_data
        cell_ids = rownames(meta_data)
        
        xy_coordinates[[CT_]] -> xy_coordinates
        rownames(xy_coordinates) <- xy_coordinates$cell
        xy_coordinates[cell_ids, ] -> xy_coordinates
        cbind(meta_data, xy_coordinates[, c("x", "y")]) -> df_plot
        df_plot %>% 
          ggplot(aes(x = x, y = y, color = predicted.id)) +
          geom_point(size = 0.5) + 
          geom_polygon(data = SCN_polygon[[CT_]], aes(x = X, y = Y), fill = NA, color = "black") + 
          theme_classic() + 
          ggtitle(CT_)
      })
  } -> p_list
patchwork::wrap_plots(p_list, ncol = 3, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SCN_region_celltype_spatial_distribution.pdf")


####
library(tidyverse)
library(Seurat)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

#list.files(path = "~/Dropbox/Spatial_transcriptomics_new/rda", pattern = "VIP_enriched_regions\\.csv", full.names = T) %>% 
list.files(path = "~/Dropbox/Spatial_transcriptomics_new/rda", pattern = "SON_PVN_regions\\.csv", full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)_SON.+", "\\1", .) %>% gsub("_", "-", .)) %>% 
  purrr::map(function(file_){
    read_csv(file_, skip = 2) -> df_
    df_ %>% dplyr::filter(grepl("SON", Selection)) -> df_
  })-> SON_polygon

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(group_){
    xy_coordinates[[group_]] -> xy_coordinates
    unique(SON_polygon[[group_]]$Selection) %>% 
      "names<-"(.,.) %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(region_, group__){
        SON_polygon[[group_]] %>% dplyr::filter(Selection == region_) %>% dplyr::select(X,Y) -> df_
        cells_se = sp::point.in.polygon(
          point.x = xy_coordinates$x,
          point.y = xy_coordinates$y,
          pol.x = df_$X,
          pol.y = df_$Y) > 0
        xy_coordinates[cells_se, ] -> df_1
        df_1 %>% dplyr::mutate(group = group__) %>% 
          dplyr::mutate(group = case_when(
            grepl("AVP", group) ~ "Shell",
            TRUE ~ "Core"
          )) -> df_1
        df_1 %>% dplyr::select(-x, -y) -> df_1
      }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_final
    print(sprintf("%s, n_cells = %s", group_, length(df_final$cell)))
    return(df_final)
  }) %>% do.call(rbind, .) -> SON_cells
SON_cells %>% dplyr::filter(cell %in% colnames(xenium_obj_filtered)) -> SON_cells

xenium_obj_filtered@meta.data[SON_cells$cell,]$predicted.id %>% 
  table() %>% 
  as.data.frame() %>% "colnames<-"(., c("celltype", "n_cells")) %>% 
  dplyr::mutate(Percent = n_cells/sum(n_cells)*100) %>% 
  dplyr::arrange(desc(n_cells)) %>% 
  dplyr::mutate(celltype = factor(celltype, levels = c("Glutamatergic", "GABAergic", "Glial", "Vascular", "Immune"))) %>% 
  ggplot(aes(y = n_cells, x = celltype, fill = celltype)) + 
  geom_bar(stat="identity") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggrepel::geom_text_repel(aes(label = paste0(round(Percent,2), "%")), 
                           position = position_stack(vjust = 1), color = "black") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SON_region_celltype_composition.pdf")

xenium_obj_filtered@meta.data[SON_cells$cell,]$cluster %>% 
  table() %>% 
  as.data.frame() %>% "colnames<-"(., c("celltype", "n_cells")) %>% 
  dplyr::mutate(Percent = n_cells/sum(n_cells)*100) %>% 
  dplyr::arrange(desc(n_cells)) %>% 
  dplyr::mutate(celltype = factor(celltype, levels = c("1", "4", "10", "2", "7", "9", "13", "3", "5", "8", "6", "11", "12"))) %>% 
  ggplot(aes(x = celltype, y = n_cells, fill = celltype)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("1" = "#00bf7d", "4" = "#00bf7d", "10" = "#00bf7d", "2" = "#f8766d", "7" = "#f8766d", "9" = "#f8766d"
                               , "13" = "#f8766d", "3" = "#a3a500", "5" = "#a3a500", "8" = "#a3a500", "6" = "#e76bf3", 
                               "11" = "#00b0f6", "12" = "#00b0f6")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggrepel::geom_text_repel(aes(label = paste0(round(Percent,2), "%")), 
                           position = position_stack(vjust = 1), color = "black") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/SON_region_cluster_composition.pdf")

list.files(path = "~/Dropbox/Spatial_transcriptomics_new/rda", pattern = "SON_PVN_regions\\.csv", full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)_SON.+", "\\1", .) %>% gsub("_", "-", .)) %>% 
  purrr::map(function(file_){
    read_csv(file_, skip = 2) -> df_
    df_ %>% dplyr::filter(grepl("PVN", Selection)) -> df_
  })-> PVN_polygon

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(group_){
    xy_coordinates[[group_]] -> xy_coordinates
    unique(PVN_polygon[[group_]]$Selection) %>% 
      "names<-"(.,.) %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(region_, group__){
        PVN_polygon[[group_]] %>% dplyr::filter(Selection == region_) %>% dplyr::select(X,Y) -> df_
        cells_se = sp::point.in.polygon(
          point.x = xy_coordinates$x,
          point.y = xy_coordinates$y,
          pol.x = df_$X,
          pol.y = df_$Y) > 0
        xy_coordinates[cells_se, ] -> df_1
        df_1 %>% dplyr::mutate(group = group__) %>% 
          dplyr::mutate(group = case_when(
            grepl("AVP", group) ~ "Shell",
            TRUE ~ "Core"
          )) -> df_1
        df_1 %>% dplyr::select(-x, -y) -> df_1
      }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_final
    print(sprintf("%s, n_cells = %s", group_, length(df_final$cell)))
    return(df_final)
  }) %>% do.call(rbind, .) -> PVN_cells
PVN_cells %>% dplyr::filter(cell %in% colnames(xenium_obj_filtered)) -> PVN_cells

xenium_obj_filtered@meta.data[PVN_cells$cell,]$predicted.id %>% 
  table() %>% 
  as.data.frame() %>% "colnames<-"(., c("celltype", "n_cells")) %>% 
  dplyr::mutate(Percent = n_cells/sum(n_cells)*100) %>% 
  dplyr::arrange(desc(n_cells)) %>% 
  dplyr::mutate(celltype = factor(celltype, levels = c("Glutamatergic", "GABAergic", "Glial", "Vascular", "Immune"))) %>% 
  ggplot(aes(y = n_cells, x = celltype, fill = celltype)) + 
  geom_bar(stat="identity") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggrepel::geom_text_repel(aes(label = paste0(round(Percent,2), "%")), 
                           position = position_stack(vjust = 1), color = "black") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/PVN_region_celltype_composition.pdf")

xenium_obj_filtered@meta.data[PVN_cells$cell,]$cluster %>% 
  table() %>% 
  as.data.frame() %>% "colnames<-"(., c("celltype", "n_cells")) %>% 
  dplyr::mutate(Percent = n_cells/sum(n_cells)*100) %>% 
  dplyr::arrange(desc(n_cells)) %>% 
  dplyr::mutate(celltype = factor(celltype, levels = c("1", "4", "10", "2", "7", "9", "13", "3", "5", "8", "6", "11", "12"))) %>% 
  ggplot(aes(x = celltype, y = n_cells, fill = celltype)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("1" = "#00bf7d", "4" = "#00bf7d", "10" = "#00bf7d", "2" = "#f8766d", "7" = "#f8766d", "9" = "#f8766d"
                               , "13" = "#f8766d", "3" = "#a3a500", "5" = "#a3a500", "8" = "#a3a500", "6" = "#e76bf3", 
                               "11" = "#00b0f6", "12" = "#00b0f6")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggrepel::geom_text_repel(aes(label = paste0(round(Percent,2), "%")), 
                           position = position_stack(vjust = 1), color = "black") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/PVN_region_cluster_composition.pdf")
