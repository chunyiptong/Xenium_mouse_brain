library(tidyverse)

list.files("~/Dropbox/Spatial_transcriptomics_new/rda", pattern = "region_annotated\\.csv", recursive = T, full.names = T) %>% 
  gtools::mixedsort() %>% 
  "names<-"(., gsub("/.+/(.+)_region.+", "\\1", .) %>% gsub("_", "-", .)) %>% 
  purrr::map(function(file_){
    read_csv(file_) -> df_
  }) -> region_annotated_list

names(region_annotated_list) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(CT_){
    region_annotated_list[[CT_]] -> region_annotated
    region_annotated %>% dplyr::mutate(cell_id = sprintf("%s_%s", CT_, cell_id)) %>% dplyr::select(cell_id, acronym, struct_id) -> region_annotated
  }) -> region_annotated_list

save(region_annotated_list, file = "~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-region-annotation.rda")

library(jsonlite)

## 1. Load ontology JSON (Mouse, P56, structure_graph = 1)
url <- "http://api.brain-map.org/api/v2/structure_graph_download/1.json"
onto <- fromJSON(url, simplifyVector = FALSE)

# The root node is the first element
root <- onto$msg[[1]]

## 2. Recursively flatten the tree into a data.frame
flatten_structs <- function(node) {
  parent_id <- if (is.null(node$parent_structure_id)) NA_integer_ else node$parent_structure_id
  
  df <- data.frame(
    id                  = node$id,
    acronym             = node$acronym,
    name                = node$name,
    parent_structure_id = parent_id,
    stringsAsFactors    = FALSE
  )
  
  if (length(node$children) == 0) {
    return(df)
  }
  
  children_df <- do.call(
    rbind,
    lapply(node$children, flatten_structs)
  )
  
  rbind(df, children_df)
}

structs <- flatten_structs(root)
head(structs)
save(structs, file = "~/Dropbox/Spatial_transcriptomics_new/rda/Allen_brain_ontology_structs.rda")

get_parents <- function(acronym, structs_df = structs) {
  row <- structs_df[structs_df$acronym == acronym, ]
  if (nrow(row) == 0) stop("Acronym not found: ", acronym)
  
  parents <- list()
  parent_id <- row$parent_structure_id
  
  while (!is.na(parent_id) && parent_id != 0) {
    parent_row <- structs_df[structs_df$id == parent_id, ]
    if (nrow(parent_row) == 0) break
    parents <- c(parents, list(parent_row))
    parent_id <- parent_row$parent_structure_id
  }
  
  do.call(rbind, parents)  # closest parent first
}

safe_fun <- possibly(get_parents, otherwise = NULL)

region_annotated_list %>% 
  purrr::map(function(df_){
    df_$acronym
  }) %>% unlist() %>% 
  unique() %>% 
#  .[30] %>% 
  purrr::map(function(acronym_){
    onto = c("CTX", "CNU", "TH", "HY", "fiber tracts", "VS")
    if ( acronym_ %in% onto){
      parent_ = acronym_
    }else{
      safe_fun(acronym_) %>% .$acronym -> acronyms_
      safe_fun(acronym_)      
      acronyms_ %>% {.[. %in% onto]} -> parent_
    }
    
    if (is.null(parent_)){
      parent_ = NA
    }
    data.frame(acronym = acronym_, parent = parent_)
  }) %>% do.call(rbind, .) %>% 
  dplyr::arrange(parent) -> region_parent_mapping

region_annotated_list %>% 
  purrr::map(function(df_){
    left_join(x = df_, y = region_parent_mapping, by = "acronym") -> df_ 
  }) -> region_annotated_list_1

save(region_annotated_list, region_annotated_list_1, file = "~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-region-annotation.rda")

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

names(region_annotated_list_1) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(CT_){
    region_annotated_list_1[[CT_]] -> df_
    xy_coordinates[[CT_]] -> df_1
    left_join(x = df_, y = df_1, by = c("cell_id" = "cell"))
  }) -> region_annotated_list_1

region_annotated_list_1 %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, CT_){
    print(CT_)
    df_ %>%
      dplyr::filter(!is.na(parent)) %>% 
      dplyr::mutate(parent = factor(parent, levels = c("CTX", "CNU", "TH", "HY", "VS", "fiber tracts"))) %>% 
      ggplot(aes(x = x, y = y, color = parent)) + 
      ggrastr::geom_point_rast(size = 0.2, alpha = 0.8) + 
      theme_classic() +
      ggtitle(CT_)
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 3, guides = "collect") + patchwork::plot_annotation(title = "Region annotation based on Allen Brain Atlas ontology") -> p_all
p_all
ggsave(plot = p_all, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/annotation-allen-ontology.pdf", width = 19.8, height = 11.9) #Figure S3A
rm(list=ls())

library(tidyverse)
library(Seurat)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-region-annotation.rda")
rm(region_annotated_list)

xenium_obj_filtered@meta.data %>% rownames_to_column("cell_id") %>% 
  {
    df_ = .
    region_annotated_list_1 %>% 
      purrr::map(function(df_){
        df_
      }) %>% do.call(rbind, .) -> df_1
    
    left_join(x = df_, y = df_1, by = "cell_id")
  } %>% column_to_rownames("cell_id") -> meta_data_1

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")
xy_coordinates %>% 
  do.call(rbind, .) %>% 
  "rownames<-"(., NULL) -> xy_coordinates

meta_data_1 %>% rownames_to_column("cell") %>% 
  left_join(x = ., y = xy_coordinates, by = "cell") %>% 
  column_to_rownames("cell") -> meta_data_1

head(meta_data_1)

table(meta_data_1$orig.ident, meta_data_1$parent)[c(4,5,6,1,2,3), ]
table(meta_data_1$acronym, meta_data_1$orig.ident)[, c(4,5,6,1,2,3)]

xenium_obj_filtered = NormalizeData(xenium_obj_filtered, assay = "Xenium", normalization.method = "RC")
normalized_counts = GetAssayData(xenium_obj_filtered, assay = "Xenium", layer = "data")
rm(xenium_obj_filtered, xy_coordinates);gc()
rm(region_annotated_list_1);gc()

TH_markers = c("Tcf7l2")
TH_markers %in% rownames(normalized_counts)

names(region_annotated_list_1) %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  purrr::map(function(Samples_){
    meta_data_1 %>% dplyr::filter(orig.ident == Samples_) -> df_
    cells = rownames(df_)
    mat_ = normalized_counts[TH_markers, cells, drop=F] %>% t()
    mat_ = mat_ %>% as.data.frame()
    mat_ %>% "colnames<-"(., "relative_expr") -> mat_
    mat_ = cbind(mat_, df_[cells, c("x", "y")])
    mat_$relative_expr = scales::rescale(mat_$relative_expr, to = c(0, 1))
    mat_ %>% dplyr::arrange(relative_expr) -> mat_
    
    mat_ %>% 
      ggplot(aes(x = x, y = y , color = relative_expr)) +
      ggrastr::geom_point_rast(size = 0.1, alpha = 0.5) +
      scale_color_continuous(low = "lightgrey", high = "red") + 
      theme_classic() + 
      ggtitle(sprintf("Tcf7l2 in %s", Samples_)) -> p
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 3, guides = "collect") -> p_all
ggsave(plot = p_all, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/TH-marker-featureplot.pdf", width = 18.5, height = 10.3)
rm(p_list, p_all);gc()

plot_markers = function(marker, quantile_=1){
  names(region_annotated_list_1) %>% 
    "names<-"(.,.) %>% 
    purrr::map(function(Samples_){
      meta_data_1 %>% dplyr::filter(orig.ident == Samples_) -> df_
      cells = rownames(df_)
      mat_ = normalized_counts[marker, cells, drop=F] %>% t()
      mat_ = mat_ %>% as.data.frame()
      mat_ %>% "colnames<-"(., "relative_expr") -> mat_
      mat_ = cbind(mat_, df_[cells, c("x", "y")])
      mat_$relative_expr %>% {.[. != 0]} %>% quantile(., quantile_) -> thresh_
#      mat_ %>% dplyr::filter(relative_expr <= thresh_) -> mat_
      mat_ %>% dplyr::mutate(relative_expr = case_when(
        relative_expr >= thresh_ ~ thresh_,
        TRUE ~ relative_expr
      )) -> mat_
      mat_$relative_expr = scales::rescale(mat_$relative_expr, to = c(0, 1))
      mat_ %>% dplyr::arrange(relative_expr) -> mat_
      
      mat_ %>% 
        ggplot(aes(x = x, y = y , color = relative_expr)) +
        ggrastr::geom_point_rast(size = 0.1, alpha = 0.5) +
        scale_color_continuous(low = "lightgrey", high = "red") + 
        theme_classic() + 
        ggtitle(sprintf("%s in %s", marker, Samples_)) -> p
    }) -> p_list
  patchwork::wrap_plots(p_list, ncol = 3, guides = "collect") -> p_all  
  return(p_all)
}

#CTX marker

pdf("~/Dropbox/Spatial_transcriptomics_new/github/plots/CTX-markers-featureplot.pdf", width = 18.5, height = 10.3)
plot_markers(marker = "Slc17a7", quantile_ = 0.95)
plot_markers("Cux1", quantile_ = 0.95) #L2-4
plot_markers("Cux2", quantile_ = 0.95) #L2-4
plot_markers("Satb2", quantile_ = 0.95) #L2-4
plot_markers("Rorb", quantile_ = 0.95) #L4
plot_markers("Bcl11b", quantile_ = 0.95) #Deep cortex and STR
plot_markers("Tbr1", quantile_ = 0.95) #Deep cortex and STR
dev.off()

#CNU
pdf("~/Dropbox/Spatial_transcriptomics_new/github/plots/CNU-markers-featureplot.pdf", width = 18.5, height = 10.3)
plot_markers("Ppp1r1b", quantile_ = 0.95)
plot_markers("Pde10a", quantile_ = 0.95)
plot_markers("Foxp1", quantile_ = 0.95)
plot_markers("Trpc4", quantile_ = 0.95)
plot_markers("Prkcd", quantile_ = 0.95)
plot_markers("Nkx2-1", quantile_ = 0.95)
plot_markers("Lhx6", quantile_ = 0.95)
plot_markers("Foxp2", quantile_ = 0.95)
plot_markers("Npas1", quantile_ = 0.95)
dev.off()

#TH
pdf("~/Dropbox/Spatial_transcriptomics_new/github/plots/TH-markers-featureplot.pdf", width = 18.5, height = 10.3)
plot_markers(marker = "Tcf7l2", quantile_ = 0.95)
plot_markers(marker = "Gbx2", quantile_ = 0.95)
plot_markers(marker = "Esrrg", quantile_ = 0.95)
plot_markers(marker = "Six3", quantile_ = 0.95)
#plot_markers(marker = "Slc17a6", quantile_ = 0.95)
dev.off()

#HY
pdf("~/Dropbox/Spatial_transcriptomics_new/github/plots/HY-markers-featureplot.pdf", width = 18.5, height = 10.3)
plot_markers("Otp", quantile_ = 0.95)
plot_markers("Cartpt", quantile_ = 0.95)
plot_markers("Avp", quantile_ = 0.95)
plot_markers("Oxt", quantile_ = 0.95)
plot_markers("Crh", quantile_ = 0.95)
plot_markers("Trh", quantile_ = 0.95)
#plot_markers("Sst", quantile_ = 0.95)
plot_markers("Oxt", quantile_ = 0.95)
plot_markers("Vip", quantile_ = 0.95)
dev.off()

#VS
pdf("~/Dropbox/Spatial_transcriptomics_new/github/plots/VS-markers-featureplot.pdf", width = 18.5, height = 10.3)
plot_markers("Foxj1", quantile_ = 0.95)
plot_markers("Pifo", quantile_ = 0.95)
dev.off()

#Fiber tracts
pdf("~/Dropbox/Spatial_transcriptomics_new/github/plots/Fiber-markers-featureplot.pdf", width = 18.5, height = 10.3)
plot_markers("Cnp", quantile_ = 0.95)
plot_markers("Mag", quantile_ = 0.95)
plot_markers("Mog", quantile_ = 0.95)
dev.off()

library(tidyverse)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")
c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  .[1:4] %>% 
  purrr::map(function(group_){
    list.files("~/Downloads", pattern = "onto_coordinates\\.csv", full.names = T) %>% 
      {.[grepl(group_, .)]} -> file_
    polygon_ = read_csv(file_, skip = 2)
    xy_coordinates[[group_]] -> xy_coordinates
    unique(polygon_$Selection) %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(region_){
        polygon_ %>% dplyr::filter(Selection == region_) %>% dplyr::select(X,Y) -> df_
        cells_se = sp::point.in.polygon(
          point.x = xy_coordinates$x,
          point.y = xy_coordinates$y,
          pol.x = df_$X,
          pol.y = df_$Y) > 0
        xy_coordinates[cells_se, ] -> df_1
        data.frame(cell = df_1$cell, region = region_, group = group_)
      }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) %>% 
      dplyr::distinct(cell, .keep_all = T) %>% 
      dplyr::mutate(upper_onto = gsub("(.+)_\\d+", "\\1", region)) %>% 
      dplyr::arrange(upper_onto) -> df_region
    left_join(x = df_region, y = xy_coordinates, by = c("cell" = "cell")) -> df_region
  }) -> df_region_list


df_region_list %>% 
  purrr::map(function(df_){
#    head(df_)
    df_ %>% 
      ggplot(aes(x = x, y = y, color = upper_onto)) + 
      ggrastr::geom_point_rast(size = 0.1, alpha = 0.5) + 
      theme_classic()
  }) -> p_list

p_list[[1]]
p_list[[2]]
p_list[[3]]
p_list[[4]]

#####
library(tidyverse)
library(sf)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")
c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
#  .[1:4] %>% 
  purrr::map(function(group_){
    list.files("~/Downloads/Brain_region_annotation", pattern = "onto_coordinates\\.csv", full.names = T) %>% 
      {.[grepl(group_, .)]} -> file_
    polygon_ = read_csv(file_, skip = 2)
    xy_coordinates[[group_]] -> xy_coordinates
    unique(polygon_$Selection) %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(region_){
        polygon_ %>% dplyr::filter(Selection == region_) %>% dplyr::select(X,Y) -> df_
      })
  }) -> polygon_list

polygon_list %>% 
#  .[4] %>% 
  purrr::map(function(list_){
    list_ %>% 
      purrr::map(function(xy_){
        colnames(xy_) = c("x", "y")
        xy_ = xy_ %>% as.matrix()
        st_polygon(list(xy_))
      }) -> polygon_list
    names(polygon_list) -> region_names
    region_names %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      purrr::map(function(name_){
        polygon_list[[name_]] -> polygon_
        setdiff(region_names, name_) %>% 
          purrr::map(function(other_regions){
            polygon_list[[other_regions]] -> other_polygon_
            st_difference(polygon_, other_polygon_) ->> polygon_
          }) %>% .[[length(.)]] -> polygon_cleaned
        if (st_is_empty(polygon_cleaned)){
          polygon_cleaned = polygon_list[[name_]]
        }else{
          polygon_cleaned = polygon_cleaned
        }
        return(polygon_cleaned)
      })
  }) -> polygon_cleaned_list


polygon_cleaned_list %>% 
#  .[4] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(list_, CT_){
    list_ %>% 
#      .[3] %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(polygon_, region_){
        polygon_[[1]] %>% as.data.frame() %>% "colnames<-"(c("X", "Y")) -> df_
        xy_coordinates[[CT_]] -> xy_coordinates
        cells_se = sp::point.in.polygon(
          point.x = xy_coordinates$x,
          point.y = xy_coordinates$y,
          pol.x = df_$X,
          pol.y = df_$Y) > 0
          xy_coordinates[cells_se, ] -> df_1
          data.frame(cell = df_1$cell, region = region_, group = CT_, x = df_1$x, y = df_1$y)
      }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) %>% 
      dplyr::mutate(upper_onto = gsub("(.+)_\\d+", "\\1", region))
  }) -> df_region_list

polygon_cleaned_list %>% 
#    .[4] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(list_, CT_){
    xy_coordinates[[CT_]] -> xy_coordinates
    list_ %>% 
      #      .[3] %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(polygon_, region_){
        polygon_ %>% 
          purrr::map(function(polygon_1){
            polygon_1 %>% as.data.frame() %>% "colnames<-"(c("X", "Y")) -> df_
            cells_se = sp::point.in.polygon(
              point.x = xy_coordinates$x,
              point.y = xy_coordinates$y,
              pol.x = df_$X,
              pol.y = df_$Y) > 0
            xy_coordinates[cells_se, ] -> df_1
          }) %>% do.call(rbind, .) -> df_1
        df_1$cell %>% table() %>% {.[.==1]} %>% names() -> cell_se
        df_1 = df_1 %>% dplyr::filter(cell %in% cell_se)
        data.frame(cell = df_1$cell, region = region_, group = CT_, x = df_1$x, y = df_1$y)
      }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) %>% 
      dplyr::mutate(upper_onto = gsub("(.+)_\\d+", "\\1", region))
  }) -> df_region_list
save(df_region_list, file = "~/Downloads/Brain_region_annotation/df_region_list.rda")

df_region_list %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, CT_){
    #    head(df_)
    df_ %>% 
      ggplot(aes(x = x, y = y, color = upper_onto)) + 
      ggrastr::geom_point_rast(size = 0.1, alpha = 0.5) + 
      theme_classic() + 
      ggtitle(CT_)
  }) -> p_list

p_list[[1]]
p_list[[2]]
p_list[[3]]
p_list[[4]]

patchwork::wrap_plots(p_list, ncol = 3, guides = "collect") -> p
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/annotation-allen-ontology_manual.pdf")

rm(list=ls())

library(Seurat)
library(tidyverse)
load("~/Downloads/Brain_region_annotation/df_region_list.rda")
df_region_list %>% 
  purrr::map(function(df_){
    df_
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) %>% 
  dplyr::select(cell, upper_onto) -> df_regions_all
dim(df_regions_all)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
xenium_obj_filtered@meta.data %>% 
  rownames_to_column("cell") %>% 
  {
    df_ = .
    left_join(x = df_, y = df_regions_all, by = "cell") -> df_
    df_ %>% column_to_rownames("cell") -> df_
    df_ %>% dplyr::filter(!is.na(upper_onto)) -> df_
    df_
  } -> meta.data

cell_to_keep = meta.data %>% rownames()
xenium_obj_filtered = subset(xenium_obj_filtered, cells = cell_to_keep)
gc()

all(colnames(xenium_obj_filtered) == rownames(meta.data))
xenium_obj_filtered@meta.data <- meta.data

Idents(xenium_obj_filtered) = xenium_obj_filtered$upper_onto
FindAllMarkers(xenium_obj_filtered, only.pos = TRUE) -> markers_regions

levels(markers_regions$cluster) %>% 
  .[c(1,4,5,3,2)] %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(cluster_){
    markers_regions %>% dplyr::filter(cluster == cluster_) -> df_
    head(df_, 20)
  })

#genes_selected = c("Slc17a7", "Ppp1r1b", "Tcf7l2", "Hap1", "Mog")
genes_selected = c("Slc17a7", "Ppp1r1b", "Tcf7l2", "Irs4", "Mog")

unique(xenium_obj_filtered@meta.data$upper_onto) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(upper_onto_){
    xenium_obj_filtered@meta.data %>% dplyr::filter(upper_onto == upper_onto_) -> meta_data
    data.frame(cell_type = upper_onto_, nCount = meta_data$nCount_Xenium)
  }) %>% do.call(rbind, .) %>% 
  "rownames<-"(., NULL) %>% 
  dplyr::mutate(cell_type = factor(cell_type, levels = rev(c("CTX", "CNU", "TH", "HY", "Fiber")))) %>% 
  ggplot(aes(y = cell_type, x = nCount, fill = cell_type)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) + 
  theme_classic() + 
  theme(legend.position = "none") -> p_nCount
p_nCount + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_nCount_1

unique(xenium_obj_filtered@meta.data$upper_onto) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(upper_onto_){
    xenium_obj_filtered@meta.data %>% dplyr::filter(upper_onto == upper_onto_) -> meta_data
    data.frame(cell_type = upper_onto_, nFeature = meta_data$nFeature_Xenium)
  }) %>% do.call(rbind, .) %>% 
  "rownames<-"(., NULL) %>% 
  dplyr::mutate(cell_type = factor(cell_type, levels = rev(c("CTX", "CNU", "TH", "HY", "Fiber")))) %>% 
  ggplot(aes(y = cell_type, x = nFeature, fill = cell_type)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) + 
  theme_classic() + 
  theme(legend.position = "none") -> p_nFeature
p_nFeature + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_nFeature_1

xenium_obj_filtered@meta.data -> meta_data
mat_selected = GetAssayData(object = xenium_obj_filtered, assay = "Xenium", layer = "counts")[genes_selected, ]
dim(mat_selected)

unique(meta_data$upper_onto) %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(upper_onto_){
    meta_data %>% dplyr::filter(upper_onto == upper_onto_) -> df_
    df_$nCount_Xenium %>% sum()
  }) %>% unlist() %>% median() -> median_

unique(meta_data$upper_onto) %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(upper_onto_){
    meta_data %>% dplyr::filter(upper_onto == upper_onto_) -> df_
    df_$nCount_Xenium %>% sum() -> lib_size
    cells_selected = rownames(df_)
    mat_selected[, cells_selected] %>% rowSums() %>% {(./lib_size)*median_} %>% 
      as.data.frame() %>% "colnames<-"(., "expr") %>% rownames_to_column("gene") %>% 
      dplyr::mutate(cell_type = upper_onto_, .before = 1) -> df_1
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_expr

unique(meta_data$upper_onto) %>% 
  "names<-"(.,.) %>% 
  #    .[1] %>% 
  purrr::map(function(upper_onto_){
    meta_data %>% dplyr::filter(upper_onto == upper_onto_) -> df_
    cells_selected = rownames(df_)
    n_cells = nrow(df_)
    mat_ = mat_selected[, cells_selected] 
    (mat_ != 0) %>% rowSums() %>% {./n_cells} %>% 
      as.data.frame() %>% "colnames<-"(., "non_zero_prop") %>% 
      rownames_to_column("gene") %>% 
      dplyr::mutate(cell_type = upper_onto_, .before = 1) -> df_1
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_non_zero_prop

full_join(x = df_expr, y = df_non_zero_prop, by = c("cell_type", "gene")) -> df_final
unique(df_final$gene) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    df_final %>% dplyr::filter(gene == gene_) %>% 
      dplyr::mutate(scaled_expr = scales::rescale(expr, to = c(0,1))) %>% 
      dplyr::mutate(non_zero_prop = non_zero_prop*100)
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_plot

df_plot$cell_type = factor(df_plot$cell_type, levels = rev(c("CTX", "CNU", "TH", "HY", "Fiber")))
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
ggsave(plot = p_dotplot_2, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/brain_region_marker_genes_dotplot.pdf", width = 12.1, height = 4.68) #Figure S3B

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    unique(xenium_obj_filtered@meta.data$upper_onto) %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(upper_onto_){
        xenium_obj_filtered@meta.data %>% dplyr::filter(upper_onto == upper_onto_, orig.ident == sample_) -> meta_data
        data.frame(cell_type = upper_onto_, nCount = meta_data$nCount_Xenium)
      }) %>% do.call(rbind, .) %>% 
      "rownames<-"(., NULL) %>% 
      dplyr::mutate(cell_type = factor(cell_type, levels = rev(c("CTX", "CNU", "TH", "HY", "Fiber")))) %>% 
      ggplot(aes(y = cell_type, x = nCount, fill = cell_type)) + 
      geom_boxplot(outlier.shape = NA, width = 0.3) + 
      theme_classic() + 
      theme(legend.position = "none") -> p_nCount
  }) -> p_nCount_list

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    unique(xenium_obj_filtered@meta.data$upper_onto) %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(upper_onto_){
        xenium_obj_filtered@meta.data %>% dplyr::filter(upper_onto == upper_onto_, orig.ident == sample_) -> meta_data
        data.frame(cell_type = upper_onto_, nFeature = meta_data$nFeature_Xenium)
      }) %>% do.call(rbind, .) %>% 
      "rownames<-"(., NULL) %>% 
      dplyr::mutate(cell_type = factor(cell_type, levels = rev(c("CTX", "CNU", "TH", "HY", "Fiber")))) %>% 
      ggplot(aes(y = cell_type, x = nFeature, fill = cell_type)) + 
      geom_boxplot(outlier.shape = NA, width = 0.3) + 
      theme_classic() + 
      theme(legend.position = "none") -> p_nFeature
    p_nFeature + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_nFeature_1    
  }) -> p_nFeature_list_1

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    unique(meta_data$upper_onto) %>% 
      "names<-"(.,.) %>% 
      #  .[1] %>% 
      purrr::map(function(upper_onto_){
        meta_data %>% dplyr::filter(upper_onto == upper_onto_, orig.ident == sample_) -> df_
        df_$nCount_Xenium %>% sum()
      }) %>% unlist() %>% median() -> median_ 
  }) -> median_list

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    unique(meta_data$upper_onto) %>% 
      "names<-"(.,.) %>% 
      #  .[1] %>% 
      purrr::map(function(upper_onto_){
        meta_data %>% dplyr::filter(upper_onto == upper_onto_, orig.ident == sample_) -> df_
        df_$nCount_Xenium %>% sum() -> lib_size
        cells_selected = rownames(df_)
        mat_selected[, cells_selected] %>% rowSums() %>% {(./lib_size)*median_} %>% 
          as.data.frame() %>% "colnames<-"(., "expr") %>% rownames_to_column("gene") %>% 
          dplyr::mutate(cell_type = upper_onto_, .before = 1) -> df_1
      }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_expr    
  }) -> df_expr_list

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    unique(meta_data$upper_onto) %>% 
      "names<-"(.,.) %>% 
      #    .[1] %>% 
      purrr::map(function(upper_onto_){
        meta_data %>% dplyr::filter(upper_onto == upper_onto_, orig.ident == sample_) -> df_
        cells_selected = rownames(df_)
        n_cells = nrow(df_)
        mat_ = mat_selected[, cells_selected] 
        (mat_ != 0) %>% rowSums() %>% {./n_cells} %>% 
          as.data.frame() %>% "colnames<-"(., "non_zero_prop") %>% 
          rownames_to_column("gene") %>% 
          dplyr::mutate(cell_type = upper_onto_, .before = 1) -> df_1
      }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_non_zero_prop    
  }) -> df_non_zero_prop_list

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    df_expr = df_expr_list[[sample_]]
    df_non_zero_prop = df_non_zero_prop_list[[sample_]]
    full_join(x = df_expr, y = df_non_zero_prop, by = c("cell_type", "gene")) -> df_final
  }) -> df_final_list

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    df_final = df_final_list[[sample_]]
    unique(df_final$gene) %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(gene_){
        df_final %>% dplyr::filter(gene == gene_) %>% 
          dplyr::mutate(scaled_expr = scales::rescale(expr, to = c(0,1))) %>% 
          dplyr::mutate(non_zero_prop = non_zero_prop*100)
      }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_plot    
  }) -> df_plot_list

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  purrr::map(function(sample_){
    df_plot = df_plot_list[[sample_]]
    df_plot$cell_type = factor(df_plot$cell_type, levels = rev(c("CTX", "CNU", "TH", "HY", "Fiber")))
    df_plot$gene = factor(df_plot$gene, levels = genes_selected)
    df_plot %>% 
      ggplot(aes(x = gene, y = cell_type)) +
      geom_point(aes(size = non_zero_prop, color = scaled_expr)) + 
      scale_size(range = c(0,15), limits = c(0,100)) + 
      scale_color_gradient(low = "grey", high = "blue") + 
      theme_classic() -> p_dotplot
    p_dotplot + ylab(NULL) + theme(axis.text.y = element_blank()) -> p_dotplot_1
    
    p_nCount = p_nCount_list[[sample_]]
    p_nFeature_1 = p_nFeature_list_1[[sample_]]
    
    patchwork::wrap_plots(p_nCount, p_nFeature_1, p_dotplot_1, nrow = 1) + 
      patchwork::plot_layout(widths = c(1,1,3)) -> p_dotplot_2
    p_dotplot_2 + patchwork::plot_annotation(title = sample_)
  }) -> p_dotplot_2_list

pdf(file = "~/Dropbox/Spatial_transcriptomics_new/github/plots/brain_region_marker_genes_dotplot_per_sample.pdf", width = 13.76, height = 6.32)
p_dotplot_2_list
dev.off()

table(xenium_obj_filtered@meta.data$orig.ident, xenium_obj_filtered@meta.data$upper_onto)