library(tidyverse)
library(Seurat)

list.dirs("~/Dropbox/Spatial_transcriptomics_new/data_proseg_processed", recursive = F, full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)", "\\1", .)) %>% 
  gtools::mixedsort() %>% rev() -> dat_dir_

meta_data_list = list()
dat_dir_ %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(file_, group_){
    print(group_)
    fov_ = gsub("-", ".", group_)
    cell_df_path = sprintf("%s/cells.csv.gz", file_)
    output_sc = sprintf("%s/%s-sc.rda", "~/Dropbox/Spatial_transcriptomics_new/sc/segmented", group_)
    
    xenium.obj <- LoadXenium(file_, fov = fov_, mols.qv.threshold = 20, segmentations = "cell") 
    read_csv(cell_df_path, col_names = T) -> cell_df
    cell_df %>% dplyr::select(cell_id, x_centroid, y_centroid, cell_area) -> cell_df
    xenium.obj@meta.data %>% 
      rownames_to_column("cell_id") %>% 
      left_join(x = ., y = cell_df, by = "cell_id") %>% 
      column_to_rownames("cell_id") -> meta_data
    xenium.obj@meta.data <- meta_data
    meta_data_list[[group_]] <<- meta_data
    save(xenium.obj, file = output_sc)
    gc()
  })

save(meta_data_list, file = "~/Dropbox/Spatial_transcriptomics_new/sc/segmented/meta_data_list.rda")
rm(list=ls())

library(Seurat)
library(tidyverse)
list.files("~/Dropbox/Spatial_transcriptomics_new/sc/segmented", pattern = "CT.+sc\\.rda", full.names = T) %>% 
  gtools::mixedsort() %>% rev() %>% 
  "names<-"(., gsub("/.+/(CT-\\d+-\\d+?)-sc\\.rda", "\\1", .)) -> sc_files

sc_files %>% 
#  .[1:2] %>% 
  map(function(file_){
    print(file_)
    env <- new.env()
    load(file_, envir = env)
    return(env$xenium.obj)
  }) -> xenium.obj.list
names(xenium.obj.list) -> cell_ids_

options(future.globals.maxSize = Inf)
xen_merged <- merge(
  x = xenium.obj.list[[1]],
  y = xenium.obj.list[-1],
  add.cell.ids = cell_ids_,
  project = "MergedXenium",
  merge.data = TRUE
)
colnames(xen_merged) %>% gsub("(CT-\\d+-\\d+)_.+", "\\1", .) %>% gsub("-", ".", .) -> xen_merged$fov
xen_merged@meta.data %>% head()

rm(xenium.obj.list);gc()
xen_merged[["Xenium"]] <- JoinLayers(xen_merged[["Xenium"]])
xen_merged@assays$Xenium@layers$counts %>% dim()

save(xen_merged, file = "~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc.rda")

library(tidyverse)
library(Seurat)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc.rda")

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(CT_){
    print(CT_)
    CT_ = gsub("-", ".", CT_)
    xy_coordinates = GetTissueCoordinates(xen_merged, image = CT_)
  }) -> xy_coordinates
save(xy_coordinates, file = "~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")
