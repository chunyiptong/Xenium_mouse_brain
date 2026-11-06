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
ggsave(plot = p, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/annotation-allen-ontology_manual.pdf") #Figure 1C & Figure S3A