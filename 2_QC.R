library(Seurat)
library(tidyverse)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc.rda")

xen_merged$fov %>% 
  unique() %>% 
  purrr::map(function(fov_){
    xen_merged@meta.data %>% dplyr::filter(fov == fov_) %>% 
      dplyr::select(nCount_Xenium, nFeature_Xenium, cell_area) -> df_
    df_ %>% dplyr::mutate(group = gsub("\\.", "-", fov_)) -> df_
  }) %>% do.call(rbind, .) -> df_

library(ggbreak)
df_$group = factor(df_$group, levels = c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3"))
df_ %>% 
  ggplot(aes(x = group, y = nCount_Xenium, fill = group)) + 
  geom_violin() + 
  geom_boxplot(outliers = F, width = 0.5) + 
  scale_y_break(breaks = c(7500,20000), scales = 0.5) + 
  theme_classic() -> p_nCount_pre_QC

df_ %>% 
  ggplot(aes(x = group, y = nFeature_Xenium, fill = group)) + 
  geom_violin() + 
  geom_boxplot(outliers = F, width = 0.5) + 
  theme_classic() -> p_nFeature_pre_QC

df_ %>% 
  ggplot(aes(x = group, y = cell_area, fill = group)) + 
  geom_violin() + 
  geom_boxplot(outliers = F, width = 0.5) + 
  ggbreak::scale_y_break(breaks = c(1000,3000), scales = 0.5) +
  theme_classic() -> p_cell_area_pre_QC


levels(df_$group) %>% 
  purrr::map(function(group_){
    df_ %>% dplyr::filter(group == group_) -> df_
    # filter nCount by median+3*medain_absolute_deviation
    lower_threshold_nCount = 100
    upper_threshold_nCount = median(df_$nCount_Xenium) + 3*mad(df_$nCount_Xenium)
    lower_threshold_nFeature = 50
    upper_threshold_nFeature = median(df_$nFeature_Xenium) + 3*mad(df_$nFeature_Xenium)
    upper_threshold_cell_area = median(df_$cell_area) + 3*mad(df_$cell_area)
    sprintf("For %s, nCount_Xenium threshold: %s - %s; nFeature_Xenium threshold: %s - %s, cell area threshold: %s", group_, lower_threshold_nCount, upper_threshold_nCount, lower_threshold_nFeature, upper_threshold_nFeature, upper_threshold_cell_area) %>% print()
    df_ %>% dplyr::filter(nCount_Xenium > lower_threshold_nCount & nCount_Xenium < upper_threshold_nCount) %>%
      dplyr::filter(nFeature_Xenium > lower_threshold_nFeature & nFeature_Xenium < upper_threshold_nFeature) -> df_
    df_ %>% dplyr::filter(cell_area < upper_threshold_cell_area) -> df_
  }) %>% do.call(rbind, .) -> df_1

df_1 %>% 
  ggplot(aes(x = group, y = nCount_Xenium, fill = group)) + 
  geom_violin() + 
  geom_boxplot(outliers = F, width = 0.5) + 
  theme_classic() -> p_nCount_after_QC

df_1 %>% 
  ggplot(aes(x = group, y = nFeature_Xenium, fill = group)) + 
  geom_violin() + 
  geom_boxplot(outliers = F, width = 0.5) + 
  theme_classic() -> p_nFeature_after_QC

df_1 %>% 
  ggplot(aes(x = group, y = cell_area, fill = group)) + 
  geom_violin() + 
  geom_boxplot(outliers = F, width = 0.5) + 
  theme_classic() -> p_cell_area_after_QC

patchwork::wrap_plots(p_nCount_after_QC, p_nFeature_after_QC, p_cell_area_after_QC, ncol = 3, guides = "collect") #Figure S1

rbind(
df_$group %>% table() %>% as.list() %>% as.data.frame(),
df_1$group %>% table() %>% as.list() %>% as.data.frame()
) %>% "rownames<-"(., c("pre-QC", "after-QC")) %>% 
  rowwise() %>% 
  dplyr::mutate(total = sum(c_across(everything()))) %>% 
  ungroup() %>% 
  dplyr::mutate(across(everything(), function(x){
    max_ = max(x)
    percentage_ = (x/max_)*100
    return(sprintf("%s (%.1f%%)", x, percentage_))
  })) %>% dplyr::mutate(group = c("Pre-QC", "After-QC"), .before = 1) -> table_QC 
  
rownames(df_1) -> cells_select
xen_merged_filtered <- subset(xen_merged, cells = cells_select)
save(xen_merged_filtered, file = "~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered.rda")

