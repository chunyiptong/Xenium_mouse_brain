#ml GCC/13.2.0 OpenMPI/4.1.6 R/4.4.1
.libPaths("/scratch/user/tongtong1993/project/R/4.4")
library(Seurat)
library(tidyverse)
library(future)
library(glmGamPoi)

load("~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/wmb_merged_ref_sc.rda")
gc()
allen.wmb.ref = seurat_obj
rm(seurat_obj);gc()

DefaultAssay(allen.wmb.ref) <- "SCT"
allen.wmb.ref[["RNA"]] <- NULL  # drop RNA assay entirely                                                                                                                         
gc()

load("~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-mini.rda")
xenium_obj_filtered[["Xenium"]] = NULL

feat <- intersect(                                                                       
		  rownames(GetAssayData(allen.wmb.ref, "SCT", "scale.data")),
		  rownames(GetAssayData(xenium_obj_filtered,  "SCT", "scale.data"))
		  )
gc()

options(future.globals.maxSize = Inf)
anchors <- FindTransferAnchors(
			       reference = allen.wmb.ref, 
				 query = xenium_obj_filtered,
				 normalization.method = "SCT", reduction = "rpca",
				 reference.reduction  = "pca", 
				 dims = 1:30, 
				 features = feat,
				 recompute.residuals  = FALSE                                                             
			       )
gc()
save(anchors, file = "~/Desktop/Research/Analysis/spatial_transcriptome/External_data/Allen_Brain_Cell_Atlas/Allen_wmb_to_xenium_anchors.rda")

allen.wmb.ref@meta.data -> meta.data
rm(allen.wmb.ref);gc()
xenium_obj_filtered <- TransferData(
		       anchorset        = anchors,
		         refdata          = meta.data$custom_class,
		         prediction.assay = FALSE,
			   store.weights    = FALSE,
			   weight.reduction = "pcaproject",
			     dims             = 1:30,
			     query            = xenium_obj_filtered
		     )
gc()
rpca_cell_annotation = xenium_obj_filtered@meta.data
save(rpca_cell_annotation, file = "~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/rpca_cell_annotation_wmb.rda")
rm(list=ls())

library(tidyverse)
library(Seurat)

load("~/Downloads/rpca_cell_annotation_wmb.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-mini.rda")

unique(rpca_cell_annotation$predicted.id) %>% 
  purrr::map(function(predicted.id_){
    rpca_cell_annotation %>% dplyr::filter(predicted.id == predicted.id_) -> df_
    df_ %>% 
      ggplot(aes(x = predicted.id.score)) + 
      geom_histogram(bins = 30) + 
      theme_classic() + 
      ggtitle(sprintf("%s n=%s", predicted.id_, nrow(df_)))
  }) -> p_list_score
length(p_list_score)
patchwork::wrap_plots(p_list_score, ncol = 5) -> p
ggsave(plot = p, filename = "~/Downloads/rpca_cell_annotation_wmb_score_histogram.pdf", width = 15.28, height = 2.47) #Figure S2D

(rpca_cell_annotation %>% dplyr::filter(predicted.id.score >= 0.8) %>% nrow())/
(rpca_cell_annotation %>% nrow())

cells_to_keep = rpca_cell_annotation %>% dplyr::filter(predicted.id.score >= 0.8) %>% rownames()
subset(xenium_obj_filtered, cells = cells_to_keep) -> xenium_obj_filtered

DefaultAssay(xenium_obj_filtered) <- "SCT"
xenium_obj_filtered <- RunPCA(xenium_obj_filtered, npcs = 30, verbose = T)
pca_mat <- Embeddings(xenium_obj_filtered, "pca")[, 1:30]
meta.data = xenium_obj_filtered@meta.data
save(pca_mat, meta.data, file="~/Downloads/merged-sc-filtered-normalized-matched-nucleus-pca-mini.rda")
rm(list=ls())
gc()

library(tidyverse)
library(uwot)
library(RSpectra)

load("~/Downloads/merged-sc-filtered-normalized-matched-nucleus-pca-mini.rda")
load("~/Downloads/rpca_cell_annotation_wmb.rda")

rownames(pca_mat) -> cells_to_keep
rpca_cell_annotation[cells_to_keep, ] -> rpca_cell_annotation

umap_embedding <- uwot::umap(
  X = pca_mat,
  n_neighbors  = 30,
  n_components = 2,
  metric       = "cosine",
  n_threads    = 6,
  verbose      = TRUE,
  init = "spectral"
)
umap_embedding %>% 
  as.data.frame() %>% 
  "colnames<-"(., c("UMAP_1", "UMAP_2")) -> umap_embedding_df
cbind(umap_embedding_df, rpca_cell_annotation[,c("predicted.id", "predicted.id.score")]) %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = predicted.id)) + 
  ggrastr::geom_point_rast(size = 0.1) + 
  theme_classic() -> p_umap

db = dbscan::dbscan((umap_embedding_df %>% "colnames<-"(., c("x", "y"))), eps = 0.8, minPts = 20)

umap_embedding_df %>% dplyr::mutate(cluster = as.character(db$cluster)) %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  ggrastr::geom_point_rast(size = 0.1) + 
  theme_classic() -> p_umap_cluster

cbind(umap_embedding_df, rpca_cell_annotation[,c("predicted.id", "predicted.id.score")]) %>% 
  dplyr::mutate(cluster = as.character(db$cluster)) %>% 
  {
    df_ = .
#    table(df_$cluster, df_$predicted.id)
    df_$cluster %>% unique() %>% 
      purrr::map(function(cluster_){
        df_ %>% dplyr::filter(cluster == cluster_) -> df_1
        major_cell_type = df_1$predicted.id %>% table() %>% sort(decreasing = T) %>% names() %>% .[1]
        df_1 %>% dplyr::filter(predicted.id == major_cell_type) -> df_1
        if (nrow(df_1) >= 100){
          return(df_1)
        }else{
          return(NULL)
        }
      }) %>% do.call(rbind, .)
  } -> umap_embedding_df_cleaned
table(umap_embedding_df_cleaned$cluster)
table(umap_embedding_df_cleaned$predicted.id)

umap_embedding_df_cleaned %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = predicted.id)) + 
  ggrastr::geom_point_rast(size = 0.1) + 
  theme_classic() -> p_umap_cleaned


umap_embedding_2 <- uwot::umap2(
  X = pca_mat,
  n_neighbors  = 30,
  n_components = 2,
  metric       = "cosine",
  n_threads    = 6,
  verbose      = TRUE,
  init = "pca"
)
umap_embedding_2 %>% 
  as.data.frame() %>% 
  "colnames<-"(., c("UMAP_1", "UMAP_2")) -> umap_embedding_2_df

cbind(umap_embedding_2_df, rpca_cell_annotation[,c("predicted.id", "predicted.id.score")]) %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = predicted.id)) + 
  ggrastr::geom_point_rast(size = 0.1) + 
  theme_classic() -> p_umap_2

rm(db, meta.data, p_umap_2, p_umap_cluster, rpca_cell_annotation, umap_embedding_2, umap_embedding_2_df, cells_to_keep);gc()

library(Seurat)
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-mini.rda")
cells_to_keep = rownames(pca_mat)
subset(xenium_obj_filtered, cells = cells_to_keep) -> xenium_obj_filtered
  
DefaultAssay(xenium_obj_filtered)

xenium_obj_filtered[["pca"]] <- CreateDimReducObject(
  embeddings = pca_mat,
  key        = "PC_",
  assay      = DefaultAssay(xenium_obj_filtered)
)

xenium_obj_filtered[["umap"]] <- CreateDimReducObject(
  embeddings = umap_embedding, 
  key = "UMAP_",
  assay = DefaultAssay(xenium_obj_filtered)
)

load("~/Downloads/rpca_cell_annotation_wmb.rda")
left_join(
x = xenium_obj_filtered@meta.data %>% dplyr::select(-matches("res|clusters")) %>% rownames_to_column("cell_id"), 
y = rpca_cell_annotation %>% rownames_to_column("cell_id") %>% dplyr::select(cell_id, predicted.id, predicted.id.score),
by = "cell_id"
) %>% column_to_rownames("cell_id") -> meta_data

all(colnames(xenium_obj_filtered) == rownames(meta_data))
xenium_obj_filtered@meta.data <- meta_data
rm(meta_data)

DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "predicted.id", label = T) + ggtitle("Celltypes") -> p_umap_annotated

rownames(umap_embedding_df_cleaned) -> cells_to_keep
cells_to_keep %>% length()

xenium_obj_filtered = subset(xenium_obj_filtered, cells = cells_to_keep)
gc()
DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "predicted.id", label = T) + ggtitle("Celltypes") -> p_umap_annotated_cleaned
ggsave(plot = p_umap_annotated_cleaned, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/umap_wmb_label_transfer_cleaned_xenium.pdf") #Figure 1B

factor(umap_embedding_df_cleaned[colnames(xenium_obj_filtered), "cluster"], levels = as.character(seq(1,13))) -> cluster_
xenium_obj_filtered$cluster = cluster_

DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "cluster", label = T) -> p_umap_annotated_cluster

DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "orig.ident", label = T) + ggtitle("Samples") -> p_umap_annotated_sample
ggsave(plot = p_umap_annotated_sample, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/umap_wmb_label_transfer_cleaned_xenium_by_sample.pdf") #Figure S2B

table(xenium_obj_filtered$orig.ident)[c(4,5,6,1,2,3)]
save(xenium_obj_filtered, file = "~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")

###Further delete minor clusters
library(Seurat)
library(tidyverse)
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")

xenium_obj_filtered@meta.data$cluster %>% levels()
xenium_obj_filtered@meta.data$cluster %>% table() %>% sort() %>% .[5:length(.)] %>% names() -> clusters_to_keep
xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% clusters_to_keep) %>% rownames() -> cells_to_keep

xenium_obj_filtered = subset(xenium_obj_filtered, cells = cells_to_keep)
cluster_ = xenium_obj_filtered@meta.data$cluster
as.character(seq_along(levels(cluster_))) -> levels(cluster_)
DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "cluster", label = T) + 
DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "predicted.id", label = T) + 
DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "orig.ident", label = T)

table(xenium_obj_filtered$orig.ident) %>% .[c(4,5,6,1,2,3)]
save(xenium_obj_filtered, file = "/Volumes/Crucial_X8/Spatial//merged-sc-filtered-normalized-matched-nucleus-annotated-mini-cleaned.rda")


#Check minor clusters ----
library(tidyverse)
library(Seurat)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")

DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "cluster", label = T) -> p_umap_annotated_cluster
table(xenium_obj_filtered$cluster, xenium_obj_filtered$predicted.id)

C11_immune = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("11")) %>% rownames()
C13_GABA = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("13")) %>% rownames()
C4_Gluta = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("4")) %>% rownames()
C9_GABA = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("9")) %>% rownames()
C8_Glial = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("8")) %>% rownames()
C3_Glial = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("3")) %>% rownames()
C5_Glial = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("5")) %>% rownames()
C7_GABA = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("7")) %>% rownames()
C2_GABA = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("2")) %>% rownames()
C1_Gluta = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("1")) %>% rownames()
C10_Gluta = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("10")) %>% rownames()
C6_Vascular = xenium_obj_filtered@meta.data %>% dplyr::filter(cluster %in% c("6")) %>% rownames()

#where are this minor cluster cells located? 

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

xy_coordinates %>% 
  purrr::map(function(df_){
    df_ %>% dplyr::filter(cell %in% colnames(xenium_obj_filtered)) -> df_
  }) -> xy_coordinates

list(
  C11_immune = C11_immune,
  C13_GABA = C13_GABA,
  C4_Gluta = C4_Gluta,
  C9_GABA = C9_GABA,
  C8_Glial = C8_Glial,
  C3_Glial = C3_Glial,
  C5_Glial = C5_Glial,
  C7_GABA = C7_GABA,
  C2_GABA = C2_GABA, 
  C1_Gluta = C1_Gluta,
  C10_Gluta = C10_Gluta,
  C6_Vascular = C6_Vascular
) %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(cells, cluster_){
    names(xy_coordinates) %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      purrr::map(function(sample_){
        xy_coordinates[[sample_]] -> df_
        df_ %>% dplyr::mutate(group = case_when(
          cell %in% cells ~ cluster_, 
          TRUE ~ "Other"
        )) -> df_
        df_ %>% dplyr::arrange(desc(group)) -> df_
        df_ %>% 
          ggplot(aes(x = x, y = y, color = group)) + 
          geom_point(size = 0.5, alpha = 0.5) + 
          ggtitle(sprintf("%s: %s", cluster_, sample_)) -> p
      }) -> p_list
    patchwork::wrap_plots(p_list, ncol = 3, guides = "collect")
  }) -> p_list

Idents(xenium_obj_filtered) <- xenium_obj_filtered$cluster
markers = FindAllMarkers(xenium_obj_filtered, only.pos = T)

levels(markers$cluster) %>% 
  "names<-"(., sprintf("C_%s", .)) %>% 
  purrr::map(function(cluster_){
    markers %>% dplyr::filter(cluster == cluster_) -> df_
  }) -> markers_list
save(markers_list, file = "~/Downloads/cluster_markers_list.rda")

markers_list %>% 
  purrr::map(function(df_){
    head(df_, 20)
  })

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
rm(class_)

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
  )) -> wmb_cell_annotation

table(wmb_cell_annotation$custom_class)

wmb_cell_annotation %>% 
  column_to_rownames("cell_id") %>% 
  .[C11_immune, ] %>% 
  .$custom_class %>% table() %>% 
  {(./sum(.))*100}

wmb_cell_annotation %>% 
  column_to_rownames("cell_id") %>% 
  .[C11_immune, ] %>% 
  {
    df_ = .
    df_$custom_class %>% unique() %>% sort() -> levels_
    df_$custom_class = factor(df_$custom_class, levels = levels_)
    
    df_ %>% 
      ggplot(aes(x = custom_class, y = class_bootstrapping_probability, color = custom_class)) + 
      geom_jitter(size = 0.5) + 
      theme_classic() + 
      theme(legend.position = "none") + 
      ggtitle("C11 Immune")
  } -> p_C11_immune

wmb_cell_annotation %>% 
  column_to_rownames("cell_id") %>% 
  .[C13_GABA, ] %>% 
  .$custom_class %>% table() %>% 
  {(./sum(.))*100}

wmb_cell_annotation %>% 
  column_to_rownames("cell_id") %>% 
  .[C13_GABA, ] %>% 
  {
    df_ = .
    df_$custom_class %>% unique() %>% sort() -> levels_
    df_$custom_class = factor(df_$custom_class, levels = levels_)
    
    df_ %>% 
      ggplot(aes(x = custom_class, y = class_bootstrapping_probability, color = custom_class)) + 
      geom_jitter(size = 0.5) + 
      theme_classic() + 
      theme(legend.position = "none") + 
      ggtitle("C13 GABAergic")
  } -> p_C13_GABA

wmb_cell_annotation %>% 
  column_to_rownames("cell_id") %>% 
  .[C4_Gluta, ] %>% 
  .$custom_class %>% table() %>% 
  {(./sum(.))*100}

wmb_cell_annotation %>% 
  column_to_rownames("cell_id") %>% 
  .[C4_Gluta, ] %>% 
  {
    df_ = .
    df_$custom_class %>% unique() %>% sort() -> levels_
    df_$custom_class = factor(df_$custom_class, levels = levels_)
    
    df_ %>% 
      ggplot(aes(x = custom_class, y = class_bootstrapping_probability, color = custom_class)) + 
      geom_jitter(size = 0.5) + 
      theme_classic() + 
      theme(legend.position = "none") + 
      ggtitle("C4 Glutamatergic")
  } -> p_C4_Gluta

patchwork::wrap_plots(p_C4_Gluta, p_C11_immune, p_C13_GABA, ncol = 1) -> p1
patchwork::wrap_plots(p_umap_annotated_cluster, p1, ncol = 2) -> p2
ggsave(plot = p2, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/umap_minor_clusters.pdf")
rm(p_C11_immune, p_C13_GABA, p_C4_Gluta, p_umap_annotated_cluster, p1, p2, C11_immune, C13_GABA, C4_Gluta);gc()

table(xenium_obj_filtered$cluster)

cbind(
  xenium_obj_filtered@meta.data["predicted.id"],
  wmb_cell_annotation %>% column_to_rownames("cell_id") %>% .[colnames(xenium_obj_filtered), ] %>% .["custom_class"]
) %>% 
  {
    df_ = .
    table(df_$predicted.id, df_$custom_class) %>% 
      as.matrix.data.frame() %>% .[,-5] %>% 
      apply(., 1, function(x){max(x)}) %>% 
      sum()
  }
# If rpca wmb >= 0.8 and mapmycells (same celltype), 525392 cells retained

cbind(
  xenium_obj_filtered@meta.data["predicted.id"],
  wmb_cell_annotation %>% column_to_rownames("cell_id") %>% .[colnames(xenium_obj_filtered), ] %>% .["custom_class"]
) %>% dplyr::filter(predicted.id == custom_class) %>% 
  rownames() -> cells_to_keep

xenium_obj_filtered = subset(xenium_obj_filtered, cells = cells_to_keep)
gc()

DefaultAssay(xenium_obj_filtered) <- "SCT"
xenium_obj_filtered <- RunPCA(xenium_obj_filtered, npcs = 30, verbose = T)
xenium_obj_filtered <- RunUMAP(xenium_obj_filtered, dims = 1:30, verbose = T)

DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "predicted.id", label = T) -> p_umap_celltype
ggsave(plot = p_umap_celltype, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/umap_wmb_label_transfer_cleaned_xenium_1.pdf")
DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "orig.ident", label = T) -> p_umap_samples
ggsave(plot = p_umap_samples, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/umap_wmb_label_transfer_cleaned_xenium_by_sample_1.pdf")

save(xenium_obj_filtered, file="/Volumes/Crucial_X8/Spatial/merged-sc-filtered-normalized-matched-nucleus-annotated-1-mini.rda")
rm(p_umap_samples, p_umap_celltype);gc()
####

cbind(xenium_obj_filtered@meta.data, 
      wmb_cell_annotation %>% column_to_rownames("cell_id") %>% .[colnames(xenium_obj_filtered), ]) %>% 
  dplyr::filter(class_bootstrapping_probability >= 0.5) %>% 
  dim()
#484665 13

cbind(xenium_obj_filtered@meta.data, 
      wmb_cell_annotation %>% column_to_rownames("cell_id") %>% .[colnames(xenium_obj_filtered), ]) %>% 
  dplyr::filter(class_bootstrapping_probability >= 0.5) %>% 
  rownames() -> cells_to_keep

subset(xenium_obj_filtered, cells = cells_to_keep) -> xenium_obj_filtered
DefaultAssay(xenium_obj_filtered) <- "SCT"
xenium_obj_filtered <- RunPCA(xenium_obj_filtered, npcs = 30, verbose = T)
xenium_obj_filtered <- RunUMAP(xenium_obj_filtered, dims = 1:30, verbose = T)

DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "predicted.id", label = T) -> p_umap_celltype
ggsave(plot = p_umap_celltype, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/umap_wmb_label_transfer_cleaned_xenium_2.pdf")
DimPlot(xenium_obj_filtered, reduction = "umap", group.by = "orig.ident", label = T) -> p_umap_samples
ggsave(plot = p_umap_samples, filename = "~/Dropbox/Spatial_transcriptomics_new/github/plots/umap_wmb_label_transfer_cleaned_xenium_by_sample_2.pdf")

save(xenium_obj_filtered, file="/Volumes/Crucial_X8/Spatial/merged-sc-filtered-normalized-matched-nucleus-annotated-2-mini.rda")
