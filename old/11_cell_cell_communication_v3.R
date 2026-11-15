library(tidyverse)
library(Seurat)
library(SpatialCellChat)
library(CellChat)
library(ComplexHeatmap)
library(RColorBrewer)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates-aligned.rda")

list.files(path = "~/Dropbox/Spatial_transcriptomics_new/rda/Brain_region_annotation/SCN_SON_PVN", pattern = "v5\\.csv", full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)_regions.+", "\\1", .) %>% gsub("_", "-", .)) %>% 
  #  .[1] %>% 
  purrr::map(function(file_){
    read_csv(file_, skip = 2) -> df_
    df_$Selection %>% unique() %>% {.[grepl("^SCN_\\d+$", .)]} -> Selection_
    df_ = df_ %>% dplyr::filter(Selection %in% Selection_)
  }) -> SCN_polygon

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
        df_1 %>% dplyr::select(-x, -y) -> df_1
      }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_final
    print(sprintf("%s, n_cells = %s", group_, length(df_final$cell)))
    return(df_final)
  }) %>% do.call(rbind, .) -> SCN_cells
SCN_cells %>% dplyr::filter(cell %in% colnames(xenium_obj_filtered)) -> SCN_cells

xenium_obj_filtered = subset(xenium_obj_filtered, cells = SCN_cells$cell)
gc()

xenium_obj_filtered@meta.data %>% dplyr::filter(predicted.id == "GABAergic") %>% rownames() -> GABA_cells
xenium_obj_filtered_GABA = subset(xenium_obj_filtered, cells = GABA_cells)
dim(xenium_obj_filtered_GABA)
xenium_obj_filtered_GABA@meta.data %>% head()

xenium_obj_filtered_GABA_split = SplitObject(xenium_obj_filtered_GABA, split.by = "orig.ident")
for (i in 1:length(xenium_obj_filtered_GABA_split)) {
  xenium_obj_filtered_GABA_split[[i]] <- SCTransform(xenium_obj_filtered_GABA_split[[i]], 
                                             vars.to.regress = c("nCount_Xenium", "nFeature_Xenium"), assay="Xenium", verbose = TRUE, return.only.var.genes = FALSE)
}

v.common <- rownames(xenium_obj_filtered_GABA_split[[1]]$SCT@data)
for (  i in 1:length(xenium_obj_filtered_GABA_split)){
  if (i == 1){
    print(sprintf("condtion %s number of gene in SCT matrix: %s", i, length(v.common))) 
  }else{
    v.common <- intersect(v.common, rownames(xenium_obj_filtered_GABA_split[[i]]$SCT@data))
    print(sprintf("condtion %s number of gene in SCT matrix: %s", i, length(rownames(xenium_obj_filtered_GABA_split[[i]]$SCT@data))))
  }
} 
length( v.common );
cat( "length(v.common)=", length(v.common), "\n" )

v.inte.features <- vector()
if ( 1 ){                                     # find variable features/genes common to all
  v.var.common <- xenium_obj_filtered_GABA_split[[1]]@assays$SCT@var.features
  print(sprintf("condition %s number of SCT var.features: %s", 1, length(v.var.common)))
  for (i in 2:length(xenium_obj_filtered_GABA_split)){
    v.var.common <- c(v.var.common, xenium_obj_filtered_GABA_split[[i]]@assays$SCT@var.features)
    print(sprintf("condition %s number of SCT var.features: %s", i, length(xenium_obj_filtered_GABA_split[[i]]@assays$SCT@var.features)))
  }
  length(unique(v.var.common))
  v.var.common <- names(which(table(v.var.common) == length(xenium_obj_filtered_GABA_split)))
  
  cat("length(v.var.common)=", length(v.var.common), "\n")
  length(v.var.common)
  v.inte.features <- v.var.common
} else {
  v.inte.features <- SelectIntegrationFeatures(object.list = xenium_obj_filtered_GABA_split, nfeatures = 3000)
}

nDims = 50
k.anchor.use = 5
k.filter.use = 199
k.score.use = 30
#k.weight.use = 50
k.weight.use = 80
xenium_obj_filtered_GABA_split.anchors <- FindIntegrationAnchors(object.list = xenium_obj_filtered_GABA_split, dims = 1:nDims, 
                                                         assay=rep( "SCT", length(xenium_obj_filtered_GABA_split)), 
                                                         anchor.features=v.inte.features, k.anchor=k.anchor.use, 
                                                         k.filter=k.filter.use, k.score=k.score.use, verbose=TRUE )
xenium_obj_filter_GABA_integrated = IntegrateData(anchorset = xenium_obj_filtered_GABA_split.anchors, dims = 1:nDims, features.to.integrate=v.common, 
                                   k.weight = k.weight.use, verbose = T, new.assay.name = "XeniumIntegrated")

VariableFeatures(xenium_obj_filter_GABA_integrated, assay="XeniumIntegrated") <- v.inte.features
xenium_obj_filter_GABA_integrated@assays$XeniumIntegrated@key <- "xeniumintegrated_"

DefaultAssay(xenium_obj_filter_GABA_integrated) <- "XeniumIntegrated"
xenium_obj_filter_GABA_integrated <- ScaleData(xenium_obj_filter_GABA_integrated, verbose=TRUE, assay="XeniumIntegrated", features=rownames(xenium_obj_filter_GABA_integrated$XeniumIntegrated@data))


xenium_obj_filter_GABA_integrated <- RunPCA(xenium_obj_filter_GABA_integrated, assay="XeniumIntegrated", npcs = nDims, verbose = FALSE, 
                             reduction.name = "xeniumintegrated_pca", 
                             reduction.key = "xeniumintegratedpc_")
xenium_obj_filter_GABA_integrated <- RunUMAP(xenium_obj_filter_GABA_integrated, reduction = "xeniumintegrated_pca", dims = 1:nDims, 
                              reduction.name = "xeniumintegrated_umap", 
                              reduction.key = "xeniumintegratedumap_", 
                              assay = "XeniumIntegrated")
xenium_obj_filter_GABA_integrated <- FindNeighbors(xenium_obj_filter_GABA_integrated, assay="XeniumIntegrated", reduction = "xeniumintegrated_pca", dims = 1:nDims, force.recalc=T, 
                                    graph.name = c("xeniumintegrated_nn", "xeniumintegrated_snn"))
xenium_obj_filter_GABA_integrated <- FindClusters(xenium_obj_filter_GABA_integrated, assay="XeniumIntegrated", resolution = 0.3, graph.name = "xeniumintegrated_snn")
xenium_obj_filter_GABA_integrated@meta.data %>% head()
xenium_obj_filter_GABA_integrated$xeniumintegrated_snn_res.0.3 -> xenium_obj_filter_GABA_integrated$seurat_clusters

DimPlot(xenium_obj_filter_GABA_integrated, group.by = "seurat_clusters", label = TRUE, reduction = "xeniumintegrated_umap", repel = TRUE)
DimPlot(xenium_obj_filter_GABA_integrated, group.by = "orig.ident", label = TRUE, reduction = "xeniumintegrated_umap", repel = TRUE)

Idents(xenium_obj_filter_GABA_integrated) <- xenium_obj_filter_GABA_integrated$seurat_clusters

all(colnames(xenium_obj_filtered_GABA) == colnames(xenium_obj_filter_GABA_integrated))
xenium_obj_filtered_GABA$integrated_cluster = xenium_obj_filter_GABA_integrated$xeniumintegrated_snn_res.0.3

markers_SCN_neurons = FindAllMarkers(xenium_obj_filtered_GABA, only.pos = FALSE, group.by = "integrated_cluster")
levels(markers_SCN_neurons$cluster) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(cluster_){
    markers_SCN_neurons %>% dplyr::filter(cluster == cluster_) %>% 
      dplyr::arrange(p_val_adj) %>% 
      dplyr::filter(avg_log2FC > 0, p_val_adj < 0.05) -> df_
  }) -> markers_SCN_neurons_list

markers_SCN_neurons_list %>% purrr::map(function(df_){
  df_$gene
})

markers_SCN_neurons_list$`0` %>% dplyr::filter(gene %in% c("Vip", "Grp", "Nms", "Avp", "C1ql3", "Bdnf", "Cck", "lrs4"))
markers_SCN_neurons_list$`1` %>% dplyr::filter(gene %in% c("Vip", "Grp", "Nms", "Avp", "C1ql3", "Bdnf", "Cck", "lrs4"))

xenium_obj_filter_GABA_integrated@meta.data %>% dplyr::filter(!(xeniumintegrated_snn_res.0.3 %in% c("2", "3"))) %>% rownames() -> cells_
xenium_obj_filter_GABA_integrated = subset(xenium_obj_filter_GABA_integrated, cell = cells_)
dim(xenium_obj_filter_GABA_integrated)
xenium_obj_filter_GABA_integrated <- ScaleData(xenium_obj_filter_GABA_integrated, verbose=TRUE, assay="XeniumIntegrated", features=rownames(xenium_obj_filter_GABA_integrated$XeniumIntegrated@data))


xenium_obj_filter_GABA_integrated <- RunPCA(xenium_obj_filter_GABA_integrated, assay="XeniumIntegrated", npcs = nDims, verbose = FALSE, 
                                            reduction.name = "xeniumintegrated_pca", 
                                            reduction.key = "xeniumintegratedpc_")
xenium_obj_filter_GABA_integrated <- RunUMAP(xenium_obj_filter_GABA_integrated, reduction = "xeniumintegrated_pca", dims = 1:nDims, 
                                             reduction.name = "xeniumintegrated_umap", 
                                             reduction.key = "xeniumintegratedumap_", 
                                             assay = "XeniumIntegrated")
xenium_obj_filter_GABA_integrated <- FindNeighbors(xenium_obj_filter_GABA_integrated, assay="XeniumIntegrated", reduction = "xeniumintegrated_pca", dims = 1:nDims, force.recalc=T, 
                                                   graph.name = c("xeniumintegrated_nn", "xeniumintegrated_snn"))

xenium_obj_filter_GABA_integrated <- FindClusters(xenium_obj_filter_GABA_integrated, assay="XeniumIntegrated", resolution = 0.3, graph.name = "xeniumintegrated_snn")
xenium_obj_filter_GABA_integrated@meta.data %>% head()
xenium_obj_filter_GABA_integrated$xeniumintegrated_snn_res.0.5 -> xenium_obj_filter_GABA_integrated$seurat_clusters

DimPlot(xenium_obj_filter_GABA_integrated, group.by = "seurat_clusters", label = TRUE, reduction = "xeniumintegrated_umap", repel = TRUE)
DimPlot(xenium_obj_filter_GABA_integrated, group.by = "orig.ident", label = TRUE, reduction = "xeniumintegrated_umap", repel = TRUE)
save(xenium_obj_filter_GABA_integrated, file = "~/Downloads/xenium_obj_filter_GABA_integrated.rda")

xenium_obj_filtered_GABA = subset(xenium_obj_filtered_GABA, cell = cells_)
all(colnames(xenium_obj_filtered_GABA) == colnames(xenium_obj_filter_GABA_integrated))
xenium_obj_filtered_GABA$integrated_cluster = xenium_obj_filter_GABA_integrated$seurat_clusters

markers_SCN_neurons = FindAllMarkers(xenium_obj_filtered_GABA, only.pos = FALSE, group.by = "integrated_cluster")
markers_SCN_neurons %>% dplyr::filter(cluster == "2", p_val_adj < 0.05)

table(xenium_obj_filtered_GABA$orig.ident, xenium_obj_filtered_GABA$integrated_cluster)

levels(markers_SCN_neurons$cluster) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(cluster_){
    markers_SCN_neurons %>% dplyr::filter(cluster == cluster_) %>% 
      dplyr::arrange(p_val_adj) %>% 
      dplyr::filter(, p_val_adj < 0.05, avg_log2FC > 0) -> df_
  }) -> markers_SCN_neurons_list
markers_SCN_neurons_list
markers_SCN_neurons_list = markers_SCN_neurons_list[c(1:4)]

markers_SCN_neurons %>% dplyr::filter(cluster == "2")


markers_SCN_neurons_list$`0` %>% dplyr::filter(gene %in% c("Vip", "Grp", "Nms", "Avp", "C1ql3", "Bdnf", "Cck", "lrs4"))
markers_SCN_neurons_list$`1` %>% dplyr::filter(gene %in% c("Vip", "Grp", "Nms", "Avp", "C1ql3", "Bdnf", "Cck", "lrs4"))
markers_SCN_neurons_list$`2` %>% dplyr::filter(gene %in% c("Vip", "Grp", "Nms", "Avp", "C1ql3", "Bdnf", "Cck", "lrs4"))
markers_SCN_neurons_list$`3` %>% dplyr::filter(gene %in% c("Vip", "Grp", "Nms", "Avp", "C1ql3", "Bdnf", "Cck", "lrs4"))

xenium_obj_filtered_GABA@meta.data %>% dplyr::filter(integrated_cluster %in% c("0", "1", "3")) %>%rownames() -> cells_
xenium_obj_filtered_GABA = subset(xenium_obj_filtered_GABA, cell = cells_)
xenium_obj_filtered_GABA@meta.data %>% 
  dplyr::mutate(annot = case_when(
    integrated_cluster == "0" ~ "Avp_GABA",
    integrated_cluster == "1" ~ "Vip_Nmns_Clql3_GABA",
    integrated_cluster == "3" ~ "Grp_GABA",
    TRUE ~ NA
  )) %>% drop_na() %>% rownames_to_column("cell") %>% 
  dplyr::select(cell, annot) -> GABA_subclass

save(GABA_subclass, file = "~/Downloads/GABA_subclass.rda")

#### ----
load("~/Downloads/GABA_subclass.rda")
left_join(
  x = xenium_obj_filtered@meta.data %>% rownames_to_column("cell"),
  y = GABA_subclass, 
  by = "cell"
) %>% dplyr::filter(!((predicted.id == "GABAergic")&(is.na(annot)))) %>% 
  dplyr::mutate(annot = case_when(
    is.na(annot) ~ predicted.id,
    TRUE ~ annot
  )) %>% 
  dplyr::mutate(predicted.id = annot) -> meta_tmp
meta_tmp %>% column_to_rownames("cell") -> meta_tmp

xenium_obj_filtered = subset(xenium_obj_filtered, cell = rownames(meta_tmp))
xenium_obj_filtered@meta.data <- meta_tmp

DefaultAssay(xenium_obj_filtered)

DefaultAssay(xenium_obj_filtered) <- "Xenium"
GetAssayData(xenium_obj_filtered, assay = "Xenium", layer = "counts") -> mat_
rowSums(mat_) %>% {.[.!=0]} %>% names() -> genes_expressed
rm(mat_);gc()
length(genes_expressed)

dim(xenium_obj_filtered)
xenium_obj_filtered = xenium_obj_filtered[genes_expressed, ]
dim(xenium_obj_filtered)

xenium_obj_filtered@meta.data %>% dplyr::mutate(celltype = case_when(
  cluster == "5" ~ "Oligodendrocyte", 
  cluster == "3" ~ "Astrocyte", 
  cluster == "8" ~ "Ependymal", 
  TRUE ~ predicted.id
)) -> xenium_obj_filtered@meta.data
xenium_obj_filtered$celltype %>% table()

xenium_obj_filtered$celltype %>% table() %>% {.[. > 100]} %>% names() -> celltype_kept
xenium_obj_filtered@meta.data %>% dplyr::filter(celltype %in% celltype_kept) %>% rownames() -> cell_kept
xenium_obj_filtered = subset(xenium_obj_filtered, cells = cell_kept)

Idents(xenium_obj_filtered) = xenium_obj_filtered$celltype

xenium_obj_filtered@meta.data %>% dplyr::filter(grepl("CT-8", orig.ident)) %>% rownames() -> CT8_cells
xenium_obj_filtered@meta.data %>% dplyr::filter(grepl("CT-20", orig.ident)) %>% rownames() -> CT20_cells

xenium_obj_filtered_CT8 = subset(xenium_obj_filtered, cells = CT8_cells)
xenium_obj_filtered_CT20 = subset(xenium_obj_filtered, cells = CT20_cells)
rm(xenium_obj_filtered)
gc()

xenium_obj_filtered_CT8@meta.data %>% 
  dplyr::mutate(celltype = case_when(
    grepl("Vip", celltype) ~ "Vip_Nms_C1ql3_GABA",
    TRUE ~ celltype
  )) -> xenium_obj_filtered_CT8@meta.data
table(xenium_obj_filtered_CT8$celltype)
xenium_obj_filtered_CT8@meta.data %>% 
  {
    df_ = .
    data.frame(labels = df_$celltype, samples = df_$orig.ident) -> df_1
    rownames(df_1) = rownames(df_)
    
    df_1$labels = factor(df_1$labels, levels = c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular"))
    df_1$samples = factor(df_1$samples, levels = c("CT-8-1", "CT-8-2", "CT-8-3"))
    
    df_1$labels %>% unique()
    df_1$samples %>% unique()
    
    df_1
  } -> meta_CT8

xenium_obj_filtered_CT20@meta.data %>% 
  dplyr::mutate(celltype = case_when(
    grepl("Vip", celltype) ~ "Vip_Nms_C1ql3_GABA",
    TRUE ~ celltype
  )) -> xenium_obj_filtered_CT20@meta.data
table(xenium_obj_filtered_CT20$celltype)
xenium_obj_filtered_CT20@meta.data %>% 
  {
    df_ = .
    data.frame(labels = df_$celltype, samples = df_$orig.ident) -> df_1
    rownames(df_1) = rownames(df_)
    
    df_1$labels = factor(df_1$labels, levels = c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular"))
    df_1$samples = factor(df_1$samples, levels = c("CT-20-1", "CT-20-2", "CT-20-3"))
    
    df_1$labels %>% unique()
    df_1$samples %>% unique()
    
    df_1
  } -> meta_CT20

coords = do.call(rbind, xy_coordinates) %>% "rownames<-"(., NULL) %>% column_to_rownames("cell")
coords.use_CT8 = coords[colnames(xenium_obj_filtered_CT8), ]
coords.use_CT20 = coords[colnames(xenium_obj_filtered_CT20), ]

options(future.globals.maxSize = 8 * 1024^3)  # 4 GiB

DefaultAssay(xenium_obj_filtered_CT8)
xenium_obj_filtered_CT8 <- NormalizeData(
  xenium_obj_filtered_CT8,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

DefaultAssay(xenium_obj_filtered_CT20)
xenium_obj_filtered_CT20 <- NormalizeData(
  xenium_obj_filtered_CT20,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

conversion.factor = 1
spot.size = 10 # use the typical human cell size
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

rbind(spatial.factors, spatial.factors, spatial.factors) %>% "rownames<-"(c("CT-8-1", "CT-8-2", "CT-8-3")) -> spatial.factors_CT8
rbind(spatial.factors, spatial.factors, spatial.factors) %>% "rownames<-"(c("CT-20-1", "CT-20-2", "CT-20-3")) -> spatial.factors_CT20


cellchat_CT8 <- SpatialCellChat::createSpatialCellChat(
  object = xenium_obj_filtered_CT8,
  group.by = "labels",
  assay = "Xenium",
  datatype = "spatial",
  coordinates = coords.use_CT8,
  spatial.factors = spatial.factors_CT8, 
  meta = meta_CT8
)

cellchat_CT20 <- SpatialCellChat::createSpatialCellChat(
  object = xenium_obj_filtered_CT20,
  group.by = "labels",
  assay = "Xenium",
  datatype = "spatial",
  coordinates = coords.use_CT20,
  spatial.factors = spatial.factors_CT20, 
  meta = meta_CT20
)

CellChatDB <- CellChatDB.mouse
cellchat_CT8@DB <- CellChatDB
cellchat_CT20@DB <- CellChatDB

cellchat_CT8
cellchat_CT20

cellchat_CT8 = SpatialCellChat::subsetData(cellchat_CT8, features = genes_expressed)
cellchat_CT20 = SpatialCellChat::subsetData(cellchat_CT20, features = genes_expressed)

cellchat_CT8 <- CellChat::identifyOverExpressedGenes(cellchat_CT8, do.DE = FALSE, min.cells = 10)
cellchat_CT20 <- CellChat::identifyOverExpressedGenes(cellchat_CT20, do.DE = FALSE, min.cells = 10)

cellchat_CT8 <- SpatialCellChat::identifyOverExpressedInteractions(cellchat_CT8, variable.both = F)
cellchat_CT20 <- SpatialCellChat::identifyOverExpressedInteractions(cellchat_CT20, variable.both = F)

future::plan("multisession", workers = 8) 
cellchat_CT8 <- SpatialCellChat::computeCommunProb(
  cellchat_CT8, 
  distance.use = T,
  scale.distance = 100
)
cellchat_CT20 <- SpatialCellChat::computeCommunProb(
  cellchat_CT20, 
  distance.use = T,
  scale.distance = 100
)

cellchat_CT8 <- SpatialCellChat::filterProbability(cellchat_CT8)
cellchat_CT20 <- SpatialCellChat::filterProbability(cellchat_CT20)

cellchat_CT8 <- SpatialCellChat::filterCommunication(
  cellchat_CT8, 
  min.cells = 3,
  min.links = 3,
  min.cells.sr = 3
)
cellchat_CT20 <- SpatialCellChat::filterCommunication(
  cellchat_CT20, 
  min.cells = 3,
  min.links = 3,
  min.cells.sr = 3
)

cellchat_CT8 <- SpatialCellChat::computeAvgCommunProb(
  cellchat_CT8,
  nboot = 100
)
cellchat_CT20 <- SpatialCellChat::computeAvgCommunProb(
  cellchat_CT20,
  nboot = 100
)

cellchat_CT8 <- SpatialCellChat::computeCommunProbPathway(cellchat_CT8)
cellchat_CT20 <- SpatialCellChat::computeCommunProbPathway(cellchat_CT20)

cellchat_CT8 <- SpatialCellChat::aggregateNet(cellchat_CT8)
cellchat_CT20 <- SpatialCellChat::aggregateNet(cellchat_CT20)

save(cellchat_CT8, cellchat_CT20, file = "~/Downloads/merged_cellchat_results_v3.RData")

set.seed(1234)
future.seed = TRUE
cellchat_CT8 <- CellChat::netAnalysis_computeCentrality(cellchat_CT8, slot.name = "netP")
cellchat_CT20 <- CellChat::netAnalysis_computeCentrality(cellchat_CT20, slot.name = "netP")

CellChatDB.mouse$interaction$pathway_name %>% unique() -> pathways
SpatialCellChat::subsetCommunication(cellchat_CT8, thresh = 1.1) -> communication_df_CT8
SpatialCellChat::subsetCommunication(cellchat_CT20, thresh = 1.1) -> communication_df_CT20
cellchat_CT8@netP$pathways -> CT8_netP
cellchat_CT20@netP$pathways -> CT20_netP

####
safe_subsetCommunication_netP <- function(cellchat, netP) {
  empty_df <- data.frame(
    source = character(),
    target = character(),
    pathway_name = character(),
    prob = numeric(),
    pval = numeric(),
    stringsAsFactors = FALSE
  )
  tryCatch(
    {
      df <- subsetCommunication(
        cellchat,
        slot.name = "netP",
        signaling = netP
      )
      return(df)
    },
    error = function(e) {
      message("subsetCommunication failed for: ", netP)
      message("Error message: ", e$message)
      return(empty_df)
    }
  )
}

draw_heatmap = function(cellchat, netP, pval_ = 1){
  #  subsetCommunication(cellchat, slot.name = "netP", signaling = netP) -> df_
  safe_subsetCommunication_netP(cellchat = cellchat, netP = netP) -> df_
  celltype = c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular")
  source = celltype %>% purrr::map(function(x){rep(x, length(celltype))}) %>% unlist()
  data.frame(source = source, target = celltype) -> df_1
  left_join(x = df_1, y = df_, by = c("source", "target")) %>% 
    dplyr::select(source, target, prob, pval) %>% 
    dplyr::mutate(prob = case_when(
      is.na(prob) ~ 0, 
      TRUE ~ prob
    )) %>% 
    dplyr::mutate(pval = case_when(
      is.na(pval) ~ 1, 
      TRUE ~ pval
    )) -> df_plot
  
  df_plot$source = factor(df_plot$source, levels = rev(celltype))
  df_plot$target = factor(df_plot$target, levels = celltype)
  
  df_plot %>% 
    dplyr::mutate(prob = case_when(
      pval >= pval_ ~ 0,
      TRUE ~ prob
    )) -> df_plot
  
  df_plot %>% 
    dplyr::mutate(pval_1 = case_when(
      pval < pval_ ~ sprintf("p=%s", round(pval,3)),
      TRUE ~ ""
    )) -> df_plot
  
  max_prob <- max(df_plot$prob, na.rm = TRUE)
  
  df_plot %>% 
    ggplot(aes(x = target, y = source, fill = prob)) +
#    geom_tile(color = "black") + 
    geom_tile() + 
    geom_text(aes(label = pval_1), size = 3, color = "black") +
    scale_fill_gradient(
      low = "white",
      high = "red",
      limits = c(0, ifelse(max_prob == 0, 1, max_prob)),
      na.value = "white"
    ) + 
#    theme_bw() + 
    theme_classic() + 
    ggtitle(netP) + 
    labs(fill = "communication\nprobability") + 
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )
    )
}

plot_dotplot = function(cellchat, signaling, pval_){
  res = CellChat::extractEnrichedLR(cellchat, signaling = signaling, geneLR.return = T, enriched.only = F, thresh = pval_)
  res$geneLR -> pathway_genes
  res = CellChat::extractEnrichedLR(cellchat, signaling = signaling, geneLR.return = T, enriched.only = T, thresh = pval_)
  if (nrow(res$pairLR) == 0){
    enriched_genes = c("")
  }else{
    enriched_genes = res$geneLR
  }
  
  cellchat@meta -> meta
  cellchat@data[pathway_genes, ] -> mat
  celltypes = c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular")
  
  celltypes %>% 
    "names<-"(.,.) %>% 
    #    .[1] %>% 
    purrr::map(function(celltype_){
      meta %>% dplyr::filter(labels == celltype_) %>% rownames() -> cells_
      mat[,cells_] -> mat
      rowMeans(mat) -> norm_expr
      apply(mat, 1, function(x){
        non_zero_cell = sum(x != 0)
        non_zero_prop = (non_zero_cell/length(cells_))*100
      }) -> non_zero_prop
      data.frame(expr = norm_expr, percent = non_zero_prop) %>% 
        rownames_to_column("gene") %>% dplyr::mutate(celltype = celltype_)
    }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_plot
  df_plot$celltype = factor(df_plot$celltype, levels = celltypes)
  df_plot$gene = factor(df_plot$gene, levels = pathway_genes)
  df_plot %>% 
    ggplot(aes(x = celltype, y = gene, size = percent, color = expr)) + 
    geom_point() + 
    scale_y_discrete(
      labels = function(x) {
        ifelse(x %in% enriched_genes, paste0("<b>", x, "</b>"), x)
      }
    ) +
    scale_color_gradient(low = "white", high = "red") + 
    theme_bw() + 
    theme(
      axis.text.y = ggtext::element_markdown()
    )
}

plot_violin = function(signaling, pval_){
  list(CT8 = cellchat_CT8, CT20 = cellchat_CT20) %>% 
    #    .[1] %>% 
    purrr::map2(.x=.,.y=names(.),.f=function(cellchat, timepoint){
      res = CellChat::extractEnrichedLR(cellchat, signaling = signaling, geneLR.return = T, enriched.only = F, thresh = pval_)
      res$geneLR -> pathway_genes
      res$geneLR ->> pathway_genes
      res = CellChat::extractEnrichedLR(cellchat, signaling = signaling, geneLR.return = T, enriched.only = T, thresh = pval_)
      if (nrow(res$pairLR) == 0){
        enriched_genes = c("")
        enriched_genes <<- c("")
      }else{
        enriched_genes = res$geneLR
        enriched_genes <<- res$geneLR
      }
      
      cellchat@meta -> meta
      cellchat@data[pathway_genes, ] -> mat
      celltypes = c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular")
      
      celltypes %>% 
        "names<-"(.,.) %>% 
        #    .[1] %>% 
        purrr::map(function(celltype_){
          meta %>% dplyr::filter(labels == celltype_) %>% rownames() -> cells_
          mat[,cells_] -> mat
          mat %>% as.data.frame() %>% 
            rownames_to_column("gene") %>% 
            pivot_longer(cols = -gene, names_to = "cell", values_to = "expr") %>% 
            dplyr::mutate(celltype = celltype_, timepoint = timepoint)
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_plot
  pathway_genes %>% 
    #    .[1] %>% 
    "names<-"(.,.) %>% 
    purrr::map(function(gene_){
      df_plot %>% dplyr::filter(gene == gene_) -> df_plot
      celltypes = c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular")
      df_plot$celltype = factor(df_plot$celltype, levels = celltypes)
      df_plot$timepoint = factor(df_plot$timepoint, levels = c("CT8", "CT20"))
      df_plot %>% 
        ggplot(aes(x = celltype, y = expr, fill = timepoint)) + 
        geom_violin() + 
        theme_classic() + 
        ylab(gene_) +
        theme(
          axis.title.y = element_text(
            angle = 0,
            vjust = 0.5,
            hjust = 1
          )
        ) -> p
      if (gene_ %in% enriched_genes){
        p = p + 
          theme(
            axis.title.y = element_text(face = "bold")
          )
#        print(gene_)
      }
      
      p + 
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
        ) -> p
      
      return(p)
    }) -> p_list
  patchwork::wrap_plots(p_list, ncol = 1, guides = "collect")
}

#GABA
subsetCommunication(cellchat_CT8, slot.name = "netP", signaling = "GABA-A") 
subsetCommunication(cellchat_CT20, slot.name = "netP", signaling = "GABA-A")

draw_heatmap(cellchat = cellchat_CT8, netP = "GABA-A", pval_ = 0.01) + ggtitle("GABA-A (CT-8)") -> p1
draw_heatmap(cellchat = cellchat_CT20, netP = "GABA-A", pval_ = 0.01) + ggtitle("GABA-A (CT-20)") -> p2
draw_heatmap(cellchat = cellchat_CT8, netP = "GABA-B", pval_ = 0.01) + ggtitle("GABA-B (CT-8)") -> p3
draw_heatmap(cellchat = cellchat_CT20, netP = "GABA-B", pval_ = 0.01) + ggtitle("GABA-B (CT-20)") -> p4
list(
  p1,p2,p3,p4
) %>% 
  purrr::map(function(p_){
    p_$data$prob
  }) %>% unlist() %>% max() -> max_prob

i = 0
list(
  p1,p2,p3,p4
) %>% purrr::map(function(p_){
  i <<- i + 1
  p_ + 
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, max_prob),
    na.value = "white"
  ) -> p
  if (i < 3){
    p + theme(axis.text.x = element_blank()) + xlab(NULL) -> p
    print(i)
  }
  return(p)
}) %>% patchwork::wrap_plots(.,ncol = 2, guides = "collect") -> p_GABA_heatmap

plot_violin(signaling = "GABA-A", pval_ = 0.01) -> p
plot_violin(signaling = "GABA-B", pval_ = 0.01) -> p1
1:2 %>% 
  purrr::map(function(i){
    p[[i]] -> p
    if (i != 2){
      p +
        theme(
          axis.text.x = element_blank()
        ) + 
        xlab(NULL)
    } -> p
    return(p)
  }) %>% patchwork::wrap_plots(., ncol = 1, guides = "collect") + 
  patchwork::plot_annotation(title = "GABA transportors expression") -> p_GABA_violoin_Slc

4:5 %>% 
  purrr::map(function(i){
    p[[i]] -> p
    if (i != 5){
      p +
        theme(
          axis.text.x = element_blank()
        ) + 
        xlab(NULL)
    } -> p
    return(p)
  }) %>% patchwork::wrap_plots(., ncol = 1, guides = "collect") + 
  patchwork::plot_annotation(title = "GABA enzymes expression")-> p_GABA_violoin_GAD

c(6,8,10,11,12) %>% 
  purrr::map(function(i){
    p[[i]] -> p
    if (i != 12){
      p +
        theme(
          axis.text.x = element_blank()
        ) + 
        xlab(NULL)
    } -> p
    return(p)
  }) %>% patchwork::wrap_plots(., ncol = 1, guides = "collect") + 
  patchwork::plot_annotation(title = "GABA-A receptors expression")-> p_GABA_violoin_GABA_A_receptor


c(6,7) %>% 
  purrr::map(function(i){
    p1[[i]] -> p
    if (i != 7){
      p +
        theme(
          axis.text.x = element_blank()
        ) + 
        xlab(NULL)
    } -> p
    return(p)
  }) %>% patchwork::wrap_plots(., ncol = 1, guides = "collect") + 
  patchwork::plot_annotation(title = "GABA-B receptors expression")-> p_GABA_violoin_GABA_B_receptor



calc_distance <- function(x1, y1, x2, y2) {
  
  sqrt((x2 - x1)^2 + (y2 - y1)^2)
  
}
cellchat_CT8@netP$prob.cell[,, "GABA-B"] %>% 
  as.data.frame() %>% 
  rownames_to_column("source") %>% 
  pivot_longer(cols = -source, names_to = "target", values_to = "prob") %>% {
    df_prob = .
    cellchat_CT8@meta %>% rownames_to_column("cell") %>% "colnames<-"(., c("cell", "source_labels", "source_samples")) -> df_meta_source
    cellchat_CT8@meta %>% rownames_to_column("cell") %>% "colnames<-"(., c("cell", "target_labels", "target_samples")) -> df_meta_target
    
    left_join(x = df_prob, y = df_meta_source, by = c("source"="cell")) %>% 
      left_join(x = ., y = df_meta_target, by = c("target"="cell")) -> df_
    
    cellchat_CT8@images$result.computeCellDistance$d.spatial %>% 
      as.data.frame() %>% 
      "rownames<-"(., df_meta_source$cell) %>% 
      "colnames<-"(., df_meta_target$cell) -> df_dist
    df_dist %>% 
      rownames_to_column("source") %>% 
      pivot_longer(cols = -source, names_to = "target", values_to = "dist") %>% 
      distinct() -> df_dist
    
    cellchat_CT8@images$coordinates %>% "colnames<-"(., c("source_x", "source_y")) %>% rownames_to_column("cell") -> df_coord_source
    cellchat_CT8@images$coordinates %>% "colnames<-"(., c("target_x", "target_y")) %>% rownames_to_column("cell") -> df_coord_target
    
    left_join(x = df_, y = df_coord_source, by = c("source"="cell")) %>% 
      left_join(x = ., y = df_coord_target, by = c("target"="cell")) -> df_
    
    left_join(x = df_, y = df_dist, by = c("source", "target")) -> df_
    
    df_ %>% dplyr::filter(source_samples != target_samples) %>% dplyr::arrange(desc(prob))
  }
calc_distance(x1 = 4581, y1 = 887, x2 = 4646, y2 = 619)


netVisual_heatmap(object = cellchat_CT8, signaling = "GABA-A")
netVisual_heatmap(object = cellchat_CT20, signaling = "EPHA")
plot_violin(signaling = "GABA-A", pval_ = 0.01)
#SST
subsetCommunication(cellchat_CT8, slot.name = "netP", signaling = "SOMATOSTATIN")
subsetCommunication(cellchat_CT20, slot.name = "netP", signaling = "SOMATOSTATIN")
draw_heatmap(cellchat = cellchat_CT8, netP = "SOMATOSTATIN", pval_ = 0.01)
draw_heatmap(cellchat = cellchat_CT20, netP = "SOMATOSTATIN", pval_ = 0.01)
netVisual_heatmap(object = cellchat_CT20, signaling = "SOMATOSTATIN")
netVisual_heatmap(object = cellchat_CT8, signaling = "SOMATOSTATIN")

subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "SOMATOSTATIN")
subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "SOMATOSTATIN")

#EPHA/EPHB
subsetCommunication(cellchat_CT8, slot.name = "netP", signaling = "EPHA", thresh = 0.01)
subsetCommunication(cellchat_CT20, slot.name = "netP", signaling = "EPHA", thresh = 0.01)
subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "EPHA")
subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "EPHA")
draw_heatmap(cellchat = cellchat_CT8, netP = "EPHA", pval_ = 0.01)
draw_heatmap(cellchat = cellchat_CT20, netP = "EPHA", pval_ = 0.01)
netVisual_heatmap(object = cellchat_CT8, signaling = "EPHA")
netVisual_heatmap(object = cellchat_CT20, signaling = "EPHA")
plot_dotplot(cellchat = cellchat_CT8, signaling = "EPHA", pval_ = 0.01)
plot_dotplot(cellchat = cellchat_CT20, signaling = "EPHA", pval_ = 0.01)
plot_violin(signaling = "EPHA", pval_ = 0.01)


rbind(
communication_df_CT8 %>% 
  dplyr::filter(pval < 0.05) %>% 
  dplyr::mutate(timepoint = "CT8"),
communication_df_CT20 %>% 
  dplyr::filter(pval < 0.05) %>% 
  dplyr::mutate(timepoint = "CT20")
) %>% 
  {
    df_ = .
    df_ %>% dplyr::mutate(group = sprintf("%s_%s", ligand, receptor)) -> df_
    df_$group %>% unique() %>% 
      "names<-"(.,.) %>% 
      .[1:5] %>% 
      purrr::map(function(group_){
        df_ %>% dplyr::filter(group == group_)
      })
  }

