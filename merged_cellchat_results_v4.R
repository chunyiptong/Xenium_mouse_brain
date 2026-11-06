library(tidyverse)
library(Seurat)
#library(SpatialCellChat)
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

GetAssayData(object = xenium_obj_filtered, assay = "Xenium", layer = "counts") -> mat_
mat_ %>% as.data.frame() %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6})) -> mat_norm

xenium_obj_filtered@meta.data$annot %>% 
  unique() %>% {.[grepl("GABA",.)]} %>% 
  purrr::map(function(annot_){
    xenium_obj_filtered@meta.data %>% dplyr::filter(annot == annot_) %>% rownames() -> cells_
    genes = c("Avp", "Vip", "Nms", "C1ql3", "Grp")
    mat_norm[genes, cells_] -> mat_norm 
    mat_norm %>% rowMeans() -> mean_expr
    apply(mat_norm, 1, function(x){
      non_zero_cells = (x != 0) %>% sum()
      non_zero_prop = non_zero_cells/length(cells_)
    }) -> non_zero_prop
    data.frame(gene = genes, norm_expr = mean_expr, percent = non_zero_prop, celltype = annot_)
  }) %>% do.call(rbind, .) -> df_plot

df_plot %>% 
  group_by(gene) %>% 
  group_map(function(df_,y){
    df_ %>% dplyr::mutate(scaled_expr = scales::rescale(x = norm_expr, to = c(0,5)))
  },.keep = T) %>% do.call(rbind, .) %>% 
  dplyr::mutate(gene = factor(gene, levels = c("Avp", "Vip", "Nms", "C1ql3", "Grp"))) %>% 
  ggplot(aes(x = gene, y = celltype, size = percent, color = scaled_expr)) + 
  geom_point() + 
  scale_color_gradient(low = "blue", high = "red") + 
  theme_classic() -> p

ggsave(plot = p, filename = "~/Downloads/GABA_neuron_subtype_dotplot.pdf", height = 4.16, width = 4.64)
  

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

cellchat_CT8 <- CellChat::createCellChat(
  object = xenium_obj_filtered_CT8,
  group.by = "labels",
  assay = "Xenium",
  datatype = "spatial",
  coordinates = coords.use_CT8,
  spatial.factors = spatial.factors_CT8, 
  meta = meta_CT8
)

cellchat_CT20 <- CellChat::createCellChat(
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

cellchat_CT8 = CellChat::subsetData(cellchat_CT8, features = genes_expressed)
cellchat_CT20 = CellChat::subsetData(cellchat_CT20, features = genes_expressed)

cellchat_CT8 <- CellChat::identifyOverExpressedGenes(cellchat_CT8, do.DE = FALSE, min.cells = 10)
cellchat_CT20 <- CellChat::identifyOverExpressedGenes(cellchat_CT20, do.DE = FALSE, min.cells = 10)

cellchat_CT8 <- CellChat::identifyOverExpressedInteractions(cellchat_CT8, variable.both = F)
cellchat_CT20 <- CellChat::identifyOverExpressedInteractions(cellchat_CT20, variable.both = F)

future::plan("multisession", workers = 8) 
cellchat_CT8 <- CellChat::computeCommunProb(
  cellchat_CT8, 
  distance.use = T,
  contact.range = 10,
  interaction.range = 250,
  scale.distance = 4
)
cellchat_CT20 <- CellChat::computeCommunProb(
  cellchat_CT20, 
  distance.use = T,
  contact.range = 10,
  interaction.range = 250,
  scale.distance = 4
)

cellchat_CT8 <- CellChat::filterCommunication(
  cellchat_CT8, 
  min.cells = 5,
)
cellchat_CT20 <- CellChat::filterCommunication(
  cellchat_CT20, 
  min.cells = 5,
)

cellchat_CT8 <- CellChat::computeCommunProbPathway(cellchat_CT8)
cellchat_CT20 <- CellChat::computeCommunProbPathway(cellchat_CT20)

cellchat_CT8 <- CellChat::aggregateNet(cellchat_CT8)
cellchat_CT20 <- CellChat::aggregateNet(cellchat_CT20)

save(cellchat_CT8, cellchat_CT20, file = "~/Downloads/merged_cellchat_results_v4.RData")

set.seed(1234)
future.seed = TRUE
cellchat_CT8 <- CellChat::netAnalysis_computeCentrality(cellchat_CT8, slot.name = "netP")
cellchat_CT20 <- CellChat::netAnalysis_computeCentrality(cellchat_CT20, slot.name = "netP")

#Functions ----
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

safe_subsetCommunication_net <- function(cellchat, netP) {
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
        slot.name = "net",
        pairLR.use = data.frame(interaction_name = netP)
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
  
#  df_plot %>% 
#    dplyr::mutate(prob = case_when(
#      pval >= pval_ ~ 0,
#      TRUE ~ prob
#    )) -> df_plot
  
  df_plot %>% 
    dplyr::mutate(pval_1 = case_when(
#      pval < pval_ ~ sprintf("p=%s", round(pval,3)),
      pval < pval_ ~ "*",
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

draw_heatmap_net = function(cellchat, netP, pval_ = 1){
  safe_subsetCommunication_net(cellchat = cellchat, netP = netP) -> df_
  
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
  
#  df_plot %>% 
#    dplyr::mutate(prob = case_when(
#      pval >= pval_ ~ 0,
#      TRUE ~ prob
#    )) -> df_plot
  
  df_plot %>% 
    dplyr::mutate(pval_1 = case_when(
#      pval < pval_ ~ sprintf("p=%s", round(pval,3)),
      pval < pval_ ~ "*",
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
          mat %>% as.matrix() %>% 
            as.data.frame() %>% 
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

#custom heatmap function ----
source("~/Downloads/netVisual_heatmap_custom.R")
#analsyis ----
pdf(file = "~/Downloads/communication_count_heatmap.pdf", width = 5.75, height = 4.94)
netVisual_heatmap(cellchat_CT8, measure = "count", color.heatmap = "Reds", slot.name = "net")
netVisual_heatmap_custom(cellchat_CT8, measure = "count", color.heatmap = "Reds", slot.name = "net", heatmap.max = 60)
cellchat_CT8@net$count %>% 
  as.data.frame() %>% 
  rownames_to_column("source") %>% 
  pivot_longer(cols = -source, names_to = "target", values_to = "interaction") %>% 
  dplyr::mutate(source = case_when(
    grepl("GABA", source) ~ "GABAergic",
    TRUE ~ source
  )) %>% 
  dplyr::mutate(target = case_when(
    grepl("GABA", target) ~ "GABAergic",
    TRUE ~ target
  )) %>% 
  group_by(source, target) %>% 
  summarise(interaction = sum(interaction)) %>% 
  dplyr::arrange(desc(interaction)) -> df_plot

df_plot %>% 
  dplyr::mutate(communication = sprintf("%s -> %s", source, target)) %>% 
  {
    df_ = .
    df_$communication = factor(df_$communication, levels = rev(df_$communication))
    df_
  } %>% 
  ggplot(aes(y = communication, x = interaction)) + 
  geom_col(fill = "red") + 
  theme_classic() + 
  xlab("Number of interactions") + 
  ggtitle("CT8") -> p1

netVisual_heatmap(cellchat_CT20, measure = "count", color.heatmap = "Blues", slot.name = "net")
netVisual_heatmap_custom(cellchat_CT20, measure = "count", color.heatmap = "Reds", slot.name = "net", heatmap.max = 60)
dev.off()
cellchat_CT20@net$count %>% 
  as.data.frame() %>% 
  rownames_to_column("source") %>% 
  pivot_longer(cols = -source, names_to = "target", values_to = "interaction") %>% 
  dplyr::mutate(source = case_when(
    grepl("GABA", source) ~ "GABAergic",
    TRUE ~ source
  )) %>% 
  dplyr::mutate(target = case_when(
    grepl("GABA", target) ~ "GABAergic",
    TRUE ~ target
  )) %>% 
  group_by(source, target) %>% 
  summarise(interaction = sum(interaction)) %>% 
  dplyr::arrange(desc(interaction)) -> df_plot

df_plot %>% 
  dplyr::mutate(communication = sprintf("%s -> %s", source, target)) %>% 
  {
    df_ = .
    df_$communication = factor(df_$communication, levels = rev(df_$communication))
    df_
  } %>% 
  ggplot(aes(y = communication, x = interaction)) + 
  geom_col(fill = "blue") + 
  theme_classic() + 
  xlab("Number of interactions") + 
  ggtitle("CT20") -> p2

rbind(
p1$data %>% as.data.frame() %>% dplyr::select(interaction, communication) %>% dplyr::mutate(timepoint = "CT8"),
p2$data %>% as.data.frame() %>% dplyr::select(interaction, communication) %>% dplyr::mutate(timepoint = "CT20")
) %>% dplyr::mutate(timepoint = factor(timepoint, levels = c("CT8", "CT20"))) %>% 
  {
    df_ = .
    df_ %>% group_by(communication) %>% summarise(interaction = sum(interaction)) %>% dplyr::arrange(desc(interaction)) %>% .$communication -> levels_
    df_ %>% 
      dplyr::mutate(communication = factor(communication, levels = rev(levels_))) %>% 
      ggplot(aes(x = interaction, y = communication, group = timepoint, fill = timepoint)) + 
      geom_col() + 
      theme_classic() + 
      xlab("number of communications")
  } -> p3
ggsave(plot = p3, width = 5.75, height = 4.94, filename = "~/Downloads/communication_count_barplot.pdf")

list(
`CT8 (LR signaling pathways)` = cellchat_CT8@netP$pathways,
`CT20 (LR signaling pathways)` = cellchat_CT20@netP$pathways
) %>% ggvenn::ggvenn(fill_color = c("#e37365", "#4cb3bb"), text_size = 8)

list(
  `CT8 (LRs)` = cellchat_CT8@net$LRs,
  `CT20 (LRs)` = cellchat_CT20@net$LRs
) %>% ggvenn::ggvenn(fill_color = c("#e37365", "#4cb3bb"), text_size = 6) -> p4
ggsave(plot = p4, width = 5.75, height = 4.94, filename = "~/Downloads/communication_count_venn.pdf")

#GABA ----
subsetCommunication(cellchat_CT8, slot.name = "netP", signaling = "GABA-A") 
subsetCommunication(cellchat_CT20, slot.name = "netP", signaling = "GABA-A")
subsetCommunication(cellchat_CT8, slot.name = "netP", signaling = "GABA-B") 
subsetCommunication(cellchat_CT20, slot.name = "netP", signaling = "GABA-B")

c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(celltype_){
    subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "GABA-A") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(source == celltype_) %>% .$ligand %>% unique() %>% {.[grepl("SLC32A1",.)]} %>% paste(., collapse = " ; ") -> ligand
    subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "GABA-A") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(target == celltype_) %>% .$receptor %>% unique() %>% paste(., collapse = " ; ") -> receptor
    data.frame(celltype = celltype_, ligand = ligand, receptor = receptor)
  }) %>% do.call(rbind ,.) %>% dplyr::mutate(CT = "CT8") %>% "rownames<-"(.,NULL)

c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(celltype_){
    subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "GABA-A") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(source == celltype_) %>% .$ligand %>% unique() %>% {.[grepl("SLC32A1",.)]} %>% paste(., collapse = " ; ") -> ligand
    subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "GABA-A") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(target == celltype_) %>% .$receptor %>% unique() %>% paste(., collapse = " ; ") -> receptor
    data.frame(celltype = celltype_, ligand = ligand, receptor = receptor)
  }) %>% do.call(rbind ,.) %>% dplyr::mutate(CT = "CT20") %>% "rownames<-"(.,NULL)

c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(celltype_){
    subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "GABA-B") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(source == celltype_) %>% .$ligand %>% unique() %>% {.[grepl("SLC32A1",.)]} %>% paste(., collapse = " ; ") -> ligand
    subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "GABA-B") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(target == celltype_) %>% .$receptor %>% unique() %>% paste(., collapse = " ; ") -> receptor
    data.frame(celltype = celltype_, ligand = ligand, receptor = receptor)
  }) %>% do.call(rbind ,.) %>% dplyr::mutate(CT = "CT8") %>% "rownames<-"(.,NULL)

c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(celltype_){
    subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "GABA-B") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(source == celltype_) %>% .$ligand %>% unique() %>% {.[grepl("SLC32A1",.)]} %>% paste(., collapse = " ; ") -> ligand
    subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "GABA-B") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(target == celltype_) %>% .$receptor %>% unique() %>% paste(., collapse = " ; ") -> receptor
    data.frame(celltype = celltype_, ligand = ligand, receptor = receptor)
  }) %>% do.call(rbind ,.) %>% dplyr::mutate(CT = "CT20") %>% "rownames<-"(.,NULL)

c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(celltype_){
    subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "VIP") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(source == celltype_) %>% .$ligand %>% unique()  %>% paste(., collapse = " ; ") -> ligand
    subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "VIP") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(target == celltype_) %>% .$receptor %>% unique() %>% paste(., collapse = " ; ") -> receptor
    data.frame(celltype = celltype_, ligand = ligand, receptor = receptor)
  }) %>% do.call(rbind ,.) %>% dplyr::mutate(CT = "CT8") %>% "rownames<-"(.,NULL)

c("Vip_Nms_C1ql3_GABA", "Avp_GABA", "Grp_GABA", "Astrocyte", "Oligodendrocyte", "Vascular") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(celltype_){
    subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "VIP") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(source == celltype_) %>% .$ligand %>% unique() %>% paste(., collapse = " ; ") -> ligand
    subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "VIP") %>% dplyr::filter(pval <= 0.01) %>% 
      dplyr::filter(target == celltype_) %>% .$receptor %>% unique() %>% paste(., collapse = " ; ") -> receptor
    data.frame(celltype = celltype_, ligand = ligand, receptor = receptor)
  }) %>% do.call(rbind ,.) %>% dplyr::mutate(CT = "CT20") %>% "rownames<-"(.,NULL)


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
ggsave(plot = p_GABA_heatmap, width = 8.96, height = 5.45, filename = "~/Downloads/GABA_A_B.heatmap_V2.pdf")

plot_violin(signaling = "GABA-A", pval_ = 0.01) -> p
plot_violin(signaling = "GABA-B", pval_ = 0.01) -> p1

c(4,5,1) %>% 
  purrr::map(function(i){
    p[[i]] -> p
    if (i != 6){
      p +
        theme(
          axis.text.x = element_blank()
        ) + 
        xlab(NULL)
    } -> p
    return(p)
  }) %>% patchwork::wrap_plots(., ncol = 1, guides = "collect") + 
  patchwork::plot_annotation(title = "GABA LR pathway (enzymes and transporters)") -> p_GABA_1

c(6,10,11,12) %>% 
  purrr::map(function(i){
    p[[i]] -> p
    if (i != 13){
      p +
        theme(
          axis.text.x = element_blank()
        ) + 
        xlab(NULL)
    } -> p
    return(p)
  }) %>% patchwork::wrap_plots(., ncol = 1, guides = "collect") + 
  patchwork::plot_annotation(title = "GABA-A LR pathway (receptors)") -> p_GABA_2

c(6:7) %>% 
  purrr::map(function(i){
    p1[[i]] -> p
    if (i != 8){
      p +
        theme(
          axis.text.x = element_blank()
        ) + 
        xlab(NULL)
    } -> p
    return(p)
  }) %>% patchwork::wrap_plots(., ncol = 1, guides = "collect") + 
  patchwork::plot_annotation(title = "GABA-B LR pathway (receptors)") -> p_GABA_3

patchwork::wrap_plots(p_GABA_1, p_GABA_2, p_GABA_3, ncol = 1, heights = c(3,4,2)) -> p
ggsave(plot = p, filename = "~/Downloads/GABA_LR_pathway_gene_expression.pdf", height = 7.68, width = 6.90)

c(1, 4:6, 11, 10, 12) %>% 
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
  patchwork::plot_annotation(title = "GABA-A LR pathway (gene expression)") -> p_GABA

c(1, 4:7) %>% 
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
  patchwork::plot_annotation(title = "GABA-B LR pathway (gene expression)") -> p_GABA_1

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

c(6,10,11,12) %>% 
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


rbind(
  netAnalysis_contribution(cellchat_CT8, signaling = "GABA-A", thresh = 0.01)$data %>% dplyr::mutate(timepoint = "CT8") %>% dplyr::filter(contribution != 0),
  netAnalysis_contribution(cellchat_CT20, signaling = "GABA-A", thresh = 0.01)$data %>% dplyr::mutate(timepoint = "CT20") %>% dplyr::filter(contribution != 0)
) %>% dplyr::mutate(timepoint = factor(timepoint, levels = c("CT8", "CT20"))) %>% 
  group_by(name) %>% 
  group_map(function(df_, y){
    df_ %>% dplyr::mutate(total_contribution = sum(contribution)) -> df_
    df_ %>% dplyr::mutate(contribution = case_when(
      timepoint == "CT8" ~ -contribution,
      TRUE ~ contribution
    ))
  }, .keep = T) %>% do.call(rbind, .) %>% 
  dplyr::arrange(desc(total_contribution)) %>% 
  head(12) %>% 
  dplyr::mutate(name = factor(name, levels = rev(unique(name)))) %>% 
  ggplot(aes(x = contribution, y = name, fill = timepoint)) + 
  geom_col() + 
  theme_classic() + 
  ylab(NULL) + 
  ggtitle("GABA-A LRs contribution") -> p1

rbind(
  netAnalysis_contribution(cellchat_CT8, signaling = "GABA-B", thresh = 0.01)$data %>% dplyr::mutate(timepoint = "CT8") %>% dplyr::filter(contribution != 0),
  netAnalysis_contribution(cellchat_CT20, signaling = "GABA-B", thresh = 0.01)$data %>% dplyr::mutate(timepoint = "CT20") %>% dplyr::filter(contribution != 0)
) %>% dplyr::mutate(timepoint = factor(timepoint, levels = c("CT8", "CT20"))) %>% 
  group_by(name) %>% 
  group_map(function(df_, y){
    df_ %>% dplyr::mutate(total_contribution = sum(contribution)) -> df_
    df_ %>% dplyr::mutate(contribution = case_when(
      timepoint == "CT8" ~ -contribution,
      TRUE ~ contribution
    ))
  }, .keep = T) %>% do.call(rbind, .) %>% 
  dplyr::arrange(desc(total_contribution)) %>% 
  head(8) %>% 
  dplyr::mutate(name = factor(name, levels = rev(unique(name)))) %>% 
  ggplot(aes(x = contribution, y = name, fill = timepoint)) + 
  geom_col() + 
  theme_classic() + 
  ylab(NULL) + 
  ggtitle("GABA-B LRs contribution") -> p2

patchwork::wrap_plots(p1,p2,ncol = 2, guides = "collect")

rbind(
netAnalysis_contribution(cellchat_CT8, signaling = "GABA-A", thresh = 0.01)$data %>% dplyr::mutate(timepoint = "CT8"),
netAnalysis_contribution(cellchat_CT20, signaling = "GABA-A", thresh = 0.01)$data %>% dplyr::mutate(timepoint = "CT20")
) %>% dplyr::mutate(timepoint = factor(timepoint, levels = c("CT8", "CT20"))) %>% 
  dplyr::arrange(desc(contribution)) %>% 
  head(10) %>% 
  dplyr::mutate(name = sprintf("%s (%s)", name, timepoint)) %>% 
  dplyr::mutate(name = factor(name, levels = rev(name))) %>% 
  ggplot(aes(x = contribution, y = name, fill = timepoint)) + 
    geom_col() + 
  theme_classic() + 
  ylab(NULL) + 
  ggtitle("GABA-A LRs contribution") -> p1

rbind(
  netAnalysis_contribution(cellchat_CT8, signaling = "GABA-B", thresh = 0.01)$data %>% dplyr::mutate(timepoint = "CT8"),
  netAnalysis_contribution(cellchat_CT20, signaling = "GABA-B", thresh = 0.01)$data %>% dplyr::mutate(timepoint = "CT20")
) %>% dplyr::mutate(timepoint = factor(timepoint, levels = c("CT8", "CT20"))) %>% 
  dplyr::arrange(desc(contribution)) %>% 
  head(10) %>% 
  dplyr::mutate(name = sprintf("%s (%s)", name, timepoint)) %>% 
  dplyr::mutate(name = factor(name, levels = rev(name))) %>% 
  ggplot(aes(x = contribution, y = name, fill = timepoint)) + 
  geom_col() + 
  theme_classic() + 
  ylab(NULL) + 
  ggtitle("GABA-B LRs contribution") -> p2

patchwork::wrap_plots(p1,p2,ncol = 2, guides = "collect")

#rbind(
#subsetCommunication(cellchat_CT8, slot.name = "netP", signaling = "GABA-A") %>% dplyr::filter(pval <= 0.01) %>% dplyr::mutate(timepoint = "CT8"),
#subsetCommunication(cellchat_CT20, slot.name = "netP", signaling = "GABA-A") %>% dplyr::filter(pval <= 0.01) %>% dplyr::mutate(timepoint = "CT20"),
#subsetCommunication(cellchat_CT8, slot.name = "netP", signaling = "GABA-B")  %>% dplyr::filter(pval <= 0.01) %>% dplyr::mutate(timepoint = "CT8"),
#subsetCommunication(cellchat_CT20, slot.name = "netP", signaling = "GABA-B") %>% dplyr::filter(pval <= 0.01) %>% dplyr::mutate(timepoint = "CT20")
#) %>% dplyr::select(pathway_name, prob, timepoint) %>% 
#  dplyr::mutate(timepoint = factor(timepoint, levels = c("CT8", "CT20"))) %>% 
#  dplyr::mutate(pathway_name = factor(pathway_name, levels = c("GABA-A", "GABA-B"))) %>% 
#  ggplot(aes(x = pathway_name, y = prob, fill = timepoint)) + 
#  geom_boxplot() + 
#  theme_classic() -> p
#
#p + 
#  ggpubr::stat_compare_means(
#    aes(group = timepoint),
#    method = "t.test",
#    method.args = list(alternative = "greater"),
#    label = "p.format"
#  ) -> p1

rbind(
  subsetCommunication(cellchat_CT8, slot.name = "netP", signaling = "GABA-A") %>% dplyr::filter(pval <= 0.01) %>% dplyr::mutate(timepoint = "CT8"),
  subsetCommunication(cellchat_CT20, slot.name = "netP", signaling = "GABA-A") %>% dplyr::filter(pval <= 0.01) %>% dplyr::mutate(timepoint = "CT20")
) %>% 
  {
    df_ = .
    df_ %>% dplyr::filter(timepoint == "CT8") -> df_CT8
    df_ %>% dplyr::filter(timepoint == "CT20") -> df_CT20
    
    full_join(x = df_CT8, y = df_CT20[,c("source", "target")], by = c("source", "target")) %>% 
      dplyr::mutate(pathway_name = "GABA-A", timepoint = "CT8") %>% 
      dplyr::mutate(prob = case_when(
        is.na(prob) ~ 0,
        TRUE ~ prob
      )) -> df_CT8_1
    
    full_join(x = df_CT20, y = df_CT8[,c("source", "target")], by = c("source", "target")) %>% 
      dplyr::mutate(pathway_name = "GABA-A", timepoint = "CT20") %>% 
      dplyr::mutate(prob = case_when(
        is.na(prob) ~ 0,
        TRUE ~ prob
      )) -> df_CT20_1
    
    rbind(df_CT8_1, df_CT20_1)
  } %>% 
  dplyr::select(pathway_name, source, target, prob, timepoint) %>% 
  dplyr::arrange(source, target) %>% 
  dplyr::mutate(timepoint = factor(timepoint, levels = c("CT8", "CT20"))) %>% 
  dplyr::mutate(pathway_name = factor(pathway_name, levels = c("GABA-A"))) %>% 
  ggplot(aes(x = pathway_name, y = prob, fill = timepoint)) + 
  geom_boxplot() + 
  theme_classic() -> p

p + 
  ggpubr::stat_compare_means(
    aes(group = timepoint),
    method = "t.test",
    method.args = list(alternative = "greater"),
    paired = T,
    label = "p.format"
  ) -> p_GABA_A

rbind(
  subsetCommunication(cellchat_CT8, slot.name = "netP", signaling = "GABA-B") %>% dplyr::filter(pval <= 0.01) %>% dplyr::mutate(timepoint = "CT8"),
  subsetCommunication(cellchat_CT20, slot.name = "netP", signaling = "GABA-B") %>% dplyr::filter(pval <= 0.01) %>% dplyr::mutate(timepoint = "CT20")
) %>% 
  {
    df_ = .
    df_ %>% dplyr::filter(timepoint == "CT8") -> df_CT8
    df_ %>% dplyr::filter(timepoint == "CT20") -> df_CT20
    
    full_join(x = df_CT8, y = df_CT20[,c("source", "target")], by = c("source", "target")) %>% 
      dplyr::mutate(pathway_name = "GABA-B", timepoint = "CT8") %>% 
      dplyr::mutate(prob = case_when(
        is.na(prob) ~ 0,
        TRUE ~ prob
      )) -> df_CT8_1
    
    full_join(x = df_CT20, y = df_CT8[,c("source", "target")], by = c("source", "target")) %>% 
      dplyr::mutate(pathway_name = "GABA-B", timepoint = "CT20") %>% 
      dplyr::mutate(prob = case_when(
        is.na(prob) ~ 0,
        TRUE ~ prob
      )) -> df_CT20_1
    
    rbind(df_CT8_1, df_CT20_1)
  } %>% 
  dplyr::select(pathway_name, source, target, prob, timepoint) %>% 
  dplyr::arrange(source, target) %>% 
  dplyr::mutate(timepoint = factor(timepoint, levels = c("CT8", "CT20"))) %>% 
  dplyr::mutate(pathway_name = factor(pathway_name, levels = c("GABA-B"))) %>% 
  ggplot(aes(x = pathway_name, y = prob, fill = timepoint)) + 
  geom_boxplot() + 
  theme_classic() -> p

p + 
  ggpubr::stat_compare_means(
    aes(group = timepoint),
    method = "t.test",
    method.args = list(alternative = "greater"),
    paired = T,
    label = "p.format"
  ) -> p_GABA_B

patchwork::wrap_plots(p_GABA_A, p_GABA_B, ncol = 2, guides = "collect")
  


p$data %>% 
  {
    df_ = .
    df_ %>% dplyr::filter(pathway_name == "GABA-A", timepoint == "CT8") %>% .$prob -> CT8_prob
    df_ %>% dplyr::filter(pathway_name == "GABA-A", timepoint == "CT20") %>% .$prob -> CT20_prob
    t.test(x = CT20_prob, y = CT8_prob, alternative = "greater", paired = F)
  }

p$data %>% 
  {
    df_ = .
    df_ %>% dplyr::filter(pathway_name == "GABA-B", timepoint == "CT8") %>% .$prob -> CT8_prob
    df_ %>% dplyr::filter(pathway_name == "GABA-B", timepoint == "CT20") %>% .$prob -> CT20_prob
    t.test(x = CT20_prob, y = CT8_prob, alternative = "greater", paired = F)
  }

#Other pathways ----
list(
CT8 = cellchat_CT8@netP$pathways,
CT20 = cellchat_CT20@netP$pathways
) %>% ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    list(
      CT8_only = df_ %>% dplyr::filter(CT8 == TRUE, CT20 == FALSE) %>% .$key,
      CT20_only = df_ %>% dplyr::filter(CT20 == TRUE, CT8 == FALSE) %>% .$key,
      Both = df_ %>% dplyr::filter(CT20 == TRUE, CT8 == TRUE) %>% .$key
    )
  } -> pathway_list

list(
  CT8 = cellchat_CT8@net$LRs,
  CT20 = cellchat_CT20@net$LRs
) %>% ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    list(
      CT8_only = df_ %>% dplyr::filter(CT8 == TRUE, CT20 == FALSE) %>% .$key,
      CT20_only = df_ %>% dplyr::filter(CT20 == TRUE, CT8 == FALSE) %>% .$key,
      Both = df_ %>% dplyr::filter(CT20 == TRUE, CT8 == TRUE) %>% .$key
    )
  } -> LR_list

pathway_list %>% 
  purrr::map(function(pathways_){
    pathways_ %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(pathway_){
        draw_heatmap(cellchat = cellchat_CT8, netP = pathway_, pval_ = 0.01) + ggtitle(sprintf("%s (CT-8)", pathway_)) -> p1
        draw_heatmap(cellchat = cellchat_CT20, netP = pathway_, pval_ = 0.01) + ggtitle(sprintf("%s (CT-20)", pathway_)) -> p2        
        patchwork::wrap_plots(p1,p2, ncol = 2, guides = "collect")
      })
  }) -> p_list

pdf("~/Downloads/Cellchat_pathmay_both.pdf", width = 9.94, height = 4.18)
p_list$Both
dev.off()

#VIP ----
subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "VIP")
subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "VIP")
netVisual_heatmap(object = cellchat_CT8, signaling = "VIP", slot.name = "netP") -> p1
netVisual_heatmap(object = cellchat_CT20, signaling = "VIP", slot.name = "netP") -> p2
draw_heatmap(cellchat = cellchat_CT8, netP = "VIP", pval_ = 0.01) + ggtitle("VIP (CT-8)") -> p1
draw_heatmap(cellchat = cellchat_CT20, netP = "VIP", pval_ = 0.01) + ggtitle("VIP (CT-20)") -> p2

list(
  p1,p2
) %>% 
  purrr::map(function(p_){
    p_$data$prob
  }) %>% unlist() %>% max() -> max_prob

i = 0
list(
  p1,p2
) %>% purrr::map(function(p_){
  i <<- i + 1
  p_ + scale_fill_gradient(
      low = "white",
      high = "red",
      limits = c(0, max_prob),
      na.value = "white"
    ) -> p
#  if (i < 3){
#    p + theme(axis.text.x = element_blank()) + xlab(NULL) -> p
#    print(i)
#  }
  return(p)
}) %>% patchwork::wrap_plots(.,ncol = 2, guides = "collect") #-> p_GABA_heatmap

plot_violin(signaling = "VIP", pval_ = 0.01) -> p
c(1,3) %>% 
  purrr::map(function(i){
    p[[i]] -> p
    if (i != 3){
      p +
        theme(
          axis.text.x = element_blank()
        ) + 
        xlab(NULL)
    } -> p
    return(p)
  }) %>% patchwork::wrap_plots(., ncol = 1, guides = "collect") + 
  patchwork::plot_annotation(title = "VIP LR pathway (gene expression)") -> p_violin_VIP
ggsave(plot = p_violin_VIP, filename = "~/Downloads/VIP_LR_pathway_gene_expression.pdf", width = 7.34, height = 4.14)


rbind(
  subsetCommunication(cellchat_CT8, slot.name = "netP", signaling = "VIP") %>% dplyr::filter(pval <= 0.01) %>% dplyr::mutate(timepoint = "CT8"),
  subsetCommunication(cellchat_CT20, slot.name = "netP", signaling = "VIP") %>% dplyr::filter(pval <= 0.01) %>% dplyr::mutate(timepoint = "CT20")
) %>% 
  {
    df_ = .
    df_ %>% dplyr::filter(timepoint == "CT8") -> df_CT8
    df_ %>% dplyr::filter(timepoint == "CT20") -> df_CT20
    
    full_join(x = df_CT8, y = df_CT20[,c("source", "target")], by = c("source", "target")) %>% 
      dplyr::mutate(pathway_name = "VIP", timepoint = "CT8") %>% 
      dplyr::mutate(prob = case_when(
        is.na(prob) ~ 0,
        TRUE ~ prob
      )) -> df_CT8_1
    
    full_join(x = df_CT20, y = df_CT8[,c("source", "target")], by = c("source", "target")) %>% 
      dplyr::mutate(pathway_name = "VIP", timepoint = "CT20") %>% 
      dplyr::mutate(prob = case_when(
        is.na(prob) ~ 0,
        TRUE ~ prob
      )) -> df_CT20_1
    
    rbind(df_CT8_1, df_CT20_1)
  } %>% 
  dplyr::select(pathway_name, source, target, prob, timepoint) %>% 
  dplyr::arrange(source, target) %>% 
  dplyr::mutate(timepoint = factor(timepoint, levels = c("CT8", "CT20"))) %>% 
  dplyr::mutate(pathway_name = factor(pathway_name, levels = c("VIP"))) %>% 
  ggplot(aes(x = pathway_name, y = prob, fill = timepoint)) + 
  geom_boxplot(outliers = TRUE) + 
  theme_classic() -> p

p + 
  ggpubr::stat_compare_means(
    aes(group = timepoint),
    method = "t.test",
    method.args = list(alternative = "greater"),
    paired = T,
    label = "p.format"
  ) -> p_VIP

rbind(
  netAnalysis_contribution(cellchat_CT8, signaling = "VIP", thresh = 0.01)$data %>% dplyr::mutate(timepoint = "CT8"),
  netAnalysis_contribution(cellchat_CT20, signaling = "VIP", thresh = 0.01)$data %>% dplyr::mutate(timepoint = "CT20")
) %>% dplyr::mutate(timepoint = factor(timepoint, levels = c("CT8", "CT20"))) %>% 
  dplyr::arrange(desc(contribution)) %>% 
  head(2) %>% 
  dplyr::mutate(name = sprintf("%s (%s)", name, timepoint)) %>% 
  dplyr::mutate(name = factor(name, levels = rev(name))) %>% 
  ggplot(aes(y = contribution, x = name, fill = timepoint)) + 
  geom_col() + 
  theme_classic() + 
  ylab(NULL) + 
  xlab(NULL) + 
  ggtitle("VIP LRs contribution") + 
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  ) -> p_VIP

##
subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "OPIOID")
draw_heatmap(cellchat = cellchat_CT20, netP = "OPIOID", pval_ = 0.01) + ggtitle("OPIOID (CT-20)")
plot_violin(signaling = "OPIOID", pval_ = 0.01) -> p
c(2,7) %>% 
  purrr::map(function(i){
    p[[i]] -> p
    if (i != 7){
      p +
        theme(
          axis.text.x = element_blank()
        ) + 
        xlab(NULL)
    } -> p
    return(p)
  }) %>% patchwork::wrap_plots(., ncol = 1, guides = "collect") + 
  patchwork::plot_annotation(title = "OPIOID pathway") -> p_violin_OPIOID
netAnalysis_contribution(cellchat_CT20, signaling = "OPIOID", thresh = 0.01)

subsetCommunication(cellchat_CT20, slot.name = "net", signaling = "NRG")
draw_heatmap(cellchat = cellchat_CT20, netP = "NRG", pval_ = 0.01) + ggtitle("NRG (CT-20)")
netAnalysis_contribution(cellchat_CT20, signaling = "NRG", thresh = 0.01)
plot_violin(signaling = "NRG", pval_ = 0.01) -> p
c(1,3) %>% 
  purrr::map(function(i){
    p[[i]] -> p
    if (i != 3){
      p +
        theme(
          axis.text.x = element_blank()
        ) + 
        xlab(NULL)
    } -> p
    return(p)
  }) %>% patchwork::wrap_plots(., ncol = 1, guides = "collect") + 
  patchwork::plot_annotation(title = "NRG pathway") -> p_violin_NRG

draw_heatmap(cellchat = cellchat_CT8, netP = "Netrin", pval_ = 0.01) + ggtitle("Netrin (CT8)")
netAnalysis_contribution(cellchat_CT8, signaling = "Netrin", thresh = 0.01)
plot_violin(signaling = "Netrin", pval_ = 0.01) -> p
c(1,4,5) %>% 
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
  patchwork::plot_annotation(title = "NRG pathway") -> p_violin_Netrin

subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "SEMA5")
subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "CXCL")
subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "Netrin")
subsetCommunication(cellchat_CT8, slot.name = "net", signaling = "GALECTIN")

pathway_list %>% 
#  .[1] %>% 
  purrr::map(function(pathways){
    pathways %>% 
#      .[1] %>% 
      purrr::map(function(pathway_){
        CellChatDB.mouse$interaction %>% dplyr::filter(pathway_name == pathway_) %>% .$annotation %>% unique()
      }) %>% unlist()
  })

subsetCommunication(cellchat_CT8, slot.name = "netP", thresh = 0.05)
subsetCommunication(cellchat_CT20, slot.name = "netP", thresh = 1.1) %>% dplyr::filter(pathway_name == "AVP")
subsetCommunication(cellchat_CT20, slot.name = "net", thresh = 1.1) %>% dplyr::filter(ligand == "Grp")

subsetCommunication(cellchat_CT8, slot.name = "net", thresh = 1.1) %>% dplyr::arrange(desc(pval))

##
rbind(
subsetCommunication(object = cellchat_CT8, slot.name = "netP") %>% dplyr::mutate(timepoint = "CT8"),
subsetCommunication(object = cellchat_CT8, slot.name = "netP") %>% dplyr::mutate(timepoint = "CT20")
) %>% 
  dplyr::filter(pval <= 0.01) %>% 
  dplyr::filter(grepl("GABA", source), target == "Astrocyte") %>% 
  {
    df_ = .
    df_$pathway_name %>% 
      unique() %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(pathway_name_){
        rbind(
          subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT8"),
          subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT20")
        ) %>% dplyr::filter(pathway_name == pathway_name_) -> df_1
      })
  }

list(
  draw_heatmap_net(cellchat = cellchat_CT8, netP = "HBEGF_EGFR", pval_ = 0.01) + ggtitle("HBEGF_EGFR (CT8)"),
  draw_heatmap_net(cellchat = cellchat_CT20, netP = "HBEGF_EGFR", pval_ = 0.01) + ggtitle("HBEGF_EGFR (CT20)")
) %>% patchwork::wrap_plots(., ncol = 2)


#GABA to oligodendrocytes
rbind(
  subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT8"),
  subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT20")
) %>% 
  {
    df_ = .
    df_$interaction_name %>% unique() %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(interaction_name_){
        df_ %>% dplyr::filter(interaction_name == interaction_name_) -> df_
        df_$source %>% unique() -> source
        df_$target %>% unique() -> target
        if ((all(grepl("GABA", source))&(all(target == "Oligodendrocyte")))){
          return(df_)
        }else{
          return(NULL)
        }
      })
  }


rbind(
  subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT8"),
  subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT20")
) %>% dplyr::filter(grepl("GABA", source), target == "Oligodendrocyte") %>% 
  .$interaction_name %>% unique() %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(interaction_name_){
    rbind(
      subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT8"),
      subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT20")
    ) %>% 
      dplyr::filter(interaction_name == interaction_name_)
  })


#Glial to neurons
rbind(
  subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT8"),
  subsetCommunication(object = cellchat_CT20, slot.name = "net") %>% dplyr::mutate(timepoint = "CT20")
) %>% dplyr::filter((source %in% c("Oligodendrocyte", "Astrocyte", "Vascular")), grepl("GABA", target)) %>% 
  .$interaction_name %>% unique() %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(interaction_name_){
    rbind(
      subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT8"),
      subsetCommunication(object = cellchat_CT20, slot.name = "net") %>% dplyr::mutate(timepoint = "CT20")
    ) %>% 
      dplyr::filter(interaction_name == interaction_name_) -> df_
    df_$target %>% grepl("GABA", .) %>% all() -> se_
    (df_$source %in% c("Oligodendrocyte", "Astrocyte", "Vascular")) %>% all() -> se_1
    se_2 = all(c(se_, se_1))
    if (se_2){
      df_ = df_
    }else{
      df_ = NULL
    }
    return(df_)
  }) %>% do.call(rbind, .) -> df_

df_$interaction_name %>% 
  unique() %>% 
  "names<-"(.,.) %>%  
  purrr::map(function(interaction_name_){
    list(
      draw_heatmap_net(cellchat = cellchat_CT8, netP = interaction_name_, pval_ = 0.01) + ggtitle(sprintf("%s (CT8)", interaction_name_)),
      draw_heatmap_net(cellchat = cellchat_CT20, netP = interaction_name_, pval_ = 0.01) + ggtitle(sprintf("%s (CT20)", interaction_name_))
    ) %>% patchwork::wrap_plots(., ncol = 2)    
  }) -> p_list

p_list$SEMA3C_NRP2_PLXNA1
p_list$SEMA3D_NRP2_PLXNA1
p_list$SEMA6D_PLXNA1
p_list$MDK_ALK

#Neuron to non-neuron
rbind(
  subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT8"),
  subsetCommunication(object = cellchat_CT20, slot.name = "net") %>% dplyr::mutate(timepoint = "CT20")
) %>% dplyr::filter((target %in% c("Oligodendrocyte", "Astrocyte", "Vascular")), grepl("GABA", source)) %>% 
  .$interaction_name %>% unique() %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(interaction_name_){
    rbind(
      subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT8"),
      subsetCommunication(object = cellchat_CT20, slot.name = "net") %>% dplyr::mutate(timepoint = "CT20")
    ) %>% 
      dplyr::filter(interaction_name == interaction_name_) -> df_
    df_$source %>% grepl("GABA", .) %>% all() -> se_
    (df_$target %in% c("Oligodendrocyte", "Astrocyte", "Vascular")) %>% all() -> se_1
    se_2 = all(c(se_, se_1))
    if (se_2){
      df_ = df_
    }else{
      df_ = NULL
    }
    return(df_)
  }) %>% do.call(rbind, .) -> df_

df_$interaction_name %>% 
  unique() %>% 
  "names<-"(.,.) %>%  
  purrr::map(function(interaction_name_){
    list(
      draw_heatmap_net(cellchat = cellchat_CT8, netP = interaction_name_, pval_ = 0.01) + ggtitle(sprintf("%s (CT8)", interaction_name_)),
      draw_heatmap_net(cellchat = cellchat_CT20, netP = interaction_name_, pval_ = 0.01) + ggtitle(sprintf("%s (CT20)", interaction_name_))
    ) %>% patchwork::wrap_plots(., ncol = 2)    
  }) -> p_list

df_ %>% dplyr::filter(grepl("GABA", source), target == "Oligodendrocyte")
df_ %>% dplyr::filter(grepl("GABA", source), target == "Vascular")

p_list$HBEGF_EGFR
p_list$NRG1_ERBB4
p_list$PDGFA_PDGFRA
p_list$CDH12_ITGA1_ITGB1

rbind(
  subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT8"),
  subsetCommunication(object = cellchat_CT20, slot.name = "net") %>% dplyr::mutate(timepoint = "CT20")
) %>% dplyr::filter(interaction_name == "HBEGF_EGFR_ERBB2")

rbind(
  subsetCommunication(object = cellchat_CT8, slot.name = "net") %>% dplyr::mutate(timepoint = "CT8"),
  subsetCommunication(object = cellchat_CT20, slot.name = "net") %>% dplyr::mutate(timepoint = "CT20")
) %>% 
  group_by(interaction_name) %>% 
  group_map(function(df_,y){
    if ((nrow(df_) == 2)&(length(unique(df_$timepoint)) == 2)){
      df_
    }else{
      df_ = NULL
    }
  }, .keep = T) %>% do.call(rbind, .) %>% as.data.frame()

draw_heatmap(cellchat = cellchat_CT20, netP = "OPIOID") + ggtitle("OPIOID (CT20)")
netAnalysis_contribution(cellchat_CT20, signaling = "OPIOID", thresh = 0.01)
draw_heatmap(cellchat = cellchat_CT8, netP = "CXCL") + ggtitle("CXCL (CT8)")
netAnalysis_contribution(cellchat_CT8, signaling = "CXCL", thresh = 0.01)

###UMAP ----
library(Seurat)
library(tidyverse)
load("~/Downloads/xenium_obj_filter_GABA_integrated.rda")
xenium_obj_filter_GABA_integrated@meta.data %>% dplyr::filter(seurat_clusters %in% c("0", "1", "3")) %>% rownames() -> cells

DimPlot(xenium_obj_filter_GABA_integrated, group.by = "seurat_clusters", label = TRUE, reduction = "xeniumintegrated_umap", repel = TRUE, cells = cells) + 
  scale_color_manual(
    name = "Cell type",
    values = c(
      "0" = "#F8766D",
      "1" = "#00BA38",
      "3" = "#619CFF"
    ),
    labels = c(
      "0" = "Avp_GABA",
      "1" = "Vip_Nms_C1ql3_GABA",
      "3" = "GRP_GABA"
    )
  ) -> p_UMAP
ggsave(plot = p_UMAP, filename = "~/Downloads/GABA_subtype_umap.pdf", width = 6.58, height = 4.16)

rm(list=ls())
gc()

#Supplementary ----
subsetCommunication(cellchat_CT8)$pval %>% summary()

subsetCommunication(cellchat_CT8, slot.name = "netP") %>% 
  {
    df_ = .
    df_$pathway_name %>% 
      unique() %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(pathway_name_){
        df_ %>% dplyr::filter(pathway_name == pathway_name_) -> df_1
        data.frame(pathway_name = pathway_name_, mean_prob = mean(df_1$prob), sum_prob = sum(df_1$prob))
      }) %>% do.call(rbind, .) %>% 
      dplyr::arrange(desc(sum_prob)) %>% dplyr::mutate(group = "CT8")
  } -> CT8_df

subsetCommunication(cellchat_CT20, slot.name = "netP") %>% 
  {
    df_ = .
    df_$pathway_name %>% 
      unique() %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(pathway_name_){
        df_ %>% dplyr::filter(pathway_name == pathway_name_) -> df_1
        data.frame(pathway_name = pathway_name_, mean_prob = mean(df_1$prob), sum_prob = sum(df_1$prob))
      }) %>% do.call(rbind, .) %>% 
      dplyr::arrange(desc(sum_prob)) %>% dplyr::mutate(group = "CT20")
  } -> CT20_df

intersect(CT8_df$pathway_name, CT20_df$pathway_name) %>% 
  {
    intersected_pathway = .
    rbind(
    CT8_df[intersected_pathway, ],
    CT20_df[intersected_pathway, ]
    ) -> df_
    intersected_pathway %>% 
      purrr::map(function(intersected_pathway_){
        df_ %>% dplyr::filter(pathway_name == intersected_pathway_) -> df_
        data.frame(pathway_name = intersected_pathway_, sum_prob = sum(df_$sum_prob))
      }) %>% do.call(rbind, .) %>% dplyr::arrange(desc(sum_prob)) -> df_1
    df_1$pathway_name[1:20] -> pathways_se
    df_ %>% dplyr::filter(pathway_name %in% pathways_se) -> df_plot
    df_plot$pathway_name = factor(df_plot$pathway_name, levels = rev(pathways_se))
    df_plot$group = factor(df_plot$group, levels = c("CT8", "CT20"))
    df_plot %>% dplyr::mutate(sum_prob = case_when(
      group == "CT8" ~ -sum_prob,
      TRUE ~ sum_prob
    )) -> df_plot
    df_plot %>% 
      ggplot(aes(y = pathway_name, x = sum_prob, fill = group)) + 
      geom_col() + 
      theme_classic() + 
      ylab("LR signaling pathway") + 
      xlab("Community probability") +
      ggtitle("Top 20 LR signaling pathways at CT8 and CT20")
  }

object.list <- list(CT8 = cellchat_CT8,
                    CT20 = cellchat_CT20)

cellchat <- mergeCellChat(object.list,
                          add.names = names(object.list))

rankNet(
  cellchat,
  mode = "comparison",
  stacked = F
) -> p

rankNet(
  cellchat,
  mode = "comparison",
  stacked = T 
) -> p
ggsave(plot = p, filename = "~/Downloads/cellchat_ranknet.pdf", width = 8.33, height = 6.59)


rankNet(
  cellchat,
  mode = "comparison",
  stacked = T,
  do.stat = T,
  return.data = T
) -> df_ 

df_$signaling.contribution

gaba_ct8 <- cellchat_CT8@netP$prob[, , "GABA-A"]
gaba_ct20 <- cellchat_CT20@netP$prob[, , "GABA-A"]

df_CT8 <- as.data.frame(as.table(gaba_ct8)) %>%
  mutate(Group = "CT8")

df_CT20 <- as.data.frame(as.table(gaba_ct20)) %>%
  mutate(Group = "CT20")
df_plot <- bind_rows(df_CT8, df_CT20)

ggplot(df_plot %>% dplyr::filter(Freq != 0),
       aes(Group, Freq, fill = Group)) +
  geom_boxplot(outlier.shape = NA)

wilcox.test(x = df_plot %>% dplyr::filter(Group == "CT8") %>% .$Freq, y = df_plot %>% dplyr::filter(Group == "CT20") %>% .$Freq, paired = T)


gaba_ct8 <- cellchat_CT8@netP$prob[, , "GABA-B"]
gaba_ct20 <- cellchat_CT20@netP$prob[, , "GABA-B"]

df_CT8 <- as.data.frame(as.table(gaba_ct8)) %>%
  mutate(Group = "CT8")

df_CT20 <- as.data.frame(as.table(gaba_ct20)) %>%
  mutate(Group = "CT20")
df_plot <- bind_rows(df_CT8, df_CT20)

ggplot(df_plot %>% dplyr::filter(Freq != 0),
       aes(Group, Freq, fill = Group)) +
  geom_boxplot(outlier.shape = NA)

wilcox.test(x = df_plot %>% dplyr::filter(Group == "CT8") %>% .$Freq, y = df_plot %>% dplyr::filter(Group == "CT20") %>% .$Freq, paired = T)


vip_ct8 <- cellchat_CT8@netP$prob[, , "VIP"]
vip_ct20 <- cellchat_CT20@netP$prob[, , "VIP"]

df_CT8 <- as.data.frame(as.table(vip_ct8)) %>%
  mutate(Group = "CT8")

df_CT20 <- as.data.frame(as.table(vip_ct20)) %>%
  mutate(Group = "CT20")
df_plot <- bind_rows(df_CT8, df_CT20)

wilcox.test(x = df_plot %>% dplyr::filter(Group == "CT8") %>% .$Freq, y = df_plot %>% dplyr::filter(Group == "CT20") %>% .$Freq, paired = T)

ggplot(df_plot %>% dplyr::filter(Freq != 0),
       aes(Group, Freq, fill = Group)) +
  geom_boxplot(outlier.shape = NA)
