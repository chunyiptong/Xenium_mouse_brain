#devtools::install_github('YMa-lab/spMOCA')
library(spMOCA)
library(Seurat)
library(tidyverse)
library(WGCNA)
library(SPARK)

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
#scn_core_cells = select_cells(region__ = "SCN_core")
#scn_shell_cells = select_cells(region__ = "SCN_shell")

xenium_obj_filtered = subset(xenium_obj_filtered, cells = scn_cells)
GetAssayData(xenium_obj_filtered, assay = "Xenium", layer = "counts") -> counts_matrix
rowSums(counts_matrix) %>% {.[.!=0]} %>% names() -> gene_keep
counts_matrix = counts_matrix %>% as.data.frame()
counts_matrix = counts_matrix[gene_keep, ]
dim(counts_matrix)

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    scn_cells %>% {.[grepl(sample_, .)]} -> cells_
    counts_matrix[,cells_]
  }) -> count_matrix_per_sample

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    scn_cells %>% {.[grepl(sample_, .)]} -> cells_
    xy_coordinates[[sample_]] %>% column_to_rownames("cell") %>% .[cells_, ]
  }) -> xy_per_sample

load("~/Downloads/df_spatiotemporal_all_2.rda")

df_spatiotemporal_all_2 %>% 
  {
    df_ = .
    df_ %>% dplyr::filter(`CT-8_cauchy_padj_sparkx` < 0.05) %>% .$gene -> CT8_SVG
    df_ %>% dplyr::filter(`CT-20_cauchy_padj_sparkx` < 0.05) %>% .$gene -> CT20_SVG
    list(CT8 = CT8_SVG,
         CT20 = CT20_SVG)
  } -> spark_gene_list

moduleDetectionspMOCA_custom <- function(
    object,
    module.method = "WGCNA",
#    minClusterSize = 30,
    minClusterSize = 10,
    deepSplit = 4,
    pamRespectsDendro = FALSE,
    cutHeight = NULL,
#    power = 1
    power = 4 
) {
  REst <- object@corr_est
  
  if (module.method == "WGCNA") {
    cat("## Performing WGCNA ...\n")
    
    corr.mat <- round(REst, 6)
    corr.mat[lower.tri(corr.mat)] <- t(corr.mat)[lower.tri(corr.mat)]
    
    adj.mat <- WGCNA::adjacency.fromSimilarity(
      corr.mat,
      power = power,
      type = "signed"
    )
    
    TOM <- WGCNA::TOMsimilarity(adj.mat)
    dissTOM <- 1 - TOM
    
    rownames(dissTOM) <- colnames(dissTOM) <- rownames(corr.mat)
    
    hclust_dist <- hclust(as.dist(dissTOM), method = "average")
    
    wcgna_label <- dynamicTreeCut::cutreeDynamic(
      dendro = hclust_dist,
      distM = dissTOM,
      deepSplit = deepSplit,
      pamRespectsDendro = pamRespectsDendro,
      minClusterSize = minClusterSize,
      method = "hybrid",
      cutHeight = cutHeight
    )
    
    gene.module.all <- WGCNA::intramodularConnectivity(
      (corr.mat + 1) / 2,
      colors = wcgna_label
    )
    
    gene.module.all$k <- wcgna_label
    gene.module.all$gene <- object@feature_select
    
    object@network <- gene.module.all
  }
  
  return(object)
}

moduleDetectionspMOCA_custom_1 <- function(
    object,
    module.method = "WGCNA",
    minClusterSize = 15,
    deepSplit = 4,
    cutHeight = NULL,
    power = 1,
    networkType = "signed",
    pamRespectsDendro = FALSE
) {
  REst <- object@corr_est
  
  if (module.method == "WGCNA") {
    cat("## Performing WGCNA ...\n")
    
    corr.mat <- round(REst, 6)
    corr.mat[lower.tri(corr.mat)] <- t(corr.mat)[lower.tri(corr.mat)]
    
    adj.mat <- WGCNA::adjacency.fromSimilarity(
      corr.mat,
      power = power,
      type = networkType
    )
    
    TOM <- WGCNA::TOMsimilarity(adj.mat)
    dissTOM <- 1 - TOM
    
    rownames(dissTOM) <- colnames(dissTOM) <- rownames(corr.mat)
    
    hclust_dist <- hclust(as.dist(dissTOM), method = "average")
    
    wcgna_label <- dynamicTreeCut::cutreeDynamic(
      dendro = hclust_dist,
      distM = dissTOM,
      method = "hybrid",
      deepSplit = deepSplit,
      pamRespectsDendro = pamRespectsDendro,
      minClusterSize = minClusterSize,
      cutHeight = cutHeight
    )
    
    gene.module.all <- WGCNA::intramodularConnectivity(
      (corr.mat + 1) / 2,
      colors = wcgna_label
    )
    
    gene.module.all$k <- wcgna_label
    gene.module.all$gene <- object@feature_select
    
    object@network <- gene.module.all
  }
  
  return(object)
}

run_spMOCA = function(count_data, location, all_gene, gene_spark){
  ##### Create spMOCA object
  object = createSPMOCAObject(spatial_countMat = count_data,
                              spatial_location = location,
                              feature_list = all_gene)
  
  #### Normalization
  object = normalizeSPMOCA(object)
  
  #### Feature Selection
  object = featureSelectionSPMOCA(object,
                                  feature.selection = NULL,
                                  feature.list = gene_spark)

  
  object = createSpatialKernel(object)
  object = estimatespMOCA(object)
#  object = moduleDetectionspMOCA(object)
  object = moduleDetectionspMOCA_custom(object)
#  object = moduleDetectionspMOCA_custom_1(object)
  
  ms_prop = calModuleScore(object)
  res_list = list(object = object, ms_prop = ms_prop)
  return(res_list)
}

all_genes = counts_matrix %>% rownames()
CT8_1_spMOCA_res = run_spMOCA(count_data = count_matrix_per_sample$`CT-8-1`, location = xy_per_sample$`CT-8-1`, all_gene = all_genes, gene_spark = spark_gene_list$CT8)
CT8_1_spMOCA_res$object@network$k %>% table()

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    timepoint_ = gsub("(CT-\\d+?)-.+", "\\1", sample_) %>% gsub("-", "", .)
    spMOCA_res = run_spMOCA(count_data = count_matrix_per_sample[[sample_]], location = xy_per_sample[[sample_]], all_gene = all_genes, gene_spark = spark_gene_list[[timepoint_]])
  }) -> spMOCA_res_list
save(spMOCA_res_list, file = "~/Downloads/spMOCA_res_list.rda")

#spMOCA_res_list %>% 
##  .[1] %>% 
#  purrr::map(function(res_){
#    res_$object@network$k %>% unique() %>% sort() %>% 
#      "names<-"(., sprintf("MOD%s", .)) %>% 
#      purrr::map(function(mod_){
#        res_$object@network %>% dplyr::filter(k == mod_) %>% 
#          dplyr::arrange(desc(kWithin)) %>% dplyr::mutate(idx = 1:nrow(.)) -> df_
#        df_ %>% dplyr::filter(gene %in% c("Vip", "Pcsk1", "Lhx1"))
#      })
#  })

#spMOCA_res_list %>% 
#  #  .[1] %>% 
#  purrr::map(function(res_){
#    res_$object@network$k %>% table()
#  })

load("~/Downloads/DEG_dfs.rda")

c(1:2) %>% 
  "names<-"(., sprintf("MOD%s", .)) %>% 
  purrr::map(function(mod_){
    spMOCA_res_list %>% 
      .[1:3] %>% 
      purrr::map(function(res_){
        res_$object@network -> df_
        df_ %>% 
          dplyr::filter(k == mod_) -> df_
        df_ %>% 
          dplyr::select(gene, kWithin) -> df_
      }) %>% purrr::reduce(., inner_join, by = "gene") -> df_
    df_ %>% column_to_rownames("gene") %>% 
      apply(., 1, mean) %>% unname() -> mean_KWithin
    data.frame(gene = df_$gene, mean_KWithin = mean_KWithin) %>% 
      dplyr::arrange(desc(mean_KWithin)) %>% 
      dplyr::filter(gene %in% DEG_df_SCN$gene)
  }) -> CT8_MOD

c(1:2) %>% 
  "names<-"(., sprintf("MOD%s", .)) %>% 
  purrr::map(function(mod_){
    spMOCA_res_list %>% 
      .[4:6] %>% 
      purrr::map(function(res_){
        res_$object@network -> df_
        df_ %>% 
          dplyr::filter(k == mod_) -> df_
        df_ %>% 
          dplyr::select(gene, kWithin) -> df_
      }) %>% purrr::reduce(., inner_join, by = "gene") -> df_
    df_ %>% column_to_rownames("gene") %>% 
      apply(., 1, mean) %>% unname() -> mean_KWithin
    data.frame(gene = df_$gene, mean_KWithin = mean_KWithin) %>% 
      dplyr::arrange(desc(mean_KWithin)) %>% 
      dplyr::filter(gene %in% DEG_df_SCN$gene)
  }) -> CT20_MOD

load("~/Downloads/New_SVG.rda")
CT8_MOD %>% 
  purrr::map(function(genes_){
    genes_ %>% dplyr::filter(gene %in% CT8_SVG_1) %>% dplyr::arrange(desc(mean_KWithin))
  }) -> CT8_MOD

CT20_MOD %>% 
  purrr::map(function(genes_){
    genes_ %>% dplyr::filter(gene %in% CT20_SVG_1) %>% dplyr::arrange(desc(mean_KWithin))
  }) -> CT20_MOD

CT8_MOD %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, MOD_){
    round(nrow(df_)*0.05, 0) -> row_
    df_[1:row_, ] %>% dplyr::mutate(MOD = sprintf("CT8_%s", MOD_))
  }) %>% do.call(rbind, .) %>% 
  dplyr::mutate(
    gene_id = row_number(),
    angle = 90 - 360 * (gene_id - 0.5) / n(),
    hjust = ifelse(angle < -90, 1, 0),
    angle_label = ifelse(angle < -90, angle + 180, angle)
  ) %>% 
  dplyr::mutate(gene = factor(gene, levels = gene)) %>% 
  ggplot(aes(x = gene, y = 1, color = MOD)) + 
  geom_point(size = 2) +
  geom_text(
    aes(
      y = 1.15,
      label = gene,
      angle = angle_label,
      hjust = hjust
    ),
    size = 3,
    show.legend = FALSE
  ) +
  coord_polar(start = 0) +
  ylim(0, 1.4) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) + 
  scale_color_manual(values = c("CT8_MOD1" = "#f66a62", "CT8_MOD2" = "#02b7bc")) + 
  ggtitle("CT8 top 5% module hub genes") -> CT8_top5_percent_hub_genes

CT8_top5_percent_hub_genes + 
  scale_color_manual(
    values = c("CT8_MOD1" = "#f66a62", "CT8_MOD2" = "#02b7bc"), 
    labels = c("CT8_MOD1" = "CT8 Neuronal module", "CT8_MOD2" = "CT8 Astrocyte module")
  ) -> CT8_top5_percent_hub_genes


CT20_MOD %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, MOD_){
    round(nrow(df_)*0.05, 0) -> row_
    df_[1:row_, ] %>% dplyr::mutate(MOD = sprintf("CT20_%s", MOD_))
  }) %>% do.call(rbind, .) %>% 
  dplyr::mutate(
    gene_id = row_number(),
    angle = 90 - 360 * (gene_id - 0.5) / n(),
    hjust = ifelse(angle < -90, 1, 0),
    angle_label = ifelse(angle < -90, angle + 180, angle)
  ) %>% 
  dplyr::mutate(gene = factor(gene, levels = gene)) %>% 
  ggplot(aes(x = gene, y = 1, color = MOD)) + 
  geom_point(size = 2) + 
  geom_text(
    aes(
      y = 1.15,
      label = gene,
      angle = angle_label,
      hjust = hjust
    ),
    size = 3,
    show.legend = FALSE
  ) +
  coord_polar(start = 0) +
  ylim(0, 1.4) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) + 
  scale_color_manual(values = c("CT20_MOD1" = "#e39500", "CT20_MOD2" = "#c071ff")) + 
  ggtitle("CT20 top 5% module hub genes") -> CT20_top5_percent_hub_genes

CT20_top5_percent_hub_genes + 
  scale_color_manual(
    values = c("CT20_MOD1" = "#e39500", "CT20_MOD2" = "#c071ff"), 
    labels = c("CT20_MOD1" = "CT20 Neuronal module", "CT20_MOD2" = "CT20 Astrocyte module")
  ) -> CT20_top5_percent_hub_genes


patchwork::wrap_plots(CT8_top5_percent_hub_genes, CT20_top5_percent_hub_genes, ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom") &
  guides(color = guide_legend(nrow = 1)) -> p_top_5_percent_hub_genes

ggsave(plot = p_top_5_percent_hub_genes, filename = "~/Downloads/p_top_5_percent_hub_genes_V2.pdf", width = 9.09, height = 5.96) #Figure S9A

#CT8_MOD %>% 
#  purrr::map(function(df_){
#    df_ %>% dplyr::filter(gene %in% c("Vip", "Pcsk1"))
#  })
#
#CT20_MOD %>% 
#  purrr::map(function(df_){
#    df_ %>% dplyr::filter(gene %in% c("Vip", "Pcsk1"))
#  })

list(
  CT8 = CT8_MOD,
  CT20 = CT20_MOD
) %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(MOD, CT_){
    MOD %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(df_, MOD_){
        data.frame(n_genes = nrow(df_), CT = CT_, MOD = sprintf("%s_%s", CT_, MOD_))
      }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) %>% 
  dplyr::mutate(CT = factor(CT, levels = c("CT8", "CT20"))) %>% 
  ggplot(aes(x = CT, y = n_genes, fill = MOD, group = MOD)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  ylab("Number of genes") + xlab("Timepoint") + 
  theme_classic() + 
  scale_fill_manual(values = c("CT8_MOD1" = "#f66a62", "CT8_MOD2" = "#02b7bc", "CT20_MOD1" = "#e39500", "CT20_MOD2" = "#c071ff")) + 
  geom_text(
    aes(label = n_genes),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 4
  ) -> p_barplot
p_barplot + scale_fill_manual(values = c("CT8_MOD1" = "#f66a62", "CT8_MOD2" = "#02b7bc", "CT20_MOD1" = "#e39500", "CT20_MOD2" = "#c071ff"), 
                              labels = c("CT8_MOD1" = "CT8 Neuronal module", "CT8_MOD2" = "CT8 Astrocyte module", "CT20_MOD1" = "CT20 Neuronal module", "CT20_MOD2" = "CT20 Astrocyte module")) -> p_barplot_1
ggsave(plot = p_barplot_1, filename = "~/Downloads/p_barplot_1_V2.pdf", width = 8.08, height = 4.74) #Figure S9B

ggvenn_list = list()
list(
  CT8 = CT8_MOD,
  CT20 = CT20_MOD
) %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(MOD, CT_){
    MOD %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(df_, MOD_){
        data.frame(n_genes = nrow(df_), CT = CT_, MOD = sprintf("%s_%s", CT_, MOD_))
        ggvenn_list[[sprintf("%s_%s", CT_, MOD_)]] <<- df_$gene
      }) 
  })
ggvenn_list %>% 
  ggvenn::ggvenn(fill_color = c("#f66a62", "#02b7bc", "#e39500", "#c071ff"), fill_alpha = 0.7) -> p_ggvenn

patchwork::wrap_plots(p_barplot, p_ggvenn, ncol = 1)

patchwork::wrap_plots(
ggvenn_list[c("CT8_MOD1", "CT20_MOD1")] %>% "names<-"(., c("CT8 Neuronal module", "CT20 Neuranal module")) %>% ggvenn::ggvenn(fill_color = c("#f66a62", "#e39500"), fill_alpha = 0.7, text_size = 6),
ggvenn_list[c("CT8_MOD2", "CT20_MOD2")] %>% "names<-"(., c("CT8 Glial module", "CT20 Glial module")) %>% ggvenn::ggvenn(fill_color = c("#02b7bc", "#c071ff"), fill_alpha = 0.7, text_size = 6), 
ncol = 1) -> p_venn #Figure S9C


spMOCA_res_list$`CT-8-1`$ms_prop
spMOCA_res_list$`CT-8-1`$object@spatial_location
spMOCA_res_list %>% 
#    .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(res_, sample_){
    res_$ms_prop %>% 
      as.data.frame() -> df_
#    xy_per_sample[[sample_]][rownames(df_), ] -> df_1
    res_$object@spatial_location -> df_1
    cbind(df_, df_1) -> df_2
    df_2 %>% 
      ggplot(aes(x = x, y = y, color = Mod1)) + 
      geom_point(alpha = 0.7) + 
      scale_color_gradient(low = "grey", high = "red", name = "Neuronal\nodule score") + 
      theme_void() + 
      ggtitle(sample_)-> p_mod1
    df_2 %>% 
      ggplot(aes(x = x, y = y, color = Mod2)) + 
      geom_point(alpha = 0.7) + 
      scale_color_gradient(low = "grey", high = "red", name = "Glial\nmoduel score") + 
      theme_void() + 
      ggtitle(sample_)-> p_mod2
    patchwork::wrap_plots(p_mod1, p_mod2, ncol = 1)
  }) %>% patchwork::wrap_plots(., ncol = 6) -> p_module_score

ggvenn_list %>% 
  ggvenn::list_to_data_frame() %>% {
    df_ = .
    list(
    CT8_MOD1_only = df_ %>% dplyr::filter(CT8_MOD1 == TRUE, CT8_MOD2 == FALSE, CT20_MOD1 == FALSE, CT20_MOD2 == FALSE) %>% .$key,
    CT8_MOD2_only = df_ %>% dplyr::filter(CT8_MOD1 == FALSE, CT8_MOD2 == TRUE, CT20_MOD1 == FALSE, CT20_MOD2 == FALSE) %>% .$key,
    CT20_MOD1_only = df_ %>% dplyr::filter(CT8_MOD1 == FALSE, CT8_MOD2 == FALSE, CT20_MOD1 == TRUE, CT20_MOD2 == FALSE) %>% .$key,
    CT20_MOD2_only = df_ %>% dplyr::filter(CT8_MOD1 == FALSE, CT8_MOD2 == FALSE, CT20_MOD1 == FALSE, CT20_MOD2 == TRUE) %>% .$key,
    MOD1_common = df_ %>% dplyr::filter(CT8_MOD1 == TRUE, CT8_MOD2 == FALSE, CT20_MOD1 == TRUE, CT20_MOD2 == FALSE) %>% .$key,
    MOD2_common = df_ %>% dplyr::filter(CT8_MOD1 == FALSE, CT8_MOD2 == TRUE, CT20_MOD1 == FALSE, CT20_MOD2 == TRUE) %>% .$key
    )
  } -> gene_list

gene_list %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(genes_, group_){
    data.frame(idx = 1:length(genes_)) -> df_
    df_[[group_]] <- genes_
    return(df_)
  }) %>% purrr::reduce(., full_join, by = "idx") %>% 
  column_to_rownames("idx") %>% 
  write_csv(file = "~/Dropbox/Spatial_transcriptomics_new/rda/SCN_coexpressed_genes.csv")

#2. GO Analysis ----
library(clusterProfiler)
library(org.Mm.eg.db)

read.table("~/Dropbox/Spatial_transcriptomics_new/rda/GRCm39_mart_export.txt", header = T, sep = "\t", stringsAsFactors = F) -> gene.ref
gene.ref %>% dplyr::filter(Gene.name != "") %>% .[,c(1,3)] %>% 
  distinct() -> gene.ref
readRDS("~/Dropbox/singulomics/github_rda/gene.ref.rds") -> gene.ref.old

gene_list %>% 
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
}) -> res_go_analysis

res_go_analysis %>% 
  purrr::map(function(list_){
    list_ %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(res_, test_){
        if (test_ == "kk"){
          res_ = res_
        }else{
          print("run simplify")
          res_ <- clusterProfiler::simplify(
            res_,
            cutoff = 0.5,
            by = "p.adjust",
            select_fun = min
          )                  
        }
        return(res_)
      })
  }) -> res_go_analysis_simplify

save(res_go_analysis, res_go_analysis_simplify, file = "~/Downloads/coexpressed_genes_GO_analysis.rda")

res_go_analysis_simplify %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(list_, group_){
    list_ %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(res_, test_){
        res_@result %>% dplyr::filter(p.adjust < 0.05) %>% head(5) %>% 
          dplyr::select(Description, p.adjust) %>% 
          dplyr::mutate(method = test_)
      }) %>% do.call(rbind, .) -> df_plot
    df_plot %>% 
      dplyr::mutate(Description = factor(Description, levels = rev(Description)), method = factor(method, levels = c("ego_CC", "ego_BP", "ego_MF", "kk"))) -> df_plot
    df_plot %>% 
      ggplot(aes(x = -log10(p.adjust), y = Description, fill = method)) + 
      geom_col() + 
      theme_classic() + 
      ggtitle(group_)
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 2, guides = "collect")

#MOD1_only
#Circaidan
res_go_analysis_simplify$CT8_MOD1_only$kk@result %>% dplyr::filter(Description == "Circadian rhythm") %>% 
  {
    df_ = .
    entrez_id_ = strsplit(df_$geneID, "/") %>% unlist()
    AnnotationDbi::select(org.Mm.eg.db, keys = entrez_id_, 
                          columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                          keytype = "ENTREZID")    
  } %>% .$SYMBOL
#[1] "Fbxw11" "Clock"  "Csnk1e" "Prkaa2" "Per1"  

res_go_analysis_simplify$CT20_MOD1_only$kk@result %>% dplyr::filter(Description == "Circadian entrainment") %>% 
  {
    df_ = .
    entrez_id_ = strsplit(df_$geneID, "/") %>% unlist()
    AnnotationDbi::select(org.Mm.eg.db, keys = entrez_id_, 
                          columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                          keytype = "ENTREZID")    
  } %>% .$SYMBOL
#[1] "Fos"    "Mtnr1a" "Kcnj5"  "Per2"   "Prkg2"  "Ryr2"   "Ryr3"   "Adcy1"  "Mapk3"  "Grin2b" "Gnb1"   "Gria4" 

#MOD1_common
res_go_analysis_simplify$MOD1_common$kk@result %>% dplyr::filter(Description == "GABAergic synapse") %>% 
  {
    df_ = .
    entrez_id_ = strsplit(df_$geneID, "/") %>% unlist()
    AnnotationDbi::select(org.Mm.eg.db, keys = entrez_id_, 
                          columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                          keytype = "ENTREZID")    
  } %>% .$SYMBOL


#MOD2 
res_go_analysis_simplify$MOD2_common$ego_BP@result %>% dplyr::filter(Description == "glial cell differentiation") %>% 
  .$geneID %>% str_split("/") %>% unlist()

res_go_analysis_simplify$CT8_MOD2_only$ego_BP@result %>% dplyr::filter(Description == "glial cell differentiation") %>% 
  .$geneID %>% str_split("/") %>% unlist()

res_go_analysis_simplify$CT20_MOD2_only$ego_BP@result %>% dplyr::filter(Description == "neural precursor cell proliferation") %>% 
  .$geneID %>% str_split("/") %>% unlist()

res_go_analysis_simplify$CT20_MOD2_only$ego_BP@result %>% dplyr::filter(Description == "cell fate commitment") %>% 
  .$geneID %>% str_split("/") %>% unlist()

#Permutation ----
#load spMOC_res_list, CT20_MOD, CT8MOD
read_csv("~/Downloads/new_spatiotemporal_genes_df.csv") -> new_spatiotemporal_genes_df

CT8_MOD$MOD1 %>% dplyr::filter(gene %in% c("Lhx1", "Vip", "Pcsk1"))
CT20_MOD$MOD1 %>% dplyr::filter(gene %in% c("Lhx1", "Vip", "Pcsk1"))

permutation_test = function(MOD_df, genes, n_permutation){
  MOD_df %>% dplyr::filter(gene %in% genes) -> se_gene_df
  print(se_gene_df)
  se_gene_df$mean_KWithin -> Kwithin_se
  1:n_permutation %>% 
    purrr::map(function(i){
      set.seed(i)
      sample(MOD_df$mean_KWithin, length(genes), replace = F) -> permuted_Kwithin
      mean(permuted_Kwithin) >  mean(Kwithin_se)
    }) %>% unlist() -> test_result
  sum(test_result)/n_permutation
}

permutation_test_1 = function(MOD_df, genes, n_permutation){
  MOD_df = MOD_df %>% dplyr::mutate(rank = 1:nrow(.))
  MOD_df %>% dplyr::filter(gene %in% genes) -> se_gene_df
  print(se_gene_df)
  se_gene_df$rank -> rank_se
  rank_se
  1:n_permutation %>% 
    purrr::map(function(i){
      set.seed(i)
      sample(MOD_df$rank, length(genes), replace = F) -> permuted_rank
      mean(permuted_rank) <  mean(rank_se)
    }) %>% unlist() -> test_result
  sum(test_result)/n_permutation
}

permutation_test(MOD_df = CT8_MOD$MOD1, genes = c("Lhx1", "Vip", "Pcsk1"), n_permutation = 10000)  #p=0.0245
permutation_test(MOD_df = CT20_MOD$MOD1, genes = c("Lhx1", "Vip", "Pcsk1"), n_permutation = 10000) #p=0.0525

permutation_test(MOD_df = CT8_MOD$MOD1, genes = c("Lhx1", "Vip"), n_permutation = 10000)  #p=0.0127
permutation_test(MOD_df = CT20_MOD$MOD1, genes = c("Lhx1", "Vip"), n_permutation = 10000) #p=0.016

permutation_test(MOD_df = CT8_MOD$MOD1, genes = c("Vip", "Pcsk1"), n_permutation = 10000)  #p=0.0687
permutation_test(MOD_df = CT20_MOD$MOD1, genes = c("Vip", "Pcsk1"), n_permutation = 10000) #p=0.2078

permutation_test(MOD_df = CT8_MOD$MOD1, genes = c("Vip", "Pcsk1"), n_permutation = 10000)  #p=0.0687
permutation_test(MOD_df = CT20_MOD$MOD1, genes = c("Vip", "Pcsk1"), n_permutation = 10000) #p=0.2078

permutation_test(MOD_df = CT8_MOD$MOD1, genes = c("Per2", "Csnk1e"), n_permutation = 10000) #p=0.2078
permutation_test(MOD_df = CT20_MOD$MOD1, genes = c("Per1", "Csnk1e"), n_permutation = 10000) #p=0.2078

CT8_MOD$MOD1$mean_KWithin %>% density() %>% plot()
CT8_MOD$MOD1$mean_KWithin %>% summary()
CT8_MOD$MOD1 %>% dplyr::filter(gene %in% c("Per1", "Csnk1e", "Per2"))

CT8_MOD$MOD1 %>% dplyr::filter(gene %in% c("Lhx1", "Vip", "Pcsk1"))
CT8_MOD$MOD1$mean_KWithin %>% summary()
CT20_MOD$MOD1 %>% dplyr::filter(gene %in% c("Lhx1", "Vip", "Pcsk1"))
CT20_MOD$MOD1$mean_KWithin %>% summary()


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
permutaion_test_2 = function(genes, n_permutation, mod){
  spMOCA_res_list %>% 
    purrr::map2(.x=.,.y=names(.),.f=function(res_, sample_){
      res_$object@network %>% dplyr::filter(k == mod) -> df_
      res_$object@network[genes, ] -> res_df
      res_df$kWithin %>% mean() -> mean_k_within
      1:n_permutation %>% 
        purrr::map(function(i){
          set.seed(i)
          sample(df_$kWithin, length(genes), replace = F) %>% mean() -> permutated_k_within
        }) %>% unlist() -> permutated_k_within
      sum(permutated_k_within > mean_k_within)/n_permutation
    }) %>% as.data.frame() %>% 
    pivot_longer(cols = everything(), names_to = "sample", values_to = "pval") %>% 
    dplyr::mutate(CT = c(rep("CT8",3), rep("CT20", 3))) -> df_1
  
  df_1 %>% 
    group_by(CT) %>% 
    summarise(cauchy_pval = cauchy.test(pval), mean_pval = mean(pval)) %>% 
    ungroup()
}

spMOCA_res_list %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(res_, sample_){
    res_$object@network %>% dplyr::filter(gene %in% c("Per1", "Csnk1e", "Per2") ) -> df_
  })



permutaion_test_2(genes = c("Lhx1", "Vip", "Pcsk1"), n_permutation = 10000, mod = 1)
permutaion_test_2(genes = c("Lhx1", "Vip"), n_permutation = 10000, mod = 1)
permutaion_test_2(genes = c("Vip", "Pcsk1"), n_permutation = 10000, mod = 1)

permutaion_test_2(genes = c("Per2", "Csnk1e"), n_permutation = 10000, mod = 1)


#gene x gene co-expression plot ----
library(tidyverse)
library(Seurat)

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
scn_core_cells = select_cells(region__ = "SCN_core")
scn_shell_cells = select_cells(region__ = "SCN_shell")

cell_keep = c(scn_cells, scn_core_cells, scn_shell_cells) %>% unique()
xenium_obj_filtered = subset(xenium_obj_filtered, cells = cell_keep)

GetAssayData(xenium_obj_filtered, assay = "Xenium", layer = "counts") -> counts_matrix

counts_matrix %>% colSums() %>% median() -> median_
counts_matrix %>% as.data.frame() %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*median_})) -> norm_counts_matrix
#xenium_obj_filtered@meta.data[rownames(df_),] %>% .$predicted.id -> df_$celltype

library(GGally)
library(hexbin)

hex_plot <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_hex(bins = 30) +
    scale_fill_viridis_c(
      trans = "log10",
      name = "Count"
    ) +
    theme_classic()
}

cor_p_panel <- function(data, mapping, method = "spearman", ...) {
  
  x <- GGally::eval_data_col(data, mapping$x)
  
  y <- GGally::eval_data_col(data, mapping$y)
  
  keep <- complete.cases(x, y)
  
  test <- cor.test(
    x[keep],
    y[keep],
    method = method,
    exact = FALSE
  )
  
  r_text <- if (method == "spearman") {
    paste0("\u03c1 = ", sprintf("%.2f", unname(test$estimate)))
  } else {
    paste0("r = ", sprintf("%.2f", unname(test$estimate)))
  }
  
  p_text <- if (test$p.value < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", sprintf("%.3f", test$p.value))
  }
  
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.58,
      label = r_text,
      size = 5
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 0.40,
      label = p_text,
      size = 4
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

norm_counts_matrix[c("Lhx1" ,"Vip", "Pcsk1"), scn_core_cells] %>% 
  t() %>% 
  as.data.frame() -> data

ggpairs(
  data,
  columns = c("Lhx1", "Vip", "Pcsk1"),
  lower = list(continuous = hex_plot),
  upper = list(continuous = 
                 wrap(cor_p_panel, method = "spearman")
                      ),
  diag = list(continuous = "densityDiag")
) -> p

ggsave(plot = p, filename = "~/Downloads/pairwise_correlation.pdf", width = 4.98, height = 4.03)

#Wen's scRNA-seq
list.files("~/Dropbox/Spatial_transcriptomics_new/external_data/Wen_mouse_SCN_circadian/GSE117295_RAW", full.names = T) %>% 
  .[grep("CT", .)] %>% 
  "names<-"(., gsub("/.+/.+(CT\\d+)\\.csv\\.gz", "\\1", .)) %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(file_, CT_){
    print(file_)
    read.csv(file = file_, row.names = 1, header = T) -> df_
  }) -> mat_list

mat_list$CT14[1:5,1:5]
mat_list %>% 
  purrr::map(function(df_){
    df_ %>% rownames_to_column("gene")
  }) %>% purrr::reduce(., inner_join, by = "gene") -> mat
mat %>% column_to_rownames("gene") -> mat
c("Lhx1", "Vip", "Pcsk1") %in% rownames(mat)

colSums(mat) %>% median() -> median_
mat %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*median_})) -> mat_norm

mat_norm[c("Lhx1", "Vip", "Pcsk1"), ] %>% 
  t() %>% 
  as.data.frame() -> data_scRNA

ggpairs(
  data_scRNA,
  columns = c("Lhx1", "Vip", "Pcsk1"),
  lower = list(continuous = hex_plot),
  upper = list(continuous = 
                 wrap(cor_p_panel, method = "spearman")
  ),
  diag = list(continuous = "densityDiag")
) -> p_scRNA

ggsave(plot = p_scRNA, filename = "~/Downloads/pairwise_correlation_scRNA.pdf", width = 4.98, height = 4.03)

#New analysis ----
library(tidyverse)
library(ggpubr)

load("~/Downloads/spMOCA_res_list.rda")
load("~/Downloads/DEG_dfs.rda")
load("~/Downloads/df_combinedPval_all_samples_sparkx.Rdata") ##2.2 SVG results ----

c(1:2) %>% 
  "names<-"(., sprintf("MOD%s", .)) %>% 
  purrr::map(function(mod_){
    spMOCA_res_list %>% 
      .[1:3] %>% 
      purrr::map(function(res_){
        res_$object@network -> df_
        df_ %>% 
          dplyr::filter(k == mod_) -> df_
        df_ %>% 
          dplyr::select(gene, kWithin) -> df_
      }) %>% purrr::reduce(., inner_join, by = "gene") -> df_
    df_ %>% column_to_rownames("gene") %>% 
      apply(., 1, mean) %>% unname() -> mean_KWithin
    data.frame(gene = df_$gene, mean_KWithin = mean_KWithin) %>% 
      dplyr::arrange(desc(mean_KWithin))
  }) -> CT8_MOD

c(1:2) %>% 
  "names<-"(., sprintf("MOD%s", .)) %>% 
  purrr::map(function(mod_){
    spMOCA_res_list %>% 
      .[4:6] %>% 
      purrr::map(function(res_){
        res_$object@network -> df_
        df_ %>% 
          dplyr::filter(k == mod_) -> df_
        df_ %>% 
          dplyr::select(gene, kWithin) -> df_
      }) %>% purrr::reduce(., inner_join, by = "gene") -> df_
    df_ %>% column_to_rownames("gene") %>% 
      apply(., 1, mean) %>% unname() -> mean_KWithin
    data.frame(gene = df_$gene, mean_KWithin = mean_KWithin) %>% 
      dplyr::arrange(desc(mean_KWithin))
  }) -> CT20_MOD

#counts_matrix %>% rowSums() %>% {.[.!=0]} %>% names() -> gene_selected
DEG_df_SCN$gene -> gene_selected
df_combinedPval_all_samples_sparkx %>% 
  {
    df_ = .
    df_ %>% dplyr::select(gene, `CT-20_cauchy_pval`, `CT-8_cauchy_pval`) -> df_
    setdiff(gene_selected, df_$gene) -> genes_
    tibble(gene = genes_, `CT-20_cauchy_pval` = 1, `CT-8_cauchy_pval` = 1) -> df_1
    #    df_1 %>% head()
    rbind(df_, df_1) -> df_2
    df_2 %>% dplyr::mutate(`CT-20_cauchy_padj` = p.adjust(`CT-20_cauchy_pval`, method = "BH")) %>% 
      dplyr::mutate(`CT-8_cauchy_padj` = p.adjust(`CT-8_cauchy_pval`, method = "BH"))  -> df_2
    CT8_SVG = df_2 %>% dplyr::filter(`CT-8_cauchy_padj` < 0.05) %>% .$gene
    CT20_SVG = df_2 %>% dplyr::filter(`CT-20_cauchy_padj` < 0.05) %>% .$gene
    list(
      CT8_SVG = CT8_SVG,
      CT20_SVG = CT20_SVG
    ) ->> CT_SVG
    ggvenn::ggvenn(
      list(
        CT8_SVG = CT8_SVG,
        CT20_SVG = CT20_SVG
      )
    )
  } -> CT_sepcfic_SVG_venn ##2.3CT_sepcfic_SVG_venn ----

load("~/Downloads/cauchy_scn_core_vs_shell.Rdata")
cauchy_scn_core_vs_shell %>%
  rowwise() %>% 
  dplyr::mutate(`CT-8_avg_log2FC` = mean(c(`CT-8-1_avg_log2FC`, `CT-8-2_avg_log2FC`, `CT-8-3_avg_log2FC`))) %>% 
  dplyr::mutate(`CT-20_avg_log2FC` = mean(c(`CT-20-1_avg_log2FC`, `CT-20-2_avg_log2FC`, `CT-20-3_avg_log2FC`))) %>% 
  as.data.frame() -> cauchy_scn_core_vs_shell

cauchy_scn_core_vs_shell %>% 
  dplyr::mutate(SVG_CT8_group = case_when(
    (CT8_cauchy_padj < 0.05)&(`CT-8_avg_log2FC` > 0) ~ "SCN core",
    (CT8_cauchy_padj < 0.05)&(`CT-8_avg_log2FC` < 0) ~ "SCN shell",
    TRUE ~ "None"
  )) %>% 
  dplyr::mutate(SVG_CT20_group = case_when(
    (CT20_cauchy_padj < 0.05)&(`CT-20_avg_log2FC` > 0) ~ "SCN core",
    (CT20_cauchy_padj < 0.05)&(`CT-20_avg_log2FC` < 0) ~ "SCN shell",
    TRUE ~ "None"
  )) -> cauchy_scn_core_vs_shell

CT8_MOD %>% 
  purrr::map(function(df_){
#    gene_se = CT_SVG$CT8_SVG
#    df_ %>% dplyr::filter(gene %in% gene_se) -> df_1
#    sprintf("df_ = %s, df_1 = %s", nrow(df_), nrow(df_1))
    cauchy_scn_core_vs_shell = cauchy_scn_core_vs_shell %>% dplyr::select(gene, SVG_CT8_group, SVG_CT20_group)
    left_join(x = df_, y = cauchy_scn_core_vs_shell, by = "gene") -> df_1
    df_1 %>% dplyr::filter(SVG_CT8_group != "None") -> df_1
    df_1$SVG_CT8_group = factor(df_1$SVG_CT8_group, levels = c("SCN core", "SCN shell"))
    
    comparisons <- list(
      c("SCN core", "SCN shell")
    )
    
#    df_1 = df_1[1:100, ]
    
    df_1 %>% 
      ggplot(aes(x = SVG_CT8_group, y = mean_KWithin, fill = SVG_CT8_group)) + 
      geom_boxplot() + 
      stat_compare_means(
        comparisons = comparisons,
        method = "t.test",
        label = "p.format"
      ) + 
      theme_classic()
  }) %>% patchwork::wrap_plots(., ncol = 2, guides = "collect") -> p_CT8

CT20_MOD %>% 
  purrr::map(function(df_){
#    gene_se = CT_SVG$CT20_SVG
#    df_ %>% dplyr::filter(gene %in% gene_se) -> df_1
#    sprintf("df_ = %s, df_1 = %s", nrow(df_), nrow(df_1))
    cauchy_scn_core_vs_shell = cauchy_scn_core_vs_shell %>% dplyr::select(gene, SVG_CT8_group, SVG_CT20_group)
    left_join(x = df_, y = cauchy_scn_core_vs_shell, by = "gene") -> df_1
    df_1 %>% dplyr::filter(SVG_CT20_group != "None") -> df_1
    df_1$SVG_CT20_group = factor(df_1$SVG_CT20_group, levels = c("SCN core", "SCN shell"))
    
    comparisons <- list(
      c("SCN core", "SCN shell")
    )
    
#    df_1 = df_1[1:100, ]
    
    df_1 %>% 
      ggplot(aes(x = SVG_CT20_group, y = mean_KWithin, fill = SVG_CT20_group)) + 
      geom_boxplot() + 
      stat_compare_means(
        comparisons = comparisons,
        method = "t.test",
        label = "p.format"
      ) + 
      theme_classic()
  }) %>% patchwork::wrap_plots(., ncol = 2, guides = "collect") -> p_CT20

patchwork::wrap_plots(p_CT8[[1]], p_CT20[[1]], guides = "collect") -> p_core
patchwork::wrap_plots(p_CT8[[2]], p_CT20[[2]]) -> p_shell
ggsave(plot = p_core, filename = "~/Downloads/p_core.pdf", width = 7.04, height = 4.44)
ggsave(plot = p_shell, filename = "~/Downloads/p_shell.pdf", width = 7.04, height = 4.44)

p_CT8[[1]]$data$SVG_CT8_group %>% table()
p_CT8[[2]]$data$SVG_CT8_group %>% table()

p_CT20[[1]]$data$SVG_CT20_group %>% table()
p_CT20[[2]]$data$SVG_CT20_group %>% table()

p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN core") %>% head(10)
p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN shell") %>% head(10)

p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN core") %>% head(10)
p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN shell") %>% head(10)

p_CT8[[2]]$data %>% dplyr::filter(SVG_CT8_group == "SCN core") %>% head(10) %>% .$gene
p_CT8[[2]]$data %>% dplyr::filter(SVG_CT8_group == "SCN shell") %>% head(10) %>% .$gene

p_CT20[[2]]$data %>% dplyr::filter(SVG_CT20_group == "SCN core") %>% head(10) %>% .$gene
p_CT20[[2]]$data %>% dplyr::filter(SVG_CT20_group == "SCN shell") %>% head(10)

p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN core")
p_CT20[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN core")

list(
  CT8 = p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN core") %>% .$gene,
  CT20 = p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN core") %>% .$gene
) %>% ggvenn::ggvenn() -> p_ggvenn_core

list(
  CT8 = p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN shell") %>% .$gene,
  CT20 = p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN shell") %>% .$gene
) %>% ggvenn::ggvenn() -> p_ggvenn_shell

list(
  CT8 = p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN core") %>% .$gene,
  CT20 = p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN core") %>% .$gene
) %>% ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    CT8_only = df_ %>% dplyr::filter(CT8 == T, CT20 == F) %>% .$key
    CT20_only = df_ %>% dplyr::filter(CT8 == F, CT20 == T) %>% .$key
    Both = df_ %>% dplyr::filter(CT8 == T, CT20 == T) %>% .$key

    p_CT8[[1]]$data %>% dplyr::mutate(scaled_score = scales::rescale(x = mean_KWithin, to = c(0,1))) %>% dplyr::filter(gene %in% CT8_only) %>% .$scaled_score -> CT8_only_score
    p_CT20[[1]]$data %>% dplyr::mutate(scaled_score = scales::rescale(x = mean_KWithin, to = c(0,1))) %>% dplyr::filter(gene %in% CT20_only) %>% .$scaled_score -> CT20_only_score
    
    ((p_CT8[[1]]$data %>% dplyr::mutate(scaled_score = scales::rescale(x = mean_KWithin, to = c(0,1))) %>% dplyr::filter(gene %in% Both) %>% .$scaled_score) + 
        (p_CT20[[1]]$data %>% dplyr::mutate(scaled_score = scales::rescale(x = mean_KWithin, to = c(0,1))) %>% dplyr::filter(gene %in% Both) %>% .$scaled_score))/2 -> Both_score
    
    t.test(x = CT8_only_score, y = CT20_only_score)
    t.test(x = CT8_only_score, y = Both_score)
    t.test(x = CT20_only_score, y = Both_score)
    
    rbind(data.frame(score = CT8_only_score, group = "CT8_only"),
          data.frame(score = CT20_only_score, group = "CT20_only"),
          data.frame(score = Both_score, group = "Both")) -> df_plot
    
    comparisons <- list(
      c("Both", "CT8_only"),
      c("Both", "CT20_only"),
      c("CT8_only", "CT20_only")
    )
    
    df_plot %>% 
      ggplot(aes(x = group, y = score, fill = group)) + 
      geom_boxplot() + 
      stat_compare_means(
        comparisons = comparisons,
        method = "t.test",
        label = "p.format"
      ) + 
      theme_classic() + 
      ylab("scaled hub score") + 
      ggtitle("SCN core")
  }

list(
  CT8 = p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN shell") %>% .$gene,
  CT20 = p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN shell") %>% .$gene
) %>% ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    CT8_only = df_ %>% dplyr::filter(CT8 == T, CT20 == F) %>% .$key
    CT20_only = df_ %>% dplyr::filter(CT8 == F, CT20 == T) %>% .$key
    Both = df_ %>% dplyr::filter(CT8 == T, CT20 == T) %>% .$key
    
    p_CT8[[1]]$data %>% dplyr::mutate(scaled_score = scales::rescale(x = mean_KWithin, to = c(0,1))) %>% dplyr::filter(gene %in% CT8_only) %>% .$scaled_score -> CT8_only_score
    p_CT20[[1]]$data %>% dplyr::mutate(scaled_score = scales::rescale(x = mean_KWithin, to = c(0,1))) %>% dplyr::filter(gene %in% CT20_only) %>% .$scaled_score -> CT20_only_score
    
    ((p_CT8[[1]]$data %>% dplyr::mutate(scaled_score = scales::rescale(x = mean_KWithin, to = c(0,1))) %>% dplyr::filter(gene %in% Both) %>% .$scaled_score) + 
        (p_CT20[[1]]$data %>% dplyr::mutate(scaled_score = scales::rescale(x = mean_KWithin, to = c(0,1))) %>% dplyr::filter(gene %in% Both) %>% .$scaled_score))/2 -> Both_score
    
    t.test(x = CT8_only_score, y = CT20_only_score)
    t.test(x = CT8_only_score, y = Both_score)
    t.test(x = CT20_only_score, y = Both_score)
    
    rbind(data.frame(score = CT8_only_score, group = "CT8_only"),
          data.frame(score = CT20_only_score, group = "CT20_only"),
          data.frame(score = Both_score, group = "Both")) -> df_plot
    
    comparisons <- list(
      c("Both", "CT8_only"),
      c("Both", "CT20_only"),
      c("CT8_only", "CT20_only")
    )
    
    df_plot %>% 
      ggplot(aes(x = group, y = score, fill = group)) + 
      geom_boxplot() + 
      stat_compare_means(
        comparisons = comparisons,
        method = "t.test",
        label = "p.format"
      ) + 
      theme_classic() + 
      ylab("scaled hub score") + 
      ggtitle("SCN shell")
  }

list(
  CT8 = p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN core") %>% .$gene,
  CT20 = p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN core") %>% .$gene
) %>% ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    list(
      CT8_only = df_ %>% dplyr::filter(CT8 == T, CT20 == F) %>% .$key,
      CT20_only = df_ %>% dplyr::filter(CT8 == F, CT20 == T) %>% .$key,
      Both = df_ %>% dplyr::filter(CT8 == T, CT20 == T) %>% .$key
    ) -> list_
    
    full_join(
    x = p_CT8[[1]]$data %>% dplyr::filter(gene %in% list_$Both) %>% dplyr::mutate(rank = 1:nrow(.)) %>% 
      dplyr::select(gene, rank),
    y = p_CT20[[1]]$data %>% dplyr::filter(gene %in% list_$Both) %>% dplyr::mutate(rank = 1:nrow(.)) %>% 
      dplyr::select(gene, rank),
    by = "gene"
    ) %>% 
      rowwise() %>% dplyr::mutate(rank = (rank.x+rank.y)/2) %>% 
      dplyr::arrange(rank) %>% head(10) %>% .$gene -> Both_top_10
    CT8_top_10 = p_CT8[[1]]$data %>% dplyr::filter(gene %in% list_$CT8_only) %>% head(10) %>% .$gene
    CT20_top_10 = p_CT20[[1]]$data %>% dplyr::filter(gene %in% list_$CT20_only) %>% head(10) %>% .$gene
    
    rbind(
      data.frame(gene = Both_top_10, group = "Both") %>% dplyr::mutate(rank = nrow(.):1),
      data.frame(gene = CT8_top_10, group = "CT8_only") %>% dplyr::mutate(rank = nrow(.):1),
      data.frame(gene = CT20_top_10, group = "CT20_only") %>% dplyr::mutate(rank = nrow(.):1)
    )
  } -> df_plot_SCN_core

df_plot_SCN_core %>% 
  {
    df_plot = .
    df_plot$gene = factor(df_plot$gene, levels = df_plot$gene)
    
    df_plot %>% 
      dplyr::mutate(
        gene_id = row_number(),
        angle = 90 - 360 * (gene_id - 0.5) / n(),
        hjust = ifelse(angle < -90, 1, 0),
        angle_label = ifelse(angle < -90, angle + 180, angle)
      ) -> df_plot
    
    df_plot %>% 
      ggplot(aes(x = gene, y = 1, color = group, size = rank)) + 
      geom_point() + 
      geom_text(
        aes(
          y = 1.15,
          label = gene,
          angle = angle_label,
          hjust = hjust
        ),
        size = 3,
        show.legend = FALSE
      ) +
      coord_polar(start = 0) +
      ylim(0, 1.4) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)
      ) + 
      ggtitle("SCN core top ranked co-expressed genes")
  } -> p1

list(
  CT8 = p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN shell") %>% .$gene,
  CT20 = p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN shell") %>% .$gene
) %>% ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    list(
      CT8_only = df_ %>% dplyr::filter(CT8 == T, CT20 == F) %>% .$key,
      CT20_only = df_ %>% dplyr::filter(CT8 == F, CT20 == T) %>% .$key,
      Both = df_ %>% dplyr::filter(CT8 == T, CT20 == T) %>% .$key
    ) -> list_
    
    full_join(
      x = p_CT8[[1]]$data %>% dplyr::filter(gene %in% list_$Both) %>% dplyr::mutate(rank = 1:nrow(.)) %>% 
        dplyr::select(gene, rank),
      y = p_CT20[[1]]$data %>% dplyr::filter(gene %in% list_$Both) %>% dplyr::mutate(rank = 1:nrow(.)) %>% 
        dplyr::select(gene, rank),
      by = "gene"
    ) %>% 
      rowwise() %>% dplyr::mutate(rank = (rank.x+rank.y)/2) %>% 
      dplyr::arrange(rank) %>% head(10) %>% .$gene -> Both_top_10
    CT8_top_10 = p_CT8[[1]]$data %>% dplyr::filter(gene %in% list_$CT8_only) %>% head(10) %>% .$gene
    CT20_top_10 = p_CT20[[1]]$data %>% dplyr::filter(gene %in% list_$CT20_only) %>% head(10) %>% .$gene
    
    rbind(
      data.frame(gene = Both_top_10, group = "Both") %>% dplyr::mutate(rank = nrow(.):1),
      data.frame(gene = CT8_top_10, group = "CT8_only") %>% dplyr::mutate(rank = nrow(.):1),
      data.frame(gene = CT20_top_10, group = "CT20_only") %>% dplyr::mutate(rank = nrow(.):1)
    )
  } -> df_plot_SCN_shell

df_plot_SCN_shell %>% 
  {
    df_plot = .
    df_plot$gene = factor(df_plot$gene, levels = df_plot$gene)
    
    df_plot %>% 
      dplyr::mutate(
        gene_id = row_number(),
        angle = 90 - 360 * (gene_id - 0.5) / n(),
        hjust = ifelse(angle < -90, 1, 0),
        angle_label = ifelse(angle < -90, angle + 180, angle)
      ) -> df_plot
    
    df_plot %>% 
      ggplot(aes(x = gene, y = 1, color = group, size = rank)) + 
      geom_point() + 
      geom_text(
        aes(
          y = 1.15,
          label = gene,
          angle = angle_label,
          hjust = hjust
        ),
        size = 3,
        show.legend = FALSE
      ) +
      coord_polar(start = 0) +
      ylim(0, 1.4) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)
      ) + 
      ggtitle("SCN shell top ranked co-expressed genes")
  } -> p2

patchwork::wrap_plots(p1,p2)


list(
CT8 = p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN core") %>% .$gene,
CT20 = p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN core") %>% .$gene
) %>% ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    list(
      CT8_only = df_ %>% dplyr::filter(CT8 == T, CT20 == F) %>% .$key,
      CT20_only = df_ %>% dplyr::filter(CT8 == F, CT20 == T) %>% .$key,
      Both = df_ %>% dplyr::filter(CT8 == T, CT20 == T) %>% .$key
    ) %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(genes, group_){
        data.frame(gene = genes, group = group_)
      }) %>% do.call(rbind, .) %>% "rownames<-"(.,NULL) -> df_plot
    
    df_plot$gene = factor(df_plot$gene, levels = df_plot$gene)
    
    df_plot %>% 
      dplyr::mutate(
        gene_id = row_number(),
        angle = 90 - 360 * (gene_id - 0.5) / n(),
        hjust = ifelse(angle < -90, 1, 0),
        angle_label = ifelse(angle < -90, angle + 180, angle)
      ) -> df_plot
    
    df_plot %>% 
    ggplot(aes(x = gene, y = 1, color = group)) + 
      geom_point(size = 2) + 
      geom_text(
        aes(
          y = 1.15,
          label = gene,
          angle = angle_label,
          hjust = hjust
        ),
        size = 3,
        show.legend = FALSE
      ) +
      coord_polar(start = 0) +
      ylim(0, 1.4) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)
      ) + 
      ggtitle("SCN core top 10 co-expressed genes")
  } -> p_SCN_core

list(
  CT8 = p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN shell") %>% head(10) %>% .$gene,
  CT20 = p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN shell") %>% head(10) %>% .$gene
) %>% ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    list(
      CT8_only = df_ %>% dplyr::filter(CT8 == T, CT20 == F) %>% .$key,
      CT20_only = df_ %>% dplyr::filter(CT8 == F, CT20 == T) %>% .$key,
      Both = df_ %>% dplyr::filter(CT8 == T, CT20 == T) %>% .$key
    ) %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(genes, group_){
        data.frame(gene = genes, group = group_)
      }) %>% do.call(rbind, .) %>% "rownames<-"(.,NULL) -> df_plot
    
    df_plot$gene = factor(df_plot$gene, levels = df_plot$gene)
    
    df_plot %>% 
      dplyr::mutate(
        gene_id = row_number(),
        angle = 90 - 360 * (gene_id - 0.5) / n(),
        hjust = ifelse(angle < -90, 1, 0),
        angle_label = ifelse(angle < -90, angle + 180, angle)
      ) -> df_plot
    
    df_plot %>% 
      ggplot(aes(x = gene, y = 1, color = group)) + 
      geom_point(size = 2) + 
      geom_text(
        aes(
          y = 1.15,
          label = gene,
          angle = angle_label,
          hjust = hjust
        ),
        size = 3,
        show.legend = FALSE
      ) +
      coord_polar(start = 0) +
      ylim(0, 1.4) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)
      ) + 
      ggtitle("SCN shell top 10 co-expressed genes")
  } -> p_SCN_shell

patchwork::wrap_plots(p_SCN_core, p_SCN_shell, guides = "collect") -> p_SCN_coexpression
ggsave(plot = p_SCN_coexpression, filename = "~/Downloads/p_SCN_coexpression.pdf", width = 12.42, height = 5.10)

p_CT8[[1]]$data %>% dplyr::filter(SVG_CT8_group == "SCN shell") %>% head(10) %>% .$gene
p_CT20[[1]]$data %>% dplyr::filter(SVG_CT20_group == "SCN shell") %>% head(10) %>% .$gene




list(
  CT8 = p_CT8[[1]]$data %>% head(10),
  CT20 = p_CT20[[1]]$data %>% head(10)
) %>% 
  purrr::map(function(df_){
    df_$gene
  })

p_CT8[[2]]$data %>% head(10)
p_CT20[[2]]$data %>% head(10)

p_CT8
p_CT8[[1]]$data

p_CT20[[2]]$data
