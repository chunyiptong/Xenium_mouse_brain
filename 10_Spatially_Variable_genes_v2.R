library(tidyverse)
library(Seurat)
library(SPARK)

# 1. load data ----
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-annotated-mini.rda")
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")

meta_data = xenium_obj_filtered@meta.data

list.files(path = "~/Dropbox/Spatial_transcriptomics_new/rda/Brain_region_annotation/SCN_SON_PVN", pattern = "v5\\.csv", full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)_regions.+", "\\1", .) %>% gsub("_", "-", .)) %>% 
  #    .[1] %>% 
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

GetAssayData(xenium_obj_filtered, assay = "Xenium", layer = "counts") -> counts_matrix
rm(xenium_obj_filtered);gc()

names(xy_coordinates) %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    SCN_cells$cell %>% {.[grepl(sample_, .)]} -> cells_
    xy_coordinates[[sample_]] %>% dplyr::filter(cell %in% cells_)
  }) -> xy_coordinates

xy_coordinates %>% purrr::map(function(df_){nrow(df_)})
SCN_cells$cell %>% {gsub("(CT.+?)_.+", "\\1", .)} %>% table()

#2. Run SPARK-X ----
run_sparkx = function(counts_matrix, xy_coordintates, cells_list, Sample){
  xy_coordinates[[Sample]] %>% column_to_rownames("cell") %>% .[cells_list[[Sample]], ] %>% as.matrix() -> locus_in_
  counts_matrix[,rownames(locus_in_)] -> count_in_
  sx = SPARK::sparkx(count_in = count_in_, locus_in = locus_in_, numCores = 4)
  return(sx)
}

c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(sample_){
    SCN_cells$cell %>% {.[grepl(sample_, .)]}
  }) -> cells_list


c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(Sample_){
    sx = run_sparkx(counts_matrix = counts_matrix, xy_coordintates = xy_coordinates, cells_list = cells_list, Sample = Sample_)
  }) -> sparkx_res_list
save(sparkx_res_list, file = "~/Downloads/sparkx_res_list.rda")

sparkx_res_list %>% 
  purrr::map(function(df_){
    df_$res_mtest[c("Vip", "Avp"), ]
  })

counts_matrix %>% colSums() %>% median() -> median_
counts_matrix %>% as.data.frame() %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*median_})) -> norm_counts_matrix
save(counts_matrix, norm_counts_matrix, xy_coordinates, file = "~/Downloads/SCN_mats.rda")
sparkx_res_list %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, sample_){
    df_$res_mtest %>% dplyr::arrange(combinedPval) %>% .[1:10,] -> df_
    df_res = df_
    rownames(df_) %>% 
      "names<-"(.,.) %>% 
      #      .[1] %>% 
      purrr::map(function(gene_){
        norm_counts_matrix[gene_, ] %>% 
          t() %>% 
          as.data.frame() %>% "colnames<-"(., "normalized_expr") %>% 
          rownames_to_column("cell") -> df_
        xy_coordinates[[sample_]] -> df_xy
        inner_join(x = df_, y = df_xy, by = "cell") -> df_plot
        
        pval = df_res[gene_, ]$adjustedPval %>% sprintf("%.2e", .)
        
        df_plot %>% ggplot(aes(x = x, y = y, color = normalized_expr)) + 
          geom_point(size = 1) + 
          scale_color_gradient(low = "lightgrey", high = "red") + 
          theme_void() + 
          labs(color = "Normalized\nExpression") +
          ggtitle(sprintf("%s: p.adj=%s", gene_, pval))
      }) %>% patchwork::wrap_plots(., ncol = 4) + 
      patchwork::plot_annotation(title = sample_, theme = theme(plot.title = element_text(hjust = 0.5)))
  }) -> p_list ##2.1 per sample top 10 SVGs -----

p_list$`CT-8-1`
p_list$`CT-8-2`
p_list$`CT-8-3`
p_list$`CT-20-1`
p_list$`CT-20-2`
p_list$`CT-20-3`

sparkx_res_list %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, sample_){
    df_ = df_$res_mtest
    df_ = df_ %>% dplyr::filter(adjustedPval < 0.05)
    gene_ = rownames(df_)
  }) %>% ggvenn::list_to_data_frame() %>% 
  dplyr::mutate(across(!matches("key"), function(x){x %>% as.integer()})) %>%
  as.data.frame() %>%
  UpSetR::upset(sets = rev(c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3")), keep.order = T, order.by = "freq",
                #                main.bar.color = color_list[[direction_]],
                #                sets.bar.color = color_list[[direction_]], 
                #                matrix.color = color_list[[direction_]] -> p_upset
  ) -> p_upset
p_upset

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

sparkx_res_list_per_CT = list()
c("CT-8-1", "CT-8-2", "CT-8-3") %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(sample_){
    sparkx_res_list[[sample_]]$res_mtest %>% 
      dplyr::select(combinedPval) %>% 
      rownames_to_column("gene")
  }) %>% purrr::reduce(., inner_join, by = "gene") %>% 
  column_to_rownames("gene") %>% 
  apply(., 1, function(x){cauchy.test(x)}) %>% 
  as.data.frame() %>% 
  "colnames<-"(., "cauchy_pval") %>% 
  dplyr::mutate(cauchy_padj = p.adjust(cauchy_pval, method = "BH")) %>% 
  dplyr::arrange(cauchy_padj) -> sparkx_res_list_per_CT[["CT-8"]]

c("CT-20-1", "CT-20-2", "CT-20-3") %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(sample_){
    sparkx_res_list[[sample_]]$res_mtest %>% 
      dplyr::select(combinedPval) %>% 
      rownames_to_column("gene")
  }) %>% purrr::reduce(., inner_join, by = "gene") %>% 
  column_to_rownames("gene") %>% 
  apply(., 1, function(x){cauchy.test(x)}) %>% 
  as.data.frame() %>% 
  "colnames<-"(., "cauchy_pval") %>% 
  dplyr::mutate(cauchy_padj = p.adjust(cauchy_pval, method = "BH")) %>% 
  dplyr::arrange(cauchy_padj) -> sparkx_res_list_per_CT[["CT-20"]]

sparkx_res_list %>% 
  purrr::map2(.x=., .y=names(.), .f=function(res_, sample_){
    res_$res_mtest %>% 
      #      dplyr::select(combinedPval) %>% 
      "colnames<-"(., sprintf("%s_%s", sample_, colnames(.))) %>% 
      rownames_to_column("gene")
  }) %>% purrr::reduce(., inner_join, by = "gene") -> df_combinedPval_all_samples

sparkx_res_list_per_CT %>% 
  purrr::map2(.x=.,.y=names(.), .f=function(df_, CT_){
    df_ %>% 
      "colnames<-"(., sprintf("%s_%s", CT_, colnames(df_))) %>% 
      rownames_to_column("gene")
  }) %>% purrr::reduce(., inner_join, by = "gene") %>% 
  inner_join(x = df_combinedPval_all_samples, y = ., by ="gene") -> df_combinedPval_all_samples_sparkx

save(df_combinedPval_all_samples_sparkx, file = "~/Downloads/df_combinedPval_all_samples_sparkx.Rdata") 
##2.2 SVG results ----

counts_matrix %>% rowSums() %>% {.[.!=0]} %>% names() -> gene_selected
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
  } -> CT_sepcfic_SVG_venn #Figure 3A before filtering
##2.3CT_sepcfic_SVG_venn ----

CT_SVG %>% 
  ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    list(
    CT20_only = df_ %>% dplyr::filter(CT20_SVG == TRUE, CT8_SVG == FALSE) %>% .$key,
    CT8_only = df_ %>% dplyr::filter(CT8_SVG == TRUE, CT20_SVG == FALSE) %>% .$key,
    Both = df_ %>% dplyr::filter(CT8_SVG == TRUE, CT20_SVG == TRUE) %>% .$key
    )
  } -> gene_list

load("~/Downloads/DEG_dfs.rda")
DEG_df_SCN %>% 
  dplyr::filter(p_adj.loc < 0.05, abs(logFC) > log2(1.1)) %>% 
  {
    df_ = .
    df_$Sig %>% table() %>% 
      names() %>% 
      "names<-"(., c("CT20_TVG", "CT8_TVG")) %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(Sig_, group_){
        df_ %>% dplyr::filter(Sig == Sig_) %>% .$gene -> genes_
        data.frame(gene = genes_, group_TVG = group_)
      }) %>% do.call(rbind, .) %>% 
      "rownames<-"(., NULL) -> df_TVG
    gene_list %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(genes_, group_){
        data.frame(gene = genes_, group_SVG = group_)
      }) %>% do.call(rbind, .) %>% 
      "rownames<-"(., NULL) %>% 
      dplyr::mutate(group_SVG = case_when(
        group_SVG == "CT20_only" ~ "CT20_SVG_only",
        group_SVG == "CT8_only" ~ "CT8_SVG_only",
        group_SVG == "Both" ~ "Both_SVG"
      )) -> df_SVG
    left_join(df_SVG, df_TVG, by = "gene") %>% 
      dplyr::mutate(group_TVG = case_when(
        is.na(group_TVG) ~ "Not_TVG",
        TRUE ~ group_TVG
      )) -> df_plot
    table(df_plot$group_SVG, df_plot$group_TVG) %>% 
      as.data.frame() %>% 
      "colnames<-"(., c("SVG", "TVG", "Freq")) -> df_plot
    levels(df_plot$SVG) %>% 
      purrr::map(function(x_){
        df_plot %>% dplyr::filter(SVG == x_) -> df_
        round((df_$Freq/sum(df_$Freq))*100, 2) -> percent_
        df_ %>% dplyr::mutate(label = sprintf("%s\n(%s%%)", Freq, percent_))
      }) %>% do.call(rbind, .) -> df_plot
    
    df_plot$SVG = factor(df_plot$SVG, levels = c("CT8_SVG_only", "Both_SVG", "CT20_SVG_only"))
    df_plot$TVG = factor(df_plot$TVG, levels = c("CT8_TVG", "CT20_TVG", "Not_TVG"))
    df_plot %>% 
      ggplot(aes(x = SVG, y = Freq, fill = TVG)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_classic() + 
      geom_text(aes(label = label), position = position_dodge(width = 0.9), vjust = 0.5) + 
      scale_fill_manual(values = c("CT8_TVG" = "#f86a62", "CT20_TVG" = "#06b7bc", "Not_TVG" = "grey")) + 
      ylab("Number of genes")
  } -> p
ggsave(plot = p, filename = "~/Downloads/SVG_TVG_barplot.pdf", width = 8.24, height = 4.10)

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

cauchy_scn_core_vs_shell
CT_SVG %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(genes_, group_){
    CT_ = gsub("(CT\\d+)_.+", "\\1", group_)
    cauchy_scn_core_vs_shell %>% dplyr::filter(gene %in% genes_) %>% 
    dplyr::select(gene, SVG_CT8_group, SVG_CT20_group) %>% 
      pivot_longer(cols = -gene, names_to = "group", values_to = "region") -> df_
    df_ %>% dplyr::mutate(CT = gsub(".+_(CT\\d+)_.+", "\\1", group)) -> df_
    df_ %>% dplyr::filter(CT == CT_) -> df_
    df_$region %>% table() %>% 
      as.data.frame() %>% 
      "colnames<-"(., c("region", "freq")) %>% 
      dplyr::mutate(group = group_) %>% 
      dplyr::mutate(label = sprintf("%s\n(%s%%)", freq, round((freq/sum(freq))*100, 2))) -> df_
  }) %>% do.call(rbind, .) %>% {
    df_ = .
    rownames(df_) = NULL
    df_$region = factor(df_$region, levels = c("SCN core", "SCN shell", "None"))
    df_$group = factor(df_$group, levels = c("CT8_SVG", "CT20_SVG"))

    df_ %>% 
      ggplot(aes(x = group, y = freq, fill = region)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_classic() + 
      geom_text(aes(label = label), position = position_dodge(width = 0.9), vjust = 0.5) + 
      scale_fill_manual(values = c("SCN shell" = "#f86a62", "SCN core" = "#06b7bc", "None" = "grey")) + 
      ylab("Number of genes")
  } -> p

ggsave(plot = p, filename = "~/Downloads/SVG_TVG_barplot_1.pdf", width = 8.24, height = 4.10)


cauchy_scn_core_vs_shell %>% dplyr::filter(gene %in% gene_list$Both) %>% 
  dplyr::filter(gene %in% (DEG_df_SCN %>% dplyr::filter(Sig == "Not_significant") %>% .$gene)) %>% 
  dplyr::filter(SVG_CT8_group != SVG_CT20_group) %>% 
  dplyr::filter((SVG_CT8_group != "None")&(SVG_CT20_group != "None")) %>% {
    df_ = .
    DEG_df_SCN %>% dplyr::filter(gene %in% df_$gene) %>% dplyr::filter(Sig == "Not_significant") %>% .$gene -> genes_
    df_ %>% dplyr::filter(gene %in% genes_) #%>% .$gene
  }
#"Slc39a6" "Gad2"    "Fgf12"  

cauchy_scn_core_vs_shell %>% dplyr::filter(gene %in% gene_list$CT20_only) %>% 
  dplyr::filter(gene %in% (DEG_df_SCN %>% dplyr::filter(Sig == "Not_significant") %>% .$gene)) %>% 
  dplyr::filter(SVG_CT20_group != "None", SVG_CT8_group == "None") %>% 
  dplyr::select(gene, CT8_cauchy_padj, CT20_cauchy_padj, `CT-8_avg_log2FC`, `CT-20_avg_log2FC`, SVG_CT8_group, SVG_CT20_group) %>% .$gene -> SVG_CT20

DEG_df_SCN$gene -> se_genes
cauchy_scn_core_vs_shell %>% 
  dplyr::filter(gene %in% se_genes) %>% 
  dplyr::select(gene, CT8_cauchy_pval, CT8_cauchy_padj, CT20_cauchy_pval, `CT20_cauchy_padj`, `CT-8_avg_log2FC`, `CT-20_avg_log2FC`, SVG_CT8_group, SVG_CT20_group) %>% 
  dplyr::left_join(x = ., y = DEG_df_SCN[,c("gene", "Sig")], by = "gene") %>% 
  dplyr::mutate(SVG_Sig = case_when(
    gene %in% gene_list$CT8_only ~ "CT8_only",
    gene %in% gene_list$CT20_only ~ "CT20_only",
    gene %in% gene_list$Both ~ "Both"
  )) -> df_plot

df_plot %>% 
  dplyr::filter(SVG_Sig %in% c("CT8_only", "Both")) %>% 
  dplyr::mutate(TVG_Sig = case_when(
    grepl("enriched", Sig) ~ "TVG", 
    TRUE ~ "Not_TVG"
  )) %>% 
  {
    df_plot = .
    
    rbind(
    df_plot %>% dplyr::filter(SVG_CT8_group == "SCN core") %>% dplyr::arrange(CT8_cauchy_pval, desc(`CT-8_avg_log2FC`)) %>% head(10),
    df_plot %>% dplyr::filter(SVG_CT8_group == "SCN shell") %>% dplyr::arrange(CT8_cauchy_pval, desc(abs(`CT-8_avg_log2FC`))) %>% head(10)
    ) -> top10
    
    df_plot %>% 
      ggplot(aes(x = `CT-8_avg_log2FC`, y = -log10(CT8_cauchy_padj), color = TVG_Sig)) + 
      geom_point(alpha = 0.7) + 
      geom_text_repel(
        data = top10,
        aes(
          x = `CT-8_avg_log2FC`,
          y = -log10(CT8_cauchy_padj),
          label = gene
        ),
        size = 3,
        box.padding = 0.4,
        point.padding = 0.3,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 123,
        show.legend = FALSE
      ) + 
      theme_classic() + 
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + 
      ggtitle("CT8") + 
      ylab("-log10(Core vs Shell adj. P)") + 
      xlab("log2FC")
  } -> p1


df_plot %>% 
  dplyr::filter(SVG_Sig %in% c("CT20_only", "Both")) %>% 
  dplyr::mutate(TVG_Sig = case_when(
    grepl("enriched", Sig) ~ "TVG", 
    TRUE ~ "Not_TVG"
  )) %>% 
  {
    df_plot = .
    
    rbind(
      df_plot %>% dplyr::filter(SVG_CT20_group == "SCN core") %>% dplyr::arrange(CT20_cauchy_pval, desc(`CT-20_avg_log2FC`)) %>% head(10),
      df_plot %>% dplyr::filter(SVG_CT20_group == "SCN shell") %>% dplyr::arrange(CT20_cauchy_pval, desc(abs(`CT-20_avg_log2FC`))) %>% head(10)
    ) -> top10
    
    df_plot %>% 
      ggplot(aes(x = `CT-20_avg_log2FC`, y = -log10(CT20_cauchy_padj), color = TVG_Sig)) + 
      geom_point(alpha = 0.7) + 
      geom_text_repel(
        data = top10,
        aes(
          x = `CT-20_avg_log2FC`,
          y = -log10(CT20_cauchy_padj),
          label = gene
        ),
        size = 3,
        box.padding = 0.4,
        point.padding = 0.3,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 123,
        show.legend = FALSE
      ) + 
      theme_classic() + 
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + 
      ggtitle("CT20") + 
      ylab("-log10(Core vs Shell adj. P)") + 
      xlab("log2FC")
  } -> p2

patchwork::wrap_plots(p1, p2, ncol = 2, guides = "collect") -> p3

DEG_df_SCN %>% dplyr::mutate(TVG = case_when(
  grepl("enriched", Sig) ~ "Yes",
  TRUE ~ "No"
)) %>% dplyr::select(gene,  TVG) -> df_1


cauchy_scn_core_vs_shell %>% 
  dplyr::filter(gene %in% se_genes) %>% 
  dplyr::filter(gene %in% (gene_list %>% unlist() %>% unname() %>% unique())) %>% 
  left_join(x = ., y = df_1, by = "gene") %>% 
  dplyr::filter(TVG == "Yes") %>% 
  {
    df_plot = .
    df_plot = df_plot %>% dplyr::mutate(Sig = case_when(
      (CT8_cauchy_padj < 0.05 )&(CT20_cauchy_padj < 0.05)&(abs(`CT-8_avg_log2FC`) > 0.2)&(abs(`CT-20_avg_log2FC`) > 0.2) ~ "Sig",
      TRUE ~ "Not_Sig"
    ))
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC` < 0, `CT-20_avg_log2FC` > 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC`)), desc(abs(`CT-20_avg_log2FC`))) %>% 
      head(10) -> df_1
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC` < 0, `CT-20_avg_log2FC` < 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC`)), desc(abs(`CT-20_avg_log2FC`))) %>% 
      head(10) -> df_2
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC` > 0, `CT-20_avg_log2FC` > 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC`)), desc(abs(`CT-20_avg_log2FC`))) %>% 
      head(10) -> df_3
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC` > 0, `CT-20_avg_log2FC` < 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC`)), desc(abs(`CT-20_avg_log2FC`))) %>% 
      head(10) -> df_4
    rbind(df_1, df_2, df_3, df_4) -> top10
    
    top10
    
    df_plot %>% 
      ggplot(aes(x = `CT-8_avg_log2FC`, y = `CT-20_avg_log2FC`, color = Sig)) + 
      geom_point(alpha = 0.5) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + 
      theme_classic() + 
      geom_text_repel(
        data = top10,
        aes(
          x = `CT-8_avg_log2FC`,
          y = `CT-20_avg_log2FC`,
          label = gene, 
          color = Sig
        ),
        size = 3,
        box.padding = 0.4,
        point.padding = 0.3,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 123,
        show.legend = F,
        inherit.aes =T 
      ) + 
      xlab("Core vs shell log2FC (CT8)") + 
      ylab("log2FC (CT20)") + 
      ggtitle("Spatiotemporal genes (SVG + TVG) enriched in core/shell")
    
  } -> p4

ggsave(plot = p4, filename = "~/Downloads/Spatiotemporal_TVG_and_SVG_geom_point.pdf", width = 9.35, height = 5.35)
save(p4, file = "~/Downloads/Spatiotemporal_TVG_and_SVG_geom_point.rdata")
p4$data %>% dplyr::filter(Sig == "Sig") %>% dim()

list(
SVG = gene_list %>% unlist() %>% unname() %>% unique(),
TVG = DEG_df_SCN %>% dplyr::filter(Sig != "Not_significant") %>% .$gene
) %>% ggvenn::ggvenn() -> p5
ggsave(plot = p5, filename = "~/Downloads/TVG_vs_SVG.pdf", width = 5.79, height = 3.93)



p4$data %>% dplyr::filter(SVG_CT8_group != "None", SVG_CT20_group != "None", TVG == "Yes", Sig == "Sig") %>% 
  {
    df_ = .
    DEG_df_SCN %>% dplyr::select(gene, Sig) -> df_1
    left_join(x = df_, y = df_1, by = "gene") %>% 
      dplyr::select(gene, SVG_CT8_group, SVG_CT20_group, Sig.y) %>% 
      "colnames<-"(., c("gene", "SVG_CT8_group", "SVG_CT20_group", "TVG_Sig")) %>% 
      dplyr::arrange(desc(TVG_Sig)) -> df_2
    
    df_2 %>% dplyr::filter(TVG_Sig == "CT8_enriched") %>% 
      dplyr::arrange(SVG_CT8_group, SVG_CT20_group) -> df_CT8
    
    df_2 %>% dplyr::filter(TVG_Sig == "CT20_enriched") %>% 
      dplyr::arrange(SVG_CT8_group, SVG_CT20_group) -> df_CT20
    
    rbind(df_CT8, df_CT20)
  }


c("Slc39a6", "Gad2", "Fgf12", "Vip", "Rgs16", "Sema5a", "Rpa1", "Grin2c") %>% 
  "names<-"(.,.) %>% 
  #      .[1] %>% 
  purrr::map(function(gene_){
    c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(sample_){
        norm_counts_matrix[gene_, ] %>% 
          t() %>% 
          as.data.frame() %>% "colnames<-"(., "normalized_expr") %>% 
          rownames_to_column("cell") -> df_
        xy_coordinates[[sample_]] -> df_xy
        inner_join(x = df_, y = df_xy, by = "cell") -> df_plot
        
        df_combinedPval_all_samples_sparkx %>% dplyr::filter(gene == gene_) %>% 
          dplyr::select(matches(sample_)) %>% .[1,2] %>% sprintf("%.2e", .) -> pval
        
        #            pval = sparkx_res_list[[sample_]]$res_mtest[gene_, ]$adjustedPval %>% sprintf("%.2e", .)
        
        df_plot %>% 
          dplyr::arrange(normalized_expr) -> df_plot
        upper = quantile(df_plot$normalized_expr, 0.99, na.rm = TRUE)
        df_plot %>% 
          ggplot(aes(x = x, y = y, color = normalized_expr, size = normalized_expr)) + 
          #              geom_point(size = 1) + 
          geom_point(alpha = 0.8) + 
          scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, upper), oob = scales::squish) + 
          scale_size_continuous(range = c(2, 5)) +
          theme_void() + 
          labs(color = "Normalized\nExpression") +
          ggtitle(sprintf("%s, %s: p.adj=%s", sample_, gene_, pval)) + 
          theme(legend.position = "none")
      }) %>% patchwork::wrap_plots(., ncol = 3) -> p
  }) -> p_list

SVG_CT20 %>% 
  "names<-"(.,.) %>% 
  #      .[1] %>% 
  purrr::map(function(gene_){
    c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(sample_){
        norm_counts_matrix[gene_, ] %>% 
          t() %>% 
          as.data.frame() %>% "colnames<-"(., "normalized_expr") %>% 
          rownames_to_column("cell") -> df_
        xy_coordinates[[sample_]] -> df_xy
        inner_join(x = df_, y = df_xy, by = "cell") -> df_plot
        
        df_combinedPval_all_samples_sparkx %>% dplyr::filter(gene == gene_) %>% 
          dplyr::select(matches(sample_)) %>% .[1,2] %>% sprintf("%.2e", .) -> pval
        
        #            pval = sparkx_res_list[[sample_]]$res_mtest[gene_, ]$adjustedPval %>% sprintf("%.2e", .)
        
        df_plot %>% 
          dplyr::arrange(normalized_expr) -> df_plot
        upper = quantile(df_plot$normalized_expr, 0.99, na.rm = TRUE)
        df_plot %>% 
          ggplot(aes(x = x, y = y, color = normalized_expr, size = normalized_expr)) + 
          #              geom_point(size = 1) + 
          geom_point(alpha = 0.8) + 
          scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, upper), oob = scales::squish) + 
          scale_size_continuous(range = c(2, 5)) +
          theme_void() + 
          labs(color = "Normalized\nExpression") +
          ggtitle(sprintf("%s, %s: p.adj=%s", sample_, gene_, pval)) + 
          theme(legend.position = "none")
      }) %>% patchwork::wrap_plots(., ncol = 3) -> p
  }) -> p_list
pdf("~/Downloads/p_list.pdf", width = 14.56, height = 7.61)
p_list
dev.off()
#Etv5, Adra1b, Calnt18, Fos, Chrm3, Ptprf


c("Prkaa2", "Nrarp") %in% gene_list$CT8_only
c("Chst8", "Rab3b") %in% gene_list$CT20_only
c("Vip", "Gpr101") %in% gene_list$Both

gene_list %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(genes_, group_){
    df_combinedPval_all_samples_sparkx %>% dplyr::filter(gene %in% genes_) -> df_
    df_ %>% dplyr::select(gene, `CT-8_cauchy_padj`, `CT-20_cauchy_padj`) -> df_
    if (group_ == "CT8_only"){
      df_ %>% dplyr::arrange(`CT-8_cauchy_padj`) %>% dplyr::filter(`CT-20_cauchy_padj` > 0.3) -> df_
    }else if (group_ == "CT20_only"){
      df_ %>% dplyr::arrange(`CT-20_cauchy_padj`) %>% dplyr::filter(`CT-8_cauchy_padj` > 0.3) -> df_
    }
    
    df_ = head(df_, 20)
    return(df_$gene)
  }) -> gene_list_2

gene_list_1 = list(
  CT8 = c("Prkaa2", "Nrarp"), 
  CT20 = c("Chst8", "Rab3b"),
  Both = c("Vip", "Gpr101")
)

gene_list_1 %>% 
  #  .[1] %>% 
  purrr::map(function(genes_){
    genes_ %>% 
      "names<-"(.,.) %>% 
      #      .[1] %>% 
      purrr::map(function(gene_){
        c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
          "names<-"(.,.) %>% 
          purrr::map(function(sample_){
            norm_counts_matrix[gene_, ] %>% 
              t() %>% 
              as.data.frame() %>% "colnames<-"(., "normalized_expr") %>% 
              rownames_to_column("cell") -> df_
            xy_coordinates[[sample_]] -> df_xy
            inner_join(x = df_, y = df_xy, by = "cell") -> df_plot
            
            df_combinedPval_all_samples_sparkx %>% dplyr::filter(gene == gene_) %>% 
              dplyr::select(matches(sample_)) %>% .[1,2] %>% sprintf("%.2e", .) -> pval
            
#            pval = sparkx_res_list[[sample_]]$res_mtest[gene_, ]$adjustedPval %>% sprintf("%.2e", .)
            
            df_plot %>% 
              dplyr::arrange(normalized_expr) -> df_plot
            upper = quantile(df_plot$normalized_expr, 0.99, na.rm = TRUE)
            df_plot %>% 
              ggplot(aes(x = x, y = y, color = normalized_expr, size = normalized_expr)) + 
#              geom_point(size = 1) + 
              geom_point(alpha = 0.8) + 
              scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, upper), oob = scales::squish) + 
              scale_size_continuous(range = c(2, 5)) +
              theme_void() + 
              labs(color = "Normalized\nExpression") +
              ggtitle(sprintf("%s, %s: p.adj=%s", sample_, gene_, pval)) + 
              theme(legend.position = "none")
          }) %>% patchwork::wrap_plots(., ncol = 3) -> p
      })
  }) -> p_list_1

p_list_1$CT8$Prkaa2 #Figure 3C
p_list_1$CT8$His2h2bb #Figure 3C
  
p_list_1$CT20$Chst8 #Figure 3C
p_list_1$CT20$Rab3b #Figure 3C

p_list_1$Both$Vip #Figure 3B
p_list_1$Both$Gpr101 #Figure 3B

p_list_1_1 = c(p_list_1$CT8, p_list_1$CT20, p_list_1$Both)
#pdf("~/Downloads/SPARK_X_CT8_CT20_Both_genes.pdf", width = 13.06, height = 6.23)
pdf("~/Downloads/SPARK_X_CT8_CT20_Both_genes_v2.pdf", width = 9.8, height = 5.57)
p_list_1_1
dev.off()

gene_list_2 %>% 
  #  .[1] %>% 
  purrr::map(function(genes_){
    genes_ %>% 
      "names<-"(.,.) %>% 
      #      .[1] %>% 
      purrr::map(function(gene_){
        c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
          "names<-"(.,.) %>% 
          purrr::map(function(sample_){
            norm_counts_matrix[gene_, ] %>% 
              t() %>% 
              as.data.frame() %>% "colnames<-"(., "normalized_expr") %>% 
              rownames_to_column("cell") -> df_
            xy_coordinates[[sample_]] -> df_xy
            inner_join(x = df_, y = df_xy, by = "cell") -> df_plot
            
            df_combinedPval_all_samples_sparkx %>% dplyr::filter(gene == gene_) %>% 
              dplyr::select(matches(sample_)) %>% .[1,2] %>% sprintf("%.2e", .) -> pval
            
            #            pval = sparkx_res_list[[sample_]]$res_mtest[gene_, ]$adjustedPval %>% sprintf("%.2e", .)
            
            df_plot %>% 
              dplyr::arrange(normalized_expr) -> df_plot
            upper = quantile(df_plot$normalized_expr, 0.99, na.rm = TRUE)
            df_plot %>% 
              ggplot(aes(x = x, y = y, color = normalized_expr, size = normalized_expr)) + 
              #              geom_point(size = 1) + 
              geom_point(alpha = 0.8) + 
              scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, upper), oob = scales::squish) + 
              scale_size_continuous(range = c(2, 5)) +
              theme_void() + 
              labs(color = "Normalized\nExpression") +
              ggtitle(sprintf("%s, %s: p.adj=%s", sample_, gene_, pval)) + 
              theme(legend.position = "none")
          }) %>% patchwork::wrap_plots(., ncol = 3) -> p
      })
  }) -> p_list_2

p_list_2$CT8_only$Hist2h2bb
p_list_2$CT8_only$Rnaset2a
p_list_2$CT8_only$Klc1
p_list_2$CT8_only$Gad1
