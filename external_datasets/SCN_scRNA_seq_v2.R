library(tidyverse)

list.files("~/Dropbox/Spatial_transcriptomics_new/external_data/Wen_mouse_SCN_circadian/GSE117295_RAW", full.names = T) %>% 
  .[grep("CT", .)] %>% 
  "names<-"(., gsub("/.+/.+(CT\\d+)\\.csv\\.gz", "\\1", .)) %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(file_, CT_){
    print(file_)
    read.csv(file = file_, row.names = 1, header = T) -> df_
  }) -> mat_list

mat_list %>% 
  purrr::map(function(df_){
    df_ %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6})) -> df_
  }) -> mat_norm_list
rm(mat_list);gc()

mat_norm_list$CT14
mat_norm_list %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_,CT_){
    rowMeans(df_) %>% as.data.frame() %>% "colnames<-"(.,CT_) %>% rownames_to_column("gene")
  }) %>% purrr::reduce(., inner_join, by = "gene") -> mat_norm_mean

se <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
}

mat_norm_list %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_,CT_){
    apply(df_, 1, se) %>% as.data.frame() %>% "colnames<-"(.,CT_) %>% rownames_to_column("gene")
  }) %>% purrr::reduce(., inner_join, by = "gene") -> mat_norm_se

mat_norm_mean_1 = mat_norm_mean %>% column_to_rownames("gene")
mat_norm_mean_1 %>% colnames() %>% gsub("CT(\\d+)", "\\1", .) %>% as.numeric() -> ZT_
colnames(mat_norm_mean_1) %>% gsub("CT", "ZT", .) %>% sprintf("%s_REP1", .) -> colnames(mat_norm_mean_1)
mat_norm_mean_1 %>% 
  rownames_to_column("Gene") %>% 
  write.csv(., "~/Downloads/temp.csv", row.names = F, quote = F)

library(MetaCycle)
library(HarmonicRegression)

JTK.i=meta2d(infile = "~/Downloads/temp.csv", filestyle = "csv",
             outputFile = FALSE, timepoints = ZT_, cycMethod = c("JTK"),
             maxper = 24, minper = 20, combinePvalue = 'bonferroni', parallelize = T, nCores = 8)$meta

save(mat_norm_mean, mat_norm_se, JTK.i, file = "~/Downloads/scRNA_dat.rda")


c("Csnk1e", "Calb1", "Slc39a6", "Lhx1", "Rgs16", "Grin2c", "Vip", "Vipr2") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    mat_norm_mean %>% dplyr::filter(gene == gene_) %>% 
      dplyr::select(-gene) %>% 
      pivot_longer(cols = everything(), names_to = "CT", values_to = "norm_expr") -> df_expr
    
    mat_norm_se %>% dplyr::filter(gene == gene_) %>% 
      dplyr::select(-gene) %>% 
      pivot_longer(cols = everything(), names_to = "CT", values_to = "se") -> df_se
    
    full_join(x = df_expr, y = df_se, by = "CT") %>% 
      dplyr::mutate(CT = gsub("CT(\\d+)", "\\1", CT) %>% as.numeric()) -> df_plot
    
    night_df <- data.frame(
      xmin = c(14, 36),
      xmax = c(24, 48)
    )
    
    
    df_plot %>% 
      ggplot(aes(x = CT, y = norm_expr)) + 
      geom_rect(
        data = night_df,
        aes(xmin = xmin, xmax = xmax),
        ymin = -Inf,
        ymax = Inf,
        inherit.aes = FALSE,
        fill = "grey40",
        alpha = 0.15,
      ) +
      geom_line() + 
      geom_ribbon(aes(ymin = norm_expr - se, ymax = norm_expr + se), 
                  alpha = 0.2) + 
      theme_classic() + 
      ylab("Normalized expression") -> p
    
    JTK.i %>% dplyr::filter(CycID == gene_) %>% .$JTK_pvalue %>% {round(., 1)} -> pval_
    
    p + ggtitle(
      sprintf("(SCN scRNA-seq) %s: JTK pval = %s", gene_, pval_)
    )
    
    p + scale_x_continuous(breaks = seq(14,58,4)) + ggtitle(gene_)
  }) -> p_list

compare_day_night = function(gene_, day_CT, night_CT){
  mat_norm_list[[day_CT]][gene_, ] %>% as.numeric() -> day_expr
  mat_norm_list[[night_CT]][gene_, ] %>% as.numeric() -> night_expr
  t.test(x = day_expr, y = night_expr)
}

load("~/Downloads/df_spatiotemporal_all_2.rda")
compare_day_night_adj = function(day_CT, night_CT){
  c(rownames(mat_norm_list[[day_CT]]), 
    rownames(mat_norm_list[[night_CT]]),
    df_spatiotemporal_all_2$gene) %>% table() %>% 
    {.[. == 3]} %>% names() -> genes_
  i = 0
  genes_ %>% 
#    .[1:10] %>% 
    purrr::map(function(gene__){
      i <<- i + 1
      print(i)
      mat_norm_list[[day_CT]][gene__, ] %>% as.numeric() -> day_expr
      mat_norm_list[[night_CT]][gene__, ] %>% as.numeric() -> night_expr
      t.test(x = day_expr, y = night_expr)$p.value -> pval
      data.frame(gene = gene__, pval = pval)
    }) %>% do.call(rbind, .) -> df_
  df_ %>% dplyr::mutate(BH = p.adjust(pval, method = "BH"))
}

compare_day_night_adj(day_CT = "CT30", night_CT = "CT18") -> df_day_night_CT30_CT18
compare_day_night_adj(day_CT = "CT54", night_CT = "CT42") -> df_day_night_CT54_CT42
compare_day_night_adj(day_CT = "CT30", night_CT = "CT42") -> df_day_night_CT30_CT42

mat_norm_mean %>% pivot_longer(cols = -gene, names_to = "CT", values_to = "norm_expr") %>% 
  dplyr::mutate(day_night = case_when(
    CT %in% c("CT30", "CT54") ~ "day",
    CT %in% c("CT18", "CT42") ~ "night",
    TRUE ~ NA
  )) %>% dplyr::filter(!is.na(day_night)) -> mat_norm_mean_long

save(df_day_night_CT30_CT18, df_day_night_CT54_CT42, df_day_night_CT30_CT42, mat_norm_mean_long, file = "~/Downloads/df_day_night.rda")

compare_day_night(gene_ = "Calb1", day_CT = "CT30", night_CT = "CT18")
df_day_night_CT30_CT18 %>% dplyr::filter(gene == "Calb1")
compare_day_night(gene_ = "Csnk1e", day_CT = "CT30", night_CT = "CT18")
df_day_night_CT30_CT18 %>% dplyr::filter(gene == "Csnk1e")
compare_day_night(gene_ = "Slc39a6", day_CT = "CT54", night_CT = "CT42")
df_day_night_CT54_CT42 %>% dplyr::filter(gene == "Slc39a6")

compare_day_night(gene_ = "Lhx1", day_CT = "CT30", night_CT = "CT18")
df_day_night_CT30_CT18 %>% dplyr::filter(gene == "Lhx1")
compare_day_night(gene_ = "Rgs16", day_CT = "CT54", night_CT = "CT42")
df_day_night_CT54_CT42 %>% dplyr::filter(gene == "Rgs16")
compare_day_night(gene_ = "Il33", day_CT = "CT30", night_CT = "CT18")
df_day_night_CT30_CT18 %>% dplyr::filter(gene == "Il33")

redraw_figure = function(data, CT_start, CT_end, gene_){
  data %>% dplyr::filter(CT <= CT_end, CT >= CT_start) -> df_plot
  night_df <- data.frame(
    xmin = c(CT_start),
    xmax = c((CT_start+CT_end)/2)
  )
  
  df_plot %>% 
    ggplot(aes(x = CT, y = norm_expr)) + 
    geom_rect(
      data = night_df,
      aes(xmin = xmin, xmax = xmax),
      ymin = -Inf,
      ymax = Inf,
      inherit.aes = FALSE,
      fill = "grey40",
      alpha = 0.15,
    ) +
    geom_line() + 
    geom_ribbon(aes(ymin = norm_expr - se, ymax = norm_expr + se), 
                alpha = 0.2) + 
    geom_point() + 
    theme_classic() + 
    ylab("Normalized expression") -> p
  
  JTK.i %>% dplyr::filter(CycID == gene_) %>% .$JTK_pvalue %>% {round(., 1)} -> pval_
  
  p + ggtitle(
    sprintf("(SCN scRNA-seq) %s: JTK pval = %s", gene_, pval_)
  )
  
  p + scale_x_continuous(breaks = seq(14,58,4)) + ggtitle(gene_) -> p1
}

compare_day_night_boxplot = function(gene_, day_CT, night_CT){
  mat_norm_list[[day_CT]][gene_, ] %>% as.numeric() -> day_expr
  mat_norm_list[[night_CT]][gene_, ] %>% as.numeric() -> night_expr
  t.test(x = day_expr, y = night_expr)$p.value %>% sprintf("%.2e", .) -> pval
  
  rbind(
    data.frame(CT = day_CT, expr = day_expr),
    data.frame(CT = night_CT, expr = night_expr) 
  ) -> df_plot
  
  df_plot %>% 
#    dplyr::filter(expr != 0) %>% 
    ggplot(aes(x = CT, y = expr)) + 
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.12,
      linewidth = 0.8
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 3
    ) + 
    theme_classic() + 
    ggtitle(sprintf("P=%s", pval))
}

p_list$Calb1
redraw_figure(data = p_list$Calb1$data, CT_start = 38, CT_end = 58, gene_ = "Calb1") -> p_Calb1_1
compare_day_night_boxplot(gene_ = "Calb1", day_CT = "CT54", night_CT = "CT42") -> p_Calb1_2
patchwork::wrap_plots(p_Calb1_1, p_Calb1_2, ncol = 2, widths = c(3,1)) -> p_Calb1_3

redraw_figure(data = p_list$Calb1$data, CT_start = 14, CT_end = 34, gene_ = "Calb1") -> p_Calb1_1
compare_day_night_boxplot(gene_ = "Calb1", day_CT = "CT30", night_CT = "CT18") -> p_Calb1_2
patchwork::wrap_plots(p_Calb1_1, p_Calb1_2, ncol = 2, widths = c(3,1)) -> p_Calb1_3

p_list$Csnk1e
redraw_figure(data = p_list$Csnk1e$data, CT_start = 14, CT_end = 34, gene_ = "Csnk1e") -> p_Csnk1e_1
compare_day_night_boxplot(gene_ = "Csnk1e", day_CT = "CT30", night_CT = "CT18") -> p_Csnk1e_2
patchwork::wrap_plots(p_Csnk1e_1, p_Csnk1e_2, ncol = 2, widths = c(3,1)) -> p_Csnk1e_3

p_list$Slc39a6
redraw_figure(data = p_list$Slc39a6$data, CT_start = 38, CT_end = 58, gene_ = "Slc39a6") -> p_Slc39a6_1
compare_day_night_boxplot(gene_ = "Slc39a6", day_CT = "CT54", night_CT = "CT42") -> p_Slc39a6_2
patchwork::wrap_plots(p_Slc39a6_1, p_Slc39a6_2, ncol = 2, widths = c(3,1)) -> p_Slc39a6_3

p_list$Lhx1
redraw_figure(data = p_list$Lhx1$data, CT_start = 14, CT_end = 34, gene_ = "Lhx") -> p_Lhx_1
compare_day_night_boxplot(gene_ = "Lhx1", day_CT = "CT30", night_CT = "CT18") -> p_Lhx_2
patchwork::wrap_plots(p_Lhx_1, p_Lhx_2, ncol = 2, widths = c(3,1)) -> p_Lhx_3

p_list$Rgs16
redraw_figure(data = p_list$Rgs16$data, CT_start = 38, CT_end = 58, gene_ = "Rgs16") -> p_Rgs16_1
compare_day_night_boxplot(gene_ = "Rgs16", day_CT = "CT54", night_CT = "CT42") -> p_Rgs16_2
patchwork::wrap_plots(p_Rgs16_1, p_Rgs16_2, ncol = 2, widths = c(3,1)) -> p_Rgs16_3


p_list$Grin2c
redraw_figure(data = p_list$Grin2c$data, CT_start = 14, CT_end = 34, gene_ = "Grin2c") -> p_Grin2c_1
compare_day_night_boxplot(gene_ = "Grin2c", day_CT = "CT30", night_CT = "CT18") -> p_Grin2c_2
patchwork::wrap_plots(p_Grin2c_1, p_Grin2c_2, ncol = 2, widths = c(3,1)) -> p_Grin2c_3

ggsave(plot = p_Csnk1e_3, width = 6.92, height = 4.17, filename = "~/Downloads/p_Csnk1e_3.pdf")
ggsave(plot = p_Calb1_3, width = 6.92, height = 4.17, filename = "~/Downloads/p_Calb1_3.pdf")
ggsave(plot = p_Slc39a6_3, width = 6.92, height = 4.17, filename = "~/Downloads/p_Slc39a6_3.pdf")
ggsave(plot = p_Lhx_3, width = 6.92, height = 4.17, filename = "~/Downloads/p_Lhx1_3.pdf")
ggsave(plot = p_Rgs16_3, width = 6.92, height = 4.17, filename = "~/Downloads/p_Rgs16_3.pdf")
ggsave(plot = p_Grin2c_3, width = 6.92, height = 4.17, filename = "~/Downloads/p_Grin2c_3.pdf")

pdf("~/Downloads/scRNA_circadian_testing.pdf", height = 4.02, width = 4.06)
p_list
dev.off()

###
load("~/Downloads/scRNA_dat.rda")
colnames(mat_norm_mean)[-1] %>% gsub("CT(\\d+)", "\\1", .) %>% as.numeric() %>% {ifelse(.>24, .-24, .)} %>% {ifelse(.>24, .-24, .)} %>% 
  {c("gene", sprintf("CT%s", .))} -> new_colname

day_night_compare = function(gene_){
  mat_norm_mean %>% dplyr::filter(gene == gene_) %>% "colnames<-"(., new_colname) %>% 
    pivot_longer(-gene, names_to = "CT", values_to = "expr") -> df_
#    df_ %>% dplyr::filter(CT %in% c("CT6")) %>% dplyr::mutate(group = "day") -> df_day
#    df_ %>% dplyr::filter(CT %in% c("CT18")) %>% dplyr::mutate(group = "night") -> df_night
    df_ %>% dplyr::filter(CT %in% c("CT10")) %>% dplyr::mutate(group = "day") -> df_day
    df_ %>% dplyr::filter(CT %in% c("CT22")) %>% dplyr::mutate(group = "night") -> df_night
    t.test(x = df_day$expr, y = df_night$expr)
    
    rbind(df_day, df_night) %>% 
      ggplot(aes(x = group, y = expr, color = group)) + 
      geom_point() + 
      theme_classic() + 
      ggtitle(gene_)
}

day_night_compare(gene = "Csnk1e")
day_night_compare(gene = "Calb1")
day_night_compare(gene = "Slc39a6")
day_night_compare(gene = "Lhx1")
day_night_compare(gene = "Rgs16")
day_night_compare(gene = "Grin2c")

JTK.i %>% dplyr::filter(CycID == "Lhx1")
JTK.i %>% dplyr::filter(CycID == "Rgs16")
JTK.i %>% dplyr::filter(CycID == "Grin2c")

JTK.i %>% dplyr::filter(CycID == "Csnk1e")
JTK.i %>% dplyr::filter(CycID == "Calb1")
JTK.i %>% dplyr::filter(CycID == "Slc39a6")


### New analysis----
list.files("~/Dropbox/Spatial_transcriptomics_new/external_data/Wen_mouse_SCN_circadian/GSE117295_RAW", full.names = T) %>% 
  .[grep("CT", .)] %>% 
  "names<-"(., gsub("/.+/.+(CT\\d+)\\.csv\\.gz", "\\1", .)) %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(file_, CT_){
    print(file_)
    read.csv(file = file_, row.names = 1, header = T) -> df_
  }) -> mat_list

barcode_list = list()
list.files("~/Dropbox/singulomics/github_rda/GSE117295_RAW_Mouse_SCN/cluster_vs_barcode", pattern = "\\.tsv", full.names = T) %>% .[2:3] %>% 
  purrr::map(function(files_){
    read.table(file = files_, header = T, sep = "\t", stringsAsFactors = F) -> df_
    colnames(df_) <- c("barcode", "celltypes")
    df_$celltypes %>% unique() %>% sort() -> celltypes_
    celltypes_ %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(celltype_){
        df_ %>% dplyr::filter(celltypes == celltype_) %>% .$barcode -> barcodes_
        barcode_list[[celltype_]] <<- barcodes_
      })
  })

barcode_list_1 = list()
barcode_list[c("Neurons1", "Neurons2", "Neurons3", "Neurons4", "Neurons5", "Neurons6")] %>% unlist() %>% unname() %>% unique() -> barcode_list_1[["Neurons"]]

barcode_list_1 %>% 
#  .[1] %>% 
  purrr::map(function(barcodes_){
    mat_list %>% 
#      .[1] %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(df_,CT_){
        intersect(colnames(df_), barcodes_) -> barcodes_
        df_ = df_[, barcodes_]
        df_ %>% rowSums() %>% as.data.frame() %>% "colnames<-"(., CT_) -> df_
        df_ %>% rownames_to_column("gene") -> df_
      }) %>% purrr::reduce(., inner_join, by = "gene") -> df_
    df_ %>% column_to_rownames("gene") -> df_
    return(df_)
  }) -> scn_celltype_raw_counts

scn_celltype_raw_counts %>% 
  purrr::map(function(df_){
    df_ %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6}))
  }) -> scn_celltype_norm_counts

barcode_list_2 = list()
barcode_list$`vip/nms` %>% unlist() %>% unname() %>% unique() -> barcode_list_2[["Vip_nms"]]

barcode_list_2 %>% 
  #  .[1] %>% 
  purrr::map(function(barcodes_){
    mat_list %>% 
      #      .[1] %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(df_,CT_){
        intersect(colnames(df_), barcodes_) -> barcodes_
        df_ = df_[, barcodes_]
        df_ %>% rowSums() %>% as.data.frame() %>% "colnames<-"(., CT_) -> df_
        df_ %>% rownames_to_column("gene") -> df_
      }) %>% purrr::reduce(., inner_join, by = "gene") -> df_
    df_ %>% column_to_rownames("gene") -> df_
    return(df_)
  }) -> scn_vip_nms_raw_counts

scn_vip_nms_raw_counts %>% 
  purrr::map(function(df_){
    df_ %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6}))
  }) -> scn_vip_nms_norm_counts

barcode_list_3 = list()
barcode_list$`avp/nms` %>% unlist() %>% unname() %>% unique() -> barcode_list_3[["Avp_nms"]]

barcode_list_3 %>% 
  #  .[1] %>% 
  purrr::map(function(barcodes_){
    mat_list %>% 
      #      .[1] %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(df_,CT_){
        intersect(colnames(df_), barcodes_) -> barcodes_
        df_ = df_[, barcodes_]
        df_ %>% rowSums() %>% as.data.frame() %>% "colnames<-"(., CT_) -> df_
        df_ %>% rownames_to_column("gene") -> df_
      }) %>% purrr::reduce(., inner_join, by = "gene") -> df_
    df_ %>% column_to_rownames("gene") -> df_
    return(df_)
  }) -> scn_avp_nms_raw_counts

scn_avp_nms_raw_counts %>% 
  purrr::map(function(df_){
    df_ %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6}))
  }) -> scn_avp_nms_norm_counts


library(MetaCycle)
library(HarmonicRegression)

scn_celltype_norm_counts$Neurons %>% "colnames<-"(., gsub("CT", "ZT", colnames(.))) %>% 
  {
    df_ = .
    ZT_ = colnames(df_) %>% gsub("ZT(\\d+)", "\\1", .) %>% as.numeric()
    sprintf("%s_REP1", colnames(df_)) -> colnames(df_)
    df_ %>% rownames_to_column("gene") %>% 
      write.csv(., "~/Downloads/temp.csv", row.names = F, quote = F)
    JTK.i=meta2d(infile = "~/Downloads/temp.csv", filestyle = "csv",
                 outputFile = FALSE, timepoints = ZT_, cycMethod = c("JTK"),
                 maxper = 24, minper = 24, combinePvalue = 'bonferroni', parallelize = T, nCores = 8)
    JTK.i
  } -> JTK_results

scn_vip_nms_norm_counts$Vip_nms %>% "colnames<-"(., gsub("CT", "ZT", colnames(.))) %>% 
  {
    df_ = .
    ZT_ = colnames(df_) %>% gsub("ZT(\\d+)", "\\1", .) %>% as.numeric()
    sprintf("%s_REP1", colnames(df_)) -> colnames(df_)
    df_ %>% rownames_to_column("gene") %>% 
      write.csv(., "~/Downloads/temp.csv", row.names = F, quote = F)
    JTK.i=meta2d(infile = "~/Downloads/temp.csv", filestyle = "csv",
                 outputFile = FALSE, timepoints = ZT_, cycMethod = c("JTK"),
                 maxper = 24, minper = 24, combinePvalue = 'bonferroni', parallelize = T, nCores = 8)
    JTK.i
  } -> JTK_results_vip_nms

scn_avp_nms_norm_counts$Avp_nms %>% "colnames<-"(., gsub("CT", "ZT", colnames(.))) %>% 
  {
    df_ = .
    ZT_ = colnames(df_) %>% gsub("ZT(\\d+)", "\\1", .) %>% as.numeric()
    sprintf("%s_REP1", colnames(df_)) -> colnames(df_)
    df_ %>% rownames_to_column("gene") %>% 
      write.csv(., "~/Downloads/temp.csv", row.names = F, quote = F)
    JTK.i=meta2d(infile = "~/Downloads/temp.csv", filestyle = "csv",
                 outputFile = FALSE, timepoints = ZT_, cycMethod = c("JTK"),
                 maxper = 24, minper = 24, combinePvalue = 'bonferroni', parallelize = T, nCores = 8)
    JTK.i
  } -> JTK_results_avp_nms

save(barcode_list, barcode_list_1, JTK_results, scn_celltype_norm_counts, scn_celltype_raw_counts, file = "~/Downloads/SCN_circadian_analysis.rda")

JTK_results$meta %>% dim()
JTK_results$JTK %>% dplyr::filter(CycID %in% c("Arntl", "Per1", "Per2", "Nr1d1", "Nr1d2", "Dbp", "Tef", "Ciart"))

JTK_results$JTK %>% dplyr::filter(CycID == "Gad1")
JTK_results$JTK %>% dplyr::filter(CycID == "Calb1")
JTK_results$JTK %>% dplyr::filter(CycID == "Slc39a6")

JTK_results$JTK %>% dplyr::filter(CycID == "Lhx1")
JTK_results$JTK %>% dplyr::filter(CycID == "Il33")
JTK_results$JTK %>% dplyr::filter(CycID == "Rgs16")

draw_lineplot = function(gene){
  scn_celltype_norm_counts$Neurons[gene, ] %>% 
    pivot_longer(cols = everything(), names_to = "CT", values_to = "TPM") %>% 
    dplyr::mutate(CT = gsub("CT", "", CT) %>% as.numeric()) -> df_plot
  df_plot %>% 
    ggplot(aes(x = CT, y = TPM)) + 
    geom_rect(data = data.frame(xmin = c(14, 36), xmax = c(24, 48), ymin = -Inf, ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey", alpha = 0.5, inherit.aes = FALSE) + 
    geom_line() + 
    geom_point() + 
    scale_x_continuous(breaks = seq(14,58,4)) + 
    theme_classic() -> p
  JTK_results$JTK %>% dplyr::filter(CycID == gene) %>% .$ADJ.P %>% sprintf("%.2e", .) -> p.adj
  p + ggtitle(sprintf("%s (adj.p = %s)", gene, p.adj))
}

draw_lineplot(gene = "Gad1") #Figure 3F
draw_lineplot(gene = "Calb1") #Figure 3G
draw_lineplot(gene = "Slc39a6") #Figure 3H
draw_lineplot(gene = "Per2")
draw_lineplot(gene = "Rgs16")
draw_lineplot(gene = "Lhx1")

c("Lhx1") %>% 
  purrr::map(function(gene_){
    scn_vip_nms_norm_counts$Vip_nms[gene_, ] %>% 
      pivot_longer(cols = everything(), names_to = "CT", values_to = "TPM") %>% 
      dplyr::mutate(CT = gsub("CT", "", CT) %>% as.numeric()) %>% dplyr::mutate(group = "vip_nms") -> df_plot 
    scn_avp_nms_norm_counts$Avp_nms[gene_, ] %>% 
      pivot_longer(cols = everything(), names_to = "CT", values_to = "TPM") %>% 
      dplyr::mutate(CT = gsub("CT", "", CT) %>% as.numeric()) %>% dplyr::mutate(group = "avp_nms") -> df_plot_1
    rbind(df_plot, df_plot_1) -> df_plot
    df_plot %>% 
      ggplot(aes(x = CT, y = TPM, group = group, color = group)) + 
      geom_rect(data = data.frame(xmin = c(14, 36), xmax = c(24, 48), ymin = -Inf, ymax = Inf),
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "grey", alpha = 0.5, inherit.aes = FALSE) + 
      geom_line() + 
      geom_point() + 
      scale_x_continuous(breaks = seq(14,58,4)) + 
      theme_classic() -> p
    JTK_results_vip_nms$JTK %>% dplyr::filter(CycID == gene_) %>% .$ADJ.P %>% sprintf("%.2e", .) -> p.adj_vip_nms
    JTK_results_avp_nms$JTK %>% dplyr::filter(CycID == gene_) %>% .$ADJ.P %>% sprintf("%.2e", .) -> p.adj_avp_nms
    p + ggtitle(sprintf("%s (vip_nms adj.p = %s) (avp_nms adj.p = %s)", gene_, p.adj_vip_nms, p.adj_avp_nms))
  }) #Figure 4E

c("Rgs16") %>% 
  purrr::map(function(gene_){
    scn_vip_nms_norm_counts$Vip_nms[gene_, ] %>% 
      pivot_longer(cols = everything(), names_to = "CT", values_to = "TPM") %>% 
      dplyr::mutate(CT = gsub("CT", "", CT) %>% as.numeric()) %>% dplyr::mutate(group = "vip_nms") -> df_plot 
    scn_avp_nms_norm_counts$Avp_nms[gene_, ] %>% 
      pivot_longer(cols = everything(), names_to = "CT", values_to = "TPM") %>% 
      dplyr::mutate(CT = gsub("CT", "", CT) %>% as.numeric()) %>% dplyr::mutate(group = "avp_nms") -> df_plot_1
    rbind(df_plot, df_plot_1) -> df_plot
    df_plot %>% 
      ggplot(aes(x = CT, y = TPM, group = group, color = group)) + 
      geom_rect(data = data.frame(xmin = c(14, 36), xmax = c(24, 48), ymin = -Inf, ymax = Inf),
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "grey", alpha = 0.5, inherit.aes = FALSE) + 
      geom_line() + 
      geom_point() + 
      scale_x_continuous(breaks = seq(14,58,4)) + 
      theme_classic() -> p
    JTK_results_vip_nms$JTK %>% dplyr::filter(CycID == gene_) %>% .$ADJ.P %>% sprintf("%.2e", .) -> p.adj_vip_nms
    JTK_results_avp_nms$JTK %>% dplyr::filter(CycID == gene_) %>% .$ADJ.P %>% sprintf("%.2e", .) -> p.adj_avp_nms
    p + ggtitle(sprintf("%s (vip_nms adj.p = %s) (avp_nms adj.p = %s)", gene_, p.adj_vip_nms, p.adj_avp_nms))
  }) #Figure 4F

c("Per2") %>% 
  purrr::map(function(gene_){
    scn_vip_nms_norm_counts$Vip_nms[gene_, ] %>% 
      pivot_longer(cols = everything(), names_to = "CT", values_to = "TPM") %>% 
      dplyr::mutate(CT = gsub("CT", "", CT) %>% as.numeric()) %>% dplyr::mutate(group = "vip_nms") -> df_plot 
    scn_avp_nms_norm_counts$Avp_nms[gene_, ] %>% 
      pivot_longer(cols = everything(), names_to = "CT", values_to = "TPM") %>% 
      dplyr::mutate(CT = gsub("CT", "", CT) %>% as.numeric()) %>% dplyr::mutate(group = "avp_nms") -> df_plot_1
    rbind(df_plot, df_plot_1) -> df_plot
    df_plot %>% 
      ggplot(aes(x = CT, y = TPM, group = group, color = group)) + 
      geom_rect(data = data.frame(xmin = c(14, 36), xmax = c(24, 48), ymin = -Inf, ymax = Inf),
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "grey", alpha = 0.5, inherit.aes = FALSE) + 
      geom_line() + 
      geom_point() + 
      scale_x_continuous(breaks = seq(14,58,4)) + 
      theme_classic() -> p
    JTK_results_vip_nms$JTK %>% dplyr::filter(CycID == gene_) %>% .$ADJ.P %>% sprintf("%.2e", .) -> p.adj_vip_nms
    JTK_results_avp_nms$JTK %>% dplyr::filter(CycID == gene_) %>% .$ADJ.P %>% sprintf("%.2e", .) -> p.adj_avp_nms
    p + ggtitle(sprintf("%s (vip_nms adj.p = %s) (avp_nms adj.p = %s)", gene_, p.adj_vip_nms, p.adj_avp_nms))
  }) #Figure 4C

save(JTK_results_vip_nms, JTK_results_avp_nms, scn_avp_nms_norm_counts, 
     scn_avp_nms_raw_counts, scn_vip_nms_norm_counts, scn_vip_nms_raw_counts, 
     file = "~/Downloads/SCN_circadian_analysis_1.rda")
