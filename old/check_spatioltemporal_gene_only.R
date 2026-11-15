library(tidyverse)
load("~/Downloads/df_spatiotemporal_all_2.rda")
df_spatiotemporal_all_2 %>% dplyr::filter(spatiotemp.padj_spatiotemporal < 0.05) %>% 
  dplyr::filter(p_adj_TVG > 0.05, `CT-8_cauchy_padj_sparkx` > 0.05, `CT-20_cauchy_padj_sparkx` > 0.05) %>% 
  dplyr::arrange(spatiotemp.padj_spatiotemporal) %>% 
  .$gene -> spatiotemporal_gene_only

source("~/Dropbox/Spatial_transcriptomics_new/github/segmented/functions.R")
spatiotemporal_gene_only %>% 
#  "names<-"(.,.) %>% 
#  .[1] %>% 
  purrr::map(function(gene_){
    norm_counts_matrix[gene_, ] %>% {.!=0} %>% sum() %>% {./ncol(norm_counts_matrix)} %>% {.*100} -> percent_
    data.frame(gene = gene_, percent = percent_)
  }) %>% do.call(rbind, .) -> df_

df_$percent %>% density() %>% plot()

apply(norm_counts_matrix, 1, function(x){(sum(x!=0)/length(x))*100}) -> gene_percent
gene_percent %>% as.data.frame() %>% rownames_to_column("gene") %>% "colnames<-"(., c("gene", "percent")) -> gene_percent_df
gene_percent_df %>% dplyr::filter(percent > 10) %>% dim()

df_spatiotemporal_all_2 %>% dplyr::filter(`CT-8_cauchy_padj_sparkx` < 0.05) %>% .$gene -> CT8_SVG
df_spatiotemporal_all_2 %>% dplyr::filter(`CT-20_cauchy_padj_sparkx` < 0.05) %>% .$gene -> CT20_SVG
df_spatiotemporal_all_2 %>% dplyr::filter(p_adj_TVG < 0.05) %>% .$gene -> TVG
df_spatiotemporal_all_2 %>% dplyr::filter(spatiotemp.padj_spatiotemporal < 0.05) %>% .$gene -> STVG
gene_percent_df %>% dplyr::filter(percent > 5) %>% .$gene -> HESG

length(TVG) #436

CT8_SVG %in% HESG %>% sum() #819
length(CT8_SVG) #896
CT8_SVG[(CT8_SVG %in% HESG)] -> CT8_SVG_1

CT20_SVG %in% HESG %>% sum() #737
length(CT20_SVG) #844
CT20_SVG[(CT20_SVG %in% HESG)] -> CT20_SVG_1

STVG %in% HESG %>% sum() #176
length(STVG) #521
STVG[(STVG %in% HESG)] -> STVG_1

list(
  CT8_SVG = CT8_SVG_1,
  CT20_SVG = CT20_SVG_1,
  TVG = TVG,
  STVG = STVG_1
) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::mutate(across(!matches("key"), function(x){x %>% as.integer()})) %>%
  as.data.frame() %>%
  UpSetR::upset(sets = rev(c("TVG", "CT8_SVG", "CT20_SVG", "STVG")), keep.order = T, order.by = "freq"#,
#                main.bar.color = color_list[[direction_]],
#                sets.bar.color = color_list[[direction_]], 
#                matrix.color = color_list[[direction_]])
) -> p_upset #Figure 4A


#Find gene change it's spatial location across times
list(
  CT8_SVG = CT8_SVG_1,
  CT20_SVG = CT20_SVG_1,
  TVG = TVG,
  STVG = STVG_1
) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::filter(TVG == F, CT8_SVG == T, CT20_SVG == T, STVG == T) %>% 
  .$key %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    map_feature_plot(gene_ = gene_) + 
      patchwork::plot_annotation(
        title = gene_
      )
  }) -> p_list

pdf("~/Downloads/p_list_100.pdf", width = 15.14, height = 7.73)
p_list
dev.off()



list(
  CT8_SVG = CT8_SVG_1,
  CT20_SVG = CT20_SVG_1,
  TVG = TVG) %>% 
  ggvenn::ggvenn()

list(
  CT8_SVG = CT8_SVG_1,
  CT20_SVG = CT20_SVG_1
) %>% ggvenn::ggvenn() #Figure 3A after filtering

save(CT8_SVG_1, CT20_SVG_1, file = "~/Downloads/New_SVG.rda")
save(TVG, file = "~/Downloads/New_TVG.rda")
save(STVG_1, file = "~/Downloads/New_STVG.rda")

df_spatiotemporal_all_2 %>% dplyr::filter(gene == "Map1b") %>% dplyr::select(`CT-8_cauchy_padj_sparkx`, 
                                                                             `CT-20_cauchy_padj_sparkx`, 
                                                                             p_adj_TVG, 
                                                                             spatiotemp.padj_spatiotemporal, 
                                                                             CT8_cauchy_padj_core_vs_shell, 
                                                                             CT20_cauchy_padj_core_vs_shell)
patchwork::wrap_plots(map_feature_plot(gene_ = "Map1b"), DEG_plot(gene_ = "Map1b"), ncol = 2, widths = c(3,1))

df_spatiotemporal_all_2 %>% dplyr::filter(gene == "Taf12") %>% dplyr::select(`CT-8_cauchy_padj_sparkx`, 
                                                                             `CT-20_cauchy_padj_sparkx`, 
                                                                             p_adj_TVG, 
                                                                             spatiotemp.padj_spatiotemporal, 
                                                                             CT8_cauchy_padj_core_vs_shell, 
                                                                             CT20_cauchy_padj_core_vs_shell)
patchwork::wrap_plots(map_feature_plot(gene_ = "Taf12"), DEG_plot(gene_ = "Taf12"), ncol = 2, widths = c(3,1))

Prpf8

dim(norm_counts_matrix)
#SVG but not TVG
norm_counts_matrix["Csnk1e", ] %>% {.!=0} %>% sum() %>% {./ncol(norm_counts_matrix)} %>% {.*100}
norm_counts_matrix["Calb1", ] %>% {.!=0} %>% sum() %>% {./ncol(norm_counts_matrix)} %>% {.*100}
norm_counts_matrix["Slc39a6", ] %>% {.!=0} %>% sum() %>% {./ncol(norm_counts_matrix)} %>% {.*100}
#SVG and TVG
norm_counts_matrix["Lhx1", ] %>% {.!=0} %>% sum() %>% {./ncol(norm_counts_matrix)} %>% {.*100}
norm_counts_matrix["Rgs16", ] %>% {.!=0} %>% sum() %>% {./ncol(norm_counts_matrix)} %>% {.*100}
norm_counts_matrix["Grin2c", ] %>% {.!=0} %>% sum() %>% {./ncol(norm_counts_matrix)} %>% {.*100}


df_spatiotemporal_all_2 %>% 
  column_to_rownames("gene") %>% 
  .[c("Csnk1e", "Calb1", "Slc39a6", "Lhx1", "Rgs16", "Grin2c"), ] %>% 
  dplyr::select(`CT-8_cauchy_padj_sparkx`, 
               `CT-20_cauchy_padj_sparkx`, 
               p_adj_TVG, 
               spatiotemp.padj_spatiotemporal, 
               CT8_cauchy_padj_core_vs_shell, 
               CT20_cauchy_padj_core_vs_shell)

#gene alwasys enriched in shell
gene_se = c("Kcng2", "Il4ra", "Grin2c", "Il33", "Hopx", "Dlx1", "Zic4", "F3", "Gck", "Khdrbs2")

cbind(
df_spatiotemporal_all_2 %>% 
  column_to_rownames("gene") %>% 
  .[gene_se, ] %>% 
  dplyr::select(`CT-8_cauchy_padj_sparkx`, 
                `CT-20_cauchy_padj_sparkx`, 
                p_adj_TVG, 
                logFC_TVG,
                spatiotemp.padj_spatiotemporal, 
                CT8_cauchy_padj_core_vs_shell, 
                CT20_cauchy_padj_core_vs_shell),
gene_percent_df %>% column_to_rownames("gene") %>% .[gene_se, ] %>% as.data.frame() %>% "colnames<-"(., c("percent"))
) %>% dplyr::filter(spatiotemp.padj_spatiotemporal < 0.05) %>% 
  dplyr::mutate(TVG_group = case_when(
    logFC_TVG > 0 ~ "CT20_enriched",
    TRUE ~ "CT8_enriched"
  ))
patchwork::wrap_plots(map_feature_plot(gene_ = "Il33"), DEG_plot(gene_ = "Il33"), ncol = 2, widths = c(3,1))
patchwork::wrap_plots(map_feature_plot(gene_ = "Gck"), DEG_plot(gene_ = "Gck"), ncol = 2, widths = c(3,1))

map_feature_plot(gene_ = "Slc9a3r1")

load("~/Downloads/df_day_night.rda")
t_test_day_night = function(gene_){
  mat_norm_mean_long %>% dplyr::filter(gene == gene_) -> df_
  day = df_ %>% dplyr::filter(day_night == "day") %>% .$norm_expr
  night = df_ %>% dplyr::filter(day_night == "night") %>% .$norm_expr
  t.test(day, night)
}

map_feature_plot(gene_ = "Gad1") -> p1 #Figure 3F
DEG_plot(gene_ = "Gad1") -> p2 #Figure 3F
scRNA_lineplot(gene_ = "Gad1") -> p3
scRNA_lineplot_1(data = p3$data, CT_start = 26, CT_end = 46, gene_ = "Gad1")
scRNA_lineplot_2_new(data = p3$data, CT_start = 26, CT_end = 46, gene_ = "Gad1") -> p3

map_feature_plot(gene_ = "Per1") -> p1
DEG_plot(gene_ = "Per1") -> p2
scRNA_lineplot(gene_ = "Per1") -> p3
scRNA_lineplot_1(data = p3$data, CT_start = 38, CT_end = 58, gene_ = "Per1")
scRNA_lineplot_2(data = p3$data, CT_start = 38, CT_end = 58, gene_ = "Per1") -> p3

patchwork::wrap_plots(p1,p2,p3, ncol = 3)
df_day_night_CT54_CT42 %>% dplyr::filter(gene == "Per1")

map_feature_plot(gene_ = "Per2") -> p1 #Figure 4C
DEG_plot(gene_ = "Per2") -> p2 #Figure S10A
scRNA_lineplot(gene_ = "Per2") -> p3
scRNA_lineplot_1(data = p3$data, CT_start = 38, CT_end = 58, gene_ = "Per2")
scRNA_lineplot_2(data = p3$data, CT_start = 38, CT_end = 58, gene_ = "Per2") -> p3

patchwork::wrap_plots(p1,p2,p3, ncol = 3)
df_day_night_CT54_CT42 %>% dplyr::filter(gene == "Per2")

map_feature_plot(gene_ = "Lhx1") -> p1 #Figure 4E
DEG_plot(gene_ = "Lhx1") -> p2 #Figure S10B

map_feature_plot(gene_ = "Rgs16") -> p1 #Figure 4F
DEG_plot(gene_ = "Rgs16") -> p2 #Figure S10C

map_feature_plot(gene_ = "Il33") -> p1
DEG_plot(gene_ = "Il33") -> p2
scRNA_lineplot(gene_ = "Il33") -> p3
scRNA_lineplot_1(data = p3$data, CT_start = 14, CT_end = 34, gene_ = "Il33")
scRNA_lineplot_2(data = p3$data, CT_start = 14, CT_end = 34, gene_ = "Il33") -> p3

patchwork::wrap_plots(p1,p2,p3, ncol = 3)

scRNA_lineplot(gene_ = "Gck") -> p_Gck
scRNA_lineplot_1(data = p_Gck$data, CT_start = 14, CT_end = 34, gene_ = "Gck")
scRNA_lineplot_2(data = p_Gck$data, CT_start = 14, CT_end = 34, gene_ = "Gck")

scRNA_lineplot(gene_ = "Csnk1e") -> p_Csnk1e
scRNA_lineplot_1(data = p_Csnk1e$data, CT_start = 14, CT_end = 34, gene_ = "Csnk1e")
scRNA_lineplot_2(data = p_Csnk1e$data, CT_start = 14, CT_end = 34, gene_ = "Csnk1e") -> p_Csnk1e
ggsave(plot = p_Csnk1e, "~/Downloads/p_Csnk1e_4.pdf", width = 8.21, height = 4.27)


scRNA_lineplot(gene_ = "Calb1") -> p_Calb1 #Figure 3G
scRNA_lineplot_1(data = p_Calb1$data, CT_start = 14, CT_end = 34, gene_ = "Calb1")
scRNA_lineplot_2(data = p_Calb1$data, CT_start = 14, CT_end = 34, gene_ = "Calb1") -> p_Calb1
ggsave(plot = p_Calb1, "~/Downloads/p_Calb1_4.pdf", width = 8.21, height = 4.27)


scRNA_lineplot(gene_ = "Slc39a6") -> p_Slc39a6 #Figure 3H
scRNA_lineplot_1(data = p_Slc39a6$data, CT_start = 38, CT_end = 58, gene_ = "Slc39a6")
scRNA_lineplot_2(data = p_Slc39a6$data, CT_start = 38, CT_end = 58, gene_ = "Slc39a6") -> p_Slc39a6
ggsave(plot = p_Slc39a6, "~/Downloads/p_Slc39a6_4.pdf", width = 8.21, height = 4.27)

scRNA_lineplot(gene_ = "Lhx1") -> p_Lhx1
scRNA_lineplot_1(data = p_Lhx1$data, CT_start = 14, CT_end = 34, gene_ = "Lhx1")
scRNA_lineplot_2(data = p_Lhx1$data, CT_start = 14, CT_end = 34, gene_ = "Lhx1") -> p_Lhx1
ggsave(plot = p_Lhx1, "~/Downloads/p_Lhx1_4.pdf", width = 8.21, height = 4.27)

scRNA_lineplot(gene_ = "Rgs16") -> p_Rgs16
scRNA_lineplot_1(data = p_Rgs16$data, CT_start = 38, CT_end = 58, gene_ = "Rgs16")
scRNA_lineplot_2(data = p_Rgs16$data, CT_start = 38, CT_end = 58, gene_ = "Rgs16") -> p_Rgs16
ggsave(plot = p_Rgs16, "~/Downloads/p_Rgs16_4.pdf", width = 8.21, height = 4.27)

scRNA_lineplot(gene_ = "Grin2c") -> p_Grin2c
scRNA_lineplot_1(data = p_Grin2c$data, CT_start = 14, CT_end = 34, gene_ = "Grin2c")
scRNA_lineplot_2(data = p_Grin2c$data, CT_start = 14, CT_end = 34, gene_ = "Grin2c") -> p_Grin2c
ggsave(plot = p_Grin2c, "~/Downloads/p_Grin2c_4.pdf", width = 8.21, height = 4.27)

###spatialtemporal analysis ----
library(tidyverse)
load("~/Downloads/df_spatiotemporal_all_2.rda")
load("~/Downloads/New_SVG.rda")
load("~/Downloads/New_TVG.rda")
load("~/Downloads/New_STVG.rda")

list(
  CT8_SVG = CT8_SVG_1,
  CT20_SVG = CT20_SVG_1,
  TVG = TVG,
  STVG = STVG_1
) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::mutate(across(!matches("key"), function(x){x %>% as.integer()})) %>%
  as.data.frame() %>%
  UpSetR::upset(sets = rev(c("TVG", "CT8_SVG", "CT20_SVG", "STVG")), keep.order = T, order.by = "freq") -> p_upset

pdf("~/Downloads/upset_plot.pdf", width = 7.65, height = 4.86)
p_upset
dev.off()


list(
  CT8_SVG = CT8_SVG_1,
  CT20_SVG = CT20_SVG_1,
  TVG = TVG,
  STVG = STVG_1
) %>% ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    df_ = df_ %>% dplyr::filter(STVG == T, TVG == T)
    df_ %>% dplyr::filter(((CT8_SVG==T)|(CT20_SVG==T))) -> df_1
    dim(df_1)
    df_1$key
  } -> STVG_TVG_SVG_intersection

df_spatiotemporal_all_2 %>% dplyr::filter(gene %in% STVG_TVG_SVG_intersection) %>% 
  dplyr::select(gene, CT8_cauchy_padj_core_vs_shell, CT20_cauchy_padj_core_vs_shell, `CT-8_avg_log2FC_core_vs_shell`, `CT-20_avg_log2FC_core_vs_shell`) %>% 
  dplyr::mutate(Sig = case_when(
    (CT8_cauchy_padj_core_vs_shell<0.05)&(CT20_cauchy_padj_core_vs_shell<0.05) ~ "Sig",
    TRUE ~ "Not Sig"
  )) -> df_plot

library(ggrepel)
df_plot %>% 
  {
    df_plot = .
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC_core_vs_shell` < 0, `CT-20_avg_log2FC_core_vs_shell` > 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC_core_vs_shell`)), desc(abs(`CT-20_avg_log2FC_core_vs_shell`))) %>% 
      head(10) -> df_1
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC_core_vs_shell` < 0, `CT-20_avg_log2FC_core_vs_shell` < 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC_core_vs_shell`)), desc(abs(`CT-20_avg_log2FC_core_vs_shell`))) %>% 
      head(10) -> df_2
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC_core_vs_shell` > 0, `CT-20_avg_log2FC_core_vs_shell` > 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC_core_vs_shell`)), desc(abs(`CT-20_avg_log2FC_core_vs_shell`))) %>% 
      head(10) -> df_3
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC_core_vs_shell` > 0, `CT-20_avg_log2FC_core_vs_shell` < 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC_core_vs_shell`)), desc(abs(`CT-20_avg_log2FC_core_vs_shell`))) %>% 
      head(10) -> df_4
    
    rbind(df_1, df_2, df_3, df_4) -> top10
    top10
    
    df_plot %>% 
      ggplot(aes(x = `CT-8_avg_log2FC_core_vs_shell`, y = `CT-20_avg_log2FC_core_vs_shell`, color = Sig)) + 
      geom_point(alpha = 0.5) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + 
      theme_classic() + 
      geom_text_repel(
        data = top10,
        aes(
          x = `CT-8_avg_log2FC_core_vs_shell`,
          y = `CT-20_avg_log2FC_core_vs_shell`,
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
      ggtitle("Spatiotemporal genes (SVG + TVG + STVG) enriched in core/shell")
  } -> p #Figure 4B

ggsave(plot = p, filename = "~/Downloads/Spatiotemporal_TVG_and_SVG_geom_point.pdf", width = 11.41, height = 6.36)

p$data %>% dplyr::filter(Sig == "Sig") %>% dim()
p$data %>% dplyr::filter(Sig == "Sig") %>% 
  {
    df_ = .
    df_$gene -> genes_
    df_spatiotemporal_all_2 %>% dplyr::filter(gene %in% genes_) %>% dplyr::select(gene, logFC_TVG) -> df_1
    full_join(df_, df_1, by = "gene") %>% 
      dplyr::mutate(CT8_SVG = case_when(
        `CT-8_avg_log2FC_core_vs_shell` > 0 ~ "Core", 
        TRUE ~ "Shell"
      )) %>% 
      dplyr::mutate(CT20_SVG = case_when(
        `CT-20_avg_log2FC_core_vs_shell` > 0 ~ "Core", 
        TRUE ~ "Shell"
      )) %>% 
      dplyr::mutate(TVG = case_when(
        logFC_TVG > 0 ~ "CT20",
        TRUE ~ "CT8"
      )) -> df_2
    c("CT8", "CT20") %>% 
      purrr::map(function(CT_){
        df_2 %>% dplyr::filter(TVG == CT_) %>% 
          dplyr::arrange(CT8_SVG, CT20_SVG) %>% 
          dplyr::select(gene, CT8_SVG, CT20_SVG, TVG)
      })

  }

list(
  CT8_SVG = CT8_SVG_1,
  CT20_SVG = CT20_SVG_1,
  TVG = TVG) %>% 
  ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    df_ %>% dplyr::filter(TVG == F) -> df_
    
    df_ %>% dplyr::filter(CT8_SVG == T, CT20_SVG == F) %>% .$key -> CT8_SVG
    df_ %>% dplyr::filter(CT8_SVG == F, CT20_SVG == T) %>% .$key -> CT20_SVG
    df_ %>% dplyr::filter(CT8_SVG == T, CT20_SVG == T) %>% .$key -> Both
    list(
      CT8_only = CT8_SVG,
      CT20_only = CT20_SVG,
      Both = Both
    ) %>% 
#      .[2] %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(genes_,group_){
        df_spatiotemporal_all_2 %>% dplyr::filter(gene %in% genes_) %>% 
          dplyr::select(gene, `CT-8_cauchy_padj_sparkx`, `CT-20_cauchy_padj_sparkx`) -> df_
        if (group_ == "CT8_only"){
          df_ %>% dplyr::arrange(`CT-8_cauchy_padj_sparkx`, desc(`CT-20_cauchy_padj_sparkx`)) %>% dplyr::mutate(group = group_) -> df_1
        }else if (group_ == "CT20_only"){
          df_ %>% dplyr::arrange(`CT-20_cauchy_padj_sparkx`, desc(`CT-8_cauchy_padj_sparkx`)) %>% dplyr::mutate(group = group_) -> df_1
        }else{
          df_ %>% dplyr::arrange(`CT-8_cauchy_padj_sparkx`, `CT-20_cauchy_padj_sparkx`) %>% dplyr::mutate(group = group_) -> df_1
        }
        df_spatiotemporal_all_2 %>% 
          dplyr::select(gene, CT8_cauchy_padj_core_vs_shell, CT20_cauchy_padj_core_vs_shell, 
                        `CT-8_avg_log2FC_core_vs_shell`, `CT-20_avg_log2FC_core_vs_shell`) %>% 
          dplyr::mutate(CT8_SVG_group = case_when(
            (CT8_cauchy_padj_core_vs_shell < 0.05)&(`CT-8_avg_log2FC_core_vs_shell` > 0) ~ "Core",
            (CT8_cauchy_padj_core_vs_shell < 0.05)&(`CT-8_avg_log2FC_core_vs_shell` < 0) ~ "Shell",
            TRUE ~ "None"
          )) %>% 
          dplyr::mutate(CT20_SVG_group = case_when(
            (CT20_cauchy_padj_core_vs_shell < 0.05)&(`CT-20_avg_log2FC_core_vs_shell` > 0) ~ "Core",
            (CT20_cauchy_padj_core_vs_shell < 0.05)&(`CT-20_avg_log2FC_core_vs_shell` < 0) ~ "Shell",
            TRUE ~ "None"
          )) -> df_2
        left_join(df_1, df_2, by = "gene") -> df_3
#        df_3 %>% dplyr::filter(!((CT8_SVG_group=="None")&(CT20_SVG_group=="None"))) -> df_4
        df_3 -> df_4
        df_4 %>% dplyr::select(gene, CT8_SVG_group,  CT20_SVG_group) -> df_5
        df_spatiotemporal_all_2 %>% dplyr::select(gene, `CT-8_avg_log2FC_core_vs_shell`, `CT-20_avg_log2FC_core_vs_shell`) -> df_6
        left_join(df_5, df_6, by = "gene")
      })
  } -> plot_dat_list

plot_dat_list$CT8_only %>% 
  dplyr::mutate(Sig = case_when(
    (CT8_SVG_group=="None")&(CT20_SVG_group=="None") ~ "Not Sig",
      TRUE ~ "Sig"
  )) %>% 
  {
    df_plot = .
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC_core_vs_shell` < 0, `CT-20_avg_log2FC_core_vs_shell` > 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC_core_vs_shell`)), desc(abs(`CT-20_avg_log2FC_core_vs_shell`))) %>% 
      head(10) -> df_1
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC_core_vs_shell` < 0, `CT-20_avg_log2FC_core_vs_shell` < 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC_core_vs_shell`)), desc(abs(`CT-20_avg_log2FC_core_vs_shell`))) %>% 
      head(10) -> df_2
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC_core_vs_shell` > 0, `CT-20_avg_log2FC_core_vs_shell` > 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC_core_vs_shell`)), desc(abs(`CT-20_avg_log2FC_core_vs_shell`))) %>% 
      head(10) -> df_3
    
    df_plot %>% 
      dplyr::filter(Sig == "Sig") %>% 
      dplyr::filter(`CT-8_avg_log2FC_core_vs_shell` > 0, `CT-20_avg_log2FC_core_vs_shell` < 0) %>% 
      dplyr::arrange(desc(abs(`CT-8_avg_log2FC_core_vs_shell`)), desc(abs(`CT-20_avg_log2FC_core_vs_shell`))) %>% 
      head(10) -> df_4
    
    rbind(df_1, df_2, df_3, df_4) -> top10
    top10
    
    df_plot %>% 
      ggplot(aes(x = `CT-8_avg_log2FC_core_vs_shell`, y = `CT-20_avg_log2FC_core_vs_shell`, color = Sig)) + 
      geom_point(alpha = 0.5) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + 
      theme_classic() + 
      geom_text_repel(
        data = top10,
        aes(
          x = `CT-8_avg_log2FC_core_vs_shell`,
          y = `CT-20_avg_log2FC_core_vs_shell`,
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
      ggtitle("Spatiotemporal genes (SVG + TVG + STVG) enriched in core/shell")
  } -> p

p$data %>% dplyr::filter(gene == "Csnk1e")

plot_dat_list$CT8_only %>% 
  dplyr::filter(CT8_SVG_group != "None", CT20_SVG_group == "None") %>% 
  dplyr::arrange(desc(abs(`CT-8_avg_log2FC_core_vs_shell`))) %>% 
  head(29) -> df_plot

rbind(df_plot, 
plot_dat_list$Both %>% 
  dplyr::filter(CT20_SVG_group != "None", CT8_SVG_group != "None") %>% 
  dplyr::filter(gene %in% c("Slc39a6"))
) %>% dplyr::arrange(desc(abs(`CT-8_avg_log2FC_core_vs_shell`))) -> df_plot

df_plot %>% dplyr::arrange(`CT-8_avg_log2FC_core_vs_shell`) -> df_plot

df_plot %>% {
  df_plot = .
  df_plot$gene = factor(df_plot$gene, rev(df_plot$gene))
  df_plot %>% 
#    ggplot(aes(x = abs(`CT-8_avg_log2FC_core_vs_shell`), y = gene, fill = CT8_SVG_group)) + 
    ggplot(aes(x = `CT-8_avg_log2FC_core_vs_shell`, y = gene, fill = CT8_SVG_group)) + 
    geom_col() + 
    theme_classic()
} -> p_CT8_only

plot_dat_list$CT20_only %>% 
  dplyr::filter(CT20_SVG_group != "None", CT8_SVG_group == "None") %>% 
  dplyr::arrange(desc(abs(`CT-20_avg_log2FC_core_vs_shell`))) %>% 
  head(29) -> df_plot

rbind(df_plot, 
      plot_dat_list$Both %>% 
        dplyr::filter(CT20_SVG_group != "None", CT8_SVG_group != "None") %>% 
        dplyr::filter(gene %in% c("Slc39a6"))
) %>% dplyr::arrange(desc(abs(`CT-20_avg_log2FC_core_vs_shell`))) -> df_plot

df_plot %>% dplyr::arrange(`CT-20_avg_log2FC_core_vs_shell`) -> df_plot

df_plot %>% {
  df_plot = .
  df_plot$gene = factor(df_plot$gene, rev(df_plot$gene))
  df_plot %>% 
#    ggplot(aes(x = abs(`CT-20_avg_log2FC_core_vs_shell`), y = gene, fill = CT20_SVG_group)) + 
    ggplot(aes(x = `CT-20_avg_log2FC_core_vs_shell`, y = gene, fill = CT20_SVG_group)) + 
    geom_col() + 
    theme_classic()
} -> p_CT20_only

list(CT8_SVG = CT8_SVG_1, 
     TVG = TVG) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::filter(CT8_SVG == T) %>%
  .$TVG %>% {ifelse(. == T, "TVG", "Not_TVG")} %>% 
  table() %>% 
  as.data.frame() %>% 
  "colnames<-"(., c("group", "freq")) %>% 
  dplyr::mutate(percent = round((freq/sum(freq))*100, 2)) %>% 
  dplyr::mutate(Label = sprintf("%s\n(%s)", freq, percent)) %>% 
  ggplot(aes(x = "", y = freq, fill = group)) + 
  geom_col(width = 1, color = "white") + 
  coord_polar(theta = "y") + 
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) + 
  theme_void() -> p1 #Figure 3D

list(CT20_SVG = CT20_SVG_1, 
     TVG = TVG) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::filter(CT20_SVG == T) %>%
  .$TVG %>% {ifelse(. == T, "TVG", "Not_TVG")} %>% 
  table() %>% 
  as.data.frame() %>% 
  "colnames<-"(., c("group", "freq")) %>% 
  dplyr::mutate(percent = round((freq/sum(freq))*100, 2)) %>% 
  dplyr::mutate(Label = sprintf("%s\n(%s)", freq, percent)) %>% 
  ggplot(aes(x = "", y = freq, fill = group)) + 
  geom_col(width = 1, color = "white") + 
  coord_polar(theta = "y") + 
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) + 
  theme_void() -> p2 #Figure 3D

patchwork::wrap_plots(p1, p2, ncol = 2, guides = "collect")


list(SVG = c(CT8_SVG_1, CT20_SVG_1) %>% unique(), 
     TVG = TVG) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::filter(SVG == T) %>%
  .$TVG %>% {ifelse(. == T, "TVG", "Not_TVG")} %>% 
  table() %>% 
  as.data.frame() %>% 
  "colnames<-"(., c("group", "freq")) %>% 
  dplyr::mutate(percet = round((freq/sum(freq))*100, 2))


plot_dat_list$Both %>% 
  dplyr::filter(CT20_SVG_group != "None", CT8_SVG_group != "None") %>% 
  rowwise() %>% 
  dplyr::mutate(log2fc = mean(c(abs(`CT-8_avg_log2FC_core_vs_shell`), abs(`CT-20_avg_log2FC_core_vs_shell`)))) %>% 
  dplyr::arrange(desc(log2fc)) %>% 
  ungroup() %>% as.data.frame() %>% 
  {
    df_plot = .
    df_plot_1 = df_plot %>% head(27)
    df_plot_2 = df_plot %>% dplyr::filter(CT8_SVG_group != CT20_SVG_group)
    rbind(
      df_plot_1,
      df_plot_2
    ) -> df_plot_3
    df_plot_3 %>% dplyr::mutate(group = sprintf("%s -> %s", CT8_SVG_group, CT20_SVG_group)) %>% 
      dplyr::mutate(gene = factor(gene, levels = rev(gene))) %>% 
      ggplot(aes(x = log2fc, y = gene, fill = group)) + 
      geom_col() + 
      theme_classic()
  } -> p_both
  

patchwork::wrap_plots(p_CT8_only, p_CT20_only, p_both, ncol = 3) -> p
ggsave(plot = p, filename = "~/Downloads/SVG_not_TVG.pdf", width = 16.53, height = 5.60) #Figure 3E

map_feature_plot(gene_ = "Adra1b")

map_feature_plot(gene_ = "Gad1") #plot 4 and 4
DEG_plot(gene_ = "Gad1")
scRNA_lineplot(gene_ = "Gad1") -> p
scRNA_lineplot_3(data = p$data, gene_ = "Gad1", CT_start = 26, CT_end = 46)

map_feature_plot(gene_ = "Galnt18") #plot 4 and 4

DEG_plot(gene_ = "Gad1")
scRNA_lineplot(gene_ = "Gad1") -> p
scRNA_lineplot_3(data = p$data, gene_ = "Gad1", CT_start = 26, CT_end = 46)


list(
  CT8_SVG = CT8_SVG_1,
  CT20_SVG = CT20_SVG_1,
  TVG = TVG,
  STVG = STVG_1) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::filter(CT8_SVG == F, CT20_SVG == F, TVG == F, STVG == T) %>% 
  left_join(x = ., y = df_spatiotemporal_all_2, by = c("key"="gene")) %>% 
  dplyr::select(key, spatiotemp.padj_spatiotemporal) %>% 
  left_join(x = ., y = gene_percent_df, by = c("key"="gene")) %>% 
  dplyr::arrange(spatiotemp.padj_spatiotemporal, percent) -> df_

#STVG only
df_$key %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    DEG_plot(gene_ = gene_)
  }) -> p_list

pdf("~/Downloads/p_list_new.pdf", width = 13.79, height = 5.60)
p_list
dev.off()


df_$key %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    map_feature_plot(gene_ = gene_)
  }) -> p_list_1

p_list$Numb
p_list_1$Numb

p_list$Kdm6a
p_list_1$Kdm6a

p_list$


df_spatiotemporal_all_2 %>% dplyr::filter(gene == "Map1b")
map_feature_plot(gene_ = "Map1b") -> p
patchwork::wrap_plots(p[[3]], p[[4]], ncol = 2, guides = "collect") -> p_1
DEG_plot(gene_ = "Map1b") -> p_2
patchwork::wrap_plots(p_1, p_2, ncol = 2, widths = c(2,1)) -> p_Map1b

#Myh9
map_feature_plot(gene_ = "Myh9") -> p
patchwork::wrap_plots(p[[2]], p[[4]], ncol = 2, guides = "collect") -> p_1
DEG_plot(gene_ = "Myh9") -> p_2
patchwork::wrap_plots(p_1, p_2, ncol = 2, widths = c(2,1)) -> p_Myh9

#Cacna1g
map_feature_plot(gene_ = "Cacna1g") -> p #Figure S11A
patchwork::wrap_plots(p[[1]], p[[6]], ncol = 2, guides = "collect") -> p_1
DEG_plot(gene_ = "Cacna1g") -> p_2 #Figure S11A
patchwork::wrap_plots(p_1, p_2, ncol = 2, widths = c(2,1)) -> p_Cacna1g

#Smarcc1
#map_feature_plot(gene_ = "Smarcc1") -> p
#patchwork::wrap_plots(p[[1]], p[[6]], ncol = 2, guides = "collect") -> p_1
#DEG_plot(gene_ = "Smarcc1") -> p_2
#patchwork::wrap_plots(p_1, p_2, ncol = 2, widths = c(2,1)) -> p_Cacna1g

#Grm8
map_feature_plot(gene_ = "Grm8") -> p
patchwork::wrap_plots(p[[3]], p[[6]], ncol = 2, guides = "collect") -> p_1
DEG_plot(gene_ = "Grm8") -> p_2
patchwork::wrap_plots(p_1, p_2, ncol = 2, widths = c(2,1)) -> p_Grm8

#Rbms1
map_feature_plot(gene_ = "Rbms1") -> p
patchwork::wrap_plots(p[[2]], p[[5]], ncol = 2, guides = "collect") -> p_1
DEG_plot(gene_ = "Rbms1") -> p_2
patchwork::wrap_plots(p_1, p_2, ncol = 2, widths = c(2,1)) -> p_Rbms1

#Gigyf2
map_feature_plot(gene_ = "Gigyf2") -> p
patchwork::wrap_plots(p[[1]], p[[5]], ncol = 2, guides = "collect") -> p_1
DEG_plot(gene_ = "Gigyf2") -> p_2
patchwork::wrap_plots(p_1, p_2, ncol = 2, widths = c(2,1)) -> p_Gigyf2

pdf("~/Downloads/STVG_only.pdf", width = 15.6, height = 4.80)
p_Map1b
p_Myh9
p_Cacna1g
p_Grm8
p_Rbms1
p_Gigyf2
dev.off()

