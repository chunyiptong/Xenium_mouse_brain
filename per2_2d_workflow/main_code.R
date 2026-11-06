library(tidyverse)

# Generate ImageJ points files for each gene and time point ----
source("~/Downloads/20261010/per2_2d_workflow/write_imagej_points.R")
list.files("~/Downloads/20261010/vip-cck-per2-cor", pattern = "points\\.tsv", recursive = T, full.names = T) -> files_
c("CT2", "CT6", "CT10", "CT14", "CT18", "CT22") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(CT_){
    files_ %>% .[grepl(sprintf("/%s/" ,CT_), .)] -> files_
    c("Per2", "Cck", "Vip") %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(gene_){
        files_ %>% .[grepl(sprintf("%s" ,gene_), .)] -> file_ 
        read.table(file_, header = T, sep = "\t") -> points
        gsub("\\.tsv", ".roi", file_) -> out_
        write_imagej_points(points, out_)
      })
  })

# Determine if smFISh spot is inside the SCN ----
library(sp)
library(RImageJROI)

source("/Users/chunyiptong/Downloads/20261010/per2_2d_workflow/classify_spots_scn.R")
c("CT2", "CT6", "CT10", "CT14", "CT18", "CT22") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(CT_){
    files_ %>% .[grepl(sprintf("/%s/" ,CT_), .)] -> files_
    c("Per2", "Cck", "Vip") %>% 
      "names<-"(.,.) %>% 
      purrr::map(function(gene_){
        files_ %>% .[grepl(sprintf("%s" ,gene_), .)] -> file_ 
        gsub("(/.+)/.+/.+", "\\1", file_) -> dir_
        list.files(dir_, pattern = "\\ROIs.zip$", recursive = F, full.names = T) -> roi_file_
        gsub("(/.+/.+)/.+", "\\1", file_) -> dir_
        results = classify_spots_scn(spots_file = file_, roi_file = roi_file_, output_dir = dir_, gene = gene_)
      })
  }) -> results_list

library(png)
source("/Users/chunyiptong/Downloads/20261010/per2_2d_workflow/assign_spots_to_mask_v3.R")
results_list %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(list_, CT_){
    list_ %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(list__, gene_){
        list__$spots -> spots
        list.files("~/Downloads/20261010/vip-cck-per2-cor", pattern = "cp_masks.png", recursive = T, full.names = T) %>% 
          .[grepl(sprintf("/%s/", CT_), .)] -> mask_file
        list.files("~/Downloads/20261010/vip-cck-per2-cor", pattern = "ROIs\\.zip", recursive = T, full.names = T) %>% 
          .[grepl(sprintf("/%s/", CT_), .)] -> roi_file
        assigned <- assign_spots_to_mask_v3(
          spots = spots,
          mask_file = mask_file,
          roi_file = roi_file,
          image_id = CT_,
          mask_kind = "nuclei"
        )
      })
  }) -> results_list_1

results_list_1 %>% 
  purrr::map(function(list_){
    list_ %>% 
      purrr::map2(.x=.,.y=names(.),.f=function(list__, gene_){
        list__$spots %>% dplyr::filter(in_SCN == T, mask_assignment == "inside_nucleus") -> df_
        df_$cell_id %>% table() %>% as.data.frame() %>% "colnames<-"(., c("cell", gene_))
      }) %>% purrr::reduce(., full_join, by = "cell") -> df_1
    df_1[is.na(df_1)] <- 0
    df_1 %>% column_to_rownames("cell") %>% t() %>% as.data.frame() -> df_1
    df_1
  }) -> smFISH_SCN_count_mat

smFISH_SCN_count_mat %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, CT_){
    df_ %>% t() %>% as.data.frame() -> df_
    df_ %>% 
      ggplot(aes(x = Vip, y = Per2)) + 
      geom_count(alpha = 0.6) +
      scale_size_area(max_size = 10, name = "Number of nucleus", limits = c(0,400), breaks = c(50, 100, 200, 300, 400)) +
      theme_classic() + 
      ggtitle(CT_) + 
      ylab("Per2 (Count per nucleus)") + 
      xlab("Vip (Count per nucleus)")
  }) %>% patchwork::wrap_plots(., ncol = 3, guides = "collect")

smFISH_SCN_count_mat %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, CT_){
    df_ %>% t() %>% as.data.frame() -> df_
    df_ %>% 
      ggplot(aes(x = Cck, y = Per2)) + 
      geom_count(alpha = 0.6) +
      scale_size_area(max_size = 10, name = "Number of cells", limits = c(0,400), breaks = c(50, 100, 200, 300, 400)) +
      theme_classic() + 
      ggtitle(CT_) + 
      ylab("Per2 (Count per cell)") + 
      xlab("Cck (Count per cell)")
  }) %>% patchwork::wrap_plots(., ncol = 3, guides = "collect")


smFISH_SCN_count_mat %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_,CT_){
    CT_ = gsub("CT", "", CT_) %>% as.numeric()
    df_ %>% t() %>% as.data.frame() %>% dplyr::mutate(CT = CT_)
  }) %>% do.call(rbind, .) -> dat_

factor(dat_$CT) -> dat_$CT

fit <- MASS::glm.nb(
  Per2 ~ Vip * CT,
  data = dat_
)

fit_pois <- glm(
  Per2 ~ Vip * CT,
  family = poisson(link = "log"),
  data = dat_
)

AIC(fit_pois, fit)

summary(fit)
pred <- dat_ %>%
  group_by(CT) %>%
  reframe(
    Vip = seq(min(Vip), max(Vip), length.out = 100)
  )
pred$Per2_pred <- predict(
  fit,
  newdata = pred,
  type = "response"
)

ggplot(pred, aes(Vip, Per2_pred, colour = CT, group = CT)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Vip count per nucleus",
    y = "Predicted mean Per2 count per nucleus",
    colour = "CT"
  ) +
  theme_classic(base_size = 14)

results_list_1$CT22$Per2$spots

c("Per2", "Vip") %>% 
  "names<-"(.,.) %>% 
#    .[1] %>% 
  purrr::map(function(gene_){
    results_list %>% 
#      .[2] %>% 
      purrr::map2(.x=., .y=names(.),.f=function(list_, CT_){
        list_[[gene_]]$spots -> df_
        df_ %>% 
          ggplot(aes(x = log(intensity), y = score)) + 
          geom_point() + 
          theme_classic() + 
          ggtitle(sprintf("%s (%s)", gene_, CT_))
      }) %>% patchwork::wrap_plots(., ncol = 3)
  }) -> p_list

c("Per2", "Vip") %>% 
  "names<-"(.,.) %>% 
  #    .[1] %>% 
  purrr::map(function(gene_){
    results_list %>% 
      #      .[2] %>% 
      purrr::map2(.x=., .y=names(.),.f=function(list_, CT_){
        list_[[gene_]]$spots -> df_
        df_ %>% 
          ggplot(aes(x = log(intensity))) + 
          geom_density() + 
          theme_classic() + 
          ggtitle(sprintf("%s (%s)", gene_, CT_))
      }) %>% patchwork::wrap_plots(., ncol = 3)
  }) -> p_list_1

p_list_1$Vip
    

c("Per2", "Vip") %>% 
#  .[1] %>% 
  purrr::map(function(gene_){
    results_list %>% 
      purrr::map2(.x=., .y=names(.),.f=function(list_, CT_){
        list_[[gene_]]$counts %>% dplyr::filter(region == "SCN_union") %>% .$n_spots -> n_spots
        data.frame(CT = CT_, n_spots = n_spots)
      }) %>% do.call(rbind, .) -> df_
    df_ %>% dplyr::mutate(CT = gsub("CT", "", CT) %>% as.numeric()) -> df_
    df_ %>% 
      ggplot(aes(x = CT, y = n_spots)) + 
      geom_line() + 
      geom_point() +
      scale_x_continuous(breaks = seq(2, 22, by = 4)) + 
      theme_classic() + 
      ggtitle(gene_)
  }) %>% patchwork::wrap_plots(., ncol = 2)

#New Analysis ----
source("~/Downloads/20261010/per2_2d_workflow/count_nuclei_in_scn.R")

c("CT2", "CT6", "CT10", "CT14", "CT18", "CT22") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(CT_){
    list.files("~/Downloads/20261010", pattern = "ROIs\\.zip", full.names = T, recursive = T) %>% .[grepl(sprintf("/%s/", CT_), .)] -> roi
    list.files("~/Downloads/20261010", pattern = "cp_masks\\.png", full.names = T, recursive = T) %>% .[grepl(sprintf("/%s/", CT_), .)] -> mask
    count_nuclei_in_scn(mask_file = mask, roi_file = roi)
  }) -> scn_n_nuclei

scn_n_nuclei$CT2$counts

list.files("~/Downloads/20261010/vip-cck-per2-cor", pattern = "cp_masks.png", recursive = T, full.names = T) -> mask_file
results_list_1$CT2$Per2$spots$mask_assignment %>% table()

results_list_1 %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(list_, CT_){
    list_$Per2$spots %>% dplyr::filter(mask_assignment == "inside_nucleus") %>% nrow() -> n_transcript_per2
    list_$Vip$spots %>% dplyr::filter(mask_assignment == "inside_nucleus") %>% nrow() -> n_transcript_vip
    n_transcript_per2_vip = n_transcript_per2 + n_transcript_vip
    scn_n_nuclei[[CT_]]$counts %>% dplyr::filter(region == "total") %>% .$n_nuclei -> n_nuclei
    data.frame(Per2 = n_transcript_per2/n_nuclei, Vip = n_transcript_vip/n_nuclei, Per2_Vip = n_transcript_per2_vip/n_nuclei)
  })

seq(2,22,4) %>% 
  "names<-"(., sprintf("CT%s", .)) %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(ct, ct_){
    dat_ %>% dplyr::filter(CT == ct) %>% dplyr::select(Per2, Vip) -> dat_
    dat_ %>% dplyr::filter(Per2 != 0, Vip != 0) %>% sum() -> n_transcript_Per2_Vip
    dat_$Per2 %>% sum() -> n_transcript_Per2
    dat_$Vip %>% sum() -> n_transcript_Vip
    
    scn_n_nuclei[[ct_]]$counts %>% dplyr::filter(region == "total") %>% .$n_nuclei -> n_nuclei
    data.frame(Per2 = n_transcript_Per2, Vip = n_transcript_Vip, Per2_Vip = n_transcript_Per2_Vip) -> df_
    df_/n_nuclei -> df_
    df_$CT = ct
    df_
  }) %>% do.call(rbind, .) -> df_plot

df_plot %>% pivot_longer(cols = -CT, names_to = "gene", values_to = "expr") %>% 
  dplyr::mutate(gene = factor(gene, levels = c("Vip", "Per2", "Per2_Vip"))) %>% 
  ggplot(aes(x = CT, y = expr, group = gene, fill = gene)) + 
  geom_col(position = "dodge") + 
  ylab("RNA concentration") + 
  scale_fill_manual(values = c(
    Vip      = "#F8766D",
    Per2     = "#00BA38",
    Per2_Vip = "#B8860B"
  )) + 
  scale_x_continuous(breaks = seq(2,22,4)) + 
  theme_classic()


seq(2,22,4) %>% 
  "names<-"(., sprintf("CT%s", .)) %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(ct, ct_){
    dat_ %>% dplyr::filter(CT == ct) %>% dplyr::select(Per2, Vip) -> dat_
    dat_ %>% dplyr::filter(Per2 != 0, Vip != 0) %>% nrow() -> n_nucleus_Per2_Vip
    dat_ %>% dplyr::filter(Per2 != 0) %>% nrow() -> n_nucleus_Per2
    dat_ %>% dplyr::filter(Vip != 0) %>% nrow() -> n_nucleus_Vip
#    
    scn_n_nuclei[[ct_]]$counts %>% dplyr::filter(region == "total") %>% .$n_nuclei -> n_total_nuclei
    data.frame(Per2 = n_nucleus_Per2, Vip = n_nucleus_Vip, Per2_Vip = n_nucleus_Per2_Vip) -> df_
    (df_/n_total_nuclei)*100 -> df_
    df_$CT = ct
    df_
  }) %>% do.call(rbind, .) -> df_plot

df_plot %>% pivot_longer(cols = -CT, names_to = "gene", values_to = "expr") %>% 
  dplyr::mutate(gene = factor(gene, levels = c("Vip", "Per2", "Per2_Vip"))) %>% 
  ggplot(aes(x = CT, y = expr, group = gene, fill = gene)) + 
  geom_col(position = "dodge") + 
  ylab("RNA concentration") + 
  scale_fill_manual(values = c(
    Vip      = "#F8766D",
    Per2     = "#00BA38",
    Per2_Vip = "#B8860B"
  )) + 
  scale_x_continuous(breaks = seq(2,22,4)) + 
  theme_classic()


results_list_1$CT2$Per2$cells
dat_ %>% rownames_to_column("cell") %>% dplyr::mutate(cell = gsub(".+\\.(.+)", "\\1", cell)) %>% dplyr::mutate(CT = as.character(CT)) -> dat_new
results_list_1 %>% 
  purrr::map(function(list_){
    list_[["Per2"]]$cell %>% #dplyr::filter(in_SCN == T) %>% 
      dplyr::select(cell_id, in_SCN)
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) %>% 
  {
    df_scn_cells = .
    full_join(x = dat_new, y = df_scn_cells, by = c("cell" = "cell_id")) %>% 
      dplyr::filter(in_SCN == T) -> df_
    df_[is.na(df_)] <- 0
    summary(df_)
    df_$CT = as.numeric(df_$CT)
    df_ %>% dplyr::select(-in_SCN) %>% column_to_rownames("cell")
  } -> dat_new

seq(2,22,4) %>% 
  "names<-"(., sprintf("CT%s", .)) %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(ct, ct_){
    dat_new %>% dplyr::filter(CT == ct) %>% dplyr::select(Per2, Vip) -> dat_
    dat_ %>% sum() -> n_transcript_Per2_Vip
    dat_$Per2 %>% sum() -> n_transcript_Per2
    dat_$Vip %>% sum() -> n_transcript_Vip
    
    scn_n_nuclei[[ct_]]$counts %>% dplyr::filter(region == "total") %>% .$n_nuclei -> n_nuclei
    data.frame(Per2 = n_transcript_Per2, Vip = n_transcript_Vip, Per2_Vip = n_transcript_Per2_Vip) -> df_
    df_/n_nuclei -> df_
    df_$CT = ct
    df_
  }) %>% do.call(rbind, .) -> df_plot

df_plot %>% pivot_longer(cols = -CT, names_to = "gene", values_to = "expr") %>% 
  dplyr::mutate(gene = factor(gene, levels = c("Vip", "Per2", "Per2_Vip"))) %>% 
  ggplot(aes(x = CT, y = expr, group = gene, fill = gene)) + 
  geom_col(position = "dodge") + 
  ylab("RNA concentration") + 
  scale_fill_manual(values = c(
    Vip      = "#F8766D",
    Per2     = "#00BA38",
    Per2_Vip = "#B8860B"
  )) + 
  scale_x_continuous(breaks = seq(2,22,4)) + 
  theme_classic() + 
  ylab("Number of RNA transcript per nucleus")

seq(2,22,4) %>% 
  "names<-"(., sprintf("CT%s", .)) %>% 
#  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(ct, ct_){
    dat_new %>% dplyr::filter(CT == ct) %>% dplyr::select(Per2, Vip) -> dat_
    
    dat_ %>% dplyr::filter(Per2 != 0, Vip != 0) -> dat_vip_per2
    dat_vip_per2 %>% sum() -> vip_per2_n_transcript
    nrow(dat_vip_per2) -> vip_per2_n_nucleus
    vip_per2_concen = vip_per2_n_transcript/vip_per2_n_nucleus
    
    dat_ %>% dplyr::filter(Per2 != 0) -> dat_per2
    dat_per2$Per2 %>% sum() -> per2_n_transcript
    nrow(dat_per2) -> per2_n_nucleus
    per2_concen = per2_n_transcript/per2_n_nucleus
    
    dat_ %>% dplyr::filter(Vip != 0) -> dat_vip
    dat_vip$Vip %>% sum() -> vip_n_transcript
    nrow(dat_vip) -> vip_n_nucleus
    vip_concen = vip_n_transcript/vip_n_nucleus
    
   data.frame(Per2 = per2_concen, Vip = vip_concen, Per2_Vip = vip_per2_concen) -> df_
    df_$CT = ct
    df_
  }) %>% do.call(rbind, .) -> df_plot

df_plot %>% pivot_longer(cols = -CT, names_to = "gene", values_to = "expr") %>% 
  dplyr::mutate(gene = factor(gene, levels = c("Vip", "Per2", "Per2_Vip"))) %>% 
  ggplot(aes(x = CT, y = expr, group = gene, fill = gene)) + 
  geom_col(position = "dodge") + 
  ylab("RNA concentration") + 
  scale_fill_manual(values = c(
    Vip      = "#F8766D",
    Per2     = "#00BA38",
    Per2_Vip = "#B8860B"
  )) + 
  scale_x_continuous(breaks = seq(2,22,4)) + 
  theme_classic() + 
  ylab("Number of RNA transcript per active expressing nucleus")


seq(2,22,4) %>% 
  "names<-"(., sprintf("CT%s", .)) %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(ct, ct_){
    dat_new %>% dplyr::filter(CT == ct) %>% dplyr::select(Per2, Vip) -> dat_
    
    dat_ %>% dplyr::filter(Vip != 0) -> dat_vip_pos
    dat_vip_pos$Per2 %>% sum() -> vip_pos_per2_n_transcript
    nrow(dat_vip_pos) -> vip_pos_per2_n_nucleus
    vip_pos_per2_concen = vip_pos_per2_n_transcript/vip_pos_per2_n_nucleus
    
    dat_ %>% dplyr::filter(Vip == 0) -> dat_vip_neg
    dat_vip_neg$Per2 %>% sum() -> vip_neg_per2_n_transcript
    nrow(dat_vip_neg) -> vip_neg_per2_n_nucleus
    vip_neg_per2_concen = vip_neg_per2_n_transcript/vip_neg_per2_n_nucleus
    
    data.frame(Per2 = c(vip_pos_per2_concen, vip_neg_per2_concen), Vip = c("pos", "neg")) -> df_
    df_$CT = ct
    df_
  }) %>% do.call(rbind, .) -> df_plot

df_plot %>% 
  dplyr::mutate(Vip = factor(Vip, levels = c("pos", "neg"))) %>% 
  ggplot(aes(x = CT, y = Per2, group = Vip, fill = Vip)) + 
  geom_col(position = "dodge") + 
  ylab("RNA concentration") + 
  scale_x_continuous(breaks = seq(2,22,4)) + 
  theme_classic() + 
  ylab("Mean Per2 transcript per nucleus") -> p_1

seq(2,22,4) %>% 
  "names<-"(., sprintf("CT%s", .)) %>% 
  purrr::map2(.x=., .y=names(.), .f=function(ct, ct_){
    
    dat_new %>% 
      dplyr::filter(CT == ct) %>% 
      dplyr::select(Per2, Vip) -> dat_
    
    # Vip positive nuclei
    dat_ %>% dplyr::filter(Vip != 0) -> dat_vip_pos
    vip_pos_mean = mean(dat_vip_pos$Per2)
    vip_pos_sd = sd(dat_vip_pos$Per2)
    
    # Vip negative nuclei
    dat_ %>% dplyr::filter(Vip == 0) -> dat_vip_neg
    vip_neg_mean = mean(dat_vip_neg$Per2)
    vip_neg_sd = sd(dat_vip_neg$Per2)
    
    data.frame(
      Per2 = c(vip_pos_mean, vip_neg_mean),
      SD = c(vip_pos_sd, vip_neg_sd),
      Vip = c("pos", "neg")
    ) -> df_
    
    df_$CT = ct
    df_
  }) %>% 
  do.call(rbind, .) -> df_plot

# Plot Mean ± SD
df_plot %>% 
  dplyr::mutate(Vip = factor(Vip, levels = c("pos", "neg"))) %>% 
  ggplot(aes(x = CT, y = Per2, group = Vip, fill = Vip)) + 
  geom_col(
    position = "dodge") + 
  geom_errorbar(
    aes(
      ymin = pmax(Per2 - SD, 0),
      ymax = Per2 + SD
    ),
    position = "dodge"
  ) +
  scale_x_continuous(breaks = seq(2,22,4)) + 
  theme_classic() + 
  ylab("Mean Per2 transcripts per nucleus") -> p_1 #Figure 4D

p_1

dat_new %>%
  filter(CT %in% seq(2, 22, 4)) %>%
  select(CT, Per2, Vip) %>%
  mutate(
    Vip_group = if_else(Vip > 0, "pos", "neg"),
    Vip_group = factor(Vip_group, levels = c("pos", "neg")),
    CT = factor(CT, levels = seq(2, 22, 4))
  ) -> df_plot

# Boxplot
df_plot %>%
  ggplot(aes(x = CT, y = Per2, fill = Vip_group)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.alpha = 0.2
  ) +
  scale_x_discrete(labels = function(x) paste0("CT", x)) +
  labs(
    x = "Circadian Time",
    y = "Per2 transcripts per nucleus",
    fill = "Vip"
  ) +
  theme_classic() -> p_2

p_2$data %>% 
  dplyr::filter(CT == 2) %>% 
  {
    df_ = .
    df_ %>% dplyr::filter(Vip_group == "pos") %>% dplyr::pull(Per2) -> pos
    df_ %>% dplyr::filter(Vip_group == "neg") %>% dplyr::pull(Per2) -> neg
    t.test(pos, neg)
  }



seq(2,22,4) %>% 
  "names<-"(., sprintf("CT%s", .)) %>% 
  #  .[1] %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(ct, ct_){
    dat_new %>% dplyr::filter(CT == ct) %>% dplyr::select(Per2, Vip) -> dat_
    
    dat_ %>% dplyr::filter(Vip != 0) -> dat_vip_pos
    dat_vip_pos$Per2 %>% sum() -> vip_pos_per2_n_transcript
    nrow(dat_vip_pos) -> vip_pos_per2_n_nucleus
    vip_pos_per2_concen = vip_pos_per2_n_transcript/vip_pos_per2_n_nucleus
    
    dat_ %>% dplyr::filter(Vip == 0) -> dat_vip_neg
    dat_vip_neg$Per2 %>% sum() -> vip_neg_per2_n_transcript
    nrow(dat_vip_neg) -> vip_neg_per2_n_nucleus
    vip_neg_per2_concen = vip_neg_per2_n_transcript/vip_neg_per2_n_nucleus
    
    data.frame(Per2 = c(vip_pos_per2_concen, vip_neg_per2_concen), Vip = c("pos", "neg")) -> df_
    df_$CT = ct
    df_
  }) %>% do.call(rbind, .) -> df_plot

results_list_1$CT2$Per2$spots %>% head()
results_list_1 %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(list_, CT_){
    list_[["Vip"]]$spots %>% dplyr::filter(in_SCN == T) -> df_
    df_plot = df_ %>% dplyr::mutate(y_std = scales::rescale(y_px, to = c(0,1)))
#    df_plot = df_plot %>% dplyr::filter(score >= 0.6)
#    df_plot = df_plot %>% dplyr::filter(intensity >= 500)
    df_plot %>% 
      ggplot(aes(x = -y_std, y = intensity)) + 
#      geom_density() + 
      geom_point(alpha = 0.3) + 
      ggtitle(CT_)
  }) %>% patchwork::wrap_plots(., ncol = 3)

results_list_1$CT22$Vip$spots %>% 
  dplyr::filter(in_SCN  == T) %>% 
  ggplot(aes(x = x_px, y = -y_px)) + 
  geom_point()

### new new-----

rotate_xy <- function(x, y, angle_deg, center = c(0, 0)) {
  stopifnot(
    is.numeric(x), is.numeric(y), length(x) == length(y),
    is.numeric(angle_deg), length(angle_deg) == 1L,
    is.finite(angle_deg),
    is.numeric(center), length(center) == 2L, all(is.finite(center))
  )
  
  theta <- angle_deg * pi / 180
  
  # 將旋轉中心移至原點
  dx <- x - center[1]
  dy <- y - center[2]
  
  # 旋轉，再移回原本位置
  data.frame(
    x_rot = dx * cos(theta) - dy * sin(theta) + center[1],
    y_rot = dx * sin(theta) + dy * cos(theta) + center[2]
  )
}

dat_ %>% rownames_to_column("cell") %>% 
  dplyr::mutate(cell = gsub(".+\\.(.+)", "\\1", cell)) -> dat_1

levels(dat_1$CT) %>% 
  purrr::map(function(CT_){
    ct_ = sprintf("CT%s", CT_)
    results_list_1[[ct_]]$Per2$cells %>% dplyr::filter(in_SCN == T) %>% dplyr::select(cell_id, centroid_x_px, centroid_y_px) -> all_scn_cells_df_
    dat_1 %>% dplyr::filter(CT == CT_) -> df_
    full_join(x = df_, y = all_scn_cells_df_, by = c("cell"="cell_id")) -> df_
    df_ %>% dplyr::filter(!is.na(centroid_x_px)) %>% dplyr::filter(!is.na(centroid_y_px)) -> df_
    df_$CT = CT_
    df_[is.na(df_)] <- 0
    
    if (ct_ == "CT10"){
      print("hi")
      center = c(mean(df_$centroid_x_px), mean(df_$centroid_y_px))
      
      df_[c("centroid_x_px", "centroid_y_px")] <- rotate_xy(
        x = df_$centroid_x_px,
        y = df_$centroid_y_px,
        angle_deg = 15,
        center = center
      )
    }
    
    df_ %>% dplyr::mutate(x_std = scales::rescale(x = centroid_x_px, to = c(0,1))) %>% 
      dplyr::mutate(y_std = scales::rescale(x = centroid_y_px, to = c(0,1))) -> df_
    

    
    df_plot = df_
    df_ %>% 
      ggplot(aes(y = Vip, x = -y_std)) + 
#      geom_point()
      geom_smooth()
  }) %>% patchwork::wrap_plots(., ncol = 3)
head(dat_1)
results_list_1$CT2$Per2$cells %>% dplyr::filter(cell_id == "CT2_n1820")
