library(tidyverse)

#load('~/Downloads/SCN_count_matrix.rda')
#rm(counts_matrix, df_count_mat_list, df_count_mat_norm_list)
#load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-xy-coordinates.rda")
#
#load("~/Downloads/SCN_core_shell_matrix_norm.Rdata")
#SCN_core_shell_matrix %>% 
#  purrr::map2(.x=.,.y=names(.),.f=function(df_, name_){
#    df_ %>% colSums()
#  }) %>% unlist() %>% median() -> median_
#
#SCN_core_shell_matrix %>% 
#  purrr::map2(.x=.,.y=names(.),.f=function(df_, name_){
#    colnames(df_) = sprintf("%s_%s", name_, colnames(df_))
#    df_ %>% rownames_to_column("gene")
#  }) %>% purrr::reduce(., inner_join, by = "gene") %>% 
#  column_to_rownames("gene") %>% 
#  dplyr::mutate(across(everything(), function(x){(x/sum(x))*median_})) -> SCN_core_shell_matrix_norm
#rm(median_, SCN_core_shell_matrix)
#
#list.files(path = "~/Dropbox/Spatial_transcriptomics_new/rda/Brain_region_annotation/SCN_SON_PVN", pattern = "_regions_v5\\.csv", full.names = T) %>% 
#  "names<-"(., gsub("/.+/(.+)_regions.+", "\\1", .) %>% gsub("_", "-", .)) %>% 
#  purrr::map(function(file_){
#    read_csv(file_, skip = 2) -> df_
#  })-> polygon_
#
#save(norm_counts_matrix, polygon_, SCN_core_shell_matrix_norm, xy_coordinates, file = "~/Downloads/SCN_feature_plot_dat.rda")

load("~/Downloads/SCN_feature_plot_dat.rda")
map_feature_plot = function(gene_){
  norm_counts_matrix[gene_, ] %>% as.numeric() -> gene_expr
  upper = quantile(gene_expr, 0.99, na.rm = TRUE)
  c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") %>% 
    "names<-"(.,.) %>% 
    purrr::map(function(sample_){
      norm_counts_matrix[gene_, ] %>% 
        t() %>% 
        as.data.frame() %>% "colnames<-"(., "normalized_expr") %>% 
        rownames_to_column("cell") -> df_
      xy_coordinates[[sample_]] -> df_xy
      inner_join(x = df_, y = df_xy, by = "cell") -> df_plot
      
      print(upper)
      df_plot %>% 
        dplyr::arrange(normalized_expr) -> df_plot
      
      df_plot %>% dplyr::mutate(normalized_expr_1 = case_when(
        normalized_expr > upper ~ upper,
        TRUE ~ normalized_expr
      )) -> df_plot
      #        upper = quantile(df_plot$normalized_expr, 0.99, na.rm = TRUE)
      df_plot %>% 
        ggplot(aes(x = x, y = y, color = normalized_expr, size = normalized_expr_1)) + 
        #          geom_point(size = 1) + 
        geom_point() + 
        scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, upper), oob = scales::squish) + 
        scale_size_continuous(range = c(1,4)) +
        theme_void() + 
        labs(color = "Normalized\nExpression") +
        ggtitle(sprintf("%s", sample_)) -> p
      
      polygon_[[sample_]]$Selection %>% unique() %>% {.[grepl("shell|core", .)]} -> regions_
      for (region_ in regions_){
        polygon_[[sample_]] %>% dplyr::filter(Selection == region_) -> polygon_sub
        p = p + geom_polygon(data = polygon_sub, aes(x = X, y = Y), fill = NA, color = "black", linewidth = 0.3, inherit.aes = F)
      }
      
      return(p)
      
    }) %>% patchwork::wrap_plots(., ncol = 3, guides = "collect") -> p 
  return(p)
}

DEG_plot = function(gene_){
  SCN_core_shell_matrix_norm[gene_, ] %>% 
    pivot_longer(cols = everything(), names_to = "group", values_to = "norm_expr") %>% 
    dplyr::mutate(region = gsub("(.+)_CT.+", "\\1", group)) %>% 
    dplyr::mutate(CT = gsub("(.+)_(CT-\\d+?)-.+", "\\2", group)) %>% 
    dplyr::mutate(region = gsub("SCN_", "", region)) %>% 
    dplyr::mutate(CT = gsub("-", "", CT)) %>% 
    dplyr::mutate(group_1 = sprintf("%s_%s", CT, region)) -> df_plot
  
  df_plot$group_1 = factor(df_plot$group_1, levels = c("CT8_core", "CT8_shell", "CT20_core", "CT20_shell"))
  df_plot$CT = factor(df_plot$CT, levels = c("CT8", "CT20"))
  
  df_plot %>% 
    ggplot(aes(x = group_1, y = norm_expr, shape = region, color = CT)) + 
    geom_point(size = 4, alpha = 0.7) + 
    stat_summary(
      aes(group = group_1, color = CT),
      fun = mean,
      geom = "crossbar",
      width = 0.4,
      #        color = "black",
      linewidth = 0.4
    ) +
    theme_classic() + 
    ggtitle(gene_) -> p
  return(p)
}

#load("~/Downloads/scRNA_dat.rda")
#rm(JTKi)
#scRNA_mat_norm_mean = mat_norm_mean
#scRNA_mat_norm_se = mat_norm_se
#rm(mat_norm_mean, mat_norm_se)
#gc()
#save(scRNA_mat_norm_mean, scRNA_mat_norm_se, file = "~/Downloads/scRNA_SCN_line_plot_dat.rda")
#rm(scRNA_mat_norm_mean, scRNA_mat_norm_se)

#load("~/Downloads/scRNA_SCN_line_plot_dat.rda")
#list.files("~/Dropbox/Spatial_transcriptomics_new/external_data/Wen_mouse_SCN_circadian/GSE117295_RAW", full.names = T) %>% 
#  .[grep("CT", .)] %>% 
#  "names<-"(., gsub("/.+/.+(CT\\d+)\\.csv\\.gz", "\\1", .)) %>% 
#  #  .[1] %>% 
#  purrr::map2(.x=.,.y=names(.),.f=function(file_, CT_){
#    print(file_)
#    read.csv(file = file_, row.names = 1, header = T) -> df_
#  }) -> mat_list
#
#mat_list %>% 
#  purrr::map(function(df_){
#    df_ %>% dplyr::mutate(across(everything(), function(x){(x/sum(x))*10^6})) -> df_
#  }) -> mat_norm_list
#rm(mat_list);gc()
#save(scRNA_mat_norm_mean, scRNA_mat_norm_se, mat_norm_list, file = "~/Downloads/scRNA_SCN_line_plot_dat.rda")

load("~/Downloads/scRNA_SCN_line_plot_dat.rda")
scRNA_lineplot = function(gene_){
  
  scRNA_mat_norm_mean %>% dplyr::filter(gene == gene_) %>% 
    dplyr::select(-gene) %>% 
    pivot_longer(cols = everything(), names_to = "CT", values_to = "norm_expr") -> df_expr
  
  scRNA_mat_norm_se %>% dplyr::filter(gene == gene_) %>% 
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
    geom_point() + 
    geom_ribbon(aes(ymin = norm_expr - se, ymax = norm_expr + se), 
                alpha = 0.2) + 
    theme_classic() + 
    ylab("Normalized expression") -> p
  
  p + scale_x_continuous(breaks = seq(14,58,4)) -> p
  return(p)
}

scRNA_lineplot_1 = function(data, CT_start, CT_end, gene_){
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
  
  p + scale_x_continuous(breaks = seq(14,58,4)) + ggtitle(gene_) -> p1
  return(p1)
}


scRNA_lineplot_2 = function(data, CT_start, CT_end, gene_){
  data %>% dplyr::filter(CT <= CT_end, CT >= CT_start) -> df_plot
  df_plot %>% dplyr::mutate(CT_1 = case_when(
    CT > CT_end-12 ~ CT-24,
    TRUE ~ CT
  )) -> df_plot
  
  CT_1_min = min(df_plot$CT_1)
  if (CT_1_min > 24){
    df_plot$CT_1 = df_plot$CT_1 - 24
  }
  
  night_df <- data.frame(
    xmin = 12,
    xmax = 22
  )
  
  df_plot %>% 
    ggplot(aes(x = CT_1, y = norm_expr)) + 
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
  
#  p + scale_x_continuous(breaks = seq(2,22,4)) -> p1
  
  df_plot %>% dplyr::filter(CT_1 == 6) %>% .$CT %>% sprintf("CT%s", .) -> CT_day
  df_plot %>% dplyr::filter(CT_1 == 18) %>% .$CT %>% sprintf("CT%s", .) -> CT_night
  compare_day_night(gene = gene_, day_CT = CT_day, night_CT = CT_night) -> pval
  p + scale_x_continuous(breaks = seq(2,22,4)) + ggtitle(sprintf("%s: CT6 vs CT18 pval=%s", gene_, pval)) -> p1
  return(p1)
}

scRNA_lineplot_2_new = function(data, CT_start, CT_end, gene_){
  data %>% dplyr::filter(CT <= CT_end, CT >= CT_start) -> df_plot
  df_plot$CT_1 = df_plot$CT
  while (min(df_plot$CT_1) > 0){
    df_plot$CT_1 = df_plot$CT_1 - 24
  }
  df_plot %>% dplyr::mutate(CT_1 = case_when(
    CT_1 < 0 ~ CT_1 + 24, 
    TRUE ~ CT_1
  )) -> df_plot
  
  night_df <- data.frame(
    xmin = 12,
    xmax = 22
  )
  
  df_plot %>% 
    ggplot(aes(x = CT_1, y = norm_expr)) + 
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
  
  #  p + scale_x_continuous(breaks = seq(2,22,4)) -> p1
  
  df_plot %>% dplyr::filter(CT_1 == 6) %>% .$CT %>% sprintf("CT%s", .) -> CT_day
  df_plot %>% dplyr::filter(CT_1 == 18) %>% .$CT %>% sprintf("CT%s", .) -> CT_night
  compare_day_night(gene = gene_, day_CT = CT_day, night_CT = CT_night) -> pval
  p + scale_x_continuous(breaks = seq(2,22,4)) + ggtitle(sprintf("%s: CT6 vs CT18 pval=%s", gene_, pval)) -> p1
  return(p1)
}

compare_day_night = function(gene, day_CT, night_CT){
  mat_norm_list[[day_CT]][gene, ] %>% as.numeric() -> day_expr
  mat_norm_list[[night_CT]][gene, ] %>% as.numeric() -> night_expr
  t.test(x = day_expr, y = night_expr)$p.value %>% sprintf("%.2e", .) -> pval
  return(pval)
}

scRNA_lineplot_3 = function(data, CT_start, CT_end, gene_){
  data %>% dplyr::filter(CT <= CT_end, CT >= CT_start) -> df_plot
  
  
  function_1 = function(CT_){
    while (CT_ > 24){
      CT_ = CT_ - 24
    }
    return(CT_)
  }
  
  df_plot$CT %>% 
    purrr::map(function(x){
      function_1(CT_ = x)
    }) %>% unlist() -> df_plot$CT_1
  
  night_df <- data.frame(
    xmin = 12,
    xmax = 22
  )
  
  df_plot %>% 
    ggplot(aes(x = CT_1, y = norm_expr)) + 
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
  
  #  p + scale_x_continuous(breaks = seq(2,22,4)) -> p1
  
  df_plot %>% dplyr::filter(CT_1 == 6) %>% .$CT %>% sprintf("CT%s", .) -> CT_day
  df_plot %>% dplyr::filter(CT_1 == 18) %>% .$CT %>% sprintf("CT%s", .) -> CT_night
  compare_day_night(gene = gene_, day_CT = CT_day, night_CT = CT_night) -> pval
  p + scale_x_continuous(breaks = seq(2,22,4)) + ggtitle(sprintf("%s: CT6 vs CT18 pval=%s", gene_, pval)) -> p1
  return(p1)
}
