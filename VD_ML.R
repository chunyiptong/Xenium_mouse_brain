library(tidyverse)

source("~/Dropbox/Spatial_transcriptomics_new/github/segmented/functions.R")
load("~/Downloads/New_SVG.rda")
load("~/Downloads/New_TVG.rda")
load("~/Downloads/New_STVG.rda")

colnames(norm_counts_matrix) -> scn_cells
xy_coordinates %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_,CT_){
#    df_ %>% head()
    df_ %>% dplyr::filter(cell %in% scn_cells) -> df_
#    dim(df_)
    df_ %>% dplyr::mutate(sample = CT_) -> df_
#    head(df_)
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> xy_coordinates_scn_df_

standardize_coords_by_slice <- function(x, y, slice_id, center = FALSE) {
  # x, y: numeric vectors
  # slice_id: factor or character vector (length n)
  # center: if TRUE, also center to [-0.5, 0.5]
  
  stopifnot(length(x) == length(y))
  stopifnot(length(x) == length(slice_id))
  
  slice_id <- as.factor(slice_id)
  n <- length(x)
  
  x_std <- numeric(n)
  y_std <- numeric(n)
  
  for (s in levels(slice_id)) {
    idx <- which(slice_id == s)
    
    x_s <- x[idx]
    y_s <- y[idx]
    
    # handle constant edge case
    range_x <- max(x_s) - min(x_s)
    range_y <- max(y_s) - min(y_s)
    
    if (range_x == 0) {
      x_std[idx] <- 0.5
    } else {
      x_std[idx] <- (x_s - min(x_s)) / range_x
    }
    
    if (range_y == 0) {
      y_std[idx] <- 0.5
    } else {
      y_std[idx] <- (y_s - min(y_s)) / range_y
    }
    
    if (center) {
      x_std[idx] <- x_std[idx] - 0.5
      y_std[idx] <- y_std[idx] - 0.5
    }
  }
  
  data.frame(
    x_std = x_std,
    y_std = y_std
  )
}

coordinate_aligned=standardize_coords_by_slice(xy_coordinates_scn_df_$x, 
                                               xy_coordinates_scn_df_$y,
                                               xy_coordinates_scn_df_$sample,
                                               center=TRUE)
coordinate_aligned %>% dim() # 4915 scn cells
cbind(xy_coordinates_scn_df_, coordinate_aligned) -> xy_coordinates_scn_df_

extract_expr = function(gene_, xy_df_){
  cells_ = xy_df_$cell
  norm_counts_matrix[gene_,cells_] %>% as.numeric() -> expr_
  xy_df_ %>% dplyr::mutate(expr = expr_)
}

draw_line_VD = function(df){
  df %>% 
    dplyr::mutate(CT = factor(CT, levels = c("CT-8", "CT-20"))) -> df_plot
  
  df_plot %>% 
    ggplot(aes(x = y_std, y = pred_expr_VD, color = CT)) + 
    geom_line() + 
#    geom_point(data = df_plot %>% dplyr::filter(expr != 0), aes(x = y_std, y = log2(expr+1), color = CT), alpha = 0.1) +
#    geom_smooth(se = F) + 
    #    geom_density_2d() + 
    theme_classic() + 
    scale_x_continuous(
      breaks = c(-0.5, 0.5),
      labels = c("Ventral", "Dorsal")
    ) +
    labs(
      x = NULL,
      y = "log(normalized expression + 1)"
    )
}

draw_line_ML = function(df){
  df %>% 
    dplyr::mutate(CT = factor(CT, levels = c("CT-8", "CT-20"))) -> df_plot
  
  df_plot %>% 
    ggplot(aes(x = x_std, y = pred_expr_ML, color = CT)) + 
    geom_line() + 
    geom_point(data = df_plot %>% dplyr::filter(expr != 0), aes(x = x_std, y = log2(expr+1), color = CT), alpha = 0.1) +
#    geom_smooth(se = F) + 
    #    geom_density_2d() + 
    theme_classic() + 
    scale_x_continuous(
      breaks = c(0, 0.5),
      labels = c("Medial", "Lateral")
    ) +
    labs(
      x = NULL,
      y = "log(normalized expression + 1)"
    )
}

draw_VD_MD = function(gene){
  extract_expr(gene_ = gene, xy_df_ = xy_coordinates_scn_df_) -> df_
  df_$x_std = abs(df_$x_std)
  df_ %>% dplyr::filter(grepl("CT-8", sample)) -> df_CT8
  df_ %>% dplyr::filter(grepl("CT-20", sample)) -> df_CT20
  
  fit_quad_CT8_VD <- lm(log2(expr + 1) ~ y_std + I(y_std^2), data = df_CT8)
  fit_quad_CT20_VD <- lm(log2(expr + 1) ~ y_std + I(y_std^2), data = df_CT20)
  pred_CT8_VD = predict(fit_quad_CT8_VD, newdata = data.frame(y_std = df_CT8$y_std))
  pred_CT20_VD = predict(fit_quad_CT20_VD, newdata = data.frame(y_std = df_CT20$y_std))
  rbind(df_CT8, df_CT20) %>% dplyr::mutate(pred_expr_VD = c(pred_CT8_VD, pred_CT20_VD)) %>% dplyr::mutate(CT = gsub("(CT-\\d+?)-.+", "\\1", sample))-> df_
  
  fit_quad_CT8_ML <- lm(log2(expr + 1) ~ x_std + I(x_std^2), data = df_CT8)
  fit_quad_CT20_ML <- lm(log2(expr + 1) ~ x_std + I(x_std^2), data = df_CT20)
  pred_CT8_ML = predict(fit_quad_CT8_ML, newdata = data.frame(x_std = df_CT8$x_std))
  pred_CT20_ML = predict(fit_quad_CT20_ML, newdata = data.frame(x_std = df_CT20$x_std))
  df_ %>% dplyr::mutate(pred_expr_ML = c(pred_CT8_ML, pred_CT20_ML)) -> df_
  
  draw_line_VD(df = df_) -> p1
  draw_line_ML(df = df_) -> p2
  patchwork::wrap_plots(p1, p2, ncol = 2, guides = "collect") + patchwork::plot_annotation(title = gene)  
}


draw_VD_MD_1 <- function(gene, spline_df = 5, return_details = FALSE) {
  
  # Extract expression and coordinates
  df_ <- extract_expr(
    gene_ = gene,
    xy_df_ = xy_coordinates_scn_df_
  )
  
  df_$x_std <- abs(df_$x_std)
  
  df_ <- df_ %>%
    dplyr::filter(
      is.finite(x_std),
      is.finite(y_std),
      is.finite(expr),
      expr >= 0
    ) %>%
    dplyr::mutate(
      CT = dplyr::case_when(
        grepl("^CT-8-", sample) ~ "CT-8",
        grepl("^CT-20-", sample) ~ "CT-20",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(CT))
  
  df_CT8 <- df_ %>%
    dplyr::filter(CT == "CT-8")
  
  df_CT20 <- df_ %>%
    dplyr::filter(CT == "CT-20")
  
  # Helper function for spline fitting
  fit_spline <- function(dat, coord) {
    
    stats::smooth.spline(
      x = dat[[coord]],
      y = log2(dat$expr + 1),
      df = spline_df#,
#      cv = NA
    )
  }
  
  # 1. Ventral-Dorsal spline models
  fit_CT8_VD <- fit_spline(df_CT8, "y_std")
  fit_CT20_VD <- fit_spline(df_CT20, "y_std")
  
  # 2. Medial-Lateral spline models
  fit_CT8_ML <- fit_spline(df_CT8, "x_std")
  fit_CT20_ML <- fit_spline(df_CT20, "x_std")
  
  # 3. Predict fitted expression
  df_CT8 <- df_CT8 %>%
    dplyr::mutate(
      pred_expr_VD = predict(
        fit_CT8_VD, x = y_std
      )$y,
      pred_expr_ML = predict(
        fit_CT8_ML, x = x_std
      )$y
    )
  
  df_CT20 <- df_CT20 %>%
    dplyr::mutate(
      pred_expr_VD = predict(
        fit_CT20_VD, x = y_std
      )$y,
      pred_expr_ML = predict(
        fit_CT20_ML, x = x_std
      )$y
    )
  
  # 4. Combine data
  df_ <- dplyr::bind_rows(df_CT8, df_CT20)
  
  # 5. Plot
  p1 <- draw_line_VD(df = df_)
  p2 <- draw_line_ML(df = df_)
  
  final_plot <- patchwork::wrap_plots(
    p1, p2,
    ncol = 2,
    guides = "collect"
  ) +
    patchwork::plot_annotation(title = gene)
  
  # Optional: return fitted models and data
  if (return_details) {
    return(list(
      plot = final_plot,
      data = df_,
      models = list(
        CT8_VD = fit_CT8_VD,
        CT20_VD = fit_CT20_VD,
        CT8_ML = fit_CT8_ML,
        CT20_ML = fit_CT20_ML
      )
    ))
  }
  
  return(final_plot)
}


draw_VD_MD_2 <- function(
    gene,
    spline_df = 5,
    return_details = FALSE,
    n_grid = 201
) {
  
  library(dplyr)
  library(ggplot2)
  
  # 1. Extract expression and coordinates
  df_ <- extract_expr(
    gene_ = gene,
    xy_df_ = xy_coordinates_scn_df_
  )
  
  df_ <- df_ %>%
    mutate(
      x_std = abs(x_std),
      CT = case_when(
        grepl("^CT-8-", sample) ~ "CT-8",
        grepl("^CT-20-", sample) ~ "CT-20",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      is.finite(x_std),
      is.finite(y_std),
      is.finite(expr),
      expr >= 0,
      !is.na(CT)
    ) %>%
    mutate(
      CT = factor(CT, levels = c("CT-8", "CT-20"))
    )
  
  if (!all(c("CT-8", "CT-20") %in% df_$CT)) {
    stop("Both CT-8 and CT-20 are required.")
  }
  
  # 2. Split by biological replicate
  sample_list <- split(df_, df_$sample)
  
  # 3. Fit replicate-level splines on each axis
  fit_by_sample <- function(coord, axis) {
    
    # Determine common spatial range across samples
    bounds <- df_ %>%
      group_by(sample) %>%
      summarise(
#        lower = quantile(.data[[coord]], 0),
#        upper = quantile(.data[[coord]], 1),
        lower = min(.data[[coord]]),
        upper = max(.data[[coord]]),
        .groups = "drop"
      )
#    print(bounds)
    
    lower <- max(bounds$lower)
    upper <- min(bounds$upper)
    
    if (lower >= upper) {
      stop(paste("No common spatial range for", axis))
    }
    
    grid <- seq(lower, upper, length.out = n_grid)
#    print(grid)
    
    # Fit one spline per sample
    models <- lapply(sample_list, function(dat) {
      
      stats::smooth.spline(
        x = dat[[coord]],
        y = log2(dat$expr + 1),
        df = spline_df
      )
      
    })
    
    # Predict each replicate on the same grid
    predictions <- bind_rows(
      lapply(names(sample_list), function(s) {
        
        dat <- sample_list[[s]]
        
        data.frame(
          sample = s,
          CT = as.character(dat$CT[1]),
          axis = axis,
          pos = grid,
          fitted = predict(
            models[[s]],
            x = grid
          )$y
        )
        
      })
    )
    
    # Mean and SD across biological replicates
    summary <- predictions %>%
      group_by(axis, CT, pos) %>%
      summarise(
        mean_pred = mean(fitted),
        sd_pred = sd(fitted),
        n_reps = n_distinct(sample),
        .groups = "drop"
      ) %>%
      mutate(
        CT = factor(CT, levels = c("CT-8", "CT-20"))
      )
    
    list(
      predictions = predictions,
      summary = summary,
      models = models,
      grid = grid
    )
  }
  
  # 4. Replicate-level spline fitting
  VD <- fit_by_sample("y_std", "VD")
  ML <- fit_by_sample("x_std", "ML")
  
  # 5. Plotting helper
  draw_axis <- function(res, axis) {
    
    pred <- res$predictions
    summ <- res$summary
    
    pred$CT <- factor(
      pred$CT,
      levels = c("CT-8", "CT-20")
    )
    
    p <- ggplot() +
      
      # Mean +/- SD ribbon
      geom_ribbon(
        data = summ,
        aes(
          x = pos,
          ymin = mean_pred - sd_pred,
          ymax = mean_pred + sd_pred,
          fill = CT,
          group = CT
        ),
        alpha = 0.18,
        color = NA
      ) +
      
#      # Individual replicate curves
#      geom_line(
#        data = pred,
#        aes(
#          x = pos,
#          y = fitted,
#          color = CT,
#          group = sample
#        ),
#        linewidth = 0.5,
#        alpha = 0.4
#      ) +
      
      # Mean fitted curve
      geom_line(
        data = summ,
        aes(
          x = pos,
          y = mean_pred,
          color = CT,
          group = CT
        ),
        linewidth = 1.2
      ) +
      
      theme_classic() +
      guides(fill = "none") +
      labs(
        x = NULL,
        y = "log2(normalized expression + 1)",
        color = "CT"
      )
    
    if (axis == "VD") {
      
      p <- p +
        scale_x_continuous(
          limits = c(-0.5, 0.5),
          breaks = c(-0.5, 0.5),
          labels = c("Ventral", "Dorsal")
        )
      
    } else {
      
      p <- p +
        scale_x_continuous(
          limits = range(summ$pos),
          breaks = range(summ$pos),
          labels = c("Medial", "Lateral")
        )
    }
    
    return(p)
  }
  
  # 6. Draw both axes
  p1 <- draw_axis(VD, "VD")
  p2 <- draw_axis(ML, "ML")
  
  final_plot <- patchwork::wrap_plots(
    p1, p2,
    ncol = 2,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = gene
    )
  
  # 7. Return details if requested
  if (return_details) {
    
    # Fitted values at original cell coordinates
    df_fitted <- bind_rows(
      lapply(names(sample_list), function(s) {
        
        dat <- sample_list[[s]]
        
        dat$pred_expr_VD <- predict(
          VD$models[[s]], x = dat$y_std
        )$y
        
        dat$pred_expr_ML <- predict(
          ML$models[[s]], x = dat$x_std
        )$y
        
        dat
      })
    )
    
    return(list(
      plot = final_plot,
      data = df_fitted,
      predictions = bind_rows(
        VD$predictions,
        ML$predictions
      ),
      summary = bind_rows(
        VD$summary,
        ML$summary
      ),
      models = list(
        VD = VD$models,
        ML = ML$models
      )
    ))
  }
  
  return(final_plot)
}



draw_VD_MD_2(gene = "Vip", spline_df = 5)#[[1]]$data
draw_VD_MD_1(gene = "Gad1")
draw_VD_MD(gene = "Calb1")
draw_VD_MD(gene = "Slc39a6")
draw_VD_MD(gene = "Per2") 
draw_VD_MD(gene = "Lhx1")
draw_VD_MD_1(gene = "Rgs16")

draw_VD_MD_2(gene = "Cacna1g", spline_df = 3)
draw_VD_MD_2(gene = "Kdm6a", spline_df = 3)
draw_VD_MD_2(gene = "Vhl", spline_df = 3)

draw_VD_MD_2(gene = "Per2", spline_df = 4) #Figure 4C
draw_VD_MD_2(gene = "Lhx1", spline_df = 5) #Figure 4E
draw_VD_MD_2(gene = "Rgs16", spline_df = 4) #Figure 4F


list(
  CT8_SVG = CT8_SVG_1,
  CT20_SVG = CT20_SVG_1,
  TVG = TVG,
  STVG = STVG_1
) %>% 
  ggvenn::list_to_data_frame() -> SVG_TVG_STVG_df


SVG_TVG_STVG_df %>% dplyr::filter(TVG == F, CT8_SVG == F, CT20_SVG == F, STVG == T) %>% .$key %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    draw_VD_MD(gene = gene_) 
  }) -> p_list

pdf("~/Downloads/p_list_101.pdf", width = 10.51, height = 4.15)
p_list
dev.off()

map_feature_plot(gene_ = "Cacna1g")
draw_VD_MD(gene = "Cacna1g")

map_feature_plot(gene_ = "Kdm6a")
draw_VD_MD(gene = "Kdm6a")
DEG_plot(gene_ = "Kdm6a")


map_feature_plot(gene_ = "Vhl")
draw_VD_MD(gene = "Vhl")
DEG_plot(gene_ = "Vhl")

SVG_TVG_STVG_df %>% dplyr::filter(TVG == T, CT8_SVG == T, CT20_SVG == T, STVG == T) %>% .$key %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    draw_VD_MD(gene = gene_) 
  }) -> p_list

pdf("~/Downloads/p_list_102.pdf", width = 10.51, height = 4.15)
p_list
dev.off()

source("~/Dropbox/Spatial_transcriptomics_new/github/segmented/classify_stvg_axis.R")
draw_VD_MD(gene = "Rgs16")[[2]]$data -> df_
prepare_curve(df = df_, ct = "CT-8", axis = "ML") -> df_CT8
prepare_curve(df = df_, ct = "CT-20", axis = "ML") -> df_CT20
classify_shape(pos = df_CT8$pos, fitted = df_CT8$fitted)
classify_shape(pos = df_CT20$pos, fitted = df_CT20$fitted)

STVG_1 %>% length()
STVG_1 %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(gene_){
    draw_VD_MD_1(gene = gene_)    
  }) -> p_list_STVG

pdf("~/Downloads/p_list_STVGs.pdf", width = 10.51, height = 4.15)
p_list_STVG
dev.off()

#p_list_STVG[c("Rgs16", "Vip")] %>% 
p_list_STVG %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(p,gene_){
    p[[1]]$data -> df_VD
    p[[2]]$data -> df_ML
    
    prepare_curve(df = df_VD, ct = "CT-8", axis = "VD") -> df_VD_CT8
    prepare_curve(df = df_VD, ct = "CT-20", axis = "VD") -> df_VD_CT20
    
    classify_shape(pos = df_VD_CT8$pos, fitted = df_VD_CT8$fitted) -> df_VD_CT8_shape
    classify_shape(pos = df_VD_CT20$pos, fitted = df_VD_CT20$fitted) -> df_VD_CT20_shape
    
    prepare_curve(df = df_ML, ct = "CT-8", axis = "ML") -> df_ML_CT8
    prepare_curve(df = df_ML, ct = "CT-20", axis = "ML") -> df_ML_CT20
    
    classify_shape(pos = df_ML_CT8$pos, fitted = df_ML_CT8$fitted) -> df_ML_CT8_shape
    classify_shape(pos = df_ML_CT20$pos, fitted = df_ML_CT20$fitted) -> df_ML_CT20_shape
    
    sprintf("VD: %s (CT8), %s (CT20); ML %s (CT8), %s (CT20)", df_VD_CT8_shape, df_VD_CT20_shape, df_ML_CT8_shape, df_ML_CT20_shape)
    data.frame(gene = gene_, VD_CT8_shape = df_VD_CT8_shape, VD_CT20_shape = df_VD_CT20_shape, ML_CT8_shape = df_ML_CT8_shape, ML_CT20_shape = df_ML_CT20_shape)
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> df_shape


