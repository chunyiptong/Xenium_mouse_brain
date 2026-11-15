library(tidyverse)

rpkm <- read.csv(
  gzfile('~/Downloads/LCM_seq/GSE118403_rpkm.txt.gz'),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Keep gene annotation columns separately
gene_anno <- rpkm[, 1:4]

# Expression matrix starts after the first 4 annotation columns
expr <- as.matrix(rpkm[, -(1:4)])
rownames(expr) <- rpkm$gene_id
storage.mode(expr) <- "numeric"

dim(gene_anno)
dim(expr)

all(rownames(expr) == gene_anno$gene_id)
rownames(expr) <- gene_anno$symbol
expr[1:5,]

read.csv(file = "~/Downloads/LCM_seq/label_coord_final_180522.csv", stringsAsFactors = F) -> meta_data
meta_data %>% arrange(y_re) -> meta_data

extract_expr = function(gene){
  meta_data$sample %>% 
    purrr::map(function(sample_){
      sample__ = gsub("-", "_", sample_)
     print(sample__)
     expr[gene, sample__] -> expr_
     meta_data %>% dplyr::filter(sample == sample_) %>% .$y_re -> y_re
     data.frame(expr = expr_, y_re = y_re)
    }) %>% do.call(rbind, .)
}

draw_line = function(df){
  df %>% 
    ggplot(aes(x = -y_re, y = pred_expr)) + 
    geom_line(color = "red") + 
    geom_point(data = df, aes(y = log2(expr+1), x = -y_re), alpha = 0.3) + 
#    geom_smooth(color = "red") + 
#    geom_density_2d() + 
    theme_classic() + 
    ylab("log(normalized expression + 1)")
}

#extract_expr(gene = "Vip") -> df_
#fit_linear <- lm(log2(expr + 1) ~ y_re, data = df_)
#pred = predict(fit_linear, newdata = data.frame(y_re = df_$y_re))
#df_ %>% dplyr::mutate(pred_expr = pred) -> df_
#draw_line(df = df_)

extract_expr(gene = "Vip") -> df_
fit_quad <- lm(log2(expr + 1) ~ y_re + I(y_re^2), data = df_)
pred = predict(fit_quad, newdata = data.frame(y_re = df_$y_re))
df_ %>% dplyr::mutate(pred_expr = pred) -> df_
draw_line(df = df_)

extract_expr(gene = "Gpr101") -> df_
fit_quad <- lm(log2(expr + 1) ~ y_re + I(y_re^2), data = df_)
pred = predict(fit_quad, newdata = data.frame(y_re = df_$y_re))
df_ %>% dplyr::mutate(pred_expr = pred) -> df_
draw_line(df = df_)

extract_expr(gene = "Vip") -> df_
# Fit smoothing spline
fit_spline <- smooth.spline(
  x = df_$y_re,
  y = log2(df_$expr + 1),
  df = 3
)
# Predict fitted expression
pred <- predict(
  fit_spline,
  x = df_$y_re
)$y
# Add predictions to dataframe
df_ %>%
  dplyr::mutate(pred_expr = pred) -> df_
# Plot
draw_line(df = df_) #Figure 3B


extract_expr(gene = "Gpr101") -> df_
# Fit smoothing spline
fit_spline <- smooth.spline(
  x = df_$y_re,
  y = log2(df_$expr + 1),
  df = 3
)
# Predict fitted expression
pred <- predict(
  fit_spline,
  x = df_$y_re
)$y
# Add predictions to dataframe
df_ %>%
  dplyr::mutate(pred_expr = pred) -> df_
# Plot
draw_line(df = df_) #Figure 3B

