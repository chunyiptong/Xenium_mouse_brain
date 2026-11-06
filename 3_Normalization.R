library(Seurat)
library(tidyverse)
library(glmGamPoi)
library(future)

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered.rda")
cm <- GetAssayData(xen_merged_filtered, assay = "Xenium", slot = "counts")
rm(xen_merged_filtered);gc()

mini <- CreateSeuratObject(counts = cm, assay = "Xenium", project = "mini")
rm(cm);gc()

save(mini, file = "~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-mini.rda")
####

rm(list=ls())
library(Seurat)
library(tidyverse)
library(glmGamPoi)
library(future)
options(future.globals.maxSize = Inf)
#plan("multisession", workers = 1)
plan(sequential)
gc()

load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-mini.rda")
gc()

mini <- SCTransform(
  mini,
  assay = "Xenium",
  vst.flavor = "v2",
  method = "glmGamPoi",
  do.correct.umi = FALSE,
  return.only.var.genes = TRUE,     # keeps memory down
  residual.features = NULL,         # residuals only for var genes
  conserve.memory = TRUE,
  verbose = TRUE,
  bin_size = 100
)
gc()
save(mini, file="~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-mini.rda")

rm(list=ls())
library(Seurat)
library(tidyverse)
library(glmGamPoi)
library(future)
load("~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-mini.rda")
mini <- RunPCA(mini, npcs = 30, features = rownames(mini))
mini <- RunUMAP(mini, dims = 1:30)
mini <- FindNeighbors(mini, reduction = "pca", dims = 1:30)
mini <- FindClusters(mini, resolution = 0.3)
DimPlot(mini, group.by = "orig.ident")
save(mini, file="~/Dropbox/Spatial_transcriptomics_new/sc/segmented/merged-sc-filtered-normalized-mini.rda")
