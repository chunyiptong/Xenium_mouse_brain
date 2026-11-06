#ml GCC/13.2.0 OpenMPI/4.1.6 R/4.4.1

.libPaths("/scratch/user/tongtong1993/project/R/4.4")
library(Seurat)
library(tidyverse)
library(future)
library(glmGamPoi)
options(future.globals.maxSize = Inf)
plan(sequential)

meta_data = read_csv("~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/metadata/WMB-10X/20241115/views/cell_metadata_with_cluster_annotation.csv")
meta_data %>% as.data.frame() %>% "rownames<-"(., .$cell_label) -> meta_data

class_ = meta_data$class %>% unique()
class_list = list(
		  Glutamatergic = class_[grepl("Glut", class_)],
		  GABAergic = class_[grepl("GABA", class_)],
		  Other_neuronal = class_[grepl("21|22", class_)],
		  Glial = class_[grepl("30|31|32", class_)],
		  Vascular = class_[grepl("33", class_)],
		  Immune = class_[grepl("34", class_)]
		  )

list.files("~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices", pattern="\\.rds", full=T, recursive=T) -> rds_
rds_ %>% {.[grepl("Isocortex|CTxsp|OLF|TH|HY|STR|PAL", .)]} -> rds_
gsub("/.+/(.+)\\.rds", "\\1", rds_) -> names(rds_)

make_seurat_obj = function(raw_counts, meta_data){
	counts_sparse = readRDS(raw_counts)
	seurat_obj <- CreateSeuratObject(counts = counts_sparse)
	rm(counts_sparse);gc()

	all(colnames(seurat_obj) %in% rownames(meta_data))
	cbind(seurat_obj@meta.data, meta_data[colnames(seurat_obj), c("feature_matrix_label", "class")]) -> meta_data_1

	meta_data_1 %>%
		dplyr::mutate(custom_class = case_when(
						       class %in% class_list$Glutamatergic ~ "Glutamatergic",
						       class %in% class_list$GABAergic ~ "GABAergic",
						       class %in% class_list$Other_neuronal ~ "Other_neuronal",
						       class %in% class_list$Glial ~ "Glial",
						       class %in% class_list$Vascular ~ "Vascular",
						       class %in% class_list$Immune ~ "Immune",
						       TRUE ~ class
						       )) -> meta_data_1 
	seurat_obj@meta.data <- meta_data_1

	seurat_obj <- SCTransform(
				  seurat_obj,
				  assay = "RNA",
				  vst.flavor = "v2",
				  method = "glmGamPoi",
				  do.correct.umi = FALSE,
				  return.only.var.genes = TRUE,     # keeps memory down
				  residual.features = NULL,         # residuals only for var genes
				  conserve.memory = TRUE,
				  verbose = TRUE,
				  bin_size = 50
				 )
	DefaultAssay(seurat_obj) <- "SCT"
	seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
	seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
	return(seurat_obj)
}

names(rds_) -> rds_names_

make_seurat_obj(raw_counts=rds_["WMB-10Xv2-HY-raw"], meta_data=meta_data) -> WMB_10Xv2_HY_sc
DimPlot(WMB_10Xv2_HY_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv2-HY-raw.png")
save(WMB_10Xv2_HY_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv2_HY_sc")
rm(WMB_10Xv2_HY_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv2-Isocortex-1-raw"], meta_data=meta_data) -> WMB_10Xv2_Isocortex_1_sc
DimPlot(WMB_10Xv2_Isocortex_1_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv2-Isocortex-1-raw.png")
save(WMB_10Xv2_Isocortex_1_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv2_Isocortex-1_sc")
rm(WMB_10Xv2_Isocortex_1_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv2-Isocortex-2-raw"], meta_data=meta_data) -> WMB_10Xv2_Isocortex_2_sc
DimPlot(WMB_10Xv2_Isocortex_2_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv2-Isocortex-2-raw.png")
save(WMB_10Xv2_Isocortex_2_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv2_Isocortex-2_sc")
rm(WMB_10Xv2_Isocortex_2_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv2-Isocortex-3-raw"], meta_data=meta_data) -> WMB_10Xv2_Isocortex_3_sc
DimPlot(WMB_10Xv2_Isocortex_3_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv2-Isocortex-3-raw.png")
save(WMB_10Xv2_Isocortex_3_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv2_Isocortex-3_sc")
rm(WMB_10Xv2_Isocortex_3_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv2-Isocortex-4-raw"], meta_data=meta_data) -> WMB_10Xv2_Isocortex_4_sc
DimPlot(WMB_10Xv2_Isocortex_4_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv2-Isocortex-4-raw.png")
save(WMB_10Xv2_Isocortex_4_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv2_Isocortex-4_sc")
rm(WMB_10Xv2_Isocortex_4_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv2-OLF-raw"], meta_data=meta_data) -> WMB_10Xv2_OLF_sc
DimPlot(WMB_10Xv2_OLF_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv2-OLF-raw.png")
save(WMB_10Xv2_OLF_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv2_OLF_sc")
rm(WMB_10Xv2_OLF_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv2-TH-raw"], meta_data=meta_data) -> WMB_10Xv2_TH_sc
DimPlot(WMB_10Xv2_TH_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv2-TH-raw.png")
save(WMB_10Xv2_TH_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv2_TH_sc")
rm(WMB_10Xv2_TH_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv3-HY-raw"], meta_data=meta_data) -> WMB_10Xv3_HY_sc
DimPlot(WMB_10Xv3_HY_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv3-HY-raw.png")
save(WMB_10Xv3_HY_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv3_HY_sc")
rm(WMB_10Xv3_HY_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv3-Isocortex-1-raw"], meta_data=meta_data) -> WMB_10Xv3_Isocortex_1_sc
DimPlot(WMB_10Xv3_Isocortex_1_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv3-Isocortex-1-raw.png")
save(WMB_10Xv3_Isocortex_1_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv3_Isocortex-1_sc")
rm(WMB_10Xv3_Isocortex_1_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv3-Isocortex-2-raw"], meta_data=meta_data) -> WMB_10Xv3_Isocortex_2_sc
DimPlot(WMB_10Xv3_Isocortex_2_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv3-Isocortex-2-raw.png")
save(WMB_10Xv3_Isocortex_2_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv3_Isocortex-2_sc")
rm(WMB_10Xv3_Isocortex_2_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv3-OLF-raw"], meta_data=meta_data) -> WMB_10Xv3_OLF_sc
DimPlot(WMB_10Xv3_OLF_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv3-OLF-raw.png")
save(WMB_10Xv3_OLF_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv3_OLF_sc")
rm(WMB_10Xv3_OLF_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv3-PAL-raw"], meta_data=meta_data) -> WMB_10Xv3_PAL_sc
DimPlot(WMB_10Xv3_PAL_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv3-PAL-raw.png")
save(WMB_10Xv3_PAL_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv3_PAL_sc")
rm(WMB_10Xv3_PAL_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv3-STR-raw"], meta_data=meta_data) -> WMB_10Xv3_STR_sc
DimPlot(WMB_10Xv3_STR_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv3-STR-raw.png")
save(WMB_10Xv3_STR_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv3_STR_sc")
rm(WMB_10Xv3_STR_sc);gc()

make_seurat_obj(raw_counts=rds_["WMB-10Xv3-TH-raw"], meta_data=meta_data) -> WMB_10Xv3_TH_sc
DimPlot(WMB_10Xv3_TH_sc, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/WMB-10Xv3-TH-raw.png")
save(WMB_10Xv3_TH_sc, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/sc/WMB_10Xv3_TH_sc")
rm(WMB_10Xv3_TH_sc);gc()
rm(list=ls())

.libPaths("/scratch/user/tongtong1993/project/R/4.4")
library(Seurat)
library(tidyverse)
library(future)
library(glmGamPoi)
options(future.globals.maxSize = Inf)
plan(sequential)


meta_data = read_csv("~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/metadata/WMB-10X/20241115/views/cell_metadata_with_cluster_annotation.csv")
meta_data %>% as.data.frame() %>% "rownames<-"(., .$cell_label) -> meta_data

class_ = meta_data$class %>% unique()
class_list = list(
		  Glutamatergic = class_[grepl("Glut", class_)],
		  GABAergic = class_[grepl("GABA", class_)],
		  Other_neuronal = class_[grepl("21|22", class_)],
		  Glial = class_[grepl("30|31|32", class_)],
		  Vascular = class_[grepl("33", class_)],
		  Immune = class_[grepl("34", class_)]
		  )

list.files("~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/WMB-10Xv3/20230630", pattern="\\.rds", full=T, recursive=F) -> rds_
rds_ %>% {.[grepl("Isocortex|OLF|STR|PAL|TH|HY", .)]} -> rds_
gsub("/.+/(.+)\\.rds", "\\1", rds_) -> names(rds_)

rds_ %>% 
	purrr::map(function(rds){
			   print(rds)
			   mat_ = readRDS(rds)
			   print(dim(mat_))
			   return(mat_)
		  }) -> rds_list

rds_list %>% 
	purrr::map(function(mat_){
			   rownames(mat_) -> genes_
		  }) %>% purrr::reduce(., intersect) -> intersected_genes_

#rds_list %>% 
#	purrr::map(function(mat_){
#			   colnames(mat_) -> cells_
#		  }) %>% unlist() -> cells_

rds_list %>% 
	purrr::map(function(mat_){
			   mat_[intersected_genes_, ] -> mat_
		  }) %>% do.call(cbind, .) -> merged_mat_
dim(merged_mat_)
rm(rds_list);gc()

seurat_obj <- CreateSeuratObject(counts = merged_mat_)
rm(merged_mat_);gc()

all(colnames(seurat_obj) %in% rownames(meta_data))
cbind(seurat_obj@meta.data, meta_data[colnames(seurat_obj), c("feature_matrix_label", "class")]) -> meta_data_1

meta_data_1 %>%
	dplyr::mutate(custom_class = case_when(
					       class %in% class_list$Glutamatergic ~ "Glutamatergic",
					       class %in% class_list$GABAergic ~ "GABAergic",
					       class %in% class_list$Other_neuronal ~ "Other_neuronal",
					       class %in% class_list$Glial ~ "Glial",
					       class %in% class_list$Vascular ~ "Vascular",
					       class %in% class_list$Immune ~ "Immune",
					       TRUE ~ class
					       )) -> meta_data_1 
seurat_obj@meta.data <- meta_data_1

cells_keep = seurat_obj@meta.data %>% dplyr::filter(custom_class != "Other_neuronal") %>% rownames()
seurat_obj = subset(seurat_obj, cells=cells_keep)

seurat_obj <- SCTransform(
			  seurat_obj,
			  assay = "RNA",
			  vst.flavor = "v2",
			  method = "glmGamPoi",
			  do.correct.umi = FALSE,
			  return.only.var.genes = TRUE,     # keeps memory down
			  residual.features = NULL,         # residuals only for var genes
			  conserve.memory = TRUE,
			  verbose = TRUE,
			  bin_size = 50
			 )
DefaultAssay(seurat_obj) <- "SCT"
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
save(seurat_obj, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/wmb_merged_ref_sc.rda")
rm(list=ls())

.libPaths("/scratch/user/tongtong1993/project/R/4.4")
library(Seurat)
library(future)

load("~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/wmb_merged_ref_sc.rda")
plan(multicore, workers = 12)
pca_mat <- Embeddings(seurat_obj, "pca")[, 1:30]
save(pca_mat, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/wmb_merged_ref_sc_pca.rda")
rm(list=ls())


.libPaths("/scratch/user/tongtong1993/project/R/4.4")
library(uwot)
library(future)
#library(RSpectra) #Don't attach this pkg if you want to use irlba

set.seed(42)
load("~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/wmb_merged_ref_sc_pca.rda")

umap_embedding <- uwot::umap(
			     X = pca_mat,
			     n_neighbors  = 30,
			     n_components = 2,
			     metric       = "cosine",
			     n_threads    = 12,
			     verbose      = TRUE,
			     init = "spectral" 
			     )
save(umap_embedding, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/wmb_merged_ref_sc_umap_irlba.rda")

umap_embedding_2 <- uwot::umap2(
			     X = pca_mat,
			     n_neighbors  = 30,
			     n_components = 2,
			     metric       = "cosine",
			     n_threads    = 12,
			     verbose      = TRUE,
			     init = "pca" 
			     )
save(umap_embedding_2, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/wmb_merged_ref_sc_umap2_pca.rda")

library(Seurat)
load("~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/wmb_merged_ref_sc.rda")
seurat_obj[["umap"]] <- CreateDimReducObject(
					     embeddings = umap_embedding_2,
					     key        = "UMAP_",
					     assay      = DefaultAssay(seurat_obj)
					     )

DimPlot(seurat_obj, reduction="umap", group.by="custom_class", label=T) -> p_umap
ggplot2::ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/wmb_umap.png")
ggplot2::ggsave(plot=p_umap, filename="~/Desktop/Research/Analysis/spatial_transcriptome/figures/wmb_umap.pdf") #Figure S2A

save(seurat_obj, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/mapmycells/raw_data/expression_matrices/wmb_merged_ref_sc.rda")
