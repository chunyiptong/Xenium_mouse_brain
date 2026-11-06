#ml GCC/13.2.0 OpenMPI/4.1.6 HDF5/1.14.3 R/4.4.1 UDUNITS/2.2.28 GDAL/3.9.0

.libPaths("/scratch/user/tongtong1993/project/R/4.4")

library(arrow); library(dplyr); library(Matrix); library(Seurat); library(sf)

# CT-8-1 ----
NUC_PATH = "~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058138__CT-8-1__20250604__231706/nucleus_boundaries.parquet"

df <- read_parquet(NUC_PATH)

build_polygon <- function(g) {
	  # keep order as given; drop consecutive duplicate vertices if any
	  g <- g |> distinct(vertex_x, vertex_y, .keep_all = TRUE)
  ring <- as.matrix(g[, c("vertex_x","vertex_y")])
    if (!all(ring[1,] == ring[nrow(ring),])) ring <- rbind(ring, ring[1,])  # close ring
    sf::st_polygon(list(ring))
    
}

g <- df |> group_by(cell_id, label_id)
keys   <- dplyr::group_keys(g)                         # IDs for each nucleus
polys  <- group_split(g) |> lapply(build_polygon)
sfc    <- sf::st_sfc(polys)                      # microns; no CRS
cent   <- sf::st_centroid(sfc)
xy     <- sf::st_coordinates(cent)

centroids <- bind_cols(keys, tibble(nucleus_x = xy[,1], nucleus_y = xy[,2]))
centroids[["nucleus_poly"]] <- polys
cell_keep = centroids$cell_id %>% table() %>% {.[.==1]} %>% names()
length(cell_keep) 
centroids %>% dplyr::filter(cell_id %in% cell_keep) %>% dplyr::select(-label_id) -> centroids
nucleus_centroids = centroids
dim(nucleus_centroids)

rm(cent, centroids, df, g, keys, sfc, xy, cell_keep, NUC_PATH, polys);gc()

cb = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058138__CT-8-1__20250604__231706/proseg/outs/cell_boundaries.parquet", col_select = c("cell_id", "vertex_x", "vertex_y"))

cells <- cb %>%
	  group_by(cell_id) %>%
	  summarise(
		    geometry = list({
			          r <- as.matrix(cbind(vertex_x, vertex_y))
				        if (!all(r[1, ] == r[nrow(r), ])) r <- rbind(r, r[1, ])  # close ring
				        sf::st_polygon(list(r))   # sfg
					    
		    }),
		        .groups = "drop"
			  
		    ) %>%
	    sf::st_as_sf(sf_column_name = "geometry")

pts   <- sf::st_as_sf(nucleus_centroids, coords = c("nucleus_x","nucleus_y"), remove = FALSE)
match <- sf::st_join(pts, cells["cell_id"], join = sf::st_within, left = FALSE)
match %>% as_tibble() -> match
colnames(match) = c("nucleus_id", "nucleus_x", "nucleus_y", "nucleus_poly", "nucleus_geometry", "cell_id")

match$nucleus_id %>% table() %>% .[.==1] %>% names() -> unique_uncleus
match %>% dplyr::filter(nucleus_id %in% unique_uncleus) -> match
match$cell_id %>% table() %>% .[.==1] %>% names() -> unique_cell_id
match %>% dplyr::filter(cell_id %in% unique_cell_id) -> match
nrow(match)

rm(cb, cells, nucleus_centroids, pts, unique_cell_id, unique_uncleus); gc()

tp = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058138__CT-8-1__20250604__231706/proseg/outs/transcripts.parquet")

tp %>% dplyr::filter(cell_id != "UNASSIGNED") -> tp
tp %>% dplyr::filter(cell_id %in% match$cell_id) -> tp
tp %>% dplyr::filter(qv >= 20) -> tp

all(tp$cell_id %in% match$cell_id)
gc()

sfc <- sf::st_sfc(match[["nucleus_poly"]])
nucleus_poly = sf::st_as_sf(sfc)
rm(sfc);gc()
tp %>% dplyr::select(cell_id, transcript_id, feature_name, x_location, y_location) -> transcripts
transcripts <- sf::st_as_sf(transcripts, coords = c("x_location", "y_location"), crs = NA)
nucleus_poly[["nucleus_id"]] = match$nucleus_id
nucleus_poly[["cell_id"]] = match$cell_id
rm(match, tp);gc()
transcript_in_nucleus <- sf::st_join(transcripts, nucleus_poly, join = sf::st_within, left = FALSE)
gc()
rm(nucleus_poly, transcripts);gc()

transcript_in_nucleus %>% dplyr::filter(cell_id.x == cell_id.y) -> transcript_in_nucleus
transcript_in_nucleus %>% as_tibble() %>% .[, c(2,3,5,6)] %>% "colnames<-"(., c("transcript_id", "feature_name", "nucleus_id", "cell_id")) -> transcript_in_nucleus
transcript_in_nucleus %>% count(cell_id, feature_name, name = "n") -> transcript_in_nucleus

cells <- unique(transcript_in_nucleus$cell_id); genes <- unique(transcript_in_nucleus$feature_name)

i <- match(transcript_in_nucleus$cell_id,  cells)
j <- match(transcript_in_nucleus$feature_name, genes)
mat <- sparseMatrix(i = i, j = j, x = transcript_in_nucleus$n,
                    dims = c(length(cells), length(genes)),
                    dimnames = list(cells, genes)) # cells x genes
mat = t(mat)
save(mat, file = "~/Desktop/Research/Analysis/spatial_transcriptome/rda/CT-8-1_RNA_counts_nucleus_mat.rda")
rm(list=ls())

#CT-8-2 ----
.libPaths("/scratch/user/tongtong1993/project/R/4.4")

library(arrow); library(dplyr); library(Matrix); library(Seurat); library(sf)

NUC_PATH = "~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058138__CT-8-2__20250604__231707/nucleus_boundaries.parquet"

df <- read_parquet(NUC_PATH)

build_polygon <- function(g) {
		  # keep order as given; drop consecutive duplicate vertices if any
		  g <- g |> distinct(vertex_x, vertex_y, .keep_all = TRUE)
  ring <- as.matrix(g[, c("vertex_x","vertex_y")])
      if (!all(ring[1,] == ring[nrow(ring),])) ring <- rbind(ring, ring[1,])  # close ring
      sf::st_polygon(list(ring))
          
      
}

g <- df |> group_by(cell_id, label_id)
keys   <- dplyr::group_keys(g)                         # IDs for each nucleus
polys  <- group_split(g) |> lapply(build_polygon)
sfc    <- sf::st_sfc(polys)                      # microns; no CRS
cent   <- sf::st_centroid(sfc)
xy     <- sf::st_coordinates(cent)

centroids <- bind_cols(keys, tibble(nucleus_x = xy[,1], nucleus_y = xy[,2]))
centroids[["nucleus_poly"]] <- polys
cell_keep = centroids$cell_id %>% table() %>% {.[.==1]} %>% names()
length(cell_keep) 
centroids %>% dplyr::filter(cell_id %in% cell_keep) %>% dplyr::select(-label_id) -> centroids
nucleus_centroids = centroids
dim(nucleus_centroids)

rm(cent, centroids, df, g, keys, sfc, xy, cell_keep, NUC_PATH, polys);gc()

cb = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058138__CT-8-2__20250604__231707/proseg/outs/cell_boundaries.parquet", col_select=c("cell_id", "vertex_x", "vertex_y"))

cells <- cb %>%
	  group_by(cell_id) %>%
		  summarise(
			    geometry = list({
				  r <- as.matrix(cbind(vertex_x, vertex_y))
				if (!all(r[1, ] == r[nrow(r), ])) r <- rbind(r, r[1, ])  # close ring
				sf::st_polygon(list(r))   # sfg
					    }), .groups = "drop"
			    ) %>%
			    sf::st_as_sf(sf_column_name = "geometry")

pts   <- sf::st_as_sf(nucleus_centroids, coords = c("nucleus_x","nucleus_y"), remove = FALSE)
match <- sf::st_join(pts, cells["cell_id"], join = sf::st_within, left = FALSE)
match %>% as_tibble() -> match
colnames(match) = c("nucleus_id", "nucleus_x", "nucleus_y", "nucleus_poly", "nucleus_geometry", "cell_id")

match$nucleus_id %>% table() %>% .[.==1] %>% names() -> unique_uncleus
match %>% dplyr::filter(nucleus_id %in% unique_uncleus) -> match
match$cell_id %>% table() %>% .[.==1] %>% names() -> unique_cell_id
match %>% dplyr::filter(cell_id %in% unique_cell_id) -> match
nrow(match)

rm(cb, cells, nucleus_centroids, pts, unique_cell_id, unique_uncleus); gc()

tp = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058138__CT-8-2__20250604__231707/proseg/outs/transcripts.parquet")

tp %>% dplyr::filter(cell_id != "UNASSIGNED") -> tp
tp %>% dplyr::filter(cell_id %in% match$cell_id) -> tp
tp %>% dplyr::filter(qv >= 20) -> tp

all(tp$cell_id %in% match$cell_id)
gc()

sfc <- sf::st_sfc(match[["nucleus_poly"]])
nucleus_poly = sf::st_as_sf(sfc)
rm(sfc);gc()
tp %>% dplyr::select(cell_id, transcript_id, feature_name, x_location, y_location) -> transcripts
transcripts <- sf::st_as_sf(transcripts, coords = c("x_location", "y_location"), crs = NA)
nucleus_poly[["nucleus_id"]] = match$nucleus_id
nucleus_poly[["cell_id"]] = match$cell_id
rm(match, tp);gc()
transcript_in_nucleus <- sf::st_join(transcripts, nucleus_poly, join = sf::st_within, left = FALSE)
gc()
rm(nucleus_poly, transcripts);gc()

transcript_in_nucleus %>% dplyr::filter(cell_id.x == cell_id.y) -> transcript_in_nucleus
transcript_in_nucleus %>% as_tibble() %>% .[, c(2,3,5,6)] %>% "colnames<-"(., c("transcript_id", "feature_name", "nucleus_id", "cell_id")) -> transcript_in_nucleus
transcript_in_nucleus %>% count(cell_id, feature_name, name = "n") -> transcript_in_nucleus

cells <- unique(transcript_in_nucleus$cell_id); genes <- unique(transcript_in_nucleus$feature_name)

i <- match(transcript_in_nucleus$cell_id,  cells)
j <- match(transcript_in_nucleus$feature_name, genes)
mat <- sparseMatrix(i = i, j = j, x = transcript_in_nucleus$n,
		                        dims = c(length(cells), length(genes)),
					    dimnames = list(cells, genes)) # cells x genes
mat = t(mat)
save(mat, file = "~/Desktop/Research/Analysis/spatial_transcriptome/rda/CT-8-2_RNA_counts_nucleus_mat.rda")
rm(list=ls())

#CT-8-3 ----
.libPaths("/scratch/user/tongtong1993/project/R/4.4")

library(arrow); library(dplyr); library(Matrix); library(Seurat); library(sf)

NUC_PATH = "~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058138__CT-8-3__20250604__231707/nucleus_boundaries.parquet"

df <- read_parquet(NUC_PATH)

build_polygon <- function(g) {
		  # keep order as given; drop consecutive duplicate vertices if any
		  g <- g |> distinct(vertex_x, vertex_y, .keep_all = TRUE)
  ring <- as.matrix(g[, c("vertex_x","vertex_y")])
      if (!all(ring[1,] == ring[nrow(ring),])) ring <- rbind(ring, ring[1,])  # close ring
      sf::st_polygon(list(ring))
          
      
}

g <- df |> group_by(cell_id, label_id)
keys   <- dplyr::group_keys(g)                         # IDs for each nucleus
polys  <- group_split(g) |> lapply(build_polygon)
sfc    <- sf::st_sfc(polys)                      # microns; no CRS
cent   <- sf::st_centroid(sfc)
xy     <- sf::st_coordinates(cent)

centroids <- bind_cols(keys, tibble(nucleus_x = xy[,1], nucleus_y = xy[,2]))
centroids[["nucleus_poly"]] <- polys
cell_keep = centroids$cell_id %>% table() %>% {.[.==1]} %>% names()
length(cell_keep) 
centroids %>% dplyr::filter(cell_id %in% cell_keep) %>% dplyr::select(-label_id) -> centroids
nucleus_centroids = centroids
dim(nucleus_centroids)

rm(cent, centroids, df, g, keys, sfc, xy, cell_keep, NUC_PATH, polys);gc()

cb = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058138__CT-8-3__20250604__231707/proseg/outs/cell_boundaries.parquet", col_select=c("cell_id", "vertex_x", "vertex_y"))

cells <- cb %>%
	  group_by(cell_id) %>%
		  summarise(
			    geometry = list({
				  r <- as.matrix(cbind(vertex_x, vertex_y))
				if (!all(r[1, ] == r[nrow(r), ])) r <- rbind(r, r[1, ])  # close ring
				sf::st_polygon(list(r))   # sfg
					    }), .groups = "drop"
			    ) %>%
			    sf::st_as_sf(sf_column_name = "geometry")

pts   <- sf::st_as_sf(nucleus_centroids, coords = c("nucleus_x","nucleus_y"), remove = FALSE)
match <- sf::st_join(pts, cells["cell_id"], join = sf::st_within, left = FALSE)
match %>% as_tibble() -> match
colnames(match) = c("nucleus_id", "nucleus_x", "nucleus_y", "nucleus_poly", "nucleus_geometry", "cell_id")

match$nucleus_id %>% table() %>% .[.==1] %>% names() -> unique_uncleus
match %>% dplyr::filter(nucleus_id %in% unique_uncleus) -> match
match$cell_id %>% table() %>% .[.==1] %>% names() -> unique_cell_id
match %>% dplyr::filter(cell_id %in% unique_cell_id) -> match
nrow(match)

rm(cb, cells, nucleus_centroids, pts, unique_cell_id, unique_uncleus); gc()

tp = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058138__CT-8-3__20250604__231707/proseg/outs/transcripts.parquet")

tp %>% dplyr::filter(cell_id != "UNASSIGNED") -> tp
tp %>% dplyr::filter(cell_id %in% match$cell_id) -> tp
tp %>% dplyr::filter(qv >= 20) -> tp

all(tp$cell_id %in% match$cell_id)
gc()

sfc <- sf::st_sfc(match[["nucleus_poly"]])
nucleus_poly = sf::st_as_sf(sfc)
rm(sfc);gc()
tp %>% dplyr::select(cell_id, transcript_id, feature_name, x_location, y_location) -> transcripts
transcripts <- sf::st_as_sf(transcripts, coords = c("x_location", "y_location"), crs = NA)
nucleus_poly[["nucleus_id"]] = match$nucleus_id
nucleus_poly[["cell_id"]] = match$cell_id
rm(match, tp);gc()
transcript_in_nucleus <- sf::st_join(transcripts, nucleus_poly, join = sf::st_within, left = FALSE)
gc()
rm(nucleus_poly, transcripts);gc()

transcript_in_nucleus %>% dplyr::filter(cell_id.x == cell_id.y) -> transcript_in_nucleus
transcript_in_nucleus %>% as_tibble() %>% .[, c(2,3,5,6)] %>% "colnames<-"(., c("transcript_id", "feature_name", "nucleus_id", "cell_id")) -> transcript_in_nucleus
transcript_in_nucleus %>% count(cell_id, feature_name, name = "n") -> transcript_in_nucleus

cells <- unique(transcript_in_nucleus$cell_id); genes <- unique(transcript_in_nucleus$feature_name)

i <- match(transcript_in_nucleus$cell_id,  cells)
j <- match(transcript_in_nucleus$feature_name, genes)
mat <- sparseMatrix(i = i, j = j, x = transcript_in_nucleus$n,
		                        dims = c(length(cells), length(genes)),
					    dimnames = list(cells, genes)) # cells x genes
mat = t(mat)
save(mat, file = "~/Desktop/Research/Analysis/spatial_transcriptome/rda/CT-8-3_RNA_counts_nucleus_mat.rda")
rm(list=ls())

#CT-20-1 ----
.libPaths("/scratch/user/tongtong1993/project/R/4.4")

library(arrow); library(dplyr); library(Matrix); library(Seurat); library(sf)

NUC_PATH = "~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058140__CT-20-1__20250604__231707/nucleus_boundaries.parquet"

df <- read_parquet(NUC_PATH)

build_polygon <- function(g) {
		  # keep order as given; drop consecutive duplicate vertices if any
		  g <- g |> distinct(vertex_x, vertex_y, .keep_all = TRUE)
  ring <- as.matrix(g[, c("vertex_x","vertex_y")])
      if (!all(ring[1,] == ring[nrow(ring),])) ring <- rbind(ring, ring[1,])  # close ring
      sf::st_polygon(list(ring))
          
      
}

g <- df |> group_by(cell_id, label_id)
keys   <- dplyr::group_keys(g)                         # IDs for each nucleus
polys  <- group_split(g) |> lapply(build_polygon)
sfc    <- sf::st_sfc(polys)                      # microns; no CRS
cent   <- sf::st_centroid(sfc)
xy     <- sf::st_coordinates(cent)

centroids <- bind_cols(keys, tibble(nucleus_x = xy[,1], nucleus_y = xy[,2]))
centroids[["nucleus_poly"]] <- polys
cell_keep = centroids$cell_id %>% table() %>% {.[.==1]} %>% names()
length(cell_keep) 
centroids %>% dplyr::filter(cell_id %in% cell_keep) %>% dplyr::select(-label_id) -> centroids
nucleus_centroids = centroids
dim(nucleus_centroids)

rm(cent, centroids, df, g, keys, sfc, xy, cell_keep, NUC_PATH, polys);gc()

cb = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058140__CT-20-1__20250604__231707/proseg/outs/cell_boundaries.parquet", col_select=c("cell_id", "vertex_x", "vertex_y"))

cells <- cb %>%
	  group_by(cell_id) %>%
		  summarise(
			    geometry = list({
				  r <- as.matrix(cbind(vertex_x, vertex_y))
				if (!all(r[1, ] == r[nrow(r), ])) r <- rbind(r, r[1, ])  # close ring
				sf::st_polygon(list(r))   # sfg
					    }), .groups = "drop"
			    ) %>%
			    sf::st_as_sf(sf_column_name = "geometry")

pts   <- sf::st_as_sf(nucleus_centroids, coords = c("nucleus_x","nucleus_y"), remove = FALSE)
match <- sf::st_join(pts, cells["cell_id"], join = sf::st_within, left = FALSE)
match %>% as_tibble() -> match
colnames(match) = c("nucleus_id", "nucleus_x", "nucleus_y", "nucleus_poly", "nucleus_geometry", "cell_id")

match$nucleus_id %>% table() %>% .[.==1] %>% names() -> unique_uncleus
match %>% dplyr::filter(nucleus_id %in% unique_uncleus) -> match
match$cell_id %>% table() %>% .[.==1] %>% names() -> unique_cell_id
match %>% dplyr::filter(cell_id %in% unique_cell_id) -> match
nrow(match)

rm(cb, cells, nucleus_centroids, pts, unique_cell_id, unique_uncleus); gc()

tp = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058140__CT-20-1__20250604__231707/proseg/outs/transcripts.parquet")

tp %>% dplyr::filter(cell_id != "UNASSIGNED") -> tp
tp %>% dplyr::filter(cell_id %in% match$cell_id) -> tp
tp %>% dplyr::filter(qv >= 20) -> tp

all(tp$cell_id %in% match$cell_id)
gc()

sfc <- sf::st_sfc(match[["nucleus_poly"]])
nucleus_poly = sf::st_as_sf(sfc)
rm(sfc);gc()
tp %>% dplyr::select(cell_id, transcript_id, feature_name, x_location, y_location) -> transcripts
transcripts <- sf::st_as_sf(transcripts, coords = c("x_location", "y_location"), crs = NA)
nucleus_poly[["nucleus_id"]] = match$nucleus_id
nucleus_poly[["cell_id"]] = match$cell_id
rm(match, tp);gc()
transcript_in_nucleus <- sf::st_join(transcripts, nucleus_poly, join = sf::st_within, left = FALSE)
gc()
rm(nucleus_poly, transcripts);gc()

transcript_in_nucleus %>% dplyr::filter(cell_id.x == cell_id.y) -> transcript_in_nucleus
transcript_in_nucleus %>% as_tibble() %>% .[, c(2,3,5,6)] %>% "colnames<-"(., c("transcript_id", "feature_name", "nucleus_id", "cell_id")) -> transcript_in_nucleus
transcript_in_nucleus %>% count(cell_id, feature_name, name = "n") -> transcript_in_nucleus

cells <- unique(transcript_in_nucleus$cell_id); genes <- unique(transcript_in_nucleus$feature_name)

i <- match(transcript_in_nucleus$cell_id,  cells)
j <- match(transcript_in_nucleus$feature_name, genes)
mat <- sparseMatrix(i = i, j = j, x = transcript_in_nucleus$n,
		                        dims = c(length(cells), length(genes)),
					    dimnames = list(cells, genes)) # cells x genes
mat = t(mat)
save(mat, file = "~/Desktop/Research/Analysis/spatial_transcriptome/rda/CT-20-1_RNA_counts_nucleus_mat.rda")
rm(list=ls())

#CT-20-2 ----
.libPaths("/scratch/user/tongtong1993/project/R/4.4")

library(arrow); library(dplyr); library(Matrix); library(Seurat); library(sf)

NUC_PATH = "~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058140__CT-20-2__20250604__231707/nucleus_boundaries.parquet"

df <- read_parquet(NUC_PATH)

build_polygon <- function(g) {
		  # keep order as given; drop consecutive duplicate vertices if any
		  g <- g |> distinct(vertex_x, vertex_y, .keep_all = TRUE)
  ring <- as.matrix(g[, c("vertex_x","vertex_y")])
      if (!all(ring[1,] == ring[nrow(ring),])) ring <- rbind(ring, ring[1,])  # close ring
      sf::st_polygon(list(ring))
          
      
}

g <- df |> group_by(cell_id, label_id)
keys   <- dplyr::group_keys(g)                         # IDs for each nucleus
polys  <- group_split(g) |> lapply(build_polygon)
sfc    <- sf::st_sfc(polys)                      # microns; no CRS
cent   <- sf::st_centroid(sfc)
xy     <- sf::st_coordinates(cent)

centroids <- bind_cols(keys, tibble(nucleus_x = xy[,1], nucleus_y = xy[,2]))
centroids[["nucleus_poly"]] <- polys
cell_keep = centroids$cell_id %>% table() %>% {.[.==1]} %>% names()
length(cell_keep) 
centroids %>% dplyr::filter(cell_id %in% cell_keep) %>% dplyr::select(-label_id) -> centroids
nucleus_centroids = centroids
dim(nucleus_centroids)

rm(cent, centroids, df, g, keys, sfc, xy, cell_keep, NUC_PATH, polys);gc()

cb = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058140__CT-20-2__20250604__231707/proseg/outs/cell_boundaries.parquet", col_select=c("cell_id", "vertex_x", "vertex_y"))

cells <- cb %>%
	  group_by(cell_id) %>%
		  summarise(
			    geometry = list({
				  r <- as.matrix(cbind(vertex_x, vertex_y))
				if (!all(r[1, ] == r[nrow(r), ])) r <- rbind(r, r[1, ])  # close ring
				sf::st_polygon(list(r))   # sfg
					    }), .groups = "drop"
			    ) %>%
			    sf::st_as_sf(sf_column_name = "geometry")

pts   <- sf::st_as_sf(nucleus_centroids, coords = c("nucleus_x","nucleus_y"), remove = FALSE)
match <- sf::st_join(pts, cells["cell_id"], join = sf::st_within, left = FALSE)
match %>% as_tibble() -> match
colnames(match) = c("nucleus_id", "nucleus_x", "nucleus_y", "nucleus_poly", "nucleus_geometry", "cell_id")

match$nucleus_id %>% table() %>% .[.==1] %>% names() -> unique_uncleus
match %>% dplyr::filter(nucleus_id %in% unique_uncleus) -> match
match$cell_id %>% table() %>% .[.==1] %>% names() -> unique_cell_id
match %>% dplyr::filter(cell_id %in% unique_cell_id) -> match
nrow(match)

rm(cb, cells, nucleus_centroids, pts, unique_cell_id, unique_uncleus); gc()

tp = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058140__CT-20-2__20250604__231707/proseg/outs/transcripts.parquet")

tp %>% dplyr::filter(cell_id != "UNASSIGNED") -> tp
tp %>% dplyr::filter(cell_id %in% match$cell_id) -> tp
tp %>% dplyr::filter(qv >= 20) -> tp

all(tp$cell_id %in% match$cell_id)
gc()

sfc <- sf::st_sfc(match[["nucleus_poly"]])
nucleus_poly = sf::st_as_sf(sfc)
rm(sfc);gc()
tp %>% dplyr::select(cell_id, transcript_id, feature_name, x_location, y_location) -> transcripts
transcripts <- sf::st_as_sf(transcripts, coords = c("x_location", "y_location"), crs = NA)
nucleus_poly[["nucleus_id"]] = match$nucleus_id
nucleus_poly[["cell_id"]] = match$cell_id
rm(match, tp);gc()
transcript_in_nucleus <- sf::st_join(transcripts, nucleus_poly, join = sf::st_within, left = FALSE)
gc()
rm(nucleus_poly, transcripts);gc()

transcript_in_nucleus %>% dplyr::filter(cell_id.x == cell_id.y) -> transcript_in_nucleus
transcript_in_nucleus %>% as_tibble() %>% .[, c(2,3,5,6)] %>% "colnames<-"(., c("transcript_id", "feature_name", "nucleus_id", "cell_id")) -> transcript_in_nucleus
transcript_in_nucleus %>% count(cell_id, feature_name, name = "n") -> transcript_in_nucleus

cells <- unique(transcript_in_nucleus$cell_id); genes <- unique(transcript_in_nucleus$feature_name)

i <- match(transcript_in_nucleus$cell_id,  cells)
j <- match(transcript_in_nucleus$feature_name, genes)
mat <- sparseMatrix(i = i, j = j, x = transcript_in_nucleus$n,
		                        dims = c(length(cells), length(genes)),
					    dimnames = list(cells, genes)) # cells x genes
mat = t(mat)
save(mat, file = "~/Desktop/Research/Analysis/spatial_transcriptome/rda/CT-20-2_RNA_counts_nucleus_mat.rda")
rm(list=ls())

#CT-20-3 ----
.libPaths("/scratch/user/tongtong1993/project/R/4.4")

library(arrow); library(dplyr); library(Matrix); library(Seurat); library(sf)

NUC_PATH = "~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058140__CT-20-3__20250604__231707/nucleus_boundaries.parquet"

df <- read_parquet(NUC_PATH)

build_polygon <- function(g) {
		  # keep order as given; drop consecutive duplicate vertices if any
		  g <- g |> distinct(vertex_x, vertex_y, .keep_all = TRUE)
  ring <- as.matrix(g[, c("vertex_x","vertex_y")])
      if (!all(ring[1,] == ring[nrow(ring),])) ring <- rbind(ring, ring[1,])  # close ring
      sf::st_polygon(list(ring))
          
      
}

g <- df |> group_by(cell_id, label_id)
keys   <- dplyr::group_keys(g)                         # IDs for each nucleus
polys  <- group_split(g) |> lapply(build_polygon)
sfc    <- sf::st_sfc(polys)                      # microns; no CRS
cent   <- sf::st_centroid(sfc)
xy     <- sf::st_coordinates(cent)

centroids <- bind_cols(keys, tibble(nucleus_x = xy[,1], nucleus_y = xy[,2]))
centroids[["nucleus_poly"]] <- polys
cell_keep = centroids$cell_id %>% table() %>% {.[.==1]} %>% names()
length(cell_keep) 
centroids %>% dplyr::filter(cell_id %in% cell_keep) %>% dplyr::select(-label_id) -> centroids
nucleus_centroids = centroids
dim(nucleus_centroids)

rm(cent, centroids, df, g, keys, sfc, xy, cell_keep, NUC_PATH, polys);gc()

cb = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058140__CT-20-3__20250604__231707/proseg/outs/cell_boundaries.parquet", col_select=c("cell_id", "vertex_x", "vertex_y"))

cells <- cb %>%
	  group_by(cell_id) %>%
		  summarise(
			    geometry = list({
				  r <- as.matrix(cbind(vertex_x, vertex_y))
				if (!all(r[1, ] == r[nrow(r), ])) r <- rbind(r, r[1, ])  # close ring
				sf::st_polygon(list(r))   # sfg
					    }), .groups = "drop"
			    ) %>%
			    sf::st_as_sf(sf_column_name = "geometry")

pts   <- sf::st_as_sf(nucleus_centroids, coords = c("nucleus_x","nucleus_y"), remove = FALSE)
match <- sf::st_join(pts, cells["cell_id"], join = sf::st_within, left = FALSE)
match %>% as_tibble() -> match
colnames(match) = c("nucleus_id", "nucleus_x", "nucleus_y", "nucleus_poly", "nucleus_geometry", "cell_id")

match$nucleus_id %>% table() %>% .[.==1] %>% names() -> unique_uncleus
match %>% dplyr::filter(nucleus_id %in% unique_uncleus) -> match
match$cell_id %>% table() %>% .[.==1] %>% names() -> unique_cell_id
match %>% dplyr::filter(cell_id %in% unique_cell_id) -> match
nrow(match)

rm(cb, cells, nucleus_centroids, pts, unique_cell_id, unique_uncleus); gc()

tp = read_parquet("~/Desktop/Research/Analysis/spatial_transcriptome/data/2025-06-11_Two_time_points_data/20250604__231552__UTSWMC-Yuchao-MouseBrain/output-XETG00124R1__0058140__CT-20-3__20250604__231707/proseg/outs/transcripts.parquet")

tp %>% dplyr::filter(cell_id != "UNASSIGNED") -> tp
tp %>% dplyr::filter(cell_id %in% match$cell_id) -> tp
tp %>% dplyr::filter(qv >= 20) -> tp

all(tp$cell_id %in% match$cell_id)
gc()

sfc <- sf::st_sfc(match[["nucleus_poly"]])
nucleus_poly = sf::st_as_sf(sfc)
rm(sfc);gc()
tp %>% dplyr::select(cell_id, transcript_id, feature_name, x_location, y_location) -> transcripts
transcripts <- sf::st_as_sf(transcripts, coords = c("x_location", "y_location"), crs = NA)
nucleus_poly[["nucleus_id"]] = match$nucleus_id
nucleus_poly[["cell_id"]] = match$cell_id
rm(match, tp);gc()
transcript_in_nucleus <- sf::st_join(transcripts, nucleus_poly, join = sf::st_within, left = FALSE)
gc()
rm(nucleus_poly, transcripts);gc()

transcript_in_nucleus %>% dplyr::filter(cell_id.x == cell_id.y) -> transcript_in_nucleus
transcript_in_nucleus %>% as_tibble() %>% .[, c(2,3,5,6)] %>% "colnames<-"(., c("transcript_id", "feature_name", "nucleus_id", "cell_id")) -> transcript_in_nucleus
transcript_in_nucleus %>% count(cell_id, feature_name, name = "n") -> transcript_in_nucleus

cells <- unique(transcript_in_nucleus$cell_id); genes <- unique(transcript_in_nucleus$feature_name)

i <- match(transcript_in_nucleus$cell_id,  cells)
j <- match(transcript_in_nucleus$feature_name, genes)
mat <- sparseMatrix(i = i, j = j, x = transcript_in_nucleus$n,
		                        dims = c(length(cells), length(genes)),
					    dimnames = list(cells, genes)) # cells x genes
mat = t(mat)
save(mat, file = "~/Desktop/Research/Analysis/spatial_transcriptome/rda/CT-20-3_RNA_counts_nucleus_mat.rda")
rm(list=ls())

.libPaths("/scratch/user/tongtong1993/project/R/4.4")
library(tidyverse)
library(Seurat)

nucleus_RNA_count_list = list.files(path = "~/Desktop/Research/Analysis/spatial_transcriptome/rda", pattern="RNA_counts_nucleus_mat\\.rda", full=T) %>% gtools::mixedsort() %>% rev()
c("CT-8-1", "CT-8-2", "CT-8-3", "CT-20-1", "CT-20-2", "CT-20-3") -> names(nucleus_RNA_count_list)

load(nucleus_RNA_count_list["CT-8-1"])
CT_8_1_nuc_RNA_count = mat;rm(mat);gc()
load(nucleus_RNA_count_list["CT-8-2"])
CT_8_2_nuc_RNA_count = mat;rm(mat);gc()
load(nucleus_RNA_count_list["CT-8-3"])
CT_8_3_nuc_RNA_count = mat;rm(mat);gc()
load(nucleus_RNA_count_list["CT-20-1"])
CT_20_1_nuc_RNA_count = mat;rm(mat);gc()
load(nucleus_RNA_count_list["CT-20-2"])
CT_20_2_nuc_RNA_count = mat;rm(mat);gc()
load(nucleus_RNA_count_list["CT-20-3"])
CT_20_3_nuc_RNA_count = mat;rm(mat);gc()
rm(nucleus_RNA_count_list)

colnames(CT_8_1_nuc_RNA_count) <- sprintf("CT-8-1_%s", colnames(CT_8_1_nuc_RNA_count))
colnames(CT_8_2_nuc_RNA_count) <- sprintf("CT-8-2_%s", colnames(CT_8_2_nuc_RNA_count))
colnames(CT_8_3_nuc_RNA_count) <- sprintf("CT-8-3_%s", colnames(CT_8_3_nuc_RNA_count))
colnames(CT_20_1_nuc_RNA_count) <- sprintf("CT-20-1_%s", colnames(CT_20_1_nuc_RNA_count))
colnames(CT_20_2_nuc_RNA_count) <- sprintf("CT-20-2_%s", colnames(CT_20_2_nuc_RNA_count))
colnames(CT_20_3_nuc_RNA_count) <- sprintf("CT-20-3_%s", colnames(CT_20_3_nuc_RNA_count))

load("~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/merged-sc-filtered-normalized-mini.rda")
xenium_obj_filtered = mini; rm(mini) ;gc()
rownames(xenium_obj_filtered) -> genes_

nrow(CT_8_1_nuc_RNA_count);nrow(CT_8_2_nuc_RNA_count);nrow(CT_8_3_nuc_RNA_count)
nrow(CT_20_1_nuc_RNA_count);nrow(CT_20_2_nuc_RNA_count);nrow(CT_20_3_nuc_RNA_count)

CT_8_1_nuc_RNA_count[genes_,] -> CT_8_1_nuc_RNA_count
CT_8_2_nuc_RNA_count[genes_,] -> CT_8_2_nuc_RNA_count
CT_8_3_nuc_RNA_count[genes_,] -> CT_8_3_nuc_RNA_count
CT_20_1_nuc_RNA_count[genes_,] -> CT_20_1_nuc_RNA_count
CT_20_2_nuc_RNA_count[genes_,] -> CT_20_2_nuc_RNA_count
CT_20_3_nuc_RNA_count[genes_,] -> CT_20_3_nuc_RNA_count

nuc_RNA_count = cbind(CT_8_1_nuc_RNA_count, CT_8_2_nuc_RNA_count, CT_8_3_nuc_RNA_count, CT_20_1_nuc_RNA_count, CT_20_2_nuc_RNA_count, CT_20_3_nuc_RNA_count)
rm(CT_8_1_nuc_RNA_count, CT_8_2_nuc_RNA_count, CT_8_3_nuc_RNA_count, CT_20_1_nuc_RNA_count, CT_20_2_nuc_RNA_count, CT_20_3_nuc_RNA_count);gc()

cells_matched_ = match(colnames(xenium_obj_filtered), colnames(nuc_RNA_count))
cells_matched_ = cells_matched_[!is.na(cells_matched_)]
nuc_RNA_count = nuc_RNA_count[, cells_matched_]
all(colnames(nuc_RNA_count) %in% colnames(xenium_obj_filtered))

xenium_obj_filtered = subset(xenium_obj_filtered, cells = colnames(nuc_RNA_count))
all(colnames(nuc_RNA_count) == colnames(xenium_obj_filtered))
all(rownames(nuc_RNA_count) == rownames(xenium_obj_filtered))
rm(cells_matched_, genes_)

cell_RNA_count = GetAssayData(object=xenium_obj_filtered, assay="Xenium", layer="counts")

nuc_RNA_count = t(nuc_RNA_count)
cell_RNA_count = t(cell_RNA_count)
nuc_cell_ratio = nuc_RNA_count

for (j in seq_len(ncol(nuc_RNA_count))) {
  print(j)
  idx  <- (nuc_RNA_count@p[j] + 1) : nuc_RNA_count@p[j + 1]   # positions in @x for column j
    if (length(idx) == 0) next
    rows <- nuc_RNA_count@i[idx] + 1
      den  <- cell_RNA_count[rows, j]
      # choose behavior for den==0:
      good <- den != 0
	print(all(good))
        nuc_cell_ratio@x[idx[good]] <- nuc_cell_ratio@x[idx[good]] / den[good]
	  nuc_cell_ratio@x[idx[!good]] <- 0                  # drop undefined divisions
}
nuc_cell_ratio <- Matrix::drop0(nuc_cell_ratio)

save(nuc_RNA_count, nuc_cell_ratio, file="~/Desktop/Research/Analysis/spatial_transcriptome/rda/nuc_cell_ratio.mat")
rm(den, good, idx, j, rows);gc()
rm(cell_RNA_count);gc()

nuc_RNA_count = t(nuc_RNA_count)
nuc_cell_ratio = t(nuc_cell_ratio)
dim(nuc_RNA_count);dim(nuc_cell_ratio)

nuc_RNA_assay <- CreateAssayObject(counts = nuc_RNA_count)
xenium_obj_filtered[["Xenium_nuc"]] <- nuc_RNA_assay
rm(nuc_RNA_count);gc()
rm(nuc_RNA_assay);gc()

nuc_cell_ratio_assay <- CreateAssayObject(data = nuc_cell_ratio)
xenium_obj_filtered[["Xenium_nuc_cell_ratio"]] <- nuc_cell_ratio_assay
rm(nuc_cell_ratio, nuc_cell_ratio_assay);gc()

DefaultAssay(xenium_obj_filtered) <- "SCT"
xenium_obj_filtered = RunPCA(xenium_obj_filtered, npcs = 30)
xenium_obj_filtered = RunUMAP(xenium_obj_filtered, dims = 1:30)

save(xenium_obj_filtered, file="~/Desktop/Research/Analysis/spatial_transcriptome/sc/segmented/merged-sc-filtered-normalized-matched-nucleus-mini.rda")
