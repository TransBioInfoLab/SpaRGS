# run a small example for sGSVA
## use a small subsets from breast Visium dataset

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(msigdbr)
  library(GSVA); library(spdep); library(Matrix)
  library(matrixStats); library(jsonlite); library(readr)
  library(ggplot2); library(FNN)
})

source("R/sgsva_score.R")
source("R/spatial_permutation.R")

set.seed(1)

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
breast_results_dir <- file.path(project_dir, "results", "breast_visium")
results_dir <- file.path(project_dir, "results", "example")
data_dir <- file.path(project_dir, "data", "breast_visium")

matrixdir <- file.path(data_dir, "1160920F_filtered_count_matrix")
imgdir    <- file.path(data_dir, "1160920F_spatial")
metadir   <- file.path(data_dir, "1160920F_metadata.csv")

# check input
if (!dir.exists(data_dir)) {
  stop(
    "Breast Visium data directory not found: ", data_dir, "\n",
    "Please download the public 10x Visium breast cancer dataset ",
    "and place it in data/breast_visium/."
  )
}

# ---------------------------
# load data 
# ---------------------------
meta <- read.csv(metadir)

rna_data <- Read10X(
  matrixdir,
  gene.column = 1,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

# create seurat object from data dgCMatrix
bc_st <- CreateSeuratObject(counts = rna_data,
                            assay = "Spatial"
)
# read Visium image file
image <- Read10X_Image(imgdir,
                       image.name = "tissue_lowres_image.png",
                       filter.matrix = TRUE
)
image <- image[Cells(x = bc_st)]
DefaultAssay(object = image) <- "Spatial"
# add Visium data to image slot
bc_st@images$slice <- image

# add annotation to spatial metadata
meta$barcode <- meta$X
bc_st@meta.data$barcode <- rownames(bc_st@meta.data)
bc_st@meta.data <- bc_st@meta.data %>%
  dplyr::left_join(meta, by="barcode") %>%
  dplyr::select(-X)
rownames(bc_st@meta.data) <- bc_st@meta.data$barcode

# ---------------------------
# QC and normalization
# ---------------------------
DefaultAssay(bc_st) <- "Spatial"

bc_st[["percent.mt"]] <- PercentageFeatureSet(
  bc_st,
  assay = "Spatial",
  pattern = "^MT-"
)

qc_count_var <- paste0("nCount_", DefaultAssay(bc_st))
qc_feature_var <- paste0("nFeature_", DefaultAssay(bc_st))

required_qc_cols <- c(qc_count_var, qc_feature_var, "percent.mt")
missing_qc_cols <- setdiff(required_qc_cols, colnames(bc_st@meta.data))

if (length(missing_qc_cols) > 0) {
  stop(
    "Missing expected QC columns: ",
    paste(missing_qc_cols, collapse = ", ")
  )
}

keep_spots_qc <- rownames(bc_st@meta.data)[
  bc_st@meta.data[[qc_feature_var]] > 1000 &
    bc_st@meta.data[[qc_count_var]] < 30000 &
    bc_st@meta.data[["percent.mt"]] < 9
]

bc_subset <- subset(bc_st, cells = keep_spots_qc)
spatial_counts <- GetAssayData(bc_subset, assay = "Spatial", layer = "counts")
keep_genes <- rownames(spatial_counts)[Matrix::rowSums(spatial_counts > 0) >= 3]
bc_subset <- subset(bc_subset, features = keep_genes)

# SCTransform normalization
bc_norm <- SCTransform(bc_subset, assay = "Spatial", verbose = FALSE)
# clustering
bc_obj <- RunPCA(bc_norm, assay = "SCT", verbose = FALSE)
bc_obj <- FindNeighbors(bc_obj, reduction = "pca", dims = 1:30)
bc_obj <- FindClusters(bc_obj, verbose = FALSE, resolution = 0.5)
bc_obj <- RunUMAP(bc_obj, reduction = "pca", dims = 1:30)

# Compute local enrichment scores
expr_matrix <- as.matrix(SeuratObject::LayerData(bc_obj, assay = "SCT", layer = "data"))
nonzero_counts <- rowSums(expr_matrix > 0)
expr_filtered <- expr_matrix[nonzero_counts >= floor(0.05 * ncol(expr_matrix)), ,drop = FALSE]

# spatial coordinates
spatial_coords <- GetTissueCoordinates(bc_obj) %>%
  dplyr::select(x, y)

sf <- fromJSON(file.path(imgdir, "scalefactors_json.json"))
spot_diam_px <- sf$spot_diameter_fullres
spot_diam_um <- 65
um_per_px <- spot_diam_um / spot_diam_px
coords_um <- spatial_coords * um_per_px

# diagnostic
nn_um <- FNN::get.knn(as.matrix(coords_um), k = 1)$nn.dist[, 1]
median_nn_um <- median(nn_um[nn_um > 0 & is.finite(nn_um)], na.rm = TRUE)

message("Estimated microns per pixel: ", signif(um_per_px, 4))
message("Median nearest-neighbor distance after conversion: ", round(median_nn_um, 2), " um")

# KEGG gene sets
kegg_gene_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcollection = "CP:KEGG_LEGACY"
)
gene_sets <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)


# ---------------------------
# get a small subset
# ---------------------------
# choose a few representative pathways
gene_sets_small <- gene_sets[c(
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_GLUTATHIONE_METABOLISM",
  "KEGG_FATTY_ACID_METABOLISM"
)]


# ----------------------------------
# run sGSVA scoring and inference
# ----------------------------------
sigma_global <- round(estimate_sigma_knn(coords_um))
sgsva_res <- run_sgsva(expr_filtered, coords_um, gene_sets_small, sigma_global = sigma_global, 
                       mode = "delaunay", cor_min = 0.1, weight_clip = 0.99, min_genes = 15, max_n = 25, one_sided = "center_0.5")
sgsva_scores <- sgsva_res$scores
z_sgsva <- sgsva_res$z_score
balanced_sets <- sgsva_res$gene_sets
write.csv(as.data.frame(z_sgsva), file.path(results_dir, "minimal_example_visium_sgsva_zscores.csv"))


block_res <- spatial_permutation(z_sgsva, coords = coords_um, balanced_sets, 
                                 mode = "delaunay", sigma_global = sigma_global, one_sided = "center_0.5",
                                 perm_methods = c("hex_block"), hex_cellsize = 200)
hex_df <- block_res$permutation_results$hex_block_results$hex_block_summary
write.csv(hex_df, file.path(results_dir, "minimal_example_visium_hex_block_results.csv"), row.names = FALSE)
print(hex_df)

mini_res <- list(
  sgsva_results = sgsva_res, hex_block_permutation = block_res
)
saveRDS(mini_res, file.path(results_dir, "minimal_example_visium_all_results.rds"))

