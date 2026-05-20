# sGSVA analysis on breast cancer Visium data

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
data_dir <- file.path(project_dir, "data", "breast_visium")
results_dir <- file.path(project_dir, "results", "breast_visium")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

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
bc_st[["percent.mt"]] <- PercentageFeatureSet(bc_st, pattern = "^MT-")

bc_subset <- subset(
  bc_st,
  subset = nFeature_RNA > 1000 &
    nCount_RNA < 30000 &
    percent.mt < 9
)
keep_genes <- rownames(bc_subset)[rowSums(GetAssayData(bc_subset, layer = "counts") > 0) >= 3]

# Filter lowly expressed genes
bc_subset <- subset(bc_subset, features = keep_genes)

# SCTransform normalization
bc_norm <- SCTransform(bc_subset, assay = "Spatial", verbose = FALSE)
# clustering
bc_obj <- RunPCA(bc_norm, assay = "SCT", verbose = FALSE)
bc_obj <- FindNeighbors(bc_obj, reduction = "pca", dims = 1:30)
bc_obj <- FindClusters(bc_obj, verbose = FALSE, resolution = 0.5)
bc_obj <- RunUMAP(bc_obj, reduction = "pca", dims = 1:30)

# Compute local enrichment scores
expr_matrix <- as.matrix(bc_obj[["SCT"]]@data)
nonzero_counts <- rowSums(expr_matrix > 0)
expr_filtered <- expr_matrix[nonzero_counts >= floor(0.05 * ncol(expr_matrix)), ]

saveRDS(bc_obj, file.path(results_dir, "bc_obj.rds"))

# ---------------------------
# coordinates
# ---------------------------
spatial_coords <- GetTissueCoordinates(bc_obj) %>%
  dplyr::select(x, y)

sf <- fromJSON(file.path(imgdir, "scalefactors_json.json"))
spot_diam_px <- sf$spot_diameter_fullres
spot_diam_um <- 55
um_per_px <- spot_diam_um / spot_diam_px
coords_um <- spatial_coords * um_per_px

# ---------------------------
# gene sets
# ---------------------------
kegg_gene_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcollection = "CP:KEGG_LEGACY"
)
gene_sets <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)

# ---------------------------
# GSVA score
# ---------------------------
t0 <- proc.time()[3]
params <- gsvaParam(expr_filtered, gene_sets)
gsva_scores <- gsva(params, verbose = FALSE)
rt <- proc.time()[3] - t0

# ---------------------------
# sGSVA score
# ---------------------------
sigma_global <- round(estimate_sigma_knn(coords_um)) # 127
t0 <- proc.time()[3]
sgsva_res <- run_sgsva(expr_filtered, coords_um, gene_sets, sigma_global = sigma_global, 
                       mode = "delaunay", one_sided = "center_0.5")
rt <- proc.time()[3] - t0
saveRDS(sgsva_res, file.path(results_dir, "sgsva_score_res.rds"))
sgsva_scores <- sgsva_res$scores
z_sgsva <- sgsva_res$z_score
balanced_sets <- sgsva_res$gene_sets

write.csv(
  sgsva_scores,
  file.path(results_dir, "sgsva_scores.csv"),
  row.names = TRUE
)

write.csv(
  gsva_scores,
  file.path(results_dir, "gsva_scores.csv"),
  row.names = TRUE
)


# ----------------------------------
# statistical tests across regions
# ----------------------------------
# statistical tests for mean z-score across regions
t_sgsva <- as.data.frame(t(z_sgsva))
t_sgsva$barcode <- rownames(t_sgsva)

merge_dat <- merge(t_sgsva, bc_obj@meta.data[,c("barcode","Classification")], by="barcode")
results <- lapply(rownames(z_sgsva), function(pathway) {
  df <- merge_dat %>% dplyr::select(Classification, !!sym(pathway))
  colnames(df)[2] <- "score"
  
  # Compute mean z per region
  region_means <- df %>% 
    group_by(Classification) %>% 
    summarise(mean_z = mean(score, na.rm=TRUE)) %>%
    tidyr::pivot_wider(names_from = Classification, values_from = mean_z, names_prefix = "mean_z_")
  # ANOVA or Kruskal-Wallis test
  kw <- kruskal.test(score ~ Classification, data = df)  # Use a non-parametric test
  # Combine
  cbind(
    pathway = pathway,
    pval = kw$p.value,
    region_means
  )
})
# Combine all rows into one data frame
stats_table <- bind_rows(results) %>%
  mutate(fdr = p.adjust(pval, method = "BH")) %>%
  relocate(fdr, .after=pval)
write.csv(stats_table, file.path(results_dir, "region_overlays_meanz_fdr.csv"))


# ---------------------------
# spatial permutation
# ---------------------------
sgsva_res <- spatial_permutation(z_sgsva, coords = coords_um, balanced_sets, 
                                 mode = "delaunay", sigma_global = sigma_global, 
                                 one_sided = "center_0.5", perm_methods = c("moran_mc","geary_mc"))
saveRDS(sgsva_res, file.path(results_dir, "moranmc_res.rds"))

mc_df <- sgsva_res$permutation_results$moran_geary_results$moran_geary_summary
out_csv <- file.path(results_dir, "moran_geary_pathway_summary.csv")
readr::write_csv(mc_df, out_csv)

block_sizes <- c(150, 200, 300)
block_res <- lapply(block_sizes, function(b) {
  spatial_permutation(z_sgsva, coords = coords_um, balanced_sets, 
                      mode = "delaunay", sigma_global = sigma_global, one_sided = "center_0.5",
                      perm_methods = c("hex_block"), hex_cellsize = b)
})
saveRDS(block_res, file.path(results_dir, "block_res.rds"))

hex150_df <- block_res[[1]]$permutation_results$hex_block_results$hex_block_summary
readr::write_csv(hex150_df, file.path(results_dir, "block150_pathway_summary.csv"))

hex200_df <- block_res[[2]]$permutation_results$hex_block_results$hex_block_summary
readr::write_csv(hex200_df, file.path(results_dir, "block200_pathway_summary.csv"))

hex300_df <- block_res[[3]]$permutation_results$hex_block_results$hex_block_summary
readr::write_csv(hex300_df, file.path(results_dir, "block300_pathway_summary.csv"))


