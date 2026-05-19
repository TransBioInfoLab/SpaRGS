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

# require files
required_files <- c(
  file.path(breast_results_dir, "bc_obj.rds"),
  file.path(imgdir, "scalefactors_json.json")
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input files for the minimal example:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nRun the Breast Visium preprocessing script first or place the required files in the expected folders."
  )
}

# ---------------------------
# load data 
# ---------------------------
bc_obj <- readRDS(file.path(breast_results_dir, "bc_obj.rds"))

# Compute local enrichment scores
expr_matrix <- as.matrix(bc_obj[["SCT"]]@data)
nonzero_counts <- rowSums(expr_matrix > 0)
expr_filtered <- expr_matrix[nonzero_counts >= floor(0.05 * ncol(expr_matrix)), ]

# spatial coordinates
spatial_coords <- GetTissueCoordinates(bc_obj) %>%
  dplyr::select(x, y)

sf <- fromJSON(file.path(imgdir, "scalefactors_json.json"))
spot_diam_px <- sf$spot_diameter_fullres
spot_diam_um <- 55
um_per_px <- spot_diam_um / spot_diam_px
coords_um <- spatial_coords * um_per_px

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
# choose a manageable number of spots
n_example <- 200

# use metadata from the QC-passed object
meta_example <- bc_obj@meta.data
meta_example$barcode <- rownames(meta_example)

# keep only barcodes that exist in expr_filtered
meta_example <- meta_example[meta_example$barcode %in% colnames(expr_filtered), , drop = FALSE]

# stratified sampling by Classification if available
if ("Classification" %in% colnames(meta_example)) {
  sampled_barcodes <- meta_example %>%
    dplyr::group_by(Classification) %>%
    dplyr::group_modify(~ .x %>% dplyr::slice_sample(n = min(20, nrow(.x)))) %>%
    dplyr::ungroup() %>%
    dplyr::pull(barcode)
  
  if (length(sampled_barcodes) > n_example) {
    sampled_barcodes <- sample(sampled_barcodes, n_example)
  } else if (length(sampled_barcodes) < n_example) {
    remaining <- setdiff(meta_example$barcode, sampled_barcodes)
    extra <- sample(remaining, min(n_example - length(sampled_barcodes), length(remaining)))
    sampled_barcodes <- c(sampled_barcodes, extra)
  }
} else {
  sampled_barcodes <- sample(colnames(expr_filtered), min(n_example, ncol(expr_filtered)))
}

# subset expression
expr_small <- expr_filtered[, sampled_barcodes, drop = FALSE]

# subset coordinates and reorder to match expression columns
coords_df <- as.data.frame(coords_um)
coords_df$barcode <- rownames(coords_df)
coords_small <- coords_df[match(sampled_barcodes, coords_df$barcode), c("x", "y"), drop = FALSE]
rownames(coords_small) <- sampled_barcodes

# choose a few representative pathways
gene_sets_small <- gene_sets[c(
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_GLUTATHIONE_METABOLISM",
  "KEGG_FATTY_ACID_METABOLISM"
)]

# subset metadata
meta_small <- meta_example[match(sampled_barcodes, meta_example$barcode), , drop = FALSE]

# save minimal example
mini_exm <- list(
  expr_small = expr_small, coords_small = coords_small,
  meta_small = meta_small, gene_sets_small = gene_sets_small
)
saveRDS(mini_exm, file.path(results_dir, "minimal_example_breast_visium.rds"))


# ----------------------------------
# run sGSVA scoring and inference
# ----------------------------------
sigma_global <- estimate_sigma_knn(coords_small)
sgsva_res <- run_sgsva(expr_small, coords_small, gene_sets_small, 
                       sigma_global = sigma_global, one_sided = "center_0.5")
sgsva_scores <- sgsva_res$scores
z_sgsva <- sgsva_res$z_score
balanced_sets <- sgsva_res$gene_sets
print(sgsva_scores)

block_res <- spatial_permutation(z_sgsva, coords = coords_small, balanced_sets, 
                                 mode = "delaunay", sigma_global = sigma_global, one_sided = "center_0.5",
                                 perm_methods = c("hex_block"), hex_cellsize = 300)
hex_df <- block_res$hex_block_results$hex_block_summary
print(hex_df)

mini_res <- list(
  sgsva_results = sgsva_res, hex_block_permutation = block_res
)

write.csv(as.data.frame(z_sgsva), file.path(results_dir, "minimal_example_visium_sgsva_z_scores.csv"))
write.csv(hex_df, file.path(results_dir, "minimal_example_visium_hex_block_results.csv"), row.names = FALSE)
saveRDS(mini_res, file.path(results_dir, "minimal_example_visium_all_results.rds"))

