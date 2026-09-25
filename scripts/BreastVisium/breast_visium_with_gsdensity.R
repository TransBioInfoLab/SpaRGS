# spatial pathway comparison on breast cancer Visium data

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(dplyr); library(tidyr)
  library(purrr); library(tibble); library(readr); library(stringr)
  library(msigdbr); library(GSVA); library(spdep); library(Matrix)
  library(matrixStats); library(jsonlite); library(ggplot2); library(ggspatial)
  library(patchwork); library(scales); library(FNN); library(sf)
  library(gsdensity); library(CelliD); library(here); library(fs)
})

set.seed(1)

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
data_dir <- file.path(project_dir, "data", "breast_visium")
results_dir <- file.path(project_dir, "results", "breast_visium")
qc_dir <- file.path(results_dir, "reproducibility")
supp_dir <- file.path(project_dir, "results", "supplementary_data")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

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

source(file.path(project_dir, "R", "sgsva_score.R"))
source(file.path(project_dir, "R", "spatial_permutation.R"))

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

saveRDS(bc_obj, file.path(results_dir, "bc_obj.rds"))

# Compute local enrichment scores
expr_matrix <- as.matrix(SeuratObject::LayerData(bc_obj, assay = "SCT", layer = "data"))
nonzero_counts <- rowSums(expr_matrix > 0)
expr_filtered <- expr_matrix[nonzero_counts >= floor(0.05 * ncol(expr_matrix)), ,drop = FALSE]

# ---------------------------
# coordinates
# ---------------------------
spatial_coords <- GetTissueCoordinates(bc_obj)
spatial_coords <- spatial_coords[colnames(expr_filtered), c("x", "y")]
stopifnot(identical(rownames(spatial_coords), colnames(expr_filtered)))

sf <- fromJSON(file.path(imgdir, "scalefactors_json.json"))
spot_diam_px <- sf$spot_diameter_fullres
spot_diam_um <- 65
um_per_px <- spot_diam_um / spot_diam_px
coords_um <- spatial_coords * um_per_px

saveRDS(
  list(expr_filtered=expr_filtered, coords_um=coords_um),
  file.path(results_dir, "breast_visium_sct_data.rds")
)

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
# generate runtime table
runtime_results <- tibble::tibble(
  method = character(),
  runtime = numeric()
)

t0 <- proc.time()[3]
params <- gsvaParam(expr_filtered, gene_sets)
gsva_scores <- gsva(params, verbose = FALSE)
rt <- proc.time()[3] - t0
runtime_results <- bind_rows(runtime_results,
                             tibble(method = "GSVA", runtime = rt))

# ---------------------------
# sGSVA score
# ---------------------------
sigma_global <- round(estimate_sigma_knn(coords_um)) # 150
t0 <- proc.time()[3]
sgsva_res <- run_sgsva(expr_filtered, coords_um, gene_sets, sigma_global = sigma_global, 
                       mode = "delaunay", one_sided = "center_0.5")
rt <- proc.time()[3] - t0
runtime_results <- bind_rows(runtime_results,
                             tibble(method = "sGSVA", runtime = rt))

sgsva_scores <- sgsva_res$scores
z_sgsva <- sgsva_res$z_score
balanced_sets <- sgsva_res$gene_sets

# ---------------------------
# GSVA-smoothed score
# ---------------------------
t0 <- proc.time()[3]
smoothed_expr <- smooth_expr(expr_filtered, coords_um, sigma = sigma_global, mode = "delaunay")
params_smoothed <- gsvaParam(as.matrix(smoothed_expr), gene_sets)
gsva_smoothed_scores <- gsva(params_smoothed, verbose = FALSE)
rt <- proc.time()[3] - t0
runtime_results <- bind_rows(runtime_results,
                             tibble(method = "GSVA-smoothed", runtime = rt))

# ---------------------------
# GSDensity
# ---------------------------
# prepare matched Seurat object for GSDensity
genes_use <- rownames(expr_filtered)
spots_use <- colnames(expr_filtered)
stopifnot(identical(spots_use, rownames(coords_um)))

# Create a separate object so that GSDensity does not modify bc_obj
bc_gsdensity <- bc_obj
DefaultAssay(bc_gsdensity) <- "Spatial"

bc_gsdensity <- subset(bc_gsdensity, cells = spots_use, features = genes_use)
  
stopifnot(
  identical(colnames(bc_gsdensity), spots_use),
  identical(rownames(bc_gsdensity), genes_use)
)

# adjust original function to solve the compatibility issue
compute_mca_v5 <- function(
    object,
    dims.use = 1:10,
    genes.use = rownames(object),
    assay = Seurat::DefaultAssay(object),
    verbose = FALSE
) {
  # Preserve the original compute.mca() logic:
  # 1. intersect genes
  # 2. NormalizeData
  # 3. RunMCA
  # 4. combine gene and spot embeddings
  
  Seurat::DefaultAssay(object) <- assay
  
  genes.use.1 <- intersect(genes.use, rownames(object))
  
  # Same normalization used by the original compute.mca()
  object <- Seurat::NormalizeData(object = object)
  
  # Seurat v5-compatible replacement for:
  # GetAssayData(object, slot = "data")
  expression_data <- SeuratObject::LayerData(
    object = object, assay = assay, layer = "data")
  
  # Subset before converting to a dense matrix
  expression_data <- expression_data[genes.use.1,,drop = FALSE]
  
  # CelliD's matrix method is compatible because it does not access
  # a Seurat slot/layer internally.
  mca_result <- CelliD::RunMCA(
    X = as.matrix(expression_data),
    nmcs = max(dims.use),
    reduction.name = "MCA"
  )
  
  genes.use.update <- intersect(genes.use.1,
    rownames(mca_result$featuresCoordinates))
  
  # Same returned object as gsdensity::compute.mca():
  # gene coordinates followed by spot coordinates
  coembed <- rbind(
    mca_result$featuresCoordinates[genes.use.update, dims.use, drop = FALSE],
    mca_result$cellsCoordinates[, dims.use, drop = FALSE]
  )
  
  return(coembed)
}

# MCA co-embedding
set.seed(1)

t0 <- proc.time()[["elapsed"]]

ce <- compute_mca_v5(
  object = bc_gsdensity, dims.use = 1:10,
  genes.use = genes_use, assay = "Spatial"
)
genes_in_mca <- intersect(rownames(ce), rownames(bc_gsdensity))

# Keep gene sets with at least three genes in the MCA graph.
# This matches the default gene.set.cutoff used by GSDensity.
gene_sets_gsdensity <- lapply(
  gene_sets,
  function(gs) intersect(unique(gs), genes_in_mca)
)

gene_sets_gsdensity <- gene_sets_gsdensity[
  lengths(gene_sets_gsdensity) >= 3
]
message("GSDensity retained ", length(gene_sets_gsdensity), " of ",
  length(gene_sets), " pathways.")

# Build nearest-neighbor graph in MCA space.
el <- gsdensity::compute.nn.edges(coembed = ce, nn.use = 300)

# propagation probabilities for all retained pathways
gsdensity_spot_scores <- gsdensity::run.rwr.list(
  el = el,
  gene_set_list = gene_sets_gsdensity,
  cells = spots_use,
  restart = 0.75
)
rt <- proc.time()[3] - t0
runtime_results <- bind_rows(runtime_results, tibble(
  method = "GSDensity", runtime = rt
))

gsdensity_spot_scores <- as.matrix(gsdensity_spot_scores)
gsdensity_spot_scores <- gsdensity_spot_scores[spots_use,,drop = FALSE]

# Convert to pathways × spots, matching GSVA and sGSVA.
gsdensity_scores <- t(gsdensity_spot_scores)

# save all results
write.csv(gsdensity_scores, file.path(results_dir, "gsdensity_scores.csv"), row.names = TRUE)
write.csv(runtime_results, file.path(results_dir, "pathway_methods_runtime.csv"), row.names = FALSE)

saveRDS(sgsva_res, file.path(results_dir, "sgsva_score_res.rds"))
write.csv(sgsva_scores, file.path(results_dir, "sgsva_scores.csv"), row.names = TRUE)
write.csv(gsva_scores, file.path(results_dir, "gsva_scores.csv"), row.names = TRUE)
write.csv(gsva_smoothed_scores, file.path(results_dir, "gsva_smoothed_scores.csv"), row.names = TRUE)


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
# global permutation
t0 <- proc.time()[3]
sgsva_res <- spatial_permutation(z_sgsva, coords = coords_um, balanced_sets, 
                                 mode = "delaunay", sigma_global = sigma_global, 
                                 one_sided = "center_0.5", perm_methods = c("moran_mc","geary_mc"))
rt <- proc.time()[3] - t0 
print(rt) # 14542.04 
saveRDS(sgsva_res, file.path(results_dir, "moranmc_res.rds"))

mc_df <- sgsva_res$permutation_results$moran_geary_results$moran_geary_summary
out_csv <- file.path(results_dir, "moran_geary_pathway_summary.csv")
readr::write_csv(mc_df, out_csv)

# export output
mc_export <- mc_df %>%
  dplyr::mutate(
    permutation_scheme = "global",
    one_sided = "center_0.5",
    seed_permutation = 124
  ) %>%
  dplyr::select(-obs_geary_flipped, -z_geary_flipped, 
                -null_mean_geary_flipped, -null_sd_geary_flipped)

# hex-block permutation
block_sizes <- c(150, 200, 300)
block_res <- lapply(block_sizes, function(b) {
  t0 <- proc.time()[3]
  spatial_permutation(z_sgsva, coords = coords_um, balanced_sets, 
                      mode = "delaunay", sigma_global = sigma_global, one_sided = "center_0.5",
                      perm_methods = c("hex_block"), hex_cellsize = b)
  rt <- proc.time()[3] - t0 #19364.64
})
saveRDS(block_res, file.path(results_dir, "block_res.rds"))

hex150_df <- block_res[[1]]$permutation_results$hex_block_results$hex_block_summary
readr::write_csv(hex150_df, file.path(results_dir, "block150_pathway_summary.csv"))

hex200_df <- block_res[[2]]$permutation_results$hex_block_results$hex_block_summary
readr::write_csv(hex200_df, file.path(results_dir, "block200_pathway_summary.csv"))

hex300_df <- block_res[[3]]$permutation_results$hex_block_results$hex_block_summary
readr::write_csv(hex300_df, file.path(results_dir, "block300_pathway_summary.csv"))

# export output
hex150_export <- hex150_df %>%
  dplyr::mutate(
    permutation_scheme = "hex-block",
    block_size_um = block_size,
    one_sided = "center_0.5",
    seed_permutation = 124
  ) %>%
  dplyr::select(-obs_geary_flipped, -z_geary_flipped, 
                -null_mean_geary_flipped, -null_sd_geary_flipped, -block_size)
hex200_export <- hex200_df %>%
  dplyr::mutate(
    permutation_scheme = "hex-block",
    block_size_um = block_size,
    one_sided = "center_0.5",
    seed_permutation = 124
  ) %>%
  dplyr::select(-obs_geary_flipped, -z_geary_flipped, 
                -null_mean_geary_flipped, -null_sd_geary_flipped, -block_size)
hex300_export <- hex300_df %>%
  dplyr::mutate(
    permutation_scheme = "hex-block",
    block_size_um = block_size,
    one_sided = "center_0.5",
    seed_permutation = 124
  ) %>%
  dplyr::select(-obs_geary_flipped, -z_geary_flipped, 
                -null_mean_geary_flipped, -null_sd_geary_flipped, -block_size)


# ---------------------------
# Gene-set coverage
# ---------------------------
normalize_pathway_names <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  sub("^KEGG_", "", x)
}

gene_set_coverage_base <- tibble(
  pathway_original = names(gene_sets),
  pathway = normalize_pathway_names(names(gene_sets)),
  n_genes_original = lengths(gene_sets),
  n_genes_in_expr_filtered = lengths(lapply(gene_sets, function(gs) {
    intersect(unique(gs), rownames(expr_filtered))
  }))
)

score_pathway_sets <- list(
  sGSVA = normalize_pathway_names(rownames(sgsva_scores)),
  GSVA = normalize_pathway_names(rownames(gsva_scores)),
  `GSVA-smoothed` = normalize_pathway_names(rownames(gsva_smoothed_scores)),
  GSDensity = normalize_pathway_names(rownames(gsdensity_scores))
)

gene_set_coverage_main <- purrr::imap_dfr(
  score_pathway_sets,
  function(pathways_scored, method_name) {
    gene_set_coverage_base %>%
      mutate(
        method = method_name,
        scored_by_method = pathway %in% pathways_scored,
        n_pathways_scored_by_method = sum(scored_by_method),
        .before = 1
      )
  }
)

readr::write_csv(
  gene_set_coverage_main,
  file.path(qc_dir, "breast_gene_set_coverage_main_methods.csv")
)


################################################################################
# Method-wise Moran's I and Geary's C comparisons
# Collect score matrices
score_list <- list(
  sGSVA = as.matrix(sgsva_scores),
  GSVA = as.matrix(gsva_scores),
  GSVA_smoothed = as.matrix(gsva_smoothed_scores),
  GSDensity = as.matrix(gsdensity_scores)
)

# standardize barcodes
standardize_barcodes <- function(x) {
  gsub(".", "-", x, fixed = TRUE)
}

score_list <- lapply(score_list, function(mat) {
  rownames(mat) <- normalize_pathway_names(rownames(mat))
  colnames(mat) <- standardize_barcodes(colnames(mat))
  mat
})

# Restrict all methods to the same spots and shared pathways
common_spots <- Reduce(intersect,
  c(lapply(score_list, colnames),list(rownames(coords_um))))

shared_pathways <- Reduce(intersect, lapply(score_list, rownames))

coords_common <- as.data.frame(coords_um[common_spots, , drop = FALSE])
coords_common <- coords_common[, 1:2, drop = FALSE]

score_list <- lapply(score_list, function(mat) {
  mat[shared_pathways, common_spots, drop = FALSE]
})

message("Common spots: ", length(common_spots))
message("Shared pathways: ", length(shared_pathways))

# Build one common Delaunay graph for every method
nb <- spdep::tri2nb(as.matrix(coords_common))
lw <- spdep::nb2listw(
  nb,
  style = "W",
  zero.policy = TRUE
)

# Compute observed Moran's I and Geary's C
spatial_stats <- imap_dfr(score_list, function(mat, method_name) {
  map_dfr(shared_pathways, function(pathway_name) {
    x <- as.numeric(mat[pathway_name, ])
    
    if (all(!is.finite(x)) || sd(x, na.rm = TRUE) == 0) {
      return(tibble(
        method = method_name,
        pathway = pathway_name,
        moran_i = NA_real_,
        geary_c = NA_real_
      ))
    }
    
    moran_value <- spdep::moran.test(x, lw, zero.policy = TRUE)$estimate[["Moran I statistic"]]
    geary_value <- spdep::geary.test(x, lw, zero.policy = TRUE)$estimate[["Geary C statistic"]]
    
    tibble(
      method = method_name,
      pathway = pathway_name,
      moran_i = moran_value,
      geary_c = geary_value
    )
  })
})

# save outputs
readr::write_csv(
  spatial_stats,
  file.path(results_dir, "pathway_level_moran_geary_shared_pathways.csv")
)


# add input harmonization
main_input_harmonization <- tibble(
  method = c("sGSVA", "GSVA", "GSVA-smoothed", "GSDensity"),
  count_source = c(
    "Same breast Visium Spatial assay",
    "Same breast Visium Spatial assay",
    "Same breast Visium Spatial assay",
    "Same breast Visium Spatial assay"
  ),
  expression_input = c(
    "Filtered SCT data matrix",
    "Filtered SCT data matrix",
    "Gaussian-smoothed filtered SCT data matrix",
    "Method-specific Normalize followed by MCA"
  ),
  retained_spots_before_shared_filter = c(
    ncol(sgsva_scores),
    ncol(gsva_scores),
    ncol(gsva_smoothed_scores),
    ncol(gsdensity_scores)
  ),
  retained_pathways_before_shared_filter = c(
    nrow(sgsva_scores),
    nrow(gsva_scores),
    nrow(gsva_smoothed_scores),
    nrow(gsdensity_scores)
  ),
  common_spots_used_for_spatial_stats = length(common_spots),
  shared_pathways_used_for_spatial_stats = length(shared_pathways),
  spatial_graph_for_moran_geary = "Common Delaunay graph from shared spot coordinates",
  score_scaling_for_maps = "Pathway-wise z-standardization for visualization only",
  notes = c(
    "Data-driven signed pathway scoring",
    "Classical pathway score baseline",
    "Smoothing-only baseline",
    "MCA/RWR workflow with Seurat v5 compatibility wrapper"
  )
)

readr::write_csv(
  main_input_harmonization,
  file.path(qc_dir, "breast_input_harmonization_main_methods.csv")
)



################################################################################
# save package versions and session info
writeLines(
  capture.output(sessionInfo()),
  con = file.path(qc_dir, paste0(tools::file_path_sans_ext(basename(commandArgs(trailingOnly = FALSE)[1])), "_sessionInfo.txt"))
)

packages_to_record <- c(
  "Seurat", "SeuratObject", "GSVA", "gsdensity", "CelliD",
  "spacexr", "spdep", "msigdbr", "Matrix", "dplyr",
  "tidyr", "purrr", "ggplot2", "readr", "jsonlite"
)

package_versions <- tibble(
  package = packages_to_record,
  version = purrr::map_chr(
    packages_to_record,
    ~ {
      if (requireNamespace(.x, quietly = TRUE)) {
        as.character(utils::packageVersion(.x))
      } else {
        NA_character_
      }
    }
  )
)

readr::write_csv(
  package_versions,
  file.path(qc_dir, "breast_benchmark_package_versions.csv")
)


library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "SData3_histology_associations")
writeData(wb, "SData3_histology_associations", stats_table)
addWorksheet(wb, "SData3_global_permutation")
writeData(wb, "SData3_global_permutation", mc_df)
addWorksheet(wb, "SData3_block150_permutation")
writeData(wb, "SData3_block150_permutation", hex150_df)
addWorksheet(wb, "SData3_block200_permutation")
writeData(wb, "SData3_block200_permutation", hex200_df)
addWorksheet(wb, "SData3_block300_permutation")
writeData(wb, "SData3_block300_permutation", hex300_df)
addWorksheet(wb, "SData3_geneset_coverage")
writeData(wb, "SData3_geneset_coverage", gene_set_coverage_main)
saveWorkbook(wb, file.path(supp_dir, "Supplementary Data 3.xlsx"), overwrite = TRUE)


