### Real data showcase2 for DLPFC
suppressPackageStartupMessages ({
  library(spatialLIBD); library(SpatialExperiment)
  library(scater); library(scran); library(ggplot2)
  library(viridis); library(spdep); library(msigdbr)
  library(dplyr); library(purrr) ; library(FNN)
  library(GSVA); library(org.Hs.eg.db)
})
source("R/sgsva_score.R")
source("R/spatial_permutation.R")
set.seed(1)

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
results_dir <- file.path(project_dir, "results", "DLPFC")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------
# Load data and gene sets
# ---------------------------
spe <- fetch_data(type="spe")
sample_ids <- unique(spe$sample_id)
spe_list <- lapply(sample_ids, function(id) spe[, spe$sample_id == id])
names(spe_list) <- sample_ids

# extract KEGG pathway sets
kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
gene_sets <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)

# map ensembl to gene symbol
all_genes <- rownames(counts(spe))
ensembl_ids <- gsub("\\..*", "", all_genes)
symbols <- mapIds(org.Hs.eg.db,
                  keys = ensembl_ids,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first") # map to HGNC symbols
map_table <- tibble(ensembl = ensembl_ids, symbol = symbols) %>%
  filter(!is.na(symbol) & symbol != "") %>%
  distinct(ensembl, .keep_all = TRUE)

# helper function to convert one matrix
convert_rows_to_symbols <- function(expr, map_table) {
  rn <- gsub("\\..*", "", rownames(expr))
  matched <- map_table$symbol[match(rn, map_table$ensembl)]
  expr <- expr[!is.na(matched) & matched != "", ]
  matched <- matched[!is.na(matched) & matched != ""]
  rownames(expr) <- matched
  expr <- expr[!duplicated(rownames(expr)), ]
  expr
}

# Helper function to extract and preprocess per sample
preprocess_per_sample <- function(spe_sub, map_table) {
  # Convert Ensembl → gene symbol
  counts_mat <- as.matrix(assay(spe_sub, "counts"))
  counts_mat <- convert_rows_to_symbols(counts_mat, map_table)
  
  # Filter lowly expressed genes
  keep <- rowSums(counts_mat > 0) >= floor(0.05 * ncol(counts_mat))
  expr <- counts_mat[keep, ]
  
  # Normalize to logCPM
  lib <- colSums(expr)
  expr_logCPM <- log1p(t(t(expr) / pmax(lib, 1)) * 1e6)
  
  # Extract spatial coordinates
  coords <- as.data.frame(spatialCoords(spe_sub))[, c("pxl_col_in_fullres", "pxl_row_in_fullres")]
  colnames(coords) <- c("x", "y")
  
  list(
    sample_id = unique(spe_sub$sample_id),
    expr = expr_logCPM,
    coords = coords,
    counts_mat = counts_mat
  )
}
dlpfc_preprocessed <- purrr::map(spe_list, preprocess_per_sample, map_table = map_table)
saveRDS(dlpfc_preprocessed, file.path(results_dir, "DLPFC_preprocess.rds"))



# ----------------------------
# get sGSVA and GSVA score
# ----------------------------
# Compute base σ for each sample
sigma_table <- sapply(unique(spe$sample_id), function(id) {
  coords <- as.data.frame(spatialCoords(spe[, spe$sample_id == id]))
  estimate_sigma_knn(coords)
})
sigma_table # all sigmas are around 207

# run sGSVA and GSVA across 12 samples
full_res <- list()
for (i in 1:length(sample_ids)) {
  id <- sample_ids[i]
  # GSVA
  params <- gsvaParam(
    dlpfc_preprocessed[[id]]$expr, gene_sets
  )
  gsva_scores <- gsva(params, verbose = FALSE)
  # sGSVA
  t0 <- proc.time()[3]
  sgsva_res <- run_sgsva(dlpfc_preprocessed[[id]]$expr, 
                         dlpfc_preprocessed[[id]]$coords, gene_sets,
                         sigma_global = 207, cor_min = 0.2,
                         mode = "delaunay", one_sided = "center_0.5")
  rt <- proc.time()[3] - t0
  full_res[[id]] <- list(
    gsva_scores = gsva_scores,
    sgsva_res = sgsva_res,
    runtime = rt
  )
}
saveRDS(full_res, file.path(results_dir, "full_results_12_samples.rds"))

# record spots, runtimes, and pathways per section
spot <- NULL
runtime <- NULL
num_path <- NULL
for (i in 1:length(sample_ids)){
  id <- sample_ids[i]
  spot <- c(spot, dim(dlpfc_preprocessed[[id]]$coords)[1])
  runtime <- c(runtime, full_res[[id]]$runtime)
  
  sgsva_res <- full_res[[id]]$sgsva_res
  balanced_sets <- sgsva_res$gene_sets
  num_path <- c(num_path, length(balanced_sets))
}
#print(spot)
#print(runtime)
#print(num_path)



# ---------------------------
# spatial permutation
# ---------------------------
for (i in 1:length(sample_ids)) {
  id <- sample_ids[i]
  # get sGSVA scores
  sgsva_res <- full_res[[id]]$sgsva_res
  z_sgsva <- sgsva_res$z_score
  balanced_sets <- sgsva_res$gene_sets
  
  # get coords
  idx <- spe$sample_id == id
  coords <- as.data.frame(spatialCoords(spe)[idx, ])
  colnames(coords)[1:2] <- c("x","y")
  
  sgsva_res <- spatial_permutation(z_sgsva, coords, balanced_sets, id, mode = "delaunay",
                                   perm_methods = c("moran_mc","geary_mc"))
  saveRDS(sgsva_res, file.path(results_dir, paste0(id,"_moranmc_res.rds")))
  
  block_sizes <- c(150, 200, 300)
  block_res <- lapply(block_sizes, function(b) {
    spatial_permutation(z_sgsva, coords, balanced_sets, id, mode = "delaunay",
                        perm_methods = c("hex_block"), hex_cellsize = b)
  })
  saveRDS(block_res, file.path(results_dir, paste0(id,"_block_res.rds")))
}

# meta-analysis across 12 samples
mc_all <- NULL
hex150_all <- NULL
hex200_all <- NULL
hex300_all <- NULL
for (i in 1:length(sample_ids)) {
  id <- sample_ids[i]
  sgsva_res <- readRDS(file.path(results_dir, paste0(id,"_moranmc_res.rds")))
  block_res <- readRDS(file.path(results_dir, paste0(id,"_block_res.rds")))
  
  mc_df <- sgsva_res$moran_geary_results$moran_geary_summary
  mc_df$sample_id <- id
  mc_all <- rbind(mc_all, mc_df)
  hex150_df <- block_res[[1]]$hex_block_results$hex_block_summary
  hex150_df$sample_id <- id
  hex150_all <- rbind(hex150_all, hex150_df)
  hex200_df <- block_res[[2]]$hex_block_results$hex_block_summary
  hex200_df$sample_id <- id
  hex200_all <- rbind(hex200_all, hex200_df)
  hex300_df <- block_res[[3]]$hex_block_results$hex_block_summary
  hex300_df$sample_id <- id
  hex300_all <- rbind(hex300_all, hex300_df)
}

mc_all$one_sided <- "center_0.5"
hex150_all$one_sided <- "center_0.5"
hex200_all$one_sided <- "center_0.5"
hex300_all$one_sided <- "center_0.5"
write.csv(mc_all, file.path(outdir, "mc_meta_stouffer_res_all_sections.csv"))
write.csv(hex150_all, file.path(outdir, "hex150_meta_stouffer_res_all_sections.csv"))
write.csv(hex200_all, file.path(outdir, "hex200_meta_stouffer_res_all_sections.csv"))
write.csv(hex300_all, file.path(outdir, "hex300_meta_stouffer_res_all_sections.csv"))


# Stouffer meta 
mc_meta <- mc_all %>%
  group_by(pathway) %>%
  summarize(
    k = n(),
    n_sig = sum(p_adj_moran < 0.05, na.rm = TRUE),
    median_obs_moran = median(obs_moran),
    median_z_moran = median(z_moran),
    median_obs_geary = median(obs_geary),
    median_z_geary = median(z_geary),
    z_stouffer = sum(qnorm(1 - p_moran / 2)) / sqrt(n()),
    p_meta = 2 * pnorm(-abs(z_stouffer)),
    .groups = "drop"
  ) %>%
  mutate(fdr_meta = p.adjust(p_meta, method = "BH")) %>%
  arrange(desc(median_obs_moran))
write.csv(mc_meta, file.path(outdir, "mc_meta_stouffer_res.csv"))

hex150_meta <- hex150_all %>%
  group_by(pathway) %>%
  summarize(
    k = n(),
    n_sig = sum(p_adj_moran < 0.05, na.rm = TRUE),
    median_obs_moran = median(obs_moran),
    median_z_moran = median(z_moran),
    median_obs_geary = median(obs_geary),
    median_z_geary = median(z_geary),
    z_stouffer = sum(qnorm(1 - p_moran / 2)) / sqrt(n()),
    p_meta = 2 * pnorm(-abs(z_stouffer)),
    .groups = "drop"
  ) %>%
  mutate(fdr_meta = p.adjust(p_meta, method = "BH")) %>%
  arrange(desc(median_obs_moran))
write.csv(hex150_meta, file.path(outdir, "hex150_meta_stouffer_res.csv"))

hex200_meta <- hex200_all %>%
  group_by(pathway) %>%
  summarize(
    k = n(),
    n_sig = sum(p_adj_moran < 0.05, na.rm = TRUE),
    median_obs_moran = median(obs_moran),
    median_z_moran = median(z_moran),
    median_obs_geary = median(obs_geary),
    median_z_geary = median(z_geary),
    z_stouffer = sum(qnorm(1 - p_moran / 2)) / sqrt(n()),
    p_meta = 2 * pnorm(-abs(z_stouffer)),
    .groups = "drop"
  ) %>%
  mutate(fdr_meta = p.adjust(p_meta, method = "BH")) %>%
  arrange(desc(median_obs_moran))
write.csv(hex200_meta, file.path(outdir, "hex200_meta_stouffer_res.csv"))

hex300_meta <- hex300_all %>%
  group_by(pathway) %>%
  summarize(
    k = n(),
    n_sig = sum(p_adj_moran < 0.05, na.rm = TRUE),
    median_obs_moran = median(obs_moran),
    median_z_moran = median(z_moran),
    median_obs_geary = median(obs_geary),
    median_z_geary = median(z_geary),
    z_stouffer = sum(qnorm(1 - p_moran / 2)) / sqrt(n()),
    p_meta = 2 * pnorm(-abs(z_stouffer)),
    .groups = "drop"
  ) %>%
  mutate(fdr_meta = p.adjust(p_meta, method = "BH")) %>%
  arrange(desc(median_obs_moran))
write.csv(hex300_meta, file.path(outdir, "hex300_meta_stouffer_res.csv"))


# save combined hex block permutation summary
for (i in 1:length(sample_ids)) {
  id <- sample_ids[i]
  block_res <- readRDS(file.path(results_dir, paste0(id,"_block_res.rds")))
  hex150_df <- block_res[[1]]$hex_block_results$hex_block_summary
  hex200_df <- block_res[[2]]$hex_block_results$hex_block_summary
  hex300_df <- block_res[[3]]$hex_block_results$hex_block_summary
  
  combined_df <- rbind(hex150_df, hex200_df, hex300_df)
  
  write.csv(combined_df, file.path(results_dir, paste0("combined_hex_block_summary_",id,".csv")))
}

