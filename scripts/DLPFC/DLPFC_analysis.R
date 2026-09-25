### DLPFC analysis for final manuscript
### Section-level / donor-aware summaries; no Stouffer meta-analysis

suppressPackageStartupMessages({
  library(spatialLIBD)
  library(SpatialExperiment)
  library(scater)
  library(scran)
  library(ggplot2)
  library(viridis)
  library(spdep)
  library(msigdbr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(FNN)
  library(GSVA)
  library(org.Hs.eg.db)
  library(here)
})

# ---------------------------
# Paths and reproducibility
# ---------------------------
project_dir <- "."
results_dir <- file.path(project_dir, "results", "DLPFC")
supp_dir <- file.path(project_dir, "results", "supplementary_data")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_dir, "R", "sgsva_score.R"))
source(file.path(project_dir, "R", "spatial_permutation.R"))

set.seed(1)

# ---------------------------
# Load data and gene sets
# ---------------------------
spe <- fetch_data(type = "spe")
sample_ids <- unique(as.character(spe$sample_id))
spe_list <- lapply(sample_ids, function(id) spe[, as.character(spe$sample_id) == id])
names(spe_list) <- sample_ids

kegg_gene_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcollection = "CP:KEGG_LEGACY"
)
gene_sets <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)

# ---------------------------
# Section-to-donor mapping
# ---------------------------
# Prefer donor/subject metadata in spatialLIBD when available. The fallback
# reflects the known grouping of the 12 sections into three donors.

infer_donor_map <- function(spe, sample_ids) {
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  if (!"sample_id" %in% colnames(cd)) {
    cd$sample_id <- as.character(spe$sample_id)
  }
  cd$sample_id <- as.character(cd$sample_id)

  candidate_cols <- intersect(
    c("subject", "subject_id", "donor", "donor_id", "brnum"),
    colnames(cd)
  )

  for (nm in candidate_cols) {
    tmp <- cd %>%
      dplyr::select(sample_id, donor_id = all_of(nm)) %>%
      dplyr::filter(!is.na(donor_id), donor_id != "") %>%
      dplyr::distinct()

    n_per_section <- tmp %>% dplyr::count(sample_id)
    if (all(sample_ids %in% tmp$sample_id) && all(n_per_section$n == 1)) {
      out <- setNames(as.character(tmp$donor_id), tmp$sample_id)
      return(out[sample_ids])
    }
  }

  fallback <- c(
    "151507" = "Donor1", "151508" = "Donor1",
    "151509" = "Donor1", "151510" = "Donor1",
    "151669" = "Donor2", "151670" = "Donor2",
    "151671" = "Donor2", "151672" = "Donor2",
    "151673" = "Donor3", "151674" = "Donor3",
    "151675" = "Donor3", "151676" = "Donor3"
  )

  if (!all(sample_ids %in% names(fallback))) {
    stop("Could not infer donor mapping from spatialLIBD metadata and fallback mapping is incomplete.")
  }
  fallback[sample_ids]
}

donor_by_section <- infer_donor_map(spe, sample_ids)

write.csv(
  data.frame(
    sample_id = sample_ids,
    donor_id = unname(donor_by_section[sample_ids])
  ),
  file.path(results_dir, "DLPFC_section_donor_map.csv"),
  row.names = FALSE
)

# ---------------------------
# Map Ensembl IDs to symbols
# ---------------------------
all_genes <- rownames(counts(spe))
ensembl_ids <- gsub("\\..*", "", all_genes)
symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
map_table <- tibble(ensembl = ensembl_ids, symbol = symbols) %>%
  filter(!is.na(symbol) & symbol != "") %>%
  distinct(ensembl, .keep_all = TRUE)

# Preserve the manuscript analysis behavior: keep the first row when multiple
# Ensembl IDs map to the same HGNC symbol.
convert_rows_to_symbols <- function(expr, map_table) {
  rn <- gsub("\\..*", "", rownames(expr))
  matched <- map_table$symbol[match(rn, map_table$ensembl)]
  keep_map <- !is.na(matched) & matched != ""
  expr <- expr[keep_map, , drop = FALSE]
  matched <- matched[keep_map]
  rownames(expr) <- matched
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  expr
}

# ---------------------------
# Pixel-to-micron conversion
# ---------------------------
estimate_um_per_px_from_grid <- function(coords_px, spot_spacing_um = 100) {
  nn <- FNN::get.knn(as.matrix(coords_px[, c("x", "y")]), k = 1)$nn.dist[, 1]
  median_nn_px <- median(nn[is.finite(nn) & nn > 0], na.rm = TRUE)

  data.frame(
    median_nn_px = median_nn_px,
    spot_spacing_um_used = spot_spacing_um,
    um_per_px = spot_spacing_um / median_nn_px
  )
}

preprocess_per_sample <- function(spe_sub, map_table) {
  counts_mat <- as.matrix(assay(spe_sub, "counts"))
  counts_mat <- convert_rows_to_symbols(counts_mat, map_table)

  keep <- rowSums(counts_mat > 0) >= floor(0.05 * ncol(counts_mat))
  expr <- counts_mat[keep, , drop = FALSE]

  lib <- colSums(expr)
  expr_logCPM <- log1p(t(t(expr) / pmax(lib, 1)) * 1e6)

  coords_px <- as.data.frame(spatialCoords(spe_sub))[, c(
    "pxl_col_in_fullres",
    "pxl_row_in_fullres"
  ), drop = FALSE]
  colnames(coords_px) <- c("x", "y")

  scale_info <- estimate_um_per_px_from_grid(coords_px)
  coords_um <- coords_px * scale_info$um_per_px

  # Explicitly align coordinate rows to expression columns.
  stopifnot(all(colnames(expr_logCPM) %in% rownames(coords_um)))
  coords_px <- coords_px[colnames(expr_logCPM), , drop = FALSE]
  coords_um <- coords_um[colnames(expr_logCPM), , drop = FALSE]

  list(
    sample_id = unique(as.character(spe_sub$sample_id)),
    expr = expr_logCPM,
    coords_px = coords_px,
    coords_um = coords_um,
    um_per_px = scale_info$um_per_px,
    median_nn_px = scale_info$median_nn_px,
    counts_mat = counts_mat
  )
}

dlpfc_preprocessed <- purrr::map(
  spe_list,
  preprocess_per_sample,
  map_table = map_table
)
saveRDS(dlpfc_preprocessed, file.path(results_dir, "DLPFC_preprocess.rds"))

# ---------------------------
# Score GSVA and sGSVA
# ---------------------------
sigma_table_um <- sapply(sample_ids, function(id) {
  estimate_sigma_knn(dlpfc_preprocessed[[id]]$coords_um)
})

write.csv(
  data.frame(
    sample_id = names(sigma_table_um),
    donor_id = unname(donor_by_section[names(sigma_table_um)]),
    sigma_knn_um = as.numeric(sigma_table_um),
    sigma_used_um = 151,
    um_per_px = sapply(dlpfc_preprocessed, `[[`, "um_per_px"),
    median_nn_px = sapply(dlpfc_preprocessed, `[[`, "median_nn_px")
  ),
  file.path(results_dir, "DLPFC_coordinate_scale_summary.csv"),
  row.names = FALSE
)

full_res <- list()
for (id in sample_ids) {
  message("Scoring section ", id, " ...")

  t0 <- proc.time()[3]
  params <- gsvaParam(dlpfc_preprocessed[[id]]$expr, gene_sets)
  gsva_scores <- gsva(params, verbose = FALSE)
  rt_gsva <- proc.time()[3] - t0

  t0 <- proc.time()[3]
  sgsva_res <- run_sgsva(
    dlpfc_preprocessed[[id]]$expr,
    dlpfc_preprocessed[[id]]$coords_um,
    gene_sets,
    sigma_global = 151,
    cor_min = 0.2,
    mode = "delaunay",
    one_sided = "center_0.5"
  )
  rt_sgsva <- proc.time()[3] - t0

  full_res[[id]] <- list(
    gsva_scores = gsva_scores,
    sgsva_res = sgsva_res,
    runtime = list(sgsva = rt_sgsva, gsva = rt_gsva)
  )
}
saveRDS(full_res, file.path(results_dir, "full_results_12_samples.rds"))

run_summary <- bind_rows(lapply(sample_ids, function(id) {
  data.frame(
    sample_id = id,
    donor_id = unname(donor_by_section[id]),
    n_spots = ncol(dlpfc_preprocessed[[id]]$expr),
    n_sgsva_pathways = length(full_res[[id]]$sgsva_res$gene_sets),
    runtime_sgsva_s = as.numeric(full_res[[id]]$runtime$sgsva),
    runtime_gsva_s = as.numeric(full_res[[id]]$runtime$gsva)
  )
}))
write.csv(run_summary, file.path(results_dir, "DLPFC_run_summary.csv"), row.names = FALSE)

# ---------------------------
# Spatial permutation analyses
# ---------------------------
for (id in rev(sample_ids)) {
  message("Spatial inference for section ", id, " ...")

  sgsva_res <- full_res[[id]]$sgsva_res
  z_sgsva <- sgsva_res$z_score
  balanced_sets <- sgsva_res$gene_sets

  coords <- dlpfc_preprocessed[[id]]$coords_um
  colnames(coords)[1:2] <- c("x", "y")
  
  t0 <- proc.time()[3]
  global_res <- spatial_permutation(
    z_sgsva, coords, balanced_sets,
    mode = "delaunay",
    sigma_global = 151,
    cor_min = 0.2,
    one_sided = "center_0.5",
    perm_methods = c("moran_mc", "geary_mc")
  )
  rt <- proc.time()[3]-t0
  print(rt)
  saveRDS(global_res, file.path(results_dir, paste0(id, "_moranmc_res.rds")))

  block_sizes <- c(150, 200, 300)
  block_res <- lapply(block_sizes, function(b) {
    
    message(b, "hex-block permutation")
    
    t0 <- proc.time()[3]
    spatial_permutation(
      z_sgsva, coords, balanced_sets,
      mode = "delaunay",
      sigma_global = 151,
      cor_min = 0.2,
      one_sided = "center_0.5",
      perm_methods = "hex_block",
      hex_cellsize = b
    )
    rt <- proc.time()[3]-t0
    print(rt)
  })
  names(block_res) <- paste0("block", block_sizes)
  saveRDS(block_res, file.path(results_dir, paste0(id, "_block_res.rds")))
}

# -------------------------------------------------------------------------
# FINAL DLPFC CROSS-SECTION SUMMARIES
# No Stouffer meta-analysis: summarize section recurrence, donor coverage,
# and section-level effect sizes.
# -------------------------------------------------------------------------
standardize_spatial_summary <- function(
    df,
    sample_id,
    donor_id,
    permutation_scheme,
    block_size_um = NA_real_) {
  
  out <- df
  out$sample_id = sample_id
  out$donor_id = donor_id
  out$permutation_scheme = permutation_scheme
  out$block_size_um = block_size_um
  out$moran_fdr <- p.adjust(out$p_moran, method = "BH")
  out$geary_fdr <- p.adjust(out$p_geary, method = "BH")
  out$one_sided <- "center_0.5"
  out$seed_permutation <- 124
  out %>%
    dplyr::select(-obs_geary_flipped, -z_geary_flipped,
                  -null_mean_geary_flipped, -null_sd_geary_flipped)
}

extract_global_summary <- function(x) {
  if (!is.null(x$permutation_results$moran_geary_results$moran_geary_summary)) {
    return(x$permutation_results$moran_geary_results$moran_geary_summary)
  }
  if (!is.null(x$moran_geary_results$moran_geary_summary)) {
    return(x$moran_geary_results$moran_geary_summary)
  }
  stop("Could not locate moran_geary_summary in global spatial_permutation output.")
}

extract_hex_summary <- function(x) {
  if (!is.null(x$permutation_results$hex_block_results$hex_block_summary)) {
    return(x$permutation_results$hex_block_results$hex_block_summary)
  }
  if (!is.null(x$hex_block_results$hex_block_summary)) {
    return(x$hex_block_results$hex_block_summary)
  }
  stop("Could not locate hex_block_summary in block spatial_permutation output.")
}

global_section_list <- list()
hex_section_list <- list()
block_sizes <- c(150, 200, 300)

for (id in sample_ids) {
  donor_id <- unname(donor_by_section[id])

  global_res <- readRDS(file.path(results_dir, paste0(id, "_moranmc_res.rds")))
  global_df <- extract_global_summary(global_res)
  global_section_list[[id]] <- standardize_spatial_summary(
    global_df,
    sample_id = id,
    donor_id = donor_id,
    permutation_scheme = "global",
    block_size_um = NA_real_
  )

  block_res <- readRDS(file.path(results_dir, paste0(id, "_block_res.rds")))
  for (j in seq_along(block_sizes)) {
    b <- block_sizes[j]
    hex_df <- extract_hex_summary(block_res[[j]])
    hex_section_list[[paste(id, b, sep = "_")]] <- standardize_spatial_summary(
      hex_df,
      sample_id = id,
      donor_id = donor_id,
      permutation_scheme = "hex_block",
      block_size_um = b
    )
  }
}

global_section_stats <- bind_rows(global_section_list)
hex_section_stats <- bind_rows(hex_section_list)

write.csv(
  global_section_stats,
  file.path(results_dir, "DLPFC_section_level_global_mc.csv"),
  row.names = FALSE
)
write.csv(
  hex_section_stats,
  file.path(results_dir, "DLPFC_section_level_hex_block.csv"),
  row.names = FALSE
)

# Per-section combined hex files are retained as convenient downstream inputs.
for (id in sample_ids) {
  tmp <- hex_section_stats %>% filter(sample_id == id)
  write.csv(
    tmp,
    file.path(results_dir, paste0("combined_hex_block_summary_", id, ".csv")),
    row.names = FALSE
  )
}

summarise_recurrence <- function(df, alpha = 0.05) {
  df %>%
    group_by(block_size_um, pathway) %>%
    summarise(
      n_sections_tested = n_distinct(sample_id[is.finite(obs_moran)]),
      n_sections_sig_moran = n_distinct(sample_id[!is.na(moran_fdr) & moran_fdr < alpha]),
      n_sections_sig_geary = n_distinct(sample_id[!is.na(geary_fdr) & geary_fdr < alpha]),
      n_donors_tested = n_distinct(donor_id[is.finite(obs_moran)]),
      n_donors_sig_moran = n_distinct(donor_id[!is.na(moran_fdr) & moran_fdr < alpha]),
      n_donors_sig_geary = n_distinct(donor_id[!is.na(geary_fdr) & geary_fdr < alpha]),
      median_moran_I = median(obs_moran, na.rm = TRUE),
      q1_moran_I = quantile(obs_moran, 0.25, na.rm = TRUE, names = FALSE),
      q3_moran_I = quantile(obs_moran, 0.75, na.rm = TRUE, names = FALSE),
      median_geary_C = median(obs_geary, na.rm = TRUE),
      q1_geary_C = quantile(obs_geary, 0.25, na.rm = TRUE, names = FALSE),
      q3_geary_C = quantile(obs_geary, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}

hex_recurrence <- summarise_recurrence(hex_section_stats)
write.csv(
  hex_recurrence,
  file.path(results_dir, "DLPFC_hex_block_recurrence_summary.csv"),
  row.names = FALSE
)

for (b in block_sizes) {
  write.csv(
    hex_recurrence %>% filter(block_size_um == b),
    file.path(results_dir, paste0("DLPFC_hex", b, "_recurrence_summary.csv")),
    row.names = FALSE
  )
}

# Donor-stratified summaries make the within-donor structure explicit.
donor_recurrence <- hex_section_stats %>%
  group_by(block_size_um, pathway, donor_id) %>%
  summarise(
    n_sections_tested = n_distinct(sample_id[is.finite(obs_moran)]),
    n_sections_sig_moran = n_distinct(sample_id[!is.na(moran_fdr) & moran_fdr < 0.05]),
    n_sections_sig_geary = n_distinct(sample_id[!is.na(geary_fdr) & geary_fdr < 0.05]),
    median_moran_I = median(obs_moran, na.rm = TRUE),
    median_geary_C = median(obs_geary, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(
  donor_recurrence,
  file.path(results_dir, "DLPFC_donor_stratified_recurrence_summary.csv"),
  row.names = FALSE
)

# Main Fig. 5G source: rank first by section recurrence, then donor coverage,
# then median section-level Moran's I. This does not pool p-values.
top20_hex200 <- hex_recurrence %>%
  filter(block_size_um == 200) %>%
  arrange(
    desc(n_sections_sig_moran),
    desc(n_donors_sig_moran),
    desc(median_moran_I),
    median_geary_C
  ) %>%
  slice_head(n = 20)

write.csv(
  top20_hex200,
  file.path(results_dir, "DLPFC_top20_recurrent_pathways_hex200.csv"),
  row.names = FALSE
)


library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "SData4_global_statistics")
writeData(wb, "SData4_global_statistics", global_section_stats)
addWorksheet(wb, "SData4_hexblock_statistics")
writeData(wb, "SData4_hexblock_statistics", hex_section_stats)
addWorksheet(wb, "SData4_hex_recurrence")
writeData(wb, "SData4_hex_recurrence", hex_recurrence)
addWorksheet(wb, "SData4_donor_recurrence")
writeData(wb, "SData4_donor_recurrence", donor_recurrence)
addWorksheet(wb, "SData4_runtime")
writeData(wb, "SData4_runtime", run_summary)
saveWorkbook(wb, file.path(supp_dir, "Supplementary Data 4.xlsx"), overwrite = TRUE)





