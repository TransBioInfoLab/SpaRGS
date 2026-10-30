# making breast Visium table
library(dplyr)
library(readr)
library(purrr)

# paths
project_dir <- "."
results_dir <- file.path(project_dir, "results", "breast_visium")
table_dir <- file.path(project_dir, "results", "tables")
pilot_dir <- file.path(results_dir, "enrichmap_pilot")
anchoring_dir <- file.path(results_dir, "celltype_anchoring")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# --------------------------------------
# load outputs
# --------------------------------------
# original data
bc_obj <- readRDS(file.path(results_dir, "bc_obj.rds"))

visium_data <- readRDS(file.path(results_dir, "breast_visium_sct_data.rds"))
expr_filtered <- visium_data$expr_filtered
coords_um <- visium_data$coords_um

# read other pathway comparators
spa_res <- readRDS(file.path(results_dir, "spargs_score_res.rds"))
spa_scores <- spa_res$scores
balanced_sets <- spa_res$gene_sets
gsva_scores <- read.csv(file.path(results_dir, "gsva_scores.csv"),
                        row.names = 1,
                        check.names = FALSE)
gsva_smoothed_scores <- read.csv(file.path(results_dir, "gsva_smoothed_scores.csv"),
                                 row.names = 1,
                                 check.names = FALSE)
gsdensity_scores <- read.csv(file.path(results_dir, "gsdensity_scores.csv"),
                             row.names = 1,
                             check.names = FALSE)
# evaluate breast visium EnrichMap output
em_default <- read.csv(
  file.path(pilot_dir, "breast_visium_enrichmap_default_corrected_all_scores.csv"),
  row.names = 1,
  check.names = FALSE
) 
em_smooth <- read.csv(
  file.path(pilot_dir, "breast_visium_enrichmap_smooth_only_all_scores.csv"),
  row.names = 1,
  check.names = FALSE
)

# runtime
runtime_four_results <- read.csv(file.path(results_dir, "pathway_methods_runtime.csv"))
runtime_enrichmap <- read.csv(file.path(pilot_dir, "breast_visium_enrichmap_all_runtime_summary.csv"))

# spatial statistics
spatial_stats <- read.csv(file.path(results_dir, "pathway_level_moran_geary_shared_pathways.csv"))
spatial_stats_all <- read.csv(file.path(pilot_dir, "breast_visium_enrichmap_pilot_pathway_spatial_stats.csv"))

# cell-type anchoring
six_method_anchor_clean <- read.csv(file.path(results_dir, "enrichmap_pilot", "supplement_exploratory_six_method_anchoring.csv"))

# prespecified biological anchors
prespecified_results <- read.csv(file.path(anchoring_dir, "anchoring_prespecified_pairs.csv"))


# --------------------------------------
# Generate Main Table
# --------------------------------------
# Main Moran/Geary summary: median and IQR
moran_geary_summary_iqr <- spatial_stats %>%
  group_by(method) %>%
  summarise(
    n_pathways = n(),
    median_moran = median(moran_i, na.rm = TRUE),
    q1_moran = quantile(moran_i, 0.25, na.rm = TRUE),
    q3_moran = quantile(moran_i, 0.75, na.rm = TRUE),
    median_geary = median(geary_c, na.rm = TRUE),
    q1_geary = quantile(geary_c, 0.25, na.rm = TRUE),
    q3_geary = quantile(geary_c, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    moran_display = sprintf(
      "%.3f (%.3f-%.3f)",
      median_moran, q1_moran, q3_moran
    ),
    geary_display = sprintf(
      "%.3f (%.3f-%.3f)",
      median_geary, q1_geary, q3_geary
    )
  )

# save outputs
readr::write_csv(
  moran_geary_summary_iqr,
  file.path(table_dir, "main_table2_breast_moran_geary_summary.csv")
)


# --------------------------------------
# Supplementary tables
# --------------------------------------
normalize_pathway_names <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  sub("^KEGG_", "", x)
}

# runtime table
runtime_enrichmap <- runtime_enrichmap %>%
  dplyr::mutate(
    method = ifelse(mode=="default_corrected", "EnrichMap-default", "EnrichMap-smooth"),
    runtime = runtime_seconds) %>%
  dplyr::select(method, runtime)
runtime_tab <- rbind(runtime_four_results, runtime_enrichmap)

# get gene-set coverage table
exploratory_score_list <- list(
  SpaRGS = as.matrix(spa_scores),
  GSVA = as.matrix(gsva_scores),
  `GSVA-smoothed` = as.matrix(gsva_smoothed_scores),
  GSDensity = as.matrix(gsdensity_scores),
  `EnrichMap-default` = as.matrix(em_default),
  `EnrichMap-smooth` = as.matrix(em_smooth)
)

exploratory_score_list <- lapply(exploratory_score_list, function(mat) {
  rownames(mat) <- normalize_pathway_names(rownames(mat))
  mat
})
exploratory_common_spots <- Reduce(
  intersect,
  c(lapply(exploratory_score_list, colnames), list(rownames(coords_um)))
)
exploratory_shared_pathways <- Reduce(
  intersect,
  lapply(exploratory_score_list, rownames)
)

enrichmap_exploratory_coverage <- purrr::imap_dfr(
  exploratory_score_list,
  function(mat, method_name) {
    tibble(
      method = method_name,
      pathways_available = nrow(mat),
      spots_available = ncol(mat),
      shared_pathways_for_exploratory_comparison = length(exploratory_shared_pathways),
      common_spots_for_exploratory_comparison = length(exploratory_common_spots)
    )
  }
)

coverage_runtime_tab <- enrichmap_exploratory_coverage %>%
  left_join(runtime_tab, by="method")
write.csv(coverage_runtime_tab, 
          file.path(table_dir, "supplement_Table_S7_breast_coverage_runtime.csv"))


# spargs orientation diagnostics
kegg_gene_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcollection = "CP:KEGG_LEGACY"
)
gene_sets <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)
pathways <- c("KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
              "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
              "KEGG_GLUTATHIONE_METABOLISM",
              "KEGG_FATTY_ACID_METABOLISM")

check_final_orientation <- function(
    pathway,
    spa_scores,
    expr_matrix,
    gene_sets
) {
  genes <- intersect(gene_sets[[pathway]], rownames(expr_matrix))
  
  pathway_mean <- colMeans(expr_matrix[genes, , drop = FALSE],na.rm = TRUE)
  
  common_spots <- intersect(colnames(spa_scores), names(pathway_mean))
  
  cor(
    as.numeric(spa_scores[pathway, common_spots]),
    pathway_mean[common_spots],
    method = "spearman",
    use = "complete.obs"
  )
}

res_orien <- sapply(pathways, function(p) 
  check_final_orientation(p, spa_scores, expr_filtered, gene_sets))

diag_tab <- data.frame(
  pathway = gsub("KEGG_", "", names(res_orien)), 
  cor_expr = as.numeric(res_orien),
  signed_rho = NA_real_, abs_rho = NA_real_,
  pos_genes = NA_real_, neg_genes = NA_real_
)
for (i in seq_len(length(res_orien))) {
  pathway <- diag_tab$pathway[i]
  diag_tab$signed_rho[i] <- prespecified_results$rho[prespecified_results$method == "SpaRGS" 
                                                     & prespecified_results$pathway == pathway]
  diag_tab$abs_rho[i] <- prespecified_results$abs_rho[prespecified_results$method == "SpaRGS" 
                                                      & prespecified_results$pathway == pathway]
  diag_tab$pos_genes[i] <- length(balanced_sets[[names(res_orien)[i]]]$U)
  diag_tab$neg_genes[i] <- length(balanced_sets[[names(res_orien)[i]]]$D)
}

write_csv(diag_tab, 
          file.path(table_dir, "supplement_Table_S8_breast_spargs_orientation_diagnostics.csv"))


# Paired bootstrap CIs for Moran/Geary differences
bootstrap_paired_spatial_diff <- function(
    dat,
    comparator,
    reference = "SpaRGS",
    B = 5000,
    seed = 1
) {
  set.seed(seed)
  
  wide <- dat %>%
    filter(method %in% c(reference, comparator)) %>%
    dplyr::select(pathway, method, moran_i, geary_c) %>%
    tidyr::pivot_wider(
      names_from = method,
      values_from = c(moran_i, geary_c)
    )
  
  moran_ref <- paste0("moran_i_", reference)
  moran_cmp <- paste0("moran_i_", comparator)
  geary_ref <- paste0("geary_c_", reference)
  geary_cmp <- paste0("geary_c_", comparator)
  
  wide <- wide %>%
    filter(
      is.finite(.data[[moran_ref]]),
      is.finite(.data[[moran_cmp]]),
      is.finite(.data[[geary_ref]]),
      is.finite(.data[[geary_cmp]])
    ) %>%
    mutate(
      # Positive values favor SpaRGS
      delta_moran = .data[[moran_ref]] - .data[[moran_cmp]],
      delta_geary = .data[[geary_cmp]] - .data[[geary_ref]]
    )
  
  n <- nrow(wide)
  
  obs <- tibble(
    comparator = comparator,
    n_pathways = n,
    median_delta_moran = median(wide$delta_moran, na.rm = TRUE),
    median_delta_geary = median(wide$delta_geary, na.rm = TRUE)
  )
  
  boot <- purrr::map_dfr(seq_len(B), function(b) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    
    tibble(
      bootstrap_id = b,
      delta_moran = median(wide$delta_moran[idx], na.rm = TRUE),
      delta_geary = median(wide$delta_geary[idx], na.rm = TRUE)
    )
  })
  
  ci <- boot %>%
    summarise(
      moran_ci_lower = quantile(delta_moran, 0.025, na.rm = TRUE),
      moran_ci_upper = quantile(delta_moran, 0.975, na.rm = TRUE),
      geary_ci_lower = quantile(delta_geary, 0.025, na.rm = TRUE),
      geary_ci_upper = quantile(delta_geary, 0.975, na.rm = TRUE)
    )
  
  bind_cols(obs, ci)
}

spatial_comparators <- c("GSVA", "GSVA_smoothed", "GSDensity")
paired_moran_geary_bootstrap <- purrr::map_dfr(
  spatial_comparators,
  ~ bootstrap_paired_spatial_diff(
    dat = spatial_stats,
    comparator = .x,
    reference = "SpaRGS",
    B = 5000,
    seed = 1
  )
)

write_csv(
  paired_moran_geary_bootstrap,
  file.path(table_dir, "supplement_Table_S9_breast_paired_moran_geary_bootstrap_ci.csv")
)


# Exploratory EnrichMap Moran/Geary summary 
moran_geary_summary_iqr <- spatial_stats_all %>%
  group_by(method) %>%
  summarise(
    n_pathways = n(),
    median_moran = median(moran_i, na.rm = TRUE),
    q1_moran = quantile(moran_i, 0.25, na.rm = TRUE),
    q3_moran = quantile(moran_i, 0.75, na.rm = TRUE),
    median_geary = median(geary_c, na.rm = TRUE),
    q1_geary = quantile(geary_c, 0.25, na.rm = TRUE),
    q3_geary = quantile(geary_c, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    moran_display = sprintf(
      "%.3f (%.3f-%.3f)",
      median_moran, q1_moran, q3_moran
    ),
    geary_display = sprintf(
      "%.3f (%.3f-%.3f)",
      median_geary, q1_geary, q3_geary
    )
  )
write_csv(
  moran_geary_summary_iqr,
  file.path(table_dir, "supplement_Table_S10_breast_exploratory_enrichmap_spatial_summary.csv")
)


# Exploratory six-method pathway-cell-type spatial concordance
# Wide table for report/supplement
six_method_anchor_wide <- six_method_anchor_clean %>%
  dplyr::select(method, anchor_label_display, signed_rho) %>%
  tidyr::pivot_wider(
    names_from = anchor_label_display,
    values_from = signed_rho
  )

write_csv(
  six_method_anchor_wide,
  file.path(table_dir, "supplement_Table_S11_breast_exploratory_cell_type_concordance.csv")
)





