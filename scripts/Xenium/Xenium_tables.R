# making Xenium tables
library(dplyr)
library(hexbin)

# paths
project_dir <- "."
results_dir <- file.path(project_dir, "results", "xenium")
table_dir <- file.path(project_dir, "results", "tables")

source(file.path(project_dir,"R/sgsva_score.R"))

# generate spatial autocorrelation summary
moran_geary_all <- read.csv(file.path(results_dir, "moran_geary_three_methods.csv"))
sum_moran_geary <- moran_geary_all %>%
  group_by(method) %>%
  summarise(
    n_pathways = n(),
    median_moran = median(moran_I, na.rm = TRUE),
    q1_moran = quantile(moran_I, 0.25, na.rm = TRUE),
    q3_moran = quantile(moran_I, 0.75, na.rm = TRUE),
    median_geary = median(geary_C, na.rm = TRUE),
    q1_geary = quantile(geary_C, 0.25, na.rm = TRUE),
    q3_geary = quantile(geary_C, 0.75, na.rm = TRUE),
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
write.csv(sum_moran_geary, file.path(table_dir, "supplement_Table_S12A_xenium_median_autocorrelation_summary.csv"))

bootstrap_paired_spatial_diff <- function(
    dat,
    comparator,
    reference = "sGSVA",
    B = 5000,
    seed = 1
) {
  set.seed(seed)
  
  wide <- dat %>%
    filter(method %in% c(reference, comparator)) %>%
    dplyr::select(pathway, method, moran_I, geary_C) %>%
    tidyr::pivot_wider(
      names_from = method,
      values_from = c(moran_I, geary_C)
    )
  
  moran_ref <- paste0("moran_I_", reference)
  moran_cmp <- paste0("moran_I_", comparator)
  geary_ref <- paste0("geary_C_", reference)
  geary_cmp <- paste0("geary_C_", comparator)
  
  wide <- wide %>%
    filter(
      is.finite(.data[[moran_ref]]),
      is.finite(.data[[moran_cmp]]),
      is.finite(.data[[geary_ref]]),
      is.finite(.data[[geary_cmp]])
    ) %>%
    mutate(
      delta_moran = .data[[moran_ref]] - .data[[moran_cmp]],
      delta_geary = .data[[geary_cmp]] - .data[[geary_ref]]
    )
  
  n <- nrow(wide)
  
  boot <- purrr::map_dfr(seq_len(B), function(b) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    
    tibble::tibble(
      bootstrap_id = b,
      delta_moran = median(wide$delta_moran[idx], na.rm = TRUE),
      delta_geary = median(wide$delta_geary[idx], na.rm = TRUE)
    )
  })
  
  tibble::tibble(
    comparator = comparator,
    n_pathways = n,
    median_delta_moran = median(wide$delta_moran, na.rm = TRUE),
    moran_ci_lower = quantile(boot$delta_moran, 0.025, na.rm = TRUE),
    moran_ci_upper = quantile(boot$delta_moran, 0.975, na.rm = TRUE),
    median_delta_geary = median(wide$delta_geary, na.rm = TRUE),
    geary_ci_lower = quantile(boot$delta_geary, 0.025, na.rm = TRUE),
    geary_ci_upper = quantile(boot$delta_geary, 0.975, na.rm = TRUE)
  )
}

xenium_spatial_bootstrap <- bind_rows(
  bootstrap_paired_spatial_diff(moran_geary_all, "GSVA"),
  bootstrap_paired_spatial_diff(moran_geary_all, "GSVA_smoothed")
)

write.csv(
  xenium_spatial_bootstrap,
  file.path(table_dir, "supplement_Table_S12B_xenium_paired_moran_geary_bootstrap.csv"),
  row.names = FALSE
)


# Xenium pathway-module gene overlap and non-overlap concordance sensitivity
sensitivity_long <- read.csv(file.path(results_dir, "xenium_representative_nonoverlap_sensitivity_long.csv"))
overlap_df <- read.csv(file.path(results_dir, "xenium_representative_pathway_module_gene_overlap.csv"))
# wide summary table
sensitivity_wide <- sensitivity_long %>%
  tidyr::pivot_wider(
    id_cols = c(anchor_label, pathway, module, method, n_bins),
    names_from = module_version,
    values_from = c(rho, abs_rho),
    names_glue = "{.value}_{module_version}"
  ) %>%
  dplyr::left_join(
    overlap_df %>%
      dplyr::select(
        anchor_label,
        n_pathway_genes,
        n_module_genes,
        n_overlap_genes,
        n_nonoverlap_module_genes,
        overlap_fraction_pathway,
        overlap_fraction_module,
        overlap_genes
      ),
    by = "anchor_label"
  ) %>%
  dplyr::mutate(
    delta_abs_rho_nonoverlap_minus_original =
      `abs_rho_Non-overlap module` - `abs_rho_Original module`
  )

write.csv(
  sensitivity_wide,
  file.path(table_dir, "supplement_Table_S13_xenium_gene_overlap_and_nonoverlap_correlation.csv"),
  row.names = FALSE
)






