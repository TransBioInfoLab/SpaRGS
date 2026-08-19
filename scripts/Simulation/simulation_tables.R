# making simulation tables
library(dplyr)

# paths
project_dir <- "."
results_dir <- file.path(project_dir, "results", "simulation")
table_dir <- file.path(project_dir, "results", "tables")

# load save outputs
region_raw <- read.csv(file.path(results_dir, "evaluate_metrics_all_replicates_data_driven_plus_gsdensity.csv"))
region_raw$X <- NULL


# ------------------------------------------------------------
# Retained-matched simulation summary
# ------------------------------------------------------------
# Identify sGSVA-retained pathway-replicate-scenario keys from main benchmark
retained_keys <- region_raw %>%
  filter(methods == "sGSVA", retained == TRUE, !is.na(AUROC)) %>%
  distinct(pathway, replicate, ell, sigma_f, up_size, down_size)

region_retained_matched <- region_raw %>%
  semi_join(
    retained_keys,
    by = c("pathway", "replicate", "ell", "sigma_f", "up_size", "down_size")
  ) %>%
  mutate(
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0 ~ "One-sided (25/0)",
      TRUE ~ "Other"
    )
  )

# get comparator summary results
spatial_comp_summary <- region_retained_matched %>%
  group_by(methods) %>%
  summarise(
    n_evaluable = sum(!is.na(AUROC)),
    mean_auroc = mean(AUROC, na.rm = TRUE),
    sd_auroc = sd(AUROC, na.rm = TRUE),
    mean_abs_spearman = mean(spearman, na.rm = TRUE),
    sd_abs_spearman = sd(spearman, na.rm = TRUE),
    mean_signed_spearman = mean(spearman_signed, na.rm = TRUE),
    sd_signed_spearman = sd(spearman_signed, na.rm = TRUE),
    .groups = "drop"
  )

runtime_comp_summary <- region_raw %>%
  distinct(methods, replicate, ell, sigma_f, up_size, down_size, times) %>%
  group_by(methods) %>%
  summarise(
    mean_runtime = mean(times, na.rm = TRUE),
    sd_runtime = sd(times, na.rm = TRUE),
    .groups = "drop"
  )
main_simulation_summary <- spatial_comp_summary %>%
  left_join(runtime_comp_summary, by = "methods")

write.csv(
  main_simulation_summary,
  file.path(table_dir, "main_table1_simulation_retained_matched_summary.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# Paired bootstrap confidence intervals for retained-matched simulation benchmark
# ------------------------------------------------------------
bootstrap_paired_ci <- function(
    dat,
    metric,
    reference_method = "sGSVA",
    cluster_cols = c("replicate", "ell", "sigma_f", "up_size", "down_size"),
    B = 1000,
    seed = 1
) {
  set.seed(seed)
  
  dat <- dat %>% filter(!is.na(.data[[metric]]))
  
  dat$cluster_id <- do.call(
    interaction,
    c(dat[, cluster_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  )
  
  clusters <- unique(dat$cluster_id)
  
  observed_means <- dat %>%
    group_by(methods) %>%
    summarise(
      mean_metric = mean(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    )
  
  ref_observed <- observed_means %>%
    filter(methods == reference_method) %>%
    pull(mean_metric)
  
  observed_diff <- observed_means %>%
    filter(methods != reference_method) %>%
    transmute(
      comparator = methods,
      observed_reference = ref_observed,
      observed_comparator = mean_metric,
      observed_difference = ref_observed - mean_metric
    )
  
  boot_diffs <- purrr::map_dfr(seq_len(B), function(b) {
    sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
    
    boot_dat <- purrr::map_dfr(seq_along(sampled_clusters), function(i) {
      dat %>%
        filter(cluster_id == sampled_clusters[i]) %>%
        mutate(boot_cluster = i)
    })
    
    boot_means <- boot_dat %>%
      group_by(methods) %>%
      summarise(
        mean_metric = mean(.data[[metric]], na.rm = TRUE),
        .groups = "drop"
      )
    
    ref_mean <- boot_means %>%
      filter(methods == reference_method) %>%
      pull(mean_metric)
    
    boot_means %>%
      filter(methods != reference_method) %>%
      transmute(
        bootstrap_id = b,
        comparator = methods,
        difference = ref_mean - mean_metric
      )
  })
  
  ci <- boot_diffs %>%
    group_by(comparator) %>%
    summarise(
      ci_lower = quantile(difference, 0.025, na.rm = TRUE),
      ci_upper = quantile(difference, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  observed_diff %>%
    left_join(ci, by = "comparator") %>%
    mutate(metric = metric, .before = 1) %>%
    arrange(desc(observed_difference))
}

# retained_matched should be the 745-case matched benchmark dataset
boot_auc <- bootstrap_paired_ci(
  dat = region_retained_matched,
  metric = "AUROC",
  B = 1000,
  seed = 1
)

boot_abs_spearman <- bootstrap_paired_ci(
  dat = region_retained_matched,
  metric = "spearman",
  B = 1000,
  seed = 1
)

paired_bootstrap_summary <- bind_rows(
  boot_auc,
  boot_abs_spearman
)

write_csv(
  paired_bootstrap_summary,
  file.path(outdir, "supplement_Table_S2_simulation_paired_bootstrap_ci.csv")
)


# ------------------------------------------------------------
# Retained vs non-retained sGSVA cases
# ------------------------------------------------------------
region_full_query <- region_raw %>%
  mutate(
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0  ~ "One-sided (25/0)",
      TRUE ~ "Other"
    ),
    evaluable = retained == TRUE & !is.na(AUROC)
  )

key_cols <- c(
  "pathway", "replicate", "ell", "sigma_f", "up_size", "down_size"
)

sgsva_status <- region_full_query %>%
  filter(methods == "sGSVA") %>%
  distinct(
    across(all_of(key_cols)),
    set_type,
    sgsva_retained = evaluable
  )

# Retention rate by signal
sgsva_retention_by_signal <- sgsva_status %>%
  group_by(sigma_f) %>%
  summarise(
    n_queried = n(),
    n_retained = sum(sgsva_retained, na.rm = TRUE),
    n_nonretained = sum(!sgsva_retained, na.rm = TRUE),
    retention_rate = n_retained / n_queried,
    .groups = "drop"
  ) 

write_csv(
  sgsva_retention_by_signal,
  file.path(table_dir, "supplement_Table_S3A_simulation_sgsva_retention_by_signal.csv")
)

# Use other methods as a rough difficulty check on the same cases
competitor_case_difficulty <- region_full_query %>%
  filter(methods != "sGSVA") %>%
  group_by(across(all_of(key_cols))) %>%
  summarise(
    competitor_mean_auroc = mean(AUROC, na.rm = TRUE),
    competitor_mean_abs_spearman = mean(spearman, na.rm = TRUE),
    .groups = "drop"
  )

retained_vs_nonretained_summary <- sgsva_status %>%
  left_join(competitor_case_difficulty, by = key_cols) %>%
  group_by(sgsva_retained) %>%
  summarise(
    n_cases = n(),
    prop_low_signal = mean(sigma_f == 0.5, na.rm = TRUE),
    prop_moderate_signal = mean(sigma_f == 1.0, na.rm = TRUE),
    prop_high_signal = mean(sigma_f == 1.5, na.rm = TRUE),
    prop_balanced = mean(set_type == "Balanced (25/25)", na.rm = TRUE),
    prop_onesided= mean(set_type == "One-sided (25/0)", na.rm = TRUE),
    mean_competitor_auroc = mean(competitor_mean_auroc, na.rm = TRUE),
    mean_competitor_abs_spearman = mean(competitor_mean_abs_spearman, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  retained_vs_nonretained_summary,
  file.path(table_dir, "supplement_Table_S3B_simulation_sgsva_retained_vs_nonretained.csv")
)

# simulation signed-orientation summary
simulation_orientation_summary <- region_retained_matched %>%
  filter(methods == "sGSVA", retained == TRUE, !is.na(spearman_signed)) %>%
  group_by(set_type, sigma_f) %>%
  summarise(
    n_retained= n(),
    n_positive_signed_spearman = sum(spearman_signed > 0, na.rm = TRUE),
    prop_positive_signed_spearman = mean(spearman_signed > 0, na.rm = TRUE),
    mean_signed_spearman = mean(spearman_signed, na.rm = TRUE),
    median_signed_spearman = median(spearman_signed, na.rm = TRUE),
    mean_abs_spearman = mean(spearman, na.rm = TRUE),
    median_abs_spearman = median(spearman, na.rm = TRUE),
    .groups = "drop"
  ) 
write_csv(
  simulation_orientation_summary,
  file.path(table_dir, "supplement_Table_S3C_simulation_signed_orientation_summary.csv")
)


# ------------------------------------------------------------
# Full-query simulation summary
# ------------------------------------------------------------
# Purpose:
# This summarizes all intended pathway-replicate-scenario queries,
# not only the cases retained by data-driven sGSVA.
#
# The retained-matched table remains the primary paired accuracy benchmark.
# This full-query table is for transparency about evaluability/coverage.

# Overall full-query summary by method
full_query_summary <- region_full_query %>%
  group_by(methods) %>%
  summarise(
    n_queried = n(),
    n_evaluable = sum(evaluable, na.rm = TRUE),
    evaluable_rate = n_evaluable / n_queried,
    
    mean_auroc = mean(AUROC, na.rm = TRUE),
    sd_auroc = sd(AUROC, na.rm = TRUE),
    mean_abs_spearman = mean(spearman, na.rm = TRUE),
    sd_abs_spearman = sd(spearman, na.rm = TRUE),
    mean_signed_spearman = mean(spearman_signed, na.rm = TRUE),
    sd_signed_spearman = sd(spearman_signed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_auroc))


# ------------------------------------------------------------
# Conservative full-query sensitivity for data-driven sGSVA
# ------------------------------------------------------------
# Non-retained sGSVA cases are assigned:
# AUROC = 0.5 and absolute Spearman = 0.
# This is a sensitivity analysis, not the primary benchmark.

region_full_query_conservative <- region_full_query %>%
  mutate(
    AUROC_conservative = case_when(
      methods == "sGSVA" & !evaluable ~ 0.5,
      TRUE ~ AUROC
    ),
    spearman_conservative = case_when(
      methods == "sGSVA" & !evaluable ~ 0,
      TRUE ~ spearman
    ),
    spearman_signed_conservative = case_when(
      methods == "sGSVA" & !evaluable ~ 0,
      TRUE ~ spearman_signed
    )
  )

conservative_full_query_summary <- region_full_query_conservative %>%
  group_by(methods) %>%
  summarise(
    n_queried = n(),
    n_evaluable = sum(evaluable, na.rm = TRUE),
    evaluable_rate = n_evaluable / n_queried,
    
    mean_auroc_conservative = mean(AUROC_conservative, na.rm = TRUE),
    mean_abs_spearman_conservative = mean(spearman_conservative, na.rm = TRUE),
    mean_signed_spearman_conservative = mean(spearman_signed_conservative, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(desc(mean_auroc_conservative))

full_query_dat <- full_query_summary %>%
  left_join(conservative_full_query_summary %>%
              dplyr::select(methods, mean_auroc_conservative,mean_abs_spearman_conservative, mean_signed_spearman_conservative),
            by = "methods")
write.csv(full_query_dat,
          file.path(table_dir, "supplement_Table_S4_simulation_full_query_summary.csv"))



# ------------------------------------------------------------
# GSDensity parameter sensitivity check
# ------------------------------------------------------------
gsdensity_param_res <- read.csv(file.path(results_dir, "gsdensity_parameter_sensitivity_pathway_level.csv"))

# Overall sensitivity summary
gsdensity_param_summary <- gsdensity_param_res %>%
  group_by(setting_id, n_mca_dims, nn_use, restart) %>%
  summarise(
    n_queried = n(),
    n_evaluable = sum(evaluable, na.rm = TRUE),
    evaluable_rate = n_evaluable / n_queried,
    
    mean_auroc = mean(AUROC, na.rm = TRUE),
    sd_auroc = sd(AUROC, na.rm = TRUE),
    
    mean_abs_spearman = mean(spearman, na.rm = TRUE),
    sd_abs_spearman = sd(spearman, na.rm = TRUE),
    
    mean_signed_spearman = mean(spearman_signed, na.rm = TRUE),
    sd_signed_spearman = sd(spearman_signed, na.rm = TRUE),
    
    mean_runtime = mean(times, na.rm = TRUE),
    sd_runtime = sd(times, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(desc(mean_auroc))

write.csv(
  gsdensity_param_summary,
  file.path(table_dir, "supplement_Table_S5_simulation_gsdensity_parameter_sensitivity.csv"),
  row.names = FALSE
)

