# making robustness tables
library(dplyr)
library(msigdbr)
library(tidyr)

# paths
project_dir <- "."
results_dir <- file.path(project_dir, "results", "robustness")
table_dir <- file.path(project_dir, "results", "tables")

# load save output
res_rho <- read.csv(file.path(results_dir, "abl_correlation_threshold.csv"))
res_rho$X <- NULL

# ---- gene sets ( MSigDB KEGG legacy ) ----
msig <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
gene_sets <- split(msig$gene_symbol, msig$gs_name)

# summary retained pathways
rho_levels <- c(0, 0.05, 0.10, 0.20)

rho_sum <- res_rho %>%
  mutate(
    set_type = case_when(
      true_up == 25 & true_down == 25 ~ "Balanced (25/25)",
      true_up == 25 & true_down == 0  ~ "One-sided (25/0)",
      TRUE ~ NA_character_
    ),
    set_type = factor(
      set_type,
      levels = c("Balanced (25/25)", "One-sided (25/0)")
    ),
    rho = as.numeric(rho)
  ) %>%
  filter(
    !is.na(set_type),
    rho %in% rho_levels
  )

# pathway retention
retention_summary <- rho_sum %>%
  group_by(set_type, rho) %>%
  summarise(
    n_kept = n_distinct(pathway),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    set_type = c("Balanced (25/25)", "One-sided (25/0)"),
    rho = rho_levels,
    fill = list(n_kept = 0)
  ) %>%
  mutate(
    n_total = length(gene_sets),
    retention_rate = n_kept / n_total
  ) 

# AUROC among retained pathways
auroc_summary <- rho_sum %>%
  group_by(set_type, rho) %>%
  summarise(
    n_evaluable = sum(!is.na(AUROC_data)),
    mean_auroc = ifelse(
      n_evaluable > 0,
      mean(AUROC_data, na.rm = TRUE),
      NA_real_
    ),
    sd_auroc = ifelse(
      n_evaluable > 1,
      sd(AUROC_data, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(-n_evaluable)

# runtime
runtime_summary <- rho_sum %>%
  dplyr::select(set_type, rho, times_data_driven) %>%
  dplyr::distinct() %>%
  dplyr::group_by(set_type, rho) %>%
  dplyr::summarise(
    runtime_seconds = dplyr::first(times_data_driven),
    .groups = "drop"
  )

# combine to table s3
table_s3 <- retention_summary %>%
  left_join(
    auroc_summary,
    by = c("set_type", "rho")
  ) %>%
  left_join(
    runtime_summary,
    by = c("set_type", "rho")
  ) %>%
  arrange(set_type, rho)

write.csv(table_s3, 
          file.path(table_dir, "supplement_Table_S6_robustness_correlation_retains.csv"),
          row.names = FALSE)
