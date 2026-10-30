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
    rho = as.numeric(rho),
    retained = as.logical(retained)
  ) %>%
  filter(
    !is.na(set_type),
    rho %in% rho_levels
  )

# pathway retention
retention_summary <- rho_sum %>%
  group_by(set_type, rho) %>%
  summarise(
   n_kept = sum(retained, na.rm=TRUE),
   n_total = n(),
   retention_rate = n_kept/n_total,
   .groups = "drop"
  ) 

# AUROC among retained pathways
auroc_summary <- rho_sum %>%
  filter(retained, !is.na(AUROC)) %>%
  group_by(set_type, rho) %>%
  summarise(
    mean_auroc = mean(AUROC, na.rm = TRUE),
    sd_auroc = sd(AUROC, na.rm = TRUE),
    .groups = "drop"
  )

# full-query penalized AUROC (no-retained AUROC=0.5)
penalized_summary <- rho_sum %>%
  mutate(
    AUROC_penalized = ifelse(retained & !is.na(AUROC), AUROC, 0.5)
  ) %>%
  group_by(set_type, rho) %>%
  summarise(
    mean_auroc_penalized = mean(AUROC_penalized),
    .groups = "drop"
  )
  
# runtime
runtime_summary <- rho_sum %>%
  dplyr::select(replicate, ell, sigma_f, set_type, rho, runtime_seconds) %>%
  distinct() %>%
  group_by(set_type, rho) %>%
  summarise(
    runtime_seconds = mean(runtime_seconds, na.rm = TRUE),
    n_runtime_runs = n(),
    .groups = "drop"
  )

# Combine main robustness table
table_s6 <- retention_summary %>%
  left_join(
    auroc_summary,
    by = c("set_type", "rho")
  ) %>%
  left_join(
    penalized_summary,
    by = c("set_type", "rho")
  ) %>%
  left_join(
    runtime_summary,
    by = c("set_type", "rho")
  ) %>%
  arrange(set_type, rho)

write.csv(table_s6, 
          file.path(table_dir, "supplement_Table_S6_robustness_correlation_retains.csv"),
          row.names = FALSE)
