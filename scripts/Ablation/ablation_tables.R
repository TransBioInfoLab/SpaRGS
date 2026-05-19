# making ablation tables
library(dplyr)
library(msigdbr)

# paths
project_dir <- "."
results_dir <- file.path(project_dir, "results", "ablation")
table_dir <- file.path(project_dir, "results", "tables")

# load save output
res_rho <- read.csv(file.path(results_dir, "abl_correlation_threshold.csv"))
res_rho$X <- NULL

# ---- gene sets ( MSigDB KEGG legacy ) ----
msig <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
gene_sets <- split(msig$gene_symbol, msig$gs_name)

# summary retained pathways
rho_sum <- res_rho %>%
  mutate(
    set_type = case_when(
      true_up==25 & true_down==25 ~ "Balanced (25/25)",
      true_up==25 & true_down==0 ~ "One-sided (25/0)"
    )
  )
retain_tbl <- rho_sum %>%
  group_by(set_type, rho) %>%
  summarise(
    n_kept = n_distinct(pathway),
    .groups = "drop"
  ) %>%
  group_by(set_type) %>%
  mutate(
    n_total = length(gene_sets),
    retain_path = n_kept / length(gene_sets)
  ) %>%
  ungroup() %>%
  arrange(set_type, rho)
write.csv(retain_tbl, file.path(table_dir, "supplement_Table_S5_ablation_correlation_retains.csv"))
