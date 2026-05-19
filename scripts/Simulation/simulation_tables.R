# making simulation tables
library(dplyr)

# paths
project_dir <- "."
results_dir <- file.path(project_dir, "results", "simulation")
table_dir <- file.path(project_dir, "results", "tables")

# load save outputs
region_raw <- read.csv(file.path(results_dir, "AUROC_spearman_runtime_all_replicates_oracle.csv"))
region_raw$X <- NULL
region_sum <- read.csv(file.path(results_dir, "auroc_aggregate_sum_oracle.csv"))
region_sum$X <- NULL
cor_sum <- read.csv(file.path(results_dir, "spearman_aggregate_sum_oracle.csv"))
cor_sum$X <- NULL

# generate summary tables
# Table 1
scen_summary <- region_raw %>%
  mutate(set_type = if_else(down_size == 0, "One-sided (25/0)", "Balanced (25/25)")) %>%
  group_by(methods, set_type, ell, sigma_f) %>%
  summarise(
    auroc = mean(AUROC, na.rm = TRUE),
    spear = mean(spearman, na.rm = TRUE),
    runtime = mean(times, na.rm = TRUE),
    .groups = "drop"
  )
sum_table <- scen_summary %>%
  group_by(methods) %>%
  summarise(
    n_scen = n(),
    mean_auroc = mean(auroc, na.rm = TRUE),
    sd_auroc   = sd(auroc, na.rm = TRUE),
    mean_spearman = mean(spear, na.rm = TRUE),
    sd_spearman   = sd(spear, na.rm = TRUE),
    mean_runtime  = mean(runtime, na.rm = TRUE),
    .groups = "drop"
  )
writexl::write_xlsx(sum_table, file.path(table_dir, "table1_simulation_summary.xlsx")) 

# supplementary Table S1-4
# AUROC tables
auroc_sum <- region_sum %>%
  mutate(
    set_type = case_when(
      up_size==25 & down_size==25 ~ "Balanced (25/25)",
      up_size==25 & down_size==0 ~ "One-sided (25/0)"
    )
  )
auroc_balanced <- auroc_sum %>%
  filter(set_type=="Balanced (25/25)") %>%
  group_by(methods, ell) %>%
  summarise(mean_auroc = mean(mean_auroc, na.rm=TRUE))
write.csv(auroc_balanced, file.path(table_dir, "supplement_Table_S1_simulation_auroc_balanced.csv"))
auroc_onesided <- auroc_sum %>%
  filter(set_type=="One-sided (25/0)") %>%
  group_by(methods, ell) %>%
  summarise(mean_auroc = mean(mean_auroc, na.rm=TRUE))
write.csv(auroc_onesided, file.path(table_dir, "supplement_Table_S2_simulation_auroc_onesided.csv"))

# spearman tables
cor_sum <- cor_sum %>%
  mutate(
    set_type = case_when(
      up_size==25 & down_size==25 ~ "Balanced (25/25)",
      up_size==25 & down_size==0 ~ "One-sided (25/0)"
    )
  )
cor_balanced <- cor_sum %>%
  filter(set_type=="Balanced (25/25)") %>%
  group_by(methods, ell) %>%
  summarise(mean_cor = mean(mean_cor, na.rm=TRUE))
write.csv(cor_balanced, file.path(table_dir, "supplement_Table_S3_simulation_spearman_balanced.csv"))
cor_onesided <- cor_sum %>%
  filter(set_type=="One-sided (25/0)") %>%
  group_by(methods, ell) %>%
  summarise(mean_cor = mean(mean_cor, na.rm=TRUE))
write.csv(cor_onesided, file.path(table_dir, "supplement_Table_S4_simulation_spearman_onesided.csv"))



