# making breast Visium table
library(dplyr)

# load results
project_dir <- "."
results_dir <- file.path(project_dir, "results", "breast_visium")
table_dir <- file.path(project_dir, "results", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(file.path(results_dir, "moran_geary_pathway_summary.csv"))
df$X <- NULL
anchor_dat <- read.csv(file.path(results_dir, "correlation_score_vs_celltype.csv"))
anchor_dat$X <- NULL

# generate top10 summary table
top10_df <- df %>%
  slice_max(z_moran, n=10) %>%
  mutate(p_moran="< 0.001",
         p_geary="< 0.001",
         B = 1000) %>%
  dplyr::select(pathway, obs_moran, z_moran, p_moran, p_adj_moran, B) %>%
  left_join(anchor_dat, by="pathway")
write.csv(top10_df, file.path(table_dir, "supplement_Table_S6_breast_top10.csv"))
