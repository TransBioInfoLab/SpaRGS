# making DLPFC tables
library(dplyr)

# path
project_dir <- "."
results_dir <- file.path(project_dir, "results", "DLPFC")
table_dir <- file.path(project_dir, "results", "tables")

# generate top20 summary meta-analysis table
df <- read.csv(file.path(results_dir, "hex200_meta_stouffer_res.csv"))
df$X <- NULL
df_top20 <- df %>%
  arrange(desc(n_sig), desc(median_obs_moran), median_obs_geary) %>%
  slice_head(n = 20)
write.csv(df_top20, file.path(table_dir, "supplement_Table_S7_DLPFC_top20.csv"))
