# making Xenium tables
library(dplyr)
library(hexbin)

source("R/sgsva_score.R")

# Required source files for Xenium table generation
required_files <- c(
  "results/xenium/sgsva_score_res.rds",
  "results/xenium/gsva_score.csv",
  "results/xenium/gsva_smooth_score.csv",
  "results/xenium/xenium_cell_type_meta_data.csv",
  "results/xenium/xenium_coordinates.csv",
  "results/xenium/moran_geary_three_methods.csv"
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "\nMissing Xenium source files required for table generation:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    "\n\nThese large intermediate files are not included in the lightweight GitHub repository.\n",
    "To regenerate them, run:\n\n",
    "  Rscript scripts/Xenium/Xenium_analysis.R\n",
    "After these files are generated, rerun:\n\n",
    "  Rscript scripts/Xenium/Xenium_tables.R\n\n",
    "Alternatively, use the archived full analysis outputs if provided separately.\n",
    call. = FALSE
  )
}

# paths
project_dir <- "."
results_dir <- file.path(project_dir, "results", "xenium")
table_dir <- file.path(project_dir, "results", "tables")

## generate Table 2 for spatial autocorrelation summary
moran_geary_all <- read.csv(file.path(results_dir, "moran_geary_three_methods.csv"))
sum_moran_geary <- moran_geary_all %>%
  group_by(method) %>%
  summarize(
    median_moran = median(moran_I),
    min_moran = min(moran_I),
    max_moran = max(moran_I),
    median_geary = median(geary_C),
    min_geary = min(geary_C),
    max_geary = max(geary_C)
  )
write.csv(sum_moran_geary, file.path(table_dir, "table2_xenium_median_autocorrelation_summary.csv"))


# Built hex bins
coords <- read.csv(file.path(results_dir, "xenium_coordinates.csv"))
coords_df <- data.frame(
  cell_id = coords$X,
  x = coords[, "x"],
  y = coords[, "y"]
)

# choose hex size ≈ 25–60 µm
bin_size <- 40
# hexbin wants number of bins across x; approximate from range / size
xrange <- range(coords_df$x)
yrange <- range(coords_df$y)
n_x <- ceiling((xrange[2] - xrange[1]) / bin_size)

hb <- hexbin(x = coords_df$x, y = coords_df$y, xbins = n_x, IDs = TRUE)
coords_df$bin_id <- hb@cID        # each cell gets a bin index
bin_centers <- hcell2xy(hb)       # bin centers, rows correspond to cID values
cell_ids_hex <- sort(unique(hb@cell))
bin_df <- data.frame(
  bin_id = cell_ids_hex,
  bx     = bin_centers$x,
  by     = bin_centers$y
)

# read meta data
meta <- read.csv(file.path(results_dir, "xenium_cell_type_meta_data.csv"))
rownames(meta) <- meta$X
meta$X <- NULL
ct_cols <- grep("^GG_", colnames(meta), value = TRUE)

# sgsva
sgsva_res <- readRDS(file.path(results_dir, "sgsva_score_res.rds"))
z_sgsva <- sgsva_res$z_score
# gsva
gsva_scores <- read.csv(file.path(results_dir, "gsva_score.csv"))
rownames(gsva_scores) <- gsva_scores$X
gsva_scores$X <- NULL
z_gsva <- z_normalization(as.matrix(gsva_scores))
colnames(z_gsva) <- gsub("\\.","-",colnames(z_gsva))
# gsva-smoothed
gsva_smooth <- read.csv(file.path(results_dir, "gsva_smooth_score.csv"))
rownames(gsva_smooth) <- gsva_smooth$X
gsva_smooth$X <- NULL
z_gsva_smooth <- z_normalization(as.matrix(gsva_smooth))
colnames(z_gsva_smooth) <- gsub("\\.","-",colnames(z_gsva_smooth))

# make sure order matches
common_cells <- intersect(colnames(z_sgsva), rownames(meta))
z_sgsva_sub  <- z_sgsva[, common_cells, drop = FALSE]
z_gsva_sub   <- z_gsva[,  common_cells, drop = FALSE]
z_smooth_sub <- z_gsva_smooth[, common_cells, drop = FALSE]

# helper function
# aggregate to bin level
aggregate_bin_level <- function(pathway) {
  df <- data.frame(
    cell_id = common_cells,
    bin_id  = meta$bin_id,
    sgsva   = as.numeric(z_sgsva_sub[pathway, ]),
    gsva    = as.numeric(z_gsva_sub[pathway, ]),
    smooth  = as.numeric(z_smooth_sub[pathway, ])
  )
  
  # add GG_* columns
  df <- cbind(df, meta[common_cells, ct_cols, drop = FALSE])
  bin_scores <- df %>%
    group_by(bin_id) %>%
    summarise(
      sgsva_mean  = mean(sgsva,  na.rm = TRUE),
      gsva_mean   = mean(gsva,   na.rm = TRUE),
      smooth_mean = mean(smooth, na.rm = TRUE),
      across(all_of(ct_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
      n_cells     = n(),
      .groups     = "drop"
    ) %>%
    left_join(bin_df, by = "bin_id")
  
  bin_scores
}

bin_anchoring_stats <- function(bin_scores, ct_col) {
  x  <- bin_scores[[ct_col]]
  ys <- list(
    sGSVA        = bin_scores$sgsva_mean,
    GSVA         = bin_scores$gsva_mean,
    GSVA_smoothed  = bin_scores$smooth_mean
  )
  
  out <- lapply(names(ys), function(m) {
    y   <- ys[[m]]
    rho <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
    data.frame(
      method = m,
      rho    = rho,
      abs_rho = abs(rho)
    )
  })
  do.call(rbind, out)
}

## generate supplement Table S8
bin_scores_chemokine <- aggregate_bin_level("KEGG_CHEMOKINE_SIGNALING_PATHWAY")
# check chemokine pathway with multiple local fractions
anchors <- c("GG_T.cells_mean",
             "GG_Myeloid_mean",
             "GG_Fibroblasts_mean",
             "GG_Endothelial.cells_mean",
             "GG_total_immune_mean")

chem_anchoring <- do.call(rbind, lapply(anchors, function(ct_col) {
  stats <- bin_anchoring_stats(bin_scores_chemokine, ct_col)
  stats$anchor <- ct_col
  stats
}))
write.csv(chem_anchoring, file.path(table_dir, "supplement_Table_S8_xenium_chemokine_anchoring_summary.csv"))







