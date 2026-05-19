# making Xenium figures
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(grid)
library(hexbin)

source("R/sgsva_score.R")
source("R/spatial_permutation.R")

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
data_dir <- file.path(project_dir, "data", "xenium")
results_dir <- file.path(project_dir, "results", "xenium")
fig_dir <- file.path(project_dir, "results", "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------
# load saved outputs
# ---------------------------
# sGSVA
sgsva_res <- readRDS(file.path(results_dir, "sgsva_score_res.rds"))
z_sgsva <- sgsva_res$z_score
balanced_sets <- sgsva_res$gene_sets

# GSVA
gsva_scores <- read.csv(file.path(results_dir, "gsva_score.csv"))
rownames(gsva_scores) <- gsva_scores$X
gsva_scores$X <- NULL
z_gsva <- z_normalization(as.matrix(gsva_scores))
colnames(z_gsva) <- gsub("\\.","-",colnames(z_gsva))

# GSVA-smoothed
gsva_smooth <- read.csv(file.path(results_dir, "gsva_smooth_score.csv"))
rownames(gsva_smooth) <- gsva_smooth$X
gsva_smooth$X <- NULL
z_gsva_smooth <- z_normalization(as.matrix(gsva_smooth))
colnames(z_gsva_smooth) <- gsub("\\.","-",colnames(z_gsva_smooth))

# Xenium meta data with cell type proportions
meta <- read.csv(file.path(results_dir, "xenium_cell_type_meta_data.csv"))
rownames(meta) <- meta$X
meta$X <- NULL
colnames(meta) <- gsub("\\.", " ", colnames(meta))
# Xenium coordinates
coords <- read.csv(file.path(results_dir, "xenium_coordinates.csv"))
rownames(coords) <- coords$X
coords$X <- NULL

# ---------------------------
# Built hex bins
# ---------------------------
coords_df <- data.frame(
  cell_id = rownames(coords),
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

# make sure order matches
common_cells <- intersect(colnames(z_sgsva), coords_df$cell_id)
coords_sub   <- coords_df[match(common_cells, coords_df$cell_id), ]
z_sgsva_sub  <- z_sgsva[, common_cells, drop = FALSE]
z_gsva_sub   <- z_gsva[,  common_cells, drop = FALSE]
z_smooth_sub <- z_gsva_smooth[, common_cells, drop = FALSE]


# ---------------------------
# helper functions
# ---------------------------
# aggregate to bin level
aggregate_bin_level <- function(pathway, ct_cols) {
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

# bin-level scatter plots
plot_bin_scatter <- function(bin_scores, ct_col, ct_label, pathway_name) {
  df <- bin_scores %>%
    transmute(
      ct = .data[[ct_col]],
      sGSVA       = sgsva_mean,
      GSVA        = gsva_mean,
      GSVA_smoothed = smooth_mean
    ) %>%
    tidyr::pivot_longer(
      cols      = c("sGSVA", "GSVA", "GSVA_smoothed"),
      names_to  = "method",
      values_to = "score"
    ) %>%
    filter(!is.na(ct), !is.na(score))
  
  ggplot(df, aes(x = ct, y = score)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    ggpubr::stat_cor(method = "spearman",
                     label.x.npc = "left",
                     label.y.npc = "top",
                     size = 4) +
    facet_wrap(~ method, ncol = 3) +
    labs(
      x = ct_label,
      y = paste0("Bin-level z-score (", pathway_name, ")"),
      title = paste0(pathway_name, " vs ", ct_label, " (bin level)")
    ) +
    theme_minimal()
}

plot_bin_maps <- function(
    bin_scores, ct_col, ct_title, bin_size = 1.8,
    pathway_quantile_clip = 0.99, ct_quantile_clip = 0.99
) {
  pathway_vals <- c(bin_scores$sgsva_mean,
                    bin_scores$gsva_mean,
                    bin_scores$smooth_mean)
  pathway_vals <- pathway_vals[is.finite(pathway_vals)]
  
  # symmetric clipping for pathway z-scores
  lim_path <- quantile(abs(pathway_vals), probs = pathway_quantile_clip, na.rm = TRUE)
  lim_path <- max(1, lim_path)  # avoid too-small limits
  
  # clip cell-type score for cleaner display
  ct_vals <- bin_scores[[ct_col]]
  ct_vals <- ct_vals[is.finite(ct_vals)]
  lim_ct <- quantile(ct_vals, probs = c(1 - ct_quantile_clip, ct_quantile_clip), na.rm = TRUE)
  
  map_theme <- theme_void(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.height = unit(0.55, "cm"),
      legend.key.width  = unit(0.35, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    )
  
  p_ct <- ggplot(bin_scores, aes(x = bx, y = by, color = .data[[ct_col]])) +
    geom_point(size = bin_size, shape = 16, stroke = 0) +
    coord_equal() +
    scale_color_gradient(
      low = "grey92",
      high = "#D7301F",
      limits = lim_ct,
      oob = squish,
      name = "Cell-type\nscore"
    ) +
    labs(title = ct_title) +
    map_theme
  
  pathway_scale <- scale_color_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#D7301F",
    midpoint = 0,
    limits = c(-lim_path, lim_path),
    oob = squish,
    name = "z-score"
  )
  
  p_sg <- ggplot(bin_scores, aes(x = bx, y = by, color = sgsva_mean)) +
    geom_point(size = bin_size, shape = 16, stroke = 0) +
    coord_equal() +
    pathway_scale +
    labs(title = "sGSVA") +
    map_theme
  p_gs <- ggplot(bin_scores, aes(x = bx, y = by, color = gsva_mean)) +
    geom_point(size = bin_size, shape = 16, stroke = 0) +
    coord_equal() +
    pathway_scale +
    labs(title = "GSVA") +
    map_theme
  p_sm <- ggplot(bin_scores, aes(x = bx, y = by, color = smooth_mean)) +
    geom_point(size = bin_size, shape = 16, stroke = 0) +
    coord_equal() +
    pathway_scale +
    labs(title = "GSVA-smoothed") +
    map_theme
  
  combined_plot <-
    (p_ct | p_sg | p_gs | p_sm) +
    plot_layout(guides = "collect") 
  return(combined_plot)
}


# -----------------------------
# Xenium main figures
# -----------------------------
# aggregate bin-level scores
bin_scores_chemokine <- aggregate_bin_level("KEGG_CHEMOKINE_SIGNALING_PATHWAY", ct_col = "GG_T cells")
bin_scores_cellcycle <- aggregate_bin_level("KEGG_CELL_CYCLE", ct_col = "GG_Breast cancer")

# sGSVA vs GSVA vs GSVA-smoothed three maps
map_chemokine_main <- plot_bin_maps(bin_scores_chemokine, ct_col = "GG_T cells_mean",
                                   ct_title  = "T-cell module (bin mean)")
map_cellcycle_main <- plot_bin_maps(bin_scores_cellcycle, ct_col = "GG_Breast cancer_mean",
                                   ct_title = "Breast cancer module (bin mean)")
row1_titled <- map_cellcycle_main + plot_annotation(
  title = "KEGG Cell Cycle",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
)
row2_titled <- map_chemokine_main + plot_annotation(
  title = "KEGG Chemokine Signaling Pathway",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
)
comb_map <- wrap_plots(
  wrap_elements(panel = row1_titled),
  wrap_elements(panel = row2_titled),
  ncol = 1,
  heights = c(1, 1)
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))
ggsave(filename=file.path(fig_dir, "Fig6_Xenium_main.pdf"),
       plot=comb_map, width = 12, height = 10)



# -----------------------------
# Xenium supplementary figures
# -----------------------------
# bin-level scatter plots
scatter_cellcycle <- plot_bin_scatter(bin_scores_cellcycle, ct_col = "GG_Breast cancer_mean",
                                      ct_label = "Breast cancer module",
                                      pathway_name = "KEGG cell cycle")
scatter_chemokine <- plot_bin_scatter(bin_scores_chemokine, ct_col = "GG_T cells_mean",
                                      ct_label = "T-cell module",
                                      pathway_name = "KEGG chemokine signaling pathway")
comb_scatter <- scatter_cellcycle/scatter_chemokine +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")
ggsave(filename=file.path(fig_dir, "supplement_Fig_S11_Xenium_scatter.pdf"),
       plot=comb_scatter, width = 12, height = 12)

