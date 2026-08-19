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
xenium_data <- readRDS(file.path(results_dir, "xenium_process_data.rds"))
coords <- xenium_data$coords

ct_cols <- grep("^GG_", colnames(meta), value=TRUE)
resolve_col <- function(pattern, cols) {
  hit <- grep(pattern, cols, value = TRUE)
  if (length(hit) != 1) {
    stop("Could not resolve column for pattern: ", pattern,
         "\nHits: ", paste(hit, collapse = ", "))
  }
  hit
}


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
      name = "Module score"
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


# pathway-cell-type module concordance heatmap
prespecified_xenium_pairs <- tibble::tribble(
  ~pathway,                              ~anchor_pattern,       ~anchor_label,
  "KEGG_CELL_CYCLE",                     "Breast cancer",       "Cell cycle - breast cancer",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",    "T cells",             "Chemokine - T cells",
  "KEGG_VEGF_SIGNALING_PATHWAY",         "Endothelial",         "VEGF - endothelial",
  "KEGG_FOCAL_ADHESION",                 "Fibroblasts",         "Focal adhesion - fibroblasts"
)

compute_pair_concordance <- function(pathway, anchor_pattern) {
  ct_col <- resolve_col(anchor_pattern, ct_cols)
  ct_mean <- paste0(ct_col, "_mean")
  
  bin_scores <- aggregate_bin_level(pathway, ct_col) %>%
    mutate(sGSVA = sgsva_mean,
           GSVA = gsva_mean,
           `GSVA-smoothed` = smooth_mean)
  
  method_cols <- c("sGSVA", "GSVA", "GSVA-smoothed")
  
  bind_rows(lapply(method_cols, function(m) {
    rho <- suppressWarnings(
      cor(
        bin_scores[[m]],
        bin_scores[[ct_mean]],
        method = "spearman",
        use = "complete.obs"
      )
    )
    
    tibble::tibble(
      pathway = pathway,
      anchor = anchor_pattern,
      method = m,
      signed_rho = rho,
      abs_rho = abs(rho),
      n_bins = sum(is.finite(bin_scores[[m]]) & is.finite(bin_scores[[ct_mean]]))
    )
  }))
}

xenium_concordance <- bind_rows(
  lapply(seq_len(nrow(prespecified_xenium_pairs)), function(i) {
    out <- compute_pair_concordance(
      prespecified_xenium_pairs$pathway[i],
      prespecified_xenium_pairs$anchor_pattern[i]
    )
    out$anchor_label <- prespecified_xenium_pairs$anchor_label[i]
    out
  })
)

p_xenium_concordance <- xenium_concordance %>%
  mutate(
    method = factor(method, levels = c("sGSVA", "GSVA", "GSVA-smoothed")),
    anchor_label = factor(anchor_label, levels = rev(unique(anchor_label))),
    signed_label = sprintf("%+.2f", signed_rho)
  ) %>%
  ggplot(aes(x = method, y = anchor_label, fill = abs_rho)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = signed_label), size = 3.8) +
  scale_fill_gradient(
    limits = c(0, 1),
    low = "white",
    high = "steelblue",
    oob = scales::squish,
    name = expression("|Spearman " * rho * "|")
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Magnitude of pathway-module spatial concordance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 4)
    ),
    plot.margin = margin(4, 4, 4, 4)
  )


# Panel title themes
row_title_theme <- theme(
  plot.title = element_text(
    size = 16,
    face = "bold",
    hjust = 0.5,
    margin = margin(b = 6)
  )
)

# Row titles
row_cellcycle_clean <- map_cellcycle_main +
  patchwork::plot_annotation(
    title = "KEGG Cell Cycle",
    theme = row_title_theme
  )

row_chemokine_clean <- map_chemokine_main +
  patchwork::plot_annotation(
    title = "KEGG Chemokine Signaling Pathway",
    theme = row_title_theme
  )
heatmap_row <- (
  patchwork::plot_spacer() |
    p_xenium_concordance |
    patchwork::plot_spacer()
) +
  patchwork::plot_layout(
    widths = c(0.60, 2.80, 0.60)
  )

fig_xenium_main <- (
  patchwork::wrap_elements(full = heatmap_row) /
    patchwork::wrap_elements(full = row_cellcycle_clean) /
    patchwork::wrap_elements(full = row_chemokine_clean)
) +
  patchwork::plot_layout(
    heights = c(0.90, 1.15, 1.15)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.tag = element_text(
      size = 18,
      face = "bold"
    ),
    plot.tag.position = c(0.01, 0.98),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(filename=file.path(fig_dir, "Fig6_Xenium_main.pdf"),
       plot=fig_xenium_main, width = 12, height = 13)

