# making robustness figures
library(patchwork)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

source("R/spatial_permutation.R")
source("R/simulation.R")
source("R/simulation_scoring.R")
source("R/evaluation_metrics.R")

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
results_dir <- file.path(project_dir, "results", "robustness")
fig_dir <- file.path(project_dir, "results", "figures")

# ---------------------------
# load saved outputs
# ---------------------------
# combine bandwidth outputs from oracle and data-driven directionality
bandwidth_oracle <- read.csv(file.path(results_dir, "abl_bandwidths_aggregate_sum_oracle.csv"))
bandwidth_oracle$X <- NULL
bandwidth_data <- read.csv(file.path(results_dir, "abl_bandwidths_aggregate_sum_data_driven.csv"))
bandwidth_data$X <- NULL
bandwidth_comb <- bandwidth_oracle %>%
  full_join(bandwidth_data, by=c("pathway","ell","sigma_f","set_type","sigma_param","sigma_global")) %>%
  dplyr::rename(n_oracle=n.x, mean_auroc_oracle=mean_auroc.x,
                se_auroc_oracle=se_auroc.x, mean_spearman_oracle=mean_spearman.x,
                se_spearman_oracle=se_spearman.x, mean_times_oracle=mean_times.x,
                se_times_oracle=se_times.x, n_data=n.y,
                mean_auroc_data=mean_auroc.y, se_auroc_data=se_auroc.y,
                mean_spearman_data=mean_spearman.y, se_spearman_data=se_spearman.y,
                mean_times_data=mean_times.y, se_times_data=se_times.y)

# combine weight clipping outputs from oracle and data-driven directionality
clip_oracle <- read.csv(file.path(results_dir, "abl_weight_clipping_aggregate_sum_oracle.csv"))
clip_oracle$X <- NULL
clip_data <- read.csv(file.path(results_dir, "abl_weight_clipping_aggregate_sum_data_driven.csv"))
clip_data$X <- NULL
clip_comb <- clip_oracle %>%
  full_join(clip_data, by=c("pathway","ell","sigma_f","set_type","weight_clip")) %>%
  dplyr::rename(n_oracle=n.x, mean_auroc_oracle=mean_auroc.x,
                se_auroc_oracle=se_auroc.x, mean_times_oracle=mean_times.x,
                se_times_oracle=se_times.x, mean_moran_oracle=mean_moran.x,
                se_moran_oracle=se_moran.x, 
                n_data=n.y, mean_auroc_data=mean_auroc.y,
                se_auroc_data=se_auroc.y, mean_times_data=mean_times.y,
                se_times_data=se_times.y, mean_moran_data=mean_moran.y,
                se_moran_data=se_moran.y) 

# spatial graph output
res_graph <- read.csv(file.path(results_dir, "abl_spatial_graph_moran_oracle.csv"))
res_graph$X <- NULL

# block stability output
res_block <- read.csv(file.path(results_dir, "abl_block_stability_autocorrelation_oracle.csv"))
res_block$X <- NULL

# -------------------------------
# ablation main figures
# -------------------------------
# Shared theme for manuscript figures
theme_manuscript <- function(base_size = 15) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.9, "lines"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# bandwidth plot
bandwidth_comb <- bandwidth_oracle %>%
  filter(!(sigma_f == 0.5 & set_type == "One-sided (25/0)")) %>%
  full_join(
    bandwidth_data,
    by = c("pathway", "ell", "sigma_f", "set_type", "sigma_param", "sigma_global")
  ) %>%
  dplyr::rename(
    mean_auroc_oracle = mean_auroc.x,
    mean_auroc_data = mean_auroc.y
  ) %>%
  mutate(
    auc_diff = mean_auroc_oracle - mean_auroc_data,
    ell = factor(
      ell,
      levels = c(50, 200, 400),
      labels = c(expression(ell == 50), expression(ell == 200), expression(ell == 400))
    ),
    sigma_param = factor(
      sigma_param,
      levels = c("1/sqrt(2)*ell", "1*ell", "1.5*ell"),
      labels = c(expression(ell / sqrt(2)), expression(ell), expression(1.5 * ell))
    )
  ) %>%
  filter(!is.na(auc_diff))
summary(bandwidth_comb$auc_diff)

p_band <- ggplot(
  bandwidth_comb,
  aes(x = factor(sigma_f), y = auc_diff, fill = sigma_param)) +
  geom_boxplot(
    width = 0.65,
    alpha = 0.85,
    outlier.size = 0.7,
    linewidth = 0.35,
    position = position_dodge(width = 0.75)
  ) +
  geom_hline(yintercept = 0, linetype = 2, color = "firebrick", linewidth = 0.4) +
  facet_grid(set_type ~ ell, labeller = label_parsed) +
  scale_fill_discrete(
    name = expression(sigma~"relative to"~ell),
    labels = c(expression(ell / sqrt(2)), expression(1 %*% ell), expression(1.5 %*% ell))
  ) +
  labs(
    x = expression("Signal amplitude " * sigma[f]),
    y = expression(Delta * "AUROC: oracle - data-driven"),
    title = "Oracle-data-driven gap across bandwidths"
  ) +
  theme_manuscript(base_size = 15) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


# weight clipping plot
clip_comb <- clip_oracle %>%
  filter(!(sigma_f == 0.5 & set_type == "One-sided (25/0)")) %>%
  full_join(
    clip_data,
    by = c("pathway", "ell", "sigma_f", "set_type", "weight_clip")
  ) %>%
  dplyr::rename(
    mean_auroc_oracle = mean_auroc.x,
    mean_auroc_data = mean_auroc.y
  ) %>%
  mutate(
    auc_diff = mean_auroc_oracle - mean_auroc_data,
    ell = factor(
      ell,
      levels = c(50, 200, 400),
      labels = c(expression(ell == 50), expression(ell == 200), expression(ell == 400))
    )
  ) %>%
  filter(!is.na(auc_diff))
summary(clip_comb$auc_diff)

p_clip <- ggplot(
  clip_comb, aes(x = factor(sigma_f), y = auc_diff, fill = factor(weight_clip))) +
  geom_boxplot(
    width = 0.65,
    alpha = 0.85,
    outlier.size = 0.7,
    linewidth = 0.35,
    position = position_dodge(width = 0.75)
  ) +
  geom_hline(yintercept = 0, linetype = 2, color = "firebrick", linewidth = 0.4) +
  facet_grid(set_type ~ ell, labeller = label_parsed) +
  labs(
    x = expression("Signal amplitude " * sigma[f]),
    y = expression(Delta * "AUROC: oracle - data-driven"),
    fill = "CV weight\nclipping",
    title = "Oracle-data-driven gap across weight clippings"
  ) +
  theme_manuscript(base_size = 15) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


# spatial graph plot
# helper function
get_graph_ablation_df <- function(dat) {
  ell_vals <- sort(unique(dat$ell))
  map_dfr(ell_vals, function(l) {
    sub <- dat %>% filter(ell == l)
    
    wide <- sub %>%
      dplyr::select(graph, pathway, z_moran_hex) %>%
      tidyr::pivot_wider(names_from = graph, values_from = z_moran_hex)
    rank_mat <- apply(wide[, -1], 2, rank, ties.method = "average", na.last = "keep")
    cor_mat <- cor(rank_mat, method = "spearman", use = "pairwise.complete.obs")
    
    as.data.frame(cor_mat) %>%
      tibble::rownames_to_column("row_graph") %>%
      tidyr::pivot_longer(
        -row_graph,
        names_to = "col_graph",
        values_to = "value"
      ) %>%
      mutate(metric = "Rank correlation", ell = l)
  })
}
# rank correlation heatmap across ells
dat <- res_graph %>%
  mutate(
    graph = case_when(
      mode == "delaunay" ~ "Delaunay",
      mode == "knn" & k == 6 ~ "KNN-6",
      mode == "knn" & k == 8 ~ "KNN-8",
      mode == "knn" & k == 12 ~ "KNN-12"
    ),
    graph = as.factor(graph)
  )
ablation_df <- get_graph_ablation_df(dat)
p_graph <- ablation_df %>%
  mutate(
    row_graph = factor(row_graph, levels = c("Delaunay", "KNN-6", "KNN-8", "KNN-12")),
    col_graph = factor(col_graph, levels = c("Delaunay", "KNN-6", "KNN-8", "KNN-12")),
    ell = factor(
      ell,
      levels = c(50, 200, 400),
      labels = c(expression(ell == 50), expression(ell == 200), expression(ell == 400))
    )
  ) %>%
  ggplot(aes(x = col_graph, y = row_graph, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 4.8, fontface = "bold") +
  scale_fill_viridis_c(
    limits = c(0, 1),
    name = expression("Spearman " * rho)
  ) +
  coord_equal() +
  facet_grid(. ~ ell, labeller = label_parsed) +
  labs(
    x = NULL,
    y = NULL,
    title = "Graph sensitivity"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    panel.spacing.x = unit(1.1, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

# block stability plot
# Rank correlation heatmap
wide <- res_block %>%
  dplyr::select(pathway, block_size, z_moran_hex) %>%
  tidyr::pivot_wider(names_from = block_size, values_from = z_moran_hex, names_glue = "Block {block_size}")
# Convert to ranks
rank_mat <- apply(wide[,-1], 2, rank, ties.method = "average", na.last = "keep")
cor_mat <- cor(rank_mat, method = "spearman", use = "pairwise.complete.obs")
p_R <- as.data.frame(cor_mat) %>%
  tibble::rownames_to_column("row_block") %>%
  tidyr::pivot_longer(
    -row_block,
    names_to = "col_block",
    values_to = "value"
  ) %>%
  mutate(
    row_block = factor(row_block, levels = c("Block 150", "Block 200", "Block 300")),
    col_block = factor(col_block, levels = c("Block 150", "Block 200", "Block 300"))
  ) %>%
  ggplot(aes(x = col_block, y = row_block, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 5.2, fontface = "bold") +
  scale_fill_viridis_c(
    limits = c(0, 1),
    name = expression("Spearman " * rho)
  ) +
  coord_equal() +
  labs(
    x = NULL,
    y = NULL,
    title = "Hex-block sensitivity"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(5, 5, 5, 5)
  )

# combine into main figure
combined <- (p_band | p_clip) / (p_graph | p_R) +
  plot_annotation(tag_levels = "A") +
  plot_layout(
    heights = c(1.05, 0.9),
    guides = "collect"
  ) &
  theme(
    legend.position = "right",
    plot.tag = element_text(size = 20, face = "bold")
  )
ggsave(
  filename = file.path(fig_dir, "Fig3_ablation_robustness.pdf"),
  plot = combined,
  width = 15,
  height = 12
)


# -----------------------------------
# ablation supplementary figures
# -----------------------------------
# ============================================================
# Supplementary Figure S1
# Direct sensitivity of data-driven sGSVA to bandwidth and
# CV weight-clipping choices
#
# Delta AUROC = alternative setting - default setting
# Default bandwidth: ell / sqrt(2)
# Default weight clipping: 0.99
# ============================================================

# ------------------------------------------------------------
# A. Bandwidth sensitivity relative to default
# ------------------------------------------------------------

bandwidth_default <- "1/sqrt(2)*ell"
bandwidth_default_df <- bandwidth_data %>%
  filter(sigma_param == bandwidth_default) %>%
  dplyr::select(pathway, ell, sigma_f, set_type,
    default_auroc = mean_auroc
  ) %>%
  distinct()

bandwidth_sens <- bandwidth_data %>%
  dplyr::select(pathway, ell, sigma_f, set_type, sigma_param, mean_auroc) %>%
  left_join(
    bandwidth_default_df,
    by = c("pathway", "ell", "sigma_f", "set_type")
  ) %>%
  mutate(
    delta_from_default = mean_auroc - default_auroc,
    set_type = factor(
      set_type,
      levels = c("Balanced (25/25)", "One-sided (25/0)")
    ),
    ell_label = factor(
      ell,
      levels = c(50, 200, 400),
      labels = c(expression(ell == 50), expression(ell == 200), expression(ell == 400))
    ),
    sigma_param_plot = factor(
      sigma_param,
      levels = c("1*ell", "1.5*ell"),
      labels = c(expression(ell), expression(1.5 * ell)
      )
    )
  ) %>%
  # Do not plot the default setting because its delta is always zero
  filter(
    sigma_param != bandwidth_default,
    !is.na(delta_from_default),
    !is.na(sigma_param_plot)
  )

p_band_sens <- ggplot(
  bandwidth_sens,
  aes(x = factor(sigma_f), y = delta_from_default, fill = sigma_param_plot)) +
  geom_hline(yintercept = 0, linetype = 2, color = "firebrick", linewidth = 0.4) +
  geom_boxplot(
    width = 0.65, alpha = 0.85, outlier.size = 0.7, linewidth = 0.35,
    position = position_dodge(width = 0.75),
    na.rm = TRUE
  ) +
  facet_grid(
    set_type ~ ell_label,
    labeller = label_parsed
  ) +
  labs(
    x = expression("Signal amplitude " * sigma[f]),
    y = expression(
      Delta * "AUROC (alternative " - " default)"
    ),
    fill = expression(sigma ~ "relative to" ~ ell),
    title = "Sensitivity to smoothing bandwidth"
  ) +
  theme_manuscript(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )


# ------------------------------------------------------------
# B. CV weight-clipping sensitivity relative to default
# ------------------------------------------------------------

weight_clip_default <- 0.99
clip_default_df <- clip_data %>%
  filter(dplyr::near(weight_clip, weight_clip_default)) %>%
  dplyr::select(pathway, ell, sigma_f, set_type,
    default_auroc = mean_auroc
  ) %>%
  distinct()

clip_sens <- clip_data %>%
  dplyr::select(pathway, ell, sigma_f, set_type, weight_clip, mean_auroc) %>%
  left_join(
    clip_default_df,
    by = c("pathway", "ell", "sigma_f", "set_type")
  ) %>%
  mutate(
    delta_from_default = mean_auroc - default_auroc,
    set_type = factor(
      set_type,
      levels = c("Balanced (25/25)", "One-sided (25/0)")
    ),
    ell_label = factor(
      ell,
      levels = c(50, 200, 400),
      labels = c(expression(ell == 50), expression(ell == 200), expression(ell == 400))
    ),
    weight_clip_plot = factor(
      weight_clip,
      levels = c(0.95, 0.995),
      labels = c("0.95", "0.995")
    )
  ) %>%
  # Do not plot the default setting because its delta is always zero
  filter(
    !dplyr::near(weight_clip, weight_clip_default),
    !is.na(delta_from_default),
    !is.na(weight_clip_plot)
  )

p_clip_sens <- ggplot(clip_sens,
  aes(x = factor(sigma_f), y = delta_from_default, fill = weight_clip_plot)) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    color = "firebrick",
    linewidth = 0.4
  ) +
  geom_boxplot(
    width = 0.65,
    alpha = 0.85,
    outlier.size = 0.7,
    linewidth = 0.35,
    position = position_dodge(width = 0.75),
    na.rm = TRUE
  ) +
  facet_grid(
    set_type ~ ell_label,
    labeller = label_parsed
  ) +
  labs(
    x = expression("Signal amplitude " * sigma[f]),
    y = expression(
      Delta * "AUROC (alternative " - " default)"
    ),
    fill = "CV weight-clipping\nquantile",
    title = "Sensitivity to CV weight clipping"
  ) +
  theme_manuscript(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )


# ------------------------------------------------------------
# Combine bandwidth and clipping panels
# ------------------------------------------------------------
comb_sens <- p_band_sens / p_clip_sens +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1), guides = "keep") &
  theme(
    plot.tag = element_text(size = 18, face = "bold")
  )

ggsave(
  filename = file.path(
    fig_dir,
    "supplement_Fig_S1_robustness_bandwidth_clipping_sensitivity.pdf"
  ),
  plot = comb_sens,
  width = 12,
  height = 12
)


# ============================================================
# supplementary figure S2 - hex-block stability (Jaccard)
# ============================================================

# helper function
get_top_k_jaccard <- function(df, rank_column = z_moran_hex, k_values = c(10, 20, 30)) {
  block_levels <- sort(unique(df$block_size))
  out_list <- list()
  
  for (k in k_values) {
    top_k_sets <- lapply(block_levels, function(b) {
      df %>% filter(block_size == b) %>%
        slice_max({{rank_column}}, n = k, with_ties = FALSE) %>%
        pull(pathway)
    })
    names(top_k_sets) <- as.character(block_levels)
    
    jaccard_mat <- outer(block_levels, block_levels, Vectorize(function(i, j) {
      A <- top_k_sets[[as.character(i)]]
      B <- top_k_sets[[as.character(j)]]
      length(intersect(A, B)) / length(union(A, B))
    }))
    colnames(jaccard_mat) <- rownames(jaccard_mat) <- paste0("block", block_levels)
    out_list[[paste0("Top", k)]] <- jaccard_mat
  }
  return(out_list)
}

# plot Jaccard heatmaps for Moran's I
jaccard_list_moran <- get_top_k_jaccard(res_block, rank_column = z_moran_hex)
jaccard_map <- list()
for (name in names(jaccard_list_moran)) {
  jacc_mat <- jaccard_list_moran[[name]]
  jac_df <- as.data.frame(jacc_mat) %>%
    tibble::rownames_to_column("row_graph") %>%
    tidyr::pivot_longer(-row_graph, names_to = "col_graph",
                 values_to = "value")
  p_J <- jac_df %>%
    ggplot(aes(x = col_graph, y = row_graph, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value), color = value > 0.5), size = 3) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none")+
    scale_fill_viridis_c(limits = c(0, 1)) +
    coord_equal() +
    labs(x = NULL, y = NULL, fill = "Value", title = paste0("Jaccard ", name, " overlap across hex-block sizes")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.spacing.x = unit(1, "lines"),
          panel.spacing.y = unit(0.5, "lines"))
  jaccard_map[[name]] <- p_J
}

comb_jac <- (
  jaccard_map$Top10 + 
    jaccard_map$Top20 +
    jaccard_map$Top30
) +
  plot_layout(
    ncol = 2, widths = c(1, 1),
    heights = c(1, 1), guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    plot.tag = element_text(size = 16, face = "bold")
  )
ggsave(filename=file.path(fig_dir, "supplement_Fig_S2_robustness_hexblock_jaccard.pdf"),
       plot=comb_jac, width = 12, height = 10)






