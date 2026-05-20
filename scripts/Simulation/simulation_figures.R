# making simulation figures
library(patchwork)
library(dplyr)
library(ggplot2)

source("R/spatial_permutation.R")
source("R/simulation.R")
source("R/simulation_scoring.R")
source("R/evaluation_metrics.R")


# Required source files for Simulation figure generation
required_files <- c(
  "results/simulation/AUROC_spearman_runtime_all_replicates_oracle.csv",
  "results/simulation/auroc_aggregate_sum_oracle.csv",
  "results/simulation/spearman_aggregate_sum_oracle.csv",
  "results/simulation/delongtest_aggregate_sum_oracle.csv",
  "results/simulation/moran_geary_aggregate_sum_oracle.csv",
  "results/simulation/simulation_replicate1.rds"
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "\nMissing Simulation source files required for figure generation:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    "\n\nThese large intermediate files are not included in the lightweight GitHub repository.\n",
    "To regenerate them, run:\n\n",
    "  Rscript scripts/Simulation/simulation_analysis.R\n",
    "After these files are generated, rerun:\n\n",
    "  Rscript scripts/Simulation/simulation_figures.R\n\n",
    "Alternatively, use the archived full analysis outputs if provided separately.\n",
    call. = FALSE
  )
}

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
results_dir <- file.path(project_dir, "results", "simulation")
fig_dir <- file.path(project_dir, "results", "figures")

# ---------------------------
# load saved outputs
# ---------------------------
region_raw <- read.csv(file.path(results_dir, "AUROC_spearman_runtime_all_replicates_oracle.csv"))
region_raw$X <- NULL
region_sum <- read.csv(file.path(results_dir, "auroc_aggregate_sum_oracle.csv"))
region_sum$X <- NULL
cor_sum <- read.csv(file.path(results_dir, "spearman_aggregate_sum_oracle.csv"))
cor_sum$X <- NULL
delong_sum <- read.csv(file.path(results_dir, "delongtest_aggregate_sum_oracle.csv"))
delong_sum$X <- NULL
moran_sum <- read.csv(file.path(results_dir, "moran_geary_aggregate_sum_oracle.csv"))
moran_sum$X <- NULL

sim_rep <- readRDS(file.path(results_dir, "simulation_replicate1.rds"))

# ---- gene sets ( MSigDB KEGG legacy ) ----
msig <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
gene_sets <- split(msig$gene_symbol, msig$gs_name)

sel_paths <- c("KEGG_CHEMOKINE_SIGNALING_PATHWAY",
               "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
               "KEGG_FATTY_ACID_METABOLISM")
# build spatial weights ( Delaunay triangulation )
coords <- make_hex_grid(nx=70, ny=70, pitch=50)
nb <- tri2nb(coords)
listw <- nb2listw(nb, style = "W")
islands <- which(card(nb)==0)
length(islands)

# -------------------------------
# simulation main figures
# -------------------------------
method_levels <- c(
  "GSVA",
  "GSVA_smoothed",
  "sGSVA",
  "ssGSEA",
  "AUCell",
  "UCell",
  "AddModuleScore"
)

method_labels <- c(
  "GSVA" = "GSVA",
  "GSVA_smoothed" = "GSVA-smoothed",
  "sGSVA" = "sGSVA",
  "ssGSEA" = "ssGSEA",
  "AUCell" = "AUCell",
  "UCell" = "UCell",
  "AddModuleScore" = "AddModuleScore"
)

## Optional: make sigma labels cleaner
sigma_labeller <- as_labeller(function(x) paste0(expression(sigma[f]), " = ", x))

## Shared theme for readability at final print size
theme_fig2 <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = "grey50"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.45, "cm")
  )

# AUROC plot
auroc_plot <- region_sum %>%
  mutate(
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0  ~ "One-sided (25/0)",
      TRUE ~ NA_character_
    ),
    set_type = factor(set_type, levels = c("Balanced (25/25)", "One-sided (25/0)")),
    methods = factor(methods, levels = method_levels)
  )
p1 <- ggplot(auroc_plot, aes(x = factor(ell), y = mean_auroc, fill = methods)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.35) +
  facet_grid(
    set_type ~ sigma_f,
    labeller = labeller(sigma_f = label_both)
  ) +
  scale_fill_discrete(labels = method_labels, drop = FALSE) +
  labs(
    x = expression(ell ~ "(" * mu * "m)"),
    y = "AUROC",
    fill = "Methods",
    title = "Simulation AUROC performance"
  ) +
  theme_fig2 +
  theme(legend.position = "right")

# spearman plot
cor_plot <- cor_sum %>%
  mutate(
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0  ~ "One-sided (25/0)",
      TRUE ~ NA_character_
    ),
    set_type = factor(set_type, levels = c("Balanced (25/25)", "One-sided (25/0)")),
    methods = factor(methods, levels = method_levels)
  )
p2 <- ggplot(cor_plot, aes(x = factor(ell), y = mean_cor, fill = methods)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.35) +
  facet_grid(
    set_type ~ sigma_f,
    labeller = labeller(sigma_f = label_both)
  ) +
  scale_fill_discrete(labels = method_labels, drop = FALSE) +
  coord_cartesian(ylim = c(-0.2, 1.0)) +
  labs(
    x = expression(ell ~ "(" * mu * "m)"),
    y = "Map-to-truth Spearman correlation",
    fill = "Methods",
    title = "Simulation map-to-truth Spearman"
  ) +
  theme_fig2 +
  theme(legend.position = "right")

# running time plot
time_plot <- region_raw %>%
  mutate(
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0  ~ "One-sided (25/0)",
      TRUE ~ NA_character_
    ),
    methods = factor(methods, levels = method_levels)
  )
p3 <- ggplot(time_plot, aes(x = methods, y = times, fill = methods)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.35) +
  scale_y_log10() +
  scale_x_discrete(labels = method_labels, drop = FALSE) +
  scale_fill_discrete(labels = method_labels, drop = FALSE) +
  labs(
    x = "Methods",
    y = expression("Runtime (log"[10] * " scale)"),
    title = "Runtime by method"
  ) +
  theme_fig2 +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 11),
    legend.position = "none"
  )

# DeLong test plot
delong_plot_df <- delong_sum %>%
  filter(method1 == "sGSVA") %>%
  mutate(
    set_type = factor(set_type, levels = c("Balanced (25/25)", "One-sided (25/0)")),
    method2 = factor(method2, levels = method_levels),
    comparison = paste0("sGSVA vs ", method_labels[as.character(method2)]),
    comparison = factor(
      comparison,
      levels = paste0("sGSVA vs ", method_labels[method_levels[method_levels != "sGSVA"]])
    ),
    smoothing_baseline = ifelse(method2 == "GSVA_smoothed", "GSVA-smoothed", "Other methods")
  ) %>%
  group_by(ell, sigma_f, set_type, method1, method2, comparison, smoothing_baseline) %>%
  summarise(
    n_path = dplyr::n(),
    mu = mean(auc_diff_mean, na.rm = TRUE),
    se = sd(auc_diff_mean, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
delong_plot <- ggplot(
  delong_plot_df,
  aes(
    x = factor(ell),
    y = mu,
    group = comparison,
    color = comparison,
    linewidth = smoothing_baseline
  )
) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
  geom_line() +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = mu - 1.96 * se, ymax = mu + 1.96 * se),
    width = 0.15,
    linewidth = 0.45
  ) +
  facet_grid(
    set_type ~ sigma_f,
    labeller = labeller(sigma_f = label_both)
  ) +
  scale_color_discrete(name = "Method pairs", drop = FALSE) +
  scale_linewidth_manual(
    values = c("GSVA-smoothed" = 1.15, "Other methods" = 0.55),
    guide = "none"
  ) +
  labs(
    x = expression(ell ~ "(" * mu * "m)"),
    y = expression(Delta * "AUROC relative to sGSVA"),
    title = expression(Delta * "AUROC relative to sGSVA")
  ) +
  theme_fig2 +
  theme(legend.position = "right")
  
combined <- (p1 | p2) / (delong_plot | p3) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(1, 0.8), guides = "collect")
combined & theme(legend.position = "bottom")
ggsave(
  filename = file.path(fig_dir, "Fig2_simulation_performance.pdf"),
  plot = combined,
  width = 20,
  height = 15,
  device = cairo_pdf
)


# -----------------------------------
# simulation supplementary figures
# -----------------------------------
# supplementary figure S1 for Spatial Autocorrelation Calibration
pathway_labels <- c(
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY" = "B cell receptor",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY" = "Chemokine signaling",
  "KEGG_FATTY_ACID_METABOLISM" = "Fatty acid metabolism"
)
moran_plot_df <- moran_sum %>%
  mutate(
    pathway_short = recode(pathway, !!!pathway_labels),
    pathway_short = factor(
      pathway_short,
      levels = c("B cell receptor", "Chemokine signaling", "Fatty acid metabolism")
    ),
    set_type = factor(set_type, levels = c("Balanced (25/25)", "One-sided (25/0)"))
  )
theme_supp <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = "grey50"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Moran's I plot
p_moran <- ggplot(
  moran_plot_df,
  aes(x = pathway_short, y = z_moran, fill = pathway_short)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.6, width = 0.6, linewidth = 0.35) +
  coord_flip() +
  facet_wrap(~ set_type, ncol = 1) +
  labs(
    x = NULL,
    y = "Moran's I z-score",
    title = "Moran's I"
  ) +
  theme_supp
p_geary <- ggplot(
  moran_plot_df,
  aes(x = pathway_short, y = z_geary, fill = pathway_short)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.6, width = 0.6, linewidth = 0.35) +
  coord_flip() +
  facet_wrap(~ set_type, ncol = 1) +
  labs(
    x = NULL,
    y = "Geary's C z-score",
    title = "Geary's C"
  ) +
  theme_supp
comb_auto <- (p_moran | p_geary) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )
ggsave(
  filename = file.path(fig_dir, "supplement_Fig_S1_simulation_autocorrelation.pdf"),
  plot = comb_auto,
  width = 12,
  height = 7,
  device = cairo_pdf
)

# supplementary figure S2 for Permutation Null Calibration
qq_perm_vs_normal <- function(perm_stats, stat_name = "Moran's I", path_name) {
  m <- length(perm_stats)
  z <- scale(perm_stats)[,1]
  z_sorted <- sort(z)
  u <- (1:m)/(m+1)
  theo <- qnorm(u)
  # Binomial pointwise 95% band via Beta
  u_lo <- qbeta(0.025, 1:m, rev(1:m) + 1)
  u_hi <- qbeta(0.975, 1:m, rev(1:m) + 1)
  band <- data.frame(theo = theo, lo = qnorm(u_lo), hi = qnorm(u_hi))
  df   <- data.frame(theo = theo, samp = z_sorted)
  
  p <- ggplot() +
    geom_ribbon(data = band, aes(x = theo, ymin = lo, ymax = hi), fill = "grey85") +
    geom_abline(slope = 1, intercept = 0, color = "grey40") +
    geom_point(data = df, aes(x = theo, y = samp), size = 0.6, alpha = 0.7) +
    theme_bw() + labs(x = "Theoretical N(0,1) quantiles", y = "Permutation z-quantiles",
                      title = paste("QQ plot for", path_name, "-", stat_name))
  p
}

# get simulation results at chosen scenario (ell=200, sigma_f=0.5, 25/25)
scn <- sim_rep[[8]]$scn
sim <- sim_rep[[8]]$sim
expr <- as.matrix(sim$counts)

# build spatial weights
coords <- sim$coords
nb <- tri2nb(coords)
listw <- nb2listw(nb, style = "W")
islands <- which(card(nb)==0)
length(islands)

# get sgsva scores
sgsva_res <- run_sgsva(sim$counts, sim$coords, sim$per_path, 
                       sigma_global = floor(scn$ell/sqrt(2)), directionality = "oracle")
sgsva_scores <- sgsva_res$z_score
hid <- make_hex_ids(coords, cellsize = scn$ell)

mc_perm_null_res <- list()
hex_perm_null_res <- list()
for (p in seq_along(sel_paths)) {
  pathway <- sel_paths[p]
  x <- sgsva_scores[pathway,]
  
  set.seed(123)
  # get moran.mc results
  out_moran_mc <- perm_with_early_stop(stat_moran, x, listw)
  obs_mc <- out_moran_mc$obs
  # get hex-block results
  out_moran_hex <- perm_with_early_stop(stat_moran_hex, x, listw, hid)
  obs_hex <- out_moran_hex$obs
  
  # store permutation results
  mc_perm_null_res[[pathway]] <- list(
    moran_perm = out_moran_mc$perm_stats,
    obs = obs_mc
  )
  hex_perm_null_res[[pathway]] <- list(
    moran_perm = out_moran_hex$perm_stats,
    obc = obs_hex
  )
}

qq_plot <- list()
p_names <- c("chemokine", "Bcell", "fatty_acid")
for (p in seq_len(length(sel_paths))) {
  # moran.mc permutation
  moran_perm <- mc_perm_null_res[[p]]$moran_perm
  obs <- mc_perm_null_res[[p]]$obs
  qq_mc <- qq_perm_vs_normal(moran_perm, stat_name = "moran.mc", path_name = p_names[p])
  
  # hex-block permutation
  moran_perm <- hex_perm_null_res[[p]]$moran_perm
  obs <- hex_perm_null_res[[p]]$obc
  qq_hex <- qq_perm_vs_normal(moran_perm, stat_name = "hex block", path_name = p_names[p])
  
  qq_plot[[p]] <- list(qq_mc=qq_mc, qq_hex=qq_hex)
}
comb_qq <- (qq_plot[[1]]$qq_mc|qq_plot[[2]]$qq_mc|qq_plot[[3]]$qq_mc) /
  (qq_plot[[1]]$qq_hex|qq_plot[[2]]$qq_hex|qq_plot[[3]]$qq_hex)
ggsave(filename = file.path(fig_dir, "supplement_Fig_S2_simulation_perm_null_qqplot.pdf"),
      plot=comb_qq, width=12, height=15)

# supplementary figure S3 for Qualitative Spatial Recovery
path_name <- "Bcell"
path_id <- sel_paths[2]
set_type <- c("balanced", "one-sided")

map_list <- list()
for (set in set_type) {
  set.seed(1)
  if (set=="balanced") {
    sim <- simulate_counts(coords, gene_sets, sel_paths, ell=200, sigma_f=0.5, up=25, down=25)
  } else if(set=="one-sided") {
    sim <- simulate_counts(coords, gene_sets, sel_paths, ell=200, sigma_f=0.5, up=25, down=0)
  }
  
  # sGSVA
  sgsva_res <- run_sgsva(sim$counts, sim$coords, sim$per_path, 
                         directionality = "oracle", one_sided = "center_0.5")
  sgsva_scores <- sgsva_res$z_score
  
  # locCPM and filter the expression matrix
  counts <- sim$counts
  G <- nrow(counts); N <- ncol(counts)
  lib <- colSums(counts)+1
  lcp <- log1p(t(t(counts)/lib*1e6)) #logCPM
  expr <- as.matrix(lcp)
  keep <- rowSums(expr > 0) >= floor(0.05 * ncol(expr))
  expr <- expr[keep, , drop = FALSE]
  # GSVA
  gsva_scores <- run_gsva(expr, sim$per_path)
  z_gsva <- z_normalization(as.matrix(gsva_scores))
  
  df <- data.frame(
    x = sim$coords$x,
    y = sim$coords$y,
    sgsva = sgsva_scores[path_id, ],
    gsva = z_gsva[path_id, ],
    truth = sim$fields[[2]]
  )
  # --- unified scale for sGSVA + GSVA ---
  lim_methods <- max(abs(c(df$sgsva, df$gsva)), na.rm = TRUE)
  
  p1 <- ggplot(df, aes(x, y, color = sgsva)) +
    geom_point() + scale_color_viridis_c() + 
    scale_fill_gradient2(limits = c(-lim_methods, lim_methods)) +
    ggtitle(paste0("sGSVA (", set, ")")) + theme(legend.position = "none")
  p2 <- ggplot(df, aes(x, y, color = gsva)) +
    geom_point() + scale_color_viridis_c() + 
    scale_fill_gradient2(limits = c(-lim_methods, lim_methods)) +
    ggtitle(paste0("GSVA (", set, ")")) + labs(color = "scores")
  p3 <- ggplot(df, aes(x, y, color = truth)) +
    geom_point() + scale_color_viridis_c() + ggtitle(paste0("True field (", set, ")"))
  map_list[[set]] <- list(p1=p1, p2=p2, p3=p3)
}
comb_map <- (map_list$balanced$p1|map_list$balanced$p2|map_list$balanced$p3) /
  (map_list$`one-sided`$p1|map_list$`one-sided`$p2|map_list$`one-sided`$p3)
ggsave(filename = file.path(fig_dir, "supplement_Fig_S3_simulation_BCell_maps.pdf"),
       plot=comb_map, width=12, height=15)




