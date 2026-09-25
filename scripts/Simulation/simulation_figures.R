# making simulation figures
library(patchwork)
library(dplyr)
library(ggplot2)
library(scales)

# ---------------------------
# paths and reproducibility
# ---------------------------
project_dir <- "."
results_dir <- file.path(project_dir, "results", "simulation")
fig_dir <- file.path(project_dir, "results", "figures")

source(file.path(project_dir,"R/spatial_permutation.R"))
source(file.path(project_dir,"R/simulation.R"))
source(file.path(project_dir,"R/simulation_scoring.R"))
source(file.path(project_dir,"R/evaluation_metrics.R"))

# ---------------------------
# load saved outputs
# ---------------------------
region_raw <- read.csv(file.path(results_dir, "evaluate_metrics_all_replicates_data_driven_plus_gsdensity.csv"))
region_raw$X <- NULL

sim_res <- readRDS(file.path(results_dir, "rep_simulation.rds"))

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
  "sGSVA",
  "GSVA",
  "GSVA-smoothed",
  "GSDensity",
  "ssGSEA",
  "AUCell",
  "UCell",
  "AddModuleScore"
)

# create plotting data
retained_keys <- region_raw %>%
  filter(methods == "sGSVA", retained == TRUE, !is.na(AUROC)) %>%
  distinct(pathway, replicate, ell, sigma_f, up_size, down_size)
region_retained_matched <- region_raw %>%
  semi_join(
    retained_keys,
    by = c("pathway", "replicate", "ell", "sigma_f", "up_size", "down_size")
  )

spatial_methods_plot <- region_retained_matched %>%
  mutate(
    methods = recode(methods,
                     "GSVA_smoothed" = "GSVA-smoothed"),
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0  ~ "One-sided (25/0)",
      TRUE ~ NA_character_
    ),
    set_type = factor(set_type, levels = c("Balanced (25/25)", "One-sided (25/0)")),
    methods = factor(methods, levels = method_levels)
  )

# AUROC plot
p_auc <- ggplot(spatial_methods_plot, aes(x = methods, y = AUROC, fill = methods)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.35) +
  facet_wrap(~ set_type, nrow = 1) +
  labs(
    x = NULL,
    y = "Orientation-invariant AUROC",
    title = "Region-level recovery"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# spearman plot
p_abs_spearman <- ggplot(spatial_methods_plot, aes(x = methods, y = spearman, fill = methods)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.35) +
  facet_wrap(~ set_type, nrow = 1) +
  labs(
    x = NULL,
    y = "Absolute Spearman correlation",
    title = "Map-to-truth recovery"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# retention heatmap
retention_sum <- region_raw %>%
  filter(methods == "sGSVA") %>%
  mutate(
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0  ~ "One-sided (25/0)",
      TRUE ~ "Other"
    ),
    set_type = factor(set_type, levels = c("Balanced (25/25)", "One-sided (25/0)")),
    ell = factor(ell, levels = c(50, 200, 400)),
    sigma_f = factor(sigma_f, levels = c(0.5, 1.0, 1.5))
  ) %>%
  group_by(set_type, ell, sigma_f) %>%
  summarise(
    n_total = n(),
    n_retained = sum(retained, na.rm = TRUE),
    retention_rate = n_retained / n_total,
    .groups = "drop"
  ) %>%
  mutate(
    retention_label = ifelse(retention_rate == 0 | retention_rate == 1,
                             paste0(round(100 * retention_rate), "%"),
                             paste0(round(100 * retention_rate, 1), "%")),
    label_color = ifelse(retention_rate < 0.55, "white", "black")
  )

p_retention <- ggplot(retention_sum, aes(x = ell, y = sigma_f, fill = retention_rate)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = retention_label, color = label_color),
            size = 6, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ set_type, nrow = 1) +
  scale_color_identity() +
  scale_fill_gradientn(
    colours = c("#0b1f3a", "#355f8c", "#5aa3e1", "#dceefb"),
    values = scales::rescale(c(0, 0.35, 0.7, 1)),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = percent_format(accuracy = 1),
    name = "Retention"
  ) +
  labs(
    x = expression("Spatial scale " * ell * " (" * mu * "m)"),
    y = expression("Signal amplitude " * sigma[f]),
    title = "Data-driven sGSVA retention"
  )

fig_simulation_main <- (p_auc | p_abs_spearman) / p_retention +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  file.path(fig_dir, "Fig2_simulation_performance.pdf"),
  fig_simulation_main,
  width = 12,
  height = 8
)
