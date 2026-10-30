### DLPFC figures for final manuscript
### Uses section-level recurrence, donor coverage, and effect-size summaries.
### No Stouffer meta-analysis is used.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(patchwork)
  library(spatialLIBD)
  library(SpatialExperiment)
  library(stringr)
  library(scales)
  library(viridis)
  library(reshape2)
  library(here)
})

project_dir <- "."
results_dir <- file.path(project_dir, "results", "DLPFC")
fig_dir <- file.path(project_dir, "results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_dir, "R", "spargs_score.R"))

required_files <- c(
  file.path(results_dir, "full_results_12_samples.rds"),
  file.path(results_dir, "DLPFC_preprocess.rds"),
  file.path(results_dir, "DLPFC_section_level_hex_block.csv"),
  file.path(results_dir, "DLPFC_hex_block_recurrence_summary.csv"),
  file.path(results_dir, "DLPFC_top20_recurrent_pathways_hex200.csv"),
  file.path(results_dir, "DLPFC_blocksize_effect_size_summary.csv")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "\nMissing DLPFC source files required for figure generation:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    "\n\nRegenerate them with:\n",
    "  Rscript scripts/DLPFC/DLPFC_analysis.R\n",
    call. = FALSE
  )
}

# ---------------------------
# Load data and outputs
# ---------------------------
spe <- fetch_data(type = "spe")
sample_ids <- unique(as.character(spe$sample_id))

full_res <- readRDS(file.path(results_dir, "full_results_12_samples.rds"))
dlpfc_preprocessed <- readRDS(file.path(results_dir, "DLPFC_preprocess.rds"))
hex_section_stats <- read.csv(
  file.path(results_dir, "DLPFC_section_level_hex_block.csv"),
  stringsAsFactors = FALSE
)
hex_recurrence <- read.csv(
  file.path(results_dir, "DLPFC_hex_block_recurrence_summary.csv"),
  stringsAsFactors = FALSE
)
top20_hex200 <- read.csv(
  file.path(results_dir, "DLPFC_top20_recurrent_pathways_hex200.csv"),
  stringsAsFactors = FALSE
)

# Representative section
id <- "151673"
spa_res <- full_res[[id]]$spargs_res
gsva_scores <- full_res[[id]]$gsva_scores
z_gsva <- z_normalization(as.matrix(gsva_scores))
z_spa <- spa_res$z_score

# ---------------------------
# Helpers
# ---------------------------
pretty_kegg_name <- function(x) {
  x %>%
    stringr::str_replace("^KEGG_", "") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_title()
}

pretty_kegg_short <- function(x) {
  dplyr::case_when(
    x == "KEGG_CALCIUM_SIGNALING_PATHWAY" ~ "Calcium signaling",
    x == "KEGG_AXON_GUIDANCE" ~ "Axon guidance",
    x == "KEGG_LONG_TERM_POTENTIATION" ~ "Long-term potentiation",
    x == "KEGG_LONG_TERM_DEPRESSION" ~ "Long-term depression",
    x == "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION" ~ "Neuroactive ligand-receptor interaction",
    x == "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS" ~ "Fc gamma R-mediated phagocytosis",
    x == "KEGG_EPITHELIAL_CELL_SIGNALING_IN_HELICOBACTER_PYLORI_INFECTION" ~
      "Epithelial cell signaling in H. pylori infection",
    TRUE ~ pretty_kegg_name(x)
  )
}

get_spatial_map <- function(z_spa, z_gsva, coords, pathway) {
  stopifnot(pathway %in% rownames(z_spa))
  stopifnot(pathway %in% rownames(z_gsva))
  stopifnot(all(coords$spot %in% colnames(z_spa)))
  stopifnot(all(coords$spot %in% colnames(z_gsva)))

  coords <- coords %>%
    mutate(
      score_spa = as.numeric(z_spa[pathway, spot]),
      score_gsva = as.numeric(z_gsva[pathway, spot])
    )

  lim_max <- quantile(
    abs(c(coords$score_spa, coords$score_gsva)),
    probs = 0.98,
    na.rm = TRUE
  )
  lims <- c(-lim_max, lim_max)
  pathway_title <- pretty_kegg_short(pathway)

  base_map <- function(score_col, title) {
    ggplot(coords, aes(x = x, y = y, color = .data[[score_col]])) +
      geom_point(size = 1.6) +
      scale_color_gradientn(
        colors = c("#313695", "#ffffbf", "#a50026"),
        limits = lims,
        oob = scales::squish,
        name = "z-score"
      ) +
      scale_y_reverse() +
      coord_fixed() +
      theme_void() +
      ggtitle(title) +
      theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
      )
  }

  p1 <- base_map("score_spa", paste0(pathway_title, " - SpaRGS"))
  p2 <- base_map("score_gsva", paste0(pathway_title, " - GSVA"))
  p1 + p2 + patchwork::plot_layout(ncol = 2, guides = "collect")
}

plot_violin <- function(z_spa, z_gsva, coords, pathway) {
  stopifnot(pathway %in% rownames(z_spa))
  stopifnot(pathway %in% rownames(z_gsva))
  stopifnot(all(coords$spot %in% colnames(z_spa)))
  stopifnot(all(coords$spot %in% colnames(z_gsva)))

  layer_levels <- c("Layer1", "Layer2", "Layer3", "Layer4", "Layer5", "Layer6", "WM")

  coords_long <- coords %>%
    mutate(
      score_spa = as.numeric(z_spa[pathway, spot]),
      score_gsva = as.numeric(z_gsva[pathway, spot]),
      layer = factor(as.character(layer), levels = layer_levels)
    ) %>%
    filter(!is.na(layer)) %>%
    dplyr::select(layer, score_spa, score_gsva) %>%
    pivot_longer(
      cols = starts_with("score_"),
      names_to = "method",
      values_to = "score"
    ) %>%
    mutate(
      method = ifelse(method == "score_spa", "SpaRGS", "GSVA"),
      method = factor(method, levels = c("SpaRGS", "GSVA"))
    )

  ggplot(coords_long, aes(x = layer, y = score, fill = method)) +
    geom_violin(width = 0.9, position = position_dodge(width = 0.8), trim = TRUE) +
    geom_boxplot(
      width = 0.2,
      position = position_dodge(width = 0.8),
      outlier.size = 0.3
    ) +
    scale_fill_manual(values = c("SpaRGS" = "#377eb8", "GSVA" = "#e41a1c")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    ggtitle(paste0(pretty_kegg_short(pathway), " - layer score distributions")) +
    ylab("z-score") +
    xlab("Layer")
}

# rank correlation per sample
get_rank_cor <- function(df_sample, stat = "z_moran") {
  mat <- df_sample %>%
    dplyr::select(pathway, block_size_um, !!sym(stat)) %>%
    tidyr::pivot_wider(names_from = block_size_um, values_from = !!sym(stat), names_glue = "block{block_size_um}") %>%
    tibble::column_to_rownames("pathway")
  
  # rank by pathway within each column
  rank_mat <- apply(mat, 2, rank, ties.method = "average", na.last = "keep")
  cor(rank_mat, use = "pairwise.complete.obs", method = "spearman")
}

# Jaccard between two character vectors
.jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
# Build a square Jaccard matrix for one sample
jaccard_matrix <- function(df_sample, topk = 30, stat_col = "z_moran") {
  # Top-k pathway sets per block size
  sets <- df_sample %>%
    group_by(block_size_um) %>%
    slice_max(.data[[stat_col]], n = topk, with_ties = FALSE) %>%
    summarize(set = list(pathway), .groups = "drop")
  
  bs <- sets$block_size_um
  block_levels <- sort(unique(bs))
  n  <- length(bs)
  J  <- matrix(NA_real_, n, n, dimnames = list(bs, bs))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    J[i, j] <- .jac(sets$set[[i]], sets$set[[j]])
  }
  colnames(J) <- rownames(J) <- paste0("block", block_levels)
  J
}

extract_hex_component <- function(x) {
  if (!is.null(x$permutation_results$hex_block_results)) {
    return(x$permutation_results$hex_block_results)
  }
  if (!is.null(x$hex_block_results)) return(x$hex_block_results)
  stop("Could not locate hex_block_results in block permutation output.")
}

extract_mc_component <- function(x) {
  if (!is.null(x$permutation_results$moran_geary_results)) {
    return(x$permutation_results$moran_geary_results)
  }
  if (!is.null(x$moran_geary_results)) return(x$moran_geary_results)
  stop("Could not locate moran_geary_results in global permutation output.")
}

# ---------------------------
# Main Figure 5
# ---------------------------
coords <- as.data.frame(dlpfc_preprocessed[[id]]$coords_um)
colnames(coords)[1:2] <- c("x", "y")
coords$spot <- rownames(coords)

idx <- as.character(spe$sample_id) == id
layer_df <- data.frame(
  spot = colnames(spe)[idx],
  layer = spe$layer_guess_reordered[idx],
  stringsAsFactors = FALSE
)
coords <- coords %>% left_join(layer_df, by = "spot")

stopifnot(all(coords$spot %in% colnames(z_spa)))
stopifnot(all(coords$spot %in% colnames(z_gsva)))

cal_map <- get_spatial_map(
  z_spa, z_gsva, coords,
  pathway = "KEGG_CALCIUM_SIGNALING_PATHWAY"
)
cal_violin <- plot_violin(
  z_spa, z_gsva, coords,
  pathway = "KEGG_CALCIUM_SIGNALING_PATHWAY"
)
axon_map <- get_spatial_map(
  z_spa, z_gsva, coords,
  pathway = "KEGG_AXON_GUIDANCE"
)
axon_violin <- plot_violin(
  z_spa, z_gsva, coords,
  pathway = "KEGG_AXON_GUIDANCE"
)

# Panel G: leading recurrent pathways at 200 µm.
df_plot <- top20_hex200 %>%
  arrange(
    n_sections_sig_moran,
    n_donors_sig_moran,
    median_moran_I
  ) %>%
  mutate(
    pathway_lab = str_wrap(pretty_kegg_short(pathway), width = 28),
    repl_lab = paste0(
      n_sections_sig_moran, "/", n_sections_tested,
      "; ", n_donors_sig_moran, "/", n_donors_tested
    ),
    pathway_lab = factor(pathway_lab, levels = pathway_lab)
  )

p_hex200 <- ggplot(df_plot, aes(x = pathway_lab, y = median_moran_I)) +
  geom_col(fill = "grey40", width = 0.75) +
  geom_text(
    aes(label = repl_lab),
    hjust = -0.10,
    color = "firebrick",
    size = 3.1
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    name = "Median section-level Moran's I",
    limits = c(0, max(df_plot$median_moran_I, na.rm = TRUE) + 0.12),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Pathway",
    title = "Recurrent spatial pathways across 12 DLPFC sections from three donors (200-µm hex-block)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.margin = margin(5.5, 65, 5.5, 5.5),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 11)
  )

# Panels H-I: block-size stability
R_list <- list()
J_list <- list()
topks <- c(10, 20, 30)
for (i in 1:length(sample_ids)) {
  id <- sample_ids[i]
  combined_df <- read.csv(file.path(results_dir, paste0("combined_hex_block_summary_",id,".csv")))
  
  # rank correlation heatmap
  cor_mat <- get_rank_cor(combined_df)
  R_list[[id]] <- cor_mat
  
  # Jaccard heatmap for top 10/20/30
  for (k in topks) {
    jaccard_topk <- jaccard_matrix(combined_df, topk = k)
    J_list[[paste0("top",k)]][[id]] <- jaccard_topk
  }
}

# aggregate mean across samples
mean_matrix <- function(mats) {
  Reduce(`+`, mats) / length(mats)
}
J_top30_mean <- mean_matrix(J_list$top30)
df <- reshape2::melt(J_top30_mean)
colnames(df) <- c("Row", "Col", "Value")
p_J <- ggplot(df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile(color = "grey70") +
  geom_text(aes(label = sprintf("%.2f", Value), color = Value > 0.6), size = 5.2, fontface = "bold") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_viridis_c(limits = c(0, 1), name = "Value") +
  scale_y_discrete(limits = rev) +   # flips to match pheatmap orientation
  labs(title = "Top30 Jaccard overlap (block-size sensitivity)", x = "", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )

R_mean <- mean_matrix(R_list)
df <- reshape2::melt(R_mean)
colnames(df) <- c("Row", "Col", "Value")
p_R <- ggplot(df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile(color = "grey70") +
  geom_text(aes(label = sprintf("%.2f", Value), color = Value > 0.6), size = 5.2, fontface = "bold") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_viridis_c(limits = c(0, 1), name = "Value") +
  scale_y_discrete(limits = rev) +   # flips to match pheatmap orientation
  labs(title = "Rank correlation (block-size sensitivity)", x = "", y = "") +
  theme_minimal(base_size = 14)  +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )

comb_fig <-
  (cal_map | cal_violin) /
  (axon_map | axon_violin) /
  (p_hex200 | ((p_J + theme(legend.position = "none")) | p_R)) +
  plot_annotation(tag_levels = "A") +
  plot_layout(
    heights = c(1, 1, 1.15),
    widths = c(1.35, 1),
    guides = "collect"
  )

ggsave(
  filename = file.path(fig_dir, "Fig5_DLPFC_main.pdf"),
  plot = comb_fig,
  width = 25,
  height = 15
)

# -----------------------------
# Supplementary Figure S5
# Additional DLPFC pathways in section 151673
# -----------------------------
oxidat_map <- get_spatial_map(
  z_spa, z_gsva, coords,
  pathway = "KEGG_OXIDATIVE_PHOSPHORYLATION"
)
chemo_map <- get_spatial_map(
  z_spa, z_gsva, coords,
  pathway = "KEGG_CHEMOKINE_SIGNALING_PATHWAY"
)
oxidat_violin <- plot_violin(
  z_spa, z_gsva, coords,
  pathway = "KEGG_OXIDATIVE_PHOSPHORYLATION"
)
chemo_violin <- plot_violin(
  z_spa, z_gsva, coords,
  pathway = "KEGG_CHEMOKINE_SIGNALING_PATHWAY"
)

comb_map <-
  (oxidat_map | oxidat_violin) /
  (chemo_map | chemo_violin) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1, 0.9), guides = "collect")

ggsave(
  filename = file.path(fig_dir, "supplement_Fig_S5_DLPFC_151673_map_violin.pdf"),
  plot = comb_map,
  width = 20,
  height = 10
)

# -----------------------------
# Supplementary Figure S6
# Calcium signaling across all 12 sections from three donors
# -----------------------------
calcium_pathway <- "KEGG_CALCIUM_SIGNALING_PATHWAY"
all_calcium_scores <- unlist(lapply(sample_ids, function(section_id) {
  z <- full_res[[section_id]]$sgsva_res$z_score
  if (!calcium_pathway %in% rownames(z)) return(numeric(0))
  as.numeric(z[calcium_pathway, ])
}))
cal_lim <- quantile(abs(all_calcium_scores), 0.98, na.rm = TRUE)

map_list <- lapply(sample_ids, function(section_id) {
  z <- full_res[[section_id]]$sgsva_res$z_score
  stopifnot(calcium_pathway %in% rownames(z))

  coords_i <- as.data.frame(dlpfc_preprocessed[[section_id]]$coords_um)
  colnames(coords_i)[1:2] <- c("x", "y")
  coords_i$spot <- rownames(coords_i)
  stopifnot(all(coords_i$spot %in% colnames(z)))
  coords_i$score_sgsva <- as.numeric(z[calcium_pathway, coords_i$spot])

  donor_id <- unique(
    hex_section_stats$donor_id[hex_section_stats$sample_id == section_id]
  )
  donor_id <- donor_id[!is.na(donor_id)][1]

  ggplot(coords_i, aes(x, y, color = score_sgsva)) +
    geom_point(size = 1.5) +
    scale_color_gradientn(
      colors = c("#313695", "#ffffbf", "#a50026"),
      limits = c(-cal_lim, cal_lim),
      oob = scales::squish,
      name = "z-score"
    ) +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    ggtitle(paste0(section_id, " (", donor_id, ")")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
})

all_map <- wrap_plots(map_list, ncol = 4, guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = file.path(fig_dir, "supplement_Fig_S6_DLPFC_CALCIUM_maps.pdf"),
  plot = all_map,
  width = 20,
  height = 15
)

# -----------------------------
# Supplementary Figure S7
# Null distributions in representative section 151673
# -----------------------------
id <- "151673"
block_res <- readRDS(file.path(results_dir, paste0(id, "_block_res.rds")))
global_res <- readRDS(file.path(results_dir, paste0(id, "_moranmc_res.rds")))

top4_moran <- c(
  "KEGG_CALCIUM_SIGNALING_PATHWAY",
  "KEGG_AXON_GUIDANCE",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY"
)

make_hex_hist_list <- function(block_obj, block_size_um, pathways) {
  comp <- extract_hex_component(block_obj)
  summary_df <- comp$hex_block_summary
  null_res <- comp$hex_perm_null_results

  out <- list()
  for (nm in pathways) {
    df <- data.frame(val = null_res[[nm]]$moran_perm)
    observed <- summary_df$obs_moran[summary_df$pathway == nm]
    p <- ggplot(df, aes(x = val)) +
      geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
      geom_vline(
        xintercept = observed,
        color = "firebrick",
        linetype = 2,
        linewidth = 1
      ) +
      theme_bw() +
      labs(
        x = "Moran's I (permuted)",
        y = "Count",
        title = paste0(tolower(gsub("_", " ", sub("^KEGG_", "", nm))),
                       " (hex-block ", block_size_um, " µm)")
      )
    out[[nm]] <- p
  }
  out
}

hist_block150 <- make_hex_hist_list(block_res[[1]], 150, top4_moran)
hist_block200 <- make_hex_hist_list(block_res[[2]], 200, top4_moran)
hist_block300 <- make_hex_hist_list(block_res[[3]], 300, top4_moran)

mc_comp <- extract_mc_component(global_res)
mc_results_df <- mc_comp$moran_geary_summary
mc_perm_null_results <- mc_comp$mc_perm_null_results
hist_mc <- list()
for (nm in top4_moran) {
  df <- data.frame(val = mc_perm_null_results[[nm]]$moran_perm)
  observed <- mc_results_df$obs_moran[mc_results_df$pathway == nm]
  hist_mc[[nm]] <- ggplot(df, aes(x = val)) +
    geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
    geom_vline(
      xintercept = observed,
      color = "firebrick",
      linetype = 2,
      linewidth = 1
    ) +
    theme_bw() +
    labs(
      x = "Moran's I (permuted)",
      y = "Count",
      title = paste0(tolower(gsub("_", " ", sub("^KEGG_", "", nm))), " (global permutation)")
    )
}

comb_hist <-
  wrap_plots(hist_block150[top4_moran], ncol = 4) /
  wrap_plots(hist_block200[top4_moran], ncol = 4) /
  wrap_plots(hist_block300[top4_moran], ncol = 4) /
  wrap_plots(hist_mc[top4_moran], ncol = 4) +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = file.path(fig_dir, "supplement_Fig_S7_DLPFC_151673_null_hist.pdf"),
  plot = comb_hist,
  width = 20,
  height = 20
)
