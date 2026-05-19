# making DLPFC figures
suppressPackageStartupMessages ({
  library(ggplot2); library(dplyr); library(purrr);
  library(patchwork); library(gridExtra); library(spatialLIBD); 
  library(SpatialExperiment); library(stringr); library(scales); 
  library(stringr); library(viridis)
})
source("R/sgsva_score.R")

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
results_dir <- file.path(project_dir, "results", "DLPFC")
fig_dir <- file.path(project_dir, "results", "figures")

# ---------------------------
# load saved outputs
# ---------------------------
# load DLPFC dataset
spe <- fetch_data(type="spe")
sample_ids <- unique(spe$sample_id)
spe_list <- lapply(sample_ids, function(id) spe[, spe$sample_id == id])
names(spe_list) <- sample_ids

# load full score results
full_res <- readRDS(file.path(results_dir, "full_results_12_samples.rds"))

# load meta analysis results
hex150_meta <- read.csv(file.path(results_dir, "hex150_meta_stouffer_res.csv"))
hex150_meta$X <- NULL
hex200_meta <- read.csv(file.path(results_dir, "hex200_meta_stouffer_res.csv"))
hex200_meta$X <- NULL
hex300_meta <- read.csv(file.path(results_dir, "hex300_meta_stouffer_res.csv"))
hex300_meta$X <- NULL
mc_meta <- read.csv(file.path(results_dir, "mc_meta_stouffer_res.csv"))
mc_meta$X <- NULL

# select a representative sample
id <- "151673"
# get GSVA, sGSVA scores
sgsva_res <- full_res[[id]]$sgsva_res
gsva_scores <- full_res[[id]]$gsva_scores
z_gsva <- z_normalization(as.matrix(gsva_scores))
z_sgsva <- sgsva_res$z_score
balanced_sets <- sgsva_res$gene_sets

# ---------------------------
# helper functions
# ---------------------------
pretty_kegg_name <- function(x) {
  x %>%
    str_replace("^KEGG_", "") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}
pretty_kegg_short <- function(x) {
  case_when(
    x == "KEGG_CALCIUM_SIGNALING_PATHWAY" ~ "Calcium signaling",
    x == "KEGG_AXON_GUIDANCE" ~ "Axon guidance",
    x == "KEGG_LONG_TERM_POTENTIATION" ~ "Long-term potentiation",
    x == "KEGG_LONG_TERM_DEPRESSION" ~ "Long-term depression",
    x == "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION" ~ "Neuroactive ligand–receptor interaction",
    x == "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS" ~ "Fc gamma R-mediated phagocytosis",
    x == "KEGG_EPITHELIAL_CELL_SIGNALING_IN_HELICOBACTER_PYLORI_INFECTION" ~ 
      "Epithelial cell signaling in H. pylori infection",
    TRUE ~ pretty_kegg_name(x)
  )
}

# generate spatial map
get_spatial_map <- function(z_sgsva, z_gsva, coords, pathway) {
  coords$score_sgsva <- as.numeric(z_sgsva[pathway, ])
  coords$score_gsva <- as.numeric(z_gsva[pathway, ])
  lims <- quantile(c(coords$score_sgsva, coords$score_gsva), probs = c(0.02, 0.98), na.rm = TRUE)
  pathway_title <- pretty_kegg_short(pathway)
  
  p1 <- ggplot(coords, aes(x, y, color = pmin(pmax(score_sgsva, lims[1]), lims[2]))) +
    geom_point(size = 1.6) +
    viridis::scale_color_viridis(option = "C", limits = lims, name = "sGSVA") +
    scale_y_reverse() + coord_fixed() + theme_void() +
    ggtitle(paste0(pathway_title, "-sGSVA"))
  
  p2 <- ggplot(coords, aes(x, y, color = pmin(pmax(score_gsva, lims[1]), lims[2]))) +
    geom_point(size = 1.6) +
    viridis::scale_color_viridis(option = "C", limits = lims, name = "GSVA") +
    scale_y_reverse() + coord_fixed() + theme_void() +
    ggtitle(paste0(pathway_title, "-GSVA"))
  p1 + p2 + plot_layout(ncol = 2)
}

# violin plots by layer
plot_violin <- function(z_sgsva, z_gsva, coords, pathway) {
  coords$score_sgsva <- as.numeric(z_sgsva[pathway,])
  coords$score_gsva <- as.numeric(z_gsva[pathway,])
  coords_long <- coords %>%
    dplyr::select(layer, score_sgsva, score_gsva) %>%
    tidyr::pivot_longer(cols = starts_with("score_"),
                        names_to = "method",
                        values_to = "score") %>%
    mutate(method = ifelse(method == "score_sgsva", "sGSVA", "GSVA"),
           layer = factor(layer, levels = sort(unique(layer))))
  
  pathway_title <- pretty_kegg_short(pathway)
  p <- ggplot(coords_long, aes(x = layer, y = score, fill = method)) +
    geom_violin(width = 0.9, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), outlier.size = 0.3) +
    scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          panel.grid.major.x = element_blank()) +
    ggtitle(paste0(pathway_title, " – Layer score distributions")) +
    ylab("Score") + xlab("Layer")
  return(p)
}

# rank correlation per sample
get_rank_cor <- function(df_sample, stat = "z_moran") {
  mat <- df_sample %>%
    dplyr::select(pathway, block_size, !!sym(stat)) %>%
    tidyr::pivot_wider(names_from = block_size, values_from = !!sym(stat), names_glue = "block{block_size}") %>%
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
    group_by(block_size) %>%
    slice_max(.data[[stat_col]], n = topk, with_ties = FALSE) %>%
    summarize(set = list(pathway), .groups = "drop")
  
  bs <- sets$block_size
  block_levels <- sort(unique(bs))
  n  <- length(bs)
  J  <- matrix(NA_real_, n, n, dimnames = list(bs, bs))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    J[i, j] <- .jac(sets$set[[i]], sets$set[[j]])
  }
  colnames(J) <- rownames(J) <- paste0("block", block_levels)
  J
}


# -----------------------------
# DLPFC main figures
# -----------------------------
# generate spatial maps/violin plots
id <- "151673"
idx <- spe$sample_id == id
coords <- as.data.frame(spatialCoords(spe)[idx, ])
colnames(coords)[1:2] <- c("x","y")   # usually fullres pixel columns/rows
coords$spot <- colnames(spe)[idx]
coords$layer <- spe$layer_guess_reordered[idx]

cal_map <- get_spatial_map(z_sgsva, z_gsva, coords, pathway="KEGG_CALCIUM_SIGNALING_PATHWAY")
cal_violin <- plot_violin(z_sgsva, z_gsva, coords, pathway="KEGG_CALCIUM_SIGNALING_PATHWAY")
axon_map <- get_spatial_map(z_sgsva, z_gsva, coords, pathway="KEGG_AXON_GUIDANCE")
axon_violin <- plot_violin(z_sgsva, z_gsva, coords, pathway="KEGG_AXON_GUIDANCE")

# meta-analysis bar plot
df_plot <- hex200_meta %>%
  arrange(desc(n_sig), desc(median_obs_moran), median_obs_geary) %>%
  slice_head(n = 20) %>%
  mutate(
    pathway_lab = pretty_kegg_short(pathway),
    pathway_lab = str_wrap(pathway_lab, width = 28),
    repl_lab = paste0(n_sig, "/", k),
    pathway_lab = factor(pathway_lab, levels = rev(pathway_lab))
  )
p_hex200 <- ggplot(
  df_plot,
  aes(x = pathway_lab, y = median_obs_moran)
) +
  geom_col(fill = "grey40", width = 0.75) +
  geom_text(
    aes(label = repl_lab),
    hjust = -0.15,
    color = "red",
    size = 3.3
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    name = "Median section-level Moran's I",
    limits = c(0, max(df_plot$median_obs_moran, na.rm = TRUE) + 0.08),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Pathway",
    title = "Reproducible spatial pathways across 12 DLPFC sections (hex-block 200)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.margin = margin(5.5, 45, 5.5, 5.5),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )


# cross-section stability
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

# combine into one figure
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
# DLPFC supplementary figures
# -----------------------------
# supplementary figure S1 for 151673
oxidat_map <- get_spatial_map(z_sgsva, z_gsva, coords, pathway="KEGG_OXIDATIVE_PHOSPHORYLATION")
chemo_map <- get_spatial_map(z_sgsva, z_gsva, coords, pathway="KEGG_CHEMOKINE_SIGNALING_PATHWAY")
oxidat_violin <- plot_violin(z_sgsva, z_gsva, coords, pathway="KEGG_OXIDATIVE_PHOSPHORYLATION")
chemo_violin <- plot_violin(z_sgsva, z_gsva, coords, pathway="KEGG_CHEMOKINE_SIGNALING_PATHWAY")
comb_map <- (oxidat_map | oxidat_violin) / (chemo_map | oxidat_violin) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1, 0.9), guides = "collect")
ggsave(
  filename = file.path(fig_dir, "supplement_Fig_S8_DLPFC_151673_map_violin.pdf"),
  plot = comb_map, width = 20, height = 10
)

# supplementary figure S2 for all CALCIUM map
map_list <- list()
for (i in 1:length(sample_ids)) {
  id <- sample_ids[i]
  # get sGSVA scores
  sgsva_res <- full_res[[id]]$sgsva_res
  z_sgsva <- sgsva_res$z_score
  balanced_sets <- sgsva_res$gene_sets
  
  # generate spatial maps
  idx <- spe$sample_id == id
  coords <- as.data.frame(spatialCoords(spe)[idx, ])
  colnames(coords)[1:2] <- c("x","y")   # usually fullres pixel columns/rows
  coords$spot <- colnames(spe)[idx]
  coords$layer <- spe$layer_guess_reordered[idx]
  coords$score_sgsva <- as.numeric(z_sgsva["KEGG_CALCIUM_SIGNALING_PATHWAY", ])
  
  p <- ggplot(coords, aes(x, y, color = score_sgsva)) +
    geom_point(size = 1.6) +
    scale_color_viridis(option = "C", name = "sGSVA") +
    scale_y_reverse() + coord_fixed() + theme_void() +
    ggtitle(paste0(id, "-sGSVA"))
  map_list[[i]] <- p
}
all_map <- ((map_list[[1]] + theme(legend.position="none"))|(map_list[[2]] + theme(legend.position="none"))
            |(map_list[[3]] + theme(legend.position="none"))|(map_list[[4]] + theme(legend.position="none"))) /
  ((map_list[[5]] + theme(legend.position="none"))|(map_list[[6]] + theme(legend.position="none"))
   |(map_list[[7]] + theme(legend.position="none"))|(map_list[[8]] + theme(legend.position="none"))) /
  ((map_list[[9]] + theme(legend.position="none"))|(map_list[[10]] + theme(legend.position="none"))|
     (map_list[[11]] + theme(legend.position="none"))|map_list[[12]]) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1, 0.9), guides = "collect")
ggsave(
  filename = file.path(fig_dir, "supplement_Fig_S9_DLPFC_CALCIUM_maps.pdf"),
  plot = all_map, width = 20, height = 15
)

# supplementary figure S3 for 151673 null histograms
id <- "151673"
block_res <- readRDS(file.path(results_dir,paste0(id, "_block_res.rds")))
top4_moran <- c("KEGG_CALCIUM_SIGNALING_PATHWAY", "KEGG_AXON_GUIDANCE",
                "KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_CHEMOKINE_SIGNALING_PATHWAY")
hex_results_df <- block_res[[1]]$hex_block_results$hex_block_summary
hex_perm_null_results <- block_res[[1]]$hex_block_results$hex_perm_null_results
hist_block150 <- list()
for (nm in top4_moran) {
  df <- data.frame(val = hex_perm_null_results[[nm]]$moran_perm)
  path_name <- stringr::str_replace_all(nm, "^KEGG_", "")
  path_name <- tolower(gsub("_", " ", path_name))
  p <- ggplot(df, aes(x = val)) +
    geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
    geom_vline(xintercept = hex_results_df$obs_moran[hex_results_df$pathway==nm], 
               color = "firebrick", linetype = 2, linewidth = 1) +
    theme_bw() +
    labs(x = "Moran's I (permuted)", y = "Count",
         title = paste0(path_name, " (hex-block 150)"))
  hist_block150[[nm]] <- p
}

hex_results_df <- block_res[[2]]$hex_block_results$hex_block_summary
hex_perm_null_results <- block_res[[2]]$hex_block_results$hex_perm_null_results
hist_block200 <- list()
for (nm in top4_moran) {
  df <- data.frame(val = hex_perm_null_results[[nm]]$moran_perm)
  path_name <- stringr::str_replace_all(nm, "^KEGG_", "")
  path_name <- tolower(gsub("_", " ", path_name))
  p <- ggplot(df, aes(x = val)) +
    geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
    geom_vline(xintercept = hex_results_df$obs_moran[hex_results_df$pathway==nm], 
               color = "firebrick", linetype = 2, linewidth = 1) +
    theme_bw() +
    labs(x = "Moran's I (permuted)", y = "Count",
         title = paste0(path_name, " (hex-block 200)"))
  hist_block200[[nm]] <- p
}

hex_results_df <- block_res[[3]]$hex_block_results$hex_block_summary
hex_perm_null_results <- block_res[[3]]$hex_block_results$hex_perm_null_results
hist_block300 <- list()
for (nm in top4_moran) {
  df <- data.frame(val = hex_perm_null_results[[nm]]$moran_perm)
  path_name <- stringr::str_replace_all(nm, "^KEGG_", "")
  path_name <- tolower(gsub("_", " ", path_name))
  p <- ggplot(df, aes(x = val)) +
    geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
    geom_vline(xintercept = hex_results_df$obs_moran[hex_results_df$pathway==nm], 
               color = "firebrick", linetype = 2, linewidth = 1) +
    theme_bw() +
    labs(x = "Moran's I (permuted)", y = "Count",
         title = paste0(path_name, " (hex-block 300)"))
  hist_block300[[nm]] <- p
}

sgsva_res <- readRDS(file.path(results_dir, paste0(id, "_moranmc_res.rds")))
mc_results_df <- sgsva_res$moran_geary_results$moran_geary_summary
mc_perm_null_results <- sgsva_res$moran_geary_results$mc_perm_null_results
hist_mc <- list()
for (nm in top4_moran) {
  df <- data.frame(val = mc_perm_null_results[[nm]]$moran_perm)
  path_name <- stringr::str_replace_all(nm, "^KEGG_", "")
  path_name <- tolower(gsub("_", " ", path_name))
  p <- ggplot(df, aes(x = val)) +
    geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
    geom_vline(xintercept = mc_results_df$obs_moran[mc_results_df$pathway==nm], 
               color = "firebrick", linetype = 2, linewidth = 1) +
    theme_bw() +
    labs(x = "Moran's I (permuted)", y = "Count",
         title = paste0(path_name, " (moran.mc)"))
  hist_mc[[nm]] <- p
}

comb_hist <- (hist_block150$KEGG_CALCIUM_SIGNALING_PATHWAY|hist_block150$KEGG_AXON_GUIDANCE|
                hist_block150$KEGG_OXIDATIVE_PHOSPHORYLATION|hist_block150$KEGG_CHEMOKINE_SIGNALING_PATHWAY) /
  (hist_block200$KEGG_CALCIUM_SIGNALING_PATHWAY|hist_block200$KEGG_AXON_GUIDANCE|
     hist_block200$KEGG_OXIDATIVE_PHOSPHORYLATION|hist_block200$KEGG_CHEMOKINE_SIGNALING_PATHWAY) /
  (hist_block300$KEGG_CALCIUM_SIGNALING_PATHWAY|hist_block300$KEGG_AXON_GUIDANCE|
     hist_block300$KEGG_OXIDATIVE_PHOSPHORYLATION|hist_block300$KEGG_CHEMOKINE_SIGNALING_PATHWAY) /
  (hist_mc$KEGG_CALCIUM_SIGNALING_PATHWAY|hist_mc$KEGG_AXON_GUIDANCE|
     hist_mc$KEGG_OXIDATIVE_PHOSPHORYLATION|hist_mc$KEGG_CHEMOKINE_SIGNALING_PATHWAY) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1, 0.9), guides = "collect")
ggsave(
  filename = file.path(fig_dir, "supplement_Fig_S10_DLPFC_151673_null_hist.pdf"),
  plot = comb_hist, width = 30, height = 30
)

