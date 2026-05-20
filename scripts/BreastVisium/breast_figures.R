# making breast Visium figures
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2)
  library(patchwork); library(gridExtra);library(ggspatial)
  library(ggpubr)
})

source("R/sgsva_score.R")
source("R/spatial_permutation.R")


# Required source files for Breast Visium figure generation
required_files <- c(
  "results/breast_visium/bc_obj.rds",
  "results/breast_visium/block_res.rds",
  "results/breast_visium/moranmc_res.rds",
  "results/breast_visium/sgsva_scores.csv",
  "results/breast_visium/gsva_scores.csv"
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "\nMissing Breast Visium source files required for figure generation:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    "\n\nThese large intermediate files are not included in the lightweight GitHub repository.\n",
    "To regenerate them, run:\n\n",
    "  Rscript scripts/BreastVisium/breast_analysis.R\n",
    "  Rscript scripts/BreastVisium/breast_celltype_anchoring.R\n",
    "After these files are generated, rerun:\n\n",
    "  Rscript scripts/BreastVisium/breast_figures.R\n\n",
    "Alternatively, use the archived full analysis outputs if provided separately.\n",
    call. = FALSE
  )
}

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
results_dir <- file.path(project_dir, "results", "breast_visium")
fig_dir     <- file.path(project_dir, "results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# load saved outputs
# ---------------------------
block_res <- readRDS(file.path(results_dir, "block_res.rds"))
hex_results_df <- block_res[[2]]$permutation_results$hex_block_results$hex_block_summary
hex_perm_null_results <- block_res[[2]]$permutation_results$hex_block_results$hex_perm_null_results
sgsva_res <- readRDS(file.path(results_dir, "moranmc_res.rds"))

anchor_merge_dat <- read.csv(file.path(results_dir, "anchoring_sgsva_scores_merged.csv"))

bc_obj <- readRDS(file.path(results_dir, "bc_obj.rds"))
sgsva_scores <- read.csv(file.path(results_dir, "sgsva_scores.csv"), row.names = 1, check.names = FALSE)
z_sgsva <- z_normalization(as.matrix(sgsva_scores))
gsva_scores <- read.csv(file.path(results_dir, "gsva_scores.csv"), row.names = 1, check.names = FALSE)
z_gsva <- z_normalization(as.matrix(gsva_scores))

# Convert pathway scores to a matrix (sGSVA)
gsva_smooth_matrix <- as.matrix(z_sgsva)
rownames(gsva_smooth_matrix) <- gsub("_", "-", rownames(gsva_smooth_matrix))
colnames(gsva_smooth_matrix) <- gsub(".", "-", colnames(gsva_smooth_matrix), fixed = TRUE)

# Create a new assay in the Seurat object
bc_obj[["sGSVA"]] <- CreateAssayObject(data = gsva_smooth_matrix)
DefaultAssay(bc_obj) <- "sGSVA"
# Set an empty meta.features dataframe with the same number of rows as pathways
bc_obj[["sGSVA"]]@meta.features <- data.frame(
  row.names = rownames(gsva_smooth_matrix)
)

# Convert pathway scores to a matrix (GSVA)
gsva_score_matrix <- as.matrix(z_gsva)
rownames(gsva_score_matrix) <- gsub("_", "-", rownames(gsva_score_matrix))
colnames(gsva_score_matrix) <- gsub(".", "-", colnames(gsva_score_matrix), fixed = TRUE)

# Create a new assay in the Seurat object
bc_obj[["GSVA"]] <- CreateAssayObject(data = gsva_score_matrix)
DefaultAssay(bc_obj) <- "GSVA"
# Set an empty meta.features dataframe with the same number of rows as pathways
bc_obj[["GSVA"]]@meta.features <- data.frame(
  row.names = rownames(gsva_score_matrix)
)


# ---------------------------
# helper functions
# ---------------------------
pretty_pathway_short <- function(pathway) {
  dplyr::case_when(
    grepl("B-CELL-RECEPTOR", pathway) ~ "BCR signaling",
    grepl("T-CELL-RECEPTOR", pathway) ~ "TCR signaling",
    grepl("GLUTATHIONE", pathway) ~ "Glutathione metabolism",
    grepl("FATTY-ACID", pathway) ~ "Fatty acid metabolism",
    TRUE ~ gsub("-", " ", pathway)
  )
}
# violin plots by annotation
plot_violin <- function(spatial_obj, pathway, annotation_col = "Classification") {
  
  # Fetch pathway scores from both assays
  DefaultAssay(spatial_obj) <- "sGSVA"
  sgsva_scores <- FetchData(spatial_obj, vars = pathway, layer = "data")
  DefaultAssay(spatial_obj) <- "GSVA"
  gsva_scores  <- FetchData(spatial_obj, vars = pathway, layer = "data")
  
  # Add annotation column
  ann <- spatial_obj@meta.data[[annotation_col]]
  
  # Combine into one data frame
  df <- data.frame(
    Annotation = ann,
    sGSVA = as.numeric(sgsva_scores[[1]]),
    GSVA  = as.numeric(gsva_scores[[1]])
  ) %>%
    tidyr::pivot_longer(cols = c("sGSVA", "GSVA"), names_to = "Method", values_to = "Score")
  
  # --- biologically defined grouping for main figure ---
  df <- df %>%
    mutate(
      Annotation_group = dplyr::case_when(
        Annotation == "Lymphocytes" ~ "Lymphocytes",
        Annotation == "Stroma" ~ "Stroma",
        Annotation %in% c(
          "Invasive cancer + stroma + lymphocytes",
          "Normal glands + lymphocytes",
          "DCIS",
          "Cancer trapped in lymphocyte aggregation"
        ) ~ "Epithelial/tumor",
        Annotation %in% c("Adipose tissue", "Artefact") ~ "Other",
        TRUE ~ "Other"
      )
    )
  
  # order by size
  df$Annotation_group <- factor(
    df$Annotation_group,
    levels = c("Epithelial/tumor", "Stroma", "Lymphocytes", "Other")
  )
  # Set consistent y-axis limits
  lims <- range(df$Score, na.rm = TRUE)
  
  # Plot combined violin
  p <- ggplot(df, aes(x = Annotation_group, y = Score, fill = Method)) +
    geom_violin(position = position_dodge(width = 0.8), trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.11, position = position_dodge(width = 0.8), outlier.shape = NA, color = "black", alpha = 0.6) +
    scale_fill_manual(values = c("GSVA" = "#377eb8", "sGSVA" = "#e41a1c")) +
    coord_cartesian(ylim = lims) +
    labs(
      title = pretty_pathway_short(pathway),
      x = NULL,
      y = "Pathway Score"
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  return(p)
}

# generate spatial maps
get_spatial_map <- function(spatial_obj, method = "sGSVA", pathway, lims) {
  DefaultAssay(spatial_obj) <- method
  p <- SpatialFeaturePlot(spatial_obj, features = pathway, pt.size.factor = 3) +
    ggtitle(paste0(method, ": ", pretty_pathway_short(pathway))) +
    scale_fill_gradientn(colors = c("#313695", "#ffffbf", "#a50026"),
                         limits = lims, oob = scales::squish) +
    annotation_scale(location = "bl", width_hint = 0.35, text_cex = 0.55) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.45, "cm")
    )
  return(p)
}
get_score_limits <- function(obj, assay, features) {
  DefaultAssay(obj) <- assay
  v <- GetAssayData(obj, layer = "data")[features, ]
  c(min(v, na.rm = TRUE), max(v, na.rm = TRUE))
}


# -----------------------------
# breast Visium main figures
# -----------------------------

# Example: save a few histograms for Moran's I
path_short <- function(pathway) {
  dplyr::case_when(
    grepl("B_CELL_RECEPTOR", pathway) ~ "BCR signaling",
    grepl("T_CELL_RECEPTOR", pathway) ~ "TCR signaling",
    grepl("GLUTATHIONE", pathway) ~ "Glutathione metabolism",
    grepl("FATTY_ACID", pathway) ~ "Fatty acid metabolism",
    TRUE ~ gsub("-", " ", pathway)
  )
}
hist_plot <- list()
path_moran <- c("KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY", 
                "KEGG_GLUTATHIONE_METABOLISM")
for (nm in path_moran) {
  df <- data.frame(val = hex_perm_null_results[[nm]]$moran_perm)
  p <- ggplot(df, aes(x = val)) +
    geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
    geom_vline(xintercept = hex_results_df$obs_moran[hex_results_df$pathway==nm], color = "firebrick", linetype = 2, linewidth = 1) +
    theme_bw() +
    labs(x = "Moran's I (permuted)", y = "Count",
         title = path_short(nm))
  hist_plot[[nm]] <- p
}

# violin plots
bcell_violin <- plot_violin(bc_obj, pathway = "KEGG-B-CELL-RECEPTOR-SIGNALING-PATHWAY")
gluta_violin <- plot_violin(bc_obj, pathway = "KEGG-GLUTATHIONE-METABOLISM")

# KEGG-B-CELL-RECEPTOR-SIGNALING-PATHWAY maps
lim_sgsva <- get_score_limits(bc_obj, "sGSVA", "KEGG-B-CELL-RECEPTOR-SIGNALING-PATHWAY")
lim_gsva  <- get_score_limits(bc_obj, "GSVA", "KEGG-B-CELL-RECEPTOR-SIGNALING-PATHWAY")
lims <- range(lim_sgsva, lim_gsva)
bcell_sgsva <- get_spatial_map(bc_obj, method = "sGSVA", pathway = "KEGG-B-CELL-RECEPTOR-SIGNALING-PATHWAY", lims)
bcell_gsva <- get_spatial_map(bc_obj, method = "GSVA", pathway = "KEGG-B-CELL-RECEPTOR-SIGNALING-PATHWAY", lims)
# KEGG-GLUTATHIONE-METABOLISM maps
lim_sgsva <- get_score_limits(bc_obj, "sGSVA", "KEGG-GLUTATHIONE-METABOLISM")
lim_gsva  <- get_score_limits(bc_obj, "GSVA", "KEGG-GLUTATHIONE-METABOLISM")
lims <- range(lim_sgsva, lim_gsva)
gluta_sgsva <- get_spatial_map(bc_obj, method = "sGSVA", pathway = "KEGG-GLUTATHIONE-METABOLISM", lims)
gluta_gsva <- get_spatial_map(bc_obj, method = "GSVA", pathway = "KEGG-GLUTATHIONE-METABOLISM", lims)

# scatterplots of BCR and glutathione
scatter_BCR <- ggplot(anchor_merge_dat, aes(x = B.cells, y = KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  ggpubr::stat_cor(method = "spearman", 
                   label.x = 0.05, label.y = max(anchor_merge_dat$KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY, na.rm=TRUE) + 0.5,
                   size = 5) +
  labs(x = "B-cell proportion", y = "sGSVA z-score (BCR Pathway)",
       title = "BCR Pathway vs B-cell") +
  theme_minimal()
scatter_Gluta <- ggplot(anchor_merge_dat, aes(x = Cancer.Epithelial, y = KEGG_GLUTATHIONE_METABOLISM)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  ggpubr::stat_cor(method = "spearman", 
                   label.x = 0.05, label.y = max(anchor_merge_dat$KEGG_GLUTATHIONE_METABOLISM, na.rm=TRUE) + 0.5,
                   size = 5) +
  labs(x = "Cancer Epithelial proportion", y = "sGSVA z-score (Glutathione metabolism)",
       title = "Glutathione metabolism vs Epithelial") +
  theme_minimal()

# combine figures
comb_fig <- (bcell_sgsva | bcell_gsva | bcell_violin) /
  ((gluta_sgsva+ theme(legend.position="none")) | (gluta_gsva+ theme(legend.position="none")) | gluta_violin) /
  (hist_plot$KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY | hist_plot$KEGG_GLUTATHIONE_METABOLISM |
     scatter_BCR | scatter_Gluta) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1, 0.9), guides = "collect")
ggsave(file.path(fig_dir, "Fig4_BreastVisium_main.pdf"),
       comb_fig, width = 16, height = 12)


# ------------------------------------
# breast Visium supplementary figures
# ------------------------------------
# supplementary figure S1 - Additional pathway spatial maps
# KEGG-T-CELL-RECEPTOR-SIGNALING-PATHWAY maps
lim_sgsva <- get_score_limits(bc_obj, "sGSVA", "KEGG-T-CELL-RECEPTOR-SIGNALING-PATHWAY")
lim_gsva  <- get_score_limits(bc_obj, "GSVA", "KEGG-T-CELL-RECEPTOR-SIGNALING-PATHWAY")
lims <- range(lim_sgsva, lim_gsva)
tcell_sgsva <- get_spatial_map(bc_obj, method = "sGSVA", pathway = "KEGG-T-CELL-RECEPTOR-SIGNALING-PATHWAY", lims)
tcell_gsva <- get_spatial_map(bc_obj, method = "GSVA", pathway = "KEGG-T-CELL-RECEPTOR-SIGNALING-PATHWAY", lims)
# KEGG-FATTY-ACID-METABOLISM maps
lim_sgsva <- get_score_limits(bc_obj, "sGSVA", "KEGG-FATTY-ACID-METABOLISM")
lim_gsva  <- get_score_limits(bc_obj, "GSVA", "KEGG-FATTY-ACID-METABOLISM")
lims <- range(lim_sgsva, lim_gsva)
fatty_sgsva <- get_spatial_map(bc_obj, method = "sGSVA", pathway = "KEGG-FATTY-ACID-METABOLISM", lims)
fatty_gsva <- get_spatial_map(bc_obj, method = "GSVA", pathway = "KEGG-FATTY-ACID-METABOLISM", lims)

# violin plots
tcell_violin <- plot_violin(bc_obj, pathway = "KEGG-T-CELL-RECEPTOR-SIGNALING-PATHWAY")
fatty_violin <- plot_violin(bc_obj, pathway = "KEGG-FATTY-ACID-METABOLISM")

# moran null hist
hex_results_df <- block_res[[2]]$permutation_results$hex_block_results$hex_block_summary
hex_perm_null_results <- block_res[[2]]$permutation_results$hex_block_results$hex_perm_null_results
hist_plot <- list()
top4_moran <- c("KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", 
                "KEGG_FATTY_ACID_METABOLISM")
for (nm in top4_moran) {
  df <- data.frame(val = hex_perm_null_results[[nm]]$moran_perm)
  p <- ggplot(df, aes(x = val)) +
    geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
    geom_vline(xintercept = hex_results_df$obs_moran[hex_results_df$pathway==nm], color = "firebrick", linetype = 2, linewidth = 1) +
    theme_bw() +
    labs(x = "Moran's I (permuted)", y = "Count",
         title = path_short(nm))
  hist_plot[[nm]] <- p
}

# scatterplots of TCR and Fatty acid
scatter_TCR <- ggplot(anchor_merge_dat, aes(x = T.cells, y = KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  ggpubr::stat_cor(method = "spearman", 
                   label.x = 0.05, label.y = max(anchor_merge_dat$KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY, na.rm=TRUE) + 0.5,
                   size = 5) +
  labs(x = "T-cell proportion", y = "sGSVA z-score (TCR Pathway)",
       title = "TCR Pathway vs T-cell") +
  theme_minimal()
scatter_Fatty <- ggplot(anchor_merge_dat, aes(x = Cancer.Epithelial, y = KEGG_FATTY_ACID_METABOLISM)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  ggpubr::stat_cor(method = "spearman", 
                   label.x = 0.05, label.y = max(anchor_merge_dat$KEGG_FATTY_ACID_METABOLISM, na.rm=TRUE) + 0.5,
                   size = 5) +
  labs(x = "Cancer Epithelial proportion", y = "sGSVA z-score (Fatty acid metabolism)",
       title = "Fatty acid metabolism vs Epithelial") +
  theme_minimal()

comb_map <- (tcell_sgsva | tcell_gsva | tcell_violin) /
  ((fatty_sgsva+theme(legend.position="none")) | (fatty_gsva+theme(legend.position="none")) | fatty_violin) /
  (hist_plot$KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY | hist_plot$KEGG_FATTY_ACID_METABOLISM |
     scatter_TCR | scatter_Fatty) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1, 0.9), guides = "collect")
ggsave(file.path(fig_dir, "Supplement_Fig_S6_BreastVisium_additional_pathways_results.pdf"),
       comb_map, width = 16, height = 12)

# supplementary figure S2 — Full spatial permutation diagnostics
pathway <- "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY"
methods <- c("hex-block 150", "hex-block 200", "hex-block 300", "Global permutation")
plot_null_hist <- function(stat = c("moran", "geary"), method_index, method_label) {
  
  stat <- match.arg(stat)
  
  if (grepl("hex-block", method_label)) {
    
    hex_results_df <- block_res[[method_index]]$permutation_results$hex_block_results$hex_block_summary
    hex_perm_null_results <- block_res[[method_index]]$permutation_results$hex_block_results$hex_perm_null_results
    if (stat == "moran") {
      null_vals <- hex_perm_null_results[[pathway]]$moran_perm
      obs_val <- hex_results_df$obs_moran[hex_results_df$pathway == pathway]
      x_lab <- "Moran's I (permuted)"
    } else {
      null_vals <- hex_perm_null_results[[pathway]]$geary_perm
      obs_val <- hex_results_df$obs_geary[hex_results_df$pathway == pathway]
      x_lab <- "Geary's C (permuted)"
    }
    
  } else {
    mc_results_df <- sgsva_res$permutation_results$moran_geary_results$moran_geary_summary
    mc_perm_null_results <- sgsva_res$permutation_results$moran_geary_results$mc_perm_null_results
    
    if (stat == "moran") {
      null_vals <- mc_perm_null_results[[pathway]]$moran_perm
      obs_val <- mc_results_df$obs_moran[mc_results_df$pathway == pathway]
      x_lab <- "Moran's I (permuted)"
    } else {
      null_vals <- mc_perm_null_results[[pathway]]$geary_perm
      obs_val <- mc_results_df$obs_geary[mc_results_df$pathway == pathway]
      x_lab <- "Geary's C (permuted)"
    }
  }
  
  df <- data.frame(val = null_vals)
  ggplot(df, aes(x = val)) +
    geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
    geom_vline(xintercept = obs_val, color = "firebrick", linetype = 2, linewidth = 0.8) +
    theme_bw(base_size = 8) +
    labs(
      x = x_lab,
      y = "Count",
      title = method_label
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      panel.grid.minor = element_blank()
    )
}
moran_list <- list(
  plot_null_hist("moran", 1, "Hex-block 150 µm"),
  plot_null_hist("moran", 2, "Hex-block 200 µm"),
  plot_null_hist("moran", 3, "Hex-block 300 µm"),
  plot_null_hist("moran", 4, "Global permutation")
)

geary_list <- list(
  plot_null_hist("geary", 1, "Hex-block 150 µm"),
  plot_null_hist("geary", 2, "Hex-block 200 µm"),
  plot_null_hist("geary", 3, "Hex-block 300 µm"),
  plot_null_hist("geary", 4, "Global permutation")
)

comb_hist <- 
  (moran_list[[1]] | moran_list[[2]] | moran_list[[3]] | moran_list[[4]]) /
  (geary_list[[1]] | geary_list[[2]] | geary_list[[3]] | geary_list[[4]]) +
  plot_annotation(
    title = "KEGG B cell receptor signaling pathway",
    subtitle = "Permutation null distributions for Moran's I and Geary's C",
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )   
ggsave(file.path(fig_dir, "Supplement_Fig_S7_BreastVisium_perm_null_hist.pdf"),
       comb_hist, width = 12, height = 6, units = "in")







