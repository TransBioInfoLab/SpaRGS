# ============================================================
# Breast Visium figures
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggspatial)
  library(scales)
  library(here)
})

# ------------------------------------------------------------
# Paths and reproducibility
# ------------------------------------------------------------

project_dir <- "."
results_dir   <- file.path(project_dir, "results", "breast_visium")
anchoring_dir <- file.path(results_dir, "celltype_anchoring")
pilot_dir     <- file.path(results_dir, "enrichmap_pilot")
fig_dir       <- file.path(project_dir, "results", "figures")

source(file.path(project_dir,"R/spargs_score.R"))

# ---------------------------
# load saved outputs
# ---------------------------
bc_obj <- readRDS(file.path(results_dir, "bc_obj.rds"))
prespecified_results <- read.csv(file.path(anchoring_dir, "anchoring_prespecified_pairs.csv"))

# main comparator scores
spa_scores <- read.csv(file.path(results_dir, "spargs_scores.csv"), row.names = 1)
gsva_scores <- read.csv(file.path(results_dir, "gsva_scores.csv"), row.names = 1)
gsva_smoothed_scores <- read.csv(file.path(results_dir, "gsva_smoothed_scores.csv"), row.names = 1)
gsdensity_scores <- read.csv(file.path(results_dir, "gsdensity_scores.csv"), row.names = 1)

# EnrichMap output
em_default <- read.csv(
  file.path(pilot_dir, "breast_visium_enrichmap_default_corrected_all_scores.csv"),
  row.names = 1,
  check.names = FALSE
) 
em_smooth <- read.csv(
  file.path(pilot_dir, "breast_visium_enrichmap_smooth_only_all_scores.csv"),
  row.names = 1,
  check.names = FALSE
)


# ------------------------------------------------
# Name harmonization
# ------------------------------------------------
# standardize barcodes
standardize_barcodes <- function(x) {
  gsub(".", "-", x, fixed = TRUE)
}

# normalize pathway
normalize_pathway_names <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  sub("^KEGG_", "", x)
}

method_labels <- c(
  SpaRGS = "SpaRGS",
  GSVA = "GSVA",
  GSVA_smoothed = "GSVA-smoothed",
  GSDensity = "GSDensity",
  EnrichMap_default = "EnrichMap-default",
  EnrichMap_smooth = "EnrichMap-smooth"
)

pathway_labels <- c(
  B_CELL_RECEPTOR_SIGNALING_PATHWAY = "B-cell receptor signaling",
  FATTY_ACID_METABOLISM = "Fatty-acid metabolism"
)


# ------------------------------------------------------------
# Prepare score matrices
# ------------------------------------------------------------

score_list_all <- list(
  SpaRGS = as.matrix(spa_scores),
  GSVA = as.matrix(gsva_scores),
  GSVA_smoothed = as.matrix(gsva_smoothed_scores),
  GSDensity = as.matrix(gsdensity_scores),
  EnrichMap_default = as.matrix(em_default),
  EnrichMap_smooth = as.matrix(em_smooth)
)

score_list_all <- lapply(score_list_all, function(mat) {
  mat <- as.matrix(mat)
  
  rownames(mat) <- normalize_pathway_names(rownames(mat))
  colnames(mat) <- standardize_barcodes(colnames(mat))
  
  if (anyDuplicated(rownames(mat))) {
    stop("Duplicated normalized pathway names were detected.")
  }
  
  mat
})

# Confirm that every method contains all Visium spots
spot_check <- vapply(
  score_list_all,
  function(mat) {
    all(colnames(bc_obj) %in% colnames(mat))
  },
  logical(1)
)

if (!all(spot_check)) {
  stop(
    "Some score matrices are missing Visium spots: ",
    paste(names(spot_check)[!spot_check], collapse = ", ")
  )
}

# Z-standardize without changing score polarity
z_scores_all <- lapply(score_list_all, z_normalization)

# Add each matrix once as a display assay
for (method in names(z_scores_all)) {
  mat <- z_scores_all[[method]]
  
  mat <- mat[, colnames(bc_obj), drop = FALSE]
  rownames(mat) <- gsub("_", "-", rownames(mat))
  
  assay_name <- paste0(method, "_display")
  
  bc_obj[[assay_name]] <- CreateAssayObject(data = mat)
}


# ------------------------------------------------------------
# Generate one spatial map
# ------------------------------------------------------------
get_spatial_map_main <- function(
    spatial_obj,
    assay_name,
    pathway,
    limits
) {
  DefaultAssay(spatial_obj) <- assay_name
  
  method_key <- sub("_display$", "", assay_name)
  method_title <- unname(method_labels[method_key])
  
  SpatialFeaturePlot(
    spatial_obj,
    features = pathway,
    pt.size.factor = 3
  ) +
    ggtitle(method_title) +
    scale_fill_gradientn(
      colors = c("#313695", "#ffffbf", "#a50026"),
      limits = limits,
      oob = scales::squish
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.height = grid::unit(1.4, "cm"),
      legend.key.width = grid::unit(0.28, "cm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# ------------------------------------------------------------
# Make one horizontal four-method pathway row
# ------------------------------------------------------------
make_four_method_row <- function(
    spatial_obj,
    assay_names,
    pathway,
    pathway_title
) {
  values_all <- unlist(
    lapply(assay_names, function(assay_name) {
      DefaultAssay(spatial_obj) <- assay_name
      FetchData(spatial_obj, vars = pathway)[, 1]
    })
  )
  
  # Same color limits across methods within this pathway
  limit_max <- unname(
    quantile(abs(values_all), probs = 0.99, na.rm = TRUE)
  )
  
  shared_limits <- c(-limit_max, limit_max)
  
  map_list <- lapply(
    assay_names,
    function(assay_name) {
      get_spatial_map_main(
        spatial_obj = spatial_obj,
        assay_name = assay_name,
        pathway = pathway,
        limits = shared_limits
      )
    }
  )
  
  patchwork::wrap_plots(
    map_list, nrow = 1, guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = pathway_title,
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 13,
          margin = margin(b = 4)
        )
      )
    ) &
    theme(legend.position = "right")
}

# ------------------------------------------------------
# breast Visium main figures
# ------------------------------------------------------
# representative spatial maps
main_assays <- c(
  "SpaRGS_display",
  "GSVA_display",
  "GSVA_smoothed_display",
  "GSDensity_display"
)

p_bcr_main <- make_four_method_row(
  spatial_obj = bc_obj,
  assay_names = main_assays,
  pathway = "B-CELL-RECEPTOR-SIGNALING-PATHWAY",
  pathway_title = "B-cell receptor signaling"
) 
p_fatty_main <- make_four_method_row(
  spatial_obj = bc_obj,
  assay_names = main_assays,
  pathway = "FATTY-ACID-METABOLISM",
  pathway_title = "Fatty-acid metabolism"
) 

# Pathway-cell-type spatial concordance heatmap
method_order <- c(
  "SpaRGS",
  "GSVA",
  "GSVA-smoothed",
  "GSDensity"
)

anchor_label_short <- c(
  "BCR signaling vs B cells" =
    "BCR - B cells",
  "TCR signaling vs T cells" =
    "TCR - T cells",
  "Glutathione metabolism vs cancer epithelial" =
    "Glutathione - cancer epithelial",
  "Fatty-acid metabolism vs cancer epithelial" =
    "Fatty acid - cancer epithelial"
)

anchor_order <- c(
  "BCR - B cells",
  "TCR - T cells",
  "Glutathione - cancer epithelial",
  "Fatty acid - cancer epithelial"
)

anchoring_plot_df <- prespecified_results %>%
  mutate(
    method = recode(
      method,
      "GSVA_smoothed" = "GSVA-smoothed",
      .default = method
    ),
    method = factor(method, levels = method_order),
    
    anchor_short = recode(anchor_label, !!!anchor_label_short),
    anchor_short = factor(anchor_short, levels = rev(anchor_order)),
    
    signed_rho = as.numeric(rho),
    abs_rho = abs(signed_rho),
    signed_label = sprintf("%+.2f", signed_rho)
  )

p_anchor <- ggplot(anchoring_plot_df,
  aes(x = method, y = anchor_short, fill = abs_rho)) +
  geom_tile(color = "white", linewidth = 0.8,
    width = 0.98, height = 0.98) +
  geom_text(aes(label = signed_label), size = 4) +
  scale_fill_gradient(
    limits = c(0, 1),
    low = "white",
    high = "steelblue",
    oob = scales::squish,
    name = expression("|Spearman " * rho * "|")
  ) +
  scale_x_discrete(
    labels = c("SpaRGS", "GSVA", "GSVA-smoothed", "GSDensity")
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Magnitude of pathway-cell-type spatial concordance"
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(5, "cm"),
      barheight = grid::unit(0.35, "cm")
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 9
    ),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0,
      size = 13
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 9
    ),
    
    legend.position = "bottom",
    legend.margin = margin(t = -3),
    plot.margin = margin(4, 4, 4, 4)
  )

# combine three figures
heatmap_row <- (
  patchwork::plot_spacer() |
    p_anchor |
    patchwork::plot_spacer()
) +
  patchwork::plot_layout(
    widths = c(0.75, 2.5, 0.75)
  )

panel_A <- patchwork::wrap_elements(
  full = patchwork::patchworkGrob(heatmap_row)
)
panel_B <- patchwork::wrap_elements(
  full = patchwork::patchworkGrob(p_bcr_main)
)
panel_C <- patchwork::wrap_elements(
  full = patchwork::patchworkGrob(p_fatty_main)
)

fig_breast_main <- (panel_A / panel_B / panel_C) +
  patchwork::plot_layout(
    heights = c(0.95, 1.15, 1.15)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    plot.tag.position = c(0.01, 0.98)
  )
ggsave(
  filename = file.path(fig_dir, "Fig4_BreastVisium_main.pdf"),
  plot = fig_breast_main,
  width = 10,
  height = 10
)



# --------------------------------------------------------
# breast Visium supplementart figures
# --------------------------------------------------------
# ------------------------------------------------------------
# Make one multi-method pathway panel
# ------------------------------------------------------------
make_pathway_panel <- function(
    spatial_obj,
    assay_names,
    pathway,
    ncol = 3,
    pathway_title = NULL,
    quantile_limit = 0.99
) {
  # Check assay availability
  missing_assays <- setdiff(assay_names, names(spatial_obj@assays))
  if (length(missing_assays) > 0) {
    stop(
      "These assays are missing from the Seurat object: ",
      paste(missing_assays, collapse = ", ")
    )
  }
  
  # Check pathway availability in each assay
  pathway_check <- vapply(
    assay_names,
    function(assay_name) {
      pathway %in% rownames(spatial_obj[[assay_name]])
    },
    logical(1)
  )
  
  if (!all(pathway_check)) {
    stop(
      "Pathway not found in these assays: ",
      paste(assay_names[!pathway_check], collapse = ", "),
      "\nPathway requested: ", pathway
    )
  }
  
  # Extract all values across selected assays to define common color limits
  values_all <- unlist(
    lapply(assay_names, function(assay_name) {
      DefaultAssay(spatial_obj) <- assay_name
      FetchData(spatial_obj, vars = pathway)[, 1]
    })
  )
  
  limit_max <- unname(
    quantile(abs(values_all), probs = quantile_limit, na.rm = TRUE)
  )
  
  shared_limits <- c(-limit_max, limit_max)
  
  # Generate one map per method
  map_list <- lapply(
    assay_names,
    function(assay_name) {
      get_spatial_map_main(
        spatial_obj = spatial_obj,
        assay_name = assay_name,
        pathway = pathway,
        limits = shared_limits
      )
    }
  )
  
  # Default title
  if (is.null(pathway_title)) {
    pathway_title <- pathway
  }
  
  # Combine maps
  patchwork::wrap_plots(
    map_list,
    ncol = ncol,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = pathway_title,
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 15,
          margin = margin(b = 6)
        )
      )
    ) &
    theme(
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA)
    )
}
supp_assays <- c(
  "SpaRGS_display",
  "GSVA_display",
  "GSVA_smoothed_display",
  "GSDensity_display",
  "EnrichMap_default_display",
  "EnrichMap_smooth_display"
)

p_bcr_supp <- make_pathway_panel(
  spatial_obj = bc_obj,
  assay_names = supp_assays,
  pathway = "B-CELL-RECEPTOR-SIGNALING-PATHWAY",
  ncol = 3,
  pathway_title = "B-cell receptor signaling"
)

p_fatty_supp <- make_pathway_panel(
  spatial_obj = bc_obj,
  assay_names = supp_assays,
  pathway = "FATTY-ACID-METABOLISM",
  ncol = 3,
  pathway_title = "Fatty-acid metabolism"
)

ggsave(
  file.path(fig_dir, "supplement_Fig_S3_breast_BCR_six_methods.pdf"),
  p_bcr_supp,
  width = 11,
  height = 7.5
)

ggsave(
  file.path(fig_dir, "supplement_Fig_S4_breast_fatty_acid_six_methods.pdf"),
  p_fatty_supp,
  width = 11,
  height = 7.5
)






