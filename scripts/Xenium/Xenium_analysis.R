### Real data case3 for Xenium FFPE
suppressPackageStartupMessages ({
  library(Seurat); library(Matrix);library(org.Hs.eg.db)
  library(gridExtra); library(jsonlite); library(ggplot2)
  library(viridis); library(spdep); library(msigdbr)
  library(dplyr); library(purrr) ; library(FNN)
  library(GSVA); library(hexbin); library(readr); 
  library(tidyr); library(stringr)
})

set.seed(1)

# ---------------------------
# paths and reproducibility
# ---------------------------
project_dir <- "."
data_dir <- file.path(project_dir, "data", "xenium")
results_dir <- file.path(project_dir, "results", "xenium")
fig_dir <- file.path(project_dir, "results", "figures")
table_dir <- file.path(project_dir, "results", "tables")
supp_dir <- file.path(project_dir, "results", "supplementary_data")

source(file.path(project_dir, "R/sgsva_score.R"))
source(file.path(project_dir, "R/spatial_permutation.R"))


# ---------------------------
# load data
# ---------------------------
# load data
xe <- ReadXenium(
  data.dir = data_dir,
  outs     = "matrix",        
  type     = "centroids",     
  mols.qv.threshold = 20      
)
counts <- xe$matrix$`Gene Expression`
centroids <- xe$centroids

# extract KEGG pathway sets
kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
gene_sets <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)


# ---------------------------
# Preprocess data
# ---------------------------
# panel-aware gene-set filtering
panel_filter_sets <- function(counts,
                              sets,
                              min_genes = 15,
                              min_coverage = 0.30) {
  # Panel gene universe from Xenium
  panel_genes <- rownames(counts)
  panel_key <- panel_genes
  names(panel_key) <- panel_genes
  sets_upper <- lapply(sets, unique)
  
  # For each pathway: intersect with panel and compute coverage
  stats_list <- vector("list", length(sets))
  names(stats_list) <- names(sets)
  
  sets_on_panel <- vector("list", length(sets))
  names(sets_on_panel) <- names(sets)
  
  for (nm in names(sets)) {
    gs_all   <- sets_upper[[nm]]  # full gene set
    n_total  <- length(gs_all)
    
    in_panel <- gs_all %in% panel_key # which of these are on the panel
    n_panel  <- sum(in_panel)
    
    # Map back to original panel gene symbols (respect panel’s case)
    if (n_panel > 0) {
      matched_panel <- panel_genes[match(gs_all[in_panel], panel_key)]
      matched_panel <- unique(matched_panel)
    } else {
      matched_panel <- character(0)
    }
    
    sets_on_panel[[nm]] <- matched_panel
    
    stats_list[[nm]] <- tibble(
      pathway        = nm,
      n_total_genes  = n_total,
      n_panel_genes  = n_panel,
      coverage       = ifelse(n_total > 0, n_panel / n_total, NA_real_)
    )
  }
  
  stats <- bind_rows(stats_list) %>%
    mutate(
      keep = (n_panel_genes >= min_genes) | (coverage >= min_coverage)
    )
  
  sets_filtered <- sets_on_panel[stats$pathway[stats$keep]]
  
  list(
    sets_filtered = sets_filtered,  # only panel genes, only kept pathways
    stats         = stats,          # table for report (coverage%)
    sets_on_panel = sets_on_panel   # all pathways, restricted to panel genes
  )
}

pf <- panel_filter_sets(
  counts       = counts,
  sets         = gene_sets,
  min_genes    = 15,
  min_coverage = 0.30
)
kegg_sets_xenium <- pf$sets_filtered
genes_used <- unique(unlist(kegg_sets_xenium))
counts_sub <- counts[intersect(rownames(counts), genes_used), , drop = FALSE]

# coverage QC table/plots
stats <- pf$stats %>%
  arrange(desc(keep), desc(n_panel_genes))
write.csv(stats, file.path(results_dir, "xenium_kegg_panel_coverage.csv"))

# histogram of QC
coverage_table <- stats %>%
  mutate(
    coverage_pct = coverage*100,
    status = ifelse(keep==TRUE, "Retained", "Filtered")
  )
status_cols <- c("Filtered" = "#E8A19A", "Retained" = "#7FC4C9")
gg_coverage_hist <- ggplot(coverage_table,
                           aes(x = coverage_pct, fill = status)) +
  geom_histogram(binwidth = 5, color = "black", alpha = 0.7) +
  geom_vline(xintercept = 30, linetype = "dashed") +
  scale_fill_manual(values = status_cols, name = "Pathway status") +
  labs(
    x = "Panel coverage per pathway (%)",
    y = "Number of pathways",
    title = "Xenium panel coverage of gene sets"
  ) +
  theme_classic()

# Cell-level QC filtering (based on the quantiles)
cell_total_counts <- Matrix::colSums(counts_sub)
cell_n_features   <- Matrix::colSums(counts_sub > 0)
keep_cells <- (cell_total_counts >= 30) & (cell_n_features >= 15)

# Scatter: nCount vs nFeature, colored by keep/drop
qc_df <- data.frame(
  cell_id    = colnames(counts_sub),
  nCount     = as.numeric(cell_total_counts),
  nFeature   = as.numeric(cell_n_features),
  keep_cell  = ifelse(
    (cell_total_counts >= 30) & (cell_n_features >= 15),"Retained","Filtered")
)
cell_qc_plot <- ggplot(qc_df, aes(x = nCount, y = nFeature, color = keep_cell)) +
  geom_point(alpha = 0.3, size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  geom_vline(xintercept = 30, linetype = "dashed") +
  geom_hline(yintercept = 15, linetype = "dashed") +
  scale_color_manual(values = status_cols, name = "Cell status") +
  labs(x = "Total transcripts per cell (log10)",
       y = "Detected genes per cell (log10)",
       title = "Xenium FFPE cell-level QC") +
  theme_minimal()

# keep cells
counts_qc <- counts_sub[, keep_cells, drop = FALSE]
centroids <- centroids[keep_cells, , drop = FALSE]
cell_total_counts <- cell_total_counts[keep_cells]
cell_n_features   <- cell_n_features[keep_cells]

n_cells_before <- ncol(counts_sub)
n_cells_after  <- ncol(counts_qc)    # 272,469
frac_kept      <- n_cells_after / n_cells_before # 39%
summary_df <- data.frame(
  n_cells_before = n_cells_before,
  n_cells_after  = n_cells_after,
  frac_kept      = frac_kept
)

# Gene-level filtering
gene_n_cells <- Matrix::rowSums(counts_qc > 0)
keep_genes   <- gene_n_cells >= 20
counts_qc <- counts_qc[keep_genes, , drop = FALSE]

# Histogram of log10 prevalence
gene_qc_df <- data.frame(
  gene        = rownames(counts_sub),
  n_cells     = as.numeric(gene_n_cells),
  keep_gene20 = ifelse(gene_n_cells >= 20, "Retained", "Filtered")
)
gene_qc_plot <- ggplot(gene_qc_df, aes(x = n_cells, fill = keep_gene20)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_x_log10() +
  scale_fill_manual(values = status_cols, name = "Gene status") +
  labs(x = "Number of cells with nonzero counts (log10)",
       y = "Genes",
       title = "Gene-level prevalence and QC threshold") +
  theme_minimal()

genes_panel    <- rownames(counts_sub)
genes_after_qc <- rownames(counts_qc)
# genes per pathway before vs after QC
path_qc <- lapply(kegg_sets_xenium, function(gs) {
  before <- length(intersect(gs, genes_panel))
  after  <- length(intersect(gs, genes_after_qc))
  c(before = before, after = after)
})
path_qc_df <- do.call(rbind, path_qc)
head(path_qc_df)

coords <- as.data.frame(centroids[, c("x", "y")])
rownames(coords) <- centroids$cell

# metadata for downstream QC / plotting
meta <- data.frame(
  cell_id          = colnames(counts_qc),
  nCount           = as.numeric(cell_total_counts),
  nFeature         = as.numeric(cell_n_features),
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$cell_id

# Update gene sets to final universe
final_genes <- rownames(counts_qc)
kegg_sets_final <- lapply(kegg_sets_xenium, function(gs) intersect(gs, final_genes))

final_coverage <- tibble::tibble(
  pathway = names(kegg_sets_final),
  n_panel_genes_before_gene_qc = lengths(kegg_sets_xenium[names(kegg_sets_final)]),
  n_final_genes = lengths(kegg_sets_final)
) %>%
  mutate(
    keep_final = n_final_genes >= 0
  )
kegg_sets_final <- kegg_sets_final[final_coverage$keep_final]

write.csv(
  final_coverage,
  file.path(results_dir, "xenium_kegg_final_gene_coverage.csv"),
  row.names = FALSE
)

# save results for downstream analysis
xenium_process_data <- list(
  counts_qc = counts_qc, coords = coords,
  kegg_sets_final = kegg_sets_final
)
saveRDS(xenium_process_data, file.path(results_dir, "xenium_process_data.rds"))

# supplementary S3 for Xenium preprocessing
comb_qc <- (gg_coverage_hist) /
  (cell_qc_plot|gene_qc_plot) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")
ggsave(filename=file.path(fig_dir, "supplement_Fig_S11_Xenium_preprocessing.pdf"),
       plot=comb_qc, width = 10, height = 10)



# ---------------------------------
# computing three pathway scores
# ---------------------------------
# generate runtime table
runtime_results <- tibble::tibble(
  method = character(),
  runtime = numeric()
)

# sGSVA
optimal_sigma <- estimate_sigma_knn(coords, k=8) # 24
lib <- colSums(counts_qc)+1
lcp <- log1p(t(t(counts_qc)/lib*1e6)) #logCPM
expr <- as.matrix(lcp)

t0 <- proc.time()[3]
sgsva_res <- run_sgsva(expr, coords, kegg_sets_final,
                       sigma_global = 24, cor_min = 0.1, 
                       mode = "knn", k=8, one_sided = "center_0.5")
rt <- proc.time()[3] - t0
runtime_results <- bind_rows(runtime_results,
                             tibble(method = "sGSVA", runtime = rt))

saveRDS(sgsva_res, file.path(results_dir, "sgsva_score_res.rds"))

sgsva_res <- readRDS(file.path(results_dir, "sgsva_score_res.rds"))
z_sgsva <- sgsva_res$z_score
balanced_sets <- sgsva_res$gene_sets
orient <- sgsva_res$directionality_diagnostics
write.csv(orient, 
          file.path(results_dir, "xenium_sgsva_orientation_diagnostics.csv"), row.names = FALSE)

# GSVA
lib <- colSums(counts_qc)+1
lcp <- log1p(t(t(counts_qc)/lib*1e6)) #logCPM
expr <- as.matrix(lcp)
keep <- rowSums(expr > 0) >= floor(0 * ncol(expr))
expr <- expr[keep, , drop = FALSE]

t0 <- proc.time()[3]
params <- gsvaParam(
  expr, kegg_sets_final
)
gsva_scores <- gsva(params, verbose = FALSE)
rt <- proc.time()[3] - t0
runtime_results <- bind_rows(runtime_results,
                             tibble(method = "GSVA", runtime = rt))

write.csv(gsva_scores, file.path(results_dir, "gsva_score.csv"))

gsva_scores <- read.csv(file.path(results_dir, "gsva_score.csv"))
rownames(gsva_scores) <- gsva_scores$X
gsva_scores$X <- NULL
z_gsva <- z_normalization(as.matrix(gsva_scores))
colnames(z_gsva) <- gsub("\\.","-",colnames(z_gsva))

# GSVA on smoothed
t0 <- proc.time()[3]
smoothed <- smooth_expr(expr, coords, sigma=24, mode="knn", k=8)
smoothed <- as.matrix(smoothed)
colnames(smoothed) <- colnames(expr)
params <- gsvaParam(
  smoothed, kegg_sets_final
)
gsva_smooth <- gsva(params, verbose = FALSE)
rt <- proc.time()[3] - t0
runtime_results <- bind_rows(runtime_results,
                             tibble(method = "GSVA-smoothed", runtime = rt))
write.csv(runtime_results, file.path(results_dir, "Xenium_runtime_summary.csv"))
write.csv(gsva_smooth, file.path(results_dir, "gsva_smooth_score.csv"))

gsva_smooth <- read.csv(file.path(results_dir, "gsva_smooth_score.csv"))
rownames(gsva_smooth) <- gsva_smooth$X
gsva_smooth$X <- NULL
z_gsva_smooth <- z_normalization(as.matrix(gsva_smooth))
colnames(z_gsva_smooth) <- gsub("\\.","-",colnames(z_gsva_smooth))


# ---------------------------------
# summary Xenium coverage
# ---------------------------------
# Collect score matrices
score_list <- list(
  sGSVA = as.matrix(sgsva_res$scores),
  GSVA = as.matrix(gsva_scores),
  GSVA_smoothed = as.matrix(gsva_smooth)
)
# Restrict all methods to the same spots and shared pathways
common_spots <- Reduce(intersect,
                       c(lapply(score_list, colnames),list(rownames(coords))))
shared_pathways <- Reduce(intersect, lapply(score_list, rownames))

# supplementary xenium coverage summary
xenium_coverage_sum <- tibble::tibble(
  quantity = c(
    "Cells before QC",
    "Cells after QC",
    "Genes in Xenium count matrix",
    "Genes represented in retained panel-filtered pathways",
    "Genes retained after final cell/gene QC",
    "KEGG pathways queried",
    "KEGG pathways passing initial panel filter",
    "KEGG pathways passing final QC gene filter",
    "Pathways retained by data-driven sGSVA",
    "Pathways scored by GSVA",
    "Pathways scored by GSVA-smoothed",
    "Shared pathways across three methods"
  ),
  value = c(
    ncol(counts),
    ncol(counts_qc),
    nrow(counts),
    nrow(counts_sub),
    nrow(counts_qc),
    length(gene_sets),
    length(kegg_sets_xenium),
    length(kegg_sets_final),
    nrow(sgsva_res$scores),
    nrow(gsva_scores),
    nrow(gsva_smooth),
    length(shared_pathways)
  )
)
write.csv(xenium_coverage_sum, file.path(table_dir, "supplement_Table_S11_xenium_coverage_summary.csv"))



# ============================================================
# Xenium pathway-module concordance analyses
# ============================================================
# ---------------------------
# marker score panel
# ---------------------------
gene_groups <- read.csv(file.path(data_dir, "Xenium_V1_FFPE_gene_groups.csv"))
# keep only genes on the Xenium panel
panel_genes <- rownames(counts_qc)  # genes × cells Xenium matrix after QC
marker_list <- gene_groups %>%
  filter(group!="") %>%
  filter(gene %in% panel_genes) %>%
  group_by(group) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")
marker_sets <- setNames(marker_list$genes, marker_list$group)

# creat a seurat object for counts
xen <- CreateSeuratObject(counts = counts_qc, assay = "Xenium")
# add coordinates to meta.data
xen@meta.data$x <- coords[colnames(counts_qc), "x"]
xen@meta.data$y <- coords[colnames(counts_qc), "y"]
xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 1e4)

# score each cell type
xen <- AddModuleScore(
  object   = xen,
  features = marker_sets,
  name     = "GG_",
  ctrl     = 25,   # fewer control genes per bin
  nbin     = 24
)
gg_cols <- grep("^GG_", colnames(xen@meta.data), value = TRUE)
colnames(xen@meta.data)[match(gg_cols, colnames(xen@meta.data))] <-
  paste0("GG_", names(marker_sets))
gg_cols <- paste0("GG_", names(marker_sets))

# assign each cell with best gene group
meta <- xen@meta.data
score_mat <- as.matrix(meta[, gg_cols])
rownames(score_mat) <- rownames(meta)
max_idx   <- max.col(score_mat, ties.method = "first")
best_group <- names(marker_sets)[max_idx]

# choose max weight
max_score <- apply(score_mat, 1, max)
score_threshold <- 0
best_group[max_score < score_threshold] <- "Unassigned"
meta$gene_group_label <- best_group
xen@meta.data <- meta
table(meta$gene_group_label)


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

# bin-level cell type fractions
meta <- xen@meta.data[common_cells, , drop = FALSE]
meta$bin_id <- coords_sub$bin_id

# myeloid: macrophages + monocytes + dendritic
meta$GG_Myeloid <- rowMeans(
  meta[, c("GG_Macrophages", "GG_Monocytes", "GG_Dendritic cells")],
  na.rm = TRUE
)
# total immune: T, B, NK, dendritic, macrophages, monocytes, mast
immune_cols <- c("GG_T cells", "GG_B cells", "GG_NK cells",
                 "GG_Dendritic cells", "GG_Macrophages",
                 "GG_Monocytes", "GG_Mast cells")
meta$GG_total_immune <- rowMeans(meta[, immune_cols], na.rm = TRUE)
ct_cols <- grep("^GG_", colnames(meta), value = TRUE)
write.csv(meta, file.path(results_dir, "xenium_cell_type_meta_data.csv"))


# ---------------------------
# bin-level correlation
# ---------------------------
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
bin_scores_chemokine <- aggregate_bin_level("KEGG_CHEMOKINE_SIGNALING_PATHWAY")

# bin-level correlations
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
# check chemokine pathway with multiple local fractions
anchors <- c("GG_T cells_mean",
             "GG_Myeloid_mean",
             "GG_Fibroblasts_mean",
             "GG_Endothelial cells_mean",
             "GG_total_immune_mean")

chem_anchoring <- do.call(rbind, lapply(anchors, function(ct_col) {
  stats <- bin_anchoring_stats(bin_scores_chemokine, ct_col)
  stats$anchor <- ct_col
  stats
}))
write.csv(chem_anchoring, file.path(results_dir, "chemokine_anchoring_summary.csv"))

# correlation summary report for all pathways and cell types
pathways <- rownames(z_sgsva_sub)
all_stats <- list()
for (pw in pathways) {
  bin_scores <- aggregate_bin_level(pw)
  for (ct_col in ct_cols) {
    cell_type_name <- paste0(ct_col, "_mean")
    stats <- bin_anchoring_stats(bin_scores, cell_type_name)
    stats$pathway  <- pw
    stats$cell_type <- ct_col
    
    all_stats[[length(all_stats) + 1]] <- stats
  }
}
anchoring_summary <- do.call(rbind, all_stats)
write.csv(anchoring_summary,
          file = file.path(results_dir, "xenium_bin_level_anchoring_summary.csv"),
          row.names = FALSE)



# ============================================================
# Xenium pathway-module gene overlap and non-overlap sensitivity
# ============================================================
# This checks whether the main pathway-module concordance results
# are driven by shared genes between the KEGG pathway and marker module.
# It is focused on the four representative pathway-module pairs used in Fig. 6.

# ---------------------------
# Define representative pairs
# ---------------------------
xenium_pairs <- tibble::tribble(
  ~anchor_label,                  ~pathway,                            ~module,
  "Cell cycle - breast cancer",   "KEGG_CELL_CYCLE",                   "Breast cancer",
  "Chemokine - T cells",          "KEGG_CHEMOKINE_SIGNALING_PATHWAY",  "T cells",
  "VEGF - endothelial",           "KEGG_VEGF_SIGNALING_PATHWAY",       "Endothelial cells",
  "Focal adhesion - fibroblasts", "KEGG_FOCAL_ADHESION",               "Fibroblasts"
)

# ---------------------------
# Gene-overlap table
# ---------------------------
overlap_list <- list()

for (i in seq_len(nrow(xenium_pairs))) {
  pw <- xenium_pairs$pathway[i]
  mod <- xenium_pairs$module[i]
  
  pathway_genes <- intersect(kegg_sets_final[[pw]], rownames(counts_qc))
  module_genes  <- intersect(marker_sets[[mod]], rownames(counts_qc))
  overlap_genes <- intersect(pathway_genes, module_genes)
  nonoverlap_module_genes <- setdiff(module_genes, pathway_genes)
  
  overlap_list[[i]] <- data.frame(
    anchor_label = xenium_pairs$anchor_label[i],
    pathway = pw,
    module = mod,
    n_pathway_genes = length(pathway_genes),
    n_module_genes = length(module_genes),
    n_overlap_genes = length(overlap_genes),
    n_nonoverlap_module_genes = length(nonoverlap_module_genes),
    overlap_fraction_pathway = length(overlap_genes) / length(pathway_genes),
    overlap_fraction_module = length(overlap_genes) / length(module_genes),
    overlap_genes = paste(overlap_genes, collapse = ";"),
    nonoverlap_module_genes = paste(nonoverlap_module_genes, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

overlap_df <- dplyr::bind_rows(overlap_list)
write.csv(
  overlap_df,
  file.path(results_dir, "xenium_representative_pathway_module_gene_overlap.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------
# Recompute marker-module scores after removing overlapping genes
# ---------------------------------------------------------------

min_nonoverlap_genes <- 2

# Fixed-control implementation
make_nonoverlap_score_fixed_control <- function(
    xen,
    original_col,
    module_genes,
    nonoverlap_genes,
    expr_module,
    cells
) {
  
  module_genes <- intersect(module_genes, rownames(expr_module))
  nonoverlap_genes <- intersect(nonoverlap_genes, rownames(expr_module))
  
  if (length(module_genes) == 0)
    stop("No original module genes found.")
  
  if (length(nonoverlap_genes) == 0)
    stop("No non-overlap module genes remain.")
  
  # Original feature-gene mean
  original_feature_mean <- Matrix::colMeans(
    expr_module[module_genes, cells, drop = FALSE])
  
  # Original AddModuleScore
  original_score <- as.numeric(xen@meta.data[cells, original_col])
  names(original_score) <- cells
  
  # Recover the exact original control contribution
  original_control_mean <- original_feature_mean - original_score
  
  # Feature mean after removing overlap
  nonoverlap_feature_mean <- Matrix::colMeans(
    expr_module[nonoverlap_genes, cells, drop = FALSE])
  
  # Keep control contribution FIXED
  nonoverlap_score <- nonoverlap_feature_mean - original_control_mean
  
  list(
    original_score = original_score,
    nonoverlap_score = nonoverlap_score,
    original_feature_mean = original_feature_mean,
    nonoverlap_feature_mean = nonoverlap_feature_mean,
    original_control_mean = original_control_mean
  )
}

expr_module <- GetAssayData(xen, assay = "Xenium", layer = "data")

fixed_scores <- list()
for (i in seq_len(nrow(xenium_pairs))) {
  
  pw  <- xenium_pairs$pathway[i]
  mod <- xenium_pairs$module[i]
  
  pathway_genes <- intersect(kegg_sets_final[[pw]], rownames(expr_module))
  module_genes <- intersect(marker_sets[[mod]], rownames(expr_module))
  nonoverlap_genes <- setdiff(module_genes, pathway_genes)
  
  if (length(nonoverlap_genes) < min_nonoverlap_genes)
    next
  
  original_col <- paste0("GG_", mod)
  
  res <- make_nonoverlap_score_fixed_control(
    xen = xen,
    original_col = original_col,
    module_genes = module_genes,
    nonoverlap_genes = nonoverlap_genes,
    expr_module = expr_module,
    cells = common_cells
  )
  
  new_col <- paste0("NO_FIXED_", i)
  
  meta[[new_col]] <- res$nonoverlap_score
  
  fixed_scores[[paste0("pair", i)]] <- list(
    score_col = new_col,
    original_genes = module_genes,
    nonoverlap_genes = nonoverlap_genes,
    overlap_genes = intersect(module_genes, pathway_genes)
  )
}


# ---------------------------------------------------------
# Helper function for one pathway and one module-score column
# ---------------------------------------------------------
get_pair_bin_cor <- function(pathway, module_col) {
  
  df <- data.frame(
    cell_id = common_cells,
    bin_id  = meta$bin_id,
    module  = as.numeric(meta[common_cells, module_col]),
    sgsva   = as.numeric(z_sgsva_sub[pathway, common_cells]),
    gsva    = as.numeric(z_gsva_sub[pathway, common_cells]),
    smooth  = as.numeric(z_smooth_sub[pathway, common_cells])
  )
  
  bin_df_pair <- df %>%
    group_by(bin_id) %>%
    summarise(
      module_mean = mean(module, na.rm = TRUE),
      sgsva_mean  = mean(sgsva, na.rm = TRUE),
      gsva_mean   = mean(gsva, na.rm = TRUE),
      smooth_mean = mean(smooth, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    )
  
  out <- data.frame(
    method = c("sGSVA", "GSVA", "GSVA-smoothed"),
    rho = c(
      suppressWarnings(cor(bin_df_pair$module_mean, bin_df_pair$sgsva_mean,
                           method = "spearman", use = "complete.obs")),
      suppressWarnings(cor(bin_df_pair$module_mean, bin_df_pair$gsva_mean,
                           method = "spearman", use = "complete.obs")),
      suppressWarnings(cor(bin_df_pair$module_mean, bin_df_pair$smooth_mean,
                           method = "spearman", use = "complete.obs"))
    )
  )
  
  out$abs_rho <- abs(out$rho)
  out$n_bins <- nrow(bin_df_pair)
  out
}

# ---------------------------
# Original-module and non-overlap-module sensitivity
# ---------------------------
sens_list <- list()

for (i in seq_len(nrow(xenium_pairs))) {
  
  anchor_label <- xenium_pairs$anchor_label[i]
  pw <- xenium_pairs$pathway[i]
  mod <- xenium_pairs$module[i]
  
  # Original marker-module score column, e.g. GG_T cells
  original_col <- paste0("GG_", mod)
  
  # Original module correlation
  original_stats <- get_pair_bin_cor(
    pathway = pw,
    module_col = original_col
  )
  original_stats$anchor_label <- anchor_label
  original_stats$pathway <- pw
  original_stats$module <- mod
  original_stats$module_version <- "Original module"
  
  sens_list[[length(sens_list) + 1]] <- original_stats
  
  # Non-overlap module correlation
  pair_id <- paste0("pair", i)
  
  no_col <- fixed_scores[[pair_id]]$score_col
  nonoverlap_stats <- get_pair_bin_cor(
    pathway = pw,
    module_col = no_col
  )
  nonoverlap_stats$anchor_label <- anchor_label
  nonoverlap_stats$pathway <- pw
  nonoverlap_stats$module <- mod
  nonoverlap_stats$module_version <- "Non-overlap module"
  
  sens_list[[length(sens_list) + 1]] <- nonoverlap_stats
}

sensitivity_long <- dplyr::bind_rows(sens_list) %>%
  dplyr::select(
    anchor_label,
    pathway,
    module,
    module_version,
    method,
    rho,
    abs_rho,
    n_bins
  ) %>%
  dplyr::arrange(anchor_label, method, module_version)

write.csv(
  sensitivity_long,
  file.path(results_dir, "xenium_representative_nonoverlap_sensitivity_long.csv"),
  row.names = FALSE
)


# ---------------------------
# spatial permutation
# ---------------------------
# list of bins that actually have cells
bin_ids <- sort(unique(coords_df$bin_id))
bin_to_cells <- split(coords_df$cell_id, coords_df$bin_id)

# aggregate: pathways × bins matrix of mean sGSVA per bin
bin_scores_sgsva <- sapply(bin_ids, function(b) {
  cells_b <- bin_to_cells[[as.character(b)]]
  rowMeans(z_sgsva_sub[, cells_b, drop = FALSE], na.rm = TRUE)
})

# sapply returns pathways × length(bin_ids)
dim(bin_scores_sgsva)
colnames(bin_scores_sgsva) <- paste0("bin_", bin_ids)
rownames(bin_scores_sgsva) <- rownames(z_sgsva_sub)

# bin centroids
centers <- hcell2xy(hb)
cell_ids_hex <- sort(unique(hb@cell))
bin_df <- data.frame(
  bin_id = cell_ids_hex,
  bx     = centers$x,
  by     = centers$y
)
# subset to used bins in the same order as bin_ids
bin_df_sub <- bin_df[match(bin_ids, bin_df$bin_id), ]
coords_bin <- as.matrix(bin_df_sub[, c("bx", "by")])
rownames(coords_bin) <- paste0("bin_", bin_ids)

# global permutation
t0 <- proc.time()[3]
sgsva_res <- spatial_permutation(bin_scores_sgsva, coords_bin, balanced_sets,
                                 mode = "knn", k=8, sigma_global = 24, 
                                 one_sided = "center_0.5",
                                 perm_methods = c("moran_mc","geary_mc"))
rt <- proc.time()[3]-t0
print(rt) # rt = 242,914 seconds (67.5 hours)
mc_df <- sgsva_res$moran_geary_results$moran_geary_summary
out_csv <- file.path(results_dir, "moran_geary_pathway_summary.csv")
readr::write_csv(mc_df, out_csv)

# export output
mc_export <- mc_df %>%
  dplyr::mutate(
    permutation_scheme = "global",
    one_sided = "center_0.5",
    seed_permutation = 124
  ) %>%
  dplyr::select(-obs_geary_flipped, -z_geary_flipped, 
                -null_mean_geary_flipped, -null_sd_geary_flipped)

# hex-block permutation
block_sizes <- c(150, 200, 300)
block_res <- lapply(block_sizes, function(b) {
  spatial_permutation(bin_scores_sgsva, coords_bin, balanced_sets,
                      mode = "knn", k=8, sigma_global = 24, one_sided = "center_0.5",
                      perm_methods = c("hex_block"), hex_cellsize = b)
})
hex150_df <- block_res[[1]]$hex_block_results$hex_block_summary
readr::write_csv(hex150_df, file.path(results_dir, "block150_pathway_summary.csv"))

hex200_df <- block_res[[2]]$hex_block_results$hex_block_summary
readr::write_csv(hex200_df, file.path(results_dir, "block200_pathway_summary.csv"))

hex300_df <- block_res[[3]]$hex_block_results$hex_block_summary
readr::write_csv(hex300_df, file.path(results_dir, "block300_pathway_summary.csv"))

# export output
hex150_export <- hex150_df %>%
  dplyr::mutate(
    permutation_scheme = "hex-block",
    block_size_um = block_size,
    one_sided = "center_0.5",
    seed_permutation = 124
  ) %>%
  dplyr::select(-obs_geary_flipped, -z_geary_flipped, 
                -null_mean_geary_flipped, -null_sd_geary_flipped, -block_size)
hex200_export <- hex200_df %>%
  dplyr::mutate(
    permutation_scheme = "hex-block",
    block_size_um = block_size,
    one_sided = "center_0.5",
    seed_permutation = 124
  ) %>%
  dplyr::select(-obs_geary_flipped, -z_geary_flipped, 
                -null_mean_geary_flipped, -null_sd_geary_flipped, -block_size)
hex300_export <- hex300_df %>%
  dplyr::mutate(
    permutation_scheme = "hex-block",
    block_size_um = block_size,
    one_sided = "center_0.5",
    seed_permutation = 124
  ) %>%
  dplyr::select(-obs_geary_flipped, -z_geary_flipped, 
                -null_mean_geary_flipped, -null_sd_geary_flipped, -block_size)


# Moran's I/Geary's C effect sizes across sGSVA vs GSVA vs GSVA-smoothed
pathways_sub <- rownames(z_sgsva_sub)
bin_scores_list <- lapply(pathways_sub, function(pw) {
  aggregate_bin_level(pw)   # one bin_scores table per pathway
})
names(bin_scores_list) <- pathways_sub

knn_obj <- spdep::knearneigh(coords_bin, k = 8)
nb <- spdep::knn2nb(knn_obj)
lw_bin <- nb2listw(nb, style = "W")
moran_geary_all <- do.call(rbind, lapply(pathways_sub, function(pw){
  bs <- bin_scores_list[[pw]]
  vals <- list(
    sGSVA = bs$sgsva_mean,
    GSVA = bs$gsva_mean,
    GSVA_smoothed = bs$smooth_mean
  )
  
  out <- lapply(names(vals), function(m) {
    v <- vals[[m]]
    
    mi <- moran.test(v, lw_bin, alternative = "greater")
    gi <- geary.test(v, lw_bin, alternative = "greater")
    
    data.frame(
      pathway = pw,
      method  = m,
      moran_I = mi$estimate[["Moran I statistic"]],
      moran_p = mi$p.value,
      geary_C = gi$estimate[["Geary C statistic"]],
      geary_p = gi$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}))
write.csv(moran_geary_all, file.path(results_dir, "moran_geary_three_methods.csv"))


library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "SData5_panel_coverage")
writeData(wb, "SData5_panel_coverage", stats)
addWorksheet(wb, "SData5_sgsva_diagnostics")
writeData(wb, "SData5_sgsva_diagnostics", orient)
addWorksheet(wb, "SData5_chemokine_sensitivity")
writeData(wb, "SData5_chemokine_sensitivity", chem_anchoring)
addWorksheet(wb, "SData5_all_module_concordance")
writeData(wb, "SData5_all_module_concordance", anchoring_summary)
addWorksheet(wb, "SData5_global_permutation")
writeData(wb, "SData5_global_permutation", mc_export)
addWorksheet(wb, "SData5_block150_permutation")
writeData(wb, "SData5_block150_permutation", hex150_export)
addWorksheet(wb, "SData5_block200_permutation")
writeData(wb, "SData5_block200_permutation", hex200_export)
addWorksheet(wb, "SData5_block300_permutation")
writeData(wb, "SData5_block300_permutation", hex300_export)
addWorksheet(wb, "SData5_runtime")
writeData(wb, "SData5_runtime", runtime_results)
saveWorkbook(wb, file.path(supp_dir, "Supplementary Data 5.xlsx"), overwrite = TRUE)


