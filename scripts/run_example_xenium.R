# run a small example for sGSVA
## use a small subsets from Xenium FFPE dataset

suppressPackageStartupMessages ({
  library(Seurat); library(Matrix); library(jsonlite); library(ggplot2)
  library(spdep); library(msigdbr); library(dplyr); 
  library(purrr) ; library(FNN); library(hexbin)
})

source("R/sgsva_score.R")
source("R/spatial_permutation.R")

set.seed(1)

# ---------------------------
# Paths and settings
# ---------------------------
project_dir <- "."
data_dir    <- file.path(project_dir, "data", "xenium")
out_dir     <- file.path(project_dir, "results", "example")

n_cells_example <- 5000
bin_size_um     <- 40
k_knn           <- 8
cor_min         <- 0.10

target_pathways <- c(
  "KEGG_CELL_CYCLE",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_ERBB_SIGNALING_PATHWAY",
  "KEGG_FOCAL_ADHESION"
)


# ---------------------------
# Helper functions
# ---------------------------
z_scale_rows <- function(mat) {
  mat <- as.matrix(mat)
  t(scale(t(mat)))
}

panel_filter_sets <- function(counts, sets, min_genes = 15, min_coverage = 0.30) {
  panel_genes <- rownames(counts)
  
  stats_list <- vector("list", length(sets))
  names(stats_list) <- names(sets)
  
  sets_on_panel <- vector("list", length(sets))
  names(sets_on_panel) <- names(sets)
  
  for (nm in names(sets)) {
    gs_all <- unique(sets[[nm]])
    n_total <- length(gs_all)
    
    matched <- intersect(gs_all, panel_genes)
    n_panel <- length(matched)
    
    sets_on_panel[[nm]] <- matched
    
    stats_list[[nm]] <- data.frame(
      pathway = nm,
      n_total_genes = n_total,
      n_panel_genes = n_panel,
      coverage = ifelse(n_total > 0, n_panel / n_total, NA_real_),
      stringsAsFactors = FALSE
    )
  }
  
  stats <- bind_rows(stats_list) %>%
    mutate(keep = (n_panel_genes >= min_genes) | (coverage >= min_coverage))
  
  sets_filtered <- sets_on_panel[stats$pathway[stats$keep]]
  
  list(
    sets_filtered = sets_filtered,
    stats = stats,
    sets_on_panel = sets_on_panel
  )
}
  
aggregate_to_bins <- function(score_mat, coords_df, bin_size_um = 40) {
  xrange <- range(coords_df$x, na.rm = TRUE)
  n_x <- max(1, ceiling((xrange[2] - xrange[1]) / bin_size_um))
  
  hb <- hexbin(
    x = coords_df$x,
    y = coords_df$y,
    xbins = n_x,
    IDs = TRUE
  )
  
  coords_df$bin_id <- hb@cID
  bin_centers <- hcell2xy(hb)
  used_bins <- sort(unique(hb@cell))
  
  bin_df <- data.frame(
    bin_id = used_bins,
    bx = bin_centers$x,
    by = bin_centers$y
  )
  
  bin_to_cells <- split(coords_df$cell_id, coords_df$bin_id)
  
  bin_score_mat <- sapply(used_bins, function(b) {
    cells_b <- bin_to_cells[[as.character(b)]]
    rowMeans(score_mat[, cells_b, drop = FALSE], na.rm = TRUE)
  })
  
  colnames(bin_score_mat) <- paste0("bin_", used_bins)
  rownames(bin_score_mat) <- rownames(score_mat)
  
  coords_bin <- as.matrix(bin_df[match(used_bins, bin_df$bin_id), c("bx", "by")])
  rownames(coords_bin) <- colnames(bin_score_mat)
  
  list(
    bin_scores = bin_score_mat,
    bin_coords = coords_bin,
    cell_bins = coords_df,
    bin_df = bin_df
  )
}


# ---------------------------
# Check inputs
# ---------------------------
if (!dir.exists(data_dir)) {
  stop(
    "Xenium data directory not found: ", data_dir, "\n",
    "Please download the public 10x Genomics Xenium FFPE breast cancer dataset ",
    "and place it in data/xenium/."
  )
}

message("Loading Xenium data from: ", data_dir)


# ---------------------------
# Load Xenium data
# ---------------------------
xe <- ReadXenium(
  data.dir = data_dir,
  outs = "matrix",
  type = "centroids",
  mols.qv.threshold = 20
)

counts <- xe$matrix$`Gene Expression`
centroids <- xe$centroids

if (!all(c("x", "y") %in% colnames(centroids))) {
  stop("Centroid table must contain columns named x and y.")
}

if (!"cell" %in% colnames(centroids)) {
  centroids$cell <- rownames(centroids)
}


# ---------------------------
# KEGG gene sets and panel-aware filtering
# ---------------------------
kegg_tbl <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcollection = "CP:KEGG_LEGACY"
)
gene_sets <- split(kegg_tbl$gene_symbol, kegg_tbl$gs_name)

pf <- panel_filter_sets(
  counts = counts,
  sets = gene_sets,
  min_genes = 15,
  min_coverage = 0.30
)

# Keep target pathways if available; otherwise use first retained pathways.
sets_available <- pf$sets_filtered
target_keep <- intersect(target_pathways, names(sets_available))

if (length(target_keep) == 0) {
  warning("Target pathways were not retained after panel filtering. Using first 3 retained KEGG pathways.")
  target_keep <- head(names(sets_available), 3)
}

gene_sets_example <- sets_available[target_keep]
genes_used <- unique(unlist(gene_sets_example))
counts_sub <- counts[intersect(rownames(counts), genes_used), , drop = FALSE]


# ---------------------------
# Cell and gene QC
# ---------------------------
cell_total_counts <- Matrix::colSums(counts_sub)
cell_n_features   <- Matrix::colSums(counts_sub > 0)

keep_cells <- (cell_total_counts >= 30) & (cell_n_features >= 15)

counts_qc <- counts_sub[, keep_cells, drop = FALSE]
centroids_qc <- centroids[match(colnames(counts_qc), centroids$cell), , drop = FALSE]

gene_n_cells <- Matrix::rowSums(counts_qc > 0)
keep_genes <- gene_n_cells >= 20

counts_qc <- counts_qc[keep_genes, , drop = FALSE]
gene_sets_example <- lapply(gene_sets_example, function(gs) intersect(gs, rownames(counts_qc)))
gene_sets_example <- gene_sets_example[lengths(gene_sets_example) > 0]

# use a contiguous spatial subset
select_spatial_window_cells <- function(centroids, n = 5000, seed = 1) {
  set.seed(seed)
  center_id <- sample(seq_len(nrow(centroids)), 1)
  center <- as.numeric(centroids[center_id, c("x", "y")])
  d <- sqrt((centroids$x - center[1])^2 + (centroids$y - center[2])^2)
  centroids$cell[order(d)][seq_len(min(n, nrow(centroids)))]
}

n_sample <- min(n_cells_example, ncol(counts_qc))
sample_cells <- select_spatial_window_cells(
  centroids_qc,
  n = n_cells_example
)

counts_ex <- counts_qc[, sample_cells, drop = FALSE]
centroids_ex <- centroids_qc[match(sample_cells, centroids_qc$cell), , drop = FALSE]
coords <- data.frame(
  x = centroids_ex$x,
  y = centroids_ex$y,
  row.names = sample_cells
)

# recompute gene sets after sampling
gene_n_cells_ex <- Matrix::rowSums(counts_ex > 0)
counts_ex <- counts_ex[gene_n_cells_ex > 0, , drop = FALSE]

gene_sets_example <- lapply(gene_sets_example, function(gs) intersect(gs, rownames(counts_ex)))
gene_sets_example <- gene_sets_example[lengths(gene_sets_example) >= 2]

if (length(gene_sets_example) == 0) {
  stop("No example pathways retained after spatial subsetting. Increase n_cells_example or choose a different tissue window.")
}


# ---------------------------
# Normalize expression
# ---------------------------
lib <- Matrix::colSums(counts_ex) + 1
expr <- log1p(t(t(counts_ex) / lib * 1e6))
expr <- as.matrix(expr)


# ---------------------------
# Run sGSVA
# ---------------------------
message("Running minimal sGSVA Xenium example on ", ncol(expr), " cells...")

optimal_sigma <- round(estimate_sigma_knn(coords, k=8))
sgsva_res <- run_sgsva(
  expr, coords, gene_sets_example,
  sigma_global = optimal_sigma,
  cor_min = cor_min,
  mode = "knn",
  k = k_knn,
  one_sided = "center_0.5"
)
z_sgsva <- sgsva_res$z_score
write.csv(
  z_sgsva,
  file.path(out_dir, "minimal_example_xenium_sgsva_zscores.csv")
)


# ---------------------------
# Aggregate to 40 um hex-bin pseudospots
# ---------------------------
coords_df <- data.frame(
  cell_id = colnames(z_sgsva),
  x = coords[colnames(z_sgsva), "x"],
  y = coords[colnames(z_sgsva), "y"],
  stringsAsFactors = FALSE
)

bin_res <- aggregate_to_bins(
  score_mat = z_sgsva,
  coords_df = coords_df,
  bin_size_um = bin_size_um
)

write.csv(bin_res$bin_scores,file.path(out_dir, "minimal_example_xenium_bin_level_sgsva_scores.csv"))

write.csv(
  data.frame(bin_id = rownames(bin_res$bin_coords), bin_res$bin_coords),
  file.path(out_dir, "minimal_example_xenium_bin_coordinates.csv"),
  row.names = FALSE
)


# ----------------------------------
# spatial permutation testing on bin-level scores
# ----------------------------------
# run moran.mc/geary.mc
perm_res <- spatial_permutation(z_sgsva = bin_res$bin_scores, 
                                coords = bin_res$bin_coords, 
                                balanced_sets = sgsva_res$gene_sets,
                                 mode = "knn", k=8, sigma_global = optimal_sigma, 
                                 cor_min = cor_min, one_sided = "center_0.5",
                                 perm_methods = c("moran_mc","geary_mc"))
mc_df <- perm_res$permutation_results$moran_geary_results$moran_geary_summary
write.csv(mc_df, file.path(out_dir, "minimal_example_xenium_hex_block_results_csv"), row.names = FALSE)

mini_res <- list(
  sgsva_results = sgsva_res, hex_block_permutation = perm_res
  )

saveRDS(mini_res,file.path(out_dir, "minimal_example_xenium_all_result.rds"))
