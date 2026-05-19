### Real data case3 for Xenium FFPE
suppressPackageStartupMessages ({
  library(Seurat); library(Matrix);library(org.Hs.eg.db)
  library(gridExtra); library(jsonlite); library(ggplot2)
  library(viridis); library(spdep); library(msigdbr)
  library(dplyr); library(purrr) ; library(FNN)
  library(GSVA); library(hexbin)
})

source("R/sgsva_score.R")
source("R/spatial_permutation.R")

set.seed(1)

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
data_dir <- file.path(project_dir, "data", "xenium")
results_dir <- file.path(project_dir, "results", "xenium")
fig_dir <- file.path(project_dir, "results", "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


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
#write.csv(coords, file.path(results_dir, "xenium_coordinates.csv"))

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
kegg_sets_final <- kegg_sets_final[lengths(kegg_sets_final) > 0]

# supplementary S3 for Xenium preprocessing
comb_qc <- (gg_coverage_hist) /
  (cell_qc_plot|gene_qc_plot) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")
ggsave(filename=file.path(fig_dir, "supplement_Fig_S12_Xenium_preprocessing.pdf"),
       plot=comb_qc, width = 10, height = 10)


# ---------------------------------
# computing three pathway scores
# ---------------------------------

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
saveRDS(sgsva_res, file.path(results_dir, "sgsva_score_res.rds"))

sgsva_res <- readRDS(file.path(results_dir, "sgsva_score_res.rds"))
z_sgsva <- sgsva_res$z_score
balanced_sets <- sgsva_res$gene_sets

# GSVA
lib <- colSums(counts_qc)+1
lcp <- log1p(t(t(counts_qc)/lib*1e6)) #logCPM
expr <- as.matrix(lcp)
keep <- rowSums(expr > 0) >= floor(0 * ncol(expr))
expr <- expr[keep, , drop = FALSE]
params <- gsvaParam(
  expr, kegg_sets_final
)
gsva_scores <- gsva(params, verbose = FALSE)
write.csv(gsva_scores, file.path(results_dir, "gsva_score.csv"))

gsva_scores <- read.csv(file.path(results_dir, "gsva_score.csv"))
rownames(gsva_scores) <- gsva_scores$X
gsva_scores$X <- NULL
z_gsva <- z_normalization(as.matrix(gsva_scores))
colnames(z_gsva) <- gsub("\\.","-",colnames(z_gsva))

# GSVA on smoothed
smoothed <- smooth_expr(expr, coords, sigma=24, mode="knn", k=8)
smoothed <- as.matrix(smoothed)
colnames(smoothed) <- colnames(expr)
params <- gsvaParam(
  smoothed, kegg_sets_final
)
gsva_smooth <- gsva(params, verbose = FALSE)
write.csv(gsva_smooth, file.path(results_dir, "gsva_smooth_score.csv"))

gsva_smooth <- read.csv(file.path(results_dir, "gsva_smooth_score.csv"))
rownames(gsva_smooth) <- gsva_smooth$X
gsva_smooth$X <- NULL
z_gsva_smooth <- z_normalization(as.matrix(gsva_smooth))
colnames(z_gsva_smooth) <- gsub("\\.","-",colnames(z_gsva_smooth))


# ---------------------------
# cell type anchoring
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

# run moran.mc/geary.mc
sgsva_res <- spatial_permutation(bin_scores_sgsva, coords_bin, balanced_sets,
                                 mode = "knn", k=8, sigma_global = 24, 
                                 one_sided = "center_0.5",
                                 perm_methods = c("moran_mc","geary_mc"))
mc_df <- sgsva_res$permutation_results$moran_geary_results$moran_geary_summary
out_csv <- file.path(results_dir, "moran_geary_pathway_summary.csv")
readr::write_csv(mc_df, out_csv)

block_sizes <- c(150, 200, 300)
block_res <- lapply(block_sizes, function(b) {
  spatial_permutation(bin_scores_sgsva, coords_bin, balanced_sets,
                      mode = "knn", k=8, sigma_global = 24, one_sided = "center_0.5",
                      perm_methods = c("hex_block"), hex_cellsize = b)
})
hex150_df <- block_res[[1]]$permutation_results$hex_block_results$hex_block_summary
readr::write_csv(hex150_df, file.path(results_dir, "block150_pathway_summary.csv"))

hex200_df <- block_res[[2]]$permutation_results$hex_block_results$hex_block_summary
readr::write_csv(hex200_df, file.path(results_dir, "block200_pathway_summary.csv"))

hex300_df <- block_res[[3]]$permutation_results$hex_block_results$hex_block_summary
readr::write_csv(hex300_df, file.path(results_dir, "block300_pathway_summary.csv"))

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






