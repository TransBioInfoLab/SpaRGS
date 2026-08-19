# Cell-type anchoring for breast cancer Visium
# Methods: sGSVA, GSVA, GSVA-smoothed, and GSDensity

suppressPackageStartupMessages({
  library(Matrix); library(spacexr); library(Seurat)
  library(SeuratObject); library(dplyr); library(tidyr)
  library(purrr); library(ggplot2); library(readr)
  library(stringr); library(tibble); library(here); library(fs)
})

set.seed(1)

# -------------------------------------------------------------------------
# 1. Paths and user-editable settings
# -------------------------------------------------------------------------
project_dir <- "."
data_dir <- file.path(project_dir, "data", "breast_visium")
results_dir <- file.path(project_dir, "results", "breast_visium")
anchoring_dir <- file.path(results_dir, "celltype_anchoring")
table_dir <- file.path(project_dir, "results", "tables")
dir.create(anchoring_dir, recursive = TRUE, showWarnings = FALSE)

bc_object_file <- file.path(results_dir, "bc_obj.rds")

# Update any filename that differs in your project.
score_files <- c(
  "sGSVA" = file.path(results_dir, "sgsva_scores.csv"),
  "GSVA" = file.path(results_dir, "gsva_scores.csv"),
  "GSVA-smoothed" = file.path(results_dir, "gsva_smoothed_scores.csv"),
  "GSDensity" = file.path(results_dir, "gsdensity_scores.csv")
)

# TNBC single-cell reference
sample_id <- "CID4495"
reference_dir <- file.path(data_dir, "GSE176078", "extracted", sample_id)
metadata_file <- file.path(reference_dir, "metadata.csv")
matrix_file <- file.path(reference_dir, "count_matrix_sparse.mtx")
genes_file <- file.path(reference_dir, "count_matrix_genes.tsv")
barcodes_file <- file.path(reference_dir, "count_matrix_barcodes.tsv")

# Saved RCTD output. Once generated, it will be reused.
rctd_file <- file.path(anchoring_dir, "rctd_object.rds")
rctd_weights_file <- file.path(anchoring_dir, "rctd_cell_type_proportions.rds")

max_cores <- 4

# gene sets
kegg_gene_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcollection = "CP:KEGG_LEGACY"
)
gene_sets <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)

# If automatic cell-type matching is ambiguous, enter exact RCTD column names here.
# Leave as "" to use pattern matching.
cell_type_overrides <- c(
  B_cell = "",
  T_cell = "",
  Cancer_epithelial = ""
)

# Prespecified pathway-cell-type anchors.
anchor_pairs <- tribble(
  ~pathway,                                  ~cell_type_key,       ~anchor_label,
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",  "B_cell",             "BCR signaling vs B cells",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",  "T_cell",             "TCR signaling vs T cells",
  "KEGG_GLUTATHIONE_METABOLISM",              "Cancer_epithelial",  "Glutathione metabolism vs cancer epithelial",
  "KEGG_FATTY_ACID_METABOLISM",               "Cancer_epithelial",  "Fatty-acid metabolism vs cancer epithelial"
)

# -------------------------------------------------------------------------
# 2. Helper functions
# -------------------------------------------------------------------------
normalize_pathway_names <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  # Remove the KEGG prefix so methods using KEGG-prefixed and unprefixed
  # pathway names can still be matched.
  x <- sub("^KEGG_", "", x)
  x
}

collapse_duplicate_pathways <- function(mat) {
  if (!anyDuplicated(rownames(mat))) {
    return(mat)
  }

  pathway_index <- split(seq_len(nrow(mat)), rownames(mat))
  out <- lapply(pathway_index, function(idx) {
    colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- names(pathway_index)
  colnames(out) <- colnames(mat)
  out
}

read_score_matrix <- function(file, spot_ids, method_name) {
  if (!file.exists(file)) {
    stop("Score file not found for ", method_name, ": ", file)
  }

  mat <- read.csv(file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  mat <- as.matrix(mat)
  suppressWarnings(storage.mode(mat) <- "numeric")

  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(method_name, " score matrix must have pathway and spot names.")
  }

  row_overlap <- sum(rownames(mat) %in% spot_ids)
  col_overlap <- sum(colnames(mat) %in% spot_ids)

  # Convert to pathway x spot orientation.
  if (row_overlap > col_overlap) {
    mat <- t(mat)
    col_overlap <- sum(colnames(mat) %in% spot_ids)
  }

  original_names <- rownames(mat)
  rownames(mat) <- normalize_pathway_names(original_names)
  mat <- collapse_duplicate_pathways(mat)

  message(
    method_name, ": ", nrow(mat), " pathways; ",
    col_overlap, " matching spots."
  )

  mat
}

get_spatial_coordinates <- function(seurat_obj, spot_ids) {
  coords_raw <- as.data.frame(GetTissueCoordinates(seurat_obj))

  if (!all(spot_ids %in% rownames(coords_raw))) {
    missing_ids <- setdiff(spot_ids, rownames(coords_raw))
    stop(
      "Coordinates are missing for ", length(missing_ids),
      " spatial-count barcodes."
    )
  }

  coords_raw <- coords_raw[spot_ids, , drop = FALSE]

  if (all(c("x", "y") %in% colnames(coords_raw))) {
    coords <- coords_raw[, c("x", "y"), drop = FALSE]
  } else if (all(c("imagecol", "imagerow") %in% colnames(coords_raw))) {
    coords <- coords_raw[, c("imagecol", "imagerow"), drop = FALSE]
    colnames(coords) <- c("x", "y")
  } else {
    stop(
      "Could not identify coordinate columns. Available columns: ",
      paste(colnames(coords_raw), collapse = ", ")
    )
  }
  coords
}

resolve_cell_type <- function(key, available_names, overrides) {
  override <- overrides[[key]]

  if (!is.null(override) && nzchar(override)) {
    if (!override %in% available_names) {
      stop(
        "Cell-type override '", override, "' for ", key,
        " is not present in the RCTD weights. Available names:\n",
        paste(available_names, collapse = ", ")
      )
    }
    return(override)
  }

  patterns <- switch(
    key,
    B_cell = c("^B[ _-]?CELLS?$", "B[ _-]?CELL"),
    T_cell = c("^T[ _-]?CELLS?$", "T[ _-]?CELL"),
    Cancer_epithelial = c(
      "CANCER.*EPITHELIAL",
      "EPITHELIAL.*CANCER",
      "MALIGNANT",
      "TUMOU?R.*EPITHELIAL"
    ),
    stop("Unknown cell_type_key: ", key)
  )

  hits <- unique(unlist(lapply(patterns, function(pattern) {
    grep(pattern, available_names, ignore.case = TRUE, value = TRUE)
  })))

  if (length(hits) != 1L) {
    stop(
      "Could not uniquely resolve ", key, ". Candidate matches: ",
      ifelse(length(hits) == 0L, "none", paste(hits, collapse = ", ")),
      "\nAvailable RCTD cell types:\n",
      paste(available_names, collapse = ", "),
      "\nSet an exact name in cell_type_overrides."
    )
  }

  hits
}

safe_spearman <- function(x, y, min_n = 10L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < min_n || sd(x) == 0 || sd(y) == 0) {
    return(tibble(
      n = length(x),
      rho = NA_real_,
      p_value = NA_real_
    ))
  }

  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))

  tibble(
    n = length(x),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}

z_standardize <- function(x) {
  if (all(!is.finite(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

# -------------------------------------------------------------------------
# 3. Load Visium object
# -------------------------------------------------------------------------
if (!file.exists(bc_object_file)) {
  stop("Breast Visium Seurat object not found: ", bc_object_file)
}

bc_obj <- readRDS(bc_object_file)

spatial_counts <- GetAssayData(
  bc_obj,
  assay = "Spatial",
  layer = "counts"
)

spot_ids <- colnames(spatial_counts)
coords <- get_spatial_coordinates(bc_obj, spot_ids)
stopifnot(identical(rownames(coords), spot_ids))

# -------------------------------------------------------------------------
# 4. Run or load RCTD
# -------------------------------------------------------------------------
if (file.exists(rctd_weights_file)) {
  message("Loading existing RCTD cell-type proportions.")
  cell_type_proportions <- readRDS(rctd_weights_file)
} else {
  required_reference_files <- c(
    metadata_file, matrix_file, genes_file, barcodes_file
  )
  missing_reference_files <- required_reference_files[
    !file.exists(required_reference_files)
  ]

  if (length(missing_reference_files) > 0L) {
    stop(
      "Missing single-cell reference files:\n",
      paste(missing_reference_files, collapse = "\n")
    )
  }

  message("Loading single-cell reference counts and metadata.")
  sc_counts <- readMM(matrix_file)
  sc_barcodes <- read.delim(barcodes_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  sc_genes <- read.delim(genes_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  sc_metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  rownames(sc_counts) <- sc_genes
  colnames(sc_counts) <- sc_barcodes

  if (!"celltype_major" %in% colnames(sc_metadata)) {
    stop("Single-cell metadata does not contain 'celltype_major'.")
  }

  common_reference_cells <- intersect(
    colnames(sc_counts),
    rownames(sc_metadata)
  )

  if (length(common_reference_cells) == 0L) {
    stop("No reference-cell barcodes match between counts and metadata.")
  }

  sc_counts <- sc_counts[, common_reference_cells, drop = FALSE]
  cell_types <- sc_metadata[common_reference_cells, "celltype_major", drop = TRUE]

  keep_reference_cells <- !is.na(cell_types) & nzchar(cell_types)
  sc_counts <- sc_counts[, keep_reference_cells, drop = FALSE]
  cell_types <- factor(cell_types[keep_reference_cells])
  names(cell_types) <- colnames(sc_counts)

  # RCTD requires raw counts. SCTransform, PCA, clustering, and UMAP are
  # not needed for constructing the reference.
  reference <- Reference(
    sc_counts,
    cell_types,
    nUMI = Matrix::colSums(sc_counts)
  )

  puck <- SpatialRNA(
    coords,
    spatial_counts,
    nUMI = Matrix::colSums(spatial_counts)
  )

  message("Running RCTD.")
  rctd <- create.RCTD(
    puck,
    reference,
    max_cores = max_cores
  )
  rctd <- run.RCTD(rctd)

  saveRDS(rctd, rctd_file)

  cell_type_proportions <- as.data.frame(rctd@results$weights)
  saveRDS(cell_type_proportions, rctd_weights_file)
  write_csv(
    rownames_to_column(cell_type_proportions, "barcode"),
    file.path(anchoring_dir, "rctd_cell_type_proportions.csv")
  )
}

cell_type_proportions <- as.data.frame(cell_type_proportions)

if (is.null(rownames(cell_type_proportions))) {
  stop("RCTD cell-type proportions must have spot barcodes as row names.")
}

message(
  "RCTD cell types: ",
  paste(colnames(cell_type_proportions), collapse = ", ")
)

# -------------------------------------------------------------------------
# 5. Load all four pathway-score matrices
# -------------------------------------------------------------------------
score_list <- imap(score_files, function(file, method_name) {
  read_score_matrix(file, spot_ids, method_name)
})

# Use the same spots for every method and RCTD.
common_spots <- rownames(cell_type_proportions)
for (mat in score_list) {
  common_spots <- common_spots[common_spots %in% colnames(mat)]
}
common_spots <- common_spots[common_spots %in% spot_ids]

message("Common spots used for all methods: ", length(common_spots))

cell_type_proportions <- cell_type_proportions[common_spots,,drop = FALSE]
score_list <- lapply(score_list, function(mat) {
  mat[, common_spots, drop = FALSE]
})

stopifnot(all(vapply(
  score_list,
  function(mat) identical(colnames(mat), common_spots),
  logical(1)
)))

# Strict shared-pathway set for paired head-to-head summaries.
shared_pathways <- Reduce(intersect, lapply(score_list, rownames))

write_lines(
  shared_pathways,
  file.path(anchoring_dir, "shared_pathways_all_methods.txt")
)

coverage_table <- imap_dfr(score_list, function(mat, method_name) {
  tibble(
    method = method_name,
    pathways_available = nrow(mat),
    common_spots = ncol(mat),
    pathways_shared_by_all_methods = length(shared_pathways)
  )
})

write_csv(
  coverage_table,
  file.path(anchoring_dir, "anchoring_method_coverage.csv")
)

# -------------------------------------------------------------------------
# 6. Resolve RCTD cell-type names for the prespecified anchors
# -------------------------------------------------------------------------
resolved_cell_types <- c(
  B_cell = resolve_cell_type(
    "B_cell",
    colnames(cell_type_proportions),
    cell_type_overrides
  ),
  T_cell = resolve_cell_type(
    "T_cell",
    colnames(cell_type_proportions),
    cell_type_overrides
  ),
  Cancer_epithelial = resolve_cell_type(
    "Cancer_epithelial",
    colnames(cell_type_proportions),
    cell_type_overrides
  )
)

anchor_pairs <- anchor_pairs %>%
  mutate(
    pathway_norm = normalize_pathway_names(pathway),
    cell_type = unname(resolved_cell_types[cell_type_key])
  )

write_csv(
  anchor_pairs,
  file.path(anchoring_dir, "prespecified_anchor_definitions.csv")
)

# -------------------------------------------------------------------------
# 7. Compute all pathway-cell-type correlations
# -------------------------------------------------------------------------
all_available_results <- imap_dfr(score_list, function(mat, method_name) {
  pathway_names <- rownames(mat)
  cell_type_names <- colnames(cell_type_proportions)

  result_list <- vector(
    "list",
    length(pathway_names) * length(cell_type_names)
  )
  counter <- 1L

  for (pathway_name in pathway_names) {
    score_vector <- as.numeric(mat[pathway_name, common_spots])

    for (cell_type_name in cell_type_names) {
      abundance_vector <- cell_type_proportions[common_spots, cell_type_name]

      test_result <- safe_spearman(score_vector, abundance_vector)

      result_list[[counter]] <- tibble(
        method = method_name,
        pathway = pathway_name,
        cell_type = cell_type_name,
        n = test_result$n,
        rho = test_result$rho,
        abs_rho = abs(test_result$rho),
        p_value = test_result$p_value
      )
      counter <- counter + 1L
    }
  }

  bind_rows(result_list)
}) %>%
  group_by(method) %>%
  mutate(fdr_within_method = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write_csv(
  all_available_results,
  file.path(
    anchoring_dir,
    "anchoring_all_available_pathways_celltypes.csv"
  )
)

shared_results <- all_available_results %>%
  filter(pathway %in% shared_pathways) %>%
  group_by(method) %>%
  mutate(fdr_shared_within_method = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write_csv(
  shared_results,
  file.path(
    anchoring_dir,
    "anchoring_shared_pathways_all_methods.csv"
  )
)

# -------------------------------------------------------------------------
# 8. Compute the prespecified biological anchors for every method
# -------------------------------------------------------------------------
prespecified_results <- imap_dfr(score_list, function(mat, method_name) {
  pmap_dfr(anchor_pairs, function(
      pathway, cell_type_key, anchor_label, pathway_norm, cell_type
  ) {
    if (!pathway_norm %in% rownames(mat)) {
      return(tibble(
        method = method_name,
        pathway = pathway_norm,
        cell_type_key = cell_type_key,
        cell_type = cell_type,
        anchor_label = anchor_label,
        status = "pathway_not_available",
        n = NA_integer_,
        rho = NA_real_,
        abs_rho = NA_real_,
        p_value = NA_real_
      ))
    }

    score_vector <- as.numeric(mat[pathway_norm, common_spots])
    abundance_vector <- cell_type_proportions[common_spots, cell_type]

    test_result <- safe_spearman(score_vector, abundance_vector)

    tibble(
      method = method_name,
      pathway = pathway_norm,
      cell_type_key = cell_type_key,
      cell_type = cell_type,
      anchor_label = anchor_label,
      status = "evaluated",
      n = test_result$n,
      rho = test_result$rho,
      abs_rho = abs(test_result$rho),
      p_value = test_result$p_value
    )
  })
}) %>%
  group_by(method) %>%
  mutate(fdr_prespecified_within_method = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write_csv(
  prespecified_results,
  file.path(anchoring_dir, "anchoring_prespecified_pairs.csv")
)

# -------------------------------------------------------------------------
# 9. Figures
# -------------------------------------------------------------------------
anchoring_plot_df <- prespecified_results %>%
  mutate(
    method = factor(
      method,
      levels = c("GSVA", "GSVA-smoothed", "GSDensity", "sGSVA")
    ),
    anchor_label = factor(
      anchor_label,
      levels = c(
        "BCR signaling vs B cells",
        "TCR signaling vs T cells",
        "Glutathione metabolism vs cancer epithelial",
        "Fatty-acid metabolism vs cancer epithelial"
      )
    ),
    signed_label = sprintf("%+.2f", rho)
  )

p_anchor <- ggplot(
  anchoring_plot_df,
  aes(x = method, y = anchor_label, fill = abs_rho)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = signed_label), size = 3.8) +
  scale_fill_gradient(
    limits = c(0, 1),
    low = "white",
    high = "steelblue",
    name = expression("|Spearman " * rho * "|")
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Magnitude of pathway-cell-type spatial concordance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(anchoring_dir, "main_anchoring_prespecified_heatmap.pdf"),
  p_anchor,
  width = 9,
  height = 5.5
)


# Faceted scatterplot using standardized pathway scores only for display.
scatter_data <- imap_dfr(score_list, function(mat, method_name) {
  pmap_dfr(anchor_pairs, function(
      pathway,
      cell_type_key,
      anchor_label,
      pathway_norm,
      cell_type
  ) {
    if (!pathway_norm %in% rownames(mat)) {
      return(tibble())
    }

    tibble(
      barcode = common_spots,
      method = method_name,
      anchor_label = anchor_label,
      abundance = as.numeric(
        cell_type_proportions[common_spots, cell_type]
      ),
      pathway_z = z_standardize(
        as.numeric(mat[pathway_norm, common_spots])
      )
    )
  })
})

scatter_annotations <- prespecified_results %>%
  filter(status == "evaluated") %>%
  mutate(
    annotation = paste0(
      "rho = ", sprintf("%.2f", rho),
      "\nFDR = ", formatC(
        fdr_prespecified_within_method,
        format = "g",
        digits = 2
      )
    )
  )

p_scatter <- ggplot(
  scatter_data,
  aes(x = abundance, y = pathway_z)
) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_grid(
    anchor_label ~ method,
    scales = "free_x"
  ) +
  labs(
    x = "RCTD-estimated cell-type proportion",
    y = "Standardized pathway score",
    title = "Pathway-cell-type anchoring by method"
  ) +
  theme_bw()

ggsave(
  file.path(anchoring_dir, "anchoring_prespecified_scatterplots.pdf"),
  p_scatter,
  width = 14,
  height = 11
)

# -------------------------------------------------------------------------
# 10. Save a compact summary
# -------------------------------------------------------------------------
saveRDS(
  list(
    score_files = score_files,
    common_spots = common_spots,
    shared_pathways = shared_pathways,
    resolved_cell_types = resolved_cell_types,
    coverage = coverage_table,
    all_available_results = all_available_results,
    shared_results = shared_results,
    prespecified_results = prespecified_results,
    session_info = capture.output(sessionInfo())
  ),
  file.path(anchoring_dir, "breast_celltype_anchoring_results.rds")
)

message("Cell-type anchoring analysis completed.")
message("Outputs written to: ", anchoring_dir)


# -------------------------------------------------------------------------
# 11. checking orientation
# -------------------------------------------------------------------------
# check orientation
expr_matrix <- as.matrix(bc_obj[["SCT"]]@data)
nonzero_counts <- rowSums(expr_matrix > 0)
expr_filtered <- expr_matrix[nonzero_counts >= floor(0.05 * ncol(expr_matrix)), ]

pathways <- anchor_pairs$pathway
sgsva_res <- readRDS(file.path(results_dir, "sgsva_score_res.rds"))
sgsva_scores <- sgsva_res$scores
balanced_sets <- sgsva_res$gene_sets

check_final_orientation <- function(
    pathway,
    sgsva_scores,
    expr_matrix,
    gene_sets
) {
  genes <- intersect(gene_sets[[pathway]], rownames(expr_matrix))
  
  pathway_mean <- colMeans(expr_matrix[genes, , drop = FALSE],na.rm = TRUE)
  
  common_spots <- intersect(colnames(sgsva_scores), names(pathway_mean))
  
  cor(
    as.numeric(sgsva_scores[pathway, common_spots]),
    pathway_mean[common_spots],
    method = "spearman",
    use = "complete.obs"
  )
}

res_orien <- sapply(pathways, function(p) 
  check_final_orientation(p, sgsva_scores, expr_filtered, gene_sets))

diag_tab <- data.frame(
  pathway = gsub("KEGG_", "", names(res_orien)), 
  cor_expr = as.numeric(res_orien),
  signed_rho = NA_real_, abs_rho = NA_real_,
  pos_genes = NA_real_, neg_genes = NA_real_
)
for (i in seq_len(length(res_orien))) {
  pathway <- diag_tab$pathway[i]
  diag_tab$signed_rho[i] <- prespecified_results$rho[prespecified_results$method == "sGSVA" 
                                         & prespecified_results$pathway == pathway]
  diag_tab$abs_rho[i] <- prespecified_results$abs_rho[prespecified_results$method == "sGSVA" 
                                                  & prespecified_results$pathway == pathway]
  diag_tab$pos_genes[i] <- length(balanced_sets[[names(res_orien)[i]]]$U)
  diag_tab$neg_genes[i] <- length(balanced_sets[[names(res_orien)[i]]]$D)
}
