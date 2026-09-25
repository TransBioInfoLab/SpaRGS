### Export EnrichMap inputs for breast Visium data

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(dplyr)
  library(tidyr); library(purrr); library(tibble)
  library(msigdbr); library(spdep); library(Matrix)
  library(ggplot2); library(ggspatial); library(patchwork)
  library(scales); library(matrixStats); library(jsonlite)
  library(readr); library(stringr); library(here); library(fs)
})

set.seed(1)

# ---------------------------
# paths and reproducibility
# ---------------------------
project_dir <- "."
results_dir <- file.path(project_dir, "results", "breast_visium")
pilot_dir <- file.path("results", "breast_visium", "enrichmap_pilot")
dir.create(pilot_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_dir, "R/sgsva_score.R"))
source(file.path(project_dir, "R/spatial_permutation.R"))


# ---------------------------
# load data
# ---------------------------
bc_obj <- readRDS(file.path(results_dir, "bc_obj.rds"))

visium_data <- readRDS(file.path(results_dir, "breast_visium_sct_data.rds"))
expr_filtered <- visium_data$expr_filtered
coords_um <- visium_data$coords_um

# ---------------------------
# gene sets
# ---------------------------
kegg_gene_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcollection = "CP:KEGG_LEGACY"
)
gene_sets <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)

# add input metadata
enrichmap_input_metadata <- tibble(
  dataset = "breast_visium",
  expression_matrix = "expr_filtered",
  expression_source = "Filtered SCT data matrix from the same breast Visium object",
  n_spots = ncol(expr_filtered),
  n_genes = nrow(expr_filtered),
  coordinate_source = "coords_um from Visium tissue coordinates converted to microns",
  min_genes_per_gene_set = 2,
  note = "Used for exploratory breast-only EnrichMap supplement"
)

readr::write_csv(
  enrichmap_input_metadata,
  file.path(pilot_dir, "breast_visium_enrichmap_input_metadata.csv")
)


################################################################################
# export breast visium 

normalize_pathway_names <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

export_enrichmap_input <- function(expr_mat, coords, gene_sets,
                                   min_genes = 2, outdir, prefix) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  expr_mat <- as.matrix(expr_mat)
  
  # Make sure coordinates align with expression columns
  common_spots <- intersect(colnames(expr_mat), rownames(coords))
  expr_mat <- expr_mat[, common_spots, drop = FALSE]
  coords <- coords[common_spots, , drop = FALSE]
  
  # AnnData expects spots x genes
  expr_out <- as.data.frame(t(expr_mat))
  expr_out <- tibble::rownames_to_column(expr_out, "spot_id")
  
  coords_out <- as.data.frame(coords[, 1:2, drop = FALSE])
  colnames(coords_out)[1:2] <- c("x", "y")
  coords_out <- tibble::rownames_to_column(coords_out, "spot_id")
  
  # intersect all gene sets with available genes
  gene_sets_filtered <- lapply(gene_sets, function(g) {
    intersect(g, rownames(expr_mat))
  })
  
  # EnrichMap scoring function requires at least two genes from a signature to be present in dataset
  gene_sets_filtered <- gene_sets_filtered[lengths(gene_sets_filtered) >= min_genes]
  
  write_csv(
    expr_out,
    file.path(outdir, paste0(prefix, "_expr_spots_by_genes.csv"))
  )
  
  write_csv(
    coords_out,
    file.path(outdir, paste0(prefix, "_coords.csv"))
  )
  
  write_json(
    gene_sets_filtered,
    file.path(outdir, paste0(prefix, "_gene_sets.json")),
    pretty = TRUE,
    auto_unbox = TRUE
  )
  
  coverage <- tibble(
    pathway_original = names(gene_sets),
    pathway = normalize_pathway_names(names(gene_sets)),
    n_genes_original = lengths(gene_sets),
    n_genes_in_expr = lengths(lapply(gene_sets, function(g) {
      intersect(unique(g), rownames(expr_mat))
    })),
    passed_enrichmap_min_genes = n_genes_in_expr >= min_genes
  )
  
  write_csv(
    coverage,
    file.path(outdir, paste0(prefix, "_enrichmap_gene_set_coverage.csv"))
  )
  
  message("Exported ", nrow(expr_out), " spots.")
  message("Exported ", ncol(expr_out) - 1, " genes.")
  message("Exported ", length(gene_sets_filtered), " gene sets with >= ", min_genes, " genes.")
}

export_enrichmap_input(
  expr_mat = expr_filtered,
  coords = coords_um,
  gene_sets = gene_sets,
  outdir = pilot_dir,
  prefix = "breast_visium"
)
message("Enrichmap scoring runs in 01_enrichmap_breast_visium_pilot.py file")
