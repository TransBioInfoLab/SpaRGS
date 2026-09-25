#!/usr/bin/env Rscript

# ============================================================
# Minimal reproducibility example for sGSVA
#
# This script is a lightweight smoke test for a clean sGSVA
# installation. To improve cross-platform reproducibility, it
# uses a frozen small simulated input stored in the repository
# rather than regenerating the Gaussian field and counts at run time.
#
# IMPORTANT:
# - This is NOT the manuscript benchmark.
# - It is intended only to verify that the released sGSVA code
#   runs end-to-end from a clean clone on fixed inputs.
# - Run from the root of the sGSVA repository:
#
#     Rscript scripts/run_minimal_example.R
#
# Expected repository files:
#     R/sgsva_score.R
#     R/spatial_permutation.R
#     data/example/minimal_example_input.rds
# ============================================================

set.seed(1)

# ------------------------------------------------------------
# 1. Check repository location and packages
# ------------------------------------------------------------

input_file <- file.path("data", "example", "minimal_example_input.rds")

required_files <- c(
  file.path("R", "sgsva_score.R"),
  file.path("R", "spatial_permutation.R"),
  input_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required repository files were not found:\n  ",
    paste(missing_files, collapse = "\n  "),
    "\n\nRun this script from the root of the sGSVA repository:\n",
    "  Rscript scripts/run_minimal_example.R\n\n",
    "The frozen example input should be committed at:\n  ",
    input_file
  )
}

required_packages <- c(
  "Matrix",
  "matrixStats",
  "FNN",
  "spdep"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running the example."
  )
}

suppressPackageStartupMessages({
  library(Matrix)
  library(spdep)
})

source(file.path("R", "sgsva_score.R"))
source(file.path("R", "spatial_permutation.R"))

out_dir <- file.path("results", "example")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load the frozen minimal-example input
# ------------------------------------------------------------

example_input <- readRDS(input_file)

required_input_names <- c(
  "counts",
  "coords",
  "gene_sets",
  "fields",
  "per_path",
  "selected_pathways",
  "parameters"
)

missing_input_names <- setdiff(required_input_names, names(example_input))

if (length(missing_input_names) > 0L) {
  stop(
    "The frozen example input is missing required component(s): ",
    paste(missing_input_names, collapse = ", ")
  )
}

counts <- example_input$counts
coords <- as.data.frame(example_input$coords)
gene_sets <- example_input$gene_sets
truth_fields <- example_input$fields
true_sets <- example_input$per_path
sel_paths <- example_input$selected_pathways
params <- example_input$parameters

# Preserve spot names and verify alignment.
if (is.null(rownames(coords))) {
  rownames(coords) <- colnames(counts)
}

if (is.null(colnames(counts)) || is.null(rownames(counts))) {
  stop("Frozen count matrix must contain gene and spot names.")
}

if (!identical(rownames(coords), colnames(counts))) {
  stop(
    "Spot ordering mismatch between frozen coordinates and count matrix."
  )
}

ell <- params$ell
sigma_f <- params$sigma_f
up <- params$up
down <- params$down
sigma_global <- params$sigma_global

if (is.null(sigma_global) || !is.finite(sigma_global)) {
  sigma_global <- floor(ell / sqrt(2))
}

input_md5 <- unname(tools::md5sum(input_file))

message("Running fixed-input minimal sGSVA example...")
message("Frozen input: ", input_file)
message("Input MD5: ", input_md5)
message("Spatial locations: ", nrow(coords))
message("Input genes: ", nrow(counts))
message("Input KEGG pathways: ", length(sel_paths))
message("Scenario used to generate frozen input: ell = ", ell,
        " um; sigma_f = ", sigma_f)
message("Balanced pathway design: ", up, " up / ", down, " down")
message("sGSVA smoothing bandwidth: ", sigma_global, " um")

# Save the exact frozen pathway membership for auditability.
example_gene_sets_df <- do.call(
  rbind,
  lapply(sel_paths, function(p) {
    data.frame(
      pathway = p,
      gene = gene_sets[[p]],
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  example_gene_sets_df,
  file.path(out_dir, "minimal_example_kegg_gene_sets.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Normalize as logCPM, matching the simulation workflow
# ------------------------------------------------------------

lib <- Matrix::colSums(counts) + 1
lcp <- log1p(t(t(counts) / lib * 1e6))
expr <- as.matrix(lcp)
keep <- rowSums(expr > 0) >= floor(0.05 * ncol(expr))
expr <- expr[keep, , drop = FALSE]

message("Genes after filtering: ", nrow(expr))

# ------------------------------------------------------------
# 4. Run the released sGSVA implementation
# ------------------------------------------------------------

sgsva_res <- run_sgsva(
  expr = expr,
  coords = coords,
  sets = gene_sets[sel_paths],
  mode = "delaunay",
  sigma_global = sigma_global,
  weight_clip = 0.99,
  cor_min = 0.10,
  min_genes = 15,
  max_n = 25,
  one_sided = "center_0.5"
)

z_sgsva <- as.matrix(sgsva_res$z_score)

# ------------------------------------------------------------
# 5. Validate the score output
# ------------------------------------------------------------

if (nrow(z_sgsva) == 0L) {
  stop(
    "No pathways were retained by data-driven sGSVA in the ",
    "minimal example."
  )
}

if (ncol(z_sgsva) != nrow(coords)) {
  stop(
    "Score matrix columns do not match the number of spatial locations."
  )
}

if (!all(is.finite(z_sgsva))) {
  stop("Non-finite values were found in the sGSVA score matrix.")
}

message(
  "sGSVA scoring completed successfully: ",
  nrow(z_sgsva),
  " retained pathways x ",
  ncol(z_sgsva),
  " spatial locations."
)

# ------------------------------------------------------------
# 6. Compare inferred maps with the known simulated fields
# ------------------------------------------------------------

retained_paths <- intersect(rownames(z_sgsva), names(truth_fields))

truth_summary <- do.call(
  rbind,
  lapply(retained_paths, function(p) {
    rho_signed <- suppressWarnings(
      stats::cor(
        as.numeric(z_sgsva[p, ]),
        as.numeric(truth_fields[[p]]),
        method = "spearman"
      )
    )

    data.frame(
      pathway = p,
      signed_spearman = rho_signed,
      absolute_spearman = abs(rho_signed),
      stringsAsFactors = FALSE
    )
  })
)

if (is.null(truth_summary) || nrow(truth_summary) == 0L) {
  stop("No retained pathways could be matched to simulated truth.")
}

# ------------------------------------------------------------
# 7. Lightweight spatial-inference smoke test
# ------------------------------------------------------------
#
# Use only 100 permutations so the clean-install example remains
# fast. These are NOT the manuscript inference settings.

perm_res <- spatial_permutation(
  z_sgsva = z_sgsva,
  coords = coords,
  balanced_sets = sgsva_res$gene_sets,
  mode = "delaunay",
  sigma_global = sigma_global,
  weight_clip = 0.99,
  cor_min = 0.10,
  one_sided = "center_0.5",
  perm_methods = c("moran_mc", "geary_mc"),
  n_perm_min = 100,
  n_perm_max = 100,
  early_stop = TRUE,
  alpha = 0.05,
  seed_permutation = 124
)

spatial_summary <- perm_res$permutation_results$moran_geary_results$moran_geary_summary

if (is.null(spatial_summary) || nrow(spatial_summary) == 0L) {
  stop("Spatial permutation output was empty.")
}

# ------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------
write.csv(
  data.frame(
    spot_id = rownames(coords),
    coords,
    row.names = NULL
  ),
  file.path(out_dir, "minimal_example_coordinates.csv"),
  row.names = FALSE
)

write.csv(
  expr,
  file.path(out_dir, "minimal_example_logCPM.csv")
)

write.csv(
  z_sgsva,
  file.path(out_dir, "minimal_example_sgsva_zscores.csv")
)

write.csv(
  truth_summary,
  file.path(out_dir, "minimal_example_truth_recovery.csv"),
  row.names = FALSE
)

write.csv(
  spatial_summary,
  file.path(out_dir, "minimal_example_global_permutation_results.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    normalized_expression = expr,
    selected_pathways = sel_paths,
    sgsva_results = sgsva_res,
    truth_recovery = truth_summary,
    spatial_permutation = perm_res
  ),
  file.path(out_dir, "minimal_example_all_results.rds")
)

summary_lines <- c(
  "sGSVA fixed-input minimal reproducibility example",
  paste0("R version: ", R.version.string),
  paste0("Frozen input: ", input_file),
  paste0("Frozen input MD5: ", input_md5),
  paste0("Spatial locations: ", nrow(coords)),
  paste0("Input genes: ", nrow(counts)),
  paste0("Genes after filtering: ", nrow(expr)),
  paste0("Input KEGG pathways: ", length(sel_paths)),
  paste0("Retained pathways: ", nrow(z_sgsva)),
  paste0("ell (um): ", ell),
  paste0("sigma_f: ", sigma_f),
  paste0("Pathway composition: ", up, " up / ", down, " down"),
  paste0("Scoring bandwidth (um): ", sigma_global),
  "Spatial inference smoke test: global Moran's I and Geary's C",
  "Example permutations per statistic: 100",
  paste0(
    "Retained pathway names: ",
    paste(rownames(z_sgsva), collapse = ", ")
  )
)

writeLines(
  summary_lines,
  file.path(out_dir, "minimal_example_summary.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "minimal_example_session_info.txt")
)

message("")
message("Truth-recovery summary:")
print(truth_summary)

message("")
message("Minimal sGSVA example completed successfully.")
message("Outputs written to: ", normalizePath(out_dir, mustWork = FALSE))
