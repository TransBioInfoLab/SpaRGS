### Evaluate EnrichMap in breast Visium data

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
qc_dir <- file.path(results_dir, "reproducibility")
supp_dir <- file.path(project_dir, "results", "supplementary_data")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_dir, "R/sgsva_score.R"))
source(file.path(project_dir, "R/spatial_permutation.R"))


################################################################################
# evaluate breast visium EnrichMap output
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

enrichmap_coverage <- read.csv(file.path(pilot_dir, "breast_visium_enrichmap_gene_set_coverage.csv"))

# read other pathway comparators
sgsva_scores <- read.csv(file.path(results_dir, "sgsva_scores.csv"),
                         row.names = 1,
                         check.names = FALSE)
gsva_scores <- read.csv(file.path(results_dir, "gsva_scores.csv"),
                        row.names = 1,
                        check.names = FALSE)
gsva_smoothed_scores <- read.csv(file.path(results_dir, "gsva_smoothed_scores.csv"),
                                 row.names = 1,
                                 check.names = FALSE)
gsdensity_scores <- read.csv(file.path(results_dir, "gsdensity_scores.csv"),
                             row.names = 1,
                             check.names = FALSE)

# get gene-set coverage table
exploratory_score_list <- list(
  sGSVA = as.matrix(sgsva_scores),
  GSVA = as.matrix(gsva_scores),
  `GSVA-smoothed` = as.matrix(gsva_smoothed_scores),
  GSDensity = as.matrix(gsdensity_scores),
  `EnrichMap-default` = as.matrix(em_default),
  `EnrichMap-smooth-only` = as.matrix(em_smooth)
)

# Optional pathway-name harmonization.
normalize_pathway_names <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  sub("^KEGG_", "", x)
}

exploratory_score_list <- lapply(exploratory_score_list, function(mat) {
  rownames(mat) <- normalize_pathway_names(rownames(mat))
  mat
})

exploratory_common_spots <- Reduce(
  intersect,
  c(lapply(exploratory_score_list, colnames), list(rownames(coords_um)))
)

exploratory_shared_pathways <- Reduce(
  intersect,
  lapply(exploratory_score_list, rownames)
)

enrichmap_exploratory_coverage <- purrr::imap_dfr(
  exploratory_score_list,
  function(mat, method_name) {
    tibble(
      method = method_name,
      pathways_available = nrow(mat),
      spots_available = ncol(mat),
      shared_pathways_for_exploratory_comparison = length(exploratory_shared_pathways),
      common_spots_for_exploratory_comparison = length(exploratory_common_spots)
    )
  }
)

readr::write_csv(
  enrichmap_exploratory_coverage,
  file.path(pilot_dir, "supplement_enrichmap_exploratory_coverage.csv")
)


################################################################################
# Method-wise Moran's I and Geary's C comparisons
# Collect score matrices
score_list <- list(
  sGSVA = as.matrix(sgsva_scores),
  GSVA = as.matrix(gsva_scores),
  GSVA_smoothed = as.matrix(gsva_smoothed_scores),
  GSDensity = as.matrix(gsdensity_scores),
  EnrichMap_default = as.matrix(em_default),
  EnrichMap_smooth = as.matrix(em_smooth)
)

score_list <- lapply(score_list, function(mat) {
  rownames(mat) <- normalize_pathway_names(rownames(mat))
  mat
})

# Restrict all methods to the same spots and shared pathways
common_spots <- Reduce(intersect,
                       c(lapply(score_list, colnames),list(rownames(coords_um))))

shared_pathways <- Reduce(intersect, lapply(score_list, rownames))

coords_common <- as.data.frame(coords_um[common_spots, , drop = FALSE])
coords_common <- coords_common[, 1:2, drop = FALSE]

score_list <- lapply(score_list, function(mat) {
  mat[shared_pathways, common_spots, drop = FALSE]
})

message("Common spots: ", length(common_spots))
message("Shared pathways: ", length(shared_pathways))

# Build one common Delaunay graph for every method
nb <- spdep::tri2nb(as.matrix(coords_common))
lw <- spdep::nb2listw(
  nb,
  style = "W",
  zero.policy = TRUE
)

# Compute observed Moran's I and Geary's C
spatial_stats <- imap_dfr(score_list, function(mat, method_name) {
  map_dfr(shared_pathways, function(pathway_name) {
    x <- as.numeric(mat[pathway_name, ])
    
    if (all(!is.finite(x)) || sd(x, na.rm = TRUE) == 0) {
      return(tibble(
        method = method_name,
        pathway = pathway_name,
        moran_i = NA_real_,
        geary_c = NA_real_
      ))
    }
    
    moran_value <- spdep::moran.test(x, lw, zero.policy = TRUE)$estimate[["Moran I statistic"]]
    geary_value <- spdep::geary.test(x, lw, zero.policy = TRUE)$estimate[["Geary C statistic"]]
    
    tibble(
      method = method_name,
      pathway = pathway_name,
      moran_i = moran_value,
      geary_c = geary_value
    )
  })
})

# save enrichmap pilot
readr::write_csv(
  spatial_stats,
  file.path(pilot_dir, "breast_visium_enrichmap_pilot_pathway_spatial_stats.csv")
)



################################################################################
# get cell type proportions
cell_type_proportions <- read.csv(file.path(results_dir, "celltype_anchoring", "rctd_cell_type_proportions.csv"))
rownames(cell_type_proportions) <- cell_type_proportions$barcode
cell_type_proportions$barcode <- NULL

# evaluate enrichmap anchoring
anchor_pairs_eval <- tribble(
  ~pathway,                                  ~cell_type_key,       ~anchor_label,
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",  "B_cell",             "BCR signaling vs B cells",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",  "T_cell",             "TCR signaling vs T cells",
  "KEGG_GLUTATHIONE_METABOLISM",              "Cancer_epithelial",  "Glutathione metabolism vs cancer epithelial",
  "KEGG_FATTY_ACID_METABOLISM",               "Cancer_epithelial",  "Fatty-acid metabolism vs cancer epithelial"
) %>%
  mutate(pathway_norm = normalize_pathway_names(pathway))

resolve_cell_type_simple <- function(key, available_names) {
  patterns <- switch(
    key,
    B_cell = c("^B[._ -]?CELLS?$", "B[._ -]?CELL"),
    T_cell = c("^T[._ -]?CELLS?$", "T[._ -]?CELL"),
    Cancer_epithelial = c(
      "CANCER.*EPITHELIAL",
      "EPITHELIAL.*CANCER",
      "MALIGNANT",
      "TUMOU?R.*EPITHELIAL"
    ),
    stop("Unknown key: ", key)
  )
  
  hits <- unique(unlist(lapply(patterns, function(p) {
    grep(p, available_names, ignore.case = TRUE, value = TRUE)
  })))
  
  if (length(hits) != 1L) {
    stop(
      "Could not uniquely resolve ", key, ". Hits: ",
      paste(hits, collapse = ", "),
      "\nAvailable: ",
      paste(available_names, collapse = ", ")
    )
  }
  
  hits
}

resolved_cell_types <- c(
  B_cell = resolve_cell_type_simple("B_cell", colnames(cell_type_proportions)),
  T_cell = resolve_cell_type_simple("T_cell", colnames(cell_type_proportions)),
  Cancer_epithelial = resolve_cell_type_simple("Cancer_epithelial", colnames(cell_type_proportions))
)

evaluate_anchoring <- function(score_mat, method_name) {
  score_mat <- as.matrix(score_mat)
  rownames(score_mat) <- normalize_pathway_names(rownames(score_mat))
  
  common_spots <- intersect(colnames(score_mat), rownames(cell_type_proportions))
  
  dplyr::bind_rows(lapply(seq_len(nrow(anchor_pairs_eval)), function(i) {
    p <- anchor_pairs_eval$pathway_norm[i]
    ct <- resolved_cell_types[anchor_pairs_eval$cell_type_key[i]]
    
    if (!p %in% rownames(score_mat)) {
      return(tibble(
        method = method_name,
        anchor_label = anchor_pairs_eval$anchor_label[i],
        pathway = p,
        cell_type = ct,
        n = NA_integer_,
        signed_rho = NA_real_,
        abs_rho = NA_real_,
        status = "pathway_not_available"
      ))
    }
    
    x <- as.numeric(score_mat[p, common_spots])
    y <- as.numeric(cell_type_proportions[common_spots, ct])
    
    keep <- is.finite(x) & is.finite(y)
    rho <- if (sum(keep) >= 10 && sd(x[keep]) > 0 && sd(y[keep]) > 0) {
      suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
    } else {
      NA_real_
    }
    
    tibble(
      method = method_name,
      anchor_label = anchor_pairs_eval$anchor_label[i],
      pathway = p,
      cell_type = ct,
      n = sum(keep),
      signed_rho = rho,
      abs_rho = abs(rho),
      status = "evaluated"
    )
  }))
}

em_anchor <- bind_rows(
  evaluate_anchoring(em_default, "EnrichMap-default"),
  evaluate_anchoring(em_smooth, "EnrichMap-smooth")
)

write_csv(
  em_anchor,
  file.path(pilot_dir, "breast_visium_enrichmap_anchoring.csv")
)


# combine six methods anchoring
main_anchor <- read_csv(
  file.path(results_dir, "celltype_anchoring", "anchoring_prespecified_pairs.csv"),
  show_col_types = FALSE
)

# Standardize column names if needed
main_anchor_clean <- main_anchor %>%
  transmute(
    method = method,
    anchor_label = paste(pathway, cell_type, sep = " vs "),
    pathway = pathway,
    cell_type = cell_type,
    n = n,
    signed_rho = rho,
    abs_rho = abs(rho),
    source = "primary_four_method_benchmark"
  )

em_anchor_clean <- em_anchor %>%
  transmute(
    method = method,
    anchor_label = anchor_label,
    pathway = pathway,
    cell_type = cell_type,
    n = n,
    signed_rho = signed_rho,
    abs_rho = abs_rho,
    source = "exploratory_enrichmap_pilot"
  )

six_method_anchor <- bind_rows(
  main_anchor_clean,
  em_anchor_clean
) %>%
  distinct(method, anchor_label, .keep_all = TRUE)

# Standardize six-method anchoring labels
canonical_anchor_map <- tibble::tribble(
  ~pathway_norm,                         ~anchor_id,        ~pathway_display,          ~cell_type_display,      ~anchor_label_display,                          ~anchor_order,
  "B_CELL_RECEPTOR_SIGNALING_PATHWAY",   "BCR_Bcell",       "BCR signaling",           "B cells",               "BCR signaling vs B cells",                         1,
  "T_CELL_RECEPTOR_SIGNALING_PATHWAY",   "TCR_Tcell",       "TCR signaling",           "T cells",               "TCR signaling vs T cells",                         2,
  "GLUTATHIONE_METABOLISM",              "GSH_epithelial",  "Glutathione metabolism",  "Cancer epithelial",     "Glutathione metabolism vs cancer epithelial",      3,
  "FATTY_ACID_METABOLISM",               "FA_epithelial",   "Fatty-acid metabolism",   "Cancer epithelial",     "Fatty-acid metabolism vs cancer epithelial",       4
)

method_order <- c(
  "sGSVA",
  "GSVA",
  "GSVA-smoothed",
  "GSDensity",
  "EnrichMap-default",
  "EnrichMap-smooth"
)

six_method_anchor_clean <- six_method_anchor %>%
  dplyr::mutate(
    pathway_norm = normalize_pathway_names(pathway)
  ) %>%
  left_join(canonical_anchor_map, by = "pathway_norm") %>%
  dplyr::mutate(
    method = factor(method, levels = method_order),
    anchor_label_display = factor(anchor_label_display, levels = canonical_anchor_map$anchor_label_display),
    signed_rho = as.numeric(signed_rho),
    abs_rho = abs(signed_rho)
  ) %>%
  arrange(method, anchor_order) %>%
  dplyr::select(method, anchor_label_display,
                pathway_display, cell_type_display, n,
                signed_rho, abs_rho, source)

write_csv(
  six_method_anchor_clean,
  file.path(results_dir, "enrichmap_pilot", "supplement_exploratory_six_method_anchoring.csv")
)


# supplementary data
library(openxlsx)
xlsx_file <- file.path(supp_dir, "Supplementary Data 3.xlsx")
wb <- loadWorkbook(xlsx_file)
addWorksheet(wb, "SData3_spatial_statistics")
writeData(wb, "SData3_spatial_statistics", spatial_stats)
addWorksheet(wb, "SData3_celltype_concordance")
writeData(wb, "SData3_celltype_concordance", cell_type_proportions)
addWorksheet(wb, "SData3_enrichmap_coverage")
writeData(wb, "SData3_enrichmap_coverage", enrichmap_coverage)
saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)





