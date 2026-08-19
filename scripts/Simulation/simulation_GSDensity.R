# Simulation benchmarking
suppressPackageStartupMessages ({
  library(Matrix) ; library(matrixStats) ; library(MASS)
  library(dplyr) ; library(tidyr) ; library(purrr) ; 
  library(readr); library(spdep) ; library(sf) ; 
  library(pROC); library(msigdbr) ; library(ggplot2); 
  library(Seurat); library(patchwork); library(gsdensity);
  library(here); library(fs)
})

set.seed(1)
source("R/simulation.R")
source("R/simulation_scoring.R")
source("R/evaluation_metrics.R")

project_dir <- "."
outdir <- file.path(project_dir, "results", "simulation")
repdir <- file.path(outdir, "replicates")
qcdir <- file.path(outdir, "reproducibility")
supp_dir <- file.path(outdir, "results", "supplementary_data")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(repdir, recursive = TRUE, showWarnings = FALSE)
dir.create(qcdir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

sim_res <- readRDS(file.path(outdir, "rep_simulation.rds"))

# gene sets ( MSigDB KEGG legacy )
msig <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
gene_sets <- split(msig$gene_symbol, msig$gs_name)

coords <- make_hex_grid(nx=70, ny=70, pitch=50)
N <- nrow(coords)

# build spatial weights ( Delaunay triangulation )
nb <- tri2nb(coords)
listw <- nb2listw(nb, style = "W")
islands <- which(card(nb)==0)
length(islands)

# representative pathways
sel_paths <- c("KEGG_CHEMOKINE_SIGNALING_PATHWAY",
               "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
               "KEGG_FATTY_ACID_METABOLISM")

# simulaion scenarios
sizes <- tibble::tibble(up_size = c(25, 25),
                        down_size = c(0, 25))
grid <- tidyr::crossing(
  ell = c(50, 200, 400),
  sigma_f = c(0.5, 1.0, 1.5),
  sizes
)
replicates <- 20

# ------------------------------------------------------------
# Shared utilities
# ------------------------------------------------------------
prepare_expr <- function(counts, floor_frac = 0.05) {
  lib <- Matrix::colSums(counts) + 1
  lcp <- log1p(t(t(counts) / lib * 1e6))
  expr <- as.matrix(lcp)
  keep <- rowSums(expr > 0) >= floor(floor_frac * ncol(expr))
  expr[keep, , drop = FALSE]
}

eval_metrics_safe <- function(score, method, scn, sim, rt, q = 0.1,
                              expected_paths = names(sim$fields)) {
  score <- as.matrix(score)
  
  out <- data.frame(
    pathway = expected_paths, methods = method,
    AUROC = NA_real_, AUROC_fixed = NA_real_,
    spearman = NA_real_, spearman_signed = NA_real_,
    retained = expected_paths %in% rownames(score),
    times = rt, ell = scn$ell,
    sigma_f = scn$sigma_f, up_size = scn$up_size,
    down_size = scn$down_size,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(expected_paths)) {
    p <- expected_paths[i]
    if (!p %in% rownames(score)) next
    
    field <- sim$fields[[p]]
    truth <- top_mask(field, q = q)
    s <- as.numeric(score[p, ])
    
    auc_fixed <- score_auc_fixed(truth, s)
    auc_invariant <- score_auc_invariant(truth, s)
    out$AUROC[i] <- auc_invariant
    out$AUROC_fixed[i] <- auc_fixed
    
    rho <- suppressWarnings(cor(field, s, method = "spearman", use = "pairwise.complete.obs"))
    out$spearman_signed[i] <- rho
    out$spearman[i] <- abs(rho)
  }
  
  out
}

# ------------------------------------------------------------
# GSDensity wrapper
# ------------------------------------------------------------

# adjust original function to solve the compatibility issue
compute_mca_v5 <- function(
    object,
    dims.use = 1:10,
    genes.use = rownames(object),
    assay = Seurat::DefaultAssay(object),
    scale.factor = 1e6,
    verbose = FALSE
) {
  # Preserve the original compute.mca() logic:
  # 1. intersect genes
  # 2. NormalizeData
  # 3. RunMCA
  # 4. combine gene and spot embeddings
  
  Seurat::DefaultAssay(object) <- assay
  
  genes.use.1 <- intersect(genes.use, rownames(object))
  
  # Same normalization used by the original compute.mca()
  object <- Seurat::NormalizeData(object = object, assay = assay,
                                  normalization.method = "LogNormalize", 
                                  scale.factor = scale.factor)
  
  # Seurat v5-compatible replacement for:
  # GetAssayData(object, slot = "data")
  expression_data <- SeuratObject::LayerData(
    object = object, assay = assay, layer = "data")
  
  # Subset before converting to a dense matrix
  expression_data <- expression_data[genes.use.1,,drop = FALSE]
  
  # CelliD's matrix method is compatible because it does not access
  # a Seurat slot/layer internally.
  mca_result <- CelliD::RunMCA(
    X = as.matrix(expression_data),
    nmcs = max(dims.use),
    reduction.name = "MCA"
  )
  
  genes.use.update <- intersect(genes.use.1,
                                rownames(mca_result$featuresCoordinates))
  
  # Same returned object as gsdensity::compute.mca():
  # gene coordinates followed by spot coordinates
  coembed <- rbind(
    mca_result$featuresCoordinates[genes.use.update, dims.use, drop = FALSE],
    mca_result$cellsCoordinates[, dims.use, drop = FALSE]
  )
  
  return(coembed)
}

run_gsdensity_scores <- function(counts,
                                 sets,
                                 n_mca_dims = NULL,
                                 nn_use = 300,
                                 restart = 0.75,
                                 scale_factor = 1e6,
                                 verbose = FALSE) {
  
  seu <- Seurat::CreateSeuratObject(
    counts = as(counts, "dgCMatrix"),
    assay = "RNA",
    project = "simulation_gsdensity"
  )
  
  dims_use <- if (is.null(n_mca_dims)) 1:10 else seq_len(n_mca_dims)
  if (verbose) message("Computing GSDensity MCA co-embedding...")
  ce <- compute_mca_v5(
    object = seu, dims.use = dims_use, assay = "RNA", scale.factor = scale_factor
  )
  
  genes_in_mca <- intersect(rownames(ce), rownames(seu))
  sets_use <- lapply(
    sets,
    function(x) intersect(unique(x), genes_in_mca)
  )
  
  # match with the default gene.set.cutoff
  sets_use <- sets_use[lengths(sets_use) >= 3]
  
  if (length(sets_use) == 0) {
    return(matrix(
      NA_real_,
      nrow = 0,
      ncol = ncol(counts),
      dimnames = list(character(0), colnames(counts))
    ))
  }
  
  cells <- colnames(seu)
  
  if (verbose) message("Building GSDensity nearest-neighbor graph...")
  el <- gsdensity::compute.nn.edges(coembed = ce, nn.use = nn_use)
  
  score_mat <- gsdensity::run.rwr.list(
    el = el, gene_set_list = sets_use, cells = cells, restart = restart
  )
  
  gsdensity_score <- t(as.matrix(score_mat))
  z_normalization(gsdensity_score)
}


# ------------------------------------------------------------
# Spatial-method supplemental benchmark
# ------------------------------------------------------------
run_gsdensity_benchmark <- function(sim_res, gene_sets, sel_paths, replicates = replicates) {
  all_res <- list()
  
  for (r in seq_len(replicates)) {
    message("Replicate ", r)
    sim_rep <- sim_res[[paste0("replicate", r)]]
    rep_res <- list()
    
    for (g in seq_along(sim_rep)) {
      message("  Scenario ", g)
      
      scn <- sim_rep[[g]]$scn
      sim <- sim_rep[[g]]$sim
      sets_use <- gene_sets[sel_paths]
      
      scenario_res <- NULL
      t0 <- proc.time()[3]
      gsden_scores <- tryCatch(
        run_gsdensity_scores(
          counts = sim$counts,
          sets = sets_use,
          nn_use = 300,
          restart = 0.75
        ),
        error = function(e) {
          warning("GSDensity failed at replicate ", r, ", scenario ", g, ": ",
                  conditionMessage(e))
          matrix(
            NA_real_,
            nrow = 0,
            ncol = ncol(sim$counts),
            dimnames = list(character(0), colnames(sim$counts))
          )
        }
      )
      rt <- proc.time()[3] - t0
      
      scenario_res <- eval_metrics_safe(
        score = gsden_scores, method = "GSDensity", scn = scn,
        sim = sim, rt = rt, expected_paths = names(sim$fields)
      )
      
      rep_res[[g]] <- scenario_res
    }
    
    rep_df <- bind_rows(rep_res) %>% mutate(replicate = r)
    
    write.csv(
      rep_df,
      file.path(repdir, paste0("replicate", r, "_gsdensity_evaluation.csv")),
      row.names = FALSE
    )
    
    all_res[[r]] <- rep_df
  }
  
  all_df <- bind_rows(all_res)
  
  write.csv(
    all_df,
    file.path(outdir, "gsdensity_evaluation_all_replicates.csv"),
    row.names = FALSE
  )
  all_df
}


# ---- run GSDensity analysis ----
gsdensity_res <- run_gsdensity_benchmark(
  sim_res = sim_res,
  gene_sets = gene_sets,
  sel_paths = sel_paths,
  replicates = 20
)
write.csv(gsdensity_res, file.path(outdir, "evaluate_metrics_gsdensity.csv"))

# integrate GSDensity with existing methods
region_raw <- read.csv(file.path(outdir, "evaluate_metrics_all_replicates_data_driven.csv"))
region_raw$X <- NULL

region_with_gsdensity <- bind_rows(region_raw, gsdensity_res)
region_with_gsdensity <- region_with_gsdensity %>%
  mutate(
    methods = recode(
      methods,
      "GSVA_smoothed" = "GSVA-smoothed",
      "gsdensity" = "GSDensity"
    )
  )
write.csv(region_with_gsdensity, file.path(outdir, "evaluate_metrics_all_replicates_data_driven_plus_gsdensity.csv"))

# check availability and failure
gsdensity_availability <- gsdensity_res %>%
  mutate(
    evaluable = !is.na(AUROC),
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0 ~ "One-sided (25/0)"
    )
  ) %>%
  group_by(set_type, ell, sigma_f) %>%
  summarise(
    n_total = n(),
    n_evaluable = sum(evaluable),
    evaluable_rate = mean(evaluable),
    .groups = "drop"
  )
write.csv(gsdensity_availability, file.path(outdir, "gsdensity_availability.csv"))


# ------------------------------------------------------------
# Full-query simulation summary
# ------------------------------------------------------------
# Purpose:
# This summarizes all intended pathway-replicate-scenario queries,
# not only the cases retained by data-driven sGSVA.
#
# The retained-matched table remains the primary paired accuracy benchmark.
# This full-query table is for transparency about evaluability/coverage.

region_full_query <- region_with_gsdensity %>%
  mutate(
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0  ~ "One-sided (25/0)",
      TRUE ~ "Other"
    ),
    evaluable = retained == TRUE & !is.na(AUROC)
  )

# Overall full-query summary by method
full_query_summary <- region_full_query %>%
  group_by(methods) %>%
  summarise(
    n_queried = n(),
    n_evaluable = sum(evaluable, na.rm = TRUE),
    evaluable_rate = n_evaluable / n_queried,
    
    mean_auroc = mean(AUROC, na.rm = TRUE),
    sd_auroc = sd(AUROC, na.rm = TRUE),
    mean_abs_spearman = mean(spearman, na.rm = TRUE),
    sd_abs_spearman = sd(spearman, na.rm = TRUE),
    mean_signed_spearman = mean(spearman_signed, na.rm = TRUE),
    sd_signed_spearman = sd(spearman_signed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_auroc))

write.csv(
  full_query_summary,
  file.path(outdir, "supplement_simulation_full_query_simulation_summary.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# Retained-matched simulation summary
# ------------------------------------------------------------
# Identify sGSVA-retained pathway-replicate-scenario keys from main benchmark
retained_keys <- region_with_gsdensity %>%
  filter(methods == "sGSVA", retained == TRUE, !is.na(AUROC)) %>%
  distinct(pathway, replicate, ell, sigma_f, up_size, down_size)

region_retained_matched <- region_with_gsdensity %>%
  semi_join(
    retained_keys,
    by = c("pathway", "replicate", "ell", "sigma_f", "up_size", "down_size")
  ) %>%
  mutate(
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0 ~ "One-sided (25/0)",
      TRUE ~ "Other"
    )
  )


# ------------------------------------------------------------
# Retained vs non-retained sGSVA cases
# ------------------------------------------------------------
key_cols <- c(
  "pathway", "replicate", "ell", "sigma_f", "up_size", "down_size"
)

sgsva_status <- region_full_query %>%
  filter(methods == "sGSVA") %>%
  distinct(
    across(all_of(key_cols)),
    set_type,
    sgsva_retained = evaluable
  )

# Compare scenario composition
retained_composition <- sgsva_status %>%
  dplyr::count(sgsva_retained, set_type, sigma_f, ell) %>%
  dplyr::group_by(sgsva_retained) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

write_csv(
  retained_composition,
  file.path(outdir, "sgsva_retained_vs_nonretained_composition.csv")
)

# Retention rate by all main simulation scenario factors
sgsva_retention_by_scenario <- sgsva_status %>%
  group_by(set_type, sigma_f, ell) %>%
  summarise(
    n_queried = n(),
    n_retained = sum(sgsva_retained, na.rm = TRUE),
    n_nonretained = sum(!sgsva_retained, na.rm = TRUE),
    retention_rate = n_retained / n_queried,
    .groups = "drop"
  ) %>%
  arrange(set_type, sigma_f, ell)

write_csv(
  sgsva_retention_by_scenario,
  file.path(outdir, "supplement_sgsva_retention_by_settype_signal_ell.csv")
)

# Use other methods as a rough difficulty check on the same cases
competitor_case_difficulty <- region_full_query %>%
  filter(methods != "sGSVA") %>%
  group_by(across(all_of(key_cols))) %>%
  summarise(
    competitor_mean_auroc = mean(AUROC, na.rm = TRUE),
    competitor_mean_abs_spearman = mean(spearman, na.rm = TRUE),
    .groups = "drop"
  )

retained_vs_nonretained_summary <- sgsva_status %>%
  left_join(competitor_case_difficulty, by = key_cols) %>%
  group_by(sgsva_retained) %>%
  summarise(
    n_cases = n(),
    prop_low_signal = mean(sigma_f == 0.5, na.rm = TRUE),
    prop_moderate_signal = mean(sigma_f == 1.0, na.rm = TRUE),
    prop_high_signal = mean(sigma_f == 1.5, na.rm = TRUE),
    prop_balanced = mean(set_type == "Balanced (25/25)", na.rm = TRUE),
    prop_onesided= mean(set_type == "One-sided (25/0)", na.rm = TRUE),
    mean_competitor_auroc = mean(competitor_mean_auroc, na.rm = TRUE),
    mean_competitor_abs_spearman = mean(competitor_mean_abs_spearman, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  retained_vs_nonretained_summary,
  file.path(outdir, "sgsva_retained_vs_nonretained_summary.csv")
)


# ------------------------------------------------------------
# GSDensity parameter sensitivity check
# ------------------------------------------------------------
# Purpose:
# Check whether GSDensity results are driven by one arbitrary parameter choice.
#
# Baseline setting used in the main benchmark:
# n_mca_dims = 10, nn_use = 300, restart = 0.75
#
# Sensitivity settings vary one component at a time.

gsdensity_param_grid <- tibble::tribble(
  ~setting_id,        ~n_mca_dims, ~nn_use, ~restart,
  "baseline",         10,          300,     0.75,
  "mca_dims_5",        5,          300,     0.75,
  "mca_dims_20",      20,          300,     0.75,
  "nn_100",           10,          100,     0.75,
  "nn_500",           10,          500,     0.75,
  "restart_0.50",     10,          300,     0.50,
  "restart_0.90",     10,          300,     0.90
)

run_gsdensity_parameter_check <- function(
    sim_res,
    gene_sets,
    sel_paths,
    param_grid,
    replicates_use = 1:5
) {
  all_setting_res <- list()
  
  for (pp in seq_len(nrow(param_grid))) {
    
    setting <- param_grid[pp, ]
    
    message(
      "Running GSDensity sensitivity setting: ",
      setting$setting_id,
      " | dims=", setting$n_mca_dims,
      " | nn=", setting$nn_use,
      " | restart=", setting$restart
    )
    
    setting_res <- list()
    
    for (r in replicates_use) {
      
      message("  Replicate ", r)
      
      sim_rep <- sim_res[[paste0("replicate", r)]]
      rep_res <- list()
      
      for (g in seq_along(sim_rep)) {
        
        message("    Scenario ", g)
        
        scn <- sim_rep[[g]]$scn
        sim <- sim_rep[[g]]$sim
        sets_use <- gene_sets[sel_paths]
        
        t0 <- proc.time()[3]
        gsden_scores <- tryCatch(
          run_gsdensity_scores(
            counts = sim$counts,
            sets = sets_use,
            n_mca_dims = setting$n_mca_dims,
            nn_use = setting$nn_use,
            restart = setting$restart
          ),
          error = function(e) {
            warning(
              "GSDensity sensitivity failed at setting ",
              setting$setting_id,
              ", replicate ", r,
              ", scenario ", g,
              ": ",
              conditionMessage(e)
            )
            
            matrix(
              NA_real_,
              nrow = 0,
              ncol = ncol(sim$counts),
              dimnames = list(character(0), colnames(sim$counts))
            )
          }
        )
        rt <- proc.time()[3] - t0
        
        scenario_res <- eval_metrics_safe(
          score = gsden_scores,
          method = "GSDensity",
          scn = scn,
          sim = sim,
          rt = rt,
          expected_paths = names(sim$fields)
        ) %>%
          mutate(
            setting_id = setting$setting_id,
            n_mca_dims = setting$n_mca_dims,
            nn_use = setting$nn_use,
            restart = setting$restart,
            replicate = r,
            scenario_id = g,
            set_type = case_when(
              up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
              up_size == 25 & down_size == 0  ~ "One-sided (25/0)",
              TRUE ~ "Other"
            ),
            evaluable = retained == TRUE & !is.na(AUROC)
          )
        
        rep_res[[g]] <- scenario_res
      }
      
      setting_res[[as.character(r)]] <- bind_rows(rep_res)
    }
    
    all_setting_res[[setting$setting_id]] <- bind_rows(setting_res)
  }
  
  bind_rows(all_setting_res)
}

# Run small sensitivity check
gsdensity_param_res <- run_gsdensity_parameter_check(
  sim_res = sim_res,
  gene_sets = gene_sets,
  sel_paths = sel_paths,
  param_grid = gsdensity_param_grid,
  replicates_use = 1:5
)

write.csv(
  gsdensity_param_res,
  file.path(outdir, "gsdensity_parameter_sensitivity_pathway_level.csv"),
  row.names = FALSE
)

# Plot GSDensity parameter sensitivity
p_gsdensity_param_auc <- ggplot(
  gsdensity_param_res, aes(x = setting_id, y = AUROC)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = NULL,
    y = "AUROC",
    title = "GSDensity parameter sensitivity: AUROC"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  filename = file.path(outdir, "gsdensity_parameter_sensitivity_auroc.pdf"),
  plot = p_gsdensity_param_auc,
  width = 7,
  height = 4
)

p_gsdensity_param_spearman <- ggplot(gsdensity_param_res, aes(x = setting_id, y = spearman)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = NULL,
    y = "Absolute Spearman correlation",
    title = "GSDensity parameter sensitivity: map-to-truth recovery"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  filename = file.path(outdir, "gsdensity_parameter_sensitivity_spearman.pdf"),
  plot = p_gsdensity_param_spearman,
  width = 7,
  height = 4
)



# ------------------------------------------------------------
# Reproducibility metadata
# ------------------------------------------------------------

simulation_design <- tibble::tibble(
  item = c(
    "n_replicates",
    "n_spots",
    "grid_nx",
    "grid_ny",
    "pitch",
    "ell_values",
    "sigma_f_values",
    "pathway_types",
    "selected_pathways",
    "smoothing_sigma",
    "primary_metric",
    "map_recovery_metric",
    "retained_matched_definition"
  ),
  value = c(
    as.character(replicates),
    as.character(nrow(coords)),
    "70",
    "70",
    "50",
    paste(c(50, 200, 400), collapse = ", "),
    paste(c(0.5, 1.0, 1.5), collapse = ", "),
    "Balanced 25/25 and one-sided 25/0",
    paste(sel_paths, collapse = "; "),
    "floor(ell / sqrt(2))",
    "AUROC for top 10% active-region recovery",
    "absolute Spearman correlation with latent spatial field",
    "all methods restricted to pathway-replicate-scenario cases retained by data-driven sGSVA"
  )
)

readr::write_csv(
  simulation_design,
  file.path(qc_dir, "simulation_design_metadata.csv")
)

package_versions <- tibble::tibble(
  package = c(
    "Matrix", "matrixStats", "MASS", "dplyr", "tidyr", "purrr",
    "readr", "GSVA", "AUCell", "UCell", "spdep", "sf",
    "pROC", "msigdbr", "ggplot2", "Seurat", "SeuratObject",
    "gsdensity", "CelliD"
  ),
  version = purrr::map_chr(
    package,
    ~ if (requireNamespace(.x, quietly = TRUE)) {
      as.character(utils::packageVersion(.x))
    } else {
      NA_character_
    }
  )
)

readr::write_csv(
  package_versions,
  file.path(qc_dir, "simulation_package_versions.csv")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(qc_dir, "simulation_sessionInfo.txt")
)


library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "SData1_case_level_metrics")
writeData(wb, "SData1_case_level_metrics", region_with_gsdensity)
addWorksheet(wb, "SData1_gsdensity_sensitivity")
writeData(wb, "SData1_gsdensity_sensitivity", gsdensity_param_res)
saveWorkbook(wb, file.path(supp_dir, "Supplementary Data 1.xlsx"), overwrite = TRUE)







