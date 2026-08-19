# Simulation benchmarking
suppressPackageStartupMessages ({
  library(Matrix) ; library(matrixStats) ; library(MASS)
  library(dplyr) ; library(tidyr) ; library(purrr) ; 
  library(readr); library(GSVA) ; library(AUCell) ; 
  library(UCell); library(spdep) ; library(sf) ; 
  library(pROC); library(msigdbr) ; library(ggplot2); 
  library(Seurat); library(patchwork)
})

set.seed(1)

source("R/spatial_permutation.R")
source("R/simulation.R")
source("R/simulation_scoring.R")
source("R/evaluation_metrics.R")

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
results_dir <- file.path(project_dir, "results", "simulation")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------
# generate simulation 
# ---------------------------
# ---- gene sets ( MSigDB KEGG legacy ) ----
msig <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
gene_sets <- split(msig$gene_symbol, msig$gs_name)

# pick a handful for demo ; in the real run , use all sets passing filters
sel_paths <- c("KEGG_CHEMOKINE_SIGNALING_PATHWAY",
               "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
               "KEGG_FATTY_ACID_METABOLISM")

coords <- make_hex_grid(nx=70, ny=70, pitch=50)
N <- nrow(coords)

# build spatial weights ( Delaunay triangulation )
nb <- tri2nb(coords)
listw <- nb2listw(nb, style = "W")
islands <- which(card(nb)==0)
length(islands)

# --- full grid scenario ---
sizes <- tibble::tibble(up_size = c(25, 25),
                        down_size = c(0, 25))
grid <- tidyr::crossing(
  ell = c(50, 200, 400),
  sigma_f = c(0.5, 1.0, 1.5),
  sizes
)
replicates <- 20
method_pairs <- list(c("sGSVA","GSVA"), c("sGSVA", "GSVA_smoothed"),
                     c("sGSVA", "AUCell"), c("sGSVA", "UCell"), 
                     c("sGSVA", "ssGSEA"), c("sGSVA", "AddModuleScore"))

# generate simulation for multiple scenarios across replicates
sim_res <- NULL
for (r in 1:replicates) {
  rep <- NULL
  for (g in seq_len(dim(grid)[1])) {
    # get simulation results
    scn <- grid[g,]
    sim <- simulate_counts(coords, gene_sets, sel_paths,
                           up=scn$up_size, down=scn$down_size, ell=scn$ell, sigma_f=scn$sigma_f,
                           mean_lib=8000, sd_lib=0.6)
    rep[[g]] <- list(scn=scn, sim=sim)
  }
  sim_res[[paste0("replicate",r)]] <- rep
}
saveRDS(sim_res, file.path(results_dir, "rep_simulation.rds"))



# ---------------------------
# benchmarking
# ---------------------------
# AUROC results across replicates
res_auc <- NULL
for (r in 1:replicates) {
  sim_rep <- sim_res[[paste0("replicate",r)]]
  for (g in seq_len(length(sim_rep))) {
    # get simulation results
    scn <- sim_rep[[g]]$scn
    sim <- sim_rep[[g]]$sim
    
    # compute sGSVA scores
    t0 <- proc.time()[3]
    sgsva_res <- run_sgsva(sim$counts, sim$coords, gene_sets[sel_paths], 
                           sigma_global = floor(scn$ell/sqrt(2)), 
                           directionality = "data-driven", one_sided = "center_0.5")
    rt <- proc.time()[3] - t0
    sgsva_scores <- sgsva_res$z_score
    temp <- eval_metrics(score = sgsva_scores, method = "sGSVA", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    
    # locCPM and filter the expression matrix
    counts <- sim$counts
    G <- nrow(counts); N <- ncol(counts)
    lib <- colSums(counts)+1
    lcp <- log1p(t(t(counts)/lib*1e6)) #logCPM
    expr <- as.matrix(lcp)
    keep <- rowSums(expr > 0) >= floor(0.05 * ncol(expr))
    expr <- expr[keep, , drop = FALSE]
    # get smoothed expression
    smoothed <- smooth_expr(expr, sim$coords, sigma=floor(scn$ell/sqrt(2)), mode="delaunay")
    smoothed <- as.matrix(smoothed)
    colnames(smoothed) <- colnames(expr)
    #GSVA
    t0 <- proc.time()[3]
    gsva_scores <- run_gsva(expr, gene_sets[sel_paths])
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = gsva_scores, method = "GSVA", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    #GSVA on smoothed
    t0 <- proc.time()[3]
    gsva_smoothed <- run_gsva(smoothed, gene_sets[sel_paths])
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = gsva_smoothed, method = "GSVA_smoothed", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    #AUCell
    t0 <- proc.time()[3]
    aucell_scores <- run_aucell(expr, gene_sets[sel_paths])
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = aucell_scores, method = "AUCell", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    #UCell
    t0 <- proc.time()[3]
    ucell_scores <- run_ucell(expr, gene_sets[sel_paths])
    rt <- proc.time()[3] - t0
    rownames(ucell_scores) <- sel_paths
    temp <- eval_metrics(score = ucell_scores, method = "UCell", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    #ssGSEA
    t0 <- proc.time()[3]
    ssgsea_scores <- run_ssgsea(expr, gene_sets[sel_paths])
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = ssgsea_scores, method = "ssGSEA", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    #AddModuleScore
    t0 <- proc.time()[3]
    addmod_scores <- run_addmod(sim$counts, gene_sets[sel_paths]) # input counts
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = addmod_scores, method = "AddModuleScore", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    res_auc$replicate <- r
  }
  #write.csv(res_auc, file.path(base_dir, paste0("replicate",r,"_auroc.csv")))
}


# Moran's I results across replicates
res_spatial <- NULL
for (r in 1:replicates) {
  sim_rep <- sim_res[[paste0("replicate",r)]]
  for (g in seq_len(length(sim_rep))) {
    # get simulation results
    scn <- sim_rep[[g]]$scn
    sim <- sim_rep[[g]]$sim
    sel_paths <- names(sim$per_path)
    
    # get sgsva scores
    sgsva_res <- run_sgsva(sim$counts, sim$coords, gene_sets[sel_paths], 
                           sigma_global = floor(scn$ell/sqrt(2)), 
                           directionality = "data-driven", one_sided = "center_0.5")
    sgsva_scores <- sgsva_res$z_score
    hid <- make_hex_ids(coords, cellsize = scn$ell)
    df <- vector("list", length(sel_paths))
    
    for (p in seq_along(sel_paths)) {
      pathway <- sel_paths[p]
      x <- sgsva_scores[pathway,]
      
      set.seed(123)
      # get moran.mc/geary.mc results
      out_moran_mc <- perm_with_early_stop(stat_moran, x, listw)
      out_geary_mc <- perm_with_early_stop(stat_geary, x, listw)
      
      # convert results into data.frame
      df[[p]] <- data.frame(
        pathway = pathway,
        obs_moran = out_moran_mc$obs, z_moran = out_moran_mc$z, 
        p_moran = out_moran_mc$p, 
        obs_geary = -out_geary_mc$obs, z_geary = -out_geary_mc$z, 
        p_geary = out_geary_mc$p, 
        block_size = scn$ell, ell = scn$ell, sigma_f = scn$sigma_f,
        up_size = scn$up_size, down_size = scn$down_size,
        sigma_global = floor(scn$ell/sqrt(2)),
        stringsAsFactors = FALSE
      )
    }
    temp_spatial <- do.call(rbind, df)
    rownames(temp_spatial) <- NULL
    temp_spatial$p_moran <- ifelse(temp_spatial$p_moran < 0.001, "< 0.001", temp_spatial$p_moran)
    temp_spatial$p_geary <- ifelse(temp_spatial$p_geary < 0.001, "< 0.001", temp_spatial$p_geary)
    res_spatial <- rbind(res_spatial, temp_spatial)
  }
  #write.csv(res_spatial, file.path(results_dir, paste0("replicate",r,"_moran_geary.csv")))
}



# --- organize summary results ---
# get metrics result
region_raw <- res_auc
write.csv(region_raw, file.path(results_dir, "evaluate_metrics_all_replicates_data_driven.csv"))

# signed vs absolute spearman orientation diagnostic
# Identify pathway-replicate-scenario pairs retained by data-driven sGSVA
retained_keys <- region_raw %>%
  filter(methods == "sGSVA", retained == TRUE, !is.na(AUROC)) %>%
  distinct(pathway, replicate, ell, sigma_f, up_size, down_size)

# Retained-matched dataset: all methods evaluated on sGSVA-retained pathway-replicate pairs
region_retained_matched <- region_raw %>%
  semi_join(
    retained_keys,
    by = c("pathway", "replicate", "ell", "sigma_f", "up_size", "down_size")
  ) %>%
  mutate(
    set_type = ifelse(down_size == 0, "One-sided (25/0)", "Balanced (25/25)")
  )

orientation_diag <- region_retained_matched %>%
  group_by(methods, set_type, ell, sigma_f) %>%
  summarise(
    n_evaluable = sum(!is.na(spearman_signed)),
    mean_signed_spearman = mean(spearman_signed, na.rm = TRUE),
    mean_abs_spearman = mean(spearman, na.rm = TRUE),
    negative_spearman_rate = mean(spearman_signed < 0, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(orientation_diag, 
          file.path(results_dir, "simulation_signed_spearman_orientation_diagnostic.csv"))


# get Moran's I/Geary's C
moran_raw <- res_spatial
write.csv(moran_raw, file.path(results_dir, "spatial_test_all_replicates_data_driven.csv"))

moran_sum <- moran_raw %>%
  mutate(
    set_type = case_when(
      up_size==25 & down_size==25 ~ "Balanced (25/25)",
      up_size==25 & down_size==0 ~ "One-sided (25/0)"
    )
  ) %>%
  group_by(pathway, ell, sigma_f, set_type, sigma_global) %>%
  summarise(
    n_retained = sum(retained),
    obs_moran = mean(obs_moran, na.rm=TRUE),
    z_moran = mean(z_moran, na.rm=TRUE),
    obs_geary = mean(obs_geary, na.rm=TRUE),
    z_geary = mean(z_geary, na.rm=TRUE)
  )
write.csv(moran_sum, file.path(results_dir, "spatial_test_aggregate_sum_data_driven.csv"))
