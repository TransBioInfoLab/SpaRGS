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
rep1 <- sim_res[[paste0("replicate",1)]] # use replicate1 as representation
saveRDS(rep1, file.path(results_dir, "simulation_replicate1.rds"))



# ---------------------------
# benchmarking
# ---------------------------
# AUROC/Delong test results across replicates
res_auc <- NULL
res_delong <- NULL
for (r in 1:replicates) {
  sim_rep <- sim_res[[paste0("replicate",r)]]
  for (g in seq_len(length(sim_rep))) {
    # get simulation results
    scn <- sim_rep[[g]]$scn
    sim <- sim_rep[[g]]$sim
    
    # compute sGSVA scores
    t0 <- proc.time()[3]
    sgsva_res <- run_sgsva(sim$counts, sim$coords, sim$per_path, 
                           sigma_global = floor(scn$ell/sqrt(2)), 
                           directionality = "oracle", one_sided = "center_0.5")
    sgsva_scores <- sgsva_res$z_score
    rt <- proc.time()[3] - t0
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
    gsva_scores <- run_gsva(expr, sim$per_path)
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = gsva_scores, method = "GSVA", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    #GSVA on smoothed
    t0 <- proc.time()[3]
    gsva_smoothed <- run_gsva(smoothed, sim$per_path)
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = gsva_smoothed, method = "GSVA_on_smoothed_expression", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    #AUCell
    t0 <- proc.time()[3]
    aucell_scores <- run_aucell(expr, sim$per_path)
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = aucell_scores, method = "AUCell", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    #UCell
    t0 <- proc.time()[3]
    ucell_scores <- run_ucell(expr, sim$per_path)
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = ucell_scores, method = "UCell", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    #ssGSEA
    t0 <- proc.time()[3]
    ssgsea_scores <- run_ssgsea(expr, sim$per_path)
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = ssgsea_scores, method = "ssGSEA", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    #AddModuleScore
    t0 <- proc.time()[3]
    addmod_scores <- run_addmod(sim$counts, sim$per_path) # input counts
    rt <- proc.time()[3] - t0
    temp <- eval_metrics(score = addmod_scores, method = "AddModuleScore", scn, sim, rt)
    res_auc <- rbind(res_auc, temp)
    
    # perform Delong tests per pathway
    scores <- list(sGSVA=sgsva_scores, GSVA=gsva_scores, GSVA_on_smoothed_expression=gsva_smoothed,
                   AUCell=aucell_scores, UCell=ucell_scores, ssGSEA=ssgsea_scores, AddModuleScore=addmod_scores)
    temp_de <- bind_rows(lapply(method_pairs, function(methods) {
      delong_per_pathway(sim$fields, scores, methods, q = 0.10, scn)
    }))
    res_delong <- rbind(res_delong, temp_de)
  }
  #write.csv(res_auc, file.path(results_dir, paste0("replicate",r,"_auroc.csv")))
  #write.csv(res_delong, file.path(results_dir, paste0("replicate",r,"_delong_test.csv")))
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
    sgsva_res <- run_sgsva(sim$counts, sim$coords, sim$per_path, 
                           sigma_global = floor(scn$ell/sqrt(2)), 
                           directionality = "oracle", one_sided = "center_0.5")
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
      
      # get hex-block results
      out_moran_hex <- perm_with_early_stop(stat_moran_hex, x, listw, hid)
      out_geary_hex <- perm_with_early_stop(stat_geary_hex, x, listw, hid)
      
      # convert results into data.frame
      df[[p]] <- data.frame(
        pathway = pathway,
        obs_moran = out_moran_mc$obs, z_moran = out_moran_mc$z, 
        p_moran = out_moran_mc$p, z_moran_hex = out_moran_hex$z, 
        p_moran_hex = out_moran_hex$p, 
        obs_geary = -out_geary_mc$obs, z_geary = -out_geary_mc$z, 
        p_geary = out_geary_mc$p, z_geary_hex = -out_geary_hex$z,
        p_geary_hex = out_geary_hex$p,
        block_size = scn$ell, ell = scn$ell, sigma_f = scn$sigma_f,
        up_size = scn$up_size, down_size = scn$down_size,
        sigma_global = floor(scn$ell/sqrt(2)),
        stringsAsFactors = FALSE
      )
    }
    temp_spatial <- do.call(rbind, df)
    rownames(temp_spatial) <- NULL
    temp_spatial$p_moran <- ifelse(temp_spatial$p_moran < 0.001, "< 0.001", temp_spatial$p_moran)
    temp_spatial$p_moran_hex <- ifelse(temp_spatial$p_moran_hex < 0.001, "< 0.001", temp_spatial$p_moran_hex)
    temp_spatial$p_geary <- ifelse(temp_spatial$p_geary < 0.001, "< 0.001", temp_spatial$p_geary)
    temp_spatial$p_geary_hex <- ifelse(temp_spatial$p_geary_hex < 0.001, "< 0.001", temp_spatial$p_geary_hex)
    res_spatial <- rbind(res_spatial, temp_spatial)
  }
  #write.csv(res_spatial, file.path(results_dir, paste0("replicate",r,"_moran_geary.csv")))
}

# get metrics result
region_raw <- res_auc
region_raw <- region_raw %>%
  dplyr::mutate(methods=ifelse(methods=="GSVA_on_smoothed_expression", 
                               "GSVA_smoothed", methods))
write.csv(region_raw, file.path(results_dir, "AUROC_spearman_runtime_all_replicates_oracle.csv"))

# AUROC summary results
region_sum <- region_raw %>%
  group_by(pathway, methods, ell, sigma_f, up_size, down_size) %>%
  summarise(
    n = n(),
    mean_auroc = mean(AUROC, na.rm=TRUE),
    se_auroc = sd(AUROC, na.rm=TRUE)/sqrt(n)
  )
write.csv(region_sum, file.path(results_dir, "auroc_aggregate_sum_oracle.csv"))

# spearman correlation result
cor_sum <- region_raw %>%
  group_by(pathway, methods, ell, sigma_f, up_size, down_size) %>%
  summarise(
    n = n(),
    mean_cor = mean(spearman, na.rm=TRUE),
    se_cor = sd(spearman, na.rm=TRUE)/sqrt(n)
  )
write.csv(cor_sum, file.path(results_dir, "spearman_aggregate_sum_oracle.csv"))

# get running time
time_sum <- region_raw %>%
  mutate(
    set_type = case_when(
      up_size==25 & down_size==25 ~ "Balanced (25/25)",
      up_size==25 & down_size==0 ~ "One-sided (25/0)"
    )
  ) %>%
  group_by(methods, ell, sigma_f, set_type) %>%
  summarise(
    n = n(),
    mean_times = mean(times, na.rm=TRUE),
    se_times = sd(times, na.rm=TRUE)/sqrt(n)
  ) %>%
  mutate(
    ci_low  = mean_times - 1.96 * se_times,
    ci_high = mean_times + 1.96 * se_times
  )
write.csv(time_sum, file.path(results_dir, "time_aggregate_sum_oracle.csv"))


# get delong test results
delong_raw <- res_delong
delong_raw <- delong_raw %>%
  mutate(method2 = ifelse(method2=="GSVA_on_smoothed_expression",
                          "GSVA_smoothed", method2))
write.csv(delong_raw, file.path(results_dir, "Delongtest_all_replicates_oracle.csv"))

# summary results across replicates
delong_sum <- delong_raw %>%
  mutate(
    set_type = case_when(
      up_size == 25 & down_size == 25 ~ "Balanced (25/25)",
      up_size == 25 & down_size == 0  ~ "One-sided (25/0)"
    )
  ) %>%
  group_by(pathway, method1, method2, set_type, ell, sigma_f) %>%
  summarise(
    n_rep      = n(),
    auc1_mean  = mean(auc1, na.rm = TRUE),
    auc2_mean  = mean(auc2, na.rm = TRUE),
    auc_diff_mean = mean(auc_diff, na.rm = TRUE),
    auc_diff_se   = sd(auc_diff, na.rm = TRUE) / sqrt(n_rep),
    p_value     = stouffer_two_sided(p_delong, auc_diff),
    .groups    = "drop"
  )
write.csv(delong_sum, file.path(results_dir, "delongtest_aggregate_sum_oracle.csv"))


# get Moran's I
moran_raw <- res_spatial
write.csv(moran_raw, file.path(results_dir, "Moran_geary_all_replicates_oracle.csv"))

moran_sum <- moran_raw %>%
  mutate(
    set_type = case_when(
      up_size==25 & down_size==25 ~ "Balanced (25/25)",
      up_size==25 & down_size==0 ~ "One-sided (25/0)"
    )
  ) %>%
  group_by(pathway, block_size, ell, sigma_f, set_type, sigma_global) %>%
  summarise(
    obs_moran = mean(obs_moran, na.rm=TRUE),
    z_moran = mean(z_moran, na.rm=TRUE),
    z_moran_hex = mean(z_moran_hex, na.rm=TRUE),
    obs_geary = mean(obs_geary, na.rm=TRUE),
    z_geary = mean(z_geary, na.rm=TRUE),
    z_geary_hex = mean(z_geary_hex, na.rm=TRUE),
  )
write.csv(moran_sum, file.path(results_dir, "moran_geary_aggregate_sum_oracle.csv"))
