# Ablation analysis
suppressPackageStartupMessages ({
  library(Matrix) ; library(matrixStats) ; library(MASS)
  library(dplyr) ; library(tidyr) ; library(purrr) ; 
  library(readr); library(spdep) ; library(sf) ; 
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
results_dir <- file.path(project_dir, "results", "ablation")
sim_dir <- file.path(project_dir, "results", "simulation")
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

# build spatial weights (Delaunay triangulation)
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

# simulation for two set types with all KEGG pathways
sel_paths <- names(gene_sets)
sizes <- tibble::tibble(up_size = c(25, 25),
                        down_size = c(0, 25))
grid <- tidyr::crossing(sizes)
sim_set_type <- NULL

for (g in seq_len(dim(grid)[1])) {
  up <- grid$up_size[g]
  down <- grid$down_size[g]
  sim <- simulate_counts(coords, gene_sets, sel_paths, 
                         up = up, down = down) # default ell, sigma_f
  sim_set_type[[paste0("up",up,"_down",down)]] <- sim
}
#saveRDS(sim_set_type, file.path(results_dir, "sim_set_type.rds"))

# simulation for varied ell with all KEGG pathways
ell <- c(50, 200, 400)
sel_paths <- names(gene_sets)
sim_block_res <- NULL
for (l in ell) {
  sim <- simulate_counts(coords, gene_sets, sel_paths, ell=l)
  sim_block_res[[paste0("block",l)]] <- sim
}
#saveRDS(sim_block_res, file.path(results_dir, "sim_block_res.rds"))

# read all simulations
sim_set_type <- readRDS(file.path(results_dir, "sim_set_type.rds"))
sim_res <- readRDS(file.path(sim_dir, "rep_simulation.rds"))
sim_block_res <- readRDS(file.path(results_dir, "sim_block_res.rds"))


# ---------------------------
# ablation analysis
# ---------------------------
# correlation threshold for directionality
rho <- c(0, 0.05, 0.1, 0.2)
res_rho <- NULL
for (g in seq_len(dim(grid)[1])) {
  up <- grid$up_size[g]
  down <- grid$down_size[g]
  sim <- sim_set_type[[paste0("up",up,"_down",down)]]
  
  t0 <- proc.time()[3]
  sgsva_res <- run_sgsva(sim$counts, sim$coords, sim$per_path, 
                         directionality = "oracle", one_sided = "center_0.5")
  rt0 <- proc.time()[3] - t0
  sgsva_oracle <- sgsva_res$z_score
  
  for (r in rho) {
    t0 <- proc.time()[3]
    sgsva_res <- run_sgsva(sim$counts, sim$coords, gene_sets, one_sided = "center_0.5",
                           directionality = "data-driven", cor_min = r)
    rt <- proc.time()[3] - t0
    sgsva_scores <- sgsva_res$z_score
    final_sets <- sgsva_res$gene_sets
    
    sel_paths <- rownames(sgsva_scores)
    temp <- data.frame(
      pathway = sel_paths, rho = r,
      AUROC_data_driven = NA_real_, AUROC_oracle = NA_real_,
      auc_diff = NA_real_, p_delong = NA_real_,
      up_size = NA_real_, down_size = NA_real_,
      true_up = up, true_down = down,
      times_data_driven = rt, times_oracle = rt0,
      reduced_paths = length(gene_sets)-length(final_sets),
      stringsAsFactors = FALSE
    )
    for (i in seq_along(sel_paths)) {
      p <- sel_paths[i]
      field <- sim$fields[[p]]
      truth <- top_mask(field, q=0.1)
      r1 <- pROC::roc(truth, sgsva_scores[p,], quiet = TRUE)
      r2 <- pROC::roc(truth, sgsva_oracle[p,], quiet = TRUE)
      tt <- pROC::roc.test(r1, r2, method = "delong", paired = TRUE)
      
      temp$AUROC_data_driven[i] <- as.numeric(pROC::auc(r1))
      temp$AUROC_oracle[i] <- as.numeric(pROC::auc(r2))
      temp$auc_diff[i] <- as.numeric(pROC::auc(r1) - pROC::auc(r2))
      temp$p_delong[i] <- tt$p.value
      temp$up_size[i] <- length(final_sets[[p]]$U)
      temp$down_size[i] <- length(final_sets[[p]]$D)
    }
    res_rho <- rbind(res_rho, temp)
  }
}
write.csv(res_rho, file.path(results_dir, "abl_correlation_threshold.csv"))

# bandwidth sensitivity 
param <- c(1/sqrt(2), 1, 1.5) # (1/sqrt(2), 1, 1.5)*ell
res_sigma_oracle <- NULL
res_sigma_data <- NULL
for (r in 1:replicates) {
  sim_rep <- sim_res[[paste0("replicate",r)]]
  for (g in seq_len(length(sim_rep))) {
    scn <- sim_rep[[g]]$scn
    sim <- sim_rep[[g]]$sim
    for (p in param) {
      sigma_global <- floor(scn$ell*p)
      path_names <- names(sim$per_path)
      # get sgsva scores
      t0 <- proc.time()[3]
      sgsva_res_oracle <- run_sgsva(sim$counts, sim$coords, sim$per_path, sigma_global, 
                             directionality = "oracle", one_sided = "center_0.5")
      rt <- proc.time()[3] - t0
      sgsva_scores <- sgsva_res_oracle$z_score
      sets <- sgsva_res_oracle$gene_sets
      temp <- eval_metrics_abl(score = sgsva_scores, sets, method = "sGSVA", scn, sim, rt)
      temp$sigma_param <- p
      temp$sigma_global <- sigma_global
      res_sigma_oracle <- rbind(res_sigma_oracle, temp)
      
      # get sgsva scores
      t0 <- proc.time()[3]
      sgsva_res_data <- run_sgsva(sim$counts, sim$coords, gene_sets[path_names], sigma_global, 
                             directionality = "data-driven", one_sided = "center_0.5")
      rt <- proc.time()[3] - t0
      sgsva_scores <- sgsva_res_data$z_score
      sets <- sgsva_res_data $gene_sets
      temp <- eval_metrics_abl(score = sgsva_scores, sets, method = "sGSVA", scn, sim, rt)
      temp$sigma_param <- p
      temp$sigma_global <- sigma_global
      res_sigma_data <- rbind(res_sigma_data, temp)
    }
  }
}
write.csv(res_sigma_oracle, file.path(results_dir, "abl_bandwidths_all_replicates_oracle.csv"))
write.csv(res_sigma_data, file.path(outdir, "abl_bandwidths_all_replicates_data_driven.csv"))

bandwidth_sum_oracle <- res_sigma_oracle %>%
  filter(!is.na(AUROC)) %>%
  mutate(
    set_type = case_when(
      true_up==25 & true_down==25 ~ "Balanced (25/25)",
      true_up==25 & true_down==0 ~ "One-sided (25/0)"
    ),
    sigma_param = case_when(
      sigma_param==1.5 ~ "1.5*ell",
      sigma_param==1 ~ "1*ell",
      TRUE ~ "1/sqrt(2)*ell"
    )
  )  %>%
  group_by(pathway, ell, sigma_f, set_type, sigma_param, sigma_global) %>%
  summarise(
    n = n(),
    mean_auroc = mean(AUROC, na.rm=TRUE),
    se_auroc = sd(AUROC, na.rm=TRUE)/sqrt(n),
    mean_spearman = mean(spearman, na.rm=TRUE),
    se_spearman = sd(spearman, na.rm=TRUE)/sqrt(n),
    mean_times = mean(times, na.rm=TRUE),
    se_times = sd(times, na.rm=TRUE)/sqrt(n),
  )
bandwidth_sum_oracle$sigma_param <- factor(bandwidth_sum_oracle$sigma_param, levels = c("1/sqrt(2)*ell", "1*ell", "1.5*ell"))
write.csv(bandwidth_sum, file.path(results_dir, "abl_bandwidths_aggregate_sum_oracle.csv"))

bandwidth_sum_data <- res_sigma_data %>%
  filter(!is.na(AUROC)) %>%
  mutate(
    set_type = case_when(
      true_up==25 & true_down==25 ~ "Balanced (25/25)",
      true_up==25 & true_down==0 ~ "One-sided (25/0)"
    ),
    sigma_param = case_when(
      sigma_param==1.5 ~ "1.5*ell",
      sigma_param==1 ~ "1*ell",
      TRUE ~ "1/sqrt(2)*ell"
    )
  )  %>%
  group_by(pathway, ell, sigma_f, set_type, sigma_param, sigma_global) %>%
  summarise(
    n = n(),
    mean_auroc = mean(AUROC, na.rm=TRUE),
    se_auroc = sd(AUROC, na.rm=TRUE)/sqrt(n),
    mean_spearman = mean(spearman, na.rm=TRUE),
    se_spearman = sd(spearman, na.rm=TRUE)/sqrt(n),
    mean_times = mean(times, na.rm=TRUE),
    se_times = sd(times, na.rm=TRUE)/sqrt(n),
  )
bandwidth_sum_data$sigma_param <- factor(bandwidth_sum_data$sigma_param, levels = c("1/sqrt(2)*ell", "1*ell", "1.5*ell"))
write.csv(bandwidth_sum_data, file.path(results_dir, "abl_bandwidths_aggregate_sum_data_driven.csv"))


# CV weight clipping sensitivity
# helper function
get_temp <- function(sets, rt, scn, c, sel_paths, sgsva_scores) {
  if (length(sets)==0) {
    temp <- data.frame(
      pathway = "no pathway passed filter",
      AUROC = NA_real_, MoranI = NA_real_,
      times = rt, ell = scn$ell,
      sigma_global = floor(scn$ell/sqrt(2)),
      sigma_f = scn$sigma_f, weight_clip = c,
      up_size = 0, down_size = 0,
      true_up = scn$up_size, true_down = scn$down_size,
      stringsAsFactors = FALSE
    )
  } else {
    temp <- data.frame(
      pathway = sel_paths,
      AUROC = NA_real_, MoranI = NA_real_,
      times = rt, ell = scn$ell,
      sigma_global = floor(scn$ell/sqrt(2)),
      sigma_f = scn$sigma_f, weight_clip = c,
      up_size = NA_real_, down_size = NA_real_,
      true_up = scn$up_size, true_down = scn$down_size,
      stringsAsFactors = FALSE
    )
    for (i in seq_along(sel_paths)) {
      p <- sel_paths[i]
      field <- sim$fields[[p]]
      truth <- top_mask(field, q=0.1)
      temp$AUROC[i] <- score_auc(truth, sgsva_scores[p,])
      temp$MoranI[i] <- moran.test(sgsva_scores[p,], listw)$estimate[1]
      temp$up_size[i] <- length(sets[[p]]$U)
      temp$down_size[i] <- length(sets[[p]]$D)
    }
  }
  temp
}

clip_param <- c(0.95, 0.99, 0.995)
res_clip_oracle <- NULL
res_clip_data <- NULL
for (r in 1:replicates) {
  sim_rep <- sim_res[[paste0("replicate",r)]]
  for (g in seq_len(length(sim_rep))) {
    scn <- sim_rep[[g]]$scn
    sim <- sim_rep[[g]]$sim
    for (c in clip_param) {
      path_names <- names(sim$per_path)
      # get sgsva scores
      t0 <- proc.time()[3]
      sgsva_res_oracle <- run_sgsva(sim$counts, sim$coords, sim$per_path, weight_clip = c, 
                                sigma_global = floor(scn$ell/sqrt(2)), directionality = "oracle", one_sided = "center_0.5")
      rt <- proc.time()[3] - t0
      sgsva_scores <- sgsva_res_oracle$z_score
      sets <- sgsva_res_oracle$gene_sets
      sel_paths <- rownames(sgsva_scores)
      temp_oracle <- get_temp(sets, rt, scn, c, sel_paths, sgsva_scores)
      
      # get sgsva scores
      t0 <- proc.time()[3]
      sgsva_res_data <- run_sgsva(sim$counts, sim$coords, gene_sets[path_names], weight_clip = c, 
                             sigma_global = floor(scn$ell/sqrt(2)), directionality = "data-driven", one_sided = "center_0.5")
      rt <- proc.time()[3] - t0
      sgsva_scores <- sgsva_res_data$z_score
      sets <- sgsva_res_data$gene_sets
      sel_paths <- rownames(sgsva_scores)
      temp_data <- get_temp(sets, rt, scn, c, sel_paths, sgsva_scores)
      
      res_clip_oracle <- rbind(res_clip_oracle, temp_oracle)
      res_clip_data <- rbind(res_clip_data, temp_data)
    }
  }
}

clip_sum_oracle <- res_clip_oracle %>%
  filter(!is.na(AUROC)) %>%
  mutate(
    set_type = case_when(
      true_up==25 & true_down==25 ~ "Balanced (25/25)",
      true_up==25 & true_down==0 ~ "One-sided (25/0)"
    )
  ) %>%
  group_by(pathway, ell, sigma_f, set_type, weight_clip) %>%
  summarise(
    n = n(),
    mean_auroc = mean(AUROC, na.rm=TRUE),
    se_auroc = sd(AUROC, na.rm=TRUE)/sqrt(n),
    mean_times = mean(times, na.rm=TRUE),
    se_times = sd(times, na.rm=TRUE)/sqrt(n),
    mean_moran = mean(MoranI, na.rm=TRUE),
    se_moran = sd(MoranI, na.rm=TRUE)/sqrt(n)
  )
write.csv(clip_sum_oracle, file.path(results_dir, "abl_weight_clipping_aggregate_sum_oracle.csv"))

clip_sum_data <- res_clip_data %>%
  filter(!is.na(AUROC)) %>%
  mutate(
    set_type = case_when(
      true_up==25 & true_down==25 ~ "Balanced (25/25)",
      true_up==25 & true_down==0 ~ "One-sided (25/0)"
    )
  ) %>%
  group_by(pathway, ell, sigma_f, set_type, weight_clip) %>%
  summarise(
    n = n(),
    mean_auroc = mean(AUROC, na.rm=TRUE),
    se_auroc = sd(AUROC, na.rm=TRUE)/sqrt(n),
    mean_times = mean(times, na.rm=TRUE),
    se_times = sd(times, na.rm=TRUE)/sqrt(n),
    mean_moran = mean(MoranI, na.rm=TRUE),
    se_moran = sd(MoranI, na.rm=TRUE)/sqrt(n)
  )
write.csv(clip_sum_data, file.path(results_dir, "abl_weight_clipping_aggregate_sum_data_driven.csv"))


# Spatial graph ablation (delaunay vs KNN (k=6,8,12))
mode_df <- data.frame(
  mode = c("delaunay", "knn", "knn", "knn"),
  k = c(0, 6, 8, 12)
)
res_graph <- NULL
for (l in ell) {
  sim <- sim_block_res[[paste0("block",l)]]
  
  for (i in 1:dim(mode_df)[1]) {
    mode <- mode_df$mode[i]
    k <- mode_df$k[i]
    
    # get sgsva scores
    sgsva_res <- run_sgsva(sim$counts, sim$coords, sim$per_path,
                           mode = mode, k = k, sigma_global = floor(l/sqrt(2)), 
                           directionality = "oracle", one_sided = "center_0.5")
    sgsva_scores <- sgsva_res$z_score
    
    hid <- make_hex_ids(coords, cellsize = l)
    df <- vector("list", length(sel_paths))
    
    for (p in seq_along(sel_paths)) {
      pathway <- sel_paths[p]
      x <- sgsva_scores[pathway,]
      
      set.seed(123)
      # get hex-block results
      out_moran_hex <- perm_with_early_stop(stat_moran_hex, x, listw, hid)
      # convert results into data.frame
      df[[p]] <- data.frame(
        pathway = pathway, obs_moran = out_moran_hex$obs,
        z_moran_hex = out_moran_hex$z, p_moran_hex = out_moran_hex$p,
        mode = mode, k=k,
        block_size = l, ell = l,
        sigma_global = floor(l/sqrt(2)),
        stringsAsFactors = FALSE
      )
    }
    temp_spatial <- do.call(rbind, df)
    rownames(temp_spatial) <- NULL
    res_graph <- rbind(res_graph, temp_spatial)
  }
}
write.csv(res_graph, file.path(results_dir, "abl_spatial_graph_moran_oracle.csv"))


# hex block stability
l <- 200
block_sizes <- c(150, 200, 300)
sim <- sim_block_res[[paste0("block",l)]]
# get sgsva scores
sgsva_res <- run_sgsva(sim$counts, sim$coords, sim$per_path, sigma_global = floor(l/sqrt(2)), 
                       directionality = "oracle", one_sided = "center_0.5")
sgsva_scores <- sgsva_res$z_score
res_block <- NULL
for (b in block_sizes) {
  hid <- make_hex_ids(coords, cellsize = b)
  df <- vector("list", length(sel_paths))
  
  set.seed(123)
  for (p in seq_along(sel_paths)) {
    pathway <- sel_paths[p]
    x <- sgsva_scores[pathway,]
    
    # get hex-block results
    out_moran_hex <- perm_with_early_stop(stat_moran_hex, x, listw, hid)
    out_geary_hex <- perm_with_early_stop(stat_geary_hex, x, listw, hid)
    # convert results into data.frame
    df[[p]] <- data.frame(
      pathway = pathway,
      obs_moran = out_moran_hex$obs, z_moran_hex = out_moran_hex$z, 
      p_moran_hex = out_moran_hex$p, obs_geary = -out_geary_hex$obs, 
      z_geary_hex = -out_geary_hex$z, p_geary_hex = out_geary_hex$p,
      block_size = b, ell = l,
      sigma_global = floor(l/sqrt(2)),
      stringsAsFactors = FALSE
    )
  }
  temp_spatial <- do.call(rbind, df)
  rownames(temp_spatial) <- NULL
  res_block <- rbind(res_block, temp_spatial)
}
write.csv(res_block, file.path(results_dir, "abl_block_stability_autocorrelation_oracle.csv"))


