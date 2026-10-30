# Robustness analysis
suppressPackageStartupMessages ({
  library(Matrix) ; library(matrixStats) ; library(MASS)
  library(dplyr) ; library(tidyr) ; library(purrr) ; 
  library(readr); library(spdep) ; library(sf) ; 
  library(pROC); library(msigdbr) ; library(ggplot2); 
  library(Seurat); library(patchwork)
})

set.seed(1)

# ---------------------------
# paths and reproducibility
# ---------------------------
project_dir <- "."
results_dir <- file.path(project_dir, "results", "robustness")
sim_dir <- file.path(project_dir, "results", "simulation")
supp_dir <- file.path(project_dir, "results", "supplementary_data")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_dir, "R/spatial_permutation.R"))
source(file.path(project_dir, "R/simulation.R"))
source(file.path(project_dir, "R/simulation_scoring.R"))
source(file.path(project_dir, "R/evaluation_metrics.R"))


# ---------------------------
# generate simulation 
# ---------------------------
# ---- gene sets ( MSigDB KEGG legacy ) ----
msig <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
gene_sets <- split(msig$gene_symbol, msig$gs_name)

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

ell <- c(50, 200, 400)
sel_paths <- names(gene_sets)

# simulation for varied ell with all KEGG pathways
sim_block_res <- NULL
base_seed <- 2000L
for (i in seq_along(ell)) {
  set.seed(base_seed + i)
  l <- ell[i]
  sim <- simulate_counts(coords, gene_sets, sel_paths, ell=l)
  sim_block_res[[paste0("block",l)]] <- sim
}
saveRDS(sim_block_res, file.path(results_dir, "sim_block_res.rds"))

# read all simulations
sim_block_res <- readRDS(file.path(results_dir, "sim_block_res.rds"))
sim_res <- readRDS(file.path(sim_dir, "rep_simulation.rds"))

# ---------------------------
# ablation analysis
# ---------------------------
# correlation threshold for directionality
rho_levels <- c(0, 0.05, 0.10, 0.20)
res_rho <- NULL
for (r in seq_len(replicates)) {
  
  sim_rep <- sim_res[[paste0("replicate", r)]]
  
  for (g in seq_along(sim_rep)) {
    
    scn <- sim_rep[[g]]$scn
    sim <- sim_rep[[g]]$sim
    
    # locCPM and filter the expression matrix
    counts <- sim$counts
    G <- nrow(counts); N <- ncol(counts)
    lib <- colSums(counts)+1
    lcp <- log1p(t(t(counts)/lib*1e6)) #logCPM
    expr <- as.matrix(lcp)
    keep <- rowSums(expr > 0) >= floor(0.05 * ncol(expr))
    expr <- expr[keep, , drop = FALSE]
    path_names <- names(sim$per_path)
    
    for (rho in rho_levels) {
      
      t0 <- proc.time()[3]
      spa_res <- run_spargs(expr, sim$coords, gene_sets[path_names],
        mode = "delaunay", sigma_global = floor(scn$ell / sqrt(2)),
        directionality = "data-driven", cor_min = rho, one_sided = "center_0.5"
      )
      runtime <- proc.time()[3] - t0
      
      spa_scores <- spa_res$z_score
      retained_paths <- rownames(spa_scores)
      
      # Retention should be recorded even when a pathway is not retained
      for (p in path_names) {
        retained <- p %in% retained_paths
        auc <- NA_real_
        
        if (retained) {
          truth <- top_mask(sim$fields[[p]], q = 0.1)
          auc <- score_auc_invariant(truth, spa_scores[p, ])
        }
        
        res_rho <- rbind(
          res_rho,
          data.frame(
            replicate = r, pathway = p, ell = scn$ell,
            sigma_f = scn$sigma_f, true_up = scn$up_size,
            true_down = scn$down_size, rho = rho,
            retained = retained, AUROC = auc,
            runtime_seconds = runtime,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}
write.csv(res_rho, file.path(results_dir, "abl_correlation_threshold.csv"), row.names = FALSE)

# bandwidth sensitivity 
param <- c(1/sqrt(2), 1, 1.5) # (1/sqrt(2), 1, 1.5)*ell
res_sigma_oracle <- NULL
res_sigma_data <- NULL
for (r in 1:replicates) {
  sim_rep <- sim_res[[paste0("replicate",r)]]
  for (g in seq_len(length(sim_rep))) {
    scn <- sim_rep[[g]]$scn
    sim <- sim_rep[[g]]$sim
    
    # locCPM and filter the expression matrix
    counts <- sim$counts
    G <- nrow(counts); N <- ncol(counts)
    lib <- colSums(counts)+1
    lcp <- log1p(t(t(counts)/lib*1e6)) #logCPM
    expr <- as.matrix(lcp)
    keep <- rowSums(expr > 0) >= floor(0.05 * ncol(expr))
    expr <- expr[keep, , drop = FALSE]
    
    for (p in param) {
      sigma_global <- floor(scn$ell*p)
      path_names <- names(sim$per_path)
      
      # get spargs scores
      t0 <- proc.time()[3]
      spa_res_oracle <- run_spargs(expr, sim$coords, sim$per_path, mode = "delaunay",
                                   sigma_global = sigma_global, directionality = "oracle", one_sided = "center_0.5")
      rt <- proc.time()[3] - t0
      spa_scores <- spa_res_oracle$z_score
      sets <- spa_res_oracle$gene_sets
      temp <- eval_metrics_abl(score = spa_scores, sets, method = "SpaRGS", scn, sim, rt)
      temp$sigma_param <- p
      temp$sigma_global <- sigma_global
      res_sigma_oracle <- rbind(res_sigma_oracle, temp)
      
      # get spargs scores
      t0 <- proc.time()[3]
      spa_res_data <- run_spargs(expr, sim$coords, gene_sets[path_names], mode = "delaunay",
                                 sigma_global = sigma_global, directionality = "data-driven", one_sided = "center_0.5")
      rt <- proc.time()[3] - t0
      spa_scores <- spa_res_data$z_score
      sets <- spa_res_data$gene_sets
      temp <- eval_metrics_abl(score = spa_scores, sets, method = "SpaRGS", scn, sim, rt)
      temp$sigma_param <- p
      temp$sigma_global <- sigma_global
      res_sigma_data <- rbind(res_sigma_data, temp)
    }
  }
}
write.csv(res_sigma_oracle, file.path(results_dir, "abl_bandwidths_all_replicates_oracle.csv"), row.names = FALSE)
write.csv(res_sigma_data, file.path(results_dir, "abl_bandwidths_all_replicates_data_driven.csv"), row.names = FALSE)

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
    mean_spearman_signed = mean(spearman_signed, na.rm=TRUE),
    se_spearman_signed = sd(spearman_signed, na.rm=TRUE)/sqrt(n),
    mean_times = mean(times, na.rm=TRUE),
    se_times = sd(times, na.rm=TRUE)/sqrt(n),
    .groups = "drop"
  )
bandwidth_sum_oracle$sigma_param <- factor(bandwidth_sum_oracle$sigma_param, levels = c("1/sqrt(2)*ell", "1*ell", "1.5*ell"))
write.csv(bandwidth_sum_oracle, file.path(results_dir, "abl_bandwidths_aggregate_sum_oracle.csv"), row.names = FALSE)

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
    mean_spearman_signed = mean(spearman_signed, na.rm=TRUE),
    se_spearman_signed = sd(spearman_signed, na.rm=TRUE)/sqrt(n),
    mean_times = mean(times, na.rm=TRUE),
    se_times = sd(times, na.rm=TRUE)/sqrt(n),
    .groups = "drop"
  )
bandwidth_sum_data$sigma_param <- factor(bandwidth_sum_data$sigma_param, levels = c("1/sqrt(2)*ell", "1*ell", "1.5*ell"))
write.csv(bandwidth_sum_data, file.path(results_dir, "abl_bandwidths_aggregate_sum_data_driven.csv"), row.names = FALSE)


# CV weight clipping sensitivity
# helper function
get_temp <- function(sets, rt, scn, c, sel_paths, spa_scores, sim, listw) {
  if (length(sets)==0) {
    temp <- data.frame(
      pathway = "no pathway passed filter",
      AUROC = NA_real_, MoranI = NA_real_,
      spearman = NA_real_, spearman_signed = NA_real_,
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
      spearman = NA_real_, spearman_signed = NA_real_,
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
      temp$AUROC[i] <- score_auc_invariant(truth, spa_scores[p,])
      temp$MoranI[i] <- moran.test(spa_scores[p,], listw)$estimate[1]
      rho_signed <- suppressWarnings(
        cor(field, spa_scores[p, ], method = "spearman", use = "pairwise.complete.obs")
      )
      temp$spearman_signed[i] <- rho_signed
      temp$spearman[i] <- abs(rho_signed)
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
    
    # locCPM and filter the expression matrix
    counts <- sim$counts
    G <- nrow(counts); N <- ncol(counts)
    lib <- colSums(counts)+1
    lcp <- log1p(t(t(counts)/lib*1e6)) #logCPM
    expr <- as.matrix(lcp)
    keep <- rowSums(expr > 0) >= floor(0.05 * ncol(expr))
    expr <- expr[keep, , drop = FALSE]
    for (c in clip_param) {
      path_names <- names(sim$per_path)
      
      # get spargs scores
      t0 <- proc.time()[3]
      spa_res_oracle <- run_spargs(expr, sim$coords, sim$per_path, weight_clip = c, mode = "delaunay",
                                    sigma_global = floor(scn$ell/sqrt(2)), directionality = "oracle", one_sided = "center_0.5")
      rt <- proc.time()[3] - t0
      spa_scores <- spa_res_oracle$z_score
      sets <- spa_res_oracle$gene_sets
      sel_paths <- rownames(spa_scores)
      temp_oracle <- get_temp(sets, rt, scn, c, sel_paths, spa_scores, sim, listw)
      
      # get spargs scores
      t0 <- proc.time()[3]
      spa_res_data <- run_spargs(expr, sim$coords, gene_sets[path_names], weight_clip = c, mode = "delaunay",
                                  sigma_global = floor(scn$ell/sqrt(2)), directionality = "data-driven", one_sided = "center_0.5")
      rt <- proc.time()[3] - t0
      spa_scores <- spa_res_data$z_score
      sets <- spa_res_data$gene_sets
      sel_paths <- rownames(spa_scores)
      temp_data <- get_temp(sets, rt, scn, c, sel_paths, spa_scores, sim, listw)
      
      res_clip_oracle <- rbind(res_clip_oracle, temp_oracle)
      res_clip_data <- rbind(res_clip_data, temp_data)
    }
  }
}
write.csv(res_clip_oracle, file.path(results_dir, "abl_weight_clipping_all_replicates_oracle.csv"), row.names = FALSE)
write.csv(res_clip_data, file.path(results_dir, "abl_weight_clipping_all_replicates_data_driven.csv"), row.names = FALSE)


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
    mean_spearman = mean(spearman, na.rm=TRUE),
    se_spearman = sd(spearman, na.rm=TRUE)/sqrt(n),
    mean_spearman_signed = mean(spearman_signed, na.rm=TRUE),
    se_spearman_signed = sd(spearman_signed, na.rm=TRUE)/sqrt(n),
    mean_times = mean(times, na.rm=TRUE),
    se_times = sd(times, na.rm=TRUE)/sqrt(n),
    mean_moran = mean(MoranI, na.rm=TRUE),
    se_moran = sd(MoranI, na.rm=TRUE)/sqrt(n),
    .groups = "drop"
  )
write.csv(clip_sum_oracle, file.path(results_dir, "abl_weight_clipping_aggregate_sum_oracle.csv"), row.names = FALSE)

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
    mean_spearman = mean(spearman, na.rm=TRUE),
    se_spearman = sd(spearman, na.rm=TRUE)/sqrt(n),
    mean_spearman_signed = mean(spearman_signed, na.rm=TRUE),
    se_spearman_signed = sd(spearman_signed, na.rm=TRUE)/sqrt(n),
    mean_times = mean(times, na.rm=TRUE),
    se_times = sd(times, na.rm=TRUE)/sqrt(n),
    mean_moran = mean(MoranI, na.rm=TRUE),
    se_moran = sd(MoranI, na.rm=TRUE)/sqrt(n),
    .groups = "drop"
  )
write.csv(clip_sum_data, file.path(results_dir, "abl_weight_clipping_aggregate_sum_data_driven.csv"), row.names = FALSE)


# Spatial graph ablation (delaunay vs KNN (k=6,8,12))
mode_df <- data.frame(
  mode = c("delaunay", "knn", "knn", "knn"),
  k = c(0, 6, 8, 12)
)
res_graph <- NULL
for (l in ell) {
  sim <- sim_block_res[[paste0("block",l)]]
  test_nb <- spdep::tri2nb(as.matrix(sim$coords))
  test_listw <- spdep::nb2listw(test_nb, style = "W")
  hid <- make_hex_ids(sim$coords, cellsize = l)
  
  # locCPM and filter the expression matrix
  counts <- sim$counts
  G <- nrow(counts); N <- ncol(counts)
  lib <- colSums(counts)+1
  lcp <- log1p(t(t(counts)/lib*1e6)) #logCPM
  expr <- as.matrix(lcp)
  keep <- rowSums(expr > 0) >= floor(0.05 * ncol(expr))
  expr <- expr[keep, , drop = FALSE]
  
  for (i in 1:dim(mode_df)[1]) {
    mode <- mode_df$mode[i]
    k <- mode_df$k[i]
    
    # get spargs scores (use oracle to isolate graph effects)
    spa_res <- run_spargs(expr, sim$coords, sim$per_path,
                           mode = mode, k = k, sigma_global = floor(l/sqrt(2)), 
                           directionality = "oracle", one_sided = "center_0.5")
    spa_scores <- spa_res$z_score
    sel_paths <- names(sim$per_path)
    
    df <- vector("list", length(sel_paths))
    
    for (p in seq_along(sel_paths)) {
      pathway <- sel_paths[p]
      x <- spa_scores[pathway,]
      
      obs_moran <- stat_moran(x, test_listw)
      
      df[[p]] <- data.frame(
        pathway = pathway,
        obs_moran = obs_moran,
        mode = mode,
        k = k, ell = l,
        sigma_global = floor(l / sqrt(2)),
        stringsAsFactors = FALSE
      )
    }
    temp_spatial <- do.call(rbind, df)
    rownames(temp_spatial) <- NULL
    res_graph <- rbind(res_graph, temp_spatial)
  }
}
write.csv(res_graph, 
          file.path(results_dir, "abl_spatial_graph_moran_oracle.csv"), row.names = FALSE)


# hex block stability
l <- 200
block_sizes <- c(150, 200, 300)
sim <- sim_block_res[[paste0("block",l)]]
# locCPM and filter the expression matrix
counts <- sim$counts
G <- nrow(counts); N <- ncol(counts)
lib <- colSums(counts)+1
lcp <- log1p(t(t(counts)/lib*1e6)) #logCPM
expr <- as.matrix(lcp)
keep <- rowSums(expr > 0) >= floor(0.05 * ncol(expr))
expr <- expr[keep, , drop = FALSE]
# get spargs scores
spa_res <- run_spargs(expr, sim$coords, sim$per_path, mode = "delaunay",
                      sigma_global = floor(l/sqrt(2)), directionality = "oracle", one_sided = "center_0.5")
spa_scores <- spa_res$z_score
sel_paths <- names(sim$per_path)
test_nb <- spdep::tri2nb(as.matrix(sim$coords))
test_listw <- spdep::nb2listw(test_nb, style = "W")

res_block <- NULL
for (b in block_sizes) {
  hid <- make_hex_ids(sim$coords, cellsize = b)
  df <- vector("list", length(sel_paths))

  for (p in seq_along(sel_paths)) {
    pathway <- sel_paths[p]
    x <- spa_scores[pathway,]
    
    set.seed(123 + p)
    # get hex-block results
    out_moran_hex <- perm_with_early_stop(stat_moran_hex, x, test_listw, hid)
    
    # convert results into data.frame
    df[[p]] <- data.frame(
      pathway = pathway,
      obs_moran = out_moran_hex$obs, z_moran_hex = out_moran_hex$z, 
      p_moran_hex = out_moran_hex$p,
      block_size = b, ell = l,
      sigma_global = floor(l/sqrt(2)),
      stringsAsFactors = FALSE
    )
  }
  temp_spatial <- do.call(rbind, df)
  rownames(temp_spatial) <- NULL
  res_block <- rbind(res_block, temp_spatial)
}
write.csv(res_block, 
          file.path(results_dir, "abl_block_stability_autocorrelation_oracle.csv"), row.names = FALSE)


# ---------------------------
# reproducibility
# ---------------------------
writeLines(
  capture.output(sessionInfo()),
  con = file.path(results_dir, "Robustness_sessionInfo.txt")
)


library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "SData2_directionality")
writeData(wb, "SData2_directionality", res_rho)
addWorksheet(wb, "SData2_bandwidth_data_driven")
writeData(wb, "SData2_bandwidth_data_driven", res_sigma_data)
addWorksheet(wb, "SData2_bandwidth_oracle")
writeData(wb, "SData2_bandwidth_oracle", res_sigma_oracle)
addWorksheet(wb, "SData2_clipping_data_driven")
writeData(wb, "SData2_clipping_data_driven", res_clip_data)
addWorksheet(wb, "SData2_clipping_oracle")
writeData(wb, "SData2_clipping_oracle", res_clip_oracle)
addWorksheet(wb, "SData2_graph")
writeData(wb, "SData2_graph", res_graph)
addWorksheet(wb, "SData2_hex_block")
writeData(wb, "SData2_hex_block", res_block)
saveWorkbook(wb, file.path(supp_dir, "Supplementary Data 2.xlsx"), overwrite = TRUE)



