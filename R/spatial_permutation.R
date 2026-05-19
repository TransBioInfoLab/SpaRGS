# R/spatial_permutation.R

### --- spatial permutation function ---
spatial_permutation <- function(
    z_sgsva, coords, balanced_sets,
    # sgsva parameters
    mode = c("delaunay","knn"), k=6, # mode for neighborhood
    sigma_global = 200, # global smoothing bandwidth (um in lowres coords)
    weight_clip = 0.99, # CV weight clipping quantile
    cor_min = 0.10, # abs(Spearman) threshold for directionality
    # scoring
    one_sided = c("off", "center_0.5", "center_median"),
    # permutation
    perm_methods = c("moran_mc", "geary_mc", "hex_block"),
    n_perm_max = 10000, n_perm_min = 1000, B = 10000,
    early_stop = TRUE, alpha = 0.05, 
    path_names = NULL, # option to choose pathway histograms
    hex_cellsize = 200, # um for block permutation
    plot_dir = NULL, # plot folder
    seed_permutation = 124 # RNG seed
){
  mode <- match.arg(mode)
  one_sided <- match.arg(one_sided)
  set.seed(seed_permutation)
  
  make_plots <- !is.null(plot_dir) && nzchar(plot_dir)
  if (make_plots) {
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # ----- Build spatial neighbors on BIN centroids -----
  coords <- as.data.frame(coords)
  stopifnot(nrow(coords) == ncol(z_sgsva))
  
  if (mode == "delaunay") {
    nb <- spdep::tri2nb(coords)
  } else if (mode == "knn") {
    knn_obj <- spdep::knearneigh(coords, k = k)
    nb <- spdep::knn2nb(knn_obj)
  }
  
  # drop isolated bins with no neighbors
  iso <- which(spdep::card(nb) == 0L)
  if (length(iso) > 0L) {
    message("Dropping ", length(iso), " island spots with no neighbors.")
    keep_idx <- setdiff(seq_len(nrow(coords)), iso)
    coords <- coords[keep_idx, , drop=FALSE]
    z_sgsva <- z_sgsva[, keep_idx, drop=FALSE]
    if (mode == "delaunay") {
      nb <- spdep::tri2nb(coords)
    } else {
      knn_obj <- spdep::knearneigh(coords_bin, k = k)
      nb <- spdep::knn2nb(knn_obj)
    }
  }
  lw <- spdep::nb2listw(nb, style="W")
  
  # --- Spatial permutation ---
  scores <- z_sgsva
  permutation_results <- list()
  if (early_stop==TRUE) {
    n_pathways <- nrow(scores)
    
    if (("moran_mc" %in% perm_methods) & ("geary_mc" %in% perm_methods)) {
      mc_pathway_results <- vector("list", nrow(scores))
      mc_perm_null_results <- list()
      
      for (p in 1:n_pathways) {
        pathway <- rownames(scores)[p]
        x <- as.numeric(scores[pathway,])
        
        # get moran.mc/geary.mc results
        set.seed(seed_permutation)
        out_moran_mc <- perm_with_early_stop(stat_moran, x, lw)
        out_geary_mc <- perm_with_early_stop(stat_geary, x, lw)
        
        # convert results into data.frame
        mc_pathway_results[[p]] <- data.frame(
          pathway = pathway,
          obs_moran = out_moran_mc$obs, z_moran = out_moran_mc$z,
          p_moran = out_moran_mc$p, p_adj_moran = NA_real_,
          null_mean_moran = out_moran_mc$null_mean,
          null_sd_moran = out_moran_mc$null_sd,
          t_moran = out_moran_mc$t, b_moran = out_moran_mc$b,
          stop_reason_moran = out_moran_mc$stop_reason,
          obs_geary = -out_geary_mc$obs, z_geary = -out_geary_mc$z,
          p_geary = out_geary_mc$p, p_adj_geary = NA_real_,
          null_mean_geary = -out_geary_mc$null_mean,
          null_sd_geary = out_geary_mc$null_sd,
          obs_geary_flipped = out_geary_mc$obs, z_geary_flipped = out_geary_mc$z,
          null_mean_geary_flipped = out_geary_mc$null_mean,
          null_sd_geary_flipped = out_geary_mc$null_sd,
          t_geary = out_geary_mc$t, b_geary = out_geary_mc$b,
          stop_reason_geary = out_geary_mc$stop_reason,
          up = length(balanced_sets[[pathway]]$U),
          down = length(balanced_sets[[pathway]]$D),
          cv_clip = weight_clip, global_sigma = sigma_global,
          graph = mode, cor_min = cor_min,
          stringsAsFactors = FALSE
        )
        
        # store permutation results
        mc_perm_null_results[[pathway]] <- list(
          moran_perm = out_moran_mc$perm_stats,
          geary_perm = -out_geary_mc$perm_stats
        )
      }
      
      # Combine all results into a single data frame
      mc_results_df <- do.call(rbind, mc_pathway_results)
      # Adjust p-values (FDR)
      mc_results_df$p_adj_moran <- p.adjust(mc_results_df$p_moran, method = "BH")
      mc_results_df$p_adj_geary <- p.adjust(mc_results_df$p_geary, method = "BH")
      rownames(mc_results_df) <- NULL
      
      permutation_results[["moran_geary_results"]] <- list(
        moran_geary_summary = mc_results_df,
        mc_perm_null_results = mc_perm_null_results
      )
      
      if (make_plots) {
        # Example: save a few histograms for Moran's I
        hist_plot <- list()
        if (is.null(path_names)){
          top4_moran <- mc_results_df$pathway[order(mc_results_df$obs_moran, decreasing=TRUE)[1:4]]
        } else {
          top4_moran <- path_names
        }
        for (nm in top4_moran) {
          df <- data.frame(val = mc_perm_null_results[[nm]]$moran_perm)
          p <- ggplot(df, aes(x = val)) +
            geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
            geom_vline(xintercept = mc_results_df$obs_moran[mc_results_df$pathway==nm], color = "firebrick", linetype = 2, linewidth = 1) +
            theme_bw() +
            labs(x = "Moran's I (permuted)", y = "Count",
                 title = nm)
          hist_plot[[nm]] <- p
        }
        pdf(file.path(plot_dir, "mc_perm_null_hist_moran.pdf"), width=12, height=8)
        gridExtra::grid.arrange(grobs=hist_plot, ncol=2)
        dev.off()
        
        # Example: save a few histograms for Geary's C
        hist_plot <- list()
        if (is.null(path_names)) {
          top4_geary <- mc_results_df$pathway[order(mc_results_df$obs_geary, decreasing=FALSE)[1:4]]
        } else {
          top4_geary <- path_names
        }
        for (nm in top4_geary) {
          df <- data.frame(val = mc_perm_null_results[[nm]]$geary_perm)
          p <- ggplot(df, aes(x = val)) +
            geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
            geom_vline(xintercept = mc_results_df$obs_geary[mc_results_df$pathway==nm], color = "firebrick", linetype = 2, linewidth = 1) +
            theme_bw() +
            labs(x = "Geary's C (permuted)", y = "Count",
                 title = nm)
          hist_plot[[nm]] <- p
        }
        pdf(file.path(plot_dir, "mc_perm_null_hist_geary.pdf"), width=12, height=8)
        gridExtra::grid.arrange(grobs=hist_plot, ncol=2)
        dev.off()
      }
    }
    
    if ("hex_block" %in% perm_methods) {
      hex_pathway_results <- vector("list", nrow(scores))
      hex_perm_null_results <- list()
      hex_id <- make_hex_ids(coords, cellsize = hex_cellsize)
      
      for (p in 1:n_pathways) {
        pathway <- rownames(scores)[p]
        x <- as.numeric(scores[pathway,])
        
        # get hex-block results
        set.seed(seed_permutation)
        out_moran_hex <- perm_with_early_stop(stat_moran_hex, x, lw, hex_id)
        out_geary_hex <- perm_with_early_stop(stat_geary_hex, x, lw, hex_id)
        
        # convert results into data.frame
        hex_pathway_results[[p]] <- data.frame(
          pathway = pathway,
          obs_moran = out_moran_hex$obs, z_moran = out_moran_hex$z,
          p_moran = out_moran_hex$p, p_adj_moran = NA_real_,
          null_mean_moran = out_moran_hex$null_mean,
          null_sd_moran = out_moran_hex$null_sd,
          t_moran = out_moran_hex$t, b_moran = out_moran_hex$b,
          stop_reason_moran = out_moran_hex$stop_reason,
          obs_geary = -out_geary_hex$obs, z_geary = -out_geary_hex$z,
          p_geary = out_geary_hex$p, p_adj_geary = NA_real_,
          null_mean_geary = -out_geary_hex$null_mean,
          null_sd_geary = out_geary_hex$null_sd,
          obs_geary_flipped = out_geary_hex$obs, z_geary_flipped = out_geary_hex$z,
          null_mean_geary_flipped = out_geary_hex$null_mean,
          null_sd_geary_flipped = out_geary_hex$null_sd,
          t_geary = out_geary_hex$t, b_geary = out_geary_hex$b,
          stop_reason_geary = out_geary_hex$stop_reason,
          block_size = hex_cellsize,
          up = length(balanced_sets[[pathway]]$U),
          down = length(balanced_sets[[pathway]]$D),
          cv_clip = weight_clip, global_sigma = sigma_global,
          graph = mode, cor_min = cor_min,
          stringsAsFactors = FALSE
        )
        
        # store permutation results
        hex_perm_null_results[[pathway]] <- list(
          moran_perm = out_moran_hex$perm_stats,
          geary_perm = -out_geary_hex$perm_stats
        )
        
      }
      # Combine all results into a single data frame
      hex_results_df <- do.call(rbind, hex_pathway_results)
      # Adjust p-values (FDR)
      hex_results_df$p_adj_moran <- p.adjust(hex_results_df$p_moran, method = "BH")
      hex_results_df$p_adj_geary <- p.adjust(hex_results_df$p_geary, method = "BH")
      rownames(hex_results_df) <- NULL
      
      permutation_results[["hex_block_results"]] <- list(
        hex_block_summary = hex_results_df,
        hex_perm_null_results = hex_perm_null_results
      )
      
      if (make_plots) {
        # Example: save a few histograms for Moran's I
        hist_plot <- list()
        if (is.null(path_names)) {
          top4_moran <- hex_results_df$pathway[order(hex_results_df$obs_moran, decreasing=TRUE)[1:4]]
        } else {
          top4_moran <- path_names
        }
        for (nm in top4_moran) {
          df <- data.frame(val = hex_perm_null_results[[nm]]$moran_perm)
          p <- ggplot(df, aes(x = val)) +
            geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
            geom_vline(xintercept = hex_results_df$obs_moran[hex_results_df$pathway==nm], color = "firebrick", linetype = 2, linewidth = 1) +
            theme_bw() +
            labs(x = "Moran's I (permuted)", y = "Count",
                 title = nm)
          hist_plot[[nm]] <- p
        }
        pdf(file.path(plot_dir, paste0("block", hex_cellsize, "_null_hist_moran.pdf")), width=12, height=8)
        gridExtra::grid.arrange(grobs=hist_plot, ncol=2)
        dev.off()
        
        # Example: save a few histograms for Geary's C
        hist_plot <- list()
        if (is.null(path_names)) {
          top4_geary <- hex_results_df$pathway[order(hex_results_df$obs_geary, decreasing=FALSE)[1:4]]
        } else {
          top4_geary <- path_names
        }
        for (nm in top4_geary) {
          df <- data.frame(val = hex_perm_null_results[[nm]]$geary_perm)
          p <- ggplot(df, aes(x = val)) +
            geom_histogram(bins = 40, fill = "grey80", color = "grey50") +
            geom_vline(xintercept = hex_results_df$obs_geary[hex_results_df$pathway==nm], color = "firebrick", linetype = 2, linewidth = 1) +
            theme_bw() +
            labs(x = "Geary's C (permuted)", y = "Count",
                 title = nm)
          hist_plot[[nm]] <- p
        }
        pdf(file.path(plot_dir, paste0("block", hex_cellsize, "_null_hist_geary.pdf")), width=12, height=8)
        gridExtra::grid.arrange(grobs=hist_plot, ncol=2)
        dev.off()
      }
    }
    
  } else {
    
    # no early stopping
    if (("moran_mc" %in% perm_methods) & ("geary_mc" %in% perm_methods)) {
      set.seed(seed_permutation)
      mc_perm_res <- run_moran_geary(scores, w = lw, n_perm = B,
                                     balanced_sets, pathway_sigmas)
      mc_results_df <- mc_perm_res$res
      mc_perm_null_results <- mc_perm_res$perm_null_res
      
      permutation_results[["moran_geary_results"]] <- list(
        moran_geary_summary = mc_results_df,
        mc_perm_null_results = mc_perm_null_results
      )
      
    }
    
    if ("hex_block" %in% perm_methods) {
      set.seed(seed_permutation)
      hex_perm_res <- hex_block_perm(scores, lw, coords, sets = balanced_sets, 
                                     block_size = hex_cellsize, B = B)
      hex_results_df <- hex_perm_res$res
      hex_perm_null_results <- hex_perm_res$perm_null_res
      
      permutation_results[["hex_block_results"]] <- list(
        hex_block_summary = hex_results_df,
        hex_perm_null_results = hex_perm_null_results
      )
      
    }
  }
  
  return(permutation_results)
}

# --- helper function ---
# make a hexagonal grid
make_hex_ids <- function(coords, cellsize = 200) {
  pts <- st_as_sf(coords, coords = c("x","y"), crs = NA)
  hex <- st_make_grid(pts, cellsize = cellsize, square = FALSE)
  hid <- st_within(pts, hex)
  sapply(hid, function(x) if (length(x) == 0) NA_integer_ else x[1])
}

# hex block permutation
hex_block_perm <- function(z, lw, coords, sets, block_size=200, 
                           B = 1000, weight_clip=0.99, sigma_global=207) {
  
  hex_id <- make_hex_ids(coords, cellsize = block_size)
  path_names <- rownames(z)
  res <- data.frame(
    pathway = path_names,
    obs_moran = NA_real_, z_moran = NA_real_,
    p_moran = NA_real_, p_adj_moran = NA_real_,
    null_mean_moran = NA_real_, null_sd_moran = NA_real_, 
    obs_geary = NA_real_, z_geary = NA_real_,
    p_geary = NA_real_, p_adj_geary = NA_real_,
    null_mean_geary = NA_real_, null_sd_geary = NA_real_,
    up = NA_integer_, down = NA_integer_,
    block_size = block_size, B = B,
    cv_clip = weight_clip, sigma_global = sigma_global,
    optimal_sigma = NA_real_,
    stringsAsFactors = FALSE
  )
  
  null_hist <- vector("list", length(path_names))
  names(null_hist) <- path_names
  
  for (i in seq_along(path_names)) {
    pth <- path_names[i];s <- as.numeric(z[pth, ])
    
    obs_I <- moran.test(s, lw)$estimate[1]
    obs_C <- geary.test(s, lw)$estimate[1]
    
    perm_I <- replicate(B, {
      x <- s
      for (b in unique(na.omit(hex_id))) {
        idx <- which(hex_id == b)
        if (length(idx) > 1) x[idx] <- sample(x[idx])
      }
      moran.test(x, lw)$estimate[1]
    })
    perm_C <- replicate(B, {
      x <- s
      for (b in unique(na.omit(hex_id))) {
        idx <- which(hex_id == b)
        if (length(idx) > 1) x[idx] <- sample(x[idx])
      }
      geary.test(x, lw)$estimate[1]
    })
    
    mu_moran <- mean(perm_I); sd_moran <- sd(perm_I)
    z_moran <- if (sd_moran > 0) (obs_I - mu_moran) / sd_moran else NA_real_
    mu_geary <- mean(perm_C); sd_geary <- sd(perm_C)
    z_geary <- if (sd_geary > 0) (obs_C - mu_geary) / sd_geary else NA_real_
    
    # Unbiased p-values
    res$obs_moran[i] <- obs_I; res$z_moran[i] <- z_moran
    res$p_moran[i] <- mean(c(obs_I, perm_I)>=obs_I)
    res$null_mean_moran[i] <- mu_moran
    res$null_sd_moran[i] <- sd_moran
    res$obs_geary[i] <- obs_C; res$z_geary[i] <- z_geary
    res$null_mean_geary[i] <- mu_geary
    res$null_sd_geary[i] <- sd_geary
    res$p_geary[i] <- mean(c(obs_C, perm_C)<=obs_C)
    res$up[i] <- length(sets[[pth]]$U)
    res$down[i] <- length(sets[[pth]]$D)
    res$optimal_sigma[i] <- pathway_sigmas[[pth]]$median_sigma
    
    null_hist[[pth]] <- list(moran_perm = perm_I, geary_perm = perm_C)
  }
  
  # Adjust p-values (FDR)
  pval_moran <- p.adjust(res$p_moran, method = "BH")
  pval_geary <- p.adjust(res$p_geary, method = "BH")
  res$p_adj_moran <- ifelse(pval_moran < 1/(B+1), paste0("p < ", formatC(1/(B+1), format = "e", digits = 1)), pval_moran)
  res$p_adj_geary <- ifelse(pval_geary < 1/(B+1), paste0("p < ", formatC(1/(B+1), format = "e", digits = 1)), pval_geary)
  
  return(
    list(res = res, perm_null_res = null_hist)
  )
}

# moran.mc/geary.mc
run_moran_geary <- function(scores, w, B = 1000,
                            sets, weight_clip=0.99, sigma_global=207) {
  # Extract scores for selected pathways
  n_pathways <- nrow(scores)
  paths <- rownames(scores)
  n_perm <- B
  null_hist <- vector("list", length(paths))
  names(null_hist) <- paths
  
  # Initialize result table
  results <- data.frame(
    pathway = paths,
    obs_moran = NA_real_, z_moran = NA_real_,
    p_moran = NA_real_, p_adj_moran = NA_real_,
    null_mean_moran = NA_real_, null_sd_moran = NA_real_,
    obs_geary = NA_real_, z_geary = NA_real_,  
    p_geary = NA_real_, p_adj_geary = NA_real_,
    null_mean_geary = NA_real_, null_sd_geary = NA_real_,
    up = NA_real_, down = NA_real_,
    B = n_perm, cv_clip = weight_clip,
    global_sigma = sigma_global, optimal_sigma = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Loop over each pathway
  for (i in seq_len(n_pathways)) {
    vec <- scores[i, ]
    pth <- paths[i]
    
    ## Moran's I
    mc_moran <- spdep::moran.mc(vec, w, nsim = n_perm, alternative = "greater")
    obs_moran <- as.numeric(mc_moran$statistic)
    sim_moran <- as.numeric(mc_moran$res)
    p_moran <- (sum(sim_moran >= obs_moran) + 1) / (n_perm + 1)
    
    results$obs_moran[i] <- obs_moran
    results$null_mean_moran[i] <- mean(sim_moran)
    results$null_sd_moran[i] <- sd(sim_moran)
    results$z_moran[i] <- if (sd(sim_moran) > 0) {
      (obs_moran - mean(sim_moran)) / sd(sim_moran)
    } else NA_real_
    results$p_moran[i] <- p_moran
    
    ## Geary's C
    mc_geary <- spdep::geary.mc(vec, w, nsim = n_perm, alternative = "less")
    obs_geary <- as.numeric(mc_geary$statistic)
    sim_geary <- as.numeric(mc_geary$res)
    p_geary <- (sum(sim_geary <= obs_geary) + 1) / (n_perm + 1)
    
    results$obs_geary[i] <- obs_geary
    results$null_mean_geary[i] <- mean(sim_geary)
    results$null_sd_geary[i] <- sd(sim_geary)
    results$z_geary[i] <- if (sd(sim_geary) > 0) {
      (obs_geary - mean(sim_geary)) / sd(sim_geary)
    } else NA_real_
    results$p_geary[i] <- p_geary
    
    results$up[i] <- length(sets[[pth]]$U)
    results$down[i] <- length(sets[[pth]]$D)
    results$optimal_sigma[i] <- pathway_sigmas[[pth]]$median_sigma
    
    null_hist[[pth]] <- list(
      perm_moran = sim_moran, perm_geary = sim_geary
    )
  }
  
  # Adjust p-values (FDR)
  pval_moran <- p.adjust(results$p_moran, method = "BH")
  pval_geary <- p.adjust(results$p_geary, method = "BH")
  results$p_adj_moran <- ifelse(pval_moran < 1/(B+1), paste0("p < ", formatC(1/(B+1), format = "e", digits = 1)), pval_moran)
  results$p_adj_geary <- ifelse(pval_geary < 1/(B+1), paste0("p < ", formatC(1/(B+1), format = "e", digits = 1)), pval_geary)
  
  return(list(
    res = results, perm_null_res = null_hist
  ))
}

# --- implementation of early stopping ---
# Clopper-Pearson bounds
clopper_bounds <- function(k, t, alpha = 0.05) {
  list(
    lower = qbeta(alpha / 2, k + 1, t - k + 1),
    upper = qbeta(1 - alpha / 2, k + 1, t - k + 1)
  )
}
# permute within blocks (identity for non-block)
permute_within_blocks <- function(x, blocks) {
  if (is.null(blocks)) {
    return(sample(x))
  } else {
    perm <- x
    for (b in unique(blocks)) {
      idx <- which(blocks==b)
      if (length(idx)>1) { perm[idx] <- sample(perm[idx]) }
    }
    return(perm)
  }
}
# general permutation test with early stopping
perm_with_early_stop <- function(stat_fun, x, listw, blocks = NULL, 
                                 t_min = 1000, t_max = 10000, alpha = 0.05) {
  alpha_decisive <- alpha/10
  obs <- stat_fun(x, listw)
  b <- 0; t <- 0; stats <- numeric(t_max)
  stop_reason <- NA_real_
  
  while(t < t_max) {
    t <- t+1
    x_perm <- permute_within_blocks(x, blocks)
    s_perm <- stat_fun(x_perm, listw)
    stats[t] <- s_perm
    b <- b + as.integer(s_perm >= obs) # for Geary's C, use wrapper to flip sign
    
    if (t >= t_min) {
      cb <- clopper_bounds(b, t, alpha)
      if (cb$lower > alpha) {
        stop_reason <- "non_significant" 
        break} 
      if (cb$upper < alpha_decisive) {
        stop_reason <- "significant" 
        break}
    }
  }
  
  p_hat <- (b+1) / (t+1)
  mu <- mean(stats[1:t]); sd <- sd(stats[1:t])
  z <- if (sd > 0) (obs - mu) / sd else NA_real_
  
  return(list(
    obs = obs, p = p_hat,
    t = t, b = b, 
    stop_reason = ifelse (t==t_max, "reach the maximum", stop_reason),
    null_mean = mu, null_sd = sd,
    z = z, perm_stats = stats[1:t]
  ))
}

# Moran’s I wrapper (right-tailed)
stat_moran <- function(x, listw) moran.test(x, listw)$estimate[1]
# Geary’s C wrapper (left-tailed → flip sign to reuse logic)
stat_geary <- function(x, listw) -geary.test(x, listw)$estimate[1]

# Moran with hex-block
stat_moran_hex <- function(x, w) moran.test(x, w)$estimate[1]
# Geary with hex-block (flipped)
stat_geary_hex <- function(x, w) -geary.test(x, w)$estimate[1]
