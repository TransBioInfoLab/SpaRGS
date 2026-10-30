# evaluation metrics functions

# metrics : AUROC recovery of top -q true regions
top_mask <- function(v, q=0.1) as.integer(v >= quantile(v, probs=1-q))

# Explicit AUROC calculation
score_auc_fixed <- function(truth, score) {
  
  keep <- is.finite(truth) & is.finite(score)
  if (length(unique(truth[keep])) < 2) {
    return(NA_real_)
  }
  
  # Fixed orientation: higher scores are expected to indicate active regions.
  roc_obj <- pROC::roc(
    response = truth[keep],
    predictor = score[keep],
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  as.numeric(pROC::auc(roc_obj))
}

score_auc_invariant <- function(truth, score) {
  auc_fixed <- score_auc_fixed(truth, score)
  
  if (!is.finite(auc_fixed)) {
    return(NA_real_)
  }
  max(auc_fixed, 1 - auc_fixed)
}

# make AUROC and spearman metrics for benchmarking with oracle directionality
eval_metrics <- function(score, method, scn, sim, rt, q=0.1) {
  truth_paths <- names(sim$fields)
  temp <- data.frame(
    pathway = truth_paths, methods = method,
    AUROC = NA_real_, AUROC_fixed = NA_real_,
    spearman = NA_real_, spearman_signed = NA_real_, 
    retained = FALSE, times = unname(rt), 
    ell = scn$ell, sigma_f = scn$sigma_f,
    up_size = scn$up_size, down_size = scn$down_size,
    stringsAsFactors = FALSE
  )
  
  # Handle case where one method has no valid scores
  no_score <- is.null(score) || (length(score) == 1 && identical(as.numeric(score), 0))
  if (no_score) {
    return(temp)
  }
  
  scored_paths <- intersect(truth_paths, rownames(score))
  temp$retained[temp$pathway %in% scored_paths] <- TRUE
  
  for (p in scored_paths) {
    idx <- temp$pathway == p
    field <- sim$fields[[p]]
    truth <- top_mask(field, q=q)
    
    auc_fixed <- score_auc_fixed(truth, score[p,])
    auc_invariant <- score_auc_invariant(truth, score[p,])
    temp$AUROC[idx] <- auc_invariant
    temp$AUROC_fixed[idx] <- auc_fixed
    
    rho <- cor(field, score[p, ], method = "spearman", use = "pairwise.complete.obs")
    temp$spearman_signed[idx] <- rho          # keep for auditing sign behaviour
    temp$spearman[idx]        <- abs(rho)      # report this: pattern recovery
  }
  return(temp)
}

# make AUROC and spearman metrics for ablation with data-driven directionality
eval_metrics_abl <- function(score, sets, method, scn, sim, rt, q=0.1) {
  if (length(sets)==0) {
    temp <- data.frame(
      pathway = "no pathways passed filter", 
      methods = method, AUROC = NA_real_, 
      spearman = NA_real_, spearman_signed = NA_real_,
      times = rt, ell = scn$ell, sigma_f = scn$sigma_f,
      up_size = 0, down_size = 0,
      true_up = scn$up_size, true_down = scn$down_size,
      stringsAsFactors = FALSE
    )
    return(temp)
  }
  sel_paths <- rownames(score)
  temp <- data.frame(
    pathway = sel_paths, methods = method,
    AUROC = NA_real_, spearman = NA_real_, 
    spearman_signed = NA_real_, times = rt,
    ell = scn$ell, sigma_f = scn$sigma_f,
    up_size = NA_real_, down_size = NA_real_,
    true_up = scn$up_size, true_down = scn$down_size,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(sel_paths)) {
    p <- sel_paths[i]
    field <- sim$fields[[p]]
    truth <- top_mask(field, q=q)
    temp$AUROC[i] <- score_auc_invariant(truth, score[p,])
    rho <- cor(field, score[p, ], method = "spearman", use = "pairwise")
    temp$spearman[i] <- abs(rho)
    temp$spearman_signed[i] <- rho
    temp$up_size[i] <- length(sets[[p]]$U)
    temp$down_size[i] <- length(sets[[p]]$D)
  }
  temp
}

# aggregate results
read_with_rep <- function(files, pattern_rep = "replicate(\\d+)") {
  if (!length(files)) return(NULL)
  bind_rows(lapply(files, function(f) {
    repx <- suppressWarnings(as.integer(stringr::str_match(basename(f), pattern_rep)[,2]))
    df <- readr::read_csv(f, show_col_types = FALSE)
    if (!"replicate" %in% names(df)) df$replicate <- repx
    df$.src_file <- f
    df
  }))
}

# two-sided Stouffer combiner
stouffer_two_sided <- function(p, effect, w = NULL) {
  stopifnot(length(p) == length(effect))
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  if (is.null(w)) w <- rep(1, length(p))
  z <- qnorm(p/2, lower.tail = FALSE)*sign(effect) # two-sided
  zc <- sum(w * z) / sqrt(sum(w^2))
  2 * pnorm(abs(zc), lower.tail = FALSE)
}

