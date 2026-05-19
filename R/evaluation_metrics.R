# evaluation metrics functions

# metrics : AUROC recovery of top -q true regions
top_mask <- function(v, q=0.1) as.integer(v >= quantile(v, probs=1-q))
score_auc <- function(truth, score) {
  roc <- pROC::roc(truth, as.numeric(score), quiet=TRUE)
  as.numeric(pROC::auc(roc))
}

# DeLong tests
delong_per_pathway <- function(fields, scores, methods = c("sGSVA","GSVA"), q = 0.10, scn) {
  stopifnot(all(methods %in% names(scores)))
  paths <- names(fields)
  
  bind_rows(lapply(paths, function(p) {
    y <- as.integer(fields[[p]] >= quantile(fields[[p]], 1 - q))  # labels
    s1 <- as.numeric(scores[[methods[1]]][p,])
    s2 <- as.numeric(scores[[methods[2]]][p,])
    
    r1 <- pROC::roc(y, s1, quiet = TRUE)
    r2 <- pROC::roc(y, s2, quiet = TRUE)
    tt <- pROC::roc.test(r1, r2, method = "delong", paired = TRUE)
    
    data.frame(
      pathway  = p,
      method1  = methods[1], method2  = methods[2],
      auc1 = as.numeric(pROC::auc(r1)), auc2 = as.numeric(pROC::auc(r2)),
      auc_diff = as.numeric(pROC::auc(r1) - pROC::auc(r2)),
      p_delong = tt$p.value,
      ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
      ell = scn$ell, sigma_f = scn$sigma_f,
      up_size = scn$up_size, down_size = scn$down_size)
  }))
}

# make AUROC and spearman metrics for benchmarking with oracle directionality
eval_metrics <- function(score, method, scn, sim, rt, q=0.1) {
  sel_paths <- rownames(score)
  temp <- data.frame(
    pathway = sel_paths, methods = rep(method, length(sel_paths)),
    AUROC = NA_real_, spearman = NA_real_,
    times = rep(rt, length(sel_paths)),
    ell = rep(scn$ell, length(sel_paths)), 
    sigma_f = rep(scn$sigma_f, length(sel_paths)),
    up_size = rep(scn$up_size, length(sel_paths)), 
    down_size = rep(scn$down_size, length(sel_paths)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(sel_paths)) {
    p <- sel_paths[i]
    field <- sim$fields[[p]]
    truth <- top_mask(field, q=q)
    temp$AUROC[i] <- score_auc(truth, score[p,])
    temp$spearman[i] <- cor(field, score[p,], method = "spearman", use = "pairwise")
  }
  temp
}

# make AUROC and spearman metrics for ablation with data-driven directionality
eval_metrics_abl <- function(score, sets, method, scn, sim, rt, q=0.1) {
  if (length(sets)==0) {
    temp <- data.frame(
      pathway = "no pathways passed filter", 
      methods = method, AUROC = NA_real_, spearman = NA_real_,
      times = rt, ell = scn$ell, sigma_f = scn$sigma_f,
      up_size = 0, down_size = 0,
      true_up = scn$up_size, true_down = scn$down_size,
      stringsAsFactors = FALSE
    )
    return(temp)
  }
  sel_paths <- rownames(score)
  temp <- data.frame(
    pathway = sel_paths, methods = rep(method, length(sel_paths)),
    AUROC = NA_real_, spearman = NA_real_,
    times = rep(rt, length(sel_paths)),
    ell = rep(scn$ell, length(sel_paths)), 
    sigma_f = rep(scn$sigma_f, length(sel_paths)),
    up_size = NA_real_, down_size = NA_real_,
    true_up = rep(scn$up_size, length(sel_paths)), 
    true_down = rep(scn$down_size, length(sel_paths)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(sel_paths)) {
    p <- sel_paths[i]
    field <- sim$fields[[p]]
    truth <- top_mask(field, q=q)
    temp$AUROC[i] <- score_auc(truth, score[p,])
    temp$spearman[i] <- cor(field, score[p,], method = "spearman", use = "pairwise")
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

