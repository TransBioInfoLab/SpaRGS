# R/sgsva_score.R

# set sigma proportional to the median inter-neighbor distance
estimate_sigma_knn <- function(coords, k = 6, scale_factor = 1.5) {
  knn <- FNN::get.knn(as.matrix(coords), k = k)
  # knn$nn.dist is an n × k matrix of neighbor distances
  d_med <- median(knn$nn.dist, na.rm = TRUE)
  scale_factor * d_med
}

# --- set sGSVA function ---
run_sgsva <- function(expr, coords, sets,
                      # smoothing/weighting
                      no_smoothing = FALSE, # remove smoothing
                      mode = c("delaunay","knn"), k=6, # mode for neighborhood
                      sigma_global = 200, # global smoothing bandwidth (um in lowres coords)
                      weight_clip = 0.99, # CV weight clipping quantile
                      cor_min = 0.10, # abs(Spearman) threshold for directionality
                      min_genes = 15,  # minimum U+D
                      max_n = 25, # maximum balance cap
                      one_sided = c("off", "center_0.5", "center_median") # one-sided options
) {
  # revised version from previous run_sgsva function based on seurat object, Visium data
  # minimal stand -in: Gaussian - smooth logCPM then rank - based scoring (U minus D)
  mode <- match.arg(mode)
  one_sided <- match.arg(one_sided)
  stopifnot(ncol(expr) == nrow(coords))
  
  # --- Smoothing ---
  if (no_smoothing == FALSE) {
    smoothed <- smooth_expr(expr, coords, sigma = sigma_global, mode=mode, k=k)
    smoothed <- as.matrix(smoothed)
    colnames(smoothed) <- colnames(expr)
  } else {smoothed <- expr}
  
  # --- Weights ---
  w <- gene_weights(expr, clip = weight_clip)
  
  # --- Directionality + balancing ---
  corr_res <- derive_direction_correlation(expr, sets,
                                           min_cor = cor_min, min_genes = min_genes)
  corr_sets <- corr_res$sets # drop pathways that failed min_genes
  sets <- lapply(corr_sets, function(s) balance_sets(s$U, s$D, w, max_n = max_n))
  
  if (length(sets) == 0) {
    z_score <- matrix(
      NA_real_,
      nrow = 0,
      ncol = ncol(expr),
      dimnames = list(character(0), colnames(expr))
    )
    return(list(
      z_score = z_score,
      gene_sets = sets,
      directionality_diagnostics = corr_res$diagnostics
    ))
  }
  
  # --- Rank scaling ---
  ranks <- apply(smoothed, 2, function(x) rank(x, ties.method="average")/length(x))
  rownames(ranks) <- rownames(smoothed)
  
  # --- sGSVA scoring + z-normalization ---
  scores <- signed_gsva(smoothed, sets, w, ranks, one_sided=one_sided)
  z <- z_normalization(as.matrix(scores))
  
  return(list(
    scores = scores, 
    z_score = z,
    gene_sets = sets,
    directionality_diagnostics = corr_res$diagnostics
  ))
}


# --- helper functions ---
# Gaussian kernel smoothing
smooth_expr <- function(expr, coords, sigma=15,
                        mode = c("delaunay", "knn"), k=6) {
  n <- nrow(coords)
  coords <- as.matrix(coords)
  if (mode=="delaunay") {
    nb <- spdep::tri2nb(coords)
  } else if (mode=="knn") {
    nb <- spdep::knn2nb(spdep::knearneigh(coords, k=k))
  }
  
  # compute Gaussian kernel weights as sparse triplets
  I <- integer(0); J <- integer(0); VAL <- numeric(0)
  for (i in seq_len(n)) {
    ngb <- nb[[i]]
    if (length(ngb)==0L) {next}
    
    dist_sq <- rowSums((coords[ngb, , drop = FALSE] - matrix(coords[i, ], nrow=length(ngb), ncol=2, byrow=TRUE))^2)
    w <- exp(-dist_sq / (2*sigma^2))
    I <- c(I, rep(i, length(ngb)))
    J <- c(J, ngb)
    VAL <- c(VAL, w)
  }
  
  # build sparse weight W
  W <- Matrix::sparseMatrix(i=I, j=J, x=VAL, dims=c(n,n))
  rs <- Matrix::rowSums(W); W <- Matrix::Diagonal(n, 1/pmax(rs, 1e-12)) %*% W
  expr <- as(expr, "dgCMatrix")
  gene_names <- rownames(expr)
  spot_names <- colnames(expr)
  
  expr_smooth <- expr %*% Matrix::t(W)
  rownames(expr_smooth) <- gene_names
  colnames(expr_smooth) <- spot_names
  return(expr_smooth)
}

# gene weights (coefficient of variation)
gene_weights <- function(expr, clip=0.99) {
  mu <- rowMeans(expr)
  sd <- matrixStats::rowSds(expr)
  cv <- sd / pmax(mu, 0.1)
  w <- pmin(cv, stats::quantile(cv, probs=clip, na.rm=TRUE))
  names(w) <- rownames(expr); w
}

# apply correlation-based directionality
derive_direction_correlation <- function(expr, gene_sets, min_cor = 0.1, min_genes = 15) {
  kept <- list()
  diag <- vector("list", length(gene_sets))
  names(diag) <- names(gene_sets)
  
  for (nm in names(gene_sets)){
    genes <- intersect(gene_sets[[nm]], rownames(expr))
    out <- data.frame(
      pathway = nm,
      n_genes_input = length(genes),
      n_up = 0L, n_down = 0L,
      n_selected = 0L, retained = FALSE,
      pc1_anchor_cor = NA_real_,
      min_cor = min_cor, min_genes = min_genes,
      stringsAsFactors = FALSE
    )
    
    if (length(genes)<2){ 
      diag[[nm]] <- out
      next 
    }
    X <- t(expr[genes, , drop=FALSE])
    pc <- prcomp(X, center=TRUE, scale.=FALSE)
    eig <- pc$x[,1]
    
    ## --- SIGN ANCHORING (new) -------------------------------------------------
    anchor <- rowMeans(X)                       # mean pathway expression / location
    a_cor  <- suppressWarnings(cor(eig, anchor, method = "spearman"))
    if (!is.na(a_cor) && a_cor < 0) {
      eig <- -eig
    }
    
    cors <- apply(X, 2, function(z) suppressWarnings(cor(z, eig, method="spearman")))
    pos <- names(cors)[cors >= min_cor]
    neg <- names(cors)[cors <= -min_cor]
    U <- intersect(pos, rownames(expr))
    D <- intersect(neg, rownames(expr))
    
    out$n_up <- length(U)
    out$n_down <- length(D)
    out$n_selected <- length(U) + length(D)
    out$pc1_anchor_cor <- a_cor
    out$retained <- out$n_selected >= min_genes
    diag[[nm]] <- out
    
    if (length(U) + length(D) < min_genes) next
    kept[[nm]] <- list(U=U, D=D)
  }
  list(sets=kept, diagnostics = dplyr::bind_rows(diag))
}

# Balancing helper function
balance_sets <- function(U, D, w, max_n = 25) {
  U <- intersect(U, names(w)); D <- intersect(D, names(w))
  if (length(U)>0) U <- head(U[order(-w[U])], max_n)
  if (length(D)>0) D <- head(D[order(-w[D])], max_n)
  list(U = U, D = D)
}

# calculate signed, weighted GSVA
signed_gsva <- function(expr, split_sets, w, ranks, 
                        one_sided = c("off", "center_0.5", "center_median")){
  n_p <- length(split_sets); n_s <- ncol(expr)
  out <- matrix(0, nrow=n_p, ncol=n_s,
                dimnames=list(names(split_sets), colnames(expr)))
  
  for (p in seq_len(n_p)){
    nm <- names(split_sets)[p]
    U <- intersect(split_sets[[nm]]$U, rownames(expr))
    D <- intersect(split_sets[[nm]]$D, rownames(expr))
    if (length(D)==0 || length(U)==0){
      genes <- union(U, D)
      for (i in 1:n_s){
        R <- ranks[genes,i]; ww <- w[genes]
        sc <- sum(ww*R, na.rm=TRUE)/sum(abs(ww), na.rm=TRUE)
        if (one_sided == "center_0.5") {
          out[p,i] <- sum(ww*(R-0.5), na.rm=TRUE)/sum(abs(ww), na.rm=TRUE)
        } else if (one_sided == "center_median") {
          base <- median(ranks[genes, i], na.rm=TRUE)
          out[p,i] <- sc-base
        } else if (one_sided == "off") {
          out[p,i] <- sc
        }
      }
    } else {
      for (i in 1:n_s){
        Rp <- ranks[U,i]; Rm <- ranks[D,i]
        wp <- w[U]; wm <- w[D]
        num <- sum(wp*Rp, na.rm=TRUE) - sum(wm*Rm, na.rm=TRUE)
        den <- sum(abs(wp), na.rm=TRUE) + sum(abs(wm), na.rm=TRUE)
        out[p,i] <- ifelse(den>0, num/den, 0)
      }
    }
  }
  out
}

# Apply z-normalization to sGSVA scores
z_normalization <- function(scores) {
  n_pathways <- nrow(scores)
  n_spots <- ncol(scores)
  z_scores <- matrix(0, nrow = n_pathways, ncol = n_spots)
  rownames(z_scores) <- rownames(scores)
  colnames(z_scores) <- colnames(scores)
  
  for (p in 1:n_pathways) {
    mean_p <- mean(scores[p, ], na.rm = TRUE)
    sd_p <- sd(scores[p, ], na.rm = TRUE)
    if (sd_p > 0) {
      z_scores[p, ] <- (scores[p, ] - mean_p) / sd_p
    } else {
      z_scores[p, ] <- scores[p, ]  # Retain original scores if SD is zero
    }
  }
  return(z_scores)
}
