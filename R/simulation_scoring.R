# pathway methods for simulation


# --- Helper functions (as used above) ---

# Gaussian kernel smoothing
smooth_expr <- function(expr, coords, sigma=15, 
                        mode=c("delaunay","knn"), k=6) {
  mode <- match.arg(mode)
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

signed_score <- function(expr, split_sets, w, ranks, 
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
          base <- median(ranks[, i], na.rm=TRUE)
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

# Apply z-normalization to SpaRGS scores
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


# ---- SpaRGS ----
# Requirements :
# - inputs : counts ( genes x spots ), coords (x,y), gene_sets ( list )
# - returns : score matrix ( pathways x spots )
run_spargs <- function(expr, coords, sets, 
                      # smoothing/weighting
                      no_smoothing = FALSE, # remove smoothing
                      mode = c("delaunay","knn"), k = 6, # method for getting neighborhood
                      sigma_global = 150, # global smoothing bandwidth
                      floor_frac = 0.05, # min non-zero fraction for gene filtering
                      weight_clip = 0.99, # CV weight clipping quantile
                      # correlation direactionality
                      directionality = c("oracle", "data-driven"),
                      cor_min = 0.10, # abs(Spearman) threshold for directionality
                      min_genes = 15,  # minimum U+D
                      max_n = 25, # maximum balance cap
                      one_sided = c("off", "center_0.5", "center_median") # one-sided options
) {
  directionality <- match.arg(directionality)
  one_sided <- match.arg(one_sided)
  mode <- match.arg(mode)
  
  # Gene-set coverage relative to the supplied pathway definitions
  n_input <- vapply(sets, function(gs) length(unique(gs)), integer(1))
  n_available <- vapply(
    sets, function(gs) length(intersect(unique(gs), rownames(expr))), integer(1)
  )
  
  pathway_coverage <- data.frame(
    pathway = names(sets),
    n_genes_input = unname(n_input),
    n_genes_available = unname(n_available),
    coverage_fraction = ifelse(n_input > 0, n_available / n_input, NA_real_),
    row.names = NULL
  )
  
  # --- Smoothing ---
  if (no_smoothing == FALSE) {
    smoothed <- smooth_expr(expr, coords, sigma = sigma_global, mode=mode, k=k)
    smoothed <- as.matrix(smoothed)
    colnames(smoothed) <- colnames(expr)
  } else {smoothed <- expr}
  
  # --- Weights ---
  w <- gene_weights(expr, clip = weight_clip)
  
  # --- Directionality + balancing ---
  if (directionality == "data-driven") {
    corr_res <- derive_direction_correlation(expr, sets,
                                             min_cor = cor_min, min_genes = min_genes)
    corr_sets <- corr_res$sets # drop pathways that failed min_genes
    sets <- lapply(corr_sets, function(s) balance_sets(s$U, s$D, w, max_n = max_n))
  } else if (directionality == "oracle") {
    sets <- sets
  }
  
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
      weights = w,
      pathway_coverage = pathway_coverage
    ))
  }
  
  # Final number of selected genes after directional filtering and caps
  pathway_coverage$n_selected <- vapply(
    pathway_coverage$pathway,
    function(p) {if (!p %in% names(sets)) return(0L) 
      length(sets[[p]]$U) + length(sets[[p]]$D)
    }, integer(1)
  )
  
  # --- Rank scaling ---
  ranks <- apply(smoothed, 2, function(x) rank(x, ties.method="average")/length(x))
  rownames(ranks) <- rownames(smoothed)
  
  # --- scoring + z-normalization ---
  scores <- signed_score(smoothed, sets, w, ranks, one_sided = one_sided)
  z <- z_normalization(as.matrix(scores))
  
  return(list(
    scores = scores, 
    z_score = z,
    gene_sets = sets,
    weights = w,
    pathway_coverage = pathway_coverage
  ))
}

# ---- GSVA ----
run_gsva <- function(expr, sets) {
  params <- gsvaParam(
    expr, sets
  )
  gsva_scores <- gsva(params, verbose = FALSE)
  return(gsva_scores)
}

# ---- AUCell ----
run_aucell <- function(expr, sets) {
  rk <- AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
  AUCell::getAUC(AUCell_calcAUC(sets, rk))
}

# ---- UCell ----
run_ucell <- function(expr, sets) {
  score <- UCell::ScoreSignatures_UCell(as.matrix(expr), features = sets)
  ucell_scores <- t(score)
  rownames(ucell_scores) <- names(sets)
  ucell_scores
}

# ---- ssGSEA ----
run_ssgsea <- function(expr, sets) {
  params <- ssgseaParam(
    expr, sets
  )
  gsva_scores <- gsva(params, verbose = FALSE)
  return(gsva_scores)
}

# ---- AddModuleScore ----
run_addmod <- function(counts, sets, nbin = 12, ctrl = 10) {
  # Convert to Seurat object temporarily
  expr_full <- as(as.matrix(counts), "dgCMatrix")
  seu <- Seurat::CreateSeuratObject(counts = expr_full)
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e6)
  seu <- Seurat::AddModuleScore(object=seu, features = sets, 
                                name = "pathway", nbin = nbin, ctrl = ctrl)
  # Extract matrix of scores
  addmod_scores <- as.matrix(seu@meta.data[, grep("pathway", colnames(seu@meta.data))])
  rownames(addmod_scores) <- colnames(expr_full)  # ensure sample alignment
  colnames(addmod_scores) <- names(sets)
  addmod_scores <- t(addmod_scores)
  return(addmod_scores)
}