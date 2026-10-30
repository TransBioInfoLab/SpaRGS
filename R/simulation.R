# R/simulation.R

# ---- geometry : use a realistic layout ( import from a real slide ) or make a grid ----
make_hex_grid <- function(nx=70, ny=70, pitch=50){
  # regular hex-like layout in microns
  xs <- rep(seq(0, by=pitch, length.out=nx), each=ny)
  ys <- rep(seq(0, by=pitch*sqrt(3)/2, length.out=ny), nx)
  # offset for even columns
  xs_off <- rep(c(0, pitch/2), length.out=nx)
  ys <- ys + rep(xs_off, each=ny)
  data.frame(x=xs, y=ys)
}

# ---- simulate latent spatial field via Gaussian kernel ----
rfield <- function(coords, ell=200, sigma_f=1.0) {
  D <- as.matrix(dist(coords))
  K <- sigma_f^2 * exp(-D/ell) + 1e-6*diag(nrow(D))
  as.numeric(mvrnorm(1, mu=rep(0, nrow(D)), Sigma=K))
}

# ---- simulate counts with pathway - driven loadings ----
simulate_counts <- function(coords, gene_sets, sel_paths, 
                            up=25, down=25, ell=200, sigma_f=1.0,
                            mean_lib=8000, sd_lib=0.6, cap=5000) {
  N <- nrow(coords)
  # make gene universe
  G <- sort(unique(unlist(gene_sets[sel_paths])))
  G <- unique(G)[1:min(cap, length(unique(G)))] # cap for seed
  # choose up/down subsets within each pathway
  pd <- lapply(sel_paths, function(p) {
    genes <- intersect(gene_sets[[p]], G)
    # No genes from this pathway remain in the simulated universe
    if (length(genes) == 0) {
      return(list(U = character(0), D = character(0)))
    }
    if (length(genes)<(up + down)) genes <- rep(genes, length.out=(up + down))
    U <- sample(genes, up, replace=FALSE)
    avail <- setdiff(genes, U)
    if (down > 0) {
      if (length(avail) == 0) {
        D <- character(0)
      } else {
        if (length(avail) < down) {
          avail <- rep(avail, length.out = down)
        }
        D <- sample(avail, down, replace = FALSE)
      }
    } else {D <- character(0)}
    
    list(U=U, D=D)
  })
  names(pd) <- sel_paths
  
  # latent fields per pathway
  fields <- lapply(sel_paths, function(p) rfield(coords, ell=ell, sigma_f=sigma_f))
  names(fields) <- sel_paths
  
  # per-gene baselines and dispersions
  mu_g <- rlnorm(length(G), meanlog=log(0.02), sdlog=1.0) # baseline relative means
  names(mu_g) <- G
  theta_g <- runif(length(G), min=0.5, max=2.0); names(theta_g) <- G
  
  # gene loadings; positive for up, negative for down
  beta <- rep(0, length(G)); names(beta) <- G
  for (p in sel_paths) {
    U <- intersect(pd[[p]]$U, G)
    D <- intersect(pd[[p]]$D, G)
    beta[U] <- beta[U] + rlnorm(length(U), -1.2, 0.5)
    beta[D] <- beta[D] - rlnorm(length(D), -1.2, 0.5)
  }
  
  # library sizes
  L <- round(rlnorm(nrow(coords), log(mean_lib), sd_lib))
  # small spatial library - size gradient ( optional confounder )
  gx <- scale(coords$x); gy <- scale(coords$y)
  gamma_i <- as.numeric(0.2*gx-0.2*gy)
  
  # construct mean matrix (N x G) and sample NB counts
  # log mu_{ig} = log (mu_g) + (sum_p beta_g^{(p)} * f_i^{( p)}) + gamma _i
  log_mu <- outer(rep(1,N), log(mu_g))
  for (p in sel_paths) {
    U <- pd[[p]]$U
    D <- pd[[p]]$D
    genes <- intersect(c(U, D), names(mu_g))
    if (length(genes) > 0) {
      log_mu[, genes] <- log_mu[, genes] + outer(fields[[p]], beta[genes])
    }
  }
  log_mu <- log_mu + gamma_i
  # clip for stability
  log_mu <- pmin(pmax(log_mu, -20), 20)
  MU <- exp(log_mu)
  MU <- MU*L # scale by library size
  
  X <- matrix(0, nrow=length(G), ncol=N, dimnames=list(G, paste0("spot", seq_len(N))))
  for (j in seq_along(G)) {
    mu <- MU[,j]
    th <- theta_g[j]
    X[j,] <- rnbinom(N, mu=mu, size=th)
  }
  list(counts = Matrix(X, sparse = TRUE),
       fields = fields,
       per_path = pd,
       coords = coords)
}

