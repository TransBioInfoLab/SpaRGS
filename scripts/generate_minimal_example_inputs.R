#!/usr/bin/env Rscript

# ============================================================
# ONE-TIME helper to create the frozen input used by
# scripts/run_minimal_example.R.
#
# Users performing the clean-install test should NOT run this script.
# ============================================================

set.seed(1)

required_packages <- c("Matrix", "MASS", "msigdbr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

make_hex_grid_small <- function(nx = 12, ny = 12, pitch = 50) {
  xs <- rep(seq(0, by = pitch, length.out = nx), each = ny)
  ys <- rep(seq(0, by = pitch * sqrt(3) / 2, length.out = ny), nx)
  
  # Preserve the offset convention used in R/simulation.R.
  xs_off <- rep(c(0, pitch / 2), length.out = nx)
  ys <- ys + rep(xs_off, each = ny)
  
  data.frame(x = xs, y = ys)
}

simulate_spatial_field <- function(coords, ell = 200, sigma_f = 1.5) {
  D <- as.matrix(stats::dist(coords))
  K <- sigma_f^2 * exp(-D / ell) + 1e-6 * diag(nrow(D))
  
  as.numeric(
    MASS::mvrnorm(
      n = 1,
      mu = rep(0, nrow(D)),
      Sigma = K
    )
  )
}

simulate_counts_small <- function(
    coords,
    gene_sets,
    sel_paths,
    up = 25,
    down = 25,
    ell = 200,
    sigma_f = 1.5,
    mean_lib = 8000,
    sd_lib = 0.6,
    seed = 1) {
  
  set.seed(seed)
  
  N <- nrow(coords)
  G <- sort(unique(unlist(gene_sets[sel_paths])))
  
  if (length(G) == 0L) {
    stop("No genes were found for the selected KEGG pathways.")
  }
  
  # Require enough unique genes for the balanced example.
  n_unique <- vapply(
    gene_sets[sel_paths],
    function(gs) length(unique(intersect(gs, G))),
    integer(1)
  )
  
  if (any(n_unique < (up + down))) {
    bad <- names(n_unique)[n_unique < (up + down)]
    stop(
      "The following example pathway(s) contain fewer than ",
      up + down,
      " unique genes in the current msigdbr release: ",
      paste(bad, collapse = ", "),
      ".\nChoose pathways with at least ",
      up + down,
      " genes or reduce up/down for the smoke test."
    )
  }
  
  # Select simulated up/down genes within each pathway.
  per_path <- lapply(sel_paths, function(p) {
    genes <- unique(intersect(gene_sets[[p]], G))
    U <- sample(genes, up, replace = FALSE)
    D <- sample(setdiff(genes, U), down, replace = FALSE)
    list(U = U, D = D)
  })
  names(per_path) <- sel_paths
  
  # One latent spatial field per pathway.
  fields <- lapply(
    sel_paths,
    function(p) simulate_spatial_field(coords, ell = ell, sigma_f = sigma_f)
  )
  names(fields) <- sel_paths
  
  # Baseline relative means and NB dispersions.
  mu_g <- stats::rlnorm(length(G), meanlog = log(0.02), sdlog = 1.0)
  names(mu_g) <- G
  
  theta_g <- stats::runif(length(G), min = 0.5, max = 2.0)
  names(theta_g) <- G
  
  beta_by_path <- lapply(sel_paths, function(p) {
    b <- setNames(rep(0, length(G)), G)
    
    U <- per_path[[p]]$U
    D <- per_path[[p]]$D
    
    b[U] <- stats::rlnorm(length(U), -1.2, 0.5)
    b[D] <- -stats::rlnorm(length(D), -1.2, 0.5)
    
    b
  })
  names(beta_by_path) <- sel_paths
  
  L <- round(
    stats::rlnorm(N, meanlog = log(mean_lib), sdlog = sd_lib)
  )
  
  gx <- scale(coords$x)
  gy <- scale(coords$y)
  gamma_i <- as.numeric(0.2 * gx - 0.2 * gy)
  
  # Construct mean matrix.
  log_mu <- outer(rep(1, N), log(mu_g))
  
  for (p in sel_paths) {
    b <- beta_by_path[[p]]
    genes <- names(b)[b != 0]
    
    if (length(genes) > 0L) {
      log_mu[, genes] <-
        log_mu[, genes, drop = FALSE] +
        outer(fields[[p]], b[genes])
    }
  }
  
  log_mu <- log_mu + gamma_i
  log_mu <- pmin(pmax(log_mu, -20), 20)
  
  MU <- exp(log_mu) * L
  
  X <- matrix(
    0,
    nrow = length(G),
    ncol = N,
    dimnames = list(
      G,
      paste0("spot", seq_len(N))
    )
  )
  
  for (j in seq_along(G)) {
    X[j, ] <- stats::rnbinom(N, mu = MU[, j], size = theta_g[j])
  }
  
  list(
    counts = Matrix::Matrix(X, sparse = TRUE),
    fields = fields,
    per_path = per_path,
    coords = coords
  )
}

msig <- tryCatch(
  msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:KEGG_LEGACY"
  ),
  error = function(e) {
    msigdbr::msigdbr(
      species = "Homo sapiens",
      category = "C2",
      subcollection = "CP:KEGG_LEGACY"
    )
  }
)

gene_sets <- split(msig$gene_symbol, msig$gs_name)
gene_sets <- lapply(gene_sets, unique)

# Representative immune, signaling, and metabolic/neuronal programs
# used elsewhere in the sGSVA project.
sel_paths <- c(
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_OXIDATIVE_PHOSPHORYLATION"
)

missing_paths <- setdiff(sel_paths, names(gene_sets))

if (length(missing_paths) > 0L) {
  stop(
    "The following KEGG Legacy pathways were not found in the ",
    "installed msigdbr release: ",
    paste(missing_paths, collapse = ", ")
  )
}

coords <- make_hex_grid_small(nx = 12, ny = 12, pitch = 50)
rownames(coords) <- paste0("spot", seq_len(nrow(coords)))

# Use one strong-signal balanced setting from the manuscript
# simulation design. This is chosen for a stable smoke test.
ell <- 200
sigma_f <- 1.5
up <- 25
down <- 25
mean_lib <- 8000
sd_lib <- 0.6
sigma_global <- floor(ell / sqrt(2))

sim <- simulate_counts_small(
  coords = coords,
  gene_sets = gene_sets,
  sel_paths = sel_paths,
  up = up,
  down = down,
  ell = ell,
  sigma_f = sigma_f,
  mean_lib = mean_lib,
  sd_lib = sd_lib,
  seed = 1
)

example_input <- list(
  counts = sim$counts,
  coords = sim$coords,
  gene_sets = gene_sets[sel_paths],
  fields = sim$fields,
  per_path = sim$per_path,
  selected_pathways = sel_paths,
  parameters = list(
    ell = ell,
    sigma_f = sigma_f,
    up = up,
    down = down,
    mean_lib = mean_lib,
    sd_lib = sd_lib,
    simulation_seed = 1,
    sigma_global = sigma_global,
    msigdbr_version = as.character(utils::packageVersion("msigdbr")),
    R_version = R.version.string
  )
)

out_dir <- file.path("data", "example")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "minimal_example_input.rds")

saveRDS(example_input, out_file, version = 3)

cat("Frozen minimal-example input written to:\n  ", out_file, "\n", sep = "")
cat("MD5: ", unname(tools::md5sum(out_file)), "\n", sep = "")
cat("Genes: ", nrow(sim$counts), "\n", sep = "")
cat("Spots: ", ncol(sim$counts), "\n", sep = "")
cat("Pathways: ", paste(sel_paths, collapse = "; "), "\n", sep = "")



