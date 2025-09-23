# qtomography_full.R
# All-in-one: general Born sampling, counts, log-likelihood, gradient,
# projection to density via eigen + simplex, and PGD-MLE solver.
# 所有功能合一：Born 抽样、计数、似然/梯度、单纯形谱投影到密度矩阵、PGD-MLE。

# -------------------------
# Utilities
# -------------------------
Trace <- function(A) sum(diag(A))

hs_inner <- function(A, B) {
  # Hilbert–Schmidt inner product <A,B> = Re tr(A^H B)
  Re(Trace(Conj(t(A)) %*% B))
}

is_density_matrix <- function(rho, atol=1e-10) {
  rhoH <- (rho + Conj(t(rho))) / 2
  ev   <- eigen(rhoH, symmetric = TRUE, only.values = TRUE)$values
  herm <- max(Mod(rho - rhoH)) < atol
  psd  <- min(Re(ev)) >= -1e-12
  tr1  <- Mod(Trace(rho) - 1) < 1e-10
  herm && psd && tr1
}

# -------------------------
# Simplex projection for eigenvalues  (Duchi et al., 2008)
# -------------------------
proj_simplex <- function(d) {
  d <- as.numeric(d); n <- length(d)
  u <- sort(d, decreasing = TRUE)
  css <- cumsum(u); j <- seq_len(n)
  crit <- u - (css - 1) / j
  rho <- max(j[crit > 0])
  tau <- (css[rho] - 1) / rho
  lam <- pmax(d - tau, 0)
  s <- sum(lam)
  if (s <= 0) lam[] <- 1/n else lam <- lam / s
  lam
}

# -------------------------
# Correct projection to the density set D = {rho ⪰ 0, tr rho = 1}
# -------------------------
project_to_density <- function(Y) {
  YH <- (Y + Conj(t(Y))) / 2
  ev <- eigen(YH, symmetric = TRUE)
  d  <- Re(ev$values)
  U  <- ev$vectors
  lam <- proj_simplex(d)
  rho <- U %*% diag(lam) %*% Conj(t(U))
  (rho + Conj(t(rho))) / 2
}

# -------------------------
# Observables and spectral data
# -------------------------
spectral_projectors <- function(O) {
  # For Hermitian O: O = sum_r lam_r P_r with P_r = |xi_r><xi_r|
  ev <- eigen(O, symmetric = TRUE)
  lams <- Re(ev$values); U <- ev$vectors
  Ps <- vector("list", length(lams))
  for (k in seq_along(lams)) {
    u <- matrix(U[,k], ncol=1)
    Ps[[k]] <- u %*% Conj(t(u))
  }
  list(lams = lams, Ps = Ps)
}

observables_to_spectral <- function(observables) {
  lapply(observables, spectral_projectors)
}

build_observables_pauli <- function() {
  # Pauli X, Y, Z (2x2)
  X <- matrix(c(0,1,1,0), 2, 2)
  Y <- matrix(c(0,-1i,1i,0), 2, 2)
  Z <- matrix(c(1,0,0,-1), 2, 2)
  list(X, Y, Z)
}

# -------------------------
# Born rule probabilities and data generation
# -------------------------
born_probs <- function(rho, Ps) {
  # p_r = tr(rho P_r), normalized
  p <- sapply(Ps, function(P) Re(Trace(rho %*% P)))
  p[p < 0] <- 0
  s <- sum(p)
  if (s <= 0) rep(1/length(p), length(p)) else p/s
}

generate_data <- function(N, rho, specs, a_seq=NULL, seed=0L) {
  # Generate (X_n, a_n, b_n) with exact Born rule via projectors Ps
  set.seed(seed)
  K <- length(specs)
  if (is.null(a_seq)) {
    a <- rep(seq_len(K), length.out=N)  # cycle 1..K
  } else {
    a <- as.integer(a_seq); stopifnot(length(a) == N)
  }
  X <- numeric(N)
  b <- integer(N)
  for (n in seq_len(N)) {
    i <- a[n]
    Ps_i   <- specs[[i]]$Ps
    lams_i <- specs[[i]]$lams
    p <- born_probs(rho, Ps_i)                    # Born rule
    r <- sample.int(length(Ps_i), size=1, prob=p) # outcome index
    b[n] <- r
    X[n] <- lams_i[r]                             # measured eigenvalue
  }
  list(X=X, a=a, b=b)
}

counts_from_ab <- function(a, b, specs) {
  # Count N_{i,r} from (a_n, b_n) only
  K <- length(specs)
  counts <- vector("list", K)
  for (i in seq_len(K)) counts[[i]] <- integer(length(specs[[i]]$Ps))
  for (n in seq_along(a)) {
    i <- a[n]; r <- b[n]
    counts[[i]][r] <- counts[[i]][r] + 1L
  }
  counts
}

# -------------------------
# General log-likelihood and gradient (matrix form)
# -------------------------
loglik_rho <- function(rho, specs, counts, eps=1e-12) {
  ll <- 0
  for (i in seq_along(specs)) {
    Ps_i <- specs[[i]]$Ps
    for (r in seq_along(Ps_i)) {
      trip <- Re(Trace(rho %*% Ps_i[[r]]))
      ll <- ll + counts[[i]][r] * log(max(trip, eps))
    }
  }
  as.numeric(ll)
}

grad_loglik_rho <- function(rho, specs, counts, eps=1e-12) {
  d <- nrow(rho)
  G <- matrix(0+0i, d, d)
  for (i in seq_along(specs)) {
    Ps_i <- specs[[i]]$Ps
    for (r in seq_along(Ps_i)) {
      trip <- Re(Trace(rho %*% Ps_i[[r]]))
      w <- counts[[i]][r] / max(trip, eps)
      G <- G + w * Ps_i[[r]]
    }
  }
  (G + Conj(t(G))) / 2
}

# -------------------------
# PGD-MLE solver (projected gradient ascent on rho)
# -------------------------
pgd_mle_rho <- function(specs, counts,
                        rho_init = NULL,
                        step0 = 1.0, bt_beta = 0.5,   # backtracking params
                        bt_c = 1e-4,                  # Armijo-like slope factor (mild)
                        maxit = 1000L, tol = 1e-8, tmin = 1e-12,
                        verbose = TRUE) {
  d <- nrow(specs[[1]]$Ps[[1]])
  if (is.null(rho_init)) rho <- diag(d)/d else rho <- rho_init
  rho <- project_to_density(rho)   # ensure feasibility
  ll <- loglik_rho(rho, specs, counts)
  if (verbose) cat(sprintf("init ll = %.6f\n", ll))
  hist <- numeric(0)
  
  for (it in seq_len(maxit)) {
    G <- grad_loglik_rho(rho, specs, counts)
    t <- step0
    improved <- FALSE
    # backtracking to ensure monotone increase
    while (t >= tmin) {
      cand <- project_to_density(rho + t * G)
      ll_new <- loglik_rho(cand, specs, counts)
      # Armijo-like condition with HS inner product
      sufficient <- ll_new >= ll + bt_c * t * hs_inner(G, cand - rho)
      if (isTRUE(sufficient) || ll_new > ll) {
        improved <- TRUE; break
      }
      t <- t * bt_beta
    }
    if (!improved) {
      if (verbose) cat("No further improvement; stopping.\n")
      break
    }
    diff_norm <- sqrt(hs_inner(cand - rho, cand - rho))
    rho <- cand; ll <- ll_new; hist <- c(hist, ll)
    if (verbose) cat(sprintf("it=%d  ll=%.6f  step=%.3e  ||Δρ||_F=%.3e\n", it, ll, t, diff_norm))
    if (diff_norm < tol) break
  }
  list(rho_hat = rho, loglik = ll, history = hist, iters = length(hist))
}

# -------------------------
# Optional: Bloch vector from rho (Pauli case) — for reporting only
# -------------------------
bloch_from_rho <- function(rho) {
  X <- matrix(c(0,1,1,0), 2, 2)
  Y <- matrix(c(0,-1i,1i,0), 2, 2)
  Z <- matrix(c(1,0,0,-1), 2, 2)
  c(Re(Trace(rho %*% X)),
    Re(Trace(rho %*% Y)),
    Re(Trace(rho %*% Z)))
}

# -------------------------
# Minimal example / 最小示例
# -------------------------
if (sys.nframe() == 0) {
  # 1) Build observables (Pauli) and spectral data
  observables <- build_observables_pauli()
  specs <- observables_to_spectral(observables)
  
  # 2) Make a random Hermitian and project to a valid rho_true
  set.seed(123)
  M <- matrix(rnorm(4) + 1i*rnorm(4), 2, 2)
  M <- (M + Conj(t(M))) / 2
  rho_true <- project_to_density(M)
  
  # 3) Generate data with a_n cycling 1,2,3,...
  N <- 5000L
  dat <- generate_data(N, rho_true, specs, seed = 42)
  counts <- counts_from_ab(dat$a, dat$b, specs)
  
  # 4) MLE by PGD
  fit <- pgd_mle_rho(specs, counts,
                     rho_init = diag(2)/2,
                     step0 = 1.0, bt_beta = 0.5, bt_c = 1e-4,
                     maxit = 500, tol = 1e-9, verbose = TRUE)
  
  rho_hat <- fit$rho_hat
  cat("\n---- Results ----\n")
  cat("trace(rho_hat) =", Re(Trace(rho_hat)),
      "  min eig =", min(Re(eigen((rho_hat+Conj(t(rho_hat)))/2, symmetric=TRUE)$values)), "\n")
  cat("Frobenius ||rho_hat - rho_true|| =", 
      sqrt(hs_inner(rho_hat - rho_true, rho_hat - rho_true)), "\n")
  
  # Optional: compare Bloch vectors (Pauli case)
  cat("Bloch (true): ", round(bloch_from_rho(rho_true), 4), "\n")
  cat("Bloch (hat) : ", round(bloch_from_rho(rho_hat),  4), "\n")
}
