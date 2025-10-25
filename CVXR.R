## ------------------------------------------------------------
## Bloch-basis tomography via theta-MLE with CVXR
## ------------------------------------------------------------
suppressPackageStartupMessages({
  library(CVXR)
  library(Matrix)
})

## ---------- Utilities ----------
traceC <- function(A) Re(sum(diag(A)))                # complex-safe trace
hermitianize <- function(A) (A + Conj(t(A)))/2        # enforce Hermiticity

## Real embedding:  C^{N x N} --> R^{2N x 2N}
##   M ↦ [[Re M, -Im M],
##         [Im M,  Re M]]
real_embed <- function(M) {
  R <- Re(M); I <- Im(M)
  rbind(cbind(R, -I),
        cbind(I,  R))
}

## Elementary matrix |j><k|
E_jk <- function(j, k, N) {
  M <- matrix(0+0i, N, N); M[j, k] <- 1+0i; M
}

## ---------- 1) SU(N) Bloch-vector basis {u_jk, v_jk, w_l} ----------
## Matches Kimura (Eqs. 8a–8c) and your Section~\ref{sec:bloch_basis}.
build_suN_basis <- function(N) {
  stopifnot(N >= 2)
  sigmas <- list()
  labels <- character()
  add <- function(M, lab) { sigmas[[length(sigmas)+1]] <<- hermitianize(M); labels <<- c(labels, lab) }
  
  ## Off-diagonal pairs: 1 <= j < k <= N
  for (j in 1:(N-1)) for (k in (j+1):N) {
    U_jk <- E_jk(j,k,N) + E_jk(k,j,N)                      # \hat u_{jk}
    V_jk <- -1i*(E_jk(j,k,N) - E_jk(k,j,N))                # \hat v_{jk}
    add(U_jk, sprintf("u_%d_%d", j, k))
    add(V_jk, sprintf("v_%d_%d", j, k))
  }
  ## Diagonal traceless set: \hat w_l, 1 <= l <= N-1
  for (l in 1:(N-1)) {
    coef <- sqrt(2/(l*(l+1)))
    diag_vec <- c(rep(1, l), -l, rep(0, N-l-1))
    W_l <- diag(diag_vec) * coef                             # \hat w_l
    add(W_l, sprintf("w_%d", l))
  }
  
  d <- length(sigmas)                                       # should be N^2 - 1
  stopifnot(d == N^2 - 1)
  
  ## Normalize to Tr(σ_i σ_j) = 2 δ_ij  (already satisfied by this construction)
  ## Optionally check numerically for small N:
  # G <- matrix(0, d, d)
  # for (i in 1:d) for (j in 1:d) G[i,j] <- traceC(sigmas[[i]] %*% sigmas[[j]])
  # print(round(G[1:min(d,8),1:min(d,8)], 6))
  
  list(sigmas = sigmas, labels = labels, N = N, d = d)
}

## ---------- 2) Stable spectral projectors for each B_a ----------
## Input: Hermitian B, tolerance 'tol' for grouping eigenvalues (relative).
## Output: list(projectors = {Q_b}, lambdas = eigenvalues (unique), r = count)
spectral_projectors <- function(B, tol = 1e-10) {
  N <- nrow(B)
  H <- hermitianize(B)
  ee <- eigen(H, symmetric = TRUE)             # Hermitian -> stable, real eigenvalues
  vals <- ee$values
  vecs <- ee$vectors
  
  ## Group eigenvalues by relative tolerance
  ord <- order(vals); vals <- vals[ord]; vecs <- vecs[, ord, drop = FALSE]
  groups <- list()
  g_start <- 1
  for (i in 2:length(vals)) {
    if (abs(vals[i] - vals[g_start]) > tol * max(1, abs(vals[g_start]))) {
      groups[[length(groups)+1]] <- g_start:(i-1)
      g_start <- i
    }
  }
  groups[[length(groups)+1]] <- g_start:length(vals)
  
  Qlist <- vector("list", length(groups))
  lambdas <- sapply(groups, function(idx) mean(vals[idx]))
  for (b in seq_along(groups)) {
    idx <- groups[[b]]
    V   <- vecs[, idx, drop = FALSE]
    ## Projector onto the eigenspace: Q = V V^†
    Qb  <- V %*% Conj(t(V))
    Qlist[[b]] <- hermitianize(Qb)
  }
  list(projectors = Qlist, lambdas = lambdas, r = length(Qlist))
}

## Build measurement list: B_a = sigma_a; then eigen-decompose each
build_measurements_from_basis <- function(sigmas, eig_tol = 1e-10) {
  d <- length(sigmas)
  Q_list   <- vector("list", d)       # each entry is list of Q_{a,b}
  lam_list <- vector("list", d)       # lambdas for each 'a'
  r_vec    <- integer(d)
  for (a in 1:d) {
    sp <- spectral_projectors(sigmas[[a]], tol = eig_tol)
    Q_list[[a]]   <- sp$projectors
    lam_list[[a]] <- sp$lambdas
    r_vec[a]      <- sp$r
  }
  ## Create a handy (a,b)->row index map
  ab_row <- vector("list", d)
  row <- 0L
  for (a in 1:d) {
    ab_row[[a]] <- (row + 1L):(row + r_vec[a]); row <- row + r_vec[a]
  }
  ab_df <- do.call(rbind, lapply(1:d, function(a) {
    data.frame(a = a, b = seq_len(r_vec[a]), row = ab_row[[a]], lambda = lam_list[[a]])
  }))
  list(Q_list = Q_list, lam_list = lam_list, r_vec = r_vec, ab_row = ab_row, ab_df = ab_df)
}

## ---------- 3) Born probabilities and sampling ----------
## ρ(θ) = I/N + 1/2 Σ_j θ_j σ_j   (complex N x N)
rho_of_theta <- function(theta, sigmas, N) {
  R <- diag(N)/N
  if (length(theta) > 0) {
    for (j in seq_along(sigmas)) R <- R + 0.5 * theta[j] * sigmas[[j]]
  }
  hermitianize(R)
}

## p_{a,b}(ρ) = Tr(ρ Q_{a,b})
born_probs <- function(rho, Q_list) {
  sapply(Q_list, function(Qb) sapply(Qb, function(Q) {
    val <- traceC(rho %*% Q); as.numeric(Re(val))
  }))
  ## returns a matrix with columns 'a', but lengths differ across 'a'; we will use list form instead:
}

born_probs_list <- function(rho, Q_list) {
  lapply(Q_list, function(Qb) {
    p <- vapply(Qb, function(Q) as.numeric(Re(traceC(rho %*% Q))), numeric(1))
    ## Normalize defensively (tolerant to tiny numeric drift):
    p <- p / sum(p)
    p
  })
}

## Sequential schedule: a_1=1,...,a_d=d,a_{d+1}=1,...
sequential_actions <- function(d, n) ((0:(n-1)) %% d) + 1L

## Sample b_n from Born rule for each chosen a_n
sample_b_sequence <- function(a_seq, prob_by_a, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(a_seq)
  b <- integer(n)
  for (m in 1:n) {
    a <- a_seq[m]
    b[m] <- sample.int(length(prob_by_a[[a]]), size = 1, prob = prob_by_a[[a]])
  }
  b
}

## ---------- 4) Affine probabilities p_{a,b}(θ) = c_ab + 1/2 * S_ab %*% θ ----------
## S_ab[m,j] = Tr(σ_j Q_{a,b}), c_ab[m] = (1/N) Tr(Q_{a,b})
build_Sab_cab <- function(sigmas, Q_list, N, ab_df) {
  d <- length(sigmas)
  M <- nrow(ab_df)
  S_ab <- matrix(0, M, d)
  c_ab <- numeric(M)
  for (m in 1:M) {
    a <- ab_df$a[m]; b <- ab_df$b[m]
    Q <- Q_list[[a]][[b]]
    c_ab[m] <- traceC(Q) / N
    for (j in 1:d) S_ab[m, j] <- traceC(sigmas[[j]] %*% Q)   # real by construction
  }
  list(S_ab = S_ab, c_ab = c_ab)
}

## Optional: numeric NLL for diagnostics
nll_numeric <- function(theta, S_ab, c_ab, Nab, eps = 1e-12) {
  p <- as.numeric(c_ab + 0.5 * (S_ab %*% theta))
  if (any(p < -1e-8)) return(Inf)
  p <- pmax(p, eps)
  -sum(Nab * log(p))
}

## ---------- 5) CVXR MLE in θ with spectrahedral constraint ----------
## A(θ) = real_embed( ρ(θ) )  ⪰ 0
## p(θ) = c_ab + 0.5 * S_ab %*% θ
fit_theta_cvxr <- function(N, sigmas, S_ab, c_ab, Nab,
                           solver = c("MOSEK","ECOS","SCS"), eps_log = 1e-12, verbose = FALSE) {
  d  <- length(sigmas)
  M  <- length(c_ab)
  theta <- Variable(d)
  ## Affine probabilities
  S_const <- as.matrix(S_ab)
  p <- c_ab + 0.5 * (S_const %*% theta)
  ## Real-embedded PSD constraint on rho(θ)
  SI <- diag(2*N)   # real-embedded identity as dense matrix
  Slist <- lapply(sigmas, real_embed)
  A_affine <- SI / N
  for (j in 1:d) A_affine <- A_affine + 0.5 * theta[j] * Slist[[j]]
  rho <- Variable(2 * N, 2 * N, PSD = TRUE)
  
  ## Objective and constraints
  obj <- -sum_entries(Nab * log(p + eps_log))
  constr <- list(rho == A_affine, p >= 0)  # rho PSD via variable, p>=0 aids numerics
  prob <- Problem(Minimize(obj), constr)
  
  solver <- match.arg(solver)
  ## Fallback chain if a premium solver isn't available
  if (solver == "MOSEK" && !requireNamespace("Rmosek", quietly = TRUE)) solver <- "ECOS"
  if (solver == "ECOS"  && !requireNamespace("ECOSolveR", quietly = TRUE)) solver <- "SCS"
  
  res <- solve(prob, solver = solver, verbose = verbose,
               feastol = 1e-8, reltol = 1e-8, abstol = 1e-8,
               max_iters = 5000)
  list(theta_hat = drop(res$getValue(theta)),
       status = res$status,
       value  = res$value)
}

## ---------- Convenience: random density, true theta, counts ----------
random_density <- function(N, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(rnorm(N*N) + 1i*rnorm(N*N), N, N)
  Y <- X %*% Conj(t(X))
  Y / traceC(Y)
}

theta_from_rho <- function(rho, sigmas) {
  sapply(sigmas, function(S) as.numeric(Re(traceC(rho %*% S))))
}

## Build counts Nab for the (a,b) cells
counts_from_ab <- function(a_seq, b_seq, ab_row, M) {
  Nab <- numeric(M)
  for (m in 1:length(a_seq)) {
    a <- a_seq[m]; b <- b_seq[m]
    r <- ab_row[[a]][b]
    Nab[r] <- Nab[r] + 1
  }
  Nab
}

## ---------- End-to-end demo wrapper ----------
run_demo <- function(N = 3, n_per_d = 10, eig_tol = 1e-10,
                     solver = c("MOSEK","ECOS","SCS"), seed = 123,
                     rho_true = NULL, theta_true = NULL) {
  solver <- match.arg(solver)
  cat(sprintf("Building SU(%d) basis ...\n", N))
  bas <- build_suN_basis(N); sigmas <- bas$sigmas; d <- bas$d
  
  cat("Building measurement list B_a = sigma_a and spectral projectors ...\n")
  meas <- build_measurements_from_basis(sigmas, eig_tol = eig_tol)
  Q_list <- meas$Q_list; ab_df <- meas$ab_df; ab_row <- meas$ab_row
  M <- nrow(ab_df)
  
  cat("Preparing affine model p_ab(theta) = c_ab + 0.5 * S_ab %*% theta ...\n")
  SC <- build_Sab_cab(sigmas, Q_list, N, ab_df)
  S_ab <- SC$S_ab; c_ab <- SC$c_ab
  
  cat("Preparing ground-truth state ...\n")
  ## Derive a consistent (rho_true, theta_true) pair from the supplied inputs.
  if (!is.null(theta_true)) {
    rho_true <- rho_of_theta(theta_true, sigmas, N)
  } else if (!is.null(rho_true)) {
    theta_true <- theta_from_rho(rho_true, sigmas)
  } else {
    rho_true   <- random_density(N, seed = seed)
    theta_true <- theta_from_rho(rho_true, sigmas)
  }
  
  a_seq <- sequential_actions(d, n = d*n_per_d)
  p_list_true <- born_probs_list(rho_true, Q_list)
  b_seq <- sample_b_sequence(a_seq, p_list_true, seed = seed+1L)
  Nab   <- counts_from_ab(a_seq, b_seq, ab_row, M)
  
  cat("Solving theta-MLE with CVXR ...\n")
  fit <- fit_theta_cvxr(N, sigmas, S_ab, c_ab, Nab, solver = solver)
  theta_hat <- fit$theta_hat
  
  cat(sprintf("Solver status: %s,  optimal value: %.6f\n", fit$status, fit$value))
  rel_err <- sqrt(sum((theta_hat - theta_true)^2)) / max(1e-9, sqrt(sum(theta_true^2)))
  cat(sprintf("||theta_hat - theta_true||_2 / ||theta_true||_2 = %.3e\n", rel_err))
  
  list(N=N, d=d, ab_df=ab_df, Nab=Nab,
       a_seq=a_seq, b_seq=b_seq,
       theta_true=theta_true, theta_hat=theta_hat,
       rho_true=rho_true,
       S_ab=S_ab, c_ab=c_ab, solver_status=fit$status)
}

## --------------- Example run ---------------
## Execute a demo only when run as a script, to keep sourcing lightweight.
if (sys.nframe() == 0) {
  out <- run_demo(N = 2, n_per_d = 10, solver = "SCS", seed = 42)
  str(out, max.level = 1)
}
