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
  
  list(N=N, d=d, ab_df=ab_df, ab_row=ab_row, Nab=Nab,
       a_seq=a_seq, b_seq=b_seq,
       theta_true=theta_true, theta_hat=theta_hat,
       rho_true=rho_true, sigmas=sigmas, Q_list=Q_list,
       S_ab=S_ab, c_ab=c_ab, solver_status=fit$status)
}

## --------------- Example run ---------------
## Execute a demo only when run as a script, to keep sourcing lightweight.
if (sys.nframe() == 0) {
  out <- run_demo(N = 2, n_per_d = 10, solver = "SCS", seed = 42)
  str(out, max.level = 1)
}


###########6.4 

## ---------- 6) Fisher information for a single observation (Eq. 6.4) ----------
## Implements: I_a(theta) = (1/4) * sum_b s(a,b) s(a,b)^T / p_{a,b}(theta)
## Where: s(a,b) are rows of S_a, and p_{a,b}(theta) = c_{a,b} + 0.5 * S_a(b,.) %*% theta

## Helper: pick the (a,b) block rows from (S_ab, c_ab) using ab_df
get_Sc_for_setting <- function(a, S_ab, c_ab, ab_df) {
  idx <- which(ab_df$a == a)
  if (length(idx) == 0L) stop(sprintf("No rows found for setting a = %d.", a))
  list(S_a = S_ab[idx, , drop = FALSE],
       c_a = c_ab[idx],
       rows = idx)
}

## ---- (6.4) Single-setting Fisher information I_a(theta)
fisher_info_setting <- function(theta, a, S_ab, c_ab, ab_df, eps = 1e-12, check_positive = TRUE) {
  sc <- get_Sc_for_setting(a, S_ab, c_ab, ab_df)
  S_a <- sc$S_a                                    # r_a x d matrix with entries s_j(a,b)
  p_a <- sc$c_a + 0.5 * as.numeric(S_a %*% theta)  # length r_a

  if (check_positive && any(p_a <= 0)) {
    stop(sprintf("p_{a,b}(theta) has nonpositive entries at setting a = %d. "+
                 "Consider projecting theta into the feasible set or increase eps.", a))
  }
  p_a <- pmax(p_a, eps)                            # numerical safety
  ## Matrix form: I_a = (1/4) * t(S_a) %*% diag(1/p_a) %*% S_a
  Ia <- 0.25 * crossprod(S_a, S_a * (1 / p_a))     # row-wise scaling by 1/p_a
  Ia
}

## ---- Convenience: Fisher info for all settings (list)
fisher_info_by_setting <- function(theta, S_ab, c_ab, ab_df, eps = 1e-12, check_positive = TRUE) {
  A <- sort(unique(ab_df$a))
  setNames(lapply(A, function(a)
    fisher_info_setting(theta, a, S_ab, c_ab, ab_df, eps, check_positive)), A)
}

## ---- Weighted average Fisher info over settings: sum_a pi[a] * I_a(theta)
## pi can be a probability vector over settings, or counts normalized by sum.
fisher_info_weighted <- function(theta, pi, S_ab, c_ab, ab_df, eps = 1e-12, check_positive = TRUE) {
  A <- sort(unique(ab_df$a))
  if (length(pi) != length(A)) stop("Length of 'pi' must match number of unique settings.")
  if (any(pi < 0)) stop("'pi' must be nonnegative.")
  if (sum(pi) == 0) stop("Sum of 'pi' must be positive.")
  pi <- pi / sum(pi)

  I_list <- fisher_info_by_setting(theta, S_ab, c_ab, ab_df, eps, check_positive)
  Reduce(`+`, Map(function(w, Ia) w * Ia, pi, I_list))
}

## ---- Helper: counts by setting from Nab
counts_by_setting <- function(Nab, ab_df) {
  A <- sort(unique(ab_df$a))
  vapply(A, function(a) sum(Nab[ab_df$a == a]), numeric(1))
}

## ---- Helper: total Fisher information given counts
total_fisher_info <- function(theta, counts, S_ab, c_ab, ab_df,
                              eps = 1e-12, check_positive = TRUE) {
  A <- sort(unique(ab_df$a))
  counts_full <- numeric(length(A))
  names(counts_full) <- as.character(A)
  if (length(counts) > 0) {
    if (is.null(names(counts)) && length(counts) == length(A)) {
      names(counts) <- as.character(A)
    }
    overlap <- intersect(names(counts), names(counts_full))
    counts_full[overlap] <- counts[overlap]
  }
  I_list <- fisher_info_by_setting(theta, S_ab, c_ab, ab_df,
                                   eps = eps, check_positive = check_positive)
  d <- ncol(S_ab)
  I_total <- matrix(0, d, d)
  for (nm in names(I_list)) {
    n_a <- counts_full[[nm]]
    if (n_a > 0) I_total <- I_total + n_a * I_list[[nm]]
  }
  list(I_total = I_total, I_list = I_list, counts = counts_full)
}

## ---- Helpers for design criteria
safe_logdet <- function(M, ridge = 1e-8) {
  d <- nrow(M)
  M_reg <- M + diag(ridge, d)
  det_res <- suppressWarnings(determinant(M_reg, logarithm = TRUE))
  if (det_res$sign <= 0) return(-Inf)
  as.numeric(det_res$modulus)
}

safe_trace_inverse <- function(M, ridge = 1e-8) {
  d <- nrow(M)
  M_reg <- M + diag(ridge, d)
  Minv <- tryCatch(solve(M_reg), error = function(e) NULL)
  if (is.null(Minv)) return(Inf)
  sum(diag(Minv))
}

first_order_d_score <- function(I_total, Ia, ridge = 1e-8) {
  d <- nrow(I_total)
  M_reg <- I_total + diag(ridge, d)
  Minv <- tryCatch(solve(M_reg), error = function(e) NULL)
  if (is.null(Minv)) return(-Inf)
  sum(Minv * Ia)
}

first_order_a_score <- function(I_total, Ia, ridge = 1e-8) {
  d <- nrow(I_total)
  M_reg <- I_total + diag(ridge, d)
  Minv <- tryCatch(solve(M_reg), error = function(e) NULL)
  if (is.null(Minv)) return(-Inf)
  sum(diag(Minv %*% Ia %*% Minv))
}

select_next_setting <- function(theta, S_ab, c_ab, ab_df, counts = NULL,
                                criterion = c("D", "A"),
                                method = c("exact", "first-order"),
                                ridge = 1e-8,
                                eps = 1e-12,
                                check_positive = TRUE) {
  criterion <- match.arg(criterion)
  method <- match.arg(method)
  tf <- total_fisher_info(theta, counts, S_ab, c_ab, ab_df,
                          eps = eps, check_positive = check_positive)
  I_total <- tf$I_total
  I_list <- tf$I_list
  scores <- numeric(length(I_list))
  names(scores) <- names(I_list)
  method_used <- method
  for (nm in names(I_list)) {
    Ia <- I_list[[nm]]
    if (method == "exact") {
      if (criterion == "D") {
        scores[[nm]] <- safe_logdet(I_total + Ia, ridge = ridge)
      } else {
        scores[[nm]] <- -safe_trace_inverse(I_total + Ia, ridge = ridge)
      }
    } else {
      if (criterion == "D") {
        scores[[nm]] <- first_order_d_score(I_total, Ia, ridge = ridge)
      } else {
        scores[[nm]] <- first_order_a_score(I_total, Ia, ridge = ridge)
      }
    }
  }
  if (method == "exact" && !any(is.finite(scores))) {
    scores <- rep(-Inf, length(I_list))
    names(scores) <- names(I_list)
    method_used <- "first-order (fallback)"
    for (nm in names(I_list)) {
      Ia <- I_list[[nm]]
      if (criterion == "D") {
        scores[[nm]] <- first_order_d_score(I_total, Ia, ridge = ridge)
      } else {
        scores[[nm]] <- first_order_a_score(I_total, Ia, ridge = ridge)
      }
    }
  }
  best_nm <- names(scores)[which.max(scores)]
  list(setting = as.integer(best_nm), scores = scores, I_total = I_total, I_list = I_list,
       method_used = method_used, criterion = criterion)
}

## ---- (Optional) Observed (per-outcome) information J_{a,b}(theta)
## J_{a,b} = (1/4) * s(a,b) s(a,b)^T / p_{a,b}(theta)^2   (Eq. (8.1) for reference)
observed_info_one_outcome <- function(theta, a, b, S_ab, c_ab, ab_df, eps = 1e-12, check_positive = TRUE) {
  sc <- get_Sc_for_setting(a, S_ab, c_ab, ab_df)
  if (b < 1L || b > nrow(sc$S_a)) stop(sprintf("b = %d out of range for setting a = %d.", b, a))
  s_ab <- matrix(sc$S_a[b, ], ncol = 1)                    # d x 1
  p_ab <- sc$c_a[b] + 0.5 * sum(sc$S_a[b, ] * theta)
  if (check_positive && p_ab <= 0) stop("p_{a,b}(theta) <= 0.")
  p_ab <- max(p_ab, eps)
  (0.25 / (p_ab^2)) * (s_ab %*% t(s_ab))
}

## ---- Average Fisher information for a given measurement schedule
average_fisher_info_from_sequence <- function(a_seq,
                                              S_ab, c_ab, ab_df,
                                              theta = NULL,
                                              b_seq = NULL,
                                              N = NULL,
                                              sigmas = NULL,
                                              ab_row = NULL,
                                              solver = c("MOSEK","ECOS","SCS"),
                                              eps = 1e-12,
                                              check_positive = TRUE,
                                              eps_log = 1e-12,
                                              verbose = FALSE) {
  if (length(a_seq) == 0L) stop("'a_seq' must contain at least one measurement setting.")
  solver <- match.arg(solver)
  theta_use <- theta
  theta_source <- "input"
  solver_status <- NULL
  solver_value  <- NA_real_
  if (is.null(theta_use)) {
    if (is.null(b_seq)) stop("Provide 'b_seq' when 'theta' is NULL.")
    if (length(a_seq) != length(b_seq)) stop("'a_seq' and 'b_seq' must have the same length.")
    if (is.null(N) || is.null(sigmas) || is.null(ab_row)) {
      stop("To obtain theta via MLE, supply 'N', 'sigmas', and 'ab_row'.")
    }
    Nab <- counts_from_ab(a_seq, b_seq, ab_row, nrow(ab_df))
    fit <- fit_theta_cvxr(N, sigmas, S_ab, c_ab, Nab,
                          solver = solver, eps_log = eps_log, verbose = verbose)
    theta_use <- fit$theta_hat
    theta_source <- "MLE"
    solver_status <- fit$status
    solver_value  <- fit$value
  }
  A <- sort(unique(ab_df$a))
  idx <- match(a_seq, A)
  if (any(is.na(idx))) stop("Some entries in 'a_seq' are not present in 'ab_df$a'.")
  counts <- tabulate(idx, nbins = length(A))
  if (sum(counts) == 0) stop("'a_seq' produced zero counts; check input.")
  weights <- counts / sum(counts)
  names(weights) <- as.character(A)
  I_list <- fisher_info_by_setting(theta_use, S_ab, c_ab, ab_df,
                                   eps = eps, check_positive = check_positive)
  I_avg <- fisher_info_weighted(theta_use, counts, S_ab, c_ab, ab_df,
                                eps = eps, check_positive = check_positive)
  list(theta = theta_use,
       theta_source = theta_source,
       fisher_by_setting = I_list,
       fisher_average = I_avg,
       weights = weights,
       counts = counts,
       solver_status = solver_status,
       solver_value = solver_value)
}

born_probs_from_theta <- function(theta, sigmas, Q_list, N) {
  rho <- rho_of_theta(theta, sigmas, N)
  born_probs_list(rho, Q_list)
}

increment_ab_count <- function(Nab, a, b, ab_row) {
  if (a < 1L || a > length(ab_row)) stop(sprintf("'a' = %d is out of range.", a))
  rows_ab <- ab_row[[a]]
  if (b < 1L || b > length(rows_ab)) stop(sprintf("'b' = %d is out of range for setting a = %d.", b, a))
  idx <- rows_ab[b]
  Nab[idx] <- Nab[idx] + 1
  Nab
}

adaptive_design_sequence <- function(steps,
                                     N, sigmas, Q_list,
                                     S_ab, c_ab, ab_df, ab_row,
                                     Nab_init = NULL,
                                     a_seq_init = integer(0),
                                     b_seq_init = integer(0),
                                     theta_init = NULL,
                                     theta_sampling = NULL,
                                     rho_sampling = NULL,
                                     criterion = c("D", "A"),
                                     method = c("exact", "first-order"),
                                     solver = c("MOSEK","ECOS","SCS"),
                                     ridge = 1e-8,
                                     eps = 1e-12,
                                     eps_log = 1e-12,
                                     verbose = FALSE,
                                     seed = NULL,
                                     check_positive = TRUE) {
  if (steps < 0L) stop("'steps' must be nonnegative.")
  solver <- match.arg(solver)
  criterion <- match.arg(criterion)
  method <- match.arg(method)
  M <- nrow(ab_df)
  if (!is.null(Nab_init)) {
    if (length(Nab_init) != M) stop("'Nab_init' must have length equal to nrow(ab_df).")
    Nab_current <- as.numeric(Nab_init)
  } else if (length(a_seq_init) > 0) {
    if (length(a_seq_init) != length(b_seq_init)) stop("'a_seq_init' and 'b_seq_init' must have the same length.")
    Nab_current <- counts_from_ab(a_seq_init, b_seq_init, ab_row, M)
  } else {
    Nab_current <- numeric(M)
  }
  d <- length(sigmas)
  theta_current <- NULL
  theta_source <- "input"
  if (!is.null(theta_init)) {
    theta_current <- as.numeric(theta_init)
  } else if (sum(Nab_current) > 0) {
    fit0 <- fit_theta_cvxr(N, sigmas, S_ab, c_ab, Nab_current,
                           solver = solver, eps_log = eps_log, verbose = verbose)
    theta_current <- fit0$theta_hat
    theta_source <- sprintf("MLE (status: %s)", fit0$status)
  } else {
    theta_current <- rep(0, d)
    theta_source <- "zero"
  }
  theta_initial <- theta_current
  rho_sampling_use <- NULL
  if (!is.null(rho_sampling)) {
    rho_sampling_use <- rho_sampling
  } else if (!is.null(theta_sampling)) {
    rho_sampling_use <- rho_of_theta(theta_sampling, sigmas, N)
  }
  history <- vector("list", steps)
  if (!is.null(seed)) set.seed(seed)
  for (step in seq_len(steps)) {
    counts_current <- counts_by_setting(Nab_current, ab_df)
    sel <- select_next_setting(theta_current, S_ab, c_ab, ab_df,
                               counts = counts_current,
                               criterion = criterion,
                               method = method,
                               ridge = ridge,
                               eps = eps,
                               check_positive = check_positive)
    a_next <- sel$setting
    rho_for_sampling <- if (!is.null(rho_sampling_use)) {
      rho_sampling_use
    } else {
      rho_of_theta(theta_current, sigmas, N)
    }
    prob_list <- born_probs_list(rho_for_sampling, Q_list)
    probs_a <- prob_list[[a_next]]
    b_next <- sample.int(length(probs_a), size = 1, prob = probs_a)
    Nab_current <- increment_ab_count(Nab_current, a_next, b_next, ab_row)
    fit <- fit_theta_cvxr(N, sigmas, S_ab, c_ab, Nab_current,
                          solver = solver, eps_log = eps_log, verbose = verbose)
    theta_current <- fit$theta_hat
    history[[step]] <- list(
      step = step,
      criterion = sel$criterion,
      method_requested = method,
      method_used = sel$method_used,
      a_next = a_next,
      b_next = b_next,
      probs = probs_a,
      counts_prior = counts_current,
      scores = sel$scores,
      I_total_prior = sel$I_total,
      theta = theta_current,
      solver_status = fit$status,
      solver_value = fit$value
    )
  }
  theta_matrix <- if (steps > 0) {
    do.call(rbind, lapply(history, function(h) h$theta))
  } else {
    NULL
  }
  theta_path <- if (steps > 0) {
    rbind(theta_initial, theta_matrix)
  } else {
    matrix(theta_initial, nrow = 1)
  }
  list(
    Nab = Nab_current,
    counts_per_setting = counts_by_setting(Nab_current, ab_df),
    theta = theta_current,
    theta_initial = theta_initial,
    theta_source_initial = theta_source,
    theta_path = theta_path,
    a_path = if (steps > 0) vapply(history, `[[`, integer(1), "a_next") else integer(0),
    b_path = if (steps > 0) vapply(history, `[[`, integer(1), "b_next") else integer(0),
    criterion = criterion,
    method_requested = method,
    history = history,
    method_used_last = if (steps > 0) history[[steps]]$method_used else method
  )
}

## --------------- Section 6.4 demo ---------------
if (sys.nframe() == 0) {
  cat("\nSection 6.4 example: Fisher information via Eq. (6.4)\n")
  fisher_demo <- average_fisher_info_from_sequence(
    a_seq = out$a_seq,
    b_seq = out$b_seq,
    S_ab = out$S_ab,
    c_ab = out$c_ab,
    ab_df = out$ab_df,
    theta = out$theta_hat,
    N = out$N,
    sigmas = out$sigmas,
    ab_row = out$ab_row,
    solver = "SCS")
  for (nm in names(fisher_demo$fisher_by_setting)) {
    cat(sprintf("Setting a = %s -> I_a(theta_hat):\n", nm))
    print(round(fisher_demo$fisher_by_setting[[nm]], 6))
  }
  cat("Average Fisher information for provided schedule:\n")
  print(round(fisher_demo$fisher_average, 6))

  summarize_design <- function(label, res) {
    cat(sprintf("\n%s\n", label))
    cat(sprintf("  a-path: %s\n", paste(res$a_path, collapse = ", ")))
    cat(sprintf("  b-path: %s\n", paste(res$b_path, collapse = ", ")))
    cat(sprintf("  theta_hat (last): %s\n", paste(round(res$theta, 4), collapse = ", ")))
    cat("  counts per setting:\n")
    print(res$counts_per_setting)
    cat(sprintf("  method used last: %s\n", res$method_used_last))
  }

  cat("\nAdaptive design preview (5 additional shots)\n")
  design_D_exact <- adaptive_design_sequence(
    steps = 5,
    N = out$N,
    sigmas = out$sigmas,
    Q_list = out$Q_list,
    S_ab = out$S_ab,
    c_ab = out$c_ab,
    ab_df = out$ab_df,
    ab_row = out$ab_row,
    Nab_init = out$Nab,
    theta_init = out$theta_hat,
    rho_sampling = out$rho_true,
    criterion = "D",
    method = "exact",
    solver = "SCS",
    seed = 100)
  summarize_design("D-optimal (exact)", design_D_exact)

  design_D_first <- adaptive_design_sequence(
    steps = 5,
    N = out$N,
    sigmas = out$sigmas,
    Q_list = out$Q_list,
    S_ab = out$S_ab,
    c_ab = out$c_ab,
    ab_df = out$ab_df,
    ab_row = out$ab_row,
    Nab_init = out$Nab,
    theta_init = out$theta_hat,
    rho_sampling = out$rho_true,
    criterion = "D",
    method = "first-order",
    solver = "SCS",
    seed = 101)
  summarize_design("D-optimal (first-order approximation)", design_D_first)

  design_A_exact <- adaptive_design_sequence(
    steps = 5,
    N = out$N,
    sigmas = out$sigmas,
    Q_list = out$Q_list,
    S_ab = out$S_ab,
    c_ab = out$c_ab,
    ab_df = out$ab_df,
    ab_row = out$ab_row,
    Nab_init = out$Nab,
    theta_init = out$theta_hat,
    rho_sampling = out$rho_true,
    criterion = "A",
    method = "exact",
    solver = "SCS",
    seed = 102)
  summarize_design("A-optimal (exact)", design_A_exact)

  design_A_first <- adaptive_design_sequence(
    steps = 5,
    N = out$N,
    sigmas = out$sigmas,
    Q_list = out$Q_list,
    S_ab = out$S_ab,
    c_ab = out$c_ab,
    ab_df = out$ab_df,
    ab_row = out$ab_row,
    Nab_init = out$Nab,
    theta_init = out$theta_hat,
    rho_sampling = out$rho_true,
    criterion = "A",
    method = "first-order",
    solver = "SCS",
    seed = 103)
  summarize_design("A-optimal (first-order approximation)", design_A_first)
}
