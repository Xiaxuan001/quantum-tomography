## Utilities ---------------------------------------------------------------

# Standard matrix unit |j><k| in C^{N x N}
E_jk <- function(N, j, k) {
  M <- matrix(0+0i, nrow = N, ncol = N)
  M[j, k] <- 1+0i
  M
}

# Complex trace (works for Hermitian and general complex)
ctrace <- function(A) sum(diag(A))

# Hilbert–Schmidt inner product tr(A^\dagger B)
hs_inner <- function(A, B) ctrace(Conj(t(A)) %*% B)

## SU(N) basis construction: {u_jk, v_jk, w_l} ----------------------------PDF 7 Bloch-vector basis and admissible set 

# Build {sigma_i} with the normalization tr(sigma_i sigma_j) = 2 delta_{ij}
# and tr(sigma_i) = 0, exactly matching the requested definitions.
build_suN_basis <- function(N, verify = TRUE) {  #for basis of density matrix  build sigma
  stopifnot(N >= 2)
  
  sigmas  <- list()
  names_v <- character(0)
  
  # Off-diagonals: u_{jk}, v_{jk},  1 <= j < k <= N
  for (j in 1:(N-1)) {
    for (k in (j+1):N) {
      U <- E_jk(N, j, k) + E_jk(N, k, j)                 # |j><k| + |k><j|
      V <- -1i * (E_jk(N, j, k) - E_jk(N, k, j))         # -i(|j><k| - |k><j|)
      
      sigmas[[length(sigmas) + 1]] <- U
      names_v <- c(names_v, sprintf("u_%d_%d", j, k))
      
      sigmas[[length(sigmas) + 1]] <- V
      names_v <- c(names_v, sprintf("v_%d_%d", j, k))
    }
  }
  
  # Diagonals: w_l for l = 1, ..., N-1
  for (l in 1:(N-1)) {
    # Your definition:
    # w_l = sqrt(2/(l(l+1))) * sum_{m=1}^l (|m><m| - |l+1><l+1|)
    W <- matrix(0+0i, N, N)
    for (m in 1:l) {
      W <- W + E_jk(N, m, m) - E_jk(N, l+1, l+1)
    }
    W <- sqrt(2 / (l*(l+1))) * W
    
    sigmas[[length(sigmas) + 1]] <- W
    names_v <- c(names_v, sprintf("w_%d", l))
  }
  
  names(sigmas) <- names_v
  
  if (verify) {
    # Check tracelessness
    tr0 <- vapply(sigmas, function(S) isTRUE(all.equal(as.numeric(Re(ctrace(S))), 0, tolerance=1e-10)),
                  logical(1))
    if (!all(tr0)) stop("Some sigma have nonzero trace.")
    
    # Check orthonormality: tr(sigma_i sigma_j) = 2 delta_{ij}
    d <- length(sigmas)
    # For speed, check on a sample and diagonals; set verify = FALSE to skip fully.
    diag_ok <- vapply(sigmas, function(S) {
      val <- hs_inner(S, S)
      isTRUE(all.equal(as.numeric(Re(val)), 2, tolerance = 1e-10)) &&
        isTRUE(all.equal(as.numeric(Im(val)), 0, tolerance = 1e-10))
    }, logical(1))
    if (!all(diag_ok)) stop("Some sigma are not normalized to tr(S S)=2.")
    
    # (Optional) spot check a few off-diagonal inner products
    if (d >= 4) {
      idx <- utils::head(which(diag_ok), 4)
      for (i in idx) for (j in idx) if (i != j) {
        val <- hs_inner(sigmas[[i]], sigmas[[j]])
        if (abs(val) > 1e-10) stop("Off-diagonal inner product not ~0.")
      }
    }
  }
  
  sigmas
}


gram <- function(sigmas) { 
  n <- length(sigmas)
  G <- matrix(0+0i, n, n)
  for (i in seq_len(n))
    for (j in seq_len(n))
      G[i, j] <- hs_inner(sigmas[[i]], sigmas[[j]])
  G
}

gram(build_suN_basis(2))
gram(build_suN_basis(3))
