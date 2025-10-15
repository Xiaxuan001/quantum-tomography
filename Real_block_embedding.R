library(Matrix)

# Complex -> real block embedding
# Returns a (2N x 2N) real symmetric Matrix suitable for CVXR's PSD cone.
real_embed <- function(A, sparse = TRUE) {
  ReA <- Re(A)
  ImA <- Im(A)
  
  top    <- cbind(ReA, -ImA)
  bottom <- cbind(ImA,  ReA)
  M      <- rbind(top, bottom)
  
  if (sparse) {
    M <- Matrix(M, sparse = TRUE)
  } else {
    M <- Matrix(M, sparse = FALSE)
  }
  # Force symmetry numerically (Hermitian A -> symmetric embedding)
  M <- forceSymmetric((M + t(M)) / 2)
  M
}

# Convert a list of complex sigma_j to a list of embedded real matrices
# Slist[[j]] = Phi(sigma_j)
embed_sigma_list <- function(sigma_list, sparse = TRUE) {
  lapply(sigma_list, function(S) real_embed(S, sparse = sparse))
}

# Real-embedded identity: Phi(I_N) = I_{2N}
real_embed_identity <- function(N, sparse = TRUE) {
  if (sparse) Diagonal(2*N) else Matrix(diag(2*N), sparse = FALSE)
}

