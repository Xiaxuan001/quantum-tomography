# Flatten a nested list of projectors into a single list
flatten_list <- function(x) {
  if (!is.list(x)) return(list(x))
  unlist(x, recursive = TRUE)
}

#**** Build c_ab (vector) and S_ab (matrix with rows s(a,b)^T) from complex projectors Q_{a,b}
# Qs may be a list-of-lists by setting or a flat list; each Q must be N x N Hermitian projector.
build_linear_prob_parts <- function(Qs, sigma_list, N) {
  Qflat <- flatten_list(Qs)
  M <- length(Qflat)
  d <- length(sigma_list)
  
  c_ab <- numeric(M)
  S_ab <- matrix(0, nrow = M, ncol = d)
  
  for (m in 1:M) {
    Q <- Qflat[[m]]
    c_ab[m] <- as.numeric(Re(ctrace(Q)) / N)
    for (j in 1:d) {
      S_ab[m, j] <- as.numeric(Re(ctrace(sigma_list[[j]] %*% Q)))
    }
  }
  list(c_ab = c_ab, S_ab = S_ab)
}
