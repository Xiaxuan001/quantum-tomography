proj_simplex <- function(d) {
  d <- as.numeric(d)                                      # EN: ensure numeric; 中文：转数值向量
  n <- length(d)                                          # EN: dimension; 中文：维度
  u <- sort(d, decreasing = TRUE)                         # EN: sort desc; 中文：降序排序
  css <- cumsum(u)                                        # EN: cumulative sums; 中文：前缀和
  # EN: rho = max{ j : u_j - (css_j - 1)/j > 0 }.
  # 中文：ρ = 满足 u_j - (css_j - 1)/j > 0 的最大 j。
  j <- seq_len(n)
  crit <- u - (css - 1) / j
  rho <- max(j[crit > 0])
  tau <- (css[rho] - 1) / rho                             # EN: threshold τ; 中文：阈值 τ
  lam <- pmax(d - tau, 0)                                 # EN: λ = max(d - τ, 0); 中文：λ = max(d-τ,0)
  s <- sum(lam)
  if (s <= 0) lam[] <- 1/n else lam <- lam / s           # EN: renormalize; 中文：重标使和为1
  lam
}