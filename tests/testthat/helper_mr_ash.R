# Shared helpers for the mr.ash tests
# (test_mr_ash.R + test_mr_ash_equivalence.R).
# testthat auto-sources helper_*.R before running any test file.
#
# Naming map (old -> unified):
#   make_data            (test_mr_ash.R)            -> mr_ash_data()
#   setup_mr_ash_test    (test_mr_ash_equivalence.R)-> mr_ash_data()
#   derive_summary_stats (test_mr_ash_equivalence.R)-> derive_summary_stats()
#
# mr_ash_data() is the canonical generator: it reproduces setup_mr_ash_test's
# exact construction (so the numeric-equivalence tolerances still hold) and
# also returns beta_true and k for tests that need them.

# Generate centered X, y and a default mr.ash prior (sa2 / pi0 / sigma2_init).
mr_ash_data <- function(n = 100, p = 50, k = 5, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = FALSE)
  beta_true <- rep(0, p)
  causal <- sample(1:p, k)
  beta_true[causal] <- rnorm(k, sd = 2)
  y <- c(X %*% beta_true + rnorm(n))
  y <- y - mean(y)

  # Prior matching mr.ash defaults.
  sa2 <- c(0, (2^((1:19) / 20) - 1)^2)
  w <- colSums(X^2)
  sa2 <- sa2 / median(w) * n
  K <- length(sa2)
  pi0 <- rep(1 / K, K)
  sigma2_init <- c(var(y))

  list(X = X, y = y, n = n, p = p, k = k,
       sa2 = sa2, K = K, pi0 = pi0, sigma2_init = sigma2_init,
       beta_true = beta_true)
}

# Derive summary statistics from individual-level data, matching the PVE
# adjustment in mr.ash.rss (uses n-2 df for shat):
#   bhat_j = X_j'y / X_j'X_j,  shat_j = sqrt(RSS_j / ((n-2) * X_j'X_j)),
#   R = cor(X),  var_y = var(y).
derive_summary_stats <- function(X, y) {
  n <- nrow(X)
  p <- ncol(X)
  bhat <- sapply(1:p, function(j) sum(X[, j] * y) / sum(X[, j]^2))
  shat <- sapply(1:p, function(j) {
    resid <- y - X[, j] * bhat[j]
    sqrt(sum(resid^2) / ((n - 2) * sum(X[, j]^2)))
  })
  R_mat <- cor(X)
  var_y <- c(var(y))
  list(bhat = bhat, shat = shat, R = R_mat, var_y = var_y, n = n)
}
