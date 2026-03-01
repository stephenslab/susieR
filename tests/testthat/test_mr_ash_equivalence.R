# =============================================================================
# Test: mr.ash vs mr.ash.rss Equivalence
# =============================================================================
#
# Verifies that mr.ash (individual-level data) and mr.ash.rss (summary
# statistics) produce equivalent results when fed the same data.
#
# The key mathematical relationship is:
#   mr.ash model: beta_j ~ sum_k pi_k * N(0, sigma2 * sa2[k])
#   mr.ash.rss reconstructs X'X, X'y from (bhat, shat, R, var_y, n)
#   and uses s0 as the prior variance scale (multiplied by sigma2_e internally).
#
# Summary statistics are derived from individual data as:
#   bhat_j = X_j'y / X_j'X_j  (univariate OLS)
#   shat_j = sqrt(RSS_j / ((n-2) * X_j'X_j))  (standard error with n-2 df)
#   R = cor(X)
#   var_y = var(y)
# =============================================================================

# Helper: derive summary statistics from individual data
# Uses n-2 df for shat to match the PVE adjustment in mr.ash.rss
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

# Helper: generate test data and prior
setup_mr_ash_test <- function(n = 100, p = 50, k = 5, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = FALSE)
  beta_true <- rep(0, p)
  causal <- sample(1:p, k)
  beta_true[causal] <- rnorm(k, sd = 2)
  y <- c(X %*% beta_true + rnorm(n))
  y <- y - mean(y)

  # Prior matching mr.ash defaults
  sa2 <- c(0, (2^((1:19) / 20) - 1)^2)
  w <- colSums(X^2)
  sa2 <- sa2 / median(w) * n
  K <- length(sa2)
  pi0 <- rep(1 / K, K)
  sigma2_init <- c(var(y))

  list(
    X = X, y = y, n = n, p = p,
    sa2 = sa2, K = K, pi0 = pi0, sigma2_init = sigma2_init
  )
}

test_that("mr.ash and mr.ash.rss produce identical beta with fixed sigma and pi", {
  d <- setup_mr_ash_test(n = 100, p = 50, k = 5, seed = 42)

  fit_ind <- mr.ash(d$X, d$y,
    sa2 = d$sa2, pi = d$pi0, sigma2 = d$sigma2_init,
    intercept = FALSE, standardize = FALSE,
    update.sigma2 = FALSE, update.pi = FALSE,
    max.iter = 100, verbose = FALSE
  )

  ss <- derive_summary_stats(d$X, d$y)
  fit_rss <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R,
    var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 100,
    update_w0 = FALSE, update_sigma = FALSE
  )

  # Should match to near-machine precision (ignore dim attributes from Armadillo)
  expect_equal(c(fit_rss$beta), c(fit_ind$beta), tolerance = 1e-10)
  expect_equal(c(fit_rss$sigma2), c(fit_ind$sigma2), tolerance = 1e-10)
  expect_equal(c(fit_rss$pi), c(fit_ind$pi), tolerance = 1e-10)
})

test_that("mr.ash and mr.ash.rss agree with sigma2 updates enabled", {
  d <- setup_mr_ash_test(n = 100, p = 50, k = 5, seed = 42)

  fit_ind <- mr.ash(d$X, d$y,
    sa2 = d$sa2, pi = d$pi0, sigma2 = d$sigma2_init,
    intercept = FALSE, standardize = FALSE,
    update.sigma2 = TRUE, update.pi = FALSE,
    max.iter = 200, verbose = FALSE
  )

  ss <- derive_summary_stats(d$X, d$y)
  fit_rss <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R,
    var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 200,
    update_w0 = FALSE, update_sigma = TRUE
  )

  expect_equal(c(fit_rss$beta), c(fit_ind$beta), tolerance = 1e-3)
  expect_equal(c(fit_rss$sigma2), c(fit_ind$sigma2), tolerance = 1e-4)
})

test_that("mr.ash and mr.ash.rss agree with full EM (sigma + pi updates)", {
  d <- setup_mr_ash_test(n = 100, p = 50, k = 5, seed = 42)

  fit_ind <- mr.ash(d$X, d$y,
    sa2 = d$sa2, pi = d$pi0, sigma2 = d$sigma2_init,
    intercept = FALSE, standardize = FALSE,
    update.sigma2 = TRUE, update.pi = TRUE,
    max.iter = 200, verbose = FALSE
  )

  ss <- derive_summary_stats(d$X, d$y)
  fit_rss <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R,
    var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 200,
    update_w0 = TRUE, update_sigma = TRUE
  )

  expect_equal(c(fit_rss$beta), c(fit_ind$beta), tolerance = 1e-3)
  expect_equal(c(fit_rss$sigma2), c(fit_ind$sigma2), tolerance = 1e-3)
  expect_equal(c(fit_rss$pi), c(fit_ind$pi), tolerance = 1e-2)
})

test_that("mr.ash.rss output format matches mr.ash", {
  d <- setup_mr_ash_test(n = 80, p = 20, k = 3, seed = 123)

  ss <- derive_summary_stats(d$X, d$y)
  fit_rss <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R,
    var_y = ss$var_y, n = ss$n,
    sigma2_e = 1.0, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 100,
    update_w0 = FALSE, update_sigma = FALSE
  )

  # Check that mr.ash-compatible fields exist and have correct types
  expect_true(is.numeric(fit_rss$beta))
  expect_true(is.numeric(fit_rss$sigma2))
  expect_true(is.numeric(fit_rss$pi))
  expect_true(is.integer(fit_rss$iter))
  expect_true(is.numeric(fit_rss$varobj))

  # Check dimensions
  expect_length(fit_rss$beta, d$p)
  expect_length(fit_rss$sigma2, 1)
  expect_length(fit_rss$pi, d$K)
  expect_true(fit_rss$iter > 0)
  expect_true(length(fit_rss$varobj) > 0)
  expect_true(length(fit_rss$varobj) <= 100)

  # Original RSS-specific fields also present
  expect_true(!is.null(fit_rss$mu1))
  expect_true(!is.null(fit_rss$sigma2_1))
  expect_true(!is.null(fit_rss$w1))
  expect_true(!is.null(fit_rss$sigma2_e))
  expect_true(!is.null(fit_rss$w0))
})

test_that("mr.ash and mr.ash.rss agree on different data sizes", {
  # Test with a wider range of n, p combinations
  for (params in list(
    list(n = 80, p = 20, k = 3, seed = 100),
    list(n = 200, p = 30, k = 5, seed = 200)
  )) {
    d <- setup_mr_ash_test(
      n = params$n, p = params$p,
      k = params$k, seed = params$seed
    )

    fit_ind <- mr.ash(d$X, d$y,
      sa2 = d$sa2, pi = d$pi0, sigma2 = d$sigma2_init,
      intercept = FALSE, standardize = FALSE,
      update.sigma2 = FALSE, update.pi = FALSE,
      max.iter = 50, verbose = FALSE
    )

    ss <- derive_summary_stats(d$X, d$y)
    fit_rss <- mr.ash.rss(
      bhat = ss$bhat, shat = ss$shat, R = ss$R,
      var_y = ss$var_y, n = ss$n,
      sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
      tol = 1e-4, max_iter = 50,
      update_w0 = FALSE, update_sigma = FALSE
    )

    expect_equal(c(fit_rss$beta), c(fit_ind$beta), tolerance = 1e-10,
      label = sprintf("beta (n=%d, p=%d)", params$n, params$p)
    )
  }
})
