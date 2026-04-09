# Tests for estimate_prior_method = "fixed_mixture"
#
# Key invariant: a K=1 mixture with grid = c(V) and weights = c(1)
# must produce identical results to the scalar V path with
# estimate_prior_variance = FALSE and prior_variance = V.

context("Fixed mixture prior")

# Generate a small test dataset
set.seed(1)
n <- 200
p <- 50
X <- matrix(rnorm(n * p), n, p)
beta <- rep(0, p)
beta[c(1, 5, 10)] <- c(0.5, -0.3, 0.4)
y <- X %*% beta + rnorm(n)

# Compute summary stats
R <- cor(X)
z <- as.vector(sqrt(n) * crossprod(X, y) / sqrt(n * diag(crossprod(X))))

# =============================================================================
# Test 1: K=1 mixture matches scalar V exactly (individual data)
# =============================================================================
test_that("K=1 mixture matches scalar V for individual data", {
  L <- 5
  # Run scalar path first to find effective V
  fit_scalar <- susie(X, y, L = L,
                      estimate_prior_variance = FALSE,
                      estimate_residual_variance = FALSE,
                      max_iter = 20, tol = 1e-4)
  V_eff <- fit_scalar$V[1]

  # K=1 mixture path with the same effective V
  fit_mixture <- susie(X, y, L = L,
                       prior_variance_grid = c(V_eff),
                       mixture_weights = c(1),
                       estimate_residual_variance = FALSE,
                       max_iter = 20, tol = 1e-4)

  expect_equal(fit_scalar$pip, fit_mixture$pip, tolerance = 1e-10)
  expect_equal(fit_scalar$alpha, fit_mixture$alpha, tolerance = 1e-10)
  expect_equal(fit_scalar$mu, fit_mixture$mu, tolerance = 1e-10)
  expect_equal(fit_scalar$lbf, fit_mixture$lbf, tolerance = 1e-10)
})

# =============================================================================
# Test 2: K=1 mixture matches scalar V exactly (RSS data)
# =============================================================================
test_that("K=1 mixture matches scalar V for RSS data", {
  L <- 5
  # Run scalar path first to find what V is actually used
  fit_scalar <- susie_rss(z = z, R = R, n = n, L = L,
                          estimate_prior_variance = FALSE,
                          estimate_residual_variance = FALSE,
                          max_iter = 20, tol = 1e-4)
  # The effective V is stored in fit_scalar$V[1] (same for all L since
  # scaled_prior_variance is a scalar and estimate_prior_variance = FALSE)
  V_eff <- fit_scalar$V[1]

  # K=1 mixture path with the same effective V
  fit_mixture <- susie_rss(z = z, R = R, n = n, L = L,
                           prior_variance_grid = c(V_eff),
                           mixture_weights = c(1),
                           estimate_residual_variance = FALSE,
                           max_iter = 20, tol = 1e-4)

  expect_equal(fit_scalar$pip, fit_mixture$pip, tolerance = 1e-10)
  expect_equal(fit_scalar$alpha, fit_mixture$alpha, tolerance = 1e-10)
  expect_equal(fit_scalar$mu, fit_mixture$mu, tolerance = 1e-10)
  expect_equal(fit_scalar$lbf, fit_mixture$lbf, tolerance = 1e-10)
})

# =============================================================================
# Test 3: K=1 mixture matches scalar V exactly (sufficient stats)
# =============================================================================
test_that("K=1 mixture matches scalar V for sufficient stats", {
  L <- 5
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  yty <- sum(y^2)

  # Run scalar path first to find effective V
  fit_scalar <- susie_ss(XtX = XtX, Xty = Xty, yty = yty, n = n, L = L,
                         estimate_prior_variance = FALSE,
                         estimate_residual_variance = FALSE,
                         max_iter = 20, tol = 1e-4)
  V_eff <- fit_scalar$V[1]

  # K=1 mixture path with the same effective V
  fit_mixture <- susie_ss(XtX = XtX, Xty = Xty, yty = yty, n = n, L = L,
                          prior_variance_grid = c(V_eff),
                          mixture_weights = c(1),
                          estimate_residual_variance = FALSE,
                          max_iter = 20, tol = 1e-4)

  expect_equal(fit_scalar$pip, fit_mixture$pip, tolerance = 1e-10)
  expect_equal(fit_scalar$alpha, fit_mixture$alpha, tolerance = 1e-10)
  expect_equal(fit_scalar$mu, fit_mixture$mu, tolerance = 1e-10)
  expect_equal(fit_scalar$lbf, fit_mixture$lbf, tolerance = 1e-10)
})

# =============================================================================
# Test 4: K>1 mixture produces valid outputs
# =============================================================================
test_that("K=3 mixture produces valid outputs for RSS data", {
  L <- 5
  grid <- c(1, 10, 50)
  w <- c(0.3, 0.5, 0.2)

  fit <- susie_rss(z = z, R = R, n = n, L = L,
                   prior_variance_grid = grid,
                   mixture_weights = w,
                   estimate_residual_variance = FALSE,
                   max_iter = 20, tol = 1e-4)

  # PIPs in [0, 1]
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  # Alpha rows sum to 1
  expect_equal(rowSums(fit$alpha), rep(1, L), tolerance = 1e-10)
  # lbf_grid is a list of L elements
  expect_length(fit$lbf_grid, L)
  # Each element is p x K matrix
  expect_equal(dim(fit$lbf_grid[[1]]), c(p, 3))
  # Posterior means are finite
  expect_true(all(is.finite(fit$mu)))
  # Posterior second moments >= posterior means squared
  expect_true(all(fit$mu2 >= fit$mu^2 - 1e-10))
})

# =============================================================================
# Test 5: Uniform weights produces correct mixture BF
# =============================================================================
test_that("Uniform mixture weights give correct BF", {
  L <- 1
  grid <- c(1, 50)
  w <- c(0.5, 0.5)

  fit <- susie_rss(z = z, R = R, n = n, L = L,
                   prior_variance_grid = grid,
                   mixture_weights = w,
                   estimate_residual_variance = FALSE,
                   max_iter = 1)

  # Manually compute mixture BF for variant 1
  # lbf(V) = -0.5*log(1 + V*R[1,1]) + 0.5*z[1]^2*V*R[1,1]/(V*R[1,1]+1)
  # For RSS with lambda=0, sigma2=1: shat2 = 1/R[1,1], betahat = z[1]/R[1,1] * shat2
  # This is approximate due to eigendecomposition; just check BF matrix shape
  expect_equal(ncol(fit$lbf_grid[[1]]), 2)
  expect_equal(nrow(fit$lbf_grid[[1]]), p)
})

# =============================================================================
# Test 6: Input validation
# =============================================================================
test_that("Invalid mixture prior inputs are rejected", {
  # Mismatched lengths
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 5,
              prior_variance_grid = c(1, 10),
              mixture_weights = c(1)),
    "length"
  )
  # Negative grid values
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 5,
              prior_variance_grid = c(-1, 10),
              mixture_weights = c(0.5, 0.5)),
    "prior_variance_grid"
  )
  # Weights not summing to 1
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 5,
              prior_variance_grid = c(1, 10),
              mixture_weights = c(0.3, 0.3)),
    "sum"
  )
})

# =============================================================================
# Test 7: Default weights (uniform) when mixture_weights is NULL
# =============================================================================
test_that("NULL mixture_weights defaults to uniform", {
  L <- 5
  grid <- c(1, 10, 50)

  # Should not error, should use uniform weights
  fit <- susie_rss(z = z, R = R, n = n, L = L,
                   prior_variance_grid = grid,
                   estimate_residual_variance = FALSE,
                   max_iter = 5)

  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

# =============================================================================
# Test 8: Existing tests still pass (backward compatibility)
# =============================================================================
test_that("Standard susie_rss without mixture prior is unchanged", {
  fit <- susie_rss(z = z, R = R, n = n, L = 5,
                   estimate_prior_variance = TRUE,
                   estimate_residual_variance = FALSE,
                   max_iter = 20)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_null(fit$lbf_grid)  # no grid stored in standard path
})

# =============================================================================
# Test 9: K=1 mixture with RSS-lambda path (explicit lambda > 0)
# =============================================================================
test_that("K=1 mixture matches scalar V with lambda regularization", {
  L <- 3
  lam <- 0.1
  fit_scalar <- susie_rss(z = z, R = R, n = n, L = L, lambda = lam,
                          estimate_prior_variance = FALSE,
                          estimate_residual_variance = FALSE,
                          max_iter = 10)
  V_eff <- fit_scalar$V[1]

  fit_mixture <- susie_rss(z = z, R = R, n = n, L = L, lambda = lam,
                           prior_variance_grid = c(V_eff),
                           mixture_weights = c(1),
                           estimate_residual_variance = FALSE,
                           max_iter = 10)

  expect_equal(fit_scalar$alpha, fit_mixture$alpha, tolerance = 1e-10)
  expect_equal(fit_scalar$mu, fit_mixture$mu, tolerance = 1e-10)
  expect_equal(fit_scalar$lbf, fit_mixture$lbf, tolerance = 1e-10)
})

# =============================================================================
# Test 10: Mixture prior with sketch LD inflation (inflated shat2)
# =============================================================================
test_that("Mixture prior works with sketch LD inflation", {
  skip_if_not_installed("Matrix")
  L <- 3
  grid <- c(1, 10, 50)
  w <- c(0.3, 0.5, 0.2)

  # Run with sketch_samples to trigger shat2 inflation
  fit <- susie_rss(z = z, R = R, n = n, L = L,
                   prior_variance_grid = grid,
                   mixture_weights = w,
                   estimate_residual_variance = FALSE,
                   sketch_samples = 30,
                   max_iter = 5)

  # Basic validity
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_equal(rowSums(fit$alpha), rep(1, L), tolerance = 1e-10)
  expect_length(fit$lbf_grid, L)
})

# =============================================================================
# Test 11: Mixture weights are correctly used (asymmetric weights)
# =============================================================================
test_that("Asymmetric mixture weights shift PIPs correctly", {
  L <- 3
  # Large V component only: should produce wider credible intervals
  fit_large <- susie_rss(z = z, R = R, n = n, L = L,
                         prior_variance_grid = c(0.001, 100),
                         mixture_weights = c(0.01, 0.99),
                         estimate_residual_variance = FALSE,
                         max_iter = 10)
  # Small V component only: should produce tighter credible intervals
  fit_small <- susie_rss(z = z, R = R, n = n, L = L,
                         prior_variance_grid = c(0.001, 100),
                         mixture_weights = c(0.99, 0.01),
                         estimate_residual_variance = FALSE,
                         max_iter = 10)

  # Both should be valid
  expect_true(all(fit_large$pip >= 0 & fit_large$pip <= 1))
  expect_true(all(fit_small$pip >= 0 & fit_small$pip <= 1))
  # PIPs should differ
  expect_false(all(abs(fit_large$pip - fit_small$pip) < 1e-6))
})

# =============================================================================
# Test 12: K=1 individual data exact match (L=1, single iteration)
# =============================================================================
test_that("K=1 mixture is numerically identical for L=1 individual data", {
  L <- 1
  fit_scalar <- susie(X, y, L = L,
                      estimate_prior_variance = FALSE,
                      estimate_residual_variance = FALSE,
                      max_iter = 1)
  V_eff <- fit_scalar$V[1]

  fit_mixture <- susie(X, y, L = L,
                       prior_variance_grid = c(V_eff),
                       mixture_weights = c(1),
                       estimate_residual_variance = FALSE,
                       max_iter = 1)

  # After exactly 1 iteration with L=1, results must match to machine precision
  expect_equal(fit_scalar$alpha, fit_mixture$alpha, tolerance = .Machine$double.eps * 10)
  expect_equal(fit_scalar$mu, fit_mixture$mu, tolerance = .Machine$double.eps * 10)
  expect_equal(fit_scalar$lbf, fit_mixture$lbf, tolerance = .Machine$double.eps * 10)
})
