# Test that susie_rss gives equivalent results for raw, centered,
# and standardized X inputs.

context("X centering equivalence")

test_that("Full-rank X: raw, centered, standardized give same results", {
  set.seed(1)
  n <- 200
  p <- 50

  # Generate X with non-zero column means (raw)
  X_raw <- matrix(rnorm(n * p, mean = 5, sd = 2), n, p)
  X_raw[, 1:3] <- X_raw[, 1:3] + 10  # make some columns have large means

  # Centered and standardized versions
  X_centered <- scale(X_raw, center = TRUE, scale = FALSE)
  X_standardized <- scale(X_raw, center = TRUE, scale = TRUE)

  # Generate z-scores from centered X
  beta_true <- rep(0, p)
  beta_true[c(5, 15, 30)] <- c(0.5, -0.3, 0.4)
  y <- X_centered %*% beta_true + rnorm(n)
  z <- as.vector(sqrt(n) * cor(X_centered, y))

  # Fit with all three forms of X (full-rank path: nrow >= ncol)
  fit_raw  <- susie_rss(z = z, X = X_raw, n = n, max_iter = 50)
  fit_cent <- susie_rss(z = z, X = X_centered, n = n, max_iter = 50)
  fit_std  <- susie_rss(z = z, X = X_standardized, n = n, max_iter = 50)

  # All should produce identical results (full-rank path uses safe_cor
  # which centers internally, but we also center before calling safe_cor)
  expect_equal(fit_raw$elbo, fit_cent$elbo, tolerance = 1e-10)
  expect_equal(fit_raw$elbo, fit_std$elbo, tolerance = 1e-10)
  expect_equal(fit_raw$pip, fit_cent$pip, tolerance = 1e-10)
  expect_equal(fit_raw$pip, fit_std$pip, tolerance = 1e-10)
})

test_that("Low-rank X: raw and centered give identical results", {
  set.seed(2)
  n <- 200
  p <- 500
  B <- 100  # sketch size < p, triggers low-rank path

  # Generate full X and z-scores
  X_full <- matrix(rnorm(n * p), n, p)
  beta_true <- rep(0, p)
  beta_true[c(10, 50, 200)] <- c(0.5, -0.3, 0.4)
  y <- X_full %*% beta_true + rnorm(n)
  z <- as.vector(sqrt(n) * cor(X_full, y))

  # Create a sketch matrix (B x p, B < p)
  S <- matrix(rnorm(B * n) / sqrt(B), B, n)
  X_sketch <- S %*% X_full  # B x p sketch

  # Add offset to create "raw" version with non-zero column means
  X_sketch_raw <- X_sketch + 10
  # Manually center (avoid scale() attributes)
  X_sketch_centered <- X_sketch_raw - rep(colMeans(X_sketch_raw), each = B)

  # Fit with both forms (low-rank path: nrow < ncol)
  fit_raw  <- susie_rss(z = z, X = X_sketch_raw, n = n, L = 5, max_iter = 50)
  fit_cent <- susie_rss(z = z, X = X_sketch_centered, n = n, L = 5, max_iter = 50)

  # Raw and centered should give identical results
  expect_equal(fit_raw$elbo, fit_cent$elbo, tolerance = 1e-10)
  expect_equal(fit_raw$pip, fit_cent$pip, tolerance = 1e-10)
})

test_that("Low-rank X: raw vs centered give same results (no n)", {
  set.seed(3)
  p <- 300
  B <- 80

  # Generate sketch matrix with non-zero means
  X_raw <- matrix(rnorm(B * p, mean = 3), B, p)
  X_centered <- X_raw - rep(colMeans(X_raw), each = B)

  z <- rnorm(p)
  z[c(5, 100)] <- c(4, -3.5)

  fit_raw  <- susie_rss(z = z, X = X_raw, L = 5, max_iter = 30)
  fit_cent <- susie_rss(z = z, X = X_centered, L = 5, max_iter = 30)

  expect_equal(fit_raw$elbo, fit_cent$elbo, tolerance = 1e-10)
  expect_equal(fit_raw$pip, fit_cent$pip, tolerance = 1e-10)
})

test_that("Low-rank X with lambda: raw and centered give same results", {
  set.seed(4)
  p <- 200
  B <- 50

  X_raw <- matrix(rnorm(B * p, mean = 2, sd = 3), B, p)
  X_centered <- X_raw - rep(colMeans(X_raw), each = B)

  z <- rnorm(p)
  z[c(10, 80)] <- c(5, -4)

  fit_raw  <- susie_rss(z = z, X = X_raw, lambda = 0.1, L = 5, max_iter = 30)
  fit_cent <- susie_rss(z = z, X = X_centered, lambda = 0.1, L = 5, max_iter = 30)

  expect_equal(fit_raw$elbo, fit_cent$elbo, tolerance = 1e-10)
  expect_equal(fit_raw$pip, fit_cent$pip, tolerance = 1e-10)
})

test_that("Low-rank X: large offset does not break fitting", {
  # Verify that adding a large constant offset to X columns does not
  # break the fitting (because centering removes it).
  set.seed(5)
  p <- 200
  B <- 80

  X_centered <- matrix(rnorm(B * p), B, p)
  X_offset <- X_centered + 1000  # huge offset

  z <- rnorm(p)
  z[c(10, 80)] <- c(5, -4)

  # Both should give identical results
  fit_cent <- susie_rss(z = z, X = X_centered, n = 500, L = 5, max_iter = 50)
  fit_off  <- susie_rss(z = z, X = X_offset, n = 500, L = 5, max_iter = 50)

  # Tolerance accounts for floating-point cancellation when subtracting
  # the large offset (1000 * machine_eps ≈ 2e-13, amplified by iterations)
  expect_equal(fit_cent$elbo, fit_off$elbo, tolerance = 1e-6)
  expect_equal(fit_cent$pip, fit_off$pip, tolerance = 1e-6)
})
