context("X centering equivalence")

# ---- full-rank path (nrow >= ncol) ------------------------------------------

test_that("Full-rank X: raw, centered, and standardized give identical results", {
  set.seed(1)
  n <- 200
  p <- 50

  X_raw <- matrix(rnorm(n * p, mean = 5, sd = 2), n, p)
  X_raw[, 1:3] <- X_raw[, 1:3] + 10

  X_centered     <- scale(X_raw, center = TRUE,  scale = FALSE)
  X_standardized <- scale(X_raw, center = TRUE,  scale = TRUE)

  beta_true <- rep(0, p)
  beta_true[c(5, 15, 30)] <- c(0.5, -0.3, 0.4)
  y <- X_centered %*% beta_true + rnorm(n)
  z <- as.vector(sqrt(n) * cor(X_centered, y))

  fit_raw  <- suppressWarnings(susie_rss(z = z, X = X_raw,          n = n, max_iter = 50))
  fit_cent <- suppressWarnings(susie_rss(z = z, X = X_centered,     n = n, max_iter = 50))
  fit_std  <- suppressWarnings(susie_rss(z = z, X = X_standardized, n = n, max_iter = 50))

  expect_equal(fit_raw$elbo, fit_cent$elbo, tolerance = 1e-8)
  expect_equal(fit_raw$elbo, fit_std$elbo,  tolerance = 1e-8)
  expect_equal(fit_raw$pip,  fit_cent$pip,  tolerance = 1e-8)
  expect_equal(fit_raw$pip,  fit_std$pip,   tolerance = 1e-8)
})

# ---- low-rank path (nrow < ncol) --------------------------------------------

test_that("Low-rank X with n: raw and centered give identical results", {
  set.seed(2)
  n <- 200
  p <- 500
  B <- 100

  X_full <- matrix(rnorm(n * p), n, p)
  beta_true <- rep(0, p)
  beta_true[c(10, 50, 200)] <- c(0.5, -0.3, 0.4)
  y <- X_full %*% beta_true + rnorm(n)
  z <- as.vector(sqrt(n) * cor(X_full, y))

  S <- matrix(rnorm(B * n) / sqrt(B), B, n)
  X_ref_raw      <- S %*% X_full + 10
  X_ref_centered <- X_ref_raw - rep(colMeans(X_ref_raw), each = B)

  fit_raw  <- suppressWarnings(susie_rss(z = z, X = X_ref_raw,      n = n, L = 5, max_iter = 50))
  fit_cent <- suppressWarnings(susie_rss(z = z, X = X_ref_centered, n = n, L = 5, max_iter = 50))

  expect_equal(fit_raw$elbo, fit_cent$elbo, tolerance = 1e-8)
  expect_equal(fit_raw$pip,  fit_cent$pip,  tolerance = 1e-8)
})

test_that("Low-rank X without n: raw and centered give identical results", {
  set.seed(3)
  p <- 300
  B <- 80

  X_raw      <- matrix(rnorm(B * p, mean = 3), B, p)
  X_centered <- X_raw - rep(colMeans(X_raw), each = B)

  z <- rnorm(p)
  z[c(5, 100)] <- c(4, -3.5)

  fit_raw  <- suppressWarnings(susie_rss(z = z, X = X_raw,      L = 5, max_iter = 30))
  fit_cent <- suppressWarnings(susie_rss(z = z, X = X_centered, L = 5, max_iter = 30))

  expect_equal(fit_raw$elbo, fit_cent$elbo, tolerance = 1e-8)
  expect_equal(fit_raw$pip,  fit_cent$pip,  tolerance = 1e-8)
})

test_that("Low-rank X with lambda: raw and centered give identical results", {
  set.seed(4)
  p <- 200
  B <- 50

  X_raw      <- matrix(rnorm(B * p, mean = 2, sd = 3), B, p)
  X_centered <- X_raw - rep(colMeans(X_raw), each = B)

  z <- rnorm(p)
  z[c(10, 80)] <- c(5, -4)

  fit_raw  <- suppressWarnings(susie_rss_lambda(z = z, X = X_raw,      lambda = 0.1, L = 5, max_iter = 30))
  fit_cent <- suppressWarnings(susie_rss_lambda(z = z, X = X_centered, lambda = 0.1, L = 5, max_iter = 30))

  expect_equal(fit_raw$elbo, fit_cent$elbo, tolerance = 1e-8)
  expect_equal(fit_raw$pip,  fit_cent$pip,  tolerance = 1e-8)
})

test_that("Low-rank X: large column offset does not affect fitted results", {
  set.seed(5)
  p <- 200
  B <- 80

  X_centered <- matrix(rnorm(B * p), B, p)
  X_offset   <- X_centered + 1000

  z <- rnorm(p)
  z[c(10, 80)] <- c(5, -4)

  fit_cent <- suppressWarnings(susie_rss(z = z, X = X_centered, n = 500, L = 5, max_iter = 50))
  fit_off  <- suppressWarnings(susie_rss(z = z, X = X_offset,   n = 500, L = 5, max_iter = 50))

  # Tolerance accounts for floating-point cancellation when subtracting 1000
  expect_equal(fit_cent$elbo, fit_off$elbo, tolerance = 1e-6)
  expect_equal(fit_cent$pip,  fit_off$pip,  tolerance = 1e-6)
})
