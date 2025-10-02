# Source helper functions
source(file.path("..", "helpers", "helper_reference.R"), local = TRUE)

context("susie_rss reference comparison")

# =============================================================================
# REFERENCE TESTS FOR susie_rss()
# =============================================================================
#
# These tests compare the new susie_rss() implementation against the reference
# package susie_rss() from stephenslab/susieR@1f9166c
#
# Tests standard RSS (no lambda parameter)

# =============================================================================
# Standard RSS Tests
# =============================================================================

test_that("susie_rss() matches reference with default parameters", {
  skip_if_no_reference()

  set.seed(1)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute z-scores and R using standard approach
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10)
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with bhat and shat", {
  skip_if_no_reference()

  set.seed(2)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute bhat, shat and R
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  bhat <- ss$betahat
  shat <- ss$sebetahat

  args <- list(bhat = bhat, shat = shat, R = R, n = n, L = 10)
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference without n provided", {
  skip_if_no_reference()

  set.seed(3)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute z-scores and R
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Note: n not provided - uses large n approximation
  args <- list(z = z, R = R, L = 10)
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with estimate_residual_variance=TRUE", {
  skip_if_no_reference()

  set.seed(4)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute z-scores and R
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_residual_variance = TRUE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with estimate_prior_variance=FALSE", {
  skip_if_no_reference()

  set.seed(5)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute z-scores and R
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_prior_variance = FALSE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with different L values", {
  skip_if_no_reference()

  set.seed(6)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute z-scores and R
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Test L=1
  args1 <- list(z = z, R = R, n = n, L = 1)
  compare_to_reference("susie_rss", args1, tolerance = 1e-5)

  # Test L=5
  args5 <- list(z = z, R = R, n = n, L = 5)
  compare_to_reference("susie_rss", args5, tolerance = 1e-5)

  # Test L=20
  args20 <- list(z = z, R = R, n = n, L = 20)
  compare_to_reference("susie_rss", args20, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with different estimate_prior_method", {
  skip_if_no_reference()

  set.seed(7)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute z-scores and R
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Test "EM"
  args_em <- list(z = z, R = R, n = n, L = 10, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args_em, tolerance = 1e-5)

  # Test "simple"
  args_simple <- list(z = z, R = R, n = n, L = 10, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args_simple, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with prior_weights", {
  skip_if_no_reference()

  set.seed(8)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute z-scores and R
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Uniform prior weights
  prior_weights <- rep(1/p, p)
  args <- list(z = z, R = R, n = n, L = 10, prior_weights = prior_weights)
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with var_y", {
  skip_if_no_reference()

  set.seed(9)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute bhat, shat and R
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  bhat <- ss$betahat
  shat <- ss$sebetahat
  var_y <- var(y)

  args <- list(bhat = bhat, shat = shat, R = R, n = n, L = 10, var_y = var_y)
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with different coverage", {
  skip_if_no_reference()

  set.seed(10)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute z-scores and R
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, coverage = 0.99)
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})
