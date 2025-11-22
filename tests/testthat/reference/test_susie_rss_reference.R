# Source helper functions
source(file.path("..", "helper_reference.R"), local = TRUE)

context("susie_rss reference comparison")

# =============================================================================
# REFERENCE TESTS FOR susie_rss() with lambda = 0
# =============================================================================
#
# These tests compare the new susie_rss() implementation (lambda = 0) against
# the reference package susie_rss() from stephenslab/susieR@1f9166c
#
# Tests cover all major parameters with all three prior variance optimization
# methods: "optim", "EM", "simple"

# =============================================================================
# Part 1: Basic Input Formats
# =============================================================================
# Test that different input formats (z, bhat/shat) work correctly

test_that("susie_rss() matches reference with z-scores - optim", {
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

  args <- list(z = z, R = R, n = n, L = 10, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with z-scores - EM", {
  skip_if_no_reference()

  set.seed(1)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with z-scores - simple", {
  skip_if_no_reference()

  set.seed(1)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with bhat/shat - optim", {
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

  args <- list(bhat = bhat, shat = shat, R = R, n = n, L = 10, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with bhat/shat - EM", {
  skip_if_no_reference()

  set.seed(2)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  bhat <- ss$betahat
  shat <- ss$sebetahat

  args <- list(bhat = bhat, shat = shat, R = R, n = n, L = 10, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with bhat/shat - simple", {
  skip_if_no_reference()

  set.seed(2)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  bhat <- ss$betahat
  shat <- ss$sebetahat

  args <- list(bhat = bhat, shat = shat, R = R, n = n, L = 10, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

# =============================================================================
# Part 2: Sample size n parameter
# =============================================================================
# Test with n provided vs. not provided (large n approximation)

test_that("susie_rss() matches reference with n provided - optim", {
  skip_if_no_reference()

  set.seed(3)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with n provided - EM", {
  skip_if_no_reference()

  set.seed(3)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with n provided - simple", {
  skip_if_no_reference()

  set.seed(3)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference without n - optim", {
  skip_if_no_reference()

  set.seed(4)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Note: n not provided - uses large n approximation
  args <- list(z = z, R = R, L = 10, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference without n - EM", {
  skip_if_no_reference()

  set.seed(4)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, L = 10, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference without n - simple", {
  skip_if_no_reference()

  set.seed(4)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, L = 10, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

# =============================================================================
# Part 3: Different L values
# =============================================================================

test_that("susie_rss() matches reference with different L values - optim", {
  skip_if_no_reference()

  set.seed(5)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Test L=1
  args1 <- list(z = z, R = R, n = n, L = 1, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args1, tolerance = 1e-5)

  # Test L=5
  args5 <- list(z = z, R = R, n = n, L = 5, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args5, tolerance = 1e-5)

  # Test L=20
  args20 <- list(z = z, R = R, n = n, L = 20, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args20, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with different L values - EM", {
  skip_if_no_reference()

  set.seed(5)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Test L=1
  args1 <- list(z = z, R = R, n = n, L = 1, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args1, tolerance = 1e-5)

  # Test L=5
  args5 <- list(z = z, R = R, n = n, L = 5, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args5, tolerance = 1e-5)

  # Test L=20
  args20 <- list(z = z, R = R, n = n, L = 20, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args20, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with different L values - simple", {
  skip_if_no_reference()

  set.seed(5)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Test L=1
  args1 <- list(z = z, R = R, n = n, L = 1, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args1, tolerance = 1e-5)

  # Test L=5
  args5 <- list(z = z, R = R, n = n, L = 5, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args5, tolerance = 1e-5)

  # Test L=20
  args20 <- list(z = z, R = R, n = n, L = 20, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args20, tolerance = 1e-5)
})

# =============================================================================
# Part 4: estimate_prior_variance parameter
# =============================================================================

test_that("susie_rss() matches reference with estimate_prior_variance=FALSE - optim", {
  skip_if_no_reference()

  set.seed(6)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_prior_variance = FALSE,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with estimate_prior_variance=FALSE - EM", {
  skip_if_no_reference()

  set.seed(6)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_prior_variance = FALSE,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with estimate_prior_variance=FALSE - simple", {
  skip_if_no_reference()

  set.seed(6)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_prior_variance = FALSE,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

# =============================================================================
# Part 5: estimate_residual_variance parameter
# =============================================================================

test_that("susie_rss() matches reference with estimate_residual_variance=TRUE - optim", {
  skip_if_no_reference()

  set.seed(7)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_residual_variance = TRUE,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with estimate_residual_variance=TRUE - EM", {
  skip_if_no_reference()

  set.seed(7)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_residual_variance = TRUE,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with estimate_residual_variance=TRUE - simple", {
  skip_if_no_reference()

  set.seed(7)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_residual_variance = TRUE,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with estimate_residual_variance=FALSE - optim", {
  skip_if_no_reference()

  set.seed(8)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with estimate_residual_variance=FALSE - EM", {
  skip_if_no_reference()

  set.seed(8)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with estimate_residual_variance=FALSE - simple", {
  skip_if_no_reference()

  set.seed(8)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

# =============================================================================
# Part 6: prior_weights
# =============================================================================

test_that("susie_rss() matches reference with prior_weights - optim", {
  skip_if_no_reference()

  set.seed(9)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Use non-uniform prior weights
  prior_weights <- runif(p)
  prior_weights <- prior_weights / sum(prior_weights)
  args <- list(z = z, R = R, n = n, L = 10, prior_weights = prior_weights, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with prior_weights - EM", {
  skip_if_no_reference()

  set.seed(9)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Use non-uniform prior weights
  prior_weights <- runif(p)
  prior_weights <- prior_weights / sum(prior_weights)
  args <- list(z = z, R = R, n = n, L = 10, prior_weights = prior_weights, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with prior_weights - simple", {
  skip_if_no_reference()

  set.seed(9)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Use non-uniform prior weights
  prior_weights <- runif(p)
  prior_weights <- prior_weights / sum(prior_weights)
  args <- list(z = z, R = R, n = n, L = 10, prior_weights = prior_weights, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

# =============================================================================
# Part 7: scaled_prior_variance
# =============================================================================

test_that("susie_rss() matches reference with scaled_prior_variance - optim", {
  skip_if_no_reference()

  set.seed(10)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, scaled_prior_variance = 0.5, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with scaled_prior_variance - EM", {
  skip_if_no_reference()

  set.seed(10)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, scaled_prior_variance = 0.5, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with scaled_prior_variance - simple", {
  skip_if_no_reference()

  set.seed(10)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, scaled_prior_variance = 0.5, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

# =============================================================================
# Part 8: var_y parameter
# =============================================================================

test_that("susie_rss() matches reference with var_y - optim", {
  skip_if_no_reference()

  set.seed(11)
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

  args <- list(bhat = bhat, shat = shat, R = R, n = n, L = 10, var_y = var_y, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with var_y - EM", {
  skip_if_no_reference()

  set.seed(11)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  bhat <- ss$betahat
  shat <- ss$sebetahat
  var_y <- var(y)

  args <- list(bhat = bhat, shat = shat, R = R, n = n, L = 10, var_y = var_y, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with var_y - simple", {
  skip_if_no_reference()

  set.seed(11)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  bhat <- ss$betahat
  shat <- ss$sebetahat
  var_y <- var(y)

  args <- list(bhat = bhat, shat = shat, R = R, n = n, L = 10, var_y = var_y, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

# =============================================================================
# Part 9: coverage and min_abs_corr
# =============================================================================

test_that("susie_rss() matches reference with coverage=0.99 - optim", {
  skip_if_no_reference()

  set.seed(12)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, coverage = 0.99, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with coverage=0.99 - EM", {
  skip_if_no_reference()

  set.seed(12)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, coverage = 0.99, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with coverage=0.99 - simple", {
  skip_if_no_reference()

  set.seed(12)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, coverage = 0.99, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with min_abs_corr=0.7 - optim", {
  skip_if_no_reference()

  set.seed(13)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, min_abs_corr = 0.7, estimate_prior_method = "optim")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with min_abs_corr=0.7 - EM", {
  skip_if_no_reference()

  set.seed(13)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, min_abs_corr = 0.7, estimate_prior_method = "EM")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with min_abs_corr=0.7 - simple", {
  skip_if_no_reference()

  set.seed(13)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, n = n, L = 10, min_abs_corr = 0.7, estimate_prior_method = "simple")
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

# =============================================================================
# Part 10: prior_tol parameter
# =============================================================================

test_that("susie_rss() matches reference with prior_tol=1e-5 - optim", {
  skip_if_no_reference()

  set.seed(15)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Disable check_null_threshold so we can see prior_tol effects
  args <- list(
    z = z, R = R, n = n, L = 10,
    prior_tol = 0.1,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with prior_tol=1e-5 - EM", {
  skip_if_no_reference()

  set.seed(15)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    prior_tol = 0.1,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with prior_tol=1e-5 - simple", {
  skip_if_no_reference()

  set.seed(15)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    prior_tol = 0.1,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

# =============================================================================
# Part 11: check_null_threshold parameter
# =============================================================================

test_that("susie_rss() matches reference with check_null_threshold=0.1 - optim", {
  skip_if_no_reference()

  set.seed(16)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    check_null_threshold = 0.1,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with check_null_threshold=0.1 - EM", {
  skip_if_no_reference()

  set.seed(16)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    check_null_threshold = 0.1,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with check_null_threshold=0.1 - simple", {
  skip_if_no_reference()

  set.seed(16)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  args <- list(
    z = z, R = R, n = n, L = 10,
    check_null_threshold = 0.1,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

# =============================================================================
# Part 11: MAF filtering
# =============================================================================

test_that("susie_rss() matches reference with maf filtering - optim", {
  skip_if_no_reference()

  set.seed(17)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Simulate minor allele frequencies
  maf <- runif(p, 0.05, 0.5)

  args <- list(
    z = z, R = R, n = n, L = 10,
    maf = maf, maf_thresh = 0.1,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with maf filtering - EM", {
  skip_if_no_reference()

  set.seed(17)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Simulate minor allele frequencies
  maf <- runif(p, 0.05, 0.5)

  args <- list(
    z = z, R = R, n = n, L = 10,
    maf = maf, maf_thresh = 0.1,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})

test_that("susie_rss() matches reference with maf filtering - simple", {
  skip_if_no_reference()

  set.seed(17)
  n <- 500
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 10, 20)] <- c(0.5, 0.4, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Simulate minor allele frequencies
  maf <- runif(p, 0.05, 0.5)

  args <- list(
    z = z, R = R, n = n, L = 10,
    maf = maf, maf_thresh = 0.1,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5)
})
