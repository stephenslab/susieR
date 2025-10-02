# Source helper functions
source(file.path("..", "helpers", "helper_reference.R"), local = TRUE)

context("susie_rss with lambda reference comparison")

# =============================================================================
# REFERENCE TESTS FOR susie_rss() WITH LAMBDA > 0
# =============================================================================
#
# These tests compare the new susie_rss(lambda > 0) implementation against the
# reference package susie_rss_lambda() from stephenslab/susieR@1f9166c
#
# NOTE: Default differences between implementations:
# - Reference susie_rss_lambda(): estimate_residual_variance = TRUE (default)
# - New susie_rss(): estimate_residual_variance = FALSE (default)
# Tests explicitly set this parameter to ensure fair comparison.

# =============================================================================
# Part 1: Tests with estimate_residual_variance = TRUE
# =============================================================================

test_that("susie_rss(lambda=1e-5) matches reference with estimate_residual_variance=TRUE", {
  skip_if_no_reference()

  set.seed(11)
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

  # Both with estimate_residual_variance = TRUE
  args <- list(z = z, R = R, L = 10, lambda = 1e-5,
               estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda=0.1) matches reference with estimate_residual_variance=TRUE", {
  skip_if_no_reference()

  set.seed(12)
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

  args <- list(z = z, R = R, L = 10, lambda = 0.1,
               estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda=0.5) matches reference with estimate_residual_variance=TRUE", {
  skip_if_no_reference()

  set.seed(13)
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

  args <- list(z = z, R = R, L = 10, lambda = 0.5,
               estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda > 0) matches reference with different L values (estimate_residual_variance=TRUE)", {
  skip_if_no_reference()

  set.seed(14)
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
  args1 <- list(z = z, R = R, L = 1, lambda = 0.1,
                estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args1, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")

  # Test L=5
  args5 <- list(z = z, R = R, L = 5, lambda = 0.1,
                estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args5, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda > 0) matches reference with prior_variance (estimate_residual_variance=TRUE)", {
  skip_if_no_reference()

  set.seed(15)
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

  args <- list(z = z, R = R, L = 10, lambda = 0.1, prior_variance = 100,
               estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda > 0) matches reference with different estimate_prior_method (estimate_residual_variance=TRUE)", {
  skip_if_no_reference()

  set.seed(16)
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
  args_em <- list(z = z, R = R, L = 10, lambda = 0.1, estimate_prior_method = "EM",
                  estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args_em, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")

  # Test "simple"
  args_simple <- list(z = z, R = R, L = 10, lambda = 0.1, estimate_prior_method = "simple",
                      estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args_simple, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda > 0) matches reference with intercept_value (estimate_residual_variance=TRUE)", {
  skip_if_no_reference()

  set.seed(17)
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

  args <- list(z = z, R = R, L = 10, lambda = 0.1, intercept_value = 0.5,
               estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda > 0) matches reference with coverage (estimate_residual_variance=TRUE)", {
  skip_if_no_reference()

  set.seed(18)
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

  args <- list(z = z, R = R, L = 10, lambda = 0.1, coverage = 0.99,
               estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 2: Tests with estimate_residual_variance = FALSE
# =============================================================================

test_that("susie_rss(lambda=1e-5) matches reference with estimate_residual_variance=FALSE", {
  skip_if_no_reference()

  set.seed(21)
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

  # Both with estimate_residual_variance = FALSE
  args <- list(z = z, R = R, L = 10, lambda = 1e-5,
               estimate_residual_variance = FALSE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda=0.1) matches reference with estimate_residual_variance=FALSE", {
  skip_if_no_reference()

  set.seed(22)
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

  args <- list(z = z, R = R, L = 10, lambda = 0.1,
               estimate_residual_variance = FALSE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda=0.5) matches reference with estimate_residual_variance=FALSE", {
  skip_if_no_reference()

  set.seed(23)
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

  args <- list(z = z, R = R, L = 10, lambda = 0.5,
               estimate_residual_variance = FALSE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda > 0) matches reference with estimate_prior_variance=FALSE and estimate_residual_variance=FALSE", {
  skip_if_no_reference()

  set.seed(24)
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

  args <- list(z = z, R = R, L = 10, lambda = 0.2,
               estimate_prior_variance = FALSE,
               estimate_residual_variance = FALSE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda > 0) matches reference with different L values (estimate_residual_variance=FALSE)", {
  skip_if_no_reference()

  set.seed(25)
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
  args1 <- list(z = z, R = R, L = 1, lambda = 0.1,
                estimate_residual_variance = FALSE)
  compare_to_reference("susie_rss", args1, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")

  # Test L=5
  args5 <- list(z = z, R = R, L = 5, lambda = 0.1,
                estimate_residual_variance = FALSE)
  compare_to_reference("susie_rss", args5, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss(lambda > 0) matches reference with coverage (estimate_residual_variance=FALSE)", {
  skip_if_no_reference()

  set.seed(26)
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

  args <- list(z = z, R = R, L = 10, lambda = 0.1, coverage = 0.99,
               estimate_residual_variance = FALSE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5,
                       ref_func_name = "susie_rss_lambda")
})
