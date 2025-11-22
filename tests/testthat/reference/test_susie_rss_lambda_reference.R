# Source helper functions
source(file.path("..", "helper_reference.R"), local = TRUE)

context("susie_rss with lambda reference comparison")

# =============================================================================
# REFERENCE TESTS FOR susie_rss() with lambda > 0
# =============================================================================
#
# These tests compare susie_rss(lambda > 0) against susie_rss_lambda()
# from stephenslab/susieR@1f9166c
#

# =============================================================================
# Part 1: Different lambda values
# =============================================================================

test_that("susie_rss() matches reference with lambda=1e-5 - optim", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with lambda=1e-5 - EM", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with lambda=1e-5 - simple", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with lambda=0.1 - optim", {
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
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, L = 10, lambda = 0.1, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with lambda=0.1 - EM", {
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
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, L = 10, lambda = 0.1, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with lambda=0.1 - simple", {
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
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, L = 10, lambda = 0.1, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with lambda=0.5 - optim", {
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

  args <- list(z = z, R = R, L = 10, lambda = 0.5, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with lambda=0.5 - EM", {
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

  args <- list(z = z, R = R, L = 10, lambda = 0.5, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with lambda=0.5 - simple", {
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

  args <- list(z = z, R = R, L = 10, lambda = 0.5, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 2: Different L values
# =============================================================================

test_that("susie_rss() matches reference with different L values - optim", {
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

  # Test L=1
  args1 <- list(z = z, R = R, L = 1, lambda = 1e-5, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args1, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")

  # Test L=5
  args5 <- list(z = z, R = R, L = 5, lambda = 1e-5, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args5, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")

  # Test L=20
  args20 <- list(z = z, R = R, L = 20, lambda = 1e-5, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args20, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with different L values - EM", {
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

  # Test L=1
  args1 <- list(z = z, R = R, L = 1, lambda = 1e-5, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args1, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")

  # Test L=5
  args5 <- list(z = z, R = R, L = 5, lambda = 1e-5, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args5, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")

  # Test L=20
  args20 <- list(z = z, R = R, L = 20, lambda = 1e-5, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args20, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with different L values - simple", {
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

  # Test L=1
  args1 <- list(z = z, R = R, L = 1, lambda = 1e-5, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args1, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")

  # Test L=5
  args5 <- list(z = z, R = R, L = 5, lambda = 1e-5, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args5, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")

  # Test L=20
  args20 <- list(z = z, R = R, L = 20, lambda = 1e-5, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args20, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 3: estimate_prior_variance parameter
# =============================================================================

test_that("susie_rss() matches reference with estimate_prior_variance=FALSE - optim", {
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

  args <- list(
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = TRUE,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with estimate_prior_variance=FALSE - EM", {
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

  args <- list(
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = TRUE,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with estimate_prior_variance=FALSE - simple", {
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

  args <- list(
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = TRUE,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 4: estimate_residual_variance parameter
# =============================================================================

test_that("susie_rss() matches reference with estimate_residual_variance=FALSE - optim", {
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
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_residual_variance = FALSE,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with estimate_residual_variance=FALSE - EM", {
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
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_residual_variance = FALSE,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with estimate_residual_variance=FALSE - simple", {
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
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_residual_variance = FALSE,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with residual_variance fixed - optim", {
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
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with residual_variance fixed - EM", {
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
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with residual_variance fixed - simple", {
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
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 5: prior_variance parameter
# =============================================================================

test_that("susie_rss() matches reference with prior_variance=100 - optim", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, prior_variance = 100, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with prior_variance=100 - EM", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, prior_variance = 100, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with prior_variance=100 - simple", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, prior_variance = 100, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
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
  args <- list(z = z, R = R, L = 10, lambda = 1e-5, prior_weights = prior_weights, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
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
  args <- list(z = z, R = R, L = 10, lambda = 1e-5, prior_weights = prior_weights, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
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
  args <- list(z = z, R = R, L = 10, lambda = 1e-5, prior_weights = prior_weights, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 7: maf filtering
# =============================================================================

test_that("susie_rss() matches reference with maf filtering - optim", {
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

  # Simulate minor allele frequencies
  maf <- runif(p, 0.05, 0.5)

  args <- list(
    z = z, R = R, L = 10, lambda = 1e-5,
    maf = maf, maf_thresh = 0.1,
    estimate_prior_method = "optim",
    estimate_residual_variance = TRUE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with maf filtering - EM", {
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

  # Simulate minor allele frequencies
  maf <- runif(p, 0.05, 0.5)

  args <- list(
    z = z, R = R, L = 10, lambda = 1e-5,
    maf = maf, maf_thresh = 0.1,
    estimate_prior_method = "EM",
    estimate_residual_variance = TRUE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with maf filtering - simple", {
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

  # Simulate minor allele frequencies
  maf <- runif(p, 0.05, 0.5)

  args <- list(
    z = z, R = R, L = 10, lambda = 1e-5,
    maf = maf, maf_thresh = 0.1,
    estimate_prior_method = "simple",
    estimate_residual_variance = TRUE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 8: coverage and min_abs_corr
# =============================================================================

test_that("susie_rss() matches reference with coverage=0.99 - optim", {
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
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, coverage = 0.99, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with coverage=0.99 - EM", {
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
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, coverage = 0.99, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with coverage=0.99 - simple", {
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
  z <- with(ss, betahat / sebetahat)

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, coverage = 0.99, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with min_abs_corr=0.7 - optim", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, min_abs_corr = 0.7, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with min_abs_corr=0.7 - EM", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, min_abs_corr = 0.7, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with min_abs_corr=0.7 - simple", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, min_abs_corr = 0.7, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 9: prior_tol parameter
# =============================================================================

test_that("susie_rss() matches reference with prior_tol=1e-5 - optim", {
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

  # Disable check_null_threshold so we can see prior_tol effects
  args <- list(
    z = z, R = R, L = 10, lambda = 1e-5,
    prior_tol = 0.1,
    estimate_prior_method = "optim",
    estimate_residual_variance = TRUE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with prior_tol=1e-5 - EM", {
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

  args <- list(
    z = z, R = R, L = 10, lambda = 1e-5,
    prior_tol = 0.1,
    estimate_prior_method = "EM",
    estimate_residual_variance = TRUE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with prior_tol=1e-5 - simple", {
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

  args <- list(
    z = z, R = R, L = 10, lambda = 1e-5,
    prior_tol = 0.1,
    estimate_prior_method = "simple",
    estimate_residual_variance = TRUE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 10: check_null_threshold parameter
# =============================================================================

test_that("susie_rss() matches reference with check_null_threshold=0.1 - optim", {
  skip_if_no_reference()

  set.seed(14)
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
    z = z, R = R, L = 10, lambda = 1e-5,
    check_null_threshold = 0.1,
    estimate_prior_method = "optim",
    estimate_residual_variance = TRUE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with check_null_threshold=0.1 - EM", {
  skip_if_no_reference()

  set.seed(14)
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
    z = z, R = R, L = 10, lambda = 1e-5,
    check_null_threshold = 0.1,
    estimate_prior_method = "EM",
    estimate_residual_variance = TRUE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with check_null_threshold=0.1 - simple", {
  skip_if_no_reference()

  set.seed(14)
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
    z = z, R = R, L = 10, lambda = 1e-5,
    check_null_threshold = 0.1,
    estimate_prior_method = "simple",
    estimate_residual_variance = TRUE
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 11: intercept_value parameter
# =============================================================================

test_that("susie_rss() matches reference with intercept_value=0.5 - optim", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, intercept_value = 0.5, estimate_prior_method = "optim", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with intercept_value=0.5 - EM", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, intercept_value = 0.5, estimate_prior_method = "EM", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with intercept_value=0.5 - simple", {
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

  args <- list(z = z, R = R, L = 10, lambda = 1e-5, intercept_value = 0.5, estimate_prior_method = "simple", estimate_residual_variance = TRUE)
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

# =============================================================================
# Part 12: Combined parameter tests
# =============================================================================

test_that("susie_rss() matches reference with combined params - optim", {
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

  # Test combination: estimate_prior_variance=FALSE, estimate_residual_variance=FALSE
  args <- list(
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with combined params - EM", {
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
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})

test_that("susie_rss() matches reference with combined params - simple", {
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
    z = z, R = R, L = 10, lambda = 1e-5,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_rss", args, tolerance = 1e-5, ref_func_name = "susie_rss_lambda")
})
