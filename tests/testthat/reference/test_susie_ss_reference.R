# Source helper functions
source(file.path("..", "helpers", "helper_reference.R"), local = TRUE)

context("susie_ss reference comparison")

# =============================================================================
# REFERENCE TESTS FOR susie_ss()
# =============================================================================
#
# These tests compare the new susie_ss() implementation against the reference
# susie_suff_stat() from stephenslab/susieR@1f9166c

test_that("susie_ss() matches reference with default parameters", {
  skip_if_no_reference()

  # Generate test data
  set.seed(1)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute sufficient statistics
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  # Test with defaults
  args <- list(
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    L = 10
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with X_colmeans and y_mean", {
  skip_if_no_reference()

  set.seed(2)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Compute sufficient statistics
  X_colmeans <- colMeans(X)
  y_mean <- mean(y)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - y_mean
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  args <- list(
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    L = 10,
    X_colmeans = X_colmeans,
    y_mean = y_mean
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with standardize=FALSE", {
  skip_if_no_reference()

  set.seed(3)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  args <- list(
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    L = 10,
    standardize = FALSE
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with estimate_prior_variance=FALSE", {
  skip_if_no_reference()

  set.seed(4)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  args <- list(
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    L = 10,
    estimate_prior_variance = FALSE
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with estimate_residual_variance=FALSE", {
  skip_if_no_reference()

  set.seed(5)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  args <- list(
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    L = 10,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with different L values", {
  skip_if_no_reference()

  set.seed(6)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  # Test L=1
  args1 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 1)
  compare_to_reference("susie_ss", args1, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test L=5
  args5 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 5)
  compare_to_reference("susie_ss", args5, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test L=20
  args20 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 20)
  compare_to_reference("susie_ss", args20, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with different estimate_prior_method", {
  skip_if_no_reference()

  set.seed(7)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  # Test "EM"
  args_em <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args_em, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test "simple"
  args_simple <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args_simple, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with prior_weights", {
  skip_if_no_reference()

  set.seed(8)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  # Uniform prior weights
  prior_weights <- rep(1/p, p)
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    prior_weights = prior_weights
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with different scaled_prior_variance", {
  skip_if_no_reference()

  set.seed(9)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  # Test with scaled_prior_variance = 0.5
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    scaled_prior_variance = 0.5
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with combined parameter variations", {
  skip_if_no_reference()

  set.seed(10)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  # Test combination: standardize=FALSE, estimate_prior_variance=FALSE
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    standardize = FALSE,
    estimate_prior_variance = FALSE
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with different coverage values", {
  skip_if_no_reference()

  set.seed(11)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  # Test coverage = 0.99
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    coverage = 0.99
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with different min_abs_corr values", {
  skip_if_no_reference()

  set.seed(12)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - mean(y)
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  # Test min_abs_corr = 0.7
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    min_abs_corr = 0.7
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})
