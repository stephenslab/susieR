# Source helper functions
source(file.path("..", "helpers", "helper_reference.R"), local = TRUE)

context("susie_ss reference comparison")

# =============================================================================
# REFERENCE TESTS FOR susie_ss()
# =============================================================================
#
# These tests compare the new susie_ss() implementation against the reference
# susie_suff_stat() from stephenslab/susieR@1f9166c
#
# Tests cover all major parameters and their combinations with all three
# prior variance optimization methods: "optim", "EM", "simple"

# =============================================================================
# Part 1: Basic Parameter Tests
# =============================================================================

test_that("susie_ss() matches reference with default parameters - optim", {
  skip_if_no_reference()

  set.seed(1)
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with default parameters - EM", {
  skip_if_no_reference()

  set.seed(1)
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with default parameters - simple", {
  skip_if_no_reference()

  set.seed(1)
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 2: X_colmeans and y_mean (intercept estimation)
# =============================================================================

test_that("susie_ss() matches reference with X_colmeans and y_mean - optim", {
  skip_if_no_reference()

  set.seed(2)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_colmeans <- colMeans(X)
  y_mean <- mean(y)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - y_mean
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    X_colmeans = X_colmeans, y_mean = y_mean,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with X_colmeans and y_mean - EM", {
  skip_if_no_reference()

  set.seed(2)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_colmeans <- colMeans(X)
  y_mean <- mean(y)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - y_mean
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    X_colmeans = X_colmeans, y_mean = y_mean,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with X_colmeans and y_mean - simple", {
  skip_if_no_reference()

  set.seed(2)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  X_colmeans <- colMeans(X)
  y_mean <- mean(y)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  y_centered <- y - y_mean
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)
  yty <- sum(y_centered^2)

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    X_colmeans = X_colmeans, y_mean = y_mean,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 3: standardize parameter
# =============================================================================

test_that("susie_ss() matches reference with standardize=FALSE - optim", {
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    standardize = FALSE, estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with standardize=FALSE - EM", {
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    standardize = FALSE, estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with standardize=FALSE - simple", {
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    standardize = FALSE, estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 4: estimate_prior_variance=FALSE
# =============================================================================

test_that("susie_ss() matches reference with estimate_prior_variance=FALSE - optim", {
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_variance = FALSE, estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with estimate_prior_variance=FALSE - EM", {
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_variance = FALSE, estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with estimate_prior_variance=FALSE - simple", {
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_variance = FALSE, estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 5: estimate_residual_variance parameter
# =============================================================================

test_that("susie_ss() matches reference with estimate_residual_variance=FALSE - optim", {
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_residual_variance = FALSE, residual_variance = 1.0,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with estimate_residual_variance=FALSE - EM", {
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_residual_variance = FALSE, residual_variance = 1.0,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with estimate_residual_variance=FALSE - simple", {
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_residual_variance = FALSE, residual_variance = 1.0,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 6: Different L values
# =============================================================================

test_that("susie_ss() matches reference with different L values - optim", {
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
  args1 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 1, estimate_prior_method = "optim")
  compare_to_reference("susie_ss", args1, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test L=5
  args5 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 5, estimate_prior_method = "optim")
  compare_to_reference("susie_ss", args5, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test L=20
  args20 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 20, estimate_prior_method = "optim")
  compare_to_reference("susie_ss", args20, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with different L values - EM", {
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
  args1 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 1, estimate_prior_method = "EM")
  compare_to_reference("susie_ss", args1, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test L=5
  args5 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 5, estimate_prior_method = "EM")
  compare_to_reference("susie_ss", args5, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test L=20
  args20 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 20, estimate_prior_method = "EM")
  compare_to_reference("susie_ss", args20, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with different L values - simple", {
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
  args1 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 1, estimate_prior_method = "simple")
  compare_to_reference("susie_ss", args1, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test L=5
  args5 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 5, estimate_prior_method = "simple")
  compare_to_reference("susie_ss", args5, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test L=20
  args20 <- list(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 20, estimate_prior_method = "simple")
  compare_to_reference("susie_ss", args20, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 7: prior_weights
# =============================================================================

test_that("susie_ss() matches reference with prior_weights - optim", {
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

  # Use non-uniform prior weights
  prior_weights <- runif(p)
  prior_weights <- prior_weights / sum(prior_weights)
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    prior_weights = prior_weights, estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with prior_weights - EM", {
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

  # Use non-uniform prior weights
  prior_weights <- runif(p)
  prior_weights <- prior_weights / sum(prior_weights)
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    prior_weights = prior_weights, estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with prior_weights - simple", {
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

  # Use non-uniform prior weights
  prior_weights <- runif(p)
  prior_weights <- prior_weights / sum(prior_weights)
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    prior_weights = prior_weights, estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 8: scaled_prior_variance
# =============================================================================

test_that("susie_ss() matches reference with scaled_prior_variance - optim", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    scaled_prior_variance = 0.5, estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with scaled_prior_variance - EM", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    scaled_prior_variance = 0.5, estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with scaled_prior_variance - simple", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    scaled_prior_variance = 0.5, estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 9: coverage, min_abs_corr, and n_purity
# =============================================================================

test_that("susie_ss() matches reference with coverage=0.99 - optim", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    coverage = 0.99, estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with coverage=0.99 - EM", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    coverage = 0.99, estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with coverage=0.99 - simple", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    coverage = 0.99, estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with min_abs_corr=0.7 - optim", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    min_abs_corr = 0.7, estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with min_abs_corr=0.7 - EM", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    min_abs_corr = 0.7, estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with min_abs_corr=0.7 - simple", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    min_abs_corr = 0.7, estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with n_purity=3 - optim", {
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

  # Set seed before calling to ensure same variant selection
  set.seed(999)
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    n_purity = 3, estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with n_purity=3 - EM", {
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

  # Set seed before calling to ensure same variant selection
  set.seed(999)
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    n_purity = 3, estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with n_purity=3 - simple", {
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

  # Set seed before calling to ensure same variant selection
  set.seed(999)
  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    n_purity = 3, estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 10: Combined parameter variations
# =============================================================================

test_that("susie_ss() matches reference with combined parameters - optim", {
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

  # Test combination: standardize=FALSE, estimate_prior_variance=FALSE
  args1 <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    standardize = FALSE, estimate_prior_variance = FALSE,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args1, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test combination: estimate_prior_variance=FALSE, estimate_residual_variance=FALSE
  args2 <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_variance = FALSE, estimate_residual_variance = FALSE,
    residual_variance = 1.0, estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args2, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with combined parameters - EM", {
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

  # Test combination: standardize=FALSE, estimate_prior_variance=FALSE
  args1 <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    standardize = FALSE, estimate_prior_variance = FALSE,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args1, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test combination: estimate_prior_variance=FALSE, estimate_residual_variance=FALSE
  args2 <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_variance = FALSE, estimate_residual_variance = FALSE,
    residual_variance = 1.0, estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args2, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with combined parameters - simple", {
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

  # Test combination: standardize=FALSE, estimate_prior_variance=FALSE
  args1 <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    standardize = FALSE, estimate_prior_variance = FALSE,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args1, tolerance = 1e-5, ref_func_name = "susie_suff_stat")

  # Test combination: estimate_prior_variance=FALSE, estimate_residual_variance=FALSE
  args2 <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    estimate_prior_variance = FALSE, estimate_residual_variance = FALSE,
    residual_variance = 1.0, estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args2, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 11: prior_tol parameter
# =============================================================================

test_that("susie_ss() matches reference with prior_tol=1e-5 - optim", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    prior_tol = 0.1,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with prior_tol=0.1 - EM", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    prior_tol = 0.1,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with prior_tol=0.1 - simple", {
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

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    prior_tol = 0.1,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 12: check_null_threshold parameter
# =============================================================================

test_that("susie_ss() matches reference with check_null_threshold=0.1 - optim", {
  skip_if_no_reference()

  set.seed(13)
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    check_null_threshold = 0.1, estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with check_null_threshold=0.1 - EM", {
  skip_if_no_reference()

  set.seed(13)
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    check_null_threshold = 0.1, estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with check_null_threshold=0.1 - simple", {
  skip_if_no_reference()

  set.seed(13)
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
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    check_null_threshold = 0.1, estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

# =============================================================================
# Part 13: maf and maf_thresh parameters
# =============================================================================

test_that("susie_ss() matches reference with maf filtering - optim", {
  skip_if_no_reference()

  set.seed(14)
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

  # Simulate minor allele frequencies
  maf <- runif(p, 0.05, 0.5)

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    maf = maf, maf_thresh = 0.1,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with maf filtering - EM", {
  skip_if_no_reference()

  set.seed(14)
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

  # Simulate minor allele frequencies
  maf <- runif(p, 0.05, 0.5)

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    maf = maf, maf_thresh = 0.1,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})

test_that("susie_ss() matches reference with maf filtering - simple", {
  skip_if_no_reference()

  set.seed(14)
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

  # Simulate minor allele frequencies
  maf <- runif(p, 0.05, 0.5)

  args <- list(
    XtX = XtX, Xty = Xty, yty = yty, n = n, L = 10,
    maf = maf, maf_thresh = 0.1,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie_ss", args, tolerance = 1e-5, ref_func_name = "susie_suff_stat")
})
