# Source helper functions
source(file.path("..", "helper_reference.R"), local = TRUE)

context("susie reference comparison")

# =============================================================================
# REFERENCE TESTS FOR susie()
# =============================================================================
#
# These functions compare the new susieR implementation against the reference
# package (stephenslab/susieR@1f9166c) to ensure results are identical.
#
# Tests cover all major parameters and their combinations with all three
# prior variance optimization methods: "optim", "EM", "simple"

# =============================================================================
# Part 1: Basic Parameter Tests
# =============================================================================

test_that("susie() matches reference with default parameters - optim", {
  skip_if_no_reference()

  set.seed(1)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with default parameters - EM", {
  skip_if_no_reference()

  set.seed(1)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, estimate_prior_method = "EM")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with default parameters - simple", {
  skip_if_no_reference()

  set.seed(1)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, estimate_prior_method = "simple")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 2: standardize parameter
# =============================================================================

test_that("susie() matches reference with standardize=FALSE - optim", {
  skip_if_no_reference()

  set.seed(2)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, standardize = FALSE, estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with standardize=FALSE - EM", {
  skip_if_no_reference()

  set.seed(2)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, standardize = FALSE, estimate_prior_method = "EM")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with standardize=FALSE - simple", {
  skip_if_no_reference()

  set.seed(2)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, standardize = FALSE, estimate_prior_method = "simple")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 3: intercept parameter
# =============================================================================

test_that("susie() matches reference with intercept=FALSE - optim", {
  skip_if_no_reference()

  set.seed(3)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, intercept = FALSE, estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with intercept=FALSE - EM", {
  skip_if_no_reference()

  set.seed(3)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, intercept = FALSE, estimate_prior_method = "EM")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with intercept=FALSE - simple", {
  skip_if_no_reference()

  set.seed(3)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, intercept = FALSE, estimate_prior_method = "simple")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 4: estimate_prior_variance=FALSE
# =============================================================================

test_that("susie() matches reference with estimate_prior_variance=FALSE", {
  skip_if_no_reference()

  set.seed(4)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # When estimate_prior_variance=FALSE, the method doesn't matter
  args <- list(X = X, y = y, L = 10, estimate_prior_variance = FALSE)
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 5: estimate_residual_variance parameter
# =============================================================================

test_that("susie() matches reference with estimate_residual_variance=FALSE - optim", {
  skip_if_no_reference()

  set.seed(5)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with estimate_residual_variance=FALSE - EM", {
  skip_if_no_reference()

  set.seed(5)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with estimate_residual_variance=FALSE - simple", {
  skip_if_no_reference()

  set.seed(5)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 6: Sparse matrix input
# =============================================================================

test_that("susie() matches reference with sparse matrix input - optim", {
  skip_if_no_reference()

  set.seed(6)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  X_sparse <- Matrix::Matrix(X, sparse = TRUE)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X_sparse, y = y, L = 10, estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with sparse matrix input - EM", {
  skip_if_no_reference()

  set.seed(6)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  X_sparse <- Matrix::Matrix(X, sparse = TRUE)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X_sparse, y = y, L = 10, estimate_prior_method = "EM")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with sparse matrix input - simple", {
  skip_if_no_reference()

  set.seed(6)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  X_sparse <- Matrix::Matrix(X, sparse = TRUE)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X_sparse, y = y, L = 10, estimate_prior_method = "simple")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 7: Different L values
# =============================================================================

test_that("susie() matches reference with different L values - optim", {
  skip_if_no_reference()

  set.seed(7)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Test L=1
  args1 <- list(X = X, y = y, L = 1, estimate_prior_method = "optim")
  compare_to_reference("susie", args1, tolerance = 1e-5)

  # Test L=5
  args5 <- list(X = X, y = y, L = 5, estimate_prior_method = "optim")
  compare_to_reference("susie", args5, tolerance = 1e-5)

  # Test L=20
  args20 <- list(X = X, y = y, L = 20, estimate_prior_method = "optim")
  compare_to_reference("susie", args20, tolerance = 1e-5)
})

test_that("susie() matches reference with different L values - EM", {
  skip_if_no_reference()

  set.seed(7)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Test L=1
  args1 <- list(X = X, y = y, L = 1, estimate_prior_method = "EM")
  compare_to_reference("susie", args1, tolerance = 1e-5)

  # Test L=5
  args5 <- list(X = X, y = y, L = 5, estimate_prior_method = "EM")
  compare_to_reference("susie", args5, tolerance = 1e-5)

  # Test L=20
  args20 <- list(X = X, y = y, L = 20, estimate_prior_method = "EM")
  compare_to_reference("susie", args20, tolerance = 1e-5)
})

test_that("susie() matches reference with different L values - simple", {
  skip_if_no_reference()

  set.seed(7)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Test L=1
  args1 <- list(X = X, y = y, L = 1, estimate_prior_method = "simple")
  compare_to_reference("susie", args1, tolerance = 1e-5)

  # Test L=5
  args5 <- list(X = X, y = y, L = 5, estimate_prior_method = "simple")
  compare_to_reference("susie", args5, tolerance = 1e-5)

  # Test L=20
  args20 <- list(X = X, y = y, L = 20, estimate_prior_method = "simple")
  compare_to_reference("susie", args20, tolerance = 1e-5)
})

# =============================================================================
# Part 8: prior_weights
# =============================================================================

test_that("susie() matches reference with prior_weights - optim", {
  skip_if_no_reference()

  set.seed(8)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Use non-uniform prior weights
  prior_weights <- runif(p)
  prior_weights <- prior_weights / sum(prior_weights)
  args <- list(X = X, y = y, L = 10, prior_weights = prior_weights, estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with prior_weights - EM", {
  skip_if_no_reference()

  set.seed(8)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Use non-uniform prior weights
  prior_weights <- runif(p)
  prior_weights <- prior_weights / sum(prior_weights)
  args <- list(X = X, y = y, L = 10, prior_weights = prior_weights, estimate_prior_method = "EM")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with prior_weights - simple", {
  skip_if_no_reference()

  set.seed(8)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Use non-uniform prior weights
  prior_weights <- runif(p)
  prior_weights <- prior_weights / sum(prior_weights)
  args <- list(X = X, y = y, L = 10, prior_weights = prior_weights, estimate_prior_method = "simple")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 9: scaled_prior_variance
# =============================================================================

test_that("susie() matches reference with scaled_prior_variance - optim", {
  skip_if_no_reference()

  set.seed(9)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, scaled_prior_variance = 0.5, estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with scaled_prior_variance - EM", {
  skip_if_no_reference()

  set.seed(9)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, scaled_prior_variance = 0.5, estimate_prior_method = "EM")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with scaled_prior_variance - simple", {
  skip_if_no_reference()

  set.seed(9)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, scaled_prior_variance = 0.5, estimate_prior_method = "simple")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 10: coverage and min_abs_corr
# =============================================================================

test_that("susie() matches reference with coverage=0.99 - optim", {
  skip_if_no_reference()

  set.seed(10)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, coverage = 0.99, estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with coverage=0.99 - EM", {
  skip_if_no_reference()

  set.seed(10)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, coverage = 0.99, estimate_prior_method = "EM")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with coverage=0.99 - simple", {
  skip_if_no_reference()

  set.seed(10)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, coverage = 0.99, estimate_prior_method = "simple")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with min_abs_corr=0.7 - optim", {
  skip_if_no_reference()

  set.seed(11)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, min_abs_corr = 0.7, estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with min_abs_corr=0.7 - EM", {
  skip_if_no_reference()

  set.seed(11)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, min_abs_corr = 0.7, estimate_prior_method = "EM")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with min_abs_corr=0.7 - simple", {
  skip_if_no_reference()

  set.seed(11)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L = 10, min_abs_corr = 0.7, estimate_prior_method = "simple")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 11: Combined parameter variations
# =============================================================================

test_that("susie() matches reference with combined parameters - optim", {
  skip_if_no_reference()

  set.seed(13)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Test combination: standardize=FALSE, intercept=FALSE
  args1 <- list(X = X, y = y, L = 10, standardize = FALSE, intercept = FALSE, estimate_prior_method = "optim")
  compare_to_reference("susie", args1, tolerance = 1e-5)

  # Test combination: estimate_prior_variance=FALSE, estimate_residual_variance=FALSE
  args2 <- list(
    X = X, y = y, L = 10,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie", args2, tolerance = 1e-5)
})

test_that("susie() matches reference with combined parameters - EM", {
  skip_if_no_reference()

  set.seed(13)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Test combination: standardize=FALSE, intercept=FALSE
  args1 <- list(X = X, y = y, L = 10, standardize = FALSE, intercept = FALSE, estimate_prior_method = "EM")
  compare_to_reference("susie", args1, tolerance = 1e-5)

  # Test combination: estimate_prior_variance=FALSE, estimate_residual_variance=FALSE
  args2 <- list(
    X = X, y = y, L = 10,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie", args2, tolerance = 1e-5)
})

test_that("susie() matches reference with combined parameters - simple", {
  skip_if_no_reference()

  set.seed(13)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Test combination: standardize=FALSE, intercept=FALSE
  args1 <- list(X = X, y = y, L = 10, standardize = FALSE, intercept = FALSE, estimate_prior_method = "simple")
  compare_to_reference("susie", args1, tolerance = 1e-5)

  # Test combination: estimate_prior_variance=FALSE, estimate_residual_variance=FALSE
  args2 <- list(
    X = X, y = y, L = 10,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = FALSE,
    residual_variance = 1.0,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie", args2, tolerance = 1e-5)
})

# =============================================================================
# Part 12: prior_tol parameter
# =============================================================================

test_that("susie() matches reference with prior_tol=0.1 - optim", {
  skip_if_no_reference()

  set.seed(14)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    prior_tol = 0.1,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with prior_tol=0.1 - EM", {
  skip_if_no_reference()

  set.seed(14)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    prior_tol = 0.1,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with prior_tol=0.1 - simple", {
  skip_if_no_reference()

  set.seed(14)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    prior_tol = 0.1,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 13: check_null_threshold parameter
# =============================================================================

test_that("susie() matches reference with check_null_threshold=0.1 - optim", {
  skip_if_no_reference()

  set.seed(15)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    check_null_threshold = 0.1,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with check_null_threshold=0.1 - EM", {
  skip_if_no_reference()

  set.seed(15)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    check_null_threshold = 0.1,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with check_null_threshold=0.1 - simple", {
  skip_if_no_reference()

  set.seed(15)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    check_null_threshold = 0.1,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with check_null_threshold=0.5 - optim", {
  skip_if_no_reference()

  set.seed(16)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    check_null_threshold = 0.5,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with check_null_threshold=0.5 - EM", {
  skip_if_no_reference()

  set.seed(16)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    check_null_threshold = 0.5,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with check_null_threshold=0.5 - simple", {
  skip_if_no_reference()

  set.seed(16)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(
    X = X, y = y, L = 10,
    check_null_threshold = 0.5,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 14: residual_variance bounds
# =============================================================================

test_that("susie() matches reference with residual_variance_upperbound - optim", {
  skip_if_no_reference()

  set.seed(17)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Set upperbound lower than natural variance (~1.07) to ensure it's binding
  args <- list(
    X = X, y = y, L = 10,
    residual_variance_upperbound = 0.8,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with residual_variance_upperbound - EM", {
  skip_if_no_reference()

  set.seed(17)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Set upperbound lower than natural variance (~1.07) to ensure it's binding
  args <- list(
    X = X, y = y, L = 10,
    residual_variance_upperbound = 0.8,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with residual_variance_upperbound - simple", {
  skip_if_no_reference()

  set.seed(17)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Set upperbound lower than natural variance (~1.07) to ensure it's binding
  args <- list(
    X = X, y = y, L = 10,
    residual_variance_upperbound = 0.8,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with residual_variance_lowerbound - optim", {
  skip_if_no_reference()

  set.seed(18)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Set lowerbound higher than natural variance (~1.07) to ensure it's binding
  args <- list(
    X = X, y = y, L = 10,
    residual_variance_lowerbound = 1.5,
    estimate_prior_method = "optim"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with residual_variance_lowerbound - EM", {
  skip_if_no_reference()

  set.seed(18)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Set lowerbound higher than natural variance (~1.07) to ensure it's binding
  args <- list(
    X = X, y = y, L = 10,
    residual_variance_lowerbound = 1.5,
    estimate_prior_method = "EM"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with residual_variance_lowerbound - simple", {
  skip_if_no_reference()

  set.seed(18)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Set lowerbound higher than natural variance (~1.07) to ensure it's binding
  args <- list(
    X = X, y = y, L = 10,
    residual_variance_lowerbound = 1.5,
    estimate_prior_method = "simple"
  )
  compare_to_reference("susie", args, tolerance = 1e-5)
})

# =============================================================================
# Part 15: na.rm parameter
# =============================================================================

test_that("susie() matches reference with na.rm=TRUE - optim", {
  skip_if_no_reference()

  set.seed(19)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Introduce NA values
  y[c(1, 25, 50)] <- NA

  args <- list(X = X, y = y, L = 10, na.rm = TRUE, estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with na.rm=TRUE - EM", {
  skip_if_no_reference()

  set.seed(19)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Introduce NA values
  y[c(1, 25, 50)] <- NA

  args <- list(X = X, y = y, L = 10, na.rm = TRUE, estimate_prior_method = "EM")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with na.rm=TRUE - simple", {
  skip_if_no_reference()

  set.seed(19)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Introduce NA values
  y[c(1, 25, 50)] <- NA

  args <- list(X = X, y = y, L = 10, na.rm = TRUE, estimate_prior_method = "simple")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with na.rm=TRUE and single NA - optim", {
  skip_if_no_reference()

  set.seed(20)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Single NA (the bug report case)
  y[1] <- NA

  args <- list(X = X, y = y, L = 10, na.rm = TRUE, estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with na.rm=TRUE and standardize=FALSE", {
  skip_if_no_reference()

  set.seed(21)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Introduce NA values
  y[c(5, 10, 15)] <- NA

  args <- list(X = X, y = y, L = 10, na.rm = TRUE, standardize = FALSE,
               estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

test_that("susie() matches reference with na.rm=TRUE and intercept=FALSE", {
  skip_if_no_reference()

  set.seed(22)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Introduce NA values
  y[c(10, 20, 30)] <- NA

  args <- list(X = X, y = y, L = 10, na.rm = TRUE, intercept = FALSE,
               estimate_prior_method = "optim")
  compare_to_reference("susie", args, tolerance = 1e-5)
})

