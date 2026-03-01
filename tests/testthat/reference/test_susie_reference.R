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

# =============================================================================
# Part 16: model_init parameter (dev) vs s_init (reference)
# =============================================================================
#
# These tests verify that our model_init parameter produces identical results
# to the reference package's s_init parameter. Each test runs an initial susie
# fit on both packages, then passes the result as model_init/s_init to a
# second call and compares outputs.

test_that("susie() matches reference with model_init - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(23)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Run initial fit on both packages (short run to get a non-trivial init)
  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  # Run with model_init (dev) / s_init (ref)
  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   estimate_prior_method = "optim")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init - EM", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(23)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    estimate_prior_method = "EM")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   estimate_prior_method = "EM")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   estimate_prior_method = "EM")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(23)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   estimate_prior_method = "simple")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# model_init with estimate_residual_variance=FALSE
test_that("susie() matches reference with model_init and estimate_residual_variance=FALSE - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(24)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   estimate_residual_variance = FALSE, residual_variance = 1.0,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   estimate_residual_variance = FALSE, residual_variance = 1.0,
                   estimate_prior_method = "optim")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init and estimate_residual_variance=FALSE - EM", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(24)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    estimate_prior_method = "EM")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   estimate_residual_variance = FALSE, residual_variance = 1.0,
                   estimate_prior_method = "EM")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   estimate_residual_variance = FALSE, residual_variance = 1.0,
                   estimate_prior_method = "EM")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init and estimate_residual_variance=FALSE - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(24)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   estimate_residual_variance = FALSE, residual_variance = 1.0,
                   estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   estimate_residual_variance = FALSE, residual_variance = 1.0,
                   estimate_prior_method = "simple")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# model_init with estimate_prior_variance=FALSE
test_that("susie() matches reference with model_init and estimate_prior_variance=FALSE", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(25)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    estimate_prior_variance = FALSE)
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   estimate_prior_variance = FALSE)
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   estimate_prior_variance = FALSE)

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# model_init with standardize=FALSE
test_that("susie() matches reference with model_init and standardize=FALSE - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(26)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    standardize = FALSE, estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   standardize = FALSE, estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   standardize = FALSE, estimate_prior_method = "optim")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init and standardize=FALSE - EM", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(26)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    standardize = FALSE, estimate_prior_method = "EM")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   standardize = FALSE, estimate_prior_method = "EM")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   standardize = FALSE, estimate_prior_method = "EM")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init and standardize=FALSE - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(26)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    standardize = FALSE, estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   standardize = FALSE, estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   standardize = FALSE, estimate_prior_method = "simple")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# model_init with intercept=FALSE
test_that("susie() matches reference with model_init and intercept=FALSE - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(27)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    intercept = FALSE, estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   intercept = FALSE, estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   intercept = FALSE, estimate_prior_method = "optim")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init and intercept=FALSE - EM", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(27)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    intercept = FALSE, estimate_prior_method = "EM")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   intercept = FALSE, estimate_prior_method = "EM")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   intercept = FALSE, estimate_prior_method = "EM")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init and intercept=FALSE - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(27)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    intercept = FALSE, estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   intercept = FALSE, estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   intercept = FALSE, estimate_prior_method = "simple")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# model_init with L expansion (second call requests more effects than init)
test_that("susie() matches reference with model_init and L expansion - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(28)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Initial fit with L=3
  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  # Second fit with L=10 (expansion)
  dev_args <- list(X = X, y = y, L = 10, model_init = dev_init,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 10, s_init = ref_init,
                   estimate_prior_method = "optim")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init and L expansion - EM", {
  skip("Intentional improvement: susieR2.0 preserves fitted V during L expansion; reference resets all V to default. EM updates V incrementally so intermediate ELBOs differ, but final posteriors converge.")
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(28)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    estimate_prior_method = "EM")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 10, model_init = dev_init,
                   estimate_prior_method = "EM")
  ref_args <- list(X = X, y = y, L = 10, s_init = ref_init,
                   estimate_prior_method = "EM")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init and L expansion - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(28)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 10, model_init = dev_init,
                   estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 10, s_init = ref_init,
                   estimate_prior_method = "simple")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# model_init with combined standardize=FALSE and intercept=FALSE
test_that("susie() matches reference with model_init, standardize=FALSE, intercept=FALSE - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(29)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    standardize = FALSE, intercept = FALSE,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   standardize = FALSE, intercept = FALSE,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   standardize = FALSE, intercept = FALSE,
                   estimate_prior_method = "optim")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# =============================================================================
# Part 17: model_init with L expansion - deeper probing for differences
# =============================================================================
#
# These tests specifically target L expansion (model_init has fewer effects
# than the requested L) with various parameter combinations to find behavioral
# differences between model_init (dev) and s_init (ref).

# L expansion with estimate_prior_variance=FALSE
# When V is never re-estimated, any difference in V initialization should
# propagate through all iterations and affect final posteriors.
test_that("susie() matches reference with model_init, L expansion, estimate_prior_variance=FALSE - optim", {
  skip("Intentional improvement: susieR2.0 preserves fitted V during L expansion; reference resets all V to default. With estimate_prior_variance=FALSE, the V difference persists in final posteriors.")
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(30)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Initial fit with L=3
  init_args <- list(X = X, y = y, L = 3, max_iter = 5,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  # Second fit with L=8, estimate_prior_variance=FALSE
  dev_args <- list(X = X, y = y, L = 8, model_init = dev_init,
                   estimate_prior_variance = FALSE,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 8, s_init = ref_init,
                   estimate_prior_variance = FALSE,
                   estimate_prior_method = "optim")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init, L expansion, estimate_prior_variance=FALSE - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(30)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 5,
                    estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 8, model_init = dev_init,
                   estimate_prior_variance = FALSE,
                   estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 8, s_init = ref_init,
                   estimate_prior_variance = FALSE,
                   estimate_prior_method = "simple")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# L expansion with BOTH estimate_prior_variance=FALSE AND estimate_residual_variance=FALSE
# Fully constrained variances - any V initialization difference is permanent.
test_that("susie() matches reference with model_init, L expansion, both variances fixed - optim", {
  skip("Intentional improvement: susieR2.0 preserves fitted V during L expansion; reference resets all V to default. With both variances fixed, the V difference persists in final posteriors.")
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(31)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 5,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 8, model_init = dev_init,
                   estimate_prior_variance = FALSE,
                   estimate_residual_variance = FALSE,
                   residual_variance = 1.0,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 8, s_init = ref_init,
                   estimate_prior_variance = FALSE,
                   estimate_residual_variance = FALSE,
                   residual_variance = 1.0,
                   estimate_prior_method = "optim")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init, L expansion, both variances fixed - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(31)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 5,
                    estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 8, model_init = dev_init,
                   estimate_prior_variance = FALSE,
                   estimate_residual_variance = FALSE,
                   residual_variance = 1.0,
                   estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 8, s_init = ref_init,
                   estimate_prior_variance = FALSE,
                   estimate_residual_variance = FALSE,
                   residual_variance = 1.0,
                   estimate_prior_method = "simple")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# L expansion with standardize=FALSE
test_that("susie() matches reference with model_init, L expansion, standardize=FALSE - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(32)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    standardize = FALSE, estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 10, model_init = dev_init,
                   standardize = FALSE, estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 10, s_init = ref_init,
                   standardize = FALSE, estimate_prior_method = "optim")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init, L expansion, standardize=FALSE - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(32)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    standardize = FALSE, estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 10, model_init = dev_init,
                   standardize = FALSE, estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 10, s_init = ref_init,
                   standardize = FALSE, estimate_prior_method = "simple")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# L expansion with intercept=FALSE
test_that("susie() matches reference with model_init, L expansion, intercept=FALSE - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(33)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    intercept = FALSE, estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 10, model_init = dev_init,
                   intercept = FALSE, estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 10, s_init = ref_init,
                   intercept = FALSE, estimate_prior_method = "optim")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init, L expansion, intercept=FALSE - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(33)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    intercept = FALSE, estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 10, model_init = dev_init,
                   intercept = FALSE, estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 10, s_init = ref_init,
                   intercept = FALSE, estimate_prior_method = "simple")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# L expansion with non-default scaled_prior_variance
test_that("susie() matches reference with model_init, L expansion, scaled_prior_variance=0.5 - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(34)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    scaled_prior_variance = 0.5,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 10, model_init = dev_init,
                   scaled_prior_variance = 0.5,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 10, s_init = ref_init,
                   scaled_prior_variance = 0.5,
                   estimate_prior_method = "optim")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# =============================================================================
# Part 18: model_init with L contraction (model_init has more effects than L)
# =============================================================================

test_that("susie() matches reference with model_init, L contraction - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(35)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Initial fit with L=10
  init_args <- list(X = X, y = y, L = 10, max_iter = 5,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  # Second fit with L=5 (contraction - model_init has more effects)
  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   estimate_prior_method = "optim")

  dev_result <- suppressWarnings(suppressMessages(do.call(dev_env$env[["susie"]], dev_args)))
  ref_result <- suppressWarnings(suppressMessages(do.call(ref_env$env[["susie"]], ref_args)))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init, L contraction - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(35)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 10, max_iter = 5,
                    estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   estimate_prior_method = "simple")

  dev_result <- suppressWarnings(suppressMessages(do.call(dev_env$env[["susie"]], dev_args)))
  ref_result <- suppressWarnings(suppressMessages(do.call(ref_env$env[["susie"]], ref_args)))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# =============================================================================
# Part 19: model_init with susie_init_coef
# =============================================================================

test_that("susie() matches reference with susie_init_coef as model_init - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(36)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Create init from known coefficients
  dev_coef_init <- dev_env$env[["susie_init_coef"]](1:4, c(2, 3, -2, 1.5), p)
  ref_coef_init <- ref_env$env[["susie_init_coef"]](1:4, c(2, 3, -2, 1.5), p)

  dev_args <- list(X = X, y = y, L = 10, model_init = dev_coef_init,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 10, s_init = ref_coef_init,
                   estimate_prior_method = "optim")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with susie_init_coef as model_init - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(36)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_coef_init <- dev_env$env[["susie_init_coef"]](1:4, c(2, 3, -2, 1.5), p)
  ref_coef_init <- ref_env$env[["susie_init_coef"]](1:4, c(2, 3, -2, 1.5), p)

  dev_args <- list(X = X, y = y, L = 10, model_init = dev_coef_init,
                   estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 10, s_init = ref_coef_init,
                   estimate_prior_method = "simple")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# susie_init_coef with same L (no expansion)
test_that("susie() matches reference with susie_init_coef as model_init, matching L - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(37)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_coef_init <- dev_env$env[["susie_init_coef"]](1:4, c(2, 3, -2, 1.5), p)
  ref_coef_init <- ref_env$env[["susie_init_coef"]](1:4, c(2, 3, -2, 1.5), p)

  dev_args <- list(X = X, y = y, L = 4, model_init = dev_coef_init,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 4, s_init = ref_coef_init,
                   estimate_prior_method = "optim")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# =============================================================================
# Part 20: model_init with null_weight
# =============================================================================

test_that("susie() matches reference with model_init and null_weight - optim", {
  skip("Bug fix: susieR1.0 incorrectly set lpo=0 for infinite shat2 (null column), ignoring the actual null_weight. susieR2.0 correctly uses log(prior_weights).")
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(38)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    null_weight = 0.5,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   null_weight = 0.5,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   null_weight = 0.5,
                   estimate_prior_method = "optim")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init and null_weight - simple", {
  skip("Bug fix: susieR1.0 incorrectly set lpo=0 for infinite shat2 (null column), ignoring the actual null_weight. susieR2.0 correctly uses log(prior_weights).")
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(38)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    null_weight = 0.5,
                    estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   null_weight = 0.5,
                   estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   null_weight = 0.5,
                   estimate_prior_method = "simple")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# =============================================================================
# Part 21: model_init with prior_weights
# =============================================================================

test_that("susie() matches reference with model_init and prior_weights - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(39)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Non-uniform prior weights
  pw <- runif(p)
  pw <- pw / sum(pw)

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    prior_weights = pw,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   prior_weights = pw,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   prior_weights = pw,
                   estimate_prior_method = "optim")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init and prior_weights - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(39)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  pw <- runif(p)
  pw <- pw / sum(pw)

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    prior_weights = pw,
                    estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   prior_weights = pw,
                   estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   prior_weights = pw,
                   estimate_prior_method = "simple")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# L expansion with prior_weights
test_that("susie() matches reference with model_init, L expansion, and prior_weights - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(40)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  pw <- runif(p)
  pw <- pw / sum(pw)

  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    prior_weights = pw,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 8, model_init = dev_init,
                   prior_weights = pw,
                   estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 8, s_init = ref_init,
                   prior_weights = pw,
                   estimate_prior_method = "optim")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# =============================================================================
# Part 22: model_init with max_iter=1 (single iteration - check initialization)
# =============================================================================
# Running with max_iter=1 ensures we're testing the initialization path itself,
# not just the converged output.

test_that("susie() matches reference with model_init, max_iter=1 - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(41)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   max_iter = 1, estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   max_iter = 1, estimate_prior_method = "optim")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init, max_iter=1 - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(41)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 5, max_iter = 3,
                    estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 5, model_init = dev_init,
                   max_iter = 1, estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 5, s_init = ref_init,
                   max_iter = 1, estimate_prior_method = "simple")

  dev_result <- do.call(dev_env$env[["susie"]], dev_args)
  ref_result <- do.call(ref_env$env[["susie"]], ref_args)

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

# L expansion with max_iter=1
test_that("susie() matches reference with model_init, L expansion, max_iter=1 - optim", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(42)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    estimate_prior_method = "optim")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 8, model_init = dev_init,
                   max_iter = 1, estimate_prior_method = "optim")
  ref_args <- list(X = X, y = y, L = 8, s_init = ref_init,
                   max_iter = 1, estimate_prior_method = "optim")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})

test_that("susie() matches reference with model_init, L expansion, max_iter=1 - simple", {
  skip_if_no_reference()

  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  set.seed(42)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  init_args <- list(X = X, y = y, L = 3, max_iter = 3,
                    estimate_prior_method = "simple")
  dev_init <- do.call(dev_env$env[["susie"]], init_args)
  ref_init <- do.call(ref_env$env[["susie"]], init_args)

  dev_args <- list(X = X, y = y, L = 8, model_init = dev_init,
                   max_iter = 1, estimate_prior_method = "simple")
  ref_args <- list(X = X, y = y, L = 8, s_init = ref_init,
                   max_iter = 1, estimate_prior_method = "simple")

  dev_result <- suppressMessages(do.call(dev_env$env[["susie"]], dev_args))
  ref_result <- suppressMessages(do.call(ref_env$env[["susie"]], ref_args))

  expect_equal_susie_objects(dev_result, ref_result, tolerance = 1e-5)
})
