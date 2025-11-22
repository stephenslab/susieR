# Source helper functions
source(file.path("..", "helper_reference.R"), local = TRUE)

context("susie_auto reference comparison")

# =============================================================================
# REFERENCE TESTS FOR susie_auto()
# =============================================================================
#
# These tests compare the new susieR implementation against the reference
# package (stephenslab/susieR@1f9166c) to ensure results are identical.
#
# =============================================================================
# Part 1: Basic Tests with Default Parameters
# =============================================================================

test_that("susie_auto() matches reference with default parameters", {
  skip_if_no_reference()

  set.seed(1)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

test_that("susie_auto() matches reference with L_init=2", {
  skip_if_no_reference()

  set.seed(2)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L_init = 2)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

test_that("susie_auto() matches reference with L_init=5, L_max=10", {
  skip_if_no_reference()

  set.seed(3)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:3] <- c(2, -1.5, 1)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L_init = 5, L_max = 10)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

# =============================================================================
# Part 2: standardize Parameter
# =============================================================================

test_that("susie_auto() matches reference with standardize=FALSE", {
  skip_if_no_reference()

  set.seed(4)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, standardize = FALSE)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

test_that("susie_auto() matches reference with standardize=TRUE", {
  skip_if_no_reference()

  set.seed(5)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, standardize = TRUE)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

# =============================================================================
# Part 3: intercept Parameter
# =============================================================================

test_that("susie_auto() matches reference with intercept=FALSE", {
  skip_if_no_reference()

  set.seed(6)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, intercept = FALSE)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

test_that("susie_auto() matches reference with intercept=TRUE", {
  skip_if_no_reference()

  set.seed(7)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, intercept = TRUE)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

# =============================================================================
# Part 4: Tolerance Parameters
# =============================================================================

test_that("susie_auto() matches reference with init_tol=0.1", {
  skip_if_no_reference()

  set.seed(8)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, init_tol = 0.1)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

test_that("susie_auto() matches reference with tol=1e-3", {
  skip_if_no_reference()

  set.seed(9)
  n <- 100
  p <- 1000
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, tol = 1e-3)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

# =============================================================================
# Part 5: max_iter Parameter
# =============================================================================

test_that("susie_auto() matches reference with max_iter=50", {
  skip_if_no_reference()

  set.seed(10)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, max_iter = 50)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

test_that("susie_auto() matches reference with max_iter=200", {
  skip_if_no_reference()

  set.seed(11)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, max_iter = 200)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

# =============================================================================
# Part 6: Combined Parameter Tests
# =============================================================================

test_that("susie_auto() matches reference with standardize=FALSE, intercept=FALSE", {
  skip_if_no_reference()

  set.seed(12)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, standardize = FALSE, intercept = FALSE)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

test_that("susie_auto() matches reference with L_init=2, L_max=8, init_tol=0.5", {
  skip_if_no_reference()

  set.seed(13)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L_init = 2, L_max = 8, init_tol = 0.5)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

# =============================================================================
# Part 7: Edge Cases
# =============================================================================

test_that("susie_auto() matches reference with sparse signal", {
  skip_if_no_reference()

  set.seed(14)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[5] <- 3  # Only one effect
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

test_that("susie_auto() matches reference with dense signal", {
  skip_if_no_reference()

  set.seed(15)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:10] <- rnorm(10)  # Ten effects
  y <- as.vector(X %*% beta + rnorm(n))

  args <- list(X = X, y = y, L_init = 5, L_max = 20)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})

test_that("susie_auto() matches reference with high noise", {
  skip_if_no_reference()

  set.seed(16)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n, sd = 3))  # High noise

  args <- list(X = X, y = y)
  compare_to_reference("susie_auto", args, tolerance = 1e-5)
})
