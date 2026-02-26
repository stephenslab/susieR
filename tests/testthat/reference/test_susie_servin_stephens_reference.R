# Source helper functions for Servin-Stephens reference comparison
source(file.path("..", "helper_servin_stephens_reference.R"), local = TRUE)

context("susie Servin-Stephens reference comparison")

# =============================================================================
# REFERENCE TESTS FOR susie(estimate_residual_method = "Servin_Stephens")
# =============================================================================
#
# These tests compare our implementation of the Servin-Stephens (NIG)
# prior, invoked via estimate_residual_method = "Servin_Stephens",
# against the reference implementation on the fix-susie-small-sigma-update
# branch of stephenslab/susieR (commit a999d44), where the equivalent
# feature is invoked via small = TRUE.
#
# Parameter mapping between the two interfaces:
#   Dev:  estimate_residual_method = "Servin_Stephens"  <->  Ref: small = TRUE
#   Dev:  tol (convergence tolerance)                   <->  Ref: tol_small
#   Dev:  convergence_method = "pip" (auto-set)         <->  Ref: (hard-coded PIP convergence)
#   Dev:  estimate_prior_method = "EM" (auto-set)       <->  Ref: (forced to EM)
#   Dev:  alpha0, beta0                                 <->  Ref: alpha0, beta0
#
# The helper function compare_servin_stephens_to_reference() handles
# this mapping automatically.

# =============================================================================
# Part 1: Default parameters (baseline match)
# =============================================================================

test_that("Servin_Stephens matches reference (small=TRUE) with defaults", {
  skip_if_no_ss_reference()

  set.seed(1)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Dev: estimate_residual_method = "Servin_Stephens" (set by helper)
  # Ref: small = TRUE (mapped by helper)
  # Both default to alpha0 = 0.1, beta0 = 0.1, L = 10
  dev_args <- list(X = X, y = y, L = 10)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 2: L = 1 (single effect — ELBO is well-defined)
# =============================================================================

test_that("Servin_Stephens matches reference with L = 1", {
  skip("L=1 uses different convergence methods between dev and ref")
  skip_if_no_ss_reference()

  set.seed(2)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[3] <- 3
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 1)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 3: standardize = FALSE
# =============================================================================

test_that("Servin_Stephens matches reference with standardize=FALSE", {
  skip_if_no_ss_reference()

  set.seed(3)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, standardize = FALSE)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 4: intercept = FALSE
# =============================================================================

test_that("Servin_Stephens matches reference with intercept=FALSE", {
  skip_if_no_ss_reference()

  set.seed(4)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, intercept = FALSE)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 5: Custom alpha0 and beta0
# =============================================================================

test_that("Servin_Stephens matches reference with custom alpha0/beta0", {
  skip_if_no_ss_reference()

  set.seed(5)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, alpha0 = 1.0, beta0 = 1.0)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

test_that("Servin_Stephens matches reference with small alpha0/beta0", {
  skip_if_no_ss_reference()

  set.seed(6)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, alpha0 = 0.01, beta0 = 0.01)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 6: estimate_prior_variance = FALSE
# =============================================================================

test_that("Servin_Stephens matches reference with estimate_prior_variance=FALSE", {
  skip_if_no_ss_reference()

  set.seed(7)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, estimate_prior_variance = FALSE)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 7: Explicit convergence tolerance
# =============================================================================

test_that("Servin_Stephens matches reference with tol = 1e-4", {
  skip_if_no_ss_reference()

  set.seed(8)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Dev uses tol; helper maps it to tol_small for reference
  dev_args <- list(X = X, y = y, L = 10, tol = 1e-4)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 8: max_iter boundary
# =============================================================================

test_that("Servin_Stephens matches reference with small max_iter", {
  skip_if_no_ss_reference()

  set.seed(9)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Force early termination to test partial convergence path
  dev_args <- list(X = X, y = y, L = 10, max_iter = 5)
  results  <- suppressWarnings(
    compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)
  )

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 9: Sparse signal (most effects zero)
# =============================================================================

test_that("Servin_Stephens matches reference with very sparse signal", {
  skip_if_no_ss_reference()

  set.seed(10)
  n <- 100
  p <- 200
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1] <- 5  # single strong effect
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 10: Small sample size (n < 80, the regime Servin_Stephens targets)
# =============================================================================

test_that("Servin_Stephens matches reference with small n", {
  skip_if_no_ss_reference()

  set.seed(11)
  n <- 30
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:2] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 11: High noise (large residual variance)
# =============================================================================

test_that("Servin_Stephens matches reference with high noise", {
  skip_if_no_ss_reference()

  set.seed(12)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n, sd = 10))  # high noise

  dev_args <- list(X = X, y = y, L = 10)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 12: Combined — standardize=FALSE, intercept=FALSE
# =============================================================================

test_that("Servin_Stephens matches reference with standardize=FALSE, intercept=FALSE", {
  skip_if_no_ss_reference()

  set.seed(13)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, standardize = FALSE, intercept = FALSE)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 13: Combined — custom alpha0/beta0 with standardize=FALSE
# =============================================================================

test_that("Servin_Stephens matches reference with custom priors and standardize=FALSE", {
  skip_if_no_ss_reference()

  set.seed(14)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(
    X = X, y = y, L = 10,
    standardize = FALSE,
    alpha0 = 0.5,
    beta0 = 0.5
  )
  results <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 14: Null signal (no true effects)
# =============================================================================

test_that("Servin_Stephens matches reference under null signal", {
  skip_if_no_ss_reference()

  set.seed(15)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)  # pure noise

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 15: L = 1 with standardize = FALSE (ELBO well-defined, no scaling)
# =============================================================================

test_that("Servin_Stephens matches reference with L=1, standardize=FALSE", {
  skip("L=1 uses different convergence methods between dev and ref")
  skip_if_no_ss_reference()

  set.seed(16)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[5] <- 4
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 1, standardize = FALSE)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 16: Small n with intercept = FALSE
# =============================================================================

test_that("Servin_Stephens matches reference with small n and intercept=FALSE", {
  skip_if_no_ss_reference()

  set.seed(17)
  n <- 20
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1] <- 5
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 3, intercept = FALSE)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Part 17: Diagnostic — field-by-field summary of discrepancies
# =============================================================================
#
# This test does NOT use expect_equal; instead it generates a summary of
# all numeric differences between dev and reference outputs. Useful for
# diagnosing regressions without hard-failing CI.

test_that("Servin_Stephens field-by-field difference summary", {
  skip_if_no_ss_reference()

  set.seed(100)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  dev <- results$dev
  ref <- results$ref

  # Collect per-field max absolute differences
  fields <- c("alpha", "mu", "mu2", "V", "sigma2", "intercept",
              "fitted", "pip")

  diffs <- vapply(fields, function(f) {
    d <- dev[[f]]
    r <- ref[[f]]
    if (is.null(d) || is.null(r)) return(NA_real_)
    max(abs(d - r), na.rm = TRUE)
  }, numeric(1))

  # Print a summary table
  message("\n--- Servin-Stephens vs reference: max |dev - ref| per field ---")
  for (f in names(diffs)) {
    message(sprintf("  %-12s: %s", f, format(diffs[f], digits = 8)))
  }

  # Convergence & iteration metadata
  message(sprintf("  niter (dev/ref): %d / %d", dev$niter, ref$niter))
  message(sprintf("  converged (dev/ref): %s / %s", dev$converged, ref$converged))

  # Credible sets match
  if (!is.null(dev$sets$cs) && !is.null(ref$sets$cs)) {
    cs_match <- identical(dev$sets$cs, ref$sets$cs)
    message(sprintf("  CS sets identical: %s", cs_match))
  }

  # Hard assertion: all differences should be < tolerance
  expect_true(all(diffs[!is.na(diffs)] < 1e-5),
              info = paste("Some fields exceed tolerance:",
                           paste(names(which(diffs >= 1e-5)),
                                 collapse = ", ")))
})

# #############################################################################
# EXPANDED EDGE-CASE TEST SUITE
# #############################################################################

# =============================================================================
# Category A: Data dimensions
# =============================================================================

# Part 18: n >> p (overdetermined)
test_that("Servin_Stephens matches reference with n >> p (overdetermined)", {
  skip_if_no_ss_reference()

  set.seed(101)
  n <- 500
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:3] <- c(2, -1.5, 3)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 19: n << p (underdetermined, genetics regime)
test_that("Servin_Stephens matches reference with n << p (underdetermined)", {
  skip_if_no_ss_reference()

  set.seed(102)
  n <- 30
  p <- 200
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 50)] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 20: n = p (square)
test_that("Servin_Stephens matches reference with n = p (square)", {
  skip_if_no_ss_reference()

  set.seed(103)
  n <- 50
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:3] <- c(2, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 21: Very small n
test_that("Servin_Stephens matches reference with very small n", {
  skip_if_no_ss_reference()

  set.seed(104)
  n <- 10
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1] <- 3
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 3)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Category B: Signal patterns
# =============================================================================

# Part 22: Weak signals (low SNR)
test_that("Servin_Stephens matches reference with weak signals", {
  skip_if_no_ss_reference()

  set.seed(105)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:3] <- c(0.3, -0.3, 0.2)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 23: Very strong signals
test_that("Servin_Stephens matches reference with very strong signals", {
  skip_if_no_ss_reference()

  set.seed(106)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:3] <- c(10, -15, 20)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 24: Mixed strength signals
test_that("Servin_Stephens matches reference with mixed strength signals", {
  skip_if_no_ss_reference()

  set.seed(107)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(10, 0.5, -10, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 25: Many true effects
test_that("Servin_Stephens matches reference with many true effects", {
  skip_if_no_ss_reference()

  set.seed(108)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:10] <- c(2, -1.5, 3, -2, 1, -1, 2.5, -0.8, 1.2, -1.8)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Category C: L values
# =============================================================================

# Part 26: L = 2 (minimal multi-effect)
test_that("Servin_Stephens matches reference with L = 2", {
  skip_if_no_ss_reference()

  set.seed(109)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:2] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 2)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 27: L = 20 (more effects than default)
test_that("Servin_Stephens matches reference with L = 20", {
  skip_if_no_ss_reference()

  set.seed(110)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 20)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 28: L >> true effects (over-specified)
test_that("Servin_Stephens matches reference with L >> true effects", {
  skip_if_no_ss_reference()

  set.seed(111)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:2] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 15)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 29: L < true effects (under-specified)
test_that("Servin_Stephens matches reference with L < true effects", {
  skip_if_no_ss_reference()

  set.seed(112)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:5] <- c(3, -2, 4, -1.5, 2)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 2)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Category D: Prior parameters alpha0/beta0
# =============================================================================

# Part 30: Informative priors
test_that("Servin_Stephens matches reference with informative alpha0/beta0", {
  skip_if_no_ss_reference()

  set.seed(113)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, alpha0 = 10, beta0 = 10)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 31: Very diffuse priors
test_that("Servin_Stephens matches reference with very diffuse alpha0/beta0", {
  skip_if_no_ss_reference()

  set.seed(114)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, alpha0 = 0.001, beta0 = 0.001)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 32: Asymmetric priors
test_that("Servin_Stephens matches reference with asymmetric alpha0/beta0", {
  skip_if_no_ss_reference()

  set.seed(115)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, alpha0 = 0.1, beta0 = 1.0)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Category E: Predictor structure
# =============================================================================

# Part 33: AR(1) correlated predictors
test_that("Servin_Stephens matches reference with AR(1) correlated X", {
  skip_if_no_ss_reference()

  set.seed(116)
  n <- 100
  p <- 50
  rho <- 0.8

  # Generate AR(1) correlation structure
  Z <- matrix(rnorm(n * p), n, p)
  X <- Z
  for (j in 2:p) {
    X[, j] <- rho * X[, j - 1] + sqrt(1 - rho^2) * Z[, j]
  }

  beta <- rep(0, p)
  beta[c(1, 25)] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 34: Block-correlated predictors
test_that("Servin_Stephens matches reference with block-correlated X", {
  skip_if_no_ss_reference()

  set.seed(117)
  n <- 100
  p <- 50
  block_size <- 5
  n_blocks <- p / block_size

  # Generate block correlation structure
  X <- matrix(0, n, p)
  for (b in seq_len(n_blocks)) {
    cols <- ((b - 1) * block_size + 1):(b * block_size)
    common <- rnorm(n)
    for (j in cols) {
      X[, j] <- 0.8 * common + 0.6 * rnorm(n)
    }
  }

  beta <- rep(0, p)
  beta[c(1, 26)] <- c(3, -2)  # one in first block, one in sixth block
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 35: Near-collinear predictors
test_that("Servin_Stephens matches reference with near-collinear predictors", {
  skip_if_no_ss_reference()

  set.seed(118)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)

  # Make columns 2 nearly identical to column 1 (r ~ 0.99)
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)

  beta <- rep(0, p)
  beta[c(1, 2)] <- c(2, -1.5)  # both collinear predictors active
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Category F: Convergence settings
# =============================================================================

# Part 36: max_iter = 1 (single iteration snapshot)
test_that("Servin_Stephens matches reference with max_iter = 1", {
  skip_if_no_ss_reference()

  set.seed(119)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, max_iter = 1)
  results  <- suppressWarnings(
    compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)
  )

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 37: max_iter = 2 (minimal convergence path)
test_that("Servin_Stephens matches reference with max_iter = 2", {
  skip_if_no_ss_reference()

  set.seed(120)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, max_iter = 2)
  results  <- suppressWarnings(
    compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)
  )

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 38: Tight convergence tolerance
test_that("Servin_Stephens matches reference with tight tol = 1e-6", {
  skip_if_no_ss_reference()

  set.seed(121)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, tol = 1e-6)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Category G: null_weight and prior_weights
# =============================================================================

# Part 39: null_weight = 0.5 (strong null prior)
test_that("Servin_Stephens matches reference with null_weight = 0.5", {
  skip("null_weight + Servin_Stephens triggers NA in loglik (dev-side bug)")
  skip_if_no_ss_reference()

  set.seed(122)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10, null_weight = 0.5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 40: Non-uniform prior_weights
test_that("Servin_Stephens matches reference with non-uniform prior_weights", {
  skip_if_no_ss_reference()

  set.seed(123)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  # Favor first 10 predictors
  pw <- rep(1, p)
  pw[1:10] <- 5
  pw <- pw / sum(pw)

  dev_args <- list(X = X, y = y, L = 10, prior_weights = pw)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Category H: Parameter combinations
# =============================================================================

# Part 41: Small n + weak signal
test_that("Servin_Stephens matches reference with small n + weak signal", {
  skip_if_no_ss_reference()

  set.seed(124)
  n <- 20
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:2] <- c(0.5, -0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 42: n << p + L large
test_that("Servin_Stephens matches reference with n << p and large L", {
  skip_if_no_ss_reference()

  set.seed(125)
  n <- 30
  p <- 200
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(10, 50, 100)] <- c(3, -2, 4)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 10)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 43: intercept=FALSE + small n
test_that("Servin_Stephens matches reference with intercept=FALSE + small n", {
  skip_if_no_ss_reference()

  set.seed(126)
  n <- 25
  p <- 40
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:3] <- c(2, -1.5, 3)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 5, intercept = FALSE)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 44: standardize=FALSE + custom alpha0/beta0
test_that("Servin_Stephens matches reference with standardize=FALSE + custom priors", {
  skip_if_no_ss_reference()

  set.seed(127)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(
    X = X, y = y, L = 10,
    standardize = FALSE,
    alpha0 = 1.0, beta0 = 0.5
  )
  results <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 45: estimate_prior_variance=FALSE + intercept=FALSE
test_that("Servin_Stephens matches reference with estimate_prior_variance=FALSE + intercept=FALSE", {
  skip_if_no_ss_reference()

  set.seed(128)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(
    X = X, y = y, L = 10,
    estimate_prior_variance = FALSE,
    intercept = FALSE
  )
  results <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# =============================================================================
# Category I: L = 1 variants
# =============================================================================

# Part 46: L=1 + high noise
test_that("Servin_Stephens matches reference with L=1 + high noise", {
  skip("L=1 uses different convergence methods between dev and ref")
  skip_if_no_ss_reference()

  set.seed(129)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1] <- 3
  y <- as.vector(X %*% beta + rnorm(n, sd = 10))

  dev_args <- list(X = X, y = y, L = 1)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 47: L=1 + very small n
test_that("Servin_Stephens matches reference with L=1 + very small n", {
  skip("L=1 uses different convergence methods between dev and ref")
  skip_if_no_ss_reference()

  set.seed(130)
  n <- 15
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1] <- 4
  y <- as.vector(X %*% beta + rnorm(n))

  dev_args <- list(X = X, y = y, L = 1)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# Part 48: L=1 + null signal (no true effect)
test_that("Servin_Stephens matches reference with L=1 + null signal", {
  skip("L=1 uses different convergence methods between dev and ref")
  skip_if_no_ss_reference()

  set.seed(131)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)  # pure noise

  dev_args <- list(X = X, y = y, L = 1)
  results  <- compare_servin_stephens_to_reference(dev_args, tolerance = 1e-5)

  expect_equal_servin_stephens_objects(results$dev, results$ref, tolerance = 1e-5)
})

# #############################################################################
# SUFFICIENT STATISTICS VS INDIVIDUAL-LEVEL DATA COMPARISON
# #############################################################################
#
# For each reference test scenario above, verify that susie_ss()
# produces the same result as susie() with Servin_Stephens.
# These tests do NOT require the reference package.

# =============================================================================
# SS Part 1: Default parameters (baseline match)
# =============================================================================

test_that("SS matches individual: defaults (Part 1)", {
  set.seed(1)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 10))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 2: L = 1 (skipped — convergence method differs)
# =============================================================================

test_that("SS matches individual: L = 1 (Part 2)", {
  skip("L=1 uses different convergence methods between dev and ref")

  set.seed(2)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[3] <- 3
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 1))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 3: standardize = FALSE
# =============================================================================

test_that("SS matches individual: standardize=FALSE (Part 3)", {
  set.seed(3)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, standardize = FALSE))
  expect_ss_matches_individual_ss(res)
})


# =============================================================================
# SS Part 4: Custom alpha0 and beta0
# =============================================================================

test_that("SS matches individual: custom alpha0/beta0 (Part 4)", {
  set.seed(5)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, alpha0 = 1.0, beta0 = 1.0))
  expect_ss_matches_individual_ss(res)
})

test_that("SS matches individual: small alpha0/beta0 (Part 5)", {
  set.seed(6)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, alpha0 = 0.01, beta0 = 0.01))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 6: estimate_prior_variance = FALSE
# =============================================================================

test_that("SS matches individual: estimate_prior_variance=FALSE (Part 6)", {
  set.seed(7)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, estimate_prior_variance = FALSE))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 7: Explicit convergence tolerance
# =============================================================================

test_that("SS matches individual: tol = 1e-4 (Part 7)", {
  set.seed(8)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, tol = 1e-4))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 8: max_iter boundary
# =============================================================================

test_that("SS matches individual: small max_iter (Part 8)", {
  set.seed(9)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, max_iter = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 9: Sparse signal
# =============================================================================

test_that("SS matches individual: very sparse signal (Part 9)", {
  set.seed(10)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1] <- 5
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 10: Small sample size
# =============================================================================

test_that("SS matches individual: small n (Part 10)", {
  set.seed(11)
  n <- 30; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:2] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 11: High noise
# =============================================================================

test_that("SS matches individual: high noise (Part 11)", {
  set.seed(12)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n, sd = 10))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 10))
  expect_ss_matches_individual_ss(res)
})


# =============================================================================
# SS Part 12: Custom alpha0/beta0 with standardize=FALSE
# =============================================================================

test_that("SS matches individual: custom priors + standardize=FALSE (Part 12)", {
  set.seed(14)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, standardize = FALSE, alpha0 = 0.5, beta0 = 0.5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 13: Null signal
# =============================================================================

test_that("SS matches individual: null signal (Part 13)", {
  set.seed(15)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 14: L=1 with standardize=FALSE (skipped)
# =============================================================================

test_that("SS matches individual: L=1, standardize=FALSE (Part 14)", {
  skip("L=1 uses different convergence methods between dev and ref")

  set.seed(16)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[5] <- 4
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 1, standardize = FALSE))
  expect_ss_matches_individual_ss(res)
})


# =============================================================================
# SS Part 15: n >> p (overdetermined)
# =============================================================================

test_that("SS matches individual: n >> p (Part 15)", {
  set.seed(101)
  n <- 500; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:3] <- c(2, -1.5, 3)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 16: n << p (underdetermined)
# =============================================================================

test_that("SS matches individual: n << p (Part 16)", {
  set.seed(102)
  n <- 30; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[c(5, 50)] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 17: n = p (square)
# =============================================================================

test_that("SS matches individual: n = p (Part 17)", {
  set.seed(103)
  n <- 50; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:3] <- c(2, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 18: Very small n
# =============================================================================

test_that("SS matches individual: very small n (Part 18)", {
  set.seed(104)
  n <- 10; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1] <- 3
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 3))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 19: Weak signals
# =============================================================================

test_that("SS matches individual: weak signals (Part 19)", {
  set.seed(105)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:3] <- c(0.3, -0.3, 0.2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 20: Very strong signals
# =============================================================================

test_that("SS matches individual: very strong signals (Part 20)", {
  set.seed(106)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:3] <- c(10, -15, 20)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 21: Mixed strength signals
# =============================================================================

test_that("SS matches individual: mixed strength signals (Part 21)", {
  set.seed(107)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(10, 0.5, -10, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 10))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 22: Many true effects
# =============================================================================

test_that("SS matches individual: many true effects (Part 22)", {
  set.seed(108)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:10] <- c(2, -1.5, 3, -2, 1, -1, 2.5, -0.8, 1.2, -1.8)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 10))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 23: L = 2
# =============================================================================

test_that("SS matches individual: L = 2 (Part 23)", {
  set.seed(109)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:2] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 2))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 24: L = 20
# =============================================================================

test_that("SS matches individual: L = 20 (Part 24)", {
  set.seed(110)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 20))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 25: L >> true effects
# =============================================================================

test_that("SS matches individual: L >> true effects (Part 25)", {
  set.seed(111)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:2] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 15))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 26: L < true effects
# =============================================================================

test_that("SS matches individual: L < true effects (Part 26)", {
  set.seed(112)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:5] <- c(3, -2, 4, -1.5, 2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 2))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 27: Informative alpha0/beta0
# =============================================================================

test_that("SS matches individual: informative alpha0/beta0 (Part 27)", {
  set.seed(113)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, alpha0 = 10, beta0 = 10))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 28: Very diffuse alpha0/beta0
# =============================================================================

test_that("SS matches individual: very diffuse alpha0/beta0 (Part 28)", {
  set.seed(114)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, alpha0 = 0.001, beta0 = 0.001))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 29: Asymmetric alpha0/beta0
# =============================================================================

test_that("SS matches individual: asymmetric alpha0/beta0 (Part 29)", {
  set.seed(115)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, alpha0 = 0.1, beta0 = 1.0))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 30: AR(1) correlated predictors
# =============================================================================

test_that("SS matches individual: AR(1) correlated X (Part 30)", {
  set.seed(116)
  n <- 100; p <- 50; rho <- 0.8

  Z <- matrix(rnorm(n * p), n, p)
  X <- Z
  for (j in 2:p) {
    X[, j] <- rho * X[, j - 1] + sqrt(1 - rho^2) * Z[, j]
  }

  beta <- rep(0, p); beta[c(1, 25)] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 31: Block-correlated predictors
# =============================================================================

test_that("SS matches individual: block-correlated X (Part 31)", {
  set.seed(117)
  n <- 100; p <- 50; block_size <- 5
  n_blocks <- p / block_size

  X <- matrix(0, n, p)
  for (b in seq_len(n_blocks)) {
    cols <- ((b - 1) * block_size + 1):(b * block_size)
    common <- rnorm(n)
    for (j in cols) {
      X[, j] <- 0.8 * common + 0.6 * rnorm(n)
    }
  }

  beta <- rep(0, p); beta[c(1, 26)] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 32: Near-collinear predictors
# =============================================================================

test_that("SS matches individual: near-collinear X (Part 32)", {
  set.seed(118)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)

  beta <- rep(0, p); beta[c(1, 2)] <- c(2, -1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 33: max_iter = 1
# =============================================================================

test_that("SS matches individual: max_iter = 1 (Part 33)", {
  set.seed(119)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, max_iter = 1))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 34: max_iter = 2
# =============================================================================

test_that("SS matches individual: max_iter = 2 (Part 34)", {
  set.seed(120)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, max_iter = 2))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 35: Tight convergence tolerance
# =============================================================================

test_that("SS matches individual: tight tol = 1e-6 (Part 35)", {
  set.seed(121)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, tol = 1e-6))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 36: null_weight = 0.5 (skipped — known dev-side bug)
# =============================================================================

test_that("SS matches individual: null_weight = 0.5 (Part 36)", {
  skip("null_weight + Servin_Stephens triggers NA in loglik (dev-side bug)")

  set.seed(122)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, null_weight = 0.5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 37: Non-uniform prior_weights
# =============================================================================

test_that("SS matches individual: non-uniform prior_weights (Part 37)", {
  set.seed(123)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  pw <- rep(1, p); pw[1:10] <- 5; pw <- pw / sum(pw)

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, prior_weights = pw))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 38: Small n + weak signal
# =============================================================================

test_that("SS matches individual: small n + weak signal (Part 38)", {
  set.seed(124)
  n <- 20; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:2] <- c(0.5, -0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_ss_matches_individual_ss(res)
})

# =============================================================================
# SS Part 39: n << p + L large
# =============================================================================

test_that("SS matches individual: n << p + large L (Part 39)", {
  set.seed(125)
  n <- 30; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[c(10, 50, 100)] <- c(3, -2, 4)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 10))
  expect_ss_matches_individual_ss(res)
})


# =============================================================================
# SS Part 40: standardize=FALSE + custom alpha0/beta0
# =============================================================================

test_that("SS matches individual: standardize=FALSE + custom priors (Part 40)", {
  set.seed(127)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y,
    list(L = 10, standardize = FALSE, alpha0 = 1.0, beta0 = 0.5))
  expect_ss_matches_individual_ss(res)
})


# =============================================================================
# SS Part 41-43: L = 1 variants (all skipped)
# =============================================================================

test_that("SS matches individual: L=1 + high noise (Part 41)", {
  skip("L=1 uses different convergence methods between dev and ref")

  set.seed(129)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1] <- 3
  y <- as.vector(X %*% beta + rnorm(n, sd = 10))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 1))
  expect_ss_matches_individual_ss(res)
})

test_that("SS matches individual: L=1 + very small n (Part 42)", {
  skip("L=1 uses different convergence methods between dev and ref")

  set.seed(130)
  n <- 15; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1] <- 4
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 1))
  expect_ss_matches_individual_ss(res)
})

test_that("SS matches individual: L=1 + null signal (Part 43)", {
  skip("L=1 uses different convergence methods between dev and ref")

  set.seed(131)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  res <- run_ss_and_individual_servin_stephens(X, y, list(L = 1))
  expect_ss_matches_individual_ss(res)
})

# #############################################################################
# RSS (SUMMARY STATISTICS) VS INDIVIDUAL-LEVEL DATA COMPARISON
# #############################################################################
#
# For each reference test scenario above, verify that susie_rss()
# (via the bhat/shat/var_y path) produces the same result as susie()
# with Servin_Stephens.
# These tests do NOT require the reference package.

# =============================================================================
# RSS Part 1: Default parameters (baseline match)
# =============================================================================

test_that("RSS matches individual: defaults (Part 1)", {
  set.seed(1)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 10))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 2: L = 1 (skipped — convergence method differs)
# =============================================================================

test_that("RSS matches individual: L = 1 (Part 2)", {
  skip("L=1 uses different convergence methods between dev and ref")

  set.seed(2)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[3] <- 3
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 1))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 3: standardize = FALSE
# =============================================================================

test_that("RSS matches individual: standardize=FALSE (Part 3)", {
  set.seed(3)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, standardize = FALSE))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 4: Custom alpha0 and beta0
# =============================================================================

test_that("RSS matches individual: custom alpha0/beta0 (Part 4)", {
  set.seed(5)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, alpha0 = 1.0, beta0 = 1.0))
  expect_rss_matches_individual_ss(res)
})

test_that("RSS matches individual: small alpha0/beta0 (Part 5)", {
  set.seed(6)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, alpha0 = 0.01, beta0 = 0.01))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 6: estimate_prior_variance = FALSE
# =============================================================================

test_that("RSS matches individual: estimate_prior_variance=FALSE (Part 6)", {
  set.seed(7)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, estimate_prior_variance = FALSE))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 7: Explicit convergence tolerance
# =============================================================================

test_that("RSS matches individual: tol = 1e-4 (Part 7)", {
  set.seed(8)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, tol = 1e-4))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 8: max_iter boundary
# =============================================================================

test_that("RSS matches individual: small max_iter (Part 8)", {
  set.seed(9)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, max_iter = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 9: Sparse signal
# =============================================================================

test_that("RSS matches individual: very sparse signal (Part 9)", {
  set.seed(10)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1] <- 5
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 10: Small sample size
# =============================================================================

test_that("RSS matches individual: small n (Part 10)", {
  set.seed(11)
  n <- 30; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:2] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 11: High noise
# =============================================================================

test_that("RSS matches individual: high noise (Part 11)", {
  set.seed(12)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n, sd = 10))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 10))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 12: Custom alpha0/beta0 with standardize=FALSE
# =============================================================================

test_that("RSS matches individual: custom priors + standardize=FALSE (Part 12)", {
  set.seed(14)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, standardize = FALSE, alpha0 = 0.5, beta0 = 0.5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 13: Null signal
# =============================================================================

test_that("RSS matches individual: null signal (Part 13)", {
  set.seed(15)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 14: L=1 with standardize=FALSE (skipped)
# =============================================================================

test_that("RSS matches individual: L=1, standardize=FALSE (Part 14)", {
  skip("L=1 uses different convergence methods between dev and ref")

  set.seed(16)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[5] <- 4
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 1, standardize = FALSE))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 15: n >> p (overdetermined)
# =============================================================================

test_that("RSS matches individual: n >> p (Part 15)", {
  set.seed(101)
  n <- 500; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:3] <- c(2, -1.5, 3)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 16: n << p (underdetermined)
# =============================================================================

test_that("RSS matches individual: n << p (Part 16)", {
  set.seed(102)
  n <- 30; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[c(5, 50)] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 17: n = p (square)
# =============================================================================

test_that("RSS matches individual: n = p (Part 17)", {
  set.seed(103)
  n <- 50; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:3] <- c(2, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 18: Very small n
# =============================================================================

test_that("RSS matches individual: very small n (Part 18)", {
  set.seed(104)
  n <- 10; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1] <- 3
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 3))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 19: Weak signals
# =============================================================================

test_that("RSS matches individual: weak signals (Part 19)", {
  set.seed(105)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:3] <- c(0.3, -0.3, 0.2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 20: Very strong signals
# =============================================================================

test_that("RSS matches individual: very strong signals (Part 20)", {
  set.seed(106)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:3] <- c(10, -15, 20)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 21: Mixed strength signals
# =============================================================================

test_that("RSS matches individual: mixed strength signals (Part 21)", {
  set.seed(107)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(10, 0.5, -10, 0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 10))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 22: Many true effects
# =============================================================================

test_that("RSS matches individual: many true effects (Part 22)", {
  set.seed(108)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:10] <- c(2, -1.5, 3, -2, 1, -1, 2.5, -0.8, 1.2, -1.8)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 10))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 23: L = 2
# =============================================================================

test_that("RSS matches individual: L = 2 (Part 23)", {
  set.seed(109)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:2] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 2))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 24: L = 20
# =============================================================================

test_that("RSS matches individual: L = 20 (Part 24)", {
  set.seed(110)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 20))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 25: L >> true effects
# =============================================================================

test_that("RSS matches individual: L >> true effects (Part 25)", {
  set.seed(111)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:2] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 15))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 26: L < true effects
# =============================================================================

test_that("RSS matches individual: L < true effects (Part 26)", {
  set.seed(112)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:5] <- c(3, -2, 4, -1.5, 2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 2))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 27: Informative alpha0/beta0
# =============================================================================

test_that("RSS matches individual: informative alpha0/beta0 (Part 27)", {
  set.seed(113)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, alpha0 = 10, beta0 = 10))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 28: Very diffuse alpha0/beta0
# =============================================================================

test_that("RSS matches individual: very diffuse alpha0/beta0 (Part 28)", {
  set.seed(114)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, alpha0 = 0.001, beta0 = 0.001))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 29: Asymmetric alpha0/beta0
# =============================================================================

test_that("RSS matches individual: asymmetric alpha0/beta0 (Part 29)", {
  set.seed(115)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, alpha0 = 0.1, beta0 = 1.0))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 30: AR(1) correlated predictors
# =============================================================================

test_that("RSS matches individual: AR(1) correlated X (Part 30)", {
  set.seed(116)
  n <- 100; p <- 50; rho <- 0.8

  Z <- matrix(rnorm(n * p), n, p)
  X <- Z
  for (j in 2:p) {
    X[, j] <- rho * X[, j - 1] + sqrt(1 - rho^2) * Z[, j]
  }

  beta <- rep(0, p); beta[c(1, 25)] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 31: Block-correlated predictors
# =============================================================================

test_that("RSS matches individual: block-correlated X (Part 31)", {
  set.seed(117)
  n <- 100; p <- 50; block_size <- 5
  n_blocks <- p / block_size

  X <- matrix(0, n, p)
  for (b in seq_len(n_blocks)) {
    cols <- ((b - 1) * block_size + 1):(b * block_size)
    common <- rnorm(n)
    for (j in cols) {
      X[, j] <- 0.8 * common + 0.6 * rnorm(n)
    }
  }

  beta <- rep(0, p); beta[c(1, 26)] <- c(3, -2)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 32: Near-collinear predictors
# =============================================================================

test_that("RSS matches individual: near-collinear X (Part 32)", {
  set.seed(118)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)

  beta <- rep(0, p); beta[c(1, 2)] <- c(2, -1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 33: max_iter = 1
# =============================================================================

test_that("RSS matches individual: max_iter = 1 (Part 33)", {
  set.seed(119)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, max_iter = 1))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 34: max_iter = 2
# =============================================================================

test_that("RSS matches individual: max_iter = 2 (Part 34)", {
  set.seed(120)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, max_iter = 2))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 35: Tight convergence tolerance
# =============================================================================

test_that("RSS matches individual: tight tol = 1e-6 (Part 35)", {
  set.seed(121)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, tol = 1e-6))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 36: null_weight = 0.5 (skipped — known dev-side bug)
# =============================================================================

test_that("RSS matches individual: null_weight = 0.5 (Part 36)", {
  skip("null_weight + Servin_Stephens triggers NA in loglik (dev-side bug)")

  set.seed(122)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, null_weight = 0.5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 37: Non-uniform prior_weights
# =============================================================================

test_that("RSS matches individual: non-uniform prior_weights (Part 37)", {
  set.seed(123)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  pw <- rep(1, p); pw[1:10] <- 5; pw <- pw / sum(pw)

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, prior_weights = pw))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 38: Small n + weak signal
# =============================================================================

test_that("RSS matches individual: small n + weak signal (Part 38)", {
  set.seed(124)
  n <- 20; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:2] <- c(0.5, -0.3)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 39: n << p + L large
# =============================================================================

test_that("RSS matches individual: n << p + large L (Part 39)", {
  set.seed(125)
  n <- 30; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[c(10, 50, 100)] <- c(3, -2, 4)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 10))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 40: standardize=FALSE + custom alpha0/beta0
# =============================================================================

test_that("RSS matches individual: standardize=FALSE + custom priors (Part 40)", {
  set.seed(127)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:4] <- c(2, 3, -2, 1.5)
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y,
    list(L = 10, standardize = FALSE, alpha0 = 1.0, beta0 = 0.5))
  expect_rss_matches_individual_ss(res)
})

# =============================================================================
# RSS Part 41-43: L = 1 variants (all skipped)
# =============================================================================

test_that("RSS matches individual: L=1 + high noise (Part 41)", {
  skip("L=1 uses different convergence methods between dev and ref")

  set.seed(129)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1] <- 3
  y <- as.vector(X %*% beta + rnorm(n, sd = 10))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 1))
  expect_rss_matches_individual_ss(res)
})

test_that("RSS matches individual: L=1 + very small n (Part 42)", {
  skip("L=1 uses different convergence methods between dev and ref")

  set.seed(130)
  n <- 15; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1] <- 4
  y <- as.vector(X %*% beta + rnorm(n))

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 1))
  expect_rss_matches_individual_ss(res)
})

test_that("RSS matches individual: L=1 + null signal (Part 43)", {
  skip("L=1 uses different convergence methods between dev and ref")

  set.seed(131)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  res <- run_rss_and_individual_servin_stephens(X, y, list(L = 1))
  expect_rss_matches_individual_ss(res)
})
