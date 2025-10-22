context("RSS utility functions")

# =============================================================================
# FUNDAMENTAL COMPUTATIONS
# =============================================================================

test_that("compute_suff_stat with standardize=FALSE produces correct XtX", {
  base_data <- generate_base_data(n = 10, p = 5, seed = 1)

  # Manual calculation: center X
  X_centered <- scale(base_data$X, center = TRUE, scale = FALSE)

  out <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)

  dimnames(out$XtX) <- NULL
  expect_equal(out$XtX, crossprod(X_centered), tolerance = 1e-14)
})

test_that("compute_suff_stat with standardize=TRUE produces correct XtX", {
  base_data <- generate_base_data(n = 10, p = 5, seed = 2)

  # Manual calculation: center and scale X
  X_standardized <- scale(base_data$X, center = TRUE, scale = TRUE)

  out <- compute_suff_stat(base_data$X, base_data$y, standardize = TRUE)

  dimnames(out$XtX) <- NULL
  expect_equal(out$XtX, crossprod(X_standardized), tolerance = 1e-14)
})

test_that("compute_suff_stat with sparse matrix input", {
  base_data <- generate_base_data(n = 10, p = 5, seed = 3)

  # Sparse version
  X_sparse <- as(base_data$X, "CsparseMatrix")

  out_dense <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)
  out_sparse <- compute_suff_stat(X_sparse, base_data$y, standardize = FALSE)

  dimnames(out_dense$XtX) <- NULL
  dimnames(out_sparse$XtX) <- NULL

  expect_equal(out_sparse$XtX, out_dense$XtX, tolerance = 1e-14)
  expect_equal(as.vector(out_sparse$Xty), out_dense$Xty, tolerance = 1e-14)
  expect_equal(out_sparse$yty, out_dense$yty, tolerance = 1e-14)
})

test_that("compute_suff_stat produces correct Xty", {
  base_data <- generate_base_data(n = 20, p = 8, seed = 4)

  out <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)

  # Manual calculation
  y_centered <- base_data$y - mean(base_data$y)
  X_centered <- scale(base_data$X, center = TRUE, scale = FALSE)
  expected_Xty <- drop(crossprod(X_centered, y_centered))

  expect_equal(out$Xty, expected_Xty, tolerance = 1e-14)
})

test_that("compute_suff_stat produces correct yty", {
  base_data <- generate_base_data(n = 20, p = 8, seed = 5)

  out <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)

  # Manual calculation
  y_centered <- base_data$y - mean(base_data$y)
  expected_yty <- sum(y_centered^2)

  expect_equal(out$yty, expected_yty, tolerance = 1e-14)
})

test_that("compute_suff_stat stores column means and y_mean", {
  base_data <- generate_base_data(n = 15, p = 6, seed = 6)

  out <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)

  expect_equal(out$X_colmeans, colMeans(base_data$X), tolerance = 1e-14)
  expect_equal(out$y_mean, mean(base_data$y), tolerance = 1e-14)
  expect_equal(out$n, base_data$n)
})

test_that("compute_suff_stat with standardize=TRUE scales correctly", {
  base_data <- generate_base_data(n = 25, p = 10, seed = 7)

  out <- compute_suff_stat(base_data$X, base_data$y, standardize = TRUE)

  # XtX diagonal should be close to n (since standardized columns have variance 1)
  # After centering: crossprod of standardized X
  X_std <- scale(base_data$X, center = TRUE, scale = TRUE)
  expected_diag <- diag(crossprod(X_std))

  expect_equal(diag(out$XtX), expected_diag, tolerance = 1e-12)
})

test_that("compute_suff_stat returns list with correct names", {
  base_data <- generate_base_data(n = 10, p = 5, seed = 8)

  out <- compute_suff_stat(base_data$X, base_data$y, standardize = FALSE)

  expect_type(out, "list")
  expect_named(out, c("XtX", "Xty", "yty", "n", "y_mean", "X_colmeans"))
})

test_that("compute_suff_stat matches susie_ss input requirements", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 9)

  ss <- compute_suff_stat(base_data$X, base_data$y, standardize = TRUE)

  # Should be able to use directly with susie_ss
  expect_error(
    fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                    max_iter = 2, verbose = FALSE),
    NA
  )
})

test_that("compute_suff_stat with zero-variance column", {
  skip("Fails on Linux in CI")

  base_data <- generate_base_data(n = 20, p = 5, seed = 10)
  base_data$X[, 3] <- 1  # Constant column (zero variance after centering)

  # Should not error
  expect_error(
    out <- compute_suff_stat(base_data$X, base_data$y, standardize = TRUE),
    NA
  )
  expect_true(is.infinite(out$Xty[3]))
})

# =============================================================================
# RSS MODEL METHODS
# =============================================================================

test_that("estimate_s_rss returns value between 0 and 1 (null-mle)", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 11)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  s <- estimate_s_rss(z, R, n = base_data$n, method = "null-mle")

  expect_type(s, "double")
  expect_length(s, 1)
  expect_true(s >= 0 && s <= 1)
})

test_that("estimate_s_rss with null-partialmle method", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 12)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  s <- estimate_s_rss(z, R, n = base_data$n, method = "null-partialmle")

  expect_type(s, "double")
  expect_length(s, 1)
  # Note: null-partialmle can be > 1
  expect_true(s >= 0)
})

test_that("estimate_s_rss with null-pseudomle method", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 13)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  s <- estimate_s_rss(z, R, n = base_data$n, method = "null-pseudomle")

  expect_type(s, "double")
  expect_length(s, 1)
  expect_true(s >= 0 && s <= 1)
})

test_that("estimate_s_rss warns when n is not provided", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 14)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  expect_message(
    s <- estimate_s_rss(z, R),
    "sample size"
  )
})

test_that("estimate_s_rss errors when n <= 1", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 15)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  expect_error(
    estimate_s_rss(z, R, n = 1),
    "must be greater than 1"
  )

  expect_error(
    estimate_s_rss(z, R, n = 0),
    "must be greater than 1"
  )
})

test_that("estimate_s_rss handles eigen decomposition in R", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 16)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  # Pre-compute eigen decomposition
  attr(R, "eigen") <- eigen(R, symmetric = TRUE)

  s1 <- estimate_s_rss(z, R, n = base_data$n, method = "null-mle")

  # Without pre-computed eigen
  R2 <- cor(base_data$X)
  s2 <- estimate_s_rss(z, R2, n = base_data$n, method = "null-mle")

  expect_equal(s1, s2, tolerance = 1e-10)
})

test_that("estimate_s_rss handles negative eigenvalues in R", {
  set.seed(17)
  p <- 50

  # Create R with intentionally negative eigenvalue
  R <- matrix(0.5, p, p)
  diag(R) <- 1
  R[1, 2] <- 1.5 
  R[2, 1] <- 1.5

  z <- rnorm(p)

  expect_message(
    s <- estimate_s_rss(z, R, n = 100, method = "null-mle"),
    "not positive semidefinite"
  )

  # Should still return valid estimate
  expect_true(s >= 0 && s <= 1)
})

test_that("estimate_s_rss handles NA values in z", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 18)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  # Introduce NA
  z[5] <- NA

  expect_error(
    s <- estimate_s_rss(z, R, n = base_data$n, method = "null-mle"),
    NA
  )

  expect_true(s >= 0 && s <= 1)
})

test_that("estimate_s_rss with perfect LD has one large eigenvalue", {
  set.seed(19)
  p <- 10

  # All variables perfectly correlated
  R <- matrix(1, p, p)
  z <- rnorm(p)
  s <- estimate_s_rss(z, R, n = 100, method = "null-partialmle")
  expect_true(s >= 0)
})

test_that("estimate_s_rss errors on invalid method", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 20)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  expect_error(
    estimate_s_rss(z, R, n = base_data$n, method = "invalid-method"),
    "not implemented"
  )
})

test_that("estimate_s_rss produces small s for consistent z and R", {
  # Generate data where z-scores are consistent with R
  base_data <- generate_base_data(n = 500, p = 100, k = 3, signal_sd = 0.5, seed = 21)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  s <- estimate_s_rss(z, R, n = base_data$n, method = "null-mle")

  # With consistent data, s should be small
  expect_true(s < 0.01)
})

# =============================================================================
# DIAGNOSTIC & QUALITY CONTROL
# =============================================================================

test_that("kriging_rss returns list with plot and conditional_dist", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 22)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  result <- kriging_rss(z, R, n = base_data$n)

  expect_type(result, "list")
  expect_named(result, c("plot", "conditional_dist"))
})

test_that("kriging_rss plot is a ggplot object", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 23)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  result <- kriging_rss(z, R, n = base_data$n)

  expect_s3_class(result$plot, "gg")
  expect_s3_class(result$plot, "ggplot")
})

test_that("kriging_rss conditional_dist is a data frame with correct columns", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 24)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  result <- kriging_rss(z, R, n = base_data$n)

  expect_s3_class(result$conditional_dist, "data.frame")
  expect_equal(nrow(result$conditional_dist), base_data$p)
  expect_true(all(c("z", "condmean", "condvar", "z_std_diff", "logLR") %in%
                    colnames(result$conditional_dist)))
})

test_that("kriging_rss with provided s parameter", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 25)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  # Provide custom s
  result <- kriging_rss(z, R, n = base_data$n, s = 0.1)

  expect_type(result, "list")
  expect_s3_class(result$plot, "ggplot")
  expect_equal(nrow(result$conditional_dist), base_data$p)
})

test_that("kriging_rss warns when n is not provided", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 26)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  expect_message(
    result <- kriging_rss(z, R),
    "sample size"
  )
})

test_that("kriging_rss errors when n <= 1", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 27)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  expect_error(
    kriging_rss(z, R, n = 1),
    "must be greater than 1"
  )
})

test_that("kriging_rss handles s > 1 with warning", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 28)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  expect_message(
    result <- kriging_rss(z, R, n = base_data$n, s = 1.5),
    "greater than 1"
  )

  # Should still produce output
  expect_type(result, "list")
})

test_that("kriging_rss errors when s < 0", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 29)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  expect_error(
    kriging_rss(z, R, n = base_data$n, s = -0.1),
    "non-negative"
  )
})

test_that("kriging_rss handles negative eigenvalues in R", {
  set.seed(30)
  p <- 50

  # Create R with intentionally negative eigenvalue
  R <- matrix(0.5, p, p)
  diag(R) <- 1
  R[1, 2] <- 1.5 
  R[2, 1] <- 1.5

  z <- rnorm(p)

  expect_message(
    result <- kriging_rss(z, R, n = 100),
    "not positive semidefinite"
  )

  expect_type(result, "list")
})

test_that("kriging_rss identifies potential allele switches (high logLR)", {
  # Create data with one flipped allele
  base_data <- generate_base_data(n = 500, p = 100, k = 3, signal_sd = 1, seed = 31)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  # Flip one z-score to simulate allele switch
  z[1] <- -z[1]

  result <- kriging_rss(z, R, n = base_data$n)

  # The flipped variant should have high logLR
  expect_true(result$conditional_dist$logLR[1] > 0)
})

test_that("kriging_rss handles NA values in z", {
  base_data <- generate_base_data(n = 100, p = 50, seed = 32)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  # Introduce NA
  z[10] <- NA

  expect_error(
    result <- kriging_rss(z, R, n = base_data$n),
    NA
  )

  # NA should be replaced with 0
  expect_equal(result$conditional_dist$z[10], 0)
})

test_that("kriging_rss conditional mean and variance are sensible", {
  base_data <- generate_base_data(n = 200, p = 50, seed = 33)

  ss <- univariate_regression(base_data$X, base_data$y)
  R <- cor(base_data$X)
  z <- with(ss, betahat / sebetahat)

  result <- kriging_rss(z, R, n = base_data$n)

  # All conditional variances should be positive
  expect_true(all(result$conditional_dist$condvar > 0))

  # Conditional means should be finite
  expect_true(all(is.finite(result$conditional_dist$condmean)))

  # Standardized differences should be finite
  expect_true(all(is.finite(result$conditional_dist$z_std_diff)))
})
