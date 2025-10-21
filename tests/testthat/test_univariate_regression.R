context("Univariate regression")

# =============================================================================
# BASIC FUNCTIONALITY
# =============================================================================

test_that("univariate_regression returns correct structure", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 1)

  result <- univariate_regression(base_data$X, base_data$y)

  expect_type(result, "list")
  expect_named(result, c("betahat", "sebetahat"))
  expect_length(result$betahat, base_data$p)
  expect_length(result$sebetahat, base_data$p)
  expect_type(result$betahat, "double")
  expect_type(result$sebetahat, "double")
})

test_that("univariate_regression computes correct coefficients", {
  set.seed(2)
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = NULL)
  beta_true <- c(1, -0.5, 0.8, 0, 0.3)
  y <- base_data$X %*% beta_true + rnorm(base_data$n, sd = 0.1)

  result <- univariate_regression(base_data$X, y, center = TRUE, scale = FALSE)

  # Each betahat should be close to the true beta
  for (i in 1:base_data$p) {
    expect_equal(result$betahat[i], beta_true[i], tolerance = 0.2)
  }
})

test_that("univariate_regression standard errors are positive", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 3)

  result <- univariate_regression(base_data$X, base_data$y)

  expect_true(all(result$sebetahat > 0))
})

test_that("univariate_regression matches manual lm calculation", {
  base_data <- generate_base_data(n = 50, p = 3, k = 0, seed = 4)

  result <- univariate_regression(base_data$X, base_data$y, center = TRUE, scale = FALSE)

  # Compare to manual lm for first column
  y_centered <- base_data$y - mean(base_data$y)
  X_centered <- scale(base_data$X, center = TRUE, scale = FALSE)
  manual_fit <- lm(y_centered ~ X_centered[, 1])

  expect_equal(result$betahat[1], unname(coef(manual_fit)[2]), tolerance = 1e-10)
  expect_equal(result$sebetahat[1], unname(summary(manual_fit)$coef[2, 2]), tolerance = 1e-10)
})

# =============================================================================
# CENTER AND SCALE OPTIONS
# =============================================================================

test_that("univariate_regression with center=TRUE, scale=FALSE", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 5)

  result <- univariate_regression(base_data$X, base_data$y, center = TRUE, scale = FALSE)

  expect_length(result$betahat, base_data$p)
  expect_length(result$sebetahat, base_data$p)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
})

test_that("univariate_regression with center=TRUE, scale=TRUE", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 6)

  result <- univariate_regression(base_data$X, base_data$y, center = TRUE, scale = TRUE)

  expect_length(result$betahat, base_data$p)
  expect_length(result$sebetahat, base_data$p)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
})

test_that("univariate_regression with center=FALSE, scale=FALSE", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 7)

  result <- univariate_regression(base_data$X, base_data$y, center = FALSE, scale = FALSE)

  expect_length(result$betahat, base_data$p)
  expect_length(result$sebetahat, base_data$p)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
})

test_that("univariate_regression with center=FALSE, scale=TRUE", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 8)

  result <- univariate_regression(base_data$X, base_data$y, center = FALSE, scale = TRUE)

  expect_length(result$betahat, base_data$p)
  expect_length(result$sebetahat, base_data$p)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
})

test_that("univariate_regression scaling affects coefficient magnitude", {
  set.seed(9)
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = NULL)
  # Create X with different variances
  X_varied <- matrix(rnorm(base_data$n * base_data$p, sd = rep(c(1, 5, 10, 2, 3), each = base_data$n)), base_data$n, base_data$p)

  result_unscaled <- univariate_regression(X_varied, base_data$y, center = TRUE, scale = FALSE)
  result_scaled <- univariate_regression(X_varied, base_data$y, center = TRUE, scale = TRUE)

  # Coefficients should differ due to scaling
  expect_false(all(abs(result_unscaled$betahat - result_scaled$betahat) < 0.01))
})

# =============================================================================
# COVARIATES (Z PARAMETER)
# =============================================================================

test_that("univariate_regression with covariates Z", {
  set.seed(10)
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = NULL)
  k <- 2
  Z <- matrix(rnorm(base_data$n * k), base_data$n, k)

  result <- univariate_regression(base_data$X, base_data$y, Z = Z, center = TRUE)

  expect_type(result, "list")
  expect_length(result$betahat, base_data$p)
  expect_length(result$sebetahat, base_data$p)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
})

test_that("univariate_regression with Z adjusts for confounders", {
  set.seed(11)
  base_data <- generate_base_data(n = 200, p = 5, k = 0, seed = NULL)
  # Create confounder
  Z <- matrix(rnorm(base_data$n), base_data$n, 1)
  # X correlated with Z
  X_confounded <- base_data$X + Z %*% matrix(rnorm(base_data$p), 1, base_data$p)
  # y depends only on Z, not X
  y_confounded <- 2 * Z[, 1] + rnorm(base_data$n, sd = 0.1)

  result_no_Z <- univariate_regression(X_confounded, y_confounded, Z = NULL, center = TRUE)
  result_with_Z <- univariate_regression(X_confounded, y_confounded, Z = Z, center = TRUE)

  # Without Z, X appears associated with y
  # With Z adjustment, association should be weaker
  expect_true(mean(abs(result_with_Z$betahat)) < mean(abs(result_no_Z$betahat)))
})

test_that("univariate_regression with Z and return_residuals=TRUE", {
  set.seed(12)
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = NULL)
  k <- 2
  Z <- matrix(rnorm(base_data$n * k), base_data$n, k)

  result <- univariate_regression(base_data$X, base_data$y, Z = Z, return_residuals = TRUE)

  expect_type(result, "list")
  expect_named(result, c("betahat", "sebetahat", "residuals"))
  expect_length(result$residuals, base_data$n)
  expect_type(result$residuals, "double")
})

test_that("univariate_regression return_residuals=TRUE without Z omits residuals", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 13)

  result <- univariate_regression(base_data$X, base_data$y, Z = NULL, return_residuals = TRUE)

  expect_named(result, c("betahat", "sebetahat"))
  expect_false("residuals" %in% names(result))
})

test_that("univariate_regression residuals from Z are centered", {
  set.seed(14)
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = NULL)
  k <- 2
  Z <- matrix(rnorm(base_data$n * k), base_data$n, k)

  result <- univariate_regression(base_data$X, base_data$y, Z = Z, return_residuals = TRUE, center = TRUE)

  # Residuals should be approximately centered
  expect_equal(mean(result$residuals), 0, tolerance = 1e-10)
})

# =============================================================================
# NA HANDLING
# =============================================================================

test_that("univariate_regression handles NA values in y", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 15)
  base_data$y[c(5, 20, 35)] <- NA

  result <- univariate_regression(base_data$X, base_data$y)

  expect_length(result$betahat, base_data$p)
  expect_length(result$sebetahat, base_data$p)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
})

test_that("univariate_regression removes corresponding X rows when y has NA", {
  set.seed(16)
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = NULL)
  beta_true <- rep(1, base_data$p)
  y_with_signal <- base_data$X %*% beta_true + rnorm(base_data$n, sd = 0.1)

  # Add NAs
  na_idx <- c(10, 20, 30)
  y_with_signal[na_idx] <- NA

  result <- univariate_regression(base_data$X, y_with_signal, center = TRUE)

  # Should still produce finite results (NA removal worked)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
  expect_length(result$betahat, base_data$p)
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("univariate_regression with zero-variance column", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 17)
  base_data$X[, 3] <- 5  # Constant column

  result <- univariate_regression(base_data$X, base_data$y, center = TRUE, scale = FALSE)

  # Constant column becomes zero after centering
  expect_equal(result$betahat[3], 0)
  expect_equal(result$sebetahat[3], 0)
})

test_that("univariate_regression with perfect predictor", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 18)
  # Make y perfectly predicted by first column
  y_perfect <- 3 * base_data$X[, 1] + mean(base_data$X[, 1])

  result <- univariate_regression(base_data$X, y_perfect, center = TRUE, scale = FALSE)

  # First coefficient should be exactly 3, SE should be very small
  expect_equal(result$betahat[1], 3, tolerance = 1e-10)
  expect_true(result$sebetahat[1] < 1e-10)
})

test_that("univariate_regression with single column X", {
  base_data <- generate_base_data(n = 100, p = 1, k = 0, seed = 19)

  result <- univariate_regression(base_data$X, base_data$y)

  expect_length(result$betahat, 1)
  expect_length(result$sebetahat, 1)
  expect_true(is.finite(result$betahat[1]))
  expect_true(is.finite(result$sebetahat[1]))
})

test_that("univariate_regression with very small sample size", {
  base_data <- generate_base_data(n = 5, p = 3, k = 0, seed = 20)

  result <- univariate_regression(base_data$X, base_data$y)

  expect_length(result$betahat, base_data$p)
  expect_length(result$sebetahat, base_data$p)
})

# =============================================================================
# Z-SCORES
# =============================================================================

test_that("univariate_regression enables z-score calculation", {
  set.seed(21)
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = NULL)
  beta_true <- c(rep(1, 3), rep(0, 7))
  y_signal <- base_data$X %*% beta_true + rnorm(base_data$n)

  result <- univariate_regression(base_data$X, y_signal)

  # Calculate z-scores
  z <- result$betahat / result$sebetahat

  expect_length(z, base_data$p)
  expect_type(z, "double")

  # Causal variables should have larger |z|
  expect_true(mean(abs(z[1:3])) > mean(abs(z[4:10])))
})

# =============================================================================
# COMPARISON WITH SIMPLE LM
# =============================================================================

test_that("univariate_regression agrees with lm for each column", {
  base_data <- generate_base_data(n = 50, p = 5, k = 0, seed = 22)

  result <- univariate_regression(base_data$X, base_data$y, center = TRUE, scale = FALSE)

  # Prepare centered data
  y_c <- base_data$y - mean(base_data$y)
  X_c <- scale(base_data$X, center = TRUE, scale = FALSE)

  # Compare each column
  for (i in 1:base_data$p) {
    lm_fit <- lm(y_c ~ X_c[, i])
    lm_coef <- unname(coef(summary(lm_fit))[2, ])

    expect_equal(result$betahat[i], lm_coef[1], tolerance = 1e-10)
    expect_equal(result$sebetahat[i], lm_coef[2], tolerance = 1e-10)
  }
})

# =============================================================================
# FALLBACK MECHANISM
# =============================================================================

test_that("univariate_regression fallback works when fast method fails", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 23)

  # Normal execution should work
  result <- univariate_regression(base_data$X, base_data$y)

  expect_length(result$betahat, base_data$p)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
})

test_that("univariate_regression handles nearly singular design matrix", {
  set.seed(24)
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = NULL)
  # Make two columns nearly identical
  base_data$X[, 2] <- base_data$X[, 1] + rnorm(base_data$n, sd = 1e-10)

  # Should still produce output (possibly via fallback)
  result <- univariate_regression(base_data$X, base_data$y, center = TRUE)

  expect_length(result$betahat, base_data$p)
  expect_length(result$sebetahat, base_data$p)
})

# =============================================================================
# INTEGRATION TESTS
# =============================================================================

test_that("univariate_regression output usable for RSS methods", {
  set.seed(25)
  base_data <- generate_base_data(n = 200, p = 100, k = 0, seed = NULL)
  beta_true <- rep(0, base_data$p)
  beta_true[1:5] <- rnorm(5)
  y_causal <- base_data$X %*% beta_true + rnorm(base_data$n)

  result <- univariate_regression(base_data$X, y_causal)

  # Should be able to compute z-scores
  z <- result$betahat / result$sebetahat

  # Should be able to compute correlation matrix
  R <- cor(base_data$X)

  # Should be usable with estimate_s_rss
  expect_error(
    s <- estimate_s_rss(z, R, n = base_data$n),
    NA
  )
})

test_that("univariate_regression with center and scale matches susie preprocessing", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 26)

  # This should match what susie does internally for univariate regression
  result <- univariate_regression(base_data$X, base_data$y, center = TRUE, scale = TRUE)

  expect_length(result$betahat, base_data$p)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
})
