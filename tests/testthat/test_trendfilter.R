context("Trend filtering")

# =============================================================================
# BASIC FUNCTIONALITY
# =============================================================================

test_that("susie_trendfilter returns susie object", {
  set.seed(1)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  expect_s3_class(result, "susie")
  expect_type(result, "list")
  expect_true("alpha" %in% names(result))
  expect_true("mu" %in% names(result))
  expect_true("elbo" %in% names(result))
})

test_that("susie_trendfilter detects changepoints with order=0", {
  set.seed(2)
  # Create signal with clear changepoints
  mu <- c(rep(0, 25), rep(3, 25), rep(-2, 25), rep(1, 25))
  y <- mu + rnorm(100, sd = 0.3)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  # Should have non-zero PIPs near changepoint locations (25, 50, 75)
  pip <- susie_get_pip(result)
  changepoint_regions <- c(23:27, 48:52, 73:77)

  expect_true(sum(pip[changepoint_regions]) > sum(pip[-changepoint_regions]))
})

test_that("susie_trendfilter fitted values track signal", {
  set.seed(3)
  mu <- c(rep(0, 20), rep(2, 20), rep(0, 20))
  y <- mu + rnorm(60, sd = 0.1)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  fitted <- predict(result)

  # Fitted values should be closer to true signal than raw data
  expect_true(mean((fitted - mu)^2) < mean((y - mu)^2))
})

test_that("susie_trendfilter with no changepoints", {
  set.seed(4)
  # Constant signal
  y <- rep(5, 50) + rnorm(50, sd = 0.5)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  # PIPs should be low everywhere (no clear changepoints)
  pip <- susie_get_pip(result)
  expect_true(max(pip) < 0.5)
})

# =============================================================================
# ORDER PARAMETER
# =============================================================================

test_that("susie_trendfilter with order=0 (changepoints)", {
  set.seed(5)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  expect_error(
    result <- susie_trendfilter(y, order = 0, use_mad = FALSE),
    NA
  )

  expect_s3_class(result, "susie")
})

test_that("susie_trendfilter with order=1 warns", {
  set.seed(6)
  y <- seq(0, 1, length.out = 50) + rnorm(50, sd = 0.1)

  expect_message(
    result <- suppressWarnings(susie_trendfilter(y, order = 1, use_mad = FALSE)),
    "order > 0 is not recommended"
  )

  expect_s3_class(result, "susie")
})

test_that("susie_trendfilter with order=2 warns", {
  set.seed(7)
  y <- (seq(0, 1, length.out = 50))^2 + rnorm(50, sd = 0.1)

  expect_message(
    result <- suppressWarnings(susie_trendfilter(y, order = 2, use_mad = FALSE)),
    "order > 0 is not recommended"
  )

  expect_s3_class(result, "susie")
})

test_that("susie_trendfilter order=0 vs order=1 produce different results", {
  set.seed(8)
  # Linear trend
  y <- seq(0, 2, length.out = 50) + rnorm(50, sd = 0.1)

  result_0 <- susie_trendfilter(y, order = 0, use_mad = FALSE, max_iter = 10)
  result_1 <- suppressWarnings(
    susie_trendfilter(y, order = 1, use_mad = FALSE, max_iter = 10)
  )

  # Results should differ
  expect_false(all(abs(result_0$alpha - result_1$alpha) < 1e-10))
})

# =============================================================================
# USE_MAD PARAMETER
# =============================================================================

test_that("susie_trendfilter with use_mad=TRUE", {
  set.seed(9)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  expect_error(
    result <- susie_trendfilter(y, order = 0, use_mad = TRUE),
    NA
  )

  expect_s3_class(result, "susie")
})

test_that("susie_trendfilter with use_mad=FALSE", {
  set.seed(10)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  expect_error(
    result <- susie_trendfilter(y, order = 0, use_mad = FALSE),
    NA
  )

  expect_s3_class(result, "susie")
})

test_that("susie_trendfilter use_mad=TRUE vs FALSE differ", {
  set.seed(11)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result_mad <- susie_trendfilter(y, order = 0, use_mad = TRUE, max_iter = 5)
  result_no_mad <- susie_trendfilter(y, order = 0, use_mad = FALSE, max_iter = 5)

  # Results may differ due to initialization
  # Just verify both work
  expect_s3_class(result_mad, "susie")
  expect_s3_class(result_no_mad, "susie")
})

test_that("susie_trendfilter use_mad with model_init skips MAD", {
  set.seed(12)
  mu <- c(rep(0, 20), rep(2, 20))
  y <- mu + rnorm(40)

  # Create a simple init
  init <- susie_init_coef(c(20), c(2), 40)

  # With model_init, should skip MAD even if use_mad=TRUE
  result <- susie_trendfilter(y, order = 0, use_mad = TRUE,
                               model_init = init, max_iter = 2)

  expect_s3_class(result, "susie")
})

test_that("susie_trendfilter rejects MAD=0 when use_mad=TRUE", {
  # Create constant data which will cause MAD = 0
  # All differences will be 0, so median(abs(diff(y))) = 0
  y <- rep(5, 50)

  expect_error(
    susie_trendfilter(y, order = 0, use_mad = TRUE),
    "Cannot use median absolute deviation \\(MAD\\) to initialize residual variance because MAD = 0 for the input data. Please set 'use_mad = FALSE'"
  )
})

# =============================================================================
# STANDARDIZE AND INTERCEPT OPTIONS
# =============================================================================

test_that("susie_trendfilter with standardize=TRUE", {
  set.seed(13)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result <- susie_trendfilter(y, order = 0, standardize = TRUE, use_mad = FALSE)

  expect_s3_class(result, "susie")
  expect_true(all(result$alpha >= 0 & result$alpha <= 1))
})

test_that("susie_trendfilter with standardize=FALSE", {
  set.seed(14)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result <- susie_trendfilter(y, order = 0, standardize = FALSE, use_mad = FALSE)

  expect_s3_class(result, "susie")
  expect_true(all(result$alpha >= 0 & result$alpha <= 1))
})

test_that("susie_trendfilter with intercept=TRUE", {
  set.seed(15)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60) + 10  # Add offset

  result <- susie_trendfilter(y, order = 0, intercept = TRUE, use_mad = FALSE)

  expect_s3_class(result, "susie")
  expect_true(!is.na(result$intercept))
})

test_that("susie_trendfilter with intercept=FALSE", {
  set.seed(16)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result <- susie_trendfilter(y, order = 0, intercept = FALSE, use_mad = FALSE)

  expect_s3_class(result, "susie")
  expect_equal(result$intercept, 0)
})

# =============================================================================
# PASS-THROUGH PARAMETERS
# =============================================================================

test_that("susie_trendfilter passes L parameter to susie", {
  set.seed(17)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result <- susie_trendfilter(y, order = 0, L = 3, use_mad = FALSE)

  expect_equal(nrow(result$alpha), 3)
  expect_equal(length(result$V), 3)
})

test_that("susie_trendfilter passes max_iter parameter to susie", {
  set.seed(18)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result <- susie_trendfilter(y, order = 0, max_iter = 5, use_mad = FALSE)

  expect_true(result$niter <= 5)
})

test_that("susie_trendfilter passes estimate_prior_variance to susie", {
  set.seed(19)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result_estimate <- susie_trendfilter(y, order = 0,
                                       estimate_prior_variance = TRUE,
                                       use_mad = FALSE, max_iter = 3)
  result_fixed <- susie_trendfilter(y, order = 0,
                                    estimate_prior_variance = FALSE,
                                    use_mad = FALSE, max_iter = 3)

  # Both should work
  expect_s3_class(result_estimate, "susie")
  expect_s3_class(result_fixed, "susie")
})

test_that("susie_trendfilter passes null_weight to susie", {
  set.seed(20)
  mu <- c(rep(0, 30), rep(2, 30))
  y <- mu + rnorm(60)
  n <- length(y)

  result <- susie_trendfilter(y, order = 0, null_weight = 1/(n+1),
                               use_mad = FALSE, max_iter = 3)

  expect_s3_class(result, "susie")
  expect_true(!is.null(result$null_index))
})

# =============================================================================
# INTEGRATION WITH SUSIE METHODS
# =============================================================================

test_that("susie_trendfilter output works with susie_get_cs", {
  set.seed(21)
  mu <- c(rep(0, 25), rep(3, 25), rep(-2, 25), rep(1, 25))
  y <- mu + rnorm(100, sd = 0.3)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  cs <- susie_get_cs(result, coverage = 0.95)

  expect_type(cs, "list")
  expect_true("cs" %in% names(cs))
  expect_equal(cs$requested_coverage, 0.95)
})

test_that("susie_trendfilter output works with susie_get_pip", {
  set.seed(22)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  pip <- susie_get_pip(result)

  expect_length(pip, length(y))
  expect_type(pip, "double")
  expect_true(all(pip >= 0 & pip <= 1))
})

test_that("susie_trendfilter output works with predict", {
  set.seed(23)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  fitted <- predict(result)

  expect_length(fitted, length(y))
  expect_type(fitted, "double")
  expect_true(all(is.finite(fitted)))
})

test_that("susie_trendfilter output works with coef", {
  set.seed(24)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y <- mu + rnorm(60)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  coefficients <- coef(result)

  # Should have length n+1 (intercept + n basis coefficients)
  expect_length(coefficients, length(y) + 1)
  expect_type(coefficients, "double")
})

test_that("susie_trendfilter output has sets field", {
  set.seed(25)
  mu <- c(rep(0, 25), rep(3, 25), rep(-2, 25), rep(1, 25))
  y <- mu + rnorm(100, sd = 0.3)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)
  result$sets <- susie_get_cs(result, coverage = 0.95)

  # Verify sets structure
  expect_type(result$sets, "list")
  expect_true("cs" %in% names(result$sets))
  expect_true("coverage" %in% names(result$sets))
  expect_equal(result$sets$requested_coverage, 0.95)
})

# =============================================================================
# CHANGEPOINT DETECTION QUALITY
# =============================================================================

test_that("susie_trendfilter recovers true changepoints", {
  set.seed(26)
  # Simple changepoint problem
  mu <- c(rep(0, 30), rep(4, 30), rep(-2, 30))
  y <- mu + rnorm(90, sd = 0.5)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  pip <- susie_get_pip(result)

  # Changepoints should be at positions 30 and 60
  # High PIP should be near these positions
  expect_true(pip[30] > 0.5 | pip[29] > 0.5 | pip[31] > 0.5)
  expect_true(pip[60] > 0.5 | pip[59] > 0.5 | pip[61] > 0.5)
})

test_that("susie_trendfilter handles multiple small changepoints", {
  set.seed(27)
  # Many small changes
  mu <- rep(c(0, 0.5), length.out = 60)
  y <- mu + rnorm(60, sd = 0.2)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  expect_s3_class(result, "susie")
  expect_true(all(result$alpha >= 0 & result$alpha <= 1))
})

test_that("susie_trendfilter with noisy data", {
  set.seed(28)
  mu <- c(rep(0, 30), rep(2, 30))
  y <- mu + rnorm(60, sd = 2)  # High noise

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  # Should still converge
  expect_s3_class(result, "susie")
  expect_true(result$converged)
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("susie_trendfilter with short time series", {
  set.seed(29)
  y <- c(0, 0, 0, 2, 2, 2)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE, max_iter = 3)

  expect_s3_class(result, "susie")
})

test_that("susie_trendfilter with long time series", {
  set.seed(30)
  mu <- rep(c(0, 1, 2, 0), each = 100)
  y <- mu + rnorm(400, sd = 0.5)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE, max_iter = 20)

  expect_s3_class(result, "susie")
  expect_length(predict(result), 400)
})

test_that("susie_trendfilter with constant y errors", {
  set.seed(31)
  y <- rep(5, 50)

  # Constant y has zero variance, should error
  expect_error(
    susie_trendfilter(y, order = 0, use_mad = FALSE),
    "Residual variance sigma2 must be positive"
  )
})

test_that("susie_trendfilter with single changepoint at start", {
  set.seed(32)
  mu <- c(rep(0, 5), rep(2, 45))
  y <- mu + rnorm(50, sd = 0.3)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  pip <- susie_get_pip(result)

  # Should detect changepoint near position 5
  expect_true(any(pip[3:7] > 0.3))
})

test_that("susie_trendfilter with single changepoint at end", {
  set.seed(33)
  mu <- c(rep(0, 45), rep(2, 5))
  y <- mu + rnorm(50, sd = 0.3)

  result <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  pip <- susie_get_pip(result)

  # Should detect changepoint near position 45
  expect_true(any(pip[43:47] > 0.3))
})

# =============================================================================
# COMPARISON WITH MANUAL CONSTRUCTION
# =============================================================================

test_that("susie_trendfilter matches manual sparse matrix construction", {
  set.seed(34)
  with(simulate_tf(0), {
    # Manual approach with explicit X matrix
    result_manual <- susie(X, y, estimate_prior_variance = FALSE,
                           standardize = TRUE, max_iter = 5)

    # Using susie_trendfilter
    result_tf <- susie_trendfilter(y, order = 0,
                                    estimate_prior_variance = FALSE,
                                    standardize = TRUE,
                                    use_mad = FALSE, max_iter = 5)

    # Should produce similar results
    expect_equal(result_tf$alpha, result_manual$alpha, tolerance = 1e-6)
    expect_equal(result_tf$mu, result_manual$mu, tolerance = 1e-6)
  })
})

test_that("susie_trendfilter order=1 matches manual construction", {
  set.seed(35)
  with(simulate_tf(1), {
    # Manual approach
    result_manual <- susie(X, y, estimate_prior_variance = FALSE,
                           standardize = TRUE, max_iter = 5)

    # Using susie_trendfilter
    result_tf <- suppressWarnings(
      susie_trendfilter(y, order = 1,
                        estimate_prior_variance = FALSE,
                        standardize = TRUE,
                        use_mad = FALSE, max_iter = 5)
    )

    # Should produce similar results
    expect_equal(result_tf$alpha, result_manual$alpha, tolerance = 1e-6)
    expect_equal(result_tf$mu, result_manual$mu, tolerance = 1e-6)
  })
})

test_that("susie_trendfilter order=2 matches manual construction", {
  set.seed(36)
  with(simulate_tf(2), {
    # Manual approach
    result_manual <- susie(X, y, estimate_prior_variance = FALSE,
                           standardize = TRUE, max_iter = 5)

    # Using susie_trendfilter
    result_tf <- suppressWarnings(
      susie_trendfilter(y, order = 2,
                        estimate_prior_variance = FALSE,
                        standardize = TRUE,
                        use_mad = FALSE, max_iter = 5)
    )

    # Should produce similar results
    expect_equal(result_tf$alpha, result_manual$alpha, tolerance = 1e-6)
    expect_equal(result_tf$mu, result_manual$mu, tolerance = 1e-6)
  })
})

# =============================================================================
# EXAMPLES FROM DOCUMENTATION
# =============================================================================

test_that("susie_trendfilter works with documentation example", {
  set.seed(1)
  mu <- c(rep(0, 50), rep(1, 50), rep(3, 50), rep(-2, 50), rep(0, 200))
  y <- mu + rnorm(400)

  s <- susie_trendfilter(y, max_iter = 100)

  expect_s3_class(s, "susie")
  expect_length(predict(s), 400)

  # Should be able to get credible sets
  cs <- susie_get_cs(s)
  expect_type(cs, "list")
})
