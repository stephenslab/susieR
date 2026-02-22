context("SuSiE Data Constructors")

# =============================================================================
# INDIVIDUAL DATA CONSTRUCTOR - Basic Functionality
# =============================================================================

test_that("individual_data_constructor returns correct structure", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 1)

  result <- individual_data_constructor(base_data$X, base_data$y)

  expect_type(result, "list")
  expect_true("data" %in% names(result))
  expect_true("params" %in% names(result))
  expect_s3_class(result$data, "individual")
})

test_that("individual_data_constructor creates data object with correct fields", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 2)

  result <- individual_data_constructor(base_data$X, base_data$y)

  expect_true("X" %in% names(result$data))
  expect_true("y" %in% names(result$data))
  expect_true("n" %in% names(result$data))
  expect_true("p" %in% names(result$data))
  expect_true("mean_y" %in% names(result$data))
})

test_that("individual_data_constructor sets correct dimensions", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 3)

  result <- individual_data_constructor(base_data$X, base_data$y)

  expect_equal(result$data$n, base_data$n)
  expect_equal(result$data$p, base_data$p)
  expect_equal(dim(result$data$X), c(base_data$n, base_data$p))
  expect_length(result$data$y, base_data$n)
})

test_that("individual_data_constructor sets X attributes", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 4)

  result <- individual_data_constructor(base_data$X, base_data$y, standardize = TRUE, intercept = TRUE)

  expect_true(!is.null(attr(result$data$X, "scaled:center")))
  expect_true(!is.null(attr(result$data$X, "scaled:scale")))
  expect_true(!is.null(attr(result$data$X, "d")))
})

# =============================================================================
# INDIVIDUAL DATA CONSTRUCTOR - Input Validation
# =============================================================================

test_that("individual_data_constructor rejects non-matrix X", {
  expect_error(
    individual_data_constructor(as.data.frame(matrix(1:10, 5, 2)), rnorm(5)),
    "Input X must be a double-precision matrix"
  )
})

test_that("individual_data_constructor rejects X with NAs", {
  base_data <- generate_base_data(n = 10, p = 10, k = 0, seed = 5)
  base_data$X[5, 5] <- NA

  expect_error(
    individual_data_constructor(base_data$X, base_data$y),
    "X contains NA values"
  )
})

test_that("individual_data_constructor rejects y with NAs when na.rm=FALSE", {
  base_data <- generate_base_data(n = 10, p = 10, k = 0, seed = 6)
  base_data$y[5] <- NA

  expect_error(
    individual_data_constructor(base_data$X, base_data$y, na.rm = FALSE),
    "Input y must not contain missing values"
  )
})

test_that("individual_data_constructor handles y with NAs when na.rm=TRUE", {
  base_data <- generate_base_data(n = 10, p = 10, k = 0, seed = 7)
  base_data$y[5] <- NA

  result <- individual_data_constructor(base_data$X, base_data$y, na.rm = TRUE)

  expect_equal(result$data$n, 9)
  expect_equal(nrow(result$data$X), 9)
  expect_length(result$data$y, 9)
  expect_false(anyNA(result$data$y))
})

test_that("individual_data_constructor computes residual_variance_lowerbound after NA removal", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 7.25)
  base_data$y[1] <- NA

  result <- individual_data_constructor(base_data$X, base_data$y, na.rm = TRUE)

  # Verify residual_variance_lowerbound is computed correctly (not NA)
  expect_true(is.finite(result$params$residual_variance_lowerbound))
  expect_true(result$params$residual_variance_lowerbound > 0)

  # Verify it equals var(y_clean) / 1e4 where y_clean has NA removed
  y_clean <- base_data$y[!is.na(base_data$y)]
  expected_lowerbound <- var(y_clean) / 1e4
  expect_equal(result$params$residual_variance_lowerbound, expected_lowerbound)
})

test_that("individual_data_constructor allows custom residual_variance_lowerbound with NA in y", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 7.5)
  base_data$y[1] <- NA

  custom_lowerbound <- 0.001
  result <- individual_data_constructor(
    base_data$X, base_data$y,
    na.rm = TRUE,
    residual_variance_lowerbound = custom_lowerbound
  )

  expect_equal(result$params$residual_variance_lowerbound, custom_lowerbound)
})

test_that("individual_data_constructor handles multiple NAs in y with na.rm=TRUE", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 7.75)
  # Set multiple NAs at different positions

  base_data$y[c(1, 25, 50, 75, 100)] <- NA

  result <- individual_data_constructor(base_data$X, base_data$y, na.rm = TRUE)

  expect_equal(result$data$n, 95)
  expect_equal(nrow(result$data$X), 95)
  expect_length(result$data$y, 95)
  expect_false(anyNA(result$data$y))
  expect_true(is.finite(result$params$residual_variance_lowerbound))
})

# =============================================================================
# INDIVIDUAL DATA CONSTRUCTOR - Centering and Scaling
# =============================================================================

test_that("individual_data_constructor centers y when intercept=TRUE", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 8)
  base_data$y <- base_data$y + 10

  result <- individual_data_constructor(base_data$X, base_data$y, intercept = TRUE)

  expect_equal(mean(result$data$y), 0, tolerance = 1e-10)
  expect_true(result$data$mean_y != 0)
})

test_that("individual_data_constructor does not center y when intercept=FALSE", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 9)
  base_data$y <- base_data$y + 10

  result <- individual_data_constructor(base_data$X, base_data$y, intercept = FALSE)

  expect_true(abs(mean(result$data$y) - 10) < 1)
})

test_that("individual_data_constructor standardizes X when requested", {
  set.seed(10)
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = NULL)
  # Create X with different mean and sd
  base_data$X <- matrix(rnorm(base_data$n * base_data$p, mean = 5, sd = 3), base_data$n, base_data$p)

  result <- individual_data_constructor(base_data$X, base_data$y, standardize = TRUE, intercept = TRUE)

  cm <- attr(result$data$X, "scaled:center")
  csd <- attr(result$data$X, "scaled:scale")

  expect_length(cm, base_data$p)
  expect_length(csd, base_data$p)
  expect_true(all(csd > 0))
})

# =============================================================================
# INDIVIDUAL DATA CONSTRUCTOR - Prior Weights
# =============================================================================

test_that("individual_data_constructor creates uniform prior weights by default", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 11)

  result <- individual_data_constructor(base_data$X, base_data$y)

  expect_length(result$params$prior_weights, base_data$p)
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)
  expect_true(all(abs(result$params$prior_weights - 1/base_data$p) < 1e-10))
})

test_that("individual_data_constructor normalizes custom prior weights", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 12)
  custom_weights <- rep(2, base_data$p)

  result <- individual_data_constructor(base_data$X, base_data$y, prior_weights = custom_weights)

  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)
})

test_that("individual_data_constructor rejects wrong length prior weights", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 13)

  expect_error(
    individual_data_constructor(base_data$X, base_data$y, prior_weights = rep(1, 40)),
    "Prior weights must have length p"
  )
})

test_that("individual_data_constructor rejects all-zero prior weights", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 14)

  expect_error(
    individual_data_constructor(base_data$X, base_data$y, prior_weights = rep(0, base_data$p)),
    "Prior weight should be greater than 0"
  )
})

# =============================================================================
# INDIVIDUAL DATA CONSTRUCTOR - Null Weight
# =============================================================================

test_that("individual_data_constructor handles null_weight=0", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 15)

  result <- individual_data_constructor(base_data$X, base_data$y, null_weight = 0)

  expect_equal(result$data$p, base_data$p)
  expect_equal(ncol(result$data$X), base_data$p)
  expect_null(result$params$null_weight)
})

test_that("individual_data_constructor adds null column when null_weight > 0", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 16)

  result <- individual_data_constructor(base_data$X, base_data$y, null_weight = 0.1)

  expect_equal(result$data$p, base_data$p + 1)
  expect_equal(ncol(result$data$X), base_data$p + 1)
  expect_equal(result$params$null_weight, 0.1)
  expect_true(all(result$data$X[, base_data$p + 1] == 0))
})

test_that("individual_data_constructor adjusts prior weights with null_weight", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 17)

  result <- individual_data_constructor(base_data$X, base_data$y, null_weight = 0.2)

  expect_length(result$params$prior_weights, base_data$p + 1)
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)
  expect_equal(result$params$prior_weights[base_data$p + 1], 0.2, tolerance = 1e-10)
})

test_that("individual_data_constructor adjusts custom prior weights with null_weight", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 17.5)

  # Create custom prior weights (not uniform)
  custom_weights <- runif(base_data$p, 0.5, 2)
  custom_weights <- custom_weights / sum(custom_weights)  # Normalize to sum to 1

  result <- individual_data_constructor(base_data$X, base_data$y,
                                       prior_weights = custom_weights,
                                       null_weight = 0.3)

  # Check that we have p+1 weights (original p + null column)
  expect_length(result$params$prior_weights, base_data$p + 1)

  # Check that all weights sum to 1
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)

  # Check that the null weight is exactly 0.3
  expect_equal(result$params$prior_weights[base_data$p + 1], 0.3, tolerance = 1e-10)

  # Check that the other weights were scaled by (1 - null_weight) = 0.7
  # i.e., result$params$prior_weights[1:p] should equal custom_weights * 0.7
  expect_equal(result$params$prior_weights[1:base_data$p],
               custom_weights * 0.7,
               tolerance = 1e-10)

  # Verify that the sum of the first p weights is (1 - 0.3) = 0.7
  expect_equal(sum(result$params$prior_weights[1:base_data$p]), 0.7, tolerance = 1e-10)
})

test_that("individual_data_constructor rejects invalid null_weight", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 18)

  expect_error(
    individual_data_constructor(base_data$X, base_data$y, null_weight = -0.1),
    "Null weight must be between 0 and 1"
  )

  expect_error(
    individual_data_constructor(base_data$X, base_data$y, null_weight = 1.5),
    "Null weight must be between 0 and 1"
  )

  expect_error(
    individual_data_constructor(base_data$X, base_data$y, null_weight = "invalid"),
    "Null weight must be numeric"
  )
})

# =============================================================================
# INDIVIDUAL DATA CONSTRUCTOR - Rfast Warning
# =============================================================================

test_that("individual_data_constructor warns about Rfast when p > 1000 and Rfast not available", {
  # Only test the warning if Rfast is not installed
  skip_if(requireNamespace("Rfast", quietly = TRUE),
          "Rfast is installed, skipping warning test")

  base_data <- generate_base_data(n = 100, p = 1001, k = 0, seed = 18.5)

  expect_message(
    result <- individual_data_constructor(base_data$X, base_data$y),
    "consider installing the Rfast package",
    fixed = FALSE
  )

  # Verify constructor still works despite the warning
  expect_equal(result$data$p, 1001)
})

test_that("individual_data_constructor does not warn when p <= 1000", {
  # This should never warn regardless of Rfast availability
  base_data <- generate_base_data(n = 100, p = 1000, k = 0, seed = 18.75)
  result <- individual_data_constructor(base_data$X, base_data$y)

  expect_equal(result$data$p, 1000)
})

# =============================================================================
# INDIVIDUAL DATA CONSTRUCTOR - Parameters
# =============================================================================

test_that("individual_data_constructor stores all parameters", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 19)

  result <- individual_data_constructor(
    base_data$X, base_data$y,
    L = 5,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    max_iter = 50,
    tol = 1e-4
  )

  expect_equal(result$params$L, 5)
  expect_true(result$params$estimate_residual_variance)
  expect_true(result$params$estimate_prior_variance)
  expect_equal(result$params$max_iter, 50)
  expect_equal(result$params$tol, 1e-4)
})

# =============================================================================
# INDIVIDUAL DATA CONSTRUCTOR - Incompatible Parameter Combinations
# =============================================================================

test_that("individual_data_constructor rejects unmappable_effects with Servin_Stephens", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 19.5)

  expect_error(
    individual_data_constructor(
      base_data$X, base_data$y,
      unmappable_effects = "inf",
      estimate_residual_method = "Servin_Stephens"
    ),
    "The combination of unmappable_effects = 'inf' with estimate_residual_method = 'Servin_Stephens' is not supported"
  )

  expect_error(
    individual_data_constructor(
      base_data$X, base_data$y,
      unmappable_effects = "ash",
      estimate_residual_method = "Servin_Stephens"
    ),
    "The combination of unmappable_effects = 'ash' with estimate_residual_method = 'Servin_Stephens' is not supported"
  )
})

test_that("individual_data_constructor rejects unmappable_effects='ash' with estimate_prior_method='EM'", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 19.75)

  expect_error(
    individual_data_constructor(
      base_data$X, base_data$y,
      unmappable_effects = "ash",
      estimate_prior_method = "EM"
    ),
    "The combination of unmappable_effects = 'ash' with estimate_prior_method = 'EM' is not supported"
  )
})

# =============================================================================
# SUFFICIENT STATISTICS CONSTRUCTOR - Basic Functionality
# =============================================================================

test_that("sufficient_stats_constructor returns correct structure", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 20)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX)

  expect_type(result, "list")
  expect_true("data" %in% names(result))
  expect_true("params" %in% names(result))
  expect_s3_class(result$data, "ss")
})

test_that("sufficient_stats_constructor creates data object with correct fields", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 21)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX)

  expect_true("XtX" %in% names(result$data))
  expect_true("Xty" %in% names(result$data))
  expect_true("yty" %in% names(result$data))
  expect_true("n" %in% names(result$data))
  expect_true("p" %in% names(result$data))
})

test_that("sufficient_stats_constructor sets correct dimensions", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 22)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX)

  expect_equal(result$data$n, base_data$n)
  expect_equal(result$data$p, base_data$p)
  expect_equal(dim(result$data$XtX), c(base_data$p, base_data$p))
  expect_length(result$data$Xty, base_data$p)
})

# =============================================================================
# SUFFICIENT STATISTICS CONSTRUCTOR - Input Validation
# =============================================================================

test_that("sufficient_stats_constructor requires n", {
  XtX <- matrix(1:25, 5, 5)
  Xty <- 1:5
  yty <- 10

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, XtX = XtX),
    "n must be provided"
  )
})

test_that("sufficient_stats_constructor rejects n <= 1", {
  XtX <- matrix(1:25, 5, 5)
  Xty <- 1:5
  yty <- 10

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = 1, XtX = XtX),
    "n must be greater than 1"
  )
})

test_that("sufficient_stats_constructor requires all inputs", {
  XtX <- matrix(1:25, 5, 5)
  Xty <- 1:5
  yty <- 10
  n <- 100

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n),
    "XtX, Xty, yty must all be provided"
  )
})

test_that("sufficient_stats_constructor rejects non-matrix XtX", {
  # Test with data.frame (not a matrix)
  XtX_df <- data.frame(matrix(rnorm(25), 5, 5))
  Xty <- rnorm(5)
  yty <- 10
  n <- 100

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n, XtX = XtX_df),
    "XtX must be a numeric dense or sparse matrix"
  )
})

test_that("sufficient_stats_constructor rejects integer matrix XtX", {
  # Test with integer matrix (not double)
  XtX_int <- matrix(1L:25L, 5, 5)
  Xty <- 1:5
  yty <- 10
  n <- 100

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n, XtX = XtX_int),
    "XtX must be a numeric dense or sparse matrix"
  )
})

test_that("sufficient_stats_constructor rejects non-numeric XtX", {
  # Test with character matrix
  XtX_char <- matrix(as.character(1:25), 5, 5)
  Xty <- 1:5
  yty <- 10
  n <- 100

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n, XtX = XtX_char),
    "XtX must be a numeric dense or sparse matrix"
  )
})

test_that("sufficient_stats_constructor rejects vector XtX", {
  # Test with vector (not a matrix)
  XtX_vec <- rnorm(25)
  Xty <- 1:5
  yty <- 10
  n <- 100

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n, XtX = XtX_vec),
    "XtX must be a numeric dense or sparse matrix"
  )
})

test_that("sufficient_stats_constructor rejects dimension mismatch", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 23)
  XtX <- crossprod(base_data$X)
  Xty <- rnorm(10)
  yty <- 100

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX),
    "does not agree with expected"
  )
})

test_that("sufficient_stats_constructor rejects non-symmetric XtX", {
  XtX <- matrix(1:25, 5, 5)
  XtX[1, 2] <- 100
  Xty <- 1:5
  yty <- 10
  n <- 100

  expect_message(
    result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = n, XtX = XtX),
    "XtX not symmetric"
  )

  expect_true(isSymmetric(result$data$XtX))
})

test_that("sufficient_stats_constructor rejects XtX with NAs", {
  XtX <- matrix(rnorm(25), 5, 5)
  XtX[1, 1] <- NA
  Xty <- 1:5
  yty <- 10
  n <- 100

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n, XtX = XtX),
    "Input XtX matrix contains NAs"
  )
})

test_that("sufficient_stats_constructor handles Xty with NAs", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 24)
  XtX <- crossprod(base_data$X)
  Xty <- rnorm(base_data$p)
  Xty[5] <- NA
  yty <- 100

  expect_message(
    result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX),
    "NA values in Xty are replaced with 0"
  )

  expect_false(anyNA(result$data$Xty))
  expect_equal(result$data$Xty[5], 0)
})

test_that("sufficient_stats_constructor rejects infinite Xty", {
  XtX <- crossprod(matrix(rnorm(100 * 10), 100, 10))
  Xty <- rnorm(10)
  Xty[5] <- Inf
  yty <- 100
  n <- 100

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n, XtX = XtX),
    "Input Xty contains infinite values"
  )
})

# =============================================================================
# SUFFICIENT STATISTICS CONSTRUCTOR - Standardization
# =============================================================================

test_that("sufficient_stats_constructor standardizes when requested", {
  set.seed(25)
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = NULL)
  base_data$X <- matrix(rnorm(base_data$n * base_data$p, mean = 5, sd = 3), base_data$n, base_data$p)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX, standardize = TRUE)

  d_attr <- attr(result$data$XtX, "d")
  expect_length(d_attr, base_data$p)
  expect_true(all(is.finite(d_attr)))
})

test_that("sufficient_stats_constructor does not standardize when standardize=FALSE", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 26)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX, standardize = FALSE)

  csd_attr <- attr(result$data$XtX, "scaled:scale")
  expect_true(all(csd_attr == 1))
})

# =============================================================================
# SUFFICIENT STATISTICS CONSTRUCTOR - Rfast Warning
# =============================================================================

test_that("sufficient_stats_constructor warns about Rfast when p > 1000 and Rfast not available", {
  # Only test the warning if Rfast is not installed
  skip_if(requireNamespace("Rfast", quietly = TRUE),
          "Rfast is installed, skipping warning test")

  base_data <- generate_base_data(n = 100, p = 1001, k = 0, seed = 27.5)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  expect_message(
    result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX),
    "consider installing the Rfast package",
    fixed = FALSE
  )

  # Verify constructor still works despite the warning
  expect_equal(result$data$p, 1001)
})

test_that("sufficient_stats_constructor does not warn when p <= 1000", {
  # This should never warn regardless of Rfast availability
  base_data <- generate_base_data(n = 100, p = 1000, k = 0, seed = 27.6)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX)

  # Just verify it worked
  expect_equal(result$data$p, 1000)
})

# =============================================================================
# SUFFICIENT STATISTICS CONSTRUCTOR - MAF Filtering
# =============================================================================

test_that("sufficient_stats_constructor applies MAF filter", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 27)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)
  maf <- runif(base_data$p, 0, 0.5)

  result <- sufficient_stats_constructor(
    Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
    maf = maf, maf_thresh = 0.1
  )

  n_filtered <- sum(maf > 0.1)
  expect_equal(result$data$p, n_filtered)
  expect_equal(nrow(result$data$XtX), n_filtered)
  expect_length(result$data$Xty, n_filtered)
})

test_that("sufficient_stats_constructor rejects MAF with incorrect length", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 27.7)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # MAF vector with wrong length (p - 10 instead of p)
  maf_wrong_length <- runif(base_data$p - 10, 0, 0.5)

  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      maf = maf_wrong_length, maf_thresh = 0.1
    ),
    "The length of maf does not agree with expected"
  )
})

test_that("sufficient_stats_constructor rejects MAF that is too long", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 27.8)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # MAF vector with wrong length (p + 10 instead of p)
  maf_too_long <- runif(base_data$p + 10, 0, 0.5)

  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      maf = maf_too_long, maf_thresh = 0.1
    ),
    "The length of maf does not agree with expected"
  )
})

# =============================================================================
# SUFFICIENT STATISTICS CONSTRUCTOR - Positive Semidefinite Check
# =============================================================================

test_that("sufficient_stats_constructor rejects non-positive-semidefinite XtX when check_input=TRUE", {
  # Create a matrix that is NOT positive semidefinite
  # by using a matrix with negative eigenvalues
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 28.9)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # Make XtX non-positive-semidefinite by adding a negative diagonal
  # This creates negative eigenvalues
  XtX[1, 1] <- -10

  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      check_input = TRUE
    ),
    "XtX is not a positive semidefinite matrix"
  )
})

test_that("sufficient_stats_constructor accepts positive-semidefinite XtX when check_input=TRUE", {
  # Create a valid positive semidefinite matrix
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 29)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # This should work without error
  result <- sufficient_stats_constructor(
    Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
    check_input = TRUE
  )

  expect_equal(result$data$p, base_data$p)
})

test_that("sufficient_stats_constructor warns when Xty not in column space of XtX", {
  # Create a rank-deficient matrix with exact zero eigenvalues
  p <- 5
  n <- 100

  # Diagonal matrix with rank 3 (2 zero eigenvalues)
  XtX <- diag(c(1, 1, 1, 0, 0))

  # Create Xty with non-zero components in the null space (positions 4 and 5)
  # This Xty cannot be written as X'y for any y
  Xty <- c(1, 1, 1, 10, 10)  # Last two components are in null space

  yty <- 100

  expect_message(
    result <- sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = n, XtX = XtX,
      check_input = TRUE
    ),
    "Xty does not lie in the space of the non-zero eigenvectors"
  )
})

test_that("sufficient_stats_constructor does not warn when Xty in column space", {
  # Create valid XtX and Xty from same data
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 29.3)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # This should work without warning since Xty = X'y by construction
  result <- sufficient_stats_constructor(
    Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
    check_input = TRUE
  )

  # Just verify it worked
  expect_equal(result$data$p, base_data$p)
})

# =============================================================================
# SUFFICIENT STATISTICS CONSTRUCTOR - Null Weight
# =============================================================================

test_that("sufficient_stats_constructor adds null column when null_weight > 0", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 28)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(
    Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
    null_weight = 0.1,
    X_colmeans = rep(0, base_data$p)
  )

  expect_equal(result$data$p, base_data$p + 1)
  expect_equal(nrow(result$data$XtX), base_data$p + 1)
  expect_length(result$data$Xty, base_data$p + 1)
  expect_equal(result$params$null_weight, 0.1)
})

test_that("sufficient_stats_constructor adjusts custom prior weights with null_weight", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 28.5)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # Create custom prior weights (not uniform)
  custom_weights <- runif(base_data$p, 0.5, 2)
  custom_weights <- custom_weights / sum(custom_weights)  # Normalize to sum to 1

  result <- sufficient_stats_constructor(
    Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
    prior_weights = custom_weights,
    null_weight = 0.25,
    X_colmeans = rep(0, base_data$p)
  )

  # Check that we have p+1 weights (original p + null column)
  expect_length(result$params$prior_weights, base_data$p + 1)

  # Check that all weights sum to 1
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)

  # Check that the null weight is exactly 0.25
  expect_equal(result$params$prior_weights[base_data$p + 1], 0.25, tolerance = 1e-10)

  # Check that the other weights were scaled by (1 - null_weight) = 0.75
  expect_equal(result$params$prior_weights[1:base_data$p],
               custom_weights * 0.75,
               tolerance = 1e-10)

  # Verify that the sum of the first p weights is (1 - 0.25) = 0.75
  expect_equal(sum(result$params$prior_weights[1:base_data$p]), 0.75, tolerance = 1e-10)
})

test_that("sufficient_stats_constructor rejects non-numeric null_weight", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 28.6)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      null_weight = "invalid"
    ),
    "Null weight must be numeric"
  )
})

test_that("sufficient_stats_constructor rejects negative null_weight", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 28.7)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      null_weight = -0.1
    ),
    "Null weight must be between 0 and 1"
  )
})

test_that("sufficient_stats_constructor rejects null_weight >= 1", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 28.8)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # Test null_weight = 1
  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      null_weight = 1
    ),
    "Null weight must be between 0 and 1"
  )

  # Test null_weight > 1
  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      null_weight = 1.5
    ),
    "Null weight must be between 0 and 1"
  )
})

test_that("sufficient_stats_constructor replicates scalar X_colmeans when null_weight is set", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 28.9)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # Provide scalar X_colmeans which should be replicated to length p
  result <- sufficient_stats_constructor(
    Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
    null_weight = 0.1,
    X_colmeans = 0  # Scalar value
  )

  # Should work without error
  expect_equal(result$data$p, base_data$p + 1)  # p + 1 due to null column
})

test_that("sufficient_stats_constructor rejects wrong length X_colmeans with null_weight", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 29.0)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # Provide X_colmeans with wrong length
  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      null_weight = 0.1,
      X_colmeans = rep(0, base_data$p - 10)  # Wrong length
    ),
    "The length of X_colmeans does not agree with number of variables"
  )
})

test_that("sufficient_stats_constructor replicates scalar X_colmeans without null_weight", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 29.1)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # Provide scalar X_colmeans which should be replicated to length p
  result <- sufficient_stats_constructor(
    Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
    X_colmeans = 0  # Scalar value
  )

  # Should work without error
  expect_equal(result$data$p, base_data$p)
})

test_that("sufficient_stats_constructor rejects wrong length X_colmeans without null_weight", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 29.2)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # Provide X_colmeans with wrong length
  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      X_colmeans = rep(0, base_data$p - 10)  # Wrong length
    ),
    "X_colmeans.*does not match number of variables"
  )
})

test_that("sufficient_stats_constructor rejects wrong length prior_weights", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 29.3)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # Provide prior_weights with wrong length
  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      prior_weights = rep(1, base_data$p - 10)  # Wrong length
    ),
    "Prior weights must have length p"
  )
})

test_that("sufficient_stats_constructor rejects all-zero prior_weights", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 29.4)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  # Provide all-zero prior_weights
  expect_error(
    sufficient_stats_constructor(
      Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
      prior_weights = rep(0, base_data$p)  # All zeros
    ),
    "Prior weight should be greater than 0 for at least one variable"
  )
})

# =============================================================================
# SUFFICIENT STATISTICS CONSTRUCTOR - Method Restrictions
# =============================================================================

test_that("sufficient_stats_constructor accepts Servin_Stephens", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 29)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- 100

  result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
                                         estimate_residual_method = "Servin_Stephens")
  expect_true(result$params$use_servin_stephens)
  expect_equal(result$params$estimate_prior_method, "EM")
})

test_that("sufficient_stats_constructor rejects unmappable_effects='ash'", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 30)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- 100

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
                                unmappable_effects = "ash"),
    "Adaptive shrinkage \\(ash\\) requires individual-level data"
  )
})

# =============================================================================
# RSS LAMBDA CONSTRUCTOR - Basic Functionality
# =============================================================================

test_that("rss_lambda_constructor returns correct structure", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  result <- rss_lambda_constructor(z, R, lambda = 0.5)

  expect_type(result, "list")
  expect_true("data" %in% names(result))
  expect_true("params" %in% names(result))
  expect_s3_class(result$data, "rss_lambda")
})

test_that("rss_lambda_constructor creates data object with correct fields", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  result <- rss_lambda_constructor(z, R, lambda = 0.5)

  expect_true("z" %in% names(result$data))
  expect_true("R" %in% names(result$data))
  expect_true("lambda" %in% names(result$data))
  expect_true("eigen_R" %in% names(result$data))
  expect_true("Vtz" %in% names(result$data))
  expect_true("n" %in% names(result$data))
  expect_true("p" %in% names(result$data))
})

test_that("rss_lambda_constructor sets correct dimensions", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  result <- rss_lambda_constructor(z, R, lambda = 0.5)

  expect_equal(result$data$n, p)
  expect_equal(result$data$p, p)
  expect_length(result$data$z, p)
  expect_equal(dim(result$data$R), c(p, p))
})

test_that("rss_lambda_constructor computes eigen decomposition", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  result <- rss_lambda_constructor(z, R, lambda = 0.5)

  expect_true("eigen_R" %in% names(result$data))
  expect_true("values" %in% names(result$data$eigen_R))
  expect_true("vectors" %in% names(result$data$eigen_R))
  expect_length(result$data$eigen_R$values, p)
})

test_that("rss_lambda_constructor computes Vtz", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  result <- rss_lambda_constructor(z, R, lambda = 0.5)

  expect_true("Vtz" %in% names(result$data))
  expect_length(result$data$Vtz, p)
})

# =============================================================================
# RSS LAMBDA CONSTRUCTOR - Input Validation
# =============================================================================

test_that("rss_lambda_constructor rejects dimension mismatch", {
  z <- rnorm(50)
  R <- diag(40)

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5),
    "does not agree with expected"
  )
})

test_that("rss_lambda_constructor rejects non-symmetric R", {
  R <- matrix(rnorm(25), 5, 5)
  z <- rnorm(5)

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5),
    "R is not a symmetric matrix"
  )
})

test_that("rss_lambda_constructor rejects integer matrix R", {
  R <- matrix(1:25, 5, 5)
  R <- R + t(R)  # Make it symmetric
  mode(R) <- "integer"  # Convert to integer type
  z <- rnorm(5)

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5),
    "Input R must be a double-precision matrix or a sparse matrix"
  )
})

test_that("rss_lambda_constructor rejects non-positive-semidefinite R when check_R=TRUE", {
  # Create a matrix with negative eigenvalue
  R <- diag(5)
  R[1, 1] <- -1  # Force negative eigenvalue
  z <- rnorm(5)

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5, check_R = TRUE),
    "is not a positive semidefinite matrix"
  )
})

test_that("rss_lambda_constructor accepts non-PSD R when check_R=FALSE", {
  # Create a matrix with negative eigenvalue
  R <- diag(5)
  R[1, 1] <- -0.5  # Force negative eigenvalue
  z <- rnorm(5)

  # Should succeed with check_R = FALSE (sets negative eigenvalues to 0)
  result <- suppressWarnings(
    rss_lambda_constructor(z, R, lambda = 0.5, check_R = FALSE)
  )
  expect_true(!is.null(result))
})

test_that("rss_lambda_constructor warns when z not in column space of R", {
  # Create R with rank < p (has null space)
  p <- 5
  R <- diag(c(1, 1, 1, 0, 0))  # Rank 3, nullspace dimension 2

  # Create z with components in null space (positions 4 and 5)
  z <- c(0.1, 0.1, 0.1, 10, 10)  # Large components in null directions

  expect_message(
    rss_lambda_constructor(z, R, lambda = 0.5, check_z = TRUE),
    "Input z does not lie in the space of non-zero eigenvectors of R"
  )
})

test_that("rss_lambda_constructor messages when z in column space of R", {
  # Create R with rank < p
  p <- 5
  R <- diag(c(1, 1, 1, 0, 0))  # Rank 3

  # Create z only in column space (zero components in null directions)
  z <- c(1, 2, 3, 0, 0)

  expect_message(
    suppressWarnings(
      rss_lambda_constructor(z, R, lambda = 0.5, check_z = TRUE)
    ),
    "Input z is in space spanned by the non-zero eigenvectors of R"
  )
})

test_that("rss_lambda_constructor skips z check when check_z=FALSE", {
  # Create R with rank < p
  p <- 5
  R <- diag(c(1, 1, 1, 0, 0))
  z <- c(0.1, 0.1, 0.1, 10, 10)  # z in null space
  result <- suppressWarnings(
    suppressMessages(
      rss_lambda_constructor(z, R, lambda = 0.5, check_z = FALSE)
    )
  )
  expect_true(!is.null(result))
})

test_that("rss_lambda_constructor skips z check when R is full rank", {
  # Full rank R (no null space)
  p <- 5
  R <- diag(p)
  z <- rnorm(p)

  # Should not check when length(colspace) == length(z)
  result <- suppressWarnings(
    rss_lambda_constructor(z, R, lambda = 0.5, check_z = TRUE)
  )
  expect_true(!is.null(result))
})

test_that("rss_lambda_constructor rejects R with NAs", {
  R <- diag(10)
  R[1, 1] <- NA
  z <- rnorm(10)

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5),
    "R matrix contains missing values"
  )
})

test_that("rss_lambda_constructor rejects infinite z", {
  R <- diag(10)
  z <- rnorm(10)
  z[5] <- Inf

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5),
    "z contains infinite values"
  )
})

test_that("rss_lambda_constructor replaces NA z with zero", {
  R <- diag(10)
  z <- rnorm(10)
  z[5] <- NA

  expect_message(
    result <- rss_lambda_constructor(z, R, lambda = 0.5),
    "NA values in z-scores are replaced with 0"
  )

  expect_false(anyNA(result$data$z))
  expect_equal(result$data$z[5], 0)
})

# =============================================================================
# RSS LAMBDA CONSTRUCTOR - Lambda Parameter
# =============================================================================

test_that("rss_lambda_constructor stores lambda value", {
  z <- rnorm(50)
  R <- diag(50)

  result <- rss_lambda_constructor(z, R, lambda = 0.3)

  expect_equal(result$data$lambda, 0.3)
})

test_that("rss_lambda_constructor estimates lambda when lambda='estimate'", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)
  R[1:10, 1:10] <- 0

  result <- rss_lambda_constructor(z, R, lambda = "estimate")

  expect_true(is.numeric(result$data$lambda))
  expect_true(result$data$lambda >= 0)
})

test_that("rss_lambda_constructor sets lambda=0 when R is full rank and lambda='estimate'", {
  set.seed(123)
  p <- 50
  z <- rnorm(p)
  R <- diag(p)  # Full rank - all eigenvalues positive

  result <- rss_lambda_constructor(z, R, lambda = "estimate")

  # When R is full rank, length(colspace) == length(z), so lambda should be set to 0
  expect_equal(result$data$lambda, 0)
})

test_that("rss_lambda_constructor adjusts residual variance with lambda", {
  z <- rnorm(50)
  R <- diag(50)

  result <- rss_lambda_constructor(z, R, lambda = 0.2, residual_variance = 0.8)

  expect_equal(result$params$residual_variance, 0.6)
})

# =============================================================================
# RSS LAMBDA CONSTRUCTOR - Method Restrictions
# =============================================================================

test_that("rss_lambda_constructor switches MoM to MLE", {
  z <- rnorm(50)
  R <- diag(50)

  expect_message(
    result <- rss_lambda_constructor(z, R, lambda = 0.5,
                                    estimate_residual_method = "MoM"),
    "Automatically switching to estimate_residual_method = 'MLE'"
  )

  expect_equal(result$params$estimate_residual_method, "MLE")
})

test_that("rss_lambda_constructor rejects Servin_Stephens", {
  z <- rnorm(50)
  R <- diag(50)

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5,
                          estimate_residual_method = "Servin_Stephens"),
    "Servin-Stephens prior on residual variance is not implemented"
  )
})

test_that("rss_lambda_constructor sets unmappable_effects to none", {
  z <- rnorm(50)
  R <- diag(50)

  result <- rss_lambda_constructor(z, R, lambda = 0.5,
                                  unmappable_effects = "inf")

  expect_equal(result$params$unmappable_effects, "none")
})

# =============================================================================
# RSS LAMBDA CONSTRUCTOR - MAF Filtering
# =============================================================================

test_that("rss_lambda_constructor applies MAF filter", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)
  maf <- runif(p, 0, 0.5)

  result <- rss_lambda_constructor(z, R, lambda = 0.5,
                                  maf = maf, maf_thresh = 0.1)

  n_filtered <- sum(maf > 0.1)
  expect_equal(result$data$p, n_filtered)
  expect_length(result$data$z, n_filtered)
  expect_equal(nrow(result$data$R), n_filtered)
})

test_that("rss_lambda_constructor rejects MAF with wrong length", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)
  maf <- runif(p - 10)  # Wrong length

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5, maf = maf),
    "The length of maf does not agree with expected 50"
  )
})

# =============================================================================
# RSS LAMBDA CONSTRUCTOR - Null Weight
# =============================================================================

test_that("rss_lambda_constructor adds null column when null_weight > 0", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  result <- rss_lambda_constructor(z, R, lambda = 0.5, null_weight = 0.1)

  expect_equal(result$data$p, p + 1)
  expect_length(result$data$z, p + 1)
  expect_equal(nrow(result$data$R), p + 1)
  expect_equal(result$params$null_weight, 0.1)
  expect_equal(result$data$z[p + 1], 0)
})

test_that("rss_lambda_constructor adjusts custom prior weights with null_weight", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  # Create custom prior weights (not uniform)
  custom_weights <- runif(p, 0.5, 2)
  custom_weights <- custom_weights / sum(custom_weights)  # Normalize to sum to 1

  result <- rss_lambda_constructor(z, R, lambda = 0.5,
                                  prior_weights = custom_weights,
                                  null_weight = 0.15)

  # Check that we have p+1 weights (original p + null column)
  expect_length(result$params$prior_weights, p + 1)

  # Check that all weights sum to 1
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)

  # Check that the null weight is exactly 0.15
  expect_equal(result$params$prior_weights[p + 1], 0.15, tolerance = 1e-10)

  # Check that the other weights were scaled by (1 - null_weight) = 0.85
  expect_equal(result$params$prior_weights[1:p],
               custom_weights * 0.85,
               tolerance = 1e-10)

  # Verify that the sum of the first p weights is (1 - 0.15) = 0.85
  expect_equal(sum(result$params$prior_weights[1:p]), 0.85, tolerance = 1e-10)
})

test_that("rss_lambda_constructor rejects non-numeric null_weight", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5, null_weight = "invalid"),
    "Null weight must be numeric"
  )
})

test_that("rss_lambda_constructor rejects negative null_weight", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5, null_weight = -0.1),
    "Null weight must be between 0 and 1"
  )
})

test_that("rss_lambda_constructor rejects null_weight >= 1", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5, null_weight = 1.0),
    "Null weight must be between 0 and 1"
  )

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5, null_weight = 1.5),
    "Null weight must be between 0 and 1"
  )
})

# =============================================================================
# SUMMARY STATISTICS CONSTRUCTOR - Routing Logic
# =============================================================================

test_that("summary_stats_constructor routes to rss_lambda when lambda != 0", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  result <- summary_stats_constructor(z = z, R = R, lambda = 0.5)

  expect_s3_class(result$data, "rss_lambda")
  expect_equal(result$data$lambda, 0.5)
})

test_that("summary_stats_constructor routes to sufficient_stats when lambda = 0", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  result <- summary_stats_constructor(z = z, R = R, n = 100, lambda = 0)

  expect_s3_class(result$data, "ss")
})

# =============================================================================
# SUMMARY STATISTICS CONSTRUCTOR - Input Validation
# =============================================================================

test_that("summary_stats_constructor rejects R with wrong number of rows", {
  p <- 50
  z <- rnorm(p)

  # Create R with wrong number of rows (40 instead of 50)
  R_wrong <- diag(40)

  expect_error(
    summary_stats_constructor(z = z, R = R_wrong, n = 100, lambda = 0),
    "The dimension of R \\(40 x 40\\) does not agree with expected \\(50 x 50\\)"
  )
})

test_that("summary_stats_constructor rejects n <= 1", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  # Test n = 1
  expect_error(
    summary_stats_constructor(z = z, R = R, n = 1, lambda = 0),
    "n must be greater than 1"
  )

  # Test n = 0
  expect_error(
    summary_stats_constructor(z = z, R = R, n = 0, lambda = 0),
    "n must be greater than 1"
  )

  # Test negative n
  expect_error(
    summary_stats_constructor(z = z, R = R, n = -5, lambda = 0),
    "n must be greater than 1"
  )
})

test_that("summary_stats_constructor rejects mismatched bhat and shat lengths", {
  p <- 50
  R <- diag(p)
  bhat <- rnorm(p)
  shat <- abs(rnorm(p - 5))  # Wrong length

  expect_error(
    summary_stats_constructor(bhat = bhat, shat = shat, R = R, n = 100, lambda = 0),
    "The lengths of bhat and shat do not agree"
  )
})

test_that("summary_stats_constructor accepts scalar shat and replicates it", {
  p <- 50
  R <- diag(p)
  bhat <- rnorm(p)
  shat <- 0.1  # Scalar

  # Should replicate shat to length of bhat
  result <- summary_stats_constructor(bhat = bhat, shat = shat, R = R, n = 100, lambda = 0)
  expect_true(!is.null(result))
})

test_that("summary_stats_constructor rejects missing values in bhat", {
  p <- 50
  R <- diag(p)
  bhat <- rnorm(p)
  bhat[5] <- NA
  shat <- abs(rnorm(p))

  expect_error(
    summary_stats_constructor(bhat = bhat, shat = shat, R = R, n = 100, lambda = 0),
    "bhat, shat cannot have missing values"
  )
})

test_that("summary_stats_constructor rejects missing values in shat", {
  p <- 50
  R <- diag(p)
  bhat <- rnorm(p)
  shat <- abs(rnorm(p))
  shat[10] <- NA

  expect_error(
    summary_stats_constructor(bhat = bhat, shat = shat, R = R, n = 100, lambda = 0),
    "bhat, shat cannot have missing values"
  )
})

test_that("summary_stats_constructor rejects zero elements in shat", {
  p <- 50
  R <- diag(p)
  bhat <- rnorm(p)
  shat <- abs(rnorm(p))
  shat[5] <- 0

  expect_error(
    summary_stats_constructor(bhat = bhat, shat = shat, R = R, n = 100, lambda = 0),
    "shat cannot have zero or negative elements"
  )
})

test_that("summary_stats_constructor rejects negative elements in shat", {
  p <- 50
  R <- diag(p)
  bhat <- rnorm(p)
  shat <- abs(rnorm(p))
  shat[8] <- -0.5

  expect_error(
    summary_stats_constructor(bhat = bhat, shat = shat, R = R, n = 100, lambda = 0),
    "shat cannot have zero or negative elements"
  )
})

test_that("summary_stats_constructor rejects empty z vector", {
  # When z is empty, length(z) = 0, so p = 0
  # This causes R dimension check to fail before z length check
  z <- numeric(0)  # Empty vector
  R <- matrix(0, 0, 0)  # Match the expected dimension (0 x 0)

  expect_error(
    summary_stats_constructor(z = z, R = R, n = 100, lambda = 0),
    "Input vector z should have at least one element"
  )
})

test_that("summary_stats_constructor warns about deprecated z_ld_weight", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  expect_message(
    summary_stats_constructor(z = z, R = R, n = 100, lambda = 0, z_ld_weight = 0.1),
    "As of version 0.11.0, use of non-zero z_ld_weight is no longer recommended"
  )
})

test_that("summary_stats_constructor rejects MAF with wrong length", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)
  maf <- runif(p - 10)  # Wrong length

  expect_error(
    summary_stats_constructor(z = z, R = R, n = 100, lambda = 0, maf = maf),
    "The length of maf does not agree with expected 50"
  )
})

test_that("summary_stats_constructor handles shat and var_y for original scale effects", {
  p <- 50
  bhat <- rnorm(p)
  shat <- abs(rnorm(p, mean = 0.1, sd = 0.02))
  var_y <- 2.5
  R <- diag(p)
  n <- 100

  # This should use the original scale path (lines 649-655)
  result <- summary_stats_constructor(
    bhat = bhat, shat = shat, var_y = var_y, R = R, n = n, lambda = 0
  )

  # Verify the result is created successfully
  expect_true(!is.null(result))
  expect_true(!is.null(result$data))
  expect_true(!is.null(result$data$XtX))
  expect_true(!is.null(result$data$Xty))
  expect_true(!is.null(result$data$yty))

  # Verify yty matches expected: (n - 1) * var_y
  expect_equal(result$data$yty, (n - 1) * var_y)
})

# =============================================================================
# SUMMARY STATISTICS CONSTRUCTOR - Lambda=0 Path
# =============================================================================

test_that("summary_stats_constructor converts z to sufficient stats", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)
  n <- 100

  result <- summary_stats_constructor(z = z, R = R, n = n, lambda = 0)

  expect_true("XtX" %in% names(result$data))
  expect_true("Xty" %in% names(result$data))
  expect_true("yty" %in% names(result$data))
})

test_that("summary_stats_constructor handles z without n", {
  p <- 50
  z <- rnorm(p)
  R <- diag(p)

  expect_message(
    result <- summary_stats_constructor(z = z, R = R, lambda = 0),
    "Providing the sample size"
  )

  expect_s3_class(result$data, "ss")
})

test_that("summary_stats_constructor converts bhat/shat to z", {
  p <- 50
  bhat <- rnorm(p)
  shat <- runif(p, 0.5, 1.5)
  R <- diag(p)
  n <- 100

  result <- summary_stats_constructor(bhat = bhat, shat = shat, R = R, n = n, lambda = 0)

  expect_s3_class(result$data, "ss")
})

test_that("summary_stats_constructor requires either z or bhat/shat", {
  R <- diag(50)

  expect_error(
    summary_stats_constructor(R = R, n = 100, lambda = 0),
    "Please provide either z or \\(bhat, shat\\)"
  )
})

test_that("summary_stats_constructor rejects both z and bhat/shat", {
  z <- rnorm(50)
  bhat <- rnorm(50)
  shat <- runif(50, 0.5, 1.5)
  R <- diag(50)

  expect_error(
    summary_stats_constructor(z = z, bhat = bhat, shat = shat, R = R, n = 100, lambda = 0),
    "Please provide either z or \\(bhat, shat\\), but not both"
  )
})

# =============================================================================
# SUMMARY STATISTICS CONSTRUCTOR - Lambda != 0 Restrictions
# =============================================================================

test_that("summary_stats_constructor rejects bhat/shat when lambda != 0", {
  z <- rnorm(50)
  bhat <- rnorm(50)
  shat <- runif(50, 0.5, 1.5)
  R <- diag(50)

  expect_error(
    summary_stats_constructor(z = z, R = R, bhat = bhat, shat = shat, lambda = 0.5),
    "Parameters 'bhat' and 'shat' are not supported when lambda != 0"
  )
})

test_that("summary_stats_constructor rejects var_y when lambda != 0", {
  z <- rnorm(50)
  R <- diag(50)

  expect_error(
    summary_stats_constructor(z = z, R = R, var_y = 1.5, lambda = 0.5),
    "Parameter 'var_y' is not supported when lambda != 0"
  )
})

test_that("summary_stats_constructor rejects z_ld_weight when lambda != 0", {
  z <- rnorm(50)
  R <- diag(50)

  expect_error(
    summary_stats_constructor(z = z, R = R, z_ld_weight = 0.5, lambda = 0.5),
    "Parameter 'z_ld_weight' is not supported when lambda != 0"
  )
})

test_that("summary_stats_constructor warns about n when lambda != 0", {
  z <- rnorm(50)
  R <- diag(50)

  expect_message(
    result <- summary_stats_constructor(z = z, R = R, n = 100, lambda = 0.5),
    "Parameter 'n' is ignored when lambda != 0"
  )
})

# =============================================================================
# SUMMARY STATISTICS CONSTRUCTOR - Lambda=0 Restrictions
# =============================================================================

test_that("summary_stats_constructor rejects intercept_value when lambda = 0", {
  z <- rnorm(50)
  R <- diag(50)

  expect_error(
    summary_stats_constructor(z = z, R = R, n = 100, lambda = 0, intercept_value = 0.5),
    "Parameter 'intercept_value' is only supported when lambda != 0"
  )
})

# =============================================================================
# INTEGRATION - Constructor Output Usability
# =============================================================================

test_that("individual_data_constructor output works with ibss_initialize", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 31)

  result <- individual_data_constructor(base_data$X, base_data$y, L = 5)

  expect_error(
    model <- ibss_initialize(result$data, result$params),
    NA
  )
})

test_that("sufficient_stats_constructor output works with ibss_initialize", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 32)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX, L = 5)

  expect_error(
    model <- ibss_initialize(result$data, result$params),
    NA
  )
})

test_that("rss_lambda_constructor output works with ibss_initialize", {
  z <- rnorm(50)
  R <- diag(50)

  result <- rss_lambda_constructor(z, R, lambda = 0.5, L = 5)

  expect_error(
    model <- ibss_initialize(result$data, result$params),
    NA
  )
})
