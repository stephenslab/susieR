devtools::load_all(".")


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
# SUFFICIENT STATISTICS CONSTRUCTOR - Basic Functionality
# =============================================================================

test_that("sufficient_stats_constructor returns correct structure", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 20)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(XtX, Xty, yty, base_data$n)

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

  result <- sufficient_stats_constructor(XtX, Xty, yty, base_data$n)

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

  result <- sufficient_stats_constructor(XtX, Xty, yty, base_data$n)

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
    sufficient_stats_constructor(XtX, Xty, yty),
    "n must be provided"
  )
})

test_that("sufficient_stats_constructor rejects n <= 1", {
  XtX <- matrix(1:25, 5, 5)
  Xty <- 1:5
  yty <- 10

  expect_error(
    sufficient_stats_constructor(XtX, Xty, yty, n = 1),
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

test_that("sufficient_stats_constructor rejects dimension mismatch", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 23)
  XtX <- crossprod(base_data$X)
  Xty <- rnorm(10)
  yty <- 100

  expect_error(
    sufficient_stats_constructor(XtX, Xty, yty, base_data$n),
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
    result <- sufficient_stats_constructor(XtX, Xty, yty, n),
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
    sufficient_stats_constructor(XtX, Xty, yty, n),
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
    result <- sufficient_stats_constructor(XtX, Xty, yty, base_data$n),
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
    sufficient_stats_constructor(XtX, Xty, yty, n),
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

  result <- sufficient_stats_constructor(XtX, Xty, yty, n = base_data$n, standardize = TRUE)

  d_attr <- attr(result$data$XtX, "d")
  expect_length(d_attr, base_data$p)
  expect_true(all(is.finite(d_attr)))
})

test_that("sufficient_stats_constructor does not standardize when standardize=FALSE", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 26)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(XtX, Xty, yty, n = base_data$n, standardize = FALSE)

  csd_attr <- attr(result$data$XtX, "scaled:scale")
  expect_true(all(csd_attr == 1))
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
    XtX, Xty, yty, n = base_data$n,
    maf = maf, maf_thresh = 0.1
  )

  n_filtered <- sum(maf > 0.1)
  expect_equal(result$data$p, n_filtered)
  expect_equal(nrow(result$data$XtX), n_filtered)
  expect_length(result$data$Xty, n_filtered)
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
    XtX, Xty, yty, n = base_data$n,
    null_weight = 0.1,
    X_colmeans = rep(0, base_data$p)
  )

  expect_equal(result$data$p, base_data$p + 1)
  expect_equal(nrow(result$data$XtX), base_data$p + 1)
  expect_length(result$data$Xty, base_data$p + 1)
  expect_equal(result$params$null_weight, 0.1)
})

# =============================================================================
# SUFFICIENT STATISTICS CONSTRUCTOR - Method Restrictions
# =============================================================================

test_that("sufficient_stats_constructor rejects Servin_Stephens", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 29)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- 100

  expect_error(
    sufficient_stats_constructor(XtX, Xty, yty, n = base_data$n,
                                estimate_residual_method = "Servin_Stephens"),
    "Small sample correction not implemented for SS/RSS data"
  )
})

test_that("sufficient_stats_constructor rejects unmappable_effects='ash'", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 30)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- 100

  expect_error(
    sufficient_stats_constructor(XtX, Xty, yty, n = base_data$n,
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
    "Servin Stephens small sample correction is not implemented"
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

  result <- sufficient_stats_constructor(XtX, Xty, yty, n = base_data$n, L = 5)

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
