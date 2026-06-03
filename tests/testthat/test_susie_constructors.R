context("SuSiE Data Constructors")

# =============================================================================
# INDIVIDUAL DATA CONSTRUCTOR
# =============================================================================

test_that("individual_data_constructor returns correct structure", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 1)
  result <- individual_data_constructor(base_data$X, base_data$y)

  expect_type(result, "list")
  expect_named(result, c("data", "params"), ignore.order = TRUE)
  expect_s3_class(result$data, "individual")
  expect_true(all(c("X", "y", "n", "p", "mean_y") %in% names(result$data)))
})

test_that("individual_data_constructor sets correct dimensions", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 3)
  result <- individual_data_constructor(base_data$X, base_data$y)

  expect_equal(result$data$n, base_data$n)
  expect_equal(result$data$p, base_data$p)
  expect_equal(dim(result$data$X), c(base_data$n, base_data$p))
  expect_length(result$data$y, base_data$n)
})

test_that("individual_data_constructor sets X standardization attributes", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 4)
  result <- individual_data_constructor(base_data$X, base_data$y, standardize = TRUE, intercept = TRUE)

  expect_length(attr(result$data$X, "scaled:center"), base_data$p)
  expect_length(attr(result$data$X, "scaled:scale"), base_data$p)
  expect_length(attr(result$data$X, "d"), base_data$p)
  expect_true(all(attr(result$data$X, "scaled:scale") > 0))
})

# --- Input validation ---

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

test_that("individual_data_constructor removes NA rows in y when na.rm=TRUE", {
  base_data <- generate_base_data(n = 10, p = 10, k = 0, seed = 7)
  base_data$y[5] <- NA

  result <- individual_data_constructor(base_data$X, base_data$y, na.rm = TRUE)

  expect_equal(result$data$n, 9)
  expect_equal(nrow(result$data$X), 9)
  expect_length(result$data$y, 9)
  expect_false(anyNA(result$data$y))
})

test_that("individual_data_constructor computes residual_variance_lowerbound after NA removal", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 42)
  base_data$y[1] <- NA

  result <- individual_data_constructor(base_data$X, base_data$y, na.rm = TRUE)

  y_clean <- base_data$y[!is.na(base_data$y)]
  expect_equal(result$params$residual_variance_lowerbound, var(y_clean) / 1e4,
               tolerance = 1e-8)
})

test_that("individual_data_constructor respects custom residual_variance_lowerbound with NA in y", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 43)
  base_data$y[1] <- NA

  result <- individual_data_constructor(
    base_data$X, base_data$y,
    na.rm = TRUE,
    residual_variance_lowerbound = 0.001
  )

  expect_equal(result$params$residual_variance_lowerbound, 0.001)
})

test_that("individual_data_constructor handles multiple NAs in y with na.rm=TRUE", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 44)
  base_data$y[c(1, 25, 50, 75, 100)] <- NA

  result <- individual_data_constructor(base_data$X, base_data$y, na.rm = TRUE)

  expect_equal(result$data$n, 95)
  expect_equal(nrow(result$data$X), 95)
  expect_length(result$data$y, 95)
  expect_false(anyNA(result$data$y))
  expect_true(is.finite(result$params$residual_variance_lowerbound))
})

# --- Centering and scaling ---

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
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p, mean = 5, sd = 3), n, p)
  y <- rnorm(n)

  result <- individual_data_constructor(X, y, standardize = TRUE, intercept = TRUE)

  cm  <- attr(result$data$X, "scaled:center")
  csd <- attr(result$data$X, "scaled:scale")
  expect_length(cm, p)
  expect_length(csd, p)
  expect_true(all(csd > 0))
})

# --- Prior weights ---

test_that("individual_data_constructor creates uniform prior weights by default", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 11)
  result <- individual_data_constructor(base_data$X, base_data$y)

  expect_length(result$params$prior_weights, base_data$p)
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)
  expect_equal(result$params$prior_weights, rep(1 / base_data$p, base_data$p),
               tolerance = 1e-10)
})

test_that("individual_data_constructor normalizes custom prior weights", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 12)
  result <- individual_data_constructor(base_data$X, base_data$y,
                                        prior_weights = rep(2, base_data$p))

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

# --- Null weight ---

test_that("individual_data_constructor handles null_weight=0 (no null column)", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 15)
  result <- individual_data_constructor(base_data$X, base_data$y, null_weight = 0)

  expect_equal(result$data$p, base_data$p)
  expect_equal(ncol(result$data$X), base_data$p)
  expect_null(result$params$null_weight)
})

test_that("individual_data_constructor adds null column and adjusts weights when null_weight > 0", {
  set.seed(16)
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = NULL)

  result <- individual_data_constructor(base_data$X, base_data$y, null_weight = 0.1)

  expect_equal(result$data$p, base_data$p + 1)
  expect_equal(ncol(result$data$X), base_data$p + 1)
  expect_equal(result$params$null_weight, 0.1)
  expect_true(all(result$data$X[, base_data$p + 1] == 0))
  expect_length(result$params$prior_weights, base_data$p + 1)
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)
  expect_equal(result$params$prior_weights[base_data$p + 1], 0.1, tolerance = 1e-10)
})

test_that("individual_data_constructor scales custom prior weights with null_weight", {
  set.seed(17)
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = NULL)
  custom_weights <- runif(base_data$p, 0.5, 2)
  custom_weights <- custom_weights / sum(custom_weights)

  result <- individual_data_constructor(base_data$X, base_data$y,
                                        prior_weights = custom_weights,
                                        null_weight = 0.3)

  expect_length(result$params$prior_weights, base_data$p + 1)
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)
  expect_equal(result$params$prior_weights[base_data$p + 1], 0.3, tolerance = 1e-10)
  expect_equal(result$params$prior_weights[1:base_data$p], custom_weights * 0.7,
               tolerance = 1e-10)
})

test_that("individual_data_constructor rejects invalid null_weight values", {
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

# --- Rfast warning ---

test_that("individual_data_constructor warns about Rfast when p > 1000 and Rfast not available", {
  skip_if(requireNamespace("Rfast", quietly = TRUE),
          "Rfast is installed, skipping warning test")

  base_data <- generate_base_data(n = 100, p = 1001, k = 0, seed = 19)

  expect_message(
    result <- individual_data_constructor(base_data$X, base_data$y),
    "consider installing the Rfast package",
    fixed = FALSE
  )
  expect_equal(result$data$p, 1001)
})

test_that("individual_data_constructor does not warn when p <= 1000", {
  base_data <- generate_base_data(n = 100, p = 1000, k = 0, seed = 20)
  result <- individual_data_constructor(base_data$X, base_data$y)

  expect_equal(result$data$p, 1000)
})

# --- Parameters ---

test_that("individual_data_constructor stores all parameters", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 21)

  result <- individual_data_constructor(
    base_data$X, base_data$y,
    L = 5, estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE, max_iter = 50, tol = 1e-4
  )

  expect_equal(result$params$L, 5)
  expect_true(result$params$estimate_residual_variance)
  expect_true(result$params$estimate_prior_variance)
  expect_equal(result$params$max_iter, 50)
  expect_equal(result$params$tol, 1e-4)
})

# --- Incompatible combinations ---

test_that("individual_data_constructor rejects unmappable_effects with NIG", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 22)

  for (ue in c("inf", "ash")) {
    expect_error(
      individual_data_constructor(
        base_data$X, base_data$y,
        unmappable_effects = ue,
        estimate_residual_method = "NIG"
      ),
      sprintf("The combination of unmappable_effects = '%s' with estimate_residual_method = 'NIG' is not supported", ue)
    )
  }
})

test_that("individual_data_constructor rejects unmappable_effects='ash' with estimate_prior_method='EM'", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 23)

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
# SUFFICIENT STATISTICS CONSTRUCTOR
# =============================================================================

test_that("sufficient_stats_constructor returns correct structure with correct dimensions", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 30)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  result <- sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX)

  expect_type(result, "list")
  expect_named(result, c("data", "params"), ignore.order = TRUE)
  expect_s3_class(result$data, "ss")
  expect_true(all(c("XtX", "Xty", "yty", "n", "p") %in% names(result$data)))
  expect_equal(result$data$n, base_data$n)
  expect_equal(result$data$p, base_data$p)
  expect_equal(dim(result$data$XtX), c(base_data$p, base_data$p))
  expect_length(result$data$Xty, base_data$p)
})

# --- Input validation ---

test_that("sufficient_stats_constructor requires n and rejects n <= 1", {
  XtX <- matrix(1:25, 5, 5); Xty <- 1:5; yty <- 10

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, XtX = XtX),
    "n must be provided"
  )
  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = 1, XtX = XtX),
    "n must be greater than 1"
  )
})

test_that("sufficient_stats_constructor requires XtX, Xty, yty together", {
  expect_error(
    sufficient_stats_constructor(Xty = 1:5, yty = 10, n = 100),
    "XtX, Xty, yty must all be provided"
  )
})

test_that("sufficient_stats_constructor rejects invalid XtX types", {
  Xty <- rnorm(5); yty <- 10; n <- 100

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n,
                                 XtX = data.frame(matrix(rnorm(25), 5, 5))),
    "XtX must be a numeric dense or sparse matrix"
  )
  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n,
                                 XtX = matrix(1L:25L, 5, 5)),
    "XtX must be a numeric dense or sparse matrix"
  )
  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n,
                                 XtX = matrix(as.character(1:25), 5, 5)),
    "XtX must be a numeric dense or sparse matrix"
  )
  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = yty, n = n, XtX = rnorm(25)),
    "XtX must be a numeric dense or sparse matrix"
  )
})

test_that("sufficient_stats_constructor rejects dimension mismatch between XtX and Xty", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 31)
  XtX <- crossprod(base_data$X)

  expect_error(
    sufficient_stats_constructor(Xty = rnorm(10), yty = 100, n = base_data$n, XtX = XtX),
    "does not agree with expected"
  )
})

test_that("sufficient_stats_constructor symmetrizes non-symmetric XtX with message", {
  XtX <- matrix(1:25, 5, 5)
  XtX[1, 2] <- 100

  expect_message(
    result <- sufficient_stats_constructor(Xty = 1:5, yty = 10, n = 100, XtX = XtX),
    "XtX is not symmetric"
  )
  expect_true(isSymmetric(result$data$XtX))
})

test_that("sufficient_stats_constructor rejects XtX with NAs", {
  XtX <- matrix(rnorm(25), 5, 5)
  XtX[1, 1] <- NA

  expect_error(
    sufficient_stats_constructor(Xty = 1:5, yty = 10, n = 100, XtX = XtX),
    "Input XtX matrix contains NAs"
  )
})

test_that("sufficient_stats_constructor replaces NA Xty with zero via message", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 32)
  XtX <- crossprod(base_data$X)
  Xty <- rnorm(base_data$p)
  Xty[5] <- NA

  expect_message(
    result <- sufficient_stats_constructor(Xty = Xty, yty = 100, n = base_data$n, XtX = XtX),
    "NA values in Xty are replaced with 0"
  )
  expect_false(anyNA(result$data$Xty))
  expect_equal(result$data$Xty[5], 0)
})

test_that("sufficient_stats_constructor rejects infinite Xty", {
  XtX <- crossprod(matrix(rnorm(100 * 10), 100, 10))
  Xty <- rnorm(10)
  Xty[5] <- Inf

  expect_error(
    sufficient_stats_constructor(Xty = Xty, yty = 100, n = 100, XtX = XtX),
    "Input Xty contains infinite values"
  )
})

# --- Standardization ---

test_that("sufficient_stats_constructor standardizes when standardize=TRUE", {
  set.seed(33)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p, mean = 5, sd = 3), n, p)
  y <- rnorm(n)
  ss <- compute_summary_stats(X, y)

  result <- sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n,
                                         XtX = ss$XtX, standardize = TRUE)

  d_attr <- attr(result$data$XtX, "d")
  expect_length(d_attr, p)
  expect_true(all(is.finite(d_attr)))
})

test_that("sufficient_stats_constructor sets csd=1 when standardize=FALSE", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 34)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  result <- sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n,
                                         XtX = ss$XtX, standardize = FALSE)

  expect_true(all(attr(result$data$XtX, "scaled:scale") == 1))
})

# --- Rfast warning ---

test_that("sufficient_stats_constructor warns about Rfast when p > 1000", {
  skip_if(requireNamespace("Rfast", quietly = TRUE),
          "Rfast is installed, skipping warning test")

  base_data <- generate_base_data(n = 100, p = 1001, k = 0, seed = 35)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  expect_message(
    result <- sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty,
                                           n = ss$n, XtX = ss$XtX),
    "consider installing the Rfast package",
    fixed = FALSE
  )
  expect_equal(result$data$p, 1001)
})

test_that("sufficient_stats_constructor does not warn when p <= 1000", {
  base_data <- generate_base_data(n = 100, p = 1000, k = 0, seed = 36)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  result <- sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty,
                                         n = ss$n, XtX = ss$XtX)
  expect_equal(result$data$p, 1000)
})

# --- MAF filtering ---

test_that("sufficient_stats_constructor applies MAF filter", {
  set.seed(37)
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = NULL)
  ss <- compute_summary_stats(base_data$X, base_data$y)
  maf <- runif(base_data$p, 0, 0.5)

  result <- sufficient_stats_constructor(
    Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
    maf = maf, maf_thresh = 0.1
  )

  n_filtered <- sum(maf > 0.1)
  expect_equal(result$data$p, n_filtered)
  expect_equal(nrow(result$data$XtX), n_filtered)
  expect_length(result$data$Xty, n_filtered)
})

test_that("sufficient_stats_constructor rejects MAF with incorrect length", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 38)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 maf = runif(base_data$p - 10), maf_thresh = 0.1),
    "The length of maf does not agree with expected"
  )
  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 maf = runif(base_data$p + 10), maf_thresh = 0.1),
    "The length of maf does not agree with expected"
  )
})

# --- Positive semidefinite check ---

test_that("sufficient_stats_constructor rejects non-PSD XtX when check_input=TRUE", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 39)
  XtX <- crossprod(base_data$X)
  XtX[1, 1] <- -10

  expect_error(
    sufficient_stats_constructor(
      Xty = crossprod(base_data$X, base_data$y), yty = sum(base_data$y^2),
      n = base_data$n, XtX = XtX, check_input = TRUE
    ),
    "XtX is not a positive semidefinite matrix"
  )
})

test_that("sufficient_stats_constructor accepts PSD XtX when check_input=TRUE", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 40)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  result <- sufficient_stats_constructor(
    Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
    check_input = TRUE
  )
  expect_equal(result$data$p, base_data$p)
})

test_that("sufficient_stats_constructor warns when Xty not in column space of XtX", {
  XtX  <- diag(c(1, 1, 1, 0, 0))
  Xty  <- c(1, 1, 1, 10, 10)

  expect_message(
    sufficient_stats_constructor(Xty = Xty, yty = 100, n = 100, XtX = XtX,
                                 check_input = TRUE),
    "Xty does not lie in the space of the non-zero eigenvectors"
  )
})

test_that("sufficient_stats_constructor does not warn when Xty in column space of XtX", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 41)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  result <- sufficient_stats_constructor(
    Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
    check_input = TRUE
  )
  expect_equal(result$data$p, base_data$p)
})

# --- Null weight ---

test_that("sufficient_stats_constructor adds null column when null_weight > 0", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 50)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  result <- sufficient_stats_constructor(
    Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
    null_weight = 0.1, X_colmeans = rep(0, base_data$p)
  )

  expect_equal(result$data$p, base_data$p + 1)
  expect_equal(nrow(result$data$XtX), base_data$p + 1)
  expect_length(result$data$Xty, base_data$p + 1)
  expect_equal(result$params$null_weight, 0.1)
})

test_that("sufficient_stats_constructor scales custom prior weights with null_weight", {
  set.seed(51)
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = NULL)
  ss <- compute_summary_stats(base_data$X, base_data$y)
  custom_weights <- runif(base_data$p, 0.5, 2)
  custom_weights <- custom_weights / sum(custom_weights)

  result <- sufficient_stats_constructor(
    Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
    prior_weights = custom_weights, null_weight = 0.25,
    X_colmeans = rep(0, base_data$p)
  )

  expect_length(result$params$prior_weights, base_data$p + 1)
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)
  expect_equal(result$params$prior_weights[base_data$p + 1], 0.25, tolerance = 1e-10)
  expect_equal(result$params$prior_weights[1:base_data$p], custom_weights * 0.75,
               tolerance = 1e-10)
})

test_that("sufficient_stats_constructor rejects invalid null_weight values", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 52)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 null_weight = "invalid"),
    "Null weight must be numeric"
  )
  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 null_weight = -0.1),
    "Null weight must be between 0 and 1"
  )
  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 null_weight = 1),
    "Null weight must be between 0 and 1"
  )
  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 null_weight = 1.5),
    "Null weight must be between 0 and 1"
  )
})

test_that("sufficient_stats_constructor replicates scalar X_colmeans", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 53)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  result_with_nw <- sufficient_stats_constructor(
    Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
    null_weight = 0.1, X_colmeans = 0
  )
  expect_equal(result_with_nw$data$p, base_data$p + 1)

  result_no_nw <- sufficient_stats_constructor(
    Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
    X_colmeans = 0
  )
  expect_equal(result_no_nw$data$p, base_data$p)
})

test_that("sufficient_stats_constructor rejects wrong length X_colmeans", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 54)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 null_weight = 0.1,
                                 X_colmeans = rep(0, base_data$p - 10)),
    "The length of X_colmeans does not agree with number of variables"
  )
  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 X_colmeans = rep(0, base_data$p - 10)),
    "X_colmeans.*does not match number of variables"
  )
})

test_that("sufficient_stats_constructor rejects wrong length or all-zero prior_weights", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 55)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 prior_weights = rep(1, base_data$p - 10)),
    "Prior weights must have length p"
  )
  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 prior_weights = rep(0, base_data$p)),
    "Prior weight should be greater than 0 for at least one variable"
  )
})

# --- Method restrictions ---

test_that("sufficient_stats_constructor accepts NIG residual method", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 56)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  result <- sufficient_stats_constructor(Xty = ss$Xty, yty = 100, n = ss$n,
                                         XtX = ss$XtX,
                                         estimate_residual_method = "NIG")
  expect_true(result$params$use_NIG)
  expect_equal(result$params$estimate_prior_method, "EM")
})

test_that("sufficient_stats_constructor accepts unmappable_effects='ash'", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 57)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  result <- sufficient_stats_constructor(Xty = ss$Xty, yty = 100, n = ss$n,
                                         XtX = ss$XtX,
                                         unmappable_effects = "ash")
  expect_s3_class(result$data, "ss")
  expect_equal(result$params$unmappable_effects, "ash")
})

# --- X-path branches ---

test_that("sufficient_stats_constructor errors on ncol(X)/Xty dimension mismatch", {
  set.seed(58)
  n <- 30; p_true <- 10
  X <- matrix(rnorm(n * p_true), n, p_true)

  expect_error(
    sufficient_stats_constructor(X = X, Xty = rnorm(7), yty = 100, n = n),
    "does not agree with the length of Xty"
  )
})

test_that("sufficient_stats_constructor X-path adds null column when null_weight > 0", {
  set.seed(59)
  n <- 30; p <- 10
  X <- matrix(rnorm(n * p), n, p)

  result <- sufficient_stats_constructor(
    X = X, Xty = rnorm(p), yty = 100, n = n,
    null_weight = 0.1, X_colmeans = rep(0, p)
  )

  expect_equal(result$data$p, p + 1)
  expect_equal(result$params$null_weight, 0.1)
})

test_that("sufficient_stats_constructor X-path sets csd=1 when standardize=FALSE", {
  set.seed(60)
  n <- 30; p <- 10
  X <- matrix(rnorm(n * p), n, p)

  result <- sufficient_stats_constructor(
    X = X, Xty = rnorm(p), yty = 100, n = n,
    standardize = FALSE
  )

  expect_true(all(attr(result$data$X, "scaled:scale") == 1))
  expect_equal(result$data$p, p)
})

# =============================================================================
# RSS LAMBDA CONSTRUCTOR
# =============================================================================

test_that("rss_lambda_constructor returns correct structure with all required fields", {
  set.seed(70)
  p <- 50; z <- rnorm(p); R <- diag(p)

  result <- rss_lambda_constructor(z, R, lambda = 0.5)

  expect_type(result, "list")
  expect_named(result, c("data", "params"), ignore.order = TRUE)
  expect_s3_class(result$data, "rss_lambda")
  expect_true(all(c("z", "R", "lambda", "eigen_R", "Vtz", "n", "p") %in%
                    names(result$data)))
})

test_that("rss_lambda_constructor stores n and dimensions correctly", {
  set.seed(71)
  p <- 50; z <- rnorm(p); R <- diag(p)

  result <- rss_lambda_constructor(z, R, lambda = 0.5)
  expect_true(is.na(result$data$n))
  expect_equal(result$data$p, p)
  expect_length(result$data$z, p)
  expect_equal(dim(result$data$R), c(p, p))

  result_n <- rss_lambda_constructor(z, R, lambda = 0.5, n = 1000)
  expect_equal(result_n$data$n, 1000L)
  adj <- (1000 - 1) / (z^2 + 1000 - 2)
  expect_equal(result_n$data$z, sqrt(adj) * z)
})

test_that("rss_lambda_constructor computes eigen decomposition and Vtz", {
  set.seed(72)
  p <- 50; z <- rnorm(p); R <- diag(p)

  result <- rss_lambda_constructor(z, R, lambda = 0.5)

  expect_true(all(c("values", "vectors") %in% names(result$data$eigen_R)))
  expect_length(result$data$eigen_R$values, p)
  expect_length(result$data$Vtz, p)
})

# --- Input validation ---

test_that("rss_lambda_constructor rejects dimension mismatch between z and R", {
  expect_error(
    rss_lambda_constructor(rnorm(50), diag(40), lambda = 0.5),
    "does not agree with expected"
  )
})

test_that("rss_lambda_constructor errors on non-symmetric / non-PSD R", {
  R_nonsym <- matrix(rnorm(25), 5, 5)
  z <- rnorm(5)

  expect_error(
    rss_lambda_constructor(z, R_nonsym, lambda = 0.5),
    "not a positive semidefinite matrix|R is not a symmetric matrix"
  )
})

test_that("rss_lambda_constructor accepts R with various dimname configurations", {
  set.seed(73)
  p <- 10; z <- rnorm(p)

  R_mismatch <- diag(p)
  rownames(R_mismatch) <- paste0("row_", 1:p)
  colnames(R_mismatch) <- paste0("col_", 1:p)
  expect_equal(rss_lambda_constructor(z, R_mismatch, lambda = 0.5)$data$p, p)

  R_matching <- diag(p)
  rownames(R_matching) <- paste0("SNP", 1:p)
  colnames(R_matching) <- paste0("SNP", 1:p)
  expect_equal(rss_lambda_constructor(z, R_matching, lambda = 0.5)$data$p, p)

  R_none <- diag(p)
  expect_equal(rss_lambda_constructor(z, R_none, lambda = 0.5)$data$p, p)
})

test_that("rss_lambda_constructor rejects integer matrix R", {
  R <- matrix(1:25, 5, 5)
  R <- R + t(R)
  mode(R) <- "integer"

  expect_error(
    rss_lambda_constructor(rnorm(5), R, lambda = 0.5),
    "Input R must be a double-precision matrix or a sparse matrix"
  )
})

test_that("rss_lambda_constructor rejects non-PSD R when check_R=TRUE", {
  R <- diag(5)
  R[1, 1] <- -1

  expect_error(
    rss_lambda_constructor(rnorm(5), R, lambda = 0.5, check_R = TRUE),
    "is not a positive semidefinite matrix"
  )
})

test_that("rss_lambda_constructor accepts non-PSD R when check_R=FALSE", {
  R <- diag(5)
  R[1, 1] <- -0.5

  result <- suppressWarnings(
    rss_lambda_constructor(rnorm(5), R, lambda = 0.5, check_R = FALSE)
  )
  expect_s3_class(result$data, "rss_lambda")
})

test_that("rss_lambda_constructor warns when z not in column space of R", {
  R <- diag(c(1, 1, 1, 0, 0))
  z <- c(0.1, 0.1, 0.1, 10, 10)

  expect_message(
    rss_lambda_constructor(z, R, lambda = 0.5, check_z = TRUE),
    "Input z does not lie in the space of non-zero eigenvectors of R"
  )
})

test_that("rss_lambda_constructor messages when z in column space of R", {
  R <- diag(c(1, 1, 1, 0, 0))
  z <- c(1, 2, 3, 0, 0)

  expect_message(
    suppressWarnings(rss_lambda_constructor(z, R, lambda = 0.5, check_z = TRUE)),
    "Input z is in space spanned by the non-zero eigenvectors of R"
  )
})

test_that("rss_lambda_constructor skips z check when check_z=FALSE or R is full rank", {
  R_rank_deficient <- diag(c(1, 1, 1, 0, 0))
  z_null_space <- c(0.1, 0.1, 0.1, 10, 10)

  result <- suppressWarnings(suppressMessages(
    rss_lambda_constructor(z_null_space, R_rank_deficient, lambda = 0.5, check_z = FALSE)
  ))
  expect_s3_class(result$data, "rss_lambda")

  R_full <- diag(5)
  result2 <- suppressWarnings(
    rss_lambda_constructor(rnorm(5), R_full, lambda = 0.5, check_z = TRUE)
  )
  expect_s3_class(result2$data, "rss_lambda")
})

test_that("rss_lambda_constructor rejects R with NAs", {
  R <- diag(10)
  R[1, 1] <- NA

  expect_error(
    rss_lambda_constructor(rnorm(10), R, lambda = 0.5),
    "R matrix contains missing values"
  )
})

test_that("rss_lambda_constructor rejects infinite z", {
  z <- rnorm(10)
  z[5] <- Inf

  expect_error(
    rss_lambda_constructor(z, diag(10), lambda = 0.5),
    "z contains infinite values"
  )
})

test_that("rss_lambda_constructor replaces NA z with zero via message", {
  z <- rnorm(10)
  z[5] <- NA

  expect_message(
    result <- rss_lambda_constructor(z, diag(10), lambda = 0.5),
    "NA values in z-scores are replaced with 0"
  )
  expect_false(anyNA(result$data$z))
  expect_equal(result$data$z[5], 0)
})

# --- Lambda parameter ---

test_that("rss_lambda_constructor stores lambda and adjusts residual_variance", {
  set.seed(74)
  z <- rnorm(50); R <- diag(50)

  result_stored <- rss_lambda_constructor(z, R, lambda = 0.3)
  expect_equal(result_stored$data$lambda, 0.3)

  result_adj <- rss_lambda_constructor(z, R, lambda = 0.2, residual_variance = 0.8)
  expect_equal(result_adj$params$residual_variance, 0.6)
})

test_that("rss_lambda_constructor estimates lambda when lambda='estimate'", {
  set.seed(75)
  p <- 50; z <- rnorm(p); R <- diag(p)
  R[1:10, 1:10] <- 0

  result <- rss_lambda_constructor(z, R, lambda = "estimate")

  expect_true(is.numeric(result$data$lambda))
  expect_true(result$data$lambda >= 0)
})

test_that("rss_lambda_constructor sets lambda=0 when R is full rank and lambda='estimate'", {
  set.seed(123)
  z <- rnorm(50); R <- diag(50)

  result <- rss_lambda_constructor(z, R, lambda = "estimate")
  expect_equal(result$data$lambda, 0)
})

# --- Method restrictions ---

test_that("rss_lambda_constructor rejects non-MLE residual variance methods", {
  z <- rnorm(50); R <- diag(50)

  for (method in c("MoM", "NIG")) {
    expect_error(
      rss_lambda_constructor(z, R, lambda = 0.5, estimate_residual_method = method),
      "RSS-lambda supports estimate_residual_method"
    )
  }
})

test_that("rss_lambda_constructor does not expose unmappable_effects", {
  expect_error(
    rss_lambda_constructor(rnorm(50), diag(50), lambda = 0.5, unmappable_effects = "inf"),
    "unused argument"
  )
})

# --- MAF filtering ---

test_that("rss_lambda_constructor applies MAF filter", {
  set.seed(76)
  p <- 50; z <- rnorm(p); R <- diag(p)
  maf <- runif(p, 0, 0.5)

  result <- rss_lambda_constructor(z, R, lambda = 0.5, maf = maf, maf_thresh = 0.1)

  n_filtered <- sum(maf > 0.1)
  expect_equal(result$data$p, n_filtered)
  expect_length(result$data$z, n_filtered)
  expect_equal(nrow(result$data$R), n_filtered)
})

test_that("rss_lambda_constructor rejects MAF with wrong length", {
  p <- 50

  expect_error(
    rss_lambda_constructor(rnorm(p), diag(p), lambda = 0.5, maf = runif(p - 10)),
    "The length of maf does not agree with expected 50"
  )
})

# --- Null weight ---

test_that("rss_lambda_constructor adds null column and adjusts weights when null_weight > 0", {
  set.seed(77)
  p <- 50; z <- rnorm(p); R <- diag(p)

  result <- rss_lambda_constructor(z, R, lambda = 0.5, null_weight = 0.1)

  expect_equal(result$data$p, p + 1)
  expect_length(result$data$z, p + 1)
  expect_equal(nrow(result$data$R), p + 1)
  expect_equal(result$params$null_weight, 0.1)
  expect_equal(result$data$z[p + 1], 0)
})

test_that("rss_lambda_constructor scales custom prior weights with null_weight", {
  set.seed(78)
  p <- 50; z <- rnorm(p); R <- diag(p)
  custom_weights <- runif(p, 0.5, 2)
  custom_weights <- custom_weights / sum(custom_weights)

  result <- rss_lambda_constructor(z, R, lambda = 0.5,
                                   prior_weights = custom_weights,
                                   null_weight = 0.15)

  expect_length(result$params$prior_weights, p + 1)
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-10)
  expect_equal(result$params$prior_weights[p + 1], 0.15, tolerance = 1e-10)
  expect_equal(result$params$prior_weights[1:p], custom_weights * 0.85,
               tolerance = 1e-10)
})

test_that("rss_lambda_constructor rejects invalid null_weight values", {
  p <- 50; z <- rnorm(p); R <- diag(p)

  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5, null_weight = "invalid"),
    "Null weight must be numeric"
  )
  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5, null_weight = -0.1),
    "Null weight must be between 0 and 1"
  )
  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5, null_weight = 1.0),
    "Null weight must be between 0 and 1"
  )
  expect_error(
    rss_lambda_constructor(z, R, lambda = 0.5, null_weight = 1.5),
    "Null weight must be between 0 and 1"
  )
})

# --- R-path and X-path branches ---

test_that("rss_lambda_constructor errors when neither R nor X is provided", {
  set.seed(79)
  expect_error(
    rss_lambda_constructor(z = rnorm(10), lambda = 0.1),
    "Please provide either R \\(correlation matrix\\) or X \\(factor matrix\\)."
  )
})

test_that("rss_lambda_constructor errors when both R and X are provided", {
  set.seed(80)
  p <- 8; n <- 100
  expect_error(
    rss_lambda_constructor(z = rnorm(p), R = diag(p),
                           X = matrix(rnorm(n * p), n, p), lambda = 0.1),
    "Please provide either R or X, but not both."
  )
})

test_that("rss_lambda_constructor X-path errors when ncol(X) != length(z)", {
  set.seed(81)
  expect_error(
    rss_lambda_constructor(z = rnorm(10), X = matrix(rnorm(100 * 7), 100, 7), lambda = 0.1),
    "The number of columns of X"
  )
})

test_that("rss_lambda_constructor X-path applies MAF filter to columns", {
  set.seed(82)
  p <- 15; n <- 100
  X <- matrix(rnorm(n * p), n, p)
  maf <- c(rep(0.01, 5), rep(0.30, 10))
  n_keep <- sum(maf > 0.05)

  result <- rss_lambda_constructor(z = rnorm(p), X = X, lambda = 0.1,
                                   maf = maf, maf_thresh = 0.05, max_iter = 5)
  expect_equal(result$data$p, n_keep)
})

test_that("rss_lambda_constructor X-path adds null column when null_weight > 0", {
  set.seed(83)
  p <- 10; n <- 80
  X <- matrix(rnorm(n * p), n, p)

  result <- rss_lambda_constructor(z = rnorm(p), X = X, lambda = 0.1,
                                   null_weight = 0.15, max_iter = 5)
  expect_equal(result$data$p, p + 1L)
  expect_equal(result$params$null_weight, 0.15)
})

test_that("rss_lambda_constructor maf + prior_weights subsetting preserves normalization", {
  set.seed(84)
  p <- 20
  prior_weights <- runif(p); prior_weights <- prior_weights / sum(prior_weights)
  maf <- c(rep(0.01, 8), rep(0.3, 12))
  n_keep <- sum(maf > 0.05)

  result <- rss_lambda_constructor(
    z = rnorm(p), R = diag(p), lambda = 0.3,
    maf = maf, maf_thresh = 0.05,
    prior_weights = prior_weights
  )

  expect_equal(result$data$p, n_keep)
  expect_length(result$params$prior_weights, n_keep)
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-8)
})

test_that("rss_lambda_constructor null_weight R-path adds null column", {
  set.seed(85)
  p <- 15

  result <- rss_lambda_constructor(z = rnorm(p), R = diag(p), lambda = 0.2,
                                   null_weight = 0.15)

  expect_equal(result$data$p, p + 1L)
  expect_equal(result$data$z[p + 1L], 0)
  expect_equal(result$params$null_weight, 0.15)
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-8)
})

test_that("rss_lambda_constructor estimate_residual_variance=TRUE is stored", {
  set.seed(86)
  result <- rss_lambda_constructor(z = rnorm(15), R = diag(15), lambda = 0.2,
                                   estimate_residual_variance = TRUE)

  expect_true(result$params$estimate_residual_variance)
})

test_that("rss_lambda_constructor maf+prior_weights+null_weight chain preserves normalization", {
  set.seed(87)
  p <- 12
  maf <- c(rep(0.01, 4), rep(0.4, 8))
  n_keep <- sum(maf > 0.05)

  result <- rss_lambda_constructor(
    z = rnorm(p), R = diag(p), lambda = 0.1,
    prior_weights = rep(2, p),
    maf = maf, maf_thresh = 0.05,
    null_weight = 0.1
  )

  expect_equal(result$data$p, n_keep + 1L)
  expect_equal(sum(result$params$prior_weights), 1, tolerance = 1e-8)
})

# =============================================================================
# SUMMARY STATISTICS CONSTRUCTOR
# =============================================================================

test_that("summary_stats_constructor routes to sufficient_stats data class", {
  set.seed(90)
  p <- 50; z <- rnorm(p); R <- diag(p)

  result <- summary_stats_constructor(z = z, R = R, n = 100)

  expect_s3_class(result$data, "ss")
  expect_true(all(c("XtX", "Xty", "yty") %in% names(result$data)))
})

# --- Input validation ---

test_that("summary_stats_constructor rejects R with wrong dimensions", {
  expect_error(
    summary_stats_constructor(z = rnorm(50), R = diag(40), n = 100),
    "The dimension of R \\(40 x 40\\) does not agree with expected \\(50 x 50\\)"
  )
})

test_that("summary_stats_constructor rejects n <= 1", {
  z <- rnorm(50); R <- diag(50)

  for (bad_n in c(1, 0, -5)) {
    expect_error(
      summary_stats_constructor(z = z, R = R, n = bad_n),
      "n must be greater than 1"
    )
  }
})

test_that("summary_stats_constructor rejects mismatched bhat/shat lengths", {
  expect_error(
    summary_stats_constructor(bhat = rnorm(50), shat = abs(rnorm(45)),
                              R = diag(50), n = 100),
    "The lengths of bhat and shat do not agree"
  )
})

test_that("summary_stats_constructor accepts scalar shat and replicates it", {
  set.seed(91)
  result <- summary_stats_constructor(bhat = rnorm(50), shat = 0.1,
                                      R = diag(50), n = 100)
  expect_s3_class(result$data, "ss")
})

test_that("summary_stats_constructor rejects missing values in bhat or shat", {
  p <- 50; R <- diag(p)
  bhat <- rnorm(p); shat <- abs(rnorm(p))

  bhat_na <- bhat; bhat_na[5] <- NA
  expect_error(
    summary_stats_constructor(bhat = bhat_na, shat = shat, R = R, n = 100),
    "bhat, shat cannot have missing values"
  )

  shat_na <- shat; shat_na[10] <- NA
  expect_error(
    summary_stats_constructor(bhat = bhat, shat = shat_na, R = R, n = 100),
    "bhat, shat cannot have missing values"
  )
})

test_that("summary_stats_constructor rejects zero or negative elements in shat", {
  p <- 50; R <- diag(p); bhat <- rnorm(p); shat <- abs(rnorm(p))

  shat_zero <- shat; shat_zero[5] <- 0
  expect_error(
    summary_stats_constructor(bhat = bhat, shat = shat_zero, R = R, n = 100),
    "shat cannot have zero or negative elements"
  )

  shat_neg <- shat; shat_neg[8] <- -0.5
  expect_error(
    summary_stats_constructor(bhat = bhat, shat = shat_neg, R = R, n = 100),
    "shat cannot have zero or negative elements"
  )
})

test_that("summary_stats_constructor rejects empty z vector", {
  expect_error(
    summary_stats_constructor(z = numeric(0), R = matrix(0, 0, 0), n = 100),
    "Input vector z should have at least one element"
  )
})

test_that("summary_stats_constructor rejects MAF with wrong length", {
  expect_error(
    summary_stats_constructor(z = rnorm(50), R = diag(50), n = 100,
                              maf = runif(40)),
    "The length of maf does not agree with expected 50"
  )
})

test_that("summary_stats_constructor validates R-mismatch diagnostic parameters", {
  set.seed(92)
  p <- 30; z <- rnorm(p); R <- diag(p)

  expect_error(summary_stats_constructor(z = z, R = R, n = 100, eig_delta_rel = -1),
               "eig_delta_rel")
  expect_error(summary_stats_constructor(z = z, R = R, n = 100,
                                         eig_delta_rel = c(0.1, 0.2)),
               "eig_delta_rel")
  expect_error(summary_stats_constructor(z = z, R = R, n = 100, eig_delta_abs = -0.5),
               "eig_delta_abs")
  expect_error(summary_stats_constructor(z = z, R = R, n = 100, artifact_threshold = -0.1),
               "artifact_threshold")
  expect_error(summary_stats_constructor(z = z, R = R, n = 100, artifact_threshold = 1.1),
               "artifact_threshold")
  expect_error(summary_stats_constructor(z = z, R = R, n = 100, R_sensitivity_threshold = -1),
               "R_sensitivity_threshold")
  expect_error(summary_stats_constructor(z = z, R = R, n = 100, R_sensitivity_threshold = Inf),
               "R_sensitivity_threshold")
})

test_that("summary_stats_constructor rejects when both R and X provided", {
  set.seed(93)
  p <- 30
  expect_error(
    summary_stats_constructor(z = rnorm(p), R = diag(p),
                              X = matrix(rnorm(100 * p), 100, p), n = 100),
    "Please provide either R or X, but not both"
  )
})

test_that("summary_stats_constructor rejects when neither R nor X provided", {
  expect_error(
    summary_stats_constructor(z = rnorm(30), n = 100),
    "Please provide either R \\(correlation matrix\\) or X \\(factor matrix\\)"
  )
})

test_that("summary_stats_constructor requires either z or bhat/shat but not both", {
  R <- diag(50)
  expect_error(
    summary_stats_constructor(R = R, n = 100),
    "Please provide either z or \\(bhat, shat\\)"
  )
  expect_error(
    summary_stats_constructor(z = rnorm(50), bhat = rnorm(50), shat = runif(50, 0.5, 1.5),
                              R = R, n = 100),
    "Please provide either z or \\(bhat, shat\\), but not both"
  )
})

test_that("summary_stats_constructor: NIG residual method requires valid n", {
  set.seed(94)
  p <- 8; z <- rnorm(p); R <- diag(p)

  expect_error(
    summary_stats_constructor(z = z, R = R, estimate_residual_method = "NIG", n = NULL),
    "requires a valid sample"
  )
})

test_that("summary_stats_constructor: R_finite = TRUE requires X", {
  set.seed(95)
  p <- 8

  expect_error(
    summary_stats_constructor(z = rnorm(p), R = diag(p), n = 100, R_finite = TRUE),
    "R_finite = TRUE requires X"
  )
})

test_that("summary_stats_constructor: single-panel X must be numeric matrix", {
  set.seed(96)
  p <- 8

  X_char <- matrix(letters[1:(3 * p)], nrow = 3, ncol = p)
  expect_error(
    summary_stats_constructor(z = rnorm(p), X = X_char, n = 100),
    "X must be a numeric matrix."
  )
})

test_that("summary_stats_constructor errors when ncol(X) != length(z) for low-rank X", {
  set.seed(97)
  p_z <- 15; p_X <- 12; n_X <- 8
  X <- matrix(rnorm(n_X * p_X), n_X, p_X)

  expect_error(
    summary_stats_constructor(z = rnorm(p_z), X = X, n = 200),
    "number of columns of X"
  )
})

# --- Lambda=0 / sufficient stats path ---

test_that("summary_stats_constructor handles z without n with informational message", {
  set.seed(98)
  p <- 50; z <- rnorm(p); R <- diag(p)

  expect_message(
    result <- summary_stats_constructor(z = z, R = R),
    "Providing the sample size"
  )
  expect_s3_class(result$data, "ss")
})

test_that("summary_stats_constructor converts bhat/shat to z then to ss", {
  set.seed(99)
  p <- 50

  result <- summary_stats_constructor(bhat = rnorm(p), shat = runif(p, 0.5, 1.5),
                                      R = diag(p), n = 100)
  expect_s3_class(result$data, "ss")
})

test_that("summary_stats_constructor handles shat and var_y for original scale effects", {
  set.seed(100)
  p <- 50
  bhat <- rnorm(p)
  shat <- abs(rnorm(p, mean = 0.1, sd = 0.02))
  var_y <- 2.5; n <- 100

  result <- summary_stats_constructor(bhat = bhat, shat = shat, var_y = var_y,
                                      R = diag(p), n = n)

  expect_s3_class(result$data, "ss")
  expect_equal(result$data$yty, (n - 1) * var_y)
})

test_that("summary_stats_constructor applies MAF filter to bhat/shat", {
  set.seed(101)
  p <- 20
  bhat <- rnorm(p); shat <- abs(rnorm(p, mean = 0.3, sd = 0.05))
  maf <- c(rep(0.01, 10), rep(0.3, 10))

  result <- summary_stats_constructor(bhat = bhat, shat = shat, R = diag(p), n = 300,
                                      maf = maf, maf_thresh = 0.05)

  expect_equal(result$data$p, sum(maf > 0.05))
  expect_s3_class(result$data, "ss")
})

test_that("summary_stats_constructor applies MAF filter to low-rank X", {
  set.seed(102)
  p <- 20; n_rows <- 8
  X <- matrix(rnorm(n_rows * p), n_rows, p)
  maf <- c(rep(0.01, 5), rep(0.30, 15))
  n_keep <- sum(maf > 0.05)

  result <- summary_stats_constructor(z = rnorm(p), X = X, n = 200,
                                      maf = maf, maf_thresh = 0.05)
  expect_equal(result$data$p, n_keep)
})

test_that("summary_stats_constructor warns for ash + z-only input (standardized scale)", {
  set.seed(103)
  p <- 15

  expect_message(
    result <- summary_stats_constructor(z = rnorm(p), R = diag(p), n = 200,
                                        unmappable_effects = "ash"),
    "standardized scale"
  )
  expect_s3_class(result$data, "ss")
})

# --- R_mismatch branches ---

test_that("summary_stats_constructor disables sigma^2 estimation when R_mismatch is active", {
  set.seed(104)
  p <- 30

  expect_message(
    result <- summary_stats_constructor(
      z = rnorm(p), R = diag(p), n = 500,
      R_mismatch = "eb_no_init",
      estimate_residual_variance = TRUE,
      verbose = FALSE
    ),
    "incompatible with"
  )
  expect_false(result$params$estimate_residual_variance)
})

# --- X-based paths ---

test_that("summary_stats_constructor converts full-rank X to R", {
  set.seed(105)
  n <- 100; p <- 30
  X <- scale(matrix(rnorm(n * p), n, p), TRUE, TRUE)

  result <- summary_stats_constructor(z = rnorm(p), X = X, n = n)
  expect_s3_class(result$data, "ss")
})

test_that("summary_stats_constructor keeps low-rank X without conversion", {
  set.seed(106)
  n <- 30; p <- 80
  X <- scale(matrix(rnorm(n * p), n, p), TRUE, TRUE)

  result <- summary_stats_constructor(z = rnorm(p), X = X, n = n)
  expect_type(result, "list")
})

test_that("summary_stats_constructor emits warning for low-rank X with var_y/shat", {
  set.seed(107)
  n <- 20; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  bhat <- rnorm(p)
  shat <- abs(rnorm(p, mean = 0.3, sd = 0.05))

  expect_message(
    result <- summary_stats_constructor(bhat = bhat, shat = shat, var_y = 2.0,
                                        X = X, n = n),
    "low-rank factor matrix"
  )
  expect_s3_class(result$data, "ss")
})

# --- Multi-panel ---

test_that("summary_stats_constructor handles multi-panel R input", {
  set.seed(108)
  p <- 25; z <- rnorm(p)
  Ra <- cor(matrix(rnorm(500 * p), 500, p))
  Rb <- cor(matrix(rnorm(500 * p), 500, p))

  result <- summary_stats_constructor(
    z = z, R = list(Ra, Rb), n = 500, L = 5, max_iter = 20, verbose = FALSE
  )
  expect_type(result, "list")
  expect_length(result$multi_panel_meta$fits, 2)
  expect_length(result$multi_panel_meta$elbos, 2)
})

test_that("summary_stats_constructor multi-panel rejects bhat/shat and var_y", {
  set.seed(109)
  p <- 6; z <- rnorm(p); Rlist <- list(diag(p), diag(p))

  expect_error(
    summary_stats_constructor(z = z, R = Rlist, n = 100,
                              bhat = rnorm(p), shat = rep(1, p)),
    "multi-panel"
  )
  expect_error(
    summary_stats_constructor(z = z, R = Rlist, n = 100, var_y = 1),
    "multi-panel"
  )
})

test_that("summary_stats_constructor multi-panel requires n", {
  set.seed(110)
  p <- 6; z <- rnorm(p); Rlist <- list(diag(p), diag(p))

  expect_error(
    summary_stats_constructor(z = z, R = Rlist, n = NULL),
    "Sample size 'n' is required"
  )
})

test_that("summary_stats_constructor multi-panel validates R list elements", {
  set.seed(111)
  p <- 6; z <- rnorm(p)

  expect_error(
    summary_stats_constructor(z = z, R = list(diag(p), matrix(1, p, p - 1)), n = 100),
    "must be square"
  )
  expect_error(
    summary_stats_constructor(z = z, R = list(diag(p), "notamatrix"), n = 100),
    "numeric matrix"
  )
  expect_error(
    summary_stats_constructor(z = z, X = list(diag(p),
                                               as.data.frame(matrix(rnorm(100 * p), 100, p))),
                              n = 100),
    "Each element of X list must be a numeric matrix."
  )
})

# =============================================================================
# SS MIXTURE CONSTRUCTOR
# =============================================================================

test_that("ss_mixture_constructor runs 2-panel fit and returns ss_mixture class", {
  set.seed(120)
  p <- 15; n <- 300
  z <- rnorm(p)
  R1 <- cor(matrix(rnorm(n * p), n, p))
  R2 <- cor(matrix(rnorm(n * p), n, p))

  result <- ss_mixture_constructor(z = z, R = list(R1, R2), n = n, L = 3,
                                   max_iter = 5, verbose = FALSE)

  expect_s3_class(result$data, "ss_mixture")
  expect_true("omega_init" %in% names(result$data))
  expect_equal(result$data$K, 2L)
  expect_length(result$data$panel_R, 2L)
})

test_that("ss_mixture_constructor accepts explicit init_panel parameter", {
  set.seed(121)
  p <- 25; n <- 500
  z <- rnorm(p)
  Ra <- cor(matrix(rnorm(n * p), n, p))
  Rb <- cor(matrix(rnorm(n * p), n, p))

  res_default  <- ss_mixture_constructor(z = z, R = list(Ra, Rb), n = n, L = 5,
                                         max_iter = 5, verbose = FALSE)
  res_explicit <- ss_mixture_constructor(z = z, R = list(Ra, Rb), n = n, L = 5,
                                         max_iter = 5, verbose = FALSE, init_panel = 2)

  expect_s3_class(res_default$data,  "ss_mixture")
  expect_s3_class(res_explicit$data, "ss_mixture")
})

test_that("ss_mixture_constructor applies maf filter + null_weight + estimate_residual_variance", {
  set.seed(122)
  p <- 20; n <- 500
  z <- rnorm(p)
  Ra <- cor(matrix(rnorm(n * p), n, p))
  Rb <- cor(matrix(rnorm(n * p), n, p))
  maf <- c(rep(0.01, 5), rep(0.25, 15))
  n_keep <- sum(maf > 0.05)

  result <- ss_mixture_constructor(
    z = z, R = list(Ra, Rb), n = n, L = 3,
    maf = maf, maf_thresh = 0.05,
    null_weight = 0.1, estimate_residual_variance = TRUE,
    max_iter = 5, verbose = FALSE
  )

  expect_equal(result$data$p, n_keep + 1L)
  expect_equal(result$params$null_weight, 0.1)
  expect_true(result$params$estimate_residual_variance)
})

# --- Input validation ---

test_that("ss_mixture_constructor rejects invalid constructor inputs", {
  set.seed(123)
  p <- 10; n <- 300; z <- rnorm(p)
  R_ok <- list(diag(p), diag(p))

  expect_error(ss_mixture_constructor(z = z, R = R_ok, n = NULL),
               "Sample size 'n' is required for multi-panel mode.")
  expect_error(ss_mixture_constructor(z = z, R = R_ok, n = n,
                                      R_sensitivity_threshold = -1),
               "R_sensitivity_threshold must be a single nonnegative finite numeric.")
  expect_error(ss_mixture_constructor(z = NULL, R = R_ok, n = n),
               "Multi-panel mode requires z-scores.")
  expect_error(ss_mixture_constructor(z = z, R = R_ok,
                                      X = list(matrix(rnorm(n * p), n, p),
                                               matrix(rnorm(n * p), n, p)),
                                      n = n),
               "Please provide either R or X, but not both.")
  expect_error(ss_mixture_constructor(z = z, R = list(), n = n),
               "Multi-panel input must contain at least one panel.")
})

test_that("ss_mixture_constructor validates R panel types and dimensions", {
  set.seed(124)
  p <- 6; n <- 200; z <- rnorm(p)

  R_char <- matrix(letters[seq_len(p * p)], p, p)
  expect_error(ss_mixture_constructor(z = z, R = list(R_char), n = n),
               "Each element of R list must be a numeric matrix.")

  R_wrong_dim <- matrix(rnorm(81), 9, 9)
  R_wrong_dim <- cov2cor(R_wrong_dim %*% t(R_wrong_dim) + diag(9))
  expect_error(ss_mixture_constructor(z = z, R = list(R_wrong_dim), n = n),
               "Each element of R list must have dimension length\\(z\\) by length\\(z\\).")
})

test_that("ss_mixture_constructor validates X panel types and dimensions", {
  set.seed(125)
  p <- 6; n <- 100; z <- rnorm(p)

  X_char <- matrix(as.character(rnorm(n * p)), n, p)
  expect_error(ss_mixture_constructor(z = z, X = list(X_char), n = n),
               "Each element of X list must be a numeric matrix.")

  X_wrong_ncol <- matrix(rnorm(n * 8), n, 8)
  expect_error(ss_mixture_constructor(z = z, X = list(X_wrong_ncol), n = n),
               "Each element of X list must have length\\(z\\) columns.")
})

test_that("ss_mixture_constructor rejects MAF length mismatch", {
  set.seed(126)
  p <- 10; n <- 300; z <- rnorm(p)

  expect_error(
    ss_mixture_constructor(z = z, R = list(diag(p), diag(p)), n = n,
                           maf = c(0.1, 0.2)),
    "The length of maf does not agree with expected"
  )
})

test_that("ss_mixture_constructor rejects infinite z and replaces NA z with zero", {
  set.seed(127)
  p <- 8; n <- 200
  R <- list(diag(p), diag(p))

  z_inf <- rnorm(p); z_inf[2] <- Inf
  expect_error(ss_mixture_constructor(z = z_inf, R = R, n = n),
               "z contains infinite values.")

  z_na <- rnorm(p); z_na[3] <- NA
  expect_message(
    result <- ss_mixture_constructor(z = z_na, R = R, n = n, L = 2, max_iter = 2),
    "z-score"
  )
  expect_s3_class(result$data, "ss_mixture")
})

test_that("ss_mixture_constructor X-list path runs and applies maf + null_weight", {
  set.seed(128)
  p <- 12; n <- 200; z <- rnorm(p)
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p), n, p)

  result_x <- ss_mixture_constructor(z = z, X = list(X1, X2), n = n, L = 2,
                                     max_iter = 2)
  expect_s3_class(result_x$data, "ss_mixture")
  expect_equal(result_x$data$K, 2L)

  maf <- c(rep(0.01, 4), rep(0.25, 8))
  pw  <- rep(1 / p, p)
  n_keep <- sum(maf > 0.05)

  result_maf <- ss_mixture_constructor(
    z = z, X = list(X1, X2), n = n, L = 2,
    maf = maf, maf_thresh = 0.05, prior_weights = pw, max_iter = 2
  )
  expect_equal(result_maf$data$p, n_keep)
  expect_length(result_maf$params$prior_weights, n_keep)

  result_nw <- ss_mixture_constructor(z = z, X = list(X1, X2), n = n, L = 2,
                                      null_weight = 0.1, max_iter = 2)
  expect_equal(result_nw$data$p, p + 1L)
  expect_equal(result_nw$params$null_weight, 0.1)
})

test_that("ss_mixture_constructor R_mismatch != none sets R_finite_B = Inf", {
  set.seed(129)
  p <- 8; n <- 300
  R1 <- cor(matrix(rnorm(n * p), n, p))
  R2 <- cor(matrix(rnorm(n * p), n, p))

  result <- ss_mixture_constructor(z = rnorm(p), R = list(R1, R2), n = n, L = 2,
                                   R_mismatch = "eb", max_iter = 2)
  expect_s3_class(result$data, "ss_mixture")
  expect_equal(result$data$R_finite_B, Inf)
})

# =============================================================================
# CONSTRUCTOR HELPER FUNCTIONS
# =============================================================================

test_that("resolve_model_init validates mutual exclusivity and deprecation", {
  expect_error(
    resolve_model_init(model_init = list(a = 1), s_init = list(b = 2)),
    "Cannot specify both"
  )
  expect_message(
    res <- resolve_model_init(model_init = NULL, s_init = list(b = 2)),
    "s_init is deprecated"
  )
  expect_equal(res, list(b = 2))
  expect_equal(resolve_model_init(model_init = list(a = 1), s_init = NULL),
               list(a = 1))
})

test_that("normalize_null_weight validates type, range, and zero case", {
  expect_error(normalize_null_weight("x", NULL, 5), "Null weight must be numeric")
  expect_error(normalize_null_weight(c(0.1, 0.2), NULL, 5), "Null weight must be numeric")
  expect_error(normalize_null_weight(NA_real_, NULL, 5), "Null weight must be numeric")
  expect_error(normalize_null_weight(-0.1, NULL, 5), "between 0 and 1")
  expect_error(normalize_null_weight(1, NULL, 5), "between 0 and 1")

  res_zero <- normalize_null_weight(0, NULL, 5)
  expect_null(res_zero$null_weight)
  expect_false(res_zero$add_null)
})

test_that("normalize_null_weight builds prior_weights with and without supplied weights", {
  res1 <- normalize_null_weight(0.1, NULL, 4)
  expect_true(res1$add_null)
  expect_length(res1$prior_weights, 5)
  expect_equal(res1$prior_weights[5], 0.1)
  expect_equal(res1$prior_weights[1:4], rep((1 - 0.1) / 4, 4))

  pw <- c(0.4, 0.3, 0.2, 0.1)
  res2 <- normalize_null_weight(0.2, pw, 4)
  expect_equal(res2$prior_weights, c(pw * 0.8, 0.2))
})

test_that("normalize_prior_weights validates length, positivity, and NULL default", {
  expect_error(normalize_prior_weights(rep(1, 4), 5), "length p")
  expect_error(normalize_prior_weights(rep(0, 5), 5), "greater than 0")
  expect_equal(normalize_prior_weights(NULL, 5), rep(1 / 5, 5))
})

test_that("normalize_summary_stats_input validates input exclusivity and types", {
  expect_error(normalize_summary_stats_input(z = rnorm(5), bhat = rnorm(5)), "not both")
  expect_error(normalize_summary_stats_input(z = NULL, bhat = rnorm(5), shat = NULL),
               "either z or")
  expect_error(normalize_summary_stats_input(bhat = letters[1:3], shat = rep(1, 3)),
               "must be numeric")
  expect_error(normalize_summary_stats_input(bhat = rnorm(4), shat = rep(1, 3)),
               "do not agree")
  expect_error(normalize_summary_stats_input(bhat = c(1, NA, 3), shat = rep(1, 3)),
               "missing values")
  expect_error(normalize_summary_stats_input(bhat = rnorm(3), shat = c(1, 0, 1)),
               "zero or negative")
  expect_error(normalize_summary_stats_input(z = letters[1:3]), "z must be numeric")
  expect_error(normalize_summary_stats_input(z = numeric(0)), "at least one element")
})

test_that("normalize_summary_stats_input recycles scalar shat and carries names", {
  res_scalar <- normalize_summary_stats_input(bhat = c(2, 4, 6), shat = 2)
  expect_equal(res_scalar$z, c(1, 2, 3))

  bhat_named <- setNames(rnorm(8), paste0("SNP", seq_len(8)))
  shat_named <- abs(rnorm(8, mean = 0.5, sd = 0.1))
  res_names <- normalize_summary_stats_input(bhat = bhat_named, shat = shat_named)
  expect_equal(names(res_names$z), names(bhat_named))
  expect_equal(unname(res_names$z), as.numeric(bhat_named) / as.numeric(shat_named),
               tolerance = 1e-8)
})

test_that("summary_stats_working_quantities errors when pve_adjustment missing", {
  expect_error(
    summary_stats_working_quantities(z = rnorm(5), n = 100, shat = rep(0.1, 5),
                                     var_y = 1, pve_adjustment = NULL),
    "PVE adjustment is required"
  )
})

# =============================================================================
# INTEGRATION - Constructor output usability
# =============================================================================

test_that("individual_data_constructor output works with ibss_initialize", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 200)
  result <- individual_data_constructor(base_data$X, base_data$y, L = 5)

  expect_no_error(ibss_initialize(result$data, result$params))
})

test_that("sufficient_stats_constructor output works with ibss_initialize", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 201)
  ss <- compute_summary_stats(base_data$X, base_data$y)

  result <- sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n,
                                         XtX = ss$XtX, L = 5)

  expect_no_error(ibss_initialize(result$data, result$params))
})

test_that("rss_lambda_constructor output works with ibss_initialize", {
  set.seed(202)
  result <- rss_lambda_constructor(rnorm(50), diag(50), lambda = 0.5, L = 5)

  expect_no_error(ibss_initialize(result$data, result$params))
})

test_that("sufficient_stats_constructor maf filtering subsets prior_weights consistently", {
  p <- 6; n <- 100
  X <- matrix(rnorm(n * p), n, p)
  ss <- compute_summary_stats(X, rnorm(n))
  prior_weights <- c(0.10, 0.20, 0.30, 0.25, 0.10, 0.05)
  maf  <- c(0.20, 0.02, 0.30, 0.25, 0.01, 0.40)
  keep <- maf > 0.05

  result <- sufficient_stats_constructor(
    Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX, L = 1,
    prior_weights = prior_weights, maf = maf, maf_thresh = 0.05
  )

  expect_equal(result$params$prior_weights,
               prior_weights[keep] / sum(prior_weights[keep]))
})

test_that("ss_mixture_constructor maf filtering subsets prior_weights consistently", {
  p <- 6; n <- 100
  z <- rnorm(p); R <- list(diag(p), diag(p))
  prior_weights <- c(0.10, 0.20, 0.30, 0.25, 0.10, 0.05)
  maf  <- c(0.20, 0.02, 0.30, 0.25, 0.01, 0.40)
  keep <- maf > 0.05

  result <- ss_mixture_constructor(z = z, R = R, n = n, L = 1,
                                   prior_weights = prior_weights,
                                   maf = maf, maf_thresh = 0.05)

  expect_equal(result$params$prior_weights,
               prior_weights[keep] / sum(prior_weights[keep]))
})

test_that("rss_lambda_constructor maf filtering subsets prior_weights consistently", {
  p <- 6
  z <- rnorm(p); R <- diag(p)
  prior_weights <- c(0.10, 0.20, 0.30, 0.25, 0.10, 0.05)
  maf  <- c(0.20, 0.02, 0.30, 0.25, 0.01, 0.40)
  keep <- maf > 0.05

  result <- rss_lambda_constructor(z = z, R = R, lambda = 0.5, L = 1,
                                   prior_weights = prior_weights,
                                   maf = maf, maf_thresh = 0.05)

  expect_equal(result$params$prior_weights,
               prior_weights[keep] / sum(prior_weights[keep]))
})

test_that("susie_rss with init_only = TRUE works on multi-panel input", {
  set.seed(203)
  p <- 25
  z <- rnorm(p)
  Ra <- cor(matrix(rnorm(500 * p), 500, p))
  Rb <- cor(matrix(rnorm(500 * p), 500, p))

  expect_no_error(
    res <- susie_rss(z = z, R = list(Ra, Rb), n = 500, L = 5,
                     init_only = TRUE, max_iter = 20, verbose = FALSE)
  )
  expect_length(res$multi_panel_meta$fits, 2)
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("individual_data_constructor handles p=1 single variable", {
  set.seed(300)
  n <- 50
  X <- matrix(rnorm(n), n, 1)
  y <- rnorm(n)

  result <- individual_data_constructor(X, y)

  expect_equal(result$data$p, 1)
  expect_length(result$params$prior_weights, 1)
  expect_equal(result$params$prior_weights, 1)
})

test_that("sufficient_stats_constructor handles p=1 single variable", {
  set.seed(301)
  n <- 50; p <- 1
  X <- matrix(rnorm(n), n, p)
  y <- rnorm(n)
  ss <- compute_summary_stats(X, y)

  result <- sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n,
                                         XtX = ss$XtX)
  expect_equal(result$data$p, 1)
})

test_that("sufficient_stats_constructor maf filter that removes all variables errors gracefully", {
  set.seed(302)
  p <- 5; n <- 50
  ss <- compute_summary_stats(matrix(rnorm(n * p), n, p), rnorm(n))
  maf_all_low <- rep(0.01, p)

  expect_error(
    sufficient_stats_constructor(Xty = ss$Xty, yty = ss$yty, n = ss$n, XtX = ss$XtX,
                                 maf = maf_all_low, maf_thresh = 0.05),
    regexp = "."
  )
})

test_that("individual_data_constructor null_weight at boundary 0 leaves p unchanged", {
  base_data <- generate_base_data(n = 50, p = 10, k = 0, seed = 303)
  result <- individual_data_constructor(base_data$X, base_data$y, null_weight = 0)

  expect_equal(result$data$p, 10)
  expect_null(result$params$null_weight)
  expect_length(result$params$prior_weights, 10)
})
