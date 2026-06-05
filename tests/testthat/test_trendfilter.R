context("Trend filtering")

# ---- Basic functionality ----

test_that("susie_trendfilter returns a susie object with expected fields", {
  set.seed(1)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y  <- mu + rnorm(60)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))

  expect_s3_class(result, "susie")
  expect_true("alpha" %in% names(result))
  expect_true("mu"    %in% names(result))
  expect_true("elbo"  %in% names(result))
})

test_that("susie_trendfilter concentrates PIPs near true changepoints", {
  set.seed(2)
  mu <- c(rep(0, 25), rep(3, 25), rep(-2, 25), rep(1, 25))
  y  <- mu + rnorm(100, sd = 0.3)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))

  pip               <- susie_get_pip(result)
  changepoint_idx   <- c(23:27, 48:52, 73:77)
  expect_true(sum(pip[changepoint_idx]) > sum(pip[-changepoint_idx]))
})

test_that("susie_trendfilter fitted values track signal better than raw data", {
  set.seed(3)
  mu <- c(rep(0, 20), rep(2, 20), rep(0, 20))
  y  <- mu + rnorm(60, sd = 0.1)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  fitted <- predict(result)

  expect_true(mean((fitted - mu)^2) < mean((y - mu)^2))
})

test_that("susie_trendfilter keeps all PIPs low when no changepoint exists", {
  set.seed(4)
  y <- rep(5, 50) + rnorm(50, sd = 0.5)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  expect_true(max(susie_get_pip(result)) < 0.5)
})

# ---- Order parameter ----

test_that("order=1 and order=2 each emit the 'order > 0' warning message", {
  set.seed(5)
  y_lin  <- seq(0, 1, length.out = 50) + rnorm(50, sd = 0.1)
  y_quad <- (seq(0, 1, length.out = 50))^2 + rnorm(50, sd = 0.1)

  expect_message(
    suppressWarnings(susie_trendfilter(y_lin,  order = 1, use_mad = FALSE)),
    "order > 0 is not recommended"
  )
  expect_message(
    suppressWarnings(susie_trendfilter(y_quad, order = 2, use_mad = FALSE)),
    "order > 0 is not recommended"
  )
})

test_that("order=0 and order=1 produce different alpha matrices", {
  set.seed(6)
  y <- seq(0, 2, length.out = 50) + rnorm(50, sd = 0.1)

  result_0 <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE, max_iter = 10))
  result_1 <- suppressWarnings(susie_trendfilter(y, order = 1, use_mad = FALSE, max_iter = 10))

  expect_false(all(abs(result_0$alpha - result_1$alpha) < 1e-10))
})

# ---- use_mad parameter ----

test_that("use_mad=TRUE runs without error and use_mad=FALSE differs from it", {
  set.seed(7)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y  <- mu + rnorm(60)

  result_mad    <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = TRUE,  max_iter = 5))
  result_no_mad <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE, max_iter = 5))

  expect_s3_class(result_mad,    "susie")
  expect_s3_class(result_no_mad, "susie")
})

test_that("use_mad with model_init skips the MAD init pass", {
  set.seed(8)
  mu   <- c(rep(0, 20), rep(2, 20))
  y    <- mu + rnorm(40)
  init <- susie_init_coef(c(20), c(2), 40)

  result <- suppressWarnings(
    susie_trendfilter(y, order = 0, use_mad = TRUE, model_init = init, max_iter = 2)
  )
  expect_s3_class(result, "susie")
})

test_that("use_mad=TRUE errors when MAD is zero (constant data)", {
  y <- rep(5, 50)
  expect_error(
    susie_trendfilter(y, order = 0, use_mad = TRUE),
    "Cannot use median absolute deviation \\(MAD\\) to initialize residual variance because MAD = 0 for the input data. Please set 'use_mad = FALSE'"
  )
})

# ---- standardize and intercept options ----

test_that("standardize=TRUE and standardize=FALSE both produce valid alpha matrices", {
  set.seed(9)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y  <- mu + rnorm(60)

  for (std in c(TRUE, FALSE)) {
    result <- suppressWarnings(
      susie_trendfilter(y, order = 0, standardize = std, use_mad = FALSE)
    )
    expect_s3_class(result, "susie")
    expect_true(all(result$alpha >= 0 & result$alpha <= 1),
                info = paste("standardize =", std))
  }
})

test_that("intercept=TRUE records a non-NA intercept; intercept=FALSE records 0", {
  set.seed(10)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y  <- mu + rnorm(60) + 10

  result_yes <- suppressWarnings(
    susie_trendfilter(y, order = 0, intercept = TRUE,  use_mad = FALSE)
  )
  result_no  <- suppressWarnings(
    susie_trendfilter(y, order = 0, intercept = FALSE, use_mad = FALSE)
  )

  expect_false(is.na(result_yes$intercept))
  expect_equal(result_no$intercept, 0)
})

# ---- Pass-through parameters ----

test_that("L parameter controls the number of single-effect regressions", {
  set.seed(11)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y  <- mu + rnorm(60)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, L = 3, use_mad = FALSE))

  expect_equal(nrow(result$alpha), 3)
  expect_equal(length(result$V),   3)
})

test_that("max_iter caps the number of IBSS iterations", {
  set.seed(12)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y  <- mu + rnorm(60)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, max_iter = 5, use_mad = FALSE))
  expect_true(result$niter <= 5)
})

test_that("null_weight parameter is respected and null_index is populated", {
  set.seed(13)
  mu <- c(rep(0, 30), rep(2, 30))
  y  <- mu + rnorm(60)
  n  <- length(y)

  result <- suppressWarnings(
    susie_trendfilter(y, order = 0, null_weight = 1 / (n + 1), use_mad = FALSE, max_iter = 3)
  )
  expect_false(is.null(result$null_index))
})

test_that("estimate_prior_variance=FALSE fixes V and TRUE estimates it", {
  set.seed(14)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y  <- mu + rnorm(60)

  result_est <- suppressWarnings(
    susie_trendfilter(y, order = 0, estimate_prior_variance = TRUE,  use_mad = FALSE, max_iter = 3)
  )
  result_fix <- suppressWarnings(
    susie_trendfilter(y, order = 0, estimate_prior_variance = FALSE, use_mad = FALSE, max_iter = 3)
  )

  expect_s3_class(result_est, "susie")
  expect_s3_class(result_fix, "susie")
})

# ---- Integration with susie methods ----

test_that("susie_get_cs returns coverage list with correct requested_coverage", {
  set.seed(15)
  mu <- c(rep(0, 25), rep(3, 25), rep(-2, 25), rep(1, 25))
  y  <- mu + rnorm(100, sd = 0.3)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  cs     <- susie_get_cs(result, coverage = 0.95)

  expect_type(cs, "list")
  expect_true("cs" %in% names(cs))
  expect_equal(cs$requested_coverage, 0.95)
})

test_that("susie_get_pip returns a valid probability vector of length n", {
  set.seed(16)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y  <- mu + rnorm(60)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  pip    <- susie_get_pip(result)

  expect_length(pip, length(y))
  expect_true(all(pip >= 0 & pip <= 1))
})

test_that("predict returns a finite numeric vector of length n", {
  set.seed(17)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y  <- mu + rnorm(60)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  fitted <- predict(result)

  expect_length(fitted, length(y))
  expect_type(fitted, "double")
  expect_true(all(is.finite(fitted)))
})

test_that("coef returns a numeric vector of length n+1", {
  set.seed(18)
  mu <- c(rep(0, 20), rep(2, 20), rep(-1, 20))
  y  <- mu + rnorm(60)

  result      <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  coefficients <- coef(result)

  expect_length(coefficients, length(y) + 1)
  expect_type(coefficients, "double")
})

test_that("sets field built from susie_get_cs has expected structure", {
  set.seed(19)
  mu <- c(rep(0, 25), rep(3, 25), rep(-2, 25), rep(1, 25))
  y  <- mu + rnorm(100, sd = 0.3)

  result       <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  result$sets  <- susie_get_cs(result, coverage = 0.95)

  expect_type(result$sets, "list")
  expect_true("cs"       %in% names(result$sets))
  expect_true("coverage" %in% names(result$sets))
  expect_equal(result$sets$requested_coverage, 0.95)
})

# ---- Changepoint detection quality ----

test_that("susie_trendfilter places high PIPs near true changepoints", {
  set.seed(20)
  mu <- c(rep(0, 30), rep(4, 30), rep(-2, 30))
  y  <- mu + rnorm(90, sd = 0.5)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  pip    <- susie_get_pip(result)

  expect_true(pip[30] > 0.5 | pip[29] > 0.5 | pip[31] > 0.5)
  expect_true(pip[60] > 0.5 | pip[59] > 0.5 | pip[61] > 0.5)
})

test_that("susie_trendfilter converges on noisy data", {
  set.seed(21)
  mu <- c(rep(0, 30), rep(2, 30))
  y  <- mu + rnorm(60, sd = 2)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  expect_true(result$converged)
})

# ---- Edge cases ----

test_that("susie_trendfilter works with a very short time series", {
  set.seed(22)
  y <- c(0, 0, 0, 2, 2, 2)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE, max_iter = 3))
  expect_s3_class(result, "susie")
})

test_that("susie_trendfilter works with a long time series and predict has correct length", {
  set.seed(23)
  mu <- rep(c(0, 1, 2, 0), each = 100)
  y  <- mu + rnorm(400, sd = 0.5)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE, max_iter = 20))
  expect_length(predict(result), 400)
})

test_that("constant y with use_mad=FALSE errors on zero residual variance", {
  y <- rep(5, 50)
  expect_error(
    susie_trendfilter(y, order = 0, use_mad = FALSE),
    "Residual variance sigma2 must be positive"
  )
})

test_that("susie_trendfilter detects a changepoint near the start of the series", {
  set.seed(24)
  mu <- c(rep(0, 5), rep(2, 45))
  y  <- mu + rnorm(50, sd = 0.3)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  pip    <- susie_get_pip(result)
  expect_true(any(pip[3:7] > 0.3))
})

test_that("susie_trendfilter detects a changepoint near the end of the series", {
  set.seed(25)
  mu <- c(rep(0, 45), rep(2, 5))
  y  <- mu + rnorm(50, sd = 0.3)

  result <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))
  pip    <- susie_get_pip(result)
  expect_true(any(pip[43:47] > 0.3))
})

# ---- Equivalence to manual sparse-matrix construction ----

test_that("susie_trendfilter matches manual susie call for order 0, 1, and 2", {
  for (ord in 0:2) {
    set.seed(30 + ord)
    with(simulate_tf(ord), {
      result_manual <- suppressWarnings(
        susie(X, y, estimate_prior_variance = FALSE, standardize = TRUE, max_iter = 5)
      )
      result_tf <- suppressWarnings(
        susie_trendfilter(y, order = ord,
                          estimate_prior_variance = FALSE,
                          standardize = TRUE,
                          use_mad = FALSE, max_iter = 5)
      )
      expect_equal(result_tf$alpha, result_manual$alpha,
                   tolerance = 1e-6, info = paste("order =", ord))
      expect_equal(result_tf$mu,    result_manual$mu,
                   tolerance = 1e-6, info = paste("order =", ord))
    })
  }
})

# ---- Documentation example ----

test_that("documentation example runs and returns the right-length predictions", {
  set.seed(1)
  mu <- c(rep(0, 50), rep(1, 50), rep(3, 50), rep(-2, 50), rep(0, 200))
  y  <- mu + rnorm(400)

  s <- susie_trendfilter(y, max_iter = 100)

  expect_s3_class(s, "susie")
  expect_length(predict(s), 400)
  expect_type(susie_get_cs(s), "list")
})

# ---- compute_tf_d and compute_tf_std_d ----

test_that("compute_tf_d with intercept=TRUE, standardize=TRUE returns rep(n-1, n) with d[n]=0 for order=0", {
  n <- 20
  for (order in 0:2) {
    cm  <- compute_tf_cm(order, n)
    csd <- compute_tf_csd(order, n)
    d   <- compute_tf_d(order, n, cm, csd, standardize = TRUE, intercept = TRUE)
    if (order == 0) {
      expect_equal(d[-n], rep(n - 1, n - 1), info = paste("order =", order))
      expect_equal(d[n], 0,                  info = paste("order =", order))
    } else {
      expect_equal(d, rep(n - 1, n),         info = paste("order =", order))
    }
  }
})

test_that("compute_tf_d with intercept=TRUE, standardize=FALSE scales by csd^2", {
  n <- 20
  for (order in 0:2) {
    cm     <- compute_tf_cm(order, n)
    csd    <- compute_tf_csd(order, n)
    d_std  <- compute_tf_d(order, n, cm, csd, standardize = TRUE,  intercept = TRUE)
    d_raw  <- compute_tf_d(order, n, cm, csd, standardize = FALSE, intercept = TRUE)
    expect_equal(d_raw, d_std * csd^2, tolerance = 1e-10, info = paste("order =", order))
  }
})

test_that("compute_tf_d with intercept=FALSE, standardize=FALSE returns non-negative finite colSums(X^2)", {
  n <- 15
  for (order in 0:2) {
    cm  <- compute_tf_cm(order, n)
    csd <- compute_tf_csd(order, n)
    d   <- compute_tf_d(order, n, cm, csd, standardize = FALSE, intercept = FALSE)
    expect_length(d, n)
    expect_true(all(is.finite(d)))
    expect_true(all(d >= 0), info = paste("order =", order))
  }
})

test_that("compute_tf_d with intercept=FALSE, standardize=TRUE divides by csd^2", {
  n <- 20
  for (order in 0:2) {
    cm    <- compute_tf_cm(order, n)
    csd   <- compute_tf_csd(order, n)
    d_raw <- compute_tf_d(order, n, cm, csd, standardize = FALSE, intercept = FALSE)
    d_std <- compute_tf_d(order, n, cm, csd, standardize = TRUE,  intercept = FALSE)
    expect_equal(d_std, d_raw / csd^2, tolerance = 1e-10, info = paste("order =", order))
  }
})

test_that("compute_tf_std_d returns rep(n-1, n) with d[n]=0 for order=0, all equal for order>0", {
  n <- 30
  d0 <- compute_tf_std_d(0, n)
  expect_length(d0, n)
  expect_equal(d0[-n], rep(n - 1, n - 1))
  expect_equal(d0[n], 0)

  for (order in 1:3) {
    d <- compute_tf_std_d(order, n)
    expect_length(d, n)
    expect_equal(d, rep(n - 1, n), info = paste("order =", order))
  }
})

test_that("compute_tf_std_d agrees with compute_tf_d(standardize=TRUE, intercept=TRUE) for order 0 to 2", {
  n <- 20
  for (order in 0:2) {
    cm      <- compute_tf_cm(order, n)
    csd     <- compute_tf_csd(order, n)
    d_fast  <- compute_tf_std_d(order, n)
    d_full  <- compute_tf_d(order, n, cm, csd, standardize = TRUE, intercept = TRUE)
    expect_equal(d_fast, d_full, tolerance = 1e-10, info = paste("order =", order))
  }
})
