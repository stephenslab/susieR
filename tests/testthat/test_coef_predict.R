context("coef and predict S3 methods")

# =============================================================================
# COEF.SUSIE - EXTRACT COEFFICIENTS
# =============================================================================

test_that("coef.susie returns correct format", {
  set.seed(1)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  coefs <- coef(fit)

  expect_length(coefs, dat$p + 1)
  expect_type(coefs, "double")
  expect_named(coefs, NULL)
})

test_that("coef.susie includes intercept as first element", {
  set.seed(2)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  coefs <- coef(fit)

  expect_equal(coefs[1], fit$intercept)
})

test_that("coef.susie computes coefficients correctly", {
  set.seed(3)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  coefs <- coef(fit)

  expected <- c(fit$intercept,
                colSums(fit$alpha * fit$mu) / fit$X_column_scale_factors)
  expect_equal(coefs, expected)
})

test_that("coef.susie handles intercept=FALSE", {
  set.seed(4)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, intercept = FALSE, verbose = FALSE)

  coefs <- coef(fit)

  expect_equal(coefs[1], 0)
})

test_that("coef.susie handles standardize=FALSE", {
  set.seed(5)
  dat <- simulate_regression(n = 100, p = 50, k = 3, center = FALSE, scale = FALSE)
  fit <- susie(dat$X, dat$y, L = 5, standardize = FALSE, verbose = FALSE)

  coefs <- coef(fit)

  expect_length(coefs, dat$p + 1)
  expect_type(coefs, "double")
})

# =============================================================================
# PREDICT.SUSIE - MAKE PREDICTIONS
# =============================================================================

test_that("predict.susie with type='coefficients' returns coef", {
  set.seed(6)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pred_coef <- predict(fit, type = "coefficients")
  expected_coef <- coef(fit)

  expect_equal(pred_coef, expected_coef)
})

test_that("predict.susie with type='coefficients' errors if newx provided", {
  set.seed(7)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  newx <- matrix(rnorm(10 * 50), 10, 50)

  expect_error(
    predict(fit, newx = newx, type = "coefficients"),
    "Do not supply newx"
  )
})

test_that("predict.susie without newx returns fitted values", {
  set.seed(8)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pred <- predict(fit, type = "response")

  expect_equal(pred, fit$fitted)
  expect_length(pred, dat$n)
})

test_that("predict.susie with newx computes predictions correctly", {
  set.seed(9)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  newx <- matrix(rnorm(20 * 50), 20, 50)
  newx <- scale(newx, center = TRUE, scale = TRUE)

  pred <- predict(fit, newx = newx, type = "response")

  expect_length(pred, 20)
  expect_type(pred, "double")

  coefs <- coef(fit)
  expected <- drop(fit$intercept + newx %*% coefs[-1])
  expect_equal(pred, expected)
})

test_that("predict.susie handles intercept=FALSE with newx", {
  set.seed(10)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, intercept = FALSE, verbose = FALSE)

  newx <- matrix(rnorm(20 * 50), 20, 50)
  newx <- scale(newx, center = TRUE, scale = TRUE)

  pred <- predict(fit, newx = newx)

  coefs <- coef(fit)
  expected <- drop(newx %*% coefs[-1])
  expect_equal(pred, expected)
})

test_that("predict.susie with NA intercept warns and uses 0", {
  set.seed(11)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_summary_stats(dat$X, dat$y)
  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)

  expect_true(is.na(fit$intercept))

  newx <- matrix(rnorm(20 * 50), 20, 50)
  newx <- scale(newx, center = TRUE, scale = TRUE)

  expect_message(
    pred <- predict(fit, newx = newx),
    "intercept = 0"
  )

  coefs <- coef(fit)
  expected <- drop(newx %*% coefs[-1])
  expect_equal(pred, expected)
})

test_that("predict.susie default type is response", {
  set.seed(12)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pred1 <- predict(fit)
  pred2 <- predict(fit, type = "response")

  expect_equal(pred1, pred2)
})

test_that("predict.susie works with single new observation", {
  set.seed(13)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  newx <- matrix(rnorm(50), 1, 50)
  newx <- scale(newx, center = TRUE, scale = TRUE)

  pred <- predict(fit, newx = newx)

  expect_length(pred, 1)
  expect_type(pred, "double")
})

test_that("predict.susie handles matrix vs data.frame newx", {
  set.seed(14)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  newx_mat <- matrix(rnorm(20 * 50), 20, 50)
  newx_mat <- scale(newx_mat, center = TRUE, scale = TRUE)
  newx_df <- as.data.frame(newx_mat)

  pred_mat <- predict(fit, newx = newx_mat)
  pred_df <- predict(fit, newx = as.matrix(newx_df))

  expect_equal(pred_mat, pred_df)
})

# =============================================================================
# INTEGRATION - COEF & PREDICT CONSISTENCY
# =============================================================================

test_that("coef and predict work together consistently", {
  set.seed(15)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  coefs <- coef(fit)
  newx <- matrix(rnorm(10 * 50), 10, 50)
  newx <- scale(newx, center = TRUE, scale = TRUE)

  pred_via_predict <- predict(fit, newx = newx)
  pred_via_coef <- drop(coefs[1] + newx %*% coefs[-1])

  expect_equal(pred_via_predict, pred_via_coef)
})

test_that("coef.susie with all V=0 returns zero coefficients", {
  set.seed(16)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V <- rep(0, 5)
  fit$alpha <- matrix(1/dat$p, 5, dat$p)
  fit$mu <- matrix(0, 5, dat$p)

  coefs <- coef(fit)

  expect_equal(coefs[-1], rep(0, dat$p))
})

test_that("predict.susie with standardize=FALSE", {
  set.seed(17)
  dat <- simulate_regression(n = 100, p = 50, k = 3, center = FALSE, scale = FALSE)
  fit <- susie(dat$X, dat$y, L = 5, standardize = FALSE, verbose = FALSE)

  pred <- predict(fit)

  expect_equal(pred, fit$fitted)
  expect_length(pred, dat$n)
})
