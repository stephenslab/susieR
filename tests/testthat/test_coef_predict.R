context("coef and predict S3 methods")

# ---- coef.susie ----

test_that("coef.susie returns a double vector of length p+1 with no names", {
  set.seed(1)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  coefs <- coef(fit)
  expect_length(coefs, dat$p + 1)
  expect_type(coefs, "double")
  expect_named(coefs, NULL)
})

test_that("coef.susie first element equals fit$intercept", {
  set.seed(2)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  expect_equal(coef(fit)[1], fit$intercept, tolerance = 1e-8)
})

test_that("coef.susie computes (alpha*mu)/X_column_scale_factors correctly", {
  set.seed(3)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  expected <- c(fit$intercept,
                colSums(fit$alpha * fit$mu) / fit$X_column_scale_factors)
  expect_equal(coef(fit), expected, tolerance = 1e-8)
})

test_that("coef.susie returns intercept=0 when intercept=FALSE", {
  set.seed(4)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, intercept = FALSE, verbose = FALSE)
  expect_equal(coef(fit)[1], 0, tolerance = 1e-8)
})

test_that("coef.susie works with standardize=FALSE", {
  set.seed(5)
  dat <- simulate_regression(n = 100, p = 50, k = 3, center = FALSE, scale = FALSE)
  fit <- susie(dat$X, dat$y, L = 5, standardize = FALSE, verbose = FALSE)
  coefs <- coef(fit)
  expect_length(coefs, dat$p + 1)
  expect_type(coefs, "double")
})

test_that("coef.susie returns zero non-intercept coefficients when all mu=0", {
  set.seed(16)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$V   <- rep(0, 5)
  fit$alpha <- matrix(1 / dat$p, 5, dat$p)
  fit$mu  <- matrix(0, 5, dat$p)
  expect_equal(coef(fit)[-1], rep(0, dat$p), tolerance = 1e-8)
})

test_that("coef.susie includes theta contribution when unmappable_effects='inf'", {
  set.seed(401)
  n <- 60; p <- 25
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:3] <- 1
  y <- as.vector(X %*% beta + rnorm(n))
  fit <- suppressWarnings(suppressMessages(
    susie(X, y, L = 5, unmappable_effects = "inf", max_iter = 5, verbose = FALSE)
  ))
  expect_false(is.null(fit$theta))
  expect_length(fit$theta, p)
  coefs <- coef(fit)
  expect_length(coefs, p + 1)
  expected <- c(fit$intercept,
                (colSums(fit$alpha * fit$mu) + fit$theta) / fit$X_column_scale_factors)
  expect_equal(coefs, expected, tolerance = 1e-8)
})

test_that("coef.susie without theta differs from coef with theta (unmappable_effects='none' vs 'inf')", {
  set.seed(402)
  n <- 60; p <- 25
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[1:3] <- 1
  y <- as.vector(X %*% beta + rnorm(n))
  fit_inf <- suppressWarnings(suppressMessages(
    susie(X, y, L = 5, unmappable_effects = "inf", max_iter = 5, verbose = FALSE)
  ))
  fit_none <- suppressWarnings(
    susie(X, y, L = 5, unmappable_effects = "none", max_iter = 5, verbose = FALSE)
  )
  expect_length(coef(fit_inf),  p + 1)
  expect_length(coef(fit_none), p + 1)
  expect_false(is.null(fit_inf$theta))
  expect_null(fit_none$theta)
})

# ---- predict.susie ----

test_that("predict.susie with type='coefficients' returns coef(fit)", {
  set.seed(6)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  expect_equal(predict(fit, type = "coefficients"), coef(fit), tolerance = 1e-8)
})

test_that("predict.susie type='coefficients' errors when newx is supplied", {
  set.seed(7)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  newx <- matrix(rnorm(10 * 50), 10, 50)
  expect_error(predict(fit, newx = newx, type = "coefficients"), "Do not supply newx")
})

test_that("predict.susie without newx returns fit$fitted", {
  set.seed(8)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  pred <- predict(fit, type = "response")
  expect_equal(pred, fit$fitted, tolerance = 1e-8)
  expect_length(pred, dat$n)
})

test_that("predict.susie with newx computes intercept + newx %*% coefs[-1]", {
  set.seed(9)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  newx <- scale(matrix(rnorm(20 * 50), 20, 50))
  pred <- predict(fit, newx = newx, type = "response")
  expected <- drop(fit$intercept + newx %*% coef(fit)[-1])
  expect_equal(pred, expected, tolerance = 1e-8)
  expect_length(pred, 20)
})

test_that("predict.susie with intercept=FALSE omits intercept term", {
  set.seed(10)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, intercept = FALSE, verbose = FALSE)
  newx <- scale(matrix(rnorm(20 * 50), 20, 50))
  pred <- predict(fit, newx = newx)
  expected <- drop(newx %*% coef(fit)[-1])
  expect_equal(pred, expected, tolerance = 1e-8)
})

test_that("predict.susie with NA intercept warns and substitutes 0", {
  set.seed(11)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss  <- compute_summary_stats(dat$X, dat$y)
  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)
  expect_true(is.na(fit$intercept))
  newx <- scale(matrix(rnorm(20 * 50), 20, 50))
  expect_message(
    pred <- predict(fit, newx = newx),
    "intercept = 0"
  )
  expected <- drop(newx %*% coef(fit)[-1])
  expect_equal(pred, expected, tolerance = 1e-8)
})

test_that("predict.susie default type is 'response'", {
  set.seed(12)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  expect_equal(predict(fit), predict(fit, type = "response"), tolerance = 1e-8)
})

test_that("predict.susie handles a single new observation (1 x p newx)", {
  set.seed(13)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  newx <- scale(matrix(rnorm(50), 1, 50))
  pred <- predict(fit, newx = newx)
  expect_length(pred, 1)
  expect_type(pred, "double")
})

test_that("predict.susie gives identical results for matrix and data.frame newx", {
  set.seed(14)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  newx_mat <- scale(matrix(rnorm(20 * 50), 20, 50))
  pred_mat <- predict(fit, newx = newx_mat)
  pred_df  <- predict(fit, newx = as.matrix(as.data.frame(newx_mat)))
  expect_equal(pred_mat, pred_df, tolerance = 1e-8)
})

test_that("predict.susie works with standardize=FALSE fit", {
  set.seed(17)
  dat <- simulate_regression(n = 100, p = 50, k = 3, center = FALSE, scale = FALSE)
  fit <- susie(dat$X, dat$y, L = 5, standardize = FALSE, verbose = FALSE)
  pred <- predict(fit)
  expect_equal(pred, fit$fitted, tolerance = 1e-8)
  expect_length(pred, dat$n)
})

# ---- coef / predict consistency ----

test_that("predict(newx) equals intercept + newx %*% coef[-1]", {
  set.seed(15)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  coefs <- coef(fit)
  newx  <- scale(matrix(rnorm(10 * 50), 10, 50))
  expect_equal(predict(fit, newx = newx),
               drop(coefs[1] + newx %*% coefs[-1]),
               tolerance = 1e-8)
})
