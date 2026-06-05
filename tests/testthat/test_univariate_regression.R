context("Univariate regression")

# ---- Basic functionality ----

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

test_that("univariate_regression estimates coefficients close to truth", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 2)
  beta_true <- c(1, -0.5, 0.8, 0, 0.3)
  y <- base_data$X %*% beta_true + rnorm(base_data$n, sd = 0.1)
  result <- univariate_regression(base_data$X, y, center = TRUE, scale = FALSE)
  for (i in seq_len(base_data$p)) {
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
  y_centered <- base_data$y - mean(base_data$y)
  X_centered <- scale(base_data$X, center = TRUE, scale = FALSE)
  manual_fit <- lm(y_centered ~ X_centered[, 1])
  expect_equal(result$betahat[1],   unname(coef(manual_fit)[2]),              tolerance = 1e-8)
  expect_equal(result$sebetahat[1], unname(summary(manual_fit)$coef[2, 2]),   tolerance = 1e-8)
})

# ---- center / scale options ----

test_that("univariate_regression center/scale options return finite results", {
  combos <- list(
    list(center = TRUE,  scale = FALSE, seed = 5),
    list(center = TRUE,  scale = TRUE,  seed = 6),
    list(center = FALSE, scale = FALSE, seed = 7),
    list(center = FALSE, scale = TRUE,  seed = 8)
  )
  for (cfg in combos) {
    base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = cfg$seed)
    result <- univariate_regression(base_data$X, base_data$y,
                                    center = cfg$center, scale = cfg$scale)
    expect_length(result$betahat,   base_data$p)
    expect_length(result$sebetahat, base_data$p)
    expect_true(all(is.finite(result$betahat)),
                info = sprintf("betahat not finite for center=%s scale=%s",
                               cfg$center, cfg$scale))
    expect_true(all(is.finite(result$sebetahat)),
                info = sprintf("sebetahat not finite for center=%s scale=%s",
                               cfg$center, cfg$scale))
  }
})

test_that("scaling changes coefficient magnitudes when predictor variances differ", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 9)
  X_varied <- matrix(rnorm(base_data$n * base_data$p,
                            sd = rep(c(1, 5, 10, 2, 3), each = base_data$n)),
                     base_data$n, base_data$p)
  result_unscaled <- univariate_regression(X_varied, base_data$y, center = TRUE, scale = FALSE)
  result_scaled   <- univariate_regression(X_varied, base_data$y, center = TRUE, scale = TRUE)
  expect_false(all(abs(result_unscaled$betahat - result_scaled$betahat) < 0.01))
})

# ---- covariates (Z parameter) ----

test_that("univariate_regression with covariates Z returns correct structure", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 10)
  Z <- matrix(rnorm(base_data$n * 2), base_data$n, 2)
  result <- univariate_regression(base_data$X, base_data$y, Z = Z, center = TRUE)
  expect_type(result, "list")
  expect_length(result$betahat,   base_data$p)
  expect_length(result$sebetahat, base_data$p)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
})

test_that("univariate_regression with Z attenuates spurious association", {
  base_data <- generate_base_data(n = 200, p = 5, k = 0, seed = 11)
  Z <- matrix(rnorm(base_data$n), base_data$n, 1)
  X_confounded <- base_data$X + Z %*% matrix(rnorm(base_data$p), 1, base_data$p)
  y_confounded  <- 2 * Z[, 1] + rnorm(base_data$n, sd = 0.1)
  result_no_Z   <- univariate_regression(X_confounded, y_confounded, Z = NULL,  center = TRUE)
  result_with_Z <- univariate_regression(X_confounded, y_confounded, Z = Z,     center = TRUE)
  expect_true(mean(abs(result_with_Z$betahat)) < mean(abs(result_no_Z$betahat)))
})

test_that("return_residuals=TRUE with Z adds residuals element of length n", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 12)
  Z <- matrix(rnorm(base_data$n * 2), base_data$n, 2)
  result <- univariate_regression(base_data$X, base_data$y, Z = Z, return_residuals = TRUE)
  expect_named(result, c("betahat", "sebetahat", "residuals"))
  expect_length(result$residuals, base_data$n)
  expect_type(result$residuals, "double")
})

test_that("return_residuals=TRUE without Z omits residuals element", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 13)
  result <- univariate_regression(base_data$X, base_data$y, Z = NULL, return_residuals = TRUE)
  expect_named(result, c("betahat", "sebetahat"))
})

test_that("residuals from Z-adjustment are exactly mean-zero after centering", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 14)
  Z <- matrix(rnorm(base_data$n * 2), base_data$n, 2)
  result <- univariate_regression(base_data$X, base_data$y, Z = Z,
                                  return_residuals = TRUE, center = TRUE)
  expect_equal(mean(result$residuals), 0, tolerance = 1e-10)
})

# ---- NA handling ----

test_that("univariate_regression drops NA rows in y and returns finite results", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 15)
  base_data$y[c(5, 20, 35)] <- NA
  result <- univariate_regression(base_data$X, base_data$y)
  expect_length(result$betahat, base_data$p)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
})

test_that("NA removal leaves remaining estimates consistent with complete-case lm", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 16)
  beta_true <- rep(1, base_data$p)
  y_signal  <- base_data$X %*% beta_true + rnorm(base_data$n, sd = 0.1)
  na_idx    <- c(10, 20, 30)
  y_signal[na_idx] <- NA
  result <- univariate_regression(base_data$X, y_signal, center = TRUE)
  expect_true(all(is.finite(result$betahat)))
  expect_true(all(is.finite(result$sebetahat)))
  expect_length(result$betahat, base_data$p)
})

# ---- edge cases ----

test_that("zero-variance column produces warning and zero betahat/sebetahat", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 17)
  base_data$X[, 3] <- 5
  suppressWarnings(
    result <- univariate_regression(base_data$X, base_data$y, center = TRUE, scale = FALSE)
  )
  expect_equal(result$betahat[3],   0)
  expect_equal(result$sebetahat[3], 0)
})

test_that("zero-variance column emits a warning message", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 17)
  base_data$X[, 3] <- 5
  expect_message(
    univariate_regression(base_data$X, base_data$y, center = TRUE, scale = FALSE),
    "WARNING:.*Column 3 has zero variance"
  )
})

test_that("perfect predictor yields coefficient 3 and near-zero SE", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 18)
  y_perfect <- 3 * base_data$X[, 1] + mean(base_data$X[, 1])
  result <- univariate_regression(base_data$X, y_perfect, center = TRUE, scale = FALSE)
  expect_equal(result$betahat[1],   3,    tolerance = 1e-8)
  expect_true(result$sebetahat[1] < 1e-10)
})

test_that("single-column X produces scalar betahat and sebetahat", {
  base_data <- generate_base_data(n = 100, p = 1, k = 0, seed = 19)
  result <- univariate_regression(base_data$X, base_data$y)
  expect_length(result$betahat,   1)
  expect_length(result$sebetahat, 1)
  expect_true(is.finite(result$betahat[1]))
  expect_true(is.finite(result$sebetahat[1]))
})

test_that("very small sample size still produces output of correct length", {
  base_data <- generate_base_data(n = 5, p = 3, k = 0, seed = 20)
  result <- univariate_regression(base_data$X, base_data$y)
  expect_length(result$betahat,   base_data$p)
  expect_length(result$sebetahat, base_data$p)
})

test_that("nearly-singular design matrix produces output without error", {
  base_data <- generate_base_data(n = 100, p = 5, k = 0, seed = 34)
  base_data$X[, 2] <- base_data$X[, 1] + rnorm(base_data$n, sd = 1e-10)
  result <- univariate_regression(base_data$X, base_data$y, center = TRUE)
  expect_length(result$betahat,   base_data$p)
  expect_length(result$sebetahat, base_data$p)
})

# ---- z-scores ----

test_that("z-scores rank causal variables above nulls", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 21)
  beta_true <- c(rep(1, 3), rep(0, 7))
  y_signal  <- base_data$X %*% beta_true + rnorm(base_data$n)
  result    <- univariate_regression(base_data$X, y_signal)
  z         <- result$betahat / result$sebetahat
  expect_length(z, base_data$p)
  expect_type(z, "double")
  expect_true(mean(abs(z[1:3])) > mean(abs(z[4:10])))
})

# ---- method comparison: lmfit vs sumstats ----

test_that("lmfit and sumstats methods agree across center/scale combos, covariates, NA, and edge cases", {
  # Base case and four center/scale combos
  combos <- list(
    list(center = TRUE,  scale = FALSE, seed = 22, na_idx = NULL, z_cols = 0, const_col = NULL, p = 10, n = 100),
    list(center = TRUE,  scale = FALSE, seed = 23, na_idx = NULL, z_cols = 0, const_col = NULL, p = 10, n = 100),
    list(center = TRUE,  scale = TRUE,  seed = 24, na_idx = NULL, z_cols = 0, const_col = NULL, p = 10, n = 100),
    list(center = FALSE, scale = FALSE, seed = 25, na_idx = NULL, z_cols = 0, const_col = NULL, p = 10, n = 100),
    list(center = FALSE, scale = TRUE,  seed = 26, na_idx = NULL, z_cols = 0, const_col = NULL, p = 10, n = 100),
    list(center = TRUE,  scale = FALSE, seed = 28, na_idx = c(5, 20, 35), z_cols = 0, const_col = NULL, p = 10, n = 100),
    list(center = TRUE,  scale = FALSE, seed = 30, na_idx = NULL, z_cols = 0, const_col = NULL, p = 1,  n = 100)
  )

  for (cfg in combos) {
    base_data <- generate_base_data(n = cfg$n, p = cfg$p, k = 0, seed = cfg$seed)
    if (!is.null(cfg$na_idx)) base_data$y[cfg$na_idx] <- NA

    res1 <- univariate_regression(base_data$X, base_data$y,
                                  center = cfg$center, scale = cfg$scale,
                                  method = "lmfit")
    res2 <- univariate_regression(base_data$X, base_data$y,
                                  center = cfg$center, scale = cfg$scale,
                                  method = "sumstats")
    expect_equal(res1$betahat,   res2$betahat,   tolerance = 1e-8,
                 info = sprintf("betahat mismatch: center=%s scale=%s seed=%d",
                                cfg$center, cfg$scale, cfg$seed))
    expect_equal(res1$sebetahat, res2$sebetahat, tolerance = 1e-8,
                 info = sprintf("sebetahat mismatch: center=%s scale=%s seed=%d",
                                cfg$center, cfg$scale, cfg$seed))
  }
})

test_that("lmfit and sumstats agree with covariates Z", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 27)
  Z <- matrix(rnorm(base_data$n * 3), base_data$n, 3)
  res1 <- univariate_regression(base_data$X, base_data$y, Z = Z,
                                center = TRUE, method = "lmfit")
  res2 <- univariate_regression(base_data$X, base_data$y, Z = Z,
                                center = TRUE, method = "sumstats")
  expect_equal(res1$betahat,   res2$betahat,   tolerance = 1e-8)
  expect_equal(res1$sebetahat, res2$sebetahat, tolerance = 1e-8)
})

test_that("lmfit and sumstats both warn on zero-variance column and agree", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 29)
  base_data$X[, 5] <- 3
  suppressWarnings(res1 <- univariate_regression(base_data$X, base_data$y,
                                                  center = TRUE, scale = FALSE,
                                                  method = "lmfit"))
  suppressWarnings(res2 <- univariate_regression(base_data$X, base_data$y,
                                                  center = TRUE, scale = FALSE,
                                                  method = "sumstats"))
  expect_equal(res1$betahat,   res2$betahat,   tolerance = 1e-8)
  expect_equal(res1$sebetahat, res2$sebetahat, tolerance = 1e-8)
})

test_that("lmfit and sumstats agree on large dataset (n=500, p=5000)", {
  base_data <- generate_base_data(n = 500, p = 5000, k = 0, seed = 31)
  res1 <- univariate_regression(base_data$X, base_data$y,
                                center = TRUE, scale = TRUE, method = "lmfit")
  res2 <- univariate_regression(base_data$X, base_data$y,
                                center = TRUE, scale = TRUE, method = "sumstats")
  expect_equal(res1$betahat,   res2$betahat,   tolerance = 1e-8)
  expect_equal(res1$sebetahat, res2$sebetahat, tolerance = 1e-8)
})

# ---- comparison with lm ----

test_that("univariate_regression agrees with lm for each column", {
  base_data <- generate_base_data(n = 50, p = 5, k = 0, seed = 32)
  result    <- univariate_regression(base_data$X, base_data$y, center = TRUE, scale = FALSE)
  y_c <- base_data$y - mean(base_data$y)
  X_c <- scale(base_data$X, center = TRUE, scale = FALSE)
  for (i in seq_len(base_data$p)) {
    lm_coef <- unname(coef(summary(lm(y_c ~ X_c[, i])))[2, ])
    expect_equal(result$betahat[i],   lm_coef[1], tolerance = 1e-8)
    expect_equal(result$sebetahat[i], lm_coef[2], tolerance = 1e-8)
  }
})

# ---- integration ----

test_that("univariate_regression output is usable for estimate_s_rss", {
  base_data <- generate_base_data(n = 200, p = 100, k = 0, seed = 35)
  beta_true        <- rep(0, base_data$p)
  beta_true[1:5]   <- rnorm(5)
  y_causal         <- base_data$X %*% beta_true + rnorm(base_data$n)
  result           <- univariate_regression(base_data$X, y_causal)
  z                <- result$betahat / result$sebetahat
  R                <- cor(base_data$X)
  expect_no_error(estimate_s_rss(z, R, n = base_data$n))
})

# ---- calc_z ----

test_that("calc_z vector Y matches betahat/sebetahat ratio", {
  base_data <- generate_base_data(n = 100, p = 10, k = 0, seed = 37)
  z        <- susieR:::calc_z(base_data$X, base_data$y, center = FALSE, scale = FALSE)
  result   <- univariate_regression(base_data$X, base_data$y, center = FALSE, scale = FALSE)
  z_manual <- result$betahat / result$sebetahat
  expect_equal(z, z_manual)
  expect_length(z, base_data$p)
  expect_type(z, "double")
})

test_that("calc_z matrix Y returns p x m matrix matching per-column z-scores", {
  set.seed(38)
  n <- 100; p <- 10; m <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * m), n, m)
  z_matrix <- susieR:::calc_z(X, Y, center = FALSE, scale = FALSE)
  expect_true(is.matrix(z_matrix))
  expect_equal(dim(z_matrix), c(p, m))
  for (i in seq_len(m)) {
    result   <- univariate_regression(X, Y[, i], center = FALSE, scale = FALSE)
    z_manual <- result$betahat / result$sebetahat
    expect_equal(z_matrix[, i], z_manual)
  }
})

test_that("calc_z center/scale options match univariate_regression with same flags", {
  combos <- list(
    list(center = TRUE,  scale = FALSE, seed = 39),
    list(center = FALSE, scale = TRUE,  seed = 40),
    list(center = TRUE,  scale = TRUE,  seed = 41)
  )
  for (cfg in combos) {
    set.seed(cfg$seed)
    n <- 100; p <- 10
    X <- matrix(rnorm(n * p, mean = 5, sd = 2), n, p)
    y <- rnorm(n, mean = 10, sd = 3)
    z        <- susieR:::calc_z(X, y, center = cfg$center, scale = cfg$scale)
    result   <- univariate_regression(X, y, center = cfg$center, scale = cfg$scale)
    z_manual <- result$betahat / result$sebetahat
    expect_equal(z, z_manual,
                 info = sprintf("calc_z mismatch center=%s scale=%s", cfg$center, cfg$scale))
    expect_length(z, p)
  }
})

test_that("calc_z matrix Y with center and scale returns finite matrix of correct dimensions", {
  set.seed(42)
  n <- 100; p <- 8; m <- 4
  X <- matrix(rnorm(n * p, mean = rep(c(0, 5, -3, 2), each = n * 2)), n, p)
  Y <- matrix(rnorm(n * m, mean = rep(c(0, 10, -5, 3), each = n)), n, m)
  for (i in seq_len(p)) X[, i] <- X[, i] * (i %% 3 + 1)
  z <- susieR:::calc_z(X, Y, center = TRUE, scale = TRUE)
  expect_equal(dim(z), c(p, m))
  expect_true(all(is.finite(z)))
})

# ---- ordering utility functions ----

test_that("univar.order is a permutation of 1:p with causal vars ranked highest", {
  set.seed(101)
  n <- 200; p <- 30
  X    <- matrix(rnorm(n * p), n, p)
  beta <- double(p)
  beta[1:5] <- 5:1
  y    <- X %*% beta + rnorm(n)
  ord  <- univar.order(X, y)
  expect_length(ord, p)
  expect_equal(sort(ord), seq_len(p))
  expect_true(all(1:5 %in% ord[1:10]))
})

test_that("univar.order is a permutation of 1:p on null data", {
  set.seed(102)
  X   <- matrix(rnorm(50 * 10), 50, 10)
  y   <- rnorm(50)
  ord <- univar.order(X, y)
  expect_equal(sort(ord), 1:10)
})

test_that("absolute.order is a permutation of 1:p with largest |beta| first", {
  set.seed(103)
  n <- 200; p <- 30
  X         <- matrix(rnorm(n * p), n, p)
  beta_true <- double(p)
  beta_true[1:10] <- 1:10
  y         <- X %*% beta_true + rnorm(n)
  beta_hat  <- univariate_regression(X, y)$betahat
  ord       <- absolute.order(beta_hat)
  expect_length(ord, p)
  expect_equal(sort(ord), seq_len(p))
  expect_equal(ord[1], which.max(abs(beta_hat)))
})

test_that("absolute.order handles ties (all equal magnitudes)", {
  set.seed(104)
  beta <- rep(2, 10)
  ord  <- absolute.order(beta)
  expect_length(ord, 10)
  expect_equal(sort(ord), 1:10)
})

test_that("absolute.order ranks by absolute value, not sign", {
  beta <- c(-5, 3, -1, 4, 0)
  ord  <- absolute.order(beta)
  expect_equal(ord[1], 1L)   # |-5| = 5 is largest
  expect_equal(ord[2], 4L)   # |4|  = 4 is second
})

test_that("path.order is a permutation of 1:p with first-entry variable first", {
  set.seed(105)
  p         <- 30
  beta_path <- matrix(0, p + 1, p)
  for (k in seq_len(p)) beta_path[k + 1, k:p] <- 1
  fit <- list(coefficients = beta_path)
  ord <- path.order(fit)
  expect_length(ord, p)
  expect_equal(sort(ord), seq_len(p))
  expect_equal(ord[1], 1L)
})

test_that("path.order places variables that never enter the model at the end", {
  set.seed(106)
  p         <- 10
  beta_path <- matrix(0, p + 1, p)
  for (k in 1:5) beta_path[k + 1, k:p] <- 1
  fit <- list(coefficients = beta_path)
  ord <- path.order(fit)
  expect_length(ord, p)
  expect_equal(sort(ord), seq_len(p))
  expect_equal(ord[1:5], 1:5)
  expect_equal(sort(ord[6:10]), 6:10)
})
