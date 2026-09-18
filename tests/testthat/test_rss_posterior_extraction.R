# Posterior extraction on RSS / SER-fallback fits, where pecotmr saves fits with
# X_column_scale_factors and intercept absent.

test_that("assume_unit_scale returns 1 (with a hint) on NULL, else unchanged", {
  expect_equal(susieSlide:::assume_unit_scale(c(1, 2, 3)), c(1, 2, 3))
  expect_message(v <- susieSlide:::assume_unit_scale(NULL), "assuming 1")
  expect_equal(v, 1)
})

test_that("scale stripped: getters are full-length and equal the untrimmed fit, hint fires", {
  set.seed(1); p <- 30; z <- rnorm(p); z[5] <- 6; z[20] <- -5; R <- diag(p)
  fit <- suppressMessages(susie_rss(z = z, R = R, n = 2000))
  f0 <- fit; f0$X_column_scale_factors <- NULL
  expect_message(pm0 <- susie_get_posterior_mean(f0), "assuming 1")
  psd0 <- suppressMessages(susie_get_posterior_sd(f0))
  expect_length(pm0, p)
  expect_length(psd0, p)
  expect_equal(unname(pm0), unname(susie_get_posterior_mean(fit)))
  expect_equal(unname(psd0), unname(susie_get_posterior_sd(fit)))
  # and it matches colSums(alpha*mu) (scale-1 ground truth) to zero difference
  expect_equal(unname(pm0), unname(colSums(fit$alpha * fit$mu)))
})

test_that("intercept stripped: coef length p+1 and coef[1] == 0 (fixes the off-by-one)", {
  set.seed(3); p <- 20; z <- rnorm(p); z[2] <- 5; R <- diag(p)
  fit <- suppressMessages(susie_ser(z = z, n = 1000))
  f0 <- fit; f0$X_column_scale_factors <- NULL; f0$intercept <- NULL
  cf <- suppressMessages(coef(f0))
  expect_length(cf, ncol(f0$mu) + 1L)
  expect_equal(cf[[1]], 0)
  expect_equal(unname(cf[-1]), unname(colSums(f0$alpha * f0$mu)))
})

test_that("individual-level fit: coef unchanged (regression guard, no hint)", {
  set.seed(4); n <- 200; p <- 20
  X <- matrix(rnorm(n * p), n, p); X[, 1] <- X[, 1] * 5  # unequal column SDs
  y <- X[, 3] * 0.5 + rnorm(n)
  fit <- suppressMessages(susie_additive(X, y, L = 5))
  expect_false(is.null(fit$X_column_scale_factors))  # present -> guard is a no-op
  cf <- coef(fit)                                     # no "assuming 1" hint
  expect_length(cf, p + 1L)
})

test_that("NULL intercept / NULL scale guards do not error", {
  set.seed(5); p <- 12; z <- rnorm(p); z[4] <- 5; R <- diag(p)
  fit <- suppressMessages(susie_rss(z = z, R = R, n = 800))
  fit$intercept <- NULL; fit$X_column_scale_factors <- NULL
  expect_error(suppressMessages(coef(fit)), NA)                       # isTRUE(is.na(NULL)) ok
  expect_error(suppressMessages(susie_get_posterior_mean(fit)), NA)
  expect_error(suppressMessages(susie_get_posterior_samples(fit, 3)), NA)
})
