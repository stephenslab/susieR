context("Greedy-L outer loop in susie_workhorse")

# Shared fixtures: two real effects, low noise — easy signal for greedy-L tests.
.make_greedy_data <- function() {
  set.seed(42)
  N <- 200; J <- 100
  X <- matrix(rnorm(N * J), N, J)
  beta <- numeric(J); beta[c(10, 30)] <- c(2.5, -1.8)
  y <- X %*% beta + rnorm(N, sd = 0.3)
  list(X = X, y = y, N = N, J = J)
}

# ---- NULL L_greedy is bit-identical to fixed-L ----

test_that("L_greedy = NULL produces bit-identical output to fixed-L susie", {
  d <- .make_greedy_data()
  fit_fixed <- susie(d$X, d$y, L = 5)

  obj <- susie(d$X, d$y, L = 5, init_only = TRUE)
  obj$params$L_greedy <- NULL
  fit_direct <- susie_workhorse(obj$data, obj$params)

  expect_equal(fit_direct$alpha, fit_fixed$alpha, tolerance = 0)
  expect_equal(fit_direct$lbf,   fit_fixed$lbf,   tolerance = 0)
  expect_equal(fit_direct$elbo,  fit_fixed$elbo,   tolerance = 0)
})

# ---- greedy growth: L grows in multiples, capped at params$L ----

test_that("L_greedy grows in steps of L_greedy, capped at params$L", {
  d <- .make_greedy_data()
  obj <- susie(d$X, d$y, L = 10, init_only = TRUE)
  obj$params$L_greedy <- 3
  obj$params$lbf_min  <- 0.1
  fit <- susie_workhorse(obj$data, obj$params)

  final_L <- nrow(fit$alpha)
  expect_true(final_L %in% c(3, 6, 9, 10))
  expect_lte(final_L, 9)   # K=2 real effects saturates well before L=10
})

test_that("L_greedy >= K_true saturates in a single round", {
  d <- .make_greedy_data()
  obj <- susie(d$X, d$y, L = 12, init_only = TRUE)
  obj$params$L_greedy <- 6
  obj$params$lbf_min  <- 0.1
  fit <- susie_workhorse(obj$data, obj$params)

  expect_identical(nrow(fit$alpha), 6L)
})

test_that("K_true > L_greedy keeps growing past the first round", {
  set.seed(7)
  N <- 200; J <- 100
  Xh <- matrix(rnorm(N * J), N, J)
  bh <- numeric(J); bh[c(5, 20, 45, 70)] <- c(2.5, -2.0, 1.8, -1.5)
  yh <- Xh %*% bh + rnorm(N, sd = 0.3)

  obj <- susie(Xh, yh, L = 10, init_only = TRUE)
  obj$params$L_greedy <- 3
  obj$params$lbf_min  <- 0.1
  fit <- susie_workhorse(obj$data, obj$params)

  expect_gte(nrow(fit$alpha), 6)
  expect_lte(nrow(fit$alpha), 10)
})

test_that("L_greedy = L stops after one round at L", {
  d <- .make_greedy_data()
  obj <- susie(d$X, d$y, L = 3, init_only = TRUE)
  obj$params$L_greedy <- 3
  fit <- susie_workhorse(obj$data, obj$params)

  expect_identical(nrow(fit$alpha), 3L)
})

# ---- L=1 edge case ----

test_that("L_greedy = 1 produces a valid single-effect model", {
  d <- .make_greedy_data()
  obj <- susie(d$X, d$y, L = 5, init_only = TRUE)
  obj$params$L_greedy <- 1
  obj$params$lbf_min  <- 0.1
  fit <- susie_workhorse(obj$data, obj$params)

  expect_gte(nrow(fit$alpha), 1L)
  expect_equal(rowSums(fit$alpha), rep(1, nrow(fit$alpha)), tolerance = 1e-10)
})

# ---- exposed through susie / susie_ss / susie_rss interfaces ----

test_that("L_greedy is exposed through susie, susie_ss, and susie_rss interfaces", {
  d <- .make_greedy_data()

  fit <- susie(d$X, d$y, L = 8, L_greedy = 3, greedy_lbf_cutoff = 0.1,
               verbose = FALSE)
  expect_lte(nrow(fit$alpha), 8)
  expect_gte(nrow(fit$alpha), 3)

  y_vec <- drop(d$y)
  ss <- compute_suff_stat(d$X, y_vec, standardize = TRUE)
  fit_ss <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 8,
                     L_greedy = 3, greedy_lbf_cutoff = 0.1,
                     verbose = FALSE)
  expect_lte(nrow(fit_ss$alpha), 8)
  expect_gte(nrow(fit_ss$alpha), 3)

  z <- as.vector(crossprod(scale(d$X), drop(scale(y_vec))) / sqrt(nrow(d$X) - 1))
  R <- cor(d$X)
  fit_rss <- susie_rss(z = z, R = R, n = nrow(d$X), L = 8,
                       L_greedy = 3, greedy_lbf_cutoff = 0.1,
                       verbose = FALSE)
  expect_lte(nrow(fit_rss$alpha), 8)
  expect_gte(nrow(fit_rss$alpha), 3)
})
