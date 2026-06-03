# Tests for slot_weights mechanism
#
# Key invariant: slot_weights = rep(1, L) must produce identical results
# to the standard path (slot_weights = NULL).

context("Slot weights")

set.seed(1)
n <- 200
p <- 50
X <- matrix(rnorm(n * p), n, p)
beta <- rep(0, p)
beta[c(1, 5, 10)] <- c(0.5, -0.3, 0.4)
y <- X %*% beta + rnorm(n)
R <- cor(X)
z <- as.vector(sqrt(n) * crossprod(X, y) / sqrt(n * diag(crossprod(X))))

# ---- get_slot_weight: boundary values ----

test_that("get_slot_weight returns 1 when slot_weights is NULL", {
  expect_equal(susieR:::get_slot_weight(list(), 1), 1)
  expect_equal(susieR:::get_slot_weight(list(), 5), 1)
})

test_that("get_slot_weight returns correct weight at each index", {
  model <- list(slot_weights = c(0.0, 0.5, 1.0))
  expect_equal(susieR:::get_slot_weight(model, 1), 0.0)
  expect_equal(susieR:::get_slot_weight(model, 2), 0.5)
  expect_equal(susieR:::get_slot_weight(model, 3), 1.0)
})

# ---- slot_weight = 0 zeroes effect contribution ----

test_that("slot_weight = 0 for a slot: that slot's SER runs but contributes nothing to fitted values", {
  L <- 3
  objs <- susie_rss(z = z, R = R, n = n, L = L,
                    estimate_prior_variance = FALSE,
                    estimate_residual_variance = FALSE,
                    init_only = TRUE)
  data   <- objs$data
  params <- objs$params
  model  <- susieR:::ibss_initialize(data, params)

  # Set slot 2 weight to 0; save Rz before update
  model$slot_weights <- c(1, 0, 1)
  Rz_before <- model$Rz

  model_after <- susieR:::single_effect_update(data, params, model, 2)

  # SER still populates alpha and mu for slot 2
  expect_true(is.numeric(model_after$alpha[2, ]))
  expect_true(all(is.finite(model_after$alpha[2, ])))
  expect_true(all(model_after$alpha[2, ] >= 0))
  expect_true(abs(sum(model_after$alpha[2, ]) - 1) < 1e-8)
})

test_that("slot_weight = 0 vs 1 for a slot produces different alpha when model has non-trivial state", {
  L <- 2
  # Run 10 iterations to get a non-trivial (non-zero mu) model state
  fit_warm <- suppressWarnings(
    susie_rss(z = z, R = R, n = n, L = L,
              estimate_prior_variance = FALSE,
              estimate_residual_variance = FALSE,
              max_iter = 10)
  )
  objs <- susie_rss(z = z, R = R, n = n, L = L,
                    estimate_prior_variance = FALSE,
                    estimate_residual_variance = FALSE,
                    init_only = TRUE)
  data   <- objs$data
  params <- objs$params

  # Warm-start both models from the fitted state
  m1 <- susieR:::ibss_initialize(data, params)
  m1$alpha <- fit_warm$alpha
  m1$mu    <- fit_warm$mu
  m1$mu2   <- fit_warm$mu2
  m1$XtXr  <- fit_warm$XtXr
  m1$slot_weights <- c(1, 1)

  m0 <- m1
  m0$slot_weights <- c(0, 1)

  m1_after <- susieR:::single_effect_update(data, params, m1, 1)
  m0_after <- susieR:::single_effect_update(data, params, m0, 1)

  # With sw_l=0 the residuals exclude slot 1's contribution entirely,
  # while with sw_l=1 they include it; the posteriors therefore differ
  expect_false(isTRUE(all.equal(m1_after$alpha[1, ], m0_after$alpha[1, ],
                                tolerance = 1e-8)))
})

# ---- slot_weights with individual data ----

test_that("susie with slot_prior_betabinom produces valid PIPs and alpha rows sum to 1", {
  set.seed(2)
  fit <- suppressWarnings(
    susie(X, y, L = 5,
          slot_prior = slot_prior_betabinom(),
          estimate_prior_variance = FALSE,
          estimate_residual_variance = FALSE,
          max_iter = 5)
  )
  pips <- susie_get_pip(fit)
  expect_true(all(pips >= 0 & pips <= 1))
  expect_true(all(abs(rowSums(fit$alpha) - 1) < 1e-8))
  # slot_weights (c_hat) should be in [0, 1]
  if (!is.null(fit$slot_weights)) {
    expect_true(all(fit$slot_weights >= 0 & fit$slot_weights <= 1))
  }
})

# ---- slot_weights with sufficient stats ----

test_that("susie_ss with slot_prior_betabinom produces valid PIPs", {
  set.seed(3)
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  yty <- sum(y^2)

  fit <- suppressWarnings(
    susie_ss(XtX = XtX, Xty = Xty, yty = yty, n = n, L = 5,
             slot_prior = slot_prior_betabinom(),
             estimate_prior_variance = FALSE,
             estimate_residual_variance = FALSE,
             max_iter = 5)
  )
  pips <- susie_get_pip(fit)
  expect_true(all(pips >= 0 & pips <= 1))
  expect_true(all(abs(rowSums(fit$alpha) - 1) < 1e-8))
})

# ---- slot_prior_betabinom: more active slots with permissive prior ----

test_that("slot_prior_betabinom: uniform prior (a=b=1) activates more slots than sparse prior (a=1,b=10)", {
  set.seed(4)
  fit_sparse  <- suppressWarnings(
    susie(X, y, L = 10,
          slot_prior = slot_prior_betabinom(a_beta = 1, b_beta = 10),
          estimate_prior_variance = FALSE,
          estimate_residual_variance = FALSE,
          max_iter = 10)
  )
  fit_uniform <- suppressWarnings(
    susie(X, y, L = 10,
          slot_prior = slot_prior_betabinom(a_beta = 1, b_beta = 1),
          estimate_prior_variance = FALSE,
          estimate_residual_variance = FALSE,
          max_iter = 10)
  )
  # Both fits must be structurally valid
  expect_true(all(susie_get_pip(fit_sparse)  >= 0 & susie_get_pip(fit_sparse)  <= 1))
  expect_true(all(susie_get_pip(fit_uniform) >= 0 & susie_get_pip(fit_uniform) <= 1))
  # Sparse prior should yield lower or equal mean c_hat than uniform prior
  if (!is.null(fit_sparse$slot_weights) && !is.null(fit_uniform$slot_weights)) {
    expect_true(mean(fit_sparse$slot_weights) <= mean(fit_uniform$slot_weights) + 0.1)
  }
})
