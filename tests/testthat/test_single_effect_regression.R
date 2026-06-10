context("Single Effect Regression")

# ---- single_effect_regression structure and validity ----

test_that("single_effect_regression returns required fields with correct dimensions", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_type(result, "list")
  for (field in c("alpha", "mu", "mu2", "lbf", "lbf_variable", "V")) {
    expect_true(field %in% names(result), label = paste("field", field))
  }
  p <- setup$data$p
  expect_length(result$alpha[l, ], p)
  expect_length(result$mu[l, ], p)
  expect_length(result$mu2[l, ], p)
  expect_length(result$lbf_variable[l, ], p)
  expect_length(result$lbf[l], 1)
  expect_length(result$V[l], 1)
  expect_equal(sum(result$alpha[l, ]), 1, tolerance = 1e-10)
  expect_true(all(result$alpha[l, ] >= 0 & result$alpha[l, ] <= 1))
  expect_gte(result$V[l], 0)
  expect_true(is.finite(result$V[l]))
})

# ---- estimate_prior_method variants ----

test_that("single_effect_regression runs all valid estimate_prior_method values", {
  methods <- c("optim", "EM", "simple", "none")
  for (method in methods) {
    setup <- setup_individual_data(n = 100, p = 50, L = 5)
    setup$params$estimate_prior_method <- method
    l <- 1
    setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

    result <- single_effect_regression(setup$data, setup$params, setup$model, l)

    expect_gte(result$V[l], 0, label = paste("V >= 0 for method", method))
    expect_equal(sum(result$alpha[l, ]), 1, tolerance = 1e-10,
                 label = paste("alpha sums to 1 for method", method))
  }
})

test_that("single_effect_regression rejects invalid estimate_prior_method", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$estimate_prior_method <- "invalid_method"
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  expect_error(
    single_effect_regression(setup$data, setup$params, setup$model, l),
    "Invalid option for estimate_prior_method: invalid_method"
  )
})

test_that("single_effect_regression fixed_mixture path sets V to weighted grid mean", {
  setup <- setup_individual_data(n = 100, p = 20, L = 2)
  setup$params$estimate_prior_method <- "fixed_mixture"
  grid <- c(0.1, 0.5, 1.0)
  w    <- c(1/3, 1/3, 1/3)
  setup$params$prior_variance_grid <- grid
  setup$params$mixture_weights <- w
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_equal(sum(result$alpha[l, ]), 1, tolerance = 1e-10)
  expect_equal(result$V[l], sum(w * grid), tolerance = 1e-8)
  expect_true(!is.null(result$lbf_grid))
})

# ---- single_effect_update ----

test_that("single_effect_update produces valid posterior and Xr for all effects", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  for (l in 1:setup$params$L) {
    updated <- single_effect_update(setup$data, setup$params, setup$model, l)
    expect_equal(sum(updated$alpha[l, ]), 1, tolerance = 1e-10,
                 label = paste("alpha sums to 1 for l =", l))
    expect_gte(updated$V[l], 0, label = paste("V >= 0 for l =", l))
    expect_gte(updated$KL[l], -1e-6, label = paste("KL >= 0 for l =", l))
    expect_true("lbf" %in% names(updated))
    expect_true("lbf_variable" %in% names(updated))
    expect_true("Xr" %in% names(updated))
  }
})

# ---- posterior moment properties ----

test_that("SER posterior variance is non-negative (E[b^2] - (E[b])^2 >= 0)", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  post_second_moment <- sum(result$alpha[l, ] * result$mu2[l, ])
  post_mean_squared  <- sum(result$alpha[l, ] * result$mu[l, ])^2
  expect_gte(post_second_moment - post_mean_squared, -1e-10)
})

test_that("SER log Bayes factors and posterior moments are all finite", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_true(all(is.finite(result$lbf_variable[l, ])))
  expect_true(is.finite(result$lbf[l]))
  expect_true(all(is.finite(result$mu[l, ])))
  expect_true(all(is.finite(result$mu2[l, ])))
  expect_true(all(result$mu2[l, ] >= 0))
})

# ---- signal detection ----

test_that("SER with strong signal estimates positive V", {
  set.seed(123)
  n <- 100; p <- 50
  base_data <- generate_base_data(n, p, k = 1, signal_sd = 10, seed = NULL)
  X <- set_X_attributes(base_data$X, center = TRUE, scale = TRUE)
  y <- base_data$y - mean(base_data$y)
  data   <- structure(list(X = X, y = y, n = n, p = p, mean_y = mean(base_data$y)),
                      class = "individual")
  params <- create_base_params(L = 1, p = p, additional_params = list(
    estimate_prior_method = "optim", use_NIG = FALSE, check_null_threshold = 0.1))
  model  <- create_base_model(L = 1, p = p, n = n, X_attr = attr(X, "d"))
  model  <- compute_residuals.individual(data, params, model, 1)
  result <- single_effect_regression(data, params, model, 1)

  expect_gt(result$V, 0.1)
})

test_that("SER with no signal estimates V = 0 (null prior wins)", {
  set.seed(456)
  n <- 100; p <- 50
  base_data <- generate_base_data(n, p, k = 0, seed = NULL)
  X <- set_X_attributes(base_data$X, center = TRUE, scale = TRUE)
  y <- base_data$y - mean(base_data$y)
  data   <- structure(list(X = X, y = y, n = n, p = p, mean_y = mean(base_data$y)),
                      class = "individual")
  params <- create_base_params(L = 1, p = p, additional_params = list(
    estimate_prior_method = "optim", use_NIG = FALSE, check_null_threshold = 0.1))
  model  <- create_base_model(L = 1, p = p, n = n, X_attr = attr(X, "d"))
  model  <- compute_residuals.individual(data, params, model, 1)
  result <- single_effect_regression(data, params, model, 1)

  expect_equal(result$V, 0, tolerance = 1e-10)
})

# ---- edge case: single variable ----

test_that("SER with p=1 assigns all probability to the only variable", {
  setup <- setup_individual_data(n = 100, p = 1, L = 1)
  l <- 1
  setup$model$alpha         <- matrix(1, 1, 1)
  setup$model$mu            <- matrix(0, 1, 1)
  setup$model$mu2           <- matrix(0, 1, 1)
  setup$model$V             <- 1
  setup$model$pi            <- 1
  setup$model$predictor_weights <- attr(setup$data$X, "d")
  setup$model$lbf           <- 0
  setup$model$lbf_variable  <- matrix(0, 1, 1)
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_length(result$alpha, 1)
  expect_equal(result$alpha[1], 1, tolerance = 1e-10)
  expect_gte(result$V, 0)
})

# ---- gaussian_ser_lbf ----

test_that("gaussian_ser_lbf matches Wakefield ABF formula", {
  betahat <- 0.5; shat2 <- 0.1; V <- 1.0
  expected <- -0.5 * log(1 + V / shat2) +
    0.5 * betahat^2 * V / (shat2 * (V + shat2))
  expect_equal(gaussian_ser_lbf(betahat, shat2, V), expected, tolerance = 1e-8)
})

test_that("gaussian_ser_lbf returns 0 for non-finite betahat or shat2", {
  expect_equal(gaussian_ser_lbf(Inf,  0.1, 1.0), 0)
  expect_equal(gaussian_ser_lbf(NA,   0.1, 1.0), 0)
  expect_equal(gaussian_ser_lbf(0.5,  NA,  1.0), 0)
  expect_equal(gaussian_ser_lbf(0.5, -Inf, 1.0), 0)
})

test_that("gaussian_ser_lbf clips shat2 to machine epsilon to avoid -Inf", {
  lbf <- gaussian_ser_lbf(0.5, 0, 1.0)
  expect_true(is.finite(lbf))
  expect_gt(lbf, 0)
})

# ---- gaussian_ser_moments ----

test_that("gaussian_ser_moments matches conjugate normal posterior", {
  betahat <- 0.5; shat2 <- 0.1; V <- 1.0
  post_var  <- V * shat2 / (V + shat2)
  post_mean <- post_var / shat2 * betahat
  m <- gaussian_ser_moments(betahat, shat2, V)
  expect_equal(m$post_mean,  post_mean,              tolerance = 1e-8)
  expect_equal(m$post_mean2, post_var + post_mean^2, tolerance = 1e-8)
})

test_that("gaussian_ser_moments zeros out non-finite entries", {
  m <- gaussian_ser_moments(c(0.5, Inf, NA), c(0.1, 0.2, 0.3), V = 1.0)
  expect_equal(m$post_mean[2:3],  c(0, 0))
  expect_equal(m$post_mean2[2:3], c(0, 0))
})

# ---- gaussian_ser_posterior_e_loglik ----

test_that("gaussian_ser_posterior_e_loglik computes expected log-likelihood", {
  alpha <- c(0.3, 0.2, 0.5)
  mu    <- c(0.1, -0.2, 0.0)
  mu2   <- c(0.1^2 + 0.01, 0.2^2 + 0.02, 0.03)
  bh    <- c(0.5, -1.2, 0.0)
  sh2   <- c(0.1, 0.2, 0.3)
  e     <- gaussian_ser_posterior_e_loglik(alpha, mu, mu2, bh, sh2)
  expected <- -0.5 * sum((-2 * (alpha * mu) * bh + alpha * mu2) / sh2)
  expect_equal(e, expected, tolerance = 1e-8)
})

test_that("gaussian_ser_posterior_e_loglik ignores non-finite shat2 entries", {
  alpha <- c(0.5, 0.5)
  mu    <- c(0.1, 0.2)
  mu2   <- c(0.02, 0.05)
  bh    <- c(0.5, 1.0)
  sh2   <- c(0.1, Inf)
  e_full  <- gaussian_ser_posterior_e_loglik(alpha, mu, mu2, bh, sh2)
  e_first <- gaussian_ser_posterior_e_loglik(alpha[1], mu[1], mu2[1], bh[1], sh2[1])
  expect_equal(e_full, e_first, tolerance = 1e-8)
})

# ---- optimize_scalar_prior_variance ----

test_that("optimize_scalar_prior_variance optim path finds exact minimum", {
  # neg_loglik is minimized at V=2; optim (Brent) should recover it
  V <- optimize_scalar_prior_variance(
    V_init               = 0.5,
    estimate_prior_method = "optim",
    neg_loglik_fn        = function(Vp) (Vp - 2)^2,
    loglik_fn            = function(Vp) -(Vp - 2)^2,
    optim_init           = 0.5,
    optim_bounds         = c(0, 5),
    optim_scale          = "linear",
    check_null_threshold = 0
  )
  expect_equal(V, 2, tolerance = 1e-4)
})

test_that("optimize_scalar_prior_variance optim reverts to V_init when optimum is worse", {
  # Parabola centered at 5 but V_param_init=1: optim finds 5, which is worse
  # than V=1 on neg_loglik
  # neg_loglik(1) = 0 (minimum), neg_loglik(5) = 16
  V <- optimize_scalar_prior_variance(
    V_init               = 1.0,
    estimate_prior_method = "optim",
    neg_loglik_fn        = function(Vp) (Vp - 1)^2,
    loglik_fn            = function(Vp) -(Vp - 1)^2,
    optim_init           = 4.9,          # starts near 5 to trip the revert
    optim_bounds         = c(4, 6),      # optimizes in [4,6], best there is 4
    optim_scale          = "linear",
    check_null_threshold = 0
  )
  # best in [4,6] is 4; neg_loglik(4)=9 > neg_loglik(1)=0 -> reverts to 1
  expect_equal(V, 1.0, tolerance = 1e-6)
})

test_that("optimize_scalar_prior_variance resets V to 0 when null is better", {
  # loglik(0) + 0 >= loglik(V) => null wins
  V <- optimize_scalar_prior_variance(
    V_init               = 1.0,
    estimate_prior_method = "optim",
    neg_loglik_fn        = function(Vp) (Vp - 1)^2,
    loglik_fn            = function(Vp) -(Vp - 1)^2,  # loglik(0)=-1, loglik(1)=0
    optim_init           = 0.5,
    optim_bounds         = c(0, 5),
    optim_scale          = "linear",
    check_null_threshold = 2             # loglik(0)+2=-1+2=1 >= loglik(1)=0 -> null wins
  )
  expect_equal(V, 0)
})

test_that("optimize_scalar_prior_variance uniroot finds root on linear scale", {
  # Gradient of (Vp-2)^2 is zero at Vp=2; uniroot recovers it
  V <- optimize_scalar_prior_variance(
    V_init               = 0.5,
    estimate_prior_method = "uniroot",
    neg_loglik_fn        = function(Vp) (Vp - 2)^2,
    loglik_fn            = function(Vp) -(Vp - 2)^2,
    optim_init           = 1,
    optim_bounds         = c(0, 5),
    optim_scale          = "linear",
    check_null_threshold = 0
  )
  expect_equal(V, 2, tolerance = 1e-4)
})

test_that("optimize_scalar_prior_variance uniroot falls back to V_init on error (linear scale)", {
  # Monotone neg_loglik has same-sign gradient throughout the interval, so
  # uniroot errors; error handler returns V_init on linear scale.
  # loglik(Vp) = Vp means loglik(0)=0 < loglik(2)=2, so null-check does not fire.
  V <- optimize_scalar_prior_variance(
    V_init               = 2.0,
    estimate_prior_method = "uniroot",
    neg_loglik_fn        = function(Vp) -Vp,
    loglik_fn            = function(Vp) Vp,
    optim_init           = 2.0,
    optim_bounds         = c(1, 3),
    optim_scale          = "linear",
    check_null_threshold = 0
  )
  expect_equal(V, 2.0, tolerance = 1e-8)
})

test_that("optimize_scalar_prior_variance uniroot resets V to 0 when null model is better", {
  # Error fallback on log scale returns log(V_init); after exp(), V_new = V_init.
  # loglik_fn(0) = 0 >= loglik_fn(V_init) when loglik_fn = -Vp -> null wins.
  V <- optimize_scalar_prior_variance(
    V_init               = 3,
    estimate_prior_method = "uniroot",
    neg_loglik_fn        = function(Vp) Vp,
    loglik_fn            = function(Vp) -Vp,
    optim_init           = 1,
    optim_bounds         = c(1, 2),
    optim_scale          = "log",
    check_null_threshold = 0
  )
  expect_equal(V, 0)
})

test_that("optimize_scalar_prior_variance uniroot reverts V_new to V_init when root is worse", {
  # neg_loglik: Gaussian bump with mode at Vp=0 (log scale); neg_loglik(0) >> neg_loglik(log(3))
  # -> the found root is at ~0 but is WORSE than V_param_init=log(3), so V_new reverts to V_init=3.
  neg_ll <- function(Vp) 10 * exp(-0.5 * Vp^2) + 0.01
  loglik_fn <- function(Vp) -neg_ll(Vp)
  V <- optimize_scalar_prior_variance(
    V_init               = 3.0,
    estimate_prior_method = "uniroot",
    neg_loglik_fn        = neg_ll,
    loglik_fn            = loglik_fn,
    optim_init           = log(3.0),
    optim_bounds         = c(-5, 5),
    optim_scale          = "log",
    check_null_threshold = 0
  )
  # loglik(0)=-10.01, loglik(3)=-0.12 -> no null-reset; V_new reverted to 3
  expect_equal(V, 3.0, tolerance = 1e-6)
})

test_that("optimize_scalar_prior_variance rejects invalid estimate_prior_method", {
  expect_error(
    optimize_scalar_prior_variance(
      V_init               = 1,
      estimate_prior_method = "bogus_method",
      neg_loglik_fn        = function(x) x^2,
      loglik_fn            = function(x) -x^2,
      optim_init           = 0,
      optim_bounds         = c(-5, 5),
      optim_scale          = "linear"
    ),
    "Invalid option for estimate_prior_method: bogus_method"
  )
})

# ---- uniroot method via single_effect_regression (integration) ----

test_that("single_effect_regression with estimate_prior_method='uniroot' gives valid output", {
  set.seed(2001)
  setup <- setup_individual_data(n = 100, p = 30, L = 3)
  setup$params$estimate_prior_method <- "uniroot"
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_gte(result$V[l], 0)
  expect_true(is.finite(result$V[l]))
  expect_equal(sum(result$alpha[l, ]), 1, tolerance = 1e-10)
})

test_that("single_effect_regression uniroot estimates positive V for strong signal", {
  set.seed(2002)
  n <- 100; p <- 30
  base_data <- generate_base_data(n, p, k = 1, signal_sd = 10, seed = NULL)
  X <- set_X_attributes(base_data$X, center = TRUE, scale = TRUE)
  y <- base_data$y - mean(base_data$y)
  data   <- structure(list(X = X, y = y, n = n, p = p, mean_y = mean(base_data$y)),
                      class = "individual")
  params <- create_base_params(L = 1, p = p, additional_params = list(
    estimate_prior_method = "uniroot", use_NIG = FALSE, check_null_threshold = 0.1))
  model  <- create_base_model(L = 1, p = p, n = n, X_attr = attr(X, "d"))
  model  <- compute_residuals.individual(data, params, model, 1)

  result <- single_effect_regression(data, params, model, 1)

  expect_gt(result$V, 0)
  expect_true(is.finite(result$V))
})
