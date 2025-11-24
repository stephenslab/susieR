context("Single Effect Regression")

# =============================================================================
# SINGLE_EFFECT_REGRESSION - Returns Correct Structure
# =============================================================================

test_that("single_effect_regression returns correct structure", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_type(result, "list")
  expect_true("alpha" %in% names(result))
  expect_true("mu" %in% names(result))
  expect_true("mu2" %in% names(result))
  expect_true("lbf" %in% names(result))
  expect_true("lbf_variable" %in% names(result))
  expect_true("V" %in% names(result))
  expect_length(result$alpha[l, ], setup$data$p)
  expect_length(result$mu[l, ], setup$data$p)
  expect_length(result$mu2[l, ], setup$data$p)
  expect_length(result$lbf_variable[l, ], setup$data$p)
  expect_length(result$lbf[l], 1)
  expect_length(result$V[l], 1)
})

# =============================================================================
# SINGLE_EFFECT_REGRESSION - Alpha Sums to 1
# =============================================================================

test_that("single_effect_regression alpha is valid probability distribution", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_equal(sum(result$alpha[l, ]), 1, tolerance = 1e-10)
  expect_true(all(result$alpha[l, ] >= 0 & result$alpha[l, ] <= 1))
})

# =============================================================================
# SINGLE_EFFECT_REGRESSION - V Non-negative
# =============================================================================

test_that("single_effect_regression V is non-negative and finite", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_true(result$V[l] >= 0)
  expect_true(is.finite(result$V[l]))
})

# =============================================================================
# SINGLE_EFFECT_REGRESSION - Different Estimation Methods
# =============================================================================

test_that("single_effect_regression works with estimate_prior_method='optim'", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$estimate_prior_method <- "optim"
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_true(result$V[l] >= 0)
  expect_equal(sum(result$alpha[l, ]), 1, tolerance = 1e-10)
})

test_that("single_effect_regression works with estimate_prior_method='EM'", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$estimate_prior_method <- "EM"
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_true(result$V[l] >= 0)
  expect_equal(sum(result$alpha[l, ]), 1, tolerance = 1e-10)
})

test_that("single_effect_regression works with estimate_prior_method='simple'", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$estimate_prior_method <- "simple"
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_true(result$V[l] >= 0)
  expect_equal(sum(result$alpha[l, ]), 1, tolerance = 1e-10)
})

test_that("single_effect_regression works with estimate_prior_method='none'", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$estimate_prior_method <- "none"
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_true(result$V[l] >= 0)
  expect_equal(sum(result$alpha[l, ]), 1, tolerance = 1e-10)
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

# =============================================================================
# SINGLE_EFFECT_UPDATE
# =============================================================================

test_that("single_effect_update updates all model components", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1

  updated_model <- single_effect_update(setup$data, setup$params, setup$model, l)

  expect_equal(sum(updated_model$alpha[l, ]), 1, tolerance = 1e-10)
  expect_true(updated_model$V[l] >= 0)
  expect_true(updated_model$KL[l] >= -1e-6)
  expect_true("lbf" %in% names(updated_model))
  expect_true("lbf_variable" %in% names(updated_model))
  expect_true("Xr" %in% names(updated_model))
})

test_that("single_effect_update maintains valid probability constraints", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1

  updated_model <- single_effect_update(setup$data, setup$params, setup$model, l)

  expect_equal(sum(updated_model$alpha[l, ]), 1, tolerance = 1e-10)
  expect_true(all(updated_model$alpha[l, ] >= 0))
  expect_true(all(updated_model$alpha[l, ] <= 1))
})

test_that("single_effect_update works for all effects l=1,...,L", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)

  for (l in 1:setup$params$L) {
    updated_model <- single_effect_update(setup$data, setup$params, setup$model, l)
    expect_equal(sum(updated_model$alpha[l, ]), 1, tolerance = 1e-10)
    expect_true(updated_model$V[l] >= 0)
  }
})

# =============================================================================
# MATHEMATICAL PROPERTIES
# =============================================================================

test_that("SER variance decomposition", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  post_second_moment <- sum(result$alpha[l, ] * result$mu2[l, ])
  post_mean_squared <- (sum(result$alpha[l, ] * result$mu[l, ]))^2
  post_var <- post_second_moment - post_mean_squared

  expect_true(post_var >= -1e-10)
})

test_that("SER log Bayes factors are finite", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_true(all(is.finite(result$lbf_variable[l, ])))
  expect_true(is.finite(result$lbf[l]))
})

test_that("SER posterior moments are finite", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  l <- 1
  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_true(all(is.finite(result$mu[l, ])))
  expect_true(all(is.finite(result$mu2[l, ])))
  expect_true(all(result$mu2[l, ] >= 0))
})

# =============================================================================
# SIGNAL DETECTION
# =============================================================================

test_that("SER with strong signal has large V", {
  set.seed(123)
  n <- 100
  p <- 50

  base_data <- generate_base_data(n, p, k = 1, signal_sd = 10, seed = NULL)
  X <- set_X_attributes(base_data$X, center = TRUE, scale = TRUE)
  y <- base_data$y - mean(base_data$y)

  data <- structure(
    list(X = X, y = y, n = n, p = p, mean_y = mean(base_data$y)),
    class = "individual"
  )

  params <- create_base_params(L = 1, p = p, additional_params = list(
    estimate_prior_method = "optim",
    use_servin_stephens = FALSE,
    check_null_threshold = 0.1
  ))

  model <- create_base_model(L = 1, p = p, n = n, X_attr = attr(X, "d"))

  model <- compute_residuals.individual(data, params, model, 1)
  result <- single_effect_regression(data, params, model, 1)

  expect_true(result$V > 0.1)
})

test_that("SER with no signal has V close to 0", {
  set.seed(456)
  n <- 100
  p <- 50

  base_data <- generate_base_data(n, p, k = 0, seed = NULL)
  X <- set_X_attributes(base_data$X, center = TRUE, scale = TRUE)
  y <- base_data$y - mean(base_data$y)

  data <- structure(
    list(X = X, y = y, n = n, p = p, mean_y = mean(base_data$y)),
    class = "individual"
  )

  params <- create_base_params(L = 1, p = p, additional_params = list(
    estimate_prior_method = "optim",
    use_servin_stephens = FALSE,
    check_null_threshold = 0.1
  ))

  model <- create_base_model(L = 1, p = p, n = n, X_attr = attr(X, "d"))

  model <- compute_residuals.individual(data, params, model, 1)
  result <- single_effect_regression(data, params, model, 1)

  expect_equal(result$V, 0, tolerance = 1e-10)
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("SER handles single variable (p=1)", {
  setup <- setup_individual_data(n = 100, p = 1, L = 1)
  l <- 1
  setup$model$alpha <- matrix(1, 1, 1)
  setup$model$mu <- matrix(0, 1, 1)
  setup$model$mu2 <- matrix(0, 1, 1)
  setup$model$V <- 1
  setup$model$pi <- 1
  setup$model$predictor_weights <- attr(setup$data$X, "d")
  setup$model$lbf <- 0
  setup$model$lbf_variable <- matrix(0, 1, 1)

  setup$model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  result <- single_effect_regression(setup$data, setup$params, setup$model, l)

  expect_length(result$alpha, 1)
  expect_equal(result$alpha[1], 1)
  expect_true(result$V >= 0)
})
