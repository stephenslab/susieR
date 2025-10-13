devtools::load_all(".")
source(file.path("..", "helpers", "helper_testthat.R"), local = TRUE)

context("S3 methods for individual data class")

# =============================================================================
# DATA INITIALIZATION & CONFIGURATION
# =============================================================================

test_that("configure_data.individual returns data when unmappable_effects='none'", {
  setup <- setup_individual_data()
  setup$params$unmappable_effects <- "none"

  result <- configure_data.individual(setup$data, setup$params)

  expect_true("individual" %in% class(result))
})

test_that("get_var_y.individual computes variance of y", {
  setup <- setup_individual_data()

  var_y <- get_var_y.individual(setup$data)

  expect_type(var_y, "double")
  expect_length(var_y, 1)
  expect_true(var_y > 0)
  expect_equal(var_y, var(setup$data$y))
})

# =============================================================================
# MODEL INITIALIZATION & SETUP
# =============================================================================

test_that("initialize_susie_model.individual creates model with predictor_weights", {
  setup <- setup_individual_data()
  var_y <- var(setup$data$y)

  model <- initialize_susie_model.individual(setup$data, setup$params, var_y)

  expect_true("predictor_weights" %in% names(model))
  expect_length(model$predictor_weights, setup$data$p)
  expect_equal(model$predictor_weights, attr(setup$data$X, "d"))
})

test_that("initialize_fitted.individual creates Xr", {
  setup <- setup_individual_data()

  mat_init <- list(
    alpha = setup$model$alpha,
    mu = setup$model$mu
  )

  fitted <- initialize_fitted.individual(setup$data, mat_init)

  expect_true("Xr" %in% names(fitted))
  expect_length(fitted$Xr, setup$data$n)
})

test_that("validate_prior.individual delegates to default method", {
  setup <- setup_individual_data()

  result <- validate_prior.individual(setup$data, setup$params, setup$model)

  expect_type(result, "logical")
})

test_that("track_ibss_fit.individual delegates to default method", {
  setup <- setup_individual_data()
  tracking <- list()
  iter <- 1
  elbo <- -100

  result <- track_ibss_fit.individual(setup$data, setup$params, setup$model,
                                      tracking, iter, elbo)

  expect_type(result, "list")
})

# =============================================================================
# SINGLE EFFECT REGRESSION & ELBO
# =============================================================================

test_that("compute_residuals.individual computes residuals correctly", {
  setup <- setup_individual_data()
  l <- 1

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)

  expect_true("residuals" %in% names(model))
  expect_true("fitted_without_l" %in% names(model))
  expect_true("raw_residuals" %in% names(model))
  expect_true("residual_variance" %in% names(model))

  expect_length(model$residuals, setup$data$p)
  expect_length(model$raw_residuals, setup$data$n)
})

test_that("compute_ser_statistics.individual computes betahat and shat2", {
  setup <- setup_individual_data()
  l <- 1

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  ser_stats <- compute_ser_statistics.individual(setup$data, setup$params, model, l)

  expect_true("betahat" %in% names(ser_stats))
  expect_true("shat2" %in% names(ser_stats))
  expect_true("optim_init" %in% names(ser_stats))
  expect_true("optim_bounds" %in% names(ser_stats))
  expect_true("optim_scale" %in% names(ser_stats))

  expect_length(ser_stats$betahat, setup$data$p)
  expect_length(ser_stats$shat2, setup$data$p)
  expect_true(all(ser_stats$shat2 > 0))
})

test_that("calculate_posterior_moments.individual computes posterior correctly", {
  setup <- setup_individual_data()
  l <- 1
  V <- 1.0

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  moments <- calculate_posterior_moments.individual(setup$data, setup$params, model, V)

  expect_true("post_mean" %in% names(moments))
  expect_true("post_mean2" %in% names(moments))
  expect_true("post_var" %in% names(moments))

  expect_length(moments$post_mean, setup$data$p)
  expect_length(moments$post_mean2, setup$data$p)
  expect_length(moments$post_var, setup$data$p)

  expect_true(all(moments$post_var >= 0))
  expect_true(all(moments$post_mean2 >= moments$post_mean^2 - 1e-10))
})

test_that("calculate_posterior_moments.individual handles V=0", {
  setup <- setup_individual_data()
  l <- 1
  V <- 0

  setup$params$use_servin_stephens <- TRUE
  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  moments <- calculate_posterior_moments.individual(setup$data, setup$params, model, V)

  expect_equal(moments$post_mean, rep(0, setup$data$p))
  expect_equal(moments$post_mean2, rep(0, setup$data$p))
  expect_equal(moments$post_var, rep(0, setup$data$p))
})

test_that("compute_kl.individual computes KL divergence", {
  setup <- setup_individual_data()
  l <- 1

  setup$model$lbf <- rep(0, setup$params$L)
  setup$model$alpha[l, ] <- rep(1/setup$data$p, setup$data$p)
  setup$model$mu[l, ] <- rnorm(setup$data$p, sd = 0.1)
  setup$model$mu2[l, ] <- setup$model$mu[l, ]^2 + 0.1

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  kl <- compute_kl.individual(setup$data, setup$params, model, l)

  expect_type(kl, "double")
  expect_length(kl, 1)
})

test_that("get_ER2.individual computes expected squared residuals", {
  setup <- setup_individual_data()

  er2 <- get_ER2.individual(setup$data, setup$model)

  expect_type(er2, "double")
  expect_length(er2, 1)
  expect_true(er2 >= 0)
})

test_that("Eloglik.individual computes expected log-likelihood", {
  setup <- setup_individual_data()

  e_loglik <- Eloglik.individual(setup$data, setup$model)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
})

test_that("loglik.individual computes log Bayes factors", {
  setup <- setup_individual_data()
  l <- 1
  V <- 1.0

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  ser_stats <- compute_ser_statistics.individual(setup$data, setup$params, model, l)
  result <- loglik.individual(setup$data, setup$params, model, V, ser_stats)

  expect_true("lbf" %in% names(result))
  expect_true("lbf_model" %in% names(result))
  expect_true("alpha" %in% names(result))
  expect_true("gradient" %in% names(result))

  expect_length(result$lbf, setup$data$p)
  expect_length(result$alpha, setup$data$p)

  expect_true(all(result$alpha >= 0))
  expect_true(abs(sum(result$alpha) - 1) < 1e-10)
})

test_that("neg_loglik.individual returns negative log-likelihood", {
  setup <- setup_individual_data()
  l <- 1
  V_param <- log(1.0)

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  ser_stats <- compute_ser_statistics.individual(setup$data, setup$params, model, l)
  neg_ll <- neg_loglik.individual(setup$data, setup$params, model, V_param, ser_stats)

  expect_type(neg_ll, "double")
  expect_length(neg_ll, 1)
})

test_that("SER_posterior_e_loglik.individual computes expected log-likelihood", {
  setup <- setup_individual_data()
  l <- 1

  setup$model$alpha[l, ] <- rep(1/setup$data$p, setup$data$p)
  setup$model$mu[l, ] <- rnorm(setup$data$p)
  setup$model$mu2[l, ] <- setup$model$mu[l, ]^2 + 0.1

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  e_loglik <- SER_posterior_e_loglik.individual(setup$data, setup$params, model, l)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
})

# =============================================================================
# MODEL UPDATES & FITTING
# =============================================================================

test_that("update_fitted_values.individual updates Xr", {
  setup <- setup_individual_data()
  l <- 1

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  setup$model$fitted_without_l <- model$fitted_without_l

  updated_model <- update_fitted_values.individual(setup$data, setup$params, setup$model, l)

  expect_true("Xr" %in% names(updated_model))
  expect_length(updated_model$Xr, setup$data$n)
})

test_that("update_variance_components.individual delegates to default method", {
  setup <- setup_individual_data()

  result <- update_variance_components.individual(setup$data, setup$params, setup$model)

  expect_type(result, "list")
})

test_that("update_derived_quantities.individual delegates to default method", {
  setup <- setup_individual_data()

  result <- update_derived_quantities.individual(setup$data, setup$params, setup$model)

  expect_type(result, "list")
})

# =============================================================================
# OUTPUT GENERATION & POST-PROCESSING
# =============================================================================

test_that("get_scale_factors.individual returns column scale factors", {
  setup <- setup_individual_data()

  scales <- get_scale_factors.individual(setup$data, setup$params)

  expect_length(scales, setup$data$p)
  expect_true(all(scales > 0))
  expect_equal(scales, attr(setup$data$X, "scaled:scale"))
})

test_that("get_intercept.individual computes intercept when intercept=TRUE", {
  setup <- setup_individual_data()
  setup$params$intercept <- TRUE

  intercept <- get_intercept.individual(setup$data, setup$params, setup$model)

  expect_type(intercept, "double")
  expect_length(intercept, 1)
})

test_that("get_intercept.individual returns 0 when intercept=FALSE", {
  setup <- setup_individual_data()
  setup$params$intercept <- FALSE

  intercept <- get_intercept.individual(setup$data, setup$params, setup$model)

  expect_equal(intercept, 0)
})

test_that("get_fitted.individual returns fitted values with correct length", {
  setup <- setup_individual_data()

  fitted <- get_fitted.individual(setup$data, setup$params, setup$model)

  expect_length(fitted, setup$data$n)
  expect_type(fitted, "double")
})

test_that("get_fitted.individual adds intercept when intercept=TRUE", {
  setup <- setup_individual_data()
  setup$params$intercept <- TRUE
  setup$data$mean_y <- 5.0

  fitted <- get_fitted.individual(setup$data, setup$params, setup$model)

  expect_true(any(fitted != setup$model$Xr))
})

test_that("get_fitted.individual does not add intercept when intercept=FALSE", {
  setup <- setup_individual_data()
  setup$params$intercept <- FALSE
  setup$data$mean_y <- 0

  fitted <- get_fitted.individual(setup$data, setup$params, setup$model)

  expect_equal(fitted, drop(setup$model$Xr))
})

test_that("get_cs.individual returns NULL when coverage is NULL", {
  setup <- setup_individual_data()
  setup$params$coverage <- NULL

  cs <- get_cs.individual(setup$data, setup$params, setup$model)

  expect_null(cs)
})

test_that("get_cs.individual returns NULL when min_abs_corr is NULL", {
  setup <- setup_individual_data()
  setup$params$min_abs_corr <- NULL

  cs <- get_cs.individual(setup$data, setup$params, setup$model)

  expect_null(cs)
})

test_that("get_variable_names.individual assigns variable names to model", {
  setup <- setup_individual_data()
  colnames(setup$data$X) <- paste0("var", 1:setup$data$p)
  setup$model$pip <- rep(0.1, setup$data$p)
  setup$model$null_weight <- NULL
  setup$model$alpha <- matrix(0, 5, setup$data$p)
  setup$model$mu <- matrix(0, 5, setup$data$p)
  setup$model$mu2 <- matrix(0, 5, setup$data$p)
  setup$model$lbf_variable <- matrix(0, 5, setup$data$p)

  model_with_names <- get_variable_names.individual(setup$data, setup$model)

  expect_true(all(grepl("var", colnames(model_with_names$alpha))))
  expect_true(all(grepl("var", colnames(model_with_names$mu))))
  expect_true(all(grepl("var", colnames(model_with_names$mu2))))
  expect_true(all(grepl("var", names(model_with_names$pip))))
})

test_that("get_zscore.individual computes z-scores", {
  setup <- setup_individual_data()
  setup$params$compute_univariate_zscore <- TRUE

  z <- get_zscore.individual(setup$data, setup$params, setup$model)

  expect_length(z, setup$data$p)
  expect_type(z, "double")
})

test_that("get_zscore.individual handles null_weight", {
  setup <- setup_individual_data()
  setup$params$compute_univariate_zscore <- TRUE
  setup$model$null_weight <- 0.1

  setup$data$X <- cbind(setup$data$X, 0)

  z <- get_zscore.individual(setup$data, setup$params, setup$model)

  expect_length(z, setup$data$p)
})

test_that("get_zscore.individual returns default when compute_univariate_zscore=FALSE", {
  setup <- setup_individual_data()
  setup$params$compute_univariate_zscore <- FALSE

  z <- get_zscore.individual(setup$data, setup$params, setup$model)

  expect_null(z)
})

test_that("cleanup_model.individual removes temporary fields", {
  setup <- setup_individual_data()

  setup$model$raw_residuals <- rnorm(setup$data$n)
  setup$model$residuals <- rnorm(setup$data$p)

  cleaned <- cleanup_model.individual(setup$data, setup$params, setup$model)

  expect_false("raw_residuals" %in% names(cleaned))
})