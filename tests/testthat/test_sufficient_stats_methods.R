context("S3 methods for sufficient statistics (ss) data class")

# =============================================================================
# DATA INITIALIZATION & CONFIGURATION
# =============================================================================

test_that("configure_data.ss returns data when unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")

  result <- configure_data.ss(setup$data, setup$params)

  expect_true("ss" %in% class(result))
  expect_false("eigen_values" %in% names(result))
})

test_that("configure_data.ss adds eigen decomposition for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")

  # Remove eigen components to test they get added
  setup$data$eigen_values <- NULL
  setup$data$eigen_vectors <- NULL
  setup$data$VtXty <- NULL

  result <- configure_data.ss(setup$data, setup$params)

  expect_true("eigen_values" %in% names(result))
  expect_true("eigen_vectors" %in% names(result))
  expect_true("VtXty" %in% names(result))
})

test_that("configure_data.ss requires individual data for unmappable_effects='ash'", {
  setup <- setup_ss_data(unmappable_effects = "none")
  setup$params$unmappable_effects <- "ash"

  # Remove X and y (individual data) to test error
  setup$data$X <- NULL
  setup$data$y <- NULL

  expect_error(
    configure_data.ss(setup$data, setup$params),
    "Adaptive shrinkage \\(ash\\) requires individual-level data"
  )
})

test_that("get_var_y.ss computes variance of y", {
  setup <- setup_ss_data()

  var_y <- get_var_y.ss(setup$data)

  expect_type(var_y, "double")
  expect_length(var_y, 1)
  expect_true(var_y > 0)
  expect_equal(var_y, setup$data$yty / (setup$data$n - 1))
})

# =============================================================================
# MODEL INITIALIZATION & SETUP
# =============================================================================

test_that("initialize_susie_model.ss creates model with predictor_weights (none)", {
  setup <- setup_ss_data(unmappable_effects = "none")
  var_y <- var(setup$data$yty / setup$data$n)

  model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  expect_true("predictor_weights" %in% names(model))
  expect_length(model$predictor_weights, setup$data$p)
  expect_equal(model$predictor_weights, attr(setup$data$XtX, "d"))
})

test_that("initialize_susie_model.ss initializes omega quantities for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  var_y <- setup$data$yty / (setup$data$n - 1)

  model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  expect_true("omega_var" %in% names(model))
  expect_true("predictor_weights" %in% names(model))
  expect_true("XtOmegay" %in% names(model))
  expect_true("tau2" %in% names(model))
  expect_true("theta" %in% names(model))

  expect_equal(model$tau2, 0)
  expect_equal(model$theta, rep(0, setup$data$p))
})

test_that("initialize_fitted.ss creates XtXr", {
  setup <- setup_ss_data()

  mat_init <- list(
    alpha = setup$model$alpha,
    mu = setup$model$mu
  )

  fitted <- initialize_fitted.ss(setup$data, mat_init)

  expect_true("XtXr" %in% names(fitted))
  expect_length(fitted$XtXr, setup$data$p)
})

test_that("validate_prior.ss checks prior variance", {
  setup <- setup_ss_data()
  setup$params$check_prior <- TRUE

  # Should not error for reasonable prior variance
  expect_error(
    validate_prior.ss(setup$data, setup$params, setup$model),
    NA
  )
})

test_that("track_ibss_fit.ss delegates to default when unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")
  tracking <- list()
  iter <- 1
  elbo <- -100

  result <- track_ibss_fit.ss(setup$data, setup$params, setup$model,
                               tracking, iter, elbo)

  expect_type(result, "list")
})

test_that("track_ibss_fit.ss tracks tau2 for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  setup$params$track_fit <- TRUE
  tracking <- list()
  iter <- 1
  elbo <- -100

  result <- track_ibss_fit.ss(setup$data, setup$params, setup$model,
                               tracking, iter, elbo)

  expect_true("tau2" %in% names(result[[1]]))
  expect_equal(result[[1]]$tau2, setup$model$tau2)
})

# =============================================================================
# SINGLE EFFECT REGRESSION & ELBO
# =============================================================================

test_that("compute_residuals.ss computes residuals for unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")
  l <- 1

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)

  expect_true("residuals" %in% names(model))
  expect_true("fitted_without_l" %in% names(model))
  expect_true("residual_variance" %in% names(model))

  expect_length(model$residuals, setup$data$p)
  expect_equal(model$residual_variance, setup$model$sigma2)
})

test_that("compute_residuals.ss computes omega-weighted residuals for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  l <- 1

  # Initialize omega quantities first
  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)

  expect_true("residuals" %in% names(model))
  expect_true("predictor_weights" %in% names(model))
  expect_true("residual_variance" %in% names(model))

  expect_length(model$residuals, setup$data$p)
  expect_equal(model$residual_variance, 1)  
})

test_that("compute_ser_statistics.ss computes betahat and shat2 for unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")
  l <- 1

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  ser_stats <- compute_ser_statistics.ss(setup$data, setup$params, model, l)

  expect_true("betahat" %in% names(ser_stats))
  expect_true("shat2" %in% names(ser_stats))
  expect_true("optim_init" %in% names(ser_stats))
  expect_true("optim_bounds" %in% names(ser_stats))
  expect_true("optim_scale" %in% names(ser_stats))

  expect_length(ser_stats$betahat, setup$data$p)
  expect_length(ser_stats$shat2, setup$data$p)
  expect_equal(ser_stats$optim_scale, "log")
  expect_equal(ser_stats$optim_bounds, c(-30, 15))
})

test_that("compute_ser_statistics.ss uses linear scale for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  l <- 1

  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)
  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  ser_stats <- compute_ser_statistics.ss(setup$data, setup$params, model, l)

  expect_equal(ser_stats$optim_scale, "linear")
  expect_equal(ser_stats$optim_bounds, c(0, 1))
  expect_equal(ser_stats$optim_init, model$V[l])
})

test_that("SER_posterior_e_loglik.ss computes expected log-likelihood for unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")
  l <- 1

  setup$model$alpha[l, ] <- rep(1/setup$data$p, setup$data$p)
  setup$model$mu[l, ] <- rnorm(setup$data$p)
  setup$model$mu2[l, ] <- setup$model$mu[l, ]^2 + 0.1

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  e_loglik <- SER_posterior_e_loglik.ss(setup$data, setup$params, model, l)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
})

test_that("SER_posterior_e_loglik.ss uses omega-weighted likelihood for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  l <- 1

  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)
  setup$model$alpha[l, ] <- rep(1/setup$data$p, setup$data$p)
  setup$model$mu[l, ] <- rnorm(setup$data$p, sd = 0.01)
  setup$model$mu2[l, ] <- setup$model$mu[l, ]^2 + 0.01

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  e_loglik <- SER_posterior_e_loglik.ss(setup$data, setup$params, model, l)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
})

test_that("calculate_posterior_moments.ss computes posterior correctly", {
  setup <- setup_ss_data()
  l <- 1
  V <- 1.0

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  moments <- calculate_posterior_moments.ss(setup$data, setup$params, model, V)

  expect_true("post_mean" %in% names(moments))
  expect_true("post_mean2" %in% names(moments))
  expect_true("post_var" %in% names(moments))

  expect_length(moments$post_mean, setup$data$p)
  expect_length(moments$post_mean2, setup$data$p)
  expect_length(moments$post_var, setup$data$p)

  expect_true(all(moments$post_var >= 0))
  expect_true(all(moments$post_mean2 >= moments$post_mean^2 - 1e-10))
})

test_that("compute_kl.ss delegates to default method", {
  setup <- setup_ss_data()
  l <- 1

  setup$model$lbf <- rep(0, setup$params$L)
  setup$model$alpha[l, ] <- rep(1/setup$data$p, setup$data$p)
  setup$model$mu[l, ] <- rnorm(setup$data$p, sd = 0.1)
  setup$model$mu2[l, ] <- setup$model$mu[l, ]^2 + 0.1

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  kl <- compute_kl.ss(setup$data, setup$params, model, l)

  expect_type(kl, "double")
  expect_length(kl, 1)
})

test_that("get_ER2.ss computes expected squared residuals", {
  setup <- setup_ss_data()

  er2 <- get_ER2.ss(setup$data, setup$model)

  expect_type(er2, "double")
  expect_length(er2, 1)
  expect_true(er2 >= 0)
})

test_that("Eloglik.ss computes expected log-likelihood", {
  setup <- setup_ss_data()

  e_loglik <- Eloglik.ss(setup$data, setup$model)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
})

test_that("loglik.ss computes log Bayes factors", {
  setup <- setup_ss_data()
  l <- 1
  V <- 1.0

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  ser_stats <- compute_ser_statistics.ss(setup$data, setup$params, model, l)
  result <- loglik.ss(setup$data, setup$params, model, V, ser_stats)

  expect_true("lbf" %in% names(result))
  expect_true("lbf_model" %in% names(result))
  expect_true("alpha" %in% names(result))
  expect_true("gradient" %in% names(result))

  expect_length(result$lbf, setup$data$p)
  expect_length(result$alpha, setup$data$p)

  expect_true(all(result$alpha >= 0))
  expect_true(abs(sum(result$alpha) - 1) < 1e-10)
})

test_that("neg_loglik.ss returns negative log-likelihood for unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")
  l <- 1
  V_param <- log(1.0)

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  ser_stats <- compute_ser_statistics.ss(setup$data, setup$params, model, l)
  neg_ll <- neg_loglik.ss(setup$data, setup$params, model, V_param, ser_stats)

  expect_type(neg_ll, "double")
  expect_length(neg_ll, 1)
})

test_that("neg_loglik.ss uses unmappable objective for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  l <- 1
  V_param <- 0.5  # Linear scale

  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)
  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  ser_stats <- compute_ser_statistics.ss(setup$data, setup$params, model, l)
  neg_ll <- neg_loglik.ss(setup$data, setup$params, model, V_param, ser_stats)

  expect_type(neg_ll, "double")
  expect_length(neg_ll, 1)
})

# =============================================================================
# MODEL UPDATES & FITTING
# =============================================================================

test_that("update_fitted_values.ss updates XtXr for unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")
  l <- 1

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  setup$model$fitted_without_l <- model$fitted_without_l

  updated_model <- update_fitted_values.ss(setup$data, setup$params, setup$model, l)

  expect_true("XtXr" %in% names(updated_model))
  expect_length(updated_model$XtXr, setup$data$p)
})

test_that("update_fitted_values.ss includes theta for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  l <- 1

  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  updated_model <- update_fitted_values.ss(setup$data, setup$params, setup$model, l)

  expect_true("XtXr" %in% names(updated_model))
  expect_length(updated_model$XtXr, setup$data$p)
})

test_that("update_variance_components.ss delegates to default for unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")

  result <- update_variance_components.ss(setup$data, setup$params, setup$model)

  expect_type(result, "list")
  expect_true("sigma2" %in% names(result))
})

test_that("update_derived_quantities.ss delegates to default for unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")

  result <- update_derived_quantities.ss(setup$data, setup$params, setup$model)

  expect_type(result, "list")
})

test_that("update_derived_quantities.ss updates omega quantities for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")

  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  result <- update_derived_quantities.ss(setup$data, setup$params, setup$model)

  expect_true("omega_var" %in% names(result))
  expect_true("predictor_weights" %in% names(result))
  expect_true("XtOmegay" %in% names(result))
  expect_true("XtXr" %in% names(result))
})

# =============================================================================
# OUTPUT GENERATION & POST-PROCESSING
# =============================================================================

test_that("get_scale_factors.ss returns column scale factors", {
  setup <- setup_ss_data()

  scales <- get_scale_factors.ss(setup$data, setup$params)

  expect_length(scales, setup$data$p)
  expect_true(all(scales > 0))
  expect_equal(scales, attr(setup$data$XtX, "scaled:scale"))
})

test_that("get_intercept.ss computes intercept", {
  setup <- setup_ss_data()
  setup$params$intercept <- TRUE

  intercept <- get_intercept.ss(setup$data, setup$params, setup$model)

  expect_type(intercept, "double")
  expect_length(intercept, 1)
})

test_that("get_fitted.ss delegates to default method", {
  setup <- setup_ss_data()

  fitted <- get_fitted.ss(setup$data, setup$params, setup$model)

  # Default method returns NULL for SS data
  expect_null(fitted)
})

test_that("get_cs.ss returns NULL when coverage is NULL", {
  setup <- setup_ss_data()
  setup$params$coverage <- NULL

  cs <- get_cs.ss(setup$data, setup$params, setup$model)

  expect_null(cs)
})

test_that("get_cs.ss returns NULL when min_abs_corr is NULL", {
  setup <- setup_ss_data()
  setup$params$min_abs_corr <- NULL

  cs <- get_cs.ss(setup$data, setup$params, setup$model)

  expect_null(cs)
})

test_that("get_cs.ss computes correlation from XtX when diagonal not standardized", {
  setup <- setup_ss_data()

  # Make diagonal not 0 or 1
  diag(setup$data$XtX) <- diag(setup$data$XtX) * 1.5

  # Add strong signal to create credible set
  setup$model$alpha[1, 1] <- 0.95
  setup$model$alpha[1, -1] <- 0.05 / (setup$data$p - 1)

  cs <- get_cs.ss(setup$data, setup$params, setup$model)

  # May or may not find CS, but should not error
  expect_true(is.null(cs) || is.list(cs))
})

test_that("get_variable_names.ss assigns variable names to model", {
  setup <- setup_ss_data()
  colnames(setup$data$XtX) <- paste0("var", 1:setup$data$p)
  setup$model$pip <- rep(0.1, setup$data$p)
  setup$model$null_weight <- NULL
  setup$model$alpha <- matrix(0, 5, setup$data$p)
  setup$model$mu <- matrix(0, 5, setup$data$p)
  setup$model$mu2 <- matrix(0, 5, setup$data$p)
  setup$model$lbf_variable <- matrix(0, 5, setup$data$p)

  model_with_names <- get_variable_names.ss(setup$data, setup$model)

  expect_true(all(grepl("var", colnames(model_with_names$alpha))))
  expect_true(all(grepl("var", colnames(model_with_names$mu))))
  expect_true(all(grepl("var", colnames(model_with_names$mu2))))
  expect_true(all(grepl("var", names(model_with_names$pip))))
})

test_that("get_zscore.ss delegates to default method", {
  setup <- setup_ss_data()
  setup$params$compute_univariate_zscore <- TRUE

  z <- get_zscore.ss(setup$data, setup$params, setup$model)

  expect_true(is.null(z) || is.numeric(z))
})

test_that("cleanup_model.ss removes temporary fields for unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")

  setup$model$residuals <- rnorm(setup$data$p)

  cleaned <- cleanup_model.ss(setup$data, setup$params, setup$model)

  expect_false("residuals" %in% names(cleaned))
})

test_that("cleanup_model.ss removes omega fields for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")

  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)
  setup$model$residuals <- rnorm(setup$data$p)

  cleaned <- cleanup_model.ss(setup$data, setup$params, setup$model)

  expect_false("omega_var" %in% names(cleaned))
  expect_false("XtOmegay" %in% names(cleaned))
  expect_false("residuals" %in% names(cleaned))
})