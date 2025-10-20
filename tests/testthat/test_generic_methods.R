devtools::load_all(".")

context("Generic methods infrastructure")

# =============================================================================
# GENERIC EXISTENCE
# =============================================================================

test_that("all core generics are defined", {
  # Data initialization
  expect_true(exists("configure_data", mode = "function"))
  expect_true(exists("get_var_y", mode = "function"))

  # Model initialization
  expect_true(exists("initialize_susie_model", mode = "function"))
  expect_true(exists("initialize_fitted", mode = "function"))
  expect_true(exists("validate_prior", mode = "function"))
  expect_true(exists("track_ibss_fit", mode = "function"))

  # Single effect regression
  expect_true(exists("compute_residuals", mode = "function"))
  expect_true(exists("compute_ser_statistics", mode = "function"))
  expect_true(exists("SER_posterior_e_loglik", mode = "function"))
  expect_true(exists("calculate_posterior_moments", mode = "function"))
  expect_true(exists("compute_kl", mode = "function"))
  expect_true(exists("get_ER2", mode = "function"))
  expect_true(exists("Eloglik", mode = "function"))
  expect_true(exists("loglik", mode = "function"))
  expect_true(exists("neg_loglik", mode = "function"))

  # Model updates
  expect_true(exists("update_fitted_values", mode = "function"))
  expect_true(exists("update_variance_components", mode = "function"))
  expect_true(exists("update_derived_quantities", mode = "function"))

  # Output generation
  expect_true(exists("get_scale_factors", mode = "function"))
  expect_true(exists("get_intercept", mode = "function"))
  expect_true(exists("get_fitted", mode = "function"))
  expect_true(exists("get_cs", mode = "function"))
  expect_true(exists("get_variable_names", mode = "function"))
  expect_true(exists("get_zscore", mode = "function"))
  expect_true(exists("cleanup_model", mode = "function"))
})

# =============================================================================
# METHOD DISPATCH
# =============================================================================

test_that("methods exist for all three data types", {
  classes <- c("individual", "ss", "rss_lambda")

  # Core generics that all data types must implement
  key_generics <- c(
    "configure_data",
    "get_var_y",
    "initialize_susie_model",
    "initialize_fitted",
    "compute_residuals",
    "compute_ser_statistics",
    "SER_posterior_e_loglik",
    "calculate_posterior_moments",
    "get_ER2",
    "Eloglik",
    "loglik",
    "neg_loglik",
    "update_fitted_values",
    "update_variance_components",
    "get_scale_factors",
    "get_intercept",
    "get_fitted",
    "get_cs",
    "get_variable_names",
    "cleanup_model"
  )

  for (generic in key_generics) {
    for (cls in classes) {
      method_name <- paste0(generic, ".", cls)
      expect_true(exists(method_name, mode = "function"),
                  info = paste("Missing method:", method_name))
    }
  }
})

test_that("default methods exist for optional generics", {
  default_methods <- c(
    "configure_data.default",
    "validate_prior.default",
    "track_ibss_fit.default",
    "compute_kl.default",
    "update_variance_components.default",
    "update_derived_quantities.default",
    "get_fitted.default",
    "get_zscore.default",
    "cleanup_model.default"
  )

  for (method in default_methods) {
    expect_true(exists(method, mode = "function"),
                info = paste("Missing default method:", method))
  }
})

# =============================================================================
# DEFAULT METHOD BEHAVIOR
# =============================================================================

test_that("default methods have sensible fallback behavior", {
  data <- structure(list(n = 50, p = 10), class = "test_class")
  params <- list(track_fit = FALSE)
  model <- list(alpha = matrix(1/10, 3, 10), V = c(0.1, 0.2, 0.3), sigma2 = 1)

  # configure_data.default returns data unchanged
  expect_identical(configure_data.default(data, params), data)

  # validate_prior.default returns TRUE
  expect_true(validate_prior.default(data, params, model))

  # update_derived_quantities.default returns model unchanged
  expect_identical(update_derived_quantities.default(data, params, model), model)

  # get_fitted.default and get_zscore.default return NULL
  expect_null(get_fitted.default(data, params, model))
  expect_null(get_zscore.default(data, params, model))
})

test_that("track_ibss_fit.default maintains convergence tracking", {
  data <- structure(list(), class = "test_class")
  params <- list(track_fit = FALSE)
  model <- list(alpha = matrix(1/10, 3, 10), V = c(0.1, 0.2, 0.3), sigma2 = 1)
  tracking <- list()

  # Should initialize convergence tracking
  result <- track_ibss_fit.default(data, params, model, tracking, iter = 1, elbo = c(-Inf))
  expect_true("convergence" %in% names(result))
  expect_true(is.matrix(result$convergence$prev_alpha))

  # Should update convergence tracking
  result2 <- track_ibss_fit.default(data, params, model, result, iter = 2, elbo = c(-Inf, 100))
  expect_equal(result2$convergence$prev_elbo, 100)
})

test_that("cleanup_model.default removes temporary fields", {
  data <- structure(list(), class = "test_class")
  params <- list()
  model <- list(
    alpha = matrix(1/10, 3, 10),
    mu = matrix(0, 3, 10),
    sigma2 = 1,
    V = c(0.1, 0.2, 0.3),
    # Temporary fields to remove
    null_weight = 0,
    predictor_weights = rep(1/10, 10),
    residuals = rnorm(50),
    fitted_without_l = rnorm(50)
  )

  result <- cleanup_model.default(data, params, model)

  # Keep core fields
  expect_true("alpha" %in% names(result))
  expect_true("mu" %in% names(result))
  expect_true("sigma2" %in% names(result))
  expect_true("V" %in% names(result))

  # Remove temporary fields
  expect_false("null_weight" %in% names(result))
  expect_false("predictor_weights" %in% names(result))
  expect_false("residuals" %in% names(result))
  expect_false("fitted_without_l" %in% names(result))
})
