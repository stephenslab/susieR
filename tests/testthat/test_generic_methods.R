context("Generic methods infrastructure")

# ---- Generic existence ----

test_that("all core generics are defined as functions", {
  generics <- c(
    # Data initialization
    "configure_data", "get_var_y",
    # Model initialization
    "initialize_susie_model", "initialize_fitted", "validate_prior",
    "track_ibss_fit",
    # Single effect regression
    "compute_residuals", "compute_ser_statistics",
    "SER_posterior_e_loglik", "calculate_posterior_moments",
    "compute_kl", "get_ER2", "Eloglik", "loglik", "neg_loglik",
    # Model updates
    "update_fitted_values", "update_variance_components",
    "update_derived_quantities",
    # Output generation
    "get_scale_factors", "get_intercept", "get_fitted", "get_cs",
    "get_variable_names", "get_zscore", "cleanup_model"
  )
  for (g in generics) {
    expect_true(exists(g, mode = "function"), info = paste("Missing generic:", g))
  }
})

# ---- Method dispatch by data type ----

test_that("concrete methods exist for individual, ss, and rss_lambda data classes", {
  classes <- c("individual", "ss", "rss_lambda")
  key_generics <- c(
    "configure_data", "get_var_y", "initialize_susie_model",
    "initialize_fitted", "compute_residuals", "compute_ser_statistics",
    "SER_posterior_e_loglik", "calculate_posterior_moments",
    "get_ER2", "Eloglik", "loglik", "neg_loglik",
    "update_fitted_values", "update_variance_components",
    "get_scale_factors", "get_intercept", "get_fitted", "get_cs",
    "get_variable_names", "cleanup_model"
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
    "configure_data.default", "validate_prior.default",
    "track_ibss_fit.default", "compute_kl.default",
    "update_variance_components.default",
    "update_derived_quantities.default",
    "get_fitted.default", "get_zscore.default",
    "cleanup_model.default"
  )
  for (method in default_methods) {
    expect_true(exists(method, mode = "function"),
                info = paste("Missing default method:", method))
  }
})

# ---- Default method behavior ----

test_that("configure_data.default returns data unchanged", {
  data   <- structure(list(n = 50, p = 10), class = "test_class")
  params <- list(track_fit = FALSE)
  expect_identical(configure_data.default(data, params), data)
})

test_that("validate_prior.default returns TRUE", {
  data   <- structure(list(n = 50, p = 10), class = "test_class")
  params <- list()
  model  <- list(alpha = matrix(1/10, 3, 10), V = c(0.1, 0.2, 0.3), sigma2 = 1)
  expect_true(validate_prior.default(data, params, model))
})

test_that("update_derived_quantities.default returns model unchanged", {
  data   <- structure(list(n = 50, p = 10), class = "test_class")
  params <- list()
  model  <- list(alpha = matrix(1/10, 3, 10), V = c(0.1, 0.2, 0.3), sigma2 = 1)
  expect_identical(update_derived_quantities.default(data, params, model), model)
})

test_that("get_fitted.default and get_zscore.default return NULL", {
  data   <- structure(list(n = 50, p = 10), class = "test_class")
  params <- list()
  model  <- list(alpha = matrix(1/10, 3, 10))
  expect_null(get_fitted.default(data, params, model))
  expect_null(get_zscore.default(data, params, model))
})

test_that("track_ibss_fit.default stores iteration snapshots and respects track_fit=FALSE", {
  data    <- structure(list(), class = "test_class")
  model   <- list(alpha = matrix(1/10, 3, 10), V = c(0.1, 0.2, 0.3), sigma2 = 1)
  tracking <- list()

  # track_fit = TRUE: first snapshot is stored
  params1 <- list(track_fit = TRUE)
  result1 <- track_ibss_fit.default(data, params1, model, tracking, iter = 1, elbo = c(-Inf))
  expect_true(is.list(result1[[1]]))
  expect_s3_class(result1[[1]]$alpha,     "data.frame")
  expect_s3_class(result1[[1]]$effect,    "data.frame")
  expect_s3_class(result1[[1]]$iteration, "data.frame")
  expect_equal(result1[[1]]$iteration$sigma2, 1)

  # Second call appends a second snapshot
  result2 <- track_ibss_fit.default(data, params1, model, result1, iter = 2, elbo = c(-Inf, 100))
  expect_equal(length(result2), 2)
  expect_equal(result2[[2]]$iteration$iteration, 1)

  # track_fit = FALSE: tracking list stays empty
  params2 <- list(track_fit = FALSE)
  result3 <- track_ibss_fit.default(data, params2, model, list(), iter = 1, elbo = c(-Inf))
  expect_equal(length(result3), 0)
})

test_that("cleanup_model.default removes temporary fields and keeps core fields", {
  data  <- structure(list(), class = "test_class")
  model <- list(
    alpha            = matrix(1/10, 3, 10),
    mu               = matrix(0, 3, 10),
    sigma2           = 1,
    V                = c(0.1, 0.2, 0.3),
    null_weight      = 0,
    predictor_weights = rep(1/10, 10),
    residuals        = rnorm(50),
    fitted_without_l = rnorm(50),
    runtime          = list(prev_elbo = -100, prev_alpha = matrix(1/10, 3, 10),
                            prev_pip_diff = 0.01)
  )

  result <- cleanup_model.default(data, params = list(), model)

  expect_true("alpha"  %in% names(result))
  expect_true("mu"     %in% names(result))
  expect_true("sigma2" %in% names(result))
  expect_true("V"      %in% names(result))

  expect_false("null_weight"       %in% names(result))
  expect_false("predictor_weights" %in% names(result))
  expect_false("residuals"         %in% names(result))
  expect_false("fitted_without_l"  %in% names(result))
  expect_false("runtime"           %in% names(result))
})

# ---- Default method error dispatch (stop() branches) ----

test_that("required generics stop() with informative message for unsupported class", {
  unsupported <- structure(list(n = 50, p = 10), class = "unsupported_class")
  params   <- list()
  model    <- list(alpha = matrix(1/10, 5, 10), sigma2 = 1,
                   lbf_variable = matrix(0, 5, 10))
  ser_stats <- list(betahat = rnorm(10), shat2 = rep(1, 10))

  cases <- list(
    list(fn = "get_var_y.default",
         args = list(unsupported),
         pat = "get_var_y: no method for class 'unsupported_class'"),

    list(fn = "initialize_susie_model.default",
         args = list(unsupported, params),
         pat = "initialize_susie_model: no method for class 'unsupported_class'"),

    list(fn = "initialize_fitted.default",
         args = list(unsupported, matrix(0, 5, 10)),
         pat = "initialize_fitted: no method for class 'unsupported_class'"),

    list(fn = "compute_residuals.default",
         args = list(unsupported, params, model, 1L),
         pat = "compute_residuals: no method for class 'unsupported_class'"),

    list(fn = "compute_ser_statistics.default",
         args = list(unsupported, params,
                     modifyList(model, list(residuals = rnorm(50))), 1L),
         pat = "compute_ser_statistics: no method for class 'unsupported_class'"),

    list(fn = "SER_posterior_e_loglik.default",
         args = list(unsupported, params, model, 1L),
         pat = "SER_posterior_e_loglik: no method for class 'unsupported_class'"),

    list(fn = "calculate_posterior_moments.default",
         args = list(unsupported, params, model, 1.0),
         pat = "calculate_posterior_moments: no method for class 'unsupported_class'"),

    list(fn = "get_ER2.default",
         args = list(unsupported, model),
         pat = "get_ER2: no method for class 'unsupported_class'"),

    list(fn = "Eloglik.default",
         args = list(unsupported, model),
         pat = "Eloglik: no method for class 'unsupported_class'"),

    list(fn = "loglik.default",
         args = list(unsupported, params, model, 1.0, ser_stats),
         pat = "loglik: no method for class 'unsupported_class'"),

    list(fn = "neg_loglik.default",
         args = list(unsupported, params, model, 0.0, ser_stats),
         pat = "neg_loglik: no method for class 'unsupported_class'"),

    list(fn = "update_fitted_values.default",
         args = list(unsupported, params,
                     modifyList(model, list(mu = matrix(0, 5, 10))), 1L),
         pat = "update_fitted_values: no method for class 'unsupported_class'"),

    list(fn = "get_scale_factors.default",
         args = list(unsupported, params),
         pat = "get_scale_factors: no method for class 'unsupported_class'"),

    list(fn = "get_intercept.default",
         args = list(unsupported, params,
                     modifyList(model, list(mu = matrix(0, 5, 10)))),
         pat = "get_intercept: no method for class 'unsupported_class'"),

    list(fn = "get_cs.default",
         args = list(unsupported,
                     modifyList(params, list(coverage = 0.95, min_abs_corr = 0.5)),
                     model),
         pat = "get_cs: no method for class 'unsupported_class'"),

    list(fn = "get_variable_names.default",
         args = list(unsupported, model),
         pat = "get_variable_names: no method for class 'unsupported_class'")
  )

  for (case in cases) {
    # Resolve the unexported .default method from the namespace so this works
    # under both devtools::load_all() and installed-package testing (R CMD check).
    fn   <- getFromNamespace(case$fn, "susieR")
    args <- case$args
    pat  <- case$pat
    expect_error(
      do.call(fn, args),
      pat,
      info = paste("Expected stop() in", case$fn)
    )
  }
})

# ---- Verbose row generics ----

test_that("format_sigma2_summary.default returns sprintf '%.4f' of scalar sigma2", {
  model  <- list(sigma2 = 1.2345678)
  result <- format_sigma2_summary.default(model)
  expect_equal(result, sprintf("%.4f", 1.2345678))
  expect_type(result, "character")
  expect_length(result, 1L)
})

test_that("format_extra_diag.default returns empty string", {
  expect_identical(format_extra_diag.default(list()), "")
})

# ---- cleanup_extra_fields generic ----

test_that("cleanup_extra_fields.default returns character(0)", {
  expect_identical(cleanup_extra_fields.default(list()), character(0))
})

test_that("cleanup_model.default strips only standard temp fields, preserves others", {
  model <- list(
    null_weight      = 0.5,
    runtime          = list(prev_elbo = -Inf),
    fitted_without_l = NA,
    keep_me          = 42
  )
  out <- cleanup_model.default(list(), params = list(), model = model)
  expect_null(out$null_weight)
  expect_null(out$runtime)
  expect_null(out$fitted_without_l)
  expect_equal(out$keep_me, 42)
})

# ---- get_objective.default na.rm = TRUE branch ----

test_that("get_objective.default skips NA entries in KL when computing standard ELBO", {
  # Use a real individual data object so Eloglik.individual is dispatched.
  # Inject NA into model$KL to exercise the na.rm = TRUE code path.
  set.seed(7)
  setup <- setup_individual_data(n = 50, p = 20, L = 5)
  setup$params$max_iter <- 5
  setup$params$tol      <- 1e-2

  # Run one iteration to get a properly initialised model with fitted_without_l etc.
  model_warm <- suppressWarnings(
    susie_workhorse(setup$data, setup$params)
  )

  # Build a minimal model compatible with get_objective: needs the fields that
  # Eloglik.individual and the standard branch expect.
  # predictor_weights = attr(data$X, "d") is the per-predictor squared norm,
  # stripped by cleanup_model but required by get_ER2.individual.
  model_na <- list(
    alpha             = model_warm$alpha,
    mu                = model_warm$mu,
    mu2               = model_warm$mu2,
    V                 = model_warm$V,
    sigma2            = model_warm$sigma2,
    # Inject NA in two positions to exercise na.rm = TRUE
    KL                = replace(model_warm$KL, c(2, 4), NA_real_),
    lbf               = model_warm$lbf,
    lbf_variable      = model_warm$lbf_variable,
    Xr                = model_warm$fitted,
    predictor_weights = attr(setup$data$X, "d")
  )

  params_std <- modifyList(setup$params, list(
    unmappable_effects = "none",
    use_NIG            = FALSE
  ))

  obj_na <- get_objective(setup$data, params_std, model_na)

  # The same call with KL NAs replaced by zero should differ only by those KL terms
  model_zero_kl         <- model_na
  model_zero_kl$KL      <- replace(model_warm$KL, c(2, 4), 0)
  obj_zero              <- get_objective(setup$data, params_std, model_zero_kl)

  # Both should be finite (na.rm prevented NA propagation)
  expect_true(is.finite(obj_na))
  expect_true(is.finite(obj_zero))
  # The difference equals the sum of the zeroed-out KL terms
  expected_diff <- sum(model_warm$KL[c(2, 4)], na.rm = TRUE)
  expect_equal(obj_zero - obj_na, expected_diff, tolerance = 1e-8)
})
