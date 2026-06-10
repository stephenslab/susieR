context("S3 methods for individual data class")

# ---- Data Initialization & Configuration ----

test_that("configure_data.individual returns data when unmappable_effects='none'", {
  setup <- setup_individual_data()
  setup$params$unmappable_effects <- "none"

  result <- configure_data.individual(setup$data, setup$params)

  expect_s3_class(result, "individual")
})

test_that("get_var_y.individual equals var(y)", {
  setup <- setup_individual_data()

  var_y <- get_var_y.individual(setup$data)

  expect_equal(var_y, var(setup$data$y), tolerance = 1e-15)
})

# ---- Model Initialization & Setup ----

test_that("initialize_susie_model.individual predictor_weights equal attr(X, 'd')", {
  setup <- setup_individual_data()
  var_y <- var(setup$data$y)

  model <- initialize_susie_model.individual(setup$data, setup$params, var_y)

  expect_equal(model$predictor_weights, attr(setup$data$X, "d"), tolerance = 1e-15)
  expect_length(model$predictor_weights, setup$data$p)
})

test_that("initialize_fitted.individual produces Xr of length n", {
  setup <- setup_individual_data()

  mat_init <- list(alpha = setup$model$alpha, mu = setup$model$mu)
  fitted <- initialize_fitted.individual(setup$data, mat_init)

  expect_length(fitted$Xr, setup$data$n)
  expect_true(is.numeric(fitted$Xr))
})

test_that("validate_prior.individual returns a logical scalar", {
  setup <- setup_individual_data()
  result <- validate_prior.individual(setup$data, setup$params, setup$model)
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("track_ibss_fit.individual returns a list", {
  setup <- setup_individual_data()
  result <- track_ibss_fit.individual(setup$data, setup$params, setup$model,
                                      list(), 1, -100)
  expect_type(result, "list")
})

# ---- Single Effect Regression & ELBO ----

test_that("compute_residuals.individual produces correctly-named residual components", {
  setup <- setup_individual_data()

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l = 1)

  expect_length(model$residuals, setup$data$p)
  expect_length(model$raw_residuals, setup$data$n)
  expect_length(model$fitted_without_l, setup$data$n)
  expect_true(is.numeric(model$residual_variance) && model$residual_variance > 0)
})

test_that("compute_ser_statistics.individual produces positive shat2", {
  setup <- setup_individual_data()

  model      <- compute_residuals.individual(setup$data, setup$params, setup$model, l = 1)
  ser_stats  <- compute_ser_statistics.individual(setup$data, setup$params, model, l = 1)

  expect_length(ser_stats$betahat, setup$data$p)
  expect_length(ser_stats$shat2, setup$data$p)
  expect_true(all(ser_stats$shat2 > 0))
  expect_true("optim_init" %in% names(ser_stats))
  expect_true("optim_bounds" %in% names(ser_stats))
  expect_true("optim_scale" %in% names(ser_stats))
})

test_that("compute_ser_statistics.individual shat2_inflation branch multiplies shat2", {
  set.seed(305)
  setup  <- setup_individual_data(n = 80, p = 10)
  var_y  <- get_var_y.individual(setup$data)
  model  <- initialize_susie_model.individual(setup$data, setup$params, var_y)
  model$residuals        <- as.vector(crossprod(setup$data$X, setup$data$y) /
                                        model$predictor_weights)
  model$residual_variance <- var_y
  model$V                <- rep(0.2, setup$params$L)
  model$shat2_inflation  <- rep(1.2, setup$data$p)

  result         <- compute_ser_statistics.individual(setup$data, setup$params, model, 1)
  baseline_shat2 <- var_y / model$predictor_weights

  expect_equal(result$shat2, baseline_shat2 * 1.2, tolerance = 1e-8)
})

test_that("calculate_posterior_moments.individual mu2 >= mu^2 (positive variance)", {
  setup <- setup_individual_data()

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l = 1)
  model <- calculate_posterior_moments.individual(setup$data, setup$params, model, V = 1.0, l = 1)

  expect_length(model$mu[1, ],  setup$data$p)
  expect_length(model$mu2[1, ], setup$data$p)
  expect_true(all(model$mu2[1, ] >= model$mu[1, ]^2 - 1e-10))
})

test_that("calculate_posterior_moments.individual V=0 zeroes mu and mu2", {
  setup <- setup_individual_data()
  setup$params$use_NIG <- TRUE

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l = 1)
  model <- calculate_posterior_moments.individual(setup$data, setup$params, model, V = 0, l = 1)

  expect_equal(model$mu[1, ],  rep(0, setup$data$p))
  expect_equal(model$mu2[1, ], rep(0, setup$data$p))
})

test_that("calculate_posterior_moments.individual shat2_inflation branch returns finite mu", {
  set.seed(306)
  setup  <- setup_individual_data(n = 80, p = 10)
  var_y  <- get_var_y.individual(setup$data)
  model  <- initialize_susie_model.individual(setup$data, setup$params, var_y)
  model$residuals        <- as.vector(crossprod(setup$data$X, setup$data$y) /
                                        model$predictor_weights)
  model$residual_variance <- var_y
  model$shat2_inflation  <- rep(1.3, setup$data$p)

  result <- calculate_posterior_moments.individual(setup$data, setup$params, model, V = 0.2, l = 1)

  expect_true(all(is.finite(result$mu[1, ])))
  expect_true(all(is.finite(result$mu2[1, ])))
})

test_that("compute_kl.individual returns a finite scalar KL after a full SER step", {
  setup <- setup_individual_data()
  l <- 1

  # Run a full SER step so that lbf / alpha / mu are mutually consistent.
  model     <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  ser_stats <- compute_ser_statistics.individual(setup$data, setup$params, model, l)
  model     <- loglik.individual(setup$data, setup$params, model, V = 1.0, ser_stats, l)
  model     <- calculate_posterior_moments.individual(setup$data, setup$params, model, V = 1.0, l)
  model     <- compute_kl.individual(setup$data, setup$params, model, l)

  expect_length(model$KL[l], 1)
  expect_true(is.finite(model$KL[l]))
})

test_that("get_ER2.individual returns a non-negative scalar", {
  setup <- setup_individual_data()

  er2 <- get_ER2.individual(setup$data, setup$model)

  expect_length(er2, 1)
  expect_true(er2 >= 0)
})

test_that("Eloglik.individual returns a finite scalar", {
  setup <- setup_individual_data()

  e_loglik <- Eloglik.individual(setup$data, setup$model)

  expect_length(e_loglik, 1)
  expect_true(is.finite(e_loglik))
})

test_that("loglik.individual alpha sums to 1 and lbf is numeric", {
  setup <- setup_individual_data()

  model     <- compute_residuals.individual(setup$data, setup$params, setup$model, l = 1)
  ser_stats <- compute_ser_statistics.individual(setup$data, setup$params, model, l = 1)
  model     <- loglik.individual(setup$data, setup$params, model, V = 1.0, ser_stats, l = 1)

  expect_equal(sum(model$alpha[1, ]), 1, tolerance = 1e-10)
  expect_true(all(model$alpha[1, ] >= 0))
  expect_true(is.numeric(model$lbf[1]))
  expect_length(model$lbf_variable[1, ], setup$data$p)
})

test_that("neg_loglik.individual returns a finite scalar", {
  setup <- setup_individual_data()

  model     <- compute_residuals.individual(setup$data, setup$params, setup$model, l = 1)
  ser_stats <- compute_ser_statistics.individual(setup$data, setup$params, model, l = 1)
  neg_ll    <- neg_loglik.individual(setup$data, setup$params, model, log(1.0), ser_stats)

  expect_length(neg_ll, 1)
  expect_true(is.finite(neg_ll))
})

test_that("SER_posterior_e_loglik.individual returns a finite scalar", {
  set.seed(2)
  setup <- setup_individual_data()
  l <- 1

  setup$model$alpha[l, ] <- rep(1 / setup$data$p, setup$data$p)
  setup$model$mu[l, ]    <- rnorm(setup$data$p)
  setup$model$mu2[l, ]   <- setup$model$mu[l, ]^2 + 0.1

  model    <- compute_residuals.individual(setup$data, setup$params, setup$model, l)
  e_loglik <- SER_posterior_e_loglik.individual(setup$data, setup$params, model, l)

  expect_length(e_loglik, 1)
  expect_true(is.finite(e_loglik))
})

# ---- Model Updates & Fitting ----

test_that("update_fitted_values.individual produces Xr of length n", {
  setup <- setup_individual_data()

  model <- compute_residuals.individual(setup$data, setup$params, setup$model, l = 1)
  setup$model$fitted_without_l <- model$fitted_without_l
  updated <- update_fitted_values.individual(setup$data, setup$params, setup$model, l = 1)

  expect_length(updated$Xr, setup$data$n)
  expect_true(is.numeric(updated$Xr))
})

test_that("update_variance_components.individual returns a list", {
  setup  <- setup_individual_data()
  result <- update_variance_components.individual(setup$data, setup$params, setup$model)
  expect_type(result, "list")
})

test_that("update_derived_quantities.individual returns a list", {
  setup  <- setup_individual_data()
  result <- update_derived_quantities.individual(setup$data, setup$params, setup$model)
  expect_type(result, "list")
})

# ---- Variance Component Branches ----

test_that("update_variance_components.individual inf+MLE path returns positive sigma2", {
  set.seed(303)
  n <- 120; p <- 12; L <- 2
  X    <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[4] <- 2
  y    <- as.vector(X %*% beta + rnorm(n, sd = 0.5))

  obj        <- individual_data_constructor(X, y, L = L,
                                            unmappable_effects = "inf",
                                            estimate_residual_method = "MoM",
                                            convergence_method = "pip",
                                            max_iter = 5)
  data       <- obj$data
  params_mom <- obj$params
  var_y      <- get_var_y.individual(data)
  model      <- initialize_susie_model.individual(data, params_mom, var_y)

  for (i in seq_len(3)) {
    model <- ibss_fit(data, params_mom, model)
    model <- update_model_variance(data, params_mom, model)
  }

  params_mle                          <- params_mom
  params_mle$estimate_residual_method <- "MLE"

  result <- update_variance_components.individual(data, params_mle, model)

  expect_true(is.numeric(result$sigma2) && result$sigma2 > 0)
  expect_true("tau2"  %in% names(result))
  expect_true("theta" %in% names(result))
})

test_that("update_variance_components.individual inf+MoM path sets tau2", {
  set.seed(303)
  n <- 120; p <- 25
  X    <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[7] <- 2
  y    <- as.vector(X %*% beta + rnorm(n, sd = 0.5))

  fit <- suppressWarnings(susie(X, y, L = 3,
                                unmappable_effects = "inf",
                                estimate_residual_method = "MoM",
                                convergence_method = "pip",
                                max_iter = 10, verbose = FALSE))

  expect_s3_class(fit, "susie")
  expect_true(!is.null(fit$tau2) && is.numeric(fit$tau2))
})

test_that("update_variance_components.individual ash path preserves positive sigma2", {
  set.seed(304)
  n <- 120; p <- 25
  X    <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[c(2, 15)] <- c(1, -1)
  y    <- as.vector(X %*% beta + rnorm(n, sd = 0.5))

  fit <- suppressWarnings(susie(X, y, L = 3,
                                unmappable_effects = "ash",
                                max_iter = 10, verbose = FALSE))

  expect_s3_class(fit, "susie")
  expect_true(is.numeric(fit$sigma2) && fit$sigma2 > 0)
})

test_that("calculate_posterior_moments.individual NIG path with V>0 returns a susie object", {
  set.seed(301)
  n    <- 100; p <- 20
  X    <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[5] <- 2
  y    <- as.vector(X %*% beta + rnorm(n, sd = 0.5))

  fit <- suppressWarnings(susie(X, y, L = 1,
                                estimate_residual_method = "NIG",
                                alpha0 = 2, beta0 = 1,
                                max_iter = 10, verbose = FALSE))

  expect_s3_class(fit, "susie")
  expect_true(is.numeric(fit$mu) && all(is.finite(fit$mu)))
})

test_that("compute_ser_statistics.individual inf path (optim_bounds = c(0,1)) runs correctly", {
  set.seed(302)
  n    <- 120; p <- 15
  X    <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[c(3, 9)] <- c(1.5, -1.5)
  y    <- as.vector(X %*% beta + rnorm(n, sd = 0.5))

  fit <- suppressWarnings(susie(X, y, L = 3,
                                unmappable_effects = "inf",
                                estimate_residual_method = "MoM",
                                convergence_method = "pip",
                                max_iter = 10, verbose = FALSE))

  expect_s3_class(fit, "susie")
  expect_true(is.numeric(fit$sigma2) && fit$sigma2 > 0)
  expect_true(!is.null(fit$tau2))
  expect_true(!is.null(fit$theta))
})

# ---- Output Generation & Post-Processing ----

test_that("get_scale_factors.individual returns column scale factors matching attr(X, 'scaled:scale')", {
  setup  <- setup_individual_data()
  scales <- get_scale_factors.individual(setup$data, setup$params)

  expect_equal(scales, attr(setup$data$X, "scaled:scale"), tolerance = 1e-15)
  expect_length(scales, setup$data$p)
  expect_true(all(scales > 0))
})

test_that("get_intercept.individual returns a finite scalar when intercept=TRUE", {
  setup                  <- setup_individual_data()
  setup$params$intercept <- TRUE

  intercept <- get_intercept.individual(setup$data, setup$params, setup$model)

  expect_length(intercept, 1)
  expect_true(is.finite(intercept))
})

test_that("get_intercept.individual returns 0 when intercept=FALSE", {
  setup                  <- setup_individual_data()
  setup$params$intercept <- FALSE

  expect_equal(get_intercept.individual(setup$data, setup$params, setup$model), 0)
})

test_that("get_fitted.individual incorporates intercept when intercept=TRUE", {
  setup                  <- setup_individual_data()
  setup$params$intercept <- TRUE
  setup$data$mean_y      <- 5.0

  fitted <- get_fitted.individual(setup$data, setup$params, setup$model)

  expect_length(fitted, setup$data$n)
  expect_true(any(fitted != setup$model$Xr))
})

test_that("get_fitted.individual equals Xr when intercept=FALSE and mean_y=0", {
  setup                  <- setup_individual_data()
  setup$params$intercept <- FALSE
  setup$data$mean_y      <- 0

  fitted <- get_fitted.individual(setup$data, setup$params, setup$model)

  expect_equal(fitted, drop(setup$model$Xr))
})

test_that("get_cs.individual returns NULL when coverage is NULL", {
  setup                  <- setup_individual_data()
  setup$params$coverage  <- NULL

  expect_null(get_cs.individual(setup$data, setup$params, setup$model))
})

test_that("get_cs.individual returns NULL when min_abs_corr is NULL", {
  setup                     <- setup_individual_data()
  setup$params$min_abs_corr <- NULL

  expect_null(get_cs.individual(setup$data, setup$params, setup$model))
})

test_that("get_variable_names.individual assigns X column names to model matrices", {
  setup                  <- setup_individual_data()
  colnames(setup$data$X) <- paste0("var", seq_len(setup$data$p))
  setup$model$pip           <- rep(0.1, setup$data$p)
  setup$model$null_weight    <- NULL
  setup$model$alpha         <- matrix(0, 5, setup$data$p)
  setup$model$mu            <- matrix(0, 5, setup$data$p)
  setup$model$mu2           <- matrix(0, 5, setup$data$p)
  setup$model$lbf_variable  <- matrix(0, 5, setup$data$p)

  model_with_names <- get_variable_names.individual(setup$data, setup$model)

  expect_true(all(grepl("^var", colnames(model_with_names$alpha))))
  expect_true(all(grepl("^var", colnames(model_with_names$mu))))
  expect_true(all(grepl("^var", names(model_with_names$pip))))
})

test_that("get_zscore.individual returns p-length numeric vector when enabled", {
  setup                             <- setup_individual_data()
  setup$params$compute_univariate_zscore <- TRUE

  z <- get_zscore.individual(setup$data, setup$params, setup$model)

  expect_length(z, setup$data$p)
  expect_true(is.numeric(z))
})

test_that("get_zscore.individual respects null_weight (drops last column)", {
  setup                             <- setup_individual_data()
  setup$params$compute_univariate_zscore <- TRUE
  setup$model$null_weight           <- 0.1
  setup$data$X                      <- cbind(setup$data$X, 0)

  z <- get_zscore.individual(setup$data, setup$params, setup$model)

  expect_length(z, setup$data$p)
})

test_that("get_zscore.individual returns NULL when compute_univariate_zscore=FALSE", {
  setup                             <- setup_individual_data()
  setup$params$compute_univariate_zscore <- FALSE

  expect_null(get_zscore.individual(setup$data, setup$params, setup$model))
})

test_that("get_zscore.individual warns and still computes z-scores for sparse X", {
  setup                             <- setup_individual_data()
  setup$params$compute_univariate_zscore <- TRUE
  setup$data$X <- Matrix::Matrix(setup$data$X, sparse = TRUE)

  expect_message(
    z <- get_zscore.individual(setup$data, setup$params, setup$model),
    "Calculation of univariate regression z-scores is not implemented specifically"
  )
  expect_length(z, setup$data$p)
  expect_true(is.numeric(z))
})

test_that("cleanup_model.individual removes raw_residuals from the model", {
  setup                    <- setup_individual_data()
  setup$model$raw_residuals <- rnorm(setup$data$n)
  setup$model$residuals     <- rnorm(setup$data$p)

  cleaned <- cleanup_model.individual(setup$data, setup$params, setup$model)

  expect_false("raw_residuals" %in% names(cleaned))
})
