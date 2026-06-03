context("S3 methods for rss_lambda data class")

# =============================================================================
# DATA INITIALIZATION & CONFIGURATION
# =============================================================================

test_that("configure_data.rss_lambda returns configured data object", {
  dat <- setup_rss_lambda_data(seed = 1)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(
    L = 5,
    prior_variance = 0.2,
    residual_variance = 1.0
  )

  configured <- configure_data.rss_lambda(data, params)

  expect_s3_class(configured, "rss_lambda")
  expect_length(configured$z, dat$p)
  expect_equal(dim(configured$R), c(dat$p, dat$p))
  expect_equal(configured$lambda, dat$lambda)
})

test_that("get_var_y.rss_lambda returns 1", {
  dat <- setup_rss_lambda_data(seed = 2)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  var_y <- get_var_y.rss_lambda(result$data)

  expect_equal(var_y, 1)
  expect_type(var_y, "double")
  expect_length(var_y, 1)
})

# =============================================================================
# MODEL INITIALIZATION & SETUP
# =============================================================================

test_that("initialize_susie_model.rss_lambda creates valid model", {
  dat <- setup_rss_lambda_data(seed = 3)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(
    L = 5,
    prior_variance = 0.2,
    residual_variance = 1.0,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE
  )

  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)

  expect_type(model, "list")
  expect_equal(dim(model$alpha), c(5, dat$p))
  expect_equal(dim(model$mu), c(5, dat$p))
  expect_equal(dim(model$SinvRj), c(dat$p, dat$p))
  expect_length(model$RjSinvRj, dat$p)
  expect_equal(dim(model$mu2), c(5, dat$p))
})

test_that("validate_prior.rss_lambda delegates to default method", {
  dat <- setup_rss_lambda_data(seed = 21)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)

  result <- validate_prior.rss_lambda(data, params, model)

  expect_type(result, "logical")
})

test_that("track_ibss_fit.rss_lambda delegates to default method", {
  dat <- setup_rss_lambda_data(seed = 22)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, track_fit = TRUE,
                 scaled_prior_variance = 0.2, residual_variance = 1.0,
                 prior_weights = rep(1/dat$p, dat$p))
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)

  result <- track_ibss_fit.rss_lambda(data, params, model, list(), iter = 1, elbo = -100)

  expect_type(result, "list")
})

# =============================================================================
# SINGLE EFFECT REGRESSION & ELBO
# =============================================================================

test_that("initialize_fitted.rss_lambda creates Rz", {
  dat <- setup_rss_lambda_data(seed = 4)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data

  mat_init <- list(
    alpha = matrix(1/dat$p, nrow = 5, ncol = dat$p),
    mu = matrix(0, nrow = 5, ncol = dat$p)
  )

  fitted <- initialize_fitted.rss_lambda(data, mat_init)

  expect_type(fitted, "list")
  expect_true("Rz" %in% names(fitted))
  expect_length(fitted$Rz, dat$p)
  expect_type(fitted$Rz, "double")
})

test_that("compute_residuals.rss_lambda computes correct residuals", {
  dat <- setup_rss_lambda_data(seed = 5)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))

  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  expect_true("residuals" %in% names(model))
  expect_true("fitted_without_l" %in% names(model))
  expect_length(model$residuals, dat$p)
  expect_length(model$fitted_without_l, dat$p)
  expect_equal(model$residual_variance, 1)
})

test_that("compute_ser_statistics.rss_lambda computes shat2 and optim params", {
  dat <- setup_rss_lambda_data(seed = 6)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  ser_stats <- compute_ser_statistics.rss_lambda(data, params, model, l = 1)

  expect_type(ser_stats, "list")
  expect_true("shat2" %in% names(ser_stats))
  expect_true("optim_init" %in% names(ser_stats))
  expect_true("optim_bounds" %in% names(ser_stats))
  expect_true("optim_scale" %in% names(ser_stats))
  expect_length(ser_stats$shat2, dat$p)
  expect_true(all(ser_stats$shat2 > 0))
  expect_equal(ser_stats$optim_scale, "log")
})

test_that("SER_posterior_e_loglik.rss_lambda computes expected log-likelihood", {
  dat <- setup_rss_lambda_data(seed = 7)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  e_loglik <- SER_posterior_e_loglik.rss_lambda(data, params, model, l = 1)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
  expect_true(is.finite(e_loglik))
})

test_that("compute_kl.rss_lambda delegates to default method", {
  set.seed(23)
  dat <- setup_rss_lambda_data(seed = 23)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))

  l <- 1
  model$lbf <- rep(0, params$L)
  model$alpha[l, ] <- rep(1/dat$p, dat$p)
  model$mu[l, ] <- rnorm(dat$p, sd = 0.1)
  model$mu2[l, ] <- model$mu[l, ]^2 + 0.1

  model <- compute_residuals.rss_lambda(data, params, model, l)
  model <- compute_kl.rss_lambda(data, params, model, l)

  expect_type(model$KL[l], "double")
  expect_length(model$KL[l], 1)
  expect_true(is.finite(model$KL[l]))
})

test_that("calculate_posterior_moments.rss_lambda computes moments", {
  dat <- setup_rss_lambda_data(seed = 8)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  V <- 0.2
  l <- 1
  model <- calculate_posterior_moments.rss_lambda(data, params, model, V, l)

  expect_length(model$mu[l, ], dat$p)
  expect_length(model$mu2[l, ], dat$p)
  post_var <- model$mu2[l, ] - model$mu[l, ]^2
  expect_true(all(post_var > -1e-10))
  expect_equal(model$mu2[l, ], post_var + model$mu[l, ]^2, tolerance = 1e-15)
})

test_that("Eloglik.rss_lambda computes expected log-likelihood", {
  dat <- setup_rss_lambda_data(seed = 24)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model <- precompute_rss_lambda_terms(data, model)

  e_loglik <- Eloglik.rss_lambda(data, model)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
  expect_true(is.finite(e_loglik))
})

test_that("loglik.rss_lambda computes log Bayes factors with valid posterior", {
  dat <- setup_rss_lambda_data(seed = 25)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0,
                 prior_weights = rep(1/dat$p, dat$p))
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  V <- 0.2
  l <- 1
  ser_stats <- compute_ser_statistics.rss_lambda(data, params, model, l = l)
  model <- loglik.rss_lambda(data, params, model, V, ser_stats, l)

  expect_length(model$lbf_variable[l, ], dat$p)
  expect_length(model$alpha[l, ], dat$p)
  expect_true(all(model$alpha[l, ] >= 0))
  expect_equal(sum(model$alpha[l, ]), 1, tolerance = 1e-10)
  expect_true(is.numeric(model$lbf[l]))
  expect_true(is.finite(model$lbf[l]))
})

test_that("neg_loglik.rss_lambda returns finite negative log-likelihood", {
  dat <- setup_rss_lambda_data(seed = 26)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0,
                 prior_weights = rep(1/dat$p, dat$p))
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  V_param <- log(1.0)
  ser_stats <- compute_ser_statistics.rss_lambda(data, params, model, l = 1)
  neg_ll <- neg_loglik.rss_lambda(data, params, model, V_param, ser_stats)

  expect_type(neg_ll, "double")
  expect_length(neg_ll, 1)
  expect_true(is.finite(neg_ll))
})

test_that("get_ER2.rss_lambda returns non-negative finite scalar", {
  dat <- setup_rss_lambda_data(seed = 27)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model <- precompute_rss_lambda_terms(data, model)

  er2 <- get_ER2.rss_lambda(data, model)

  expect_type(er2, "double")
  expect_length(er2, 1)
  expect_true(er2 >= 0)
  expect_true(is.finite(er2))
})

# =============================================================================
# MODEL UPDATES & FITTING
# =============================================================================

test_that("update_fitted_values.rss_lambda updates Rz correctly", {
  dat <- setup_rss_lambda_data(seed = 28)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  model <- update_fitted_values.rss_lambda(data, params, model, l = 1)

  expect_true("Rz" %in% names(model))
  expect_length(model$Rz, dat$p)
  expect_type(model$Rz, "double")
})

test_that("update_variance_components.rss_lambda estimates sigma2 within bounds", {
  dat <- setup_rss_lambda_data(seed = 29)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0,
                 estimate_residual_variance = TRUE)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model <- precompute_rss_lambda_terms(data, model)

  variance_update <- update_variance_components.rss_lambda(data, params, model)

  expect_type(variance_update, "list")
  expect_true("sigma2" %in% names(variance_update))
  expect_type(variance_update$sigma2, "double")
  expect_length(variance_update$sigma2, 1)
  expect_true(variance_update$sigma2 > 0)
  expect_true(variance_update$sigma2 <= 1 - dat$lambda)
})

test_that("update_variance_components.rss_lambda returns empty list when estimate_residual_variance=FALSE", {
  set.seed(403)
  n <- 200; p <- 12
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  ss <- compute_suff_stat(X, y, standardize = TRUE)
  ur <- univariate_regression(X, y)
  z  <- with(ur, betahat / sebetahat)
  R  <- cov2cor(ss$XtX); R <- (R + t(R)) / 2

  ctor <- rss_lambda_constructor(z = z, R = R, n = n, L = 2, lambda = 0.3,
                                 estimate_residual_variance = FALSE)
  data   <- ctor$data
  params <- ctor$params
  var_y  <- get_var_y.rss_lambda(data)
  model  <- initialize_susie_model.rss_lambda(data, params, var_y)

  result <- update_variance_components.rss_lambda(data, params, model)

  expect_type(result, "list")
  expect_length(result, 0)
})

test_that("update_derived_quantities.rss_lambda updates SinvRj and RjSinvRj", {
  dat <- setup_rss_lambda_data(seed = 12)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$sigma2 <- 0.5

  updated_model <- update_derived_quantities.rss_lambda(data, params, model)

  expect_true("SinvRj" %in% names(updated_model))
  expect_true("RjSinvRj" %in% names(updated_model))
  expect_equal(dim(updated_model$SinvRj), c(dat$p, dat$p))
  expect_length(updated_model$RjSinvRj, dat$p)
  expect_true(all(updated_model$RjSinvRj > 0))
})

# =============================================================================
# OUTPUT GENERATION & POST-PROCESSING
# =============================================================================

test_that("get_scale_factors.rss_lambda returns vector of 1s", {
  dat <- setup_rss_lambda_data(seed = 13)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  scale_factors <- get_scale_factors.rss_lambda(result$data, list())

  expect_type(scale_factors, "double")
  expect_length(scale_factors, dat$p)
  expect_equal(scale_factors, rep(1, dat$p))
})

test_that("get_intercept.rss_lambda returns intercept_value from data", {
  dat <- setup_rss_lambda_data(seed = 14)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)

  intercept <- get_intercept.rss_lambda(data, params, model)

  expect_type(intercept, "double")
  expect_length(intercept, 1)
  expect_equal(intercept, data$intercept_value)
})

test_that("get_fitted.rss_lambda returns NULL", {
  dat <- setup_rss_lambda_data(seed = 30)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)

  expect_null(get_fitted.rss_lambda(data, params, model))
})

test_that("get_cs.rss_lambda returns NULL when coverage or min_abs_corr is NULL", {
  for (null_param in c("coverage", "min_abs_corr")) {
    dat <- setup_rss_lambda_data(seed = 31)

    result <- rss_lambda_constructor(
      z = dat$z,
      R = dat$R,
      lambda = dat$lambda,
      n = dat$n
    )

    data <- result$data
    params <- list(L = 5, coverage = 0.95, min_abs_corr = 0.5,
                   scaled_prior_variance = 0.2, residual_variance = 1.0,
                   prior_weights = rep(1/dat$p, dat$p))
    params[[null_param]] <- NULL

    var_y <- get_var_y.rss_lambda(data)
    model <- initialize_susie_model.rss_lambda(data, params, var_y)

    expect_null(get_cs.rss_lambda(data, params, model),
                info = paste(null_param, "= NULL should return NULL CS"))
  }
})

test_that("get_cs.rss_lambda uses correlation from R matrix", {
  dat <- setup_rss_lambda_data(seed = 33)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, coverage = 0.95, min_abs_corr = 0.5, n_purity = 100,
                 scaled_prior_variance = 0.2, residual_variance = 1.0,
                 prior_weights = rep(1/dat$p, dat$p))
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$alpha[1, 1] <- 0.95
  model$alpha[1, -1] <- 0.05 / (dat$p - 1)

  cs <- get_cs.rss_lambda(data, params, model)

  expect_true(is.null(cs) || is.list(cs))
})

test_that("get_variable_names.rss_lambda assigns variable names to model", {
  dat <- setup_rss_lambda_data(seed = 34)

  z_named <- dat$z
  names(z_named) <- paste0("var", 1:dat$p)

  result <- rss_lambda_constructor(
    z = z_named,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0,
                 prior_weights = rep(1/dat$p, dat$p))
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$pip <- rep(0.1, dat$p)
  model$null_weight <- NULL
  model$alpha <- matrix(0, 5, dat$p)
  model$mu <- matrix(0, 5, dat$p)
  model$mu2 <- matrix(0, 5, dat$p)
  model$lbf_variable <- matrix(0, 5, dat$p)

  model_with_names <- get_variable_names.rss_lambda(data, model)

  expect_true(all(grepl("var", colnames(model_with_names$alpha))))
  expect_true(all(grepl("var", colnames(model_with_names$mu))))
  expect_true(all(grepl("var", colnames(model_with_names$mu2))))
  expect_true(all(grepl("var", names(model_with_names$pip))))
})

test_that("get_zscore.rss_lambda returns NULL", {
  dat <- setup_rss_lambda_data(seed = 35)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, compute_univariate_zscore = TRUE,
                 scaled_prior_variance = 0.2, residual_variance = 1.0,
                 prior_weights = rep(1/dat$p, dat$p))
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)

  expect_null(get_zscore.rss_lambda(data, params, model))
})

test_that("cleanup_model.rss_lambda removes temporary fields and retains essential ones", {
  dat <- setup_rss_lambda_data(seed = 38)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data
  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0)
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- rep(0, dat$p)
  model$Z <- matrix(0, 5, dat$p)
  model$zbar <- rep(0, dat$p)
  model$diag_postb2 <- rep(0, dat$p)

  cleaned <- cleanup_model.rss_lambda(data, params, model)

  expect_false("SinvRj" %in% names(cleaned))
  expect_false("RjSinvRj" %in% names(cleaned))
  expect_false("Rz" %in% names(cleaned))
  expect_false("Z" %in% names(cleaned))
  expect_false("zbar" %in% names(cleaned))
  expect_false("diag_postb2" %in% names(cleaned))
  expect_true("alpha" %in% names(cleaned))
  expect_true("mu" %in% names(cleaned))
  expect_true("mu2" %in% names(cleaned))
})

# =============================================================================
# FINITE-REFERENCE R INFLATION TESTS
# =============================================================================

test_that("compute_ser_statistics.rss_lambda returns finite betahat", {
  dat <- setup_rss_lambda_data(seed = 40)

  data <- dat$data
  params <- dat$params
  model <- dat$model
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  ser_stats <- compute_ser_statistics.rss_lambda(data, params, model, l = 1)

  expect_true("betahat" %in% names(ser_stats))
  expect_length(ser_stats$betahat, dat$p)
  expect_true(all(is.finite(ser_stats$betahat)))
})

test_that("compute_residuals.rss_lambda does not set shat2_inflation when R_finite_B is absent", {
  dat <- setup_rss_lambda_data(seed = 42)

  data <- dat$data
  params <- dat$params
  model <- dat$model
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))

  expect_null(data$R_finite_B)
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)
  expect_null(model$shat2_inflation)
})

# =============================================================================
# R vs X INPUT PATH AGREEMENT
# =============================================================================

test_that("R and X input paths produce numerically identical results", {
  set.seed(50)
  p <- 30
  n <- 500
  X_full <- matrix(rnorm(n * p), n, p)
  X_full <- scale(X_full, center = TRUE, scale = TRUE)
  y <- X_full[, 1] * 0.5 + rnorm(n)
  input_ss <- compute_suff_stat(X_full, y, standardize = TRUE)
  R <- cov2cor(input_ss$XtX)
  R <- (R + t(R)) / 2
  ss <- univariate_regression(X_full, y)
  z <- ss$betahat / ss$sebetahat

  res_R <- rss_lambda_constructor(z = z, R = R, lambda = 0.1, n = n)
  res_X <- rss_lambda_constructor(z = z, X = X_full, lambda = 0.1, n = n)

  expect_equal(res_R$data$eigen_R$values, res_X$data$eigen_R$values, tolerance = 1e-6)

  var_y_R <- get_var_y.rss_lambda(res_R$data)
  model_R <- initialize_susie_model.rss_lambda(res_R$data, res_R$params, var_y_R)
  model_R$Rz <- as.vector(R %*% colSums(model_R$alpha * model_R$mu))

  var_y_X <- get_var_y.rss_lambda(res_X$data)
  model_X <- initialize_susie_model.rss_lambda(res_X$data, res_X$params, var_y_X)
  model_X$Rz <- as.vector(compute_Rv(res_X$data, colSums(model_X$alpha * model_X$mu)))

  expect_equal(model_R$RjSinvRj, model_X$RjSinvRj, tolerance = 1e-6)

  model_R <- compute_residuals.rss_lambda(res_R$data, res_R$params, model_R, l = 1)
  model_X <- compute_residuals.rss_lambda(res_X$data, res_X$params, model_X, l = 1)
  expect_equal(model_R$residuals, model_X$residuals, tolerance = 1e-6)

  stats_R <- compute_ser_statistics.rss_lambda(res_R$data, res_R$params, model_R, l = 1)
  stats_X <- compute_ser_statistics.rss_lambda(res_X$data, res_X$params, model_X, l = 1)
  expect_equal(stats_R$betahat, stats_X$betahat, tolerance = 1e-6)
  expect_equal(stats_R$shat2, stats_X$shat2, tolerance = 1e-6)
})

# =============================================================================
# END-TO-END susie_rss_lambda
# =============================================================================

test_that("susie_rss_lambda defaults max_iter to 50 with a hint", {
  set.seed(510)
  p <- 10
  n <- 200
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  expect_message(
    obj <- susie_rss_lambda(z = z, R = R, n = n, L = 2,
                            lambda = 0.1, init_only = TRUE,
                            verbose = FALSE),
    "Setting max_iter = 50 for the SuSiE RSS-lambda model"
  )
  expect_equal(obj$params$max_iter, 50)

  obj2 <- susie_rss_lambda(z = z, R = R, n = n, L = 2,
                           lambda = 0.1, max_iter = 7,
                           init_only = TRUE, verbose = FALSE)
  expect_equal(obj2$params$max_iter, 7)
})

test_that("susie_rss_lambda with lambda > 0 converges and recovers signal", {
  set.seed(51)
  p <- 50
  n <- 2000
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = TRUE)
  beta <- rep(0, p)
  beta[1] <- 0.5
  beta[10] <- -0.3
  y <- drop(X %*% beta + rnorm(n))
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  R <- cov2cor(input_ss$XtX)
  R <- (R + t(R)) / 2
  ss <- univariate_regression(X, y)
  z <- ss$betahat / ss$sebetahat

  fit <- suppressWarnings(
    susie_rss_lambda(z = z, R = R, lambda = 0.1, n = n, L = 5,
                     max_iter = 50, verbose = FALSE)
  )
  expect_true(fit$converged)
  expect_true(is.finite(fit$elbo[length(fit$elbo)]))
  expect_true(fit$pip[1] > 0.5)
})

test_that("susie_rss_lambda rejects unsupported arguments", {
  set.seed(511)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  expect_error(
    susie_rss_lambda(z = z, R = R, n = n, L = 3, lambda = 0.1,
                     R_finite = 5000, max_iter = 2, verbose = FALSE),
    "unused argument"
  )
  expect_error(
    susie_rss_lambda(z = z, R = R, n = n, L = 3, lambda = 0.1,
                     R_mismatch = "eb_no_init", max_iter = 2, verbose = FALSE),
    "unused argument"
  )
  expect_error(
    susie_rss_lambda(z = z, R = R, n = n, L = 3, lambda = 0.1,
                     R_mismatch = "eb", max_iter = 2, verbose = FALSE),
    "unused argument"
  )
  expect_error(
    susie_rss_lambda(z = z, X = list(X, X), n = n, L = 3, lambda = 0.1,
                     max_iter = 2, verbose = FALSE),
    "single X matrix"
  )
  expect_error(
    susie_rss_lambda(z = z, R = list(R, R), n = n, L = 3, lambda = 0.1,
                     max_iter = 2, verbose = FALSE),
    "single R matrix"
  )
  expect_error(
    susie_rss_lambda(z = z, R = R, n = n, L = 3, lambda = 0.1,
                     estimate_residual_method = "MoM",
                     max_iter = 2, verbose = FALSE),
    "MLE"
  )
})

test_that("susie_rss_lambda with estimate_residual_variance=TRUE bounds sigma2 by 1-lambda", {
  set.seed(401)
  n <- 300; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[c(4, 12)] <- c(1, -1)
  y <- as.vector(X %*% beta + rnorm(n, sd = 0.8))

  ss <- compute_suff_stat(X, y, standardize = TRUE)
  ur <- univariate_regression(X, y)
  z  <- with(ur, betahat / sebetahat)
  R  <- cov2cor(ss$XtX)
  R  <- (R + t(R)) / 2

  lambda_val <- 0.3
  fit <- suppressWarnings(
    susie_rss_lambda(z = z, R = R, n = n, L = 3,
                     lambda = lambda_val,
                     estimate_residual_variance = TRUE,
                     estimate_residual_method = "MLE",
                     max_iter = 15, verbose = FALSE)
  )

  expect_s3_class(fit, "susie")
  expect_true(is.numeric(fit$sigma2) && fit$sigma2 > 0)
  expect_true(fit$sigma2 <= 1 - lambda_val + 1e-8)
})

test_that("susie_rss_lambda sigma2 bound holds for lambda = 0.5", {
  set.seed(402)
  n <- 250; p <- 15
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  ss <- compute_suff_stat(X, y, standardize = TRUE)
  ur <- univariate_regression(X, y)
  z  <- with(ur, betahat / sebetahat)
  R  <- cov2cor(ss$XtX); R <- (R + t(R)) / 2

  lambda_val <- 0.5
  fit <- suppressWarnings(
    susie_rss_lambda(z = z, R = R, n = n, L = 2,
                     lambda = lambda_val,
                     estimate_residual_variance = TRUE,
                     max_iter = 10, verbose = FALSE)
  )

  expect_true(fit$sigma2 >= 0)
  expect_true(fit$sigma2 <= 1 - lambda_val + 1e-8)
})

# =============================================================================
# SS-PATH R_MISMATCH REGRESSION TESTS
# =============================================================================

test_that("R_mismatch works with and without finite-reference input", {
  set.seed(511)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit_inf <- susie_rss(z = z, R = R, n = n, L = 3, R_mismatch = "eb_no_init",
                       max_iter = 2, verbose = FALSE)
  expect_equal(fit_inf$R_finite_diagnostics$B, Inf)
  expect_length(fit_inf$R_finite_diagnostics$lambda_bias, 1)
  expect_length(fit_inf$R_finite_diagnostics$B_corrected, 1)
  expect_true(fit_inf$R_finite_diagnostics$lambda_bias >= 0)
  expect_equal(fit_inf$R_finite_diagnostics$B_corrected,
               1 / fit_inf$R_finite_diagnostics$lambda_bias,
               tolerance = 1e-12)

  fit_false <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = FALSE,
                         R_mismatch = "eb_no_init", max_iter = 2, verbose = FALSE)
  expect_equal(fit_false$R_finite_diagnostics$B, Inf)

  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, R_finite = 10000, R_mismatch = "mle",
              max_iter = 2, verbose = FALSE),
    "should be one of"
  )

  expect_message(
    fit_warn <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = 10000,
                          R_mismatch = "eb_no_init", estimate_residual_variance = TRUE,
                          max_iter = 2, verbose = FALSE),
    "incompatible with"
  )

  fit <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = 10000,
                   R_mismatch = "eb_no_init", max_iter = 2, verbose = FALSE)
  expect_length(fit$R_finite_diagnostics$lambda_bias, 1)
  expect_length(fit$R_finite_diagnostics$B_corrected, 1)
  expect_true(fit$R_finite_diagnostics$lambda_bias >= 0)
  R_finiteB <- fit$R_finite_diagnostics$B
  if (fit$R_finite_diagnostics$lambda_bias > 0) {
    expect_true(fit$R_finite_diagnostics$B_corrected < R_finiteB)
  } else {
    expect_equal(fit$R_finite_diagnostics$B_corrected, R_finiteB)
  }
})

test_that("R_mismatch = 'none' is identical to no-R_mismatch call", {
  set.seed(913)
  p <- 25
  n <- 2000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)
  z[3] <- 4

  fit_none <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
                        R_mismatch = "none", max_iter = 5, verbose = FALSE)
  fit_default <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
                           max_iter = 5, verbose = FALSE)
  expect_equal(fit_none$pip, fit_default$pip, tolerance = 1e-12)
  expect_equal(fit_none$alpha, fit_default$alpha, tolerance = 1e-12)
  expect_null(fit_none$lambda_bias)
})

test_that("Fisher SE zero-mask sends near-boundary estimates to 0 under null", {
  set.seed(7)
  p <- 50
  n <- 5000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)
  fit <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = 10000,
                   R_mismatch = "eb_no_init", max_iter = 5, verbose = FALSE)
  lb <- fit$R_finite_diagnostics$lambda_bias
  expect_true(all(lb == 0 | lb > 1e-6),
              info = "Fisher zero-mask must leave no values in the (0, 1e-6) gap")
})

test_that("In-sample LD yields lambda_bias = 0", {
  set.seed(2024)
  p <- 30
  n <- 4000
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = TRUE)
  beta <- rep(0, p); beta[5] <- 0.4
  y <- drop(X %*% beta + rnorm(n))
  ss <- compute_suff_stat(X, y, standardize = TRUE)
  R <- cov2cor(ss$XtX)
  z <- ss$XtX %*% beta / sqrt(diag(ss$XtX)) + rnorm(p)
  z <- as.numeric(z)
  fit <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
                   R_mismatch = "eb_no_init", max_iter = 8, verbose = FALSE)
  expect_true(all(fit$R_finite_diagnostics$lambda_bias == 0),
              info = "In-sample LD must produce lambda_bias = 0")
})

test_that("In-sample LD with multiple sparse signals does not inflate lambda_bias", {
  set.seed(44)
  n <- 1000
  p <- 120
  rho <- 0.95
  Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
  X <- scale(X, center = TRUE, scale = TRUE)
  beta <- rep(0, p)
  causal <- c(20, 60, 100)
  beta[causal] <- c(0.18, -0.20, 0.22)
  y <- drop(X %*% beta + rnorm(n))
  z <- calc_z(X, y, center = TRUE, scale = FALSE)

  fit <- suppressWarnings(
    susie_rss(z = z, X = X, n = n, L = 6, R_finite = TRUE,
              R_mismatch = "eb_no_init", max_iter = 50, verbose = FALSE)
  )

  expect_true(max(fit$R_finite_diagnostics$lambda_bias) < 0.01,
              info = "In-sample LD should not estimate large population mismatch")
  expect_gt(max(fit$pip[causal]), 0.5)
})

test_that("R_mismatch = 'mle' is rejected at susie_rss and summary_stats_constructor", {
  set.seed(31)
  p <- 20; n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X); z <- rnorm(p)
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, R_finite = 10000, R_mismatch = "mle",
              max_iter = 1, verbose = FALSE),
    "should be one of"
  )
  expect_error(
    summary_stats_constructor(z = z, R = R, n = n, L = 3, R_finite = 10000,
                              R_mismatch = "mle"),
    "should be one of"
  )
})

test_that("Large R_finite limit reduces to pure-drift estimator", {
  set.seed(11)
  p <- 30; n <- 4000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  beta <- rep(0, p); beta[1] <- 0.6
  z <- as.numeric(R %*% beta * sqrt(n) + rnorm(p))
  fit <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = 1e12,
                   R_mismatch = "eb_no_init", max_iter = 5, verbose = FALSE)
  lb <- fit$R_finite_diagnostics$lambda_bias
  bc <- fit$R_finite_diagnostics$B_corrected
  expect_false(is.null(lb))
  expect_true(all(is.finite(lb)) && all(lb >= 0))
  # Where a drift correction is active, B_corrected -> 1/lambda_bias as R_finite -> Inf.
  active <- lb > 0
  if (any(active)) {
    expect_equal(bc[active], 1 / lb[active], tolerance = 1e-6,
                 info = "B_corrected -> 1/lambda_bias as R_finite -> Inf")
  }
})

test_that("tau_j^2 is monotone non-decreasing in lambda_bias", {
  s     <- c(0.5, 1.5, 3.0, 0.0)
  sigma2 <- 1.2
  B     <- 1000
  tau2 <- function(lambda) sigma2 + (1 / B + lambda) * s
  expect_true(all(tau2(0.05) >= tau2(0)))
  expect_true(all(tau2(0.5)  >= tau2(0.05)))
  expect_equal(tau2(0)[s == 0], rep(sigma2, sum(s == 0)))
})

test_that("Wakefield ABF and signal-squared BF form agree when there is no inflation", {
  dat <- setup_rss_lambda_data(seed = 53)
  data <- dat$data
  params <- dat$params
  model <- dat$model
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  ser_stats <- compute_ser_statistics.rss_lambda(data, params, model, l = 1)

  V <- 0.2
  shat2 <- pmax(ser_stats$shat2, .Machine$double.eps)
  lbf_wakefield <- -0.5 * log(1 + V / shat2) +
    0.5 * ser_stats$betahat^2 * V / (shat2 * (V + shat2))

  signal <- as.vector(crossprod(model$SinvRj, model$residuals))
  RjSinvRj <- model$RjSinvRj
  lbf_original <- -0.5 * log(1 + V * RjSinvRj) +
    0.5 * V * signal^2 / (1 + V * RjSinvRj)

  expect_equal(lbf_wakefield, lbf_original, tolerance = 1e-8)
})

# =============================================================================
# SS vs RSS-LAMBDA CROSS-PATH AGREEMENT
# =============================================================================

test_that("SS and RSS-lambda paths agree with small lambda", {
  set.seed(200)
  p <- 50; n <- 2000
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = TRUE)
  beta <- rep(0, p)
  beta[1] <- 0.5; beta[10] <- -0.3
  y <- drop(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  R <- cov2cor(input_ss$XtX); R <- (R + t(R)) / 2
  ss <- univariate_regression(X, y)
  z <- ss$betahat / ss$sebetahat

  fit_ss <- suppressWarnings(
    susie_rss(z = z, R = R, n = n, L = 5, max_iter = 100, verbose = FALSE)
  )
  fit_rss <- suppressWarnings(
    susie_rss_lambda(z = z, R = R, n = n, L = 5, lambda = 1e-6,
                     max_iter = 100, verbose = FALSE)
  )

  expect_true(fit_ss$converged)
  expect_true(fit_rss$converged)
  expect_equal(fit_ss$alpha, fit_rss$alpha, tolerance = 1e-4)
  expect_equal(fit_ss$pip, fit_rss$pip, tolerance = 1e-4)
})

# =============================================================================
# MULTI-PANEL HELPER UTILITIES
# =============================================================================

test_that("form_X_meta combines panels correctly", {
  set.seed(42)
  p <- 10
  X1 <- matrix(rnorm(50 * p), 50, p)
  X2 <- matrix(rnorm(30 * p), 30, p)
  omega <- c(0.6, 0.4)

  X_meta <- form_X_meta(list(X1, X2), omega)

  expect_equal(nrow(X_meta), 80)
  expect_equal(ncol(X_meta), p)
  expect_equal(X_meta[1:50, ], sqrt(0.6) * X1)
  expect_equal(X_meta[51:80, ], sqrt(0.4) * X2)
})

test_that("eigen_from_X recovers eigendecomposition of X'X", {
  set.seed(43)
  p <- 20
  X <- matrix(rnorm(100 * p), 100, p)
  R <- crossprod(X)
  eigen_R_direct <- eigen(R, symmetric = TRUE)

  eigen_R_svd <- eigen_from_X(X, p)

  expect_equal(eigen_R_svd$values, eigen_R_direct$values, tolerance = 1e-10)
  for (j in seq_len(p)) {
    inner <- abs(sum(eigen_R_svd$vectors[, j] * eigen_R_direct$vectors[, j]))
    expect_gt(inner, 0.99)
  }
})

test_that("eval_omega_eloglik_reduced matches pure R reference", {
  set.seed(44)
  p <- 50
  K <- 2

  X1 <- matrix(rnorm(15 * p), 15, p)
  X2 <- matrix(rnorm(10 * p), 10, p)
  X_list <- list(X1, X2)
  panel_R <- list(crossprod(X1), crossprod(X2))

  z <- rnorm(p)
  zbar <- rnorm(p) * 0.1
  diag_postb2 <- abs(rnorm(p)) * 0.01
  L <- 3
  Z <- matrix(rnorm(L * p) * 0.05, L, p)
  sigma2 <- 0.9
  lambda <- 0.01
  omega <- c(0.7, 0.3)

  val_R <- susieR:::eval_omega_eloglik_R(panel_R, omega, z, zbar, diag_postb2,
                                          Z, sigma2, lambda, K, p)

  cache <- susieR:::precompute_omega_cache(X_list, z)
  iter_cache <- susieR:::precompute_omega_iteration(cache, zbar, diag_postb2, Z)
  val_reduced <- susieR:::eval_omega_eloglik_reduced(cache, omega, iter_cache,
                                                      sigma2, lambda, K, p)

  expect_equal(val_R, val_reduced, tolerance = 1e-6)
})

test_that("eval_omega_eloglik is concave in omega", {
  set.seed(45)
  p <- 20
  K <- 2

  X1 <- matrix(rnorm(60 * p), 60, p)
  X2 <- matrix(rnorm(50 * p), 50, p)
  panel_R <- list(crossprod(X1) / 60, crossprod(X2) / 50)

  z <- rnorm(p)
  zbar <- rnorm(p) * 0.1
  diag_postb2 <- abs(rnorm(p)) * 0.01
  Z <- matrix(rnorm(2 * p) * 0.05, 2, p)

  eloglik <- function(w1) {
    susieR:::eval_omega_eloglik_R(panel_R, c(w1, 1 - w1), z, zbar,
                                   diag_postb2, Z, 0.9, 0.01, K, p)
  }

  vals <- sapply(seq(0, 1, 0.1), eloglik)
  for (i in 1:(length(vals) - 2)) {
    expect_gte(vals[i + 1], (vals[i] + vals[i + 2]) / 2 - 1e-8)
  }
})

test_that("accessor helpers fall through for single panel", {
  dat <- setup_rss_lambda_data(seed = 50)
  model <- dat$model

  expect_null(model$eigen_R)
  eigen_R <- get_eigen_R(dat$data, model)
  expect_equal(eigen_R$values, dat$data$eigen_R$values)

  expect_null(model$Vtz)
  Vtz <- get_Vtz(dat$data, model)
  expect_equal(Vtz, dat$data$Vtz)
})

test_that(".omega_tol has expected fields with positive values", {
  tol <- susieR:::.omega_tol
  expect_true(is.list(tol))
  expect_true("convergence" %in% names(tol))
  expect_true("grid_spacing" %in% names(tol))
  expect_true("fw_stop" %in% names(tol))
  expect_true("fw_max_iter" %in% names(tol))
  expect_true(tol$convergence > 0)
  expect_true(tol$grid_spacing > 0 && tol$grid_spacing < 1)
  expect_true(tol$fw_stop > 0)
  expect_true(tol$fw_max_iter >= 1L)
})

test_that("eigen_from_reduced recovers full eigendecomposition", {
  set.seed(55)
  p <- 30; B1 <- 40; B2 <- 35
  X1 <- matrix(rnorm(B1 * p), B1, p)
  X2 <- matrix(rnorm(B2 * p), B2, p)
  X_list <- lapply(list(X1, X2), susieR:::standardize_X)
  z <- rnorm(p)

  cache <- susieR:::precompute_omega_cache(X_list, z)

  omega <- c(0.7, 0.3)
  eig_reduced <- susieR:::eigen_from_reduced(cache, omega, K = 2, p = p)

  R_omega <- omega[1] * crossprod(X_list[[1]]) + omega[2] * crossprod(X_list[[2]])
  R_omega <- 0.5 * (R_omega + t(R_omega))
  eig_direct <- eigen(R_omega, symmetric = TRUE)

  r <- cache$r
  expect_equal(eig_reduced$values[1:r], eig_direct$values[1:r], tolerance = 1e-8)

  overlap <- abs(crossprod(eig_reduced$vectors[, 1:r], eig_direct$vectors[, 1:r]))
  expect_true(all(apply(overlap, 1, max) > 1 - 1e-8))
})

test_that("optimize_omega handles vertex optimum", {
  set.seed(57)
  p <- 25; B <- 100

  X1 <- matrix(rnorm(B * p), B, p)
  X_list <- lapply(list(X1, matrix(rnorm(B * p), B, p)), susieR:::standardize_X)
  z <- rnorm(p)

  R1 <- crossprod(X_list[[1]])
  R2 <- diag(p)
  panel_R <- list(R1, R2)

  zbar <- rnorm(p) * 0.1
  diag_postb2 <- abs(rnorm(p)) * 0.01
  Z <- matrix(rnorm(2 * p) * 0.05, 2, p)

  eval_fn <- function(omega_vec) {
    susieR:::eval_omega_eloglik_R(panel_R, omega_vec, z, zbar,
                                   diag_postb2, Z, 0.9, 0.1, 2, p)
  }

  result <- susieR:::optimize_omega(eval_fn, c(0.5, 0.5), K = 2)

  expect_equal(sum(result$omega), 1, tolerance = 1e-10)
  expect_true(all(result$omega >= -1e-10))
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("rss_lambda handles lambda near boundary: lambda = 0 and lambda near 1", {
  set.seed(601)
  p <- 20; n <- 500
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p); z[1] <- 4

  fit0 <- suppressWarnings(
    susie_rss_lambda(z = z, R = R, n = n, L = 3, lambda = 0,
                     max_iter = 10, verbose = FALSE)
  )
  expect_s3_class(fit0, "susie")
  expect_true(all(is.finite(fit0$pip)))

  fit_hi <- suppressWarnings(
    susie_rss_lambda(z = z, R = R, n = n, L = 3, lambda = 0.99,
                     max_iter = 10, verbose = FALSE)
  )
  expect_s3_class(fit_hi, "susie")
  expect_true(all(is.finite(fit_hi$pip)))
})

test_that("rss_lambda handles p = 1", {
  set.seed(602)
  n <- 200
  x <- rnorm(n)
  y <- x * 2 + rnorm(n)
  ss <- univariate_regression(matrix(x, ncol = 1), y)
  z <- ss$betahat / ss$sebetahat
  R <- matrix(1, 1, 1)

  fit <- suppressWarnings(
    susie_rss_lambda(z = z, R = R, n = n, L = 1, lambda = 0.1,
                     max_iter = 10, verbose = FALSE)
  )
  expect_s3_class(fit, "susie")
  expect_length(fit$pip, 1)
  expect_true(fit$pip[1] > 0.5)
})
