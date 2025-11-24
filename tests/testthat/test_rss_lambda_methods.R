context("S3 methods for rss_lambda data class")

# =============================================================================
# DATA INITIALIZATION & CONFIGURATION
# =============================================================================

test_that("configure_data.rss_lambda returns configured data object", {
  dat <- setup_rss_lambda_data(seed = 1)

  # Create rss_lambda data object
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
  expect_true(!is.null(configured$z))
  expect_true(!is.null(configured$R))
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

  data <- result$data

  var_y <- get_var_y.rss_lambda(data)

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
  expect_true(!is.null(model$alpha))
  expect_true(!is.null(model$mu))
  expect_true(!is.null(model$mu2))
  expect_true(!is.null(model$SinvRj))
  expect_true(!is.null(model$RjSinvRj))

  # Check dimensions
  expect_equal(dim(model$alpha), c(5, dat$p))
  expect_equal(dim(model$mu), c(5, dat$p))
  expect_equal(dim(model$SinvRj), c(dat$p, dat$p))
  expect_length(model$RjSinvRj, dat$p)
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
  params <- list(
    L = 5,
    prior_variance = 0.2,
    residual_variance = 1.0
  )
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

  tracking <- list()
  iter <- 1
  elbo <- -100

  result <- track_ibss_fit.rss_lambda(data, params, model, tracking, iter, elbo)

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

  # Create minimal mat_init
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

  # Add Rz to model
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))

  # Compute residuals for effect 1
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

  # Set up for KL computation
  l <- 1
  model$lbf <- rep(0, params$L)
  model$alpha[l, ] <- rep(1/dat$p, dat$p)
  model$mu[l, ] <- rnorm(dat$p, sd = 0.1)
  model$mu2[l, ] <- model$mu[l, ]^2 + 0.1

  model <- compute_residuals.rss_lambda(data, params, model, l)
  model <- compute_kl.rss_lambda(data, params, model, l)

  expect_type(model$KL[l], "double")
  expect_length(model$KL[l], 1)
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

  # Variance should be positive
  post_var <- model$mu2[l, ] - model$mu[l, ]^2
  expect_true(all(post_var > -1e-10))

  # post_mean2 = post_var + post_mean^2
  expect_equal(model$mu2[l, ], post_var + model$mu[l, ]^2)
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
  # Precompute cached terms needed for Eloglik
  model <- precompute_rss_lambda_terms(data, model)

  e_loglik <- Eloglik.rss_lambda(data, model)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
  expect_true(is.finite(e_loglik))
})


test_that("loglik.rss_lambda computes log Bayes factors", {
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
  expect_true(abs(sum(model$alpha[l, ]) - 1) < 1e-10)
  expect_true(is.numeric(model$lbf[l]))
})


test_that("neg_loglik.rss_lambda returns negative log-likelihood", {
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

  V_param <- log(1.0)  # Log scale
  ser_stats <- compute_ser_statistics.rss_lambda(data, params, model, l = 1)
  neg_ll <- neg_loglik.rss_lambda(data, params, model, V_param, ser_stats)

  expect_type(neg_ll, "double")
  expect_length(neg_ll, 1)
})

test_that("get_ER2.rss_lambda computes expected squared residuals", {
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
  # Precompute cached terms needed for get_ER2
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

  # Update fitted values for effect 1
  old_Rz <- model$Rz
  model <- update_fitted_values.rss_lambda(data, params, model, l = 1)

  expect_true("Rz" %in% names(model))
  expect_length(model$Rz, dat$p)
  expect_type(model$Rz, "double")
})


test_that("update_variance_components.rss_lambda estimates sigma2", {
  dat <- setup_rss_lambda_data(seed = 29)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data

  params <- list(L = 5, scaled_prior_variance = 0.2, residual_variance = 1.0, estimate_residual_variance = TRUE)

  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  # Precompute cached terms
  model <- precompute_rss_lambda_terms(data, model)


  variance_update <- update_variance_components.rss_lambda(data, params, model)

  expect_type(variance_update, "list")
  expect_true("sigma2" %in% names(variance_update))
  expect_type(variance_update$sigma2, "double")
  expect_length(variance_update$sigma2, 1)
  expect_true(variance_update$sigma2 > 0)
  expect_true(variance_update$sigma2 <= 1 - dat$lambda)  # Upper bound
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

  # Change sigma2
  model$sigma2 <- 0.5

  # Update derived quantities
  updated_model <- update_derived_quantities.rss_lambda(data, params, model)

  expect_true("SinvRj" %in% names(updated_model))
  expect_true("RjSinvRj" %in% names(updated_model))
  expect_equal(dim(updated_model$SinvRj), c(dat$p, dat$p))
  expect_length(updated_model$RjSinvRj, dat$p)
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

  data <- result$data

  params <- list()
  scale_factors <- get_scale_factors.rss_lambda(data, params)

  expect_type(scale_factors, "double")
  expect_length(scale_factors, dat$p)
  expect_equal(scale_factors, rep(1, dat$p))
})


test_that("get_intercept.rss_lambda returns intercept_value", {
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


test_that("get_fitted.rss_lambda delegates to default method", {
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

  fitted <- get_fitted.rss_lambda(data, params, model)

  # Default method returns NULL for RSS data
  expect_null(fitted)
})


test_that("get_cs.rss_lambda returns NULL when coverage is NULL", {
  dat <- setup_rss_lambda_data(seed = 31)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data

  params <- list(L = 5, coverage = NULL, min_abs_corr = 0.5,
                 scaled_prior_variance = 0.2, residual_variance = 1.0,
                 prior_weights = rep(1/dat$p, dat$p))
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)

  cs <- get_cs.rss_lambda(data, params, model)

  expect_null(cs)
})

test_that("get_cs.rss_lambda returns NULL when min_abs_corr is NULL", {
  dat <- setup_rss_lambda_data(seed = 32)

  result <- rss_lambda_constructor(
    z = dat$z,
    R = dat$R,
    lambda = dat$lambda,
    n = dat$n
  )

  data <- result$data

  params <- list(L = 5, coverage = 0.95, min_abs_corr = NULL,
                 scaled_prior_variance = 0.2, residual_variance = 1.0,
                 prior_weights = rep(1/dat$p, dat$p))
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)

  cs <- get_cs.rss_lambda(data, params, model)

  expect_null(cs)
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

  # Add strong signal to create credible set
  model$alpha[1, 1] <- 0.95
  model$alpha[1, -1] <- 0.05 / (dat$p - 1)

  cs <- get_cs.rss_lambda(data, params, model)

  # May or may not find CS, but should not error
  expect_true(is.null(cs) || is.list(cs))
})


test_that("get_variable_names.rss_lambda assigns variable names to model", {
  dat <- setup_rss_lambda_data(seed = 34)

  # Create named z-scores
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


test_that("get_zscore.rss_lambda delegates to default method", {
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

  z <- get_zscore.rss_lambda(data, params, model)

  # Default returns NULL
  expect_null(z)
})

test_that("cleanup_model.rss_lambda removes temporary fields", {
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

  # Cleanup model
  cleaned <- cleanup_model.rss_lambda(data, params, model)

  # Check that temporary fields are removed
  expect_false("SinvRj" %in% names(cleaned))
  expect_false("RjSinvRj" %in% names(cleaned))
  expect_false("Rz" %in% names(cleaned))
  expect_false("Z" %in% names(cleaned))
  expect_false("zbar" %in% names(cleaned))
  expect_false("diag_postb2" %in% names(cleaned))

  # Check that essential fields remain
  expect_true("alpha" %in% names(cleaned))
  expect_true("mu" %in% names(cleaned))
  expect_true("mu2" %in% names(cleaned))
})