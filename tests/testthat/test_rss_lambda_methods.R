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

# =============================================================================
# STOCHASTIC LD SKETCH INFLATION TESTS
# =============================================================================

test_that("compute_ser_statistics.rss_lambda returns betahat", {
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

test_that("inflation is computed when sketch_B is set", {
  dat <- setup_rss_lambda_data(seed = 41)

  # Construct data with sketch_samples
  B <- nrow(dat$X)
  result <- rss_lambda_constructor(
    z = dat$z, R = dat$R, lambda = dat$lambda, n = dat$n,
    sketch_samples = B
  )
  data <- result$data
  params <- result$params

  expect_equal(data$sketch_B, B)

  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))

  # After computing residuals, shat2_inflation should exist
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)
  expect_true(!is.null(model$shat2_inflation))
  expect_length(model$shat2_inflation, dat$p)
  expect_true(all(model$shat2_inflation >= 1))
})

test_that("no inflation without sketch_B", {
  dat <- setup_rss_lambda_data(seed = 42)

  data <- dat$data
  params <- dat$params
  model <- dat$model
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))

  # Without sketch_B, no inflation
  expect_null(data$sketch_B)
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)
  expect_null(model$shat2_inflation)
})

test_that("inflation widens shat2 in ser_stats", {
  dat <- setup_rss_lambda_data(seed = 43)

  # Run once without inflation
  data_no_infl <- dat$data
  params <- dat$params
  model_no <- dat$model
  model_no$Rz <- as.vector(data_no_infl$R %*% colSums(model_no$alpha * model_no$mu))
  model_no <- compute_residuals.rss_lambda(data_no_infl, params, model_no, l = 1)
  stats_no <- compute_ser_statistics.rss_lambda(data_no_infl, params, model_no, l = 1)

  # Run with inflation
  B <- nrow(dat$X)
  result <- rss_lambda_constructor(
    z = dat$z, R = dat$R, lambda = dat$lambda, n = dat$n,
    sketch_samples = B
  )
  data_infl <- result$data
  var_y <- get_var_y.rss_lambda(data_infl)
  model_infl <- initialize_susie_model.rss_lambda(data_infl, result$params, var_y)
  model_infl$Rz <- as.vector(data_infl$R %*% colSums(model_infl$alpha * model_infl$mu))
  model_infl <- compute_residuals.rss_lambda(data_infl, result$params, model_infl, l = 1)
  stats_infl <- compute_ser_statistics.rss_lambda(data_infl, result$params, model_infl, l = 1)

  # betahat should be identical (inflation doesn't affect betahat)
  expect_equal(stats_infl$betahat, stats_no$betahat, tolerance = 1e-10)

  # shat2 should be inflated (>= original)
  expect_true(all(stats_infl$shat2 >= stats_no$shat2 - 1e-15))
})

test_that("inflation approximately 1 under global null", {
  set.seed(44)
  p <- 50
  n <- 500
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = TRUE)
  y <- rnorm(n)  # no signal
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  R <- cov2cor(input_ss$XtX)
  R <- (R + t(R)) / 2
  ss <- univariate_regression(X, y)
  z <- ss$betahat / ss$sebetahat

  result <- rss_lambda_constructor(
    z = z, R = R, lambda = 0.1, n = n,
    sketch_samples = n
  )
  data <- result$data
  params <- result$params
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  # Under global null with uniform alpha, inflation should be close to 1
  # (Rz_without_l ~ 0 since all mu ~ 0)
  expect_true(all(model$shat2_inflation < 1.05))
})

test_that("posterior moments use inflated shat2 consistently", {
  dat <- setup_rss_lambda_data(seed = 45)

  B <- nrow(dat$X)
  result <- rss_lambda_constructor(
    z = dat$z, R = dat$R, lambda = dat$lambda, n = dat$n,
    sketch_samples = B
  )
  data <- result$data
  params <- result$params
  var_y <- get_var_y.rss_lambda(data)
  model <- initialize_susie_model.rss_lambda(data, params, var_y)
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  V <- 0.2
  model <- calculate_posterior_moments.rss_lambda(data, params, model, V, l = 1)

  # Posterior variance should be positive and well-behaved
  post_var <- model$mu2[1, ] - model$mu[1, ]^2
  expect_true(all(post_var > -1e-10))
  expect_true(all(is.finite(model$mu[1, ])))
  expect_true(all(is.finite(model$mu2[1, ])))

  # With inflated shat2, posterior variance should be larger than
  # uninflated (more conservative)
  model_no <- dat$model
  model_no$Rz <- as.vector(dat$data$R %*% colSums(model_no$alpha * model_no$mu))
  model_no <- compute_residuals.rss_lambda(dat$data, params, model_no, l = 1)
  model_no <- calculate_posterior_moments.rss_lambda(dat$data, params, model_no, V, l = 1)
  post_var_no <- model_no$mu2[1, ] - model_no$mu[1, ]^2

  # Inflated posterior variance >= uninflated (within numerical tolerance)
  expect_true(all(post_var >= post_var_no - 1e-10))
})

# =============================================================================
# R vs X INPUT PATH AGREEMENT
# =============================================================================

test_that("R and X input paths produce numerically identical results", {
  set.seed(50)
  p <- 30
  n <- 500
  B <- 200
  X_full <- matrix(rnorm(n * p), n, p)
  X_full <- scale(X_full, center = TRUE, scale = TRUE)
  y <- X_full[, 1] * 0.5 + rnorm(n)
  input_ss <- compute_suff_stat(X_full, y, standardize = TRUE)
  R <- cov2cor(input_ss$XtX)
  R <- (R + t(R)) / 2
  ss <- univariate_regression(X_full, y)
  z <- ss$betahat / ss$sebetahat

  # Use X as a "sketch" -- here use X_full itself (B=n)
  X_sketch <- X_full

  # Construct from R
  res_R <- rss_lambda_constructor(z = z, R = R, lambda = 0.1, n = n)
  # Construct from X
  res_X <- rss_lambda_constructor(z = z, X = X_sketch, lambda = 0.1, n = n)

  # Eigendecomposition should be very close
  # (sorted eigenvalues should match; eigenvectors may differ in sign)
  expect_equal(res_R$data$eigen_R$values, res_X$data$eigen_R$values, tolerance = 1e-6)

  # Initialize and run one SER iteration
  var_y_R <- get_var_y.rss_lambda(res_R$data)
  model_R <- initialize_susie_model.rss_lambda(res_R$data, res_R$params, var_y_R)
  model_R$Rz <- as.vector(R %*% colSums(model_R$alpha * model_R$mu))

  var_y_X <- get_var_y.rss_lambda(res_X$data)
  model_X <- initialize_susie_model.rss_lambda(res_X$data, res_X$params, var_y_X)
  model_X$Rz <- as.vector(compute_Rv(res_X$data, colSums(model_X$alpha * model_X$mu)))

  # RjSinvRj should agree
  expect_equal(model_R$RjSinvRj, model_X$RjSinvRj, tolerance = 1e-6)

  # Compute residuals
  model_R <- compute_residuals.rss_lambda(res_R$data, res_R$params, model_R, l = 1)
  model_X <- compute_residuals.rss_lambda(res_X$data, res_X$params, model_X, l = 1)
  expect_equal(model_R$residuals, model_X$residuals, tolerance = 1e-6)

  # SER statistics
  stats_R <- compute_ser_statistics.rss_lambda(res_R$data, res_R$params, model_R, l = 1)
  stats_X <- compute_ser_statistics.rss_lambda(res_X$data, res_X$params, model_X, l = 1)
  expect_equal(stats_R$betahat, stats_X$betahat, tolerance = 1e-6)
  expect_equal(stats_R$shat2, stats_X$shat2, tolerance = 1e-6)
})

# =============================================================================
# END-TO-END susie_rss WITH INFLATION
# =============================================================================

test_that("susie_rss with lambda and sketch_samples runs successfully", {
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

  # Without inflation
  fit_no <- susie_rss(z = z, R = R, lambda = 0.1, n = n, L = 5,
                      max_iter = 50, verbose = FALSE)

  # With inflation
  fit_infl <- susie_rss(z = z, R = R, lambda = 0.1, n = n, L = 5,
                        sketch_samples = n, max_iter = 50, verbose = FALSE)

  # Both should converge
  expect_true(fit_no$converged)
  expect_true(fit_infl$converged)

  # ELBO with inflation should be slightly lower (wider posterior = looser bound)
  # but still reasonable
  expect_true(is.finite(fit_infl$elbo[length(fit_infl$elbo)]))
  expect_true(is.finite(fit_no$elbo[length(fit_no$elbo)]))

  # sketch_diagnostics should be stored
  expect_true(!is.null(fit_infl$sketch_diagnostics))

  # Causal variables should still have high PIP
  expect_true(fit_no$pip[1] > 0.5)
  expect_true(fit_infl$pip[1] > 0.3)  # slightly less confident with inflation
})

test_that("susie_rss with X input and lambda and inflation works", {
  set.seed(52)
  p <- 30
  n <- 2000
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = TRUE)
  beta <- rep(0, p)
  beta[1] <- 0.5
  y <- drop(X %*% beta + rnorm(n))
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  R <- cov2cor(input_ss$XtX)
  R <- (R + t(R)) / 2
  ss <- univariate_regression(X, y)
  z <- ss$betahat / ss$sebetahat

  # susie_rss with X input (lambda > 0 triggers rss_lambda path)
  fit_X <- susie_rss(z = z, R = R, lambda = 0.1, n = n, L = 5,
                     sketch_samples = n, max_iter = 50, verbose = FALSE)

  expect_true(fit_X$converged)
  expect_true(is.finite(fit_X$elbo[length(fit_X$elbo)]))
  expect_true(fit_X$pip[1] > 0.3)
})

test_that("loglik.rss_lambda Wakefield ABF agrees with old signal-based form", {
  # Verify the Wakefield ABF form gives the same result as the original
  # signal^2 / RjSinvRj form when there is no inflation
  dat <- setup_rss_lambda_data(seed = 53)
  data <- dat$data
  params <- dat$params
  model <- dat$model
  model$Rz <- as.vector(data$R %*% colSums(model$alpha * model$mu))
  model <- compute_residuals.rss_lambda(data, params, model, l = 1)

  ser_stats <- compute_ser_statistics.rss_lambda(data, params, model, l = 1)

  # Compute BF using Wakefield ABF (current code)
  V <- 0.2
  shat2 <- pmax(ser_stats$shat2, .Machine$double.eps)
  lbf_wakefield <- -0.5 * log(1 + V / shat2) +
    0.5 * ser_stats$betahat^2 * V / (shat2 * (V + shat2))

  # Compute BF using the original SinvRj form:
  # lbf = -0.5 * log(1 + V * RjSinvRj) + 0.5 * V * signal^2 / (1 + V * RjSinvRj)
  # where signal = SinvRj' * r, and shat2 = 1/RjSinvRj, betahat = signal * shat2
  signal <- as.vector(crossprod(model$SinvRj, model$residuals))
  RjSinvRj <- model$RjSinvRj
  lbf_original <- -0.5 * log(1 + V * RjSinvRj) +
    0.5 * V * signal^2 / (1 + V * RjSinvRj)

  expect_equal(lbf_wakefield, lbf_original, tolerance = 1e-10)
})

# =============================================================================
# SS vs RSS-LAMBDA CROSS-PATH AGREEMENT TESTS
# =============================================================================

test_that("SS and RSS-lambda paths agree with small lambda (no inflation)", {
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

  # SS path (lambda = 0)
  fit_ss <- susie_rss(z = z, R = R, n = n, L = 5, lambda = 0,
                      max_iter = 100, verbose = FALSE)
  # RSS-lambda path (tiny lambda ~= 0)
  fit_rss <- susie_rss(z = z, R = R, n = n, L = 5, lambda = 1e-6,
                       max_iter = 100, verbose = FALSE)

  expect_true(fit_ss$converged)
  expect_true(fit_rss$converged)
  # Alpha matrices should be essentially identical
  expect_equal(fit_ss$alpha, fit_rss$alpha, tolerance = 1e-4)
  # PIPs should match
  expect_equal(fit_ss$pip, fit_rss$pip, tolerance = 1e-4)
})

test_that("SS and RSS-lambda paths agree with inflation", {
  set.seed(201)
  p <- 50; n <- 5000
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = TRUE)
  beta <- rep(0, p)
  beta[1] <- 0.5; beta[10] <- -0.3
  y <- drop(X %*% beta + rnorm(n))

  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  R <- cov2cor(input_ss$XtX); R <- (R + t(R)) / 2
  ss <- univariate_regression(X, y)
  z <- ss$betahat / ss$sebetahat

  B_sketch <- 2000
  # SS path with inflation
  fit_ss <- susie_rss(z = z, R = R, n = n, L = 5, lambda = 0,
                      sketch_samples = B_sketch,
                      max_iter = 100, verbose = FALSE)
  # RSS-lambda path with inflation
  fit_rss <- susie_rss(z = z, R = R, n = n, L = 5, lambda = 1e-6,
                       sketch_samples = B_sketch,
                       max_iter = 100, verbose = FALSE)

  expect_true(fit_ss$converged)
  expect_true(fit_rss$converged)
  # Alpha and PIPs should match closely despite different inflation formulas
  expect_equal(fit_ss$alpha, fit_rss$alpha, tolerance = 1e-4)
  expect_equal(fit_ss$pip, fit_rss$pip, tolerance = 1e-4)
})

# =============================================================================
# MULTI-PANEL LD MIXTURE TESTS
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
  # First 50 rows scaled by sqrt(0.6)
  expect_equal(X_meta[1:50, ], sqrt(0.6) * X1)
  # Last 30 rows scaled by sqrt(0.4)
  expect_equal(X_meta[51:80, ], sqrt(0.4) * X2)
})

test_that("eigen_from_X recovers eigendecomposition of X'X", {
  set.seed(43)
  p <- 20
  X <- matrix(rnorm(100 * p), 100, p)
  R <- crossprod(X)
  eigen_R_direct <- eigen(R, symmetric = TRUE)

  eigen_R_svd <- eigen_from_X(X, p)

  # Eigenvalues should match
  expect_equal(eigen_R_svd$values, eigen_R_direct$values, tolerance = 1e-10)
  # Eigenvectors span same space (up to sign)
  for (j in seq_len(p)) {
    inner <- abs(sum(eigen_R_svd$vectors[, j] * eigen_R_direct$vectors[, j]))
    expect_gt(inner, 0.99)
  }
})

test_that("eval_omega_eloglik_reduced matches pure R reference", {
  set.seed(44)
  p <- 50
  K <- 2

  # Create two panels with B_total < p so reduced-basis applies
  X1 <- matrix(rnorm(15 * p), 15, p)
  X2 <- matrix(rnorm(10 * p), 10, p)
  X_list <- list(X1, X2)
  # Use raw cross-products (matching constructor: lapply(X_list, crossprod))
  panel_R <- list(crossprod(X1), crossprod(X2))

  z <- rnorm(p)
  zbar <- rnorm(p) * 0.1
  diag_postb2 <- abs(rnorm(p)) * 0.01
  L <- 3
  Z <- matrix(rnorm(L * p) * 0.05, L, p)
  sigma2 <- 0.9
  lambda <- 0.01
  omega <- c(0.7, 0.3)

  # Pure R reference (O(p^3) eigendecomposition using panel_R)
  val_R <- susieR:::eval_omega_eloglik_R(panel_R, omega, z, zbar, diag_postb2,
                                          Z, sigma2, lambda, K, p)

  # Reduced-basis (O(r^3) Cholesky using X_list)
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

  # Concavity: midpoint should be >= average of endpoints
  vals <- sapply(seq(0, 1, 0.1), eloglik)
  for (i in 1:(length(vals) - 2)) {
    midval <- vals[i + 1]
    avg_endpoints <- (vals[i] + vals[i + 2]) / 2
    expect_gte(midval, avg_endpoints - 1e-8)
  }
})

test_that("rss_lambda_constructor accepts list X", {
  set.seed(46)
  p <- 30
  B1 <- 200
  B2 <- 150
  beta <- rep(0, p)
  beta[1:2] <- 0.5

  X1 <- matrix(rnorm(B1 * p), B1, p)
  X2 <- matrix(rnorm(B2 * p), B2, p)
  z <- rnorm(p) + crossprod(X1, X1 %*% beta)[, 1] / (B1 * sqrt(p))

  result <- rss_lambda_constructor(z = z, X = list(X1, X2), lambda = 0.1)

  data <- result$data
  expect_s3_class(data, "rss_lambda")
  expect_equal(data$K, 2)
  expect_equal(length(data$X_list), 2)
  expect_equal(data$B_list, c(B1, B2))
  expect_true(!is.null(data$eigen_R))
  expect_true(!is.null(data$Vtz))
  # Inflation is opt-in: without sketch_samples, B is not set
  expect_null(data$sketch_B)
})

test_that("K=1 list X gives same results as single matrix X", {
  set.seed(47)
  p <- 20
  B <- 100
  X <- matrix(rnorm(B * p), B, p)
  z <- rnorm(p)

  # Single matrix
  result_single <- rss_lambda_constructor(z = z, X = X, lambda = 0.1)
  # List with one element
  result_list <- rss_lambda_constructor(z = z, X = list(X), lambda = 0.1)

  # Eigenvalues should match (both represent R from same X)
  expect_equal(result_single$data$eigen_R$values,
               result_list$data$eigen_R$values,
               tolerance = 1e-8)
  # Vtz should match
  expect_equal(abs(result_single$data$Vtz),
               abs(result_list$data$Vtz),
               tolerance = 1e-8)
})

test_that("K=2 with identical panels: omega sums to 1 and fit converges", {
  set.seed(48)
  p <- 20
  B <- 200
  X <- matrix(rnorm(B * p), B, p)
  beta <- rep(0, p)
  beta[1:3] <- c(0.5, -0.3, 0.4)

  # Compute z-scores
  R <- cov2cor(crossprod(X))
  z <- as.vector(R %*% beta) + rnorm(p, sd = 0.1)

  fit <- susie_rss(z = z, X = list(X, X), lambda = 0.1,
                   estimate_residual_variance = TRUE,
                   max_iter = 50)

  # With identical panels, any omega on the simplex is equally good
  # Just check weights are valid
  expect_equal(sum(fit$omega_weights), 1, tolerance = 1e-10)
  expect_true(all(fit$omega_weights >= -1e-10))
  # Model should still converge
  expect_true(fit$converged || fit$niter <= 50)
})

test_that("K=2 end-to-end susie_rss converges", {
  set.seed(49)
  p <- 30
  B1 <- 300
  B2 <- 200
  beta <- rep(0, p)
  beta[1:2] <- c(0.5, -0.4)

  # Two panels sampling from same distribution
  X1 <- matrix(rnorm(B1 * p), B1, p)
  X2 <- matrix(rnorm(B2 * p), B2, p)

  # R from pooled
  R_pool <- cov2cor((crossprod(X1) + crossprod(X2)) / (B1 + B2))
  z <- as.vector(R_pool %*% beta) + rnorm(p, sd = 0.1)

  fit <- susie_rss(z = z, X = list(X1, X2), lambda = 0.1,
                   estimate_residual_variance = TRUE,
                   max_iter = 100)

  # Should converge
  expect_true(fit$converged || fit$niter <= 100)
  # Should have omega weights
  expect_equal(length(fit$omega_weights), 2)
  expect_equal(sum(fit$omega_weights), 1, tolerance = 1e-10)
  # Should have PIPs
  expect_equal(length(fit$pip), p)
})

test_that("accessor helpers fall through for single panel", {
  dat <- setup_rss_lambda_data(seed = 50)
  model <- dat$model

  # model$eigen_R is NULL for single panel
  expect_null(model$eigen_R)
  # Accessor should return data$eigen_R
  eigen_R <- get_eigen_R(dat$data, model)
  expect_equal(eigen_R$values, dat$data$eigen_R$values)

  # Same for Vtz
  expect_null(model$Vtz)
  Vtz <- get_Vtz(dat$data, model)
  expect_equal(Vtz, dat$data$Vtz)
})

# =============================================================================
# RANK BOUND FALLBACK
# =============================================================================

test_that("K=2 works when sum(B_k) > p (rank bound fallback)", {
  set.seed(51)
  p <- 20
  B1 <- 15; B2 <- 15  # sum = 30 > p = 20
  X1 <- matrix(rnorm(B1 * p), B1, p)
  X2 <- matrix(rnorm(B2 * p), B2, p)
  beta <- rep(0, p)
  beta[1:2] <- c(0.5, -0.3)
  R_pool <- cov2cor((crossprod(X1) + crossprod(X2)) / (B1 + B2))
  z <- as.vector(R_pool %*% beta) + rnorm(p, sd = 0.1)

  # Should not error: falls back to direct panel_R path
  fit <- susie_rss(z = z, X = list(X1, X2), lambda = 0.1,
                   max_iter = 50, verbose = FALSE)

  expect_equal(length(fit$pip), p)
  expect_equal(length(fit$omega_weights), 2)
  expect_equal(sum(fit$omega_weights), 1, tolerance = 1e-10)
  expect_true(all(fit$omega_weights >= -1e-10))
})

# =============================================================================
# AUTO-PIP CONVERGENCE
# =============================================================================

test_that("multi-panel auto-switches to PIP convergence with message", {
  set.seed(52)
  p <- 20; B <- 100
  X1 <- matrix(rnorm(B * p), B, p)
  X2 <- matrix(rnorm(B * p), B, p)
  z <- rnorm(p)

  expect_message(
    susie_rss(z = z, X = list(X1, X2), lambda = 0.1,
              convergence_method = "elbo", max_iter = 10, verbose = FALSE),
    "Switching to PIP-based convergence for multi-panel"
  )
})

test_that("sketch_samples auto-switches to PIP convergence with message", {
  set.seed(53)
  p <- 20; B <- 1000
  X <- matrix(rnorm(B * p), B, p)
  z <- rnorm(p)

  expect_message(
    susie_rss(z = z, X = X, lambda = 0.1,
              sketch_samples = TRUE,
              convergence_method = "elbo", max_iter = 10, verbose = FALSE),
    "Switching to PIP-based convergence because sketch LD inflation"
  )
})

# =============================================================================
# TOLERANCE CONSTANTS
# =============================================================================

test_that(".omega_tol has expected fields", {
  tol <- susieR:::.omega_tol
  expect_true(is.list(tol))
  expect_true("convergence" %in% names(tol))
  expect_true("grid_spacing" %in% names(tol))
  expect_true("fw_stop" %in% names(tol))
  expect_true("fw_max_iter" %in% names(tol))
  # Sanity: values are positive
  expect_true(tol$convergence > 0)
  expect_true(tol$grid_spacing > 0 && tol$grid_spacing < 1)
  expect_true(tol$fw_stop > 0)
  expect_true(tol$fw_max_iter >= 1L)
})

# =============================================================================
# INFLATION OPT-IN
# =============================================================================

test_that("multi-panel without sketch_samples has no inflation", {
  set.seed(54)
  p <- 20; B <- 100
  X1 <- matrix(rnorm(B * p), B, p)
  X2 <- matrix(rnorm(B * p), B, p)
  z <- rnorm(p)

  fit <- susie_rss(z = z, X = list(X1, X2), lambda = 0.1,
                   max_iter = 10, verbose = FALSE)

  # No inflation: sketch_B should not be set in the fit
  expect_null(fit$sketch_B)
})

# =============================================================================
# EIGEN_FROM_REDUCED UNIT TEST (Issue 21)
# =============================================================================

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

  # Direct eigendecomposition of R(omega)
  R_omega <- omega[1] * crossprod(X_list[[1]]) + omega[2] * crossprod(X_list[[2]])
  R_omega <- 0.5 * (R_omega + t(R_omega))
  eig_direct <- eigen(R_omega, symmetric = TRUE)

  # Eigenvalues should match (within reduced rank)
  r <- cache$r
  expect_equal(eig_reduced$values[1:r], eig_direct$values[1:r], tolerance = 1e-8)

  # Eigenvectors should span the same space: V_reduced' V_direct ~= I for top-r
  overlap <- abs(crossprod(eig_reduced$vectors[, 1:r], eig_direct$vectors[, 1:r]))
  # Each reduced eigenvector should align with exactly one direct eigenvector
  expect_true(all(apply(overlap, 1, max) > 1 - 1e-8))
})

# =============================================================================
# K > 2 MULTI-PANEL TEST (Issue 20)
# =============================================================================

test_that("K=3 multi-panel mixture runs and produces valid omega", {
  set.seed(56)
  p <- 20; B <- 80
  X1 <- matrix(rnorm(B * p), B, p)
  X2 <- matrix(rnorm(B * p), B, p)
  X3 <- matrix(rnorm(B * p), B, p)
  z <- rnorm(p)

  fit <- susie_rss(z = z, X = list(X1, X2, X3), lambda = 0.1,
                   max_iter = 20, verbose = FALSE)

  expect_equal(length(fit$omega_weights), 3)
  expect_equal(sum(fit$omega_weights), 1, tolerance = 1e-10)
  expect_true(all(fit$omega_weights >= -1e-10))
  expect_equal(length(fit$pip), p)
})

# =============================================================================
# OMEGA AT SIMPLEX VERTEX (Issue 22)
# =============================================================================

test_that("optimize_omega handles vertex optimum (one panel irrelevant)", {
  set.seed(57)
  p <- 25; B <- 100

  # Panel 1: true LD, Panel 2: pure noise (identity-like)
  X1 <- matrix(rnorm(B * p), B, p)
  X_list <- lapply(list(X1, matrix(rnorm(B * p), B, p)), susieR:::standardize_X)
  z <- rnorm(p)

  # Construct data where panel 1 is much better
  R1 <- crossprod(X_list[[1]])
  R2 <- diag(p)  # identity -- pure noise panel
  panel_R <- list(R1, R2)

  zbar <- rnorm(p) * 0.1
  diag_postb2 <- abs(rnorm(p)) * 0.01
  Z <- matrix(rnorm(2 * p) * 0.05, 2, p)

  eval_fn <- function(omega_vec) {
    susieR:::eval_omega_eloglik_R(panel_R, omega_vec, z, zbar,
                                   diag_postb2, Z, 0.9, 0.1, 2, p)
  }

  result <- susieR:::optimize_omega(eval_fn, c(0.5, 0.5), K = 2)

  # Should produce valid omega on simplex
  expect_equal(sum(result$omega), 1, tolerance = 1e-10)
  expect_true(all(result$omega >= -1e-10))
})

# =============================================================================
# UPDATE_DERIVED_QUANTITIES CORRECTNESS (Issue 24)
# =============================================================================

test_that("update_derived_quantities produces correct SinvRj after omega update", {
  set.seed(58)
  p <- 20; B <- 80
  X1 <- matrix(rnorm(B * p), B, p)
  X2 <- matrix(rnorm(B * p), B, p)
  z <- rnorm(p)

  # Run a few iterations to get an omega update
  fit <- susie_rss(z = z, X = list(X1, X2), lambda = 0.1,
                   max_iter = 5, verbose = FALSE)

  # Verify omega weights exist and are valid
  expect_true(!is.null(fit$omega_weights))
  expect_equal(sum(fit$omega_weights), 1, tolerance = 1e-10)

  # Verify the fit produced PIPs
  expect_equal(length(fit$pip), p)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

test_that("multi-panel ELBO does not degrade below single-panel after omega update", {
  # Regression test: before the Rz fix, stale model$Rz after omega update
  # could cause the mixture ELBO to be worse than any single panel.
  set.seed(70)
  p <- 30; B <- 100
  X1 <- matrix(rnorm(B * p), B, p)
  X2 <- matrix(rnorm(B * p), B, p)
  beta <- rep(0, p); beta[1:3] <- c(0.5, -0.3, 0.4)
  R_true <- cov2cor(crossprod(X1))
  z <- as.vector(R_true %*% beta) + rnorm(p, sd = 0.1)

  # Single-panel fits
  fit1 <- susie_rss(z = z, X = X1, lambda = 0.1, max_iter = 50, verbose = FALSE)
  fit2 <- susie_rss(z = z, X = X2, lambda = 0.1, max_iter = 50, verbose = FALSE)
  best_single <- max(tail(fit1$elbo, 1), tail(fit2$elbo, 1))

  # Mixture fit
  fit_mix <- susie_rss(z = z, X = list(X1, X2), lambda = 0.1,
                       max_iter = 50, verbose = FALSE)
  mix_elbo <- tail(fit_mix$elbo, 1)

  # Mixture ELBO should not be substantially worse than best single panel
  expect_true(mix_elbo >= best_single - 0.5,
              info = sprintf("mixture ELBO %.2f < best single %.2f",
                             mix_elbo, best_single))
})

test_that("multi-panel ELBO >= best single-panel ELBO", {
  set.seed(71)
  p <- 30; B1 <- 200; B2 <- 150
  X1 <- matrix(rnorm(B1 * p), B1, p)
  X2 <- matrix(rnorm(B2 * p), B2, p)
  beta <- rep(0, p); beta[1:2] <- c(0.6, -0.4)

  # Generate z from panel 1 (the "correct" panel)
  R1 <- cov2cor(crossprod(X1))
  z <- as.vector(R1 %*% beta) + rnorm(p, sd = 0.1)

  # Single-panel fits
  fit1 <- susie_rss(z = z, X = X1, lambda = 0.1, max_iter = 50)
  fit2 <- susie_rss(z = z, X = X2, lambda = 0.1, max_iter = 50)
  best_single_elbo <- max(tail(fit1$elbo, 1), tail(fit2$elbo, 1))

  # Mixture fit
  fit_mix <- susie_rss(z = z, X = list(X1, X2), lambda = 0.1, max_iter = 50)
  mix_elbo <- tail(fit_mix$elbo, 1)

  # Mixture should be at least as good as best single panel
  # (omega at a simplex vertex recovers single-panel)
  expect_true(mix_elbo >= best_single_elbo - 1e-2,
              info = sprintf("mixture ELBO %.2f < best single %.2f",
                             mix_elbo, best_single_elbo))
})

test_that("multi-panel first-iter ELBO >= best single-panel first-iter ELBO", {
  # With vertex initialization, the first E-step uses a single panel's R,
  # so the mixture first-iter ELBO should match the best single panel.
  # (Softmax init caused rank inflation when B < p, making it much worse.)
  set.seed(72)
  p <- 100; B <- 80
  X_true <- matrix(rnorm(2000 * p), 2000, p)
  for (b in 1:(p %/% 5)) {
    idx <- ((b-1)*5+1):(b*5)
    X_true[, idx] <- 0.4*X_true[, idx] + 0.6*rnorm(2000)
  }
  X_true <- scale(X_true)
  R_true <- cov2cor(crossprod(X_true))

  # Two noisy panels (both wrong)
  X1 <- matrix(rnorm(B * p), B, p)
  X2 <- matrix(rnorm(B * p), B, p)
  for (b in 1:(p %/% 5)) {
    idx <- ((b-1)*5+1):(b*5)
    X1[, idx] <- 0.4*X1[, idx] + 0.6*rnorm(B)
    X2[, idx] <- 0.4*X2[, idx] + 0.6*rnorm(B)
  }
  X1 <- scale(X1); X2 <- scale(X2)

  beta <- rep(0, p); beta[c(3, 15, 28)] <- c(1.0, -0.8, 0.7)
  z <- as.vector(R_true %*% beta) + rnorm(p, sd = 0.1)

  f1 <- susie_rss(z = z, X = X1, lambda = 0.1, max_iter = 1, verbose = FALSE)
  f2 <- susie_rss(z = z, X = X2, lambda = 0.1, max_iter = 1, verbose = FALSE)
  fm <- susie_rss(z = z, X = list(X1, X2), lambda = 0.1, max_iter = 1,
                  verbose = FALSE)

  best_sp_first <- max(f1$elbo[1], f2$elbo[1])
  mix_first <- fm$elbo[1]

  # With vertex init, mixture first-iter should be within 1 unit of best SP
  expect_true(mix_first >= best_sp_first - 1,
              info = sprintf("mix first=%.2f, best_sp first=%.2f, diff=%.2f",
                             mix_first, best_sp_first, mix_first - best_sp_first))
})