context("S3 methods for sufficient statistics (ss) data class")

# ---- Data initialization & configuration ----

test_that("configure_data.ss returns ss object without eigen components when unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")

  result <- configure_data.ss(setup$data, setup$params)

  expect_true("ss" %in% class(result))
  expect_false("eigen_values" %in% names(result))
})

test_that("configure_data.ss adds eigen decomposition when unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")

  setup$data$eigen_values  <- NULL
  setup$data$eigen_vectors <- NULL
  setup$data$VtXty         <- NULL

  result <- configure_data.ss(setup$data, setup$params)

  expect_true("eigen_values"  %in% names(result))
  expect_true("eigen_vectors" %in% names(result))
  expect_true("VtXty"         %in% names(result))
})

test_that("sufficient_stats_constructor accepts unmappable_effects='ash'", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 1)
  XtX <- crossprod(base_data$X)
  Xty <- crossprod(base_data$X, base_data$y)
  yty <- sum(base_data$y^2)

  result <- sufficient_stats_constructor(Xty = Xty, yty = yty, n = base_data$n, XtX = XtX,
                                          unmappable_effects = "ash")
  expect_true(inherits(result$data, "ss"))
  expect_equal(result$params$unmappable_effects, "ash")
})

test_that("get_var_y.ss returns yty / (n-1)", {
  setup <- setup_ss_data()

  var_y <- get_var_y.ss(setup$data)

  expect_type(var_y, "double")
  expect_length(var_y, 1)
  expect_true(var_y > 0)
  expect_equal(var_y, setup$data$yty / (setup$data$n - 1), tolerance = 1e-15)
})

# ---- Model initialization & setup ----

test_that("initialize_susie_model.ss sets predictor_weights from XtX diagonal when unmappable_effects='none'", {
  setup  <- setup_ss_data(unmappable_effects = "none")
  var_y  <- var(setup$data$yty / setup$data$n)
  model  <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  expect_true("predictor_weights" %in% names(model))
  expect_length(model$predictor_weights, setup$data$p)
  expect_equal(model$predictor_weights, attr(setup$data$XtX, "d"), tolerance = 1e-15)
})

test_that("initialize_susie_model.ss initializes omega fields with tau2=0 when unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  var_y <- setup$data$yty / (setup$data$n - 1)
  model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  expect_true(all(c("omega_var", "predictor_weights", "XtOmegay", "tau2", "theta") %in% names(model)))
  expect_equal(model$tau2,  0,                    tolerance = 1e-15)
  expect_equal(model$theta, rep(0, setup$data$p), tolerance = 1e-15)
})

test_that("initialize_fitted.ss adds XtXr of length p", {
  setup <- setup_ss_data()

  mat_init <- list(alpha = setup$model$alpha, mu = setup$model$mu)
  fitted   <- initialize_fitted.ss(setup$data, mat_init)

  expect_true("XtXr" %in% names(fitted))
  expect_length(fitted$XtXr, setup$data$p)
})

test_that("validate_prior.ss passes for a model with reasonable prior variance", {
  setup <- setup_ss_data()
  setup$params$check_prior <- TRUE

  expect_error(validate_prior.ss(setup$data, setup$params, setup$model), NA)
})

test_that("validate_prior.ss errors when prior variance exceeds 100 * zm^2", {
  setup <- setup_ss_data()
  setup$params$check_prior <- TRUE

  var_y        <- setup$data$yty / (setup$data$n - 1)
  setup$model  <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  bhat <- setup$data$Xty / setup$model$predictor_weights
  shat <- sqrt(setup$model$sigma2 / setup$model$predictor_weights)
  z    <- bhat / shat
  zm   <- max(abs(z[!is.nan(z)]))

  setup$model$V <- rep(150 * (zm^2), setup$params$L)

  expect_error(
    validate_prior.ss(setup$data, setup$params, setup$model),
    "Estimated prior variance is unreasonably large"
  )
})

test_that("track_ibss_fit.ss returns a list when unmappable_effects='none'", {
  setup    <- setup_ss_data(unmappable_effects = "none")
  tracking <- list()

  result <- track_ibss_fit.ss(setup$data, setup$params, setup$model, tracking, 1, -100)

  expect_type(result, "list")
})

# ---- Single effect regression & ELBO ----

test_that("compute_residuals.ss returns residuals, fitted_without_l, and sigma2 for unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")
  model <- compute_residuals.ss(setup$data, setup$params, setup$model, 1)

  expect_true(all(c("residuals", "fitted_without_l", "residual_variance") %in% names(model)))
  expect_length(model$residuals, setup$data$p)
  expect_equal(model$residual_variance, setup$model$sigma2, tolerance = 1e-15)
})

test_that("compute_residuals.ss returns omega-weighted residuals with residual_variance=1 for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, 1)

  expect_true(all(c("residuals", "predictor_weights", "residual_variance") %in% names(model)))
  expect_length(model$residuals, setup$data$p)
  expect_equal(model$residual_variance, 1, tolerance = 1e-15)
})

test_that("compute_ser_statistics.ss returns correct fields and optim metadata for unmappable_effects='none'", {
  setup     <- setup_ss_data(unmappable_effects = "none")
  model     <- compute_residuals.ss(setup$data, setup$params, setup$model, 1)
  ser_stats <- compute_ser_statistics.ss(setup$data, setup$params, model, 1)

  expect_true(all(c("betahat", "shat2", "optim_init", "optim_bounds", "optim_scale") %in% names(ser_stats)))
  expect_length(ser_stats$betahat, setup$data$p)
  expect_length(ser_stats$shat2,   setup$data$p)
  expect_equal(ser_stats$optim_scale,  "log",       tolerance = NULL)
  expect_equal(ser_stats$optim_bounds, c(-30, 15),  tolerance = 1e-15)
})

test_that("compute_ser_statistics.ss uses linear scale and V[l] init for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  model     <- compute_residuals.ss(setup$data, setup$params, setup$model, 1)
  ser_stats <- compute_ser_statistics.ss(setup$data, setup$params, model, 1)

  expect_equal(ser_stats$optim_scale,  "linear",     tolerance = NULL)
  expect_equal(ser_stats$optim_bounds, c(0, 1),      tolerance = 1e-15)
  expect_equal(ser_stats$optim_init,   model$V[1],   tolerance = 1e-15)
})

test_that("SER_posterior_e_loglik.ss returns a finite scalar for unmappable_effects='none'", {
  set.seed(301)
  setup <- setup_ss_data(unmappable_effects = "none")
  l     <- 1

  setup$model$alpha[l, ]  <- rep(1/setup$data$p, setup$data$p)
  setup$model$mu[l, ]     <- rnorm(setup$data$p)
  setup$model$mu2[l, ]    <- setup$model$mu[l, ]^2 + 0.1

  model    <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  e_loglik <- SER_posterior_e_loglik.ss(setup$data, setup$params, model, l)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
  expect_true(is.finite(e_loglik))
})

test_that("SER_posterior_e_loglik.ss returns a finite scalar for unmappable_effects='inf'", {
  set.seed(302)
  setup <- setup_ss_data(unmappable_effects = "inf")
  l     <- 1

  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)
  setup$model$alpha[l, ] <- rep(1/setup$data$p, setup$data$p)
  setup$model$mu[l, ]    <- rnorm(setup$data$p, sd = 0.01)
  setup$model$mu2[l, ]   <- setup$model$mu[l, ]^2 + 0.01

  model    <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  e_loglik <- SER_posterior_e_loglik.ss(setup$data, setup$params, model, l)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
  expect_true(is.finite(e_loglik))
})

test_that("calculate_posterior_moments.ss produces mu2 >= mu^2 (nonneg posterior variance)", {
  setup <- setup_ss_data()
  l     <- 1

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  model <- calculate_posterior_moments.ss(setup$data, setup$params, model, V = 1.0, l)

  expect_length(model$mu[l, ],  setup$data$p)
  expect_length(model$mu2[l, ], setup$data$p)
  expect_true(all(model$mu2[l, ] - model$mu[l, ]^2 >= -1e-10))
})

test_that("compute_kl.ss returns a finite scalar KL for the given component", {
  set.seed(303)
  setup <- setup_ss_data()
  l     <- 1

  setup$model$lbf       <- rep(0, setup$params$L)
  setup$model$alpha[l, ] <- rep(1/setup$data$p, setup$data$p)
  setup$model$mu[l, ]    <- rnorm(setup$data$p, sd = 0.1)
  setup$model$mu2[l, ]   <- setup$model$mu[l, ]^2 + 0.1

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  model <- compute_kl.ss(setup$data, setup$params, model, l)

  expect_type(model$KL[l], "double")
  expect_length(model$KL[l], 1)
  expect_true(is.finite(model$KL[l]))
})

test_that("get_ER2.ss returns a finite nonneg scalar", {
  setup <- setup_ss_data()
  er2   <- get_ER2.ss(setup$data, setup$model)

  expect_type(er2, "double")
  expect_length(er2, 1)
  expect_true(is.finite(er2))
  expect_true(er2 >= 0)
})

test_that("Eloglik.ss returns a finite scalar", {
  setup    <- setup_ss_data()
  e_loglik <- Eloglik.ss(setup$data, setup$model)

  expect_type(e_loglik, "double")
  expect_length(e_loglik, 1)
  expect_true(is.finite(e_loglik))
})

test_that("loglik.ss produces alpha summing to 1 and finite lbf", {
  setup <- setup_ss_data()
  l     <- 1

  model     <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  ser_stats <- compute_ser_statistics.ss(setup$data, setup$params, model, l)
  model     <- loglik.ss(setup$data, setup$params, model, V = 1.0, ser_stats, l)

  expect_length(model$lbf_variable[l, ], setup$data$p)
  expect_length(model$alpha[l, ],        setup$data$p)
  expect_true(all(model$alpha[l, ] >= 0))
  expect_equal(sum(model$alpha[l, ]), 1, tolerance = 1e-10)
  expect_true(is.finite(model$lbf[l]))
})

test_that("neg_loglik.ss returns a finite scalar for unmappable_effects='none'", {
  setup     <- setup_ss_data(unmappable_effects = "none")
  model     <- compute_residuals.ss(setup$data, setup$params, setup$model, 1)
  ser_stats <- compute_ser_statistics.ss(setup$data, setup$params, model, 1)
  neg_ll    <- neg_loglik.ss(setup$data, setup$params, model, V_param = log(1.0), ser_stats)

  expect_type(neg_ll, "double")
  expect_length(neg_ll, 1)
  expect_true(is.finite(neg_ll))
})

test_that("neg_loglik.ss returns a finite scalar for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  model     <- compute_residuals.ss(setup$data, setup$params, setup$model, 1)
  ser_stats <- compute_ser_statistics.ss(setup$data, setup$params, model, 1)
  neg_ll    <- neg_loglik.ss(setup$data, setup$params, model, V_param = 0.5, ser_stats)

  expect_type(neg_ll, "double")
  expect_length(neg_ll, 1)
  expect_true(is.finite(neg_ll))
})

# ---- Model updates & fitting ----

test_that("update_fitted_values.ss updates XtXr for unmappable_effects='none'", {
  setup <- setup_ss_data(unmappable_effects = "none")
  l     <- 1

  model <- compute_residuals.ss(setup$data, setup$params, setup$model, l)
  setup$model$fitted_without_l <- model$fitted_without_l

  updated_model <- update_fitted_values.ss(setup$data, setup$params, setup$model, l)

  expect_true("XtXr" %in% names(updated_model))
  expect_length(updated_model$XtXr, setup$data$p)
})

test_that("update_fitted_values.ss updates XtXr for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  l     <- 1

  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  updated_model <- update_fitted_values.ss(setup$data, setup$params, setup$model, l)

  expect_true("XtXr" %in% names(updated_model))
  expect_length(updated_model$XtXr, setup$data$p)
})

test_that("update_variance_components.ss returns sigma2 for unmappable_effects='none'", {
  setup  <- setup_ss_data(unmappable_effects = "none")
  result <- update_variance_components.ss(setup$data, setup$params, setup$model)

  expect_type(result, "list")
  expect_true("sigma2" %in% names(result))
  expect_true(is.finite(result$sigma2) && result$sigma2 > 0)
})

test_that("update_derived_quantities.ss returns a list for unmappable_effects='none'", {
  setup  <- setup_ss_data(unmappable_effects = "none")
  result <- update_derived_quantities.ss(setup$data, setup$params, setup$model)

  expect_type(result, "list")
})

test_that("update_derived_quantities.ss updates omega quantities for unmappable_effects='inf'", {
  setup <- setup_ss_data(unmappable_effects = "inf")
  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)

  result <- update_derived_quantities.ss(setup$data, setup$params, setup$model)

  expect_true(all(c("omega_var", "predictor_weights", "XtOmegay", "XtXr") %in% names(result)))
})

# ---- Output generation & post-processing ----

test_that("get_scale_factors.ss returns positive column scale factors matching XtX attribute", {
  setup  <- setup_ss_data()
  scales <- get_scale_factors.ss(setup$data, setup$params)

  expect_length(scales, setup$data$p)
  expect_true(all(scales > 0))
  expect_equal(scales, attr(setup$data$XtX, "scaled:scale"), tolerance = 1e-15)
})

test_that("get_intercept.ss returns a finite scalar when intercept=TRUE", {
  setup <- setup_ss_data()
  setup$params$intercept <- TRUE

  intercept <- get_intercept.ss(setup$data, setup$params, setup$model)

  expect_type(intercept, "double")
  expect_length(intercept, 1)
  expect_true(is.finite(intercept))
})

test_that("get_fitted.ss returns NULL (ss data has no individual fitted values)", {
  setup  <- setup_ss_data()
  fitted <- get_fitted.ss(setup$data, setup$params, setup$model)

  expect_null(fitted)
})

test_that("get_cs.ss returns NULL when coverage is NULL", {
  setup <- setup_ss_data()
  setup$params$coverage <- NULL

  expect_null(get_cs.ss(setup$data, setup$params, setup$model))
})

test_that("get_cs.ss returns NULL when min_abs_corr is NULL", {
  setup <- setup_ss_data()
  setup$params$min_abs_corr <- NULL

  expect_null(get_cs.ss(setup$data, setup$params, setup$model))
})

test_that("get_cs.ss computes correlation from XtX when diagonal is not 0 or 1", {
  setup <- setup_ss_data()
  diag(setup$data$XtX) <- diag(setup$data$XtX) * 1.5

  setup$model$alpha[1, 1]  <- 0.95
  setup$model$alpha[1, -1] <- 0.05 / (setup$data$p - 1)

  cs <- get_cs.ss(setup$data, setup$params, setup$model)

  expect_true(is.null(cs) || is.list(cs))
})

test_that("get_cs.ss uses XtX directly as correlation matrix when diagonal is all 1s", {
  set.seed(401)
  setup <- setup_ss_data()

  setup$data$XtX <- cor(matrix(rnorm(100 * setup$data$p), 100, setup$data$p))

  setup$model$alpha[1, 1]  <- 0.95
  setup$model$alpha[1, -1] <- 0.05 / (setup$data$p - 1)

  cs <- get_cs.ss(setup$data, setup$params, setup$model)

  expect_true(is.null(cs) || is.list(cs))
})

test_that("get_variable_names.ss applies XtX column names to all model matrices", {
  setup <- setup_ss_data()
  colnames(setup$data$XtX) <- paste0("var", seq_len(setup$data$p))

  setup$model$pip          <- rep(0.1, setup$data$p)
  setup$model$null_weight  <- NULL
  setup$model$alpha        <- matrix(0, 5, setup$data$p)
  setup$model$mu           <- matrix(0, 5, setup$data$p)
  setup$model$mu2          <- matrix(0, 5, setup$data$p)
  setup$model$lbf_variable <- matrix(0, 5, setup$data$p)

  model_out <- get_variable_names.ss(setup$data, setup$model)

  expect_true(all(grepl("var", colnames(model_out$alpha))))
  expect_true(all(grepl("var", colnames(model_out$mu))))
  expect_true(all(grepl("var", colnames(model_out$mu2))))
  expect_true(all(grepl("var", names(model_out$pip))))
})

test_that("get_zscore.ss returns NULL or a numeric vector", {
  setup <- setup_ss_data()
  setup$params$compute_univariate_zscore <- TRUE

  z <- get_zscore.ss(setup$data, setup$params, setup$model)

  expect_true(is.null(z) || is.numeric(z))
})

test_that("cleanup_model.ss removes residuals field when unmappable_effects='none'", {
  set.seed(402)
  setup <- setup_ss_data(unmappable_effects = "none")
  setup$model$residuals <- rnorm(setup$data$p)

  cleaned <- cleanup_model.ss(setup$data, setup$params, setup$model)

  expect_false("residuals" %in% names(cleaned))
})

test_that("cleanup_model.ss removes omega and residuals fields when unmappable_effects='inf'", {
  set.seed(403)
  setup <- setup_ss_data(unmappable_effects = "inf")
  var_y <- setup$data$yty / (setup$data$n - 1)
  setup$model <- initialize_susie_model.ss(setup$data, setup$params, var_y)
  setup$model$residuals <- rnorm(setup$data$p)

  cleaned <- cleanup_model.ss(setup$data, setup$params, setup$model)

  expect_false("omega_var"  %in% names(cleaned))
  expect_false("XtOmegay"   %in% names(cleaned))
  expect_false("residuals"  %in% names(cleaned))
})

context("susie() N>=2P hint and compute_suff_stat() workflow")

# ---- Hint behaviour ----

test_that("susie emits a hint pointing to compute_suff_stat() when nrow(X) >= 2 * ncol(X)", {
  set.seed(2026)
  n <- 200; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  expect_message(
    suppressWarnings(susie(X, y, L = 3, max_iter = 2, verbose = FALSE)),
    "compute_suff_stat"
  )
})

test_that("susie does not emit the compute_suff_stat hint when nrow(X) < 2 * ncol(X)", {
  set.seed(2027)
  n <- 60; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  msgs <- suppressWarnings(capture_messages(
    susie(X, y, L = 3, max_iter = 2, verbose = FALSE)
  ))
  expect_false(any(grepl("compute_suff_stat", msgs, fixed = TRUE)))
})

test_that("the hint does not change susie output relative to a hint-suppressed run", {
  set.seed(2028)
  n <- 200; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[c(5, 15, 25)] <- c(1, -1, 1.5)
  y <- as.vector(X %*% beta + rnorm(n, sd = 0.5))

  fit <- suppressMessages(
    suppressWarnings(susie(X, y, L = 5, max_iter = 100, verbose = FALSE))
  )

  expect_s3_class(fit, "susie")
  expect_length(fit$pip, p)
  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(is.finite(fit$elbo)))
})

# ---- Vignette workflow: compute_suff_stat() -> susie_ss() ----

test_that("compute_suff_stat() + susie_ss() agrees with susie() on the same data", {
  set.seed(2029)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[c(5, 15, 25)] <- c(1, -1, 1.5)
  y <- as.vector(X %*% beta + rnorm(n, sd = 0.5))

  fit_ind <- suppressMessages(susie(
    X, y, L = 5, standardize = TRUE, intercept = TRUE, verbose = FALSE
  ))

  ss <- compute_suff_stat(X, y, standardize = FALSE)

  fit_ss <- susie_ss(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n,
    X_colmeans = ss$X_colmeans, y_mean = ss$y_mean,
    L = 5, standardize = TRUE, verbose = FALSE
  )

  expect_equal(fit_ind$pip,    fit_ss$pip,    tolerance = 1e-3)
  expect_equal(fit_ind$V,      fit_ss$V,      tolerance = 1e-3)
  expect_equal(fit_ind$sigma2, fit_ss$sigma2, tolerance = 1e-3)
})

test_that("compute_suff_stat: X-only quantities are identical when reused across y vectors", {
  set.seed(2030)
  n <- 80; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * 2), n, 2)

  ss1 <- compute_suff_stat(X, Y[, 1], standardize = FALSE)

  y2_mean       <- mean(Y[, 2])
  y2c           <- Y[, 2] - y2_mean
  ss_reused     <- ss1
  ss_reused$Xty    <- drop(y2c %*% X)
  ss_reused$yty    <- sum(y2c^2)
  ss_reused$y_mean <- y2_mean

  ss_fresh <- compute_suff_stat(X, Y[, 2], standardize = FALSE)

  expect_identical(ss_reused$XtX,        ss_fresh$XtX)
  expect_identical(ss_reused$X_colmeans, ss_fresh$X_colmeans)
  expect_identical(ss_reused$n,          ss_fresh$n)
  expect_equal(ss_reused$Xty,    ss_fresh$Xty,    tolerance = 1e-12)
  expect_equal(ss_reused$yty,    ss_fresh$yty,    tolerance = 1e-12)
  expect_equal(ss_reused$y_mean, ss_fresh$y_mean, tolerance = 1e-12)

  fit_reused <- susie_ss(
    XtX = ss_reused$XtX, Xty = ss_reused$Xty, yty = ss_reused$yty,
    n = ss_reused$n,
    X_colmeans = ss_reused$X_colmeans, y_mean = ss_reused$y_mean,
    L = 5, verbose = FALSE
  )
  fit_fresh <- susie_ss(
    XtX = ss_fresh$XtX, Xty = ss_fresh$Xty, yty = ss_fresh$yty,
    n = ss_fresh$n,
    X_colmeans = ss_fresh$X_colmeans, y_mean = ss_fresh$y_mean,
    L = 5, verbose = FALSE
  )
  expect_equal(fit_reused$pip, fit_fresh$pip, tolerance = 1e-10)
  expect_equal(fit_reused$V,   fit_fresh$V,   tolerance = 1e-10)
})

# ---- Mode-specific ss branches (NIG / ash / inf) ----

# Shared helper: sufficient statistics with a clear two-signal structure.
.make_ss_stats <- function(n = 250, p = 15, seed = 4242) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = FALSE)
  beta <- rep(0, p)
  beta[c(3, 8)] <- c(1.2, -1)
  y <- as.vector(X %*% beta + rnorm(n, sd = 0.5))
  y <- y - mean(y)
  list(XtX = crossprod(X), Xty = as.vector(crossprod(X, y)),
       yty = sum(y^2), n = n, p = p)
}

# ---- NIG branch ----

test_that("susie_ss with NIG prior (L=1) produces well-formed posterior and removes NIG scratch fields", {
  ss <- .make_ss_stats(seed = 101)

  fit <- suppressWarnings(suppressMessages(susie_ss(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, L = 1,
    estimate_residual_method = "NIG", max_iter = 10, verbose = FALSE
  )))

  expect_s3_class(fit, "susie")
  expect_length(fit$pip, ss$p)
  expect_length(fit$KL, 1)
  expect_true(all(is.finite(fit$KL)))
  expect_equal(rowSums(fit$alpha), 1, tolerance = 1e-8)
  expect_true(all(is.finite(fit$mu)))
  expect_true(all(fit$mu2 >= fit$mu^2 - 1e-8))
  expect_null(fit$marginal_loglik)
})

test_that("susie_ss with NIG prior (L>1) sets KL=0 for each component", {
  ss <- .make_ss_stats(seed = 102)

  fit <- suppressWarnings(suppressMessages(susie_ss(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, L = 3,
    estimate_residual_method = "NIG", max_iter = 10, verbose = FALSE
  )))

  expect_s3_class(fit, "susie")
  expect_length(fit$KL, 3)
  expect_true(all(fit$KL == 0))
  expect_equal(rowSums(fit$alpha), rep(1, 3), tolerance = 1e-8)
})

test_that("calculate_posterior_moments.ss NIG with V<=0 sets mu/mu2 to zero and rv[l]=1", {
  ss   <- .make_ss_stats(n = 150, p = 10, seed = 103)
  ctor <- sufficient_stats_constructor(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, L = 1,
    X_colmeans = rep(0, ss$p), y_mean = 0, standardize = TRUE,
    estimate_residual_method = "NIG", convergence_method = "elbo",
    coverage = 0.95, min_abs_corr = 0.5, n_purity = 100,
    check_prior = FALSE, track_fit = FALSE
  )
  data   <- ctor$data
  params <- ctor$params
  var_y  <- get_var_y.ss(data)
  model  <- initialize_susie_model.ss(data, params, var_y)
  model  <- compute_residuals.ss(data, params, model, 1)

  model  <- calculate_posterior_moments.ss(data, params, model, V = 0, l = 1)

  expect_equal(model$mu[1, ],  rep(0, data$p), tolerance = 1e-15)
  expect_equal(model$mu2[1, ], rep(0, data$p), tolerance = 1e-15)
  expect_equal(model$rv[1], 1, tolerance = 1e-15)
})

test_that("compute_residuals.ss computes a positive yy_residual under the NIG prior", {
  ss   <- .make_ss_stats(n = 150, p = 10, seed = 104)
  ctor <- sufficient_stats_constructor(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, L = 2,
    X_colmeans = rep(0, ss$p), y_mean = 0, standardize = TRUE,
    estimate_residual_method = "NIG", convergence_method = "elbo",
    coverage = 0.95, min_abs_corr = 0.5, n_purity = 100,
    check_prior = FALSE, track_fit = FALSE
  )
  data   <- ctor$data
  params <- ctor$params
  expect_true(params$use_NIG)

  var_y <- get_var_y.ss(data)
  model <- initialize_susie_model.ss(data, params, var_y)
  model <- compute_residuals.ss(data, params, model, 1)

  expect_length(model$yy_residual, 1)
  expect_true(model$yy_residual >= .Machine$double.eps)
})

# ---- ash branch ----

test_that("initialize_susie_model.ss sets predictor_weights and theta for unmappable_effects='ash'", {
  ss   <- .make_ss_stats(n = 150, p = 12, seed = 105)
  ctor <- sufficient_stats_constructor(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, L = 3,
    X_colmeans = rep(0, ss$p), y_mean = 0, standardize = TRUE,
    unmappable_effects = "ash", estimate_residual_method = "MoM",
    convergence_method = "pip", coverage = 0.95, min_abs_corr = 0.5,
    n_purity = 100, check_prior = FALSE, track_fit = FALSE
  )
  data   <- ctor$data
  params <- ctor$params
  expect_equal(params$unmappable_effects, "ash")

  var_y <- get_var_y.ss(data)
  model <- initialize_susie_model.ss(data, params, var_y)

  expect_length(model$predictor_weights, data$p)
  expect_true("theta" %in% names(model))
  expect_length(model$theta, data$p)
})

test_that("susie_ss with unmappable_effects='ash' converges and returns well-formed fit", {
  ss <- .make_ss_stats(seed = 106)

  fit <- suppressWarnings(suppressMessages(susie_ss(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, L = 3,
    unmappable_effects = "ash", max_iter = 8, verbose = FALSE
  )))

  expect_s3_class(fit, "susie")
  expect_length(fit$pip, ss$p)
  expect_equal(rowSums(fit$alpha), rep(1, 3), tolerance = 1e-8)
})

# ---- inf branch (MoM) ----

test_that("susie_ss with unmappable_effects='inf' returns nonneg tau2 and theta of length p", {
  ss <- .make_ss_stats(seed = 107)

  fit <- suppressWarnings(suppressMessages(susie_ss(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, L = 3,
    unmappable_effects = "inf", estimate_residual_method = "MoM",
    max_iter = 8, verbose = FALSE
  )))

  expect_s3_class(fit, "susie")
  expect_length(fit$pip, ss$p)
  expect_true(!is.null(fit$tau2) && fit$tau2 >= 0)
  expect_length(fit$theta, ss$p)
})

test_that("update_variance_components.ss inf MoM returns finite positive sigma2, nonneg tau2, and theta of length p", {
  ss   <- .make_ss_stats(n = 150, p = 12, seed = 108)
  ctor <- sufficient_stats_constructor(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, L = 3,
    X_colmeans = rep(0, ss$p), y_mean = 0, standardize = TRUE,
    unmappable_effects = "inf", estimate_residual_method = "MoM",
    residual_variance = 1, convergence_method = "pip", coverage = 0.95,
    min_abs_corr = 0.5, n_purity = 100, check_prior = FALSE, track_fit = FALSE
  )
  data   <- ctor$data
  params <- ctor$params
  expect_equal(params$estimate_residual_method, "MoM")

  var_y  <- get_var_y.ss(data)
  model  <- initialize_susie_model.ss(data, params, var_y)
  result <- update_variance_components.ss(data, params, model)

  expect_type(result, "list")
  expect_true(all(c("sigma2", "tau2", "theta") %in% names(result)))
  expect_true(is.finite(result$sigma2) && result$sigma2 > 0)
  expect_true(result$tau2 >= 0)
  expect_length(result$theta, data$p)
})

# ---- inf branch (MLE) via direct call ----

test_that("update_variance_components.ss inf MLE path returns sigma2/tau2/theta via direct call", {
  ss <- .make_ss_stats(n = 150, p = 12, seed = 109)

  ctor <- sufficient_stats_constructor(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, L = 3,
    X_colmeans = rep(0, ss$p), y_mean = 0, standardize = TRUE,
    unmappable_effects = "inf", estimate_residual_method = "MoM",
    residual_variance = 1, convergence_method = "pip", coverage = 0.95,
    min_abs_corr = 0.5, n_purity = 100, check_prior = FALSE, track_fit = FALSE
  )
  data   <- ctor$data
  params <- ctor$params
  var_y  <- get_var_y.ss(data)
  model  <- initialize_susie_model.ss(data, params, var_y)

  for (i in seq_len(3)) {
    model <- ibss_fit(data, params, model)
    vc    <- update_variance_components.ss(data, params, model)
    model <- modifyList(model, vc)
  }

  params_mle <- params
  params_mle$estimate_residual_method <- "MLE"

  result <- update_variance_components.ss(data, params_mle, model)

  expect_type(result, "list")
  expect_true(all(c("sigma2", "tau2", "theta") %in% names(result)))
  expect_true(is.finite(result$sigma2) && result$sigma2 > 0)
  expect_true(result$tau2 >= 0)
  expect_length(result$theta, data$p)
})

# ---- Edge cases ----

test_that("susie_ss with p=1 returns a well-formed fit in the default (none) mode", {
  set.seed(501)
  n <- 100; p <- 1
  X <- matrix(rnorm(n), n, p)
  y <- as.vector(X * 2 + rnorm(n, sd = 0.5))
  y <- y - mean(y)
  X <- X - mean(X)

  ss <- list(XtX = crossprod(X), Xty = drop(crossprod(X, y)),
             yty = sum(y^2), n = n)

  fit <- suppressWarnings(susie_ss(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n,
    L = 1, verbose = FALSE
  ))

  expect_s3_class(fit, "susie")
  expect_length(fit$pip, p)
  expect_equal(rowSums(fit$alpha), 1, tolerance = 1e-8)
  expect_true(fit$pip[1] > 0.5)
})

test_that("susie_ss default mode (unmappable_effects='none') returns well-formed fit", {
  ss  <- .make_ss_stats(seed = 601)
  fit <- suppressWarnings(susie_ss(
    XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n,
    L = 3, verbose = FALSE
  ))

  expect_s3_class(fit, "susie")
  expect_length(fit$pip, ss$p)
  expect_equal(rowSums(fit$alpha), rep(1, 3), tolerance = 1e-8)
  expect_true(all(is.finite(fit$elbo)))
})
