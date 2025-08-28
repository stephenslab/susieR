susie <- function(X, y, L = min(10, ncol(X)),
                  scaled_prior_variance = 0.2,
                  residual_variance = NULL,
                  prior_weights = NULL,
                  null_weight = 0,
                  standardize = TRUE,
                  intercept = TRUE,
                  estimate_residual_variance = TRUE,
                  estimate_residual_method = c("MLE", "MoM", "Servin_Stephens"),
                  estimate_prior_variance = TRUE,
                  estimate_prior_method = c("optim", "EM", "simple"),
                  unmappable_effects = c("none", "inf", "ash"),
                  check_null_threshold = 0,
                  prior_tol = 1e-9,
                  residual_variance_upperbound = Inf,
                  model_init = NULL,
                  coverage = 0.95,
                  min_abs_corr = 0.5,
                  compute_univariate_zscore = FALSE,
                  na.rm = FALSE,
                  max_iter = 100,
                  tol = 1e-3,
                  convergence_method = c("elbo", "pip"),
                  verbose = FALSE,
                  track_fit = FALSE,
                  residual_variance_lowerbound = var(drop(y)) / 1e4,
                  refine = FALSE,
                  n_purity = 100,
                  alpha0 = 0,
                  beta0 = 0) {
  
  # Validate method arguments
  unmappable_effects       <- match.arg(unmappable_effects)
  estimate_prior_method    <- match.arg(estimate_prior_method)
  estimate_residual_method <- match.arg(estimate_residual_method)
  convergence_method       <- match.arg(convergence_method)

  # Construct data object
  data <- individual_data_constructor(
    X, y, intercept, standardize, na.rm,
    prior_weights, null_weight, unmappable_effects,
    estimate_residual_method, convergence_method,
    estimate_prior_method, alpha0, beta0
  )

  model <- susie_workhorse(
    data, L, intercept, standardize, scaled_prior_variance,
    residual_variance, data$prior_weights, data$null_weight,
    model_init, estimate_prior_variance, data$estimate_prior_method,
    check_null_threshold, estimate_residual_variance,
    estimate_residual_method,
    residual_variance_lowerbound, residual_variance_upperbound,
    max_iter, tol, verbose, track_fit, coverage, min_abs_corr,
    prior_tol, n_purity, compute_univariate_zscore,
    check_prior = FALSE, data$convergence_method,
    alpha0 = alpha0, beta0 = beta0, refine = refine
  )

  return(model)
}

susie_ss <- function(XtX, Xty, yty, n,
                     L = min(10, ncol(XtX)),
                     X_colmeans = NA, y_mean = NA,
                     maf = NULL, maf_thresh = 0,
                     check_input = FALSE,
                     r_tol = 1e-8,
                     standardize = TRUE,
                     scaled_prior_variance = 0.2,
                     residual_variance = NULL,
                     prior_weights = NULL,
                     null_weight = 0,
                     model_init = NULL,
                     estimate_residual_variance = TRUE,
                     estimate_residual_method = c("MLE", "MoM", "Servin_Stephens"),
                     residual_variance_lowerbound = 0,
                     residual_variance_upperbound = Inf,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = c("optim", "EM", "simple"),
                     unmappable_effects = c("none", "inf"),
                     check_null_threshold = 0,
                     prior_tol = 1e-9,
                     max_iter = 100,
                     tol = 1e-3,
                     convergence_method = c("elbo", "pip"),
                     coverage = 0.95,
                     min_abs_corr = 0.5,
                     n_purity = 100,
                     verbose = FALSE,
                     track_fit = FALSE,
                     check_prior = FALSE,
                     ...) {
  
  # Validate method arguments
  unmappable_effects       <- match.arg(unmappable_effects)
  estimate_prior_method    <- match.arg(estimate_prior_method)
  estimate_residual_method <- match.arg(estimate_residual_method)
  convergence_method       <- match.arg(convergence_method)

  # Construct data object
  data <- sufficient_stats_constructor(
    XtX, Xty, yty, n, X_colmeans, y_mean,
    maf, maf_thresh, standardize, r_tol, check_input,
    prior_weights, null_weight, unmappable_effects, estimate_residual_method,
    convergence_method
  )

  # Run SuSiE workhorse
  model <- susie_workhorse(
    data, L,
    intercept = FALSE, standardize, scaled_prior_variance,
    residual_variance, data$prior_weights, data$null_weight,
    model_init, estimate_prior_variance, estimate_prior_method,
    check_null_threshold, estimate_residual_variance,
    estimate_residual_method,
    residual_variance_lowerbound, residual_variance_upperbound,
    max_iter, tol, verbose, track_fit, coverage, min_abs_corr,
    prior_tol, n_purity, compute_univariate_zscore = FALSE,
    check_prior, data$convergence_method
  )

  return(model)
}

susie_rss <- function(z = NULL, R, n = NULL,
                      bhat = NULL, shat = NULL, var_y = NULL,
                      L = min(10, ncol(R)),
                      lambda = 0,
                      maf = NULL,
                      maf_thresh = 0,
                      z_ld_weight = 0,
                      prior_variance = 50,
                      scaled_prior_variance = 0.2,
                      residual_variance = NULL,
                      prior_weights = NULL,
                      null_weight = 0,
                      standardize = TRUE,
                      intercept_value = 0,
                      estimate_residual_variance = FALSE,
                      estimate_residual_method = c("MLE", "MoM"),
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = c("optim", "EM", "simple"),
                      unmappable_effects = c("none", "inf"),
                      check_null_threshold = 0,
                      prior_tol = 1e-9,
                      residual_variance_lowerbound = 0,
                      residual_variance_upperbound = Inf,
                      s_init = NULL,
                      coverage = 0.95,
                      min_abs_corr = 0.5,
                      max_iter = 100,
                      tol = 1e-3,
                      convergence_method = c("elbo", "pip"),
                      verbose = FALSE,
                      track_fit = FALSE,
                      check_input = FALSE,
                      check_prior = TRUE,
                      check_R = TRUE,
                      check_z = FALSE,
                      n_purity = 100,
                      r_tol = 1e-8) {
  
  # Validate method arguments
  unmappable_effects       <- match.arg(unmappable_effects)
  estimate_prior_method    <- match.arg(estimate_prior_method)
  estimate_residual_method <- match.arg(estimate_residual_method)
  convergence_method       <- match.arg(convergence_method)

  # Construct data object
  data <- summary_stats_constructor(
    z, R, n, bhat, shat, var_y, lambda,
    maf, maf_thresh, z_ld_weight,
    prior_weights, null_weight, unmappable_effects,
    standardize, check_input,
    check_R, check_z, r_tol,
    prior_variance, scaled_prior_variance, intercept_value,
    estimate_residual_variance, estimate_residual_method,
    convergence_method, check_prior, residual_variance_upperbound
  )

  # Run SuSiE workhorse
  model <- susie_workhorse(
    data, L,
    intercept = FALSE, data$standardize, data$scaled_prior_variance,
    residual_variance, data$prior_weights, data$null_weight,
    s_init, estimate_prior_variance, estimate_prior_method,
    check_null_threshold, estimate_residual_variance,
    estimate_residual_method, residual_variance_lowerbound,
    data$residual_variance_upperbound, max_iter, tol, verbose,
    track_fit, coverage, min_abs_corr, prior_tol, n_purity,
    compute_univariate_zscore = FALSE, data$check_prior,
    convergence_method
  )

  return(model)
}
