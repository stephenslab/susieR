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

  # Construct data and params objects with ALL parameters
  susie_objects <- individual_data_constructor(
    X, y, L, scaled_prior_variance, residual_variance,
    prior_weights, null_weight, standardize, intercept,
    estimate_residual_variance, estimate_residual_method,
    estimate_prior_variance, estimate_prior_method,
    unmappable_effects, check_null_threshold, prior_tol,
    residual_variance_upperbound, model_init, coverage,
    min_abs_corr, compute_univariate_zscore, na.rm,
    max_iter, tol, convergence_method, verbose, track_fit,
    residual_variance_lowerbound, refine, n_purity,
    alpha0, beta0
  )
  data <- susie_objects$data
  params <- susie_objects$params

  model <- susie_workhorse(data, params)

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

  # Construct data and params objects
  susie_objects <- sufficient_stats_constructor(
    XtX, Xty, yty, n, L, X_colmeans, y_mean,
    maf = NULL, maf_thresh = 0, check_input, r_tol, standardize,
    scaled_prior_variance, residual_variance, prior_weights, null_weight,
    model_init, estimate_residual_variance, estimate_residual_method,
    residual_variance_lowerbound, residual_variance_upperbound,
    estimate_prior_variance, estimate_prior_method, unmappable_effects,
    check_null_threshold, prior_tol, max_iter, tol, convergence_method,
    coverage, min_abs_corr, n_purity, verbose, track_fit, check_prior
  )
  data <- susie_objects$data
  params <- susie_objects$params

  # Run SuSiE workhorse
  model <- susie_workhorse(data, params)

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

  # Construct data and params objects with ALL parameters
  susie_objects <- summary_stats_constructor(
    z = z, R = R, n = n, bhat = bhat, shat = shat, var_y = var_y,
    L = L, lambda = lambda, maf = maf, maf_thresh = maf_thresh,
    z_ld_weight = z_ld_weight, prior_variance = prior_variance,
    scaled_prior_variance = scaled_prior_variance,
    residual_variance = residual_variance, prior_weights = prior_weights,
    null_weight = null_weight, standardize = standardize,
    intercept_value = intercept_value,
    estimate_residual_variance = estimate_residual_variance,
    estimate_residual_method = estimate_residual_method,
    estimate_prior_variance = estimate_prior_variance,
    estimate_prior_method = estimate_prior_method,
    unmappable_effects = unmappable_effects,
    check_null_threshold = check_null_threshold, prior_tol = prior_tol,
    residual_variance_lowerbound = residual_variance_lowerbound,
    residual_variance_upperbound = residual_variance_upperbound,
    model_init = s_init, coverage = coverage, min_abs_corr = min_abs_corr,
    max_iter = max_iter, tol = tol, convergence_method = convergence_method,
    verbose = verbose, track_fit = track_fit, check_input = check_input,
    check_prior = check_prior, check_R = check_R, check_z = check_z,
    n_purity = n_purity, r_tol = r_tol,
    compute_univariate_zscore = compute_univariate_zscore
  )
  data <- susie_objects$data
  params <- susie_objects$params

  # Run SuSiE workhorse using new pattern for RSS
  model <- susie_workhorse(data, params)

  return(model)
}
