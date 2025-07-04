susie_new <- function(X, y, L = min(10,ncol(X)),
                  scaled_prior_variance = 0.2,
                  residual_variance = NULL,
                  prior_weights = NULL,
                  null_weight = 0,
                  standardize = TRUE,
                  intercept = TRUE,
                  estimate_residual_variance = TRUE,
                  estimate_prior_variance = TRUE,
                  estimate_prior_method = c("optim", "EM", "simple"),
                  non_sparse_method = c("none", "inf", "ash"),
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
                  verbose = FALSE,
                  track_fit = FALSE,
                  residual_variance_lowerbound = var(drop(y))/1e4,
                  refine = FALSE,
                  n_purity = 100,
                  check_prior = FALSE) {

  # Estimate Non-sparse Effects Method
  non_sparse_method <- match.arg(non_sparse_method)

  # Construct Data Object
  data <- susie_constructor(X, y, intercept, standardize,
                            na.rm, prior_weights, null_weight,
                            non_sparse_method = non_sparse_method)

  # Estimate Prior Variance Method
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Run SuSiE. TODO: Another option would be to store model inference parameters
  # like the prior estimation or residual/inf estimation method in a list in the
  # model object for later.
  model <- susie_engine(data, L, intercept, standardize, scaled_prior_variance,
                        residual_variance, prior_weights, null_weight,
                        model_init, estimate_prior_variance,
                        estimate_prior_method,
                        check_null_threshold,
                        estimate_residual_variance,
                        residual_variance_lowerbound,
                        residual_variance_upperbound,
                        max_iter, tol, verbose, track_fit,
                        coverage, min_abs_corr,
                        prior_tol, n_purity, compute_univariate_zscore,
                        check_prior, refine)

  return(model)
}

susie_ss <- function(XtX, Xty, yty, n,
                     L = min(10, ncol(XtX)),
                     X_colmeans = NA, y_mean = NA,
                     maf = NULL, maf_thresh = 0,
                     check_input                = FALSE,
                     r_tol                      = 1e-8,
                     standardize                = TRUE,
                     scaled_prior_variance      = 0.2,
                     residual_variance          = NULL,
                     prior_weights              = NULL,
                     null_weight                = 0,
                     model_init                 = NULL,
                     estimate_residual_variance = TRUE,
                     residual_variance_lowerbound = 0,
                     residual_variance_upperbound = Inf,
                     estimate_prior_variance    = TRUE,
                     estimate_prior_method      = c("optim","EM","simple"),
                     non_sparse_method          = c("none", "inf", "ash"),
                     check_null_threshold       = 0,
                     prior_tol                  = 1e-9,
                     max_iter                   = 100,
                     tol                        = 1e-3,
                     coverage                   = 0.95,
                     min_abs_corr               = 0.5,
                     n_purity                   = 100,
                     verbose                    = FALSE,
                     track_fit                  = FALSE,
                     check_prior                = FALSE,
                     refine                     = FALSE, ...) {

  # Estimate Non-sparse Effects Method
  non_sparse_method <- match.arg(non_sparse_method)

  # Create data object
  data <- susie_ss_constructor(XtX, Xty, yty, n,
                               X_colmeans    = X_colmeans,
                               y_mean        = y_mean,
                               maf           = maf,
                               maf_thresh    = maf_thresh,
                               standardize   = standardize,
                               r_tol         = r_tol,
                               check_input   = check_input,
                               prior_weights = prior_weights,
                               null_weight   = null_weight,
                               non_sparse_method = non_sparse_method)

  # Estimate Prior Variance Method
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Run SuSiE
  model <- susie_engine(data, L, intercept, standardize, scaled_prior_variance,
                        residual_variance, prior_weights, null_weight,
                        model_init, estimate_prior_variance,
                        estimate_prior_method,
                        check_null_threshold,
                        estimate_residual_variance,
                        residual_variance_lowerbound,
                        residual_variance_upperbound,
                        max_iter, tol, verbose, track_fit,
                        coverage, min_abs_corr,
                        prior_tol, n_purity,
                        compute_univariate_zscore, check_prior)


}
