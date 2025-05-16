susie_new = function (X, y, L = min(10,ncol(X)),
                  scaled_prior_variance = 0.2,
                  residual_variance = NULL,
                  prior_weights = NULL,
                  null_weight = 0,
                  standardize = TRUE,
                  intercept = TRUE,
                  estimate_residual_variance = TRUE,
                  estimate_prior_variance = TRUE,
                  estimate_prior_method = c("optim", "EM", "simple"),
                  check_null_threshold = 0,
                  prior_tol = 1e-9,
                  residual_variance_upperbound = Inf,
                  model_init = NULL,
                  coverage = 0.95,
                  min_abs_corr = 0.5,
                  #median_abs_corr = NULL,
                  compute_univariate_zscore = FALSE,
                  na.rm = FALSE,
                  max_iter = 100,
                  tol = 1e-3,
                  verbose = FALSE,
                  track_fit = FALSE,
                  residual_variance_lowerbound = var(drop(y))/1e4,
                  refine = FALSE,
                  n_purity = 100) {

  # Estimate Prior Variance Method
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Handle null weight
  if (is.numeric(null_weight) && null_weight == 0)
    null_weight <- NULL

  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight) || null_weight <= 0 || null_weight >= 1)
      stop("null_weight must be between 0 and 1 (exclusive)")

    if (is.null(prior_weights))
      prior_weights <- rep(1 / ncol(X), ncol(X))

    # rescale existing prior_weights and append null weight
    prior_weights <- c(prior_weights * (1 - null_weight), null_weight)

    # add the extra 0 column to X
    X <- cbind(X, 0)
  }

  data <- susie_constructor(X, y, intercept, standardize, na.rm)

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
                        prior_tol, n_purity, compute_univariate_zscore)

  return(model)
}

susie_ss_new = function(XtX, Xty, yty, n,
                        X_colmeans       = NA,
                                    y_mean           = NA,
                                    maf              = NULL,
                                    maf_thresh       = 0,
                                    check_input      = FALSE,
                                    r_tol            = 1e-8,
                                    standardize      = TRUE,
                                    L                        = min(10, ncol(XtX)),
                                    scaled_prior_variance    = 0.2,
                                    residual_variance        = NULL,
                                    prior_weights            = NULL,
                                    null_weight              = 0,
                                    estimate_residual_variance = TRUE,
                                    estimate_prior_variance    = TRUE,
                                    estimate_prior_method      = c("optim","EM","simple"),
                                    check_null_threshold      = 0,
                                    prior_tol                 = 1e-9,
                                    max_iter                  = 100,
                                    tol                       = 1e-3,
                                    refine                    = FALSE,
                                    coverage                  = 0.95,
                                    min_abs_corr              = 0.5,
                                    median_abs_corr           = NULL,
                                    n_purity                  = 100,
                                    verbose                   = FALSE,
                                    track_fit                 = FALSE,
                                    ...) {

  # Estimate prior variance method
  estimate_prior_method <- match.arg(estimate_prior_method)

  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (is.null(prior_weights))
      prior_weights = c(rep(1/ncol(XtX)*(1-null_weight),ncol(XtX)),null_weight)
    else
      prior_weights = c(prior_weights*(1 - null_weight),null_weight)
    XtX = cbind(rbind(XtX,0),0)
    Xty = c(Xty,0)
    if (length(X_colmeans) == 1)
      X_colmeans = rep(X_colmeans,p)
    if (length(X_colmeans) != p)
      stop("The length of X_colmeans does not agree with number of variables")
  }

  data <- susie_ss_constructor(XtX, Xty, yty, n,
                               X_colmeans  = X_colmeans,
                               y_mean      = y_mean,
                               maf         = maf,
                               maf_thresh  = maf_thresh,
                               standardize = standardize,
                               r_tol       = r_tol,
                               check_input = check_input)

  # need to add engine with init fit etc.

  # model <- ibss_initialize(data, L = 10, scaled_prior_variance = 0.2,
  #                 residual_variance = NULL, prior_weights = NULL,
  #                 null_weight = 0, model_init = NULL)  # FOR TESTING PURPOSES


}


# Individual data works! --> now verify ss
# dat <- readRDS("/Users/alexmccreight/Columbia/Research/SuSiE-ASH/simple_data.rds")
# res <- susie_new(X = dat$X, y = dat$y, L = 10)
# str(res)

# dat <- readRDS("/Users/alexmccreight/Columbia/Research/SuSiE-ASH/simple_data.rds")
# X <- dat$X
# y <- dat$y
# XtX <- crossprod(X)
# Xty <- crossprod(X, y)
# yty <- as.numeric(crossprod(y))
# n   <- length(y)
