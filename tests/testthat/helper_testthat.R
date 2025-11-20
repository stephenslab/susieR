# =============================================================================
# HELPER FUNCTIONS FOR UNIT TESTS
# =============================================================================
#
# This file provides helper functions for testing the susieR package. These
# functions are automatically loaded by testthat when running tests.
#
# CONTENTS:
#
# 1. DATA SIMULATION FUNCTIONS
#    - simulate()              : Legacy simulation (sparse/dense matrices)
#    - simulate_tf()           : Simulate trend filtering data
#    - simulate_regression()   : Simulate linear regression with causal effects
#
# 2. DATA SETUP FUNCTIONS (Constructor-based)
#    - setup_individual_data()  : Create 'individual' class test data
#    - setup_ss_data()          : Create 'ss' class test data
#    - setup_rss_lambda_data()  : Create 'rss_lambda' class test data
#
# 3. CUSTOM EXPECTATION FUNCTIONS
#    - expect_equal_susie_*()   : Compare susie objects (individual/ss/rss)
#    - expect_equal_SER_*()     : Compare single effect regression results
#
# 4. UTILITY FUNCTIONS
#    - set_X_attributes()       : Set standardization attributes on X matrix
#    - compute_summary_stats()  : Compute XtX, Xty, yty from X, y
#    - create_sparsity_mat()    : Create sparse matrix with given sparsity
#
# USAGE NOTES:
#
# - All simulation functions use set.seed() internally for reproducibility
# - Setup functions return list(data, params, model) ready for testing
# - Custom expectation functions handle tolerance and class-specific fields
# - For new tests, prefer simulate_regression() over legacy simulate()
# - Setup functions use actual constructors to ensure correct initialization
#
# =============================================================================

# -----------------------------------------------------------------------------
# UTILITY FUNCTIONS
# -----------------------------------------------------------------------------

#' Set standardization attributes on matrix X
#'
#' Sets three attributes on the input matrix: `scaled:center` (column means),
#' `scaled:scale` (column standard deviations), and `d` (column sums of
#' squared standardized values). These attributes are used by SuSiE algorithm
#' to efficiently handle standardized data.
#'
#' @param X An n by p data matrix (dense, sparse, or trend filtering matrix)
#' @param center Logical; if TRUE, center by column means
#' @param scale Logical; if TRUE, scale by column standard deviations
#' @return X with three attributes set: `scaled:center`, `scaled:scale`, and `d`
#' @keywords internal
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
set_X_attributes <- function(X, center = TRUE, scale = TRUE) {
    
  # if X is a trend filtering matrix
  if (!is.null(attr(X,"matrix.type"))) {
    order = attr(X,"order")
    n = ncol(X)
    
    # Set three attributes for X.
    attr(X,"scaled:center") = compute_tf_cm(order,n)
    attr(X,"scaled:scale") = compute_tf_csd(order,n)
    attr(X,"d") = compute_tf_d(order,n,attr(X,"scaled:center"),
                               attr(X,"scaled:scale"),scale,center)
    if (!center)
      attr(X,"scaled:center") = rep(0,n)
    if (!scale)
      attr(X,"scaled:scale") = rep(1,n)
  } else {
      
    # If X is either a dense or sparse ordinary matrix.
    # Get column means.
    cm = colMeans(X,na.rm = TRUE)
    
    # Get column standard deviations.
    csd = compute_colSds(X)
    
    # Set sd = 1 when the column has variance 0.
    csd[csd == 0] = 1
    if (!center)
      cm = rep(0,length = length(cm))
    if (!scale) 
      csd = rep(1,length = length(cm))

    # Ah, this code is very inefficient because the matrix becomes
    # dense!
    # 
    #   X.std = as.matrix(X)
    #   X.std = (t(X.std) - cm)/csd
    #   attr(X,"d") = rowSums(X.std * X.std)
    #
    # Set three attributes for X.
    n = nrow(X)
    d = n*colMeans(X)^2 + (n-1)*compute_colSds(X)^2
    d = (d - n*cm^2)/csd^2
    attr(X,"d") = d
    attr(X,"scaled:center") = cm
    attr(X,"scaled:scale") = csd
  }
  return(X)
}

#' Create sparse matrix with specified sparsity level
#'
#' Generates a binary matrix with a specified proportion of non-zero entries.
#' Used for testing sparse matrix operations.
#'
#' @param sparsity Proportion of zero entries (between 0 and 1)
#' @param n Number of rows
#' @param p Number of columns
#' @return Binary matrix with (1-sparsity)*n*p non-zero entries
#' @keywords internal
create_sparsity_mat <- function(sparsity, n, p) {
  nonzero <- round(n * p * (1 - sparsity))
  nonzero.idx <- sample(n * p, nonzero)
  mat <- numeric(n * p)
  mat[nonzero.idx] <- 1
  mat <- matrix(mat, nrow = n, ncol = p)
  return(mat)
}

# -----------------------------------------------------------------------------
# DATA SIMULATION FUNCTIONS
# -----------------------------------------------------------------------------

#' Simulate trend filtering data
#'
#' Generates synthetic data for testing trend filtering functionality.
#' Creates piecewise constant (order=0), linear (order=1), or quadratic
#' (order=2) signals with noise.
#'
#' @param order Trend filtering order (0, 1, or 2)
#' @return List with X (trend filtering matrix) and y (response vector)
#' @keywords internal
simulate_tf <- function(order) {
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(2)
  n = 50
  D = diag(-1, n)
  for (i in 1:(n-1)){
    D[i, i+1] = 1
  }
  if (order==0) {
    beta = c(rep(0,5),rep(1,5),rep(3,5),rep(-2,5),rep(0,30))
    y = beta + rnorm(n)
    X = solve(D)
  } else if (order==1) {
    beta = numeric(n)
    for (i in 1:n){
      if (i <= 5){
        beta[i] = 0.001*i + 2
      } else if (i <= 15){
        beta[i] = 5*0.001*i + 1.6
      } else{
        beta[i] = 6.1 - 10*0.001*i
      }
    }
    y = beta + rnorm(n)
    X = solve(D%*%D)
  } else if (order==2) {
    beta = numeric(n)
    for (i in 1:n){
      if (i <= 5){
        beta[i] = (0.001*i)^2
      } else if (i <= 35){
        beta[i] = -5*(0.001*i)^2 + 0.06
      } else{
        beta[i] = 3*(0.001*i)^2 - 3.86
      }
    }
    y = beta + rnorm(n)
    X = solve(D%*%D%*%D)
  }
  return(list(X=X, y=y))
}

# -----------------------------------------------------------------------------
# CUSTOM EXPECTATION FUNCTIONS
# -----------------------------------------------------------------------------

expect_equal_susie_update = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal(new.res$alpha, original.res$alpha, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu, original.res$mu, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu2, original.res$mu2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$Xr, original.res$Xr, scale = 1, tolerance = tolerance)
  expect_equal(new.res$KL, original.res$KL, scale = 1, tolerance = tolerance)
  expect_equal(new.res$sigma2, original.res$sigma2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$V, original.res$V, scale = 1, tolerance = tolerance)
}

expect_equal_susie_suff_stat_update = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal(new.res$alpha, original.res$alpha, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu, original.res$mu, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu2, original.res$mu2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$XtXr, original.res$XtXr, scale = 1, tolerance = tolerance)
  expect_equal(new.res$KL, original.res$KL, scale = 1, tolerance = tolerance)
  expect_equal(new.res$sigma2, original.res$sigma2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$V, original.res$V, scale = 1, tolerance = tolerance)
}

expect_equal_susie_rss_update = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal(new.res$alpha, original.res$alpha, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu, original.res$mu, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu2, original.res$mu2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$Rz, original.res$Rz, scale = 1, tolerance = tolerance)
  expect_equal(new.res$KL, original.res$KL, scale = 1, tolerance = tolerance)
  expect_equal(new.res$sigma2, original.res$sigma2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$V, original.res$V, scale = 1, tolerance = tolerance)
}

expect_equal_SER = function(new.res, original.res){
  expect_equal(new.res$alpha, original.res$alpha)
  expect_equal(new.res$mu, original.res$mu)
  expect_equal(new.res$mu2, original.res$mu2)
  expect_equal(new.res$lbf, original.res$lbf)
  expect_equal(new.res$V, original.res$V)
  expect_equal(new.res$loglik, original.res$loglik)
}

expect_equal_SER_suff_stat = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal(new.res$alpha, original.res$alpha, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu, original.res$mu, scale = 1, tolerance = tolerance)
  expect_equal(new.res$mu2, original.res$mu2, scale = 1, tolerance = tolerance)
  expect_equal(new.res$lbf, original.res$lbf, scale = 1, tolerance = tolerance)
  expect_equal(new.res$V, original.res$V, scale = 1, tolerance = tolerance)
  expect_equal(new.res$lbf_model, original.res$lbf_model, scale = 1, tolerance = tolerance)
}

expect_equal_susie = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal_susie_update(new.res, original.res, tolerance = tolerance)
  expect_equal(new.res$elbo, original.res$elbo, scale = 1, tolerance = tolerance)
  expect_equal(new.res$niter, original.res$niter, scale = 1, tolerance = tolerance)
  expect_equal(new.res$intercept, original.res$intercept, scale = 1, tolerance = tolerance)
  expect_equal(new.res$fitted, original.res$fitted, scale = 1, tolerance = tolerance)
  expect_equal(new.res$X_column_scale_factors, original.res$X_column_scale_factors, scale = 1, tolerance = tolerance)
}

expect_equal_susie_suff_stat = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal_susie_suff_stat_update(new.res, original.res, tolerance = tolerance)
  expect_equal(new.res$elbo, original.res$elbo, scale = 1, tolerance = tolerance)
  expect_equal(new.res$niter, original.res$niter, scale = 1, tolerance = tolerance)
  expect_equal(new.res$intercept, original.res$intercept, scale = 1, tolerance = tolerance)
  expect_equal(new.res$Xtfitted, original.res$Xtfitted, scale = 1, tolerance = tolerance)
}

expect_equal_susie_rss = function(new.res, original.res, tolerance = .Machine$double.eps^0.5){
  expect_equal_susie_rss_update(new.res, original.res, scale = 1, tolerance = tolerance)
  expect_equal(new.res$elbo, original.res$elbo, scale = 1, tolerance = tolerance)
  expect_equal(new.res$niter, original.res$niter, scale = 1, tolerance = tolerance)
  expect_equal(new.res$intercept, original.res$intercept, scale = 1, tolerance = tolerance)
  expect_equal(new.res$Rz, original.res$Rz, scale = 1, tolerance = tolerance)
}

#' Unified dispatcher for comparing susie objects
#'
#' Automatically detects the type of susie object and calls the appropriate
#' comparison function. This simplifies test code and ensures correct
#' comparison based on object structure.
#'
#' @param new.res New susie result object
#' @param original.res Original susie result object to compare against
#' @param tolerance Numerical tolerance for comparisons (default: sqrt(.Machine$double.eps))
#'
#' @details
#' Detects object type by checking for class-specific fields:
#' - Individual data: has 'Xr' field
#' - Sufficient stats: has 'XtXr' field
#' - RSS/RSS lambda: has 'Rz' field
#'
#' @examples
#' # Automatically handles all susie object types
#' fit1 <- susie(X, y, L = 5)
#' fit2 <- susie(X, y, L = 5)
#' expect_equal_susie_objects(fit1, fit2)
#'
expect_equal_susie_objects <- function(new.res, original.res,
                                       tolerance = .Machine$double.eps^0.5) {
  # Detect type based on class-specific fields
  if (!is.null(new.res$Xr) && !is.null(original.res$Xr)) {
    # Individual data (has Xr)
    expect_equal_susie(new.res, original.res, tolerance = tolerance)
  } else if (!is.null(new.res$XtXr) && !is.null(original.res$XtXr)) {
    # Sufficient stats (has XtXr)
    expect_equal_susie_suff_stat(new.res, original.res, tolerance = tolerance)
  } else if (!is.null(new.res$Rz) && !is.null(original.res$Rz)) {
    # RSS or RSS lambda (has Rz)
    expect_equal_susie_rss(new.res, original.res, tolerance = tolerance)
  } else {
    stop("Cannot determine susie object type. Unknown structure.")
  }
}

#' Unified dispatcher for comparing SER results
#'
#' Automatically detects the type of single effect regression result and calls
#' the appropriate comparison function.
#'
#' @param new.res New SER result object
#' @param original.res Original SER result object to compare against
#' @param tolerance Numerical tolerance for comparisons (default: sqrt(.Machine$double.eps))
#'
#' @details
#' Detects result type by checking for 'lbf_model' field (sufficient stats only)
#'
expect_equal_SER_objects <- function(new.res, original.res,
                                     tolerance = .Machine$double.eps^0.5) {
  # Detect type based on presence of lbf_model (sufficient stats specific)
  if (!is.null(new.res$lbf_model) && !is.null(original.res$lbf_model)) {
    # Sufficient stats SER
    expect_equal_SER_suff_stat(new.res, original.res, tolerance = tolerance)
  } else {
    # Individual data SER (default)
    expect_equal_SER(new.res, original.res)
  }
}

# -----------------------------------------------------------------------------
# BASE HELPER FUNCTIONS (Internal - reduce duplication)
# -----------------------------------------------------------------------------

#' Generate base regression data for testing
#'
#' Creates random X matrix and y vector for use in test setup functions.
#' This function encapsulates the common data generation pattern used across
#' multiple setup functions to reduce code duplication.
#'
#' @param n Sample size
#' @param p Number of variables
#' @param k Number of causal variables (if 0, generates random y)
#' @param signal_sd Standard deviation of effect sizes for causal variables
#' @param seed Random seed (if NULL, no seed is set)
#' @return List with X, y, and optionally beta and causal_idx
#' @keywords internal
generate_base_data <- function(n, p, k = 0, signal_sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)

  if (k > 0) {
    # Generate data with known causal structure
    beta <- rep(0, p)
    causal_idx <- sort(sample(1:p, k))
    beta[causal_idx] <- rnorm(k, mean = 0, sd = signal_sd)
    y <- as.vector(X %*% beta + rnorm(n))
    return(list(X = X, y = y, n = n, p = p, beta = beta, causal_idx = causal_idx))
  } else {
    # Generate random y (no causal structure)
    y <- rnorm(n)
    return(list(X = X, y = y, n = n, p = p))
  }
}

#' Create base model structure
#'
#' Creates the common model list structure used across all data types.
#' This function encapsulates the shared model initialization pattern.
#'
#' @param L Number of single effects
#' @param p Number of variables
#' @param n Number of samples (for individual data, adds Xr field)
#' @param X_attr Optional attributes from X (for predictor_weights)
#' @return List with alpha, mu, mu2, V, sigma2, pi, lbf, lbf_variable, KL
#' @keywords internal
create_base_model <- function(L, p, n = NULL, X_attr = NULL) {
  model <- list(
    alpha = matrix(1 / p, L, p),
    mu = matrix(0, L, p),
    mu2 = matrix(0, L, p),
    V = rep(1, L),
    sigma2 = 1,
    pi = rep(1 / p, p),
    lbf = rep(0, L),
    lbf_variable = matrix(0, L, p),
    KL = rep(0, L),
    null_weight = 0
  )

  # Add class-specific fields
  if (!is.null(X_attr)) {
    model$predictor_weights <- X_attr
  }

  if (!is.null(n)) {
    model$Xr <- rep(0, n)  # For individual data
  }

  return(model)
}

#' Create standard parameter list
#'
#' Creates the default params list used by setup functions. Provides consistent
#' defaults that can be overridden by specific setup functions.
#'
#' @param L Number of single effects
#' @param p Number of variables
#' @param unmappable_effects One of "none", "inf", or "ash"
#' @param additional_params Named list of additional/override parameters
#' @return List of parameters for SuSiE fitting
#' @keywords internal
create_base_params <- function(L, p, unmappable_effects = "none",
                               additional_params = list()) {
  params <- list(
    L = L,
    intercept = TRUE,
    standardize = TRUE,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    unmappable_effects = unmappable_effects,
    use_servin_stephens = FALSE,
    compute_univariate_zscore = TRUE,
    coverage = 0.95,
    min_abs_corr = 0.5,
    n_purity = 100,
    check_null_threshold = 0.1,
    scaled_prior_variance = 0.2,
    prior_weights = rep(1 / p, p),
    null_weight = 0,
    residual_variance = NULL,
    track_fit = FALSE,
    prior_tol = 1e-9,
    max_iter = 100,
    tol = 1e-3,
    convergence_method = "elbo",
    verbose = FALSE,
    refine = FALSE,
    model_init = NULL
  )

  # Override with additional params if provided
  if (length(additional_params) > 0) {
    params[names(additional_params)] <- additional_params
  }

  return(params)
}

# -----------------------------------------------------------------------------
# DATA SETUP FUNCTIONS (Constructor-based)
# -----------------------------------------------------------------------------

#' Setup individual-level data for testing
#'
#' Creates a complete test setup with individual-level data (X, y matrices),
#' parameters, and an initialized model. This is the primary setup function
#' for testing individual data methods.
#'
#' @param n Sample size
#' @param p Number of variables
#' @param L Number of single effects
#' @param seed Random seed for reproducibility
#' @return List with data (class: individual), params, and model
#' @keywords internal
#' @examples
#' # Internal use in tests
#' setup <- setup_individual_data(n = 100, p = 50, L = 5)
#' fit <- susie(setup$data$X, setup$data$y, L = setup$params$L)
setup_individual_data <- function(n = 100, p = 50, L = 5, seed = 42) {
  # Use base helper for data generation
  base_data <- generate_base_data(n, p, k = 0, seed = seed)
  X <- base_data$X
  y <- base_data$y

  X <- set_X_attributes(X, center = TRUE, scale = TRUE)
  mean_y <- mean(y)
  y <- y - mean_y

  data <- structure(
    list(
      X = X,
      y = y,
      n = n,
      p = p,
      mean_y = mean_y
    ),
    class = "individual"
  )

  # Use base helper for standard params
  params <- create_base_params(L, p, unmappable_effects = "none")

  # Use base helper for model, then add individual-specific fields
  model <- create_base_model(L, p, n = n, X_attr = attr(X, "d"))

  list(data = data, params = params, model = model)
}

#' Setup sufficient statistics data with unmappable_effects support
#'
#' Creates a complete test setup with sufficient statistics (XtX, Xty, yty),
#' parameters, and an initialized model. Supports unmappable effects testing.
#'
#' @param n Number of samples
#' @param p Number of variables
#' @param L Number of single effects
#' @param seed Random seed
#' @param unmappable_effects One of "none" or "inf"
#' @return List with data (class: ss), params, and model
#' @keywords internal
setup_ss_data <- function(n = 100, p = 50, L = 5, seed = 42,
                          unmappable_effects = "none") {
  # Use base helper for data generation
  base_data <- generate_base_data(n, p, k = 0, seed = seed)
  X <- base_data$X
  y <- base_data$y

  # Center and scale X like the constructor does
  X_colmeans <- colMeans(X)
  X <- sweep(X, 2, X_colmeans)
  y_mean <- mean(y)
  y <- y - y_mean

  # Compute sufficient statistics
  XtX <- crossprod(X)
  Xty <- as.vector(crossprod(X, y))
  yty <- sum(y^2)

  # Use the actual constructor like susie_ss does
  # This ensures proper setup including eigen decomposition for unmappable effects
  susie_objects <- sufficient_stats_constructor(
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    L = L,
    X_colmeans = X_colmeans,
    y_mean = y_mean,
    standardize = TRUE,
    unmappable_effects = unmappable_effects,
    residual_variance = 1,  # Set initial residual variance
    estimate_residual_method = if (unmappable_effects != "none") "MoM" else "MLE",
    convergence_method = if (unmappable_effects != "none") "pip" else "elbo",
    coverage = 0.95,
    min_abs_corr = 0.5,
    n_purity = 100,
    check_prior = FALSE,
    track_fit = FALSE
  )

  data <- susie_objects$data
  params <- susie_objects$params

  # Use base helper for model, add ss-specific fields
  model <- create_base_model(L, data$p, n = NULL, X_attr = attr(data$XtX, "d"))
  model$XtXr <- rep(0, data$p)  # ss-specific field

  # Add unmappable components if needed
  if (unmappable_effects == "inf") {
    model$tau2 <- 0
    model$theta <- rep(0, data$p)
  }

  list(data = data, params = params, model = model)
}

#' Setup RSS lambda test data
#'
#' Creates a complete test setup for RSS with correlated residuals (lambda > 0).
#' Generates data with known causal structure, computes z-scores and correlation
#' matrix, and initializes model using the rss_lambda constructor.
#'
#' @param n Number of samples
#' @param p Number of variables
#' @param k Number of causal variables
#' @param lambda Lambda parameter (residual correlation, between 0 and 1)
#' @param signal_sd Standard deviation of causal effects
#' @param seed Random seed
#' @param L Number of single effects
#' @return List with X, y, beta, causal_idx, z, R, n, p, k, lambda, data, params, model
#' @keywords internal
setup_rss_lambda_data <- function(n = 500, p = 50, k = 3, lambda = 0.5,
                                  signal_sd = 0.5, seed = NULL, L = 5) {
  # Use base helper for data generation with causal structure
  base_data <- generate_base_data(n, p, k = k, signal_sd = signal_sd, seed = seed)
  X <- base_data$X
  y <- base_data$y
  beta <- base_data$beta
  causal_idx <- base_data$causal_idx

  # Compute sufficient statistics and z-scores
  input_ss <- compute_suff_stat(X, y, standardize = TRUE)
  ss <- univariate_regression(X, y)
  R <- with(input_ss, cov2cor(XtX))
  R <- (R + t(R)) / 2
  z <- with(ss, betahat / sebetahat)

  # Build data and params using constructor
  constructor_result <- rss_lambda_constructor(z = z, R = R, lambda = lambda, n = n, L = L)

  # Initialize model properly
  var_y <- get_var_y.rss_lambda(constructor_result$data)
  model <- initialize_susie_model.rss_lambda(constructor_result$data, constructor_result$params, var_y)

  list(
    X = X,
    y = y,
    beta = beta,
    causal_idx = causal_idx,
    z = z,
    R = R,
    n = n,
    p = p,
    k = k,
    lambda = lambda,
    data = constructor_result$data,
    params = constructor_result$params,
    model = model
  )
}

# -----------------------------------------------------------------------------
# ADDITIONAL SIMULATION FUNCTIONS
# -----------------------------------------------------------------------------

#' Simulate simple regression data with known causal variables
#'
#' @param n Sample size
#' @param p Number of variables
#' @param k Number of causal variables
#' @param signal_sd Standard deviation of effect sizes
#' @param noise_sd Standard deviation of noise
#' @param center Whether to center X and y
#' @param scale Whether to scale X to unit variance
#' @return List with X, y, beta, causal_idx
simulate_regression <- function(n = 100, p = 50, k = 3,
                               signal_sd = 1, noise_sd = 1,
                               center = TRUE, scale = TRUE) {
  # Generate independent X
  X <- matrix(rnorm(n * p), n, p)

  # Optionally standardize X
  if (center || scale) {
    X <- scale(X, center = center, scale = scale)
  }

  # Generate causal effects
  beta <- rep(0, p)
  causal_idx <- sort(sample(1:p, k))
  beta[causal_idx] <- rnorm(k, mean = 0, sd = signal_sd)

  # Generate y
  y <- drop(X %*% beta + rnorm(n, sd = noise_sd))

  # Optionally center y
  if (center) {
    y <- y - mean(y)
  }

  list(
    X = X,
    y = y,
    beta = beta,
    causal_idx = causal_idx,
    n = n,
    p = p,
    k = k
  )
}

#' Compute summary statistics (XtX, Xty, yty) from X and y
#'
#' @param X n x p matrix
#' @param y n vector
#' @return List with XtX, Xty, yty, n
compute_summary_stats <- function(X, y) {
  list(
    XtX = crossprod(X),
    Xty = drop(crossprod(X, y)),
    yty = sum(y^2),
    n = length(y)
  )
}

#' Create model with credible sets for refinement testing
#'
#' Generates synthetic data with known causal structure, fits SuSiE model,
#' and returns both the model and the data/params objects needed for
#' refinement testing. Used primarily by test_refinement.R.
#'
#' @param n Sample size
#' @param p Number of variables
#' @param L Number of single effects
#' @param n_causal Number of causal variables to simulate
#' @param seed Random seed for reproducibility
#' @param run_susie Logical; if TRUE, runs susie and returns model, otherwise just returns data
#' @return List with model, data, params, X, y, beta, causal_idx
#' @keywords internal
create_model_with_cs <- function(n = 100, p = 50, L = 5, n_causal = 3,
                                  seed = 42, run_susie = TRUE) {
  set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = TRUE)

  beta <- rep(0, p)
  causal_idx <- sample(1:p, n_causal)
  beta[causal_idx] <- rnorm(n_causal, sd = 1)

  y <- as.vector(X %*% beta + rnorm(n, sd = 0.5))

  if (run_susie) {
    model <- susie(X, y, L = L, verbose = FALSE)

    constructor_result <- individual_data_constructor(
      X = X, y = y, L = L,
      standardize = TRUE, intercept = TRUE,
      estimate_residual_method = "MLE",
      convergence_method = "elbo",
      coverage = 0.95, min_abs_corr = 0.5,
      n_purity = 100,
      track_fit = FALSE
    )

    return(list(
      model = model,
      data = constructor_result$data,
      params = constructor_result$params,
      X = X,
      y = y,
      beta = beta,
      causal_idx = causal_idx
    ))
  } else {
    return(list(
      X = X,
      y = y,
      beta = beta,
      causal_idx = causal_idx
    ))
  }
}
