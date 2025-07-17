# Individual-level data backend methods

# Initialize fitted values
initialize_fitted.individual <- function(data, alpha, mu) {
  return(list(Xr = compute_Xb(data$X, colSums(alpha * mu))))
}

# Initialize susie model
initialize_susie_model.individual <- function(data, L, scaled_prior_variance, var_y,
                                           residual_variance, prior_weights, ...) {
  return(initialize_matrices(data$p, L, scaled_prior_variance, var_y,
                                residual_variance, prior_weights))
}

# Get variance of y
get_var_y.individual <- function(data, ...) {
  return(var(drop(data$y)))
}

# Configure individual data for specified method
configure_data.individual <- function(data, non_sparse_method) {
  if (non_sparse_method == "none") {
    return(data)
  } else {
    # Convert to sufficient statistics for non-sparse methods
    warning("Individual-level data converted to sufficient statistics for non-sparse methods")
    return(convert_individual_to_ss_non_sparse(data, non_sparse_method))
  }
}

# Convert individual data to ss with non-sparse components.
convert_individual_to_ss_non_sparse <- function(individual_data, non_sparse_method) {

  # Extract components from individual data
  X <- individual_data$X
  y <- individual_data$y
  n <- individual_data$n
  p <- individual_data$p
  mean_y <- individual_data$mean_y

  # Compute sufficient statistics
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  yty <- sum(y^2)

  # Get column means and scaling from attributes
  X_colmeans <- attr(X, "scaled:center")

  # Create sufficient statistics data object with multiple classes
  class_vector <- c(paste0("ss_", non_sparse_method), "ss")

  ss_data <- structure(list(
    XtX        = XtX,
    Xty        = Xty,
    yty        = yty,
    n          = n,
    p          = p,
    X_colmeans = X_colmeans,
    y_mean     = mean_y),
    class = class_vector)

  # Copy attributes from X to XtX
  attr(ss_data$XtX, "d") <- attr(X, "d")
  attr(ss_data$XtX, "scaled:scale") <- attr(X, "scaled:scale")

  # Add eigen decomposition for non-sparse methods
  ss_data <- add_eigen_decomposition(ss_data)

  return(ss_data)
}

# Extract core parameters across iterations
extract_core.individual <- function(data, model, tracking, iter, track_fit, ...) {
  if (isTRUE(track_fit)) {
    tracking[[iter]] <- list(alpha = model$alpha,
                             niter = iter,
                             V = model$V,
                             sigma2 = model$sigma2)
  }
  return(tracking)
}

# Validate prior variance
validate_prior.individual <- function(data, model, check_prior, ...) {
  invisible(TRUE)
}

# Posterior expected log-likelihood for single effect regression
SER_posterior_e_loglik.individual <- function(data, model, R, Eb, Eb2) {
  return(-0.5 * data$n * log(2 * pi * model$sigma2) -
         0.5 / model$sigma2 * (sum(R * R) - 2 * sum(R * compute_Xb(data$X, Eb)) +
                               sum(attr(data$X, "d") * Eb2)))
}

# Expected squared residuals
get_ER2.individual <- function(data, model) {
  Xr_L <- compute_MXt(model$alpha * model$mu, data$X)
  postb2 <- model$alpha * model$mu2
  return(sum((data$y - model$Xr)^2) - sum(Xr_L^2) + sum(attr(data$X, "d") * t(postb2)))
}

# Single Effect Update
single_effect_update.individual <- function(
    data, model, l,
    optimize_V, check_null_threshold) {

  # Remove lth effect
  model$Xr <- model$Xr - compute_Xb(data$X, model$alpha[l, ] * model$mu[l, ])

  # Compute Residuals
  R <- data$y - model$Xr
  XtR <- crossprod(data$X, R)
  d <- attr(data$X, "d")

  res <- single_effect_regression(
    Xty                  = XtR,
    dXtX                 = d,
    V                    = model$V[l],
    residual_variance    = model$sigma2,
    prior_weights        = model$pi,
    optimize_V           = optimize_V,
    check_null_threshold = check_null_threshold)

  # log-likelihood term using current residual vector (not available in ss)
  res$loglik <- res$lbf_model +
    sum(dnorm(R, 0, sqrt(model$sigma2), log = TRUE))

  res$KL <- -res$loglik +
    SER_posterior_e_loglik(data, model, R,
                           Eb  = res$alpha * res$mu,
                           Eb2 = res$alpha * res$mu2)

  # Update alpha and mu for adding effect back
  model$alpha[l,]         <- res$alpha
  model$mu[l, ]           <- res$mu
  model$mu2[l, ]          <- res$mu2
  model$V[l]              <- res$V
  model$lbf[l]            <- res$lbf_model
  model$lbf_variable[l, ] <- res$lbf
  model$KL[l]             <- res$KL

  model$Xr <- model$Xr + compute_Xb(data$X, model$alpha[l, ] * model$mu[l, ])

  return(model)
}

# Get column scale factors
get_scale_factors.individual <- function(data){
  return(attr(data$X,"scaled:scale"))
}

# Get intercept
get_intercept.individual <- function(data, model, intercept){
  if (intercept) {
    return(data$mean_y - sum(attr(data$X, "scaled:center") *
                                 (colSums(model$alpha * model$mu) / attr(data$X, "scaled:scale"))))
  } else {
    return(0)
  }
}

# Get Fitted Values
get_fitted.individual <- function(data, model, intercept){
  if(intercept) {
    fitted = model$Xr + data$mean_y
  } else{
    fitted = model$Xr
  }

  fitted = drop(fitted)
  names(fitted) = `if`(is.null(names(data$y)),rownames(data$X),names(data$y))

  return(fitted)
}

# Get Credible Sets
get_cs.individual <- function(data, model, coverage, min_abs_corr, n_purity){

  if (is.null(coverage) || is.null(min_abs_corr)) return(NULL)

  return(susie_get_cs(model,coverage = coverage,X = data$X,
                      min_abs_corr = min_abs_corr,
                      n_purity = n_purity))
}


# Get Variable Names
get_variable_names.individual <- function(data, model, null_weight){
  variable_names <- colnames(data$X)
  return(assign_names(model, variable_names, null_weight, data$p))
}

# Get univariate z-score
get_zscore.individual <- function(data, model, compute_univariate_zscore,
                                  intercept, standardize, null_weight){

  if(isFALSE(compute_univariate_zscore)) return(NULL)

  X <- data$X

  if (!is.matrix(X))
    warning("Calculation of univariate regression z-scores is not ",
            "implemented specifically for sparse or trend filtering ",
            "matrices, so this step may be slow if the matrix is large; ",
            "to skip this step set compute_univariate_zscore = FALSE")
  if (!is.null(null_weight) && null_weight != 0)
    X = X[,1:(ncol(X) - 1)]

  return(calc_z(X,data$y,center = intercept,scale = standardize))

}

# Update variance components for individual data (standard approach)
update_variance_components.individual <- function(data, model) {
  sigma2 <- est_residual_variance(data, model)
  return(list(sigma2 = sigma2, tausq = NULL))
}

# Update derived quantities for individual data (no-op for standard)
update_derived_quantities.individual <- function(data, model) {
  return(data)  # No changes needed for individual data
}

# Check convergence for individual data (uses ELBO)
check_convergence.individual <- function(data, model_prev, model_current, elbo_prev, elbo_current, tol) {
  # Standard ELBO-based convergence (uses pre-computed ELBO values)
  # Handle case where elbo_prev is NA (first iteration)
  if (is.na(elbo_prev)) {
    return(FALSE)  # Cannot converge on first iteration
  }
  return(elbo_current - elbo_prev < tol)
}

# Update variance before convergence check for individual data
update_variance_before_convergence.individual <- function(data) {
  # Standard behavior: Check convergence first, then update variance
  return(FALSE)
}

# Handle convergence and variance updates for individual data (standard behavior)
handle_convergence_and_variance.individual <- function(data, model, model_prev, elbo_prev, elbo_current,
                                                       tol, estimate_residual_variance,
                                                       residual_variance_lowerbound, residual_variance_upperbound) {
  # Standard: Check convergence first, then update variance
  converged <- check_convergence(data, model_prev, model, elbo_prev, elbo_current, tol)

  if (!converged && estimate_residual_variance) {
    result <- update_model_variance(data, model, residual_variance_lowerbound, residual_variance_upperbound)
    data <- result$data
    model <- result$model
  }

  return(list(data = data, model = model, converged = converged))
}
