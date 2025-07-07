### SuSiE sufficient statistic data backend functions ###

# Initialize Fitted Values
initialize_fitted.ss <- function(data, alpha, mu){
  return(list(XtXr = data$XtX %*% colSums(alpha * mu)))
}

# Initialize Matrices
initialize_matrices.ss <- function(data, L, scaled_prior_variance, var_y,
                                   residual_variance, prior_weights, ...){
  p <- data$p

  mat_init <- list(
    alpha        = matrix(1 / p, L, p),
    mu           = matrix(0,     L, p),
    mu2          = matrix(0,     L, p),
    V            = rep(scaled_prior_variance * var_y, L),
    KL           = rep(as.numeric(NA), L),
    lbf          = rep(as.numeric(NA), L),
    lbf_variable = matrix(as.numeric(NA), L, p),
    sigma2       = residual_variance,
    pi           = prior_weights
  )

  return(mat_init)
}

# Get variance of y
get_var_y.ss <- function(data, ...) {
  return(data$yty / (data$n - 1))
}

# Extract core parameters of a susie fit across iterations
susie_extract_core.ss <- function(data, model, tracking, iter, track_fit, ...){
  if (isTRUE(track_fit)) {
    tracking[[iter]] <- list(alpha = model$alpha,
                             niter = iter,
                             V = model$V,
                             sigma2 = model$sigma2)
  }
  return(tracking)
}

# Validate Prior Variance
validate_prior.ss <- function(data, model, check_prior, ...) {
  if (isTRUE(check_prior)) {
    if (is.null(data$zm)) {
      d    <- attr(data$XtX, "d")
      bhat <- data$Xty / d
      shat <- sqrt(model$sigma2 / d)
      z    <- bhat / shat
      data$zm <- max(abs(z[!is.nan(z)]))
    }
    if (any(model$V > 100 * (data$zm^2))) {
      stop(
        "Estimated prior variance is unreasonably large.\n",
        "This usually caused by mismatch between the summary statistics and the LD Matrix.\n",
        "Please check the input."
      )
    }
  }
  invisible(TRUE)
}

# Posterior expected log-likelihood for a single effect regression
SER_posterior_e_loglik.ss = function (data, model, XtR, Eb, Eb2)
  return(-0.5/model$sigma2 * (-2*sum(Eb*XtR) + sum(attr(data$XtX,"d") * as.vector(Eb2))))

# Expected Squared Residuals
get_ER2.ss <- function (data, model) {
  B = model$alpha * model$mu
  XB2 = sum((B %*% data$XtX) * B)
  betabar = colSums(B)
  d = attr(data$XtX,"d")
  postb2 = model$alpha * model$mu2 # Posterior second moment.
  return(data$yty - 2*sum(betabar * data$Xty) + sum(betabar * (data$XtX %*% betabar)) -
           XB2 + sum(d * t(postb2)))
}

# Expected loglikelihood for a susie fit.
Eloglik.ss <- function (data, model) {
  return(-data$n/2*log(2*pi*model$sigma2) - 1/(2*model$sigma2) * get_ER2(data, model))
}

# Get Objective
get_objective.ss <- function (data, model) {
  objective <- Eloglik(data, model) - sum(model$KL)
  if(is.infinite(objective)){
    stop("get_objective.ss() produced an infinite ELBO value")
  }
  return(objective)
}

# Estimate Residual Variance
est_residual_variance.ss <- function(data, model){
  resid_var <- (1/data$n)*get_ER2(data,model)
  if(resid_var < 0){
    stop("est_residual_variance.ss() failed: the estimated value is negative")
  }
  return(resid_var)
}

# Single Effect Update
single_effect_update.ss <- function(
    data, model, l,
    optimize_V, check_null_threshold) {

  # Remove lth effect
  model$XtXr <- model$XtXr - data$XtX %*% (model$alpha[l, ] * model$mu[l, ])

  # Compute Residuals
  XtR <- data$Xty - model$XtXr
  d <- attr(data$XtX, "d")

  res <- single_effect_regression(
    Xty                  = XtR,
    dXtX                 = d,
    V                    = model$V[l],
    residual_variance    = model$sigma2,
    prior_weights        = model$pi,
    optimize_V           = optimize_V,
    check_null_threshold = check_null_threshold)

  res$KL <- -res$lbf_model +
    SER_posterior_e_loglik(data, model, XtR,
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

  model$XtXr <- model$XtXr + data$XtX %*% (model$alpha[l, ] * model$mu[l, ])

  return(model)
}

# Get column scale factors
get_scale_factors.ss <- function(data){
  return(attr(data$XtX,"scaled:scale"))
}

# Get intercept
get_intercept.ss <- function(data, model, ...){
  return(sum(data$X_colmeans * (colSums(model$alpha * model$mu)/model$X_column_scale_factors)))
}

# Get Fitted Values
get_fitted.ss <- function(data, model, ...){
  return(NULL)
}

# Get Credible Sets
get_cs.ss <- function(data, model, coverage, min_abs_corr, n_purity){

  if (is.null(coverage) || is.null(min_abs_corr)) return(NULL)

  if(any(!(diag(data$XtX) %in% c(0,1)))){
    return(susie_get_cs(model,coverage = coverage,
                        Xcorr = muffled_cov2cor(data$XtX),
                        min_abs_corr = min_abs_corr,
                        check_symmetric = FALSE,
                        n_purity = n_purity))
  }else{
    return(susie_get_cs(model,coverage = coverage,
                        Xcorr = data$XtX,
                        min_abs_corr = min_abs_corr,
                        check_symmetric = FALSE,
                        n_purity = n_purity))
  }
}

# Get PIP
get_pip.ss <- function(data, model, coverage, min_abs_corr, prior_tol){

  if (is.null(coverage) || is.null(min_abs_corr)) return(NULL)

  return(susie_get_pip(model,prune_by_cs = FALSE,prior_tol = prior_tol))

}

# Get Variable Names
get_variable_names.ss <- function(data, model, null_weight){
  if (!is.null(colnames(data$XtX))) {
    variable_names = colnames(data$XtX)
    if (!is.null(null_weight)) {
      variable_names[length(variable_names)] = "null"
      names(model$pip) = variable_names[-data$p]
    } else{
    names(model$pip)             = variable_names
    }
    colnames(model$alpha)        = variable_names
    colnames(model$mu)           = variable_names
    colnames(model$mu2)          = variable_names
    colnames(model$lbf_variable) = variable_names
  }
  return(model)
}

# Get univariate z-score
get_zscore.ss <- function(data, model, ...) {
  return(NULL)
}

# Configure ss data for specified method  
configure_data.ss <- function(data, non_sparse_method) {
  if (non_sparse_method == "none") {
    return(data)  # No changes needed
  } else {
    # Add non-sparse class and eigen decomposition
    class(data) <- c(paste0("ss_", non_sparse_method), "ss")
    data <- add_eigen_decomposition(data)
    return(data)
  }
}

# Add eigen decomposition to ss objects (for non-sparse methods)
add_eigen_decomposition.ss <- function(data) {
  # Compute eigen decomposition of correlation matrix
  eigen_decomp <- compute_eigen_decomposition(data$XtX, data$n)

  # Add eigen components to data object
  data$eigen_vectors <- eigen_decomp$V
  data$eigen_values  <- eigen_decomp$Dsq
  data$VtXty         <- t(eigen_decomp$V) %*% data$Xty  # Compute VtXty

  # Initialize derived quantities for non-sparse methods
  # These will be updated when variance components change
  sigmasq <- 1  # Default initial value
  tausq <- 0    # Default initial value
  var <- tausq * data$eigen_values + sigmasq
  data$var <- var
  data$diagXtOmegaX <- rowSums(sweep(data$eigen_vectors^2, 2, (data$eigen_values / var), `*`))
  data$XtOmegay <- data$eigen_vectors %*% (data$VtXty / var)

  return(data)
}

# Update variance components for ss data (standard approach)
update_variance_components.ss <- function(data, model) {
  sigma2 <- est_residual_variance(data, model)
  return(list(sigma2 = sigma2, tausq = NULL))
}

# Update derived quantities for ss data (no-op for standard)
update_derived_quantities.ss <- function(data, model) {
  return(data)  # No changes needed for standard ss data
}

# Check convergence for ss data (uses ELBO)
check_convergence.ss <- function(data, model_prev, model_current, elbo_prev, elbo_current, tol) {
  # Standard ELBO-based convergence (uses pre-computed ELBO values)
  return(elbo_current - elbo_prev < tol)
}

# Update variance before convergence check for ss data
update_variance_before_convergence.ss <- function(data) {
  # Standard behavior: Check convergence first, then update variance
  return(FALSE)
}

# Handle convergence and variance updates for ss data (standard behavior)
handle_convergence_and_variance.ss <- function(data, model, model_prev, elbo_prev, elbo_current, 
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

