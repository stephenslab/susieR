### SuSiE Individual-level data backend functions ###

# Initialize Fitted Values
initialize_fitted.individual <- function(data, alpha, mu){
  return(list(Xr = compute_Xb(data$X, colSums(alpha * mu))))
}

# Get variance of y
get_var_y.individual <- function(data, ...) {
  return(var(drop(data$y)))
}

# Extract core parameters of a susie fit across iterations
susie_extract_core.individual <- function(data, model, tracking, iter, track_fit, ...){
  if (isTRUE(track_fit)) {
    tracking[[iter]] <- list(alpha = model$alpha,
                             niter = iter,
                             V = model$V,
                             sigma2 = model$sigma2)
  }
  return(tracking)
}

# Validate Prior Variance
validate_prior.individual <- function(data, model, check_prior, ...) {
  invisible(TRUE)
}

# Posterior expected log-likelihood for a single effect regression
SER_posterior_e_loglik.individual <- function (data, model, R, Eb, Eb2) {
  return(-0.5*data$n*log(2*pi*model$sigma2) - 0.5/model$sigma2*(sum(R*R)
                                                                - 2*sum(R*compute_Xb(data$X,Eb))
                                                                + sum(attr(data$X,"d") * Eb2)))}

# Expected Squared Residuals
get_ER2.individual <- function (data, model) {
  Xr_L = compute_MXt(model$alpha * model$mu, data$X) # L by N matrix
  postb2 = model$alpha * model$mu2 # Posterior second moment.
  return(sum((data$y - model$Xr)^2) - sum(Xr_L^2) + sum(attr(data$X,"d") * t(postb2)))
}


# Expected loglikelihood for a susie fit.
Eloglik.individual <- function (data, model) {
  return(-(data$n/2) * log(2*pi*model$sigma2) - (1/(2*model$sigma2)) * get_ER2(data, model))
}

# Get Objective
get_objective.individual <- function (data, model) {
  objective <- Eloglik(data,model) - sum(model$KL)
  if(is.infinite(objective)){
    stop("get_objective.individual() produced an infinite ELBO value")
  }
  return(objective)
}

# Estimate Residual Variance
est_residual_variance.individual <- function(data, model){
  resid_var <- (1/data$n)*get_ER2(data,model)
  if(resid_var < 0){
    stop("est_residual_variance.individual() failed: the estimated value is negative")
  }
  return(resid_var) #TODO: Instead of stopping, we could give warning() and run MoM estimator instead.
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
    return(data$mean_y - sum(attr(data$X,"scaled:center") *
                                 (colSums(model$alpha * model$mu)/attr(data$X,"scaled:scale"))))
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

# Get PIP
get_pip.individual <- function(data, model, coverage, min_abs_corr, prior_tol){

  if (is.null(coverage) || is.null(min_abs_corr)) return(NULL)

  return(susie_get_pip(model,prune_by_cs = FALSE,prior_tol = prior_tol))

}

# Get Variable Names
get_variable_names.individual <- function(data, model, null_weight){
  if (!is.null(colnames(data$X))) {
    variable_names = colnames(data$X)
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
