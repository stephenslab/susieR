### SuSiE sufficient statistic data backend functions ###

# Initialize Fitted Values
initialize_fitted.ss <- function(data, alpha, mu){
  return(list(XtXr = data$XtX %*% colSums(alpha * mu)))
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

  if(any(!(diag(XtX) %in% c(0,1)))){
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
