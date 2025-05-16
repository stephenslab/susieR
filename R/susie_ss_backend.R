### SuSiE sufficient statistic data backend functions ###

# Posterior expected log-likelihood for a single effect regression
SER_posterior_e_loglik.ss = function (data, model, XtR, Eb, Eb2)
  return(-0.5/model$sigma2 * (-2*sum(Eb*XtR) + sum(attr(data$XtX,"d") * as.vector(Eb2))))

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
