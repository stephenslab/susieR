#' @title Bayesian single-effect linear regression of Y on X.
#' @details Performs single-effect linear regression of Y on X. That is, this function
#' fits the regression model Y= Xb + e, where elements of e are iid N(0,residual_variance) and the
#' b is a p vector of effects to be estimated.
#' The assumption is that b has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=V).
#' Only the summary statistcs t(X)Y and diagonal elements of t(X)X are avialable.
#' @param Xty a p vector
#' @param dXtX a p vector, diagonal elements of t(X)X
#' @param V the prior variance
#' @param residual_variance the residual variance
#' @param prior_weights a p vector of prior weights
#' @param optimize_V indicating the method to optimize V (by maximum likelihood)
#' @return a list with elements: \cr
#' \item{alpha}{vector of posterior inclusion probabilities. ie alpha[i] is posterior probability that
#'  that b[i] is non-zero}
#' \item{mu}{vector of posterior means (conditional on inclusion)}
#' \item{mu2}{vector of posterior second moments (conditional on inclusion)}
#' \item{lbf}{vector of log Bayes factors for each variable}
#' \item{V}{the prior variance (after optimization, if optimize_V is TRUE)}
#' \item{lbf_model}{(scalar) the loglikelihood for the total model minus the log-likelihood for the null model}
#'
#' @importFrom stats uniroot
#' @importFrom stats optim
#'
single_effect_regression_ss = function(Xty,dXtX,V=1,residual_variance=1,prior_weights=NULL,optimize_V=c("none", "optim", "EM")){
  optimize_V = match.arg(optimize_V)
  betahat = (1/dXtX) * Xty
  shat2 = residual_variance/dXtX
  if (is.null(prior_weights))
    prior_weights = rep(1/length(dXtX), length(dXtX))

  if(optimize_V=="optim") V=optimize_prior_variance(optimize_V, betahat, shat2, prior_weights, alpha=NULL, post_mean2=NULL)

  lbf = dnorm(betahat,0,sqrt(V+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  #log(bf) on each SNP

  lbf[shat2==Inf] = 0 # deal with special case of infinite shat2 (eg happens if X does not vary)

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w is proportional to BF, but subtract max for numerical stability
  # posterior prob on each SNP
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w

  post_var = (1/V + dXtX/residual_variance)^(-1) # posterior variance
  post_mean = (1/residual_variance) * post_var * Xty
  post_mean2 = post_var + post_mean^2 # second moment
  lbf_model = maxlbf + log(weighted_sum_w) #analogue of loglik in the non-summary case

  if(optimize_V=="EM") V=optimize_prior_variance(optimize_V, betahat, shat2, prior_weights, alpha, post_mean2)

  return(list(alpha=alpha,mu=post_mean,mu2 = post_mean2,lbf=lbf, V=V, lbf_model=lbf_model))
}


