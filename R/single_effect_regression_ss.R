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
#' @param optimize_V boolean indicating whether to optimize V (by maximum likelihood)
#' @param optimV_method the method to estimate V, 'optim', 'EM' or 'uniroot'
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
#'
single_effect_regression_ss = function(Xty,dXtX,V=1,residual_variance=1,prior_weights=NULL,optimize_V=FALSE, optimV_method = "optim", niter){
  betahat = (1/dXtX) * Xty
  shat2 = residual_variance/dXtX
  if (is.null(prior_weights))
    prior_weights = rep(1/length(dXtX), length(dXtX))

  if(optimize_V && optimV_method == "uniroot"){
    V = est_V_uniroot(betahat, shat2, prior_weights)
    if(loglik(0,betahat,shat2,prior_weights) >= loglik(V,betahat,shat2,prior_weights)){
      V=0 # set V exactly 0 if that beats the numerical value
    }
  }

  if(optimize_V && optimV_method=="optim"){
    lV = optim(par=log(max(c(betahat^2-shat2, 1), na.rm = TRUE)), fn=neg.loglik.logscale,
               gr = negloglik.grad.logscale, betahat=betahat, shat2=shat2, prior_weights = prior_weights,
               method='Brent', lower = -10, upper = 15)$par
    V = exp(lV)
    if(loglik(0,betahat,shat2,prior_weights) >= loglik(V,betahat,shat2,prior_weights)){
      V=0 # set V exactly 0 if that beats the numerical value
    }
  }


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
  if(optimize_V && optimV_method == "EM"){
    # if(niter <= 5 || niter %% 3)
    V = sum(alpha*post_mean2)
    if(loglik(0,betahat,shat2,prior_weights) >= loglik(V,betahat,shat2,prior_weights)){
      V=0 # set V exactly 0 if that beats the numerical value
    }
  }
  return(list(alpha=alpha,mu=post_mean,mu2 = post_mean2,lbf=lbf, V=V, lbf_model=lbf_model))
}
