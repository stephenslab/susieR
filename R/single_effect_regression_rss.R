#' @title Bayesian single-effect linear regression of Y on X.
#' @details Performs single-effect linear regression of Y on X. That is, this function
#' fits the regression model Y= Xb + e, where elements of e are iid N(0,residual_variance) and the
#' b is a p vector of effects to be estimated.
#' The assumption is that b has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=V).
#' Only the summary statistcs t(X)Y and diagonal elements of t(X)X are avialable.
#' @param z a p vector
#' @param Sigma residual_var * R + lambda I
#' @param V the prior variance
#' @param residual_variance the residual variance
#' @param prior_weights a p vector of prior weights
#' @param optimize_V boolean indicating whether to optimize V (by maximum likelihood)
#' @return a list with elements: \cr
#' \item{alpha}{vector of posterior inclusion probabilities. ie alpha[i] is posterior probability that
#'  that b[i] is non-zero}
#' \item{mu}{vector of posterior means (conditional on inclusion)}
#' \item{mu2}{vector of posterior second moments (conditional on inclusion)}
#' \item{lbf}{vector of log Bayes factors for each variable}
#' \item{V}{the prior variance (after optimization, if optimize_V is TRUE)}
#' \item{lbf_model}{(scalar) the loglikelihood for the total model minus the log-likelihood for the null model}
#'
#'
single_effect_regression_rss = function(z,Sigma,V=1,residual_variance=1,prior_weights=NULL,optimize_V=c("none", "optim", "EM")){
  p = length(z)
  shat2 = diag(Sigma)
  if (is.null(prior_weights))
    prior_weights = rep(1/p, p)

  if(optimize_V=="optim"){
    V=optimize_prior_variance_rss(optimize_V, z, Sigma, prior_weights, alpha=NULL, post_mean2=NULL)
  }

  lbf = sapply(1:p, function(j){
    -0.5 * log(1+V*attr(Sigma, 'RjSinvRj')[j]) +
      0.5 * (V/(1+V*attr(Sigma, 'RjSinvRj')[j])) * sum(attr(Sigma,'SinvRj')[[j]] * z)^2
  })
  #log(bf) on each SNP

  lbf[shat2==Inf] = 0 # deal with special case of infinite shat2 (eg happens if X does not vary)

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w is proportional to BF, but subtract max for numerical stability
  # posterior prob on each SNP
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w

  post_var = (attr(Sigma, 'RjSinvRj') + 1/V)^(-1) # posterior variance
  post_mean = sapply(1:p, function(j) (post_var[j]) * sum(attr(Sigma,'SinvRj')[[j]]* z))
  post_mean2 = post_var + post_mean^2 # second moment
  lbf_model = maxlbf + log(weighted_sum_w) #analogue of loglik in the non-summary case

  if(optimize_V=="EM"){
    V=optimize_prior_variance_rss(optimize_V, z, Sigma, prior_weights, alpha, post_mean2)
  }

  return(list(alpha=alpha,mu=post_mean,mu2 = post_mean2,lbf=lbf, V=V, lbf_model=lbf_model))
}

loglik_rss = function(V,z,Sigma,prior_weights) {
  p = length(z)

  lbf = sapply(1:p, function(j){
    -0.5 * log(1+V*attr(Sigma, 'RjSinvRj')[j]) +
      0.5 * (V/(1+V*attr(Sigma, 'RjSinvRj')[j])) * sum(attr(Sigma,'SinvRj')[[j]] * z)^2
  })
  #log(bf) on each SNP

  # lbf[shat2==Inf] = 0 # deal with special case of infinite shat2 (eg happens if X does not vary)

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w =BF/BFmax
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  return(log(weighted_sum_w)+ maxlbf)
}

neg.loglik_z.logscale_rss = function(lV,z,Sigma,prior_weights){
  return(-loglik_rss(exp(lV),z,Sigma,prior_weights))
}

optimize_prior_variance_rss = function(optimize_V, z, Sigma, prior_weights, alpha=NULL, post_mean2=NULL){
  if(optimize_V=="optim"){
    lV = optim(par=log(max(c(z^2, 1), na.rm = TRUE)),
               fn=neg.loglik_z.logscale_rss,
               z=z, Sigma = Sigma, prior_weights = prior_weights,
               method='Brent', lower = -10, upper = 15)$par
    V = exp(lV)
  }else if(optimize_V=="EM"){
    V = sum(alpha*post_mean2)
  }else stop('Invalid option for `optimize_V` method')
  if(loglik_rss(0,z,Sigma,prior_weights) >= loglik_rss(V,z,Sigma,prior_weights)) V=0 # set V exactly 0 if that beats the numerical value
  return(V)
}
