#' @title Bayesian single-effect linear regression of Y on X
#' @details Performs single-effect linear regression of Y on X. That is, this function
#' fits the regression model Y= Xb + e, where elements of e are iid N(0,s2) and the
#' b is a p vector of effects to be estimated.
#' The assumption is that b has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=sa2*s2).
#' @param Y an n vector
#' @param X an n by p matrix of covariates
#' @param sa2 the scaled prior variance (so prior variance is sa2*s2)
#' @param s2 the residual variance
#' @return a list with elements: \cr
#' \item{alpha}{vector of posterior inclusion probabilities. ie alpha[i] is posterior probability that
#'  that b[i] is non-zero}
#' \item{mu}{vector of posterior means (conditional on inclusion)}
#' \item{mu2}{vector of posterior second moments (conditional on inclusion)}
#' \item{bf}{vector of Bayes factors for each variable}
single_effect_regression = function(Y,X,sa2=1,s2=1){
  d = colSums(X^2)
  sa2 = s2*sa2 # scale by residual variance
  betahat = (1/d) * t(X) %*% Y
  shat2 = s2/d

  lbf = dnorm(betahat,0,sqrt(sa2+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  #log(bf) on each SNP

  maxlbf = max(lbf)
  w = exp(lbf-max(lbf)) # w is proportional to BF, but subtract max for numerical stability
  alpha = w/sum(w) # posterior prob on each SNP

  post_var = (1/sa2 + d/s2)^(-1) # posterior variance
  post_mean = (d/s2) * post_var * betahat
  post_mean2 = post_var + post_mean^2 # second moment

  return(list(alpha=alpha,mu=post_mean,mu2 = post_mean2,lbf=lbf))
}
