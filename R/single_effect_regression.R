#' @title Bayesian single-effect linear regression of Y on X
#' @details Performs single-effect linear regression of Y on X. That is, this function
#' fits the regression model Y= Xb + e, where elements of e are iid N(0,s2) and the
#' b is a p vector of effects to be estimated.
#' The assumption is that b has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=V).
#' @param Y an n vector
#' @param X an n by p matrix of covariates
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
#' \item{lbf_model}{log Bayes factor for the single effect regression}
#' \item{V}{the prior variance (after optimization, if optimize_V is TRUE)}
#' \item{loglik}{The log-likelihood p(Y|X,V)}
#'
#' @importFrom stats uniroot
#' @importFrom stats optim
#' @importFrom Matrix colSums
#'
single_effect_regression = function(Y,X,V,residual_variance=1,prior_weights=NULL, optimize_V=c("none", "optim", "EM")){
  optimize_V = match.arg(optimize_V)
  Xty = compute_Xty(X, Y)
  betahat = (1/attr(X, "d")) * Xty
  shat2 = residual_variance/attr(X, "d")
  if (is.null(prior_weights))
    prior_weights = rep(1/ncol(X), ncol(X))

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
  post_var = (1/V + attr(X, "d")/residual_variance)^(-1) # posterior variance
  post_mean = (1/residual_variance) * post_var * Xty
  post_mean2 = post_var + post_mean^2 # second moment
  # BF for single effect model
  lbf_model = maxlbf + log(weighted_sum_w)
  loglik = lbf_model + sum(dnorm(Y,0,sqrt(residual_variance),log=TRUE))

  if(optimize_V=="EM") V=optimize_prior_variance(optimize_V, betahat, shat2, prior_weights, alpha, post_mean2)

  return(list(alpha=alpha,mu=post_mean,mu2 = post_mean2,lbf=lbf,lbf_model=lbf_model,V=V,loglik=loglik))
}

#' @title estimate prior variance
#' @keywords internal
est_V_uniroot = function(betahat, shat2, prior_weights){
  V.u=uniroot(negloglik.grad.logscale,c(-10,10),extendInt = "upX",betahat=betahat,shat2=shat2,prior_weights=prior_weights)
  V = exp(V.u$root)
  return(V)
}

optimize_prior_variance = function(optimize_V, betahat, shat2, prior_weights, alpha=NULL, post_mean2=NULL){
  if(optimize_V=="optim"){
    lV = optim(par=log(max(c(betahat^2-shat2, 1), na.rm = TRUE)), fn=neg.loglik.logscale, betahat=betahat, shat2=shat2, prior_weights = prior_weights, method='Brent', lower = -10, upper = 15)$par
    V = exp(lV)
  }else if(optimize_V=="uniroot"){
    V = est_V_uniroot(betahat, shat2, prior_weights)
  }else if(optimize_V=="EM"){
    V = sum(alpha*post_mean2)
  }else stop('Invalid option for `optimize_V` method')
  if(loglik(0,betahat,shat2,prior_weights) >= loglik(V,betahat,shat2,prior_weights)) V=0 # set V exactly 0 if that beats the numerical value
  return(V)
}

# In these functions s2 represents residual_variance and shat2 is an estimate of it

#' The log likelihood function for SER model (based on summary data betahat, shat2), as a function of prior variance V
#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik = function(V,betahat,shat2,prior_weights) {

  lbf = dnorm(betahat,0,sqrt(V+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  #log(bf) on each SNP

  lbf[shat2==Inf] = 0 # deal with special case of infinite shat2 (eg happens if X does not vary)

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w =BF/BFmax
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  return(log(weighted_sum_w)+ maxlbf)
}

neg.loglik.logscale = function(lV,betahat,shat2,prior_weights){
  return(-loglik(exp(lV),betahat,shat2,prior_weights))
}

#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik.grad = function(V,betahat,shat2,prior_weights) {

  lbf = dnorm(betahat,0,sqrt(V+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  #log(bf) on each SNP

  lbf[shat2==Inf] = 0 # deal with special case of infinite shat2 (eg happens if X does not vary)

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w =BF/BFmax
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  sum(alpha*lbf.grad(V,shat2,betahat^2/shat2))
}

# define loglikelihood and gradient as function of lV:=log(V)
# to improve numerical optimization
#negloglik.logscale = function(lV, Y,X,s2){-loglik(exp(lV),Y,X,s2)}
negloglik.grad.logscale = function(lV,betahat,shat2,prior_weights){-exp(lV)*loglik.grad(exp(lV),betahat,shat2,prior_weights)}

#
# numDeriv::grad(negloglik.logscale,0, X =X, Y=Y,s2=s2)
# negloglik.grad.logscale(0,Y,X,s2)
#
# numDeriv::grad(loglik, 0.1, X =X, Y=Y,s2=s2)
# loglik.grad(0.1,X,Y,s2)
#
# numDeriv::grad(loglik, 1, X =X, Y=Y,s2=s2)
# loglik.grad(1,X,Y,s2)

# set.seed(1)
# n = 1000
# p = 1000
# beta = rep(0,p)
# beta[1] = 1
# X = matrix(rnorm(n*p),nrow=n,ncol=p)
# Y = X %*% beta + rnorm(n)
# s2 = 1
# optim(par=0,fn=negloglik.logscale,gr = negloglik.grad.logscale, X=X,Y=Y,s2=s2,method="BFGS")


# vector of gradients of logBF_j for each j, with respect to prior variance V
lbf.grad = function(V,shat2,T2){
  l = 0.5* (1/(V+shat2)) * ((shat2/(V+shat2))*T2-1)
  l[is.nan(l)] = 0
  return(l)
}

lbf = function(V,shat2,T2){
  l = 0.5*log(shat2/(V+shat2)) + 0.5*T2*(V/(V+shat2))
  l[is.nan(l)] = 0
  return(l)
}





