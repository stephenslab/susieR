#' @title Bayesian single-effect linear regression of Y on X
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
#' @param optimize_V boolean indicating whether to optimize V (by maximum likelihood)
#' @return a list with elements: \cr
#' \item{alpha}{vector of posterior inclusion probabilities. ie alpha[i] is posterior probability that
#'  that b[i] is non-zero}
#' \item{mu}{vector of posterior means (conditional on inclusion)}
#' \item{mu2}{vector of posterior second moments (conditional on inclusion)}
#' \item{lbf}{vector of log Bayes factors for each variable}
#' \item{V}{the prior variance (after optimization, if optimize_V is TRUE)}
#' \item{logBF}{(scalar) the loglikelihood for the total model minus the log-likelihood for the null model}
single_effect_regression_ss = function(Xty,dXtX,V=1,residual_variance=1,optimize_V=FALSE){
  betahat = (1/dXtX) * Xty
  shat2 = residual_variance/dXtX

  if(optimize_V){
    if(loglik.grad_ss(0,Xty,dXtX,residual_variance)<0){
      V=0
    } else {
      #V.o = optim(par=log(V),fn=negloglik.logscale,gr = negloglik.grad.logscale, X=X,Y=Y,s2=s2,method="BFGS")
      #if(V.o$convergence!=0){
      #  warning("optimization over prior variance failed to converge")
      #}
      V.u=uniroot(negloglik.grad.logscale_ss,c(-10,10),extendInt = "upX",Xty=Xty,dXtX=dXtX,s2=residual_variance)
      V = exp(V.u$root)
    }
  }

  lbf = dnorm(betahat,0,sqrt(V+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  #log(bf) on each SNP

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w is proportional to BF, but subtract max for numerical stability
  alpha = w/sum(w) # posterior prob on each SNP

  post_var = (1/V + dXtX/residual_variance)^(-1) # posterior variance
  post_mean = (dXtX/residual_variance) * post_var * betahat
  post_mean2 = post_var + post_mean^2 # second moment
  logBF = maxlbf + log(mean(w)) #analogue of loglik in the non-summary case

  return(list(alpha=alpha,mu=post_mean,mu2 = post_mean2,lbf=lbf, V=V, logBF = logBF))
}


loglik.grad_ss = function(V,Xty,dXtX,s2){
  betahat = (1/dXtX) * Xty
  shat2 = s2/dXtX

  lbf = dnorm(betahat,0,sqrt(V+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  #log(bf) on each SNP

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w =BF/BFmax
  alpha = w/sum(w)
  sum(alpha*lbf.grad(V,shat2,betahat^2/shat2))
}

# define loglikelihood and gradient as function of lV:=log(V)
# to improve numerical optimization
# negloglik.logscale_ss = function(lV, Xty,dXtX,s2){-loglik(exp(lV),Xty,dXtX,s2)}
negloglik.grad.logscale_ss = function(lV,Xty,dXtX,s2){-exp(lV)*loglik.grad_ss(exp(lV),Xty,dXtX,s2)}
