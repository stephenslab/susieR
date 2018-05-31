#' @title Bayesian single-effect linear regression of Y on X
#' @details Performs single-effect linear regression of Y on X. That is, this function
#' fits the regression model Y= Xb + e, where elements of e are iid N(0,s2) and the
#' b is a p vector of effects to be estimated.
#' The assumption is that b has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=sa2*s2).
#' Only the summary statistcs t(X)Y and diagonal elements of t(X)X are avialable.
#' @param Xty an p vector
#' @param dXtX an p vector, diagonal elements of t(X)X
#' @param sa2 the scaled prior variance (so prior variance is sa2*s2)
#' @param s2 the residual variance
#' @return a list with elements: \cr
#' \item{alpha}{vector of posterior inclusion probabilities. ie alpha[i] is posterior probability that
#'  that b[i] is non-zero}
#' \item{mu}{vector of posterior means (conditional on inclusion)}
#' \item{mu2}{vector of posterior second moments (conditional on inclusion)}
#' \item{bf}{vector of Bayes factors for each variable}
single_effect_regression_ss = function(Xty,dXtX,sa2=1,s2=1,optimize_sa2=FALSE){
  d = dXtX
  V = s2*sa2 # scale by residual variance
  betahat = (1/d) * Xty
  shat2 = s2/d

  if(optimize_sa2){
    if(loglik.grad_ss(0,Xty,dXtX,s2)<0){
      V=0
    } else {
      #V.o = optim(par=log(V),fn=negloglik.logscale,gr = negloglik.grad.logscale, X=X,Y=Y,s2=s2,method="BFGS")
      #if(V.o$convergence!=0){
      #  warning("optimization over prior variance failed to converge")
      #}
      V.u=uniroot(negloglik.grad.logscale_ss,c(-10,10),extendInt = "upX",Xty=Xty,dXtX=d,s2=s2)
      V = exp(V.u$root)
    }
  }

  lbf = dnorm(betahat,0,sqrt(V+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  #log(bf) on each SNP

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w is proportional to BF, but subtract max for numerical stability
  alpha = w/sum(w) # posterior prob on each SNP

  post_var = (1/V + d/s2)^(-1) # posterior variance
  post_mean = (d/s2) * post_var * betahat
  post_mean2 = post_var + post_mean^2 # second moment
  # loglik = maxlbf + log(mean(w)) + sum(dnorm(Y,0,sqrt(s2),log=TRUE))

  return(list(alpha=alpha,mu=post_mean,mu2 = post_mean2,lbf=lbf,sa2=V/s2))
}


loglik.grad_ss = function(V,Xty,dXtX,s2){
  d = dXtX
  betahat = (1/d) * Xty
  shat2 = s2/d

  lbf = dnorm(betahat,0,sqrt(V+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  #log(bf) on each SNP

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w =BF/BFmax
  alpha = w/sum(w)
  sum(alpha*lbf.grad(V,shat2,betahat^2/shat2))
}

# define loglikelihood and gradient as function of lV:=log(V)
# to improve numerical optimization
# negloglik.logscale_ss = function(lV, XtY,dXtX,s2){-loglik(exp(lV),XtY,dXtX,s2)}
negloglik.grad.logscale_ss = function(lV,Xty,dXtX,s2){-exp(lV)*loglik.grad_ss(exp(lV),Xty,dXtX,s2)}
