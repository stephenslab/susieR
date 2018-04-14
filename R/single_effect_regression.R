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
single_effect_regression = function(Y,X,sa2=1,s2=1,optimize_sa2=FALSE){
  d = colSums(X^2)
  V = s2*sa2 # scale by residual variance
  betahat = (1/d) * t(X) %*% Y
  shat2 = s2/d

  if(optimize_sa2){
    if(loglik.grad(0,Y,X,s2)<0){
      V=0
    } else {
      V.o = optim(par=V,fn=negloglik.logscale,gr = negloglik.grad.logscale, X=X,Y=Y,s2=s2,method="BFGS")
      V = exp(V.o$par)
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

  return(list(alpha=alpha,mu=post_mean,mu2 = post_mean2,lbf=lbf,sa2=V/s2))
}

# compute loglik (up to a constant that depends on s2 but not on sa2)
loglik = function(V,Y,X,s2){
  d = colSums(X^2)
  betahat = (1/d) * t(X) %*% Y
  shat2 = s2/d

  lbf = dnorm(betahat,0,sqrt(V+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  #log(bf) on each SNP

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w =BF/BFmax
  return(maxlbf + log(mean(w)))
}

loglik.grad = function(V,Y,X,s2){
  d = colSums(X^2)
  betahat = (1/d) * t(X) %*% Y
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
negloglik.logscale = function(lV, Y,X,s2){-loglik(exp(lV),Y,X,s2)}
negloglik.grad.logscale = function(lV,Y,X,s2){-exp(lV)*loglik.grad(exp(lV),Y,X,s2)}

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
  0.5* (1/(V+shat2)) * ((shat2/(V+shat2))*T2-1)
}

lbf = function(V,shat2,T2){
  0.5*log(shat2/(V+shat2)) + 0.5*T2*(V/(V+shat2))
}







#
# em_update_prior_variance_single_regression = function(Y,X,sa2,s2){
#   d = colSums(X^2)
#   sa2 = s2*sa2 # scale by residual variance
#   betahat = (1/d) * t(X) %*% Y
#   shat2 = s2/d
#
#   lbf = dnorm(betahat,0,sqrt(sa2+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
#   #log(bf) on each SNP
#
#   maxlbf = max(lbf)
#   w = exp(lbf-max(lbf)) # w is proportional to BF, but subtract max for numerical stability
#   alpha = w/sum(w) # posterior prob on each SNP
#
#   t2 = betahat^2/shat2 # t-squared
#   # if(min(d)==max(d)){
#   #       vmax = shat2[1]*(sum(alpha*t2)-1)
#   # }
#   #
#   # this is the derivative of the complete data log-likelihood wrt v = prior variance
#   cdll_negloglik = function(v){
#     lbf = dnorm(betahat,0,sqrt(v+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
#     return(-sum(alpha*lbf))
#   }
#
#   cdll_negloglik.logscale = function(v){
#     lbf = dnorm(betahat,0,sqrt(exp(v)+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
#     return(-sum(alpha*lbf))
#   }
#
#
#   cdll_negderiv = function(v){-0.5*sum(alpha * (1/(v+shat2)) * ((shat2/(shat2+v))*t2 -1))}
#   v_upper = max(shat2*(t2-1)) # upper bound on vhat
#
#
#
#
#   if(v_upper>0 && cdll_deriv(0)>0){
#     v_init = v_upper/2
#     v_opt = optim(log(vinit), cdll_negloglik.logscale, method="BFGS")
#   #  v_opt = uniroot(cdll_deriv,interval = c(0,v_upper), extendInt = "downX") #$root
#   } else{}
#     v_opt = 0
#   }
# }
