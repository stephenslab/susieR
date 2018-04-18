#this is for scaled prior in which effect prior variance is sa * sigma
# s is a susie fit
# a list with elements alpha, mu, sigma2, sa2
elbo_fn = function(X,Y,s){
  L = nrow(s$alpha)
  n = nrow(X)
  p = ncol(X)
  d = colSums(X*X) #note this is currently being computed 2 times - here and in Eloglik; could be made more efficient by avoiding this

  ss <- s$mu2 - s$mu^2 # posterior variance (conditional on inclusion)
  postb2 = s$alpha * s$mu2 # posterior second moment of b
  Ell = Eloglik(X,Y,s)

  sub = s$alpha>0 # this is to avoid taking 0 * log(0) in next line
  KL1 = sum(s$alpha[sub] * log(s$alpha[sub]/(1/p)))

  KL2 = - 0.5* rowSums(s$alpha * (1 + log(ss)-log(s$sigma2*s$sa2)))
  + 0.5 * rowSums(postb2/(s$sigma2*s$sa2))
  KL2[s$sa2==0] = 0 # deal with 0 prior

  return(Ell - KL1 - sum(KL2))
}

#' @title Get objective function from data and susie fit object.
#'
#' @param data A flash data object.
#'
#' @param f A flash fit object.
#'
#' @export
#'
susie_get_objective = function(X, Y, s) {
  return(Eloglik(X,Y,s)-sum(s$KL))
}



#' @title expected loglikelihood for a susie fit
Eloglik = function(X,Y,s){
  n = nrow(X)
  p = ncol(X)
  result =  -(n/2) * log(2*pi* s$sigma2) - (1/(2*s$sigma2)) * get_ER2(X,Y,s)
  return(result)
}


get_ER2 = function(X,Y,s){
  Xr = (s$alpha*s$mu) %*% t(X)
  Xrsum = colSums(Xr)

  d = colSums(X*X)
  postb2 = s$alpha * s$mu2 #posterior second moment

  return(sum((Y-Xrsum)^2) - sum(Xr^2) + sum(d*t(postb2)))
}


#' @title posterior expected loglikelihood for a single effect regression
#' @param X an n by p matrix of covariates
#' @param Y an n vector of regression outcome
#' @param s2 the residual variance
#' @param Eb the posterior mean of b (p vector) (alpha * mu)
#' @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
SER_posterior_e_loglik = function(X,Y,s2,Eb,Eb2){
  n = nrow(X)
  -0.5*n*log(2*pi*s2)  - (0.5/s2) * (sum(Y*Y) - 2*sum(Y*(X %*% Eb)) + sum(t(X^2)*as.vector(Eb2)))
}
