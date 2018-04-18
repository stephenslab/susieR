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

# expected squared residuals
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
