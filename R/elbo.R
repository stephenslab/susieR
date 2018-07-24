#' @title Get objective function from data and susie fit object.
#'
#' @param data A flash data object.
#'
#' @param f A flash fit object.
#'
#' @export
#'
susie_get_objective = function(X, X.sparse, Y, cm, csd, s){
  return(Eloglik(X, X.sparse, Y, cm, csd, s)-sum(s$KL))
}

#' @title expected loglikelihood for a susie fit
Eloglik = function(X, X.sparse, Y, cm, csd, s){
  n = nrow(X)
  p = ncol(X)
  result =  -(n/2) * log(2*pi* s$sigma2) - (1/(2*s$sigma2)) * get_ER2(X, X.sparse, Y, cm, csd, s)
  return(result)
}

#' @title expected squared residuals
get_ER2 = function(X, X.sparse, Y, cm, csd, s){
  M = s$alpha*s$mu
  Xr = compute_sparse_MtX(M, X.sparse, cm, csd)
  Xrsum = colSums(Xr)

  d = colSums(X*X)
  postb2 = s$alpha * s$mu2 #posterior second moment

  return(sum((Y-Xrsum)^2) - sum(Xr^2) + sum(d*t(postb2)))
}

#' @title posterior expected loglikelihood for a single effect regression
#' @param X an n by p matrix of covariates, scaled 
#' @param X.sparse an n by p matrix of covariates, unscaled
#' @param Y an n vector of regression outcome
#' @param cm a p vector of column means
#' @param csd a p vector of column standard deviations
#' @param s2 the residual variance
#' @param Eb the posterior mean of b (p vector) (alpha * mu)
#' @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
SER_posterior_e_loglik = function(X, X.sparse,Y,cm,csd, s2,Eb,Eb2){
  n = nrow(X)
  XEb = compute_sparse_Xy(X.sparse, Eb, cm, csd)
  -0.5*n*log(2*pi*s2)  - (0.5/s2) * (sum(Y*Y) - 2*sum(Y*XEb) + sum(t(X^2)*as.vector(Eb2)))
}
