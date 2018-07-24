#' @title Estimate residual variance
#' @param X an n by p matrix of covariates, scaled 
#' @param X.sparse an n by p matrix of covariates, unscaled
#' @param Y an n vector of regression outcome
#' @param cm a p vector of column means
#' @param csd a p vector of column standard deviations
#' @param s a susie fit
estimate_residual_variance = function(X, X.sparse, Y, cm, csd, s){
  n = nrow(X)
  return( (1/n)* get_ER2(X, X.sparse, Y, cm, csd, s) )
}


