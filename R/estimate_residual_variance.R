#' @title Estimate residual variance
#' @param X an n by p matrix of covariantes
#' @param Y an n vector of data
#' @param s a susie fit
estimate_residual_variance = function(X,Y, s){
  n = nrow(X)
  return( (1/n)* get_ER2(X,Y,s) )
}


