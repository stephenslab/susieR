#' @title Estimate residual variance
#' @param X an n by p matrix of covariantes
#' @param Y an n vector of data
#' @param s a susie fit
estimate_residual_variance = function(X,Y,s){
  n = nrow(X)
  return( (1/n)* get_ER2(X,Y,s) )
}


#' @title Estimate residual variance for summary statistics
#' @param XtX an n by p matrix of covariantes
#' @param Xty an n vector of data
#' @param s a susie fit
#' @param var_y variance of y
#' @param n sample size
estimate_residual_variance_ss = function(XtX,Xty,s,var_y,n){
  return( (1/n)* get_ER2_ss(XtX,Xty,s,var_y,n) )
}
