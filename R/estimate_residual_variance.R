#' @title Estimate residual variance
#' @param X an n by p matrix of covariantes
#' @param Y an n vector of data
#' @param s a susie fit
estimate_residual_variance = function(X,Y, s){
  # R = get_R(X,Y,s) #residuals
  # d = colSums(X^2)
  # post_var = s$alpha*s$mu2 - (s$alpha*s$mu)^2
  # V = colSums(post_var)
  # n = nrow(X)
  return( (1/n)* get_ER2(X,Y,s) )
}


