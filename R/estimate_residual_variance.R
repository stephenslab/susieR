# @title Estimate residual variance
# @param X an n by p matrix of covariantes
# @param y an n vector of data
# @param s a susie fit
estimate_residual_variance = function (X, y, s) {
  n = nrow(X)
  return((1/n)*get_ER2(X,y,s))
}

# @title Estimate residual variance for summary statistics
# @param XtX a p by p matrix
# @param Xty a p vector
# @param s a susie fit
# @param yty a scaler, y'y, where y is centered to have mean 0
# @param n sample size
estimate_residual_variance_ss = function (XtX, Xty, s, yty, n)
  (1/n)*get_ER2_ss(XtX,Xty,s,yty)
