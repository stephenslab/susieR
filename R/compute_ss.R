#' @title Compute the needed summary statistics for `susie_ss` from the 'raw' data
#' @param X an n by p matrix of covariates
#' @param y an n vector
#' @param standardize logical flag (default=TRUE) for whether to standardize columns of X to unit variance prior to fitting.
#' @return a list with the input for `susie_ss`
#' @export
compute_ss = function(X, y, standardize = TRUE){
  y = y - mean(y)
  is.sparse = !(is.matrix(X))
  X = susieR:::safe_colScale(as.matrix(X), center=TRUE, scale = standardize)

  XtX = crossprod(X)
  if(is.sparse){
    XtX = as(XtX, 'dgCMatrix')
  }
  Xty = c(y %*% X)
  n = length(y)
  vary = sum(y^2)/n

  return(list(XtX = XtX, Xty = Xty, vary = vary, n = n))
}
