#' @title Compute sufficient statistics for susie.
#' 
#' @param X An n by p matrix of covariates.
#' 
#' @param y An n vector.
#' 
#' @param standardize Logical flag indicating whether to standardize
#'   columns of X to unit variance prior to fitting.
#' 
#' @return A list of sufficient statistics.
#' 
#' @importFrom methods as
#' 
#' @export
#' 
compute_ss = function(X, y, standardize = TRUE) {
  y = y - mean(y)
  is.sparse = !is.matrix(X)
  X = set_X_attributes(as.matrix(X),center=TRUE,scale = standardize)
  X = t((t(X) - attr(X,"scaled:center"))/attr(X,"scaled:scale"))
  XtX = crossprod(X)
  if(is.sparse)
    XtX = as(XtX,"dgCMatrix")
  Xty = c(y %*% X)
  n = length(y)
  yty = sum(y^2)
  return(list(XtX = XtX,Xty = Xty,yty = yty,n = n))
}
