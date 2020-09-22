# @title sets three attributes for matrix X
# @param X an n by p data matrix that can be either a trend filtering
#   matrix or a regular dense/sparse matrix
# @param center boolean indicating centered by column means or not
# @param scale boolean indicating scaled by column standard deviations or not
# @return X with three attributes e.g. `attr(X, 'scaled:center') is a
# p vector of column means of X if center=TRUE, a p vector of zeros
# otherwise. 'attr(X, 'scaled:scale') is a p vector of columan standard
# deviations of X if scale=TRUE, a p vector of 1s otherwise. 'attr(X,
# 'd') is a p vector of column sums of X.standardized^2,' where
# X.standardized is the matrix X centered by attr(X, 'scaled:center')
# and scaled by attr(X, 'scaled:scale').
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
set_X_attributes = function(X, center = TRUE, scale = TRUE) {
    
  # if X is a trend filtering matrix
  if (!is.null(attr(X,"matrix.type"))) {
    order = attr(X,"order")
    n = ncol(X)
    
    # Set three attributes for X.
    attr(X,"scaled:center") = compute_tf_cm(order,n)
    attr(X,"scaled:scale") = compute_tf_csd(order,n)
    attr(X,"d") = compute_tf_d(order,n,attr(X,"scaled:center"),
                               attr(X,"scaled:scale"),scale,center)
    if (!center)
      attr(X,"scaled:center") = rep(0,n)
    if (!scale)
      attr(X,"scaled:scale") = rep(1,n)
  } else {
      
    # If X is either a dense or sparse ordinary matrix.
    # Get column means.
    cm = colMeans(X,na.rm = TRUE)
    
    # Get column standard deviations.
    csd = compute_colSds(X)
    
    # Set sd = 1 when the column has variance 0.
    csd[csd == 0] = 1
    if (!center)
      cm = rep(0,length = length(cm))
    if (!scale) 
      csd = rep(1,length = length(cm))
    X.std = (t(X) - cm)/csd
    
    # Set three attributes for X.
    attr(X,"d") = rowSums(X.std * X.std)
    attr(X,"scaled:center") = cm
    attr(X,"scaled:scale") = csd
  }
  return(X)
}

# @title computes column standard deviations for any type of matrix
# @details This should give the same result as matrixStats::colSds(X),
#   but allows for sparse matrices as well as dense ones.
# @param X an n by p matrix of any type, e.g. sparse, dense.
# @return a p vector of column standard deviations.
#
#' @importFrom Matrix colSums
compute_colSds = function(X) {
  n = nrow(X)
  return(sqrt((colSums(X^2)/n - (colSums(X)/n)^2)*(n/(n-1))))
}
