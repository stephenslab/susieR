# @title computes column standard deviations for any type of matrix
# @details This should give the same result as matrixStats::colSds(X),
#   but allows for sparse matrices as well as dense ones.
# @param X an n by p matrix of any type, e.g. sparse, dense.
# @return a p vector of column standard deviations.
#
#' @importFrom Matrix colSums
#' @importFrom matrixStats colSds
compute_colSds = function(X) {
  if (is.matrix(X))
    y = colSds(X)
  else {
    n = nrow(X)
    y = sqrt((colSums(X^2)/n - (colSums(X)/n)^2)*(n/(n-1)))
  }
  return(y)
}
