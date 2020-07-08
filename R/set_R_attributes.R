# @title sets the attributes for the R matrix
# @param R a p by p LD matrix
# @param r_tol tolerance level for eigen value check of positive
#   semidefinite matrix of R.
# @return R with attribute e.g., attr(R, 'eigenR') is the eigen
#   decomposition of R.
set_R_attributes = function (R, r_tol) {
  if (is.null(attr(R,"eigen")))
    eigenR = eigen(R,symmetric = TRUE)
  else
    eigenR = attr(R,"eigen")

  # drop small eigenvalues
  eigenR$values[abs(eigenR$values) < r_tol] = 0
  if(any(eigenR$values < 0)) {
    eigenR$values[eigenR$values < 0] = 0
    warning("Negative eigenvalues are set to 0")
  }
  res = eigenR$vectors %*% (t(eigenR$vectors) * eigenR$values)

  attr(res,"eigen") = eigenR
  attr(res,"d") = diag(res)
  attr(res,"scaled:scale") = rep(1,length = nrow(R))
  return(res)
}
