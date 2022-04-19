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

  # Drop small eigenvalues.
  eigenR$values[abs(eigenR$values) < r_tol] = 0
  if(any(eigenR$values < 0)) {
    min_lambda = min(eigenR$values)
    eigenR$values[eigenR$values < 0] = 0
    warning_message(paste0("The input correlation matrix has negative eigenvalues ",
                   "(smallest one is ", min_lambda, "). The correlation ",
                   "matrix is adjusted such that these negative eigenvalues ",
                   "are now zeros. You can ignore this message, only if you ",
                   "believe the negative eigenvalue is result of numerical ",
                   "rounding errors."))
  }
  res = eigenR$vectors %*% (t(eigenR$vectors) * eigenR$values)

  attr(res,"eigen") = eigenR
  attr(res,"d") = diag(res)
  attr(res,"scaled:scale") = rep(1,length = nrow(R))
  return(res)
}
