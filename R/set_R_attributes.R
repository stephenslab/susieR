#' @title sets the attributes for the X matrix
#' @param R a p by p LD matrix
#' @param expected_dim the expected dimension for R
#' @param r_tol tolerance level for eigen value check of positive semidefinite matrix of R.
#' @param z a p vector of z scores
#' @return R with two attributes e.g.
#'         attr(R, 'det') is the determinant of R. It is 1 if R is not full rank.
#'         attr(R, 'ztRinvz') is t(z)R^{-1}z. We use pseudoinverse of R when R is not invertible.

set_R_attributes = function(R, expected_dim, r_tol, z) {
  svdR <- svd(R)
  eigenvalues <- svdR$d
  eigenvalues[abs(eigenvalues) < r_tol] <- 0

  if(all(eigenvalues > 0)){
    attr(R, 'det') = prod(eigenvalues)
  }else{
    attr(R, 'det') = 1
  }

  if(!missing(z)){
    Dinv = numeric(expected_dim)
    Dinv[eigenvalues != 0] = 1/(eigenvalues[eigenvalues!=0])
    attr(R, 'ztRinvz') <- sum(z*(svdR$v %*% (Dinv * crossprod(svdR$u, z))))
  }

  return(R)
}
