#' @title sets the attributes for the R matrix
#' @param R a p by p LD matrix
#' @param z a p vector of z scores
#' @return R with attribute e.g.
#'         attr(R, 'eigenR') is the eigen decomposition of R.

set_R_attributes = function(R, expected_dim) {
  if(nrow(R) != expected_dim) {
    stop(paste0('The dimension of R (', nrow(R), ' by ', nrow(R), ') does not agree with expected (', expected_dim, ' by ', expected_dim, ')'))
  }
  if(!is_symmetric_matrix(R)){
    stop('R is not a symmetric matrix.')
  }

  X0 = diag(R) == 0
  # convert any input R to correlation matrix
  # if R has 0 colums and rows, cov2cor produces NaN and warning
  R = muffled_cov2cor(R)
  # change the columns and rows with NaN to 0
  if(sum(X0) > 0){
    R[X0, ] = R[,X0] = 0
  }

  attr(R, 'eigenR') <- eigen(R, symmetric = TRUE)

  return(R)
}
