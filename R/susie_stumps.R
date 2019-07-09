#' @title Applies susie to perform non-linear regression using "stumps" (trees
#' with  a  single split)
#'
#' @details TBD
#'
#' @param y an n vector of observations
#'
#' @param X an n by p matrix of covariates
#'
#' @param ... other parameters to pass to susie
#' @return a "susie" object; see ?susie for details
#' @examples
#' TBD
#' @export
susie_stumps = function(X, y, ...){
  n = length(y)
  attr(X, "matrix.type") = "stumps_matrix"
  s = susie(X=X, Y=y, ...)
  return(s)
}
