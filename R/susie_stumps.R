#' @title Applies susie to perform non-linear regression using "stumps" (trees
#' with  a  single split)
#'
#' @details TBD
#'
#' @param X an n by p matrix of covariates
#'
#' @param Y an n vector of observations
#'
#' @param include_linear logical indicating whether to include the linear terms of X as well as the non-parametric stumps
#'
#' @param ... other parameters to pass to susie
#' @return a "susie" object; see ?susie for details
#' @examples
#' TBD
#' @export
susie_stumps = function(X, Y, include_linear = TRUE, ...){
  xl = make_stumps_matrix(X, include_linear)
  s = susie(X=xl, Y=Y, standardize=FALSE,...)
  return(s)
}

#' if Xtrain is supplied then use that to define the breaks of
#' the stumps
#' @export
make_stumps_matrix=function(X, include_linear, Xtrain=NULL){
  if(is.null(Xtrain)){Xtrain = X}

  xl=list() # initialize
  if(include_linear){ #include X as a regular matrix first
    attr(X,"nrow") <- nrow(X)
    attr(X,"ncol") <- ncol(X)
    attr(X,"scaled:center") <- rep(0,ncol(X))
    attr(X,"scaled:scale") <- rep(1,ncol(X))
    xl=c(xl,list(X))
  }

  for(i in 1:ncol(X)){xl= c(xl,list(susieR:::make_tfg_matrix(X[,i],Xtrain[,i])))}
  return(xl)
}

#' @title Applies susie to perform non-linear regression using "stumps" (trees
#' with  a  single split). Old version (original beta implementation)
#'
#' @details TBD
#'
#' @param X an n by p matrix of covariates
#'
#' @param Y an n vector of observations
#'
#' @param ... other parameters to pass to susie
#' @return a "susie" object; see ?susie for details
#' @examples
#' TBD
#' @export
susie_stumps_old = function(X, Y, ...){
  attr(X, "matrix.type") = "stumps_matrix"
  s = susie(X=X, Y=Y, ...)
  return(s)
}
