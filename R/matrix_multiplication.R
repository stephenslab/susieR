#' @title Computes standardized.X \%*\% b
#' @param X an n by p matrix with three attributes: attr(X, 'scaled:center'), attr(X, 'scaled:scale'), and attr(X, 'd')
#' @param b a p vector
#' @return an n vector
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
#' @keywords internal
compute_Xb = function(X, b){
  cm = attr(X, 'scaled:center')
  csd = attr(X, 'scaled:scale')
  #scale Xb
  #when X is a trend filtering matrix or stumps matrix, use special matrix mult
  if (is.tfmatrix(X))
    scaled.Xb <- compute_tf_Xb(attr(X, 'order'), b/csd)
  else if(is.stumps_matrix(X)){
    scaled.Xb <- compute_stumps_Xb(attr(X, 'Xord'),  b/csd)
  } else
  #when X is an ordinary sparse/dense matrix
     scaled.Xb <- tcrossprod(X, t(b/csd))
  #center Xb
  Xb <- scaled.Xb - sum(cm*b/csd)
  return(as.numeric(Xb))
}

#' @title Computes t(standardized.X)\%*\%y using sparse multiplication trick
#' @param X an n by p unstandardized matrix with three attributes: attr(X, 'scaled:center'), attr(X, 'scaled:scale'), and attr(X, 'd')
#' @param y an n vector
#' @return a p vector
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
compute_Xty = function(X, y){
  cm = attr(X, 'scaled:center')
  csd = attr(X, 'scaled:scale')
  ytX <- crossprod(y, X)
  #scale Xty
  #when X is a trend filtering matrix
  if (is.tfmatrix(X))
    scaled.Xty <- compute_tf_Xty(attr(X, 'order'),y)/csd
  else if(is.stumps_matrix(X))
    scaled.Xty <- compute_stumps_Xty(attr(X,'Xord'),y)/csd
  #when X is an ordinary sparse/dense matrix
  else
    scaled.Xty <- t(ytX/csd)
  #center Xty
  centered.scaled.Xty <- scaled.Xty - cm/csd * sum(y)
  return(as.numeric(centered.scaled.Xty))
}

#' @title Computes M\%*\%t(standardized.X) using sparse multiplication trick
#' @param M a L by p matrix
#' @param X an n by p unstandardized matrix with three attributes: attr(X, 'scaled:center'), attr(X, 'scaled:scale'), and attr(X, 'd')
#' @return a L by n matrix
#' @importFrom Matrix t
compute_MXt = function(M, X){
  cm = attr(X, 'scaled:center')
  csd = attr(X, 'scaled:scale')
  #when X is a trend filtering matrix
  if (is.tfmatrix(X) | is.stumps_matrix(X)) {
    return(as.matrix(t(apply(M,1,function(b) compute_Xb(X, b)))))
  }
  #when X is an ordinary sparse/dense matrix
  else return(as.matrix(tcrossprod(M,sweep(X,2,csd,"/")) - drop(tcrossprod(M, t(cm/csd)))))

    # This should be the same as
    #
    #   t(apply(M, 1, function(b) compute_Xb(X, b))))
    #
    # as well as
    #
    #   M %*% (t(X)/csd) - drop(tcrossprod(M,t(cm/csd)))
    #
    # but should be more memory-efficient.
}



