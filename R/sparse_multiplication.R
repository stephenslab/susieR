# @title Compute scaled.X %*% y using sparse multiplication
# @param X is a scaled dense matrix or an unscaled sparse matrix 
# @param y a p vector
# @param cm a p vector of column means 
# @param csd a p vector of column standard deviations
#
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
compute_sparse_Xy = function(X, y){
  cm = attr(X, 'scaled:center')
  csd = attr(X, 'scaled:scale')
  if(is.matrix(X)){
    return(X%*%y)
  }else{
    #scale Xy
    scaled.X  <- t(t(X)/csd)
    scaled.Xy <- tcrossprod(scaled.X,t(y))
    #center Xy
    Xy <- scaled.Xy - sum(cm*y/csd) 
    return(as.numeric(Xy))
  }
}

# @title Compute t(scaled.X)%*%y using sparse multiplication
# @param X is a scaled dense matrix or an unscaled sparse matrix 
# @param y an n vector
# @param cm a p vector of column means 
# @param csd a p vector of column standard deviations
#
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
compute_sparse_Xty = function(X, y){
  cm = attr(X, 'scaled:center')
  csd = attr(X, 'scaled:scale')
  if(is.matrix(X)){
    return(t(X)%*%y)
  }else{
    Xty        <- crossprod(X, y)
    scaled.Xty <- t(t(Xty)/csd)
    centered.scaled.Xty <- scaled.Xty - cm/csd * sum(y)     
    return(as.numeric(centered.scaled.Xty))
  }
}

# @title Compute M%*%t(scaled.X) using sparse multiplication
# @param M a L by p matrix
# @param X is a scaled dense matrix or an unscaled sparse matrix 
# @param cm a p vector of column means 
# @param csd a p vector of column standard deviations 
compute_sparse_MtX = function(M, X){
  cm = attr(X, 'scaled:center')
  csd = attr(X, 'scaled:scale')
  if(is.matrix(X)){
    return(M %*% t(X))
  }else{
    return(t(apply(M, 1, function(y) compute_sparse_Xy(X, y))))
  }
}
