#' @title problem1: compute scaled.X %*% y using sparse multiplication
#' @param X.sparse an nonscaled n by p sparse matrix
#' @param y a p vector
#' @param cm a p vector of column means 
#' @param csd a p vector of column standard deviations 
compute_sparse_Xy = function(X.sparse, y, cm, csd){
  if(is.matrix(X.sparse)){
    return(X.sparse%*%y)
  }else{
    #scale Xy
    scaled.X = t(t(X.sparse)/csd)
    scaled.Xy = tcrossprod(scaled.X, t(y))
    #center Xy
    Xy = scaled.Xy - sum(cm*y/csd) 
    return(as.numeric(Xy))
  }
}

#' @title problem2: compute t(scaled.X)%*%y using sparse multiplication
#' @param X.sparse an nonscaled n by p sparse matrix
#' @param y an n vector
#' @param cm a p vector of column means 
#' @param csd a p vector of column standard deviations 
compute_sparse_Xty = function(X.sparse, y, cm, csd){
  if(is.matrix(X.sparse)){
    return(t(X.sparse)%*%y)
  }else{
    Xty = crossprod(X.sparse, y)
    scaled.Xty = t(t(Xty)/csd)
    centered.scaled.Xty =scaled.Xty - cm/csd * sum(y)     
    return(as.numeric(centered.scaled.Xty))
  }
}

#' @title problem3: compute M%*%t(scaled.X) using sparse multiplication
#' @param M a L by p matrix
#' @param X.sparse an nonscaled n by p sparse matrix
#' @param cm a p vector of column means 
#' @param csd a p vector of column standard deviations 
compute_sparse_MtX = function(M, X.sparse, cm, csd){
  if(is.matrix(X.sparse)){
    return(M%*%t(X.sparse))
  }else{
    return(t(apply(M, 1, function(y) compute_sparse_Xy(X.sparse, y, cm, csd))))
  }
}