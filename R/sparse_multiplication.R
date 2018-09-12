# @title Compute scaled.X %*% y using sparse multiplication
# @param X is a scaled dense matrix or an unscaled sparse matrix 
# @param y a p vector
# @param cm a p vector of column means 
# @param csd a p vector of column standard deviations
#
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
compute_Xy = function(X, y){
  if (is.matrix(X)) {
    return(tcrossprod(X,t(y)))
  } else {
    #cm = attr(X, 'scaled:center')
    #csd = attr(X, 'scaled:scale')
    ##scale Xy
    #scaled.X  <- t(t(X)/csd)
    #scaled.Xy <- tcrossprod(scaled.X,t(y))
    ##center Xy
    #Xy <- scaled.Xy - sum(cm*y/csd) 
    #return(as.numeric(Xy))
    return(tcrossprod(attr(X, 'scaled.X'),t(y)))
 }
  #return(attr(X, 'scaled.X')%*%y)
}

# @title Compute t(scaled.X)%*%y using sparse multiplication
# @param X is a scaled dense matrix or an unscaled sparse matrix 
# @param y an n vector
# @param cm a p vector of column means 
# @param csd a p vector of column standard deviations
#
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
compute_Xty = function(X, y){
  if (is.matrix(X)) {
    return(crossprod(X,y))
  } else {
  cm = attr(X, 'scaled:center')
  csd = attr(X, 'scaled:scale')
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
compute_MtX = function(M, X){
  return(t(apply(M, 1, function(y) compute_Xy(X, y))))
}

# @title Compute square of a scaled X
# @param X is a scaled dense X, or an unscaled sparse X with scaled.X as one attribute
compute_X2 = function(X){
  if(is.matrix(X)){
    return(X*X)
  }else{
    scaled.X = attr(X, 'scaled.X')
    return(scaled.X*scaled.X) 
  }
}
  
  
  
  
  
