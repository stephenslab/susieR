#' @title Compute scaled.X \%*\% b using sparse multiplication
#' @param X is a scaled dense matrix or an unscaled sparse matrix
#' @param b a p vector
#' @return an n vector
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
#' @keywords internal
compute_Xb = function(X, b){
  if (is.matrix(X)) { #when X is a dense matrix
    return(tcrossprod(X,t(b))) #tcrossprod(A,B) performs A%*%t(B) but faster
  } else {
    ## given larger p since using many transpose here.
    cm = attr(X, 'scaled:center')
    csd = attr(X, 'scaled:scale')
    #scale Xb
    #when trend filtering
    if (!is.null(attr(X,"matrix.type"))) scaled.Xb <- compute_tf_Xb(attr(X, 'order'), b/csd)
    #when X is sparse
    else scaled.Xb <- tcrossprod(X, t(b/csd))
    #center Xb
    Xb <- scaled.Xb - sum(cm*b/csd)
    return(as.numeric(Xb))
 }
}

#' @title Compute t(scaled.X)\%*\%y using sparse multiplication
#' @param X is a scaled dense matrix or an unscaled sparse matrix
#' @param y an n vector
#' @return a p vector
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
compute_Xty = function(X, y){
  if (is.matrix(X)) { #when X is a dense matrix
    return(crossprod(X,y)) #corssprod(A,B) performs t(A)%*%B but faster
  } else { #when X is sparse matrix
  cm = attr(X, 'scaled:center')
  csd = attr(X, 'scaled:scale')
  ytX <- crossprod(y, X)
  #when trend filtering 
  if (!is.null(attr(X,"matrix.type"))) scaled.Xty <- compute_tf_Xty(attr(X, 'order'),y)/csd
  #when X is sparse
  else scaled.Xty <- t(ytX/csd)
  centered.scaled.Xty <- scaled.Xty - cm/csd * sum(y)
  return(as.numeric(centered.scaled.Xty))
  }
}

#' @title Compute M\%*\%t(scaled.X) using sparse multiplication
#' @param M a L by p matrix
#' @param X is a scaled dense matrix or an unscaled sparse matrix
#' @return a L by n matrix
#' @importFrom Matrix t
compute_MXt = function(M, X){
  if(is.matrix(X)){
    return(tcrossprod(M,X))
  }else{
    cm = attr(X, 'scaled:center')
    csd = attr(X, 'scaled:scale')
    #when trend filtring
    if (!is.null(attr(X,"matrix.type"))) {
      return(as.matrix(t(apply(M,1,function(b) compute_Xb(X, b)))))
    }
    #when X is sparse
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
}



