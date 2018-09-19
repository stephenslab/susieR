#' @title Compute scaled.X \%*\% b using sparse multiplication
#' @param X is a scaled dense matrix or an unscaled sparse matrix
#' @param b a p vector
#' @return an n vector
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
compute_Xb = function(X, b){
  if (is.matrix(X)) { #when X is a dense matrix
    return(tcrossprod(X,t(b))) #tcrossprod(A,B) performs A%*%t(B) but faster
  } else {
    ## FIXME: sparse matrix multiplication advantage will be diluted
    ## given larger p since using many transpose here.
    cm = attr(X, 'scaled:center')
    csd = attr(X, 'scaled:scale')
    #scale Xb
    scaled.Xb <- tcrossprod(X, t(b/csd))
    #center Xb
    Xb <- scaled.Xb - sum(cm*b/csd)
    return(as.numeric(Xb))
    #return(tcrossprod(attr(X, 'scaled.X'),t(b)))
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
  ytX        <- crossprod(y, X)
  scaled.Xty <- t(ytX/csd)
  centered.scaled.Xty <- scaled.Xty - cm/csd * sum(y)
  return(as.numeric(centered.scaled.Xty))
  }
}

#' @title Compute M\%*\%t(scaled.X) using sparse multiplication
#' @param M a L by p matrix
#' @param X is a scaled dense matrix or an unscaled sparse matrix
#' @return a L by n matrix
compute_MXt = function(M, X){
  return(t(apply(M, 1, function(b) compute_Xb(X, b))))
}
