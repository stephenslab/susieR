#' @title extract regression coefficients from susie fit
#' @param object a susie fit
#' @return a p+1 vector, the first element being an intercept, and the remaining p elements being estimated regression coefficients
#' @export coef.susie
#' @export
coef.susie = function(object, ...){
  s <- object
  c(s$intercept,colSums(s$alpha*s$mu)/s$X_column_scale_factors)
}

#' @title predict future observations or extract coefficients from susie fit
#' @param object a susie fit
#' @param newx a new value for X at which to do predictions
#' @param type if this is coefficients, then calls coef.susie
#' @importFrom stats coef
#' @export predict.susie
#' @export
predict.susie = function(object,newx = NULL,
                         type=c("response","coefficients"),...) {
  s <- object
  type <- match.arg(type)
  if (type=="coefficients"){
    if(!missing(newx)){
      stop("Do not supply newx when predicting coefficients")
    }
    return(coef(s))
  }

  if(missing(newx)){return(s$fitted)}

  return(drop(s$intercept + newx %*% coef(s)[-1]))
}

#' @title predict future observations from susie_stumps fit
#' @param x a susie fit
#' @param newx a new value for X at which to do predictions
#' @param Xtrain the value of X used to train
#' @param include_linear boolean for whether linear term was included in training
#' @export
susie_stumps_predict = function(s,newx,Xtrain,include_linear=TRUE){
  xl = make_stumps_matrix(newx,include_linear,Xtrain)
  compute_Xb(xl,coef(s)[-1])+s$intercept
  # need to think how best to do this - need to create a matrix with same breaks
  # as Xtrain, but with data from Xtest
}

