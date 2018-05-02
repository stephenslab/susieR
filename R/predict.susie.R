#' @title extract regression coefficients from susie fit
#' @param s a susie fit
#' @return a p vector of estimated regression coefficients
#' @method coef susie
#' @export
coef.susie = function(s){
  colSums(s$alpha*s$mu)
}

#' @title predict future observations or extract coefficients from susie fit
#' @param s a susie fit
#' @param newx a new value for X at which to do predictions
#' @param type if this is coefficients, then calls coef.susie
#' @method predict susie
#' @export
predict.susie = function(s,newx = NULL,type=c("response","coefficients")){
  type <- match.arg(type)
  if (type=="coefficients"){
    if(!missing(newx)){
      stop("Do not supply newx when predicting coefficients")
    }
    return(coef(s))
  }

  if(missing(newx)){return(s$Xr)}

  return(newx %*% coef(s))
}
