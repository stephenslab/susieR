#' @param s a susie fit
#' @return a p vector of estimated regression coefficients
#' @export
coef.susie = function(s){
  colSums(s$alpha*s$mu)
}

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
