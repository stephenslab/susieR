#' @title Extract regression coefficients from susie fit
#' 
#' @param object A susie fit.
#'
#' @param \dots Additional arguments passed to the default \code{coef}
#'   method.
#' 
#' @return A p+1 vector, the first element being an intercept, and the
#'   remaining p elements being estimated regression coefficients.
#'
#' @importFrom stats coef
#'
#' @method coef susie
#' @export coef.susie
#' @export
#' 
coef.susie = function (object, ...) {
  s = object
  return(c(s$intercept,colSums(s$alpha*s$mu)/s$X_column_scale_factors))
}

#' @title Predict outcomes or extract coefficients from susie fit.
#' 
#' @param object A susie fit.
#' 
#' @param newx A new value for X at which to do predictions.
#' 
#' @param type The type of output. For \code{type = "response"},
#'   predicted or fitted outcomes are returned; for \code{type =
#'   "coefficients"}, the estimated coefficients are returned.
#'
#' @return For \code{type = "response"}, predicted or fitted outcomes
#'   are returned; for \code{type = "coefficients"}, the estimated
#'   coefficients are returned.
#' 
#' @importFrom stats coef
#'
#' @method predict susie
#' @export predict.susie
#' @export
#' 
predict.susie = function (object, newx = NULL,
                          type = c("response","coefficients"), ...) {
  s = object
  type = match.arg(type)
  if (type == "coefficients") {
    if (!missing(newx))
      stop("Do not supply newx when predicting coefficients")
    return(coef(s))
  }
  if (missing(newx))
    return(s$fitted)
  return(drop(s$intercept + newx %*% coef(s)[-1]))
}
