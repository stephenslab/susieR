#' @title Perform univariate regression between each column of X and y.
#'   Remove covariates if Z is not NULL.
#'
#' @importFrom stats lm
#' @importFrom stats .lm.fit
#' @importFrom stats coef
#' @importFrom stats summary.lm
#'
#' @export
#' 
univariate_regression = function (X, y, Z = NULL, center = TRUE,
                                  scale = FALSE, return_residuals = FALSE) {
  y_na = which(is.na(y))
  if (length(y_na)) {
    X = X[-y_na,]
    y = y[-y_na]
  }
  if (center) {
    y = y - mean(y)
    X = scale(X,center = TRUE,scale = scale)
  } else 
    X = scale(X,center = FALSE,scale = scale)
  X[is.nan(X)] = 0
  if (!is.null(Z)) {
    if (center)
      Z = scale(Z,center = TRUE,scale = FALSE)
    y = .lm.fit(Z,y)$residuals
  }
  output = try(do.call(rbind,
                       lapply(1:ncol(X), function (i) {
                         g = .lm.fit(cbind(1,X[,i]),y)
                         return(c(coef(g)[2],calc_stderr(cbind(1,X[,i]),
                                                         g$residuals)[2]))
                       })),
               silent = TRUE)
  
  # Exception occurs, fall back to a safer but slower calculation.
  if (inherits(output,"try-error")) {
    output = matrix(0,ncol(X),2)
    for (i in 1:ncol(X)) {
      fit = summary(lm(y ~ X[,i]))$coef
      if (nrow(fit) == 2)
        output[i,] = as.vector(summary(lm(y ~ X[,i]))$coef[2,1:2])
      else
        output[i,] = c(0,0)
    }
  }
  if (return_residuals) 
    return(list(betahat = output[,1],sebetahat = output[,2],residuals = y))
  else
    return(list(betahat = output[,1],sebetahat = output[,2]))
}

# Computes the z-scores (t-statistics) for association between Y and
# each column of X.
calc_z = function (X, Y, center = FALSE, scale = FALSE) {
  univariate_z = function(X,Y,center,scale) {
    out = univariate_regression(X,Y,center = center,scale = scale)
    return(out$betahat/out$sebetahat)
  }
  if (is.null(dim(Y)))
    return(univariate_z(X,Y,center,scale))
  else
    return(do.call(cbind,lapply(1:ncol(Y),
                                function(i) univariate_z(X,Y[,i],
                                                         center = center,
                                                         scale = scale))))
}
