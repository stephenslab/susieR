#' @title Perform Univariate Linear Regression Separately for Columns of X
#' 
#' @description This function performs the univariate linear
#'   regression y ~ x separately for each column x of X. Each regression
#'   is implemented using \code{.lm.fit()}. The estimated effect size
#'   and stardard error for each variable are outputted.
#' 
#' @param X n by p matrix of regressors.
#' 
#' @param y n-vector of response variables.
#' 
#' @param Z Optional n by k matrix of covariates to be included in all
#'   regresions. If Z is not \code{NULL}, the linear effects of
#'   covariates are removed from y first, and the resulting residuals
#'   are used in place of y.
#' 
#' @param center If \code{center = TRUE}, center X, y and Z.
#' 
#' @param scale If \code{scale = TRUE}, scale X, y and Z.
#' 
#' @param return_residuals Whether or not to output the residuals if Z
#'   is not \code{NULL}.
#'
#' @return A list with two vectors containing the least-squares
#'   estimates of the coefficients (\code{betahat}) and their standard
#'   errors (\code{sebetahat}). Optionally, and only when a matrix of
#'   covariates \code{Z} is provided, a third vector \code{residuals}
#'   containing the residuals is returned.
#' 
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#' res = univariate_regression(X,y)
#' plot(res$betahat/res$sebetahat)
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
      Z = scale(Z,center = TRUE,scale = scale)
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
  if (return_residuals && !is.null(Z)) 
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
