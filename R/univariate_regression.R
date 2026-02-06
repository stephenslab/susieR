#' @title Perform Univariate Linear Regression Separately for Columns of X
#' 
#' @description This function performs the univariate linear
#'   regression y ~ x separately for each column x of X. The estimated effect size
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
#' @param method Either \dQuote{sumstats} (faster implementation) or
#'   \dQuote{lmfit} (uses \code{\link[stats]{.lm.fit}}).
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
                                  scale = FALSE, return_residuals = FALSE,
                                  method = c("lmfit", "sumstats")) {
  method <- match.arg(method)
  y_na <- which(is.na(y))
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

  # fast implementation: computes X'X and X'y without forming X
  if (method == "sumstats") {
    output <- try({
      n  <- length(y)
      sy <- sum(y)
      yy <- sum(y * y)
      p  <- ncol(X)
      res <- matrix(NA_real_, nrow = p, ncol = 2)
      
      for (i in seq_len(p)) {
        x   <- X[, i]
        sx  <- sum(x)
        sxx <- sum(x * x)
        sxy <- sum(x * y)

        # XtX and Xty for [1, x]
        # XtX = [[ n,  sx ],
        #        [ sx, sxx]]
        detXtX <- n * sxx - sx * sx
        if (!is.finite(detXtX) || detXtX <= 0) {
          warning_message("Column ", i, " has zero variance after centering/scaling")
          res[i, ] <- c(0, 0)  # constant/degenerate column
          next
        }

        XtX <- matrix(c(n, sx, sx, sxx), nrow = 2, ncol = 2)
        Xty <- c(sy, sxy)

        # Solve (XtX) beta = Xty via Cholesky
        R    <- chol(XtX)                              # XtX = R^T R
        beta <- backsolve(R, forwardsolve(t(R), Xty))  # slope is beta[2]

        # RSS = y'y - 2 beta^T X'y + beta^T XtX beta (no need to form
        # residuals)
        rss <- yy - 2 * sum(beta * Xty) +
               as.numeric(crossprod(beta, XtX %*% beta))
        sigma2 <- rss / (n - 2)         # p = 2 (intercept + slope)

        # Var(beta) = sigma2 * (XtX)^{-1}; se(slope) = sqrt( ... [2,2] )
        XtX_inv   <- chol2inv(R)
        se_slope  <- sqrt(sigma2 * XtX_inv[2, 2])
        res[i, ] <- c(beta[2], se_slope)
      }
      res
    }, silent = TRUE)
  } else {
    # original .lm.fit-based implementation
    output = try(do.call(rbind,
                        lapply(1:ncol(X), function (i) {
                          g = .lm.fit(cbind(1,X[,i]),y)
                          return(c(coef(g)[2],calc_stderr(cbind(1,X[,i]),
                                                          g$residuals)[2]))
                        })),
                silent = TRUE)
  }

  # Exception occurs, fall back to a safer but slower calculation.
  if (inherits(output,"try-error")) {
    output = matrix(0,ncol(X),2)
    for (i in 1:ncol(X)) {
      fit = summary(lm(y ~ X[,i]))$coef
      if (nrow(fit) == 2)
        output[i,] = as.vector(summary(lm(y ~ X[,i]))$coef[2,1:2])
      else {
        warning_message("Column ", i, " has zero variance after centering/scaling")
        output[i,] = c(0,0)
      }
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


# ----------------------------------------------------------------------
# Some miscellaneuous auxiliary functions are listed below.
# Some functions are directly copied from varbvs,
# https://github.com/pcarbo/varbvs
# ----------------------------------------------------------------------

# Remove covariate effects Regresses Z out from X and y; that is, X
# and y are projected into the space orthogonal to Z.
#' 
#' @importFrom Matrix forceSymmetric
#'
remove_covariate <- function (X, y, Z, standardize = FALSE, intercept = TRUE) {
  
  # check if Z is null and intercept = FALSE
  if (is.null(Z) & (intercept == FALSE)) {
    return(list(X = X, y = y, Z = Z,
                ZtZiZX = rep(0,dim(X)[2]), ZtZiZy = 0))
  }
  
  # redefine y
  y = c(as.double(y))
  n = length(y)
  
  # add intercept if intercept = TRUE
  if (intercept) {
    if (is.null(Z))
      Z <- matrix(1,n,1)
    else
      Z <- cbind(1,Z)
  }
  
  if (ncol(Z) == 1) {
    ZtZ         = forceSymmetric(crossprod(Z))       # (Z^T Z) symmetric
    ZtZiZy      = as.vector(solve(ZtZ,c(y %*% Z)))   # (Z^T Z)^{-1} Z^T y
    ZtZiZX      = as.matrix(solve(ZtZ,t(Z) %*% X))   # (Z^T Z)^{-1} Z^T X
    X           = scale(X, center = intercept, scale = standardize)
    alpha       = mean(y)
    y           = y - alpha
    
  } else {
    ZtZ         = forceSymmetric(crossprod(Z))       # (Z^T Z) symmetric
    ZtZiZy      = as.vector(solve(ZtZ,c(y %*% Z)))   # (Z^T Z)^{-1} Z^T y
    ZtZiZX      = as.matrix(solve(ZtZ,t(Z) %*% X))   # (Z^T Z)^{-1} Z^T X
    
    #   y = y - Z (Z^T Z)^{-1} Z^T y
    #   X = X - Z (Z^T Z)^{-1} Z^T X  
    y     = y - c(Z %*% ZtZiZy)
    X     = X - Z %*% ZtZiZX
  }
  
  return(list(X = X, y = y, Z = Z,
              ZtZiZX = ZtZiZX, ZtZiZy = ZtZiZy))
}

#' @title Ordering of Predictors from Univariate Regression
#' 
#' @description This function extracts the ordering of the predictors
#'   according to the coefficients estimated in a basic univariate
#'   regression; in particular, the predictors are ordered in decreasing
#'   order by magnitude of the univariate regression coefficient
#'   estimate.
#' 
#' @param X An input design matrix. This may be centered and/or
#'   standardized prior to calling function.
#' 
#' @param y A vector of response variables.
#'
#' @return An ordering of the predictors.
#' 
#' @examples
#' ### generate synthetic data
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#' 
#' univ.order = univar.order(X,y)
#' 
#' @export
#' 
univar.order = function(X, y) {
  colnorm = c(colMeans(X^2))
  return (order(abs(c(t(X) %*% y) / colnorm), decreasing = TRUE))
}

#' @title Ordering of Predictors from Coefficient Estimates 
#' 
#' @param beta A vector of estimated regression coefficients.
#' 
#' @description This function orders the predictors by decreasing
#'   order of the magnitude of the estimated regression coefficient.
#'
#' @return An ordering of the predictors.
#' 
#' @examples
#' ### generate synthetic data
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#' 
#' ### glmnet fit
#' library(glmnet)
#' beta.lasso = coef(cv.glmnet(X, y))[-1]
#' lasso.order = absolute.order(beta.lasso)
#' 
#' ### ncvreg fit
#' library(ncvreg)
#' beta.scad = c(coef(cv.ncvreg(X, y))[-1])
#' scad.order = absolute.order(beta.scad)
#' 
#' @export
#' 
absolute.order = function (beta) {
  abs_order = c(order(abs(beta), decreasing = TRUE))
  return (abs_order)
}

#' @title Ordering of Predictors by Regularization Path
#' 
#' @param fit The output of a function such as \code{glmnet} from the
#'   \code{glmnet} package or \code{ncvreg} from the \code{ncvfeg} that
#'   estimates a "regularization path" for all predictors.
#' 
#' @description This function determines an ordering of the predictors
#'  based on the regularization path of the penalized regression; in
#'   particular, the predictors are ordered based on the order in which
#'   the coefficients are included in the model as the penalty strength
#'   decreases.
#' 
#' @return An ordering of the predictors.
#' 
#' @examples
#' ### generate synthetic data
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#' 
#' ### glmnet fit
#' library(glmnet)
#' fit.lasso = glmnet(X, y)
#' lasso.order = path.order(fit.lasso)
#' 
#' ### ncvreg fit
#' library(ncvreg)
#' fit.scad = ncvreg(X, y)
#' scad.order = path.order(fit.scad)
#'
#' @export
#' 
path.order = function (fit) {
  beta_path = coef(fit)[-1,]
  K = dim(beta_path)[2]
  path_order = c()
  for (k in 1:K) {
    crt_path = which(beta_path[,k] != 0)
    if (length(crt_path) != 0 & length(path_order) == 0) {
      path_order = c(path_order, crt_path)
    } else if(length(crt_path) != 0) {
      path_order = c(path_order, crt_path[-which(crt_path %in% path_order)] )
    }
  }
  path_order = unname(path_order)
  index_order = c(path_order, seq(1,dim(beta_path)[1])[-path_order])
  return (index_order)
}
