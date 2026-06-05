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

#' @rdname univariate_regression
#' @export
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


#' @title Per-Position Marginal OLS Regression of `Y` on Each Column of `X`
#'
#' @description Computes the marginal OLS regression coefficient and
#'   standard error for each `(X column, Y column)` pair, treating
#'   the regressions as independent. `X` is assumed column-centred
#'   (no intercept term in the per-pair regression); each `Y`
#'   column is treated independently. Returns the J x T matrices
#'   `Bhat` and `Shat`.
#'
#' Used internally by single-effect-regression style routines that
#' need a per-position marginal estimate. Vectorised across columns
#' of `Y` so callers can pass either a numeric vector (T = 1) or a
#' numeric matrix (T > 1) without looping at the call site.
#'
#' @param X numeric matrix `n x J`, expected column-centred.
#' @param Y numeric matrix `n x T` or numeric vector of length `n`.
#'   When a vector, is treated as a one-column matrix.
#' @param predictor_weights optional numeric vector of length `J`
#'   giving `colSums(X^2)`. Computed internally when `NULL`.
#'   Callers that have this cached on the data object pass it
#'   through to avoid recomputation.
#' @param sigma2 optional numeric scalar giving a known residual
#'   variance. When supplied, `Shat[j, t] = sqrt(sigma2 /
#'   predictor_weights[j])` (single-effect-residual form). When
#'   `NULL`, `Shat` is the per-pair empirical residual standard
#'   error: for each `(j, t)` pair, the sample SD of `Y[, t] -
#'   X[, j] * Bhat[j, t]` divided by `sqrt(n - 1)`. The latter
#'   matches the form used by data-driven prior init routines
#'   (e.g., for fitting a normal-mixture prior via `ashr::ash`).
#'
#' @return list with elements `Bhat` (`J x T`) and `Shat` (`J x T`).
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(50 * 5), 50, 5)
#' X <- scale(X, center = TRUE, scale = FALSE)
#' Y <- matrix(rnorm(50 * 3), 50, 3)
#' out <- compute_marginal_bhat_shat(X, Y)
#' dim(out$Bhat)   # 5 x 3
#' dim(out$Shat)   # 5 x 3
#'
#' @importFrom Rfast colVars
#' @export
compute_marginal_bhat_shat <- function(X, Y,
                                       predictor_weights = NULL,
                                       sigma2 = NULL) {
  if (is.null(dim(Y))) {
    Y <- matrix(Y, ncol = 1)
  }
  n <- nrow(Y)
  J <- ncol(X)
  T_y <- ncol(Y)

  if (is.null(predictor_weights)) {
    predictor_weights <- colSums(X^2)
  }

  Bhat <- crossprod(X, Y) / predictor_weights      # J x T

  if (!is.null(sigma2)) {
    Shat <- matrix(sqrt(sigma2 / predictor_weights), nrow = J, ncol = T_y)
  } else {
    Shat <- vapply(
      seq_len(T_y),
      function(t) Rfast::colVars(Y[, t] - sweep(X, 2, Bhat[, t], "*")),
      numeric(J)
    )
    if (!is.matrix(Shat)) Shat <- matrix(Shat, nrow = J, ncol = T_y) # nocov — vapply always returns matrix
    Shat <- sqrt(pmax(Shat, 1e-64)) / sqrt(n - 1)
  }

  list(Bhat = Bhat, Shat = Shat)
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
    # Still scale the columns of X when standardize = TRUE, so callers that
    # rely on attr(X, "scaled:scale") to map coefficients back to the
    # original scale (e.g. mr.ash) do not divide by a missing attribute.
    if (standardize)
      X = scale(X, center = FALSE, scale = TRUE)
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
#' ### order predictors by magnitude of univariate regression coefficient
#' beta.hat    = univariate_regression(X,y)$betahat
#' order       = absolute.order(beta.hat)
#'
#' @export
#'
absolute.order = function (beta) {
  abs_order = c(order(abs(beta), decreasing = TRUE))
  return (abs_order)
}

#' @title Ordering of Predictors by Regularization Path
#' 
#' @param fit A fit object whose \code{coef()} method returns a matrix of
#'   coefficients with the intercept in the first row and one column per
#'   penalty strength (as produced by typical penalized-regression
#'   implementations).
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
#' p           = 30
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#'
#' ### build a minimal example 'fit' object with the same structure as a
#' ### fit from a penalized regression: a coefficient matrix with the
#' ### intercept in row 1 and one column per (decreasing) penalty value.
#' beta_path   = matrix(0, p + 1, p)
#' for (k in 1:p) beta_path[k + 1, k:p] = 1
#' fit         = list(coefficients = beta_path)
#' order       = path.order(fit)
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
