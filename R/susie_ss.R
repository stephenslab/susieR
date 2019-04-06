#' @title Bayesian sum of single-effect (susie) linear regression using summary stat
#' @details Performs sum of single-effect (susie) linear regression of Y on X when
#' only summary statistics are available. The summary data required are
#' the p by p matrix X'X, the p vector X'Y, and the sample variance of Y, or (1/n)Y'Y. The summary stats should come from the same individuals.
#' Both the columns of X and the vector Y
#' should be centered to have mean 0 before
#' computing these summary statistics; you may also want to scale each column of X and Y to have variance 1 (see examples).
#' This function fits the regression model Y = sum_l Xb_l + e, where elements of e are iid N(0,var=residual_variance) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=scaled_prior_variance*var(Y)).
#' @param XtX a p by p matrix, X'X, where columns of X are centered to have mean 0
#' @param Xty a p vector, X'Y, where columns of X are centered and Y is centered to have mean 0
#' @param yty a scaler, Y'Y, where Y is centered to have mean 0
#' @param n sample size
#' @param L maximum number of non-zero effects
#' @param type sufficient or z scores
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param residual_variance the residual variance (defaults to variance of Y)
#' @param estimate_residual_variance indicates whether to estimate residual variance
#' @param estimate_prior_variance indicates whether to estimate prior
#' @param estimate_prior_method The method used for estimating prior variance, 'optim' or 'EM'
#' @param r_tol tolerance level for eigen value check of positive semidefinite matrix of R.
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param standardize logical flag (default=TRUE) for whether to standardize columns of X to unit variance prior to fitting.
#' Note that `scaled_prior_variance` specifies the prior on the coefficients of X *after* standardization (if performed).
#' If you do not standardize you may need
#' to think more carefully about specifying
#' `scaled_prior_variance`. Whatever the value of standardize, the coefficients (returned from `coef`) are for X on the original input scale.
#' @param max_iter maximum number of iterations to perform
#' @param s_init a previous susie fit with which to initialize
#' @param intercept_value a value to assign to the intercept (since the intercept cannot be estimated from centered summary data). This
#' value will be used by coef.susie() to assign an intercept value, for consistency with the non-summary-statistic version of this function \code{susie}.
#' Set to NULL if you want coef.susie() not to include an intercept term (and so only return a p vector).
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' @param tol convergence tolerance based on alpha
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param skip_checks whether to skip the checks for XtX and Xty
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{mu}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{mu2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{XtXr}{an p vector of t(X) times fitted values, the fitted values equal to X times colSums(alpha*mu))}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#'
#' @examples
#' set.seed(1)
#' n    <- 1000
#' p    <- 1000
#' beta <- rep(0,p)
#' beta[1:4] <- 1
#' X        <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' y        <- c(X %*% beta + rnorm(n))
#' input_ss <- compute_ss(X,y,standardize = TRUE)
#' res      <- with(input_ss,susie_ss(XtX,Xty,yty,n))
#' coef(res)
#'
#' @export
susie_ss = function(XtX, Xty, yty, n, L=10,
                    scaled_prior_variance=0.2,
                    residual_variance=NULL,
                    r_tol = 1e-08,
                    prior_weights = NULL, null_weight = NULL,
                    standardize = TRUE,
                    estimate_residual_variance = TRUE,
                    estimate_prior_variance = TRUE,
                    estimate_prior_method = c("optim","EM"),
                    max_iter=100,s_init = NULL, intercept_value=0,
                    coverage=0.95, min_abs_corr=0.5,
                    tol=1e-3, verbose=FALSE, track_fit = FALSE, skip_checks=FALSE){
  # Process input estimate_prior_method.
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Check input XtX.
  if (!(is.double(XtX) & is.matrix(XtX)) & !inherits(XtX,"CsparseMatrix"))
    stop("Input X must be a double-precision matrix, or a sparse matrix.")
  if (is.numeric(null_weight) && null_weight == 0) null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight<0 || null_weight>=1)
      stop('Null weight must be between 0 and 1')
    if (is.null(prior_weights))
      prior_weights = c(rep(1/ncol(XtX)*(1-null_weight), ncol(XtX)), null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight), null_weight)
    XtX = cbind(rbind(XtX, 0),0)
    Xty = c(Xty, 0)
  }

  p = ncol(XtX)

  if(standardize){
    dXtX = diag(XtX)
    csd = sqrt(dXtX/(n-1))
    csd[csd == 0] = 1
    XtX = (1/csd) * t((1/csd) * XtX)
    Xty = (1/csd) * Xty
  }else{
    csd = rep(1, length = p)
  }
  attr(XtX, "d") <- diag(XtX)
  attr(XtX, "scaled:scale") <- csd

  attr(XtX, 'eigen') <- eigen(XtX, symmetric = TRUE)
  attr(XtX, 'eigen')$values[abs(attr(XtX, 'eigen')$values) < r_tol] <- 0

  # check whether Xty in space spanned by the non-zero eigenvectors of XtX
  if(!skip_checks){
    A = attr(XtX, 'eigen')$vectors[,attr(XtX, 'eigen')$values!=0]
    in_space = all.equal(as.vector(A%*%solve(crossprod(A)) %*% crossprod(A, Xty)), Xty)
    if(in_space!=TRUE){
      warning('Xty does not lie in the space of non-zero eigenvectors of XtX')
    }
  }

  # initialize susie fit
  s = init_setup(0,p,L,scaled_prior_variance,residual_variance,
                 prior_weights,null_weight,yty/(n-1))
  s$Xr = NULL; s$XtXr = rep(0,p)

  if (!missing(s_init)) {
    s = modifyList(s, s_init)
    s = init_finalize(s, X=XtX)
    s$XtXr = s$Xr
    s$Xr = NULL
  } else {
    s = init_finalize(s)
  }

  # intialize elbo to NA
  elbo = rep(NA,max_iter+1)
  elbo[1] = -Inf;
  tracking = list()

  # alpha_new = s$alpha
  for(i in 1:max_iter){
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    # alpha_old = alpha_new
    s = update_each_effect_ss(XtX, Xty, s, estimate_prior_variance,estimate_prior_method)
    # alpha_new = s$alpha

    if(verbose){
      print(paste0("objective:",get_objective_ss(XtX, Xty, s, yty, n)))
    }
    if(estimate_residual_variance){
      est_sigma2 = estimate_residual_variance_ss(XtX,Xty,s,yty,n)
      if(est_sigma2 < 0){
        stop('Estimating residual variance failed: the estimated value is negative')
      }
      s$sigma2 = est_sigma2
      if(verbose){
        print(paste0("objective:",get_objective_ss(XtX, Xty, s, yty, n)))
      }
    }
    elbo[i+1] = get_objective_ss(XtX, Xty, s, yty, n)

    # if(max(abs(alpha_new - alpha_old)) < tol) break;

    if((elbo[i+1]-elbo[i])<tol) {
      s$converged = TRUE
      break;
    }
  }
  elbo = elbo[2:(i+1)] # Remove first (infinite) entry, and trailing NAs.
  s$elbo <- elbo
  s$niter <- i

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in", max_iter, "iterations!"))
    s$converged = FALSE
  }

  s$intercept = intercept_value
  s$Xtfitted = s$XtXr

  s$X_column_scale_factors = attr(XtX,"scaled:scale")

  if (track_fit)
    s$trace = tracking

  ## SuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    R = muffled_cov2cor(XtX)
    s$sets = susie_get_cs(s, coverage=coverage, Xcorr=R, min_abs_corr=min_abs_corr)
    s$pip = susie_get_pip(s, prune_by_cs = FALSE)
  }

  return(s)
}

#' @title Check input covariance / correlation matrix
#' @param R a p by p matrix of X'X, covariance matrix or correlation matrix
#' @return a verified matrix
#' @keywords internal

check_r_matrix <- function(R, expected_dim, r_tol) {
  n = nrow(R)
  if(n != expected_dim) {
    stop(paste0('The dimension of R (', n, ' by ', n, ') does not agree with expected (', expected_dim, ' by ', expected_dim, ')'))
  }
  if(!is_symmetric_matrix(R)){
    stop('R is not a symmetric matrix.')
  }

  E <- tryCatch(suppressWarnings(chol(R, pivot = TRUE, tol=r_tol)),error = function(e) FALSE)
  if (is.logical(E)) {
    stop('R is not a positive semidefinite matrix.')
  }

  if(any(is.na(R))){
    warning('R matrix contains NA, it is replaced with 0.')
    isnaR = is.na(R)
    naind = which(rowSums(isnaR) == ncol(R)-1)
    R[,naind] = 0
    R[naind,] = 0
    R[isnaR] = 0
  }
  X0 = diag(R) == 0
  # convert any input R to correlation matrix
  # if R has 0 colums and rows, cov2cor produces NaN and warning
  R = muffled_cov2cor(R)
  # change the columns and rows with NaN to 0
  if(sum(X0) > 0){
    R[X0, ] = R[,X0] = 0
  }
  return(R)
}

#' @title Summary statistics version of SuSiE on betahat, the corresponding standard error, and correlation (or covariance) matrix
#' @param bhat a p vector of estimated effects.
#' @param shat a p vector of corresponding standard errors.
#' @param R a p by p symmetric and positive semidefinite matrix. It can be X'X, covariance matrix or correlation matrix.
#' @param n sample size.
#' @param var_y the (sample) variance of y. If it is unknown, the coefficients (returned from `coef`) are on the standardized X, y scale.
#' @param r_tol tolerance level for eigen value check of positive semidefinite matrix of R.
#' @param L maximum number of non-zero effects.
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times)
#' @param estimate_residual_variance indicates whether to estimate residual variance
#' @param estimate_prior_variance indicates whether to estimate prior
#' @param estimate_prior_method The method used for estimating prior variance, 'optim' or 'EM'
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param standardize logical flag (default=TRUE) for whether to standardize columns of X to unit variance prior to fitting. It is useful when `var_y` is given.
#' Note that `scaled_prior_variance` specifies the prior on the coefficients of X *after* standardization (if performed).
#' If you do not standardize you may need
#' to think more carefully about specifying `scaled_prior_variance`.
#' Whatever the value of standardize, the coefficients (returned from `coef`) are on the bhat, shat scale.
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param ... further arguments to be passed to \code{\link{susie_ss}}
#' @return a susie fit
#'
#' @export
susie_bhat = function(bhat, shat, R, n, var_y = 1, r_tol = 1e-08,
                      L=10,
                      scaled_prior_variance=0.2,
                      estimate_residual_variance = TRUE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = c("optim", "EM"),
                      prior_weights = NULL, null_weight = NULL,
                      standardize = TRUE,
                      coverage=0.95, min_abs_corr=0.5,
                      verbose=FALSE, track_fit = FALSE, ...){
  if(missing(n)) {
    stop('n must be provided')
  }
  if(length(shat) == 1) {
    shat = rep(shat, length(bhat))
  }
  if(length(bhat) != length(shat)) {
    stop('The length of bhat does not agree with length of shat.')
  }
  if(anyNA(bhat) || anyNA(shat)){
    stop('The input summary statistics have missing value.')
  }
  #
  that = bhat/shat
  that[is.na(that)] = 0
  R2 = that^2/(that^2 + n-2)
  sigma2 = (n-1)*(1-R2)/(n-2)
  #
  R = check_r_matrix(R, length(bhat), r_tol)
  #
  if(missing(var_y)) {
    XtX = (n-1)*R
    Xty = sqrt(sigma2) * sqrt(n-1) * that
  } else {
    XtXdiag = var_y*sigma2/(shat^2)
    Xty = that * var_y* sigma2/shat
    XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
  }

  estimate_prior_method = match.arg(estimate_prior_method)
  susie_ss(XtX = XtX, Xty = Xty, yty = var_y * (n-1), n = n, L = L,
           scaled_prior_variance = scaled_prior_variance,
           estimate_residual_variance = estimate_residual_variance,
           estimate_prior_variance = estimate_prior_variance,
           estimate_prior_method = estimate_prior_method,
           prior_weights = prior_weights, null_weight = null_weight,
           standardize = standardize,
           coverage=coverage, min_abs_corr=min_abs_corr,
           verbose=verbose, track_fit = track_fit, ...)
}


