#' @title Bayesian sum of single-effect (susie) linear regression using sufficient summary stat
#' @details Performs sum of single-effect (susie) linear regression of y on X when
#' only sufficient summary statistics are available. The sufficient summary data required are EITHER
#' the p vector bhat, the p vector shat, the p by p symmetric and positive semidefinite correlation
#' (or covariance) matrix R, the sample size n, the variance of y; OR
#' the p by p matrix X'X, the p vector X'y, the sum of squares of y (y'y) and the sample size.
#' The sufficient summary stats should come from the same individuals.
#' Both the columns of X and the vector y
#' should be centered to have mean 0 before
#' computing these summary statistics; you may also want to scale each column of X and y to have variance 1 (see examples).
#' This function fits the regression model y = sum_l Xb_l + e, where elements of e are iid N(0,var=residual_variance) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=scaled_prior_variance*y'y/(n-1)).
#' @param bhat a p vector of estimated effects.
#' @param shat a p vector of corresponding standard errors.
#' @param R a p by p symmetric and positive semidefinite matrix. It can be X'X, covariance matrix (X'X/(n-1)) or correlation matrix.
#' It should from the same samples used to compute `bhat` and `shat`. Using out of sample matrix may produce unreliable results.
#' @param n sample size.
#' @param var_y the (sample) variance of y, defined as y'y/(n-1) . If it is unknown, the coefficients (returned from `coef`) are on the standardized X, y scale.
#' @param XtX a p by p matrix, X'X, where columns of X are centered to have mean 0
#' @param Xty a p vector, X'y, where columns of X are centered and y is centered to have mean 0
#' @param yty a scaler, y'y, where y is centered to have mean 0
#' @param maf minor allele frequency; to be used along with `maf_thresh` to filter input summary statistics
#' @param maf_thresh variants having MAF smaller than this threshold will be filtered out
#' @param L maximum number of non-zero effects
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times)
#' @param residual_variance the residual variance (defaults to variance of y)
#' @param estimate_residual_variance indicates whether to estimate residual variance
#' @param estimate_prior_variance indicates whether to estimate prior
#' @param estimate_prior_method The method used for estimating prior variance.
#' @param check_null_threshold when prior variance is estimated, compare the
#' estimate with the null and set prior variance to null (zero) unless the log-likelihood
#' using the estimate is larger than that of null by this threshold. For example,
#' you can set it to 0.1 to nudge the estimate towards zero. Default is 0. Notice that setting it to non-zero
#' may lead to decreasing ELBO in some cases.
#' @param prior_tol when prior variance is estimated, compare the estimated value to this tol at the end of
#' the analysis and exclude a single effect from PIP computation if the estimated prior variance is smaller than it.
#' @param r_tol tolerance level for eigen value check of positive semidefinite matrix of R.
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param standardize logical flag (default=TRUE) for whether to adjust XtX and Xty such that they are computed from column standardized X, prior to fitting.
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
#' @param check_input whether to perform checks on XtX and Xty, the checks are:
#'
#' 1. Check whether XtX is positive semidefinite
#'
#' 2. Check whether Xty in space spanned by the non-zero eigenvectors of XtX
#'
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
#' ss <- susieR:::univariate_regression(X, y)
#' R <- with(input_ss, cov2cor(XtX))
#' res1      <- with(input_ss,susie_suff_stat(XtX = XtX,Xty = Xty, yty = yty,n=n))
#' coef(res1)
#' res2      <- with(ss,susie_suff_stat(bhat = betahat, shat = sebetahat, R = R, n=n, var_y = var(y)))
#' coef(res2)
#'
#' @export
susie_suff_stat = function(bhat, shat, R, n, var_y = 1,
                           XtX, Xty, yty, maf=NULL, maf_thresh=0,
                           L=10,
                           scaled_prior_variance=0.2,
                           residual_variance=NULL,
                           estimate_residual_variance = TRUE,
                           estimate_prior_variance = TRUE,
                           estimate_prior_method = c("optim","EM","simple"),
                           check_null_threshold=0, prior_tol = 1E-9,
                           r_tol = 1e-08,
                           prior_weights = NULL, null_weight = NULL,
                           standardize = TRUE,
                           max_iter=100,s_init = NULL, intercept_value=0,
                           coverage=0.95, min_abs_corr=0.5,
                           tol=1e-3, verbose=FALSE, track_fit = FALSE, check_input=FALSE){
  # Process input estimate_prior_method.
  estimate_prior_method <- match.arg(estimate_prior_method)

  if(missing(n)) {
    stop('n must be provided')
  }

  # Check sufficient stat input
  ## check bhat, shat
  missing_bhat = c(missing(bhat), missing(shat), missing(R))
  ## check XtX, Xty
  missing_XtX = c(missing(XtX), missing(Xty), missing(yty))

  if(all(missing_bhat) & all(missing_XtX)){
    stop('Please provide EITHER "bhat, shat, R, n, var(y)" OR "XtX, Xty, yty, n".')
  }
  if(any(missing_bhat) & any(missing_XtX)){
    stop('Please provide EITHER "bhat, shat, R, n, var(y)" OR "XtX, Xty, yty, n".')
  }
  if(all(missing_bhat) & any(missing_XtX)){
    stop('Please provide all "XtX, Xty, yty, n"')
  }
  if(all(missing_XtX) & any(missing_bhat)){
    stop('Please provide all "bhat, shat, R, n, var(y)"')
  }
  if((!any(missing_XtX)) & (!all(missing_bhat))){
    warning('Only use information from "XtX, Xty, yty, n"')
  }
  if((!any(missing_bhat))){
    if(!all(missing_XtX)){
      warning('Only use information from "bhat, shat, R, n, var(y)"')
    }

    # Compute XtX, Xty, yty from bhat, shat, R, n, var(y)
    if(length(shat) == 1) {
      shat = rep(shat, length(bhat))
    }
    if(length(bhat) != length(shat)) {
      stop('The length of bhat does not agree with length of shat.')
    }

    if(anyNA(bhat) || anyNA(shat)){
      stop('The input summary statistics have missing value.')
    }
    if(any(shat==0)){
      stop('shat contains zero.')
    }
    #
    that = bhat/shat
    that[is.na(that)] = 0
    R2 = that^2/(that^2 + n-2)
    sigma2 = (n-1)*(1-R2)/(n-2)

    # convert any input R to correlation matrix
    # if R has 0 colums and rows, cov2cor produces NaN and warning
    X0 = diag(R) == 0
    R = muffled_cov2cor(R)
    # change the columns and rows with NaN to 0
    if(sum(X0) > 0){
      R[X0, ] = R[,X0] = 0
    }

    #
    if(missing(var_y)) {
      XtX = (n-1)*R
      Xty = sqrt(sigma2) * sqrt(n-1) * that
    }else{
      XtXdiag = var_y*sigma2/(shat^2)
      Xty = that * var_y* sigma2/shat
      XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
    }
    yty = var_y * (n-1)
  }


  # Check input XtX.
  if(ncol(XtX) != length(Xty)) {
    stop(paste0('The dimension of XtX (', nrow(XtX), ' by ', ncol(XtX), ') does not agree with expected (', length(Xty), ' by ', length(Xty), ')'))
  }
  if(!is_symmetric_matrix(XtX)){
    stop('XtX is not a symmetric matrix.')
  }
  # MAF filter
  if(!is.null(maf)){
    if(length(maf) != length(Xty)){
      stop(paste0('The length of maf does not agree with expected ', length(Xty)))
    }
    id = which(maf > maf_thresh)
    XtX = XtX[id, id]
    Xty = Xty[id]
  }

  if(any(is.infinite(Xty))){
    stop('Xty contains infinite value.')
  }
  if (!(is.double(XtX) & is.matrix(XtX)) & !inherits(XtX,"CsparseMatrix"))
    stop("Input X must be a double-precision matrix, or a sparse matrix.")
  if(any(is.na(XtX)))
    stop('XtX matrix contains NA.')

  if(check_input){
    # Check whether XtX is positive semidefinite
    semi_pd = check_semi_pd(XtX, r_tol)
    if(semi_pd$status == FALSE){
      stop('XtX is not a positive semidefinite matrix.')
    }
    # Check whether Xty in space spanned by the non-zero eigenvectors of XtX
    proj = check_projection(semi_pd$matrix, Xty)
    if(proj$status == FALSE)
      warning('Xty does not lie in the space of non-zero eigenvectors of XtX')
  }

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
    s = update_each_effect_ss(XtX, Xty, s, estimate_prior_variance,estimate_prior_method,check_null_threshold)
    # alpha_new = s$alpha

    if(verbose){
      print(paste0("objective:",get_objective_ss(XtX, Xty, s, yty, n)))
    }
    # compute objective before updating residual variance
    # because part of the objective s$kl has already been computed
    # under the residual variance before the update
    elbo[i+1] = get_objective_ss(XtX, Xty, s, yty, n)
    # if(max(abs(alpha_new - alpha_old)) < tol) break;
    if((elbo[i+1]-elbo[i])<tol) {
      s$converged = TRUE
      break;
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
    s$pip = susie_get_pip(s, prune_by_cs = FALSE, prior_tol = prior_tol)
  }

  return(s)
}

#' @title Check whether A is positive semidefinite
#' @param A a symmetric matrix
#' @return a list of result: \cr
#' \item{matrix}{The matrix with eigen decomposition}
#' \item{status}{whether A is positive semidefinite}
#' \item{eigenvalues}{eigenvalues of A truncated by r_tol}
#' @export
check_semi_pd <- function(A, tol){
  attr(A, 'eigen') = eigen(A, symmetric = TRUE)

  eigenvalues = attr(A, 'eigen')$values
  eigenvalues[abs(eigenvalues) < tol] <- 0

  # E <- tryCatch(suppressWarnings(chol(R, pivot = TRUE, tol=r_tol)),error = function(e) FALSE)
  # if (is.logical(E)) {
  #   stop('R is not a positive semidefinite matrix.')
  # }

  return(list(matrix = A, status = !any(eigenvalues < 0), eigenvalues = eigenvalues))
}

#' @title Check whether b in space spanned by the non-zero eigenvectors of A
#' @param A a p by p matrix
#' @param b a length p vector
#' @return a list of result: \cr
#' \item{status}{whether b in space spanned by the non-zero eigenvectors of A}
#' \item{msg}{msg gives the difference between the projected b and b if status is FALSE}
#' @export
check_projection <- function(A, b){
  if(is.null(attr(A, 'eigen')))
    attr(A, 'eigen') = eigen(A, symmetric = TRUE)

  B = attr(A, 'eigen')$vectors[,attr(A, 'eigen')$values > -sqrt(.Machine$double.eps)]
  msg = all.equal(as.vector(B %*% crossprod(B, b)), b, check.names=FALSE)

  if (length(msg)==1 && msg == TRUE)
    return(list(status=T, msg=NA))
  else
    return(list(status=F, msg=msg))
}
