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
#' @param n sample size
#' @param var_y the (sample) variance of the vector Y
#' @param L maximum number of non-zero effects
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param residual_variance the residual variance (defaults to variance of Y)
#' @param estimate_prior_variance indicates whether to estimate prior (currently not recommended as not working as well)
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
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% confidence interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a confidence set.
#' @param tol convergence tolerance based on alpha
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
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
#' res      <- with(input_ss,susie_ss(XtX,Xty,n,vary))
#' coef(res)
#'
#' @export
susie_ss = function(XtX, Xty, n, var_y = 1, L=10,
                    scaled_prior_variance=0.2,
                    residual_variance=NULL,
                    prior_weights = NULL, null_weight = NULL,
                    standardize = TRUE,
                    estimate_prior_variance = FALSE,
                    max_iter=100,s_init = NULL, intercept_value=0,
                    coverage=0.95, min_abs_corr=0.5,
                    tol=1e-4, verbose=FALSE, track_fit = FALSE){
  # Check input XtX.
  if (!(is.double(XtX) & is.matrix(XtX)) & !inherits(XtX,"CsparseMatrix"))
    stop("Input X must be a double-precision matrix, or a sparse matrix.")
  if (is.numeric(null_weight) && null_weight == 0) null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight<0 || null_weight>=1)
      stop('Null weight must be between 0 and 1')
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(XtX)*(1-null_weight), ncol(XtX)), null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight), null_weight)
    XtX = cbind(rbind(XtX, 0),0)
    Xty = c(Xty, 0)
  }

  p = ncol(XtX)

  if(standardize){
    dXtX = diag(XtX)
    dXtX[dXtX==0] = 1
    csd = sqrt(dXtX/(n-1))
    XtX = (1/csd) * t((1/csd) * XtX)
    Xty = (1/csd) * Xty
  }else{
    csd = rep(1, length = p)
  }
  attr(XtX, "d") <- diag(XtX)
  attr(XtX, "scaled:scale") <- csd

  # initialize susie fit
  s = init_setup(n,p,L,scaled_prior_variance,residual_variance,
                 prior_weights,null_weight,as.numeric(var_y))
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

  alpha_new = s$alpha
  for(i in 1:max_iter){
    if (track_fit)
      tracking[[i]] = s
    alpha_old = alpha_new
    s = update_each_effect_ss(XtX, Xty, s, estimate_prior_variance)
    alpha_new = s$alpha

    elbo[i+1] = susie_get_objective_ss(XtX, Xty, s, var_y, n)

    if(verbose){
      print(paste0("objective:",elbo[i+1]))
    }
    if(max(abs(alpha_new - alpha_old)) < tol) break;

    # if((elbo[i+1]-elbo[i])<tol) break;
  }
  elbo = elbo[1:(i+1)] #remove trailing NAs
  s$elbo <- elbo
  s$niter <- i

  s$intercept = intercept_value
  s$Xtfitted = s$XtXr

  s$X_column_scale_factors = attr(XtX,"scaled:scale")

  if (track_fit)
    s$trace = tracking

  ## SuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    R = muffled_cov2cor(XtX)
    s$sets = susie_get_CS(s, coverage=coverage, Xcorr=R, min_abs_corr=min_abs_corr)
    s$pip = susie_get_PIP(s,s$sets$cs_index)
  }

  return(s)
}

#' @title Apply susie with z scores and correlation matrix R
#' @param z a p vector of z scores.
#' @param R a p by p symmetric and positive semidefinite matrix
#' @param n sample size
#' @param var_y the (sample) variance of the y
#' @param R_type the type of the matrix R, one of 'Cor', 'XtX' or 'Cov'
#' @param L maximum number of non-zero effects
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param estimate_prior_variance indicates whether to estimate prior (currently not recommended as not working as well)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% confidence interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a confidence set.
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param ... further arguments to be passed to \code{\link{susie_ss}}
#' @return a susie fit
#'
#' @export
susie_z = function(z, R, n, var_y = 1, R_type = c('Cor', 'Cov', 'XtX'),
                   L=10, scaled_prior_variance=0.2,
                   estimate_prior_variance = TRUE,
                   prior_weights = NULL, null_weight = NULL,
                   coverage=0.95, min_abs_corr=0.5,
                   verbose=FALSE, track_fit = FALSE, ...){
  if(nrow(R) != length(z)){
    stop('The dimension of R does not agree with length of z.')
  }

  if(!isSymmetric(R)){
    stop('R is not a symmetric matrix.')
  }
  eigenvalues <- eigen(R, only.values = TRUE)$values
  eigenvalues[abs(eigenvalues) < 1e-08] <- 0
  if(any(eigenvalues < 0)){
    stop('R is not a positive semidefinite matrix.')
  }

  R_type = match.arg(R_type)
  if(R_type != 'XtX'){
    XtX = (n-1)*R
  }else{
    XtX = R
  }
  susie_ss(XtX = XtX, Xty = z*sqrt(diag(XtX)*var_y),
           n = n, var_y = var_y, L = L,
           scaled_prior_variance = scaled_prior_variance,
           estimate_prior_variance = estimate_prior_variance,
           prior_weights = prior_weights, null_weight = null_weight,
           coverage=coverage, min_abs_corr=min_abs_corr,
           verbose=verbose, track_fit = track_fit, ...)
}

#' @title Apply susie with bhat, shat and correlation matrix R
#' @param bhat a p vector of estimated effects.
#' @param shat a p vector of corresponding standard errors.
#' @param R a p by p symmetric and positive semidefinite matrix
#' @param n sample size
#' @param var_y the (sample) variance of the y
#' @param R_type the type of the matrix R, one of 'Cor', 'XtX' or 'Cov'
#' @param L maximum number of non-zero effects
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times )
#' @param estimate_prior_variance indicates whether to estimate prior (currently not recommended as not working as well)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% confidence interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a confidence set.
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @param ... further arguments to be passed to \code{\link{susie_ss}}
#' @return a susie fit
#'
#' @export
susie_bhat = function(bhat, shat, R, n, var_y = 1, R_type = c('Cor', 'Cov', 'XtX'),
                      L=10, standardize = TRUE,
                      scaled_prior_variance=0.2,
                      estimate_prior_variance = FALSE,
                      prior_weights = NULL, null_weight = NULL,
                      coverage=0.95, min_abs_corr=0.5,
                      verbose=FALSE, track_fit = FALSE, ...){
  if(length(bhat) != length(shat)){
    stop('The length of bhat does not agree with length of shat.')
  }
  if(nrow(R) != length(bhat)){
    stop('The dimension of R does not agree with length of bhat.')
  }

  R_type = match.arg(R_type)
  if(R_type != 'Cor'){
    R = cov2cor(R)
  }

  if(!isSymmetric(R)){
    stop('R is not a symmetric matrix.')
  }
  eigenvalues <- eigen(R, only.values = TRUE)$values
  eigenvalues[abs(eigenvalues) < 1e-08] <- 0
  if(any(eigenvalues < 0)){
    stop('R is not a positive semidefinite matrix.')
  }

  that = bhat/shat

  R2 = that^2/(that^2 + n-2)
  sigma2 = var_y*(n-1)*(1-R2)/(n-2)

  if(standardize){
    XtX = (n-1)*R
    Xty = sqrt(sigma2) * sqrt(n-1) * that
  }else{
    XtXdiag = sigma2/(shat^2)
    Xty = bhat * XtXdiag
    XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
  }

  susie_ss(XtX = XtX, Xty = Xty, n = n, var_y = var_y, L = L,
           scaled_prior_variance = scaled_prior_variance,
           estimate_prior_variance = estimate_prior_variance,
           prior_weights = prior_weights, null_weight = null_weight,
           coverage=coverage, min_abs_corr=min_abs_corr,
           verbose=verbose, track_fit = track_fit, ...)
}
