#' @title Bayesian sum of single-effect (susie) linear regression of Y on X
#' @details Performs sum of single-effect (susie) linear regression of Y on X.
#' That is, this function
#' fits the regression model Y= sum_l Xb_l + e, where elements of e are iid N(0,residual_variance) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero. The prior on the non-zero element is N(0,var=var(Y)*scaled_prior_variance).
#' @param X an n by p matrix of covariates
#' @param Y an n vector
#' @param L maximum number of non-zero effects
#' @param scaled_prior_variance the scaled prior variance (vector of length L, or scalar. In latter case gets repeated L times). The prior variance on each non-zero element of b is set to be var(Y)*scaled_prior_variance.
#' @param residual_variance the residual variance (defaults to variance of Y)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param null_weight probability of no effect, for each single effect model
#' @param standardize logical flag (default=TRUE) for whether to standardize columns of X to unit variance prior to fitting.
#' Note that `scaled_prior_variance` specifies the prior on the coefficients of X *after* standardization (if performed).
#' If you do not standardize you may need
#' to think more carefully about specifying
#' `scaled_prior_variance`. Whatever the value of standardize, the coefficients (returned from `coef`) are for X on the original input scale.
#' Any column of X that has zero variance is not standardized, but left as is.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE). The latter is generally not recommended.
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param estimate_residual_variance indicates whether to estimate residual variance
#' @param estimate_prior_variance indicates whether to estimate prior (currently not recommended as not fully tested and assessed)
#' @param s_init a previous susie fit with which to initialize
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% confidence interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a confidence set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' @param verbose if true outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves current values of all iterations
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{mu}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{mu2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{Xr}{an n vector of fitted values, equal to X times colSums(alpha*mu))}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#' \item{elbo}{a vector of values of elbo achieved (objective function)}
#' \item{sets}{a list of `cs`, `purity` and selected `cs_index`}
#' \item{pip}{a vector of posterior inclusion probability}
#' \item{z}{a vector of univariate z-scores}
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1] = 1
#' beta[2] = 1
#' beta[3] = 1
#' beta[4] = 1
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' y = X %*% beta + rnorm(n)
#' res =susie(X,y,L=10)
#' coef(res)
#' plot(y,predict(res))
#' @export
susie = function(X,Y,L=10,scaled_prior_variance=0.2,residual_variance=NULL,
                 prior_weights=NULL, null_weight=NULL,
                 standardize=TRUE,intercept=TRUE,
                 estimate_residual_variance=TRUE,estimate_prior_variance = FALSE,
                 s_init = NULL,coverage=0.95,min_abs_corr=0.5,
                 max_iter=100,tol=1e-3,
                 verbose=FALSE,track_fit=FALSE) {
  # Check input X.
  if (!(is.double(X) & is.matrix(X)) & !is(X, 'CsparseMatrix'))
    stop("Input X must be a double-precision matrix, or a sparse matrix.")
  if (!missing(null_weight)) {
    if (null_weight<=0 || null_weight>=1)
      stop('Null weight must be between 0 and 1')
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(X)*(1-null_weight), ncol(X)), null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight), null_weight)
    X = cbind(X,0)
  }
  p = ncol(X)
  n = nrow(X)
  mean_y = mean(Y)
  # center and scale input.
  if(intercept){
    Y = Y-mean_y
  }
  X = safe_colScale(X,center=intercept, scale=standardize)
  # initialize susie fit
  s = init_setup(n,p,L,scaled_prior_variance,residual_variance,prior_weights,as.numeric(var(Y)))
  if (!missing(s_init)) {
    s = modifyList(s, s_init)
    s = init_finalize(s, X=X)
  } else {
    s = init_finalize(s)
  }

  #initialize elbo to NA
  elbo = rep(NA,max_iter+1)
  elbo[1] = -Inf;
  tracking = list()


  for(i in 1:max_iter){
    #s = add_null_effect(s,0)
    if (track_fit)
      tracking[[i]] = s
    s = update_each_effect(X, Y, s, estimate_prior_variance)
    if(verbose){
        print(paste0("objective:",susie_get_objective(X,Y,s)))
    }
    if(estimate_residual_variance){
      s$sigma2 = estimate_residual_variance(X,Y,s)
      if(verbose){
        print(paste0("objective:",susie_get_objective(X,Y,s)))
      }
    }
    #s = remove_null_effects(s)

    elbo[i+1] = susie_get_objective(X,Y,s)
    if((elbo[i+1]-elbo[i])<tol) break;
  }
  elbo = elbo[1:(i+1)] #remove trailing NAs
  s$elbo <- elbo
  s$niter <- i

  if(intercept){
    s$intercept = mean_y - sum(attr(X,"scaled:center")* (colSums(s$alpha*s$mu)/attr(X,"scaled:scale")))# estimate intercept (unshrunk)
    s$fitted = s$Xr + mean_y
  } else {
    s$intercept = 0
    s$fitted = s$Xr
  }

  if (track_fit)
    s$trace = tracking

  ## report z-scores from univariate regression
  s$z = calc_z(X,Y)
  ## SuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_CS(s, coverage=coverage, X=X, min_abs_corr=min_abs_corr)
    s$pip = susie_get_PIP(s,s$sets$cs_index)
  }
  ## for prediction
  s$X_column_scale_factors = attr(X,"scaled:scale")
  return(s)
}
