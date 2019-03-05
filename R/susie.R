#' @title SUm of Single Effects (SuSiE) Regression
#' 
#' @description Performs Bayesian multiple linear regression of Y on
#'   X. That is, this function fits the regression model \eqn{Y = sum_l
#'   Xb_l + e}, where elements of e are \emph{i.i.d.} normal with zero
#'   mean and variance \code{residual_variance}, and \eqn{sum_l b_l} is
#'   a vector of length p representing the effects to be estimated.  The
#'   SuSiE assumption is that each b_l has exactly one non-zero
#'   element. The prior on the non-zero element is normal with zero mean
#'   and variance \code{var(Y)*scaled_prior_variance}.
#'
#'   The model is fitted using the "Iterative Bayesian Stepwise
#'   Selection" (IBSS) algorithm.
#'
#' @param X An n by p matrix of covariates.
#' 
#' @param Y A vector of length n.
#' 
#' @param L Number of components (nonzero elements) in the SuSiE
#'   regression model. If \code{L} is larger than the number of
#'   covariates, p, \code{L} is set to p.
#' 
#' @param scaled_prior_variance The scaled prior variance. This is
#'   either a scalar, or a vector of length \code{L}. The prior variance
#'   on each non-zero element of b is set to be
#'   \code{var(Y)*scaled_prior_variance}.
#' 
#' @param residual_variance The variance of the residual. By default,
#'   it is set to \code{var(Y)}.
#' 
#' @param prior_weights A vector of length p, in which each entry
#'   gives the prior probability that corresponding column of X has a
#'   nonzero effect on the outcome, Y.
#' 
#' @param null_weight Prior probability of no effect, for each single
#'   effect model.
#' 
#' @param standardize If \code{standardize = TRUE}, standardize the
#'   columns of X to unit variance prior to fitting. Note that
#'   `scaled_prior_variance` specifies the prior on the coefficients of
#'   X \emph{after} standardization (if it is performed). If you do not
#'   standardize, you may need to think more carefully about specifying
#'   \code{scaled_prior_variance}. Whatever your choice, the
#'   coefficients returned by \code{coef} are given for \code{X} on the
#'   original input scale. Any column of \code{X} that has zero variance is
#'   not standardized, but left as is.
#' 
#' @param intercept If \code{intercept = TRUE}, the intercept is
#'   fitted; otherwise, it is set to zero. Setting \code{intercept =
#'   FALSE} is generally not recommended.
#' 
#' @param estimate_residual_variance If
#'   \code{estimate_residual_variance = TRUE}, the variance of the
#'   residual is estimated.
#' 
#' @param estimate_prior_variance If \code{estimate_prior_variance =
#' TRUE}, the prior variance is estimated father than fixed.
#' 
#' @param optimV_method The method used for estimating prior
#' variance. Posible settings are "EM" and "optim".
#' 
#' @param s_init a previous susie fit with which to initialize
#' 
#' @param coverage A number between 0 and 1 specifying the coverage of
#'  the estimated confidence sets.
#' 
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' 
#' @param compute_univariate_zscore if true, outputs z-score from per variable univariate regression
#' 
#' @param max_iter maximum number of iterations to perform
#' 
#' @param tol convergence tolerance
#' 
#' @param verbose If true outputs some progress messages
#' 
#' @param track_fit  \code{trace} to output that saves current values of all iterations
#' 
#' @return A susie fit, which is a list with some or all of the following elements\cr
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
#'
#' @references
#'
#' G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2018). A simple
#' new approach to variable selection in regression, with application
#' to genetic fine-mapping. \emph{bioRxiv}
#' \url{https://doi.org/10.1101/501114}.
#'
#' @examples
#' 
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' y = X %*% beta + rnorm(n)
#' res = susie(X,y,L=10)
#' coef(res)
#' plot(y,predict(res))
#'
#' @importFrom stats var
#' @importFrom utils modifyList
#'
#' @export
#'
susie <- function(X,Y,L = min(10,ncol(X)),scaled_prior_variance = 0.2,
                 residual_variance=NULL,
                 prior_weights=NULL, null_weight=NULL,
                 standardize=TRUE,intercept=TRUE,
                 estimate_residual_variance=TRUE,
                 estimate_prior_variance = FALSE,
                 optimV_method = 'optim',
                 s_init = NULL,coverage=0.95,min_abs_corr=0.5,
                 compute_univariate_zscore = FALSE,
                 max_iter=100,tol=1e-3,
                 verbose=FALSE,track_fit=FALSE) {
    
  # Check input X.
  if (!(is.double(X) & is.matrix(X)) & !inherits(X,"CsparseMatrix") & is.null(attr(X,"matrix.type")))
    stop("Input X must be a double-precision matrix, or a sparse matrix, or a trend filtering matrix.")
  if (is.numeric(null_weight) && null_weight == 0) null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight<0 || null_weight>=1)
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
  X = set_X_attributes(X,center=intercept, scale=standardize)
  # initialize susie fit
  s = init_setup(n,p,L,scaled_prior_variance,residual_variance,prior_weights,null_weight,as.numeric(var(Y)))
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
      tracking[[i]] = susie_slim(s)
    s = update_each_effect(X, Y, s, estimate_prior_variance,optimV_method)
    if(verbose){
        print(paste0("objective:",get_objective(X,Y,s)))
    }
    if(estimate_residual_variance){
      s$sigma2 = estimate_residual_variance(X,Y,s)
      if(verbose){
        print(paste0("objective:",get_objective(X,Y,s)))
      }
    }
    #s = remove_null_effects(s)

    elbo[i+1] = get_objective(X,Y,s)
    if((elbo[i+1]-elbo[i])<tol) {
      s$converged = TRUE
      break;
    }
  }
  elbo = elbo[1:(i+1)] #remove trailing NAs
  s$elbo <- elbo
  s$niter <- i

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in", max_iter, "iterations!"))
    s$converged = FALSE
  }

  if(intercept){
    s$intercept = mean_y - sum(attr(X,"scaled:center")* (colSums(s$alpha*s$mu)/attr(X,"scaled:scale")))# estimate intercept (unshrunk)
    s$fitted = s$Xr + mean_y
  } else {
    s$intercept = 0
    s$fitted = s$Xr
  }

  if (track_fit)
    s$trace = tracking

  ## SuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s, coverage=coverage, X=X, min_abs_corr=min_abs_corr)
    s$pip = susie_get_pip(s, prune_by_cs = FALSE)
  }
  ## report z-scores from univariate regression
  if (compute_univariate_zscore) {
    if (!is.null(null_weight) && null_weight != 0) {
      X = X[,1:(ncol(X)-1)]
    }
    s$z = calc_z(X,Y,centered=intercept)
  }
  ## for prediction
  s$X_column_scale_factors = attr(X,"scaled:scale")
  return(s)
}
