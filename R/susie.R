#' @title SUm of Single Effects (SuSiE) Regression
#'
#' @description Performs Bayesian multiple linear regression of Y on
#'   X; that is, this function fits the regression model \eqn{Y = sum_l
#'   Xb_l + e}, where elements of e are \emph{i.i.d.} normal with zero
#'   mean and variance \code{residual_variance}, and \eqn{sum_l b_l} is
#'   a vector of length p representing the effects to be estimated. The
#'   SuSiE assumption is that each b_l has exactly one non-zero
#'   element. The prior on the non-zero element is normal with zero mean
#'   and variance \code{var(Y)*scaled_prior_variance}.
#'
#'   The model is fitted using the "Iterative Bayesian Stepwise
#'   Selection" (IBSS) algorithm.
#'
#'   See also \code{susie_trendfilter} for applying susie to non-parametric regression, particularly changepoint problems
#'
#' @param X An n by p matrix of covariates.
#'
#' @param Y A vector of length n.
#'
#' @param L Number of components (nonzero elements) in the SuSiE
#'   regression model. If \code{L} is larger than the number of
#'   covariate (p), \code{L} is set to p.
#'
#' @param scaled_prior_variance The scaled prior variance. This is
#'    either a scalar, or a vector of length \code{L}. The prior variance
#'   of each non-zero element of b is set to
#'   \code{var(Y)*scaled_prior_variance}. If
#'   \code{estimate_prior_variance = TRUE}, this input provides the
#'   initial estimates of the prior variances.
#'
#' @param residual_variance The variance of the residual. If
#'   \code{estimate_residual_variance = TRUE}, this value provides the
#'   initial estimate of the residual variance. By default, it is set to
#'   \code{var(Y)}.
#'
#' @param prior_weights A vector of length p, in which each entry
#'   gives the prior probability that corresponding column of X has a
#'   nonzero effect on the outcome, Y.
#'
#' @param null_weight Prior probability of no effect (a number between
#'   0 and 1, and cannot be exactly 1).
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
#'   residual is estimated separately for each of the \code{L}
#'   components, and \code{scaled_prior_variance} is used as an initial
#'   estimate of the variances.
#'
#' @param estimate_prior_variance If \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated father than fixed.
#'
#' @param estimate_prior_method The method used for estimating prior
#'   variance.
#'
#' @param s_init A previous susie fit with which to initialize.
#'
#' @param coverage A number between 0 and 1 specifying the coverage of
#'  the estimated confidence sets.
#'
#' @param min_abs_corr Minimum of absolute value of correlation
#'   allowed in a credible set. The default, 0.5, corresponds to squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetics studies.
#'
#' @param compute_univariate_zscore If \code{compute_univariate_zscore
#'   = TRUE}, the univariate regression z-scores are outputted for each
#'   variable.
#'
#' @param max_iter Maximum number of iterations of the IBSS fitting
#'   procedure.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   "ELBO" (this is the objective function to be maximized), is less
#'   than \code{tol}.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's
#'   progress and a summary of the optimization settings are printed to
#'   the console.
#'
#' @param track_fit If \code{track_fit = TRUE}, an object \code{trace}
#'   is also returned containing detailed information about the
#'   estimates at each iteration of the IBSS fitting procedure.
#'
#' @return A \code{"susie"} object with some or all of the following
#'   elements:
#'
#' \item{alpha}{An L by p matrix of posterior inclusion probabilites.}
#'
#' \item{mu}{An L by p matrix of posterior means, conditional on
#'   inclusion.}
#'
#' \item{mu2}{An L by p matrix of posterior second moments,
#'   conditional on inclusion.}
#'
#' \item{Xr}{An vector of length n, equal to \code{X \%*\% colSums(alpha
#'   * mu)}.}
#'
#' \item{intercept}{The intercept (fixed or estimated).}
#'
#' \item{sigma2}{Residual variance (fixed or estimated).}
#'
#' \item{V}{Prior variance of the non-zero elements of b, equal to
#'   \code{scaled_prior_variance * var(Y)}.}
#'
#' \item{elbo}{The value of the variational lower bound, or "ELBO"
#'   (the objective function to be maximized), achieved at each
#'   iteration of the IBSS fitting procedure.}
#'
#' \item{fitted}{Vector of length n containing the "fitted" values of
#'   the outcome.}
#'
#' \item{sets}{Credible sets estimated from model fit; see
#'   \code{\link{susie_get_cs}} for details.}
#'
#' \item{pip}{A vector of length p giving the (marginal) posterior
#'   inclusion probabilities for all p covariates.}
#'
#' \item{z}{A vector of univariate z-scores.}
#'
#' \item{niter}{Number of IBSS iterations that were run.}
#'
#' \item{converged}{\code{TRUE} or \code{FALSE}, indicating whether
#'   the IBSS converged to a solution within the chosen tolerance
#'   level.}
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
                 estimate_prior_variance = TRUE,
                 estimate_prior_method = c("optim","EM"),
                 s_init = NULL,coverage=0.95,min_abs_corr=0.5,
                 compute_univariate_zscore = FALSE,
                 max_iter=100,tol=1e-3,
                 verbose=FALSE,track_fit=FALSE) {

  # Process input estimate_prior_method.
  estimate_prior_method <- match.arg(estimate_prior_method)

  # Check input X.
  if (!(is.double(X) & is.matrix(X)) & !inherits(X,"CsparseMatrix") & is.null(attr(X,"matrix.type")))
    stop("Input X must be a double-precision matrix, or a sparse matrix, or a trend filtering matrix.")
  if (is.numeric(null_weight) && null_weight == 0) null_weight = NULL
  if (!is.null(null_weight) && is.null(attr(X, "matrix.type"))) {
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
  s = init_setup(n,p,L,scaled_prior_variance,residual_variance,prior_weights,null_weight,as.numeric(var(Y)),standardize)
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
    s = update_each_effect(X, Y, s, estimate_prior_variance,estimate_prior_method)
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
  elbo = elbo[2:(i+1)] # Remove first (infinite) entry, and trailing NAs.
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
  s$fitted <- drop(s$fitted)

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
