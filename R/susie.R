#' @rdname susie
#'
#' @title Sum of Single Effects (SuSiE) Regression
#'
#' @description Performs Bayesian multiple linear regression of Y on
#'   X; that is, this function fits the regression model \eqn{Y = \sum_l
#'   X b_{l=1}^L + e}, where elements of e are \emph{i.i.d.} normal with
#'   zero mean and variance \code{residual_variance}, and
#'   \eqn{\sum_{l=1}^L b_l} is a vector of length p representing the
#'   effects to be estimated. The \dQuote{susie assumption} is that each
#'   \eqn{b_l} has exactly one non-zero element. The prior on the
#'   non-zero element is normal with zero mean and variance \code{var(Y)
#'   * scaled_prior_variance}. The model is fitted using the
#'   \dQuote{Iterative Bayesian Stepwise Selection} (IBSS) algorithm.
#'   See also \code{\link{susie_trendfilter}} for applying susie to
#'   non-parametric regression, particularly changepoint problems.
#'
#' @details \code{susie_suff_stat} performs sum of single-effect
#' linear regression with summary statistics. The required summary
#' data are either: \code{bhat}, \code{shat}, the p by p symmetric,
#' positive semidefinite correlation (or covariance) matrix \code{R},
#' the sample size \code{n}, and the variance of y; or the p by p
#' matrix \eqn{X'X}, the p-vector \eqn{X'y}, the sum of squares
#' \eqn{y'y}, and the sample size \code{n}. The summary statistics
#' should come from the same individuals. Both the columns of X and
#' the vector y should be centered to have mean zero before computing
#' these summary statistics; you may also want to scale each column of
#' X and y to have variance 1 (see examples).
#' 
#' \code{susie_rss} performs sum of single-effect linear regression
#' with z scores; all posterior calculations are for z-scores. The
#' required summary data are the p by p correlation matrix, \code{R},
#' and the p-vector \code{z}. The summary stats should come from the
#' same individuals (samples).
#' 
#' susie_auto is an attempt to automate reliable running of susie even
#' on hard problems. It implements a three-stage strategy for each L:
#' first, fit susie with very small residual error; next, estimate
#' residual error; finally, estimate the prior variance. If the last
#' step estimates some prior variances to be zero, stop. Otherwise,
#' double L, and repeat. Initial runs are performed with relaxed
#' tolerance; the final run is performed using the default susie
#' tolerance.
#'
#' @param X An n by p matrix of covariates.
#'
#' @param Y The observed responses, a vector of length n.
#'
#' @param L Number of components (nonzero coefficients) in the susie
#'   regression model. If L is larger than the number of covariates, p,
#'   L is set to p.
#'
#' @param scaled_prior_variance The scaled prior variance. This is
#'   either a scalar or a vector of length \code{L}. The prior variance
#'   of each non-zero element of b is set to \code{var(Y) *
#'   scaled_prior_variance}. If \code{estimate_prior_variance = TRUE},
#'   this provides initial estimates of the prior variances.
#'
#' @param residual_variance Variance of the residual. If
#'   \code{estimate_residual_variance = TRUE}, this value provides the
#'   initial estimate of the residual variance. By default, it is
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
#'   columns of X (or XtX and Xty) to unit variance prior to
#'   fitting. Note that \code{scaled_prior_variance} specifies the prior
#'   on the coefficients of X \emph{after} standardization (if it is
#'   performed). If you do not standardize, you may need to think more
#'   carefully about specifying \code{scaled_prior_variance}. Whatever
#'   your choice, the coefficients returned by \code{coef} are given for
#'   \code{X} on the original input scale. Any column of \code{X} that
#'   has zero variance is not standardized.
#'
#' @param intercept If \code{intercept = TRUE}, the intercept is
#'   fitted; it \code{intercept = FALSE}, the intercept is set to
#'   zero. Setting \code{intercept = FALSE} is generally not
#'   recommended.
#'
#' @param estimate_residual_variance If
#'   \code{estimate_residual_variance = TRUE}, the residual variance is
#'   estimated, using \code{residual_variance} as an initial value. If
#'   \code{estimate_residual_variance = FALSE}, the residual variance is
#'   fixed to the value supplied by \code{residual_variance}.
#'
#' @param estimate_prior_variance If \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated (this is a separate
#'   parameter for each of the L effects). If provided,
#'   \code{scaled_prior_variance} is then used as an initial value for
#'   the optimization. When \code{estimate_prior_variance = FALSE}, the
#'   prior variance for each of the L effects is determined by the
#'   value supplied to \code{scaled_prior_variance}.
#'
#' @param estimate_prior_method The method used for estimating prior
#'   variance. When \code{estimate_prior_method = "simple"} is used, the
#'   likelihood at the specified prior variance is compared to the
#'   likelihood at a variance of zero, and the setting with the larger
#'   likelihood is retained.
#'
#' @param check_null_threshold When the prior variance is estimated,
#'   compare the estimate with the null, and set the prior variance to
#'   zero unless the log-likelihood using the estimate is larger by this
#'   threshold amount. For example, if you set
#'   \code{check_null_threshold = 0.1}, this will "nudge" the estimate
#'   towards zero when the difference in log-likelihoods is small. A
#'   note of caution that setting this to a value greater than zero may
#'   lead the IBSS fitting procedure to occasionally decrease the ELBO.
#'
#' @param prior_tol When the prior variance is estimated, compare the
#'   estimated value to \code{prior_tol} at the end of the computation,
#'   and exclude a single effect from PIP computation if the estimated
#'   prior variance is smaller than this tolerance value.
#'
#' @param residual_variance_upperbound Upper limit on the estimated
#'   residual variance. It is only relevant when
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @param s_init A previous susie fit with which to initialize.
#'
#' @param coverage A number between 0 and 1 specifying the
#'   \dQuote{coverage} of the estimated confidence sets.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param compute_univariate_zscore If \code{compute_univariate_zscore
#'   = TRUE}, the univariate regression z-scores are outputted for each
#'   variable.
#'
#' @param na.rm Drop any missing values in Y from both X and Y.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#'   and a summary of the optimization settings, are printed to the
#'   console.
#'
#' @param track_fit If \code{track_fit = TRUE}, \code{trace}
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
#' \item{Xr}{A vector of length n, equal to \code{X \%*\% colSums(alpha
#'   * mu)}.}
#'
#' \item{intercept}{Intercept (fixed or estimated).}
#'
#' \item{sigma2}{Residual variance (fixed or estimated).}
#'
#' \item{V}{Prior variance of the non-zero elements of b, equal to
#'   \code{scaled_prior_variance * var(Y)}.}
#'
#' \item{elbo}{The value of the variational lower bound, or
#'   \dQuote{ELBO} (objective function to be maximized), achieved at
#'   each iteration of the IBSS fitting procedure.}
#'
#' \item{fitted}{Vector of length n containing the fitted values of
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
#' \item{niter}{Number of IBSS iterations that were performed.}
#'
#' \item{converged}{\code{TRUE} or \code{FALSE} indicating whether
#'   the IBSS converged to a solution within the chosen tolerance
#'   level.}
#'
#' \code{susie_suff_stat} returns also outputs:
#'
#' \item{XtXr}{A p-vector of \code{t(X)} times the fitted values,
#'   \code{X \%*\% colSums(alpha*mu)}.}
#' 
#' \code{susie_rss} also outputs:
#' 
#' \item{Rr}{An p-vector of \code{t(X)} times fitted values, \code{X
#'   \%*\% colSums(alpha*mu)}.}
#'
#' @references
#'
#' G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple
#'   new approach to variable selection in regression, with application
#'   to genetic fine-mapping. \emph{bioRxiv}
#'   \url{https://doi.org/10.1101/501114}.
#'
#' @examples
#'
#' # susie example.
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#' res1 = susie(X,y,L = 10)
#' plot(beta,coef(res1)[-1])
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#' plot(y,predict(res1))
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#'
#' # susie_suff_stat example.
#' input_ss = compute_ss(X,y,standardize = TRUE)
#' res2 = with(input_ss,
#'             susie_suff_stat(XtX = XtX,Xty = Xty,yty = yty,n = n,L = 10))
#' plot(coef(res1)[-1],coef(res2)[-1])
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#'
#' # susie_rss example.
#' ss   <- susieR:::univariate_regression(X,y)
#' R    <- with(input_ss,cov2cor(XtX))
#' zhat <- with(ss,betahat/sebetahat)
#' res3 <- susie_rss(zhat,R,L = 10)
#' plot(coef(res1)[-1]/ss$sebetahat,coef(res3)[-1])
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#'
#' # susie_auto example.
#' 
#' @importFrom stats var
#' @importFrom utils modifyList
#'
#' @export
#'
susie <- function (X,Y,L = min(10,ncol(X)),
                   scaled_prior_variance = 0.2,
                   residual_variance = NULL,
                   prior_weights = NULL,
                   null_weight = NULL,
                   standardize = TRUE,
                   intercept = TRUE,
                   estimate_residual_variance = TRUE,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = c("optim", "EM", "simple"),
                   check_null_threshold = 0,
                   prior_tol = 1e-9,
                   residual_variance_upperbound = Inf,
                   s_init = NULL,
                   coverage = 0.95,
                   min_abs_corr = 0.5,
                   compute_univariate_zscore = FALSE,
                   na.rm = FALSE,
                   max_iter = 100,
                   tol = 1e-3,
                   verbose = FALSE,
                   track_fit = FALSE) {

  # Process input estimate_prior_method.
  estimate_prior_method = match.arg(estimate_prior_method)

  # Check input X.
  if (!(is.double(X) & is.matrix(X)) &
      !inherits(X,"CsparseMatrix") &
      is.null(attr(X,"matrix.type")))
    stop("Input X must be a double-precision matrix, or a sparse matrix, or ",
         "a trend filtering matrix")
  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight) && is.null(attr(X,"matrix.type"))) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(X) * (1 - null_weight),ncol(X)),null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight),null_weight)
    X = cbind(X,0)
  }
  if (any(is.na(X)))
    stop("Input X must not contain missing values")
  if (any(is.na(Y))) {
    if (na.rm) {
      samples_kept = which(!is.na(Y))
      Y = Y[samples_kept]
      X = X[samples_kept,]
    } else 
      stop("Input Y must not contain missing values")
  }
  
  # Check input Y.
  p = ncol(X)
  n = nrow(X)
  mean_y = mean(Y)

  # Center and scale input.
  if (intercept)
    Y = Y - mean_y
  X = set_X_attributes(X,center = intercept,scale = standardize)
  
  # Initialize susie fit.
  s = init_setup(n,p,L,scaled_prior_variance,residual_variance,prior_weights,
                 null_weight,as.numeric(var(Y)),standardize)
  if (!missing(s_init)) {
    s = modifyList(s,s_init)
    s = init_finalize(s,X = X)
  } else 
    s = init_finalize(s)

  # Initialize elbo to NA.
  elbo = rep(as.numeric(NA),max_iter + 1)
  elbo[1] = -Inf;
  tracking = list()

  for (i in 1:max_iter) {
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    s = update_each_effect(X,Y,s,estimate_prior_variance,estimate_prior_method,
                           check_null_threshold)
    if (verbose)
      print(paste0("objective:",get_objective(X,Y,s)))
    
    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance before the update.
    elbo[i+1] = get_objective(X,Y,s)
    if ((elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }
    if (estimate_residual_variance) {
      s$sigma2 = estimate_residual_variance(X,Y,s)
      if (s$sigma2 > residual_variance_upperbound)
        s$sigma2 = residual_variance_upperbound
      if (verbose)
        print(paste0("objective:",get_objective(X,Y,s)))
    }
  }

  # Remove first (infinite) entry, and trailing NAs.
  elbo = elbo[2:(i+1)] 
  s$elbo = elbo
  s$niter = i

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }

  if (intercept) {

    # Estimate unshrunk intercept.
    s$intercept = mean_y - sum(attr(X,"scaled:center") *
      (colSums(s$alpha * s$mu)/attr(X,"scaled:scale"))) 
    s$fitted = s$Xr + mean_y
  } else {
    s$intercept = 0
    s$fitted = s$Xr
  }
  s$fitted = drop(s$fitted)

  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s,coverage = coverage,X = X,
                          min_abs_corr = min_abs_corr)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }
  
  # report z-scores from univariate regression.
  if (compute_univariate_zscore) {
    if (!is.null(null_weight) && null_weight != 0) 
      X = X[,1:(ncol(X) - 1)]
    s$z = calc_z(X,Y,center = intercept,scale = standardize)
  }
  
  # For prediction.
  s$X_column_scale_factors = attr(X,"scaled:scale")
  return(s)
}
