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
#' @param residual_variance_lowerbound Lower limit on the estimated
#'   residual variance. It is only relevant when
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @param refine If \code{refine = TRUE}, we use a procedure to help
#'   SuSiE get out of local optimum.
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
#' \item{lbf}{log-Bayes Factor for each single effect.}
#'
#' \item{lbf_variable}{log-Bayes Factor for each variable and single effect.}
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
#' @references
#'
#' G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple
#'   new approach to variable selection in regression, with application
#'   to genetic fine-mapping. \emph{Journal of the Royal Statistical
#'   Society, Series B} \url{https://doi.org/10.1101/501114}.
#'
#' @seealso \code{\link{susie_rss}}
#' 
#' @examples
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
                   track_fit = FALSE,
                   residual_variance_lowerbound = var(drop(Y))/1e4,
                   refine = FALSE) {

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
  if (!missing(s_init) && !is.null(s_init)) {
    if (!inherits(s_init,"susie"))
      stop("s_init should be a susie object")
    if (max(s_init$alpha) > 1 || min(s_init$alpha) < 0)
      stop("s_init$alpha has invalid values outside range [0,1]; please ",
           "check your input")
    # First, remove effects with s_init$V = 0
    s_init = susie_prune_single_effects(s_init, verbose=FALSE)
    # Then prune or expand
    s_init = susie_prune_single_effects(s_init, L, s$V, verbose)
    s = modifyList(s,s_init)
    s = init_finalize(s,X = X)
  } else {
    s = init_finalize(s)
  }
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
      s$sigma2 = pmax(residual_variance_lowerbound,
                      estimate_residual_variance(X,Y,s))
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
  names(s$fitted) = `if`(is.null(names(Y)), rownames(X), names(Y))

  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s,coverage = coverage,X = X,
                          min_abs_corr = min_abs_corr)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }

  if (!is.null(colnames(X))) {
    variable_names = colnames(X)
    if (!is.null(null_weight))
      variable_names = c("null", variable_names)
    colnames(s$alpha) = variable_names
    colnames(s$mu) = variable_names
    colnames(s$mu2) = variable_names
    colnames(s$lbf_variable) = variable_names
    names(s$pip) = variable_names
  }
  # report z-scores from univariate regression.
  if (compute_univariate_zscore) {
    if (!is.null(null_weight) && null_weight != 0)
      X = X[,1:(ncol(X) - 1)]
    s$z = calc_z(X,Y,center = intercept,scale = standardize)
  }

  # For prediction.
  s$X_column_scale_factors = attr(X,"scaled:scale")

  if(refine){
    if(!is.null(null_weight) && null_weight!=0 && !compute_univariate_zscore){
      ## if null_weight is specified, and the extra 0 column is not removed from compute_univariate_zscore,
      ## we remove it here
      X = X[,1:(ncol(X) - 1)]
    }
    conti = TRUE
    while(conti){
      m = list()
      for(cs in 1:length(s$sets$cs)){
        if(!missing(s_init) && !is.null(s_init)){
          warning('The given s_init is not used in refinement.')
        }
        pw = rep(1, ncol(X))
        pw[s$sets$cs[[cs]]] = 0
        s2 = susie(X,Y,L = L,
                   scaled_prior_variance = scaled_prior_variance,residual_variance = residual_variance,
                   prior_weights = pw,s_init = NULL,
                   null_weight = null_weight,standardize = standardize,intercept = intercept,
                   estimate_residual_variance = estimate_residual_variance,
                   estimate_prior_variance = estimate_prior_variance,
                   estimate_prior_method = estimate_prior_method,
                   check_null_threshold = check_null_threshold,
                   prior_tol = prior_tol,coverage = coverage,
                   residual_variance_upperbound = residual_variance_upperbound,
                   min_abs_corr = min_abs_corr,compute_univariate_zscore = FALSE,
                   na.rm = na.rm,max_iter = max_iter,tol = tol,
                   verbose = FALSE,track_fit = FALSE,residual_variance_lowerbound = var(drop(Y))/1e4,
                   refine = FALSE)
        sinit2 = s2[c('alpha', 'mu', 'mu2')]
        class(sinit2) = 'susie'
        s3 = susie(X,Y,L = L,
                   scaled_prior_variance = scaled_prior_variance,residual_variance = residual_variance,
                   prior_weights = NULL, s_init = sinit2,
                   null_weight = null_weight,standardize = standardize,intercept = intercept,
                   estimate_residual_variance = estimate_residual_variance,
                   estimate_prior_variance = estimate_prior_variance,
                   estimate_prior_method = estimate_prior_method,
                   check_null_threshold = check_null_threshold,
                   prior_tol = prior_tol,coverage = coverage,
                   residual_variance_upperbound = residual_variance_upperbound,
                   min_abs_corr = min_abs_corr,compute_univariate_zscore = FALSE,
                   na.rm = na.rm,max_iter = max_iter,tol = tol,
                   verbose = FALSE,track_fit = FALSE,residual_variance_lowerbound = var(drop(Y))/1e4,
                   refine = FALSE)
        m = c(m, list(s3))
      }
      elbo = sapply(m, function(x) susie_get_objective(x))
      if((max(elbo) - susie_get_objective(s)) <= 0){
        conti=FALSE
      }else{
        s = m[[which.max(elbo)]]
      }
    }
  }
  return(s)
}
