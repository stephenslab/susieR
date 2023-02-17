#' @rdname susie
#'
#' @title Sum of Single Effects (SuSiE) Regression
#'
#' @description Performs a sparse Bayesian multiple linear regression
#'   of y on X, using the "Sum of Single Effects" model from Wang et al
#'   (2020). In brief, this function fits the regression model \eqn{y =
#'   \mu + X b + e}, where elements of \eqn{e} are \emph{i.i.d.} normal
#'   with zero mean and variance \code{residual_variance}, \eqn{\mu} is
#'   an intercept term and \eqn{b} is a vector of length p representing
#'   the effects to be estimated. The \dQuote{susie assumption} is that
#'   \eqn{b = \sum_{l=1}^L b_l} where each \eqn{b_l} is a vector of
#'   length p with exactly one non-zero element. The prior on the
#'   non-zero element is normal with zero mean and variance \code{var(y)
#'   * scaled_prior_variance}. The value of \code{L} is fixed, and
#'   should be chosen to provide a reasonable upper bound on the number
#'   of non-zero effects to be detected. Typically, the hyperparameters
#'   \code{residual_variance} and \code{scaled_prior_variance} will be
#'   estimated during model fitting, although they can also be fixed as
#'   specified by the user. See functions \code{\link{susie_get_cs}} and
#'   other functions of form \code{susie_get_*} to extract the most
#'   commonly-used results from a susie fit.
#'
#' @details The function \code{susie} implements the IBSS algorithm
#' from Wang et al (2020). The option \code{refine = TRUE} implements
#' an additional step to help reduce problems caused by convergence of
#' the IBSS algorithm to poor local optima (which is rare in our
#' experience, but can provide misleading results when it occurs). The
#' refinement step incurs additional computational expense that
#' increases with the number of CSs found in the initial run.
#'
#' The function \code{susie_suff_stat} implements essentially the same
#' algorithms, but using sufficient statistics. (The statistics are
#' sufficient for the regression coefficients \eqn{b}, but not for the
#' intercept \eqn{\mu}; see below for how the intercept is treated.)
#' If the sufficient statistics are computed correctly then the
#' results from \code{susie_suff_stat} should be the same as (or very
#' similar to) \code{susie}, although runtimes will differ as
#' discussed below. The sufficient statistics are the sample
#' size \code{n}, and then the p by p matrix \eqn{X'X}, the p-vector
#' \eqn{X'y}, and the sum of squared y values \eqn{y'y}, all computed
#' after centering the columns of \eqn{X} and the vector \eqn{y} to
#' have mean 0; these can be computed using \code{compute_suff_stat}.
#'
#' The handling of the intercept term in \code{susie_suff_stat} needs
#' some additional explanation. Computing the summary data after
#' centering \code{X} and \code{y} effectively ensures that the
#' resulting posterior quantities for \eqn{b} allow for an intercept
#' in the model; however, the actual value of the intercept cannot be
#' estimated from these centered data. To estimate the intercept term
#' the user must also provide the column means of \eqn{X} and the mean
#' of \eqn{y} (\code{X_colmeans} and \code{y_mean}). If these are not
#' provided, they are treated as \code{NA}, which results in the
#' intercept being \code{NA}. If for some reason you prefer to have
#' the intercept be 0 instead of \code{NA} then set
#' \code{X_colmeans = 0,y_mean = 0}.
#'
#' For completeness, we note that if \code{susie_suff_stat} is run on
#' \eqn{X'X, X'y, y'y} computed \emph{without} centering \eqn{X} and
#' \eqn{y}, and with \code{X_colmeans = 0,y_mean = 0}, this is
#' equivalent to \code{susie} applied to \eqn{X, y} with
#' \code{intercept = FALSE} (although results may differ due to
#' different initializations of \code{residual_variance} and
#' \code{scaled_prior_variance}). However, this usage is not
#' recommended for for most situations.
#'
#' The computational complexity of \code{susie} is \eqn{O(npL)} per
#' iteration, whereas \code{susie_suff_stat} is \eqn{O(p^2L)} per
#' iteration (not including the cost of computing the sufficient
#' statistics, which is dominated by the \eqn{O(np^2)} cost of
#' computing \eqn{X'X}). Because of the cost of computing \eqn{X'X},
#' \code{susie} will usually be faster. However, if \eqn{n >> p},
#' and/or if \eqn{X'X} is already computed, then
#' \code{susie_suff_stat} may be faster.
#'
#' @param X An n by p matrix of covariates.
#'
#' @param y The observed responses, a vector of length n.
#'
#' @param L Maximum number of non-zero effects in the susie
#'   regression model. If L is larger than the number of covariates, p,
#'   L is set to p.
#'
#' @param scaled_prior_variance The prior variance, divided by
#'   \code{var(y)} (or by \code{(1/(n-1))yty} for
#'   \code{susie_suff_stat}); that is, the prior variance of each
#'   non-zero element of b is \code{var(y) * scaled_prior_variance}. The
#'   value provided should be either a scalar or a vector of length
#'   \code{L}. If \code{estimate_prior_variance = TRUE}, this provides
#'   initial estimates of the prior variances.
#'
#' @param residual_variance Variance of the residual. If
#'   \code{estimate_residual_variance = TRUE}, this value provides the
#'   initial estimate of the residual variance. By default, it is set to
#'   \code{var(y)} in \code{susie} and \code{(1/(n-1))yty} in
#'   \code{susie_suff_stat}.
#'
#' @param prior_weights A vector of length p, in which each entry
#'   gives the prior probability that corresponding column of X has a
#'   nonzero effect on the outcome, y.
#'
#' @param null_weight Prior probability of no effect (a number between
#'   0 and 1, and cannot be exactly 1).
#'
#' @param standardize If \code{standardize = TRUE}, standardize the
#'   columns of X to unit variance prior to fitting (or equivalently
#'   standardize XtX and Xty to have the same effect). Note that
#'   \code{scaled_prior_variance} specifies the prior on the
#'   coefficients of X \emph{after} standardization (if it is
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
#' @param na.rm Drop any missing values in y from both X and y.
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
#' @param refine If \code{refine = TRUE}, then an additional
#'  iterative refinement procedure is used, after the IBSS algorithm,
#'  to check and escape from local optima (see details).
#'
#' @param n_purity Passed as argument \code{n_purity} to
#'   \code{\link{susie_get_cs}}.
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
#'   \code{scaled_prior_variance * var(y)}.}
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
#' G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple
#' new approach to variable selection in regression, with application
#' to genetic fine-mapping. \emph{Journal of the Royal Statistical
#' Society, Series B} \bold{82}, 1273-1300 \doi{10.1101/501114}.
#'
#' Y. Zou, P. Carbonetto, G. Wang, G and M. Stephens
#' (2022). Fine-mapping from summary data with the \dQuote{Sum of
#' Single Effects} model. \emph{PLoS Genetics} \bold{18},
#' e1010299. \doi{10.1371/journal.pgen.1010299}.
#'
#' @seealso \code{\link{susie_get_cs}} and other \code{susie_get_*}
#' functions for extracting results; \code{\link{susie_trendfilter}} for
#' applying the SuSiE model to non-parametric regression, particularly
#' changepoint problems, and \code{\link{susie_rss}} for applying the
#' SuSiE model when one only has access to limited summary statistics
#' related to \eqn{X} and \eqn{y} (typically in genetic applications).
#'
#' @examples
#' # susie example
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#' res1 = susie(X,y,L = 10)
#' susie_get_cs(res1) # extract credible sets from fit
#' plot(beta,coef(res1)[-1])
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#' plot(y,predict(res1))
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#'
#' # susie_suff_stat example
#' input_ss = compute_suff_stat(X,y)
#' res2 = with(input_ss,
#'             susie_suff_stat(XtX = XtX,Xty = Xty,yty = yty,n = n,
#'                             X_colmeans = X_colmeans,y_mean = y_mean,L = 10))
#' plot(coef(res1),coef(res2))
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#'
#' @importFrom stats var
#' @importFrom utils modifyList
#'
#' @export
#'
susie = function (X,y,L = min(10,ncol(X)),
                   scaled_prior_variance = 0.2,
                   residual_variance = NULL,
                   prior_weights = NULL,
                   null_weight = 0,
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
                   residual_variance_lowerbound = var(drop(y))/1e4,
                   refine = FALSE,
                   n_purity = 100) {

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
  if (anyNA(X))
    stop("Input X must not contain missing values")
  if (anyNA(y)) {
    if (na.rm) {
      samples_kept = which(!is.na(y))
      y = y[samples_kept]
      X = X[samples_kept,]
    } else
      stop("Input y must not contain missing values")
  }
  p = ncol(X)
  if (p > 1000 & !requireNamespace("Rfast",quietly = TRUE))
    warning_message("For an X with many columns, please consider installing",
                    "the Rfast package for more efficient credible set (CS)",
                    "calculations.", style='hint')

  # Check input y.
  n = nrow(X)
  mean_y = mean(y)

  # Center and scale input.
  if (intercept)
    y = y - mean_y

  # Set three attributes for matrix X: attr(X,'scaled:center') is a
  # p-vector of column means of X if center=TRUE, a p vector of zeros
  # otherwise; 'attr(X,'scaled:scale') is a p-vector of column
  # standard deviations of X if scale=TRUE, a p vector of ones
  # otherwise; 'attr(X,'d') is a p-vector of column sums of
  # X.standardized^2,' where X.standardized is the matrix X centered
  # by attr(X,'scaled:center') and scaled by attr(X,'scaled:scale').
  out = compute_colstats(X,center = intercept,scale = standardize)
  attr(X,"scaled:center") = out$cm
  attr(X,"scaled:scale") = out$csd
  attr(X,"d") = out$d

  # Initialize susie fit.
  s = init_setup(n,p,L,scaled_prior_variance,residual_variance,prior_weights,
                 null_weight,as.numeric(var(y)),standardize)
  if (!missing(s_init) && !is.null(s_init)) {
    if (!inherits(s_init,"susie"))
      stop("s_init should be a susie object")
    if (max(s_init$alpha) > 1 || min(s_init$alpha) < 0)
      stop("s_init$alpha has invalid values outside range [0,1]; please ",
           "check your input")

    # First, remove effects with s_init$V = 0
    s_init = susie_prune_single_effects(s_init)
    num_effects = nrow(s_init$alpha)
    if(missing(L)){
      L = num_effects
    }else if(min(p,L) < num_effects){
      warning_message(paste("Specified number of effects L =",min(p,L),
                    "is smaller than the number of effects",num_effects,
                    "in input SuSiE model. The SuSiE model will have",
                    num_effects,"effects."))
      L = num_effects
    }
    # expand s_init if L > num_effects.
    s_init = susie_prune_single_effects(s_init, min(p, L), s$V)
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
    s = update_each_effect(X,y,s,estimate_prior_variance,estimate_prior_method,
                           check_null_threshold)
    if (verbose)
      print(paste0("objective:",get_objective(X,y,s)))

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance before the update.
    elbo[i+1] = get_objective(X,y,s)
    if ((elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }

    if (estimate_residual_variance) {
      s$sigma2 = pmax(residual_variance_lowerbound,
                      estimate_residual_variance(X,y,s))
      if (s$sigma2 > residual_variance_upperbound)
        s$sigma2 = residual_variance_upperbound
      if (verbose)
        print(paste0("objective:",get_objective(X,y,s)))
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
  names(s$fitted) = `if`(is.null(names(y)),rownames(X),names(y))

  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s,coverage = coverage,X = X,
                          min_abs_corr = min_abs_corr,
                          n_purity = n_purity)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }

  if (!is.null(colnames(X))) {
    variable_names = colnames(X)
    if (!is.null(null_weight)) {
      variable_names[length(variable_names)] = "null"
      names(s$pip) = variable_names[-p]
    } else
      names(s$pip)    = variable_names
    colnames(s$alpha) = variable_names
    colnames(s$mu)    = variable_names
    colnames(s$mu2)   = variable_names
    colnames(s$lbf_variable) = variable_names
  }
  # report z-scores from univariate regression.
  if (compute_univariate_zscore) {
    if (!is.matrix(X))
      warning("Calculation of univariate regression z-scores is not ",
              "implemented specifically for sparse or trend filtering ",
              "matrices, so this step may be slow if the matrix is large; ",
              "to skip this step set compute_univariate_zscore = FALSE")
    if (!is.null(null_weight) && null_weight != 0)
      X = X[,1:(ncol(X) - 1)]
    s$z = calc_z(X,y,center = intercept,scale = standardize)
  }

  # For prediction.
  s$X_column_scale_factors = attr(X,"scaled:scale")

  if (refine) {
    if (!missing(s_init) && !is.null(s_init))
      warning("The given s_init is not used in refinement")
    if (!is.null(null_weight) && null_weight != 0) {

      ## if null_weight is specified, we compute the original prior_weight
      pw_s = s$pi[-s$null_index]/(1-null_weight)
      if (!compute_univariate_zscore)

        ## if null_weight is specified, and the extra 0 column is not
        ## removed from compute_univariate_zscore, we remove it here
        X = X[,1:(ncol(X) - 1)]
    } else
      pw_s = s$pi
    conti = TRUE
    while (conti & length(s$sets$cs)>0) {
      m = list()
      for(cs in 1:length(s$sets$cs)){
        pw_cs = pw_s
        pw_cs[s$sets$cs[[cs]]] = 0
        if(all(pw_cs == 0)){
          break
        }
        s2 = susie(X,y,L = L,scaled_prior_variance = scaled_prior_variance,
            residual_variance = residual_variance,
            prior_weights = pw_cs, s_init = NULL,null_weight = null_weight,
            standardize = standardize,intercept = intercept,
            estimate_residual_variance = estimate_residual_variance,
            estimate_prior_variance = estimate_prior_variance,
            estimate_prior_method = estimate_prior_method,
            check_null_threshold = check_null_threshold,
            prior_tol = prior_tol,coverage = coverage,
            residual_variance_upperbound = residual_variance_upperbound,
            min_abs_corr = min_abs_corr,
            compute_univariate_zscore = compute_univariate_zscore,
            na.rm = na.rm,max_iter = max_iter,tol = tol,verbose = FALSE,
            track_fit = FALSE,residual_variance_lowerbound = var(drop(y))/1e4,
            refine = FALSE)
        sinit2 = s2[c("alpha","mu","mu2")]
        class(sinit2) = "susie"
        s3 = susie(X,y,L = L,scaled_prior_variance = scaled_prior_variance,
            residual_variance = residual_variance,prior_weights = pw_s,
            s_init = sinit2,null_weight = null_weight,
            standardize = standardize,intercept = intercept,
            estimate_residual_variance = estimate_residual_variance,
            estimate_prior_variance = estimate_prior_variance,
            estimate_prior_method = estimate_prior_method,
            check_null_threshold = check_null_threshold,
            prior_tol = prior_tol,coverage = coverage,
            residual_variance_upperbound = residual_variance_upperbound,
            min_abs_corr = min_abs_corr,
            compute_univariate_zscore = compute_univariate_zscore,
            na.rm = na.rm,max_iter = max_iter,tol = tol,verbose = FALSE,
            track_fit = FALSE,residual_variance_lowerbound = var(drop(y))/1e4,
            refine = FALSE)
        m = c(m,list(s3))
      }
      if(length(m) == 0){
        conti = FALSE
      }else{
        elbo = sapply(m,function(x) susie_get_objective(x))
        if ((max(elbo) - susie_get_objective(s)) <= 0)
          conti = FALSE
        else
          s = m[[which.max(elbo)]]
      }
    }
  }
  return(s)
}
