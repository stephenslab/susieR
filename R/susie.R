# =============================================================================
# SuSiE WITH INDIVIDUAL-LEVEL DATA
# =============================================================================

#' @title Sum of Single Effects (SuSiE) Regression
#'
#' @description Performs a sparse Bayesian multiple linear regression
#' of y on X, using the "Sum of Single Effects" model from Wang et al
#' (2020). In brief, this function fits the regression model \eqn{y =
#' \mu + X b + e}, where elements of \eqn{e} are \emph{i.i.d.} normal
#' with zero mean and variance \code{residual_variance}, \eqn{\mu} is
#' an intercept term and \eqn{b} is a vector of length p representing
#' the effects to be estimated. The \dQuote{susie assumption} is that
#' \eqn{b = \sum_{l=1}^L b_l} where each \eqn{b_l} is a vector of
#' length p with exactly one non-zero element. The prior on the
#' non-zero element is normal with zero mean and variance \code{var(y)
#' * scaled_prior_variance}. The value of \code{L} is fixed, and
#' should be chosen to provide a reasonable upper bound on the number
#' of non-zero effects to be detected. Typically, the hyperparameters
#' \code{residual_variance} and \code{scaled_prior_variance} will be
#' estimated during model fitting, although they can also be fixed as
#' specified by the user. See functions \code{\link{susie_get_cs}} and
#' other functions of form \code{susie_get_*} to extract the most
#' commonly-used results from a susie fit.
#'
#' #' @details The function \code{susie} implements the IBSS algorithm
#' from Wang et al (2020). The option \code{refine = TRUE} implements
#' an additional step to help reduce problems caused by convergence of
#' the IBSS algorithm to poor local optima (which is rare in our
#' experience, but can provide misleading results when it occurs). The
#' refinement step incurs additional computational expense that
#' increases with the number of CSs found in the initial run.
#'
#' The function \code{susie_ss} implements essentially the same
#' algorithms, but using sufficient statistics. (The statistics are
#' sufficient for the regression coefficients \eqn{b}, but not for the
#' intercept \eqn{\mu}; see below for how the intercept is treated.)
#' If the sufficient statistics are computed correctly then the
#' results from \code{susie_ss} should be the same as (or very
#' similar to) \code{susie}, although runtimes will differ as
#' discussed below. The sufficient statistics are the sample
#' size \code{n}, and then the p by p matrix \eqn{X'X}, the p-vector
#' \eqn{X'y}, and the sum of squared y values \eqn{y'y}, all computed
#' after centering the columns of \eqn{X} and the vector \eqn{y} to
#' have mean 0; these can be computed using \code{compute_suff_stat}.
#'
#' The handling of the intercept term in \code{susie_ss} needs
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
#' For completeness, we note that if \code{susie_ss} is run on
#' \eqn{X'X, X'y, y'y} computed \emph{without} centering \eqn{X} and
#' \eqn{y}, and with \code{X_colmeans = 0,y_mean = 0}, this is
#' equivalent to \code{susie} applied to \eqn{X, y} with
#' \code{intercept = FALSE} (although results may differ due to
#' different initializations of \code{residual_variance} and
#' \code{scaled_prior_variance}). However, this usage is not
#' recommended for for most situations.
#'
#' The computational complexity of \code{susie} is \eqn{O(npL)} per
#' iteration, whereas \code{susie_ss} is \eqn{O(p^2L)} per
#' iteration (not including the cost of computing the sufficient
#' statistics, which is dominated by the \eqn{O(np^2)} cost of
#' computing \eqn{X'X}). Because of the cost of computing \eqn{X'X},
#' \code{susie} will usually be faster. However, if \eqn{n >> p},
#' and/or if \eqn{X'X} is already computed, then
#' \code{susie_ss} may be faster.
#'
#' @param X An n by p matrix of covariates.
#'
#' @param y The observed responses, a vector of length n.
#'
#' @param L Maximum number of non-zero effects in the model. If L is larger than
#' the number of covariates, p, L is set to p.
#'
#' @param scaled_prior_variance The prior variance, divided by
#'   \code{var(y)} (or by \code{(1/(n-1))yty} for
#'   \code{susie_ss}); that is, the prior variance of each
#'   non-zero element of b is \code{var(y) * scaled_prior_variance}. The
#'   value provided should be either a scalar or a vector of length
#'   \code{L}. If \code{estimate_prior_variance = TRUE}, this provides
#'   initial estimates of the prior variances.
#'
#' @param residual_variance Variance of the residual. If
#'   \code{estimate_residual_variance = TRUE}, this value provides the
#'   initial estimate of the residual variance. By default, it is set to
#'   \code{var(y)} in \code{susie} and \code{(1/(n-1))yty} in
#'   \code{susie_ss}.
#'
#' @param prior_weights A vector of length p, in which each entry
#'   gives the prior probability that corresponding column of X has a
#'   nonzero effect on the outcome, y.
#'
#' @param null_weight Prior probability of no effect (a number between 0 and 1,
#' and cannot be exactly 1).
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
#' @param estimate_residual_method The method used for estimating residual variance.
#'   For the original SuSiE model, "MLE" and "MoM" estimation is equivalent, but for
#'   the infinitesimal model, "MoM" is more stable. We recommend using "Servin_Stephens"
#'   when n < 80 for improved coverage, although it is currently only implemented
#'   for individual-level data.
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
#' @param unmappable_effects The method for modeling unmappable effects:
#'   "none", "inf", "ash".
#'
#' @param check_null_threshold When the prior variance is estimated,
#'   compare the estimate with the null, and set the prior variance to
#'   zero unless the log-likelihood using the estimate is larger by this
#'   threshold amount. For example, if you set
#'   \code{check_null_threshold = 0.1}, this will "nudge" the estimate
#'   towards zero when the difference in log-likelihoods is small. A
#'   note of caution that setting this to a value greater than zero may
#'   lead the IBSS fitting procedure to occasionally decrease the ELBO. This
#'   setting is disabled when using \code{unmappable_effects = "inf"} or
#'   \code{unmappable_effects = "ash"}.
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
#' @param model_init A previous susie fit with which to initialize.
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
#' @param tol tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure..
#'
#' @param convergence_method When \code{converge_method = "elbo"} the fitting
#'   procedure halts when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}. When \code{converge_method = "pip"} the fitting
#'   procedure halts when the maximum absolute difference in \code{alpha} is less
#'   than \code{tol}.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#'  a summary of the optimization settings, and refinement progress (if
#'  \code{refine = TRUE}) are printed to the console.
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
#' @param alpha0 Numerical parameter for the NIG prior when using
#' \code{estimate_residual_method = "Servin_Stephens}.
#'
#' @param beta0 Numerical parameter for the NIG prior when using
#' \code{estimate_residual_method = "Servin_Stephens}.
#'
#' @return A \code{"susie"} object with some or all of the following elements:
#'
#' \item{alpha}{An L by p matrix of posterior inclusion probabilities.}
#'
#' \item{mu}{An L by p matrix of posterior means, conditional on inclusion.}
#'
#' \item{mu2}{An L by p matrix of posterior second moments, conditional on
#'   inclusion.}
#'
#' \item{Xr}{A vector of length n, equal to \code{X \%*\% colSums(alpha * mu)}.}
#'
#' \item{lbf}{Log-Bayes Factor for each single effect.}
#'
#' \item{lbf_variable}{Log-Bayes Factor for each variable and single effect.}
#'
#' \item{intercept}{Intercept (fixed or estimated).}
#'
#' \item{sigma2}{Residual variance (fixed or estimated).}
#'
#' \item{V}{Prior variance of the non-zero elements of b.}
#'
#' \item{elbo}{The variational lower bound (or ELBO) achieved at each iteration.}
#'
#' \item{fitted}{Vector of length n containing the fitted values.}
#'
#' \item{sets}{Credible sets estimated from model fit.}
#'
#' \item{pip}{A vector of length p giving the marginal posterior inclusion
#'   probabilities.}
#'
#' \item{z}{A vector of univariate z-scores.}
#'
#' \item{niter}{Number of IBSS iterations performed.}
#'
#' \item{converged}{\code{TRUE} or \code{FALSE} indicating whether
#'   the IBSS converged to a solution within the chosen tolerance
#'   level.}
#'
#' \item{theta}{If \code{unmappable_effects = "inf"} or
#'   \code{unmappable_effects = "ash"}, then \code{theta} is a p-vector of posterior
#'   means for the unmappable effects.
#'
#' \item{tau2}{If \code{unmappable_effects = "inf"} or
#'   \code{unmappable_effects = "ash"}, then \code{tau2} is the unmappable variance.
#'
#' @export
susie <- function(X, y, L = min(10, ncol(X)),
                  scaled_prior_variance = 0.2,
                  residual_variance = NULL,
                  prior_weights = NULL,
                  null_weight = 0,
                  standardize = TRUE,
                  intercept = TRUE,
                  estimate_residual_variance = TRUE,
                  estimate_residual_method = c("MLE", "MoM", "Servin_Stephens"),
                  estimate_prior_variance = TRUE,
                  estimate_prior_method = c("optim", "EM", "simple"),
                  unmappable_effects = c("none", "inf", "ash"),
                  check_null_threshold = 0,
                  prior_tol = 1e-9,
                  residual_variance_upperbound = Inf,
                  model_init = NULL,
                  coverage = 0.95,
                  min_abs_corr = 0.5,
                  compute_univariate_zscore = FALSE,
                  na.rm = FALSE,
                  max_iter = 100,
                  tol = 1e-3,
                  convergence_method = c("elbo", "pip"),
                  verbose = FALSE,
                  track_fit = FALSE,
                  residual_variance_lowerbound = var(drop(y)) / 1e4,
                  refine = FALSE,
                  n_purity = 100,
                  alpha0 = 0,
                  beta0 = 0) {

  # Validate method arguments
  unmappable_effects       <- match.arg(unmappable_effects)
  estimate_prior_method    <- match.arg(estimate_prior_method)
  estimate_residual_method <- match.arg(estimate_residual_method)
  convergence_method       <- match.arg(convergence_method)

  # Construct data and params objects
  susie_objects <- individual_data_constructor(
    X, y, L, scaled_prior_variance, residual_variance,
    prior_weights, null_weight, standardize, intercept,
    estimate_residual_variance, estimate_residual_method,
    estimate_prior_variance, estimate_prior_method,
    unmappable_effects, check_null_threshold, prior_tol,
    residual_variance_upperbound, model_init, coverage,
    min_abs_corr, compute_univariate_zscore, na.rm,
    max_iter, tol, convergence_method, verbose, track_fit,
    residual_variance_lowerbound, refine, n_purity,
    alpha0, beta0
  )

  # Run main SuSiE algorithm
  model <- susie_workhorse(susie_objects$data, susie_objects$params)

  return(model)
}

# =============================================================================
# SuSiE WITH SUFFICIENT STATISTICS
# =============================================================================

#' @title SuSiE using Sufficient Statistics
#'
#' @description Performs SuSiE regression using sufficient statistics (XtX, Xty,
#' yty, n) instead of individual-level data (X, y).
#'
#' @param XtX A p by p matrix, X'X, with columns of X centered to have mean zero.
#'
#' @param Xty A p-vector, X'y, with y and columns of X centered to have mean zero.
#'
#' @param yty A scalar, y'y, with y centered to have mean zero.
#'
#' @param n The sample size.
#'
#' @param X_colmeans A p-vector of column means of \code{X}. If both
#'   \code{X_colmeans} and \code{y_mean} are provided, the intercept
#'   is estimated; otherwise, the intercept is NA.
#'
#' @param y_mean A scalar containing the mean of \code{y}. If both
#'   \code{X_colmeans} and \code{y_mean} are provided, the intercept
#'   is estimated; otherwise, the intercept is NA.
#'
#' @param maf A p-vector of minor allele frequencies; to be used along with
#'   \code{maf_thresh} to filter input summary statistics.
#'
#' @param maf_thresh Variants with MAF smaller than this threshold are not used.
#'
#' @param check_input If \code{check_input = TRUE}, \code{susie_ss} performs
#'   additional checks on \code{XtX} and \code{Xty}. The checks are:
#'   (1) check that \code{XtX} is positive semidefinite; (2) check that
#'   \code{Xty} is in the space spanned by the non-zero eigenvectors of \code{XtX}.
#'
#' @param r_tol Tolerance level for eigenvalue check of positive semidefinite
#'   matrix \code{XtX}.
#'
#' @param check_prior If \code{check_prior = TRUE}, it checks if the
#'   estimated prior variance becomes unreasonably large (comparing with
#'   10 * max(abs(z))^2).
#'
#' @export
susie_ss <- function(XtX, Xty, yty, n,
                     L = min(10, ncol(XtX)),
                     X_colmeans = NA, y_mean = NA,
                     maf = NULL, maf_thresh = 0,
                     check_input = FALSE,
                     r_tol = 1e-8,
                     standardize = TRUE,
                     scaled_prior_variance = 0.2,
                     residual_variance = NULL,
                     prior_weights = NULL,
                     null_weight = 0,
                     model_init = NULL,
                     estimate_residual_variance = TRUE,
                     estimate_residual_method = c("MLE", "MoM", "Servin_Stephens"),
                     residual_variance_lowerbound = 0,
                     residual_variance_upperbound = Inf,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = c("optim", "EM", "simple"),
                     unmappable_effects = c("none", "inf"),
                     check_null_threshold = 0,
                     prior_tol = 1e-9,
                     max_iter = 100,
                     tol = 1e-3,
                     convergence_method = c("elbo", "pip"),
                     coverage = 0.95,
                     min_abs_corr = 0.5,
                     n_purity = 100,
                     verbose = FALSE,
                     track_fit = FALSE,
                     check_prior = FALSE,
                     refine = FALSE,
                     ...) {

  # Validate method arguments
  unmappable_effects       <- match.arg(unmappable_effects)
  estimate_prior_method    <- match.arg(estimate_prior_method)
  estimate_residual_method <- match.arg(estimate_residual_method)
  convergence_method       <- match.arg(convergence_method)

  # Construct data and params objects
  susie_objects <- sufficient_stats_constructor(
    XtX, Xty, yty, n, L, X_colmeans, y_mean,
    maf = NULL, maf_thresh = 0, check_input, r_tol, standardize,
    scaled_prior_variance, residual_variance, prior_weights, null_weight,
    model_init, estimate_residual_variance, estimate_residual_method,
    residual_variance_lowerbound, residual_variance_upperbound,
    estimate_prior_variance, estimate_prior_method, unmappable_effects,
    check_null_threshold, prior_tol, max_iter, tol, convergence_method,
    coverage, min_abs_corr, n_purity, verbose, track_fit, check_prior,
    refine
  )

  # Run main SuSiE algorithm
  model <- susie_workhorse(susie_objects$data, susie_objects$params)

  return(model)
}

# =============================================================================
# SuSiE WITH REGRESSION SUMMARY STATISTICS
# =============================================================================

#' @title SuSiE with Regression Summary Statistics (RSS)
#'
#' @description Performs SuSiE regression using z-scores and correlation matrix.
#' Supports both standard RSS (lambda = 0) and RSS with correlated errors (lambda > 0).
#'
#' @param z A p-vector of z-scores.
#'
#' @param R A p by p correlation matrix.
#'
#' @param n The sample size, not required but recommended.
#'
#' @param bhat Alternative summary data giving the estimated effects
#'   (a vector of length p). This, together with \code{shat}, may be
#'   provided instead of \code{z}.
#'
#' @param shat Alternative summary data giving the standard errors of
#'   the estimated effects (a vector of length p). This, together with
#'   \code{bhat}, may be provided instead of \code{z}.
#'
#' @param lambda Regularization parameter for correlated errors. When
#' \code{lambda} > 0, you cannot use \code{unmappable_effects} methods.
#'
#' @param z_ld_weight This parameter is included for backwards
#'   compatibility with previous versions of the function, but it is no
#'   longer recommended to set this to a non-zero value. When
#'   \code{z_ld_weight > 0}, the matrix \code{R} is adjusted to be
#'   \code{cov2cor((1-w)*R + w*tcrossprod(z))}, where \code{w =
#'   z_ld_weight}.
#'
#' @param estimate_residual_variance The default is FALSE, the
#'   residual variance is fixed to 1 or variance of y. If the in-sample
#'   LD matrix is provided, we recommend setting
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @param check_R If TRUE, check that R is positive semidefinite.
#'
#' @param check_z If TRUE, check that z lies in column space of R.
#'
#' @export
susie_rss <- function(z = NULL, R, n = NULL,
                      bhat = NULL, shat = NULL, var_y = NULL,
                      L = min(10, ncol(R)),
                      lambda = 0,
                      maf = NULL,
                      maf_thresh = 0,
                      z_ld_weight = 0,
                      prior_variance = 50,
                      scaled_prior_variance = 0.2,
                      residual_variance = NULL,
                      prior_weights = NULL,
                      null_weight = 0,
                      standardize = TRUE,
                      intercept_value = 0,
                      estimate_residual_variance = FALSE,
                      estimate_residual_method = c("MLE", "MoM"),
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = c("optim", "EM", "simple"),
                      unmappable_effects = c("none", "inf"),
                      check_null_threshold = 0,
                      prior_tol = 1e-9,
                      residual_variance_lowerbound = 0,
                      residual_variance_upperbound = Inf,
                      s_init = NULL,
                      coverage = 0.95,
                      min_abs_corr = 0.5,
                      max_iter = 100,
                      tol = 1e-3,
                      convergence_method = c("elbo", "pip"),
                      verbose = FALSE,
                      track_fit = FALSE,
                      check_input = FALSE,
                      check_prior = TRUE,
                      check_R = TRUE,
                      check_z = FALSE,
                      n_purity = 100,
                      r_tol = 1e-8,
                      refine = FALSE) {

  # Validate method arguments
  unmappable_effects       <- match.arg(unmappable_effects)
  estimate_prior_method    <- match.arg(estimate_prior_method)
  estimate_residual_method <- match.arg(estimate_residual_method)
  convergence_method       <- match.arg(convergence_method)

  # Construct data and params objects with ALL parameters
  susie_objects <- summary_stats_constructor(
    z, R, n, bhat, shat, var_y, L, lambda, maf, maf_thresh,
    z_ld_weight, prior_variance, scaled_prior_variance,
    residual_variance, prior_weights, null_weight, standardize,
    intercept_value, estimate_residual_variance, estimate_residual_method,
    estimate_prior_variance, estimate_prior_method,
    unmappable_effects, check_null_threshold, prior_tol,
    residual_variance_lowerbound, residual_variance_upperbound,
    model_init = s_init, coverage, min_abs_corr,
    max_iter, tol, convergence_method, verbose, track_fit, check_input,
    check_prior, check_R, check_z, n_purity, r_tol, compute_univariate_zscore,
    refine
  )

  # Run main SuSiE algorithm
  model <- susie_workhorse(susie_objects$data, susie_objects$params)

  return(model)
}
