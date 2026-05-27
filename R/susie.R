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
#'   nonzero effect on the outcome, y. The weights are internally
#'   normalized to sum to 1. When \code{NULL} (the default), uniform
#'   prior weights are used (each variable is assigned probability
#'   \code{1/p}).
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
#'   the infinitesimal model, "MoM" is more stable. We recommend using "NIG"
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
#'   likelihood is retained. When \code{prior_variance_grid} is provided,
#'   this is automatically set to \code{"fixed_mixture"}.
#'
#' @param prior_variance_grid Numeric vector of K prior variances defining
#'   a mixture-of-normals prior on effect sizes. When provided, the SER
#'   evaluates Bayes factors at each grid point and forms a mixture BF
#'   weighted by \code{mixture_weights}. This bypasses the scalar prior
#'   variance optimization. Default is \code{NULL} (standard scalar V path).
#'
#' @param mixture_weights Numeric vector of K non-negative weights summing
#'   to 1, giving the mixture proportions for the variance grid. Default is
#'   \code{NULL}, which uses uniform weights when \code{prior_variance_grid}
#'   is provided.
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
#' @param s_init Deprecated alias for \code{model_init}.
#'
#' @param coverage A number between 0 and 1 specifying the
#'   \dQuote{coverage} of the estimated confidence sets.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies. This "purity" filter is
#'   applied to the CSs reported in the fit object, so the CS list
#'   returned here may be a subset of the one produced by calling
#'   \code{\link{susie_get_cs}} on the same fit without passing
#'   \code{X} or \code{Xcorr} (in which case the purity filter is
#'   skipped).
#'
#' @param compute_univariate_zscore If \code{compute_univariate_zscore
#'   = TRUE}, the univariate regression z-scores are outputted for each
#'   variable.
#'
#' @param na.rm Drop any missing values in y from both X and y.
#'
#' @param max_iter Maximum number of IBSS iterations to perform. For
#'   \code{susie_rss()} and \code{susie_rss_lambda()}, \code{NULL} uses
#'   \code{50} with a hint; other interfaces use \code{100}.
#'
#' @param L_greedy Integer or \code{NULL}. When non-\code{NULL}, run a
#'   greedy outer loop that grows the number of effects from
#'   \code{L_greedy} up to \code{L} in linear steps until the fit
#'   saturates. The default \code{NULL} runs the usual fixed-\code{L}
#'   fit.
#'
#' @param greedy_lbf_cutoff Numeric saturation threshold for the
#'   \code{L_greedy} outer loop. Default is 0.1.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. When \code{NULL}, the
#'   default is \code{1e-4}; for \code{estimate_residual_method = "NIG"},
#'   the default is tightened to \code{1e-6}.
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
#' @param track_fit If \code{track_fit = TRUE}, a compact
#'   \code{susie_track} object is returned in \code{trace}, containing
#'   alpha history, effect summaries and available diagnostics at each
#'   iteration of the IBSS fitting procedure.
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
#'   \code{estimate_residual_method = "NIG"}. Defaults to
#'   \code{1/sqrt(n)}, where \code{n} is the sample size. When calling
#'   \code{susie_rss} with NIG, \code{n} must be supplied; otherwise
#'   validation errors.
#'
#' @param beta0 Numerical parameter for the NIG prior when using
#'   \code{estimate_residual_method = "NIG"}. Defaults to
#'   \code{1/sqrt(n)}, where \code{n} is the sample size. When calling
#'   \code{susie_rss} with NIG, \code{n} must be supplied; otherwise
#'   validation errors.
#'
#' @param slot_prior Optional slot activity prior created by
#'   \code{\link{slot_prior_betabinom}} or \code{\link{slot_prior_poisson}}.
#'   Use \code{slot_prior_betabinom(a_beta, b_beta)} for the usual
#'   single-locus setting; it places a Beta-Binomial prior on the
#'   number of active effects and gives an adaptive multiplicity
#'   correction. Use \code{slot_prior_poisson(C, nu)} when you want a
#'   Gamma-Poisson prior centered on an expected number \code{C} of
#'   active effects. When supplied, each single-effect slot has an
#'   estimated activity probability \code{c_hat}; fitted values and
#'   PIPs are weighted by these activity probabilities, and convergence
#'   is checked using \code{convergence_method = "pip"}.
#'
#' @param init_only Logical. If \code{TRUE}, return a list with
#'   \code{data} and \code{params} objects without running the IBSS
#'   algorithm. Used by packages like susieAnn that implement their own
#'   outer loop around SuSiE's building blocks. Default is \code{FALSE}.
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
#'   means for the unmappable effects.}
#'
#' \item{tau2}{If \code{unmappable_effects = "inf"} or
#'   \code{unmappable_effects = "ash"}, then \code{tau2} is the unmappable variance.}
#'
#' @importFrom stats var
#' 
#' @export
#' 
susie <- function(X, y, L = min(10, ncol(X)),
                  scaled_prior_variance = 0.2,
                  residual_variance = NULL,
                  prior_weights = NULL,
                  null_weight = 0,
                  standardize = TRUE,
                  intercept = TRUE,
                  estimate_residual_variance = TRUE,
                  estimate_residual_method = c("MoM", "MLE", "NIG"),
                  estimate_prior_variance = TRUE,
                  estimate_prior_method = c("optim", "EM", "simple"),
                  prior_variance_grid = NULL,
                  mixture_weights = NULL,
                  unmappable_effects = c("none", "inf", "ash", "ash_filter_archived"),
                  check_null_threshold = 0,
                  prior_tol = 1e-9,
                  residual_variance_upperbound = Inf,
                  model_init = NULL,
                  s_init = NULL,
                  coverage = 0.95,
                  min_abs_corr = 0.5,
                  compute_univariate_zscore = FALSE,
                  na.rm = FALSE,
                  max_iter = 100,
                  L_greedy = NULL,
                  greedy_lbf_cutoff = 0.1,
                  tol = NULL,
                  convergence_method = c("elbo", "pip"),
                  verbose = FALSE,
                  track_fit = FALSE,
                  residual_variance_lowerbound = NULL,
                  refine = FALSE,
                  n_purity = 100,
                  alpha0 = NULL,
                  beta0 = NULL,
                  init_only = FALSE,
                  slot_prior = NULL) {

  # Validate method arguments
  unmappable_effects       <- match.arg(unmappable_effects)
  estimate_residual_method <- match.arg(estimate_residual_method)
  convergence_method       <- match.arg(convergence_method)
  mp <- resolve_mixture_prior(estimate_prior_method, estimate_prior_variance,
                              prior_variance_grid, mixture_weights)
  estimate_prior_method   <- mp$estimate_prior_method
  estimate_prior_variance <- mp$estimate_prior_variance
  prior_variance_grid     <- mp$prior_variance_grid
  mixture_weights         <- mp$mixture_weights

  # Construct data and params objects
  susie_objects <- individual_data_constructor(
    X, y, L, scaled_prior_variance, residual_variance,
    prior_weights, null_weight, standardize, intercept,
    estimate_residual_variance, estimate_residual_method,
    estimate_prior_variance, estimate_prior_method,
    prior_variance_grid, mixture_weights,
    unmappable_effects, check_null_threshold, prior_tol,
    residual_variance_upperbound, model_init, s_init, coverage,
    min_abs_corr, compute_univariate_zscore, na.rm,
    max_iter, tol, convergence_method, verbose, track_fit,
    residual_variance_lowerbound, refine, n_purity,
    alpha0, beta0, slot_prior, L_greedy, greedy_lbf_cutoff
  )

  # Return data and params without fitting if init_only is TRUE.
  # The caller is responsible for calling ibss_initialize() on these.
  if (init_only) {
    return(susie_objects)
  }

  # Run main SuSiE algorithm
  model <- susie_workhorse(susie_objects$data, susie_objects$params)

  return(model)
}

# =============================================================================
# SuSiE WITH SUFFICIENT STATISTICS
# =============================================================================

#' @title SuSiE using Sufficient Statistics
#'
#' @inheritParams susie
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
#' 
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
                     s_init = NULL,
                     estimate_residual_variance = TRUE,
                     estimate_residual_method = c("MoM", "MLE", "NIG"),
                     residual_variance_lowerbound = 0,
                     residual_variance_upperbound = Inf,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = c("optim", "EM", "simple"),
                     prior_variance_grid = NULL,
                     mixture_weights = NULL,
                     unmappable_effects = c("none", "inf", "ash", "ash_filter_archived"),
                     check_null_threshold = 0,
                     prior_tol = 1e-9,
                     max_iter = 100,
                     L_greedy = NULL,
                     greedy_lbf_cutoff = 0.1,
                     tol = NULL,
                     convergence_method = c("elbo", "pip"),
                     coverage = 0.95,
                     min_abs_corr = 0.5,
                     n_purity = 100,
                     verbose = FALSE,
                     track_fit = FALSE,
                     check_prior = FALSE,
                     refine = FALSE,
                     alpha0 = NULL,
                     beta0 = NULL,
                     slot_prior = NULL) {

  # Validate method arguments
  unmappable_effects       <- match.arg(unmappable_effects)
  estimate_residual_method <- match.arg(estimate_residual_method)
  convergence_method       <- match.arg(convergence_method)
  mp <- resolve_mixture_prior(estimate_prior_method, estimate_prior_variance,
                              prior_variance_grid, mixture_weights)
  estimate_prior_method   <- mp$estimate_prior_method
  estimate_prior_variance <- mp$estimate_prior_variance
  prior_variance_grid     <- mp$prior_variance_grid
  mixture_weights         <- mp$mixture_weights

  # Construct data and params objects
  susie_objects <- sufficient_stats_constructor(
    Xty = Xty, yty = yty, n = n, XtX = XtX,
    L = L, X_colmeans = X_colmeans, y_mean = y_mean,
    maf = maf, maf_thresh = maf_thresh,
    check_input = check_input, r_tol = r_tol, standardize = standardize,
    scaled_prior_variance = scaled_prior_variance,
    residual_variance = residual_variance,
    prior_weights = prior_weights, null_weight = null_weight,
    model_init = model_init, s_init = s_init,
    estimate_residual_variance = estimate_residual_variance,
    estimate_residual_method = estimate_residual_method,
    residual_variance_lowerbound = residual_variance_lowerbound,
    residual_variance_upperbound = residual_variance_upperbound,
    estimate_prior_variance = estimate_prior_variance,
    estimate_prior_method = estimate_prior_method,
    prior_variance_grid = prior_variance_grid,
    mixture_weights = mixture_weights,
    unmappable_effects = unmappable_effects,
    check_null_threshold = check_null_threshold, prior_tol = prior_tol,
    max_iter = max_iter, tol = tol, convergence_method = convergence_method,
    coverage = coverage, min_abs_corr = min_abs_corr, n_purity = n_purity,
    verbose = verbose, track_fit = track_fit, check_prior = check_prior,
    refine = refine, alpha0 = alpha0, beta0 = beta0,
    slot_prior = slot_prior, L_greedy = L_greedy,
    greedy_lbf_cutoff = greedy_lbf_cutoff
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
#' @inheritParams susie_ss
#' 
#' @description Performs SuSiE regression using z-scores and correlation matrix.
#' This is the sufficient-statistics RSS interface. For the specialized
#' regularized eigendecomposition likelihood with \code{lambda > 0}, use
#' \code{\link{susie_rss_lambda}}.
#'
#' @param z A p-vector of z-scores.
#'
#' @param R A p by p correlation matrix. Exactly one of \code{R} or
#'   \code{X} must be provided.
#'
#' @param n The sample size, not required but recommended.
#'
#' @param X A factor matrix (B x p) such that \code{R = crossprod(X) /
#'   nrow(X)} approximates the R (correlation) matrix. When
#'   \code{nrow(X) >= ncol(X)}, the correlation matrix \code{R} is
#'   formed explicitly and the standard path is used. When
#'   \code{nrow(X) < ncol(X)}, a low-rank path is used that avoids
#'   forming the p x p matrix, reducing per-iteration cost from
#'   O(Lp^2) to O(LBp). Columns of \code{X} are standardized
#'   internally.
#'
#' @param bhat Alternative summary data giving the estimated effects
#'   (a vector of length p). This, together with \code{shat}, may be
#'   provided instead of \code{z}.
#'
#' @param shat Alternative summary data giving the standard errors of
#'   the estimated effects (a vector of length p). This, together with
#'   \code{bhat}, may be provided instead of \code{z}.
#'
#' @param var_y The sample variance of y, defined as \eqn{y'y/(n-1)}.
#'   When the sample variance is not provided, the coefficients
#'   (returned from \code{coef}) are computed on the
#'   \dQuote{standardized} X, y scale.
#' 
#' @param estimate_residual_variance The default is FALSE, the
#'   residual variance is fixed to 1 or variance of y. If the in-sample
#'   R matrix is provided, we recommend setting
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @param R_finite Controls variance inflation to account
#'   for estimating the R matrix from a finite reference panel. Accepts four
#'   types of input:
#'   \describe{
#'     \item{\code{NULL} (default)}{The R matrix is treated as trusted, and no
#'       finite-reference variance inflation is applied. With
#'       \code{estimate_residual_variance = TRUE}, this keeps the in-sample
#'       R warning active.}
#'     \item{\code{FALSE}}{No finite-reference variance inflation is applied,
#'       and the in-sample R warning for
#'       \code{estimate_residual_variance = TRUE} is silenced. Use this only
#'       when \code{R} is the in-sample correlation matrix.}
#'     \item{\code{TRUE}}{Infer the reference sample size B from the input
#'       \code{X}. Sets \code{B = nrow(X)} for single-panel input,
#'       or one B per panel for multi-panel
#'       input. Requires \code{X} to be provided (errors if only
#'       \code{R} is given, since B cannot be inferred).}
#'     \item{Number}{Explicit reference sample size B.}
#'   }
#'   When active, this dynamically inflates the null variance of each
#'   variable's score statistic at every IBSS iteration to account for
#'   finite-reference uncertainty in the Single Effect Regression (SER).
#'   When provided, the output includes a
#'   \code{R_finite_diagnostics} element with per-region and
#'   per-variable quality metrics.
#'
#' @param R_mismatch R-mismatch correction mode. \code{"none"} (default) is
#'   off. \code{"eb"} is the recommended empirical Bayes procedure for
#'   mismatch correction described in Sun et al. (2026+). It updates a
#'   region-level variance component during the IBSS iterations and reports a
#'   QC score (\code{Q_art}) that extends the Zou et al. (2022) column-space
#'   check from the input summary vector to the fitted residual after
#'   correction. It warns when that residual still projects onto near-null
#'   directions of the supplied \code{R}, and auto-disables
#'   \code{estimate_residual_variance} with a warning.
#'
#' @param R_mismatch_method Estimator for the region-level
#'   \code{lambda_bias} variance component when \code{R_mismatch != "none"}.
#'   \code{"mle"} (default) maximizes the working Gaussian likelihood.
#'   \code{"map"} uses a half-Cauchy MAP estimator.
#'
#' @param eig_delta_rel,eig_delta_abs Cutoffs for "low-eigenvalue"
#'   directions of \code{R} used by the QC diagnostic when
#'   \code{R_mismatch != "none"}. Default \code{eig_delta_rel = 1e-3},
#'   \code{eig_delta_abs = 0}; the threshold is
#'   \code{max(eig_delta_abs, eig_delta_rel * max_eigenvalue(R))}. Tighter
#'   (smaller) values flag fewer regions.
#'
#' @param artifact_threshold Flag threshold on the QC score \code{Q_art}
#'   (a fraction in [0, 1]). Default \code{0.1}; flag fires when
#'   \code{Q_art > artifact_threshold}. Heuristic, not a calibrated test.
#'
#' @param R_sensitivity_threshold Flag threshold for the credible-set
#'   Bayes-factor attenuation diagnostic. Default \code{log(20)}; flag fires
#'   when a credible set contains a variable whose nominal log BF exceeds its
#'   R-adjusted log BF by at least this amount.
#'
#' @param init_only Logical. If \code{TRUE}, return a list with
#'   \code{data} and \code{params} objects without running the IBSS
#'   algorithm. Default is \code{FALSE}.
#'
#' @return In addition to the standard \code{"susie"} output (see
#'   \code{\link{susie}}), the returned object may contain:
#'
#' \item{R_finite_diagnostics}{A list of diagnostics for the
#'   R-uncertainty correction (only present when \code{R_finite} is provided
#'   or \code{R_mismatch != "none"}), containing:
#'   \code{B} (the reference sample size);
#'   \code{p} (number of variables);
#'   \code{effective_rank} (debiased \eqn{\tilde{r} = p^2 / \|R\|_F^2});
#'   \code{r_over_B} (\eqn{\tilde{r}/B}, one number per region; values
#'     \eqn{\le 0.2} indicate the reference panel is adequate);
#'   \code{Rhat_diag_deviation} (\eqn{|\hat{R}_{jj} - 1|}, one number
#'     per variable);
#'   \code{lambda_bias} (region-level scalar on the default
#'     \code{lambda = 0} sufficient-statistics path when
#'     \code{R_mismatch != "none"});
#'   \code{B_corrected} (effective reference sample size after the
#'     R-mismatch correction, \eqn{1/(1/B + \lambda_{\mathrm{bias}})};
#'     substantially
#'     smaller than the input \code{B} flags a dominant population
#'     mismatch component);
#'   \code{per_variable_penalty} (final-iteration
#'     \eqn{v_j / \sigma^2 = \tau_j^2 / \sigma^2 - 1}, one number per
#'     variable; values \eqn{\le 0.2} indicate minimal power loss,
#'     values \eqn{\gg 1} flag variables where the correction is doing
#'     heavy lifting).}
#'
#' @export
#' 
susie_rss <- function(z = NULL, R = NULL, n = NULL,
                      X = NULL,
                      bhat = NULL, shat = NULL, var_y = NULL,
                      L = min(10, if (is.list(R) && !is.matrix(R)) ncol(R[[1]])
                               else if (!is.null(R)) ncol(R)
                               else if (is.list(X) && !is.matrix(X)) ncol(X[[1]])
                               else ncol(X)),
                      maf = NULL,
                      maf_thresh = 0,
                      scaled_prior_variance = 0.2,
                      residual_variance = NULL,
                      prior_weights = NULL,
                      null_weight = 0,
                      standardize = TRUE,
                      estimate_residual_variance = FALSE,
                      estimate_residual_method = c("MoM", "MLE", "NIG"),
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = c("optim", "EM", "simple"),
                      prior_variance_grid = NULL,
                      mixture_weights = NULL,
                      unmappable_effects = c("none", "inf", "ash", "ash_filter_archived"),
                      check_null_threshold = 0,
                      prior_tol = 1e-9,
                      residual_variance_lowerbound = 0,
                      residual_variance_upperbound = Inf,
                      model_init = NULL,
                      s_init = NULL,
                      coverage = 0.95,
                      min_abs_corr = 0.5,
                      max_iter = NULL,
                      L_greedy = NULL,
                      greedy_lbf_cutoff = 0.1,
                      tol = NULL,
                      convergence_method = c("elbo", "pip"),
                      verbose = FALSE,
                      track_fit = FALSE,
                      check_input = FALSE,
                      check_prior = TRUE,
                      n_purity = 100,
                      r_tol = 1e-8,
                      refine = FALSE,
                      R_finite = NULL,
                      R_mismatch = c("none", "eb"),
                      R_mismatch_method = c("mle", "map"),
                      eig_delta_rel = 1e-3,
                      eig_delta_abs = 0,
                      artifact_threshold = 0.1,
                      R_sensitivity_threshold = log(20),
                      alpha0 = NULL,
                      beta0 = NULL,
                      init_only = FALSE,
                      slot_prior = NULL) {

  # Validate method arguments
  unmappable_effects       <- match.arg(unmappable_effects)
  estimate_residual_method <- match.arg(estimate_residual_method)
  convergence_method       <- match.arg(convergence_method)
  if (length(R_mismatch) > 1L)
    R_mismatch <- R_mismatch[1L]
  R_mismatch <- match.arg(R_mismatch,
                          c("none", "eb", "eb_ser_init",
                            "eb_force_init", "eb_no_init"))
  R_mismatch_method        <- match.arg(R_mismatch_method)
  mp <- resolve_mixture_prior(estimate_prior_method, estimate_prior_variance,
                              prior_variance_grid, mixture_weights)
  estimate_prior_method   <- mp$estimate_prior_method
  estimate_prior_variance <- mp$estimate_prior_variance
  prior_variance_grid     <- mp$prior_variance_grid
  mixture_weights         <- mp$mixture_weights

  susie_objects <- summary_stats_constructor(
    z = z, R = R, X = X, n = n,
    bhat = bhat, shat = shat, var_y = var_y,
    L = L, maf = maf, maf_thresh = maf_thresh,
    scaled_prior_variance = scaled_prior_variance,
    residual_variance = residual_variance,
    prior_weights = prior_weights, null_weight = null_weight,
    standardize = standardize,
    estimate_residual_variance = estimate_residual_variance,
    estimate_residual_method = estimate_residual_method,
    estimate_prior_variance = estimate_prior_variance,
    estimate_prior_method = estimate_prior_method,
    prior_variance_grid = prior_variance_grid,
    mixture_weights = mixture_weights,
    unmappable_effects = unmappable_effects,
    check_null_threshold = check_null_threshold, prior_tol = prior_tol,
    residual_variance_lowerbound = residual_variance_lowerbound,
    residual_variance_upperbound = residual_variance_upperbound,
    model_init = model_init, s_init = s_init,
    coverage = coverage, min_abs_corr = min_abs_corr,
    max_iter = max_iter, tol = tol, convergence_method = convergence_method,
    verbose = verbose, track_fit = track_fit, check_input = check_input,
    check_prior = check_prior,
    n_purity = n_purity, r_tol = r_tol, refine = refine,
    R_finite = R_finite,
    R_mismatch = R_mismatch,
    R_mismatch_method = R_mismatch_method,
    eig_delta_rel = eig_delta_rel, eig_delta_abs = eig_delta_abs,
    artifact_threshold = artifact_threshold,
    R_sensitivity_threshold = R_sensitivity_threshold,
    alpha0 = alpha0, beta0 = beta0,
    slot_prior = slot_prior, L_greedy = L_greedy,
    greedy_lbf_cutoff = greedy_lbf_cutoff
  )

  # Return constructed data and params without running IBSS (for susieAnn
  # and other packages that implement their own outer loop). The caller
  # is responsible for calling ibss_initialize() on the returned objects.
  if (init_only) {
    return(susie_objects)
  }

  # Run main SuSiE algorithm
  model <- susie_workhorse(susie_objects$data, susie_objects$params)

  # Attach multi-panel metadata when the constructor ran sub-fits to pick
  # the mixture-init panel. Always return the mixture result; users can
  # pick a single-panel fit from `model$single_panel_fits` if they prefer.
  mp_meta <- susie_objects$multi_panel_meta
  if (!is.null(mp_meta)) {
    if (verbose) {
      mix_elbo <- tail(model$elbo, 1)
      best_k <- which.max(mp_meta$elbos)
      omega_str <- paste(round(model$omega_weights, 3), collapse = ", ")
      message(sprintf(
        "Multi-panel: mixture ELBO = %.2f (omega = %s), best single-panel ELBO = %.2f (panel %d).",
        mix_elbo, omega_str, mp_meta$elbos[best_k], best_k))
    }
    model$single_panel_fits <- mp_meta$fits
  }

  return(model)
}

#' Sum of Single Effects Regression using the RSS-lambda likelihood
#'
#' @description Specialized interface for the regularized eigendecomposition
#' RSS likelihood of Zou et al. (2022). This path accepts a single reference
#' matrix or a single factor matrix and does not support multi-panel mixture,
#' finite-reference inflation, or R-mismatch correction.
#'
#' @inheritParams susie_rss
#'
#' @param X Optional factor matrix (B by p) such that
#'   \code{R = crossprod(X) / nrow(X)} approximates the correlation
#'   matrix \code{R}. The RSS-lambda path always builds the eigen
#'   decomposition via SVD of standardized \code{X}; there is no
#'   separate low-rank branch. Provide either \code{R} or \code{X}, not
#'   both.
#' @param lambda Regularization parameter for the RSS-lambda likelihood.
#'   Must be supplied. \code{lambda = "estimate"} estimates lambda from
#'   the null-space residual.
#' @param prior_variance Prior variance for each non-zero effect on the
#'   z-score scale. Replaces \code{scaled_prior_variance} from
#'   \code{\link{susie_rss}}. Default \code{50}.
#' @param residual_variance Residual variance on the RSS-lambda scale. If
#'   \code{NULL}, it is initialized to \code{1 - lambda}; if supplied, the
#'   working residual variance is \code{residual_variance - lambda}.
#' @param intercept_value Intercept used by the RSS-lambda likelihood.
#'   Default \code{0}.
#' @param estimate_residual_variance If \code{TRUE}, estimate the
#'   RSS-lambda residual variance by maximum likelihood.
#' @param estimate_residual_method Variance-component estimator. The
#'   RSS-lambda path supports \code{"MLE"} only; any other value errors.
#' @param estimate_prior_variance If \code{estimate_prior_variance = TRUE},
#'   the prior variance is estimated (a separate parameter for each of
#'   the L effects). When \code{TRUE}, \code{prior_variance} provides the
#'   initial value; when \code{FALSE}, it is held fixed.
#' @param check_null_threshold When the prior variance is estimated,
#'   compare its likelihood to the likelihood at zero and use zero
#'   unless the larger value exceeds it by at least
#'   \code{check_null_threshold}. \code{0} (default) takes the larger
#'   likelihood at face value.
#' @param check_R If TRUE, verify that \code{R} is positive semidefinite.
#' @param check_z If TRUE, verify that \code{z} lies in the column space
#'   of \code{R}.
#'
#' @return A \code{"susie"} fit (or, with \code{init_only = TRUE}, the
#'   constructed data and params objects).
#'
#' @export
susie_rss_lambda <- function(z = NULL, R = NULL, n = NULL,
                             X = NULL,
                             L = min(10, if (!is.null(R)) ncol(R)
                                      else ncol(X)),
                             lambda,
                             maf = NULL,
                             maf_thresh = 0,
                             prior_variance = 50,
                             residual_variance = NULL,
                             prior_weights = NULL,
                             null_weight = 0,
                             intercept_value = 0,
                             estimate_residual_variance = FALSE,
                             estimate_residual_method = "MLE",
                             estimate_prior_variance = TRUE,
                             estimate_prior_method = c("optim", "EM", "simple"),
                             prior_variance_grid = NULL,
                             mixture_weights = NULL,
                             check_null_threshold = 0,
                             prior_tol = 1e-9,
                             residual_variance_lowerbound = 0,
                             model_init = NULL,
                             coverage = 0.95,
                             min_abs_corr = 0.5,
                             max_iter = NULL,
                             L_greedy = NULL,
                             greedy_lbf_cutoff = 0.1,
                             tol = NULL,
                             convergence_method = c("elbo", "pip"),
                             verbose = FALSE,
                             track_fit = FALSE,
                             check_prior = TRUE,
                             check_R = TRUE,
                             check_z = FALSE,
                             n_purity = 100,
                             r_tol = 1e-8,
                             refine = FALSE,
                             init_only = FALSE,
                             slot_prior = NULL) {
  if (missing(lambda))
    stop("susie_rss_lambda() requires lambda.")

  convergence_method       <- match.arg(convergence_method)
  mp <- resolve_mixture_prior(estimate_prior_method, estimate_prior_variance,
                              prior_variance_grid, mixture_weights)
  estimate_prior_method   <- mp$estimate_prior_method
  estimate_prior_variance <- mp$estimate_prior_variance
  prior_variance_grid     <- mp$prior_variance_grid
  mixture_weights         <- mp$mixture_weights

  susie_objects <- rss_lambda_constructor(
    z = z, R = R, X = X, n = n,
    L = L, lambda = lambda, maf = maf, maf_thresh = maf_thresh,
    prior_variance = prior_variance,
    residual_variance = residual_variance,
    prior_weights = prior_weights, null_weight = null_weight,
    intercept_value = intercept_value,
    estimate_residual_variance = estimate_residual_variance,
    estimate_residual_method = estimate_residual_method,
    estimate_prior_variance = estimate_prior_variance,
    estimate_prior_method = estimate_prior_method,
    prior_variance_grid = prior_variance_grid,
    mixture_weights = mixture_weights,
    check_null_threshold = check_null_threshold, prior_tol = prior_tol,
    residual_variance_lowerbound = residual_variance_lowerbound,
    model_init = model_init, coverage = coverage, min_abs_corr = min_abs_corr,
    max_iter = max_iter, tol = tol, convergence_method = convergence_method,
    verbose = verbose, track_fit = track_fit,
    check_prior = check_prior, check_R = check_R, check_z = check_z,
    n_purity = n_purity, r_tol = r_tol, refine = refine,
    slot_prior = slot_prior, L_greedy = L_greedy,
    greedy_lbf_cutoff = greedy_lbf_cutoff
  )

  if (init_only)
    return(susie_objects)

  susie_workhorse(susie_objects$data, susie_objects$params)
}
