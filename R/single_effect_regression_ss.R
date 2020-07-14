#' @title Bayesian single-effect linear regression of y on X using summary data
#' 
#' @description Performs single-effect linear regression of y on X, in
#'   which only the summary statistcs \eqn{X^T y} and diagonal elements
#'   of \eqn{X^T X} are available. That is, this function fits the
#'   regression model eqn{y = Xb + e}, where elements of e are
#'   i.i.d. \eqn{N(0,s2)}, and the b is a p-vector of effects to be
#'   estimated. The assumption is that b has exactly one non-zero
#'   element, with all elements equally likely to be non-zero. The prior
#'   on the non-zero element is \eqn{N(0,V)}.
#' 
#' @param Xty A p vector.
#' 
#' @param dXtX A p vector containing the diagonal elements of
#' \code{t(X) \%*\% X}.
#' 
#' @param V The prior variance.
#' 
#' @param residual_variance The residual variance.
#' 
#' @param prior_weights A p vector of prior weights.
#' 
#' @param optimize_V Specifies the method to optimize V (by maximum
#' likelihood).
#' 
#' @param check_null_threshold A sccalar threshold on the log-scale to
#'   compare the likelihood between current estimate and zero (the
#'   null).
#' 
#' @return A list with the following elements:
#' 
#' \item{alpha}{Vector of posterior inclusion probabilities. i.e.,
#'   \code{alpha[i]} is posterior probability that that the ith
#'   coefficient is non-zero.}
#' 
#' \item{mu}{Vector of posterior means (conditional on inclusion).}
#' 
#' \item{mu2}{Vector of posterior second moments (conditional on
#'   inclusion).}
#' 
#' \item{lbf}{Vector of log-Bayes factors for each variable.}
#' 
#' \item{V}{The prior variance (after optimization if
#'   \code{optimize_V != "none"}).}
#' 
#' \item{lbf_model}{The log-likelihood for the total model minus the
#'   log-likelihood for the null model.}
#'
#' @importFrom stats uniroot
#' @importFrom stats optim
#'
#' @keywords internal
#' 
single_effect_regression_ss =
  function (Xty, dXtX, V = 1, residual_variance = 1, prior_weights = NULL,
            optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
            check_null_threshold = 0) {
  optimize_V = match.arg(optimize_V)
  betahat = (1/dXtX) * Xty
  shat2 = residual_variance/dXtX
  if (is.null(prior_weights))
    prior_weights = rep(1/length(dXtX),length(dXtX))

  if(optimize_V != "EM" && optimize_V != "none")
    V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
                                alpha = NULL,post_mean2 = NULL,V_init = V,
                                check_null_threshold = check_null_threshold)

  # log(bf) for each SNP.
  lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
        dnorm(betahat,0,sqrt(shat2),log = TRUE)

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0 

  # w is proportional to BF, but subtract max for numerical stability
  # posterior prob on each SNP.
  maxlbf = max(lbf)
  w = exp(lbf - maxlbf) 
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w

  post_var = (1/V + dXtX/residual_variance)^(-1) # Posterior variance.
  post_mean = (1/residual_variance) * post_var * Xty
  post_mean2 = post_var + post_mean^2 # Second moment.
  lbf_model = maxlbf + log(weighted_sum_w) # Analogue of loglik in the
                                           # non-summary case.

  if (optimize_V == "EM")
    V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,alpha,
          post_mean2,check_null_threshold = check_null_threshold)
  return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
              V = V,lbf_model = lbf_model))
}


