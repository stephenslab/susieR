#' @rdname single_effect_regression
#' 
#' @param Xty A p-vector.
#' 
#' @param dXtX A p-vector containing the diagonal elements of
#'   \code{crossprod(X)}.
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

  if (optimize_V != "EM" && optimize_V != "none")
    V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
                                alpha = NULL,post_mean2 = NULL,V_init = V,
                                check_null_threshold = check_null_threshold)

  # log(po) = log(BF * prior) for each SNP
  lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
        dnorm(betahat,0,sqrt(shat2),log = TRUE)
  lpo = lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0 
  lpo[is.infinite(shat2)] = 0
  maxlpo = max(lpo)

  # w is proportional to
  #
  #   posterior odds = BF * prior,
  #
  # but subtract max for numerical stability.
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  post_var = (1/V + dXtX/residual_variance)^(-1) # Posterior variance.
  post_mean = (1/residual_variance) * post_var * Xty
  post_mean2 = post_var + post_mean^2 # Second moment.
  lbf_model = maxlpo + log(weighted_sum_w) # Analogue of loglik in the
                                           # non-summary case.

  if (optimize_V == "EM")
    V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,alpha,
          post_mean2,check_null_threshold = check_null_threshold)
  return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
              V = V,lbf_model = lbf_model))
}


