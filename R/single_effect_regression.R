#' @rdname single_effect_regression
#'
#' @title Bayesian single-effect linear regression
#'
#' @description These methods fit the regression model \eqn{y = Xb +
#'   e}, where elements of e are \emph{i.i.d.}  \eqn{N(0,s^2)}, and b is
#'   a p-vector of effects to be estimated. The assumption is that b has
#'   exactly one non-zero element, with all elements equally likely to
#'   be non-zero. The prior on the coefficient of the non-zero element
#'   is \eqn{N(0,V)}.
#'
#' @details \code{single_effect_regression_ss} performs single-effect
#' linear regression with summary data, in which only the statistcs
#' \eqn{X^Ty} and diagonal elements of \eqn{X^TX} are provided to the
#' method.
#'
#' \code{single_effect_regression_rss} performs single-effect linear
#' regression with z scores. That is, this function fits the
#' regression model \eqn{z = R*b + e}, where e is \eqn{N(0,Sigma)},
#' \eqn{Sigma = residual_var*R + lambda*I}, and the b is a p-vector of
#' effects to be estimated. The assumption is that b has exactly one
#' non-zero element, with all elements equally likely to be non-zero.
#' The prior on the non-zero element is \eqn{N(0,V)}. The required
#' summary data are the p-vector \code{z} and the p by p matrix
#' \code{Sigma}. The summary statistics should come from the same
#' individuals.
#'
#' @param y An n-vector.
#'
#' @param X An n by p matrix of covariates.
#'
#' @param V A scalar giving the (initial) prior variance
#'
#' @param residual_variance The residual variance.
#'
#' @param prior_weights A p-vector of prior weights.
#'
#' @param optimize_V The optimization method to use for fitting the
#'   prior variance.
#'
#' @param check_null_threshold Scalar specifying threshold on the
#'   log-scale to compare likelihood between current estimate and zero
#'   the null.
#'
#' @return A list with the following elements:
#'
#' \item{alpha}{Vector of posterior inclusion probabilities;
#'   \code{alpha[i]} is posterior probability that the ith coefficient
#'   is non-zero.}
#'
#' \item{mu}{Vector of posterior means (conditional on inclusion).}
#'
#' \item{mu2}{Vector of posterior second moments (conditional on
#'   inclusion).}
#'
#' \item{lbf}{Vector of log-Bayes factors for each variable.}
#'
#' \item{lbf_model}{Log-Bayes factor for the single effect regression.}
#'
#' \code{single_effect_regression} and \code{single_effect_regression_ss}
#' additionally output:
#'
#' \item{V}{Prior variance (after optimization if \code{optimize_V !=
#'   "none"}).}
#'
#' \item{loglik}{The log-likelihood, \eqn{\log p(y | X, V)}.}
#'
#' @importFrom stats dnorm
#' @importFrom stats uniroot
#' @importFrom stats optim
#' @importFrom Matrix colSums
#'
#' @keywords internal
#'
single_effect_regression =
  function (y, X, V, residual_variance = 1, prior_weights = NULL,
            optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
            check_null_threshold = 0) {
  optimize_V = match.arg(optimize_V)
  Xty = compute_Xty(X,y)
  betahat = (1/attr(X,"d")) * Xty
  shat2 = residual_variance/attr(X,"d")
  if (is.null(prior_weights))
    prior_weights = rep(1/ncol(X),ncol(X))
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

  # Posterior prob for each SNP.
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  post_var = (1/V + attr(X,"d")/residual_variance)^(-1) # Posterior variance.
  post_mean = (1/residual_variance) * post_var * Xty
  post_mean2 = post_var + post_mean^2 # Second moment.

  # BF for single effect model.
  lbf_model = maxlpo + log(weighted_sum_w)
  loglik = lbf_model + sum(dnorm(y,0,sqrt(residual_variance),log = TRUE))

  if(optimize_V == "EM")
    V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
                                alpha,post_mean2,
                                check_null_threshold = check_null_threshold)

  return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
              lbf_model = lbf_model,V = V,loglik = loglik))
}

# Estimate prior variance.
est_V_uniroot = function (betahat, shat2, prior_weights) {
  V.u = uniroot(negloglik.grad.logscale,c(-10,10),extendInt = "upX",
                betahat = betahat,shat2 = shat2,prior_weights = prior_weights)
  return(exp(V.u$root))
}

optimize_prior_variance = function (optimize_V, betahat, shat2, prior_weights,
                                    alpha = NULL, post_mean2 = NULL,
                                    V_init = NULL, check_null_threshold = 0) {
  V = V_init
  if (optimize_V != "simple") {
    if(optimize_V == "optim") {
      lV = optim(par = log(max(c(betahat^2-shat2,1),na.rm = TRUE)),
          fn = neg.loglik.logscale,betahat = betahat,shat2 = shat2,
          prior_weights = prior_weights,method = "Brent",lower = -30,
          upper = 15)$par
      ## if the estimated one is worse than current one, don't change it.
      if(neg.loglik.logscale(lV, betahat = betahat,shat2 = shat2,prior_weights = prior_weights) >
         neg.loglik.logscale(log(V), betahat = betahat,
                             shat2 = shat2,prior_weights = prior_weights)){
        lV = log(V)
      }
      V = exp(lV)
    } else if (optimize_V == "uniroot")
      V = est_V_uniroot(betahat,shat2,prior_weights)
    else if (optimize_V == "EM")
      V = sum(alpha * post_mean2)
    else
     stop("Invalid option for optimize_V method")
  }

  # Set V exactly 0 if that beats the numerical value by
  # check_null_threshold in loglik. check_null_threshold = 0.1 is
  # exp(0.1) = 1.1 on likelihood scale; it means that for parsimony
  # reasons we set estiate of V to zero, if its numerical estimate is
  # only "negligibly" different from zero. We use a likelihood ratio
  # of exp(check_null_threshold) to define "negligible" in this
  # context. This is fairly modest condition compared to, say, a
  # formal LRT with p-value 0.05. But the idea is to be lenient to
  # non-zeros estimates unless they are indeed small enough to be
  # neglible. See more intuition at
  # https://stephens999.github.io/fiveMinuteStats/LR_and_BF.html
  if (loglik(0,betahat,shat2,prior_weights) +
      check_null_threshold >= loglik(V,betahat,shat2,prior_weights))
    V = 0
  return(V)
}

# In these functions, s2 represents residual_variance, and shat2 is an
# estimate of it.

# The log likelihood function for SER model (based on summary data
# betahat, shat2) as a function of prior variance V.
#
#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik = function (V, betahat, shat2, prior_weights) {

  #log(bf) for each SNP
  lbf = dnorm(betahat,0,sqrt(V+shat2),log = TRUE) -
        dnorm(betahat,0,sqrt(shat2),log = TRUE)
  lpo = lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0
  lpo[is.infinite(shat2)] = 0

  maxlpo = max(lpo)
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)
  return(log(weighted_sum_w) + maxlpo)
}

neg.loglik.logscale = function(lV,betahat,shat2,prior_weights)
  -loglik(exp(lV),betahat,shat2,prior_weights)

#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik.grad = function(V, betahat, shat2, prior_weights) {

  # log(bf) for each SNP.
  lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
        dnorm(betahat,0,sqrt(shat2),log = TRUE)
  lpo = lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0
  lpo[is.infinite(shat2)] = 0

  maxlpo = max(lpo)
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  return(sum(alpha * lbf.grad(V,shat2,betahat^2/shat2)))
}

# Define loglikelihood and gradient as function of lV:=log(V)
# to improve numerical optimization
negloglik.grad.logscale = function (lV, betahat, shat2, prior_weights)
  -exp(lV) * loglik.grad(exp(lV),betahat,shat2,prior_weights)

# Vector of gradients of logBF_j for each j, with respect to prior
# variance V.
lbf.grad = function (V, shat2, T2) {
  l = 0.5*(1/(V + shat2)) * ((shat2/(V + shat2))*T2 - 1)
  l[is.nan(l)] = 0
  return(l)
}

lbf = function (V, shat2, T2) {
  l = 0.5*log(shat2/(V + shat2)) + 0.5*T2*(V/(V + shat2))
  l[is.nan(l)] = 0
  return(l)
}
