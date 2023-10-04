#' @rdname single_effect_regression
#' 
#' @param z A p-vector of z scores.
#' 
#' @param Sigma \code{residual_var*R + lambda*I}
#' 
#' @keywords internal
#' 
single_effect_regression_rss =
  function (z, Sigma, V = 1, prior_weights = NULL,
            optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
            check_null_threshold = 0) {
  p = length(z)
  shat2 = 1/attr(Sigma,"RjSinvRj")
  if (is.null(prior_weights))
    prior_weights = rep(1/p,p)

  if (optimize_V != "EM" && optimize_V != "none") 
    V = optimize_prior_variance_rss(optimize_V,z,Sigma,prior_weights,
                                    alpha = NULL,post_mean2 = NULL,V_init = V,
                                    check_null_threshold=check_null_threshold)

  # log(po) = log(BF * prior) for each SNP
  lbf = sapply(1:p, function(j)
    -0.5 * log(1 + (V/shat2[j])) +
     0.5 * (V/(1 + (V/shat2[j]))) * sum(attr(Sigma,"SinvRj")[,j] * z)^2
  )
  lpo = lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does not
  # vary).
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
  post_var = (attr(Sigma,"RjSinvRj") + 1/V)^(-1) # Posterior variance.
  post_mean = sapply(1:p,function(j) (post_var[j]) *
              sum(attr(Sigma,"SinvRj")[,j] * z))
  post_mean2 = post_var + post_mean^2 # Second moment.
  lbf_model = maxlpo + log(weighted_sum_w) # Analogue of loglik in the
                                           # non-summary case.

  if (optimize_V=="EM") 
    V = optimize_prior_variance_rss(optimize_V,z,Sigma,prior_weights,
        alpha,post_mean2,check_null_threshold = check_null_threshold)
  
  return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
              V = V,lbf_model = lbf_model))
}

loglik_rss = function (V, z, Sigma, prior_weights) {
  p = length(z)
  shat2 = 1/attr(Sigma,"RjSinvRj")

  # log(bf) for each SNP.
  lbf = sapply(1:p,function (j)
    -0.5 * log(1 + (V/shat2[j])) +
     0.5 * (V/(1 + (V/shat2[j]))) * sum(attr(Sigma,"SinvRj")[,j] * z)^2)
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

neg.loglik_z.logscale_rss = function (lV, z, Sigma, prior_weights)
  -loglik_rss(exp(lV),z,Sigma,prior_weights)

optimize_prior_variance_rss = function (optimize_V, z, Sigma, prior_weights,
                                        alpha = NULL, post_mean2 = NULL,
                                        V_init = NULL,
                                        check_null_threshold = 0) {
  V = V_init
  if (optimize_V != "simple") {
    if(optimize_V == "optim") {
      lV = optim(par = log(max(c((colSums(attr(Sigma,"SinvRj") * z)^2) -
                       (1/attr(Sigma,"RjSinvRj")),1e-6),na.rm = TRUE)),
                 fn = neg.loglik_z.logscale_rss,z = z,Sigma = Sigma,
                 prior_weights = prior_weights,method = "Brent",
                 lower = -30,upper = 15)$par
      ## if the estimated one is worse than current one, don't change it.
      if(neg.loglik_z.logscale_rss(lV, z = z,Sigma = Sigma,prior_weights = prior_weights) > 
         neg.loglik_z.logscale_rss(log(V), z = z,Sigma = Sigma,prior_weights = prior_weights)){
        lV = log(V)
      }
      V = exp(lV)
    } else if (optimize_V == "EM")
      V = sum(alpha * post_mean2)
    else
      stop("Invalid option for optimize_V")
  }
  
  # Set V exactly 0 if that beats the numerical value. By
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
  if (loglik_rss(0,z,Sigma,prior_weights) + check_null_threshold >=
      loglik_rss(V,z,Sigma,prior_weights))
    V = 0
  return(V)
}
