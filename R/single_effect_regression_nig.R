# Normal-inverse-gamma SER, with one convention throughout:
# b | sigma2, gamma=j ~ N(0, sigma2 * V), sigma2 ~ IG(alpha0/2,beta0/2).
# V is a variance RATIO, never an absolute variance or a standard deviation.
# n is the residual dimension after integrating out unpenalized covariates.
single_effect_regression_nig = function (xty, d, yty, n, V,
                                        prior_weights = NULL,
                                        optimize_V = "none",
                                        check_null_threshold = 0,
                                        alpha0 = 0, beta0 = 0) {
  scalar_nonnegative = function (x)
    is.numeric(x) && length(x) == 1 && is.finite(x) && x >= 0
  if (!scalar_nonnegative(V) || !scalar_nonnegative(alpha0) ||
      !scalar_nonnegative(beta0) || !scalar_nonnegative(check_null_threshold))
    stop("V, alpha0, beta0 and check_null_threshold must be finite nonnegative scalars")
  if (!scalar_nonnegative(n) || n + alpha0 <= 2)
    stop("NIG posterior second moments require n_effective + alpha0 > 2")
  if (!scalar_nonnegative(yty) || yty + beta0 <= 0)
    stop("NIG requires a positive yty + beta0")
  if (!length(d) || length(d) != length(xty) ||
      any(!is.finite(d)) || any(d < 0) || any(!is.finite(xty)))
    stop("Invalid NIG sufficient statistics")
  if (is.null(prior_weights))
    prior_weights = rep(1/length(d),length(d))
  if (length(prior_weights) != length(d) ||
      any(!is.finite(prior_weights)) || any(prior_weights < 0) ||
      sum(prior_weights) <= 0)
    stop("NIG prior_weights must be nonnegative with positive sum")
  prior_weights = prior_weights/sum(prior_weights)
  optimize_V = match.arg(optimize_V,c("none","optim","EM","simple","uniroot"))
  if (optimize_V == "uniroot")
    stop("Use optim, EM or simple for NIG prior estimation")

  # Orthogonal decomposition avoids cancellation near a perfect fit.
  explained = numeric(length(d))
  active = d > 0
  explained[active] = (xty[active]/sqrt(d[active]))^2
  if (any(explained > yty + 1e-8 * max(yty,.Machine$double.xmin)) ||
      any(xty[!active] != 0))
    stop("Inconsistent NIG sufficient statistics")
  rss = pmax(yty - explained,0)
  a = (n + alpha0)/2
  evaluate = function (v) {
    omega = if (v == 0) numeric(length(d)) else 1/(d + 1/v)
    mu = omega * xty
    rate2 = beta0 + rss + explained/(1 + v*d)
    if (any(rate2 <= 0))
      stop("Degenerate NIG posterior; use a positive beta0")
    lbf = -0.5*log1p(v*d) - a*log(rate2/(beta0 + yty))
    # Zero-information columns retain their prior odds, including null_weight.
    lbf[!active] = 0
    lpo = log(prior_weights) + lbf
    m = max(lpo)
    w = exp(lpo - m)
    alpha = w/sum(w)
    list(alpha = alpha,mu = mu,omega = omega,rate2 = rate2,
         lbf = lbf,lbf_model = m + log(sum(w)))
  }
  fit = evaluate(V)
  if (optimize_V == "optim" && any(active & prior_weights > 0)) {
    # Search on log(1 + V * mean(d)); retain zero and the initial value.
    # A coarse grid brackets possible separated modes of a mixture likelihood.
    dscale = mean(d[active & prior_weights > 0])
    upper = max(30,log1p(V*dscale) + 2)
    grid = seq(0,upper,length.out = 61)
    score = vapply(grid,function(t) evaluate(expm1(t)/dscale)$lbf_model,0.0)
    peaks = which(score >= c(-Inf,head(score,-1)) &
                  score >= c(utils::tail(score,-1),-Inf))
    candidates = c(0,V,expm1(grid[peaks])/dscale)
    for (k in peaks) {
      opt = optimize(function(t) -evaluate(expm1(t)/dscale)$lbf_model,
                     c(grid[max(1,k-1)],grid[min(length(grid),k+1)]))
      candidates = c(candidates,expm1(opt$minimum)/dscale)
    }
    scores = vapply(candidates,function(v) evaluate(v)$lbf_model,0.0)
    V = candidates[which.max(scores)]
    fit = evaluate(V)
  } else if (optimize_V == "EM") {
    # E[b^2/sigma2 | gamma=j] = omega_j + mu_j^2 E[1/sigma2 | j].
    # Recompute the posterior at the updated V, keeping outputs consistent.
    V = sum(fit$alpha * (fit$omega + fit$mu^2*(2*a)/fit$rate2))
    fit = evaluate(V)
  }
  if (optimize_V != "none" && fit$lbf_model <= check_null_threshold) {
    V = 0
    fit = evaluate(V)
  }
  sigma2 = fit$rate2/(2*a - 2)
  fit$mu2 = fit$mu^2 + fit$omega*sigma2
  fit$V = V
  fit$sigma2 = sum(fit$alpha*sigma2)
  fit$shape = a
  fit$rate = fit$rate2/2
  fit$scale2 = fit$omega*fit$rate/a
  # An improper common prior permits BF ratios, not absolute evidence.
  fit$loglik = if (alpha0 > 0 && beta0 > 0)
    fit$lbf_model - n/2*log(2*pi) + alpha0/2*log(beta0/2) -
      lgamma(alpha0/2) + lgamma(a) - a*log((beta0 + yty)/2)
  else NA_real_
  return(fit)
}
