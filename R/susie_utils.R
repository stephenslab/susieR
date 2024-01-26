#' @rdname susie_get_methods
#'
#' @title Inferences From Fitted SuSiE Model
#'
#' @description These functions access basic properties or draw
#'   inferences from a fitted susie model.
#'
#' @param res A susie fit, typically an output from
#'   \code{\link{susie}} or one of its variants. For
#'   \code{susie_get_pip} and \code{susie_get_cs}, this may instead be
#'   the posterior inclusion probability matrix, \code{alpha}.
#'
#' @param last_only If \code{last_only = FALSE}, return the ELBO from
#'   all iterations; otherwise return the ELBO from the last iteration
#'   only.
#'
#' @param warning_tol Warn if ELBO is decreasing by this
#'   tolerance level.
#'
#' @return \code{susie_get_objective} returns the evidence lower bound
#' (ELBO) achieved by the fitted susie model and, optionally, at each
#' iteration of the IBSS fitting procedure.
#'
#' \code{susie_get_residual_variance} returns the (estimated or
#' fixed) residual variance parameter.
#'
#' \code{susie_get_prior_variance} returns the (estimated or fixed)
#' prior variance parameters.
#'
#' \code{susie_get_posterior_mean} returns the posterior mean for the
#' regression coefficients of the fitted susie model.
#'
#' \code{susie_get_posterior_sd} returns the posterior standard
#' deviation for coefficients of the fitted susie model.
#'
#' \code{susie_get_niter} returns the number of model fitting
#' iterations performed.
#'
#' \code{susie_get_pip} returns a vector containing the posterior
#' inclusion probabilities (PIPs) for all variables.
#'
#' \code{susie_get_lfsr} returns a vector containing the average lfsr
#' across variables for each single-effect, weighted by the posterior
#' inclusion probability (alpha).
#'
#' \code{susie_get_posterior_samples} returns a list containing the
#' effect sizes samples and causal status with two components: \code{b},
#' an \code{num_variables} x \code{num_samples} matrix of effect
#' sizes; \code{gamma}, an \code{num_variables} x \code{num_samples}
#' matrix of causal status random draws.
#'
#' \code{susie_get_cs} returns credible sets (CSs) from a susie fit,
#' as well as summaries of correlation among the variables included in
#' each CS. If desired, one can filter out CSs that do not meet a
#' specified \dQuote{purity} threshold; to do this, either \code{X} or
#' \code{Xcorr} must be supplied. It returns a list with the following
#' elements:
#'
#' \item{cs}{A list in which each list element is a vector containing
#'   the indices of the variables in the CS.}
#'
#' \item{coverage}{The nominal coverage specified for each CS.}
#'
#' \item{purity}{If \code{X} or \code{Xcorr} iis provided), the
#'   purity of each CS.}
#'
#' \item{cs_index}{If \code{X} or \code{Xcorr} is provided) the index
#'   (number between 1 and L) of each reported CS in the supplied susie
#'   fit.}
#'
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#' s = susie(X,y,L = 10)
#' susie_get_objective(s)
#' susie_get_objective(s, last_only=FALSE)
#' susie_get_residual_variance(s)
#' susie_get_prior_variance(s)
#' susie_get_posterior_mean(s)
#' susie_get_posterior_sd(s)
#' susie_get_niter(s)
#' susie_get_pip(s)
#' susie_get_lfsr(s)
#'
#' @export
#'
susie_get_objective = function (res, last_only = TRUE, warning_tol = 1e-6) {
  if (!all(diff(res$elbo) >= (-1*warning_tol)))
    warning("Objective is decreasing")
  if (last_only)
    return(res$elbo[length(res$elbo)])
  else
    return(res$elbo)
}

#' @rdname susie_get_methods
#'
#' @export
#'
susie_get_posterior_mean = function (res, prior_tol = 1e-9) {

  # Drop the single-effects with estimated prior of zero.
  if (is.numeric(res$V))
    include_idx = which(res$V > prior_tol)
  else
    include_idx = 1:nrow(res$alpha)

  # Now extract relevant rows from alpha matrix.
  if (length(include_idx) > 0)
    return(colSums((res$alpha*res$mu)[include_idx,,drop=FALSE])/
           res$X_column_scale_factors)
  else
    return(numeric(ncol(res$mu)))
}

#' @rdname susie_get_methods
#'
#' @export
#'
susie_get_posterior_sd = function (res, prior_tol = 1e-9) {

  # Drop the single-effects with estimated prior of zero.
  if (is.numeric(res$V))
    include_idx = which(res$V > prior_tol)
  else
    include_idx = 1:nrow(res$alpha)

  # Now extract relevant rows from alpha matrix.
  if (length(include_idx) > 0)
    return(sqrt(colSums((res$alpha * res$mu2 -
                         (res$alpha*res$mu)^2)[include_idx,,drop=FALSE]))/
           (res$X_column_scale_factors))
  else
    return(numeric(ncol(res$mu)))
}

#' @rdname susie_get_methods
#'
#' @export
#'
susie_get_niter = function (res)
  res$niter

#' @rdname susie_get_methods
#'
#' @export
#'
susie_get_prior_variance = function (res)
  res$V

#' @rdname susie_get_methods
#'
#' @export
#'
susie_get_residual_variance = function (res)
  res$sigma2

#' @rdname susie_get_methods
#'
#' @importFrom stats pnorm
#'
#' @export
#'
susie_get_lfsr = function (res) {
  pos_prob = pnorm(0,mean = t(res$mu),sd = sqrt(res$mu2 - res$mu^2))
  neg_prob = 1 - pos_prob
  return(1 - rowSums(res$alpha * t(pmax(pos_prob,neg_prob))))
}

#' @rdname susie_get_methods
#'
#' @param susie_fit A susie fit, an output from \code{\link{susie}}.
#'
#' @param num_samples The number of draws from the posterior
#'   distribution.
#'
#' @importFrom stats rmultinom
#' @importFrom stats rnorm
#'
#' @export
#'
susie_get_posterior_samples = function (susie_fit, num_samples) {

  # Remove effects having estimated prior variance equals zero.
  if (is.numeric(susie_fit$V))
    include_idx = which(susie_fit$V > 1e-9)
  else
    include_idx = 1:nrow(susie_fit$alpha)

  posterior_mean = sweep(susie_fit$mu,2,susie_fit$X_column_scale_factors,"/")
  posterior_sd = sweep(sqrt(susie_fit$mu2 - (susie_fit$mu)^2),2,
                       susie_fit$X_column_scale_factors,"/")

  pip = susie_fit$alpha
  L = nrow(pip)
  num_snps = ncol(pip)
  b_samples = matrix(as.numeric(NA),num_snps,num_samples)
  gamma_samples = matrix(as.numeric(NA),num_snps,num_samples)
  for (sample_i in 1:num_samples) {
    b = 0
    if (length(include_idx) > 0) {
      for (l in include_idx) {
        gamma_l = rmultinom(1,1,pip[l,])
        effect_size = rnorm(1,mean = posterior_mean[l,which(gamma_l != 0)],
                            sd = posterior_sd[l,which(gamma_l != 0)])
        b_l = gamma_l * effect_size
        b = b + b_l
      }
    }
    b_samples[, sample_i] = b
    gamma_samples[, sample_i] = as.numeric(b != 0)
  }
  return(list(b = b_samples,gamma = gamma_samples))
}

#' @rdname susie_get_methods
#'
#' @param X n by p matrix of values of the p variables (covariates) in
#'   n samples. When provided, correlation between variables will be
#'   computed and used to remove CSs whose minimum correlation among
#'   variables is smaller than \code{min_abs_corr}.
#'
#' @param Xcorr p by p matrix of correlations between variables
#'   (covariates). When provided, it will be used to remove CSs whose
#'   minimum correlation among variables is smaller than
#'   \code{min_abs_corr}.
#'
#' @param coverage A number between 0 and 1 specifying desired
#'   coverage of each CS.
#'
#' @param min_abs_corr A "purity" threshold for the CS. Any CS that
#'   contains a pair of variables with correlation less than this
#'   threshold will be filtered out and not reported.
#'
#' @param dedup If \code{dedup = TRUE}, remove duplicate CSs.
#'
#' @param squared If \code{squared = TRUE}, report min, mean and
#' median of squared correlation instead of the absolute correlation.
#'
#' @param check_symmetric If \code{check_symmetric = TRUE}, perform a
#'   check for symmetry of matrix \code{Xcorr} when \code{Xcorr} is
#'   provided (not \code{NULL}).
#'
#' @param n_purity The maximum number of credible set (CS) variables
#'   used in calculating the correlation (\dQuote{purity})
#'   statistics. When the number of variables included in the CS is
#'   greater than this number, the CS variables are randomly subsampled.
#'
#' @param use_rfast Use the Rfast package for the purity calculations.
#'   By default \code{use_rfast = TRUE} if the Rfast package is
#'   installed.
#'
#'
#' @export
#'
susie_get_cs = function (res, X = NULL, Xcorr = NULL, coverage = 0.95,
                         min_abs_corr = 0.5, dedup = TRUE, squared = FALSE,
                         check_symmetric = TRUE, n_purity = 100, use_rfast) {
  if (!is.null(X) && !is.null(Xcorr))
    stop("Only one of X or Xcorr should be specified")
  if (check_symmetric) {
    if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr)) {
      warning_message("Xcorr is not symmetric; forcing Xcorr to be symmetric",
                  "by replacing Xcorr with (Xcorr + t(Xcorr))/2")
      Xcorr = Xcorr + t(Xcorr)
      Xcorr = Xcorr/2
    }
  }

  null_index = 0
  include_idx = rep(TRUE,nrow(res$alpha))
  if (!is.null(res$null_index)) null_index = res$null_index
  if (is.numeric(res$V)) include_idx = res$V > 1e-9
  # L x P binary matrix.
  status = in_CS(res$alpha,coverage)

  # L-list of CS positions.
  cs = lapply(1:nrow(status),function(i) which(status[i,]!=0))
  claimed_coverage = sapply(1:length(cs),
                            function (i) sum(res$alpha[i,][cs[[i]]]))
  include_idx = include_idx * (lapply(cs,length) > 0)

  # FIXME: see issue 21
  # https://github.com/stephenslab/susieR/issues/21
  if (dedup)
    include_idx = include_idx * (!duplicated(cs))
  include_idx = as.logical(include_idx)
  if (sum(include_idx) == 0)
    return(list(cs = NULL,
                coverage = NULL,
                requested_coverage = coverage))
  cs = cs[include_idx]
  claimed_coverage = claimed_coverage[include_idx]

  # Compute and filter by "purity".
  if (missing(use_rfast))
    use_rfast = requireNamespace("Rfast",quietly = TRUE)
  if (is.null(Xcorr) && is.null(X)) {
    names(cs) = paste0("L",which(include_idx))
    return(list(cs = cs,
                coverage = claimed_coverage,
                requested_coverage = coverage))
  } else {
    purity = NULL
    for (i in 1:length(cs)) {
      if (null_index > 0 && null_index %in% cs[[i]])
        purity = rbind(purity,c(-9,-9,-9))
      else
        purity =
          rbind(purity,
            matrix(get_purity(cs[[i]],X,Xcorr,squared,n_purity,use_rfast),1,3))
    }
    purity = as.data.frame(purity)
    if (squared)
      colnames(purity) = c("min.sq.corr","mean.sq.corr","median.sq.corr")
    else
      colnames(purity) = c("min.abs.corr","mean.abs.corr","median.abs.corr")
    threshold = ifelse(squared,min_abs_corr^2,min_abs_corr)
    is_pure = which(purity[,1] >= threshold)
    if (length(is_pure) > 0) {
      cs        = cs[is_pure]
      purity    = purity[is_pure,]
      row_names = paste0("L",which(include_idx)[is_pure])
      names(cs) = row_names
      rownames(purity) = row_names

      # Re-order CS list and purity rows based on purity.
      ordering = order(purity[,1],decreasing = TRUE)
      return(list(cs       = cs[ordering],
                  purity   = purity[ordering,],
                  cs_index = which(include_idx)[is_pure[ordering]],
                  coverage = claimed_coverage[ordering],
                  requested_coverage=coverage))
    } else
      return(list(cs = NULL,coverage = NULL,requested_coverage = coverage))
  }
}

#' @title Get Correlations Between CSs, using Variable with Maximum PIP From Each CS
#'
#' @description This function evaluates the correlation between single effect
#'   CSs. It is not part of the SuSiE inference. Rather, it is designed as
#'   a diagnostic tool to assess how correlated the reported CS are.
#'
#' @param model A SuSiE fit, typically an output from
#'   \code{\link{susie}} or one of its variants.
#'
#' @param X n by p matrix of values of the p variables (covariates) in
#'   n samples. When provided, correlation between variables will be
#'   computed and used to remove CSs whose minimum correlation among
#'   variables is smaller than \code{min_abs_corr}.
#'
#' @param Xcorr p by p matrix of correlations between variables
#'   (covariates). When provided, it will be used to remove CSs whose
#'   minimum correlation among variables is smaller than
#'   \code{min_abs_corr}.
#'
#' @param max When \code{max = FAFLSE}, return a matrix of CS
#'   correlations. When \code{max = TRUE}, return only the maximum
#'   absolute correlation among all pairs of correlations.
#'
#' @return A matrix of correlations between CSs, or the maximum
#'   absolute correlation when \code{max = TRUE}.
#'
#' @export
#'
get_cs_correlation = function (model, X = NULL, Xcorr = NULL, max = FALSE) {
  if (is.null(model$sets$cs) || length(model$sets$cs) == 1) return(NA)
  if (!is.null(X) && !is.null(Xcorr))
    stop("Only one of X or Xcorr should be specified")
  if (is.null(Xcorr) && is.null(X))
    stop("One of X or Xcorr must be specified")
  if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr)) {
    warning_message("Xcorr is not symmetric; forcing Xcorr to be symmetric",
                "by replacing Xcorr with (Xcorr + t(Xcorr))/2")
    Xcorr = Xcorr + t(Xcorr)
    Xcorr = Xcorr/2
  }
  # Get index for the best PIP per CS
  max_pip_idx = sapply(model$sets$cs, function(cs) cs[which.max(model$pip[cs])])
  if (is.null(Xcorr)) {
    X_sub = X[,max_pip_idx]
    cs_corr = muffled_corr(as.matrix(X_sub))
  } else {
    cs_corr = Xcorr[max_pip_idx, max_pip_idx]
  }
  if (max) {
    cs_corr = max(abs(cs_corr[upper.tri(cs_corr)]))
  }
  rownames(cs_corr) = colnames(cs_corr) = names(model$sets$cs)
  return(cs_corr)
}

#' @rdname susie_get_methods
#'
#' @param prune_by_cs Whether or not to ignore single effects not in
#'   a reported CS when calculating PIP.
#'
#' @param prior_tol Filter out effects having estimated prior variance
#'   smaller than this threshold.
#'
#' @export
#'
susie_get_pip = function (res, prune_by_cs = FALSE, prior_tol = 1e-9) {

  if (inherits(res,"susie")) {

    # Drop null weight columns.
    if (!is.null(res$null_index) && res$null_index > 0)
      res$alpha = res$alpha[,-res$null_index,drop=FALSE]

    # Drop the single-effects with estimated prior of zero.
    if (is.numeric(res$V))
      include_idx = which(res$V > prior_tol)
    else
      include_idx = 1:nrow(res$alpha)

    # Only consider variables in reported CS.
    # This is not what we do in the SuSiE paper.
    # So by default prune_by_cs = FALSE means we do not run the
    # following code.
    if (!is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = intersect(include_idx,res$sets$cs_index)
    if (is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = numeric(0)

    # now extract relevant rows from alpha matrix
    if (length(include_idx) > 0)
      res = res$alpha[include_idx,,drop = FALSE]
    else
      res = matrix(0,1,ncol(res$alpha))
  }

  return(as.vector(1 - apply(1 - res,2,prod)))
}

# Find how many variables in the CS.
# x is a probability vector.
#' @keywords internal
n_in_CS_x = function (x, coverage = 0.9)
  sum(cumsum(sort(x,decreasing = TRUE)) < coverage) + 1

# Return binary vector indicating if each point is in CS.
# x is a probability vector.
#' @keywords internal
in_CS_x = function (x, coverage = 0.9) {
  n = n_in_CS_x(x,coverage)
  o = order(x,decreasing = TRUE)
  result = rep(0,length(x))
  result[o[1:n]] = 1
  return(result)
}

# Returns an l-by-p binary matrix indicating which variables are in
# susie credible sets.
#' @keywords internal
in_CS = function (res, coverage = 0.9) {
  if (inherits(res,"susie"))
    res = res$alpha
  return(t(apply(res,1,function(x) in_CS_x(x,coverage))))
}

#' @keywords internal
n_in_CS = function(res, coverage = 0.9) {
  if (inherits(res,"susie"))
    res = res$alpha
  return(apply(res,1,function(x) n_in_CS_x(x,coverage)))
}

# Subsample and compute min, mean, median and max abs corr.
#
#' @importFrom stats median
get_purity = function (pos, X, Xcorr, squared = FALSE, n = 100,
                       use_rfast) {
  if (missing(use_rfast))
    use_rfast = requireNamespace("Rfast",quietly = TRUE)
  if (use_rfast) {
    get_upper_tri = Rfast::upper_tri
    get_median    = Rfast::med
  } else {
    get_upper_tri = function (R) R[upper.tri(R)]
    get_median    = stats::median
  }
  if (length(pos) == 1)
    return(c(1,1,1))
  else {

    # Subsample the columns if necessary.
    if (length(pos) > n)
      pos = sample(pos,n)

    if (is.null(Xcorr)) {
      X_sub = X[,pos]
      X_sub = as.matrix(X_sub)
      value = abs(get_upper_tri(muffled_corr(X_sub)))
    } else
      value = abs(get_upper_tri(Xcorr[pos,pos]))
    if (squared)
      value = value^2
    return(c(min(value),
             sum(value)/length(value),
             get_median(value)))
  }
}

# Correlation function with specified warning muffled.
#
#' @importFrom stats cor
muffled_corr = function (x)
  withCallingHandlers(cor(x),
                      warning = function(w) {
                        if (grepl("the standard deviation is zero",w$message))
                          invokeRestart("muffleWarning")
                      })

# cov2cor function with specified warning muffled.
#
#' @importFrom stats cov2cor
#' @keywords internal
muffled_cov2cor = function (x)
  withCallingHandlers(cov2cor(x),
    warning = function(w) {
      if (grepl("had 0 or NA entries; non-finite result is doubtful",
                w$message))
          invokeRestart("muffleWarning")
      })

# Check for symmetric matrix.
#' @keywords internal
is_symmetric_matrix = function (x) {
  if (requireNamespace("Rfast",quietly = TRUE))
    return(Rfast::is.symmetric(x))
  else
    return(isSymmetric(x))
}

# Compute standard error for regression coef.
# S = (X'X)^-1 \Sigma
#' @keywords internal
calc_stderr = function (X, residuals)
  sqrt(diag(sum(residuals^2)/(nrow(X) - 2) * chol2inv(chol(crossprod(X)))))

# Return residuals of Y after removing the linear effects of the susie
# model.
#
#' @importFrom stats coef
#' @keywords internal
get_R = function (X, Y, s)
  Y - X %*% coef(s)

# Slim the result of fitted susie model.
#' @keywords internal
susie_slim = function (res)
  list(alpha = res$alpha,niter = res$niter,V = res$V,sigma2 = res$sigma2)

# Prune single effects to given number L in susie model object.
#' @keywords internal
susie_prune_single_effects = function (s,L = 0,V = NULL) {
  num_effects = nrow(s$alpha)
  if (L == 0) {
    # Filtering will be based on non-zero elements in s$V.
    if (!is.null(s$V))
      L = length(which(s$V > 0))
    else
      L = num_effects
  }
  if (L == num_effects) {
    s$sets = NULL
    return(s)
  }
  if (!is.null(s$sets$cs_index))
    effects_rank = c(s$sets$cs_index,setdiff(1:num_effects,s$sets$cs_index))
  else
    effects_rank = 1:num_effects

  if (L > num_effects) {
    message(paste("Specified number of effects L =",L,
                  "is greater the number of effects",num_effects,
                  "in input SuSiE model. The SuSiE model will be expanded",
                  "to have",L,"effects."))

    s$alpha = rbind(s$alpha[effects_rank,],
                    matrix(1/ncol(s$alpha),L - num_effects,ncol(s$alpha)))
    for (n in c("mu","mu2","lbf_variable"))
      if (!is.null(s[[n]]))
        s[[n]] = rbind(s[[n]][effects_rank,],
                       matrix(0,L - num_effects,ncol(s[[n]])))
    for (n in c("KL", "lbf"))
      if (!is.null(s[[n]]))
        s[[n]] = c(s[[n]][effects_rank],rep(NA, L-num_effects))
    if (!is.null(V)) {
      if (length(V) > 1)
        V[1:num_effects] = s$V[effects_rank]
      else V = rep(V,L)
    }
    s$V = V
  }
  s$sets = NULL
  return(s)
}

#' @title Estimate s in \code{susie_rss} Model Using Regularized LD
#'
#' @description The estimated s gives information about the
#'   consistency between the z scores and LD matrix. A larger \eqn{s}
#'   means there is a strong inconsistency between z scores and LD
#'   matrix. The \dQuote{null-mle} method obtains mle of \eqn{s} under
#'   \eqn{z | R ~ N(0,(1-s)R + s I)}, \eqn{0 < s < 1}. The
#'   \dQuote{null-partialmle} method obtains mle of \eqn{s} under
#'   \eqn{U^T z | R ~ N(0,s I)}, in which \eqn{U} is a matrix containing
#'   the of eigenvectors that span the null space of R; that is, the
#'   eigenvectors corresponding to zero eigenvalues of R. The estimated
#'   \eqn{s} from \dQuote{null-partialmle} could be greater than 1. The
#'   \dQuote{null-pseudomle} method obtains mle of \eqn{s} under
#'   pseudolikelihood \eqn{L(s) = \prod_{j=1}^{p} p(z_j | z_{-j}, s,
#'   R)}, \eqn{0 < s < 1}.
#'
#' @param z A p-vector of z scores.
#'
#' @param R A p by p symmetric, positive semidefinite correlation
#'   matrix.
#'
#' @param n The sample size. (Optional, but highly recommended.)
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param method a string specifies the method to estimate \eqn{s}.
#'
#' @return A number between 0 and 1.
#'
#' @examples
#' set.seed(1)
#' n = 500
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 0.01
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#' input_ss = compute_suff_stat(X,y,standardize = TRUE)
#' ss = univariate_regression(X,y)
#' R = cor(X)
#' attr(R,"eigen") = eigen(R,symmetric = TRUE)
#' zhat = with(ss,betahat/sebetahat)
#'
#' # Estimate s using the unadjusted z-scores.
#' s0 = estimate_s_rss(zhat,R)
#' 
#' # Estimate s using the adjusted z-scores.
#' s1 = estimate_s_rss(zhat,R,n)
#'
#' @importFrom stats dnorm
#' @importFrom stats optim
#'
#' @export
#'
estimate_s_rss = function (z, R, n, r_tol = 1e-08, method = "null-mle") {

  # Check and process input arguments z, R.
  z[is.na(z)] = 0
  if (is.null(attr(R,"eigen")))
    attr(R,"eigen") = eigen(R,symmetric = TRUE)
  eigenld = attr(R,"eigen")
  if (any(eigenld$values < -r_tol))
    warning_message("The matrix R is not positive semidefinite. Negative ",
            "eigenvalues are set to zero")
  eigenld$values[eigenld$values < r_tol] = 0

  # Check input n, and adjust the z-scores if n is provided.
  if (missing(n))
    warning_message("Providing the sample size (n), or even a rough estimate of n, ",
            "is highly recommended. Without n, the implicit assumption is ",
            "n is large (Inf) and the effect sizes are small (close to zero).")
  else if (n <= 1)
    stop("n must be greater than 1")
  if (!missing(n)) {
    sigma2 = (n-1)/(z^2 + n - 2)
    z = sqrt(sigma2) * z
  }

  if (method == "null-mle") {
    negloglikelihood = function(s, ztv, d) {
      0.5 * sum(log((1 - s) * d + s)) +
      0.5 * tcrossprod(ztv / ((1 - s) * d + s), ztv)
    }
    s = optim(0.5,fn = negloglikelihood,ztv = crossprod(z, eigenld$vectors),
              d = eigenld$values,method = "Brent",lower = 0,upper = 1)$par
  } else if (method == "null-partialmle") {
    colspace = which(eigenld$values > 0)
    if (length(colspace) == length(z))
      s = 0
    else{
      znull = crossprod(eigenld$vectors[,-colspace], z) # U2^T z
      s = sum(znull^2)/length(znull)
    }
  } else if (method == "null-pseudomle") {
    pseudolikelihood = function(s, z, eigenld) {
      precision = eigenld$vectors %*% (t(eigenld$vectors) *
                  (1/((1-s)*eigenld$values + s)))
      postmean = rep(0,length(z))
      postvar = rep(0,length(z))
      for (i in 1:length(z)) {
        postmean[i] = -(1/precision[i,i])*precision[i,-i] %*% z[-i]
        postvar[i] = 1/precision[i,i]
      }
      return(-sum(dnorm(z,mean = postmean,sd = sqrt(postvar),log = TRUE)))
    }
    s = optim(0.5,fn = pseudolikelihood,z = z,eigenld = eigenld,
              method = "Brent",lower = 0,upper = 1)$par
  }
  else
    stop("The method is not implemented")
  return(s)
}

#' @title Compute Distribution of z-scores of Variant j Given Other z-scores, and Detect Possible Allele Switch Issue
#'
#' @description Under the null, the rss model with regularized LD
#'   matrix is \eqn{z|R,s ~ N(0, (1-s)R + s I))}. We use a mixture of
#'   normals to model the conditional distribution of z_j given other z
#'   scores, \eqn{z_j | z_{-j}, R, s ~ \sum_{k=1}^{K} \pi_k
#'   N(-\Omega_{j,-j} z_{-j}/\Omega_{jj}, \sigma_{k}^2/\Omega_{jj})},
#'   \eqn{\Omega = ((1-s)R + sI)^{-1}}, \eqn{\sigma_1, ..., \sigma_k}
#'   is a grid of fixed positive numbers. We estimate the mixture
#'   weights \eqn{\pi}  We detect the possible allele switch issue
#'   using likelihood ratio for each variant.
#'
#' @param z A p-vector of z scores.
#'
#' @param R A p by p symmetric, positive semidefinite correlation
#'   matrix.
#'
#' @param n The sample size. (Optional, but highly recommended.)
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param s an estimated s from \code{estimate_s_rss}
#'
#' @return a list containing a ggplot2 plot object and a table. The plot
#'   compares observed z score vs the expected value. The possible allele
#'   switched variants are labeled as red points (log LR > 2 and abs(z) > 2).
#'   The table summarizes the conditional distribution for each variant
#'   and the likelihood ratio test. The table has the following columns:
#'   the observed z scores, the conditional expectation, the conditional
#'   variance, the standardized differences between the observed z score
#'   and expected value, the log likelihood ratio statistics.
#'
#' @importFrom stats dnorm
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes_string
#' @importFrom mixsqp mixsqp
#'
#' @examples
#' # See also the vignette, "Diagnostic for fine-mapping with summary
#' # statistics."
#' set.seed(1)
#' n = 500
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 0.01
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#' ss = univariate_regression(X,y)
#' R = cor(X)
#' attr(R,"eigen") = eigen(R,symmetric = TRUE)
#' zhat = with(ss,betahat/sebetahat)
#' cond_dist = kriging_rss(zhat,R,n = n)
#' cond_dist$plot
#'
#' @export
#'
kriging_rss = function (z, R, n, r_tol = 1e-08,
                        s = estimate_s_rss(z,R,n,r_tol,method = "null-mle")) {

  # Check and process input arguments z, R.
  z[is.na(z)] = 0
  if (is.null(attr(R,"eigen")))
    attr(R,"eigen") = eigen(R,symmetric = TRUE)
  eigenld = attr(R,"eigen")
  if (any(eigenld$values < -r_tol))
    warning_message("The matrix R is not positive semidefinite. Negative ",
            "eigenvalues are set to zero.")
  eigenld$values[eigenld$values < r_tol] = 0

  # Check and progress input argument s.
  force(s)
  if (s > 1) {
    warning_message("The given s is greater than 1. We replace it with 0.8.")
    s = 0.8
  } else if (s < 0)
    stop("The s must be non-negative")
  
  # Check input n, and adjust the z-scores if n is provided.
  if ((!missing(n)) && (n <= 1))
    stop("n must be greater than 1")
  if (missing(n))
    warning_message("Providing the sample size (n), or even a rough estimate of n, ",
            "is highly recommended. Without n, the implicit assumption is ",
            "n is large (Inf) and the effect sizes are small (close to zero).")
  else {
    sigma2 = (n-1)/(z^2 + n - 2)
    z = sqrt(sigma2) * z
  }
  
  dinv = 1/((1-s) * eigenld$values + s)
  dinv[is.infinite(dinv)] = 0
  precision = eigenld$vectors %*% (t(eigenld$vectors) * dinv)
  condmean = rep(0,length(z))
  condvar = rep(0,length(z))
  for (i in 1:length(z)) {
    condmean[i] = -(1/precision[i,i]) * precision[i,-i] %*% z[-i]
    condvar[i]  = 1/precision[i,i]
  }
  z_std_diff = (z-condmean)/sqrt(condvar)

  # obtain grid
  a_min = 0.8
  if (max(z_std_diff^2) < 1)
    a_max = 2
  else
    a_max = 2*sqrt(max(z_std_diff^2))
  npoint = ceiling(log2(a_max/a_min)/log2(1.05))
  a_grid = 1.05^(seq(-npoint,0)) * a_max

  # compute likelihood
  sd_mtx      = outer(sqrt(condvar),a_grid)
  matrix_llik = dnorm(z - condmean,sd = sd_mtx,log = TRUE)
  lfactors    = apply(matrix_llik,1,max)
  matrix_llik = matrix_llik - lfactors

  # estimate weight
  w = mixsqp(matrix_llik,log = TRUE,control = list(verbose = FALSE))$x

  # Compute denominators in likelihood ratios.
  logl0mix = drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors

  # Compute numerators in likelihood ratios.
  matrix_llik = dnorm(z + condmean,sd = sd_mtx,log = TRUE)
  lfactors    = apply(matrix_llik,1,max)
  matrix_llik = matrix_llik - lfactors
  logl1mix    = drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors

  # Compute (log) likelihood ratios.
  logLRmix = logl1mix - logl0mix

  z          = drop(z)
  z_std_diff = drop(z_std_diff)
  res = data.frame(z = z,
                   condmean = condmean,
                   condvar = condvar,
                   z_std_diff = z_std_diff,
                   logLR = logLRmix)
  p = ggplot(res,aes_string(y = "z",x = "condmean")) +
        geom_point() +
        labs(y = "Observed z scores", x = "Expected value") +
        geom_abline(intercept = 0, slope = 1) +
        theme_bw()
  idx = which(logLRmix > 2 & abs(z) > 2)
  if (length(idx) > 0)
    p = p + geom_point(data = res[idx,],
                       aes_string(y = "z", x = "condmean"),col = "red")
  return(list(plot = p,conditional_dist = res))
}

# Compute the column means of X, the column standard deviations of X,
# and rowSums(Y^2), where Y is the centered and/or scaled version of
# X.
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
compute_colstats = function (X, center = TRUE, scale = TRUE) {
  n = nrow(X)
  p = ncol(X)
  if (!is.null(attr(X,"matrix.type"))) {

    # X is a trend filtering matrix.
    cm  = compute_tf_cm(attr(X,"order"),p)
    csd = compute_tf_csd(attr(X,"order"),p)
    d   = compute_tf_d(attr(X,"order"),p,cm,csd,scale,center)
    if (!center)
      cm = rep(0,p)
    if (!scale)
      csd = rep(1,p)
  } else {

    # X is an ordinary dense or sparse matrix. Set sd = 1 when the
    # column has variance 0.
    if (center)
      cm = colMeans(X,na.rm = TRUE)
    else
      cm = rep(0,p)
    if (scale) {
      csd = compute_colSds(X)
      csd[csd == 0] = 1
    } else
      csd = rep(1,p)

    # These two lines of code should give the same result as
    #
    #   Y = (t(X) - cm)/csd
    #   d = rowSums(Y^2)
    #
    # for all four combinations of "center" and "scale", but do so
    # without having to modify X, or create copies of X in memory. In
    # particular the first line should be equivalent to colSums(X^2).
    d = n*colMeans(X)^2 + (n-1)*compute_colSds(X)^2
    d = (d - n*cm^2)/csd^2
  }

  return(list(cm = cm,csd = csd,d = d))
}

# @title computes column standard deviations for any type of matrix
# @details This should give the same result as matrixStats::colSds(X),
#   but allows for sparse matrices as well as dense ones.
# @param X an n by p matrix of any type, e.g. sparse, dense.
# @return a p vector of column standard deviations.
#
#' @importFrom matrixStats colSds
#' @importFrom Matrix summary
compute_colSds = function(X) {
  if (is.matrix(X))
    return(colSds(X))
  else {
    n = nrow(X)
    Y = apply_nonzeros(X,function (u) u^2)
    d = colMeans(Y) - colMeans(X)^2
    return(sqrt(d*n/(n-1)))
  }
}

# @title Check whether A is positive semidefinite
# @param A a symmetric matrix
# @return a list of result:
# \item{matrix}{The matrix with eigen decomposition}
# \item{status}{whether A is positive semidefinite}
# \item{eigenvalues}{eigenvalues of A truncated by r_tol}
check_semi_pd = function (A, tol) {
  attr(A,"eigen") = eigen(A,symmetric = TRUE)
  v = attr(A,"eigen")$values
  v[abs(v) < tol] = 0
  return(list(matrix      = A,
              status      = !any(v < 0),
              eigenvalues = v))
}

# @title Check whether b is in space spanned by the non-zero eigenvectors
#   of A
# @param A a p by p matrix
# @param b a length p vector
# @return a list of result:
# \item{status}{whether b in space spanned by the non-zero
#  eigenvectors of A}
# \item{msg}{msg gives the difference between the projected b and b if
#   status is FALSE}
check_projection = function (A, b) {
  if (is.null(attr(A,"eigen")))
    attr(A,"eigen") = eigen(A,symmetric = TRUE)
  v = attr(A,"eigen")$values
  B = attr(A,"eigen")$vectors[,v > .Machine$double.eps]
  msg = all.equal(as.vector(B %*% crossprod(B,b)),as.vector(b),
                  check.names = FALSE)
  if (!is.character(msg))
    return(list(status = TRUE,msg = NA))
  else
    return(list(status = FALSE,msg = msg))
}

# @title Utility function to display warning messages as they occur
# @param ... warning message
# @param style either "warning" or "hint"
#'@importFrom crayon combine_styles
warning_message = function(..., style=c("warning", "hint")) {
  style = match.arg(style)
  if (style=="warning" && getOption("warn")>=0) {
    alert <- combine_styles("bold", "underline", "red")	
    message(alert("WARNING:"), " ", ...)
  } else {
    alert <- combine_styles("bold", "underline", "magenta")	
    message(alert("HINT:"), " ", ...)
  }
}

# Apply operation f to all nonzeros of a sparse matrix.
#
#' @importFrom Matrix sparseMatrix
#' 
apply_nonzeros <- function (X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i,j = d$j,x = f(d$x),dims = dim(X)))
}
