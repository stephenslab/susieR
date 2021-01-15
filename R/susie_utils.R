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
    return(colSums((res$alpha*res$mu)[include_idx,])/
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
                         (res$alpha*res$mu)^2)[include_idx,]))/
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
susie_get_posterior_samples <- function (susie_fit, num_samples) {
    
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
#' @export
#'
susie_get_cs = function (res, X = NULL, Xcorr = NULL, coverage = 0.95,
                         min_abs_corr = 0.5, dedup = TRUE, squared = FALSE) {
  if (!is.null(X) && !is.null(Xcorr))
    stop("Only one of X or Xcorr should be specified")
  if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr))
    stop("Xcorr matrix must be symmetric")
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
  if (is.null(Xcorr) && is.null(X)) {
    names(cs) = paste0("L",which(include_idx))
    return(list(cs = cs,
                coverage = claimed_coverage,
                requested_coverage = coverage))
  } else {
    purity = data.frame(do.call(rbind,lapply(1:length(cs),function (i) {
              if (null_index > 0 && null_index %in% cs[[i]])
                c(-9,-9,-9)
              else
                get_purity(cs[[i]],X,Xcorr,squared)
             })))
    if (squared)
      colnames(purity) = c("min.sq.corr","mean.sq.corr","median.sq.corr")
    else
      colnames(purity) = c("min.abs.corr","mean.abs.corr","median.abs.corr")
    threshold = ifelse(squared,min_abs_corr^2,min_abs_corr)
    is_pure = which(purity[,1] >= threshold)
    if (length(is_pure) > 0) {
      cs = cs[is_pure]
      purity = purity[is_pure,]
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
      return(list(cs = NULL,coverage = NULL, requested_coverage = coverage))
  }
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
get_purity = function(pos, X, Xcorr, squared = FALSE, n = 100) {
  if (length(pos) == 1)
    c(1,1,1)
  else {
    if (length(pos) > n)
      pos = sample(pos, n)
    if (is.null(Xcorr)) {
      X_sub = X[,pos]
      if (length(pos) > n) {

        # Remove identical columns.
        pos_rm = sapply(1:ncol(X_sub),
                       function(i) all(abs(X_sub[,i] - mean(X_sub[,i])) <
                                       .Machine$double.eps^0.5))
        if (length(pos_rm))
          X_sub = X_sub[,-pos_rm]
      }
      value = abs(muffled_corr(as.matrix(X_sub)))
    } else
      value = abs(Xcorr[pos,pos])
    if (squared)
      value = value^2
    return(c(min(value,na.rm = TRUE),
             mean(value,na.rm = TRUE),
             median(value,na.rm = TRUE)))
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
  res = isSymmetric(x)
  if (!res)
    res = isSymmetric(unname(x))
  return(res)
}

# Compute standard error for regression coef.
# S = (X'X)^-1 \Sigma
#' @keywords internal
calc_stderr = function (X, residuals)
  sqrt(diag(sum(residuals^2)/(nrow(X) - 2) * chol2inv(chol(t(X) %*% X))))

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
susie_prune_single_effects = function (s,L = 0,V = NULL,verbose = FALSE) {
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
  if (verbose) 
    warning(paste("Specified number of effects L =",L,
                  "does not match the number of effects",num_effects,
                  "in input SuSiE model. It will be",
                  ifelse(L < num_effects,"pruned","expanded"),"to have",L,
                  "effects."))
  if (L < num_effects) {
    for (n in c("alpha","mu","mu2","lbf_variable"))
      if (!is.null(s[[n]]))
        s[[n]] = s[[n]][effects_rank,][1:L,]
    for (n in c("KL","lbf","V"))
      if (!is.null(s[[n]]))
        s[[n]] = s[[n]][effects_rank][1:L]
  } else {
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
