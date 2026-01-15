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
#' n <- 1000
#' p <- 1000
#' beta <- rep(0, p)
#' beta[1:4] <- 1
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' X <- scale(X, center = TRUE, scale = TRUE)
#' y <- drop(X %*% beta + rnorm(n))
#' s <- susie(X, y, L = 10)
#' susie_get_objective(s)
#' susie_get_objective(s, last_only = FALSE)
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
susie_get_objective <- function(res, last_only = TRUE, warning_tol = 1e-6) {
  if (!all(diff(res$elbo) >= (-1 * warning_tol))) {
    warning_message("Objective is decreasing")
  }
  if (last_only) {
    return(res$elbo[length(res$elbo)])
  } else {
    return(res$elbo)
  }
}

#' @rdname susie_get_methods
#'
#' @export
#'
susie_get_posterior_mean <- function(res, prior_tol = 1e-9) {
  # Drop the single-effects with estimated prior of zero.
  if (is.numeric(res$V)) {
    include_idx <- which(res$V > prior_tol)
  } else {
    include_idx <- 1:nrow(res$alpha)
  }

  # Now extract relevant rows from alpha matrix.
  if (length(include_idx) > 0) {
    return(colSums((res$alpha * res$mu)[include_idx, , drop = FALSE]) /
      res$X_column_scale_factors)
  } else {
    return(numeric(ncol(res$mu)))
  }
}

#' @rdname susie_get_methods
#'
#' @export
#'
susie_get_posterior_sd <- function(res, prior_tol = 1e-9) {
  # Drop the single-effects with estimated prior of zero.
  if (is.numeric(res$V)) {
    include_idx <- which(res$V > prior_tol)
  } else {
    include_idx <- 1:nrow(res$alpha)
  }

  # Now extract relevant rows from alpha matrix.
  if (length(include_idx) > 0) {
    return(sqrt(colSums((res$alpha * res$mu2 -
      (res$alpha * res$mu)^2)[include_idx, , drop = FALSE])) /
      (res$X_column_scale_factors))
  } else {
    return(numeric(ncol(res$mu)))
  }
}

#' @rdname susie_get_methods
#'
#' @export
#'
susie_get_niter <- function(res) {
  res$niter
}

#' @rdname susie_get_methods
#'
#' @export
#'
susie_get_prior_variance <- function(res) {
  res$V
}

#' @rdname susie_get_methods
#'
#' @export
#'
susie_get_residual_variance <- function(res) {
  res$sigma2
}

#' @rdname susie_get_methods
#'
#' @importFrom stats pnorm
#'
#' @export
#'
susie_get_lfsr <- function(res) {
  pos_prob <- pnorm(0, mean = t(res$mu), sd = sqrt(res$mu2 - res$mu^2))
  neg_prob <- 1 - pos_prob
  return(1 - rowSums(res$alpha * t(pmax(pos_prob, neg_prob))))
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
susie_get_posterior_samples <- function(susie_fit, num_samples) {
  # Remove effects having estimated prior variance equals zero.
  if (is.numeric(susie_fit$V)) {
    include_idx <- which(susie_fit$V > 1e-9)
  } else {
    include_idx <- 1:nrow(susie_fit$alpha)
  }

  posterior_mean <- sweep(susie_fit$mu, 2, susie_fit$X_column_scale_factors, "/")
  posterior_sd <- sweep(
    sqrt(susie_fit$mu2 - (susie_fit$mu)^2), 2,
    susie_fit$X_column_scale_factors, "/"
  )

  pip <- susie_fit$alpha
  L <- nrow(pip)
  num_snps <- ncol(pip)
  b_samples <- matrix(as.numeric(NA), num_snps, num_samples)
  gamma_samples <- matrix(as.numeric(NA), num_snps, num_samples)
  for (sample_i in 1:num_samples) {
    b <- 0
    if (length(include_idx) > 0) {
      for (l in include_idx) {
        gamma_l <- rmultinom(1, 1, pip[l, ])
        effect_size <- rnorm(1,
          mean = posterior_mean[l, which(gamma_l != 0)],
          sd = posterior_sd[l, which(gamma_l != 0)]
        )
        b_l <- gamma_l * effect_size
        b <- b + b_l
      }
    }
    b_samples[, sample_i] <- b
    gamma_samples[, sample_i] <- as.numeric(b != 0)
  }
  return(list(b = b_samples, gamma = gamma_samples))
}

#' @rdname susie_get_methods
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
#'   median of squared correlation instead of the absolute correlation.
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
#' @param ld_extend_threshold Threshold for extending CS by LD (default 0.99).
#'   Variants with |correlation| > threshold with any CS member are added.
#'   Set to NULL to disable LD extension. Requires Xcorr (would not work if only X is provided).
#'
#' @export
#'
susie_get_cs <- function(res, X = NULL, Xcorr = NULL, coverage = 0.95,
                         min_abs_corr = 0.5, dedup = TRUE, squared = FALSE,
                         check_symmetric = TRUE, n_purity = 100, use_rfast,
                         ld_extend_threshold = 0.99) {
  if (!is.null(X) && !is.null(Xcorr)) {
    stop("Only one of X or Xcorr should be specified")
  }
  if (check_symmetric) {
    if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr)) {
      warning_message(
        "Xcorr is not symmetric; forcing Xcorr to be symmetric",
        "by replacing Xcorr with (Xcorr + t(Xcorr))/2"
      )
      Xcorr <- Xcorr + t(Xcorr)
      Xcorr <- Xcorr / 2
    }
  }

  null_index <- 0
  include_idx <- rep(TRUE, nrow(res$alpha))
  if (!is.null(res$null_index)) null_index <- res$null_index
  if (is.numeric(res$V)) include_idx <- res$V > 1e-9
  
  # L x P binary matrix
  status <- in_CS(res$alpha, coverage)

  # L-list of CS positions
  cs <- lapply(1:nrow(status), function(i) which(status[i, ] != 0))
  claimed_coverage <- sapply(
    1:length(cs),
    function(i) sum(res$alpha[i, ][cs[[i]]])
  )
  include_idx <- include_idx * (lapply(cs, length) > 0)

  # FIXME: see issue 21
  # https://github.com/stephenslab/susieR/issues/21
  if (dedup) {
    include_idx <- include_idx * (!duplicated(cs))
  }
  include_idx <- as.logical(include_idx)
  if (sum(include_idx) == 0) {
    return(list(
      cs = NULL,
      coverage = NULL,
      requested_coverage = coverage
    ))
  }
  
  cs <- cs[include_idx]
  claimed_coverage <- claimed_coverage[include_idx]
  # Track which original effects these correspond to
  effect_indices <- which(include_idx)

  # Compute and filter by "purity"
  if (missing(use_rfast)) {
    use_rfast <- requireNamespace("Rfast", quietly = TRUE)
  }
  
  # If no correlation info, return without purity or LD extension
  if (is.null(Xcorr) && is.null(X)) {
    names(cs) <- paste0("L", effect_indices)
    return(list(
      cs = cs,
      coverage = claimed_coverage,
      requested_coverage = coverage
    ))
  }
  
  # Extend CS by LD if threshold is set and Xcorr is available
  # Note: LD extension requires Xcorr; if only X is provided, skip extension
  # (X may be sparse, and computing full Xcorr is expensive/infeasible)
  if (!is.null(ld_extend_threshold) && !is.null(Xcorr)) {
    for (i in 1:length(cs)) {
      cs_idx <- cs[[i]]
      # Find variants in tight LD with any CS member
      ld_with_cs <- abs(Xcorr[cs_idx, , drop = FALSE]) > ld_extend_threshold
      in_tight_ld <- which(colSums(ld_with_cs) > 0)
      # Extend CS
      cs[[i]] <- sort(unique(c(cs_idx, in_tight_ld)))
      # Update coverage for extended CS
      claimed_coverage[i] <- sum(res$alpha[effect_indices[i], cs[[i]]])
    }
  }

  # Compute purity for each CS
  purity <- NULL
  for (i in 1:length(cs)) {
    if (null_index > 0 && null_index %in% cs[[i]]) {
      purity <- rbind(purity, c(-9, -9, -9))
    } else {
      purity <- rbind(
        purity,
        matrix(get_purity(cs[[i]], X, Xcorr, squared, n_purity, use_rfast), 1, 3)
      )
    }
  }
  purity <- as.data.frame(purity)
  if (squared) {
    colnames(purity) <- c("min.sq.corr", "mean.sq.corr", "median.sq.corr")
  } else {
    colnames(purity) <- c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
  }
  
  threshold <- ifelse(squared, min_abs_corr^2, min_abs_corr)
  is_pure <- which(purity[, 1] >= threshold)
  
  if (length(is_pure) > 0) {
    cs <- cs[is_pure]
    purity <- purity[is_pure, , drop = FALSE]
    claimed_coverage <- claimed_coverage[is_pure]
    effect_indices <- effect_indices[is_pure]
    
    row_names <- paste0("L", effect_indices)
    names(cs) <- row_names
    rownames(purity) <- row_names
    
    # Re-order CS list and purity rows based on purity
    ordering <- order(purity[, 1], decreasing = TRUE)
    return(list(
      cs = cs[ordering],
      purity = purity[ordering, , drop = FALSE],
      cs_index = effect_indices[ordering],
      coverage = claimed_coverage[ordering],
      requested_coverage = coverage
    ))
  } else {
    return(list(cs = NULL, coverage = NULL, requested_coverage = coverage))
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
get_cs_correlation <- function(model, X = NULL, Xcorr = NULL, max = FALSE) {
  if (is.null(model$sets$cs) || length(model$sets$cs) == 1) {
    return(NA)
  }
  if (!is.null(X) && !is.null(Xcorr)) {
    stop("Only one of X or Xcorr should be specified")
  }
  if (is.null(Xcorr) && is.null(X)) {
    stop("One of X or Xcorr must be specified")
  }
  if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr)) {
    warning_message(
      "Xcorr is not symmetric; forcing Xcorr to be symmetric",
      "by replacing Xcorr with (Xcorr + t(Xcorr))/2"
    )
    Xcorr <- Xcorr + t(Xcorr)
    Xcorr <- Xcorr / 2
  }
  # Get index for the best PIP per CS
  max_pip_idx <- sapply(model$sets$cs, function(cs) cs[which.max(model$pip[cs])])
  if (is.null(Xcorr)) {
    X_sub <- X[, max_pip_idx]
    cs_corr <- safe_cor(as.matrix(X_sub))
  } else {
    cs_corr <- Xcorr[max_pip_idx, max_pip_idx]
  }
  if (max) {
    cs_corr <- max(abs(cs_corr[upper.tri(cs_corr)]))
  } else {
    rownames(cs_corr) <- colnames(cs_corr) <- names(model$sets$cs)
  }
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
susie_get_pip <- function(res, prune_by_cs = FALSE, prior_tol = 1e-9) {
  if (inherits(res, "susie")) {
    # Drop null weight columns.
    if (!is.null(res$null_index) && res$null_index > 0) {
      res$alpha <- res$alpha[, -res$null_index, drop = FALSE]
    }

    # Drop the single-effects with estimated prior of zero.
    if (is.numeric(res$V)) {
      include_idx <- which(res$V > prior_tol)
    } else {
      include_idx <- 1:nrow(res$alpha)
    }

    # Only consider variables in reported CS.
    # This is not what we do in the SuSiE paper.
    # So by default prune_by_cs = FALSE means we do not run the
    # following code.
    if (!is.null(res$sets$cs_index) && prune_by_cs) {
      include_idx <- intersect(include_idx, res$sets$cs_index)
    }
    if (is.null(res$sets$cs_index) && prune_by_cs) {
      include_idx <- numeric(0)
    }

    # now extract relevant rows from alpha matrix
    if (length(include_idx) > 0) {
      res <- res$alpha[include_idx, , drop = FALSE]
    } else {
      res <- matrix(0, 1, ncol(res$alpha))
    }
  }

  return(as.vector(1 - apply(1 - res, 2, prod)))
}

#' @title Initialize a susie object using regression coefficients
#'
#' @param coef_index An L-vector containing the the indices of the
#'   nonzero coefficients.
#'
#' @param coef_value An L-vector containing initial coefficient
#' estimates.
#'
#' @param p A scalar giving the number of variables.
#'
#' @return A list with elements \code{alpha}, \code{mu} and \code{mu2}
#'   to be used by \code{susie}.
#'
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[sample(1:1000,4)] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#'
#' # Initialize susie to ground-truth coefficients.
#' s = susie_init_coef(which(beta != 0),beta[beta != 0],length(beta))
#' res = susie(X,y,L = 10,model_init=s)
#'
#' @export
#'
susie_init_coef = function (coef_index, coef_value, p) {
  L = length(coef_index)
  if (L <= 0)
    stop("Need at least one non-zero effect")
  if (!all(coef_value != 0))
    stop("Input coef_value must be non-zero for all its elements")
  if (L != length(coef_value))
    stop("Inputs coef_index and coef_value must of the same length")
  if (max(coef_index) > p)
    stop("Input coef_index exceeds the boundary of p")
  alpha = matrix(0,nrow = L,ncol = p)
  mu = matrix(0,nrow = L,ncol = p)
  for(i in 1:L){
    alpha[i,coef_index[i]] = 1
    mu[i,coef_index[i]] = coef_value[i]
  }
  # Add prior variances V field (required by new architecture)
  V = rep(0.2, L)  # Default scaled prior variance for each effect
  out = list(alpha = alpha, mu = mu, mu2 = mu*mu, V = V)
  class(out) = c("susie","list")
  return(out)
}
