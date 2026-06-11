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
#' @param min_abs_corr A "purity" threshold for the CS, applied to the
#'   minimum absolute correlation among the CS variables: a CS is filtered
#'   out unless this minimum is at least \code{min_abs_corr}. Set to
#'   \code{NULL} to disable this clause. This filter is only applied when
#'   \code{X} or \code{Xcorr} is provided; otherwise it is ignored and a
#'   warning is issued.
#'
#' @param median_abs_corr An optional second purity threshold applied to
#'   the \emph{median} absolute correlation among the CS variables. The
#'   default, \code{NULL}, leaves it off. When both \code{min_abs_corr} and
#'   \code{median_abs_corr} are non-\code{NULL}, they are combined with OR:
#'   a CS is kept if it clears \emph{either} threshold (its min
#'   \eqn{\ge} \code{min_abs_corr} \strong{or} its median \eqn{\ge}
#'   \code{median_abs_corr}).
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
#' @param n_purity Maximum number of CS variables used for purity when
#'   correlations are computed from \code{X}. The default, \code{"auto"}, uses
#'   a resource-aware cap; a positive number gives a fixed cap; a negative
#'   number uses all CS variables. \code{Xcorr} inputs always use all variables.
#'
#' @param use_rfast Use the Rfast package for the purity calculations.
#'   By default \code{use_rfast = TRUE} if the Rfast package is
#'   installed.
#'
#' @param cs_extension_corr Either \code{NULL} or a single number between 0
#'   and 1. If non-\code{NULL}, each credible set is extended to include every
#'   variable whose absolute correlation with a credible-set member exceeds this
#'   threshold, pulling in near-perfectly correlated proxies (variables
#'   statistically indistinguishable from the selected ones). Works from either
#'   \code{X} or \code{Xcorr}. The default, \code{NULL}, disables correlation-based
#'   extension; \code{0.99} is the recommended value when extension is desired.
#'
#' @export
#'
susie_get_cs <- function(res, X = NULL, Xcorr = NULL, coverage = 0.95,
                         min_abs_corr = 0.5, dedup = TRUE, squared = FALSE,
                         check_symmetric = TRUE, n_purity = "auto",
                         use_rfast = NULL, median_abs_corr = NULL,
                         cs_extension_corr = NULL) {
  if (!is.null(X) && !is.null(Xcorr)) {
    stop("Only one of X or Xcorr should be specified")
  }
  for (nm in c("min_abs_corr", "median_abs_corr", "cs_extension_corr")) {
    v <- get(nm)
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1 || !is.finite(v) ||
                        v < 0 || v > 1)) {
      stop(nm, " must be NULL or a single numeric value in [0, 1].")
    }
  }
  if (is.null(X) && is.null(Xcorr)) {
    warning_message(
      "Neither X nor Xcorr was provided; purity filtering is skipped ",
      "and min_abs_corr will have no effect. Pass X or Xcorr to enable ",
      "the purity filter."
    )
  }
  if (check_symmetric) {
    if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr)) {
      Xcorr <- symmetrize_warned(Xcorr, "Xcorr")
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
  if (is.null(use_rfast)) {
    use_rfast <- requireNamespace("Rfast", quietly = TRUE)
  }
  
  # If no correlation info, return without purity or extension
  if (is.null(Xcorr) && is.null(X)) {
    names(cs) <- paste0("L", effect_indices)
    return(list(
      cs = cs,
      coverage = claimed_coverage,
      requested_coverage = coverage
    ))
  }
  
  # Extend CS by correlation if requested (works from either X or Xcorr; see
  # extend_cs_by_correlation). Recompute claimed coverage over the extended sets.
  if (!is.null(cs_extension_corr) && (!is.null(Xcorr) || !is.null(X))) {
    cs <- extend_cs_by_correlation(cs, X, Xcorr, cs_extension_corr, null_index)
    claimed_coverage <- vapply(seq_along(cs), function(i)
      sum(res$alpha[effect_indices[i], cs[[i]]]), numeric(1))
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
  
  if (!is.null(min_abs_corr) && !is.null(median_abs_corr)) {
    warning_message(
      "Both min_abs_corr and median_abs_corr are set; a credible set is kept ",
      "if it clears EITHER threshold (OR), which retains more sets than either ",
      "filter alone.",
      style = "hint")
  }
  # A CS is kept if it clears any active purity threshold (OR-linked).
  # min_abs_corr gates column 1 (min |corr|); median_abs_corr gates column 3
  # (median |corr|). Either may be NULL to disable that clause.
  keep <- rep(FALSE, nrow(purity))
  active <- FALSE
  if (!is.null(min_abs_corr)) {
    keep <- keep | (purity[, 1] >= if (squared) min_abs_corr^2 else min_abs_corr)
    active <- TRUE
  }
  if (!is.null(median_abs_corr)) {
    keep <- keep | (purity[, 3] >= if (squared) median_abs_corr^2 else median_abs_corr)
    active <- TRUE
  }
  if (!active) keep <- purity[, 1] > -1   # no threshold: keep all but the null CS (purity = -9)
  is_pure <- which(keep)
  
  if (length(is_pure) > 0) {
    cs <- cs[is_pure]
    purity <- purity[is_pure, , drop = FALSE]
    claimed_coverage <- claimed_coverage[is_pure]
    effect_indices <- effect_indices[is_pure]
    
    row_names <- paste0("L", effect_indices)
    names(cs) <- row_names
    rownames(purity) <- row_names
    
    # Re-order CS list and purity rows by the primary active statistic
    # (min |corr| if on, else median |corr|).
    order_col <- if (!is.null(min_abs_corr)) 1L
                 else if (!is.null(median_abs_corr)) 3L else 1L
    ordering <- order(purity[, order_col], decreasing = TRUE)
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

# Extend each credible set to include variables with |corr| > threshold
# with any current member. Works from a precomputed correlation matrix
# (Xcorr) or directly from the data matrix X; in the latter case the needed
# correlations are computed on demand as crossprod(standardize_X(X)), which
# equals cor(X) but never materializes the full p x p matrix (only an m x p
# block per credible set). The null index, if any, is never pulled in. Returns
# the (possibly extended) cs list in the same order; the caller recomputes
# claimed coverage.
#' @keywords internal
extend_cs_by_correlation <- function(cs, X, Xcorr, threshold, null_index = 0) {
  if (is.null(threshold) || length(cs) == 0 || (is.null(X) && is.null(Xcorr)))
    return(cs)

  Xs <- NULL
  if (is.null(Xcorr)) {
    if (inherits(X, "sparseMatrix")) {
      warning_message(
        "LD extension from a sparse X would require densifying the matrix; ",
        "skipping extension. Pass Xcorr (or a dense X) to enable it.",
        style = "hint"
      )
      return(cs)
    }
    # standardize_X(X) so that crossprod(Xs[, a], Xs) equals cor(X)[a, ].
    Xs <- standardize_X(X)
  }

  for (i in seq_along(cs)) {
    cs_idx <- cs[[i]]
    if (!is.null(Xcorr)) {
      ld_with_cs <- abs(Xcorr[cs_idx, , drop = FALSE]) > threshold
    } else {
      ld_with_cs <- abs(crossprod(Xs[, cs_idx, drop = FALSE], Xs)) > threshold
    }
    in_tight_ld <- which(colSums(ld_with_cs) > 0)
    if (null_index > 0)
      in_tight_ld <- setdiff(in_tight_ld, null_index)
    cs[[i]] <- sort(unique(c(cs_idx, in_tight_ld)))
  }
  cs
}

#' @rdname susie_get_methods
#'
#' @description \code{susie_get_cs_attainable} returns credible sets after
#'   a post-hoc filter based on \emph{attainable coverage} and per-effect
#'   entropy. Use this as a fallback when no LD reference is available
#'   for the standard purity filter in \code{susie_get_cs}.
#'
#' @param ethres Entropy threshold expressed as an effective number of
#'   variables: an effect is dropped when its attainable-projection
#'   entropy is at least \code{log(ethres)}. Defaults to
#'   \code{max(100, 0.1 * p)} where \code{p} is the number of variables.
#'
#' @details For each effect \eqn{l}, the \emph{attainable coverage} is
#'   computed by zeroing every entry of \code{res$alpha} that is not the
#'   column maximum, then summing across variables. Intuitively, it is
#'   the coverage effect \eqn{l} could attain if it did not have to
#'   share probability mass with other effects. Effects with attainable
#'   coverage at most \code{coverage}, or with attainable-projection
#'   entropy at least \code{log(ethres)}, are dropped before delegating
#'   to \code{susie_get_cs}. Any \code{X} or \code{Xcorr} passed through
#'   \code{...} is ignored, since the procedure is intended for the
#'   no-LD setting.
#'
#' @export
#'
susie_get_cs_attainable <- function(res, coverage = 0.95, ethres = NULL,
                                    ...) {
  p <- ncol(res$alpha)
  if (is.null(ethres))
    ethres <- max(100, 0.1 * p)

  cs_entropy <- function(y) {
    if (!any(y > 0)) return(Inf)
    y <- y / sum(y)
    -sum(y[y > 0] * log(y[y > 0]))
  }

  col_max <- apply(res$alpha, 2, max)
  alpha_attainable <- res$alpha *
    (res$alpha == matrix(col_max, nrow(res$alpha), p, byrow = TRUE))

  coverage_attained <- rowSums(alpha_attainable)
  entropy_vals <- apply(alpha_attainable, 1, cs_entropy)
  keep_effects <- which(coverage_attained > coverage &
                          entropy_vals < log(ethres))

  if (length(keep_effects) == 0)
    return(list(cs = NULL, coverage = NULL, requested_coverage = coverage))

  res_filtered <- res
  res_filtered$alpha <- res$alpha[keep_effects, , drop = FALSE]
  if (is.numeric(res$V) && length(res$V) == nrow(res$alpha))
    res_filtered$V <- res$V[keep_effects]

  dot_args <- list(...)
  dot_args$X <- NULL
  dot_args$Xcorr <- NULL

  out <- do.call(susie_get_cs,
                 c(list(res_filtered, coverage = coverage), dot_args))

  # Remap effect indices: susie_get_cs returned indices into res_filtered;
  # translate them back to the original L positions.
  if (!is.null(out$cs)) {
    filtered_idx <- as.integer(sub("^L", "", names(out$cs)))
    original_idx <- keep_effects[filtered_idx]
    names(out$cs) <- paste0("L", original_idx)
    if (!is.null(out$cs_index)) # nocov
      out$cs_index <- original_idx # nocov
    if (!is.null(out$purity)) # nocov
      rownames(out$purity) <- paste0("L", original_idx) # nocov
  }

  out
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
    Xcorr <- symmetrize_warned(Xcorr, "Xcorr")
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

    # Extract slot weights (c_hat) if available for Gamma-Poisson weighting.
    # PIP_j = 1 - prod_l(1 - c_hat[l] * alpha[l,j])
    slot_wt <- res$slot_weights

    # now extract relevant rows from alpha matrix
    if (length(include_idx) > 0) {
      res_alpha <- res$alpha[include_idx, , drop = FALSE]
      if (!is.null(slot_wt)) {
        slot_wt <- slot_wt[include_idx]
      }
    } else {
      res_alpha <- matrix(0, 1, ncol(res$alpha))
      slot_wt <- NULL
    }
    res <- res_alpha
  }

  # c_hat-weighted PIPs when slot_weights are available
  if (exists("slot_wt", inherits = FALSE) && !is.null(slot_wt)) {
    weighted_alpha <- sweep(res, 1, slot_wt, `*`)
    return(as.vector(1 - apply(1 - weighted_alpha, 2, prod)))
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
  out = list(alpha = alpha, mu = mu, mu2 = mu*mu)
  class(out) = c("susie","list")
  return(out)
}
