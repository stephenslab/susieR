# =============================================================================
# FUNDAMENTAL BUILDING BLOCKS
#
# Basic mathematical operations and utilities that serve as dependencies
# for other functions. These include matrix operations, statistical computations,
# and general-purpose helper functions.
#
# Functions: warning_message, safe_cor, safe_cov2cor, is_symmetric_matrix,
# apply_nonzeros, compute_colSds, compute_colstats
# =============================================================================

# Report R process memory usage (GB). Uses gc() which is cheap.
#' @keywords internal
mem_used_gb <- function() {
  gc_info <- gc(verbose = FALSE, reset = FALSE)
  sum(gc_info[, "(Mb)"]) / 1024
}

# Format prior variance vector: show non-zero values, summarize zeros.
# E.g., "[1.23e-01, 5.34e-02, 0 x 3]"
#' @keywords internal
format_V_summary <- function(V) {
  n_na   <- sum(is.na(V))
  V_nona <- V[!is.na(V)]
  n_zero <- sum(V_nona == 0)
  nz <- V_nona[V_nona != 0]
  parts <- sprintf("%.2e", nz)
  if (n_zero > 0) parts <- c(parts, sprintf("0 x %d", n_zero))
  if (n_na > 0)   parts <- c(parts, sprintf("NA x %d", n_na))
  paste0("[", paste(parts, collapse = ", "), "]")
}

# Format slot activity (c_hat) summary for verbose output.
# Shows per-slot c_hat values and lbf when active, empty string when not.
# E.g., ", c_hat=[0.99,0.87,0.12,0.45], lbf=[12.3,8.1,0.0,-0.2], C_hat=2.4"
#' @keywords internal
format_chat_summary <- function(model) {
  if (is.null(model$slot_weights)) return("")
  sw <- model$slot_weights
  lbf <- model$lbf
  lbf[is.na(lbf)] <- 0
  chat_vals <- paste(sprintf("%.2f", sw), collapse = ",")
  lbf_vals <- paste(sprintf("%.1f", lbf), collapse = ",")
  n_active <- sum(sw > 0.5)
  sprintf(", c_hat=[%s], lbf=[%s], C_hat=%.1f(%d>0.5)",
          chat_vals, lbf_vals, sum(sw), n_active)
}

# Utility function to display warning messages as they occur
#' @importFrom crayon combine_styles
#' @keywords internal
warning_message <- function(..., style = c("warning", "hint")) {
  style <- match.arg(style)
  if (style == "warning" && getOption("warn") >= 0) {
    alert <- combine_styles("bold", "underline", "red")
    message(alert("WARNING:"), " ", ...)
  } else {
    alert <- combine_styles("bold", "underline", "magenta")
    message(alert("HINT:"), " ", ...)
  }
}

# Auto-switch convergence_method to "pip". Emits a hint explaining why.
# Used at every site where the chosen method makes ELBO ill-defined or
# unreliable. Returns "pip" so callers can do
# `convergence_method <- force_pip(...)`. Variadic, like
# `warning_message()`.
#' @keywords internal
force_pip <- function(...) {
  warning_message("Switching to PIP-based convergence: ", ..., ".",
                  style = "hint")
  "pip"
}

# Standard "please supply n" hint used by every RSS entry that allows
# n = NULL.
#' @keywords internal
hint_sample_size_recommended <- function() {
  warning_message(
    "Providing the sample size (n), or even a rough estimate of n, ",
    "is highly recommended. Without n, the implicit assumption is ",
    "n is large (Inf) and the effect sizes are small (close to zero).",
    style = "hint")
}

# Symmetrize a matrix and warn. Returns (M + t(M)) / 2.
#' @keywords internal
symmetrize_warned <- function(M, name) {
  warning_message(name, " is not symmetric; using (", name, " + t(",
                  name, "))/2.")
  (M + t(M)) / 2
}

# Replace NAs in a numeric vector with 0 and warn.
#' @keywords internal
replace_na_zero_warned <- function(v, name) {
  warning_message("NA values in ", name, " are replaced with 0.")
  v[is.na(v)] <- 0
  v
}

# Warn that R is not positive semidefinite (callers clamp the negative
# eigenvalues themselves, since they hold the local eigenvalue handle).
#' @keywords internal
warn_R_not_psd <- function() {
  warning_message(
    "The matrix R is not positive semidefinite. Negative ",
    "eigenvalues are set to zero.")
}

#' Converts covariance matrix to correlation matrix
#' Constant variables (zero variance) get correlation 0 with others, 1 with self
#'
#' @param V Covariance matrix
#' @return Correlation matrix
#' @keywords internal
safe_cov2cor <- function(V) {
  d <- sqrt(diag(V))
  d_inv <- 1 / d
  d_inv[d == 0] <- 0
  R <- V * outer(d_inv, d_inv)
  diag(R) <- 1
  R
}

#' Computes correlation matrix from data matrix
#' Handles constant columns without warnings - returns 0 correlation for constant cols
#' Uses Rfast::cora when available (much faster for large matrices), falls back
#' to crossprod-based computation otherwise.
#'
#' @param X Data matrix (n x p)
#' @return Correlation matrix (p x p)
#' @keywords internal
safe_cor <- function(X) {
  n <- nrow(X)
  cm <- colMeans(X)
  css <- colSums(X^2) - n * cm^2   # column sum of squares (centered)
  has_const <- any(css == 0)

  # Fast path: use Rfast::cora when available and no constant columns
  if (!has_const && requireNamespace("Rfast", quietly = TRUE)) {
    return(Rfast::cora(X))
  }

  # Fallback: manual crossprod, handling constant columns
  X_centered <- X - rep(cm, each = n)
  sds <- sqrt(css / n)
  sds_inv <- 1 / sds
  sds_inv[sds == 0] <- 0
  X_scaled <- X_centered * rep(sds_inv, each = n)
  R <- crossprod(X_scaled) / n
  diag(R) <- 1
  R
}

# Standardize X so that X'X equals the correlation matrix R.
# Centers columns, divides by column sd, then scales by 1/sqrt(n)
# so that crossprod(X_out) = cor(X_in). Constant columns are zeroed.
# This is the in-place analog of safe_cor: safe_cor(X) == crossprod(standardize_X(X)).
#' @keywords internal
standardize_X <- function(X) {
  n <- nrow(X)
  cm <- colMeans(X)
  X <- X - rep(cm, each = n)
  css <- colSums(X^2)
  sds <- sqrt(css / n)
  sds[sds < .Machine$double.eps] <- 1  # constant columns: avoid 0/0
  X <- X * rep(1 / (sds * sqrt(n)), each = n)
  X
}

#' Check for symmetric matrix
#'
#' @param x A matrix to check
#' @return Logical indicating if x is symmetric
#' @export
#' @keywords internal
is_symmetric_matrix <- function(x) {
  if (is.matrix(x) && is.numeric(x) && !isS4(x) &&
      requireNamespace("Rfast", quietly = TRUE)) {
    return(Rfast::is.symmetric(x))
  } else {
    return(Matrix::isSymmetric(x, check.attributes = FALSE))
  }
}

# Apply operation f to all nonzeros of a sparse matrix.
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix summary
#' @keywords internal
apply_nonzeros <- function(X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i, j = d$j, x = f(d$x), dims = dim(X)))
}

# Computes column standard deviations for any type of matrix
# This should give the same result as matrixStats::colSds(X),
# but allows for sparse matrices as well as dense ones.
#' @importFrom matrixStats colSds
#' @importFrom Matrix summary
#' @keywords internal
compute_colSds <- function(X) {
  if (is.matrix(X)) {
    return(colSds(X))
  } else {
    n <- nrow(X)
    Y <- apply_nonzeros(X, function(u) u^2)
    d <- colMeans(Y) - colMeans(X)^2
    return(sqrt(d * n / (n - 1)))
  }
}

# Compute the column means of X, the column standard deviations of X,
# and rowSums(Y^2), where Y is the centered and/or scaled version of
# X.
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @keywords internal
compute_colstats <- function(X, center = TRUE, scale = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  if (!is.null(attr(X, "matrix.type"))) {
    # X is a trend filtering matrix.
    cm <- compute_tf_cm(attr(X, "order"), p)
    csd <- compute_tf_csd(attr(X, "order"), p)
    d <- compute_tf_d(attr(X, "order"), p, cm, csd, scale, center)
    if (!center) {
      cm <- rep(0, p)
    }
    if (!scale) {
      csd <- rep(1, p)
    }
  } else {
    # X is an ordinary dense or sparse matrix. Set sd = 1 when the
    # column has variance 0.
    if (center) {
      cm <- colMeans(X, na.rm = TRUE)
    } else {
      cm <- rep(0, p)
    }
    if (scale) {
      csd <- compute_colSds(X)
      csd[csd == 0] <- 1
    } else {
      csd <- rep(1, p)
    }

    # These two lines of code should give the same result as
    #
    #   Y = (t(X) - cm)/csd
    #   d = rowSums(Y^2)
    #
    # for all four combinations of "center" and "scale", but do so
    # without having to modify X, or create copies of X in memory. In
    # particular the first line should be equivalent to colSums(X^2).
    d <- n * colMeans(X)^2 + (n - 1) * compute_colSds(X)^2
    d <- (d - n * cm^2) / csd^2
  }

  return(list(cm = cm, csd = csd, d = d))
}

# Compute standard error for regression coef.
# S = (X'X)^-1 \Sigma
calc_stderr = function (X, residuals)
  sqrt(diag(sum(residuals^2)/(nrow(X) - 2) * chol2inv(chol(crossprod(X)))))

#' Check alpha/PIP fixed-point or short-cycle convergence
#'
#' Uses one tolerance for both marginal PIPs and alpha. Lag 1 is ordinary
#' convergence; larger lags detect a periodic orbit and average alpha over it.
#'
#' @keywords internal
check_alpha_pip_cycle_convergence <- function(data, params, model) {
  tol <- params$tol
  cycle_window <- if (!is.null(params$pip_stall_window))
    params$pip_stall_window
  else
    5
  cycle_window <- max(1L, as.integer(cycle_window))
  prior_tol <- if (!is.null(params$prior_tol)) params$prior_tol else 1e-9

  alpha_pip <- function(alpha) {
    tmp <- model
    tmp$alpha <- alpha
    if (inherits(tmp, "susie"))
      susie_get_pip(tmp, prior_tol = prior_tol)
    else
      susie_get_pip(alpha, prior_tol = prior_tol)
  }

  current_alpha <- model$alpha
  current_pip <- alpha_pip(current_alpha)

  alpha_history <- model$runtime$alpha_history
  pip_history <- model$runtime$pip_history
  if (is.null(alpha_history) || is.null(pip_history)) {
    alpha_history <- list(model$runtime$prev_alpha)
    pip_history <- list(alpha_pip(model$runtime$prev_alpha))
  }

  state_diff <- function(alpha_old, pip_old)
    max(max(abs(current_alpha - alpha_old)), max(abs(current_pip - pip_old)))

  max_lag <- min(cycle_window, length(alpha_history))
  lag_diff <- state_diff(alpha_history[[length(alpha_history)]],
                         pip_history[[length(pip_history)]])
  reason <- NULL

  if (lag_diff < tol) {
    reason <- "alpha_pip_fixed_point"
  } else if (max_lag >= 2) {
    for (lag in 2:max_lag) {
      idx <- length(alpha_history) - lag + 1
      lag_diff <- state_diff(alpha_history[[idx]], pip_history[[idx]])
      if (lag_diff < tol) {
        reason <- paste0("alpha_pip_cycle_", lag)
        cycle_alpha <- c(tail(alpha_history, lag - 1), list(current_alpha))
        model$alpha <- Reduce(`+`, cycle_alpha) / lag
        current_alpha <- model$alpha
        current_pip <- alpha_pip(current_alpha)
        break
      }
    }
  }

  model$converged <- !is.null(reason)
  model$convergence_reason <- reason
  model$runtime$pip_diff <- lag_diff
  model$runtime$alpha_history <- c(alpha_history, list(current_alpha))
  model$runtime$pip_history <- c(pip_history, list(current_pip))
  if (length(model$runtime$alpha_history) > cycle_window) {
    keep <- seq.int(length(model$runtime$alpha_history) - cycle_window + 1,
                    length(model$runtime$alpha_history))
    model$runtime$alpha_history <- model$runtime$alpha_history[keep]
    model$runtime$pip_history <- model$runtime$pip_history[keep]
  }
  return(model)
}

# =============================================================================
# DATA PROCESSING & VALIDATION
#
# Functions for input validation, data conversion between formats, and
# preprocessing operations. These ensure data integrity and compatibility
# across different SuSiE data types.
#
# Functions: check_semi_pd, check_projection, validate_init,
# convert_individual_to_ss, extract_prior_weights, reconstruct_full_weights,
# validate_and_override_params
# =============================================================================

# Check whether A is positive semidefinite
#' @keywords internal
check_semi_pd <- function(A, tol) {
  attr(A, "eigen") <- eigen(A, symmetric = TRUE)
  v <- attr(A, "eigen")$values
  v[abs(v) < tol] <- 0
  return(list(
    matrix = A,
    status = !any(v < 0),
    eigenvalues = v
  ))
}

# Fraction of b projected onto eigenvectors with eigenvalues below a threshold.
#' @keywords internal
low_eigen_projection_fraction <- function(eigen_A, b, eigen_tol,
                                          norm_floor = 0) {
  d <- eigen_A$values
  V <- eigen_A$vectors
  low <- which(d <= eigen_tol)
  b_norm <- sum(b^2)
  if (length(low) == 0L || b_norm <= norm_floor) {
    return(list(fraction = 0, evaluable = FALSE,
                low_eigen_count = length(low), eigen_tol = eigen_tol))
  }
  coeff <- as.numeric(crossprod(V[, low, drop = FALSE], b))
  list(fraction = sum(coeff^2) / b_norm, evaluable = TRUE,
       low_eigen_count = length(low), eigen_tol = eigen_tol)
}

# Check whether b is in space spanned by the non-zero eigenvectors of A.
# A may be a matrix, optionally with attr(A, "eigen"), or an eigen() list.
#' @keywords internal
check_projection <- function(A, b, tol = .Machine$double.eps,
                             projection_tol = tol, norm_floor = 0) {
  eig <- if (is.list(A) && !is.null(A$values) && !is.null(A$vectors)) {
    A
  } else {
    if (is.null(attr(A, "eigen")))
      attr(A, "eigen") <- eigen(A, symmetric = TRUE)
    attr(A, "eigen")
  }
  proj <- low_eigen_projection_fraction(eig, b, tol, norm_floor)
  status <- !proj$evaluable || proj$fraction <= projection_tol
  msg <- if (status) NA else
    paste0("relative squared projection ", signif(proj$fraction, 3),
           " exceeds ", signif(projection_tol, 3))
  c(list(status = status, msg = msg), proj)
}

# Validate Model Initialization Object
#' @keywords internal
validate_init <- function(data, params) {
  if (!inherits(params$model_init, "susie")) {
    stop("model_init must be a 'susie' object")
  }

  # Assign values from initialized model
  L       <- params$L
  alpha   <- params$model_init$alpha
  mu      <- params$model_init$mu
  mu2     <- params$model_init$mu2
  V       <- params$model_init$V
  sigma2  <- params$model_init$sigma2
  pi_w    <- params$model_init$pi
  null_id <- params$model_init$null_index

  # Verify no NA/Inf values in alpha
  if (any(!is.finite(alpha))) {
    stop("model_init$alpha contains NA/Inf values")
  }

  # Verify no NA/Inf values in mu
  if (any(!is.finite(mu))) {
    stop("model_init$mu contains NA/Inf values")
  }

  # Verify no NA/Inf values in mu2
  if (any(!is.finite(mu2))) {
    stop("model_init$mu2 contains NA/Inf values")
  }

  # Only check V if it exists
  if (!is.null(V)) {
    # Verify no NA/Inf values in V
    if (any(!is.finite(V))) {
      stop("model_init$V contains NA/Inf values")
    }
  }

  # Only check sigma2 if it exists
  if (!is.null(sigma2)) {
    # Verify no NA/Inf values in sigma2
    if (any(!is.finite(sigma2))) {
      stop("model_init$sigma2 contains NA/Inf")
    }
  }

  # Only check pi_w if it exists
  if (!is.null(pi_w)) {
    # Verify no NA/Inf values in prior weights
    if (any(!is.finite(pi_w))) {
      stop("model_init$pi contains NA/Inf")
    }
  }

  # Verify alpha is matrix
  if (!is.matrix(alpha)) {
    stop("model_init$alpha must be a matrix")
  }

  # Verify alpha values are between [0,1]
  if (max(alpha) > 1 || min(alpha) < 0) {
    stop(
      "model_init$alpha has invalid values outside range [0,1]; please ",
      "check your input"
    )
  }

  # Verify mu & mu2 dimensions match alpha
  if (!all(dim(mu) == dim(alpha))) {
    stop("model_init$mu and model_init$alpha dimensions do not match")
  }
  if (!all(dim(mu2) == dim(alpha))) {
    stop("model_init$mu2 and model_init$alpha dimensions do not match")
  }

  # Only validate V dimensions and values if V exists
  if (!is.null(V)) {
    # Verify V & alpha dimensions agree
    if (length(V) != nrow(alpha)) {
      stop(
        "length(model_init$V) (", length(V), ") does not equal nrow(model_init$alpha) (",
        nrow(alpha), ")"
      )
    }

    # Verify V is numeric and non-negative
    if (!is.numeric(V)) {
      stop("model_init$V must be numeric")
    }
    if (any(V < 0)) {
      stop("model_init$V has at least one negative value")
    }
  }

  # Verify sigma2 is numeric and non-negative if it exists
  if (!is.null(sigma2)) {
    if (!is.numeric(sigma2)) {
      stop("model_init$sigma2 must be numeric")
    }
    if (sigma2 < 0) {
      stop("model_init$sigma2 is negative")
    }
  }

  # Verify prior weight properties if they exist
  if (!is.null(pi_w)) {
    if (length(pi_w) != ncol(alpha)) {
      stop(
        "model_init$pi should have the same length as the number of columns",
        " in model_init$alpha"
      )
    }
  }

  invisible(params$model_init)
}

# Convert individual data to ss with unmappable effects components.
#' @keywords internal
convert_individual_to_ss <- function(data, params) {
  # Compute sufficient statistics
  XtX <- compute_XtX(data$X)
  Xty <- compute_Xty(data$X, data$y)
  yty <- sum(data$y^2)

  # Get column means and scaling from attributes
  X_colmeans <- attr(data$X, "scaled:center")

  # Create sufficient statistics data object
  ss_data <- structure(
    list(
      XtX = XtX,
      X = NULL,
      Xty = Xty,
      yty = yty,
      n = data$n,
      p = data$p,
      X_colmeans = X_colmeans,
      y_mean = data$mean_y
    ),
    class = "ss"
  )

  # Set attributes on XtX from individual X
  attr(ss_data$XtX, "d") <- attr(data$X, "d")
  attr(ss_data$XtX, "scaled:scale") <- attr(data$X, "scaled:scale")

  # Add eigen decomposition for unmappable effects methods
  ss_data <- add_eigen_decomposition(ss_data, params, data)

  return(ss_data)
}

# Extract non-null prior weights from a model
#' @keywords internal
extract_prior_weights <- function(model, null_weight = NULL) {
  # Use model's null_weight if not provided (backwards compatibility)
  if (is.null(null_weight)) {
    null_weight <- model$null_weight
  }
  
  if (!is.null(null_weight) && null_weight != 0 && !is.null(model$null_index) && model$null_index != 0) {
    # Extract non-null prior weights and rescale
    pw_s <- model$pi[-model$null_index] / (1 - null_weight)
  } else {
    pw_s <- model$pi
  }
  return(pw_s)
}

# Reconstruct full prior weights with null weight handling
#' @keywords internal
reconstruct_full_weights <- function(non_null_weights, null_weight) {
  if (!is.null(null_weight) && null_weight != 0) {
    # Reconstruct full prior weights including null component
    full_weights <- c(non_null_weights * (1 - null_weight), null_weight)
  } else {
    full_weights <- non_null_weights
  }
  # Normalize to sum to 1
  return(full_weights / sum(full_weights))
}


# Validate and Override Parameters
#' @keywords internal
validate_and_override_params <- function(params) {

  # Validate prior tolerance threshold
  if (!is.numeric(params$prior_tol) || length(params$prior_tol) != 1) {
    stop("prior_tol must be a numeric scalar.")
  }
  if (params$prior_tol < 0) {
    stop("prior_tol must be non-negative.")
  }

  tol <- params$tol
  use_NIG_method <- identical(params$estimate_residual_method, "NIG")
  if (is.null(tol) && use_NIG_method)
    warning_message("tol was not specified; setting tol = 1e-6 for ",
                    "the NIG prior to improve convergence.",
                    style = "hint")
  params$tol <- if (is.null(tol)) {
    if (use_NIG_method) 1e-6 else 1e-4
  } else {
    tol
  }
  if (!is.numeric(params$tol) || length(params$tol) != 1 ||
      is.na(params$tol) || !is.finite(params$tol) ||
      params$tol < 0) {
    stop("tol must be a non-negative numeric scalar.")
  }

  # Validate greedy-L parameters.
  if (!is.null(params$L_greedy)) {
    if (!is.numeric(params$L_greedy) || length(params$L_greedy) != 1 ||
        is.na(params$L_greedy) || !is.finite(params$L_greedy) ||
        params$L_greedy < 1 ||
        params$L_greedy != as.integer(params$L_greedy)) {
      stop("L_greedy must be NULL or a positive integer.")
    }
    params$L_greedy <- as.integer(params$L_greedy)
    if (params$L_greedy > params$L) {
      warning_message("L_greedy is greater than L; using L instead.",
                      style = "hint")
      params$L_greedy <- as.integer(params$L)
    }
  }
  if (!is.numeric(params$greedy_lbf_cutoff) ||
      length(params$greedy_lbf_cutoff) != 1 ||
      is.na(params$greedy_lbf_cutoff) ||
      !is.finite(params$greedy_lbf_cutoff)) {
    stop("greedy_lbf_cutoff must be a numeric scalar.")
  }

  # Validate residual_variance_upperbound
  if (!is.numeric(params$residual_variance_upperbound) || length(params$residual_variance_upperbound) != 1) {
    stop("residual_variance_upperbound must be a numeric scalar.")
  }
  if (params$residual_variance_upperbound <= 0) {
    stop("residual_variance_upperbound must be positive.")
  }

  # Validate scaled prior variance
  if (!is.numeric(params$scaled_prior_variance) || any(params$scaled_prior_variance < 0)) {
    stop("Scaled prior variance should be positive number.")
  }
  spv_len <- length(params$scaled_prior_variance)
  if (spv_len != 1 && spv_len != params$L) {
    stop("scaled_prior_variance must be a scalar or a vector of length L.")
  }

  # Validate unmappable_effects
  # "ash_filter_archived" is a hidden option for internal diagnostics/archiving
  # of the original SuSiE-ASH filter-based masking heuristic.
  valid_unmappable <- c("none", "inf", "ash", "ash_filter_archived")
  if (!params$unmappable_effects %in% valid_unmappable) {
    stop("unmappable_effects must be one of 'none', 'inf', or 'ash'.")
  }

  # Auto-create Beta-Binomial slot prior for ash mode (identifiability)
  # and report defaults as a single hint. Users who want an
  # uninformative prior can pass slot_prior_betabinom(1, 1).
  if (params$unmappable_effects == "ash" && is.null(params$slot_prior)) {
    params$slot_prior <- slot_prior_betabinom()
  }
  if (!is.null(params$slot_prior) && inherits(params$slot_prior, "slot_prior_betabinom") &&
      isTRUE(params$slot_prior$ab_was_default)) {
    sp <- params$slot_prior
    n_active <- round(sp$a_beta / (sp$a_beta + sp$b_beta) * params$L)
    warning_message(
      "For SuSiE-ash, slot_prior was not specified; using default ",
      "Beta(a=", sp$a_beta, ", b=", sp$b_beta, "), ",
      "roughly expecting ~", n_active, " of ", params$L,
      " slots to be active. Set slot_prior = slot_prior_betabinom(a, b) ",
      "to override; for an uninformative prior use slot_prior_betabinom(1, 1).",
      style = "hint")
  }
  # Report nu default after the C message so the user sees C first.
  if (!is.null(params$slot_prior) && !inherits(params$slot_prior, "slot_prior_betabinom") &&
      isTRUE(params$slot_prior$nu_was_default)) {
    sp <- params$slot_prior
    sd_mu <- sp$C / sqrt(sp$nu)
    warning_message(
      "Overdispersion parameter nu not specified, using default nu = ",
      sp$nu, ". With C = ", sp$C,
      " this implies sd(mu) = ", round(sd_mu, 2),
      ", so the number of active effects ranges roughly from ",
      round(max(0, sp$C - 2 * sd_mu), 1), " to ",
      round(sp$C + 2 * sd_mu, 1),
      " around the prior mean. Set nu explicitly to change this behavior.",
      style = "hint")
  }

  # Override convergence method for unmappable effects or slot_prior.
  # The ELBO is not well-defined when slot_weights != 1 (c_hat active)
  # or when unmappable effects modify the residual structure.
  needs_pip <- params$unmappable_effects != "none" || !is.null(params$slot_prior)
  if (needs_pip && params$convergence_method != "pip") {
    reason <- if (params$unmappable_effects != "none")
      "unmappable effects models (inf/ash) do not have a well-defined ELBO"
    else
      paste("slot activity model modifies fitted values by slot weights,",
            "making the standard ELBO invalid")
    params$convergence_method <- force_pip(reason)
  }

  # Check for incompatible parameter combinations
  if (!is.null(params$refine) && params$refine && params$unmappable_effects != "none") {
    stop("Refinement is not supported with unmappable effects (inf/ash) as it relies on ELBO, ",
         "which is not well-defined for these models. Please set refine = FALSE.")
  }

  # Override prior estimation method when estimation is disabled,
  # unless using a fixed mixture prior (which does not estimate V but
  # still needs the mixture BF computation path)
  if (!params$estimate_prior_variance &&
      params$estimate_prior_method != "fixed_mixture") {
    params$estimate_prior_method <- "none"
  }

  # Handle NIG parameters for small sample correction
  if (use_NIG_method) {
    params$use_NIG <- TRUE

    # Require a valid sample size n. The default alpha0/beta0 scale as
    # 1/sqrt(n), so n must be a positive finite scalar. susie() infers
    # this from nrow(X); susie_ss()'s constructor enforces it; susie_rss()
    # allows n = NULL by default, so users who select NIG must also
    # supply `n` explicitly.
    if (is.null(params$n) || !is.numeric(params$n) || length(params$n) != 1 ||
        !is.finite(params$n) || params$n < 1) {
      stop("estimate_residual_method = \"NIG\" requires a valid sample ",
           "size `n` (got n = ", paste(params$n, collapse = ""), "). ",
           "susie() infers n from nrow(X); for susie_ss() and susie_rss(), ",
           "pass `n` explicitly.")
    }

    if (is.null(params$alpha0))
      params$alpha0 <- 1 / sqrt(params$n)
    if (is.null(params$beta0))
      params$beta0 <- 1 / sqrt(params$n)

    # Validate NIG prior parameters: both must be strictly positive for a proper
    # Inverse-Gamma prior. Otherwise compute_null_loglik_NIG() evaluates
    # lgamma(alpha0 / 2) at <= 0 and the marginal log-likelihood (and ELBO)
    # become Inf or NaN.
    if (!is.numeric(params$alpha0) || length(params$alpha0) != 1 ||
        !is.finite(params$alpha0) || params$alpha0 <= 0 ||
        !is.numeric(params$beta0)  || length(params$beta0)  != 1 ||
        !is.finite(params$beta0)  || params$beta0  <= 0) {
      stop("estimate_residual_method = \"NIG\" requires ",
           "alpha0 > 0 and beta0 > 0 (proper Inverse-Gamma prior). ",
           "Got alpha0 = ", params$alpha0, ", beta0 = ", params$beta0, ". ",
           "The default is alpha0 = beta0 = 1/sqrt(n).")
    }

    # The NIG prior inherently estimates residual variance (integrates out sigma^2).
    # If estimate_residual_variance is FALSE, override it -- the user chose a method
    # that estimates sigma^2 by design. To suppress this warning, explicitly set
    # estimate_residual_variance = TRUE in the function call.
    if (!isTRUE(params$estimate_residual_variance)) {
      warning_message("NIG prior integrates out residual variance, ",
                      "implying estimate_residual_variance = TRUE. ",
                      "Setting estimate_residual_variance = TRUE. ",
                      "To suppress this hint, explicitly set ",
                      "estimate_residual_variance = TRUE in the function call.",
                      style = "hint")
      params$estimate_residual_variance <- TRUE
    }

    # Override convergence method only when L > 1
    if (params$L > 1 && params$convergence_method != "pip") {
      params$convergence_method <- force_pip("NIG method with L > 1")
    }

    # Override prior variance estimation method (only when estimation is enabled)
    if (params$estimate_prior_variance && params$estimate_prior_method != "EM") {
      warning_message("NIG method works better with EM. ",
                      "Setting estimate_prior_method='EM'.",
                      style = "hint")
      params$estimate_prior_method <- "EM"
    }
  } else {
    params$use_NIG <- FALSE
    params$alpha0 <- NULL
    params$beta0 <- NULL
  }

  return(params)
}

# =============================================================================
# MODEL INITIALIZATION
#
# Functions that set up initial model states, create model matrices,
# and handle model configuration. These prepare the SuSiE model object
# for iterative fitting.
#
# Functions: initialize_null_index, assign_names,
# adjust_L, prune_single_effects, add_null_effect
# =============================================================================

# Initialize Null Index
#' @keywords internal
initialize_null_index <- function(data, model) {
  if (is.null(model$null_weight) || model$null_weight == 0) {
    null_idx <- 0
  } else {
    null_idx <- data$p
  }
  return(null_idx)
}

# Helper function to assign variable names to model components
#' @keywords internal
assign_names <- function(data, model, variable_names) {
  if (!is.null(variable_names)) {
    if (!is.null(model$null_weight) && model$null_weight != 0 && !is.null(model$null_index) && model$null_index != 0) {
      variable_names[length(variable_names)] <- "null"
      names(model$pip) <- variable_names[-data$p]
    } else {
      names(model$pip) <- variable_names
    }
    colnames(model$alpha)        <- variable_names
    colnames(model$mu)           <- variable_names
    colnames(model$mu2)          <- variable_names
    colnames(model$lbf_variable) <- variable_names
  }
  return(model)
}

# Expand scaled_prior_variance into a length-L vector of prior variances.
# Accepts either a scalar (recycled to length L) or a length-L vector
# (used as-is, one prior variance per single-effect slot).
#' @keywords internal
expand_scaled_prior_variance <- function(scaled_prior_variance, var_y, L) {
  if (length(scaled_prior_variance) == 1) {
    rep(scaled_prior_variance * var_y, L)
  } else {
    scaled_prior_variance * var_y
  }
}

# Adjust the number of effects
#' @keywords internal
adjust_L <- function(params, model_init_pruned, var_y) {
  num_effects <- nrow(model_init_pruned$alpha)
  L <- params$L

  if (num_effects > L) {
    warning_message(paste0(
      "Requested L = ", L,
      " is smaller than the ", num_effects,
      " effects in model_init after pruning; ",
      "using L = ", num_effects, " instead."
    ), style = "hint")
    L <- num_effects
  }

  V <- expand_scaled_prior_variance(params$scaled_prior_variance, var_y, L)
  model_init <- extract_model_init_fields(
    model_init_pruned, L = L, V = V,
    estimate_residual_variance = params$estimate_residual_variance)

  return(list(model_init = model_init, L = L))
}

# Extract the fitted state needed to initialize a new SuSiE fit.
#' @keywords internal
extract_model_init_fields <- function(model_init, L = 0, V = NULL,
                                      estimate_residual_variance = TRUE) {
  model_init <- prune_single_effects(model_init, L = L, V = V)
  keep <- c("alpha", "mu", "mu2", "V")
  if (isTRUE(estimate_residual_variance))
    keep <- c(keep, "sigma2")
  model_init[intersect(keep, names(model_init))]
}

# Prune single effects to given number L in susie model object.
#' @keywords internal
prune_single_effects <- function(model_init, L = 0, V = NULL) {
  num_effects <- nrow(model_init$alpha)
  if (L == 0) {
    # Filtering will be based on non-zero elements in model_init$V.
    if (!is.null(model_init$V)) {
      L <- length(which(model_init$V > 0))
    } else {
      L <- num_effects
    }
  }
  if (L == num_effects) {
    model_init$sets <- NULL
    return(model_init)
  }
  if (!is.null(model_init$sets$cs_index)) {
    effects_rank <- c(model_init$sets$cs_index, setdiff(1:num_effects, model_init$sets$cs_index))
  } else {
    effects_rank <- 1:num_effects
  }

  if (L > num_effects) {
    message(paste(
      "Specified number of effects L =", L,
      "is greater the number of effects", num_effects,
      "in input SuSiE model. The SuSiE model will be expanded",
      "to have", L, "effects."
    ))

    model_init$alpha <- rbind(
      model_init$alpha[effects_rank, ],
      matrix(1 / ncol(model_init$alpha), L - num_effects, ncol(model_init$alpha))
    )
    for (n in c("mu", "mu2", "lbf_variable")) {
      if (!is.null(model_init[[n]])) {
        model_init[[n]] <- rbind(
          model_init[[n]][effects_rank, ],
          matrix(0, L - num_effects, ncol(model_init[[n]]))
        )
      }
    }
    for (n in c("KL", "lbf")) {
      if (!is.null(model_init[[n]])) {
        model_init[[n]] <- c(model_init[[n]][effects_rank], rep(NA, L - num_effects))
      }
    }
    if (!is.null(V)) {
      if (length(V) > 1) {
        if (!is.null(model_init$V)) {
          V[1:num_effects] <- model_init$V[effects_rank]
        }
      } else {
        V <- rep(V, L)
      }
    }
    model_init$V <- V
  }
  model_init$sets <- NULL
  return(model_init)
}

# Add a null effect to the model object
#' @keywords internal
add_null_effect <- function(model_init, V) {
  p                       <- ncol(model_init$alpha)
  model_init$alpha        <- rbind(model_init$alpha, 1 / p)
  model_init$mu           <- rbind(model_init$mu, rep(0, p))
  model_init$mu2          <- rbind(model_init$mu2, rep(0, p))
  model_init$lbf_variable <- rbind(model_init$lbf_variable, rep(0, p))
  model_init$V            <- c(model_init$V, V)
  return(model_init)
}

# =============================================================================
# MATRIX-VECTOR PRODUCT HELPERS
#
# Unified helpers for predictor-matrix-times-vector operations across
# SS (XtX) and RSS-lambda (R) data types. These dispatch on what's available
# on the data object: data$X (low-rank factor), data$XtX, or data$R.
# When data$X is stored (Bxp, B < p), the two-step product X'(Xv) avoids
# forming the pxp matrix, reducing cost from O(p^2) to O(Bp).
#
# Functions: compute_Rv, compute_BR
# =============================================================================

# Compute R*v product: X'(Xv), XtX*v, or R*v.
# For multi-panel rss_lambda, pass Rv_matrix = model$X_meta to use the
# current R(omega) factor instead of data$X.
#' @keywords internal
compute_Rv <- function(data, v, Rv_matrix = NULL) {
  if (!is.null(Rv_matrix)) {
    return(as.vector(crossprod(Rv_matrix, Rv_matrix %*% v)))
  }
  if (!is.null(data$X)) {
    return(as.vector(crossprod(data$X, data$X %*% v)))
  } else if (!is.null(data$XtX)) {
    return(as.vector(data$XtX %*% v))
  } else if (!is.null(data$R)) {
    return(as.vector(data$R %*% v))
  }
  stop("No predictor matrix available on data object.")
}

# Compute B_mat %*% predictor-matrix: (Lxp) times (pxp) -> (Lxp)
# Used in get_ER2.ss for the quadratic form B %*% XtX
#' @keywords internal
compute_BR <- function(data, B_mat) {
  if (!is.null(data$X)) {
    return((B_mat %*% t(data$X)) %*% data$X)
  } else if (!is.null(data$XtX)) {
    return(B_mat %*% data$XtX)
  }
  stop("No predictor matrix available for compute_BR.")
}

# =============================================================================
# CORE ALGORITHM COMPONENTS
#
# Key computational functions that implement the mathematical core of the
# SuSiE algorithm. These handle eigen decompositions, posterior computations,
# and log Bayes factor calculations.
#
# Functions: compute_eigen_decomposition, add_eigen_decomposition,
# compute_omega_quantities, scale_design_matrix, compute_theta_blup,
# lbf_stabilization, compute_posterior_weights, compute_lbf_gradient
# =============================================================================

# Compute eigenvalue decomposition for unmappable methods
# When X (low-rank factor) is available, uses thin SVD (O(pB^2)) instead
# of eigen decomposition of XtX (O(p^3)).
#' @keywords internal
compute_eigen_decomposition <- function(XtX, n, X = NULL) {
  if (!is.null(X)) {
    # Thin SVD: V is p x r where r = min(n,p).  No null-space padding -
    # null-space components contribute 0 to every eigenspace sum in the
    # SuSiE-inf math, and storing them as p x p was the memory wall.
    sv  <- svd(X, nu = 0)
    Dsq <- pmax(sv$d^2, 0)
    idx <- order(Dsq, decreasing = TRUE)
    return(list(V = sv$v[, idx, drop = FALSE], Dsq = Dsq[idx], VtXty = NULL))
  }

  LD  <- XtX / n
  eig <- eigen(LD, symmetric = TRUE)
  idx <- order(eig$values, decreasing = TRUE)

  list(
    V     = eig$vectors[, idx],
    Dsq   = pmax(eig$values[idx] * n, 0),
    VtXty = NULL
  )
}

# Add eigen decomposition to a data object for unmappable_effects = "inf".
# Two input shapes are supported:
#   (a) Individual-level data: data$X present, data$XtX absent.  Thin SVD of
#       standardized X; stores p x r eigen factors with NO null-space padding.
#       Skips the O(np^2) XtX formation entirely.
#   (b) Sufficient-statistics data: data$XtX present.  Existing behavior.
#' @keywords internal
add_eigen_decomposition <- function(data, params, individual_data = NULL) {

  if (!is.null(data$X) && is.null(data$XtX)) {
    # (a) Individual-data path.  SVD operates on the standardized X so that
    # V/D match (scaled X)'(scaled X), which is the operator the rest of the
    # SuSiE-inf eigenspace math assumes.  X_std is materialized only for the
    # SVD call and then released.
    cm    <- attr(data$X, "scaled:center")
    csd   <- attr(data$X, "scaled:scale")
    X_std <- scale_design_matrix(data$X, cm, csd)
    sv    <- svd(X_std)
    rm(X_std)

    r <- length(sv$d)
    data$eigen_vectors    <- sv$v                                # p x r
    data$eigen_vectors_sq <- sv$v^2                              # p x r (cached)
    data$eigen_values     <- sv$d^2                              # length r
    data$Uty              <- as.vector(crossprod(sv$u, data$y))  # length r
    data$VtXty            <- sv$d * data$Uty                     # = V' Xty
    data$rank             <- r

    data$Xty <- compute_Xty(data$X, data$y)
    data$yty <- sum(data$y^2)
    return(data)
  }

  # (b) Sufficient-statistics path.
  eigen_decomp <- compute_eigen_decomposition(data$XtX, data$n, X = data$X)
  data$eigen_vectors    <- eigen_decomp$V
  data$eigen_vectors_sq <- eigen_decomp$V^2
  data$eigen_values     <- eigen_decomp$Dsq
  data$VtXty            <- t(eigen_decomp$V) %*% data$Xty   # match reference summation order
  data$rank             <- length(eigen_decomp$Dsq)
  return(data)
}

#' Scale design matrix using centering and scaling parameters
#'
#' Applies column-wise centering and scaling to match the space used by
#' compute_XtX() and compute_Xty() for unmappable effects methods.
#'
#' @param X Matrix to scale (n x p)
#' @param center Vector of column means to subtract (length p), or NULL
#' @param scale Vector of column SDs to divide by (length p), or NULL
#'
#' @return Scaled matrix with centered and scaled columns
#'
#' @keywords internal
scale_design_matrix <- function(X, center = NULL, scale = NULL) {
  if (is.null(center)) center <- rep(0, ncol(X))
  if (is.null(scale)) scale <- rep(1, ncol(X))

  X_centered <- sweep(X, 2, center, "-")
  X_scaled <- sweep(X_centered, 2, scale, "/")

  return(X_scaled)
}

# Compute (scaled X)'(scaled X) v using the cached thin-SVD factors.
# Cost is O(p * r) instead of O(p^2) or O(n p).  Used in the SuSiE-inf
# individual path where compute_Rv with raw data$X would return the
# unscaled X'X v -- here the eigen factors come from svd(X_std), so
# V (D^2 V' v) is the standardized operator that matches compute_XtX().
#' @keywords internal
compute_XtXv_eigen <- function(data, v) {
  Vtv <- crossprod(data$eigen_vectors, v)
  as.vector(data$eigen_vectors %*% (data$eigen_values * Vtv))
}

# Compute Omega-weighted quantities for unmappable effects methods.
# diagXtOmegaX[j] = sum_k V[j,k]^2 * (lambda_k / omega_var_k) is a length-p
# vector.  Uses cached V^2 (eigen_vectors_sq) to skip the inline V^2
# computation but keeps the rowSums(sweep(...)) summation order so output
# is bit-identical to the upstream reference (which forms V^2 inline).
#' @keywords internal
compute_omega_quantities <- function(data, tau2, sigma2) {
  omega_var    <- tau2 * data$eigen_values + sigma2
  V_sq         <- if (!is.null(data$eigen_vectors_sq)) data$eigen_vectors_sq
                  else data$eigen_vectors^2
  diagXtOmegaX <- rowSums(sweep(V_sq, 2, data$eigen_values / omega_var, `*`))
  list(omega_var = omega_var, diagXtOmegaX = diagXtOmegaX)
}

# Compute unmappable effects coefficient vector using BLUP.
# Reads model$omega_var and model$XtOmegay (cached at iteration boundaries by
# initialize_susie_model and update_derived_quantities).  Falls back to
# computing them on the fly if the caches are absent.
#' @keywords internal
compute_theta_blup <- function(data, model) {
  if (is.null(model$omega_var) || is.null(model$XtOmegay)) {
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega_var <- omega_res$omega_var
    XtOmegay  <- as.vector(data$eigen_vectors %*% (data$VtXty / omega_var))
  } else {
    omega_var <- model$omega_var
    XtOmegay  <- model$XtOmegay
  }

  b         <- colSums(model$mu * model$alpha)
  XtOmegaXb <- as.vector(data$eigen_vectors %*%
                           ((t(data$eigen_vectors) %*% b) *
                              data$eigen_values / omega_var))
  XtOmegar  <- XtOmegay - XtOmegaXb
  model$tau2 * XtOmegar
}

# Stabilize log Bayes factors and compute log posterior odds
#' @keywords internal
lbf_stabilization <- function(lbf, prior_weights, shat2) {
  lpo <- lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # When shat2 is infinite, set lbf=0 and lpo to prior (no information from data)
  infinite_idx      <- is.infinite(shat2)
  lbf[infinite_idx] <- 0
  lpo[infinite_idx] <- log(prior_weights[infinite_idx] + sqrt(.Machine$double.eps))

  return(list(lbf = lbf, lpo = lpo))
}

# Compute alpha and lbf for each effect
#' @keywords internal
compute_posterior_weights <- function(lpo) {

  w_weighted     <- exp(lpo - max(lpo))
  weighted_sum_w <- sum(w_weighted)
  alpha          <- w_weighted / weighted_sum_w

  return(list(
    alpha     = alpha,
    lbf_model = log(weighted_sum_w) + max(lpo)
  ))
}

# Compute gradient for prior variance optimization
#' @keywords internal
compute_lbf_gradient <- function(alpha, betahat, shat2, V, use_NIG = FALSE) {
  # No gradient computation for NIG prior
  if (use_NIG) {
    return(NULL)
  }

  T2 <- betahat^2 / shat2
  grad_components <- 0.5 * (1 / (V + shat2)) * ((shat2 / (V + shat2)) * T2 - 1)
  grad_components[is.nan(grad_components)] <- 0
  gradient <- sum(alpha * grad_components)
  return(gradient)
}

# =============================================================================
# VARIANCE ESTIMATION
#
# Functions specifically for estimating variance components using different
# methods (MLE, MoM, NIG). These handle both standard SuSiE
# and unmappable effects models.
#
# Functions: mom_unmappable, mle_unmappable, create_ash_grid,
# compute_lbf_NIG_univariate, compute_lbf_NIG,
# compute_posterior_moments_NIG, update_prior_variance_NIG_EM,
# compute_kl_NIG, inv_gamma_factor, compute_null_loglik_NIG,
# compute_marginal_loglik, est_residual_variance, update_model_variance
# =============================================================================

# Method of Moments variance estimation for unmappable effects methods
#' @keywords internal
mom_unmappable <- function(data, params, model, omega, tau2, est_tau2 = TRUE, est_sigma2 = TRUE) {
  L <- nrow(model$mu)

  A <- matrix(0, nrow = 2, ncol = 2)
  A[1, 1] <- data$n
  A[1, 2] <- sum(data$eigen_values)
  A[2, 1] <- A[1, 2]
  A[2, 2] <- sum(data$eigen_values^2)

  # Compute diag(V'MV).  Final accumulator: rowSums(sweep(t(V)^2, 2, tmpD, *))
  # = sum_j V[j,k]^2 tmpD[j] for each k.  Same summation order as upstream
  # reference, but with V^2 already cached so we skip the inline t(V)^2.
  b <- colSums(model$mu * model$alpha)
  Vtb <- crossprod(data$eigen_vectors, b)
  diagVtMV <- Vtb^2
  tmpD <- rep(0, data$p)

  for (l in seq_len(L)) {
    bl <- model$mu[l, ] * model$alpha[l, ]
    Vtbl <- crossprod(data$eigen_vectors, bl)
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + model$alpha[l, ] * (model$mu[l, ]^2 + 1 / omega[l, ])
  }

  diagVtMV <- diagVtMV + rowSums(sweep(t(data$eigen_vectors_sq), 2, tmpD, `*`))

  # Compute x
  x <- rep(0, 2)
  x[1] <- data$yty - 2 * sum(b * data$Xty) + sum(data$eigen_values * diagVtMV)
  x[2] <- sum(data$Xty^2) - 2 * sum(Vtb * data$VtXty * data$eigen_values) +
    sum(data$eigen_values^2 * diagVtMV)

  if (est_tau2) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigma2 <- sol[1]
      tau2   <- sol[2]
    } else {
      sigma2 <- x[1] / data$n
      tau2   <- 0
    }
    if (params$verbose) {
      message(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigma2, tau2))
    }
  } else if (est_sigma2) {
    sigma2 <- (x[1] - A[1, 2] * tau2) / data$n
    if (params$verbose) {
      message(sprintf("Update sigma^2 to %f\n", sigma2))
    }
  }
  return(list(sigma2 = sigma2, tau2 = tau2))
}

# MLE variance estimation for unmappable effects
#' @keywords internal
mle_unmappable <- function(data, params, model, omega, est_tau2 = TRUE, est_sigma2 = TRUE) {
  L <- nrow(model$alpha)
  r <- length(data$eigen_values)   # rank of the stored eigenspace

  # Set default ranges
  sigma2_range <- c(0.2 * data$yty / data$n, 1.2 * data$yty / data$n)
  tau2_range   <- c(1e-12, 1.2 * data$yty / (data$n * data$p))

  # Compute diag(V'MV) (same form as mom_unmappable -- preserves the upstream
  # summation order while using the cached V^2).
  b        <- colSums(model$mu * model$alpha)
  Vtb      <- crossprod(data$eigen_vectors, b)
  diagVtMV <- Vtb^2
  tmpD     <- rep(0, data$p)

  for (l in seq_len(L)) {
    bl       <- model$mu[l, ] * model$alpha[l, ]
    Vtbl     <- crossprod(data$eigen_vectors, bl)
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD     <- tmpD + model$alpha[l, ] * (model$mu[l, ]^2 + 1 / omega[l, ])
  }

  diagVtMV <- diagVtMV + rowSums(sweep(t(data$eigen_vectors_sq), 2, tmpD, `*`))

  # Negative ELBO as function of x = (sigma^2, tau^2).
  # Use (n - r) instead of (n - p): when eigen_values is stored without
  # null-space padding (length r < p), the dropped null-space components
  # would each contribute 0.5 log(sigma^2), and (n - p) over-counts them.
  # For the legacy full-eigen path r = p, so behavior is unchanged.
  f <- function(x) {
    sigma2_val <- x[1]
    tau2_val   <- x[2]
    var_val    <- tau2_val * data$eigen_values + sigma2_val

    0.5 * (data$n - r) * log(sigma2_val) + 0.5 / sigma2_val * data$yty +
      sum(0.5 * log(var_val) -
            0.5 * tau2_val / sigma2_val * data$VtXty^2 / var_val -
            Vtb * data$VtXty / var_val +
            0.5 * data$eigen_values / var_val * diagVtMV)
  }

  # Negative ELBO for sigma^2 only (when tau^2 is fixed)
  g <- function(sigma2_val) {
    f(c(sigma2_val, model$tau2))
  }

  # Initialize with current values
  sigma2 <- model$sigma2
  tau2 <- model$tau2

  if (est_tau2) {
    # Optimize both sigma^2 and tau^2
    res <- optim(
      par    = c(model$sigma2, model$tau2),
      fn     = f,
      method = "L-BFGS-B",
      lower  = c(sigma2_range[1], tau2_range[1]),
      upper  = c(sigma2_range[2], tau2_range[2])
    )

    if (res$convergence == 0) {
      sigma2 <- res$par[1]
      tau2   <- res$par[2]
      if (params$verbose) {
        message(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigma2, tau2))
      }
    } else {
      warning_message("MLE optimization failed to converge; keeping previous parameters")
    }
  } else if (est_sigma2) {
    # Optimize only sigma^2
    res <- optim(
      par    = model$sigma2,
      fn     = g,
      method = "L-BFGS-B",
      lower  = sigma2_range[1],
      upper  = sigma2_range[2]
    )

    if (res$convergence == 0) {
      sigma2 <- res$par
      if (params$verbose) {
        message(sprintf("Update sigma^2 to %f\n", sigma2))
      }
    } else {
      warning_message("MLE optimization failed to converge; keeping previous parameters")
    }
  }

  return(list(sigma2 = sigma2, tau2 = tau2))
}

# Extract NIG sufficient statistics from model, regardless of data type
# This is the ONLY function that needs to know whether we have individual or SS data.
# All other NIG functions work with (yy, sxy, tau) uniformly.
#' @keywords internal
get_nig_sufficient_stats <- function(data, model) {
  if (!is.null(model$raw_residuals)) {
    # Individual data path: compute from the SER crossproducts. Using cor()
    # would silently center X and the residuals, which is wrong when the
    # caller sets intercept = FALSE.
    yy  <- sum(model$raw_residuals^2)
    sxy <- model$residuals / sqrt(model$predictor_weights * yy)
    sxy <- pmin(pmax(sxy, -1), 1)
    tau <- 1
  } else {
    # SS/RSS path: use pre-computed quantities
    yy  <- model$yy_residual
    sxy <- model$residuals / sqrt(model$predictor_weights * yy)
    # Clamp sxy to [-1, 1]: with approximate R from a finite reference,
    # Cauchy-Schwarz may be violated numerically, giving |sxy| > 1.
    # This would make rss = yy*(1 - r0*sxy^2) negative, producing NaN in log BF.
    sxy <- pmin(pmax(sxy, -1), 1)
    tau <- if (!is.null(model$shat2_inflation)) model$shat2_inflation else 1
  }
  list(yy = yy, sxy = sxy, tau = tau)
}

# Inputs for compute_kl_NIG: single-IG projection of the SER posterior.
# a_post = (alpha0+n)/2; b_post = (beta0 + sum_j alpha_j RSS_j)/2 with
# RSS_j = yy*(1 - r0_j*sxy_j^2); s_j_sq = r0_j * tau / xx_j.
#' @keywords internal
nig_kl_inputs <- function(data, params, model, l) {
  nig_ss  <- get_nig_sufficient_stats(data, model)
  V_l     <- model$V[l]
  r0      <- V_l / (V_l + nig_ss$tau / model$predictor_weights)
  s_j_sq  <- r0 * nig_ss$tau / model$predictor_weights

  alpha_l <- model$alpha[l, ]
  rss_avg <- nig_ss$yy * (1 - sum(alpha_l * r0 * nig_ss$sxy^2))

  list(
    a_post = (params$alpha0 + data$n) / 2,
    b_post = (params$beta0 + rss_avg) / 2,
    s_j_sq = s_j_sq
  )
}

# Compute log Bayes factor for NIG prior (univariate form on raw x, y)
#' @keywords internal
compute_lbf_NIG_univariate <- function(x, y, s0, alpha0 = 0, beta0 = 0) {
  x <- x - mean(x)
  y <- y - mean(y)
  n <- length(x)
  xx <- sum(x * x)
  xy <- sum(x * y)
  yy <- sum(y * y)
  r0 <- s0 / (s0 + 1 / xx)
  sxy <- xy / sqrt(xx * yy)
  ratio <- (beta0 + yy * (1 - r0 * sxy^2)) / (beta0 + yy)
  return((log(1 - r0) - (n + alpha0) * log(ratio)) / 2)
}

# Compute log Bayes factors under Normal-Inverse-Gamma (NIG) prior
#' @keywords internal
compute_lbf_NIG <- function(n, xx, xy, yy, sxy, s0, a0, b0, tau = 1) {
  r0 <- s0 / (s0 + tau / xx)
  rss <- yy * (1 - r0 * sxy^2)

  # Update inverse-gamma parameters
  a1 <- a0 + n
  b1 <- b0 + rss

  # Compute log Bayes factor for each variable
  lbf <- -(log(1 + s0 * xx / tau) + a1 * log(b1 / (b0 + yy))) / 2

  return(lbf)
}

# Compute posterior moments under Normal-Inverse-Gamma (NIG) prior
#' @keywords internal
compute_posterior_moments_NIG <- function(n, xx, xy, yy, sxy, s0, a0, b0, tau = 1) {
  r0 <- s0 / (s0 + tau / xx)
  rss <- yy * (1 - r0 * sxy^2)

  # Update inverse-gamma parameters
  a1 <- a0 + n
  b1 <- b0 + rss

  # Compute least-squares estimate for each variable
  bhat <- xy / xx

  # Compute posterior mean
  post_mean <- r0 * bhat

  # Compute posterior variance
  post_var <- b1 / (a1 - 2) * r0 * tau / xx

  # Posterior mean of residual variance under IG((a0+n)/2, (b0+RSS)/2)
  rv <- (b1 / 2) / (a1 / 2 - 1)

  return(list(
    post_mean  = post_mean,
    post_mean2 = post_var + post_mean^2,
    post_var   = post_var,
    s_j_sq     = r0 * tau / xx,
    rv         = rv
  ))
}

# EM update for prior variance under Normal-Inverse-Gamma (NIG) prior
#' @keywords internal
update_prior_variance_NIG_EM <- function(n, xx, xy, yy, sxy, pip, s0, a0, b0, tau = 1) {
  r0   <- s0 / (s0 + tau / xx)
  rss  <- yy * (1 - r0 * sxy^2)

  # Update inverse-gamma parameters
  a1   <- a0 + n
  b1   <- b0 + rss

  # Compute posterior mean and variance component
  bhat <- xy / xx
  post_mean  <- r0 * bhat
  post_var   <- r0 * tau / xx

  u  <- gamma(1/2) / beta(a1/2, 1/2)
  mb <- post_mean * sqrt(2 / b1) * u
  vb <- post_var + post_mean^2 * 2 / b1 * (1 / beta(a1/2, 1) - u^2)

  return(sum(pip * (vb + mb^2)))
}

# KL divergence for the NIG variational form. KL_beta uses the Gaussian-
# Gaussian KL with shared sigma^2 scaling (sigma^2 log-dets cancel; only
# mu^2/(sigma^2 V) survives). KL_sigma2 is the closed form for IG||IG.
#' @keywords internal
compute_kl_NIG <- function(alpha, post_mean, post_mean2, pi, V, a0, b0, a_post, b_post,
                           s_j_sq) {
  eps <- .Machine$double.eps

  # KL for categorical assignment q(gamma) || p(gamma)
  KL_gamma <- sum(alpha * (log(pmax(alpha, eps)) - log(pmax(pi, eps))))

  # KL for b given sigma^2, integrated over q(sigma^2)
  E_inv_sigma2 <- a_post / b_post
  KL_beta <- 0.5 * sum(alpha * (
    log(V) - log(pmax(s_j_sq, eps)) + s_j_sq / V +
      post_mean^2 * E_inv_sigma2 / V - 1
  ))

  # KL between IG posterior and IG prior (closed form)
  KL_sigma2 <- lgamma(a0) - lgamma(a_post) +
    a0 * log(b_post / b0) +
    (a_post - a0) * digamma(a_post) -
    a_post + (a_post * b0) / b_post

  return(as.numeric(KL_gamma + KL_beta + KL_sigma2))
}

# Compute log-normalizing factor for the IG(a,b) distribution
#' @keywords internal
inv_gamma_factor <- function(a, b) {
  return(a * log(b) - lgamma(a))
}

# Compute null log-likelihood under NIG prior
#' @keywords internal
compute_null_loglik_NIG <- function(n, yy, a0, b0, use_NIG = FALSE) {
  # No null log-likelihood for non-NIG prior
  if (!use_NIG) {
    return(NULL)
  }

  return(-n * log(2 * pi) / 2 +
         inv_gamma_factor(a0 / 2, b0 / 2) -
         inv_gamma_factor((a0 + n) / 2, (b0 + yy) / 2))
}

# Compute marginal log-likelihood for single effect regression
#' @keywords internal
compute_marginal_loglik <- function(lbf_model, n, yy, a0, b0, use_NIG = FALSE) {
  # No marginal log-likelihood computation for non-NIG prior
  if (!use_NIG) {
    return(NULL)
  }

  ll0 <- compute_null_loglik_NIG(n, yy, a0, b0, use_NIG = TRUE)
  return(lbf_model + ll0)
}

# Estimate residual variance
#' @keywords internal
est_residual_variance <- function(data, model) {
  resid_var <- (1 / data$n) * get_ER2(data, model)
  if (resid_var < 0) {
    stop("est_residual_variance() failed: the estimated value is negative")
  }
  return(resid_var)
}


# =============================================================================
# SuSiE-ASH SHARED UTILITIES
#
# Functions shared between individual-level and summary-statistics ash paths.
# The masking logic is data-type-agnostic (only needs correlation matrix),
# while the mr.ash fitting is dispatched to either mr.ash (individual) or
# mr.ash.rss (summary stats).
#
# Design: individual-level data uses mr.ash directly, SS data uses mr.ash.rss.
# The init/masking/cleanup logic is shared since it depends only on model
# dimensions (n, p, L), not on data representation.
# =============================================================================

# Initialize ash tracking fields on a model object
#
# Shared between individual and SS model initialization.
# Individual models store X_theta (n-vector); SS models store XtX_theta (p-vector).
#
# @param model Model object to augment
# @param n Number of observations (used only by individual path for X_theta)
# @param p Number of predictors
# @param L Number of single effects
# @param is_individual Whether this is individual-level data
#
# @return Model with ash tracking fields added
#
# @keywords internal
init_ash_fields <- function(model, n, p, L, is_individual = FALSE) {
  model$tau2     <- 0
  model$theta    <- rep(0, p)
  model$ash_iter <- 0
  model$ash_pi   <- NULL
  model$ash_s0   <- NULL

  if (is_individual) {
    model$X_theta <- rep(0, n)
  } else {
    model$XtX_theta <- rep(0, p)
  }

  return(model)
}

init_ash_fields_filter_archived <- function(model, n, p, L, is_individual = FALSE) {
  model <- init_ash_fields(model, n, p, L, is_individual)

  # Additional tracking fields for the archived filter-based masking
  model$masked                 <- rep(FALSE, p)
  model$diffuse_iter_count     <- rep(0, L)
  model$prev_sentinel          <- rep(0, L)
  model$unmask_candidate_iters <- rep(0, p)
  model$ever_unmasked          <- rep(FALSE, p)
  model$force_exposed_iter     <- rep(0, p)
  model$ever_diffuse           <- rep(0, L)
  model$second_chance_used     <- rep(FALSE, p)
  model$prev_case              <- rep(0, L)

  return(model)
}

# Update ash variance components: coordinate SuSiE and Mr.ASH
#
# Manages the interaction between SuSiE's sparse effects (beta) and
# Mr.ASH's dense effects (theta). Uses three zones to prevent Mr.ASH
# from absorbing signals SuSiE is fine-mapping:
#
#   SUBTRACT: Trusted effects are fully removed from Mr.ASH residuals.
#     Mr.ASH cannot see or absorb these confirmed signals.
#   MASK: Emerging effects have theta forced to 0 in their LD
#     neighborhoods. Mr.ASH cannot absorb signal here, giving SuSiE
#     time to resolve these effects.
#   EXPOSE: Positions not covered by subtraction or masking. Mr.ASH
#     absorbs signal freely. LD spillover and false positives that
#     get flagged and unmasked end up here.
#
# Spillover is detected in two ways:
#   Across-SER: Two effects whose sentinels are in tight LD (|r| > 0.9)
#     are colliding on the same signal. Both get flagged.
#   Within-SER: A sentinel that jumps to a distant variant (|r| < 0.5)
#     indicates the effect is not well-localized. Gets flagged.
#
# Flagged effects are initially masked (protected). After unmask_delay
# iterations, the OLD sentinel's LD neighborhood is selectively exposed
# to Mr.ASH, allowing absorption of spillover at the confusion point
# while the effect's current position stays protected.
#
# @param data Data object (individual or SS).
# @param model Current SuSiE model.
# @param params Parameters object.
#
# @return List of updated model fields (sigma2, tau2, theta, etc.)
#   to be merged into model via modifyList.
#
# @keywords internal
## V0-faithful three-case classification for BB+ash filter
## Replaces the body of update_ash_variance_components()
## Key change: standard purity (not effect_purity) for case classification
## + force_mask for diffuse slots' sentinel LD (from V0)

update_ash_variance_components <- function(data, model, params) {

  # BB+ash filter: V0's 3-tier classifier (diffuse / uncertain / confident)
  # with c_hat marginalization of CASE 3 subtraction.
  # Output: b_confident (subtract from mr.ash residual) and mask.

  # --- User parameters ---
  purity_threshold <- if (!is.null(params$purity_threshold)) params$purity_threshold else 0.5
  pip_threshold    <- if (!is.null(params$pip_threshold))    params$pip_threshold    else 0.1
  ld_threshold     <- if (!is.null(params$ld_threshold))     params$ld_threshold     else 0.5

  # --- Internal constants ---
  diffuse_purity             <- 0.1   # purity below this = CASE 1
  cs_threshold               <- 0.9   # working CS coverage
  neighborhood_pip_threshold <- 0.4   # LD-spread mask threshold
  collision_threshold        <- 0.9   # strong LD = same signal
  tight_ld_threshold         <- 0.95  # WTE exposure region
  diffuse_iter_count         <- 2L    # CASE 2 stable iters before WTE
  second_chance_wait         <- 3L    # iters after expose before re-mask
  delayed_unmask_iter        <- 2L    # iters unwanted before unmask

  is_individual <- inherits(data, "individual")
  L <- nrow(model$alpha)
  p <- ncol(model$alpha)
  if (is.null(model$ash_iter)) model$ash_iter <- 0L
  model$ash_iter <- model$ash_iter + 1L

  # Lazy-init state
  if (is.null(model$prev_case))              model$prev_case <- rep(0L, L)
  if (is.null(model$prev_sentinel))          model$prev_sentinel <- rep(0L, L)
  if (is.null(model$ever_diffuse))           model$ever_diffuse <- rep(0L, L)
  if (is.null(model$diffuse_iter_count))     model$diffuse_iter_count <- rep(0L, L)
  if (is.null(model$masked))                 model$masked <- rep(FALSE, p)
  if (is.null(model$ever_unmasked))          model$ever_unmasked <- rep(FALSE, p)
  if (is.null(model$unmask_candidate_iters)) model$unmask_candidate_iters <- rep(0L, p)
  if (is.null(model$force_exposed_iter))     model$force_exposed_iter <- rep(0L, p)
  if (is.null(model$second_chance_used))     model$second_chance_used <- rep(FALSE, p)

  xcorr_result <- get_xcorr(data)
  Xcorr <- xcorr_result$Xcorr
  data <- xcorr_result$data

  # c_hat = BB posterior slot-active weight (fallback 1 = V0).
  c_hat <- if (!is.null(model$slot_weights)) model$slot_weights else rep(1, L)

  # First pass: sentinels, purity, activity
  sentinels     <- apply(model$alpha, 1, which.max)
  effect_purity <- rep(1.0, L)
  is_active     <- rep(FALSE, L)
  for (l in 1:L) {
    a <- model$alpha[l, ]
    is_active[l] <- (max(a) - min(a)) >= 5e-5
    alpha_order <- order(a, decreasing = TRUE)
    cs_size <- min(sum(cumsum(a[alpha_order]) <= cs_threshold) + 1L, p)
    if (cs_size > 1) {
      cs_indices <- alpha_order[1:cs_size]
      R_cs <- abs(Xcorr[cs_indices, cs_indices])
      effect_purity[l] <- min(R_cs[upper.tri(R_cs)])
    }
  }

  # Collision: sentinels of two active slots in strong LD -> bump ever_diffuse
  current_collision <- rep(FALSE, L)
  for (l in which(is_active)) {
    others <- setdiff(which(is_active), l)
    if (length(others) == 0) next
    if (any(abs(Xcorr[sentinels[l], sentinels[others]]) > collision_threshold)) {
      current_collision[l] <- TRUE
      model$ever_diffuse[l] <- model$ever_diffuse[l] + 1L
    }
  }

  b_confident     <- rep(0, p)
  alpha_protected <- matrix(0, nrow = L, ncol = p)
  force_unmask    <- rep(FALSE, p)
  force_mask      <- rep(FALSE, p)

  # Second pass: classify slots
  current_case <- rep(0L, L)
  for (l in 1:L) {
    purity   <- effect_purity[l]
    sentinel <- sentinels[l]

    # Sentinel-jump reset: if sentinel moved to a non-tight-LD position, reset
    # the CASE 2 stability counter so WTE doesn't fire on a newly-picked signal.
    if (sentinel != model$prev_sentinel[l] && model$prev_sentinel[l] > 0L) {
      if (abs(Xcorr[sentinel, model$prev_sentinel[l]]) < tight_ld_threshold) {
        model$diffuse_iter_count[l] <- 0L
      }
    }

    can_be_confident <- (purity >= purity_threshold) && (model$ever_diffuse[l] == 0L)

    if (purity < diffuse_purity) {
      # CASE 1: diffuse. Narrow protection + force_mask on sentinel LD.
      current_case[l] <- 1L
      model$diffuse_iter_count[l] <- 0L
      moderate_ld <- abs(Xcorr[sentinel, ]) > ld_threshold
      to_protect  <- moderate_ld | (model$alpha[l, ] > 5 / p)
      alpha_protected[l, to_protect] <- model$alpha[l, to_protect]
      force_mask <- force_mask | moderate_ld
    } else if (!can_be_confident) {
      # CASE 2: uncertain. Collision -> reset counter only (no alpha protection).
      # Otherwise full alpha protection + wait-then-expose.
      current_case[l] <- 2L
      if (current_collision[l]) {
        model$diffuse_iter_count[l] <- 0L
      } else {
        model$diffuse_iter_count[l] <- model$diffuse_iter_count[l] + 1L
        alpha_protected[l, ] <- model$alpha[l, ]
        if (model$diffuse_iter_count[l] >= diffuse_iter_count) {
          tight_ld <- abs(Xcorr[sentinel, ]) > tight_ld_threshold
          expose   <- tight_ld & !model$second_chance_used
          newly    <- expose & (model$force_exposed_iter == 0L)
          model$force_exposed_iter[newly] <- model$ash_iter
          if (any(newly)) model$diffuse_iter_count[l] <- 0L
          alpha_protected[l, expose] <- 0
          force_unmask <- force_unmask | expose
        }
      }
    } else {
      # CASE 3: confident. Full protection + c_hat-weighted subtraction.
      current_case[l] <- 3L
      model$diffuse_iter_count[l] <- 0L
      alpha_protected[l, ] <- model$alpha[l, ]
      b_confident <- b_confident + c_hat[l] * model$alpha[l, ] * model$mu[l, ]
    }
  }

  # Oscillation: slot flipping between CASE 2 and CASE 3 is unstable -> mark
  # sticky-diffuse, and if it landed on CASE 3 this iter, reverse the
  # subtraction we just added (slot is not trustworthy yet).
  oscillated    <- model$prev_case != 0L & current_case != 0L &
                   ((model$prev_case == 2L & current_case == 3L) |
                    (model$prev_case == 3L & current_case == 2L))
  unstable_case3 <- oscillated & (current_case == 3L)
  model$ever_diffuse[oscillated] <- model$ever_diffuse[oscillated] + 1L
  for (l in which(unstable_case3)) {
    b_confident <- b_confident - c_hat[l] * model$alpha[l, ] * model$mu[l, ]
  }
  model$prev_case     <- current_case
  model$prev_sentinel <- sentinels

  # Mask: PIP-based union, with per-position persistence
  pip_protected    <- susie_get_pip(alpha_protected)
  LD_adj           <- abs(Xcorr) > ld_threshold
  neighborhood_pip <- as.vector(LD_adj %*% pip_protected)
  want_masked <- (neighborhood_pip > neighborhood_pip_threshold) |
                 (pip_protected    > pip_threshold) |
                 force_mask

  # Delayed unmask: count iters a masked position is no longer wanted
  reset_idx <- want_masked | !model$masked
  model$unmask_candidate_iters[!reset_idx] <- model$unmask_candidate_iters[!reset_idx] + 1L
  model$unmask_candidate_iters[reset_idx]  <- 0L
  ready_to_unmask <- model$masked &
    ((model$unmask_candidate_iters >= delayed_unmask_iter & !model$ever_unmasked) | force_unmask)

  model$ever_unmasked[ready_to_unmask] <- TRUE
  masked <- (model$masked | want_masked) & !ready_to_unmask & !model$ever_unmasked

  # Second chance
  should_restore <- (model$force_exposed_iter > 0L) &
                    (model$ash_iter - model$force_exposed_iter >= second_chance_wait) &
                    !model$second_chance_used
  if (any(should_restore)) {
    model$second_chance_used[should_restore] <- TRUE
    model$force_exposed_iter[should_restore] <- 0L
    model$ever_unmasked[should_restore]      <- FALSE
    masked[should_restore]                   <- TRUE
  }
  model$masked <- masked
  mask <- masked

  # Mr.ASH fit
  .skip_mrash <- getOption("susie.skip_mrash", FALSE)
  model$theta[mask] <- 0
  if (.skip_mrash) {
    theta_new <- model$theta; theta_new[mask] <- 0
    ash_result <- list(
      beta = theta_new,
      sigma2 = if (!is.null(model$sigma2)) model$sigma2 else 1,
      pi = model$ash_pi,
      tau2 = if (!is.null(model$tau2)) model$tau2 else 0)
  } else {
    convtol <- if (model$ash_iter < 2) 1e-3 else 1e-4
    if (is_individual) {
      ash_result <- compute_ash_from_individual_data(
        data$X, data$y, b_confident, model, params, convtol)
    } else {
      ash_result <- compute_ash_from_summary_stats(
        data, b_confident, model, params, convtol)
    }
    theta_new <- ash_result$beta
    theta_new[mask] <- 0
  }

  # Diagnostic
  .ash_debug <- TRUE
  if (.ash_debug) {
    model$ever_uncertain <- model$ever_diffuse > 0
    diag_df <- diagnose_bb_ash_iter(
      model, Xcorr, mask, b_confident,
      sentinels, current_collision,
      current_case == 3L, current_case == 3L,
      which(current_case == 2L), which(model$V > 0), c_hat,
      list(beta = theta_new, sigma2 = ash_result$sigma2, pi = ash_result$pi),
      p,
      high_chat = integer(0), low_chat = integer(0),
      collision_threshold = collision_threshold,
      purity_threshold = ld_threshold,
      masking_threshold = ld_threshold,
      nPIP_threshold = pip_threshold,
      c_hat_excess_threshold = NA,
      alpha_entropy_threshold = NA,
      slot_prior = params$slot_prior,
      mask_smoothness = effect_purity,
      mask_amount = c_hat,
      mask_concentration = current_case,
      mask_burnin = model$ever_diffuse,
      mask_spread_pip_at_sent = pip_protected[sentinels],
      mask_pip_prot_at_sent = pip_protected[sentinels])
    if (is.null(model$.diag_env)) model$.diag_env <- new.env(parent = emptyenv())
    if (is.null(model$.diag_env$history)) model$.diag_env$history <- list()
    model$.diag_env$history[[length(model$.diag_env$history) + 1]] <- diag_df
  }

  result <- list(
    sigma2                 = ash_result$sigma2,
    tau2                   = ash_result$tau2,
    theta                  = theta_new,
    ash_pi                 = ash_result$pi,
    ash_iter               = model$ash_iter,
    prev_case              = model$prev_case,
    prev_sentinel          = model$prev_sentinel,
    ever_diffuse           = model$ever_diffuse,
    diffuse_iter_count     = model$diffuse_iter_count,
    masked                 = model$masked,
    ever_unmasked          = model$ever_unmasked,
    unmask_candidate_iters = model$unmask_candidate_iters,
    force_exposed_iter     = model$force_exposed_iter,
    second_chance_used     = model$second_chance_used,
    .diag_env              = model$.diag_env
  )

  if (is_individual) {
    result$X_theta <- as.vector(data$X %*% theta_new)
  } else {
    result$XtX_theta <- as.vector(compute_Rv(data, theta_new))
  }

  return(result)
}

# Archived: Original filter-based ash variance component update
#
# Performs the full ash update cycle with the Diffusion-Aware masking
# heuristic: get correlation matrix, compute masking (3-case classification,
# WaitThenExpose, collision detection, nPIP masking), fit mr.ash, mask theta.
# This is the original SuSiE-ASH implementation with 9+ tuning parameters
# and per-effect/per-variant tracking arrays.
#
# Kept for internal diagnostics/archiving via unmappable_effects="ash_filter_archived".
# The default "ash" path now uses update_ash_variance_components() (c_hat + 3 LD rules).
#
# @param data Data object (individual or SS)
# @param model Current SuSiE model
# @param params Parameters object
#
# @return List with sigma2, tau2, theta, fitted theta, and all tracking fields
#
# @keywords internal
update_ash_variance_components_filter_archived <- function(data, model, params) {
  is_individual <- inherits(data, "individual")

  # Step 1: Get correlation matrix (cached after first call)
  xcorr_result <- get_xcorr(data)
  Xcorr <- xcorr_result$Xcorr
  data <- xcorr_result$data

  # Step 2: Compute masking (shared 3-case classification)
  mask_result <- compute_ash_masking(Xcorr, model, params)
  b_confident <- mask_result$b_confident
  masked <- mask_result$masked
  model <- mask_result$model

  # Step 3: Fit Mr.ASH (dispatch to individual or SS backend)
  # Set options(susie.skip_mrash = TRUE) to diagnose without mr.ash.
  .skip_mrash <- getOption("susie.skip_mrash", FALSE)
  if (.skip_mrash) {
    p_v0 <- ncol(model$alpha)
    mrash_output <- list(
      beta = model$theta, sigma2 = if (!is.null(model$sigma2)) model$sigma2 else 1,
      pi = model$ash_pi, tau2 = if (!is.null(model$tau2)) model$tau2 else 0,
      sa2 = if (!is.null(model$ash_s0)) model$ash_s0 else 0)
    mrash_output$beta[masked] <- 0
  } else {
    convtol <- if (model$ash_iter < 2) 1e-3 else 1e-4
    if (is_individual) {
      mrash_output <- compute_ash_from_individual_data(
        data$X, data$y, b_confident, model, params, convtol
      )
    } else {
      mrash_output <- compute_ash_from_summary_stats(
        data, b_confident, model, params, convtol
      )
    }
  }

  # Step 4: Zero out theta for masked variants
  theta_new <- mrash_output$beta
  theta_new[masked] <- 0

  # V0 diagnostic: capture data.frame and accumulate on model
  .ash_debug <- TRUE
  if (.ash_debug) {
    diag_df <- diagnose_ash_filter_archived_iter(
      model, Xcorr, masked, b_confident,
      mask_result$sentinels, mask_result$effect_purity,
      mask_result$current_case, mask_result$current_collision,
      mrash_output)
    if (is.null(model$.diag_env)) model$.diag_env <- new.env(parent = emptyenv())
    if (is.null(model$.diag_env$history)) model$.diag_env$history <- list()
    model$.diag_env$history[[length(model$.diag_env$history) + 1]] <- diag_df
  }

  # Step 5: Compute fitted theta (data-representation-specific)
  result <- list(
    sigma2              = mrash_output$sigma2,
    tau2                = mrash_output$tau2,
    theta               = theta_new,
    ash_pi              = mrash_output$pi,
    sa2                 = mrash_output$sa2,
    ash_iter            = model$ash_iter,
    diffuse_iter_count  = model$diffuse_iter_count,
    prev_sentinel       = model$prev_sentinel,
    masked              = masked,
    unmask_candidate_iters = model$unmask_candidate_iters,
    ever_unmasked       = model$ever_unmasked,
    force_exposed_iter  = model$force_exposed_iter,
    second_chance_used  = model$second_chance_used,
    ever_diffuse        = model$ever_diffuse,
    prev_case           = model$prev_case,
    .diag_env           = model$.diag_env
  )

  if (is_individual) {
    result$X_theta <- as.vector(data$X %*% theta_new)
  } else {
    result$XtX_theta <- as.vector(compute_Rv(data, theta_new))
  }

  return(result)
}

# Remove ash-specific runtime fields from model
#
# Shared between cleanup_model.individual() and cleanup_model.ss().
#
# @param model Model object
#
# @return Model with ash runtime fields removed
#
# @keywords internal
cleanup_ash_fields <- function(model) {
  # Remove internal tracking fields from the new ash path.
  # Keep: tau2, theta, ash_pi (user-visible results)
  for (field in c("X_theta", "XtX_theta", "ash_iter", "ash_s0",
                   "ever_uncertain", "prev_sentinel")) {
    model[[field]] <- NULL
  }
  return(model)
}

cleanup_ash_fields_filter_archived <- function(model) {
  # Remove internal tracking fields from the archived filter path.
  for (field in c("X_theta", "XtX_theta", "ash_iter", "ash_s0",
                   "masked", "diffuse_iter_count", "prev_sentinel",
                   "unmask_candidate_iters", "ever_unmasked",
                   "force_exposed_iter", "ever_diffuse",
                   "second_chance_used", "prev_case")) {
    model[[field]] <- NULL
  }
  return(model)
}

# Get or compute correlation matrix for ash masking
#
# For SS data: derives from XtX via safe_cov2cor (cheap scaling).
# For individual data: computes cor(X) and caches it on the data object.
#
# @param data Data object
#
# @return List with Xcorr and (possibly updated) data object
#
# @keywords internal
get_xcorr <- function(data) {
  # Check for cached correlation matrix
  if (!is.null(data$Xcorr_cache)) {
    return(list(Xcorr = data$Xcorr_cache, data = data))
  }

  if (!is.null(data$XtX)) {
    # SS path: derive from XtX
    if (any(!(diag(data$XtX) %in% c(0, 1)))) {
      Xcorr <- safe_cov2cor(data$XtX)
    } else {
      Xcorr <- data$XtX
    }
  } else if (!is.null(data$X)) {
    # Individual path: compute correlation from X
    Xcorr <- safe_cor(data$X)
  } else {
    stop("Cannot compute correlation matrix: data has neither XtX nor X")
  }

  # Cache for future iterations
  data$Xcorr_cache <- Xcorr
  return(list(Xcorr = Xcorr, data = data))
}

# Compute ash masking: classify effects, determine protection/masking
#
# KEY INSIGHT: Protect SuSiE's sparse effects from Mr.ASH absorption,
# but let Mr.ASH absorb unmappable and unreliable signals.
#
# Two types of diffusion (both indicate unreliable effects):
#   1. WITHIN-EFFECT: Low purity - spread across variants not in tight LD
#   2. CROSS-EFFECT: Sentinel collision - multiple effects compete for
#      same position (composite signal, not clean single causal)
#
# Classification into three cases:
#   CASE 1 (diffuse): purity < 0.1 - protect neighborhood loosely
#   CASE 2 (uncertain): low purity OR ever_diffuse - expose to Mr.ASH
#   CASE 3 (confident): purity >= 0.5 AND never diffuse - subtract from residuals
#
# Cross-effect diffusion tracking:
#   - Detect via sentinel collision (sentinels in tight LD across effects)
#   - Mark effect as ever_diffuse (sticky, effect-level)
#   - ever_diffuse effects get zero protection permanently
#   - Real signals survive Mr.ASH competition; composites get absorbed
#
# Low purity (non-diffuse) effects:
#   - Use wait-then-expose mechanism
#   - After diffuse_iter_count stable iterations, expose sentinel neighborhood
#   - If Mr.ASH absorbs it, the signal was synthetic
#   - Second-chance mechanism: after wait period, restore masking to check
#
# @param Xcorr Correlation matrix (p x p)
# @param model Current SuSiE model (with alpha, mu, tracking fields)
# @param params Parameters object
#
# @return List with b_confident, masked, and updated model tracking fields
#
# @keywords internal
compute_ash_masking <- function(Xcorr, model, params) {
  # --- Protection thresholds ---
  neighborhood_pip_threshold <- if (!is.null(params$neighborhood_pip_threshold)) params$neighborhood_pip_threshold else 0.4
  direct_pip_threshold <- if (!is.null(params$direct_pip_threshold)) params$direct_pip_threshold else 0.1
  signal_separation_ld <- if (!is.null(params$signal_separation_ld)) params$signal_separation_ld else 0.5

  # --- Purity thresholds ---
  cs_threshold <- if (!is.null(params$working_cs_threshold)) params$working_cs_threshold else 0.9
  cs_formation_threshold <- if (!is.null(params$cs_formation_threshold)) params$cs_formation_threshold else 0.1
  purity_threshold <- if (!is.null(params$purity_threshold)) params$purity_threshold else 0.5

  # --- LD thresholds for collision and exposure ---
  collision_threshold <- if (!is.null(params$collision_threshold)) params$collision_threshold else 0.9
  tight_ld_threshold <- if (!is.null(params$tight_ld_threshold)) params$tight_ld_threshold else 0.95

  # --- Iteration counters for CASE 2 ---
  diffuse_iter_count <- if (!is.null(params$diffuse_iter_count)) params$diffuse_iter_count else 2
  track_sentinel <- if (!is.null(params$track_sentinel)) params$track_sentinel else TRUE

  # --- Second chance mechanism ---
  second_chance_wait <- if (!is.null(params$second_chance_wait)) params$second_chance_wait else 3

  # --- Unmasking stability ---
  delayed_unmask_iter <- 2

  L <- nrow(model$alpha)
  p <- ncol(model$alpha)
  model$ash_iter <- model$ash_iter + 1

  # =========================================================================
  # First pass: Compute sentinels and purity
  # =========================================================================
  sentinels <- apply(model$alpha, 1, which.max)
  effect_purity <- rep(NA, L)

  for (l in 1:L) {
    alpha_order <- order(model$alpha[l,], decreasing = TRUE)
    cumsum_alpha <- cumsum(model$alpha[l, alpha_order])
    cs_size <- sum(cumsum_alpha <= cs_threshold) + 1
    cs_indices <- alpha_order[1:min(cs_size, p)]
    effect_purity[l] <- get_purity(cs_indices, X = NULL, Xcorr = Xcorr, use_rfast = FALSE)[1]
  }

  # Detect current collision and update ever_diffuse
  current_collision <- rep(FALSE, L)
  for (l in 1:L) {
    if (max(model$alpha[l,]) - min(model$alpha[l,]) < 5e-5) next
    sentinel_l <- sentinels[l]
    for (other_l in (1:L)[-l]) {
      if (max(model$alpha[other_l,]) - min(model$alpha[other_l,]) < 5e-5) next
      if (abs(Xcorr[sentinel_l, sentinels[other_l]]) > collision_threshold) {
        current_collision[l] <- TRUE
      }
    }
    model$ever_diffuse[l] <- model$ever_diffuse[l] + current_collision[l]
  }

  # Initialize per-iteration outputs
  b_confident <- rep(0, p)
  alpha_protected <- matrix(0, nrow = L, ncol = p)
  force_unmask <- rep(FALSE, p)
  force_mask <- rep(FALSE, p)

  # =========================================================================
  # Second pass: Classify effects and determine protection
  # =========================================================================
  current_case <- rep(0, L)

  for (l in 1:L) {
    purity <- effect_purity[l]
    sentinel <- sentinels[l]

    if (track_sentinel && sentinel != model$prev_sentinel[l] && model$prev_sentinel[l] > 0) {
      if (abs(Xcorr[sentinel, model$prev_sentinel[l]]) < tight_ld_threshold) {
        model$diffuse_iter_count[l] <- 0
      }
    }

    is_ever_diffuse <- model$ever_diffuse[l] > 0
    can_be_confident <- (purity >= purity_threshold) && !is_ever_diffuse

    if (purity < cs_formation_threshold) {
      # CASE 1: Diffuse within effect
      current_case[l] <- 1
      model$diffuse_iter_count[l] <- 0
      moderate_ld_with_sentinel <- abs(Xcorr[sentinel,]) > signal_separation_ld
      meaningful_alpha <- model$alpha[l,] > 5/p
      to_protect <- moderate_ld_with_sentinel | meaningful_alpha
      alpha_protected[l, to_protect] <- model$alpha[l, to_protect]
      force_mask <- force_mask | moderate_ld_with_sentinel
    } else if (!can_be_confident) {
      # CASE 2: Uncertain
      current_case[l] <- 2
      if (current_collision[l]) {
        model$diffuse_iter_count[l] <- 0
      } else {
        model$diffuse_iter_count[l] <- model$diffuse_iter_count[l] + 1
        if (model$diffuse_iter_count[l] >= diffuse_iter_count) {
          tight_ld_with_sentinel <- abs(Xcorr[sentinel,]) > tight_ld_threshold
          newly_exposed <- tight_ld_with_sentinel &
                          !model$second_chance_used &
                          (model$force_exposed_iter == 0)
          model$force_exposed_iter[newly_exposed] <- model$ash_iter
          if (any(newly_exposed)) {
            model$diffuse_iter_count[l] <- 0
          }
          alpha_protected[l,] <- model$alpha[l,]
          expose_positions <- tight_ld_with_sentinel & !model$second_chance_used
          alpha_protected[l, expose_positions] <- 0
          force_unmask <- force_unmask | expose_positions
        } else {
          alpha_protected[l,] <- model$alpha[l,]
        }
      }
    } else {
      # CASE 3: Confident
      current_case[l] <- 3
      model$diffuse_iter_count[l] <- 0
      b_confident <- b_confident + model$alpha[l,] * model$mu[l,]
      alpha_protected[l,] <- model$alpha[l,]
    }

    model$prev_sentinel[l] <- sentinel
  }

  # =========================================================================
  # Oscillation detection
  # =========================================================================
  for (l in 1:L) {
    prev <- model$prev_case[l]
    curr <- current_case[l]
    if (curr == 0 || prev == 0) next
    if ((prev == 2 && curr == 3) || (prev == 3 && curr == 2)) {
      model$ever_diffuse[l] <- model$ever_diffuse[l] + 1
      if (curr == 3) {
        b_confident <- b_confident - model$alpha[l,] * model$mu[l,]
      }
    }
  }
  model$prev_case <- current_case

  # =========================================================================
  # Masking logic
  # =========================================================================
  pip_protected <- susie_get_pip(alpha_protected)

  LD_adj <- abs(Xcorr) > signal_separation_ld
  neighborhood_pip <- as.vector(LD_adj %*% pip_protected)
  want_masked <- (neighborhood_pip > neighborhood_pip_threshold) |
                 (pip_protected > direct_pip_threshold) |
                 force_mask

  dont_want_mask <- !want_masked
  model$unmask_candidate_iters[model$masked & dont_want_mask] <-
    model$unmask_candidate_iters[model$masked & dont_want_mask] + 1
  model$unmask_candidate_iters[want_masked | !model$masked] <- 0

  ready_to_unmask <- (model$masked &
                     (model$unmask_candidate_iters >= delayed_unmask_iter) &
                     !model$ever_unmasked) |
                     (model$masked & force_unmask)

  model$ever_unmasked[ready_to_unmask] <- TRUE
  masked <- (model$masked | want_masked) & !ready_to_unmask & !model$ever_unmasked

  # Second chance
  waited_long_enough <- (model$force_exposed_iter > 0) &
                        (model$ash_iter - model$force_exposed_iter) >= second_chance_wait
  should_restore <- waited_long_enough & !model$second_chance_used
  if (any(should_restore)) {
    model$second_chance_used[should_restore] <- TRUE
    model$force_exposed_iter[should_restore] <- 0
    model$ever_unmasked[should_restore] <- FALSE
    masked[should_restore] <- TRUE
  }

  list(
    b_confident = b_confident,
    masked = masked,
    model = model,
    sentinels = sentinels,
    effect_purity = effect_purity,
    current_case = current_case,
    current_collision = current_collision
  )
}

# Run Mr.ASH on individual-level data
#
# Computes residuals from raw X, y and calls mr.ash directly.
#
# @param X Design matrix (n x p)
# @param y Response vector (n)
# @param b_confident Vector of confident effects to subtract from residuals
# @param model Current SuSiE model
# @param params Parameters object
# @param convtol Convergence tolerance for mr.ash
#
# @return List with beta, sigma2, pi, sa2, tau2
#
# @keywords internal
compute_ash_from_individual_data <- function(X, y, b_confident, model, params, convtol = 1e-4) {
  residuals <- y - X %*% b_confident

  mrash_output <- mr.ash(
    X             = X,
    y             = residuals,
    intercept     = FALSE,
    standardize   = FALSE,
    sigma2        = model$sigma2,
    update.sigma2 = params$estimate_residual_variance,
    beta.init     = model$theta,
    pi            = model$ash_pi,
    tol           = list(convtol = convtol, epstol = 1e-12),
    verbose       = params$verbose,
    max.iter      = 1000
  )

  list(
    beta = mrash_output$beta,
    sigma2 = mrash_output$sigma2,
    pi = mrash_output$pi,
    sa2 = mrash_output$data$sa2,
    tau2 = sum(mrash_output$data$sa2 * mrash_output$pi) * mrash_output$sigma2
  )
}

# Run Mr.ASH using summary statistics (via mr.ash.rss)
#
# Computes residual summary statistics from sufficient statistics and
# calls mr.ash.rss. This enables the "ash" unmappable effects mode
# for susie_ss() and susie_rss() without requiring raw X and y.
#
# @param data Data object (must have $XtX, $Xty, $yty, $n)
# @param b_confident Vector of confident effects to subtract from residuals
# @param model Current SuSiE model
# @param params Parameters object
# @param convtol Convergence tolerance for mr.ash.rss
#
# @return List with beta, sigma2, pi, sa2, tau2
#
# @keywords internal
compute_ash_from_summary_stats <- function(data, b_confident, model, params, convtol = 1e-4) {
  # Compute residual sufficient statistics: r = y - X*b_confident
  # X'r = X'y - X'X * b_confident
  Xtr <- as.vector(data$Xty - data$XtX %*% b_confident)
  # r'r = y'y - 2*b'*X'y + b'*X'X*b
  rtr <- data$yty - 2 * sum(b_confident * data$Xty) +
         as.numeric(t(b_confident) %*% data$XtX %*% b_confident)

  XtXdiag <- diag(data$XtX)
  bhat <- Xtr / XtXdiag
  # Use n-2 df to match PVE adjustment in mr.ash.rss
  shat <- sqrt(pmax(0, (rtr - Xtr^2 / XtXdiag) / ((data$n - 2) * XtXdiag)))
  R_mat <- safe_cov2cor(data$XtX)
  var_r <- rtr / (data$n - 1)

  # Default prior grid (matching mr.ash defaults)
  if (is.null(model$ash_s0)) {
    sa2 <- (2^((0:24) / 25) - 1)^2
    model$ash_s0 <- sa2 / median(XtXdiag) * data$n
  }
  K <- length(model$ash_s0)
  if (is.null(model$ash_pi)) model$ash_pi <- rep(1 / K, K)

  fit <- mr.ash.rss(
    bhat = bhat, shat = shat, R = R_mat,
    var_y = var_r, n = data$n,
    sigma2_e = model$sigma2, s0 = model$ash_s0, w0 = model$ash_pi,
    mu1_init = model$theta,
    tol = convtol, max_iter = 1000,
    update_w0 = TRUE,
    update_sigma = params$estimate_residual_variance
  )

  list(
    beta = fit$beta,
    sigma2 = fit$sigma2,
    pi = fit$pi,
    sa2 = model$ash_s0,
    tau2 = sum(model$ash_s0 * fit$pi) * fit$sigma2
  )
}

# Run final Mr.ASH pass after SuSiE convergence
#
# After SuSiE converges, run one final Mr.ASH pass with ALL SuSiE effects
# removed (unmasked) to get the final unmappable effects estimate.
# Individual-level data uses mr.ash directly; SS data uses mr.ash.rss.
#
# @param data Data object
# @param params Parameters object
# @param model Converged SuSiE model
#
# @return Model with updated theta, tau2, ash_pi, and fitted theta values
#
# @keywords internal
run_final_ash_pass <- function(data, params, model) {
  b_susie <- colSums(model$alpha * model$mu)
  is_individual <- inherits(data, "individual")

  # Dispatch to individual (mr.ash) or SS (mr.ash.rss) backend
  if (is_individual) {
    mrash_output <- compute_ash_from_individual_data(
      data$X, data$y, b_susie, model, params
    )
  } else {
    mrash_output <- compute_ash_from_summary_stats(data, b_susie, model, params)
  }

  # Update model (shared across both paths)
  model$theta  <- mrash_output$beta
  model$tau2   <- mrash_output$tau2
  model$ash_pi <- mrash_output$pi

  # Compute fitted theta (data-representation-specific)
  if (is_individual) {
    model$X_theta <- as.vector(data$X %*% model$theta)
  } else {
    model$XtX_theta <- compute_Rv(data, model$theta)
  }

  return(model)
}

# Compute ELBO for infinitesimal effects model
#' @keywords internal
compute_elbo_inf <- function(alpha, mu, omega, lbf, sigma2, tau2, n, p,
                             eigen_vectors, eigen_values, VtXty, yty,
                             eigen_vectors_sq = NULL) {
  L <- nrow(mu)
  r <- length(eigen_values)        # rank of the stored eigenspace
  # Use cached V^2 (p x r) if available to skip the inline V^2; otherwise
  # form it.  Summation order matches the upstream reference.
  if (is.null(eigen_vectors_sq))
    eigen_vectors_sq <- eigen_vectors^2

  b <- colSums(mu * alpha)
  Vtb <- t(eigen_vectors) %*% b
  diagVtMV <- Vtb^2
  tmpD <- rep(0, p)

  for (l in seq_len(L)) {
    bl <- mu[l, ] * alpha[l, ]
    Vtbl <- t(eigen_vectors) %*% bl
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + alpha[l, ] * (mu[l, ]^2 + 1 / omega[l, ])
  }

  diagVtMV <- diagVtMV + rowSums(sweep(t(eigen_vectors_sq), 2, tmpD, `*`))

  # Compute variance
  var <- tau2 * eigen_values + sigma2

  # Compute negative ELBO.  (n - r) replaces the original (n - p): when the
  # eigenspace is stored without null-space padding (length r < p), the
  # dropped null-space components would each contribute 0.5 log(sigma^2),
  # so (n - p) double-counts them.  Equivalent for full-eigen (r = p).
  neg_elbo <- 0.5 * (n - r) * log(sigma2) + 0.5 / sigma2 * yty +
    sum(0.5 * log(var) -
          0.5 * tau2 / sigma2 * VtXty^2 / var -
          Vtb * VtXty / var +
          0.5 * eigen_values / var * diagVtMV)

  elbo <- -neg_elbo

  return(elbo)
}

# =============================================================================
# CREDIBLE SETS & POST-PROCESSING
#
# Functions for generating final output including credible sets, posterior
# inclusion probabilities, and summary statistics. These process the fitted
# model into interpretable results.
#
# Functions: n_in_CS_x, in_CS_x, n_in_CS, in_CS, get_purity
# =============================================================================

# Find how many variables in the CS.
# x is a probability vector.
#' @keywords internal
n_in_CS_x <- function(x, coverage = 0.9) {
  sum(cumsum(sort(x, decreasing = TRUE)) < coverage) + 1
}

# Return binary vector indicating if each point is in CS.
# x is a probability vector.
#' @keywords internal
in_CS_x <- function(x, coverage = 0.9) {
  n <- n_in_CS_x(x, coverage)
  o <- order(x, decreasing = TRUE)
  result <- rep(0, length(x))
  result[o[1:n]] <- 1
  return(result)
}

# Returns an l-by-p binary matrix indicating which variables are in
# susie credible sets.
#' @keywords internal
in_CS <- function(res, coverage = 0.9) {
  if (inherits(res, "susie")) {
    res <- res$alpha
  }
  return(t(apply(res, 1, function(x) in_CS_x(x, coverage))))
}

#' @keywords internal
n_in_CS <- function(res, coverage = 0.9) {
  if (inherits(res, "susie")) {
    res <- res$alpha
  }
  return(apply(res, 1, function(x) n_in_CS_x(x, coverage)))
}

# Subsample and compute min, mean, median and max abs corr.
#' @importFrom stats median
#' @keywords internal
get_purity <- function(pos, X, Xcorr, squared = FALSE, n = 100,
                       use_rfast = NULL) {
  if (is.null(use_rfast)) {
    use_rfast <- requireNamespace("Rfast", quietly = TRUE)
  }
  if (use_rfast) {
    get_upper_tri <- Rfast::upper_tri
    get_median <- Rfast::med
  } else {
    get_upper_tri <- function(R) R[upper.tri(R)]
    get_median <- median
  }
  if (length(pos) == 1) {
    return(c(1, 1, 1))
  } else {
    if (is.null(Xcorr)) {
      if (length(pos) > n) {
        pos <- sample(pos, n)
      }
      X_sub <- X[, pos]
      X_sub <- as.matrix(X_sub)
      value <- abs(get_upper_tri(safe_cor(X_sub)))
    } else {
      value <- abs(get_upper_tri(Xcorr[pos, pos]))
    }
    if (squared) {
      value <- value^2
    }
    result <- c(
      min(value),
      sum(value) / length(value),
      get_median(value)
    )
    if (any(is.na(result) | is.nan(result))) {
      stop("get_purity returned NaN/NA. Check for constant columns or data issues.")
    }
    return(result)
  }
}
