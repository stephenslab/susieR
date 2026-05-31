# =============================================================================
# CONSTRUCTOR HELPERS
# =============================================================================

#' @keywords internal
#' @noRd
resolve_model_init <- function(model_init, s_init) {
  if (!is.null(s_init)) {
    if (!is.null(model_init))
      stop("Cannot specify both 's_init' and 'model_init'.")
    warning_message("s_init is deprecated and will be removed in a future ",
                    "version of susieR. Please use model_init instead.",
                    style = "hint")
    model_init <- s_init
  }
  model_init
}

#' @keywords internal
#' @noRd
normalize_null_weight <- function(null_weight, prior_weights, p) {
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric.")
    if (length(null_weight) != 1L)
      stop("Null weight must be numeric.")
    if (is.na(null_weight))
      stop("Null weight must be numeric.")
    if (null_weight == 0)
      null_weight <- NULL
  }

  if (!is.null(null_weight)) {
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1.")

    if (is.null(prior_weights)) {
      prior_weights <- c(rep((1 - null_weight) / p, p), null_weight)
    } else {
      prior_weights <- c(prior_weights * (1 - null_weight), null_weight)
    }
  }

  list(null_weight = null_weight,
       prior_weights = prior_weights,
       add_null = !is.null(null_weight))
}

#' @keywords internal
#' @noRd
normalize_prior_weights <- function(prior_weights, p) {
  if (is.null(prior_weights))
    prior_weights <- rep(1 / p, p)
  if (length(prior_weights) != p)
    stop("Prior weights must have length p.")
  if (all(prior_weights == 0))
    stop("Prior weight should be greater than 0 for at least one variable.")
  prior_weights / sum(prior_weights)
}

#' @keywords internal
#' @noRd
normalize_summary_stats_input <- function(z = NULL, bhat = NULL,
                                          shat = NULL) {
  if (!is.null(z) && (!is.null(bhat) || !is.null(shat))) {
    stop("Please provide either z or (bhat, shat), but not both.")
  }
  if (is.null(z) && (is.null(bhat) || is.null(shat)))
    stop("Please provide either z or (bhat, shat).")

  if (is.null(z)) {
    if (!is.numeric(bhat) || !is.numeric(shat))
      stop("bhat and shat must be numeric.")
    variable_names <- names(bhat)
    bhat <- as.numeric(bhat)
    if (length(shat) == 1L)
      shat <- rep(shat, length(bhat))
    else
      shat <- as.numeric(shat)
    if (length(bhat) != length(shat))
      stop("The lengths of bhat and shat do not agree.")
    if (anyNA(bhat) || anyNA(shat))
      stop("bhat, shat cannot have missing values.")
    if (any(shat <= 0))
      stop("shat cannot have zero or negative elements.")
    z <- bhat / shat
    if (!is.null(variable_names))
      names(z) <- variable_names
  } else {
    if (!is.numeric(z))
      stop("z must be numeric.")
    variable_names <- names(z)
    z <- as.numeric(z)
    if (!is.null(variable_names))
      names(z) <- variable_names
  }

  if (length(z) < 1L)
    stop("Input vector z should have at least one element.")

  list(z = z, bhat = bhat, shat = shat, variable_names = variable_names)
}

#' @keywords internal
#' @noRd
apply_pve_adjustment <- function(z, n = NULL) {
  if (is.null(z) || is.null(n) || n <= 1) {
    return(list(z = z, adjustment = NULL, adjusted = FALSE))
  }
  adj <- (n - 1) / (z^2 + n - 2)
  list(z = sqrt(adj) * z, adjustment = adj, adjusted = TRUE)
}

#' @keywords internal
#' @noRd
summary_stats_working_quantities <- function(z, n = NULL, shat = NULL,
                                             var_y = NULL,
                                             pve_adjustment = NULL,
                                             prior_variance = 50,
                                             scaled_prior_variance = 0.2) {
  p <- length(z)

  if (is.null(n)) {
    return(list(path = "z_only", n = 2, Xty = z, yty = 1,
                scaled_prior_variance = prior_variance,
                ser_betahat = z, ser_shat2 = rep(1, p), ser_sigma2 = 1,
                ser_V_init = prior_variance,
                ser_scale_factors = rep(1, p)))
  }

  nm1 <- n - 1
  yty <- nm1 * if (!is.null(var_y)) var_y else 1
  y_var <- yty / nm1

  if (!is.null(shat) && !is.null(var_y)) {
    if (is.null(pve_adjustment))
      stop("PVE adjustment is required with shat and var_y.")
    XtXdiag <- var_y * pve_adjustment / (shat^2)
    Xty <- z * sqrt(pve_adjustment) * var_y / shat
    scale_factors <- sqrt(XtXdiag / nm1)
    ser_betahat <- (1 / nm1) * (Xty / scale_factors)
    return(list(path = "original_scale", n = n, Xty = z *
                  sqrt(pve_adjustment) * var_y / shat,
                yty = yty, scaled_prior_variance = scaled_prior_variance,
                XtXdiag = XtXdiag,
                ser_betahat = ser_betahat,
                ser_shat2 = rep(y_var / nm1, p), ser_sigma2 = y_var,
                ser_V_init = scaled_prior_variance * y_var,
                ser_scale_factors = scale_factors))
  }

  Xty <- sqrt(nm1) * z
  list(path = "standardized_scale", n = n, Xty = Xty,
       yty = yty, scaled_prior_variance = scaled_prior_variance,
       XtX_multiplier = nm1, X_multiplier = sqrt(nm1),
       ser_betahat = (1 / nm1) * Xty,
       ser_shat2 = rep(y_var / nm1, p),
       ser_sigma2 = y_var, ser_V_init = scaled_prior_variance * y_var,
       ser_scale_factors = rep(1, p))
}

# =============================================================================
# INDIVIDUAL-LEVEL DATA CONSTRUCTOR
#
# Constructs data and params objects for SuSiE from individual-level data (X, y).
# Handles data preprocessing, parameter validation, and object creation.
# =============================================================================
#'
#' @return A list containing:
#' \item{data}{A processed list containing X and y matrices with appropriate scaling
#' attributes and sample dimensions}
#' \item{params}{Validated params object with all input algorithm parameters}
#'
#' @keywords internal
#' @importFrom stats var
#' @noRd
individual_data_constructor <- function(X, y, L = min(10, ncol(X)),
                                        scaled_prior_variance = 0.2,
                                        residual_variance = NULL,
                                        prior_weights = NULL,
                                        null_weight = 0,
                                        standardize = TRUE,
                                        intercept = TRUE,
                                        estimate_residual_variance = TRUE,
                                        estimate_residual_method = "MoM",
                                        estimate_prior_variance = TRUE,
                                        estimate_prior_method = "optim",
                                        prior_variance_grid = NULL,
                                        mixture_weights = NULL,
                                        unmappable_effects = "none",
                                        check_null_threshold = 0,
                                        prior_tol = 1e-9,
                                        residual_variance_upperbound = Inf,
                                        model_init = NULL,
                                        s_init = NULL,
                                        coverage = 0.95,
                                        min_abs_corr = 0.5,
                                        compute_univariate_zscore = FALSE,
                                        na.rm = FALSE,
                                        max_iter = 100,
                                        tol = NULL,
                                        convergence_method = "elbo",
                                        verbose = FALSE,
                                        track_fit = FALSE,
                                        residual_variance_lowerbound = NULL,
                                        refine = FALSE,
                                        n_purity = 100,
                                        alpha0 = NULL,
                                        beta0 = NULL,
                                        slot_prior = NULL,
                                        L_greedy = NULL,
                                        greedy_lbf_cutoff = 0.1) {

  model_init <- resolve_model_init(model_init, s_init)

  # Validate input X
  if (!(is.double(X) & is.matrix(X)) &
      !inherits(X, "sparseMatrix") &
      is.null(attr(X, "matrix.type"))) {
    stop("Input X must be a double-precision matrix, or a sparse matrix, or ",
         "a trend filtering matrix.")
  }

  if (anyNA(X)) {
    stop("X contains NA values.")
  }

  # When n >> p, holding X in memory at every IBSS iteration is wasteful;
  # nudge the user toward the sufficient-statistics path. The hint is
  # advisory only and does not gate on `verbose` (matches the style of the
  # constructor's other hints, e.g. the susie_rss max_iter default).
  if (is.matrix(X) && nrow(X) >= 2 * ncol(X)) {
    warning_message(
      "nrow(X) = ", nrow(X), " >= 2 * ncol(X) = ", 2 * ncol(X), ". ",
      "Consider precomputing sufficient statistics with compute_suff_stat() ",
      "and fitting with susie_ss() instead -- this avoids holding X in ",
      "memory at every iteration and lets you reuse XtX across multiple y.",
      style = "hint"
    )
  }

  # Constant column check for regular matrix
  if (is.null(attr(X, "matrix.type")) || attr(X, "matrix.type") != "tfmatrix") {
    col_vars <- apply(X, 2, var)
    const_cols <- which(col_vars == 0 | is.na(col_vars))
    if (length(const_cols) > 0) {
      warning_message(sprintf("X contains %d constant columns (first few cols: %s).",
                 length(const_cols), paste(head(const_cols, 10), collapse = ", ")))
    }
  }

  # Handle missing values in y
  if (anyNA(y)) {
    if (na.rm) {
      samples_kept <- which(!is.na(y))
      y <- y[samples_kept]
      X <- X[samples_kept, , drop = FALSE]
    } else {
      stop("Input y must not contain missing values.")
    }
  }

  # Set residual_variance_lowerbound
  if (is.null(residual_variance_lowerbound)) {
    residual_variance_lowerbound <- var(drop(y)) / 1e4
  }

  mean_y <- mean(y)

  # Force required preprocessing for unmappable effects methods
  if (unmappable_effects != "none") {
    if (!intercept) {
      warning_message("Unmappable effects methods require centered data. ",
                      "Setting intercept=TRUE.", style = "hint")
      intercept <- TRUE
    }
    if (!standardize) {
      warning_message("Unmappable effects methods require scaled data. ",
                      "Setting standardize=TRUE.", style = "hint")
      standardize <- TRUE
    }
  }

  # Check for incompatible parameter combination
  if (unmappable_effects != "none" &&
      estimate_residual_method == "NIG") {
    stop("The combination of unmappable_effects = '", unmappable_effects,
         "' with estimate_residual_method = 'NIG' is not supported. ",
         "Please use estimate_residual_method = 'MoM' or 'MLE' instead.")
  }

  # Check for incompatible parameter combination
  if (unmappable_effects %in% c("ash", "ash_filter_archived") && estimate_prior_method == "EM") {
    stop("The combination of unmappable_effects = 'ash' with ",
         "estimate_prior_method = 'EM' is not supported. ",
         "Please use estimate_prior_method = 'optim' instead.")
  }

  nw <- normalize_null_weight(null_weight, prior_weights, ncol(X))
  null_weight <- nw$null_weight
  prior_weights <- nw$prior_weights

  if (nw$add_null) {
    # add the extra 0 column to X
    X <- cbind(X, 0)
  }

  # Store dimensions
  n <- nrow(X)
  p <- ncol(X)

  prior_weights <- normalize_prior_weights(prior_weights, p)

  # nocov start
  if (p > 1000 & !requireNamespace("Rfast", quietly = TRUE)) {
    warning_message("For an X with many columns, please consider installing ",
                    "the Rfast package for more efficient credible set (CS) ",
                    "calculations.",
                    style = "hint")
  }
  # nocov end

  # Center y if intercept is included
  if (intercept) {
    y <- y - mean_y
  }

  # Compute and set X matrix attributes
  out <- compute_colstats(X, center = intercept, scale = standardize)
  attr(X, "scaled:center") <- out$cm
  attr(X, "scaled:scale") <- out$csd
  attr(X, "d") <- out$d

  # Create params object with all input parameters
  params_object <- list(
    L = L,
    scaled_prior_variance = scaled_prior_variance,
    residual_variance = residual_variance,
    prior_weights = prior_weights,
    null_weight = null_weight,
    estimate_residual_variance = estimate_residual_variance,
    estimate_residual_method = estimate_residual_method,
    estimate_prior_variance = estimate_prior_variance,
    estimate_prior_method = estimate_prior_method,
    prior_variance_grid = prior_variance_grid,
    mixture_weights = mixture_weights,
    unmappable_effects = unmappable_effects,
    check_null_threshold = check_null_threshold,
    prior_tol = prior_tol,
    residual_variance_upperbound = residual_variance_upperbound,
    model_init = model_init,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    compute_univariate_zscore = compute_univariate_zscore,
    max_iter = max_iter,
    tol = tol,
    convergence_method = convergence_method,
    verbose = verbose,
    track_fit = track_fit,
    residual_variance_lowerbound = residual_variance_lowerbound,
    refine = refine,
    n_purity = n_purity,
    alpha0 = alpha0,
    beta0 = beta0,
    n = n,
    use_NIG = FALSE,  # Will be set by validation function
    intercept = intercept,
    standardize = standardize,
    slot_prior = slot_prior,
    L_greedy = L_greedy,
    greedy_lbf_cutoff = greedy_lbf_cutoff
  )

  # Validate and apply parameter overrides
  params_object <- validate_and_override_params(params_object)
  data_object <- structure(
    list(
      X = X,
      y = y,
      mean_y = mean_y,
      n = n,
      p = p
    ),
    class = "individual"
  )

  # Configure data object based on params
  data_object <- configure_data(data_object, params_object)

  return(list(data = data_object, params = params_object))
}

# =============================================================================
# SUFFICIENT STATISTICS DATA CONSTRUCTOR
#
# Constructs data and params objects for SuSiE from sufficient statistics (XtX, Xty, yty).
# Handles data preprocessing, parameter validation, and object creation.
# =============================================================================
#'
#' @return A list containing:
#' \item{data}{A processed list containing XtX, Xty, yty matrices with appropriate scaling
#' attributes and sample dimensions}
#' \item{params}{Validated params object with all input algorithm parameters}
#'
#' @keywords internal
#' @noRd
sufficient_stats_constructor <- function(Xty, yty, n,
                                         XtX = NULL, X = NULL,
                                         L = min(10, if (!is.null(XtX)) ncol(XtX) else ncol(X)),
                                         X_colmeans = NA, y_mean = NA,
                                         maf = NULL, maf_thresh = 0,
                                         check_input = FALSE,
                                         r_tol = 1e-8,
                                         standardize = TRUE,
                                         scaled_prior_variance = 0.2,
                                         residual_variance = NULL,
                                         prior_weights = NULL,
                                         null_weight = 0,
                                         model_init = NULL,
                                         s_init = NULL,
                                         estimate_residual_variance = TRUE,
                                         estimate_residual_method = "MoM",
                                         residual_variance_lowerbound = 0,
                                         residual_variance_upperbound = Inf,
                                         estimate_prior_variance = TRUE,
                                         estimate_prior_method = "optim",
                                         prior_variance_grid = NULL,
                                         mixture_weights = NULL,
                                         unmappable_effects = "none",
                                         check_null_threshold = 0,
                                         prior_tol = 1e-9,
                                         max_iter = 100,
                                         tol = NULL,
                                         convergence_method = "elbo",
                                         coverage = 0.95,
                                         min_abs_corr = 0.5,
                                         n_purity = 100,
                                         verbose = FALSE,
                                         track_fit = FALSE,
                                         check_prior = FALSE,
                                         refine = FALSE,
                                         alpha0 = NULL,
                                         beta0 = NULL,
                                         slot_prior = NULL,
                                         L_greedy = NULL,
                                         greedy_lbf_cutoff = 0.1) {

  model_init <- resolve_model_init(model_init, s_init)

  # Validate required inputs
  if (missing(n)) {
    stop("n must be provided.")
  }
  if (n <= 1) {
    stop("n must be greater than 1.")
  }

  if (is.null(X)) {
    # XtX path: validate XtX
    if (is.null(XtX) || missing(Xty) || missing(yty)) {
      stop("XtX, Xty, yty must all be provided.")
    }

    if (!(is.double(XtX) && is.matrix(XtX)) &&
        !inherits(XtX, "sparseMatrix")) {
      stop("XtX must be a numeric dense or sparse matrix.")
    }

    if (ncol(XtX) != length(Xty)) {
      stop(paste0(
        "The dimension of XtX (", nrow(XtX), " by ", ncol(XtX),
        ") does not agree with expected (", length(Xty), " by ",
        length(Xty), ")."
      ))
    }

    # nocov start
    if (ncol(XtX) > 1000 & !requireNamespace("Rfast", quietly = TRUE)) {
      warning_message("For large R or large XtX, consider installing the ",
                      "Rfast package for better performance.",
                      style = "hint")
    }
    # nocov end

    # Ensure XtX is symmetric
    if (!is_symmetric_matrix(XtX)) {
      XtX <- symmetrize_warned(XtX, "XtX")
    }

    # Apply MAF filter if provided
    if (!is.null(maf)) {
      if (length(maf) != length(Xty)) {
        stop(paste("The length of maf does not agree with expected", length(Xty), "."))
      }
      id <- which(maf > maf_thresh)
      XtX <- XtX[id, id]
      Xty <- Xty[id]
      if (!is.null(prior_weights))
        prior_weights <- prior_weights[id]
    }

    # Additional validation
    if (anyNA(XtX)) {
      stop("Input XtX matrix contains NAs.")
    }

    # Positive-semidefinite check
    if (check_input) {
      semi_pd <- check_semi_pd(XtX, r_tol)
      if (!semi_pd$status) {
        stop("XtX is not a positive semidefinite matrix.")
      }

      # Check whether Xty lies in space spanned by non-zero eigenvectors of XtX
      proj <- check_projection(semi_pd$matrix, Xty)
      if (!proj$status) {
        warning_message("Xty does not lie in the space of the non-zero eigenvectors ",
                        "of XtX.")
      }
    }
  } else {
    # X low-rank path: validate X
    if (ncol(X) != length(Xty)) {
      stop(paste0(
        "The number of columns of X (", ncol(X),
        ") does not agree with the length of Xty (", length(Xty), ")."
      ))
    }
  }

  # Common validation for Xty
  if (any(is.infinite(Xty))) {
    stop("Input Xty contains infinite values.")
  }
  if (anyNA(Xty)) {
    Xty <- replace_na_zero_warned(Xty, "Xty")
  }

  # Define p before null_weight handling
  p <- if (!is.null(XtX)) ncol(XtX) else ncol(X)

  nw <- normalize_null_weight(null_weight, prior_weights, p)
  null_weight <- nw$null_weight
  prior_weights <- nw$prior_weights

  if (nw$add_null) {
    if (!is.null(XtX)) {
      XtX <- cbind(rbind(XtX, 0), 0)
    }
    if (!is.null(X)) {
      X <- cbind(X, 0)
    }
    Xty <- c(Xty, 0)
    if (length(X_colmeans) == 1) {
      X_colmeans <- rep(X_colmeans, p)
    }
    if (length(X_colmeans) != p) {
      stop("The length of X_colmeans does not agree with number of variables.")
    }
    # Add 0 for null column
    X_colmeans <- c(X_colmeans, 0)
    # Update p after adding null column
    p <- p + 1
  }

  prior_weights <- normalize_prior_weights(prior_weights, p)

  # Standardize if requested
  if (!is.null(X)) {
    # Low-rank X path: standardize columns of X
    if (standardize) {
      dXtX <- colSums(X^2)
      csd <- sqrt(dXtX / (n - 1))
      csd[csd == 0] <- 1
      X <- t(t(X) / csd)
      Xty <- Xty / csd
    } else {
      csd <- rep(1, length = p)
    }
    attr(X, "d") <- colSums(X^2)
    attr(X, "scaled:scale") <- csd
    colnames(X) <- names(Xty)
  } else {
    # XtX path: standardize XtX
    if (standardize) {
      dXtX <- diag(XtX)
      csd <- sqrt(dXtX / (n - 1))
      csd[csd == 0] <- 1
      XtX <- t((1 / csd) * XtX) / csd
      Xty <- Xty / csd
    } else {
      csd <- rep(1, length = p)
    }
    attr(XtX, "d") <- diag(XtX)
    attr(XtX, "scaled:scale") <- csd
  }

  if (length(X_colmeans) == 1) {
    X_colmeans <- rep(X_colmeans, p)
  }
  if (length(X_colmeans) != p) {
    stop(
      "`X_colmeans` length (", length(X_colmeans),
      ") does not match number of variables (", p, ")."
    )
  }

  # Create params object with all input parameters
  params_object <- list(
    L = L,
    scaled_prior_variance = scaled_prior_variance,
    residual_variance = residual_variance,
    prior_weights = prior_weights,
    null_weight = null_weight,
    estimate_residual_variance = estimate_residual_variance,
    estimate_residual_method = estimate_residual_method,
    estimate_prior_variance = estimate_prior_variance,
    estimate_prior_method = estimate_prior_method,
    prior_variance_grid = prior_variance_grid,
    mixture_weights = mixture_weights,
    unmappable_effects = unmappable_effects,
    check_null_threshold = check_null_threshold,
    prior_tol = prior_tol,
    residual_variance_upperbound = residual_variance_upperbound,
    model_init = model_init,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    compute_univariate_zscore = FALSE,  # SS doesn't support univariate zscore
    max_iter = max_iter,
    tol = tol,
    convergence_method = convergence_method,
    verbose = verbose,
    track_fit = track_fit,
    residual_variance_lowerbound = residual_variance_lowerbound,
    refine = refine,
    n_purity = n_purity,
    alpha0 = alpha0,
    beta0 = beta0,
    n = n,
    use_NIG = FALSE,
    intercept = FALSE,  # SS always uses intercept = FALSE
    standardize = standardize,
    check_prior = check_prior,
    slot_prior = slot_prior,
    L_greedy = L_greedy,
    greedy_lbf_cutoff = greedy_lbf_cutoff
  )

  # Validate and apply parameter overrides
  params_object <- validate_and_override_params(params_object)

  # Assemble data object
  data_object <- structure(
    list(
      XtX = XtX,
      X = X,
      Xty = Xty,
      yty = yty,
      n = n,
      p = p,
      X_colmeans = X_colmeans,
      y_mean = y_mean
    ),
    class = "ss"
  )

  # Configure data object based on params
  data_object <- configure_data(data_object, params_object)

  return(list(data = data_object, params = params_object))
}

# =============================================================================
# SUMMARY STATISTICS (RSS) DATA CONSTRUCTOR
#
# Constructs data and params objects for SuSiE from summary statistics
# (z-scores, R matrix, or multi-panel R/X reference).
# =============================================================================
#'
#' @return A list containing:
#' \item{data}{A processed list containing converted matrices with appropriate scaling
#' attributes and sample dimensions}
#' \item{params}{Validated params object with all input algorithm parameters}
#'
#' @keywords internal
#' @noRd
summary_stats_constructor <- function(z = NULL, R = NULL, X = NULL,
                                      n = NULL, bhat = NULL,
                                      shat = NULL, var_y = NULL,
                                      L = min(10, if (!is.null(R)) ncol(R) else ncol(X)),
                                      maf = NULL,
                                      maf_thresh = 0,
                                      prior_variance = 50,
                                      scaled_prior_variance = 0.2,
                                      residual_variance = NULL,
                                      prior_weights = NULL,
                                      null_weight = 0,
                                      standardize = TRUE,
                                      estimate_residual_variance = FALSE,
                                      estimate_residual_method = "MoM",
                                      estimate_prior_variance = TRUE,
                                      estimate_prior_method = "optim",
                                      prior_variance_grid = NULL,
                                      mixture_weights = NULL,
                                      unmappable_effects = "none",
                                      check_null_threshold = 0,
                                      prior_tol = 1e-9,
                                      residual_variance_lowerbound = 0,
                                      residual_variance_upperbound = Inf,
                                      model_init = NULL,
                                      s_init = NULL,
                                      coverage = 0.95,
                                      min_abs_corr = 0.5,
                                      max_iter = 100,
                                      tol = NULL,
                                      convergence_method = "elbo",
                                      verbose = FALSE,
                                      track_fit = FALSE,
                                      check_input = FALSE,
                                      check_prior = TRUE,
                                      n_purity = 100,
                                      r_tol = 1e-8,
                                      refine = FALSE,
                                      R_finite = NULL,
                                      R_mismatch = "none",
                                      R_mismatch_method = "mle",
                                      eig_delta_rel = 1e-3,
                                      eig_delta_abs = 0,
                                      artifact_threshold = 0.1,
                                      R_sensitivity_threshold = log(20),
                                      alpha0 = NULL,
                                      beta0 = NULL,
                                      slot_prior = NULL,
                                      L_greedy = NULL,
                                      greedy_lbf_cutoff = 0.1) {

  # Validate: exactly one of R or X must be provided
  if (is.null(R) && is.null(X))
    stop("Please provide either R (correlation matrix) or X (factor matrix).")
  if (!is.null(R) && !is.null(X))
    stop("Please provide either R or X, but not both.")

  # Default max_iter for susie_rss-style summary-statistics fitting:
  # callers (such as susie_rss) that pass max_iter = NULL get the SuSiE-RSS
  # default of 50 plus a one-time hint. Direct constructor callers keep
  # the constructor's own `max_iter = 100` default and never see the hint.
  if (is.null(max_iter)) {
    max_iter <- 50
    warning_message("Setting max_iter = 50 for the SuSiE RSS model because ",
                    "slow convergence is often a sign of unstable summary-statistics ",
                    "fitting. To disable this message, explicitly set max_iter = 50 ",
                    "or another value in the susie_rss() call.",
                    style = "hint")
  }

  model_init <- resolve_model_init(model_init, s_init)

  # NIG prior requires an explicit sample size n: the default alpha0/beta0
  # scale as 1/sqrt(n) and the NIG marginal likelihood depends on n. Without
  # n, summary_stats_constructor bumps n to 2 internally (see below), which
  # would silently corrupt the NIG posterior. Reject early with a clear error.
  if (estimate_residual_method == "NIG" &&
      (is.null(n) || !is.numeric(n) || length(n) != 1 ||
       !is.finite(n) || n < 1)) {
    stop("estimate_residual_method = \"NIG\" requires a valid sample ",
         "size `n` (got n = ", paste(n, collapse = ""), "). ",
         "For susie_rss(), pass `n` explicitly.")
  }

  # PVE-adjusted z-scores: shrink large z toward zero to account for
  # winner's curse. Applied to ALL paths when sample size is available.
  # Guard: z may be NULL when bhat/shat are provided (converted later).
  # Snapshot the user-supplied z before mutation so the multi-panel
  # sub-fit dispatch (below) can hand the *unadjusted* z to per-panel
  # constructor calls — they will apply the PVE adjustment themselves
  # and double-applying it would silently corrupt Xty.
  z_user <- z
  pve <- apply_pve_adjustment(z, n)
  z <- pve$z
  adj <- pve$adjustment
  pve_adjusted <- pve$adjusted

  is_multipanel <- (is.list(X) && !is.matrix(X)) ||
                   (is.list(R) && !is.matrix(R))
  R_mismatch <- match.arg(R_mismatch, c("none", "eb", "eb_ser_init", "eb_force_init", "eb_no_init"))
  R_mismatch_method <- match.arg(R_mismatch_method, c("mle", "map"))
  if (!is.numeric(R_sensitivity_threshold) ||
      length(R_sensitivity_threshold) != 1L ||
      !is.finite(R_sensitivity_threshold) ||
      R_sensitivity_threshold < 0)
    stop("R_sensitivity_threshold must be a single nonnegative finite numeric.")
  if (!is.numeric(eig_delta_rel) || length(eig_delta_rel) != 1L ||
      eig_delta_rel < 0)
    stop("eig_delta_rel must be a single nonnegative numeric.")
  if (!is.numeric(eig_delta_abs) || length(eig_delta_abs) != 1L ||
      eig_delta_abs < 0)
    stop("eig_delta_abs must be a single nonnegative numeric.")
  if (!is.numeric(artifact_threshold) || length(artifact_threshold) != 1L ||
      artifact_threshold < 0 || artifact_threshold > 1)
    stop("artifact_threshold must be a single numeric in [0, 1].")
  R_finite_explicit_false <- identical(R_finite, FALSE)
  if (isTRUE(R_finite) && is.null(X))
    stop("R_finite = TRUE requires X input. When using precomputed R, ",
         "provide the reference sample size explicitly.")
  R_finite <- resolve_R_finite(R_finite, if (!is.null(X)) X else R,
                               is_multipanel)

  # sigma^2 and lambda_bias both inflate the residual variance and are
  # only weakly jointly identified, so sigma^2 is fixed when R_mismatch is
  # active.
    if (R_mismatch != "none" && isTRUE(estimate_residual_variance)) {
      warning_message(
        "Joint estimation of sigma^2 and lambda_bias is weakly identified; ",
        "sigma^2 and lambda_bias can trade off, especially early in EM. ",
        "Consider track_fit=TRUE to monitor, or R_mismatch_method='map' ",
        "to penalize lambda_bias.",
        style = "hint"
      )
    }
                               
  if (is_multipanel) {
    if (!is.null(bhat) || !is.null(shat)) {
      stop("Parameters 'bhat' and 'shat' are not supported in the ",
           "multi-panel summary-statistics path. ",
           "Please provide z-scores directly.")
    }
    if (!is.null(var_y))
      stop("Parameter 'var_y' is not supported in the multi-panel path.")

    # n required for multi-panel sub-fits and ss_mixture_constructor.
    if (is.null(n))
      stop("Sample size 'n' is required for multi-panel mode.")

    # Validate panel-list elements (lighter than ss_mixture's per-p check;
    # this catches type errors before the per-panel sub-fits run, so that
    # the user sees a list-element error instead of a confusing per-panel
    # constructor error).
    if (!is.null(R)) {
      for (k in seq_along(R)) {
        if (!is.matrix(R[[k]]) || !is.numeric(R[[k]]))
          stop("Each element of R list must be a numeric matrix.")
        if (nrow(R[[k]]) != ncol(R[[k]]))
          stop("Each element of R list must be square.")
      }
    } else {
      for (k in seq_along(X)) {
        if (!is.matrix(X[[k]]) || !is.numeric(X[[k]]))
          stop("Each element of X list must be a numeric matrix.")
      }
      # Center each panel before sub-fits so that crossprod(Xk) gives a
      # covariance-like quantity, matching ss_mixture_constructor's
      # downstream expectation.
      X <- lapply(X, function(Xk) {
        cm <- colMeans(Xk)
        if (max(abs(cm)) > 1e-10 * max(abs(Xk)))
          Xk <- t(t(Xk) - cm)
        Xk
      })
    }

    # Capture parent args BEFORE force_pip so sub-fits inherit the user's
    # original convergence_method (today's match.call-based recursion does
    # the same via the captured user call).
    parent_args <- mget(names(formals()))
    parent_args$z <- z_user                     # restore pre-PVE z (sub-fits PVE-adjust on their own)
    parent_args$max_iter <- max_iter
    parent_args$model_init <- model_init
    parent_args$R_finite <- R_finite
    parent_args$X <- X
    parent_args$estimate_residual_variance <- estimate_residual_variance

    if (convergence_method[1] == "elbo") {
      convergence_method <- force_pip(
        "multi-panel mixture weight updates change R(omega) each ",
        "iteration, which prevents ELBO monotonicity")
    }

    panel_arg <- if (!is.null(R)) "R" else "X"
    panels <- if (panel_arg == "R") R else X
    sub <- pick_init_panel_via_subfits(panels, panel_arg, parent_args)

    result <- ss_mixture_constructor(
      z = z, R = R, X = X, n = n, L = L,
      maf = maf, maf_thresh = maf_thresh,
      scaled_prior_variance = scaled_prior_variance,
      residual_variance = residual_variance,
      prior_weights = prior_weights, null_weight = null_weight,
      standardize = standardize,
      estimate_residual_variance = estimate_residual_variance,
      estimate_residual_method = estimate_residual_method,
      estimate_prior_variance = estimate_prior_variance,
      estimate_prior_method = estimate_prior_method,
      prior_variance_grid = prior_variance_grid,
      mixture_weights = mixture_weights,
      unmappable_effects = unmappable_effects,
      check_null_threshold = check_null_threshold, prior_tol = prior_tol,
      residual_variance_lowerbound = residual_variance_lowerbound,
      residual_variance_upperbound = residual_variance_upperbound,
      model_init = model_init, coverage = coverage,
      min_abs_corr = min_abs_corr, max_iter = max_iter, tol = tol,
      convergence_method = convergence_method, verbose = verbose,
      track_fit = track_fit, check_input = check_input,
      check_prior = check_prior, n_purity = n_purity,
      r_tol = r_tol, refine = refine, R_finite = R_finite,
      R_mismatch = R_mismatch, R_mismatch_method = R_mismatch_method,
      eig_delta_rel = eig_delta_rel,
      eig_delta_abs = eig_delta_abs, artifact_threshold = artifact_threshold,
      R_sensitivity_threshold = R_sensitivity_threshold,
      alpha0 = alpha0, beta0 = beta0, slot_prior = slot_prior,
      L_greedy = L_greedy, greedy_lbf_cutoff = greedy_lbf_cutoff,
      init_panel = sub$idx
    )
    # Stash sub-fit metadata so susie_rss can attach it to the workhorse
    # output (single_panel_fits) and emit the verbose ELBO summary.
    result$multi_panel_meta <- list(fits = sub$fits, elbos = sub$elbos)
    return(result)
  }

  # Single-panel X handling: validate, center, and convert to R when the
  # full-rank or var_y/shat path requires the correlation matrix.
  if (!is.null(X)) {
    if (!is.matrix(X) || !is.numeric(X))
      stop("X must be a numeric matrix.")
    cm <- colMeans(X)
    if (max(abs(cm)) > 1e-10 * max(abs(X)))
      X <- t(t(X) - cm)
    needs_R <- !is.null(var_y) && !is.null(shat)
    if (needs_R && nrow(X) < ncol(X)) {
      warning_message(
        "X is provided as a low-rank factor matrix, but var_y/shat ",
        "requires the full correlation matrix R. Forming ",
        "R = cov2cor(crossprod(X)/nrow(X)) and using the standard path.",
        style = "hint")
    }
    if (nrow(X) >= ncol(X) || needs_R) {
      R <- safe_cor(X)
      X <- NULL
    }
  }

  # Auto-switch to PIP convergence when the SER likelihood is dynamically
  # inflated by finite-reference R uncertainty or R-mismatch correction.
  if ((!is.null(R_finite) || R_mismatch != "none") &&
      convergence_method[1] == "elbo") {
    convergence_method <- force_pip(
      "R uncertainty inflation modifies per-variant SER ",
      "likelihoods, which prevents a consistent model-level ELBO")
  }

  # Issue warning for estimate_residual_variance if TRUE, unless the caller
  # explicitly sets R_finite = FALSE to confirm that finite-reference
  # correction is not needed.
  if (estimate_residual_variance && !R_finite_explicit_false) {
    warning_message("estimate_residual_variance = TRUE is not recommended ",
            "unless R is the \"in-sample\" R matrix; that is, the ",
            "correlation matrix obtained using the exact same data ",
            "matrix X that was used for the other summary statistics. ",
            "If R is in-sample, set R_finite = FALSE to silence this ",
            "warning. When covariates are included in the univariate ",
            "regressions that produced the summary statistics, also ",
            "consider removing these effects from X before computing R.")
  }


  # For SuSiE-ash with summary statistics, recommend providing bhat/shat/var_y
  # for best agreement with individual-level analysis. The z+R-only path
  # operates on a standardized scale (var_y=1) and may give different results.
  if (unmappable_effects %in% c("ash", "ash_filter_archived") && is.null(bhat) && is.null(var_y)) {
    warning_message("SuSiE-ash with z-scores and R only operates on a ",
            "standardized scale. For best agreement with ",
            "individual-level analysis, provide bhat, shat, and ",
            "var_y instead of z-scores.", style = "hint")
  }

  summary_input <- normalize_summary_stats_input(z = z, bhat = bhat,
                                                 shat = shat)
  z <- summary_input$z
  bhat <- summary_input$bhat
  shat <- summary_input$shat
  p <- length(z)

  # Check dimensions of R or X
  if (!is.null(R)) {
    if (nrow(R) != p) {
      stop(paste0(
        "The dimension of R (", nrow(R), " x ", ncol(R), ") does not ",
        "agree with expected (", p, " x ", p, ")."
      ))
    }
  } else if (!is.null(X)) {
    if (ncol(X) != p) {
      stop(paste0(
        "The number of columns of X (", ncol(X), ") does not ",
        "agree with expected (", p, ")."
      ))
    }
  }

  # Check input n
  if (!is.null(n)) {
    if (n <= 1) {
      stop("n must be greater than 1.")
    }
  }

  z[is.na(z)] <- 0

  # Apply PVE adjustment if not already done (when z was computed from bhat/shat)
  if (!pve_adjusted && !is.null(n) && n > 1) {
    pve <- apply_pve_adjustment(z, n)
    z <- pve$z
    adj <- pve$adjustment
    pve_adjusted <- pve$adjusted
  }

  # MAF filter (after z-scores are computed)
  if (!is.null(maf)) {
    if (length(maf) != length(z)) {
      stop(paste0("The length of maf does not agree with expected ", length(z)))
    }
    id <- which(maf > maf_thresh)
    if (!is.null(R)) R <- R[id, id]
    if (!is.null(X)) X <- X[, id, drop = FALSE]
    z <- z[id]
    if (!is.null(bhat)) bhat <- bhat[id]
    if (!is.null(shat)) shat <- shat[id]
    if (!is.null(adj)) adj <- adj[id]
    if (!is.null(prior_weights))
      prior_weights <- prior_weights[id]
    # Update p after filtering
    p <- length(z)
  }

  # Standardize X so X'X = R (correlation matrix). The model assumes
  # column-standardized X; without this, X'X/B gives sample covariance
  # which != correlation when columns have different variances.
  if (!is.null(X)) {
    X <- standardize_X(X)
  }

  R_mismatch <- match.arg(R_mismatch, c("none", "eb", "eb_ser_init", "eb_force_init", "eb_no_init"))
  R_finite_B <- R_finite
  if (R_mismatch != "none" && is.null(R_finite_B))
    R_finite_B <- Inf

  # R diagnostics (static, computed once at initialization).
  # X is standardized (X'X = R) at this point.
  R_finite_diagnostics <- NULL
  if (!is.null(R_finite_B)) {
    R_finite_diagnostics <- compute_R_finite_diagnostics(
      X = X, R = R, B = R_finite_B, p = length(z),
      x_is_standardized = TRUE)
  }

  # Cache eigen(R) for the Q_art QC diagnostic whenever R_mismatch is active.
  # Works for both R-input and X-input: after standardize_X, crossprod(X) == R.
  # Reuses the attr(R, "eigen") convention when the caller pre-computed it.
  eigen_R_cache <- NULL
  if (R_mismatch != "none") {
    eigen_R_cache <- if (!is.null(R)) attr(R, "eigen") else NULL
    if (is.null(eigen_R_cache)) {
      R_for_eigen <- if (!is.null(R)) R else crossprod(X)
      eigen_R_cache <- eigen(R_for_eigen, symmetric = TRUE)
    }
  }

  # Convert to sufficient statistics format
  working <- summary_stats_working_quantities(
    z = z, n = n, shat = shat, var_y = var_y, pve_adjustment = adj,
    prior_variance = prior_variance,
    scaled_prior_variance = scaled_prior_variance)
  XtX <- NULL
  if (is.null(n)) {
    # Sample size not provided - use unadjusted z-scores
    hint_sample_size_recommended()
    if (!is.null(R)) {
      XtX <- R
    }
    # X path: X'X = R already after standardize_X, no further scaling needed.
    Xty <- working$Xty
    yty <- working$yty
    n <- working$n
    scaled_prior_variance <- working$scaled_prior_variance
  } else if (working$path == "original_scale") {
    # Sample size provided - use PVE-adjusted z-scores
    # var_y and shat provided - effects on original scale (R path only)
    XtX <- t(R * sqrt(working$XtXdiag)) * sqrt(working$XtXdiag)
    XtX <- (XtX + t(XtX)) / 2
    Xty <- working$Xty
    yty <- working$yty
  } else {
    # Effects on standardized X scale.
    if (!is.null(R)) {
      XtX <- working$XtX_multiplier * R
    } else {
      # X path: X'X = R after standardize_X, scale to X'X = (n-1)*R
      X <- X * working$X_multiplier
    }
    Xty <- working$Xty
    yty <- working$yty
  }

  # Use sufficient_stats_constructor with ALL parameters
  result <- sufficient_stats_constructor(
    Xty = Xty, yty = yty, n = n, XtX = XtX, X = X,
    L = L, X_colmeans = NA, y_mean = NA,
    maf = NULL, maf_thresh = 0, check_input = check_input,
    r_tol = r_tol, standardize = standardize,
    scaled_prior_variance = scaled_prior_variance,
    residual_variance = residual_variance, prior_weights = prior_weights,
    null_weight = null_weight, model_init = model_init,
    estimate_residual_variance = estimate_residual_variance,
    estimate_residual_method = estimate_residual_method,
    residual_variance_lowerbound = residual_variance_lowerbound,
    residual_variance_upperbound = residual_variance_upperbound,
    estimate_prior_variance = estimate_prior_variance,
    estimate_prior_method = estimate_prior_method,
    prior_variance_grid = prior_variance_grid,
    mixture_weights = mixture_weights,
    unmappable_effects = unmappable_effects,
    check_null_threshold = check_null_threshold, prior_tol = prior_tol,
    max_iter = max_iter, tol = tol, convergence_method = convergence_method,
    coverage = coverage, min_abs_corr = min_abs_corr, n_purity = n_purity,
    verbose = verbose, track_fit = track_fit, check_prior = check_prior,
    refine = refine, alpha0 = alpha0, beta0 = beta0,
    slot_prior = slot_prior, L_greedy = L_greedy,
    greedy_lbf_cutoff = greedy_lbf_cutoff
  )

  # Attach R-uncertainty metadata. For R_mismatch without finite-reference
  # input, R_finite_B = Inf gives the B^{-1} = 0 limit of the same model.
  if (!is.null(R_finite_B)) {
    result$data$R_finite_B <- R_finite_B
    result$data$R_finite_diagnostics <- R_finite_diagnostics
    result$data$R_mismatch <- R_mismatch
  }

  # eigen(R) cache for Q_art diagnostic.
  if (!is.null(eigen_R_cache))
    result$data$eigen_R <- eigen_R_cache

  # Attach R-mismatch params consumed by R/rss_mismatch.R.
  result$params$R_mismatch <- R_mismatch
  result$params$R_mismatch_method <- R_mismatch_method
  result$params$eig_delta_rel <- eig_delta_rel
  result$params$eig_delta_abs <- eig_delta_abs
  result$params$artifact_threshold <- artifact_threshold
  result$params$R_sensitivity_threshold <- R_sensitivity_threshold

  return(result)
}

# =============================================================================
# SS MULTI-PANEL MIXTURE DATA CONSTRUCTOR
# =============================================================================
#'
#' @keywords internal
#' @noRd
ss_mixture_constructor <- function(z, R = NULL, X = NULL, n,
                                   L = min(10, if (!is.null(R)) ncol(R[[1]])
                                           else ncol(X[[1]])),
                                   maf = NULL, maf_thresh = 0,
                                   scaled_prior_variance = 0.2,
                                   residual_variance = NULL,
                                   prior_weights = NULL,
                                   null_weight = 0,
                                   standardize = TRUE,
                                   estimate_residual_variance = FALSE,
                                   estimate_residual_method = "MoM",
                                   estimate_prior_variance = TRUE,
                                   estimate_prior_method = "optim",
                                   prior_variance_grid = NULL,
                                   mixture_weights = NULL,
                                   unmappable_effects = "none",
                                   check_null_threshold = 0,
                                   prior_tol = 1e-9,
                                   residual_variance_lowerbound = 0,
                                   residual_variance_upperbound = Inf,
                                   model_init = NULL,
                                   coverage = 0.95,
                                   min_abs_corr = 0.5,
                                   max_iter = 100,
                                   tol = NULL,
                                   convergence_method = "pip",
                                   verbose = FALSE,
                                   track_fit = FALSE,
                                   check_input = FALSE,
                                   check_prior = TRUE,
                                   n_purity = 100,
                                   r_tol = 1e-8,
                                   refine = FALSE,
                                   R_finite = NULL,
                                   R_mismatch = "none",
                                   R_mismatch_method = "mle",
                                   eig_delta_rel = 1e-3,
                                   eig_delta_abs = 0,
                                   artifact_threshold = 0.1,
                                   R_sensitivity_threshold = log(20),
                                   alpha0 = NULL,
                                   beta0 = NULL,
                                   slot_prior = NULL,
                                   L_greedy = NULL,
                                   greedy_lbf_cutoff = 0.1,
                                   init_panel = NULL) {
  if (is.null(n) || !is.numeric(n) || length(n) != 1 || n <= 1)
    stop("Sample size 'n' is required for multi-panel mode.")
  R_mismatch <- match.arg(R_mismatch, c("none", "eb", "eb_ser_init", "eb_force_init", "eb_no_init"))
  R_mismatch_method <- match.arg(R_mismatch_method, c("mle", "map"))
  if (!is.numeric(R_sensitivity_threshold) ||
      length(R_sensitivity_threshold) != 1L ||
      !is.finite(R_sensitivity_threshold) ||
      R_sensitivity_threshold < 0)
    stop("R_sensitivity_threshold must be a single nonnegative finite numeric.")
  if (is.null(z))
    stop("Multi-panel mode requires z-scores.")
  if (!is.null(R) && !is.null(X))
    stop("Please provide either R or X, but not both.")

  use_R <- !is.null(R)
  panels <- if (use_R) R else X
  K <- length(panels)
  if (K < 1)
    stop("Multi-panel input must contain at least one panel.")
  p <- length(z)

  if (use_R) {
    for (k in seq_len(K)) {
      if (!is.matrix(R[[k]]) || !is.numeric(R[[k]]))
        stop("Each element of R list must be a numeric matrix.")
      if (nrow(R[[k]]) != p || ncol(R[[k]]) != p)
        stop("Each element of R list must have dimension length(z) by length(z).")
      if (!is_symmetric_matrix(R[[k]]))
        R[[k]] <- (R[[k]] + t(R[[k]])) / 2
    }
    panel_R <- lapply(R, safe_cov2cor)
    X_list <- NULL
    B_list <- R_finite
    if (is.null(init_panel)) init_panel <- attr(R, ".init_panel")
    omega_cache <- NULL
  } else {
    for (k in seq_len(K)) {
      if (!is.matrix(X[[k]]) || !is.numeric(X[[k]]))
        stop("Each element of X list must be a numeric matrix.")
      if (ncol(X[[k]]) != p)
        stop("Each element of X list must have length(z) columns.")
    }
    X_list <- lapply(X, standardize_X)
    panel_R <- lapply(X_list, function(Xk) cov2cor(crossprod(Xk)))
    B_list <- if (is.null(R_finite)) NULL else R_finite
    if (is.null(init_panel)) init_panel <- attr(X, ".init_panel")
    omega_cache <- if (sum(vapply(X_list, nrow, integer(1))) < p)
                     precompute_omega_cache(X_list, z) else NULL
  }

  if (!is.null(maf)) {
    if (length(maf) != p)
      stop(paste0("The length of maf does not agree with expected ", p, "."))
    id <- which(maf > maf_thresh)
    z <- z[id]
    panel_R <- lapply(panel_R, function(Rk) Rk[id, id, drop = FALSE])
    if (!is.null(X_list))
      X_list <- lapply(X_list, function(Xk) Xk[, id, drop = FALSE])
    if (!is.null(prior_weights))
      prior_weights <- prior_weights[id]
    p <- length(z)
  }

  if (any(is.infinite(z)))
    stop("z contains infinite values.")
  if (anyNA(z)) {
    z <- replace_na_zero_warned(z, "z-scores")
  }

  nw <- normalize_null_weight(null_weight, prior_weights, p)
  null_weight <- nw$null_weight
  prior_weights <- nw$prior_weights

  if (nw$add_null) {
    panel_R <- lapply(panel_R, function(Rk) cbind(rbind(Rk, 0), 0))
    if (!is.null(X_list))
      X_list <- lapply(X_list, function(Xk) cbind(Xk, 0))
    z <- c(z, 0)
    p <- p + 1L
  }

  prior_weights <- normalize_prior_weights(prior_weights, p)

  k_best <- if (!is.null(init_panel)) init_panel else 1L
  omega_init <- rep(0, K)
  omega_init[k_best] <- 1
  R_init <- Reduce("+", Map(function(w, Rk) w * Rk, omega_init, panel_R))
  R_init <- 0.5 * (R_init + t(R_init))

  R_finite_B <- NULL
  R_finite_diagnostics <- NULL
  if (!is.null(R_finite)) {
    B_list <- as.numeric(R_finite)
    R_finite_B <- 1 / sum(omega_init^2 / B_list)
  }
  if (R_mismatch != "none" && is.null(R_finite_B))
    R_finite_B <- Inf
  if (!is.null(R_finite_B)) {
    R_finite_diagnostics <- compute_R_finite_diagnostics(
      R = R_init, B = R_finite_B, p = p)
  }

  nm1 <- n - 1
  XtX <- nm1 * R_init
  X_ss <- NULL
  if (!is.null(X_list)) {
    X_ss <- form_X_meta(X_list, omega_init) * sqrt(nm1)
    attr(X_ss, "d") <- rep(nm1, p)
    attr(X_ss, "scaled:scale") <- rep(1, p)
    XtX <- NULL
  }

  params_object <- list(
    L = L,
    scaled_prior_variance = scaled_prior_variance,
    residual_variance = residual_variance,
    prior_weights = prior_weights,
    null_weight = null_weight,
    estimate_residual_variance = estimate_residual_variance,
    estimate_residual_method = estimate_residual_method,
    residual_variance_lowerbound = residual_variance_lowerbound,
    residual_variance_upperbound = residual_variance_upperbound,
    estimate_prior_variance = estimate_prior_variance,
    estimate_prior_method = estimate_prior_method,
    prior_variance_grid = prior_variance_grid,
    mixture_weights = mixture_weights,
    unmappable_effects = unmappable_effects,
    check_null_threshold = check_null_threshold,
    prior_tol = prior_tol,
    max_iter = max_iter,
    tol = tol,
    convergence_method = convergence_method,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    compute_univariate_zscore = FALSE,
    verbose = verbose,
    track_fit = track_fit,
    check_prior = check_prior,
    refine = refine,
    n_purity = n_purity,
    alpha0 = alpha0,
    beta0 = beta0,
    n = n,
    use_NIG = estimate_residual_method == "NIG",
    intercept = FALSE,
    standardize = standardize,
    model_init = model_init,
    slot_prior = slot_prior,
    L_greedy = L_greedy,
    greedy_lbf_cutoff = greedy_lbf_cutoff,
    R_mismatch = R_mismatch,
    R_mismatch_method = R_mismatch_method,
    eig_delta_rel = eig_delta_rel,
    eig_delta_abs = eig_delta_abs,
    artifact_threshold = artifact_threshold,
    R_sensitivity_threshold = R_sensitivity_threshold
  )
  params_object <- validate_and_override_params(params_object)

  data_object <- structure(
    list(
      X = X_ss, XtX = XtX,
      Xty = sqrt(nm1) * z, yty = nm1,
      n = n, p = p,
      X_colmeans = rep(0, p), y_mean = 0,
      nm1 = nm1, z = z, lambda = 0,
      R_finite_B = R_finite_B,
      R_finite_diagnostics = R_finite_diagnostics,
      R_mismatch = R_mismatch,
      X_list_std = X_list, B_list = B_list,
      K = K, panel_R = panel_R, omega_cache = omega_cache,
      omega_init = omega_init
    ),
    class = c("ss_mixture", "ss")
  )
  if (R_mismatch != "none")
    data_object$eigen_R <- eigen(R_init, symmetric = TRUE)

  list(data = data_object, params = params_object)
}

# =============================================================================
# RSS LAMBDA DATA CONSTRUCTOR
#
# Constructs data and params objects for SuSiE from RSS data using eigendecomposition
# (lambda >= 0).
# Handles eigen decomposition, MAF filtering, and specialized RSS-lambda preprocessing.
# =============================================================================
#'
#' @return A list containing:
#' \item{data}{A processed list containing z-scores, R matrix, eigen decomposition,
#' and RSS-lambda specific fields}
#' \item{params}{Validated params object with all input algorithm parameters}
#'
#' @keywords internal
#' @noRd
rss_lambda_constructor <- function(z, R = NULL, X = NULL, n = NULL,
                                   L = min(10, if (!is.null(R)) ncol(R) else ncol(X)),
                                   lambda = 0,
                                   maf = NULL,
                                   maf_thresh = 0,
                                   prior_variance = 50,
                                   residual_variance = NULL,
                                   prior_weights = NULL,
                                   null_weight = 0,
                                   intercept_value = 0,
                                   estimate_residual_variance = FALSE,
                                   estimate_residual_method = "MLE",
                                   estimate_prior_variance = TRUE,
                                   estimate_prior_method = "optim",
                                   prior_variance_grid = NULL,
                                   mixture_weights = NULL,
                                   check_null_threshold = 0,
                                   prior_tol = 1e-9,
                                   residual_variance_lowerbound = 0,
                                   model_init = NULL,
                                   coverage = 0.95,
                                   min_abs_corr = 0.5,
                                   max_iter = 100,
                                   tol = NULL,
                                   convergence_method = "elbo",
                                   verbose = FALSE,
                                   track_fit = FALSE,
                                   check_prior = TRUE,
                                   check_R = TRUE,
                                   check_z = FALSE,
                                   n_purity = 100,
                                   r_tol = 1e-8,
                                   refine = FALSE,
                                   slot_prior = NULL,
                                   L_greedy = NULL,
                                   greedy_lbf_cutoff = 0.1) {

  # Validate: exactly one of R or X must be provided.
  if (is.null(R) && is.null(X))
    stop("Please provide either R (correlation matrix) or X (factor matrix).")
  if (!is.null(R) && !is.null(X))
    stop("Please provide either R or X, but not both.")
  if (!identical(estimate_residual_method, "MLE")) {
    stop("RSS-lambda supports estimate_residual_method = \"MLE\" only.")
  }
  if (is.list(R) && !is.matrix(R))
    stop("rss_lambda_constructor() accepts only a single R matrix.")
  if (is.list(X) && !is.matrix(X))
    stop("rss_lambda_constructor() accepts only a single X matrix.")

  # Default max_iter for susie_rss_lambda-style fitting (see notes in
  # summary_stats_constructor). NULL from `susie_rss_lambda` triggers the
  # default and hint; direct constructor callers keep `max_iter = 100`.
  if (is.null(max_iter)) {
    max_iter <- 50
    warning_message("Setting max_iter = 50 for the SuSiE RSS-lambda model because ",
                    "slow convergence is often a sign of unstable summary-statistics ",
                    "fitting. To disable this message, explicitly set max_iter = 50 ",
                    "or another value in the susie_rss_lambda() call.",
                    style = "hint")
  }

  # PVE-adjust z when sample size is provided. Shrinks large z toward
  # zero to account for winner's curse. Same form as the SS path
  # (summary_stats_constructor); skipped when n is unavailable.
  if (!is.null(z) && !is.null(n) && is.numeric(n) && length(n) == 1 &&
      is.finite(n) && n > 1)
    z <- apply_pve_adjustment(z, n)$z

  if (is.null(X)) {
    # R path: validate R
    if (is.null(R))
      stop("Please provide either R or X for rss_lambda_constructor.")
    if (nrow(R) != length(z)) {
      stop(paste0(
        "The dimension of correlation matrix (", nrow(R), " by ",
        ncol(R), ") does not agree with expected (", length(z), " by ",
        length(z), ")."
      ))
    }
    if (!is_symmetric_matrix(R)) {
      R <- symmetrize_warned(R, "R")
    }
    if (!(is.double(R) & is.matrix(R)) & !inherits(R, "sparseMatrix")) {
      stop("Input R must be a double-precision matrix or a sparse matrix.")
    }
  } else {
    # Single-panel X path: validate X
    if (ncol(X) != length(z)) {
      stop(paste0(
        "The number of columns of X (", ncol(X),
        ") does not agree with expected (", length(z), ")."
      ))
    }
  }

  # MAF filter
  if (!is.null(maf)) {
    if (length(maf) != length(z)) {
      stop(paste0("The length of maf does not agree with expected ", length(z), "."))
    }
    id <- which(maf > maf_thresh)
    if (!is.null(R)) R <- R[id, id]
    if (!is.null(X)) X <- X[, id, drop = FALSE]
    z <- z[id]
    if (!is.null(prior_weights))
      prior_weights <- prior_weights[id]
  }

  if (any(is.infinite(z))) {
    stop("z contains infinite values.")
  }

  # Check for NAs
  if (!is.null(R) && anyNA(R)) {
    stop("R matrix contains missing values.")
  }

  # Replace NAs in z with zero
  if (anyNA(z)) {
    z <- replace_na_zero_warned(z, "z-scores")
  }

  p_cur <- if (!is.null(R)) ncol(R) else ncol(X)
  nw <- normalize_null_weight(null_weight, prior_weights, p_cur)
  null_weight <- nw$null_weight
  prior_weights <- nw$prior_weights

  if (nw$add_null) {
    if (!is.null(R)) R <- cbind(rbind(R, 0), 0)
    if (!is.null(X)) X <- cbind(X, 0)
    z <- c(z, 0)
  }

  # Determine p and set prior weights
  p <- if (!is.null(R)) ncol(R) else ncol(X)

  prior_weights <- normalize_prior_weights(prior_weights, p)

  # Eigen decomposition: from R or SVD of X
  if (!is.null(X)) {
    # Single-panel: standardize so X'X = R, then SVD
    X <- standardize_X(X)
    eigen_R <- eigen_from_X(X, p)
  } else {
    eigen_R <- eigen(R, symmetric = TRUE)
  }

  if (is.null(X) && check_R && any(eigen_R$values < -r_tol)) {
    stop(paste0(
      "The correlation matrix (", nrow(R), " by ", ncol(R),
      ") is not a positive semidefinite matrix. ",
      "The smallest eigenvalue is ", min(eigen_R$values),
      ". You can bypass this by \"check_R = FALSE\" which instead ",
      "sets negative eigenvalues to 0 to allow for continued ",
      "computations."
    ))
  }

  # Check whether z in space spanned by the non-zero eigenvectors of R
  if (is.null(X) && check_z) {
    proj <- check_projection(eigen_R, z, tol = r_tol,
                             projection_tol = r_tol)
    if (proj$low_eigen_count > 0) {
      if (!proj$status) {
        warning_message("Input z does not lie in the space of non-zero eigenvectors of R.")
      } else {
        message("Input z is in space spanned by the non-zero eigenvectors of R.\n")
      }
    }
  }

  # Set negative eigenvalues to zero
  eigen_R$values[eigen_R$values < r_tol] <- 0

  # Precompute V'z
  Vtz <- crossprod(eigen_R$vectors, z)

  # Compute Null-space z-score norm: ||z||^2 - ||V'z||^2.
  z_null_norm2 <- max(sum(z^2) - sum(Vtz^2), 0)

  # Handle lambda estimation
  if (identical(lambda, "estimate")) {
    colspace <- which(eigen_R$values > 0)
    if (length(colspace) == length(z)) {
      lambda <- 0
    } else {
      znull <- crossprod(eigen_R$vectors[, -colspace], z)
      lambda <- sum(znull^2) / length(znull)
    }
  }

  if (is.null(residual_variance)) {
    residual_variance <- 1 - lambda
  } else {
    residual_variance <- residual_variance - lambda
  }

  # Create params object with ALL algorithm parameters
  params_object <- list(
    L = L,
    scaled_prior_variance = prior_variance, # Use unscaled prior_variance for RSS-lambda
    residual_variance = residual_variance,
    prior_weights = prior_weights,
    null_weight = null_weight,
    estimate_residual_variance = estimate_residual_variance,
    estimate_residual_method = estimate_residual_method,
    estimate_prior_variance = estimate_prior_variance,
    estimate_prior_method = estimate_prior_method,
    prior_variance_grid = prior_variance_grid,
    mixture_weights = mixture_weights,
    unmappable_effects = "none", # RSS-lambda doesn't support unmappable effects
    check_null_threshold = check_null_threshold,
    prior_tol = prior_tol,
    residual_variance_upperbound = 1, # RSS constraint
    model_init = model_init,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    compute_univariate_zscore = FALSE,
    max_iter = max_iter,
    tol = tol,
    convergence_method = convergence_method,
    verbose = verbose,
    track_fit = track_fit,
    residual_variance_lowerbound = residual_variance_lowerbound,
    refine = refine,
    n_purity = n_purity,
    alpha0 = 0,  # RSS doesn't support NIG
    beta0 = 0,   # RSS doesn't support NIG
    n = n,
    use_NIG = FALSE,
    intercept = FALSE,  # RSS always uses intercept = FALSE
    standardize = FALSE, # Never standardize RSS-lambda
    check_prior = check_prior,
    slot_prior = slot_prior,
    L_greedy = L_greedy,
    greedy_lbf_cutoff = greedy_lbf_cutoff
  )

  # Validate params
  params_object <- validate_and_override_params(params_object)

  # Create data object with RSS-lambda specific fields. n is the GWAS
  # sample size (used by the PVE adjustment above and by any downstream
  # consumer that needs to know the GWAS size); we store NA_integer_ when
  # the caller did not supply it. p (the number of variants) is always
  # length(z).
  data_object <- structure(
    list(
      z = z,
      R = R,
      X = X,
      n = if (is.null(n)) NA_integer_ else as.integer(n),
      p = length(z),
      lambda = lambda,
      intercept_value = intercept_value,
      r_tol = r_tol,
      prior_variance = prior_variance,
      eigen_R = eigen_R,
      Vtz = Vtz,
      z_null_norm2 = z_null_norm2
    ),
    class = "rss_lambda"
  )

  return(list(data = data_object, params = params_object))
}
