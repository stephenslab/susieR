#' Individual-level Data Constructor
#'
#' Constructs a data object for SuSiE from individual-level data (X, y).
#' This internal function prepares and validates the input data for use
#' in the SuSiE algorithm.
#'
#' @param X An n by p matrix of covariates (or sparse matrix, or trend filtering matrix)
#' @param y An n-vector of response values
#' @param intercept If TRUE, center X and y (default TRUE)
#' @param standardize If TRUE, scale X to have unit variance (default TRUE)
#' @param na.rm If TRUE, remove samples with missing values in y (default FALSE)
#' @param prior_weights A p-vector of prior inclusion probabilities (default NULL for uniform)
#' @param null_weight Weight for the null component (default 0, no null)
#' @param unmappable_effects Method for handling unmappable effects: "none", "inf", or "ash"
#'
#' @return A list with class "individual" containing:
#' \item{X}{The (possibly centered/scaled) covariate matrix}
#' \item{y}{The (possibly centered) response vector}
#' \item{mean_y}{Original mean of y}
#' \item{n}{Number of samples}
#' \item{p}{Number of variables}
#' \item{prior_weights}{Normalized prior weights}
#' \item{null_weight}{Weight for null component}
#' \item{unmappable_effects}{Method for unmappable effects}
#' Plus additional fields added by configure_data()
#'
#' @keywords internal
#' @noRd
individual_data_constructor <- function(X, y, intercept = TRUE, standardize = TRUE,
                                        na.rm = FALSE, prior_weights = NULL,
                                        null_weight = 0, unmappable_effects = "none",
                                        estimate_residual_method = "MLE",
                                        convergence_method = "elbo",
                                        estimate_prior_method = "optim",
                                        alpha0 = 0, beta0 = 0, refine = FALSE) {

  # Validate input X
  if (!(is.double(X) & is.matrix(X)) &
      !inherits(X, "CsparseMatrix") &
      is.null(attr(X, "matrix.type"))) {
    stop("Input X must be a double-precision matrix, or a sparse matrix, or ",
         "a trend filtering matrix.")
  }

  if (anyNA(X)) {
    stop("X contains NA values.")
  }

  # Handle missing values in y
  if (anyNA(y)) {
    if (na.rm) {
      samples_kept <- which(!is.na(y))
      y <- y[samples_kept]
      X <- X[samples_kept, ]
    } else {
      stop("Input y must not contain missing values.")
    }
  }
  mean_y <- mean(y)

  # Handle null weights
  if (is.numeric(null_weight) && null_weight == 0) {
    null_weight <- NULL
  }

  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight)) {
      stop("Null weight must be numeric.")
    }
    if (null_weight < 0 || null_weight >= 1) {
      stop("Null weight must be between 0 and 1.")
    }

    if (is.null(prior_weights)) {
      prior_weights <- c(rep(1 / ncol(X) * (1 - null_weight), ncol(X)), null_weight)
    } else {
      prior_weights <- c(prior_weights * (1 - null_weight), null_weight)
    }

    # add the extra 0 column to X
    X <- cbind(X, 0)
  }

  # Store dimensions
  n <- nrow(X)
  p <- ncol(X)

  # Validate and normalize prior_weights
  if (!is.null(prior_weights)) {
    if (length(prior_weights) != p) {
      stop("Prior weights must have length p.")
    }
    if (all(prior_weights == 0)) {
      stop("Prior weight should be greater than 0 for at least one variable.")
    }
    prior_weights <- prior_weights / sum(prior_weights)
  }
  if (p > 1000 & !requireNamespace("Rfast", quietly = TRUE)) {
    warning_message("For an X with many columns, please consider installing ",
                    "the Rfast package for more efficient credible set (CS) ",
                    "calculations.",
                    style = "hint")
  }

  # Center y if intercept is included
  if (intercept) {
    y <- y - mean_y
  }

  # Compute and set X matrix attributes
  out <- compute_colstats(X, center = intercept, scale = standardize)
  attr(X, "scaled:center") <- out$cm
  attr(X, "scaled:scale") <- out$csd
  attr(X, "d") <- out$d

  # Override convergence method for unmappable effects
  if (unmappable_effects != "none") {
    warning("Unmappable effects models (inf/ash) do not have a well defined ELBO and require PIP convergence. ",
            "Setting convergence_method='pip'.\n")
    convergence_method <- "pip"
  }

  # Check for incompatible parameter combinations
  if (refine && unmappable_effects != "none") {
    stop("Refinement is not supported with unmappable effects (inf/ash) as it relies on ELBO, ",
         "which is not well-defined for these models. Please set refine = FALSE.")
  }

  # Handle Servin_Stephens parameters for small sample correction
  if (estimate_residual_method == "Servin_Stephens") {
    use_servin_stephens <- TRUE

    # Override convergence method
    if (convergence_method != "pip") {
      warning("Servin_Stephens method requires PIP convergence. Setting convergence_method='pip'.\n")
      convergence_method <- "pip"
    }

    # Override prior variance estimation method
    if (estimate_prior_method != "EM") {
      warning("Servin_Stephens method works better with EM. Setting estimate_prior_method='EM'.\n")
      estimate_prior_method <- "EM"
    }
  } else {
    use_servin_stephens <- FALSE
    alpha0 <- NULL
    beta0 <- NULL
  }

  data_object <- structure(
    list(
      X = X,
      y = y,
      mean_y = mean_y,
      n = n,
      p = p,
      prior_weights = prior_weights,
      null_weight = null_weight,
      unmappable_effects = unmappable_effects,
      use_servin_stephens = use_servin_stephens,
      alpha0 = alpha0,
      beta0 = beta0,
      convergence_method = convergence_method,
      estimate_prior_method = estimate_prior_method
    ),
    class = "individual"
  )

  # Configure data object for specified unmappable effects method
  data_object <- configure_data(data_object)

  return(data_object)
}

#' Sufficient Statistics Data Constructor
#'
#' Constructs a data object for SuSiE from sufficient statistics.
#' This internal function prepares and validates sufficient statistics
#' for use in the SuSiE algorithm without requiring individual-level data.
#'
#' @param XtX A p by p matrix of X'X
#' @param Xty A p-vector of X'y
#' @param yty A scalar of y'y
#' @param n Sample size
#' @param X_colmeans Column means of X (default NA)
#' @param y_mean Mean of y (default NA)
#' @param maf Minor allele frequencies (for genetic data, default NULL)
#' @param maf_thresh MAF threshold for filtering (default 0)
#' @param standardize If TRUE, standardize variables (default TRUE)
#' @param r_tol Tolerance for positive semi-definite check (default 1e-8)
#' @param check_input If TRUE, perform additional input validation (default FALSE)
#' @param prior_weights A p-vector of prior inclusion probabilities (default NULL for uniform)
#' @param null_weight Weight for the null component (default 0, no null)
#' @param unmappable_effects Method for handling unmappable effects: "none", "inf", or "ash"
#'
#' @return A list with class "ss" containing:
#' \item{XtX}{The (possibly standardized) X'X matrix}
#' \item{Xty}{The (possibly standardized) X'y vector}
#' \item{yty}{The y'y value}
#' \item{n}{Sample size}
#' \item{p}{Number of variables}
#' \item{prior_weights}{Normalized prior weights}
#' \item{null_weight}{Weight for null component}
#' \item{unmappable_effects}{Method for unmappable effects}
#' Plus additional fields added by configure_data()
#'
#' @keywords internal
#' @noRd
sufficient_stats_constructor <- function(XtX, Xty, yty, n, X_colmeans = NA,
                                         y_mean = NA, maf = NULL, maf_thresh = 0,
                                         standardize = TRUE, r_tol = 1e-8,
                                         check_input = FALSE, prior_weights = NULL,
                                         null_weight = 0, unmappable_effects = "none",
                                         estimate_residual_method = "MLE",
                                         convergence_method = "elbo",
                                         scaled_prior_variance = NULL,
                                         check_prior = NULL,
                                         residual_variance_upperbound = NULL) {

  # Validate required inputs
  if (missing(n)) {
    stop("n must be provided.")
  }
  if (n <= 1) {
    stop("n must be greater than 1.")
  }
  if (missing(XtX) || missing(Xty) || missing(yty)) {
    stop("XtX, Xty, yty must all be provided.")
  }

  if (!(is.double(XtX) && is.matrix(XtX)) &&
      !inherits(XtX, "CsparseMatrix")) {
    stop("XtX must be a numeric dense or sparse matrix.")
  }

  if (ncol(XtX) != length(Xty)) {
    stop(paste0(
      "The dimension of XtX (", nrow(XtX), " by ", ncol(XtX),
      ") does not agree with expected (", length(Xty), " by ",
      length(Xty), ")."
    ))
  }

  if (ncol(XtX) > 1000 & !requireNamespace("Rfast", quietly = TRUE)) {
    warning_message("For large R or large XtX, consider installing the ",
                    "Rfast package for better performance.",
                    style = "hint")
  }

  # Ensure XtX is symmetric
  if (!is_symmetric_matrix(XtX)) {
    warning("XtX not symmetric; using (XtX + t(XtX))/2.\n")
    XtX <- (XtX + t(XtX)) / 2
  }

  # Apply MAF filter if provided
  if (!is.null(maf)) {
    if (length(maf) != length(Xty)) {
      stop(paste("The length of maf does not agree with expected", length(Xty), "."))
    }
    id <- which(maf > maf_thresh)
    XtX <- XtX[id, id]
    Xty <- Xty[id]
  }

  # Additional validation
  if (any(is.infinite(Xty))) {
    stop("Input Xty contains infinite values.")
  }
  if (!(is.double(XtX) & is.matrix(XtX)) & !inherits(XtX, "CsparseMatrix")) {
    stop("Input XtX must be a double-precision matrix, or a sparse matrix.")
  }
  if (anyNA(XtX)) {
    stop("Input XtX matrix contains NAs.")
  }

  # Replace NAs in Xty with zeros
  if (anyNA(Xty)) {
    warning_message("NA values in Xty are replaced with 0.")
    Xty[is.na(Xty)] <- 0
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

  # Handle null weights
  if (is.numeric(null_weight) && null_weight == 0) {
    null_weight <- NULL
  }

  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight)) {
      stop("Null weight must be numeric.")
    }
    if (null_weight < 0 || null_weight >= 1) {
      stop("Null weight must be between 0 and 1.")
    }
    if (is.null(prior_weights)) {
      prior_weights <- c(rep(1 / ncol(XtX) * (1 - null_weight), ncol(XtX)), null_weight)
    } else {
      prior_weights <- c(prior_weights * (1 - null_weight), null_weight)
    }
    XtX <- cbind(rbind(XtX, 0), 0)
    Xty <- c(Xty, 0)
    if (length(X_colmeans) == 1) {
      X_colmeans <- rep(X_colmeans, p)
    }
    if (length(X_colmeans) != p) {
      stop("The length of X_colmeans does not agree with number of variables.")
    }
  }

  # Define p
  p <- ncol(XtX)

  # Validate and normalize prior_weights
  if (!is.null(prior_weights)) {
    if (length(prior_weights) != p) {
      stop("Prior weights must have length p.")
    }
    if (all(prior_weights == 0)) {
      stop("Prior weight should be greater than 0 for at least one variable.")
    }
    prior_weights <- prior_weights / sum(prior_weights)
  }

  # Standardize if requested
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

  if (length(X_colmeans) == 1) {
    X_colmeans <- rep(X_colmeans, p)
  }
  if (length(X_colmeans) != p) {
    stop(
      "`X_colmeans` length (", length(X_colmeans),
      ") does not match number of variables (", p, ")."
    )
  }

  # Handle Servin_Stephens parameters for small sample correction
  if (estimate_residual_method == "Servin_Stephens") {
    stop("Small sample correction not implemented for SS/RSS data.")
  }

  # Ash requires individual-level data and cannot work with sufficient statistics alone
  if (unmappable_effects == "ash") {
    stop("Adaptive shrinkage (ash) requires individual-level data and cannot be used with sufficient statistics. ",
         "Please provide X and y instead of XtX, Xty, and yty.")
  }

  # Assemble data object
  data_object <- structure(
    list(
      XtX = XtX,
      Xty = Xty,
      yty = yty,
      n = n,
      p = p,
      X_colmeans = X_colmeans,
      y_mean = y_mean,
      prior_weights = prior_weights,
      null_weight = null_weight,
      unmappable_effects = unmappable_effects,
      convergence_method = convergence_method,
      use_servin_stephens = FALSE,
      scaled_prior_variance = scaled_prior_variance,
      standardize = standardize,
      check_prior = check_prior,
      residual_variance_upperbound = residual_variance_upperbound
    ),
    class = "ss"
  )

  # Configure data object for specified unmappable effects method
  data_object <- configure_data(data_object)

  return(data_object)
}

#' Summary Statistics (RSS) Data Constructor
#'
#' Constructs a data object for SuSiE from summary statistics (RSS format).
#' This internal function converts z-scores and correlation matrix to sufficient
#' statistics format for use in the SuSiE algorithm.
#'
#' @param z A p-vector of z-scores
#' @param R A p by p correlation matrix
#' @param n Sample size (optional but highly recommended)
#' @param bhat Alternative summary data giving estimated effects (length p)
#' @param shat Alternative summary data giving standard errors (length p)
#' @param var_y Sample variance of y, defined as y'y/(n-1)
#' @param z_ld_weight Deprecated parameter for backwards compatibility
#' @param prior_weights A p-vector of prior inclusion probabilities (default NULL for uniform)
#' @param null_weight Weight for the null component (default 0, no null)
#' @param unmappable_effects Method for handling unmappable effects: "none", "inf", or "ash"
#' @param standardize If TRUE, scale to unit variance (default TRUE)
#' @param check_input If TRUE, perform additional input validation (default FALSE)
#' @param r_tol Tolerance for positive semi-definite check (default 1e-8)
#'
#' @return A list with class "ss" containing sufficient statistics
#'
#' @keywords internal
#' @noRd
summary_stats_constructor <- function(z = NULL, R, n = NULL, bhat = NULL,
                                      shat = NULL, var_y = NULL, lambda = 0,
                                      maf = NULL, maf_thresh = 0, z_ld_weight = 0,
                                      prior_weights = NULL, null_weight = 0,
                                      unmappable_effects = "none",
                                      standardize = TRUE, check_input = FALSE,
                                      check_R = TRUE, check_z = FALSE,
                                      r_tol = 1e-8, prior_variance = 50,
                                      scaled_prior_variance = 0.2,
                                      intercept_value = 0,
                                      estimate_residual_variance = FALSE,
                                      estimate_residual_method = "MLE",
                                      convergence_method = "elbo",
                                      check_prior = TRUE,
                                      residual_variance_upperbound = Inf) {

  # Check if this should use RSS-lambda path
  if (lambda != 0) {
    # Parameter validation for RSS-lambda
    if (!is.null(bhat) || !is.null(shat)) {
      stop("Parameters 'bhat' and 'shat' are not supported when lambda != 0. ",
           "Please provide z-scores directly.")
    }

    if (!is.null(var_y)) {
      stop("Parameter 'var_y' is not supported when lambda != 0.")
    }

    if (z_ld_weight != 0) {
      stop("Parameter 'z_ld_weight' is not supported when lambda != 0.")
    }

    if (!is.null(n)) {
      warning("Parameter 'n' is ignored when lambda != 0.\n")
    }

    # Delegate to rss_lambda_constructor (which will apply overrides)
    return(rss_lambda_constructor(
      z = z, R = R, maf = maf, maf_thresh = maf_thresh,
      lambda = lambda, prior_weights = prior_weights,
      null_weight = null_weight, check_R = check_R,
      check_z = check_z, r_tol = r_tol,
      prior_variance = prior_variance,
      intercept_value = intercept_value,
      estimate_residual_method = estimate_residual_method,
      convergence_method = convergence_method
    ))
  }

  # Parameter validation for standard RSS (lambda = 0)
  if (!is.null(maf) || maf_thresh != 0) {
    stop("Parameters 'maf' and 'maf_thresh' are only supported when lambda != 0.")
  }

  if (intercept_value != 0) {
    stop("Parameter 'intercept_value' is only supported when lambda != 0.")
  }

  # Issue warning for estimate_residual_variance if TRUE
  if (estimate_residual_variance && lambda == 0) {
    warning("For estimate_residual_variance = TRUE, please check ",
            "that R is the \"in-sample\" LD matrix; that is, the ",
            "correlation matrix obtained using the exact same data ",
            "matrix X that was used for the other summary ",
            "statistics. Also note, when covariates are included in ",
            "the univariate regressions that produced the summary ",
            "statistics, also consider removing these effects from ",
            "X before computing R.\n")
  }

  # Check input R
  if (is.null(z) && !is.null(bhat)) {
    p <- length(bhat)
  } else if (!is.null(z)) {
    p <- length(z)
  } else {
    stop("Please provide either z or (bhat, shat).")
  }

  if (nrow(R) != p) {
    stop(paste0(
      "The dimension of R (", nrow(R), " x ", ncol(R), ") does not ",
      "agree with expected (", p, " x ", p, ")."
    ))
  }

  # Check input n
  if (!is.null(n)) {
    if (n <= 1) {
      stop("n must be greater than 1.")
    }
  }

  # Check inputs z, bhat and shat
  if (sum(c(is.null(z), is.null(bhat) || is.null(shat))) != 1) {
    stop("Please provide either z or (bhat, shat), but not both.")
  }
  if (is.null(z)) {
    if (length(shat) == 1) {
      shat <- rep(shat, length(bhat))
    }
    if (length(bhat) != length(shat)) {
      stop("The lengths of bhat and shat do not agree.")
    }
    if (anyNA(bhat) || anyNA(shat)) {
      stop("bhat, shat cannot have missing values.")
    }
    if (any(shat <= 0)) {
      stop("shat cannot have zero or negative elements.")
    }
    z <- bhat / shat
  }
  if (length(z) < 1) {
    stop("Input vector z should have at least one element.")
  }
  z[is.na(z)] <- 0

  # When n is provided, compute the PVE-adjusted z-scores
  if (!is.null(n)) {
    adj <- (n - 1) / (z^2 + n - 2)
    z <- sqrt(adj) * z
  }

  # Modify R by z_ld_weight (deprecated but maintained for compatibility)
  if (z_ld_weight > 0) {
    warning("As of version 0.11.0, use of non-zero z_ld_weight is no longer ",
            "recommended.\n")
    R <- muffled_cov2cor((1 - z_ld_weight) * R + z_ld_weight * tcrossprod(z))
    R <- (R + t(R)) / 2
  }

  # Convert to sufficient statistics format
  if (is.null(n)) {
    # Sample size not provided - use unadjusted z-scores
    warning("Providing the sample size (n), or even a rough estimate of n, ",
            "is highly recommended. Without n, the implicit assumption is ",
            "n is large (Inf) and the effect sizes are small (close to zero).\n")
    XtX <- R
    Xty <- z
    yty <- 1
    n <- 2 # Arbitrary choice to ensure var(y) = yty/(n-1) = 1
  } else {
    # Sample size provided - use PVE-adjusted z-scores
    if (!is.null(shat) && !is.null(var_y)) {
      # var_y and shat provided - effects on original scale
      XtXdiag <- var_y * adj / (shat^2)
      XtX <- t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
      XtX <- (XtX + t(XtX)) / 2
      Xty <- z * sqrt(adj) * var_y / shat
      yty <- (n - 1) * var_y
    } else {
      # Effects on standardized X, y scale
      XtX <- (n - 1) * R
      Xty <- sqrt(n - 1) * z
      yty <- (n - 1) * (if (!is.null(var_y)) var_y else 1)
    }
  }

  # Use sufficient_stats_constructor to handle the rest
  data_object <- sufficient_stats_constructor(
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    X_colmeans = NA,
    y_mean = NA,
    standardize = standardize,
    r_tol = r_tol,
    check_input = check_input,
    prior_weights = prior_weights,
    null_weight = null_weight,
    unmappable_effects = unmappable_effects,
    estimate_residual_method = estimate_residual_method,
    convergence_method = convergence_method,
    scaled_prior_variance = scaled_prior_variance,
    check_prior = check_prior,
    residual_variance_upperbound = residual_variance_upperbound
  )

  return(data_object)
}

#' Constructor for RSS Lambda Data
#'
#' Constructs a data object for RSS data with lambda regularization.
#'
#' @param z A p-vector of z-scores
#' @param R A p by p correlation matrix
#' @param maf A p-vector of minor allele frequencies
#' @param maf_thresh MAF threshold for filtering (default 0)
#' @param lambda Regularization parameter for correlated errors (default 0)
#' @param prior_weights A p-vector of prior inclusion probabilities (default NULL for uniform)
#' @param null_weight Weight for the null component (default 0, no null)
#' @param check_R If TRUE, check that R is positive semi-definite (default TRUE)
#' @param check_z If TRUE, check that z lies in column space of R (default FALSE)
#' @param r_tol Tolerance for eigenvalue checks (default 1e-8)
#' @param prior_variance Prior variance on effect sizes (default 50)
#' @param intercept_value Intercept value for the model (default 0)
#'
#' @return A list with class "rss_lambda" containing:
#' \item{z}{The z-scores}
#' \item{R}{The correlation matrix}
#' \item{lambda}{The regularization parameter}
#' \item{eigen_R}{Eigen decomposition of R}
#' Plus additional fields for the model
#'
#' @keywords internal
#' @noRd
rss_lambda_constructor <- function(z, R, maf = NULL, maf_thresh = 0,
                                   lambda = 0, prior_weights = NULL,
                                   null_weight = 0, check_R = TRUE,
                                   check_z = FALSE, r_tol = 1e-8,
                                   prior_variance = 50, intercept_value = 0,
                                   estimate_residual_method = "MLE",
                                   convergence_method = "elbo") {

  # Validate that MoM is not requested for RSS with lambda != 0
  if (estimate_residual_method == "MoM") {
    stop("Method of Moments (MoM) variance estimation is not implemented for RSS with lambda > 0. ",
         "Please use estimate_residual_method = 'MLE' instead.")
  }

  if (estimate_residual_method == "Servin_Stephens") {
    stop("Servin Stephens small sample correction is not implemented for RSS with lambda > 0. ",
         "Please use estimate_residual_method = 'MLE' instead.")
  }

  # Check input R
  if (nrow(R) != length(z)) {
    stop(paste0(
      "The dimension of correlation matrix (", nrow(R), " by ",
      ncol(R), ") does not agree with expected (", length(z), " by ",
      length(z), ")."
    ))
  }
  if (!isSymmetric(R)) {
    stop("R is not a symmetric matrix.")
  }
  if (!(is.double(R) & is.matrix(R)) & !inherits(R, "CsparseMatrix")) {
    stop("Input R must be a double-precision matrix or a sparse matrix.")
  }

  # MAF filter
  if (!is.null(maf)) {
    if (length(maf) != length(z)) {
      stop(paste0("The length of maf does not agree with expected ", length(z), "."))
    }
    id <- which(maf > maf_thresh)
    R <- R[id, id]
    z <- z[id]
  }

  if (any(is.infinite(z))) {
    stop("z contains infinite values.")
  }

  # Check for NAs in R
  if (anyNA(R)) {
    stop("R matrix contains missing values.")
  }

  # Replace NAs in z with zero
  if (anyNA(z)) {
    warning("NA values in z-scores are replaced with 0.\n")
    z[is.na(z)] <- 0
  }

  # Handle null weight
  if (is.numeric(null_weight) && null_weight == 0) {
    null_weight <- NULL
  }
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight)) {
      stop("Null weight must be numeric.")
    }
    if (null_weight < 0 || null_weight >= 1) {
      stop("Null weight must be between 0 and 1.")
    }
    if (is.null(prior_weights)) {
      prior_weights <- c(rep(1 / ncol(R) * (1 - null_weight), ncol(R)), null_weight)
    } else {
      prior_weights <- c(prior_weights * (1 - null_weight), null_weight)
    }
    R <- cbind(rbind(R, 0), 0)
    z <- c(z, 0)
  }

  # Eigen decomposition for R
  p <- ncol(R)
  eigen_R <- eigen(R, symmetric = TRUE)

  if (check_R && any(eigen_R$values < -r_tol)) {
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
  if (check_z) {
    # Project z onto null space of R
    colspace <- which(eigen_R$values > r_tol)
    if (length(colspace) < length(z)) {
      znull <- crossprod(eigen_R$vectors[, -colspace], z)
      if (sum(znull^2) > r_tol * sum(z^2)) {
        warning("Input z does not lie in the space of non-zero eigenvectors of R.\n")
      } else {
        message("Input z is in space spanned by the non-zero eigenvectors of R.\n")
      }
    }
  }

  # Set negative eigenvalues to zero
  eigen_R$values[eigen_R$values < r_tol] <- 0

  # Handle lambda estimation
  if (identical(lambda, "estimate")) {
    colspace <- which(eigen_R$values > 0)
    if (length(colspace) == length(z)) {
      lambda <- 0
    } else {
      znull <- crossprod(eigen_R$vectors[, -colspace], z) # U2^T z
      lambda <- sum(znull^2) / length(znull)
    }
  }

  # Create data object with RSS-lambda specific parameter overrides
  data_object <- structure(
    list(
      z = z,
      R = R,
      n = length(z),
      p = length(z),
      prior_weights = prior_weights,
      null_weight = null_weight,
      lambda = lambda,
      intercept_value = intercept_value,
      r_tol = r_tol,
      prior_variance = prior_variance,
      eigen_R = eigen_R,
      convergence_method = convergence_method,
      # RSS-lambda specific parameter overrides
      scaled_prior_variance = prior_variance, # Use unscaled for RSS-lambda
      standardize = FALSE, # Never standardize RSS-lambda
      residual_variance_upperbound = 1, # RSS constraint
      check_prior = FALSE # Skip prior check for RSS-lambda
    ),
    class = "rss_lambda"
  )

  return(data_object)
}
