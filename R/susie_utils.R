# =============================================================================
# FUNDAMENTAL BUILDING BLOCKS
#
# Basic mathematical operations and utilities that serve as dependencies
# for other functions. These include matrix operations, statistical computations,
# and general-purpose helper functions.
#
# Functions: warning_message, muffled_corr, muffled_cov2cor, is_symmetric_matrix,
# apply_nonzeros, compute_colSds, compute_colstats
# =============================================================================

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

# Correlation function with specified warning muffled.
#' @importFrom stats cor
#' @keywords internal
muffled_corr <- function(x) {
  withCallingHandlers(cor(x), warning = function(w) {
    if (grepl("the standard deviation is zero", w$message)) {
      invokeRestart("muffleWarning")}
    })
}

# cov2cor function with specified warning muffled.
#' @importFrom stats cov2cor
#' @keywords internal
muffled_cov2cor <- function(x) {
  withCallingHandlers(cov2cor(x), warning = function(w) {
    if (grepl("had.*(0|non-positive).*NA entries.*result.*(dubious|doubtful)", w$message)){
      invokeRestart("muffleWarning")}
    })
}

# Check for symmetric matrix.
#' @keywords internal
is_symmetric_matrix <- function(x) {
  if (is.matrix(x) && is.numeric(x) && !isS4(x) &&
      requireNamespace("Rfast", quietly = TRUE)) {
    return(Rfast::is.symmetric(x))
  } else {
    return(Matrix::isSymmetric(x))
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

# Check whether b is in space spanned by the non-zero eigenvectors of A
#' @keywords internal
check_projection <- function(A, b) {
  if (is.null(attr(A, "eigen"))) {
    attr(A, "eigen") <- eigen(A, symmetric = TRUE)
  }
  v <- attr(A, "eigen")$values
  B <- attr(A, "eigen")$vectors[, v > .Machine$double.eps]
  msg <- all.equal(as.vector(B %*% crossprod(B, b)), as.vector(b),
                   check.names = FALSE
  )
  if (!is.character(msg)) {
    return(list(status = TRUE, msg = NA))
  } else {
    return(list(status = FALSE, msg = msg))
  }
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
      Xty = Xty,
      yty = yty,
      n = data$n,
      p = data$p,
      X_colmeans = X_colmeans,
      y_mean = data$mean_y
    ),
    class = "ss"
  )

  # Copy attributes from X to XtX
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

  # Validate unmappable_effects
  if (!params$unmappable_effects %in% c("none", "inf", "ash")) {
    stop("unmappable_effects must be one of 'none', 'inf', or 'ash'.")
  }

  # Override convergence method for unmappable effects
  if (params$unmappable_effects != "none") {
    if (params$convergence_method != "pip") {
      warning_message("Unmappable effects models (inf/ash) do not have a well defined ELBO and require PIP convergence. ",
              "Setting convergence_method='pip'.")
      params$convergence_method <- "pip"
    }
  }

  # Check for incompatible parameter combinations
  if (!is.null(params$refine) && params$refine && params$unmappable_effects != "none") {
    stop("Refinement is not supported with unmappable effects (inf/ash) as it relies on ELBO, ",
         "which is not well-defined for these models. Please set refine = FALSE.")
  }

  # Override prior estimation method when estimation is disabled
  if (!params$estimate_prior_variance) {
    params$estimate_prior_method <- "none"
  }

  # Handle Servin_Stephens parameters for small sample correction
  if (params$estimate_residual_method == "Servin_Stephens") {
    params$use_servin_stephens <- TRUE

    # Override convergence method only when L > 1
    if (params$L > 1 && params$convergence_method != "pip") {
      warning_message("Servin_Stephens method with L > 1 requires PIP convergence. Setting convergence_method='pip'.")
      params$convergence_method <- "pip"
    }

    # Override prior variance estimation method
    if (params$estimate_prior_method != "EM") {
      warning_message("Servin_Stephens method works better with EM. Setting estimate_prior_method='EM'.")
      params$estimate_prior_method <- "EM"
    }
  } else {
    params$use_servin_stephens <- FALSE
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
# Functions: initialize_matrices, initialize_null_index, assign_names,
# adjust_L, prune_single_effects, add_null_effect
# =============================================================================

# Initialize core susie model object with default parameter matrices
#' @keywords internal
initialize_matrices <- function(data, params, var_y) {
  L <- params$L
  mat_init <- list(
    alpha             = matrix(1 / data$p, L, data$p),
    mu                = matrix(0, L, data$p),
    mu2               = matrix(0, L, data$p),
    V                 = rep(params$scaled_prior_variance * var_y, L),
    KL                = rep(as.numeric(NA), L),
    lbf               = rep(as.numeric(NA), L),
    lbf_variable      = matrix(as.numeric(NA), L, data$p),
    sigma2            = params$residual_variance,
    pi                = params$prior_weights,
    null_weight       = params$null_weight,
    predictor_weights = rep(as.numeric(NA), data$p)
  )

  return(mat_init)
}

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
    ))
    L <- num_effects
  }

  V <- rep(params$scaled_prior_variance * var_y, L)
  model_init <- prune_single_effects(model_init_pruned, L = L, V = V)

  return(list(model_init = model_init, L = L))
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
        V[1:num_effects] <- model_init$V[effects_rank]
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
#' @keywords internal
compute_eigen_decomposition <- function(XtX, n) {
  LD  <- XtX / n
  eig <- eigen(LD, symmetric = TRUE)
  idx <- order(eig$values, decreasing = TRUE)

  list(
    V     = eig$vectors[, idx],
    Dsq   = pmax(eig$values[idx] * n, 0),
    VtXty = NULL
  )
}

# Add eigen decomposition to ss data objects for unmappable methods
#' @keywords internal
add_eigen_decomposition <- function(data, params, individual_data = NULL) {
  # Compute eigen decomposition
  eigen_decomp <- compute_eigen_decomposition(data$XtX, data$n)

  # Append eigen components to data object
  data$eigen_vectors <- eigen_decomp$V
  data$eigen_values  <- eigen_decomp$Dsq
  data$VtXty         <- t(eigen_decomp$V) %*% data$Xty

  if (params$unmappable_effects == "ash") {
    if (is.null(individual_data)) {
      stop("Adaptive shrinkage (ash) requires individual-level data")
    }

    X_scaled <- scale_design_matrix(
      individual_data$X,
      center = attr(individual_data$X, "scaled:center"),
      scale = attr(individual_data$X, "scaled:scale")
    )

    data$X    <- X_scaled               
    data$y    <- individual_data$y       
    data$VtXt <- t(data$eigen_vectors) %*% t(X_scaled)
  }

  return(data)
}

#' Scale design matrix using centering and scaling parameters
#'
#' Applies column-wise centering and scaling to match the space used by
#' compute_XtX() and compute_Xty() for unmappable effects methods.
#'
#' @param X Matrix to scale (n × p)
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

# Compute Omega-weighted quantities for unmappable effects methods
#' @keywords internal
compute_omega_quantities <- function(data, tau2, sigma2) {
  # Compute variance in eigen space
  omega_var <- tau2 * data$eigen_values + sigma2

  # Compute diagonal of X'OmegaX
  diagXtOmegaX <- rowSums(sweep(data$eigen_vectors^2, 2,
                                (data$eigen_values / omega_var), `*`))

  return(list(
    omega_var    = omega_var,
    diagXtOmegaX = diagXtOmegaX
  ))
}

# Compute unmappable effects coefficient vector using BLUP
#' @keywords internal
compute_theta_blup <- function(data, model) {
  # Calculate diagXtOmegaX, diagonal variances, and Beta
  omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
  b         <- colSums(model$mu * model$alpha)

  # Compute XtOmegaXb, XtOmegay, and XtOmegar
  XtOmegaXb <- data$eigen_vectors %*% ((t(data$eigen_vectors) %*% b) *
                                         data$eigen_values / omega_res$omega_var)
  XtOmegay  <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)
  XtOmegar  <- XtOmegay - XtOmegaXb

  # Compute theta
  theta     <- model$tau2 * XtOmegar

  return(theta)
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
compute_lbf_gradient <- function(alpha, betahat, shat2, V, use_servin_stephens = FALSE) {
  # No gradient computation for Servin-Stephens prior
  if (use_servin_stephens) {
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
# methods (MLE, MoM, Servin-Stephens). These handle both standard SuSiE
# and unmappable effects models.
#
# Functions: mom_unmappable, mle_unmappable, create_ash_grid,
# compute_lbf_servin_stephens, posterior_mean_servin_stephens,
# posterior_var_servin_stephens, compute_stats_NIG, update_prior_variance_NIG_EM,
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

  # Compute diag(V'MV)
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

  diagVtMV <- diagVtMV + rowSums(sweep(t(data$eigen_vectors)^2, 2, tmpD, `*`))

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

  # Set default ranges
  sigma2_range <- c(0.2 * data$yty / data$n, 1.2 * data$yty / data$n)
  tau2_range   <- c(1e-12, 1.2 * data$yty / (data$n * data$p))

  # Compute diag(V'MV)
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

  diagVtMV <- diagVtMV + rowSums(sweep(t(data$eigen_vectors)^2, 2, tmpD, `*`))

  # Negative ELBO as function of x = (sigma^2, tau^2)
  f <- function(x) {
    sigma2_val <- x[1]
    tau2_val   <- x[2]
    var_val    <- tau2_val * data$eigen_values + sigma2_val

    0.5 * (data$n - data$p) * log(sigma2_val) + 0.5 / sigma2_val * data$yty +
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

#' Create variance grid for mr.ash
#' @keywords internal
create_ash_grid <- function(data, grid_length = 25, exponent = 2) {
  n <- data$n
  w <- colSums(data$X^2)
  sa2 <- (2^((0:(grid_length - 1)) / grid_length) - 1)^exponent * n / median(w)
  return(sa2)
}

#' Diagnostic function for SuSiE-ash iteration
#' @keywords internal
diagnose_susie_ash_iter <- function(data, model, params, Xcorr, mrash_output, 
                                     residuals, b_confident, effect_purity,
                                     pip_protected, neighborhood_pip, contested,
                                     sigma2_new, tau2_new, sa2) {
  
  L <- nrow(model$alpha)
  p <- ncol(model$alpha)
  purity_threshold <- if (!is.null(params$purity_threshold)) params$purity_threshold else 0.5
  
  theta_sum2_before <- sum(mrash_output$beta^2)
  theta_max_before <- max(abs(mrash_output$beta))
  n_high_purity <- sum(effect_purity >= purity_threshold, na.rm = TRUE)
  n_low_purity <- sum(!is.na(effect_purity)) - n_high_purity
  b_full <- colSums(model$alpha * model$mu)
  grid_scale <- data$n / median(colSums(data$X^2))
  pi_null <- mrash_output$pi[1]
  pi_nonnull <- 1 - pi_null
  theta_new <- mrash_output$beta
  theta_new[contested] <- 0
  
  cat(sprintf("SuSiE-ash iter %d:\n", model$ash_iter))
  cat(sprintf("  Purity: %d high (>=%.2f), %d low | current: %s\n", 
              n_high_purity, purity_threshold, n_low_purity,
              paste(round(effect_purity[!is.na(effect_purity)], 2), collapse=", ")))
  cat(sprintf("  Effect lbf: %s\n", paste(round(model$lbf, 1), collapse=", ")))
  cat(sprintf("  Residuals: var(y-Xb_conf)=%.4f | b_full²=%.2e, b_conf²=%.2e\n",
              var(as.vector(residuals)), sum(b_full^2), sum(b_confident^2)))
  cat(sprintf("  Grid: scale=%.3f (n/med(w)) | sa2 range=[%.2e, %.2e]\n",
              grid_scale, min(sa2), max(sa2)))
  cat(sprintf("  Mr.ASH: sigma2=%.4f, tau2=%.2e | pi_null=%.1f%%, pi_nonnull=%.1f%%\n",
              sigma2_new, tau2_new, pi_null*100, pi_nonnull*100))
  cat(sprintf("  Theta before: sum²=%.2e, max|θ|=%.4f | after exclusion: sum²=%.2e\n",
              theta_sum2_before, theta_max_before, sum(theta_new^2)))
  cat(sprintf("  LD-exclusion: %d total contested\n", sum(contested)))
  cat(sprintf("  Pre-MrASH: var(residuals)=%.6f, sum(residuals)=%.6f, residuals[1:3]=%.4f,%.4f,%.4f\n",
            var(as.vector(residuals)), sum(residuals), residuals[1], residuals[2], residuals[3]))

  # Top Mr.ASH signals before zeroing
  top_theta_idx <- order(abs(mrash_output$beta), decreasing = TRUE)[1:10]
  cat(sprintf("  Top 10 |theta| before zeroing:\n"))
  for (k in 1:10) {
    idx <- top_theta_idx[k]
    cat(sprintf("    %d: theta=%.4f, contested=%s, pip_prot=%.3f, nbhd_pip=%.3f\n",
                idx, mrash_output$beta[idx], contested[idx], 
                pip_protected[idx], neighborhood_pip[idx]))
  }
  
  # Contested variants with non-trivial Mr.ASH signal
  contested_with_signal <- which(contested & abs(mrash_output$beta) > 0.001)
  if (length(contested_with_signal) > 0) {
    cat(sprintf("  Contested variants with |theta|>0.001 (%d total):\n", length(contested_with_signal)))
    top_contested <- contested_with_signal[order(abs(mrash_output$beta[contested_with_signal]), decreasing=TRUE)]
    for (k in 1:min(10, length(top_contested))) {
      idx <- top_contested[k]
      cat(sprintf("    %d: theta=%.4f, nbhd_pip=%.3f\n",
                  idx, mrash_output$beta[idx], neighborhood_pip[idx]))
    }
  }
  
  # For each active effect, show nearby theta and potentially missed variants
  cat(sprintf("  Mr.ASH near CS sentinels:\n"))
  for (l in 1:L) {
    max_alpha_l <- max(model$alpha[l,])
    if (max_alpha_l - min(model$alpha[l,]) < 1E-6) next
    sentinel_l <- which.max(model$alpha[l,])
    
    # Variants in moderate LD with sentinel
    near_sentinel <- which(abs(Xcorr[sentinel_l,]) > 0.3)
    if (length(near_sentinel) > 0) {
      theta_near <- mrash_output$beta[near_sentinel]
      top_near <- near_sentinel[order(abs(theta_near), decreasing = TRUE)[1:min(5, length(near_sentinel))]]
      cat(sprintf("    Effect %d (sentinel=%d): nearby theta = ", l, sentinel_l))
      cat(paste0(top_near, "(", round(mrash_output$beta[top_near], 4), ", r=", 
                 round(Xcorr[sentinel_l, top_near], 2), ")", collapse=", "))
      cat("\n")
    }
  }

  # Check marginal z-score at variants with Mr.ASH signal but weak SuSiE
  cat(sprintf("  Marginal z-scores in residuals for top theta variants:\n"))
  for (k in 1:min(5, length(top_theta_idx))) {
    idx <- top_theta_idx[k]
    x_j <- data$X[, idx]
    z_resid <- as.numeric(crossprod(x_j, residuals)) / sqrt(sum(x_j^2) * sigma2_new)
    cat(sprintf("    %d: z=%.2f, theta=%.4f, contested=%s\n",
                idx, z_resid, mrash_output$beta[idx], contested[idx]))
  }
  
  # Effect summary
  for (l in 1:L) {
    max_alpha_l <- max(model$alpha[l,])
    if (max_alpha_l - min(model$alpha[l,]) < 1E-6) next
    sentinel_l <- which.max(model$alpha[l,])
    cat(sprintf("  Effect %d: max_alpha=%.3f, purity=%.2f, min_purity=%.2f, sentinel=%d, lbf=%.1f\n",
                l, max_alpha_l, effect_purity[l], model$min_purity[l], sentinel_l, model$lbf[l]))
  }
}

# Compute log Bayes factor for Servin and Stephens prior
#' @keywords internal
compute_lbf_servin_stephens <- function(x, y, s0, alpha0 = 0, beta0 = 0) {
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

# Posterior mean for Servin and Stephens prior using sufficient statistics
#' @keywords internal
posterior_mean_servin_stephens <- function(xtx, xty, s0_t = 1) {
  omega <- (xtx + (1 / s0_t^2))^(-1)
  b_bar <- omega %*% xty
  return(b_bar)
}

# Posterior variance for Servin and Stephens prior using sufficient statistics
#' @keywords internal
posterior_var_servin_stephens <- function(xtx, xty, yty, n, s0_t = 1) {

  # If prior variance is too small, return 0.
  if (s0_t < 1e-5) {
    return(list(post_var = 0, beta1 = 0))
  }

  omega <- (xtx + (1 / s0_t^2))^(-1)
  b_bar <- omega %*% xty
  beta1 <- (yty - b_bar * (omega^(-1)) * b_bar)
  post_var_up <- 0.5 * (yty - b_bar * (omega^(-1)) * b_bar)
  post_var_down <- 0.5 * (n * (1 / omega))
  post_var <- omega * (post_var_up / post_var_down) * n / (n - 2)

  return(list(post_var = post_var, beta1 = beta1))
}

# Compute the (log) Bayes factors and additional statistics under Normal-Inverse-Gamma (NIG) prior
#' @keywords internal
compute_stats_NIG <- function(n, xx, xy, yy, sxy, s0, a0, b0) {

  r0 <- s0 / (s0 + 1 / xx)
  rss <- yy * (1 - r0 * sxy^2)

  # Update inverse-gamma parameters
  a1 <- a0 + n
  b1 <- b0 + rss

  # Compute log Bayes factor for each variable
  lbf <- -(log(1 + s0 * xx) + a1 * log(b1 / (b0 + yy))) / 2

  # Compute least-squares estimate for each variable
  bhat <- xy / xx

  # Compute posterior mean
  post_mean <- r0 * bhat

  # Compute posterior variance
  post_var <- b1 / (a1 - 2) * r0 / xx

  # Compute posterior mode of residual variance
  rv <- (b1 / 2) / (a1 / 2 - 1)

  return(list(
    lbf        = lbf,
    post_mean  = post_mean,
    post_mean2 = post_var + post_mean^2,
    post_var   = post_var,
    rv         = rv
  ))
}

# Compute log Bayes factors under Normal-Inverse-Gamma (NIG) prior
#' @keywords internal
compute_lbf_NIG <- function(n, xx, xy, yy, sxy, s0, a0, b0) {
  r0 <- s0 / (s0 + 1 / xx)
  rss <- yy * (1 - r0 * sxy^2)

  # Update inverse-gamma parameters
  a1 <- a0 + n
  b1 <- b0 + rss

  # Compute log Bayes factor for each variable
  lbf <- -(log(1 + s0 * xx) + a1 * log(b1 / (b0 + yy))) / 2

  return(lbf)
}

# Compute posterior moments under Normal-Inverse-Gamma (NIG) prior
#' @keywords internal
compute_posterior_moments_NIG <- function(n, xx, xy, yy, sxy, s0, a0, b0) {
  r0 <- s0 / (s0 + 1 / xx)
  rss <- yy * (1 - r0 * sxy^2)

  # Update inverse-gamma parameters
  a1 <- a0 + n
  b1 <- b0 + rss

  # Compute least-squares estimate for each variable
  bhat <- xy / xx

  # Compute posterior mean
  post_mean <- r0 * bhat

  # Compute posterior variance
  post_var <- b1 / (a1 - 2) * r0 / xx

  # Compute posterior mode of residual variance
  rv <- (b1 / 2) / (a1 / 2 - 1)

  return(list(
    post_mean  = post_mean,
    post_mean2 = post_var + post_mean^2,
    post_var   = post_var,
    rv         = rv
  ))
}

# EM update for prior variance under Normal-Inverse-Gamma (NIG) prior
#' @keywords internal
update_prior_variance_NIG_EM <- function(n, xx, xy, yy, sxy, pip, s0, a0, b0) {
  r0   <- s0 / (s0 + 1 / xx)
  rss  <- yy * (1 - r0 * sxy^2)

  # Update inverse-gamma parameters
  a1   <- a0 + n
  b1   <- b0 + rss

  # Compute posterior mean and variance component
  bhat <- xy / xx
  post_mean  <- r0 * bhat
  post_var   <- r0 / xx

  u  <- gamma(1/2) / beta(a1/2, 1/2)
  mb <- post_mean * sqrt(2 / b1) * u
  vb <- post_var + post_mean^2 * 2 / b1 * (1 / beta(a1/2, 1) - u^2)

  return(sum(pip * (vb + mb^2)))
}

# Compute KL divergence for Normal-Inverse-Gamma (NIG) prior
#' @keywords internal
compute_kl_NIG <- function(alpha, post_mean, post_mean2, pi, V, a0, b0, a_post, b_post) {
  eps <- .Machine$double.eps

  # Posterior variance from second moment
  post_var <- pmax(post_mean2 - post_mean^2, eps)

  # Prior precision (tau2 = 1/V)
  tau2 <- 1 / V

  # KL for gamma
  KL_gamma <- sum(alpha * (log(pmax(alpha, eps)) - log(pmax(pi, eps))))

  # Expectations under posterior q(sigma^2) ~ IG(a_post, b_post)
  E_log_sigma2 <- digamma(a_post) - log(b_post)
  E_inv_sigma2 <- a_post / b_post

  # KL divergence for beta given sigma^2
  KL_beta <- 0.5 * sum(alpha * (
    E_log_sigma2 + log(tau2) - log(pmax(post_var, eps)) +
      E_inv_sigma2 * (post_var + post_mean^2) / tau2 - 1
  ))

  # KL divergence between IG posterior and IG prior
  KL_sigma2 <- lgamma(a0) - lgamma(a_post) +
    a0 * log(b_post / b0) +
    (a_post - a0) * digamma(a_post) -
    a_post + (a_post * b0) / b_post

  # Total KL divergence
  KL_total <- KL_gamma + KL_beta + KL_sigma2

  return(as.numeric(KL_total))
}

# Compute log-normalizing factor for the IG(a,b) distribution
#' @keywords internal
inv_gamma_factor <- function(a, b) {
  return(a * log(b) - lgamma(a))
}

# Compute null log-likelihood under NIG prior
#' @keywords internal
compute_null_loglik_NIG <- function(n, yy, a0, b0, use_servin_stephens = FALSE) {
  # No null log-likelihood for non-Servin-Stephens prior
  if (!use_servin_stephens) {
    return(NULL)
  }

  return(-n * log(2 * pi) / 2 +
         inv_gamma_factor(a0 / 2, b0 / 2) -
         inv_gamma_factor((a0 + n) / 2, (b0 + yy) / 2))
}

# Compute marginal log-likelihood for single effect regression
#' @keywords internal
compute_marginal_loglik <- function(lbf_model, n, yy, a0, b0, use_servin_stephens = FALSE) {
  # No marginal log-likelihood computation for non-Servin-Stephens prior
  if (!use_servin_stephens) {
    return(NULL)
  }

  ll0 <- compute_null_loglik_NIG(n, yy, a0, b0, use_servin_stephens = TRUE)
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

# Helper function to update variance components and derived quantities
#' @keywords internal
#' @importFrom utils modifyList
update_model_variance <- function(data, params, model) {
  # Update variance components
  variance_result <- update_variance_components(data, params, model)
  model           <- modifyList(model, variance_result)

  # Apply bounds to residual variance
  model$sigma2    <- min(max(model$sigma2, params$residual_variance_lowerbound),
                         params$residual_variance_upperbound)

  # Update derived quantities after variance component changes
  model           <- update_derived_quantities(data, params, model)

  return(model)
}

# =============================================================================
# CONVERGENCE & OPTIMIZATION
#
# Functions related to iteration control, convergence checking, and objective
# function computation. These control the iterative fitting process and
# determine when the algorithm has converged.
#
# Functions: check_convergence, get_objective, compute_elbo_inf
# =============================================================================

# Check convergence
#' @keywords internal
check_convergence <- function(params, model, elbo, iter, tracking) {
  # Skip convergence check on first iteration
  if(iter == 1) {
    return(FALSE)
  }

  # Calculate difference in ELBO values
  ELBO_diff   <- elbo[iter + 1] - tracking$convergence$prev_elbo
  ELBO_failed <- is.na(ELBO_diff) || is.infinite(ELBO_diff)

  if (params$convergence_method == "pip" || ELBO_failed) {
    # Fallback to PIP-based convergence if ELBO calculation fails
    if (ELBO_failed && params$convergence_method == "elbo") {
      warning_message(paste0("Iteration ", iter, " produced an NA/infinite ELBO
                             value. Using pip-based convergence this iteration."))
    }

    # For Servin-Stephens prior, require at least 3 iterations and average convergence
    # over 2 consecutive iterations for more stable convergence
    if (!is.null(params$use_servin_stephens) && params$use_servin_stephens) {
      if (iter <= 2) {
        return(FALSE)  # Require at least 3 iterations
      }

      # Current iteration PIP difference
      current_diff <- max(abs(tracking$convergence$prev_alpha - model$alpha))

      # Average with previous iteration's difference if available
      if (!is.null(tracking$convergence$prev_pip_diff)) {
        avg_diff <- (current_diff + tracking$convergence$prev_pip_diff) / 2
      } else {
        avg_diff <- current_diff
      }

      # Store current diff for next iteration
      tracking$convergence$prev_pip_diff <- current_diff

      if (params$verbose) {
        message("max |change in PIP| (avg): ", format(avg_diff, digits = 6))
      }

      return(avg_diff < params$tol)
    } else {
      # Standard PIP convergence
      PIP_diff <- max(abs(tracking$convergence$prev_alpha - model$alpha))

      if (params$verbose) {
        message("max |change in PIP|: ", format(PIP_diff, digits = 6))
      }

      return(PIP_diff < params$tol)
    }
  }

  if (params$verbose) {
    message("ELBO: ", format(elbo[iter + 1], digits = 6))
  }

  return(ELBO_diff < params$tol)
}

# Objective function (ELBO)
#' @keywords internal
get_objective <- function(data, params, model) {
  if (!is.null(params$unmappable_effects) && params$unmappable_effects == "inf") {
    # Compute omega
    L         <- nrow(model$alpha)
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega     <- matrix(0, L, data$p)

    for (l in seq_len(L)) {
      omega[l, ] <- omega_res$diagXtOmegaX + 1 / model$V[l]
    }

    # Compute total ELBO for infinitesimal effects model
    objective <- compute_elbo_inf(
      model$alpha, model$mu, omega, model$lbf,
      model$sigma2, model$tau2, data$n, data$p,
      data$eigen_vectors, data$eigen_values,
      data$VtXty, data$yty
    )
  } else if (params$use_servin_stephens && nrow(model$alpha) == 1) {
    objective <- model$marginal_loglik[1]
  } else {
    # Standard ELBO computation
    objective <- Eloglik(data, model) - sum(model$KL)
  }

  if (is.infinite(objective)) {
    stop("get_objective() produced an infinite ELBO value")
  }
  return(objective)
}

# Compute ELBO for infinitesimal effects model
#' @keywords internal
compute_elbo_inf <- function(alpha, mu, omega, lbf, sigma2, tau2, n, p,
                             eigen_vectors, eigen_values, VtXty, yty) {
  L <- nrow(mu)

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

  diagVtMV <- diagVtMV + rowSums(sweep(t(eigen_vectors)^2, 2, tmpD, `*`))

  # Compute variance
  var <- tau2 * eigen_values + sigma2

  # Compute negative ELBO
  neg_elbo <- 0.5 * (n - p) * log(sigma2) + 0.5 / sigma2 * yty +
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
# Functions: n_in_CS_x, in_CS_x, n_in_CS, in_CS, get_purity, get_tracking
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
                       use_rfast) {
  if (missing(use_rfast)) {
    use_rfast <- requireNamespace("Rfast", quietly = TRUE)
  }
  if (use_rfast) {
    get_upper_tri <- Rfast::upper_tri
    get_median <- Rfast::med
  } else {
    get_upper_tri <- function(R) R[upper.tri(R)]
    get_median <- stats::median
  }
  if (length(pos) == 1) {
    return(c(1, 1, 1))
  } else {
    # Subsample the columns if necessary.
    if (length(pos) > n) {
      pos <- sample(pos, n)
    }

    if (is.null(Xcorr)) {
      X_sub <- X[, pos]
      X_sub <- as.matrix(X_sub)
      value <- abs(get_upper_tri(muffled_corr(X_sub)))
    } else {
      value <- abs(get_upper_tri(Xcorr[pos, pos]))
    }
    if (squared) {
      value <- value^2
    }
    return(c(
      min(value),
      sum(value) / length(value),
      get_median(value)
    ))
  }
}

# Clean tracking object for output by removing convergence data
#' @keywords internal
get_tracking <- function(tracking) {
  clean_tracking <- tracking
  clean_tracking$convergence <- NULL
  return(clean_tracking)
}
