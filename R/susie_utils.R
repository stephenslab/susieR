# Compute eigenvalue decomposition
#' @keywords internal
compute_eigen_decomposition <- function(XtX, n) {
  LD <- XtX / n
  eig <- eigen(LD, symmetric = TRUE)
  idx <- order(eig$values, decreasing = TRUE)

  list(
    V = eig$vectors[, idx],
    Dsq = pmax(eig$values[idx] * n, 0),
    VtXty = NULL
  )
}

# Call mr.ash for adaptive shrinkage unmappable effects
#' @keywords internal
call_mr_ash <- function(X, y, sigma2, tau2, alpha, mu,
                        V, Dsq, VtXt,
                        K_length = 10,
                        update_ash_sigma = FALSE,
                        custom_est_sa2 = NULL) {

  # Create variance grid for mr.ash
  if (is.null(custom_est_sa2)) {
    est_sa2 <- 100 * tau2 * (seq(0, 1, length.out = K_length))^2
  } else {
    est_sa2 <- custom_est_sa2
  }

  # Compute required quantities for mr.ash
  var <- tau2 * Dsq + sigma2
  diagXtOmegaX <- rowSums(sweep(V^2, 2, Dsq / var, `*`))

  # Compute XtOmega for mr.ash
  # TODO: This will become a bottleneck for larger matrices. Switch to Rcpp.
  XtOmega <- V %*% sweep(VtXt, 1, 1/var, `*`)

  # Call mr.ash. Kept hard coded for now to show this is a modified version of the
  # original mr.ash package.
  mrash_output <- mr.ash.alpha.mccreight::mr.ash(
    X = X,
    y = y,
    sa2 = est_sa2,
    intercept = FALSE,
    standardize = FALSE,
    sigma2 = sigma2,
    update.sigma2 = update_ash_sigma,
    diagXtOmegaX = diagXtOmegaX,
    XtOmega = XtOmega,
    V = V,
    tausq = tau2,
    sum_Dsq = sum(Dsq),
    Dsq = Dsq,
    VtXt = VtXt
  )

  # Extract updated variance components
  sigma2_new <- mrash_output$sigma2
  tau2_new <- sum(est_sa2 * mrash_output$pi)

  return(list(
    sigma2 = sigma2_new,
    tau2 = tau2_new,
    theta = mrash_output$beta,
    pi = mrash_output$pi,
    est_sa2 = est_sa2
  ))
}

# Method of Moments variance estimation for unmappable effects methods
#' @keywords internal
mom_unmappable <- function(alpha, mu, omega, sigma2, tau2, n, V, Dsq, VtXty, Xty, yty,
                           est_sigma2, est_tau2, verbose) {
  L <- nrow(mu)
  p <- ncol(mu)

  A <- matrix(0, nrow = 2, ncol = 2)
  A[1, 1] <- n
  A[1, 2] <- sum(Dsq)
  A[2, 1] <- A[1, 2]
  A[2, 2] <- sum(Dsq^2)

  # Compute diag(V'MV)
  b <- colSums(mu * alpha)
  Vtb <- t(V) %*% b
  diagVtMV <- Vtb^2
  tmpD <- rep(0, p)

  for (l in seq_len(L)) {
    bl <- mu[l, ] * alpha[l, ]
    Vtbl <- t(V) %*% bl
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + alpha[l, ] * (mu[l, ]^2 + 1 / omega[l, ])
  }

  diagVtMV <- diagVtMV + rowSums(sweep(t(V)^2, 2, tmpD, `*`))

  # Compute x
  x <- rep(0, 2)
  x[1] <- yty - 2 * sum(b * Xty) + sum(Dsq * diagVtMV)
  x[2] <- sum(Xty^2) - 2 * sum(Vtb * VtXty * Dsq) + sum(Dsq^2 * diagVtMV)

  if (est_tau2) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigma2 <- sol[1]
      tau2 <- sol[2]
    } else {
      sigma2 <- x[1] / n
      tau2 <- 0
    }
    if (verbose) {
      cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigma2, tau2))
    }
  } else if (est_sigma2) {
    sigma2 <- (x[1] - A[1, 2] * tau2) / n
    if (verbose) {
      cat(sprintf("Update sigma^2 to %f\n", sigma2))
    }
  }
  return(list(sigma2 = sigma2, tau2 = tau2))
}

# Compute theta using BLUP
#' @keywords internal
compute_theta_blup <- function(data, model) {
  # Calculate diagXtOmegaX, diagonal variances, and Beta
  omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
  b <- colSums(model$mu * model$alpha)

  XtOmegaXb <- data$eigen_vectors %*% ((t(data$eigen_vectors) %*% b) * data$eigen_values / omega_res$omega_var)
  XtOmegay <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)
  XtOmegar <- XtOmegay - XtOmegaXb

  theta <- model$tau2 * XtOmegar

  return(theta)
}

# Compute Omega-weighted quantities for unmappable effects methods
#' @keywords internal
compute_omega_quantities <- function(data, tau2, sigma2) {
  # Compute variance in eigen space
  omega_var <- tau2 * data$eigen_values + sigma2

  # Compute diagonal of X'OmegaX
  diagXtOmegaX <- rowSums(sweep(data$eigen_vectors^2, 2, (data$eigen_values / omega_var), `*`))

  return(list(
    omega_var = omega_var,
    diagXtOmegaX = diagXtOmegaX
  ))
}

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

# Correlation function with specified warning muffled.
#' @importFrom stats cor
#' @keywords internal
muffled_corr <- function(x) {
  withCallingHandlers(cor(x),
    warning = function(w) {
      if (grepl("the standard deviation is zero", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# cov2cor function with specified warning muffled.
#' @importFrom stats cov2cor
#' @keywords internal
muffled_cov2cor <- function(x) {
  withCallingHandlers(cov2cor(x),
    warning = function(w) {
      if (grepl(
        "had 0 or NA entries; non-finite result is doubtful",
        w$message
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Check for symmetric matrix.
#' @keywords internal
is_symmetric_matrix <- function(x) {
  if (requireNamespace("Rfast", quietly = TRUE)) {
    return(Rfast::is.symmetric(x))
  } else {
    return(Matrix::isSymmetric(x))
  }
}

# Slim the result of fitted susie model.
#' @keywords internal
susie_slim <- function(res) {
  list(alpha = res$alpha, niter = res$niter, V = res$V, sigma2 = res$sigma2)
}

# Prune single effects to given number L in susie model object.
#' @keywords internal
susie_prune_single_effects <- function(s, L = 0, V = NULL) {
  num_effects <- nrow(s$alpha)
  if (L == 0) {
    # Filtering will be based on non-zero elements in s$V.
    if (!is.null(s$V)) {
      L <- length(which(s$V > 0))
    } else {
      L <- num_effects
    }
  }
  if (L == num_effects) {
    s$sets <- NULL
    return(s)
  }
  if (!is.null(s$sets$cs_index)) {
    effects_rank <- c(s$sets$cs_index, setdiff(1:num_effects, s$sets$cs_index))
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

    s$alpha <- rbind(
      s$alpha[effects_rank, ],
      matrix(1 / ncol(s$alpha), L - num_effects, ncol(s$alpha))
    )
    for (n in c("mu", "mu2", "lbf_variable")) {
      if (!is.null(s[[n]])) {
        s[[n]] <- rbind(
          s[[n]][effects_rank, ],
          matrix(0, L - num_effects, ncol(s[[n]]))
        )
      }
    }
    for (n in c("KL", "lbf")) {
      if (!is.null(s[[n]])) {
        s[[n]] <- c(s[[n]][effects_rank], rep(NA, L - num_effects))
      }
    }
    if (!is.null(V)) {
      if (length(V) > 1) {
        V[1:num_effects] <- s$V[effects_rank]
      } else {
        V <- rep(V, L)
      }
    }
    s$V <- V
  }
  s$sets <- NULL
  return(s)
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

# computes column standard deviations for any type of matrix
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

# Check whether b is in space spanned by the non-zero eigenvectors
# of A
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

# Apply operation f to all nonzeros of a sparse matrix.
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix summary
#' @keywords internal
apply_nonzeros <- function(X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i, j = d$j, x = f(d$x), dims = dim(X)))
}

# Validate Model Initialization Object
#' @keywords internal
validate_init <- function(model_init, L, null_weight) {
  # Check if model_init is a susie object
  if (!inherits(model_init, "susie")) {
    stop("model_init must be a 'susie' object")
  }

  alpha <- model_init$alpha
  mu <- model_init$mu
  mu2 <- model_init$mu2
  V <- model_init$V
  sigma2 <- model_init$sigma2
  pi_w <- model_init$pi
  null_id <- model_init$null_index


  # TODO: Fix this check
  # if(null_id > 0 && is.null(null_weight) || null_weight == 0)
  #   stop("There is a mistmatch in null_weight between the initalization object",
  #        " and the current call. Please make them consistent.")

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

  # Verify no NA/Inf values in V
  if (any(!is.finite(V))) {
    stop("model_init$V contains NA/Inf values")
  }

  # Verify no NA/Inf values in sigma2
  if (any(!is.finite(sigma2))) {
    stop("model_init$sigma2 contains NA/Inf")
  }

  # Verify no NA/Inf values in prior weights
  if (any(!is.finite(pi_w))) {
    stop("model_init$pi contains NA/Inf")
  }

  # Verify alpha is matrix
  if (!is.matrix(alpha)) {
    stop("model_init$alpha must be a matrix")
  }

  # Verify alpha values are between [0,1]
  if (max(model_init$alpha) > 1 || min(model_init$alpha) < 0) {
    stop(
      "model_init$alpha has invalid values outside range [0,1]; please ",
      "check your input"
    )
  }

  # Verify alpha dimensions for number of requested effects
  if (nrow(model_init$alpha) > L) {
    stop("model_init has more effects than requested L")
  }

  # Verify mu & mu2 dimensions match alpha
  if (!all(dim(mu) == dim(alpha))) {
    stop("model_init$mu and model_init$alpha dimensions do not match")
  }
  if (!all(dim(mu2) == dim(alpha))) {
    stop("model_init$mu2 and model_init$alpha dimensions do not match")
  }

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

  # Verify sigma2 is numeric and non-negative
  if (!is.numeric(sigma2)) {
    stop("model_init$sigma2 must be numeric")
  }
  if (sigma2 < 0) {
    stop("model_init$sigma2 is negative")
  }

  # Verify prior weight properties
  if (length(pi_w) != ncol(alpha)) {
    stop(
      "model_init$pi should have the same length as the number of columns",
      " in model_init$alpha"
    )
  }
  # TODO: fix this check. is this a floating point difference ? maybe set a tol ?
  # if (sum(pi_w) != 1)
  #   stop("model_init$pi must sum to one")

  invisible(model_init)
}

# Adjust the number of effects
#' @keywords internal
adjust_L <- function(model_init, L, V) {
  num_effects <- nrow(model_init$alpha)
  if (num_effects > L) {
    warning(paste0(
      "Requested L = ", L,
      " is smaller than the ", num_effects,
      " effects in model_init after pruning; ",
      "using L = ", num_effects, " instead."
    ))
    L <- num_effects
  }

  model_init <- susie_prune_single_effects(model_init, L = L, V = V)

  return(list(model_init = model_init, L = L))
}

# Initialize Null Index
#' @keywords internal
initialize_null_index <- function(null_weight, p) {
  if (is.null(null_weight) || null_weight == 0) {
    null_idx <- 0
  } else {
    null_idx <- p
  }
  return(null_idx)
}

# Helper function to assign variable names to model components
#' @keywords internal
assign_names <- function(model, variable_names, null_weight, p) {
  if (!is.null(variable_names)) {
    if (!is.null(null_weight)) {
      variable_names[length(variable_names)] <- "null"
      names(model$pip) <- variable_names[-p]
    } else {
      names(model$pip) <- variable_names
    }
    colnames(model$alpha) <- variable_names
    colnames(model$mu) <- variable_names
    colnames(model$mu2) <- variable_names
    colnames(model$lbf_variable) <- variable_names
  }
  return(model)
}

# Helper function to update variance components and derived quantities
#' @keywords internal
update_model_variance <- function(data, model, lowerbound, upperbound, estimate_method = "MLE") {
  variance_result <- update_variance_components(data, model, estimate_method)
  model$sigma2 <- max(lowerbound, variance_result$sigma2)
  model$sigma2 <- min(model$sigma2, upperbound)

  # Update additional variance components if they exist
  if (!is.null(variance_result$tau2)) {
    model$tau2 <- variance_result$tau2
  }

  # Handle unmappable effects outputs (theta for both inf and ash)
  if (!is.null(variance_result$theta)) {
    model$theta <- variance_result$theta
  }

  # Handle ash weights
  if (!is.null(variance_result$ash_pi)) {
    model$ash_pi <- variance_result$ash_pi
  }

  # Update derived quantities after variance component changes
  data <- update_derived_quantities(data, model)

  # Update fitted values to include theta if it exists
  if (!is.null(model$theta)) {
    b <- colSums(model$alpha * model$mu)
    model$XtXr <- data$XtX %*% (b + model$theta)
  }

  # Transfer RSS lambda specific updates from data to model
  if (!is.null(data$SinvRj_temp)) {
    model$SinvRj <- data$SinvRj_temp
    model$RjSinvRj <- data$RjSinvRj_temp
    data$SinvRj_temp <- NULL
    data$RjSinvRj_temp <- NULL
  }

  return(list(data = data, model = model))
}

# Get posterior inclusion probabilities
#' @keywords internal
get_pip <- function(data, model, coverage, min_abs_corr, prior_tol) {

}


# Objective function (ELBO)
#' @keywords internal
get_objective <- function(data, model, verbose = FALSE) {
  if (!is.null(data$unmappable_effects) && data$unmappable_effects == "inf") {
    # Compute omega
    L <- nrow(model$alpha)
    omega_res <- compute_omega_quantities(data, model$tau2, model$sigma2)
    omega <- matrix(0, L, data$p)
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
  } else {
    # Standard ELBO computation
    objective <- Eloglik(data, model) - sum(model$KL)
  }

  if (is.infinite(objective)) {
    stop("get_objective() produced an infinite ELBO value")
  }
  if (verbose) {
    print(paste0("objective:", objective))
  }
  return(objective)
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

# Initialize core susie model object with default parameter matrices
#' @keywords internal
initialize_matrices <- function(p, L, scaled_prior_variance, var_y, residual_variance,
                                prior_weights, include_unmappable = FALSE) {
  mat_init <- list(
    alpha = matrix(1 / p, L, p),
    mu = matrix(0, L, p),
    mu2 = matrix(0, L, p),
    V = rep(scaled_prior_variance * var_y, L),
    KL = rep(as.numeric(NA), L),
    lbf = rep(as.numeric(NA), L),
    lbf_variable = matrix(as.numeric(NA), L, p),
    sigma2 = residual_variance,
    pi = prior_weights
  )

  # Add unmappable effects specific components
  if (include_unmappable) {
    mat_init$tau2 <- 0
    mat_init$theta <- rep(0, p)
  }

  return(mat_init)
}

# MLE variance estimation for unmappable effects
#' @keywords internal
mle_unmappable <- function(alpha, mu, omega, sigma2, tau2, n,
                           eigen_vectors, eigen_values, VtXty, yty,
                           est_sigma2 = TRUE, est_tau2 = TRUE,
                           sigma2_range = NULL, tau2_range = NULL,
                           verbose = FALSE) {
  L <- nrow(mu)
  p <- ncol(mu)

  # Set default ranges if not provided
  if (is.null(sigma2_range)) {
    sigma2_range <- c(0.2 * yty / n, 1.2 * yty / n)
  }
  if (is.null(tau2_range)) {
    tau2_range <- c(1e-12, 1.2 * yty / (n * p))
  }

  # Compute diag(V'MV)
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

  # Negative ELBO as function of x = (sigma^2, tau^2)
  f <- function(x) {
    sigma2_val <- x[1]
    tau2_val <- x[2]
    var_val <- tau2_val * eigen_values + sigma2_val

    0.5 * (n - p) * log(sigma2_val) + 0.5 / sigma2_val * yty +
      sum(0.5 * log(var_val) -
        0.5 * tau2_val / sigma2_val * VtXty^2 / var_val -
        Vtb * VtXty / var_val +
        0.5 * eigen_values / var_val * diagVtMV)
  }

  # Negative ELBO for sigma^2 only (when tau^2 is fixed)
  g <- function(sigma2_val) {
    f(c(sigma2_val, tau2))
  }

  if (est_tau2) {
    # Optimize both sigma^2 and tau^2
    res <- optim(
      par = c(sigma2, tau2),
      fn = f,
      method = "L-BFGS-B",
      lower = c(sigma2_range[1], tau2_range[1]),
      upper = c(sigma2_range[2], tau2_range[2])
    )

    if (res$convergence == 0) {
      sigma2 <- res$par[1]
      tau2 <- res$par[2]
      if (verbose) {
        cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigma2, tau2))
      }
    } else {
      warning("MLE optimization failed to converge; keeping previous parameters")
    }
  } else if (est_sigma2) {
    # Optimize only sigma^2
    res <- optim(
      par = sigma2,
      fn = g,
      method = "L-BFGS-B",
      lower = sigma2_range[1],
      upper = sigma2_range[2]
    )

    if (res$convergence == 0) {
      sigma2 <- res$par
      if (verbose) {
        cat(sprintf("Update sigma^2 to %f\n", sigma2))
      }
    } else {
      warning("MLE optimization failed to converge; keeping previous parameters")
    }
  }

  return(list(sigma2 = sigma2, tau2 = tau2))
}

# Compute ELBO for infinitesimal effects model
#' @keywords internal
compute_elbo_inf <- function(alpha, mu, omega, lbf, sigma2, tau2, n, p,
                             eigen_vectors, eigen_values, VtXty, yty) {
  L <- nrow(mu)

  # Compute b = sum_l alpha_l * mu_l
  b <- colSums(mu * alpha)
  Vtb <- t(eigen_vectors) %*% b

  # Compute E[theta'theta] diagonal
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

  # Compute negative ELBO (matching original formulation)
  neg_elbo <- 0.5 * (n - p) * log(sigma2) + 0.5 / sigma2 * yty +
    sum(0.5 * log(var) -
      0.5 * tau2 / sigma2 * VtXty^2 / var -
      Vtb * VtXty / var +
      0.5 * eigen_values / var * diagVtMV)

  # Total ELBO includes both infinitesimal and sparse components
  elbo_infinitesimal <- -neg_elbo
  elbo_sparse <- sum(lbf)
  elbo_total <- elbo_infinitesimal + elbo_sparse

  return(elbo_total)
}

# @title sets the attributes for the R matrix
# @param R a p by p LD matrix
# @param r_tol tolerance level for eigen value check of positive
#   semidefinite matrix of R.
# @return R with attribute e.g., attr(R, 'eigenR') is the eigen
#   decomposition of R.
set_R_attributes <- function(R, r_tol) {
  if (is.null(attr(R, "eigen"))) {
    eigenR <- eigen(R, symmetric = TRUE)
  } else {
    eigenR <- attr(R, "eigen")
  }

  # Drop small eigenvalues.
  eigenR$values[abs(eigenR$values) < r_tol] <- 0
  if (any(eigenR$values < 0)) {
    min_lambda <- min(eigenR$values)
    eigenR$values[eigenR$values < 0] <- 0
    warning_message(paste0(
      "The input correlation matrix has negative eigenvalues ",
      "(smallest one is ", min_lambda, "). The correlation ",
      "matrix is adjusted such that these negative eigenvalues ",
      "are now zeros. You can ignore this message, only if you ",
      "believe the negative eigenvalue is result of numerical ",
      "rounding errors."
    ))
  }
  res <- eigenR$vectors %*% (t(eigenR$vectors) * eigenR$values)

  attr(res, "eigen") <- eigenR
  attr(res, "d") <- diag(res)
  attr(res, "scaled:scale") <- rep(1, length = nrow(R))
  return(res)
}

remove_null_effects <- function(s) {
  null_indices <- (s$V == 0)
  s$alpha <- s$alpha[!null_indices, , drop = FALSE]
  s$mu <- s$mu[!null_indices, , drop = FALSE]
  s$mu2 <- s$mu2[!null_indices, , drop = FALSE]
  s$lbf_variable <- s$lbf_variable[!null_indices, , drop = FALSE]
  s$V <- s$V[!null_indices, drop = FALSE]
  return(s)
}

add_null_effect <- function(s, V) {
  p <- ncol(s$alpha)
  s$alpha <- rbind(s$alpha, 1 / p)
  s$mu <- rbind(s$mu, rep(0, p))
  s$mu2 <- rbind(s$mu2, rep(0, p))
  s$lbf_variable <- rbind(s$lbf_variable, rep(0, p))
  s$V <- c(s$V, V)
  return(s)
}

# Servin and Stephens prior helper functions
#' @keywords internal

# Compute log Bayes factor for Servin and Stephens prior
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
posterior_mean_servin_stephens <- function(xtx, xty, s0_t = 1) {
  omega <- (xtx + (1 / s0_t^2))^(-1)
  b_bar <- omega %*% xty
  return(b_bar)
}

# Posterior variance for Servin and Stephens prior using sufficient statistics
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

# Convert individual data to ss with unmappable effects components.
#' @keywords internal
convert_individual_to_ss_unmappable <- function(data) {
  # Compute sufficient statistics
  XtX <- crossprod(data$X)
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
      y_mean = data$mean_y,
      prior_weights = data$prior_weights,
      null_weight = data$null_weight,
      unmappable_effects = data$unmappable_effects,
      convergence_method = data$convergence_method,
      use_servin_stephens = FALSE
    ),
    class = "ss"
  )

  # Copy attributes from X to XtX
  attr(ss_data$XtX, "d") <- attr(data$X, "d")
  attr(ss_data$XtX, "scaled:scale") <- attr(data$X, "scaled:scale")

  # Add eigen decomposition for unmappable effects methods
  ss_data <- add_eigen_decomposition(ss_data, data)

  return(ss_data)
}

# Check convergence
#' @keywords internal
check_convergence <- function(model_prev, model_current, elbo, tol, convergence_method, iter) {
  # Skip convergence check on first iteration
  if(iter == 1) {
    return(FALSE)
  }

  ELBO_diff <- elbo[iter + 1] - elbo[iter]
  PIP_diff <- max(abs(model_prev$alpha - model_current$alpha))

  # If ELBO calculation produces NA/Inf value, fallback to PIP-based convergence
  if((is.na(ELBO_diff) || is.infinite(ELBO_diff)) && convergence_method == "elbo"){
    warning(paste0("Iteration ", iter, " produced an NA/infinite ELBO value. Using pip-based convergence this iteration."))
    convergence_method <- "pip"
  }

  if (convergence_method == "pip") {
    return(PIP_diff < tol)
  } else {
    return(ELBO_diff < tol)
  }
}

# Stabilize log Bayes factors and compute log posterior odds
#' @keywords internal
lbf_stabilization <- function(lbf, prior_weights, shat2 = NULL) {

  # Add numerical stability to prior weights
  lpo <- lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Handle special case of infinite shat2 (e.g., when variable doesn't vary)
  infinite_idx <- is.infinite(shat2)
  lbf[infinite_idx] <- 0
  lpo[infinite_idx] <- 0

  return(list(lbf = lbf, lpo = lpo))
}

# Compute alpha and lbf for each effect
#' @keywords internal
compute_posterior_weights <- function(lpo) {

  w_weighted <- exp(lpo - max(lpo))
  weighted_sum_w <- sum(w_weighted)
  alpha <- w_weighted / weighted_sum_w

  return(list(
    alpha = alpha,
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

# Add eigen decomposition to ss objects (for unmappable effects methods)
#' @keywords internal
add_eigen_decomposition <- function(data, individual_data) {
  # Compute eigen decomposition of correlation matrix
  eigen_decomp <- compute_eigen_decomposition(data$XtX, data$n)

  # Add eigen components to data object
  data$eigen_vectors <- eigen_decomp$V
  data$eigen_values <- eigen_decomp$Dsq
  data$VtXty <- t(eigen_decomp$V) %*% data$Xty

  # SuSiE.ash requires the original X matrix, y vector, and VtXt
  if (data$unmappable_effects == "ash") {
    data$X <- individual_data$X
    data$y <- individual_data$y
    data$VtXt <- t(data$eigen_vectors) %*% t(individual_data$X)
  }

  # Precompute diagXtOmegaX and XtOmegay using initial values of sigma2 and tau2
  omega_res <- compute_omega_quantities(data, tau2 = 0, sigma2 = 1)
  data$omega_var <- omega_res$omega_var
  data$diagXtOmegaX <- omega_res$diagXtOmegaX
  data$XtOmegay <- data$eigen_vectors %*% (data$VtXty / omega_res$omega_var)

  return(data)
}
