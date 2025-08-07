# internal utility functions

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

# Method of Moments variance estimation for unmappable effects methods
#' @keywords internal
MoM <- function(alpha, mu, omega, sigma2, tau2, n, V, Dsq, VtXty, Xty, yty,
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
  alpha <- model$alpha
  mu <- model$mu
  sigma2 <- model$sigma2
  tau2 <- if (is.null(model$tau2)) 0 else model$tau2

  b <- colSums(mu * alpha)

  var <- tau2 * data$eigen_values + sigma2
  XtOmegaXb <- data$eigen_vectors %*% ((t(data$eigen_vectors) %*% b) * data$eigen_values / var)

  XtOmegar <- data$XtOmegay - XtOmegaXb

  theta <- tau2 * XtOmegar

  return(theta)
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
#' @importFrom stats median
#' @keywords internal
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
    if (length(pos) > n) {
      pos = sample(pos,n)
    }

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
#' @importFrom stats cor
#' @keywords internal
muffled_corr = function (x)
  withCallingHandlers(cor(x),
                      warning = function(w) {
                        if (grepl("the standard deviation is zero",w$message))
                          invokeRestart("muffleWarning")
                      })

# cov2cor function with specified warning muffled.
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
    return(Matrix::isSymmetric(x))
}

# Compute standard error for regression coef.
# S = (X'X)^-1 \Sigma
#' @keywords internal
calc_stderr = function (X, residuals)
  sqrt(diag(sum(residuals^2)/(nrow(X) - 2) * chol2inv(chol(crossprod(X)))))

# Return residuals of Y after removing the linear effects of the susie
# model.
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

# Compute the column means of X, the column standard deviations of X,
# and rowSums(Y^2), where Y is the centered and/or scaled version of
# X.
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @keywords internal
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

# computes column standard deviations for any type of matrix
# This should give the same result as matrixStats::colSds(X),
# but allows for sparse matrices as well as dense ones.
#' @importFrom matrixStats colSds
#' @importFrom Matrix summary
#' @keywords internal
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

# Check whether A is positive semidefinite
#' @keywords internal
check_semi_pd = function (A, tol) {
  attr(A,"eigen") = eigen(A,symmetric = TRUE)
  v = attr(A,"eigen")$values
  v[abs(v) < tol] = 0
  return(list(matrix      = A,
              status      = !any(v < 0),
              eigenvalues = v))
}

# Check whether b is in space spanned by the non-zero eigenvectors
# of A
#' @keywords internal
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

# Utility function to display warning messages as they occur
#' @importFrom crayon combine_styles
#' @keywords internal
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
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix summary
#' @keywords internal
apply_nonzeros <- function (X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i,j = d$j,x = f(d$x),dims = dim(X)))
}

# Validate Model Initialization Object
#' @keywords internal
validate_init <- function(model_init, L, null_weight){

  # Check if model_init is a susie object
  if (!inherits(model_init, "susie"))
    stop("model_init must be a 'susie' object")

  alpha   <- model_init$alpha
  mu      <- model_init$mu
  mu2     <- model_init$mu2
  V       <- model_init$V
  sigma2  <- model_init$sigma2
  pi_w    <- model_init$pi
  null_id <- model_init$null_index


  # TODO: Fix this check
  # if(null_id > 0 && is.null(null_weight) || null_weight == 0)
  #   stop("There is a mistmatch in null_weight between the initalization object",
  #        " and the current call. Please make them consistent.")

  # Verify no NA/Inf values in alpha
  if (any(!is.finite(alpha)))
    stop("model_init$alpha contains NA/Inf values")

  # Verify no NA/Inf values in mu
  if (any(!is.finite(mu)))
    stop("model_init$mu contains NA/Inf values")

  # Verify no NA/Inf values in mu2
  if (any(!is.finite(mu2)))
    stop("model_init$mu2 contains NA/Inf values")

  # Verify no NA/Inf values in V
  if (any(!is.finite(V)))
    stop("model_init$V contains NA/Inf values")

  # Verify no NA/Inf values in sigma2
  if (any(!is.finite(sigma2)))
    stop("model_init$sigma2 contains NA/Inf")

  # Verify no NA/Inf values in prior weights
  if (any(!is.finite(pi_w)))
    stop("model_init$pi contains NA/Inf")

  # Verify alpha is matrix
  if (!is.matrix(alpha))
    stop("model_init$alpha must be a matrix")

  # Verify alpha values are between [0,1]
  if (max(model_init$alpha) > 1 || min(model_init$alpha) < 0)
    stop("model_init$alpha has invalid values outside range [0,1]; please ",
         "check your input")

  # Verify alpha dimensions for number of requested effects
  if (nrow(model_init$alpha) > L)
    stop("model_init has more effects than requested L")

  # Verify mu & mu2 dimensions match alpha
  if (!all(dim(mu)  == dim(alpha)))
    stop("model_init$mu and model_init$alpha dimensions do not match")
  if (!all(dim(mu2) == dim(alpha)))
    stop("model_init$mu2 and model_init$alpha dimensions do not match")

  # Verify V & alpha dimensions agree
  if (length(V) != nrow(alpha))
    stop("length(model_init$V) (", length(V), ") does not equal nrow(model_init$alpha) (",
         nrow(alpha), ")")

  # Verify V is numeric and non-negative
  if(!is.numeric(V))
    stop("model_init$V must be numeric")
  if(any(V < 0))
    stop("model_init$V has at least one negative value")

  # Verify sigma2 is numeric and non-negative
  if(!is.numeric(sigma2))
    stop("model_init$sigma2 must be numeric")
  if(sigma2 < 0)
    stop("model_init$sigma2 is negative")

  # Verify prior weight properties
  if(length(pi_w) != ncol(alpha))
    stop("model_init$pi should have the same length as the number of columns",
         " in model_init$alpha")
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
    warning(paste0("Requested L = ", L,
                   " is smaller than the ", num_effects,
                   " effects in model_init after pruning; ",
                   "using L = ", num_effects, " instead."))
    L <- num_effects
  }

  model_init <- susie_prune_single_effects(model_init, L = L, V = V)

  return(list(model_init = model_init, L = L))
}

# Initialize Null Index
#' @keywords internal
initialize_null_index <- function(null_weight, p){
  if(is.null(null_weight) || null_weight == 0){
    null_idx = 0
  }else{
    null_idx = p
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
update_model_variance <- function(data, model, lowerbound, upperbound) {
  variance_result <- update_variance_components(data, model)
  model$sigma2 <- max(lowerbound, variance_result$sigma2)
  model$sigma2 <- min(model$sigma2, upperbound)

  # Update additional variance components if they exist
  if (!is.null(variance_result$tau2)) {
    model$tau2 <- variance_result$tau2
  }

  # Update derived quantities after variance component changes
  data <- update_derived_quantities(data, model)

  # Transfer theta from data to model if computed (for unmappable effects methods)
  if (!is.null(data$theta)) {
    model$theta <- data$theta

    # Update fitted values to include theta: XtXr = XtX %*% (b + theta)
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
  if (is.null(coverage) || is.null(min_abs_corr)) return(NULL)
  return(susie_get_pip(model, prune_by_cs = FALSE, prior_tol = prior_tol))
}


# Objective function (ELBO)
#' @keywords internal
get_objective <- function(data, model, verbose = FALSE) {
  objective <- Eloglik(data, model) - sum(model$KL)
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
  if(resid_var < 0) {
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

