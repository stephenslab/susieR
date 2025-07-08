# Individual-level data constructor
susie_constructor <- function(X, y,
                              intercept = TRUE,
                              standardize = TRUE,
                              na.rm = FALSE,
                              prior_weights = NULL,
                              null_weight = 0,
                              non_sparse_method = "none") {
  # Validate input X
  if (!(is.double(X) & is.matrix(X)) &
      !inherits(X, "CsparseMatrix") &
      is.null(attr(X, "matrix.type")))
    stop("Input X must be a double-precision matrix, or a sparse matrix, or ",
         "a trend filtering matrix")

  if (anyNA(X))
    stop("X contains NA values")

  # Handle missing values in y
  if (anyNA(y)) {
    if (na.rm) {
      samples_kept <- which(!is.na(y))
      y <- y[samples_kept]
      X <- X[samples_kept, ]
    } else {
      stop("Input y must not contain missing values")
    }
  }
  mean_y <- mean(y)

  # Handle null weights
  if (is.numeric(null_weight) && null_weight == 0)
    null_weight <- NULL

  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    
    if (is.null(prior_weights))
      prior_weights <- c(rep(1/ncol(X) * (1 - null_weight), ncol(X)), null_weight)
    else
      prior_weights <- c(prior_weights * (1 - null_weight), null_weight)

    # add the extra 0 column to X
    X <- cbind(X, 0)
  }

  # Store dimensions
  n <- nrow(X)
  p <- ncol(X)
  if (p > 1000 & !requireNamespace("Rfast", quietly = TRUE))
    warning_message("For an X with many columns, please consider installing ",
                    "the Rfast package for more efficient credible set (CS) ",
                    "calculations.", style = "hint")

  # Center y if intercept is included
  if (intercept)
    y <- y - mean_y

  # Compute and set X matrix attributes
  out <- compute_colstats(X, center = intercept, scale = standardize)
  attr(X, "scaled:center") <- out$cm
  attr(X, "scaled:scale") <- out$csd
  attr(X, "d") <- out$d

  data_object <- structure(list(
    X = X,
    y = y,
    mean_y = mean_y,
    n = n,
    p = p,
    prior_weights = prior_weights,
    null_weight = null_weight),
    class = "individual")

  # Configure data object for specified non-sparse method
  data_object <- configure_data(data_object, non_sparse_method)

  return(data_object)
}

# Sufficient statistics data constructor
susie_ss_constructor <- function(XtX, Xty, yty, n,
                                 X_colmeans = NA, y_mean = NA, maf = NULL,
                                 maf_thresh = 0, standardize = TRUE,
                                 r_tol = 1e-8, check_input = FALSE,
                                 prior_weights = NULL, null_weight = 0,
                                 non_sparse_method = "none") {

  # Validate required inputs
  if (missing(n))
    stop("n must be provided")
  if (n <= 1)
    stop("n must be greater than 1")
  if (missing(XtX) || missing(Xty) || missing(yty))
    stop("XtX, Xty, yty must all be provided")

  if (!(is.double(XtX) && is.matrix(XtX)) &&
      !inherits(XtX, "CsparseMatrix"))
    stop("XtX must be a numeric dense or sparse matrix")

  if (ncol(XtX) != length(Xty))
    stop(paste0("The dimension of XtX (",nrow(XtX)," by ",ncol(XtX),
                ") does not agree with expected (",length(Xty)," by ",
                length(Xty),")"))

  if (ncol(XtX) > 1000 & !requireNamespace("Rfast", quietly = TRUE))
    warning_message("For large R or large XtX, consider installing the ",
                    "Rfast package for better performance.", style = "hint")

  # Ensure XtX is symmetric
  if (!is_symmetric_matrix(XtX)) {
    warning("XtX not symmetric; using (XtX + t(XtX))/2")
    XtX <- (XtX + t(XtX))/2
  }

  # Apply MAF filter if provided
  if (!is.null(maf)) {
    if (length(maf) != length(Xty))
      stop(paste("The length of maf does not agree with expected", length(Xty)))
    id <- which(maf > maf_thresh)
    XtX <- XtX[id, id]
    Xty <- Xty[id]
  }

  # Additional validation
  if (any(is.infinite(Xty)))
    stop("Input Xty contains infinite values")
  if (!(is.double(XtX) & is.matrix(XtX)) & !inherits(XtX, "CsparseMatrix"))
    stop("Input XtX must be a double-precision matrix, or a sparse matrix")
  if (anyNA(XtX))
    stop("Input XtX matrix contains NAs")

  # Replace NAs in Xty with zeros
  if (anyNA(Xty)) {
    warning_message("NA values in Xty are replaced with 0")
    Xty[is.na(Xty)] <- 0
  }

  # Positive-semidefinite check
  if (check_input) {
    semi_pd <- check_semi_pd(XtX, r_tol)
    if (!semi_pd$status)
      stop("XtX is not a positive semidefinite matrix")

    # Check whether Xty lies in space spanned by non-zero eigenvectors of XtX
    proj <- check_projection(semi_pd$matrix, Xty)
    if (!proj$status)
      warning_message("Xty does not lie in the space of the non-zero eigenvectors ",
                      "of XtX")
  }

  # Handle null weights
  if (is.numeric(null_weight) && null_weight == 0)
    null_weight <- NULL

  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (is.null(prior_weights))
      prior_weights <- c(rep(1/ncol(XtX)*(1-null_weight), ncol(XtX)), null_weight)
    else
      prior_weights <- c(prior_weights*(1 - null_weight), null_weight)
    XtX <- cbind(rbind(XtX, 0), 0)
    Xty <- c(Xty, 0)
    if (length(X_colmeans) == 1)
      X_colmeans <- rep(X_colmeans, p)
    if (length(X_colmeans) != p)
      stop("The length of X_colmeans does not agree with number of variables")
  }

  # Define p
  p <- ncol(XtX)

  # Standardize if requested
  if (standardize) {
    dXtX <- diag(XtX)
    csd <- sqrt(dXtX/(n-1))
    csd[csd == 0] <- 1
    XtX <- t((1/csd) * XtX) / csd
    Xty <- Xty / csd
  } else {
    csd <- rep(1, length = p)
  }

  attr(XtX, "d") <- diag(XtX)
  attr(XtX, "scaled:scale") <- csd

  if (length(X_colmeans) == 1)
    X_colmeans <- rep(X_colmeans, p)
  if (length(X_colmeans) != p)
    stop("`X_colmeans` length (", length(X_colmeans),
         ") does not match number of variables (", p, ")")

  # Assemble data object
  data_object <- structure(list(
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    p = p,
    X_colmeans = X_colmeans,
    y_mean = y_mean,
    prior_weights = prior_weights,
    null_weight = null_weight),
    class = "ss")

  # Configure data object for specified non-sparse method
  data_object <- configure_data(data_object, non_sparse_method)

  return(data_object)
}
