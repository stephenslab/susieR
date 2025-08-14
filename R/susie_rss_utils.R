#' @title Compute sufficient statistics for input to \code{susie_suff_stat}
#'
#' @description Computes the sufficient statistics \eqn{X'X, X'y, y'y}
#'   and \eqn{n} after centering (and possibly standardizing) the
#'   columns of \eqn{X} and centering \eqn{y} to have mean zero. We also
#'   store the column means of \eqn{X} and mean of \eqn{y}.
#'
#' @param X An n by p matrix of covariates.
#'
#' @param y An n vector.
#'
#' @param standardize Logical flag indicating whether to standardize
#'   columns of X to unit variance prior to computing summary data
#'
#' @return A list of sufficient statistics (\code{XtX, Xty, yty, n})
#'   and \code{X_colmeans}, \code{y_mean}.
#'
#' @importFrom methods as
#' @importFrom Matrix colMeans
#' @importFrom Matrix crossprod
#'
#' @examples
#' data(N2finemapping)
#' ss <- compute_suff_stat(N2finemapping$X, N2finemapping$Y[, 1])
#'
#' @export
#'
compute_suff_stat <- function(X, y, standardize = FALSE) {
  y_mean <- mean(y)
  y <- y - y_mean
  n <- nrow(X)
  mu <- colMeans(X)
  s <- compute_colSds(X)
  Xty <- drop(y %*% X)
  XtX <- crossprod(X)
  XtX <- as.matrix(XtX)
  XtX <- XtX - n * tcrossprod(mu)
  if (standardize) {
    XtX <- XtX / s
    XtX <- t(XtX)
    XtX <- XtX / s
    Xty <- Xty / s
  }
  n <- length(y)
  yty <- sum(y^2)
  return(list(
    XtX = XtX, Xty = Xty, yty = yty, n = n,
    y_mean = y_mean, X_colmeans = mu
  ))
}

#' @title Estimate s in \code{susie_rss} Model Using Regularized LD
#'
#' @description The estimated s gives information about the
#'   consistency between the z scores and LD matrix. A larger \eqn{s}
#'   means there is a strong inconsistency between z scores and LD
#'   matrix. The \dQuote{null-mle} method obtains mle of \eqn{s} under
#'   \eqn{z | R ~ N(0,(1-s)R + s I)}, \eqn{0 < s < 1}. The
#'   \dQuote{null-partialmle} method obtains mle of \eqn{s} under
#'   \eqn{U^T z | R ~ N(0,s I)}, in which \eqn{U} is a matrix containing
#'   the of eigenvectors that span the null space of R; that is, the
#'   eigenvectors corresponding to zero eigenvalues of R. The estimated
#'   \eqn{s} from \dQuote{null-partialmle} could be greater than 1. The
#'   \dQuote{null-pseudomle} method obtains mle of \eqn{s} under
#'   pseudolikelihood \eqn{L(s) = \prod_{j=1}^{p} p(z_j | z_{-j}, s,
#'   R)}, \eqn{0 < s < 1}.
#'
#' @param z A p-vector of z scores.
#'
#' @param R A p by p symmetric, positive semidefinite correlation
#'   matrix.
#'
#' @param n The sample size. (Optional, but highly recommended.)
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param method a string specifies the method to estimate \eqn{s}.
#'
#' @return A number between 0 and 1.
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' p <- 1000
#' beta <- rep(0, p)
#' beta[1:4] <- 0.01
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' X <- scale(X, center = TRUE, scale = TRUE)
#' y <- drop(X %*% beta + rnorm(n))
#' input_ss <- compute_suff_stat(X, y, standardize = TRUE)
#' ss <- univariate_regression(X, y)
#' R <- cor(X)
#' attr(R, "eigen") <- eigen(R, symmetric = TRUE)
#' zhat <- with(ss, betahat / sebetahat)
#'
#' # Estimate s using the unadjusted z-scores.
#' s0 <- estimate_s_rss(zhat, R)
#'
#' # Estimate s using the adjusted z-scores.
#' s1 <- estimate_s_rss(zhat, R, n)
#'
#' @importFrom stats dnorm
#' @importFrom stats optim
#'
#' @export
#'
estimate_s_rss <- function(z, R, n, r_tol = 1e-08, method = "null-mle") {
  # Check and process input arguments z, R.
  z[is.na(z)] <- 0
  if (is.null(attr(R, "eigen"))) {
    attr(R, "eigen") <- eigen(R, symmetric = TRUE)
  }
  eigenld <- attr(R, "eigen")
  if (any(eigenld$values < -r_tol)) {
    warning_message(
      "The matrix R is not positive semidefinite. Negative ",
      "eigenvalues are set to zero"
    )
  }
  eigenld$values[eigenld$values < r_tol] <- 0

  # Check input n, and adjust the z-scores if n is provided.
  if (missing(n)) {
    warning_message(
      "Providing the sample size (n), or even a rough estimate of n, ",
      "is highly recommended. Without n, the implicit assumption is ",
      "n is large (Inf) and the effect sizes are small (close to zero)."
    )
  } else if (n <= 1) {
    stop("n must be greater than 1")
  }
  if (!missing(n)) {
    sigma2 <- (n - 1) / (z^2 + n - 2)
    z <- sqrt(sigma2) * z
  }

  if (method == "null-mle") {
    negloglikelihood <- function(s, ztv, d) {
      0.5 * sum(log((1 - s) * d + s)) +
        0.5 * tcrossprod(ztv / ((1 - s) * d + s), ztv)
    }
    s <- optim(0.5,
      fn = negloglikelihood, ztv = crossprod(z, eigenld$vectors),
      d = eigenld$values, method = "Brent", lower = 0, upper = 1
    )$par
  } else if (method == "null-partialmle") {
    colspace <- which(eigenld$values > 0)
    if (length(colspace) == length(z)) {
      s <- 0
    } else {
      znull <- crossprod(eigenld$vectors[, -colspace], z) # U2^T z
      s <- sum(znull^2) / length(znull)
    }
  } else if (method == "null-pseudomle") {
    pseudolikelihood <- function(s, z, eigenld) {
      precision <- eigenld$vectors %*% (t(eigenld$vectors) *
        (1 / ((1 - s) * eigenld$values + s)))
      postmean <- rep(0, length(z))
      postvar <- rep(0, length(z))
      for (i in 1:length(z)) {
        postmean[i] <- -(1 / precision[i, i]) * precision[i, -i] %*% z[-i]
        postvar[i] <- 1 / precision[i, i]
      }
      return(-sum(dnorm(z, mean = postmean, sd = sqrt(postvar), log = TRUE)))
    }
    s <- optim(0.5,
      fn = pseudolikelihood, z = z, eigenld = eigenld,
      method = "Brent", lower = 0, upper = 1
    )$par
  } else {
    stop("The method is not implemented")
  }
  return(s)
}

#' @title Compute Distribution of z-scores of Variant j Given Other z-scores, and Detect Possible Allele Switch Issue
#'
#' @description Under the null, the rss model with regularized LD
#'   matrix is \eqn{z|R,s ~ N(0, (1-s)R + s I))}. We use a mixture of
#'   normals to model the conditional distribution of z_j given other z
#'   scores, \eqn{z_j | z_{-j}, R, s ~ \sum_{k=1}^{K} \pi_k
#'   N(-\Omega_{j,-j} z_{-j}/\Omega_{jj}, \sigma_{k}^2/\Omega_{jj})},
#'   \eqn{\Omega = ((1-s)R + sI)^{-1}}, \eqn{\sigma_1, ..., \sigma_k}
#'   is a grid of fixed positive numbers. We estimate the mixture
#'   weights \eqn{\pi}  We detect the possible allele switch issue
#'   using likelihood ratio for each variant.
#'
#' @param z A p-vector of z scores.
#'
#' @param R A p by p symmetric, positive semidefinite correlation
#'   matrix.
#'
#' @param n The sample size. (Optional, but highly recommended.)
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param s an estimated s from \code{estimate_s_rss}
#'
#' @return a list containing a ggplot2 plot object and a table. The plot
#'   compares observed z score vs the expected value. The possible allele
#'   switched variants are labeled as red points (log LR > 2 and abs(z) > 2).
#'   The table summarizes the conditional distribution for each variant
#'   and the likelihood ratio test. The table has the following columns:
#'   the observed z scores, the conditional expectation, the conditional
#'   variance, the standardized differences between the observed z score
#'   and expected value, the log likelihood ratio statistics.
#'
#' @importFrom stats dnorm
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes_string
#' @importFrom mixsqp mixsqp
#'
#' @examples
#' # See also the vignette, "Diagnostic for fine-mapping with summary
#' # statistics."
#' set.seed(1)
#' n <- 500
#' p <- 1000
#' beta <- rep(0, p)
#' beta[1:4] <- 0.01
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' X <- scale(X, center = TRUE, scale = TRUE)
#' y <- drop(X %*% beta + rnorm(n))
#' ss <- univariate_regression(X, y)
#' R <- cor(X)
#' attr(R, "eigen") <- eigen(R, symmetric = TRUE)
#' zhat <- with(ss, betahat / sebetahat)
#' cond_dist <- kriging_rss(zhat, R, n = n)
#' cond_dist$plot
#'
#' @export
#'
kriging_rss <- function(z, R, n, r_tol = 1e-08,
                        s = estimate_s_rss(z, R, n, r_tol, method = "null-mle")) {
  # Check and process input arguments z, R.
  z[is.na(z)] <- 0
  if (is.null(attr(R, "eigen"))) {
    attr(R, "eigen") <- eigen(R, symmetric = TRUE)
  }
  eigenld <- attr(R, "eigen")
  if (any(eigenld$values < -r_tol)) {
    warning_message(
      "The matrix R is not positive semidefinite. Negative ",
      "eigenvalues are set to zero."
    )
  }
  eigenld$values[eigenld$values < r_tol] <- 0

  # Check and progress input argument s.
  force(s)
  if (s > 1) {
    warning_message("The given s is greater than 1. We replace it with 0.8.")
    s <- 0.8
  } else if (s < 0) {
    stop("The s must be non-negative")
  }

  # Check input n, and adjust the z-scores if n is provided.
  if ((!missing(n)) && (n <= 1)) {
    stop("n must be greater than 1")
  }
  if (missing(n)) {
    warning_message(
      "Providing the sample size (n), or even a rough estimate of n, ",
      "is highly recommended. Without n, the implicit assumption is ",
      "n is large (Inf) and the effect sizes are small (close to zero)."
    )
  } else {
    sigma2 <- (n - 1) / (z^2 + n - 2)
    z <- sqrt(sigma2) * z
  }

  dinv <- 1 / ((1 - s) * eigenld$values + s)
  dinv[is.infinite(dinv)] <- 0
  precision <- eigenld$vectors %*% (t(eigenld$vectors) * dinv)
  condmean <- rep(0, length(z))
  condvar <- rep(0, length(z))
  for (i in 1:length(z)) {
    condmean[i] <- -(1 / precision[i, i]) * precision[i, -i] %*% z[-i]
    condvar[i] <- 1 / precision[i, i]
  }
  z_std_diff <- (z - condmean) / sqrt(condvar)

  # obtain grid
  a_min <- 0.8
  if (max(z_std_diff^2) < 1) {
    a_max <- 2
  } else {
    a_max <- 2 * sqrt(max(z_std_diff^2))
  }
  npoint <- ceiling(log2(a_max / a_min) / log2(1.05))
  a_grid <- 1.05^(seq(-npoint, 0)) * a_max

  # compute likelihood
  sd_mtx <- outer(sqrt(condvar), a_grid)
  matrix_llik <- dnorm(z - condmean, sd = sd_mtx, log = TRUE)
  lfactors <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors

  # estimate weight
  w <- mixsqp(matrix_llik, log = TRUE, control = list(verbose = FALSE))$x

  # Compute denominators in likelihood ratios.
  logl0mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors

  # Compute numerators in likelihood ratios.
  matrix_llik <- dnorm(z + condmean, sd = sd_mtx, log = TRUE)
  lfactors <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  logl1mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors

  # Compute (log) likelihood ratios.
  logLRmix <- logl1mix - logl0mix

  z <- drop(z)
  z_std_diff <- drop(z_std_diff)
  res <- data.frame(
    z = z,
    condmean = condmean,
    condvar = condvar,
    z_std_diff = z_std_diff,
    logLR = logLRmix
  )
  p <- ggplot(res, aes_string(y = "z", x = "condmean")) +
    geom_point() +
    labs(y = "Observed z scores", x = "Expected value") +
    geom_abline(intercept = 0, slope = 1) +
    theme_bw()
  idx <- which(logLRmix > 2 & abs(z) > 2)
  if (length(idx) > 0) {
    p <- p + geom_point(
      data = res[idx, ],
      aes_string(y = "z", x = "condmean"), col = "red"
    )
  }
  return(list(plot = p, conditional_dist = res))
}
