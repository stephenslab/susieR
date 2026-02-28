#' @title Bayesian Multiple Regression with Mixture-of-Normals Prior (RSS)
#'
#' @description This function performs Bayesian multiple regression with a
#'   mixture-of-normals prior using summary statistics (RSS: Regression with
#'   Summary Statistics). It uses a C++ implementation for efficient computation.
#'
#' @param bhat Numeric vector of observed effect sizes (standardized).
#' @param shat Numeric vector of standard errors of effect sizes.
#' @param R Numeric matrix of the correlation matrix.
#' @param var_y Numeric value of the variance of the outcome.
#'   If NULL, it is set to Inf (effects on standardized scale).
#' @param n Integer value of the sample size.
#' @param sigma2_e Numeric value of the error variance.
#' @param s0 Numeric vector of prior variances for the mixture components.
#' @param w0 Numeric vector of prior weights for the mixture components.
#' @param mu1_init Numeric vector of initial values for the posterior mean of
#'   the coefficients. Default is \code{numeric(0)} (initialize to zero).
#' @param tol Numeric value of the convergence tolerance. Default is 1e-8.
#' @param max_iter Integer value of the maximum number of iterations.
#'   Default is 1e5.
#' @param z Numeric vector of Z-scores. If not provided, computed as
#'   \code{bhat / shat}.
#' @param update_w0 Logical value indicating whether to update the mixture
#'   weights. Default is TRUE.
#' @param update_sigma Logical value indicating whether to update the error
#'   variance. Default is TRUE.
#' @param compute_ELBO Logical value indicating whether to compute the
#'   Evidence Lower Bound (ELBO). Default is TRUE.
#' @param standardize Logical value indicating whether to standardize the
#'   input data. Default is FALSE.
#' @param ncpu An integer specifying the number of CPU cores to use for
#'   parallel computation. Default is 1.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{mu1}{Numeric vector of the posterior mean of the coefficients.}
#'   \item{sigma2_1}{Numeric vector of the posterior variance of the coefficients.}
#'   \item{w1}{Numeric matrix of the posterior assignment probabilities.}
#'   \item{sigma2_e}{Numeric value of the error variance.}
#'   \item{w0}{Numeric vector of the mixture weights.}
#'   \item{ELBO}{Numeric value of the Evidence Lower Bound (if \code{compute_ELBO = TRUE}).}
#' }
#'
#' @examples
#' # Generate example data
#' set.seed(985115)
#' n <- 350
#' p <- 16
#' sigmasq_error <- 0.5
#' zeroes <- rbinom(p, 1, 0.6)
#' beta.true <- rnorm(p, 1, sd = 4)
#' beta.true[zeroes] <- 0
#'
#' X <- cbind(matrix(rnorm(n * p), nrow = n))
#' X <- scale(X, center = TRUE, scale = FALSE)
#' y <- X %*% matrix(beta.true, ncol = 1) + rnorm(n, 0, sqrt(sigmasq_error))
#' y <- scale(y, center = TRUE, scale = FALSE)
#'
#' # Set the prior
#' K <- 9
#' sigma0 <- c(0.001, .1, .5, 1, 5, 10, 20, 30, .005)
#' omega0 <- rep(1 / K, K)
#'
#' # Calculate summary statistics
#' b.hat <- sapply(1:p, function(j) {
#'   summary(lm(y ~ X[, j]))$coefficients[-1, 1]
#' })
#' s.hat <- sapply(1:p, function(j) {
#'   summary(lm(y ~ X[, j]))$coefficients[-1, 2]
#' })
#' R.hat <- cor(X)
#' var_y <- var(y)
#' sigmasq_init <- 1.5
#'
#' # Run mr.ash.rss
#' out <- mr.ash.rss(b.hat, s.hat,
#'   R = R.hat, var_y = var_y, n = n,
#'   sigma2_e = sigmasq_init, s0 = sigma0, w0 = omega0,
#'   mu1_init = rep(0, ncol(X)), tol = 1e-8, max_iter = 1e5,
#'   update_w0 = TRUE, update_sigma = TRUE, compute_ELBO = TRUE,
#'   standardize = FALSE
#' )
#'
#' @export
mr.ash.rss <- function(bhat, shat, R, var_y, n,
                       sigma2_e, s0, w0, mu1_init = numeric(0),
                       tol = 1e-8, max_iter = 1e5, z = numeric(0),
                       update_w0 = TRUE, update_sigma = TRUE,
                       compute_ELBO = TRUE, standardize = FALSE, ncpu = 1L) {
  # Check if ncpu is greater than 0 and is an integer
  if (ncpu <= 0 || !is.integer(ncpu)) {
    stop("ncpu must be a positive integer.")
  }

  if (is.null(var_y)) var_y <- Inf
  if (identical(z, numeric(0))) z <- bhat / shat
  result <- rcpp_mr_ash_rss(
    bhat = bhat, shat = shat, z = z, R = R,
    var_y = var_y, n = n, sigma2_e = sigma2_e,
    s0 = s0, w0 = w0, mu1_init = mu1_init,
    tol = tol, max_iter = max_iter,
    update_w0 = update_w0, update_sigma = update_sigma,
    compute_ELBO = compute_ELBO, standardize = standardize,
    ncpus = ncpu
  )

  return(result)
}
