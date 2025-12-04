context("Main susie interface functions")

# =============================================================================
# SUSIE() - BASIC FUNCTIONALITY
# =============================================================================

test_that("susie returns valid susie object", {
  set.seed(1)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_true("alpha" %in% names(fit))
  expect_true("mu" %in% names(fit))
  expect_true("mu2" %in% names(fit))
  expect_true("V" %in% names(fit))
  expect_true("sigma2" %in% names(fit))
  expect_true("pip" %in% names(fit))
  expect_true("sets" %in% names(fit))
  expect_true("elbo" %in% names(fit))
})

test_that("susie has correct dimensions", {
  set.seed(2)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_equal(dim(fit$alpha), c(5, 50))
  expect_equal(dim(fit$mu), c(5, 50))
  expect_equal(dim(fit$mu2), c(5, 50))
  expect_length(fit$V, 5)
  expect_length(fit$pip, 50)
  expect_length(fit$fitted, 100)
})

test_that("susie maintains valid probability distributions", {
  set.seed(3)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Alpha rows sum to 1
  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)

  # Alpha values are valid probabilities
  expect_true(all(fit$alpha >= 0 & fit$alpha <= 1))

  # PIPs are valid probabilities
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

test_that("susie ELBO is monotonically increasing", {
  set.seed(4)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  elbo_diff <- diff(fit$elbo)
  expect_true(all(elbo_diff > -1e-6))
})

test_that("susie converges within max_iter", {
  set.seed(5)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, max_iter = 100, verbose = FALSE)

  expect_true(fit$niter <= 100)
  expect_true("converged" %in% names(fit))
})

# =============================================================================
# SUSIE() - PARAMETER HANDLING
# =============================================================================

test_that("susie respects L parameter", {
  set.seed(6)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_L3 <- susie(dat$X, dat$y, L = 3, verbose = FALSE)
  fit_L7 <- susie(dat$X, dat$y, L = 7, verbose = FALSE)

  expect_equal(nrow(fit_L3$alpha), 3)
  expect_equal(nrow(fit_L7$alpha), 7)
})

test_that("susie adjusts L when L > p", {
  set.seed(7)
  dat <- simulate_regression(n = 100, p = 20, k = 3)

  fit <- susie(dat$X, dat$y, L = 50, verbose = FALSE)

  expect_equal(nrow(fit$alpha), 20)
})

test_that("susie handles standardize parameter", {
  set.seed(8)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_std <- susie(dat$X, dat$y, L = 5, standardize = TRUE, verbose = FALSE)
  fit_nostd <- susie(dat$X, dat$y, L = 5, standardize = FALSE, verbose = FALSE)

  expect_s3_class(fit_std, "susie")
  expect_s3_class(fit_nostd, "susie")
  expect_true(all(fit_std$X_column_scale_factors > 0))
})

test_that("susie handles intercept parameter", {
  set.seed(9)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_int <- susie(dat$X, dat$y, L = 5, intercept = TRUE, verbose = FALSE)
  fit_noint <- susie(dat$X, dat$y, L = 5, intercept = FALSE, verbose = FALSE)

  expect_true(is.finite(fit_int$intercept))
  expect_equal(fit_noint$intercept, 0)
})

test_that("susie handles prior_weights parameter", {
  set.seed(10)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  # Uniform weights
  fit_uniform <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Custom weights (favor first 10 variables)
  custom_weights <- c(rep(10, 10), rep(1, 40))
  fit_custom <- susie(dat$X, dat$y, L = 5, prior_weights = custom_weights, verbose = FALSE)

  expect_s3_class(fit_uniform, "susie")
  expect_s3_class(fit_custom, "susie")
})

test_that("susie handles null_weight parameter", {
  set.seed(11)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_nonull <- susie(dat$X, dat$y, L = 5, null_weight = 0, verbose = FALSE)
  fit_null <- susie(dat$X, dat$y, L = 5, null_weight = 0.1, verbose = FALSE)

  expect_equal(ncol(fit_nonull$alpha), 50)
  expect_equal(ncol(fit_null$alpha), 51)
})

test_that("susie handles estimate_residual_variance parameter", {
  set.seed(12)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fixed_sigma2 <- 1.5
  fit_fixed <- susie(dat$X, dat$y, L = 5,
                     estimate_residual_variance = FALSE,
                     residual_variance = fixed_sigma2,
                     verbose = FALSE)

  fit_est <- susie(dat$X, dat$y, L = 5,
                   estimate_residual_variance = TRUE,
                   verbose = FALSE)

  expect_equal(fit_fixed$sigma2, fixed_sigma2)
  expect_true(fit_est$sigma2 != var(dat$y))
})

test_that("susie handles estimate_prior_variance parameter", {
  set.seed(13)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_fixed <- susie(dat$X, dat$y, L = 5,
                     estimate_prior_variance = FALSE,
                     scaled_prior_variance = 0.5,
                     verbose = FALSE)

  fit_est <- susie(dat$X, dat$y, L = 5,
                   estimate_prior_variance = TRUE,
                   verbose = FALSE)

  expect_s3_class(fit_fixed, "susie")
  expect_s3_class(fit_est, "susie")
})

test_that("susie handles convergence_method parameter", {
  set.seed(14)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_elbo <- susie(dat$X, dat$y, L = 5, convergence_method = "elbo", verbose = FALSE)
  fit_pip <- susie(dat$X, dat$y, L = 5, convergence_method = "pip", verbose = FALSE)

  expect_s3_class(fit_elbo, "susie")
  expect_s3_class(fit_pip, "susie")
})

test_that("susie handles compute_univariate_zscore parameter", {
  set.seed(15)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_noz <- susie(dat$X, dat$y, L = 5, compute_univariate_zscore = FALSE, verbose = FALSE)
  fit_z <- susie(dat$X, dat$y, L = 5, compute_univariate_zscore = TRUE, verbose = FALSE)

  expect_null(fit_noz$z)
  expect_true(!is.null(fit_z$z))
  expect_length(fit_z$z, 50)
})

# =============================================================================
# SUSIE() - VARIANCE ESTIMATION METHODS
# =============================================================================

test_that("susie handles estimate_residual_method = MoM", {
  set.seed(16)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5,
               estimate_residual_method = "MoM",
               verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_true(fit$sigma2 > 0)
})

test_that("susie handles estimate_residual_method = MLE", {
  set.seed(17)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5,
               estimate_residual_method = "MLE",
               verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_true(fit$sigma2 > 0)
})

test_that("susie handles estimate_residual_method = Servin_Stephens", {
  set.seed(18)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 1,
               estimate_residual_method = "Servin_Stephens",
               verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_true(fit$sigma2 > 0)
})

test_that("susie handles estimate_prior_method options", {
  set.seed(19)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_optim <- susie(dat$X, dat$y, L = 5, estimate_prior_method = "optim", verbose = FALSE)
  fit_em <- susie(dat$X, dat$y, L = 5, estimate_prior_method = "EM", verbose = FALSE)
  fit_simple <- susie(dat$X, dat$y, L = 5, estimate_prior_method = "simple", verbose = FALSE)

  expect_s3_class(fit_optim, "susie")
  expect_s3_class(fit_em, "susie")
  expect_s3_class(fit_simple, "susie")
})

# =============================================================================
# SUSIE() - UNMAPPABLE EFFECTS
# =============================================================================

test_that("susie handles unmappable_effects = none", {
  set.seed(20)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, unmappable_effects = "none", verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_false("theta" %in% names(fit))
  expect_false("tau2" %in% names(fit))
})

test_that("susie handles unmappable_effects = inf", {
  set.seed(21)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, unmappable_effects = "inf", verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_true("theta" %in% names(fit))
  expect_true("tau2" %in% names(fit))
  expect_length(fit$theta, 50)
})

test_that("susie handles unmappable_effects = ash", {
  skip_if_not_installed("mr.ash.alpha")

  set.seed(22)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, unmappable_effects = "ash", verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_true("theta" %in% names(fit))
  expect_true("tau2" %in% names(fit))
})

# =============================================================================
# SUSIE() - SIGNAL RECOVERY
# =============================================================================

test_that("susie identifies true causal variables", {
  set.seed(23)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)

  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  # Top PIPs should include true causal variables
  top_vars <- order(fit$pip, decreasing = TRUE)[1:10]
  overlap <- length(intersect(top_vars, dat$causal_idx))

  expect_true(overlap >= 1)
})

test_that("susie maintains low PIPs for null variables", {
  set.seed(24)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)

  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  null_vars <- setdiff(1:dat$p, dat$causal_idx)
  null_pips <- fit$pip[null_vars]

  expect_true(median(null_pips) < 0.3)
})

# =============================================================================
# SUSIE() - EDGE CASES
# =============================================================================

test_that("susie handles L = 1", {
  set.seed(25)
  dat <- simulate_regression(n = 100, p = 50, k = 1)

  fit <- susie(dat$X, dat$y, L = 1, verbose = FALSE)

  expect_equal(nrow(fit$alpha), 1)
  expect_equal(sum(fit$alpha), 1, tolerance = 1e-10)
})

test_that("susie handles small p", {
  set.seed(26)
  dat <- simulate_regression(n = 100, p = 5, k = 2)

  fit <- susie(dat$X, dat$y, L = 3, verbose = FALSE)

  expect_equal(ncol(fit$alpha), 5)
})

test_that("susie errors on NA values when na.rm = FALSE", {
  set.seed(28)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  dat$y[5] <- NA

  expect_error(
    susie(dat$X, dat$y, L = 5, na.rm = FALSE, verbose = FALSE),
    "must not contain missing values"
  )
})

# =============================================================================
# SUSIE_SS() - BASIC FUNCTIONALITY
# =============================================================================

test_that("susie_ss returns valid susie object", {
  set.seed(29)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  X_centered <- scale(dat$X, center = TRUE, scale = FALSE)
  y_centered <- dat$y - mean(dat$y)

  XtX <- crossprod(X_centered)
  Xty <- as.vector(crossprod(X_centered, y_centered))
  yty <- sum(y_centered^2)

  fit <- susie_ss(XtX, Xty, yty, n = 100, L = 5, verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_true("alpha" %in% names(fit))
  expect_true("mu" %in% names(fit))
  expect_true("V" %in% names(fit))
  expect_true("sigma2" %in% names(fit))
  expect_true("pip" %in% names(fit))
})

test_that("susie_ss has correct dimensions", {
  set.seed(30)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)

  expect_equal(dim(fit$alpha), c(5, 50))
  expect_equal(dim(fit$mu), c(5, 50))
  expect_length(fit$V, 5)
  expect_length(fit$pip, 50)
})

test_that("susie_ss maintains valid probability distributions", {
  set.seed(31)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)

  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(fit$alpha >= 0 & fit$alpha <= 1))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

test_that("susie_ss ELBO is monotonically increasing", {
  set.seed(32)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)

  elbo_diff <- diff(fit$elbo)
  expect_true(all(elbo_diff > -1e-6))
})

# =============================================================================
# SUSIE_SS() - CONSISTENCY WITH SUSIE()
# =============================================================================

test_that("susie_ss agrees with susie on same data", {
  set.seed(33)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  # Fit with individual data
  fit_ind <- susie(dat$X, dat$y, L = 5, standardize = TRUE, verbose = FALSE)

  # Compute sufficient statistics
  X_centered <- scale(dat$X, center = TRUE, scale = FALSE)
  y_centered <- dat$y - mean(dat$y)
  XtX <- crossprod(X_centered)
  Xty <- as.vector(crossprod(X_centered, y_centered))
  yty <- sum(y_centered^2)

  # Fit with sufficient statistics
  fit_ss <- susie_ss(XtX, Xty, yty, n = 100, L = 5,
                     X_colmeans = colMeans(dat$X),
                     y_mean = mean(dat$y),
                     standardize = TRUE, verbose = FALSE)

  # Results should be very similar
  expect_equal(fit_ind$pip, fit_ss$pip, tolerance = 1e-3)
  expect_equal(fit_ind$V, fit_ss$V, tolerance = 1e-3)
  expect_equal(fit_ind$sigma2, fit_ss$sigma2, tolerance = 1e-3)
})

# =============================================================================
# SUSIE_SS() - PARAMETER HANDLING
# =============================================================================

test_that("susie_ss handles X_colmeans and y_mean for intercept", {
  set.seed(34)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit_noint <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)

  fit_int <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                      X_colmeans = colMeans(dat$X),
                      y_mean = mean(dat$y),
                      verbose = FALSE)

  expect_true(is.na(fit_noint$intercept))
  expect_true(is.finite(fit_int$intercept))
})

test_that("susie_ss handles maf filtering", {
  set.seed(35)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  maf <- runif(50, 0, 0.5)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                  maf = maf, maf_thresh = 0.1,
                  verbose = FALSE)

  n_filtered <- sum(maf > 0.1)
  expect_equal(ncol(fit$alpha), n_filtered)
})

test_that("susie_ss handles check_input parameter", {
  set.seed(36)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                  check_input = TRUE,
                  verbose = FALSE)

  expect_s3_class(fit, "susie")
})

test_that("susie_ss handles unmappable_effects = inf", {
  set.seed(37)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                  unmappable_effects = "inf",
                  verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_true("theta" %in% names(fit))
  expect_true("tau2" %in% names(fit))
})

# =============================================================================
# SUSIE_RSS() - BASIC FUNCTIONALITY (lambda = 0)
# =============================================================================

test_that("susie_rss with lambda = 0 returns valid susie object", {
  set.seed(39)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                   lambda = 0, verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_true("alpha" %in% names(fit))
  expect_true("mu" %in% names(fit))
  expect_true("V" %in% names(fit))
  expect_true("pip" %in% names(fit))
})

test_that("susie_rss with lambda = 0 has correct dimensions", {
  set.seed(40)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                   lambda = 0, verbose = FALSE)

  expect_equal(dim(fit$alpha), c(5, 50))
  expect_equal(dim(fit$mu), c(5, 50))
  expect_length(fit$V, 5)
  expect_length(fit$pip, 50)
})

test_that("susie_rss with lambda = 0 maintains valid probability distributions", {
  set.seed(41)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                   lambda = 0, verbose = FALSE)

  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(fit$alpha >= 0 & fit$alpha <= 1))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

test_that("susie_rss with lambda = 0 accepts bhat and shat instead of z", {
  set.seed(42)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  univar <- univariate_regression(dat$X, dat$y)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(bhat = univar$betahat, shat = univar$sebetahat,
                   R = R, n = 100, L = 5, lambda = 0, verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_equal(dim(fit$alpha), c(5, 50))
})

test_that("susie_rss with lambda = 0 handles maf filtering", {
  set.seed(43)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))
  maf <- runif(50, 0, 0.5)

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                   lambda = 0, maf = maf, maf_thresh = 0.1,
                   verbose = FALSE)

  n_filtered <- sum(maf > 0.1)
  expect_equal(ncol(fit$alpha), n_filtered)
})

# =============================================================================
# SUSIE_RSS() - BASIC FUNCTIONALITY (lambda > 0)
# =============================================================================

test_that("susie_rss with lambda > 0 returns valid susie object", {
  set.seed(44)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  fit <- susie_rss(z = setup$z, R = setup$R, L = 5,
                   lambda = 1e-5, verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_true("alpha" %in% names(fit))
  expect_true("mu" %in% names(fit))
  expect_true("V" %in% names(fit))
  expect_true("pip" %in% names(fit))
})

test_that("susie_rss with lambda > 0 has correct dimensions", {
  set.seed(45)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  fit <- susie_rss(z = setup$z, R = setup$R, L = 5,
                   lambda = 1e-5, verbose = FALSE)

  expect_equal(dim(fit$alpha), c(5, 50))
  expect_equal(dim(fit$mu), c(5, 50))
  expect_length(fit$V, 5)
  expect_length(fit$pip, 50)
})

test_that("susie_rss with lambda > 0 maintains valid probability distributions", {
  set.seed(46)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  fit <- susie_rss(z = setup$z, R = setup$R, L = 5,
                   lambda = 1e-5, verbose = FALSE)

  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(fit$alpha >= 0 & fit$alpha <= 1))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

test_that("susie_rss with lambda > 0 handles maf filtering", {
  set.seed(47)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)
  maf <- runif(50, 0, 0.5)

  fit <- susie_rss(z = setup$z, R = setup$R, L = 5,
                   lambda = 1e-5, maf = maf, maf_thresh = 0.1,
                   verbose = FALSE)

  n_filtered <- sum(maf > 0.1)
  expect_equal(ncol(fit$alpha), n_filtered)
})

test_that("susie_rss with lambda > 0 ELBO is monotonically increasing", {
  set.seed(48)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  fit <- susie_rss(z = setup$z, R = setup$R, L = 5,
                   lambda = 1e-5, verbose = FALSE)

  elbo_diff <- diff(fit$elbo)
  expect_true(all(elbo_diff > -1e-6))
})

test_that("susie_rss with lambda > 0 identifies causal variables", {
  set.seed(49)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5,
                                  signal_sd = 1, seed = NULL)

  fit <- susie_rss(z = setup$z, R = setup$R, L = 10,
                   lambda = 1e-5, verbose = FALSE)

  # Top PIPs should include at least one true causal variable
  top_vars <- order(fit$pip, decreasing = TRUE)[1:10]
  overlap <- length(intersect(top_vars, setup$causal_idx))

  expect_true(overlap >= 1)
})

# =============================================================================
# SUSIE_RSS() - LAMBDA PARAMETER HANDLING
# =============================================================================

test_that("susie_rss switches data type based on lambda", {
  set.seed(50)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  # lambda = 0 should use sufficient statistics
  fit_lambda0 <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                           lambda = 0, verbose = FALSE)

  # lambda > 0 should use rss_lambda class
  fit_lambda_pos <- susie_rss(z = z_scores, R = R, L = 5,
                              lambda = 1e-5, verbose = FALSE)

  expect_s3_class(fit_lambda0, "susie")
  expect_s3_class(fit_lambda_pos, "susie")
})

test_that("susie_rss with lambda > 0 ignores n parameter", {
  set.seed(51)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  expect_message(
    susie_rss(z = setup$z, R = setup$R, n = 100, L = 5,
              lambda = 1e-5, verbose = FALSE),
    "Parameter 'n' is ignored when lambda != 0"
  )
})

test_that("susie_rss with lambda > 0 rejects bhat/shat", {
  set.seed(52)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  expect_error(
    susie_rss(z = setup$z, R = setup$R, L = 5,
              bhat = rnorm(50), shat = runif(50, 0.5, 1),
              lambda = 1e-5, verbose = FALSE),
    "Parameters 'bhat' and 'shat' are not supported when lambda != 0"
  )
})

test_that("susie_rss with lambda > 0 rejects var_y", {
  set.seed(53)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  expect_error(
    susie_rss(z = setup$z, R = setup$R, L = 5,
              var_y = 1.5, lambda = 1e-5, verbose = FALSE),
    "Parameter 'var_y' is not supported when lambda != 0"
  )
})

test_that("susie_rss with lambda = 0 rejects intercept_value", {
  set.seed(54)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  expect_error(
    susie_rss(z = z_scores, R = R, n = 100, L = 5,
              lambda = 0, intercept_value = 0.5, verbose = FALSE),
    "Parameter 'intercept_value' is only supported when lambda != 0"
  )
})

# =============================================================================
# SUSIE_RSS() - INPUT VALIDATION
# =============================================================================

test_that("susie_rss requires either z or bhat/shat", {
  R <- diag(50)

  expect_error(
    susie_rss(R = R, n = 100, L = 5, lambda = 0, verbose = FALSE),
    "Please provide either z or \\(bhat, shat\\)"
  )
})

test_that("susie_rss rejects both z and bhat/shat", {
  z <- rnorm(50)
  bhat <- rnorm(50)
  shat <- runif(50, 0.5, 1)
  R <- diag(50)

  expect_error(
    susie_rss(z = z, bhat = bhat, shat = shat, R = R, n = 100,
              L = 5, lambda = 0, verbose = FALSE),
    "Please provide either z or \\(bhat, shat\\), but not both"
  )
})

test_that("susie_rss handles check_R parameter", {
  set.seed(55)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                   lambda = 0, check_R = TRUE, verbose = FALSE)

  expect_s3_class(fit, "susie")
})

# =============================================================================
# INTEGRATION TESTS - CROSS-METHOD COMPARISONS
# =============================================================================

test_that("susie, susie_ss, and susie_rss give similar PIPs", {
  set.seed(56)
  dat <- simulate_regression(n = 100, p = 50, k = 3, signal_sd = 2)

  # Fit with individual data
  fit_ind <- susie(dat$X, dat$y, L = 5, standardize = TRUE,
                   intercept = TRUE, verbose = FALSE)

  # Fit with sufficient statistics
  X_centered <- scale(dat$X, center = TRUE, scale = FALSE)
  y_centered <- dat$y - mean(dat$y)
  XtX <- crossprod(X_centered)
  Xty <- as.vector(crossprod(X_centered, y_centered))
  yty <- sum(y_centered^2)

  fit_ss <- susie_ss(XtX, Xty, yty, n = 100, L = 5,
                     X_colmeans = colMeans(dat$X),
                     y_mean = mean(dat$y),
                     standardize = TRUE, verbose = FALSE)

  # Fit with RSS (lambda = 0)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit_rss <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                       lambda = 0, verbose = FALSE, estimate_residual_variance = TRUE)

  # PIPs should be very similar
  expect_equal(fit_ind$pip, fit_ss$pip, tolerance = 1e-3)
  expect_equal(fit_ind$pip, fit_rss$pip, tolerance = 1e-2)
})

test_that("All three interfaces find credible sets", {
  set.seed(57)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)

  fit_ind <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  ss <- compute_summary_stats(dat$X, dat$y)
  fit_ss <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 10, verbose = FALSE)

  ss_full <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss_full, cov2cor(XtX))
  fit_rss <- susie_rss(z = z_scores, R = R, n = 200, L = 10,
                       lambda = 0, verbose = FALSE)

  # At least one method should find credible sets
  has_cs <- (!is.null(fit_ind$sets$cs) && length(fit_ind$sets$cs) > 0) ||
            (!is.null(fit_ss$sets$cs) && length(fit_ss$sets$cs) > 0) ||
            (!is.null(fit_rss$sets$cs) && length(fit_rss$sets$cs) > 0)

  expect_true(has_cs)
})

# =============================================================================
# REFINE PARAMETER
# =============================================================================

test_that("susie handles refine = TRUE", {
  set.seed(58)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)

  fit_norefine <- susie(dat$X, dat$y, L = 10, refine = FALSE, verbose = FALSE)
  fit_refine <- susie(dat$X, dat$y, L = 10, refine = TRUE, verbose = FALSE)

  expect_s3_class(fit_norefine, "susie")
  expect_s3_class(fit_refine, "susie")

  # Refined model should have equal or better ELBO
  elbo_norefine <- susie_get_objective(fit_norefine, last_only = TRUE)
  elbo_refine <- susie_get_objective(fit_refine, last_only = TRUE)

  expect_true(elbo_refine >= elbo_norefine - 1e-6)
})

test_that("susie_ss handles refine = TRUE", {
  set.seed(59)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 10,
                  refine = TRUE, verbose = FALSE)

  expect_s3_class(fit, "susie")
})

test_that("susie_rss handles refine = TRUE", {
  set.seed(60)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 200, L = 10,
                   lambda = 0, refine = TRUE, verbose = FALSE)

  expect_s3_class(fit, "susie")
})

# =============================================================================
# TRACK_FIT PARAMETER
# =============================================================================

test_that("susie handles track_fit = TRUE", {
  set.seed(61)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, track_fit = TRUE, verbose = FALSE)

  expect_true("trace" %in% names(fit))
  expect_type(fit$trace, "list")
})

test_that("susie_ss handles track_fit = TRUE", {
  set.seed(62)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                  track_fit = TRUE, verbose = FALSE)

  expect_true("trace" %in% names(fit))
  expect_type(fit$trace, "list")
})

test_that("susie_rss handles track_fit = TRUE", {
  set.seed(63)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                   lambda = 0, track_fit = TRUE, verbose = FALSE)

  expect_true("trace" %in% names(fit))
  expect_type(fit$trace, "list")
})

# =============================================================================
# VERBOSE PARAMETER
# =============================================================================

test_that("susie verbose output works", {
  set.seed(64)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  expect_message(
    susie(dat$X, dat$y, L = 5, verbose = TRUE),
    "ELBO:"
  )
})

test_that("susie_ss verbose output works", {
  set.seed(65)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  expect_message(
    susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = TRUE),
    "ELBO:"
  )
})

test_that("susie_rss verbose output works", {
  set.seed(66)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  expect_message(
    susie_rss(z = z_scores, R = R, n = 100, L = 5,
              lambda = 0, verbose = TRUE),
    "ELBO:"
  )
})

# =============================================================================
# MATHEMATICAL PROPERTIES - ALL INTERFACES
# =============================================================================

test_that("All interfaces maintain non-negative prior variances", {
  set.seed(67)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_ind <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  expect_true(all(fit_ind$V >= 0))

  ss <- compute_summary_stats(dat$X, dat$y)
  fit_ss <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)
  expect_true(all(fit_ss$V >= 0))

  ss_full <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss_full, cov2cor(XtX))
  fit_rss <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                       lambda = 0, verbose = FALSE)
  expect_true(all(fit_rss$V >= 0))
})

test_that("All interfaces maintain positive residual variance", {
  set.seed(68)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_ind <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  expect_true(fit_ind$sigma2 > 0)

  ss <- compute_summary_stats(dat$X, dat$y)
  fit_ss <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)
  expect_true(fit_ss$sigma2 > 0)

  ss_full <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss_full, cov2cor(XtX))
  fit_rss <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                       lambda = 0, verbose = FALSE)
  expect_true(fit_rss$sigma2 > 0)
})

# =============================================================================
# OUTPUT COMPATIBILITY
# =============================================================================

test_that("All interfaces produce output compatible with susie_get functions", {
  set.seed(69)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  # Test susie
  fit_ind <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  expect_length(susie_get_pip(fit_ind), 50)
  expect_equal(susie_get_objective(fit_ind, last_only = TRUE),
               fit_ind$elbo[length(fit_ind$elbo)])

  # Test susie_ss
  ss <- compute_summary_stats(dat$X, dat$y)
  fit_ss <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)
  expect_length(susie_get_pip(fit_ss), 50)
  expect_equal(susie_get_objective(fit_ss, last_only = TRUE),
               fit_ss$elbo[length(fit_ss$elbo)])

  # Test susie_rss
  ss_full <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss_full, cov2cor(XtX))
  fit_rss <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                       lambda = 0, verbose = FALSE)
  expect_length(susie_get_pip(fit_rss), 50)
  expect_equal(susie_get_objective(fit_rss, last_only = TRUE),
               fit_rss$elbo[length(fit_rss$elbo)])
})
