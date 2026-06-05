context("Main susie interface functions")

# =============================================================================
# SUSIE() - BASIC FUNCTIONALITY
# =============================================================================

test_that("susie returns valid susie object with expected fields", {
  set.seed(1)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_s3_class(fit, "susie")
  for (field in c("alpha", "mu", "mu2", "V", "sigma2", "pip", "sets", "elbo")) {
    expect_true(field %in% names(fit), info = paste("missing field:", field))
  }
})

test_that("susie output has correct dimensions", {
  set.seed(2)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_equal(dim(fit$alpha), c(5, 50))
  expect_equal(dim(fit$mu),    c(5, 50))
  expect_equal(dim(fit$mu2),   c(5, 50))
  expect_length(fit$V,       5)
  expect_length(fit$pip,    50)
  expect_length(fit$fitted, 100)
})

test_that("susie alpha rows are valid probability distributions", {
  set.seed(3)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(fit$alpha >= 0 & fit$alpha <= 1))
  expect_true(all(fit$pip   >= 0 & fit$pip   <= 1))
})

test_that("susie ELBO is monotonically increasing", {
  set.seed(4)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_true(all(diff(fit$elbo) > -1e-6))
})

test_that("susie ELBO is monotonically increasing for NIG residuals with L = 1 (data_small)", {
  set.seed(1)
  data(data_small)
  fit <- suppressWarnings(
    susie(data_small$X, data_small$y, L = 1,
          estimate_residual_method = "NIG",
          alpha0 = 0.1, beta0 = 0.1, tol = 1e-6,
          verbose = FALSE)
  )
  expect_true(all(diff(fit$elbo) >= 0))
})

test_that("susie niter is at most max_iter and converged field is present", {
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

test_that("susie caps L at p when L > p", {
  set.seed(7)
  dat <- simulate_regression(n = 100, p = 20, k = 3)

  fit <- susie(dat$X, dat$y, L = 50, verbose = FALSE)

  expect_equal(nrow(fit$alpha), 20)
})

test_that("susie standardize=TRUE sets positive X_column_scale_factors; both variants produce valid output", {
  set.seed(8)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_std   <- susie(dat$X, dat$y, L = 5, standardize = TRUE,  verbose = FALSE)
  fit_nostd <- susie(dat$X, dat$y, L = 5, standardize = FALSE, verbose = FALSE)

  expect_true(all(fit_std$X_column_scale_factors > 0))
  expect_equal(rowSums(fit_std$alpha),   rep(1, 5), tolerance = 1e-10)
  expect_equal(rowSums(fit_nostd$alpha), rep(1, 5), tolerance = 1e-10)
})

test_that("susie intercept=TRUE gives finite intercept; intercept=FALSE gives 0", {
  set.seed(9)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_int   <- susie(dat$X, dat$y, L = 5, intercept = TRUE,  verbose = FALSE)
  fit_noint <- susie(dat$X, dat$y, L = 5, intercept = FALSE, verbose = FALSE)

  expect_true(is.finite(fit_int$intercept))
  expect_equal(fit_noint$intercept, 0)
})

test_that("susie prior_weights shifts posterior mass toward favored variables", {
  set.seed(10)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  # Heavy weight on first 10 variables
  custom_weights <- c(rep(10, 10), rep(1, 40))
  fit_uniform <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit_custom  <- susie(dat$X, dat$y, L = 5, prior_weights = custom_weights, verbose = FALSE)

  # Mean PIP in favored region should be higher under custom weights
  expect_true(mean(fit_custom$pip[1:10]) >= mean(fit_uniform$pip[1:10]))
})

test_that("susie null_weight=0 gives p alpha columns; null_weight>0 gives p+1", {
  set.seed(11)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_nonull <- susie(dat$X, dat$y, L = 5, null_weight = 0,   verbose = FALSE)
  fit_null   <- susie(dat$X, dat$y, L = 5, null_weight = 0.1, verbose = FALSE)

  expect_equal(ncol(fit_nonull$alpha), 50)
  expect_equal(ncol(fit_null$alpha),   51)
})

test_that("susie estimate_residual_variance=FALSE fixes sigma2 at specified value", {
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
  expect_true(fit_est$sigma2 > 0)
  expect_false(isTRUE(all.equal(fit_est$sigma2, fixed_sigma2, tolerance = 1e-3)))
})

test_that("susie estimate_prior_variance=FALSE uses fixed scaled_prior_variance", {
  set.seed(13)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_fixed <- susie(dat$X, dat$y, L = 5,
                     estimate_prior_variance = FALSE,
                     scaled_prior_variance = 0.5,
                     verbose = FALSE)
  fit_est <- susie(dat$X, dat$y, L = 5,
                   estimate_prior_variance = TRUE,
                   verbose = FALSE)

  # Fixed: all V components should equal 0.5 * var(y) (scaled back)
  expect_true(all(fit_fixed$V > 0))
  # Estimated prior should differ from fixed
  expect_false(isTRUE(all.equal(fit_fixed$V, fit_est$V, tolerance = 1e-3)))
})

test_that("susie convergence_method pip and elbo both produce valid output", {
  set.seed(14)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_elbo <- susie(dat$X, dat$y, L = 5, convergence_method = "elbo", verbose = FALSE)
  fit_pip  <- susie(dat$X, dat$y, L = 5, convergence_method = "pip",  verbose = FALSE)

  expect_equal(rowSums(fit_elbo$alpha), rep(1, 5), tolerance = 1e-10)
  expect_equal(rowSums(fit_pip$alpha),  rep(1, 5), tolerance = 1e-10)
  expect_true("converged" %in% names(fit_elbo))
  expect_true("converged" %in% names(fit_pip))
})

test_that("susie compute_univariate_zscore=TRUE attaches z of length p; FALSE gives NULL", {
  set.seed(15)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_noz <- susie(dat$X, dat$y, L = 5, compute_univariate_zscore = FALSE, verbose = FALSE)
  fit_z   <- susie(dat$X, dat$y, L = 5, compute_univariate_zscore = TRUE,  verbose = FALSE)

  expect_null(fit_noz$z)
  expect_length(fit_z$z, 50)
})

# =============================================================================
# SUSIE() - VARIANCE ESTIMATION METHODS
# =============================================================================

test_that("susie estimate_residual_method MoM produces positive sigma2", {
  set.seed(16)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5,
               estimate_residual_method = "MoM",
               verbose = FALSE)

  expect_true(fit$sigma2 > 0)
  expect_true(all(diff(fit$elbo) > -1e-6))
})

test_that("susie estimate_residual_method MLE produces positive sigma2", {
  set.seed(17)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5,
               estimate_residual_method = "MLE",
               verbose = FALSE)

  expect_true(fit$sigma2 > 0)
  expect_true(all(diff(fit$elbo) > -1e-6))
})

test_that("susie estimate_residual_method NIG (L=1) produces positive sigma2 and monotone ELBO", {
  set.seed(18)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 1,
               estimate_residual_method = "NIG",
               verbose = FALSE)

  expect_true(fit$sigma2 > 0)
  expect_true(all(diff(fit$elbo) > -1e-6))
})

test_that("susie NIG rejects improper Inverse-Gamma priors (alpha0=0 or beta0=0)", {
  set.seed(18)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  expect_error(
    susie(dat$X, dat$y, L = 1,
          min_abs_corr = 0, check_null_threshold = -1000,
          estimate_residual_method = "NIG",
          alpha0 = 0, beta0 = 0.5,
          verbose = FALSE),
    "alpha0 > 0 and beta0 > 0"
  )

  expect_error(
    susie(dat$X, dat$y, L = 1,
          min_abs_corr = 0, check_null_threshold = -1000,
          estimate_residual_method = "NIG",
          alpha0 = 0, beta0 = 0,
          verbose = FALSE),
    "alpha0 > 0 and beta0 > 0"
  )

  expect_error(
    susie(dat$X, dat$y, L = 2,
          min_abs_corr = 0, check_null_threshold = -1000,
          estimate_residual_method = "NIG",
          alpha0 = 0, beta0 = 0.5,
          verbose = FALSE),
    "alpha0 > 0 and beta0 > 0"
  )

  # Valid alpha0/beta0 should still work
  fit_l1 <- susie(dat$X, dat$y, L = 1,
                  min_abs_corr = 0, check_null_threshold = -1000,
                  estimate_residual_method = "NIG",
                  alpha0 = 0.1, beta0 = 0.1,
                  verbose = FALSE)
  expect_s3_class(fit_l1, "susie")
  expect_true(fit_l1$sigma2 > 0)

  fit_l2 <- susie(dat$X, dat$y, L = 2,
                  min_abs_corr = 0, check_null_threshold = -1000,
                  estimate_residual_method = "NIG",
                  alpha0 = 0.1, beta0 = 0.1,
                  verbose = FALSE)
  expect_s3_class(fit_l2, "susie")
  expect_true(fit_l2$sigma2 > 0)
})

test_that("susie estimate_prior_method options all produce valid ELBO-monotone fits", {
  set.seed(19)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  for (method in c("optim", "EM", "simple")) {
    fit <- susie(dat$X, dat$y, L = 5, estimate_prior_method = method, verbose = FALSE)
    expect_true(all(diff(fit$elbo) > -1e-6), info = paste("method =", method))
    expect_true(all(fit$V >= 0),             info = paste("method =", method))
  }
})

# =============================================================================
# SUSIE() - UNMAPPABLE EFFECTS
# =============================================================================

test_that("susie unmappable_effects=none does not add theta/tau2 fields", {
  set.seed(20)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, unmappable_effects = "none", verbose = FALSE)

  expect_false("theta" %in% names(fit))
  expect_false("tau2"  %in% names(fit))
})

test_that("susie unmappable_effects=inf adds theta/tau2 of correct length", {
  set.seed(21)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, unmappable_effects = "inf", verbose = FALSE)

  expect_true("theta"  %in% names(fit))
  expect_true("tau2"   %in% names(fit))
  expect_false("omega_weights" %in% names(fit))
  expect_length(fit$theta, 50)
  expect_true(is.finite(fit$tau2))
})

test_that("susie drops unmappable fields when SuSiE-inf model_init used with none", {
  set.seed(21)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_inf <- susie(dat$X, dat$y, L = 5, unmappable_effects = "inf",
                   refine = FALSE, verbose = FALSE)
  fit <- susie(dat$X, dat$y, L = 5, unmappable_effects = "none",
               model_init = fit_inf, verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_false("theta"         %in% names(fit))
  expect_false("tau2"          %in% names(fit))
  expect_false("omega_weights" %in% names(fit))
})

test_that("susie unmappable_effects=ash adds theta/tau2", {
  set.seed(22)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, unmappable_effects = "ash", verbose = FALSE)

  expect_true("theta" %in% names(fit))
  expect_true("tau2"  %in% names(fit))
  expect_length(fit$theta, 50)
})

# =============================================================================
# SUSIE() - SIGNAL RECOVERY
# =============================================================================

test_that("susie top-10 PIPs include at least one true causal variable", {
  set.seed(23)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)

  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  top_vars <- order(fit$pip, decreasing = TRUE)[1:10]
  expect_true(length(intersect(top_vars, dat$causal_idx)) >= 1)
})

test_that("susie median null-variable PIP is below 0.3 on strong signal data", {
  set.seed(24)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)

  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  null_vars <- setdiff(1:dat$p, dat$causal_idx)
  expect_true(median(fit$pip[null_vars]) < 0.3)
})

# =============================================================================
# SUSIE() - EDGE CASES
# =============================================================================

test_that("susie L=1 gives single-row alpha summing to 1", {
  set.seed(25)
  dat <- simulate_regression(n = 100, p = 50, k = 1)

  fit <- susie(dat$X, dat$y, L = 1, verbose = FALSE)

  expect_equal(nrow(fit$alpha), 1)
  expect_equal(sum(fit$alpha), 1, tolerance = 1e-10)
})

test_that("susie handles small p (p=5)", {
  set.seed(26)
  dat <- simulate_regression(n = 100, p = 5, k = 2)

  fit <- susie(dat$X, dat$y, L = 3, verbose = FALSE)

  expect_equal(ncol(fit$alpha), 5)
  expect_equal(rowSums(fit$alpha), rep(1, 3), tolerance = 1e-10)
})

test_that("susie near-zero signal y produces valid output with small PIPs", {
  set.seed(27)
  X <- matrix(rnorm(100 * 30), 100, 30)
  y <- rnorm(100, sd = 1e-4)

  fit <- suppressWarnings(susie(X, y, L = 5, verbose = FALSE))

  expect_s3_class(fit, "susie")
  expect_true(all(fit$pip < 0.5))
  expect_true(fit$sigma2 > 0)
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

test_that("susie_ss returns valid susie object with expected fields", {
  set.seed(29)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  X_centered <- scale(dat$X, center = TRUE, scale = FALSE)
  y_centered <- dat$y - mean(dat$y)

  XtX <- crossprod(X_centered)
  Xty <- as.vector(crossprod(X_centered, y_centered))
  yty <- sum(y_centered^2)

  fit <- susie_ss(XtX, Xty, yty, n = 100, L = 5, verbose = FALSE)

  expect_s3_class(fit, "susie")
  for (field in c("alpha", "mu", "V", "sigma2", "pip")) {
    expect_true(field %in% names(fit), info = paste("missing field:", field))
  }
})

test_that("susie_ss output has correct dimensions", {
  set.seed(30)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)

  expect_equal(dim(fit$alpha), c(5, 50))
  expect_equal(dim(fit$mu),    c(5, 50))
  expect_length(fit$V,   5)
  expect_length(fit$pip, 50)
})

test_that("susie_ss alpha rows are valid probability distributions", {
  set.seed(31)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)

  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(fit$alpha >= 0 & fit$alpha <= 1))
  expect_true(all(fit$pip   >= 0 & fit$pip   <= 1))
})

test_that("susie_ss ELBO is monotonically increasing", {
  set.seed(32)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)

  expect_true(all(diff(fit$elbo) > -1e-6))
})

# =============================================================================
# SUSIE_SS() - CONSISTENCY WITH SUSIE()
# =============================================================================

test_that("susie_ss PIPs, V, and sigma2 agree with susie on the same data", {
  set.seed(33)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_ind <- susie(dat$X, dat$y, L = 5, standardize = TRUE, verbose = FALSE)

  X_centered <- scale(dat$X, center = TRUE, scale = FALSE)
  y_centered <- dat$y - mean(dat$y)
  XtX <- crossprod(X_centered)
  Xty <- as.vector(crossprod(X_centered, y_centered))
  yty <- sum(y_centered^2)

  fit_ss <- susie_ss(XtX, Xty, yty, n = 100, L = 5,
                     X_colmeans = colMeans(dat$X),
                     y_mean = mean(dat$y),
                     standardize = TRUE, verbose = FALSE)

  expect_equal(fit_ind$pip,    fit_ss$pip,    tolerance = 1e-3)
  expect_equal(fit_ind$V,      fit_ss$V,      tolerance = 1e-3)
  expect_equal(fit_ind$sigma2, fit_ss$sigma2, tolerance = 1e-3)
})

# =============================================================================
# SUSIE_SS() - PARAMETER HANDLING
# =============================================================================

test_that("susie_ss X_colmeans/y_mean gives finite intercept; without gives NA intercept", {
  set.seed(34)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit_noint <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)
  fit_int   <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                        X_colmeans = colMeans(dat$X),
                        y_mean     = mean(dat$y),
                        verbose = FALSE)

  expect_true(is.na(fit_noint$intercept))
  expect_true(is.finite(fit_int$intercept))
})

test_that("susie_ss maf filtering keeps only variables above maf_thresh", {
  set.seed(35)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  maf <- runif(50, 0, 0.5)
  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                  maf = maf, maf_thresh = 0.1,
                  verbose = FALSE)

  expect_equal(ncol(fit$alpha), sum(maf > 0.1))
})

test_that("susie_ss check_input=TRUE runs without error", {
  set.seed(36)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                  check_input = TRUE,
                  verbose = FALSE)

  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)
})

test_that("susie_ss unmappable_effects=inf adds theta/tau2", {
  set.seed(37)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                  unmappable_effects = "inf",
                  verbose = FALSE)

  expect_true("theta" %in% names(fit))
  expect_true("tau2"  %in% names(fit))
  expect_length(fit$theta, 50)
})

# =============================================================================
# SUSIE_RSS() - BASIC FUNCTIONALITY (lambda = 0)
# =============================================================================

test_that("susie_rss (lambda=0) returns valid susie object with expected fields", {
  set.seed(39)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5, verbose = FALSE)

  expect_s3_class(fit, "susie")
  for (field in c("alpha", "mu", "V", "pip")) {
    expect_true(field %in% names(fit), info = paste("missing field:", field))
  }
})

test_that("susie_rss (lambda=0) output has correct dimensions", {
  set.seed(40)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5, verbose = FALSE)

  expect_equal(dim(fit$alpha), c(5, 50))
  expect_equal(dim(fit$mu),    c(5, 50))
  expect_length(fit$V,   5)
  expect_length(fit$pip, 50)
})

test_that("susie_rss (lambda=0) alpha rows are valid probability distributions", {
  set.seed(41)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5, verbose = FALSE)

  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(fit$alpha >= 0 & fit$alpha <= 1))
  expect_true(all(fit$pip   >= 0 & fit$pip   <= 1))
})

test_that("susie_rss (lambda=0) accepts bhat/shat instead of z", {
  set.seed(42)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  univar <- univariate_regression(dat$X, dat$y)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(bhat = univar$betahat, shat = univar$sebetahat,
                   R = R, n = 100, L = 5, verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_equal(dim(fit$alpha), c(5, 50))
  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)
})

test_that("susie_rss (lambda=0) maf filtering keeps only variables above maf_thresh", {
  set.seed(43)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))
  maf <- runif(50, 0, 0.5)

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                   maf = maf, maf_thresh = 0.1,
                   verbose = FALSE)

  expect_equal(ncol(fit$alpha), sum(maf > 0.1))
})

# =============================================================================
# SUSIE_RSS_LAMBDA() - BASIC FUNCTIONALITY
# =============================================================================

test_that("susie_rss_lambda (lambda>0) returns valid susie object with expected fields", {
  set.seed(44)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  fit <- susie_rss_lambda(z = setup$z, R = setup$R, L = 5,
                          lambda = 1e-5, verbose = FALSE)

  expect_s3_class(fit, "susie")
  for (field in c("alpha", "mu", "V", "pip")) {
    expect_true(field %in% names(fit), info = paste("missing field:", field))
  }
})

test_that("susie_rss_lambda (lambda>0) output has correct dimensions", {
  set.seed(45)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  fit <- susie_rss_lambda(z = setup$z, R = setup$R, L = 5,
                          lambda = 1e-5, verbose = FALSE)

  expect_equal(dim(fit$alpha), c(5, 50))
  expect_equal(dim(fit$mu),    c(5, 50))
  expect_length(fit$V,   5)
  expect_length(fit$pip, 50)
})

test_that("susie_rss_lambda (lambda>0) alpha rows are valid probability distributions", {
  set.seed(46)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  fit <- susie_rss_lambda(z = setup$z, R = setup$R, L = 5,
                          lambda = 1e-5, verbose = FALSE)

  expect_equal(rowSums(fit$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(fit$alpha >= 0 & fit$alpha <= 1))
  expect_true(all(fit$pip   >= 0 & fit$pip   <= 1))
})

test_that("susie_rss_lambda (lambda>0) maf filtering keeps only variables above maf_thresh", {
  set.seed(47)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)
  maf <- runif(50, 0, 0.5)

  fit <- susie_rss_lambda(z = setup$z, R = setup$R, L = 5,
                          lambda = 1e-5, maf = maf, maf_thresh = 0.1,
                          verbose = FALSE)

  expect_equal(ncol(fit$alpha), sum(maf > 0.1))
})

test_that("susie_rss_lambda (lambda>0) ELBO is monotonically increasing", {
  set.seed(48)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  fit <- susie_rss_lambda(z = setup$z, R = setup$R, L = 5,
                          lambda = 1e-5, verbose = FALSE)

  expect_true(all(diff(fit$elbo) > -1e-6))
})

test_that("susie_rss_lambda (lambda>0) top PIPs include at least one true causal variable", {
  set.seed(49)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5,
                                 signal_sd = 1, seed = NULL)

  fit <- susie_rss_lambda(z = setup$z, R = setup$R, L = 10,
                          lambda = 1e-5, verbose = FALSE)

  top_vars <- order(fit$pip, decreasing = TRUE)[1:10]
  expect_true(length(intersect(top_vars, setup$causal_idx)) >= 1)
})

# =============================================================================
# SUSIE_RSS() / SUSIE_RSS_LAMBDA() - SPLIT API
# =============================================================================

test_that("susie_rss and susie_rss_lambda are separate interfaces; susie_rss rejects lambda", {
  set.seed(50)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit_lambda0   <- susie_rss(z = z_scores, R = R, n = 100, L = 5, verbose = FALSE)
  fit_lambda_pos <- susie_rss_lambda(z = z_scores, R = R, L = 5,
                                     lambda = 1e-5, verbose = FALSE)

  expect_s3_class(fit_lambda0,    "susie")
  expect_s3_class(fit_lambda_pos, "susie")
  expect_error(
    susie_rss(z = z_scores, R = R, n = 100, L = 5, lambda = 1e-5, verbose = FALSE),
    "unused argument"
  )
})

test_that("susie_rss_lambda stores n and applies PVE z-score adjustment", {
  set.seed(51)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  init_no_n <- susie_rss_lambda(z = setup$z, R = setup$R, L = 5,
                                lambda = 1e-5, init_only = TRUE)
  init_n    <- susie_rss_lambda(z = setup$z, R = setup$R, n = 100, L = 5,
                                lambda = 1e-5, init_only = TRUE)
  adj <- (100 - 1) / (setup$z^2 + 100 - 2)

  expect_true(is.na(init_no_n$data$n))
  expect_equal(init_n$data$n, 100L)
  expect_equal(init_no_n$data$z, setup$z)
  expect_equal(init_n$data$z, sqrt(adj) * setup$z)
})

test_that("susie_rss_lambda rejects non-MLE residual variance methods", {
  set.seed(510)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  expect_error(
    susie_rss_lambda(z = setup$z, R = setup$R, L = 5, lambda = 1e-5,
                     estimate_residual_method = "MoM", verbose = FALSE),
    "MLE"
  )
  expect_error(
    susie_rss_lambda(z = setup$z, R = setup$R, L = 5, lambda = 1e-5,
                     estimate_residual_method = "NIG", verbose = FALSE),
    "MLE"
  )
})

test_that("R_finite=FALSE silences in-sample R residual variance warning", {
  z <- c(2, -1, 0.5)
  R <- diag(3)

  expect_message(
    susie_rss(z = z, R = R, n = 100, L = 1, max_iter = 50,
              estimate_residual_variance = TRUE, init_only = TRUE),
    "not recommended unless R is the"
  )

  expect_no_message(
    susie_rss(z = z, R = R, n = 100, L = 1, max_iter = 50,
              estimate_residual_variance = TRUE, R_finite = FALSE,
              init_only = TRUE)
  )
})

test_that("susie_rss_lambda does not expose bhat/shat arguments", {
  set.seed(52)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  expect_error(
    susie_rss_lambda(z = setup$z, R = setup$R, L = 5,
        bhat = rnorm(50), shat = runif(50, 0.5, 1),
        lambda = 1e-5, verbose = FALSE),
    "unused argument"
  )
})

test_that("susie_rss_lambda does not expose var_y argument", {
  set.seed(53)
  setup <- setup_rss_lambda_data(n = 500, p = 50, k = 3, lambda = 1e-5, seed = NULL)

  expect_error(
    susie_rss_lambda(z = setup$z, R = setup$R, L = 5,
        var_y = 1.5, lambda = 1e-5, verbose = FALSE),
    "unused argument"
  )
})

test_that("susie_rss does not expose intercept_value argument", {
  set.seed(54)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  expect_error(
    susie_rss(z = z_scores, R = R, n = 100, L = 5, intercept_value = 0.5, verbose = FALSE),
    "unused argument"
  )
})

# =============================================================================
# SUSIE_RSS() - INPUT VALIDATION
# =============================================================================

test_that("susie_rss errors when neither z nor bhat/shat is provided", {
  R <- diag(50)

  expect_error(
    susie_rss(R = R, n = 100, L = 5, verbose = FALSE),
    "Please provide either z or \\(bhat, shat\\)"
  )
})

test_that("susie_rss errors when both z and bhat/shat are provided", {
  z    <- rnorm(50)
  bhat <- rnorm(50)
  shat <- runif(50, 0.5, 1)
  R    <- diag(50)

  expect_error(
    susie_rss(z = z, bhat = bhat, shat = shat, R = R, n = 100, L = 5, verbose = FALSE),
    "Please provide either z or \\(bhat, shat\\), but not both"
  )
})

test_that("susie_rss does not expose check_R argument", {
  set.seed(55)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  expect_error(
    susie_rss(z = z_scores, R = R, n = 100, L = 5, check_R = TRUE, verbose = FALSE),
    "unused argument"
  )
})

# =============================================================================
# INTEGRATION TESTS - CROSS-METHOD COMPARISONS
# =============================================================================

test_that("susie, susie_ss, and susie_rss give similar PIPs on the same data", {
  set.seed(56)
  dat <- simulate_regression(n = 100, p = 50, k = 3, signal_sd = 2)

  fit_ind <- susie(dat$X, dat$y, L = 5, standardize = TRUE,
                   intercept = TRUE, verbose = FALSE)

  X_centered <- scale(dat$X, center = TRUE, scale = FALSE)
  y_centered <- dat$y - mean(dat$y)
  XtX <- crossprod(X_centered)
  Xty <- as.vector(crossprod(X_centered, y_centered))
  yty <- sum(y_centered^2)

  fit_ss <- susie_ss(XtX, Xty, yty, n = 100, L = 5,
                     X_colmeans = colMeans(dat$X),
                     y_mean = mean(dat$y),
                     standardize = TRUE, verbose = FALSE)

  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))
  fit_rss <- susie_rss(z = z_scores, R = R, n = 100, L = 5,
                       verbose = FALSE, estimate_residual_variance = TRUE)

  expect_equal(fit_ind$pip, fit_ss$pip,  tolerance = 1e-3)
  expect_equal(fit_ind$pip, fit_rss$pip, tolerance = 1e-2)
})

test_that("at least one of the three interfaces finds credible sets on strong signal data", {
  set.seed(57)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)

  fit_ind <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  ss <- compute_summary_stats(dat$X, dat$y)
  fit_ss <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 10, verbose = FALSE)

  ss_full <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss_full, cov2cor(XtX))
  fit_rss <- susie_rss(z = z_scores, R = R, n = 200, L = 10, verbose = FALSE)

  has_cs <- (!is.null(fit_ind$sets$cs) && length(fit_ind$sets$cs) > 0) ||
            (!is.null(fit_ss$sets$cs)  && length(fit_ss$sets$cs)  > 0) ||
            (!is.null(fit_rss$sets$cs) && length(fit_rss$sets$cs) > 0)
  expect_true(has_cs)
})

# =============================================================================
# REFINE PARAMETER
# =============================================================================

test_that("susie refine=TRUE produces ELBO >= refine=FALSE (up to tolerance)", {
  set.seed(58)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)

  fit_norefine <- susie(dat$X, dat$y, L = 10, refine = FALSE, verbose = FALSE)
  fit_refine   <- susie(dat$X, dat$y, L = 10, refine = TRUE,  verbose = FALSE)

  elbo_norefine <- susie_get_objective(fit_norefine, last_only = TRUE)
  elbo_refine   <- susie_get_objective(fit_refine,   last_only = TRUE)

  expect_true(elbo_refine >= elbo_norefine - 1e-6)
})

test_that("susie_ss handles refine=TRUE", {
  set.seed(59)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 10,
                  refine = TRUE, verbose = FALSE)

  expect_equal(rowSums(fit$alpha), rep(1, 10), tolerance = 1e-10)
})

test_that("susie_rss handles refine=TRUE", {
  set.seed(60)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 200, L = 10, refine = TRUE, verbose = FALSE)

  expect_equal(rowSums(fit$alpha), rep(1, 10), tolerance = 1e-10)
})

# =============================================================================
# TRACK_FIT PARAMETER
# =============================================================================

test_that("susie track_fit=TRUE attaches a susie_track trace object", {
  set.seed(61)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, track_fit = TRUE, verbose = FALSE)

  expect_true("trace" %in% names(fit))
  expect_s3_class(fit$trace, "susie_track")
})

test_that("susie_ss track_fit=TRUE attaches a susie_track trace object", {
  set.seed(62)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5,
                  track_fit = TRUE, verbose = FALSE)

  expect_true("trace" %in% names(fit))
  expect_s3_class(fit$trace, "susie_track")
})

test_that("susie_rss track_fit=TRUE attaches a susie_track trace object", {
  set.seed(63)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  fit <- susie_rss(z = z_scores, R = R, n = 100, L = 5, track_fit = TRUE, verbose = FALSE)

  expect_true("trace" %in% names(fit))
  expect_s3_class(fit$trace, "susie_track")
})

# =============================================================================
# VERBOSE OUTPUT
# =============================================================================

test_that("susie verbose=TRUE emits ELBO messages", {
  set.seed(64)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  expect_message(
    susie(dat$X, dat$y, L = 5, verbose = TRUE),
    "ELBO"
  )
})

test_that("susie_ss verbose=TRUE emits ELBO messages", {
  set.seed(65)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  expect_message(
    susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = TRUE),
    "ELBO"
  )
})

test_that("susie_rss verbose=TRUE emits ELBO messages", {
  set.seed(66)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss, cov2cor(XtX))

  expect_message(
    susie_rss(z = z_scores, R = R, n = 100, L = 5, verbose = TRUE),
    "ELBO"
  )
})

# =============================================================================
# MATHEMATICAL PROPERTIES - ALL INTERFACES
# =============================================================================

test_that("all interfaces maintain non-negative prior variances", {
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
  fit_rss <- susie_rss(z = z_scores, R = R, n = 100, L = 5, verbose = FALSE)
  expect_true(all(fit_rss$V >= 0))
})

test_that("all interfaces maintain positive residual variance", {
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
  fit_rss <- susie_rss(z = z_scores, R = R, n = 100, L = 5, verbose = FALSE)
  expect_true(fit_rss$sigma2 > 0)
})

# =============================================================================
# OUTPUT COMPATIBILITY WITH SUSIE_GET FUNCTIONS
# =============================================================================

test_that("all interfaces produce output compatible with susie_get functions", {
  set.seed(69)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit_ind <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  expect_length(susie_get_pip(fit_ind), 50)
  expect_equal(susie_get_objective(fit_ind, last_only = TRUE),
               fit_ind$elbo[length(fit_ind$elbo)], tolerance = 1e-8)

  ss <- compute_summary_stats(dat$X, dat$y)
  fit_ss <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)
  expect_length(susie_get_pip(fit_ss), 50)
  expect_equal(susie_get_objective(fit_ss, last_only = TRUE),
               fit_ss$elbo[length(fit_ss$elbo)], tolerance = 1e-8)

  ss_full <- compute_suff_stat(dat$X, dat$y, standardize = TRUE)
  z_scores <- with(univariate_regression(dat$X, dat$y), betahat / sebetahat)
  R <- with(ss_full, cov2cor(XtX))
  fit_rss <- susie_rss(z = z_scores, R = R, n = 100, L = 5, verbose = FALSE)
  expect_length(susie_get_pip(fit_rss), 50)
  expect_equal(susie_get_objective(fit_rss, last_only = TRUE),
               fit_rss$elbo[length(fit_rss$elbo)], tolerance = 1e-8)
})

# =============================================================================
# MULTI-PANEL SUSIE_RSS
# =============================================================================

test_that("susie_rss multi-panel verbose=TRUE emits Multi-panel message", {
  set.seed(101)
  n <- 200; p <- 20
  X <- scale(matrix(rnorm(n * p), n, p), center = TRUE, scale = FALSE)
  beta <- rep(0, p); beta[c(3, 10)] <- c(1, -1)
  y <- as.vector(X %*% beta + rnorm(n, sd = 0.5))
  ss_full <- compute_suff_stat(X, y, standardize = TRUE)
  z_scores <- with(univariate_regression(X, y), betahat / sebetahat)
  R_single <- with(ss_full, cov2cor(XtX))
  R_list   <- list(R_single, R_single)

  expect_message(
    suppressWarnings(
      susie_rss(z = z_scores, R = R_list, n = n, L = 3,
                verbose = TRUE, max_iter = 5)
    ),
    "Multi-panel"
  )
})

test_that("susie_rss multi-panel attaches single_panel_fits regardless of verbose", {
  set.seed(102)
  n <- 200; p <- 20
  X <- scale(matrix(rnorm(n * p), n, p), center = TRUE, scale = FALSE)
  z <- as.vector(scale(rnorm(p)))
  R_single <- cov2cor(crossprod(X) / (n - 1))
  R_list   <- list(R_single, R_single)

  fit <- suppressWarnings(
    susie_rss(z = z, R = R_list, n = n, L = 3, verbose = FALSE, max_iter = 5)
  )
  expect_true("single_panel_fits" %in% names(fit))
})

test_that("susie_rss_lambda errors when lambda argument is missing", {
  p <- 10
  z <- rnorm(p)
  R <- diag(p)

  expect_error(
    susie_rss_lambda(z = z, R = R),
    "requires lambda"
  )
})
