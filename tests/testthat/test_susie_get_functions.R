devtools::load_all(".")


context("susie_get_* functions")

# =============================================================================
# Get Model Information
# =============================================================================

test_that("susie_get_objective returns last ELBO when last_only=TRUE", {
  set.seed(1)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  obj <- susie_get_objective(fit, last_only = TRUE)

  expect_type(obj, "double")
  expect_length(obj, 1)
  expect_equal(obj, fit$elbo[length(fit$elbo)])
})

test_that("susie_get_objective returns full ELBO vector when last_only=FALSE", {
  set.seed(2)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  obj <- susie_get_objective(fit, last_only = FALSE)

  expect_type(obj, "double")
  expect_equal(length(obj), fit$niter)
  expect_equal(obj, fit$elbo)
})

test_that("susie_get_objective detects ELBO decrease", {
  set.seed(3)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$elbo <- c(-100, -90, -95, -85)

  expect_message(
    susie_get_objective(fit, warning_tol = 1e-6),
    "Objective is decreasing"
  )
})

# =============================================================================
# Get Posterior Quantities
# =============================================================================

test_that("susie_get_posterior_mean computes correctly", {
  set.seed(4)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pm <- susie_get_posterior_mean(fit)

  # Should return p-length vector
  expect_length(pm, dat$p)
  expect_type(pm, "double")

  # Manual calculation
  expected <- colSums(fit$alpha * fit$mu) / fit$X_column_scale_factors
  expect_equal(pm, expected)
})

test_that("susie_get_posterior_mean filters effects with V < prior_tol", {
  set.seed(5)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Set some V values to zero
  fit$V[c(1, 3)] <- 0

  pm <- susie_get_posterior_mean(fit, prior_tol = 1e-9)

  # Only effects 2, 4, 5 should contribute
  expected <- colSums((fit$alpha * fit$mu)[c(2, 4, 5), , drop = FALSE]) /
              fit$X_column_scale_factors
  expect_equal(pm, expected)
})

test_that("susie_get_posterior_mean returns zeros when all V=0", {
  set.seed(6)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Set all V to zero
  fit$V <- rep(0, 5)

  pm <- susie_get_posterior_mean(fit)

  expect_length(pm, dat$p)
  expect_equal(pm, rep(0, dat$p))
})

test_that("susie_get_posterior_sd computes correctly", {
  set.seed(7)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  psd <- susie_get_posterior_sd(fit)

  # Should return p-length vector
  expect_length(psd, dat$p)
  expect_type(psd, "double")
  expect_true(all(psd >= 0))  # SD must be non-negative

  # Manual calculation
  expected <- sqrt(colSums(fit$alpha * fit$mu2 - (fit$alpha * fit$mu)^2)) /
              fit$X_column_scale_factors
  expect_equal(psd, expected)
})

test_that("susie_get_posterior_sd filters effects with V < prior_tol", {
  set.seed(8)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Set some V values to zero
  fit$V[c(2, 4)] <- 0

  psd <- susie_get_posterior_sd(fit, prior_tol = 1e-9)

  # Only effects 1, 3, 5 should contribute
  expected <- sqrt(colSums((fit$alpha * fit$mu2 -
                           (fit$alpha * fit$mu)^2)[c(1, 3, 5), , drop = FALSE])) /
              fit$X_column_scale_factors
  expect_equal(psd, expected)
})

test_that("susie_get_niter returns correct iteration count", {
  set.seed(9)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, max_iter = 50, verbose = FALSE)

  niter <- susie_get_niter(fit)

  expect_type(niter, "integer")
  expect_equal(niter, fit$niter)
  expect_equal(niter, length(fit$elbo))
})

test_that("susie_get_prior_variance returns V", {
  set.seed(10)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  V <- susie_get_prior_variance(fit)

  expect_equal(V, fit$V)
  expect_length(V, 5)
  expect_true(all(V >= 0))
})

test_that("susie_get_residual_variance returns sigma2", {
  set.seed(11)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  sigma2 <- susie_get_residual_variance(fit)

  expect_type(sigma2, "double")
  expect_length(sigma2, 1)
  expect_equal(sigma2, fit$sigma2)
  expect_true(sigma2 > 0)
})

test_that("susie_get_lfsr computes local false sign rate", {
  set.seed(17)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  lfsr <- susie_get_lfsr(fit)

  # Should return L-length vector (one per effect)
  expect_length(lfsr, 5)
  expect_type(lfsr, "double")

  # LFSR should be in [0, 1]
  expect_true(all(lfsr >= 0 & lfsr <= 1))
})

test_that("susie_get_posterior_samples generates samples", {
  set.seed(18)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  num_samples <- 100
  samples <- susie_get_posterior_samples(fit, num_samples = num_samples)

  # Should return list with b and gamma
  expect_type(samples, "list")
  expect_named(samples, c("b", "gamma"))

  # Check dimensions
  expect_equal(dim(samples$b), c(dat$p, num_samples))
  expect_equal(dim(samples$gamma), c(dat$p, num_samples))

  # Gamma should be binary
  expect_true(all(samples$gamma %in% c(0, 1)))

  # b should be non-zero only where gamma is 1
  for (i in 1:num_samples) {
    expect_true(all((samples$b[, i] != 0) == (samples$gamma[, i] == 1)))
  }
})

test_that("susie_get_posterior_samples filters effects with V < 1e-9", {
  set.seed(19)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Set all V to zero (no effects)
  fit$V <- rep(0, 5)

  samples <- susie_get_posterior_samples(fit, num_samples = 50)

  # With all V=0, all samples should be zero
  expect_true(all(samples$b == 0))
  expect_true(all(samples$gamma == 0))
})

# =============================================================================
# Get Credible Sets and Correlations
# =============================================================================

test_that("susie_get_cs identifies credible sets", {
  set.seed(20)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  cs <- susie_get_cs(fit, coverage = 0.95)

  # Should return list with cs, coverage, requested_coverage
  expect_type(cs, "list")
  expect_true("cs" %in% names(cs))
  expect_true("coverage" %in% names(cs))
  expect_true("requested_coverage" %in% names(cs))

  expect_equal(cs$requested_coverage, 0.95)

  # If CS found, check structure
  if (!is.null(cs$cs)) {
    expect_type(cs$cs, "list")
    expect_true(all(sapply(cs$cs, is.numeric)))
    expect_equal(length(cs$cs), length(cs$coverage))
  }
})

test_that("susie_get_cs filters by purity when X provided", {
  set.seed(21)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  cs_with_purity <- susie_get_cs(fit, X = dat$X, min_abs_corr = 0.5, coverage = 0.95)

  # Should have purity and cs_index fields when X provided
  if (!is.null(cs_with_purity$cs)) {
    expect_true("purity" %in% names(cs_with_purity))
    expect_true("cs_index" %in% names(cs_with_purity))

    # Purity should be data frame with min, mean, median
    expect_s3_class(cs_with_purity$purity, "data.frame")
    expect_true(all(c("min.abs.corr", "mean.abs.corr", "median.abs.corr") %in%
                    colnames(cs_with_purity$purity)))

    # All purity values should be >= min_abs_corr
    expect_true(all(cs_with_purity$purity$min.abs.corr >= 0.5))
  }
})

test_that("susie_get_cs handles dedup parameter", {
  set.seed(22)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  cs_dedup <- susie_get_cs(fit, coverage = 0.95, dedup = TRUE)
  cs_no_dedup <- susie_get_cs(fit, coverage = 0.95, dedup = FALSE)

  # With dedup=TRUE, should have <= CS than without
  n_cs_dedup <- if (is.null(cs_dedup$cs)) 0 else length(cs_dedup$cs)
  n_cs_no_dedup <- if (is.null(cs_no_dedup$cs)) 0 else length(cs_no_dedup$cs)
  expect_true(n_cs_dedup <= n_cs_no_dedup)
})

test_that("get_cs_correlation computes correlations between CS", {
  set.seed(23)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 1) {
    cs_corr <- get_cs_correlation(fit, X = dat$X)

    expect_true(is.matrix(cs_corr))
    expect_equal(nrow(cs_corr), length(fit$sets$cs))
    expect_equal(ncol(cs_corr), length(fit$sets$cs))

    expect_equal(as.numeric(diag(cs_corr)), rep(1, nrow(cs_corr)))
  } else {
    skip("No multiple CS found for correlation test")
  }
})

test_that("get_cs_correlation with Xcorr instead of X", {
  set.seed(24)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 1) {
    Xcorr <- cor(dat$X)
    cs_corr <- get_cs_correlation(fit, Xcorr = Xcorr)

    expect_true(is.matrix(cs_corr))
    expect_equal(nrow(cs_corr), length(fit$sets$cs))
  } else {
    skip("No multiple CS found for Xcorr test")
  }
})

# =============================================================================
# Get PIPs and Related Functions
# =============================================================================

test_that("susie_get_pip computes PIPs correctly", {
  set.seed(12)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pip <- susie_get_pip(fit)

  # Should return p-length vector
  expect_length(pip, dat$p)
  expect_type(pip, "double")

  # All PIPs should be in [0, 1]
  expect_true(all(pip >= 0 & pip <= 1))

  # Manual calculation: 1 - prod(1 - alpha)
  expected <- 1 - apply(1 - fit$alpha, 2, prod)
  expect_equal(pip, expected)
})

test_that("susie_get_pip handles null_index correctly", {
  set.seed(13)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, null_weight = 0.1, verbose = FALSE)

  pip <- susie_get_pip(fit)

  expect_length(pip, dat$p)
  expect_true(all(pip >= 0 & pip <= 1))

  if (!is.null(fit$null_index) && fit$null_index > 0) {
    expect_true(ncol(fit$alpha) == dat$p + 1)
  }
})

test_that("susie_get_pip filters by prior_tol", {
  set.seed(14)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Set some V to zero
  fit$V[c(1, 5)] <- 0

  pip <- susie_get_pip(fit, prior_tol = 1e-9)

  # Only effects 2, 3, 4 should contribute
  expected <- 1 - apply(1 - fit$alpha[c(2, 3, 4), , drop = FALSE], 2, prod)
  expect_equal(pip, expected)
})

test_that("susie_get_pip with prune_by_cs filters to CS effects only", {
  set.seed(15)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  # Get CS
  fit$sets <- susie_get_cs(fit, coverage = 0.95)

  pip_pruned <- susie_get_pip(fit, prune_by_cs = TRUE)

  # Should still return p-length vector
  expect_length(pip_pruned, dat$p)
  expect_true(all(pip_pruned >= 0 & pip_pruned <= 1))

  # If there are CS, pruned PIPs should be different from unpruned
  if (!is.null(fit$sets$cs_index)) {
    pip_full <- susie_get_pip(fit, prune_by_cs = FALSE)
    # At least some should differ (unless all effects are in CS)
    expect_true(any(pip_pruned != pip_full) || length(fit$sets$cs_index) == nrow(fit$alpha))
  }
})

test_that("susie_get_pip returns zeros when no CS and prune_by_cs=TRUE", {
  set.seed(16)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Force no credible sets
  fit$sets <- list(cs = NULL, cs_index = NULL)

  pip <- susie_get_pip(fit, prune_by_cs = TRUE)

  expect_length(pip, dat$p)
  # When no CS, all PIPs should be zero
  expect_true(all(pip == 0))
})

# =============================================================================
# Initialization Functions
# =============================================================================

test_that("susie_init_coef creates valid initialization object", {
  p <- 100
  coef_index <- c(5, 20, 45, 80)
  coef_value <- c(1.5, -2.0, 0.8, -1.2)

  init <- susie_init_coef(coef_index, coef_value, p)

  # Should return susie object
  expect_s3_class(init, "susie")
  expect_type(init, "list")

  # Should have required fields
  expect_true(all(c("alpha", "mu", "mu2", "V") %in% names(init)))

  # Check dimensions
  L <- length(coef_index)
  expect_equal(dim(init$alpha), c(L, p))
  expect_equal(dim(init$mu), c(L, p))
  expect_equal(dim(init$mu2), c(L, p))
  expect_length(init$V, L)
})

test_that("susie_init_coef sets alpha correctly", {
  p <- 50
  coef_index <- c(10, 25, 40)
  coef_value <- c(1.0, 2.0, 3.0)

  init <- susie_init_coef(coef_index, coef_value, p)

  # Alpha should be indicator matrix
  for (i in seq_along(coef_index)) {
    expect_equal(init$alpha[i, coef_index[i]], 1)
    expect_equal(sum(init$alpha[i, ]), 1)  # Each row sums to 1
    expect_equal(sum(init$alpha[i, -coef_index[i]]), 0)  # All others are 0
  }
})

test_that("susie_init_coef sets mu and mu2 correctly", {
  p <- 50
  coef_index <- c(10, 25, 40)
  coef_value <- c(1.5, -2.0, 0.8)

  init <- susie_init_coef(coef_index, coef_value, p)

  # Mu should have coef_value at coef_index
  for (i in seq_along(coef_index)) {
    expect_equal(init$mu[i, coef_index[i]], coef_value[i])
    expect_equal(sum(init$mu[i, -coef_index[i]]), 0)  # All others are 0
  }

  # mu2 should equal mu^2
  expect_equal(init$mu2, init$mu * init$mu)
})

test_that("susie_init_coef errors on invalid inputs", {
  # No effects
  expect_error(
    susie_init_coef(integer(0), numeric(0), 100),
    "Need at least one non-zero effect"
  )

  # Zero coefficient value
  expect_error(
    susie_init_coef(c(1, 5), c(1.0, 0.0), 100),
    "Input coef_value must be non-zero for all its elements"
  )

  # Mismatched lengths
  expect_error(
    susie_init_coef(c(1, 5, 10), c(1.0, 2.0), 100),
    "Inputs coef_index and coef_value must of the same length"
  )

  # Index out of bounds
  expect_error(
    susie_init_coef(c(1, 5, 150), c(1.0, 2.0, 3.0), 100),
    "Input coef_index exceeds the boundary of p"
  )
})

test_that("susie_init_coef works with susie", {
  set.seed(25)
  n <- 100
  p <- 50

  # Create data with known true effects
  dat <- simulate_regression(n = n, p = p, k = 3)
  true_coef_idx <- which(dat$beta != 0)
  true_coef_val <- dat$beta[true_coef_idx]

  # Initialize with true coefficients
  init <- susie_init_coef(true_coef_idx, true_coef_val, p)

  # Fit susie with initialization
  fit <- susie(dat$X, dat$y, L = 10, model_init = init, verbose = FALSE)

  # Should return valid susie fit
  expect_s3_class(fit, "susie")
  expect_true(!is.null(fit$alpha))
  expect_true(!is.null(fit$mu))
  expect_true(!is.null(fit$elbo))
})
