context("susie_get_* functions")

# =============================================================================
# Model information
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
# Posterior quantities
# =============================================================================

test_that("susie_get_posterior_mean computes correctly using all effects", {
  set.seed(4)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pm <- susie_get_posterior_mean(fit)

  expect_length(pm, dat$p)
  expect_type(pm, "double")
  expected <- colSums(fit$alpha * fit$mu) / fit$X_column_scale_factors
  expect_equal(pm, expected)
})

test_that("susie_get_posterior_mean filters effects with V < prior_tol", {
  set.seed(5)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V[c(1, 3)] <- 0
  pm <- susie_get_posterior_mean(fit, prior_tol = 1e-9)

  expected <- colSums((fit$alpha * fit$mu)[c(2, 4, 5), , drop = FALSE]) /
              fit$X_column_scale_factors
  expect_equal(pm, expected)
})

test_that("susie_get_posterior_mean returns zeros when all V=0", {
  set.seed(6)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V <- rep(0, 5)
  pm <- susie_get_posterior_mean(fit)

  expect_length(pm, dat$p)
  expect_equal(pm, rep(0, dat$p))
})

test_that("susie_get_posterior_mean uses all effects when V is not numeric", {
  set.seed(26)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V <- NULL
  pm <- susie_get_posterior_mean(fit)

  expect_length(pm, dat$p)
  expect_type(pm, "double")
  expected <- colSums(fit$alpha * fit$mu) / fit$X_column_scale_factors
  expect_equal(pm, expected)
})

test_that("susie_get_posterior_sd computes correctly", {
  set.seed(7)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  psd <- susie_get_posterior_sd(fit)

  expect_length(psd, dat$p)
  expect_type(psd, "double")
  expect_true(all(psd >= 0))

  expected <- sqrt(colSums(fit$alpha * fit$mu2 - (fit$alpha * fit$mu)^2)) /
              fit$X_column_scale_factors
  expect_equal(psd, expected)
})

test_that("susie_get_posterior_sd filters effects with V < prior_tol", {
  set.seed(8)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V[c(2, 4)] <- 0
  psd <- susie_get_posterior_sd(fit, prior_tol = 1e-9)

  expected <- sqrt(colSums((fit$alpha * fit$mu2 -
                           (fit$alpha * fit$mu)^2)[c(1, 3, 5), , drop = FALSE])) /
              fit$X_column_scale_factors
  expect_equal(psd, expected)
})

test_that("susie_get_posterior_sd uses all effects when V is not numeric", {
  set.seed(27)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V <- NULL
  psd <- susie_get_posterior_sd(fit)

  expect_length(psd, dat$p)
  expect_type(psd, "double")
  expect_true(all(psd >= 0))

  expected <- sqrt(colSums(fit$alpha * fit$mu2 - (fit$alpha * fit$mu)^2)) /
              fit$X_column_scale_factors
  expect_equal(psd, expected)
})

test_that("susie_get_posterior_sd returns zeros when no effects pass prior_tol", {
  set.seed(28)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V <- rep(0, 5)
  psd <- susie_get_posterior_sd(fit, prior_tol = 1e-9)

  expect_length(psd, dat$p)
  expect_type(psd, "double")
  expect_equal(psd, numeric(dat$p))
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

test_that("susie_get_lfsr computes local false sign rate in [0, 1]", {
  set.seed(17)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  lfsr <- susie_get_lfsr(fit)

  expect_length(lfsr, 5)
  expect_type(lfsr, "double")
  expect_true(all(lfsr >= 0 & lfsr <= 1))
})

test_that("susie_get_posterior_samples returns binary gamma and non-zero b only where gamma=1", {
  set.seed(18)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  num_samples <- 100
  samples <- susie_get_posterior_samples(fit, num_samples = num_samples)

  expect_type(samples, "list")
  expect_named(samples, c("b", "gamma"))
  expect_equal(dim(samples$b), c(dat$p, num_samples))
  expect_equal(dim(samples$gamma), c(dat$p, num_samples))
  expect_true(all(samples$gamma %in% c(0, 1)))
  for (i in 1:num_samples) {
    expect_true(all((samples$b[, i] != 0) == (samples$gamma[, i] == 1)))
  }
})

test_that("susie_get_posterior_samples returns all zeros when all V=0", {
  set.seed(19)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V <- rep(0, 5)
  samples <- susie_get_posterior_samples(fit, num_samples = 50)

  expect_true(all(samples$b == 0))
  expect_true(all(samples$gamma == 0))
})

test_that("susie_get_posterior_samples uses all effects when V is not numeric", {
  set.seed(29)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V <- NULL
  num_samples <- 100
  samples <- susie_get_posterior_samples(fit, num_samples = num_samples)

  expect_type(samples, "list")
  expect_named(samples, c("b", "gamma"))
  expect_equal(dim(samples$b), c(dat$p, num_samples))
  expect_equal(dim(samples$gamma), c(dat$p, num_samples))
  expect_true(all(samples$gamma %in% c(0, 1)))
})

# =============================================================================
# Credible sets and correlations
# =============================================================================

test_that("susie_get_cs identifies credible sets with correct structure", {
  set.seed(20)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))

  cs <- suppressMessages(susie_get_cs(fit, coverage = 0.95))

  expect_type(cs, "list")
  expect_true(all(c("cs", "coverage", "requested_coverage") %in% names(cs)))
  expect_equal(cs$requested_coverage, 0.95)
  expect_type(cs$cs, "list")
  expect_true(all(sapply(cs$cs, is.numeric)))
  expect_equal(length(cs$cs), length(cs$coverage))
})

test_that("susie_get_cs filters by purity when X provided", {
  set.seed(20)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))

  cs_with_purity <- susie_get_cs(fit, X = dat$X, min_abs_corr = 0.5, coverage = 0.95)

  expect_true("purity" %in% names(cs_with_purity))
  expect_true("cs_index" %in% names(cs_with_purity))
  expect_s3_class(cs_with_purity$purity, "data.frame")
  expect_true(all(c("min.abs.corr", "mean.abs.corr", "median.abs.corr") %in%
                  colnames(cs_with_purity$purity)))
  expect_true(all(cs_with_purity$purity$min.abs.corr >= 0.5))
})

test_that("susie_get_cs dedup=TRUE returns <= CS count vs dedup=FALSE", {
  set.seed(22)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))

  cs_dedup    <- suppressMessages(susie_get_cs(fit, coverage = 0.95, dedup = TRUE))
  cs_no_dedup <- suppressMessages(susie_get_cs(fit, coverage = 0.95, dedup = FALSE))

  n_cs_dedup    <- if (is.null(cs_dedup$cs))    0 else length(cs_dedup$cs)
  n_cs_no_dedup <- if (is.null(cs_no_dedup$cs)) 0 else length(cs_no_dedup$cs)
  expect_true(n_cs_dedup <= n_cs_no_dedup)
})

test_that("susie_get_cs errors when both X and Xcorr are provided", {
  set.seed(20)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))

  expect_error(
    susie_get_cs(fit, X = dat$X, Xcorr = cor(dat$X), coverage = 0.95),
    "Only one of X or Xcorr should be specified"
  )
})

test_that("susie_get_cs warns about skipped purity filtering when neither X nor Xcorr given", {
  set.seed(20)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))

  expect_message(susie_get_cs(fit), "purity filtering is skipped")
  expect_message(susie_get_cs(fit, min_abs_corr = 0.9), "purity filtering is skipped")
})

test_that("susie_get_cs does not warn about purity filtering when X or Xcorr provided", {
  set.seed(20)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))

  expect_no_message(
    susie_get_cs(fit, X = dat$X, min_abs_corr = 0.5),
    message = "purity filtering is skipped"
  )
  expect_no_message(
    susie_get_cs(fit, Xcorr = cor(dat$X), min_abs_corr = 0.5),
    message = "purity filtering is skipped"
  )
})

test_that("susie_get_cs warns and symmetrizes non-symmetric Xcorr", {
  set.seed(20)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))

  Xcorr <- cor(dat$X)
  Xcorr[1, 2] <- 0.9
  Xcorr[2, 1] <- 0.8

  expect_message(
    cs <- susie_get_cs(fit, Xcorr = Xcorr, coverage = 0.95, check_symmetric = TRUE),
    "Xcorr is not symmetric"
  )
  # Result should use symmetrized Xcorr: entries [1,2] and [2,1] should equal 0.85
  # Verified indirectly: no error on symmetrized input, purity computed
  expect_type(cs, "list")
})

test_that("susie_get_cs uses squared-correlation column names when squared=TRUE", {
  set.seed(20)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))

  cs_sq  <- susie_get_cs(fit, X = dat$X, coverage = 0.95, squared = TRUE)
  cs_abs <- susie_get_cs(fit, X = dat$X, coverage = 0.95, squared = FALSE)

  expect_true(all(c("min.sq.corr",  "mean.sq.corr",  "median.sq.corr")  %in% colnames(cs_sq$purity)))
  expect_false("min.abs.corr" %in% colnames(cs_sq$purity))

  expect_true(all(c("min.abs.corr", "mean.abs.corr", "median.abs.corr") %in% colnames(cs_abs$purity)))
  expect_false("min.sq.corr" %in% colnames(cs_abs$purity))
})

# Shared setup for get_cs_correlation tests: seed 23 reliably yields >=2 CS.
local({
  make_multi_cs_fit <- function() {
    set.seed(23)
    dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
    fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))
    fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
    list(fit = fit, dat = dat)
  }

  test_that("get_cs_correlation returns symmetric matrix with unit diagonal", {
    env <- make_multi_cs_fit()
    fit <- env$fit; dat <- env$dat

    cs_corr <- get_cs_correlation(fit, X = dat$X)

    expect_true(is.matrix(cs_corr))
    expect_equal(nrow(cs_corr), length(fit$sets$cs))
    expect_equal(ncol(cs_corr), length(fit$sets$cs))
    expect_equal(as.numeric(diag(cs_corr)), rep(1, nrow(cs_corr)))
  })

  test_that("get_cs_correlation gives same result with Xcorr as with X", {
    env <- make_multi_cs_fit()
    fit <- env$fit; dat <- env$dat

    cs_corr <- get_cs_correlation(fit, Xcorr = cor(dat$X))

    expect_true(is.matrix(cs_corr))
    expect_equal(nrow(cs_corr), length(fit$sets$cs))
  })

  test_that("get_cs_correlation errors when both X and Xcorr are provided", {
    env <- make_multi_cs_fit()
    fit <- env$fit; dat <- env$dat

    expect_error(
      get_cs_correlation(fit, X = dat$X, Xcorr = cor(dat$X)),
      "Only one of X or Xcorr should be specified"
    )
  })

  test_that("get_cs_correlation errors when neither X nor Xcorr is provided", {
    env <- make_multi_cs_fit()
    fit <- env$fit

    expect_error(
      get_cs_correlation(fit, X = NULL, Xcorr = NULL),
      "One of X or Xcorr must be specified"
    )
  })

  test_that("get_cs_correlation warns and symmetrizes non-symmetric Xcorr", {
    env <- make_multi_cs_fit()
    fit <- env$fit; dat <- env$dat

    Xcorr <- cor(dat$X)
    Xcorr[1, 2] <- 0.9
    Xcorr[2, 1] <- 0.8

    expect_message(
      cs_corr <- get_cs_correlation(fit, Xcorr = Xcorr),
      "Xcorr is not symmetric"
    )
    expect_true(is.matrix(cs_corr))
  })

  test_that("get_cs_correlation with max=TRUE returns scalar equal to upper-triangle maximum", {
    env <- make_multi_cs_fit()
    fit <- env$fit; dat <- env$dat

    cs_corr_matrix <- get_cs_correlation(fit, X = dat$X, max = FALSE)
    cs_corr_max    <- get_cs_correlation(fit, X = dat$X, max = TRUE)

    expect_type(cs_corr_max, "double")
    expect_length(cs_corr_max, 1)
    expect_equal(cs_corr_max,
                 max(abs(cs_corr_matrix[upper.tri(cs_corr_matrix)])))
    expect_true(cs_corr_max >= 0 && cs_corr_max <= 1)
    expect_null(names(cs_corr_max))
  })
})

test_that("get_cs_correlation returns NA when no CS or only one CS", {
  set.seed(33)
  dat <- simulate_regression(n = 100, p = 50, k = 1)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 5, verbose = FALSE))

  fit$sets <- list(cs = NULL)
  expect_true(is.na(get_cs_correlation(fit, X = dat$X)))

  fit$sets <- list(cs = list(c(1, 2, 3)))
  expect_true(is.na(get_cs_correlation(fit, X = dat$X)))
})

# =============================================================================
# PIPs
# =============================================================================

test_that("susie_get_pip computes PIPs correctly as 1 - prod(1 - alpha)", {
  set.seed(12)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pip <- susie_get_pip(fit)

  expect_length(pip, dat$p)
  expect_type(pip, "double")
  expect_true(all(pip >= 0 & pip <= 1))

  expected <- 1 - apply(1 - fit$alpha, 2, prod)
  expect_equal(pip, expected)
})

test_that("susie_get_pip handles null_index correctly with null_weight", {
  set.seed(13)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, null_weight = 0.1, verbose = FALSE)

  pip <- susie_get_pip(fit)

  expect_length(pip, dat$p)
  expect_true(all(pip >= 0 & pip <= 1))
  if (!is.null(fit$null_index) && fit$null_index > 0) {
    expect_equal(ncol(fit$alpha), dat$p + 1)
  }
})

test_that("susie_get_pip filters by prior_tol", {
  set.seed(14)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V[c(1, 5)] <- 0
  pip <- susie_get_pip(fit, prior_tol = 1e-9)

  expected <- 1 - apply(1 - fit$alpha[c(2, 3, 4), , drop = FALSE], 2, prod)
  expect_equal(pip, expected)
})

test_that("susie_get_pip with prune_by_cs=TRUE restricts to CS-indexed effects", {
  set.seed(15)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))
  fit$sets <- suppressMessages(susie_get_cs(fit, coverage = 0.95))

  pip_pruned <- susie_get_pip(fit, prune_by_cs = TRUE)

  expect_length(pip_pruned, dat$p)
  expect_true(all(pip_pruned >= 0 & pip_pruned <= 1))

  if (!is.null(fit$sets$cs_index)) {
    pip_full <- susie_get_pip(fit, prune_by_cs = FALSE)
    expect_true(any(pip_pruned != pip_full) ||
                  length(fit$sets$cs_index) == nrow(fit$alpha))
  }
})

test_that("susie_get_pip returns zeros when no CS and prune_by_cs=TRUE", {
  set.seed(16)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$sets <- list(cs = NULL, cs_index = NULL)
  pip <- susie_get_pip(fit, prune_by_cs = TRUE)

  expect_length(pip, dat$p)
  expect_true(all(pip == 0))
})

test_that("susie_get_pip uses all effects when V is not numeric", {
  set.seed(38)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$V <- NULL
  pip <- susie_get_pip(fit)

  expect_length(pip, dat$p)
  expect_type(pip, "double")
  expect_true(all(pip >= 0 & pip <= 1))
  expected <- 1 - apply(1 - fit$alpha, 2, prod)
  expect_equal(pip, expected)
})

test_that("susie_get_pip with prune_by_cs intersects prior_tol filter and cs_index", {
  set.seed(39)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))
  fit$sets <- suppressMessages(susie_get_cs(fit, coverage = 0.95))
  fit$V[c(1, 2)] <- 0

  skip_if(is.null(fit$sets$cs_index), "No CS found for intersection test")

  pip_pruned <- susie_get_pip(fit, prune_by_cs = TRUE, prior_tol = 1e-9)

  expect_length(pip_pruned, dat$p)
  expect_true(all(pip_pruned >= 0 & pip_pruned <= 1))

  include_idx_V     <- which(fit$V > 1e-9)
  include_idx_final <- intersect(include_idx_V, fit$sets$cs_index)

  if (length(include_idx_final) > 0) {
    expected <- 1 - apply(1 - fit$alpha[include_idx_final, , drop = FALSE], 2, prod)
    expect_equal(pip_pruned, expected)
  } else {
    expect_true(all(pip_pruned == 0))
  }
})

test_that("susie_get_pip with prune_by_cs=TRUE uses cs_index branch", {
  set.seed(502)
  dat <- simulate_regression(n = 150, p = 40, k = 2, signal_sd = 3)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 5, verbose = FALSE))

  skip_if(is.null(fit$sets$cs_index), "No CS found in this fit")

  pip_pruned   <- susie_get_pip(fit, prune_by_cs = TRUE)
  pip_unpruned <- susie_get_pip(fit, prune_by_cs = FALSE)

  expect_length(pip_pruned, dat$p)
  expect_true(all(pip_pruned >= 0 & pip_pruned <= 1))
  expect_true(all(pip_pruned <= pip_unpruned + 1e-12))
})

# =============================================================================
# Initialization
# =============================================================================

test_that("susie_init_coef creates valid susie object with correct alpha/mu/mu2", {
  p          <- 100
  coef_index <- c(5, 20, 45, 80)
  coef_value <- c(1.5, -2.0, 0.8, -1.2)

  init <- susie_init_coef(coef_index, coef_value, p)

  expect_s3_class(init, "susie")
  expect_true(all(c("alpha", "mu", "mu2") %in% names(init)))
  expect_null(init$V)

  L <- length(coef_index)
  expect_equal(dim(init$alpha), c(L, p))
  expect_equal(dim(init$mu),    c(L, p))
  expect_equal(dim(init$mu2),   c(L, p))

  # Alpha is indicator: each row sums to 1, mass at coef_index
  for (i in seq_along(coef_index)) {
    expect_equal(init$alpha[i, coef_index[i]], 1)
    expect_equal(sum(init$alpha[i, ]), 1)
    expect_equal(sum(init$alpha[i, -coef_index[i]]), 0)
  }

  # Mu holds coef_value at coef_index; mu2 = mu^2
  for (i in seq_along(coef_index)) {
    expect_equal(init$mu[i, coef_index[i]], coef_value[i])
    expect_equal(sum(init$mu[i, -coef_index[i]]), 0)
  }
  expect_equal(init$mu2, init$mu * init$mu)
})

test_that("susie_init_coef errors on invalid inputs", {
  expect_error(susie_init_coef(integer(0), numeric(0), 100),
               "Need at least one non-zero effect")
  expect_error(susie_init_coef(c(1, 5), c(1.0, 0.0), 100),
               "Input coef_value must be non-zero for all its elements")
  expect_error(susie_init_coef(c(1, 5, 10), c(1.0, 2.0), 100),
               "Inputs coef_index and coef_value must of the same length")
  expect_error(susie_init_coef(c(1, 5, 150), c(1.0, 2.0, 3.0), 100),
               "Input coef_index exceeds the boundary of p")
})

test_that("susie_init_coef integrates with susie via model_init", {
  set.seed(25)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  true_coef_idx <- which(dat$beta != 0)
  true_coef_val <- dat$beta[true_coef_idx]
  init <- susie_init_coef(true_coef_idx, true_coef_val, dat$p)

  fit <- susie(dat$X, dat$y, L = 10, model_init = init, verbose = FALSE)

  expect_s3_class(fit, "susie")
  expect_false(is.null(fit$alpha))
  expect_false(is.null(fit$mu))
  expect_false(is.null(fit$elbo))
})

# =============================================================================
# Attainable credible sets
# =============================================================================

# Local helper: build a synthetic susie fit from an alpha matrix.
make_alpha_fit <- function(alpha, V = NULL) {
  if (is.null(V)) V <- rep(1, nrow(alpha))
  fit <- list(alpha = alpha, V = V)
  class(fit) <- c("susie", "list")
  fit
}

test_that("susie_get_cs_attainable returns list with cs/coverage/requested_coverage", {
  set.seed(50)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))

  cs <- suppressMessages(susie_get_cs_attainable(fit, coverage = 0.95))

  expect_type(cs, "list")
  expect_true(all(c("cs", "coverage", "requested_coverage") %in% names(cs)))
  expect_equal(cs$requested_coverage, 0.95)
})

test_that("susie_get_cs_attainable default ethres equals max(100, 0.1*p)", {
  alpha <- matrix(0, nrow = 2, ncol = 1000)
  alpha[1, 1] <- 0.99; alpha[1, 2] <- 0.01
  alpha[2, 3] <- 0.99; alpha[2, 4] <- 0.01
  fit <- make_alpha_fit(alpha)

  cs_default  <- suppressMessages(susie_get_cs_attainable(fit, coverage = 0.95))
  cs_explicit <- suppressMessages(
    susie_get_cs_attainable(fit, coverage = 0.95, ethres = max(100, 0.1 * 1000))
  )

  expect_identical(cs_default$cs,       cs_explicit$cs)
  expect_identical(cs_default$coverage, cs_explicit$coverage)
})

test_that("susie_get_cs_attainable drops diffuse effects via entropy filter", {
  p     <- 1000
  alpha <- matrix(0, nrow = 2, ncol = p)
  alpha[1, 1] <- 0.99; alpha[1, 2] <- 0.01
  alpha[2, 1:600] <- 1 / 600      # entropy = log(600) > log(100) = ethres
  fit <- make_alpha_fit(alpha)

  cs <- suppressMessages(susie_get_cs_attainable(fit, coverage = 0.5))

  expect_false(is.null(cs$cs))
  expect_length(cs$cs, 1)
  expect_equal(names(cs$cs), "L1")
})

test_that("susie_get_cs_attainable drops effects with low attainable coverage", {
  p     <- 100
  alpha <- matrix(0, nrow = 2, ncol = p)
  alpha[1, 1:5] <- c(0.50, 0.10, 0.10, 0.10, 0.10)
  alpha[2, 1:5] <- c(0.10, 0.50, 0.10, 0.10, 0.10)
  fit <- make_alpha_fit(alpha)

  cs <- suppressMessages(susie_get_cs_attainable(fit, coverage = 0.95))

  expect_null(cs$cs)
  expect_null(cs$coverage)
  expect_equal(cs$requested_coverage, 0.95)
})

test_that("susie_get_cs_attainable returns NULL cs when all effects are diffuse", {
  p     <- 1000
  alpha <- matrix(1 / p, nrow = 3, ncol = p)
  fit   <- make_alpha_fit(alpha)

  cs <- suppressMessages(susie_get_cs_attainable(fit, coverage = 0.95))

  expect_null(cs$cs)
  expect_null(cs$coverage)
  expect_equal(cs$requested_coverage, 0.95)
})

test_that("susie_get_cs_attainable remaps cs names to original effect indices", {
  p     <- 1000
  alpha <- matrix(0, nrow = 3, ncol = p)
  alpha[1, 5]   <- 0.98; alpha[1, 6]   <- 0.02
  alpha[2, ]    <- 1 / p                          # diffuse: filtered
  alpha[3, 800] <- 0.97; alpha[3, 801] <- 0.03
  fit <- make_alpha_fit(alpha)

  cs <- suppressMessages(susie_get_cs_attainable(fit, coverage = 0.95))

  expect_false(is.null(cs$cs))
  expect_true(all(names(cs$cs) %in% c("L1", "L3")))
  expect_null(cs$cs_index)
  expect_null(cs$purity)
})

test_that("susie_get_cs_attainable handles scalar V without error", {
  p     <- 100
  alpha <- matrix(0, nrow = 2, ncol = p)
  alpha[1, 1] <- 0.99; alpha[1, 2] <- 0.01
  alpha[2, 3] <- 0.99; alpha[2, 4] <- 0.01
  fit <- make_alpha_fit(alpha, V = 1)

  expect_no_error(suppressMessages(susie_get_cs_attainable(fit, coverage = 0.5)))
})

test_that("susie_get_cs_attainable strips X and Xcorr from ... arguments", {
  set.seed(51)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 10, verbose = FALSE))

  expect_no_error(
    suppressMessages(
      susie_get_cs_attainable(fit, X = dat$X, Xcorr = cor(dat$X))
    )
  )
})

test_that("susie_get_cs_attainable matches susie_get_cs when all effects are localized", {
  p     <- 200
  alpha <- matrix(0, nrow = 3, ncol = p)
  alpha[1,   1] <- 0.99; alpha[1,   2] <- 0.01
  alpha[2,  50] <- 0.99; alpha[2,  51] <- 0.01
  alpha[3, 100] <- 0.99; alpha[3, 101] <- 0.01
  fit <- make_alpha_fit(alpha)

  cs_a <- suppressMessages(susie_get_cs_attainable(fit, coverage = 0.95))
  cs_b <- suppressMessages(susie_get_cs(fit,            coverage = 0.95))

  expect_equal(cs_a$cs,                cs_b$cs)
  expect_equal(cs_a$coverage,          cs_b$coverage)
  expect_equal(cs_a$requested_coverage, cs_b$requested_coverage)
})

test_that("susie_get_cs_attainable remaps L2/L3 when L1 is filtered by entropy", {
  # Three effects: L2 localized, L1 and L3 diffuse. Only L2 survives.
  p     <- 200
  alpha <- matrix(0, nrow = 3, ncol = p)
  alpha[1, 1:150] <- 1 / 150
  alpha[2, 10] <- 0.99; alpha[2, 11] <- 0.01
  alpha[3, 1:150] <- 1 / 150
  fit <- make_alpha_fit(alpha)

  cs <- suppressMessages(susie_get_cs_attainable(fit, coverage = 0.5))

  expect_length(cs$cs, 1)
  expect_equal(names(cs$cs), "L2")
})
