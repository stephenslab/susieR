context("RSS R-reference mismatch (R_mismatch correction)")

# ---- API surface guards ----

test_that("susie_rss defaults max_iter to 50 with a hint", {
  set.seed(10)
  p <- 10
  n <- 200
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  expect_message(
    obj <- susie_rss(z = z, R = R, n = n, L = 2,
                     init_only = TRUE, verbose = FALSE),
    "Setting max_iter = 50 for the SuSiE RSS model"
  )
  expect_equal(obj$params$max_iter, 50)

  obj2 <- susie_rss(z = z, R = R, n = n, L = 2, max_iter = 7,
                    init_only = TRUE, verbose = FALSE)
  expect_equal(obj2$params$max_iter, 7)
})

test_that("R_mismatch = 'eb' stores Q_art diagnostics in valid ranges", {
  set.seed(11)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- suppressWarnings(susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
                   R_mismatch = "eb", max_iter = 2, verbose = FALSE))
  d <- fit$R_finite_diagnostics
  expect_false(is.null(d$Q_art))
  expect_gte(d$Q_art, 0)
  expect_lte(d$Q_art, 1)
  expect_true(is.logical(d$artifact_flag))
  expect_true(d$mode_label %in% c("normal", "warning", "conservative"))
})

test_that("R_mismatch = 'eb' with B = Inf stores finite r_over_B = 0", {
  set.seed(12)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3,
                   R_mismatch = "eb", max_iter = 2, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  expect_equal(d$B, Inf)
  expect_true(is.finite(d$r_over_B))
  expect_equal(d$r_over_B, 0)
  expect_length(d$lambda_bias, 1)
  expect_gte(d$lambda_bias, 0)
  expect_false(is.null(d$Q_art))
})

test_that("MLE lambda estimator handles zero, interior, and capped fits", {
  lam_zero <- susieR:::estimate_lambda_bias(
    r = c(0.2, -0.1), s = c(1, 2), sigma2 = 1,
    R_finite_B = Inf, method = "eb",
    R_mismatch_method = "mle")
  expect_equal(lam_zero, 0)

  lam_interior <- susieR:::estimate_lambda_bias(
    r = sqrt(1.25), s = 1, sigma2 = 1,
    R_finite_B = Inf, method = "eb",
    R_mismatch_method = "mle")
  expect_equal(lam_interior, 0.25, tolerance = 5e-5)

  lam_small <- susieR:::estimate_lambda_bias(
    r = sqrt(1.01), s = 1, sigma2 = 1,
    R_finite_B = Inf, method = "eb",
    R_mismatch_method = "mle")
  expect_equal(lam_small, 0)

  lam_large <- susieR:::estimate_lambda_bias(
    r = 4, s = 1, sigma2 = 1,
    R_finite_B = 500, method = "eb",
    R_mismatch_method = "mle")
  expect_equal(lam_large, 16 - 1 - 1 / 500, tolerance = 5e-5)
  expect_true(1 / (1 / 500 + lam_large) < 1)
})

test_that("R_mismatch_method mle vs map both yield positive B_corrected on mismatched data", {
  rho <- 0.99
  R <- matrix(c(1, -rho, -rho, 1), 2, 2)
  z <- c(-8, -8)

  fit_mle <- suppressWarnings(susie_rss(
    z = z, R = R, n = 5000, L = 1, R_mismatch = "eb_no_init",
    R_mismatch_method = "mle", max_iter = 3, verbose = FALSE))
  d_mle <- fit_mle$R_finite_diagnostics
  expect_equal(d_mle$R_mismatch_method, "mle")
  expect_gt(d_mle$lambda_bias, 1)
  expect_lt(d_mle$B_corrected, 1)

  fit_map <- suppressWarnings(susie_rss(
    z = z, R = R, n = 5000, L = 1, R_mismatch = "eb_no_init",
    R_mismatch_method = "map", max_iter = 3, verbose = FALSE))
  d_map <- fit_map$R_finite_diagnostics
  expect_equal(d_map$R_mismatch_method, "map")
  expect_lt(d_map$B_corrected, 1)
})

test_that("R_mismatch = 'eb' uses SER-protected initialization", {
  set.seed(13)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)
  z[4] <- 5

  fit <- susie_rss(z = z, R = R, n = n, L = 3,
                   R_mismatch = "eb", max_iter = 2,
                   track_fit = TRUE, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  expect_false(is.null(d$R_mismatch_init))
  expect_equal(d$R_mismatch_init$method, "ser")
  expect_false(is.null(d$R_mismatch_trace))
  expect_equal(d$R_mismatch_trace[[1]]$phase, "init_ser")
  expect_true(is.finite(d$R_mismatch_init$ld_coherence))
  expect_gte(d$R_mismatch_init$ld_coherence, 0)
  expect_lte(d$R_mismatch_init$ld_coherence, 1)
})

test_that("R_mismatch = 'eb_ser_init' preserves finite-B initialization without SER init", {
  set.seed(14)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- suppressWarnings(susie_rss(z = z, R = R, n = n, L = 3,
                   R_finite = 5000, R_mismatch = "eb_ser_init", max_iter = 2,
                   track_fit = TRUE, verbose = FALSE))
  d <- fit$R_finite_diagnostics
  expect_null(d$R_mismatch_init)
  expect_false(is.null(d$Q_art))
})

test_that("R_mismatch = 'eb_force_init' initializes even with finite R_finite", {
  set.seed(16)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3,
                   R_finite = 5000, R_mismatch = "eb_force_init",
                   max_iter = 2, track_fit = TRUE, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  expect_false(is.null(d$R_mismatch_init))
  expect_equal(d$R_mismatch_trace[[1]]$phase, "init_ser")
})

test_that("R_mismatch = 'eb' uses adaptive SER initialization with finite R_finite", {
  set.seed(19)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3,
                   R_finite = 5000, R_mismatch = "eb",
                   max_iter = 2, track_fit = TRUE, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  init <- d$R_mismatch_init
  expect_false(is.null(init))
  expect_false(is.null(d$Q_art))
  expect_equal(d$R_mismatch_trace[[1]]$phase, "init_ser")
  expect_true(is.finite(init$ld_coherence))
  expect_gte(init$ld_coherence, 0)
  expect_lte(init$ld_coherence, 1)
  expect_true(is.finite(init$lambda_bias))
  expect_false("lambda_bias_raw" %in% names(init))
  expect_false("B_corrected" %in% names(init))
})

test_that("R_mismatch = 'eb_no_init' skips SER-protected initialization", {
  set.seed(18)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3,
                   R_mismatch = "eb_no_init", max_iter = 2,
                   track_fit = TRUE, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  expect_null(d$R_mismatch_init)
  expect_false(is.null(d$Q_art))
})

test_that("Optional artifact args validate ranges", {
  set.seed(17)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
              R_mismatch = "eb_no_init", artifact_threshold = -0.1,
              max_iter = 2, verbose = FALSE),
    "artifact_threshold"
  )
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
              R_mismatch = "eb_no_init", artifact_threshold = 1.1,
              max_iter = 2, verbose = FALSE),
    "artifact_threshold"
  )
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
              R_mismatch = "eb_no_init", eig_delta_rel = -1,
              max_iter = 2, verbose = FALSE),
    "eig_delta_rel"
  )
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
              R_mismatch = "eb_no_init", R_sensitivity_threshold = -1,
              max_iter = 2, verbose = FALSE),
    "R_sensitivity_threshold"
  )
})

test_that("BF attenuation diagnostic stores nonnegative BF loss", {
  model <- list(
    alpha = matrix(c(0.8, 0.2), nrow = 1),
    pi = c(0.5, 0.5),
    shat2_inflation = c(4, 1)
  )
  ser_stats <- list(betahat = c(5, 0), shat2 = c(4, 1))
  lbf_adjusted <- susieR:::gaussian_ser_lbf(ser_stats$betahat,
                                            ser_stats$shat2, V = 1)
  out <- susieR:::record_R_bf_attenuation(model, ser_stats,
                                          lbf_adjusted, V = 1, l = 1)
  expect_equal(dim(out$R_bf_attenuation), c(1, 2))
  expect_gt(out$R_bf_attenuation[1, 1], 0)
  expect_gte(out$R_bf_attenuation[1, 2], 0)
})

test_that("BF attenuation summary flags sensitive credible sets", {
  model <- list(
    alpha = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE),
    R_bf_attenuation = matrix(c(log(25), 0, 0, log(2)),
                              nrow = 2, byrow = TRUE),
    sets = list(cs = list(L1 = 1L, L2 = 2L), cs_index = c(1L, 2L)),
    R_finite_diagnostics = list(artifact_flag = FALSE)
  )
  out <- susieR:::summarize_R_bf_attenuation(model, threshold = log(20))
  d <- out$R_finite_diagnostics
  expect_true(d$R_sensitivity_flag)
  expect_equal(d$bf_attenuation$cs_label[["L1"]], "sensitive")
  expect_equal(d$bf_attenuation$cs_label[["L2"]], "stable")
})

# ---- Region-level scalar lambda_bias on the SS path ----

test_that("SS path stores scalar lambda_bias / B_corrected (not per-slot)", {
  set.seed(101)
  p <- 25
  n <- 1500
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
                   R_mismatch = "eb_no_init", max_iter = 5, verbose = FALSE)
  expect_length(fit$R_finite_diagnostics$lambda_bias, 1)
  expect_length(fit$R_finite_diagnostics$B_corrected, 1)
  expect_gte(fit$R_finite_diagnostics$lambda_bias, 0)
  expect_equal(fit$R_finite_diagnostics$B_corrected,
               1 / (1 / fit$R_finite_diagnostics$B +
                      fit$R_finite_diagnostics$lambda_bias),
               tolerance = 1e-12)
})

# ---- Q_art unit tests ----

test_that("compute_Q_art recovers Q ~ 1 when r_fit lies in low-eigen direction", {
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  r_fit <- c(0, 0, 1)
  out <- susieR:::compute_Q_art(eig, r_fit)
  expect_equal(out$Q_art, 1, tolerance = 1e-12)
  expect_true(out$evaluable)
  expect_equal(out$low_eigen_count, 1L)
})

test_that("compute_Q_art returns Q ~ 0 when r_fit avoids low-eigen directions", {
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  r_fit <- c(1, 0.5, 0)
  out <- susieR:::compute_Q_art(eig, r_fit)
  expect_equal(out$Q_art, 0, tolerance = 1e-12)
})

test_that("compute_Q_art is non-evaluable when r_fit has negligible energy", {
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  out <- susieR:::compute_Q_art(eig, rep(0, 3))
  expect_equal(out$Q_art, 0)
  expect_false(out$evaluable)
})

test_that("compute_Q_art is non-evaluable when no low-eigenvalues exist", {
  V <- diag(3)
  d <- c(2, 1, 0.5)
  eig <- list(values = d, vectors = V)
  out <- susieR:::compute_Q_art(eig, c(1, 0, 0))
  expect_equal(out$low_eigen_count, 0L)
  expect_false(out$evaluable)
})

test_that("compute_Q_art stays in [0, 1] for typical inputs", {
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  for (r_fit in list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1),
                     c(0.5, 0.5, 0.5), c(-1, 1, -1))) {
    out <- susieR:::compute_Q_art(eig, r_fit)
    expect_gte(out$Q_art, 0)
    expect_lte(out$Q_art, 1)
  }
})

# ---- eb end-to-end smoke ----

test_that("eb on well-behaved data yields Q_art near 0 and no flag", {
  set.seed(11)
  p <- 25
  n <- 2000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)
  z[3] <- 4

  fit <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
                   R_mismatch = "eb", max_iter = 5, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  expect_lt(d$Q_art, 0.1)
  expect_false(d$artifact_flag)
  expect_equal(d$mode_label, "normal")
})

test_that("eb emits one final R warning when reliability flag triggers", {
  rho <- 0.9999
  z <- c(-8, -8)
  R <- matrix(c(1, -rho, -rho, 1), 2, 2)
  warnings <- character()
  fit <- withCallingHandlers(
    susie_rss(z = z, R = R, n = 5000, L = 1,
              R_mismatch = "eb", max_iter = 5,
              estimate_prior_variance = FALSE,
              estimate_residual_variance = FALSE, verbose = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  final_warnings <- grep("Possible summary statistics and R reference mismatch detected",
                         warnings, value = TRUE)
  expect_length(final_warnings, 1)
  expect_match(final_warnings, "fit\\$R_finite_diagnostics\\$ser_model")
  expect_match(final_warnings, "one-effect credible-set model as in Maller et al. 2012")
  expect_true(fit$R_finite_diagnostics$artifact_flag)
  expect_true(fit$R_finite_diagnostics$R_reliability_flag)
  expect_equal(fit$R_finite_diagnostics$Q_art, 1, tolerance = 1e-6)
  ser <- fit$R_finite_diagnostics$ser_model
  expect_s3_class(ser, "susie")
  expect_equal(nrow(ser$alpha), 1)
  expect_equal(as.vector(ser$alpha[1, ]), ser$pip)
  expect_named(ser$sets, c("cs", "purity", "cs_index", "coverage",
                           "requested_coverage"))
})

test_that("eb surfaces all expected Q_art and mode_label diagnostic fields", {
  set.seed(12)
  p <- 25; n <- 2000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
                   R_mismatch = "eb", max_iter = 3, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  for (fld in c("Q_art", "artifact_flag", "artifact_evaluable",
                "low_eigen_count", "low_eigen_fraction", "eig_delta",
                "mode_label", "lambda_bias", "B_corrected"))
    expect_false(is.null(d[[fld]]), info = paste("missing diagnostic:", fld))
})

test_that("eb with X-input runs and surfaces Q_art in [0, 1]", {
  set.seed(15)
  p <- 25; n <- 2000
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = TRUE)
  beta <- rep(0, p); beta[5] <- 0.4
  y <- drop(X %*% beta + rnorm(n))
  z <- as.numeric(crossprod(X, y) / sqrt(diag(crossprod(X))))
  fit <- susie_rss(z = z, X = X, n = n, L = 3, R_finite = 5000,
                   R_mismatch = "eb", max_iter = 3, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  expect_false(is.null(d$Q_art))
  expect_gte(d$Q_art, 0)
  expect_lte(d$Q_art, 1)
})

test_that("eb on multi-panel SS path with lambda=0 stores scalar lambda_bias", {
  set.seed(19)
  n <- 80
  p <- 12
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p), n, p)
  z <- rnorm(p)

  fit <- susie_rss(z = z, X = list(X1, X2), n = 1000, L = 3,
                   R_finite = TRUE, R_mismatch = "eb", max_iter = 3,
                   verbose = FALSE)
  d <- fit$R_finite_diagnostics
  expect_false(is.null(d$Q_art))
  expect_gte(d$Q_art, 0)
  expect_lte(d$Q_art, 1)
  expect_length(d$lambda_bias, 1)
})

test_that("eb handles multi-panel X input with adaptive initialization", {
  set.seed(21)
  n <- 80
  p <- 12
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p), n, p)
  z <- rnorm(p)

  fit <- susie_rss(z = z, X = list(X1, X2), n = 1000, L = 3,
                   R_finite = TRUE, R_mismatch = "eb",
                   max_iter = 3, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  expect_false(is.null(d$R_mismatch_init))
  expect_true(is.finite(d$R_mismatch_init$ld_coherence))
  expect_gte(d$R_mismatch_init$ld_coherence, 0)
  expect_lte(d$R_mismatch_init$ld_coherence, 1)
  expect_length(d$lambda_bias, 1)
  expect_length(d$B_corrected, 1)
})

test_that("R_mismatch EB initialization respects non-uniform prior_weights", {
  p <- 6
  n <- 2000
  R <- diag(p)
  z <- rep(0, p)
  z[2:3] <- 4
  prior_weights <- rep(0.01, p)
  prior_weights[3] <- 0.95
  prior_weights <- prior_weights / sum(prior_weights)

  fit <- susie_rss(z = z, R = R, n = n, L = 1,
                   prior_weights = prior_weights,
                   R_mismatch = "eb", max_iter = 2,
                   track_fit = TRUE, verbose = FALSE)
  ser <- fit$R_finite_diagnostics$ser_model

  expect_equal(fit$pi, prior_weights)
  expect_equal(ser$pi, prior_weights)
  expect_gt(ser$pip[3], 0.9)
  expect_lt(ser$pip[2], 0.1)
})

test_that("R_mismatch with maf filtering subsets prior_weights", {
  p <- 6
  n <- 2000
  R <- diag(p)
  z <- rnorm(p)
  prior_weights <- c(0.10, 0.20, 0.30, 0.25, 0.10, 0.05)
  maf <- c(0.20, 0.02, 0.30, 0.25, 0.01, 0.40)
  keep <- maf > 0.05

  fit <- susie_rss(z = z, R = R, n = n, L = 1,
                   prior_weights = prior_weights,
                   maf = maf, maf_thresh = 0.05,
                   R_mismatch = "eb", max_iter = 2,
                   verbose = FALSE)

  expect_equal(length(fit$pi), sum(keep))
  expect_equal(fit$pi, prior_weights[keep] / sum(prior_weights[keep]))
})

# ---- resolve_R_finite (direct, all branches) ----

test_that("resolve_R_finite handles NULL / FALSE / TRUE single-panel inputs", {
  expect_null(resolve_R_finite(NULL))
  expect_null(resolve_R_finite(FALSE))

  expect_error(resolve_R_finite(TRUE), "requires X input")

  X <- matrix(0, 120, 5)
  expect_equal(resolve_R_finite(TRUE, X), 120)

  expect_equal(resolve_R_finite(5000), 5000)
})

test_that("resolve_R_finite validates and expands multi-panel inputs", {
  Xl <- list(matrix(0, 50, 5), matrix(0, 70, 5))

  expect_equal(resolve_R_finite(TRUE, Xl, is_multi_panel = TRUE), c(50L, 70L))
  expect_equal(resolve_R_finite(5000, Xl, is_multi_panel = TRUE), c(5000, 5000))
  expect_equal(resolve_R_finite(c(100, 200), Xl, is_multi_panel = TRUE),
               c(100, 200))
  expect_error(
    resolve_R_finite(c(1, 2, 3), Xl, is_multi_panel = TRUE),
    "multi-panel"
  )
})

test_that("resolve_R_finite rejects invalid numeric values", {
  expect_error(resolve_R_finite(-1), "positive numeric")
  expect_error(resolve_R_finite(c(0, 1)), "positive numeric")
  expect_error(resolve_R_finite(Inf), "positive numeric")
  expect_error(resolve_R_finite(c(100, 200)), "single positive number")
})

# ---- compute_R_finite_diagnostics (direct, X / R / identity branches) ----

test_that("compute_R_finite_diagnostics handles the X (factor matrix) branch", {
  set.seed(301)
  X <- matrix(rnorm(100 * 8), 100, 8)

  d <- compute_R_finite_diagnostics(X = X, B = 100, p = 8,
                                    x_is_standardized = FALSE)
  expect_equal(d$B, 100)
  expect_equal(d$p, 8)
  expect_true(is.finite(d$effective_rank) && d$effective_rank > 0)
  expect_equal(d$r_over_B, d$effective_rank / 100, tolerance = 1e-12)
  expect_length(d$Rhat_diag_deviation, 8)

  Xs <- scale(X) / sqrt(99)
  ds <- compute_R_finite_diagnostics(X = Xs, B = 100, p = 8,
                                     x_is_standardized = TRUE)
  expect_true(is.finite(ds$R_frob_sq_debiased))
})

test_that("compute_R_finite_diagnostics handles the R branch and identity fallback", {
  set.seed(302)
  X <- matrix(rnorm(80 * 6), 80, 6)
  R <- cov2cor(crossprod(X))

  d <- compute_R_finite_diagnostics(R = R, B = 5000, p = 6)
  expect_equal(d$p, 6)
  expect_length(d$Rhat_diag_deviation, 6)
  expect_lt(max(d$Rhat_diag_deviation), 1e-8)

  d_id <- compute_R_finite_diagnostics(B = Inf, p = 6)
  expect_equal(d_id$R_frob_sq_debiased, 6)
  expect_equal(d_id$effective_rank, 6)
  expect_equal(d_id$Rhat_diag_deviation, rep(0, 6))
})

# ---- get_mixture_omega (direct, fallback branches) ----

test_that("get_mixture_omega returns NULL for non-mixture or K-less data", {
  expect_null(get_mixture_omega(structure(list(), class = "ss"), list()))
  expect_null(get_mixture_omega(structure(list(), class = "ss_mixture"),
                                list()))
})

test_that("get_mixture_omega normalizes, clamps, and falls back to uniform", {
  data <- structure(list(K = 3, omega_init = c(1, 0, 0)),
                    class = "ss_mixture")

  expect_equal(get_mixture_omega(data, list(omega = c(0.2, 0.3, 0.5))),
               c(0.2, 0.3, 0.5), tolerance = 1e-12)

  expect_equal(get_mixture_omega(data, list()), c(1, 0, 0))

  expect_equal(get_mixture_omega(data, list(omega = c(0.5, 0.5))),
               rep(1 / 3, 3), tolerance = 1e-12)

  expect_equal(get_mixture_omega(data, list(omega = c(-1, 1, 1))),
               c(0, 0.5, 0.5), tolerance = 1e-12)

  expect_equal(get_mixture_omega(data, list(omega = c(0, 0, 0))),
               rep(1 / 3, 3), tolerance = 1e-12)
})

test_that("get_mixture_omega returns uniform when both model omega and data omega_init are NULL", {
  K <- 3
  p <- 10
  data <- structure(list(K = K, p = p, nm1 = 99), class = "ss_mixture")
  model <- list(sigma2 = 1)

  result <- get_mixture_omega(data, model)
  expect_length(result, K)
  expect_equal(sum(result), 1, tolerance = 1e-12)
  expect_true(all(result == 1 / K))
})

# ---- get_R_mismatch_eigen (direct, mixture branches) ----

test_that("get_R_mismatch_eigen prefers model$eigen_R then SS data$eigen_R", {
  eg <- list(values = c(3, 2, 1), vectors = diag(3))
  expect_identical(
    get_R_mismatch_eigen(structure(list(), class = "ss"),
                         list(eigen_R = eg)), eg)
  expect_identical(
    get_R_mismatch_eigen(structure(list(eigen_R = eg), class = "ss"),
                         list()), eg)
})

test_that("get_R_mismatch_eigen returns data$eigen_R as fallback", {
  set.seed(312)
  p <- 6
  R <- diag(p)
  eig <- eigen(R, symmetric = TRUE)

  data <- structure(list(eigen_R = eig, p = p), class = "ss")
  model <- list(sigma2 = 1)

  result <- get_R_mismatch_eigen(data, model)
  expect_equal(result$values, eig$values)
})

test_that("get_R_mismatch_eigen returns NULL when eigen_R is absent everywhere", {
  p <- 4
  data  <- structure(list(eigen_R = NULL), class = "ss")
  model <- list(eigen_R = NULL)

  result <- get_R_mismatch_eigen(data, model)
  expect_null(result)
})

test_that("get_R_mismatch_eigen recovers the mixture spectrum from omega_cache", {
  set.seed(303)
  p <- 5
  r <- 3
  V_s <- qr.Q(qr(matrix(rnorm(p * r), p, r)))
  A1 <- crossprod(matrix(rnorm(r * r), r, r))
  A2 <- crossprod(matrix(rnorm(r * r), r, r))
  cache <- list(V_s = V_s, r = r, A_list = list(A1, A2))
  data <- structure(list(K = 2, p = p, omega_cache = cache,
                         omega_init = c(1, 0)),
                    class = c("ss_mixture", "ss"))

  eig <- get_R_mismatch_eigen(data, list(omega = c(0.5, 0.5)))
  expect_length(eig$values, p)
  expect_true(all(eig$values >= 0))
})

test_that("get_R_mismatch_eigen falls back to panel_R then to X for mixtures", {
  set.seed(304)
  p <- 5
  R1 <- cov2cor(crossprod(matrix(rnorm(30 * p), 30, p)))
  R2 <- cov2cor(crossprod(matrix(rnorm(30 * p), 30, p)))

  data_pr <- structure(list(K = 2, p = p, panel_R = list(R1, R2),
                            omega_init = c(1, 0)),
                       class = c("ss_mixture", "ss"))
  eig_pr <- get_R_mismatch_eigen(data_pr, list(omega = c(0.4, 0.6)))
  expect_length(eig_pr$values, p)
  expect_true(all(eig_pr$values >= 0))

  data_x <- structure(list(K = 2, p = p, X = matrix(rnorm(40 * p), 40, p),
                           nm1 = 39, omega_init = c(1, 0)),
                      class = c("ss_mixture", "ss"))
  eig_x <- get_R_mismatch_eigen(data_x, list(omega = c(0.5, 0.5)))
  expect_length(eig_x$values, p)
  expect_true(all(eig_x$values >= 0))
})

# ---- compute_ser_ld_coherence (direct, all branches) ----

test_that("compute_ser_ld_coherence returns 1 for degenerate inputs", {
  p <- 6
  data <- structure(list(p = p), class = "ss")

  expect_equal(compute_ser_ld_coherence(data, list(pi = rep(1 / p, p)),
                                        rep(0, p)), 1)
  expect_equal(compute_ser_ld_coherence(data, list(pi = rep(1 / p, p)),
                                        rep(1 / p, p)), 1)
  data_nm1 <- structure(list(nm1 = 49, p = p), class = "ss")
  expect_equal(compute_ser_ld_coherence(data_nm1, list(pi = rep(1 / p, p)),
                                        c(0.7, 0.3, 0, 0, 0, 0)), 1)
})

test_that("compute_ser_ld_coherence computes coherence in [0, 1] across R sources", {
  set.seed(305)
  p <- 6
  alpha <- c(0.7, 0.3, 0, 0, 0, 0)
  pi <- rep(1 / p, p)

  XtX <- crossprod(matrix(rnorm(50 * p), 50, p))
  d_xtx <- structure(list(XtX = XtX, nm1 = 49, p = p), class = "ss")
  c_xtx <- compute_ser_ld_coherence(d_xtx, list(pi = pi), alpha)
  expect_gte(c_xtx, 0)
  expect_lte(c_xtx, 1)

  d_x <- structure(list(X = matrix(rnorm(50 * p), 50, p), nm1 = 49, p = p),
                   class = "ss")
  c_x <- compute_ser_ld_coherence(d_x, list(pi = pi), alpha)
  expect_gte(c_x, 0)
  expect_lte(c_x, 1)

  R1 <- cov2cor(crossprod(matrix(rnorm(40 * p), 40, p)))
  R2 <- cov2cor(crossprod(matrix(rnorm(40 * p), 40, p)))
  d_mix <- structure(list(K = 2, p = p, panel_R = list(R1, R2), nm1 = 49,
                          omega_init = c(1, 0)),
                     class = c("ss_mixture", "ss"))
  c_mix <- compute_ser_ld_coherence(d_mix,
                                    list(pi = pi, omega = c(0.5, 0.5)), alpha)
  expect_gte(c_mix, 0)
  expect_lte(c_mix, 1)

  c_badpi <- compute_ser_ld_coherence(d_xtx, list(pi = c(1, 2)), alpha)
  expect_gte(c_badpi, 0)
  expect_lte(c_badpi, 1)
})

# ---- estimate_lambda_bias (MAP method and method='none') ----

test_that("estimate_lambda_bias MAP returns a positive estimate or floors to zero", {
  lam_map <- estimate_lambda_bias(r = 4, s = 1, sigma2 = 1,
                                  R_finite_B = 500, method = "eb",
                                  R_mismatch_method = "map")
  expect_gt(lam_map, 0)

  lam_zero <- estimate_lambda_bias(r = sqrt(1.001), s = 1, sigma2 = 1,
                                   R_finite_B = Inf, method = "eb",
                                   R_mismatch_method = "map")
  expect_equal(lam_zero, 0)
})

test_that("estimate_lambda_bias returns 0 when method='none'", {
  set.seed(308)
  r <- rnorm(20)
  s <- abs(rnorm(20)) + 0.01
  result <- estimate_lambda_bias(r, s, sigma2 = 0.5, R_finite_B = 1000,
                                 method = "none")
  expect_equal(result, 0)
})

# ---- compute_R_mismatch_state on the ss_mixture path (direct) ----

test_that("compute_R_mismatch_state runs on the ss_mixture path and stores Q_art", {
  set.seed(306)
  p <- 8
  n <- 1000
  make_R <- function(s) {
    set.seed(s)
    X <- matrix(rnorm(60 * p), 60, p)
    cov2cor(crossprod(scale(X)))
  }
  R1 <- make_R(1)
  R2 <- make_R(2)
  z <- rnorm(p)
  z[3] <- 4

  ctor <- ss_mixture_constructor(z = z, R = list(R1, R2), n = n, L = 2,
                                 R_mismatch = "eb", R_finite = c(3000, 4000))
  data <- ctor$data
  params <- ctor$params
  expect_false(is.null(data$eigen_R))

  var_y <- get_var_y.ss(data)
  model <- initialize_susie_model.ss(data, params, var_y)
  model$alpha <- matrix(1 / p, params$L, p)
  model$alpha[1, 3] <- 0.8
  model$alpha[1, -3] <- 0.2 / (p - 1)
  model$mu <- matrix(0.02, params$L, p)
  model$mu2 <- model$mu^2 + 0.01
  model$sigma2 <- 1
  model$XtXr <- compute_XtXv_mixture(data, model,
                                     colSums(model$alpha * model$mu))

  res <- compute_R_mismatch_state(data, params, model, phase = "sweep")

  expect_length(res$lambda_bias, 1)
  expect_gte(res$lambda_bias, 0)
  expect_equal(res$B_corrected,
               1 / (1 / data$R_finite_B + res$lambda_bias),
               tolerance = 1e-12)
  expect_false(is.null(res$Q_art))
  expect_gte(res$Q_art, 0)
  expect_lte(res$Q_art, 1)
  expect_true(res$mode_label %in% c("normal", "warning"))
  expect_true(is.logical(res$artifact_flag))
})

test_that("compute_R_mismatch_state returns early when nm1 is invalid", {
  set.seed(307)
  p <- 6
  make_R <- function(s) {
    set.seed(s)
    X <- matrix(rnorm(40 * p), 40, p)
    cov2cor(crossprod(scale(X)))
  }
  ctor <- ss_mixture_constructor(z = rnorm(p), R = list(make_R(1), make_R(2)),
                                 n = 500, L = 2, R_mismatch = "eb",
                                 R_finite = c(2000, 2000))
  data <- ctor$data
  params <- ctor$params
  var_y <- get_var_y.ss(data)
  model <- initialize_susie_model.ss(data, params, var_y)
  model$alpha <- matrix(1 / p, params$L, p)
  model$mu <- matrix(0.01, params$L, p)
  model$mu2 <- model$mu^2 + 0.01
  model$sigma2 <- 1
  model$XtXr <- compute_XtXv_mixture(data, model,
                                     colSums(model$alpha * model$mu))

  data_bad <- data
  data_bad$nm1 <- 0

  res <- compute_R_mismatch_state(data_bad, params, model, phase = "sweep")
  expect_null(res$lambda_bias)
})

test_that("compute_R_mismatch_state returns early when sigma2 is zero", {
  set.seed(310)
  p <- 8
  X <- matrix(rnorm(60 * p), 60, p)
  R <- cov2cor(crossprod(scale(X)))
  z <- rnorm(p)

  ctor <- summary_stats_constructor(
    z = z, R = R, n = 200, L = 2,
    R_mismatch = "eb", R_finite = 5000,
    convergence_method = "pip", coverage = 0.95, min_abs_corr = 0.5,
    n_purity = 100, check_prior = FALSE, track_fit = FALSE
  )
  data <- ctor$data
  params <- ctor$params
  var_y <- get_var_y.ss(data)
  model <- initialize_susie_model.ss(data, params, var_y)
  model$alpha <- matrix(1 / p, 2, p)
  model$mu    <- matrix(0.01, 2, p)
  model$mu2   <- model$mu^2 + 0.01
  model$sigma2 <- 0

  result <- compute_R_mismatch_state(data, params, model)
  expect_null(result$lambda_bias)
})

test_that("compute_R_mismatch_state returns early when data is not ss/ss_mixture", {
  p <- 6
  model <- list(
    sigma2 = 0.5,
    alpha = matrix(1 / p, 2, p),
    mu = matrix(0.01, 2, p),
    R_finite_B = 5000
  )
  data <- structure(list(K = 1, p = p, n = 50, nm1 = 49,
                         Xty = rnorm(p), eigen_R = NULL),
                    class = "individual")
  params <- list(R_mismatch = "eb")

  result <- compute_R_mismatch_state(data, params, model)
  expect_null(result$lambda_bias)
})

test_that("compute_R_mismatch_state short-circuits when R_mismatch is 'none'", {
  data <- structure(list(K = 1), class = c("ss"))
  model <- list(sigma2 = 1)
  expect_identical(compute_R_mismatch_state(data, list(), model), model)
})

test_that("compute_R_mismatch_state stops when eigen_R is unavailable", {
  set.seed(901)
  p <- 6; n <- 100
  z <- rnorm(p)
  R <- diag(p)

  data <- structure(
    list(
      Xty   = z * sqrt(n - 1),
      nm1   = n - 1,
      z     = z,
      p     = p,
      n     = n,
      lambda = 0,
      X_colmeans = rep(0, p),
      y_mean     = 0,
      XtX        = (n - 1) * R,
      R_finite_B = Inf,
      R_mismatch = "eb",
      eigen_R    = NULL,
      X          = NULL
    ),
    class = c("ss", "list")
  )

  params <- list(
    R_mismatch         = "eb",
    R_mismatch_method  = "mle",
    eig_delta_rel      = 1e-3,
    eig_delta_abs      = 0,
    artifact_threshold = 0.1
  )

  model <- list(
    sigma2   = 1,
    lambda_bias = 0,
    converged   = FALSE,
    eigen_R     = NULL,
    alpha       = matrix(1 / p, 3, p),
    mu          = matrix(0, 3, p),
    XtXr        = rep(0, p)
  )

  expect_error(
    compute_R_mismatch_state(data, params, model),
    "R_mismatch requires data\\$eigen_R"
  )
})

# ---- summarize_R_bf_attenuation (cs_index name-parsing + early returns) ----

test_that("summarize_R_bf_attenuation parses cs_index from names when absent", {
  model <- list(
    alpha = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE),
    R_bf_attenuation = matrix(c(log(25), 0, 0, log(2)),
                              nrow = 2, byrow = TRUE),
    sets = list(cs = list(L1 = 1L, L2 = 2L)),
    R_finite_diagnostics = list(artifact_flag = FALSE)
  )

  out <- summarize_R_bf_attenuation(model, threshold = log(20))
  d <- out$R_finite_diagnostics
  expect_true(d$R_sensitivity_flag)
  expect_equal(d$bf_attenuation$cs_label[["L1"]], "sensitive")
  expect_equal(d$bf_attenuation$cs_label[["L2"]], "stable")
})

test_that("summarize_R_bf_attenuation returns the model unchanged on missing inputs", {
  m1 <- list(R_finite_diagnostics = list(artifact_flag = FALSE))
  expect_identical(summarize_R_bf_attenuation(m1), m1)

  m2 <- list(R_bf_attenuation = matrix(0, 1, 1))
  expect_identical(summarize_R_bf_attenuation(m2), m2)
})

# ---- Edge cases ----

test_that("perfectly matched R (identity) yields lambda_bias ~ 0", {
  set.seed(400)
  p <- 10
  n <- 500
  R <- diag(p)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 2,
                   R_mismatch = "eb_no_init", max_iter = 5, verbose = FALSE)
  expect_equal(fit$R_finite_diagnostics$lambda_bias, 0)
})

test_that("R_mismatch = 'none' skips all mismatch machinery", {
  set.seed(401)
  p <- 10
  n <- 500
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 2,
                   R_mismatch = "none", max_iter = 5, verbose = FALSE)
  expect_null(fit$R_finite_diagnostics)
})
