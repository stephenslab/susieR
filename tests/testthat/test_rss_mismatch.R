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

test_that("R_mismatch = 'eb' runs and returns Q_art diagnostics", {
  set.seed(11)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- suppressWarnings(susie_rss(z = z, R = R, n = n, L = 3, R_finite = 5000,
                   R_mismatch = "eb", max_iter = 2, verbose = FALSE))
  d <- fit$R_finite_diagnostics
  expect_true(!is.null(d$Q_art))
  expect_true(d$Q_art >= 0 && d$Q_art <= 1)
  expect_true(is.logical(d$artifact_flag))
  expect_true(d$mode_label %in% c("normal", "warning", "conservative"))
})

test_that("R_mismatch = 'eb' supports the B = Inf limit", {
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
  expect_true(d$lambda_bias >= 0)
  expect_true(!is.null(d$Q_art))
})

test_that("MLE lambda estimator handles zero, interior, and large fits", {
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

test_that("R_mismatch_method controls lambda estimator on SS path only", {
  rho <- 0.99
  R <- matrix(c(1, -rho, -rho, 1), 2, 2)
  z <- c(-8, -8)

  fit_mle <- suppressWarnings(susie_rss(
    z = z, R = R, n = 5000, L = 1, R_mismatch = "eb_no_init",
    R_mismatch_method = "mle", max_iter = 3, verbose = FALSE))
  d_mle <- fit_mle$R_finite_diagnostics
  expect_equal(d_mle$R_mismatch_method, "mle")
  expect_true(d_mle$lambda_bias > 1)
  expect_true(d_mle$B_corrected < 1)

  fit_map <- suppressWarnings(susie_rss(
    z = z, R = R, n = 5000, L = 1, R_mismatch = "eb_no_init",
    R_mismatch_method = "map", max_iter = 3, verbose = FALSE))
  d_map <- fit_map$R_finite_diagnostics
  expect_equal(d_map$R_mismatch_method, "map")
  expect_true(d_map$B_corrected < 1)
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
  expect_true(!is.null(d$R_mismatch_init))
  expect_equal(d$R_mismatch_init$method, "ser")
  expect_true(!is.null(d$R_mismatch_trace))
  expect_equal(d$R_mismatch_trace[[1]]$phase, "init_ser")
})

test_that("R_mismatch = 'eb' skips SER-protected initialization with finite R_finite", {
  set.seed(14)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- suppressWarnings(susie_rss(z = z, R = R, n = n, L = 3,
                   R_finite = 5000, R_mismatch = "eb", max_iter = 2,
                   track_fit = TRUE, verbose = FALSE))
  d <- fit$R_finite_diagnostics
  expect_true(is.null(d$R_mismatch_init))
  expect_true(!is.null(d$Q_art))
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
  expect_true(!is.null(d$R_mismatch_init))
  expect_equal(d$R_mismatch_trace[[1]]$phase, "init_ser")
})

test_that("R_mismatch = 'eb_adaptive_init' tempers SER initialization", {
  set.seed(19)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3,
                   R_finite = 5000, R_mismatch = "eb_adaptive_init",
                   max_iter = 2, track_fit = TRUE, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  init <- d$R_mismatch_init
  expect_true(!is.null(init))
  expect_true(!is.null(d$Q_art))
  expect_equal(d$R_mismatch_trace[[1]]$phase, "init_ser")
  expect_true(is.finite(init$ld_coherence))
  expect_true(init$ld_coherence >= 0 && init$ld_coherence <= 1)
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
  expect_true(is.null(d$R_mismatch_init))
  expect_true(!is.null(d$Q_art))
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
  expect_true(out$R_bf_attenuation[1, 1] > 0)
  expect_true(out$R_bf_attenuation[1, 2] >= 0)
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
  expect_true(fit$R_finite_diagnostics$lambda_bias >= 0)
  expect_equal(fit$R_finite_diagnostics$B_corrected,
               1 / (1 / fit$R_finite_diagnostics$B +
                      fit$R_finite_diagnostics$lambda_bias),
               tolerance = 1e-12)
})

# ---- Q_art unit tests ----

test_that("compute_Q_art recovers Q ~ 1 when r_fit lies in low-eigen direction", {
  # Diagonal R with eigenvalues (2, 1, 1e-6). Default eig_delta_rel=1e-3
  # selects only the third eigenvalue.
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  r_fit <- c(0, 0, 1)  # purely in the low-eigen direction
  out <- susieR:::compute_Q_art(eig, r_fit)
  expect_equal(out$Q_art, 1, tolerance = 1e-12)
  expect_true(out$evaluable)
  expect_equal(out$low_eigen_count, 1L)
})

test_that("compute_Q_art returns Q ~ 0 when r_fit avoids low-eigen directions", {
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  r_fit <- c(1, 0.5, 0)  # fully in top two eigen directions
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
  d <- c(2, 1, 0.5)  # all > 1e-3 * 2 = 2e-3
  eig <- list(values = d, vectors = V)
  out <- susieR:::compute_Q_art(eig, c(1, 0, 0))
  expect_equal(out$low_eigen_count, 0L)
  expect_false(out$evaluable)
})

test_that("compute_Q_art is in [0, 1] for typical inputs", {
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  for (r_fit in list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1),
                     c(0.5, 0.5, 0.5), c(-1, 1, -1))) {
    out <- susieR:::compute_Q_art(eig, r_fit)
    expect_true(out$Q_art >= 0 && out$Q_art <= 1)
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

test_that("eb surfaces Q_art and mode_label diagnostics", {
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
    expect_true(!is.null(d[[fld]]), info = paste("missing diagnostic:", fld))
})

test_that("eb with X-input runs and surfaces Q_art", {
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
  expect_true(!is.null(d$Q_art))
  expect_true(d$Q_art >= 0 && d$Q_art <= 1)
})

test_that("eb works on lambda=0 multi-panel SS path", {
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
  expect_true(!is.null(d$Q_art))
  expect_true(d$Q_art >= 0 && d$Q_art <= 1)
  expect_length(d$lambda_bias, 1)
})

test_that("eb_adaptive_init handles multi-panel X input", {
  set.seed(21)
  n <- 80
  p <- 12
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p), n, p)
  z <- rnorm(p)

  fit <- susie_rss(z = z, X = list(X1, X2), n = 1000, L = 3,
                   R_finite = TRUE, R_mismatch = "eb_adaptive_init",
                   max_iter = 3, verbose = FALSE)
  d <- fit$R_finite_diagnostics
  expect_true(!is.null(d$R_mismatch_init))
  expect_true(is.finite(d$R_mismatch_init$ld_coherence))
  expect_true(d$R_mismatch_init$ld_coherence >= 0)
  expect_true(d$R_mismatch_init$ld_coherence <= 1)
  expect_length(d$lambda_bias, 1)
  expect_length(d$B_corrected, 1)
})
