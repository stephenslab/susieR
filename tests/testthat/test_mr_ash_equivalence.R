context("mr.ash vs mr.ash.rss equivalence")

# Verifies that mr.ash (individual-level data) and mr.ash.rss (summary
# statistics) agree when fed the same data, and exercises the mr.ash.rss
# argument-defaulting branches and the C++ kernels in mr_ash_rss.h.
#
# Summary statistics are derived from individual data via derive_summary_stats
# (helper_mr_ash.R): bhat_j = X_j'y / X_j'X_j, shat_j with n-2 df, R = cor(X),
# var_y = var(y). Test data comes from mr_ash_data() (helper_mr_ash.R), whose
# construction matches the priors the equivalence tolerances were tuned to.

# ---- mr.ash <-> mr.ash.rss numeric equivalence -----------------------------

test_that("mr.ash and mr.ash.rss produce identical beta with fixed sigma and pi", {
  d <- mr_ash_data(n = 100, p = 50, k = 5, seed = 42)

  fit_ind <- mr.ash(d$X, d$y,
    sa2 = d$sa2, pi = d$pi0, sigma2 = d$sigma2_init,
    intercept = FALSE, standardize = FALSE,
    update.sigma2 = FALSE, update.pi = FALSE,
    max.iter = 100, verbose = FALSE)

  ss <- derive_summary_stats(d$X, d$y)
  fit_rss <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 100, update_w0 = FALSE, update_sigma = FALSE)

  # Match to near-machine precision (drop Armadillo dim attributes via c()).
  expect_equal(c(fit_rss$beta), c(fit_ind$beta), tolerance = 1e-10)
  expect_equal(c(fit_rss$sigma2), c(fit_ind$sigma2), tolerance = 1e-10)
  expect_equal(c(fit_rss$pi), c(fit_ind$pi), tolerance = 1e-10)
})

test_that("mr.ash and mr.ash.rss agree with sigma2 updates enabled", {
  d <- mr_ash_data(n = 100, p = 50, k = 5, seed = 42)

  fit_ind <- mr.ash(d$X, d$y,
    sa2 = d$sa2, pi = d$pi0, sigma2 = d$sigma2_init,
    intercept = FALSE, standardize = FALSE,
    update.sigma2 = TRUE, update.pi = FALSE,
    max.iter = 200, verbose = FALSE)

  ss <- derive_summary_stats(d$X, d$y)
  fit_rss <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 200, update_w0 = FALSE, update_sigma = TRUE)

  expect_equal(c(fit_rss$beta), c(fit_ind$beta), tolerance = 1e-3)
  expect_equal(c(fit_rss$sigma2), c(fit_ind$sigma2), tolerance = 1e-4)
})

test_that("mr.ash and mr.ash.rss agree with full EM (sigma + pi updates)", {
  d <- mr_ash_data(n = 100, p = 50, k = 5, seed = 42)

  fit_ind <- mr.ash(d$X, d$y,
    sa2 = d$sa2, pi = d$pi0, sigma2 = d$sigma2_init,
    intercept = FALSE, standardize = FALSE,
    update.sigma2 = TRUE, update.pi = TRUE,
    max.iter = 200, verbose = FALSE)

  ss <- derive_summary_stats(d$X, d$y)
  fit_rss <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 200, update_w0 = TRUE, update_sigma = TRUE)

  expect_equal(c(fit_rss$beta), c(fit_ind$beta), tolerance = 1e-3)
  expect_equal(c(fit_rss$sigma2), c(fit_ind$sigma2), tolerance = 1e-3)
  expect_equal(c(fit_rss$pi), c(fit_ind$pi), tolerance = 1e-2)
})

test_that("mr.ash and mr.ash.rss agree across different data sizes", {
  for (params in list(
    list(n = 80,  p = 20, k = 3, seed = 100),
    list(n = 200, p = 30, k = 5, seed = 200))) {
    d <- mr_ash_data(n = params$n, p = params$p, k = params$k,
                     seed = params$seed)

    fit_ind <- mr.ash(d$X, d$y,
      sa2 = d$sa2, pi = d$pi0, sigma2 = d$sigma2_init,
      intercept = FALSE, standardize = FALSE,
      update.sigma2 = FALSE, update.pi = FALSE,
      max.iter = 50, verbose = FALSE)

    ss <- derive_summary_stats(d$X, d$y)
    fit_rss <- mr.ash.rss(
      bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
      sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
      tol = 1e-4, max_iter = 50, update_w0 = FALSE, update_sigma = FALSE)

    expect_equal(c(fit_rss$beta), c(fit_ind$beta), tolerance = 1e-10,
                 label = sprintf("beta (n=%d, p=%d)", params$n, params$p))
  }
})

# ---- mr.ash.rss output format ----------------------------------------------

test_that("mr.ash.rss returns mr.ash-compatible fields plus the RSS-specific ones", {
  d  <- mr_ash_data(n = 80, p = 20, k = 3, seed = 123)
  ss <- derive_summary_stats(d$X, d$y)
  fit_rss <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
    sigma2_e = 1.0, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 100, update_w0 = FALSE, update_sigma = FALSE)

  # mr.ash-compatible fields: correct types and dimensions.
  expect_true(is.numeric(fit_rss$beta))
  expect_true(is.numeric(fit_rss$sigma2))
  expect_true(is.numeric(fit_rss$pi))
  expect_true(is.integer(fit_rss$iter))
  expect_true(is.numeric(fit_rss$varobj))
  expect_length(fit_rss$beta, d$p)
  expect_length(fit_rss$sigma2, 1)
  expect_length(fit_rss$pi, d$K)
  expect_true(fit_rss$iter > 0)
  expect_true(length(fit_rss$varobj) > 0 && length(fit_rss$varobj) <= 100)

  # RSS-specific fields present.
  for (nm in c("mu1", "sigma2_1", "w1", "sigma2_e", "w0")) {
    expect_false(is.null(fit_rss[[nm]]), info = nm)
  }
})

# ---- Argument-defaulting branches in mr.ash.rss ----------------------------

test_that("mr.ash.rss defaults sigma2_e from var_y (or 1 when var_y is NULL/Inf)", {
  # var_y = NULL -> Inf, and sigma2_e = NULL -> defaults to 1.
  d1 <- mr_ash_data(n = 80, p = 16, k = 3, seed = 601)
  ss1 <- derive_summary_stats(d1$X, d1$y)
  fit_inf <- mr.ash.rss(
    bhat = ss1$bhat, shat = ss1$shat, R = ss1$R,
    var_y = NULL, n = ss1$n, sigma2_e = NULL,
    s0 = d1$sa2, w0 = d1$pi0,
    tol = 1e-4, max_iter = 50, update_w0 = FALSE, update_sigma = FALSE)
  expect_length(fit_inf$beta, d1$p)
  expect_true(is.finite(fit_inf$sigma2))

  # Finite var_y with sigma2_e = NULL -> sigma2_e initialised to var_y.
  d2 <- mr_ash_data(n = 80, p = 16, k = 3, seed = 602)
  ss2 <- derive_summary_stats(d2$X, d2$y)
  fit_var <- mr.ash.rss(
    bhat = ss2$bhat, shat = ss2$shat, R = ss2$R,
    var_y = ss2$var_y, n = ss2$n, sigma2_e = NULL,
    s0 = d2$sa2, w0 = d2$pi0,
    tol = 1e-4, max_iter = 50, update_w0 = FALSE, update_sigma = FALSE)
  expect_length(fit_var$beta, d2$p)
  expect_true(is.finite(fit_var$sigma2) && fit_var$sigma2 > 0)
})

test_that("mr.ash.rss accepts an explicit z and derives z = bhat/shat otherwise", {
  d  <- mr_ash_data(n = 100, p = 16, k = 3, seed = 603)
  ss <- derive_summary_stats(d$X, d$y)

  # Explicit z bypasses the bhat/shat derivation; default z = numeric(0) makes
  # the C++ derive z_use = bhat/shat. Both paths must return a valid fit and,
  # with compute_ELBO = TRUE, a finite varobj trace.
  for (z_arg in list(ss$bhat / ss$shat, numeric(0))) {
    fit <- mr.ash.rss(
      bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
      sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
      z = z_arg,
      tol = 1e-4, max_iter = 50,
      update_w0 = TRUE, update_sigma = TRUE, compute_ELBO = TRUE)
    expect_length(fit$beta, d$p)
    expect_true(all(is.finite(fit$beta)))
    expect_true(length(fit$varobj) > 0 && all(is.finite(fit$varobj)))
  }
})

# ---- compute_ELBO toggle ---------------------------------------------------

test_that("compute_ELBO controls the varobj trace and leaves beta finite either way", {
  d  <- mr_ash_data(n = 80, p = 16, k = 3, seed = 604)
  ss <- derive_summary_stats(d$X, d$y)
  args <- list(
    bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 50, update_w0 = FALSE, update_sigma = FALSE)

  fit_on  <- do.call(mr.ash.rss, c(args, list(compute_ELBO = TRUE)))
  fit_off <- do.call(mr.ash.rss, c(args, list(compute_ELBO = FALSE)))

  expect_length(fit_on$beta, d$p)
  expect_true(all(is.finite(fit_on$beta)))
  expect_true(length(fit_on$varobj) > 0 && all(is.finite(fit_on$varobj)))

  expect_length(fit_off$beta, d$p)
  expect_true(all(is.finite(fit_off$beta)))
  expect_true(is.numeric(fit_off$varobj))
})

# ---- C++ kernel paths (mr_ash_rss.h) ---------------------------------------

test_that("mr.ash.rss hits the non-convergence path when max_iter is too small", {
  # tol is tight and max_iter = 1, so the iteration cap is reached (the C++
  # prints a non-convergence note to stderr) and iter == max_iter.
  d  <- mr_ash_data(n = 100, p = 20, k = 4, seed = 605)
  ss <- derive_summary_stats(d$X, d$y)
  fit <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-12, max_iter = 1,
    update_w0 = TRUE, update_sigma = TRUE, compute_ELBO = TRUE)
  expect_length(fit$beta, d$p)
  expect_equal(as.integer(fit$iter), 1L)
})

test_that("mr.ash.rss with standardize = TRUE returns a valid, rescaled fit", {
  # Exercises the XtX standardization + elementwise rescale-back path. sigma2_1
  # is a length-p variance vector (not a p x p covariance), so the rescale is
  # elementwise and standardize = TRUE returns a valid fit.
  d  <- mr_ash_data(n = 100, p = 16, k = 3, seed = 611)
  ss <- derive_summary_stats(d$X, d$y)
  fit <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 100,
    update_w0 = TRUE, update_sigma = TRUE, compute_ELBO = TRUE,
    standardize = TRUE)
  expect_length(c(fit$beta), d$p)
  expect_true(all(is.finite(c(fit$beta))))
  expect_length(c(fit$sigma2_1), d$p)
  expect_true(all(is.finite(c(fit$sigma2_1))) && all(c(fit$sigma2_1) >= 0))
  expect_true(is.finite(fit$sigma2) && fit$sigma2 > 0)
})

test_that("mr.ash.rss constructs XtX/Xty from finite var_y on the raw scale", {
  # standardize = FALSE + finite var_y exercises the XtXdiag construction.
  d  <- mr_ash_data(n = 120, p = 16, k = 2, seed = 606)
  ss <- derive_summary_stats(d$X, d$y)
  fit <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 80,
    update_w0 = TRUE, update_sigma = TRUE, compute_ELBO = TRUE,
    standardize = FALSE)
  expect_length(fit$beta, d$p)
  expect_true(all(is.finite(fit$beta)))
  expect_true(is.finite(fit$sigma2) && fit$sigma2 > 0)
})

test_that("mr.ash.rss uses the standardized-scale XtX = (n-1)*R when var_y is Inf", {
  # var_y = NULL -> Inf selects the standardized-scale XtX path; also covers
  # the update_w0 = FALSE + update_sigma = FALSE combination.
  d  <- mr_ash_data(n = 100, p = 16, k = 3, seed = 609)
  ss <- derive_summary_stats(d$X, d$y)
  fit <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = NULL, n = ss$n,
    sigma2_e = NULL, s0 = d$sa2, w0 = d$pi0,
    tol = 1e-4, max_iter = 80,
    update_w0 = FALSE, update_sigma = FALSE, compute_ELBO = FALSE,
    standardize = FALSE)
  expect_length(fit$beta, d$p)
  expect_true(all(is.finite(fit$beta)))
})

test_that("mr.ash.rss with a single-component prior falls back to bayes_ridge", {
  # K = 1 -> bayes_mix_sufficient calls bayes_ridge_sufficient once.
  d  <- mr_ash_data(n = 100, p = 16, k = 2, seed = 610)
  ss <- derive_summary_stats(d$X, d$y)
  fit <- mr.ash.rss(
    bhat = ss$bhat, shat = ss$shat, R = ss$R, var_y = ss$var_y, n = ss$n,
    sigma2_e = d$sigma2_init, s0 = c(1.0), w0 = c(1.0),
    tol = 1e-4, max_iter = 50,
    update_w0 = FALSE, update_sigma = FALSE, compute_ELBO = TRUE)
  expect_length(fit$beta, d$p)
  expect_true(all(is.finite(fit$beta)))
  expect_true(is.numeric(fit$varobj))
})
