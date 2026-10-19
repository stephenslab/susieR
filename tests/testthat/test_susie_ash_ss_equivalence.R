# Helper: prepare individual and SS data from a common dataset
setup_susie_ash_test <- function(n = 200, p = 50, k = 5, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = TRUE)
  beta_true <- rep(0, p)
  causal <- sample(1:p, k)
  beta_true[causal] <- rnorm(k, sd = 2)
  y <- c(X %*% beta_true + rnorm(n))
  y <- y - mean(y)

  # Sufficient stats
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  yty <- sum(y^2)

  # RSS inputs
  bhat <- sapply(1:p, function(j) sum(X[, j] * y) / sum(X[, j]^2))
  shat <- sapply(1:p, function(j) {
    resid <- y - X[, j] * bhat[j]
    sqrt(sum(resid^2) / ((n - 2) * sum(X[, j]^2)))
  })
  R_mat <- susieR:::safe_cor(X)

  list(X = X, y = y, n = n, p = p,
       XtX = XtX, Xty = Xty, yty = yty,
       bhat = bhat, shat = shat, R_mat = R_mat,
       beta_true = beta_true, causal = causal)
}

test_that("susie_ss ash agrees with susie individual-level ash", {
  d <- setup_susie_ash_test(n = 200, p = 50, k = 5, seed = 42)

  fit_ind <- susie(d$X, d$y, L = 5,
    unmappable_effects = "ash_filter_archived",
    estimate_residual_variance = TRUE,
    estimate_prior_method = "optim",
    intercept = FALSE, standardize = FALSE,
    max_iter = 20, verbose = FALSE
  )

  fit_ss <- susie_ss(
    XtX = d$XtX, Xty = d$Xty, yty = d$yty, n = d$n, L = 5,
    unmappable_effects = "ash_filter_archived",
    estimate_residual_variance = TRUE,
    estimate_prior_method = "optim",
    max_iter = 20, verbose = FALSE
  )

  expect_equal(fit_ss$theta, fit_ind$theta, tolerance = 1e-4,
    label = "theta (mr.ash coefficients)")
  expect_equal(fit_ss$sigma2, fit_ind$sigma2, tolerance = 1e-4,
    label = "sigma2 (residual variance)")
  expect_true(cor(susie_get_pip(fit_ss), susie_get_pip(fit_ind)) > 0.999,
    label = "PIP correlation > 0.999")
})

test_that("susie_ss ash agrees with susie individual-level ash across data sizes", {
  for (params in list(
    list(n = 100, p = 30, k = 3, seed = 100),
    list(n = 300, p = 80, k = 8, seed = 200)
  )) {
    d <- setup_susie_ash_test(
      n = params$n, p = params$p,
      k = params$k, seed = params$seed
    )

    fit_ind <- suppressWarnings(susie(d$X, d$y, L = 5,
      unmappable_effects = "ash_filter_archived",
      estimate_residual_variance = TRUE,
      estimate_prior_method = "optim",
      intercept = FALSE, standardize = FALSE,
      max_iter = 15, verbose = FALSE
    ))

    fit_ss <- susie_ss(
      XtX = d$XtX, Xty = d$Xty, yty = d$yty, n = d$n, L = 5,
      unmappable_effects = "ash_filter_archived",
      estimate_residual_variance = TRUE,
      estimate_prior_method = "optim",
      max_iter = 15, verbose = FALSE
    )

    expect_equal(fit_ss$theta, fit_ind$theta, tolerance = 1e-4,
      label = sprintf("theta (n=%d, p=%d)", params$n, params$p))
    expect_equal(fit_ss$sigma2, fit_ind$sigma2, tolerance = 1e-4,
      label = sprintf("sigma2 (n=%d, p=%d)", params$n, params$p))
    expect_true(cor(susie_get_pip(fit_ss), susie_get_pip(fit_ind)) > 0.999,
      label = sprintf("PIP cor > 0.999 (n=%d, p=%d)", params$n, params$p))
  }
})

# ---- output structure ----

test_that("susie individual-level ash output has expected fields and dimensions", {
  d <- setup_susie_ash_test(n = 100, p = 30, k = 3, seed = 123)

  fit_ind <- susie(d$X, d$y, L = 5,
    unmappable_effects = "ash_filter_archived",
    estimate_residual_variance = TRUE,
    estimate_prior_method = "optim",
    intercept = FALSE, standardize = FALSE,
    max_iter = 10, verbose = FALSE
  )

  expect_true(is.numeric(fit_ind$theta))
  expect_true(is.numeric(fit_ind$sigma2))
  expect_true(is.numeric(fit_ind$tau2))
  expect_true(is.matrix(fit_ind$alpha))
  expect_length(fit_ind$theta, d$p)
  expect_length(fit_ind$sigma2, 1)
  expect_equal(ncol(fit_ind$alpha), d$p)
  expect_null(fit_ind$X_theta)
})

test_that("susie_ss ash output has expected fields and dimensions", {
  d <- setup_susie_ash_test(n = 100, p = 30, k = 3, seed = 123)

  fit_ss <- susie_ss(
    XtX = d$XtX, Xty = d$Xty, yty = d$yty, n = d$n, L = 5,
    unmappable_effects = "ash_filter_archived",
    estimate_residual_variance = TRUE,
    estimate_prior_method = "optim",
    max_iter = 10, verbose = FALSE
  )

  expect_true(is.numeric(fit_ss$theta))
  expect_true(is.numeric(fit_ss$sigma2))
  expect_true(is.numeric(fit_ss$tau2))
  expect_true(is.matrix(fit_ss$alpha))
  expect_length(fit_ss$theta, d$p)
  expect_length(fit_ss$sigma2, 1)
  expect_equal(ncol(fit_ss$alpha), d$p)
})

test_that("susie_rss ash works with correlation matrix input", {
  d <- setup_susie_ash_test(n = 200, p = 50, k = 5, seed = 42)

  fit_rss <- susie_rss(
    bhat = d$bhat, shat = d$shat, R = d$R_mat, n = d$n, L = 5,
    unmappable_effects = "ash_filter_archived",
    estimate_residual_variance = TRUE,
    estimate_prior_method = "optim",
    max_iter = 20, verbose = FALSE
  )

  expect_true(is.numeric(fit_rss$theta))
  expect_length(fit_rss$theta, d$p)
  expect_true(is.numeric(fit_rss$sigma2))
  expect_true(length(fit_rss$sets$cs) >= 0)
})

# ---- unstandardized X (individual-level) ----

# susie() keeps X on its original scale and stores the centering/scaling as
# attributes. The tests above pass an already-standardized X, so here X is
# given as uncentered 0/1/2 dosages (plus a copy with heterogeneous column
# scales), with a polygenic background so that theta is non-zero.
.make_raw_dosage_ash_data <- function()
  simulate_dosage_regression(n = 400, p = 60, seed = 2024, causal = c(10, 40),
                             effects = c(0.6, -0.5))

.raw_dosage_fit_args <- list(L = 5, estimate_residual_variance = TRUE,
                             max_iter = 20, verbose = FALSE)

test_that("individual-level ash does not depend on the scaling of X (ash and archived paths)", {
  d <- .make_raw_dosage_ash_data()
  # Same standardized X, column scales spanning exp(-3) to exp(3).
  X_het <- sweep(d$X, 2, exp(seq(-3, 3, length.out = d$p)), "*")

  for (method in c("ash", "ash_filter_archived")) {
    args <- c(.raw_dosage_fit_args, list(unmappable_effects = method))
    fit_raw <- suppressMessages(do.call(susie, c(list(X = d$X, y = d$y), args)))
    fit_std <- suppressMessages(do.call(susie, c(list(X = scale(d$X), y = d$y), args)))
    fit_het <- suppressMessages(do.call(susie, c(list(X = X_het, y = d$y), args)))

    expect_true(any(fit_std$theta != 0), label = paste(method, "theta is non-zero"))
    # Compare converged fits, so differences are not an artefact of max_iter.
    expect_true(fit_raw$converged, label = paste(method, "raw fit converged"))

    for (nm in c("theta", "sigma2", "tau2", "pip", "fitted")) {
      expect_equal(fit_raw[[nm]], fit_std[[nm]], tolerance = 1e-6,
        label = sprintf("%s from unstandardized X (%s)", nm, method))
      expect_equal(fit_het[[nm]], fit_std[[nm]], tolerance = 1e-6,
        label = sprintf("%s from heterogeneously scaled X (%s)", nm, method))
    }

    # coef() maps the standardized-scale effects (including theta) back to
    # the scale of X.
    expect_equal(coef(fit_raw)[-1],
                 coef(fit_std)[-1] / attr(scale(d$X), "scaled:scale"),
                 tolerance = 1e-6,
                 label = sprintf("coef from unstandardized X (%s)", method))
    expect_equal(coef(fit_het)[-1],
                 coef(fit_std)[-1] / attr(scale(X_het), "scaled:scale"),
                 tolerance = 1e-6,
                 label = sprintf("coef from heterogeneously scaled X (%s)", method))
  }
})

test_that("individual-level ash on unstandardized X agrees with susie_ss", {
  d <- .make_raw_dosage_ash_data()
  args <- c(.raw_dosage_fit_args, list(unmappable_effects = "ash"))
  fit_raw <- suppressMessages(do.call(susie, c(list(X = d$X, y = d$y), args)))

  Xs <- scale(d$X)
  yc <- d$y - mean(d$y)
  fit_ss <- suppressMessages(do.call(susie_ss,
    c(list(XtX = crossprod(Xs), Xty = c(crossprod(Xs, yc)), yty = sum(yc^2), n = d$n),
      args)))

  # mr.ash and mr.ash.rss stop at slightly different points, so theta agrees
  # in absolute (not relative) terms; PIPs agree much more tightly.
  expect_lt(max(abs(fit_raw$theta - fit_ss$theta)), 5e-4)
  expect_equal(fit_raw$sigma2, fit_ss$sigma2, tolerance = 1e-3,
    label = "sigma2 vs susie_ss")
  expect_equal(susie_get_pip(fit_raw), susie_get_pip(fit_ss), tolerance = 1e-5,
    label = "PIPs vs susie_ss")
})

test_that("susie_auto with ash does not depend on the scaling of X", {
  d <- .make_raw_dosage_ash_data()
  # mr.ash warns about its default prior-variance grid in the small-L stages.
  fit_raw <- suppressWarnings(suppressMessages(
    susie_auto(d$X, d$y, unmappable_effects = "ash", verbose = FALSE)))
  fit_std <- suppressWarnings(suppressMessages(
    susie_auto(scale(d$X), d$y, unmappable_effects = "ash", verbose = FALSE)))

  expect_lte(nrow(fit_raw$alpha), 8)
  expect_equal(nrow(fit_raw$alpha), nrow(fit_std$alpha))
  expect_equal(fit_raw$pip, fit_std$pip, tolerance = 1e-6)
})
