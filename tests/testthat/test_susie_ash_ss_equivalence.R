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
