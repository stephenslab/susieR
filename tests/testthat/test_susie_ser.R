context("susie_ser")

machine_tol <- .Machine$double.eps

# ---- core computation ----

test_that("susie_ser simple method matches direct ABF calculation", {
  z <- c(-0.5, 4, 0.3, 1.1)
  prior_weights <- c(0.1, 0.6, 0.2, 0.1)
  V <- 2

  fit <- susie_ser(z = z, prior_variance = V,
                   prior_weights = prior_weights,
                   estimate_prior_method = "simple",
                   coverage = NULL)

  lbf <- gaussian_ser_lbf(z, rep(1, length(z)), V)
  lpo <- lbf + log(prior_weights + sqrt(.Machine$double.eps))
  alpha <- exp(lpo - max(lpo))
  alpha <- alpha / sum(alpha)
  moments <- gaussian_ser_moments(z, rep(1, length(z)), V)

  expect_equal(as.numeric(fit$alpha), alpha)
  expect_equal(as.numeric(fit$mu), moments$post_mean)
  expect_equal(as.numeric(fit$mu2), moments$post_mean2)
  expect_equal(fit$V, V)
  expect_equal(fit$niter, 1)
  expect_true(fit$converged)
})

test_that("susie_ser is a single-effect model: niter=1 and converged", {
  set.seed(10)
  fit <- susie_ser(z = c(0, 5, 0, 0), estimate_prior_method = "simple",
                   coverage = NULL)
  expect_equal(fit$niter, 1L)
  expect_true(fit$converged)
  expect_equal(nrow(fit$alpha), 1L)
})

test_that("susie_ser only accepts optim and simple prior methods", {
  expect_error(
    susie_ser(z = rnorm(5), estimate_prior_method = "EM", coverage = NULL),
    "should be one of"
  )
})

# ---- PVE adjustment and z storage ----

test_that("susie_ser applies PVE adjustment when n is provided", {
  z <- c(8, 0.5, -1)
  n <- 100
  fit <- susie_ser(z = z, n = n, estimate_prior_method = "simple",
                   coverage = NULL)
  adj <- (n - 1) / (z^2 + n - 2)
  expect_equal(fit$pve_adjustment, adj)
  expect_equal(fit$z, z)
})

# ---- agreement with susie_rss on diagonal LD ----

test_that("susie_ser with z and n matches susie_rss with diagonal R", {
  set.seed(1)
  p <- 30
  n <- 500
  z <- rnorm(p)
  z[7] <- 6

  fit_ser <- susie_ser(z = z, n = n, estimate_prior_method = "optim",
                       coverage = NULL)
  fit_rss <- suppressWarnings(suppressMessages(
    susie_rss(z = z, R = diag(p), n = n, L = 1, max_iter = 1,
              estimate_prior_method = "optim", coverage = NULL,
              check_prior = FALSE)
  ))

  expect_equal(fit_ser$alpha, fit_rss$alpha, tolerance = machine_tol)
  expect_equal(fit_ser$mu, fit_rss$mu, tolerance = machine_tol)
  expect_equal(fit_ser$mu2, fit_rss$mu2, tolerance = machine_tol)
  expect_equal(fit_ser$lbf, fit_rss$lbf, tolerance = machine_tol)
  expect_equal(fit_ser$lbf_variable, fit_rss$lbf_variable,
               tolerance = machine_tol)
  expect_equal(fit_ser$V, fit_rss$V, tolerance = machine_tol)
  expect_equal(fit_ser$pip, fit_rss$pip, tolerance = machine_tol)
})

test_that("susie_ser simple method matches diagonal susie_rss", {
  set.seed(2)
  p <- 25
  n <- 400
  z <- rnorm(p)
  z[4] <- 5

  fit_ser <- susie_ser(z = z, n = n, estimate_prior_method = "simple",
                       scaled_prior_variance = 0.4, coverage = NULL)
  fit_rss <- suppressWarnings(suppressMessages(
    susie_rss(z = z, R = diag(p), n = n, L = 1, max_iter = 1,
              estimate_prior_method = "simple",
              scaled_prior_variance = 0.4,
              coverage = NULL, check_prior = FALSE)
  ))

  expect_equal(fit_ser$alpha, fit_rss$alpha, tolerance = machine_tol)
  expect_equal(fit_ser$mu, fit_rss$mu, tolerance = machine_tol)
  expect_equal(fit_ser$mu2, fit_rss$mu2, tolerance = machine_tol)
  expect_equal(fit_ser$V, fit_rss$V, tolerance = machine_tol)
})

test_that("susie_ser z-only scale matches diagonal susie_rss without n", {
  set.seed(3)
  p <- 20
  z <- rnorm(p)
  z[2] <- 4.5

  fit_ser <- susie_ser(z = z, estimate_prior_method = "optim",
                       coverage = NULL)
  fit_rss <- suppressWarnings(suppressMessages(
    susie_rss(z = z, R = diag(p), L = 1, max_iter = 1,
              estimate_prior_method = "optim", coverage = NULL,
              check_prior = FALSE)
  ))

  expect_equal(fit_ser$alpha, fit_rss$alpha, tolerance = machine_tol)
  expect_equal(fit_ser$mu, fit_rss$mu, tolerance = machine_tol)
  expect_equal(fit_ser$mu2, fit_rss$mu2, tolerance = machine_tol)
  expect_equal(fit_ser$V, fit_rss$V, tolerance = machine_tol)
})

test_that("susie_ser with z, n, and var_y matches diagonal susie_rss", {
  set.seed(6)
  p <- 22
  n <- 450
  var_y <- 2.3
  z <- rnorm(p)
  z[11] <- 5.2

  fit_ser <- susie_ser(z = z, n = n, var_y = var_y,
                       estimate_prior_method = "optim", coverage = NULL)
  fit_rss <- suppressWarnings(suppressMessages(
    susie_rss(z = z, R = diag(p), n = n, var_y = var_y,
              L = 1, max_iter = 1, estimate_prior_method = "optim",
              coverage = NULL, check_prior = FALSE)
  ))

  expect_equal(fit_ser$alpha, fit_rss$alpha, tolerance = machine_tol)
  expect_equal(fit_ser$mu, fit_rss$mu, tolerance = machine_tol)
  expect_equal(fit_ser$mu2, fit_rss$mu2, tolerance = machine_tol)
  expect_equal(fit_ser$V, fit_rss$V, tolerance = machine_tol)
  expect_equal(coef(fit_ser)[-1], coef(fit_rss)[-1],
               tolerance = machine_tol)
})

# ---- bhat/shat input path ----

test_that("susie_ser with bhat, shat, n, and var_y matches diagonal susie_rss", {
  set.seed(4)
  p <- 25
  n <- 600
  bhat <- rnorm(p, sd = 0.04)
  shat <- runif(p, 0.04, 0.12)
  bhat[9] <- 0.65
  var_y <- 1.7

  fit_ser <- susie_ser(bhat = bhat, shat = shat, n = n, var_y = var_y,
                       estimate_prior_method = "optim", coverage = NULL)
  fit_rss <- suppressWarnings(suppressMessages(
    susie_rss(bhat = bhat, shat = shat, R = diag(p), n = n,
              var_y = var_y, L = 1, max_iter = 1,
              estimate_prior_method = "optim", coverage = NULL,
              check_prior = FALSE)
  ))

  expect_equal(fit_ser$alpha, fit_rss$alpha, tolerance = machine_tol)
  expect_equal(fit_ser$mu, fit_rss$mu, tolerance = machine_tol)
  expect_equal(fit_ser$mu2, fit_rss$mu2, tolerance = machine_tol)
  expect_equal(fit_ser$V, fit_rss$V, tolerance = machine_tol)
  expect_equal(fit_ser$X_column_scale_factors,
               fit_rss$X_column_scale_factors, tolerance = machine_tol)
  expect_equal(coef(fit_ser)[-1], coef(fit_rss)[-1],
               tolerance = machine_tol)
})

test_that("susie_ser bhat/scalar-shat path agrees with equivalent z input", {
  bhat <- c(0.2, 0.8, -0.1)
  shat <- 0.2
  fit <- susie_ser(bhat = bhat, shat = shat,
                   estimate_prior_method = "simple", coverage = NULL)
  fit_z <- susie_ser(z = bhat / shat,
                     estimate_prior_method = "simple", coverage = NULL)
  expect_equal(fit$alpha, fit_z$alpha)
  expect_equal(fit$mu, fit_z$mu)
  expect_equal(fit$mu2, fit_z$mu2)
})

# ---- coverage and CS construction ----

test_that("susie_ser default coverage returns attainable CS and hint", {
  z <- c(5, 3, 0.2, -0.1)
  expect_message(
    fit <- susie_ser(z = z, estimate_prior_method = "simple",
                     coverage = 0.95),
    "Maller et al. 2012"
  )
  expect_true(all(c("cs", "coverage", "requested_coverage") %in%
                    names(fit$sets)))
  expect_equal(fit$sets$requested_coverage, 0.95)
})

test_that("susie_ser coverage NULL skips CS and hint", {
  z <- c(5, 3, 0.2, -0.1)
  expect_silent(
    fit <- susie_ser(z = z, estimate_prior_method = "simple",
                     coverage = NULL)
  )
  expect_null(fit$sets)
})

# ---- output structure ----

test_that("susie_ser output does not carry matrix-path fields", {
  set.seed(5)
  fit <- susie_ser(z = rnorm(5000), estimate_prior_method = "simple",
                   coverage = NULL)
  expect_false(any(c("R", "X", "XtX", "Rz", "XtXr") %in% names(fit)))
  expect_equal(ncol(fit$alpha), 5000)
})

test_that("susie_ser preserves variable names without inventing blank names", {
  z <- setNames(c(3, -1, 0.2), paste0("snp", 1:3))
  fit <- susie_ser(z = z, null_weight = 0.1,
                   estimate_prior_method = "simple", coverage = NULL)
  expect_equal(colnames(fit$alpha), c(names(z), "null"))
  expect_equal(names(fit$pip), names(z))

  fit_unnamed <- susie_ser(z = unname(z), null_weight = 0.1,
                           estimate_prior_method = "simple", coverage = NULL)
  expect_null(names(fit_unnamed$pip))
})

# ---- input validation ----

test_that("susie_ser validates summary-statistic inputs", {
  expect_error(susie_ser(z = 1:3, bhat = 1:3, shat = 1:3,
                         coverage = NULL),
               "either z or \\(bhat, shat\\)")
  expect_error(susie_ser(bhat = 1:3, coverage = NULL),
               "either z or \\(bhat, shat\\)")
  expect_error(susie_ser(bhat = 1:3, shat = c(1, 2),
                         coverage = NULL),
               "lengths of bhat and shat")
  expect_error(susie_ser(bhat = 1:3, shat = c(1, 0, 1),
                         coverage = NULL),
               "zero or negative")
})

test_that("susie_ser errors on invalid check_null_threshold", {
  z <- rnorm(5)
  expect_error(
    susie_ser(z = z, check_null_threshold = -1, coverage = NULL),
    "check_null_threshold must be a single nonneg"
  )
  expect_error(
    susie_ser(z = z, check_null_threshold = Inf, coverage = NULL),
    "check_null_threshold must be a single nonneg"
  )
  expect_error(
    susie_ser(z = z, check_null_threshold = c(0, 1), coverage = NULL),
    "check_null_threshold must be a single nonneg"
  )
})

test_that("susie_ser errors on invalid prior_tol", {
  z <- rnorm(5)
  expect_error(
    susie_ser(z = z, prior_tol = -1e-9, coverage = NULL),
    "prior_tol must be a single nonneg"
  )
  expect_error(
    susie_ser(z = z, prior_tol = NA_real_, coverage = NULL),
    "prior_tol must be a single nonneg"
  )
})

test_that("susie_ser errors on invalid coverage", {
  z <- rnorm(5)
  expect_error(
    susie_ser(z = z, coverage = 0),
    "coverage must be NULL or a single number between 0 and 1"
  )
  expect_error(
    susie_ser(z = z, coverage = 1),
    "coverage must be NULL or a single number between 0 and 1"
  )
  expect_error(
    susie_ser(z = z, coverage = c(0.9, 0.95)),
    "coverage must be NULL or a single number between 0 and 1"
  )
})

test_that("susie_ser errors on invalid n", {
  z <- rnorm(5)
  for (bad_n in list(1, -10, Inf, c(100, 200))) {
    expect_error(
      susie_ser(z = z, n = bad_n, coverage = NULL),
      "n must be a single number greater than 1"
    )
  }
})

test_that("susie_ser errors on invalid var_y", {
  z <- rnorm(5)
  for (bad_var_y in list(0, -1, Inf)) {
    expect_error(
      susie_ser(z = z, n = 100, var_y = bad_var_y, coverage = NULL),
      "var_y must be a single positive finite"
    )
  }
})

test_that("susie_ser errors on invalid prior_variance", {
  z <- rnorm(5)
  expect_error(
    susie_ser(z = z, prior_variance = -1, coverage = NULL),
    "prior_variance must be a single positive finite"
  )
  expect_error(
    susie_ser(z = z, prior_variance = 0, coverage = NULL),
    "prior_variance must be a single positive finite"
  )
})

test_that("susie_ser errors on invalid scaled_prior_variance", {
  z <- rnorm(5)
  expect_error(
    susie_ser(z = z, n = 100, scaled_prior_variance = 0, coverage = NULL),
    "scaled_prior_variance must be a single positive finite"
  )
})

# ---- infinite-value input errors ----

test_that("susie_ser stops on non-finite bhat/shat", {
  bhat <- c(1, Inf, 0.5)
  shat <- c(0.1, 0.1, 0.1)
  expect_error(
    susie_ser(bhat = bhat, shat = shat, coverage = NULL),
    "bhat and shat must be finite"
  )
})

test_that("susie_ser stops on infinite z", {
  expect_error(
    susie_ser(z = c(1.0, Inf, -0.5), coverage = NULL),
    "z contains infinite values"
  )
})

# ---- name_ser_output ----

test_that("name_ser_output with null_weight attaches names and excludes null from pip", {
  set.seed(201)
  z <- setNames(rnorm(4), paste0("var", 1:4))

  fit <- susie_ser(z = z, null_weight = 0.1,
                   estimate_prior_method = "simple", coverage = NULL)

  expect_equal(ncol(fit$alpha), 5L)
  expect_equal(colnames(fit$alpha)[5], "null")
  expect_length(fit$pip, 4L)
  expect_equal(names(fit$pip), paste0("var", 1:4))
  expect_equal(names(fit$z), paste0("var", 1:4))
})

test_that("name_ser_output with unnamed z leaves pip names as NULL", {
  set.seed(202)
  fit <- susie_ser(z = rnorm(5), null_weight = 0.1,
                   estimate_prior_method = "simple", coverage = NULL)

  expect_null(names(fit$pip))
  expect_null(colnames(fit$alpha))
})

test_that("name_ser_output returns model unchanged when variable_names length mismatches", {
  set.seed(203)
  z <- setNames(rnorm(5), paste0("snp", 1:5))

  fit0 <- susie_ser(z = z, estimate_prior_method = "simple", coverage = NULL)
  result <- name_ser_output(fit0, variable_names = rep("x", 4))
  expect_equal(colnames(result$alpha), colnames(fit0$alpha))
})
