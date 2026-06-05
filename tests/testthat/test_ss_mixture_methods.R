context("ss_mixture methods")

# ---- Local helper ----

# Build a small 2-panel ss_mixture data object (p variables, 2 related panels).
make_ss_mixture <- function(p = 12, n = 500, B1 = 60, B2 = 60, seed = 1,
                            L = 3, R_mismatch = "none", R_finite = NULL) {
  set.seed(seed)
  make_R <- function(s, B) {
    set.seed(s)
    X <- matrix(rnorm(B * p), B, p)
    cov2cor(crossprod(scale(X)))
  }
  R1 <- make_R(seed, B1)
  R2 <- make_R(seed + 1, B2)
  z <- rnorm(p)
  z[3] <- 4
  ctor <- ss_mixture_constructor(
    z = z, R = list(R1, R2), n = n, L = L,
    R_mismatch = R_mismatch, R_finite = R_finite
  )
  var_y <- get_var_y.ss(ctor$data)
  model <- initialize_susie_model.ss(ctor$data, ctor$params, var_y)
  model$alpha <- matrix(1 / p, L, p)
  model$alpha[1, 3] <- 0.8
  model$alpha[1, -3] <- 0.2 / (p - 1)
  model$mu <- matrix(0.02, L, p)
  model$mu2 <- model$mu^2 + 0.01
  model$sigma2 <- 1
  model$XtXr <- compute_XtXv_mixture(ctor$data, model,
                                     colSums(model$alpha * model$mu))
  list(data = ctor$data, params = ctor$params, model = model, z = z)
}

# ---- compute_XtXv_mixture ----

test_that("compute_XtXv_mixture uses panel_R branch and matches (n-1)*R(omega)*v", {
  setup <- make_ss_mixture(p = 10, seed = 11)
  data <- setup$data
  model <- setup$model

  expect_false(is.null(data$panel_R))
  expect_true(is.null(data$omega_cache))

  v <- rnorm(data$p)
  out <- compute_XtXv_mixture(data, model, v)

  expect_length(out, data$p)
  expect_true(all(is.finite(out)))

  omega <- get_mixture_omega(data, model)
  Rv_ref <- Reduce("+", Map(function(w, R) w * (R %*% v), omega, data$panel_R))
  expect_equal(out, data$nm1 * as.vector(Rv_ref), tolerance = 1e-10)
})

test_that("compute_XtXv_mixture falls back to compute_Rv when panel_R is absent", {
  setup <- make_ss_mixture(p = 8, seed = 12)
  data <- setup$data
  model <- setup$model

  set.seed(99)
  Xmeta <- matrix(rnorm(40 * data$p), 40, data$p)
  attr(Xmeta, "d") <- rep(data$nm1, data$p)
  data_x <- data
  data_x$panel_R <- NULL
  data_x$X <- Xmeta
  data_x$XtX <- NULL

  v <- rnorm(data$p)
  out <- compute_XtXv_mixture(data_x, model, v)
  expect_length(out, data$p)
  expect_equal(out, as.vector(compute_Rv(data_x, v)), tolerance = 1e-10)
})

# ---- get_ER2.ss_mixture ----

test_that("get_ER2.ss_mixture returns a finite nonnegative scalar", {
  setup <- make_ss_mixture(p = 12, seed = 13)
  er2 <- get_ER2.ss_mixture(setup$data, setup$model)
  expect_type(er2, "double")
  expect_length(er2, 1)
  expect_true(is.finite(er2))
  expect_true(er2 >= 0)
})

test_that("get_ER2.ss_mixture matches its explicit definition", {
  setup <- make_ss_mixture(p = 10, seed = 14)
  data <- setup$data
  model <- setup$model

  B       <- model$alpha * model$mu
  betabar <- colSums(B)
  postb2  <- model$alpha * model$mu2
  XtX_betabar <- compute_XtXv_mixture(data, model, betabar)
  XB2 <- 0
  for (l in seq_len(nrow(B))) {
    bl  <- B[l, ]
    XB2 <- XB2 + sum(bl * compute_XtXv_mixture(data, model, bl))
  }
  expected <- data$yty - 2 * sum(betabar * data$Xty) +
    sum(betabar * XtX_betabar) - XB2 + data$nm1 * sum(postb2)
  expect_equal(get_ER2.ss_mixture(data, model), expected, tolerance = 1e-10)
})

# ---- update_model_variance.ss_mixture (omega M-step) ----

test_that("update_model_variance.ss_mixture runs omega M-step and returns simplex omega", {
  setup <- make_ss_mixture(p = 12, seed = 15, L = 3)
  data <- setup$data
  params <- setup$params
  model <- setup$model

  model$zbar        <- colSums(model$alpha * model$mu)
  model$diag_postb2 <- colSums(model$alpha * model$mu2)
  model$Z           <- model$alpha * model$mu

  expect_equal(data$K, 2)
  expect_false(isTRUE(model$omega_converged))

  res <- update_model_variance.ss_mixture(data, params, model)

  expect_length(res$omega, data$K)
  expect_true(all(res$omega >= 0))
  expect_equal(sum(res$omega), 1, tolerance = 1e-8)
  expect_length(res$XtXr, data$p)
  expect_true(all(is.finite(res$XtXr)))
})

test_that("update_model_variance.ss_mixture updates sigma2 when estimate_residual_variance=TRUE", {
  setup <- make_ss_mixture(p = 12, seed = 16, L = 2)
  data <- setup$data
  params <- setup$params
  model <- setup$model
  params$estimate_residual_variance <- TRUE
  model$zbar        <- colSums(model$alpha * model$mu)
  model$diag_postb2 <- colSums(model$alpha * model$mu2)
  model$Z           <- model$alpha * model$mu

  res <- update_model_variance.ss_mixture(data, params, model)
  expect_true(is.finite(res$sigma2))
  expect_true(res$sigma2 > 0)
})

test_that("update_model_variance.ss_mixture skips omega M-step once omega_converged=TRUE", {
  setup <- make_ss_mixture(p = 10, seed = 17, L = 2)
  data <- setup$data
  params <- setup$params
  model <- setup$model
  params$estimate_residual_variance <- FALSE
  model$omega_converged <- TRUE
  model$omega <- get_mixture_omega(data, model)
  omega_before <- model$omega

  res <- update_model_variance.ss_mixture(data, params, model)
  expect_equal(get_mixture_omega(data, res), omega_before, tolerance = 1e-12)
})

# ---- end-to-end multi-panel tests ----

test_that("susie_rss with list of R panels dispatches through ss_mixture and returns valid fit", {
  set.seed(18)
  p <- 12; Bn <- 60
  make_R <- function(s) {
    set.seed(s)
    cov2cor(crossprod(scale(matrix(rnorm(Bn * p), Bn, p))))
  }
  z <- rnorm(p); z[4] <- 4
  fit <- suppressWarnings(susie_rss(z = z, R = list(make_R(1), make_R(2)),
                                    n = 1000, L = 3, max_iter = 5, verbose = FALSE))
  expect_s3_class(fit, "susie")
  expect_length(fit$pip, p)
  expect_false(is.null(fit$omega_weights))
  expect_equal(sum(fit$omega_weights), 1, tolerance = 1e-8)
})

test_that("susie_rss with X list builds omega_cache path and returns valid fit", {
  set.seed(19)
  p <- 30; n1 <- 10; n2 <- 8
  z <- rnorm(p); z[5] <- 3.5
  fit <- suppressWarnings(susie_rss(
    z = z, X = list(matrix(rnorm(n1 * p), n1, p), matrix(rnorm(n2 * p), n2, p)),
    n = 500, L = 2, max_iter = 5, verbose = FALSE
  ))
  expect_s3_class(fit, "susie")
  expect_length(fit$pip, p)
  expect_false(is.null(fit$omega_weights))
  expect_equal(sum(fit$omega_weights), 1, tolerance = 1e-8)
})
