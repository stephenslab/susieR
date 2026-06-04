context("Iterative Bayesian Stepwise Selection (IBSS)")

# ---- ibss_initialize: structure ----

test_that("ibss_initialize returns susie-class list with correct dimensions", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)

  expect_s3_class(model, "susie")
  expect_type(model, "list")
  expect_equal(dim(model$alpha), c(5L, 50L))
  expect_equal(dim(model$mu),    c(5L, 50L))
  expect_equal(dim(model$mu2),   c(5L, 50L))
  expect_equal(dim(model$lbf_variable), c(5L, 50L))
  expect_length(model$V,   5L)
  expect_length(model$lbf, 5L)
  expect_length(model$KL,  5L)
  expect_length(model$Xr,  100L)
  for (nm in c("alpha", "mu", "mu2", "V", "sigma2", "lbf", "lbf_variable",
               "KL", "pi", "predictor_weights", "Xr", "null_index"))
    expect_true(nm %in% names(model))
})

test_that("ibss_initialize reduces L to p when L > p", {
  setup <- setup_individual_data(n = 100, p = 10, L = 20)
  model <- ibss_initialize(setup$data, setup$params)
  expect_equal(nrow(model$alpha), 10L)
  expect_length(model$V, 10L)
})

# ---- ibss_initialize: residual variance validation ----

test_that("ibss_initialize errors on negative residual variance", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$residual_variance <- -1
  expect_error(ibss_initialize(setup$data, setup$params),
               "Residual variance sigma2 must be positive")
})

test_that("ibss_initialize errors on non-scalar residual variance", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$residual_variance <- c(1, 2)
  expect_error(ibss_initialize(setup$data, setup$params),
               "Input residual variance sigma2 must be a scalar")
})

test_that("ibss_initialize errors on non-numeric residual variance", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$residual_variance <- "one"
  expect_error(ibss_initialize(setup$data, setup$params),
               "Input residual variance sigma2 must be numeric")
})

test_that("ibss_initialize defaults sigma2 to var(y) when residual_variance is NULL", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$residual_variance <- NULL
  var_y <- var(drop(setup$data$y))
  model <- ibss_initialize(setup$data, setup$params)
  expect_equal(model$sigma2, var_y, tolerance = 1e-8)
})

test_that("ibss_initialize uses provided residual variance", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$residual_variance <- 2.5
  model <- ibss_initialize(setup$data, setup$params)
  expect_equal(model$sigma2, 2.5)
})

# ---- ibss_initialize: model_init warm-start ----

test_that("ibss_initialize works without model_init", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$model_init <- NULL
  model <- ibss_initialize(setup$data, setup$params)
  expect_equal(dim(model$alpha), c(5L, 50L))
  expect_true(all(model$alpha >= 0 & model$alpha <= 1))
  expect_true(all(is.finite(model$mu)))
})

test_that("ibss_initialize expands L when model_init has fewer effects", {
  setup <- setup_individual_data(n = 100, p = 50, L = 2)
  model_init <- ibss_initialize(setup$data, setup$params)
  model_init$V <- rep(0.5, 2)

  setup2 <- setup_individual_data(n = 100, p = 50, L = 5)
  setup2$params$model_init <- model_init
  model <- ibss_initialize(setup2$data, setup2$params)
  expect_equal(dim(model$alpha), c(5L, 50L))
})

test_that("ibss_initialize keeps all effects when model_init has more than L", {
  setup <- setup_individual_data(n = 100, p = 50, L = 6)
  model_init <- ibss_initialize(setup$data, setup$params)
  model_init$V <- rep(0.5, 6)

  setup2 <- setup_individual_data(n = 100, p = 50, L = 3)
  setup2$params$model_init <- model_init
  expect_message(
    model <- ibss_initialize(setup2$data, setup2$params),
    "using L = 6"
  )
  expect_equal(dim(model$alpha), c(6L, 50L))
})

# ---- ibss_initialize: mathematical properties ----

test_that("ibss_initialize produces valid initial state", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)

  # alpha rows sum to 1 and are valid probabilities
  expect_equal(rowSums(model$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(model$alpha >= 0 & model$alpha <= 1))

  # V non-negative and finite
  expect_true(all(model$V >= 0))
  expect_true(all(is.finite(model$V)))

  # sigma2 positive and finite
  expect_gt(model$sigma2, 0)
  expect_true(is.finite(model$sigma2))

  # KL and lbf initialized to NA before first fit
  setup$params$model_init <- NULL
  m2 <- ibss_initialize(setup$data, setup$params)
  expect_true(all(is.na(m2$KL)))
  expect_true(all(is.na(m2$lbf)))
})

# ---- ibss_initialize: data type compatibility ----

for (.dtype in c("individual", "ss", "rss_lambda")) {
  local({
    dtype <- .dtype
    test_that(paste("ibss_initialize creates correct fitted field for", dtype), {
      if (dtype == "individual") {
        setup <- setup_individual_data(n = 100, p = 50, L = 5)
        model <- ibss_initialize(setup$data, setup$params)
        expect_s3_class(model, "susie")
        expect_length(model$Xr, 100)
        expect_true(all(is.finite(model$Xr)))
      } else if (dtype == "ss") {
        setup <- setup_ss_data(n = 100, p = 50, L = 5)
        model <- ibss_initialize(setup$data, setup$params)
        expect_s3_class(model, "susie")
        expect_length(model$XtXr, 50)
        expect_true(all(is.finite(model$XtXr)))
      } else {
        setup <- setup_rss_lambda_data(n = 500, p = 50, L = 5, lambda = 0.5)
        model <- ibss_initialize(setup$data, setup$params)
        expect_s3_class(model, "susie")
        expect_length(model$Rz, 50)
        expect_true(all(is.finite(model$Rz)))
      }
    })
  })
}

# ---- ibss_initialize: null_index ----

test_that("ibss_initialize sets null_index to 0 when null_weight = 0", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$model$null_weight <- 0
  model <- ibss_initialize(setup$data, setup$params)
  expect_equal(model$null_index, 0)
})

test_that("ibss_initialize sets positive null_index when null_weight > 0", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$null_weight <- 0.5
  model <- ibss_initialize(setup$data, setup$params)
  expect_gt(model$null_index, 0)
})

# ---- ibss_fit: updates and mathematical properties ----

test_that("ibss_fit maintains valid distributions and finite values after one step", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  upd   <- ibss_fit(setup$data, setup$params, model)

  # alpha
  expect_equal(rowSums(upd$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(upd$alpha >= 0 & upd$alpha <= 1))

  # V non-negative finite
  expect_length(upd$V, 5)
  expect_true(all(upd$V >= 0))
  expect_true(all(is.finite(upd$V)))

  # lbf finite
  expect_length(upd$lbf, 5)
  expect_true(all(is.finite(upd$lbf)))

  # KL non-negative
  expect_length(upd$KL, 5)
  expect_true(all(is.finite(upd$KL)))
  expect_true(all(upd$KL >= -1e-6))

  # mu and mu2 finite
  expect_true(all(is.finite(upd$mu)))
  expect_true(all(is.finite(upd$mu2)))

  # Xr updated
  expect_length(upd$Xr, 100)
  expect_true(all(is.finite(upd$Xr)))
})

test_that("ibss_fit works with L=1", {
  setup <- setup_individual_data(n = 100, p = 50, L = 1)
  model <- ibss_initialize(setup$data, setup$params)
  upd   <- ibss_fit(setup$data, setup$params, model)

  expect_equal(dim(upd$alpha), c(1L, 50L))
  expect_equal(sum(upd$alpha), 1, tolerance = 1e-10)
})

test_that("ibss_fit works with L=0", {
  setup <- setup_individual_data(n = 100, p = 50, L = 0)
  model <- list(alpha = matrix(0, 0, 50))
  upd   <- ibss_fit(setup$data, setup$params, model)
  expect_equal(nrow(upd$alpha), 0)
})

test_that("ibss_fit is idempotent over multiple iterations", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  for (iter in 1:3) {
    model <- ibss_fit(setup$data, setup$params, model)
    expect_equal(rowSums(model$alpha), rep(1, 5), tolerance = 1e-10)
    expect_true(all(model$V >= 0))
  }
})

test_that("ibss_fit produces valid output for each data type", {
  for (dtype in c("individual", "ss", "rss_lambda")) {
    setup <- switch(dtype,
      individual  = setup_individual_data(n = 100, p = 50, L = 5),
      ss          = setup_ss_data(n = 100, p = 50, L = 5),
      rss_lambda  = setup_rss_lambda_data(n = 500, p = 50, L = 5, lambda = 0.5)
    )
    model <- ibss_initialize(setup$data, setup$params)
    upd   <- ibss_fit(setup$data, setup$params, model)
    expect_s3_class(upd, "susie")
  }
})

# ---- ibss_finalize: output fields and values ----

test_that("ibss_finalize produces complete finalized susie object", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)
  final <- ibss_finalize(setup$data, setup$params, model,
                         elbo = NULL, iter = 42L, tracking = NULL)

  expect_s3_class(final, "susie")
  expect_equal(final$niter, 42L)

  # pip valid probabilities
  expect_length(final$pip, 50)
  expect_true(all(final$pip >= 0 & final$pip <= 1))
  expect_true(all(is.finite(final$pip)))

  # sets is a list
  expect_type(final$sets, "list")

  # fitted values length n and finite
  expect_length(final$fitted, 100)
  expect_true(all(is.finite(final$fitted)))

  # intercept finite
  expect_true(is.finite(final$intercept))

  # scale factors
  expect_length(final$X_column_scale_factors, 50)

  # z-scores present and finite
  expect_true("z" %in% names(final))
  expect_length(final$z, 50)
  expect_true(all(is.finite(final$z)))
})

test_that("ibss_finalize attaches variable names to pip", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  colnames(setup$data$X) <- paste0("var", 1:50)
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)
  final <- ibss_finalize(setup$data, setup$params, model,
                         elbo = NULL, iter = 10L, tracking = NULL)
  expect_named(final$pip, paste0("var", 1:50))
})

test_that("ibss_finalize includes trace when track_fit=TRUE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$track_fit <- TRUE
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)
  tracking <- list(make_track_snapshot(model, 0L))
  final <- ibss_finalize(setup$data, setup$params, model,
                         elbo = NULL, iter = 3L, tracking = tracking)
  expect_true("trace" %in% names(final))
  expect_s3_class(final$trace, "susie_track")
})

test_that("ibss_finalize omits trace when track_fit=FALSE", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  setup$params$track_fit <- FALSE
  model <- ibss_initialize(setup$data, setup$params)
  model <- ibss_fit(setup$data, setup$params, model)
  final <- ibss_finalize(setup$data, setup$params, model,
                         elbo = NULL, iter = 3L, tracking = NULL)
  expect_false("trace" %in% names(final))
})

# ---- full IBSS pipeline ----

test_that("full IBSS pipeline produces valid susie object with correct properties", {
  setup <- setup_individual_data(n = 100, p = 50, L = 5)
  model <- ibss_initialize(setup$data, setup$params)
  for (i in 1:5)
    model <- ibss_fit(setup$data, setup$params, model)
  model <- ibss_finalize(setup$data, setup$params, model,
                         elbo = NULL, iter = 5L, tracking = NULL)

  expect_s3_class(model, "susie")
  expect_equal(model$niter, 5L)
  expect_equal(rowSums(model$alpha), rep(1, 5), tolerance = 1e-10)
  expect_true(all(model$pip >= 0 & model$pip <= 1))
  expect_true(all(model$V >= 0))
  for (nm in c("alpha", "mu", "V", "pip", "sets", "fitted", "niter"))
    expect_true(nm %in% names(model))
})

# ---- slot activity prior (c_hat) initialization ----

test_that("ibss_initialize sets up beta-binomial c_hat_state at prior mean", {
  setup <- setup_individual_data(n = 80, p = 40, L = 4)
  setup$params$slot_prior <- slot_prior_betabinom()
  model <- ibss_initialize(setup$data, setup$params)

  expect_false(is.null(model$c_hat_state))
  expect_equal(model$c_hat_state$prior_type, "betabinom")
  expect_length(model$slot_weights, 4)
  # Beta(1, 2) prior mean = 1/3
  expect_equal(model$slot_weights, rep(1/3, 4), tolerance = 1e-8)
  expect_equal(model$c_hat_state$a_beta, 1)
  expect_equal(model$c_hat_state$b_beta, 2)
})

test_that("ibss_initialize honors beta-binomial c_hat_init warm start", {
  setup <- setup_individual_data(n = 80, p = 40, L = 3)
  setup$params$slot_prior <- slot_prior_betabinom(c_hat_init = c(0.2, 0.3, 0.4))
  model <- ibss_initialize(setup$data, setup$params)
  expect_equal(model$slot_weights, c(0.2, 0.3, 0.4), tolerance = 1e-8)
})

test_that("ibss_initialize sets up gamma-poisson c_hat_state", {
  setup <- setup_ss_data(n = 100, p = 40, L = 4)
  setup$params$slot_prior <- slot_prior_poisson(C = 3, nu = 8)
  model <- ibss_initialize(setup$data, setup$params)

  expect_false(is.null(model$c_hat_state))
  expect_equal(model$c_hat_state$prior_type, "poisson")
  expect_equal(model$c_hat_state$C, 3)
  expect_equal(model$c_hat_state$nu, 8)
  expect_equal(model$c_hat_state$a_g, 8 + 3)
})

test_that("ibss_initialize honors gamma-poisson c_hat_init warm start", {
  setup <- setup_ss_data(n = 100, p = 40, L = 4)
  setup$params$slot_prior <- slot_prior_poisson(
    C = 3, nu = 8, c_hat_init = c(0.1, 0.2, 0.3, 0.4))
  model <- ibss_initialize(setup$data, setup$params)

  expect_equal(model$slot_weights, c(0.1, 0.2, 0.3, 0.4), tolerance = 1e-8)
  expect_equal(model$c_hat_state$a_g, 8 + sum(c(0.1, 0.2, 0.3, 0.4)))
})

test_that("ibss_initialize rejects invalid slot_prior object", {
  setup <- setup_individual_data(n = 80, p = 40, L = 4)
  setup$params$slot_prior <- list(not = "a real slot prior")
  expect_error(ibss_initialize(setup$data, setup$params), "slot_prior must be created by")
})

test_that("ibss_initialize computes finite weighted fitted values for c_hat models", {
  setup_ind <- setup_individual_data(n = 80, p = 40, L = 4)
  setup_ind$params$slot_prior <- slot_prior_betabinom()
  m_ind <- ibss_initialize(setup_ind$data, setup_ind$params)
  expect_true(all(is.finite(m_ind$Xr)))

  setup_ss <- setup_ss_data(n = 100, p = 40, L = 4)
  setup_ss$params$slot_prior <- slot_prior_poisson(C = 3, nu = 8)
  m_ss <- ibss_initialize(setup_ss$data, setup_ss$params)
  expect_true(all(is.finite(m_ss$XtXr)))
})

# ---- ibss_fit: slot activity (c_hat) update loop ----

test_that("ibss_fit updates beta-binomial slot weights to valid probabilities", {
  set.seed(900)
  setup <- setup_individual_data(n = 80, p = 40, L = 4)
  setup$params$slot_prior <- slot_prior_betabinom()
  model <- ibss_initialize(setup$data, setup$params)
  model$runtime <- list(prev_elbo = -Inf, prev_alpha = model$alpha)
  upd <- ibss_fit(setup$data, setup$params, model)

  expect_length(upd$slot_weights, 4)
  expect_true(all(upd$slot_weights >= 0 & upd$slot_weights <= 1))
})

test_that("ibss_fit computes positive skip threshold after beta-binomial iteration", {
  set.seed(901)
  setup <- setup_individual_data(n = 80, p = 40, L = 4)
  setup$params$slot_prior <- slot_prior_betabinom(skip_threshold_multiplier = 0.5)
  model <- ibss_initialize(setup$data, setup$params)
  model$runtime <- list(prev_elbo = -Inf, prev_alpha = model$alpha)
  upd <- ibss_fit(setup$data, setup$params, model)

  expect_gt(upd$c_hat_state$skip_threshold, 0)
})

test_that("ibss_fit performs batch a_g update satisfying nu + sum(weights) identity", {
  set.seed(902)
  setup <- setup_individual_data(n = 80, p = 40, L = 4)
  setup$params$slot_prior <- slot_prior_poisson(
    C = 3, update_schedule = "batch", skip_threshold_multiplier = 0.5)
  model <- ibss_initialize(setup$data, setup$params)
  model$runtime <- list(prev_elbo = -Inf, prev_alpha = model$alpha)
  upd <- ibss_fit(setup$data, setup$params, model)

  expect_equal(upd$c_hat_state$a_g,
               upd$c_hat_state$nu + sum(upd$slot_weights),
               tolerance = 1e-8)
  expect_gt(upd$c_hat_state$skip_threshold, 0)
})

test_that("ibss_fit skips slots below the skip threshold, preserving their alpha row", {
  set.seed(903)
  setup <- setup_individual_data(n = 80, p = 40, L = 4)
  setup$params$slot_prior <- slot_prior_poisson(
    C = 3, update_schedule = "batch", skip_threshold_multiplier = 0.5)
  model <- ibss_initialize(setup$data, setup$params)
  model$runtime <- list(prev_elbo = -Inf, prev_alpha = model$alpha)
  model <- ibss_fit(setup$data, setup$params, model)

  # Force slot 1 far below threshold so next iteration skips it
  model$slot_weights[1] <- 1e-12
  alpha_before <- model$alpha[1, ]
  upd <- ibss_fit(setup$data, setup$params, model)

  expect_equal(upd$slot_weights[1], 1e-12)
  expect_equal(upd$alpha[1, ], alpha_before)
})

test_that("ibss_fit works with c_hat on sufficient statistics data", {
  set.seed(904)
  setup <- setup_ss_data(n = 100, p = 40, L = 4)
  setup$params$slot_prior <- slot_prior_poisson(C = 3, nu = 8)
  model <- ibss_initialize(setup$data, setup$params)
  model$runtime <- list(prev_elbo = -Inf, prev_alpha = model$alpha)
  upd <- ibss_fit(setup$data, setup$params, model)

  expect_length(upd$slot_weights, 4)
  expect_true(all(is.finite(upd$XtXr)))
})

# ---- slot activity helpers ----

test_that("detect_fitted_field identifies Rz, XtXr, and Xr fields", {
  expect_equal(detect_fitted_field(list(Rz   = 1:3, alpha = matrix(1, 1, 3))), "Rz")
  expect_equal(detect_fitted_field(list(XtXr = 1:3, alpha = matrix(1, 1, 3))), "XtXr")
  expect_equal(detect_fitted_field(list(Xr   = 1:3, alpha = matrix(1, 1, 3))), "Xr")
})

test_that("detect_fitted_field errors when no fitted field is present", {
  expect_error(detect_fitted_field(list(alpha = matrix(1, 1, 3))),
               "Cannot detect fitted-values field")
})

test_that("update_c_hat treats NA lbf as zero and returns finite slot weight", {
  set.seed(905)
  setup <- setup_individual_data(n = 80, p = 40, L = 3)
  setup$params$slot_prior <- slot_prior_betabinom()
  model <- ibss_initialize(setup$data, setup$params)
  model$lbf[1] <- NA
  upd <- update_c_hat(setup$data, model, 1)

  expect_true(is.finite(upd$slot_weights[1]))
  expect_true(upd$slot_weights[1] >= 0 & upd$slot_weights[1] <= 1)
})

test_that("update_c_hat sequential gamma-poisson satisfies a_g = nu + sum(weights)", {
  set.seed(906)
  setup <- setup_ss_data(n = 100, p = 40, L = 4)
  setup$params$slot_prior <- slot_prior_poisson(
    C = 3, nu = 8, update_schedule = "sequential")
  model <- ibss_initialize(setup$data, setup$params)
  model$lbf[1] <- 2.0
  a_g_before <- model$c_hat_state$a_g
  upd <- update_c_hat(setup$data, model, 1)

  expect_equal(upd$c_hat_state$a_g,
               upd$c_hat_state$nu + sum(upd$slot_weights),
               tolerance = 1e-8)
  expect_false(isTRUE(all.equal(upd$c_hat_state$a_g, a_g_before)))
})

test_that("adjust_fitted_for_c_hat produces finite XtXr for ss data", {
  set.seed(907)
  setup <- setup_ss_data(n = 100, p = 40, L = 4)
  setup$params$slot_prior <- slot_prior_poisson(C = 3, nu = 8)
  model <- ibss_initialize(setup$data, setup$params)
  b_bar_l <- model$alpha[1, ] * model$mu[1, ]
  upd <- adjust_fitted_for_c_hat(setup$data, model, b_bar_l, 0.1)
  expect_true(all(is.finite(upd$XtXr)))
})

test_that("recompute_fitted_weighted produces finite Rz for rss_lambda data", {
  set.seed(908)
  setup <- setup_rss_lambda_data(n = 400, p = 40, L = 4, lambda = 0.5)
  setup$params$slot_prior <- slot_prior_betabinom()
  model <- ibss_initialize(setup$data, setup$params)
  recomputed <- recompute_fitted_weighted(setup$data, model)
  expect_true("Rz" %in% names(recomputed))
  expect_true(all(is.finite(recomputed$Rz)))
})

# ---- end-to-end slot prior fits ----

test_that("susie with beta-binomial slot_prior returns valid c_hat output", {
  set.seed(909)
  n <- 100; p <- 80
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[c(1, 20, 40)] <- c(1.2, -1, 0.9)
  y <- as.vector(X %*% b + rnorm(n))

  fit <- suppressMessages(suppressWarnings(
    susie(X, y, L = 8, slot_prior = slot_prior_betabinom(), verbose = FALSE)))

  expect_length(fit$c_hat, 8)
  expect_true(all(fit$c_hat >= 0 & fit$c_hat <= 1))
  expect_equal(fit$C_hat, sum(fit$c_hat), tolerance = 1e-8)
  expect_equal(fit$a_beta, 1)
  expect_equal(fit$b_beta, 2)
})

test_that("susie with unmappable_effects='ash' returns valid c_hat and tau2", {
  set.seed(910)
  n <- 100; p <- 120
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[1:2] <- 1.5
  y <- as.vector(X %*% b + rnorm(n))

  fit <- suppressMessages(suppressWarnings(
    susie(X, y, L = 6, unmappable_effects = "ash",
          slot_prior = slot_prior_betabinom(),
          verbose = FALSE, max_iter = 8)))

  expect_length(fit$c_hat, 6)
  expect_false(is.null(fit$tau2))
})

# ---- model_init warm-start path ----

test_that("ibss_initialize warm-start path produces valid fit on second run", {
  set.seed(911)
  n <- 100; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[c(5, 25)] <- c(1.5, -1.5)
  y <- as.vector(X %*% b + rnorm(n))

  fit0 <- suppressWarnings(susie(X, y, L = 3, max_iter = 5, verbose = FALSE))

  init <- list(alpha = fit0$alpha, mu = fit0$mu, mu2 = fit0$mu2)
  class(init) <- "susie"
  fit1 <- suppressWarnings(susie(X, y, L = 3, model_init = init, max_iter = 5, verbose = FALSE))

  expect_s3_class(fit1, "susie")
  expect_equal(dim(fit1$alpha), c(3L, p))
})
