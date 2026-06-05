context("slot_prior class")

# ---- Constructor defaults and fields ----

test_that("slot_prior_betabinom constructs with default a/b parameters", {
  sp <- slot_prior_betabinom()
  expect_s3_class(sp, "slot_prior_betabinom")
  expect_s3_class(sp, "slot_prior")
  expect_equal(sp$a_beta, 1)
  expect_equal(sp$b_beta, 2)
  expect_true(sp$ab_was_default)
  expect_equal(sp$update_schedule, "sequential")
  expect_null(sp$c_hat_init)
  expect_equal(sp$skip_threshold_multiplier, 0)
})

test_that("slot_prior_betabinom accepts explicit a/b and marks them non-default", {
  sp <- slot_prior_betabinom(a_beta = 2, b_beta = 3)
  expect_equal(sp$a_beta, 2)
  expect_equal(sp$b_beta, 3)
  expect_false(sp$ab_was_default)
})

test_that("slot_prior_poisson constructs with explicit nu", {
  sp <- suppressMessages(slot_prior_poisson(C = 5, nu = 8))
  expect_s3_class(sp, "slot_prior_poisson")
  expect_s3_class(sp, "slot_prior")
  expect_equal(sp$C, 5)
  expect_equal(sp$nu, 8)
  expect_false(sp$nu_was_default)
  expect_equal(sp$update_schedule, "sequential")
  expect_null(sp$c_hat_init)
  expect_equal(sp$skip_threshold_multiplier, 0)
})

test_that("slot_prior_poisson uses default nu=8 when omitted and records it", {
  sp <- slot_prior_poisson(C = 4)
  expect_true(sp$nu_was_default)
  expect_equal(sp$nu, 8)
})

# ---- Input validation ----

test_that("slot_prior_betabinom rejects invalid a_beta/b_beta", {
  expect_error(slot_prior_betabinom(a_beta = -1))
  expect_error(slot_prior_betabinom(b_beta = -1))
  expect_error(slot_prior_betabinom(a_beta = "x"))
})

test_that("slot_prior_poisson rejects invalid C and nu", {
  expect_error(slot_prior_poisson(C = -1, nu = 8))
  expect_error(slot_prior_poisson(C = "abc", nu = 8))
  expect_error(slot_prior_poisson(C = 5, nu = -1))
})

# ---- is.slot_prior ----

test_that("is.slot_prior correctly identifies slot_prior objects", {
  expect_true(is.slot_prior(slot_prior_poisson(C = 5, nu = 8)))
  expect_true(is.slot_prior(slot_prior_betabinom()))
  expect_false(is.slot_prior(list(C = 5)))
  expect_false(is.slot_prior(NULL))
})

# ---- print.slot_prior ----

test_that("print.slot_prior shows beta-binomial header and a/b values", {
  sp <- slot_prior_betabinom(a_beta = 1, b_beta = 2)
  out <- capture.output(print(sp))
  expect_true(any(grepl("beta-binomial", out, fixed = TRUE)))
  expect_true(any(grepl("a_beta", out, fixed = TRUE)))
  expect_true(any(grepl("b_beta", out, fixed = TRUE)))
})

test_that("print.slot_prior for betabinom returns the object invisibly", {
  sp <- slot_prior_betabinom()
  ret <- withVisible(print(sp))
  expect_false(ret$visible)
  expect_identical(ret$value, sp)
})

test_that("print.slot_prior for poisson shows C, nu and update_schedule", {
  sp <- slot_prior_poisson(C = 3, nu = 8)
  out <- capture.output(print(sp))
  expect_true(any(grepl("poisson", out, fixed = TRUE)))
  expect_true(any(grepl("C \\(expected", out)))
  expect_true(any(grepl("nu \\(overdispersion", out)))
  expect_true(any(grepl("update schedule", out, fixed = TRUE)))
})

test_that("print.slot_prior shows warm-start line only when c_hat_init is set", {
  sp_warm <- slot_prior_poisson(C = 3, nu = 8, c_hat_init = rep(0.3, 5))
  out_warm <- capture.output(print(sp_warm))
  expect_true(any(grepl("warm start", out_warm, fixed = TRUE)))
  expect_true(any(grepl("5-vector", out_warm, fixed = TRUE)))

  sp_cold <- slot_prior_betabinom()
  out_cold <- capture.output(print(sp_cold))
  expect_false(any(grepl("warm start", out_cold, fixed = TRUE)))
})

# ---- slot_prior_elbo ELBO contribution ----

test_that("slot_prior_elbo returns finite scalar for betabinom state", {
  L <- 5
  model <- list(
    slot_weights = rep(0.3, L),
    c_hat_state  = list(
      prior_type = "betabinom",
      a_beta = 1, b_beta = 2
    )
  )
  elbo <- slot_prior_elbo(model)
  expect_length(elbo, 1)
  expect_true(is.finite(elbo))
})

test_that("slot_prior_elbo returns finite scalar for poisson state", {
  L <- 10
  model <- list(
    slot_weights = rep(0.3, L),
    c_hat_state  = list(
      prior_type = "poisson",
      C = 3, nu = 8,
      a_g = 8 + 3, b_g = 8 / 3 + 1
    )
  )
  elbo <- slot_prior_elbo(model)
  expect_length(elbo, 1)
  expect_true(is.finite(elbo))
})

test_that("slot_prior_elbo betabinom ELBO penalizes dense configurations under sparse prior", {
  # Beta(1,2) prior (default) favors sparse configs: E[rho]=1/3
  # Sparse chat should have higher ELBO than dense chat
  make_model <- function(chat) {
    list(
      slot_weights = chat,
      c_hat_state  = list(prior_type = "betabinom", a_beta = 1, b_beta = 2)
    )
  }
  elbo_sparse <- slot_prior_elbo(make_model(rep(0.1, 10)))
  elbo_active <- slot_prior_elbo(make_model(rep(0.9, 10)))
  expect_true(is.finite(elbo_sparse))
  expect_true(is.finite(elbo_active))
  expect_gt(elbo_sparse, elbo_active)
})

# ---- susie integration: slot_prior produces c_hat ----

test_that("susie with slot_prior produces valid c_hat output", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[1:3] <- 1
  y <- X %*% b + rnorm(n)
  fit <- susie(X, y, L = 10, slot_prior = slot_prior_poisson(C = 3, nu = 8),
               verbose = FALSE)
  expect_false(is.null(fit$c_hat))
  expect_equal(length(fit$c_hat), 10)
  expect_true(all(fit$c_hat >= 0 & fit$c_hat <= 1))
  expect_false(is.null(fit$C_hat))
  expect_true(fit$C_hat > 0)
})

test_that("susie without slot_prior does not produce c_hat", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  fit <- susie(X, y, L = 5, verbose = FALSE)
  expect_null(fit$c_hat)
})

# ---- betabinom integration ----

test_that("susie with betabinom slot_prior produces valid c_hat output", {
  set.seed(2)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[1:3] <- 1
  y <- X %*% b + rnorm(n)
  fit <- susie(X, y, L = 10, slot_prior = slot_prior_betabinom(),
               verbose = FALSE)
  expect_false(is.null(fit$c_hat))
  expect_equal(length(fit$c_hat), 10)
  expect_true(all(fit$c_hat >= 0 & fit$c_hat <= 1))
})

# ---- ash model auto-creates slot_prior ----

test_that("ash model auto-creates betabinom slot_prior with message", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  expect_message(
    fit <- susie(X, y, L = 10, unmappable_effects = "ash",
                 verbose = FALSE, max_iter = 5),
    "slot_prior was not specified"
  )
  expect_false(is.null(fit$c_hat))
})

test_that("ash model with explicit slot_prior does not warn about C", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  fit <- withCallingHandlers(
    susie(X, y, L = 10, unmappable_effects = "ash",
          slot_prior = slot_prior_poisson(C = 3, nu = 8),
          verbose = FALSE, max_iter = 5),
    warning = function(w) {
      if (grepl("strongly advised", conditionMessage(w)))
        stop("Got unexpected C warning")
      invokeRestart("muffleWarning")
    }
  )
  expect_false(is.null(fit$c_hat))
})

# ---- Warm start ----

test_that("c_hat warm start does not increase iteration count", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[1:3] <- 1
  y <- X %*% b + rnorm(n)

  fit1 <- susie(X, y, L = 10, slot_prior = slot_prior_poisson(C = 3, nu = 8),
                verbose = FALSE)

  sp_warm <- slot_prior_poisson(C = 3, nu = 8, c_hat_init = fit1$c_hat)
  fit2 <- susie(X, y, L = 10, slot_prior = sp_warm,
                model_init = fit1, verbose = FALSE)

  expect_true(fit2$niter <= fit1$niter)
})

# ---- Batch vs sequential update schedules ----

test_that("batch and sequential update schedules both converge", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[1:3] <- 1
  y <- X %*% b + rnorm(n)
  fit_batch <- susie(X, y, L = 10,
                     slot_prior = slot_prior_poisson(C = 3,
                                                     update_schedule = "batch"),
                     verbose = FALSE)
  fit_seq <- susie(X, y, L = 10,
                   slot_prior = slot_prior_poisson(C = 3,
                                                   update_schedule = "sequential"),
                   verbose = FALSE)
  expect_true(fit_batch$converged)
  expect_true(fit_seq$converged)
  expect_equal(sum(fit_batch$c_hat > 0.5), sum(fit_seq$c_hat > 0.5))
})
