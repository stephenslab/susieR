context("slot_prior class")

test_that("slot_prior_binomial constructs correctly", {
  sp <- suppressMessages(slot_prior_binomial(C = 5, nu = 8))
  expect_s3_class(sp, "slot_prior_binomial")
  expect_s3_class(sp, "slot_prior")
  expect_equal(sp$C, 5)
  expect_equal(sp$nu, 8)  # default when NULL
  expect_equal(sp$update_schedule, "batch")  # binomial default
  expect_null(sp$c_hat_init)
  expect_equal(sp$skip_threshold_multiplier, 0)
})

test_that("slot_prior_poisson constructs correctly", {
  sp <- slot_prior_poisson(C = 3, nu = 10)
  expect_s3_class(sp, "slot_prior_poisson")
  expect_s3_class(sp, "slot_prior")
  expect_equal(sp$C, 3)
  expect_equal(sp$nu, 10)
  expect_equal(sp$update_schedule, "sequential")  # poisson default
})

test_that("slot_prior validates inputs", {
  expect_error(slot_prior_binomial(C = -1, nu = 8))
  expect_error(slot_prior_binomial(C = "abc", nu = 8))
  expect_error(slot_prior_poisson(C = 5, nu = -1))
})

test_that("slot_prior tracks nu_was_default", {
  sp_default <- slot_prior_poisson(C = 4)
  expect_true(sp_default$nu_was_default)
  expect_equal(sp_default$nu, 8)
  sp_explicit <- slot_prior_poisson(C = 4, nu = 8)
  expect_false(sp_explicit$nu_was_default)
  expect_equal(sp_explicit$nu, 8)
})

test_that("slot_prior_binomial default for update_schedule is batch", {
  sp <- slot_prior_binomial(C = 5, nu = 8)
  expect_equal(sp$update_schedule, "batch")
})

test_that("slot_prior_poisson default for update_schedule is sequential", {
  sp <- slot_prior_poisson(C = 5, nu = 8)
  expect_equal(sp$update_schedule, "sequential")
})

test_that("is.slot_prior works", {
  expect_true(is.slot_prior(slot_prior_binomial(C = 5, nu = 8)))
  expect_true(is.slot_prior(slot_prior_poisson(C = 5, nu = 8)))
  expect_false(is.slot_prior(list(C = 5)))
  expect_false(is.slot_prior(NULL))
})

test_that("print.slot_prior produces output", {
  expect_output(print(slot_prior_binomial(C = 5, nu = 8)), "binomial")
  expect_output(print(slot_prior_poisson(C = 3, nu = 8)), "poisson")
})

test_that("susie with slot_prior produces c_hat output", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[1:3] <- 1
  y <- X %*% b + rnorm(n)
  fit <- susie(X, y, L = 10, slot_prior = slot_prior_binomial(C = 3, nu = 8),
               verbose = FALSE)
  expect_true(!is.null(fit$c_hat))
  expect_equal(length(fit$c_hat), 10)
  expect_true(all(fit$c_hat >= 0 & fit$c_hat <= 1))
  expect_true(!is.null(fit$C_hat))
  expect_true(fit$C_hat > 0)
})

test_that("susie with binomial and poisson give similar results", {
  set.seed(42)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[1:3] <- 1
  y <- X %*% b + rnorm(n)
  fit_b <- susie(X, y, L = 10, slot_prior = slot_prior_binomial(C = 3, nu = 8),
                 verbose = FALSE)
  fit_p <- susie(X, y, L = 10, slot_prior = slot_prior_poisson(C = 3, nu = 8),
                 verbose = FALSE)
  # Both should find approximately the same effects
  expect_equal(length(fit_b$sets$cs), length(fit_p$sets$cs))
  # c_hat values should be similar (binomial correction is small for L >> C)
  expect_equal(fit_b$C_hat, fit_p$C_hat, tolerance = 1)
})

test_that("susie without slot_prior does not produce c_hat", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  fit <- susie(X, y, L = 5, verbose = FALSE)
  expect_null(fit$c_hat)
})

test_that("ash model auto-creates binomial slot_prior with warning", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  expect_message(
    fit <- susie(X, y, L = 10, unmappable_effects = "ash",
                 verbose = FALSE, max_iter = 5),
    "strongly advised"
  )
  expect_true(!is.null(fit$c_hat))
})

test_that("ash model with explicit slot_prior does not warn about C", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  # Should not produce the "strongly advised" warning about C
  # (may still produce convergence method warning, which is expected)
  fit <- withCallingHandlers(
    susie(X, y, L = 10, unmappable_effects = "ash",
          slot_prior = slot_prior_binomial(C = 3, nu = 8),
          verbose = FALSE, max_iter = 5),
    warning = function(w) {
      if (grepl("strongly advised", conditionMessage(w)))
        stop("Got unexpected C warning")
      invokeRestart("muffleWarning")
    }
  )
  expect_true(!is.null(fit$c_hat))
})

test_that("c_hat warm start works", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[1:3] <- 1
  y <- X %*% b + rnorm(n)

  # First fit
  fit1 <- susie(X, y, L = 10, slot_prior = slot_prior_binomial(C = 3, nu = 8),
                verbose = FALSE)

  # Warm start with previous c_hat
  sp_warm <- slot_prior_binomial(C = 3, nu = 8, c_hat_init = fit1$c_hat)
  fit2 <- susie(X, y, L = 10, slot_prior = sp_warm,
                model_init = fit1, verbose = FALSE)

  # Should converge immediately or in very few iterations
  expect_true(fit2$niter <= fit1$niter)
})

test_that("batch and sequential schedules both converge", {
  set.seed(1)
  n <- 100; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[1:3] <- 1
  y <- X %*% b + rnorm(n)
  fit_batch <- susie(X, y, L = 10,
                     slot_prior = slot_prior_binomial(C = 3, update_schedule = "batch"),
                     verbose = FALSE)
  fit_seq <- susie(X, y, L = 10,
                   slot_prior = slot_prior_binomial(C = 3, update_schedule = "sequential"),
                   verbose = FALSE)
  expect_true(fit_batch$converged)
  expect_true(fit_seq$converged)
  # Results should be very similar
  expect_equal(sum(fit_batch$c_hat > 0.5), sum(fit_seq$c_hat > 0.5))
})
