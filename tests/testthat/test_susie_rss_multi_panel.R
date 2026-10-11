# =============================================================================
# Multi-panel two-stage fit: omega from the joint fit, then refit on R(omega)
# =============================================================================

test_that("multi-panel R list: two-stage fit equals a single fit on sum_k omega_k R_k", {
  set.seed(20)
  p <- 15; Bn <- 80
  X1 <- matrix(rnorm(Bn * p), Bn, p)
  X2 <- matrix(rnorm(Bn * p), Bn, p)
  R1 <- cor(X1); R2 <- cor(X2)
  z <- rnorm(p); z[4] <- 5
  fit <- suppressWarnings(susie_rss(z = z, R = list(R1, R2), n = 1000,
                                    L = 3, max_iter = 20,
                                    R_finite = c(80, 120),
                                    R_mismatch = "eb_mix"))
  # omega is the joint mixture fit's omega (stage 1).
  obj <- suppressWarnings(susie_rss(z = z, R = list(R1, R2), n = 1000,
                                    L = 3, max_iter = 20,
                                    R_finite = c(80, 120),
                                    R_mismatch = "eb_mix", init_only = TRUE))
  joint <- suppressWarnings(susie_workhorse(obj$data, obj$params))
  om <- fit$omega_weights
  expect_equal(om, joint$omega_weights)
  expect_equal(sum(om), 1, tolerance = 1e-12)
  expect_false(is.null(fit$single_panel_fits))

  direct <- suppressWarnings(susie_rss(z = z, R = om[1] * R1 + om[2] * R2,
                                       n = 1000, L = 3, max_iter = 20,
                                       R_finite = 1 / sum(om^2 / c(80, 120)),
                                       R_mismatch = "eb_mix"))
  expect_equal(fit$pip, direct$pip)
  expect_equal(fit$alpha, direct$alpha)
  expect_equal(fit$R_finite_diagnostics$B, 1 / sum(om^2 / c(80, 120)))
  expect_equal(fit$R_finite_diagnostics$lambda_bias,
               direct$R_finite_diagnostics$lambda_bias)
})

test_that("multi-panel X list: weighted sketch reproduces sum_k omega_k R_k", {
  set.seed(21)
  p <- 40; B1 <- 12; B2 <- 9
  X1 <- matrix(rnorm(B1 * p, mean = 2, sd = 3), B1, p)
  X2 <- matrix(rnorm(B2 * p), B2, p)
  om <- c(0.3, 0.7)
  ref <- form_weighted_reference(list(X1, X2), "X", om)
  expect_equal(nrow(ref$X), B1 + B2)
  expect_equal(crossprod(ref$X), om[1] * cor(X1) + om[2] * cor(X2))
  # A zero weight drops that panel's rows.
  expect_equal(nrow(form_weighted_reference(list(X1, X2), "X", c(0, 1))$X), B2)

  z <- rnorm(p); z[7] <- 4.5
  fit <- suppressWarnings(susie_rss(z = z, X = list(X1, X2), n = 500, L = 2,
                                    max_iter = 20, R_finite = TRUE,
                                    R_mismatch = "eb"))
  om <- fit$omega_weights
  ref <- form_weighted_reference(list(X1, X2), "X", om)
  direct <- suppressWarnings(susie_rss(z = z, X = ref$X, n = 500, L = 2,
                                       max_iter = 20,
                                       R_finite = 1 / sum(om^2 / c(B1, B2)),
                                       R_mismatch = "eb"))
  expect_equal(fit$pip, direct$pip)
  expect_equal(fit$R_finite_diagnostics$B, 1 / sum(om^2 / c(B1, B2)))
})

test_that("multi-panel verbose message reports each reference panel's weight", {
  set.seed(22)
  p <- 15; Bn <- 80
  R1 <- cor(matrix(rnorm(Bn * p), Bn, p))
  R2 <- cor(matrix(rnorm(Bn * p), Bn, p))
  z <- rnorm(p); z[4] <- 5
  msgs <- character(0)
  fit <- withCallingHandlers(
    suppressWarnings(susie_rss(z = z, R = list(EUR = R1, AFR = R2), n = 1000,
                               L = 3, max_iter = 20, verbose = TRUE)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  w <- grep("Multi-panel weights", msgs, value = TRUE)
  expect_length(w, 1)
  expect_match(w, sprintf("the weight of reference panel 'EUR' is %.3f",
                          fit$omega_weights[1]), fixed = TRUE)
  expect_match(w, sprintf("the weight of reference panel 'AFR' is %.3f",
                          fit$omega_weights[2]), fixed = TRUE)
  expect_false(any(grepl("B_eff", msgs)))
})

test_that("multi-panel: only the fit on the mixed LD reports IBSS non-convergence", {
  set.seed(23)
  p <- 40
  R1 <- cor(matrix(rnorm(300 * p), 300, p))
  R2 <- cor(matrix(rnorm(200 * p), 200, p))
  z <- rnorm(p); z[c(5, 20)] <- 7
  for (v in c(FALSE, TRUE)) {
    msgs <- character(0); warns <- character(0)
    fit <- withCallingHandlers(
      susie_rss(z = z, R = list(R1, R2), n = 5000, L = 10, max_iter = 2,
                verbose = v),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    expect_false(fit$converged)
    expect_equal(sum(grepl("IBSS algorithm did not converge", msgs)), 1)
    expect_equal(sum(grepl("IBSS algorithm did not converge", warns)), 1)
    # The per-panel fits keep their own convergence status.
    expect_false(any(vapply(fit$single_panel_fits, `[[`, logical(1),
                            "converged")))
  }
})
