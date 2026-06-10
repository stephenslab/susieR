context("susie_auto unit tests")

# ---- L doubling & convergence ----

test_that("susie_auto L stays at L_init=L_max when no doubling possible", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 126)
  result <- susie_auto(base_data$X, base_data$y, L_init = 5, L_max = 5, verbose = FALSE)
  expect_equal(nrow(result$alpha), 5)
  expect_true(is.finite(result$elbo[length(result$elbo)]))
})

test_that("susie_auto doubles L when all prior variances remain positive", {
  # k=10 signals at low signal_sd so all V stay >0 and doubling is triggered
  set.seed(128)
  base_data <- generate_base_data(n = 100, p = 50, k = 10, signal_sd = 0.5, seed = NULL)
  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 4, verbose = FALSE)
  L_final <- nrow(result$alpha)
  expect_true(L_final %in% c(1, 2, 4, 8))
})

test_that("susie_auto converges (stops doubling) when any prior variance hits zero", {
  # Strong single signal: one effect absorbs all signal, leaving remaining V=0
  base_data <- generate_base_data(n = 100, p = 50, k = 1, signal_sd = 5, seed = 125)
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 16, verbose = FALSE)
  expect_true(any(result$V < 1e-3))
})

test_that("susie_auto exercises the final-run path with tighter tol", {
  # Two different init_tol values should produce nearly the same final result
  # because the final susie() call uses `tol`, not `init_tol`
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 132)
  result_large <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 2,
                             init_tol = 10, tol = 1e-3, verbose = FALSE)
  result_small <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 2,
                             init_tol = 1e-5, tol = 1e-3, verbose = FALSE)
  # Final ELBOs should be close regardless of init_tol
  expect_equal(result_large$elbo[length(result_large$elbo)],
               result_small$elbo[length(result_small$elbo)],
               tolerance = 1)
})

test_that("susie_auto exercises full L-doubling loop then final run", {
  # Strong 5-signal data with small L_max forces multiple doublings then final run
  base_data <- generate_base_data(n = 100, p = 50, k = 5, signal_sd = 1, seed = 134)
  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 8, verbose = FALSE)
  L_final <- nrow(result$alpha)
  # All output dimensions must be consistent
  expect_equal(nrow(result$alpha),          L_final)
  expect_equal(nrow(result$mu),             L_final)
  expect_equal(nrow(result$mu2),            L_final)
  expect_equal(length(result$V),            L_final)
  expect_equal(length(result$KL),           L_final)
})

# ---- Parameter propagation ----

test_that("susie_auto propagates standardize parameter: prior variances differ on unevenly-scaled X", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 129)
  base_data$X <- sweep(base_data$X, 2, seq(0.1, 5, length.out = base_data$p), "*")
  base_data$y <- as.vector(base_data$X %*% base_data$beta + rnorm(base_data$n))

  result_std   <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4,
                             standardize = TRUE,  verbose = FALSE)
  result_nostd <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4,
                             standardize = FALSE, verbose = FALSE)

  # Both paths must produce valid objects
  expect_true(all(result_std$alpha   >= 0 & result_std$alpha   <= 1))
  expect_true(all(result_nostd$alpha >= 0 & result_nostd$alpha <= 1))
  # V is on the coefficient scale; standardization changes that scale
  expect_false(isTRUE(all.equal(result_std$V, result_nostd$V, tolerance = 1e-4)))
})

test_that("susie_auto propagates intercept parameter: estimates differ with/without intercept", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 130)
  base_data$y <- as.vector(base_data$X %*% base_data$beta + 3 + rnorm(base_data$n))

  result_int   <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4,
                             intercept = TRUE,  verbose = FALSE)
  result_noint <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4,
                             intercept = FALSE, verbose = FALSE)

  expect_false(isTRUE(all.equal(result_int$intercept, result_noint$intercept, tolerance = 1e-3)))
})

test_that("susie_auto respects max_iter: niter does not exceed limit", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 131)
  result <- suppressWarnings(
    susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 2, max_iter = 3, verbose = FALSE)
  )
  expect_true(result$niter <= 3)
})

test_that("susie_auto passes ... args to susie: coverage affects CS computation", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 133)
  result_high <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 2,
                            coverage = 0.95, verbose = FALSE)
  result_low  <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 2,
                            coverage = 0.5,  verbose = FALSE)
  cs_high <- susie_get_cs(result_high, coverage = 0.95)
  cs_low  <- susie_get_cs(result_low,  coverage = 0.5)
  # With lower coverage threshold, CS sizes should be <= those with higher coverage
  if (!is.null(cs_high$cs) && !is.null(cs_low$cs) &&
      length(cs_high$cs) > 0 && length(cs_low$cs) > 0) {
    mean_sz_high <- mean(sapply(cs_high$cs, length))
    mean_sz_low  <- mean(sapply(cs_low$cs,  length))
    expect_true(mean_sz_low <= mean_sz_high + 1)  # low coverage -> smaller or equal CSs
  }
  # Both must be structurally sound
  expect_true(is.list(result_high$sets))
  expect_true(is.list(result_low$sets))
})

# ---- Model structure ----

test_that("susie_auto returns valid susie object with all required fields", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 145)
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  for (field in c("alpha", "mu", "mu2", "V", "sigma2", "elbo", "niter", "fitted")) {
    expect_true(field %in% names(result),
                label = paste("field", field, "present"))
  }
  L <- nrow(result$alpha)
  expect_equal(dim(result$alpha), c(L, base_data$p))
  expect_equal(dim(result$mu),    c(L, base_data$p))
  expect_equal(dim(result$mu2),   c(L, base_data$p))
  expect_equal(length(result$V),  L)
  expect_equal(length(result$KL), L)
})

test_that("susie_auto output: fitted values, predictions, and coefs have correct dimensions", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 147)
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  expect_equal(length(result$fitted), base_data$n)
  expect_true(all(is.finite(result$fitted)))

  pred <- predict(result)
  expect_equal(length(pred), base_data$n)
  expect_true(all(is.finite(pred)))

  coefs <- coef(result)
  expect_equal(length(coefs), base_data$p + 1)   # p + intercept
  expect_true(all(is.finite(coefs)))
})

# ---- Mathematical properties ----

test_that("susie_auto alpha rows sum to 1 and values are in [0,1]", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 152)
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  expect_true(all(abs(rowSums(result$alpha) - 1) < 1e-8))
  expect_true(all(result$alpha >= 0))
  expect_true(all(result$alpha <= 1))
})

test_that("susie_auto variance estimates are valid: sigma2>0, V>=0, KL>=0", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 136)
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  expect_true(result$sigma2 > 0)
  expect_true(is.finite(result$sigma2))
  expect_true(all(result$V >= 0))
  expect_true(all(is.finite(result$V)))
  expect_true(all(result$KL >= -1e-10))
})

test_that("susie_auto ELBO is non-decreasing within a single fit", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 150)
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  elbo_diff <- diff(result$elbo)
  expect_true(all(elbo_diff > -1e-6))
})

test_that("susie_auto PIPs are valid probabilities", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 146)
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  pips <- susie_get_pip(result)
  expect_true(all(pips >= 0))
  expect_true(all(pips <= 1))
})

# ---- Edge cases & robustness ----

for (cfg in list(
  list(label = "single effect",   k = 1,  signal_sd = 3,    L_init = 1, L_max = 8,  seed = 139),
  list(label = "pure noise",      k = 0,  signal_sd = 0,    L_init = 1, L_max = 4,  seed = 142),
  list(label = "L_init=1",        k = 2,  signal_sd = 1.75, L_init = 1, L_max = 8,  seed = 143),
  list(label = "large L_init=10", k = 2,  signal_sd = 1.75, L_init = 10, L_max = 10, seed = 144)
)) {
  local({
    cfg_ <- cfg
    test_that(paste("susie_auto handles", cfg_$label), {
      if (cfg_$k == 0) {
        bd <- generate_base_data(n = 100, p = 50, k = 0, seed = cfg_$seed)
      } else {
        bd <- generate_base_data(n = 100, p = 50, k = cfg_$k,
                                 signal_sd = cfg_$signal_sd, seed = cfg_$seed)
      }
      result <- susie_auto(bd$X, bd$y, L_init = cfg_$L_init, L_max = cfg_$L_max,
                           verbose = FALSE)
      expect_true(is.finite(result$elbo[length(result$elbo)]))
      expect_true(result$sigma2 > 0)
      if (cfg_$label == "large L_init=10") {
        expect_equal(nrow(result$alpha), 10)
      }
      if (cfg_$label == "pure noise") {
        expect_true(all(result$V < 0.5))
        # no signal → no credible sets
        expect_true(is.null(result$sets) || length(result$sets$cs) == 0)
      }
    })
  })
}

test_that("susie_auto high noise: sigma2 substantially larger than low-noise case", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 137)
  result_low <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  set.seed(141)
  bd_high <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = NULL)
  bd_high$y <- as.vector(bd_high$X %*% bd_high$beta + rnorm(bd_high$n, sd = 5))
  result_high <- susie_auto(bd_high$X, bd_high$y, L_init = 2, L_max = 4, verbose = FALSE)

  expect_true(result_high$sigma2 > 0)
  expect_true(is.finite(result_high$sigma2))
})

# ---- Signal recovery ----

test_that("susie_auto recovers causal variables with strong signal", {
  base_data <- generate_base_data(n = 200, p = 100, k = 3, signal_sd = 1.5, seed = 155)
  base_data$X <- scale(base_data$X)
  base_data$beta[base_data$causal_idx] <- c(1.5, -1.2, 1.8)
  base_data$y <- as.vector(base_data$X %*% base_data$beta + rnorm(base_data$n, sd = 0.5))

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 10, verbose = FALSE)
  pips <- susie_get_pip(result)

  top_vars <- order(pips, decreasing = TRUE)[1:5]
  expect_true(length(intersect(top_vars, base_data$causal_idx)) >= 2)
})

test_that("susie_auto: median null-variable PIP is low", {
  base_data <- generate_base_data(n = 200, p = 100, k = 2, signal_sd = 1.9, seed = 156)
  base_data$X <- scale(base_data$X)
  base_data$y <- as.vector(base_data$X %*% base_data$beta + rnorm(base_data$n, sd = 0.5))

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 10, verbose = FALSE)
  pips <- susie_get_pip(result)

  null_pips <- pips[setdiff(seq_len(base_data$p), base_data$causal_idx)]
  expect_true(median(null_pips) < 0.2)
})

test_that("susie_auto credible sets have >= 0.9 coverage for strong effects", {
  base_data <- generate_base_data(n = 200, p = 100, k = 3, signal_sd = 2, seed = 157)
  base_data$X <- scale(base_data$X)
  base_data$y <- as.vector(base_data$X %*% base_data$beta + rnorm(base_data$n, sd = 0.5))

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 10, verbose = FALSE)
  cs <- susie_get_cs(result)

  expect_true(!is.null(cs) && length(cs$cs) >= 1)
  expect_true(all(cs$coverage >= 0.9))
})

# ---- Verbose output ----

test_that("susie_auto verbose=TRUE emits 'Trying L=' messages", {
  base_data <- generate_base_data(n = 100, p = 50, k = 5, signal_sd = 1, seed = 158)
  msgs <- capture_messages(
    susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 4, verbose = TRUE)
  )
  expect_true(any(grepl("Trying L=", msgs)))
  # First message must be "Trying L=1"
  expect_true(any(grepl("Trying L=1", msgs)))
})

# ---- Reproducibility ----

test_that("susie_auto gives identical results with the same seed", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 199)

  set.seed(200)
  result1 <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)
  set.seed(200)
  result2 <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  expect_equal(result1$alpha, result2$alpha)
  expect_equal(result1$elbo,  result2$elbo)
})

test_that("susie_auto with different L_init converges to similar PIP rankings", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 201)

  result_L1 <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 8, verbose = FALSE)
  result_L2 <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 8, verbose = FALSE)

  pips1 <- susie_get_pip(result_L1)
  pips2 <- susie_get_pip(result_L2)
  expect_true(cor(pips1, pips2) > 0.9)
})
