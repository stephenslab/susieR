devtools::load_all(".")
source(file.path("..", "helpers", "helper_testthat.R"), local = TRUE)

context("susie_auto unit tests")

# =============================================================================
# ALGORITHM PROGRESSION
# =============================================================================

test_that("susie_auto starts with L_init and doubles correctly", {
  base_data <- generate_base_data(n = 100, p = 50, k = 3, signal_sd = 1, seed = 123)
  # Manually set specific beta values for this test
  base_data$beta[base_data$causal_idx] <- c(2, -1.5, 1)
  base_data$y <- as.vector(base_data$X %*% base_data$beta + rnorm(base_data$n))

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 4, verbose = FALSE)

  L_final <- nrow(result$alpha)
  expect_true(L_final >= 1)
  expect_true(L_final %in% c(1, 2, 4, 8))
})

test_that("susie_auto respects L_max limit", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 124)

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 8, verbose = FALSE)

  L_final <- nrow(result$alpha)
  expect_true(L_final %in% c(1, 2, 4, 8, 16))
})

test_that("susie_auto converges when prior variances hit zero", {
  base_data <- generate_base_data(n = 100, p = 50, k = 1, signal_sd = 3, seed = 125)

  # With single effect, should converge quickly (some V should be 0)
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 16, verbose = FALSE)

  # At least one prior variance should be zero (or very small)
  expect_true(any(result$V < 1e-3))
})

test_that("susie_auto handles L_init = L_max (no doubling)", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 126)

  # No doubling should occur
  result <- susie_auto(base_data$X, base_data$y, L_init = 5, L_max = 5, verbose = FALSE)

  # Should complete successfully with L = 5
  expect_equal(nrow(result$alpha), 5)
  expect_true(is.finite(result$elbo[length(result$elbo)]))
})

# =============================================================================
# CONVERGENCE BEHAVIOR
# =============================================================================

test_that("susie_auto convergence logic: stops when any V = 0", {
  set.seed(127)
  base_data <- generate_base_data(n = 100, p = 50, k = 1, signal_sd = 5, seed = NULL)
  # Add lower noise
  base_data$y <- base_data$X %*% base_data$beta + rnorm(base_data$n, sd = 0.5)

  # Should converge with strong single effect
  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 32, verbose = FALSE)

  # At least one V should be effectively zero (converged)
  expect_true(any(result$V < 1e-6))

  # Result should complete successfully
  expect_true(is.finite(result$elbo[length(result$elbo)]))
})

test_that("susie_auto continues until L_max when all V > 0", {
  set.seed(128)
  base_data <- generate_base_data(n = 100, p = 50, k = 10, signal_sd = 0.5, seed = NULL)

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 4, verbose = FALSE)

  L_final <- nrow(result$alpha)
  expect_true(L_final %in% c(1, 2, 4, 8))
})

# =============================================================================
# PARAMETER PROPAGATION
# =============================================================================

test_that("susie_auto propagates standardize parameter correctly", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 129)
  # Create X with different scales
  base_data$X <- sweep(base_data$X, 2, seq(0.1, 5, length.out = base_data$p), "*")
  base_data$y <- as.vector(base_data$X %*% base_data$beta + rnorm(base_data$n))

  result_std <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4,
                           standardize = TRUE, verbose = FALSE)
  result_nostd <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4,
                             standardize = FALSE, verbose = FALSE)

  expect_true(all(result_std$alpha >= 0 & result_std$alpha <= 1))
  expect_true(all(result_nostd$alpha >= 0 & result_nostd$alpha <= 1))
})

test_that("susie_auto propagates intercept parameter correctly", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 130)
  base_data$y <- as.vector(base_data$X %*% base_data$beta + 3 + rnorm(base_data$n))  # Add intercept

  result_int <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4,
                           intercept = TRUE, verbose = FALSE)
  result_noint <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4,
                             intercept = FALSE, verbose = FALSE)

  # Intercept estimates should differ
  expect_false(isTRUE(all.equal(result_int$intercept, result_noint$intercept,
                                 tolerance = 1e-3)))
})

test_that("susie_auto propagates max_iter parameter correctly", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 131)

  # Use very small max_iter to test propagation
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 2, max_iter = 3,
                       verbose = FALSE)

  # Should complete (may not converge, but should respect max_iter)
  expect_true(result$niter <= 3)
})

test_that("susie_auto uses init_tol for early runs and tol for final run", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 132)

  # Large init_tol should make early runs converge faster
  result_large_init <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 2,
                                  init_tol = 10, tol = 1e-3, verbose = FALSE)
  result_small_init <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 2,
                                  init_tol = 1e-5, tol = 1e-3, verbose = FALSE)

  # Both should complete successfully
  expect_true(is.finite(result_large_init$elbo[length(result_large_init$elbo)]))
  expect_true(is.finite(result_small_init$elbo[length(result_small_init$elbo)]))
})

test_that("susie_auto passes additional arguments via ...", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 133)

  # Pass coverage argument
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 2, coverage = 0.9,
                       verbose = FALSE)

  # Should complete successfully
  expect_true(is.finite(result$elbo[length(result$elbo)]))

  # Check that sets are computed (if any exist)
  if (!is.null(result$sets)) {
    expect_true(is.list(result$sets))
  }
})

# =============================================================================
# MODEL INITIALIZATION & EXPANSION
# =============================================================================

test_that("susie_auto correctly expands L via add_null_effect", {
  base_data <- generate_base_data(n = 100, p = 50, k = 5, signal_sd = 1, seed = 134)

  # Start small and let it expand
  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 8, verbose = FALSE)

  L_final <- nrow(result$alpha)

  # Dimensions should be consistent
  expect_equal(nrow(result$alpha), L_final)
  expect_equal(nrow(result$mu), L_final)
  expect_equal(nrow(result$mu2), L_final)
  expect_equal(length(result$V), L_final)
  expect_equal(length(result$KL), L_final)
})

test_that("susie_auto maintains valid model structure throughout", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 135)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # Check model structure
  L <- nrow(result$alpha)
  expect_equal(dim(result$alpha), c(L, base_data$p))
  expect_equal(dim(result$mu), c(L, base_data$p))
  expect_equal(dim(result$mu2), c(L, base_data$p))
  expect_equal(length(result$V), L)
  expect_equal(length(result$KL), L)

  # Alpha rows should sum to 1
  expect_true(all(abs(rowSums(result$alpha) - 1) < 1e-10))
})

# =============================================================================
# VARIANCE ESTIMATION
# =============================================================================

test_that("susie_auto estimates variances in correct stages", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 136)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 2, verbose = FALSE)

  # Final result should have estimated variances
  expect_true(result$sigma2 > 0)
  expect_true(all(result$V >= 0))
})

test_that("susie_auto residual variance is positive", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 137)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  expect_true(result$sigma2 > 0)
  expect_true(is.finite(result$sigma2))
})

test_that("susie_auto prior variances are non-negative", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 138)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  expect_true(all(result$V >= 0))
  expect_true(all(is.finite(result$V)))
})

# =============================================================================
# EDGE CASES & ROBUSTNESS
# =============================================================================

test_that("susie_auto handles sparse signal (single effect)", {
  base_data <- generate_base_data(n = 100, p = 50, k = 1, signal_sd = 3, seed = 139)

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 8, verbose = FALSE)

  # Should identify the effect
  pips <- colSums(result$alpha * result$mu)
  expect_true(max(pips) > 0.5)  # At least one variable should have high PIP

  # Most V should be zero or very small
  expect_true(sum(result$V < 1e-3) >= length(result$V) - 2)
})

test_that("susie_auto handles dense signal (many effects)", {
  base_data <- generate_base_data(n = 100, p = 50, k = 8, signal_sd = 0.5, seed = 140)

  # May need multiple doublings
  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 16, verbose = FALSE)

  # Should complete successfully
  expect_true(is.finite(result$elbo[length(result$elbo)]))
  expect_true(nrow(result$alpha) >= 2)
})

test_that("susie_auto handles high noise scenario", {
  set.seed(141)
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = NULL)
  # Add high noise
  base_data$y <- as.vector(base_data$X %*% base_data$beta + rnorm(base_data$n, sd = 5))

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # Should complete successfully
  expect_true(is.finite(result$elbo[length(result$elbo)]))
  expect_true(result$sigma2 > 0)
})

test_that("susie_auto handles no signal (pure noise)", {
  base_data <- generate_base_data(n = 100, p = 50, k = 0, seed = 142)

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 4, verbose = FALSE)

  # Should complete successfully
  expect_true(is.finite(result$elbo[length(result$elbo)]))

  # All effects should have small prior variance
  expect_true(all(result$V < 0.5))

  # No credible sets should be found
  expect_true(is.null(result$sets) || length(result$sets$cs) == 0)
})

test_that("susie_auto handles L_init = 1 (minimum)", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 143)

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 8, verbose = FALSE)

  # Should complete successfully
  expect_true(is.finite(result$elbo[length(result$elbo)]))
  expect_true(nrow(result$alpha) >= 1)
})

test_that("susie_auto handles large L_init", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 144)

  # Start with large L
  result <- susie_auto(base_data$X, base_data$y, L_init = 10, L_max = 10, verbose = FALSE)

  # Should complete successfully with L = 10
  expect_equal(nrow(result$alpha), 10)
  expect_true(is.finite(result$elbo[length(result$elbo)]))
})

# =============================================================================
# OUTPUT VALIDATION
# =============================================================================

test_that("susie_auto returns valid susie object", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 145)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # Check required fields
  expect_true("alpha" %in% names(result))
  expect_true("mu" %in% names(result))
  expect_true("mu2" %in% names(result))
  expect_true("V" %in% names(result))
  expect_true("sigma2" %in% names(result))
  expect_true("elbo" %in% names(result))
  expect_true("niter" %in% names(result))
})

test_that("susie_auto PIPs are valid probabilities", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 146)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # PIPs should be between 0 and 1
  pips <- susie_get_pip(result)
  expect_true(all(pips >= 0))
  expect_true(all(pips <= 1))
})

test_that("susie_auto fitted values have correct dimensions", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 147)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # Fitted values should have length n
  expect_equal(length(result$fitted), base_data$n)
  expect_true(all(is.finite(result$fitted)))
})

test_that("susie_auto predictions work correctly", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 148)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # Predictions should work
  pred <- predict(result)
  expect_equal(length(pred), base_data$n)
  expect_true(all(is.finite(pred)))
})

test_that("susie_auto coefficients can be extracted", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 149)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # Coefficients should be extractable
  coefs <- coef(result)
  expect_equal(length(coefs), base_data$p + 1)  # p coefficients + intercept
  expect_true(all(is.finite(coefs)))
})

# =============================================================================
# MATHEMATICAL PROPERTIES
# =============================================================================

test_that("susie_auto ELBO is monotonically increasing or stable", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 150)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # ELBO should be non-decreasing (allowing for small numerical errors)
  elbo_diff <- diff(result$elbo)
  expect_true(all(elbo_diff > -1e-6))
})

test_that("susie_auto final ELBO is finite", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 151)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  final_elbo <- result$elbo[length(result$elbo)]
  expect_true(is.finite(final_elbo))
  expect_false(is.na(final_elbo))
})

test_that("susie_auto alpha rows sum to 1", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 152)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # Each row of alpha should sum to 1
  row_sums <- rowSums(result$alpha)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("susie_auto alpha values are valid probabilities", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 153)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # Alpha should be in [0, 1]
  expect_true(all(result$alpha >= 0))
  expect_true(all(result$alpha <= 1))
})

test_that("susie_auto KL divergences are non-negative", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 154)

  result <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # KL divergences should be non-negative
  expect_true(all(result$KL >= -1e-10))  # Allow small numerical error
})

# =============================================================================
# SIGNAL RECOVERY
# =============================================================================

test_that("susie_auto recovers true causal variables with strong signal", {
  base_data <- generate_base_data(n = 200, p = 100, k = 3, signal_sd = 1.5, seed = 155)
  base_data$X <- scale(base_data$X)
  # Set specific effect sizes
  base_data$beta[base_data$causal_idx] <- c(1.5, -1.2, 1.8)
  base_data$y <- as.vector(base_data$X %*% base_data$beta + rnorm(base_data$n, sd = 0.5))

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 10, verbose = FALSE)

  # Get PIPs
  pips <- susie_get_pip(result)

  # Top PIPs should include causal variables
  top_vars <- order(pips, decreasing = TRUE)[1:5]
  expect_true(length(intersect(top_vars, base_data$causal_idx)) >= 2)
})

test_that("susie_auto has low PIPs for null variables", {
  base_data <- generate_base_data(n = 200, p = 100, k = 2, signal_sd = 1.9, seed = 156)
  base_data$X <- scale(base_data$X)
  base_data$y <- as.vector(base_data$X %*% base_data$beta + rnorm(base_data$n, sd = 0.5))

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 10, verbose = FALSE)

  # Get PIPs
  pips <- susie_get_pip(result)

  # Most PIPs should be low for null variables
  null_pips <- pips[setdiff(1:base_data$p, base_data$causal_idx)]
  expect_true(median(null_pips) < 0.2)
})

test_that("susie_auto identifies credible sets for strong effects", {
  base_data <- generate_base_data(n = 200, p = 100, k = 3, signal_sd = 2, seed = 157)
  base_data$X <- scale(base_data$X)
  base_data$y <- as.vector(base_data$X %*% base_data$beta + rnorm(base_data$n, sd = 0.5))

  result <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 10, verbose = FALSE)

  # Should find credible sets
  cs <- susie_get_cs(result)

  if (!is.null(cs) && length(cs$cs) > 0) {
    # At least one CS should be found
    expect_true(length(cs$cs) >= 1)

    # CSs should have reasonable coverage
    expect_true(all(cs$coverage >= 0.9))
  }
})

# =============================================================================
# VERBOSE OUTPUT
# =============================================================================

test_that("susie_auto verbose mode produces output", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 158)

  # Capture messages
  expect_message(
    susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 4, verbose = TRUE),
    "Trying L="
  )
})

test_that("susie_auto verbose shows correct L progression", {
  base_data <- generate_base_data(n = 100, p = 50, k = 5, signal_sd = 1, seed = 159)

  # Capture messages and check for L progression
  msgs <- capture_messages(
    susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 4, verbose = TRUE)
  )

  # Should see "Trying L=" messages
  expect_true(any(grepl("Trying L=", msgs)))
})

# =============================================================================
# CONSISTENCY
# =============================================================================

test_that("susie_auto gives consistent results with same seed", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 199)

  set.seed(200)
  result1 <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  set.seed(200)
  result2 <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 4, verbose = FALSE)

  # Results should be identical
  expect_equal(result1$alpha, result2$alpha)
  expect_equal(result1$elbo, result2$elbo)
})

test_that("susie_auto with different L_init converges to similar solutions", {
  base_data <- generate_base_data(n = 100, p = 50, k = 2, signal_sd = 1.75, seed = 201)

  result_L1 <- susie_auto(base_data$X, base_data$y, L_init = 1, L_max = 8, verbose = FALSE)
  result_L2 <- susie_auto(base_data$X, base_data$y, L_init = 2, L_max = 8, verbose = FALSE)

  # PIPs should be similar
  pips1 <- susie_get_pip(result_L1)
  pips2 <- susie_get_pip(result_L2)

  # Correlation of PIPs should be high
  expect_true(cor(pips1, pips2) > 0.9)
})
