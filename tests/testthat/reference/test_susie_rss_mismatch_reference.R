# Source helper functions
source(file.path("..", "helper_mismatch_reference.R"), local = TRUE)

context("susie_rss RSS-mismatch reference comparison")

# =============================================================================
# REFERENCE TESTS FOR susie_rss / susie_rss_lambda — RSS MISMATCH SURFACE
# =============================================================================
#
# These tests compare the current susie_rss / susie_rss_lambda implementation
# against the reference commit alexmccreight/susieR@2a86b39, which is the
# pre-refactor snapshot of this codebase used to verify that the thin-wrapper
# refactor of susie_rss preserves machine-precision results.
#
# Coverage:
#   1. The 8 vignette fits from vignettes/rss_mismatch.Rmd
#      (3 from the SummaryConsistency toy example,
#       4 from the rss_mismatch_example full example,
#       plus susie_rss_lambda smoke).
#   2. Multi-panel R and X (K = 2, K = 3) with various R_finite shapes.
#   3. All R_mismatch enum values × all R_finite shapes.
#   4. init_only paths (single-panel and multi-panel).
#   5. susie_rss_lambda across lambda values and L.
#
# Each test calls compare_to_mismatch_reference() which runs the same arguments
# against both the reference and dev environments and asserts machine-precision
# equality on every fit field (alpha, mu, mu2, V, sigma2, elbo, pip, niter,
# sets, R_finite_diagnostics, R-mismatch diagnostics, omega_weights, etc.).

# =============================================================================
# Part 1: Vignette toy example (SummaryConsistency dataset)
# =============================================================================

test_that("vignette toy: susie_rss(zflip, Rtoy, n) matches reference", {
  skip_if_no_mismatch_reference()
  data(SummaryConsistency, package = "susieR")
  with(SummaryConsistency, {
    args <- list(z = z, R = ldref, n = 10000)
    compare_to_mismatch_reference("susie_rss", args)
  })
})

test_that("vignette toy: susie_rss(zfix, Rtoy, n) matches reference", {
  skip_if_no_mismatch_reference()
  data(SummaryConsistency, package = "susieR")
  with(SummaryConsistency, {
    zfix <- z
    zfix[flip_id] <- -zfix[flip_id]
    args <- list(z = zfix, R = ldref, n = 10000)
    compare_to_mismatch_reference("susie_rss", args)
  })
})

test_that("vignette toy: susie_rss with R_finite = 500, R_mismatch = 'eb' matches reference", {
  skip_if_no_mismatch_reference()
  data(SummaryConsistency, package = "susieR")
  with(SummaryConsistency, {
    args <- list(z = z, R = ldref, n = 10000,
                 R_finite = 500, R_mismatch = "eb")
    compare_to_mismatch_reference("susie_rss", args)
  })
})

# =============================================================================
# Part 2: Vignette full example (rss_mismatch_example dataset)
# =============================================================================

test_that("vignette full: susie_rss baseline (no mismatch correction) matches reference", {
  skip_if_no_mismatch_reference()
  data(rss_mismatch_example, package = "susieR")
  with(rss_mismatch_example, {
    args <- list(bhat = bhat, shat = shat, R = R, n = n)
    compare_to_mismatch_reference("susie_rss", args)
  })
})

test_that("vignette full: susie_rss with R_finite = 500 matches reference", {
  skip_if_no_mismatch_reference()
  data(rss_mismatch_example, package = "susieR")
  with(rss_mismatch_example, {
    args <- list(bhat = bhat, shat = shat, R = R, n = n,
                 R_finite = 500)
    compare_to_mismatch_reference("susie_rss", args)
  })
})

test_that("vignette full: susie_rss with R_finite = 500 + R_mismatch = 'eb' matches reference", {
  skip_if_no_mismatch_reference()
  data(rss_mismatch_example, package = "susieR")
  with(rss_mismatch_example, {
    args <- list(bhat = bhat, shat = shat, R = R, n = n,
                 R_finite = 500, R_mismatch = "eb")
    compare_to_mismatch_reference("susie_rss", args)
  })
})

test_that("vignette full: susie_rss with R_mismatch = 'eb' only matches reference", {
  skip_if_no_mismatch_reference()
  data(rss_mismatch_example, package = "susieR")
  with(rss_mismatch_example, {
    args <- list(bhat = bhat, shat = shat, R = R, n = n,
                 R_mismatch = "eb")
    compare_to_mismatch_reference("susie_rss", args)
  })
})

# =============================================================================
# Part 3: R_mismatch × R_finite parameter grid (single-panel)
# =============================================================================
# Synthetic data with controlled seed so the reference and dev share an exact
# input. Covers every R_mismatch enum value × {NULL, FALSE, scalar} R_finite.

.mismatch_make_data <- function(seed = 202, n = 500, p = 80, k = 3) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, TRUE, TRUE)
  beta <- rep(0, p)
  beta[sort(sample(p, k))] <- rnorm(k, 0, 1.5)
  y <- drop(X %*% beta + rnorm(n))
  bhat <- drop(crossprod(X, y) / (n - 1))
  shat <- rep(sqrt(var(y) / (n - 1)), p)
  list(X = X, y = y,
       z = bhat / shat,
       R = cor(X),
       n = n, bhat = bhat, shat = shat)
}

for (rmm in c("eb", "eb_no_init", "eb_force_init")) {
  for (rfin_label in c("NULL", "FALSE", "5000")) {
    rfin_val <- switch(rfin_label, "NULL" = NULL, "FALSE" = FALSE,
                       "5000" = 5000)
    local({
      rmm_l <- rmm
      rfin_label_l <- rfin_label
      rfin_val_l <- rfin_val
      test_that(sprintf("R-grid: R_mismatch=%s R_finite=%s matches reference",
                        rmm_l, rfin_label_l), {
        skip_if_no_mismatch_reference()
        d <- .mismatch_make_data()
        args <- list(z = d$z, R = d$R, n = d$n, L = 10,
                     R_mismatch = rmm_l, R_finite = rfin_val_l,
                     max_iter = 50)
        compare_to_mismatch_reference("susie_rss", args)
      })
    })
  }
}

# =============================================================================
# Part 4: R_finite = TRUE on X path
# =============================================================================

test_that("X path: R_finite = TRUE matches reference", {
  skip_if_no_mismatch_reference()
  d <- .mismatch_make_data()
  args <- list(z = d$z, X = d$X, n = d$n, L = 10,
               R_finite = TRUE, max_iter = 50)
  compare_to_mismatch_reference("susie_rss", args)
})

# =============================================================================
# Part 5: Multi-panel R and X (K = 2, K = 3)
# =============================================================================

test_that("multi-panel R (K=2) matches reference", {
  skip_if_no_mismatch_reference()
  d <- .mismatch_make_data()
  set.seed(606); Rb <- cor(matrix(rnorm(500 * 80), 500, 80))
  args <- list(z = d$z, R = list(d$R, Rb), n = d$n, L = 10,
               max_iter = 30)
  compare_to_mismatch_reference("susie_rss", args)
})

test_that("multi-panel R (K=3) matches reference", {
  skip_if_no_mismatch_reference()
  d <- .mismatch_make_data()
  set.seed(606); Rb <- cor(matrix(rnorm(500 * 80), 500, 80))
  set.seed(666); Rc <- cor(matrix(rnorm(500 * 80), 500, 80))
  args <- list(z = d$z, R = list(d$R, Rb, Rc), n = d$n, L = 10,
               max_iter = 30)
  compare_to_mismatch_reference("susie_rss", args)
})

test_that("multi-panel X (K=2) matches reference", {
  skip_if_no_mismatch_reference()
  d <- .mismatch_make_data()
  set.seed(707); Xb <- scale(matrix(rnorm(500 * 80), 500, 80), TRUE, TRUE)
  args <- list(z = d$z, X = list(d$X, Xb), n = d$n, L = 10,
               max_iter = 30)
  compare_to_mismatch_reference("susie_rss", args)
})

test_that("multi-panel R with R_finite vector matches reference", {
  skip_if_no_mismatch_reference()
  d <- .mismatch_make_data()
  set.seed(606); Rb <- cor(matrix(rnorm(500 * 80), 500, 80))
  args <- list(z = d$z, R = list(d$R, Rb), n = d$n, L = 10,
               R_finite = c(5000, 4000), R_mismatch = "eb",
               max_iter = 30)
  compare_to_mismatch_reference("susie_rss", args)
})

test_that("multi-panel R with R_finite scalar broadcast matches reference", {
  skip_if_no_mismatch_reference()
  d <- .mismatch_make_data()
  set.seed(606); Rb <- cor(matrix(rnorm(500 * 80), 500, 80))
  args <- list(z = d$z, R = list(d$R, Rb), n = d$n, L = 10,
               R_finite = 4500, R_mismatch = "eb",
               max_iter = 30)
  compare_to_mismatch_reference("susie_rss", args)
})

# =============================================================================
# Part 6: init_only paths
# =============================================================================
# Single-panel init_only is supported in both reference and dev. Multi-panel
# init_only with the old reference would crash (latent bug), so it is not
# included as a reference-comparison test — its post-refactor behavior is
# tested at the constructor level in test_susie_constructors.R.

test_that("init_only single-panel R matches reference", {
  skip_if_no_mismatch_reference()
  ref_env <- load_mismatch_reference_env()
  dev_env <- load_mismatch_development_env()
  d <- .mismatch_make_data()
  args <- list(z = d$z, R = d$R, n = d$n, L = 10,
               init_only = TRUE, max_iter = 50)
  dev_res <- suppressMessages(do.call(dev_env$env$susie_rss, args))
  ref_res <- suppressMessages(do.call(ref_env$env$susie_rss, args))
  # init_only returns a data + params list, not a fit. Compare data fields.
  expect_equal(dev_res$data, ref_res$data, tolerance = 1e-12,
               info = "init_only data object differs")
  expect_equal(dev_res$params, ref_res$params, tolerance = 1e-12,
               info = "init_only params object differs")
})

# =============================================================================
# Part 7: susie_rss_lambda parameter grid
# =============================================================================

for (lam in c(0.01, 0.1, 0.5)) {
  for (L in c(5, 10)) {
    local({
      lam_l <- lam
      L_l <- L
      test_that(sprintf("susie_rss_lambda: L=%d lambda=%g matches reference",
                        L_l, lam_l), {
        skip_if_no_mismatch_reference()
        d <- .mismatch_make_data()
        args <- list(z = d$z, R = d$R, n = d$n, L = L_l,
                     lambda = lam_l, max_iter = 50)
        compare_to_mismatch_reference("susie_rss_lambda", args)
      })
    })
  }
}

test_that("susie_rss_lambda with default max_iter (triggers hint) matches reference", {
  skip_if_no_mismatch_reference()
  d <- .mismatch_make_data()
  args <- list(z = d$z, R = d$R, n = d$n, L = 5, lambda = 0.1)
  compare_to_mismatch_reference("susie_rss_lambda", args)
})

# =============================================================================
# Part 8: Default max_iter (susie_rss hint path)
# =============================================================================
# Confirms the max_iter NULL -> 50 + hint resolution that moved from the
# entry point to the constructor produces identical fits to the reference.

test_that("susie_rss with default max_iter (hint path) matches reference", {
  skip_if_no_mismatch_reference()
  d <- .mismatch_make_data()
  args <- list(z = d$z, R = d$R, n = d$n, L = 5)
  compare_to_mismatch_reference("susie_rss", args)
})

test_that("susie_rss with default max_iter on X path matches reference", {
  skip_if_no_mismatch_reference()
  d <- .mismatch_make_data()
  args <- list(z = d$z, X = d$X, n = d$n, L = 5)
  compare_to_mismatch_reference("susie_rss", args)
})
