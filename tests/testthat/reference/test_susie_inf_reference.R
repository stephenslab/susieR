# Source helper functions
source(file.path("..", "helper_inf_reference.R"), local = TRUE)

context("susie / susie_ss SuSiE-inf reference comparison")

# =============================================================================
# REFERENCE TESTS FOR SUSIE-INF (unmappable_effects = "inf")
# =============================================================================
#
# These tests compare the current SuSiE-inf implementation against the
# upstream reference commit stephenslab/susieR@f110692.  Two paths are
# exercised:
#
#   1. susie_ss(XtX, Xty, yty, n, unmappable_effects = "inf")
#   2. susie(X, y, unmappable_effects = "inf")
#
# Tolerance depends on whether the prior-variance update is closed-form
# (`estimate_prior_method = "EM"` or `"simple"`) or iterative
# (`"optim"`):
#
#   - `optim` calls L-BFGS-B inside each SER step.  Its
#     gradient-norm-gated stopping criterion is sensitive to 1e-16 input
#     noise: dev and reference can land at slightly different V values,
#     growing to ~1e-8 over the converged fit.  This is inherent to
#     iterative VI, not an algorithmic divergence.
#
#   - `EM` / `simple` use closed-form updates with no convergence check
#     inside the update step.  Dev and reference agree at machine
#     precision (~1e-14 to ~1e-15) across all iterations.
#
# MoM residual-variance estimation is fully closed-form (`solve(A, x)`),
# so MoM + EM / MoM + simple are end-to-end deterministic and pin
# machine-precision identity.  MLE uses optim on (sigma^2, tau^2) and
# carries the same iterative-noise floor as `optim` prior updates.
#
# On the individual path dev replaces the reference's eigen(XtX) detour
# with a thin SVD of standardized X.  The eigenspace math is
# mathematically identical (null-space eigencomponents contribute 0 to
# every sum) but LAPACK eigen and SVD use different summation orders, so
# the noise floor is one order of magnitude higher (~1e-13 vs ~1e-15).
#
# Fixtures are p = 1000, n = 500 with a sparse signal and optional dense
# polygenic background.  L = 5 unless noted.

# -----------------------------------------------------------------------------
# Fixture builders
# -----------------------------------------------------------------------------

make_inf_X_y <- function(n = 500, p = 1000, k_dense = 0, dense_sd = 0,
                         seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- rep(0, p)
  beta_true[sample(p, 3)] <- c(1.5, -1.2, 0.8)
  if (k_dense > 0) {
    idx <- sample(p, k_dense)
    beta_true[idx] <- beta_true[idx] + rnorm(k_dense, sd = dense_sd)
  }
  y <- as.vector(X %*% beta_true) + rnorm(n, sd = 0.5)
  list(X = X, y = y)
}

make_inf_ss <- function(X, y) {
  # Match the canonical centering/scaling that compute_XtX / compute_Xty
  # apply on-the-fly so the sufficient statistics correspond to the same
  # model the individual path fits.
  X_cm  <- colMeans(X)
  X_csd <- apply(X, 2, sd)
  X_std <- sweep(sweep(X, 2, X_cm, "-"), 2, X_csd, "/")
  y_c   <- y - mean(y)
  list(
    XtX = crossprod(X_std),
    Xty = as.vector(crossprod(X_std, y_c)),
    yty = sum(y_c^2),
    n   = nrow(X)
  )
}

# Tolerances:
#   - bit_ident_tol:       closed-form path (MoM + EM / MoM + simple),
#                          ss data: machine precision.
#   - bit_ident_indiv_tol: closed-form path, individual data.  SVD vs
#                          eigen-of-XtX routing adds ~2 orders of magnitude
#                          to ELBO-trajectory noise (per-fit-field values
#                          still match at ~1e-14).
#   - ss_mom_optim_tol:    optim path (L-BFGS-B) on ss + MoM, converged
#                          fit.  Tighter than ss_tol because MoM gives the
#                          optimizer well-conditioned moment-based gradients
#                          and the prior-variance optimizer drifts very
#                          little from the reference over 30 iters.
#   - ss_tol / indiv_tol:  optim path (L-BFGS-B), converged fit (MLE on
#                          ss, MoM/MLE on individual data).  Loose tolerance
#                          absorbs gradient-norm-gating noise.
.bit_ident_tol       <- 1e-12
.bit_ident_indiv_tol <- 1e-10
.ss_mom_optim_tol    <- 1e-9
.ss_tol              <- 1e-5
.indiv_tol           <- 1e-5

# =============================================================================
# Part 0a: 1-iter bit-identity contracts (optim path)
# =============================================================================
#
# After exactly one outer iteration, the L-BFGS-B noise hasn't yet
# amplified, so even the optim prior-method path is bit-identical to
# reference.  Pinning this confirms that the per-iter algebra (including
# Task 2's caching) has not drifted.

test_that("susie_ss + inf is bit-identical to reference at 1 iter (MoM, optim)", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 11)
  ss <- make_inf_ss(d$X, d$y)
  args <- c(ss, list(L = 5,
                     unmappable_effects        = "inf",
                     estimate_residual_method  = "MoM",
                     max_iter = 1))
  compare_to_inf_reference("susie_ss", args, tolerance = .bit_ident_tol)
})

test_that("susie_ss + inf is bit-identical to reference at 1 iter (MLE, optim)", {
  skip("SuSiE-inf + MLE not supported")
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 21)
  ss <- make_inf_ss(d$X, d$y)
  args <- c(ss, list(L = 5,
                     unmappable_effects        = "inf",
                     estimate_residual_method  = "MLE",
                     max_iter = 1))
  compare_to_inf_reference("susie_ss", args, tolerance = .bit_ident_tol)
})

# =============================================================================
# Part 0b: Machine-precision contracts for closed-form prior methods
# =============================================================================
#
# MoM + EM and MoM + simple are end-to-end deterministic (no L-BFGS-B in
# either the residual-variance or prior-variance update).  Dev and
# reference must agree to machine precision across all outer iterations.
# These tests run to convergence, not just 1 iter.

test_that("susie_ss + inf + MoM + EM is machine-precision identical to reference", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 11)
  ss <- make_inf_ss(d$X, d$y)
  args <- c(ss, list(L = 5,
                     unmappable_effects        = "inf",
                     estimate_residual_method  = "MoM",
                     estimate_prior_method     = "EM",
                     max_iter = 30, tol = 1e-4))
  compare_to_inf_reference("susie_ss", args, tolerance = .bit_ident_tol)
})

test_that("susie_ss + inf + MoM + EM (polygenic, L=10) is machine-precision identical", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 13, k_dense = 200, dense_sd = 0.20)
  ss <- make_inf_ss(d$X, d$y)
  args <- c(ss, list(L = 10,
                     unmappable_effects        = "inf",
                     estimate_residual_method  = "MoM",
                     estimate_prior_method     = "EM",
                     max_iter = 30, tol = 1e-4))
  compare_to_inf_reference("susie_ss", args, tolerance = .bit_ident_tol)
})

test_that("susie_ss + inf + MoM + simple is machine-precision identical to reference", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 11)
  ss <- make_inf_ss(d$X, d$y)
  args <- c(ss, list(L = 5,
                     unmappable_effects        = "inf",
                     estimate_residual_method  = "MoM",
                     estimate_prior_method     = "simple",
                     max_iter = 30, tol = 1e-4))
  compare_to_inf_reference("susie_ss", args, tolerance = .bit_ident_tol)
})

test_that("susie + inf + MoM + EM is machine-precision identical (individual path)", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 31)
  args <- list(X = d$X, y = d$y, L = 5,
               unmappable_effects        = "inf",
               estimate_residual_method  = "MoM",
               estimate_prior_method     = "EM",
               max_iter = 30, tol = 1e-4)
  compare_to_inf_reference("susie", args, tolerance = .bit_ident_indiv_tol)
})

test_that("susie + inf + MoM + EM (polygenic, L=10) is machine-precision identical (individual path)", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 33, k_dense = 200, dense_sd = 0.20)
  args <- list(X = d$X, y = d$y, L = 10,
               unmappable_effects        = "inf",
               estimate_residual_method  = "MoM",
               estimate_prior_method     = "EM",
               max_iter = 30, tol = 1e-4)
  compare_to_inf_reference("susie", args, tolerance = .bit_ident_indiv_tol)
})

test_that("susie + inf + MoM + simple is machine-precision identical (individual path)", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 31)
  args <- list(X = d$X, y = d$y, L = 5,
               unmappable_effects        = "inf",
               estimate_residual_method  = "MoM",
               estimate_prior_method     = "simple",
               max_iter = 30, tol = 1e-4)
  compare_to_inf_reference("susie", args, tolerance = .bit_ident_indiv_tol)
})

# =============================================================================
# Part 1: susie_ss(..., unmappable_effects = "inf") - MoM
# =============================================================================

test_that("susie_ss + inf + MoM matches reference (sparse only, L=5)", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 11)
  ss <- make_inf_ss(d$X, d$y)
  args <- c(ss, list(L = 5,
                     unmappable_effects        = "inf",
                     estimate_residual_method  = "MoM",
                     max_iter = 30, tol = 1e-4))
  compare_to_inf_reference("susie_ss", args, tolerance = .ss_mom_optim_tol)
})

test_that("susie_ss + inf + MoM matches reference (sparse + polygenic, L=5)", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 12, k_dense = 200, dense_sd = 0.10)
  ss <- make_inf_ss(d$X, d$y)
  args <- c(ss, list(L = 5,
                     unmappable_effects        = "inf",
                     estimate_residual_method  = "MoM",
                     max_iter = 30, tol = 1e-4))
  compare_to_inf_reference("susie_ss", args, tolerance = .ss_mom_optim_tol)
})

test_that("susie_ss + inf + MoM matches reference (sparse + polygenic, L=10)", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 13, k_dense = 200, dense_sd = 0.20)
  ss <- make_inf_ss(d$X, d$y)
  args <- c(ss, list(L = 10,
                     unmappable_effects        = "inf",
                     estimate_residual_method  = "MoM",
                     max_iter = 30, tol = 1e-4))
  compare_to_inf_reference("susie_ss", args, tolerance = .ss_mom_optim_tol)
})

# =============================================================================
# Part 2: susie_ss(..., unmappable_effects = "inf") - MLE
# =============================================================================

test_that("susie_ss + inf + MLE matches reference (sparse only, L=5)", {
  skip("SuSiE-inf + MLE not supported")
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 21)
  ss <- make_inf_ss(d$X, d$y)
  args <- c(ss, list(L = 5,
                     unmappable_effects        = "inf",
                     estimate_residual_method  = "MLE",
                     max_iter = 30, tol = 1e-4))
  compare_to_inf_reference("susie_ss", args, tolerance = .ss_tol)
})

test_that("susie_ss + inf + MLE matches reference (sparse + polygenic, L=5)", {
  skip("SuSiE-inf + MLE not supported")
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 22, k_dense = 200, dense_sd = 0.10)
  ss <- make_inf_ss(d$X, d$y)
  args <- c(ss, list(L = 5,
                     unmappable_effects        = "inf",
                     estimate_residual_method  = "MLE",
                     max_iter = 30, tol = 1e-4))
  compare_to_inf_reference("susie_ss", args, tolerance = .ss_tol)
})

# =============================================================================
# Part 3: susie(X, y, unmappable_effects = "inf") - individual-data MoM
# =============================================================================

test_that("susie + inf + MoM matches reference (sparse only, L=5)", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 31)
  args <- list(X = d$X, y = d$y, L = 5,
               unmappable_effects        = "inf",
               estimate_residual_method  = "MoM",
               max_iter = 30, tol = 1e-4)
  compare_to_inf_reference("susie", args, tolerance = .indiv_tol)
})

test_that("susie + inf + MoM matches reference (sparse + polygenic, L=5)", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 32, k_dense = 200, dense_sd = 0.10)
  args <- list(X = d$X, y = d$y, L = 5,
               unmappable_effects        = "inf",
               estimate_residual_method  = "MoM",
               max_iter = 30, tol = 1e-4)
  compare_to_inf_reference("susie", args, tolerance = .indiv_tol)
})

test_that("susie + inf + MoM matches reference (sparse + polygenic, L=10)", {
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 33, k_dense = 200, dense_sd = 0.20)
  args <- list(X = d$X, y = d$y, L = 10,
               unmappable_effects        = "inf",
               estimate_residual_method  = "MoM",
               max_iter = 30, tol = 1e-4)
  compare_to_inf_reference("susie", args, tolerance = .indiv_tol)
})

# =============================================================================
# Part 4: susie(X, y, unmappable_effects = "inf") - individual-data MLE
# =============================================================================

test_that("susie + inf + MLE matches reference (sparse only, L=5)", {
  skip("SuSiE-inf + MLE not supported")
  skip_if_no_inf_reference()
  d <- make_inf_X_y(seed = 41)
  args <- list(X = d$X, y = d$y, L = 5,
               unmappable_effects        = "inf",
               estimate_residual_method  = "MLE",
               max_iter = 30, tol = 1e-4)
  compare_to_inf_reference("susie", args, tolerance = .indiv_tol)
})

# MLE on a heavy polygenic individual-path fixture is unstable: optim's
# L-BFGS-B on the joint (sigma^2, tau^2) objective fails to converge
# identically in dev vs reference because the SVD vs eigen-of-XtX paths
# perturb the gradient by ~1e-12, which is enough to land on different
# non-converged points.  Skipped (matches the existing MLE convergence
# warning seen in scratch_check_tau2.R).
