# =============================================================================
# HELPER FUNCTIONS FOR SERVIN-STEPHENS REFERENCE COMPARISON
# =============================================================================
#
# These functions compare the local susieR implementation of
# estimate_residual_method = "Servin_Stephens" against the reference
# implementation on stephenslab/susieR@fix-susie-small-sigma-update
# (commit a999d44), where the equivalent feature is small = TRUE.
#
# This helper parallels helper_reference.R but targets a different
# reference commit and maps between the two parameter interfaces.
#
library(pkgload)
library(rprojroot)

# Reference package details for the Servin-Stephens comparison
.ss_ref_repo   <- "stephenslab/susieR"
.ss_ref_commit <- "a999d44"

# Cached environments (separate from helper_reference.R's globals)
.ss_ref_env         <- NULL
.ss_dev_env         <- NULL
.ss_ref_source_path <- NULL

# Get reference source for the fix-susie-small-sigma-update branch
get_ss_reference_source <- function() {
  if (!is.null(.ss_ref_source_path) && dir.exists(.ss_ref_source_path)) {
    return(.ss_ref_source_path)
  }

  ref_source <- file.path(tempdir(), "susieR_ss_reference_source")

  if (!dir.exists(ref_source)) {
    message("Downloading Servin-Stephens reference source from GitHub...")

    result <- system(sprintf("git clone -q https://github.com/%s.git %s 2>&1",
                             .ss_ref_repo, ref_source),
                     intern = FALSE)

    if (result != 0) {
      stop("Failed to clone reference package")
    }

    result <- system(sprintf("cd %s && git checkout -q %s 2>&1",
                             ref_source, .ss_ref_commit),
                     intern = FALSE)

    if (result != 0) {
      stop("Failed to checkout commit ", .ss_ref_commit)
    }

    message("\u2713 Servin-Stephens reference source downloaded")
  }

  .ss_ref_source_path <<- ref_source
  return(ref_source)
}

# Load the fix-susie-small-sigma-update reference using pkgload
load_ss_reference_env <- function() {
  if (!is.null(.ss_ref_env)) {
    return(.ss_ref_env)
  }

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required. Install with: install.packages('pkgload')")
  }

  ref_source <- get_ss_reference_source()

  message("Loading Servin-Stephens reference package with pkgload...")
  env <- pkgload::load_all(ref_source, export_all = FALSE, quiet = TRUE)

  .ss_ref_env <<- env
  return(env)
}

# Load development package using pkgload
load_ss_development_env <- function() {
  if (!is.null(.ss_dev_env)) {
    return(.ss_dev_env)
  }

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required. Install with: install.packages('pkgload')")
  }

  dev_source <- tryCatch({
    rprojroot::find_root(rprojroot::is_r_package)
  }, error = function(e) {
    normalizePath(file.path(getwd(), "../.."))
  })

  message("Loading development package with pkgload...")
  env <- pkgload::load_all(dev_source, export_all = FALSE, quiet = TRUE)

  .ss_dev_env <<- env
  return(env)
}

# Skip test if reference not available
skip_if_no_ss_reference <- function() {
  tryCatch({
    load_ss_reference_env()
    load_ss_development_env()
  }, error = function(e) {
    skip(paste("Servin-Stephens reference comparison not available:", e$message))
  })
}

# -----------------------------------------------------------------------
# compare_servin_stephens_to_reference
#
# Runs susie() with estimate_residual_method = "Servin_Stephens" on the
# development package and susie() with small = TRUE on the reference
# branch, then compares all output fields.
#
# Parameters:
#   dev_args  - named list of arguments for the development susie() call
#               (must include X and y; estimate_residual_method is set
#               automatically to "Servin_Stephens")
#   ref_args  - (optional) named list of arguments for the reference
#               susie() call. If NULL, derived from dev_args by mapping
#               estimate_residual_method -> small = TRUE and
#               tol -> tol_small.
#   tolerance - numeric tolerance for expect_equal comparisons
# -----------------------------------------------------------------------
compare_servin_stephens_to_reference <- function(dev_args,
                                                  ref_args = NULL,
                                                  tolerance = 1e-5) {
  skip_if_no_ss_reference()

  ref_env <- load_ss_reference_env()
  dev_env <- load_ss_development_env()

  # Ensure the dev call uses Servin_Stephens
  dev_args$estimate_residual_method <- "Servin_Stephens"

  # Match reference behavior: disable V null threshold check and use

  # the same convergence tolerance as the reference (tol_small = 1e-4)
  if (is.null(dev_args$check_null_threshold))
    dev_args$check_null_threshold <- -Inf
  if (is.null(dev_args$tol))
    dev_args$tol <- 1e-4

  # Build reference args by mapping interface differences
  if (is.null(ref_args)) {
    ref_args <- dev_args

    # Map estimate_residual_method -> small
    ref_args$estimate_residual_method <- NULL
    ref_args$small <- TRUE

    # Map tol -> tol_small (reference replaces tol with tol_small when small=TRUE)
    if (!is.null(ref_args$tol)) {
      ref_args$tol_small <- ref_args$tol
      ref_args$tol       <- NULL
    }

    # Remove parameters that don't exist in the reference interface
    ref_args$convergence_method <- NULL

    # The reference uses s_init instead of model_init
    if (!is.null(ref_args$model_init)) {
      ref_args$s_init     <- ref_args$model_init
      ref_args$model_init <- NULL
    }
  }

  ref_func <- ref_env$env[["susie"]]
  dev_func <- dev_env$env[["susie"]]

  if (is.null(ref_func)) stop("susie() not found in reference package")
  if (is.null(dev_func)) stop("susie() not found in development package")

  # Suppress known warnings (method override messages)
  dev_result <- suppressWarnings(do.call(dev_func, dev_args))
  ref_result <- suppressWarnings(do.call(ref_func, ref_args))

  # Return both results for custom assertions
  invisible(list(dev = dev_result, ref = ref_result))
}

# -----------------------------------------------------------------------
# expect_equal_servin_stephens_objects
#
# Deep comparison of susie objects produced under the Servin-Stephens /
# small = TRUE prior. Compares the standard fields (alpha, mu, mu2, V,
# sigma2, elbo, fitted, intercept, pip, sets) plus the NIG-specific
# rv field.
# -----------------------------------------------------------------------
expect_equal_servin_stephens_objects <- function(dev_obj, ref_obj,
                                                 tolerance = 1e-5) {

  # --- Core posterior quantities ---
  expect_equal(dev_obj$alpha, ref_obj$alpha, tolerance = tolerance,
               info = "alpha (posterior inclusion probabilities) differ")

  expect_equal(dev_obj$mu, ref_obj$mu, tolerance = tolerance,
               info = "mu (posterior means) differ")

  expect_equal(dev_obj$mu2, ref_obj$mu2, tolerance = tolerance,
               info = "mu2 (posterior second moments) differ")

  # --- Variance parameters ---
  expect_equal(dev_obj$V, ref_obj$V, tolerance = tolerance,
               info = "V (prior variance, after rv scaling) differs")

  expect_equal(dev_obj$sigma2, ref_obj$sigma2, tolerance = tolerance,
               info = "sigma2 (residual variance) differs")

  # --- Residual variance per effect (NIG-specific) ---
  if (!is.null(dev_obj$rv) && !is.null(ref_obj$rv)) {
    expect_equal(dev_obj$rv, ref_obj$rv, tolerance = tolerance,
                 info = "rv (per-effect residual variance) differs")
  }

  # --- ELBO / convergence ---
  # For L = 1 the dev package intentionally uses ELBO convergence while
  # the reference uses PIP convergence, so niter and elbo may differ.
  # Only compare these for L > 1.
  L <- nrow(dev_obj$alpha)
  if (L > 1) {
    expect_equal(dev_obj$niter, ref_obj$niter,
                 info = "Number of iterations differs")

    expect_equal(dev_obj$converged, ref_obj$converged,
                 info = "Convergence status differs")
  }

  # For L = 1 the ELBO (loglik) is well-defined; compare if both present
  # and the iteration counts match (they may differ due to convergence method)
  if (!is.null(dev_obj$elbo) && !is.null(ref_obj$elbo) &&
      !all(is.na(dev_obj$elbo)) && !all(is.na(ref_obj$elbo)) &&
      length(dev_obj$elbo) == length(ref_obj$elbo)) {
    expect_equal(dev_obj$elbo, ref_obj$elbo, tolerance = tolerance,
                 info = "ELBO values differ")
  }

  # --- Fitted values and intercept ---
  if (!is.null(dev_obj$fitted) && !is.null(ref_obj$fitted)) {
    expect_equal(dev_obj$fitted, ref_obj$fitted, tolerance = tolerance,
                 info = "Fitted values differ")
  }

  expect_equal(dev_obj$intercept, ref_obj$intercept, tolerance = tolerance,
               info = "Intercept differs")

  # --- PIPs ---
  if (!is.null(dev_obj$pip) && !is.null(ref_obj$pip)) {
    expect_equal(dev_obj$pip, ref_obj$pip, tolerance = tolerance,
                 info = "PIPs differ")
  }

  # --- Credible sets ---
  if (!is.null(dev_obj$sets) && !is.null(ref_obj$sets)) {
    expect_equal(dev_obj$sets$cs, ref_obj$sets$cs,
                 info = "Credible sets differ")
    if (!is.null(dev_obj$sets$purity) && !is.null(ref_obj$sets$purity)) {
      expect_equal(dev_obj$sets$purity, ref_obj$sets$purity,
                   tolerance = tolerance, info = "CS purity differs")
    }
    expect_equal(dev_obj$sets$coverage, ref_obj$sets$coverage,
                 tolerance = tolerance, info = "CS coverage differs")
  }

  invisible(TRUE)
}

# -----------------------------------------------------------------------
# run_ss_and_individual_servin_stephens
#
# Given X, y, and extra arguments (L, standardize, intercept, alpha0,
# beta0, etc.), runs both susie() and susie_ss() with
# estimate_residual_method = "Servin_Stephens", ensuring that the
# sufficient statistics are computed to match susie()'s internal
# preprocessing.
#
# Returns list(ind = ..., ss = ...) with both results.
# -----------------------------------------------------------------------
run_ss_and_individual_servin_stephens <- function(X, y, extra_args = list()) {
  n <- nrow(X)
  p <- ncol(X)

  # Extract preprocessing settings (defaults match susie)
  intercept   <- if (!is.null(extra_args$intercept))   extra_args$intercept   else TRUE
  standardize <- if (!is.null(extra_args$standardize)) extra_args$standardize else TRUE

  # Preprocess exactly as susie() does internally
  y_mean     <- mean(y)
  X_colmeans <- colMeans(X)

  if (intercept) {
    y_c <- y - y_mean
    X_c <- scale(X, center = TRUE, scale = FALSE)
  } else {
    y_c <- y
    X_c <- X
  }

  if (standardize) {
    csd <- apply(X, 2, sd)
    csd[csd == 0] <- 1
    X_cs <- t(t(X_c) / csd)
  } else {
    X_cs <- X_c
  }

  # Compute sufficient statistics from preprocessed data
  XtX <- crossprod(X_cs)
  Xty <- drop(t(X_cs) %*% y_c)
  yty <- sum(y_c^2)

  # Run individual-level susie
  ind_args <- c(list(X = X, y = y,
                     estimate_residual_method = "Servin_Stephens"),
                extra_args)
  res_ind <- suppressWarnings(do.call(susie, ind_args))

  # Build SS arguments: remove individual-only params, add SS-specific ones
  ss_extra <- extra_args
  ss_extra$intercept   <- NULL
  ss_extra$standardize <- NULL

  ss_args <- c(list(XtX = XtX, Xty = Xty, yty = yty, n = n,
                    estimate_residual_method = "Servin_Stephens",
                    standardize = FALSE,
                    X_colmeans = if (intercept) X_colmeans else NA,
                    y_mean     = if (intercept) y_mean else NA),
               ss_extra)
  res_ss <- suppressWarnings(do.call(susie_ss, ss_args))

  list(ind = res_ind, ss = res_ss)
}

# -----------------------------------------------------------------------
# run_rss_and_individual_servin_stephens
#
# Given X, y, and extra arguments (L, standardize, intercept, alpha0,
# beta0, etc.), runs both susie() and susie_rss() with
# estimate_residual_method = "Servin_Stephens", ensuring that the
# summary statistics (bhat, shat, R, var_y) are computed to match
# susie()'s internal preprocessing.
#
# Uses the bhat/shat/var_y input path of susie_rss(), which recovers
# exact sufficient statistics (XtX, Xty, yty) from summary statistics.
# This is necessary because the NIG prior's alpha0/beta0 break scale
# invariance, so the z-score-only path would not match.
#
# Returns list(ind = ..., rss = ...) with both results.
# -----------------------------------------------------------------------
run_rss_and_individual_servin_stephens <- function(X, y, extra_args = list()) {
  n <- nrow(X)
  p <- ncol(X)

  # Extract preprocessing settings (defaults match susie)
  intercept   <- if (!is.null(extra_args$intercept))   extra_args$intercept   else TRUE
  standardize <- if (!is.null(extra_args$standardize)) extra_args$standardize else TRUE

  # Preprocess exactly as susie() does internally
  if (intercept) {
    y_c <- y - mean(y)
    X_c <- scale(X, center = TRUE, scale = FALSE)
  } else {
    y_c <- y
    X_c <- X
  }

  if (standardize) {
    csd <- apply(X, 2, sd)
    csd[csd == 0] <- 1
    X_cs <- t(t(X_c) / csd)
  } else {
    X_cs <- X_c
  }

  # Compute correlation matrix from preprocessed data
  R <- cor(X_cs)
  R <- (R + t(R)) / 2  # ensure symmetry

  # Compute bhat/shat via univariate regression on preprocessed data
  # center=FALSE because we already centered
  ss <- univariate_regression(X_cs, y_c, center = FALSE)

  # Compute var_y
  var_y <- sum(y_c^2) / (n - 1)

  # Run individual-level susie
  ind_args <- c(list(X = X, y = y,
                     estimate_residual_method = "Servin_Stephens"),
                extra_args)
  res_ind <- suppressWarnings(do.call(susie, ind_args))

  # Build RSS arguments: remove individual-only params, add RSS-specific ones
  rss_extra <- extra_args
  rss_extra$intercept   <- NULL
  rss_extra$standardize <- NULL

  rss_args <- c(list(bhat = ss$betahat, shat = ss$sebetahat,
                     R = R, n = n, var_y = var_y,
                     estimate_residual_method = "Servin_Stephens",
                     standardize = FALSE),
                rss_extra)
  res_rss <- suppressWarnings(do.call(susie_rss, rss_args))

  list(ind = res_ind, rss = res_rss)
}

# -----------------------------------------------------------------------
# expect_rss_matches_individual_ss
#
# Deep comparison of susie objects produced by the individual-level and
# RSS interfaces under Servin-Stephens. Delegates to
# expect_ss_matches_individual_ss by mapping rss -> ss.
# -----------------------------------------------------------------------
expect_rss_matches_individual_ss <- function(res, tolerance = 1e-6) {
  expect_ss_matches_individual_ss(
    list(ind = res$ind, ss = res$rss),
    tolerance = tolerance
  )
}

# -----------------------------------------------------------------------
# expect_ss_matches_individual_ss
#
# Deep comparison of susie objects produced by the individual-level and
# sufficient-statistics interfaces under Servin-Stephens.
# -----------------------------------------------------------------------
expect_ss_matches_individual_ss <- function(res, tolerance = 1e-6) {
  ind <- res$ind
  ss  <- res$ss

  # Core posterior quantities
  expect_equal(ind$alpha, ss$alpha, tolerance = tolerance,
               info = "alpha differs between susie and susie_ss")
  expect_equal(ind$mu, ss$mu, tolerance = tolerance,
               info = "mu differs between susie and susie_ss")
  expect_equal(ind$mu2, ss$mu2, tolerance = tolerance,
               info = "mu2 differs between susie and susie_ss")
  expect_equal(ind$V, ss$V, tolerance = tolerance,
               info = "V differs between susie and susie_ss")
  expect_equal(ind$sigma2, ss$sigma2, tolerance = tolerance,
               info = "sigma2 differs between susie and susie_ss")
  expect_equal(ind$pip, ss$pip, tolerance = tolerance,
               info = "pip differs between susie and susie_ss")

  # Convergence
  expect_equal(ind$niter, ss$niter,
               info = "niter differs between susie and susie_ss")
  expect_equal(ind$converged, ss$converged,
               info = "converged differs between susie and susie_ss")

  # NIG-specific: per-effect residual variance
  if (!is.null(ind$rv) && !is.null(ss$rv)) {
    expect_equal(ind$rv, ss$rv, tolerance = tolerance,
                 info = "rv differs between susie and susie_ss")
  }

  # Credible sets
  if (!is.null(ind$sets) && !is.null(ss$sets)) {
    expect_equal(ind$sets$cs, ss$sets$cs,
                 info = "Credible sets differ between susie and susie_ss")
    if (!is.null(ind$sets$purity) && !is.null(ss$sets$purity)) {
      expect_equal(ind$sets$purity, ss$sets$purity,
                   tolerance = tolerance,
                   info = "CS purity differs between susie and susie_ss")
    }
  }

  invisible(TRUE)
}
