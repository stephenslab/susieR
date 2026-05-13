# =============================================================================
# HELPER FUNCTIONS FOR RSS MISMATCH REFERENCE COMPARISON
# =============================================================================
#
# These functions compare the current susie_rss / susie_rss_lambda
# implementation against the upstream gold-standard commit
# stephenslab/susieR@f110692 ("code refactor"), which is the head of
# the upstream master branch at the time of the susie_rss thin-wrapper
# refactor. The current public "eb" mode has intentionally changed since this
# reference commit, while the previous public behavior is still available
# internally as "eb_ser_init" for targeted comparisons.
#
# The contract: the refactored thin-wrapper susie_rss produces
# machine-precision identical output to the upstream gold-standard
# implementation across the full RSS-mismatch surface (R_mismatch,
# R_finite, multi-panel, susie_rss_lambda).
#
# This helper parallels helper_reference.R and helper_nig_reference.R
# but targets a different reference commit. It uses prefixed names
# (`.mismatch_ref_*`, `*_mismatch_*`) so it does not collide with the
# other reference helpers when testthat sources them all at startup.
#
library(pkgload)
library(rprojroot)

# Reference package details for the RSS-mismatch comparison
.mismatch_ref_repo   <- "stephenslab/susieR"
.mismatch_ref_commit <- "f110692"

# Cached environments (separate from other helpers' globals)
.mismatch_ref_env         <- NULL
.mismatch_dev_env         <- NULL
.mismatch_ref_source_path <- NULL

# Get reference source for the pinned commit
get_mismatch_reference_source <- function() {
  if (!is.null(.mismatch_ref_source_path) && dir.exists(.mismatch_ref_source_path)) {
    return(.mismatch_ref_source_path)
  }

  ref_source <- file.path(tempdir(), "susieR_mismatch_reference_source")

  if (!dir.exists(ref_source)) {
    message("Downloading RSS-mismatch reference source from GitHub...")

    result <- system(sprintf("git clone -q https://github.com/%s.git %s 2>&1",
                             .mismatch_ref_repo, ref_source),
                     intern = FALSE)
    if (result != 0) stop("Failed to clone reference package")

    result <- system(sprintf("cd %s && git checkout -q %s 2>&1",
                             ref_source, .mismatch_ref_commit),
                     intern = FALSE)
    if (result != 0) stop("Failed to checkout commit ", .mismatch_ref_commit)

    message("✓ RSS-mismatch reference source downloaded")
  }

  .mismatch_ref_source_path <<- ref_source
  return(ref_source)
}

load_mismatch_reference_env <- function() {
  if (!is.null(.mismatch_ref_env)) return(.mismatch_ref_env)
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("Package 'pkgload' is required.")

  ref_source <- get_mismatch_reference_source()
  message("Loading RSS-mismatch reference package with pkgload...")
  env <- pkgload::load_all(ref_source, export_all = FALSE, quiet = TRUE)
  .mismatch_ref_env <<- env
  env
}

load_mismatch_development_env <- function() {
  if (!is.null(.mismatch_dev_env)) return(.mismatch_dev_env)
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("Package 'pkgload' is required.")

  dev_source <- tryCatch(
    rprojroot::find_root(rprojroot::is_r_package),
    error = function(e) normalizePath(file.path(getwd(), "../.."))
  )
  message("Loading RSS-mismatch development package with pkgload...")
  env <- pkgload::load_all(dev_source, export_all = FALSE, quiet = TRUE)
  .mismatch_dev_env <<- env
  env
}

skip_if_no_mismatch_reference <- function() {
  tryCatch({
    load_mismatch_reference_env()
    load_mismatch_development_env()
  }, error = function(e) {
    skip(paste("RSS-mismatch reference comparison not available:", e$message))
  })
}

# -----------------------------------------------------------------------
# compare_to_mismatch_reference
#
# Calls a function (default susie_rss) on both the reference and dev
# environments with the same arguments, then asserts deep field-by-field
# equality at machine precision. Because the reference is the
# pre-refactor state of the same codebase, results are expected to be
# bit-identical; the default tolerance is tight (1e-12).
# -----------------------------------------------------------------------
compare_to_mismatch_reference <- function(func_name = "susie_rss",
                                          args,
                                          tolerance = 1e-12) {
  skip_if_no_mismatch_reference()

  ref_env <- load_mismatch_reference_env()
  dev_env <- load_mismatch_development_env()

  ref_func <- ref_env$env[[func_name]]
  dev_func <- dev_env$env[[func_name]]
  if (is.null(ref_func)) stop("'", func_name, "' not found in reference package")
  if (is.null(dev_func)) stop("'", func_name, "' not found in development package")

  ref_args <- args
  if (identical(ref_args$R_mismatch, "eb"))
    skip("Current R_mismatch = 'eb' has intentionally changed since the pinned reference commit.")
  else if (identical(ref_args$R_mismatch, "eb_ser_init"))
    ref_args$R_mismatch <- "eb"

  dev_result <- suppressWarnings(suppressMessages(do.call(dev_func, args)))
  ref_result <- suppressWarnings(suppressMessages(do.call(ref_func, ref_args)))

  expect_equal_susie_rss_mismatch_objects(dev_result, ref_result, tolerance)
  invisible(list(dev = dev_result, ref = ref_result))
}

# -----------------------------------------------------------------------
# expect_equal_susie_rss_mismatch_objects
#
# Deep comparison of susie_rss objects including all R-mismatch and
# R-finite specific fields. Skips fields that are NULL in both objects
# (e.g. R_finite_diagnostics when R_finite is NULL).
# -----------------------------------------------------------------------
expect_equal_susie_rss_mismatch_objects <- function(dev_obj, ref_obj,
                                                    tolerance = 1e-12) {
  # Core posterior
  expect_equal(dev_obj$alpha,  ref_obj$alpha,  tolerance = tolerance,
               info = "alpha differs")
  expect_equal(dev_obj$mu,     ref_obj$mu,     tolerance = tolerance,
               info = "mu differs")
  expect_equal(dev_obj$mu2,    ref_obj$mu2,    tolerance = tolerance,
               info = "mu2 differs")
  expect_equal(dev_obj$V,      ref_obj$V,      tolerance = tolerance,
               info = "V differs")
  expect_equal(dev_obj$sigma2, ref_obj$sigma2, tolerance = tolerance,
               info = "sigma2 differs")

  # ELBO & convergence
  expect_equal(dev_obj$elbo,      ref_obj$elbo,      tolerance = tolerance,
               info = "elbo differs")
  expect_equal(dev_obj$niter,     ref_obj$niter,     info = "niter differs")
  expect_equal(dev_obj$converged, ref_obj$converged, info = "converged differs")
  expect_equal(dev_obj$convergence_reason, ref_obj$convergence_reason,
               info = "convergence_reason differs")

  # PIPs and credible sets
  expect_equal(dev_obj$pip, ref_obj$pip, tolerance = tolerance,
               info = "pip differs")
  if (!is.null(dev_obj$sets) || !is.null(ref_obj$sets)) {
    expect_equal(dev_obj$sets$cs,        ref_obj$sets$cs,
                 info = "sets$cs differs")
    expect_equal(dev_obj$sets$cs_index,  ref_obj$sets$cs_index,
                 info = "sets$cs_index differs")
    expect_equal(dev_obj$sets$coverage,  ref_obj$sets$coverage,
                 tolerance = tolerance, info = "sets$coverage differs")
    if (!is.null(dev_obj$sets$purity) || !is.null(ref_obj$sets$purity))
      expect_equal(dev_obj$sets$purity, ref_obj$sets$purity,
                   tolerance = tolerance, info = "sets$purity differs")
  }

  # RSS-specific
  if (!is.null(dev_obj$XtXr) || !is.null(ref_obj$XtXr))
    expect_equal(dev_obj$XtXr, ref_obj$XtXr, tolerance = tolerance,
                 info = "XtXr differs")
  expect_equal(dev_obj$intercept, ref_obj$intercept, tolerance = tolerance,
               info = "intercept differs")

  # R-mismatch / R-finite specific
  for (fld in c("lambda_bias", "B_corrected", "R_finite_B", "R_mismatch_method",
                "Q_art", "artifact_evaluable", "artifact_flag",
                "low_eigen_count", "low_eigen_fraction", "eig_delta",
                "mode_label", "R_mismatch_init")) {
    if (!is.null(dev_obj[[fld]]) || !is.null(ref_obj[[fld]])) {
      expect_equal(dev_obj[[fld]], ref_obj[[fld]],
                   tolerance = tolerance,
                   info = sprintf("%s differs", fld))
    }
  }

  # R_finite_diagnostics is a nested structure when present
  if (!is.null(dev_obj$R_finite_diagnostics) ||
      !is.null(ref_obj$R_finite_diagnostics)) {
    expect_equal(dev_obj$R_finite_diagnostics,
                 ref_obj$R_finite_diagnostics,
                 tolerance = tolerance,
                 info = "R_finite_diagnostics differs")
  }

  # Multi-panel specific
  for (fld in c("omega_weights", "single_panel_fits")) {
    if (!is.null(dev_obj[[fld]]) || !is.null(ref_obj[[fld]])) {
      expect_equal(dev_obj[[fld]], ref_obj[[fld]],
                   tolerance = tolerance,
                   info = sprintf("%s differs", fld))
    }
  }

  invisible(TRUE)
}
