# =============================================================================
# HELPER FUNCTIONS FOR SUSIE-INF REFERENCE COMPARISON
# =============================================================================
#
# These functions compare the current SuSiE-inf implementation
# (`unmappable_effects = "inf"`) against the upstream gold-standard commit
# stephenslab/susieR@f110692 ("code refactor"), which is the upstream master
# at the time of the SuSiE-inf speedup work in this fork.
#
# Two paths are exercised:
#   1. susie_ss(..., unmappable_effects = "inf")    -- sufficient-statistics
#   2. susie(X, y, unmappable_effects = "inf")      -- individual-level data
#
# The ss path is algebraically identical between dev and reference (Task 2
# only rearranges intermediates), so machine-precision agreement is the
# contract.  The individual path in dev uses thin SVD of standardized X
# instead of the reference's eigen-of-XtX, so a slightly looser tolerance
# is used to absorb LAPACK summation-order differences.
#
# This helper parallels helper_mismatch_reference.R and helper_nig_reference.R.
# Prefixed names (`.inf_ref_*`, `*_inf_*`) prevent collisions when testthat
# sources all helpers at startup.
#
library(pkgload)
library(rprojroot)

# Reference package details for the SuSiE-inf comparison
.inf_ref_repo   <- "stephenslab/susieR"
.inf_ref_commit <- "f110692"

# Cached environments (separate from other helpers' globals)
.inf_ref_env         <- NULL
.inf_dev_env         <- NULL
.inf_ref_source_path <- NULL

get_inf_reference_source <- function() {
  if (!is.null(.inf_ref_source_path) && dir.exists(.inf_ref_source_path)) {
    return(.inf_ref_source_path)
  }

  ref_source <- file.path(tempdir(), "susieR_inf_reference_source")

  if (!dir.exists(ref_source)) {
    message("Downloading SuSiE-inf reference source from GitHub...")

    result <- system(sprintf("git clone -q https://github.com/%s.git %s 2>&1",
                             .inf_ref_repo, ref_source),
                     intern = FALSE)
    if (result != 0) stop("Failed to clone reference package")

    result <- system(sprintf("cd %s && git checkout -q %s 2>&1",
                             ref_source, .inf_ref_commit),
                     intern = FALSE)
    if (result != 0) stop("Failed to checkout commit ", .inf_ref_commit)

    message("Reference source downloaded")
  }

  .inf_ref_source_path <<- ref_source
  return(ref_source)
}

load_inf_reference_env <- function() {
  if (!is.null(.inf_ref_env)) return(.inf_ref_env)
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("Package 'pkgload' is required.")

  ref_source <- get_inf_reference_source()
  message("Loading SuSiE-inf reference package with pkgload...")
  env <- pkgload::load_all(ref_source, export_all = FALSE, quiet = TRUE)
  .inf_ref_env <<- env
  env
}

load_inf_development_env <- function() {
  if (!is.null(.inf_dev_env)) return(.inf_dev_env)
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("Package 'pkgload' is required.")

  dev_source <- tryCatch(
    rprojroot::find_root(rprojroot::is_r_package),
    error = function(e) normalizePath(file.path(getwd(), "../.."))
  )
  message("Loading SuSiE-inf development package with pkgload...")
  env <- pkgload::load_all(dev_source, export_all = FALSE, quiet = TRUE)
  .inf_dev_env <<- env
  env
}

skip_if_no_inf_reference <- function() {
  tryCatch({
    load_inf_reference_env()
    load_inf_development_env()
  }, error = function(e) {
    skip(paste("SuSiE-inf reference comparison not available:", e$message))
  })
}

# -----------------------------------------------------------------------
# compare_to_inf_reference
#
# Calls a function on both reference and dev environments with the same
# arguments, then asserts deep field-by-field equality including the
# inf-specific outputs (tau2, theta).
# -----------------------------------------------------------------------
compare_to_inf_reference <- function(func_name, args, tolerance = 1e-10) {
  skip_if_no_inf_reference()

  ref_env <- load_inf_reference_env()
  dev_env <- load_inf_development_env()

  ref_func <- ref_env$env[[func_name]]
  dev_func <- dev_env$env[[func_name]]
  if (is.null(ref_func)) stop("'", func_name, "' not found in reference package")
  if (is.null(dev_func)) stop("'", func_name, "' not found in development package")

  dev_result <- suppressWarnings(suppressMessages(do.call(dev_func, args)))
  ref_result <- suppressWarnings(suppressMessages(do.call(ref_func, args)))

  expect_equal_susie_inf_objects(dev_result, ref_result, tolerance)
  invisible(list(dev = dev_result, ref = ref_result))
}

# -----------------------------------------------------------------------
# expect_equal_susie_inf_objects
#
# Deep comparison of susie / susie_ss objects with unmappable_effects = "inf".
# Compares all standard fit fields (alpha, mu, mu2, V, sigma2, elbo, pip,
# niter, sets) plus the inf-specific (tau2, theta).  Skips fields that are
# NULL in both objects.
# -----------------------------------------------------------------------
expect_equal_susie_inf_objects <- function(dev_obj, ref_obj, tolerance = 1e-10) {
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

  # SuSiE-inf specific.  theta is a 1-col matrix on the ss path and a
  # vector on the individual path -- ignore the attribute, compare values.
  expect_equal(dev_obj$tau2,   ref_obj$tau2,   tolerance = tolerance,
               info = "tau2 differs")
  expect_equal(as.vector(dev_obj$theta), as.vector(ref_obj$theta),
               tolerance = tolerance, info = "theta differs")

  # ELBO & convergence (elbo is a vector; niter / converged are scalars)
  if (!is.null(dev_obj$elbo) || !is.null(ref_obj$elbo))
    expect_equal(dev_obj$elbo, ref_obj$elbo, tolerance = tolerance,
                 info = "elbo differs")
  expect_equal(dev_obj$niter,     ref_obj$niter,     info = "niter differs")
  expect_equal(dev_obj$converged, ref_obj$converged, info = "converged differs")

  # PIPs & credible sets
  expect_equal(dev_obj$pip, ref_obj$pip, tolerance = tolerance,
               info = "pip differs")
  if (!is.null(dev_obj$sets) || !is.null(ref_obj$sets)) {
    expect_equal(dev_obj$sets$cs,       ref_obj$sets$cs,
                 info = "sets$cs differs")
    expect_equal(dev_obj$sets$cs_index, ref_obj$sets$cs_index,
                 info = "sets$cs_index differs")
    expect_equal(dev_obj$sets$coverage, ref_obj$sets$coverage,
                 tolerance = tolerance, info = "sets$coverage differs")
    if (!is.null(dev_obj$sets$purity) || !is.null(ref_obj$sets$purity))
      expect_equal(dev_obj$sets$purity, ref_obj$sets$purity,
                   tolerance = tolerance, info = "sets$purity differs")
  }

  # Intercept (susie only; susie_ss objects may not carry one)
  if (!is.null(dev_obj$intercept) || !is.null(ref_obj$intercept))
    expect_equal(dev_obj$intercept, ref_obj$intercept, tolerance = tolerance,
                 info = "intercept differs")

  invisible(TRUE)
}
