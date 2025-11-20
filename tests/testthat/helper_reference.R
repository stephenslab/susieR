# =============================================================================
# HELPER FUNCTIONS FOR REFERENCE PACKAGE COMPARISON (PKGLOAD APPROACH)
# =============================================================================
#
# These functions compare the new susieR implementation against the reference
# package (stephenslab/susieR@1f9166c) to ensure results are identical.
#
# Strategy: Use pkgload to load both packages into separate environments
library(pkgload)
library(rprojroot)

# Reference package details
.ref_repo <- "stephenslab/susieR"
.ref_commit <- "1f9166c"

# Cached environments
.ref_env <- NULL
.dev_env <- NULL
.ref_source_path <- NULL

# Get reference package source (download once, cache path)
get_reference_source <- function() {
  if (!is.null(.ref_source_path) && dir.exists(.ref_source_path)) {
    return(.ref_source_path)
  }

  # Download to temp directory
  ref_source <- file.path(tempdir(), "susieR_reference_source")

  if (!dir.exists(ref_source)) {
    message("Downloading reference package source from GitHub...")

    result <- system(sprintf("git clone -q https://github.com/%s.git %s 2>&1",
                             .ref_repo, ref_source),
                     intern = FALSE)

    if (result != 0) {
      stop("Failed to clone reference package")
    }

    result <- system(sprintf("cd %s && git checkout -q %s 2>&1",
                             ref_source, .ref_commit),
                     intern = FALSE)

    if (result != 0) {
      stop("Failed to checkout commit ", .ref_commit)
    }

    message("✓ Reference source downloaded")
  }

  .ref_source_path <<- ref_source
  return(ref_source)
}

# Load reference package using pkgload
load_reference_env <- function() {
  if (!is.null(.ref_env)) {
    return(.ref_env)
  }

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required. Install with: install.packages('pkgload')")
  }

  ref_source <- get_reference_source()

  message("Loading reference package with pkgload...")
  env <- pkgload::load_all(ref_source, export_all = FALSE, quiet = TRUE)

  .ref_env <<- env
  return(env)
}

# Load development package using pkgload
load_development_env <- function() {
  if (!is.null(.dev_env)) {
    return(.dev_env)
  }

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required. Install with: install.packages('pkgload')")
  }

  # Get path to development package (current package being tested)
  # Use rprojroot to find package root
  dev_source <- tryCatch({
    rprojroot::find_root(rprojroot::is_r_package)
  }, error = function(e) {
    # Fallback: assume we're in tests/testthat
    normalizePath(file.path(getwd(), "../.."))
  })

  message("Loading development package with pkgload...")
  env <- pkgload::load_all(dev_source, export_all = FALSE, quiet = TRUE)

  .dev_env <<- env
  return(env)
}

# Skip test if reference not available
skip_if_no_reference <- function() {
  tryCatch({
    load_reference_env()
    load_development_env()
  }, error = function(e) {
    skip(paste("Reference comparison not available:", e$message))
  })
}

# Compare new implementation to reference
compare_to_reference <- function(func_name, args, tolerance = 1e-8, ref_func_name = NULL, ref_args = NULL) {
  skip_if_no_reference()

  # Load both environments
  ref_env <- load_reference_env()
  dev_env <- load_development_env()

  # If ref_func_name not specified, use same name as dev function
  if (is.null(ref_func_name)) {
    ref_func_name <- func_name
  }

  # If ref_args not specified, use same args as dev function
  if (is.null(ref_args)) {
    ref_args <- args
  }

  # Get functions from each environment
  ref_func <- ref_env$env[[ref_func_name]]
  dev_func <- dev_env$env[[func_name]]

  if (is.null(ref_func)) {
    stop("Function '", ref_func_name, "' not found in reference package")
  }

  if (is.null(dev_func)) {
    stop("Function '", func_name, "' not found in development package")
  }

  # Call both implementations (potentially with different arguments)
  dev_result <- do.call(dev_func, args)
  ref_result <- do.call(ref_func, ref_args)

  # Deep comparison of all fields
  expect_equal_susie_objects(dev_result, ref_result, tolerance)

  invisible(list(dev = dev_result, ref = ref_result))
}

# Deep comparison of susie objects
expect_equal_susie_objects <- function(dev_obj, ref_obj, tolerance = 1e-8) {

  # Core posterior quantities
  expect_equal(dev_obj$alpha, ref_obj$alpha, tolerance = tolerance,
               info = "alpha (posterior inclusion probabilities) differ")

  expect_equal(dev_obj$mu, ref_obj$mu, tolerance = tolerance,
               info = "mu (posterior means) differ")

  expect_equal(dev_obj$mu2, ref_obj$mu2, tolerance = tolerance,
               info = "mu2 (posterior second moments) differ")

  # Variance parameters
  expect_equal(dev_obj$V, ref_obj$V, tolerance = tolerance,
               info = "V (prior variance) differs")

  expect_equal(dev_obj$sigma2, ref_obj$sigma2, tolerance = tolerance,
               info = "sigma2 (residual variance) differs")

  # ELBO and convergence
  expect_equal(dev_obj$elbo, ref_obj$elbo, tolerance = tolerance,
               info = "ELBO values differ")

  expect_equal(dev_obj$niter, ref_obj$niter,
               info = "Number of iterations differs")

  expect_equal(dev_obj$converged, ref_obj$converged,
               info = "Convergence status differs")

  # Fitted values and intercept
  if (!is.null(dev_obj$fitted) && !is.null(ref_obj$fitted)) {
    expect_equal(dev_obj$fitted, ref_obj$fitted, tolerance = tolerance,
                 info = "Fitted values differ")
  }

  expect_equal(dev_obj$intercept, ref_obj$intercept, tolerance = tolerance,
               info = "Intercept differs")

  # PIPs (if present)
  if (!is.null(dev_obj$pip) && !is.null(ref_obj$pip)) {
    expect_equal(dev_obj$pip, ref_obj$pip, tolerance = tolerance,
                 info = "PIPs differ")
  }

  # Credible sets (if present)
  if (!is.null(dev_obj$sets) && !is.null(ref_obj$sets)) {
    expect_equal(dev_obj$sets$cs, ref_obj$sets$cs,
                 info = "Credible sets differ")
    if (!is.null(dev_obj$sets$purity) && !is.null(ref_obj$sets$purity)) {
      expect_equal(dev_obj$sets$purity, ref_obj$sets$purity, tolerance = tolerance,
                   info = "CS purity differs")
    }
    expect_equal(dev_obj$sets$coverage, ref_obj$sets$coverage, tolerance = tolerance,
                 info = "CS coverage differs")
  }

  invisible(TRUE)
}

# Compare susie_ss objects
expect_equal_susie_ss_objects <- function(dev_obj, ref_obj, tolerance = 1e-8) {

  # Use the same comparisons as susie objects
  expect_equal_susie_objects(dev_obj, ref_obj, tolerance)

  # Additional checks specific to sufficient statistics
  if (!is.null(dev_obj$XtXr) && !is.null(ref_obj$XtXr)) {
    expect_equal(dev_obj$XtXr, ref_obj$XtXr, tolerance = tolerance,
                 info = "XtXr differs")
  }

  invisible(TRUE)
}

# Compare susie_rss objects
expect_equal_susie_rss_objects <- function(dev_obj, ref_obj, tolerance = 1e-8) {

  # Use the same comparisons as susie objects
  expect_equal_susie_objects(dev_obj, ref_obj, tolerance)

  # Additional checks specific to RSS
  if (!is.null(dev_obj$Rz) && !is.null(ref_obj$Rz)) {
    expect_equal(dev_obj$Rz, ref_obj$Rz, tolerance = tolerance,
                 info = "Rz differs")
  }

  invisible(TRUE)
}
