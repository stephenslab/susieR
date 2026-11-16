#' Control parameters for SuSiE RSS
#'
#' @param mismatch_estimator Estimator for the region-level mismatch variance:
#'   \code{"mle"} or \code{"map"}. With \code{unmappable_effects = "inf"} or
#'   \code{estimate_residual_method = "NIG"}, which estimate the residual
#'   variance, it is always \code{"map"}.
#' @param mixture_reference_p Two-sided reference P-value used to calibrate
#'   prior odds for \code{R_mismatch = "eb_mix"}.
#' @param qc_eigen_tol_rel,qc_eigen_tol_abs Relative and absolute cutoffs for
#'   low-eigenvalue directions in the RSS mismatch diagnostic.
#' @param artifact_threshold Threshold in \code{[0, 1]} for the residual
#'   artifact score.
#' @param sensitivity_threshold Threshold for credible-set Bayes-factor
#'   attenuation.
#' @param r_tol Tolerance for positive-semidefinite matrix checks.
#' @param check_input Whether to check summary-statistic matrices for positive
#'   semidefiniteness and consistency.
#' @param check_prior Whether to check for an unreasonably large estimated
#'   prior variance relative to the marginal z-scores.
#'
#' @return A named list for the \code{control} argument of
#'   \code{\link{susie_rss}}.
#'
#' @export
susie_rss_control <- function(mismatch_estimator = c("mle", "map"),
                              mixture_reference_p = 5e-8,
                              qc_eigen_tol_rel = 1e-3,
                              qc_eigen_tol_abs = 0,
                              artifact_threshold = 0.1,
                              sensitivity_threshold = 30,
                              r_tol = 1e-8,
                              check_input = FALSE,
                              check_prior = TRUE) {
  mismatch_estimator <- match.arg(mismatch_estimator)
  mixture_reference_p <- validate_mixture_reference_p(mixture_reference_p)

  nonnegative <- list(
    qc_eigen_tol_rel = qc_eigen_tol_rel,
    qc_eigen_tol_abs = qc_eigen_tol_abs,
    sensitivity_threshold = sensitivity_threshold,
    r_tol = r_tol
  )
  invalid <- vapply(nonnegative, function(x)
    !is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0,
    logical(1))
  if (any(invalid))
    stop(names(nonnegative)[which(invalid)[1]],
         " must be a single nonnegative finite numeric value.")

  if (!is.numeric(artifact_threshold) || length(artifact_threshold) != 1L ||
      !is.finite(artifact_threshold) || artifact_threshold < 0 ||
      artifact_threshold > 1)
    stop("artifact_threshold must be a single finite numeric value in [0, 1].")

  logicals <- list(check_input = check_input, check_prior = check_prior)
  invalid <- vapply(logicals, function(x)
    !is.logical(x) || length(x) != 1L || is.na(x), logical(1))
  if (any(invalid))
    stop(names(logicals)[which(invalid)[1]], " must be TRUE or FALSE.")

  list(
    mismatch_estimator = mismatch_estimator,
    mixture_reference_p = mixture_reference_p,
    qc_eigen_tol_rel = qc_eigen_tol_rel,
    qc_eigen_tol_abs = qc_eigen_tol_abs,
    artifact_threshold = artifact_threshold,
    sensitivity_threshold = sensitivity_threshold,
    r_tol = r_tol,
    check_input = check_input,
    check_prior = check_prior
  )
}

#' @keywords internal
validate_susie_rss_control <- function(control) {
  if (!is.list(control))
    stop("control must be a list created by susie_rss_control().")
  if (length(control) &&
      (is.null(names(control)) || any(!nzchar(names(control))) ||
       anyDuplicated(names(control))))
    stop("control must be a uniquely named list.")

  unknown <- setdiff(names(control), names(formals(susie_rss_control)))
  if (length(unknown))
    stop("Unknown control parameter", if (length(unknown) > 1L) "s" else "",
         ": ", paste(unknown, collapse = ", "), ".")

  do.call(susie_rss_control, control)
}
