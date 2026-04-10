#' @title Slot Activity Prior for SuSiE
#'
#' @description Construct a prior specification for the Gamma-Poisson
#'   slot activity model. This model regularizes the number of active
#'   single effects in SuSiE, addressing identifiability between sparse
#'   and dense components.
#'
#' @param C Expected number of causal variants. Acts as the prior mean
#'   for the number of active single-effect slots. Must be positive.
#' @param nu Overdispersion parameter for the Gamma prior on the
#'   per-block causal rate. Larger values give stronger shrinkage
#'   toward C. Default 8 when not specified.
#' @param update_schedule How the Gamma shape parameter is updated
#'   during IBSS iterations. \code{"batch"} updates once per full
#'   sweep (standard CAVI with strict ELBO monotonicity guarantee).
#'   \code{"sequential"} updates after each slot (faster convergence
#'   per iteration, used by susieAnn).
#' @param c_hat_init Optional numeric L-vector of initial slot activity
#'   probabilities for warm-starting. If NULL, initialized uniformly
#'   at C/L.
#' @param skip_threshold_multiplier Multiplier for the adaptive skip
#'   threshold. Slots with c_hat below this fraction of the baseline
#'   (prior with zero signal) are skipped. Default 0 (no skipping).
#'
#' @return A list of class \code{"slot_prior"} with subclass
#'   \code{"slot_prior_poisson"} or \code{"slot_prior_binomial"}.
#'
#' @details
#' Two prior types are available:
#' \describe{
#'   \item{\code{slot_prior_poisson}}{Uses the Poisson approximation
#'     for slot indicators. This is the recommended default and is
#'     used by susieAnn for genome-wide applications.}
#'   \item{\code{slot_prior_binomial}}{Uses the exact Binomial model
#'     for slot indicators. The CAVI update includes a correction term
#'     \eqn{-\log(1 - E[\mu]/L)} that accounts for the finite slot
#'     capacity. No restriction on the ratio L/C is required, but
#'     the correction can inflate slot activity estimates when L is
#'     not much larger than C.}
#' }
#'
#' @examples
#' # Default: Poisson with C=5 active effects expected
#' slot_prior_poisson(C = 5, nu = 8)
#'
#' # Exact Binomial alternative
#' slot_prior_binomial(C = 5, nu = 8)
#'
#' # Pass to susie
#' # fit <- susie(X, y, slot_prior = slot_prior_poisson(C = 5, nu = 8))
#'
#' @export
slot_prior_poisson <- function(C, nu = NULL, update_schedule = c("sequential", "batch"),
                               c_hat_init = NULL, skip_threshold_multiplier = 0) {
  update_schedule <- match.arg(update_schedule)
  stopifnot(is.numeric(C), length(C) == 1, C > 0)
  nu_was_null <- is.null(nu)
  if (nu_was_null) nu <- 8
  stopifnot(is.numeric(nu), length(nu) == 1, nu > 0)
  structure(
    list(
      C = C,
      nu = nu,
      nu_was_default = nu_was_null,
      update_schedule = update_schedule,
      c_hat_init = c_hat_init,
      skip_threshold_multiplier = skip_threshold_multiplier
    ),
    class = c("slot_prior_poisson", "slot_prior")
  )
}

#' @rdname slot_prior_poisson
#' @export
slot_prior_binomial <- function(C, nu = NULL, update_schedule = c("batch", "sequential"),
                                c_hat_init = NULL, skip_threshold_multiplier = 0) {
  update_schedule <- match.arg(update_schedule)
  stopifnot(is.numeric(C), length(C) == 1, C > 0)
  nu_was_null <- is.null(nu)
  if (nu_was_null) nu <- 8
  stopifnot(is.numeric(nu), length(nu) == 1, nu > 0)
  structure(
    list(
      C = C,
      nu = nu,
      nu_was_default = nu_was_null,
      update_schedule = update_schedule,
      c_hat_init = c_hat_init,
      skip_threshold_multiplier = skip_threshold_multiplier
    ),
    class = c("slot_prior_binomial", "slot_prior")
  )
}

#' @export
print.slot_prior <- function(x, ...) {
  type <- if (inherits(x, "slot_prior_binomial")) "binomial" else "poisson"
  cat(sprintf("Slot activity prior (%s)\n", type))
  cat(sprintf("  C (expected causal):  %g\n", x$C))
  cat(sprintf("  nu (overdispersion):  %g\n", x$nu))
  cat(sprintf("  update schedule:      %s\n", x$update_schedule))
  if (!is.null(x$c_hat_init))
    cat(sprintf("  warm start:           %d-vector\n", length(x$c_hat_init)))
  invisible(x)
}

#' Check if an object is a slot_prior
#' @param x Object to test.
#' @return Logical.
#' @keywords internal
#' @noRd
is.slot_prior <- function(x) inherits(x, "slot_prior")
