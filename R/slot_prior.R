#' @title Slot Activity Prior for SuSiE
#'
#' @description Construct a prior specification for the slot activity
#'   model, which regularizes the number of active single effects in
#'   SuSiE. Two prior families are available: Beta-Binomial (default,
#'   recommended for single-locus) and Gamma-Poisson (recommended for
#'   genome-wide applications via susieAnn).
#'
#' @param C Expected number of causal variants. For \code{slot_prior_betabinom},
#'   this determines the Beta prior parameters as \code{a = C, b = L - C}
#'   (set internally when L is known). For Gamma-Poisson variants, it sets
#'   the prior mean for the per-block causal rate. Must be positive.
#' @param nu Overdispersion parameter for the Gamma-Poisson prior on the
#'   per-block causal rate. Not used by \code{slot_prior_betabinom}.
#'   Larger values give stronger shrinkage toward C. Default 8 when
#'   not specified.
#' @param a_beta Shape parameter for the Beta prior on inclusion
#'   probability rho. Default 1.
#' @param b_beta Shape parameter for the Beta prior on inclusion
#'   probability rho. Default 1 (uniform prior on rho).
#'   When \code{a_beta = 1} and \code{b_beta = 1}, the prior is
#'   uniform on [0,1], providing automatic multiplicity correction
#'   following Scott and Berger (2010).
#' @param update_schedule How the sufficient statistics are updated
#'   during IBSS iterations. \code{"batch"} updates once per full
#'   sweep (standard CAVI). \code{"sequential"} updates after each
#'   slot (faster convergence per iteration, used by susieAnn).
#' @param c_hat_init Optional numeric L-vector of initial slot activity
#'   probabilities for warm-starting. If NULL, initialized at the
#'   prior mean.
#' @param skip_threshold_multiplier Multiplier for the adaptive skip
#'   threshold. Slots with c_hat below this fraction of the baseline
#'   (prior with zero signal) are skipped. Default 0 (no skipping).
#'   The threshold is recomputed after each sweep from the current
#'   model state, and is set to 0 on the first sweep so all slots
#'   are evaluated at least once.
#'
#' @return A list of class \code{"slot_prior"} with the appropriate
#'   subclass.
#'
#' @details
#' Two prior types are available:
#' \describe{
#'   \item{\code{slot_prior_betabinom}}{Uses a Beta-Binomial model
#'     for slot inclusion. The inclusion probability rho is given a
#'     Beta(a_beta, b_beta) prior and integrated out analytically,
#'     yielding an adaptive multiplicity correction that penalizes
#'     less when more slots are active. This is the recommended
#'     default for single-locus applications. See Scott and Berger
#'     (2010) for the theoretical justification.}
#'   \item{\code{slot_prior_poisson}}{Uses the Gamma-Poisson model
#'     with Poisson approximation for slot indicators. Recommended
#'     for genome-wide applications via susieAnn, where C and nu
#'     are estimated across loci.}
#' }
#'
#' @references
#' Scott, J. G. and Berger, J. O. (2010). Bayes and empirical-Bayes
#' multiplicity adjustment in the variable-selection problem.
#' \emph{Annals of Statistics}, 38(5), 2587--2619.
#'
#' @examples
#' # Default: Beta-Binomial with Beta(1, 1.2) prior on inclusion probability
#' slot_prior_betabinom()
#'
#' # Gamma-Poisson for susieAnn
#' slot_prior_poisson(C = 5, nu = 8)
#'
#' # Pass to susie
#' # fit <- susie(X, y, slot_prior = slot_prior_betabinom())
#'
# Beta-Binomial slot activity model.
# rho ~ Beta(a_beta, b_beta), c_l | rho ~ Bernoulli(rho).
# Marginalizing rho gives the collapsed posterior update:
#   log(c_l/(1-c_l)) = log(a + k_{-l}) - log(b + L-1 - k_{-l}) + lbf_l
# where k_{-l} is the soft count of other active slots.
# This provides automatic multiplicity correction that adapts to the
# current number of active effects, following:
#   Scott, J.G. and Berger, J.O. (2010). Bayes and empirical-Bayes
#   multiplicity adjustment in the variable-selection problem.
#   Annals of Statistics, 38(5), 2587-2619.
#
# Default: Beta(1, 1.2) gives mild sparsity (b > 1) compared to assuming
# as many as > 50% slots active (uniform). 
#' @export
slot_prior_betabinom <- function(a_beta = NULL, b_beta = NULL,
                                 update_schedule = c("sequential", "batch"),
                                 c_hat_init = NULL,
                                 skip_threshold_multiplier = 0) {
  update_schedule <- match.arg(update_schedule)

  # Default a_beta = 1 (standard reference prior for Bernoulli).
  ab_was_default <- is.null(a_beta) && is.null(b_beta)
  if (is.null(a_beta)) a_beta <- 1
  # Default b_beta = 1.2: gives E[rho] = 1/(1+1.2) = 0.45, i.e. ~45% of
  # slots expected active. With L=10 this means ~4 to 5 active effects.
  # The density p(rho) ~ (1-rho)^0.5 is a mild sparsity preference.
  # In the collapsed Beta-Binomial update, b > 1 penalizes inactive slots
  # more, reducing residual contamination from partially-weighted effects.
  if (is.null(b_beta)) b_beta <- 1.2

  stopifnot(is.numeric(a_beta), length(a_beta) == 1, a_beta > 0)
  stopifnot(is.numeric(b_beta), length(b_beta) == 1, b_beta > 0)
  structure(
    list(
      a_beta = a_beta,
      b_beta = b_beta,
      ab_was_default = ab_was_default,
      update_schedule = update_schedule,
      c_hat_init = c_hat_init,
      skip_threshold_multiplier = skip_threshold_multiplier
    ),
    class = c("slot_prior_betabinom", "slot_prior")
  )
}

#' @rdname slot_prior_betabinom
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

# slot_prior_binomial removed (merged into slot_prior_poisson).
# Use slot_prior_betabinom (single-locus) or slot_prior_poisson (genome-wide).

#' @export
print.slot_prior <- function(x, ...) {
  type <- if (inherits(x, "slot_prior_betabinom")) "beta-binomial"
          else "poisson"
  cat(sprintf("Slot activity prior (%s)\n", type))
  if (type == "beta-binomial") {
    cat(sprintf("  a_beta:               %g\n", x$a_beta))
    cat(sprintf("  b_beta:               %g\n", x$b_beta))
  } else {
    cat(sprintf("  C (expected causal):  %g\n", x$C))
    cat(sprintf("  nu (overdispersion):  %g\n", x$nu))
  }
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
