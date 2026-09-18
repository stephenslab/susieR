#' @title Extract regression coefficients from susie fit
#'
#' @param object A susie fit.
#'
#' @param \dots Additional arguments passed to the generic \code{coef}
#'   method.
#'
#' @return A p+1 vector, the first element being an intercept, and the
#'   remaining p elements being estimated regression coefficients.
#'
#' @importFrom stats coef
#'
#' @method coef susie
#'
#' @export coef.susie
#'
#' @export
#'
coef.susie <- function(object, ...) {
  if (inherits(object, "susie_slide")) return(coef.susie_slide(object, ...))
  s <- object

  # Assume unit scale when X_column_scale_factors is absent (susie_rss() /
  # SER-fallback fits), else the divide below yields numeric(0).
  scale <- assume_unit_scale(s$X_column_scale_factors)

  # Compute mappable effects
  mappable_coef <- colSums(s$alpha * s$mu) / scale

  if (!is.null(s$theta)) {
    total_coef <- mappable_coef + s$theta / scale
  } else {
    total_coef <- mappable_coef
  }

  # The intercept must occupy the first slot, or `c(NULL, total_coef)` silently
  # returns length p and every coef(fit)[-1] misaligns against names(pip).
  # susie_rss() fits have no intercept (Zou et al. 2022); return 0 as a
  # placeholder so the length is p + 1. Use isTRUE(is.na(.)) -- a bare is.na()
  # errors on the NULL that saved fits carry.
  if (is.null(s$intercept) || isTRUE(is.na(s$intercept))) {
    warning_message("Intercept is missing or NA; returning 0. susie_rss() ",
                    "fits have no intercept (Zou et al. 2022).", style = "hint")
    icept <- 0
  } else {
    icept <- s$intercept
  }

  return(c(icept, total_coef))
}

#' @title Predict outcomes or extract coefficients from susie fit.
#'
#' @param object A susie fit.
#'
#' @param newx A new value for X at which to do predictions.
#'
#' @param type The type of output. For \code{type = "response"},
#'   predicted or fitted outcomes are returned; for \code{type =
#'   "coefficients"}, the estimated coefficients are returned.
#'
#' @param \dots Other arguments used by generic predict function. These
#'   extra arguments are not used here.
#'
#' @return For \code{type = "response"}, predicted or fitted outcomes
#'   are returned; for \code{type = "coefficients"}, the estimated
#'   coefficients are returned. If the susie fit has intercept =
#'   \code{NA} (which is common when using \code{susie_ss}) then
#'   predictions are computed using an intercept of 0, and a warning is
#'   emitted.
#'
#' @importFrom stats coef
#'
#' @method predict susie
#'
#' @export predict.susie
#'
#' @export
#'
predict.susie <- function(object, newx = NULL,
                          type = c("response", "coefficients"), ...) {
  if (inherits(object, "susie_slide"))
    return(predict.susie_slide(object, newx = newx, type = type, ...))
  s <- object
  type <- match.arg(type)
  if (type == "coefficients") {
    if (!missing(newx)) {
      stop("Do not supply newx when predicting coefficients")
    }
    return(coef(s))
  }
  if (missing(newx)) {
    return(s$fitted)
  }
  if (is.na(s$intercept)) {
    warning_message("The prediction assumes intercept = 0.",
                    style = "hint")
    return(drop(newx %*% coef(s)[-1]))
  } else {
    return(drop(s$intercept + newx %*% coef(s)[-1]))
  }
}
