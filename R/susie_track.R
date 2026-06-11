#' Make a compact SuSiE track history
#'
#' \code{make_susie_track_history} organizes the compact fitting history
#' produced by \code{track_fit = TRUE}. The alpha history is stored in sparse
#' long format rather than as dense \code{L x p} matrices.
#'
#' @param fit A fitted \code{"susie"} object with compact track snapshots.
#'
#' @return A \code{"susie_track"} object.
#'
#' @noRd
# nocov start
make_susie_track_history <- function(fit) {
  if (inherits(fit, "susie_track"))
    return(fit)
  if (!inherits(fit, "susie"))
    stop("fit must be a fitted susie object")
  if (is.null(fit$trace))
    stop("fit does not contain track information; refit with track_fit = TRUE")
  if (inherits(fit$trace, "susie_track"))
    return(fit$trace)
  if (!is_compact_track_snapshots(fit$trace))
    stop("fit$trace is not a compact SuSiE track")

  alpha <- do.call(rbind, lapply(fit$trace, `[[`, "alpha"))
  effect <- do.call(rbind, lapply(fit$trace, `[[`, "effect"))
  iteration <- do.call(rbind, lapply(fit$trace, `[[`, "iteration"))
  rownames(alpha) <- rownames(effect) <- rownames(iteration) <- NULL
  variable_names <- names(fit$pip)
  if (is.null(variable_names))
    variable_names <- colnames(fit$alpha)
  if (!is.null(variable_names)) {
    if (nrow(alpha) > 0)
      alpha$variable_name <- variable_names[alpha$variable]
    if (nrow(effect) > 0)
      effect$top_variable_name <- variable_names[effect$top_variable]
  }

  out <- list(
    alpha = alpha,
    effect = effect,
    iteration = iteration,
    diagnosis = make_track_diagnosis(fit),
    cs_sensitivity = make_track_cs_sensitivity(fit),
    meta = make_track_meta(fit)
  )
  class(out) <- c("susie_track", "list")
  out
}

#' Restore a dense alpha matrix from a compact SuSiE track history
#'
#' @param track A \code{"susie_track"} object.
#' @param iteration Iteration to restore.
#'
#' @return An \code{L x p} alpha matrix with omitted entries filled with zero.
#'
#' @noRd
restore_alpha_from_track <- function(track, iteration) {
  if (!inherits(track, "susie_track"))
    stop("track must be a susie_track object")
  if (missing(iteration) || length(iteration) != 1)
    stop("provide a single iteration to restore")

  alpha <- matrix(0, nrow = track$meta$L, ncol = track$meta$p)
  if (!is.null(track$meta$variable_names))
    colnames(alpha) <- track$meta$variable_names

  x <- track$alpha[track$alpha$iteration == iteration, , drop = FALSE]
  if (nrow(x) > 0)
    alpha[cbind(x$effect, x$variable)] <- x$alpha
  alpha
}
# nocov end

make_track_snapshot <- function(model, iteration) {
  p <- ncol(model$alpha)
  min_max_alpha <- 2 / p
  min_alpha <- 1 / p
  variable_names <- colnames(model$alpha)
  if (is.null(variable_names))
    variable_names <- as.character(seq_len(p))
  summarize_track_snapshot(model, iteration, variable_names, min_max_alpha,
                           min_alpha)
}

summarize_track_snapshot <- function(model, iteration, variable_names,
                                     min_max_alpha, min_alpha) {
  alpha <- model$alpha
  L <- nrow(alpha)
  row_max <- apply(alpha, 1, max)
  row_keep <- row_max >= min_max_alpha
  pip <- 1 - apply(1 - alpha, 2, prod)
  top <- max.col(alpha, ties.method = "first")
  top_alpha <- alpha[cbind(seq_len(L), top)]

  alpha_rows <- lapply(which(row_keep), function(l) {
    vars <- which(alpha[l, ] >= min_alpha)
    if (length(vars) == 0)
      return(NULL) # nocov
    data.frame(iteration = iteration, effect = l, variable = vars,
               variable_name = variable_names[vars],
               alpha = as.numeric(alpha[l, vars]),
               stringsAsFactors = FALSE)
  })
  alpha_df <- do.call(rbind, alpha_rows)
  if (is.null(alpha_df))
    alpha_df <- data.frame(iteration = integer(0), effect = integer(0),
                           variable = integer(0), variable_name = character(0),
                           alpha = numeric(0), stringsAsFactors = FALSE)

  effect <- data.frame(
    iteration = iteration,
    effect = seq_len(L),
    V = track_vector(model$V, L),
    sigma2 = track_scalar(model$sigma2),
    lambda_bias = track_scalar(model$lambda_bias),
    effect_lbf = track_vector(model$lbf, L),
    top_variable = top,
    top_variable_name = variable_names[top],
    top_alpha = top_alpha,
    top_lbf = track_top_lbf(model, top),
    max_lbf_variable = track_max_lbf(model, L),
    kept_effect = row_keep,
    n_alpha_above_background = rowSums(alpha >= min_alpha),
    alpha_mass_above_background = rowSums(alpha * (alpha >= min_alpha)),
    stringsAsFactors = FALSE
  )
  if (!is.null(model$slot_weights)) # nocov start
    effect$slot_weight <- track_vector(model$slot_weights, L) # nocov end
  if (!is.null(model$tau2)) # nocov start
    effect$tau2 <- track_vector(model$tau2, L) # nocov end

  iteration_df <- data.frame(
    iteration = iteration,
    sigma2 = track_scalar(model$sigma2),
    lambda_bias = track_scalar(model$lambda_bias),
    max_alpha = max(alpha),
    n_effects_kept = sum(row_keep),
    n_effects_alpha_gt_half = sum(row_max > 0.5),
    max_pip = max(pip),
    n_pip_gt_0.1 = sum(pip > 0.1),
    n_alpha_entries = nrow(alpha_df),
    stringsAsFactors = FALSE
  )

  list(alpha = alpha_df, effect = effect, iteration = iteration_df)
}

make_track_meta <- function(fit) {
  p <- length(fit$pip)
  if (is.null(p) || p == 0) # nocov start
    p <- ncol(fit$alpha)    # nocov end
  variable_names <- names(fit$pip)
  if (is.null(variable_names))
    variable_names <- colnames(fit$alpha)
  list(
    p = p,
    L = nrow(fit$alpha),
    niter = fit$niter,
    min_max_alpha = 2 / p,
    min_alpha = 1 / p,
    variable_names = variable_names,
    alpha_is_thresholded = TRUE
  )
}

is_compact_track_snapshots <- function(x) {
  is.list(x) && length(x) > 0 && all(vapply(x, function(y) {
    is.list(y) && is.data.frame(y$alpha) && is.data.frame(y$effect) &&
      is.data.frame(y$iteration)
  }, logical(1)))
}

track_top_lbf <- function(model, top) {
  if (is.null(model$lbf_variable))
    return(rep(NA_real_, length(top)))
  model$lbf_variable[cbind(seq_along(top), top)]
}

track_max_lbf <- function(model, L) {
  if (is.null(model$lbf_variable))
    return(rep(NA_real_, L))
  apply(model$lbf_variable, 1, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  })
}

track_vector <- function(x, n) {
  if (is.null(x))
    return(rep(NA_real_, n))
  x <- as.vector(x)
  if (length(x) >= n)
    return(x[seq_len(n)])
  c(x, rep(NA_real_, n - length(x))) # nocov
}

track_scalar <- function(x) {
  if (is.null(x) || length(x) != 1)
    return(NA_real_)
  as.numeric(x)
}

make_track_diagnosis <- function(fit) {
  d <- fit$R_finite_diagnostics
  if (is.null(d))
    return(NULL)

  iter_diag <- NULL
  tr <- d$R_mismatch_trace
  if (!is.null(tr) && length(tr) > 0) {
    tr <- tr[vapply(tr, function(x) is.null(x$phase) || x$phase != "init_ser",
                   logical(1))]
    if (length(tr) > 0) {
      iter_diag <- do.call(rbind, lapply(seq_along(tr), function(i) {
        x <- tr[[i]]
        data.frame(
          iteration = i,
          lambda_bias = scalar_from_list(x, "lambda_bias"),
          B_effective = scalar_from_list(x, "B_corrected"),
          B = scalar_from_list(x, "B"),
          sigma2 = scalar_from_list(x, "sigma2"),
          mean_r2 = scalar_from_list(x, "mean_r2"),
          median_r2 = scalar_from_list(x, "median_r2"),
          max_r2 = scalar_from_list(x, "max_r2"),
          mean_s = scalar_from_list(x, "mean_s"),
          median_s = scalar_from_list(x, "median_s"),
          max_s = scalar_from_list(x, "max_s"),
          cor_r2_s = scalar_from_list(x, "cor_r2_s"),
          Q_art = scalar_from_list(x, "Q_art"),
          artifact_flag = logical_from_list(x, "artifact_flag"),
          low_eigen_count = scalar_from_list(x, "low_eigen_count"),
          stringsAsFactors = FALSE)
      }))
      rownames(iter_diag) <- NULL
    }
  }

  final <- data.frame(
    lambda_bias = scalar_from_list(d, "lambda_bias"),
    B_effective = scalar_from_list(d, "B_corrected"),
    B = scalar_from_list(d, "B"),
    Q_art = scalar_from_list(d, "Q_art"),
    artifact_flag = logical_from_list(d, "artifact_flag"),
    R_sensitivity_flag = logical_from_list(d, "R_sensitivity_flag"),
    R_reliability_flag = logical_from_list(d, "R_reliability_flag"),
    low_eigen_count = scalar_from_list(d, "low_eigen_count"),
    stringsAsFactors = FALSE)

  list(iteration = iter_diag, final = final)
}

# nocov start
make_track_cs_sensitivity <- function(fit) {
  d <- fit$R_finite_diagnostics
  ba <- if (!is.null(d)) d$bf_attenuation else NULL
  if (is.null(ba) || is.null(fit$sets) || is.null(fit$sets$cs) ||
      length(fit$sets$cs) == 0)
    return(NULL)

  cs_names <- names(fit$sets$cs)
  if (is.null(cs_names))
    cs_names <- paste0("L", seq_along(fit$sets$cs))
  cs_index <- fit$sets$cs_index
  if (is.null(cs_index))
    cs_index <- as.integer(sub("^L", "", cs_names))
  threshold <- if (!is.null(ba$threshold)) ba$threshold else NA_real_
  variable_names <- names(fit$pip)

  out <- do.call(rbind, lapply(seq_along(fit$sets$cs), function(i) {
    vars <- fit$sets$cs[[i]]
    top_var <- value_by_name(ba$cs_top_variable, cs_names[i])
    top_atten <- value_by_name(ba$cs_top_attenuation, cs_names[i])
    cs_max <- value_by_name(ba$cs_max, cs_names[i])
    data.frame(
      cs = cs_names[i],
      effect = cs_index[i],
      cs_size = length(vars),
      cs_max_bf_attenuation = cs_max,
      cs_weighted_bf_attenuation = value_by_name(ba$cs_weighted, cs_names[i]),
      cs_sensitivity_label = as.character(value_by_name(ba$cs_label,
                                                        cs_names[i])),
      is_sensitive = is.finite(cs_max) && is.finite(threshold) &&
        cs_max >= threshold,
      threshold = threshold,
      top_sensitive_variable = top_var,
      top_sensitive_variable_name = if (!is.null(variable_names) &&
                                       !is.na(top_var)) variable_names[top_var]
                                    else NA_character_,
      top_sensitive_variable_attenuation = top_atten,
      stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}
# nocov end

scalar_from_list <- function(x, name) {
  if (is.null(x[[name]]) || length(x[[name]]) != 1)
    return(NA_real_) # nocov
  as.numeric(x[[name]])
}

logical_from_list <- function(x, name) {
  if (is.null(x[[name]]) || length(x[[name]]) != 1)
    return(NA)
  as.logical(x[[name]])
}

# nocov start
value_by_name <- function(x, name) {
  if (is.null(x) || length(x) == 0)
    return(NA)
  if (!is.null(names(x)) && name %in% names(x))
    return(x[[name]])
  NA
}
# nocov end
