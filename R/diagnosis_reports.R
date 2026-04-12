# Diagnostic functions for SuSiE-ash filter
#
# Per-iteration functions (called from susie_utils.R, return data.frame):
#   diagnose_bb_ash_iter()              - BB+ash code path
#   diagnose_ash_filter_archived_iter() - V0 code path
#
# Post-run helpers (called via susieR:::):
#   collect_ash_diag(fit)               - rbind all iterations into ML table
#   label_diag_truth(df, fit, causal)   - add TP/FP labels
#   add_delta_features(df)              - add per-slot change-over-iteration features
#   extract_bb_ash_features(fit, X, causal) - quick feature extraction from converged fit
#   compare_ash_methods(df1, df2)       - side-by-side comparison
#
# Usage example:
#   fit <- susie(X, y, L=10, slot_prior=slot_prior_betabinom(),
#                unmappable_effects="ash", max_iter=50)
#   df <- susieR:::collect_ash_diag(fit)
#   df <- susieR:::label_diag_truth(df, fit, causal)
#   df <- susieR:::add_delta_features(df)
#
#   # ML analysis: per-slot, per-iteration features with TP/FP labels
#   # Key columns: iter, slot, c_hat, lbf, purity, V, max_alpha,
#   #   alpha_entropy, mask_tier, collision, ever_uncertain, ...
#   # Delta columns: delta_c_hat, delta_V, delta_lbf, delta_max_alpha, ...
#   #   (change from previous iteration for the same slot)
#   #
#   # For 4-way comparison (BB+ash/V0 x mrash/no-mrash):
#   #   options(susie.skip_mrash = TRUE)  # toggle mr.ash off
#   #   fit_nomrash <- susie(...)
#   #   options(susie.skip_mrash = FALSE) # restore
#
# Data.frames accumulated on fit$.diag_env$history during the run.
# Debug flag .ash_debug in susie_utils.R (TRUE = on, never turn off).


#' BB+ash per-iteration diagnostic
#'
#' @return data.frame with one row per slot, all features
#' @keywords internal
diagnose_bb_ash_iter <- function(model, Xcorr, mask, b_confident,
                                 sentinels, sentinel_collision,
                                 is_confident_now, is_trusted,
                                 emerging_slots, active_slots, c_hat,
                                 ash_result, p,
                                 high_chat = NULL, low_chat = NULL,
                                 # Tunable parameters (captured for reproducibility)
                                 collision_threshold = 0.9,
                                 purity_threshold = 0.5,
                                 masking_threshold = 0.5,
                                 nPIP_threshold = 0.05,
                                 c_hat_mask_threshold = 0.9) {
  L <- nrow(model$alpha)
  theta_raw <- ash_result$beta
  theta_masked <- theta_raw
  theta_masked[mask] <- 0
  cs_coverage <- 0.9
  iter <- model$ash_iter

  rows <- list()
  for (l in seq_len(L)) {
    sent <- sentinels[l]
    alpha_l <- model$alpha[l, ]
    max_a <- max(alpha_l)

    # Purity and CS size
    cs_order <- order(alpha_l, decreasing = TRUE)
    cs_size <- min(which(cumsum(alpha_l[cs_order]) >= cs_coverage))
    if (cs_size > 1 && l %in% active_slots) {
      cs_idx <- cs_order[1:cs_size]
      pur <- min(abs(Xcorr[cs_idx, cs_idx]))
    } else if (cs_size == 1) {
      pur <- 1.0
    } else {
      pur <- 0.0
    }

    # Status
    status <- if (is_trusted[l]) "trusted"
              else if (is_confident_now[l] && model$ever_uncertain[l]) "conf_unc"
              else if (is_confident_now[l]) "confident"
              else if (sentinel_collision[l]) "collide"
              else if (model$V[l] == 0) "null"
              else "emerging"

    # Mask tier (key c_hat feature)
    mask_tier <- if (is_trusted[l]) "trusted"
                 else if (!is.null(high_chat) && l %in% high_chat) "high_chat"
                 else if (!is.null(low_chat) && l %in% low_chat) "low_chat"
                 else "unknown"

    # Max cross-sentinel |r| (collision strength)
    max_cross_r <- 0
    if (sent > 0 && length(active_slots) > 1) {
      other_sents <- sentinels[setdiff(active_slots, l)]
      other_sents <- other_sents[other_sents > 0]
      if (length(other_sents) > 0)
        max_cross_r <- max(abs(Xcorr[sent, other_sents]))
    }

    # Sentinel change
    prev_sent_l <- if (!is.null(model$prev_sentinel)) model$prev_sentinel[l] else 0L
    sent_changed <- (prev_sent_l > 0) && (sent != prev_sent_l)

    # Theta at sentinel
    theta_at_sent <- if (sent > 0) theta_masked[sent] else 0
    theta_raw_at_sent <- if (sent > 0) theta_raw[sent] else 0

    # Alpha entropy (low = concentrated, possibly FP)
    alpha_nz <- alpha_l[alpha_l > 1e-10]
    alpha_entropy <- -sum(alpha_nz * log(alpha_nz))

    # Number of colliding partners
    n_colliding <- 0
    if (sent > 0 && length(active_slots) > 1) {
      other_sents <- sentinels[setdiff(active_slots, l)]
      other_sents <- other_sents[other_sents > 0]
      if (length(other_sents) > 0)
        n_colliding <- sum(abs(Xcorr[sent, other_sents]) > 0.9)
    }

    # Per-slot mu properties
    mu_l <- model$mu[l, ]
    mu_at_sent <- if (sent > 0) mu_l[sent] else 0
    max_abs_mu <- max(abs(mu_l))

    rows[[l]] <- data.frame(
      method = "bb_ash", iter = iter, slot = l,
      sentinel = sent, purity = pur, V = model$V[l],
      c_hat = c_hat[l], lbf = if (!is.null(model$lbf)) model$lbf[l] else NA,
      max_alpha = max_a, cs_size = cs_size,
      alpha_entropy = alpha_entropy,
      is_confident = is_confident_now[l],
      is_trusted = is_trusted[l],
      status = status, mask_tier = mask_tier,
      was_low_chat = if (!is.null(model$was_low_chat)) model$was_low_chat[l] else FALSE,
      was_exposed = if (!is.null(model$was_exposed)) model$was_exposed[l] else FALSE,
      collision = sentinel_collision[l],
      ever_uncertain = model$ever_uncertain[l],
      n_colliding = n_colliding,
      max_cross_r = max_cross_r,
      sent_changed = sent_changed,
      prev_sentinel = prev_sent_l,
      mu_at_sent = mu_at_sent,
      max_abs_mu = max_abs_mu,
      theta_at_sent = theta_at_sent,
      theta_raw_at_sent = theta_raw_at_sent,
      mask_size = sum(mask), mask_frac = round(sum(mask) / p, 3),
      n_active = length(active_slots),
      n_trusted = sum(is_trusted),
      n_high_chat = if (!is.null(high_chat)) length(high_chat) else NA,
      n_low_chat = if (!is.null(low_chat)) length(low_chat) else NA,
      C_hat = round(sum(c_hat), 1),
      sigma2 = ash_result$sigma2,
      pi0 = if (!is.null(ash_result$pi)) ash_result$pi[1] else NA,
      theta_ss = sum(theta_masked^2),
      theta_raw_ss = sum(theta_raw^2),
      b_conf_ss = sum(b_confident^2),
      b_conf_max = max(abs(b_confident)),
      sent_masked = if (sent > 0) mask[sent] else FALSE,
      skip_mrash = getOption("susie.skip_mrash", FALSE),
      # Rule activation: which decision rules kicked in for this slot
      rule_active_gate = (l %in% active_slots),
      rule_collision = sentinel_collision[l],
      rule_jump = model$ever_uncertain[l] && !sentinel_collision[l],
      rule_trusted = is_trusted[l],
      rule_low_chat = if (!is.null(model$was_low_chat)) model$was_low_chat[l] else FALSE,
      rule_high_chat_pip = (!is.null(high_chat) && l %in% high_chat),
      rule_low_chat_sentinel = (!is.null(low_chat) && l %in% low_chat &&
        !(if (!is.null(model$was_exposed)) model$was_exposed[l] else FALSE)),
      rule_exposed = if (!is.null(model$was_exposed)) model$was_exposed[l] else FALSE,
      # Tunable parameter values
      param_collision_threshold = collision_threshold,
      param_purity_threshold = purity_threshold,
      param_masking_threshold = masking_threshold,
      param_nPIP_threshold = nPIP_threshold,
      param_c_hat_mask_threshold = c_hat_mask_threshold,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


#' V0 archived filter per-iteration diagnostic
#'
#' @return data.frame with one row per slot, all features
#' @keywords internal
diagnose_ash_filter_archived_iter <- function(model, Xcorr, masked,
                                              b_confident, sentinels,
                                              effect_purity, current_case,
                                              current_collision,
                                              mrash_output) {
  L <- nrow(model$alpha)
  p <- ncol(model$alpha)
  theta_raw <- mrash_output$beta
  theta_masked <- theta_raw
  theta_masked[masked] <- 0
  cs_coverage <- 0.9
  iter <- model$ash_iter

  is_active <- sapply(seq_len(L), function(l) {
    max(model$alpha[l, ]) - min(model$alpha[l, ]) >= 5e-5
  })

  rows <- list()
  for (l in seq_len(L)) {
    sent <- sentinels[l]
    alpha_l <- model$alpha[l, ]
    max_a <- max(alpha_l)
    cs_order <- order(alpha_l, decreasing = TRUE)
    cs_size <- min(which(cumsum(alpha_l[cs_order]) >= cs_coverage))
    case_str <- if (!is.na(current_case[l])) paste0("C", current_case[l]) else "inactive"

    # Max cross-sentinel |r|
    max_cross_r <- 0
    active_idx <- which(is_active)
    if (sent > 0 && length(active_idx) > 1) {
      other_sents <- sentinels[setdiff(active_idx, l)]
      other_sents <- other_sents[other_sents > 0]
      if (length(other_sents) > 0)
        max_cross_r <- max(abs(Xcorr[sent, other_sents]))
    }

    # Sentinel change
    prev_sent_l <- if (!is.null(model$prev_sentinel)) model$prev_sentinel[l] else 0L
    sent_changed <- (prev_sent_l > 0) && (sent != prev_sent_l)

    # Theta at sentinel
    theta_at_sent <- if (sent > 0) theta_masked[sent] else 0
    theta_raw_at_sent <- if (sent > 0) theta_raw[sent] else 0

    # Alpha entropy
    alpha_nz <- alpha_l[alpha_l > 1e-10]
    alpha_entropy <- -sum(alpha_nz * log(alpha_nz))

    # Number of colliding partners
    n_colliding <- 0
    if (sent > 0 && length(active_idx) > 1) {
      other_sents <- sentinels[setdiff(active_idx, l)]
      other_sents <- other_sents[other_sents > 0]
      if (length(other_sents) > 0)
        n_colliding <- sum(abs(Xcorr[sent, other_sents]) > 0.9)
    }

    # Per-slot mu properties
    mu_l <- model$mu[l, ]
    mu_at_sent <- if (sent > 0) mu_l[sent] else 0
    max_abs_mu <- max(abs(mu_l))

    rows[[l]] <- data.frame(
      method = "v0", iter = iter, slot = l,
      sentinel = sent, purity = effect_purity[l], V = model$V[l],
      lbf = model$lbf[l], max_alpha = max_a, cs_size = cs_size,
      alpha_entropy = alpha_entropy,
      status = case_str,
      current_collision = current_collision[l],
      ever_diffuse = model$ever_diffuse[l],
      diffuse_iter_count = if (!is.null(model$diffuse_iter_count)) model$diffuse_iter_count[l] else 0L,
      prev_case = if (!is.null(model$prev_case)) model$prev_case[l] else 0L,
      n_colliding = n_colliding,
      max_cross_r = max_cross_r,
      sent_changed = sent_changed,
      prev_sentinel = prev_sent_l,
      mu_at_sent = mu_at_sent,
      max_abs_mu = max_abs_mu,
      theta_at_sent = theta_at_sent,
      theta_raw_at_sent = theta_raw_at_sent,
      mask_size = sum(masked), mask_frac = round(sum(masked) / p, 3),
      n_active = sum(is_active),
      sigma2 = mrash_output$sigma2,
      pi0 = if (!is.null(mrash_output$pi)) mrash_output$pi[1] else NA,
      theta_ss = sum(theta_masked^2),
      theta_raw_ss = sum(theta_raw^2),
      b_conf_ss = sum(b_confident^2),
      b_conf_max = max(abs(b_confident)),
      sent_masked = if (sent > 0) masked[sent] else FALSE,
      skip_mrash = getOption("susie.skip_mrash", FALSE),
      # Rule activation: which V0 decision rules kicked in
      rule_collision = current_collision[l],
      rule_ever_diffuse = (model$ever_diffuse[l] > 0),
      rule_case1 = (!is.na(current_case[l]) && current_case[l] == 1),
      rule_case2 = (!is.na(current_case[l]) && current_case[l] == 2),
      rule_case3 = (!is.na(current_case[l]) && current_case[l] == 3),
      rule_exposure = (if (!is.null(model$diffuse_iter_count)) model$diffuse_iter_count[l] else 0) >= 2,
      rule_second_chance = if (!is.null(model$second_chance_used)) model$second_chance_used[sent] else FALSE,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


# ---- Helper functions for ML analysis ----

#' Collect diagnostic data.frames across iterations
#'
#' Call this after running susie() to rbind all per-iteration diagnostics
#' into a single ML-ready data.frame.
#'
#' @param fit SuSiE fit object (must have been run with .ash_debug = TRUE)
#' @return data.frame with nrow = L * n_ash_iters, or NULL if no diagnostics
#'
#' @examples
#' \dontrun{
#' # Full ML pipeline:
#' data(unmappable_data)
#' X <- unmappable_data$X; y <- as.vector(unmappable_data$y)
#' causal <- which(unmappable_data$beta != 0)
#'
#' fit <- susie(X, y, L=10, slot_prior=slot_prior_betabinom(),
#'              unmappable_effects="ash", max_iter=50)
#'
#' df <- susieR:::collect_ash_diag(fit)         # all iterations
#' df <- susieR:::label_diag_truth(df, fit, causal) # TP/FP labels
#' df <- susieR:::add_delta_features(df)         # temporal features
#'
#' # Inspect FP slot across iterations:
#' subset(df, cs_label == "FP", select = c(iter, slot, sentinel,
#'        c_hat, lbf, max_alpha, alpha_entropy, mask_tier, delta_c_hat))
#'
#' # 4-way comparison (BB+ash vs V0, with/without mr.ash):
#' options(susie.skip_mrash = TRUE)
#' fit_nomrash <- susie(X, y, L=10, slot_prior=slot_prior_betabinom(),
#'                      unmappable_effects="ash", max_iter=50)
#' options(susie.skip_mrash = FALSE)
#' df_nomrash <- susieR:::collect_ash_diag(fit_nomrash)
#'
#' # Decision tree analysis:
#' # library(rpart)
#' # last_iter <- df[df$iter == max(df$iter) & df$V > 0, ]
#' # tree <- rpart(cs_label ~ c_hat + lbf + max_alpha + alpha_entropy +
#' #               purity + mask_tier + delta_c_hat, data = last_iter)
#' }
#'
#' @keywords internal
collect_ash_diag <- function(fit) {
  if (!is.null(fit$.diag_env) && !is.null(fit$.diag_env$history)) {
    return(do.call(rbind, fit$.diag_env$history))
  }
  # Fallback to old .diag_history list
  if (!is.null(fit$.diag_history)) {
    return(do.call(rbind, fit$.diag_history))
  }
  return(NULL)
}


#' Label diagnostic table with ground truth TP/FP
#'
#' For each slot at the final iteration, check if its CS (if any) contains
#' a causal variant.
#'
#' @param df Diagnostic data.frame (from collect_ash_diag or single iter)
#' @param fit SuSiE fit object
#' @param causal Integer vector of causal variant indices
#' @return df with added 'cs_label' column: "TP", "FP", or "-" (no CS)
#' @keywords internal
label_diag_truth <- function(df, fit, causal) {
  cs <- fit$sets$cs
  L <- max(df$slot)
  # Map each CS to its owning slot
  cs_slot_map <- rep(NA, L)
  cs_tp_map <- rep(NA, L)
  if (length(cs) > 0) {
    for (i in seq_along(cs)) {
      sent <- cs[[i]][which.max(fit$pip[cs[[i]]])]
      owner <- which.max(fit$alpha[, sent])
      cs_slot_map[owner] <- i
      cs_tp_map[owner] <- any(cs[[i]] %in% causal)
    }
  }
  df$cs_label <- sapply(df$slot, function(l) {
    if (is.na(cs_slot_map[l])) "-"
    else if (cs_tp_map[l]) "TP"
    else "FP"
  })
  df
}


#' Add per-slot delta features (change from previous iteration)
#'
#' Computes delta_c_hat, delta_V, delta_lbf, delta_max_alpha,
#' delta_alpha_entropy, delta_purity for each slot across iterations.
#' Also adds lag1 features (previous iteration values) and
#' cumulative features (max c_hat ever, min alpha_entropy ever).
#' These temporal features help ML models detect trajectories
#' (e.g., a slot that was strong then weakened = collapse signal).
#'
#' @param df data.frame from collect_ash_diag + label_diag_truth
#' @return df with added delta_, lag1_, and cum_ columns
#'
#' @examples
#' \dontrun{
#' df <- susieR:::collect_ash_diag(fit)
#' df <- susieR:::label_diag_truth(df, fit, causal)
#' df <- susieR:::add_delta_features(df)
#' # Now df has delta_c_hat, lag1_c_hat, cum_max_c_hat, etc.
#' # Use for decision tree: rpart::rpart(cs_label ~ ., data = df_last_iter)
#' }
#'
#' @keywords internal
add_delta_features <- function(df) {
  iters <- sort(unique(df$iter))
  slots <- sort(unique(df$slot))

  # Features to compute deltas/lags for
  feat_cols <- c("c_hat", "V", "lbf", "max_alpha", "alpha_entropy", "purity",
                 "mask_size", "theta_ss")

  # Initialize new columns
  for (f in feat_cols) {
    df[[paste0("delta_", f)]] <- NA_real_
    df[[paste0("lag1_", f)]] <- NA_real_
  }
  df$cum_max_c_hat <- NA_real_
  df$cum_min_entropy <- NA_real_
  df$cum_max_lbf <- NA_real_

  for (s in slots) {
    idx <- which(df$slot == s)
    if (length(idx) < 2) next

    slot_df <- df[idx, ]
    for (f in feat_cols) {
      vals <- slot_df[[f]]
      if (all(is.na(vals))) next
      # Delta: current - previous
      delta <- c(NA, diff(vals))
      df[[paste0("delta_", f)]][idx] <- delta
      # Lag1: previous value
      lag1 <- c(NA, vals[-length(vals)])
      df[[paste0("lag1_", f)]][idx] <- lag1
    }

    # Cumulative features
    if ("c_hat" %in% names(slot_df)) {
      df$cum_max_c_hat[idx] <- cummax(ifelse(is.na(slot_df$c_hat), 0, slot_df$c_hat))
    }
    if ("alpha_entropy" %in% names(slot_df)) {
      df$cum_min_entropy[idx] <- cummin(ifelse(is.na(slot_df$alpha_entropy), Inf, slot_df$alpha_entropy))
    }
    if ("lbf" %in% names(slot_df)) {
      df$cum_max_lbf[idx] <- cummax(ifelse(is.na(slot_df$lbf), 0, slot_df$lbf))
    }
  }

  df
}


#' Compare two diagnostic runs side by side
#'
#' Takes two data.frames (from diagnose_bb_ash_iter or
#' diagnose_ash_filter_archived_iter) at the same iteration
#' and prints a side-by-side comparison.
#'
#' @param df1 First diagnostic data.frame
#' @param df2 Second diagnostic data.frame
#' @param label1 Label for first run (e.g., "BB+ash")
#' @param label2 Label for second run (e.g., "V0")
#' @keywords internal
compare_ash_methods <- function(df1, df2, label1 = "Method1", label2 = "Method2") {
  cat(sprintf("\n===== %s vs %s (iter %d) =====\n", label1, label2,
              df1$iter[1]))
  cat(sprintf("  %-20s %12s %12s\n", "Feature", label1, label2))
  cat(sprintf("  %s\n", strrep("-", 46)))
  cat(sprintf("  %-20s %12d %12d\n", "mask_size", df1$mask_size[1], df2$mask_size[1]))
  cat(sprintf("  %-20s %12d %12d\n", "n_active", df1$n_active[1], df2$n_active[1]))
  cat(sprintf("  %-20s %12.4f %12.4f\n", "sigma2", df1$sigma2[1], df2$sigma2[1]))
  cat(sprintf("  %-20s %12.2e %12.2e\n", "theta_ss", df1$theta_ss[1], df2$theta_ss[1]))
  cat(sprintf("  %-20s %12.2e %12.2e\n", "b_conf_ss", df1$b_conf_ss[1], df2$b_conf_ss[1]))

  # Per-slot comparison for active slots
  active1 <- df1[df1$V > 0 | df1$status != "null", ]
  active2 <- df2[df2$V > 0 | df2$status != "inactive", ]
  n_show <- max(nrow(active1), nrow(active2))

  cat(sprintf("\n  Per-slot:\n"))
  cat(sprintf("  %2s | %-5s %5s %6s %5s | %-5s %5s %6s %5s\n",
              "L", "Sent1", "Pur1", "V1", "Sta1", "Sent2", "Pur2", "V2", "Sta2"))
  cat(sprintf("  %s\n", strrep("-", 60)))

  for (i in seq_len(n_show)) {
    r1 <- if (i <= nrow(active1)) active1[i, ] else NULL
    r2 <- if (i <= nrow(active2)) active2[i, ] else NULL
    l <- if (!is.null(r1)) r1$slot else r2$slot
    cat(sprintf("  %2d |", l))
    if (!is.null(r1))
      cat(sprintf(" %5d %5.3f %6.4f %5s |", r1$sentinel, r1$purity, r1$V,
                  substr(r1$status, 1, 5)))
    else
      cat(sprintf(" %5s %5s %6s %5s |", "-", "-", "-", "-"))
    if (!is.null(r2))
      cat(sprintf(" %5d %5.3f %6.4f %5s", r2$sentinel, r2$purity, r2$V,
                  substr(r2$status, 1, 5)))
    else
      cat(sprintf(" %5s %5s %6s %5s", "-", "-", "-", "-"))
    cat("\n")
  }
  cat(sprintf("==========================================\n"))
}


#' Extract ML feature table from a completed BB+ash fit
#'
#' Computes per-slot features from the converged model. Call with
#' susieR:::extract_bb_ash_features(fit, X_or_Xcorr, causal).
#'
#' @param fit Completed susie fit (with slot_prior + ash)
#' @param X Design matrix (used to compute Xcorr if needed)
#' @param causal Integer vector of true causal indices (for labeling)
#' @return data.frame with one row per slot, all features + TP/FP label
#' @keywords internal
extract_bb_ash_features <- function(fit, X, causal = NULL) {
  L <- nrow(fit$alpha)
  p <- ncol(fit$alpha)
  Xcorr <- cor(X)
  cs_coverage <- 0.9

  c_hat <- if (!is.null(fit$c_hat)) fit$c_hat else
           if (!is.null(fit$slot_weights)) fit$slot_weights else rep(1, L)

  # Map CS to slots
  cs <- fit$sets$cs
  cs_owner <- rep(NA, L)
  cs_is_tp <- rep(NA, L)
  if (length(cs) > 0) {
    for (i in seq_along(cs)) {
      sent <- cs[[i]][which.max(fit$pip[cs[[i]]])]
      owner <- which.max(fit$alpha[, sent])
      cs_owner[owner] <- i
      if (!is.null(causal))
        cs_is_tp[owner] <- any(cs[[i]] %in% causal)
    }
  }

  rows <- list()
  for (l in seq_len(L)) {
    alpha_l <- fit$alpha[l, ]
    sent <- which.max(alpha_l)
    max_a <- max(alpha_l)

    # Purity
    cs_order <- order(alpha_l, decreasing = TRUE)
    cs_size <- min(which(cumsum(alpha_l[cs_order]) >= cs_coverage))
    if (cs_size > 1 && fit$V[l] > 0) {
      cs_idx <- cs_order[1:cs_size]
      pur <- min(abs(Xcorr[cs_idx, cs_idx]))
    } else if (cs_size == 1) {
      pur <- 1.0
    } else {
      pur <- 0.0
    }

    # Alpha entropy
    alpha_nz <- alpha_l[alpha_l > 1e-10]
    alpha_entropy <- -sum(alpha_nz * log(alpha_nz))

    # Max cross-sentinel |r|
    max_cross_r <- 0
    n_colliding <- 0
    active <- which(fit$V > 0)
    if (sent > 0 && length(active) > 1) {
      other_sents <- sapply(setdiff(active, l), function(ll) which.max(fit$alpha[ll, ]))
      other_sents <- other_sents[other_sents > 0]
      if (length(other_sents) > 0) {
        cross_r <- abs(Xcorr[sent, other_sents])
        max_cross_r <- max(cross_r)
        n_colliding <- sum(cross_r > 0.9)
      }
    }

    # Theta at sentinel
    theta_at_sent <- if (!is.null(fit$theta) && sent > 0) fit$theta[sent] else 0

    # CS label
    cs_label <- if (is.na(cs_owner[l])) "-"
                else if (!is.null(causal) && cs_is_tp[l]) "TP"
                else if (!is.null(causal)) "FP"
                else paste0("CS", cs_owner[l])

    rows[[l]] <- data.frame(
      slot = l, sentinel = sent, purity = pur, V = fit$V[l],
      c_hat = c_hat[l],
      lbf = if (!is.null(fit$lbf)) fit$lbf[l] else NA,
      max_alpha = max_a, cs_size = cs_size,
      alpha_entropy = alpha_entropy,
      max_cross_r = max_cross_r, n_colliding = n_colliding,
      theta_at_sent = theta_at_sent,
      n_active = length(active),
      cs_label = cs_label,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
