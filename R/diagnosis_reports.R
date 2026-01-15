#' Diagnostic function for SuSiE-ash iteration
#' @keywords internal
diagnose_susie_ash_iter <- function(data, model, Xcorr, mrash_output,
                                     effect_purity, sentinels, current_collision,
                                     alpha_protected, pip_protected, neighborhood_pip, 
                                     masked, sigma2_new, tau2_new,
                                     want_masked, ready_to_unmask, force_unmask, force_mask,
                                     purity_threshold, cs_formation_threshold,
                                     collision_ld_threshold, tight_ld_threshold,
                                     ld_threshold, diffuse_iter_count_param) {
  
  L <- nrow(model$alpha)
  p <- ncol(model$alpha)
  
  theta_raw <- mrash_output$beta
  theta_masked <- theta_raw
  theta_masked[masked] <- 0
  
  # Identify active effects (non-uniform alpha)
  is_active <- sapply(1:L, function(l) {
    max(model$alpha[l,]) - min(model$alpha[l,]) >= 1e-6
  })
  active_effects <- which(is_active)
  
  # Count effects by CASE
  case_assignment <- rep(NA, L)
  case_reason <- rep("", L)
  for (l in 1:L) {
    if (!is_active[l]) next
    is_ever_diffuse <- model$ever_diffuse[l] > 0
    can_conf <- (effect_purity[l] >= purity_threshold) && !is_ever_diffuse
    
    if (effect_purity[l] < cs_formation_threshold) {
      case_assignment[l] <- 1
      case_reason[l] <- "dif"
    } else if (!can_conf) {
      case_assignment[l] <- 2
      if (current_collision[l]) case_reason[l] <- "col"
      else if (is_ever_diffuse) case_reason[l] <- "edif"
      else case_reason[l] <- "pur"
    } else {
      case_assignment[l] <- 3
      case_reason[l] <- "ok"
    }
  }
  
  n_case1 <- sum(case_assignment == 1, na.rm = TRUE)
  n_case2 <- sum(case_assignment == 2, na.rm = TRUE)
  n_case3 <- sum(case_assignment == 3, na.rm = TRUE)
  n_ever_diffuse <- sum(model$ever_diffuse[active_effects] > 0)
  n_current_collision <- sum(current_collision[active_effects])
  
  cat(sprintf("\n========== SuSiE-ash iter %d ==========\n", model$ash_iter))
  
  # Header
  cat(sprintf("  Effects: %d C3, %d C2, %d C1 | ever_diff: %d, cur_col: %d\n", 
              n_case3, n_case2, n_case1, n_ever_diffuse, n_current_collision))
  cat(sprintf("  Mr.ASH: s2=%.4f, t2=%.2e, pi0=%.1f%% | theta: max=%.4f, SS=%.2e->%.2e\n",
              sigma2_new, tau2_new, mrash_output$pi[1]*100,
              max(abs(theta_raw)), sum(theta_raw^2), sum(theta_masked^2)))
  cat(sprintf("  Thresh: pur=%.2f, col_LD=%.2f, tight_LD=%.2f, mod_LD=%.2f\n",
              purity_threshold, collision_ld_threshold, tight_ld_threshold, ld_threshold))
  
  # Masking summary
  n_prev <- sum(model$masked)
  n_curr <- sum(masked)
  n_new <- sum(masked & !model$masked)
  n_unmask <- sum(ready_to_unmask)
  n_force_unmask <- sum(force_unmask & model$masked)
  n_force_mask <- sum(force_mask & !model$masked)
  cat(sprintf("  Mask: %d->%d (new:%d, unmask:%d) | ever_unmask=%d\n", 
              n_prev, n_curr, n_new, n_unmask, sum(model$ever_unmasked)))
  
  # Unmasking details
  if (n_unmask > 0) {
    unmask_idx <- which(ready_to_unmask)
    n_show <- min(5, length(unmask_idx))
    unmask_info <- sapply(unmask_idx[1:n_show], function(i) {
      type <- if(force_unmask[i]) "F" else "S"
      sprintf("%d[%s]", i, type)
    })
    extra <- if(length(unmask_idx) > 5) sprintf(" +%d", length(unmask_idx)-5) else ""
    cat(sprintf("  Unmasking: %s%s\n", paste(unmask_info, collapse=" "), extra))
  }
  
  # Second chance pending
  pending <- sum(model$force_exposed_iter > 0 & !model$second_chance_used)
  if (pending > 0) cat(sprintf("  2nd-chance pending: %d pos\n", pending))
  
  # Per-Effect Table
  cat(sprintf("\n  --- Effects ---\n"))
  cat(sprintf("  %2s %5s %5s %5s %5s %4s  %s\n", 
              "L", "Sent", "Pur", "LBF", "Case", "Cnt", "Flags"))
  cat(sprintf("  %s\n", strrep("-", 50)))
  
  for (l in 1:L) {
    if (!is_active[l]) next
    
    sentinel_l <- sentinels[l]
    purity_l <- effect_purity[l]
    case_l <- case_assignment[l]
    reason_l <- case_reason[l]
    
    # Counter string
    if (case_l == 2 && !current_collision[l] && model$ever_diffuse[l] == 0) {
      cnt_str <- sprintf("%d/%d", model$diffuse_iter_count[l], diffuse_iter_count_param)
    } else {
      cnt_str <- "-"
    }
    
    # Flags
    flags <- c()
    if (current_collision[l]) {
      # Find which effects this collides with
      col_partners <- c()
      for (other_l in active_effects) {
        if (other_l == l) next
        if (abs(Xcorr[sentinel_l, sentinels[other_l]]) > collision_ld_threshold) {
          col_partners <- c(col_partners, other_l)
        }
      }
      if (length(col_partners) > 0) {
        flags <- c(flags, sprintf("COL(%s)", paste(col_partners, collapse=",")))
      }
    }
    if (model$ever_diffuse[l] > 0) flags <- c(flags, sprintf("ED(%d)", model$ever_diffuse[l]))
    if (masked[sentinel_l]) flags <- c(flags, "M")
    flag_str <- if(length(flags) > 0) paste(flags, collapse=",") else ""
    
    cat(sprintf("  %2d %5d %5.2f %5.1f C%d-%s %4s  %s\n",
                l, sentinel_l, purity_l, model$lbf[l], case_l, reason_l, cnt_str, flag_str))
  }
  
  # Sentinel LD matrix
  if (length(active_effects) > 1 && length(active_effects) <= 10) {
    cat(sprintf("\n  --- Sentinel LD (col_thresh=%.2f) ---\n", collision_ld_threshold))
    
    # Header row
    cat(sprintf("  %8s", ""))
    for (l in active_effects) {
      cat(sprintf(" %6s", sprintf("E%d", l)))
    }
    cat("\n")
    
    # Data rows
    for (i in seq_along(active_effects)) {
      l_i <- active_effects[i]
      cat(sprintf("  E%d:%4d", l_i, sentinels[l_i]))
      for (j in seq_along(active_effects)) {
        l_j <- active_effects[j]
        if (i == j) {
          cat(sprintf(" %6s", "-"))
        } else {
          r_val <- Xcorr[sentinels[l_i], sentinels[l_j]]
          marker <- if (abs(r_val) > collision_ld_threshold) "*" else ""
          cat(sprintf(" %5.2f%s", r_val, marker))
        }
      }
      cat("\n")
    }
  }
  
  # Top theta signals
  cat(sprintf("\n  --- Top |theta| ---\n"))
  top_idx <- order(abs(theta_raw), decreasing = TRUE)[1:min(8, p)]
  
  cat(sprintf("  %5s %8s %2s %5s %5s  %s\n", "Pos", "theta", "St", "PIP", "nPIP", "Info"))
  
  for (idx in top_idx) {
    st <- if(masked[idx]) "M" else "."
    
    # Check if sentinel
    sent_for <- which(sentinels == idx & is_active)
    info <- if(length(sent_for) > 0) sprintf("<-E%d", sent_for[1]) else ""
    
    cat(sprintf("  %5d %+8.4f %2s %5.2f %5.2f  %s\n",
                idx, theta_raw[idx], st,
                pip_protected[idx], neighborhood_pip[idx], info))
  }
  
  # Warning: masked positions with large theta
  masked_with_signal <- which(masked & abs(theta_raw) > 0.01)
  if (length(masked_with_signal) > 0) {
    top_masked <- head(masked_with_signal[order(abs(theta_raw[masked_with_signal]), decreasing=TRUE)], 3)
    info <- sapply(top_masked, function(i) sprintf("%d(%.3f)", i, theta_raw[i]))
    extra <- if(length(masked_with_signal) > 3) sprintf(" +%d", length(masked_with_signal)-3) else ""
    cat(sprintf("\n  WARNING: Masked |theta|>0.01: %s%s\n", paste(info, collapse=", "), extra))
  }
  
  cat(sprintf("==========================================\n"))
}