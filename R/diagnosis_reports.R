#' Diagnostic function for SuSiE-ash iteration
#' @keywords internal
diagnose_susie_ash_iter <- function(data, model, params, Xcorr, mrash_output, 
                                     residuals, b_confident, effect_purity, sentinels,
                                     pip_protected, neighborhood_pip, masked,
                                     sigma2_new, tau2_new, sa2,
                                     want_masked, ready_to_unmask, force_unmask) {
  
  L <- nrow(model$alpha)
  p <- ncol(model$alpha)
  purity_threshold <- if (!is.null(params$purity_threshold)) params$purity_threshold else 0.5
  cs_formation_threshold <- if (!is.null(params$cs_formation_threshold)) params$cs_formation_threshold else 0.1
  low_purity_iter_count_threshold <- if (!is.null(params$low_purity_iter_count)) params$low_purity_iter_count else 2
  sentinel_ld_threshold <- if (!is.null(params$sentinel_ld_threshold)) params$sentinel_ld_threshold else 0.95
  
  theta_new <- mrash_output$beta
  theta_new[masked] <- 0
  
  # Detect current collision for display
  current_collision <- rep(FALSE, L)
  for (l in 1:L) {
    sentinel_l <- sentinels[l]
    other_sentinels <- sentinels[-l]
    if (length(other_sentinels) > 0) {
      current_collision[l] <- any(abs(Xcorr[sentinel_l, other_sentinels]) > sentinel_ld_threshold)
    }
  }
  
  # Count effects by CASE
  n_case1 <- n_case2 <- n_case3 <- 0
  for (l in 1:L) {
    if (max(model$alpha[l,]) - min(model$alpha[l,]) < 1E-6) next
    is_ever_diffuse <- isTRUE(model$ever_diffuse[l])
    can_conf <- (effect_purity[l] >= purity_threshold) && !is_ever_diffuse
    if (effect_purity[l] < cs_formation_threshold) n_case1 <- n_case1 + 1
    else if (!can_conf) n_case2 <- n_case2 + 1
    else n_case3 <- n_case3 + 1
  }
  
  n_ever_diffuse <- sum(sapply(1:L, function(l) isTRUE(model$ever_diffuse[l])))
  n_current_collision <- sum(current_collision)
  
  cat(sprintf("\n========== SuSiE-ash iter %d ==========\n", model$ash_iter))
  
  # === Compact header ===
  cat(sprintf("  Effects: %d CASE3, %d CASE2, %d CASE1 | ever_diffuse: %d, collision: %d\n", 
              n_case3, n_case2, n_case1, n_ever_diffuse, n_current_collision))
  cat(sprintf("  Mr.ASH: σ²=%.4f, τ²=%.2e, π_null=%.1f%% | θ: max=%.4f, Σ²=%.2e→%.2e\n",
              sigma2_new, tau2_new, mrash_output$pi[1]*100,
              max(abs(mrash_output$beta)), sum(mrash_output$beta^2), sum(theta_new^2)))
  
  # === Masking summary ===
  n_prev <- sum(model$masked)
  n_new <- sum(masked & !model$masked)
  n_unmask <- sum(ready_to_unmask)
  n_force <- sum(model$masked & force_unmask)
  cat(sprintf("  Mask: %d→%d (%+d,-%d) | ever_unmask=%d", 
              n_prev, sum(masked), n_new, n_unmask, sum(model$ever_unmasked)))
  if (n_force > 0) cat(sprintf(" | force=%d", n_force))
  cat("\n")
  
  # Unmasking details
  if (n_unmask > 0) {
    unmask_idx <- which(ready_to_unmask)[1:min(5, n_unmask)]
    unmask_info <- sapply(unmask_idx, function(i) {
      sprintf("%d[%s]", i, if(force_unmask[i]) "F" else "S")
    })
    cat(sprintf("  Unmasking: %s%s\n", paste(unmask_info, collapse=" "),
                if(n_unmask > 5) sprintf(" +%d", n_unmask-5) else ""))
  }
  
  # Second chance pending
  pending <- sum(model$force_exposed_iter > 0 & !model$second_chance_used)
  if (pending > 0) cat(sprintf("  2nd-chance pending: %d pos\n", pending))
  
  # === Per-Effect Table ===
  cat(sprintf("\n  --- Effects ---\n"))
  cat(sprintf("  %2s %5s %5s %5s %5s %4s  %s\n", 
              "L", "Sent", "Pur", "LBF", "Case", "Cnt", "Flags"))
  cat(sprintf("  %s\n", paste(rep("-", 45), collapse="")))
  
  for (l in 1:L) {
    max_alpha_l <- max(model$alpha[l,])
    if (max_alpha_l - min(model$alpha[l,]) < 1E-6) next
    
    sentinel_l <- sentinels[l]
    purity_l <- effect_purity[l]
    is_ever_diffuse <- isTRUE(model$ever_diffuse[l])
    has_collision <- current_collision[l]
    can_be_confident <- (purity_l >= purity_threshold) && !is_ever_diffuse
    
    # Determine CASE and reason
    if (purity_l < cs_formation_threshold) {
      case_str <- "C1"
      reason <- "dif"
    } else if (!can_be_confident) {
      case_str <- "C2"
      if (is_ever_diffuse) reason <- "edif"
      else reason <- "pur"
    } else {
      case_str <- "C3"
      reason <- "ok"
    }
    
    # Counter for CASE 2 (only for low purity, not ever_diffuse)
    cnt_str <- if(case_str == "C2" && !is_ever_diffuse) {
      sprintf("%d/%d", model$low_purity_iter_count[l], low_purity_iter_count_threshold)
    } else "-"
    
    # Flags
    flags <- c()
    if (has_collision) flags <- c(flags, "COL")
    if (is_ever_diffuse) flags <- c(flags, "ED")
    if (masked[sentinel_l]) flags <- c(flags, "M")
    flag_str <- paste(flags, collapse=",")
    
    cat(sprintf("  %2d %5d %5.2f %5.1f %2s-%s %4s  %s\n",
                l, sentinel_l, purity_l, model$lbf[l], case_str, reason, cnt_str, flag_str))
  }
  
  # === Top θ signals ===
  cat(sprintf("\n  --- Top |θ| ---\n"))
  top_idx <- order(abs(mrash_output$beta), decreasing = TRUE)[1:8]
  
  header <- sprintf("  %5s %8s %2s %4s %5s  %s", "Pos", "θ", "St", "Pip", "Nbhd", "Info")
  cat(header, "\n")
  
  for (idx in top_idx) {
    st <- if(masked[idx]) "M" else "."
    
    # Is it a sentinel?
    sent_for <- which(sentinels == idx)
    info <- if(length(sent_for) > 0) sprintf("←E%d", sent_for[1]) else ""
    
    cat(sprintf("  %5d %+8.4f %2s %4.2f %5.2f  %s\n",
                idx, mrash_output$beta[idx], st,
                pip_protected[idx], neighborhood_pip[idx], info))
  }
  
  # === Warning: masked with signal ===
  masked_signal <- which(masked & abs(mrash_output$beta) > 0.01)
  if (length(masked_signal) > 0) {
    top_masked <- head(masked_signal[order(abs(mrash_output$beta[masked_signal]), decreasing=TRUE)], 3)
    info <- sapply(top_masked, function(i) sprintf("%d(%.3f)", i, mrash_output$beta[i]))
    cat(sprintf("\n  ⚠ Masked |θ|>0.01: %s", paste(info, collapse=", ")))
    if (length(masked_signal) > 3) cat(sprintf(" +%d", length(masked_signal)-3))
    cat("\n")
  }
  
  cat(sprintf("==========================================\n"))
}