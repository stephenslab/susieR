#' Diagnostic function for SuSiE-ash iteration
#' @keywords internal
diagnose_susie_ash_iter <- function(data, model, params, Xcorr, mrash_output, 
                                     residuals, b_confident, effect_purity,
                                     pip_protected, neighborhood_pip, masked,
                                     sigma2_new, tau2_new, sa2,
                                     want_masked, ready_to_unmask, force_unmask) {
  
  L <- nrow(model$alpha)
  p <- ncol(model$alpha)
  purity_threshold <- if (!is.null(params$purity_threshold)) params$purity_threshold else 0.5
  cs_formation_threshold <- if (!is.null(params$cs_formation_threshold)) params$cs_formation_threshold else 0.1
  late_emergence_purity_threshold <- if (!is.null(params$late_emergence_purity_threshold)) params$late_emergence_purity_threshold else 0.8
  pip_threshold <- if (!is.null(params$pip_threshold)) params$pip_threshold else 0.5
  direct_pip_threshold <- if (!is.null(params$direct_pip_threshold)) params$direct_pip_threshold else 0.1
  unmask_iter_threshold <- if (!is.null(params$unmask_iter_threshold)) params$unmask_iter_threshold else 3
  
  theta_sum2_before <- sum(mrash_output$beta^2)
  theta_max_before <- max(abs(mrash_output$beta))
  n_high_purity <- sum(effect_purity >= purity_threshold, na.rm = TRUE)
  n_low_purity <- sum(!is.na(effect_purity)) - n_high_purity
  b_full <- colSums(model$alpha * model$mu)
  grid_scale <- data$n / median(colSums(data$X^2))
  pi_null <- mrash_output$pi[1]
  pi_nonnull <- 1 - pi_null
  theta_new <- mrash_output$beta
  theta_new[masked] <- 0
  
  cat(sprintf("\n========== SuSiE-ash iter %d ==========\n", model$ash_iter))
  
  # Purity and effect summary
  cat(sprintf("  Purity: %d high (>=%.2f), %d low | current: %s\n", 
              n_high_purity, purity_threshold, n_low_purity,
              paste(round(effect_purity[!is.na(effect_purity)], 2), collapse=", ")))
  cat(sprintf("  Effect lbf: %s\n", paste(round(model$lbf, 1), collapse=", ")))
  
  # Residual diagnostics
  cat(sprintf("  Residuals: var(y-Xb_conf)=%.4f | b_full²=%.2e, b_conf²=%.2e\n",
              var(as.vector(residuals)), sum(b_full^2), sum(b_confident^2)))
  
  # Mr.ASH output
  cat(sprintf("  Grid: scale=%.3f (n/med(w)) | sa2 range=[%.2e, %.2e]\n",
              grid_scale, min(sa2), max(sa2)))
  cat(sprintf("  Mr.ASH: sigma2=%.4f, tau2=%.2e | pi_null=%.1f%%, pi_nonnull=%.1f%%\n",
              sigma2_new, tau2_new, pi_null*100, pi_nonnull*100))
  cat(sprintf("  Variance for SuSiE: sigma2=%.4f (NOT sigma2_tilde=%.4f)\n",
              sigma2_new, sigma2_new + tau2_new))
  cat(sprintf("  Theta before masking: sum²=%.2e, max|θ|=%.4f | after: sum²=%.2e\n",
              theta_sum2_before, theta_max_before, sum(theta_new^2)))
  
  # =========================================================================
  # Masking/Un-masking diagnostics
  # =========================================================================
  n_previously_masked <- sum(model$masked)
  n_newly_masked <- sum(masked & !model$masked)
  n_unmasked_this_iter <- sum(ready_to_unmask)
  n_force_unmasked <- sum(model$masked & force_unmask)
  n_ever_unmasked <- sum(model$ever_unmasked)
  
  cat(sprintf("\n  --- Masking Summary ---\n"))
  cat(sprintf("  Previous iter: %d masked\n", n_previously_masked))
  cat(sprintf("  This iter: %d total masked (+%d new, -%d unmasked)\n", 
              sum(masked), n_newly_masked, n_unmasked_this_iter))
  if (n_force_unmasked > 0) {
    cat(sprintf("  Force-unmasked (tight LD exposure): %d positions\n", n_force_unmasked))
  }
  cat(sprintf("  Ever un-masked (one-chance used): %d positions\n", n_ever_unmasked))
  
  # Positions being un-masked this iteration
  if (sum(ready_to_unmask) > 0) {
    unmask_idx <- which(ready_to_unmask)
    cat(sprintf("  UN-MASKING %d positions this iter:\n", length(unmask_idx)))
    for (idx in unmask_idx[1:min(10, length(unmask_idx))]) {
      unmask_reason <- if (force_unmask[idx]) "force" else "stable"
      cat(sprintf("    %d: nbhd_pip=%.3f, pip_prot=%.3f, theta=%.4f, iters_waiting=%d [%s]\n",
                  idx, neighborhood_pip[idx], pip_protected[idx], 
                  mrash_output$beta[idx], model$unmask_candidate_iters[idx], unmask_reason))
    }
    if (length(unmask_idx) > 10) cat(sprintf("    ... and %d more\n", length(unmask_idx) - 10))
  }
  
  # Positions close to being un-masked
  close_to_unmask <- model$masked & !model$ever_unmasked & 
                     (model$unmask_candidate_iters > 0) &
                     (model$unmask_candidate_iters < unmask_iter_threshold)
  if (sum(close_to_unmask) > 0) {
    close_idx <- which(close_to_unmask)
    cat(sprintf("  Approaching un-mask (%d positions, need %d iters):\n", 
                length(close_idx), unmask_iter_threshold))
    for (idx in close_idx[1:min(5, length(close_idx))]) {
      cat(sprintf("    %d: iters=%d/%d, nbhd_pip=%.3f, pip_prot=%.3f\n",
                  idx, model$unmask_candidate_iters[idx], unmask_iter_threshold,
                  neighborhood_pip[idx], pip_protected[idx]))
    }
  }
  
  # =========================================================================
  # Top Mr.ASH signals
  # =========================================================================
  cat(sprintf("\n  --- Top Mr.ASH Signals ---\n"))
  top_theta_idx <- order(abs(mrash_output$beta), decreasing = TRUE)[1:10]
  cat(sprintf("  Top 10 |theta| before zeroing:\n"))
  for (k in 1:10) {
    idx <- top_theta_idx[k]
    status <- if (masked[idx]) "MASKED" else "active"
    cat(sprintf("    %d: theta=%.4f [%s], pip_prot=%.3f, nbhd_pip=%.3f\n",
                idx, mrash_output$beta[idx], status, 
                pip_protected[idx], neighborhood_pip[idx]))
  }
  
  # Masked variants with non-trivial Mr.ASH signal (potential signal loss)
  masked_with_signal <- which(masked & abs(mrash_output$beta) > 0.01)
  if (length(masked_with_signal) > 0) {
    cat(sprintf("\n  WARNING: %d masked variants with |theta|>0.01 (signal zeroed):\n", 
                length(masked_with_signal)))
    top_masked <- masked_with_signal[order(abs(mrash_output$beta[masked_with_signal]), decreasing=TRUE)]
    for (k in 1:min(10, length(top_masked))) {
      idx <- top_masked[k]
      why_masked <- if (neighborhood_pip[idx] > pip_threshold) "nbhd_pip" 
                    else if (pip_protected[idx] > direct_pip_threshold) "direct_pip"
                    else "inherited"
      cat(sprintf("    %d: theta=%.4f, nbhd_pip=%.3f, pip_prot=%.3f, reason=%s\n",
                  idx, mrash_output$beta[idx], neighborhood_pip[idx], 
                  pip_protected[idx], why_masked))
    }
  }
  
  # =========================================================================
  # Per-effect details
  # =========================================================================
  cat(sprintf("\n  --- Per-Effect Summary ---\n"))
  for (l in 1:L) {
    max_alpha_l <- max(model$alpha[l,])
    if (max_alpha_l - min(model$alpha[l,]) < 1E-6) next
    sentinel_l <- which.max(model$alpha[l,])
    
    # Determine effect state (matching actual CASE 1/2/3 logic)
    was_ever_diffuse <- model$min_purity[l] < cs_formation_threshold
    effective_threshold <- if (was_ever_diffuse) late_emergence_purity_threshold else purity_threshold
    
    state <- if (effect_purity[l] < cs_formation_threshold) "diffuse"
             else if (effect_purity[l] < effective_threshold) "low_purity"
             else "high_purity"
    
    cat(sprintf("  Effect %d [%s]: sentinel=%d, max_alpha=%.3f, purity=%.2f (min=%.2f), lbf=%.1f\n",
                l, state, sentinel_l, max_alpha_l, effect_purity[l], model$min_purity[l], model$lbf[l]))
    
    if (state == "low_purity") {
      cat(sprintf("           low_purity_count=%d/%d\n", 
                  model$low_purity_iter_count[l],
                  if (!is.null(params$low_purity_iter_count)) params$low_purity_iter_count else 2))
    }
    
    # Nearby theta for this effect
    near_sentinel <- which(abs(Xcorr[sentinel_l,]) > 0.3)
    if (length(near_sentinel) > 0) {
      theta_near <- mrash_output$beta[near_sentinel]
      top_near <- near_sentinel[order(abs(theta_near), decreasing = TRUE)[1:min(3, length(near_sentinel))]]
      near_info <- sapply(top_near, function(j) {
        sprintf("%d(θ=%.3f,r=%.2f,%s)", j, mrash_output$beta[j], 
                Xcorr[sentinel_l, j], if(masked[j]) "M" else ".")
      })
      cat(sprintf("           nearby theta: %s\n", paste(near_info, collapse=", ")))
    }
  }
  
  # =========================================================================
  # Marginal z-scores check
  # =========================================================================
  cat(sprintf("\n  --- Marginal z-scores for top theta ---\n"))
  for (k in 1:min(5, length(top_theta_idx))) {
    idx <- top_theta_idx[k]
    x_j <- data$X[, idx]
    z_resid <- as.numeric(crossprod(x_j, residuals)) / sqrt(sum(x_j^2) * sigma2_new)
    cat(sprintf("    %d: z=%.2f, theta=%.4f, masked=%s\n",
                idx, z_resid, mrash_output$beta[idx], masked[idx]))
  }
  
  cat(sprintf("==========================================\n"))
}
