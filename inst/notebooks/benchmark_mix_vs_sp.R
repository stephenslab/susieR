# Benchmark: Multi-panel mixture vs single-panel ELBO
#
# Verifies that the mixture model (K=2 panels, omega optimized) achieves
# ELBO >= best single-panel fit across a range of scenarios.
#
# True model: z ~ N(R_mix %*% beta, sigma2 * R_mix + lambda * I)
# where R_mix = w1 * R1 + w2 * R2 (true mixture).
#
# Also includes a "broken" mode (stale Rz bug) to demonstrate the fix.
#
# Usage: Rscript inst/notebooks/benchmark_mix_vs_sp.R
# Or:    source("inst/notebooks/benchmark_mix_vs_sp.R") after devtools::load_all()

# Load from the working tree to pick up uncommitted fixes.
# Run from the package root: Rscript inst/notebooks/benchmark_mix_vs_sp.R
# Or source() after devtools::load_all().
if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(susieR)
  message("NOTE: using installed susieR; run from package root for working-tree version")
}

# ---------------------------------------------------------------------------
# Helper: run susie_rss with stale-Rz bug (for "before fix" comparison)
# Monkey-patches update_derived_quantities.rss_lambda to skip Rz recomputation
# ---------------------------------------------------------------------------
run_without_rz_fix <- function(z, X_list, lambda, max_iter, ...) {
  fixed_fn <- susieR:::update_derived_quantities.rss_lambda
  broken_fn <- function(data, params, model) {
    if (!is.null(data$K) && data$K > 1 && !is.null(model$omega)) {
      if (!is.null(data$omega_cache)) {
        model$eigen_R <- susieR:::eigen_from_reduced(
          data$omega_cache, model$omega, data$K, data$p)
      } else if (!is.null(data$panel_R)) {
        R_omega <- Reduce("+", Map("*", model$omega, data$panel_R))
        R_omega <- 0.5 * (R_omega + t(R_omega))
        eig <- eigen(R_omega, symmetric = TRUE)
        eig$values <- pmax(eig$values, 0)
        model$eigen_R <- eig
      }
      model$Vtz          <- crossprod(model$eigen_R$vectors, data$z)
      model$z_null_norm2 <- max(sum(data$z^2) - sum(model$Vtz^2), 0)
      model$X_meta <- susieR:::form_X_meta(data$X_list, model$omega)
      # *** DELIBERATELY SKIP: model$Rz recomputation ***
      if (!is.null(data$stochastic_ld_B))
        model$stochastic_ld_B <- 1 / sum(model$omega^2 / data$B_list)
    }
    eigen_R <- susieR:::get_eigen_R(data, model)
    Dinv <- susieR:::compute_Dinv(model, data)
    V <- eigen_R$vectors; D <- eigen_R$values; Vt <- t(V)
    model$SinvRj   <- V %*% (Dinv * D * Vt)
    model$RjSinvRj <- colSums(Vt * (Dinv * (D^2) * Vt))
    return(model)
  }
  assignInNamespace("update_derived_quantities.rss_lambda", broken_fn, ns = "susieR")
  on.exit(assignInNamespace("update_derived_quantities.rss_lambda", fixed_fn, ns = "susieR"))
  susie_rss(z = z, X = X_list, lambda = lambda, max_iter = max_iter,
            estimate_residual_variance = TRUE, verbose = FALSE,
            multipanel_safeguard = FALSE, ...)
}

# ---------------------------------------------------------------------------
# Generate X with block-correlated LD structure
# ---------------------------------------------------------------------------
make_block_correlated_X <- function(n, p, rho, block_size = 5) {
  X <- matrix(rnorm(n * p), n, p)
  n_blocks <- p %/% block_size
  for (b in seq_len(n_blocks)) {
    idx <- ((b - 1) * block_size + 1):(b * block_size)
    shared <- rnorm(n)
    X[, idx] <- sqrt(1 - rho) * X[, idx] + sqrt(rho) * shared
  }
  X
}

# ---------------------------------------------------------------------------
# Single scenario runner
# ---------------------------------------------------------------------------
run_scenario <- function(scenario, n_reps = 30) {
  cat(sprintf("=== %s (tw=%.1f/%.1f, B=%d/%d) ===\n",
              scenario$name,
              scenario$true_w[1], scenario$true_w[2],
              scenario$B1, scenario$B2))

  p <- scenario$p
  n_signals <- scenario$n_signals
  signal_strength <- scenario$signal_strength
  lambda <- scenario$lambda
  max_iter <- scenario$max_iter

  mix_ge_sp     <- 0
  mix_better    <- 0
  safeguard_ct  <- 0
  broken_worse  <- 0
  elbo_diffs    <- numeric(n_reps)
  broken_diffs  <- numeric(n_reps)
  omega1_est    <- numeric(n_reps)

  for (rep in seq_len(n_reps)) {
    set.seed(scenario$seed_base + rep)

    X1 <- make_block_correlated_X(scenario$B1, p, scenario$rho1)
    X2 <- make_block_correlated_X(scenario$B2, p, scenario$rho2)

    R1 <- crossprod(X1) / scenario$B1
    R2 <- crossprod(X2) / scenario$B2
    R_true <- scenario$true_w[1] * R1 + scenario$true_w[2] * R2

    beta <- rep(0, p)
    if (n_signals > 0) {
      causal <- sample(p, n_signals)
      beta[causal] <- signal_strength * sample(c(-1, 1), n_signals, replace = TRUE)
    }
    z <- as.vector(R_true %*% beta) + rnorm(p, sd = sqrt(lambda))

    # Single-panel fits
    fit1 <- susie_rss(z = z, X = X1, lambda = lambda, max_iter = max_iter,
                      estimate_residual_variance = TRUE, verbose = FALSE)
    fit2 <- susie_rss(z = z, X = X2, lambda = lambda, max_iter = max_iter,
                      estimate_residual_variance = TRUE, verbose = FALSE)
    best_sp <- max(tail(fit1$elbo, 1), tail(fit2$elbo, 1))

    # Mixture fit (with all fixes + safeguard)
    fit_mix <- susie_rss(z = z, X = list(X1, X2), lambda = lambda,
                         max_iter = max_iter,
                         estimate_residual_variance = TRUE, verbose = FALSE,
                         multipanel_safeguard = TRUE)
    mix_elbo <- tail(fit_mix$elbo, 1)

    # Broken fit (stale Rz, no safeguard)
    fit_broken <- tryCatch(
      run_without_rz_fix(z, list(X1, X2), lambda, max_iter),
      error = function(e) NULL)
    broken_elbo <- if (!is.null(fit_broken)) tail(fit_broken$elbo, 1) else NA

    elbo_diffs[rep]  <- mix_elbo - best_sp
    broken_diffs[rep] <- if (!is.na(broken_elbo)) broken_elbo - best_sp else NA
    omega1_est[rep]   <- fit_mix$omega_weights[1]

    if (mix_elbo >= best_sp - 1e-6) mix_ge_sp <- mix_ge_sp + 1
    if (mix_elbo > best_sp + 0.1)   mix_better <- mix_better + 1
    if (!is.na(broken_elbo) && broken_elbo < best_sp - 0.5) broken_worse <- broken_worse + 1

    # Detect safeguard: omega at a vertex means safeguard kicked in
    if (any(fit_mix$omega_weights > 0.999)) safeguard_ct <- safeguard_ct + 1

    # Print ELBO trajectory for first 2 reps (diagnostic detail)
    if (rep <= 2) {
      cat(sprintf("  Rep %d: SP1=%.1f SP2=%.1f | MIX=%.1f (w=%.2f,%.2f) | BROKEN=%.1f\n",
                  rep,
                  tail(fit1$elbo, 1), tail(fit2$elbo, 1),
                  mix_elbo, fit_mix$omega_weights[1], fit_mix$omega_weights[2],
                  if (!is.na(broken_elbo)) broken_elbo else NaN))
      cat(sprintf("    MIX trajectory:    %s\n",
                  paste(round(head(fit_mix$elbo, 8), 1), collapse = " -> ")))
      if (!is.null(fit_broken))
        cat(sprintf("    BROKEN trajectory: %s\n",
                    paste(round(head(fit_broken$elbo, 8), 1), collapse = " -> ")))
    }
  }

  cat(sprintf("  Mix >= SP:         %d/%d (%d%%)\n",
              mix_ge_sp, n_reps, round(100 * mix_ge_sp / n_reps)))
  cat(sprintf("  Mix better (>0.1): %d/%d (%d%%)\n",
              mix_better, n_reps, round(100 * mix_better / n_reps)))
  cat(sprintf("  Safeguard (SP):    %d/%d\n", safeguard_ct, n_reps))
  cat(sprintf("  ELBO diff:         mean=%.2f min=%.2f max=%.2f\n",
              mean(elbo_diffs), min(elbo_diffs), max(elbo_diffs)))

  if (all(scenario$true_w > 0)) {
    cat(sprintf("  True omega1=%.1f, est omega1: mean=%.2f sd=%.2f\n",
                scenario$true_w[1], mean(omega1_est), sd(omega1_est)))
  }

  if (any(!is.na(broken_diffs))) {
    cat(sprintf("  BROKEN worse than SP: %d/%d\n", broken_worse, n_reps))
    cat(sprintf("  BROKEN ELBO diff:  mean=%.2f min=%.2f max=%.2f\n",
                mean(broken_diffs, na.rm = TRUE),
                min(broken_diffs, na.rm = TRUE),
                max(broken_diffs, na.rm = TRUE)))
  }

  cat("\n")
  invisible(list(
    name = scenario$name,
    mix_ge_sp = mix_ge_sp, mix_better = mix_better,
    safeguard_ct = safeguard_ct, broken_worse = broken_worse,
    elbo_diffs = elbo_diffs, broken_diffs = broken_diffs,
    omega1_est = omega1_est, n_reps = n_reps
  ))
}

# ---------------------------------------------------------------------------
# Scenario definitions
# ---------------------------------------------------------------------------
scenarios <- list(
  list(name = "two_wrong_equal",
       p = 50, B1 = 200, B2 = 200, rho1 = 0.5, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 3, signal_strength = 1.0,
       lambda = 0.1, max_iter = 100, seed_base = 10000),

  list(name = "two_wrong_asym_w",
       p = 50, B1 = 200, B2 = 200, rho1 = 0.5, rho2 = 0.3,
       true_w = c(0.3, 0.7), n_signals = 3, signal_strength = 1.0,
       lambda = 0.1, max_iter = 100, seed_base = 20000),

  list(name = "two_wrong_asym_B",
       p = 50, B1 = 100, B2 = 400, rho1 = 0.5, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 3, signal_strength = 1.0,
       lambda = 0.1, max_iter = 100, seed_base = 30000),

  list(name = "strong_signals",
       p = 50, B1 = 200, B2 = 200, rho1 = 0.5, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 3, signal_strength = 2.0,
       lambda = 0.1, max_iter = 100, seed_base = 40000),

  list(name = "many_signals",
       p = 50, B1 = 200, B2 = 200, rho1 = 0.5, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 8, signal_strength = 0.5,
       lambda = 0.1, max_iter = 100, seed_base = 50000),

  list(name = "weak_signals",
       p = 50, B1 = 200, B2 = 200, rho1 = 0.5, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 3, signal_strength = 0.2,
       lambda = 0.1, max_iter = 100, seed_base = 60000),

  list(name = "one_panel_true",
       p = 50, B1 = 200, B2 = 200, rho1 = 0.5, rho2 = 0.3,
       true_w = c(1.0, 0.0), n_signals = 3, signal_strength = 1.0,
       lambda = 0.1, max_iter = 100, seed_base = 70000)
)

# ---------------------------------------------------------------------------
# Run all scenarios
# ---------------------------------------------------------------------------
cat("================================================================\n")
cat("Benchmark: Multi-Panel Mixture vs Single Panel\n")
cat("================================================================\n\n")

results <- lapply(scenarios, run_scenario, n_reps = 30)

cat("================================================================\n")
cat("SUMMARY TABLE\n")
cat("================================================================\n\n")

# Build a data.frame for clean printing
summary_df <- data.frame(
  Scenario    = sapply(results, "[[", "name"),
  N           = sapply(results, "[[", "n_reps"),
  Mix_ge_SP   = sapply(results, "[[", "mix_ge_sp"),
  Mix_better  = sapply(results, "[[", "mix_better"),
  Safeguard   = sapply(results, "[[", "safeguard_ct"),
  FIXED_mean  = round(sapply(results, function(r) mean(r$elbo_diffs)), 2),
  FIXED_min   = round(sapply(results, function(r) min(r$elbo_diffs)), 2),
  BROKEN_mean = round(sapply(results, function(r) mean(r$broken_diffs, na.rm = TRUE)), 2),
  BROKEN_lt_SP = sapply(results, "[[", "broken_worse"),
  stringsAsFactors = FALSE
)

# Print header
cat(sprintf("%-22s %4s %8s %8s %9s %10s %10s %11s %10s\n",
            "Scenario", "N", "Mix>=SP", "Better", "Safeguard",
            "FIXED_mu", "FIXED_min", "BROKEN_mu", "BROKEN<SP"))
cat(paste(rep("-", 105), collapse = ""), "\n")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%-22s %4d %5d/%-2d %5d/%-2d %6d/%-2d %10.2f %10.2f %11.2f %7d/%-2d\n",
              r$Scenario, r$N,
              r$Mix_ge_SP, r$N,
              r$Mix_better, r$N,
              r$Safeguard, r$N,
              r$FIXED_mean, r$FIXED_min,
              r$BROKEN_mean,
              r$BROKEN_lt_SP, r$N))
}

cat("\n")
cat("Column legend:\n")
cat("  Mix>=SP    = reps where mixture ELBO >= best single-panel ELBO\n")
cat("  Better     = reps where mixture ELBO > best SP + 0.1 (meaningfully better)\n")
cat("  Safeguard  = reps where safeguard fell back to single panel\n")
cat("  FIXED_mu   = mean(mixture ELBO - best SP ELBO), with all fixes\n")
cat("  FIXED_min  = min(mixture ELBO - best SP ELBO), worst case\n")
cat("  BROKEN_mu  = mean(mixture ELBO - best SP ELBO), stale Rz bug, no safeguard\n")
cat("  BROKEN<SP  = reps where broken mixture ELBO < best SP - 0.5\n")

# Overall pass/fail
total_reps <- sum(summary_df$N)
total_pass <- sum(summary_df$Mix_ge_SP)
cat(sprintf("\nOVERALL: %d/%d reps mixture >= single panel (%d%%)\n",
            total_pass, total_reps, round(100 * total_pass / total_reps)))
