# Benchmark: Multi-panel mixture vs single-panel
#
# Metrics:
#   ELBO:  mixture ELBO vs best single-panel ELBO (should be >=)
#   FDR:   fraction of 95% CS NOT containing a causal variable
#   Power: fraction of causal variables covered by at least one 95% CS
#
# True model: z ~ N(R_mix %*% beta, sigma2 * R_mix + lambda * I)
# where R_mix = w1 * R1 + w2 * R2 (true mixture).
#
# Usage: Rscript inst/notebooks/benchmark_mix_vs_sp.R
# Or:    source("inst/notebooks/benchmark_mix_vs_sp.R") after devtools::load_all()

# Load from the working tree to pick up uncommitted fixes.
if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(susieR)
  message("NOTE: using installed susieR; run from package root for working-tree version")
}

# ---------------------------------------------------------------------------
# Generate X with block-correlated LD structure
# ---------------------------------------------------------------------------
make_block_correlated_X <- function(n, p, rho, block_size = 10) {
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
# CS-based FDR/Power
#   FDR   = fraction of 95% CS that do NOT contain any causal variable
#   Power = fraction of causal variables covered by at least one 95% CS
# ---------------------------------------------------------------------------
cs_fdr_power <- function(fit, causal_idx) {
  cs_list <- fit$sets$cs
  if (is.null(cs_list) || length(cs_list) == 0)
    return(list(fdr = 0, power = 0, n_cs = 0))
  n_cs <- length(cs_list)
  cs_hits <- sapply(cs_list, function(cs) any(cs %in% causal_idx))
  fdr <- sum(!cs_hits) / n_cs
  covered <- sapply(causal_idx, function(j)
    any(sapply(cs_list, function(cs) j %in% cs)))
  power <- mean(covered)
  list(fdr = fdr, power = power, n_cs = n_cs)
}

# ---------------------------------------------------------------------------
# Single scenario runner
# ---------------------------------------------------------------------------
run_scenario <- function(scenario, n_reps = 50) {
  cat(sprintf("=== %s (p=%d, tw=%.1f/%.1f, B=%d/%d, rho=%.1f/%.1f, L=%d, sig=%.1f) ===\n",
              scenario$name, scenario$p,
              scenario$true_w[1], scenario$true_w[2],
              scenario$B1, scenario$B2,
              scenario$rho1, scenario$rho2,
              scenario$n_signals, scenario$signal_strength))

  p <- scenario$p
  n_signals <- scenario$n_signals
  signal_strength <- scenario$signal_strength
  lambda <- scenario$lambda
  max_iter <- scenario$max_iter
  L <- scenario$L

  # Accumulators
  mix_ge_sp     <- 0
  mix_better    <- 0
  safeguard_ct  <- 0
  elbo_diffs    <- numeric(n_reps)
  omega1_est    <- numeric(n_reps)

  # Track per-iteration ELBO decreases in mixture fits
  elbo_decrease_ct <- 0  # reps with at least one ELBO decrease

  # CS-based FDR/power: rows = reps, cols = methods
  methods <- c("sp1", "sp2", "best_sp", "mix")
  fdr_mat   <- matrix(NA, n_reps, 4, dimnames = list(NULL, methods))
  power_mat <- matrix(NA, n_reps, 4, dimnames = list(NULL, methods))
  ncs_mat   <- matrix(NA, n_reps, 4, dimnames = list(NULL, methods))

  for (rep in seq_len(n_reps)) {
    set.seed(scenario$seed_base + rep)

    X1 <- make_block_correlated_X(scenario$B1, p, scenario$rho1)
    X2 <- make_block_correlated_X(scenario$B2, p, scenario$rho2)

    R1 <- crossprod(X1) / scenario$B1
    R2 <- crossprod(X2) / scenario$B2
    R_true <- scenario$true_w[1] * R1 + scenario$true_w[2] * R2

    beta <- rep(0, p)
    causal <- sample(p, n_signals)
    beta[causal] <- signal_strength * sample(c(-1, 1), n_signals, replace = TRUE)
    z <- as.vector(R_true %*% beta) + rnorm(p, sd = sqrt(lambda))

    # Single-panel fits
    fit1 <- susie_rss(z = z, X = X1, lambda = lambda, L = L,
                      max_iter = max_iter,
                      estimate_residual_variance = TRUE, verbose = FALSE)
    fit2 <- susie_rss(z = z, X = X2, lambda = lambda, L = L,
                      max_iter = max_iter,
                      estimate_residual_variance = TRUE, verbose = FALSE)
    best_sp_elbo <- max(tail(fit1$elbo, 1), tail(fit2$elbo, 1))
    best_sp_fit <- if (tail(fit1$elbo, 1) >= tail(fit2$elbo, 1)) fit1 else fit2

    # Mixture fit (with all fixes + safeguard)
    fit_mix <- susie_rss(z = z, X = list(X1, X2), lambda = lambda, L = L,
                         max_iter = max_iter,
                         estimate_residual_variance = TRUE, verbose = FALSE,
                         check_prior = FALSE)
    mix_elbo <- tail(fit_mix$elbo, 1)

    # ELBO tracking
    elbo_diffs[rep] <- mix_elbo - best_sp_elbo
    omega1_est[rep] <- fit_mix$omega_weights[1]

    if (mix_elbo >= best_sp_elbo - 1e-6) mix_ge_sp <- mix_ge_sp + 1
    if (mix_elbo > best_sp_elbo + 0.1)   mix_better <- mix_better + 1
    if (any(fit_mix$omega_weights > 0.999)) safeguard_ct <- safeguard_ct + 1  # omega collapsed to single panel

    # Check for ELBO decreases within the mixture fit
    elbo_traj <- fit_mix$elbo
    if (length(elbo_traj) > 1 && any(diff(elbo_traj) < -1e-6))
      elbo_decrease_ct <- elbo_decrease_ct + 1

    # CS-based FDR / Power
    r1 <- cs_fdr_power(fit1, causal)
    r2 <- cs_fdr_power(fit2, causal)
    r_best <- cs_fdr_power(best_sp_fit, causal)
    r_mix <- cs_fdr_power(fit_mix, causal)

    fdr_mat[rep, ]   <- c(r1$fdr, r2$fdr, r_best$fdr, r_mix$fdr)
    power_mat[rep, ] <- c(r1$power, r2$power, r_best$power, r_mix$power)
    ncs_mat[rep, ]   <- c(r1$n_cs, r2$n_cs, r_best$n_cs, r_mix$n_cs)

    # Print detail for first 2 reps
    if (rep <= 2) {
      cat(sprintf("  Rep %d: SP1=%.1f SP2=%.1f | MIX=%.1f (w=%.2f,%.2f)\n",
                  rep,
                  tail(fit1$elbo, 1), tail(fit2$elbo, 1),
                  mix_elbo, fit_mix$omega_weights[1], fit_mix$omega_weights[2]))
      cat(sprintf("    MIX: %d CS, FDR=%.2f, Power=%.2f | bestSP: %d CS, FDR=%.2f, Power=%.2f\n",
                  r_mix$n_cs, r_mix$fdr, r_mix$power,
                  r_best$n_cs, r_best$fdr, r_best$power))
      cat(sprintf("    ELBO trajectory: %s\n",
                  paste(round(head(fit_mix$elbo, 8), 1), collapse = " -> ")))
    }
  }

  # Per-scenario summary
  cat(sprintf("  ELBO:  Mix>=SP %d/%d | Better %d/%d | Collapsed %d/%d | ELBO-decrease %d/%d\n",
              mix_ge_sp, n_reps, mix_better, n_reps,
              safeguard_ct, n_reps, elbo_decrease_ct, n_reps))
  cat(sprintf("         diff: mean=%.2f min=%.2f max=%.2f\n",
              mean(elbo_diffs), min(elbo_diffs), max(elbo_diffs)))
  cat(sprintf("  CS-FDR:   SP1=%.3f  SP2=%.3f  bestSP=%.3f  MIX=%.3f\n",
              mean(fdr_mat[,"sp1"]), mean(fdr_mat[,"sp2"]),
              mean(fdr_mat[,"best_sp"]), mean(fdr_mat[,"mix"])))
  cat(sprintf("  CS-Power: SP1=%.3f  SP2=%.3f  bestSP=%.3f  MIX=%.3f\n",
              mean(power_mat[,"sp1"]), mean(power_mat[,"sp2"]),
              mean(power_mat[,"best_sp"]), mean(power_mat[,"mix"])))
  cat(sprintf("  Avg #CS:  SP1=%.1f  SP2=%.1f  bestSP=%.1f  MIX=%.1f\n",
              mean(ncs_mat[,"sp1"]), mean(ncs_mat[,"sp2"]),
              mean(ncs_mat[,"best_sp"]), mean(ncs_mat[,"mix"])))
  if (all(scenario$true_w > 0))
    cat(sprintf("  True w1=%.1f, est w1: mean=%.2f sd=%.2f\n",
                scenario$true_w[1], mean(omega1_est), sd(omega1_est)))
  cat("\n")

  invisible(list(
    name = scenario$name,
    mix_ge_sp = mix_ge_sp, mix_better = mix_better,
    safeguard_ct = safeguard_ct, elbo_decrease_ct = elbo_decrease_ct,
    elbo_diffs = elbo_diffs, omega1_est = omega1_est, n_reps = n_reps,
    fdr_mat = fdr_mat, power_mat = power_mat, ncs_mat = ncs_mat
  ))
}

# ---------------------------------------------------------------------------
# Scenario definitions
# ---------------------------------------------------------------------------
scenarios <- list(

  # --- Stress tests (larger p, stronger signals, more LD structure) ---

  list(name = "equal_mix",
       p = 200, B1 = 500, B2 = 500, rho1 = 0.7, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 5, signal_strength = 3.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 10000),

  list(name = "asym_weight",
       p = 200, B1 = 500, B2 = 500, rho1 = 0.7, rho2 = 0.3,
       true_w = c(0.3, 0.7), n_signals = 5, signal_strength = 3.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 20000),

  list(name = "asym_B",
       p = 200, B1 = 200, B2 = 800, rho1 = 0.7, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 5, signal_strength = 3.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 30000),

  list(name = "strong_dense",
       p = 200, B1 = 500, B2 = 500, rho1 = 0.7, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 8, signal_strength = 5.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 40000),

  list(name = "one_correct",
       p = 200, B1 = 500, B2 = 500, rho1 = 0.7, rho2 = 0.3,
       true_w = c(1.0, 0.0), n_signals = 5, signal_strength = 3.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 50000),

  list(name = "extreme_LD",
       p = 200, B1 = 500, B2 = 500, rho1 = 0.9, rho2 = 0.1,
       true_w = c(0.5, 0.5), n_signals = 5, signal_strength = 3.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 60000),

  list(name = "small_B",
       p = 200, B1 = 150, B2 = 150, rho1 = 0.7, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 5, signal_strength = 3.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 70000),

  list(name = "weak_signals",
       p = 200, B1 = 500, B2 = 500, rho1 = 0.7, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 5, signal_strength = 1.5,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 80000),

  # --- Smaller-scale sanity checks (from earlier benchmark) ---

  list(name = "small_equal",
       p = 50, B1 = 200, B2 = 200, rho1 = 0.5, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 3, signal_strength = 2.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 90000),

  list(name = "small_asym_w",
       p = 50, B1 = 200, B2 = 200, rho1 = 0.5, rho2 = 0.3,
       true_w = c(0.3, 0.7), n_signals = 3, signal_strength = 2.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 100000),

  list(name = "small_asym_B",
       p = 50, B1 = 100, B2 = 400, rho1 = 0.5, rho2 = 0.3,
       true_w = c(0.5, 0.5), n_signals = 3, signal_strength = 2.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 110000),

  list(name = "small_one_true",
       p = 50, B1 = 200, B2 = 200, rho1 = 0.5, rho2 = 0.3,
       true_w = c(1.0, 0.0), n_signals = 3, signal_strength = 2.0,
       lambda = 0.1, max_iter = 100, L = 10, seed_base = 120000)
)

# ---------------------------------------------------------------------------
# Run all scenarios
# ---------------------------------------------------------------------------
n_reps <- 50

cat("================================================================\n")
cat("Benchmark: Multi-Panel Mixture vs Single Panel\n")
cat(sprintf("  %d scenarios x %d reps = %d total runs\n",
            length(scenarios), n_reps, length(scenarios) * n_reps))
cat("  FDR/Power based on 95%% credible sets\n")
cat("================================================================\n\n")

results <- lapply(scenarios, run_scenario, n_reps = n_reps)

# ---------------------------------------------------------------------------
# ELBO summary table
# ---------------------------------------------------------------------------
cat("================================================================\n")
cat("ELBO TABLE\n")
cat("================================================================\n\n")
cat(sprintf("%-18s %4s %8s %8s %9s %9s %10s %10s\n",
            "Scenario", "N", "Mix>=SP", "Better", "Safegrd", "ELBOdecr",
            "diff_mean", "diff_min"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (r in results) {
  cat(sprintf("%-18s %4d %5d/%-2d %5d/%-2d %5d/%-2d %6d/%-2d %10.2f %10.2f\n",
              r$name, r$n_reps,
              r$mix_ge_sp, r$n_reps,
              r$mix_better, r$n_reps,
              r$safeguard_ct, r$n_reps,
              r$elbo_decrease_ct, r$n_reps,
              mean(r$elbo_diffs),
              min(r$elbo_diffs)))
}
total_reps <- sum(sapply(results, "[[", "n_reps"))
total_pass <- sum(sapply(results, "[[", "mix_ge_sp"))
total_decr <- sum(sapply(results, "[[", "elbo_decrease_ct"))
cat(sprintf("\nOVERALL: %d/%d mixture >= SP | %d/%d had ELBO decrease\n",
            total_pass, total_reps, total_decr, total_reps))

# ---------------------------------------------------------------------------
# FDR table (95% CS)
# ---------------------------------------------------------------------------
cat("\n\n================================================================\n")
cat("FDR TABLE (95% Credible Sets)\n")
cat("  FDR = fraction of CS not containing any causal variable\n")
cat("================================================================\n\n")
cat(sprintf("%-18s %8s %8s %8s %8s\n",
            "Scenario", "SP1", "SP2", "bestSP", "MIX"))
cat(paste(rep("-", 52), collapse = ""), "\n")
for (r in results) {
  cat(sprintf("%-18s %8.3f %8.3f %8.3f %8.3f\n",
              r$name,
              mean(r$fdr_mat[,"sp1"]),
              mean(r$fdr_mat[,"sp2"]),
              mean(r$fdr_mat[,"best_sp"]),
              mean(r$fdr_mat[,"mix"])))
}

# ---------------------------------------------------------------------------
# Power table (95% CS)
# ---------------------------------------------------------------------------
cat("\n\n================================================================\n")
cat("POWER TABLE (95% Credible Sets)\n")
cat("  Power = fraction of causal vars covered by >= 1 CS\n")
cat("================================================================\n\n")
cat(sprintf("%-18s %8s %8s %8s %8s\n",
            "Scenario", "SP1", "SP2", "bestSP", "MIX"))
cat(paste(rep("-", 52), collapse = ""), "\n")
for (r in results) {
  cat(sprintf("%-18s %8.3f %8.3f %8.3f %8.3f\n",
              r$name,
              mean(r$power_mat[,"sp1"]),
              mean(r$power_mat[,"sp2"]),
              mean(r$power_mat[,"best_sp"]),
              mean(r$power_mat[,"mix"])))
}

# ---------------------------------------------------------------------------
# Average number of CS
# ---------------------------------------------------------------------------
cat("\n\n================================================================\n")
cat("AVG NUMBER OF CS\n")
cat("================================================================\n\n")
cat(sprintf("%-18s %8s %8s %8s %8s\n",
            "Scenario", "SP1", "SP2", "bestSP", "MIX"))
cat(paste(rep("-", 52), collapse = ""), "\n")
for (r in results) {
  cat(sprintf("%-18s %8.1f %8.1f %8.1f %8.1f\n",
              r$name,
              mean(r$ncs_mat[,"sp1"]),
              mean(r$ncs_mat[,"sp2"]),
              mean(r$ncs_mat[,"best_sp"]),
              mean(r$ncs_mat[,"mix"])))
}
