# Shared helpers for the susie_post_outcome_configuration tests
# (test_post_outcome_configuration.R + test_post_outcome_configuration_summary.R).
# testthat auto-sources helper_*.R before running any test file.
#
# Naming map (old -> unified):
#   mk             / (none)            -> coloc_fit()        minimal susie fit
#   mk_mv          / (none)            -> coloc_mv_fit()     minimal mvsusie fit
#   mk_real_fit    / (none)            -> coloc_real_fit()   real susie fit w/ CS
#   mk_tuple       / make_susiex_tuple -> coloc_tuple()      hand-built susiex tuple
#   mk_post_obj    / make_post_obj     -> coloc_post_obj()   tagged result object
#   mk_coloc_row   / (single row)      -> coloc_row()        one coloc-result row
#   (none)         / make_coloc_df     -> coloc_df()         rows-list -> data.frame
#
# All coloc data.frames use the PP.H0..PP.H4 column order throughout.

# Reserved column order for a coloc pairwise result data.frame.
.coloc_cols <- c("trait1", "trait2", "l1", "l2", "hit1", "hit2",
                 "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")

# Minimal hand-built susie fit for precise branch control.
coloc_fit <- function(alpha, lbf, cs = list(L1 = 1L)) {
  structure(list(alpha = alpha, lbf_variable = lbf, sets = list(cs = cs)),
            class = "susie")
}

# Minimal hand-built mvsusie fit with a per-outcome log-BF array.
coloc_mv_fit <- function(alpha, lbf, lbf_outcome, cs = list(L1 = 1L)) {
  structure(list(alpha = alpha, lbf_variable = lbf,
                 lbf_variable_outcome = lbf_outcome,
                 sets = list(cs = cs)),
            class = "mvsusie")
}

# A real susie fit with a populated $sets$cs. signal_sd is large so we get
# genuine credible sets.
coloc_real_fit <- function(seed, n = 200, p = 60, k = 3, signal_sd = 3,
                           L = 5, coverage = 0.95) {
  set.seed(seed)
  dat <- simulate_regression(n = n, p = p, k = k, signal_sd = signal_sd)
  fit <- susie(dat$X, dat$y, L = L, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = coverage)
  fit
}

# A hand-built susiex tuple result (same fields the engine emits).
coloc_tuple <- function(trait_names, cs_indices, marginal_prob,
                        top_config_idx = 1L, prob_thresh = 0.8) {
  N   <- length(trait_names)
  cfg <- as.matrix(expand.grid(rep(list(c(0L, 1L)), N)))
  colnames(cfg) <- paste0("trait_", seq_len(N))
  cp  <- numeric(2L^N)
  cp[top_config_idx] <- 1
  list(
    cs_indices    = setNames(as.integer(cs_indices), trait_names),
    logBF_trait   = setNames(rep(0, N),              trait_names),
    configs       = cfg,
    config_prob   = cp,
    marginal_prob = setNames(marginal_prob,          trait_names),
    active        = setNames(marginal_prob >= prob_thresh, trait_names)
  )
}

# A tagged susie_post_outcome_configuration result wrapping susiex tuples
# and/or a coloc pairwise data.frame.
coloc_post_obj <- function(susiex = NULL, coloc = NULL) {
  out <- list()
  if (!is.null(susiex)) out$susiex <- susiex
  if (!is.null(coloc))  out$coloc_pairwise <- coloc
  class(out) <- c("susie_post_outcome_configuration", "list")
  out
}

# A single coloc-result row. `pp` is c(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4).
coloc_row <- function(t1, t2, pp, l1 = 1L, l2 = 1L, h1 = "r1", h2 = "r2") {
  data.frame(trait1 = t1, trait2 = t2, l1 = l1, l2 = l2,
             hit1 = h1, hit2 = h2,
             PP.H0 = pp[1], PP.H1 = pp[2], PP.H2 = pp[3],
             PP.H3 = pp[4], PP.H4 = pp[5],
             stringsAsFactors = FALSE, row.names = NULL)[, .coloc_cols]
}

# Stack a list of row specs into a coloc data.frame. Each element is a list
# with fields t1, t2, l1, l2, h1, h2 and pp = c(H0, H1, H2, H3, H4).
coloc_df <- function(rows) {
  do.call(rbind, lapply(rows, function(r) {
    coloc_row(r$t1, r$t2, r$pp, l1 = r$l1, l2 = r$l2, h1 = r$h1, h2 = r$h2)
  }))
}
