# Post-hoc causal-configuration probabilities for one or more SuSiE-class fits.
#
# Two algorithms live here, exposed through one entry point:
#
#   * SuSiEx (Nature Genetics, 2024): N-trait 2^N enumeration. Per CS tuple
#     (one CS chosen from each trait), report posterior probabilities over
#     all 2^N "which traits share the causal" patterns plus per-trait
#     marginals. Legacy reference implementation:
#     `mvf.susie.alpha::posthoc_multfsusie`.
#
#   * Coloc pairwise ABF (Wallace, 2020 / `coloc::coloc.bf_bf`): pairwise
#     H0/H1/H2/H3/H4 posteriors for every (trait, trait') pair across every
#     (CS in trait, CS in trait') pair. Implemented inline here as a
#     verbatim port of `coloc:::combine.abf` so susieR has no soft
#     dependency on coloc.
#
#   * mvSuSiE CS-level summary: for one multi-output mvsusie/mfsusie fit,
#     organize native mvSuSiE CS/effect evidence into readable tables.
#
# The public function normalises any supported input shape (single fit, list
# of fits, or a single multi-output fit treated outcome-wise) to a flat list
# of "trait views", then runs the requested algorithms against that list.
# Class-aware branches use `inherits()` and are confined to one helper.
#
# The return value is tagged with class `"susie_post_outcome_configuration"`
# so `summary()` dispatches to the pretty-printer at the bottom of this file.

#' Post-hoc causal-configuration probabilities for one or more SuSiE-class fits
#'
#' Runs one of three complementary post-hoc analyses, selected by
#' \code{method}: \code{"susiex"} (default) for the SuSiEx \eqn{2^N}
#' combinatorial enumeration, reporting posterior probabilities over
#' binary causality patterns across the \eqn{N} input traits;
#' \code{"coloc_pairwise"} for the coloc pairwise ABF, reporting the
#' five colocalisation hypothesis posteriors (H0/H1/H2/H3/H4) for every
#' pair of traits; or \code{"mvsusie"} for CS-level summaries
#' from one multi-output mvSuSiE fit. To get multiple views, call the
#' function more than once and combine.
#'
#' Three input organizations are supported:
#' \describe{
#'   \item{\code{by = "fit"}}{Each input fit contributes a single trait view.
#'     Multi-output fits (\code{mvsusie}, \code{mfsusie}) are kept whole: the
#'     trait's per-(CS, SNP) log Bayes factors are the joint composite
#'     stored on the fit as \code{lbf_variable}. Configuration enumeration
#'     loops over the cross-product \eqn{L_1 \times \dots \times L_N} of CS
#'     indices.}
#'   \item{\code{by = "outcome"}}{Multi-output fits fan out into per-outcome views,
#'     each with its own per-(CS, SNP) log Bayes factors read from
#'     \code{fit$lbf_variable_outcome} (an \eqn{L \times J \times R} or
#'     \eqn{L \times J \times M} array). All per-outcome views share the
#'     joint fit's PIP matrix and CS list, so the configuration enumeration
#'     reduces to a single index \eqn{l \in 1..L}. Single-output \code{susie}
#'     fits pass through unchanged. Requires \code{$lbf_variable_outcome} on the
#'     fit (set \code{attach_lbf_variable_outcome = TRUE} when fitting).}
#'   \item{\code{method = "mvsusie"}}{One multi-output mvSuSiE fit is used
#'     directly, and results are summarized at the CS/effect level.}
#' }
#'
#' After SuSiE has been fit, the post-hoc step can answer two related
#' questions:
#' \describe{
#'   \item{\code{"susiex"}}{SuSiEx \eqn{2^N} activation probabilities for
#'     each credible-set tuple. It asks which traits have evidence for a
#'     proposed shared CS-level event. In two traits, the both-active state
#'     uses same-SNP evidence such as
#'     \eqn{\sum_j \exp(\ell_{1j}+\ell_{2j})}, where \eqn{\ell_{tj}} is the
#'     SNP-level log Bayes factor for trait \eqn{t}.}
#'   \item{\code{"coloc_pairwise"}}{Pairwise coloc H0-H4 posteriors for
#'     every trait pair and CS pair. It asks whether same-SNP evidence is
#'     stronger than distinct-causal evidence, represented by terms such as
#'     \eqn{\sum_{j \ne k}\exp(\ell_{1j}+\ell_{2k})}; this separates H3,
#'     two distinct causal variants, from H4, one shared causal variant.
#'     Per-SNP log Bayes factors are aligned by variant name when available.}
#' }
#'
#' \subsection{Two-trait example}{
#' For two traits and one CS pair, SuSiEx has four activation states:
#' \eqn{(0,0)}, \eqn{(1,0)}, \eqn{(0,1)}, and \eqn{(1,1)}. The state
#' \eqn{(1,1)} means both traits are active for the candidate shared event;
#' it does not introduce a separate distinct-causal state.
#'
#' Coloc has five hypotheses: H0, no causal variant; H1, trait 1 only; H2,
#' trait 2 only; H3, two distinct causal variants; and H4, one shared causal
#' variant. Thus coloc splits the active two-trait case into H3 and H4,
#' whereas SuSiEx represents it as the single activation state \eqn{(1,1)}.
#' 
#' The \code{"susiex"} output is an \eqn{2}-trait post-hoc activation model
#' over CS tuples. The \code{"coloc_pairwise"} output is a pairwise
#' colocalisation model over H0-H4 for each CS pair. Both are summaries of
#' already fitted SuSiE-class models; neither refits the single-trait effects.
#' }
#' 
#' More generally, in \eqn{N}-trait analysis, the difference between
#' SuSiEx-style activation and colocalization is a difference in hypothesis
#' space:
#' \describe{
#'   \item{SuSiEx activation space}{\code{"susiex"} is an \eqn{N}-trait
#'     post-hoc activation, or meta-analysis-style, model over CS tuples. It
#'     enumerates the \eqn{2^N} binary activity patterns and asks which traits
#'     participate in a proposed shared event.}
#'   \item{Colocalization partition space}{A theoretical generalization to
#'     \eqn{N}-trait colocalization has a larger hypothesis space: inactive
#'     traits plus all ways of grouping active traits by shared causal
#'     variant, of size Bell(\eqn{N+1}) (see HyPrColoc, Figure 1). This space
#'     asks not only which traits are active, but which active traits share
#'     the same causal variant versus distinct causal variants.}
#' }
#' 
#' @param input A single fit of class \code{susie}, \code{mvsusie}, or
#'   \code{mfsusie}, OR a list of such fits.
#' @param by Either \code{"fit"} (one trait per input fit; default) or
#'   \code{"outcome"} (multi-output fits expand into per-outcome traits).
#' @param outcome_names Optional character vector of trait-view names. If
#'   \code{NULL}, names are taken from the input list when available, otherwise
#'   \code{trait_1}, \code{trait_2}, ... are used. If provided, its length must
#'   match the number of trait views after applying \code{by}.
#' @param method Character scalar; one of \code{"susiex"} (default),
#'   \code{"coloc_pairwise"}, or \code{"mvsusie"}. Pick the analysis to run;
#'   for multiple analyses, call the function more than once.
#' @param prob_thresh Per-trait marginal threshold for the convenience
#'   \code{$active} flags in the SuSiEx output. Default \code{0.8}.
#' @param single_effect_lfsr_cutoff LFSR threshold used to flag mvSuSiE
#'   outcome-specific CS evidence in the \code{"mvsusie"} summary.
#'   The LFSR is read from \code{fit$lfsr} at each CS sentinel SNP.
#'   Default \code{0.05}.
#' @param cs_only Logical. If \code{TRUE} (default) only enumerate over CSs
#'   present in each fit's \code{$sets$cs}; if \code{FALSE} loop over all L
#'   rows of \code{$alpha}. Either way, effects whose entire alpha row is
#'   zero are skipped. When \code{TRUE}, trait views without credible sets
#'   are skipped with a warning; if fewer than two trait views remain, the
#'   function returns \code{NULL}.
#' @param p1,p2,p12 Coloc per-SNP causal priors: \code{p1} for trait 1
#'   alone, \code{p2} for trait 2 alone, \code{p12} for shared causal.
#'   Defaults match \code{coloc::coloc.bf_bf}: \code{p1 = p2 = 1e-4},
#'   \code{p12 = 5e-6}. Only used when \code{"coloc_pairwise"} is in
#'   \code{method}.
#' @param ... Currently ignored.
#'
#' @return \code{NULL} if fewer than two trait views are available for
#' post-hoc analysis; otherwise a list of class
#' \code{"susie_post_outcome_configuration"} with exactly one of the following
#' components, depending on \code{method}:
#' \describe{
#'   \item{\code{$susiex}}{(when \code{method = "susiex"}) A list of length
#'     equal to the number of CS tuples considered. Each element has
#'     components \code{config_probability} (binary trait-activation table
#'     with a \code{config_prob} column) and \code{activation_summary}
#'     (CS-level post-hoc activation summary by trait).}
#'   \item{\code{$coloc_pairwise}}{(when \code{method = "coloc_pairwise"})
#'     A data.frame with one row per (trait1, trait2, l1, l2)
#'     combination, columns \code{trait1, trait2, l1, l2, hit1, hit2,
#'     PP.H0, PP.H1, PP.H2, PP.H3, PP.H4}.}
#'   \item{\code{$mvsusie}}{(when \code{method = "mvsusie"}) Summarizes one
#'     multi-output \code{mvsusie} or \code{mfsusie} fit by CS, returning
#'     \code{cs_summary}, \code{config_summary} (outcome log BF, sentinel
#'     LFSR, cutoff flag), and \code{cs_variant_summary} (CS variant PIP and
#'     outcome-specific LFSR).}
#' }
#' Pretty-print with \code{summary(out)}.
#'
#' @references
#' SuSiEx, Nature Genetics 2024 (combinatorial \eqn{2^N} enumeration).
#' Wallace, PLoS Genetics 2020 (coloc pairwise H0/H1/H2/H3/H4 ABF).
#' Foley et al., Nature Communications 2021 (HyPrColoc; multi-trait
#' colocalisation hypothesis space in Figure 1).
#'
#' @export
susie_post_outcome_configuration <- function(input,
                                             by          = c("fit", "outcome"),
                                             outcome_names = NULL,
                                             method      = c("susiex",
                                             "coloc_pairwise",
                                             "mvsusie"),
                                             prob_thresh = 0.8,
                                             single_effect_lfsr_cutoff = 0.05,
                                             cs_only     = TRUE,
                                             p1          = 1e-4,
                                             p2          = 1e-4,
                                             p12         = 5e-6,
                                             ...) {
  by     <- match.arg(by)
  method <- match.arg(method)

  if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L ||
      !is.finite(prob_thresh) || prob_thresh < 0 || prob_thresh > 1) {
    stop("`prob_thresh` must be a single numeric in [0, 1].")
  }
  if (!is.numeric(single_effect_lfsr_cutoff) ||
      length(single_effect_lfsr_cutoff) != 1L ||
      !is.finite(single_effect_lfsr_cutoff) ||
      single_effect_lfsr_cutoff < 0 || single_effect_lfsr_cutoff > 1) {
    stop("`single_effect_lfsr_cutoff` must be a single numeric in [0, 1].")
  }
  if (!is.logical(cs_only) || length(cs_only) != 1L || is.na(cs_only)) {
    stop("`cs_only` must be a single logical (TRUE or FALSE).")
  }
  for (nm in c("p1", "p2", "p12")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
        v <= 0 || v >= 1) {
      stop("`", nm, "` must be a single numeric in (0, 1).")
    }
  }

  out <- list()
  if (identical(method, "mvsusie")) {
    out$mvsusie <- mvsusie_cs_summary(input,
                                      outcome_names = outcome_names,
                                      single_effect_lfsr_cutoff =
                                        single_effect_lfsr_cutoff,
                                      cs_only       = cs_only)
    if (is.null(out$mvsusie)) return(NULL)
  } else {
    views <- normalise_to_views(input, by = by)
    if (!is.null(outcome_names)) {
      if (!is.character(outcome_names) ||
          length(outcome_names) != length(views) ||
          anyNA(outcome_names) || any(!nzchar(outcome_names))) {
        stop("`outcome_names` must be a non-empty character vector with ",
             "length equal to the number of trait views.")
      }
      for (k in seq_along(views)) {
        views[[k]]$name <- outcome_names[k]
      }
    }
    views <- filter_views_for_posthoc(views, cs_only = cs_only)
    if (is.null(views)) return(NULL)

    if (identical(method, "susiex")) {
      out$susiex <- susiex_configurations(views,
                                          by          = by,
                                          prob_thresh = prob_thresh)
    } else {
      # method == "coloc_pairwise"
      out$coloc_pairwise <- coloc_pairwise_abf(views,
                                               p1  = p1,
                                               p2  = p2,
                                               p12 = p12)
    }
  }
  attr(out, "prob_thresh") <- prob_thresh
  attr(out, "single_effect_lfsr_cutoff") <- single_effect_lfsr_cutoff
  attr(out, "method")      <- method
  class(out) <- c("susie_post_outcome_configuration", "list")
  out
}

# -----------------------------------------------------------------------------
# Input normalisation
# -----------------------------------------------------------------------------

is_susie_fit <- function(x) {
  inherits(x, "susie") || inherits(x, "mvsusie") || inherits(x, "mfsusie")
}

normalise_to_views <- function(input, by) {
  fits <- if (is_susie_fit(input)) list(input) else as.list(input)

  if (length(fits) == 0L) {
    stop("`input` must be a SuSiE-class fit or a non-empty list of fits.")
  }
  for (k in seq_along(fits)) {
    if (!is_susie_fit(fits[[k]])) {
      stop("Element ", k,
           " of `input` is not a SuSiE-class fit (`susie`, `mvsusie`, or ",
           "`mfsusie`).")
    }
  }

  raw_names <- names(fits)
  if (is.null(raw_names)) raw_names <- character(length(fits))
  default_names <- paste0("trait_", seq_along(fits))
  raw_names[!nzchar(raw_names)] <- default_names[!nzchar(raw_names)]

  views <- vector("list", 0)
  for (k in seq_along(fits)) {
    views <- c(views, expand_one_fit(fits[[k]], raw_names[k], by = by))
  }
  views
}

filter_views_for_posthoc <- function(views, cs_only) {
  if (cs_only) {
    has_cs <- lengths(lapply(views, view_cs_indices, cs_only = TRUE)) > 0L
    if (any(!has_cs)) {
      warning("Skipping trait view(s) with no credible sets for post-hoc ",
              "analysis: ", paste(vapply(views[!has_cs], function(v) v$name,
                                          character(1)), collapse = ", "),
              call. = FALSE)
    }
    views <- views[has_cs]
  }

  if (length(views) < 2L) {
    message("Fewer than two trait views ",
            if (cs_only) "with credible sets " else "",
            "for post-hoc outcome configuration analysis; returning NULL.")
    return(NULL)
  }
  views
}

expand_one_fit <- function(fit, base_name, by) {
  if (by == "fit") {
    return(list(make_view(
      name    = base_name,
      alpha   = fit$alpha,
      lbf     = fit$lbf_variable,
      sets_cs = fit$sets$cs
    )))
  }

  # by = "outcome": multi-output fits fan out; single-output fits pass
  # through as one view.
  if (inherits(fit, "mvsusie") || inherits(fit, "mfsusie")) {
    if (is.null(fit$lbf_variable_outcome)) {
      stop("Fit '", base_name, "': `by = \"outcome\"` requires `$lbf_variable_outcome` ",
           "(an L x J x R or L x J x M array) on the fit. ",
           "Refit with `attach_lbf_variable_outcome = TRUE` (the default in mfsusie / ",
           "mvsusie), or pass `by = \"fit\"` to use the joint composite log ",
           "BF instead.")
    }
    R <- dim(fit$lbf_variable_outcome)[3L]
    out_names <- dimnames(fit$lbf_variable_outcome)[[3L]]
    if (is.null(out_names)) out_names <- paste0("outcome_", seq_len(R))
    views <- vector("list", R)
    for (r in seq_len(R)) {
      views[[r]] <- make_view(
        name    = paste0(base_name, "_", out_names[r]),
        alpha   = fit$alpha,
        lbf     = fit$lbf_variable_outcome[, , r, drop = TRUE],
        sets_cs = fit$sets$cs
      )
    }
    return(views)
  }

  # Single-output `susie` under by = "outcome": same as by = "fit".
  list(make_view(
    name    = base_name,
    alpha   = fit$alpha,
    lbf     = fit$lbf_variable,
    sets_cs = fit$sets$cs
  ))
}

make_view <- function(name, alpha, lbf, sets_cs) {
  if (is.null(alpha) || is.null(lbf)) {
    stop("Trait '", name, "': both `$alpha` and `$lbf_variable` (or per-",
         "outcome lbf row) must be non-null.")
  }
  if (!is.matrix(alpha)) alpha <- as.matrix(alpha)
  if (!is.matrix(lbf))   lbf   <- as.matrix(lbf)
  if (!identical(dim(alpha), dim(lbf))) {
    stop("Trait '", name, "': `alpha` and `lbf` must have identical shape; ",
         "got ", paste(dim(alpha), collapse = "x"), " vs ",
         paste(dim(lbf), collapse = "x"), ".")
  }
  list(name = name, alpha = alpha, lbf = lbf, sets_cs = sets_cs)
}

# -----------------------------------------------------------------------------
# CS-tuple enumeration shared by both algorithms.
# -----------------------------------------------------------------------------

# Per-view CS index set, restricted to $sets$cs when cs_only = TRUE.
view_cs_indices <- function(view, cs_only) {
  L_n <- nrow(view$alpha)
  if (!cs_only) return(seq_len(L_n))
  if (is.null(view$sets_cs) ||
      (length(view$sets_cs) == 0L && is.null(attr(view$sets_cs, "cs_idx")))) {
    return(integer(0))
  }

  idx <- attr(view$sets_cs, "cs_idx")
  if (is.null(idx)) {
    # Fall back to the names of $sets$cs ("L1", "L2", ... in susieR's format).
    if (length(view$sets_cs) > 0L && !is.null(names(view$sets_cs))) {
      idx <- as.integer(sub("^L", "", names(view$sets_cs)))
    } else {
      idx <- seq_len(L_n)
    }
  }
  idx[idx >= 1L & idx <= L_n]
}

view_cs_label <- function(view, idx) {
  target <- paste0("L", idx)
  cs <- view$sets_cs
  if (is.null(cs) || is.na(idx)) return(target)

  cs_names <- names(cs)
  if (is.null(cs_names)) return(target)

  cs_idx <- attr(cs, "cs_idx")
  if (!is.null(cs_idx)) {
    pos <- match(idx, cs_idx)
    if (!is.na(pos) && pos <= length(cs_names) && nzchar(cs_names[pos])) {
      return(cs_names[pos])
    }
  }

  target
}

# Returns a list of length-N integer tuples (one CS index per view).
# Under by = "outcome" all views share CSs and we use the diagonal.
# Under by = "fit" we use the cross-product.
enumerate_cs_tuples <- function(views, by, cs_only) {
  per_view <- lapply(views, view_cs_indices, cs_only = cs_only)
  if (any(vapply(per_view, length, integer(1)) == 0L)) return(list())

  if (by == "outcome") {
    common <- Reduce(intersect, per_view)
    lapply(common, function(l) rep.int(l, length(views)))
  } else {
    grid <- expand.grid(per_view, KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(grid)), function(i) as.integer(grid[i, , drop = TRUE]))
  }
}

# -----------------------------------------------------------------------------
# SuSiEx 2^N configuration enumeration.
# -----------------------------------------------------------------------------

.view_variant_keys <- function(view) {
  alpha_names <- colnames(view$alpha)
  lbf_names   <- colnames(view$lbf)
  if (!is.null(alpha_names) && !is.null(lbf_names) &&
      !identical(alpha_names, lbf_names)) {
    stop("Trait '", view$name, "': column names of `alpha` and `lbf` ",
         "must match for variant-level post-hoc scoring.")
  }
  keys <- if (!is.null(lbf_names)) lbf_names else alpha_names
  if (is.null(keys)) keys <- paste0(".variant_", seq_len(ncol(view$lbf)))
  keys
}

.susiex_config_logbf <- function(config, tuple, views, variant_keys,
                                 variant_space_size) {
  active <- which(config == 1L)
  if (length(active) == 0L) return(0)

  common <- Reduce(intersect, variant_keys[active])
  if (length(common) == 0L) return(-Inf)

  logbf <- rep(-log(variant_space_size), length(common))
  for (n in active) {
    idx <- match(common, variant_keys[[n]])
    logbf <- logbf + views[[n]]$lbf[tuple[n], idx]
  }
  .logsum(logbf)
}

susiex_configurations <- function(views, by, prob_thresh,
                                  max_traits = 20L) {
  N <- length(views)
  if (N > max_traits) {
    stop("susiex: N = ", N, " traits exceeds the safety ceiling (",
         max_traits, "); 2^N enumeration would be too large.")
  }

  cs_tuples <- enumerate_cs_tuples(views, by = by, cs_only = TRUE)
  if (length(cs_tuples) == 0L) return(list())

  configs <- as.matrix(expand.grid(rep(list(c(0L, 1L)), N)))
  colnames(configs) <- paste0("trait_", seq_len(N))
  trait_names <- vapply(views, function(v) v$name, character(1))
  variant_keys <- lapply(views, .view_variant_keys)
  variant_space_size <- length(Reduce(union, variant_keys))

  out <- vector("list", length(cs_tuples))
  for (ti in seq_along(cs_tuples)) {
    tuple <- cs_tuples[[ti]]

    logBF_trait <- numeric(N)
    skip <- FALSE
    for (n in seq_len(N)) {
      l_n       <- tuple[n]
      alpha_row <- views[[n]]$alpha[l_n, ]
      if (all(alpha_row == 0)) { skip <- TRUE; break }
      singleton <- integer(N)
      singleton[n] <- 1L
      logBF_trait[n] <- .susiex_config_logbf(
        config = singleton,
        tuple = tuple,
        views = views,
        variant_keys = variant_keys,
        variant_space_size = variant_space_size
      )
    }
    if (skip) {
      out[[ti]] <- NULL
      next
    }

    logBF_conf    <- vapply(seq_len(nrow(configs)), function(m) {
      .susiex_config_logbf(
        config = configs[m, ],
        tuple = tuple,
        views = views,
        variant_keys = variant_keys,
        variant_space_size = variant_space_size
      )
    }, numeric(1))
    maxlog        <- max(logBF_conf)
    prob_conf     <- exp(logBF_conf - maxlog)
    prob_conf     <- prob_conf / sum(prob_conf)
    marginal_prob <- as.vector(crossprod(configs, prob_conf))
    names(logBF_trait) <- names(marginal_prob) <- trait_names
    active <- setNames(marginal_prob >= prob_thresh, trait_names)

    out[[ti]] <- list(
      cs_indices    = setNames(as.integer(tuple), trait_names),
      cs_labels     = setNames(vapply(seq_len(N), function(n) {
        view_cs_label(views[[n]], tuple[n])
      }, character(1)), trait_names),
      logBF_trait   = logBF_trait,
      configs       = configs,
      config_prob   = prob_conf,
      marginal_prob = marginal_prob,
      active        = active
    )
  }

  .organize_susiex_output(out[!vapply(out, is.null, logical(1))])
}

.organize_susiex_output <- function(susiex) {
  lapply(susiex, function(tup) {
    if (!is.list(tup)) return(tup)

    trait_names <- names(tup$cs_indices)
    if (is.null(trait_names) || any(!nzchar(trait_names))) {
      trait_names <- names(tup$marginal_prob)
    }

    config_probability <- if (is.matrix(tup$configs) &&
                              !is.null(tup$config_prob) &&
                              nrow(tup$configs) == length(tup$config_prob)) {
      config_probability <- as.data.frame(tup$configs)
      if (!is.null(trait_names) && length(trait_names) == ncol(tup$configs)) {
        colnames(config_probability) <- trait_names
      }
      config_probability$config_prob <- tup$config_prob
      config_probability
    } else NULL

    activation_summary <- if (!is.null(trait_names) &&
                              length(trait_names) > 0L) {
      cs_display <- if (!is.null(tup$cs_labels)) {
        tup$cs_labels[trait_names]
      } else {
        paste0("L", tup$cs_indices[trait_names])
      }
      vals <- rbind(
        cs_indices    = as.character(cs_display),
        logBF_trait   = as.character(tup$logBF_trait[trait_names]),
        posthoc_prob  = as.character(tup$marginal_prob[trait_names]),
        active        = as.character(tup$active[trait_names])
      )
      colnames(vals) <- trait_names
      as.data.frame(vals, check.names = FALSE, stringsAsFactors = FALSE)
    } else NULL

    out <- list(config_probability = config_probability,
                activation_summary = activation_summary)
    attr(out, "raw") <- tup
    out
  })
}


# -----------------------------------------------------------------------------
# mvSuSiE CS-level summary helpers.
# -----------------------------------------------------------------------------

as_single_mvsusie_fit <- function(input) {
  is_mv <- function(x) inherits(x, "mvsusie") || inherits(x, "mfsusie")
  if (is_mv(input)) return(input)
  if (is.list(input) && length(input) == 1L && is_mv(input[[1L]])) {
    return(input[[1L]])
  }
  stop("method = \"mvsusie\" requires one multi-output fit of class ",
       "`mvsusie` or `mfsusie`.")
}

mvsusie_names <- function(candidates, n, prefix, provided = NULL,
                          what = "names") {
  if (!is.null(provided)) {
    if (!is.character(provided) || length(provided) != n ||
        anyNA(provided) || any(!nzchar(provided))) {
      stop("`", what, "` must be a non-empty character vector of length ", n, ".")
    }
    return(provided)
  }
  for (x in candidates) {
    if (!is.null(x) && length(x) == n && all(!is.na(x)) && all(nzchar(x))) {
      return(as.character(x))
    }
  }
  paste0(prefix, seq_len(n))
}

mvsusie_lbf_outcome <- function(fit) {
  if (!is.null(fit$lbf_outcome)) {
    x <- fit$lbf_outcome
    if (!is.matrix(x)) x <- as.matrix(x)
    return(x)
  }

  if (is.null(fit$alpha) || is.null(fit$lbf_variable_outcome)) {
    stop("method = \"mvsusie\" requires `$lbf_outcome`, or both `$alpha` ",
         "and `$lbf_variable_outcome` to compute CS-level log Bayes factors.")
  }
  alpha <- fit$alpha
  arr <- fit$lbf_variable_outcome
  if (!is.matrix(alpha) || length(dim(arr)) != 3L ||
      nrow(alpha) != dim(arr)[1L] || ncol(alpha) != dim(arr)[2L]) {
    stop("For method = \"mvsusie\", `$alpha` must be L x J and ",
         "`$lbf_variable_outcome` must be L x J x R with matching L and J.")
  }

  out <- matrix(0, nrow = dim(arr)[1L], ncol = dim(arr)[3L])
  for (l in seq_len(nrow(out))) {
    out[l, ] <- as.vector(crossprod(alpha[l, ], arr[l, , ]))
  }
  colnames(out) <- dimnames(arr)[[3L]]
  out
}

mvsusie_cs_indices <- function(fit, cs_only) {
  view <- list(alpha = fit$alpha, sets_cs = fit$sets$cs)
  idx <- view_cs_indices(view, cs_only = cs_only)
  list(idx = idx,
       labels = vapply(idx, function(l) view_cs_label(view, l), character(1)))
}

mvsusie_cs_hit <- function(fit, label, idx, variable_names) {
  cs <- fit$sets$cs
  empty <- list(n_cs = 0L, hit_idx = NA_integer_, hit = NA_character_,
                maxPIP = NA_real_, variants = character(0),
                variant_idx = integer(0), pips = numeric(0))
  if (is.null(cs) || length(cs) == 0L) return(empty)

  pos <- if (!is.null(names(cs))) match(label, names(cs)) else NA_integer_
  if (is.na(pos)) {
    cs_idx <- attr(cs, "cs_idx")
    if (!is.null(cs_idx)) pos <- match(idx, cs_idx)
  }
  if (is.na(pos) || pos < 1L || pos > length(cs)) return(empty)

  members <- cs[[pos]]
  member_idx <- if (is.character(members)) {
    match(members, variable_names)
  } else {
    as.integer(members)
  }
  member_idx <- member_idx[!is.na(member_idx) &
                             member_idx >= 1L &
                             member_idx <= length(variable_names)]
  variants <- variable_names[member_idx]

  if (length(members) == 0L || is.null(fit$pip)) {
    empty$n_cs <- length(members)
    empty$variants <- variants
    empty$variant_idx <- member_idx
    empty$pips <- rep(NA_real_, length(variants))
    return(empty)
  }
  member_idx <- member_idx[member_idx <= length(fit$pip)]
  if (length(member_idx) == 0L) {
    empty$n_cs <- length(members)
    empty$variants <- variants
    empty$variant_idx <- member_idx
    empty$pips <- rep(NA_real_, length(variants))
    return(empty)
  }

  pip <- fit$pip[member_idx]
  best <- member_idx[which.max(pip)]
  list(n_cs = length(members), hit_idx = best, hit = variable_names[best],
       maxPIP = unname(fit$pip[best]), variants = variants,
       variant_idx = member_idx, pips = pip)
}

mvsusie_cs_summary <- function(input, outcome_names = NULL,
                               single_effect_lfsr_cutoff = 0.05,
                               cs_only) {
  fit <- as_single_mvsusie_fit(input)
  logBF_mat <- mvsusie_lbf_outcome(fit)
  if (nrow(logBF_mat) != nrow(fit$alpha)) {
    stop("For method = \"mvsusie\", `$lbf_outcome` must have one row per ",
         "single effect in `$alpha`.")
  }

  R <- ncol(logBF_mat)
  if (R < 2L) {
    message("Fewer than two outcomes in mvSuSiE fit; returning NULL.")
    return(NULL)
  }
  trait_names <- mvsusie_names(
    list(colnames(fit$lbf_outcome), colnames(fit$lfsr),
         fit$outcome_names,
         if (!is.null(fit$lbf_variable_outcome)) {
           dimnames(fit$lbf_variable_outcome)[[3L]]
         } else NULL),
    n = R, prefix = "outcome_", provided = outcome_names,
    what = "outcome_names"
  )

  cs <- mvsusie_cs_indices(fit, cs_only = cs_only)
  if (length(cs$idx) == 0L) {
    message("No credible sets available for mvSuSiE CS summary; ",
            "returning NULL.")
    return(NULL)
  }

  sentinel_lfsr_mat <- fit$lfsr
  if (!is.null(sentinel_lfsr_mat) && !is.matrix(sentinel_lfsr_mat)) {
    sentinel_lfsr_mat <- as.matrix(sentinel_lfsr_mat)
  }
  effect_lbf <- rep(NA_real_, nrow(fit$alpha))
  if (!is.null(fit$lbf)) {
    lbf <- as.numeric(unlist(fit$lbf, use.names = FALSE))
    effect_lbf[seq_len(min(length(lbf), length(effect_lbf)))] <- lbf
  }
  variable_names <- mvsusie_names(
    list(fit$variable_names, names(fit$pip), colnames(fit$alpha),
         colnames(fit$lbf_variable)),
    n = ncol(fit$alpha), prefix = "", what = "variable names"
  )

  out <- vector("list", length(cs$idx))
  for (i in seq_along(cs$idx)) {
    l <- cs$idx[i]
    if (all(fit$alpha[l, ] == 0)) next

    hit <- mvsusie_cs_hit(fit, label = cs$labels[i], idx = l,
                          variable_names = variable_names)
    sentinel_lfsr <- rep(NA_real_, R)
    if (!is.null(sentinel_lfsr_mat) &&
        !is.na(hit$hit_idx) &&
        nrow(sentinel_lfsr_mat) >= hit$hit_idx &&
        ncol(sentinel_lfsr_mat) == R) {
      sentinel_lfsr <- as.numeric(sentinel_lfsr_mat[hit$hit_idx, ])
    }
    sentinel_lfsr_pass <- sentinel_lfsr < single_effect_lfsr_cutoff

    cs_summary <- data.frame(
      cs            = cs$labels[i],
      single_effect = l,
      n_cs          = hit$n_cs,
      hit           = hit$hit,
      maxPIP        = hit$maxPIP,
      lbf           = effect_lbf[l],
      n_lfsr_pass   = sum(sentinel_lfsr_pass, na.rm = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    cs_variant_summary <- data.frame(
      variant = hit$variants,
      pip = unname(hit$pips),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (!is.null(sentinel_lfsr_mat) &&
        length(hit$variant_idx) > 0L &&
        length(hit$variant_idx) == length(hit$variants) &&
        nrow(sentinel_lfsr_mat) >= max(hit$variant_idx, na.rm = TRUE) &&
        ncol(sentinel_lfsr_mat) == R) {
      cs_lfsr <- sentinel_lfsr_mat[hit$variant_idx, , drop = FALSE]
      for (r in seq_len(R)) {
        cs_variant_summary[[paste0("lfsr_", trait_names[r])]] <- cs_lfsr[, r]
      }
    } else {
      for (trait in trait_names) {
        cs_variant_summary[[paste0("lfsr_", trait)]] <-
          rep(NA_real_, nrow(cs_variant_summary))
      }
    }

    config_summary <- data.frame(
      outcome = trait_names,
      lbf_outcome = as.numeric(logBF_mat[l, ]),
      sentinel_lfsr = sentinel_lfsr,
      lfsr_pass = sentinel_lfsr_pass,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    out[[i]] <- list(cs_summary = cs_summary,
                     config_summary = config_summary,
                     cs_variant_summary = cs_variant_summary)
  }

  out <- out[!vapply(out, is.null, logical(1))]
  names(out) <- vapply(out, function(x) x$cs_summary$cs[[1]], character(1))
  out
}

# -----------------------------------------------------------------------------
# Coloc pairwise ABF (verbatim port of coloc:::combine.abf).
# -----------------------------------------------------------------------------

# Numerically stable log(sum(exp(x))).
.logsum <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Numerically stable log(exp(a) - exp(b)) for a > b.
.logdiff <- function(a, b) {
  m <- max(a, b)
  m + log(exp(a - m) - exp(b - m))
}

# Compute (PP.H0, PP.H1, PP.H2, PP.H3, PP.H4) from per-SNP log-BF vectors,
# matching coloc:::combine.abf line-for-line.
combine_abf_pair <- function(l1, l2, p1, p2, p12) {
  stopifnot(length(l1) == length(l2))
  lsum    <- l1 + l2
  lH0     <- 0
  lH1     <- log(p1)  + .logsum(l1)
  lH2     <- log(p2)  + .logsum(l2)
  lH3     <- log(p1)  + log(p2) +
             .logdiff(.logsum(l1) + .logsum(l2), .logsum(lsum))
  lH4     <- log(p12) + .logsum(lsum)
  all_lH  <- c(lH0, lH1, lH2, lH3, lH4)
  pp      <- exp(all_lH - .logsum(all_lH))
  names(pp) <- paste0("PP.H", 0:4)
  pp
}

coloc_pairwise_abf <- function(views, p1, p2, p12) {
  N <- length(views)
  empty <- data.frame(trait1 = character(0), trait2 = character(0),
                      l1 = integer(0), l2 = integer(0),
                      hit1 = character(0), hit2 = character(0),
                      PP.H0 = numeric(0), PP.H1 = numeric(0),
                      PP.H2 = numeric(0), PP.H3 = numeric(0),
                      PP.H4 = numeric(0),
                      stringsAsFactors = FALSE)
  if (N < 2L) return(empty)

  trait_names <- vapply(views, function(v) v$name, character(1))
  variant_keys <- lapply(views, .view_variant_keys)
  variant_named <- vapply(views, function(v) {
    !is.null(colnames(v$alpha)) || !is.null(colnames(v$lbf))
  }, logical(1))
  rows <- list()

  for (a in seq_len(N - 1L)) {
    for (b in (a + 1L):N) {
      L1 <- view_cs_indices(views[[a]], cs_only = TRUE)
      L2 <- view_cs_indices(views[[b]], cs_only = TRUE)
      if (length(L1) == 0L || length(L2) == 0L) next

      keys_a <- variant_keys[[a]]
      keys_b <- variant_keys[[b]]
      if ((!variant_named[a] || !variant_named[b]) &&
          length(keys_a) != length(keys_b)) {
        stop("coloc_pairwise: cannot align unnamed variant sets with ",
             "different lengths for traits '", trait_names[a], "' and '",
             trait_names[b], "'.")
      }
      common <- intersect(keys_a, keys_b)
      if (length(common) == 0L) {
        stop("coloc_pairwise: no overlapping variants between traits '",
             trait_names[a], "' and '", trait_names[b], "'.")
      }
      idx_a <- match(common, keys_a)
      idx_b <- match(common, keys_b)
      hit_names <- if (variant_named[a] || variant_named[b]) {
        common
      } else {
        paste0("snp_", idx_a)
      }

      for (i in L1) {
        if (all(views[[a]]$alpha[i, ] == 0)) next
        l1_row <- views[[a]]$lbf[i, idx_a]
        for (j in L2) {
          if (all(views[[b]]$alpha[j, ] == 0)) next
          l2_row <- views[[b]]$lbf[j, idx_b]

          pp <- combine_abf_pair(l1_row, l2_row, p1 = p1, p2 = p2, p12 = p12)

          hit1 <- hit_names[which.max(l1_row)]
          hit2 <- hit_names[which.max(l2_row)]

          rows[[length(rows) + 1L]] <- data.frame(
            trait1 = trait_names[a], trait2 = trait_names[b],
            l1     = i,              l2     = j,
            hit1   = hit1,           hit2   = hit2,
            PP.H0  = pp["PP.H0"],    PP.H1  = pp["PP.H1"],
            PP.H2  = pp["PP.H2"],    PP.H3  = pp["PP.H3"],
            PP.H4  = pp["PP.H4"],
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        }
      }
    }
  }

  if (length(rows) == 0L) return(empty)
  do.call(rbind, rows)
}
# =============================================================================
# Summary / print methods for `susie_post_outcome_configuration` results.
# =============================================================================
#
# Goals:
#   * Be the one-stop pretty-printer so users almost never have to inspect
#     the raw nested list.
#   * Color-code signal vs. no-signal so the eye reads the table at a glance
#     (BOLD DARK GREEN = active / shared, YELLOW = ambiguous, DIM = below
#     threshold; coloc verdicts H4 = green/bold, H3 = magenta, H1/H2 = blue,
#     H0 = dim).
#   * Filter no-signal rows by default (signal_only = TRUE) and footer the
#     hidden count.
#   * Be robust to malformed / partial input objects: missing fields,
#     missing columns, empty lists, length-mismatched per-trait fields,
#     trait names colliding with reserved column names, etc. None of these
#     should error -- they should degrade gracefully.

# Reserved column names that the SuSiEx tidy table adds. Trait names that
# collide get a "trait_" prefix during materialisation.
.SUSIEX_RESERVED <- c("tuple", "top_pattern", "top_prob")

# Coloc PP columns. We tolerate the data.frame missing some, only enforce
# that PP.H0..PP.H4 are present (the source enforces all five).
.COLOC_PP_COLS  <- c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")
.COLOC_DISPLAY  <- c("trait1", "trait2", "l1", "l2", "hit1", "hit2")
.COLOC_LABELS   <- c("H0 no signal",
                     "H1 trait1-only",
                     "H2 trait2-only",
                     "H3 distinct",
                     "H4 shared")

#' Summarise a susie_post_outcome_configuration result
#'
#' Builds tidy tables from the nested list returned by
#' [susie_post_outcome_configuration()] and prints them with ANSI color
#' highlighting via [print.summary.susie_post_outcome_configuration()].
#' The summary itself is an S3 object: index `$susiex` and
#' `$coloc_pairwise` to grab the data.frames for downstream use.
#'
#' Color encoding (when ANSI colors are available):
#' \itemize{
#'   \item SuSiEx per-trait marginal probability: bold dark green when
#'     `>= prob_thresh` (active), yellow when in
#'     `[ambiguous_lower, prob_thresh)`, dim otherwise. The `active`
#'     logical from the raw result is encoded by color and is not shown
#'     as a separate column.
#'   \item Coloc verdict: bold dark green for H4 (shared causal), magenta
#'     for H3 (distinct causals), blue for H1 or H2 (single-trait causal),
#'     dim for H0 (no signal). The dominant PP per row is bolded.
#' }
#'
#' Robustness: this method is defensive against malformed input. Empty
#' lists, NULL components, missing fields, length-mismatched per-trait
#' vectors, trait names that collide with reserved columns
#' (`tuple`, `top_pattern`, `top_prob`), and coloc data.frames that
#' lack some optional columns (`hit1`, `hit2`) all degrade gracefully
#' rather than erroring.
#'
#' @param object Output of [susie_post_outcome_configuration()].
#' @param prob_thresh Threshold above which `marginal_prob` counts as a
#'   signal (default `0.8`).
#' @param ambiguous_lower Lower edge of the "ambiguous" band for the
#'   SuSiEx color coding: marginals in
#'   `[ambiguous_lower, prob_thresh)` are colored yellow. Default `0.5`.
#'   Set to `prob_thresh` to disable the band.
#' @param signal_only Logical. If `TRUE` (default), drop CS tuples where
#'   no trait is active and drop coloc rows whose dominant hypothesis is
#'   H0. Pass `FALSE` to keep everything.
#' @param color One of `"auto"` (default; honors [crayon::has_color()]),
#'   `TRUE` (force colors on), or `FALSE` (force them off).
#' @param ... Ignored.
#'
#' @return A list of class `"summary.susie_post_outcome_configuration"`
#' with components:
#' \describe{
#'   \item{`$susiex`}{`data.frame` (or `NULL` when no signals): one row per
#'     CS tuple. Columns: `tuple` (e.g. `"(1,1,1)"`), one numeric column
#'     per trait carrying that trait's `marginal_prob`, `top_pattern`
#'     (binary configuration string for the most-probable configuration),
#'     `top_prob` (its probability).}
#'   \item{`$coloc_pairwise`}{`data.frame` (or `NULL`): the original coloc
#'     table extended with `verdict` (named hypothesis label) and `top_pp`
#'     (the dominant PP value).}
#'   \item{`$susiex_n_total`, `$susiex_n_kept`, `$coloc_n_total`,
#'     `$coloc_n_kept`}{row counts before and after `signal_only`
#'     filtering, used by the print method to footer hidden rows.}
#'   \item{`$prob_thresh`, `$ambiguous_lower`, `$signal_only`, `$color`}{
#'     parameters echoed for the print method.}
#' }
#'
#' @seealso [susie_post_outcome_configuration()],
#'   [print.summary.susie_post_outcome_configuration()]
#'
#' @method summary susie_post_outcome_configuration
#' @export summary.susie_post_outcome_configuration
#' @export
summary.susie_post_outcome_configuration <- function(
    object,
    prob_thresh     = 0.8,
    ambiguous_lower = 0.5,
    signal_only     = TRUE,
    color           = "auto",
    ...) {
  if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L ||
      !is.finite(prob_thresh) || prob_thresh < 0 || prob_thresh > 1) {
    stop("`prob_thresh` must be a single numeric in [0, 1].")
  }
  if (!is.numeric(ambiguous_lower) || length(ambiguous_lower) != 1L ||
      !is.finite(ambiguous_lower) ||
      ambiguous_lower < 0 || ambiguous_lower > prob_thresh) {
    stop("`ambiguous_lower` must be a single numeric in [0, prob_thresh].")
  }
  if (!is.logical(signal_only) || length(signal_only) != 1L ||
      is.na(signal_only)) {
    stop("`signal_only` must be a single logical (TRUE or FALSE).")
  }
  if (!(isTRUE(color) || isFALSE(color) || identical(color, "auto"))) {
    stop("`color` must be one of TRUE, FALSE, or \"auto\".")
  }

  ses <- .summarise_susiex(object$susiex, signal_only, prob_thresh)
  cls <- .summarise_coloc(object$coloc_pairwise, signal_only)
  out <- list(
    susiex          = ses$df,
    susiex_n_total  = ses$n_total,
    susiex_n_kept   = ses$n_kept,
    coloc_pairwise  = cls$df,
    coloc_n_total   = cls$n_total,
    coloc_n_kept    = cls$n_kept,
    prob_thresh     = prob_thresh,
    ambiguous_lower = ambiguous_lower,
    signal_only     = signal_only,
    color           = color
  )
  class(out) <- c("summary.susie_post_outcome_configuration", "list")
  out
}

# Tidy `configs$susiex` (list of CS-tuple result lists) into a data.frame
# wrapped in a small list with kept/total counts so the print method can
# tell users what was hidden. Returns NULL `df` when input is empty or
# fully filtered.
#
# Defensive against per-tuple field omissions: a tuple missing
# `marginal_prob` or `config_prob` is silently skipped. Trait names that
# collide with reserved columns are prefixed with "trait_". Trait sets
# that vary across tuples are unioned.
.summarise_susiex <- function(susiex, signal_only, prob_thresh) {
  if (is.null(susiex) || !is.list(susiex) || length(susiex) == 0L) {
    return(list(df = NULL, n_total = 0L, n_kept = 0L))
  }
  n_total <- length(susiex)

  # Pull the union of trait names across all tuples (some tuples might be
  # malformed and missing fields; we just skip those).
  trait_names_all <- unique(unlist(lapply(susiex, function(tup) {
    raw <- .susiex_raw(tup)
    if (is.list(raw) && !is.null(raw$marginal_prob)) {
      names(raw$marginal_prob)
    } else character(0)
  })))
  if (length(trait_names_all) == 0L) {
    return(list(df = NULL, n_total = n_total, n_kept = 0L))
  }
  # Avoid collisions with reserved column names by prefixing.
  trait_cols <- ifelse(trait_names_all %in% .SUSIEX_RESERVED,
                       paste0("trait_", trait_names_all),
                       trait_names_all)
  trait_cols <- make.unique(trait_cols)
  names(trait_cols) <- trait_names_all   # raw -> column-name mapping

  rows <- lapply(susiex, function(tup) {
    raw <- .susiex_raw(tup)
    if (!is.list(raw) || is.null(raw$marginal_prob) ||
        is.null(raw$config_prob) || is.null(raw$configs)) {
      return(NULL)
    }
    mp <- raw$marginal_prob
    if (signal_only) {
      # Re-derive active using current prob_thresh (don't trust the stored
      # active flag, which was computed against the call-time threshold).
      if (!any(is.finite(mp) & mp >= prob_thresh)) return(NULL)
    }
    cp <- raw$config_prob
    if (length(cp) == 0L || !all(is.finite(cp))) return(NULL)
    top_idx <- which.max(cp)
    cfg     <- raw$configs
    top_pat <- if (is.matrix(cfg) && nrow(cfg) >= top_idx) {
      paste(cfg[top_idx, ], collapse = "")
    } else NA_character_
    cs_idx_str <- if (!is.null(raw$cs_indices)) {
      paste0("(", paste(raw$cs_indices, collapse = ","), ")")
    } else NA_character_

    out <- data.frame(tuple = cs_idx_str, stringsAsFactors = FALSE)
    for (raw in trait_names_all) {
      out[[trait_cols[[raw]]]] <- if (raw %in% names(mp)) {
        as.numeric(mp[[raw]])
      } else NA_real_
    }
    out$top_pattern <- top_pat
    out$top_prob    <- as.numeric(cp[top_idx])
    out
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(list(df = NULL, n_total = n_total, n_kept = 0L))
  }
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  list(df = df, n_total = n_total, n_kept = nrow(df))
}

.susiex_raw <- function(tup) {
  raw <- attr(tup, "raw", exact = TRUE)
  if (is.list(raw)) raw else tup
}

# Annotate the coloc data.frame with verdict + dominant PP, and optionally
# drop rows whose dominant hypothesis is H0. Returns the df and kept/total
# counts so the print method can footer the hidden count. Tolerates the
# input data.frame already carrying a `verdict` or `top_pp` column (we
# overwrite). Errors if any of PP.H0..PP.H4 is missing -- those columns
# define the algorithm.
.summarise_coloc <- function(df, signal_only) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return(list(df = NULL, n_total = 0L, n_kept = 0L))
  }
  missing_pp <- setdiff(.COLOC_PP_COLS, colnames(df))
  if (length(missing_pp) > 0L) {
    warning("coloc_pairwise table missing required columns: ",
            paste(missing_pp, collapse = ", "),
            "; skipping coloc summary.", call. = FALSE)
    return(list(df = NULL, n_total = nrow(df), n_kept = 0L))
  }

  pp_mat <- as.matrix(df[, .COLOC_PP_COLS, drop = FALSE])
  storage.mode(pp_mat) <- "double"
  # Rows where every PP is NA / non-finite are unscoreable; treat as H0.
  bad_row <- !apply(pp_mat, 1L, function(r) any(is.finite(r)))
  pp_mat[bad_row, ] <- 0
  pp_mat[bad_row, 1L] <- 1

  top    <- max.col(pp_mat, ties.method = "first")
  df$verdict <- .COLOC_LABELS[top]
  df$top_pp  <- pp_mat[cbind(seq_len(nrow(df)), top)]
  n_total    <- nrow(df)

  if (signal_only) {
    df <- df[top != 1L, , drop = FALSE]
    rownames(df) <- NULL
  }
  if (nrow(df) == 0L) {
    return(list(df = NULL, n_total = n_total, n_kept = 0L))
  }
  list(df = df, n_total = n_total, n_kept = nrow(df))
}

#' Print a summary.susie_post_outcome_configuration object
#'
#' Pretty-prints the tidy tables built by
#' [summary.susie_post_outcome_configuration()] with optional ANSI color
#' highlighting. See that page for the color encoding.
#'
#' @param x Output of [summary.susie_post_outcome_configuration()].
#' @param ... Ignored.
#' @return The input invisibly.
#'
#' @seealso [summary.susie_post_outcome_configuration()]
#'
#' @method print summary.susie_post_outcome_configuration
#' @export print.summary.susie_post_outcome_configuration
#' @export
#' @importFrom crayon has_color bold silver green yellow blue magenta cyan
print.summary.susie_post_outcome_configuration <- function(x, ...) {
  use_color <- isTRUE(x$color) ||
    (identical(x$color, "auto") && has_color())

  # Force-enable crayon when the caller asked for colors explicitly. Crayon
  # otherwise respects its own global `crayon.enabled` option and may strip
  # ANSI in non-tty contexts (R CMD CHECK, capture.output, knitr) even when
  # the user passed `color = TRUE`.
  if (isTRUE(x$color)) {
    old_opt <- options(crayon.enabled = TRUE)
    on.exit(options(old_opt), add = TRUE)
  }

  paint <- if (use_color) {
    function(s, style) style(s)
  } else {
    function(s, style) s
  }

  if (is.null(x$susiex) && is.null(x$coloc_pairwise)) {
    cat("susie_post_outcome_configuration: no signals to report",
        if (isTRUE(x$signal_only)) " (signal_only = TRUE)" else "",
        "\n", sep = "")
    return(invisible(x))
  }

  if (!is.null(x$susiex) && nrow(x$susiex) > 0L) {
    cat("\n",
        paint("SuSiEx: per-trait marginal P(active) per CS tuple", bold),
        "\n", sep = "")
    cat(paint(sprintf(
      "  prob_thresh = %.2f, ambiguous = [%.2f, %.2f)",
      x$prob_thresh, x$ambiguous_lower, x$prob_thresh),
      silver), "\n\n", sep = "")
    .print_susiex_table(x$susiex, x$prob_thresh, x$ambiguous_lower, use_color)
    if (isTRUE(x$signal_only) && x$susiex_n_total > x$susiex_n_kept) {
      cat(paint(sprintf(
        "  (%d/%d CS tuples hidden -- no trait above prob_thresh; pass signal_only = FALSE to show)",
        x$susiex_n_total - x$susiex_n_kept, x$susiex_n_total),
        silver), "\n", sep = "")
    }
  }

  if (!is.null(x$coloc_pairwise) && nrow(x$coloc_pairwise) > 0L) {
    cat("\n",
        paint("Coloc pairwise: dominant hypothesis per (trait, trait', l1, l2)",
              bold),
        "\n", sep = "")
    cat(paint(
      "  H0 no signal | H1 trait1-only | H2 trait2-only | H3 distinct | H4 shared",
      silver), "\n\n", sep = "")
    .print_coloc_table(x$coloc_pairwise, use_color)
    if (isTRUE(x$signal_only) && x$coloc_n_total > x$coloc_n_kept) {
      cat(paint(sprintf(
        "  (%d/%d pairs hidden -- H0 dominant; pass signal_only = FALSE to show)",
        x$coloc_n_total - x$coloc_n_kept, x$coloc_n_total),
        silver), "\n", sep = "")
    }
  }

  invisible(x)
}

# ---- table renderers -------------------------------------------------------

.print_susiex_table <- function(df, prob_thresh, ambiguous_lower, use_color) {
  trait_cols <- setdiff(colnames(df), .SUSIEX_RESERVED)

  fmt_prob <- function(p) {
    s <- if (is.na(p)) "  NA" else sprintf("%.3f", p)
    if (!use_color) return(s)
    if (is.na(p))                    silver(s)
    else if (p >= prob_thresh)       bold(green(s))
    else if (p >= ambiguous_lower)   yellow(s)
    else                              silver(s)
  }
  fmt_pat <- function(pat) {
    if (is.na(pat)) return("NA")
    if (!use_color) return(pat)
    cyan(pat)
  }

  hdr  <- c("CS tuple", trait_cols, "top pattern", "top P")
  rows <- lapply(seq_len(nrow(df)), function(i) {
    c(as.character(df$tuple[i]),
      vapply(trait_cols, function(t) fmt_prob(df[[t]][i]), character(1)),
      fmt_pat(df$top_pattern[i]),
      sprintf("%.3f", df$top_prob[i]))
  })
  .print_aligned(hdr, rows)
}

.print_coloc_table <- function(df, use_color) {
  display_present <- intersect(.COLOC_DISPLAY, colnames(df))
  pp_present      <- intersect(.COLOC_PP_COLS, colnames(df))

  pp_mat  <- as.matrix(df[, pp_present, drop = FALSE])
  storage.mode(pp_mat) <- "double"
  top_idx <- max.col(pp_mat, ties.method = "first")

  fmt_pp <- function(p, is_top) {
    s <- if (is.na(p)) "  NA" else sprintf("%.3f", p)
    if (!use_color) return(s)
    if (is.na(p))   silver(s)
    else if (is_top) bold(s)
    else            s
  }
  fmt_verdict <- function(v) {
    if (is.na(v) || !nzchar(v)) return(if (is.na(v)) "NA" else v)
    if (!use_color) return(v)
    style <- switch(
      substr(v, 1L, 2L),
      "H0" = silver,
      "H1" = blue,
      "H2" = blue,
      "H3" = magenta,
      "H4" = function(s) bold(green(s)),
      identity)
    style(v)
  }

  hdr  <- c(display_present, pp_present, "verdict")
  rows <- lapply(seq_len(nrow(df)), function(i) {
    pp_strs <- vapply(seq_along(pp_present), function(k) {
      fmt_pp(pp_mat[i, k], k == top_idx[i])
    }, character(1))
    c(vapply(display_present, function(col) {
        as.character(df[[col]][i])
      }, character(1)),
      pp_strs,
      fmt_verdict(df$verdict[i]))
  })
  .print_aligned(hdr, rows)
}

# Width-aware aligned printing. `vwidth` strips ANSI escape sequences so
# colored cells line up correctly; `pad_to` right-pads to a target width.
.print_aligned <- function(hdr, rows) {
  vwidth <- function(s) nchar(gsub("\033\\[[0-9;]*m", "", s))
  pad_to <- function(s, w) {
    pad <- max(0L, w - vwidth(s))
    paste0(s, strrep(" ", pad))
  }

  ncols <- length(hdr)
  if (length(rows) == 0L) {
    cat("  ", paste(hdr, collapse = "  "), "\n", sep = "")
    return(invisible())
  }
  widths <- vapply(seq_len(ncols), function(k) {
    body_w <- max(vapply(rows, function(r) vwidth(r[[k]]), integer(1)))
    max(vwidth(hdr[k]), body_w)
  }, integer(1))

  cat("  ",
      paste(vapply(seq_len(ncols), function(k) pad_to(hdr[k], widths[k]),
                   character(1)),
            collapse = "  "),
      "\n", sep = "")
  cat("  ",
      paste(strrep("-", widths), collapse = "  "),
      "\n", sep = "")
  for (r in rows) {
    cat("  ",
        paste(vapply(seq_len(ncols), function(k) pad_to(r[[k]], widths[k]),
                     character(1)),
              collapse = "  "),
        "\n", sep = "")
  }
  invisible()
}
