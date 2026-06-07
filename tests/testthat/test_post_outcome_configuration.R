context("susie_post_outcome_configuration engine")

# Engine-level tests for R/susie_post_outcome_configuration.R. The companion
# file test_post_outcome_configuration_summary.R covers the summary()/print()
# methods. Here we exercise the public entry point (validation, dispatch,
# return shape), input normalisation, CS-tuple enumeration, the SuSiEx 2^N
# kernel, the coloc pairwise ABF kernel, and the renderer's color / NA
# branches. devtools::load_all() is active so every internal helper is
# callable by name. Shared fit/tuple/data.frame builders live in
# helper_coloc.R (coloc_fit, coloc_mv_fit, coloc_real_fit, coloc_tuple,
# coloc_post_obj, coloc_row, coloc_df).

# ---- Public entry point: susie_post_outcome_configuration() ----------------

test_that("entry point rejects malformed prob_thresh", {
  fit <- coloc_fit(matrix(0.5, 1, 2), matrix(1, 1, 2))
  expect_error(susie_post_outcome_configuration(fit, prob_thresh = 1.1),
               "prob_thresh")
  expect_error(susie_post_outcome_configuration(fit, prob_thresh = -0.1),
               "prob_thresh")
  expect_error(susie_post_outcome_configuration(fit, prob_thresh = c(0.5, 0.6)),
               "prob_thresh")
  expect_error(susie_post_outcome_configuration(fit, prob_thresh = NA_real_),
               "prob_thresh")
  expect_error(susie_post_outcome_configuration(fit, prob_thresh = "x"),
               "prob_thresh")
})

test_that("entry point rejects malformed cs_only", {
  fit <- coloc_fit(matrix(0.5, 1, 2), matrix(1, 1, 2))
  expect_error(susie_post_outcome_configuration(fit, cs_only = NA),
               "cs_only")
  expect_error(susie_post_outcome_configuration(fit, cs_only = "yes"),
               "cs_only")
  expect_error(susie_post_outcome_configuration(fit, cs_only = c(TRUE, FALSE)),
               "cs_only")
})

test_that("entry point uses provided outcome_names for trait labels", {
  fits <- list(coloc_fit(matrix(c(0, 1), 1, 2), matrix(c(0, 4), 1, 2)),
               coloc_fit(matrix(c(1, 0), 1, 2), matrix(c(3, 0), 1, 2)))
  res_sx <- susie_post_outcome_configuration(fits, method = "susiex",
                                             outcome_names = c("MDD", "BP"))
  expect_named(res_sx$susiex[[1]]$activation_summary, c("MDD", "BP"))

  res_coloc <- susie_post_outcome_configuration(fits,
                                                method = "coloc_pairwise",
                                                outcome_names = c("MDD", "BP"))
  expect_equal(res_coloc$coloc_pairwise$trait1, "MDD")
  expect_equal(res_coloc$coloc_pairwise$trait2, "BP")
})

test_that("entry point rejects malformed outcome_names", {
  fits <- list(coloc_fit(matrix(0.5, 1, 2), matrix(1, 1, 2)),
               coloc_fit(matrix(0.5, 1, 2), matrix(1, 1, 2)))
  expect_error(susie_post_outcome_configuration(fits,
                                                outcome_names = "one"),
               "outcome_names")
  expect_error(susie_post_outcome_configuration(fits,
                                                outcome_names = c("one", NA)),
               "outcome_names")
  expect_error(susie_post_outcome_configuration(fits,
                                                outcome_names = c("one", "")),
               "outcome_names")
})

test_that("entry point skips traits with no credible sets", {
  no_cs <- structure(list(alpha = matrix(c(0, 1), 1, 2),
                          lbf_variable = matrix(c(0, 4), 1, 2)),
                     class = "susie")
  has_a <- coloc_fit(matrix(c(1, 0), 1, 2), matrix(c(3, 0), 1, 2))
  has_b <- coloc_fit(matrix(c(0, 1), 1, 2), matrix(c(0, 4), 1, 2))
  expect_warning(
    res <- susie_post_outcome_configuration(
      list(no_cs, has_a, has_b),
      method = "susiex",
      outcome_names = c("drop_me", "A", "B")),
    "drop_me"
  )
  expect_s3_class(res, "susie_post_outcome_configuration")
  expect_named(res$susiex[[1]]$activation_summary, c("A", "B"))
})

test_that("entry point returns NULL when all traits have no credible sets", {
  no_cs <- structure(list(alpha = matrix(c(0, 1), 1, 2),
                          lbf_variable = matrix(c(0, 4), 1, 2)),
                     class = "susie")
  res <- "not null"
  expect_message(
    expect_warning(
      res <- susie_post_outcome_configuration(
        list(no_cs, no_cs),
        method = "susiex",
        outcome_names = c("a", "b")),
      "a, b"
    ),
    "returning NULL"
  )
  expect_null(res)
})

test_that("entry point skips no-CS traits for coloc_pairwise", {
  no_cs <- structure(list(alpha = matrix(c(0, 1), 1, 2),
                          lbf_variable = matrix(c(0, 4), 1, 2)),
                     class = "susie")
  has_a <- coloc_fit(matrix(c(1, 0), 1, 2), matrix(c(3, 0), 1, 2))
  has_b <- coloc_fit(matrix(c(0, 1), 1, 2), matrix(c(0, 4), 1, 2))
  expect_warning(
    res <- susie_post_outcome_configuration(
      list(no_cs, has_a, has_b),
      method = "coloc_pairwise",
      outcome_names = c("drop_me", "A", "B")),
    "drop_me"
  )
  expect_s3_class(res, "susie_post_outcome_configuration")
  expect_equal(unique(res$coloc_pairwise$trait1), "A")
  expect_equal(unique(res$coloc_pairwise$trait2), "B")
})

test_that("entry point returns NULL when coloc_pairwise has fewer than two CS traits", {
  no_cs <- structure(list(alpha = matrix(c(0, 1), 1, 2),
                          lbf_variable = matrix(c(0, 4), 1, 2)),
                     class = "susie")
  has_cs <- coloc_fit(matrix(c(1, 0), 1, 2), matrix(c(3, 0), 1, 2))
  res <- "not null"
  expect_message(
    expect_warning(
      res <- susie_post_outcome_configuration(
        list(no_cs, has_cs),
        method = "coloc_pairwise",
        outcome_names = c("drop_me", "keep_me")),
      "drop_me"
    ),
    "Fewer than two trait views"
  )
  expect_null(res)
})

test_that("entry point rejects malformed p1/p2/p12 priors", {
  fit <- coloc_fit(matrix(0.5, 1, 2), matrix(1, 1, 2))
  expect_error(
    susie_post_outcome_configuration(fit, method = "coloc_pairwise", p1 = 0),
    "`p1` must be a single numeric")
  expect_error(
    susie_post_outcome_configuration(fit, method = "coloc_pairwise", p2 = 1),
    "`p2` must be a single numeric")
  expect_error(
    susie_post_outcome_configuration(fit, method = "coloc_pairwise", p12 = -1),
    "`p12` must be a single numeric")
  expect_error(
    susie_post_outcome_configuration(fit, method = "coloc_pairwise",
                                     p1 = c(1e-4, 1e-4)),
    "`p1` must be a single numeric")
  expect_error(
    susie_post_outcome_configuration(fit, method = "coloc_pairwise",
                                     p12 = NA_real_),
    "`p12` must be a single numeric")
  expect_error(
    susie_post_outcome_configuration(fit, method = "susiex",
                                     single_effect_lfsr_cutoff = NA_real_),
    "`single_effect_lfsr_cutoff` must be a single numeric")
})

test_that("method = 'susiex' returns tagged object with $susiex component", {
  fits <- list(coloc_real_fit(1), coloc_real_fit(2), coloc_real_fit(3))
  res  <- susie_post_outcome_configuration(fits, method = "susiex", by = "fit",
                                           prob_thresh = 0.8)
  expect_s3_class(res, "susie_post_outcome_configuration")
  expect_true("susiex" %in% names(res))
  expect_false("coloc_pairwise" %in% names(res))
  expect_identical(attr(res, "method"), "susiex")
  expect_equal(attr(res, "prob_thresh"), 0.8)
  expect_true(length(res$susiex) >= 1L)
  # Each tuple carries the documented fields.
  expect_named(res$susiex[[1]], c("config_probability", "activation_summary"))
  expect_named(attr(res$susiex[[1]], "raw"),
               c("cs_indices", "cs_labels", "logBF_trait", "configs",
                 "config_prob", "marginal_prob", "active"))
})

test_that("method = 'coloc_pairwise' returns tagged object with coloc df", {
  fits <- list(coloc_real_fit(1), coloc_real_fit(2))
  res  <- susie_post_outcome_configuration(fits, method = "coloc_pairwise",
                                           by = "fit")
  expect_s3_class(res, "susie_post_outcome_configuration")
  expect_true("coloc_pairwise" %in% names(res))
  expect_false("susiex" %in% names(res))
  expect_identical(attr(res, "method"), "coloc_pairwise")
  expect_s3_class(res$coloc_pairwise, "data.frame")
  expect_named(res$coloc_pairwise,
               c("trait1", "trait2", "l1", "l2", "hit1", "hit2",
                 "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"))
  pp <- as.matrix(res$coloc_pairwise[, paste0("PP.H", 0:4)])
  expect_equal(unname(rowSums(pp)), rep(1, nrow(pp)), tolerance = 1e-8)
})

test_that("method = 'coloc_pairwise' errors when named variants do not overlap", {
  f1 <- coloc_fit(matrix(1, 1, 1, dimnames = list(NULL, "a")),
                  matrix(2, 1, 1, dimnames = list(NULL, "a")))
  f2 <- coloc_fit(matrix(1, 1, 1, dimnames = list(NULL, "b")),
                  matrix(2, 1, 1, dimnames = list(NULL, "b")))
  expect_error(
    susie_post_outcome_configuration(list(f1, f2), method = "coloc_pairwise"),
    "no overlapping variants")
})

test_that("method = 'susiex' returns mvSuSiE CS summaries with SuSiEx activation", {
  A <- matrix(0, 3, 4)
  A[1, 1] <- 1; A[2, 2] <- 1; A[3, 3] <- 1
  LV <- matrix(0, 3, 4)
  lbf_out <- matrix(c(4, -1,
                      1,  2,
                      3,  0), nrow = 3, byrow = TRUE)
  colnames(lbf_out) <- c("old", "new")
  lfsr <- matrix(c(0.01, 1,
                   0.2,  0.03,
                   0.04, 0.5), nrow = 3, byrow = TRUE,
                 dimnames = list(paste0("L", 1:3), c("old", "new")))
  fit <- coloc_mv_fit(A, LV, array(0, dim = c(3, 4, 2)),
                      cs = list(L3 = c(3L, 4L), L1 = 1L),
                      lbf_outcome = lbf_out,
                      single_effect_lfsr = lfsr)
  fit$pip <- c(0.8, 0.1, 0.2, 0.9)
  fit$lfsr <- matrix(c(0.001, 0.2,
                       0.2,   0.2,
                       0.2,   0.2,
                       0.2,   0.003), nrow = 4, byrow = TRUE)
  fit$lbf <- c(10, 2, 5)
  fit$variable_names <- paste0("s", 1:4)
  fit$sets$purity <- data.frame(
    min.abs.corr = c(0.77, 0.91),
    mean.abs.corr = c(0.84, 0.91),
    median.abs.corr = c(0.84, 0.91),
    row.names = c("L3", "L1")
  )

  res <- susie_post_outcome_configuration(fit, method = "susiex")

  expect_s3_class(res, "susie_post_outcome_configuration")
  expect_true("susiex" %in% names(res))
  expect_false("mvsusie" %in% names(res))
  expect_identical(attr(res, "method"), "susiex")
  expect_type(res$susiex, "list")
  expect_named(res$susiex, c("L3", "L1"))
  expect_named(res$susiex$L3,
               c("mvsusie_cs_summary", "mvsusie_config_summary", "mvsusie_cs_variant_summary",
                 "susiex_config_probability", "susiex_activation_summary"))

  expect_s3_class(res$susiex$L3$mvsusie_cs_summary, "data.frame")
  expect_equal(res$susiex$L3$mvsusie_cs_summary$cs, "L3")
  expect_equal(res$susiex$L3$mvsusie_cs_summary$n_variant, 2L)
  expect_equal(res$susiex$L3$mvsusie_cs_summary$purity, 0.77)
  expect_equal(res$susiex$L3$mvsusie_cs_summary$hit, "s4")
  expect_equal(res$susiex$L3$mvsusie_cs_summary$maxPIP, 0.9)
  expect_equal(res$susiex$L3$mvsusie_cs_summary$lbf, 5)
  expect_equal(res$susiex$L3$mvsusie_cs_summary$n_lfsr_outcome, 1L)

  expect_s3_class(res$susiex$L3$mvsusie_config_summary, "data.frame")
  expect_named(res$susiex$L3$mvsusie_config_summary,
               c("outcome", "lbf_outcome", "sentinel_variant",
                 "sentinel_lfsr", "lfsr_pass", "lfsr_cutoff"))
  expect_equal(res$susiex$L3$mvsusie_config_summary$outcome, c("old", "new"))
  expect_equal(res$susiex$L3$mvsusie_config_summary$lbf_outcome, c(3, 0))
  expect_equal(res$susiex$L3$mvsusie_config_summary$sentinel_variant, c("s4", "s4"))
  expect_equal(res$susiex$L3$mvsusie_config_summary$sentinel_lfsr, c(0.2, 0.003))
  expect_equal(res$susiex$L3$mvsusie_config_summary$lfsr_pass, c(FALSE, TRUE))
  expect_equal(res$susiex$L3$mvsusie_config_summary$lfsr_cutoff, c(0.05, 0.05))
  expect_identical(attr(res, "single_effect_lfsr_cutoff"), 0.05)
  expect_false("activation_summary" %in% names(res$susiex$L3))
  expect_true("susiex_activation_summary" %in% names(res$susiex$L3))
  expect_false("posthoc_prob" %in% names(res$susiex$L3))
  expect_false("single_effect_lfsr" %in%
                 colnames(res$susiex$L3$mvsusie_config_summary))

  expect_s3_class(res$susiex$L3$mvsusie_cs_variant_summary, "data.frame")
  expect_named(res$susiex$L3$mvsusie_cs_variant_summary,
               c("variant", "pip", "lfsr_old", "lfsr_new"))
  expect_equal(res$susiex$L3$mvsusie_cs_variant_summary$variant, c("s3", "s4"))
  expect_equal(res$susiex$L3$mvsusie_cs_variant_summary$pip, c(0.2, 0.9))
  expect_equal(res$susiex$L3$mvsusie_cs_variant_summary$lfsr_old, c(0.2, 0.2))
  expect_equal(res$susiex$L3$mvsusie_cs_variant_summary$lfsr_new, c(0.2, 0.003))
  expect_s3_class(res$susiex$L3$susiex_config_probability, "data.frame")
  expect_named(res$susiex$L3$susiex_config_probability,
               c("old", "new", "config_prob"))
  expect_equal(nrow(res$susiex$L3$susiex_config_probability), 4L)
  expect_s3_class(res$susiex$L3$susiex_activation_summary, "data.frame")
  expect_named(res$susiex$L3$susiex_activation_summary, c("old", "new"))
  expect_equal(
    rownames(res$susiex$L3$susiex_activation_summary),
    c("cs_indices", "logBF_trait", "posthoc_prob", "active"))
  expect_equal(
    as.character(res$susiex$L3$susiex_activation_summary["cs_indices", ]),
    c("L3", "L3"))
  expect_named(attr(res$susiex$L3, "raw"),
               c("cs_indices", "cs_labels", "logBF_trait", "configs",
                 "config_prob", "marginal_prob", "active"))

  strict <- susie_post_outcome_configuration(
    fit, method = "susiex", single_effect_lfsr_cutoff = 0.01)
  expect_equal(strict$susiex$L3$mvsusie_config_summary$lfsr_pass, c(FALSE, TRUE))
})

test_that("method = 'susiex' can compute lbf_outcome from alpha-weighted outcome LBFs", {
  A <- matrix(c(0.25, 0.75,
                1.00, 0.00), nrow = 2, byrow = TRUE)
  LV <- matrix(0, 2, 2)
  arr <- array(0, dim = c(2, 2, 2),
               dimnames = list(NULL, NULL, c("t1", "t2")))
  arr[1, , 1] <- c(2, 6)
  arr[1, , 2] <- c(4, 8)
  arr[2, , 1] <- c(3, 9)
  arr[2, , 2] <- c(5, 7)
  fit <- coloc_mv_fit(A, LV, arr, cs = list(L1 = 1L, L2 = 2L))

  res <- susie_post_outcome_configuration(fit, method = "susiex")

  expect_equal(res$susiex$L1$mvsusie_config_summary$lbf_outcome, c(5, 7))
  expect_equal(res$susiex$L2$mvsusie_config_summary$lbf_outcome, c(3, 5))
  expect_true(all(c("cs", "n_variant", "purity", "hit", "maxPIP") %in%
                    colnames(res$susiex$L1$mvsusie_cs_summary)))
})

test_that("method = 'susiex' handles mfsusie with the multi-output contract", {
  fit <- coloc_mv_fit(
    alpha = matrix(c(1, 0), 1, 2),
    lbf = matrix(c(2, 0), 1, 2),
    lbf_variable_outcome = array(c(4, 0, 5, 0), dim = c(1, 2, 2),
                                 dimnames = list(NULL, c("s1", "s2"),
                                                 c("old", "new"))),
    cs = list(L1 = 1L)
  )
  class(fit) <- "mfsusie"
  fit$pip <- c(0.8, 0.1)
  fit$lfsr <- matrix(c(0.01, 0.02,
                       0.2,  0.2), nrow = 2, byrow = TRUE)
  fit$variable_names <- c("s1", "s2")

  res <- susie_post_outcome_configuration(fit, method = "susiex")

  expect_true("susiex" %in% names(res))
  expect_named(res$susiex, "L1")
  expect_named(res$susiex$L1$susiex_activation_summary, c("old", "new"))
  expect_equal(nrow(res$susiex$L1$susiex_config_probability), 4L)
})

test_that("method = 'susiex' covers input and missing-data edge cases", {
  A <- matrix(c(1, 0), 1, 2)
  LV <- matrix(0, 1, 2)
  fit <- coloc_mv_fit(A, LV, array(0, dim = c(1, 2, 2)),
                      cs = list(L1 = 1L),
                      lbf_outcome = as.data.frame(matrix(c(4, 2), 1, 2)))
  fit$pip <- c(0.7, 0.2)
  fit$lfsr <- data.frame(old = c(0.01, 0.2), new = c(0.2, 0.03))
  fit$lbf <- 4
  fit$variable_names <- c("s1", "s2")

  res <- susie_post_outcome_configuration(
    list(fit), method = "susiex", outcome_names = c("old", "new"))
  expect_equal(res$susiex$L1$mvsusie_config_summary$outcome, c("old", "new"))
  expect_equal(res$susiex$L1$mvsusie_config_summary$sentinel_lfsr, c(0.01, 0.2))
  expect_named(res$susiex$L1$susiex_activation_summary, c("old", "new"))

  expect_error(
    susie_post_outcome_configuration(fit, method = "mvsusie"),
    "'arg' should be one of")
  expect_error(
    susie_post_outcome_configuration(fit, method = "susiex",
                                     outcome_names = c("old", "")),
    "outcome_names")

  one_outcome <- fit
  one_outcome$lbf_outcome <- matrix(1, 1, 1)
  one_outcome$lfsr <- matrix(c(0.1, 0.2), 2, 1)
  expect_message(
    expect_null(susie_post_outcome_configuration(one_outcome,
                                                 method = "susiex")),
    "Fewer than two outcomes")

  no_cs <- fit
  no_cs$sets$cs <- list()
  expect_message(
    expect_null(susie_post_outcome_configuration(no_cs, method = "susiex")),
    "No credible sets")

  bad_lbf <- fit
  bad_lbf$lbf_outcome <- matrix(1, 2, 2)
  expect_error(
    susie_post_outcome_configuration(bad_lbf, method = "susiex"),
    "one row per")

  missing_lbf <- fit
  missing_lbf$lbf_outcome <- NULL
  missing_lbf$lbf_variable_outcome <- NULL
  expect_error(
    susie_post_outcome_configuration(missing_lbf, method = "susiex"),
    "requires.*lbf_outcome")

  malformed_lbf <- fit
  malformed_lbf$lbf_outcome <- NULL
  malformed_lbf$lbf_variable_outcome <- array(0, dim = c(2, 2, 2))
  expect_error(
    susie_post_outcome_configuration(malformed_lbf, method = "susiex"),
    "must be L x J")
})

test_that("method = 'susiex' covers CS member edge cases", {
  A <- matrix(c(1, 0), 1, 2)
  LV <- matrix(0, 1, 2)
  fit <- coloc_mv_fit(A, LV, array(0, dim = c(1, 2, 2)),
                      cs = list(L1 = c("s1", "s2")),
                      lbf_outcome = matrix(c(4, 2), 1, 2))
  fit$pip <- c(0.7, 0.2)
  fit$lfsr <- matrix(c(0.01, 0.2,
                       0.2,  0.03), 2, 2, byrow = TRUE)
  fit$lbf <- 4
  fit$variable_names <- c("s1", "s2")

  empty_hit <- mvsusie_cs_hit(
    structure(list(sets = list(cs = NULL)), class = "mvsusie"),
    label = "L1", idx = 1L, variable_names = c("s1", "s2"))
  expect_equal(empty_hit$n_cs, 0L)
  expect_true(is.na(empty_hit$hit))

  missing_hit <- mvsusie_cs_hit(fit, label = "missing", idx = 99L,
                                variable_names = c("s1", "s2"))
  expect_equal(missing_hit$n_cs, 0L)
  expect_true(is.na(missing_hit$hit))

  res <- susie_post_outcome_configuration(fit, method = "susiex")
  expect_equal(res$susiex$L1$mvsusie_cs_variant_summary$variant, c("s1", "s2"))
  expect_equal(res$susiex$L1$mvsusie_cs_variant_summary$pip, c(0.7, 0.2))
  expect_true(is.na(res$susiex$L1$mvsusie_cs_summary$purity))

  fit$sets$purity <- matrix(0.61, nrow = 1)
  expect_equal(mvsusie_cs_purity(fit, label = "L1", idx = 1L), 0.61)
  expect_true(is.na(mvsusie_cs_purity(fit, label = "missing", idx = 99L)))

  unnamed_cs <- fit
  unnamed_cs$sets$cs <- list(2L)
  unnamed_cs$sets$purity <- data.frame(min.abs.corr = 0.42)
  attr(unnamed_cs$sets$cs, "cs_idx") <- 1L
  res_unnamed <- susie_post_outcome_configuration(unnamed_cs,
                                                  method = "susiex")
  expect_equal(res_unnamed$susiex$L1$mvsusie_cs_summary$hit, "s2")
  expect_equal(res_unnamed$susiex$L1$mvsusie_cs_summary$purity, 0.42)

  no_pip <- fit
  no_pip$pip <- NULL
  res_no_pip <- susie_post_outcome_configuration(no_pip, method = "susiex")
  expect_true(all(is.na(res_no_pip$susiex$L1$mvsusie_cs_variant_summary$pip)))
  expect_true(all(is.na(res_no_pip$susiex$L1$mvsusie_cs_summary$hit)))

  bad_members <- fit
  bad_members$sets$cs <- list(L1 = 99L)
  res_bad_members <- susie_post_outcome_configuration(bad_members,
                                                      method = "susiex")
  expect_equal(res_bad_members$susiex$L1$mvsusie_cs_summary$n_variant, 1L)
  expect_equal(nrow(res_bad_members$susiex$L1$mvsusie_cs_variant_summary), 0L)

  zero_alpha <- coloc_mv_fit(matrix(c(0, 0,
                                      0, 1), 2, 2, byrow = TRUE),
                             matrix(0, 2, 2),
                             array(0, dim = c(2, 2, 2)),
                             cs = list(L1 = 1L, L2 = 2L),
                             lbf_outcome = matrix(c(1, 2,
                                                    3, 4), 2, 2,
                                                  byrow = TRUE))
  zero_alpha$pip <- c(0.1, 0.9)
  zero_alpha$lfsr <- fit$lfsr
  zero_alpha$variable_names <- c("s1", "s2")
  res_zero_alpha <- susie_post_outcome_configuration(zero_alpha,
                                                     method = "susiex")
  expect_named(res_zero_alpha$susiex, "L2")
})

test_that("entry point returns NULL for a single fit", {
  fit <- coloc_real_fit(7)
  res <- "not null"
  expect_message(
    res <- susie_post_outcome_configuration(fit, method = "susiex", by = "fit"),
    "Fewer than two trait views"
  )
  expect_null(res)
})

# ---- is_susie_fit() --------------------------------------------------------

test_that("is_susie_fit recognises susie / mvsusie / mfsusie only", {
  expect_true(is_susie_fit(structure(list(), class = "susie")))
  expect_true(is_susie_fit(structure(list(), class = "mvsusie")))
  expect_true(is_susie_fit(structure(list(), class = "mfsusie")))
  expect_false(is_susie_fit(structure(list(), class = "lm")))
  expect_false(is_susie_fit(list()))
  expect_false(is_susie_fit(42))
})

test_that("internal helpers cover defensive display branches", {
  expect_equal(view_cs_label(list(sets_cs = NULL), 2L), "L2")
  expect_equal(view_cs_label(list(sets_cs = list(1L)), 1L), "L1")

  views <- list(
    list(name = "a", lbf = matrix(1, 1, 1)),
    list(name = "b", lbf = matrix(1, 1, 1))
  )
  expect_equal(
    .susiex_config_logbf(c(1L, 1L), c(1L, 1L), views,
                         variant_keys = list("a_only", "b_only"),
                         variant_space_size = 2L),
    -Inf
  )

  expect_identical(.organize_susiex_output(list("plain"))[[1]], "plain")

  tup_fallback <- coloc_tuple(c("a", "b"), cs_indices = c(1, 2),
                              marginal_prob = c(0.9, 0.8))
  names(tup_fallback$cs_indices) <- NULL
  tup_fallback$cs_labels <- setNames(c("L1", "L2"), c("a", "b"))
  org_fallback <- .organize_susiex_output(list(tup_fallback))[[1]]
  expect_named(org_fallback$activation_summary, c("a", "b"))

  tup_no_labels <- coloc_tuple(c("a", "b"), cs_indices = c(1, 2),
                               marginal_prob = c(0.9, 0.8))
  tup_no_labels$cs_labels <- NULL
  org_no_labels <- .organize_susiex_output(list(tup_no_labels))[[1]]
  expect_equal(as.character(org_no_labels$activation_summary["cs_indices", ]),
               c("L1", "L2"))

  wrapped <- list(shown = FALSE)
  attr(wrapped, "raw") <- list(shown = TRUE)
  expect_equal(.susiex_raw(wrapped), list(shown = TRUE))
})

# ---- normalise_to_views() --------------------------------------------------

test_that("normalise_to_views wraps a single fit into a one-element view list", {
  fit   <- coloc_fit(matrix(0.5, 2, 3), matrix(1, 2, 3))
  views <- normalise_to_views(fit, by = "fit")
  expect_length(views, 1L)
  expect_equal(views[[1]]$name, "trait_1")
})

test_that("normalise_to_views errors on an empty list", {
  expect_error(normalise_to_views(list(), by = "fit"),
               "non-empty list")
})

test_that("normalise_to_views errors when an element is not a SuSiE fit", {
  fit <- coloc_fit(matrix(0.5, 1, 2), matrix(1, 1, 2))
  expect_error(
    normalise_to_views(list(fit, 1), by = "fit"),
    "Element 2 of `input` is not a SuSiE-class fit")
})

test_that("normalise_to_views does not require $sets$cs", {
  bare <- structure(list(alpha = matrix(0.5, 1, 2),
                         lbf_variable = matrix(1, 1, 2)),
                    class = "susie")
  views <- normalise_to_views(bare, by = "fit")
  expect_length(views, 1L)
})

test_that("normalise_to_views keeps explicit names and defaults unnamed ones", {
  fitA <- coloc_fit(matrix(0.5, 1, 2), matrix(1, 1, 2))
  fitB <- coloc_fit(matrix(0.5, 1, 2), matrix(1, 1, 2))
  # Named list -> names preserved.
  v_named <- normalise_to_views(list(gene = fitA, qtl = fitB), by = "fit")
  expect_equal(vapply(v_named, function(v) v$name, character(1)),
               c("gene", "qtl"))
  # Unnamed list -> trait_1, trait_2.
  v_unnamed <- normalise_to_views(list(fitA, fitB), by = "fit")
  expect_equal(vapply(v_unnamed, function(v) v$name, character(1)),
               c("trait_1", "trait_2"))
  # Partially named -> blank entries get the default.
  v_partial <- normalise_to_views(list(gene = fitA, fitB), by = "fit")
  expect_equal(vapply(v_partial, function(v) v$name, character(1)),
               c("gene", "trait_2"))
})

# ---- expand_one_fit() ------------------------------------------------------

test_that("expand_one_fit by = 'fit' yields one view from the joint lbf", {
  fit   <- coloc_fit(matrix(0.5, 2, 3), matrix(seq_len(6), 2, 3))
  views <- expand_one_fit(fit, "g", by = "fit")
  expect_length(views, 1L)
  expect_equal(views[[1]]$name, "g")
  expect_equal(views[[1]]$lbf, matrix(seq_len(6), 2, 3))
})

test_that("expand_one_fit by = 'outcome' fans an mvsusie into per-outcome views", {
  set.seed(41)
  L <- 2; J <- 4; R <- 3
  A   <- matrix(0, L, J); A[1, 2] <- 1; A[2, 3] <- 1
  LV  <- matrix(rnorm(L * J), L, J)
  arr <- array(seq_len(L * J * R), dim = c(L, J, R))
  dimnames(arr) <- list(NULL, NULL, c("o1", "o2", "o3"))
  fit <- coloc_mv_fit(A, LV, arr, cs = list(L1 = 1L, L2 = 2L))

  views <- expand_one_fit(fit, "g", by = "outcome")
  expect_length(views, R)
  expect_equal(vapply(views, function(v) v$name, character(1)),
               c("g_o1", "g_o2", "g_o3"))
  # The r-th view reads the r-th slice of lbf_variable_outcome.
  expect_equal(views[[2]]$lbf, arr[, , 2])
  # alpha is shared across outcome views.
  expect_equal(views[[1]]$alpha, A)
})

test_that("expand_one_fit by = 'outcome' names outcomes when array has no dimnames", {
  set.seed(42)
  L <- 2; J <- 3; R <- 2
  A   <- matrix(0, L, J); A[1, 2] <- 1; A[2, 3] <- 1
  LV  <- matrix(rnorm(L * J), L, J)
  arr <- array(rnorm(L * J * R), dim = c(L, J, R))   # no dimnames
  fit <- coloc_mv_fit(A, LV, arr, cs = list(L1 = 1L, L2 = 2L))
  views <- expand_one_fit(fit, "t", by = "outcome")
  expect_equal(vapply(views, function(v) v$name, character(1)),
               c("t_outcome_1", "t_outcome_2"))
})

test_that("expand_one_fit by = 'outcome' errors on mvsusie without lbf_variable_outcome", {
  fit <- structure(list(alpha = matrix(0.5, 1, 2),
                        lbf_variable = matrix(1, 1, 2),
                        sets = list(cs = list(L1 = 1L))),
                   class = "mvsusie")
  expect_error(expand_one_fit(fit, "g", by = "outcome"),
               "requires `\\$lbf_variable_outcome`")
})

test_that("expand_one_fit by = 'outcome' on a plain susie yields one view", {
  fit   <- coloc_fit(matrix(0.5, 2, 3), matrix(1, 2, 3))
  views <- expand_one_fit(fit, "g", by = "outcome")
  expect_length(views, 1L)
  expect_equal(views[[1]]$name, "g")
})

# ---- make_view() -----------------------------------------------------------

test_that("make_view errors when alpha or lbf is NULL", {
  expect_error(make_view("t", NULL, matrix(1, 1, 1), NULL),
               "must be non-null")
  expect_error(make_view("t", matrix(1, 1, 1), NULL, NULL),
               "must be non-null")
})

test_that("make_view coerces non-matrix alpha / lbf to matrices", {
  v <- make_view("t", c(0.4, 0.6), c(1, 2), NULL)
  expect_true(is.matrix(v$alpha))
  expect_true(is.matrix(v$lbf))
  expect_equal(dim(v$alpha), dim(v$lbf))
})

test_that("make_view errors on alpha / lbf shape mismatch", {
  expect_error(
    make_view("t", matrix(1, 2, 2), matrix(1, 2, 3), NULL),
    "identical shape")
})

# ---- view_cs_indices() -----------------------------------------------------

test_that("view_cs_indices returns all rows when cs_only = FALSE", {
  v <- list(alpha = matrix(0, 5, 3), lbf = matrix(0, 5, 3), sets_cs = NULL)
  expect_equal(view_cs_indices(v, cs_only = FALSE), 1:5)
})

test_that("view_cs_indices uses the cs_idx attribute when present", {
  sc <- structure(list(L1 = 1L), cs_idx = c(1L, 3L))
  v  <- list(alpha = matrix(0, 5, 3), lbf = matrix(0, 5, 3), sets_cs = sc)
  expect_equal(view_cs_indices(v, cs_only = TRUE), c(1L, 3L))
})

test_that("view_cs_indices falls back to L-prefixed names", {
  sc <- list(L1 = 1L, L3 = 1L)
  v  <- list(alpha = matrix(0, 5, 3), lbf = matrix(0, 5, 3), sets_cs = sc)
  expect_equal(view_cs_indices(v, cs_only = TRUE), c(1L, 3L))
})

test_that("view_cs_indices treats empty sets_cs as no CS", {
  v0 <- list(alpha = matrix(0, 3, 3), lbf = matrix(0, 3, 3), sets_cs = list())
  expect_equal(view_cs_indices(v0, cs_only = TRUE), integer(0))
})

test_that("view_cs_indices falls back to seq_len when sets_cs is unnamed", {
  vu <- list(alpha = matrix(0, 4, 3), lbf = matrix(0, 4, 3),
             sets_cs = list(1L, 1L))
  expect_equal(view_cs_indices(vu, cs_only = TRUE), 1:4)
})

test_that("view_cs_label uses original CS names", {
  v <- list(alpha = matrix(0, 5, 3), lbf = matrix(0, 5, 3),
            sets_cs = list(L2 = 1L, L4 = 1L))
  expect_equal(view_cs_label(v, 2L), "L2")
  expect_equal(view_cs_label(v, 4L), "L4")

  sc <- structure(list(custom = 1L), cs_idx = 3L)
  v_attr <- list(alpha = matrix(0, 5, 3), lbf = matrix(0, 5, 3),
                 sets_cs = sc)
  expect_equal(view_cs_label(v_attr, 3L), "custom")
})

test_that("view_cs_indices drops out-of-range CS indices", {
  sc <- structure(list(), cs_idx = c(2L, 99L, 7L))   # 99, 7 > nrow = 5
  v  <- list(alpha = matrix(0, 5, 3), lbf = matrix(0, 5, 3), sets_cs = sc)
  expect_equal(view_cs_indices(v, cs_only = TRUE), 2L)
})

# ---- enumerate_cs_tuples() -------------------------------------------------

test_that("enumerate_cs_tuples returns list() when any view has 0 CS indices", {
  # First view has only an out-of-range cs_idx -> 0 indices.
  v <- list(
    list(name = "a", alpha = matrix(0.5, 2, 3), lbf = matrix(1, 2, 3),
         sets_cs = structure(list(), cs_idx = 99L)),
    list(name = "b", alpha = matrix(0.5, 2, 3), lbf = matrix(1, 2, 3),
         sets_cs = list(L1 = 1L)))
  expect_equal(enumerate_cs_tuples(v, by = "fit", cs_only = TRUE), list())
})

test_that("enumerate_cs_tuples by = 'outcome' uses the diagonal of common CSs", {
  v <- list(
    list(name = "a", alpha = matrix(0.5, 3, 3), lbf = matrix(1, 3, 3),
         sets_cs = list(L1 = 1L, L2 = 2L)),
    list(name = "b", alpha = matrix(0.5, 3, 3), lbf = matrix(1, 3, 3),
         sets_cs = list(L1 = 1L, L2 = 2L)))
  tuples <- enumerate_cs_tuples(v, by = "outcome", cs_only = TRUE)
  expect_length(tuples, 2L)                 # common = {1, 2}
  expect_equal(tuples[[1]], c(1L, 1L))
  expect_equal(tuples[[2]], c(2L, 2L))
})

test_that("enumerate_cs_tuples by = 'fit' uses the cross product", {
  v <- list(
    list(name = "a", alpha = matrix(0.5, 2, 3), lbf = matrix(1, 2, 3),
         sets_cs = list(L1 = 1L, L2 = 2L)),
    list(name = "b", alpha = matrix(0.5, 3, 3), lbf = matrix(1, 3, 3),
         sets_cs = list(L1 = 1L, L2 = 2L, L3 = 3L)))
  tuples <- enumerate_cs_tuples(v, by = "fit", cs_only = TRUE)
  expect_length(tuples, 6L)                 # 2 x 3 cross product
  expect_true(all(vapply(tuples, length, integer(1)) == 2L))
})

# ---- susiex_configurations() -----------------------------------------------

test_that("susiex_configurations errors when N exceeds the safety ceiling", {
  v <- list(
    list(name = "a", alpha = matrix(0.5, 1, 2), lbf = matrix(1, 1, 2),
         sets_cs = list(L1 = 1L)),
    list(name = "b", alpha = matrix(0.5, 1, 2), lbf = matrix(1, 1, 2),
         sets_cs = list(L1 = 1L)))
  expect_error(
    susiex_configurations(v, by = "fit", prob_thresh = 0.8, max_traits = 1L),
    "exceeds the safety ceiling")
  # At the ceiling it is allowed.
  expect_silent(
    susiex_configurations(v, by = "fit", prob_thresh = 0.8, max_traits = 2L))
})

test_that("susiex_configurations returns list() when no usable CS tuple exists", {
  # A view with zero CS indices, and a view whose only CS row is all-zero
  # alpha, both yield an empty result.
  v_no_idx <- list(
    list(name = "a", alpha = matrix(0.5, 2, 3), lbf = matrix(1, 2, 3),
         sets_cs = structure(list(), cs_idx = 99L)),   # 0 indices
    list(name = "b", alpha = matrix(0.5, 2, 3), lbf = matrix(1, 2, 3),
         sets_cs = list(L1 = 1L)))
  expect_equal(susiex_configurations(v_no_idx, by = "fit", prob_thresh = 0.8),
               list())

  v_zero_alpha <- list(
    list(name = "a", alpha = matrix(0, 1, 3), lbf = matrix(1, 1, 3),
         sets_cs = list(L1 = 1L)),
    list(name = "b", alpha = matrix(0.5, 1, 3), lbf = matrix(1, 1, 3),
         sets_cs = list(L1 = 1L)))
  expect_equal(susiex_configurations(v_zero_alpha, by = "fit",
                                     prob_thresh = 0.8), list())
})

test_that("susiex_configurations computes variant-level config / marginal probs correctly", {
  # Two traits, each with one CS row over two aligned variants.
  # logBF_trait and config_prob use the same variant-level LBF kernel.
  v <- list(
    list(name = "g1",
         alpha = matrix(c(1, 0), 1, 2,
                        dimnames = list(NULL, c("rs1", "rs2"))),
         lbf = matrix(c(2, 0), 1, 2,
                      dimnames = list(NULL, c("rs1", "rs2"))),
         sets_cs = list(L1 = 1L)),
    list(name = "g2",
         alpha = matrix(c(1, 0), 1, 2,
                        dimnames = list(NULL, c("rs1", "rs2"))),
         lbf = matrix(c(3, 0), 1, 2,
                      dimnames = list(NULL, c("rs1", "rs2"))),
         sets_cs = list(L1 = 1L)))
  res <- susiex_configurations(v, by = "fit", prob_thresh = 0.8)
  expect_length(res, 1L)
  tup <- res[[1]]
  raw <- attr(tup, "raw")

  lbf <- list(g1 = c(rs1 = 2, rs2 = 0), g2 = c(rs1 = 3, rs2 = 0))

  # logBF_trait is the singleton-activation log BF, named by trait.
  expect_equal(raw$logBF_trait,
               c(g1 = log(sum(exp(-log(2) + lbf$g1))),
                 g2 = log(sum(exp(-log(2) + lbf$g2)))),
               tolerance = 1e-8)
  expect_named(raw$cs_indices, c("g1", "g2"))

  # config_prob is a proper distribution over 2^2 = 4 configs.
  expect_length(raw$config_prob, 4L)
  expect_equal(sum(raw$config_prob), 1, tolerance = 1e-8)
  expect_length(raw$marginal_prob, 2L)
  expect_named(raw$marginal_prob, c("g1", "g2"))

  # Independent recomputation from the variant-level SuSiEx kernel over
  # aligned variants.
  logBF_conf <- apply(raw$configs, 1, function(cfg) {
    active <- which(cfg == 1L)
    if (length(active) == 0L) return(0)
    log(sum(exp(-log(2) + Reduce("+", lbf[active]))))
  })
  p <- exp(logBF_conf - max(logBF_conf)); p <- p / sum(p)
  expect_equal(unname(raw$config_prob), p, tolerance = 1e-8)
  expect_equal(unname(raw$marginal_prob),
               as.vector(crossprod(raw$configs, p)), tolerance = 1e-8)
  expect_named(tup$config_probability,
               c("g1", "g2", "config_prob"))
  expect_equal(tup$config_probability$config_prob, raw$config_prob,
               tolerance = 1e-8)
  expect_equal(rownames(tup$activation_summary),
               c("cs_indices", "logBF_trait", "posthoc_prob", "active"))
  expect_named(tup$activation_summary, c("g1", "g2"))
  expect_equal(unname(unlist(tup$activation_summary["cs_indices", ])),
               c("L1", "L1"))

  # active = marginal >= prob_thresh.
  expect_equal(raw$active, raw$marginal_prob >= 0.8)
})

test_that("susiex_configurations does not call disjoint strong signals shared", {
  variants <- c("rs1", "rs2", "rs3")
  make_view <- function(name, peak) {
    alpha <- matrix(0, 1, 3, dimnames = list(NULL, variants))
    lbf   <- matrix(-20, 1, 3, dimnames = list(NULL, variants))
    alpha[1, peak] <- 1
    lbf[1, peak]   <- 12
    list(name = name, alpha = alpha, lbf = lbf,
         sets_cs = list(L1 = match(peak, variants)))
  }
  v <- list(make_view("t1", "rs1"),
            make_view("t2", "rs2"),
            make_view("t3", "rs3"))

  tup <- susiex_configurations(v, by = "fit", prob_thresh = 0.8)[[1]]
  raw <- attr(tup, "raw")
  all_active <- which(rowSums(raw$configs) == 3L)

  expect_lt(raw$config_prob[all_active], 1e-8)
  expect_true(all(raw$marginal_prob < 0.8))
  expect_false(any(raw$active))
})

test_that("susiex_configurations aligns variant-level LBFs by column names", {
  alpha1 <- matrix(c(1, 0), 1, 2, dimnames = list(NULL, c("rs1", "rs2")))
  lbf1   <- matrix(c(8, 0), 1, 2, dimnames = list(NULL, c("rs1", "rs2")))

  alpha2_same <- matrix(c(1, 0), 1, 2, dimnames = list(NULL, c("rs1", "rs2")))
  lbf2_same   <- matrix(c(8, 0), 1, 2, dimnames = list(NULL, c("rs1", "rs2")))
  alpha2_rev  <- matrix(c(0, 1), 1, 2, dimnames = list(NULL, c("rs2", "rs1")))
  lbf2_rev    <- matrix(c(0, 8), 1, 2, dimnames = list(NULL, c("rs2", "rs1")))

  v_same <- list(
    list(name = "g1", alpha = alpha1, lbf = lbf1, sets_cs = list(L1 = 1L)),
    list(name = "g2", alpha = alpha2_same, lbf = lbf2_same,
         sets_cs = list(L1 = 1L)))
  v_rev <- list(
    list(name = "g1", alpha = alpha1, lbf = lbf1, sets_cs = list(L1 = 1L)),
    list(name = "g2", alpha = alpha2_rev, lbf = lbf2_rev,
         sets_cs = list(L1 = 1L)))

  same <- susiex_configurations(v_same, by = "fit", prob_thresh = 0.8)[[1]]
  rev  <- susiex_configurations(v_rev,  by = "fit", prob_thresh = 0.8)[[1]]

  expect_equal(attr(rev, "raw")$config_prob, attr(same, "raw")$config_prob,
               tolerance = 1e-8)
  expect_equal(attr(rev, "raw")$marginal_prob,
               attr(same, "raw")$marginal_prob, tolerance = 1e-8)
})

test_that("susiex_configurations rejects mismatched alpha / LBF names", {
  v <- list(
    list(name = "bad",
         alpha = matrix(c(1, 0), 1, 2,
                        dimnames = list(NULL, c("rs1", "rs2"))),
         lbf = matrix(c(2, 0), 1, 2,
                      dimnames = list(NULL, c("rs2", "rs1"))),
         sets_cs = list(L1 = 1L)),
    list(name = "ok",
         alpha = matrix(c(1, 0), 1, 2,
                        dimnames = list(NULL, c("rs1", "rs2"))),
         lbf = matrix(c(3, 0), 1, 2,
                      dimnames = list(NULL, c("rs1", "rs2"))),
         sets_cs = list(L1 = 1L)))

  expect_error(
    susiex_configurations(v, by = "fit", prob_thresh = 0.8),
    "column names of `alpha` and `lbf` must match")
})

test_that("susiex_configurations is reachable through the public entry point", {
  fits <- list(coloc_real_fit(11), coloc_real_fit(12))
  res  <- susie_post_outcome_configuration(fits, method = "susiex", by = "fit")
  expect_true(length(res$susiex) >= 1L)
  expect_equal(sum(attr(res$susiex[[1]], "raw")$config_prob), 1,
               tolerance = 1e-8)
})

# ---- .logsum() / .logdiff() / combine_abf_pair() ---------------------------

test_that(".logsum and .logdiff match closed-form values", {
  expect_equal(.logsum(c(0, 0, 0)), log(3), tolerance = 1e-8)
  expect_equal(.logsum(log(c(1, 2, 3))), log(6), tolerance = 1e-8)
  expect_equal(.logdiff(log(5), log(2)), log(3), tolerance = 1e-8)
})

test_that("combine_abf_pair returns five named PPs that sum to 1", {
  set.seed(101)
  J  <- 20
  l1 <- rnorm(J); l2 <- rnorm(J)
  pp <- combine_abf_pair(l1, l2, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  expect_named(pp, paste0("PP.H", 0:4))
  expect_true(all(is.finite(pp)))
  expect_equal(sum(pp), 1, tolerance = 1e-8)
})

test_that("combine_abf_pair gives high PP.H4 for a shared strong signal", {
  J  <- 30
  l1 <- rep(-5, J); l1[10] <- 8
  l2 <- rep(-5, J); l2[10] <- 8           # same SNP carries both signals
  pp <- combine_abf_pair(l1, l2, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  expect_gt(pp["PP.H4"], 0.9)
  expect_equal(unname(which.max(pp)), 5L)
})

test_that("combine_abf_pair errors on length-mismatched inputs", {
  expect_error(combine_abf_pair(c(1, 2, 3), c(1, 2),
                                p1 = 1e-4, p2 = 1e-4, p12 = 5e-6))
})

# ---- coloc_pairwise_abf() --------------------------------------------------

test_that("coloc_pairwise_abf returns the empty schema for N < 2", {
  # Single trait => no pairs => the documented 0-row schema.
  v <- list(list(name = "a", alpha = matrix(0.5, 1, 3),
                 lbf = matrix(1, 1, 3), sets_cs = list(L1 = 1L)))
  df <- coloc_pairwise_abf(v, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0L)
  expect_named(df, c("trait1", "trait2", "l1", "l2", "hit1", "hit2",
                     "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"))
})

test_that("coloc_pairwise_abf produces one row per eligible CS pair", {
  # Trait A has 2 CSs, trait B has 1 CS => 2 rows.
  v <- list(
    list(name = "A",
         alpha = matrix(c(0, 1, 0,  0, 0, 1), 2, 3, byrow = TRUE),
         lbf   = matrix(c(0, 5, 0,  0, 0, 5), 2, 3, byrow = TRUE),
         sets_cs = list(L1 = 1L, L2 = 2L)),
    list(name = "B",
         alpha = matrix(c(0, 1, 0), 1, 3),
         lbf   = matrix(c(0, 5, 0), 1, 3),
         sets_cs = list(L1 = 1L)))
  df <- coloc_pairwise_abf(v, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  expect_equal(nrow(df), 2L)
  expect_equal(df$trait1, c("A", "A"))
  expect_equal(df$trait2, c("B", "B"))
  expect_equal(df$l1, c(1L, 2L))
  # No colnames on lbf -> hits reported as snp_<idx>.
  expect_equal(df$hit1, c("snp_2", "snp_3"))
  expect_equal(df$hit2, c("snp_2", "snp_2"))
})

test_that("coloc_pairwise_abf uses lbf colnames for hit labels when present", {
  cn <- c("rsX", "rsY", "rsZ")
  v <- list(
    list(name = "A", alpha = matrix(c(0, 1, 0), 1, 3),
         lbf = matrix(c(0, 5, 0), 1, 3, dimnames = list(NULL, cn)),
         sets_cs = list(L1 = 1L)),
    list(name = "B", alpha = matrix(c(0, 0, 1), 1, 3),
         lbf = matrix(c(0, 0, 5), 1, 3, dimnames = list(NULL, cn)),
         sets_cs = list(L1 = 1L)))
  df <- coloc_pairwise_abf(v, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  expect_equal(nrow(df), 1L)
  expect_equal(df$hit1, "rsY")
  expect_equal(df$hit2, "rsZ")
})

test_that("coloc_pairwise_abf aligns variant sets by name before ABF", {
  cn_a <- c("rs1", "rs2", "rs3", "rs4")
  cn_b <- c("rs3", "rs2", "rs5")
  v <- list(
    list(name = "A", alpha = matrix(c(0, 1, 0, 0), 1, 4,
                                    dimnames = list(NULL, cn_a)),
         lbf = matrix(c(0, 5, 0, 1), 1, 4, dimnames = list(NULL, cn_a)),
         sets_cs = list(L1 = 1L)),
    list(name = "B", alpha = matrix(c(0, 1, 0), 1, 3,
                                    dimnames = list(NULL, cn_b)),
         lbf = matrix(c(0, 5, 0), 1, 3, dimnames = list(NULL, cn_b)),
         sets_cs = list(L1 = 1L)))
  df <- coloc_pairwise_abf(v, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  expected <- combine_abf_pair(c(5, 0), c(5, 0),
                               p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  expect_equal(nrow(df), 1L)
  expect_equal(df$hit1, "rs2")
  expect_equal(df$hit2, "rs2")
  expect_equal(unname(as.numeric(df[1, paste0("PP.H", 0:4)])),
               unname(expected))
})

test_that("coloc_pairwise_abf errors on unequal unnamed variant sets", {
  v <- list(
    list(name = "A", alpha = matrix(c(0, 1, 0), 1, 3),
         lbf = matrix(c(0, 5, 0), 1, 3),
         sets_cs = list(L1 = 1L)),
    list(name = "B", alpha = matrix(c(0, 1), 1, 2),
         lbf = matrix(c(0, 5), 1, 2),
         sets_cs = list(L1 = 1L)))
  expect_error(coloc_pairwise_abf(v, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6),
               "cannot align unnamed variant sets")
})

test_that("coloc_pairwise_abf skips all-zero alpha rows in either trait", {
  # Both traits' only CS row is all-zero alpha => no eligible pairs.
  v_both <- list(
    list(name = "A", alpha = matrix(0, 1, 3), lbf = matrix(1, 1, 3),
         sets_cs = list(L1 = 1L)),
    list(name = "B", alpha = matrix(0, 1, 3), lbf = matrix(1, 1, 3),
         sets_cs = list(L1 = 1L)))
  df_both <- coloc_pairwise_abf(v_both, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  expect_equal(nrow(df_both), 0L)
  expect_named(df_both, c("trait1", "trait2", "l1", "l2", "hit1", "hit2",
                          "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"))

  # Trait A has one valid CS; trait B's CS 2 points at an all-zero alpha row,
  # which is skipped, leaving only the (l1=1, l2=1) pair.
  v_second <- list(
    list(name = "A", alpha = rbind(c(0.6, 0.4)), lbf = rbind(c(2.0, 1.0)),
         sets_cs = structure(list(1L), names = "L1")),
    list(name = "B", alpha = rbind(c(0.7, 0.3), c(0, 0)),
         lbf = rbind(c(1.5, 0.5), c(0, 0)),
         sets_cs = structure(list(1L, 2L), names = c("L1", "L2"))))
  df_second <- coloc_pairwise_abf(v_second, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  expect_equal(nrow(df_second), 1L)
  expect_equal(df_second$l1, 1L)
  expect_equal(df_second$l2, 1L)
})

test_that("coloc_pairwise_abf returns empty df when a view has no CS indices", {
  # Trait B has only an out-of-range cs_idx -> no CS indices -> skipped.
  v <- list(
    list(name = "A", alpha = matrix(c(0, 1, 0), 1, 3),
         lbf = matrix(c(0, 5, 0), 1, 3), sets_cs = list(L1 = 1L)),
    list(name = "B", alpha = matrix(0.5, 2, 3), lbf = matrix(1, 2, 3),
         sets_cs = structure(list(), cs_idx = 99L)))
  df <- coloc_pairwise_abf(v, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  expect_equal(nrow(df), 0L)
})

test_that("coloc_pairwise_abf is reachable through the public entry point", {
  fits <- list(coloc_real_fit(21), coloc_real_fit(22))
  res  <- susie_post_outcome_configuration(fits, method = "coloc_pairwise",
                                           by = "fit")
  expect_s3_class(res$coloc_pairwise, "data.frame")
  expect_true(nrow(res$coloc_pairwise) >= 1L)
})

# ---- by = "outcome" end-to-end with a hand-built mvsusie fit ---------------

test_that("by = 'outcome' on an mvsusie runs SuSiEx per-outcome on the diagonal", {
  set.seed(31)
  L <- 2; J <- 5; R <- 3
  A   <- matrix(0, L, J); A[1, 2] <- 1; A[2, 4] <- 1
  LV  <- matrix(rnorm(L * J), L, J)
  arr <- array(rnorm(L * J * R), dim = c(L, J, R))
  dimnames(arr) <- list(NULL, NULL, c("o1", "o2", "o3"))
  fit <- coloc_mv_fit(A, LV, arr, cs = list(L1 = 1L, L2 = 2L))

  res <- susie_post_outcome_configuration(fit, method = "susiex",
                                          by = "outcome")
  # 3 outcome views => configs over 2^3 = 8 patterns; diagonal of 2 CSs => 2.
  expect_length(res$susiex, 2L)
  expect_named(res$susiex[[1]],
               c("mvsusie_cs_summary", "mvsusie_config_summary", "mvsusie_cs_variant_summary",
                 "susiex_config_probability", "susiex_activation_summary"))
  expect_equal(nrow(res$susiex[[1]]$susiex_config_probability), 8L)
  expect_named(res$susiex[[1]]$susiex_activation_summary,
               c("o1", "o2", "o3"))
  expect_equal(unname(attr(res$susiex[[1]], "raw")$cs_indices),
               c(1L, 1L, 1L))
  expect_equal(unname(attr(res$susiex[[2]], "raw")$cs_indices),
               c(2L, 2L, 2L))
})

# ---- Renderer color / NA branches ------------------------------------------

test_that(".print_susiex_table colors active / ambiguous / below / NA cells", {
  # Marginals spanning every band: active (>=0.8), ambiguous [0.5,0.8),
  # below (<0.5), and NA.
  tup <- coloc_tuple(c("a", "b", "c", "d"), c(1, 1, 1, 1),
                     c(0.95, 0.6, 0.2, NA))
  s   <- summary(coloc_post_obj(susiex = list(tup)), color = TRUE,
                 signal_only = FALSE)
  out <- capture.output(print(s))
  expect_true(any(grepl("\033\\[", out)))   # ANSI present
  expect_true(any(grepl("NA", out)))        # NA cell rendered
})

test_that(".print_coloc_table colors H0..H4 verdicts", {
  cdf <- rbind(
    coloc_row("A", "B", c(0.90, 0.025, 0.025, 0.025, 0.025)),   # H0
    coloc_row("A", "C", c(0.025, 0.90, 0.025, 0.025, 0.025)),   # H1
    coloc_row("A", "D", c(0.025, 0.025, 0.90, 0.025, 0.025)),   # H2
    coloc_row("B", "C", c(0.025, 0.025, 0.025, 0.90, 0.025)),   # H3
    coloc_row("B", "D", c(0.025, 0.025, 0.025, 0.025, 0.90)))   # H4
  s   <- summary(coloc_post_obj(coloc = cdf), color = TRUE,
                 signal_only = FALSE)
  out <- capture.output(print(s))
  expect_true(any(grepl("\033\\[", out)))
  expect_true(any(grepl("H4 shared", out)))
  expect_true(any(grepl("H3 distinct", out)))
  expect_true(any(grepl("H0 no signal", out)))
})

test_that(".print_coloc_table renders an NA verdict / NA pp cell", {
  # Force the NA-verdict / NA-pp print branches directly.
  cdf <- coloc_row("A", "B", c(NA, NA, NA, NA, NA))
  cdf$verdict <- NA_character_
  cdf$top_pp  <- NA_real_
  out <- capture.output(.print_coloc_table(cdf, use_color = TRUE))
  expect_true(any(grepl("NA", out)))
})

test_that(".print_susiex_table renders an NA top pattern", {
  # top_pattern is NA when configs is not a usable matrix; build the rendered
  # df directly to hit fmt_pat's NA branch.
  df <- data.frame(tuple = "(1,1)", a = 0.95, b = 0.95,
                   top_pattern = NA_character_, top_prob = 0.5,
                   stringsAsFactors = FALSE)
  out <- capture.output(
    .print_susiex_table(df, prob_thresh = 0.8, ambiguous_lower = 0.5,
                        use_color = TRUE))
  expect_true(any(grepl("NA", out)))
})

test_that(".print_aligned prints just the header when there are no rows", {
  out <- capture.output(.print_aligned(c("a", "b"), list()))
  expect_length(out, 1L)
  expect_true(grepl("a", out) && grepl("b", out))
})

# ---- Residual summary-internal branches (engine-adjacent) ------------------

test_that(".summarise_susiex returns NULL df when no trait names are present", {
  # A tuple with no marginal_prob => trait_names_all is empty.
  tup <- list(cs_indices = 1L, config_prob = 1, configs = matrix(0L, 1, 1))
  res <- .summarise_susiex(list(tup), signal_only = FALSE, prob_thresh = 0.8)
  expect_null(res$df)
  expect_equal(res$n_total, 1L)
  expect_equal(res$n_kept, 0L)
})

test_that(".summarise_susiex skips tuples with non-finite config_prob", {
  tup <- list(cs_indices = setNames(1L, "a"),
              config_prob = c(NA, 1),
              configs = matrix(c(0L, 1L), 2, 1),
              marginal_prob = setNames(0.95, "a"))
  res <- .summarise_susiex(list(tup), signal_only = FALSE, prob_thresh = 0.8)
  expect_null(res$df)
  expect_equal(res$n_total, 1L)
})

test_that(".summarise_susiex returns NULL df when signal_only filters everything", {
  tup <- coloc_tuple("a", 1, 0.2)     # below threshold
  res <- .summarise_susiex(list(tup), signal_only = TRUE, prob_thresh = 0.8)
  expect_null(res$df)
  expect_equal(res$n_total, 1L)
  expect_equal(res$n_kept, 0L)
})

test_that(".summarise_coloc returns NULL df when every row is H0 under signal_only", {
  cdf <- coloc_row("A", "B", c(0.99, 0.0025, 0.0025, 0.0025, 0.0025))
  res <- .summarise_coloc(cdf, signal_only = TRUE)
  expect_null(res$df)
  expect_equal(res$n_total, 1L)
  expect_equal(res$n_kept, 0L)
})

test_that("print.summary footers hidden susiex CS tuples", {
  tups <- list(coloc_tuple(c("a", "b"), c(1, 1), c(0.95, 0.2)),  # kept (active)
               coloc_tuple(c("a", "b"), c(2, 2), c(0.3, 0.2)))   # hidden
  s <- summary(coloc_post_obj(susiex = tups), color = FALSE,
               signal_only = TRUE)
  out <- capture.output(print(s))
  expect_true(any(grepl("1/2 CS tuples hidden", out)))
})
