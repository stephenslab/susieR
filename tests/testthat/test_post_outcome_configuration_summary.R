context("susie_post_outcome_configuration summary")

# Tests for summary.susie_post_outcome_configuration and its print method.
# The numerical kernels (susiex / coloc_pairwise) are exercised in
# test_post_outcome_configuration.R; here we focus on dispatch via the class
# tag, tidy-table shape and column names, signal_only filtering with the
# kept/total bookkeeping, and defensive handling of malformed / partial
# input. Shared tuple / object / data.frame builders live in helper_coloc.R
# (coloc_tuple, coloc_post_obj, coloc_df).

# ---- Dispatch --------------------------------------------------------------

test_that("summary() dispatches on the class tag", {
  obj <- coloc_post_obj(
    susiex = list(coloc_tuple(c("a", "b"),
                              cs_indices    = c(1, 1),
                              marginal_prob = c(0.95, 0.95))))
  s <- summary(obj, color = FALSE)
  expect_s3_class(s, "summary.susie_post_outcome_configuration")
  # Print returns its input invisibly without erroring.
  out <- capture.output(p <- print(s))
  expect_identical(p, s)
  expect_true(any(grepl("SuSiEx:", out)))
})

# ---- Tidy table shape ------------------------------------------------------

test_that("susiex tidy table carries one row per CS tuple with reserved + per-trait columns", {
  tuples <- list(
    coloc_tuple(c("trait_a", "trait_b"),
                cs_indices    = c(1, 1),
                marginal_prob = c(0.95, 0.30)),
    coloc_tuple(c("trait_a", "trait_b"),
                cs_indices    = c(2, 2),
                marginal_prob = c(0.10, 0.92)))
  s <- summary(coloc_post_obj(susiex = tuples), color = FALSE,
               signal_only = FALSE)
  expect_s3_class(s$susiex, "data.frame")
  expect_equal(nrow(s$susiex), 2L)
  expect_setequal(colnames(s$susiex),
                  c("tuple", "trait_a", "trait_b", "top_pattern", "top_prob"))
  expect_equal(s$susiex$tuple, c("(1,1)", "(2,2)"))
  expect_equal(s$susiex$trait_a, c(0.95, 0.10))
})

test_that("single-trait susiex tuple yields a one-trait tidy row", {
  s <- summary(coloc_post_obj(
         susiex = list(coloc_tuple("solo", cs_indices = 1,
                                   marginal_prob = 0.95))),
       color = FALSE, signal_only = FALSE)
  expect_equal(nrow(s$susiex), 1L)
  expect_setequal(colnames(s$susiex),
                  c("tuple", "solo", "top_pattern", "top_prob"))
  expect_equal(s$susiex$tuple, "(1)")
  expect_equal(s$susiex$solo, 0.95)
})

test_that("coloc tidy table extends the input data.frame with verdict and top_pp", {
  rows <- list(
    list(t1 = "A", t2 = "B", l1 = 1, l2 = 1, h1 = "rs1", h2 = "rs1",
         pp = c(0.001, 0.001, 0.001, 0.05, 0.947)),
    list(t1 = "A", t2 = "C", l1 = 1, l2 = 1, h1 = "rs1", h2 = "rs9",
         pp = c(0.99, 0.005, 0.002, 0.002, 0.001)))   # H0 dominant
  s <- summary(coloc_post_obj(coloc = coloc_df(rows)),
               color = FALSE, signal_only = FALSE)
  expect_s3_class(s$coloc_pairwise, "data.frame")
  expect_equal(nrow(s$coloc_pairwise), 2L)
  expect_true(all(c("verdict", "top_pp") %in% colnames(s$coloc_pairwise)))
  expect_equal(s$coloc_pairwise$verdict, c("H4 shared", "H0 no signal"))
  expect_equal(s$coloc_pairwise$top_pp, c(0.947, 0.990), tolerance = 1e-9)
})

# ---- signal_only filtering + bookkeeping -----------------------------------

test_that("signal_only drops below-threshold susiex rows and counts them in n_total/n_kept", {
  tuples <- list(
    coloc_tuple(c("a", "b"),
                cs_indices    = c(1, 1),
                marginal_prob = c(0.95, 0.20)),    # signal (a active)
    coloc_tuple(c("a", "b"),
                cs_indices    = c(2, 2),
                marginal_prob = c(0.40, 0.30)))    # no signal
  s_filt <- summary(coloc_post_obj(susiex = tuples), color = FALSE,
                    signal_only = TRUE)
  s_all  <- summary(coloc_post_obj(susiex = tuples), color = FALSE,
                    signal_only = FALSE)
  expect_equal(s_filt$susiex_n_total, 2L)
  expect_equal(s_filt$susiex_n_kept,  1L)
  expect_equal(s_all$susiex_n_kept,   2L)
  expect_equal(s_filt$susiex$tuple, "(1,1)")
})

test_that("signal_only drops H0-dominant coloc rows and footers the count", {
  rows <- list(
    list(t1 = "A", t2 = "B", l1 = 1, l2 = 1, h1 = "rs1", h2 = "rs1",
         pp = c(0.99, 0.005, 0.002, 0.002, 0.001)),    # H0
    list(t1 = "A", t2 = "B", l1 = 1, l2 = 2, h1 = "rs1", h2 = "rs7",
         pp = c(0.001, 0.001, 0.001, 0.05, 0.947)))    # H4
  obj <- coloc_post_obj(coloc = coloc_df(rows))
  s   <- summary(obj, color = FALSE, signal_only = TRUE)
  expect_equal(s$coloc_n_total, 2L)
  expect_equal(s$coloc_n_kept,  1L)
  out <- capture.output(print(s))
  expect_true(any(grepl("1/2 pairs hidden", out)))
})

# ---- Defensive paths -------------------------------------------------------

test_that("summary handles entirely empty input gracefully", {
  s <- summary(coloc_post_obj(), color = FALSE)
  expect_null(s$susiex)
  expect_null(s$coloc_pairwise)
  out <- capture.output(print(s))
  expect_true(any(grepl("no signals", out)))
})

test_that("summary tolerates susiex tuples with missing fields (skips them)", {
  good   <- coloc_tuple(c("a", "b"),
                        cs_indices    = c(1, 1),
                        marginal_prob = c(0.95, 0.95))
  broken <- list(cs_indices = c(2, 2))   # missing marginal_prob etc.
  s <- summary(coloc_post_obj(susiex = list(good, broken)),
               color = FALSE, signal_only = FALSE)
  expect_equal(nrow(s$susiex), 1L)
  expect_equal(s$susiex_n_total, 2L)   # both counted in total
  expect_equal(s$susiex_n_kept,  1L)
})

test_that("summary prefixes trait names that collide with reserved column names", {
  # Trait literally named "tuple" must not clobber the CS-tuple column.
  tup <- coloc_tuple(c("tuple", "top_prob"),
                     cs_indices    = c(1, 1),
                     marginal_prob = c(0.95, 0.95))
  s <- summary(coloc_post_obj(susiex = list(tup)), color = FALSE)
  expect_true("tuple" %in% colnames(s$susiex))
  # Reserved-name traits get a "trait_" prefix.
  expect_true("trait_tuple"    %in% colnames(s$susiex))
  expect_true("trait_top_prob" %in% colnames(s$susiex))
})

test_that("summary warns and skips coloc when required PP columns are missing", {
  bad <- data.frame(trait1 = "A", trait2 = "B",
                    l1 = 1, l2 = 1, hit1 = "rs1", hit2 = "rs1",
                    PP.H0 = 0.5, PP.H1 = 0.5,
                    stringsAsFactors = FALSE)
  expect_warning(
    s <- summary(coloc_post_obj(coloc = bad), color = FALSE),
    "missing required columns")
  expect_null(s$coloc_pairwise)
  expect_equal(s$coloc_n_total, 1L)
  expect_equal(s$coloc_n_kept,  0L)
})

test_that("summary validates its arguments", {
  obj <- coloc_post_obj()
  expect_error(summary(obj, prob_thresh = 1.1),     "prob_thresh")
  expect_error(summary(obj, prob_thresh = -0.1),    "prob_thresh")
  expect_error(summary(obj, ambiguous_lower = 0.9, prob_thresh = 0.8),
               "ambiguous_lower")
  expect_error(summary(obj, signal_only = NA),      "signal_only")
  expect_error(summary(obj, color = "yes"),         "color")
})

# ---- Color toggle ----------------------------------------------------------

test_that("color toggle controls ANSI escape sequences in output", {
  obj <- coloc_post_obj(
    susiex = list(coloc_tuple(c("a", "b"),
                              cs_indices    = c(1, 1),
                              marginal_prob = c(0.95, 0.95))))
  out_plain <- capture.output(print(summary(obj, color = FALSE)))
  out_color <- capture.output(print(summary(obj, color = TRUE)))
  expect_false(any(grepl("\033\\[", out_plain)))
  expect_true(any(grepl("\033\\[", out_color)),
              info = "Forcing color = TRUE should emit at least one ANSI SGR.")
})
