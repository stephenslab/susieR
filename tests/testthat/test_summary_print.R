context("summary and print S3 methods")

# =============================================================================
# summary.susie - Summary Statistics
# =============================================================================

test_that("summary.susie creates correct structure", {
  set.seed(1)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  expect_type(summ, "list")
  expect_s3_class(summ, "summary.susie")
  expect_named(summ, c("vars", "cs"))
})

test_that("summary.susie variables data frame has correct structure", {
  set.seed(2)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  expect_s3_class(summ$vars, "data.frame")
  expect_true(all(c("variable", "variable_prob", "cs") %in% colnames(summ$vars)))
  expect_equal(nrow(summ$vars), dat$p)
})

test_that("summary.susie CS data frame has correct structure when CS exist", {
  set.seed(3)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  if (!is.null(fit$sets$cs)) {
    expect_s3_class(summ$cs, "data.frame")
    expect_true(all(c("cs", "cs_log10bf", "cs_avg_r2", "cs_min_r2", "variable") %in%
                    colnames(summ$cs)))
    expect_equal(nrow(summ$cs), length(fit$sets$cs))
  }
})

test_that("summary.susie variables sorted by PIP descending", {
  set.seed(4)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  if (!is.null(fit$sets$cs)) {
    expect_true(all(diff(summ$vars$variable_prob) <= 0))
  }
})

test_that("summary.susie cs column maps variables to credible sets", {
  set.seed(5)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  if (!is.null(fit$sets$cs)) {
    for (i in 1:length(fit$sets$cs)) {
      cs_vars <- fit$sets$cs[[i]]
      cs_idx <- fit$sets$cs_index[i]
      expect_true(all(summ$vars$cs[summ$vars$variable %in% cs_vars] == cs_idx))
    }
  }
})

test_that("summary.susie handles null_index correctly", {
  set.seed(6)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, null_weight = 0.1, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  if (!is.null(fit$null_index) && fit$null_index > 0) {
    expect_equal(nrow(summ$vars), dat$p)
  } else {
    expect_equal(nrow(summ$vars), ncol(fit$alpha))
  }
})

test_that("summary.susie errors when sets is NULL", {
  set.seed(7)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  fit$sets <- NULL

  expect_error(
    summary(fit),
    "credible set information"
  )
})

test_that("summary.susie handles no credible sets", {
  set.seed(8)
  dat <- simulate_regression(n = 100, p = 50, k = 3, signal_sd = 0.1)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95, min_abs_corr = 0.99)

  summ <- summary(fit)

  expect_null(summ$cs)
  expect_s3_class(summ$vars, "data.frame")
})

test_that("summary.susie log10BF calculation is correct", {
  set.seed(9)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  if (!is.null(summ$cs)) {
    for (i in 1:nrow(summ$cs)) {
      cs_idx <- summ$cs$cs[i]
      expected_log10bf <- fit$lbf[cs_idx] / log(10)
      expect_equal(summ$cs$cs_log10bf[i], expected_log10bf)
    }
  }
})

test_that("summary.susie r2 calculations are correct", {
  set.seed(10)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  if (!is.null(summ$cs) && !is.null(fit$sets$purity)) {
    for (i in 1:nrow(summ$cs)) {
      expect_equal(summ$cs$cs_avg_r2[i], fit$sets$purity$mean.abs.corr[i]^2)
      expect_equal(summ$cs$cs_min_r2[i], fit$sets$purity$min.abs.corr[i]^2)
    }
  }
})

# =============================================================================
# print.summary.susie - Console Output
# =============================================================================

test_that("print.summary.susie produces output", {
  set.seed(11)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  expect_message(capture.output(print(summ)), "Variables in credible sets")
  expect_message(capture.output(print(summ)), "Credible sets summary")
})

test_that("print.summary.susie handles no CS", {
  set.seed(12)
  dat <- simulate_regression(n = 100, p = 50, k = 3, signal_sd = 0.1)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95, min_abs_corr = 0.99)

  summ <- summary(fit)

  expect_message(capture.output(print(summ)))
})

test_that("print.summary.susie shows variables in CS", {
  set.seed(13)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  output <- capture.output(print(summ))
  output_text <- paste(output, collapse = "\n")

  if (!is.null(fit$sets$cs)) {
    cs_vars <- summ$vars[summ$vars$cs > 0, ]
    for (i in 1:nrow(cs_vars)) {
      expect_true(grepl(as.character(cs_vars$variable[i]), output_text))
    }
  }
})

test_that("summary and print work together", {
  set.seed(14)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  expect_message(capture.output({
    summ <- summary(fit)
    print(summ)
  }))
})

test_that("summary.susie variable column is sequential indices", {
  set.seed(15)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  expect_equal(sort(unique(summ$vars$variable)), 1:dat$p)
})

test_that("summary.susie variable_prob matches PIP", {
  set.seed(16)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  summ <- summary(fit)

  expect_equal(sort(summ$vars$variable_prob, decreasing = TRUE),
               sort(fit$pip, decreasing = TRUE))
})
