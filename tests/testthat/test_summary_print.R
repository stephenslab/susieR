context("summary and print S3 methods")

# ---- summary.susie structure ----

test_that("summary.susie returns a summary.susie list with 'vars' and 'cs'", {
  set.seed(1)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  expect_type(summ, "list")
  expect_s3_class(summ, "summary.susie")
  expect_named(summ, c("vars", "cs"))
})

test_that("summary.susie vars data frame has variable/variable_prob/cs columns with nrow == p", {
  set.seed(2)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  expect_s3_class(summ$vars, "data.frame")
  expect_true(all(c("variable", "variable_prob", "cs") %in% colnames(summ$vars)))
  expect_equal(nrow(summ$vars), dat$p)
})

test_that("summary.susie cs data frame has correct columns and one row per CS", {
  set.seed(3)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  # signal_sd=2 yields CS deterministically at this seed
  expect_s3_class(summ$cs, "data.frame")
  expect_true(all(c("cs", "cs_log10bf", "cs_avg_r2", "cs_min_r2", "variable") %in%
                    colnames(summ$cs)))
  expect_equal(nrow(summ$cs), length(fit$sets$cs))
})

test_that("summary.susie vars are sorted by PIP descending when CS exist", {
  set.seed(4)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  expect_true(all(diff(summ$vars$variable_prob) <= 0))
})

test_that("summary.susie cs column correctly maps variables to their credible set index", {
  set.seed(5)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  for (i in seq_along(fit$sets$cs)) {
    cs_vars <- fit$sets$cs[[i]]
    cs_idx  <- fit$sets$cs_index[i]
    expect_true(all(summ$vars$cs[summ$vars$variable %in% cs_vars] == cs_idx))
  }
})

test_that("summary.susie vars nrow equals p when null_weight > 0", {
  set.seed(6)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, null_weight = 0.1, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  expect_equal(nrow(summ$vars), dat$p)
})

test_that("summary.susie errors when fit$sets is NULL", {
  set.seed(7)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$sets <- NULL
  expect_error(summary(fit), "credible set information")
})

test_that("summary.susie returns NULL cs when no credible sets pass purity filter", {
  set.seed(8)
  dat <- simulate_regression(n = 100, p = 50, k = 3, signal_sd = 0.1)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95, min_abs_corr = 0.99)
  summ <- summary(fit)
  expect_null(summ$cs)
  expect_s3_class(summ$vars, "data.frame")
})

test_that("summary.susie cs_log10bf equals fit$lbf[cs_idx] / log(10)", {
  set.seed(9)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  for (i in seq_len(nrow(summ$cs))) {
    cs_idx <- summ$cs$cs[i]
    expect_equal(summ$cs$cs_log10bf[i], fit$lbf[cs_idx] / log(10), tolerance = 1e-8)
  }
})

test_that("summary.susie cs_avg_r2 and cs_min_r2 equal purity correlations squared", {
  set.seed(10)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  for (i in seq_len(nrow(summ$cs))) {
    expect_equal(summ$cs$cs_avg_r2[i], fit$sets$purity$mean.abs.corr[i]^2, tolerance = 1e-8)
    expect_equal(summ$cs$cs_min_r2[i], fit$sets$purity$min.abs.corr[i]^2, tolerance = 1e-8)
  }
})

test_that("summary.susie variable column contains all indices 1:p", {
  set.seed(15)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  expect_equal(sort(unique(summ$vars$variable)), 1:dat$p)
})

test_that("summary.susie variable_prob matches pip values", {
  set.seed(16)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  expect_equal(sort(summ$vars$variable_prob, decreasing = TRUE),
               sort(fit$pip, decreasing = TRUE),
               tolerance = 1e-8)
})

# ---- print.summary.susie ----

test_that("print.summary.susie emits 'Variables in credible sets' and 'Credible sets summary' messages", {
  set.seed(11)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  expect_message(capture.output(print(summ)), "Variables in credible sets")
  expect_message(capture.output(print(summ)), "Credible sets summary")
})

test_that("print.summary.susie still emits messages when no CS exist", {
  set.seed(12)
  dat <- simulate_regression(n = 100, p = 50, k = 3, signal_sd = 0.1)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95, min_abs_corr = 0.99)
  summ <- summary(fit)
  expect_message(capture.output(print(summ)), "Variables in credible sets")
  expect_message(capture.output(print(summ)), "Credible sets summary")
})

test_that("print.summary.susie output includes variable indices for CS members", {
  set.seed(13)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  summ <- summary(fit)
  output_text <- paste(capture.output(suppressMessages(print(summ))), collapse = "\n")
  cs_vars <- summ$vars[summ$vars$cs > 0, ]
  for (i in seq_len(nrow(cs_vars))) {
    expect_true(grepl(as.character(cs_vars$variable[i]), output_text))
  }
})

test_that("summary then print together produce messages", {
  set.seed(14)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)
  expect_message(capture.output({
    summ <- summary(fit)
    print(summ)
  }))
})
