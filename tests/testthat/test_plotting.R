context("Plotting functions")

# ---- susie_plot basic functionality ----

test_that("susie_plot produces PIP, z, z_original, and log10PIP plots without error", {
  set.seed(1)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit_no_z  <- susie(dat$X, dat$y, L = 5, compute_univariate_zscore = FALSE,
                     verbose = FALSE)
  fit_with_z <- susie(dat$X, dat$y, L = 5, compute_univariate_zscore = TRUE,
                      verbose = FALSE)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)

  # PIP and log10PIP work without z-scores
  expect_error(susie_plot(fit_no_z, "PIP"),       NA)
  expect_error(susie_plot(fit_no_z, "log10PIP"),  NA)

  # z and z_original require compute_univariate_zscore = TRUE
  expect_error(susie_plot(fit_no_z, "z"),         "z-scores are not available")
  expect_error(susie_plot(fit_no_z, "z_original"),"z-scores are not available")

  # z and z_original succeed when z-scores are present
  expect_error(susie_plot(fit_with_z, "z"),        NA)
  expect_error(susie_plot(fit_with_z, "z_original"), NA)

  # invalid y-type errors for susie objects
  expect_error(susie_plot(fit_no_z, "invalid_type"), "Need to specify")
})

test_that("susie_plot returns NULL invisibly", {
  set.seed(2)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)
  result <- susie_plot(fit, "PIP")
  expect_null(result)
})

# ---- susie_plot pos argument validation ----

test_that("susie_plot errors when pos list is missing required elements", {
  set.seed(3)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$genomic_position <- 1000 + seq_along(fit$pip)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)

  expect_error(
    susie_plot(fit, "PIP", pos = list(start = 1000, end = 1025)),
    "pos argument should be a list"
  )
  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "genomic_position", end = 1025)),
    "pos argument should be a list"
  )
  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "genomic_position", start = 1000)),
    "pos argument should be a list"
  )
})

test_that("susie_plot errors when pos$attr is not in model", {
  set.seed(4)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)
  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "nonexistent_attr", start = 1, end = 25)),
    "Cannot find attribute nonexistent_attr"
  )
})

test_that("susie_plot errors when pos$start >= pos$end", {
  set.seed(5)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$genomic_position <- 1000 + seq_along(fit$pip)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)

  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "genomic_position", start = 1025, end = 1000)),
    "Position start should be smaller than end"
  )
  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "genomic_position", start = 1000, end = 1000)),
    "Position start should be smaller than end"
  )
})

test_that("susie_plot errors when numeric pos is outside variable range", {
  set.seed(6)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)

  expect_error(
    susie_plot(fit, "PIP", pos = 1:100),        # only 50 variables
    "Provided position is outside the range"
  )
  expect_error(
    susie_plot(fit, "PIP", pos = c(0, 1, 2)),   # 0 is out of range
    "Provided position is outside the range"
  )
})

# ---- susie_plot plotting parameters ----

test_that("susie_plot optional parameters work: add_bar, numeric pos, list pos, b", {
  set.seed(7)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$genomic_position <- 1000 + seq_along(fit$pip)

  b_test <- rep(0, 50); b_test[c(10, 20, 30)] <- 1

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)

  expect_error(susie_plot(fit, "PIP", add_bar = TRUE),           NA)
  expect_error(susie_plot(fit, "PIP", pos = 1:25),               NA)
  expect_error(susie_plot(fit, "PIP",
    pos = list(attr = "genomic_position", start = 1000, end = 1025)), NA)
  expect_error(susie_plot(fit, "PIP", b = b_test),               NA)
  expect_error(susie_plot(fit, "PIP", b = dat$beta, add_bar = TRUE, add_legend = TRUE), NA)
  expect_error(susie_plot(fit, "PIP",
    main = "Custom Title", col = "blue", cex = 0.5),             NA)
  expect_error(susie_plot(fit, "PIP", xlim = c(0, 30), ylim = c(0, 1)), NA)
})

test_that("susie_plot add_legend works for all valid positions and invalid falls back to topright", {
  set.seed(8)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)

  positions <- c(TRUE, "bottomright", "bottom", "bottomleft", "left",
                 "center", "right", "topleft", "top", "topright")
  for (pos in positions) {
    expect_error(susie_plot(fit, "PIP", add_legend = pos), NA,
                 info = paste("legend position:", pos))
  }
  # invalid position silently falls back to "topright"
  expect_error(susie_plot(fit, "PIP", add_legend = "invalid_position"), NA)
})

test_that("susie_plot max_cs filters CS by size and by purity", {
  set.seed(9)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)

  # max_cs as size limit
  expect_error(susie_plot(fit, "PIP", max_cs = 2),   NA)
  expect_error(susie_plot(fit, "PIP", max_cs = 400),  NA)

  # max_cs as purity threshold (0 < max_cs < 1)
  expect_error(susie_plot(fit, "PIP", max_cs = 0.1,   add_legend = TRUE), NA)
  expect_error(susie_plot(fit, "PIP", max_cs = 0.5,   add_legend = TRUE), NA)
  expect_error(susie_plot(fit, "PIP", max_cs = 0.999, add_legend = TRUE), NA)

  # max_cs that excludes all CS (size < 1) → NULL x0/y1 path
  expect_error(susie_plot(fit, "PIP", max_cs = 1,     add_legend = TRUE), NA)
})

# ---- susie_plot vector input ----

test_that("susie_plot works with plain numeric vectors for all y types", {
  set.seed(10)
  pip_vec  <- runif(50)
  z_vec    <- rnorm(50)
  data_vec <- runif(50, 0, 10)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)

  expect_error(susie_plot(pip_vec,  "PIP"),        NA)
  expect_error(susie_plot(z_vec,    "z"),           NA)
  expect_error(susie_plot(z_vec,    "z_original"),  NA)
  expect_error(susie_plot(pip_vec,  "log10PIP"),    NA)
  expect_error(susie_plot(data_vec, "custom_data"), NA)
})

# ---- susie_plot edge cases ----

test_that("susie_plot works with a fit that has no credible sets", {
  set.seed(11)
  dat <- simulate_regression(n = 100, p = 50, k = 0)   # no signal
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)
  expect_error(susie_plot(fit, "PIP"), NA)
})

test_that("susie_plot works for a single-variable fit", {
  set.seed(12)
  dat <- simulate_regression(n = 100, p = 1, k = 1)
  fit <- susie(dat$X, dat$y, L = 1, verbose = FALSE)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)
  expect_error(susie_plot(fit, "PIP"), NA)
})

test_that("susie_plot with list pos and credible sets adjusts to new positions", {
  set.seed(13)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$genomic_position <- 1000 + seq_along(fit$pip)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)
  expect_error(
    susie_plot(fit, "PIP",
      pos = list(attr = "genomic_position", start = 1000, end = 1050),
      add_legend = TRUE),
    NA
  )
})

test_that("susie_plot PIP and log10PIP fall back to alpha[1,] when pip is NULL", {
  set.seed(14)
  p <- 20
  m <- structure(
    list(alpha = matrix(runif(p) |> (\(x) x / sum(x))(), 1, p),
         mu    = matrix(0, 1, p),
         pip   = NULL,
         sets  = NULL,
         z     = NULL),
    class = "susie"
  )

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)
  expect_error(susie_plot(m, "PIP"),      NA)
  expect_error(susie_plot(m, "log10PIP"), NA)
})

test_that("susie_plot uses cs_index when available and cs_idx as fallback", {
  set.seed(15)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  expect_true(!is.null(fit$sets$cs_index))

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)

  # cs_index present branch
  expect_error(susie_plot(fit, "PIP", add_legend = TRUE), NA)

  # cs_index absent → fallback to cs_idx
  fit_no_index <- fit
  fit_no_index$sets$cs_index <- NULL
  expect_error(susie_plot(fit_no_index, "PIP", add_legend = TRUE), NA)
})

# ---- susie_plot integration with susie_ss and susie_rss ----

test_that("susie_plot works with susie_ss and susie_rss outputs", {
  set.seed(16)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  # susie_ss
  ss_stats <- compute_summary_stats(dat$X, dat$y)
  fit_ss <- susie_ss(ss_stats$XtX, ss_stats$Xty, ss_stats$yty,
                     n = ss_stats$n, L = 5, verbose = FALSE)

  # susie_rss
  dat2 <- simulate_regression(n = 200, p = 100, k = 3)
  ur <- univariate_regression(dat2$X, dat2$y)
  z  <- with(ur, betahat / sebetahat)
  R  <- cor(dat2$X)
  fit_rss <- susie_rss(z, R, n = 200, L = 5, verbose = FALSE)

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)
  expect_error(susie_plot(fit_ss,  "PIP"), NA)
  expect_error(susie_plot(fit_rss, "PIP"), NA)
})

# ---- susie_plot_iteration ----

test_that("susie_plot_iteration uses tempdir when file_prefix is missing", {
  set.seed(17)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, track_fit = FALSE, verbose = FALSE)

  expected_path <- file.path(tempdir(), "susie_plot.pdf")
  if (file.exists(expected_path)) file.remove(expected_path)

  invisible(capture.output(
    suppressMessages(susie_plot_iteration(fit, L = 5)),
    type = "output"
  ))

  expect_true(file.exists(expected_path))
  if (file.exists(expected_path)) file.remove(expected_path)
})

test_that("susie_plot_iteration returns NULL invisibly", {
  set.seed(18)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, track_fit = FALSE, verbose = FALSE)

  temp_prefix <- tempfile("susie_invisible_")
  on.exit({
    f <- list.files(dirname(temp_prefix), pattern = basename(temp_prefix), full.names = TRUE)
    if (length(f) > 0) file.remove(f)
  })

  invisible(capture.output({
    result <- suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix))
  }, type = "output"))

  expect_null(result)
})

test_that("susie_plot_iteration with track_fit=FALSE plots only the final iteration", {
  set.seed(19)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, track_fit = FALSE, verbose = FALSE)

  temp_prefix <- tempfile("susie_notrack_")
  on.exit({
    f <- list.files(dirname(temp_prefix), pattern = basename(temp_prefix), full.names = TRUE)
    if (length(f) > 0) file.remove(f)
  })

  expect_error({
    invisible(capture.output(
      suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix)),
      type = "output"
    ))
  }, NA)
})

test_that("susie_plot_iteration with track_fit=TRUE creates PDF files", {
  set.seed(20)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- suppressWarnings(susie(dat$X, dat$y, L = 5, track_fit = TRUE, max_iter = 5, verbose = FALSE))

  temp_dir    <- tempdir()
  temp_prefix <- file.path(temp_dir, "test_susie_iter")
  on.exit({
    f <- list.files(temp_dir, pattern = "test_susie_iter", full.names = TRUE)
    if (length(f) > 0) file.remove(f)
  })

  invisible(capture.output(
    suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix)),
    type = "output"
  ))

  pdf_files <- list.files(temp_dir, pattern = "test_susie_iter.*\\.pdf$", full.names = TRUE)
  expect_true(length(pdf_files) > 0)
})

test_that("susie_plot_iteration with pos subset works", {
  set.seed(21)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, track_fit = TRUE, max_iter = 10, verbose = FALSE)

  temp_prefix <- tempfile("susie_pos_")
  on.exit({
    f <- list.files(dirname(temp_prefix), pattern = basename(temp_prefix), full.names = TRUE)
    if (length(f) > 0) file.remove(f)
  })

  expect_error({
    invisible(capture.output(
      suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix, pos = 1:25)),
      type = "output"
    ))
  }, NA)
})

test_that("susie_plot_iteration with L greater than nrow(alpha) uses available rows", {
  set.seed(22)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 3, track_fit = TRUE, verbose = FALSE)

  temp_prefix <- tempfile("susie_largeL_")
  on.exit({
    f <- list.files(dirname(temp_prefix), pattern = basename(temp_prefix), full.names = TRUE)
    if (length(f) > 0) file.remove(f)
  })

  expect_error({
    invisible(capture.output(
      suppressMessages(susie_plot_iteration(fit, L = 10, file_prefix = temp_prefix)),
      type = "output"
    ))
  }, NA)
})

test_that("susie_plot_iteration defaults L to nrow(alpha) when L is not supplied", {
  set.seed(23)
  dat <- simulate_regression(n = 80, p = 40, k = 2)
  fit <- susie(dat$X, dat$y, L = 4, track_fit = FALSE, max_iter = 5, verbose = FALSE)

  temp_prefix <- tempfile("susie_noL_")
  on.exit({
    f <- list.files(dirname(temp_prefix), pattern = basename(temp_prefix), full.names = TRUE)
    if (length(f) > 0) file.remove(f)
  })

  expect_no_error(suppressMessages(
    susie_plot_iteration(fit, file_prefix = temp_prefix)
  ))
})

test_that("susie_plot_iteration with track_fit=TRUE creates a GIF animation", {
  # GIF creation shells out to ImageMagick's `convert`; skip where it's absent.
  skip_if(Sys.which("convert") == "", "ImageMagick 'convert' not on PATH")
  set.seed(24)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, track_fit = TRUE, max_iter = 10, verbose = FALSE)

  temp_prefix <- tempfile("susie_gif_")
  on.exit({
    f <- list.files(dirname(temp_prefix), pattern = basename(temp_prefix), full.names = TRUE)
    if (length(f) > 0) file.remove(f)
  })

  expect_error({
    invisible(capture.output(
      suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix)),
      type = "output"
    ))
  }, NA)
})

test_that("susie_plot_iteration accepts a susie_track object directly", {
  # GIF creation shells out to ImageMagick's `convert`; skip where it's absent.
  skip_if(Sys.which("convert") == "", "ImageMagick 'convert' not on PATH")
  set.seed(25)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, track_fit = TRUE, max_iter = 5, verbose = FALSE)

  track <- fit$trace
  expect_s3_class(track, "susie_track")

  temp_prefix <- tempfile("susie_trackobj_")
  on.exit({
    f <- list.files(dirname(temp_prefix), pattern = basename(temp_prefix), full.names = TRUE)
    if (length(f) > 0) file.remove(f)
  })

  expect_error({
    invisible(capture.output(
      suppressMessages(susie_plot_iteration(track, L = 5, file_prefix = temp_prefix)),
      type = "output"
    ))
  }, NA)
})

test_that("susie_plot_iteration removes pre-existing gif before re-creating it", {
  set.seed(26)
  dat <- simulate_regression(n = 80, p = 40, k = 2)
  fit <- suppressWarnings(
    susie(dat$X, dat$y, L = 3, track_fit = TRUE, max_iter = 4, verbose = FALSE)
  )

  temp_prefix <- tempfile("susie_gifremove_")
  on.exit({
    f <- list.files(dirname(temp_prefix), pattern = basename(temp_prefix), full.names = TRUE)
    if (length(f) > 0) file.remove(f)
  })

  # Pre-create the .gif so the file.remove branch is exercised
  writeLines("dummy", paste0(temp_prefix, ".gif"))
  expect_true(file.exists(paste0(temp_prefix, ".gif")))

  # convert may or may not succeed; if it fails the stop() is the correct path
  tryCatch(
    suppressMessages(susie_plot_iteration(fit, L = 3, file_prefix = temp_prefix)),
    error = function(e) {
      expect_match(conditionMessage(e), "convert command failed", fixed = TRUE)
    }
  )
})

test_that("susie_plot_iteration errors on invalid model input", {
  expect_error(
    susie_plot_iteration(list(a = 1), L = 5),
    "model must be a susie fit or a susie_track object"
  )
})

test_that("susie_plot_iteration errors when trace is not a susie_track object", {
  set.seed(27)
  dat <- simulate_regression(n = 50, p = 20, k = 2)
  fit <- susie(dat$X, dat$y, L = 3, track_fit = FALSE, verbose = FALSE)
  fit$trace <- list(not_a_track = TRUE)

  expect_error(
    susie_plot_iteration(fit, L = 3),
    "model\\$trace must be a susie_track object"
  )
})

# ---- susie_plot_changepoint ----

test_that("susie_plot_changepoint returns a ggplot object with various options", {
  set.seed(28)
  mu <- c(rep(0, 25), rep(2, 25), rep(-1, 25), rep(1, 25))
  y  <- mu + rnorm(100, sd = 0.5)
  s  <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))

  result_default <- susie_plot_changepoint(s, y)
  result_color   <- susie_plot_changepoint(s, y, line_col = "red", line_size = 2)
  result_cs_col  <- susie_plot_changepoint(s, y, cs_col = "green")
  result_combo   <- susie_plot_changepoint(s, y, line_col = "red", line_size = 2,
                                           cs_col = "blue")

  expect_s3_class(result_default, "ggplot")
  expect_s3_class(result_color,   "ggplot")
  expect_s3_class(result_cs_col,  "ggplot")
  expect_s3_class(result_combo,   "ggplot")

  # returned ggplot can be further modified
  expect_error(result_default + ggplot2::ggtitle("Custom Title"), NA)
})

test_that("susie_plot_changepoint finds credible sets for strong multi-changepoint data", {
  set.seed(29)
  mu <- c(rep(0, 25), rep(2, 25), rep(-1, 25), rep(1, 25))
  y  <- mu + rnorm(100, sd = 0.3)
  s  <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))

  result <- susie_plot_changepoint(s, y)
  expect_s3_class(result, "ggplot")

  cs <- susie_get_cs(s)
  expect_true(length(cs$cs) > 0)
})

test_that("susie_plot_changepoint handles no-changepoint data gracefully", {
  set.seed(30)
  y <- rnorm(50, mean = 5, sd = 0.5)
  s <- suppressWarnings(susie_trendfilter(y, order = 0, use_mad = FALSE))

  result <- susie_plot_changepoint(s, y)
  expect_s3_class(result, "ggplot")
})

# ---- all three plot functions exercised in sequence ----

test_that("all three plot functions work in sequence on the same data", {
  set.seed(31)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  # susie_plot
  fit1 <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  pdf(tempfile()); on.exit(dev.off(), add = TRUE)
  expect_error(susie_plot(fit1, "PIP"), NA)

  # susie_plot_changepoint
  mu   <- c(rep(0, 30), rep(2, 30))
  y_cp <- mu + rnorm(60, sd = 0.3)
  fit2 <- suppressWarnings(susie_trendfilter(y_cp, order = 0, use_mad = FALSE))
  expect_error(susie_plot_changepoint(fit2, y_cp), NA)

  # susie_plot_iteration (GIF only when ImageMagick is available)
  fit3 <- suppressWarnings(
    susie(dat$X, dat$y, L = 5, track_fit = TRUE, max_iter = 5, verbose = FALSE)
  )
  temp_prefix <- tempfile("test_seq_")
  on.exit({
    f <- list.files(dirname(temp_prefix), pattern = basename(temp_prefix), full.names = TRUE)
    if (length(f) > 0) file.remove(f)
  }, add = TRUE)

  if (Sys.which("convert") != "") {
    expect_error({
      invisible(capture.output(
        suppressMessages(susie_plot_iteration(fit3, L = 5, file_prefix = temp_prefix)),
        type = "output"
      ))
    }, NA)
  }
})

test_that("susie_plot legend labels a credible set with more than one variable", {
  # Two near-identical columns land both variables in one credible set, so the
  # legend takes the size>1 branch ("L#: C=<size>/R=<purity>") rather than C=1.
  set.seed(1)
  n <- 300; p <- 6
  X <- matrix(rnorm(n * p), n, p)
  X[, 2] <- X[, 1] + rnorm(n) * 0.04
  X <- scale(X)
  beta <- rep(0, p); beta[1] <- 2
  y <- as.vector(X %*% beta + rnorm(n))
  fit <- suppressWarnings(susie(X, y, L = 3, verbose = FALSE))
  expect_gt(max(vapply(fit$sets$cs, length, integer(1))), 1)  # a multi-variable CS

  pdf(tempfile()); on.exit(dev.off(), add = TRUE)
  expect_error(susie_plot(fit, "PIP", add_legend = TRUE), NA)
})
