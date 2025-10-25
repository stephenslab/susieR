context("Plotting functions")

# =============================================================================
# SUSIE_PLOT - BASIC FUNCTIONALITY
# =============================================================================

test_that("susie_plot with PIP creates plot without error", {
  set.seed(1)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Should not error
  expect_error(
    susie_plot(fit, "PIP"),
    NA
  )
})

test_that("susie_plot with z-scores requires compute_univariate_zscore", {
  set.seed(2)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, compute_univariate_zscore = FALSE, verbose = FALSE)

  # Should error when z-scores not computed
  expect_error(
    susie_plot(fit, "z"),
    "z-scores are not available"
  )
})

test_that("susie_plot with z_original also requires z-scores", {
  set.seed(51)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, compute_univariate_zscore = FALSE, verbose = FALSE)

  # Should error when trying to plot z_original without z-scores
  expect_error(
    susie_plot(fit, "z_original"),
    "z-scores are not available"
  )
})

test_that("susie_plot with z-scores works when available", {
  set.seed(3)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, compute_univariate_zscore = TRUE, verbose = FALSE)

  # Should not error
  expect_error(
    susie_plot(fit, "z"),
    NA
  )
})

test_that("susie_plot with z_original works", {
  set.seed(4)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, compute_univariate_zscore = TRUE, verbose = FALSE)

  expect_error(
    susie_plot(fit, "z_original"),
    NA
  )
})

test_that("susie_plot with log10PIP works", {
  set.seed(5)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_error(
    susie_plot(fit, "log10PIP"),
    NA
  )
})

test_that("susie_plot with invalid y type errors", {
  set.seed(6)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_error(
    susie_plot(fit, "invalid_type"),
    "Need to specify"
  )
})

test_that("susie_plot errors when pos list missing required elements", {
  set.seed(34)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$genomic_position <- 1000 + (1:length(fit$pip))

  # Missing 'attr'
  expect_error(
    susie_plot(fit, "PIP", pos = list(start = 1000, end = 1025)),
    "pos argument should be a list"
  )

  # Missing 'start'
  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "genomic_position", end = 1025)),
    "pos argument should be a list"
  )

  # Missing 'end'
  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "genomic_position", start = 1000)),
    "pos argument should be a list"
  )
})

test_that("susie_plot errors when pos$attr not in model", {
  set.seed(35)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "nonexistent_attr", start = 1, end = 25)),
    "Cannot find attribute nonexistent_attr"
  )
})

test_that("susie_plot errors when pos$start >= pos$end", {
  set.seed(36)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$genomic_position <- 1000 + (1:length(fit$pip))

  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "genomic_position", start = 1025, end = 1000)),
    "Position start should be smaller than end"
  )

  # Equal values
  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "genomic_position", start = 1000, end = 1000)),
    "Position start should be smaller than end"
  )
})

test_that("susie_plot errors when numeric pos outside range", {
  set.seed(37)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_error(
    susie_plot(fit, "PIP", pos = 1:100),  # Only 50 variables
    "Provided position is outside the range"
  )

  expect_error(
    susie_plot(fit, "PIP", pos = c(0, 1, 2)),  # 0 is out of range
    "Provided position is outside the range"
  )
})

# =============================================================================
# SUSIE_PLOT - PARAMETERS
# =============================================================================

test_that("susie_plot with add_bar=TRUE works", {
  set.seed(7)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_error(
    susie_plot(fit, "PIP", add_bar = TRUE),
    NA
  )
})

test_that("susie_plot with add_legend=TRUE works", {
  set.seed(8)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  expect_error(
    susie_plot(fit, "PIP", add_legend = TRUE),
    NA
  )
})

test_that("susie_plot with add_legend location string works", {
  set.seed(9)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  expect_error(
    susie_plot(fit, "PIP", add_legend = "bottomright"),
    NA
  )
})

test_that("susie_plot with pos as numeric vector works", {
  set.seed(10)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Plot subset of variables
  expect_error(
    susie_plot(fit, "PIP", pos = 1:25),
    NA
  )
})

test_that("susie_plot with pos as list works", {
  set.seed(11)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)
  fit$genomic_position <- 1000 + (1:length(fit$pip))

  # Plot with custom position attribute
  expect_error(
    susie_plot(fit, "PIP", pos = list(attr = "genomic_position", start = 1000, end = 1025)),
    NA
  )
})

test_that("susie_plot with b (true effects) works", {
  set.seed(12)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Highlight true effects
  expect_error(
    susie_plot(fit, "PIP", b = dat$beta),
    NA
  )
})

test_that("susie_plot with max_cs parameter works", {
  set.seed(13)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  # Limit number of CS displayed
  expect_error(
    susie_plot(fit, "PIP", max_cs = 2),
    NA
  )
})

test_that("susie_plot with max_cs purity threshold works", {
  set.seed(38)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  # Filter by purity (< 1)
  expect_error(
    susie_plot(fit, "PIP", max_cs = 0.5, add_legend = TRUE),
    NA
  )

  # Very strict purity filter (should exclude most/all CS)
  expect_error(
    susie_plot(fit, "PIP", max_cs = 0.99),
    NA
  )

  # Very lenient purity filter
  expect_error(
    susie_plot(fit, "PIP", max_cs = 0.1),
    NA
  )
})

test_that("susie_plot with different legend positions works", {
  set.seed(39)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  positions <- c("topleft", "top", "left", "center",
                 "right", "bottomleft", "bottom")

  for (pos in positions) {
    expect_error(
      susie_plot(fit, "PIP", add_legend = pos),
      NA,
      info = paste("Failed for legend position:", pos)
    )
  }
})

test_that("susie_plot with invalid legend position defaults to topright", {
  set.seed(40)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  # Invalid position should default to "topright" (no error)
  expect_error(
    susie_plot(fit, "PIP", add_legend = "invalid_position"),
    NA
  )
})

test_that("susie_plot respects custom plotting parameters", {
  set.seed(41)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Test various plotting parameters
  expect_error(
    susie_plot(fit, "PIP", main = "Custom Title", col = "blue", cex = 0.5),
    NA
  )

  expect_error(
    susie_plot(fit, "PIP", xlim = c(0, 30), ylim = c(0, 1)),
    NA
  )
})

# =============================================================================
# SUSIE_PLOT - VECTOR INPUT
# =============================================================================

test_that("susie_plot with PIP vector input works", {
  set.seed(14)
  pip <- runif(50)

  expect_error(
    susie_plot(pip, "PIP"),
    NA
  )
})

test_that("susie_plot with z-score vector input works", {
  set.seed(15)
  z <- rnorm(50)

  expect_error(
    susie_plot(z, "z"),
    NA
  )
})

test_that("susie_plot with non-susie vector and different y types", {
  set.seed(42)

  # Test with z_original on vector
  z_vec <- rnorm(50)
  expect_error(susie_plot(z_vec, "z_original"), NA)

  # Test with log10PIP on vector
  pip_vec <- runif(50)
  expect_error(susie_plot(pip_vec, "log10PIP"), NA)

  # Test with generic data (not PIP, z, etc.)
  data_vec <- runif(50, 0, 10)
  expect_error(susie_plot(data_vec, "custom_data"), NA)
})

# =============================================================================
# SUSIE_PLOT - EDGE CASES
# =============================================================================

test_that("susie_plot with no credible sets works", {
  set.seed(16)
  dat <- simulate_regression(n = 100, p = 50, k = 0)  # No signal
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_error(
    susie_plot(fit, "PIP"),
    NA
  )
})

test_that("susie_plot with single variable works", {
  set.seed(17)
  dat <- simulate_regression(n = 100, p = 1, k = 1)
  fit <- susie(dat$X, dat$y, L = 1, verbose = FALSE)

  expect_error(
    susie_plot(fit, "PIP"),
    NA
  )
})

test_that("susie_plot returns NULL invisibly", {
  set.seed(18)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  result <- susie_plot(fit, "PIP")

  expect_null(result)
})

test_that("susie_plot with list pos and credible sets adjusts correctly", {
  set.seed(43)
  # Create data with clear signal so we get CS
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)
  fit$genomic_position <- 1000 + (1:length(fit$pip))

  # Should successfully plot with CS adjusted to new positions
  expect_error(
    susie_plot(fit, "PIP",
               pos = list(attr = "genomic_position", start = 1000, end = 1050),
               add_legend = TRUE),
    NA
  )
})

test_that("susie_plot with b parameter highlights specific positions", {
  set.seed(44)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  # Test with b having non-zero elements at specific positions
  b_test <- rep(0, 50)
  b_test[c(10, 20, 30)] <- 1
  expect_error(
    susie_plot(fit, "PIP", b = b_test),
    NA
  )

  # Test with actual beta from simulation and add_bar
  expect_error(
    susie_plot(fit, "PIP", b = dat$beta, add_bar = TRUE, add_legend = TRUE),
    NA
  )
})

test_that("susie_plot sets x0 and y1 to NULL when CS filtered by max_cs", {
  set.seed(52)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  # Get CS with purity info
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
    # Use very strict max_cs filter (size < 1) to exclude CS
    # This should trigger the else branch: x0 <- NULL; y1 <- NULL
    expect_error(
      susie_plot(fit, "PIP", max_cs = 1, add_legend = TRUE),  # Only CS with size < 1
      NA
    )

    # Also test with very high purity threshold (max_cs as purity)
    expect_error(
      susie_plot(fit, "PIP", max_cs = 0.999, add_legend = TRUE),  # Very high purity
      NA
    )
  } else {
    skip("No CS found for max_cs filter test")
  }
})

test_that("susie_plot skips CS when x0 is NULL (next statement)", {
  set.seed(53)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  # Get CS
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
    # Use max_cs to filter out large CS, causing is.null(x0) to be TRUE
    # This should trigger the next statement to skip those CS
    expect_error(
      susie_plot(fit, "PIP", max_cs = 2),  # Skip CS with > 2 variables
      NA
    )
  } else {
    skip("No CS found for next statement test")
  }
})

test_that("susie_plot uses cs_index when available (else uses cs_idx)", {
  set.seed(54)
  dat <- simulate_regression(n = 200, p = 100, k = 3, signal_sd = 2)
  fit <- susie(dat$X, dat$y, L = 10, verbose = FALSE)

  # Get CS which should populate cs_index
  fit$sets <- susie_get_cs(fit, X = dat$X, coverage = 0.95)

  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
    # When cs_index exists, should use it
    expect_true(!is.null(fit$sets$cs_index))

    # Plot with legend to see cs_index values
    expect_error(
      susie_plot(fit, "PIP", add_legend = TRUE),
      NA
    )

    # Test the else branch: remove cs_index to force use of cs_idx
    fit_no_index <- fit
    fit_no_index$sets$cs_index <- NULL

    expect_error(
      susie_plot(fit_no_index, "PIP", add_legend = TRUE),
      NA
    )
  } else {
    skip("No CS found for cs_index test")
  }
})

# =============================================================================
# SUSIE_PLOT_ITERATION
# =============================================================================

test_that("susie_plot_iteration uses tempdir when file_prefix missing", {
  set.seed(55)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, track_fit = FALSE, verbose = FALSE)

  # Don't provide file_prefix - should use tempdir()
  result <- invisible(capture.output({
    suppressMessages(susie_plot_iteration(fit, L = 5))
  }, type = "output"))

  # Check that file was created in tempdir
  expected_path <- file.path(tempdir(), "susie_plot.pdf")
  expect_true(file.exists(expected_path))

  # Clean up
  if (file.exists(expected_path)) file.remove(expected_path)
})

test_that("susie_plot_iteration with track_fit=FALSE uses final fit only", {
  set.seed(19)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  # Without track_fit
  fit <- susie(dat$X, dat$y, L = 5, track_fit = FALSE, verbose = FALSE)

  temp_prefix <- tempfile("susie_iter_no_track_")

  # Should work but only plot final iteration
  expect_error({
    invisible(capture.output({
      suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix))
    }, type = "output"))
  }, NA)

  # Clean up
  temp_files <- list.files(dirname(temp_prefix),
                           pattern = basename(temp_prefix),
                           full.names = TRUE)
  file.remove(temp_files)
})

test_that("susie_plot_iteration works with track_fit=TRUE", {
  set.seed(20)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  # With track_fit
  fit <- susie(dat$X, dat$y, L = 5, track_fit = TRUE, max_iter = 10, verbose = FALSE)

  # Create temp file for output
  temp_prefix <- tempfile("susie_iter_")

  expect_error({
    invisible(capture.output({
      suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix))
    }, type = "output"))
  }, NA)

  # Clean up temp files
  temp_files <- list.files(dirname(temp_prefix),
                           pattern = basename(temp_prefix),
                           full.names = TRUE)
  file.remove(temp_files)
})

test_that("susie_plot_iteration with pos parameter works", {
  set.seed(21)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, track_fit = TRUE, max_iter = 10, verbose = FALSE)

  temp_prefix <- tempfile("susie_iter_pos_")

  expect_error({
    invisible(capture.output({
      suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix, pos = 1:25))
    }, type = "output"))
  }, NA)

  # Clean up
  temp_files <- list.files(dirname(temp_prefix),
                           pattern = basename(temp_prefix),
                           full.names = TRUE)
  file.remove(temp_files)
})

test_that("susie_plot_iteration creates PDF files", {
  set.seed(22)
  dat <- simulate_regression(n = 100, p = 50, k = 3)

  fit <- susie(dat$X, dat$y, L = 5, track_fit = TRUE, max_iter = 5, verbose = FALSE)

  temp_dir <- tempdir()
  temp_prefix <- file.path(temp_dir, "test_susie_iter")

  invisible(capture.output({
    suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix))
  }, type = "output"))

  # Check that PDF files were created
  pdf_files <- list.files(temp_dir, pattern = "test_susie_iter.*\\.pdf$", full.names = TRUE)

  expect_true(length(pdf_files) > 0)

  # Clean up
  file.remove(pdf_files)
})

test_that("susie_plot_iteration with L greater than nrow(alpha)", {
  set.seed(45)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 3, track_fit = TRUE, verbose = FALSE)

  temp_prefix <- tempfile("test_large_L_")

  # Request L=10 when fit only has L=3
  expect_error({
    invisible(capture.output({
      suppressMessages(susie_plot_iteration(fit, L = 10, file_prefix = temp_prefix))
    }, type = "output"))
  }, NA)

  # Clean up
  temp_files <- list.files(dirname(temp_prefix),
                          pattern = basename(temp_prefix),
                          full.names = TRUE)
  if (length(temp_files) > 0) file.remove(temp_files)
})

test_that("susie_plot_iteration with few iterations", {
  set.seed(46)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  # Use max_iter=5 to ensure we have some iterations tracked
  fit <- susie(dat$X, dat$y, L = 5, track_fit = TRUE, max_iter = 5, verbose = FALSE)

  temp_prefix <- tempfile("test_few_iter_")

  expect_error({
    invisible(capture.output({
      suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix))
    }, type = "output"))
  }, NA)

  # Clean up
  temp_files <- list.files(dirname(temp_prefix),
                          pattern = basename(temp_prefix),
                          full.names = TRUE)
  if (length(temp_files) > 0) file.remove(temp_files)
})

test_that("susie_plot_iteration returns invisibly", {
  set.seed(47)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit <- susie(dat$X, dat$y, L = 5, track_fit = FALSE, verbose = FALSE)

  temp_prefix <- tempfile("test_invisible_")

  invisible(capture.output({
    result <- suppressMessages(susie_plot_iteration(fit, L = 5, file_prefix = temp_prefix))
  }, type = "output"))

  expect_null(result)

  # Clean up
  temp_files <- list.files(dirname(temp_prefix),
                          pattern = basename(temp_prefix),
                          full.names = TRUE)
  if (length(temp_files) > 0) file.remove(temp_files)
})

# =============================================================================
# SUSIE_PLOT_CHANGEPOINT
# =============================================================================

test_that("susie_plot_changepoint with basic usage works", {
  set.seed(23)
  mu <- c(rep(0, 25), rep(2, 25), rep(-1, 25), rep(1, 25))
  y <- mu + rnorm(100, sd = 0.5)

  s <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  expect_error(
    susie_plot_changepoint(s, y),
    NA
  )
})

test_that("susie_plot_changepoint with custom colors works", {
  set.seed(24)
  mu <- c(rep(0, 30), rep(2, 30))
  y <- mu + rnorm(60, sd = 0.3)

  s <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  expect_error(
    susie_plot_changepoint(s, y, line_col = "red", line_size = 2),
    NA
  )
})

test_that("susie_plot_changepoint with cs_col parameter works", {
  set.seed(25)
  mu <- c(rep(0, 30), rep(2, 30))
  y <- mu + rnorm(60, sd = 0.3)

  s <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  expect_error(
    susie_plot_changepoint(s, y, cs_col = "green"),
    NA
  )
})

test_that("susie_plot_changepoint with single changepoint works", {
  set.seed(28)
  mu <- c(rep(0, 30), rep(2, 30))
  y <- mu + rnorm(60, sd = 0.3)

  s <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  expect_error(
    susie_plot_changepoint(s, y),
    NA
  )
})

test_that("susie_plot_changepoint with no changepoints works", {
  set.seed(29)
  y <- rnorm(50, mean = 5, sd = 0.5)

  s <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  expect_error(
    susie_plot_changepoint(s, y),
    NA
  )
})

test_that("susie_plot_changepoint returns ggplot object", {
  set.seed(30)
  mu <- c(rep(0, 30), rep(2, 30))
  y <- mu + rnorm(60, sd = 0.3)

  s <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  result <- susie_plot_changepoint(s, y)

  expect_s3_class(result, "gg")
  expect_s3_class(result, "ggplot")
})

test_that("susie_plot_changepoint with multiple changepoints", {
  set.seed(48)
  # Create data with multiple clear changepoints
  mu <- c(rep(0, 25), rep(2, 25), rep(-1, 25), rep(1, 25))
  y <- mu + rnorm(100, sd = 0.3)

  s <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  result <- susie_plot_changepoint(s, y)

  # Verify it's a ggplot object
  expect_s3_class(result, "gg")

  # Check that CS were found
  cs <- susie_get_cs(s)
  expect_true(length(cs$cs) > 0)
})

test_that("susie_plot_changepoint with very strong signal", {
  set.seed(49)
  # Very clear changepoints with low noise
  mu <- c(rep(0, 30), rep(5, 30))
  y <- mu + rnorm(60, sd = 0.1)

  s <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  result <- susie_plot_changepoint(s, y, line_col = "red", line_size = 2, cs_col = "blue")

  expect_s3_class(result, "gg")
  expect_s3_class(result, "ggplot")
})

test_that("susie_plot_changepoint can be modified after creation", {
  set.seed(50)
  mu <- c(rep(0, 30), rep(2, 30))
  y <- mu + rnorm(60, sd = 0.3)

  s <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  # Use all defaults
  result <- susie_plot_changepoint(s, y)

  expect_s3_class(result, "ggplot")

  # Can add to the plot after creation
  expect_error(
    result + ggplot2::ggtitle("Custom Title"),
    NA
  )
})

# =============================================================================
# INTEGRATION TESTS
# =============================================================================

test_that("susie_plot works with susie_ss output", {
  set.seed(31)
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  ss <- compute_summary_stats(dat$X, dat$y)

  fit <- susie_ss(ss$XtX, ss$Xty, ss$yty, n = ss$n, L = 5, verbose = FALSE)

  expect_error(
    susie_plot(fit, "PIP"),
    NA
  )
})

test_that("susie_plot works with susie_rss output", {
  set.seed(32)
  n <- 200
  p <- 100
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[1:3] <- 1
  y <- X %*% beta + rnorm(n)

  ss <- univariate_regression(X, y)
  R <- cor(X)
  z <- with(ss, betahat / sebetahat)

  fit <- susie_rss(z, R, n = n, L = 5, verbose = FALSE)

  expect_error(
    susie_plot(fit, "PIP"),
    NA
  )
})

test_that("all three plot functions work in sequence", {
  set.seed(33)

  # Regular susie fit
  dat <- simulate_regression(n = 100, p = 50, k = 3)
  fit1 <- susie(dat$X, dat$y, L = 5, verbose = FALSE)

  expect_error(susie_plot(fit1, "PIP"), NA)

  # Trendfilter
  mu <- c(rep(0, 30), rep(2, 30))
  y <- mu + rnorm(60, sd = 0.3)
  fit2 <- susie_trendfilter(y, order = 0, use_mad = FALSE)

  expect_error(susie_plot_changepoint(fit2, y), NA)

  # Iteration plot
  fit3 <- susie(dat$X, dat$y, L = 5, track_fit = TRUE, max_iter = 5, verbose = FALSE)
  temp_prefix <- tempfile("test_seq_")

  expect_error({
    invisible(capture.output({
      suppressMessages(susie_plot_iteration(fit3, L = 5, file_prefix = temp_prefix))
    }, type = "output"))
  }, NA)

  # Clean up
  temp_files <- list.files(dirname(temp_prefix),
                           pattern = basename(temp_prefix),
                           full.names = TRUE)
  file.remove(temp_files)
})
