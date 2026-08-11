# Structural unit tests for the longitudinal-signal gate.
#
# Motivation: in the 8.9.2026 freeze the two-visit schedule zeroed the biomarker
# clock initializer. Every subject contributed exactly one within-subject pair at
# the same lag gap, so the gate statistic had zero lag-distance variance,
# returned 0, and fell below CLOCK_SIGNAL_MIN. The statistic was undefined, not
# zero.
#
# These tests are deliberately FAST and structural. They assert what the gate
# does, not how well registration recovers a latent quantity. The scientific
# resolution check (shift correlation and RMSE at two visits) is a pilot
# diagnostic, not a package unit test -- it costs minutes per scenario and its
# threshold is a research finding rather than a code contract. It lives in
# example/sim codes/pilot_two_visit_resolution.R.

make_visits <- function(n, lags, p = 2L, seed = 11L) {
  set.seed(seed)
  v <- data.frame(
    id = rep(seq_len(n), each = length(lags)),
    visit = rep(seq_along(lags), times = n),
    lag = rep(lags, times = n),
    X1 = 0, X2 = 0
  )
  for (k in seq_len(p)) v[[paste0("B", k)]] <- stats::rnorm(nrow(v))
  v
}

gate_statistic <- function(visits) {
  builder <- make_biomarker_kernel_builder(visits, "gaussian")
  longitudinal_kernel_signal(visits, apply_biomarker_builder(visits, builder), builder)
}

test_that("a balanced two-visit design makes the gate statistic non-evaluable", {
  expect_true(is.na(gate_statistic(make_visits(30L, c(6, 0)))))
})

test_that("a three-visit design yields an evaluable gate statistic", {
  expect_false(is.na(gate_statistic(make_visits(30L, c(6, 3, 0)))))
})

test_that("an unbalanced two-visit design is evaluable", {
  # Two visits per subject but varying lag gaps across subjects restores
  # lag-distance variance, so the statistic is defined again.
  set.seed(5)
  n <- 30L
  first_lag <- stats::runif(n, 2, 8)
  v <- data.frame(
    id = rep(seq_len(n), each = 2L),
    visit = rep(1:2, times = n),
    lag = as.vector(rbind(first_lag, 0)),
    X1 = 0, X2 = 0
  )
  v$B1 <- stats::rnorm(nrow(v))
  v$B2 <- stats::rnorm(nrow(v))
  builder <- make_biomarker_kernel_builder(v, "gaussian")
  expect_false(is.na(
    longitudinal_kernel_signal(v, apply_biomarker_builder(v, builder), builder)
  ))
})

test_that("a non-evaluable gate statistic leaves the clock enabled", {
  expect_false(clock_gate_disables(NA_real_))
  expect_false(clock_gate_disables(NA))
  expect_false(clock_gate_disables(NULL))
  expect_false(clock_gate_disables(numeric(0)))
})

test_that("an evaluable statistic below the minimum still disables the clock", {
  old <- getOption("silk.CLOCK_SIGNAL_MIN")
  on.exit(options(silk.CLOCK_SIGNAL_MIN = old), add = TRUE)
  silk_options(CLOCK_SIGNAL_MIN = 0.10)

  expect_true(clock_gate_disables(0))
  expect_true(clock_gate_disables(0.05))
  expect_false(clock_gate_disables(0.10))
  expect_false(clock_gate_disables(0.90))
})

test_that("the held-out template path also respects the non-evaluable gate", {
  # template_clock_initial_shift() returns all-zero shifts when the gate fires.
  # With a non-evaluable statistic it must NOT short-circuit, so the returned
  # shifts are not identically zero for a design with real stage signal.
  skip_if_not(exists("generate_dataset_fixed"))
  dat <- generate_dataset_fixed(80, "mean_severe", schedule_name = "m2", seed = 3)
  reg <- fit_silk_registration(dat$subjects, dat$visits, seed = 3)

  expect_true(is.na(reg$final_template$longitudinal_signal))
  shift <- template_clock_initial_shift(reg$final_template, dat$visits)
  expect_false(all(abs(shift) < 1e-8))
})

test_that("multistart records an objective for every start", {
  skip_if_not(exists("generate_dataset_fixed"))
  old <- getOption("silk.N_STARTS")
  on.exit(options(silk.N_STARTS = old), add = TRUE)
  silk_options(N_STARTS = 3L)

  dat <- generate_dataset_fixed(80, "mean_moderate", seed = 4)
  reg <- fit_silk_registration(dat$subjects, dat$visits, seed = 4)
  ms <- reg$multistart

  expect_equal(ms$n_starts, 3L)
  expect_gte(ms$n_successful_starts, 1L)
  expect_lte(ms$n_successful_starts, ms$n_starts)
  expect_true(is.finite(ms$best_objective))
  expect_true(ms$best_start %in% strsplit(ms$start_names, "|", fixed = TRUE)[[1]])
  if (ms$n_successful_starts >= 2L) {
    expect_true(is.finite(ms$objective_spread))
    expect_gte(ms$objective_spread, 0)
    expect_gte(ms$second_best_gap, 0)
  }
})
