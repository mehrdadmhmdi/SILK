# =============================================================================
# pilot_two_visit_resolution.R
#
# Scientific resolution diagnostic for sparse visit schedules. This is a PILOT
# script, not a package unit test: it takes minutes per cell, and the quantity
# it reports is a research finding whose acceptable value is not a code
# contract.
#
# What it answers. In the 8.9.2026 freeze the two-visit schedule produced
# shift correlation approximately 0.08 and estimated shift SD approximately 1.1
# against a true SD of 5, because the longitudinal-signal gate zeroed the clock
# initializer. After the gate fix the question is how much resolution a
# two-visit design actually supports.
#
# Estimand caveat. Shift and stage recovery are NOT the primary estimand; SILK
# targets the calibrated-risk functional. These statistics are latent-resolution
# diagnostics and algorithmic debugging checks only. Judge the method on
# held-out IBS and calibration, via report_two_questions.R.
#
# Sign convention. The correlation between the estimated and true shift must be
# POSITIVE. A large negative correlation is a failure, not a success, so this
# script asserts the signed value and never abs().
#
# Usage:
#   Rscript pilot_two_visit_resolution.R [n_train] [n_rep] [output_csv]
# =============================================================================

suppressPackageStartupMessages(library(SILK))
source(file.path("R", "cfg.R"))

args <- commandArgs(trailingOnly = TRUE)
N_TRAIN <- if (length(args) >= 1L) as.integer(args[[1]]) else 400L
N_REP <- if (length(args) >= 2L) as.integer(args[[2]]) else 10L
OUT <- if (length(args) >= 3L) args[[3]] else "pilot_two_visit_resolution.csv"

SCHEDULES <- c("m2", "m4", "m12")
SCENARIOS <- c("mean_severe", "mean_strong_dense", "dist_large", "weak_stage")

rows <- list()
for (scn in SCENARIOS) {
  grid <- make_shift_grid(scn)
  for (sch in SCHEDULES) {
    for (rep_id in seq_len(N_REP)) {
      seed <- 20260810L + 1000L * match(scn, SCENARIOS) + 100L * match(sch, SCHEDULES) + rep_id
      dat <- generate_dataset_fixed(N_TRAIN, scn, schedule_name = sch, seed = seed)
      reg <- try(
        fit_silk_registration(dat$subjects, dat$visits, shift_grid = grid, seed = seed),
        silent = TRUE
      )
      if (inherits(reg, "try-error")) next

      e_hat <- reg$train_stage$e_hat[match(dat$subjects$id, reg$train_stage$id)]
      e_true <- dat$subjects$eps
      ok <- is.finite(e_hat) & is.finite(e_true)
      signal <- reg$final_template$longitudinal_signal
      ms <- reg$multistart

      rows[[length(rows) + 1L]] <- data.frame(
        scenario = scn,
        schedule = sch,
        rep = rep_id,
        sigma_eps = get_scenario(scn)$sigma_eps,
        gate_statistic = if (is.null(signal)) NA_real_ else signal,
        gate_evaluable = !is.null(signal) && is.finite(signal),
        shift_correlation = if (sum(ok) > 2L && stats::sd(e_true[ok]) > 1e-12) {
          stats::cor(e_hat[ok], e_true[ok])
        } else NA_real_,
        shift_rmse = sqrt(mean((e_hat[ok] - e_true[ok])^2)),
        estimated_shift_sd = stats::sd(e_hat[ok]),
        true_shift_sd = stats::sd(e_true[ok]),
        boundary_fraction = mean(as.logical(reg$train_stage$at_boundary), na.rm = TRUE),
        median_gap_q1 = stats::median(reg$train_stage$gap_q1, na.rm = TRUE),
        n_successful_starts = if (is.null(ms)) NA_real_ else ms$n_successful_starts,
        objective_spread = if (is.null(ms)) NA_real_ else ms$objective_spread,
        second_best_gap = if (is.null(ms)) NA_real_ else ms$second_best_gap,
        stringsAsFactors = FALSE
      )
    }
  }
}

out <- do.call(rbind, rows)
utils::write.csv(out, OUT, row.names = FALSE)

agg <- stats::aggregate(
  cbind(shift_correlation, shift_rmse, estimated_shift_sd, boundary_fraction,
        objective_spread) ~ scenario + schedule,
  data = out, FUN = function(z) mean(z, na.rm = TRUE)
)
print(agg[order(agg$scenario, agg$schedule), ], row.names = FALSE)

# Signed check. weak_stage is the intended null, so it is exempt.
informative <- out[out$scenario != "weak_stage" & is.finite(out$shift_correlation), ]
if (nrow(informative)) {
  worst <- min(informative$shift_correlation)
  message(sprintf("\nMinimum signed shift correlation over informative cells: %.3f", worst))
  if (worst <= 0) {
    message("FAIL: registration is not positively associated with the true shift ",
            "in at least one informative cell. Investigate before the confirmatory run.")
  }
}
message("Wrote: ", normalizePath(OUT, mustWork = FALSE))
