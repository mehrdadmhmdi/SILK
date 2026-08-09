#!/usr/bin/env Rscript
# =============================================================================
# task.R
# One Slurm-array task = one Monte Carlo replication for one design cell.
# Outputs are prediction-focused subject-level risks plus derived metrics.
# =============================================================================

source(file.path("R", "run_simulation_comparisons.R"))
source_silk_prediction_modules()

cell <- select_task_from_env()
pred_file <- file.path(RAW_DIR, sprintf("prediction_task_%04d.csv", cell$task_id))
metric_file <- file.path(RAW_DIR, sprintf("prediction_metrics_task_%04d.csv", cell$task_id))
per_horizon_file <- file.path(RAW_DIR, sprintf("prediction_per_horizon_task_%04d.csv", cell$task_id))
paired_file <- file.path(RAW_DIR, sprintf("prediction_paired_task_%04d.csv", cell$task_id))
calibration_file <- file.path(RAW_DIR, sprintf("prediction_calibration_task_%04d.csv", cell$task_id))
status_file <- file.path(RAW_DIR, sprintf("prediction_status_task_%04d.csv", cell$task_id))
registration_file <- file.path(RAW_DIR, sprintf("registration_diagnostics_task_%04d.csv", cell$task_id))
error_file <- file.path(RAW_DIR, sprintf("prediction_failed_task_%04d.csv", cell$task_id))

res <- tryCatch(
  run_simulation_comparisons(cell),
  error = function(e) {
    data.frame(
      success = FALSE,
      task_id = cell$task_id,
      cell_id = cell$cell_id,
      phase = cell$phase,
      scenario = cell$scenario,
      schedule = cell$schedule,
      n_visits = cell$n_visits,
      n_train = cell$n_train,
      n_test = cell$n_test,
      rep = cell$rep,
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  }
)

if (is.data.frame(res) && identical(res$success[1], FALSE)) {
  write.csv(res, error_file, row.names = FALSE)
  cat("Task failed; wrote", error_file, "\n")
} else {
  write.csv(res$predictions, pred_file, row.names = FALSE)
  write.csv(res$replicate_metrics, metric_file, row.names = FALSE)
  write.csv(res$per_horizon, per_horizon_file, row.names = FALSE)
  write.csv(res$paired, paired_file, row.names = FALSE)
  write.csv(res$calibration_bins, calibration_file, row.names = FALSE)
  write.csv(res$method_status, status_file, row.names = FALSE)
  if (nrow(res$registration_diagnostics)) {
    write.csv(res$registration_diagnostics, registration_file, row.names = FALSE)
  }
  cat("Wrote", pred_file, "\n")
  cat("Wrote", metric_file, "\n")
  cat("Wrote", status_file, "\n")
  if (nrow(res$registration_diagnostics)) cat("Wrote", registration_file, "\n")
}
