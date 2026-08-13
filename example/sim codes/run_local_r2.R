#!/usr/bin/env Rscript
# =============================================================================
# run_local_r2.R
# Serial local driver for the R = 2 smoke run on the retuned DGM.
#
# Scope (2026-08-13): n_train = 400 only, visit schedules m4 and m12 only.
# Writes the same per-task CSVs as task.R into RAW_DIR, so collect.R and
# report_two_questions.R consume the output unchanged. Tasks already on disk
# are skipped, so the driver can be interrupted and restarted.
# =============================================================================

lib <- Sys.getenv("SILK_EXTRA_LIB", unset = "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))

source(file.path("R", "run_simulation_comparisons.R"))
source_silk_prediction_modules()

plan <- build_task_plan()
cat("SILK loaded from :", find.package("SILK"), "\n")
cat("tasks            :", nrow(plan), "\n")
cat("n_train grid     :", paste(unique(plan$n_train), collapse = ","), "\n")
cat("schedules        :", paste(sort(unique(plan$schedule)), collapse = ","), "\n")
cat("methods          :", length(METHOD_ORDER), "\n")
cat("raw dir          :", RAW_DIR, "\n\n")
flush.console()

only <- Sys.getenv("SILK_ONLY_TASKS", unset = "")
task_ids <- if (nzchar(only)) as.integer(trimws(strsplit(only, ",")[[1]])) else plan$task_id

started <- proc.time()[3]
done <- 0L
for (tid in task_ids) {
  cell <- plan[plan$task_id == tid, , drop = FALSE]
  metric_file <- file.path(RAW_DIR, sprintf("prediction_metrics_task_%04d.csv", tid))
  if (file.exists(metric_file)) {
    cat(sprintf("[%3d/%3d] task %4d  SKIP (already on disk)\n",
                done + 1L, length(task_ids), tid))
    done <- done + 1L
    flush.console()
    next
  }
  t0 <- proc.time()[3]
  res <- tryCatch(run_simulation_comparisons(cell), error = function(e) e)
  elapsed <- proc.time()[3] - t0

  if (inherits(res, "error")) {
    write.csv(
      data.frame(success = FALSE, task_id = tid, cell_id = cell$cell_id,
                 phase = cell$phase, scenario = cell$scenario,
                 schedule = cell$schedule, n_train = cell$n_train,
                 rep = cell$rep, error_message = conditionMessage(res),
                 stringsAsFactors = FALSE),
      file.path(RAW_DIR, sprintf("prediction_failed_task_%04d.csv", tid)),
      row.names = FALSE
    )
    cat(sprintf("[%3d/%3d] task %4d  %-20s %-4s rep%d  FAILED (%.0fs): %s\n",
                done + 1L, length(task_ids), tid, cell$scenario, cell$schedule,
                cell$rep, elapsed, conditionMessage(res)))
  } else {
    write.csv(res$predictions, file.path(RAW_DIR, sprintf("prediction_task_%04d.csv", tid)), row.names = FALSE)
    write.csv(res$replicate_metrics, metric_file, row.names = FALSE)
    write.csv(res$per_horizon, file.path(RAW_DIR, sprintf("prediction_per_horizon_task_%04d.csv", tid)), row.names = FALSE)
    write.csv(res$paired, file.path(RAW_DIR, sprintf("prediction_paired_task_%04d.csv", tid)), row.names = FALSE)
    write.csv(res$calibration_bins, file.path(RAW_DIR, sprintf("prediction_calibration_task_%04d.csv", tid)), row.names = FALSE)
    write.csv(res$method_status, file.path(RAW_DIR, sprintf("prediction_status_task_%04d.csv", tid)), row.names = FALSE)
    if (nrow(res$registration_diagnostics)) {
      write.csv(res$registration_diagnostics,
                file.path(RAW_DIR, sprintf("registration_diagnostics_task_%04d.csv", tid)),
                row.names = FALSE)
    }
    nfail <- sum(!res$method_status$success)
    cat(sprintf("[%3d/%3d] task %4d  %-20s %-4s rep%d  ok (%.0fs, %d method failures)\n",
                done + 1L, length(task_ids), tid, cell$scenario, cell$schedule,
                cell$rep, elapsed, nfail))
  }
  done <- done + 1L
  flush.console()
}

cat(sprintf("\nfinished %d tasks in %.1f min\n", done, (proc.time()[3] - started) / 60))
