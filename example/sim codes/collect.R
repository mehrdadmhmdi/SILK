#!/usr/bin/env Rscript
# =============================================================================
# collect.R
# Collect prediction-focused simulation outputs and create summaries.
# =============================================================================

source(file.path("R", "run_simulation_comparisons.R"))
source_silk_prediction_modules()

task_id_from_path <- function(path) {
  as.integer(sub(".*_task_([0-9]+)\\.csv$", "\\1", basename(path)))
}

# Optional extra raw directories (colon-separated) holding add-on method runs,
# e.g. SILK_EXTRA_RAW_DIRS="outputs_addmethods/raw". Their rows are merged into
# the data; the completeness audit below stays on the main RAW_DIR only.
EXTRA_RAW_DIRS <- Filter(nzchar, strsplit(Sys.getenv("SILK_EXTRA_RAW_DIRS", unset = ""), ":")[[1]])
RAW_DIRS_ALL <- c(RAW_DIR, EXTRA_RAW_DIRS)
list_raw <- function(pattern) {
  unlist(lapply(RAW_DIRS_ALL, function(d) list.files(d, pattern = pattern, full.names = TRUE)), use.names = FALSE)
}

metric_files <- list_raw("^prediction_metrics_task_[0-9]+\\.csv$")
per_horizon_files <- list_raw("^prediction_per_horizon_task_[0-9]+\\.csv$")
prediction_files <- list_raw("^prediction_task_[0-9]+\\.csv$")
calibration_files <- list_raw("^prediction_calibration_task_[0-9]+\\.csv$")
status_files <- list_raw("^prediction_status_task_[0-9]+\\.csv$")
failed_files <- list_raw("^prediction_failed_task_[0-9]+\\.csv$")
# Audit (missing-task detection) uses the MAIN dir only: a cell is "present"
# when its full-roster run exists in RAW_DIR; add-on dirs only add methods.
main_metric_files <- list.files(RAW_DIR, pattern = "^prediction_metrics_task_[0-9]+\\.csv$", full.names = TRUE)

if (length(metric_files) == 0L) {
  stop("No prediction_metrics_task_*.csv files found in ", RAW_DIR, call. = FALSE)
}

plan <- build_task_plan()
success_tasks <- sort(unique(task_id_from_path(main_metric_files)))
missing <- setdiff(plan$task_id, success_tasks)
failed <- if (length(failed_files)) do.call(rbind, lapply(failed_files, read.csv, stringsAsFactors = FALSE)) else data.frame()

audit <- data.frame(
  expected_tasks = nrow(plan),
  successful_tasks = length(success_tasks),
  missing_tasks = length(missing),
  failed_rows = nrow(failed),
  complete = length(missing) == 0L && nrow(failed) == 0L,
  stringsAsFactors = FALSE
)
write.csv(audit, file.path(SUMMARY_DIR, "prediction_replication_audit.csv"), row.names = FALSE)
if (length(missing) > 0L) write.csv(plan[plan$task_id %in% missing, ], file.path(SUMMARY_DIR, "prediction_missing_tasks.csv"), row.names = FALSE)
if (nrow(failed) > 0L) write.csv(failed, file.path(SUMMARY_DIR, "prediction_failed_rows.csv"), row.names = FALSE)

allow_incomplete <- identical(tolower(Sys.getenv("SILK_ALLOW_INCOMPLETE", unset = "false")), "true")
if (!audit$complete && !allow_incomplete) {
  stop("Collection incomplete. Set SILK_ALLOW_INCOMPLETE=true for provisional summaries.", call. = FALSE)
}

per_horizon <- do.call(rbind, lapply(per_horizon_files, read.csv, stringsAsFactors = FALSE))
replicate_metrics <- do.call(rbind, lapply(metric_files, read.csv, stringsAsFactors = FALSE))
summary <- prediction_method_summary(replicate_metrics)
paired <- prediction_paired_differences(replicate_metrics)
calibration_bins <- if (length(calibration_files)) {
  do.call(rbind, lapply(calibration_files, read.csv, stringsAsFactors = FALSE))
} else {
  data.frame()
}
calibration_summary <- summarize_calibration_bins(calibration_bins)
method_status <- if (length(status_files)) {
  do.call(rbind, lapply(status_files, read.csv, stringsAsFactors = FALSE))
} else {
  data.frame()
}
method_failure_summary <- if (nrow(method_status)) {
  status_true <- function(x) x %in% c(TRUE, "TRUE", "true", "True", 1, "1")
  do.call(rbind, lapply(split(method_status, method_status$method), function(z) {
    has_fit_warning <- !is.na(z$fit_warnings) & nzchar(z$fit_warnings)
    has_predict_warning <- !is.na(z$predict_warnings) & nzchar(z$predict_warnings)
    fit_ok <- status_true(z$fit_ok)
    predict_ok <- status_true(z$predict_ok)
    success <- status_true(z$success)
    data.frame(
      method = z$method[1],
      attempted_tasks = nrow(z),
      successful_tasks = sum(success),
      failed_tasks = sum(!success),
      fit_failures = sum(!fit_ok),
      prediction_failures = sum(fit_ok & !predict_ok),
      warning_rows = sum(has_fit_warning | has_predict_warning),
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame()
}

write.csv(plan, file.path(SUMMARY_DIR, "prediction_task_plan.csv"), row.names = FALSE)
write.csv(per_horizon, file.path(SUMMARY_DIR, "prediction_per_horizon_metrics.csv"), row.names = FALSE)
write.csv(replicate_metrics, file.path(SUMMARY_DIR, "prediction_replicate_metrics.csv"), row.names = FALSE)
write.csv(summary, file.path(SUMMARY_DIR, "prediction_method_summary.csv"), row.names = FALSE)
write.csv(paired, file.path(SUMMARY_DIR, "prediction_paired_differences.csv"), row.names = FALSE)
if (nrow(calibration_bins)) write.csv(calibration_bins, file.path(SUMMARY_DIR, "prediction_calibration_bins_by_task.csv"), row.names = FALSE)
if (nrow(calibration_summary)) write.csv(calibration_summary, file.path(SUMMARY_DIR, "prediction_calibration_bins.csv"), row.names = FALSE)
if (nrow(method_status)) write.csv(method_status, file.path(SUMMARY_DIR, "prediction_method_status.csv"), row.names = FALSE)
if (nrow(method_failure_summary)) write.csv(method_failure_summary, file.path(SUMMARY_DIR, "prediction_method_failure_summary.csv"), row.names = FALSE)

# Compatibility aliases for old job wrappers or notebooks that expect fixed_* names.
write.csv(replicate_metrics, file.path(SUMMARY_DIR, "fixed_all_replications.csv"), row.names = FALSE)
write.csv(summary, file.path(SUMMARY_DIR, "fixed_method_summary.csv"), row.names = FALSE)
write.csv(paired, file.path(SUMMARY_DIR, "fixed_paired_differences.csv"), row.names = FALSE)

collect_predictions <- identical(tolower(Sys.getenv("SILK_COLLECT_PREDICTIONS", unset = "false")), "true")
if (collect_predictions && length(prediction_files)) {
  predictions <- do.call(rbind, lapply(prediction_files, function(f) validate_prediction_frame(read.csv(f, stringsAsFactors = FALSE))))
  write.csv(predictions, file.path(SUMMARY_DIR, "prediction_all_rows.csv"), row.names = FALSE)
}

cat("Collected prediction simulation results\n")
cat("Successful tasks:", length(success_tasks), "of", nrow(plan), "\n")
cat("Summary:", file.path(SUMMARY_DIR, "prediction_method_summary.csv"), "\n")
