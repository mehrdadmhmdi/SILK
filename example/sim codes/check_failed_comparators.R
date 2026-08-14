#!/usr/bin/env Rscript
# Audit the dedicated repair run. This script intentionally does not merge or
# overwrite results from the original confirmatory run.

source(file.path("R", "run_simulation_comparisons.R"))
source_silk_prediction_modules()

repair_methods <- c("JM-Correct", "JM-Misspecified", "DeepSurv")
plan <- build_task_plan()
expected <- expand.grid(
  task_id = plan$task_id,
  method = repair_methods,
  stringsAsFactors = FALSE
)
expected$key <- paste(expected$task_id, expected$method, sep = "::")

status_files <- list.files(
  RAW_DIR,
  pattern = "^prediction_status_task_[0-9]+[.]csv$",
  full.names = TRUE
)
if (!length(status_files)) {
  stop("No repair status files were found under ", RAW_DIR, call. = FALSE)
}

status <- do.call(rbind, lapply(status_files, function(path) {
  out <- utils::read.csv(path, stringsAsFactors = FALSE)
  out$source_file <- basename(path)
  out
}))
status <- status[status$method %in% repair_methods, , drop = FALSE]
status$key <- paste(status$task_id, status$method, sep = "::")

duplicate_keys <- unique(status$key[duplicated(status$key)])
missing_keys <- setdiff(expected$key, status$key)
unexpected_keys <- setdiff(status$key, expected$key)
as_flag <- function(x) tolower(as.character(x)) == "true"
status$success_flag <- as_flag(status$success)
status$fit_ok_flag <- as_flag(status$fit_ok)
status$predict_ok_flag <- as_flag(status$predict_ok)

summary <- do.call(rbind, lapply(repair_methods, function(method) {
  z <- status[status$method == method, , drop = FALSE]
  data.frame(
    method = method,
    expected_tasks = nrow(plan),
    reported_tasks = nrow(z),
    successful_tasks = sum(z$success_flag, na.rm = TRUE),
    fit_failures = sum(!z$fit_ok_flag, na.rm = TRUE),
    prediction_failures = sum(z$fit_ok_flag & !z$predict_ok_flag, na.rm = TRUE),
    missing_tasks = sum(!expected$key[expected$method == method] %in% z$key),
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(
  summary,
  file.path(SUMMARY_DIR, "failed_comparator_repair_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  status,
  file.path(SUMMARY_DIR, "failed_comparator_repair_status.csv"),
  row.names = FALSE
)

task_error_files <- list.files(
  RAW_DIR,
  pattern = "^prediction_failed_task_[0-9]+[.]csv$",
  full.names = TRUE
)
cat("Failed-comparator repair audit\n")
print(summary, row.names = FALSE)

problems <- character(0)
if (length(duplicate_keys)) {
  problems <- c(problems, paste(length(duplicate_keys), "duplicate method-task keys"))
}
if (length(missing_keys)) {
  problems <- c(problems, paste(length(missing_keys), "missing method-task keys"))
}
if (length(unexpected_keys)) {
  problems <- c(problems, paste(length(unexpected_keys), "unexpected method-task keys"))
}
if (length(task_error_files)) {
  problems <- c(problems, paste(length(task_error_files), "task-level error files"))
}
if (any(!status$success_flag)) {
  problems <- c(problems, paste(sum(!status$success_flag), "reported method failures"))
}
if (length(problems)) {
  stop("Repair audit failed: ", paste(problems, collapse = "; "), call. = FALSE)
}

cat("Repair audit passed: every JM-Correct, JM-Misspecified, and DeepSurv task succeeded.\n")
