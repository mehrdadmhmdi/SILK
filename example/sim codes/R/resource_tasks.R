# =============================================================================
# resource_tasks.R
# Split the Monte Carlo task grid into resource-aware Slurm arrays.
# =============================================================================

if (!exists("build_task_plan", mode = "function")) {
  source(file.path("R", "cfg.R"))
}

RESOURCE_CLASSES <- data.frame(
  resource_class = c("standard", "medium", "high_memory", "long_time", "extreme"),
  sbatch_file = c(
    "job_array_standard.sbatch",
    "job_array_medium.sbatch",
    "job_array_high_memory.sbatch",
    "job_array_long_time.sbatch",
    "job_array_extreme.sbatch"
  ),
  cpus = c(4L, 4L, 4L, 4L, 4L),
  memory = c("12g", "16g", "24g", "20g", "32g"),
  walltime = rep("01:00:00", 5L),
  description = c(
    "Baseline cells with observed-age ML comparators.",
    "Large training size, dense visits, or JMbayes2-enabled baseline cells.",
    "Distributional biomarker cells with prior OOM risk.",
    "Severe-error cells with prior walltime risk.",
    "Combined memory/time risk from dense registration, JM, and DeepSurv."
  ),
  stringsAsFactors = FALSE
)

resource_flag_summary <- function() {
  c(
    if (isTRUE(ENABLE_JMBAYES2)) "JMbayes2 enabled",
    if (DEEPSURV_EPOCHS > 80L) paste0("DeepSurv epochs=", DEEPSURV_EPOCHS),
    if (RSF_NUM_TREES > 300L) paste0("RSF trees=", RSF_NUM_TREES),
    if (TIMEERROR_N_IMPUTE > 5L) paste0("time-error imputations=", TIMEERROR_N_IMPUTE),
    if (N_TEST > 1000L) paste0("n_test=", N_TEST)
  )
}

base_resource_class <- function(plan) {
  scenario <- plan$scenario
  high_memory <- scenario %in% c("dist_large", "dist_strong_dense") |
    (grepl("^dist_", scenario) & plan$n_visits >= 12L)
  long_time <- (
    plan$n_train >= 1000L &
      scenario %in% c("mean_severe", "mean_strong_dense", "heavy_tail", "biased_shift", "weak_stage")
  ) | (
    plan$phase == "timepoint_extension" &
      plan$schedule == "m20" &
      scenario %in% c("mean_severe", "mean_strong_dense", "heavy_tail")
  )
  medium <- plan$n_train >= 1000L |
    plan$n_visits >= 12L |
    scenario %in% c(
      "mean_severe", "mean_strong_dense", "dist_moderate", "dist_large",
      "dist_strong_dense", "heavy_tail", "asymmetric_shift", "irregular_missing"
    )

  out <- ifelse(
    high_memory & (long_time | plan$n_train >= 1000L | plan$n_visits >= 20L),
    "extreme",
    ifelse(
      high_memory,
      "high_memory",
      ifelse(long_time, "long_time", ifelse(medium, "medium", "standard"))
    )
  )
  out
}

promote_for_methods <- function(class, plan) {
  out <- class

  if (isTRUE(ENABLE_JMBAYES2)) {
    out <- ifelse(out == "standard", "medium", out)
    out <- ifelse(
      out == "medium" & (plan$n_train >= 1000L | plan$n_visits >= 12L),
      "long_time",
      out
    )
    out <- ifelse(out == "high_memory" & (plan$n_train >= 1000L | plan$n_visits >= 12L), "extreme", out)
  }

  heavier_tuning <- DEEPSURV_EPOCHS > 80L || RSF_NUM_TREES > 300L ||
    TIMEERROR_N_IMPUTE > 5L || N_TEST > 1000L
  if (heavier_tuning) {
    out <- ifelse(out == "standard", "medium", out)
    out <- ifelse(out == "medium", "long_time", out)
    out <- ifelse(out %in% c("high_memory", "long_time"), "extreme", out)
  }

  out
}

resource_class_for_plan <- function(plan) {
  promote_for_methods(base_resource_class(plan), plan)
}

resource_reason <- function(row) {
  reasons <- character(0)
  if (row$scenario %in% c("dist_large", "dist_strong_dense")) {
    reasons <- c(reasons, "distributional biomarker scenario with prior OOM risk")
  }
  if (row$scenario %in% c("mean_severe", "mean_strong_dense", "heavy_tail", "biased_shift", "weak_stage")) {
    reasons <- c(reasons, "severe-error scenario with prior walltime risk")
  }
  if (row$n_train >= 1000L) reasons <- c(reasons, "n_train >= 1000")
  if (row$n_visits >= 12L) reasons <- c(reasons, "dense visit schedule")
  if (row$n_visits >= 20L) reasons <- c(reasons, "20-visit timepoint extension")
  method_reasons <- resource_flag_summary()
  if (length(method_reasons)) reasons <- c(reasons, method_reasons)
  if (!length(reasons)) reasons <- "baseline resource class"
  paste(unique(reasons), collapse = "; ")
}

write_resource_task_lists <- function(out_dir = OUT_DIR) {
  plan <- build_task_plan()
  plan$resource_class <- resource_class_for_plan(plan)
  plan$resource_reason <- vapply(
    seq_len(nrow(plan)),
    function(i) resource_reason(plan[i, , drop = FALSE]),
    character(1)
  )
  plan <- merge(plan, RESOURCE_CLASSES, by = "resource_class", all.x = TRUE, sort = FALSE)
  plan <- plan[order(plan$task_id), , drop = FALSE]

  summary_dir <- file.path(out_dir, "summary")
  task_list_dir <- file.path(summary_dir, "resource_task_lists")
  if (!dir.exists(task_list_dir)) dir.create(task_list_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(plan, file.path(summary_dir, "prediction_resource_task_plan.csv"), row.names = FALSE)

  groups <- do.call(rbind, lapply(split(plan, plan$resource_class), function(z) {
    class <- z$resource_class[1]
    resource <- RESOURCE_CLASSES[RESOURCE_CLASSES$resource_class == class, , drop = FALSE]
    task_list_file <- file.path(task_list_dir, paste0(class, ".txt"))
    writeLines(as.character(z$task_id), task_list_file, useBytes = TRUE)
    data.frame(
      resource_class = class,
      n_tasks = nrow(z),
      sbatch_file = resource$sbatch_file,
      task_list_file = task_list_file,
      cpus = resource$cpus,
      memory = resource$memory,
      walltime = resource$walltime,
      description = resource$description,
      stringsAsFactors = FALSE
    )
  }))
  groups <- groups[match(RESOURCE_CLASSES$resource_class, groups$resource_class, nomatch = 0L), , drop = FALSE]
  write.table(
    groups,
    file.path(summary_dir, "prediction_resource_groups.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  groups
}
