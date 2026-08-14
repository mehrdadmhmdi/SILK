# =============================================================================
# run_simulation_comparisons.R
# One simulation replication with shared train/test data, shared landmarks,
# shared horizons, and a single prediction-evaluation pipeline.
# =============================================================================

source_silk_prediction_modules <- function() {
  if (!requireNamespace("SILK", quietly = TRUE)) {
    stop(
      "The SILK package is not installed in this R library. Run the cluster ",
      "package-preparation script before task.R.",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(library(SILK))
  source(file.path("R", "cfg.R"))
  source(file.path("R", "methods_common.R"))
  source(file.path("R", "methods_mmlm_recorded.R"))
  source(file.path("R", "methods_jm_recorded.R"))
  source(file.path("R", "methods_observed_ml.R"))
  source(file.path("R", "methods_timeerror_integrated.R"))
}

run_with_warnings <- function(expr) {
  warnings <- character(0)
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) structure(list(error_message = conditionMessage(e)), class = "silk_error")
  )
  list(
    value = value,
    ok = !inherits(value, "silk_error"),
    error_message = if (inherits(value, "silk_error")) value$error_message else "",
    warnings = paste(unique(warnings), collapse = " | ")
  )
}

registration_kernel_for_method <- function(method) {
  switch(
    method,
    "SILK-Gaussian" = "gaussian",
    "SILK-Laplace" = "laplace",
    "SILK-Matern32" = "matern32",
    NA_character_
  )
}

fit_one_prediction_method <- function(method, train, cell, seed_base,
                                      shared_registrations = list()) {
  require_shared_registration <- function(kernel = registration_kernel_for_method(method)) {
    shared_registration <- shared_registrations[[kernel]]
    if (inherits(shared_registration, "silk_error")) {
      stop("Shared SILK ", kernel, " registration failed: ", shared_registration$error_message,
           call. = FALSE)
    }
    if (is.null(shared_registration)) {
      stop("Shared SILK registration was not prepared.", call. = FALSE)
    }
    shared_registration
  }
  switch(
    method,
    "Recorded-Cox" = fit_recorded_cox(train$subjects),
    "Recorded-Beran" = fit_recorded_beran(train$subjects),
    "MMLM-Correct" = fit_mmlm_recorded(
      train$subjects, train$visits, specification = "correct"
    ),
    "MMLM-Misspecified" = fit_mmlm_recorded(
      train$subjects, train$visits, specification = "misspecified"
    ),
    "JM-Correct" = fit_jm_recorded(
      train$subjects, train$visits, specification = "correct"
    ),
    "JM-Misspecified" = fit_jm_recorded(
      train$subjects, train$visits, specification = "misspecified"
    ),
    "DeepSurv" = fit_deepsurv_observed(
      train$subjects, train$visits, seed = seed_base + 709L
    ),
    "RSF" = fit_rsf_observed(
      train$subjects, train$visits, seed = seed_base + 811L
    ),
    "TimeError-Integrated-Landmark" = fit_timeerror_integrated_landmark(
      train$subjects, train$visits, seed = seed_base + 503L
    ),
    "SILK-Gaussian" = fit_silk(
      train$subjects, train$visits,
      method = "SILK-Gaussian", registration = require_shared_registration()
    ),
    "SILK-Laplace" = fit_silk(
      train$subjects, train$visits,
      method = "SILK-Laplace", registration = require_shared_registration("laplace")
    ),
    "SILK-Matern32" = fit_silk(
      train$subjects, train$visits,
      method = "SILK-Matern32", registration = require_shared_registration("matern32")
    ),
    "Oracle-Cox" = fit_oracle_cox(train$subjects),
    "Oracle-Beran" = fit_oracle_beran(train$subjects),
    stop("Unknown method: ", method, call. = FALSE)
  )
}

predict_one_prediction_method <- function(method, fit, test, cell, horizons = PREDICTION_HORIZONS) {
  time_grid_setting <- paste(cell$phase, cell$scenario, cell$schedule, sep = ":")
  args <- list(
    horizons = horizons,
    fold_id = 1L,
    replicate_id = cell$rep,
    n_train_setting = cell$n_train,
    time_grid_setting = time_grid_setting
  )

  pred <- switch(
    method,
    "Recorded-Cox" = do.call(predict_recorded_cox, c(list(fit, test$subjects), args)),
    "Recorded-Beran" = do.call(predict_recorded_beran, c(list(fit, test$subjects), args)),
    "MMLM-Correct" = do.call(predict_mmlm_recorded, c(list(fit, test$subjects, test$visits), args)),
    "MMLM-Misspecified" = do.call(predict_mmlm_recorded, c(list(fit, test$subjects, test$visits), args)),
    "JM-Correct" = do.call(predict_jm_recorded, c(list(fit, test$subjects, test$visits), args)),
    "JM-Misspecified" = do.call(predict_jm_recorded, c(list(fit, test$subjects, test$visits), args)),
    "DeepSurv" = do.call(predict_deepsurv_observed, c(list(fit, test$subjects, test$visits), args)),
    "RSF" = do.call(predict_rsf_observed, c(list(fit, test$subjects, test$visits), args)),
    "TimeError-Integrated-Landmark" = do.call(predict_timeerror_integrated_landmark, c(list(fit, test$subjects, test$visits), args)),
    "SILK-Gaussian" = do.call(predict_silk, c(list(fit, test$subjects, test$visits), args)),
    "SILK-Laplace" = do.call(predict_silk, c(list(fit, test$subjects, test$visits), args)),
    "SILK-Matern32" = do.call(predict_silk, c(list(fit, test$subjects, test$visits), args)),
    "Oracle-Cox" = do.call(predict_oracle_cox, c(list(fit, test$subjects), args)),
    "Oracle-Beran" = do.call(predict_oracle_beran, c(list(fit, test$subjects), args)),
    stop("Unknown method: ", method, call. = FALSE)
  )
  validate_prediction_frame(pred)
}

method_status_row <- function(cell, method, fit_ok, predict_ok, fit_seconds, predict_seconds,
                              fit_error = "", predict_error = "", fit_warnings = "",
                              predict_warnings = "", implementation = "",
                              registration_seconds = 0,
                              registration_warnings = "") {
  data.frame(
    task_id = cell$task_id,
    cell_id = cell$cell_id,
    phase = cell$phase,
    scenario = cell$scenario,
    schedule = cell$schedule,
    n_visits = cell$n_visits,
    n_train = cell$n_train,
    n_test = cell$n_test,
    rep = cell$rep,
    method = method,
    fit_ok = fit_ok,
    predict_ok = predict_ok,
    success = fit_ok && predict_ok,
    fit_seconds = fit_seconds,
    predict_seconds = predict_seconds,
    registration_seconds = registration_seconds,
    fit_error = fit_error,
    predict_error = predict_error,
    fit_warnings = fit_warnings,
    predict_warnings = predict_warnings,
    registration_warnings = registration_warnings,
    implementation = implementation,
    stringsAsFactors = FALSE
  )
}

# Run one method end-to-end (fit + predict + status row). Safe to call in parallel.
run_one_method_full <- function(method, train, test, cell, seed_base,
                                shared_registrations = list(),
                                registration_seconds = 0,
                                registration_warnings = "") {
  fit_start <- proc.time()[3]
  fit_res <- run_with_warnings(fit_one_prediction_method(
    method, train, cell, seed_base, shared_registrations
  ))
  fit_seconds <- proc.time()[3] - fit_start
  predict_seconds <- NA_real_
  pred_res <- list(ok = FALSE, value = NULL, error_message = "Fit failed.", warnings = "")
  implementation <- ""
  pred_value <- NULL
  if (fit_res$ok) {
    implementation <- if (!is.null(fit_res$value$implementation)) fit_res$value$implementation else ""
    pred_start <- proc.time()[3]
    pred_res <- run_with_warnings(predict_one_prediction_method(method, fit_res$value, test, cell, PREDICTION_HORIZONS))
    predict_seconds <- proc.time()[3] - pred_start
    if (pred_res$ok) pred_value <- pred_res$value
  }
  status <- method_status_row(
    cell = cell, method = method,
    fit_ok = fit_res$ok, predict_ok = pred_res$ok,
    fit_seconds = fit_seconds, predict_seconds = predict_seconds,
    fit_error = fit_res$error_message, predict_error = pred_res$error_message,
    fit_warnings = fit_res$warnings, predict_warnings = pred_res$warnings,
    implementation = implementation,
    registration_seconds = registration_seconds,
    registration_warnings = registration_warnings
  )
  list(method = method, pred = pred_value, status = status)
}

# Number of CPU workers for method-level parallelism. Defaults to
# SLURM_CPUS_PER_TASK; override with SILK_METHOD_CORES (lower it to bound memory).
silk_method_cores <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("SILK_METHOD_CORES", unset = "")))
  if (is.na(n) || n < 1L) {
    n <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
    if (is.na(n) || n < 1L) n <- 1L
  }
  max(1L, min(n, length(METHOD_ORDER)))
}

# Fit/predict all methods, in parallel across CPUs (fork-based) when available.
# Each stochastic method seeds itself from seed_base, so parallel results are
# identical to sequential; only wall-clock time changes.
# Restrict to a method subset when SILK_ONLY_METHODS is set (comma-separated,
# matched against METHOD_ORDER). This lets the cluster add new methods to an
# existing run, or finish missing cells, without recomputing valid results.
methods_to_run <- function() {
  only <- Sys.getenv("SILK_ONLY_METHODS", unset = "")
  if (!nzchar(only)) return(METHOD_ORDER)
  requested <- trimws(strsplit(only, ",")[[1]])
  requested <- requested[nzchar(requested)]
  unknown <- setdiff(requested, METHOD_ORDER)
  if (length(unknown)) {
    stop("SILK_ONLY_METHODS lists unknown methods: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  intersect(METHOD_ORDER, requested)
}

run_all_methods <- function(train, test, cell, seed_base) {
  methods <- methods_to_run()
  registration_methods <- methods[!is.na(vapply(
    methods, registration_kernel_for_method, character(1)
  ))]
  required_kernels <- unique(vapply(
    registration_methods, registration_kernel_for_method, character(1)
  ))
  shared_registrations <- list()
  registration_seconds <- stats::setNames(numeric(length(required_kernels)), required_kernels)
  registration_warnings <- stats::setNames(rep("", length(required_kernels)), required_kernels)
  for (kernel in required_kernels) {
    registration_start <- proc.time()[3]
    registration_result <- run_with_warnings(fit_silk_registration(
      train$subjects, train$visits,
      shift_grid = make_shift_grid(cell$scenario),
      seed = seed_base + 211L,
      biomarker_kernel = kernel
    ))
    registration_seconds[[kernel]] <- proc.time()[3] - registration_start
    registration_warnings[[kernel]] <- registration_result$warnings
    shared_registrations[[kernel]] <- if (registration_result$ok) {
      registration_result$value
    } else {
      structure(list(error_message = registration_result$error_message), class = "silk_error")
    }
  }
  cores <- max(1L, min(silk_method_cores(), length(methods)))
  one <- function(m) tryCatch(
    {
      kernel <- registration_kernel_for_method(m)
      method_registration_seconds <- if (is.na(kernel)) 0 else registration_seconds[[kernel]]
      method_registration_warnings <- if (is.na(kernel)) "" else registration_warnings[[kernel]]
    run_one_method_full(
      m, train, test, cell, seed_base,
      shared_registrations = shared_registrations,
      registration_seconds = method_registration_seconds,
      registration_warnings = method_registration_warnings
    )
    },
    error = function(e) list(method = m, pred = NULL,
      status = method_status_row(cell = cell, method = m, fit_ok = FALSE, predict_ok = FALSE,
        fit_seconds = NA_real_, predict_seconds = NA_real_, fit_error = conditionMessage(e),
        predict_error = "method crashed", fit_warnings = "", predict_warnings = "",
        implementation = "",
        registration_seconds = 0,
        registration_warnings = "")))
  if (cores > 1L && .Platform$OS.type == "unix" && requireNamespace("parallel", quietly = TRUE)) {
    result <- parallel::mclapply(methods, one, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    result <- lapply(methods, one)
  }
  attr(result, "registrations") <- shared_registrations
  result
}

registration_metric_row <- function(registration, subjects, visits, cell,
                                    sample_role = c("train_crossfit", "test_heldout")) {
  sample_role <- match.arg(sample_role)
  if (sample_role == "train_crossfit") {
    prediction <- registration$train_stage
    prediction <- prediction[match(subjects$id, prediction$id), , drop = FALSE]
  } else {
    prediction <- predict_silk_registration(registration, visits)
    prediction <- prediction[match(subjects$id, prediction$id), , drop = FALSE]
  }
  e_hat <- as.numeric(prediction$e_hat)
  e_true <- as.numeric(subjects$eps)
  stage_hat <- as.numeric(subjects$A_obs) - e_hat
  stage_true <- as.numeric(subjects$A_star)
  finite <- is.finite(e_hat) & is.finite(e_true) & is.finite(stage_hat) & is.finite(stage_true)
  safe_cor <- function(x, y, method = "pearson") {
    if (sum(is.finite(x) & is.finite(y)) < 3L || stats::sd(x, na.rm = TRUE) < 1e-12 ||
        stats::sd(y, na.rm = TRUE) < 1e-12) return(NA_real_)
    suppressWarnings(stats::cor(x, y, method = method, use = "complete.obs"))
  }
  slope <- if (sum(finite) >= 3L && stats::sd(e_true[finite]) > 1e-12) {
    unname(stats::coef(stats::lm(e_hat[finite] ~ e_true[finite]))[2L])
  } else NA_real_
  stage_variance <- sum((stage_true[finite] - mean(stage_true[finite]))^2)
  stage_r2 <- if (sum(finite) >= 2L && stage_variance > 0) {
    1 - sum((stage_hat[finite] - stage_true[finite])^2) / stage_variance
  } else NA_real_
  recorded_rmse <- sqrt(mean((subjects$A_obs[finite] - stage_true[finite])^2))
  stage_rmse <- sqrt(mean((stage_hat[finite] - stage_true[finite])^2))
  builder <- registration$final_template$biomarker_builder

  # Optimization-stability diagnostics. ms is the full-training-set multistart;
  # fold_ms summarizes the cross-fitting folds that actually produced the
  # cross-fitted calibrated states used by the survival layer.
  ms <- registration$multistart
  if (is.null(ms)) ms <- list()
  ms_value <- function(name, default = NA_real_) {
    z <- ms[[name]]
    if (is.null(z) || !length(z)) default else z[[1L]]
  }
  fold_ms <- registration$fold_multistart
  fold_mean <- function(name) {
    if (is.null(fold_ms) || !nrow(fold_ms) || !name %in% names(fold_ms)) return(NA_real_)
    mean(as.numeric(fold_ms[[name]]), na.rm = TRUE)
  }
  data.frame(
    task_id = cell$task_id,
    phase = cell$phase,
    scenario = cell$scenario,
    schedule = cell$schedule,
    n_train = cell$n_train,
    rep = cell$rep,
    sample_role = sample_role,
    biomarker_kernel = builder$kernel,
    biomarker_bandwidth = builder$bandwidth,
    bandwidth_rule = builder$bandwidth_rule,
    clock_signal_r2 = registration$final_template$clock_signal_r2,
    longitudinal_kernel_signal = registration$final_template$longitudinal_signal,
    shift_correlation = safe_cor(e_hat, e_true),
    shift_spearman = safe_cor(e_hat, e_true, "spearman"),
    shift_rmse = sqrt(mean((e_hat[finite] - e_true[finite])^2)),
    shift_mae = mean(abs(e_hat[finite] - e_true[finite])),
    shift_bias = mean(e_hat[finite] - e_true[finite]),
    shift_slope = slope,
    shift_sd_ratio = if (stats::sd(e_true[finite]) > 1e-12) {
      stats::sd(e_hat[finite]) / stats::sd(e_true[finite])
    } else NA_real_,
    stage_rmse = stage_rmse,
    recorded_stage_rmse = recorded_rmse,
    stage_rmse_improvement = recorded_rmse - stage_rmse,
    stage_r2 = stage_r2,
    boundary_fraction = mean(as.logical(prediction$at_boundary), na.rm = TRUE),
    median_gap_q1 = stats::median(prediction$gap_q1, na.rm = TRUE),
    median_gap_q2 = stats::median(prediction$gap_q2, na.rm = TRUE),
    estimated_shift_sd = stats::sd(e_hat[finite]),
    mean_estimated_shift = mean(e_hat[finite]),
    shift_x1_correlation = safe_cor(e_hat, subjects$X1),
    n_evaluated = sum(finite),
    n_starts = ms_value("n_starts"),
    n_successful_starts = ms_value("n_successful_starts"),
    best_start = ms_value("best_start", NA_character_),
    best_objective = ms_value("best_objective"),
    objective_spread = ms_value("objective_spread"),
    second_best_gap = ms_value("second_best_gap"),
    start_names = ms_value("start_names", NA_character_),
    start_objective_values = ms_value("start_objective_values", NA_character_),
    fold_mean_successful_starts = fold_mean("n_successful_starts"),
    fold_mean_objective_spread = fold_mean("objective_spread"),
    fold_mean_second_best_gap = fold_mean("second_best_gap"),
    stringsAsFactors = FALSE
  )
}

collect_registration_diagnostics <- function(registrations, train, test, cell) {
  rows <- list()
  for (kernel in names(registrations)) {
    registration <- registrations[[kernel]]
    if (inherits(registration, "silk_error")) next
    rows[[paste(kernel, "train", sep = ":")]] <- registration_metric_row(
      registration, train$subjects, train$visits, cell, "train_crossfit"
    )
    rows[[paste(kernel, "test", sep = ":")]] <- registration_metric_row(
      registration, test$subjects, test$visits, cell, "test_heldout"
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

run_simulation_comparisons <- function(cell) {
  # Seed policy. Primary phase: unique per (task, rep). Timepoint-extension
  # phase: COMMON RANDOM NUMBERS across visit schedules -- the seed depends only
  # on (scenario, rep), so increasing the number of visits adds information to
  # the SAME subjects and outcomes. This removes the independent-sample noise
  # that otherwise makes the visit-count curves non-monotone.
  if (identical(cell$phase, "timepoint_extension")) {
    scen_idx <- match(cell$scenario, SCENARIOS_ALL)
    if (is.na(scen_idx)) scen_idx <- 0L
    seed_base <- GLOBAL_SEED + 700000000L + 100000L * scen_idx + 1000L * cell$rep
  } else {
    seed_base <- GLOBAL_SEED + 100000L * cell$task_id + 1000L * cell$rep
  }
  t0 <- proc.time()[3]
  train <- generate_dataset_fixed(cell$n_train, cell$scenario, cell$schedule, seed = seed_base + 11L)
  test <- generate_dataset_fixed(cell$n_test, cell$scenario, cell$schedule, seed = seed_base + 29L)

  results <- run_all_methods(train, test, cell, seed_base)
  registration_diagnostics <- collect_registration_diagnostics(
    attr(results, "registrations"), train, test, cell
  )

  pred_list <- list()
  status_list <- list()
  for (r in results) {
    if (!is.list(r) || is.null(r$status)) next  # skip a worker that crashed (e.g. OOM)
    if (!is.null(r$pred)) pred_list[[r$method]] <- r$pred
    status_list[[r$method]] <- r$status
  }

  if (!length(pred_list)) {
    stop("No prediction method completed successfully in task ", cell$task_id, call. = FALSE)
  }

  predictions <- do.call(rbind, pred_list)
  rownames(predictions) <- NULL
  predictions <- validate_prediction_frame(predictions)
  evaluation <- evaluate_prediction_frame(predictions)
  runtime_seconds <- proc.time()[3] - t0
  evaluation$replicate_metrics$runtime_seconds <- runtime_seconds
  method_status <- do.call(rbind, status_list)
  method_status$task_runtime_seconds <- runtime_seconds
  list(
    predictions = predictions,
    per_horizon = evaluation$per_horizon,
    replicate_metrics = evaluation$replicate_metrics,
    summary = evaluation$summary,
    paired = evaluation$paired,
    calibration_bins = evaluation$calibration_bins,
    method_status = method_status,
    registration_diagnostics = registration_diagnostics,
    runtime_seconds = runtime_seconds
  )
}
