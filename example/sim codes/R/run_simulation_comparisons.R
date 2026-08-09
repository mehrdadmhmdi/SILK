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
  source(file.path("R", "methods_landmark_recorded.R"))
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

fit_one_prediction_method <- function(method, train, cell, seed_base,
                                      shared_registration = NULL) {
  require_shared_registration <- function() {
    if (inherits(shared_registration, "silk_error")) {
      stop("Shared SILK registration failed: ", shared_registration$error_message,
           call. = FALSE)
    }
    if (is.null(shared_registration)) {
      stop("Shared SILK registration was not prepared.", call. = FALSE)
    }
    shared_registration
  }
  switch(
    method,
    "Landmark-Recorded" = fit_landmark_recorded(train$subjects, train$visits),
    "Cox-SameFeature-Recorded" = fit_same_feature_recorded_cox(train$subjects, train$visits),
    "MMLM-Recorded" = fit_mmlm_recorded(train$subjects, train$visits),
    "JM-Recorded" = fit_jm_recorded(train$subjects, train$visits),
    "Bayesian-Dynamic-Observed" = fit_bayesian_dynamic_observed(
      train$subjects, train$visits, seed = seed_base + 607L
    ),
    "DeepSurv-Observed" = fit_deepsurv_observed(
      train$subjects, train$visits, seed = seed_base + 709L
    ),
    "RSF-Observed" = fit_rsf_observed(
      train$subjects, train$visits, seed = seed_base + 811L
    ),
    "TimeError-Integrated-Landmark" = fit_timeerror_integrated_landmark(
      train$subjects, train$visits, seed = seed_base + 503L
    ),
    "SILK" = fit_silk(
      train$subjects, train$visits,
      method = "SILK", registration = require_shared_registration()
    ),
    "Beran-Recorded" = fit_beran_recorded(train$subjects),
    "Beran-SILK" = fit_beran_silk(
      train$subjects, train$visits,
      registration = require_shared_registration()
    ),
    "Beran-Oracle-Latent-Age" = fit_beran_oracle_latent_age(train$subjects),
    "Oracle-Latent-Age" = fit_oracle_latent_age(train$subjects, train$visits),
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
    "Landmark-Recorded" = do.call(predict_landmark_recorded, c(list(fit, test$subjects, test$visits), args)),
    "Cox-SameFeature-Recorded" = do.call(predict_same_feature_recorded_cox, c(list(fit, test$subjects, test$visits), args)),
    "MMLM-Recorded" = do.call(predict_mmlm_recorded, c(list(fit, test$subjects, test$visits), args)),
    "JM-Recorded" = do.call(predict_jm_recorded, c(list(fit, test$subjects, test$visits), args)),
    "Bayesian-Dynamic-Observed" = do.call(predict_bayesian_dynamic_observed, c(list(fit, test$subjects, test$visits), args)),
    "DeepSurv-Observed" = do.call(predict_deepsurv_observed, c(list(fit, test$subjects, test$visits), args)),
    "RSF-Observed" = do.call(predict_rsf_observed, c(list(fit, test$subjects, test$visits), args)),
    "TimeError-Integrated-Landmark" = do.call(predict_timeerror_integrated_landmark, c(list(fit, test$subjects, test$visits), args)),
    "SILK" = do.call(predict_silk, c(list(fit, test$subjects, test$visits), args)),
    "Beran-Recorded" = do.call(predict_beran_recorded, c(list(fit, test$subjects), args)),
    "Beran-SILK" = do.call(predict_beran_silk, c(list(fit, test$subjects, test$visits), args)),
    "Beran-Oracle-Latent-Age" = do.call(predict_beran_oracle_latent_age, c(list(fit, test$subjects), args)),
    "Oracle-Latent-Age" = do.call(predict_oracle_latent_age, c(list(fit, test$subjects, test$visits), args)),
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
                                shared_registration = NULL,
                                registration_seconds = 0,
                                registration_warnings = "") {
  fit_start <- proc.time()[3]
  fit_res <- run_with_warnings(fit_one_prediction_method(
    method, train, cell, seed_base, shared_registration
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
  registration_methods <- intersect(methods, c("SILK", "Beran-SILK"))
  shared_registration <- NULL
  registration_seconds <- 0
  registration_warnings <- ""
  if (length(registration_methods)) {
    registration_start <- proc.time()[3]
    registration_result <- run_with_warnings(fit_silk_registration(
      train$subjects, train$visits,
      shift_grid = make_shift_grid(cell$scenario),
      seed = seed_base + 211L
    ))
    registration_seconds <- proc.time()[3] - registration_start
    registration_warnings <- registration_result$warnings
    shared_registration <- if (registration_result$ok) {
      registration_result$value
    } else {
      structure(
        list(error_message = registration_result$error_message),
        class = "silk_error"
      )
    }
  }
  cores <- max(1L, min(silk_method_cores(), length(methods)))
  one <- function(m) tryCatch(
    run_one_method_full(
      m, train, test, cell, seed_base,
      shared_registration = shared_registration,
      registration_seconds = if (m %in% registration_methods) registration_seconds else 0,
      registration_warnings = if (m %in% registration_methods) registration_warnings else ""
    ),
    error = function(e) list(method = m, pred = NULL,
      status = method_status_row(cell = cell, method = m, fit_ok = FALSE, predict_ok = FALSE,
        fit_seconds = NA_real_, predict_seconds = NA_real_, fit_error = conditionMessage(e),
        predict_error = "method crashed", fit_warnings = "", predict_warnings = "",
        implementation = "",
        registration_seconds = if (m %in% registration_methods) registration_seconds else 0,
        registration_warnings = if (m %in% registration_methods) registration_warnings else "")))
  if (cores > 1L && .Platform$OS.type == "unix" && requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(methods, one, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(methods, one)
  }
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
    runtime_seconds = runtime_seconds
  )
}
