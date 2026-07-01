# =============================================================================
# 01_analysis.R
# 5-fold CV risk prediction on MACS PDS: SILK vs Landmark Cox vs MMLM
# =============================================================================
library(SILK)


macs_current_script_dir <- function(default = getwd()) {
  frames <- sys.frames()
  for (ii in rev(seq_along(frames))) {
    ofile <- frames[[ii]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      ofile_norm <- tryCatch(
        normalizePath(ofile, winslash = "/", mustWork = TRUE),
        error = function(e) normalizePath(basename(ofile), winslash = "/", mustWork = TRUE)
      )
      return(dirname(ofile_norm))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])),
                         winslash = "/", mustWork = TRUE))
  }
  normalizePath(default, winslash = "/", mustWork = TRUE)
}

source(file.path(macs_current_script_dir(), "00_setup.R"), chdir = TRUE)

# ── 1. Load MACS data ───────────────────────────────────────────────────────
subjects_raw <- read.csv(file.path(DATA_DIR, "macs_subjects.csv"))
visits_raw   <- read.csv(file.path(DATA_DIR, "macs_visits.csv"))

# Map to SILK column names
subjects <- macs_to_silk_subjects(subjects_raw)
visits   <- macs_to_silk_visits(visits_raw)

cat("MACS PDS data loaded:\n")
cat("  Subjects:", nrow(subjects), "\n")
cat("  Visits:  ", nrow(visits), "\n")
cat("  Events:  ", sum(subjects$delta),
    "(", round(100 * mean(subjects$delta), 1), "%)\n")
cat("  Biomarkers:", paste(bio_columns(visits), collapse = ", "), "\n")

# ── 2. Prediction frame for real data ────────────────────────────────────────
# The package's prediction_frame() uses D_star (simulation only).
# This version handles real censored outcomes properly.

real_prediction_frame <- function(test_subjects, landmark, horizons,
                                  risk_mat, method_name, fold_id = 1L,
                                  n_train_setting = NA_integer_,
                                  landmark_set = "last_visit",
                                  landmark_time = NA_real_) {
  n <- nrow(test_subjects)
  H <- length(horizons)
  risk_mat <- as.matrix(risk_mat)
  rows <- list()
  rr <- 1L

  for (i in seq_len(n)) {
    U_i     <- test_subjects$U[i]
    delta_i <- test_subjects$delta[i]
    for (h in seq_len(H)) {
      tau <- horizons[h]
      if (delta_i == 1 && U_i <= tau) {
        event_wh <- 1L; at_risk <- TRUE
      } else if (delta_i == 0 && U_i < tau) {
        event_wh <- 0L; at_risk <- FALSE   # censored before horizon
      } else {
        event_wh <- 0L; at_risk <- TRUE
      }
      rows[[rr]] <- data.frame(
        subject_id = test_subjects$id[i], landmark = landmark[i],
        horizon = tau, method = method_name,
        risk_pred = clip_probability(risk_mat[i, h]),
        event_within_horizon = event_wh, at_risk = at_risk,
        residual_time = U_i, event_status = delta_i,
        landmark_set = landmark_set, landmark_time = landmark_time,
        fold_id = fold_id, replicate_id = 1L,
        n_train_setting = n_train_setting,
        time_grid_setting = paste0("macs_cv:", landmark_set), stringsAsFactors = FALSE
      )
      rr <- rr + 1L
    }
  }
  do.call(rbind, rows)
}

# ── 3. Cross-validation ─────────────────────────────────────────────────────
set.seed(2024)
K           <- 5L
HORIZONS    <- silk_opt("PREDICTION_HORIZONS")
METHOD_LIST <- silk_opt("METHOD_ORDER")
landmark_sets <- macs_landmark_sets(subjects, visits)

# Extract an event-risk vector at target times from a JMbayes2 prediction object
# (mirrors the simulation's JM-Recorded comparator).
jm_event_risk <- function(pred, target_times) {
  if (is.data.frame(pred)) {
    df <- pred
  } else if (is.list(pred)) {
    dfs <- pred[vapply(pred, is.data.frame, logical(1))]
    if (!length(dfs)) stop("JMbayes2 prediction did not return a data frame.", call. = FALSE)
    df <- dfs[[1]]
  } else stop("Unsupported JMbayes2 prediction object.", call. = FALSE)
  time_col  <- intersect(c("times", "time", "Time", "A_obs_il"), names(df))
  risk_cols <- grep("risk|event|cif|prob", names(df), ignore.case = TRUE, value = TRUE)
  risk_cols <- risk_cols[vapply(df[risk_cols], is.numeric, logical(1))]
  risk_cols <- setdiff(risk_cols, time_col)
  if (!length(risk_cols)) {
    num <- names(df)[vapply(df, is.numeric, logical(1))]
    risk_cols <- setdiff(num, c(time_col, "id", "X1", "X2", "Y_obs", "status", "A_obs_il"))
  }
  if (!length(risk_cols)) stop("Could not identify the JMbayes2 event-risk column.", call. = FALSE)
  risk <- as.numeric(df[[risk_cols[1]]])
  if (length(time_col)) {
    tt  <- as.numeric(df[[time_col[1]]])
    out <- stats::approx(tt, risk, xout = target_times, rule = 2, ties = "ordered")$y
  } else if (length(risk) >= length(target_times)) {
    out <- tail(risk, length(target_times))
  } else stop("JMbayes2 prediction did not include enough event-risk values.", call. = FALSE)
  pmin(pmax(out, 0), 1)
}

fit_macs_method <- function(method, train_s, train_v, fold_seed) {
  switch(method,
    "Landmark-Recorded" = {
      x <- SILK:::landmark_covariates(train_s, train_v, clock = "recorded",
                                      include_biomarker = FALSE)
      list(method = method, fit = fit_residual_cox(train_s, x))
    },
    "MMLM-Recorded" = {
      mm <- SILK:::fit_current_value_mixed_model(train_s, train_v)
      mh <- SILK:::predict_current_value_mixed_model(mm, train_s, train_v)
      x  <- SILK:::landmark_covariates(train_s, train_v, clock = "recorded",
                                       marker = mh)
      list(method = method, marker_model = mm, fit = fit_residual_cox(train_s, x))
    },
    "JM-Recorded" = {
      if (!requireNamespace("JMbayes2", quietly = TRUE))
        stop("JM-Recorded requires the JMbayes2 package.", call. = FALSE)
      mm <- SILK:::fit_current_value_mixed_model(train_s, train_v)
      if (is.null(mm$fit))
        stop("JM-Recorded could not fit the recorded-time longitudinal mixed model.",
             call. = FALSE)
      # Absolute event/censoring time on the recorded disease-age clock:
      # landmark age plus residual follow-up (Y_obs = A_obs + U), status = delta.
      event_dat <- data.frame(
        id     = train_s$id,
        Y_obs  = train_s$A_obs + train_s$U,
        status = as.integer(train_s$delta),
        X1     = train_s$X1,
        X2     = train_s$X2
      )
      cox_fit <- survival::coxph(survival::Surv(Y_obs, status) ~ X1 + X2,
                                 data = event_dat, x = TRUE)
      jm_fit <- JMbayes2::jm(
        cox_fit, mm$fit, time_var = "A_obs_il",
        n_chains = JM_N_CHAINS, n_iter = JM_N_ITER, n_burnin = JM_N_BURNIN
      )
      list(method = method, marker_model = mm, cox_fit = cox_fit, jm_fit = jm_fit)
    },
    "SILK-LinearMMD" = fit_silk(
      train_s, train_v,
      shift_range = MACS_SHIFT_RANGE,
      feature_type = "silk",
      method = "SILK-LinearMMD",
      kernel = "linear",
      kernel_approx = "exact",
      seed = fold_seed
    ),
    "SILK" = fit_silk(
      train_s, train_v,
      shift_range = MACS_SHIFT_RANGE,
      feature_type = "silk",
      method = "SILK",
      kernel = silk_opt("REGISTRATION_KERNEL"),
      kernel_approx = silk_opt("REGISTRATION_KERNEL_APPROX"),
      rff_dim = silk_opt("KERNEL_RFF_DIM"),
      rff_seed = silk_opt("KERNEL_RFF_SEED"),
      seed = fold_seed
    ),
    stop("Unknown MACS method: ", method, call. = FALSE)
  )
}

silk_diagnostics <- function(fit, ps, fold_id, method, landmark_set, landmark_time) {
  train_diag <- fit$registration$train_stage
  train_diag$split <- "train_crossfit"
  test_diag <- ps
  test_diag$S_hat <- NA_real_
  test_diag$fold <- fold_id
  test_diag$split <- "test"
  common <- intersect(names(train_diag), names(test_diag))
  out <- rbind(train_diag[, common, drop = FALSE],
               test_diag[, common, drop = FALSE])
  out$fold_id <- fold_id
  out$method <- method
  out$landmark_set <- landmark_set
  out$landmark_time <- landmark_time
  out
}

predict_macs_method <- function(method, fit, test_s, test_v, fold_id,
                                n_train, landmark_set, landmark_time) {
  if (identical(method, "Landmark-Recorded")) {
    x <- SILK:::landmark_covariates(test_s, test_v, clock = "recorded",
                                    include_biomarker = FALSE)
    risk <- predict_residual_cox_risk(fit$fit, x, HORIZONS)
    return(list(
      predictions = real_prediction_frame(test_s, test_s$A_obs, HORIZONS, risk,
                                          method, fold_id, n_train,
                                          landmark_set, landmark_time),
      diagnostics = NULL
    ))
  }
  if (identical(method, "MMLM-Recorded")) {
    mh <- SILK:::predict_current_value_mixed_model(fit$marker_model, test_s, test_v)
    x  <- SILK:::landmark_covariates(test_s, test_v, clock = "recorded", marker = mh)
    risk <- predict_residual_cox_risk(fit$fit, x, HORIZONS)
    return(list(
      predictions = real_prediction_frame(test_s, test_s$A_obs, HORIZONS, risk,
                                          method, fold_id, n_train,
                                          landmark_set, landmark_time),
      diagnostics = NULL
    ))
  }
  if (identical(method, "JM-Recorded")) {
    b <- fit$marker_model$biomarker
    risk <- matrix(NA_real_, nrow = nrow(test_s), ncol = length(HORIZONS))
    for (i in seq_len(nrow(test_s))) {
      sv <- test_v[test_v$id == test_s$id[i], , drop = FALSE]
      if (!nrow(sv)) next
      sv$marker_value <- sv[[b]]
      sv$id     <- factor(sv$id)
      sv$Y_obs  <- test_s$A_obs[i]   # alive at the landmark; predict forward
      sv$status <- 0L
      sv$X1     <- test_s$X1[i]
      sv$X2     <- test_s$X2[i]
      target_times <- test_s$A_obs[i] + HORIZONS
      pred <- tryCatch(stats::predict(
        fit$jm_fit, newdata = sv, process = "event", times = target_times,
        control = list(cores = 1L, n_samples = JM_PRED_N_SAMPLES, return_newdata = TRUE)
      ), error = function(e) NULL)
      if (!is.null(pred)) risk[i, ] <- jm_event_risk(pred, target_times)
    }
    return(list(
      predictions = real_prediction_frame(test_s, test_s$A_obs, HORIZONS, risk,
                                          method, fold_id, n_train,
                                          landmark_set, landmark_time),
      diagnostics = NULL
    ))
  }

  grid <- fit$grid
  ps <- SILK:::predict_registration_shift(fit$registration$final_template, test_v, grid)
  stage <- test_s$A_obs[match(ps$id, test_s$id)] - ps$e_hat
  stage <- stage[match(test_s$id, ps$id)]
  ps$S_hat <- test_s$A_obs[match(ps$id, test_s$id)] - ps$e_hat
  hist <- make_history_features(test_s, test_v)
  x <- SILK:::silk_history_covariates(test_s, hist, stage)
  risk <- predict_residual_cox_risk(fit$fit, x, HORIZONS)
  list(
    predictions = real_prediction_frame(test_s, stage, HORIZONS, risk,
                                        method, fold_id, n_train,
                                        landmark_set, landmark_time),
    diagnostics = silk_diagnostics(fit, ps, fold_id, method, landmark_set, landmark_time)
  )
}

all_predictions <- list()
all_status      <- list()
all_diagnostics <- list()
counter         <- 0L
diag_counter    <- 0L

cat("\n=== ", K, "-fold cross-validation ===\n")
cat("Methods:", paste(METHOD_LIST, collapse = ", "), "\n")
cat("Landmarks:", paste(names(landmark_sets), collapse = ", "), "\n")
cat("Horizons:", paste(HORIZONS, collapse = ", "), "years\n\n")
cat("SILK kernel:", silk_opt("REGISTRATION_KERNEL"),
    "| approximation:", silk_opt("REGISTRATION_KERNEL_APPROX"),
    "| RFF dim:", silk_opt("KERNEL_RFF_DIM"), "\n\n")

for (landmark_set in names(landmark_sets)) {
  ldat <- landmark_sets[[landmark_set]]
  set_subjects <- ldat$subjects
  set_visits <- ldat$visits
  landmark_time <- ldat$landmark_time
  if (nrow(set_subjects) < K || length(unique(set_subjects$delta)) < 2L) {
    cat("Skipping", landmark_set, "because the risk set is too small or has one outcome class\n")
    next
  }
  ids <- sort(unique(set_subjects$id))
  folds <- sample(rep(seq_len(K), length.out = length(ids)))
  names(folds) <- ids
  cat("\n── Landmark set:", landmark_set, "subjects:", length(ids), "visits:", nrow(set_visits), "──\n")

  for (ff in seq_len(K)) {
    cat("Fold", ff, "/", K, "\n")
    train_ids <- ids[folds != ff]
    test_ids  <- ids[folds == ff]

    train_s <- set_subjects[set_subjects$id %in% train_ids, ]
    train_v <- set_visits[set_visits$id %in% train_ids, ]
    test_s  <- set_subjects[set_subjects$id %in% test_ids, ]
    test_v  <- set_visits[set_visits$id %in% test_ids, ]

    cat("  Train:", nrow(train_s), "subj,", nrow(train_v), "visits |",
        "Test:", nrow(test_s), "subj,", nrow(test_v), "visits\n")

    for (method in METHOD_LIST) {
      cat("  ", method, "... ")
      t0 <- proc.time()[3]
      fit_res <- tryCatch(
        list(ok = TRUE, value = fit_macs_method(method, train_s, train_v, 42L + ff)),
        error = function(e) list(ok = FALSE, error = conditionMessage(e))
      )
      fit_time <- proc.time()[3] - t0
      status_key <- paste(landmark_set, ff, method, sep = "_")

      if (!fit_res$ok) {
        cat("FIT FAILED (", fit_res$error, ")\n")
        all_status[[status_key]] <- data.frame(
          landmark_set = landmark_set, landmark_time = landmark_time,
          fold = ff, method = method, fit_ok = FALSE, predict_ok = FALSE,
          fit_time = fit_time, predict_time = NA, error = fit_res$error,
          stringsAsFactors = FALSE)
        next
      }

      pred_t0 <- proc.time()[3]
      pred_res <- tryCatch(
        list(ok = TRUE, value = predict_macs_method(
          method, fit_res$value, test_s, test_v, ff, nrow(train_s),
          landmark_set, landmark_time
        )),
        error = function(e) list(ok = FALSE, error = conditionMessage(e))
      )
      predict_time <- proc.time()[3] - pred_t0
      cat(sprintf("%.1fs fit, %.1fs predict", fit_time, predict_time))

      if (pred_res$ok) {
        cat(" done\n")
        counter <- counter + 1L
        all_predictions[[counter]] <- pred_res$value$predictions
        if (!is.null(pred_res$value$diagnostics)) {
          diag_counter <- diag_counter + 1L
          all_diagnostics[[diag_counter]] <- pred_res$value$diagnostics
        }
      } else {
        cat(" PREDICT FAILED (", pred_res$error, ")\n")
      }

      all_status[[status_key]] <- data.frame(
        landmark_set = landmark_set, landmark_time = landmark_time,
        fold = ff, method = method, fit_ok = TRUE,
        predict_ok = pred_res$ok, fit_time = fit_time,
        predict_time = predict_time,
        error = if (!pred_res$ok) pred_res$error else "",
        stringsAsFactors = FALSE)
    }
    cat("\n")
  }
}

# ── 4. Save ──────────────────────────────────────────────────────────────────
if (!length(all_predictions)) {
  stop("No MACS predictions were produced; inspect macs_method_status.csv")
}

predictions <- do.call(rbind, all_predictions)
rownames(predictions) <- NULL
status_df <- do.call(rbind, all_status)
rownames(status_df) <- NULL
diagnostics <- if (length(all_diagnostics)) do.call(rbind, all_diagnostics) else data.frame()
if (nrow(diagnostics)) rownames(diagnostics) <- NULL

write.csv(predictions, file.path(RESULTS_DIR, "macs_predictions.csv"),
          row.names = FALSE)
write.csv(status_df, file.path(RESULTS_DIR, "macs_method_status.csv"),
          row.names = FALSE)
write.csv(diagnostics, file.path(RESULTS_DIR, "macs_profile_diagnostics.csv"),
          row.names = FALSE)

cat("=== Cross-validation complete ===\n")
cat("Predictions:", nrow(predictions), "rows\n")
cat("Diagnostics:", nrow(diagnostics), "rows\n")
cat("Saved to:   ", RESULTS_DIR, "\n\n")

cat("Method status:\n")
print(table(status_df$method, status_df$predict_ok))
cat("\nAt-risk counts per method/horizon:\n")
print(xtabs(at_risk ~ method + horizon, data = predictions))
