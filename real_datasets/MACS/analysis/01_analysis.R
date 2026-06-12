# =============================================================================
# 01_analysis.R
# 5-fold CV risk prediction on MACS PDS: SILK vs Landmark Cox vs MMLM
# =============================================================================

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
                                  n_train_setting = NA_integer_) {
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
        fold_id = fold_id, replicate_id = 1L,
        n_train_setting = n_train_setting,
        time_grid_setting = "macs_cv", stringsAsFactors = FALSE
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
METHOD_LIST <- c("Landmark-Recorded", "MMLM-Recorded", "SILK")

ids   <- sort(unique(subjects$id))
folds <- sample(rep(seq_len(K), length.out = length(ids)))
names(folds) <- ids

all_predictions <- list()
all_status      <- list()
counter         <- 0L

cat("\n=== ", K, "-fold cross-validation ===\n")
cat("Methods:", paste(METHOD_LIST, collapse = ", "), "\n")
cat("Horizons:", paste(HORIZONS, collapse = ", "), "years\n\n")
cat("SILK kernel:", silk_opt("REGISTRATION_KERNEL"),
    "| approximation:", silk_opt("REGISTRATION_KERNEL_APPROX"),
    "| RFF dim:", silk_opt("KERNEL_RFF_DIM"), "\n\n")

for (ff in seq_len(K)) {
  cat("── Fold", ff, "/", K, "──\n")
  train_ids <- ids[folds != ff]
  test_ids  <- ids[folds == ff]

  train_s <- subjects[subjects$id %in% train_ids, ]
  train_v <- visits[visits$id %in% train_ids, ]
  test_s  <- subjects[subjects$id %in% test_ids, ]
  test_v  <- visits[visits$id %in% test_ids, ]

  cat("  Train:", nrow(train_s), "subj,", nrow(train_v), "visits |",
      "Test:", nrow(test_s), "subj,", nrow(test_v), "visits\n")

  for (method in METHOD_LIST) {
    cat("  ", method, "... ")
    t0 <- proc.time()[3]

    # ── FIT ──
    fit_res <- tryCatch({
      fit <- switch(method,
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
          list(method = method, marker_model = mm,
               fit = fit_residual_cox(train_s, x))
        },
        "SILK" = fit_silk(train_s, train_v,
                          shift_range = MACS_SHIFT_RANGE,
                          kernel = silk_opt("REGISTRATION_KERNEL"),
                          kernel_approx = silk_opt("REGISTRATION_KERNEL_APPROX"),
                          rff_dim = silk_opt("KERNEL_RFF_DIM"),
                          rff_seed = silk_opt("KERNEL_RFF_SEED"),
                          seed = 42L + ff)
      )
      list(ok = TRUE, value = fit)
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))

    fit_time <- proc.time()[3] - t0

    if (!fit_res$ok) {
      cat("FIT FAILED (", fit_res$error, ")\n")
      all_status[[paste0(ff, "_", method)]] <- data.frame(
        fold = ff, method = method, fit_ok = FALSE, predict_ok = FALSE,
        fit_time = fit_time, predict_time = NA, error = fit_res$error,
        stringsAsFactors = FALSE)
      next
    }

    # ── PREDICT ──
    pred_t0 <- proc.time()[3]
    pred_res <- tryCatch({
      pred <- switch(method,
        "Landmark-Recorded" = {
          x <- SILK:::landmark_covariates(test_s, test_v, clock = "recorded",
                                   include_biomarker = FALSE)
          risk <- predict_residual_cox_risk(fit_res$value$fit, x, HORIZONS)
          real_prediction_frame(test_s, test_s$A_obs, HORIZONS, risk,
                                method, fold_id = ff,
                                n_train_setting = nrow(train_s))
        },
        "MMLM-Recorded" = {
          mh <- SILK:::predict_current_value_mixed_model(fit_res$value$marker_model,
                                                   test_s, test_v)
          x  <- SILK:::landmark_covariates(test_s, test_v, clock = "recorded",
                                    marker = mh)
          risk <- predict_residual_cox_risk(fit_res$value$fit, x, HORIZONS)
          real_prediction_frame(test_s, test_s$A_obs, HORIZONS, risk,
                                method, fold_id = ff,
                                n_train_setting = nrow(train_s))
        },
        "SILK" = {
          grid <- fit_res$value$grid
          ps   <- SILK:::predict_registration_shift(
                    fit_res$value$registration$final_template, test_v, grid)
          stage <- test_s$A_obs[match(ps$id, test_s$id)] - ps$e_hat
          stage <- stage[match(test_s$id, ps$id)]
          hist  <- make_history_features(test_s, test_v)
          x     <- SILK:::silk_history_covariates(test_s, hist, stage)
          risk  <- predict_residual_cox_risk(fit_res$value$fit, x, HORIZONS)
          real_prediction_frame(test_s, stage, HORIZONS, risk,
                                method, fold_id = ff,
                                n_train_setting = nrow(train_s))
        }
      )
      list(ok = TRUE, value = pred)
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))

    predict_time <- proc.time()[3] - pred_t0
    cat(sprintf("%.1fs fit, %.1fs predict", fit_time, predict_time))

    if (pred_res$ok) {
      cat(" done\n")
      counter <- counter + 1L
      all_predictions[[counter]] <- pred_res$value
    } else {
      cat(" PREDICT FAILED (", pred_res$error, ")\n")
    }

    all_status[[paste0(ff, "_", method)]] <- data.frame(
      fold = ff, method = method, fit_ok = TRUE,
      predict_ok = pred_res$ok, fit_time = fit_time,
      predict_time = predict_time,
      error = if (!pred_res$ok) pred_res$error else "",
      stringsAsFactors = FALSE)
  }
  cat("\n")
}

# ── 4. Save ──────────────────────────────────────────────────────────────────
predictions <- do.call(rbind, all_predictions)
rownames(predictions) <- NULL
status_df <- do.call(rbind, all_status)
rownames(status_df) <- NULL

write.csv(predictions, file.path(RESULTS_DIR, "macs_predictions.csv"),
          row.names = FALSE)
write.csv(status_df, file.path(RESULTS_DIR, "macs_method_status.csv"),
          row.names = FALSE)

cat("=== Cross-validation complete ===\n")
cat("Predictions:", nrow(predictions), "rows\n")
cat("Saved to:   ", RESULTS_DIR, "\n\n")

cat("Method status:\n")
print(table(status_df$method, status_df$predict_ok))
cat("\nAt-risk counts per method/horizon:\n")
print(xtabs(at_risk ~ method + horizon, data = predictions))
