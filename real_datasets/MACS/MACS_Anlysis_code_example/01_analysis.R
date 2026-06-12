# =============================================================================
# 01_analysis.R
# Main analysis: 5-fold CV risk prediction on MACS PDS using SILK + comparators
# =============================================================================

source("00_setup.R")

# ── 1. Load MACS data ───────────────────────────────────────────────────────
subjects <- read.csv(file.path(DATA_DIR, "macs_subjects.csv"))
visits   <- read.csv(file.path(DATA_DIR, "macs_visits.csv"))

cat("Loaded MACS PDS data:\n")
cat("  Subjects:", nrow(subjects), "\n")
cat("  Visits:  ", nrow(visits), "\n")
cat("  Events:  ", sum(subjects$delta), "(", round(100 * mean(subjects$delta), 1), "%)\n")
cat("  Biomarkers:", paste(bio_columns(visits), collapse = ", "), "\n")

# ── 2. Define the real-data "scenario" for SILK ──────────────────────────────
# SILK's fit_silk() requires a scenario_name to determine the shift grid.
# We inject a custom scenario into the SCENARIOS table.
SCENARIOS <- silk_opt("SCENARIOS")
macs_row <- data.frame(
  scenario         = "macs_real",
  biomarker_signal = "mean",
  default_schedule = "m4",
  sigma_eps        = 1.2,     # estimated from SC interval mean
  eps_mean         = 0,
  eps_type         = "normal",
  n_biomarkers     = 4L,
  missing_rate     = 0,
  irregular        = FALSE,
  shift_min        = -4,
  shift_max        =  4,
  signal_amp       = 1.65,
  sigma_bio        = 0.60,
  u_bio_coef       = 0.12,
  dist_sd_base     = 0.45,
  dist_sd_slope    = 0.18,
  risk_beta_A      = 1.05,
  risk_beta_U      = 0.12,
  description      = "MACS PDS real data",
  stringsAsFactors = FALSE
)
silk_options(SCENARIOS = rbind(SCENARIOS, macs_row))
SCENARIO_NAME <- "macs_real"

# ── 3. Prediction frame for real data ────────────────────────────────────────
# The package's prediction_frame() uses D_star (true residual time) which
# only exists in simulation. For real data, we build prediction frames manually
# using observed (U, delta) with proper censoring handling.

real_prediction_frame <- function(test_subjects, landmark, horizons,
                                  risk_mat, method_name,
                                  fold_id = 1L) {
  n <- nrow(test_subjects)
  H <- length(horizons)
  risk_mat <- as.matrix(risk_mat)

  rows <- vector("list", n)
  for (i in seq_len(n)) {
    U_i     <- test_subjects$U[i]
    delta_i <- test_subjects$delta[i]

    for (h in seq_len(H)) {
      tau <- horizons[h]
      # Determine event-within-horizon and at-risk status
      if (delta_i == 1 && U_i <= tau) {
        event_wh <- 1L
        at_risk  <- TRUE
      } else if (delta_i == 0 && U_i < tau) {
        # Censored before horizon — cannot determine outcome
        event_wh <- 0L
        at_risk  <- FALSE
      } else {
        event_wh <- 0L
        at_risk  <- TRUE
      }

      rows[[length(rows) + 1L]] <- data.frame(
        subject_id           = test_subjects$id[i],
        landmark             = landmark[i],
        horizon              = tau,
        method               = method_name,
        risk_pred            = clip_probability(risk_mat[i, h]),
        event_within_horizon = event_wh,
        at_risk              = at_risk,
        fold_id              = fold_id,
        replicate_id         = 1L,
        n_train_setting      = nrow(test_subjects),
        time_grid_setting    = "macs_cv",
        stringsAsFactors     = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ── 4. Cross-validation ─────────────────────────────────────────────────────
set.seed(2024)
N_FOLDS     <- 5L
HORIZONS    <- silk_opt("PREDICTION_HORIZONS")
METHOD_LIST <- c("Landmark-Recorded", "MMLM-Recorded", "SILK")

ids     <- sort(unique(subjects$id))
n       <- length(ids)
folds   <- sample(rep(seq_len(N_FOLDS), length.out = n))
names(folds) <- ids

all_predictions <- list()
all_status      <- list()
fold_counter    <- 0L

cat("\n=== Starting", N_FOLDS, "-fold cross-validation ===\n")
cat("Methods:", paste(METHOD_LIST, collapse = ", "), "\n")
cat("Horizons:", paste(HORIZONS, collapse = ", "), "years\n\n")

for (ff in seq_len(N_FOLDS)) {
  cat("── Fold", ff, "/", N_FOLDS, "──\n")

  train_ids <- ids[folds != ff]
  test_ids  <- ids[folds == ff]

  train_subjects <- subjects[subjects$id %in% train_ids, , drop = FALSE]
  train_visits   <- visits[visits$id %in% train_ids, , drop = FALSE]
  test_subjects  <- subjects[subjects$id %in% test_ids, , drop = FALSE]
  test_visits    <- visits[visits$id %in% test_ids, , drop = FALSE]

  cat("  Train:", nrow(train_subjects), "subjects,",
      nrow(train_visits), "visits\n")
  cat("  Test: ", nrow(test_subjects), "subjects,",
      nrow(test_visits), "visits\n")

  for (method in METHOD_LIST) {
    cat("  Fitting", method, "... ")
    t0 <- proc.time()[3]

    fit_result <- tryCatch({
      fit <- switch(
        method,
        "Landmark-Recorded" = {
          x <- landmark_covariates(train_subjects, train_visits,
                                   clock = "recorded", include_biomarker = FALSE)
          list(method = method, fit = fit_residual_cox(train_subjects, x))
        },
        "MMLM-Recorded" = {
          mm <- fit_current_value_mixed_model(train_subjects, train_visits)
          marker_hat <- predict_current_value_mixed_model(mm, train_subjects, train_visits)
          x <- landmark_covariates(train_subjects, train_visits,
                                   clock = "recorded", marker = marker_hat)
          list(method = method, marker_model = mm,
               fit = fit_residual_cox(train_subjects, x))
        },
        "SILK" = fit_silk(
          train_subjects, train_visits,
          scenario_name = SCENARIO_NAME,
          seed = 42L + ff
        )
      )
      list(ok = TRUE, value = fit)
    }, error = function(e) {
      list(ok = FALSE, error = conditionMessage(e))
    })

    fit_time <- proc.time()[3] - t0

    if (!fit_result$ok) {
      cat("FAILED (", fit_result$error, ")\n")
      all_status[[paste0(ff, "_", method)]] <- data.frame(
        fold = ff, method = method, fit_ok = FALSE, predict_ok = FALSE,
        fit_time = fit_time, predict_time = NA,
        error = fit_result$error, stringsAsFactors = FALSE
      )
      next
    }

    # Predict
    pred_t0 <- proc.time()[3]
    pred_result <- tryCatch({
      pred <- switch(
        method,
        "Landmark-Recorded" = {
          x <- landmark_covariates(test_subjects, test_visits,
                                   clock = "recorded", include_biomarker = FALSE)
          risk <- predict_residual_cox_risk(fit_result$value$fit, x, HORIZONS)
          real_prediction_frame(test_subjects, test_subjects$A_obs,
                                HORIZONS, risk, method, fold_id = ff)
        },
        "MMLM-Recorded" = {
          marker_hat <- predict_current_value_mixed_model(
            fit_result$value$marker_model, test_subjects, test_visits
          )
          x <- landmark_covariates(test_subjects, test_visits,
                                   clock = "recorded", marker = marker_hat)
          risk <- predict_residual_cox_risk(fit_result$value$fit, x, HORIZONS)
          real_prediction_frame(test_subjects, test_subjects$A_obs,
                                HORIZONS, risk, method, fold_id = ff)
        },
        "SILK" = {
          grid <- make_shift_grid(SCENARIO_NAME)
          pred_shift <- predict_registration_shift(
            fit_result$value$registration$final_template,
            test_visits, grid
          )
          stage <- test_subjects$A_obs[match(pred_shift$id, test_subjects$id)] -
                   pred_shift$e_hat
          stage <- stage[match(test_subjects$id, pred_shift$id)]

          test_history <- make_history_features(test_subjects, test_visits)
          x <- silk_history_covariates(test_subjects, test_history, stage)
          risk <- predict_residual_cox_risk(fit_result$value$fit$fit, x, HORIZONS)
          real_prediction_frame(test_subjects, stage,
                                HORIZONS, risk, method, fold_id = ff)
        }
      )
      list(ok = TRUE, value = pred)
    }, error = function(e) {
      list(ok = FALSE, error = conditionMessage(e))
    })

    predict_time <- proc.time()[3] - pred_t0
    cat(sprintf("%.1fs fit, %.1fs predict",
                fit_time, predict_time))

    if (pred_result$ok) {
      cat(" ✓\n")
      fold_counter <- fold_counter + 1L
      all_predictions[[fold_counter]] <- pred_result$value
    } else {
      cat(" PREDICT FAILED (", pred_result$error, ")\n")
    }

    all_status[[paste0(ff, "_", method)]] <- data.frame(
      fold = ff, method = method,
      fit_ok = TRUE, predict_ok = pred_result$ok,
      fit_time = fit_time, predict_time = predict_time,
      error = if (!pred_result$ok) pred_result$error else "",
      stringsAsFactors = FALSE
    )
  }
  cat("\n")
}

# ── 5. Combine and save ─────────────────────────────────────────────────────
predictions <- do.call(rbind, all_predictions)
rownames(predictions) <- NULL
status_df <- do.call(rbind, all_status)
rownames(status_df) <- NULL

write.csv(predictions,
          file.path(RESULTS_DIR, "macs_predictions.csv"),
          row.names = FALSE)
write.csv(status_df,
          file.path(RESULTS_DIR, "macs_method_status.csv"),
          row.names = FALSE)

cat("=== Cross-validation complete ===\n")
cat("Total predictions:", nrow(predictions), "\n")
cat("Saved to:", file.path(RESULTS_DIR, "macs_predictions.csv"), "\n")

# ── 6. Quick summary ────────────────────────────────────────────────────────
cat("\nMethod status:\n")
print(table(status_df$method, status_df$predict_ok))

cat("\nPredictions per method:\n")
print(table(predictions$method))

cat("\nAt-risk counts per method/horizon:\n")
print(xtabs(at_risk ~ method + horizon, data = predictions))
