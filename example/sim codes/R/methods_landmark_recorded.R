# =============================================================================
# methods_landmark_recorded.R
# Landmark-Recorded: landmark Cox model on recorded time without biomarkers.
# Reference: van Houwelingen (2006) landmarking; survival::coxph for the Cox fit.
# =============================================================================

fit_landmark_recorded <- function(train_subjects, train_visits) {
  x <- landmark_covariates(
    train_subjects,
    train_visits,
    clock = "recorded",
    include_biomarker = FALSE
  )
  fit <- fit_residual_cox(train_subjects, x)
  list(method = "Landmark-Recorded", fit = fit)
}

predict_landmark_recorded <- function(fit, test_subjects, test_visits,
                                      horizons = PREDICTION_HORIZONS,
                                      fold_id = 1L, replicate_id = 1L,
                                      n_train_setting = NA,
                                      time_grid_setting = NA) {
  x <- landmark_covariates(
    test_subjects,
    test_visits,
    clock = "recorded",
    include_biomarker = FALSE
  )
  risk <- predict_residual_cox_risk(fit$fit, x, horizons)
  prediction_frame(
    test_subjects,
    landmark = test_subjects$A_obs,
    horizons = horizons,
    risk_mat = risk,
    context = method_context(
      fit$method, fold_id, replicate_id,
      n_train_setting, time_grid_setting
    )
  )
}
