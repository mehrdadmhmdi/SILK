# =============================================================================
# methods_mmlm_recorded.R
# MMLM-Recorded: mixed-model landmarking on recorded time.
# A mixed model smooths the current biomarker value only; the landmark Cox model
# uses baseline covariates, recorded landmark age, and that smoothed value.
# Reference: Rizopoulos et al. dynamic prediction work; nlme::lme and survival.
# =============================================================================

fit_mmlm_recorded <- function(train_subjects, train_visits) {
  marker_model <- fit_current_value_mixed_model(train_subjects, train_visits)
  marker_hat <- predict_current_value_mixed_model(marker_model, train_subjects, train_visits)
  x <- landmark_covariates(
    train_subjects,
    train_visits,
    clock = "recorded",
    marker = marker_hat,
    include_biomarker = TRUE
  )
  fit <- fit_residual_cox(train_subjects, x)
  list(method = "MMLM-Recorded", marker_model = marker_model, fit = fit)
}

predict_mmlm_recorded <- function(fit, test_subjects, test_visits,
                                  horizons = PREDICTION_HORIZONS,
                                  fold_id = 1L, replicate_id = 1L,
                                  n_train_setting = NA,
                                  time_grid_setting = NA) {
  marker_hat <- predict_current_value_mixed_model(fit$marker_model, test_subjects, test_visits)
  x <- landmark_covariates(
    test_subjects,
    test_visits,
    clock = "recorded",
    marker = marker_hat,
    include_biomarker = TRUE
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
