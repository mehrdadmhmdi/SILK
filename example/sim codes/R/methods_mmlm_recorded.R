# =============================================================================
# methods_mmlm_recorded.R
# MMLM: correct and deliberately misspecified mixed-model landmarking.
# A mixed model smooths the current biomarker value only; the landmark Cox model
# uses baseline covariates, recorded landmark age, and that smoothed value.
# Reference: Rizopoulos et al. dynamic prediction work; nlme::lme and survival.
# =============================================================================

fit_mmlm_recorded <- function(train_subjects, train_visits,
                              specification = c("correct", "misspecified")) {
  specification <- match.arg(specification)
  method <- if (specification == "correct") "MMLM-Correct" else "MMLM-Misspecified"
  marker_model <- fit_current_value_mixed_model(
    train_subjects, train_visits, specification = specification
  )
  marker_hat <- predict_current_value_mixed_model(marker_model, train_subjects, train_visits)
  x <- landmark_covariates(
    train_subjects,
    train_visits,
    clock = "recorded",
    marker = marker_hat,
    include_biomarker = TRUE,
    covariate_specification = specification
  )
  fit <- fit_residual_cox(train_subjects, x)
  list(
    method = method,
    specification = specification,
    marker_model = marker_model,
    fit = fit,
    implementation = paste0(
      "Mixed-model landmarking with ", specification,
      " longitudinal and event covariate specifications"
    )
  )
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
    include_biomarker = TRUE,
    covariate_specification = fit$specification
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
