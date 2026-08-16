# =============================================================================
# methods_survival_benchmarks.R
# Simulation/application-only recorded-clock and oracle benchmarks.
# These are deliberately outside the SILK package.
# =============================================================================

fit_recorded_cox <- function(train_subjects) {
  x <- age_only_covariates(train_subjects)
  list(
    method = "Recorded-Cox",
    fit = fit_age_scale_cox(train_subjects, train_subjects$A_obs, x),
    survival_time_scale = "attained_age",
    implementation = "attained-age Cox using recorded age only"
  )
}

predict_recorded_cox <- function(fit, test_subjects, horizons = NULL,
                                 fold_id = 1L, replicate_id = 1L,
                                 n_train_setting = NA,
                                 time_grid_setting = NA) {
  if (is.null(horizons)) horizons <- PREDICTION_HORIZONS
  x <- age_only_covariates(test_subjects)
  risk <- predict_age_scale_cox_risk(fit$fit, test_subjects$A_obs, x, horizons)
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

fit_recorded_beran <- function(train_subjects) {
  state <- beran_state_covariates(
    train_subjects, train_subjects$A_obs, include_x = FALSE
  )
  list(
    method = "Recorded-Beran",
    fit = fit_beran_risk(train_subjects, state),
    implementation = "Beran estimator using recorded age only"
  )
}

predict_recorded_beran <- function(fit, test_subjects,
                                   horizons = NULL,
                                   fold_id = 1L, replicate_id = 1L,
                                   n_train_setting = NA,
                                   time_grid_setting = NA) {
  if (is.null(horizons)) horizons <- PREDICTION_HORIZONS
  state <- beran_state_covariates(
    test_subjects, test_subjects$A_obs, include_x = FALSE
  )
  risk <- predict_beran_risk(fit$fit, state, horizons)
  prediction_frame(
    test_subjects,
    landmark = test_subjects$A_obs,
    horizons = horizons,
    risk_mat = risk,
    context = method_context(
      fit$method, fold_id, replicate_id, n_train_setting, time_grid_setting
    )
  )
}

fit_oracle_cox <- function(train_subjects) {
  x <- age_only_covariates(train_subjects)
  list(
    method = "Oracle-Cox",
    fit = fit_age_scale_cox(train_subjects, train_subjects$A_star, x),
    survival_time_scale = "attained_age",
    implementation = "attained-age Cox using latent age only"
  )
}

fit_oracle_beran <- function(train_subjects) {
  state <- beran_state_covariates(
    train_subjects, train_subjects$A_star, include_x = FALSE
  )
  list(
    method = "Oracle-Beran",
    fit = fit_beran_risk(train_subjects, state),
    implementation = "Beran estimator using latent age only"
  )
}

predict_oracle_cox <- function(fit, test_subjects, horizons = NULL,
                               fold_id = 1L, replicate_id = 1L,
                               n_train_setting = NA,
                               time_grid_setting = NA) {
  if (is.null(horizons)) horizons <- PREDICTION_HORIZONS
  x <- age_only_covariates(test_subjects)
  risk <- predict_age_scale_cox_risk(fit$fit, test_subjects$A_star, x, horizons)
  prediction_frame(
    test_subjects,
    landmark = test_subjects$A_star,
    horizons = horizons,
    risk_mat = risk,
    context = method_context(
      fit$method, fold_id, replicate_id,
      n_train_setting, time_grid_setting
    )
  )
}

predict_oracle_beran <- function(fit, test_subjects,
                                 horizons = NULL,
                                 fold_id = 1L, replicate_id = 1L,
                                 n_train_setting = NA,
                                 time_grid_setting = NA) {
  if (is.null(horizons)) horizons <- PREDICTION_HORIZONS
  state <- beran_state_covariates(
    test_subjects, test_subjects$A_star, include_x = FALSE
  )
  risk <- predict_beran_risk(fit$fit, state, horizons)
  prediction_frame(
    test_subjects,
    landmark = test_subjects$A_star,
    horizons = horizons,
    risk_mat = risk,
    context = method_context(
      fit$method, fold_id, replicate_id, n_train_setting, time_grid_setting
    )
  )
}
