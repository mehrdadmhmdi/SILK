# =============================================================================
# methods_timeerror_integrated.R
# Parametric recorded-time error as nuisance uncertainty, handled by Monte Carlo
# integration over shifted recorded ages. This is a sensitivity landmark model
# under a working error distribution, not a fully Bayesian posterior error model.
# =============================================================================

fit_timeerror_integrated_landmark <- function(train_subjects, train_visits,
                                              n_impute = TIMEERROR_N_IMPUTE,
                                              error_sd = TIMEERROR_SD,
                                              use_truth = TIMEERROR_USE_TRUTH,
                                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_impute <- max(1L, as.integer(n_impute))
  draws <- draw_time_error(train_subjects, n_impute, sd = error_sd, use_truth = use_truth)
  models <- vector("list", n_impute)
  for (b in seq_len(n_impute)) {
    shifted <- shift_recorded_time(train_subjects, train_visits, draws[, b])
    x <- landmark_covariates(
      shifted$subjects,
      shifted$visits,
      clock = "recorded",
      include_biomarker = FALSE
    )
    models[[b]] <- fit_residual_cox(shifted$subjects, x)
  }
  list(
    method = "TimeError-Integrated-Landmark",
    models = models,
    n_impute = n_impute,
    error_sd = error_sd,
    use_truth = use_truth,
    seed = seed
  )
}

predict_timeerror_integrated_landmark <- function(fit, test_subjects, test_visits,
                                                  horizons = PREDICTION_HORIZONS,
                                                  fold_id = 1L, replicate_id = 1L,
                                                  n_train_setting = NA,
                                                  time_grid_setting = NA) {
  if (!is.null(fit$seed)) set.seed(fit$seed + 20011L)
  draws <- draw_time_error(test_subjects, fit$n_impute, sd = fit$error_sd, use_truth = fit$use_truth)
  risk_sum <- matrix(0, nrow = nrow(test_subjects), ncol = length(horizons))
  landmark_sum <- numeric(nrow(test_subjects))
  for (b in seq_len(fit$n_impute)) {
    shifted <- shift_recorded_time(test_subjects, test_visits, draws[, b])
    x <- landmark_covariates(
      shifted$subjects,
      shifted$visits,
      clock = "recorded",
      include_biomarker = FALSE
    )
    risk_sum <- risk_sum + predict_residual_cox_risk(fit$models[[b]], x, horizons)
    landmark_sum <- landmark_sum + shifted$subjects$A_obs
  }
  prediction_frame(
    test_subjects,
    landmark = landmark_sum / fit$n_impute,
    horizons = horizons,
    risk_mat = risk_sum / fit$n_impute,
    context = method_context(
      fit$method, fold_id, replicate_id,
      n_train_setting, time_grid_setting
    )
  )
}
