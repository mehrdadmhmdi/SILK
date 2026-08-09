# =============================================================================
# silk-method.R
# fit_silk, predict_silk, and internal oracle (from methods_silk.R)
# =============================================================================

#' @keywords internal
resolve_silk_shift_grid <- function(shift_grid = NULL, shift_range = NULL) {
  if (!is.null(shift_grid)) {
    if (is.character(shift_grid) && length(shift_grid) == 1L) {
      warning(
        "Passing a simulation scenario name to fit_silk() is deprecated. ",
        "Use shift_range or shift_grid for real data, or make_shift_grid() ",
        "explicitly in simulation code.",
        call. = FALSE
      )
      return(make_shift_grid(shift_grid))
    }
    grid <- sort(unique(as.numeric(shift_grid)))
    grid <- grid[is.finite(grid)]
    if (length(grid) < 2L) {
      stop("shift_grid must contain at least two finite numeric values.", call. = FALSE)
    }
    return(grid)
  }

  if (is.null(shift_range)) {
    shift_range <- c(
      silk_opt("DEFAULT_SHIFT_GRID_MIN"),
      silk_opt("DEFAULT_SHIFT_GRID_MAX")
    )
  }
  shift_range <- as.numeric(shift_range)
  if (length(shift_range) != 2L || any(!is.finite(shift_range))) {
    stop("shift_range must be a finite numeric vector of length two.", call. = FALSE)
  }
  step <- as.numeric(silk_opt("SHIFT_GRID_STEP"))
  if (!is.finite(step) || step <= 0) {
    stop("SHIFT_GRID_STEP must be a positive finite number.", call. = FALSE)
  }
  rng <- sort(shift_range)
  seq(rng[1], rng[2], by = step)
}

#' Fit the SILK registration layer
#'
#' Estimates subject-specific origin shifts by cross-fitted registration. The
#' biomarker discrepancy is the exact RKHS loss induced by a Gaussian RBF
#' kernel on standardized biomarker vectors. This kernel is characteristic;
#' the implementation does not replace it with biomarker moments or a finite
#' biomarker feature map.
#'
#' @param train_subjects Data frame of training subjects.
#' @param train_visits Data frame of training visits with biomarker columns.
#' @param shift_grid Numeric vector of candidate origin shifts.
#' @param shift_range Numeric vector of length two used to construct the grid
#'   when \code{shift_grid} is omitted.
#' @param seed Integer or NULL. Random seed for folds and multistart fitting.
#' @return A \code{silk_registration} object containing cross-fitted training
#'   shifts and the final template used for new subjects.
#' @export
fit_silk_registration <- function(train_subjects, train_visits,
                                  shift_grid = NULL, shift_range = NULL,
                                  seed = NULL) {
  grid <- resolve_silk_shift_grid(shift_grid, shift_range)
  crossfit_registration(train_visits, train_subjects, grid = grid, seed = seed)
}

#' Predict origin shifts from a SILK registration
#'
#' @param registration Fitted object from \code{fit_silk_registration}.
#' @param new_visits Visit data for new subjects.
#' @param shift_grid Optional candidate grid. By default the fitted grid is
#'   reused.
#' @return Data frame with estimated shifts and profile-loss diagnostics.
#' @export
predict_silk_registration <- function(registration, new_visits, shift_grid = NULL) {
  if (!inherits(registration, "silk_registration") ||
      is.null(registration$final_template)) {
    stop("registration must be a fitted silk_registration object.", call. = FALSE)
  }
  grid <- if (is.null(shift_grid)) registration$grid else {
    value <- sort(unique(as.numeric(shift_grid)))
    value <- value[is.finite(value)]
    if (length(value) < 2L) stop("shift_grid must contain at least two values.", call. = FALSE)
    value
  }
  predict_registration_shift(registration$final_template, new_visits, grid)
}

#' @keywords internal
validate_silk_registration <- function(registration, train_subjects) {
  if (!inherits(registration, "silk_registration") ||
      is.null(registration$train_stage) || is.null(registration$final_template)) {
    stop("registration must be a fitted silk_registration object.", call. = FALSE)
  }
  missing_ids <- setdiff(train_subjects$id, registration$train_stage$id)
  if (length(missing_ids)) {
    stop("registration is missing training subject ids.", call. = FALSE)
  }
  registration
}

#' Fit a SILK model
#'
#' Fits the SILK (Shift-Invariant Learned Kernel) model for absolute risk
#' prediction under origin-time measurement error. For real data, the error
#' distribution is not assumed known; the user supplies only a candidate shift
#' grid or range over which registration is optimized.
#'
#' @param train_subjects Data frame of training subjects.
#' @param train_visits Data frame of training visits with biomarker columns.
#' @param shift_grid Numeric vector of candidate origin shifts. If omitted,
#'   the grid is built from \code{shift_range} and \code{SHIFT_GRID_STEP}.
#' @param seed Integer or NULL. Random seed for reproducibility.
#' @param shift_range Numeric vector of length two giving the lower and upper
#'   candidate shift values. Defaults to \code{DEFAULT_SHIFT_GRID_MIN} and
#'   \code{DEFAULT_SHIFT_GRID_MAX}.
#' @param method Character label stored in prediction outputs.
#' @param registration Optional fitted object from
#'   \code{fit_silk_registration}. Supplying it lets multiple survival layers
#'   reuse exactly the same cross-fitted registration.
#' @return A fitted SILK model object (list) for use with \code{predict_silk}.
#' @export
#' @examples
#' \dontrun{
#' dat <- generate_dataset_fixed(200, "mean_moderate", seed = 1)
#' fit <- fit_silk(dat$subjects, dat$visits, shift_range = c(-12, 12), seed = 1)
#' }
fit_silk <- function(train_subjects, train_visits, shift_grid = NULL, seed = NULL,
                     shift_range = NULL, method = "SILK", registration = NULL) {
  if (is.null(registration)) {
    registration <- fit_silk_registration(
      train_subjects, train_visits,
      shift_grid = shift_grid, shift_range = shift_range, seed = seed
    )
  } else {
    registration <- validate_silk_registration(registration, train_subjects)
    if (!is.null(shift_grid) || !is.null(shift_range)) {
      requested_grid <- resolve_silk_shift_grid(shift_grid, shift_range)
      if (!isTRUE(all.equal(requested_grid, registration$grid))) {
        stop("The supplied registration and requested shift grid differ.", call. = FALSE)
      }
    }
  }
  grid <- registration$grid
  train_history <- make_history_features(train_subjects, train_visits)
  train_stage <- registration$train_stage$S_hat[
    match(train_subjects$id, registration$train_stage$id)
  ]
  x <- silk_history_covariates(train_subjects, train_history, train_stage)
  risk_fit <- fit_residual_cox(train_subjects, x)
  list(
    method = method,
    grid = grid,
    shift_grid = grid,
    shift_range = range(grid),
    biomarker_kernel = list(
      name = "gaussian_rbf",
      characteristic = TRUE,
      bandwidth = registration$final_template$biomarker_builder$bandwidth,
      bandwidth_rule = registration$final_template$biomarker_builder$bandwidth_rule
    ),
    registration = registration,
    train_history = train_history,
    train_stage = train_stage,
    fit = risk_fit,
    implementation = registration$implementation
  )
}

#' Fit the same-feature recorded-age Cox comparator
#'
#' Uses the same biomarker-history feature vector as SILK but keeps recorded age.
#'
#' @param train_subjects Data frame of training subjects.
#' @param train_visits Data frame of training visits.
#' @return Fitted comparator object.
#' @export
fit_same_feature_recorded_cox <- function(train_subjects, train_visits) {
  x <- recorded_same_feature_covariates(train_subjects, train_visits)
  list(
    method = "Cox-SameFeature-Recorded",
    fit = fit_residual_cox(train_subjects, x)
  )
}

#' Predict from the same-feature recorded-age Cox comparator
#'
#' @param fit Fitted object from \code{fit_same_feature_recorded_cox}.
#' @param test_subjects Data frame of test subjects.
#' @param test_visits Data frame of test visits.
#' @param horizons Numeric vector of prediction horizons.
#' @param fold_id Fold identifier.
#' @param replicate_id Replicate identifier.
#' @param n_train_setting Training-size label.
#' @param time_grid_setting Time-grid label.
#' @return A prediction frame.
#' @export
predict_same_feature_recorded_cox <- function(fit, test_subjects, test_visits,
                                              horizons = NULL,
                                              fold_id = 1L, replicate_id = 1L,
                                              n_train_setting = NA,
                                              time_grid_setting = NA) {
  if (is.null(horizons)) horizons <- silk_opt("PREDICTION_HORIZONS")
  x <- recorded_same_feature_covariates(test_subjects, test_visits)
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

#' Beran and latent-age oracle survival layers
#'
#' Convenience fits and prediction methods used by the simulation study. The
#' oracle methods require latent ages and are therefore simulation-only.
#'
#' @param train_subjects Training subject data frame.
#' @param train_visits Training visit data frame.
#' @param test_subjects Test subject data frame.
#' @param test_visits Test visit data frame.
#' @param fit A fitted object from the corresponding fit function.
#' @param horizons Positive prediction horizons.
#' @param shift_grid Candidate registration-shift grid.
#' @param shift_range Two-element candidate shift range.
#' @param seed Integer or NULL.
#' @param registration Optional fitted object from
#'   \code{fit_silk_registration}; supplying it lets Cox and Beran survival
#'   layers share one registration fit.
#' @param fold_id Fold identifier for bookkeeping.
#' @param replicate_id Replicate identifier for bookkeeping.
#' @param n_train_setting Training-size label for bookkeeping.
#' @param time_grid_setting Time-grid label for bookkeeping.
#' @return A fitted survival-layer object for fit functions, or a standardized
#'   prediction frame for prediction functions.
#' @name silk_survival_layers
#' @export
fit_beran_recorded <- function(train_subjects) {
  state <- beran_state_covariates(train_subjects, train_subjects$A_obs)
  list(method = "Beran-Recorded", fit = fit_beran_risk(train_subjects, state))
}

#' @rdname silk_survival_layers
#' @export
predict_beran_recorded <- function(fit, test_subjects,
                                   horizons = NULL,
                                   fold_id = 1L, replicate_id = 1L,
                                   n_train_setting = NA,
                                   time_grid_setting = NA) {
  if (is.null(horizons)) horizons <- silk_opt("PREDICTION_HORIZONS")
  state <- beran_state_covariates(test_subjects, test_subjects$A_obs)
  risk <- predict_beran_risk(fit$fit, state, horizons)
  prediction_frame(
    test_subjects,
    landmark = test_subjects$A_obs,
    horizons = horizons,
    risk_mat = risk,
    context = method_context(fit$method, fold_id, replicate_id, n_train_setting, time_grid_setting)
  )
}

#' @rdname silk_survival_layers
#' @export
fit_beran_silk <- function(train_subjects, train_visits, shift_grid = NULL,
                           shift_range = NULL, seed = NULL, registration = NULL) {
  if (is.null(registration)) {
    registration <- fit_silk_registration(
      train_subjects, train_visits,
      shift_grid = shift_grid, shift_range = shift_range, seed = seed
    )
  } else {
    registration <- validate_silk_registration(registration, train_subjects)
    if (!is.null(shift_grid) || !is.null(shift_range)) {
      requested_grid <- resolve_silk_shift_grid(shift_grid, shift_range)
      if (!isTRUE(all.equal(requested_grid, registration$grid))) {
        stop("The supplied registration and requested shift grid differ.", call. = FALSE)
      }
    }
  }
  grid <- registration$grid
  train_stage <- registration$train_stage$S_hat[
    match(train_subjects$id, registration$train_stage$id)
  ]
  state <- beran_state_covariates(train_subjects, train_stage)
  list(
    method = "Beran-SILK",
    grid = grid,
    registration = registration,
    fit = fit_beran_risk(train_subjects, state),
    implementation = registration$implementation
  )
}

#' @rdname silk_survival_layers
#' @export
predict_beran_silk <- function(fit, test_subjects, test_visits,
                               horizons = NULL,
                               fold_id = 1L, replicate_id = 1L,
                               n_train_setting = NA,
                               time_grid_setting = NA) {
  if (is.null(horizons)) horizons <- silk_opt("PREDICTION_HORIZONS")
  ps <- predict_registration_shift(fit$registration$final_template, test_visits, fit$grid)
  stage <- test_subjects$A_obs[match(ps$id, test_subjects$id)] - ps$e_hat
  stage <- stage[match(test_subjects$id, ps$id)]
  state <- beran_state_covariates(test_subjects, stage)
  risk <- predict_beran_risk(fit$fit, state, horizons)
  prediction_frame(
    test_subjects,
    landmark = stage,
    horizons = horizons,
    risk_mat = risk,
    context = method_context(fit$method, fold_id, replicate_id, n_train_setting, time_grid_setting)
  )
}

#' Predict risk from a SILK model
#'
#' Generates absolute risk predictions for new subjects using a fitted SILK model.
#'
#' @param fit Fitted SILK model from \code{fit_silk}.
#' @param test_subjects Data frame of test subjects.
#' @param test_visits Data frame of test visits with biomarker columns.
#' @param horizons Numeric vector of prediction horizons.
#' @param fold_id Integer. Fold identifier for bookkeeping.
#' @param replicate_id Integer. Replicate identifier for bookkeeping.
#' @param n_train_setting Numeric. Training set size setting for bookkeeping.
#' @param time_grid_setting Character. Time grid setting for bookkeeping.
#' @return A prediction frame data frame.
#' @export
#' @examples
#' \dontrun{
#' dat <- generate_dataset_fixed(200, "mean_moderate", seed = 1)
#' fit <- fit_silk(dat$subjects, dat$visits, shift_range = c(-12, 12), seed = 1)
#' test <- generate_dataset_fixed(100, "mean_moderate", seed = 2)
#' pred <- predict_silk(fit, test$subjects, test$visits)
#' }
predict_silk <- function(fit, test_subjects, test_visits,
                         horizons = NULL,
                         fold_id = 1L, replicate_id = 1L,
                         n_train_setting = NA,
                         time_grid_setting = NA) {
  if (is.null(horizons)) horizons <- silk_opt("PREDICTION_HORIZONS")
  pred_shift <- predict_registration_shift(fit$registration$final_template, test_visits, fit$grid)
  stage <- test_subjects$A_obs[match(pred_shift$id, test_subjects$id)] - pred_shift$e_hat
  stage <- stage[match(test_subjects$id, pred_shift$id)]
  test_history <- make_history_features(test_subjects, test_visits)
  x <- silk_history_covariates(test_subjects, test_history, stage)
  risk <- predict_residual_cox_risk(fit$fit, x, horizons)
  prediction_frame(
    test_subjects,
    landmark = stage,
    horizons = horizons,
    risk_mat = risk,
    context = method_context(
      fit$method, fold_id, replicate_id,
      n_train_setting, time_grid_setting
    )
  )
}

#' @rdname silk_survival_layers
#' @export
fit_oracle_latent_age <- function(train_subjects, train_visits) {
  x <- landmark_covariates(train_subjects, train_visits, clock = "latent")
  fit <- fit_residual_cox(train_subjects, x)
  list(method = "Oracle-Latent-Age", fit = fit)
}

#' @rdname silk_survival_layers
#' @export
fit_beran_oracle_latent_age <- function(train_subjects) {
  state <- beran_state_covariates(train_subjects, train_subjects$A_star)
  list(method = "Beran-Oracle-Latent-Age", fit = fit_beran_risk(train_subjects, state))
}

#' @rdname silk_survival_layers
#' @export
predict_oracle_latent_age <- function(fit, test_subjects, test_visits,
                                      horizons = NULL,
                                      fold_id = 1L, replicate_id = 1L,
                                      n_train_setting = NA,
                                      time_grid_setting = NA) {
  if (is.null(horizons)) horizons <- silk_opt("PREDICTION_HORIZONS")
  x <- landmark_covariates(test_subjects, test_visits, clock = "latent")
  risk <- predict_residual_cox_risk(fit$fit, x, horizons)
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

#' @rdname silk_survival_layers
#' @export
predict_beran_oracle_latent_age <- function(fit, test_subjects,
                                            horizons = NULL,
                                            fold_id = 1L, replicate_id = 1L,
                                            n_train_setting = NA,
                                            time_grid_setting = NA) {
  if (is.null(horizons)) horizons <- silk_opt("PREDICTION_HORIZONS")
  state <- beran_state_covariates(test_subjects, test_subjects$A_star)
  risk <- predict_beran_risk(fit$fit, state, horizons)
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
