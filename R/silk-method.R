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
#' @param feature_type Registration feature map. \code{"silk"} uses the
#'   distributional SILK feature map; \code{"mean"} uses mean/path features.
#' @param method Character label stored in prediction outputs.
#' @param kernel Registration kernel. One of \code{"rbf"}, \code{"matern"},
#'   \code{"polynomial"}, or \code{"linear"}. Defaults to
#'   \code{silk_opt("REGISTRATION_KERNEL")}.
#' @param kernel_approx Kernel approximation. Use \code{"exact"} (or
#'   \code{"none"}) for the full kernel calculation, or \code{"rff"} for
#'   random Fourier features. RFF is available for RBF and Matern kernels.
#' @param rff_dim Number of random Fourier features when
#'   \code{kernel_approx = "rff"}.
#' @param rff_seed Integer or NULL. Random seed for RFF basis generation.
#' @return A fitted SILK model object (list) for use with \code{predict_silk}.
#' @export
#' @examples
#' \dontrun{
#' dat <- generate_dataset_fixed(200, "mean_moderate", seed = 1)
#' fit <- fit_silk(dat$subjects, dat$visits, shift_range = c(-12, 12), seed = 1)
#' }
fit_silk <- function(train_subjects, train_visits, shift_grid = NULL, seed = NULL,
                     shift_range = NULL, feature_type = c("silk", "mean"),
                     method = "SILK",
                     kernel = NULL, kernel_approx = NULL,
                     rff_dim = NULL, rff_seed = seed) {
  grid <- resolve_silk_shift_grid(shift_grid, shift_range)
  feature_type <- match.arg(feature_type)
  kernel_config <- registration_kernel_config(
    kernel = kernel,
    approximation = kernel_approx,
    rff_dim = rff_dim,
    rff_seed = rff_seed
  )
  cf <- crossfit_registration(
    train_visits, train_subjects,
    feature_type = feature_type,
    grid = grid,
    seed = seed,
    kernel_config = kernel_config
  )
  train_history <- make_history_features(train_subjects, train_visits)
  train_stage <- cf$train_stage$S_hat[match(train_subjects$id, cf$train_stage$id)]
  x <- silk_history_covariates(train_subjects, train_history, train_stage)
  risk_fit <- fit_residual_cox(train_subjects, x)
  list(
    method = method,
    grid = grid,
    shift_grid = grid,
    shift_range = range(grid),
    feature_type = feature_type,
    kernel = kernel_config,
    registration = cf,
    train_history = train_history,
    train_stage = train_stage,
    fit = risk_fit
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

#' @keywords internal
fit_beran_recorded <- function(train_subjects) {
  state <- beran_state_covariates(train_subjects, train_subjects$A_obs)
  list(method = "Beran-Recorded", fit = fit_beran_risk(train_subjects, state))
}

#' @keywords internal
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

#' @keywords internal
fit_beran_silk <- function(train_subjects, train_visits, shift_grid = NULL,
                           shift_range = NULL, seed = NULL) {
  grid <- resolve_silk_shift_grid(shift_grid, shift_range)
  cf <- crossfit_registration(
    train_visits, train_subjects,
    feature_type = "silk",
    grid = grid,
    seed = seed,
    kernel_config = registration_kernel_config()
  )
  train_stage <- cf$train_stage$S_hat[match(train_subjects$id, cf$train_stage$id)]
  state <- beran_state_covariates(train_subjects, train_stage)
  list(
    method = "Beran-SILK",
    grid = grid,
    registration = cf,
    fit = fit_beran_risk(train_subjects, state)
  )
}

#' @keywords internal
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

#' @keywords internal
fit_oracle_latent_age <- function(train_subjects, train_visits) {
  x <- landmark_covariates(train_subjects, train_visits, clock = "latent")
  fit <- fit_residual_cox(train_subjects, x)
  list(method = "Oracle-Latent-Age", fit = fit)
}

#' @keywords internal
fit_beran_oracle_latent_age <- function(train_subjects) {
  state <- beran_state_covariates(train_subjects, train_subjects$A_star)
  list(method = "Beran-Oracle-Latent-Age", fit = fit_beran_risk(train_subjects, state))
}

#' @keywords internal
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

#' @keywords internal
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
