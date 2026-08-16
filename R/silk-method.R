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
#' biomarker discrepancy is the exact RKHS loss induced by a supported
#' characteristic kernel on standardized biomarker vectors. Gaussian RBF,
#' Laplace, and Matern-3/2 kernels are available; none is replaced by biomarker
#' moments or a finite biomarker feature map.
#'
#' @param train_subjects Data frame of training subjects.
#' @param train_visits Data frame of training visits with biomarker columns.
#' @param shift_grid Numeric vector of candidate origin shifts.
#' @param shift_range Numeric vector of length two used to construct the grid
#'   when \code{shift_grid} is omitted.
#' @param seed Integer or NULL. Random seed for folds and multistart fitting.
#' @param biomarker_kernel Character kernel name: \code{"gaussian"},
#'   \code{"laplace"}, or \code{"matern32"}. NULL uses the package default.
#' @param data_spec Optional object from \code{silk_data_spec()} describing
#'   user column names, biomarker/covariate columns, engineered covariates,
#'   template inputs, and the shift anchor.
#' @return A \code{silk_registration} object containing cross-fitted training
#'   shifts and the final template used for new subjects.
#' @export
fit_silk_registration <- function(train_subjects, train_visits,
                                  shift_grid = NULL, shift_range = NULL,
                                  seed = NULL, biomarker_kernel = NULL,
                                  data_spec = NULL) {
  inputs <- prepare_silk_inputs(train_subjects, train_visits, data_spec)
  train_subjects <- inputs$subjects
  train_visits <- inputs$visits
  data_spec <- inputs$spec
  grid <- resolve_silk_shift_grid(shift_grid, shift_range)
  crossfit_registration(
    train_visits, train_subjects, grid = grid, seed = seed,
    biomarker_kernel = biomarker_kernel, input_spec = data_spec
  )
}

#' Predict origin shifts from a SILK registration
#'
#' @param registration Fitted object from \code{fit_silk_registration}.
#' @param new_visits Visit data for new subjects.
#' @param shift_grid Optional candidate grid. By default the fitted grid is
#'   reused.
#' @param data_spec Optional data specification. By default the specification
#'   stored in the fitted registration is reused.
#' @param new_subjects Optional subject table for new data. Supplying it is
#'   recommended when visit_lag is derived from the subject-level recorded age.
#' @return Data frame with estimated shifts and profile-loss diagnostics.
#' @export
predict_silk_registration <- function(registration, new_visits, shift_grid = NULL,
                                      data_spec = NULL, new_subjects = NULL) {
  if (!inherits(registration, "silk_registration") ||
      is.null(registration$final_template)) {
    stop("registration must be a fitted silk_registration object.", call. = FALSE)
  }
  data_spec <- if (is.null(data_spec)) registration$input_spec else data_spec
  if (!is.null(data_spec)) {
    # Registration prediction needs only the visit table.  A supplied visit_lag
    # is preferred; otherwise repeated recorded_age values in new_visits allow
    # the adapter to derive it without a subject table.
    if (!is.null(new_subjects)) {
      prepared <- prepare_silk_inputs(new_subjects, new_visits, data_spec)
    } else if (is.null(data_spec$visit_lag) || !data_spec$visit_lag %in% names(new_visits)) {
      required_age <- data_spec$recorded_age
      if (!required_age %in% names(new_visits)) {
        stop("Supply new_subjects, or include visit_lag/repeated recorded_age in new_visits.", call. = FALSE)
      }
      pseudo_subjects <- unique(new_visits[, data_spec$subject_id, drop = FALSE])
      names(pseudo_subjects) <- data_spec$subject_id
      pseudo_subjects$.__recorded_age <- new_visits[[data_spec$recorded_age]][
        match(pseudo_subjects[[data_spec$subject_id]], new_visits[[data_spec$subject_id]])
      ]
      data_spec_prediction <- data_spec
      data_spec_prediction$recorded_age <- ".__recorded_age"
      new_visits$.__recorded_age <- new_visits[[data_spec$recorded_age]]
      prepared <- prepare_silk_inputs(pseudo_subjects, new_visits, data_spec_prediction)
    } else {
      pseudo_subjects <- data.frame(
        setNames(list(unique(new_visits[[data_spec$subject_id]])), data_spec$subject_id),
        check.names = FALSE
      )
      pseudo_subjects[[data_spec$recorded_age]] <- 0
      prepared <- prepare_silk_inputs(pseudo_subjects, new_visits, data_spec)
    }
    new_visits <- prepared$visits
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
#' @param method Character label stored in prediction outputs. Defaults to
#'   \code{"SILK-Gaussian"}.
#' @param registration Optional fitted object from
#'   \code{fit_silk_registration}. Supplying it lets multiple survival layers
#'   reuse exactly the same cross-fitted registration.
#' @param biomarker_kernel Optional characteristic biomarker kernel used when
#'   \code{registration} is not supplied.
#' @param data_spec Optional object from \code{silk_data_spec()} for arbitrary
#'   user column names and engineered covariates.
#' @return A fitted SILK model object (list) for use with \code{predict_silk}.
#' @export
#' @examples
#' \dontrun{
#' dat <- generate_dataset_fixed(200, "mean_moderate", seed = 1)
#' fit <- fit_silk(dat$subjects, dat$visits, shift_range = c(-12, 12), seed = 1)
#' }
fit_silk <- function(train_subjects, train_visits, shift_grid = NULL, seed = NULL,
                     shift_range = NULL, method = "SILK-Gaussian", registration = NULL,
                     biomarker_kernel = NULL, data_spec = NULL) {
  if (is.null(data_spec) && !is.null(registration) &&
      !is.null(registration$input_spec)) {
    data_spec <- registration$input_spec
  }
  inputs <- prepare_silk_inputs(train_subjects, train_visits, data_spec)
  train_subjects <- inputs$subjects
  train_visits <- inputs$visits
  data_spec <- inputs$spec
  if (is.null(registration)) {
    registration <- fit_silk_registration(
      train_subjects, train_visits,
      shift_grid = shift_grid, shift_range = shift_range, seed = seed,
      biomarker_kernel = biomarker_kernel, data_spec = data_spec
    )
  } else {
    registration <- validate_silk_registration(registration, train_subjects)
    if (!is.null(biomarker_kernel) && !identical(
      normalize_biomarker_kernel(biomarker_kernel), registration$biomarker_kernel
    )) {
      stop("The supplied registration and requested biomarker kernel differ.", call. = FALSE)
    }
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
  x <- base_covariates(
    train_subjects,
    covariate_cols = if (!is.null(data_spec)) data_spec$covariate_cols else NULL
  )
  risk_fit <- fit_age_scale_cox(train_subjects, train_stage, x)
  builder <- registration$final_template$biomarker_builder
  list(
    method = method,
    grid = grid,
    shift_grid = grid,
    shift_range = range(grid),
    biomarker_kernel = list(
      name = builder$kernel,
      characteristic = TRUE,
      bandwidth = builder$bandwidth,
      bandwidth_rule = builder$bandwidth_rule
    ),
    registration = registration,
    train_stage = train_stage,
    fit = risk_fit,
    input_spec = data_spec,
    survival_time_scale = "attained_age",
    implementation = registration$implementation
  )
}

# Comparator benchmarks are implemented in example/sim codes/R.

# Recorded-age prediction is implemented in example/sim codes/R.

# Recorded-age Beran benchmark and prediction are implemented in
# example/sim codes/R.

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
  data_spec <- fit$input_spec
  prepared <- if (is.null(data_spec)) {
    list(subjects = test_subjects, visits = test_visits)
  } else prepare_silk_inputs(test_subjects, test_visits, data_spec)
  test_subjects <- prepared$subjects
  test_visits <- prepared$visits
  pred_shift <- predict_registration_shift(fit$registration$final_template, test_visits, fit$grid)
  stage <- test_subjects$A_obs[match(pred_shift$id, test_subjects$id)] - pred_shift$e_hat
  stage <- stage[match(test_subjects$id, pred_shift$id)]
  x <- base_covariates(
    test_subjects,
    covariate_cols = if (!is.null(data_spec)) data_spec$covariate_cols else NULL
  )
  risk <- predict_age_scale_cox_risk(fit$fit, stage, x, horizons)
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

# Oracle benchmarks are simulation-only and implemented in
# example/sim codes/R/methods_survival_benchmarks.R.
