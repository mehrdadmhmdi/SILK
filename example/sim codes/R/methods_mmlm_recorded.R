# =============================================================================
# methods_mmlm_recorded.R
# MMLM: correct and deliberately misspecified mixed-model landmarking.
# A mixed model smooths the current biomarker value only; the landmark Cox model
# uses the prespecified event covariates and that smoothed value.
# Reference: Rizopoulos et al. dynamic prediction work; nlme::lme and survival.
# =============================================================================

fit_mmlm_misspecified_event_model <- function(subjects, x,
                                               theta = MISSPECIFIED_ASSOCIATION_SHRINKAGE) {
  x <- as.matrix(x)
  required <- c("X1", "X2", "current_biomarker")
  missing <- setdiff(required, colnames(x))
  if (length(missing)) {
    stop("Misspecified MMLM event design is missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  x <- x[, required, drop = FALSE]
  keep <- apply(x, 2, stats::sd, na.rm = TRUE) > 1e-10
  x <- x[, keep, drop = FALSE]
  km <- fit_km_model(subjects$U, subjects$delta)
  if (!all(required %in% colnames(x))) {
    return(list(type = "km", km = km, keep = keep,
                center = numeric(0), scale = numeric(0)))
  }

  std <- standardize_matrix(x)
  dat <- data.frame(
    time = subjects$U,
    status = subjects$delta,
    std$train,
    check.names = FALSE
  )
  fit <- tryCatch(
    survival::coxph(
      survival::Surv(time, status) ~ X1 + X2 +
        survival::ridge(current_biomarker, theta = theta),
      data = dat,
      ties = "breslow",
      x = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || any(!is.finite(stats::coef(fit)))) {
    return(list(type = "km", km = km, keep = keep,
                center = numeric(0), scale = numeric(0)))
  }
  list(
    type = "cox_residual_ridge",
    fit = fit,
    basehaz = survival::basehaz(fit, centered = FALSE),
    keep = keep,
    center = std$center,
    scale = std$scale,
    km = km,
    theta = theta
  )
}

predict_mmlm_misspecified_event_risk <- function(model, x_new, horizons) {
  x_new <- as.matrix(x_new)
  if (length(model$keep)) x_new <- x_new[, model$keep, drop = FALSE]
  if (!identical(model$type, "cox_residual_ridge") || !ncol(x_new)) {
    return(predict_km_risk(model$km, horizons, nrow(x_new)))
  }
  x_new <- apply_standardizer(x_new, model)
  newdata <- as.data.frame(x_new, check.names = FALSE)
  lp <- as.numeric(stats::predict(
    model$fit, newdata = newdata, type = "lp", reference = "zero"
  ))
  H0 <- step_eval(model$basehaz$time, model$basehaz$hazard, horizons)
  clip_probability(1 - exp(-outer(exp(lp), H0, "*")))
}

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
  if (specification == "misspecified") {
    # The misspecified event model omits the attained-age effect together with
    # X3, X4, and X3-squared. It therefore uses only X1, X2, and a ridge-shrunk
    # current value of first-differenced B1 on the residual follow-up time scale.
    x <- x[, setdiff(colnames(x), "landmark_age"), drop = FALSE]
  }
  fit <- if (specification == "correct") {
    fit_residual_cox(train_subjects, x)
  } else {
    fit_mmlm_misspecified_event_model(train_subjects, x)
  }
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
  if (fit$specification == "misspecified") {
    x <- x[, setdiff(colnames(x), "landmark_age"), drop = FALSE]
  }
  risk <- if (fit$specification == "correct") {
    predict_residual_cox_risk(fit$fit, x, horizons)
  } else {
    predict_mmlm_misspecified_event_risk(fit$fit, x, horizons)
  }
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
