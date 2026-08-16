# =============================================================================
# utils-common.R
# Shared helpers (from methods_common.R)
# =============================================================================

#' Clip probabilities to (eps, 1-eps)
#'
#' @param x Numeric vector or matrix of probabilities.
#' @param eps Numeric. Small positive number for clipping.
#' @return Clipped probabilities with same dimensions as input.
#' @export
clip_probability <- function(x, eps = 1e-6) {
  y <- pmin(pmax(x, eps), 1 - eps)
  dim(y) <- dim(x)
  y
}

#' @keywords internal
first_biomarker_col <- function(visits) {
  bio_columns(visits)[1]
}

#' Create a method context object
#'
#' @param method Character. Method name.
#' @param fold_id Integer. Fold identifier.
#' @param replicate_id Integer. Replicate identifier.
#' @param n_train_setting Numeric. Training set size setting.
#' @param time_grid_setting Character. Time grid setting identifier.
#' @return A list with method context fields.
#' @keywords internal
method_context <- function(method, fold_id = 1L, replicate_id = 1L,
                           n_train_setting = NA, time_grid_setting = NA) {
  list(
    method = method,
    fold_id = fold_id,
    replicate_id = replicate_id,
    n_train_setting = n_train_setting,
    time_grid_setting = time_grid_setting
  )
}

#' Build a prediction frame
#'
#' Creates a standardized data frame of risk predictions across subjects and horizons.
#'
#' @param subjects Data frame with column id.  Simulation data may also carry
#'   D_star; real-data predictions return NA for the outcome label because
#'   censoring-aware evaluation belongs to the application code.
#' @param landmark Numeric vector of landmark ages.
#' @param horizons Numeric vector of prediction horizons.
#' @param risk_mat Matrix of predicted risks (subjects x horizons).
#' @param context Method context list from \code{method_context}.
#' @return Data frame with columns: subject_id, landmark, horizon, method,
#'   risk_pred, event_within_horizon, at_risk, fold_id, replicate_id,
#'   n_train_setting, time_grid_setting.
#' @export
prediction_frame <- function(subjects, landmark, horizons, risk_mat, context) {
  horizons <- as.numeric(horizons)
  risk_mat <- as.matrix(risk_mat)
  if (nrow(risk_mat) != nrow(subjects) || ncol(risk_mat) != length(horizons)) {
    stop("risk_mat must have nrow(subjects) rows and length(horizons) columns.", call. = FALSE)
  }
  event_within_horizon <- if ("D_star" %in% names(subjects)) {
    as.integer(as.vector(t(outer(subjects$D_star, horizons, "<="))))
  } else {
    # Real applications generally require IPCW evaluation and do not have
    # the simulation-only latent event time D_star.  Keep the prediction API
    # usable by returning an explicit missing outcome label rather than
    # silently treating censored observations as events/non-events.
    rep(NA_integer_, nrow(subjects) * length(horizons))
  }
  out <- data.frame(
    subject_id = rep(subjects$id, each = length(horizons)),
    landmark = rep(as.numeric(landmark), each = length(horizons)),
    horizon = rep(horizons, times = nrow(subjects)),
    method = context$method,
    risk_pred = clip_probability(as.vector(t(risk_mat))),
    event_within_horizon = event_within_horizon,
    at_risk = TRUE,
    fold_id = context$fold_id,
    replicate_id = context$replicate_id,
    n_train_setting = context$n_train_setting,
    time_grid_setting = context$time_grid_setting,
    stringsAsFactors = FALSE
  )
  out[, c(
    "subject_id", "landmark", "horizon", "method", "risk_pred",
    "event_within_horizon", "at_risk", "fold_id", "replicate_id",
    "n_train_setting", "time_grid_setting"
  ), drop = FALSE]
}

#' Standardize a matrix using training statistics
#'
#' @param train_x Numeric matrix. Training data.
#' @param test_x Numeric matrix. Test data (optional; defaults to train_x).
#' @return List with elements train, test, center, scale.
#' @keywords internal
standardize_matrix <- function(train_x, test_x = NULL) {
  train_x <- as.matrix(train_x)
  if (is.null(test_x)) test_x <- train_x
  test_x <- as.matrix(test_x)
  if (ncol(train_x) == 0L) {
    return(list(train = train_x, test = test_x, center = numeric(0), scale = numeric(0)))
  }
  center <- colMeans(train_x, na.rm = TRUE)
  scale <- apply(train_x, 2, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  train_std <- sweep(sweep(train_x, 2, center, "-"), 2, scale, "/")
  test_std <- sweep(sweep(test_x, 2, center, "-"), 2, scale, "/")
  train_std[!is.finite(train_std)] <- 0
  test_std[!is.finite(test_std)] <- 0
  list(train = train_std, test = test_std, center = center, scale = scale)
}

#' Apply standardization using stored center/scale
#' @param x Numeric matrix.
#' @param model List with center and scale vectors.
#' @return Standardized matrix.
#' @keywords internal
apply_standardizer <- function(x, model) {
  x <- as.matrix(x)
  if (ncol(x) == 0L) return(x)
  x <- sweep(sweep(x, 2, model$center, "-"), 2, model$scale, "/")
  x[!is.finite(x)] <- 0
  x
}

#' Step-function evaluation
#' @param times Numeric vector.
#' @param values Numeric vector.
#' @param grid Numeric vector.
#' @return Numeric vector.
#' @keywords internal
step_eval <- function(times, values, grid) {
  times <- as.numeric(times)
  values <- as.numeric(values)
  grid <- as.numeric(grid)
  keep <- is.finite(times) & is.finite(values)
  if (!any(keep)) return(rep(0, length(grid)))
  z <- data.frame(time = c(0, times[keep]), value = c(0, values[keep]))
  z <- z[order(z$time), , drop = FALSE]
  z <- z[!duplicated(z$time, fromLast = TRUE), , drop = FALSE]
  out <- rep(0, length(grid))
  valid_grid <- is.finite(grid)
  if (any(valid_grid)) {
    idx <- findInterval(grid[valid_grid], z$time)
    out[valid_grid] <- z$value[pmax(idx, 1L)]
  }
  out
}

#' @keywords internal
fit_km_model <- function(time, status) {
  fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  list(time = fit$time, surv = fit$surv)
}

#' @keywords internal
predict_km_risk <- function(km, horizons, n) {
  s <- step_eval_survival(km$time, km$surv, horizons)
  matrix(rep(1 - s, n), nrow = n, byrow = TRUE)
}

#' @keywords internal
make_model_frame <- function(time, status, x) {
  x <- as.matrix(x)
  df <- data.frame(time = time, status = status, x, check.names = FALSE)
  if (ncol(x) > 0L) names(df)[-(1:2)] <- paste0("x", seq_len(ncol(x)))
  df
}

#' @keywords internal
fit_age_scale_cox <- function(subjects, start_age, x) {
  x <- as.matrix(x)
  keep <- if (ncol(x) > 0L) apply(x, 2, stats::sd, na.rm = TRUE) > 1e-10 else logical(0)
  x <- x[, keep, drop = FALSE]
  km <- fit_km_model(subjects$U, subjects$delta)
  df <- data.frame(
    start = as.numeric(start_age),
    stop = as.numeric(start_age) + subjects$U,
    status = subjects$delta,
    check.names = FALSE
  )
  if (ncol(x) == 0L) {
    fit <- tryCatch(
      survival::coxph(
        survival::Surv(start, stop, status) ~ 1,
        data = df, ties = "breslow", x = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(type = "km", km = km, keep = keep, center = numeric(0), scale = numeric(0)))
    }
    bh <- survival::basehaz(fit, centered = FALSE)
    return(list(
      type = "cox_age_null", fit = fit, basehaz = bh, keep = keep,
      center = numeric(0), scale = numeric(0), km = km
    ))
  }
  std <- standardize_matrix(x)
  df <- cbind(df, as.data.frame(std$train, check.names = FALSE))
  xnames <- paste0("x", seq_len(ncol(std$train)))
  names(df)[-(1:3)] <- xnames
  form <- stats::as.formula(paste("survival::Surv(start, stop, status) ~", paste(xnames, collapse = "+")))
  fit <- tryCatch(
    survival::coxph(form, data = df, ties = "breslow", x = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) || any(!is.finite(stats::coef(fit)))) {
    return(list(type = "km", km = km, keep = keep, center = numeric(0), scale = numeric(0)))
  }
  bh <- survival::basehaz(fit, centered = FALSE)
  list(type = "cox_age", fit = fit, basehaz = bh, keep = keep,
       center = std$center, scale = std$scale, km = km)
}

#' @keywords internal
predict_age_scale_cox_risk <- function(model, landmark_age, x_new, horizons) {
  x_new <- as.matrix(x_new)
  if (length(model$keep)) x_new <- x_new[, model$keep, drop = FALSE]
  if (identical(model$type, "km")) {
    return(predict_km_risk(model$km, horizons, nrow(x_new)))
  }
  if (identical(model$type, "cox_age_null")) {
    lp <- rep(0, nrow(x_new))
  } else {
    if (ncol(x_new) == 0L) {
      stop("The fitted age-scale Cox model requires covariates.", call. = FALSE)
    }
    x_new <- apply_standardizer(x_new, model)
    beta <- as.numeric(stats::coef(model$fit))
    lp <- as.numeric(x_new %*% beta)
  }
  all_times <- sort(unique(c(as.numeric(landmark_age), as.numeric(landmark_age) + rep(horizons, each = length(landmark_age)))))
  H_all <- step_eval(model$basehaz$time, model$basehaz$hazard, all_times)
  H_lookup <- stats::setNames(H_all, as.character(all_times))
  risk <- matrix(NA_real_, nrow = nrow(x_new), ncol = length(horizons))
  for (h in seq_along(horizons)) {
    H_start <- as.numeric(H_lookup[as.character(as.numeric(landmark_age))])
    H_stop <- as.numeric(H_lookup[as.character(as.numeric(landmark_age) + horizons[h])])
    dH <- pmax(H_stop - H_start, 0)
    risk[, h] <- 1 - exp(-dH * exp(lp))
  }
  clip_probability(risk)
}

#' @keywords internal
base_covariates <- function(subjects, specification = c("correct", "misspecified"),
                             covariate_cols = NULL) {
  specification <- match.arg(specification)
  columns <- if (!is.null(covariate_cols)) {
    unique(as.character(covariate_cols))
  } else {
    unique(as.character(silk_opt("SURVIVAL_COVARIATES")))
  }
  columns <- columns[nzchar(columns)]
  missing <- setdiff(columns, names(subjects))
  if (length(missing)) {
    stop("Missing baseline covariates: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!length(columns)) return(matrix(numeric(0), nrow = nrow(subjects), ncol = 0L))
  as.matrix(subjects[, columns, drop = FALSE])
}

#' @keywords internal
current_marker_value <- function(subjects, visits, clock = c("recorded", "latent")) {
  clock <- match.arg(clock)
  b <- first_biomarker_col(visits)
  ids <- subjects$id
  out <- numeric(length(ids))
  for (i in seq_along(ids)) {
    sv <- visits[visits$id == ids[i], , drop = FALSE]
    if (!nrow(sv)) {
      out[i] <- NA_real_
      next
    }
    landmark <- if (clock == "recorded") subjects$A_obs[i] else subjects$A_star[i]
    time_col <- if (clock == "recorded") "A_obs_il" else "A_star_il"
    sv <- sv[sv[[time_col]] <= landmark + 1e-8, , drop = FALSE]
    if (!nrow(sv)) sv <- visits[visits$id == ids[i], , drop = FALSE]
    sv <- sv[order(sv[[time_col]], decreasing = TRUE), , drop = FALSE]
    out[i] <- sv[[b]][1]
  }
  out[!is.finite(out)] <- 0
  out
}

#' @keywords internal
align_covariate_columns <- function(x, column_names) {
  x <- as.data.frame(x, check.names = FALSE)
  missing <- setdiff(column_names, names(x))
  for (nm in missing) x[[nm]] <- 0
  extra <- setdiff(names(x), column_names)
  if (length(extra)) x <- x[, setdiff(names(x), extra), drop = FALSE]
  x[, column_names, drop = FALSE]
}

#' @keywords internal
draw_time_error <- function(subjects, n_impute, sd = NULL, use_truth = FALSE) {
  if (is.null(sd)) sd <- 1.5
  n <- nrow(subjects)
  if (isTRUE(use_truth) && "eps" %in% names(subjects)) {
    sd <- stats::sd(subjects$eps)
    if (!is.finite(sd) || sd <= 0) sd <- 1.5
  }
  matrix(stats::rnorm(n * n_impute, mean = 0, sd = sd), nrow = n, ncol = n_impute)
}

#' @keywords internal
shift_recorded_time <- function(subjects, visits, eps_draw) {
  subjects2 <- subjects
  visits2 <- visits
  eps_by_id <- stats::setNames(as.numeric(eps_draw), subjects$id)
  subjects2$A_obs <- subjects2$A_obs - as.numeric(eps_draw)
  for (nm in intersect(c("T_obs", "C_obs"), names(subjects2))) {
    subjects2[[nm]] <- subjects2[[nm]] - as.numeric(eps_draw)
  }
  if (!is.null(visits2)) {
    visits2$A_obs_il <- visits2$A_obs_il - as.numeric(eps_by_id[as.character(visits2$id)])
  }
  list(subjects = subjects2, visits = visits2)
}
