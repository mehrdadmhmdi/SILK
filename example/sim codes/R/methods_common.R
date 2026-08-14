# =============================================================================
# methods_common.R
# Shared helpers for prediction-focused comparison methods.
# =============================================================================

suppressPackageStartupMessages(library(survival))

clip_probability <- function(x, eps = 1e-6) {
  y <- pmin(pmax(x, eps), 1 - eps)
  dim(y) <- dim(x)
  y
}

first_biomarker_col <- function(visits) {
  bio_columns(visits)[1]
}

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

prediction_frame <- function(subjects, landmark, horizons, risk_mat, context) {
  horizons <- as.numeric(horizons)
  risk_mat <- as.matrix(risk_mat)
  if (nrow(risk_mat) != nrow(subjects) || ncol(risk_mat) != length(horizons)) {
    stop("risk_mat must have nrow(subjects) rows and length(horizons) columns.", call. = FALSE)
  }
  out <- data.frame(
    subject_id = rep(subjects$id, each = length(horizons)),
    landmark = rep(as.numeric(landmark), each = length(horizons)),
    horizon = rep(horizons, times = nrow(subjects)),
    method = context$method,
    risk_pred = clip_probability(as.vector(t(risk_mat))),
    event_within_horizon = as.integer(as.vector(t(outer(subjects$D_star, horizons, "<=")))),
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

apply_standardizer <- function(x, model) {
  x <- as.matrix(x)
  if (ncol(x) == 0L) return(x)
  x <- sweep(sweep(x, 2, model$center, "-"), 2, model$scale, "/")
  x[!is.finite(x)] <- 0
  x
}

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

step_eval_survival <- function(times, surv, grid) {
  times <- as.numeric(times)
  surv <- as.numeric(surv)
  grid <- as.numeric(grid)
  keep <- is.finite(times) & is.finite(surv)
  if (!any(keep)) return(rep(1, length(grid)))
  z <- data.frame(time = c(0, times[keep]), surv = c(1, surv[keep]))
  z <- z[order(z$time), , drop = FALSE]
  z <- z[!duplicated(z$time, fromLast = TRUE), , drop = FALSE]
  out <- rep(1, length(grid))
  valid <- is.finite(grid)
  if (any(valid)) {
    index <- findInterval(grid[valid], z$time)
    out[valid] <- z$surv[pmax(index, 1L)]
  }
  out
}

fit_km_model <- function(time, status) {
  fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  list(time = fit$time, surv = fit$surv)
}

predict_km_risk <- function(km, horizons, n) {
  s <- step_eval_survival(km$time, km$surv, horizons)
  matrix(rep(1 - s, n), nrow = n, byrow = TRUE)
}

make_model_frame <- function(time, status, x) {
  x <- as.matrix(x)
  df <- data.frame(time = time, status = status, x, check.names = FALSE)
  if (ncol(x) > 0L) names(df)[-(1:2)] <- paste0("x", seq_len(ncol(x)))
  df
}

fit_residual_cox <- function(subjects, x) {
  x <- as.matrix(x)
  keep <- if (ncol(x) > 0L) apply(x, 2, stats::sd, na.rm = TRUE) > 1e-10 else logical(0)
  x <- x[, keep, drop = FALSE]
  km <- fit_km_model(subjects$U, subjects$delta)
  if (ncol(x) == 0L) {
    return(list(type = "km", km = km, keep = keep, center = numeric(0), scale = numeric(0)))
  }
  std <- standardize_matrix(x)
  df <- make_model_frame(subjects$U, subjects$delta, std$train)
  form <- stats::as.formula(paste("Surv(time, status) ~", paste(names(df)[-(1:2)], collapse = "+")))
  fit <- tryCatch(
    survival::coxph(form, data = df, ties = "breslow", x = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) || any(!is.finite(stats::coef(fit)))) {
    return(list(type = "km", km = km, keep = keep, center = numeric(0), scale = numeric(0)))
  }
  bh <- survival::basehaz(fit, centered = FALSE)
  list(type = "cox_residual", fit = fit, basehaz = bh, keep = keep,
       center = std$center, scale = std$scale, km = km)
}

predict_residual_cox_risk <- function(model, x_new, horizons) {
  x_new <- as.matrix(x_new)
  if (length(model$keep)) x_new <- x_new[, model$keep, drop = FALSE]
  if (identical(model$type, "km") || ncol(x_new) == 0L) {
    return(predict_km_risk(model$km, horizons, nrow(x_new)))
  }
  x_new <- apply_standardizer(x_new, model)
  H0 <- step_eval(model$basehaz$time, model$basehaz$hazard, horizons)
  beta <- as.numeric(stats::coef(model$fit))
  lp <- as.numeric(x_new %*% beta)
  H <- outer(exp(lp), H0, "*")
  clip_probability(1 - exp(-H))
}

fit_beran_risk <- function(subjects, state, bandwidth = NULL) {
  state <- as.matrix(state)
  keep <- stats::complete.cases(state) &
    is.finite(subjects$U) & is.finite(subjects$delta)
  state <- state[keep, , drop = FALSE]
  time <- as.numeric(subjects$U[keep])
  status <- as.integer(subjects$delta[keep])
  std <- standardize_matrix(state)
  d <- max(1L, ncol(std$train))
  n <- max(1L, nrow(std$train))
  if (is.null(bandwidth)) bandwidth <- n^(-1 / (d + 4))
  bandwidth <- as.numeric(bandwidth)[1]
  if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- n^(-1 / (d + 4))
  list(
    type = "beran",
    time = time,
    status = status,
    state = std$train,
    center = std$center,
    scale = std$scale,
    bandwidth = bandwidth,
    km = fit_km_model(time, status)
  )
}

predict_beran_risk <- function(model, state_new, horizons) {
  state_new <- as.matrix(state_new)
  if (!identical(model$type, "beran") || !nrow(model$state)) {
    return(predict_km_risk(model$km, horizons, nrow(state_new)))
  }
  state_new <- apply_standardizer(state_new, model)
  horizons <- as.numeric(horizons)
  out <- matrix(NA_real_, nrow = nrow(state_new), ncol = length(horizons))
  event_times <- sort(unique(model$time[model$status == 1L & model$time <= max(horizons)]))
  h <- model$bandwidth
  for (ii in seq_len(nrow(state_new))) {
    z <- sweep(model$state, 2, state_new[ii, ], "-")
    w <- exp(-0.5 * rowSums(z^2) / (h^2))
    w[!is.finite(w)] <- 0
    if (sum(w) <= 1e-10 || !length(event_times)) {
      out[ii, ] <- predict_km_risk(model$km, horizons, 1L)
      next
    }
    surv <- rep(1, length(horizons))
    for (tt in event_times) {
      denom <- sum(w[model$time >= tt])
      numer <- sum(w[model$time == tt & model$status == 1L])
      if (is.finite(denom) && denom > 1e-10) {
        step <- pmin(pmax(1 - numer / denom, 0), 1)
        surv[horizons >= tt] <- surv[horizons >= tt] * step
      }
    }
    out[ii, ] <- 1 - surv
  }
  clip_probability(out)
}

fit_age_scale_cox <- function(subjects, start_age, x) {
  x <- as.matrix(x)
  keep <- if (ncol(x) > 0L) apply(x, 2, stats::sd, na.rm = TRUE) > 1e-10 else logical(0)
  x <- x[, keep, drop = FALSE]
  km <- fit_km_model(subjects$U, subjects$delta)
  if (ncol(x) == 0L) {
    return(list(type = "km", km = km, keep = keep, center = numeric(0), scale = numeric(0)))
  }
  std <- standardize_matrix(x)
  df <- data.frame(
    start = as.numeric(start_age),
    stop = as.numeric(start_age) + subjects$U,
    status = subjects$delta,
    std$train,
    check.names = FALSE
  )
  xnames <- paste0("x", seq_len(ncol(std$train)))
  names(df)[-(1:3)] <- xnames
  form <- stats::as.formula(paste("Surv(start, stop, status) ~", paste(xnames, collapse = "+")))
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

predict_age_scale_cox_risk <- function(model, landmark_age, x_new, horizons) {
  x_new <- as.matrix(x_new)
  if (length(model$keep)) x_new <- x_new[, model$keep, drop = FALSE]
  if (identical(model$type, "km") || ncol(x_new) == 0L) {
    return(predict_km_risk(model$km, horizons, nrow(x_new)))
  }
  x_new <- apply_standardizer(x_new, model)
  beta <- as.numeric(stats::coef(model$fit))
  lp <- as.numeric(x_new %*% beta)
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

base_covariates <- function(subjects, specification = c("correct", "misspecified")) {
  specification <- match.arg(specification)
  columns <- if (specification == "correct") c("X1", "X2", "X3", "X4") else c("X1", "X2")
  missing <- setdiff(columns, names(subjects))
  if (length(missing)) {
    stop("Missing baseline covariates: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- as.matrix(subjects[, columns, drop = FALSE])
  if (specification == "correct") {
    out <- cbind(out, X3_sq = as.numeric(subjects$X3)^2 - 1)
  }
  out
}

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

beran_state_covariates <- function(subjects, stage = NULL, include_x = TRUE,
                                   covariate_specification = "correct") {
  if (is.null(stage)) stage <- subjects$A_obs
  x <- cbind(landmark_age = as.numeric(stage))
  if (isTRUE(include_x)) x <- cbind(x, base_covariates(subjects, covariate_specification))
  x
}

landmark_covariates <- function(subjects, visits = NULL, clock = c("recorded", "latent"),
                                marker = NULL, include_biomarker = TRUE,
                                covariate_specification = "correct") {
  clock <- match.arg(clock)
  landmark <- if (clock == "recorded") subjects$A_obs else subjects$A_star
  x <- cbind(
    landmark_age = landmark,
    base_covariates(subjects, covariate_specification)
  )
  if (isTRUE(include_biomarker)) {
    if (is.null(marker)) {
      if (is.null(visits)) {
        stop("visits are required when include_biomarker = TRUE.", call. = FALSE)
      }
      marker <- current_marker_value(subjects, visits, clock)
    }
    x <- cbind(x, current_biomarker = marker)
  }
  x
}

observed_profile_covariates <- function(subjects, visits = NULL,
                                        include_biomarker = OBSERVED_ML_USE_CURRENT_BIOMARKER,
                                        covariate_specification = "correct") {
  if (isTRUE(include_biomarker) && !is.null(visits)) {
    return(landmark_covariates(
      subjects, visits, clock = "recorded", include_biomarker = TRUE,
      covariate_specification = covariate_specification
    ))
  }
  cbind(
    landmark_age = subjects$A_obs,
    base_covariates(subjects, covariate_specification)
  )
}

align_covariate_columns <- function(x, column_names) {
  x <- as.data.frame(x, check.names = FALSE)
  missing <- setdiff(column_names, names(x))
  for (nm in missing) x[[nm]] <- 0
  extra <- setdiff(names(x), column_names)
  if (length(extra)) x <- x[, setdiff(names(x), extra), drop = FALSE]
  x[, column_names, drop = FALSE]
}

fit_current_value_mixed_model <- function(subjects, visits,
                                          specification = c("correct", "misspecified")) {
  specification <- match.arg(specification)
  b <- first_biomarker_col(visits)
  dat <- visits
  dat$marker_value <- dat[[b]]
  if (specification == "correct") {
    missing <- setdiff(c("X1", "X2", "X3", "X4"), names(dat))
    if (length(missing)) {
      stop("Correct mixed model is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    dat$X3_sq <- dat$X3^2 - 1
  }
  dat$id <- factor(dat$id)
  fixed_formula <- if (specification == "correct") {
    marker_value ~ A_obs_il + X1 + X2 + X3 + X4 + X3_sq
  } else {
    marker_value ~ A_obs_il + X1 + X2
  }
  fit <- NULL
  if (requireNamespace("nlme", quietly = TRUE)) {
    fit <- tryCatch(
      suppressWarnings(
        nlme::lme(
          fixed_formula,
          random = ~ A_obs_il | id,
          data = dat,
          control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, returnObject = TRUE)
        )
      ),
      error = function(e) NULL
    )
  }
  fixed_effects <- NULL
  random_cov <- NULL
  residual_sd <- NA_real_
  if (!is.null(fit)) {
    fixed_effects <- tryCatch(as.numeric(nlme::fixef(fit)), error = function(e) NULL)
    if (!is.null(fixed_effects)) names(fixed_effects) <- names(nlme::fixef(fit))
    random_cov <- tryCatch(as.matrix(nlme::getVarCov(fit, type = "random.effects")), error = function(e) NULL)
    if (!is.null(random_cov) && (is.null(colnames(random_cov)) || any(!nzchar(colnames(random_cov))))) {
      rn <- c("(Intercept)", "A_obs_il")[seq_len(ncol(random_cov))]
      rownames(random_cov) <- rn
      colnames(random_cov) <- rn
    }
    residual_sd <- tryCatch(as.numeric(fit$sigma), error = function(e) NA_real_)
  }
  list(
    fit = fit,
    specification = specification,
    biomarker = b,
    fixed_effects = fixed_effects,
    random_cov = random_cov,
    residual_sd = residual_sd
  )
}

predict_current_value_mixed_model <- function(model, subjects, visits) {
  locf <- current_marker_value(subjects, visits, "recorded")
  if (is.null(model$fit) || is.null(model$fixed_effects) || is.null(model$random_cov)) return(locf)

  beta <- model$fixed_effects
  D <- model$random_cov
  sigma2 <- model$residual_sd^2
  if (!is.finite(sigma2) || sigma2 <= 0) sigma2 <- 1
  if (is.null(colnames(D))) {
    z_names <- c("(Intercept)", "A_obs_il")[seq_len(ncol(D))]
    rownames(D) <- z_names
    colnames(D) <- z_names
  }
  z_names <- colnames(D)
  b <- model$biomarker
  pred <- locf

  fixed_design <- function(a_obs, data) {
    n <- length(a_obs)
    value <- function(name) {
      if (name %in% names(data)) as.numeric(data[[name]]) else rep(0, n)
    }
    x3 <- value("X3")
    X <- cbind(
      "(Intercept)" = 1,
      A_obs_il = as.numeric(a_obs),
      X1 = value("X1"),
      X2 = value("X2"),
      X3 = x3,
      X4 = value("X4"),
      X3_sq = x3^2 - 1
    )
    X[, names(beta), drop = FALSE]
  }

  random_design <- function(a_obs) {
    Z <- cbind(
      "(Intercept)" = 1,
      A_obs_il = as.numeric(a_obs)
    )
    Z[, z_names, drop = FALSE]
  }

  for (i in seq_len(nrow(subjects))) {
    sv <- visits[visits$id == subjects$id[i], , drop = FALSE]
    sv <- sv[sv$A_obs_il <= subjects$A_obs[i] + 1e-8, , drop = FALSE]
    if (!nrow(sv)) sv <- visits[visits$id == subjects$id[i], , drop = FALSE]
    if (!nrow(sv) || !b %in% names(sv)) next

    y <- as.numeric(sv[[b]])
    X <- fixed_design(sv$A_obs_il, sv)
    Z <- random_design(sv$A_obs_il)
    keep <- is.finite(y) & apply(X, 1, function(row) all(is.finite(row))) &
      apply(Z, 1, function(row) all(is.finite(row)))
    if (!any(keep)) next

    y <- y[keep]
    X <- X[keep, , drop = FALSE]
    Z <- Z[keep, , drop = FALSE]
    eta <- as.numeric(X %*% beta)
    resid <- y - eta
    V <- Z %*% D %*% t(Z) + diag(sigma2, nrow(Z))
    bhat <- tryCatch(
      as.numeric(D %*% t(Z) %*% solve(V, resid)),
      error = function(e) rep(0, ncol(D))
    )
    X0 <- fixed_design(subjects$A_obs[i], subjects[i, , drop = FALSE])
    Z0 <- random_design(subjects$A_obs[i])
    val <- as.numeric(X0 %*% beta + Z0 %*% bhat)
    if (is.finite(val)) pred[i] <- val
  }

  pred[!is.finite(pred)] <- locf[!is.finite(pred)]
  pred
}

draw_time_error <- function(subjects, n_impute, sd = TIMEERROR_SD, use_truth = TIMEERROR_USE_TRUTH) {
  n <- nrow(subjects)
  if (isTRUE(use_truth) && "eps" %in% names(subjects)) {
    sd <- stats::sd(subjects$eps)
    if (!is.finite(sd) || sd <= 0) sd <- TIMEERROR_SD
  }
  matrix(stats::rnorm(n * n_impute, mean = 0, sd = sd), nrow = n, ncol = n_impute)
}

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
