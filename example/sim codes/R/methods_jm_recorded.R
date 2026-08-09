# =============================================================================
# methods_jm_recorded.R
# JM-Recorded: recorded-time current-value joint-model comparator.
# Uses a simple current-value association between one longitudinal biomarker and
# the event process on the recorded clock. No slope association is included.
# References: Rizopoulos (2012, 2017); JMbayes2 package; nlme and survival.
# =============================================================================

fit_jm_recorded <- function(train_subjects, train_visits) {
  if (!isTRUE(ENABLE_JMBAYES2)) {
    stop("JM-Recorded requires SILK_ENABLE_JMBAYES2=true; method marked unavailable for this run.", call. = FALSE)
  }
  if (!requireNamespace("JMbayes2", quietly = TRUE)) {
    stop("JM-Recorded requires the JMbayes2 package, which is not installed.", call. = FALSE)
  }

  event_time <- pmin(train_subjects$T_obs, train_subjects$C_obs)
  min_recorded_time <- min(c(event_time, train_subjects$A_obs, train_visits$A_obs_il), na.rm = TRUE)
  time_shift <- if (is.finite(min_recorded_time) && min_recorded_time <= 0) abs(min_recorded_time) + 1e-3 else 0
  jm_subjects <- train_subjects
  jm_visits <- train_visits
  if (time_shift > 0) {
    jm_subjects$A_obs <- jm_subjects$A_obs + time_shift
    jm_subjects$T_obs <- jm_subjects$T_obs + time_shift
    jm_subjects$C_obs <- jm_subjects$C_obs + time_shift
    jm_visits$A_obs_il <- jm_visits$A_obs_il + time_shift
  }

  marker_model <- fit_current_value_mixed_model(jm_subjects, jm_visits)
  if (is.null(marker_model$fit)) {
    stop("JM-Recorded could not fit the recorded-time longitudinal mixed model.", call. = FALSE)
  }

  event_dat <- data.frame(
    id = jm_subjects$id,
    Y_obs = pmin(jm_subjects$T_obs, jm_subjects$C_obs),
    status = jm_subjects$delta,
    X1 = jm_subjects$X1,
    X2 = jm_subjects$X2
  )
  cox_fit <- survival::coxph(survival::Surv(Y_obs, status) ~ X1 + X2, data = event_dat, x = TRUE)
  jm_fit <- JMbayes2::jm(
    cox_fit,
    marker_model$fit,
    time_var = "A_obs_il",
    n_chains = JMBAYES_N_CHAINS,
    n_iter = JMBAYES_N_ITER,
    n_burnin = JMBAYES_N_BURNIN
  )

  list(
    method = "JM-Recorded",
    marker_model = marker_model,
    cox_fit = cox_fit,
    jm_fit = jm_fit,
    time_shift = time_shift,
    implementation = "JMbayes2 current-value dynamic prediction on recorded clock"
  )
}

extract_jmbayes2_event_risk <- function(pred, target_times) {
  if (is.data.frame(pred)) {
    df <- pred
  } else if (is.list(pred)) {
    data_frames <- pred[vapply(pred, is.data.frame, logical(1))]
    if (!length(data_frames)) stop("JMbayes2 prediction did not return a data frame.", call. = FALSE)
    df <- data_frames[[1]]
  } else {
    stop("Unsupported JMbayes2 prediction object.", call. = FALSE)
  }

  time_col <- intersect(c("times", "time", "Time", "A_obs_il"), names(df))
  risk_cols <- grep("risk|event|cif|prob", names(df), ignore.case = TRUE, value = TRUE)
  risk_cols <- risk_cols[vapply(df[risk_cols], is.numeric, logical(1))]
  risk_cols <- setdiff(risk_cols, time_col)
  if (!length(risk_cols)) {
    numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    risk_cols <- setdiff(numeric_cols, c(time_col, "id", "X1", "X2", "Y_obs", "status", "A_obs_il"))
  }
  if (!length(risk_cols)) stop("Could not identify the JMbayes2 event-risk column.", call. = FALSE)
  risk <- as.numeric(df[[risk_cols[1]]])

  if (length(time_col)) {
    tt <- as.numeric(df[[time_col[1]]])
    out <- stats::approx(tt, risk, xout = target_times, rule = 2, ties = "ordered")$y
  } else if (length(risk) >= length(target_times)) {
    out <- tail(risk, length(target_times))
  } else {
    stop("JMbayes2 prediction did not include enough event-risk values.", call. = FALSE)
  }
  clip_probability(out)
}

predict_jm_recorded <- function(fit, test_subjects, test_visits,
                                horizons = PREDICTION_HORIZONS,
                                fold_id = 1L, replicate_id = 1L,
                                n_train_setting = NA,
                                time_grid_setting = NA) {
  b <- fit$marker_model$biomarker
  time_shift <- if (!is.null(fit$time_shift)) fit$time_shift else 0
  risk <- matrix(NA_real_, nrow = nrow(test_subjects), ncol = length(horizons))
  for (i in seq_len(nrow(test_subjects))) {
    sv <- test_visits[test_visits$id == test_subjects$id[i], , drop = FALSE]
    sv$marker_value <- sv[[b]]
    sv$id <- factor(sv$id)
    sv$A_obs_il <- sv$A_obs_il + time_shift
    sv$Y_obs <- test_subjects$A_obs[i] + time_shift
    sv$status <- 0
    sv$X1 <- test_subjects$X1[i]
    sv$X2 <- test_subjects$X2[i]
    target_times <- test_subjects$A_obs[i] + time_shift + horizons
    pred <- stats::predict(
      fit$jm_fit,
      newdata = sv,
      process = "event",
      times = target_times,
      control = list(cores = 1L, n_samples = JMBAYES_PRED_N_SAMPLES, return_newdata = TRUE)
    )
    risk[i, ] <- extract_jmbayes2_event_risk(pred, target_times)
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
