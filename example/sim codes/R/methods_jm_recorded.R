# =============================================================================
# methods_jm_recorded.R
# JM: correct and deliberately misspecified current-value joint models.
# Uses a simple current-value association between one longitudinal biomarker and
# the event process. Internally, time is expressed as positive follow-up time
# relative to each subject's recorded landmark. This is the usual landmark
# follow-up formulation and avoids negative-time failures in held-out subjects.
# No slope association is included.
# References: Rizopoulos (2012, 2017); JMbayes2 package; nlme and survival.
# =============================================================================

prepare_jm_followup_clock <- function(subjects, visits,
                                      landmark_time = JMBAYES_LANDMARK_TIME) {
  subject_index <- match(visits$id, subjects$id)
  if (anyNA(subject_index)) {
    stop("JM visits contain subject IDs absent from the subject table.", call. = FALSE)
  }

  landmark <- as.numeric(subjects$A_obs)
  visit_landmark <- landmark[subject_index]
  out_subjects <- subjects
  out_visits <- visits
  out_visits$A_obs_il <- as.numeric(visits$A_obs_il) - visit_landmark + landmark_time
  out_subjects$T_obs <- landmark_time + as.numeric(subjects$T_obs) - landmark
  out_subjects$C_obs <- landmark_time + as.numeric(subjects$C_obs) - landmark
  out_subjects$A_obs <- landmark_time

  internal_times <- c(
    out_visits$A_obs_il,
    out_subjects$A_obs,
    out_subjects$T_obs,
    out_subjects$C_obs
  )
  if (any(!is.finite(internal_times))) {
    stop("JM follow-up-clock transformation produced non-finite times.", call. = FALSE)
  }
  if (any(pmin(out_subjects$T_obs, out_subjects$C_obs) <= landmark_time)) {
    stop("JM requires event/censoring times strictly after the landmark.", call. = FALSE)
  }

  list(subjects = out_subjects, visits = out_visits)
}

fit_jm_longitudinal_model <- function(visits,
                                      specification = c("correct", "misspecified")) {
  specification <- match.arg(specification)
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("JM requires the nlme package.", call. = FALSE)
  }
  biomarker <- first_biomarker_col(visits)
  dat <- visits
  dat$marker_value <- marker_values_for_specification(dat, biomarker, specification)
  fixed_covariates <- if (specification == "correct") {
    c("X1", "X2", "X3", "X4", "X3_sq")
  } else {
    c("X1", "X2")
  }
  if (specification == "correct") dat$X3_sq <- dat$X3^2 - 1
  required <- c("id", "marker_value", "A_obs_il", fixed_covariates)
  dat <- dat[stats::complete.cases(dat[, required, drop = FALSE]), , drop = FALSE]
  dat$id <- droplevels(factor(dat$id))
  if (!nrow(dat) || length(unique(dat$id)) < 2L) {
    stop("JM has insufficient complete longitudinal observations.", call. = FALSE)
  }
  random_time_center <- NA_real_
  if (specification == "misspecified") {
    # Centering is an exact reparameterization of an unstructured Gaussian
    # random intercept/slope model. It avoids near-singular mixed-model systems
    # caused by using absolute follow-up-clock values (roughly 6--10) after the
    # first observation is removed by differencing.
    random_time_center <- mean(dat$A_obs_il)
    dat$A_obs_centered <- dat$A_obs_il - random_time_center
  }

  fixed_formula <- if (specification == "correct") {
    marker_value ~ A_obs_il + X1 + X2 + X3 + X4 + X3_sq
  } else {
    marker_value ~ A_obs_il + X1 + X2
  }

  fit_one <- function(random_formula) {
    tryCatch(
      suppressWarnings(nlme::lme(
        fixed_formula,
        random = random_formula,
        data = dat,
        na.action = stats::na.omit,
        control = nlme::lmeControl(
          maxIter = 100,
          msMaxIter = 100,
          niterEM = 50,
          returnObject = TRUE
        )
      )),
      error = function(e) NULL
    )
  }

  if (specification == "correct") {
    fit <- fit_one(~ A_obs_il | id)
    random_structure <- "random intercept and recorded-time slope"
    if (is.null(fit)) {
      fit <- fit_one(~ 1 | id)
      random_structure <- "random intercept fallback"
    }
  } else {
    # JMbayes2 does not support a single longitudinal outcome with only random
    # intercepts.  The supported Gaussian intercept-slope working model remains
    # deliberately wrong here because a first-differenced B1 trajectory is
    # treated as Gaussian, while X3/X4/quadratic effects are omitted from the
    # event model and the marker association is strongly shrunk toward zero.
    fit <- fit_one(~ A_obs_centered | id)
    random_structure <- "misspecified Gaussian random intercept and recorded-time slope"
    if (is.null(fit)) {
      # Retain both Gaussian random effects but remove their estimated
      # correlation as a last-resort numerical fallback.
      fit <- fit_one(list(id = nlme::pdDiag(~ A_obs_centered)))
      random_structure <- paste(
        "misspecified independent Gaussian random intercept and",
        "recorded-time slope fallback"
      )
    }
  }
  if (is.null(fit)) {
    stop("JM could not fit either longitudinal mixed model.", call. = FALSE)
  }
  list(
    fit = fit,
    biomarker = biomarker,
    random_structure = random_structure,
    specification = specification,
    fixed_covariates = fixed_covariates,
    marker_transform = if (specification == "correct") "raw" else "first_difference",
    random_time_center = random_time_center
  )
}

jm_association_priors <- function(specification = c("correct", "misspecified")) {
  specification <- match.arg(specification)
  if (specification == "correct") return(NULL)
  list(
    mean_alphas = list(0),
    Tau_alphas = list(matrix(MISSPECIFIED_ASSOCIATION_SHRINKAGE, 1L, 1L))
  )
}

fit_jm_recorded <- function(train_subjects, train_visits,
                            specification = c("correct", "misspecified")) {
  specification <- match.arg(specification)
  method <- if (specification == "correct") "JM-Correct" else "JM-Misspecified"
  if (!isTRUE(ENABLE_JMBAYES2)) {
    stop("JM requires SILK_ENABLE_JMBAYES2=true; method marked unavailable for this run.", call. = FALSE)
  }
  if (!requireNamespace("JMbayes2", quietly = TRUE)) {
    stop("JM requires the JMbayes2 package, which is not installed.", call. = FALSE)
  }

  jm_data <- prepare_jm_followup_clock(train_subjects, train_visits)
  jm_subjects <- jm_data$subjects
  jm_visits <- jm_data$visits
  marker_model <- fit_jm_longitudinal_model(jm_visits, specification)

  event_dat <- data.frame(
    id = jm_subjects$id,
    Y_obs = pmin(jm_subjects$T_obs, jm_subjects$C_obs),
    status = jm_subjects$delta,
    X1 = jm_subjects$X1,
    X2 = jm_subjects$X2,
    X3 = jm_subjects$X3,
    X4 = jm_subjects$X4,
    X3_sq = jm_subjects$X3^2 - 1
  )
  event_formula <- if (specification == "correct") {
    survival::Surv(Y_obs, status) ~ X1 + X2 + X3 + X4 + X3_sq
  } else {
    survival::Surv(Y_obs, status) ~ X1 + X2
  }
  cox_fit <- survival::coxph(event_formula, data = event_dat, x = TRUE)
  association_priors <- jm_association_priors(specification)
  jm_fit <- JMbayes2::jm(
    cox_fit,
    marker_model$fit,
    time_var = "A_obs_il",
    n_chains = JMBAYES_N_CHAINS,
    n_iter = JMBAYES_N_ITER,
    n_burnin = JMBAYES_N_BURNIN,
    priors = association_priors
  )

  list(
    method = method,
    specification = specification,
    marker_model = marker_model,
    cox_fit = cox_fit,
    jm_fit = jm_fit,
    landmark_time = JMBAYES_LANDMARK_TIME,
    implementation = paste0(
      "JMbayes2 current-value dynamic prediction with the ", specification,
      " covariate specification on a recorded-landmark follow-up clock; ",
      marker_model$random_structure
    )
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
    risk_cols <- setdiff(
      numeric_cols,
      c(time_col, "id", "X1", "X2", "X3", "X4", "X3_sq", "Y_obs", "status", "A_obs_il")
    )
  }
  if (!length(risk_cols)) stop("Could not identify the JMbayes2 event-risk column.", call. = FALSE)
  risk <- as.numeric(df[[risk_cols[1]]])

  if (length(time_col)) {
    tt <- as.numeric(df[[time_col[1]]])
    keep <- is.finite(tt) & is.finite(risk)
    tt <- tt[keep]
    risk <- risk[keep]
    if (!length(tt)) {
      stop("JMbayes2 prediction returned no finite event-risk values.", call. = FALSE)
    }
    ord <- order(tt)
    tt <- tt[ord]
    risk <- risk[ord]
    if (length(tt) == 1L) {
      out <- rep(risk, length(target_times))
    } else {
      out <- stats::approx(tt, risk, xout = target_times, rule = 2, ties = mean)$y
    }
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
  landmark_time <- if (!is.null(fit$landmark_time)) {
    as.numeric(fit$landmark_time)
  } else {
    JMBAYES_LANDMARK_TIME
  }
  risk <- matrix(NA_real_, nrow = nrow(test_subjects), ncol = length(horizons))
  for (i in seq_len(nrow(test_subjects))) {
    subject_id <- test_subjects$id[i]
    sv <- test_visits[test_visits$id == subject_id, , drop = FALSE]
    if (!nrow(sv)) {
      stop("JM has no longitudinal history for test subject ", subject_id, ".", call. = FALSE)
    }
    sv$marker_value <- marker_values_for_specification(
      sv, b, fit$specification
    )
    sv$A_obs_il <- as.numeric(sv$A_obs_il) - as.numeric(test_subjects$A_obs[i]) + landmark_time
    if (identical(fit$specification, "misspecified")) {
      sv$A_obs_centered <- sv$A_obs_il - fit$marker_model$random_time_center
    }
    sv$X1 <- test_subjects$X1[i]
    sv$X2 <- test_subjects$X2[i]
    sv$X3 <- test_subjects$X3[i]
    sv$X4 <- test_subjects$X4[i]
    sv$X3_sq <- sv$X3^2 - 1
    required_finite <- is.finite(sv$A_obs_il)
    for (nm in fit$marker_model$fixed_covariates) {
      required_finite <- required_finite & is.finite(sv[[nm]])
    }
    sv <- sv[required_finite, , drop = FALSE]
    if (!nrow(sv)) {
      stop("JM has no finite longitudinal times for test subject ", subject_id, ".", call. = FALSE)
    }
    sv <- sv[order(sv$A_obs_il), , drop = FALSE]
    finite_marker <- is.finite(sv$marker_value)
    use_y <- any(finite_marker)
    if (use_y) {
      sv <- sv[finite_marker, , drop = FALSE]
    } else {
      sv <- tail(sv, 1L)
      sv$marker_value <- 0
    }
    sv$id <- factor(sv$id)
    # JMbayes2 0.6.0 is most reliable for a single-outcome model when the
    # longitudinal and event variables share one data frame. It takes the last
    # event row per subject internally, while retaining the full marker history.
    sv$Y_obs <- landmark_time
    sv$status <- 0
    target_times <- landmark_time + horizons
    pred <- tryCatch(
      stats::predict(
        fit$jm_fit,
        newdata = sv,
        process = "event",
        times = target_times,
        control = list(
          cores = 1L,
          n_samples = JMBAYES_PRED_N_SAMPLES,
          use_Y = use_y,
          return_newdata = TRUE
        )
      ),
      error = function(e) {
        stop(
          "JM prediction failed for test subject ", subject_id,
          ": ", conditionMessage(e),
          call. = FALSE
        )
      }
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
