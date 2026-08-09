# =============================================================================
# methods_observed_ml.R
# Observed-clock machine-learning and Bayesian survival comparators.
# These methods use observed landmark age and baseline covariates by default.
# They do not estimate latent age; they are included as observed-age predictive
# benchmarks. Current biomarker level can be enabled only as an explicit
# sensitivity setting via SILK_OBSERVED_ML_USE_CURRENT_BIOMARKER=true.
# =============================================================================

fit_rsf_observed <- function(train_subjects, train_visits, seed = NULL) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("RSF-Observed requires the ranger package, which is not installed.", call. = FALSE)
  }
  x <- observed_profile_covariates(train_subjects, train_visits)
  dat <- data.frame(
    time = train_subjects$U,
    status = train_subjects$delta,
    as.data.frame(x, check.names = FALSE),
    check.names = FALSE
  )
  fit <- ranger::ranger(
    survival::Surv(time, status) ~ .,
    data = dat,
    num.trees = RSF_NUM_TREES,
    min.node.size = RSF_MIN_NODE_SIZE,
    seed = seed,
    num.threads = as.integer(Sys.getenv("SILK_RSF_THREADS", unset = "1")),
    respect.unordered.factors = "order"
  )
  list(
    method = "RSF-Observed",
    fit = fit,
    columns = colnames(x),
    implementation = "ranger::ranger survival forest on observed landmark-age profile"
  )
}

predict_rsf_observed <- function(fit, test_subjects, test_visits,
                                 horizons = PREDICTION_HORIZONS,
                                 fold_id = 1L, replicate_id = 1L,
                                 n_train_setting = NA,
                                 time_grid_setting = NA) {
  x <- observed_profile_covariates(test_subjects, test_visits)
  x <- align_covariate_columns(x, fit$columns)
  pred <- predict(fit$fit, data = x)
  if (is.null(pred$survival) || is.null(pred$unique.death.times)) {
    stop("ranger survival prediction did not return survival curves.", call. = FALSE)
  }
  risk <- t(apply(pred$survival, 1, function(s) {
    1 - step_eval_survival(pred$unique.death.times, s, horizons)
  }))
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

fit_deepsurv_observed <- function(train_subjects, train_visits, seed = NULL) {
  if (!requireNamespace("survivalmodels", quietly = TRUE)) {
    stop("DeepSurv-Observed requires the survivalmodels package, which is not installed.", call. = FALSE)
  }
  x <- observed_profile_covariates(train_subjects, train_visits)
  dat <- data.frame(
    time = train_subjects$U,
    status = train_subjects$delta,
    as.data.frame(x, check.names = FALSE),
    check.names = FALSE
  )
  if (!is.null(seed) && "set_seed" %in% getNamespaceExports("survivalmodels")) {
    survivalmodels::set_seed(seed)
  }
  fit <- survivalmodels::deepsurv(
    survival::Surv(time, status) ~ .,
    data = dat,
    epochs = DEEPSURV_EPOCHS,
    batch_size = DEEPSURV_BATCH_SIZE,
    verbose = FALSE
  )
  list(
    method = "DeepSurv-Observed",
    fit = fit,
    columns = colnames(x),
    times = sort(unique(train_subjects$U[train_subjects$delta == 1])),
    implementation = "survivalmodels::deepsurv on observed landmark-age profile"
  )
}

predict_deepsurv_observed <- function(fit, test_subjects, test_visits,
                                      horizons = PREDICTION_HORIZONS,
                                      fold_id = 1L, replicate_id = 1L,
                                      n_train_setting = NA,
                                      time_grid_setting = NA) {
  x <- observed_profile_covariates(test_subjects, test_visits)
  x <- align_covariate_columns(x, fit$columns)
  surv <- predict(fit$fit, newdata = x, type = "survival")
  surv <- as.matrix(surv)
  if (ncol(surv) != length(fit$times)) {
    stop("DeepSurv survival prediction length did not match the training event-time grid.", call. = FALSE)
  }
  risk <- t(apply(surv, 1, function(s) 1 - step_eval_survival(fit$times, s, horizons)))
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

fit_bayesian_dynamic_observed <- function(train_subjects, train_visits, seed = NULL) {
  if (!isTRUE(BAYES_DYNAMIC_ENABLE_STATIC_FALLBACK)) {
    stop(
      "Bayesian-Dynamic-Observed requires an accumulating cohort-era stream and posterior-to-prior updating; ",
      "the current DGM has no calendar-era process. Set SILK_BAYES_DYNAMIC_STATIC_FALLBACK=true ",
      "only to run a static rstanarm sensitivity model, not the dynamic-updating comparator.",
      call. = FALSE
    )
  }
  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    stop("Bayesian-Dynamic-Observed requires the rstanarm package, which is not installed.", call. = FALSE)
  }
  x <- observed_profile_covariates(train_subjects, train_visits)
  dat <- data.frame(
    time = train_subjects$U,
    status = train_subjects$delta,
    as.data.frame(x, check.names = FALSE),
    check.names = FALSE
  )
  fit <- rstanarm::stan_surv(
    survival::Surv(time, status) ~ .,
    data = dat,
    basehaz = BAYES_SURV_BASEHAZ,
    chains = BAYES_SURV_CHAINS,
    iter = BAYES_SURV_ITER,
    seed = seed,
    refresh = 0
  )
  list(
    method = "Bayesian-Dynamic-Observed",
    fit = fit,
    columns = colnames(x),
    implementation = "rstanarm::stan_surv Bayesian observed-age survival model"
  )
}

predict_bayesian_dynamic_observed <- function(fit, test_subjects, test_visits,
                                              horizons = PREDICTION_HORIZONS,
                                              fold_id = 1L, replicate_id = 1L,
                                              n_train_setting = NA,
                                              time_grid_setting = NA) {
  x <- observed_profile_covariates(test_subjects, test_visits)
  x <- align_covariate_columns(x, fit$columns)
  sf <- rstanarm::posterior_survfit(fit$fit, newdata = x, times = horizons)
  sf <- as.data.frame(sf)
  surv_col <- intersect(c("survpred", "surv", "survival"), names(sf))
  time_col <- intersect(c("time", "times"), names(sf))
  if (!length(surv_col) || !length(time_col)) {
    stop("rstanarm posterior_survfit output did not expose survival and time columns.", call. = FALSE)
  }
  sf$.row_order <- seq_len(nrow(sf))
  sf <- sf[order(sf$.row_order), , drop = FALSE]
  surv <- as.numeric(sf[[surv_col[1]]])
  if (length(surv) != nrow(test_subjects) * length(horizons)) {
    stop("rstanarm posterior_survfit output size did not match subjects times horizons.", call. = FALSE)
  }
  surv_mat <- matrix(surv, nrow = nrow(test_subjects), ncol = length(horizons), byrow = TRUE)
  prediction_frame(
    test_subjects,
    landmark = test_subjects$A_obs,
    horizons = horizons,
    risk_mat = 1 - surv_mat,
    context = method_context(
      fit$method, fold_id, replicate_id,
      n_train_setting, time_grid_setting
    )
  )
}
