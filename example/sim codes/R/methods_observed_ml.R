# =============================================================================
# methods_observed_ml.R
# Observed-clock machine-learning survival comparators.
# These methods use recorded landmark age and the same observed-history summary
# available to the primary pipeline. They do not estimate latent age and hence
# remain practical recorded-clock benchmarks.
# =============================================================================

fit_rsf_observed <- function(train_subjects, train_visits, seed = NULL) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("RSF-Observed requires the ranger package, which is not installed.", call. = FALSE)
  }
  x <- observed_profile_covariates(train_subjects, train_visits)
  keep <- apply(x, 2L, stats::sd, na.rm = TRUE) > 1e-10
  x <- x[, keep, drop = FALSE]
  standardized <- standardize_matrix(x)
  dat <- data.frame(
    time = train_subjects$U,
    status = train_subjects$delta,
    as.data.frame(standardized$train, check.names = FALSE),
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
    center = standardized$center,
    scale = standardized$scale,
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
  x <- apply_standardizer(x, fit)
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
  keep <- apply(x, 2L, stats::sd, na.rm = TRUE) > 1e-10
  x <- x[, keep, drop = FALSE]
  standardized <- standardize_matrix(x)
  dat <- data.frame(
    time = train_subjects$U,
    status = train_subjects$delta,
    as.data.frame(standardized$train, check.names = FALSE),
    check.names = FALSE
  )
  if (!is.null(seed) && "set_seed" %in% getNamespaceExports("survivalmodels")) {
    survivalmodels::set_seed(seed)
  }
  fit <- survivalmodels::deepsurv(
    survival::Surv(time, status) ~ .,
    data = dat,
    frac = DEEPSURV_VALIDATION_FRACTION,
    num_nodes = c(32L, 16L),
    dropout = 0.10,
    early_stopping = TRUE,
    best_weights = TRUE,
    patience = DEEPSURV_PATIENCE,
    epochs = DEEPSURV_EPOCHS,
    batch_size = DEEPSURV_BATCH_SIZE,
    device = "cpu",
    num_workers = 0L,
    verbose = FALSE
  )
  list(
    method = "DeepSurv-Observed",
    fit = fit,
    columns = colnames(x),
    center = standardized$center,
    scale = standardized$scale,
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
  x <- apply_standardizer(x, fit)
  # survivalmodels validates newdata as a data.frame. Matrix-preserving
  # standardization is useful internally, but must not leak across this API
  # boundary (the confirmatory run fitted successfully and then failed here).
  x <- as.data.frame(x, check.names = FALSE)
  surv <- predict(fit$fit, newdata = x, type = "survival")
  surv <- as.matrix(surv)
  if (nrow(surv) != nrow(test_subjects) && ncol(surv) == nrow(test_subjects)) {
    surv <- t(surv)
  }
  if (nrow(surv) != nrow(test_subjects)) {
    stop("DeepSurv survival prediction did not return one row per subject.", call. = FALSE)
  }
  prediction_times <- suppressWarnings(as.numeric(colnames(surv)))
  if (length(prediction_times) != ncol(surv) || any(!is.finite(prediction_times))) {
    prediction_times <- fit$times
  }
  if (ncol(surv) != length(prediction_times)) {
    stop("DeepSurv survival prediction length did not match the training event-time grid.", call. = FALSE)
  }
  risk <- t(apply(surv, 1, function(s) {
    1 - step_eval_survival(prediction_times, s, horizons)
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
