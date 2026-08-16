# =============================================================================
# 01_analysis.R
# Final fixed-landmark CV analysis of the MACS PDS
# =============================================================================

macs_current_script_dir <- function(default = getwd()) {
  frames <- sys.frames()
  for (ii in rev(seq_along(frames))) {
    ofile <- frames[[ii]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      ofile_norm <- tryCatch(
        normalizePath(ofile, winslash = "/", mustWork = TRUE),
        error = function(e) normalizePath(basename(ofile), winslash = "/", mustWork = TRUE)
      )
      return(dirname(ofile_norm))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])),
                         winslash = "/", mustWork = TRUE))
  }
  normalizePath(default, winslash = "/", mustWork = TRUE)
}

source(file.path(macs_current_script_dir(), "00_setup.R"), chdir = TRUE)

# ── 1. Load MACS data ───────────────────────────────────────────────────────
subjects_raw <- read.csv(file.path(DATA_DIR, "macs_subjects.csv"))
visits_raw   <- read.csv(file.path(DATA_DIR, "macs_visits.csv"))

# Map to SILK column names
subjects <- macs_to_silk_subjects(subjects_raw)
visits   <- macs_to_silk_visits(visits_raw)

# The final primary analysis excludes the historically imputed viral-load
# series. The retained immunologic markers are complete and therefore require
# no imputation before a landmark.
non_biomarker_columns <- setdiff(names(visits), bio_columns(visits))
missing_primary <- setdiff(MACS_PRIMARY_BIOMARKERS, names(visits))
if (length(missing_primary)) {
  stop("Missing primary biomarkers: ", paste(missing_primary, collapse = ", "),
       call. = FALSE)
}
visits <- visits[, c(non_biomarker_columns, MACS_PRIMARY_BIOMARKERS), drop = FALSE]

# MACS is an application-level configuration. The SILK package itself does
# not assume these covariates; this analysis declares them explicitly.
MACS_SILK_DATA_SPEC <- silk_data_spec(
  biomarker_cols = MACS_PRIMARY_BIOMARKERS,
  covariate_cols = c("X1", "X2"),
  template_input_covariates = c("X1", "X2", "lag"),
  clock_covariates = c("X1", "X2"),
  anchor_covariates = "X1",
  anchor_mode = "rx"
)

cat("MACS PDS data loaded:\n")
cat("  Subjects:", nrow(subjects), "\n")
cat("  Visits:  ", nrow(visits), "\n")
cat("  Events:  ", sum(subjects$delta),
    "(", round(100 * mean(subjects$delta), 1), "%)\n")
cat("  Biomarkers:", paste(bio_columns(visits), collapse = ", "), "\n")

# ── 2. Prediction frame for real data ────────────────────────────────────────
# The package's prediction_frame() uses D_star (simulation only).
# This version handles real censored outcomes properly.

real_prediction_frame <- function(test_subjects, landmark, horizons,
                                  risk_mat, method_name, fold_id = 1L,
                                  n_train_setting = NA_integer_,
                                  landmark_set = "last_visit",
                                  landmark_time = NA_real_) {
  n <- nrow(test_subjects)
  H <- length(horizons)
  risk_mat <- as.matrix(risk_mat)
  if (!identical(dim(risk_mat), c(n, H))) {
    stop("Risk predictions have the wrong subject-by-horizon dimensions.",
         call. = FALSE)
  }
  if (any(!is.finite(risk_mat))) {
    stop("Risk prediction returned missing or non-finite values.", call. = FALSE)
  }
  rows <- list()
  rr <- 1L

  for (i in seq_len(n)) {
    U_i     <- test_subjects$U[i]
    delta_i <- test_subjects$delta[i]
    for (h in seq_len(H)) {
      tau <- horizons[h]
      if (delta_i == 1 && U_i <= tau) {
        event_wh <- 1L; at_risk <- TRUE
      } else if (delta_i == 0 && U_i <= tau) {
        event_wh <- 0L; at_risk <- FALSE   # censored before horizon
      } else {
        event_wh <- 0L; at_risk <- TRUE
      }
      rows[[rr]] <- data.frame(
        subject_id = test_subjects$id[i], landmark = landmark[i],
        horizon = tau, method = method_name,
        risk_pred = clip_probability(risk_mat[i, h]),
        event_within_horizon = event_wh, at_risk = at_risk,
        residual_time = U_i, event_status = delta_i,
        landmark_set = landmark_set, landmark_time = landmark_time,
        fold_id = fold_id, replicate_id = 1L,
        n_train_setting = n_train_setting,
        time_grid_setting = paste0("macs_cv:", landmark_set), stringsAsFactors = FALSE
      )
      rr <- rr + 1L
    }
  }
  do.call(rbind, rows)
}

# ── 3. Cross-validation ─────────────────────────────────────────────────────
set.seed(2024)
K           <- as.integer(Sys.getenv("MACS_OUTER_FOLDS", unset = "5"))
if (!is.finite(K) || K < 2L) stop("MACS_OUTER_FOLDS must be at least 2.", call. = FALSE)
HORIZONS    <- silk_opt("PREDICTION_HORIZONS")
METHOD_LIST <- silk_opt("METHOD_ORDER")
method_override <- trimws(Sys.getenv("MACS_METHODS", unset = ""))
if (nzchar(method_override)) {
  METHOD_LIST <- trimws(strsplit(method_override, ",", fixed = TRUE)[[1]])
}
unknown_methods <- setdiff(METHOD_LIST, silk_opt("METHOD_ORDER"))
if (length(unknown_methods)) {
  stop("Unknown MACS_METHODS values: ", paste(unknown_methods, collapse = ", "),
       call. = FALSE)
}
development_override <- nzchar(method_override) ||
  nzchar(trimws(Sys.getenv("MACS_LANDMARKS", unset = ""))) ||
  nzchar(trimws(Sys.getenv("MACS_INNER_FOLDS", unset = ""))) ||
  !identical(K, 5L) || !isTRUE(MACS_USE_FIXED_LANDMARKS)
if (development_override &&
    !nzchar(trimws(Sys.getenv("MACS_RESULTS_DIR", unset = "")))) {
  stop(
    "Development overrides require an explicit MACS_RESULTS_DIR so that ",
    "partial runs cannot replace the final MACS results.",
    call. = FALSE
  )
}
landmark_sets <- macs_landmark_sets(subjects, visits)

make_stratified_folds <- function(subjects, number_folds, seed = 2024L) {
  set.seed(seed)
  assignment <- integer(nrow(subjects))
  for (status in sort(unique(subjects$delta))) {
    index <- which(subjects$delta == status)
    index <- sample(index, length(index))
    assignment[index] <- rep(seq_len(number_folds), length.out = length(index))
  }
  stats::setNames(assignment, subjects$id)
}

# One subject keeps the same outer fold at every landmark. This preserves the
# pairing between methods and makes landmark-to-landmark summaries auditable.
master_folds <- make_stratified_folds(subjects, K, seed = 2024L)

macs_current_marker_matrix <- function(subjects, visits, biomarkers = MACS_PRIMARY_BIOMARKERS) {
  output <- matrix(NA_real_, nrow = nrow(subjects), ncol = length(biomarkers),
                   dimnames = list(NULL, paste0("current_", biomarkers)))
  for (index in seq_len(nrow(subjects))) {
    subject_visits <- visits[visits$id == subjects$id[index] &
                               visits$A_obs_il <= subjects$A_obs[index] + 1e-8,
                             , drop = FALSE]
    if (!nrow(subject_visits)) next
    subject_visits <- subject_visits[order(subject_visits$A_obs_il, decreasing = TRUE),
                                     , drop = FALSE]
    output[index, ] <- as.numeric(subject_visits[1, biomarkers, drop = TRUE])
  }
  output[!is.finite(output)] <- 0
  output
}

macs_survival_covariates <- function(subjects, visits, stage,
                                     marker_matrix = NULL) {
  if (is.null(marker_matrix)) {
    marker_matrix <- macs_current_marker_matrix(subjects, visits)
  }
  cbind(
    landmark_age = as.numeric(stage),
    X1 = as.numeric(subjects$X1),
    X2 = as.numeric(subjects$X2),
    marker_matrix
  )
}

macs_single_marker_visits <- function(visits, biomarker) {
  non_biomarkers <- setdiff(names(visits), bio_columns(visits))
  output <- visits[, c(non_biomarkers, biomarker), drop = FALSE]
  names(output)[names(output) == biomarker] <- "B1"
  output
}

# Comparator implementation belongs to this application, not to the SILK
# package.  MACS has two baseline covariates, so its working MMLM is explicit
# about that application-level specification rather than silently requesting
# the simulation model with X3/X4.
macs_fit_current_value_mixed_model <- function(subjects, visits) {
  dat <- visits
  dat$marker_value <- dat$B1
  dat$id <- factor(dat$id)
  fit <- NULL
  if (requireNamespace("nlme", quietly = TRUE)) {
    fit <- tryCatch(
      suppressWarnings(
        nlme::lme(
          marker_value ~ A_obs_il + X1 + X2,
          random = ~ A_obs_il | id,
          data = dat,
          control = nlme::lmeControl(
            maxIter = 100, msMaxIter = 100, niterEM = 50,
            returnObject = TRUE
          )
        )
      ),
      error = function(e) NULL
    )
  }
  fixed_effects <- random_cov <- NULL
  residual_sd <- NA_real_
  if (!is.null(fit)) {
    fixed_effects <- as.numeric(nlme::fixef(fit))
    names(fixed_effects) <- names(nlme::fixef(fit))
    random_cov <- tryCatch(
      as.matrix(nlme::getVarCov(fit, type = "random.effects")),
      error = function(e) NULL
    )
    if (!is.null(random_cov) &&
        (is.null(colnames(random_cov)) || any(!nzchar(colnames(random_cov))))) {
      rn <- c("(Intercept)", "A_obs_il")[seq_len(ncol(random_cov))]
      rownames(random_cov) <- rn
      colnames(random_cov) <- rn
    }
    residual_sd <- tryCatch(as.numeric(fit$sigma), error = function(e) NA_real_)
  }
  list(fit = fit, biomarker = "B1", fixed_effects = fixed_effects,
       random_cov = random_cov, residual_sd = residual_sd)
}

macs_predict_current_value_mixed_model <- function(model, subjects, visits) {
  locf <- SILK:::current_marker_value(subjects, visits, "recorded")
  if (is.null(model$fit) || is.null(model$fixed_effects) ||
      is.null(model$random_cov)) return(locf)
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
  pred <- locf
  fixed_design <- function(a_obs, data) {
    value <- function(name) if (name %in% names(data)) {
      as.numeric(data[[name]])
    } else rep(0, length(a_obs))
    X <- cbind(
      "(Intercept)" = 1, A_obs_il = as.numeric(a_obs),
      X1 = value("X1"), X2 = value("X2")
    )
    X[, names(beta), drop = FALSE]
  }
  random_design <- function(a_obs) {
    Z <- cbind("(Intercept)" = 1, A_obs_il = as.numeric(a_obs))
    Z[, z_names, drop = FALSE]
  }
  for (i in seq_len(nrow(subjects))) {
    sv <- visits[visits$id == subjects$id[i], , drop = FALSE]
    sv <- sv[sv$A_obs_il <= subjects$A_obs[i] + 1e-8, , drop = FALSE]
    if (!nrow(sv)) sv <- visits[visits$id == subjects$id[i], , drop = FALSE]
    if (!nrow(sv)) next
    y <- as.numeric(sv$B1)
    X <- fixed_design(sv$A_obs_il, sv)
    Z <- random_design(sv$A_obs_il)
    keep <- is.finite(y) & apply(X, 1, function(row) all(is.finite(row))) &
      apply(Z, 1, function(row) all(is.finite(row)))
    if (!any(keep)) next
    y <- y[keep]; X <- X[keep, , drop = FALSE]; Z <- Z[keep, , drop = FALSE]
    resid <- y - as.numeric(X %*% beta)
    V <- Z %*% D %*% t(Z) + diag(sigma2, nrow(Z))
    bhat <- tryCatch(as.numeric(D %*% t(Z) %*% solve(V, resid)),
                     error = function(e) rep(0, ncol(D)))
    value <- as.numeric(
      fixed_design(subjects$A_obs[i], subjects[i, , drop = FALSE]) %*% beta +
        random_design(subjects$A_obs[i]) %*% bhat
    )
    if (is.finite(value)) pred[i] <- value
  }
  pred[!is.finite(pred)] <- locf[!is.finite(pred)]
  pred
}

macs_fit_residual_cox <- function(subjects, x) {
  x <- as.matrix(x)
  keep <- if (ncol(x)) apply(x, 2, stats::sd, na.rm = TRUE) > 1e-10 else logical(0)
  x <- x[, keep, drop = FALSE]
  km_fit <- survival::survfit(survival::Surv(subjects$U, subjects$delta) ~ 1)
  km <- list(time = km_fit$time, surv = km_fit$surv)
  if (!ncol(x)) return(list(type = "km", km = km, keep = keep))
  center <- colMeans(x); scale <- apply(x, 2, stats::sd)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  x_std <- sweep(sweep(x, 2, center, "-"), 2, scale, "/")
  dat <- data.frame(time = subjects$U, status = subjects$delta, x_std)
  names(dat)[-(1:2)] <- paste0("x", seq_len(ncol(x_std)))
  fit <- tryCatch(
    survival::coxph(
      stats::as.formula(paste("survival::Surv(time, status) ~",
                              paste(names(dat)[-(1:2)], collapse = "+"))),
      data = dat, ties = "breslow", x = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || any(!is.finite(stats::coef(fit)))) {
    return(list(type = "km", km = km, keep = keep))
  }
  list(type = "cox_residual", fit = fit,
       basehaz = survival::basehaz(fit, centered = FALSE), keep = keep,
       center = center, scale = scale, km = km)
}

macs_predict_residual_cox_risk <- function(model, x_new, horizons) {
  x_new <- as.matrix(x_new)
  if (length(model$keep)) x_new <- x_new[, model$keep, drop = FALSE]
  if (identical(model$type, "km") || !ncol(x_new)) {
    index <- findInterval(horizons, model$km$time)
    surv <- vapply(index, function(ii) {
      if (ii < 1L) 1 else model$km$surv[ii]
    }, numeric(1))
    return(matrix(rep(1 - surv, nrow(x_new)), nrow = nrow(x_new), byrow = TRUE))
  }
  x_new <- sweep(sweep(x_new, 2, model$center, "-"), 2, model$scale, "/")
  x_new[!is.finite(x_new)] <- 0
  H0 <- approx(model$basehaz$time, model$basehaz$hazard, horizons,
               method = "constant", rule = 2, f = 0, yleft = 0,
               yright = max(model$basehaz$hazard), ties = "ordered")$y
  lp <- as.numeric(x_new %*% stats::coef(model$fit))
  clip_probability(1 - exp(-outer(exp(lp), H0, "*")))
}

fit_multimarker_mmlm <- function(train_s, train_v) {
  models <- lapply(MACS_PRIMARY_BIOMARKERS, function(biomarker) {
    marker_visits <- macs_single_marker_visits(train_v, biomarker)
    list(
      biomarker = biomarker,
      visits_name = "B1",
      model = macs_fit_current_value_mixed_model(train_s, marker_visits)
    )
  })
  marker_matrix <- do.call(cbind, lapply(models, function(component) {
    marker_visits <- macs_single_marker_visits(train_v, component$biomarker)
    macs_predict_current_value_mixed_model(component$model, train_s, marker_visits)
  }))
  colnames(marker_matrix) <- paste0("mmlm_", MACS_PRIMARY_BIOMARKERS)
  x <- macs_survival_covariates(train_s, train_v, train_s$A_obs, marker_matrix)
  list(method = "MMLM-Recorded", marker_models = models,
       fit = macs_fit_residual_cox(train_s, x))
}

predict_multimarker_mmlm <- function(fit, test_s, test_v) {
  marker_matrix <- do.call(cbind, lapply(fit$marker_models, function(component) {
    marker_visits <- macs_single_marker_visits(test_v, component$biomarker)
    macs_predict_current_value_mixed_model(component$model, test_s, marker_visits)
  }))
  colnames(marker_matrix) <- paste0("mmlm_", MACS_PRIMARY_BIOMARKERS)
  macs_survival_covariates(test_s, test_v, test_s$A_obs, marker_matrix)
}

fit_macs_method <- function(method, train_s, train_v, fold_seed) {
  switch(method,
    "Cox-Recorded-SameFeature" = {
      x <- macs_survival_covariates(train_s, train_v, train_s$A_obs)
      list(method = method, fit = macs_fit_residual_cox(train_s, x))
    },
    "MMLM-Recorded" = fit_multimarker_mmlm(train_s, train_v),
    "SILK-Cox" = {
      registration <- fit_silk_registration(
        train_s, train_v,
        shift_range = MACS_SHIFT_RANGE,
        seed = fold_seed,
        biomarker_kernel = "gaussian",
        data_spec = MACS_SILK_DATA_SPEC
      )
      stage <- registration$train_stage$S_hat[
        match(train_s$id, registration$train_stage$id)
      ]
      x <- macs_survival_covariates(train_s, train_v, stage)
      list(method = method, registration = registration,
           grid = registration$grid, fit = macs_fit_residual_cox(train_s, x))
    },
    stop("Unknown MACS method: ", method, call. = FALSE)
  )
}

silk_diagnostics <- function(fit, ps, fold_id, method, landmark_set, landmark_time) {
  train_diag <- fit$registration$train_stage
  train_diag$split <- "train_crossfit"
  test_diag <- ps
  test_diag$S_hat <- NA_real_
  test_diag$fold <- fold_id
  test_diag$split <- "test"
  common <- intersect(names(train_diag), names(test_diag))
  out <- rbind(train_diag[, common, drop = FALSE],
               test_diag[, common, drop = FALSE])
  out$fold_id <- fold_id
  out$method <- method
  out$landmark_set <- landmark_set
  out$landmark_time <- landmark_time
  out
}

predict_macs_method <- function(method, fit, test_s, test_v, fold_id,
                                n_train, landmark_set, landmark_time) {
  if (identical(method, "Cox-Recorded-SameFeature")) {
    x <- macs_survival_covariates(test_s, test_v, test_s$A_obs)
    risk <- macs_predict_residual_cox_risk(fit$fit, x, HORIZONS)
    return(list(
      predictions = real_prediction_frame(test_s, test_s$A_obs, HORIZONS, risk,
                                          method, fold_id, n_train,
                                          landmark_set, landmark_time),
      diagnostics = NULL
    ))
  }
  if (identical(method, "MMLM-Recorded")) {
    x <- predict_multimarker_mmlm(fit, test_s, test_v)
    risk <- macs_predict_residual_cox_risk(fit$fit, x, HORIZONS)
    return(list(
      predictions = real_prediction_frame(test_s, test_s$A_obs, HORIZONS, risk,
                                          method, fold_id, n_train,
                                          landmark_set, landmark_time),
      diagnostics = NULL
    ))
  }
  grid <- fit$grid
  ps <- predict_silk_registration(fit$registration, test_v, grid)
  stage <- test_s$A_obs[match(ps$id, test_s$id)] - ps$e_hat
  stage <- stage[match(test_s$id, ps$id)]
  ps$S_hat <- test_s$A_obs[match(ps$id, test_s$id)] - ps$e_hat
  x <- macs_survival_covariates(test_s, test_v, stage)
  risk <- macs_predict_residual_cox_risk(fit$fit, x, HORIZONS)
  list(
    predictions = real_prediction_frame(test_s, stage, HORIZONS, risk,
                                        method, fold_id, n_train,
                                        landmark_set, landmark_time),
    diagnostics = silk_diagnostics(fit, ps, fold_id, method, landmark_set, landmark_time)
  )
}

all_predictions <- list()
all_status      <- list()
all_diagnostics <- list()
counter         <- 0L
diag_counter    <- 0L

cat("\n=== ", K, "-fold cross-validation ===\n")
cat("Methods:", paste(METHOD_LIST, collapse = ", "), "\n")
cat("Landmarks:", paste(names(landmark_sets), collapse = ", "), "\n")
cat("Horizons:", paste(HORIZONS, collapse = ", "), "years\n\n")
cat("SILK biomarker kernel:", silk_opt("BIOMARKER_KERNEL"),
    "| bandwidth:", silk_opt("BIOMARKER_BANDWIDTH"), "\n")
cat("Primary biomarkers:", paste(MACS_PRIMARY_BIOMARKERS, collapse = ", "),
    "(viral load excluded from the primary analysis)\n\n")

for (landmark_set in names(landmark_sets)) {
  ldat <- landmark_sets[[landmark_set]]
  set_subjects <- ldat$subjects
  set_visits <- ldat$visits
  landmark_time <- ldat$landmark_time
  if (nrow(set_subjects) < K || length(unique(set_subjects$delta)) < 2L) {
    cat("Skipping", landmark_set, "because the risk set is too small or has one outcome class\n")
    next
  }
  ids <- sort(unique(set_subjects$id))
  folds <- master_folds[as.character(ids)]
  if (anyNA(folds)) stop("Outer-fold assignment is missing for a landmark subject.")
  cat("\n── Landmark set:", landmark_set, "subjects:", length(ids), "visits:", nrow(set_visits), "──\n")

  for (ff in seq_len(K)) {
    cat("Fold", ff, "/", K, "\n")
    train_ids <- ids[folds != ff]
    test_ids  <- ids[folds == ff]

    train_s <- set_subjects[set_subjects$id %in% train_ids, ]
    train_v <- set_visits[set_visits$id %in% train_ids, ]
    test_s  <- set_subjects[set_subjects$id %in% test_ids, ]
    test_v  <- set_visits[set_visits$id %in% test_ids, ]

    cat("  Train:", nrow(train_s), "subj,", nrow(train_v), "visits |",
        "Test:", nrow(test_s), "subj,", nrow(test_v), "visits\n")

    for (method in METHOD_LIST) {
      cat("  ", method, "... ")
      t0 <- proc.time()[3]
      fit_res <- tryCatch(
        list(ok = TRUE, value = fit_macs_method(method, train_s, train_v, 42L + ff)),
        error = function(e) list(ok = FALSE, error = conditionMessage(e))
      )
      fit_time <- proc.time()[3] - t0
      status_key <- paste(landmark_set, ff, method, sep = "_")

      if (!fit_res$ok) {
        cat("FIT FAILED (", fit_res$error, ")\n")
        all_status[[status_key]] <- data.frame(
          landmark_set = landmark_set, landmark_time = landmark_time,
          fold = ff, method = method, fit_ok = FALSE, predict_ok = FALSE,
          fit_time = fit_time, predict_time = NA, error = fit_res$error,
          stringsAsFactors = FALSE)
        next
      }

      pred_t0 <- proc.time()[3]
      pred_res <- tryCatch(
        list(ok = TRUE, value = predict_macs_method(
          method, fit_res$value, test_s, test_v, ff, nrow(train_s),
          landmark_set, landmark_time
        )),
        error = function(e) list(ok = FALSE, error = conditionMessage(e))
      )
      predict_time <- proc.time()[3] - pred_t0
      cat(sprintf("%.1fs fit, %.1fs predict", fit_time, predict_time))

      if (pred_res$ok) {
        cat(" done\n")
        counter <- counter + 1L
        all_predictions[[counter]] <- pred_res$value$predictions
        if (!is.null(pred_res$value$diagnostics)) {
          diag_counter <- diag_counter + 1L
          all_diagnostics[[diag_counter]] <- pred_res$value$diagnostics
        }
      } else {
        cat(" PREDICT FAILED (", pred_res$error, ")\n")
      }

      all_status[[status_key]] <- data.frame(
        landmark_set = landmark_set, landmark_time = landmark_time,
        fold = ff, method = method, fit_ok = TRUE,
        predict_ok = pred_res$ok, fit_time = fit_time,
        predict_time = predict_time,
        error = if (!pred_res$ok) pred_res$error else "",
        stringsAsFactors = FALSE)
    }
    cat("\n")
  }
}

# ── 4. Save ──────────────────────────────────────────────────────────────────
if (!length(all_predictions)) {
  stop("No MACS predictions were produced; inspect macs_method_status.csv")
}

predictions <- do.call(rbind, all_predictions)
rownames(predictions) <- NULL
status_df <- do.call(rbind, all_status)
rownames(status_df) <- NULL
diagnostics <- if (length(all_diagnostics)) do.call(rbind, all_diagnostics) else data.frame()
if (nrow(diagnostics)) rownames(diagnostics) <- NULL

write.csv(status_df, file.path(RESULTS_DIR, "macs_method_status.csv"),
          row.names = FALSE)
if (any(!status_df$fit_ok | !status_df$predict_ok)) {
  stop(
    "The final MACS analysis is incomplete. Inspect macs_method_status.csv; ",
    "macs_predictions.csv was not replaced.",
    call. = FALSE
  )
}

prediction_file <- file.path(RESULTS_DIR, "macs_predictions.csv")
prediction_temporary <- tempfile("macs_predictions_", tmpdir = RESULTS_DIR,
                                 fileext = ".csv")
write.csv(predictions, prediction_temporary, row.names = FALSE)
if (!file.copy(prediction_temporary, prediction_file, overwrite = TRUE)) {
  unlink(prediction_temporary)
  stop("Could not publish the completed MACS prediction file.", call. = FALSE)
}
unlink(prediction_temporary)
write.csv(diagnostics, file.path(RESULTS_DIR, "macs_profile_diagnostics.csv"),
          row.names = FALSE)

analysis_specification <- data.frame(
  item = c("outer_folds", "fold_seed", "landmarks", "horizons", "methods",
           "primary_biomarkers", "viral_load_primary", "shift_range",
           "biomarker_kernel", "biomarker_bandwidth", "shift_ridge"),
  value = c(
    K, 2024, paste(MACS_FIXED_LANDMARKS, collapse = ","),
    paste(HORIZONS, collapse = ","), paste(METHOD_LIST, collapse = ","),
    paste(MACS_PRIMARY_BIOMARKERS, collapse = ","), "excluded",
    paste(MACS_SHIFT_RANGE, collapse = ","), silk_opt("BIOMARKER_KERNEL"),
    silk_opt("BIOMARKER_BANDWIDTH"), silk_opt("SHIFT_RIDGE")
  ),
  stringsAsFactors = FALSE
)
write.csv(analysis_specification,
          file.path(RESULTS_DIR, "macs_analysis_specification.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(RESULTS_DIR, "macs_session_info.txt"))

cat("=== Cross-validation complete ===\n")
cat("Predictions:", nrow(predictions), "rows\n")
cat("Diagnostics:", nrow(diagnostics), "rows\n")
cat("Saved to:   ", RESULTS_DIR, "\n\n")

cat("Method status:\n")
print(table(status_df$method, status_df$predict_ok))
cat("\nAt-risk counts per method/horizon:\n")
print(xtabs(at_risk ~ method + horizon, data = predictions))
