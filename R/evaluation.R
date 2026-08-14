# =============================================================================
# evaluation.R
# Evaluation metrics (from evaluation_prediction.R)
# =============================================================================

#' @keywords internal
PREDICTION_COLUMNS <- c(
  "subject_id", "landmark", "horizon", "method", "risk_pred",
  "event_within_horizon", "at_risk", "fold_id", "replicate_id",
  "n_train_setting", "time_grid_setting"
)

#' Validate a SILK prediction frame
#'
#' @param pred Data frame of subject-by-horizon predictions.
#' @return The validated prediction frame in canonical column order.
#' @export
validate_prediction_frame <- function(pred) {
  missing <- setdiff(PREDICTION_COLUMNS, names(pred))
  if (length(missing)) {
    stop("Prediction frame is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  pred <- pred[, PREDICTION_COLUMNS, drop = FALSE]
  pred$risk_pred <- clip_probability(pred$risk_pred)
  pred$event_within_horizon <- as.integer(pred$event_within_horizon)
  pred$at_risk <- as.logical(pred$at_risk)
  pred
}

#' Compute binary AUC via the Mann-Whitney U statistic
#'
#' @param y Integer vector of binary outcomes (0/1).
#' @param score Numeric vector of predicted scores.
#' @return Numeric scalar AUC, or NA if undefined.
#' @export
binary_auc <- function(y, score) {
  y <- as.integer(y)
  ok <- is.finite(score) & y %in% c(0L, 1L)
  y <- y[ok]
  score <- score[ok]
  # Use doubles so n1 * n0 does not overflow for large validation samples.
  n1 <- as.double(sum(y == 1L))
  n0 <- as.double(sum(y == 0L))
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Compute calibration statistics
#'
#' @param y Integer vector of binary outcomes (0/1).
#' @param p Numeric vector of predicted probabilities.
#' @return Named numeric vector with calibration_intercept, calibration_slope,
#'   and calibration_in_large.
#' @export
calibration_stats <- function(y, p) {
  y <- as.integer(y)
  p <- clip_probability(p)
  ok <- is.finite(p) & y %in% c(0L, 1L)
  y <- y[ok]
  p <- p[ok]
  if (length(unique(y)) < 2L || length(p) < 10L) {
    return(c(calibration_intercept = NA_real_, calibration_slope = NA_real_,
             calibration_in_large = mean(y) - mean(p)))
  }
  lp <- stats::qlogis(p)
  intercept_fit <- tryCatch(
    suppressWarnings(stats::glm(y ~ 1, offset = lp, family = stats::binomial())),
    error = function(e) NULL
  )
  slope_fit <- tryCatch(
    suppressWarnings(stats::glm(y ~ lp, family = stats::binomial())),
    error = function(e) NULL
  )
  c(
    calibration_intercept = if (is.null(intercept_fit)) NA_real_ else unname(stats::coef(intercept_fit)[1]),
    calibration_slope = if (is.null(slope_fit)) NA_real_ else unname(stats::coef(slope_fit)[2]),
    calibration_in_large = mean(y) - mean(p)
  )
}

#' Compute per-horizon prediction metrics
#'
#' @param pred A prediction frame data frame.
#' @return Data frame with one row per method/horizon/replicate combination.
#' @export
prediction_per_horizon_metrics <- function(pred) {
  pred <- validate_prediction_frame(pred)
  pred <- pred[pred$at_risk %in% TRUE, , drop = FALSE]
  key <- interaction(
    pred$time_grid_setting, pred$n_train_setting, pred$fold_id,
    pred$replicate_id, pred$method, pred$horizon,
    sep = "||", drop = TRUE
  )

  rows <- lapply(split(pred, key), function(z) {
    y <- z$event_within_horizon
    p <- z$risk_pred
    cal <- calibration_stats(y, p)
    data.frame(
      time_grid_setting = z$time_grid_setting[1],
      n_train_setting = z$n_train_setting[1],
      fold_id = z$fold_id[1],
      replicate_id = z$replicate_id[1],
      method = z$method[1],
      horizon = z$horizon[1],
      n = nrow(z),
      event_rate = mean(y),
      mean_risk_pred = mean(p),
      brier_score = mean((y - p)^2),
      auc = binary_auc(y, p),
      calibration_intercept = cal[["calibration_intercept"]],
      calibration_slope = cal[["calibration_slope"]],
      calibration_in_large = cal[["calibration_in_large"]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Compute replicate-level prediction metrics
#'
#' @param per_horizon Data frame from \code{prediction_per_horizon_metrics}.
#' @return Data frame with one row per method/replicate combination.
#' @export
prediction_replicate_metrics <- function(per_horizon) {
  EVAL_HORIZON <- silk_opt("EVAL_HORIZON")
  key <- interaction(
    per_horizon$time_grid_setting, per_horizon$n_train_setting,
    per_horizon$fold_id, per_horizon$replicate_id, per_horizon$method,
    sep = "||", drop = TRUE
  )
  rows <- lapply(split(per_horizon, key), function(z) {
    z <- z[order(z$horizon), , drop = FALSE]
    hidx <- which.min(abs(z$horizon - EVAL_HORIZON))
    ibs <- if (nrow(z) >= 2L) {
      trapz(z$horizon, z$brier_score) / (max(z$horizon) - min(z$horizon))
    } else {
      z$brier_score[1]
    }
    data.frame(
      time_grid_setting = z$time_grid_setting[1],
      n_train_setting = z$n_train_setting[1],
      fold_id = z$fold_id[1],
      replicate_id = z$replicate_id[1],
      method = z$method[1],
      integrated_brier_score = ibs,
      mean_brier_score = mean(z$brier_score, na.rm = TRUE),
      mean_auc = mean(z$auc, na.rm = TRUE),
      brier_score_horizon = z$brier_score[hidx],
      auc_horizon = z$auc[hidx],
      calibration_intercept = z$calibration_intercept[hidx],
      calibration_slope = z$calibration_slope[hidx],
      calibration_in_large = z$calibration_in_large[hidx],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @keywords internal
prediction_calibration_bins <- function(pred, n_bins = 10L) {
  pred <- validate_prediction_frame(pred)
  pred <- pred[pred$at_risk %in% TRUE, , drop = FALSE]
  key <- interaction(
    pred$time_grid_setting, pred$n_train_setting, pred$fold_id,
    pred$replicate_id, pred$method, pred$horizon,
    sep = "||", drop = TRUE
  )

  rows <- lapply(split(pred, key), function(z) {
    q <- stats::quantile(z$risk_pred, probs = seq(0, 1, length.out = n_bins + 1L),
                         na.rm = TRUE, type = 8)
    q <- unique(q)
    if (length(q) < 3L) return(NULL)
    z$bin <- cut(z$risk_pred, breaks = q, include.lowest = TRUE, labels = FALSE)
    do.call(rbind, lapply(split(z, z$bin), function(b) {
      data.frame(
        time_grid_setting = b$time_grid_setting[1],
        n_train_setting = b$n_train_setting[1],
        fold_id = b$fold_id[1],
        replicate_id = b$replicate_id[1],
        method = b$method[1],
        horizon = b$horizon[1],
        bin = b$bin[1],
        n = nrow(b),
        mean_predicted_risk = mean(b$risk_pred),
        observed_risk = mean(b$event_within_horizon),
        stringsAsFactors = FALSE
      )
    }))
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Summarize calibration bins across simulation replications
#'
#' @param calibration_bins Data frame returned in the calibration component of
#'   \code{evaluate_prediction_frame}.
#' @return A data frame of weighted calibration-bin summaries.
#' @export
summarize_calibration_bins <- function(calibration_bins) {
  if (is.null(calibration_bins) || !nrow(calibration_bins)) return(data.frame())
  key <- interaction(
    calibration_bins$time_grid_setting, calibration_bins$n_train_setting,
    calibration_bins$method, calibration_bins$horizon, calibration_bins$bin,
    sep = "||", drop = TRUE
  )
  rows <- lapply(split(calibration_bins, key), function(z) {
    w <- z$n / sum(z$n)
    data.frame(
      time_grid_setting = z$time_grid_setting[1],
      n_train_setting = z$n_train_setting[1],
      method = z$method[1],
      horizon = z$horizon[1],
      bin = z$bin[1],
      n = sum(z$n),
      mean_predicted_risk = sum(w * z$mean_predicted_risk),
      observed_risk = sum(w * z$observed_risk),
      n_replicates = length(unique(z$replicate_id)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @keywords internal
mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

#' @keywords internal
sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

#' Summarize prediction metrics across replicates
#'
#' @param replicate_metrics Data frame from \code{prediction_replicate_metrics}.
#' @return Data frame with mean, SD, SE, and CI for each metric by method.
#' @export
prediction_method_summary <- function(replicate_metrics) {
  METHOD_ORDER <- silk_opt("METHOD_ORDER")
  metric_names <- c(
    "integrated_brier_score", "mean_brier_score", "mean_auc",
    "brier_score_horizon", "auc_horizon", "calibration_intercept",
    "calibration_slope", "calibration_in_large"
  )
  key <- interaction(
    replicate_metrics$time_grid_setting, replicate_metrics$n_train_setting,
    replicate_metrics$method,
    sep = "||", drop = TRUE
  )
  rows <- lapply(split(replicate_metrics, key), function(z) {
    out <- data.frame(
      time_grid_setting = z$time_grid_setting[1],
      n_train_setting = z$n_train_setting[1],
      method = z$method[1],
      n_replicates = length(unique(z$replicate_id)),
      stringsAsFactors = FALSE
    )
    for (nm in metric_names) {
      vals <- z[[nm]]
      mu <- mean_or_na(vals)
      sd <- sd_or_na(vals)
      n_eff <- sum(is.finite(vals))
      se <- if (is.finite(sd) && n_eff > 1L) sd / sqrt(n_eff) else NA_real_
      out[[paste0(nm, "_mean")]] <- mu
      out[[paste0(nm, "_sd")]] <- sd
      out[[paste0(nm, "_se")]] <- se
      out[[paste0(nm, "_ci_low")]] <- if (is.finite(se)) mu - 1.96 * se else NA_real_
      out[[paste0(nm, "_ci_high")]] <- if (is.finite(se)) mu + 1.96 * se else NA_real_
    }
    out
  })
  out <- do.call(rbind, rows)
  out$method <- factor(out$method, levels = METHOD_ORDER)
  out <- out[order(out$time_grid_setting, out$n_train_setting, out$method), ]
  rownames(out) <- NULL
  out
}

#' Compute paired differences between SILK and comparators
#'
#' @param replicate_metrics Data frame from \code{prediction_replicate_metrics}.
#' @param reference Character. Reference method name (default "SILK-Gaussian").
#' @return Data frame of paired differences with CIs.
#' @export
prediction_paired_differences <- function(replicate_metrics, reference = "SILK-Gaussian") {
  metric_names <- c(
    "integrated_brier_score", "mean_brier_score", "mean_auc",
    "brier_score_horizon", "auc_horizon", "calibration_in_large"
  )
  setting_key <- interaction(
    replicate_metrics$time_grid_setting, replicate_metrics$n_train_setting,
    replicate_metrics$fold_id, sep = "||", drop = TRUE
  )
  rows <- list()
  rr <- 1L
  for (setting in split(replicate_metrics, setting_key)) {
    ref <- setting[setting$method == reference, , drop = FALSE]
    comparators <- setdiff(unique(setting$method), reference)
    for (comp in comparators) {
      zz <- setting[setting$method == comp, , drop = FALSE]
      for (met in metric_names) {
        a <- ref[, c("replicate_id", met), drop = FALSE]
        b <- zz[, c("replicate_id", met), drop = FALSE]
        names(a)[2] <- "reference"
        names(b)[2] <- "comparator_value"
        pair <- merge(a, b, by = "replicate_id")
        pair$diff <- pair$reference - pair$comparator_value
        pair <- pair[is.finite(pair$diff), , drop = FALSE]
        mu <- mean_or_na(pair$diff)
        sd <- sd_or_na(pair$diff)
        se <- if (is.finite(sd) && nrow(pair) > 1L) sd / sqrt(nrow(pair)) else NA_real_
        rows[[rr]] <- data.frame(
          time_grid_setting = setting$time_grid_setting[1],
          n_train_setting = setting$n_train_setting[1],
          fold_id = setting$fold_id[1],
          reference_method = reference,
          comparator = comp,
          metric = met,
          n_paired = nrow(pair),
          difference_mean = mu,
          difference_sd = sd,
          difference_se = se,
          difference_ci_low = if (is.finite(se)) mu - 1.96 * se else NA_real_,
          difference_ci_high = if (is.finite(se)) mu + 1.96 * se else NA_real_,
          direction = if (grepl("auc", met, fixed = TRUE)) "positive favors SILK" else "negative favors SILK",
          stringsAsFactors = FALSE
        )
        rr <- rr + 1L
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Evaluate a complete prediction frame
#'
#' Computes per-horizon metrics, replicate-level metrics, method summaries,
#' paired differences, and calibration bins from a prediction frame.
#'
#' @param pred A prediction frame data frame.
#' @return A list with elements: per_horizon, replicate_metrics, summary,
#'   paired, calibration_bins.
#' @export
#' @examples
#' \dontrun{
#' dat <- generate_dataset_fixed(200, "mean_moderate", seed = 1)
#' fit <- fit_silk(dat$subjects, dat$visits, shift_range = c(-12, 12), seed = 1)
#' test <- generate_dataset_fixed(100, "mean_moderate", seed = 2)
#' pred <- predict_silk(fit, test$subjects, test$visits)
#' eval_results <- evaluate_prediction_frame(pred)
#' }
evaluate_prediction_frame <- function(pred) {
  per_horizon <- prediction_per_horizon_metrics(pred)
  replicate_metrics <- prediction_replicate_metrics(per_horizon)
  summary <- prediction_method_summary(replicate_metrics)
  paired <- prediction_paired_differences(replicate_metrics)
  calibration_bins <- prediction_calibration_bins(pred)
  list(
    per_horizon = per_horizon,
    replicate_metrics = replicate_metrics,
    summary = summary,
    paired = paired,
    calibration_bins = calibration_bins
  )
}
