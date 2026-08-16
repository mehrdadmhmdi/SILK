# =============================================================================
# 02_results.R
# Evaluate CV predictions and produce tables + figures
# =============================================================================

macs_results_script_dir <- function(default = getwd()) {
  frames <- sys.frames()
  for (ii in rev(seq_along(frames))) {
    ofile <- frames[[ii]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(dirname(normalizePath(ofile, winslash = "/", mustWork = TRUE)))
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

source(file.path(macs_results_script_dir(), "00_setup.R"), chdir = TRUE)
library(ggplot2)

# Remove exact generated filenames from superseded exploratory reporting. The
# final analysis contains no age-correction diagnostic and no cross-cell mean.
obsolete_outputs <- c(
  file.path(RESULTS_DIR, "macs_silk_shifts.csv"),
  file.path(RESULTS_DIR, "macs_summary.csv"),
  file.path(RESULTS_DIR, "macs_predictions_previous.csv"),
  file.path(FIGURES_DIR, "fig4a_shift_distribution.pdf"),
  file.path(FIGURES_DIR, "fig4a_shift_distribution.png"),
  file.path(FIGURES_DIR, "fig4b_calibrated_landmark.pdf"),
  file.path(FIGURES_DIR, "fig4b_calibrated_landmark.png"),
  file.path(FIGURES_DIR, "fig5_km_risk_strata.pdf"),
  file.path(FIGURES_DIR, "fig5_km_risk_strata.png"),
  # Superseded generic figure names; the canonical outputs now carry the
  # manuscript's macs_* keys.
  file.path(FIGURES_DIR, "fig1_auc_by_horizon.pdf"),
  file.path(FIGURES_DIR, "fig1_auc_by_horizon.png"),
  file.path(FIGURES_DIR, "fig2_brier_by_horizon.pdf"),
  file.path(FIGURES_DIR, "fig2_brier_by_horizon.png"),
  file.path(FIGURES_DIR, "fig3_calibration.pdf"),
  file.path(FIGURES_DIR, "fig3_calibration.png"),
  file.path(FIGURES_DIR, "fig4_km_risk_strata.pdf"),
  file.path(FIGURES_DIR, "fig4_km_risk_strata.png")
)
invisible(file.remove(obsolete_outputs[file.exists(obsolete_outputs)]))

# ── 1. Load predictions ─────────────────────────────────────────────────────
predictions <- read.csv(file.path(RESULTS_DIR, "macs_predictions.csv"),
                        stringsAsFactors = FALSE)
predictions$at_risk <- as.logical(predictions$at_risk)

# Also load subjects for KM analysis
subjects_raw <- read.csv(file.path(DATA_DIR, "macs_subjects.csv"))
subjects <- macs_to_silk_subjects(subjects_raw)

METHOD_LABELS <- silk_opt("METHOD_LABELS")
UIUC_ORANGE   <- silk_opt("UIUC_ORANGE")
UIUC_BLUE     <- silk_opt("UIUC_BLUE")
HORIZONS      <- sort(unique(predictions$horizon))

cat("Loaded", nrow(predictions), "prediction rows\n")
cat("Methods:", paste(unique(predictions$method), collapse = ", "), "\n")

# ── IPCW helpers for horizon-specific real-data metrics ─────────────────────
if (!all(c("residual_time", "event_status") %in% names(predictions))) {
  pred_eval <- merge(predictions, subjects[, c("id", "U", "delta")],
                     by.x = "subject_id", by.y = "id", all.x = TRUE)
  pred_eval$residual_time <- pred_eval$U
  pred_eval$event_status <- pred_eval$delta
} else {
  pred_eval <- predictions
  pred_eval$residual_time <- as.numeric(pred_eval$residual_time)
  pred_eval$event_status <- as.integer(pred_eval$event_status)
  pred_eval$U <- pred_eval$residual_time
  pred_eval$delta <- pred_eval$event_status
}
if (!"landmark_set" %in% names(pred_eval)) pred_eval$landmark_set <- "last_visit"
if (!"landmark_time" %in% names(pred_eval)) pred_eval$landmark_time <- NA_real_

outcome_data <- unique(pred_eval[, c("subject_id", "landmark_set", "landmark_time",
                                     "residual_time", "event_status")])
outcome_data <- outcome_data[is.finite(outcome_data$residual_time), ]

make_censor_survival <- function(data) {
  fit <- survival::survfit(
    survival::Surv(residual_time, 1L - event_status) ~ 1,
    data = data
  )
  function(times, left_limit = FALSE) {
    times <- pmax(as.numeric(times), 0)
    vapply(times, function(time) {
      index <- if (isTRUE(left_limit)) {
        which(fit$time < time)
      } else {
        which(fit$time <= time)
      }
      value <- if (length(index)) fit$surv[max(index)] else 1
      pmax(value, 1e-6)
    }, numeric(1))
  }
}

ipcw_horizon <- function(U, delta, horizon, censor_survival) {
  U <- as.numeric(U)
  delta <- as.integer(delta)
  y <- rep(NA_real_, length(U))
  w <- rep(0, length(U))

  cases <- delta == 1L & U <= horizon
  controls <- U > horizon

  y[cases] <- 1
  w[cases] <- 1 / censor_survival(U[cases], left_limit = TRUE)
  y[controls] <- 0
  w[controls] <- 1 / censor_survival(rep(horizon, sum(controls)))

  data.frame(y = y, w = w, usable = w > 0)
}

weighted_auc <- function(y, p, w) {
  keep <- is.finite(y) & is.finite(p) & is.finite(w) & w > 0
  y <- y[keep]; p <- p[keep]; w <- w[keep]
  if (length(unique(y)) < 2) return(NA_real_)
  case <- y == 1
  ctrl <- y == 0
  pc <- p[case]; p0 <- p[ctrl]
  wc <- w[case]; w0 <- w[ctrl]
  denom <- sum(wc) * sum(w0)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  cmp <- outer(pc, p0, FUN = "-")
  ww <- outer(wc, w0, FUN = "*")
  sum(ww * ((cmp > 0) + 0.5 * (cmp == 0))) / denom
}

weighted_calibration_stats <- function(y, p, w) {
  keep <- is.finite(y) & is.finite(p) & is.finite(w) & w > 0
  y <- y[keep]; p <- clip_probability(p[keep]); w <- w[keep]
  out <- c(calibration_intercept = NA_real_,
           calibration_slope = NA_real_,
           calibration_in_large = NA_real_)
  if (!length(y)) return(out)
  out[["calibration_in_large"]] <- weighted.mean(y, w) - weighted.mean(p, w)
  if (length(unique(y)) < 2) return(out)
  lp <- qlogis(p)
  fit_slope <- tryCatch(
    suppressWarnings(glm(y ~ lp, family = binomial(), weights = w)),
    error = function(e) NULL)
  if (!is.null(fit_slope)) {
    out[["calibration_slope"]] <- unname(coef(fit_slope)[["lp"]])
  }
  fit_intercept <- tryCatch(
    suppressWarnings(glm(y ~ 1, family = binomial(), weights = w,
                         offset = lp)),
    error = function(e) NULL)
  if (!is.null(fit_intercept)) {
    out[["calibration_intercept"]] <- unname(coef(fit_intercept)[[1]])
  }
  out
}

# ── 2. Per-horizon metrics ───────────────────────────────────────────────────
metrics_list <- list()
for (landmark_set in sort(unique(pred_eval$landmark_set))) {
  censor_survival <- make_censor_survival(
    outcome_data[outcome_data$landmark_set == landmark_set, , drop = FALSE]
  )
  for (method in unique(pred_eval$method)) {
    for (h in HORIZONS) {
      z <- pred_eval[pred_eval$landmark_set == landmark_set &
                       pred_eval$method == method &
                       pred_eval$horizon == h, ]
      if (nrow(z) < 10) next
      hw <- ipcw_horizon(z$U, z$delta, h, censor_survival)
      y <- hw$y
      w <- hw$w
      p <- z$risk_pred

      # Per-fold metrics. The Brier score is defined whenever a fold has any
      # usable subject; the AUC additionally needs both outcome classes.
      fold_auc <- fold_brier <- numeric(0)
      for (ff in unique(z$fold_id)) {
        idx <- z$fold_id == ff & hw$usable
        if (sum(idx) > 0) {
          fold_brier <- c(fold_brier, weighted.mean((y[idx] - p[idx])^2, w[idx]))
          if (length(unique(y[idx])) == 2) {
            fold_auc <- c(fold_auc, weighted_auc(y[idx], p[idx], w[idx]))
          }
        }
      }
      cal <- weighted_calibration_stats(y, p, w)
      usable <- hw$usable

      metrics_list[[length(metrics_list) + 1]] <- data.frame(
        landmark_set = landmark_set,
        landmark_time = unique(z$landmark_time)[1],
        method = method, horizon = h, n = sum(usable),
        n_events = sum(y[usable] == 1),
        n_controls = sum(y[usable] == 0),
        n_censored_before_horizon = sum(!usable),
        event_rate   = weighted.mean(y[usable], w[usable]),
        auc_pooled   = weighted_auc(y, p, w),
        auc_mean     = if (length(fold_auc)) mean(fold_auc) else NA,
        auc_se       = if (length(fold_auc) > 1) sd(fold_auc)/sqrt(length(fold_auc)) else NA,
        brier_pooled = weighted.mean((y[usable] - p[usable])^2, w[usable]),
        brier_mean   = if (length(fold_brier)) mean(fold_brier) else NA,
        brier_se     = if (length(fold_brier) > 1) sd(fold_brier)/sqrt(length(fold_brier)) else NA,
        cal_intercept = cal[["calibration_intercept"]],
        cal_slope     = cal[["calibration_slope"]],
        cal_in_large  = cal[["calibration_in_large"]],
        stringsAsFactors = FALSE)
    }
  }
}
metrics <- do.call(rbind, metrics_list)
if (is.null(metrics) || !nrow(metrics)) {
  stop("No evaluable landmark-method-horizon cells were found.", call. = FALSE)
}
metrics$method_label <- ifelse(metrics$method %in% names(METHOD_LABELS),
                               METHOD_LABELS[metrics$method], metrics$method)
display_order <- unname(METHOD_LABELS[
  intersect(silk_opt("METHOD_ORDER"), names(METHOD_LABELS))
])
metrics$method_label <- factor(metrics$method_label, levels = display_order)

# Main figures show the primary benchmark arms; the same-feature ablation is
# retained in every results table and reported in the supplement.
figure_methods <- if (exists("MACS_PRIMARY_FIGURE_METHODS")) {
  intersect(MACS_PRIMARY_FIGURE_METHODS, unique(metrics$method))
} else unique(metrics$method)
metrics_figure <- metrics[metrics$method %in% figure_methods, , drop = FALSE]
metrics_figure$method_label <- droplevels(metrics_figure$method_label)

cat("\n=== Per-Horizon Metrics ===\n")
print(metrics[, c("landmark_set", "method_label","horizon","n","n_events","event_rate",
                   "auc_pooled","brier_pooled","cal_slope")],
      digits = 3, row.names = FALSE)

write.csv(metrics, file.path(RESULTS_DIR, "macs_metrics_per_horizon.csv"),
          row.names = FALSE)

# Complete comparator table for manuscript reporting. All prespecified
# landmark-horizon cells are retained; lower Brier scores are preferable.
comparator_table <- metrics[, c(
  "landmark_set", "landmark_time", "horizon", "method", "method_label",
  "n", "n_events", "n_controls", "n_censored_before_horizon",
  "auc_pooled", "brier_pooled", "cal_intercept", "cal_slope"
)]
method_rank <- match(comparator_table$method, silk_opt("METHOD_ORDER"))
comparator_table <- comparator_table[
  order(comparator_table$landmark_time, comparator_table$horizon, method_rank),
]
write.csv(comparator_table,
          file.path(RESULTS_DIR, "macs_comparator_table.csv"), row.names = FALSE)

# Primary benchmark contrasts: SILK against the recorded-clock time-only Cox
# (the dynamic prediction that inherits the origin error) and against
# parametric mixed-model landmarking.
silk_gain_table <- function(comparator_method, suffix) {
  primary_silk <- metrics[metrics$method == "SILK-Cox", c(
    "landmark_set", "landmark_time", "horizon", "n", "n_events",
    "auc_pooled", "brier_pooled", "cal_intercept", "cal_slope"
  )]
  comparator <- metrics[metrics$method == comparator_method, c(
    "landmark_set", "horizon", "auc_pooled", "brier_pooled",
    "cal_intercept", "cal_slope"
  )]
  if (!nrow(comparator)) return(NULL)
  names(primary_silk)[6:9] <- paste0(names(primary_silk)[6:9], "_silk")
  names(comparator)[3:6] <- paste0(names(comparator)[3:6], "_", suffix)
  comparison <- merge(primary_silk, comparator,
                      by = c("landmark_set", "horizon"), all = TRUE)
  comparison$auc_gain_silk <-
    comparison$auc_pooled_silk - comparison[[paste0("auc_pooled_", suffix)]]
  comparison$brier_gain_silk <-
    comparison[[paste0("brier_pooled_", suffix)]] - comparison$brier_pooled_silk
  comparison
}
primary_comparison <- silk_gain_table("Cox-Recorded", "recorded")
if (!is.null(primary_comparison)) {
  write.csv(primary_comparison,
            file.path(RESULTS_DIR, "macs_primary_silk_comparison.csv"),
            row.names = FALSE)
}
mmlm_comparison <- silk_gain_table("MMLM-Recorded", "mmlm")
if (!is.null(mmlm_comparison)) {
  write.csv(mmlm_comparison,
            file.path(RESULTS_DIR, "macs_silk_vs_mmlm.csv"),
            row.names = FALSE)
}

# ── 3. Profile-gap diagnostics ──────────────────────────────────────────────
diagnostics_file <- file.path(RESULTS_DIR, "macs_profile_diagnostics.csv")
if (file.exists(diagnostics_file)) {
  diagnostics <- read.csv(diagnostics_file, stringsAsFactors = FALSE)
  gap_cols <- intersect(c("gap_q1", "gap_q2", "gap"), names(diagnostics))
  if (nrow(diagnostics) && length(gap_cols)) {
    diag_summary <- aggregate(
      diagnostics[, gap_cols, drop = FALSE],
      diagnostics[, intersect(c("landmark_set", "method", "split"), names(diagnostics)),
                  drop = FALSE],
      function(x) c(median = median(x, na.rm = TRUE),
                    iqr = IQR(x, na.rm = TRUE))
    )
    diag_summary <- do.call(data.frame, diag_summary)
    names(diag_summary) <- sub("\\.", "_", names(diag_summary), fixed = TRUE)
    write.csv(diag_summary, file.path(RESULTS_DIR, "macs_profile_gap_summary.csv"),
              row.names = FALSE)
    cat("\nProfile-gap diagnostics saved to macs_profile_gap_summary.csv\n")
  }
}

# ── 4. Plotting setup ───────────────────────────────────────────────────────
method_colors <- c(
  "Recorded-clock Cox"                 = "grey50",
  "Recorded-clock Cox (same features)" = "grey70",
  "Parametric MMLM-Cox"                = UIUC_BLUE,
  "SILK-Cox (Gaussian)"                = UIUC_ORANGE
)
method_colors <- method_colors[intersect(names(method_colors), as.character(unique(metrics_figure$method_label)))]
landmark_facets <- length(unique(metrics_figure$landmark_set)) > 1L
landmark_labels <- c(
  "fixed_1y" = "1-year landmark", "fixed_2y" = "2-year landmark",
  "fixed_3y" = "3-year landmark", "fixed_5y" = "5-year landmark",
  "last_visit" = "Last observed visit"
)

# ── 5. Figure 1: AUC by horizon ─────────────────────────────────────────────
p1 <- ggplot(metrics_figure, aes(x = horizon, y = auc_pooled,
                           color = method_label, shape = method_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = HORIZONS) +
  labs(x = "Prediction horizon (years)",
       y = "Time-dependent IPCW AUC (higher is better)",
       color = "Method", shape = "Method",
       title = "Discrimination across prediction horizons") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
if (landmark_facets) {
  p1 <- p1 + facet_wrap(~ landmark_set,
                        labeller = as_labeller(landmark_labels))
}

# Figure file names match the manuscript's \includegraphics keys.
ggsave(file.path(FIGURES_DIR, "macs_auc_by_horizon.pdf"), p1, width=8, height=5.5, dpi=600)
ggsave(file.path(FIGURES_DIR, "macs_auc_by_horizon.png"), p1, width=8, height=5.5, dpi=300)

# ── 6. Figure 2: Brier score by horizon ─────────────────────────────────────
p2 <- ggplot(metrics_figure, aes(x = horizon, y = brier_pooled,
                           color = method_label, shape = method_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = HORIZONS) +
  labs(x = "Prediction horizon (years)",
       y = "IPCW Brier score (lower is better)",
       color = "Method", shape = "Method",
       title = "Prediction error across prediction horizons") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
if (landmark_facets) {
  p2 <- p2 + facet_wrap(~ landmark_set,
                        labeller = as_labeller(landmark_labels))
}

ggsave(file.path(FIGURES_DIR, "macs_brier_by_horizon.pdf"), p2, width=8, height=5.5, dpi=600)
ggsave(file.path(FIGURES_DIR, "macs_brier_by_horizon.png"), p2, width=8, height=5.5, dpi=300)

# ── 7. Figure 3: Calibration plots ──────────────────────────────────────────
cal_list <- list()
n_bins <- 5L
for (landmark_set in sort(unique(pred_eval$landmark_set))) {
  censor_survival <- make_censor_survival(
    outcome_data[outcome_data$landmark_set == landmark_set, , drop = FALSE]
  )
  # The calibration figure shows the same primary arms as the other figures.
  for (method in figure_methods) {
    for (h in HORIZONS) {
      z <- pred_eval[pred_eval$landmark_set == landmark_set &
                       pred_eval$method == method &
                       pred_eval$horizon == h, ]
      if (nrow(z) < 30) next
      hw <- ipcw_horizon(z$U, z$delta, h, censor_survival)
      z <- z[hw$usable, ]
      y <- hw$y[hw$usable]
      w <- hw$w[hw$usable]
      if (nrow(z) < 30) next
      q <- unique(quantile(z$risk_pred, seq(0, 1, length.out = n_bins+1), type=8))
      if (length(q) < 3) next
      z$bin <- cut(z$risk_pred, breaks = q, include.lowest = TRUE, labels = FALSE)
      for (b in sort(unique(z$bin))) {
        zb <- z[z$bin == b, ]
        wb <- w[z$bin == b]
        yb <- y[z$bin == b]
        cal_list[[length(cal_list)+1]] <- data.frame(
          landmark_set = landmark_set,
          method_label = ifelse(method %in% names(METHOD_LABELS),
                                METHOD_LABELS[method], method),
          horizon = h, mean_pred = weighted.mean(zb$risk_pred, wb),
          obs_rate = weighted.mean(yb, wb), stringsAsFactors = FALSE)
      }
    }
  }
}
cal_df <- if (length(cal_list)) {
  do.call(rbind, cal_list)
} else {
  data.frame(landmark_set = character(), method_label = character(), horizon = numeric(),
             mean_pred = numeric(), obs_rate = numeric())
}

calibration_horizon <- if (5 %in% cal_df$horizon) 5 else max(cal_df$horizon)
if (is.finite(calibration_horizon)) {
  cal_sub <- cal_df[cal_df$horizon == calibration_horizon, ]
  cal_sub$method_label <- droplevels(
    factor(cal_sub$method_label, levels = display_order)
  )
  cal_axis_max <- max(cal_sub$mean_pred, cal_sub$obs_rate, 0, na.rm = TRUE)
  cal_axis_max <- min(1, max(0.05, 1.05 * cal_axis_max))

  p3 <- ggplot(cal_sub, aes(x = mean_pred, y = obs_rate,
                              color = method_label, shape = method_label)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey40") +
    geom_point(size = 2.5) + geom_line(linewidth = 0.8) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(limits = c(0, cal_axis_max)) +
    scale_y_continuous(limits = c(0, cal_axis_max)) +
    coord_equal() +
    labs(x = "Predicted risk", y = "Observed risk (IPCW)",
         color = "Method", shape = "Method",
         title = paste0(calibration_horizon,
                        "-Year Risk Calibration by Landmark")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  if (length(unique(cal_sub$landmark_set)) > 1L) {
    p3 <- p3 + facet_wrap(~ landmark_set, ncol = 2,
                          labeller = as_labeller(landmark_labels))
  }

  ggsave(file.path(FIGURES_DIR, "macs_calibration.pdf"), p3,
         width=7.5, height=7, dpi=600)
  ggsave(file.path(FIGURES_DIR, "macs_calibration.png"), p3,
         width=7.5, height=7, dpi=300)
}

# ── 8. Figure 4: KM by SILK risk group ──────────────────────────────────────
target_h  <- max(HORIZONS)
silk_pred <- pred_eval[pred_eval$method == "SILK-Cox" &
                       pred_eval$horizon == target_h, ]
if (nrow(silk_pred) && length(unique(silk_pred$landmark_set)) > 1L) {
  preferred_landmark <- if ("fixed_5y" %in% silk_pred$landmark_set) "fixed_5y" else
    tail(sort(unique(silk_pred$landmark_set)), 1)
  silk_pred <- silk_pred[silk_pred$landmark_set == preferred_landmark, ]
}

if (nrow(silk_pred) > 30) {
  silk_pred$risk_group <- cut(
    silk_pred$risk_pred,
    breaks = quantile(silk_pred$risk_pred, c(0, 1/3, 2/3, 1)),
    labels = c("Low Risk", "Medium Risk", "High Risk"),
    include.lowest = TRUE
  )
  km_data <- silk_pred[, c("subject_id", "risk_group", "U", "delta")]
  km_fit  <- survival::survfit(survival::Surv(U, delta) ~ risk_group, data = km_data)

  for (ext in c("pdf", "png")) {
    if (ext == "pdf") pdf(file.path(FIGURES_DIR, "macs_km_risk_strata.pdf"), width=8, height=6)
    else png(file.path(FIGURES_DIR, "macs_km_risk_strata.png"), width=8, height=6, units="in", res=300)
    par(mar = c(5, 4.5, 3, 1))
    plot(km_fit, col = rep(UIUC_BLUE, 3),
         lty = c("dotted", "dashed", "solid"), lwd = c(1.5, 1.8, 2.2),
         xlab = "Time from landmark (years)",
         ylab = "AIDS-free survival probability",
         main = "AIDS-Free Survival by SILK-Cox Risk Group",
         xlim = c(0, min(10, max(km_data$U))))
    legend("bottomleft", legend = levels(km_data$risk_group),
           col = rep(UIUC_BLUE, 3), lty = c("dotted", "dashed", "solid"),
           lwd = c(1.5, 1.8, 2.2), bty = "n")
    dev.off()
  }

  lr <- survival::survdiff(survival::Surv(U, delta) ~ risk_group, data = km_data)
  cat("\nLog-rank test for SILK risk stratification:\n")
  print(lr)

  # Reported in the manuscript text: AIDS-free survival by SILK risk tertile.
  report_time <- min(5, max(km_data$U))
  strata_survival <- summary(km_fit, times = report_time, extend = TRUE)
  cat(sprintf("\nAIDS-free survival at %.1f years from the landmark by tertile:\n",
              report_time))
  print(data.frame(
    stratum = as.character(strata_survival$strata),
    survival = round(strata_survival$surv, 3)
  ), row.names = FALSE)
}

# ── 8b. Figure 5: registration diagnostics ──────────────────────────────────
# The manuscript reports the held-out shift distribution and the profile-gap
# diagnostic; both come from the package's registration output.
if (file.exists(diagnostics_file)) {
  diag_raw <- read.csv(diagnostics_file, stringsAsFactors = FALSE)
  held_out <- diag_raw[diag_raw$split == "test", , drop = FALSE]
  if (nrow(held_out) > 30 && "e_hat" %in% names(held_out)) {
    shift_range <- range(held_out$e_hat, na.rm = TRUE)
    pa <- ggplot(held_out, aes(x = e_hat)) +
      geom_histogram(bins = 40, fill = UIUC_ORANGE, colour = "white",
                     linewidth = 0.2) +
      geom_vline(xintercept = shift_range, linetype = "dashed",
                 colour = "grey40") +
      labs(x = "Estimated registration shift (years)", y = "Held-out subjects",
           title = "(a) Fitted shifts") +
      theme_minimal(base_size = 12)

    gap_data <- held_out[is.finite(held_out$gap_q1), , drop = FALSE]
    pb <- ggplot(gap_data, aes(x = gap_q1)) +
      geom_histogram(bins = 40, fill = UIUC_BLUE, colour = "white",
                     linewidth = 0.2) +
      geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey40") +
      labs(x = "Profile gap at one year", y = "Held-out subjects",
           title = "(b) Alignment gap") +
      theme_minimal(base_size = 12)

    combined <- NULL
    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined <- patchwork::wrap_plots(pa, pb, ncol = 2)
    }
    if (is.null(combined)) {
      # Fall back to a long-format facet so the figure does not depend on
      # an optional package.
      panel_data <- rbind(
        data.frame(panel = "(a) Fitted shifts", value = held_out$e_hat),
        data.frame(panel = "(b) Alignment gap", value = gap_data$gap_q1)
      )
      combined <- ggplot(panel_data, aes(x = value)) +
        geom_histogram(bins = 40, fill = UIUC_BLUE, colour = "white",
                       linewidth = 0.2) +
        facet_wrap(~ panel, scales = "free") +
        labs(x = NULL, y = "Held-out subjects") +
        theme_minimal(base_size = 12)
    }
    ggsave(file.path(FIGURES_DIR, "macs_registration_diagnostics.pdf"),
           combined, width = 9, height = 4, dpi = 600)
    ggsave(file.path(FIGURES_DIR, "macs_registration_diagnostics.png"),
           combined, width = 9, height = 4, dpi = 300)

    cat(sprintf(
      "\nHeld-out registration: n=%d, sd(shift)=%.2f, boundary=%.1f%%, gap<0.05=%.1f%%\n",
      nrow(held_out), sd(held_out$e_hat, na.rm = TRUE),
      100 * mean(as.logical(held_out$at_boundary), na.rm = TRUE),
      100 * mean(gap_data$gap_q1 < 0.05)
    ))
  }
}

# ── 9. Final report ─────────────────────────────────────────────────────────
cat("\n====================================================================\n")
cat("  MACS PDS Real-Data Analysis — Results Summary\n")
cat("====================================================================\n")
cat("\nData: 573 HIV seroconverters, 10,602 visits, 228 AIDS events (39.8%)\n")
cat("Design: fixed-landmark 5-fold cross-validation\n")
cat("Biomarkers: CD4 (LEU3N), CD4% (LEU3P), and CD8 (LEU2N)\n\n")

for (landmark_set in sort(unique(metrics$landmark_set))) {
  cat(sprintf("  Landmark = %s:\n", landmark_set))
  for (h in HORIZONS) {
    cat(sprintf("    Horizon = %g yr:\n", h))
    mh <- metrics[metrics$landmark_set == landmark_set & metrics$horizon == h, ]
    for (i in seq_len(nrow(mh)))
      cat(sprintf("      %-25s  AUC=%.3f  Brier=%.4f  Cal.slope=%.3f\n",
                  mh$method_label[i], mh$auc_pooled[i],
                  mh$brier_pooled[i], mh$cal_slope[i]))
  }
}
cat("\nResults:", RESULTS_DIR, "\nFigures:", FIGURES_DIR, "\n")
