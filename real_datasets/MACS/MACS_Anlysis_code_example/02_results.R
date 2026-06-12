# =============================================================================
# 02_results.R
# Evaluate predictions and produce tables + figures
# =============================================================================

source("00_setup.R")
library(ggplot2)

# ── 1. Load predictions ─────────────────────────────────────────────────────
predictions <- read.csv(file.path(RESULTS_DIR, "macs_predictions.csv"),
                        stringsAsFactors = FALSE)
predictions$at_risk <- as.logical(predictions$at_risk)

cat("Loaded", nrow(predictions), "prediction rows\n")
cat("Methods:", paste(unique(predictions$method), collapse = ", "), "\n")
cat("Horizons:", paste(sort(unique(predictions$horizon)), collapse = ", "), "\n")

METHOD_LABELS <- silk_opt("METHOD_LABELS")
UIUC_ORANGE   <- silk_opt("UIUC_ORANGE")
UIUC_BLUE     <- silk_opt("UIUC_BLUE")

# ── 2. Per-horizon metrics (only at-risk subjects) ──────────────────────────
pred_atrisk <- predictions[predictions$at_risk, , drop = FALSE]

metrics_list <- list()
for (method in unique(pred_atrisk$method)) {
  for (h in sort(unique(pred_atrisk$horizon))) {
    z <- pred_atrisk[pred_atrisk$method == method &
                     pred_atrisk$horizon == h, , drop = FALSE]
    if (nrow(z) < 10) next

    y <- z$event_within_horizon
    p <- z$risk_pred

    # Per-fold AUC and Brier
    fold_auc   <- numeric(0)
    fold_brier <- numeric(0)
    for (ff in unique(z$fold_id)) {
      zf <- z[z$fold_id == ff, , drop = FALSE]
      if (length(unique(zf$event_within_horizon)) == 2) {
        fold_auc   <- c(fold_auc, binary_auc(zf$event_within_horizon, zf$risk_pred))
        fold_brier <- c(fold_brier, mean((zf$event_within_horizon - zf$risk_pred)^2))
      }
    }

    cal <- calibration_stats(y, p)

    metrics_list[[length(metrics_list) + 1]] <- data.frame(
      method    = method,
      horizon   = h,
      n         = nrow(z),
      event_rate = mean(y),
      auc_mean  = if (length(fold_auc)) mean(fold_auc) else NA,
      auc_se    = if (length(fold_auc) > 1) sd(fold_auc) / sqrt(length(fold_auc)) else NA,
      brier_mean = if (length(fold_brier)) mean(fold_brier) else NA,
      brier_se   = if (length(fold_brier) > 1) sd(fold_brier) / sqrt(length(fold_brier)) else NA,
      brier_pooled     = mean((y - p)^2),
      auc_pooled       = binary_auc(y, p),
      cal_intercept    = cal[["calibration_intercept"]],
      cal_slope        = cal[["calibration_slope"]],
      cal_in_large     = cal[["calibration_in_large"]],
      stringsAsFactors = FALSE
    )
  }
}
metrics <- do.call(rbind, metrics_list)
rownames(metrics) <- NULL

# Pretty labels
metrics$method_label <- ifelse(
  metrics$method %in% names(METHOD_LABELS),
  METHOD_LABELS[metrics$method],
  metrics$method
)

cat("\n=== Per-Horizon Metrics ===\n")
print(metrics[, c("method_label", "horizon", "n", "event_rate",
                   "auc_pooled", "brier_pooled", "cal_slope")],
      digits = 3, row.names = FALSE)

write.csv(metrics,
          file.path(RESULTS_DIR, "macs_metrics_per_horizon.csv"),
          row.names = FALSE)

# ── 3. Summary table ────────────────────────────────────────────────────────
summary_list <- list()
for (method in unique(metrics$method)) {
  z <- metrics[metrics$method == method, , drop = FALSE]
  summary_list[[method]] <- data.frame(
    method       = method,
    method_label = z$method_label[1],
    mean_auc     = mean(z$auc_pooled, na.rm = TRUE),
    mean_brier   = mean(z$brier_pooled, na.rm = TRUE),
    mean_cal_slope     = mean(z$cal_slope, na.rm = TRUE),
    mean_cal_intercept = mean(z$cal_intercept, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

cat("\n=== Summary Across Horizons ===\n")
print(summary_df, digits = 3, row.names = FALSE)

write.csv(summary_df,
          file.path(RESULTS_DIR, "macs_summary.csv"),
          row.names = FALSE)

# ── 4. Figure 1: AUC by horizon ─────────────────────────────────────────────
method_colors <- c(
  "Landmark Cox"          = "grey50",
  "Mixed-Model Landmark"  = UIUC_BLUE,
  "SILK"                  = UIUC_ORANGE
)

p1 <- ggplot(metrics, aes(x = horizon, y = auc_pooled,
                           color = method_label, shape = method_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = auc_mean - 1.96 * auc_se,
                     ymax = auc_mean + 1.96 * auc_se),
                width = 0.15, linewidth = 0.6) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = sort(unique(metrics$horizon))) +
  labs(
    x     = "Prediction Horizon (years from landmark)",
    y     = "AUC (pooled across CV folds)",
    color = "Method", shape = "Method",
    title = "Discrimination: AUC by Prediction Horizon",
    subtitle = "MACS PDS — 573 HIV seroconverters, 5-fold CV"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(file.path(FIGURES_DIR, "fig1_auc_by_horizon.pdf"),
       p1, width = 8, height = 5.5, dpi = 600)
ggsave(file.path(FIGURES_DIR, "fig1_auc_by_horizon.png"),
       p1, width = 8, height = 5.5, dpi = 300)
cat("Saved fig1_auc_by_horizon\n")

# ── 5. Figure 2: Brier score by horizon ─────────────────────────────────────
p2 <- ggplot(metrics, aes(x = horizon, y = brier_pooled,
                           color = method_label, shape = method_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = brier_mean - 1.96 * brier_se,
                     ymax = brier_mean + 1.96 * brier_se),
                width = 0.15, linewidth = 0.6) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = sort(unique(metrics$horizon))) +
  labs(
    x     = "Prediction Horizon (years from landmark)",
    y     = "Brier Score (pooled across CV folds)",
    color = "Method", shape = "Method",
    title = "Prediction Accuracy: Brier Score by Horizon",
    subtitle = "MACS PDS — lower is better"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(file.path(FIGURES_DIR, "fig2_brier_by_horizon.pdf"),
       p2, width = 8, height = 5.5, dpi = 600)
ggsave(file.path(FIGURES_DIR, "fig2_brier_by_horizon.png"),
       p2, width = 8, height = 5.5, dpi = 300)
cat("Saved fig2_brier_by_horizon\n")

# ── 6. Figure 3: Calibration plots ──────────────────────────────────────────
cal_data <- list()
n_bins <- 10L
for (method in unique(pred_atrisk$method)) {
  for (h in sort(unique(pred_atrisk$horizon))) {
    z <- pred_atrisk[pred_atrisk$method == method &
                     pred_atrisk$horizon == h, , drop = FALSE]
    if (nrow(z) < 20) next

    q <- quantile(z$risk_pred, probs = seq(0, 1, length.out = n_bins + 1),
                  na.rm = TRUE, type = 8)
    q <- unique(q)
    if (length(q) < 3) next
    z$bin <- cut(z$risk_pred, breaks = q, include.lowest = TRUE, labels = FALSE)

    for (b in sort(unique(z$bin))) {
      zb <- z[z$bin == b, , drop = FALSE]
      cal_data[[length(cal_data) + 1]] <- data.frame(
        method       = method,
        method_label = ifelse(method %in% names(METHOD_LABELS),
                              METHOD_LABELS[method], method),
        horizon      = h,
        bin          = b,
        n            = nrow(zb),
        mean_pred    = mean(zb$risk_pred),
        obs_rate     = mean(zb$event_within_horizon),
        stringsAsFactors = FALSE
      )
    }
  }
}
cal_df <- do.call(rbind, cal_data)
rownames(cal_df) <- NULL

# Focus on key horizons
key_horizons <- c(2, 5)
key_horizons <- key_horizons[key_horizons %in% cal_df$horizon]
cal_sub <- cal_df[cal_df$horizon %in% key_horizons, , drop = FALSE]
cal_sub$horizon_label <- paste0("Horizon = ", cal_sub$horizon, " yr")

if (nrow(cal_sub) > 0) {
  p3 <- ggplot(cal_sub, aes(x = mean_pred, y = obs_rate,
                              color = method_label, shape = method_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(size = 2.5) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ horizon_label, scales = "free") +
    scale_color_manual(values = method_colors) +
    coord_equal() +
    labs(
      x     = "Predicted Risk",
      y     = "Observed Proportion",
      color = "Method", shape = "Method",
      title = "Calibration: Predicted vs Observed Risk",
      subtitle = "MACS PDS — 10-bin calibration"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")

  ggsave(file.path(FIGURES_DIR, "fig3_calibration.pdf"),
         p3, width = 10, height = 5, dpi = 600)
  ggsave(file.path(FIGURES_DIR, "fig3_calibration.png"),
         p3, width = 10, height = 5, dpi = 300)
  cat("Saved fig3_calibration\n")
}

# ── 7. Figure 4: SILK registration — estimated shift distribution ───────────
# Re-run SILK on full data to get registration estimates
cat("\nRunning SILK registration on full data for shift visualization...\n")
full_grid <- make_shift_grid(SCENARIO_NAME)

full_reg <- tryCatch({
  fit_registration_multistart(
    visits, subjects,
    feature_type = "silk",
    grid = full_grid,
    seed = 42L
  )
}, error = function(e) {
  cat("Registration failed:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(full_reg)) {
  shift_df <- data.frame(
    id     = full_reg$ids,
    e_hat  = full_reg$e_train
  )
  shift_df <- merge(shift_df, subjects[, c("id", "A_obs")], by = "id")
  shift_df$S_hat <- shift_df$A_obs - shift_df$e_hat  # calibrated landmark age

  write.csv(shift_df,
            file.path(RESULTS_DIR, "macs_silk_shifts.csv"),
            row.names = FALSE)

  p4a <- ggplot(shift_df, aes(x = e_hat)) +
    geom_histogram(fill = UIUC_ORANGE, color = "white",
                   bins = 40, alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    labs(
      x     = expression(hat(epsilon)[i] ~ "(estimated origin shift, years)"),
      y     = "Count",
      title = "SILK-Estimated Origin Shifts",
      subtitle = "MACS PDS — 573 seroconverters"
    ) +
    theme_minimal(base_size = 14)

  p4b <- ggplot(shift_df, aes(x = A_obs, y = S_hat)) +
    geom_point(alpha = 0.4, size = 1.5, color = UIUC_BLUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    labs(
      x     = expression(hat(A)[obs] ~ "(observed landmark age, years)"),
      y     = expression(hat(A)^"*" ~ "(SILK-calibrated landmark age, years)"),
      title = "Observed vs SILK-Calibrated Landmark Ages"
    ) +
    theme_minimal(base_size = 14)

  ggsave(file.path(FIGURES_DIR, "fig4a_shift_distribution.pdf"),
         p4a, width = 7, height = 4.5, dpi = 600)
  ggsave(file.path(FIGURES_DIR, "fig4a_shift_distribution.png"),
         p4a, width = 7, height = 4.5, dpi = 300)
  ggsave(file.path(FIGURES_DIR, "fig4b_calibrated_landmark.pdf"),
         p4b, width = 6, height = 5.5, dpi = 600)
  ggsave(file.path(FIGURES_DIR, "fig4b_calibrated_landmark.png"),
         p4b, width = 6, height = 5.5, dpi = 300)
  cat("Saved fig4a and fig4b\n")
}

# ── 8. Figure 5: Risk stratification ────────────────────────────────────────
# Show Kaplan-Meier curves stratified by SILK-predicted risk groups
target_horizon <- max(HORIZONS)
silk_pred <- pred_atrisk[pred_atrisk$method == "SILK" &
                         pred_atrisk$horizon == target_horizon, , drop = FALSE]

if (nrow(silk_pred) > 20) {
  silk_pred$risk_group <- cut(
    silk_pred$risk_pred,
    breaks = quantile(silk_pred$risk_pred, probs = c(0, 1/3, 2/3, 1)),
    labels = c("Low Risk", "Medium Risk", "High Risk"),
    include.lowest = TRUE
  )

  # Merge back to subjects for KM
  km_data <- merge(silk_pred[, c("subject_id", "risk_group")],
                   subjects, by.x = "subject_id", by.y = "id")
  km_data <- km_data[!duplicated(km_data$subject_id), , drop = FALSE]

  km_fit <- survival::survfit(survival::Surv(U, delta) ~ risk_group,
                              data = km_data)

  pdf(file.path(FIGURES_DIR, "fig5_km_risk_strata.pdf"), width = 8, height = 6)
  par(mar = c(5, 4.5, 3, 1))
  plot(km_fit,
       col  = c(UIUC_BLUE, "grey50", UIUC_ORANGE),
       lwd  = 2,
       xlab = "Time from Landmark (years)",
       ylab = "AIDS-Free Survival Probability",
       main = "Kaplan-Meier by SILK Risk Group",
       xlim = c(0, min(10, max(km_data$U))))
  legend("bottomleft",
         legend = levels(km_data$risk_group),
         col    = c(UIUC_BLUE, "grey50", UIUC_ORANGE),
         lwd    = 2, bty = "n")
  dev.off()

  png(file.path(FIGURES_DIR, "fig5_km_risk_strata.png"),
      width = 8, height = 6, units = "in", res = 300)
  par(mar = c(5, 4.5, 3, 1))
  plot(km_fit,
       col  = c(UIUC_BLUE, "grey50", UIUC_ORANGE),
       lwd  = 2,
       xlab = "Time from Landmark (years)",
       ylab = "AIDS-Free Survival Probability",
       main = "Kaplan-Meier by SILK Risk Group",
       xlim = c(0, min(10, max(km_data$U))))
  legend("bottomleft",
         legend = levels(km_data$risk_group),
         col    = c(UIUC_BLUE, "grey50", UIUC_ORANGE),
         lwd    = 2, bty = "n")
  dev.off()

  cat("Saved fig5_km_risk_strata\n")

  # Log-rank test
  lr <- survival::survdiff(survival::Surv(U, delta) ~ risk_group, data = km_data)
  cat("\nLog-rank test for risk stratification:\n")
  print(lr)
}

# ── 9. Print final report ───────────────────────────────────────────────────
cat("\n")
cat("====================================================================\n")
cat("  MACS PDS Real-Data Analysis — Results Summary\n")
cat("====================================================================\n")
cat("\nData: 573 HIV seroconverters, 10,602 visits, 228 AIDS events (39.8%)\n")
cat("Design: 5-fold cross-validation\n")
cat("Biomarkers: CD4, CD4%, CD8, log10(viral load)\n\n")

cat("Per-horizon metrics:\n")
for (h in sort(unique(metrics$horizon))) {
  cat(sprintf("\n  Horizon = %g years:\n", h))
  mh <- metrics[metrics$horizon == h, , drop = FALSE]
  for (i in seq_len(nrow(mh))) {
    cat(sprintf("    %-25s  AUC = %.3f  Brier = %.4f  Cal.slope = %.3f\n",
                mh$method_label[i], mh$auc_pooled[i],
                mh$brier_pooled[i], mh$cal_slope[i]))
  }
}

cat("\n\nAll results saved to:", RESULTS_DIR, "\n")
cat("All figures saved to:", FIGURES_DIR, "\n")
