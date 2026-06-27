# =============================================================================
# plotting.R
# Visualization functions (from plotting_prediction.R)
# =============================================================================

#' Method color palette
#'
#' Returns a named vector of colors for each method, suitable for ggplot2 scales.
#'
#' @return Named character vector of hex colors.
#' @export
method_palette <- function() {
  METHOD_ORDER <- silk_opt("METHOD_ORDER")
  UIUC_ORANGE <- silk_opt("UIUC_ORANGE")
  UIUC_BLUE <- silk_opt("UIUC_BLUE")
  method_colours <- c(
    "Landmark-Recorded" = "#3B7A57",
    "Cox-SameFeature-Recorded" = "#2E8B57",
    "MMLM-Recorded" = "#9467BD",
    "JM-Recorded" = "#8C564B",
    "Bayesian-Dynamic-Observed" = "#C44E52",
    "DeepSurv-Observed" = "#E377C2",
    "RSF-Observed" = "#17BECF",
    "TimeError-Integrated-Landmark" = "#2B8CBE",
    "SILK-MeanReg" = "#FDB863",
    "SILK-LinearMMD" = "#E08214",
    SILK = UIUC_ORANGE,
    "Beran-Recorded" = "#A6CEE3",
    "Beran-SILK" = "#FB9A99",
    "Beran-SILK-Linear" = "#E31A1C",
    "Beran-Oracle-Latent-Age" = "#6A3D9A",
    "Oracle-Latent-Age" = UIUC_BLUE
  )
  labels <- pretty_method(METHOD_ORDER)
  stats::setNames(unname(method_colours[METHOD_ORDER]), labels)
}

#' @keywords internal
scenario_axis_labels <- function(scenario) {
  SCENARIOS <- silk_opt("SCENARIOS")
  sc <- SCENARIOS[match(scenario, SCENARIOS$scenario), , drop = FALSE]
  error_type <- c(
    normal = "Normal",
    mixture = "Mixture",
    laplace = "Laplace",
    t3 = "T3",
    asymmetric = "Asymmetric"
  )
  et <- as.character(sc$eps_type)
  et <- ifelse(et %in% names(error_type), unname(error_type[et]), et)
  paste0(pretty_scenario(scenario), "\n", "Error SD = ", sc$sigma_eps, ", ", et)
}

#' @keywords internal
add_time_grid_columns <- function(df) {
  SCENARIO_ERROR_ORDER <- silk_opt("SCENARIO_ERROR_ORDER")
  SCENARIO_SHORT_LABELS <- silk_opt("SCENARIO_SHORT_LABELS")
  VISIT_SCHEDULES <- silk_opt("VISIT_SCHEDULES")

  if (!"time_grid_setting" %in% names(df)) return(df)
  parts <- strsplit(as.character(df$time_grid_setting), ":", fixed = TRUE)
  df$phase <- vapply(parts, function(x) if (length(x) >= 1L) x[1] else NA_character_, character(1))
  df$scenario <- vapply(parts, function(x) if (length(x) >= 2L) x[2] else NA_character_, character(1))
  df$schedule <- vapply(parts, function(x) if (length(x) >= 3L) x[3] else NA_character_, character(1))
  df$scenario_label <- factor(pretty_scenario(df$scenario), levels = pretty_scenario(SCENARIO_ERROR_ORDER))
  df$scenario_short <- factor(
    unname(SCENARIO_SHORT_LABELS[df$scenario]),
    levels = unname(SCENARIO_SHORT_LABELS[SCENARIO_ERROR_ORDER])
  )
  df$scenario_axis <- factor(scenario_axis_labels(df$scenario), levels = scenario_axis_labels(SCENARIO_ERROR_ORDER))
  df$schedule_label <- factor(pretty_schedule(df$schedule), levels = pretty_schedule(names(VISIT_SCHEDULES)))
  df
}

#' @keywords internal
add_method_columns <- function(df, method_col = "method") {
  METHOD_ORDER <- silk_opt("METHOD_ORDER")
  if (!method_col %in% names(df)) return(df)
  df <- df[df[[method_col]] %in% METHOD_ORDER, , drop = FALSE]
  label_col <- if (method_col == "method") "method_label" else paste0(method_col, "_label")
  df[[label_col]] <- factor(pretty_method(df[[method_col]]), levels = pretty_method(METHOD_ORDER))
  df
}

#' SILK ggplot2 theme
#'
#' A publication-quality ggplot2 theme for SILK visualizations.
#'
#' @param base_size Numeric. Base font size.
#' @return A ggplot2 theme object.
#' @export
theme_silk <- function(base_size = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for theme_silk", call. = FALSE)
  }
  if (is.null(base_size)) base_size <- silk_opt("FIGURE_BASE_SIZE")
  UIUC_BLUE <- silk_opt("UIUC_BLUE")
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = UIUC_BLUE, size = base_size + 4),
      plot.subtitle = ggplot2::element_text(size = base_size - 1),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 2),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 2, colour = UIUC_BLUE),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = "#B8B8B8"),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 2),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Save a ggplot to file
#'
#' @param plot A ggplot object.
#' @param out_dir Character. Output directory.
#' @param name Character. File name (without extension).
#' @param width Numeric. Figure width in inches.
#' @param height Numeric. Figure height in inches.
#' @return Invisible NULL.
#' @export
save_plot <- function(plot, out_dir, name, width = 12, height = 7) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for save_plot", call. = FALSE)
  }
  FIGURE_DPI <- silk_opt("FIGURE_DPI")
  png_file <- file.path(out_dir, paste0(name, ".png"))
  if (file.exists(png_file)) unlink(png_file)
  ggplot2::ggsave(
    png_file,
    plot,
    width = width,
    height = height,
    dpi = FIGURE_DPI,
    bg = "white"
  )
}

#' @keywords internal
rbind_fill <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, rows)
}

#' @keywords internal
wrap_labels <- function(x, width = 18) {
  vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

#' @keywords internal
auc_limits <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0.45, 1))
  lo <- max(0, floor((min(x) - 0.03) * 20) / 20)
  hi <- min(1, ceiling((max(x) + 0.03) * 20) / 20)
  if (hi - lo < 0.15) lo <- max(0, hi - 0.15)
  c(lo, hi)
}

#' @keywords internal
add_setting_label <- function(df, setting) {
  TIMEPOINT_SCHEDULES <- silk_opt("TIMEPOINT_SCHEDULES")
  if (setting == "n_train") {
    lv <- sort(unique(df$n_train_setting[is.finite(df$n_train_setting)]))
    df$setting_label <- factor(
      paste0("N = ", df$n_train_setting),
      levels = paste0("N = ", lv)
    )
  } else if (setting == "schedule") {
    df$setting_label <- factor(
      pretty_schedule(df$schedule),
      levels = pretty_schedule(TIMEPOINT_SCHEDULES)
    )
  } else {
    stop("Unknown setting: ", setting, call. = FALSE)
  }
  df
}

#' @keywords internal
design_subset <- function(df, phase, setting) {
  df <- df[df$phase == phase, , drop = FALSE]
  if (!nrow(df)) return(df)
  df <- add_setting_label(df, setting)
  df[!is.na(df$setting_label), , drop = FALSE]
}

#' @keywords internal
metric_label_for_pairs <- function(x) {
  out <- ifelse(x == "mean_auc", "Mean AUC",
                ifelse(x == "integrated_brier_score", "Integrated Brier Score", x))
  factor(out, levels = c("Integrated Brier Score", "Mean AUC"))
}

#' @keywords internal
add_scenario_plot_label <- function(df) {
  SCENARIO_PLOT_LABELS <- silk_opt("SCENARIO_PLOT_LABELS")
  SCENARIO_ERROR_ORDER <- silk_opt("SCENARIO_ERROR_ORDER")
  df$scenario_plot <- factor(
    unname(SCENARIO_PLOT_LABELS[df$scenario]),
    levels = unname(SCENARIO_PLOT_LABELS[SCENARIO_ERROR_ORDER])
  )
  df
}

#' @keywords internal
save_calibration_bias_plot <- function(ms, fig_pred_dir, src_dir, phase, setting,
                                       name, subtitle, width, height) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  UIUC_BLUE <- silk_opt("UIUC_BLUE")
  z <- design_subset(ms, phase, setting)
  z <- z[is.finite(z$calibration_in_large_mean), , drop = FALSE]
  if (!nrow(z)) return(invisible(FALSE))

  cal_source <- stats::aggregate(
    calibration_in_large_mean ~ scenario + scenario_short + method + method_label + setting_label,
    data = z,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(cal_source)[names(cal_source) == "calibration_in_large_mean"] <- "calibration_in_large"
  cal_source <- cal_source[is.finite(cal_source$calibration_in_large), , drop = FALSE]
  if (!nrow(cal_source)) return(invisible(FALSE))
  cal_source <- add_scenario_plot_label(cal_source)
  utils::write.csv(cal_source, file.path(src_dir, paste0(name, ".csv")), row.names = FALSE)

  p <- ggplot2::ggplot(
    cal_source,
    ggplot2::aes(x = .data$scenario_plot, y = .data$calibration_in_large,
                 colour = .data$method_label, group = .data$method_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = UIUC_BLUE, linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.72, alpha = 0.86) +
    ggplot2::geom_point(size = 2.15, alpha = 0.9) +
    ggplot2::facet_wrap(~ setting_label, ncol = 1) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_colour_manual(values = method_palette(), drop = FALSE) +
    ggplot2::labs(
      title = "Calibration Bias Across Error Scenarios",
      subtitle = subtitle,
      x = "Error Scenario",
      y = "Observed Minus Predicted Risk",
      colour = "Method"
    ) +
    theme_silk(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 42, hjust = 1, vjust = 1, size = 8.5, lineheight = 0.9),
      legend.position = "bottom"
    )
  save_plot(p, fig_pred_dir, name, width, height)
  invisible(TRUE)
}

#' @keywords internal
save_paired_difference_plot <- function(paired, fig_pred_dir, src_dir, phase, setting,
                                        name, subtitle, width, height) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  METHOD_ORDER <- silk_opt("METHOD_ORDER")
  UIUC_BLUE <- silk_opt("UIUC_BLUE")

  z <- paired[
    paired$metric %in% c("integrated_brier_score", "mean_auc") &
      paired$comparator != "Oracle-Latent-Age",
    ,
    drop = FALSE
  ]
  comparator_methods <- setdiff(METHOD_ORDER, c("SILK", "Oracle-Latent-Age"))
  z <- z[z$comparator %in% comparator_methods, , drop = FALSE]
  z <- design_subset(z, phase, setting)
  if (!nrow(z)) return(invisible(FALSE))
  z$metric_label <- metric_label_for_pairs(z$metric)
  comparator_levels <- pretty_method(comparator_methods)
  z$comparator_label <- factor(pretty_method(z$comparator), levels = comparator_levels)

  pair_source <- stats::aggregate(
    difference_mean ~ scenario + scenario_short + comparator + comparator_label + metric_label + setting_label,
    data = z,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  pair_source <- pair_source[is.finite(pair_source$difference_mean), , drop = FALSE]
  if (!nrow(pair_source)) return(invisible(FALSE))
  pair_source <- add_scenario_plot_label(pair_source)
  setting_levels <- levels(pair_source$setting_label)
  metric_levels <- levels(pair_source$metric_label)
  pair_source$panel_label <- factor(
    paste(pair_source$setting_label, pair_source$metric_label, sep = "\n"),
    levels = unlist(lapply(setting_levels, function(s) paste(s, metric_levels, sep = "\n")))
  )
  utils::write.csv(pair_source, file.path(src_dir, paste0(name, ".csv")), row.names = FALSE)

  p <- ggplot2::ggplot(
    pair_source,
    ggplot2::aes(x = .data$scenario_plot, y = .data$difference_mean,
                 group = .data$comparator_label, colour = .data$comparator_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = UIUC_BLUE, linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.72, alpha = 0.86) +
    ggplot2::geom_point(size = 2.1, alpha = 0.9) +
    ggplot2::facet_wrap(~ panel_label, ncol = 2, scales = "free_y") +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_colour_manual(values = method_palette()[comparator_levels], drop = TRUE) +
    ggplot2::labs(
      title = "SILK Minus Comparator Across Error Scenarios",
      subtitle = subtitle,
      x = "Error Scenario",
      y = "Mean Paired Difference",
      colour = "Comparator"
    ) +
    theme_silk(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 42, hjust = 1, vjust = 1, size = 8.5, lineheight = 0.9),
      legend.position = "bottom"
    )
  save_plot(p, fig_pred_dir, name, width, height)
  invisible(TRUE)
}

#' @keywords internal
save_oracle_gap_plot <- function(ms, fig_pred_dir, src_dir, phase, setting,
                                 name, subtitle, width, height) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  UIUC_BLUE <- silk_opt("UIUC_BLUE")

  oracle <- ms[
    ms$method == "Oracle-Latent-Age",
    c("time_grid_setting", "n_train_setting", "integrated_brier_score_mean"),
    drop = FALSE
  ]
  names(oracle)[3] <- "oracle_ibs"
  gap <- merge(ms, oracle, by = c("time_grid_setting", "n_train_setting"))
  gap <- gap[gap$method != "Oracle-Latent-Age", , drop = FALSE]
  gap <- design_subset(gap, phase, setting)
  gap$ibs_gap_to_oracle <- gap$integrated_brier_score_mean - gap$oracle_ibs
  gap <- gap[is.finite(gap$ibs_gap_to_oracle), , drop = FALSE]
  if (!nrow(gap)) return(invisible(FALSE))

  gap_source <- stats::aggregate(
    ibs_gap_to_oracle ~ scenario + scenario_short + method + method_label + setting_label,
    data = gap,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  if (!nrow(gap_source)) return(invisible(FALSE))
  gap_source <- add_scenario_plot_label(gap_source)
  utils::write.csv(gap_source, file.path(src_dir, paste0(name, ".csv")), row.names = FALSE)

  p <- ggplot2::ggplot(
    gap_source,
    ggplot2::aes(x = .data$scenario_plot, y = .data$ibs_gap_to_oracle,
                 colour = .data$method_label, group = .data$method_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = UIUC_BLUE, linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.72, alpha = 0.86) +
    ggplot2::geom_point(size = 2.1, alpha = 0.9) +
    ggplot2::facet_wrap(~ setting_label, ncol = 1) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_colour_manual(values = method_palette(), drop = TRUE) +
    ggplot2::labs(
      title = "Distance From The Latent-Age Oracle",
      subtitle = subtitle,
      x = "Error Scenario",
      y = "IBS Gap To Latent-Age Oracle",
      colour = "Method"
    ) +
    theme_silk(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 42, hjust = 1, vjust = 1, size = 8.5, lineheight = 0.9),
      legend.position = "bottom"
    )
  save_plot(p, fig_pred_dir, name, width, height)
  invisible(TRUE)
}

#' @keywords internal
make_dgm_visualization_data <- function(n = 350L) {
  SCENARIO_ERROR_ORDER <- silk_opt("SCENARIO_ERROR_ORDER")
  GLOBAL_SEED <- silk_opt("GLOBAL_SEED")

  subject_list <- vector("list", length(SCENARIO_ERROR_ORDER))
  visit_list <- vector("list", length(SCENARIO_ERROR_ORDER))
  for (k in seq_along(SCENARIO_ERROR_ORDER)) {
    scenario <- SCENARIO_ERROR_ORDER[k]
    dat <- generate_dataset_fixed(
      n = n,
      scenario_name = scenario,
      schedule_name = scenario_schedule(scenario),
      seed = GLOBAL_SEED + 90000L + k
    )
    s <- dat$subjects
    s$scenario <- scenario
    s$schedule <- dat$schedule
    s$scenario_label <- pretty_scenario(scenario)
    s$scenario_axis <- scenario_axis_labels(scenario)
    subject_list[[k]] <- s

    v <- dat$visits
    v$scenario <- scenario
    v$schedule <- dat$schedule
    v$scenario_label <- pretty_scenario(scenario)
    v$scenario_axis <- scenario_axis_labels(scenario)
    visit_list[[k]] <- v
  }
  subjects <- rbind_fill(subject_list)
  visits <- rbind_fill(visit_list)
  subjects$scenario_label <- factor(subjects$scenario_label, levels = pretty_scenario(SCENARIO_ERROR_ORDER))
  subjects$scenario_axis <- factor(subjects$scenario_axis, levels = scenario_axis_labels(SCENARIO_ERROR_ORDER))
  visits$scenario_label <- factor(visits$scenario_label, levels = pretty_scenario(SCENARIO_ERROR_ORDER))
  visits$scenario_axis <- factor(visits$scenario_axis, levels = scenario_axis_labels(SCENARIO_ERROR_ORDER))
  list(subjects = subjects, visits = visits)
}

#' Generate DGM visualization plots
#'
#' Creates publication-quality plots visualizing the data-generating mechanism
#' across all error scenarios.
#'
#' @param out_dir Character. Output directory for figures.
#' @param src_dir Character. Directory for source data CSVs.
#' @return Invisible TRUE on success.
#' @export
make_dgm_plots <- function(out_dir, src_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 is not installed; skipping DGM plots.")
    return(invisible(FALSE))
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    message("grid package is not available; skipping DGM plots.")
    return(invisible(FALSE))
  }

  SCENARIO_ERROR_ORDER <- silk_opt("SCENARIO_ERROR_ORDER")
  SCENARIOS <- silk_opt("SCENARIOS")
  UIUC_ORANGE <- silk_opt("UIUC_ORANGE")
  UIUC_BLUE <- silk_opt("UIUC_BLUE")

  dgm_dir <- file.path(out_dir, "figures", "dgm")
  dir.create(dgm_dir, recursive = TRUE, showWarnings = FALSE)

  dgm <- make_dgm_visualization_data()
  subjects <- dgm$subjects
  visits <- dgm$visits
  utils::write.csv(subjects, file.path(src_dir, "dgm_subjects_by_error_scenario.csv"), row.names = FALSE)
  utils::write.csv(visits, file.path(src_dir, "dgm_visits_by_error_scenario.csv"), row.names = FALSE)

  p_eps <- ggplot2::ggplot(subjects, ggplot2::aes(x = .data$scenario_axis, y = .data$eps)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = UIUC_BLUE, linewidth = 0.6) +
    ggplot2::geom_violin(fill = "#F9D2BA", colour = UIUC_ORANGE, linewidth = 0.5, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.16, outlier.alpha = 0.08, fill = "white", colour = UIUC_BLUE) +
    ggplot2::labs(
      title = "Origin-Error Scenarios In The Initial Data",
      subtitle = "Recorded Landmark Age Minus Latent Landmark Age, Ordered From No Error To Severe Error Regimes",
      x = "Scenario",
      y = "Origin Error"
    ) +
    theme_silk() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_plot(p_eps, dgm_dir, "DGM01_origin_error_distribution", 15, 8.5)

  key_scenarios <- intersect(
    c("no_error", "mean_moderate", "mean_severe", "mean_strong_dense",
      "dist_moderate", "dist_large", "dist_strong_dense", "heavy_tail", "asymmetric_shift"),
    SCENARIO_ERROR_ORDER
  )
  latent_obs <- subjects[subjects$scenario %in% key_scenarios, , drop = FALSE]
  p_latent <- ggplot2::ggplot(latent_obs, ggplot2::aes(x = .data$A_star, y = .data$A_obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = UIUC_BLUE, linewidth = 0.7) +
    ggplot2::geom_point(colour = UIUC_ORANGE, alpha = 0.22, size = 1.0) +
    ggplot2::facet_wrap(~ scenario_label, ncol = 3) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Recorded Versus Latent Landmark Age",
      subtitle = "The Diagonal Is The No-Error Target; Vertical Displacement Is The Origin-Time Shift",
      x = "True Latent Landmark Age",
      y = "Recorded Landmark Age"
    ) +
    theme_silk()
  save_plot(p_latent, dgm_dir, "DGM02_recorded_vs_latent_landmark_age", 13.5, 9.5)

  traj_scenarios <- c("no_error", "mean_moderate", "mean_severe", "mean_strong_dense", "dist_large", "dist_strong_dense")
  traj_scenarios <- traj_scenarios[traj_scenarios %in% SCENARIO_ERROR_ORDER]
  rep_subjects <- do.call(rbind, Map(function(scn, scenario_rank) {
    cand <- subjects[subjects$scenario == scn, , drop = FALSE]
    target <- stats::quantile(abs(cand$eps), probs = 0.90, na.rm = TRUE)
    if (!is.finite(target)) target <- 0
    cand <- cand[order(abs(abs(cand$eps) - target), decreasing = FALSE), , drop = FALSE]
    out <- cand[1, c("scenario", "id", "eps"), drop = FALSE]
    out$scenario_rank <- scenario_rank
    out
  }, traj_scenarios, seq_along(traj_scenarios)))
  traj_visits <- merge(
    visits[visits$scenario %in% traj_scenarios, , drop = FALSE],
    rep_subjects,
    by = c("scenario", "id"),
    all = FALSE,
    suffixes = c("", "_subject")
  )
  b1 <- first_biomarker_col(traj_visits)
  traj_visits <- traj_visits[order(traj_visits$scenario_rank, traj_visits$A_star_il), , drop = FALSE]
  traj_visits$facet_label <- paste0(wrap_labels(traj_visits$scenario_label, 20), "\nShift = ", round(traj_visits$eps, 1))
  facet_levels <- unique(traj_visits$facet_label)
  traj_visits$facet_label <- factor(traj_visits$facet_label, levels = facet_levels)

  latent_traj <- data.frame(
    facet_label = traj_visits$facet_label,
    age = traj_visits$A_star_il,
    value = traj_visits[[b1]],
    clock = "Latent Clock",
    stringsAsFactors = FALSE
  )
  recorded_traj <- data.frame(
    facet_label = traj_visits$facet_label,
    age = traj_visits$A_obs_il,
    value = traj_visits[[b1]],
    clock = "Recorded Clock",
    stringsAsFactors = FALSE
  )
  traj_long <- rbind(latent_traj, recorded_traj)
  traj_long$facet_label <- factor(traj_long$facet_label, levels = facet_levels)
  traj_long$clock <- factor(traj_long$clock, levels = c("Latent Clock", "Recorded Clock"))
  utils::write.csv(traj_visits, file.path(src_dir, "dgm_biomarker_trajectory_shift.csv"), row.names = FALSE)

  p_traj <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = traj_visits,
      ggplot2::aes(x = .data$A_star_il, xend = .data$A_obs_il, y = .data[[b1]], yend = .data[[b1]]),
      arrow = grid::arrow(length = grid::unit(0.11, "inches"), type = "closed"),
      colour = UIUC_ORANGE,
      alpha = 0.72,
      linewidth = 0.55
    ) +
    ggplot2::geom_line(
      data = traj_long,
      ggplot2::aes(x = .data$age, y = .data$value, colour = .data$clock, linetype = .data$clock),
      linewidth = 1.0
    ) +
    ggplot2::geom_point(
      data = traj_long,
      ggplot2::aes(x = .data$age, y = .data$value, colour = .data$clock),
      size = 2.2,
      alpha = 0.90
    ) +
    ggplot2::facet_wrap(~ facet_label, ncol = 3) +
    ggplot2::scale_colour_manual(values = c("Latent Clock" = UIUC_BLUE, "Recorded Clock" = UIUC_ORANGE)) +
    ggplot2::scale_linetype_manual(values = c("Latent Clock" = "solid", "Recorded Clock" = "dashed")) +
    ggplot2::labs(
      title = "How Origin Error Shifts A Biomarker Trajectory",
      subtitle = paste("Representative Subject Per Scenario. Arrows Move Each", b1, "Visit From Latent Age To Recorded Age."),
      x = "Age/Time Axis",
      y = paste(b1, "Value"),
      colour = "Clock",
      linetype = "Clock"
    ) +
    theme_silk()
  save_plot(p_traj, dgm_dir, "DGM03_biomarker_trajectory_shift", 14.5, 9.2)

  profile <- SCENARIOS[match(SCENARIO_ERROR_ORDER, SCENARIOS$scenario), , drop = FALSE]
  profile$scenario_axis <- factor(scenario_axis_labels(profile$scenario), levels = scenario_axis_labels(SCENARIO_ERROR_ORDER))
  utils::write.csv(profile, file.path(src_dir, "dgm_error_scenario_profile.csv"), row.names = FALSE)
  p_profile <- ggplot2::ggplot(profile, ggplot2::aes(x = .data$scenario_axis, y = .data$sigma_eps)) +
    ggplot2::geom_col(fill = UIUC_ORANGE, alpha = 0.9, width = 0.72) +
    ggplot2::geom_text(ggplot2::aes(label = .data$eps_type), angle = 90, hjust = -0.12, size = 4.8, colour = UIUC_BLUE) +
    ggplot2::labs(
      title = "Designed Origin-Error Severity Across Scenarios",
      subtitle = "Bars Show The Nominal Origin-Error Scale; Labels Show The Error Distribution Family",
      x = "Scenario",
      y = "Nominal Origin-Error Scale"
    ) +
    theme_silk() +
    ggplot2::coord_cartesian(ylim = c(0, max(profile$sigma_eps, na.rm = TRUE) + 2.2)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_plot(p_profile, dgm_dir, "DGM04_error_scenario_profile", 15, 8)

  invisible(TRUE)
}

#' Generate prediction performance plots
#'
#' Creates all prediction performance plots from summary CSV files.
#'
#' @param out_dir Character. Output directory containing summary/ subdirectory.
#' @return Invisible TRUE on success, FALSE if files not found.
#' @export
make_prediction_plots <- function(out_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 is not installed; skipping plots.")
    return(invisible(FALSE))
  }

  fig_pred_dir <- file.path(out_dir, "figures", "prediction")
  src_dir <- file.path(out_dir, "figures", "source_data")
  dir.create(fig_pred_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

  paired_file <- file.path(out_dir, "summary", "prediction_paired_differences.csv")
  method_file <- file.path(out_dir, "summary", "prediction_method_summary.csv")
  if (!file.exists(method_file) && !file.exists(paired_file)) {
    message("No prediction summary files found; skipping prediction plots.")
    make_dgm_plots(out_dir, src_dir)
    return(invisible(FALSE))
  }

  ms <- data.frame()
  if (file.exists(method_file)) {
    ms <- utils::read.csv(method_file, stringsAsFactors = FALSE)
    ms <- add_time_grid_columns(add_method_columns(ms))

    save_calibration_bias_plot(
      ms, fig_pred_dir, src_dir,
      phase = "primary", setting = "n_train",
      name = "PRED04A_calibration_bias_by_error_scenario_ntrain",
      subtitle = "Primary Design: Scenarios Are Shown Separately For Each Training-Set Size.",
      width = 16, height = 12
    )
    save_calibration_bias_plot(
      ms, fig_pred_dir, src_dir,
      phase = "timepoint_extension", setting = "schedule",
      name = "PRED04B_calibration_bias_by_error_scenario_timepoints",
      subtitle = "Timepoint Design: N = 400, With Scenarios Shown Separately For Each Visit Schedule.",
      width = 16, height = 17
    )
    save_oracle_gap_plot(
      ms, fig_pred_dir, src_dir,
      phase = "primary", setting = "n_train",
      name = "PRED07A_gap_to_latent_oracle_by_error_scenario_ntrain",
      subtitle = "Primary Design: Integrated Brier Score Gap By Training-Set Size.",
      width = 16, height = 12
    )
    save_oracle_gap_plot(
      ms, fig_pred_dir, src_dir,
      phase = "timepoint_extension", setting = "schedule",
      name = "PRED07B_gap_to_latent_oracle_by_error_scenario_timepoints",
      subtitle = "Timepoint Design: N = 400, With Integrated Brier Score Gap By Visit Schedule.",
      width = 16, height = 17
    )
  }

  paired <- data.frame()
  if (file.exists(paired_file)) {
    paired <- utils::read.csv(paired_file, stringsAsFactors = FALSE)
    paired <- add_time_grid_columns(paired)

    save_paired_difference_plot(
      paired, fig_pred_dir, src_dir,
      phase = "primary", setting = "n_train",
      name = "PRED06A_silk_minus_comparator_by_error_scenario_ntrain",
      subtitle = "Primary Design: Brier Differences Below Zero Favor SILK; AUC Differences Above Zero Favor SILK.",
      width = 18, height = 12
    )
    save_paired_difference_plot(
      paired, fig_pred_dir, src_dir,
      phase = "timepoint_extension", setting = "schedule",
      name = "PRED06B_silk_minus_comparator_by_error_scenario_timepoints",
      subtitle = "Timepoint Design With N = 400: Brier Differences Below Zero Favor SILK; AUC Differences Above Zero Favor SILK.",
      width = 18, height = 16
    )
  }

  make_dgm_plots(out_dir, src_dir)
  return(invisible(TRUE))
}
