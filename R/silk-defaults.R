# =============================================================================
# silk-defaults.R
# Package options/defaults converted from cfg.R
# =============================================================================

#' Default options for the SILK package
#'
#' Returns a named list of all default configuration values used by SILK,
#' including simulation parameters, DGM parameters, visit schedules,
#' registration tuning constants, labels, and palettes.
#'
#' @return A named list of default option values.
#' @export
silk_default_options <- function() {
  list(
    # General design
    GLOBAL_SEED = 20260530L,
    LATENT_AGE_MIN = 24,
    LATENT_AGE_MAX = 32,
    A_STAR_CENTER = 28,
    A_STAR_SD = (32 - 24) / sqrt(12),

    # Visit schedules
    VISIT_SCHEDULES = list(
      m2  = c(6, 0),
      m3  = c(6, 3, 0),
      m4  = c(6, 4, 2, 0),
      m8  = seq(6, 0, length.out = 8),
      m12 = seq(6, 0, length.out = 12),
      m20 = seq(6, 0, length.out = 20)
    ),
    PRIMARY_SCHEDULE = "m4",

    # Scenario lists
    SCENARIOS_ALL = c(
      "no_error", "mean_moderate", "mean_severe", "mean_strong_dense",
      "dist_moderate", "dist_large", "dist_strong_dense",
      "weak_stage", "biased_shift", "heavy_tail",
      "asymmetric_shift", "irregular_missing"
    ),
    TIMEPOINT_SCHEDULES = c("m2", "m4", "m8", "m12", "m20"),

    # Simulation sample sizes
    N_TRAIN_GRID = c(200L, 400L, 1000L),
    N_TEST = 1000L,
    N_REP = 50L,
    TIMEPOINT_N_TRAIN = 400L,

    # DGM parameters
    SURV_LAMBDA_D = 10.0,
    SURV_KAPPA_D = 1.35,
    SURV_BETA_A = 0.95,
    SURV_BETA_X1 = 0.20,
    SURV_BETA_X2 = 0.15,
    SURV_BETA_U = 0.25,
    CENS_RATE = 0.045,
    EVAL_HORIZON = 4.0,
    PREDICTION_HORIZONS = c(1, 2, 3, 4),

    # Biomarker controls
    N_BIO_MEAN = 4L,
    N_BIO_DIST = 20L,
    N_BIO_WEAK = 6L,
    SIGMA_BIO_MEAN = 0.70,
    SIGMA_BIO_WEAK = 1.25,
    DIST_SD_BASE = 0.45,
    DIST_SD_SLOPE = 0.18,
    DEFAULT_SIGNAL_AMP = 1.35,
    DEFAULT_U_BIO_COEF = 0.15,

    # Registration
    SHIFT_GRID_STEP = 0.2,
    DEFAULT_SHIFT_GRID_MIN = -8,
    DEFAULT_SHIFT_GRID_MAX = 8,
    H_A_BANDWIDTH = 1.00,
    H_X_BANDWIDTH = 1.50,
    REGISTRATION_KERNEL = "rbf",
    REGISTRATION_KERNEL_APPROX = "exact",
    KERNEL_RFF_DIM = 512L,
    KERNEL_RFF_SEED = 20260530L,
    KERNEL_MATERN_NU = 1.5,
    KERNEL_POLYNOMIAL_DEGREE = 2L,
    KERNEL_POLYNOMIAL_COEF0 = 1.0,
    KERNEL_POLYNOMIAL_SCALE = 0.5,
    KERNEL_POLYNOMIAL_NONNEGATIVE = TRUE,
    KERNEL_CHUNK_SIZE = 2500L,
    N_ALT_ITER_MAX = 7L,
    ALT_TOL = 1e-3,
    N_STARTS = 4L,
    RANDOM_START_SD = 0.75,
    N_FOLDS = 5L,
    ANCHOR_MODE = "rx",
    PROFILE_TEMPERATURE = 0.015,
    PROFILE_LOCAL_RADIUS = 0.8,

    # Methods
    METHODS = c(
      "Landmark-Recorded", "Cox-SameFeature-Recorded",
      "MMLM-Recorded", "JM-Recorded",
      "RSF-Observed", "DeepSurv-Observed", "Bayesian-Dynamic-Observed",
      "TimeError-Integrated-Landmark", "SILK", "SILK-LinearMMD",
      "Beran-Recorded", "Beran-SILK", "Beran-SILK-Linear",
      "Beran-Oracle-Latent-Age", "Oracle-Latent-Age"
    ),
    METHOD_ORDER = c(
      "Landmark-Recorded", "Cox-SameFeature-Recorded",
      "MMLM-Recorded", "JM-Recorded",
      "RSF-Observed", "DeepSurv-Observed", "Bayesian-Dynamic-Observed",
      "TimeError-Integrated-Landmark", "SILK", "SILK-LinearMMD",
      "Beran-Recorded", "Beran-SILK", "Beran-SILK-Linear",
      "Beran-Oracle-Latent-Age", "Oracle-Latent-Age"
    ),

    SURVIVAL_HISTORY_BIOMARKERS = 3L,
    SURVIVAL_HISTORY_KEEP = c("X1", "X2"),

    # Scenario definitions
    SCENARIOS = data.frame(
      scenario = c(
        "no_error", "mean_moderate", "mean_severe", "mean_strong_dense",
        "dist_moderate", "dist_large", "dist_strong_dense",
        "weak_stage", "biased_shift", "heavy_tail",
        "asymmetric_shift", "irregular_missing"
      ),
      biomarker_signal = c(
        "mean", "mean", "mean", "mean_strong", "distribution", "distribution",
        "distribution_strong", "weak", "mean", "mean", "mean", "mean"
      ),
      default_schedule = c(
        "m4", "m4", "m4", "m12", "m4", "m4", "m12", "m4", "m4", "m4", "m4", "m4"
      ),
      sigma_eps = c(0, 3.0, 10.0, 10.0, 3.0, 10.0, 10.0, 10.0, 6.0, 10.0, 10.0, 10.0),
      eps_mean = c(0, 0, 0, 0, 0, 0, 0, 0, 2.5, 0, 0, 0),
      eps_type = c(
        "normal", "normal", "mixture", "mixture", "normal", "mixture",
        "mixture", "mixture", "normal", "t3", "asymmetric", "mixture"
      ),
      n_biomarkers = c(4L, 4L, 4L, 5L, 20L, 20L, 40L, 6L, 4L, 4L, 4L, 4L),
      missing_rate = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25),
      irregular = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
      shift_min = c(-4, -12, -30, -30, -12, -30, -30, -30, -24, -35, -35, -30),
      shift_max = c( 4,  12,  30,  30,  12,  30,  30,  30,  30,  35,  35,  30),
      signal_amp = c(1.35, 1.65, 2.35, 4.25, 0.12, 0.12, 0.10, 0.08, 2.00, 2.00, 2.00, 2.00),
      sigma_bio = c(0.70, 0.60, 0.50, 0.22, 0.70, 0.60, 0.30, 1.45, 0.55, 0.55, 0.55, 0.65),
      u_bio_coef = c(0.15, 0.12, 0.08, 0.02, 0.12, 0.08, 0.02, 0.20, 0.10, 0.10, 0.10, 0.10),
      dist_sd_base = c(0.45, 0.45, 0.45, 0.45, 0.36, 0.22, 0.14, 0.45, 0.45, 0.45, 0.45, 0.45),
      dist_sd_slope = c(0.18, 0.18, 0.18, 0.18, 0.35, 0.85, 1.25, 0.18, 0.18, 0.18, 0.18, 0.18),
      risk_beta_A = c(0.95, 1.05, 1.65, 2.10, 1.05, 1.65, 2.10, 1.10, 1.20, 1.20, 1.20, 1.20),
      risk_beta_U = c(0.15, 0.12, 0.08, 0.03, 0.12, 0.08, 0.03, 0.18, 0.10, 0.10, 0.10, 0.10),
      description = c(
        "No origin error; calibration should not help and may add noise.",
        "Mean-stage biomarkers; moderate common origin error.",
        "Mean-stage biomarkers; severe common origin error.",
        "Favorable sanity regime: severe error, strong mean signal, dense visits, low biomarker noise.",
        "Distributional biomarkers; moderate common origin error.",
        "Distributional biomarkers; severe common origin error.",
        "Favorable distributional sanity regime: severe error, strong dispersion signal, dense visits.",
        "Weak-stage negative-control regime.",
        "Biased-shift sensitivity regime; mean-zero anchor is misspecified.",
        "Heavy-tailed common-shift regime.",
        "Asymmetric common-shift regime.",
        "Irregular visits with missing prelandmark biomarkers."
      ),
      stringsAsFactors = FALSE
    ),

    # Labels
    SCENARIO_LABELS = c(
      no_error = "No Error",
      mean_moderate = "Mean Biomarker, Moderate Error",
      mean_severe = "Mean Biomarker, Severe Error",
      mean_strong_dense = "Mean Biomarker, Severe Error, Dense Visits",
      dist_moderate = "Distribution Biomarker, Moderate Error",
      dist_large = "Distribution Biomarker, Severe Error",
      dist_strong_dense = "Distribution Biomarker, Severe Error, Dense Visits",
      weak_stage = "Weak Biomarker Signal",
      biased_shift = "Biased Origin Shift",
      heavy_tail = "Heavy-Tailed Origin Error",
      asymmetric_shift = "Asymmetric Origin Error",
      irregular_missing = "Irregular Visits With Missing Data"
    ),

    # Canonical architecture display names (registration x survival layer).
    METHOD_LABELS = c(
      "Landmark-Recorded" = "Landmark-Recorded",
      "Cox-SameFeature-Recorded" = "Recorded-Cox",
      "MMLM-Recorded" = "MMLM-Recorded",
      "JM-Recorded" = "JM-Recorded",
      "RSF-Observed" = "RSF-Observed",
      "DeepSurv-Observed" = "DeepSurv-Observed",
      "Bayesian-Dynamic-Observed" = "Bayesian-Dynamic-Observed",
      "TimeError-Integrated-Landmark" = "MI Back-Calc Landmark",
      SILK = "SILK-Gaussian-Cox",
      "SILK-LinearMMD" = "SILK-Linear-Cox",
      "Beran-Recorded" = "Recorded-Beran",
      "Beran-SILK" = "SILK-Gaussian-Beran",
      "Beran-SILK-Linear" = "SILK-Linear-Beran",
      "Beran-Oracle-Latent-Age" = "Oracle-Beran",
      "Oracle-Latent-Age" = "Oracle-Cox",
      "SILK-MeanReg" = "SILK-MeanTraj (suppl.)"
    ),

    SCHEDULE_LABELS = c(
      m2 = "2 Visits",
      m3 = "3 Visits",
      m4 = "4 Visits",
      m8 = "8 Visits",
      m12 = "12 Visits",
      m20 = "20 Visits"
    ),

    # Plotting constants
    UIUC_ORANGE = "#FF5F05",
    UIUC_BLUE = "#13294B",
    FIGURE_DPI = 600,
    FIGURE_BASE_SIZE = 18,

    SCENARIO_ERROR_ORDER = c(
      "no_error", "mean_moderate", "dist_moderate", "biased_shift",
      "mean_severe", "dist_large", "weak_stage", "heavy_tail",
      "asymmetric_shift", "irregular_missing",
      "mean_strong_dense", "dist_strong_dense"
    ),

    SCENARIO_SHORT_LABELS = c(
      no_error = "No\nError",
      mean_moderate = "Mean Signal\nModerate",
      dist_moderate = "Distribution\nModerate",
      biased_shift = "Biased\nShift",
      mean_severe = "Mean Signal\nSevere",
      dist_large = "Distribution\nSevere",
      weak_stage = "Weak\nSignal",
      heavy_tail = "Heavy-Tailed\nShift",
      asymmetric_shift = "Asymmetric\nShift",
      irregular_missing = "Irregular\nMissing",
      mean_strong_dense = "Mean Signal\nDense Severe",
      dist_strong_dense = "Distribution\nDense Severe"
    ),

    SCENARIO_PLOT_LABELS = c(
      no_error = "No\nError",
      mean_moderate = "Mean Signal\nModerate",
      dist_moderate = "Distribution\nModerate",
      biased_shift = "Biased\nShift",
      mean_severe = "Mean Signal\nSevere",
      dist_large = "Distribution\nSevere",
      weak_stage = "Weak\nSignal",
      heavy_tail = "Heavy-Tailed\nShift",
      asymmetric_shift = "Asymmetric\nShift",
      irregular_missing = "Irregular\nMissing",
      mean_strong_dense = "Mean Signal\nDense Severe",
      dist_strong_dense = "Distribution\nDense Severe"
    )
  )
}

#' Get a SILK package option
#'
#' Retrieves a single SILK package option by name.
#'
#' @param name Character string. The option name (e.g., "SHIFT_GRID_STEP").
#' @return The current value of the option.
#' @export
silk_opt <- function(name) {
  getOption(paste0("silk.", name))
}

#' Set SILK package options
#'
#' Sets one or more SILK package options.
#'
#' @param ... Named arguments where names are option names and values are the
#'   new option values.
#' @return Invisible NULL.
#' @export
#' @examples
#' silk_options(SHIFT_GRID_STEP = 0.1)
#' silk_opt("SHIFT_GRID_STEP")
silk_options <- function(...) {
  args <- list(...)
  for (nm in names(args)) {
    options(stats::setNames(list(args[[nm]]), paste0("silk.", nm)))
  }
  invisible(NULL)
}

#' Pretty-print scenario names
#'
#' @param x Character vector of scenario identifiers.
#' @return Character vector of human-readable scenario labels.
#' @export
pretty_scenario <- function(x) {
  labels <- silk_opt("SCENARIO_LABELS")
  unname(ifelse(x %in% names(labels), labels[x], gsub("_", " ", x)))
}

#' Pretty-print method names
#'
#' @param x Character vector of method identifiers.
#' @return Character vector of human-readable method labels.
#' @export
pretty_method <- function(x) {
  labels <- silk_opt("METHOD_LABELS")
  unname(ifelse(x %in% names(labels), labels[x], gsub("_", " ", x)))
}

#' Pretty-print schedule names
#'
#' @param x Character vector of schedule identifiers.
#' @return Character vector of human-readable schedule labels.
#' @export
pretty_schedule <- function(x) {
  labels <- silk_opt("SCHEDULE_LABELS")
  unname(ifelse(x %in% names(labels), labels[x], x))
}

#' Get scenario definition
#'
#' Retrieve the full row from the SCENARIOS data frame for a given scenario name.
#'
#' @param scenario_name Character string. One of the defined scenario names.
#' @return A one-row data frame with scenario parameters.
#' @export
get_scenario <- function(scenario_name) {
  scenarios <- silk_opt("SCENARIOS")
  z <- scenarios[scenarios$scenario == scenario_name, , drop = FALSE]
  if (nrow(z) != 1L) stop("Unknown scenario: ", scenario_name, call. = FALSE)
  z
}

#' Build shift grid for a scenario
#'
#' Creates the sequence of candidate shift values used by registration.
#'
#' @param scenario_name Character string. A scenario name.
#' @return Numeric vector of shift grid values.
#' @export
make_shift_grid <- function(scenario_name) {
  sc <- get_scenario(scenario_name)
  seq(sc$shift_min, sc$shift_max, by = silk_opt("SHIFT_GRID_STEP"))
}

#' Get default schedule for a scenario
#'
#' @param scenario_name Character string. A scenario name.
#' @param fallback Character string. Fallback schedule name.
#' @return Character string schedule name.
#' @keywords internal
scenario_schedule <- function(scenario_name, fallback = NULL) {
  if (is.null(fallback)) fallback <- silk_opt("PRIMARY_SCHEDULE")
  sc <- get_scenario(scenario_name)
  sch <- sc$default_schedule[1]
  if (is.na(sch) || !nzchar(sch)) sch <- fallback
  sch
}
