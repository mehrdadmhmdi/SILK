# =============================================================================
# cfg.R
# Multifaceted SILK simulation configuration, v4.
#
# Simulation target:
#   The primary target is absolute risk over fixed prediction horizons:
#   P(T <= a + horizon | T > a, information available by landmark a).
#   Recorded-age methods use recorded time; SILK uses its calibrated landmark
#   score; the latent-age oracle is simulation-only.
# =============================================================================

# ----------------------------- Paths -----------------------------------------
SIM_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
R_DIR <- file.path(SIM_ROOT, "R")
OUT_DIR_ENV <- Sys.getenv("SILK_OUT_DIR", unset = "")
OUT_DIR <- if (nzchar(OUT_DIR_ENV)) {
  normalizePath(OUT_DIR_ENV, winslash = "/", mustWork = FALSE)
} else {
  file.path(SIM_ROOT, "outputs")
}
RAW_DIR <- file.path(OUT_DIR, "raw")
SUMMARY_DIR <- file.path(OUT_DIR, "summary")
FIG_DIR <- file.path(OUT_DIR, "figures")
LOG_DIR <- file.path(OUT_DIR, "logs")
for (d in c(RAW_DIR, SUMMARY_DIR, FIG_DIR, LOG_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ----------------------------- General design --------------------------------
GLOBAL_SEED <- 20260530L
LATENT_AGE_MIN <- 24
LATENT_AGE_MAX <- 32
A_STAR_CENTER <- 28
A_STAR_SD <- (LATENT_AGE_MAX - LATENT_AGE_MIN) / sqrt(12)

# Visit times are stored as time before landmark. All schedules end at lag 0.
VISIT_SCHEDULES <- list(
  m2  = c(6, 0),
  m3  = c(6, 3, 0),
  m4  = c(6, 4, 2, 0),
  m8  = seq(6, 0, length.out = 8),
  m12 = seq(6, 0, length.out = 12),
  m20 = seq(6, 0, length.out = 20)
)

PRIMARY_SCHEDULE <- "m4"
TIMEPOINT_SCENARIOS <- SCENARIOS_ALL <- c(
  "no_error",
  "mean_moderate",
  "mean_severe",
  "mean_strong_dense",
  "dist_moderate",
  "dist_large",
  "dist_strong_dense",
  "weak_stage",
  "biased_shift",
  "heavy_tail",
  "asymmetric_shift",
  "irregular_missing"
)
TIMEPOINT_SCHEDULES <- c("m2", "m4", "m8", "m12", "m20")
TIMEPOINT_N_TRAIN <- as.integer(Sys.getenv("SILK_TP_N_TRAIN", unset = "400"))

N_TRAIN_GRID <- as.integer(strsplit(Sys.getenv("SILK_N_TRAIN_GRID", unset = "200,400,1000"), ",")[[1]])
N_TEST <- as.integer(Sys.getenv("SILK_N_TEST", unset = "1000"))
N_REP <- as.integer(Sys.getenv("SILK_N_REP", unset = "50"))
CODECHECK_N_TRAIN <- as.integer(Sys.getenv("SILK_CODECHECK_N_TRAIN", unset = "300"))
CODECHECK_N_TEST <- as.integer(Sys.getenv("SILK_CODECHECK_N_TEST", unset = "1000"))
CODECHECK_N_REP <- as.integer(Sys.getenv("SILK_CODECHECK_N_REP", unset = "2"))
CODECHECK_SCENARIOS <- c("mean_strong_dense", "dist_strong_dense", "mean_severe")

# ----------------------------- DGM parameters --------------------------------
SURV_LAMBDA_D <- 10.0
SURV_KAPPA_D <- 1.35
# Gompertz baseline on attained age. Growth 0.540 (was 0.18) puts enough of the
# prognostic signal on the corrupted time scale that the Bayes-optimal mean AUC
# is 0.808 and the recorded-versus-latent gap is 0.23 AUC, instead of 0.640 and
# 0.055. The base rate is retuned so the four-year event rate stays at 0.24.
SURV_AGE_BASE_RATE <- as.numeric(Sys.getenv("SILK_SURV_AGE_RATE", unset = "0.01107"))
SURV_AGE_GROWTH <- as.numeric(Sys.getenv("SILK_SURV_AGE_GROWTH", unset = "0.540"))
SURV_BETA_A <- 0.95
SURV_BETA_X1 <- 0.20
SURV_BETA_X2 <- 0.15
SURV_BETA_U <- 0.25
CENS_RATE <- 0.045
EVAL_HORIZON <- 4.0
PREDICTION_HORIZONS <- as.numeric(strsplit(Sys.getenv("SILK_HORIZONS", unset = "1,2,3,4"), ",")[[1]])
PREDICTION_HORIZONS <- sort(unique(PREDICTION_HORIZONS[is.finite(PREDICTION_HORIZONS) & PREDICTION_HORIZONS > 0]))
if (!length(PREDICTION_HORIZONS)) PREDICTION_HORIZONS <- EVAL_HORIZON

# Default biomarker controls; scenario-specific columns override these.
N_BIO_MEAN <- 4L
N_BIO_DIST <- 20L
N_BIO_WEAK <- 6L
SIGMA_BIO_MEAN <- 0.70
SIGMA_BIO_WEAK <- 1.25
DIST_SD_BASE <- 0.45
DIST_SD_SLOPE <- 0.18
DEFAULT_SIGNAL_AMP <- 1.35
DEFAULT_U_BIO_COEF <- 0.15

# ----------------------------- Registration ----------------------------------
SHIFT_GRID_STEP <- as.numeric(Sys.getenv("SILK_SHIFT_GRID_STEP", unset = "0.2"))
DEFAULT_SHIFT_GRID_MIN <- -8
DEFAULT_SHIFT_GRID_MAX <- 8
H_A_BANDWIDTH <- as.numeric(Sys.getenv("SILK_HA", unset = "1.00"))
H_X_BANDWIDTH <- as.numeric(Sys.getenv("SILK_HX", unset = "1.50"))
H_X2_BANDWIDTH <- as.numeric(Sys.getenv("SILK_HX2", unset = "1.00"))
H_LAG_BANDWIDTH <- as.numeric(Sys.getenv("SILK_HLAG", unset = "0.50"))
BIOMARKER_BANDWIDTH <- Sys.getenv("SILK_BIOMARKER_BW", unset = "median")
# Multiplier on the median heuristic. 1 is confirmatory; sweep a prespecified
# ladder (0.25, 0.5, 1, 2, 4) into a separate output directory to test whether
# the median rule is on the right scale rather than assuming it is not.
BIOMARKER_BANDWIDTH_SCALE <- as.numeric(
  Sys.getenv("SILK_BIOMARKER_BW_SCALE", unset = "1")
)
BIOMARKER_KERNELS <- unique(trimws(strsplit(
  Sys.getenv("SILK_BIOMARKER_KERNELS", unset = "gaussian,laplace,matern32"), ","
)[[1]]))
BIOMARKER_KERNELS <- BIOMARKER_KERNELS[nzchar(BIOMARKER_KERNELS)]
allowed_biomarker_kernels <- c("gaussian", "laplace", "matern32")
if (!length(BIOMARKER_KERNELS) || length(setdiff(BIOMARKER_KERNELS, allowed_biomarker_kernels))) {
  stop("SILK_BIOMARKER_KERNELS must be a comma-separated subset of gaussian,laplace,matern32.")
}
BIOMARKER_BANDWIDTH_MAX_POINTS <- as.integer(
  Sys.getenv("SILK_BIOMARKER_BW_MAX_POINTS", unset = "500")
)
TEMPLATE_INPUT_FEATURES <- as.integer(
  Sys.getenv("SILK_TEMPLATE_INPUT_FEATURES", unset = "16")
)
TEMPLATE_FEATURE_SEED <- as.integer(
  Sys.getenv("SILK_TEMPLATE_FEATURE_SEED", unset = "271828")
)
TEMPLATE_RIDGE_LAMBDA <- as.numeric(
  Sys.getenv("SILK_TEMPLATE_RIDGE", unset = "0.001")
)
TEMPLATE_NUMERICAL_JITTER <- as.numeric(
  Sys.getenv("SILK_TEMPLATE_JITTER", unset = "1e-10")
)
TEMPLATE_QUERY_CHUNK <- as.integer(
  Sys.getenv("SILK_TEMPLATE_QUERY_CHUNK", unset = "5000")
)
N_ALT_ITER_MAX <- as.integer(Sys.getenv("SILK_ALT_ITER", unset = "7"))
ALT_TOL <- 1e-3
# Localized profile search. Measured on a pilot (mean_severe, n = 200, m4):
# removing the localization degraded held-out shift RMSE from 2.13 to 8.62 and
# the correlation with the true shift from 0.92 to 0.11, so the restriction is
# part of the estimator and must be disclosed in the manuscript rather than
# removed. Set SILK_PROFILE_SEARCH_RADIUS=Inf for the sensitivity arm.
ALT_LOCAL_RADIUS <- as.numeric(Sys.getenv("SILK_ALT_LOCAL_RADIUS", unset = "2.0"))
N_STARTS <- as.integer(Sys.getenv("SILK_N_STARTS", unset = "3"))
RANDOM_START_SD <- 0.75
N_FOLDS <- as.integer(Sys.getenv("SILK_N_FOLDS", unset = "5"))
ANCHOR_MODE <- Sys.getenv("SILK_ANCHOR_MODE", unset = "rx") # "intercept" or "rx"
PROFILE_TEMPERATURE <- as.numeric(Sys.getenv("SILK_PROFILE_TEMP", unset = "0.015"))
PROFILE_LOCAL_RADIUS <- as.numeric(Sys.getenv("SILK_PROFILE_LOCAL_RADIUS", unset = "0.8"))
PROFILE_SEARCH_RADIUS <- as.numeric(Sys.getenv("SILK_PROFILE_SEARCH_RADIUS", unset = "2.0"))
CLOCK_SIGNAL_MIN <- as.numeric(Sys.getenv("SILK_CLOCK_SIGNAL_MIN", unset = "0.10"))

# ----------------------------- Methods ---------------------------------------
# The roster has exactly twelve methods.
#
# Biomarker-free clock benchmarks:
#   Recorded-Cox, Recorded-Beran, Oracle-Cox, Oracle-Beran.
# These use the relevant age coordinate and baseline covariates (X1, X2) only.
#
# Biomarker-using methods:
#   the three SILK registration kernels, MMLM, JM, RSF, DeepSurv, and
#   TimeError-Integrated-Landmark.
# No history-augmented Cox variants or SILK-Beran variant are included.
KERNEL_METHODS <- c(
  gaussian = "SILK-Gaussian",
  laplace = "SILK-Laplace",
  matern32 = "SILK-Matern32"
)
ACTIVE_KERNEL_METHODS <- unname(KERNEL_METHODS[BIOMARKER_KERNELS])

METHODS <- c(
  ACTIVE_KERNEL_METHODS,
  "Recorded-Cox",
  "Recorded-Beran",
  "MMLM",
  "JM",
  "RSF",
  "DeepSurv",
  "TimeError-Integrated-Landmark",
  "Oracle-Cox",
  "Oracle-Beran"
)

# Prespecified analysis roles. The reporting script must not infer these.
PRIMARY_SILK_METHOD <- "SILK-Gaussian"
PRIMARY_TRIPLE <- c(
  recorded = "Recorded-Cox",
  silk = PRIMARY_SILK_METHOD,
  oracle = "Oracle-Cox"
)

# Complete feasible procedures a practitioner could run. Prespecified, so the
# leaderboard cannot silently absorb sibling SILK kernels or oracle bounds.
PRIMARY_COMPETITORS <- c(
  "Recorded-Cox",
  "Recorded-Beran",
  "MMLM",
  "JM",
  "RSF",
  "DeepSurv",
  "TimeError-Integrated-Landmark"
)
SILK_FAMILY_METHODS <- ACTIVE_KERNEL_METHODS
ORACLE_METHODS <- c("Oracle-Cox", "Oracle-Beran")
BIOMARKER_FREE_METHODS <- c(
  "Recorded-Cox", "Recorded-Beran", "Oracle-Cox", "Oracle-Beran"
)
METHOD_ORDER <- METHODS

TIMEERROR_N_IMPUTE <- as.integer(Sys.getenv("SILK_TIMEERROR_N_IMPUTE", unset = "5"))
TIMEERROR_SD <- as.numeric(Sys.getenv("SILK_TIMEERROR_SD", unset = "1.5"))
TIMEERROR_USE_TRUTH <- identical(tolower(Sys.getenv("SILK_TIMEERROR_USE_TRUTH", unset = "false")), "true")
ENABLE_JMBAYES2 <- identical(tolower(Sys.getenv("SILK_ENABLE_JMBAYES2", unset = "true")), "true")
JMBAYES_N_CHAINS <- as.integer(Sys.getenv("SILK_JMBAYES_N_CHAINS", unset = "1"))
JMBAYES_N_ITER <- as.integer(Sys.getenv("SILK_JMBAYES_N_ITER", unset = "1000"))
JMBAYES_N_BURNIN <- as.integer(Sys.getenv("SILK_JMBAYES_N_BURNIN", unset = "500"))
JMBAYES_PRED_N_SAMPLES <- as.integer(Sys.getenv("SILK_JMBAYES_PRED_N_SAMPLES", unset = "200"))
JMBAYES_LANDMARK_TIME <- as.numeric(Sys.getenv("SILK_JMBAYES_LANDMARK_TIME", unset = "10"))
if (!is.finite(JMBAYES_LANDMARK_TIME) ||
    JMBAYES_LANDMARK_TIME <= max(unlist(VISIT_SCHEDULES, use.names = FALSE))) {
  stop(
    "SILK_JMBAYES_LANDMARK_TIME must be finite and exceed the longest pre-landmark visit lag.",
    call. = FALSE
  )
}
OBSERVED_ML_USE_CURRENT_BIOMARKER <- identical(tolower(Sys.getenv("SILK_OBSERVED_ML_USE_CURRENT_BIOMARKER", unset = "true")), "true")
RSF_NUM_TREES <- as.integer(Sys.getenv("SILK_RSF_NUM_TREES", unset = "300"))
RSF_MIN_NODE_SIZE <- as.integer(Sys.getenv("SILK_RSF_MIN_NODE_SIZE", unset = "15"))
DEEPSURV_EPOCHS <- as.integer(Sys.getenv("SILK_DEEPSURV_EPOCHS", unset = "80"))
DEEPSURV_BATCH_SIZE <- as.integer(Sys.getenv("SILK_DEEPSURV_BATCH_SIZE", unset = "128"))
DEEPSURV_VALIDATION_FRACTION <- as.numeric(
  Sys.getenv("SILK_DEEPSURV_VALIDATION_FRACTION", unset = "0.20")
)
DEEPSURV_PATIENCE <- as.integer(Sys.getenv("SILK_DEEPSURV_PATIENCE", unset = "10"))

# ----------------------------- Scenarios -------------------------------------
# Error scales are fixed through age reliability, not selected by prediction
# performance. With Var(A*) = 64/12, sigma 2.5 and 5 correspond to approximate
# reliabilities 0.46 and 0.18, respectively.
SCENARIOS <- data.frame(
  scenario = SCENARIOS_ALL,
  biomarker_signal = c(
    "mean", "mean", "mean", "mean_strong", "distribution", "distribution",
    "distribution_strong", "stage_null", "mean", "mean", "mean", "mean"
  ),
  default_schedule = c(
    "m4", "m4", "m4", "m12", "m4", "m4", "m12", "m4", "m4", "m4", "m4", "m4"
  ),
  # Offsets scaled so the recorded clock is genuinely uninformative: age
  # reliability is 0.13 at sigma_eps = 6 and 0.036 at 12, giving recorded-age
  # AUC 0.63 and 0.58 against an oracle of 0.81.
  sigma_eps = c(0, 6.0, 12.0, 12.0, 6.0, 12.0, 12.0, 12.0, 10.0, 12.0, 12.0, 12.0),
  eps_mean = c(0, 0, 0, 0, 0, 0, 0, 0, 6.0, 0, 0, 0),
  eps_type = c(
    "normal", "normal", "mixture", "mixture", "normal", "mixture",
    "mixture", "mixture", "normal", "t3", "asymmetric", "mixture"
  ),
  n_biomarkers = c(4L, 4L, 4L, 5L, 20L, 20L, 40L, 6L, 4L, 4L, 4L, 4L),
  missing_rate = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25),
  irregular = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  shift_min = c(-4, -20, -36, -36, -20, -36, -36, -36, -30, -48, -48, -36),
  shift_max = c( 4,  20,  36,  36,  20,  36,  36,  36,  42,  48,  48,  36),
  signal_amp = c(1.35, 1.65, 2.35, 4.25, 0.12, 0.12, 0.10, 0.00, 2.00, 2.00, 2.00, 2.00),
  sigma_bio = c(0.70, 0.60, 0.50, 0.22, 0.70, 0.60, 0.30, 1.45, 0.55, 0.55, 0.55, 0.65),
  u_bio_coef = rep(0, 12),
  dist_sd_base = c(0.45, 0.45, 0.45, 0.45, 0.36, 0.22, 0.14, 0.45, 0.45, 0.45, 0.45, 0.45),
  dist_sd_slope = c(0.18, 0.18, 0.18, 0.18, 0.35, 0.85, 1.25, 0.18, 0.18, 0.18, 0.18, 0.18),
  risk_age_growth = rep(SURV_AGE_GROWTH, 12),
  risk_beta_U = rep(0, 12),
  description = c(
    "No origin error; calibration should not help and may add noise.",
    "Mean-stage biomarkers; moderate common origin error.",
    "Mean-stage biomarkers; severe common origin error.",
    "Favorable sanity regime: severe error, strong mean signal, dense visits, low biomarker noise.",
    "Distributional biomarkers; moderate common origin error.",
    "Distributional biomarkers; severe common origin error.",
    "Favorable distributional sanity regime: severe error, strong dispersion signal, dense visits.",
    "Prespecified stage-null negative control: biomarkers contain no latent-age or frailty signal.",
    "Biased-shift sensitivity regime; mean-zero anchor is misspecified.",
    "Heavy-tailed common-shift regime.",
    "Asymmetric common-shift regime.",
    "Irregular visits with missing prelandmark biomarkers."
  ),
  stringsAsFactors = FALSE
)

SCENARIO_LABELS <- c(
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
)

# Canonical architecture display names (registration x survival layer).
METHOD_LABELS <- c(
  "SILK-Gaussian" = "SILK-Gaussian",
  "SILK-Laplace" = "SILK-Laplace",
  "SILK-Matern32" = "SILK-Matérn-3/2",
  "Recorded-Cox" = "Recorded-Cox",
  "Recorded-Beran" = "Recorded-Beran",
  "MMLM" = "MMLM",
  "JM" = "JM",
  "RSF" = "RSF",
  "DeepSurv" = "DeepSurv",
  "TimeError-Integrated-Landmark" = "TimeError-Integrated-Landmark",
  "Oracle-Cox" = "Oracle-Cox",
  "Oracle-Beran" = "Oracle-Beran"
)

SCHEDULE_LABELS <- c(
  m2 = "2 Visits",
  m3 = "3 Visits",
  m4 = "4 Visits",
  m8 = "8 Visits",
  m12 = "12 Visits",
  m20 = "20 Visits"
)

pretty_scenario <- function(x) unname(ifelse(x %in% names(SCENARIO_LABELS), SCENARIO_LABELS[x], gsub("_", " ", x)))
pretty_method <- function(x) unname(ifelse(x %in% names(METHOD_LABELS), METHOD_LABELS[x], gsub("_", " ", x)))
pretty_schedule <- function(x) unname(ifelse(x %in% names(SCHEDULE_LABELS), SCHEDULE_LABELS[x], x))

get_scenario <- function(scenario_name) {
  z <- SCENARIOS[SCENARIOS$scenario == scenario_name, , drop = FALSE]
  if (nrow(z) != 1L) stop("Unknown scenario: ", scenario_name, call. = FALSE)
  z
}

make_shift_grid <- function(scenario_name) {
  sc <- get_scenario(scenario_name)
  seq(sc$shift_min, sc$shift_max, by = SHIFT_GRID_STEP)
}

scenario_schedule <- function(scenario_name, fallback = PRIMARY_SCHEDULE) {
  sc <- get_scenario(scenario_name)
  sch <- sc$default_schedule[1]
  if (is.na(sch) || !nzchar(sch)) sch <- fallback
  sch
}

# The installed package is the single source of truth for data generation,
# registration, SILK survival layers, and evaluation. Mirror the simulation's
# environment-controlled design into package options once, before tasks run.
registration_cores <- Sys.getenv("SILK_N_CORES", unset = "auto")
if (grepl("^[0-9]+$", registration_cores)) {
  registration_cores <- as.integer(registration_cores)
}
silk_options(
  GLOBAL_SEED = GLOBAL_SEED,
  LATENT_AGE_MIN = LATENT_AGE_MIN,
  LATENT_AGE_MAX = LATENT_AGE_MAX,
  A_STAR_CENTER = A_STAR_CENTER,
  A_STAR_SD = A_STAR_SD,
  VISIT_SCHEDULES = VISIT_SCHEDULES,
  PRIMARY_SCHEDULE = PRIMARY_SCHEDULE,
  SCENARIOS_ALL = SCENARIOS_ALL,
  TIMEPOINT_SCHEDULES = TIMEPOINT_SCHEDULES,
  SURV_LAMBDA_D = SURV_LAMBDA_D,
  SURV_KAPPA_D = SURV_KAPPA_D,
  SURV_AGE_BASE_RATE = SURV_AGE_BASE_RATE,
  SURV_AGE_GROWTH = SURV_AGE_GROWTH,
  SURV_BETA_A = SURV_BETA_A,
  SURV_BETA_X1 = SURV_BETA_X1,
  SURV_BETA_X2 = SURV_BETA_X2,
  SURV_BETA_U = SURV_BETA_U,
  CENS_RATE = CENS_RATE,
  EVAL_HORIZON = EVAL_HORIZON,
  PREDICTION_HORIZONS = PREDICTION_HORIZONS,
  N_BIO_MEAN = N_BIO_MEAN,
  N_BIO_DIST = N_BIO_DIST,
  N_BIO_WEAK = N_BIO_WEAK,
  SIGMA_BIO_MEAN = SIGMA_BIO_MEAN,
  SIGMA_BIO_WEAK = SIGMA_BIO_WEAK,
  DIST_SD_BASE = DIST_SD_BASE,
  DIST_SD_SLOPE = DIST_SD_SLOPE,
  DEFAULT_SIGNAL_AMP = DEFAULT_SIGNAL_AMP,
  DEFAULT_U_BIO_COEF = DEFAULT_U_BIO_COEF,
  SHIFT_GRID_STEP = SHIFT_GRID_STEP,
  DEFAULT_SHIFT_GRID_MIN = DEFAULT_SHIFT_GRID_MIN,
  DEFAULT_SHIFT_GRID_MAX = DEFAULT_SHIFT_GRID_MAX,
  BIOMARKER_KERNEL = "gaussian",
  BIOMARKER_BANDWIDTH = BIOMARKER_BANDWIDTH,
  BIOMARKER_BANDWIDTH_SCALE = BIOMARKER_BANDWIDTH_SCALE,
  BIOMARKER_BANDWIDTH_MAX_POINTS = BIOMARKER_BANDWIDTH_MAX_POINTS,
  TEMPLATE_INPUT_FEATURES = TEMPLATE_INPUT_FEATURES,
  TEMPLATE_INPUT_COVARIATES = c("X1", "X2", "lag"),
  TEMPLATE_FEATURE_SEED = TEMPLATE_FEATURE_SEED,
  TEMPLATE_RIDGE_LAMBDA = TEMPLATE_RIDGE_LAMBDA,
  TEMPLATE_NUMERICAL_JITTER = TEMPLATE_NUMERICAL_JITTER,
  TEMPLATE_QUERY_CHUNK = TEMPLATE_QUERY_CHUNK,
  H_A_BANDWIDTH = H_A_BANDWIDTH,
  H_X_BANDWIDTH = H_X_BANDWIDTH,
  H_X2_BANDWIDTH = H_X2_BANDWIDTH,
  H_LAG_BANDWIDTH = H_LAG_BANDWIDTH,
  N_ALT_ITER_MAX = N_ALT_ITER_MAX,
  ALT_TOL = ALT_TOL,
  ALT_LOCAL_RADIUS = ALT_LOCAL_RADIUS,
  N_STARTS = N_STARTS,
  RANDOM_START_SD = RANDOM_START_SD,
  N_FOLDS = N_FOLDS,
  ANCHOR_MODE = ANCHOR_MODE,
  PROFILE_TEMPERATURE = PROFILE_TEMPERATURE,
  PROFILE_LOCAL_RADIUS = PROFILE_LOCAL_RADIUS,
  PROFILE_SEARCH_RADIUS = PROFILE_SEARCH_RADIUS,
  CLOCK_SIGNAL_MIN = CLOCK_SIGNAL_MIN,
  N_CORES = registration_cores,
  METHODS = METHODS,
  METHOD_ORDER = METHOD_ORDER,
  SCENARIOS = SCENARIOS,
  SCENARIO_LABELS = SCENARIO_LABELS,
  METHOD_LABELS = METHOD_LABELS,
  SCHEDULE_LABELS = SCHEDULE_LABELS
)

# ----------------------------- Task plan -------------------------------------
build_design_cells <- function() {
  primary <- do.call(rbind, lapply(SCENARIOS$scenario, function(scn) {
    data.frame(
      phase = "primary",
      scenario = scn,
      n_train = N_TRAIN_GRID,
      schedule = scenario_schedule(scn),
      stringsAsFactors = FALSE
    )
  }))

  timepoint <- do.call(rbind, lapply(SCENARIOS$scenario, function(scn) {
    data.frame(
      phase = "timepoint_extension",
      scenario = scn,
      n_train = TIMEPOINT_N_TRAIN,
      schedule = TIMEPOINT_SCHEDULES,
      stringsAsFactors = FALSE
    )
  }))

  design <- rbind(primary, timepoint)
  design$cell_id <- seq_len(nrow(design))
  design$n_test <- N_TEST
  design$n_rep <- N_REP
  design$n_visits <- vapply(design$schedule, function(x) length(VISIT_SCHEDULES[[x]]), integer(1))
  rownames(design) <- NULL
  design
}

build_task_plan <- function() {
  cells <- build_design_cells()
  plan <- cells[rep(seq_len(nrow(cells)), each = N_REP), , drop = FALSE]
  plan$rep <- rep(seq_len(N_REP), times = nrow(cells))
  plan$task_id <- seq_len(nrow(plan))
  rownames(plan) <- NULL
  plan
}

fixed_task_count <- function() nrow(build_task_plan())

select_task_from_env <- function() {
  task_id_source <- "SILK_TASK_ID"
  task_id_value <- Sys.getenv("SILK_TASK_ID", unset = "")
  if (!nzchar(task_id_value)) {
    task_id_source <- "SLURM_ARRAY_TASK_ID"
    task_id_value <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1")
  }
  task_id <- as.integer(task_id_value)
  plan <- build_task_plan()
  if (!is.finite(task_id) || task_id < 1L || task_id > nrow(plan)) {
    stop(
      "Invalid ", task_id_source, "=", task_id_value,
      "; expected 1..", nrow(plan),
      call. = FALSE
    )
  }
  plan[task_id, , drop = FALSE]
}
