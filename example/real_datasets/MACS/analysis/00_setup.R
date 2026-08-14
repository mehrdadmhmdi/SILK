# =============================================================================
# 00_setup.R
# Setup: load SILK package and configure options for MACS real-data analysis
# =============================================================================

# ── Script paths ─────────────────────────────────────────────────────────────
macs_analysis_dir <- function(default = getwd()) {
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

ANALYSIS_DIR <- macs_analysis_dir()
PACKAGE_DIR  <- normalizePath(file.path(ANALYSIS_DIR, "..", "..", "..", ".."),
                               winslash = "/", mustWork = TRUE)

# ── Dependencies ─────────────────────────────────────────────────────────────
# Do not modify the user's R library during a confirmatory analysis. Fail with
# an explicit dependency list so that the software environment can be recorded
# before the analysis is run.
required_packages <- c("survival", "nlme", "ggplot2", "pkgload")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    ". Install them before running the MACS analysis.",
    call. = FALSE
  )
}

# ── Load SILK package ────────────────────────────────────────────────────────
# During development, prefer the local checkout that owns this analysis folder.
# If pkgload is unavailable, fall back to the installed SILK package.
if (file.exists(file.path(PACKAGE_DIR, "DESCRIPTION")) &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(PACKAGE_DIR, quiet = TRUE, export_all = FALSE)
} else {
  library(SILK)
}

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_DIR    <- normalizePath(file.path(ANALYSIS_DIR, "..", "data"),
                             winslash = "/", mustWork = TRUE)
RESULTS_DIR <- Sys.getenv(
  "MACS_RESULTS_DIR", unset = file.path(ANALYSIS_DIR, "results")
)
FIGURES_DIR <- Sys.getenv(
  "MACS_FIGURES_DIR", unset = file.path(ANALYSIS_DIR, "figures")
)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# ── SILK options for MACS PDS ────────────────────────────────────────────────
# SC interval: mean ~1.2 yr, range [0, 14]. Shift grid covers plausible range.
silk_options(
  DEFAULT_SHIFT_GRID_MIN = -4,
  DEFAULT_SHIFT_GRID_MAX =  4,
  SHIFT_GRID_STEP        = 0.2,
  PREDICTION_HORIZONS    = c(1, 2, 3, 5),
  EVAL_HORIZON           = 3.0,
  N_ALT_ITER_MAX         = 7L,
  ALT_TOL                = 1e-3,
  N_STARTS               = 3L,
  RANDOM_START_SD        = 0.5,
  N_FOLDS                = as.integer(Sys.getenv("MACS_INNER_FOLDS", unset = "5")),
  ANCHOR_MODE            = "rx",
  # Parallel cross-fit folds. "auto" uses all available/allocated cores
  # (honours SILK_N_CORES, then SLURM_CPUS_PER_TASK, then detectCores()).
  # Set N_CORES = 1 or SILK_N_CORES=1 to disable. Results are identical to
  # serial because each fold is self-seeded.
  N_CORES                = "auto",
  # Prespecified stabilization. This value is not selected on the outer
  # held-out predictions.
  SHIFT_RIDGE            = as.numeric(Sys.getenv("SILK_SHIFT_RIDGE", unset = "1e-3")),
  H_A_BANDWIDTH          = 1.5,
  H_X_BANDWIDTH          = 2.0,
  BIOMARKER_KERNEL        = "gaussian",
  BIOMARKER_BANDWIDTH     = "median",
  BIOMARKER_BANDWIDTH_SCALE = 1,
  SURVIVAL_HISTORY_BIOMARKERS = 3L,
  PROFILE_TEMPERATURE    = 0.015,
  PROFILE_LOCAL_RADIUS   = 0.8,
  UIUC_ORANGE            = "#FF5F05",
  UIUC_BLUE              = "#13294B",
  FIGURE_DPI             = 600,
  FIGURE_BASE_SIZE       = 14,
  METHOD_ORDER = c(
    "SILK-Cox", "Cox-Recorded-SameFeature", "MMLM-Recorded"
  ),
  METHOD_LABELS = c(
    "SILK-Cox"                 = "SILK-Cox (Gaussian)",
    "Cox-Recorded-SameFeature" = "Recorded-clock Cox",
    "MMLM-Recorded"            = "Parametric MMLM-Cox"
  )
)

# ── Real-data registration grid ──────────────────────────────────────────────
# This is only a candidate search range. It is not an assumed error scenario.
MACS_SHIFT_RANGE <- c(-4, 4)
landmark_override <- trimws(Sys.getenv("MACS_LANDMARKS", unset = ""))
MACS_FIXED_LANDMARKS <- if (nzchar(landmark_override)) {
  as.numeric(trimws(strsplit(landmark_override, ",", fixed = TRUE)[[1]]))
} else {
  c(1, 2, 3, 5)
}
if (any(!is.finite(MACS_FIXED_LANDMARKS)) || any(MACS_FIXED_LANDMARKS <= 0)) {
  stop("MACS_LANDMARKS must be a comma-separated list of positive numbers.",
       call. = FALSE)
}
# The primary analysis uses the three complete immunologic markers. Viral load
# remains in the processed file but is excluded because its historical
# subject-median imputation can use measurements obtained after a landmark.
MACS_PRIMARY_BIOMARKERS <- c("B1", "B2", "B3")
MACS_USE_FIXED_LANDMARKS <- !identical(
  tolower(Sys.getenv("MACS_USE_FIXED_LANDMARKS", unset = "true")),
  "false"
)

# ── Column-name mapping: MACS originals → SILK internals ────────────────────
# The CSVs use descriptive MACS names; SILK functions expect generic names.
# These helpers handle the renaming.

macs_to_silk_subjects <- function(df) {
  names(df)[names(df) == "CASEID"]      <- "id"
  names(df)[names(df) == "AGE_SC"]      <- "X1"
  names(df)[names(df) == "NONWHITE"]    <- "X2"
  names(df)[names(df) == "DISEASE_AGE"] <- "A_obs"
  names(df)[names(df) == "FOLLOWUP"]    <- "U"
  names(df)[names(df) == "AIDS_EVENT"]  <- "delta"
  df
}

macs_to_silk_visits <- function(df) {
  names(df)[names(df) == "CASEID"]      <- "id"
  names(df)[names(df) == "VISIT"]       <- "visit"
  names(df)[names(df) == "LAG"]         <- "lag"
  names(df)[names(df) == "DISEASE_AGE"] <- "A_obs_il"
  names(df)[names(df) == "AGE_SC"]      <- "X1"
  names(df)[names(df) == "NONWHITE"]    <- "X2"
  names(df)[names(df) == "LEU3N"]       <- "B1"
  names(df)[names(df) == "LEU3P"]       <- "B2"
  names(df)[names(df) == "LEU2N"]       <- "B3"
  names(df)[names(df) == "LOG_VLOAD"]   <- "B4"
  df
}

macs_fixed_landmark_data <- function(subjects, visits, landmark_time) {
  landmark_time <- as.numeric(landmark_time)[1]
  total_time <- subjects$A_obs + subjects$U
  eligible <- is.finite(total_time) & total_time > landmark_time
  s <- subjects[eligible, , drop = FALSE]
  total_time <- total_time[eligible]
  if (!nrow(s)) {
    return(list(subjects = s, visits = visits[FALSE, , drop = FALSE]))
  }
  s$U <- total_time - landmark_time
  s$delta <- as.integer(s$delta == 1L)
  s$A_obs <- landmark_time
  v <- visits[visits$id %in% s$id & visits$A_obs_il <= landmark_time + 1e-8, , drop = FALSE]
  if (nrow(v)) {
    v$lag <- landmark_time - v$A_obs_il
  }
  keep_ids <- intersect(s$id, unique(v$id))
  s <- s[s$id %in% keep_ids, , drop = FALSE]
  v <- v[v$id %in% keep_ids, , drop = FALSE]
  list(subjects = s, visits = v)
}

macs_landmark_sets <- function(subjects, visits) {
  if (!isTRUE(MACS_USE_FIXED_LANDMARKS)) {
    return(list(last_visit = list(subjects = subjects, visits = visits, landmark_time = NA_real_)))
  }
  out <- list()
  for (tt in MACS_FIXED_LANDMARKS) {
    nm <- paste0("fixed_", tt, "y")
    z <- macs_fixed_landmark_data(subjects, visits, tt)
    z$landmark_time <- tt
    out[[nm]] <- z
  }
  out
}

cat("Setup complete.\n")
cat("  Package dir:", PACKAGE_DIR, "\n")
cat("  Analysis dir:", ANALYSIS_DIR, "\n")
cat("  Data dir:   ", DATA_DIR, "\n")
cat("  Results dir:", RESULTS_DIR, "\n")
cat("  Figures dir:", FIGURES_DIR, "\n")
