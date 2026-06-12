# =============================================================================
# 00_setup.R
# Setup: load the SILK package and configure options for MACS real-data analysis
# =============================================================================

# ── Dependencies ─────────────────────────────────────────────────────────────
required_packages <- c("survival", "nlme", "ggplot2", "devtools")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# ── Load SILK package from local source ──────────────────────────────────────
SILK_PKG_DIR <- file.path(
  dirname(getwd()),  # adjust if needed
  "simulation", "main_simulation"
)
if (!file.exists(file.path(SILK_PKG_DIR, "DESCRIPTION"))) {
  # Try relative to this script's location
  SILK_PKG_DIR <- normalizePath(file.path("..", "simulation", "main_simulation"))
}
cat("Loading SILK from:", SILK_PKG_DIR, "\n")
devtools::load_all(SILK_PKG_DIR, quiet = TRUE)

# ── Configure SILK options for MACS PDS ──────────────────────────────────────
# The seroconversion interval in MACS has mean ~1.2 yr, range [0, 14].
# The common origin shift eps_i is bounded by this interval.
# We set a moderate shift grid: [-4, 4] with step 0.2.

silk_options(
  # Shift grid — covers realistic SC-interval uncertainty
  DEFAULT_SHIFT_GRID_MIN = -4,
  DEFAULT_SHIFT_GRID_MAX =  4,
  SHIFT_GRID_STEP        = 0.2,

  # Prediction horizons (years from landmark)
  PREDICTION_HORIZONS = c(1, 2, 3, 5),
  EVAL_HORIZON        = 3.0,

  # Registration
  N_ALT_ITER_MAX  = 7L,
  ALT_TOL         = 1e-3,
  N_STARTS        = 3L,
  RANDOM_START_SD = 0.5,
  N_FOLDS         = 5L,
  ANCHOR_MODE     = "rx",

  # Bandwidths — wider for real data (more heterogeneity)
  H_A_BANDWIDTH = 1.5,
  H_X_BANDWIDTH = 2.0,

  # Biomarker settings (4 biomarkers: CD4, CD4%, CD8, log10VL)
  SURVIVAL_HISTORY_BIOMARKERS = 3L,

  # Profile refinement
  PROFILE_TEMPERATURE  = 0.015,
  PROFILE_LOCAL_RADIUS = 0.8,

  # Plotting
  UIUC_ORANGE = "#FF5F05",
  UIUC_BLUE   = "#13294B",
  FIGURE_DPI  = 600,
  FIGURE_BASE_SIZE = 14,

  # Method labels for real-data analysis
  METHOD_ORDER = c(
    "Landmark-Recorded", "MMLM-Recorded",
    "SILK"
  ),
  METHOD_LABELS = c(
    "Landmark-Recorded" = "Landmark Cox",
    "MMLM-Recorded"     = "Mixed-Model Landmark",
    "SILK"              = "SILK"
  )
)

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_DIR    <- normalizePath(file.path("..", "real_datasets"))
RESULTS_DIR <- file.path(getwd(), "results")
FIGURES_DIR <- file.path(getwd(), "figures")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Setup complete.\n")
cat("  Data dir:   ", DATA_DIR, "\n")
cat("  Results dir:", RESULTS_DIR, "\n")
cat("  Figures dir:", FIGURES_DIR, "\n")
