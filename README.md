# SILK: Shift-Invariant Learned Kernel for Survival Prediction

R package implementing SILK for absolute-risk prediction in the presence of
common origin-time measurement error. The registration layer uses an exact
characteristic kernel on standardized biomarker vectors: Gaussian RBF by
default, with Laplace and Matern-3/2 alternatives. The implemented RKHS loss is
not replaced by engineered biomarker moments or a finite biomarker feature map.

## Installation

Install directly from GitHub:

```r
# 1. Install the remotes package if you haven't already
install.packages("remotes")
# 2. Use the full HTTPS URL
remotes::install_git("https://github.com/mehrdadmhmdi/SILK.git", dependencies = TRUE)
```

## Data layout

SILK accepts one subject-level table and one longitudinal visit table. For
backward compatibility, the canonical columns are `id`, `A_obs`, `U`,
`delta`, `visit`, `lag`, `A_obs_il`, and `B1`, `B2`, .... For an R-ready
application, use `silk_data_spec()` instead: it maps arbitrary source column
names, selects the biomarker and survival covariates, and accepts named
functions for engineered covariates. The package does not assume or create
`X3^2 - 1`; that feature is created only when the user requests it.

## Analysis example

```r
library(SILK)

# Replace these paths with the prepared analysis tables.
train_subjects <- read.csv("train_subjects.csv")
train_visits   <- read.csv("train_visits.csv")
test_subjects  <- read.csv("test_subjects.csv")
test_visits    <- read.csv("test_visits.csv")

spec <- silk_data_spec(
  subject_id = "subject_id",
  recorded_age = "age_recorded",
  event_time = "follow_up",
  event = "status",
  visit_number = "visit_no",
  visit_age = "age_at_visit",
  visit_lag = NULL,
  biomarker_cols = c("cd4", "cd8"),
  covariate_cols = c("sex", "baseline_age"),
  engineered_covariates = list(
    baseline_age_sq = function(data) data[["baseline_age"]]^2 - 1
  ),
  template_input_covariates = c("sex", "baseline_age_sq", "lag"),
  anchor_mode = "intercept"
)

# Fit the registration once. The shift range is the scientifically plausible
# range of origin-time offsets in the application.
registration <- fit_silk_registration(
  train_subjects, train_visits,
  shift_range = c(-12, 12), seed = 1, data_spec = spec
)

# Inspect cross-fitted training shifts and profile diagnostics before fitting
# the survival layer.
head(registration$train_stage)

fit <- fit_silk(
  train_subjects, train_visits,
  registration = registration, data_spec = spec
)

# Predict absolute risk at the requested horizons.
predictions <- predict_silk(
  fit, test_subjects, test_visits,
  horizons = c(1, 2, 3, 4)
)
head(predictions)

# If outcome columns are available in the test subject table, evaluate Brier
# score, AUC, and calibration using the common prediction-frame interface.
metrics <- evaluate_prediction_frame(predictions)
metrics$summary
```

## Registration diagnostics

`fit_silk_registration()` returns cross-fitted training shifts and the final
template used for new subjects. Review boundary indicators and the profile gaps
`gap_q1` and `gap_q2`; small gaps indicate that materially different shifts
have nearly the same registration loss and should not be treated as strongly
identified estimates.

For a new dataset, `predict_silk_registration()` returns the estimated shift,
the biomarker-clock initializer, the grid minimum, boundary status, and profile
gaps without fitting or predicting the survival outcome.

## Visualization

With `ggplot2` installed:

```r
library(ggplot2)

# Publication-quality theme
p <- ggplot(predictions, aes(x = risk_pred)) +
  geom_histogram(bins = 30) +
  theme_silk()

# Method comparison palette
method_palette()
```

## Configuration

All tuning parameters are adjustable via `silk_options()`:

```r
# View current settings
silk_opt("SHIFT_GRID_STEP")

# Change registration grid resolution
silk_options(SHIFT_GRID_STEP = 0.1, N_FOLDS = 10)

# Select an exact characteristic biomarker kernel and its bandwidth rule.
silk_options(BIOMARKER_KERNEL = "gaussian")
silk_options(BIOMARKER_BANDWIDTH = "median")

# These settings control the conditional template-input feature map and its
# numerical ridge fit; they do not approximate the biomarker kernel itself.
silk_options(TEMPLATE_INPUT_FEATURES = 32,
             TEMPLATE_RIDGE_LAMBDA = 0.001,
             N_CORES = 4)

# See all defaults
str(silk_default_options())
```

## Dependencies

**Required:** survival, stats, rlang

**Optional:** ggplot2 (visualization), nlme (mixed-model biomarker prediction)

## License

MIT
