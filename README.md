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

SILK expects one subject-level table and one longitudinal visit table. The
subject table must contain `id`, recorded landmark age `A_obs`, follow-up time
`U`, event indicator `delta`, and the baseline analysis variables `X1`--`X4`.
The visit table must contain `id`, nominal visit position `visit`, time before
the landmark `lag`, recorded visit age `A_obs_il`, the same baseline variables,
and biomarker columns named `B1`, `B2`, and so on.

The package derives the centered nonlinear term `X3^2 - 1` internally for the
survival layer. Rename or construct these analysis columns before fitting if
the source dataset uses different names.

## Analysis example

```r
library(SILK)

# Replace these paths with the prepared analysis tables.
train_subjects <- read.csv("train_subjects.csv")
train_visits   <- read.csv("train_visits.csv")
test_subjects  <- read.csv("test_subjects.csv")
test_visits    <- read.csv("test_visits.csv")

# Fit the registration once. The shift range is the scientifically plausible
# range of origin-time offsets in the application.
registration <- fit_silk_registration(
  train_subjects, train_visits,
  shift_range = c(-12, 12), seed = 1
)

# Inspect cross-fitted training shifts and profile diagnostics before fitting
# the survival layer.
head(registration$train_stage)

fit <- fit_silk(
  train_subjects, train_visits,
  registration = registration
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
