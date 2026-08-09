# SILK: Shift-Invariant Learned Kernel for Survival Prediction

R package implementing SILK for absolute-risk prediction in the presence of
common origin-time measurement error. The registration layer uses an exact
Gaussian RBF kernel on standardized biomarker vectors. Because that biomarker
kernel is characteristic, the implemented RKHS loss retains the
distribution-identification premise of the method's theory; it is not replaced
by engineered biomarker moments or a finite biomarker feature map.

## Installation

Install directly from GitHub:

```r
# 1. Install the remotes package if you haven't already
install.packages("remotes")
# 2. Use the full HTTPS URL
remotes::install_git("https://github.com/mehrdadmhmdi/SILK.git", dependencies = TRUE)
```

## Quick Start

```r
library(SILK)

# Generate simulated data under moderate origin error
train_data <- generate_dataset_fixed(n = 400, scenario_name = "mean_moderate", seed = 1)
test_data  <- generate_dataset_fixed(n = 200, scenario_name = "mean_moderate", seed = 2)

# Fit SILK. The shift range is a candidate registration range, not an
# assumed error scenario.
fit <- fit_silk(train_data$subjects, train_data$visits,
                shift_range = c(-12, 12), seed = 1)

# Fit the registration once when several survival layers use the same shifts.
registration <- fit_silk_registration(
  train_data$subjects, train_data$visits,
  shift_range = c(-12, 12), seed = 1
)
fit <- fit_silk(
  train_data$subjects, train_data$visits,
  registration = registration
)

# Predict absolute risk at horizons 1, 2, 3, 4
predictions <- predict_silk(fit, test_data$subjects, test_data$visits)
head(predictions)

# Evaluate predictions (Brier score, AUC, calibration)
metrics <- evaluate_prediction_frame(predictions)
metrics$summary
```

## Comparators and Diagnostics

The package includes the manuscript comparators used to separate calibration
from ordinary history modeling. `fit_same_feature_recorded_cox()` is the
same-history-feature, recorded-age ablation. It is a diagnostic for what drives
the gain, not a second SILK kernel. `fit_beran_silk()` reuses the same
characteristic-kernel registration with a Beran survival layer; pass the same
`registration` object to both Cox and Beran fits to avoid refitting it.

Registration fits carry profile-gap diagnostics (`gap`, `gap_q1`, and
`gap_q2`) so weakly identified shift profiles can be flagged instead of treated
as routine successes. The MACS analysis scripts use fixed prospective
landmarks at 1, 2, 3, and 5 years after seroconversion and write
`macs_profile_diagnostics.csv` plus a summarized
`macs_profile_gap_summary.csv`.

## Simulation Scenarios

SILK ships with 12 simulation scenarios for `generate_dataset_fixed()`.
They are for simulation studies only. When fitting a real dataset, use
`shift_range` or `shift_grid`; the package does not assume that the
real data follow any named simulation error scenario.

| Scenario | Description |
|---|---|
| `no_error` | No origin error (calibration baseline) |
| `mean_moderate` | Mean-stage biomarkers, moderate error |
| `mean_severe` | Mean-stage biomarkers, severe error |
| `mean_strong_dense` | Strong mean signal, dense visits, severe error |
| `dist_moderate` | Distributional biomarkers, moderate error |
| `dist_large` | Distributional biomarkers, severe error |
| `dist_strong_dense` | Strong distributional signal, dense visits |
| `weak_stage` | Weak biomarker signal (negative control) |
| `biased_shift` | Biased origin shift (anchor misspecification) |
| `heavy_tail` | Heavy-tailed origin error |
| `asymmetric_shift` | Asymmetric origin error |
| `irregular_missing` | Irregular visits with missing data |

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

# The biomarker kernel is always the exact Gaussian RBF kernel. Its bandwidth
# may use the median-distance rule or a fixed positive value.
silk_options(BIOMARKER_BANDWIDTH = "median")

# These settings change only the scalar template-input kernel and numerical
# ridge fit; they never approximate the characteristic biomarker kernel.
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
