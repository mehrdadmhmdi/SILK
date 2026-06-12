# SILK: Shift-Invariant Learned Kernel for Survival Prediction

R package implementing the SILK method for absolute risk prediction in the presence of common origin-time measurement error. SILK uses biomarker trajectory registration with shift-invariant kernel features to calibrate landmark ages and improve survival prediction.

## Installation

Install directly from GitHub:

```r
# 1. Install the remotes package if you haven't already
install.packages("remotes")
# 2. Use the full HTTPS URL
remotes::install_git("https://github.com/mehrdadmhmdi/SILK.git")
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

# Fast large-sample fit with random Fourier features
fit_fast <- fit_silk(train_data$subjects, train_data$visits,
                     shift_range = c(-12, 12),
                     kernel = "rbf", kernel_approx = "rff",
                     rff_dim = 512, seed = 1)

# Predict absolute risk at horizons 1, 2, 3, 4
predictions <- predict_silk(fit, test_data$subjects, test_data$visits)
head(predictions)

# Evaluate predictions (Brier score, AUC, calibration)
metrics <- evaluate_prediction_frame(predictions)
metrics$summary
```

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

# Choose the registration kernel and optional RFF acceleration.
# Kernels: "rbf", "matern", "polynomial", "linear".
silk_options(REGISTRATION_KERNEL = "matern",
             REGISTRATION_KERNEL_APPROX = "rff",
             KERNEL_RFF_DIM = 512,
             KERNEL_MATERN_NU = 1.5)

# See all defaults
str(silk_default_options())
```

## Dependencies

**Required:** survival, stats, rlang

**Optional:** ggplot2 (visualization), nlme (mixed-model biomarker prediction)

## License

MIT
