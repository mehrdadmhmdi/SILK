# MACS PDS Real-Data Application

This folder contains a real-data application of the SILK method to the
**Multicenter AIDS Cohort Study (MACS) Public Data Set (PDS)**.

## Study Overview

The MACS PDS tracks HIV-positive men who seroconverted during follow-up.
The goal is to predict progression to AIDS from multivariate longitudinal
immunologic trajectories when the true time of HIV infection is only known to
lie between the last negative and first positive test.

The primary model, **SILK-Cox (Gaussian)**, nonparametrically registers the
longitudinal histories to a common disease clock and uses the calibrated stage
in a residual-time Cox model.

## Cohort

- **573 HIV seroconverters** with ≥ 2 post-seroconversion visits with CD4 data
- **10,602 longitudinal visits**
- **228 AIDS events** (39.8% event rate)
- Seroconversion midpoint used as estimated disease origin

## Data (`data/`)

| File | Rows | Description |
|------|------|-------------|
| `macs_subjects.csv` | 573 | One row per subject — baseline covariates and outcome |
| `macs_visits.csv` | 10,602 | Longitudinal visit-level biomarker measurements |

### `macs_subjects.csv` columns

| Column | Description |
|--------|-------------|
| `CASEID` | Subject identifier |
| `AGE_SC` | Age at seroconversion (years) |
| `NONWHITE` | Race indicator (1 = non-white, 0 = white) |
| `DISEASE_AGE` | Time from seroconversion midpoint to landmark (years) |
| `FOLLOWUP` | Residual follow-up time from landmark to event/censoring (years) |
| `AIDS_EVENT` | AIDS event indicator (1 = AIDS, 0 = censored) |

### `macs_visits.csv` columns

| Column | Description |
|--------|-------------|
| `CASEID` | Subject identifier |
| `VISIT` | Visit sequence number (1, 2, ...) |
| `LAG` | Time from visit to landmark (years, counting backward) |
| `DISEASE_AGE` | Disease age at this visit (years since SC midpoint) |
| `AGE_SC` | Age at seroconversion (years; baseline covariate) |
| `NONWHITE` | Race indicator (baseline covariate) |
| `LEU3N` | CD4+ T-cell count (cells/µL) |
| `LEU3P` | CD4+ T-cell percentage (%) |
| `LEU2N` | CD8+ T-cell count (cells/µL) |
| `LOG_VLOAD` | log₁₀(HIV viral load) — imputed with subject median where missing |

## Analysis (`analysis/`)

Three scripts define the final analysis pipeline. They resolve paths relative to
the script location, so they can be run from the package root or from
`analysis/`.

| Script | Purpose |
|--------|---------|
| `00_setup.R` | Load SILK package, configure options, define column mappings |
| `01_analysis.R` | Subject-level 5-fold cross-validated prediction for all prespecified methods |
| `02_results.R` | Landmark-specific IPCW metrics, complete comparator tables, and figures |

### Running

```r
source("example/real_datasets/MACS/analysis/01_analysis.R")
source("example/real_datasets/MACS/analysis/02_results.R")
```

The setup script loads the local SILK checkout through `pkgload::load_all()`.
SILK uses the Gaussian biomarker kernel with a median-distance bandwidth and a
prespecified candidate registration range of `[-4, 4]` years. The shift is a
registration score; it is not treated as an observed infection-time error.

The outer folds are assigned once at the subject level and reused at every
landmark and for every method. All registration, smoothing, standardization,
and model fitting are performed with training-fold data. The result script
estimates the censoring distribution separately within each landmark risk set.
It reports every prespecified landmark-horizon cell and will not replace the
final prediction file if any method fails.

### Methods compared

1. **SILK-Cox (Gaussian)** — nonparametric multivariate registration of CD4,
   CD4%, and CD8 histories, followed by residual-time Cox prediction.
2. **Recorded-clock Cox** — the primary comparator, with the same Cox
   predictors as SILK-Cox but the recorded rather than calibrated disease
   clock. This comparison isolates the contribution of registration.
3. **Parametric MMLM-Cox** — separate linear mixed models smooth the three
   biomarker histories; their current values enter a residual-time Cox model.

The primary scientific contrast is SILK-Cox versus Recorded-clock Cox.
MMLM-Cox supplies a focused parametric comparison without adding a collection
of overlapping joint-model variants.

### Primary biomarkers

The final analysis uses CD4 count (`LEU3N`), CD4 percentage (`LEU3P`), and CD8
count (`LEU2N`). Viral load remains in the processed file but is excluded from
the primary analysis because its historical subject-median imputation can use
measurements obtained after an early landmark.

### Landmarks and horizons

Risk sets are constructed at 1, 2, 3, and 5 years after the recorded
seroconversion midpoint. Prediction horizons are 1, 2, 3, and 5 years from
each landmark.

### Outputs

`macs_comparator_table.csv` contains the complete method comparison.
`macs_primary_silk_comparison.csv` gives the prespecified SILK-Cox versus
Recorded-clock Cox contrast; positive AUC and Brier gains favor SILK. Figures
use Illini Blue (`#13294B`) for SILK-Cox and Illini Orange (`#FF5F05`) for the
MMLM-Cox comparator. The single-method SILK risk-stratification figure uses
blue line types only. No age-correction diagnostic is generated.
Metrics are not averaged across landmark-horizon cells because those cells
represent different risk sets and prediction tasks.

## Notes

- The `raw_data/` directory (proprietary MACS PDS fixed-width ASCII files)
  is not included in the repository. The processed CSVs in `data/` contain
  all information needed to reproduce the analysis.
- Viral load (`LOG_VLOAD`) was not measured in early MACS visits (pre-1996).
  The processed file contains historical subject-median imputations, but viral
  load is not used by the final primary models.
- Censoring time is derived from the MACS administrative report date
  (`REPDATE` in `outcome.dat`), not from the last visit date.

## Citation

If you use this data application, please cite both the SILK method and the
MACS PDS:

> Jacobson LP, et al. "The Multicenter AIDS Cohort Study (MACS) Public
> Dataset." *NIH/NIDA*, 2007.
