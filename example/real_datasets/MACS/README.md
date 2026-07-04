# MACS PDS Real-Data Application

This folder contains a real-data application of the SILK method to the
**Multicenter AIDS Cohort Study (MACS) Public Data Set (PDS)**.

## Study Overview

The MACS PDS tracks HIV-positive men who seroconverted during follow-up.
The goal is to predict progression to AIDS using longitudinal biomarker
trajectories, where the true time of HIV infection (disease origin) is
interval-censored between the last negative and first positive test.

SILK addresses this by estimating subject-specific origin shifts, realigning
biomarker curves to a common disease clock before building a landmark
prediction model.

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

Three scripts run the full analysis pipeline. They resolve paths relative to
the script location, so they can be run from the package root or from
`analysis/`.

| Script | Purpose |
|--------|---------|
| `00_setup.R` | Load SILK package, configure options, define column mappings |
| `01_analysis.R` | 5-fold cross-validated risk prediction (Landmark Cox, MMLM, SILK) |
| `02_results.R` | Compute metrics (AUC, Brier, calibration) and produce figures |

### Running

```r
source("real_datasets/MACS/analysis/01_analysis.R")   # writes analysis/results/
source("real_datasets/MACS/analysis/02_results.R")    # writes analysis/figures/
```

The setup script prefers the local SILK checkout through `pkgload::load_all()`
when `pkgload` is installed, otherwise it falls back to the installed `SILK`
package. The MACS setup uses the RBF registration kernel with random Fourier
features (`REGISTRATION_KERNEL_APPROX = "rff"`, `KERNEL_RFF_DIM = 512`) and
fits SILK over the candidate shift range `MACS_SHIFT_RANGE = c(-4, 4)`. This
range is not a named error scenario; the real-data error distribution is
unknown and is represented only through the input data. Generated outputs
(`results/`, `figures/`) are gitignored.

### Methods compared

1. **Landmark Cox** — Standard landmark model using recorded visit times
2. **Mixed-Model Landmark (MMLM)** — Current-value joint model with recorded clock
3. **SILK** — Shift-Invariant Landmark Kernels (nonparametric registration + landmark prediction)

### Prediction horizons

1, 2, 3, and 5 years from landmark.

## Notes

- The `raw_data/` directory (proprietary MACS PDS fixed-width ASCII files)
  is not included in the repository. The processed CSVs in `data/` contain
  all information needed to reproduce the analysis.
- Viral load (`LOG_VLOAD`) was not measured in early MACS visits (pre-1996).
  Missing values (1,386 of 10,602 visits) were imputed using the
  subject-specific median; the global median was used when all visits for a
  subject were missing.
- Censoring time is derived from the MACS administrative report date
  (`REPDATE` in `outcome.dat`), not from the last visit date.

## Citation

If you use this data application, please cite both the SILK method and the
MACS PDS:

> Jacobson LP, et al. "The Multicenter AIDS Cohort Study (MACS) Public
> Dataset." *NIH/NIDA*, 2007.
