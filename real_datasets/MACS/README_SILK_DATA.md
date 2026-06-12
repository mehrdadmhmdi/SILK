# MACS PDS → SILK Data Mapping

## Source
MACS Public Data Set Release P27 (Oct 2019), 602 HIV seroconverters (STATUS=4).
After filtering for ≥2 post-seroconversion visits with CD4 data: **573 subjects, 10,602 visit records**.

## Files

| File | Rows | Description |
|------|------|-------------|
| `macs_subjects.csv` | 573 | Subject-level (one row per seroconverter) |
| `macs_visits.csv` | 10,602 | Visit-level longitudinal biomarkers |

## Column Mapping

### `macs_subjects.csv`

| SILK Column | MACS Source | Construction |
|------------|-------------|--------------|
| `id` | CASEID | Integer subject ID |
| `X1` | BIRTHDTY, V1DATY, V2DATY | Age at seroconversion = SC_midpoint − birth_year |
| `X2` | RACE (section2.dat) | Binary: 0 = White non-Hispanic, 1 = non-White |
| `A_obs` | LDATY, SC_midpoint | Disease age at landmark = last_visit_year − SC_midpoint |
| `U` | DATE1yy/REPDATEyy | Residual time from landmark to event/censoring |
| `delta` | AIDSCASE | 1 = CDC AIDS (AIDSCASE ∈ {2,3}), 0 = censored |

### `macs_visits.csv`

| SILK Column | MACS Source | Construction |
|------------|-------------|--------------|
| `id` | CASEID | Subject ID (links to subjects) |
| `visit` | sequential | Visit number within subject (1, 2, ...) |
| `lag` | derived | A_obs − A_obs_il (time from this visit to landmark) |
| `A_obs_il` | LDATY, SC_midpoint | Disease age at visit = visit_year − SC_midpoint |
| `X1` | same as subjects | Age at seroconversion |
| `X2` | same as subjects | Race indicator |
| `B1` | LEU3N | CD4+ T-cell count (cells/μl) |
| `B2` | LEU3P ÷ 10 | CD4+ percentage |
| `B3` | LEU2N | CD8+ T-cell count (cells/μl) |
| `B4` | log₁₀(VLOAD+1) | Log-transformed HIV viral load |

## Key Design Decisions

1. **Seroconversion midpoint** as estimated onset: SC_mid = (last_negative_date + first_positive_date) / 2. This is the "observed origin" A_obs, subject to common origin shift ε_i.

2. **Landmark** = each subject's last observed visit (variable landmark, matching SILK simulation where A_star varies across subjects).

3. **Residual time U**:
   - Events (delta=1): U = AIDS_diagnosis_date − last_visit_date
   - Censored (delta=0): U = REPDATE − last_visit_date (or 0.5 years if REPDATE unavailable)

4. **Missing viral load** (B4): 13% of visits lack VLOAD (not tested in early MACS years). Imputed with subject-specific median; if all missing for a subject, global median used.

5. **Seroconversion interval** (the measurement error): mean = 1.2 years, range [0, 14]. This is the source of the common origin shift that SILK is designed to handle.

## Usage in R

```r
subjects <- read.csv("macs_subjects.csv")
visits   <- read.csv("macs_visits.csv")

# These match the output format of silk::generate_dataset_fixed()
# Can be used directly with:
#   fit  <- silk::fit_silk(subjects, visits, scenario_name = "mean_moderate")
#   pred <- silk::predict_silk(fit, test_subjects, test_visits)
```

## Summary Statistics

| Variable | Mean | Median | Range |
|----------|------|--------|-------|
| N subjects | 573 | — | — |
| Events (AIDS) | 228 (39.8%) | — | — |
| A_obs (disease age at landmark) | 10.5 yr | 7.0 yr | [0, 32.5] |
| U (residual follow-up) | 1.4 yr | 0.5 yr | [0, 22.5] |
| Visits per subject | 18.5 | 13 | [2, 67] |
| B1 (CD4) | 588 | — | [0, 3015] |
| B2 (CD4%) | 26.1% | — | [0, 74] |
| B3 (CD8) | 885 | — | [14, 3627] |
| B4 (log₁₀ VL) | 3.0 | — | [0.3, 6.8] |
| X1 (age at SC) | 36.7 yr | — | [18.5, 70.5] |
| X2 (non-White) | 106 (18.5%) | — | — |
