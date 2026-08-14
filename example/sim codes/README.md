# SILK simulation pipeline

This directory is the cluster-ready implementation of the confirmatory SILK
simulation. The manuscript is not read or modified by these scripts.

## Scientific target

- Event times follow a Gompertz baseline hazard on attained latent age.
- The primary age comparison contrasts the age-only recorded and oracle Cox
  bases with SILK-Gaussian, whose calibrated age is learned from longitudinal
  biomarkers.
- Biomarkers enter the primary SILK-Cox analysis only through the registration
  layer. Recorded-clock Cox therefore cannot bypass registration using a
  hand-built biomarker summary.
- `weak_stage` is a prespecified stage-null negative control: biomarkers contain
  neither latent-age nor frailty information. Its expected result is no
  systematic registration gain; a guaranteed loss is not imposed.
- Gaussian RBF is the primary exact characteristic kernel. Exact Laplace and
  Matern-3/2 kernels are paired robustness analyses using the same data, folds,
  and Monte Carlo seeds. Linear and polynomial biomarker kernels are rejected.
- RSF and DeepSurv receive recorded age, baseline covariates, and the current
  observed biomarker.
  DeepSurv uses a fixed CPU network, a training-only validation split, early
  stopping, and no Python data-loader workers.

The simulation writes registration diagnostics for both cross-fitted training
subjects and held-out test subjects. Inspect stage RMSE and stage R-squared—not
shift correlation alone—before interpreting prediction gains. When origin
error is much more variable than latent age, even a nearly constant stage
estimate can produce a deceptively high correlation between estimated and true
shifts.

### Variables used by SILK

- The exact characteristic biomarker kernel uses all generated biomarker
  channels, `B1`--`B4`, jointly. No biomarker channel is dropped.
- Its biology-clock initializer uses one kernel on the complete longitudinal
  history (four biomarkers at every available nominal visit position), after a
  training-only ridge residualization on `X1`, `X2`, `X3`, `X4`, and centered
  `X3^2`. The neighbourhood size is selected without survival outcomes by
  leave-one-out prediction of recorded landmark age.
- The conditional RKHS template takes candidate latent age, `X1`--`X4`, and
  visit lag as inputs. Its Gaussian input kernel can represent nonlinear
  functions of `X3`; the explicit centered `X3^2` term is additionally used in
  the clock residualization and survival layer.
- The SILK survival layer uses calibrated attained age together with `X1`,
  `X2`, `X3`, `X4`, and centered `X3^2`. Thus SILK is not assigned a
  misspecified event model in this design.

The profile optimization is a local half-year refinement of the full-history
clock. This is a fixed estimator choice in the locked configuration, not an
outcome-dependent selection performed by an array task.

## Reproduce one simulation analysis in R

The package README is intentionally limited to package installation and the
analysis workflow for prepared data. Simulation-only data generation, scenario
definitions, comparator rosters, and cluster instructions are documented here.

```r
library(SILK)

train <- generate_dataset_fixed(
  n = 400, scenario_name = "mean_moderate", schedule_name = "m4", seed = 1
)
test <- generate_dataset_fixed(
  n = 1000, scenario_name = "mean_moderate", schedule_name = "m4", seed = 2
)

grid <- make_shift_grid("mean_moderate")
registration <- fit_silk_registration(
  train$subjects, train$visits,
  shift_grid = grid, seed = 1, biomarker_kernel = "gaussian"
)
fit <- fit_silk(
  train$subjects, train$visits,
  registration = registration,
  method = "SILK-Gaussian"
)
pred <- predict_silk(
  fit, test$subjects, test$visits,
  horizons = c(1, 2, 3, 4)
)
evaluate_prediction_frame(pred)$summary
```

### Simulation scenarios

| Scenario | Description |
|---|---|
| `no_error` | No origin error (calibration baseline) |
| `mean_moderate` | Mean-stage biomarkers, moderate error |
| `mean_severe` | Mean-stage biomarkers, severe error |
| `mean_strong_dense` | Strong mean signal, dense visits, severe error |
| `dist_moderate` | Distributional biomarkers, moderate error |
| `dist_large` | Distributional biomarkers, severe error |
| `dist_strong_dense` | Strong distributional signal, dense visits |
| `weak_stage` | Stage-null biomarker negative control |
| `biased_shift` | Biased origin shift (anchor misspecification) |
| `heavy_tail` | Heavy-tailed origin error |
| `asymmetric_shift` | Asymmetric origin error |
| `irregular_missing` | Irregular visits with missing data |

The named scenarios are data-generating mechanisms only. A real-data SILK fit
uses a scientifically plausible `shift_range` or explicit `shift_grid`; it does
not assume one of these scenario names.

## Cluster run

The controller installs SILK from
`https://github.com/mehrdadmhmdi/SILK.git`. Commit and push the intended package
changes before submission, or set `SILK_PACKAGE_REF` to the exact pushed commit.

```bash
cd "/path/to/SILK/git/example/sim codes"
sbatch job_validate.sbatch
```

`job_validate.sbatch` runs the complete design once (`SILK_N_REP=1`) and then
checks collection, every requested method, and the two primary reports. After
its dependent `silk-collect` job succeeds, launch the fixed 50-replication run:

```bash
sbatch job_confirmatory.sbatch
```

The confirmatory launcher refuses to run unless validation used the same
`locked_confirmatory_config.sh` file and the same installed SILK package commit.
If all validation tasks exist but its collection/reporting job stopped before
writing the completion marker, the launcher first rebuilds those summaries from
the existing raw files and proceeds only after that recheck succeeds. Do not
submit `job_full.sbatch` directly; it is the internal array controller used by
the two stage-specific launchers.
All hyperparameters are fixed in that file from a separate development design
not included in the confirmatory claims; neither launcher tunes against the
confirmatory outcomes.

### Small n=400, R=2, four-visit test

To test the complete 14-method roster across all 12 error scenarios using only
the four-visit (`m4`) design, submit:

```bash
sbatch job_test_n400_r2_m4.sbatch
```

This creates 24 simulation tasks (12 scenarios times 2 replications), all with
`n_train=400` and `n_visits=4`. It writes to
`outputs-test-n400-r2-m4/`, keeps the standard prediction horizons 1--4, runs
the collection and two-question report after the arrays finish, and does not
submit manuscript figures.

The first controller run creates a persistent CPU environment under
`.cluster-env-r44/`, installs the R dependencies and CPU-only `pycox` stack,
and installs the requested SILK GitHub commit. Later runs reuse that environment
and reinstall SILK only when the requested commit changes. Array workers only
load the prepared environment; they never clone or install packages.

Every Slurm job uses:

```bash
#SBATCH --account=bfez-delta-cpu
#SBATCH --partition=cpu

module purge
module load slurm-env/0.1
module load cray-R/4.4.0
```

Top-level submission logs are always written as `logs/<job-name>-<job-id>.log`.
Array, collection, audit, and figure jobs write `.log` files in the `logs/`
subdirectory of their run-specific output directory. The repository tracks the
top-level `logs/` directory so Slurm can open the launcher log at job start.

Arrays are submitted without a user-imposed concurrency throttle. Site/QOS
limits still apply.

The package reference and the two output locations may be changed without
changing the locked statistical configuration:

```bash
export SILK_PACKAGE_REF="<pushed-commit-or-tag>"
export SILK_VALIDATION_OUT_DIR="$PWD/outputs-validation-nrep1"
export SILK_CONFIRMATORY_OUT_DIR="$PWD/outputs-confirmatory-50rep"
sbatch job_validate.sbatch
```

### Exact method roster

The simulation accepts exactly these 14 method IDs:

| Role | Methods | Longitudinal biomarkers in survival prediction? |
|---|---|---|
| Recorded-age bases | `Recorded-Cox`, `Recorded-Beran` | No |
| SILK registration kernels | `SILK-Gaussian`, `SILK-Laplace`, `SILK-Matern32` | Yes, for registration |
| Biomarker comparators | `MMLM-Correct`, `MMLM-Misspecified`, `JM-Correct`, `JM-Misspecified`, `RSF`, `DeepSurv`, `TimeError-Integrated-Landmark` | Yes |
| Latent-age oracle bases | `Oracle-Cox`, `Oracle-Beran` | No |

The four biomarker-free methods use only their stated age coordinate; no
baseline covariate, current biomarker, biomarker history, or biomarker-derived
feature enters those fits. The three SILK methods share
the same survival layer, folds, seeds, and tuning, so the registration kernel
is the only difference among them. No rich-history Cox, same-feature recorded
Cox, SILK mean-regression, or extra Beran variant remains in the active design.
`MMLM-Correct` and `JM-Correct` include `X1`, `X2`, `X3`, `X4`, and centered
`X3^2` in both of their parametric submodels. Their misspecified partners are
identical except that both submodels retain only `X1` and `X2`.

The current comparator implementations use `B1` as their longitudinal/current
biomarker: MMLM and JM fit the `B1` trajectory, and RSF, DeepSurv, and
TimeError-Integrated-Landmark receive current `B1`. SILK alone currently uses
the joint `B1`--`B4` history. This distinction must remain explicit in reports;
calling these comparators four-biomarker methods would be inaccurate. A future
multivariate-comparator arm should be added as a separate prespecified design,
not silently substituted after inspecting outcomes.

### Required validation before the confirmatory rerun

Run the one-replication validation across the whole design first:

```bash
sbatch job_validate.sbatch
```

When its dependent collection job succeeds, run:

```bash
sbatch job_confirmatory.sbatch
```

The primary report writes the paired estimate and confidence interval for
`IBS(SILK-Gaussian) - IBS(Oracle-Cox)`. It also compares `SILK-Gaussian`
pairwise with every prespecified feasible comparator and writes
`q2_win_tie_summary.csv`.

## Repair only the failed comparators

The historical August run failed only for DeepSurv prediction and parts of JM.
DeepSurv now passes a data frame at the
`survivalmodels` prediction boundary. JM now uses a positive follow-up clock
relative to each subject's recorded landmark, filters invalid longitudinal
rows, and falls back from a random slope to a random intercept when necessary.

After copying or pulling these corrected simulation files to the cluster, run:

```bash
cd "/path/to/SILK/git/example/sim codes"
sbatch job_failed_comparators.sbatch
```

This submits one unthrottled CPU job array over the existing 4,800 design
tasks, with a two-hour limit per task, but evaluates only `JM-Correct`,
`JM-Misspecified`, and `DeepSurv`. It does not rerun SILK or any successful comparator, does
not overwrite the original results, and does not launch manuscript figures.
The corrected JM implementation is rerun for every design cell so one reported
method is not assembled from two different internal time definitions. The
dependent audit writes `failed_comparator_repair_summary.csv` and exits nonzero
if any method-task pair is missing or failed.

## Required outputs

- `raw/prediction_*`: predictions, metrics, paired contrasts, and method status.
- `raw/registration_diagnostics_*`: shift and stage recovery by kernel and task.
- `summary/prediction_method_failure_summary.csv`: comparator failure audit.
- `summary/registration_diagnostics_summary.csv`: Monte Carlo means and MCSEs
  for the registration chain.
- `summary/silk_package_provenance.txt`: repository, commit, R library, and
  Python interpreter used by the run.

Do not update the manuscript from a partial run. Require complete method status,
successful DeepSurv/JM smoke tests, stable registration diagnostics, and paired
Monte Carlo summaries first.
