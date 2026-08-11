# SILK simulation pipeline

This directory is the cluster-ready implementation of the confirmatory SILK
simulation. The manuscript is not read or modified by these scripts.

## Scientific target

- Event times follow a Gompertz baseline hazard on attained latent age.
- The primary clock-isolation comparison holds the Cox learner and baseline
  covariates fixed and changes only the time coordinate: recorded age,
  SILK-calibrated age, or latent oracle age.
- Biomarkers enter the primary SILK-Cox analysis only through the registration
  layer. Recorded-clock Cox therefore cannot bypass registration using a
  hand-built biomarker summary.
- `weak_stage` is a prespecified stage-null negative control: biomarkers contain
  neither latent-age nor frailty information. Its expected result is no
  systematic registration gain; a guaranteed loss is not imposed.
- Gaussian RBF is the primary exact characteristic kernel. Exact Laplace and
  Matern-3/2 kernels are paired robustness analyses using the same data, folds,
  and Monte Carlo seeds. Linear and polynomial biomarker kernels are rejected.
- RSF and DeepSurv receive recorded age and the full observed-history summary.
  DeepSurv uses a fixed CPU network, a training-only validation split, early
  stopping, and no Python data-loader workers.

The simulation writes registration diagnostics for both cross-fitted training
subjects and held-out test subjects. Inspect shift/stage recovery before
interpreting prediction gains.

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
All hyperparameters are fixed in that file from a separate development design
not included in the confirmatory claims; neither launcher tunes against the
confirmatory outcomes.

The first controller run creates a persistent CPU environment under
`.cluster-env-r44/`, installs the R dependencies and CPU-only `pycox` stack,
and installs the requested SILK GitHub commit. Later runs reuse that environment
and reinstall SILK only when the requested commit changes. Array workers only
load the prepared environment; they never clone or install packages.

Every Slurm job uses:

```bash
#SBATCH --account=bhmn-delta-cpu
#SBATCH --partition=cpu

module purge
module load slurm-env/0.1
module load cray-R/4.4.0
```

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

### Primary roster and the paper feature map

The manuscript's survival layer is a Cox model for **residual time** with a
fixed feature map `g(age coordinate, invariant history)`, so the age coordinate
is a covariate inside `g`, not the time axis. The primary roster is therefore
the residual-time triple, which differs only in that coordinate:

| Method | Age coordinate |
|---|---|
| `Cox-History-Recorded` | recorded `A_obs` |
| `SILK-History`, `SILK-History-Laplace`, `SILK-History-Matern32` | cross-fitted calibrated age, per kernel |
| `Oracle-History-Latent-Age` | latent `A_star` (simulation-only) |

All three kernels share this survival layer, feature map, folds and seeds, so
the kernel-sensitivity contrast varies the kernel alone.

The attained-age methods (`Cox-SameFeature-Recorded`, `SILK`, `SILK-Laplace`,
`SILK-Matern32`, `Oracle-Latent-Age`) put age on the time axis with covariates
`(X1, X2)`. They are retained as a **supplementary clock-isolation analysis**
and must be labelled as such; they are not the paper's estimator.

`SILK_INCLUDE_HISTORY_ABLATION` is obsolete and is no longer read. The locked
rerun always uses the paper feature map. Reproducing the 8.9.2026 freeze should
be done from its archived commit, not by modifying the locked rerun scripts.

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
`IBS(SILK-History) - IBS(Oracle-History-Latent-Age)`. It also compares SILK
pairwise with every primary feasible comparator, including
`Cox-History-Recorded`, and writes `q2_win_tie_summary.csv`.

## Repair only the failed comparators

The August confirmatory run failed only for `DeepSurv-Observed` prediction and
parts of `JM-Recorded`. DeepSurv now passes a data frame at the
`survivalmodels` prediction boundary. JM now uses a positive follow-up clock
relative to each subject's recorded landmark, filters invalid longitudinal
rows, and falls back from a random slope to a random intercept when necessary.

After copying or pulling these corrected simulation files to the cluster, run:

```bash
cd "/path/to/SILK/git/example/sim codes"
sbatch job_failed_comparators.sbatch
```

This submits one unthrottled CPU job array over the existing 4,800 design
tasks, with a two-hour limit per task, but evaluates only `JM-Recorded` and
`DeepSurv-Observed`. It does not rerun SILK or any successful comparator, does
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
