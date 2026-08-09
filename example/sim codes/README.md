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
sbatch job_full.sbatch
```

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

Useful overrides include:

```bash
export SILK_OUT_DIR="$PWD/outputs-confirmatory"
export SILK_PACKAGE_REF="<pushed-commit-or-tag>"
export SILK_BIOMARKER_KERNELS="gaussian,laplace,matern32"
export SILK_N_REP=50
sbatch job_full.sbatch
```

The same-feature history ablation is intentionally excluded from the primary
roster. Run it into a separate supplementary output directory with:

```bash
export SILK_INCLUDE_HISTORY_ABLATION=true
export SILK_OUT_DIR="$PWD/outputs-history-ablation"
sbatch job_full.sbatch
```

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
