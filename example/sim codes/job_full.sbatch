#!/usr/bin/env bash
#SBATCH --job-name=silk-submit
#SBATCH --account=bfez-delta-cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

module purge
module load slurm-env/0.1
module load cray-R/4.4.0

SIM_DIR="${SILK_SIM_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIM_DIR="$(cd "$SIM_DIR" && pwd)"
if [[ ! -f "$SIM_DIR/task.R" || ! -f "$SIM_DIR/cluster_prepare_package.sh" ]]; then
  echo "Simulation files were not found under $SIM_DIR" >&2
  echo "Submit from the simulation directory: cd '<repo>/example/sim codes' && sbatch job_full.sbatch" >&2
  exit 1
fi
export SILK_SIM_DIR="$SIM_DIR"
cd "$SIM_DIR"

LOCKED_CONFIG_FILE="$SIM_DIR/locked_confirmatory_config.sh"
if [[ ! -f "$LOCKED_CONFIG_FILE" ]]; then
  echo "Missing locked configuration: $LOCKED_CONFIG_FILE" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "$LOCKED_CONFIG_FILE"

# The validation and confirmatory jobs differ only in replication count and
# output directory. The smoke stage is an explicitly isolated test design:
# n=400, R=2, primary phase only, with the four-visit schedule for every
# scenario. All method hyperparameters still come from the locked file above.
SILK_RUN_STAGE="${SILK_RUN_STAGE:-validation}"
case "$SILK_RUN_STAGE" in
  validation)
    export SILK_N_REP=1
    export SILK_DESIGN_PHASES="primary,timepoint_extension"
    unset SILK_FORCE_SCHEDULE
    export SILK_OUT_DIR="${SILK_VALIDATION_OUT_DIR:-$SIM_DIR/outputs-validation-nrep1}"
    ;;
  confirmatory)
    export SILK_N_REP=50
    export SILK_DESIGN_PHASES="primary,timepoint_extension"
    unset SILK_FORCE_SCHEDULE
    export SILK_OUT_DIR="${SILK_CONFIRMATORY_OUT_DIR:-$SIM_DIR/outputs-confirmatory-50rep}"
    ;;
  smoke_n400_r2_m4)
    export SILK_N_REP=2
    export SILK_N_TRAIN_GRID="400"
    export SILK_TP_N_TRAIN="400"
    export SILK_DESIGN_PHASES="primary"
    export SILK_FORCE_SCHEDULE="m4"
    export SILK_OUT_DIR="${SILK_TEST_OUT_DIR:-$SIM_DIR/outputs-test-n400-r2-m4}"
    ;;
  *)
    echo "SILK_RUN_STAGE must be validation, confirmatory, or smoke_n400_r2_m4." >&2
    exit 1
    ;;
esac
export SILK_RUN_STAGE
export REPORT_N_REP="$SILK_N_REP"

if command -v sha256sum >/dev/null 2>&1; then
  SILK_CONFIG_SHA256="$(sha256sum "$LOCKED_CONFIG_FILE" | awk '{print $1}')"
else
  SILK_CONFIG_SHA256="$(shasum -a 256 "$LOCKED_CONFIG_FILE" | awk '{print $1}')"
fi
export SILK_CONFIG_SHA256

# When this file is submitted from a Git checkout, use
# that exact commit for the installed SILK package so the simulation scripts
# and package implementation cannot silently drift apart.
REPO_DIR="$(cd "$SIM_DIR/../.." && pwd)"
if git -C "$REPO_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  if [[ -n "$(git -C "$REPO_DIR" status --porcelain --untracked-files=normal)" ]]; then
    echo "The SILK checkout has uncommitted or untracked changes." >&2
    echo "Commit and push the locked simulation/package version before submission." >&2
    exit 1
  fi
  DEFAULT_PACKAGE_REF="$(git -C "$REPO_DIR" rev-parse HEAD)"
else
  DEFAULT_PACKAGE_REF="master"
fi

export SILK_PACKAGE_REF="${SILK_PACKAGE_REF:-$DEFAULT_PACKAGE_REF}"

if [ -n "${SILK_INCLUDE_HISTORY_ABLATION:-}" ]; then
  echo "WARNING: SILK_INCLUDE_HISTORY_ABLATION is obsolete and ignored." >&2
  echo "         Use SILK_PAPER_FEATURE_MAP instead." >&2
fi

echo "SILK $SILK_RUN_STAGE configuration"
echo "  config id:           $SILK_CONFIG_ID"
echo "  config sha256:       $SILK_CONFIG_SHA256"
echo "  tuning source:       $SILK_HYPERPARAMETER_SELECTION"
echo "  package ref:         $SILK_PACKAGE_REF"
echo "  output directory:    $SILK_OUT_DIR"
echo "  replicates:          $SILK_N_REP"
echo "  design phases:       ${SILK_DESIGN_PHASES:-primary,timepoint_extension}"
echo "  forced schedule:     ${SILK_FORCE_SCHEDULE:-none}"
echo "  kernels:             $SILK_BIOMARKER_KERNELS"
echo "  JMbayes2:            $SILK_ENABLE_JMBAYES2"
echo "  paper feature map:   $SILK_PAPER_FEATURE_MAP  (primary roster)"
echo "  profile search rad.: $SILK_PROFILE_SEARCH_RADIUS"
echo "  alt local radius:    $SILK_ALT_LOCAL_RADIUS"
echo "  registration starts: $SILK_N_STARTS"
echo "  bandwidth scale:     $SILK_BIOMARKER_BW_SCALE"
echo "  feature map:         $SILK_FEATURE_MAP"

mkdir -p "$SILK_OUT_DIR/logs" "$SILK_OUT_DIR/summary"
if [[ "$SILK_RUN_STAGE" != "confirmatory" ]]; then
  rm -f "$SILK_OUT_DIR/summary/collection_report_complete.ok"
fi

cat > "$SILK_OUT_DIR/summary/locked_configuration.txt" <<EOF
config_id=$SILK_CONFIG_ID
config_sha256=$SILK_CONFIG_SHA256
hyperparameter_selection=$SILK_HYPERPARAMETER_SELECTION
run_stage=$SILK_RUN_STAGE
n_rep=$SILK_N_REP
n_train_grid=$SILK_N_TRAIN_GRID
n_test=$SILK_N_TEST
design_phases=${SILK_DESIGN_PHASES:-primary,timepoint_extension}
forced_schedule=${SILK_FORCE_SCHEDULE:-none}
kernels=$SILK_BIOMARKER_KERNELS
biomarker_bandwidth=$SILK_BIOMARKER_BW
biomarker_bandwidth_scale=$SILK_BIOMARKER_BW_SCALE
shift_grid_step=$SILK_SHIFT_GRID_STEP
profile_search_radius=$SILK_PROFILE_SEARCH_RADIUS
alt_local_radius=$SILK_ALT_LOCAL_RADIUS
n_starts=$SILK_N_STARTS
n_folds=$SILK_N_FOLDS
feature_map=$SILK_FEATURE_MAP
paper_feature_map=$SILK_PAPER_FEATURE_MAP
package_ref=$SILK_PACKAGE_REF
EOF

if [[ "$SILK_RUN_STAGE" == "confirmatory" ]]; then
  SILK_VALIDATION_DIR="${SILK_VALIDATION_DIR:-${SILK_VALIDATION_OUT_DIR:-$SIM_DIR/outputs-validation-nrep1}}"
  export SILK_VALIDATION_DIR
  validation_config="$SILK_VALIDATION_DIR/summary/locked_configuration.txt"
  validation_audit="$SILK_VALIDATION_DIR/summary/prediction_replication_audit.csv"
  validation_failures="$SILK_VALIDATION_DIR/summary/prediction_method_failure_summary.csv"
  validation_complete="$SILK_VALIDATION_DIR/summary/collection_report_complete.ok"
  validation_q1="$SILK_VALIDATION_DIR/summary/two_questions/q1_primary.csv"
  validation_q2="$SILK_VALIDATION_DIR/summary/two_questions/q2_win_tie_summary.csv"
  if [[ ! -s "$validation_config" || ! -s "$validation_audit" ||
        ! -s "$validation_failures" || ! -s "$validation_complete" ||
        ! -s "$validation_q1" || ! -s "$validation_q2" ]]; then
    echo "The one-replication validation is incomplete under $SILK_VALIDATION_DIR." >&2
    echo "Run sbatch job_validate.sbatch and wait for its collection job to succeed." >&2
    exit 1
  fi
  validation_hash="$(awk -F= '$1 == "config_sha256" {print $2}' "$validation_config")"
  if [[ "$validation_hash" != "$SILK_CONFIG_SHA256" ]]; then
    echo "Locked configuration changed after validation; rerun job_validate.sbatch." >&2
    exit 1
  fi
  Rscript --vanilla -e '
    root <- Sys.getenv("SILK_VALIDATION_DIR")
    audit <- read.csv(file.path(root, "summary", "prediction_replication_audit.csv"))
    failures <- read.csv(file.path(root, "summary", "prediction_method_failure_summary.csv"))
    if (!isTRUE(as.logical(audit$complete[1L]))) stop("Validation task audit is incomplete.")
    if (any(failures$failed_tasks > 0L, na.rm = TRUE)) stop("At least one method failed validation.")
  '
fi

# Prepare or reuse one persistent R/Python environment. Every array worker
# receives it through --export=ALL and never installs software itself.
source "$SIM_DIR/cluster_prepare_package.sh"

if [[ "$SILK_RUN_STAGE" == "confirmatory" ]]; then
  validation_commit="$(awk -F= '$1 == "commit" {print $2}' \
    "$SILK_VALIDATION_DIR/summary/silk_package_provenance.txt")"
  if [[ -z "$validation_commit" || "$validation_commit" != "$SILK_PACKAGE_COMMIT" ]]; then
    echo "Package commit differs from the validated package; rerun validation." >&2
    exit 1
  fi
fi

Rscript -e 'source(file.path("R", "run_simulation_comparisons.R")); source_silk_prediction_modules(); source(file.path("R", "resource_tasks.R")); write_resource_task_lists()'

GROUP_FILE="$SILK_OUT_DIR/summary/prediction_resource_groups.tsv"
if [[ ! -s "$GROUP_FILE" ]]; then
  echo "Resource-group file was not created: $GROUP_FILE" >&2
  exit 1
fi

array_jobs=()
while IFS=$'\t' read -r resource_class n_tasks _sbatch task_list cpus memory walltime _description; do
  [[ "$resource_class" == "resource_class" ]] && continue
  [[ -z "$resource_class" || "$n_tasks" -lt 1 ]] && continue
  # Submit the full array without a user-imposed concurrency throttle. Slurm
  # will still enforce partition, account, QOS, and site-wide scheduling limits.
  array_spec="1-${n_tasks}"
  job_id="$(sbatch --parsable \
    --job-name="silk-${resource_class}" \
    --account=bfez-delta-cpu \
    --partition=cpu \
    --cpus-per-task="$cpus" \
    --mem="$memory" \
    --time="$walltime" \
    --array="$array_spec" \
    --output="$SILK_OUT_DIR/logs/%x-%A_%a.out" \
    --export="ALL,SILK_RESOURCE_CLASS=${resource_class},SILK_TASK_LIST_FILE=${task_list}" \
    "$SIM_DIR/job_array_worker.sh")"
  array_jobs+=("$job_id")
  echo "Submitted $resource_class: $job_id ($n_tasks tasks; $cpus CPUs; $memory; $walltime)"
done < "$GROUP_FILE"

if [[ "${#array_jobs[@]}" -eq 0 ]]; then
  echo "No simulation arrays were submitted." >&2
  exit 1
fi

dependency="$(IFS=:; echo "${array_jobs[*]}")"
check_job="$(sbatch --parsable \
  --dependency="afterany:${dependency}" \
  --export=ALL \
  --output="$SILK_OUT_DIR/logs/%x-%j.out" \
  "$SIM_DIR/job_check.sbatch")"
echo "Audit/collection job: $check_job"
if [[ "$SILK_RUN_STAGE" == "confirmatory" ]]; then
  figure_job="$(sbatch --parsable \
    --dependency="afterok:${check_job}" \
    --export=ALL \
    --output="$SILK_OUT_DIR/logs/%x-%j.out" \
    "$SIM_DIR/job_fig.sbatch")"
  echo "Figure job: $figure_job"
else
  echo "$SILK_RUN_STAGE stage: figures are not submitted."
fi
echo "Output directory: $SILK_OUT_DIR"
