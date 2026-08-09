#!/usr/bin/env bash
#SBATCH --job-name=silk-submit
#SBATCH --account=bhmn-delta-cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
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

# Confirmatory defaults. When this file is submitted from a Git checkout, use
# that exact commit for the installed SILK package so the simulation scripts
# and package implementation cannot silently drift apart. Each setting remains
# overridable by exporting it before sbatch.
REPO_DIR="$(cd "$SIM_DIR/../.." && pwd)"
if git -C "$REPO_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  DEFAULT_PACKAGE_REF="$(git -C "$REPO_DIR" rev-parse HEAD)"
else
  DEFAULT_PACKAGE_REF="master"
fi

export SILK_PACKAGE_REF="${SILK_PACKAGE_REF:-$DEFAULT_PACKAGE_REF}"
export SILK_N_REP="${SILK_N_REP:-50}"
export SILK_BIOMARKER_KERNELS="${SILK_BIOMARKER_KERNELS:-gaussian,laplace,matern32}"
export SILK_ENABLE_JMBAYES2="${SILK_ENABLE_JMBAYES2:-true}"
export SILK_INCLUDE_HISTORY_ABLATION="${SILK_INCLUDE_HISTORY_ABLATION:-false}"
export SILK_OUT_DIR="${SILK_OUT_DIR:-$SIM_DIR/outputs-confirmatory}"

echo "SILK confirmatory configuration"
echo "  package ref:      $SILK_PACKAGE_REF"
echo "  output directory: $SILK_OUT_DIR"
echo "  replicates:       $SILK_N_REP"
echo "  kernels:          $SILK_BIOMARKER_KERNELS"
echo "  JMbayes2:         $SILK_ENABLE_JMBAYES2"
echo "  history ablation: $SILK_INCLUDE_HISTORY_ABLATION"

mkdir -p "$SILK_OUT_DIR/logs" "$SILK_OUT_DIR/summary"

# Prepare or reuse one persistent R/Python environment. Every array worker
# receives it through --export=ALL and never installs software itself.
source "$SIM_DIR/cluster_prepare_package.sh"

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
    --account=bhmn-delta-cpu \
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
figure_job="$(sbatch --parsable \
  --dependency="afterok:${check_job}" \
  --export=ALL \
  --output="$SILK_OUT_DIR/logs/%x-%j.out" \
  "$SIM_DIR/job_fig.sbatch")"

echo "Audit/collection job: $check_job"
echo "Figure job: $figure_job"
echo "Output directory: $SILK_OUT_DIR"
