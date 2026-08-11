#!/usr/bin/env bash
#SBATCH --job-name=silk-repair-submit
#SBATCH --account=bhmn-delta-cpu
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
  echo "Submit from the simulation directory with: sbatch job_failed_comparators.sbatch" >&2
  exit 1
fi
export SILK_SIM_DIR="$SIM_DIR"
cd "$SIM_DIR"

REPO_DIR="$(cd "$SIM_DIR/../.." && pwd)"
if git -C "$REPO_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  DEFAULT_PACKAGE_REF="$(git -C "$REPO_DIR" rev-parse HEAD)"
else
  DEFAULT_PACKAGE_REF="master"
fi

export SILK_PACKAGE_REF="${SILK_PACKAGE_REF:-$DEFAULT_PACKAGE_REF}"
export SILK_N_REP="${SILK_N_REP:-50}"
export SILK_ENABLE_JMBAYES2=true
export SILK_ONLY_METHODS="JM-Recorded,DeepSurv-Observed"
# This job repairs two comparators in the LEGACY 8.9.2026 freeze, so it keeps
# that roster. SILK_INCLUDE_HISTORY_ABLATION is obsolete and no longer read;
# SILK_PAPER_FEATURE_MAP=false is the legacy-roster switch. Neither affects the
# two methods repaired here, but the setting is made explicit so the printed
# configuration matches what cfg.R actually does.
export SILK_PAPER_FEATURE_MAP=false
export SILK_OUT_DIR="${SILK_REPAIR_OUT_DIR:-$SIM_DIR/outputs-failed-comparators}"

echo "SILK failed-comparator repair"
echo "  methods:          $SILK_ONLY_METHODS"
echo "  package ref:      $SILK_PACKAGE_REF"
echo "  output directory: $SILK_OUT_DIR"
echo "  replicates:       $SILK_N_REP"

mkdir -p "$SILK_OUT_DIR/logs" "$SILK_OUT_DIR/summary"

# Prepare/reuse one persistent R and CPU-only Python environment. Array workers
# only load this environment; they never install software themselves.
source "$SIM_DIR/cluster_prepare_package.sh"

TASK_LIST="$SILK_OUT_DIR/summary/failed_comparator_task_ids.txt"
export SILK_REPAIR_TASK_LIST="$TASK_LIST"
Rscript --vanilla -e '
  source(file.path("R", "run_simulation_comparisons.R"))
  source_silk_prediction_modules()
  plan <- build_task_plan()
  writeLines(as.character(plan$task_id), Sys.getenv("SILK_REPAIR_TASK_LIST"), useBytes = TRUE)
  write.csv(plan, file.path(SUMMARY_DIR, "failed_comparator_task_plan.csv"), row.names = FALSE)
'

n_tasks="$(wc -l < "$TASK_LIST")"
if [[ -z "$n_tasks" || "$n_tasks" -lt 1 ]]; then
  echo "No repair tasks were written to $TASK_LIST" >&2
  exit 1
fi

# Deliberately no %N array throttle: only the site/account/QOS limits apply.
array_job="$(sbatch --parsable \
  --job-name=silk-failed-comparators \
  --account=bhmn-delta-cpu \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=16G \
  --time=02:00:00 \
  --array="1-${n_tasks}" \
  --output="$SILK_OUT_DIR/logs/%x-%A_%a.out" \
  --export="ALL,SILK_RESOURCE_CLASS=failed-comparators,SILK_TASK_LIST_FILE=${TASK_LIST}" \
  "$SIM_DIR/job_array_worker.sh")"

check_job="$(sbatch --parsable \
  --dependency="afterany:${array_job}" \
  --export=ALL \
  --output="$SILK_OUT_DIR/logs/%x-%j.out" \
  "$SIM_DIR/job_failed_comparators_check.sbatch")"

echo "Submitted repair array: $array_job ($n_tasks tasks, no user throttle)"
echo "Submitted repair audit: $check_job"
echo "Original outputs were not modified. Repair output: $SILK_OUT_DIR"
