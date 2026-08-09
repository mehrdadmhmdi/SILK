#!/usr/bin/env bash
#SBATCH --job-name=silk-submit
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

SIM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SIM_DIR"

export SILK_OUT_DIR="${SILK_OUT_DIR:-$SIM_DIR/outputs}"
export SILK_MAX_CONCURRENT="${SILK_MAX_CONCURRENT:-40}"
mkdir -p "$SILK_OUT_DIR/logs" "$SILK_OUT_DIR/summary"

# Install the exact source tree once. Every array worker receives this private
# library through --export=ALL and therefore cannot silently use a stale SILK.
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
  array_spec="1-${n_tasks}%${SILK_MAX_CONCURRENT}"
  job_id="$(sbatch --parsable \
    --job-name="silk-${resource_class}" \
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
