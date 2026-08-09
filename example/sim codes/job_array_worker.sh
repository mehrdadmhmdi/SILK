#!/usr/bin/env bash
#SBATCH --account=bhmn-delta-cpu
#SBATCH --partition=cpu

set -euo pipefail

module purge
module load slurm-env/0.1
module load cray-R/4.4.0

SIM_DIR="${SILK_SIM_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIM_DIR="$(cd "$SIM_DIR" && pwd)"
if [[ ! -f "$SIM_DIR/task.R" ]]; then
  echo "task.R was not found under $SIM_DIR; SILK_SIM_DIR is invalid." >&2
  exit 1
fi
cd "$SIM_DIR"

: "${SILK_TASK_LIST_FILE:?SILK_TASK_LIST_FILE was not exported by job_full.sbatch}"
: "${SLURM_ARRAY_TASK_ID:?This worker must run as a Slurm array task}"
: "${SILK_R_LIB:?SILK_R_LIB was not prepared by cluster_prepare_package.sh}"

mapfile -t task_ids < "$SILK_TASK_LIST_FILE"
array_index=$((SLURM_ARRAY_TASK_ID - 1))
if (( array_index < 0 || array_index >= ${#task_ids[@]} )); then
  echo "Array index $SLURM_ARRAY_TASK_ID is outside $SILK_TASK_LIST_FILE" >&2
  exit 1
fi
export SILK_TASK_ID="${task_ids[$array_index]}"

allocated_cores="${SLURM_CPUS_PER_TASK:-1}"
default_parallel_cores="$allocated_cores"
if (( default_parallel_cores > 4 )); then default_parallel_cores=4; fi

# Phase 1: independent cross-fit folds. Phase 2: independent survival methods.
# The phases are sequential, so each may safely use the same bounded core cap.
export SILK_N_CORES="${SILK_REGISTRATION_CORES:-$default_parallel_cores}"
export SILK_METHOD_CORES="${SILK_METHOD_CORES:-$default_parallel_cores}"

# Keep BLAS/OpenMP libraries single-threaded inside each R worker; otherwise a
# four-process phase can multiply into dozens of unrequested threads.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export R_LIBS="$SILK_R_LIB${R_LIBS:+:$R_LIBS}"
Rscript -e 'library(SILK); stopifnot(is.function(fit_silk_registration), identical(silk_default_options()$BIOMARKER_BANDWIDTH, "median"))'

echo "Task $SILK_TASK_ID; class ${SILK_RESOURCE_CLASS:-unclassified}; registration cores $SILK_N_CORES; method cores $SILK_METHOD_CORES"
exec Rscript "$SIM_DIR/task.R"
