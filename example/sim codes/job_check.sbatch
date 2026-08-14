#!/usr/bin/env bash
#SBATCH --job-name=silk-collect
#SBATCH --account=bfez-delta-cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.log

set -euo pipefail
module purge
module load slurm-env/0.1
module load cray-R/4.4.0

SIM_DIR="${SILK_SIM_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIM_DIR="$(cd "$SIM_DIR" && pwd)"
if [[ ! -f "$SIM_DIR/collect.R" || ! -f "$SIM_DIR/report_two_questions.R" ]]; then
  echo "Collection/reporting scripts were not found under $SIM_DIR." >&2
  exit 1
fi
cd "$SIM_DIR"
: "${SILK_R_LIB:?SILK_R_LIB is not set}"
export R_LIBS="$SILK_R_LIB${R_LIBS:+:$R_LIBS}"
export SILK_N_CORES=1
export REPORT_N_REP="${REPORT_N_REP:-${SILK_N_REP:-50}}"

Rscript "$SIM_DIR/collect.R"
Rscript "$SIM_DIR/report_two_questions.R" \
  "$SILK_OUT_DIR" \
  "$SILK_OUT_DIR/summary/two_questions" \
  --rebuild-cache

printf 'complete\n' > "$SILK_OUT_DIR/summary/collection_report_complete.ok"
