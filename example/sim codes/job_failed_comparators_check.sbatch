#!/usr/bin/env bash
#SBATCH --job-name=silk-repair-check
#SBATCH --account=bfez-delta-cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail

module purge
module load slurm-env/0.1
module load cray-R/4.4.0

SIM_DIR="${SILK_SIM_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIM_DIR="$(cd "$SIM_DIR" && pwd)"
if [[ ! -f "$SIM_DIR/check_failed_comparators.R" ]]; then
  echo "check_failed_comparators.R was not found under $SIM_DIR" >&2
  exit 1
fi
cd "$SIM_DIR"

: "${SILK_R_LIB:?SILK_R_LIB is not set}"
export R_LIBS="$SILK_R_LIB${R_LIBS:+:$R_LIBS}"
export SILK_N_CORES=1
export SILK_METHOD_CORES=1
exec Rscript --vanilla "$SIM_DIR/check_failed_comparators.R"
