#!/usr/bin/env bash
#SBATCH --job-name=silk-figures
#SBATCH --account=bfez-delta-cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00

set -euo pipefail
module purge
module load slurm-env/0.1
module load cray-R/4.4.0

SIM_DIR="${SILK_SIM_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIM_DIR="$(cd "$SIM_DIR" && pwd)"
if [[ ! -f "$SIM_DIR/fig.R" ]]; then
  echo "fig.R was not found under $SIM_DIR; SILK_SIM_DIR is invalid." >&2
  exit 1
fi
cd "$SIM_DIR"
: "${SILK_R_LIB:?SILK_R_LIB is not set}"
export R_LIBS="$SILK_R_LIB${R_LIBS:+:$R_LIBS}"
export SILK_N_CORES=1
exec Rscript "$SIM_DIR/fig.R"
