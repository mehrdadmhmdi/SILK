#!/usr/bin/env bash
#SBATCH --job-name=silk-collect
#SBATCH --account=bhmn-delta-cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00

set -euo pipefail
module purge
module load slurm-env/0.1
module load cray-R/4.4.0

SIM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SIM_DIR"
: "${SILK_R_LIB:?SILK_R_LIB is not set}"
export R_LIBS_USER="$SILK_R_LIB${R_LIBS_USER:+:$R_LIBS_USER}"
export SILK_N_CORES=1
exec Rscript "$SIM_DIR/collect.R"
