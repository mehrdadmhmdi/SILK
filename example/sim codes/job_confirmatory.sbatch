#!/usr/bin/env bash
#SBATCH --job-name=silk-confirm
#SBATCH --account=bhmn-delta-cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

module purge
module load slurm-env/0.1
module load cray-R/4.4.0

SIM_DIR="${SILK_SIM_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIM_DIR="$(cd "$SIM_DIR" && pwd)"
export SILK_SIM_DIR="$SIM_DIR"
export SILK_RUN_STAGE="confirmatory"

# job_full.sbatch refuses this launch unless the one-replication run used the
# same locked configuration and package commit and completed every method.
exec bash "$SIM_DIR/job_full.sbatch"
