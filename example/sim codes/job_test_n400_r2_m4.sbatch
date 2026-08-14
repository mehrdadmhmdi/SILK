#!/usr/bin/env bash
#SBATCH --job-name=silk-test-400-r2-m4
#SBATCH --account=bfez-delta-cpu
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

# Exact smoke-test design:
#   - all 12 error scenarios and all 12 methods
#   - n_train = 400
#   - R = 2 Monte Carlo replications
#   - primary design only, forcing the m4 (four-visit) schedule
#   - standard prediction horizons 1, 2, 3, and 4
# Results are isolated from validation and confirmatory output directories.
export SILK_RUN_STAGE="smoke_n400_r2_m4"

exec bash "$SIM_DIR/job_full.sbatch"
