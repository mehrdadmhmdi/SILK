#!/usr/bin/env bash
#SBATCH --job-name=silk-confirm
#SBATCH --account=bfez-delta-cpu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.log

set -euo pipefail

module purge
module load slurm-env/0.1
module load cray-R/4.4.0

SIM_DIR="${SILK_SIM_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
SIM_DIR="$(cd "$SIM_DIR" && pwd)"
export SILK_SIM_DIR="$SIM_DIR"

# If the simulation tasks succeeded but validation reporting stopped before its
# completion marker, repair the reports from the existing raw files first. This
# does not rerun any of the 96 validation array tasks.
LOCKED_CONFIG_FILE="$SIM_DIR/locked_confirmatory_config.sh"
if [[ ! -f "$LOCKED_CONFIG_FILE" ]]; then
  echo "Missing locked configuration: $LOCKED_CONFIG_FILE" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "$LOCKED_CONFIG_FILE"

if command -v sha256sum >/dev/null 2>&1; then
  current_config_hash="$(sha256sum "$LOCKED_CONFIG_FILE" | awk '{print $1}')"
else
  current_config_hash="$(shasum -a 256 "$LOCKED_CONFIG_FILE" | awk '{print $1}')"
fi

VALIDATION_DIR="${SILK_VALIDATION_OUT_DIR:-$SIM_DIR/outputs-validation-nrep1}"
validation_config="$VALIDATION_DIR/summary/locked_configuration.txt"
validation_provenance="$VALIDATION_DIR/summary/silk_package_provenance.txt"
validation_marker="$VALIDATION_DIR/summary/collection_report_complete.ok"
if [[ ! -s "$validation_config" || ! -s "$validation_provenance" ]]; then
  echo "Validation metadata is missing under $VALIDATION_DIR." >&2
  echo "Run sbatch job_validate.sbatch first." >&2
  exit 1
fi

validation_hash="$(awk -F= '$1 == "config_sha256" {print $2}' "$validation_config")"
if [[ "$validation_hash" != "$current_config_hash" ]]; then
  echo "The locked configuration changed after validation." >&2
  echo "Run sbatch job_validate.sbatch again before confirmation." >&2
  exit 1
fi

validation_commit="$(awk -F= '$1 == "commit" {print $2}' "$validation_provenance")"
if [[ -z "$validation_commit" ]]; then
  echo "The validated SILK package commit could not be read." >&2
  exit 1
fi
export SILK_PACKAGE_REF="$validation_commit"

if [[ ! -s "$validation_marker" ]]; then
  echo "Validation completion marker is absent."
  echo "Checking whether the existing raw tasks have the current method roster before rebuilding reports."
  export SILK_RUN_STAGE="validation"
  export SILK_N_REP=1
  export REPORT_N_REP=1
  export SILK_OUT_DIR="$VALIDATION_DIR"
  # Install/reuse the exact package commit that generated the validation data.
  # shellcheck source=/dev/null
  source "$SIM_DIR/cluster_prepare_package.sh"
  bash "$SIM_DIR/job_check.sbatch"
fi

# job_full.sbatch performs the final marker/configuration/package checks before
# submitting the 50-replication arrays.
export SILK_RUN_STAGE="confirmatory"
unset SILK_OUT_DIR
exec bash "$SIM_DIR/job_full.sbatch"
