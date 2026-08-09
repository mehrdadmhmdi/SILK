#!/usr/bin/env bash
# Remove the user-imposed concurrency throttle from already-submitted SILK
# arrays. With no arguments, discover this user's active jobs named silk-*.
# Alternatively, pass one or more array master job IDs explicitly.

set -euo pipefail

if (( $# > 0 )); then
  candidates=("$@")
else
  mapfile -t candidates < <(
    squeue --noheader --user "${USER:?USER is not set}" --format='%F|%j' |
      awk -F'|' '$2 ~ /^silk-/ && $1 ~ /^[0-9]+$/ { print $1 }' |
      sort -u
  )
fi

if (( ${#candidates[@]} == 0 )); then
  echo "No active SILK jobs were found. Pass array master job IDs explicitly:" >&2
  echo "  bash remove_array_throttle.sh 12345678 12345679" >&2
  exit 1
fi

updated=0
for raw_id in "${candidates[@]}"; do
  job_id="${raw_id%%_*}"
  job_id="${job_id%%.*}"
  if [[ ! "$job_id" =~ ^[0-9]+$ ]]; then
    echo "Skipping invalid job ID: $raw_id" >&2
    continue
  fi

  if ! before="$(scontrol show job --oneliner "$job_id" 2>/dev/null)"; then
    echo "Skipping job $job_id: it is no longer visible to scontrol." >&2
    continue
  fi
  if [[ "$before" != *"ArrayTaskId="* ]]; then
    echo "Skipping job $job_id: it is not an array master job."
    continue
  fi

  scontrol update JobId="$job_id" ArrayTaskThrottle=0
  after="$(scontrol show job --oneliner "$job_id")"
  array_spec="$(sed -n 's/.*ArrayTaskId=\([^ ]*\).*/\1/p' <<< "$after")"
  echo "Removed array throttle from job $job_id (ArrayTaskId=${array_spec:-unknown})."
  updated=$((updated + 1))
done

if (( updated == 0 )); then
  echo "No array throttles were changed." >&2
  exit 1
fi
