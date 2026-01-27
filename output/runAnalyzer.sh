#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

FCL="${FCL:-../run_analyzeEvents.fcl}"
BKG_DIR="${BKG_DIR:-/data/sbnd/background}"
BADLIST="${BADLIST:-./bad_files.txt}"

# If you ever hit "Argument list too long", lower this (e.g. 200)
CHUNK="${CHUNK:-500}"

# Load bad basenames into a set
declare -A BAD=()
while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "$line" ]] && continue
  [[ "$line" =~ ^# ]] && continue
  BAD["$line"]=1
done < "$BADLIST"

# Gather good files
good=()
for f in "$BKG_DIR"/reco2-*.root; do
  base="$(basename "$f")"
  [[ ${BAD["$base"]+x} ]] && continue
  good+=("$f")
done

# Nothing to run?
if (( ${#good[@]} == 0 )); then
  echo "No background files left after exclusions."
  exit 0
fi

# Run lar, chunking to avoid arg-length limits
i=0
while (( i < ${#good[@]} )); do
  args=()
  for (( j=0; j<CHUNK && i<${#good[@]}; j++, i++ )); do
    args+=(-s "${good[i]}")
  done
  lar -c "$FCL" "${args[@]}"
done

