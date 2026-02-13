#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

SRC_DIR="/data/sbnd/hyperons_new"
FCL="run_analyzeEvents.fcl"
BAD_LIST="bad_files.txt"

# Read bad files into an associative array for fast lookup
declare -A bad=()
while IFS= read -r line; do
  # skip empty lines and comments
  [[ -z "${line// }" ]] && continue
  [[ "${line}" =~ ^[[:space:]]*# ]] && continue
  bad["$line"]=1
done < "$BAD_LIST"

# Collect all .root files and filter out the bad ones
args=()
for f in "$SRC_DIR"/*.root; do
  [[ -n "${bad[$f]+x}" ]] && continue
  args+=(-s "$f")
done

# Safety check
if ((${#args[@]} == 0)); then
  echo "No input files left after excluding bad files."
  exit 1
fi

echo "Running lar on $(( ${#args[@]} / 2 )) files (excluded ${#bad[@]} listed bad files)."
lar -c "$FCL" "${args[@]}"

