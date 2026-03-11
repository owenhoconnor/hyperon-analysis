#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

FCL="run_analyzeEvents.fcl"
SRC_DIR="/data/sbnd/background"
NBATCH=10

LOGDIR="logs_bkg"
WORKDIR="jobs_bkg"
PROJDIR="/home/lar/ooconnor/hyperons/srcs/sbndcode/sbndcode/Hyperons"
MERGED_OUT="mergedBkg.root"
TREE_OUT="analysisOutput.root"

mkdir -p "$LOGDIR" "$WORKDIR"

# Collect all input files
files=( "$SRC_DIR"/*.root )

if ((${#files[@]} == 0)); then
  echo "No .root files found in $SRC_DIR"
  exit 1
fi

echo "Found ${#files[@]} input files"

# Split files as evenly as possible across NBATCH batches
declare -a batch_lists
for ((i=0; i<NBATCH; ++i)); do
  batch_lists[i]=""
done

for ((i=0; i<${#files[@]}; ++i)); do
  b=$(( i % NBATCH ))
  batch_lists[$b]+="${files[$i]}"$'\n'
done

run_batch () {
  local batch_id="$1"
  local jobdir="$WORKDIR/batch${batch_id}"
  local logfile="$LOGDIR/batch${batch_id}.log"
  local listfile="$jobdir/input_files.txt"

  mkdir -p "$jobdir"

  printf "%s" "${batch_lists[$((batch_id-1))]}" > "$listfile"

  # Skip empty batches if there are fewer files than NBATCH
  if [[ ! -s "$listfile" ]]; then
    echo "[$(date)] batch${batch_id} is empty, skipping" | tee "$logfile"
    return 0
  fi

  echo "[$(date)] Starting batch${batch_id}"
  (
    cd "$jobdir"

    rm -f memory.db cputime.db messages.log errors.log "$TREE_OUT"

    args=()
    while IFS= read -r f; do
      [[ -z "$f" ]] && continue
      args+=(-s "$f")
    done < input_files.txt

    lar -c "$PROJDIR/$FCL" "${args[@]}" > "$PROJDIR/$logfile" 2>&1
  )
  echo "[$(date)] Finished batch${batch_id}" >> "$logfile"
}

# Launch batches in parallel
pids=()
for ((b=1; b<=NBATCH; ++b)); do
  run_batch "$b" &
  pids+=( "$!" )
done

status=0
for pid in "${pids[@]}"; do
  wait "$pid" || status=1
done

if [[ $status -ne 0 ]]; then
  echo "[$(date)] At least one batch failed. Check logs in $LOGDIR/"
  exit 1
fi

echo "[$(date)] All batches finished successfully"

# Gather only existing analyzer outputs
outputs=()
for ((b=1; b<=NBATCH; ++b)); do
  f="$WORKDIR/batch${b}/$TREE_OUT"
  [[ -f "$f" ]] && outputs+=( "$f" )
done

if ((${#outputs[@]} == 0)); then
  echo "No $TREE_OUT files found to hadd."
  exit 1
fi

hadd -f "$MERGED_OUT" "${outputs[@]}"

echo "[$(date)] Merged analysis output written to $MERGED_OUT"
