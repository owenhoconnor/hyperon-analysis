#!/usr/bin/env bash
set -euo pipefail

FCL="run_analyzeEvents.fcl"
LOGDIR="logs"
OUTDIR="output"
WORKDIR="jobs"
PROJDIR="/home/lar/ooconnor/hyperons/srcs/sbndcode/sbndcode/Hyperons"
mkdir -p "$LOGDIR" "$OUTDIR" "$WORKDIR"

run_batch () {
  local name="$1"
  local indir="$2"
  local jobdir="$WORKDIR/$name"
  local outfile="$OUTDIR/${name}.root"
  local logfile="$LOGDIR/${name}.log"

  mkdir -p "$jobdir"

  echo "[$(date)] Starting $name from $indir"
  (
    cd "$jobdir"

    # clean up any old tracker/log sqlite files from a previous failed run
    rm -f memory.db cputime.db messages.log errors.log

    lar -c "$PROJDIR/$FCL" -s "$indir"/*.root > "$PROJDIR/$logfile" 2>&1
  )

  echo "[$(date)] Finished $name"
}

run_batch firstbatch /data/sbnd/hyperons_new/firstbatch &
pid1=$!

run_batch secondbatch /data/sbnd/hyperons_new/secondbatch &
pid2=$!

run_batch thirdbatch /data/sbnd/hyperons_new/thirdbatch &
pid3=$!

status=0

wait $pid1 || status=1
wait $pid2 || status=1
wait $pid3 || status=1

if [[ $status -ne 0 ]]; then
  echo "[$(date)] At least one batch failed. Check logs in $LOGDIR/"
  exit 1
fi

echo "[$(date)] All jobs finished"

hadd -f "$OUTDIR/hyperons.root" \
  "$WORKDIR/firstbatch/analysisOutput.root" \
  "$WORKDIR/secondbatch/analysisOutput.root" \
  "$WORKDIR/thirdbatch/analysisOutput.root"

echo "[$(date)] Merged output written to $OUTDIR/hyperons.root"
