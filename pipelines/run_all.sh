#!/usr/bin/env bash
# Driver: a quick n=20 check of the current pipeline, then full replicates sequentially.
# Usage: pipelines/run_all.sh <runs_dir> [first_seed] [last_seed] [n_ind] [jobs]
set -euo pipefail
R=$(realpath -m "$1"); A=${2:-1}; B=${3:-5}; N=${4:-400}; J=${5:-56}
REPO=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$R"
P() { echo "[$(date '+%F %T')] $*" | tee -a "$R/progress.log"; }
if [ ! -f "$R/check/.eval.done" ]; then
  P "check run (n=20, seed 102)"
  bash "$REPO/pipelines/run_rep.sh" 102 "$R/check" 20 20 > "$R/check.log" 2>&1 || { P "check run FAILED"; exit 1; }
  P "check run ok"
fi
for s in $(seq "$A" "$B"); do
  P "seed $s start"
  if bash "$REPO/pipelines/run_rep.sh" "$s" "$R/seed$s" "$N" "$J" > "$R/seed$s.log" 2>&1; then
    P "seed $s done ($(awk -F'\t' '{t+=$2} END{printf "%.1f h", t/3600}' "$R/seed$s/logs/time.tsv"))"
  else
    P "seed $s FAILED (see seed$s.log)"; exit 1
  fi
done
P "all done"
