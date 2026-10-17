#!/usr/bin/env bash
# Flat-job driver, design v2.0. Run from simulation/.
#   runs/launch.sh [-j JOBS] [-n N_REPS] [-b BATCH_SIZE] [-o OUT_DIR] [-t TIMEOUT_S] [-c CELL_REGEX]
# Builds the pending (cell, batch) list, skips batches whose .rds exists,
# runs JOBS single-threaded R processes at a time, each under `timeout`.
# Re-run the same command to resume. Works on macOS, Linux, SLURM login nodes.
# -c restricts to cells whose id matches the regex (e.g. -c N1000, -c x_classical).
# On macOS prefix with `caffeinate -i` so the machine does not sleep mid-run.
set -euo pipefail
JOBS=$(( $(sysctl -n hw.physicalcpu 2>/dev/null || nproc) - 1 ))
N_REPS=500; BATCH=25; OUT=sim_results_v2; TIMEOUT=$((6*3600)); CELLRE="."
while getopts "j:n:b:o:t:c:" opt; do case $opt in
  j) JOBS=$OPTARG;; n) N_REPS=$OPTARG;; b) BATCH=$OPTARG;; o) OUT=$OPTARG;; t) TIMEOUT=$OPTARG;; c) CELLRE=$OPTARG;;
esac; done
mkdir -p "$OUT"
# VECLIB_MAXIMUM_THREADS matters on macOS, where R links Apple's vecLib: the
# variable is read when the BLAS initialises, so pin_threads() setting it
# inside an already-started R process is too late. Exported here instead.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript -e "source('runs/grid.R'); j <- build_jobs(n_reps=$N_REPS, batch_size=$BATCH);
  f <- file.path('$OUT', j\$cell_id, sprintf('batch_%03d.rds', j\$batch));
  j <- j[!file.exists(f) & grepl('$CELLRE', j\$cell_id), ]; writeLines(paste(j\$cell_id, j\$batch), 'runs/_pending.txt');
  cat(nrow(j), 'pending jobs\n')" 2>/dev/null
[ -s runs/_pending.txt ] || { echo "nothing to do"; exit 0; }
TO=$(command -v timeout || command -v gtimeout || true)
run_one() {
  cell=$1; b=$2
  log="$OUT/_logs/${cell}__b${b}.log"; mkdir -p "$OUT/_logs"
  if [ -n "$TO" ]; then
    $TO "$TIMEOUT" Rscript runs/run_batch.R --cell="$cell" --batch="$b" --n_reps="$N_REPS" --batch_size="$BATCH" --out="$OUT" >"$log" 2>&1 || echo "FAILED $cell $b (see $log)"
  else
    Rscript runs/run_batch.R --cell="$cell" --batch="$b" --n_reps="$N_REPS" --batch_size="$BATCH" --out="$OUT" >"$log" 2>&1 || echo "FAILED $cell $b (see $log)"
  fi
  tail -1 "$log"
}
export -f run_one; export OUT N_REPS BATCH TIMEOUT TO
xargs -P "$JOBS" -L 1 bash -c 'run_one "$0" "$1"' < runs/_pending.txt | tee -a "$OUT/_progress.log"
