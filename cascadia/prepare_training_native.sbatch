#!/bin/bash
#SBATCH --job-name=casc_prep
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --qos=genome-center-grp-high-qos
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G                       # one .d's spectra at top 3000 peaks
#SBATCH --exclude=hive-as-11-2-41
#SBATCH --output=/quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native/train_v1/logs/prep_%A_%a.out
# One array task per slice of the worklist (rows task::ntasks); each run writes its own files
# (<split>/<run>.top<N>.asf, .calib.json, .json) -- one writer per file on Quobyte.
#   sbatch --array=0-39%20 prepare_training_native.sbatch <worklist.tsv> <outdir> 40
# Env: VARIANTS (default 150,1000), OWN_MIN (per-variant own-top-N coverage filter, default 0 = off),
#      DOMINANT (dominant-label margin, default 0 = off: random DIA-NN precursor per record).
set -u
B=/quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native
export PATH=/quobyte/proteomics-grp/envs/cascadia5/bin:$PATH
export PYTHONPATH=/quobyte/proteomics-grp/brett/radiant_v3/batmilk/cascadia/pyoverlay
WL=$1; OUT=$2; NT=$3; T=${SLURM_ARRAY_TASK_ID:?array task}
echo "task $T on $(hostname) $(date -Is)"
tail -n +2 "$WL" | awk -F'\t' -v t=$T -v n=$NT '(NR-1)%n==t' | while IFS=$'\t' read -r rk split grp sid rep d npre; do
  run=$(basename "$d" .d); o="$OUT/$split/$run"; mkdir -p "$OUT/$split"
  [ -s "$o.json" ] && { echo "have $o"; continue; }
  cap=1000; [ "$split" != train ] && cap=300          # val/test: fewer per run is plenty
  python3 -u -W ignore "$B/prepare_training_native.py" --report "$rep" --run "$run" --d "$d" --out "$o" \
      --variants "${VARIANTS:-150,1000}" --own-min-coverage "${OWN_MIN:-0}" --dominant-margin "${DOMINANT:-0}" --max-per-run $cap --seed $T \
      > "$o.log" 2>&1
  rc=$?
  echo "rc=$rc $split $run $(tail -1 "$o.log" | cut -c1-200)"
done
echo "task $T done $(date -Is)"
