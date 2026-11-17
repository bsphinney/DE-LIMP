#!/bin/bash
#SBATCH --job-name=casc_eval
#SBATCH --account=genome-center-grp
#SBATCH --partition=gpu-a100
#SBATCH --qos=genome-center-grp-gpu-a100-qos
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=/quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native/train_v1/logs/eval_%A_%a.out
# Score Cascadia checkpoints on held-out runs against DIA-NN. One array task per row of RUNS.tsv
# (run_name, d_path, report.parquet, calib.json). Each checkpoint runs at the peak budget it was
# trained on, with and without the run's m/z calibration:
#   sbatch --array=0-N eval_models.sbatch RUNS.tsv OUTDIR STOCK_CKPT FT150_CKPT FT1000_CKPT [BASELINE_DIR]
# Pass "-" for any checkpoint to skip its arms. BASELINE_DIR is an earlier OUTDIR whose stock arms
# are scored alongside (so the stock model need not be re-run for every new checkpoint).
set -u
B=/quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native
export PATH=/quobyte/proteomics-grp/envs/cascadia5/bin:$PATH
export PYTHONPATH=/quobyte/proteomics-grp/brett/radiant_v3/batmilk/cascadia/pyoverlay
RUNS=$1; OUT=$2; STOCK=$3; F150=$4; F1000=$5; BASE=${6:-}
IFS=$'\t' read -r N D REP CAL < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$RUNS")
TMP=${TMPDIR:-/tmp}/ce_${SLURM_JOB_ID}; mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
O="$OUT/$N"; mkdir -p "$O"
echo "run $N node $(hostname) $(date -Is)"
# run TAG "RUNNER OPTIONS" name=ckpt...   (the runner wants: [options] run.d name outdir tmpdir name=ckpt...)
run() { local tag=$1 opts=$2; shift 2
  python3 -W ignore $B/cascadia_native_run.py $opts "$D" "$N" "$O/$tag" "$TMP/$tag" "$@" > "$O/$tag.log" 2>&1
  echo "$tag exit=$? $(date -Is)"; grep -E '^\[|native augment|top_n|Error|Traceback' "$O/$tag.log" | tail -5; rm -rf "$TMP/$tag"; }
M150=""; [ "$STOCK" != - ] && M150="stock=$STOCK"; [ "$F150" != - ] && M150="$M150 ft150=$F150"
if [ -n "$M150" ]; then
  run t150    "--top-n 150"                 $M150
  run t150cal "--top-n 150 --calib $CAL"    $M150
fi
if [ "$F1000" != - ]; then
  run t1000    "--top-n 1000"              ft1000=$F1000
  run t1000cal "--top-n 1000 --calib $CAL" ft1000=$F1000
fi
arms=()
for a in stock:t150/stock stock_cal:t150cal/stock ft150:t150/ft150 ft150_cal:t150cal/ft150 \
         ft1000:t1000/ft1000 ft1000_cal:t1000cal/ft1000; do
  f=$O/${a#*:}/$N.ssl
  [ ! -s "$f" ] && [ -n "$BASE" ] && [[ ${a%%:*} == stock* ]] && f=$BASE/$N/${a#*:}/$N.ssl
  [ -s "$f" ] && arms+=("${a%%:*}=$f")
done
python3 $B/score_vs_diann.py "$REP" "$N" "${arms[@]}" | tee "$O/score.txt"
python3 $B/rank_score.py "$REP" "$N" "${arms[@]}" | tee "$O/rank.txt"
echo "done $(date -Is)"
