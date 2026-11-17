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
#   sbatch --array=0-N eval_models.sbatch RUNS.tsv OUTDIR STOCK_CKPT FT150_CKPT FT1000_CKPT
set -u
B=/quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native
export PATH=/quobyte/proteomics-grp/envs/cascadia5/bin:$PATH
export PYTHONPATH=/quobyte/proteomics-grp/brett/radiant_v3/batmilk/cascadia/pyoverlay
RUNS=$1; OUT=$2; STOCK=$3; F150=$4; F1000=$5
IFS=$'\t' read -r N D REP CAL < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$RUNS")
TMP=${TMPDIR:-/tmp}/ce_${SLURM_JOB_ID}; mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
O="$OUT/$N"; mkdir -p "$O"
echo "run $N node $(hostname) $(date -Is)"
run() { local tag=$1; shift; python3 -W ignore $B/cascadia_native_run.py "$@" "$D" "$N" "$O/$tag" "$TMP/$tag" 2>&1 | grep -E '^\[|native augment|top_n'; rm -rf "$TMP/$tag"; }
run t150      --top-n 150                stock=$STOCK ft150=$F150
run t150cal   --top-n 150  --calib "$CAL" stock=$STOCK ft150=$F150
run t1000     --top-n 1000               ft1000=$F1000
run t1000cal  --top-n 1000 --calib "$CAL" ft1000=$F1000
python3 $B/score_vs_diann.py "$REP" "$N" \
  stock=$O/t150/stock/$N.ssl stock_cal=$O/t150cal/stock/$N.ssl \
  ft150=$O/t150/ft150/$N.ssl ft150_cal=$O/t150cal/ft150/$N.ssl \
  ft1000=$O/t1000/ft1000/$N.ssl ft1000_cal=$O/t1000cal/ft1000/$N.ssl | tee "$O/score.txt"
echo "done $(date -Is)"
