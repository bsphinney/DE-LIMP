#!/bin/bash
#SBATCH --job-name=casc_train
#SBATCH --account=genome-center-grp
#SBATCH --partition=gpu-a100
#SBATCH --qos=genome-center-grp-gpu-a100-qos
#SBATCH --gres=gpu:1
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --output=/quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native/train_v1/logs/train_%j.out
# Fine-tune the stock Cascadia checkpoint on one peak-budget variant of train_v1.
#   sbatch train_native.sbatch <train_v1 dir> <top_n variant: 150|1000> [epochs]   (env JIT=ppm, default 15)
# Writes <dir>/model_top<N>/checkpoint_<date>/Cascadia-epoch=..-step=...ckpt (best 2 by Val Pep. Acc.).
# -c 10: the stock checkpoint's charge embedding is sized for max_charge 10 (inference loads it so).
set -u
D=$1; V=$2; E=${3:-10}
export PATH=/quobyte/proteomics-grp/envs/cascadia5/bin:$PATH
export PYTHONPATH=/quobyte/proteomics-grp/brett/radiant_v3/batmilk/cascadia/pyoverlay
CKPT=/quobyte/proteomics-grp/de-limp/cascadia/models/cascadia.ckpt
M=$D/model_top$V; mkdir -p "$M"; cd "$M"          # cascadia train writes checkpoint_<date>/ and logs/ here
# drift augmentation: one random +/-JIT ppm m/z offset per training record (val stays as built)
JIT=${JIT:-15}
[ -s train.asf ] || python3 /quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native/jitter_asf.py train.asf "$JIT" 0 "$D"/train/*.top$V.asf
rm -f train_gpu.hdf5 val_gpu.hdf5 train*.hdf5 val*.hdf5   # an index left by a failed start
[ -s val.asf ]   || cat "$D"/val/*.top$V.asf > val.asf
echo "variant top$V: train $(grep -c '^BEGIN IONS' train.asf), val $(grep -c '^BEGIN IONS' val.asf); node $(hostname) $(date -Is)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python3 -u /quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native/cascadia_train_wrap.py train train.asf val.asf -m "$CKPT" -b 32 -e "$E" -lr 1e-5 -c 10 2>&1 | grep --line-buffered -v 'lines/s\]$'
echo "train exit=${PIPESTATUS[0]} $(date -Is)"; ls -l checkpoint_*/ 2>/dev/null
