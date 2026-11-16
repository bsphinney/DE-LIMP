#!/bin/bash
#SBATCH --job-name=casc_train
#SBATCH --account=genome-center-grp
#SBATCH --partition=gpu-a100
#SBATCH --qos=genome-center-grp-gpu-a100-qos
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --output=/quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native/train_v1/logs/train_%j.out
# Fine-tune the stock Cascadia checkpoint on train_v1 (native-reader labelled spectra).
#   sbatch train_native.sbatch <train_v1 dir> [epochs]
set -u
D=$1; E=${2:-10}
export PATH=/quobyte/proteomics-grp/envs/cascadia5/bin:$PATH
export PYTHONPATH=/quobyte/proteomics-grp/brett/radiant_v3/batmilk/cascadia/pyoverlay
CKPT=/quobyte/proteomics-grp/de-limp/cascadia/models/cascadia.ckpt
cd "$D"                                  # cascadia train writes checkpoint_<date>/ and logs/ here
cat train/*.asf > train.asf; cat val/*.asf > val.asf
echo "train spectra $(grep -c '^BEGIN IONS' train.asf), val $(grep -c '^BEGIN IONS' val.asf); node $(hostname) $(date -Is)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
cascadia train train.asf val.asf -m "$CKPT" -b 32 -e "$E" -lr 1e-5 -c 4
echo "train exit=$? $(date -Is)"; ls -d checkpoint_*/ 2>/dev/null | tail -1
