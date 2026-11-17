#!/bin/bash
#SBATCH --job-name=im_cov
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --qos=genome-center-grp-high-qos
#SBATCH --exclude=hive-as-11-2-41
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native/logs/im_cov_%A_%a.out
B=/quobyte/proteomics-grp/brett/scratch/training_pilot_20261112/scan_lane/bruker_native
export PYTHONPATH=/quobyte/proteomics-grp/brett/radiant_v3/batmilk/cascadia/pyoverlay
IFS=$'\t' read -r N D REP CAL < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$1")
/quobyte/proteomics-grp/envs/cascadia5/bin/python3 -W ignore $B/im_cov.py "$REP" "$N" "$D" "$CAL" 800
