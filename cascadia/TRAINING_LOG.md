# Cascadia timsTOF Fine-Tuning — Training Log

## Training Run 1: Initial Fine-Tune
**Date**: March 31, 2026
**Base model**: cascadia.ckpt (pretrained on Orbitrap DIA)
**Goal**: Adapt to timsTOF HT DIA-PASEF fragmentation patterns

### Training Data

| Dataset | Files | Organism | Instrument | DIA-NN Report |
|---------|-------|----------|------------|---------------|
| HeLa affinisep/evosep | 12 | Human | timsTOF HT | ✅ |
| Berger DIA | 8 | Porcine (Sus scrofa) | timsTOF (2019) | ✅ |
| Kim | 6 | Human | timsTOF | ✅ |
| Freja | 6 | Human | timsTOF | ✅ |
| Zhao | 4 | California mouse (Peromyscus californicus) | timsTOF HT | ✅ |
| **Subtotal (available now)** | **36** | **4 species** | | |
| Liver (Oklahoma) | 60 | Bovine (Bos taurus) | timsTOF HT | Pending |
| Muscle (Oklahoma) | 60 | Bovine (Bos taurus) | timsTOF HT | Pending |
| **Total (after searches complete)** | **156** | **4 species** | | |

**Excluded**: All DDA data (Berger DDA, xlink files) — Cascadia is DIA-specific.

### Data Split
- **Training**: ~80% (~125 files) — all species represented
- **Validation**: ~10% (~16 files) — held-out subset of each species
- **Test**: ~10% (~15 files) — completely held out, never seen during training

### Data Preparation Steps

1. For each dataset with report.parquet:
   - Extract high-confidence peptide-spectrum matches (q-value < 0.01)
   - Match peptides to DIA-NN precursor IDs
   - Augment .d files to ASF format (using our native bruker_augment or mzML path)
   - Label ASF entries with matched peptide sequences

2. Combine all labeled ASF files into train/val/test splits

### Training Parameters
```
cascadia train train.asf val.asf \
    --model /quobyte/proteomics-grp/de-limp/cascadia/models/cascadia.ckpt \
    --batch-size 64 \
    --max-epochs 10 \
    --learning-rate 1e-5 \
    --width 2 \
    --max-charge 4
```

### SLURM Resources
```
Partition: gpu-a100
GPU: 1x A100 80GB
CPUs: 8
Memory: 64 GB
Time: 4 hours (estimate)
```

### Pre-Training Baseline
- Native .d path on Liver file: **3 peptides** (score > 0.5)
- mzML path on same file: **TBD** (job 11501782 running)
- These numbers establish the "before" for comparison

### Results
*(to be filled after training)*

### Notes
- Starting with 36 files (available now), will retrain with 156 when Liver/Muscle complete
- All data is timsTOF DIA-PASEF from the same instrument model
- Berger data is from 2019 (older timsTOF Pro) — different instrument generation
- Ion mobility dimension is collapsed (not used as input feature)
