# Cascadia timsTOF Fine-Tuning Plan

## Current Model Architecture

Cascadia's transformer encoder takes per-peak features as a **4-tuple**:
```
(m/z, intensity, retention_time, ms_level)
```

Ion mobility (1/K0) is **NOT** an input feature. When processing timsTOF data, the IM dimension is collapsed (summed across mobility frames) during augmentation, producing standard MS2 spectra. This means:

- The model sees the same data as from an Orbitrap (m/z + intensity)
- IM information is lost — peptides with different CCS but same m/z are merged
- This is the same approach as mzML conversion (which also collapses IM)

## Option A: Fine-Tune Without IM (Recommended First Step)

Fine-tune the existing model on timsTOF data without architectural changes. This teaches the model timsTOF-specific characteristics:

- Different fragmentation patterns (CID in PASEF vs HCD in Orbitrap)
- Different mass accuracy profiles
- Different intensity distributions
- Different noise characteristics

### Training Data Sources

**Priority 1: Own data (best instrument match)**
- HeLa affinisep/evosep files (12 files, timsTOF HT, 100 SPD)
- Ground truth: DIA-NN identifications at 1% FDR
- Location: `/quobyte/proteomics-grp/service/off_campus/University_of_Oaklahoma/Vinning-paul/Liver/`

**Priority 2: Public benchmarks**
| PRIDE ID | Description | Instrument |
|----------|-------------|------------|
| PXD014777 | Bruker timsTOF HeLa benchmark | timsTOF Pro |
| PXD038782 | MassIVE-KB timsTOF | timsTOF various |
| PXD010012 | timsTOF PASEF technical replicates | timsTOF Pro |

**Priority 3: Synthetic peptides (gold standard)**
| PRIDE ID | Description | Notes |
|----------|-------------|-------|
| PXD004947 | ProteomeTools synthetic library | Known sequences, but Orbitrap |

### Training Procedure

```bash
# 1. Prepare training data
#    - Run DIA-NN on training files (already done for HeLa)
#    - Extract peptide-spectrum matches from DIA-NN report
#    - Convert to Cascadia's training format (labeled ASF)

# 2. Split data: 80% train, 10% validation, 10% test

# 3. Fine-tune from pretrained checkpoint (transfer learning)
cascadia train \
    train_data.asf \
    val_data.asf \
    --model /quobyte/proteomics-grp/de-limp/cascadia/models/cascadia.ckpt \
    --batch-size 64 \
    --max-epochs 10 \
    --learning-rate 1e-5 \
    --output timstof_finetuned.ckpt

# 4. Evaluate on held-out test set
cascadia sequence test_data.mzML timstof_finetuned.ckpt -o test_results -t 0.5
# Compare peptide overlap with DIA-NN ground truth
```

### Expected Outcomes
- Improved peptide yield on timsTOF data (from ~3 to hundreds/thousands per file)
- Better score calibration for timsTOF fragmentation patterns
- ~1-2 hours training on single A100

### SLURM Resources
```
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
```

## Option B: Add Ion Mobility (Future — Requires Model Changes)

Extend the peak feature vector from 4-tuple to **5-tuple**:
```
(m/z, intensity, retention_time, ms_level, ion_mobility)
```

### Architecture Changes Required
1. **PeakEncoder**: Add IM encoder (FloatEncoder for continuous 1/K0 values)
2. **Combiner**: Change from `3 * d_model` to `4 * d_model` input
3. **ASF format**: Add 5th column for IM values
4. **augment.py / bruker_augment.py**: Preserve IM per-peak instead of collapsing

### Benefits
- CCS adds orthogonal information for sequence prediction
- Can distinguish co-eluting peptides with different shapes
- Better handling of diaPASEF variable windows

### Challenges
- Need to retrain from scratch (architecture change breaks checkpoint compatibility)
- Much larger training dataset needed
- No existing timsTOF + IM de novo training data
- IM values need calibration across instruments

### Recommendation
Start with Option A. If fine-tuning gives good results, evaluate whether IM integration is worth the ~10x additional effort. Option A should yield significant improvement since the main issue is likely fragmentation pattern mismatch, not missing IM information.

## Data Preparation Script (TODO)

```python
# prepare_training_data.py
# Reads DIA-NN report.parquet + raw .d files
# Extracts high-confidence peptide-spectrum matches
# Writes labeled ASF for Cascadia training

# Input: DIA-NN report.parquet, raw .d files
# Output: train.asf, val.asf, test.asf
# Filter: adj.P.Val < 0.01, peptide confidence > 0.99
```

## Timeline

1. **Week 1**: Run mzML comparison test (in progress), confirm augmentation is correct
2. **Week 2**: Prepare training data from HeLa files, download PRIDE benchmarks
3. **Week 3**: Fine-tune model (Option A), evaluate on held-out test set
4. **Week 4**: If successful, integrate fine-tuned model into DE-LIMP, benchmark against DIA-NN

## Success Criteria

- Fine-tuned model identifies >50% of DIA-NN peptides (de novo confirmation)
- >500 peptides per file at score > 0.8
- Peptide overlap between native .d and mzML paths > 90%
