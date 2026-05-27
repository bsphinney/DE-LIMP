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

### How Training Data Preparation Works

The goal: teach Cascadia "when you see THIS fragment pattern → predict THIS sequence."

**Step 1: Get ground truth labels from DIA-NN**
- Read `report.parquet` — contains all peptide identifications with sequences, precursor m/z, RT, confidence
- Filter to q-value < 0.01 (1% FDR) — only high-confidence identifications as labels

**Step 2: Match labels to raw spectra**
- Read each `.d` file with timsrust_pyo3 (native Bruker reader)
- For each DIA-NN peptide, find the MS2 spectrum that matches its precursor m/z and RT
- Also grab nearby MS1 spectra for precursor mass context

**Step 3: Write labeled ASF (Augmented Spectrum File)**
```
BEGIN IONS
PEPMASS=524.2841
CHARGE=2
RT=1523.5
SEQ=VGAHAGEYGAEALER    ← DIA-NN identification (the "answer")
234.1092  450  0.0  2   ← fragment m/z, intensity, delta_RT, MS_level
345.2183  1200  0.0  2  ← these are the actual fragment ions (the "question")
456.3274  890  0.0  2
523.8901  200  0.1  1   ← MS1 peak near precursor (level=1)
END IONS
```

During inference, `SEQ=` is a dummy placeholder. During training, it's the real peptide — the model learns to predict it from the fragment pattern.

**Step 4: Combine and split**
- Combine all labeled ASFs into train/val files
- Use different organisms for validation (tests generalization)

### Data Preparation Script
`cascadia/prepare_training_data.py` — reads DIA-NN report + raw .d files, writes labeled ASF.

```bash
python prepare_training_data.py \
    --report /path/to/report.parquet \
    --raw-dir /path/to/raw_files/ \
    --output /path/to/train.asf \
    --q-threshold 0.01
```

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
