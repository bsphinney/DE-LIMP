# HPC Paths Reference (HIVE at UC Davis)

**IMPORTANT**: Always verify paths with `ls`/`find` on the cluster before using. Do NOT rely on this file alone.

## Connection
- **Host**: hive.hpc.ucdavis.edu
- **User**: brettsp
- **SSH key**: ~/.ssh/id_ed25519
- **SLURM account**: genome-center-grp
- **Partitions**: `high` (CPU), `gpu-a100` (GPU), `publicgrp/low` (preemptible)
- **QOS**: `genome-center-grp-high-qos`, `genome-center-grp-gpu-a100-qos`
- **Per-user CPU limit**: 64 CPUs on high partition (MaxTRESPU)

## Containers (Apptainer)

| Container | Path | Notes |
|-----------|------|-------|
| **DIA-NN 2.3 (with Thermo .raw support)** | `/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif` | Has .NET runtime, reads .raw + .d + .mzML |
| DIA-NN 2.3 (Bruker only, NO .raw) | `/quobyte/proteomics-grp/apptainers/diann2.3.0.sif` | Missing dotnet — `.raw` files silently skipped |
| msconvert (ProteoWizard) | `/quobyte/proteomics-grp/apptainers/pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif` | `wine64 msconvert file.raw --mzML` |
| alphaDIA | `/quobyte/proteomics-grp/apptainers/alphadia.sif` | |
| DE-LIMP | `/quobyte/proteomics-grp/de-limp/containers/de-limp.sif` | |

**DIA-NN binary inside container**: `/diann-2.3.0/diann-linux` (NOT just `diann`)

**DIA-NN run command**:
```bash
apptainer exec --bind /quobyte:/quobyte \
  /quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif \
  /diann-2.3.0/diann-linux [flags]
```

**CRITICAL**: There are TWO different DIA-NN containers:
- `/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif` — has .NET, reads Thermo .raw
- `/quobyte/proteomics-grp/apptainers/diann2.3.0.sif` — NO .NET, .raw files fail with "dotnet: not found"

## FASTA Files

| Species | Path |
|---------|------|
| Human (HeLa) | `/quobyte/proteomics-grp/MRS/UP000005640_9606.fasta` |
| Human + contaminants | `/quobyte/proteomics-grp/MRS/UP000005640_9606_plus_universal_contam.fasta` |
| Bovine | `/quobyte/proteomics-grp/de-limp/fasta/UP000009136_bos_taurus.fasta` |
| Chicken | `/quobyte/proteomics-grp/de-limp/fasta/UP000000539_gallus_gallus.fasta` |
| Porcine | `/quobyte/proteomics-grp/de-limp/fasta/UP000008227_sus_scrofa.fasta` |

## BLAST Databases (DIAMOND)

| Database | Path |
|----------|------|
| SwissProt | `/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/uniprot_sprot` |
| TrEMBL | `/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/uniprot_trembl` |

## Storage

| Purpose | Path |
|---------|------|
| Shared group storage | `/quobyte/proteomics-grp/de-limp/` |
| Per-user output | `/quobyte/proteomics-grp/de-limp/{username}/output/` |
| Pre-staged FASTA | `/quobyte/proteomics-grp/de-limp/fasta/` |
| Downloads | `/quobyte/proteomics-grp/de-limp/downloads/` |
| Cascadia training | `/quobyte/proteomics-grp/de-limp/cascadia/training/` |
| Cascadia env | `/quobyte/proteomics-grp/envs/cascadia5/` |
| Casanovo env | `/quobyte/proteomics-grp/conda_envs/cassonovo_env/` (note typo) |
| Sage binary | `/quobyte/proteomics-grp/de-limp/cascadia/sage-v0.14.7-x86_64-unknown-linux-gnu/sage` |
| Cascadia model | `/quobyte/proteomics-grp/de-limp/cascadia/models/cascadia.ckpt` |
| Casanovo model | `/quobyte/proteomics-grp/bioinformatics_programs/casanovo_modles/casanovo_v4_2_0.ckpt` (note typo) |

## SLURM Notes
- SLURM tools need login shell: `bash -l -c '...'`
- DIA-NN is NOT a module — `module load diann` does not work
- `sacct` `.extern`/`.batch` substeps report COMPLETED even when main job failed — filter with `grep -v "\\."`
