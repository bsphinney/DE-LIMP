# HIVE Cluster Agent Bootstrap

**Purpose:** Context for Claude agents working on Brett Phinney's DE-LIMP / proteomics stack at the UC Davis HIVE HPC cluster. Read this once at session start so you know what's where and how to run things.

**Last verified:** 2026-04-21. Paths drift; `ls`-check before relying.

---

## 1. Access

| | |
|---|---|
| Host | `hive.hpc.ucdavis.edu` |
| Username | `brettsp` |
| Auth | Key-based only. Private key at `~/.ssh/id_ed25519` on the controlling machine. |
| ControlMaster socket | `/tmp/.delimp_brettsp_hive` (reuse for speed; macOS 104-byte limit on path length) |

**Connection template:**
```bash
ssh -o ControlMaster=auto -o ControlPath=/tmp/.delimp_brettsp_hive -o ControlPersist=300 \
    brettsp@hive.hpc.ucdavis.edu "<command>"
```

**SLURM commands need a login shell** — `~/.bashrc` alone doesn't populate the module path:
```bash
ssh brettsp@hive "bash -l -c 'sbatch my_job.sbatch'"
```

**Never run heavy computation on login nodes.** Login nodes are shared; intensive CPU/memory work gets flagged. Use `sbatch` or an interactive `srun` node.

---

## 2. SLURM (partitions, QOS, limits)

| Partition | Account | QOS | Use |
|---|---|---|---|
| `high` | `genome-center-grp` | `genome-center-grp-high-qos` | Priority CPU; 64-CPU per-user cap (`MaxTRESPU`) |
| `gpu-a100` | `genome-center-grp` | `genome-center-grp-gpu-a100-qos` | 1 A100 / 80 GB — use for Casanovo/Cascadia inference + training |
| `low` | `publicgrp` | `publicgrp-low-qos` | Preemptible, huge capacity; set `Requeue=1` when using |

**Auto queue switch logic**: `genome-center-grp/high` → `publicgrp/low` when pending >5 min and public has idle CPUs. Only array steps (2/4 in a 5-step parallel DIA-NN search) move; assembly steps stay on high. See `docs/QUEUE_SWITCHING.md` for details.

**Check state:**
```bash
squeue -u brettsp --format="%.10i %.20j %.10T %.10M %.10l %R"
sacct -j <jobid> --format=JobID,JobName%25,State,ExitCode,Elapsed
```

**Gotcha**: `sacct` `.extern` / `.batch` substeps report COMPLETED even when the main job failed. Filter with `grep -v "\\."` or check the top-level JobID only.

**GPU nodes**: Always `--gres=gpu:a100:1` for Casanovo inference. A6000 nodes (`--gres=gpu:1` without filter) are 2-3x slower and sometimes have kernel bugs.

---

## 3. Shared storage layout

```
/quobyte/proteomics-grp/                    Root for shared proteomics storage
├── dia-nn/
│   └── diann_2.3.0.sif                     DIA-NN 2.3 + .NET (reads Thermo .raw) — USE THIS ONE
├── apptainers/
│   ├── diann2.3.0.sif                      DIA-NN 2.3 without .NET — .raw files silently skipped. AVOID.
│   ├── msconvert.sif                       ProteoWizard (Thermo .raw → mzML)
│   ├── pwiz-skyline-...sif                 Alt msconvert (same purpose)
│   ├── alphadia.sif                        alphaDIA
│   └── dia-analyst_v0.10.5.sif             DIA-Analyst
├── de-limp/
│   ├── DE-LIMP/                            Shared repo clone
│   ├── containers/
│   │   └── de-limp.sif                     DE-LIMP Shiny app
│   ├── cascadia/
│   │   ├── models/
│   │   │   ├── cascadia.ckpt               Cascadia inference checkpoint
│   │   │   ├── cascadia_pretrained.ckpt
│   │   │   └── cascadia_astral_tuned.ckpt
│   │   ├── sage-v0.14.7-x86_64-unknown-linux-gnu/sage   Sage DDA search binary
│   │   └── training/                        Per-track training outputs + sbatch
│   ├── fasta/                               Pre-staged organism FASTAs (see §5)
│   ├── activity_log.csv                     Shared run history
│   └── users/<username>/                    Per-user dirs (output, logs, jobs)
├── bioinformatics_programs/
│   ├── blast_dbs/
│   │   ├── uniprot_sprot.dmnd              DIAMOND SwissProt
│   │   ├── uniprot_sprot_reversed.dmnd     Reversed decoy
│   │   ├── uniprot_trembl.dmnd             DIAMOND TrEMBL
│   │   └── uniprot_trembl_reversed.dmnd    Reversed decoy
│   └── casanovo_modles/                    (note typo — directory name is intentional)
│       ├── casanovo_v5_0_0.ckpt            Casanovo 5.0 (current default)
│       ├── casanovo_v4_2_0.ckpt            Casanovo 4.2
│       ├── casanovo_massivekb.ckpt         Nontryptic-trained
│       └── casanovo_nontryptic.ckpt
└── conda_envs/ and envs/                    See §4
```

**Per-user DE-LIMP storage**: `/quobyte/proteomics-grp/de-limp/users/${USER}/` — outputs, SLURM logs, sbatch scripts land here.

---

## 4. Software

### Apptainer (Singularity) containers

Always bind `/quobyte:/quobyte` so paths inside the container match host paths:

```bash
apptainer exec --bind /quobyte:/quobyte \
    /quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif \
    /diann-2.3.0/diann-linux [flags]
```

**DIA-NN binary is at `/diann-2.3.0/diann-linux` inside the container**, not just `diann`. Don't use `module load diann` — that command does not exist on this cluster.

### Modules (availability via `bash -l -c 'module avail'`)

| Module | Notes |
|---|---|
| `apptainer/latest` | Container runtime (`apptainer exec`) |
| `blast-plus/2.16.0` | NCBI BLAST+ (classic) |
| `diamond/2.1.7` | DIAMOND BLAST (use for protein-vs-protein in cascadia/casanovo) |
| `blast2go/5.2.5` | Functional annotation |
| `ncbi-rmblastn/2.14.0` | Repeatmasker variant |
| `python/3.11.9` | System Python |
| `R/4.3.3`, `R/4.4.2` | R (system); DE-LIMP uses R 4.5+ from Apptainer instead |

Load inline in sbatch: `module load diamond`

### Conda environments (`/quobyte/proteomics-grp/`)

```
/quobyte/proteomics-grp/envs/
├── cascadia5/       Cascadia current (IM-enabled, v3 bruker_patch)
├── cascadia4/       Prior Cascadia
├── cascadia/        Oldest
├── casanovo-gpu/    Casanovo with CUDA PyTorch
├── casanovo5/       Casanovo 5.0.0 with depthcharge 0.2.3 patched
├── instanovo_017_final/   Instanovo 0.1.7 (known-good lockfile)
└── datasci/         Generic DS env

/quobyte/proteomics-grp/conda_envs/
├── R4.5.1/          R 4.5.1 for scripting
├── casanovo5/       Parallel to envs/casanovo5 (same purpose, different install)
├── cassonovo_env/   (typo — do not rename, workflows reference it)
├── alphadia/        alphaDIA Python API
├── dotnet/          .NET for running DIA-NN outside the container
└── mztab/, mztab_peptide/    mzTab parsers
```

**Activation pattern inside sbatch:**
```bash
CONDA_ENV="/quobyte/proteomics-grp/conda_envs/casanovo5"
export PATH="${CONDA_ENV}/bin:$PATH"
```
(Avoids the heavy `conda activate` machinery.)

### Standalone binaries

| Tool | Path |
|---|---|
| Sage | `/quobyte/proteomics-grp/de-limp/cascadia/sage-v0.14.7-x86_64-unknown-linux-gnu/sage` |
| Sage master (source build) | `/quobyte/proteomics-grp/de-limp/cascadia/sage-master-bin/sage` |

**Sage version note**: v0.14.7 doesn't have native protein grouping — post-hoc via `sage_protein_groups.py`. Sage > DIA-NN on timsTOF DDA by ~2x PSMs; comparable on Orbitrap.

---

## 5. Reference data

### FASTAs (pre-staged at `/quobyte/proteomics-grp/de-limp/fasta/` and MRS/)

| Organism | Path |
|---|---|
| Human (UP000005640) | `/quobyte/proteomics-grp/MRS/UP000005640_9606.fasta` |
| Human + universal contaminants | `/quobyte/proteomics-grp/MRS/UP000005640_9606_plus_universal_contam.fasta` |
| Bovine (UP000009136) | `/quobyte/proteomics-grp/de-limp/fasta/UP000009136_bos_taurus.fasta` |
| Chicken (UP000000539) | `/quobyte/proteomics-grp/de-limp/fasta/UP000000539_gallus_gallus.fasta` |
| Porcine (UP000008227) | `/quobyte/proteomics-grp/de-limp/fasta/UP000008227_sus_scrofa.fasta` |

### BLAST databases

`/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/`
- `uniprot_sprot.dmnd` — DIAMOND SwissProt (target)
- `uniprot_sprot_reversed.dmnd` — decoy for TD-FDR
- `uniprot_trembl.dmnd` — larger, noisier
- `uniprot_trembl_reversed.dmnd` — decoy

---

## 6. Common workflows (with commands)

### DIA-NN search (single-file, via Apptainer)

```bash
apptainer exec --bind /quobyte:/quobyte \
    /quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif \
    /diann-2.3.0/diann-linux \
    --f /path/to/file.raw --fasta /path/to/fasta --out report.parquet \
    --threads 16 --mass-acc 10 --mass-acc-ms1 10
```

### DIA-NN parallel search (5-step)

Generated by DE-LIMP's `generate_parallel_scripts()` in `R/helpers_search.R`. Steps: (1) library prediction → (2) per-file quant array → (3) assembly → (4) second quant array → (5) cross-run report. Queue-switch on pending >5 min.

### Sage DDA search (example sbatch fragment)

```bash
#SBATCH --partition=high --account=genome-center-grp --qos=genome-center-grp-high-qos
#SBATCH --cpus-per-task=16 --mem=64G --time=4:00:00
export PATH="/quobyte/proteomics-grp/de-limp/cascadia/sage-v0.14.7-x86_64-unknown-linux-gnu:$PATH"
sage sage_config.json
python /quobyte/proteomics-grp/de-limp/cascadia/sage_protein_groups.py results.sage.parquet
```

### Casanovo inference (GPU)

```bash
#SBATCH --partition=gpu-a100 --gres=gpu:a100:1
#SBATCH --cpus-per-task=8 --mem=64G --time=4:00:00
export PATH="/quobyte/proteomics-grp/conda_envs/casanovo5/bin:$PATH"
casanovo sequence \
    --model /quobyte/proteomics-grp/bioinformatics_programs/casanovo_modles/casanovo_v5_0_0.ckpt \
    input.mgf -o output
```

**Casanovo eval gotcha**: `casanovo --evaluate` silently reports 0% due to tokenizer vocab mismatch on our data. Use `score_casanovo.py` custom scorer (in `cascadia/` subdir). Ground truth must use species-matched FASTA — human FASTA on feather spectra produces bogus 58% baseline.

### DIAMOND BLAST

```bash
module load diamond/2.1.7
diamond blastp \
    --db /quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/uniprot_sprot.dmnd \
    --query peptides.fasta --out hits.tsv \
    --outfmt 6 qseqid sseqid pident evalue stitle \
    --threads 16 --max-target-seqs 5 --evalue 1e-5
```

For de-novo peptide target-decoy FDR, run the same query against `uniprot_sprot_reversed.dmnd` and compute FDR = decoy_hits / target_hits.

### Cascadia de novo (GPU, IM-enabled)

```bash
#SBATCH --partition=gpu-a100 --gres=gpu:a100:1
export PATH="/quobyte/proteomics-grp/envs/cascadia5/bin:$PATH"
cascadia predict input.d \
    --checkpoint /quobyte/proteomics-grp/de-limp/cascadia/models/cascadia.ckpt \
    --output output.ssl
```

For native `.d` reading, Cascadia needs a routing patch to auto-detect Bruker input and use `bruker_augment.py`. See `cascadia/TRAINING_LOG.md` for status.

### Fetching log / job info

```bash
# Job stdout/stderr
cat /quobyte/proteomics-grp/de-limp/cascadia/training/casanovo_track_g/logs/track_g_<jobid>.out

# Current queue
squeue -u brettsp --format="%.10i %.20j %.10T %.10M %.10l %R" | grep -v stan-watcher

# Monitor active job output
tail -f .../logs/<jobid>.out
```

---

## 7. Gotchas (ordered by frequency of tripping)

| # | Thing | Why it matters |
|---|---|---|
| 1 | Two DIA-NN containers | `/dia-nn/diann_2.3.0.sif` (has .NET, reads `.raw`) vs `/apptainers/diann2.3.0.sif` (no .NET, `.raw` silently skipped). ALWAYS use the first. |
| 2 | DIA-NN binary path inside container | `/diann-2.3.0/diann-linux`, not `diann`. |
| 3 | SLURM needs login shell | `ssh host "bash -l -c '...'"` — otherwise `sbatch`/`sacct` not found. |
| 4 | `sacct` substeps lie | `.extern` / `.batch` report COMPLETED even when parent failed. Check top-level JobID only. |
| 5 | `module load diann` doesn't exist | DIA-NN is a container, not a module. |
| 6 | Per-user CPU cap is 64 on `high` | `MaxTRESPU`. Jobs over 64 get `InvalidQOS`. Auto-switch to `publicgrp/low` kicks in. |
| 7 | `casanovo --evaluate` is broken | Silently returns 0% on our data. Use `score_casanovo.py` custom scorer. |
| 8 | Cascadia default `temp_path` fills home quota | Always set `--temp /tmp/...` in training scripts. |
| 9 | `casanovo_modles/` directory is typo'd | Do NOT rename — scripts across the codebase reference the exact path. |
| 10 | `cassonovo_env` is also typo'd | Same as above. |
| 11 | Home directory has a quota (~10 GB) | Use `/quobyte/proteomics-grp/` for anything meaningful. Startup warning fires at >80% use. |
| 12 | ddaPASEF/diaPASEF MS1 frame reading | Use Sage's formula from `timsrust/src/converters.rs`. Formula-from-first-principles was off by 155 Da. |
| 13 | PyTorch Lightning `trainer.fit(ckpt_path=...)` hangs with non-resumable DataLoader | Use `from_pretrained()` to load weights only; accept fresh optimizer momentum. |
| 14 | DIA-NN requires `--quant-ori-names` on ALL parallel steps | Otherwise `.quant` filenames mismatch between container bind mounts. |
| 15 | Manual mass accuracy when using `--use-quant` | Auto-mode + `--use-quant` produces different results from original run. Always set `--mass-acc`/`--mass-acc-ms1`/`--window` explicitly. |

---

## 8. Cross-references

- `CLAUDE.md` — project context (R/Shiny architecture, UI patterns)
- `docs/PATTERNS.md` — R Shiny / bslib gotchas
- `docs/GOTCHAS.md` — broader gotchas table (70+ items)
- `docs/HPC_PATHS.md` — source of truth for paths (re-verified against this file)
- `docs/QUEUE_SWITCHING.md` — partition auto-switch logic
- `HPC_DEPLOYMENT.md` — Apptainer deployment walkthrough
- `WINDOWS_WSL_INSTALL.md` — WSL install (SSH key setup flows to HIVE)

---

## 9. Where to write (permission & etiquette)

- **Group-writable shared paths** (`/quobyte/proteomics-grp/de-limp/...`): fine to write to `users/<yourname>/`, `activity_log.csv`, training output dirs you own.
- **Read-only-ish**: models, containers, FASTA library, BLAST DBs — don't modify without coordination.
- **Temp compute scratch**: `/tmp/<anything>` per node, local to that node only. Use `SLURM_TMPDIR` if available, or `/tmp/casanovo_<SLURM_JOB_ID>`.
- Log files: write under the sbatch's working dir + `/logs/` subdir (DE-LIMP convention). Create with `mkdir -p` before writing or processx/SLURM will fail cryptically.

---

## 10. If in doubt

**Always verify paths before acting.** Do `ssh ... "ls /path"` or `find` before assuming a binary/container/dataset exists. Paths drift. Containers get rebuilt. Trust `ls` output over this document's contents.

If a memory claim conflicts with a `ls`, trust the `ls`.
