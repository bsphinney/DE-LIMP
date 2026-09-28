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

## DIA-NN (verified 2026-09-28)

**The skill pins DIA-NN 2.7.0**, a native Linux binary, not a container:
`/quobyte/proteomics-grp/dia-nn/build_270/diann-2.7.0/diann-linux` (`acquire_tools.sh hpc`
resolves it). Reading Thermo `.raw` with it needs a .NET runtime ≥ 8.0.17 on `DOTNET_ROOT`
(`ensure_dotnet8.sh`; the skill sets it). Earlier builds sit beside it (`build_251`,
`build_260`, `build_261`). The 2.3.0 containers below are what the DE-LIMP app's HIVE path
still uses — not the skill.

## Containers (Apptainer)

| Container | Path | Notes |
|-----------|------|-------|
| DIA-NN 2.3 (with Thermo .raw support) — DE-LIMP app, not the skill | `/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif` | Has .NET runtime, reads .raw + .d + .mzML |
| DIA-NN 2.3 (Bruker only, NO .raw) | `/quobyte/proteomics-grp/apptainers/diann2.3.0.sif` | Missing dotnet — `.raw` files silently skipped |
| msconvert (ProteoWizard) | `/quobyte/proteomics-grp/apptainers/pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif` | `wine64 msconvert file.raw --mzML` |
| alphaDIA | `/quobyte/proteomics-grp/apptainers/alphadia.sif` | |
| DE-LIMP | `/quobyte/proteomics-grp/de-limp/containers/de-limp.sif` | |

**DIA-NN binary inside the 2.3.0 container**: `/diann-2.3.0/diann-linux` (NOT just `diann`)

**DIA-NN 2.3.0 container run command** (DE-LIMP app):
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
| Human + contaminants (current, 2026-09) | `/quobyte/proteomics-grp/MRS/UP000005640_9606_plus_universal_contam_2026-09.fasta` — UniProt 2026_03 one-per-gene (20,652) + 220 `Cont_`-tagged Universal contaminants (Frankenfield 2022), built by skill 2.8.0 `fetch_fasta.py`, which removed the 161 contaminant entries that are, or cannot be told apart from, human proteins. Sidecar `.fasta.meta.json` beside it; predicted library `UP000005640_9606_plus_universal_contam_2026-09.predicted.speclib` (DIA-NN 2.7.0). Build scripts + logs: `/quobyte/proteomics-grp/claude/mrs_rebuild_2026-09-25/` |
| Human + contaminants (Sep 2025, superseded) | `/quobyte/proteomics-grp/MRS/UP000005640_9606_plus_universal_contam.fasta` + its DIA-NN 2.2.0 `.predicted.speclib` — kept only because old searches point at them. Its 381 `Cont_` entries include 153 identical to human proteins (bovine ACTB/EEF1A1/tubulins, human keratins), which DIA-NN then reports only as `Cont_` groups and keeps out of quantification. Do not use it for new searches. |
| Bovine | **none usable** — `/quobyte/proteomics-grp/de-limp/fasta/UP000009136_bos_taurus.fasta` is 0 bytes (verified 2026-09-28); build one with `fetch_fasta.py fetch --proteome UP000009136` |
| Chicken | `/quobyte/proteomics-grp/de-limp/fasta/UP000000539_gallus_gallus.fasta` |
| Porcine | `/quobyte/proteomics-grp/de-limp/fasta/UP000008227_sus_scrofa.fasta` |

`fetch_fasta.py --hive` reuses only `/quobyte/proteomics-grp/MRS/`; the `de-limp/fasta/` files
are not searched by it.

## BLAST Databases (DIAMOND)

All in `/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/` (verified 2026-09-28).
DIAMOND's `-d` takes the path with or without `.dmnd`; the file on disk has it.

| Database | Path |
|----------|------|
| SwissProt | `blast_dbs/uniprot_sprot.dmnd` (FASTA: `uniprot_sprot.fasta`; reversed decoys: `uniprot_sprot_reversed.dmnd`) |
| TrEMBL | `blast_dbs/uniprot_trembl.dmnd` (FASTA: `uniprot_trembl.fasta`; reversed decoys: `uniprot_trembl_reversed.dmnd`) |
| NCBI nr | `blast_dbs/ncbi_nr/nr.dmnd` (375 GB; taxonomy `names.dmp` / `nodes.dmp` beside it) |

## Storage

| Purpose | Path |
|---------|------|
| Shared group storage | `/quobyte/proteomics-grp/de-limp/` |
| Per-user output | `/quobyte/proteomics-grp/de-limp/{username}/output/` |
| Pre-staged FASTA | `/quobyte/proteomics-grp/de-limp/fasta/` |
| Downloads | `/quobyte/proteomics-grp/de-limp/downloads/` |
| Cascadia training | `/quobyte/proteomics-grp/de-limp/cascadia/training/` |
| Cascadia env | `/quobyte/proteomics-grp/envs/cascadia5/` |
| Casanovo v4 env | `/quobyte/proteomics-grp/conda_envs/cassonovo_env/` (typo `casso`; Casanovo 4.3.0, Python 3.10, depthcharge-ms ~0.2.x). Use with `casanovo_v4_2_0.ckpt`. |
| Casanovo v5 env | `/quobyte/proteomics-grp/conda_envs/casanovo5/` (Casanovo 5.0.0, Python 3.13, depthcharge-ms 0.4.8 with `depthcharge.tokenizers`). Required for `casanovo_v5_0_0.ckpt`. v4 env cannot load v5 ckpts (`ModuleNotFoundError: depthcharge.tokenizers`). |
| Sage binary | `/quobyte/proteomics-grp/de-limp/cascadia/sage-v0.14.7-x86_64-unknown-linux-gnu/sage` |
| Cascadia model | `/quobyte/proteomics-grp/de-limp/cascadia/models/cascadia.ckpt` |
| Casanovo model | `/quobyte/proteomics-grp/bioinformatics_programs/casanovo_modles/casanovo_v4_2_0.ckpt` (note typo) |

## SLURM Notes
- SLURM tools need login shell: `bash -l -c '...'`
- DIA-NN is NOT a module — `module load diann` does not work
- `sacct` `.extern`/`.batch` substeps report COMPLETED even when main job failed — filter with `grep -v "\\."`
