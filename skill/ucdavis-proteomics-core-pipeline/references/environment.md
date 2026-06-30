# Environment & acquisition reference

How the skill adapts to where it's running, plus FASTA resolution.

## Platform classes (`detect_env.sh`)

| class | how detected | engine acquisition | execution |
|---|---|---|---|
| `hpc` | `sbatch` on PATH, or `/quobyte/proteomics-grp` visible | reuse existing Apptainer `.sif` (DIA-NN); download Sage binary | **submit via sbatch** (never login-node) |
| `mac` | `uname` = darwin | DIA-NN via **Docker** (no native mac build); Sage native | inline |
| `linux` | everything else | native binaries | inline |

`uc_davis_hive` is true when `/quobyte/proteomics-grp` exists → enables FASTA reuse
and the existing `.sif`.

Container runtime preference: hpc→apptainer, mac→docker, linux→native.

## DIA-NN by environment

- **No native macOS build exists.** On mac you must run DIA-NN through Docker. Set
  `DIANN_DOCKER_IMAGE` to a built image, or build one from the Academia Linux zip's
  bundled Dockerfile. `acquire_tools.sh` writes a note when this is unresolved.
- **HIVE (Proteomics Core):** DIA-NN is kept under `/quobyte/proteomics-grp/dia-nn/`.
  Recent versions are **native builds**, e.g. DIA-NN **2.6.0** at
  `build_260/diann-2.6.0/diann-linux`; older 2.3.0 is a `.sif`. `acquire_tools.sh`
  resolves the **pinned** version by looking for `build_*/diann-<version>/diann-linux`
  first, then a version-matched `.sif`, and **never silently substitutes a different
  version** (reproducibility). The facility's `run_diann_*.sbatch` in that folder is
  the reference invocation. AlphaDIA is also on HIVE at
  `/quobyte/proteomics-grp/apptainers/alphadia.sif` (auto-reused).
- **Linux native:** needs glibc ≥ Linux Mint 21.2 and .NET 8. If missing, prefer
  Docker/Apptainer.
- DIA-NN reads `.raw`/`.d` natively from 2.1+.

## Version pinning (reproducibility)

`acquire_tools.sh` honors `PIN_ENGINE`/`PIN_VERSION` from the workflow bundle and
caches under `~/.proteomics-pipeline/tools/<engine>/<version>/`. Different pinned
versions coexist. The written `tools.json` records `pinned` and `versions` so the
report can state exactly what ran. **Always pass the bundle's engine+version** — a
result from "latest" is not reproducible.

## FASTA resolution (`fetch_fasta.py`)

Priority, cheapest/most-trusted first:
1. `fasta.path` override in the workflow → used verbatim (pre-staged proteome).
2. **HIVE** (`--hive`): reuse `/quobyte/proteomics-grp/MRS/` proteomes +
   contaminants instead of downloading.
3. **UniProt**: stream the reference proteome for `fasta.uniprot_proteome`
   (e.g. `UP000005640` = human), then append cRAP universal contaminants if
   `add_contaminants` is set.

The script refuses to proceed if the resolved FASTA has 0 sequences, and warns
(does not fail) if contaminants can't be fetched. Common proteome IDs:
human `UP000005640`, mouse `UP000000589`, yeast `UP000002311`, E. coli `UP000000625`.

## SLURM submission (hpc)

`run_search.py --sbatch job.sh` emits a login-node-safe script:
`--partition=high --qos=genome-center-grp-high-qos`, 64G, 12h. If that queue is
full, the DE-LIMP fallback is `publicgrp/low` (the only fallback — there is no
`genome-center-grp` LOW partition). Submit with `sbatch job.sh`, poll the
`<job>_<id>.log`, then run `run_search.py --adapt-only` for Sage/FragPipe to build
`report.parquet`.
