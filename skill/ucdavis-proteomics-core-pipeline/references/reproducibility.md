# Reproducibility contract

Every run produces a `reproducibility/` bundle. A result without one is incomplete.
This is the skill's implementation of DE-LIMP architectural rules #1 (pipeline
self-describes) and #4 (no silent gaps — `MANIFEST.txt` logs `[OK]`/`[SKIPPED]`).

## The five things that make a run reproducible

1. **Pinned validated parameters.** `fetch_workflows.py` resolves the registry to a
   git **commit SHA** and `pull --ref <sha>` fetches params at that commit, not the
   moving `main` branch. The SHA is recorded in `workflow.manifest.json.registry`
   and the run manifest. A re-run pulls byte-identical `diann.cfg`/`sage_config.json`.
2. **Pinned engine version.** `acquire_tools.sh` honors `PIN_ENGINE`/`PIN_VERSION`
   from the bundle and records resolved commands + versions in `tools.json`.
3. **Locked software environment.** `provenance.py` captures
   `environment/conda-explicit.txt` (every package pinned with URL + md5),
   `pip-freeze.txt`, and `r-sessionInfo.txt` (all R package versions). `reproduce.sh`
   rebuilds the env from the explicit lock — same packages, same versions.
4. **Recorded inputs + parameters.** Copies of the exact params file and
   `conditions.csv`; the organism taxid, instrument, contrasts, and all thresholds
   in `run_manifest.json`; sha256 of the FASTA, the search report, and DE outputs.
   Raw files get a sha256 (or, for `.d` directories / >5 GB files, a structural
   fingerprint — name+size of every member) so input drift is detectable.
5. **A runnable recipe.** `reproduce.sh` re-creates the env, re-fetches the pinned
   workflow, re-resolves the engine, rebuilds the FASTA, and re-runs search + DE
   with identical arguments. `REPRODUCE.md` is the human-readable version.

## What the orchestrator must do during the run

- **Log every command.** Append each command you execute (verbatim, full args) to
  `commands.log` and pass it via `--commands`. This is the audit trail.
- **Pass the recorded commit SHA** to `pull --ref` and into `provenance.py` (it
  reads it from the workflow manifest).
- **Pass a timestamp** (`--timestamp "$(date -u +%FT%TZ)"`) — the scripts can't read
  the clock themselves.
- **Check the bundle's `skipped` count.** If the conda lock, checksums, or
  sessionInfo were skipped, fix the cause and re-run `provenance.py`. Don't hand
  over a bundle that silently dropped a critical artifact.

## Bundle layout
```
reproducibility/
├── run_manifest.json        # full machine-readable record (the master file)
├── REPRODUCE.md             # human-readable methods + how to re-run
├── reproduce.sh             # re-creates env, re-fetches @commit, re-runs search+DE
├── MANIFEST.txt             # [OK]/[SKIPPED] capture log — read this to trust the bundle
├── environment/
│   ├── conda-explicit.txt   # fully pinned env lock (URL + md5 per package)
│   ├── pip-freeze.txt
│   ├── r-sessionInfo.txt    # R + limpa/limma/arrow/dplyr/tidyr versions
│   └── versions.txt         # engine (DIA-NN/Sage) versions + resolved commands
├── inputs/                  # exact params file, conditions.csv, workflow manifest, commands.log
└── checksums/checksums.json # sha256 / fingerprints of raw, fasta, report, DE outputs
```

## Verifying a reproduction
After `reproduce.sh` runs, compare the new `de_results/` against
`checksums/checksums.json`. DE CSVs should match bit-for-bit when the env lock,
engine version, params, and inputs all match. (DIA-NN/Sage are deterministic for a
fixed thread count + version; if you change thread count, intensities can shift
slightly — record threads in `commands.log`.)
