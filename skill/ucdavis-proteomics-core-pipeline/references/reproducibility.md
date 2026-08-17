# Reproducibility contract

Reproducibility here has **two levels**, and they answer different questions.

| | Artifact | Answers | Needs |
|---|---|---|---|
| **1. The analysis, as code** | `de_results/reproducibility_log.R` | "What did you actually do to the numbers?" | R + limpa/limma |
| **2. The whole run, pinned** | `reproducibility/` bundle | "Can I re-derive this from the raw files?" | conda env, engine build, hours |

**Level 1 is the one people read.** Level 2 is the one that makes the result
defensible. Produce both — but when you point a user at "the reproducibility", point
at the R script first.

## Level 1 — the analysis as plain R

`run_de.R` writes `reproducibility_log.R` into the DE output directory: the whole
differential-expression analysis top to bottom, with every value written out
literally — the report path, the FDR cutoff and the q-columns it was applied to, the
sample→group map with real run names, any QuantUMS pre-filter, the covariates, the
design formula, the contrasts, the significance rule. No arguments to look up, no
config to cross-reference, no skill to install:

```
Rscript reproducibility_log.R
```

It is generated **from the objects that actually ran**, not re-derived from a
parameter file, so it cannot describe a different analysis than the CSVs beside it
(DE-LIMP architectural rule #1 — the pipeline self-describes). It writes to
`de_results_rerun/` so a re-run never clobbers the original.

The same lines go into `DE-LIMP_session.rds` as `repro_log`, so dropping that session
into the DE-LIMP GUI shows this code in its Reproducibility tab.

This is deliberately the same artifact DE-LIMP's GUI produces (`app.R`
`add_to_log()` → `R/server_session.R`), so a GUI run and a skill run hand the user
the same kind of thing.

**What it does not cover:** the search. It starts from `report.parquet`. Re-running
DIA-NN/Sage identically is level 2.

## Level 2 — the pinned bundle

Every run also produces a `reproducibility/` bundle. A result without one is
incomplete. This is the skill's implementation of DE-LIMP architectural rules #1 and
#4 (no silent gaps — `MANIFEST.txt` logs `[OK]`/`[SKIPPED]`).

### The five things that make a run reproducible

1. **Parameters pinned by the skill version.** `resolve_defaults.py` derives them
   from the data type and they ship *with* the skill, so nothing is fetched at run
   time and there is no moving branch to drift. `workflow.manifest.json.registry`
   records `defaults_version`; record the skill version alongside it. Re-running the
   same skill version on the same data type reproduces the parameters exactly.
   (Before 2026-08-14 this came from a remote `workflows/` registry pinned by commit
   SHA. That registry is retired — old run records citing a SHA stay valid; see
   `workflows/README.md`.)
2. **Pinned engine version.** `acquire_tools.sh` honors `PIN_ENGINE`/`PIN_VERSION`
   from the manifest and records resolved commands + versions in `tools.json`. For
   DIA-NN it resolves the build by asset filename, so the recorded version is the
   one actually installed — never the literal string `latest`.
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

### The sequence database (`--fasta-info`)

The database is the one input that can silently differ between a run and its
"reproduction", so it is recorded explicitly rather than inferred. `fetch_fasta.py`
writes `<fasta>.meta.json` — sha256, source URL, organism + taxid, proteome ID,
database type (`content_used`), UniProt release, proteome vs contaminant sequence
counts, contaminant set + citation, and any build warnings. **Always pass it as
`provenance.py --fasta-info "$(cat search.fasta.meta.json)"`.**

`reproduce.sh` then rebuilds the database from *what actually ran*, not from the
workflow bundle's `fasta.uniprot_proteome`. This matters: the bundle holds the
workflow **default**, but the user confirms the organism at step 3 and may have
chosen a different one — regenerating from the bundle would reproduce a different
database and quietly invalidate the comparison. Without `--fasta-info`,
`reproduce.sh` falls back to the bundle and labels that step as not recorded.

Re-running later uses the *current* UniProt release, so entry counts may drift by
a few sequences. The recorded release and the FASTA sha256 in `checksums/` are
what make that drift visible instead of invisible. The same facts feed
`make_methods.py --fasta-meta` (the Methods "Sequence database" paragraph) and the
re-analysis `DIFFERENCES.md`, so an organism or database-type change shows up as a
difference rather than hiding behind an unchanged sequence count.

### What the orchestrator must do during the run

- **Log every command.** Append each command you execute (verbatim, full args) to
  `commands.log` and pass it via `--commands`. This is the audit trail.
- **Pass the recorded commit SHA** to `pull --ref` and into `provenance.py` (it
  reads it from the workflow manifest).
- **Pass a timestamp** (`--timestamp "$(date -u +%FT%TZ)"`) — the scripts can't read
  the clock themselves.
- **Check the bundle's `skipped` count.** If the conda lock, checksums, or
  sessionInfo were skipped, fix the cause and re-run `provenance.py`. Don't hand
  over a bundle that silently dropped a critical artifact.

### Bundle layout
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

### Verifying a reproduction
After `reproduce.sh` runs, compare the new `de_results/` against
`checksums/checksums.json`. DE CSVs should match bit-for-bit when the env lock,
engine version, params, and inputs all match. (DIA-NN/Sage are deterministic for a
fixed thread count + version; if you change thread count, intensities can shift
slightly — record threads in `commands.log`.)
