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
   time and there is no moving branch to drift. The skill version (`.claude-plugin/
   plugin.json`) is recorded in `environment/skill.txt` and `run_manifest.json`'s `skill`;
   `workflow.manifest.json.registry.defaults_version` is only the date of the defaults
   table, not a version. Re-running the
   same skill version on the same data type reproduces the parameters exactly.
   (Before 2026-08-14 this came from a remote `workflows/` registry pinned by commit
   SHA. That registry is retired — old run records citing a SHA stay valid; see
   `workflows/README.md`.)
2. **Pinned engine version.** `acquire_tools.sh` honors `PIN_ENGINE`/`PIN_VERSION`
   from the manifest and records resolved commands + versions in `tools.json`. The
   recorded version is the build each command actually is — read from the HIVE build
   path, `.sif` name, download asset or Docker tag for DIA-NN, from the release
   tarball for Sage, from the unpacked folder or release asset for FragPipe, from
   `pip show` for AlphaDIA, and from an image tag or `.sif` name for Radiant — never the
   literal string `latest`, and never the **request** (`PIN_VERSION`), which would make
   the record echo the manifest back. What cannot be determined is recorded as `""`; a
   `sage` found on PATH is recorded as `env` (a source, not a build). A `tools.json`
   written before this change may still say `latest`.
   The manifest's pin and `tools.json` can disagree (the FRAN pilot: manifest 2.6.1,
   2.7.0 ran). `search_provenance.json` is the record of which build ran: `version` is
   `tools.json`'s, else the version the command itself names, else `null` — **never the
   manifest's pin**, which is kept beside it in the `engine_version` record as
   `manifest_pin`, with `mismatch` (`null` when the two cannot be compared) and a
   `source` saying where `version` came from (→ `references/environment.md`, "Version
   pinning"). `version` is only ever something shaped like a version number; any other
   string in `tools.json` is kept on the record and refused there. `reproduce.sh` still re-acquires the **manifest's** pin, so when
   `engine_version.mismatch` is `true`, edit its `PIN_VERSION` to `version` before
   replaying the search.
3. **Locked software environment.** `provenance.py` captures
   `environment/conda-explicit.txt` (every package pinned with URL + md5),
   `pip-freeze.txt`, and `r-sessionInfo.txt` (all R package versions). `reproduce.sh`
   rebuilds the env from the explicit lock — same packages, same versions. All of it is
   **the environment the DE ran in**: run_de.R records `runtime` in `de_provenance.json` (R
   home, Rscript, library paths, conda env, container), `r-sessionInfo.txt` is the DE's own
   `sessionInfo.txt`, and `run_manifest.json` `environment.source` says where the env came
   from. setup.json is used only for a DE record without `runtime`, and said to be (a DE once
   ran in a separate R 4.6 env while setup.json named the default one).
4. **Recorded inputs + parameters.** Copies of the exact params file and
   `conditions.csv`; the organism taxid, instrument, contrasts, and all thresholds
   in `run_manifest.json`; sha256 of the FASTA, the search report, and DE outputs.
   Raw files get a sha256 (or, for `.d` directories / >5 GB files, a structural
   fingerprint — name+size of every member) so input drift is detectable.
5. **A runnable recipe.** `reproduce.sh` re-creates the env, re-derives the search
   defaults from the data type (`resolve_defaults.py` — they ship with the skill, nothing
   is fetched), re-resolves the engine, rebuilds the FASTA, and re-runs search + DE
   with identical arguments. `REPRODUCE.md` is the human-readable version.

### The sequence database (`--fasta-info`)

The database is the one input that can silently differ between a run and its
"reproduction", so it is recorded explicitly rather than inferred. `fetch_fasta.py`
writes `<fasta>.meta.json` — sha256, source URL, organism + taxid, proteome ID,
database type (`content_used`), UniProt release — or, for a HIVE pre-staged copy, whose
release is unknown, `staged_file` (path, sha256, file date) plus `content_inferred` /
`content_check` — proteome vs contaminant sequence counts, contaminant set + citation, and
any build warnings. **Always pass it as
`provenance.py --fasta-info "$(cat search.fasta.meta.json)"`.**

`reproduce.sh` then rebuilds the database from *what actually ran*, not from the
workflow manifest's (`input/wf/`) default. This matters: the manifest holds the
**default**, but the user confirms the organism at step 3 and may have chosen a
different one — regenerating from the default would reproduce a different database and
quietly invalidate the comparison. Without `--fasta-info`, `reproduce.sh` falls back to
the manifest and labels that step as not recorded.

The rebuild replays the flags the original database was built with, from its sidecar
(`fetch_fasta.sidecar_state()` is the one reading of which rules built it), so a replay
reproduces THAT database rather than today's corrected one. `REPRODUCE.md`'s database note
says which flag was added and why; drop the flag to get the current database instead:
- **`--enzyme <list>`** — always passed: the digestion enzyme(s) recorded as
  `digestion_enzymes_used` (they decide which protease contaminant entries stay `Cont_`
  when they match a target). A sidecar from before the field existed replays the default,
  `trypsin,lysc`.
- **`--keep-target-contaminants`** — the database predates the removal of contaminant
  entries identical to a target protein (bovine ACTB = human ACTB …; sidecar state
  `legacy`), or was built with that flag. Without it the rebuild drops entries the original
  searched (153 human `Cont_` entries in the universal set).
- **`--min-unique-peptides 0`** — the database was built by the identity rule alone (state
  `identity_only`), before near-identical contaminants (e.g. bovine EEF1A1 vs mouse) were
  also removed. A recorded threshold other than the default is replayed as recorded.

Re-running later uses the *current* UniProt release, so entry counts may drift by
a few sequences. The recorded release and the FASTA sha256 in `checksums/` are
what make that drift visible instead of invisible. The same facts feed
`make_methods.py --fasta-meta` (the Methods "Sequence database" paragraph) and the
re-analysis `DIFFERENCES.md`, so an organism or database-type change shows up as a
difference rather than hiding behind an unchanged sequence count.

### What the orchestrator must do during the run

- **Log every command.** Append each command you execute (verbatim, full args) to
  `commands.log` and pass it via `--commands`. This is the audit trail.
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
├── reproduce.sh             # re-creates env, re-derives the shipped defaults, re-runs search+DE
├── MANIFEST.txt             # [OK]/[SKIPPED] capture log — read this to trust the bundle
├── environment/
│   ├── conda-explicit.txt   # fully pinned env lock (URL + md5 per package)
│   ├── pip-freeze.txt
│   ├── r-sessionInfo.txt    # R + limpa/limma/arrow/dplyr/tidyr versions
│   └── versions.txt         # tools.json versions + resolved commands. `sage --version`
│                             # is here as `sage_self_reported_version`, with the caveat
│                             # that the binary can disagree with its own release
├── inputs/                  # exact params file, conditions.csv, workflow manifest, commands.log
└── checksums/checksums.json # sha256 / fingerprints of raw, fasta, report, DE outputs
```

### Verifying a reproduction
After `reproduce.sh` runs, compare the new `de_results/` against
`checksums/checksums.json`. DE CSVs should match bit-for-bit when the env lock,
engine version, params, inputs **and the CPU family** all match. (DIA-NN/Sage are
deterministic for a fixed thread count + version; if you change thread count, intensities can
shift slightly — record threads in `commands.log`.)

### Exact DE numbers need the same kind of CPU (DPC-Quant)

limpa's `dpcQuant` fits each protein with `optim(method = "BFGS")` at its default tolerance
(`newton.polish` is off), and the objective runs through the BLAS library. OpenBLAS picks
kernels for the CPU it runs on (AVX-512 "SkylakeX" kernels on AMD zen4, "Haswell" kernels on
zen2), and they round the last bit differently; that moves where BFGS stops, inside its
tolerance. So identical inputs, software and flags give slightly different numbers on a
different CPU family — and identical numbers on the same one.

Measured on PROT_0756 v2 (2026-09-28): v2's DE ran on a zen4 node (EPYC 9734), the re-runs on
zen2 nodes (EPYC 7532/7662), same env (R 4.6.1, limpa 1.4.0, limma 3.68.5, OpenBLAS 0.3.34), same
input md5s and command. Across 12 contrasts **|ΔlogFC| ≤ 0.0022, |Δt| ≤ 0.0064, 0 significance
calls changed**. v2's command re-run on its own zen4 node reproduced its tables exactly, and a
zen2 re-run reproduced the zen2 pre-release tables exactly. Stage by stage: readDIANN and the
filters, and `dpcCN`, were bit-identical across families; `dpcQuant` differed (protein values by
up to 0.0026, standard errors by 7e-5), and was reproduced exactly within a family; limma
(`dpcDE`, `eBayes`) given the same protein values agreed to 2e-12.

What the skill does about it:
- `run_de.R` records the machine in `de_provenance.json` `compute` and at the top of
  `sessionInfo.txt`: the CPU model, the BLAS and LAPACK libraries, `OPENBLAS_CORETYPE` (as found,
  never set) and, on a SLURM node, the node and its CPU-family feature (`cpu_family`, from
  `scontrol show node`).
- `REPRODUCE.md` names that family and says how to get it again.
- **For an exact re-analysis on HIVE, submit to the same family**: `sbatch --constraint=<family>`.
  HIVE's CPU-family features (`sinfo -o '%f'`, 2026-09-28) are `zen`, `zen2`, `zen3`, `zen4`,
  `zen5` and `icelake`; in the `high` partition, 110 zen2, 24 zen4 and 16 zen3 nodes carry them.
  `sinfo -N -o '%N %f'` lists each node's.
- **Do not pin `OPENBLAS_CORETYPE`** to force one kind of kernel: on the EPYC 9734 node,
  `OPENBLAS_CORETYPE=Haswell` (8 threads or 1) and `=Zen` all segfaulted in `dpcQuant`
  (sbatch 24180764, 24180793, 24180794).

A difference at this level is not an error in either run, and it is far below what would
change a conclusion; but "reproduces exactly" means "on the same CPU family", and the record now
says which one that was. (limpa's `newton.polish = TRUE` would polish each fit to machine
precision; that is a method choice for limpa, not the skill.)
