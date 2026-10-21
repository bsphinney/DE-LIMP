# Reproducibility contract

Reproducibility here has **two levels**, and they answer different questions.

| | Artifact | Answers | Needs |
|---|---|---|---|
| **1. The analysis, as code** | `de_results/reproducibility_log.R` | "What did you actually do to the numbers?" | R + limpa/limma |
| **2. The whole run, pinned** | `reproducibility/` bundle | "Can I re-derive this from the raw files?" | conda env, engine build, hours |

**Level 1 is the one people read.** Level 2 is the one that makes the result
defensible. Produce both — but when you point a user at "the reproducibility", point
at the R script first.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 23: Level 1 — the analysis as plain R
- 51: Level 2 — the pinned bundle
- 226: Golden rule 5: every run reproducible, at two levels
- 246: For maintainers: adding a new `--method` or analysis path
- 270: 10. The bundle (step 10)
<!-- /contents -->

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

### The sequence database (`--fasta-info-file`)

The database is the one input that can silently differ between a run and its
"reproduction", so it is recorded explicitly rather than inferred. `fetch_fasta.py`
writes `<fasta>.meta.json` — sha256, source URL, organism + taxid, proteome ID,
database type (`content_used`), UniProt release — or, for a HIVE pre-staged copy, whose
release is unknown, `staged_file` (path, sha256, file date) plus `content_inferred` /
`content_check` — proteome vs contaminant sequence counts, contaminant set + citation, and
any build warnings. **Always pass it as a file:
`provenance.py --fasta-info-file search.fasta.meta.json`** (`--fasta-info @search.fasta.meta.json`
also works). Never `--fasta-info "$(cat search.fasta.meta.json)"`: a database with thousands of
added sequences has a sidecar of megabytes, and the command then fails with "Argument list too
long". A sidecar that is not valid JSON, or not in fetch_fasta.py's shape, stops provenance.py
(exit 2) before it writes anything. The bundle records each `--add-fasta` file by path, sha256
and entry count; the per-entry records stay in the sidecar.

`reproduce.sh` then rebuilds the database from *what actually ran*, not from the
workflow manifest's (`input/wf/`) default. This matters: the manifest holds the
**default**, but the user confirms the organism at step 3 and may have chosen a
different one — regenerating from the default would reproduce a different database and
quietly invalidate the comparison. Without `--fasta-info-file`, `reproduce.sh` falls back to
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

## Golden rule 5: every run reproducible, at two levels

5. **Every run must be completely reproducible**, at two levels. This is not optional.
   - **The analysis, as code.** `run_de.R` writes `reproducibility_log.R` — the whole
     DE as flat R with every value literal, runnable with just R + limpa/limma. This
     is the artifact to point a user at; it's what they mean by "the R code".
   - **The whole run, pinned.** As you go, append every command you run (verbatim,
     with all arguments) to a `commands.log`. At the end you MUST produce the
     reproducibility bundle (step 10) capturing the skill version (from
     `.claude-plugin/plugin.json`, written to `reproducibility/environment/skill.txt`),
     exact tool + package versions, all parameters, input and output checksums, and
     a runnable `reproduce.sh`. Parameters ship with the skill, so recording the
     skill version pins them — there is no external commit to chase.

   Never describe a result without both. (DE-LIMP architectural rules #1, #4.)
   "Exactly" means on the same CPU family: DPC-Quant's numbers shift in about the third
   decimal between zen2 and zen4 nodes (de_provenance.json `compute` records the machine). To
   match an earlier HIVE run exactly, submit with `sbatch --constraint=<its family>`; never pin
   `OPENBLAS_CORETYPE` (`references/reproducibility.md`).

## For maintainers: adding a new `--method` or analysis path

`tests/test_reproducibility_contract.py` enforces one rule: **every analysis path
the skill can run must emit a `reproducibility_log.R` that reproduces it.** A new
`--method` fails CI until it is listed in `COVERED` there *and* its emitter is
verified against a real run.

That guard exists because the same defect shipped twice. `repro_script.R` has a
separate emitter per method, and both wrote the scalar `--q-cutoff` where the run
applied per-column cutoffs — so the emitted script ran clean and returned
different numbers than the analysis beside it:

| method | run | emitted script |
|---|---|---|
| dpc | 4,648 proteins / 1,859 sig | 4,624 / 1,846 |
| maxlfq | 4,526 / 1,326 | 4,444 / 1,281 |

Nothing errored either time. Fixing dpc alone is what let maxlfq ship broken.

Before adding a path: run it, execute the script it emits, and diff the results —
`DELIMP_TEST_REPORT=<report.parquet> python3 tests/test_reproducibility_contract.py
TestEndToEndReproduction`. Reading the emitted script is not enough; both bugs
were invisible on inspection and obvious on execution.

## 10. The bundle (step 10)

Inside an analysis session, **`python3 scripts/provenance.py --session <session>`** fills every
argument below from the session's own records: the files the search ran on, the engine, its
version and `tools.json` (`search_provenance.json`), the DE method, contrasts and cutoffs and the
report the DE read (`de_provenance.json`), the workflow manifest, `input/acquisition.json`, the
FASTA and its sidecar (taxid), the commands log and `logs/qc_bracket.json`. An explicit argument
wins. A value with no record stays unset and is listed under `from_session.not_found` in the
printed JSON and as a MANIFEST.txt `[NOTE]` — relay those gaps; nothing is guessed. The full
form, for a run outside a session:

Assemble the bundle that makes the whole analysis reproducible (it hashes every input file; on
HIVE's login node it moves itself to a compute node — run it in the background):

```
python3 scripts/provenance.py --outdir ./reproducibility \
  --workflow-manifest ./wf/workflow.manifest.json \
  --setup-json ~/.proteomics-pipeline/setup.json \
  --tools-json ~/.proteomics-pipeline/tools/tools.json \
  --params ./wf/<params_file> --conditions ./conditions.csv \
  --fasta ./search.fasta --fasta-info-file ./search.fasta.meta.json \
  --raw /path/to/*.d --report ./search_out/report.parquet --de-dir ./de_results \
  --engine <engine> --de-method <dpc|maxlfq> --contrasts "<...>" \
  --q-cutoff 0.01 --logfc 1.0 --adjp 0.05 \
  --acquisition <DIA|DDA> --organism-taxid <taxid> --instrument "<name>" \
  --commands ./commands.log --timestamp "$(date -u +%FT%TZ)" \
  [--qc-bracket <session>/logs/qc_bracket.json]     # step 8e: path + sha256 only, no verdict
```

This writes `reproducibility/` with `run_manifest.json`, `reproduce.sh`,
`REPRODUCE.md`, the conda lock + pip freeze + R sessionInfo + tool versions, a
`skill.txt` recording **which skill produced this and how it was installed**, copies
of the params and conditions, and sha256 checksums of inputs/outputs. This step is
**mandatory and must always run** — code + versions are not optional. The environment it
captures is **the one the DE ran in**, from run_de.R's `runtime` record in `de_provenance.json`
(R home, Rscript, library paths, conda env, container) and the DE's own `sessionInfo.txt` — not
setup.json's, which names only the env setup.sh built. A `[NOTE]` in `MANIFEST.txt` says when
setup.json names another env (e.g. a DE run in a separate R 4.6 env): relay it. A DE record with
no `runtime` (before 2.10) falls back to setup.json / PATH and says so (`[SKIPPED] DE runtime`).
Still, **check the returned `skipped` count** and, if the conda
lock / R sessionInfo / checksums were skipped, fix the cause and re-run.
