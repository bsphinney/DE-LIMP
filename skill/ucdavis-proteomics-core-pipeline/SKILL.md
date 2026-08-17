---
name: ucdavis-proteomics-core-pipeline
description: >
  Run an end-to-end proteomics search + differential expression analysis from raw
  mass-spec data. Use this whenever the user wants to "analyze my proteomics data",
  "search these raw files", "run my DIA/DDA data", "find differentially expressed
  proteins", "process this timsTOF/Astral/Orbitrap run", or points at a folder of
  .raw / .d / .mzML files and asks what's in it. Detects acquisition + instrument,
  derives the search parameters from that data type, downloads the pinned search
  engine, runs DIA-NN (DIA) or Sage (DDA), then limpa/limma DE — with full
  provenance for every parameter. Also use it to "write the LC-MS methods
  section" / "generate a publication-ready methods section with the instrument grant
  acknowledgment" from facility raw data (UC Davis Proteomics Core).
---

# Proteomics Pipeline

Take a user from **raw MS files → differentially expressed proteins** using
parameters derived from their **data type** (instrument + acquisition), never
ad-hoc numbers. Those defaults ship with this skill — `scripts/estimate_params.py`
holds the instrument → mass-accuracy table, `scripts/resolve_defaults.py` picks the
engine and version, `scripts/make_presets.py` writes FragPipe/Radiant configs. So
**the skill version fully determines the search**, and nothing is fetched at run
time. Every number carries its source; surface it rather than asserting it.

> The old remote `workflows/` registry was retired 2026-08-14 — species never
> affected a search parameter, and a mutable registry meant one skill version could
> produce different parameters on different days. See `workflows/README.md`.

All scripts are in `scripts/` next to this file. Reference detail is in
`references/` — read the relevant one before the step it covers; keep this file as
the spine.

## Golden rules (do not violate)

1. **Confirm before committing compute.** A search is multi-hour. Show the resolved
   defaults (engine + version, mass accuracy + source) with the organism/design and
   get an explicit "go" before running the engine. **One confirmation, not a menu.**
2. **Never fabricate parameters.** If a value isn't a shipped default or given by
   the user, say so — don't invent an FDR, organism, or instrument. Every default
   carries its source (`ppm_source`); quote it rather than asserting the number.
   (DE-LIMP architectural rule #2.)
3. **Never run computationally intensive work on a cluster login/head node — EVER.**
   On any cluster (HIVE/SLURM, or any other scheduler), **every** heavy step — the
   search, and any large DE / figure / conversion (e.g. msconvert) — must run as a
   **scheduled job** (`run_search.py --sbatch` → `sbatch`), not inline on the login
   node. Login nodes are shared; running compute there gets the user flagged/killed.
   The **only** things allowed on the login node are tiny orchestration commands —
   submitting jobs, `squeue`/`sacct` polling (`watch_run.sh`), and small file moves.
   If you're unsure whether a step is heavy, submit it as a job. (No SLURM but a big
   dataset locally? Warn the user it may be slow rather than hammering a shared host.)
4. **Organism decides the FASTA and nothing else**, and is **always asked, never
   assumed** — not from a folder name, not "probably human". Resolve it with
   `fetch_fasta.py resolve` and have the user confirm the proteome. It does **not**
   affect a single search parameter: mass accuracy keys on instrument, and search
   behaviour on acquisition. Both are auto-detected and confirmed. So any organism
   is supported out of the box — there is no such thing as an unsupported species.
5. **Every run must be completely reproducible**, at two levels. This is not optional.
   - **The analysis, as code.** `run_de.R` writes `reproducibility_log.R` — the whole
     DE as flat R with every value literal, runnable with just R + limpa/limma. This
     is the artifact to point a user at; it's what they mean by "the R code".
   - **The whole run, pinned.** As you go, append every command you run (verbatim,
     with all arguments) to a `commands.log`. At the end you MUST produce the
     reproducibility bundle (step 10) capturing the skill's `defaults_version`,
     exact tool + package versions, all parameters, input and output checksums, and
     a runnable `reproduce.sh`. Parameters ship with the skill, so recording the
     skill version pins them — there is no external commit to chase.

   Never describe a result without both. (DE-LIMP architectural rules #1, #4.)
6. **Assume nothing about the user's environment — this skill runs for anyone.** Most
   users are **not** UC Davis Core and **not** on HIVE: they're on macOS, Windows
   (**WSL2** recommended, but **native Windows works too** — the engines ship Windows
   builds), or generic Linux. HIVE / `/quobyte/proteomics-grp` is a Core **fast-path,
   never a requirement**. For **every** program you have the user run, give a **public
   source** anyone can obtain (download URL / `pip`/install command) and a path that
   works on *their* OS — never hardcode an internal `/quobyte/...` path without its
   public fallback. → `references/search-engines.md` ("Public program sources") +
   `references/install.md`.

## Audience: assume nothing is installed

The user may be a biologist on a fresh laptop with no R, no Python packages, no
Docker. **Do not ask them to install things by hand.** Run `setup.sh` — it
installs everything that can be installed without admin rights into one
self-contained conda env. Only one thing ever needs the user's hands (Docker
Desktop, and only for DIA-NN on macOS), and `setup.sh`/`build_diann_docker.sh`
print the exact step. Never dump a generic "please install R/limma/..." message;
the scripts handle it.

## Flow

### 0a. Access & environment — ASK THIS FIRST
**Claude Code runs locally on the user's machine.** Heavy steps either run locally or
are driven on **HIVE over SSH with the user's private key** — they do NOT run Claude
Code on HIVE. Ask two questions, then verify and pick a mode:

1. **"Do you have access to UC Davis HIVE (an account + an SSH private key)?"**
   If yes, **ask where their private key is** (e.g. `~/.ssh/id_ed25519`).
2. **"Are you a member of the UC Davis Proteomics Core?"**

Verify (don't just trust the answers):
```
bash scripts/check_access.sh <hive_user> <private_key_path>
```
Read `recommended_mode` + `facility_software_available`, then:

- **HIVE = yes → `hive_remote`:** drive HIVE over SSH from the local Claude Code
  (`export HIVE_USER=… HIVE_KEY=…`; use `bash scripts/hive_exec.sh '<cmd>'`). The search
  runs as a **SLURM job** (`run_search.py --sbatch` → `hive_exec.sh 'sbatch job.sh'`),
  never the login node. HIVE gives **compute**; the Core software is separate (next).
- **Core member = yes (with HIVE) → reuse the installed software** in
  `/quobyte/proteomics-grp/`: `acquire_tools.sh` finds the DIA-NN `.sif`,
  `fetch_fasta.py --hive` reuses pre-staged FASTAs. No rebuilding.
- **HIVE = yes but NOT Core → rebuild the toolchain in your own HIVE home** (you can't
  read the group dir). Follow `references/access.md` → "Rebuild on HIVE" (run
  `setup.sh` + `acquire_tools.sh` + `fetch_fasta.py` on HIVE via `hive_exec.sh`). **List
  these steps for the user — do not guess.**
- **Core = yes but NO HIVE →** the Core software is on HIVE, unreachable without HIVE
  access → fall back to **local**; suggest requesting a HIVE account.
- **Both = no → `local`:** `setup.sh` installs the toolchain on the user's machine
  (no admin); public engines (DIA-NN Academia, Sage). Self-contained, works fully.

Tell the user which mode you're using and why. → full runbook: `references/access.md`.

### 0. One-time setup (install everything that's missing)
```
bash scripts/setup.sh
source ~/.proteomics-pipeline/activate.sh        # puts R, python, sage on PATH
```
`setup.sh` installs micromamba (no admin), then a conda env with R + limpa +
limma + arrow + Sage + Python/pyarrow, and writes `~/.proteomics-pipeline/setup.json`.
Sourcing `activate.sh` makes every later step use those interpreters — **source it
in the same shell before running anything below** (or prefix later commands).

Read `setup.json` and **gate on `ready_for`**:
- `ready_for.de` false → DE can't run; re-run `setup.sh` and report any `notes`.
- `ready_for.dia` false → on macOS this means Docker is missing. Run
  `bash scripts/build_diann_docker.sh` and relay its instructions (install Docker
  Desktop, open it once), then continue. Don't silently fall back to a DDA engine.
- `ready_for.dda` false → Sage/R not ready, or (macOS) msconvert is unavailable for
  `.d`/`.raw` → mzML; tell the user and, if their data is DIA, route to DIA-NN.

This step is idempotent — on a machine that's already set up it just verifies and
returns in seconds. → detail: `references/install.md`.

### 0b. Detect the environment
```
bash scripts/detect_env.sh > /tmp/env.json
```
Read `platform_class` (mac|hpc|linux), `container_runtime`, `uc_davis_hive`. This
decides how tools are acquired and whether to submit via SLURM.
→ detail: `references/environment.md`.

### 0c. Is this a RESUME? — CHECK THIS BEFORE STARTING ANYTHING
A cluster search runs for hours and the user **will** close their laptop. The SLURM
jobs survive that; the conversation does not. So before doing any work, check whether
an earlier session already submitted something:
```
python3 scripts/checkpoint.py find --base <where results live>   # e.g. ~ or the project dir
python3 scripts/checkpoint.py status --session <session dir>     # live sacct re-query
```
`status` returns each stage's real state, whether the expected output file exists, and
`next_commands` — the exact commands still outstanding. If a run is found:

- **Still RUNNING** → say so, give the job ids, and offer to wait (`watch_run.sh`) rather
  than resubmitting. **Never resubmit a search that is already running** — it wastes hours
  of cluster time and the user's budget.
- **COMPLETED and the output exists** → skip the search entirely and continue from
  `next_commands` (usually DE onward).
- **FAILED or stalled** → apply the `watch_run.sh` playbook (step 7b) and resubmit only the
  incomplete steps; finished `.quant` files and the predicted library are reused.

`submit.sh` writes `RECOVERY.md` + `.recovery.json` into the session directory
automatically, so this works even if the previous session ended mid-sentence. Tell the
user, in plain words, that they can close their terminal and that the jobs will keep
going — they do not need `tmux` or `screen` for this.

### 1. Locate the raw files
Ask the user for a directory or file list if not given. Recognized: `.d` (Bruker),
`.raw` (Thermo), `.mzML[.gz]`, `.wiff` (convert first). Glob to a concrete list.

### 1b. Check for a prior analysis of this dataset
```
python3 scripts/session.py find-prior --raw /path/to/*.d
```
If a match is returned, this is a **re-analysis** of an existing dataset — note the
prior session dir; you'll pass it to `session.py init --reanalysis-of` (step 3b) so
the run nests under `<prior>/reanalysis/` and gets a `DIFFERENCES.md`, and you'll
run the Comparator at the end (step 12). If no match, it's a fresh analysis.

### 2. Detect acquisition + instrument
```
python3 scripts/detect_acquisition.py /path/to/*.d /path/to/*.raw
```
Returns per-file `acquisition` (DIA/DDA/unknown) + `confidence`, plus an overall
`instrument`. **If `needs_confirmation` is true, ask the user** before continuing
— mixed/unknown/low-confidence must not silently pick an engine.
→ detail: `references/search-engines.md`.

### 3. Ask organism + experimental design (auto-map conditions)
- **Organism** cannot be detected — you MUST ask the user, every run. **Never** take
  it from a previous session, the folder name, or an assumption that it's human.
  Nothing in the skill supplies a default organism, by design.
  Getting this wrong silently searches the wrong species and every
  downstream number is void. Resolve their answer — don't guess the proteome ID:
```
python3 scripts/fetch_fasta.py resolve --organism "<what they said>"   # or --taxid <n>
```
  Accepts a common name, a scientific name, an NCBI taxid, or a `UP…` accession.
  (`scripts/resolve_organism.py --organism …` is a thin alias for the same resolver,
  kept for older call sites — there is only one implementation.)
  Show the user `selected` (organism, proteome ID, protein count) and **confirm it**.
  If `needs_menu` is true, present `candidates` as a menu instead of auto-picking —
  it fires when the top hit isn't a clear reference proteome or several compete
  (e.g. "mouse" also matches *Myotis myotis* and mouse-ear cress). Carry the
  confirmed `taxid` into step 4 and the `proteome_id` into step 6.
- **Database type** — ask what kind of protein database they want. Don't just
  take a one-word answer if they sound unsure: **most users don't know, and that's
  expected.** Explain the trade-off in their terms and recommend — every extra
  sequence enlarges the search space, which costs sensitivity through FDR
  correction, so bigger is not better:
  - `one_per_gene` — **the default; recommend this.** One canonical protein per
    gene (human ≈ 20.6k). Right for essentially all standard expression /
    differential-abundance work, which is what almost everyone is doing.
  - `reviewed` — Swiss-Prot only (manually curated). Use for organisms where the
    reference proteome is padded with unreviewed predictions.
  - `full` — every entry incl. unreviewed TrEMBL (human ≈ 147.5k, **7× bigger**).
    Only for poorly-annotated organisms where curation is thin.
  - `*_isoforms` variants — add splice isoforms. Only if the biological question
    is *about* isoforms. Isoforms share most peptides, so they mostly add
    ambiguous protein groups rather than answers.

  If they don't know, ask what they're trying to find out (standard "which
  proteins changed between my groups?" → `one_per_gene`), then say which you're
  using and why. Don't make them learn UniProt vocabulary to answer.
- **Contaminants** — ask which set to append (they are all `Cont_`-tagged, so
  DIA-NN's `--cont-quant-exclude Cont_` keeps them out of quant either way):
  - `universal` — **default**, and what the UC Davis Core stages on HIVE. Use it
    unless the user has a reason not to.
  - sample-type matched: `cell_culture`, `mouse_tissue`, `rat_tissue`,
    `neuron_culture`, `stem_cell_culture` — offer these if the sample type is known.
  - `none` — only if the user explicitly declines. Record it; the contaminant
    anomaly check in step 10 is meaningless without contaminants in the database.
- **Conditions:** ask the user to either *tell you* the conditions in plain words
  ("the first three are control, the rest treated") **or** *upload a file* (any
  CSV/TSV with a sample column and a group column, however named). Don't make them
  hand-fill a template.

First get the real run names the conditions must map to:
```
python3 scripts/collect_conditions.py --list-runs --from-dir /path/to/raw --glob '*.d'
```
Then map their conditions onto those runs — the script does the fuzzy filename
matching so you don't guess:
```
# (a) user uploaded a file:
python3 scripts/collect_conditions.py --map conditions.csv \
    --from-dir /path/to/raw --glob '*.d' --from-file <their_file>
# (b) user described it in words: turn it into intent, then ground it against the runs:
python3 scripts/collect_conditions.py --map conditions.csv \
    --from-dir /path/to/raw --glob '*.d' \
    --mapping-json '{"groups": {"control": ["..."], "treated": ["..."]}}'
```
Read the returned `ambiguities` and **confirm every one with the user** before
proceeding — `unassigned_runs` (a raw file no condition matched), `conflicting_runs`
(a file matched to two groups), `unmatched_identifiers` (the user named something
with no matching file), `multi_match_identifiers` (one label hit several files —
usually fine, e.g. a replicate prefix), and `singleton_groups` (<2 replicates → no
within-group variance). Do **not** start a search while any run is unassigned or
conflicting. Finalize the CSV, then validate it:
```
python3 scripts/collect_conditions.py --validate conditions.csv --against report.parquet
```
(Fallback: if the user has nothing yet, `--emit-template` writes a blank
File.Name,Group sheet for them to fill.) → detail: `references/conditions.md`.

### 3b. Create the analysis session (organize all files)
**Ask the user where the results should go** — two choices:
- **In the folder with their raw data** (default, recommended — keeps results next
  to the data): pass `--raw` and omit `--base`.
- **A central location** (e.g. their Documents folder, or a path they give): pass
  `--base <that path>` (results land under `<path>/sessions/`).

Then scaffold the session and route **everything** into it:
```
# default — results live with the raw data:
python3 scripts/session.py init --name "<short study name>" --raw /path/to/*.d \
    [--reanalysis-of <prior session dir from step 1b>]
# or central, if the user chose one:
python3 scripts/session.py init --name "<short study name>" --raw /path/to/*.d \
    --base ~/Documents/DataAnalysis
```
The output's `placement` tells you which was used. **Use the printed `paths` map for
every later step** — put
`conditions.csv` and the FASTA in `paths.input_dir`, search output in
`paths.search_out`, DE results in `paths.de_dir` (= `output/tables`), the
reproducibility bundle in `paths.repro_dir`, the report in
`paths.analysis_report`, and append commands to `paths.commands_log`. Raw files
are recorded (in `input/raw_files.txt`), never copied. With `--reanalysis-of`, the
session nests under `<prior>/reanalysis/<date>_<name>/`.
→ detail: `references/outputs.md`.

### 4. Resolve the defaults for this data type (then CONFIRM once)
```
python3 scripts/resolve_defaults.py \
    --acquisition DIA --instrument "timsTOF HT" \
    [--engine diann] [--env env.json] [--fasta db.fasta] [--threads 32] --dest ./wf
```
Writes `./wf/workflow.manifest.json` (engine + pinned version, mass accuracy and its
source, DE method) and, for FragPipe/Radiant, generates the config file itself. Pass
`--env` (step 1's `detect_env.sh` output) so an engine that **cannot run on this
machine** fails here with a usable alternative rather than mid-search.

**Confirm once, then run.** State the pick in one breath — data type, engine +
version, mass accuracy + where it came from, FASTA, DE method — and get a yes:

> timsTOF HT dia-PASEF → DIA-NN 2.6.1, MS1/MS2 15 ppm (DIA-NN README), mouse
> UP000000589 + contaminants, limpa DPC-Quant + limma. Run it?

Do **not** turn this into a menu. There is exactly one confirmation before compute,
and it is this one.

- **`--engine` is honoured, never overridden.** If the user names an engine that
  isn't the default for their data type, it runs anyway and `notes` says so. Never
  silently substitute a different engine.
- **`alternatives`** lists the other engines valid for this data type. Mention them
  in one clause; don't make the user adjudicate.
- **Overrides are explicit.** `--ms1-ppm/--ms2-ppm` force a site SOP value and are
  tagged as an override in the manifest, so a run record always distinguishes an SOP
  value from the shipped default. Don't override without a stated reason.
- Parameters come from the skill, so **record the skill version**, not a registry
  commit. `registry.defaults_version` in the manifest carries it.

#### 4a. DIA: the three routes

**DIA-NN is the default. Use it unless the user asks otherwise** — don't make them
adjudicate a tooling question to get their proteins. Name the alternatives in a
clause ("FragPipe and Radiant are also available") and move on.

| Route | Engine | Instruments | Runs on a Mac? | Notes |
|---|---|---|---|---|
| **DIA-NN** — *the default* | `diann` | Thermo + Bruker | Docker only (**emulated** on Apple Silicon) | Reads `.raw` and `.d` natively. |
| **FragPipe DIA** | `fragpipe` | Thermo *or* Bruker — **different presets** | **No — no macOS build at all** | Thermo → `DIA_SpecLib_Quant`. Bruker → `DIA_SpecLib_Quant_diaPASEF` **with diaTracer**. Not interchangeable. |
| **Radiant + Fulcrum** | `radiant` | **Thermo Orbitrap only** | **Yes — natively, incl. Apple Silicon** | Container reads mzML/Parquet, so `.raw` is converted first and Bruker `.d` is refused. Needs a DIA-NN-generated library. Licence-restricted — see below. |

**Platform gating (verified against each project's released artifacts).**
`detect_env.sh` emits an `engines` map; `resolve_defaults.py --env` enforces it.

- **FragPipe has no macOS build, Intel or ARM** — 24.0 ships only a Windows
  installer and a Linux zip. On any Mac this route is unavailable; point at HIVE.
- **DIA-NN has no macOS build either** and runs only in a linux/amd64 container. On
  **Apple Silicon that means Rosetta emulation** — correct results, but far slower.
  Say so before a big cohort, and offer HIVE.
- **Radiant's image is multi-arch**, so on an M-series Mac it is the **only DIA
  engine that runs natively**. Worth naming when someone is working locally.
- **Sage** ships a native `aarch64-apple-darwin` build, so DDA on a Mac is fine.
- **diaTracer is Bruker-only.** It converts dia-PASEF `.d` into pseudo-MS/MS
  spectra; it has no role on Thermo data and is switched off in the Thermo preset
  (`diatracer.run-diatracer=false`). The two FragPipe presets also differ in MBR,
  missed cleavages, precursor charge, mass tolerance, and topN peaks — never swap
  one for the other by flipping the diaTracer flag.
- **Presets are generated, not stored per organism.** `make_presets.py` picks the
  vendor template by data type and patches only what's run-specific (FASTA, threads)
  — see step 6c. If a route can't apply to the data (Radiant on Bruker `.d`), it
  refuses with a pointer instead of producing a config no vendor validated.
- **Radiant licence:** Apache-2.0 **+ Commons Clause + mandatory grant-back**. The
  Commons Clause restricts *selling* a service whose value derives substantially
  from the software. **Say this before running it for anything fee-for-service** —
  it surfaces in `tools.json` `notes` too. Don't advise on the legal question;
  point at it and let them decide.

#### 4b. Offer to run all three and compare

After they pick, **offer the bake-off**: run two or three routes on the same files
and compare which performed better. Say roughly what it costs (each route is a full
search; Radiant is **serial** across files) so it's an informed choice.

If they accept: run each route into its **own session**, then compare with
`compare_searches.py` + `make_comparison_report.py` (step 11c). **Hold everything
except the search constant** — same FASTA, same organism, same conditions, and push
every engine's report through `run_de.R` with identical `--method`, contrasts, and
thresholds, so the search is the only variable. Name the version gap out loud:
FragPipe 24 bundles DIA-NN 1.8.2 beta 8, two majors behind what `diann_*` pins.
→ detail: `references/cross-tool-comparison.md`.

**Every organism works.** There is no "unsupported species" case any more —
parameters key on instrument and acquisition, and the FASTA comes from step 3. Never
tell a user their organism isn't covered.

⚠ **The organism always comes from step 3**, never from anything cached or inferred.
Searching rat data against a human proteome silently voids the whole run.

### 5. Acquire the pinned engine
Honor the manifest's exact version — not "latest":
```
PIN_ENGINE=diann PIN_VERSION=2.6.1 bash scripts/acquire_tools.sh <platform_class>
```
Reads/writes `~/.proteomics-pipeline/tools/tools.json`. On HIVE it reuses the
existing `.sif`; on mac it uses Docker for DIA-NN. **Read `tools.json` `notes`** —
license gates (FragPipe, Radiant) and missing-runtime warnings surface there.

Radiant is a ~3 GB image pull, so it is **not** acquired unless asked for:
```
ACQUIRE_RADIANT=1 PIN_ENGINE=radiant PIN_VERSION=2.3.3 bash scripts/acquire_tools.sh <platform_class>
```
It needs Docker or Apptainer — Seer ships no native binary. On HIVE the `.sif` must
already exist; building one pulls ~3 GB and **must run on a compute node, never the
login node** (the note in `tools.json` gives the exact `srun apptainer build` line).

### 6. Build the FASTA
Use the proteome, database type, and contaminant set the **user confirmed in
step 3**:
```
python3 scripts/fetch_fasta.py fetch --proteome <confirmed UPID> \
    --content <one_per_gene|reviewed|reviewed_isoforms|full|full_isoforms> \
    --contaminants <universal|cell_culture|...|none> \
    --out ./search.fasta [--hive]
```
Defaults are `--content one_per_gene` (the canonical set, from UniProt's
reference-proteome FTP tree) and `--contaminants universal`. **Do not pass
`--content full` casually** — for human that is 147,506 sequences instead of
20,652, a 7× larger search space. UniProt's REST `onePerGene` parameter is
silently ignored, which is why the canonical set must come from FTP.

Pass `--hive` when `uc_davis_hive` is true to reuse pre-staged proteomes +
contaminants under `/quobyte/proteomics-grp/MRS/` instead of downloading.

**Non-model organisms → NCBI.** UniProt reference proteomes cover a few thousand
species; wildlife, agricultural and other non-model organisms often have an
annotated RefSeq genome and no UniProt proteome. `resolve` searches NCBI
automatically when UniProt returns nothing and lists hits under `ncbi_candidates`:
```
python3 scripts/fetch_fasta.py fetch --ncbi-accession GCF_007827085.1 \
    --ncbi-organism "Peromyscus californicus insignis" --ncbi-taxid 42520 \
    --contaminants universal --out ./search.fasta
```
- **Confirm the assembly with the user — never auto-pick one.** It's a different
  kind of database and they need to know that's what they're searching.
- **Use the RefSeq `GCF_` accession.** A GenBank `GCA_` record usually reports no
  protein count and its download ZIP contains no `protein.faa` at all.
- **It includes every isoform** — the analogue of UniProt `full_isoforms`, *not*
  the `one_per_gene` default. The Core's Peromyscus assembly is 51,825 sequences
  against 21,877 protein-coding genes. `fetch` counts them and warns; carry that
  into the methods rather than calling it a "reference proteome".

**Read the returned JSON** (also written to `<out>.meta.json` for the repro
bundle) and act on it:
- `warnings` non-empty → tell the user before searching; a one-per-gene→full
  fallback changes the database out from under them.
- `diann_cont_quant_exclude` → pass as `--cont-quant-exclude Cont_` to DIA-NN in
  step 7 so contaminants are identified but excluded from quant + normalisation.
- Contaminants that can't be fetched are a **hard stop**, not a warning. Fix the
  source or have the user explicitly choose `--contaminants none`.
→ detail: `references/environment.md` ("FASTA").

### 6b. Estimate search parameters from the data type
Search parameters are **derived from what the data is**, not hand-maintained.
For DIA-NN and Sage, generate them (for FragPipe/Radiant see 6c — step 4 already
wrote their config):
```
python3 scripts/estimate_params.py --engine <diann|sage> \
    --acquisition <DIA|DDA> --instrument "<detected instrument>" \
    --precursor-mz-range <LO> <HI> \
    [--var-mods ox] [--overrides '<site SOP values as JSON>'] \
    --fasta-meta ./search.fasta.meta.json \
    --out ./wf/params.<cfg|json>
```
**Always pass `--precursor-mz-range`**, taking `precursor_mz_range` straight from
step 2's `detect_acquisition.py` output. Without it the range falls back to
380–980 — which on a timsTOF method acquiring 299.5–1200.5 silently discards both
tails, does not error, and is invisible in the results. The rationale tags the
value `measured …` or `FALLBACK …` so you can tell which you got; if it says
FALLBACK, say so to the user before committing to a multi-hour search.

Always pass `--fasta-meta` (step 6's sidecar): it carries the contaminant tag, so
the cfg gets `--cont-quant-exclude Cont_` and contaminants are identified but kept
out of quantification and normalisation. Both the single and 5-step parallel
DIA-NN paths read this cfg, so the flag applies to library generation *and*
analysis. With `--contaminants none` the flag is correctly absent.
The estimator keys mass tolerances on the instrument class from DIA-NN's
known-good table (Astral → MS1 4/MS2 10 ppm; timsTOF → 15 ppm; unidentified →
automatic calibration), sets DIA/DDA window mode from acquisition, and uses
standard trypsin/LFQ defaults for the rest. **It prints a `rationale` tagging
every value's provenance** — surface this to the user (and it flows into the
methods text), so a derived default is never mistaken for a confirmed setting.
Use the resulting file as `--params` below. → detail: `references/parameters.md`.

### 6c. FragPipe / Radiant: the config is generated, not stored
These two are driven by a whole config file rather than flags, so
`resolve_defaults.py` (step 4) already generated one into `./wf/` by calling:
```
python3 scripts/make_presets.py --engine <fragpipe|radiant> \
    --instrument "<detected instrument>" --acquisition DIA \
    [--fasta ./search.fasta] [--threads N] --out ./wf/<engine>.<workflow|radiantConfig>
```
Re-run it directly only if the FASTA or thread count changed after step 4.

It starts from the **vendor's own preset for that data type** and patches only what
is run-specific. That matters: FragPipe's two DIA presets differ in **20 keys**, not
one — MBR, missed cleavages, precursor charge and tolerance, topN peaks, and the
var-mod table — so flipping `diatracer.run-diatracer` on the wrong template yields a
config no vendor ever validated. Never hand-edit a preset to switch routes; pick the
right template by data type.

- **timsTOF** → `DIA_SpecLib_Quant_diaPASEF.workflow`, diaTracer **on**.
- **Thermo** → `DIA_SpecLib_Quant.workflow`, diaTracer **off** (MSFragger-DIA
  searches the spectra directly; there are no pseudo-spectra to trace).
- **Radiant on Bruker `.d` is refused**, not adapted — the container reads only
  mzML/Parquet. It also narrows Seer's flat 20/20 ppm default using the instrument,
  tagged `[DERIVED]` because that inference is ours, not Seer's.
- Mass tolerances otherwise keep the **vendor's** values: FragPipe already tunes them
  per data type, and DIA-NN's ppm table describes DIA-NN's matcher. Override only
  with a stated SOP value (`--ms1-ppm/--ms2-ppm`).

`<out>.provenance.json` records every key changed and why — it feeds the methods text.

### 7. Run the search
```
python3 scripts/run_search.py --tools ~/.proteomics-pipeline/tools/tools.json \
    --bundle ./wf/workflow.manifest.json --params ./wf/params.<cfg|json> \
    --fasta ./search.fasta --out ./search_out --files /path/to/*.d --threads 16
```
- DIA → whichever of the three routes the user chose in step 4a (**default DIA-NN**);
  DDA → Sage; FragPipe/Radiant when the bundle names them or the user asks.
- **Radiant + Fulcrum (Thermo Orbitrap DIA).** Add `--engine radiant`:
  ```
  python3 scripts/run_search.py --tools ~/.proteomics-pipeline/tools/tools.json \
      --bundle ./wf/workflow.manifest.json --params ./wf/default.radiantConfig \
      --fasta ./search.fasta --out ./search_out --files /path/to/*.raw \
      --engine radiant --threads 16
  ```
  - **Thermo only.** Bruker `.d` is refused with an error pointing at the DIA-NN or
    diaTracer route. `.raw` is converted to mzML first via the same `msconvert` path
    Sage uses — on a Mac that means no local conversion, so run it on HIVE/Linux
    (→ `references/install.md`).
  - **A spectral library is always required**, and `run_search.py` generates one with
    **DIA-NN's predictor** (`--fasta-search --predictor --gen-spec-lib --out-lib
    …tsv`) because Radiant reads DIA-NN's library TSV schema natively. Pass
    `--library` to reuse an existing DIA-NN `.tsv` library instead. This means the
    Radiant route **also needs DIA-NN acquired**.
  - ⚠ **`--libfree` is a misnomer.** In Radiant's CLI it is the same switch as
    `--no-mbr`, so it selects single-pass vs match-between-runs — it does **not**
    mean "no library". **MBR is ON by default in this skill** (`--no-mbr` to disable),
    because the DIA-NN chain builds an empirical library from all runs and re-searches
    against it; leaving Radiant single-pass makes it the only engine without cross-run
    information and understates it. Seer's shipped example TOML has `mbr = false`.
  - **On a cluster, this parallelises automatically.** Radiant's own backend is
    serial (`NotImplementedError` for parallel mode), so N files would cost N × one
    file. But the *search* is per-file independent — only rescoring/FDR/inference/
    rollup need the whole set. With >1 file and `sbatch` on PATH, `run_search.py`
    emits a **3-step chain** via `radiant_parallel.py`:
    1. DIA-NN predicts the library (1 job)
    2. **one `RadiantDIA` per file as a SLURM array** (N tasks) → `.radiantDIA` results
    3. **one Fulcrum job** rescores the whole set (`reuse_existing = true` makes it
       skip re-searching and go straight to FDR + inference + rollup)

    Submit with dependencies (the JSON output prints the exact lines), then
    `run_search.py --adapt-only --engine radiant` to build `report.parquet`. Step 2 is
    idempotent, so a preempted `publicgrp/low` task requeues without redoing work.
    Use `--no-parallel` to force the single serial search.
  - Licence: Apache-2.0 + **Commons Clause** + mandatory grant-back — see step 4a.
- **FragPipe + diaTracer (timsTOF dia-PASEF, spectrum-centric DIA).** A second DIA route,
  used when the bundle names it or the user asks for it. diaTracer traces 3-D features and
  writes pseudo-MS/MS spectra, MSFragger searches those to build a library from the data
  itself, then **DIA-NN quantifies against it**. Worth reaching for as a cross-check on the
  DIA-NN route, or when the search needs semi-tryptic / nonspecific / open-modification
  behaviour that library-free DIA-NN handles poorly.
  - **⚠ Licensing — get this right, it is easy to state wrongly.** The gate is
    **MSFragger / IonQuant / diaTracer** (one shared academic licence; commercial users
    license via fragmatics.com). The DIA-NN that FragPipe bundles is explicitly *"available
    for both academic and commercial use within FragPipe"*, so DIA-NN is **not** the
    restriction here. A commercial user is not locked out — they must license the three
    gated tools. **Never tell them the route is closed to them.**
  - **⚠ The clause that matters for a core facility:** the academic licence forbids *"any
    work product, report or deliverable by LICENSEE, for a commercial or for-profit
    organization"* — which reads on fee-for-service work for an industry client. If asked
    whether their work qualifies, send them to info@fragmatics.com and their licensing
    office. Do **not** rule on it yourself. → `references/search-engines.md`.
  - Output is DIA-NN's own, in `<workdir>/dia-quant-output/`. With FragPipe's bundled
    DIA-NN 1.8.2 Beta 8 that is **`report.tsv` — there is no `report.parquet`** unless
    someone configured DIA-NN 2.x. `adapt_fragpipe` handles either.
  - **Never point FragPipe at shared raw files directly — `run_search.py` stages symlinks
    for you.** diaTracer writes its mzML *next to the input*, so two people searching the
    same dataset (a class working from one folder) would race to write the same file, and a
    read-only share fails outright. `diatracer_stage.py` builds a per-run directory of
    symlinks; FragPipe normalizes paths without resolving them, so the output lands in
    **your** directory and the shared raw folder is never written to. This happens
    automatically for DIA + `.d` inputs; nothing to pass.
  - **Existing pseudo-spectra are reused, never regenerated.** The stager looks for the
    mzML in the staging dir *and* beside the real `.d`, and when it finds one emits the
    reuse form — the mzML as `DDA` plus the `.d` as `DIA-Quant` — so DIA-NN still
    quantifies against the real chromatograms. That saves ~20 min per file, so a re-search
    with different parameters is cheap. Report how many were reused.
  - A manifest data type that isn't one of `DIA`/`GPF-DIA`/`DIA-Quant`/`DIA-Lib`/`DDA+`/
    `DDA` is **silently treated as DDA**, quietly turning a DIA run into a DDA one.
  - FragPipe 24.0 bundles **DIA-NN 1.8.2_beta_8**, older than the 2.x the `diann_*`
    workflows pin. Expect different numbers from that alone; say so rather than
    attributing every difference to the spectrum-centric approach.
  - **Four things block this route on a fresh account — read `references/fragpipe-diatracer.md`
    BEFORE running it.** (1) the three engines are licence-gated and the tools folder must be the
    one that also contains `speclib/`; (2) MSFragger needs a **target+decoy** FASTA, which DIA-NN
    does not; (3) headless FragPipe rejects the spectral-library module until
    `fragpipe-config.bin-python` exists in `~/.config/FragPipe/.../fragpipe-ui.cache` — no CLI flag
    substitutes, run `scripts/fragpipe_bootstrap.py`; (4) diaTracer's standalone `main()` ignores its
    own documented defaults and needs all seven numeric options passed.
  - **On a cluster, parallelise the conversion automatically** with `scripts/diatracer_parallel.py`
    — one SLURM task per file instead of FragPipe's serial loop. Measured 9 files in ~25 min vs ~3 h.
    Re-stage afterwards so the manifest takes the reuse form and FragPipe skips conversion.
  - Budget roughly **20 min per file for the diaTracer conversion alone** (paper: 34 files,
    32 cores, ~19.4 min each), before MSFragger and DIA-NN. The mzML are **reusable**, so a
    re-search with different parameters skips it. Needs Java 11+, MSFragger 4.4+,
    IonQuant 1.11.18+, diaTracer 2.2.1+. The three gated jars are **not** in the FragPipe
    zip and **cannot be scripted** — a human accepts the licence once and mints a key, then
    point `FRAGPIPE_TOOLS_FOLDER` at them. → `references/search-engines.md`.
- **⚠ DIA-NN licensing (ask before a DIA run):** DIA-NN's free "Academia" build is
  licensed for **academic / non-profit use only**. If the user is in a **commercial /
  non-academic** setting, do NOT use DIA-NN — **offer AlphaDIA** (`--engine alphadia`,
  Apache-2.0, commercial-OK; the open-source DIA alternative) or Sage (MIT, also
  commercial-OK). A quick question: *"Is this academic/non-profit, or commercial?"* —
  commercial → AlphaDIA. (AlphaDIA is deep-learning based; a GPU is strongly
  recommended — best on a HIVE GPU node or a GPU workstation.)
- **⚠ AlphaDIA limitation — timsTOF whole-proteome directDIA:** AlphaDIA does **not**
  work on **timsTOF** `.d` files for **whole-proteome library-free (directDIA)**
  searches. So for timsTOF + whole-proteome + library-free, AlphaDIA is NOT an option.
  In that case: academic → DIA-NN; commercial → use a **spectral-library** approach
  (predicted/empirical library) rather than directDIA, or Sage — and tell the user the
  constraint. (Non-timsTOF instruments, or library-based timsTOF runs, are fine for AlphaDIA.)
- **DIA-NN search settings:** always use the **estimated cfg** (step 6b) — it encodes
  DIA-NN's official recommended per-instrument mass tolerances (timsTOF → MS1/MS2 15
  ppm; Orbitrap Astral → 4/10 ppm; Orbitrap by resolution 240k→4, 120k→7, 60k→10,
  30k→15). Don't override these unless the user gives a validated SOP value.
- **DIA-NN parallel — AUTOMATIC above 5 files.** `run_search.py` routes to the **5-step
  parallel chain** by itself whenever the run is DIA-NN with **more than 5 files** on a
  machine that has SLURM and a `--cfg` with **fixed** mass accuracy. You don't decide
  this; it prints `parallel routing: YES|no -- <reason>` and records the reason in
  `search_provenance.json`. This matches facility practice in DE-LIMP: per-file passes
  run as a SLURM array instead of one long single-node job.
  - It falls back to a single search, with the reason printed, when there's **no SLURM**
    (the chain is job arrays — nothing to fall back to) or when the cfg is **not
    parallel-safe**. Steps 3/5 reuse the `.quant` files, so *anything* DIA-NN
    auto-optimises per file gets stitched together inconsistently — DIA-NN's own warning
    names **mass accuracy AND scan window**. Both must be pinned:
    - **mass accuracy** — re-run `estimate_params.py` with the **real instrument** (step 6b).
    - **`--window`** — `estimate_params.py` omits it, which means auto. **Measure it**,
      never guess: on an 18-file poplar run DIA-NN inferred radius 7 for seventeen files
      and 8 for one, and the chain combined them. After step 1 has built the library:
      ```bash
      python3 scripts/probe_window.py --diann "<diann cmd>" --raw <one file> \
          --fasta <f.fasta> --lib <step1.predicted.speclib> --write-cfg <params.cfg>
      ```
      It exits as soon as DIA-NN logs `Scan window radius set to N` (during calibration),
      so it costs minutes. One file covers a set acquired with the same method.
  - Override either way with `--no-parallel` (force one job) or `--parallel-threshold N`.
  - It **generates** the chain but does not submit it. Submit `<out>/submit.sh` (over
    `hive_exec.sh` on HIVE), then watch the **step-5** job — that's the one that writes
    `report.parquet`. → `references/diann_parallel.md`.
- **On `hpc`:** add `--sbatch job.sh`, then `sbatch job.sh` (over `hive_exec.sh` for a
  remote HIVE run). Re-run with `--adapt-only` afterward for Sage/FragPipe/AlphaDIA to
  build `report.parquet`.
- Output is normalized to the **DE contract**: a DIA-NN-shaped `report.parquet`.
→ detail: `references/search-engines.md`.

### 7b. Watch the run — MANDATORY, auto-correct errors
**The moment a search starts, watch it — never leave a multi-hour search
unmonitored.** Poll with `watch_run.sh` in a loop until it finishes:
```
# local: python3 scripts/run_search.py ... run_in_background, then loop:
bash scripts/watch_run.sh --log <search log> --out <search out dir> --poll <n>
# SLURM on HIVE:
bash scripts/watch_run.sh --slurm <jobid> --log <job log> --out <out> --poll <n> --hive
# multi-step chains (diann_parallel / radiant_parallel): watch EVERY step, never just the last
bash scripts/watch_run.sh --all <search out dir> --hive
```
**Always pass `--out` and increment `--poll` each time.** `--out` is what turns "still
running" into real progress: the chain drops a `.quant` per finished file, so the watcher
reports the actual step and file count rather than guessing. `--poll` rotates the note.

**Report progress to the user on every poll — a search is hours long and silence reads as
a hang.** Give them one short line from `progress.summary` plus `doing`, and pass along
`note` (with `note_source`) so there's something to read while waiting. For example:

> Step 2 of 5 — searching each file against the predicted library. 34 of 66 files done
> (52%). *While you wait: decoy sequences exist to be wrong on purpose. How often a search
> prefers a decoy is what calibrates the false discovery rate. (Elias & Gygi 2007)*

Keep it to a couple of lines; don't dump the JSON. Never present the note as a finding
about **their** data — it is general background, and it must not be mistaken for a result.
It returns `{state, done, failed, stalled, error_class, fix}`. While `done` is false,
keep polling (sleep ~60s between polls; for long SLURM jobs use a scheduled wake-up).
Pass `--log` so it can also catch **stalls** — a job stuck in `RUNNING` with its log
frozen (a hung file; `sacct` will show RUNNING forever). **If `failed` OR `stalled` is
true, apply the `fix` and resubmit AUTONOMOUSLY — do not wait for the user.** Starting
the search needed confirmation (golden rule #1); recovering a run that's already in
flight does not. Common auto-fixes (→ `references/watcher.md` playbook):
- **stalled/hung file:** `scancel` that task, retry once on a fresh node; if it stalls
  again, **drop that one file** and continue (step 4 of the parallel chain auto-skips a
  file with no `.quant`); note the dropped file in **Data Quality Notes**.
- **failed step in the 5-step chain:** resubmit from the earliest incomplete step,
  **reusing completed `.quant` and `step1.predicted.speclib`** — never restart the whole
  cohort; broken `afterok` after a killed task → resubmit steps 3→4→5 fresh.
- OOM → raise `--mem`; timeout → raise `--time`; missing temp dir → `mkdir -p` first.

**Every search — single-shot or parallel array — must run under this watch loop.** Stop
and tell the user only after 2 failed auto-fixes of the same class (dropping 1 pathological
file out of many is success, not a failure). Report what you recovered. Only proceed to DE
once `COMPLETED` and `report.parquet` exists. → detail: `references/watcher.md`.

### 8. Differential expression
```
Rscript scripts/run_de.R --input ./search_out/report.parquet \
    --metadata conditions.csv --method <dpc|maxlfq> --outdir ./de_results
```
**`dpc` (limpa) is THE DEFAULT — use it unless the user asks otherwise or the data
cannot support it.** limpa models the detection-probability curve and quantifies from
precursor intensities directly, so it uses the whole measurement rather than a
pre-collapsed protein number. Don't switch away from it because a bundle happens to
say `maxlfq`; switch only for the two reasons below.

Use `maxlfq` when **either**:
1. **The user asks for it** — e.g. they want QuantUMS quality filtering (below), or
   to match a previous MaxLFQ-based analysis.
2. **The input cannot support limpa.** `readDIANN()` needs PRECURSOR-level columns
   (`Precursor.Id`, `Precursor.Normalised`). What matters is *which file you point at*,
   not which engine ran:

   | engine | file to use for `dpc` | verified |
   |---|---|---|
   | DIA-NN | `search_out/report.parquet` (native) | yes |
   | FragPipe | `dia-quant-output/report.tsv` + `--format tsv` | yes — its DIA route bundles DIA-NN, so this **is** a DIA-NN report, `Lib.Q.Value`/`Lib.PG.Q.Value` included |
   | Radiant | `radiant_to_delimp.py --out <x>.parquet` | yes |

   What does **not** work is the *adapted* `report.parquet` from `adapt_sage` /
   `adapt_fragpipe` / `adapt_radiant` — those deliberately collapse to one row per
   protein × run to feed the maxlfq path. `run_de.R` checks up front and names the
   precursor-level file to use instead.

Writes `DE_<method>_<contrast>.csv` + `Expression_Matrix.csv` +
`methods.txt` + `sessionInfo.txt` + `de_provenance.json` (exact R package versions) +
**`reproducibility_log.R`**.

`reproducibility_log.R` is the whole analysis as plain, flat R — every value written
out literally (report path, FDR cutoff, sample→group map, design, contrasts), runnable
with `Rscript` using nothing but R and limpa/limma. It is generated from the objects
that actually ran, so it can't drift from the CSVs next to it. **This is what to show
a user who asks "what did you do?" or "can I have the R code?"** — not the
`reproducibility/` bundle, which pins the environment and re-runs the multi-hour
search. Mention it when you hand over results.

**QuantUMS quality filtering — opt-in, `--method maxlfq`, DIA-NN input only.**
Only reach for this when the user asks. DIA-NN's
QuantUMS scores how well MS1 and MS2 agree for each measurement, and filtering on it
trades depth for quantitative reliability:

| flag | column | default | what it does |
|---|---|---|---|
| `--eq-cutoff` | `Empirical.Quality` | 0 (off) | drop precursors whose empirical quality is below this |
| `--pgq-cutoff` | `PG.MaxLFQ.Quality` | 0 (off) | drop protein groups whose MaxLFQ quality is below this |
| `--coverage-min` | — | 0.5 | require a protein in ≥ this fraction of samples before testing |

Both QuantUMS cutoffs are **off by default** — they discard real measurements, so
they are the user's call, not a silent default. Offer them when quantitative
precision matters more than depth (small fold-changes, few replicates), say roughly
how many proteins each costs, and record the choice: the counts land in
`filters_applied` and flow into `methods.txt`. They need DIA-NN input — the columns
don't exist in Sage/FragPipe/Radiant output, and the filter is skipped if absent.

`--coverage-min` is different: it is **on** at 0.5 because a protein seen in one or
two samples is an on/off observation, and letting such rows into `eBayes` destabilises
the variance moderation for *every* protein, not just those rows.

→ detail: `references/de-analysis.md`.

### 8b. Generate figures
```
Rscript scripts/make_figures.R --de-dir ./de_results --conditions ./conditions.csv \
    --outdir ./figures --adjp 0.05 --logfc 1
```
Produces publication-quality volcano (per contrast), PCA, heatmap of top proteins,
p-value distributions, and a per-sample protein-count QC plot, plus `figures.json`
(captions). These get embedded in the report.

**Significance is `adj.P.Val` alone — there is no fold-change cutoff anywhere in this
pipeline.** `--logfc` only positions a labelled reference line on the volcano. Never
re-impose it when you count, rank, or describe significant proteins, and never justify a
result by "it's ≥2-fold". The moderated t-test tested H0: log2FC = 0, so that is the only
claim the FDR covers; filtering the list on the observed fold change afterwards would
report a stronger claim than the error rate supports, because an observed fold change is
a point estimate with no uncertainty attached. (Testing a fold-change threshold properly
needs `limma::treat()`'s interval null — McCarthy & Smyth 2009 — not a post-hoc cut.) You
may report how many significant proteins also exceed the reference line, as a description
of effect size.

### 8c. Audit the results for common mistakes (surface every issue)
```
python3 scripts/audit_results.py --out AUDIT.md --conditions ./conditions.csv \
    --de-dir ./de_results --acquisition-json /tmp/acq.json --adjp 0.05 --logfc 1
```
Checks for the classic new-user pitfalls: too few/no replicates, imbalanced or
**confounded** design, **mixed acquisition or mixed instruments** in one analysis,
suspiciously low ID depth, very high missingness, contaminant dominance, and DE
results that are too-empty or implausibly-large (batch/normalization artefacts).
**Surface every `WARN` to the user, and STOP on any `FAIL`** (e.g. a group with no
replicate, or a batch confounded with the biology) until they resolve it — don't
let a new user over-interpret a broken design. The findings also go into the
report's "Audit & caveats" section. → detail: `references/audit.md`.

Then run the **biological sample-quality** check (contamination that mimics biology —
hemolysis, tissue cross-contamination, skin) — and **read `references/anomaly-checks.md`
and walk its decision tree** (the judgment branches `audit_results.py` can't automate:
intensity-distribution outliers, replicate/condition ID overlap, mass/RT/CCS shifts,
p-value-distribution shape, largest-FC-in-low-abundance, log2FC>5 artefacts, PCA sample
swaps, GSEA background; plus XL-MS / phospho / non-model branches):
```
python3 scripts/sample_quality.py --matrix ./de_results/Expression_Matrix.csv \
    --conditions ./conditions.csv --report ./search_out/report.parquet --out SAMPLE_QUALITY.md
```
**If a contamination panel is confounded with a group, STOP** — DE may be contamination,
not biology, and protein-level filtering will not fix it (→ `references/anomaly-checks.md`).
**If the sample IS keratin** (nail/hair/wool/skin/feather), pass **`--keratin-sample`** to
`sample_quality.py` **and** `audit_results.py` — keratin is the analyte there, not a
contaminant, and must not be flagged.
Flag every anomaly in the 5-part format (*what · where · why · likely cause · fix*); the
report **always** gets a **Data Quality Notes** section, even if "nothing anomalous observed."

### 8d. Conditional expert review (spawn only on triggers)
Most analyses don't need a panel. If **any** trigger fires — a statistical claim goes to a
collaborator; any `log2FC > 5` or `p < 1e-10`; a new organism/instrument/reagent; cross-tool
discordance; or the user asks — spawn **three reviewers in parallel** (one message, three
`Agent` calls; each gets the report + data paths + context): **proteomics**, **biology**
("could contamination explain this?"), **statistics**. Consolidate their findings (dedupe,
sort by severity) into a `## Expert Review Notes` section and surface every **critical**
issue to the user before finalizing. Skip for format conversion / FASTA prep / exploratory
looks with no formal report. → detail: `references/anomaly-checks.md`.

### 9. Analyze the data (you write the report)
Generate the analysis brief (it lists the figures to embed), then **do the analysis
yourself** — you are the consultant the brief addresses:
```
python3 scripts/analysis_prompt.py --out ANALYSIS_PROMPT.md \
  --de-dir ./de_results --report ./search_out/report.parquet \
  --conditions ./conditions.csv --figures-dir ./figures [--qc ./QC_Metrics.csv] \
  --engine <engine> --acquisition <DIA|DDA> --instrument "<name>" \
  --workflow-manifest ./wf/workflow.manifest.json
```
Then **read `ANALYSIS_PROMPT.md` and every data file + figure it lists, and write a
complete `AI_Analysis_Report.md`** with ALL its OUTPUT sections (Overview, QC, Key
Findings Per Comparison, Cross-Comparison Biomarkers, High-Confidence Biomarkers,
Pathway/GSEA if present, Biological Interpretation, How This Analysis Works, Methods
& Reproducibility) **plus an "Audit & caveats" section from `AUDIT.md`**, a **Data
Quality Notes** section (the 5-part anomalies from steps 8c/8d + `SAMPLE_QUALITY.md`,
even if "nothing anomalous observed"), and — if a review panel ran (8d) — an **Expert
Review Notes** section. **Embed
every figure** (`![caption](figures/<file>.png)`) with an expert interpretation of
what it shows for THIS data. Compute significant proteins, up/down splits,
cross-comparison overlaps, and lowest-CV proteins directly from the CSVs — cite
specific proteins, never fabricate. Make it thorough and expert, like the DE-LIMP
AI export. The brief takes its pipeline description from `de_provenance.json`, so
the report stays correct for whichever engine/method ran.

Then produce **both** delivery formats — a single self-contained HTML page and a
Word document. Both are required; neither is optional.

```
# 1. HTML — the DEFAULT deliverable. QC panels + figures + text in ONE file.
python3 scripts/make_analysis_html.py --session <session> \
    --title "<study name>" --out <session>/output/Analysis_Report.html

# 2. Word — for circulation and track-changes
python3 scripts/to_docx.py --in <session>/output/AI_Analysis_Report.md \
    --out <session>/output/AI_Analysis_Report.docx
```

**Point the user at the HTML first.** Every figure is inlined as a data URI, so it is
one file that opens by double-clicking in any browser on Windows, macOS or Linux, with
no network, no Word, and nothing installed. That matters because the alternative —
a report plus a `figures/` folder — silently loses every image the moment someone
copies just the report, which is exactly what people do. Tell them the filename and
that it is the whole report.

The page puts the **QC panels above the results** on purpose: a volcano plot is equally
persuasive whether or not the run was any good, so a reader who meets the biology first
has already formed a conclusion before seeing the evidence about whether to trust it.
Do not "improve" the order. It also counts significant proteins straight from the DE
CSVs rather than restating the prose, so a disagreement between the two is visible
rather than hidden.

→ detail: `references/analysis.md`.

### 9d. Publication-ready Methods section + acknowledgment
Generate a drop-in LC-MS/MS Methods section straight from the facility raw data,
with the correct UC Davis Proteomics Core instrument-grant acknowledgment:
```
python3 scripts/make_methods.py --raw /path/to/*.d \
    --fasta-meta ./search.fasta.meta.json \
    --out <session>/output/methods.md --de-dir <session>/output/tables
python3 scripts/to_docx.py --in <session>/output/methods.md \
    --out <session>/output/methods.docx
```
It extracts the acquisition parameters from the raw metadata (Bruker `.d`
`analysis.tdf`; Thermo by facility filename prefix), fills the rest from facility
defaults **tagged `[facility default — confirm]`** (the LC column defaults to a
PepSep C18 10 cm × 150 µm, 1.5 µm — override with `--lc-column`), builds a
parameter table showing the source of each value, and appends the instrument's
grant acknowledgment (Fusion Lumos → S10OD021801; Exploris 480 → S10OD026918-01A1;
timsTOF → Dr. Neil Hunter / HHMI). `--fasta-meta` writes the **Sequence database**
paragraph journals require — organism, proteome ID, UniProt release, database type,
entry count, and the contaminant library with its citation. Without it that
paragraph is left blank and tagged, never filled with a guess. **Verify the draft
against the params table and
polish the prose; keep the acknowledgment exact** (confirm wording at the source
URL). This can also be run standalone — just point `--raw` at facility data, no
search/DE needed. → detail: `references/methods.md`.

### 10. Reproducibility bundle (mandatory)
Assemble the bundle that makes the whole analysis reproducible:
```
python3 scripts/provenance.py --outdir ./reproducibility \
  --workflow-manifest ./wf/workflow.manifest.json \
  --setup-json ~/.proteomics-pipeline/setup.json \
  --tools-json ~/.proteomics-pipeline/tools/tools.json \
  --params ./wf/<params_file> --conditions ./conditions.csv \
  --fasta ./search.fasta --fasta-info "$(cat ./search.fasta.meta.json)" \
  --raw /path/to/*.d --report ./search_out/report.parquet --de-dir ./de_results \
  --engine <engine> --de-method <dpc|maxlfq> --contrasts "<...>" \
  --q-cutoff 0.01 --logfc 1.0 --adjp 0.05 \
  --acquisition <DIA|DDA> --organism-taxid <taxid> --instrument "<name>" \
  --commands ./commands.log --timestamp "$(date -u +%FT%TZ)"
```
This writes `reproducibility/` with `run_manifest.json`, `reproduce.sh`,
`REPRODUCE.md`, the conda lock + pip freeze + R sessionInfo + tool versions, a
`skill.txt` recording **which skill produced this and how it was installed**, copies
of the params and conditions, and sha256 checksums of inputs/outputs. This step is
**mandatory and must always run** — code + versions are not optional. `provenance.py`
auto-discovers the env and Rscript even if `--setup-json` is absent, so versions are
always captured; still, **check the returned `skipped` count** and, if the conda
lock / R sessionInfo / checksums were skipped, fix the cause and re-run.

### 11. Output-files report
Catalog everything the run produced so the user knows what each file is:
```
python3 scripts/make_report.py --out OUTPUT_FILES.md \
  --search-out ./search_out --de-dir ./de_results --repro ./reproducibility \
  --extra ./conditions.csv ./search.fasta ./wf ./figures ./AUDIT.md \
          ./AI_Analysis_Report.md ./AI_Analysis_Report.docx
```
`OUTPUT_FILES.md` lists every file (figures, audit, search/DE outputs, the bundle)
with its size and a plain-language description, grouped by purpose, and flags
anything unrecognized (never silently omitted).

### 11a. Comparing two searches of the same files (tool bake-off)
When the same raw files have been searched more than once — DIA-NN library-free vs
FragPipe/diaTracer, two engine versions, two parameter sets — compare the **searches**,
not just the DE tables:
```
python3 scripts/compare_searches.py --out <session>/output/search_comparison --q 0.01 \
  --search "DIA-NN:<sessA>/output/search/report.parquet" \
  --search "diaTracer:<sessB>/output/search/report.parquet"
```
Writes `SEARCH_COMPARISON.md` + CSVs: protein groups and per-run depth, how many proteins
survive in **every** run, median CV across runs, the shared/unique split with Jaccard, and
the per-run correlation of log2 intensity on shared proteins. Unique-to-one-search proteins
are listed in full.

**Do not report "engine X found more proteins" as the verdict.** Depth, completeness and CV
have to be read together — a search can identify more while quantifying them less
reproducibly or in fewer runs. Two caveats to state whenever you present this: the CV is
computed across all runs with no group labels, so it only measures precision when those runs
are replicates; and the comparison assumes protein-group identifiers are comparable, which
fails if the searches used different FASTAs or inference rules. Then run 11b on the DE output
to see whether any of it changes which proteins come out significant.

### 11b. If this is a re-analysis: compare to the original
When step 1b found a prior analysis of the same dataset, compare the two with the
Comparator (a faithful port of DE-LIMP's Run Comparator):
```
Rscript scripts/compare_analyses.R --out <session>/output/comparison \
  --adjp 0.05 --logfc 0 \
  --analysis "Original:<prior>/output/tables" \
  --analysis "Reanalysis:<this session>/output/tables"
```
It writes `COMPARISON.md` + CSVs: protein-universe overlap, the 3×3 Up/Down/NS
concordance per shared contrast, direction concordance on co-significant proteins,
and logFC correlation. (`session.py finalize` already wrote `DIFFERENCES.md` —
what *settings* changed; the Comparator shows how the *results* changed.)

### 11c. Cross-tool comparison — report it, don't just compute it

Comparing two search engines on the same raw files is a distinct deliverable from a DE
re-analysis. `compare_searches.py` compares the SEARCHES (what each found, in how many
runs, how reproducibly); `compare_analyses.R` compares the finished DE TABLES (which
proteins actually change). Run both, then build the report:
```
python3 scripts/make_comparison_report.py --out <session>/output/comparison_report \
  --search-comparison <dir from compare_searches.py> \
  --de-comparison     <dir from compare_analyses.R> \
  --qc "EngineA:<sessionA>/output/tables/QC_detected_vs_inferred.csv" \
  --qc "EngineB:<sessionB>/output/tables/QC_detected_vs_inferred.csv" \
  --instrument "<detected instrument>" --title "<A> vs <B>"
python3 scripts/to_docx.py --in <...>/COMPARISON_REPORT.md --out <...>/COMPARISON_REPORT.docx
```
**Read `references/cross-tool-comparison.md` BEFORE interpreting the output** — it holds
the design rules (both engines through the same DE pipeline; name every uncontrolled
difference) and the interpretation patterns a bake-off keeps producing: more precursors
with fewer proteins means protein *inference*, not sensitivity; a high intensity r with a
much lower logFC r means direction is trustworthy and magnitude is not; and the subtler
the contrast, the more the engine choice decides the answer. It also says to **load the
`dataviz` skill and re-run its palette validator** before changing any chart.

It emits an **interactive standalone HTML** (metric toggles, hover tooltips, a 3x3
concordance matrix, sortable tables, theme-aware, no CDN) plus a markdown skeleton in
the section order DE-LIMP's Run Comparator prompt defines — `docs/AI_PROMPTS.md` §1:
Factual Observations, Sources of Disagreement, Case for A, Case for B, Settings Audit,
Concordant Proteins, Synthesis, Follow-ups. **The generator writes every number; you
write the judgement** — fill each `TODO(model)` block from the data and never invent a
figure. Follow the prompt's guidelines: neither tool is inherently superior, every claim
carries a number or a citation, speculation is labelled as such.

**To compare engines fairly, run BOTH through the same DE method.** Push each engine's
report through `run_de.R` with identical `--method`, contrasts and thresholds, so the
search is the only variable. `run_de.R --format tsv` reads a FragPipe `report.tsv`
directly. Also state any engine-VERSION gap out loud — FragPipe 24 bundles DIA-NN
1.8.2 beta 8, two majors behind what the `diann_*` workflows pin, and that alone moves
the numbers.

### 12. Finalize the session + report to the user
```
python3 scripts/session.py finalize --dir <session> --zip
```
This tidies `output/` (tables→`tables/`, figures→`figures/`), writes the session
`README.md` (and `DIFFERENCES.md` for a re-analysis), and zips the session for easy
sharing. Then summarize: data type (instrument + acquisition), engine + **pinned
version**, mass accuracy **and its source**, the skill's `defaults_version`, FASTA
source, DE method, and per-contrast significant counts. Point them at the
**session folder** and its `README.md`, then `AI_Analysis_Report.md` (the
interpretation), `OUTPUT_FILES.md` (what every file is), `tables/methods.txt`
verbatim (the Methods paragraph — don't paraphrase), **`tables/reproducibility_log.R`
(the analysis as plain R — say this is where the code is; it is what most people
mean when they ask)**, `reproducibility/REPRODUCE.md` (the pinned recipe, for
re-running the search too), and — for a re-analysis — `DIFFERENCES.md` + the
`comparison/COMPARISON.md`.

## When something is missing
- Anything in the env (R, limpa, Sage, pyarrow) → re-run `setup.sh`; relay its
  `notes`. Never tell the biologist to "install R/limma/..." by hand.
- macOS + DIA data + no Docker → run `build_diann_docker.sh` and relay its exact
  steps (install Docker Desktop, open once). Don't silently switch to a DDA engine.
- macOS + DDA + Bruker/Thermo → msconvert is Linux-only; see `references/install.md`
  ("macOS + Sage"). Prefer DIA-NN if the data is DIA; else convert to mzML first.
- Unrecognised instrument → parameters fall back to automatic mass calibration;
  say so plainly. This is never a stop: no organism or instrument is "unsupported".
- FragPipe license (MSFragger/IonQuant) → surface the `tools.json` note and the fix.
- Acquisition/instrument unknown → ask the user; never guess into a multi-hour run.
