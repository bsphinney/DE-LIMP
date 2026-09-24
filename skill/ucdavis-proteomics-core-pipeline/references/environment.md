# Environment & acquisition reference

How the skill adapts to where it's running, plus FASTA resolution.

## Platform classes (`detect_env.sh`)

| class | how detected | engine acquisition | execution |
|---|---|---|---|
| `hpc` | `sbatch` on PATH, or `/quobyte/proteomics-grp` visible | on HIVE, reuse the Core's native DIA-NN builds (below); elsewhere download the Academia Linux build. Sage from the conda env, else its release binary | **submit via sbatch** (never login-node) |
| `mac` | `uname` = darwin | DIA-NN via **Docker** (no native mac build); Sage native | inline |
| `linux` | everything else | native binaries | inline |

`uc_davis_hive` is true when `/quobyte/proteomics-grp` exists → enables FASTA reuse
and the Core's DIA-NN builds.

Container runtime preference: hpc→apptainer, mac→docker, linux→native.

## DIA-NN by environment

- **No native macOS build exists.** On mac you must run DIA-NN through Docker. Set
  `DIANN_DOCKER_IMAGE` to a built image, or build one from the Academia Linux zip's
  bundled Dockerfile. `acquire_tools.sh` writes a note when this is unresolved.
- **HIVE (Proteomics Core):** DIA-NN is kept under `/quobyte/proteomics-grp/dia-nn/`
  as **native builds** at `build_<nnn>/diann-<version>/diann-linux` — **2.5.1, 2.6.0,
  2.6.1 and 2.7.0** (`build_270`, added 2026-09-23) — plus one older Apptainer image,
  `diann_2.3.0.sif`. The skill pins **2.7.0** (`resolve_defaults.py`). `acquire_tools.sh` resolves the
  **pinned** version by looking for `build_*/diann-<version>/diann-linux` first, then a
  version-matched `.sif`, and **never silently substitutes a different version**
  (reproducibility): a pin that is not there — as 2.7.0 was for the FRAN pilot, before
  `build_270` existed — is downloaded as the Academia Linux build into your own tools root instead. Unpinned, it
  takes the highest *version* the build paths name. Either way `tools.json` records the
  version of the build it found, never `latest`. Listing it yourself is cheap and allowed
  on the login node: `ls -d /quobyte/proteomics-grp/dia-nn/build_*/diann-*`. The
  facility's `run_diann_*.sbatch` in that folder is the reference invocation. AlphaDIA is
  also on HIVE at `/quobyte/proteomics-grp/apptainers/alphadia.sif` (auto-reused).
  (`DIANN_HIVE_DIR` overrides the build directory, and `RADIANT_HIVE_DIRS` the
  colon-separated folders searched for a Radiant `.sif` — by default
  `/quobyte/proteomics-grp/apptainers` then `/quobyte/proteomics-grp/radiant`; the tests
  use both to mock the layout.)
- **Linux native:** needs glibc ≥ Linux Mint 21.2 and .NET 8. If missing, prefer
  Docker/Apptainer.
- DIA-NN reads `.raw`/`.d` natively from 2.1+. On Linux, `.raw` also needs a **.NET 8
  runtime ≥ 8.0.17** (`ensure_dotnet8.sh`; see `references/search-engines.md`).
- **ThermoRawFileParser on HIVE:** the Core keeps TRFP 2.0.0.0 at
  `/quobyte/proteomics-grp/tools/ThermoRawFileParser/ThermoRawFileParser`, and
  `detect_acquisition.py` uses it when no parser is set, on PATH or in the pipeline env
  (`$THERMORAWFILEPARSER_SHARED` overrides the location). It is **framework-dependent**:
  it needs `Microsoft.NETCore.App` 8 **and** `Microsoft.AspNetCore.App` 8, which
  `ensure_dotnet8.sh` installs together (or `module load dotnet-core-sdk/8.0.4` for the
  parser alone — 8.0.4 is too old for DIA-NN). bioconda's `thermorawfileparser`, which
  `setup.sh` puts in the env, is self-contained and needs neither. Either parser's folder holds the
  `ThermoFisher.CommonCore.RawFileReader`/`.Data` DLLs (8.0.6, .NET 8) that
  `thermo_resolution.py` loads through pythonnet to read the Orbitrap resolution from the
  scan trailer — that needs a `Microsoft.NETCore.App` 8 root even with a self-contained
  parser. Detail: `references/install.md`.

## Version pinning (reproducibility)

`acquire_tools.sh` honors `PIN_ENGINE`/`PIN_VERSION` from the workflow bundle and
caches under `~/.proteomics-pipeline/tools/<engine>/<version>/`. Different pinned
versions coexist. The written `tools.json` records `pinned` (what was asked for) and
`versions` (the build each command **is**) so the report can state exactly what ran.
**Always pass the bundle's engine+version** — a result from "latest" is not
reproducible.

`versions` never says `latest`: an unpinned run records the number it resolved to — for
DIA-NN from the build path, `.sif` name, download asset or Docker image tag; for Sage from
the release tarball kept beside the binary, **not** `sage --version`, which prints 0.14.6
for the v0.14.7 release; for FragPipe from the unpacked folder (`fragpipe-24.0/`) or the
release asset (`FragPipe-24.0.zip`); for AlphaDIA from `pip show alphadia`; for Radiant from
an image tag or the `.sif` name. When no version can be determined it records `""` with a
note — an unpinned Radiant image (`seerbio/radiant-fulcrum:latest`) is one such case, and so
is a cached Sage whose release tarball is gone. It is never the **request**: echoing
`PIN_VERSION` back recorded the manifest's own pin as the build that ran, which
`run_search.py` could then only compare with itself.

There is a key per engine — `diann`, `sage`, `radiant`, `fragpipe`, `alphadia`. A *missing*
key is worse than an empty one: it makes every search that engine ever runs record
`version: null`, since `run_search.py` has nothing to read.

For a Radiant image on mac/linux the recorded release is the image **tag**
(`seerbio/radiant-fulcrum:2.3.3`), which is what will be pulled — but nothing here pulls it,
so the tag is the only evidence of what is inside; the note beside it says so.

The manifest and `tools.json` can still disagree — `resolve_defaults.py` pins one version,
but `tools.json` holds whatever was last acquired (the FRAN pilot: manifest 2.6.1,
`tools.json` 2.7.0). Nothing forces
them to match, so `run_search.py` compares them before submitting and records both in
`search_provenance.json` — as one `engine_version` object shaped like `scan_window` (a
`value`, and a `source` saying in words where it came from), plus top-level `version`, which
always equals `engine_version.value`:

| field | meaning |
|---|---|
| `version` = `engine_version.value` | the build that runs: `tools.json` `versions`, else the one version the **command itself** names (a build folder such as `diann-2.6.0/`, a `.sif` name, an image tag); else `null`. **Never the manifest's pin** — `fran_deposit.py` sends `version` to FRAN as the engine version it stores, so a pin nothing confirmed would be stored there as fact |
| `engine_version.source` | where `value` came from: `tools.json versions.<engine>, …`, `named by the command tools.json runs (…)`, or `unknown -- <why>` |
| `engine_version.tools_json` | `tools.json` `versions.<engine>` as written, even `latest` / `env` |
| `engine_version.named_by_command` | **every** version the command names (usually zero or one; two means the command disagrees with itself, and `value` is then `null`) |
| `engine_version.manifest_pin` | the manifest's pin as written — `null` unless the manifest's `engine.name` **is** this engine. A block naming no engine pins nothing |
| `engine_version.mismatch` | `true`/`false` when `value` and the pin are both known; **`null`** when they cannot be compared. A WARNING is printed when `true` |

Only something **shaped like a version** is read as one — `2.6.1`, `v0.14.7`, `24.0`.
Anything else is no more a build than `latest` is: `env`, `nightly`, `dev`, `n/a`, `TBD`, a
date. The check is the shape of an answer, not a list of known wrong answers, because
`tools.json` is written by a shell script whose every "cannot tell" branch is one edit away
from a new word. A refused string is kept verbatim in `engine_version.tools_json` and named
in `source` (`unknown -- tools.json records sage 'nightly', which is not a version`), so
nothing is hidden — it just never reaches `version`.

A `tools.json` whose `versions` contradicts its own command (says 2.6.1, runs
`diann-2.6.0/diann-linux`) records `version: null` with a WARNING: one is wrong and nothing
tells which. A version *directory* (`sage/0.14.6/sage`, `fragpipe/24.0/`) is not read,
because `acquire_tools.sh` names cache folders after the request, not the build inside.

A mismatch WARNING means the version the user confirmed (SKILL.md golden rule #1) is not
the one about to run: say so before submitting, and either re-acquire the pin (the WARNING
prints the command, with this `tools.json`'s platform and tools root; on macOS it rebuilds
the Docker image) or get the new version confirmed. `version: null` means the search is
recorded — and deposited to FRAN — with no engine version; for a `sage` on PATH that is
currently the only outcome.

## FASTA resolution (`fetch_fasta.py`)

### `resolve` — organism → proteome (always run this; never guess a UPID)
`fetch_fasta.py resolve --organism "mouse"` returns ranked `candidates` +
`selected` + `needs_menu` + `notes`. Accepts a common name, a scientific name, an
NCBI taxid, or a `UP…` accession (`--taxid` also works). Show the user `selected`
and confirm. When `needs_menu` is true, present the menu instead of auto-picking.

`scripts/resolve_organism.py` is a thin alias for this same resolver, kept for
older call sites. **One implementation** — don't add a third.

How it picks, and why each rule exists:
- **Curated table first.** 18 organisms a core facility actually sees map
  name/alias → *taxid*, and the proteome is then resolved live from that taxid.
  Free-text search answers some queries badly: "Escherichia coli K-12" returns
  five Non-Reference MG1655 assemblies and never surfaces the real reference
  `UP000000625` at all. Only the taxid is curated — a pinned `UP…` accession goes
  stale (the dog reference moved from `UP000002254` to `UP000805418`), so the
  table's accession is used only as an offline fallback and a staleness check,
  which is reported in `notes`.
- **Exact proteome-type match.** `"reference" in "Non Reference proteome"` is
  True, so a substring test ranks strain assemblies as references.
- **Strain qualifier stripped when name-matching.** `scientificName` is
  `Saccharomyces cerevisiae (strain ATCC 204508 / S288c)`, so the bare species
  name matched nothing and ranking fell through to protein count — which put
  *S. pastorianus* first.
- **"Excluded" proteomes dropped** unless nothing else matches.
- **Auto-picks only when unambiguous**: one reference proteome the user actually
  named. "mouse" also matches *Myotis myotis* and mouse-ear cress; "baker's yeast"
  matches two S. cerevisiae strains → menu.

Common IDs (still confirm with `resolve`): human `UP000005640`, mouse
`UP000000589`, rat `UP000002494`, yeast `UP000002311`, E. coli `UP000000625`.

### `fetch` — proteome → search FASTA
Priority, cheapest/most-trusted first:
1. `--path` override → used verbatim (pre-staged proteome).
2. **HIVE** (`--hive`): reuse `/quobyte/proteomics-grp/MRS/`. Matches only files
   whose name starts with the proteome ID and skips `*_plus_*contam*` /
   `*decoy*` / `*predicted*` variants — appending contaminants to a database that
   already contains them would duplicate them.
3. **UniProt.**

### Database type (`--content`, default `one_per_gene`)
| value | what you get | human size |
|---|---|---|
| `one_per_gene` | canonical, one protein per gene — **the default** | 20,652 |
| `reviewed` | Swiss-Prot only | ~20,400 |
| `full` | + unreviewed TrEMBL | **147,506** |
| `*_isoforms` | + splice isoforms | larger still |

**`one_per_gene` only comes from the reference-proteome FTP tree**
(`{Kingdom}/{UPID}/{UPID}_{TAXID}.fasta.gz`), because UniProt's REST
`&onePerGene=true` is **silently ignored** — verified 2026-07-29: the yeast stream
returns byte-identical output (6,067 entries) with and without it. A plain REST
`(proteome:X)` stream is therefore always the *full* set. This is why DE-LIMP
(`R/helpers_search.R`) uses FTP, and why the Core's staged HIVE database is
`UP000005640_9606.fasta` = 20,663 sequences, not 147k.

If no FTP file exists (non-reference proteome), the script warns loudly, records
the warning in its output, and falls back to the REST full set — it never swaps
databases silently. If REST also fails it exits with the reason.

### Contaminants (`--contaminants`, default `universal`)
Sets: `universal` (default; what the Core stages on HIVE), `cell_culture`,
`mouse_tissue`, `rat_tissue`, `neuron_culture`, `stem_cell_culture`, `none`.
Source order: `--contaminants-path` → HIVE (matching the *requested* set) →
a DE-LIMP checkout's `contaminants/` → the Hao lab GitHub repo (public;
Frankenfield et al. 2022, JPR 21(9):2104-2113, doi:10.1021/acs.jproteome.2c00145).

Headers are `Cont_`-tagged, so DIA-NN's `--cont-quant-exclude Cont_` keeps them
out of quantification and normalisation. The DIA-NN **Linux binary ships no
contaminant FASTA of its own** (the GUI's "Contaminants" checkbox is a
Windows-side asset), so they must be appended here. `fetch_fasta.py` reports the
tag as `diann_cont_quant_exclude`; pass the sidecar to `estimate_params.py
--fasta-meta` and the flag lands in the cfg automatically.

**Contaminant entries that ARE target proteins are removed.** Matching is by sequence, not
accession: an entry identical to a target entry, or an exact substring of one (≥ 7 aa, DIA-NN's
default minimum peptide length; I and L kept distinct), is dropped so the protein is quantified
under its own accession, and recorded under `contaminants_dropped_as_target` (sidecar) with a
warning and a methods sentence. The universal set holds 152 human-identical entries (human
keratins; bovine ACTB/EEF1A1/YWHAZ/tubulins) + 1 substring, and 31 mouse-identical ones;
left in, DIA-NN reported ACTB, EEF1A1 and KRT8 only as `Cont_` groups and
`--cont-quant-exclude` removed them from quant (a real HeLa search, 2026-09-23: 6.9% of all
intensity). The digestion enzyme(s) actually used (`fetch --enzyme`, default `trypsin,lysc`)
are kept as contaminants even when they match a target protein — they are reagents — and
recorded under `contaminants_kept_despite_target_match`; any other protease entry follows the
normal rule (S. aureus's own SspA is identical to the Glu-C entry and stays quantified on a
trypsin digest). The
auditors flag the dropped proteins as "possible contamination, kept in quantification".

**A failure to fetch contaminants is fatal, not a warning** — the GPM cRAP URL
this script used previously now 404s, and the old warn-and-continue behaviour
meant searches silently ran with no contaminants at all, which also made the
contaminant-dominance QC check meaningless. Override deliberately with
`--contaminants none` or `--allow-missing-contaminants`.

### Output
Refuses to proceed on 0 sequences, and warns if a full-set download comes back
>5% short of UniProt's declared count (truncated stream). Writes
`<out>.meta.json` alongside the FASTA — sha256, source URL, organism, taxid,
content type, UniProt release, per-part sequence counts, contaminant set +
citation, and any warnings. Pass it to `provenance.py --fasta-info` so
`reproduce.sh` rebuilds the database that was *actually searched*.

## SLURM submission (hpc)

`run_search.py --sbatch job.sh` emits a login-node-safe script (64G, 12h). **No queue is
hard-coded**: `run_search.slurm_queue()` — the one definition, also used by
`diann_parallel.py`, `radiant_parallel.py` and `diatracer_parallel.py` — asks SLURM what
the submitting user may use (`sacctmgr show assoc user=$USER`) and picks:

1. an explicit queue (pass `--partition` **and** `--account` together, plus `--qos` if the
   association has one), which skips everything below;
2. `genome-center-grp` on `high` (facility members: not preemptible);
3. `publicgrp` on `low` (everyone else: preemptible, so `#SBATCH --requeue` is added).

When an account has both, **utilisation** decides, not entitlement. `slurm_queue()` counts
the CPUs *you* already run on `high` (`squeue`) against the 64-CPU per-user cap, and the idle
CPUs on `low` (`sinfo`), and compares both with `need = min(peak, 16)`, where `peak` is what
the caller passes (16 when it passes none) — at most one array task's worth, **not** the
job's own request. It does this **once, when the scripts are generated**, not when each job
starts. Two rules:

- **A (every job):** `low` when fewer than `need` of your CPUs are free on `high` and `low`
  has at least `need` idle.
- **B (only where a caller marks a step preemption-safe):** also `low`, given the same `need`
  idle there, when fewer than `2 × need` are free on `high`, or when your usage there
  cannot be read (`squeue` cannot run).

| route | how often the queue is decided | rules |
|---|---|---|
| DIA-NN 5-step chain (`diann_parallel.py`) | **once, for all five steps**, with `need` 16 | A only |
| Radiant per-file array (`radiant_parallel.py`) | step 2 (the array) apart from steps 1 and 3 | step 2: A and B, `need = min(--threads-per-file, 16)`; steps 1 and 3: A, `need = min(--fulcrum-cpus, 16)` |
| diaTracer (`diatracer_parallel.py`) | once, `need` 16 | A only |
| single-script `--sbatch` | once per script, `need` 16 | A only |

So the DIA-NN chain's array steps (2 and 4) **always share the queue of steps 1, 3 and 5**:
they do not move to `low` on their own when `high` is merely busy or `squeue` is missing.
(`diann_parallel.py` does ask `slurm_queue()` for a separate array-step queue, but it has
already filled in partition and account by then, and an explicit pair is returned
unchanged.) The whole chain moves to `low` only under rule A. Because `need` is capped at
16, the chain's 64-CPU steps (3 and 5 by default; step 1 asks for 16) **go to `high`
whenever 16 or more of your CPUs were free there at generation, and then wait on `high`**
for the rest; they do not move to `low` because the 64 they ask for are unavailable. To put
the chain on `low`, pass `--partition low --account publicgrp` (the generator adds
`--qos=publicgrp-low-qos`).

If associations cannot be read at all the script falls back to `publicgrp/low`, **not** the
cluster default: that is `high`, which rejects a non-facility account.

Checked on HIVE 2026-09-16 (`sacctmgr show assoc` / `show qos`, `sinfo`): the default
partition is `high`; `genome-center-grp` has `high` (per-user cap 64 CPUs) and `gpu-a100`
but no `low`; `publicgrp` has `high` (8 CPUs / 128 GB **per job**) and `low` (no per-job
cap, preemptible).

`run_search.py` passes `--partition/--account/--qos` on to the job-array routes (the DIA-NN
5-step chain and the Radiant per-file array). Whether the single-script `--sbatch` paths —
DIA-NN at ≤5 files or when the chain declines, Sage, FragPipe, AlphaDIA, single-file
Radiant — honour them has changed between skill versions, so on every route **read the
`#SBATCH` header of the generated script before submitting**: it is the queue the job will
use. Submit with `sbatch job.sh` (or `bash <out>/submit.sh` for a two-job or chained
search), poll the `<job>_<id>.log`, then run `run_search.py --adapt-only` for Sage/FragPipe
to build `report.parquet`.

A DIA-NN search of more than 5 files on a SLURM host routes to the 5-step parallel chain
instead, which has no single job script: `job.sh` is not written, an existing regular file of
that name is renamed to `job.sh.stale-<time>` once the chain is generated, and
`run_search.py` exits 3. Submit `<out>/submit.sh`
(→ `references/diann_parallel.md`).
