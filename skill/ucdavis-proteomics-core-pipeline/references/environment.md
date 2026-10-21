# Environment & acquisition reference

How the skill adapts to where it's running, plus FASTA resolution.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 17: Platform classes (`detect_env.sh`)
- 30: DIA-NN by environment
- 71: The R stack: limpa >= 1.4.0 (`setup.sh`)
- 147: Version pinning (reproducibility)
- 258: FASTA resolution (`fetch_fasta.py`)
- 443: SLURM submission (hpc)
- 513: 5. Acquiring the pinned engine (step 5)
- 546: 6. Building the FASTA (step 6)
<!-- /contents -->

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
- > Superseded by "5. Acquiring the pinned engine (step 5)" below.

  **HIVE (Proteomics Core):** DIA-NN is kept under `/quobyte/proteomics-grp/dia-nn/`
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

## The R stack: limpa >= 1.4.0 (`setup.sh`)

`run_de.R --method dpc` reads the DIA-NN report with `limpa::readDIANN()`. The pinned stack is
**limpa >= 1.4.0** (Bioconductor 3.23). bioconda's only build is limpa **1.2.5**
(Bioconductor 3.22, R 4.5), and a 2.8.0 `run_de.R` died on it (HIVE sbatch 24154221:
`unused argument (annotation.columns = dpc_ann)` — 1.2.x calls that argument `extra.columns`).

- **What `setup.sh` builds.** A conda env, solved with `--override-channels -c conda-forge -c
  bioconda` (a `~/.condarc` listing `defaults` with `channel_priority: strict` otherwise left
  libmamba only defaults' R ≤ 4.3 `r-statmod`, and the solve failed: sbatch 24154212), with
  R 4.5 + bioconda `bioconductor-limma` 3.66 + `r-statmod` / `r-data.table` /
  `r-nanoparquet`. Step 2a then installs **limpa from the Bioconductor 3.23 source
  repository** into that env (`PROTEOMICS_LIMPA_BIOC` overrides the release). limpa is pure R
  (no `src/`), has no R-version floor, and every limma function it imports is in limma 3.66,
  so no compiler and no second Bioconductor stack are needed — conda-forge has `r-base` 4.6
  but no R 4.6 builds of `r-statmod`, `r-arrow` … yet. Measured on the mouse-brain IP cohort v2 (HIVE sbatch
  24176224): limpa 1.4.0 on the setup.sh env (limma 3.66) reproduces the delivered v2 tables
  (limpa 1.4.0 / limma 3.68.5) exactly — every contrast's significant set identical.
- **The check.** The last thing `setup.sh` does is verify `packageVersion("limpa") >=
  "1.4.0"`. `setup.json` → `limpa` = `{version, required, ok, source}`; when it fails,
  `notes` and stderr say `ERROR: limpa <v> ...` with the fix (re-run `setup.sh`, or an
  `install.packages()` line — on a cluster inside a job, never on a login node; HIVE's compute
  nodes reach Bioconductor and CRAN), and `setup.sh` exits 1 (`--check` only
  reports).
- **nanoparquet is checked too.** limpa reads DIA-NN's `report.parquet` with it (`nanoparquet
  package required but is not installed`), but it is only in limpa's *Suggests*, so installing
  limpa never brings it, and `create_env` lists `r-nanoparquet` only since 2.8.0 and only runs
  for a missing env — an env built before that never got it (2.8.0 report). `setup.sh` now
  checks that R can load it; if not, it installs **conda-forge's `r-nanoparquet`** into the
  existing env, **in its own call with `--freeze-installed`**: conda-forge builds it only for
  R < 4.6, and solved together with the parser / pythonnet / pandas without the freeze, an R
  4.6.1 env got the plan "r-base 4.6.1 → 4.5.3" (2.12 review; a dry run on HIVE's
  `proteomics-pipeline-r46`). Frozen, conda can only add it or fail. **R ≥ 4.6 has no conda
  build**, so there it is CRAN's source package, compiled — only when the compilers this R
  names (`R CMD config CC` / `CXX`: conda's R names conda-forge's, e.g.
  `x86_64-conda-linux-gnu-c++`) and `make` run; that is why a plain CRAN install once failed
  with "x86_64-conda-linux-gnu-c++: not found". Without them the fix first installs
  conda-forge's `c-compiler cxx-compiler make` (frozen), then the CRAN line — how HIVE's
  `proteomics-pipeline-r46` got nanoparquet 0.5.2. `setup.json` → `nanoparquet` = `{ok, fix}`
  (the fix is the command that works for THIS env's R), `ready_for.de` is false without it
  when limpa is ≥ 1.4.0 (`de_reader`, below), and
  `notes` / stderr say `ERROR: nanoparquet ...` with it; `setup.sh` exits 1 (`--check` only
  reports). The limpa fix line above leaves nanoparquet as it is. With no env built yet,
  `ready_for.dia` / `dda` can be true (an Rscript outside the env) while `de` is false, and a
  note says why: the DE runs in the env, so run `setup.sh` before the DE step.
- **An env that is still on limpa 1.2.x** (an old `setup.sh`, no internet) is not a dead
  end: `run_de.R` passes the annotation columns under whichever name this limpa's
  `readDIANN()` has (`limpa_compat.R`) and records it — `de_provenance.json` →
  `limpa_read = {limpa_version, annotation_argument, path}`. limpa 1.2.x reads `report.parquet`
  with **arrow** (`readDIANN`: `arrow::read_parquet`, RELEASE_3_22); only limpa ≥ 1.4.0 reads it
  with nanoparquet (RELEASE_3_23, where arrow left `Suggests`). So `setup.json` →
  `de_reader = {limpa, reads_parquet_with, ok}`, and `ready_for.de` is true on limpa 1.2.x
  without nanoparquet as long as arrow loads (run_de.R also filters the report with arrow on
  both methods). What such a DE loses: the pinned limpa, and run_de.R's contaminant share, which
  reads `Precursor.Quantity` with nanoparquet and falls back to `Precursor.Normalised`
  (`de_provenance.json` → `contaminants.intensity_column`). Seen on HIVE on 2026-10-16: the
  account's env (limpa 1.2.5, no nanoparquet) ran every DE that week while `setup.sh --check` said
  "no DE can run". The upgrade is still the fix: re-run `setup.sh` — 1.2.x and 1.4 do NOT give
  identical numbers (next point), so an env other analyses run on is upgraded between analyses.
- **Results depend on the limpa version, through `dpcCN()`'s row subset.** On more than
  2,000 precursors `dpcCN()` fits the detection-probability curve on 2,000 rows: in limpa
  1.2.x a SEEDED RANDOM sample (`set.seed(20250620)`, `sample.int`); in 1.4 a SYSTEMATIC one
  (rows ordered by missingness, then mean; every n/2000-th) — github.com/bioc/limpa
  RELEASE_3_22 vs RELEASE_3_23 `R/dpcCN.R`. The curve, and so every protein value and
  p-value, moves a little. Measured on the mouse-brain IP cohort v2 (83,587 precursors, 30 runs, `--block
  Mouse`, same env, sbatch 24176224): DPC (β0, β1) = (−13.200, 1.184) on 1.2.5 vs (−13.072,
  1.178) on 1.4.0; protein matrix |Δ| median 0.010, max 0.052 log2; within-mouse correlation
  0.1883 vs 0.1882; per contrast |ΔlogFC| median 0.002–0.006, max 0.064; significant sets
  Jaccard 0.958–1.000 (1.2.5 calls 0–4 more per contrast, e.g. RyR-vs-IgG Old 71 vs 68).
  Small, but not zero: `de_provenance.json` records `limpa_version` (and `packages.limpa`)
  for every run — **compare analyses only within one limpa version**, and re-run the older
  one rather than diffing across versions. On small inputs (≤ 2,000 precursors) no subset is
  taken and the two versions agree exactly.
- **R 4.6 + Bioconductor 3.23 throughout** (HIVE's `proteomics-pipeline-r46`, built with
  conda-forge `r-base` 4.6 + `BiocManager` 3.23) passes the same check.

## Version pinning (reproducibility)

> Superseded by "5. Acquiring the pinned engine (step 5)" below.

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

> Superseded by "5. Acquiring the pinned engine (step 5)" below.

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

> Superseded by "5. Acquiring the pinned engine (step 5)" below.

A mismatch WARNING means the version the user confirmed (SKILL.md golden rule #1) is not
the one about to run: say so before submitting, and either re-acquire the pin (the WARNING
prints the command, with this `tools.json`'s platform and tools root; on macOS it rebuilds
the Docker image) or get the new version confirmed. `version: null` means the search is
recorded — and deposited to FRAN — with no engine version; for a `sage` on PATH that is
currently the only outcome.

### search_provenance.json: `scan_window.mode` / `mass_acc.mode` (stable)

Every DIA-NN search's `search_provenance.json` has a top-level `scan_window` and `mass_acc`
record, each with a **`mode`**: one machine-readable word for how the value was set. The prose
(`source`, `reason`) may be reworded; **these values are stable**. FRAN ingests them — from
`fran_manifest.json`'s copy of the provenance — as a variable of its DIA-NN vs Spectronaut
comparison. A value is never renamed or reused; a new case gets a new value, added in
`diann_parallel.SCAN_WINDOW_MODES` / `MASS_ACC_MODES` (the one definition), here, and in the
test that pins them (`tests/test_probe_estale.py`, `ProvenanceModeTests`). The same records
appear under `result` where the search itself wrote one.

| `scan_window.mode` | meaning |
|---|---|
| `measured` | measured on these runs by step 1b (the 5-step chain) and pinned for every step. Written at generation as the plan; `value` is null and `value_file` (`window.txt`) holds the radius once step 1b has run |
| `fallback_auto` | step 1b's measurement failed (retried once): DIA-NN chose the radius itself, for each run. `probe_fallback` holds the `reason` |
| `pinned` | given in the cfg and passed to every step (`value`) |
| `auto` | not set, by design — a single-shot search, or DDA: DIA-NN chose the radius itself |
| `invalid` | the cfg passes a `--window` that is not one positive integer — `0` included: DIA-NN logs `scan window radius should be a positive integer` and then chooses a radius per file (the poplar run), so `0` is a rejected value, not `auto` (which means the flag was not set). What DIA-NN does with any other non-integer is unverified |
| `unknown` | the cfg could not be read |

| `mass_acc.mode` | meaning |
|---|---|
| `measured` | measured on these runs before the search (step 1b, or the single-shot search's probe): a measured level floored at the SOP, a documented level as documented. Written at generation as the plan; `value_file` (`massacc.txt`) holds the pair |
| `fallback_default` | that measurement failed (retried once): the documented level as given, the other at the facility SOP — `default` names the DEFAULT levels; `probe_fallback` holds the `reason` |
| `pinned` | given in the cfg |
| `pinned_default` | pinned by `estimate_params.py` at the facility SOP for a level that cannot be measured (DDA); `default` names the levels |
| `auto` | neither level set: DIA-NN optimised it itself |
| `partial` | one level given in the cfg: DIA-NN 2.7.0 then fixes both, the other at 20 ppm |
| `invalid` | the cfg passes a value DIA-NN does not read as a tolerance: `0` (a literal 0 ppm on the command line — "Mass accuracy will be fixed to 0 (MS2) and 0 (MS1)", 0 IDs on a 28-run Lumos search — not "automatic"), negative, not a number, or set twice with different values |
| `unknown` | the cfg could not be read |

`0` is `invalid` for **both** flags. DIA-NN's README says these settings "are set to 0, meaning
that DIA-NN will optimise them automatically" — that describes its GUI fields, which leave the
flag out at 0; a literal `0` on the command line is a rejected radius (`--window`) or a 0 ppm
tolerance (`--mass-acc`), never `auto`.

A `measured` plan whose step 1b then fails outright (a refused mass accuracy, a signal) leaves no
report, so it is never deposited as `measured`. "Was the window measured?" is `scan_window.mode
== "measured"`; "did DIA-NN choose it?" is `mode in ("fallback_auto", "auto")`.

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

> Superseded by "6. Building the FASTA (step 6)" below.

Priority, cheapest/most-trusted first:
1. `--path` override → used verbatim (pre-staged proteome). Its organism is the user's answer:
   `--organism '<scientific name>' --taxid <taxid>` (or `--organism none` for a database with
   no single organism) is required, and recorded with `organism_source: user (...)`. A name and
   a taxid the curated table says are different organisms (different genus: `Homo sapiens`
   with 10090) are refused.
2. **HIVE** (`--hive`): reuse `/quobyte/proteomics-grp/MRS/`. Matches only files
   whose name starts with the proteome ID and skips `*_plus_*contam*` /
   `*decoy*` / `*predicted*` variants — appending contaminants to a database that
   already contains them would duplicate them.
3. **UniProt.**

### Database type (`--content`, default `one_per_gene`)

> Superseded by "6. Building the FASTA (step 6)" below.

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

**Only a whole file is used.** Each download is checked three ways: the bytes against the
server's Content-Length, the gzip stream read to its end (its CRC is checked there), and the
entry count against UniProt's `geneCount` for the proteome (within 5%). Measured 2026-10-01
(release 2026_03): the FTP file's count EQUALS geneCount for human (20,652), mouse (21,860),
yeast S288c (6,066), E. coli K-12 (4,403) and Arabidopsis (27,496); the 5% only absorbs REST and
FTP serving different releases around a release day. A failed attempt is retried twice
(after 5 s and 20 s); the sidecar's `download_check` records the bytes, the count and the
attempts. A transfer that is still cut short, or a count that does not fit, **stops with exit
≠ 0** (a staff search, 2026-10-01: a truncated gzip used to fall back to the 147,520-entry full set
with exit 0, and the re-run a minute later was whole). Re-run, or pass `--content full` to
search the full set on purpose.

Only when the FTP server answers **404** (no such file: a non-reference proteome) does the
script warn loudly, record the warning in its output, and fall back to the REST full set — it
never swaps databases silently. If REST also fails it exits with the reason.

### User-supplied target sequences (`--add-fasta <file>`, repeatable)

> Superseded by "6. Building the FASTA (step 6)" below.

A bait, a tag or a construct the experiment is about (EGFP, TurboID). Each entry is a TARGET,
written after the proteome and before the contaminants, and recorded in the sidecar under
`added_sequences` (the file's path and sha256; each entry's accession, name, length and sequence
sha256). A contaminant entry most of which is the added protein — identical to it, contained in
it, or with at least `ADDED_SEQUENCE_SHARED_FRACTION` (0.5) of its own peptides also in it, in the
search's digest (I = L) — is removed, from the contaminant set and from a supplied database's own
`Cont_` entries, and listed under `contaminants_dropped_for_added_sequences` (each with
`shared_fraction` and `shared_peptides`). Measured: wild-type GFP (`Cont_P42212`) shares 21 of its
28 peptides with EGFP (75%) and has 3 of its own, so the ordinary 2-own-peptide rule kept it
beside an EGFP bait. An entry sharing fewer stays — dropping it would also hide real
contamination with it, which its own peptides show — and the peptides it shares are named under
`contaminants_sharing_peptides_with_added_sequences`, in the Methods and in a report callout, as
ambiguous between the two. The digestion enzymes in use stay
(`contaminants_kept_near_added_sequences`). A supplied database's contaminant entries are judged
under either tag (`Cont_`, and FragPipe's `contam_`). The ordinary own-peptide rule counts a
contaminant's own peptides against **every** target -- the proteome AND the added sequences -- so
one whose peptides are split between a proteome protein and a bait (one of its own left in the
whole database) is removed. Refused: text before the file's first `>` header (it would be
written onto the previous entry's sequence), an entry with a contaminant tag, no sequence, or an
accession (`protein_ids.header_accession`) used twice or already in the database. An entry
identical to a proteome entry is warned. Methods and `reproduce.sh` carry them.

### Contaminants (`--contaminants`, default `universal`)

> Superseded by "6. Building the FASTA (step 6)" below.

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
intensity).

**So are entries the search cannot tell apart from a target protein.** Identity misses the
near-identical ones: bovine EF1A1 and 1433Z each differ from the MOUSE protein by one residue,
and in a 30-run mouse-brain search (an IP cohort, DIA-NN 2.6.1) both were `Cont_`-only groups in
every run while mouse Eef1a1 had no protein group at all. `fetch` therefore also digests every
contaminant and target entry with the search's own digest (`estimate_params.DIANN_DIGEST`:
`--cut K*,R*`, 1 missed cleavage, 7–30 aa, N-terminal Met excision; I and L equal) and drops a
contaminant that shares a peptide with a target but has fewer than `--min-unique-peptides`
(default 2) peptides of its own that do not overlap one another — missed-cleavage variants of
one residue difference count once. Each record carries `n_unique_peptides` (that count),
`n_unique_peptides_all`, `n_shared_peptides` and the target sharing the most peptides. Measured
2026-09-25 on the Universal set: human UP000005640 161 dropped (8 by peptide: KRT34, KRTAP10-4,
KRTAP4-7, TBA1D_BOVIN, TPM2_BOVIN, three sheep KRB2), mouse UP000000589 41 (10 by peptide:
EF1A1_BOVIN, 1433Z_BOVIN, TBA1D_BOVIN, four human KRTAPs, three sheep KRB2). An entry with
2+ own peptides stays a contaminant (bovine ENO1, LDHB, GSN, CAP1 against mouse). The sidecar
records `min_unique_peptides` and `contaminant_digest`; a sidecar with the rule but without
`min_unique_peptides` was built by the identity rule alone, and `reproduce.sh` replays it with
`--min-unique-peptides 0`. The digestion enzyme(s) actually used (`fetch --enzyme`, default `trypsin,lysc`)
are kept as contaminants even when they match a target protein — they are reagents — and
recorded under `contaminants_kept_despite_target_match`; any other protease entry follows the
normal rule (S. aureus's own SspA is identical to the Glu-C entry and stays quantified on a
trypsin digest). The
auditors flag the dropped proteins as "possible contamination, kept in quantification".

**A keratin sample loses every keratin entry (`fetch --keratin-sample`).** For hair, wool,
feather, skin or nail, keratin is the analyte, and the keratin-family entries the rules above
leave (for human: KRT34, KRTAPs, mouse hair and sheep wool keratins) take every peptide they
share with the sample's keratins out of quantification. With `--keratin-sample`, every
keratin-family `Cont_` entry (`fetch_fasta.is_keratin_gene`: gene KRT*/KRTAP*, or a protein
name starting "Keratin" — the Universal set's 14 sheep wool keratins have no gene name; 189
entries in all) leaves the set and a supplied database's own `Cont_` entries, recorded with
`keratin_sample: true` under `contaminants_dropped_keratin_sample`. `run_search.py
--keratin-sample` refuses a database that still holds one; `reproduce.sh` replays the flag.

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
citation, and any warnings. Pass it to `provenance.py --fasta-info-file` so
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

How many CPUs each array task of the DIA-NN chain asks for is sized to the chosen queue's
per-user cap (`references/diann_parallel.md`, "CPUs per array task"): 8 each on `high` for a
large cohort, so 8 files run at once under the 64-CPU cap, instead of `--threads` each.

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

## 5. Acquiring the pinned engine (step 5)

Honor the manifest's exact version — not "latest":

```
PIN_ENGINE=diann PIN_VERSION=2.7.0 bash scripts/acquire_tools.sh <platform_class>
```

Reads/writes `~/.proteomics-pipeline/tools/tools.json`. On HIVE it reuses the Core's
native DIA-NN build for that version (downloading it if the Core has none — a download runs
itself on a compute node, so run the call in the background); on mac it uses
Docker for DIA-NN. **Read `tools.json` `notes`** — license gates (FragPipe, Radiant) and
missing-runtime warnings surface there. `tools.json` `versions` is the build it found,
never `latest` and never the version that was ASKED for (`""` when it cannot tell, `env` for
a `sage` on PATH). `run_search.py` records it as the engine version only if it is shaped
like one — `2.6.1`, `v0.14.7` — so any other word there ends as `version: null`, not as a
version deposited to FRAN. If `run_search.py`
later prints an **engine version mismatch WARNING**, the version the user confirmed is not
the one about to run — say so before submitting (both land in `search_provenance.json`). If
it prints that the record is **`version: null`**, say that too: the search will be recorded
and deposited with no engine version.

Radiant is a ~3 GB image pull, so it is **not** acquired unless asked for:

```
ACQUIRE_RADIANT=1 PIN_ENGINE=radiant PIN_VERSION=2.3.3 bash scripts/acquire_tools.sh <platform_class>
```

It needs Docker or Apptainer — Seer ships no native binary. On HIVE the `.sif` must
already exist; building one pulls ~3 GB and **must run on a compute node, never the
login node** (the note in `tools.json` gives the exact line: `on_compute.py` with the queue
spelled out, and the pull's cache in the job's `$TMPDIR`, off the 20 GB home).

## 6. Building the FASTA (step 6)

Use the proteome, database type, and contaminant set the **user confirmed in
step 3**:

```
python3 scripts/fetch_fasta.py fetch --proteome <confirmed UPID> \
    --content <one_per_gene|reviewed|reviewed_isoforms|full|full_isoforms> \
    --contaminants <universal|cell_culture|...|none> \
    [--enzyme trypsin,lysc] --keratin-sample|--no-keratin-sample --out ./search.fasta [--hive]
```

Pass the step-3 answer every time: `--no-keratin-sample` when the samples are not keratin. With
neither, the sidecar records `keratin_sample_source: default` and every downstream step tags the
"not keratin" as an assumption (`run_de.R` writes one `DEFAULT — not user-confirmed` methods line
when keratin precursors were removed on it).
`--keratin-sample` for a keratin sample (step 3): every keratin-family `Cont_` entry (gene
KRT*/KRTAP*, or a protein name starting "Keratin" — the sheep wool keratins have no gene
name) leaves the contaminant set; a `--path` database that already carries `Cont_` entries loses
its keratin ones too (still build with `--proteome`, as below).
The sidecar records `keratin_sample: true`, the removed entries under
`contaminants_dropped_keratin_sample` (each with its `source`) and
`keratin_contaminants_in_database: []`. Every other build records `keratin_sample: false` and
lists the keratin-family `Cont_` accessions its database still holds.
`--enzyme` names the digestion enzyme(s) actually used (trypsin, lysc, gluc, chymotrypsin,
aspn, argc, lysn, pepsin; comma- or slash-separated — `Trypsin/P` and `Trypsin/Lys-C` work). The default `trypsin,lysc` is the Core's usual
mix — leave it unless the user says otherwise; change it for a different digest (`--enzyme
gluc`, `--enzyme trypsin,pepsin`). Only those enzymes' contaminant entries stay `Cont_` when
they match a protein of the searched organism (their autolysis peptides stay out of quant);
any other protease is treated as sample protein — which matters for organisms that make the
enzyme themselves (S. aureus's own SspA is the Glu-C entry, quantified on a trypsin digest).
Defaults are `--content one_per_gene` (the canonical set, from UniProt's
reference-proteome FTP tree) and `--contaminants universal`. **Do not pass
`--content full` casually** — for human that is 147,506 sequences instead of
20,652, a 7× larger search space. UniProt's REST `onePerGene` parameter is
silently ignored, which is why the canonical set must come from FTP.

Pass `--hive` when `uc_davis_hive` is true to reuse pre-staged proteomes +
contaminants under `/quobyte/proteomics-grp/MRS/` instead of downloading.
A staged copy's UniProt release is unknown (it was downloaded earlier), so
`uniprot_release` stays empty; the sidecar records the copy's date and sha256 under
`staged_file`, fills organism/taxid from UniProt (or, offline, from the filename taxid —
see `organism_source`), and reports `content_inferred`, a one-per-gene check of the entry
count against UniProt's geneCount — an inference, not a verification. `--content` is not
applied to a staged file; if it matters, omit `--hive`.

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

**A database you supply (`--path`) needs its organism named** — nothing else reports it, and
Methods, the report and FRAN read it from the sidecar: `--organism '<scientific name>' --taxid
<taxid>` (the step-3 answer; `--ncbi-organism/--ncbi-taxid` work too), or `--organism none` for
a database with no single organism (several species, synthetic). Without them `fetch` stops and
asks. With `--proteome` the organism is UniProt's, and an organism flag there is refused.

**The user's own sequences — a bait (EGFP, a BirA/TurboID construct), a tag — go in with
`--add-fasta <file>`** (repeatable), on top of the usual `--proteome` build; never by
concatenating files by hand, which bypasses the contaminant rules and loses the organism and
release. They are TARGETS, written after the proteome. A contaminant entry **most of which is an
added protein** — identical, contained, or with at least half of its own peptides in one
(`fetch_fasta.ADDED_SEQUENCE_SHARED_FRACTION` = 0.5, the search's digest) — is removed, so the
added sequence keeps those peptides: the Universal set's wild-type GFP (`Cont_P42212`) shares 21
of its 28 peptides with EGFP (75%), and the ordinary 2-own-peptide rule had kept it. An entry
sharing fewer **stays** — its own peptides are the evidence of real contamination with it — and
the peptides it shares are flagged **by name** as ambiguous between the two
(`contaminants_sharing_peptides_with_added_sequences`, the Methods and a report callout). The
digestion enzymes in use stay. The sidecar records each file (path, sha256) and entry
(`added_sequences`) and the removed entries (`contaminants_dropped_for_added_sequences`);
`reproduce.sh` replays `--add-fasta`. An entry with a `Cont_` tag, an accession used twice or
already in the proteome is refused. Tell the user which contaminant entries left, and which
peptides of the added sequences are ambiguous.
**Then the bait, once the database is built** (an enrichment with added sequences): put the
FASTA sidecar in the session (`<session>/input/search.fasta.meta.json`), run
`python3 scripts/experiment_type.py propose --session <session>`, and when its `ask_bait` is not
null (it is null when the bait is already answered) ask it in its own words and record the answer:

```
python3 scripts/experiment_type.py set-bait --session <session> --bait <accession>   # the bait
python3 scripts/experiment_type.py set-bait --session <session> --bait none          # none is
```

It is kept as `bait_confirmed` (the answer, who, when) in `input/experiment_type.json`, and a
re-recorded type keeps it. Without an answer the data check would take the ONE added sequence as
the bait, unconfirmed -- never leave it there. More than 10 added sequences (a whole proteome) are
counted, not listed (`n_added_sequences`): ask for the bait by name, and `set-bait` takes its
accession, entry name or gene.

**Read the returned JSON** (also written to `<out>.meta.json` for the repro
bundle) and act on it:

- `warnings` non-empty → tell the user before searching; a one-per-gene→full
  fallback changes the database out from under them.
- `diann_cont_quant_exclude` → pass as `--cont-quant-exclude Cont_` to DIA-NN in
  step 7 so contaminants are identified but kept out of DIA-NN's normalisation. That flag
  does NOT reach the DE, which re-quantifies from the report: `run_de.R` removes them
  itself (step 8).
- `n_contaminants_dropped_as_target` > 0 is normal, not a failure. The universal
  contaminant set holds sequences identical to real proteins (human keratins; bovine
  ACTB/EEF1A1/tubulins that are residue-for-residue the human proteins) and near-identical to
  others (bovine EEF1A1 and YWHAZ each differ from the MOUSE protein by one residue). Left in,
  DIA-NN reports those proteins ONLY as `Cont_` groups, which `--cont-quant-exclude` and
  `run_de.R`'s contaminant filter (step 8) then drop from quant — ACTB, EEF1A1 and KRT8
  vanished from a real HeLa search this way, and EEF1A1 and
  YWHAZ from every run of a mouse-brain study. `fetch` removes every contaminant entry
  identical to (or contained in) a target entry, and every entry that shares peptides with a
  target but has fewer than 2 non-overlapping peptides of its own in the search's digest (I =
  L; `--min-unique-peptides`, 0 = identity rule only), lists them under
  `contaminants_dropped_as_target` with each one's `reason`, target and `n_unique_peptides`,
  and warns (expect ~161 for human: 152 identical + 1 contained + 8 by peptide; ~41 for mouse:
  31 + 10 — bovine EEF1A1, YWHAZ, TUBA1D and seven hair KRTAPs);
  the digestion enzyme(s) named by `--enzyme` are the exception and stay contaminants
  (`contaminants_kept_despite_target_match`; a warning naming a protease dropped as "not in
  --enzyme" means: if that enzyme WAS used, re-run `fetch` with it in `--enzyme`). Tell the user those proteins — skin keratins
  included — are now quantified and count toward normalisation; step 8c flags them.
  An entry with 2+ peptides of its own **stays**, even when the sample shows only the peptides it
  shares with a target (horse myoglobin `Cont_P68082` has 9 of its own; in human heart its group
  rested on human MB's peptides, one differing only by I/L): which peptides a sample shows only
  the search can say, so step 8c's `contaminant_evidence` flags such a group. `--drop-contaminant
  <accession>` (repeatable) then leaves that entry out of a rebuild: recorded in the sidecar, named
  in the Methods, replayed by `reproduce.sh`; an accession that matches no contaminant entry stops
  the build.
- A `--path` database that already contains `Cont_` entries cannot be fixed: `fetch` warns
  and lists them under `contaminants_identical_to_target_kept` — rebuild it with
  `--proteome` instead. The Core's Sep-2025 `MRS/UP000005640_9606_plus_universal_contam.fasta`
  (kept only for old searches' provenance) has this problem: do not pass it with `--path`.
  Its rebuild, `MRS/UP000005640_9606_plus_universal_contam_2026-09.fasta`, was made by this
  `fetch` (sidecar beside it, state `current`) with a DIA-NN 2.7.0
  `..._2026-09.predicted.speclib`, for searches run outside the skill. Inside the skill keep
  building the database with `fetch --proteome`.
- Contaminants that can't be fetched are a **hard stop**, not a warning. Fix the
  source or have the user explicitly choose `--contaminants none`.
