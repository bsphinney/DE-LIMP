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
  as **native builds** at `build_<nnn>/diann-<version>/diann-linux` — **2.5.1, 2.6.0 and
  2.6.1** when listed on 2026-09-16 — plus one older Apptainer image, `diann_2.3.0.sif`.
  The skill pins **2.6.1** (`resolve_defaults.py`). `acquire_tools.sh` resolves the
  **pinned** version by looking for `build_*/diann-<version>/diann-linux` first, then a
  version-matched `.sif`, and **never silently substitutes a different version**
  (reproducibility): a pin that is not there — the FRAN pilot's **2.7.0**, say — is
  downloaded as the Academia Linux build into your own tools root instead. Unpinned, it
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
for the v0.14.7 release; for Radiant from a pinned image tag or the `.sif` name. When no
version can be determined it records `""` with a note — an unpinned Radiant image
(`seerbio/radiant-fulcrum:latest`) is one such case. Two values are sources, not builds:
`"env"` (a `sage` found on PATH, e.g. the conda env's) and, in a `tools.json` written
before this was fixed, `"latest"`.

The manifest and `tools.json` can still disagree — `resolve_defaults.py` pins 2.6.1, but
`tools.json` holds whatever was last acquired (the FRAN pilot pinned 2.7.0). Nothing forces
them to match, so `run_search.py` compares them before submitting and records both in
`search_provenance.json` — as one `engine_version` object shaped like `scan_window` (a
`value`, and a `source` saying in words where it came from), plus top-level `version`, which
always equals `engine_version.value`:

| field | meaning |
|---|---|
| `version` = `engine_version.value` | the build that runs: `tools.json` `versions`, else the one version the **command itself** names (a build folder such as `diann-2.6.0/`, a `.sif` name, an image tag); else `null`. **Never the manifest's pin** — `fran_deposit.py` sends `version` to FRAN as the engine version it stores, so a pin nothing confirmed would be stored there as fact |
| `engine_version.source` | where `value` came from: `tools.json versions.<engine>, …`, `named by the command tools.json runs (…)`, or `unknown -- <why>` |
| `engine_version.tools_json` | `tools.json` `versions.<engine>` as written, even `latest` / `env` |
| `engine_version.named_by_command` | every version the command names (usually zero or one) |
| `engine_version.manifest_pin` | the manifest's pin as written (`null` if it pins another engine) |
| `engine_version.mismatch` | `true`/`false` when `value` and the pin are both known; **`null`** when they cannot be compared. A WARNING is printed when `true` |

A `tools.json` whose `versions` contradicts its own command (says 2.6.1, runs
`diann-2.6.0/diann-linux`) records `version: null` with a WARNING: one is wrong and nothing
tells which. A version *directory* (`sage/0.14.6/sage`) is not read, because
`acquire_tools.sh` names cache folders after the request, not the build inside.

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

1. an explicit `--partition` **and** `--account`, used as given;
2. `genome-center-grp` on `high` (facility members: not preemptible);
3. `publicgrp` on `low` (everyone else: preemptible, so `#SBATCH --requeue` is added).

When an account has both, **utilisation** decides, not entitlement. `slurm_queue()` counts
the CPUs *you* already run on `high` (`squeue`) against the 64-CPU per-user cap, and the idle
CPUs on `low` (`sinfo`), and compares both with `need = min(peak_cpus, 16)` — one array
task's worth, **not** the job's own request (a caller that passes no peak, like the
single-script `--sbatch` paths, gets 16):

- **any job** goes to `low` when fewer than `need` of your CPUs are free on `high` and `low`
  has at least `need` idle;
- the 5-step chain's **array steps (2 and 4)** also go to `low` (given the same `need` idle
  there) when fewer than `2 × need` are free on `high`, or when your usage there cannot be
  read (`squeue` cannot run) — a preempted task costs one file;
- **steps 1, 3 and 5** cannot restart mid-way, so only the first rule applies. Because
  `need` is capped at 16, a 64-CPU step 3 **stays on `high` while 16 or more of your CPUs
  are free there, then waits on `high`** for the rest; it does not move to `low` because
  the 64 it asked for are unavailable. To send it to `low`, pass `--partition low --account
  publicgrp` (plus `--qos` if the association has one).

If associations cannot be read at all the script falls back to `publicgrp/low`, **not** the
cluster default: that is `high`, which rejects a non-facility account.

Checked on HIVE 2026-09-16 (`sacctmgr show assoc` / `show qos`, `sinfo`): the default
partition is `high`; `genome-center-grp` has `high` (per-user cap 64 CPUs) and `gpu-a100`
but no `low`; `publicgrp` has `high` (8 CPUs / 128 GB **per job**) and `low` (no per-job
cap, preemptible).

`--partition/--account/--qos` are forwarded to the job-array routes (the DIA-NN 5-step
chain and the Radiant per-file array). The single-script `--sbatch` paths — DIA-NN at ≤5
files or when the chain declines, Sage, FragPipe, AlphaDIA, single-file Radiant — use the
detected queue; read the `#SBATCH` header they print before submitting. Submit with
`sbatch job.sh` (or `bash <out>/submit.sh` for a two-job or chained search), poll the
`<job>_<id>.log`, then run `run_search.py --adapt-only` for Sage/FragPipe to build
`report.parquet`.

A DIA-NN search of more than 5 files on a SLURM host routes to the 5-step parallel chain
instead, which has no single job script: `job.sh` is not written, an existing regular file of
that name is renamed to `job.sh.stale-<time>` once the chain is generated, and
`run_search.py` exits 3. Submit `<out>/submit.sh`
(→ `references/diann_parallel.md`).
