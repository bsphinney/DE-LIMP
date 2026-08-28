# Search engines reference

Acquisition detection, engine routing, and normalizing each engine's output to the
DE-input contract.

## Acquisition detection (`detect_acquisition.py`)

| format | how DIA/DDA is decided | instrument |
|---|---|---|
| Bruker `.d` | `analysis.tdf` SQLite: `DiaFrameMsMsInfo`/`DiaFrameMsMsWindowGroups` or `Frames.MsMsType==9` → DIA; `PasefFrameMsMsInfo`/`MsMsType==8` → DDA | `GlobalMetadata.InstrumentName` (e.g. "timsTOF Pro") |
| `.mzML[.gz]` | stream MS2 isolation windows: median width ≥3 Da over few centers → DIA; ≤2 Da, many centers → DDA | none (mzML rarely carries model reliably) |
| Thermo `.raw` | needs ThermoRawFileParser → mzML, else `unknown` | ThermoRawFileParser `metadata` model string, if present |
| `.wiff` | convert to mzML first | none |

Confidence is `high`/`medium`/`low`. **Anything not `high`, mixed, or with
disagreeing instruments sets `needs_confirmation` — ask the user.** Instrument may
be null; that's fine (matcher falls back to score 0 and the user confirms).

## Parameters

Search parameters are derived from the data type (instrument + acquisition) by
`estimate_params.py`, not hand-maintained per workflow. → `references/parameters.md`.

## Routing

Default by acquisition: **DIA → DIA-NN, DDA → Sage.** `--engine` overrides. The
bundle's `engine.name` is authoritative when present. FragPipe is opt-in only (the
bundle names it or the user asks) because MSFragger/IonQuant are license-gated.

## Licensing — pick a commercial-OK engine for non-academic users
| Engine | License | Commercial use | Notes |
|---|---|---|---|
| **DIA-NN** (Academia build) | academic / non-profit only | ❌ **no** | the free build the skill downloads is academic-only |
| **AlphaDIA** | Apache-2.0 | ✅ yes | open-source DIA alternative; deep-learning, **GPU recommended** |
| **Sage** | MIT | ✅ yes | fast, CPU; DDA + wide-window DIA |
| **FragPipe** (MSFragger/IonQuant) | license-gated | depends | requires the user's own license |

**For commercial / non-academic DIA users, route to AlphaDIA** (`--engine alphadia`)
instead of DIA-NN — same job, no license problem. Ask "academic or commercial?" before
a DIA run if it isn't already clear.

## Public program sources (anyone, any OS — cite these for non-Core / external users)

This skill is used **outside** UC Davis and the Core. Never assume the internal
`/quobyte/proteomics-grp` copy of a tool exists — **always point the user to the public
source** and a way to run it on *their* OS (macOS, Windows via WSL2 **or** native
Windows, or Linux). `setup.sh`/`acquire_tools.sh` fetch the free ones automatically.

| Program | Public source | Platforms | Install / note |
|---|---|---|---|
| **DIA-NN** (Academia) | https://github.com/vdemichev/DiaNN/releases | Windows, Linux | Windows GUI+CLI, or Linux binary; academic / non-profit only. `acquire_tools.sh` fetches the Linux build. DDA since 2.3 (`--dda`). |
| **Sage** | https://github.com/lazear/sage/releases | Win, macOS, Linux | MIT; single static binary. |
| **AlphaDIA** | https://github.com/MannLabs/alphadia | Win, macOS, Linux | Apache-2.0 (commercial-OK); `pip install alphadia`, GPU recommended. |
| **Radiant DIA + Fulcrum** (Seer) | https://github.com/seerbio/radiant-fulcrum-container · image `seerbio/radiant-fulcrum` | container only (multi-arch amd64+arm64) | **Thermo Orbitrap only** here (mzML/Parquet input). Needs a DIA-NN-generated library. **Apache-2.0 + Commons Clause + grant-back** — restricts selling a derived service. |
| **FragPipe** (MSFragger/IonQuant) | https://github.com/Nesvilab/FragPipe/releases | Win, macOS, Linux | Java GUI; MSFragger/IonQuant need the user's own (free-academic) license. |
| **ThermoRawFileParser** | https://github.com/compomics/ThermoRawFileParser/releases | Win, macOS, Linux | `.raw`→mzML; cross-platform .NET. |
| **ProteoWizard / msconvert** | https://proteowizard.sourceforge.io/ | Windows (native); Linux/macOS via Docker | vendor→mzML; Linux via the `chambm/pwiz-...` Docker image. |
| **.NET 8 runtime** | https://dotnet.microsoft.com/download/dotnet/8.0 · installer script https://dot.net/v1/dotnet-install.sh | Win, macOS, Linux | needed only so the **Linux** DIA-NN binary can read `.raw`; `ensure_dotnet8.sh` installs 8.0.latest. **Native Windows doesn't need it.** |
| **UniProt proteomes (FASTA)** | https://www.uniprot.org/proteomes/ · canonical sets: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ | web / REST / FTP | `fetch_fasta.py resolve` finds the ID from an organism name; `fetch` downloads it. One-per-gene comes from the **FTP** tree — REST `onePerGene` is silently ignored. |
| **Contaminants (Hao lab, `Cont_`-tagged)** | https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics | web | `fetch_fasta.py --contaminants <set>`. Universal + sample-type sets; JPR 2022, doi:10.1021/acs.jproteome.2c00145. (The old GPM cRAP URL now 404s — do not use it.) |

**Rule:** whenever you tell a user to run a program, give the public link **+** the exact
command **+** an OS-appropriate path — not a HIVE/`/quobyte` path they can't reach.

## Per-engine invocation & output adapter (`run_search.py`)

### DIA-NN (native contract — no adapter)
```
<cmd> --cfg <bundle diann.cfg> --f <file> [--f <file> ...] \
      --fasta <fasta> --out report.parquet --threads N [--dda]
```
`report.parquet` is already the DE contract — `run_de.R` reads it directly.

**DDA data (DIA-NN 2.6+).** DIA-NN 2.3+ can search DDA / DDA-PASEF with **`--dda`**
(`run_search.py` adds it automatically when the bundle's acquisition is `DDA`). Notes:
- `--dda` **must** be used with DDA and **must not** be used with DIA data.
- QuantUMS is auto-disabled on DDA; PTM-localisation probabilities are unreliable on DDA.
- For DDA **quant**, DIA-NN recommends extra MS1 filtering on `Ms1.Global.Q.Value`
  (< 0.0001–0.01) and `Ms1.Global.Quality` (> 0.5–0.9), optionally `Ms1.Q.Value` /
  `Averagine`. Standard DIA/DDA q-value filters still apply.
- It's officially "beta", but performs strongly — on UC Davis nail (Exploris, 67 runs)
  it matched/beat the delivered FragPipe/Scaffold result at ~equal keratin coverage.
- Default routing still sends DDA → Sage (validated); use `--engine diann` (or a
  DDA+DIA-NN workflow) to search DDA with DIA-NN.

**Reading Thermo `.raw` directly (the .NET 8 gotcha — Linux only).** On **native
Windows** DIA-NN reads `.raw` out of the box (skip this whole section). On **Linux**
(incl. WSL2 and HIVE) the 2.6 *native* binary reads `.raw` via a bundled .NET
(ThermoFisher.CommonCore) component and needs a **.NET 8 runtime ≥ 8.0.17** on PATH. Without it, DIA-NN prints
`ERROR: cannot read .raw files ... .NET Runtime 8: 8.0.17 or later` and processes **0
files**. It does a HARD ≥ 8.0.17 check: an older 8.0.x (e.g. a cluster's `8.0.4`
module) is **rejected**, and a .NET 9 runtime is **not** used (rollForward is
LatestMinor — it won't cross the major version). Fix, done for you by
**`scripts/ensure_dotnet8.sh`** (idempotent; run it on a login node — it needs
internet, compute nodes usually don't):
```
export DOTNET_ROOT="$(bash scripts/ensure_dotnet8.sh | tail -1)"   # installs 8.0.latest if missing
export PATH="$DOTNET_ROOT:$PATH"
```
`run_search.py` calls this automatically whenever DIA-NN inputs include `.raw` (inline
run, and baked into the emitted sbatch as a preamble). When it's right, DIA-NN logs
`.NET runtime found, Thermo .raw support enabled`. Alternative with **no** .NET: feed
`.mzML` (DIA-NN reads those natively) — convert `.raw`→mzML with ThermoRawFileParser /
msconvert if you don't already have them.

### Radiant DIA + Fulcrum (Seer; Thermo Orbitrap only)

Peptide-centric DIA search (Radiant) inside Seer's Fulcrum workflow engine, which
adds rescoring, FDR, protein inference and directLFQ rollup. Third of the three DIA
routes; **default to DIA-NN** unless the user asks for this one.

```
python3 scripts/run_search.py ... --engine radiant --params ./wf/default.radiantConfig
```

Everything below was verified against `seerbio/radiant-fulcrum:2.3.3` and the
Radiant source (`seerbio/radiant`, internal name **Pythia**), not inferred from docs.

**Constraints that decide whether you can use it at all**

- **Thermo Orbitrap only, in this route.** The container's search backend accepts
  **mzML or Parquet** (`radiant_fulcrum_search/search.py`), so `.raw` is converted
  with `msconvert` first and Bruker `.d` is refused with a pointer to the DIA-NN or
  diaTracer route. (Radiant's *source* does have an ion-mobility module, but that
  is not reachable through this container CLI.)
- **A spectral library is always required.** `--library` is `required=True`, and
  **`--libfree` does not mean library-free** — in the click definition it is the
  same switch as `--no-mbr` (`"--mbr/--no-mbr", "--no-libfree/--libfree"`), so it
  selects single-pass vs MBR. So **this route needs DIA-NN acquired too**.
- **⚠ DIA-NN and Radiant share no library format directly** — verified on real runs,
  and the reason `make_radiant_library.py` exists:
  1. `--out-lib foo.tsv --predictor` does **not** write `foo.tsv`. DIA-NN ignores the
     extension for predicted libraries and writes `foo.predicted.speclib`.
  2. Radiant *does* read `.speclib`, but only **format v ≥ −3**
     (`SpecLibSrc/Library.h`). DIA-NN 2.6 writes **v−11**, 2.x generally v−10, so
     `RadiantDIA` aborts with `ERROR: version is not supported11`. Every speclib on
     the UCD share is v−10/v−11.
  3. DIA-NN **cannot emit TSV at all** — asked to convert to `.tsv` it silently
     writes `.parquet`.
  4. That parquet uses **dot-separated** names (`Precursor.Mz`, `Relative.Intensity`)
     while Radiant wants concatenated ones (`PrecursorMz`, `LibraryIntensity`), so a
     straight dump is unreadable.

  `make_radiant_library.py` does predict → parquet → renamed TSV, and
  `run_search.py` calls it automatically. A human proteome library is ~50 M rows /
  ~10 GB of TSV.
- **Multi-file runs are serial in Radiant itself** — the backend raises
  `NotImplementedError` for parallel mode. On a cluster the skill works around this
  (below), so this only bites on a single machine.
- **Container only.** No native binary on any platform; multi-arch image
  (linux/amd64 + linux/arm64) so it runs natively on Apple Silicon and on HIVE.
  `acquire_tools.sh` records the runtime (`docker`/`apptainer`) and image separately
  because `run_search.py` has to inject bind mounts, and `-v` vs `--bind` differ.
  Not acquired by default (~3 GB) — set `ACQUIRE_RADIANT=1` or pin the engine.

**⚠ Licence.** Apache-2.0 as modified by a **mandatory grant-back** and the
**Commons Clause**. The Commons Clause removes the right to *Sell* — including
selling a service whose value derives substantially from the software. That is a
live question for a fee-for-service core; surface it before running and let the
institution decide. It also appears in `tools.json` `notes`.

**Cluster parallelism (`radiant_parallel.py`).** The serial limit is in Radiant's
*orchestration*, not in the science: each file is searched independently, and only
the downstream stages need every file at once. So the chain splits exactly there —
and `run_search.py` routes to it automatically when there is >1 file and `sbatch` is
on PATH:

| step | job | what it does |
|---|---|---|
| 1 | single | DIA-NN predicts the spectral library |
| 2 | **array, N tasks** | one `RadiantDIA` per mzML → `<results>/radiant-results/<name>.radiantDIA` |
| 3 | single | `fulcrum --toml-file` rescoring + FDR + inference + rollup over all files |

Step 3 does **not** re-search: `_execute_radiant()` resolves its result path as
`<output_location>/<input name>.radiantDIA` and returns it untouched when
`reuse_existing` is set, so the generated TOML carries `[overrides.search]
reuse_existing = true` and Fulcrum goes straight to the downstream stages. That
directory name is the contract between steps 2 and 3 — don't rename it.

Step 2 skips any file whose `.radiantDIA` already exists, so a preempted
`publicgrp/low` task requeues without redoing work, and step 3 refuses to rescore a
partial set. `--no-parallel` forces the single serial search.

**Feeding results to the DE-LIMP app.** `adapt_radiant` produces the skill's own DE
contract (one row per protein×run, `PG.MaxLFQ`). DE-LIMP wants something different: its
upload feeds `limpa::readDIANN()`, which is **precursor-level** and needs `Run`,
`Precursor.Id`, `Precursor.Normalised`, `Protein.Group`, `Protein.Names`, `Genes`,
`Proteotypic`, `Q.Value`, `Lib.Q.Value`, `Lib.PG.Q.Value`. Use
`radiant_to_delimp.py --results <out>/radiant_results --library <lib.tsv> --out
delimp_report.parquet` and upload that. It derives `Precursor.Id` from modified
sequence + charge, maps the `Global.*` q-values onto the `Lib.*` names, and recovers
`Protein.Names`/`Genes` from the library — Fulcrum drops them, and without the join
DE-LIMP loses gene-level annotation (GSEA, gene labels). Pass `--library` or those
columns are emitted empty rather than invented.

**Run names.** Fulcrum reports `Run` as a URI of the per-file result
(`file:///…/Sample.mzML.radiantDIA`). Both adapters reduce it to the bare `Sample`, so
a `conditions.csv` keyed on sample names matches — a single `splitext` would leave
`Sample.mzML` and silently fail group assignment.

**Adapter.** Fulcrum's `combined` output backend already uses DIA-NN-style column
names (`Run`, `Protein.Group`, `Precursor.Quantity`, `PG.Quantity`, `PG.Normalised`,
`Q.Value`, `Global.PG.Q.Value`, …) — but writes them as a **spark parquet directory**
of part-files at `<results>/fulcrum-results`, not a single file. `adapt_radiant`
reads it as a dataset, collapses PSM-level rows to one protein×run (max intensity,
best q-value), and writes the usual `report.parquet`.

### AlphaDIA (commercial-OK DIA; adapter required)
Apache-2.0 — the open-source alternative to DIA-NN for non-academic users. Library-free:
```
<cmd> -o <out> -f <raw> [-f ...] --fasta <fasta> [-c <config.yaml>]
```
- Install: `pip install alphadia` (into the conda env); `acquire_tools.sh` does this
  when AlphaDIA is pinned/requested. **GPU strongly recommended** (CPU works but is slow)
  — best on a HIVE GPU node. Reads `.raw`/`.d`/mzML; no msconvert needed.
- **Adapter:** AlphaDIA writes `pg.matrix.parquet` (protein-group × run matrix; also
  `precursors.parquet` with `raw.name`/`pg.name`/`pg.intensity`). `adapt_alphadia` melts
  the matrix → DIA-NN-shaped `report.parquet` (Run, Protein.Group, PG.MaxLFQ; Q-values
  zeroed since AlphaDIA already FDR-filtered). Like the Sage adapter, confirm it on real
  data the first time.
- **⚠ Limitation — timsTOF whole-proteome directDIA:** AlphaDIA does **not** work on
  **timsTOF** `.d` for **whole-proteome library-free (directDIA)** searches. For that
  case it is not a valid substitute for DIA-NN. Options: academic → DIA-NN; commercial →
  run AlphaDIA **library-based** (a predicted/empirical spectral library, `--library`)
  instead of directDIA, or use Sage. Non-timsTOF instruments and library-based timsTOF
  runs are unaffected.

### Sage (mzML-first; adapter required)
1. Convert `.d`/`.raw` → mzML with `msconvert` if needed (fails loudly if msconvert
   is absent and inputs aren't mzML).
2. `<cmd> <bundle sage_config.json> -f <fasta> -o <out> --parquet
   --disable-telemetry-i-dont-want-to-improve-sage <mzml...>`
3. **Adapter:** map `lfq.parquet` (protein, filename, intensity) →
   DIA-NN-shaped `report.parquet` with `Run, Protein.Group, PG.MaxLFQ` and zeroed
   Q-value columns (Sage already FDR-filtered). Q-values default to 0 so the
   MaxLFQ DE path keeps every row.
   - v0.14.x has no native protein grouping → may need `sage_protein_groups.py`
     post-hoc (DE-LIMP keeps it on HIVE). v0.15+ has IDPicker grouping.

### FragPipe + diaTracer — spectrum-centric DIA (timsTOF dia-PASEF)

Verified against FragPipe 24.0 source, the diaTracer site, and the paper. Where a claim
below is load-bearing the file that proves it is named.

**What it is.** diaTracer does 3-D (m/z, RT, ion mobility) peak tracing on dia-PASEF `.d`
data and emits **pseudo-MS/MS spectra as mzML**, which MSFragger then searches as if they
were DDA. That builds a spectral library from the data itself; DIA-NN then quantifies
against the **original `.d`**. Because identification is an ordinary database search, this
route reaches what library-free DIA-NN does poorly: semi-tryptic (N-terminomics),
nonspecific (HLA), and open / mass-offset PTM searches.

Full chain: `.d` → **diaTracer** → MSFragger (DDA mode) → MSBooster → Percolator →
Philosopher → EasyPQP (library) → **DIA-NN** (quant vs the original `.d`) → IonQuant.

Cite: K. Li et al., *diaTracer enables spectrum-centric analysis of diaPASEF proteomics
data*, Nat Commun 16, 95 (2025).

| Thing | Value | Verified from |
|---|---|---|
| FragPipe | 24.0 (2025-12-24); diaTracer arrived in 22.0 | releases API |
| Requires | MSFragger 4.4+, IonQuant 1.11.18+, diaTracer 2.2.1+, Java 11+, Python 3.9–3.11 | 24.0 release notes |
| Workflow preset | `DIA_SpecLib_Quant_diaPASEF.workflow` (`diatracer.run-diatracer=true`) | FragPipe 24.0 `workflows/` |
| Manifest data type | **`DIA`** | `InputLcmsFile.java` |
| Output dir | `<workdir>/dia-quant-output/` | `CmdDiann.java` |
| Bundled DIA-NN | 1.8.2 Beta 8 | `tools/diann/` |

Other presets that enable diaTracer: `DIA_SpecLib_Quant_Phospho_diaPASEF`,
`Nonspecific-HLA-diaPASEF`, `chemprot-ABPP-IADTB-diaPASEF`. Plain `DIA_SpecLib_Quant` and
`DIA_DIA-Umpire_SpecLib_Quant` have it **false** — those are the Orbitrap/mzML paths, so
picking one of those for `.d` data silently gives you a different method.

#### Output: expect `report.tsv`, NOT parquet

FragPipe branches on the DIA-NN version (`CmdDiann.java`):

- **Bundled DIA-NN 1.8.2 Beta 8 (the default) → `dia-quant-output/report.tsv` only. There is
  no `report.parquet`.** Also produces `MSstats.csv` (gated on DIA-NN < 2.0).
- Only if you point `--config-diann` at your own **DIA-NN 2.x** does it write
  `report.parquet` (and convert it to TSV as well).

So `adapt_fragpipe_dia()` reads whichever exists, preferring parquet and otherwise
converting the TSV. In the stock configuration the TSV path is the one that runs.
`--matrices` also yields `report.pg_matrix.tsv` / `report.pr_matrix.tsv` (MaxLFQ-normalised,
1% FDR, samples as columns) — the paper's own DE route uses `report.pg_matrix.tsv`.

`combined_protein.tsv` **also** appears, because IonQuant runs on the pseudo-MS/MS mzML
(typed `DDA`). **Do not use it as the quant of record**: it is derived from pseudo-spectra
rather than the original `.d` chromatograms, is not equivalent to the DIA-NN numbers, and
the paper never uses it. This is why `adapt_fragpipe()` tries the DIA route *first*.

#### Gotchas that will bite on a cluster

- **The mzML is written next to the `.d` file, not into `--workdir`** (`CmdDiaTracer.java`
  rewrites the input path). This is handled for you — see "Symlink staging" below — but it
  is the reason that machinery exists.
- **Never put `.d` in a parent directory name.** FragPipe builds the output name with
  Java's `String.replace`, which replaces *every* occurrence, so `/data/proj.d/run.d`
  becomes `/data/proj_diatracer.mzML/run_diatracer.mzML`.
- **A bad manifest data type does not error — it silently becomes `DDA`.** Accepted values
  are `DIA`, `GPF-DIA`, `DIA-Quant`, `DIA-Lib`, `DDA+`, `DDA`; anything else falls through
  the `default:` branch to DDA and your DIA run quietly becomes a DDA run. Use exactly `DIA`.
- **`DIA-Quant` is quantify-only** — it does not trigger diaTracer and is excluded from
  library generation. It is for reusing an existing library, not for a fresh run.
- diaTracer only fires when the type is `DIA`/`GPF-DIA`/`DIA-Lib` **and** the file ends in
  `.d`; FragPipe hard-errors if diaTracer is on and a DIA input is not `.d`.
- **First headless run needs all three config flags**: `--config-tools-folder`,
  `--config-diann`, `--config-python`. Since 22.0 MSFragger, IonQuant **and diaTracer must
  live in the same folder**. diaTracer also needs the native Bruker libraries, which ship
  in MSFragger's zip under `ext/`.
- 24.0 fixed a bug where the **log file was not saved in headless mode** — use 24.0, not
  23.x, for unattended cluster runs.
- On HIVE, an Apptainer image exists: `apptainer pull docker://fcyucn/fragpipe:latest`,
  then `apptainer shell --compat --bind /quobyte:/quobyte fragpipe_latest.sif`. Both
  `--compat` and an explicit `--bind` are required.

#### Symlink staging + reuse (`diatracer_stage.py`) — automatic

Two problems, one mechanism. `run_search.py` calls this itself for DIA + `.d` inputs.

**Problem 1 — concurrent users collide.** diaTracer derives its output path from the input
path (`CmdDiaTracer.java:138`):

```java
f.getPath().toAbsolutePath().normalize().toString().replace(".d", "_diatracer.mzML")
```

So two people searching the same shared dataset both write `<name>_diatracer.mzML` into the
same folder and race each other; a read-only share fails outright. That is the normal case
for a class working from one copy of the data.

**The fix rests on one detail: `normalize()` is not `toRealPath()`.** It collapses `.` and
`..` but does **not** resolve symlinks. So a manifest pointing at a *symlink* to the `.d`
makes diaTracer write next to the **symlink**. `diatracer_stage.py` gives each run its own
directory of symlinks (names keep the `.d` suffix — `CmdDiaTracer` requires
`endsWith(".d")`), so every user gets their own pseudo-spectra and the shared raw directory
is never written to and may be read-only.

**Problem 2 — reconversion is expensive (~20 min/file).** The stager looks for an existing
mzML both in the staging dir and beside the real `.d` (so a facility pre-conversion counts),
matching `_diatracer.mzML` / `.diatracer.mzML` case-insensitively since the spelling has
varied across versions. When one exists it emits the preset's own reuse form:

| Situation | Manifest rows |
|---|---|
| fresh | `<stage>/<name>.d` → `DIA` |
| already converted | `<mzML>` → `DDA`, **and** `<stage>/<name>.d` → `DIA-Quant` |

The `DIA-Quant` row matters: it keeps DIA-NN quantifying against the real `.d`
chromatograms while skipping the conversion. `DIA-Quant` does not trigger diaTracer, which
is exactly why it is safe here.

It refuses two things outright: an input that is not a `.d`, and a staging path with `.d`
inside a *directory* name (which Java's replace-every-occurrence would mangle).

#### Runtime

From the paper (34 dia-PASEF files, 32 logical cores, 128 GB): diaTracer **661 min total,
~19.4 min per file**; the rest of FragPipe including DIA-NN quant **87 min**. DIA-NN
library-free on the same data took 1486 min. **The pseudo-MS/MS mzML are reusable** — a
semi-tryptic re-search from existing mzML took 110 min, so conversion is a one-time cost
per cohort. No GPU is needed anywhere in this path. No published RAM minimum; `--ram` maps
to diaTracer's `-Xmx`, and real OOMs are reported, so be generous in `sbatch --mem`.

#### Licensing — read this before quoting it to a user

| Tool | Terms |
|---|---|
| MSFragger, IonQuant, **diaTracer** | Free for **academic / non-commercial / educational** use under one shared academic licence ("MSFragger + IonQuant + diaTracer Suite", The Regents of the University of Michigan — the PDF is byte-identical for all three). Commercial users license via **Fragmatics** (fragmatics.com, info@fragmatics.com). |
| **DIA-NN inside FragPipe** | *"Available for both academic and commercial use **within FragPipe**"* (`FRAGPIPE-THIRD-PARTY-LICENSES.txt`). The academic-only limit on Demichev's *standalone* DIA-NN does not apply to the bundled binary — **but the grant is scoped to use within FragPipe**, so do not treat `tools/diann/` as a commercially-usable standalone DIA-NN. |

The commercial gate on this route is MSFragger/IonQuant/diaTracer, **not** DIA-NN. A
commercial user is not locked out; they license the three gated tools via Fragmatics.
Existing University of Michigan commercial licences must be re-confirmed with Fragmatics.

**⚠ The core-facility clause — the most important line here.** The academic licence
prohibits *"any work product, report or deliverable by LICENSEE, for a commercial or
for-profit organization."* That wording reaches an **academic core delivering results to an
industry client**, which is ordinary fee-for-service work. The quote is verbatim; whether it
binds a given engagement is a legal question, not one for this skill. If a user asks whether
their fee-for-service or industry-funded work is covered, tell them to ask
info@fragmatics.com and their own licensing office. **Do not answer it for them.**

**You cannot fully script the install.** MSFragger/IonQuant/diaTracer come from a web form
requiring name, academic email, organization and licence acceptance. There is a separate
licence-key token for automated installs (Bioconda uses exactly this). So a human must
accept the licence once and mint a key; `acquire_tools.sh` can download FragPipe itself but
never the three gated jars. That is why `--config-tools-folder` is mandatory.

**Turning DIA-NN off does not remove it.** DIA-NN is also the deep-learning predictor inside
MSBooster rescoring — this preset sets `msbooster.rt-model=DIA-NN`,
`msbooster.spectra-model=DIA-NN`, `msbooster.im-model=DIA-NN`, and `CmdMSBooster.java`
passes the DIA-NN path in. Unchecking "quantify with DIA-NN" only skips the quant step; the
executable still runs during rescoring unless MSBooster is switched to Koina models. And on
**Linux there is no DIA-NN-free quant path at all**: FragPipe's Skyline backend returns null
on non-Windows (`Skyline.java`) and only ships `SkylineRunner.exe`, so DIA-NN off means IDs
and a library but no quant matrix.

**⚠ Do not swap in your own DIA-NN 2.x via `--config-diann`.** Beyond breaking the output
shape (2.x writes parquet), the DIA-NN 2.x licence adds: *"Usage of DIA-NN in the cloud or
on any computer setups that are not in your immediate possession is only permitted under
Collaborative use."* A shared cluster is arguably not "in your immediate possession". The
bundled 1.8.2 Beta 8 carries no such clause. This clause is **also worth knowing for the
`diann_*` workflows**, which pin DIA-NN 2.x and are routinely run on a shared HPC — flag it
to the user rather than deciding for them.

**Comparing it against DIA-NN.** That is the point of having both — `compare_searches.py`
(SKILL.md step 11a). FragPipe bundles DIA-NN 1.8.2 Beta 8 while the `diann_*` workflows pin
2.x, so part of any difference is engine version rather than the spectrum-centric approach.
Say so rather than crediting it all to the method.

### FragPipe (opt-in; adapter required)
- Build an `.fp-manifest` from the files + data type (DIA/DDA).
- `<cmd> --headless --workflow <.workflow> --manifest <m> --workdir <out>
  [--config-tools-folder $FRAGPIPE_TOOLS_FOLDER]`
- **Adapter:** `combined_protein.tsv` per-sample `MaxLFQ Intensity` columns →
  DIA-NN-shaped `report.parquet`.
- Needs Java 9+; MSFragger/IonQuant must already be licensed/present.
- **Reading Thermo `.raw`:** MSFragger reads `.raw` directly via
  `tools/ext/thermo/BatmassIoThermoServer.exe` — a **Windows .NET-Framework** exe that on
  **Linux needs `mono` on PATH** (the .NET 8 runtime does NOT run it — Framework ≠ Core).
  Native Windows/macOS: works out of the box. On a mono-less Linux cluster you get
  *"Could not find Batmass-IO Thermo binary"* + `Scans = 0` → either install mono
  (conda-forge `mono`, then `export PATH=<monoenv>/bin:$PATH` in the job) or pre-centroid
  to mzML. (On HIVE, mono 6.12 is installed at `/quobyte/proteomics-grp/conda_envs/mono`.)
- **Other DDA gotchas:** IonQuant MBR needs experiment groups in the manifest — set
  `ionquant.mbr=0` if you only want peptide/protein IDs; **profile-mode mzML** fail
  `CheckCentroid` (centroid during conversion, or read `.raw` so MSFragger vendor-centroids);
  the LFQ-MBR template ships **without** a `database.db-path` line (add one); and the
  headless run's exit code can read 0 even on a crash — **verify the output file exists**.
- **Which output file to verify (don't false-alarm):** the file depends on the run
  type. A **single-experiment / peptide-ID run** (`ionquant.mbr=0`, one experiment)
  writes **`psm.tsv` + `peptide.tsv` + `protein.tsv` + `ion.tsv`** — there is **no**
  `combined_peptide.tsv`/`combined_protein.tsv` (those only appear with **multi-experiment
  IonQuant LFQ**). A watcher that checks for `combined_*.tsv` will wrongly report
  "no output" on a perfectly good ID run. Success check = `peptide.tsv` (and `psm.tsv`)
  present and non-empty; use `combined_protein.tsv` only when the manifest defines ≥2
  experiment groups and you asked for the cross-run LFQ matrix.

## The adapters are the test surface

DIA-NN's path is fully wired and native. The **Sage and FragPipe → DE-contract
adapters are the part that needs real-data validation** (flagged in each
workflow's `VALIDATION.md`). Verify the protein×run matrix shape and that
intensities are in the expected (linear, pre-log) scale before trusting DE output.
