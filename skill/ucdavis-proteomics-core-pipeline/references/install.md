# Installation reference — what's automatic, what isn't

Design goal: **a biologist with a fresh laptop runs the skill and it works.**
`setup.sh` does the installing; this doc explains what it does and the one
genuine manual fallback.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 19: Installing the SKILL itself — marketplace only
- 76: What `setup.sh` installs automatically (no admin / no sudo)
- 169: The search engines
- 181: The ONE manual step: Docker on macOS (only for DIA-NN / DIA data)
- 197: macOS + Sage + Bruker/Thermo
- 209: Windows
- 240: Assume nothing about the user's computer (golden rule 7)
- 275: Step 0: the one-time setup, and what `setup.json` says
<!-- /contents -->

## Installing the SKILL itself — marketplace only

Install via the plugin marketplace, **added by its HTTPS URL**. In a terminal, one line (then
start Claude Code, or `/reload-plugins` if it is open):
```
claude plugin marketplace add https://github.com/bsphinney/DE-LIMP.git && claude plugin install ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core
```
Inside Claude Code the same is two commands, **entered one at a time** (pasted together they
fail): `/plugin marketplace add https://github.com/bsphinney/DE-LIMP.git`, then `/plugin install ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core`. Automatic updates are off for this marketplace until the
user turns them on: `/plugin` → Marketplaces → ucdavis-proteomics-core → Enable auto-update.

**Why the URL and not `bsphinney/DE-LIMP`.** Claude Code's documentation ("Host and maintain a
marketplace", read 2026-10-08) says that for the `owner/repo` shorthand it probes
`ssh -T git@github.com` and clones over SSH when the probe succeeds, and that an `https://…`
URL is cloned over HTTPS. Its plugin-commands reference says an `https://` URL ending in
`.git` is always cloned as a git repository. Without `.git`, only github.com and gitlab.com
URLs are; any other URL is fetched as a `marketplace.json`. So the `.git` is kept.
The marketplace name stays `ucdavis-proteomics-core`, so every other command is unchanged.
On 2026-10-07 a Windows laptop (Git Bash, no GitHub SSH key) with the marketplace added by the
shorthand could not refresh it: `claude plugin update` warned "marketplace not refreshed: SSH
host key is not in your known_hosts", then "SSH authentication failed", and said "already at the
latest version (2.10.0)" while main shipped 2.11.2. Its catalogue clone had an HTTPS remote all
along. **Already added by the shorthand? Keep it.** Do not remove and re-add it:
`claude plugin marketplace remove` also uninstalls the plugin. The update route below works
for it.

**To update now.** Inside the skill, step 0 does it: `bash scripts/skill_version.sh --update`
refreshes the catalogue over HTTPS when it is behind GitHub main, runs the update, and checks the
version that landed against main. It never says "up to date" while behind. By hand, in Git
Bash or a terminal, one line at a time:
```
git -C ~/.claude/plugins/marketplaces/ucdavis-proteomics-core pull --ff-only
claude plugin update ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core
```
then `/reload-plugins`. The `git pull` is what took that laptop from 2.10.0 to 2.11.2.
**A computer on 2.11.2 or older needs these two lines once, by hand.** Its own step 0 offers
only `claude plugin update` (`--update` arrived in 2.12.0), which on a laptop without a GitHub
SSH key compares against the old catalogue and stays put: staff were still on 2.11.2 on
2026-10-16, two releases later. Without
it, when Claude Code cannot refresh the catalogue, `claude plugin update` compares against the
old one, so its "already at the latest version" can be wrong. (With `CLAUDE_CONFIG_DIR` set, the
catalogue is under that folder instead of `~/.claude`.) `--update` also sets
`CLAUDE_CODE_PLUGIN_PREFER_HTTPS=1`, which Claude Code documents for machines without a GitHub
SSH key ("skip the probe and clone over HTTPS"). Whether that switch also covers the refresh of
an existing catalogue is not documented and was not tested here. The `git pull` is the route
that was seen to work.

**Do not hand-copy this directory into
`~/.claude/skills/<name>/` while the plugin is installed** — the plugin already owns
the name, so the copy is silently inert. `claude plugin list` reports it as
`@skills-dir: Not loaded — the name is already taken`: it looks installed and does
nothing, and it will drift from the released version without any error to tell you.

If you already have a hand-copied directory, remove it and reinstall from the
marketplace; check with `claude plugin list` that exactly one copy is loaded, and
that its version matches the one you expect.

## What `setup.sh` installs automatically (no admin / no sudo)

It downloads **micromamba** (a single static binary, ~5 MB) into
`~/.proteomics-pipeline/micromamba/` if no conda/mamba is already present, then
creates one conda environment — solved with `--override-channels -c conda-forge -c bioconda`, so
a `~/.condarc` listing `defaults` cannot pull in its R ≤ 4.3 builds — containing:

| package | why |
|---|---|
| `python=3.11` + `pyarrow` + `pyyaml` | the skill's Python scripts + parquet adapters |
| `r-base` (≥4.5) | the DE step |
| `bioconductor-limpa` (1.2.5, noarch) | the DPC-Quant DE pipeline — bioconda's only build, a floor: step 2a then installs **limpa ≥ 1.4.0** from the Bioconductor 3.23 source repository |
| `bioconductor-limma` (3.66) + `r-statmod` + `r-data.table` | both DE pipelines (limpa's and limma's dependencies) |
| `r-nanoparquet` | reading a DIA-NN `report.parquet` |
| `r-arrow`, `r-dplyr`, `r-tidyr` | reading parquet + the MaxLFQ matrix builder |
| `sage-proteomics` | the DDA search engine |
| `proteowizard` (msconvert) | `.d` → mzML for Sage — **Linux only on bioconda**, and built without vendor readers (`.raw` goes through ThermoRawFileParser) |
| `thermorawfileparser` (2.0.0.dev; linux-64, osx-64, osx-arm64) | reads Thermo `.raw` for `detect_acquisition.py` (step 2): acquisition, instrument and the **acquired precursor m/z range**. Without a parser every `.raw` is `unknown` and a DIA search falls back to 380–980. bioconda's build is the **self-contained** one — it needs no .NET on the machine. Installed in a **separate, non-fatal** step after the env exists, so a platform without a build cannot take R/limpa down with it |
| `pythonnet` (conda-forge, 3.1.0, noarch) | lets `thermo_resolution.py` read the **Orbitrap MS1/MS2 resolution** from each `.raw`'s scan trailer with the RawFileReader DLLs the parser ships — ThermoRawFileParser itself never outputs it, and `estimate_params.py` needs it to pin DIA-NN's documented Orbitrap tolerances. Same separate, non-fatal step as the parser; only what is missing is installed. It runs on a .NET 8 root with `Microsoft.NETCore.App` (`ensure_dotnet8.sh`'s; a self-contained parser's bundled runtime cannot host it) |
| `pandas` (conda-forge) | ad-hoc tables in Python (a user's first reach). In the same separate step — not the main solve — so an env that already exists gets it on a re-run of `setup.sh` |

> Superseded by "Step 0: the one-time setup, and what `setup.json` says" below.

`setup.sh` also runs **`ensure_dotnet8.sh`** (Linux and macOS; not with `--check`): one .NET 8
root with `Microsoft.NETCore.App` ≥ 8.0.17 + `Microsoft.AspNetCore.App` in
`~/.proteomics-pipeline/dotnet8`, which DIA-NN's `.raw` reader, a framework-dependent
ThermoRawFileParser and the resolution reader all use. It reuses a good root at once. A failure
(it needs internet the first time) is a note, not a stop; `setup.json`'s `dotnet8`
(`root`, `note`) says what happened.

**google-genai ≥ 2.25.0** (conda-forge, noarch) gives the podcast its Gemini voices
(`make_podcast.py`, SKILL.md step 9e; Gemini 3 TTS needs its `client.interactions`). A
CoreOmics submission gets the podcast by default, and on HIVE its job runs with this env's python. `setup.sh`
installs it in its own step (not with `--check`), so a failure is a note that nothing else
depends on. `setup.json`'s `podcast` (`google_genai`, `required`, `ok`) reports it.

> Superseded by "Step 0: the one-time setup, and what `setup.json` says" below.

**limpa ≥ 1.4.0** (`run_de.R` reads the report with its `readDIANN(annotation.columns =)`): step
2a installs it into the env from the Bioconductor 3.23 source repository (`PROTEOMICS_LIMPA_BIOC`
overrides the release; limpa is pure R, so no compiler). The env stays on R 4.5. The last check
of `setup.sh` is `packageVersion("limpa") >= "1.4.0"`; `setup.json` → `limpa` =
`{version, required, ok, source}`, and on failure it prints `ERROR: limpa ...` with its fix
(re-run `setup.sh`, which on a cluster runs itself on a compute node, or the exact
`install.packages()` line, on a cluster inside a job, never on a login node), and exits 1. An
env still on limpa 1.2.x runs `run_de.R` anyway (the argument was `extra.columns` there;
`limpa_compat.R`) and records it in `de_provenance.json`. → `references/environment.md`, "The R stack".

> Superseded by "Step 0: the one-time setup, and what `setup.json` says" below.

Everything lands under `~/.proteomics-pipeline/`. `source activate.sh` puts it on
PATH so `Rscript`, `python`, `sage`, `msconvert`, `ThermoRawFileParser` resolve to the env,
and exports `DOTNET_ROOT` for `ensure_dotnet8.sh`'s .NET once that exists (a `DOTNET_ROOT`
already set is kept). Nothing is installed system-wide; deleting `~/.proteomics-pipeline/`
fully uninstalls it.

`setup.sh --check` reports readiness without installing anything.

> Superseded by "Step 0: the one-time setup, and what `setup.json` says" below.

**Thermo `.raw` readiness** — `setup.json` has `ready_for.thermo_raw` and a
`thermo_raw_reader` object (`command`, `source`, `version`, `dotnet_needs`, `dotnet_root`,
`note`). It asks `detect_acquisition.py --check-reader` — the same parser search and .NET
check step 2 runs, plus one `--version`; no `.raw` is read. The parser is looked for in this
order: `$THERMORAWFILEPARSER`, PATH, the pipeline env, then shared copies in
`$THERMORAWFILEPARSER_SHARED` (default: the UC Davis Core's
`/quobyte/proteomics-grp/tools/ThermoRawFileParser/ThermoRawFileParser`; set it to your own
site's copy, or empty for none). When it is `false`, `note` is the exact fix.
`thermo_raw_reader.resolution_reader` (`ready`, `python`, `dll_dir`, `dotnet_root`, `reader`,
`note`) says whether the Orbitrap resolution can be read: it runs `thermo_resolution.py --check`
(loads .NET and the two DLLs, opens no `.raw`). It does not gate `ready_for.thermo_raw` —
without it the user is asked for the resolution instead. The DLLs are looked for in
`$THERMO_RAWFILEREADER_DIR`, then beside the parser; the interpreter is
`$THERMO_RESOLUTION_PYTHON`, else the one running step 2 if it has pythonnet, else the
pipeline env's python.

**.NET 8 (only for a framework-dependent parser, and for DIA-NN on Linux).** The Core's
shared parser, and the release's `-net8` zip run as `dotnet ThermoRawFileParser.dll`, are
*framework-dependent*: their `ThermoRawFileParser.runtimeconfig.json` needs **both**
`Microsoft.NETCore.App` 8 and `Microsoft.AspNetCore.App` 8. Missing, the parser exits 131
("You must install .NET") or 150 ("You must install or update .NET"). `detect_acquisition.py`
checks this before reading anything, and picks the first place that has both —
`$DOTNET_ROOT`, `ensure_dotnet8.sh`'s install (`$PROTEOMICS_DOTNET_DIR`, default
`~/.proteomics-pipeline/dotnet8`), `$DOTNET_CORE_SDK_ROOT` (HIVE: `module load
dotnet-core-sdk/8.0.4`), the `dotnet` on PATH, .NET's default location. With none it stops
once, up front, with the fix. `bash scripts/ensure_dotnet8.sh` (needs internet; on a cluster login
node its download moves itself to a compute node) installs both runtimes there, and adds
`AspNetCore` in place to an older NETCore-only install. The same install serves DIA-NN 2.6,
which needs `Microsoft.NETCore.App` ≥ 8.0.17 to read
`.raw` on Linux (HIVE's 8.0.4 module is fine for the parser but *not* for DIA-NN). Public
source: https://dotnet.microsoft.com/download/dotnet/8.0 (the ASP.NET Core Runtime 8.0 includes
both).

## The search engines

- **Sage** (DDA): installed by `setup.sh` from bioconda. `acquire_tools.sh` then
  finds it on PATH — no separate download.
- **DIA-NN** (DIA): license-gated, not on conda. By platform:
  - **Linux:** `acquire_tools.sh` downloads the free DIA-NN *Academia* binary.
  - **UC Davis HIVE:** Proteomics Core members reuse the Core's native DIA-NN builds
    under `/quobyte/proteomics-grp/dia-nn/`; everyone else gets the Linux download above.
  - **macOS:** *no native build exists.* See the one manual step below.
- **FragPipe** (opt-in): downloaded on demand; MSFragger/IonQuant are license-gated
  and can't be auto-downloaded — `acquire_tools.sh` says so if you opt in.

## The ONE manual step: Docker on macOS (only for DIA-NN / DIA data)

> Superseded by "Assume nothing about the user's computer (golden rule 7)" below.

DIA-NN can't run natively on a Mac. `build_diann_docker.sh` builds an image from
DIA-NN's **own official** Academia zip + Dockerfile (no third-party images), but
that needs Docker:

1. Install Docker Desktop (free): https://docs.docker.com/desktop/setup/install/mac-install/
2. Open it once so setup finishes (steady whale icon in the menu bar).
3. Re-run the analysis — `build_diann_docker.sh` builds the image and records it as
   `DIANN_DOCKER_IMAGE`; the AI continues automatically.

DDA data (Sage) and all DE steps need **no** Docker. A Mac user with DDA data that
is already in mzML never touches Docker at all.

## macOS + Sage + Bruker/Thermo

Thermo `.raw` is converted to mzML with ThermoRawFileParser, which bioconda builds for
macOS too (`setup.sh` installs it). Bruker `.d` goes to `msconvert`, which is Linux-only on
bioconda. Bioconda's recipe also builds it from ProteoWizard's `pwiz-src-without-v` source,
without the vendor readers, and on HIVE it failed at once on a `.raw` (a Core member, 2026-09-29).
Don't count on it for vendor files anywhere. For `.d` with Sage, options in order of
preference:
1. If the data is DIA, use DIA-NN (reads `.d`/`.raw` natively) — no conversion.
2. Convert to mzML elsewhere (Windows msconvert, or ProteoWizard's vendor-enabled Docker
   image) and point the skill at the `.mzML` files.

## Windows

Two supported routes — pick per the user:

1. **WSL2 (recommended, smoothest).** `wsl --install` (Ubuntu), then run the skill
   inside it and follow the Linux flow — the conda env, the bash orchestration scripts,
   Sage, and the DIA-NN Linux binary all work as-is.
2. **Native Windows.** The engines themselves ship **native Windows builds** — DIA-NN
   (Windows GUI + CLI, and it reads Thermo `.raw` natively, so the Linux-only `.NET 8`
   step is not needed), Sage, FragPipe, and ProteoWizard/`msconvert` (originally a
   Windows tool). Python + R + limpa run natively too. The skill's *orchestration
   scripts are bash*, so run them from a POSIX shell — **Git Bash** or **MSYS2** — or
   drive the engines directly per their Windows docs. Docker Desktop is a third route.

**Driving HIVE from Windows (`hive_remote`) needs none of the above** — no local Python,
R or engines; every script runs on HIVE through `hive_exec.sh`. Git Bash is enough:
- **No rsync in Git Bash.** `hive_exec.sh --put/--get` fall back to `scp` automatically
  (same key and options). The one thing scp can't express is a trailing slash on the
  source (`dir/` = "copy the contents" under rsync); that is refused with a message, so
  drop the slash.
- **`python3` may be the Microsoft Store stub** (`...\WindowsApps\python3`), which opens
  the Store instead of running. `check_access.sh` reports it as
  `local_python3.usable: false`; in `hive_remote` that is fine: Core staff's CoreOmics steps
  then run on HIVE too (`coreomics_runs_on: hive`; `references/core-submissions.md` "Python
  on this computer").
- **The HIVE username is the plain UC Davis id**, not the Windows login (`AD3+<your HIVE user>`);
  see `references/access.md` → "Windows (Git Bash)".

Public download links for every program are in `references/search-engines.md`
("Public program sources"), so a user on any OS can obtain them without HIVE access.

## Assume nothing about the user's computer (golden rule 7)

7. **Assume nothing about the user's environment — this skill runs for anyone.** Most
   users are **not** UC Davis Core and **not** on HIVE: they're on macOS, Windows
   (**WSL2** recommended, but **native Windows works too** — the engines ship Windows
   builds), or generic Linux. HIVE / `/quobyte/proteomics-grp` is a Core **fast-path,
   never a requirement**. For **every** program you have the user run, give a **public
   source** anyone can obtain (download URL / `pip`/install command) and a path that
   works on *their* OS — never hardcode an internal `/quobyte/...` path without its
   public fallback. → `references/search-engines.md` ("Public program sources") +
   `references/install.md`.


The user may be a biologist on a fresh laptop with no R, no Python packages, no
Docker. **Do not ask them to install things by hand.** Run `setup.sh` — it
installs everything that can be installed without admin rights into one
self-contained conda env. Only one thing ever needs the user's hands (Docker
Desktop, and only for DIA-NN on macOS), and `setup.sh`/`build_diann_docker.sh`
print the exact step. Never dump a generic "please install R/limma/..." message;
the scripts handle it.


- Anything in the env (R, limpa, Sage, pyarrow) → re-run `setup.sh`; relay its
  `notes`. Never tell the biologist to "install R/limma/..." by hand.
- macOS + DIA data + no Docker → run `build_diann_docker.sh` and relay its exact
  steps (install Docker Desktop, open once). Don't silently switch to a DDA engine.
- macOS + DDA + Bruker `.d` → msconvert is Linux-only on bioconda; see
  `references/install.md` ("macOS + Sage"). Prefer DIA-NN if the data is DIA; else convert
  to mzML first. Thermo `.raw` is converted with ThermoRawFileParser (bioconda has a macOS
  build).
- Unrecognised instrument → parameters fall back to automatic mass calibration;
  say so plainly. This is never a stop: no organism or instrument is "unsupported".
- FragPipe license (MSFragger/IonQuant) → surface the `tools.json` note and the fix.
- Acquisition/instrument unknown → ask the user; never guess into a multi-hour run.

## Step 0: the one-time setup, and what `setup.json` says

`setup.json` also carries **`next_steps`**: for each item that is not ready (`de`, `limpa`,
`dia`, `dda`, `thermo_raw`, `resolution_reader`) one sentence with the exact fix for this
platform and when it matters (e.g. `thermo_raw` only for Thermo `.raw` input); `{}` when all is
ready. setup.sh prints each as a `[setup] next:` line too. Do what the entry says before the step
it blocks; the bullets below are the reasoning behind each one.

**Then the one-time setup (install everything that's missing):**

```
bash scripts/setup.sh
source ~/.proteomics-pipeline/activate.sh        # puts R, python, sage on PATH
```

`setup.sh` installs micromamba (no admin), then a conda env with R + limpa +
limma + arrow + Sage + Python/pyarrow, and writes `~/.proteomics-pipeline/setup.json`.
It also provisions **.NET 8** (`scripts/ensure_dotnet8.sh`, Linux and macOS: one root in
`~/.proteomics-pipeline/dotnet8` with Microsoft.NETCore.App ≥ 8.0.17 + AspNetCore, which
DIA-NN's `.raw` reader, the Core's ThermoRawFileParser and the resolution reader all use)
and installs `thermorawfileparser`, `pythonnet` and `pandas` into the env in their own step,
so an existing env gets them on a re-run. Neither step is fatal — `setup.json`'s
`dotnet8.note` and `notes` say what did not happen.
Sourcing `activate.sh` makes every later step use those interpreters — **source it
in the same shell before running any later step** (or prefix later commands).

Read `setup.json` and **gate on `ready_for`**:

- `ready_for.de` false → DE can't run: R cannot read `report.parquet` the way this env's limpa
  does (`de_reader`: arrow for limpa 1.2.x, nanoparquet for limpa ≥ 1.4.0, arrow on both
  methods); re-run `setup.sh` and report any `notes`. `ready_for.de` is the default method's
  (`dpc`); `de_reader.maxlfq` true means `--method maxlfq` can still run (arrow + limma, no limpa). (A SILAC ratio analysis, `run_silac.R`,
  needs only limma, data.table, arrow and jsonlite in that R.)
- `limpa.ok` false (setup.sh then ends `ERROR: limpa ...` and exits 1) → the env's limpa is
  older than 1.4.0 (bioconda has only 1.2.5; setup.sh installs 1.4 from Bioconductor 3.23):
  re-run `setup.sh` (on a cluster it runs on a compute node). run_de.R still reads the report
  with limpa 1.2.x, through arrow, and records that it did (`references/environment.md`); without
  nanoparquet its contaminant share is then read from `Precursor.Normalised`. Upgrading an env
  that other analyses run on changes their limpa mid-analysis: Core staff ask the Core's director first.
- `ready_for.dia` false → on macOS this means Docker is missing. Run
  `bash scripts/build_diann_docker.sh` and relay its instructions (install Docker
  Desktop, open it once), then continue. Don't silently fall back to a DDA engine.
- `ready_for.dda` false → Sage/R not ready, or (macOS) msconvert is unavailable for
  `.d`/`.raw` → mzML; tell the user and, if their data is DIA, route to DIA-NN.
- `ready_for.thermo_raw` false **and the input has Thermo `.raw`** → fix it before step 2:
  do exactly what `thermo_raw_reader.note` says (usually re-run `setup.sh`, or `bash
  scripts/ensure_dotnet8.sh`; on a cluster both move to a compute node), then `setup.sh --check` until it is
  true. `ready_for.dia` does not cover this: without a working parser every `.raw` is
  `unknown` and a DIA search falls back to 380–980. Irrelevant for `.d`/mzML-only input.
- `thermo_raw_reader.resolution_reader.ready` false only means step 2 cannot read the
  Orbitrap resolution and the user will be asked for it; its `note` is the fix.

In `hive_remote` the toolchain lives on HIVE, so run it there instead (the laptop needs none
of it): `bash scripts/hive_exec.sh 'bash ~/proteomics-pipeline/scripts/setup.sh'`, and read
`~/.proteomics-pipeline/setup.json` on HIVE the same way. On a cluster login node `setup.sh`
(and `--check`) runs itself on a compute node (`on_compute.py`: one blocking `srun`, up to 2 h for
an install, printed first): run that call in the background and wait for it. Never run a conda,
R or .NET install on the login node by hand; HIVE's compute nodes have internet. Run it before
the first search: it installs .NET 8, so writing a `.raw` search later does not wait on a
download job. Exit 3: the compute-node run did not finish (SLURM refused it, a time limit) and
`setup.json` was not rewritten — its `ready_for` is old; read the srun lines and run it again.
Exit 2: it could not get there (no python3 there), and says how.


This step is idempotent — on a machine that's already set up it just verifies and
returns in seconds.
