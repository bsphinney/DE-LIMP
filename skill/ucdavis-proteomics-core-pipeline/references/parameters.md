# Search parameter estimation

Parameters are **derived from the data type**, not hand-maintained per workflow.
`estimate_params.py` takes the engine, acquisition, and detected instrument and
emits a complete, ready-to-run params file plus a `rationale` that tags every
value's provenance. This replaces shipping static `.cfg`/`.json` files (which go
stale and hide which values were actually chosen vs. guessed).

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 27: Where the defaults come from (the skill version pins the search)
- 37: What actually varies by data type
- 47: Known-good mass tolerances (DIA-NN README, verified 2026-06)
- 143: Sage tolerances are derived
- 151: Defaults that are NOT data-dependent (universal)
- 157: SILAC labels (`--labels silac2|silac3 --label-preset ... [--channel-norm spec|run]`)
- 167: Cysteine: which alkylating agent, or none (`--cys-alkylation`)
- 207: CPU-hours a search reserves (`cpu_hours_reserved`)
- 221: Golden rule 2: never fabricate parameters
- 228: Provenance (DE-LIMP rule #2)
- 254: Overriding a setting (the user's choice)
- 275: Golden rule 1: confirm before committing compute
- 292: 4. Resolving the defaults, and the one confirmation (step 4)
- 366: 6b. Estimating the search parameters (step 6b)
<!-- /contents -->

## Where the defaults come from (the skill version pins the search)

Take a user from **raw MS files → differentially expressed proteins** using
parameters derived from their **data type** (instrument + acquisition), never
ad-hoc numbers. Those defaults ship with this skill — `scripts/estimate_params.py`
holds the instrument → mass-accuracy table, `scripts/resolve_defaults.py` picks the
engine and version, `scripts/make_presets.py` writes FragPipe/Radiant configs. So
**the skill version fully determines the search**, and nothing is fetched at run
time. Every number carries its source; surface it rather than asserting it.

## What actually varies by data type

> Superseded by "6b. Estimating the search parameters (step 6b)" below.

Only two things genuinely depend on the data; everything else is a stable
trypsin/LFQ default:

1. **Mass tolerances** — set by the mass analyzer / resolution.
2. **DIA vs DDA window mode** — `wide_window`/`--window` and chimera handling.

## Known-good mass tolerances (DIA-NN README, verified 2026-06)

> Superseded by "6b. Estimating the search parameters (step 6b)" below.

| Instrument class | DIA-NN MS1 / MS2 | matched on |
|---|---|---|
| Orbitrap Astral | 4 / 10 ppm (assumes 240k MS1) | name contains "astral" |
| Orbitrap by resolution, per level | 240k→4, 120k→7, 60k→10, 30k→15 ppm; between tiers interpolated (tagged) | both MS1 and MS2 resolution known and inside 30k–240k |
| Bruker timsTOF (dia-PASEF / ddaPASEF) | 15 / 15 ppm | name contains "tims" |
| SCIEX TripleTOF / ZenoTOF | 20 / 20 ppm | name contains tripletof/zenotof/sciex |
| Orbitrap, **one level outside 30k–240k** (e.g. 15k MS2), the other inside | that level **measured with DIA-NN before the search**, refused unless it lands in a plausible band (MS2 3–30, MS1 1.5–25 ppm) over at least 2 agreeing runs (a one-file search: its one run; runs that disagree only warn when the median is at or below the SOP and no run is above 1.5× it, the SOP then being pinned either way), then **floored at the SOP** (`max(measured, SOP)`, `estimate_params.SOP_MASS_ACC` = MS2 20 / MS1 7 ppm): the measurement is used only where it is wider. The level with a tier keeps it and is never floored (both flags omitted from the cfg, plan `measure_with_diann`) | both resolutions known, one outside the table |
| Orbitrap, **resolution unknown** | **automatic calibration** — *not* measured: with no resolution both levels would be, and a measured MS1 is what DIA-NN warns about. Pass `--ms1-resolution`/`--ms2-resolution`, or pin the two flags | generic orbitrap names with no resolution — since 2.5.1 only when step 2 could not read it from the `.raw` scan trailer (see below) or the input is an mzML without `MS:1000800` |
| Instrument not detected | **automatic calibration** | fallback |

**Outside the table nothing is extrapolated.** The old log-log fit turned the 15,000 MS2 that
both Orbitraps in the FRAN re-search pilot acquire (Exploris 480 60k/15k and 120k/15k, Fusion
Lumos 120k/15k) into 23.3 ppm and pinned it for the whole cohort, with nothing behind the number
but a curve. Instead that level is **measured with DIA-NN** on the cohort's representative runs
before the search, and pinned for it: step 1b of the 5-step chain, or a probe between the library
and the search in a single-shot search — so a machine without SLURM measures it too
(→ `diann_parallel.md`). The measurement is *adapted from* the README's item 6 of "Changing
default settings", which reads in full: "run DIA-NN on several representative runs (best to use
any suitable empirical library, as this is the quickest) with **Unrelated runs** option checked
and review the 'Averaged recommended settings for this experiment' values reported at the end of
the log". It is not that procedure: it uses the predicted library, one DIA-NN per run, and pins
the median of what each run printed, not DIA-NN's averaged line (which differed: MS2 16 against
a median of 14 on the validation cohort).

**A level with a tier keeps it.** At 120k/15k the cfg's rationale records MS1 7 ppm under
`mass_accuracy_documented` and only MS2 is measured. Measuring MS1 too pinned 4.2 ppm, and DIA-NN
2.7.0 then warned on every pass: `the MS1 mass accuracy setting (4.2 ppm) deviates significantly
from the value recommended (7 ppm) for the Orbitrap resolution of this run (120000)`.

**Neither flag is written on its own.** DIA-NN 2.7.0 fixes BOTH levels when either is given —
`WARNING: note the mass accuracy settings used by DIA-NN, automatic optimisation will not be
performed as at least one of MS1/MS2 mass accuracies is user-provided`, and `--mass-acc-ms1 7`
alone gave `Mass accuracy will be fixed to 2e-05 (MS2) and 7e-06 (MS1)` (HIVE srun 23528991) — so
a cfg with only `--mass-acc-ms1 7` would silently search MS2 at 20 ppm. The cfg omits both; the
measurement writes both. A search that cannot measure (`--one-step`: no library exists before the
search) says so and leaves DIA-NN's first-run optimisation. The same holds for an SOP override
of one flag (below): the other level is written from the table beside it, and when the table has
no value for that level `estimate_params.py` refuses rather than write a lone flag.

**An open question, with numbers.** DIA-NN 2.7.0 also has a resolution-based MS2 value for 15k
that its README does not document: pinned at the measured 14 ppm, the search logs `WARNING: the
MS2 mass accuracy setting (14 ppm) deviates significantly from the value recommended (25 ppm)
for the Orbitrap resolution of this run (15000)`; at 25 ppm it is silent. On the validation run
the wider settings report more precursors and the narrower ones about as many protein groups
(table in `diann_parallel.md`). The skill currently pins the measurement.

Automatic calibration (instrument not detected) is DIA-NN's own default — it optimises mass
accuracy on the first run and reuses it ("use this mode for preliminary analyses only", as
DIA-NN 2.7.0 prints). We fall back to it (never to a guessed number) whenever the instrument
class can't be pinned down; such a cfg is not parallel-safe.

Each pinned flag carries **its own level's** source: a documented 120k MS1 is never tagged with
an interpolated or out-of-table MS2's provenance (it used to be — a documented tier read
"EXTRAPOLATED").

### Orbitrap resolution: where it is, and where it is not

The Orbitrap rows need the MS1/MS2 resolving power. `--from-mzml <file>` reads it from
the mzML's `mass resolving power` (`MS:1000800`) terms — **and ThermoRawFileParser does
not write that term.** Measured 2026-09-16 on HIVE with TRFP 2.0.0.0: the full indexed
mzML of an Exploris 480 run (60k / 15k) and a Fusion Lumos run (120k / 15k), ~1 GB each,
held **zero** `MS:1000800` terms, so `read_mzml_resolution()` returned `(None, None)`.
Nothing errors: the class falls to `orbitrap_generic`, mass accuracy to automatic
calibration, and — because the 5-step chain needs it pinned — a large cohort to the
slow single-shot path. TRFP's other outputs do not carry it either: its metadata JSON
`mass resolution` (`MS:1000011`) is `0.5` on every run checked (a header tolerance, not
resolving power), and its `query` JSON has no resolution attribute.

> Superseded by "4. Resolving the defaults, and the one confirmation (step 4)" below.

Where the numbers do live: the instrument method text (`Orbitrap Resolution = 60000`
then `= 15000` on that Exploris; `= 120K` then `= 15K` on the Lumos) and every scan's
trailer (`FT Resolution:` on Exploris, `Orbitrap Resolution:` on Fusion Lumos — not
`Resolution Comp. (ppm):`). **Since 2.5.1 step 2 reads the trailer itself:**
`detect_acquisition.py` runs `thermo_resolution.py` (Thermo's RawFileReader DLLs via
pythonnet) and reports `ms1_resolution`/`ms2_resolution` — verified against the embedded
instrument method on Exploris 480 (120k/15k DIA, 60k/15k DDA) and Fusion Lumos (60k/15k)
runs. Pass them with **`--resolution-source detected`**; numbers the user or operator gave
go in without it (the source defaults to `user`), and values copied from an earlier cfg
take `--resolution-source cfg`. The source is recorded — manifest `resolution.source_label`,
the rationale's `class_label` and every mass-accuracy line say "read from the raw file (scan
trailer)", "supplied by the user", etc. — so methods text never calls a typed-in number
measured. (`--ms1-res`/`--ms2-res` are accepted aliases.) Before trusting `--from-mzml`
with any converter, check the file: `grep -c 'MS:1000800' run.mzML` — `0` means it will
pin nothing.

⚠ **Automatic calibration means OMITTING `--mass-acc`/`--mass-acc-ms1`, not
setting them to 0.** `--mass-acc 0` fixes the tolerance at a literal 0 ppm — the
log reads `Mass accuracy will be fixed to 0 (MS2) and 0 (MS1)` and the search
returns **0 identifications**. The same applies to `--window 0`, which DIA-NN
rejects with `scan window radius should be a positive integer`.

## Sage tolerances are derived

Sage's docs give **no** instrument-specific tolerances, so `estimate_params.py`
derives Sage's ppm windows from the same per-instrument logic (high-res
Orbitrap/Astral → ±10 ppm fragment; timsTOF → ±20 ppm; SCIEX → ±40 ppm; unknown →
±20 ppm safe high-res default) and tags them `derived from DIA-NN per-instrument
recommendation (Sage docs give none)`. `wide_window`/`chimera` follow acquisition.

## Defaults that are NOT data-dependent (universal)
Trypsin/P, 1–2 missed cleavages, peptide length 7–30, charge 2–4, fixed
carbamidomethyl (C). Variable mods are **off by default** — the DIA-NN README
notes variable mods don't improve depth for relative quant; pass `var_mods: "ox"`
in the workflow (or `--var-mods ox`) to add Ox(M).

## SILAC labels (`--labels silac2|silac3 --label-preset ... [--channel-norm spec|run]`)

> Superseded by "6b. Estimating the search parameters (step 6b)" below.

DIA-NN plexDIA for SILAC samples (DIA only): writes `--fixed-mod SILAC,0.0,KR,label`,
`--channels` (masses from UniMod), `--original-mods` and the channel-normalisation flag, each
tagged, and the sidecar's `labels` record; never `--lib-fixed-mod` (it would label a
FASTA-predicted library twice). SILAC flags in `--overrides` are refused. The presets, the design
question and why: `references/silac.md`.

## Cysteine: which alkylating agent, or none (`--cys-alkylation`)

> Superseded by "6b. Estimating the search parameters (step 6b)" below.

Each alkylating agent leaves its own fixed modification on C, and a search for the wrong one (or
for one at all, on samples never alkylated) identifies almost no cysteine peptide: three on-bead
Core searches had 0.03–0.5% of identified peptides with a C, against 7.8% on an S-Trap prep with
DTT/IAA (staff reports, 2026-10-09). The choice is its own option (`--overrides` cannot remove a
flag), one map in `estimate_params.CYS_REAGENTS`; masses and names are UniMod's (unimod.org,
read 2026-10-15):

| `--cys-alkylation` | agent | fixed on C (UniMod) | DIA-NN cfg | Sage `static_mods` |
|---|---|---|---|---|
| not given | assumed iodo-/chloroacetamide, **DEFAULT** | Carbamidomethyl (4) +57.021464 | `--unimod4` | `{"C": 57.0215}` |
| `cam` | iodoacetamide or chloroacetamide | Carbamidomethyl (4) +57.021464 | `--unimod4` | `{"C": 57.0215}` |
| `none` | not reduced and alkylated | — | (nothing) | `{}` |
| `propionamide` | acrylamide | Propionamide (24) +71.037114 | `--fixed-mod UniMod:24,71.037114,C` | `{"C": 71.037114}` |
| `nem` | N-ethylmaleimide | Nethylmaleimide (108) +125.047679 | `--fixed-mod UniMod:108,125.047679,C` | `{"C": 125.047679}` |
| `carboxymethyl` | iodoacetic acid | Carboxymethyl (6) +58.005479 | `--fixed-mod UniMod:6,58.005479,C` | `{"C": 58.005479}` |
| `mmts` | S-methyl methanethiosulfonate | Methylthio (39) +45.987721 | `--fixed-mod UniMod:39,45.987721,C` | `{"C": 45.987721}` |
| `pyridylethyl` | 4-vinylpyridine | Pyridylethyl (31) +105.057849 | `--fixed-mod UniMod:31,105.057849,C` | `{"C": 105.057849}` |

DIA-NN's `--fixed-mod [name],[mass],[sites]` (README) applies the modification in the FASTA
digest, so the predicted library carries it; how well DIA-NN's predictor handles a modification
other than carbamidomethyl is not measured here. FragPipe (`make_presets.py --cys-alkylation`) sets
MSFragger's fixed C row to the agent's mass (0 for `none`); Radiant's DIA-NN library is predicted
with the same DIA-NN flags (`run_search.py`, `radiant_parallel.py`, `make_radiant_library.py`,
which reuses a library only when it was built with the same choice).

SKILL.md step 3 asks only for an enrichment experiment (IP / pull-down / proximity labelling, or
an on-bead digest), and then which agent. `make_methods.py` lists the fixed C modification only
when the parameters fix it, with `none` says cysteine was searched unmodified because the samples
were not reduced and alkylated, and with another agent names it. `sample_quality.py --report`
warns when a fixed C modification was searched and under 1% of the identified precursors hold a C
(`sample_quality.CYS_MIN_FRACTION`), pooled **or in any run** with at least 100 identified
precursors (`CYS_MIN_PRECURSORS_RUN`) -- a cohort half of whose samples were never alkylated stays
above 1% pooled. It names the causes and fixes: not alkylated (`none`), alkylated with another
agent than searched (the matching choice), or a protein mix naturally poor in cysteine (histone
or collagen-rich preps), which can fall under 1% though it was alkylated.

## CPU-hours a search reserves (`cpu_hours_reserved`)
Every search that writes job scripts states the most CPU-hours it can use: for each job it submits,
CPUs per task × time limit (hours) × array tasks, summed (`scripts/cpu_hours.py`, the one rule),
read from the `#SBATCH` lines actually written. It is an **upper bound** -- a job that finishes
early uses less. `run_search.py` prints it (`[run_search] cpu_hours_reserved: …` and
`cpu_hours_reserved` in its JSON and in `search_provenance.json`) for the single-shot, two-job,
5-step chain, Radiant and FragPipe routes; `diann_parallel.py` and `radiant_parallel.py` print it
in their JSON; each `submit.sh` has it in a header comment. `gpu_hours_reserved` sits beside it,
by the same rule from the GPU lines (`--gres=gpu…`, `--gpus`, `--gpus-per-task`,
`--gpus-per-node`): 0 today, since no route of the skill asks for a GPU. Every job script the skill
writes carries an `#SBATCH --time`, so for a SLURM job neither figure is null
(`tests/test_cpu_hours.py` checks every generator). The Core board approves a request under
10,000 CPU- or GPU-hours automatically and sends 10,000 or more to a person.

## Golden rule 2: never fabricate parameters

2. **Never fabricate parameters.** If a value isn't a shipped default or given by
   the user, say so — don't invent an FDR, organism, or instrument. Every default
   carries its source (`ppm_source`); quote it rather than asserting the number.
   (DE-LIMP architectural rule #2.)

## Provenance (DE-LIMP rule #2)
Every emitted value is tagged in the `rationale`:
- `data-type-default` — chosen from the instrument/acquisition (e.g. Astral → 10 ppm)
- `auto-calibration` — left to the engine because the class couldn't be pinned
- `measured with DIA-NN before the search` — an Orbitrap level with no documented DIA-NN value;
  `mass_accuracy_plan: measure_with_diann` and `mass_accuracy_documented` (in the rationale and
  at the top of the sidecar) are what `diann_parallel.mass_acc_measure_plan()` reads, for the
  chain and the single-shot search alike
- `universal trypsin/LFQ default`
- `user override: <setting> = <value>, set by Core staff; reason: <why>` — forced with
  `--overrides` (or the workflow's `param_overrides`), with `--override-by` and
  `--override-reason` (`estimate_params.override_source`, the one wording; the sidecar's
  `overrides` holds the same as fields). Never "validated SOP": an override may be a one-off
  workaround (a cohort searched at 25 ppm to cover a ~20 ppm calibration error was once
  recorded as an SOP), and the skill cannot tell. The sidecar travels to the client
  (`reproducibility/inputs/`, the session zip), so it names the role. Who set it, by name
  (`--override-by`, or the login that ran it, said to be that), is in the staff-only
  `<params>.staff.json` (`staff.py`; never delivered, never zipped), which `overrides.*.staff_record`
  names. No reason: "reason not given [ask the user -- confirm]"

> Superseded by "6b. Estimating the search parameters (step 6b)" below.

The rationale is written to `<params>.rationale.json` and is pulled into the
reproducibility bundle automatically. Surface it to the user so a derived default
is never mistaken for a confirmed setting.

## Overriding a setting (the user's choice)

> Superseded by "6b. Estimating the search parameters (step 6b)" below.

Two ways, both honored:
- `param_overrides` in the workflow.yaml, or `--overrides` (e.g. `{"--mass-acc": 8}` or
  `{"fragment_tol": {"ppm": [-15, 15]}}`) — merged on top of the estimate, tagged as a user
  override with who set it and why (above). A Sage override replaces the derived rationale
  entry it changes (`precursor_tol` drops `precursor_tol_ppm`). An override of one mass-accuracy flag is a value for that level, so
  nothing is measured: the other level is written from DIA-NN's table (`{"--mass-acc": 8}` at
  120k/15k gives `--mass-acc-ms1 7 --mass-acc 8`, plan `pinned`). When the table has no value
  for the other level (outside 30k–240k, resolution unknown, instrument not identified), the
  override is **refused** with a message naming the missing flag: written alone it would fix
  that level at 20 ppm. Give both flags, pass the resolutions, or override neither. **This is a
  hard exit**, and any cfg and `.rationale.json` already at `--out` are deleted with it, so a
  refusal never leaves a stale parameter set behind for the next step to read. A `workflow.yaml`
  whose `param_overrides` sets one mass-accuracy flag for such an instrument used to produce a
  cfg and now fails — see the changelog entry.
- Ship a full validated `params_file` in the workflow — used verbatim, estimation
  skipped entirely. Use this when a method is locked and must not move.

## Golden rule 1: confirm before committing compute

1. **Confirm before committing compute.** A search is multi-hour. Show the resolved
   defaults (engine + version, mass accuracy + source) with the organism/design and
   get an explicit "go" before running the engine. **One confirmation, not a menu.**
   After that go, the jobs that carry out the confirmed plan need no further yes, so the
   work can run unattended overnight: the chain's later steps, retries from the
   `watch_run.sh` playbook, DE, set tests, figures, the podcast. A new search with settings
   the user has not seen (another database, engine or parameters) needs their go again,
   unless they told you to carry on yourself ("run it overnight, decide as you go"): then
   log each such decision (`log_decision.py`) and go ahead. Core staff: each job is
   recorded on the board with the CPU-hours (and any GPU-hours) its scripts reserve, and the
   board approves it by itself when both are under 10,000; a job of 10,000 or more of either
   is asked for there first and waits for a person (`references/board.md` section 3).
   Sharing results with a collaborator (`core_submission.py bioshare send`) is a **second,
   separate** explicit yes — never bundled into the compute confirmation.

## 4. Resolving the defaults, and the one confirmation (step 4)

```
python3 scripts/resolve_defaults.py \
    --acquisition DIA --instrument "timsTOF HT" \
    [--ms1-resolution <MS1> (--ms2-resolution <MS2> | --ms2-analyzer <ITMS|mixed>) [--resolution-source detected]] \
    [--engine diann] [--env env.json] [--fasta db.fasta] [--threads 32] --dest ./wf
```

Writes `./wf/workflow.manifest.json` (engine + pinned version, mass accuracy and its
source, DE method) and, for FragPipe/Radiant, generates the config file itself. Pass
`--env` (step 1's `detect_env.sh` output) so an engine that **cannot run on this
machine** fails here with a usable alternative rather than mid-search.

**Orbitrap resolution — pass it, and say where it came from.** For a Thermo Orbitrap other
than the Astral, pass step 2's top-level `ms1_resolution`/`ms2_resolution` with
`--resolution-source detected` (read from the `.raw` scan trailer). Numbers the user told
you: omit `--resolution-source` (it defaults to `user`); copied from an earlier run's
configuration: `--resolution-source cfg`. (`--ms1-res`/`--ms2-res` are the same flags.)
**Ion-trap MS2:** for step 2's `ms2_ion_trap` files (common on Fusion Lumos/Eclipse OT/IT
methods) pass `--ms1-resolution <MS1> --ms2-analyzer ITMS` and **no** `--ms2-resolution` (there
is none; the two together are refused), and never ask the user for one. A file whose MS2 is
read in BOTH analyzers (step 2's `mixed_files`) is `--ms2-analyzer mixed` and gets the same
ion-trap tolerances. DIA-NN then pins
neither mass-accuracy level (there is no Orbitrap MS2 value, and either flag fixes both), so it
is left to DIA-NN's calibration; Sage uses its documented low-res fragment window (±0.4 Da —
with ±10 ppm a real Sage 0.14.7 run matched nothing on ion-trap fragments). Quote the
manifest's `ppm_source`, which says which. A cohort in which some files are ion-trap and some
Orbitrap MS2 (`resolution_mixed`) needs either `--ms2-analyzer mixed` for all of them or
separate searches. The
manifest's `resolution.source_label` and `instrument_label` carry the source — quote it as
written, and never say "read from the data" about numbers the user gave.

**`needs_confirmation: true` means a question is still open, even with exit code 0.** It is
set for an Orbitrap with no resolution on a route that uses it (DIA-NN, Radiant). Put
`ask_user` to the user as part of the one confirmation below — get the MS1/MS2 resolution
(or run `detect_acquisition.py`, which reads it from the `.raw`), then re-run this step with
it. If the user does not know, say that mass accuracy is left to DIA-NN's per-run
calibration (results then depend on file order, and a cohort that would run as the 5-step
chain runs as one single-node search), and continue.

**Confirm once, then run.** State the pick in one breath — data type, engine +
version, mass accuracy + where it came from, FASTA, DE method — and get a yes:

> timsTOF HT dia-PASEF → DIA-NN 2.7.0, MS1/MS2 15 ppm (DIA-NN README), mouse
> UP000000589 + contaminants, limpa DPC-Quant + limma. Run it?

> Exploris 480 DIA, MS1 120,000 / MS2 15,000 resolution read from the raw file (scan
> trailer) → DIA-NN 2.7.0, MS1 7 ppm (DIA-NN table), MS2 measured on representative runs
> before the search, human UP000005640 + contaminants, limpa DPC-Quant + limma. Run it?

Do **not** turn this into a menu. There is exactly one confirmation before compute,
and it is this one.

- **`--engine` is honoured, never overridden.** If the user names an engine that
  isn't the default for their data type, it runs anyway and `notes` says so. Never
  silently substitute a different engine.
- **`alternatives`** lists the other engines valid for this data type. Mention them
  in one clause; don't make the user adjudicate.
- **Overrides are explicit, and they are the user's.** `--ms1-ppm/--ms2-ppm` force a value
  the user chose. Pass `--override-by "<who chose it>"` and `--override-reason "<why, in their
  words>"` (e.g. "instrument ~20 ppm off calibration; user chose 25 ppm"): the manifest's
  `ppm_source` and `overrides` record each as a **user override** with its value, "set by Core
  staff", and why — never as an SOP, which the skill cannot know it is. The manifest and the cfg
  sidecar reach the client, so who set it, by name, goes only into the staff-only
  `<file>.staff.json` beside them. Without `--override-by` that is the login that ran it, and the
  record says so. Without a reason the record says "reason not given". Don't override without
  asking the user why, and keep names out of the reason.
- Parameters come from the skill, so **record the skill version**, not a registry
  commit: `provenance.py` reads it from `.claude-plugin/plugin.json` into
  `reproducibility/environment/skill.txt` and `run_manifest.json` (`skill.version`).
  `registry.defaults_version` in the manifest is only the date of the defaults table
  (e.g. `2026-09-23`), not a version.

## 6b. Estimating the search parameters (step 6b)

Search parameters are **derived from what the data is**, not hand-maintained.
For DIA-NN and Sage, generate them (for FragPipe/Radiant see 6c — step 4 already
wrote their config):

```
python3 scripts/estimate_params.py --engine <diann|sage> \
    --acquisition <DIA|DDA> --instrument "<detected instrument>" \
    --precursor-mz-range <LO> <HI> \
    [--ms1-resolution <MS1> (--ms2-resolution <MS2> | --ms2-analyzer <ITMS|mixed>) [--resolution-source detected]] \
    [--var-mods ox] [--cys-alkylation <agent>|none] \
    [--labels silac2|silac3 --label-preset <K4R6,K8R10|K8R10|...> [--channel-norm spec|run]] \
    [--overrides '<the user's values as JSON>' \
     --override-by "<who chose them>" --override-reason "<why, in their words>"] \
    --fasta-meta ./search.fasta.meta.json \
    --out ./wf/params.<cfg|json>
```

**`--labels` for a SILAC sample** (step 3): writes DIA-NN's plexDIA declaration into the cfg
(`--fixed-mod SILAC,0.0,KR,label`, `--channels`, `--original-mods`, the channel-normalisation
flag) and the sidecar's `labels` record, the normalisation tagged DEFAULT unless given. Never
pass these flags as `--overrides` (refused), and never add `--lib-fixed-mod`: on a library
predicted from the FASTA it labels K and R twice (`K(SILAC)(SILAC)`); both search routes refuse
such a cfg. DIA with DIA-NN only. → `references/silac.md`.
**`--cys-alkylation`** carries step 3's answer: the agent's fixed modification on C (`cam`:
`--unimod4`; another agent: `--fixed-mod UniMod:<n>,<mass>,C`; Sage: `static_mods C`), or `none`
for cysteine unmodified, and the Methods name it. The sidecar's `cys_alkylation` records the
choice and who made it. Not asked → omit it: fixed carbamidomethyl-C, tagged DEFAULT.
`sample_quality.py --report` (step 8c) warns when a fixed C modification was searched and under
1% of the identified precursors hold a C, pooled or in any one run: confirm the prep, then
re-search with `none` or the agent the lab used.
**Always pass `--precursor-mz-range`**, taking `precursor_mz_range` straight from
step 2's `detect_acquisition.py` output. Without it a DIA range falls back to
380–980 — which on a timsTOF method acquiring 299.5–1200.5 silently discards both
tails, does not error, and is invisible in the results. The rationale tags the
value `measured …` or `FALLBACK …` so you can tell which you got; if it says
FALLBACK, say so to the user before committing to a multi-hour search.
**For a DDA run `precursor_mz_range` is the MS1 survey scan range**
(`precursor_mz_range_source: ms1_survey_scan`): a DDA run picks its precursors from its
survey scans, so none lies outside it. Step 2 reads it from the Thermo scan filter
(`Full ms [350.0000-1500.0000]`), the mzML MS1 scan window, or a timsTOF `.d`'s
`MzAcqRange`. Pass it exactly as for DIA. **`estimate_params.py --engine diann --acquisition
DDA` refuses to write a cfg without it** — there is no 380–980 fallback for DDA (SET28's
Exploris surveyed 350–1500 and its PSMs spanned 360–1315 m/z, so a third would have been
cut). A DDA file with no readable range gets a NOTE (in its `reason`, and the top-level
`dda_range_note`), not a warning: Sage's DDA config does not use the range, so it does not set
`needs_confirmation`. For a DIA-NN search, ask the user for the survey scan range from the
instrument method.
A `null` on a **DIA** input means the range could not be read, the `warnings` say why, and
the search will run on the FALLBACK unless you supply the range yourself.

**DIA-NN on DDA** (`--engine diann --acquisition DDA`): the cfg carries `--dda`. Every route
reads it from the cfg — the 5-step chain included — and `run_search.py` refuses a cfg whose
`--dda` disagrees with the bundle's acquisition. **Nothing is measured for DDA**: DIA-NN logs
no scan-window radius in DDA mode, so a probe would only run to its 3600 s timeout (SET28 lost
~3 h that way). The level that has a DIA-NN tier keeps it (60k MS1 → 10 ppm); a level DIA
would measure (a 15k MS2) is pinned at the SOP (MS2 20 ppm / MS1 7 ppm) and tagged
`DEFAULT, not user-confirmed` (`mass_accuracy_default` in the sidecar). The Methods text and the
run record say so from the run's provenance ("fragment (MS2) 20 ppm (Core SOP default for this
instrument, not measured on these data)"). **Tell the user**, and offer `--overrides` if their
instrument needs another value. `--window` stays unset, and the
chain has no step 1b for DDA.

Always pass `--fasta-meta` (step 6's sidecar): it carries the contaminant tag, so
the cfg gets `--cont-quant-exclude Cont_` and contaminants are identified but kept
out of DIA-NN's normalisation and its quantities for sample proteins. Both the single and 5-step parallel
DIA-NN paths read this cfg, so the flag applies to library generation *and*
analysis. With `--contaminants none` the flag is correctly absent.
The estimator keys mass tolerances on the instrument class from DIA-NN's
known-good table (Astral → MS1 4/MS2 10 ppm; timsTOF → 15 ppm; Orbitrap by
resolution, 30k–240k; unidentified → automatic calibration). An **Orbitrap level
outside 30k–240k**, e.g. a 15k MS2, is not extrapolated: the cfg omits *both*
mass-accuracy flags (DIA-NN fixes both levels when either is given), the sidecar plans
`measure_with_diann`, and the level that has a tier keeps it
(`mass_accuracy_documented`, e.g. MS1 7 ppm at 120k). DIA-NN measures the missing level
on representative runs before the search — step 1b of the parallel chain, or a probe
inside the single-shot search job (step 7) — and the measurement is **refused unless it
is plausible**: MS2 3–30 ppm, MS1 1.5–25 ppm, from at least 2 runs agreeing to within
50% (a one-file search pins its one run; a disagreement whose median is at or below the SOP,
with no run above 1.5× it, only warns, since the SOP is pinned either way). What is pinned is then **max(measured, SOP)** (`estimate_params.SOP_MASS_ACC`: MS2
20 ppm, MS1 7 ppm) — the measurement is used only where it is WIDER than the SOP,
because a tolerance tighter than the SOP costs identifications and buys nothing.
Both numbers are recorded; a documented level is never floored. An Orbitrap whose **resolution is unknown** is *not* measured (both levels would be,
and a measured MS1 is what DIA-NN warns about): it falls to automatic calibration, and the
chain declines it. Pass the same resolution flags as step 4 — step 2's
`ms1_resolution`/`ms2_resolution` with `--resolution-source detected`, or, for files in
`orbitrap_resolution_unknown`, the user's answer with no `--resolution-source` (it defaults to
`user`); for `ms2_ion_trap` files, `--ms2-analyzer ITMS` instead of `--ms2-resolution` (step 4 —
Sage then writes `fragment_tol {"da": [-0.4, 0.4]}`, `deisotope: false` and a larger
`bucket_size`, Sage's own low-res MS/MS settings; DIA-NN omits both mass-accuracy flags).
**Pass the SAME `--ms2-analyzer` here as in step 4** — only this step changes Sage's config, and
`run_search.py` refuses a Sage search whose config does not match the manifest (a ppm fragment
window for an ion-trap MS2, or a Da window for an Orbitrap MS2) before converting or submitting
anything. With
them a 60k/15k Fusion Lumos becomes `orbitrap_untabled`: MS1 is pinned at
10 ppm from DIA-NN's table and the 15k MS2 is measured by step 1b (`measure_with_diann`), so
the 5-step chain can run. If step 2 reported `resolution_mixed`, do not pass one value for the
whole cohort. The rationale's `class_label`, `resolution.source_label` and each mass-accuracy
line then say where the numbers came from, and that is what the methods text quotes; if the
sidecar still says `needs_confirmation: true`, relay `ask_user`. Without them mass accuracy is
left to DIA-NN's per-run calibration (results then depend on file order), and a cohort
that would run as the 5-step parallel chain runs as one single-node search. If the user
does not know, say so and continue. An `--overrides` that sets only one of
the two flags gets the other level from the table, and is refused when the table has no value
for it (written alone, that flag would fix the other level at 20 ppm): relay the refusal to the
user. It sets DIA/DDA window mode from acquisition, and uses
standard trypsin/LFQ defaults for the rest. **It prints a `rationale` tagging
every value's provenance** — surface this to the user (and it flows into the
methods text), so a derived default is never mistaken for a confirmed setting.
Use the resulting file as `--params` in step 7.
