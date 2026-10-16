# Search parameter estimation

Parameters are **derived from the data type**, not hand-maintained per workflow.
`estimate_params.py` takes the engine, acquisition, and detected instrument and
emits a complete, ready-to-run params file plus a `rationale` that tags every
value's provenance. This replaces shipping static `.cfg`/`.json` files (which go
stale and hide which values were actually chosen vs. guessed).

## What actually varies by data type

Only two things genuinely depend on the data; everything else is a stable
trypsin/LFQ default:

1. **Mass tolerances** — set by the mass analyzer / resolution.
2. **DIA vs DDA window mode** — `wide_window`/`--window` and chimera handling.

## Known-good mass tolerances (DIA-NN README, verified 2026-06)

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

## Cysteine: which alkylating agent, or none (`--cys-alkylation`)
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

The rationale is written to `<params>.rationale.json` and is pulled into the
reproducibility bundle automatically. Surface it to the user so a derived default
is never mistaken for a confirmed setting.

## Overriding a setting (the user's choice)
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
