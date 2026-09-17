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
| Orbitrap, **a level with no documented tier** — resolution unknown, or outside 30k–240k (e.g. 15k MS2) | that level **measured with DIA-NN before the search**; a level with a tier keeps it (both flags omitted from the cfg, plan `measure_with_diann`) | generic orbitrap names, or resolutions outside the table |
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

## Provenance (DE-LIMP rule #2)
Every emitted value is tagged in the `rationale`:
- `data-type-default` — chosen from the instrument/acquisition (e.g. Astral → 10 ppm)
- `auto-calibration` — left to the engine because the class couldn't be pinned
- `measured with DIA-NN before the search` — an Orbitrap level with no documented DIA-NN value;
  `mass_accuracy_plan: measure_with_diann` and `mass_accuracy_documented` (in the rationale and
  at the top of the sidecar) are what `diann_parallel.mass_acc_measure_plan()` reads, for the
  chain and the single-shot search alike
- `universal trypsin/LFQ default`
- `user-override (validated SOP)` — forced via the workflow's `param_overrides`

The rationale is written to `<params>.rationale.json` and is pulled into the
reproducibility bundle automatically. Surface it to the user so a derived default
is never mistaken for a confirmed setting.

## Overriding with a validated SOP
Two ways, both honored:
- `param_overrides` in the workflow.yaml (e.g. `{"--mass-acc": 8}` or
  `{"fragment_tol": {"ppm": [-15, 15]}}`) — merged on top of the estimate, tagged
  `user-override`. An override of one mass-accuracy flag is a value for that level, so
  nothing is measured: the other level is written from DIA-NN's table (`{"--mass-acc": 8}` at
  120k/15k gives `--mass-acc-ms1 7 --mass-acc 8`, plan `pinned`). When the table has no value
  for the other level (outside 30k–240k, resolution unknown, instrument not identified), the
  override is **refused** with a message naming the missing flag: written alone it would fix
  that level at 20 ppm. Give both flags, pass the resolutions, or override neither.
- Ship a full validated `params_file` in the workflow — used verbatim, estimation
  skipped entirely. Use this when a method is locked and must not move.
