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
| Orbitrap by MS2 resolution | 240k→4, 120k→7, 60k→10, 30k→15 ppm | (when resolution is known) |
| Bruker timsTOF (dia-PASEF / ddaPASEF) | 15 / 15 ppm | name contains "tims" |
| SCIEX TripleTOF / ZenoTOF | 20 / 20 ppm | name contains tripletof/zenotof/sciex |
| Orbitrap, resolution unknown | **automatic calibration** (flags omitted) | generic orbitrap names |
| Instrument not detected | **automatic calibration** | fallback |

Automatic calibration is DIA-NN's own recommended default — it optimises mass
accuracy on the first run and reuses it. We fall back to it (never to a guessed
number) whenever the instrument class can't be pinned down.

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
`Resolution Comp. (ppm):`). Take them from the method or the instrument operator and
pass `--ms1-resolution`/`--ms2-resolution`. Before trusting `--from-mzml` with any
converter, check the file: `grep -c 'MS:1000800' run.mzML` — `0` means it will pin
nothing.

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
- `universal trypsin/LFQ default`
- `user-override (validated SOP)` — forced via the workflow's `param_overrides`

The rationale is written to `<params>.rationale.json` and is pulled into the
reproducibility bundle automatically. Surface it to the user so a derived default
is never mistaken for a confirmed setting.

## Overriding with a validated SOP
Two ways, both honored:
- `param_overrides` in the workflow.yaml (e.g. `{"--mass-acc": 8}` or
  `{"fragment_tol": {"ppm": [-15, 15]}}`) — merged on top of the estimate, tagged
  `user-override`.
- Ship a full validated `params_file` in the workflow — used verbatim, estimation
  skipped entirely. Use this when a method is locked and must not move.
