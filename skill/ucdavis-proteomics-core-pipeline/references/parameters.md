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
| Orbitrap, resolution unknown | **automatic calibration** (`--mass-acc 0`) | generic orbitrap names |
| Instrument not detected | **automatic calibration** | fallback |

Automatic calibration is DIA-NN's own recommended default — it optimises mass
accuracy on the first run and reuses it. We fall back to it (never to a guessed
number) whenever the instrument class can't be pinned down.

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
