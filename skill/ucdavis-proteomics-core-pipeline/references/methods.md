# Publication-ready Methods section (`make_methods.py`)

Generates a drop-in LC-MS/MS Methods section straight from facility raw data, plus
the correct UC Davis Proteomics Core instrument-grant acknowledgment. Can run as
part of a full analysis or **standalone** (just `--raw` at the facility data).

## What it extracts vs. defaults
- **Extracted from the raw metadata** (shown with its source in a parameter table):
  - Bruker `.d` → `analysis.tdf`: instrument (`GlobalMetadata.InstrumentName`),
    acquisition software+version, acquisition mode (`Frames.MsMsType`: 9=dia-PASEF,
    8=ddaPASEF), m/z range (`MzAcqRange*`), 1/K₀ range (`OneOverK0AcqRange*`),
    TIMS ramp/accumulation time (`Frames`), and the dia-PASEF window scheme
    (`DiaFrameMsMsWindows`: count, isolation width, collision-energy range).
  - Thermo `.raw` → identified by the **facility filename prefix** (`FL*`→Fusion
    Lumos, `Ex*`→Exploris 480); detailed parameters come from the instrument method
    (not readable here) and are flagged for confirmation.
- **Facility defaults**, every one tagged `[facility default — confirm]` so nothing
  is silently fabricated (DE-LIMP rule #2): the **LC column defaults to a PepSep C18
  10 cm × 150 µm, 1.5 µm** (override with `--lc-column`), mobile phases, the
  CaptiveSpray/PepSep emitter, and the LC system/gradient (which the raw `.d` does
  not store — supply from lab records).

With `--de-dir`, it adds a Data-processing paragraph from the run's
`de_provenance.json` (engine, FDR, quantification, thresholds, citation).

## Instrument grant acknowledgments (verified 2026-06)
Picked by instrument metadata **or** facility filename prefix, from
https://proteomics.ucdavis.edu/instrument-grant-acknowledgments:

| Instrument | Prefix | Acknowledgment |
|---|---|---|
| Orbitrap Fusion Lumos | `FL` | NIH S10 grant **S10OD021801** |
| Orbitrap Exploris 480 | `Ex` | NIH S10 grant **S10OD026918-01A1** |
| Bruker timsTOF | — | Dr. Neil Hunter / **Howard Hughes Medical Institute** |

An instrument not in this registry yields a placeholder pointing at the source URL.
**Grant wording must be exact** — the script cites the verified grant numbers and
links the source page; confirm the current wording there before publishing. To add
an instrument, extend the `ACKS` table in `make_methods.py`.

## Output + workflow
- `methods.md` — the drop-in prose (LC, MS, [Data processing], parameter table,
  Acknowledgments).
- `methods_params.json` — the extracted parameters, machine-readable.
- Render to Word with `to_docx.py --in methods.md --out methods.docx`.

The agent should **verify the draft against the parameter table and polish the
prose** (the example in `~/Documents/DataAnalysis/.../timsTOF_Methods` shows the
target quality), resolve each `[facility default — confirm]`, and keep the
acknowledgment verbatim.
