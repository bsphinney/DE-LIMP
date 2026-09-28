# Publication-ready Methods section (`make_methods.py`)

Generates a drop-in LC-MS/MS Methods section straight from facility raw data, plus
the correct UC Davis Proteomics Core instrument-grant acknowledgment. Can run as
part of a full analysis or **standalone** (just `--raw` at the facility data).

## What it extracts vs. defaults
- **Extracted from the raw metadata** (shown with its source in a parameter table):
  - Bruker `.d`, read by `bruker_method.py`. It reads small metadata files only, never the
    `analysis.tdf_bin`, and opens every sqlite file read-only and immutable.
    - `analysis.tdf`: instrument and serial, acquisition software and version, MS method name,
      acquisition date, m/z and 1/K₀ ranges, and acquisition mode (`Frames.MsMsType`: 9 =
      dia-PASEF, 8 = ddaPASEF).
    - `Frames`: polarity, TIMS ramp and accumulation times (duty cycle), frames per cycle and
      cycle time (the median MS1-to-MS1 interval).
    - `DiaFrameMsMsWindows`: window count, TIMS ramps, windows per ramp, isolation width,
      spacing and overlap, the m/z the windows cover, and each window's collision energy.
    - The method's settings, from the Properties of the first MS/MS frame: the ion source
      (named by the file's own `DisplayValueText` code table), capillary voltage, dry gas, dry
      temperature and the collision-energy ramp.
    - `<run>.m/diaSettings.diasqlite`: each window's 1/K₀ bounds.
    - `<run>.m/microTOFQImpacTemAcquisition.method`: the control software that wrote the
      method.
    - `<run>.m/hystar.method`: `ColumnInfo`, when an operator entered a column.
    - `HyStarMetadata.xml`: the LC system, its vendor and serial, and the LC method and its
      run time.
    - `SampleInfo.xml`: HyStar version and the autosampler tray (e.g. `96Evotip`).
  - The collision-energy ramp is written as a ramp only when every window's recorded energy
    matches it at the window's 1/K₀ midpoint, to within 1 eV. Otherwise the per-window range
    is written, tagged.
  - A series is checked for one method. If the runs differ in instrument, MS or LC method,
    ranges, column or window scheme, the Methods says so, above the parameter table.
  - Thermo `.raw` → identified by the **facility filename prefix** (`FL*`→Fusion
    Lumos, `Ex*`→Exploris 480); detailed parameters come from the instrument method
    (not readable here) and are flagged for confirmation.
- **The analytical column**, most authoritative first:
  1. `--lc-column`, given by the user.
  2. HyStar `ColumnInfo`, when an operator entered it.
  3. `--column-log FILE`: a CSV or JSON export of STAN's `maintenance_events`. It needs
     `event_type`, `event_date`, `column_vendor` and `column_model`, and can carry
     `instrument`, `column_serial` and `first_run`. The latest `column_change` at or before
     the first acquisition is used, and its install date is cited as the source. A change
     inside the series, a change that names no column, or a log covering several instruments
     without `--column-log-instrument` falls back to the default, with the reason printed.
     STAN keeps this log in PG Farm, which needs credentials, so the script never connects to
     it; export the rows and pass the file.
  4. Otherwise the facility's standard column, from STAN's column catalogue (PepSep MAX C18,
     10 cm × 150 µm, 1.5 µm, part 1893483), **tagged `[facility default — confirm]`**.
- **Never in a .d, always tagged**:
  - Column temperature. For a timsTOF, STAN's operator-reported 50 °C is given, tagged
    `[facility default — confirm]`.
  - The emitter. The default 20 µm CaptiveSpray emitter is tagged
    `[facility default — confirm]`.
  - Mobile phases.
  - Evotip type, loading and peptide amount.
  - The %B gradient of an Evosep method. Its fixed, named method is given instead.
  - A non-Evosep LC's gradient table, which is not parsed. It is tagged
    `[not recorded — confirm]`.
- **Tags**:
  - A value no record holds is tagged `[not recorded — confirm]`, never
    `[facility default — confirm]`.
  - A `.d` that could not be read is named on stderr and in the note under Mass spectrometry.
    When it is the only run, its values are tagged `[raw file not readable here — confirm]`.
  - The ddaPASEF precursor-selection settings are not extracted, and the text says so.

With `--params` / `--search-prov` / `--workflow-manifest`, it adds a **Database search**
paragraph and a search-parameter table (engine, the version that ran, cleavage rule, missed
cleavages, peptide length, charge and m/z range, fixed/variable modifications, mass tolerances,
precursor FDR, library mode, MBR), read by `search_record()` — the same reader
`make_deposit.py` uses for the SDRF. A value in no record prints as
`____ [not recorded — confirm]`, never as a default.

With `--submission <session>` (a CoreOmics submission attached by `submission_report.py`),
it adds a **Sample preparation** section from the record, through `prepared_by()`, the one
reading of the form. When the lab sent peptides, it says the submitting laboratory prepared
them and carries no Core-side placeholder; the form's own words (buffer, beads) go in a note
for the author, not into the prose, so nothing is added that the form does not state. When the
Core prepared them, the protocol is `____ [not recorded — confirm]`. `make_deposit.py` passes
the flag itself, and its PRIDE sample-processing protocol uses this section instead of its
TO-FILL.

With `--de-dir`, it adds a **Differential expression** paragraph from the run's
`de_provenance.json` (pipeline, quantification, the q-value filters as applied, design,
contrasts, significance rule, citation). Significance is stated as run_de.R applies it —
adj.P.Val only; |log2FC| is a volcano reference line, not a filter.

When the raw files cannot be read from where it runs (e.g. finalize away from the data), pass
`--instrument` / `--acquisition` from the session record: the Methods are then written with
those two values (sourced as the session record) and every acquisition value blank and tagged
`[raw file not readable here — confirm]`. `session.py finalize` does this automatically —
see `references/deposit.md`.

## Instrument grant acknowledgments (verified 2026-06)
The acknowledgment comes from
https://proteomics.ucdavis.edu/instrument-grant-acknowledgments.

It is picked by the **instrument name** first: the name read from the raw file, or from the
session record. The facility filename prefix (`FL`, `Ex`) is used only when the instrument is
unknown and every file is a Thermo `.raw` with that prefix. It is never used for a `.d`, because
a timsTOF run named `FLAG_IP_1.d` or `Exp3_HeLa.d` is not a Thermo run. A named instrument that
is not in the registry gets the placeholder below, not a guess from the filename.

| Instrument | Prefix | Acknowledgment |
|---|---|---|
| Orbitrap Fusion Lumos | `FL` | NIH S10 grant **S10OD021801** |
| Orbitrap Exploris 480 | `Ex` | NIH S10 grant **S10OD026918-01A1** |
| Bruker timsTOF (Pro 2 and HT) | — | Dr. Neil Hunter / **Howard Hughes Medical Institute** |

The web page names only the timsTOF Pro 2 for the HHMI acknowledgment. Brett confirmed on
2026-09-25 that it covers the timsTOF HT too.

An instrument not in this registry yields a placeholder pointing at the source URL.
**Grant wording must be exact** — the script cites the verified grant numbers and
links the source page; confirm the current wording there before publishing. To add
an instrument, extend the `ACKS` table in `make_methods.py`.

## Output + workflow
- `methods.md` — the drop-in prose (LC, MS, sequence database, [database search],
  [differential expression], parameter tables, Acknowledgments).
- `methods_params.json` — the extracted parameters, machine-readable.
- Render to Word with `to_docx.py --in methods.md --out methods.docx`.

The agent should **verify the draft against the parameter table and polish the
prose** (the example in `~/Documents/DataAnalysis/.../timsTOF_Methods` shows the
target quality), resolve each `[facility default — confirm]`, and keep the
acknowledgment verbatim.
