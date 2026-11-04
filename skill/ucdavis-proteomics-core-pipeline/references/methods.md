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
  4. Otherwise, **for a timsTOF `.d` only**, the facility's standard column, from STAN's column
     catalogue (PepSep MAX C18, 10 cm × 150 µm, 1.5 µm, part 1893483), **tagged
     `[facility default — confirm]`**. It is the Core timsTOF HT's column, so any other run (a
     Thermo `.raw`, an mzML, data from an instrument the Core does not have) gets
     `an analytical column ____ [not recorded — confirm]`: a Thermo cohort once got the Bruker
     PepSep column (2.14.2).
- **Never in a .d, always tagged** (`make_methods.column_record` and `lc_paragraph`, one rule):
  - Column temperature. **For a timsTOF `.d` only**, STAN's operator-reported 50 °C, tagged
    `[facility default — confirm]`; otherwise `____ °C [not recorded — confirm]`.
  - The emitter. **For a timsTOF `.d` only**, the default 20 µm CaptiveSpray emitter, tagged
    `[facility default — confirm]`.
  - Mobile phases (0.1% formic acid in water / in acetonitrile) and the nanospray source: the
    facility's, tagged `[facility default — confirm]`, only for a run on a Core instrument (a
    `.d` or `.raw`); an mzML or data from elsewhere gets `____ [not recorded — confirm]`.
  - Evotip type, loading and peptide amount.
  - The %B gradient of an Evosep method. Its fixed, named method is given instead.
  - A non-Evosep LC's gradient table, which is not parsed. It is tagged
    `[not recorded — confirm]`.
- **Tags.** Every one is a question for the user or a record, and **is never replaced with a
  guess** (DE-LIMP rule 2): resolve it from a record or the user's answer, or leave it tagged.
  - `[facility default — confirm]`: a Core default (above), what is usually fitted, not a record
    of this run.
  - `[not recorded — confirm]`: a value no record holds (printed as `____`), never a default.
  - `[raw file not readable here — confirm]`: a `.d` that could not be read, named on stderr
    and in the note under Mass spectrometry; when it is the only run, its values carry this tag.
  - `[vendor — confirm]`: the instrument's name matches neither vendor's words, or both
    (`make_methods.VENDOR_TAG`): an mzML from an unnamed instrument is not called Thermo.
  - `[acquired at the UC Davis Proteomics Core? — confirm]`: after the acknowledgment, when the
    instrument's name is in the Core's registry but the data are not a Core `.d`/`.raw` (an
    mzML, outside data; `ACK_OUTSIDE_TAG`). Keep the acknowledgment only if the Core acquired it.
  - `[DDA/DIA — confirm]`: no record gives the acquisition mode (order below).
  - `[Evotip type, loading protocol and peptide amount — confirm]`, and
    `[not extracted by this script — take from the MS method; confirm]` for the ddaPASEF
    precursor-selection settings, which are not extracted.

With `--fasta-meta` (fetch_fasta.py's `<fasta>.meta.json`), it adds the **Sequence database**
paragraph. Sequences added with `fetch_fasta.py --add-fasta` are named one by one only when there
are at most 10 (`ADDED_SEQUENCES_NAMED_MAX`: a bait, a tag); more are given as a count with the
file name and SHA-256, and so are the contaminant entries removed or kept beside them. The
database's build warnings never go into `methods.md`: they are printed on stderr and listed under
`database_build_warnings` in the JSON the script prints, for you to tell the user and resolve.
When `session.py finalize` writes the Methods itself (through `make_deposit.py`), it prints them
on stderr and under the same key in its own JSON; `MANIFEST.txt`, which travels in the zip, gets
only their count.

No file path reaches `methods.md`: an absolute path from any record (`/quobyte/...`, `/nfs/...`,
`/home/...`, `~/...`, `C:\...`) is cut to its file name as the file is written, and stderr shows
each cut (`<path> -> <name>`). A path starts the text or follows a space or an opening quote or
bracket, and starts with a file-system root (`make_methods.FS_ROOTS`), `~/` or a drive letter, or
has three file-name segments, so prose such as `S-Trap (ProtiFi)/Trypsin/Lys-C` or
`0.1%/80%/20%` is left alone.

With `--params` / `--search-prov` / `--workflow-manifest`, it adds a **Database search**
paragraph and a search-parameter table (engine, the version that ran, cleavage rule, missed
cleavages, peptide length, charge and m/z range, fixed/variable modifications, mass tolerances,
precursor FDR, library mode, MBR), read by `search_record()` — the same reader
`make_deposit.py` uses for the SDRF. A value in no record prints as
`____ [not recorded — confirm]`, never as a default.

**Which search, and what its record says (2.14.2):**
- `--session <S>` (with no `--search-prov`) takes the session's own search record,
  `session.search_choice`: wherever the search went (`output/search_sage/` …), and with several,
  the one the DE read. Several records and none chosen: the paragraph is **tagged blanks**, with
  why, never one search's values; `--search-prov` given more than once does the same. A
  `--search-prov` file that is not there is refused. Choose a search as SKILL.md step 12b says.
- The acquisition mode comes from `--acquisition`, then the workflow manifest's `acquisition`,
  then the session's `input/acquisition.json` (step 2); with none, `[DDA/DIA — confirm]`.
- DIA-NN's mass accuracy is "optimised automatically per run" only when the record says DIA-NN
  chose it; measured-and-pinned says so; with no record of the run it is
  `[not recorded — confirm]` (a search pinned at 10/20 ppm was once called optimised).
- Sage's settings are read from its `results.json` (the search as it ran), FragPipe's from
  `fragger.params` and the saved `fragpipe.workflow` (the Philosopher filter as run, the
  database MSFragger read).
- An identifications-only search (Sage `--ids-only`, a FragPipe workflow with no
  quantification) says so, with its FDR rule; a timsTOF `.d` searched by Sage says whether
  Sage read it directly (identifications only) or it was converted to mzML first (with the
  msconvert options); a DIA-NN `--dda` search whose second pass crashed says the first-pass
  report was used and that DIA-NN's protein MaxLFQ quantities are absent.
- The **Differential expression** paragraph states the between-run normalisation the DE record
  names (`de_provenance.json` `normalisation`; on the MaxLFQ route, the `--between-run` step,
  `references/de-analysis.md`).

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
target quality), resolve every `[… — confirm]` tag in the list above from a record or the
user's answer (never with a guess: an unresolved tag stays), and keep the acknowledgment
verbatim.
