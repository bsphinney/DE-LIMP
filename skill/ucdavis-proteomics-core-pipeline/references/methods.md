# Publication-ready Methods section (`make_methods.py`)

Generates a drop-in LC-MS/MS Methods section straight from facility raw data, plus
the correct UC Davis Proteomics Core instrument-grant acknowledgment. Can run as
part of a full analysis or **standalone** (just `--raw` at the facility data).

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 15: What it extracts vs. defaults
- 131: Instrument grant acknowledgments (verified 2026-06)
- 159: Output + workflow
- 173: 9d. The Methods section (step 9d)
<!-- /contents -->

## What it extracts vs. defaults

> Superseded by "9d. The Methods section (step 9d)" below.

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

> Superseded by "9d. The Methods section (step 9d)" below.

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

> Superseded by "9d. The Methods section (step 9d)" below.

With `--params` / `--search-prov` / `--workflow-manifest`, it adds a **Database search**
paragraph and a search-parameter table (engine, the version that ran, cleavage rule, missed
cleavages, peptide length, charge and m/z range, fixed/variable modifications, mass tolerances,
precursor FDR, library mode, MBR), read by `search_record()` — the same reader
`make_deposit.py` uses for the SDRF. A value in no record prints as
`____ [not recorded — confirm]`, never as a default.

> Superseded by "9d. The Methods section (step 9d)" below.

With `--submission <session>` (a CoreOmics submission attached by `submission_report.py`),
it adds a **Sample preparation** section from the record, through `prepared_by()`, the one
reading of the form. When the lab sent peptides, it says the submitting laboratory prepared
them and carries no Core-side placeholder; the form's own words (buffer, beads) go in a note
for the author, not into the prose, so nothing is added that the form does not state. When the
Core prepared them, the protocol is `____ [not recorded — confirm]`. `make_deposit.py` passes
the flag itself, and its PRIDE sample-processing protocol uses this section instead of its
TO-FILL.

> Superseded by "9d. The Methods section (step 9d)" below.

With `--de-dir`, it adds a **Differential expression** paragraph from the run's
`de_provenance.json` (pipeline, quantification, the q-value filters as applied, design,
contrasts, significance rule, citation). Significance is stated as run_de.R applies it —
adj.P.Val only; |log2FC| is a volcano reference line, not a filter.

> Superseded by "9d. The Methods section (step 9d)" below.

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

> Superseded by "9d. The Methods section (step 9d)" below.

| Instrument | Prefix | Acknowledgment |
|---|---|---|
| Orbitrap Fusion Lumos | `FL` | NIH S10 grant **S10OD021801** |
| Orbitrap Exploris 480 | `Ex` | NIH S10 grant **S10OD026918-01A1** |
| Bruker timsTOF (Pro 2 and HT) | — | Dr. Neil Hunter / **Howard Hughes Medical Institute** |

The web page names only the timsTOF Pro 2 for the HHMI acknowledgment. Brett confirmed on
2026-09-25 that it covers the timsTOF HT too.

> Superseded by "9d. The Methods section (step 9d)" below.

An instrument not in this registry yields a placeholder pointing at the source URL.
**Grant wording must be exact** — the script cites the verified grant numbers and
links the source page; confirm the current wording there before publishing. To add
an instrument, extend the `ACKS` table in `make_methods.py`.

## Output + workflow

> Superseded by "9d. The Methods section (step 9d)" below.

- `methods.md` — the drop-in prose (LC, MS, sequence database, [database search],
  [differential expression], parameter tables, Acknowledgments).
- `methods_params.json` — the extracted parameters, machine-readable.
- Render to Word with `to_docx.py --in methods.md --out methods.docx`.

The agent should **verify the draft against the parameter table and polish the
prose** (the example in `<your folder of instrument method files>` shows the
target quality), resolve each `[facility default — confirm]`, and keep the
acknowledgment verbatim.

## 9d. The Methods section (step 9d)

Generate a drop-in LC-MS/MS Methods section from the facility raw data and the session's own
records, with the correct UC Davis Proteomics Core instrument-grant acknowledgment:

```
python3 scripts/make_methods.py --raw /path/to/*.d \
    --fasta-meta <session>/input/search.fasta.meta.json \
    --params <session>/input/wf/<params_file> \
    --search-prov <session>/output/search/search_provenance.json \
    --workflow-manifest <session>/input/wf/workflow.manifest.json \
    --out <session>/output/methods.md --de-dir <session>/output/tables \
    [--submission <session>] \
    [--qc-bracket <session>/logs/qc_bracket.json]     # Core staff runs (step 8e)
python3 scripts/to_docx.py --in <session>/output/methods.md \
    --out <session>/output/methods.docx
```

(Step 12's `finalize` runs this for you when `output/methods.md` does not exist yet, passing
`--submission` itself when one is attached, and `--qc-bracket <session>/logs/qc_bracket.json`
when step 8e ran: an **Instrument performance** paragraph of ONE fixed sentence, the same
whatever the verdict. The verdict itself is staff-only.) `--submission` adds **Sample preparation** from
the CoreOmics record: when the lab sent peptides it says the submitting laboratory prepared
them, with no Core-side placeholder; when the Core prepared them the protocol stays a tagged
blank to fill.

What it extracts, and what it only defaults:

- It extracts the acquisition parameters from the raw metadata (Thermo by facility filename
  prefix). From a Bruker `.d` it reads the LC system, its method and run time, the ion source
  and its settings, the TIMS ramp, the dia-PASEF window scheme, the cycle time and the
  collision-energy ramp. It writes them in the order published timsTOF Methods use.
- The analytical column comes from `--lc-column`, else from HyStar's ColumnInfo, else from
  `--column-log`. That flag takes an export of STAN's column-change log, and the column in place
  at the first acquisition is used. With none of these, the facility's standard column (PepSep
  MAX C18, 10 cm × 150 µm, 1.5 µm) is printed **tagged `[facility default — confirm]`**.
  STAN keeps that log in PG Farm, which needs credentials, so the script never connects to it.
- It tags what no file records: column temperature, emitter, mobile phases, Evotip loading and
  peptide amount. STAN's 50 °C and the 20 µm emitter are given as facility defaults, tagged;
  anything else unknown is tagged `[not recorded — confirm]`.
- It builds parameter tables showing the source of each value.
- It appends the instrument's grant acknowledgment (Fusion Lumos → S10OD021801; Exploris 480 →
  S10OD026918-01A1; timsTOF → Dr. Neil Hunter / HHMI).

What each input adds:

- `--fasta-meta` writes the **Sequence database** paragraph: organism, proteome ID, UniProt
  release, database type, entry count, and the contaminant library with its citation. The
  database's build warnings are not in `methods.md`: tell the user the ones it prints under
  `database_build_warnings` (step 12's `finalize`, when it writes the Methods, prints them under
  the same key).
- `--params` / `--search-prov` / `--workflow-manifest` write the **Database search** paragraph
  from the parameters the search actually ran with: engine, the version that ran, cleavage rule,
  missed cleavages, peptide/charge/m/z ranges, modifications, mass tolerances, precursor FDR,
  library mode, and MBR.
- `--de-dir` writes the **Differential expression** paragraph from `de_provenance.json`. It states
  significance as run_de.R applies it: adj.P.Val only, with |log2FC| a volcano reference line,
  never presented as a filter.

A value in no record prints as `____ [not recorded — confirm]`, never as a guess. If the raw files
are not readable where you run it, add `--instrument "<name>" --acquisition <DIA|DDA>` from the
workflow manifest: the Methods are then written from the session record, with the acquisition
values blank and tagged. **Verify the draft against the parameter tables and polish the prose;
keep the acknowledgment exact** (confirm the wording at the source URL). This can also be run
standalone — point `--raw` at facility data; no search or DE is needed. → detail:
`references/methods.md`.
