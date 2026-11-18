# Publication-ready Methods section (`make_methods.py`)

Generates a drop-in LC-MS/MS Methods section straight from facility raw data, plus
the correct UC Davis Proteomics Core instrument-grant acknowledgment. Can run as
part of a full analysis or **standalone** (just `--raw` at the facility data).

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 15: What it extracts vs. defaults
- 135: Instrument grant acknowledgments (grant numbers checked 2026-06)
- 172: Output + workflow
- 185: 9d. The Methods section (step 9d)
<!-- /contents -->

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

## Instrument grant acknowledgments (grant numbers checked 2026-06)
The grant numbers come from
https://proteomics.ucdavis.edu/instrument-grant-acknowledgments. That page lists the grants and
names, not full sentences, so the sentences around the numbers are this skill's own wording.

**Only for data the Core acquired.** The Core's acknowledgment is printed only when the run's
record shows the Core acquired the data (`submission_report.core_acquisition`, one test): a
CoreOmics submission attached to the session, or announced at init (`session.py init
--submission`, record not attached yet), Core staff's client-work mark (`session.py
core-client`), the Core's internal work (`session.py init --internal`), or `--acquired-at-core`.
An instrument name that matches a Core instrument is not enough, because the skill is used by
other labs. Without a record the Acknowledgments section holds a tagged blank, and the Core's
LC defaults (column, emitter, column temperature, mobile phases) are tagged blanks too. A
standalone Methods (`--raw`, no session) has no record: **the agent asks the user whether the
Core acquired the data** and passes `--acquired-at-core` when they say yes (SKILL.md 9d). The
reason is recorded in `<out>_params.json` `core_acquisition`.

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
**Grant numbers must be as on the page** — the script prints the grant numbers listed there and
links the source page; confirm the current wording there before publishing. To add
an instrument, extend the `ACKS` table in `make_methods.py`.

## Output + workflow

- `methods.md` — the drop-in prose (LC, MS, sequence database, [database search],
  [differential expression], parameter tables, Acknowledgments).
- `methods_params.json` — the extracted parameters, machine-readable.
- Render to Word with `to_docx.py --in methods.md --out methods.docx`.

The agent should **verify the draft against the parameter table and polish the
prose** (the example in `<your folder of instrument method files>` shows the
target quality), resolve every `[… — confirm]` tag in the list above from a record or the
user's answer (never with a guess: an unresolved tag stays), and keep the acknowledgment
verbatim.

## 9d. The Methods section (step 9d)

Generate a drop-in LC-MS/MS Methods section from the facility raw data and the session's own
records, with the correct UC Davis Proteomics Core instrument-grant acknowledgment:

```
python3 scripts/make_methods.py --raw /path/to/*.d \
    --fasta-meta <session>/input/search.fasta.meta.json \
    --params <session>/input/wf/<params_file> \
    --session <session> \
    --workflow-manifest <session>/input/wf/workflow.manifest.json \
    --out <session>/output/methods.md --de-dir <session>/output/tables \
    [--submission <session>] \
    [--qc-bracket <session>/logs/qc_bracket.json]     # Core staff runs (step 8e)
python3 scripts/to_docx.py --in <session>/output/methods.md \
    --out <session>/output/methods.docx
```

**Which search, and what its record says (2.14.2):**
- `--session <S>` (with no `--search-prov`) takes the session's own search record,
  `session.search_choice`: wherever the search went (`output/search_sage/` …), and with several,
  the one the DE read. Several records and none chosen: the paragraph is **tagged blanks**, with
  why, never one search's values. `--search-prov <file>` names one by hand; given more than once
  it does the same as none chosen, and a file that is not there is refused. To choose one for
  `--session`, run `session.py choose-search --dir <session> --folder <its folder under output/>`
  (nothing is renamed -- the groups of a FAIMS split are valid searches), or rename a replaced one
  to `output/superseded_<name>/` before finalize (a folder whose name holds "superseded" is not a
  candidate); a Core delivery: core-submissions.md
  "12b.", "Several searches".
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
  at the first acquisition is used. With none of these, **for a timsTOF `.d` only**, the
  facility's standard column (PepSep MAX C18, 10 cm × 150 µm, 1.5 µm) is printed **tagged
  `[facility default — confirm]`**; any other run (a Thermo `.raw`, an mzML, outside data) gets
  `____ [not recorded — confirm]`.
  STAN keeps that log in PG Farm, which needs credentials, so the script never connects to it.
- It tags what no file records: column temperature, emitter, mobile phases, Evotip loading and
  peptide amount. STAN's 50 °C and the 20 µm emitter are given as facility defaults, tagged,
  **for a timsTOF `.d` only**; the facility's mobile phases only for a Core `.d` or `.raw`;
  anything else unknown is tagged `[not recorded — confirm]`.
- **The tags you will meet**, each resolved from a record or the user's answer and **never
  replaced with a guess** (an unresolved tag stays): `[facility default — confirm]`,
  `[not recorded — confirm]`, `[raw file not readable here — confirm]`, `[DDA/DIA — confirm]`,
  `[vendor — confirm]` (the instrument name names no one vendor: an mzML is not called Thermo),
  and `[acquired at the UC Davis Proteomics Core? — confirm]` after the acknowledgment (a Core
  instrument's name on data that are not a Core `.d`/`.raw`; keep the acknowledgment only if the
  Core acquired it). → the full list: **Tags** above.
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
leave the grant numbers as printed** (confirm the wording at the source URL). This can also be run
standalone — point `--raw` at facility data; no search or DE is needed. The sections above give
the detail of each part.
