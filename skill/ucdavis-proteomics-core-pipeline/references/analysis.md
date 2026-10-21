# Data analysis & reporting

The skill produces a written biological + QC interpretation of the results, plus a
catalog of every output file. The analysis prompt is ported from DE-LIMP's
"Export Prompt for Claude" (`R/server_ai.R`) so the output matches what the
DE-LIMP app produces — but here the agent IS Claude, so it does the analysis
inline instead of shipping the prompt to an external API.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 21: The flow
- 55: The CoreOmics submission (Core data)
- 68: Report sections (faithful to DE-LIMP's export prompt)
- 75: Pipeline self-description (DE-LIMP rule #1)
- 85: No fabrication (DE-LIMP rule, and the viz-data feedback)
- 94: Output-files report
- 100: Less-used flags
- 115: 9. Writing the report (step 9)
<!-- /contents -->

## The flow

> Superseded by "9. Writing the report (step 9)" below.

1. `analysis_prompt.py` writes `ANALYSIS_PROMPT.md` — the brief, parameterized by
   what's present (QC, GSEA) and by the actual engine + DE method.
2. **The agent reads the brief and the data files and writes `AI_Analysis_Report.md`.**
   It computes significant proteins, up/down splits and cross-comparison overlaps from the
   CSVs, and CVs from measured values only — citing specific proteins, never fabricating.
   Hits are tiered by per-group detection: the DE tables' `Evidence` (measured in both /
   partly inferred or missing / presence call, as `de_provenance.json` defines it), never
   by PropObs. A `not measured in either group` row was not tested in that comparison (no
   p-value): it is never a hit and not counted as tested. The brief lists every figure in `figures/figures.json` (the 2.8.0 object or an
   older bare list) and names each one make_figures.R could not draw, with its reason.
3. `make_analysis_html.py` renders it into **`Analysis_Report.html`** — the report of
   record: one self-contained page (figures inlined). Only figures in this run's
   `figures.json` are embedded. What must not depend on the writer is added by the page
   itself: Results at a glance, fixed callouts (inferred values, a contaminant-database
   risk, figures that could not be drawn), the top-protein tables and, as the last section,
   *Appendix: p-value calibration* (`qc_pvalue_panel.png`). The instrument's QC around the
   project (step 8e) is **staff-only** and never appears in the report: see
   `references/qc-bracket.md`. No Word copy of the report is made
   any more (Word mangled the figures, 2026-09-24); the Methods stay in Word
   (`methods.docx`). An older session's `AI_Analysis_Report.docx` is left in place.
   The same call writes two twins beside it:
   - **`Analysis_Report.md`** — the report as plain text with each figure's caption and
     numbers written out, for NotebookLM and other AI notebooks. `--md-out <path>` puts it
     elsewhere (default: `--out` with `.md`).
   - **`Analysis_Report.pdf`** — the HTML printed by a headless Chrome, Chromium or Edge when
     one is installed (`html_to_pdf.py`). `--no-pdf` skips it; with no browser there is no PDF
     and the output says so, and `session.py finalize` prints it later where a browser exists.
4. `make_report.py` writes `OUTPUT_FILES.md` — every output file, its size, and a
   plain-language description, grouped by purpose.

## The CoreOmics submission (Core data)

> Superseded by "9. Writing the report (step 9)" below.

When `submission_report.py attach` has stored the submission in the session,
`make_analysis_html.py` opens the report with a **Submission** section by default: the PROT id
linked to CoreOmics, PI, submitter, date, organism and UniProt as specified, the description
as written, experiment type, proteins or peptides and who prepared them, buffer, beads,
normalisation, the analysis requested, the sample sheet, and the record's Data Quality Notes.
The record is allowlisted — no email, phone, billing or internal field can reach the page.
`analysis_prompt.py --submission <session>` quotes the same record in the brief, so the
report describes the samples in the submitter's words and adds nothing they did not state.

## Report sections (faithful to DE-LIMP's export prompt)
Overview · QC Assessment (if QC present) · Key Findings Per Comparison ·
Specificity (pull-down designs) or Cross-Comparison Biomarkers + High-Confidence
Findings · Pathway/GSEA (if present) · Biological Interpretation · How This Analysis
Works (LC-MS/MS, DIA/DDA, the engine, the stats framework, key terms — written for a
biologist with no MS background) · Methods & Reproducibility · Next steps.

## Pipeline self-description (DE-LIMP rule #1)

> Superseded by "9. Writing the report (step 9)" below.

The brief takes the pipeline label, quantification, DE engine, missing-value
policy, and citation from `de_results/de_provenance.json` (written by `run_de.R`),
not from hardcoded strings. The educational + stats sections are chosen by engine
(DIA-NN vs Sage) and method (DPC-Quant vs MaxLFQ). So the report is correct for
whichever pipeline actually ran, and never contradicts `methods.txt`.

## No fabrication (DE-LIMP rule, and the viz-data feedback)

> Superseded by "9. Writing the report (step 9)" below.

Every claim must trace to a value in the attached files or to established biology.
Don't invent fold-changes, pathways, p-values, or citations. If something isn't in
the data (e.g. no GSEA was run), the corresponding section is simply omitted by the
prompt generator rather than filled with guesses.

## Output-files report
`OUTPUT_FILES.md` is generated by matching each file against a description table
(search output, DE tables, reproducibility bundle, inputs, analysis report).
Unrecognized files are listed under "Other" with a note — never silently dropped
(DE-LIMP rule #4, the MANIFEST discipline).

## Less-used flags
- `make_analysis_html.py` takes its inputs from `--session`; without one, name them:
  `--report AI_Analysis_Report.md --figures <dir> --tables <dir> [--quality SAMPLE_QUALITY.md]
  [--audit AUDIT.md]`. `--title` overrides the page title (default: the report's own `#`
  heading); `--adjp` applies only when the tables carry no `de_provenance.json`.
- `analysis_prompt.py --qc <QC_Metrics.csv> --gsea <GSEA_Results.csv>` add those sections to
  the brief; `--report-out` names the report the brief asks for (default `AI_Analysis_Report.md`).
  A `--gsea` file from gene-permutation GSEA (gseGO / fgsea) treats a set's proteins as
  independent, the flaw that keeps `geneSetTest` out of `run_sets.R`: the brief presents it as
  exploratory only, and never beside the protein-set tests as a second test
  (`references/set-tests.md`).
- `make_report.py --root <dir>`: the folder paths in OUTPUT_FILES.md are shown relative to.
- `html_to_pdf.py --html <report.html> [--pdf <out.pdf>] [--timeout 120]`: print a report by
  hand (default: beside it, `.pdf`).

## 9. Writing the report (step 9)

Generate the analysis brief (it lists the figures to embed), then **do the analysis
yourself** — you are the consultant the brief addresses:

```
python3 scripts/analysis_prompt.py --out ANALYSIS_PROMPT.md \
  --de-dir ./de_results --report ./search_out/report.parquet \
  --conditions ./conditions.csv --figures-dir ./figures [--qc ./QC_Metrics.csv] \
  --engine <engine> --acquisition <DIA|DDA> --instrument "<name>" \
  --workflow-manifest ./wf/workflow.manifest.json --session <session> [--submission <session>]
```

Always pass `--session <session>`: the brief then states the sample type as the user gave it
(`sample_type.py`), or **NOT STATED** — in which case the report names no tissue, cell type or
biofluid, and never takes one from file names.
With a CoreOmics submission attached (step 3b), pass `--submission <session>`: the brief then
quotes the record and its Data Quality Notes. Describe the samples in the submitter's words
and add nothing they did not state ("cross-linked" must not become "chemically cross-linked").
Then **read `ANALYSIS_PROMPT.md` and every data file + figure it lists, and write a
complete `AI_Analysis_Report.md`** with ALL its OUTPUT sections (Overview, QC, Key
Findings Per Comparison, Cross-Comparison Biomarkers, High-Confidence Biomarkers,
Pathway/GSEA if present, Biological Interpretation, How This Analysis Works, Methods
& Reproducibility), a **Data
Quality Notes** section (the 5-part anomalies from steps 8c/8d + `SAMPLE_QUALITY.md`,
even if "nothing anomalous observed"), and — if a review panel ran (8d) — an **Expert
Review Notes** section. **Do not write an "Audit & caveats" section:** `make_analysis_html.py`
adds it as a plain-language table from `AUDIT.json`, one row per check (what was checked, the
result, what it means for the results), in words written beside each check
(`audit_results.py`; a pipeline's own checks carry theirs). The technical `AUDIT.md` (check
codes, thresholds, file names) stays the full record. A staff report: the raw AUDIT.md in the
client report ("implausibly high ... not real biology") was unreadable to staff and clients
alike. Discuss a WARN or FAIL that bears on a finding next to that finding and in Data Quality
Notes. A report's own "Audit & caveats" section still replaces the table, so write one only
when staff ask for it, and then in plain words. **Embed
every figure** (`![caption](figures/<file>.png)`) with an expert interpretation of
what it shows for THIS data. Compute significant proteins, up/down splits,
cross-comparison overlaps, and lowest-CV proteins directly from the CSVs — cite
specific proteins, never fabricate. Make it thorough and expert, like the DE-LIMP
AI export. The brief takes its pipeline description from `de_provenance.json`, so
the report stays correct for whichever engine/method ran.

Then produce the report of record — ONE self-contained HTML page:

```
# QC panels + figures + text in ONE file: the report of record.
python3 scripts/make_analysis_html.py --session <session> \
    --title "<study name>" --out <session>/output/Analysis_Report.html
```

With a CoreOmics submission attached (step 3b) the page opens with its **Submission** section
and the header names it. `--submission <record file or fetch folder>` shows one that is not
attached; a bare `PROT_####` is not a record (the page says it could not be read) — record
the facts with `submission_report.py attach --given` instead.

**Do not make a Word copy of the report** (Brett, 2026-09-24: Word mangled the figures —
the heatmap did not display in `AI_Analysis_Report.docx`). The Methods stay in Word
(step 9d: people paste and edit them in manuscripts). An older session that already has an
`AI_Analysis_Report.docx` keeps it; just do not point anyone at it.

**Point the user at the HTML first.** Every figure is inlined as a data URI, so it is
one file that opens by double-clicking in any browser on Windows, macOS or Linux, with
no network, no Word, and nothing installed. That matters because the alternative —
a report plus a `figures/` folder — silently loses every image the moment someone
copies just the report, which is exactly what people do. Tell them the filename and
that it is the whole report.

Each `![caption](figures/x.png)` in `AI_Analysis_Report.md` becomes a numbered, embedded
figure at that spot in the text (caption from `figures.json` when there is one). Images the
report does not reference are left out, with one `WARNING` naming them — a redrawn plot's old
copy is not shipped; a referenced image that is missing shows as a visible "figure missing"
note. So reference every figure you want, where you discuss it. Only with no report at all
does the page fall back to galleries of `figures.json`'s figures. The look (header band with
the study facts, stat tiles, contents rail, figure cards with click-to-enlarge, callouts for
the audit / data-quality / expert-review sections, dark mode, print styles) lives in
`scripts/report_style.py`, shared by the skill's HTML pages — restyle there, not per page.

The same run also writes **`Analysis_Report.md`** beside the HTML: the full report as plain
text, for NotebookLM or other AI notebooks. It is rendered from the same assembled sections
(same headings, callouts as labelled blockquotes), and because those tools read text, not
images, each figure carries its caption plus a data summary from the tables (volcano: counts
and the top 5 by adj.P), with a top-20 protein table per contrast. Offer it to users who
want to "ask questions of the report".

And **`Analysis_Report.pdf`**: the same HTML printed through its print stylesheet by a
headless Chrome/Chromium/Edge (`scripts/html_to_pdf.py`) — the full report with figures, for
NotebookLM (which reads a PDF's images too) or printing. No browser (usual on HIVE) is not an
error: the run says so, and `session.py finalize` — on the laptop in hive_remote — makes it
then, recording `[OK]` or `[INFO]` (with the by-hand route: open the HTML, Print, Save as PDF)
in MANIFEST.txt.

The page puts the **QC panels above the results** on purpose: a volcano plot is equally
persuasive whether or not the run was any good, so a reader who meets the biology first
has already formed a conclusion before seeing the evidence about whether to trust it.
Do not "improve" the order. It also counts significant proteins straight from the DE
CSVs rather than restating the prose, so a disagreement between the two is visible
rather than hidden.

Core staff: board.md §3 row "results ready", and again each time a result is redone.

Once the report is final, the audio discussion (step 9e): **made by default for a CoreOmics
submission**; offered in one line for internal Core work, and to anyone else when the results go
to a collaborator or the lab.
