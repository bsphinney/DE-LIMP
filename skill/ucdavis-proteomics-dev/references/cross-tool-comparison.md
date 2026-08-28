# Comparing two search engines — the method, not just the script

`make_comparison_report.py` regenerates the artifact. This file records **how to reach a
defensible answer**, and the interpretation patterns that a bake-off keeps producing.
Written after the first real DIA-NN vs FragPipe/diaTracer comparison (9 mouse dia-PASEF
files); every pattern below was actually observed in that data.

## Get the design right before computing anything

**Run both engines through the SAME downstream pipeline.** Identical `run_de.R --method`,
identical contrasts, identical thresholds. Otherwise you are comparing search × DE
method and cannot attribute the difference. `run_de.R --format tsv` reads a FragPipe
`report.tsv` directly, and a `conditions.csv` written for the DIA-NN route works
unchanged because the `Run` column carries the original `.d` basenames.

**Name every uncontrolled difference out loud, in the report, not just in your head.**
The big one for this pair: FragPipe 24 bundles **DIA-NN 1.8.2 beta 8** while the
`diann_*` workflows pin **2.x**. Both routes end in DIA-NN, so part of any gap is two
major versions of quant engine, not the approach under test. A head-to-head that omits
this silently attributes version lag to the method. Control it with `--config-diann`
pointing at an external 2.x, or state it.

**Compare at both levels; they answer different questions.**

| Tool | Question it answers |
|---|---|
| `compare_searches.py` | What did each engine *find*, in how many runs, how reproducibly |
| `compare_analyses.R` | Which proteins actually come out *changed* |
| `QC_detected_vs_inferred.csv` (both sessions) | How much of each matrix was *measured* vs modelled |

A search-level win does not imply a DE-level win. In the reference case DIA-NN found 20%
more protein groups yet **lost** the subtle contrast on DE count.

## Interpretation patterns that recur

### More precursors, fewer proteins → protein inference, not sensitivity
diaTracer identified up to **25% more precursors** than DIA-NN on the weakest samples
while returning **fewer protein groups in every run**. That combination is not "detected
less" — it is a parsimony difference (ProteinProphet collapsing shared peptides vs
DIA-NN's own grouping). Always check precursor counts alongside protein counts before
calling one engine less sensitive. Reporting only the protein number inverts the story.

### logFC correlation vs intensity correlation → how far to trust magnitude
Intensity agreement on shared proteins was r = **0.879**; fold-change agreement was only
**0.637–0.718**. Ratios amplify disagreement because the errors do not cancel when you
divide. When intensity r is high and logFC r is much lower, the honest statement is:
**direction is trustworthy, magnitude is engine-dependent.** Report both correlations —
quoting only the intensity one overstates agreement.

### Effect size determines how much the engine matters
On a large contrast (Urea vs Beads) 708 of 1,866 DE calls overlapped. On a subtle one
(S-Trap vs Beads) only **91 of 382** did — two-thirds of that DE list depended on which
engine ran. **The weaker the effect you are chasing, the more the search engine
determines your answer.** Say this explicitly when the study's contrast is subtle; it is
the difference between "engine choice is a detail" and "engine choice is the result".

### Read the 3×3 matrix by cell type
- **Diagonal** (Up/Up, Down/Down, NS/NS) — agreement.
- **NS ↔ significant** off-diagonals — *detection* differences. Almost always where the
  disagreement lives, and the benign kind.
- **Up ↔ Down corners** — genuine direction disagreement. Should be near-empty; 9 of 708
  in the reference case. A populated corner is a red flag, not a curiosity.

### Completeness advantage tends to track sample quality
diaTracer's detected fraction was 5–7 points higher on all three urea (worst) runs and
*lower* on the strongest S-Trap runs. An empirical library built from the data helps most
where signal is scarce. Break the QC panel down by sample quality rather than reporting a
single mean — the mean hides the pattern that explains the mechanism.

### Protein-group identifiers are not perfectly comparable across tools
Different parsimony rules mean a protein can appear "unique to engine X" purely because
it was merged into a differently-named group. Some fraction of every unique list is
grouping, not detection. **Warn about this before anyone reads biology into a unique
set**, and check per-protein claims at peptide level first.

## Building the report

Follow the structure DE-LIMP's Run Comparator prompt defines
(`docs/AI_PROMPTS.md` §1) — `make_comparison_report.py` emits it as a skeleton:

1 Factual Observations · 2 Sources of Disagreement · 3 Case for A · 4 Case for B ·
5 Settings Audit · 6 Concordant Proteins · 7 Synthesis · 8 Recommended Follow-ups

Its guidelines are load-bearing, not decoration:

- **Neither tool is inherently superior.** Sections 3 and 4 exist so you argue *both*
  sides properly before section 7 weighs them. Write 3 and 4 before you know your verdict.
- **Every claim carries a number from the data or a literature citation.**
- **Label speculation** — "one possible explanation is…", never stated as fact.
- **If the evidence does not clearly favour one tool, say so.** Do not force a
  recommendation. The reference case concluded *complementary, not ranked*, and that was
  the correct answer, not a hedge.
- Section 8 wants **one follow-up probing each pipeline** — not two aimed at the same tool.

**The generator writes the numbers; you write the judgement.** Fill each `TODO(model)`
block from the data. Never invent a figure to complete a sentence.

## Charts: load the `dataviz` skill first

The HTML template's palette and mark choices are not taste — they came from the dataviz
method, and a future chart must re-derive them the same way:

- **Run the validator, don't eyeball it.**
  `node dataviz/scripts/validate_palette.js "<hex,hex>" --mode light` and again
  `--mode dark`. The two-series pair shipped in `make_comparison_report.py`
  (`#2a78d6` / `#eb6834` light, `#3987e5` / `#d95926` dark) passes all six checks in both
  modes. Changing a colour obliges re-running it.
- **Form follows the data's job.** Grouped bars for magnitude across categories;
  **dumbbell** for two paired values per sample (it shows the *gap*, which is the point,
  and needs one row per sample rather than two); a **sequential single-hue ramp** for
  counts in the 3×3 matrix, with agreement cells on the accent hue and direction-flip
  cells on the warning hue so the eye lands on the cell that matters.
- **Hover by default.** An HTML chart is interactive; every mark gets a tooltip naming
  the run, the value, and what the cell *means* ("both agree", "OPPOSITE DIRECTION").
- **Theme-aware under both scopes** — the `prefers-color-scheme` media query and the
  `data-theme` attribute, so a viewer's toggle wins either way. Redraw the SVGs on toggle;
  CSS variables alone will not repaint canvas-computed geometry.
- **Colour is never the only encoding** — legend always present for two series, direct
  value labels beside the marks.
- Self-contained: no CDN, no external fonts, inline SVG. It has to open from a file path
  on a cluster with no browser network access.

## Deliver both formats

HTML for exploration, `.docx` for circulation — `to_docx.py` handles the conversion and
resolves relative figure paths, so images actually embed rather than silently vanishing.
