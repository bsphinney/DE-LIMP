# "Show me the data for <protein>": a screen that explains where a protein's number comes from

Design, 2026-10-26. Brainstormed with Brett Phinney (UC Davis Proteomics Core) on 2026-10-23 and 2026-10-26.
Part of the next version of the UC Davis proteomics skill ("agent-built screens"). Builds on the released
skill 2.14.0 (tag `skill-v2.14.0`).

Working mockup, built on real data: `.superpowers/brainstorm/14251-1793038039/content/quant-trace-v32.html`
(serve with the superpowers visual companion, or open the file after replacing nothing: it is self-contained).
Its source, data and every extraction script: `.superpowers/brainstorm/14251-1793038039/source/` (see its README).
The mockup is the reference for layout and wording; this document is the reference for behaviour and rules.

Status (2.15): part 1 is built (`scripts/make_screen_data.R`, `scripts/make_protein_screen.py`; how to run it and
which analyses it covers: `skill/ucdavis-proteomics-core-pipeline/references/protein-screen.md`). The mockup and the
part-1 build plan named below are working files outside this repository (the mockup on the Core director's Mac, the
plan on the branch the screen was built on); they are kept as named for the record.

---

## 1. Intent

A collaborator, or a Core staff member, types **"show me the data for BASP1"** (any protein) to their AI agent.
The agent builds one screen that answers two questions:

1. **Who is this protein?** Name, gene, function, and whether the Core has seen it before.
2. **Where did its number come from?** From the fold change in the results table, down through each sample's
   protein value, each precursor's value, DIA-NN's raw quantity and confidence checks, to the actual peaks.

Brett's words: "I basically wanted a page to help me understand where all these numbers are coming from."
Success is that a reader can explain a protein's result back in their own words, including what it rests on
and how much to trust it.

The screen is built **on demand by the reader's own AI agent** from the delivered results folder (decided
2026-10-23: "ceiling only, no floor"). The static analysis report stays the record of the analysis; this screen
does not replace it, and no new built-in interactive report is added to the deliverable.

## 2. Goals, non-goals, anti-goals

**Goals**
- Trace every number on the screen to a named file and column, a line of code, or a cited paper.
- Show what a protein value rests on: which precursors, how much each one carries, which samples were inferred.
- Show the raw evidence: per-sample XICs, so a reader can see whether peaks are "intense and nice looking or
  very noisy" (Brett).
- Flag results whose conclusion depends on how the data were processed.
- Explain every term on the page, in plain language, on the page.

**Non-goals**
- Re-analysing the data, or offering live FDR or fold-change sliders (rejected 10-23 as threshold shopping).
- Replacing the static report or the results tables.
- Explaining DIA-NN's internal QuantUMS weights: DIA-NN does not report them. The screen says so.

**Anti-goals (the screen fails if it does any of these, even if it "works")**
- Shows a number that was not computed from the delivered files, UniProt or FRAN. No invented, estimated or
  "illustrative" values. Missing data is shown as missing.
- Calls a not-found peptide "inferred", or shows limpa's per-precursor guesses as if they were data.
- Reads as "limpa is wrong" or "the result is unreliable" when two pipelines merely differ in certainty.
- Answers biology from model memory. Protein descriptions come from UniProt only (see the 2026-09 fix
  "stop asking the model to recall protein biology", commit b0683c04).

## 3. The screen, top to bottom

What is visible by default, and what is folded behind a "Show details" button. The mockup implements all of it.

1. **The request**, as a chat bubble ("show me the data for BASP1"), with a one-line note that every block reads
   only this project's files plus UniProt and FRAN.
2. **Protein card.**
   - Gene, protein name, accession, UniProt ID, reviewed status, organism, length. Links: UniProt, FRAN protein
     page, FRAN gene page.
   - *What UniProt says*: function and subcellular location, verbatim from UniProt. If UniProt has no function
     text, say so (BASP1 has none).
   - *In this experiment*: fold change, adjusted p, whether it is in the report's changed list; found in n/3 per
     group; precursors and peptides used and dropped; a sequence-coverage bar (clickable peptides);
     a shared-peptide statement (section 5.6).
   - *In FRAN*: runs and searches it was seen in, FRAN's peptide count and coverage, how many of this run's
     peptides FRAN has seen, and the note that FRAN's fragment traces are cross-run averages.
   - *Checks* strip (section 6).
   - Sources line: UniProt entry version and date, FRAN fetch date, project files.
3. **The short answer**: 3 to 6 sentences written from the protein's data, each ending in a link that jumps to
   (and if needed unfolds) the part of the page that shows it. Sentence templates and their exact conditions are
   in the mockup's `storyHTML()`. Covered: the number and its call; shared peptides (if all peptides are shared);
   what the value is built from and the largest single-peptide share; how limma weighted the samples; what the
   MaxLFQ route says; DIA-NN quality scores and surprising gaps.
4. **How a protein's number is made**: one plain-language paragraph (about 250 words) covering the whole chain without dropping steps: DIA-NN and match-between-runs, QuantUMS, normalisation, the confidence checks and "not found"; limpa's model, equal weights, the detection curve, inferred values; limma's weights, borrowing across proteins and multiple-testing correction; the two check routes and the flags. Text in the mockup's `howItWorksHTML()`. Requested by Brett ("easy to understand yet does not throw away any complexity").
   Under it, folded, **the picture**: the same chain drawn with the selected protein's own numbers (DIA-NN panel
   following the selected precursor in the selected sample; limpa's precursor × sample grid and detection curve;
   limma's weights as dot sizes, borrowed degrees of freedom and BH; the two check routes along the bottom, each
   meeting limpa's result at its own flag). Mockup: `pictureHTML()`. Requested by Brett.
5. **The path strip**: one row of boxes, "the number ← this sample's protein value ← this precursor ← DIA-NN's
   raw quantity × normalisation ← the peaks", for the currently selected sample and precursor. Each box jumps to
   its step.
6. **Step 1, the fold change (limma).** Visible: the change, fold, adjusted p; a dot plot of the six protein
   values with error bars (inferred values hollow and labelled); one sentence on whether limma's weights changed
   the answer versus a plain average. Folded: *how limma weighs each sample* (section 5.4) and *two routes to a
   fold change* (section 5.5).
7. **Step 2, each sample's protein value (limpa).** Visible: the per-peptide share bars, one 100% bar per sample
   (section 5.3), and the "what decides a share" note. Folded: the precursor-by-precursor view (one line per
   precursor across samples, the protein line, the detection curve with this sample's precursors on it, and a
   receipt for the selected sample); the Sankey of the same shares, with intensity and count views; the six
   samples' run quality (section 5.8).
8. **Step 3, each precursor's value (DIA-NN).** Precursor chips (found counts; dropped precursors struck
   through); FRAN and shared-peptide lines for the selected peptide; a table per sample: DIA-NN quantity,
   × normalisation factor, = normalised, log2 used by limpa, QuantUMS quality score, detection chance
   (section 5.7), the six confidence checks with the one closest to its cutoff, and used yes/no. A sample with
   no report row says DIA-NN looked and reported no value at its cutoff, and points to its trace if one exists.
9. **Step 4, the peaks (XICs).** One panel per sample for the selected precursor: fragment traces (y ions solid,
   b ions dashed, darker = taller) above, MS1 below, DIA-NN's reported apex as a dashed line. Toggle: same
   intensity scale for all six (default) or each panel scaled to its own peak. Panel header: used / failed a
   check / trace, below 1% cutoff / no trace. Folded: what QuantUMS does with these peaks and what is not shown.
10. **Tracing a number in the report text**: the report sentence with its numbers clickable, showing the file,
   rule and step behind each (the original "where did this number come from" idea).
11. **How to read this page**: every term, grouped (the search; from precursors to a protein value; the change;
    checks and flags), each with its source, then a sources list.

Interaction: clicking a sample anywhere selects it everywhere; clicking a precursor or peptide anywhere selects
it everywhere. Folded sections stay open once opened.

## 4. Data inputs

Everything comes from the delivered results folder of one analysis, plus two public services. Names below are
as in a skill 2.14.0 session; the builder must confirm they are in what the collaborator actually receives.

| Input | Used for |
|---|---|
| `search/report.parquet` | per precursor and run: quantities, normalisation factor, quality scores, MS1 area, RT, the six q-value columns, Protein.Ids, Proteotypic, PG.MaxLFQ |
| `search/report.stats.tsv` | per run: precursors and proteins identified, normalisation instability, MS1/MS2 mass error |
| `search/xic/t*_xic/<run>.xic.parquet` | XICs (long format: `pr`, `feature`, `info`, `rt`, `value`; 12 fragments + `ms1`, about ±10 s). Written by `--xic 10 --mobilograms`, which every skill search sets |
| `search/empirical.parquet` | the library searched in the final pass: is any peptide unique to a protein? (section 5.6) |
| `tables/DE-LIMP_session.rds` | `raw_data` (precursor log2 values, exactly log2 of Precursor.Normalised), `y_protein` (values, standard errors, n.observations), `dpc_fit`, `fit` |
| `tables/DE_dpc_<contrast>.csv`, `tables/de_provenance.json` | the delivered result, cutoffs, filters, input checksum |
| the analysis report text | the sentences in block 10 |
| UniProt REST (`rest.uniprot.org/uniprotkb/<acc>.json`) | protein card |
| FRAN public API (`fran.stan-proteomics.org/api/...`), no sign-in | FRAN facts and links (`#/protein/<acc>`, `#/peptide/<seq>`, `#/gene/<gene>`) |

## 5. How each number is made (rules the builder must follow)

### 5.1 Found, not found, inferred, and match-between-runs
- **Found**: DIA-NN reported a value for the precursor in that sample that passed its confidence checks.
- **Not found**: DIA-NN looked and nothing passed. In the skill's search, match-between-runs is the
  empirical-library round trip (`diann_parallel.py` steps 3–4): a library built from every run's first pass is
  searched in every run, so DIA-NN looked for every library precursor in every sample. Never say "not measured"
  or "no measurement".
- **Inferred**: only ever a protein value in a sample where none of its precursors was found
  (n.observations = 0). Peptides are found or not found; they are never "inferred" (decided 10-23).
- DIA-NN writes XICs for more precursors than it accepts at 1% (Feb 2 of the mockup data: 35,661 traced, 27,056
  accepted; its README says "identified precursors" without a threshold). A precursor with a trace but no report
  row is labelled "trace, below 1% cutoff".

### 5.2 Values in the path
- Precursor value used by limpa = log2(Precursor.Normalised) = log2(Precursor.Quantity × Normalisation.Factor)
  (DIA-NN README; checked equal in the data). The factor differs from precursor to precursor within a run.
- Protein values, standard errors and n.observations: read from `y_protein`, never recomputed for display.

### 5.3 Shares: how much of a sample's value rests on each precursor
limpa does not report per-precursor contributions. The screen measures them with limpa's own code: for each
found precursor value, add 0.01 to it, re-run limpa's fit for that protein, and record how much each sample's
protein value moves (share = Δprotein / 0.01).
- The re-run must reproduce the delivered protein values (mockup: within 2×10⁻⁹). Replicate `dpcQuant` exactly:
  hyperparameters from `dpcQuantHyperparam` on the **full** precursor matrix sorted by protein, `prior.sd` and
  `prior.logFC` floored at 1, the protein's own `sigma` (a named per-protein vector), the same start values
  (`imputeByExpTilt`), then `peptides2ProteinBFGS`. See `source/scripts/contrib.R` and `detprob.R`.
- Per sample: shares of found precursors in that sample (own-sample derivatives); "same precursors, other
  samples" = sum of the other derivatives; "no measurement here: detection curve + limpa's prior" = 1 − total.
- Peptide share = sum over its charge states. The bars show peptides; the Sankey shows precursors.
- Label everywhere: "our calculation with limpa", not a limpa output. Finding to state: limpa weights every
  found value equally ("all precursors, in fact all observations, are weighted equally by default", limpa paper
  p.17; weights are 1 in `peptides2Proteins`), so shares depend on how many precursors were found, not on
  intensity.

### 5.4 How limma uses limpa's information
From limpa's code (`dpcDE` → `voomaLmFitWithImputation`, limpa 1.4.0):
- vooma fits a trend of each protein's spread against its average level and the mean log standard error; each
  value's weight is 1 / (predicted SD)⁴, using its own standard error. Inferred values (wide error bars) count
  less.
- If a group's fitted value rests entirely on inferred values, those values are left out of the protein's
  residual variance and degrees of freedom.
- Display the actual weights from re-running `dpcDE` on the delivered `y_protein` (it reproduces the delivered
  logFC to ~10⁻¹⁴), as each sample's share of its group's weighted mean.

### 5.5 The check routes: MaxLFQ + limma, and the Moschem et al. workflow
- Run on every analysis as an internal comparison (section 10, decision 3). Use the skill's own `run_de.R --method maxlfq` (2.14.0 `build_maxlfq.R`: same confidence filters and
  contaminant removal; log2 PG.MaxLFQ; quantile normalisation; keep proteins found in ≥ half the samples; lmFit,
  contrasts, eBayes, no weights; not-found samples left out).
- In a session that records an experiment type, `run_de.R` refuses a MaxLFQ DE without a normalisation check
  fitted for the same method (verified 10-26). The builder must produce that check for the maxlfq design, or run
  the second route as a non-deliverable comparison (section 10, question 1).
- Show both routes side by side: where the protein value comes from, not-found samples, normalisation, weights,
  samples used, change, adjusted p; plus each sample's value under both routes ("different scales; compare
  changes, not levels").
- Label it as the skill's route ("MaxLFQ + limma (our check)"). It is the route the flags compare with.
- **Third route: the Moschem et al. workflow** (added 2026-10-26, Brett: "do all 3"). Moschem JDC, de Barros BCSC,
  Serrano SMT, Chaves AFA, J Proteome Res 2025;24(8):3860–3873, doi:10.1021/acs.jproteome.5c00009, is about
  QuantUMS score filtering; the authors publish a limma workflow alongside it (github.com/41ison/limma-for-proteomics).
  Run it as published: filter `Lib.PG.Q.Value`, `Lib.Q.Value` and `PG.Q.Value` ≤ 0.01; drop contaminant rows;
  `diann::diann_matrix(id.header = "Protein.Ids", quantity.header = "Genes.MaxLFQ.Unique", proteotypic.only = TRUE,
  pg.q = 0.01)` (its default `q = 0.01` also filters precursor `Q.Value`, and it turns zero quantities into NA);
  `lmFit(log2(matrix), ~0 + group, method = "ls")`, contrasts, `eBayes`, BH; called changed at adjusted p ≤ 0.05
  and |logFC| > 0.58. The only change: their contaminant filter matches "cRAP", and the Core FASTA tags
  contaminants `Cont_`, so match `Cont_` in `Protein.Ids`. The DIA-NN R package isn't installed here; the mockup
  copies `diann_matrix`, `pivot` and `pivot_aggregate` verbatim from github.com/vdemichev/diann-rpackage
  (R/diann-R.R). Rows are keyed by `Protein.Ids`; map to `Protein.Group` through the report to line them up with
  limpa (6 rows differ, 4 map to a group twice in the mockup data). Mockup script: `scripts/moschem.R`.
- Checked against the paper's own methods (PRIDE PXD063416 data-processing protocol; the full text is paywalled): DIA-NN 2.0.2, QuantUMS high accuracy; the same three q filters; "unless otherwise mentioned, the QuantUMS scores were set to zero", i.e. no score filter by default (score cutoffs are the paper's tested variations); `diann_matrix` with `id.header = Protein.Names` (the GitHub workflow uses `Protein.Ids`); limma ls + BH. Keying by Protein.Names gives identical results here (3,836 rows vs 3,840; same numbers for BASP1, SURF6, POTEE and the whole run). The authors' QC4DIANN repo adds PG.MaxLFQ.Quality ≥ 0.75 (BASP1 adjusted p 0.093, still 0 calls); not used, since it isn't the paper's default. Whether the paper's conclusions recommend a cutoff is unchecked until someone reads the full text.
- The Moschem route's disagreements are ⚠ (screen only), never ✖: the approved ✖ rule stays tied to the skill's
  check route (my call; section 6's flag split was approved against that route).
- Mockup results (HeLa QC, Feb vs Mar): BASP1 +0.78, adjusted p 0.111 (not called); SURF6 +1.20, 0.124, from 4 of
  6 samples; POTEE and H2BC20P have no unique peptides (Genes.MaxLFQ.Unique is 0 in every run), so it can't test
  2 of limpa's 3 calls. Across the run: 3,840 rows, 3,571 testable, df.prior 2.1, 0 calls (0 even without the
  1.5-fold cut); of 3,516 proteins it and limpa both test, 27 point opposite ways (≥ 0.5 log2) and 160 fall
  outside limpa's 95% interval; 318 limpa-tested proteins can't be tested, 209 of them with a whole group inferred
  in limpa.
- PG.MaxLFQ in DIA-NN 2.x is QuantUMS's MaxLFQ-like quantity, which weights precursors by estimated noise
  (README; QuantUMS paper Supplement p.11). Say so: unlike limpa, it weights precursors.
- Show limpa's 95% interval (logFC ± qt(0.975, df.total) × logFC/t, from the delivered fit) and whether the
  MaxLFQ change falls inside it.
- Show each route's t statistic and the prior degrees of freedom limma borrowed across proteins (eBayes
  `df.prior`), and explain the certainty difference from them. In the mockup data the MaxLFQ route's t statistics
  are *larger* for BASP1 (9.99 vs 7.70) and POTEE (−11.61 vs −7.68), yet its p-values are about 40× larger,
  because it borrowed only 2.3 prior degrees of freedom against limpa's 11.7 (its protein values' noise varies
  more from protein to protein). Do not attribute the difference to inferred samples or weights unless shown:
  BASP1 and POTEE have no inferred samples and even weights.

### 5.6 Shared peptides
- A peptide is shared when DIA-NN's Protein.Ids for it lists more than this protein (the searched database).
  Name the other proteins (UniProt gene names).
- If every peptide a protein uses is shared, check the search library (`empirical.parquet`) for any target
  precursor unique to one of those proteins. If there is none, say "this data can't tell them apart; '<gene>' is
  the name DIA-NN gave the group" (mockup: POTEE; also the POTEF row in that run is mostly actin peptides).
- FRAN's sequence lookup (`/api/peptide/<seq>/proteins`) is context only: it counts isoforms and unreviewed
  entries across all of UniProt.

### 5.7 Detection chance (per precursor, per sample)
- Found value: plogis(β₀ + β₁·value), from `dpc_fit`. Descriptive only: limpa's fit uses the curve only for gaps.
- Gap: the probability of being found that limpa's own fit used: 1 − G, where G = Σ plogis(β₀ + β₁(μ + σ·zₖ),
  lower = FALSE)·wₖ over 16-point normal Gauss quadrature, μ = the fitted precursor level in that sample, σ = the
  protein's sigma. Verified 10-26: an objective rebuilt from these pieces equals limpa's value at its optimum,
  with zero gradient there.
- Show only the probability for a gap, never a value. A gap with chance ≥ 0.8 is "surprising" and pulls the
  value down; ≤ 0.4 "unsurprising".
- limpa's separate per-precursor guesses (`dpcQuantByRow`) are a different calculation that the protein values
  never used. Allowed only behind an off-by-default switch, labelled "not used for the protein value", and never
  as Sankey widths.

### 5.8 Sample health
From report.stats.tsv: precursors and proteins identified, normalisation instability (warn at ≥ 0.3, the skill's
existing line; the DIA-NN README does not define this column, and the screen says so), median MS1 mass error
before → after DIA-NN's correction.

### 5.9 XICs
Read the run's XIC parquet for the protein's precursors. Show raw intensities; never smooth, never fill. A panel
for a sample with no trace says so. The mockup's shared-scale default is what makes weak samples look weak.

## 6. Flags, and where they go

Decided by Brett, 10-26 (the "revised split"). **Provisional** until validated (see below):

| Flag | Rule | Where |
|---|---|---|
| ✖ opposite directions | limpa and MaxLFQ + limma changes have opposite signs and either is ≥ 0.5 log2 | report results table **and** this screen |
| ✖ no second opinion | limpa calls it changed (adj. p < 0.05) and MaxLFQ + limma can't test it (typically on/off proteins) | report results table **and** this screen |
| ⚠ same change, different certainty | calls differ and the MaxLFQ change is inside limpa's 95% interval | this screen only |
| ⚠ calls differ, and so do the changes | calls differ and the MaxLFQ change is outside limpa's interval | this screen only |
| ⚠ outside interval | MaxLFQ change outside limpa's 95% interval | this screen only |
| ⚠ inferred samples | any sample with n.observations = 0 | this screen only |
| ⚠ shared peptides | any peptide used is shared (stronger wording if all are) | this screen only |
| ⚠ unstable normalisation | any run with normalisation instability ≥ 0.3 | this screen only |
| ✓ / – | passed / cannot be tested | this screen only |

Every flag says the answer depends on the method, not which method is right. In the mockup data: limpa calls 3
proteins, MaxLFQ + limma none (all three are "same change, different certainty"; the difference comes from the
borrowed degrees of freedom, section 5.5); 21 of 3,581 point in opposite directions; 250 fall outside limpa's
interval; 253 limpa-tested proteins can't be tested by MaxLFQ + limma, 209 of them on/off. The earlier rule (any
call difference goes in the report table) would have flagged every hit.

**Why the split is provisional.** It keeps "calls differ" out of the client's table on the assumption that the
calls limpa makes and MaxLFQ + limma doesn't are mostly real. The mockup data cannot test that: it is the same HeLa
QC standard on two days, so there is no known answer, and two of limpa's three calls there (POTEE, H2BC20P) rest
only on shared peptides. Before the split is final, run both routes on data with known answers (for example two
species mixed at known ratios, the design Moschem et al. used) and check which route's calls, and which call
differences, match the truth. If limpa-only calls turn out to be unreliable, "calls differ" moves into the report
table.

## 7. Wording rules
- Found / not found for precursors and peptides; "inferred" only for protein values.
- Never "measured" for a found value under match-between-runs; never "no measurement" for a gap.
- Every derived number says who derived it ("our calculation with limpa", "re-run here").
- Every block ends with a provenance line naming files and columns, code, or papers.
- Defaults that flow into text are tagged as defaults (repo rule 2).
- Plain language a biologist can follow; terms defined at the bottom of the page.

## 8. Guardrails (decided 10-23, kept)
- Blocks read only the delivered results files, plus allow-listed public sources. On 2026-10-28 Brett widened the
  list from UniProt and FRAN to: UniProt, NCBI (Gene and Protein, E-utilities), HGNC, Ensembl, and the Human
  Protein Atlas, including its Blood Atlas serum/plasma concentrations. He also asked for protein info when a protein
  is not in UniProt (NCBI Protein, and UniProt's UniParc for retired accessions) and for neXtProt. neXtProt is
  retired (end of life, last release September 2023, verified 2026-10-28), so it is not used.
  The page itself makes no network requests. The builder fetches when the page is made. Every protein or gene fact is
  shown verbatim with its source, release, fetch date and link; nothing about a protein comes from the AI's memory
  (Brett: "make sure it looks up the protein info and does not rely on memory").
- Each block computes its own numbers and labels; the agent picks blocks but cannot change what they say.
- The page checks the results files against their recorded checksums before showing anything
  (`de_provenance.json` → `input_sha256` matched the report in the mockup data) and says so.
- Limit, stated honestly: nothing stops a third-party AI from drawing its own page from scratch.

## 9. Verified facts and caveats the builder can rely on
- limpa 1.4.0 code: found values enter as squared residuals with weight 1; gaps as −2·log P(missed) via the
  detection curve; priors on the protein's mean level (`prior.mean`, `prior.sd`) and on sample deviations
  (`prior.logFC`). Papers: Li, Cobbold & Smyth, bioRxiv 2025, doi:10.1101/2025.04.28.651125; Li & Smyth,
  Bioinformatics 2023;39(5):btad200.
- limpa's paper benchmarks only DIA-NN 1.6–1.8.1, before QuantUMS (from the literature check of 2026-10-26; not
  re-read since, bioRxiv was rate-limiting). This applies to every limpa result the Core delivers from DIA-NN
  searches, not just this screen: the skill pins DIA-NN 2.7.0 and never sets `--direct-quant`, so limpa always
  receives QuantUMS (precision-mode) values.
- DIA-NN may report a zero for any quantity, meaning "below the limit of reliable quantification"; its README
  recommends replacing zeros with NA before log-transforming. limpa's reader does exactly that
  (`y[y < 1e-08] <- NA`), so zeros become gaps for the detection curve. Rare in the mockup data: 7 of 181,261
  report rows.
- QuantUMS (Grossmann, Kistner … Demichev, Nat Biotechnol 2026, doi:10.1038/s41587-026-03131-2) weights a
  precursor's fragment and MS1 signals by inverse estimated variance; the weights are not reported. Its
  standalone code matches DIA-NN 1.8.2 beta; DIA-NN 2.x's version is described as improved, so details may differ.
  Quantity.Quality = 1 / (1 + estimated SD) (README).
- FRAN's fragment traces are cross-run averages (exhaustive check 2026-08-10; not re-checked since, because
  agents do not sign in to FRAN's database). Per-sample peaks must come from the project's own XIC files.
- pyarrow (this Mac, conda) cannot read DIA-NN 2.7.0 parquet files ("Repetition level histogram size
  mismatch"); R nanoparquet can.

## 9a. What delivery does not include today (found in review, verified against skill 2.14.0)
- `core_submission.py` delivers only `report.parquet`, the matrices, `report.stats.tsv`, `report.log.txt` and
  `search_provenance.json` from `search/` (`SEARCH_FILES`), and never XIC folders (`EXCLUDED_DIRS`;
  `references/core-submissions.md`: "never: XIC folders"). So a collaborator's folder has no XICs (step 4) and no
  search library (the uniqueness check in 5.6). See question 1.
- `de_provenance.json` hashes only the search report and the conditions file, not `DE-LIMP_session.rds`.

## 9b. Decisions made in review (the designer's calls; check these)
- **The paragraph is assembled, not fixed text.** Repo rule 1 forbids hard-coding what was done. Each sentence of
  "How a protein's number is made" comes from the analysis's self-description (`de_provenance.json` pipeline id,
  rollup method, filters; `search_provenance.json` engine, version, match-between-runs), with the mockup's text as
  the DIA-NN 2.x + limpa case.
- **Shares may be computed analytically.** Nudging every value is fine for this size but not for 200-run studies;
  the builder may use the fit's derivatives (implicit-function theorem on limpa's objective) instead, provided a
  nudge check on a sample of proteins agrees.
- **The page is a local HTML file** opened in a browser, not a hosted page (unpublished data). UniProt and FRAN
  facts are fetched when the page is made, stored with their fetch date, and shown with live links.
- **Blocks are separate functions with one shared loader** (checksum check, provenance footer, network
  allow-list) so later screens (linked volcano and heatmap, pathway evidence, listen-along, decided 10-23) reuse
  them. In version 1 the agent chooses the protein; the screen shows all blocks.
- **Tracing report-text numbers (block 10)** reads a machine-readable sidecar the report generator writes (each
  number with its file, rule and step), never the prose. If the protein is not mentioned, the block shows the
  report-level numbers only.
- **Checksums**: the pipeline records sha256 for every file the screen reads. A mismatch disables the affected
  blocks with a visible message; the rest still shows.
- **Older analyses**: folders from skill versions before the one that ships this screen get a short message
  saying the screen needs a newer analysis, not a partial page.

## 10. Decided after review (Brett, 2026-10-26)
1. **Every analysis ships a precomputed "screen data" folder**, built on HIVE when the analysis runs (on a compute
   node, inside the analysis job): XIC excerpts for the precursors in the final report, the library-uniqueness
   result for each protein group, shares, detection chances, limma weights, both check routes' results, and
   UniProt and FRAN snapshots with fetch dates. Any AI agent (Claude, ChatGPT, Gemini) can then build the screen
   from the delivered folder without R. `core_submission.py` must deliver this folder; it still never delivers the
   raw XIC folders or the library. Delivery size grows with run count (XICs were ~12 MB per run in the mockup
   data before excerpting); the builder should report the excerpt size per run.
2. **Version 1 scope**: DIA-NN 2.x searches, limpa (DPC-Quant) as the main method, group comparisons with any
   number of groups, one comparison shown at a time. Other designs (covariates, blocking, MaxLFQ as the main
   method, DDA/Sage, FragPipe, Radiant, SILAC, phospho, very large studies) get a reduced screen (protein card,
   precursor table, XICs) with a clear note that the limpa-specific blocks don't apply.
3. **The second route runs on every analysis as an internal comparison**, not a deliverable DE, so it adds no
   staff normalisation decision. The client's report and results tables do not change until the flags are
   validated on known-answer data (section 6). Until then the ✖ flags appear on this screen only.
4. **Three routes on the screen** (Brett, later on 2026-10-26: "do all 3"): limpa + limma (delivered), MaxLFQ +
   limma (the skill's check, which drives the flags), and the Moschem et al. workflow as published (section 5.5).
   The Moschem route also runs inside the analysis job and its results go in the screen-data folder.

## 10a. Questions still open
1. **Which known-answer dataset** validates the flag split and limpa on QuantUMS values (section 6, section 9):
   a Core two-species known-ratio run on a current instrument, or a public benchmark? **First test under way
   (Brett: go, 2026-10-26): PXD063416**, the Moschem et al. data (CC0; Exploris 480; HeLa 200 ng/µL in both mixes,
   yeast 100 vs 300 ng/µL, so truth is human log2 0, yeast +1.585; 4 isolation windows × 3 vs 3). Searched per window
   with the skill (DIA-NN 2.7.0, human + yeast one-per-gene + universal contaminants), all three routes scored against
   the truth. In the authors' own DIA-NN 2.0.2 report, yeast is only 6–7% of precursors in mix A and 11% in mix B
   (5% / 10% of signal), so normalisation barely moves the answer (quantile ≤ 0.03 log2; raw values put human at
   −0.11, DIA-NN-normalised at ~0.00), and yeast is compressed to +1.30–1.47 in every version. It cannot test on/off
   proteins or small changes, and it is not our instrument; a Core-run test is still wanted. Data:
   `/quobyte/proteomics-grp/claude/benchmarks/PXD063416/`.
2. **How the report table will show the ✖ flags**, once validated.

## 10b. What PXD063416 showed (2026-10-27/28)
Each window's session holds `output/benchmark/flags.csv`. It scores the section 6 flags against the skill's MaxLFQ
route, using the known answer (human 0, yeast +1.585).
- **✖ no second opinion holds up.** limpa was wrong on 6 of 29, 7 of 23, 8 of 31 and 5 of 24 such calls
  (2/4/6/10 m/z). For unflagged proteins it was wrong on 5 of 1,300, 16 of 1,365, 22 of 1,090 and 13 of 913. The
  flag belongs in the report table, as decided.
- **The other flags can't be judged yet, because the referee route was broken on this data.** The skill's MaxLFQ
  route adds quantile normalisation, which assumes most proteins don't change. Here about a fifth of the proteins
  rise 3-fold, so the human proteins are pushed down. The route then makes 2,473, 4,811 and 4,248 wrong calls at
  4, 6 and 10 m/z (window report, "Also found: the skill's own MaxLFQ route"), and its disagreements with limpa are
  mostly its own errors.
- **Decision (designer's call):** the referee route uses the normalisation that the analysis's own normalisation
  check chose, not a fixed quantile step. When the check says so, that means DIA-NN's normalised quantities with
  no extra quantile step. The ⚠/✖ split is then re-scored on all four windows before it is final; section 6 stays
  provisional until then.
- The Moschem route reads gene-level quantities. In a two-species sample these pool a yeast protein with its human
  namesake (235 such yeast proteins at 4 m/z). Its disagreements stay ⚠, screen only.

## 12. Build order
Each part ships working software on its own, on branch `feat/skill-protein-screen`.
1. **Screen data, core, and the page** (`docs/superpowers/plans/2026-10-28-protein-screen-1-core.md`).
   `make_screen_data.R` writes `output/screen_data/` for a finished limpa analysis: values, precursors, report rows
   and checks, shares, detection chances, limma weights, run health and checksums. `make_protein_screen.py` builds
   the local page for one protein from that folder, with the UniProt card. Blocks whose data come in a later part
   say so on the page. Both run on an existing session, so the HeLa mockup data and PXD063416 are the first tests.
2. **XICs and library uniqueness** in the screen data, and delivery (`core_submission.py` ships `screen_data/`).
3. **Check routes:** fix the referee route's normalisation (10b), add the Moschem route and the flags, and re-score
   them on PXD063416.
4. **FRAN facts, the short answer, the picture, and report-text tracing** (blocks 3, 4 and 10).
5. **The reduced screen** for other designs and methods (section 10, decision 2), and the message for older analyses.

## 11. Left to the builder
- How the agent gets the data into the page: precomputed per-protein JSON written by the pipeline, or the agent
  reading parquet on demand. Either way, the computations in section 5 must match exactly.
- Visual style, as long as the colours stay validated for colour-blind readers (mockup palette: Feb #2e86de,
  Mar #c2185b, inferred #d4880f light / #c0820f dark, selection #13a089; checked with the dataviz validator).
- Which proteins to precompute, caching, and performance.
- Exact sentence templates of the short answer, as long as each sentence is computed from the data and links to
  its evidence.
