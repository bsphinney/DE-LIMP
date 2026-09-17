# AI prompts used by DE-LIMP

**Why this file exists.** The prompts are assembled inline in R with `paste0()` from dozens of
separate string fragments, so most sentences span concatenation boundaries and **cannot be found by
searching the repo.** Searching for `"You are analyzing a proteomics cross-tool comparison"` happens
to match, because that is one whole fragment; searching for `"Structure your response as follows"`
finds nothing, because that text is split across two arguments. This file is the readable copy.

**The code is the source of truth.** This document is a transcription. If you change a prompt, update
this file in the same commit — there is a pointer comment above each builder saying so.

**To see the real prompt with your own data in it**, don't read either the source or this file: the
Run Comparator has a **"View Prompt"** button that renders the exact assembled string in a modal with
copy-to-clipboard (`R/server_comparator.R`, the `comparator_view_prompt_btn` observer). That is always
current by construction.

Notation below: `<name>` marks a value interpolated at runtime. Indented `[if ...]` marks a block that
is only included under that condition.

---

## 1. Run Comparator — in-app Gemini narrative

**Built by:** `build_gemini_comparator_prompt()` — `R/server_comparator.R`
**Sent by:** the `comparator_ai_btn` observer, via `ask_gemini_text_chat()`
**Shown to the user by:** the "View Prompt" modal

### 1.1 Full template

```text
You are analyzing a proteomics cross-tool comparison. Be objective — neither tool is inherently superior. Evaluate based on the data.

INSTRUMENT: <model>, <lc_system>, <lc_method>, <n> SPD, <n> min gradient      [if instrument metadata present; each part omitted if empty]

TOOL COMPARISON CONTEXT:
<tool_context — one of the four branches in §1.2>

COMPARISON OVERVIEW:
- Run A: <run_a_label> (<n_samples> samples)
- Run B: <run_b_label>
- Contrast: <contrast>
- Proteins shared: <n_shared>
- Concordant DE: <n_concordant>
- Discordant: <n_discordant> (<n_a_only_de> DE in A only, <n_b_only_de> DE in B only)
- Global intensity offset (median): <global_offset> log2 units
- Per-sample correlation range: <min_cor>-<max_cor>

DIA-NN LIBRARY SIZES:                                        [if a DIA-NN log was parsed for Run A]
Run A: <n> precursors (pg-level <n>), --proteoforms
Run B: <n> precursors (pg-level <n>), --proteoforms          [or "Run B: log not provided"]

SETTINGS DIFFERENCES:
<one "- <Parameter>: <Run_A> vs <Run_B>" line per row whose match is differs / structural_difference / severe;
 or "All settings identical">

TOP DISCORDANT PROTEINS (protein_id | logFC_A | logFC_B | adjP_A | adjP_B | category):
<up to 10 rows, pipe-separated; or "No discordant proteins">

MOFA2 FACTOR DECOMPOSITION:                                  [if a MOFA2 object was supplied]
- <Factor_n> is the most run-A-specific factor
  (explains <pct>% variance in Run A, <pct>% in Run B)
- Top proteins driving this factor: <top 5 features by |weight|>
- This suggests Run B's pipeline is suppressing the signal captured by this factor.

Structure your response EXACTLY as follows:

## 1. Factual Observations
List only what the data shows — concordance rates, fold-change distributions, correlation ranges, number of significant proteins per tool, global offset. No interpretation yet. Cite specific numbers from the comparison overview above.

## 2. Sources of Disagreement
For each major source of disagreement, state:
- What the difference is (quantification method, statistical model, normalization, etc.)
- How many proteins it likely affects (systematic vs protein-specific)
- The expected direction of its effect (more/fewer DE calls, shifted fold-changes, etc.)
Ground claims in proteomics literature where possible (e.g., limma empirical Bayes: Smyth 2004; MaxLFQ: Cox et al. 2014; DIA-NN: Demichev et al. 2020; Spectronaut: Bruderer et al. 2015).

## 3. Case for Run A
Argue the strongest case for trusting Run A's results for this specific experiment. What does Run A do well here? What in the data supports this?

## 4. Case for Run B
Argue the strongest case for trusting Run B's results for this specific experiment. What does Run B do well here? What in the data supports this?

## 5. Settings Audit
Review the SETTINGS DIFFERENCES table above. For each difference:
- Could this specific setting explain any of the observed discrepancies?
- Are any settings misconfigured, unusual, or potentially erroneous for this type of experiment?
- Flag anything that looks like a mistake (e.g., wrong enzyme, unexpected mass accuracy, mismatched FASTA, unusual normalization choice for the sample size).
Be specific about which discrepancies each setting could cause.

## 6. Biology of Concordant Proteins
The concordant DE proteins (significant in both pipelines with the same direction) represent the highest-confidence biological signal. Briefly characterize:
- Are they enriched in specific pathways, compartments, or functional categories?
- Do the concordant proteins suggest a coherent biological narrative for this comparison?
- Are there notable proteins in the concordant set that are well-known markers or functionally relevant?
If the protein IDs are UniProt accessions rather than gene symbols, note any you recognize.

## 7. Synthesis
Now weigh the cases from sections 3-4. Where do the pipelines agree? Where is the evidence genuinely ambiguous? If one pipeline is clearly more appropriate for this experiment design (sample size, data type), explain why with reference to the literature. If the evidence does not clearly favor one, say so explicitly — do not force a recommendation.

## 8. Recommended Follow-ups
Two concrete, actionable follow-up analyses — one investigating a potential issue in each pipeline (not both focused on the same tool). Be specific about which proteins to examine, which settings to change, or which validation to run.

GUIDELINES:
- Both pipelines are established and peer-reviewed. Do not assume one is inherently superior.
- Every claim must be supported by a specific number from the data or a literature reference.
- If a setting (e.g., Quant3) is unusual, discuss its potential impact but do not treat it as automatically invalidating results — the actual effect depends on the experiment.
- Label speculation clearly: 'This suggests...' or 'One possible explanation is...' vs stating as fact.
- If discordant proteins show a clear pattern (e.g., all low-abundance, all membrane proteins), note it. If no pattern is evident, say so.
```

### 1.2 The four `TOOL COMPARISON CONTEXT` branches

Selected on `summary_stats$source_b` — what Run B was loaded from.

**`delimp`**

```text
Both runs used DE-LIMP (DIA-NN -> limpa/limma pipeline). The comparison isolates the effect of search or analysis parameter differences on the same raw data.
```

**`spectronaut`** — the long branch, and the only one that computes statistics of its own:

```text
TOOL COMPARISON CONTEXT:
Run A: <run_a_pipeline_label>. <run_a_rollup_text> <run_a_de_engine_text>
Run B: Spectronaut <version> -- directDIA+ search -> MaxLFQ (TopN=<topn>, uses <topn> most intense peptides per protein) -> <de_test>.

QUANTIFICATION APPROACH DIFFERENCE:
Spectronaut: TopN=<topn> peptide selection for protein quantification. The TopN cap applied to ~<pct>% of proteins (average <n> peptides detected). Run A: <run_a_rollup_text> Both approaches have trade-offs: TopN may be more robust to noisy peptides; <all-precursor empirical Bayes uses more data but may include low-quality signals. | DIA-NN's PG.MaxLFQ is also a paired-ratio aggregator, but operates over the full pre-filtered precursor set.>

STATISTICAL MODEL DIFFERENCE:
Run A (limma): moderated t-test with empirical Bayes variance shrinkage — called <n_sig_a> proteins significant. Run B (<de_test>): called <n> significant at q<=0.05 (<pct>% of <n_total> total). Median |logFC| among Run B significant proteins: <n> log2. These methods differ in how they estimate per-protein variance, which can produce different significance calls even from similar fold-changes.

PRE-FILTERING DIFFERENCES:
Run B applies a Log2 Ratio Candidate Filter (<threshold> log2) BEFORE statistical testing — proteins with insufficient fold-change evidence are excluded from DE testing entirely (NaN ratios). Run B also uses <imputation> imputation. Run A: <run_a_missing_text>

NORMALIZATION DIFFERENCE:
Run A: DIA-NN RT-dependent normalization upstream + <run_a_norm_text>. Run B: <normalization>. Different normalization strategies can introduce a global intensity offset.

QUANT3 SETTING:                                              [only if "Use All MS-Level Quantities" was ON]
Spectronaut's 'Use All MS-Level Quantities' was ON. This combines MS1 and MS2 measurements, increasing the observation count per protein in the t-test. This can increase statistical power but also raises questions about observation independence (MS1 and MS2 from the same peptide are correlated). Consider whether this affects the p-value comparison.
```

**`fragpipe_analyst`**

```text
Run A: DE-LIMP — <run_a_pipeline_label>.
Run B: FragPipe + FragPipe-Analyst (MSFragger -> IonQuant MaxLFQ rollup -> limma DE).
Key structural differences:
- Search engine: DIA-NN (DIA-optimized) vs MSFragger (DDA/DIA)
- Protein rollup: <run_a_rollup_short> vs MaxLFQ (pairwise ratios)
- Normalization: <run_a_norm_short> vs IonQuant (optional VSN)
- Missingness: DIA-NN MBR vs IonQuant MBR (different implementations)
- Missing values: <run_a_missing_text> vs Perseus-style imputation (FP-Analyst default)
Both use limma for DE, so p-value calibration should be broadly comparable.
```

**`fragpipe_raw`**

```text
Run A: DE-LIMP (DIA-NN). Run B: raw FragPipe output (no DE stats). Only quantification differences can be assessed.
```

### 1.3 The pipeline-aware Run-A substitutions

Every `<run_a_*>` placeholder above is resolved from `comp_results$run_a$settings$pipeline_id`, **not
hardcoded** — this is the v3.10.0 fix for exports that described DPC-Quant on MaxLFQ-quantified data
(see the "Architectural rules" section of `CLAUDE.md`, rule 1).

| placeholder | `pipeline_id == "dpc"` (default) | `pipeline_id == "maxlfq"` |
|---|---|---|
| `run_a_rollup_short` | `DPC-Quant` | from `settings$rollup_method` |
| `run_a_rollup_text` | DPC-Quant aggregates all detected precursors with empirical Bayes weighting. | DIA-NN PG.MaxLFQ pivot to a Protein.Group x Run matrix, log2 + quantile-normalized cross-sample, then plain limma::lmFit (NA-tolerant per row). |
| `run_a_de_engine_text` | limpa::dpcDE (limma adapted for DIA proteomics) -> contrasts.fit -> eBayes. | limma::lmFit -> contrasts.fit -> eBayes (no DPC-Quant detection model). |
| `run_a_missing_text` | DPC-Quant models missing values probabilistically via a detection probability curve (not imputation), and tests all quantified proteins. | NAs left in place; limma drops them per row at fit time. Proteins entirely missing in one condition produce NA logFC and are surfaced separately as on/off calls (no imputation). |
| `run_a_norm_short` | DPC-CN cyclic loess | quantile (limma::normalizeBetweenArrays) |

---

## 2. Run Comparator — export prompt (`claude_prompt.md`)

**Built by:** `build_claude_comparator_prompt()` — `R/server_comparator.R`
**Written to:** `claude_prompt.md` in the Comparator export ZIP (export step 6)
**Purpose:** pasted into a chat alongside the ZIP's CSVs, so it describes the attachments.

```text
I am sharing a proteomics run comparison from DE-LIMP.

INSTRUMENT: <model>, <lc_system>, <lc_method>, <n> SPD, <n> min gradient      [if instrument metadata present]
COMPARISON TYPE: <tool_label — see §2.1>
EXPERIMENT: <contrast> contrast, <n_samples> samples.

METHODOLOGY NOTE:
<DPC branch or MaxLFQ branch — see §2.2>
<Spectronaut addendum — see §2.3>                            [only if source_b == "spectronaut"]

KEY FINDING: Of <n_shared> shared proteins, <n_discordant> are discordant (<n_a_only_de> DE in Run A only, <n_b_only_de> DE in Run B only).
Most common discordance pattern: <dominant hypothesis_category>.
Global intensity offset: <n> log2 units (<SYSTEMATIC BIAS DETECTED if |offset| > 0.2, else no systematic bias>).

AI PRE-ANALYSIS (<provider>, model <model>):                 [if the in-app AI narrative was generated first; provider + model from ai_source_label()]
<narrative text>

FILES ATTACHED:
- discordant_proteins.csv: All disagreements with per-protein diagnostic flags <(includes n_ratios_B column) if Spectronaut>
- de_results_combined.csv: Full DE stats from both runs side-by-side
- settings_diff.csv: Parameter comparison table
- protein_universe.csv: All proteins with tier classification
- diann_search_params.txt: DIA-NN search parameters for Run A (if available)
- precursor_summary_discordant.csv: Precursor-level data for discordant proteins (if available)
- spectronaut_run_qc.csv: Per-sample QC metrics from Spectronaut (if available)      [Spectronaut only]
- spectronaut_library_info.csv: Spectronaut library composition and version          [Spectronaut only]

QUESTIONS:
1. Based on the discordant protein patterns and tool differences, what is the most likely root cause?
2. Are there proteins where the two runs disagree biologically (not just statistically)?

   [if Spectronaut:]
3. Do low n_ratios_B proteins in discordant_proteins.csv cluster biologically, or are they random?
4. If spectronaut_run_qc.csv is present, are outlier samples driving specific discordant proteins?
5. Among proteins significant only in DE-LIMP: do they have high or low n_ratios_B in Spectronaut?

   [otherwise:]
3. If precursor_summary_discordant.csv is included, do any discordant proteins have unusual precursor characteristics that might explain the disagreement?
4. What additional information would resolve the ambiguity?
```

### 2.1 `COMPARISON TYPE` values

`run_a_short` is `DIA-NN/MaxLFQ/limma` or `DIA-NN/DPC-Quant/limma` depending on `pipeline_id`.

| `source_b` | label |
|---|---|
| `delimp` | two DE-LIMP sessions (same pipeline, different settings) |
| `spectronaut` | DE-LIMP (`<run_a_short>`) vs Spectronaut |
| `fragpipe_analyst` | DE-LIMP (`<run_a_short>`) vs FragPipe+FragPipe-Analyst (MSFragger/MaxLFQ/limma) |
| `fragpipe_raw` | DE-LIMP (`<run_a_short>`) vs FragPipe raw output (no DE stats) |
| anything else | unknown comparison |

### 2.2 `METHODOLOGY NOTE` — the two pipeline branches

**DPC-Quant** (`pipeline_id != "maxlfq"`):

```text
DE-LIMP uses DPC-Quant (Detection Probability-based Combined Quantification) from the limpa R package. DPC-Quant does NOT leave missing values, and it is NOT traditional imputation. Instead, it uses an empirical Bayes model that jointly estimates protein abundance and detection probability across all precursors. Proteins with zero detected precursors in a sample still receive abundance estimates informed by the detection probability model — these are statistically valid but less precise (higher SE). The missing_pct columns in the export reflect the RAW precursor-level missingness BEFORE DPC-Quant, not the final protein matrix (which is always complete). DE testing uses limma moderated t-tests with empirical Bayes variance shrinkage.
```

**MaxLFQ** (`pipeline_id == "maxlfq"`):

```text
DE-LIMP processed Run A through the MaxLFQ + limma pipeline (Moschem et al., J. Proteome Res. 2025; 24:3860). Precursor rows in the DIA-NN report were filtered at 1% FDR plus optional QuantUMS quality cutoffs (Empirical.Quality, PG.MaxLFQ.Quality), then aggregated to a Protein.Group x Run matrix using DIA-NN's PG.MaxLFQ values. log2 + quantile-normalized cross-sample (limma::normalizeBetweenArrays). Plain limma::lmFit + eBayes was used for DE — NAs are left in place; limma drops missing values per row at fit time. Proteins entirely missing in one condition produce NA logFC and are surfaced separately as on/off calls (no imputation, no detection-probability modelling). The missing_pct columns in the export reflect missingness in the FINAL MaxLFQ protein matrix.
```

Note the deliberate asymmetry in the closing sentence: `missing_pct` means *precursor-level
missingness before rollup* under DPC-Quant but *missingness in the final protein matrix* under MaxLFQ.
Getting that backwards would make a reader misjudge data completeness, which is why it is stated in
the prompt rather than left implicit.

### 2.3 Spectronaut addendum

```text
Spectronaut uses MaxLFQ on TopN peptides with a standard (Welch) t-test. The moderated t-test is substantially more conservative for small sample sizes, which is a key structural difference when interpreting discordant DE calls.
```

---

## 3. Design intent worth preserving

The Comparator prompts are deliberately structured against tool bias, and the structure is easy to
erode by well-meaning editing:

- **Separate "Case for Run A" and "Case for Run B" sections, before any synthesis.** Both are argued
  at full strength before section 7 weighs them.
- **Section 7 explicitly permits declining to recommend**: *"If the evidence does not clearly favor
  one, say so explicitly — do not force a recommendation."*
- **Section 8 forbids one-sided follow-ups**: *"one investigating a potential issue in each pipeline
  (not both focused on the same tool)."*
- **An unusual setting is not automatically disqualifying** — the Quant3 guideline says to discuss
  impact without treating it as invalidating.
- **Every claim needs a number or a citation**, and speculation must be labelled as such.

---

## 4. Other prompt builders in the app

Not transcribed here. Same findability caveat applies to all of them — each is assembled with
`paste0()` in R.

| area | file |
|---|---|
| AI Summary, Data Chat, HTML report | `R/server_ai.R` |
| Phospho site-level interpretation | `R/server_phospho.R`, `R/helpers_phospho.R` |
| MOFA2 factor interpretation | `R/server_mofa.R` |
| De novo / DDA | `R/server_dda.R`, `R/server_denovo.R`, `R/helpers_dda.R` |
| Proteogenomics | `R/helpers_proteogenomics.R` |
| Instrument diagnostics | `R/helpers_instrument.R` |
| Session / reproducibility / methods text | `R/server_session.R` |
