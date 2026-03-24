# Claude Analysis Drift Test — Methodology

## Purpose

The drift test monitors whether Claude's interpretation of a fixed proteomics dataset changes over time as Anthropic updates the model. This is critical for reproducibility: if an AI assistant is part of a data analysis pipeline, we need to know whether its conclusions are stable.

## How It Works

### 1. Fixed Input

Every week, we send Claude the **exact same** prompt and data:
- The DE-LIMP example dataset (HeLa proteomics, 6 comparisons)
- A structured prompt requesting analysis of DE results, QC metrics, and biological interpretation
- The prompt includes the full DE results CSV, expression matrix, and QC metrics

### 2. Structured Extraction

From Claude's free-text response, we extract quantitative metrics:

| Metric | What It Measures | How |
|--------|-----------------|-----|
| **Key Proteins** (9 total) | Does Claude identify the known biology? | Check for NAMPT, ALDH3A1, LGALS1, C3, ITGA2, PRSS23, TIMP1, GSTM3, CGA |
| **Stable Markers** (3 total) | Does Claude correctly identify non-changing proteins? | Check for MED14, TUBA1B, TBC1D17 |
| **Sections** (7 total) | Does the response cover all expected topics? | Check for: overview, QC assessment, key findings, cross-comparison, biomarker, biological interpretation, how this analysis works |
| **Genes Mentioned** | How many genes does Claude discuss? | Regex for uppercase 2-10 char words, filtered against noise list |
| **Gene Overlap (vs prev)** | Run-to-run consistency | Jaccard-like: \|intersection\| / \|previous genes\| |
| **Gene Overlap (vs baseline)** | Long-term drift from first run | Same metric but compared to the very first golden baseline |
| **Word Count** | Response verbosity | Total words in response |
| **Up/Down Mentions** | Directional language balance | Count of "upregulated/increased" vs "downregulated/decreased" |
| **P-value References** | Statistical rigor | Count of "p-value", "FDR", "q-value" mentions |
| **Fold Change References** | Quantitative specificity | Count of "fold change", "logFC" mentions |
| **Decimal Numbers** | Data citation density | Count of specific numbers (e.g., "0.05", "1.32") |
| **Hedging Language** | Scientific caution | Count of "suggest", "may", "could", "potential", "appears", "likely", "possible" |
| **Confident Language** | Definitive statements | Count of "clearly", "strongly", "significantly", "definitively", "robust" |
| **Model ID** | Which Claude model was used | Captured from API response (e.g., `claude-sonnet-4-20250514`) |
| **Input/Output Tokens** | API usage and response size | Token counts from API response metadata |

### Hedging vs Confident Language — Why We Track This

In scientific writing, **hedging** is the use of cautious, tentative language to qualify claims:
- *"These results **suggest** that NAMPT **may** play a role..."*
- *"ALDH3A1 **could** be a **potential** biomarker..."*
- *"The data **appears** to indicate..."*

**Confident language** makes stronger, more definitive assertions:
- *"NAMPT is **clearly** upregulated..."*
- *"The results **strongly** support..."*
- *"This protein is **significantly** enriched..."*

**Why it matters for drift detection:**
- A good scientific analysis should **hedge more than it asserts** — proteomics data with 3 replicates per group warrants caution
- If the hedge/confident ratio shifts dramatically between model versions, the model's "personality" may have changed
- Too little hedging = overconfident claims that could mislead researchers
- Too much hedging = vague analysis that isn't actionable
- Our March 24 baseline: 18 hedges vs 2 confident statements (9:1 ratio) — appropriately cautious

**Tracked keywords:**
| Type | Keywords |
|------|----------|
| Hedging | suggest, may, could, potential, appears, likely, possible |
| Confident | clearly, strongly, significantly, definitively, robust |

### 3. Ground Truth

The 9 **key proteins** are biologically validated markers from the example dataset:
- These proteins are consistently DE across multiple comparisons
- A competent analysis should identify all of them
- Missing any is a signal that the model may be under-exploring the data

The 3 **stable markers** are proteins with low CV and no significant DE:
- A good analysis should mention these as internal controls or housekeeping candidates
- Identifying stability is harder than identifying change

The 7 **sections** ensure Claude covers the full analytical scope:
1. **Overview** — dataset summary, sample counts, groups
2. **QC Assessment** — sample quality, outliers, completeness
3. **Key Findings** — top DE proteins, fold changes, significance
4. **Cross-Comparison** — proteins significant in multiple contrasts
5. **Biomarker** — potential biomarker candidates
6. **Biological Interpretation** — pathways, functions, mechanisms
7. **How This Analysis Works** — educational explanation of methods

### 4. Trend Assessment

With 4+ weekly baselines, the report computes:

| Trend | Method | Interpretation |
|-------|--------|----------------|
| **Gene overlap trend** | Linear regression on overlap % over runs | STABLE (±2%/run), IMPROVING (positive slope), DEGRADING (negative slope) |
| **Baseline drift** | Overlap vs first-ever run over time | CONVERGING (getting more similar to original), DRIFTING (diverging) |
| **Word count trend** | Linear regression on word count | STABLE (±20 words/run), GROWING, SHRINKING |
| **Gene count trend** | Linear regression on genes mentioned | STABLE (±1 gene/run), EXPANDING, CONTRACTING |

### 5. Gene Stability Analysis

- **Core genes**: Mentioned in every single run. These are the most robust findings — Claude always highlights them regardless of model version.
- **Frequent genes** (≥75% of runs): Reliably mentioned but occasionally dropped.
- **One-off genes**: Mentioned in exactly one run. These are stochastic — model temperature/sampling variation.

### 6. Overall Verdict (4+ runs)

| Verdict | Criteria |
|---------|----------|
| **HEALTHY** | Mean overlap ≥65%, all key proteins found, no alerts |
| **ACCEPTABLE** | No critical issues but overlap below 65% |
| **ATTENTION NEEDED** | Overlap below 50%, missing key proteins, or high word count variance |

### 7. Alerts

The report flags specific concerns:
- Missing key proteins in any run
- Gene overlap below 50% (vs previous) or 40% (vs baseline)
- Word count below 500 or above 2.5x baseline
- Gene count swing >50% between consecutive runs

## Infrastructure

- **Schedule**: GitHub Actions workflow runs every Monday at 3 AM Pacific
- **Storage**: Golden baselines saved as `.rds` files in `tests/golden/`
- **API**: Uses Claude API via `httr2` with the project's `ANTHROPIC_API_KEY` secret
- **Commit**: Workflow auto-commits new baselines back to the repo
- **Artifacts**: Each run's response also saved as a downloadable GitHub Actions artifact (90-day retention)

## File Locations

| File | Purpose |
|------|---------|
| `tests/testthat/test-claude_export.R` | Test suite — builds prompt, calls API, extracts findings, compares to golden |
| `tests/drift_report.R` | Standalone report generator — reads all golden baselines, computes trends |
| `tests/golden/claude_findings_*.rds` | Extracted metrics per run (structured R list) |
| `tests/golden/claude_response_*.md` | Full Claude response text per run |
| `tests/golden/drift_report.csv` | CSV summary of all runs |
| `.github/workflows/claude-drift-test.yml` | Weekly GitHub Actions workflow |

## Interpreting Results

**75% gene overlap is normal.** Claude's analysis involves some stochastic variation — it won't mention the exact same 28 genes every time. What matters is:

1. **Key proteins are always found** (9/9) — the core biology is stable
2. **Core genes grow over time** — more genes become consistently mentioned
3. **One-off genes stay low** — most variation is in which secondary genes get highlighted
4. **No sudden drops** — a week where overlap falls to 30% would signal a model change

**What would be concerning:**
- Key protein recall dropping below 8/9
- Gene overlap trending below 50% over 4+ weeks
- Word count doubling or halving (prompt interpretation changed)
- Sections dropping (model no longer covers expected topics)
