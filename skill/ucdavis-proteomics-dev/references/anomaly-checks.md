# Anomaly checks, expert review & analysis principles

The UC Davis Proteomics Core's data-analysis methodology, ported into the skill. Apply
it during **step 8c (audit)** and **step 9 (write the report)**. `audit_results.py` and
`sample_quality.py` automate the deterministic checks; **you apply the judgment branches
below** on the tables/figures you already have. Anything you flag goes into the report's
**Data Quality Notes** section and — if it changes a conclusion — is surfaced to the user
immediately, not buried.

## Core principles

- **Ask before assuming.** Experimental design, instrument, organism, reagents, and the
  biological question are confirmed with the user *before* analysis. Filenames are not
  authoritative ("SM2"/"SM4" could be cross-linkers, samples, or instruments). A 30-second
  clarification prevents hours of rework. (This is golden rule #1 — it applies here too.)
- **Verify tool flags against the binary, not memory.** Order of trust: `tool --help` on
  the exact binary → official docs for the exact version → source → third-party guides.
  Don't assume a flag exists because a sibling tool has it.
- **Flag anomalies as they appear**, in a **5-part format**: *what* (the specific anomaly)
  · *where* (rows/columns/files/samples) · *why it matters* (impact on the result) ·
  *likely cause* · *suggested fix*.
- **Trust observation over documentation.** If a path/convention here conflicts with what
  `ls`/`head`/the file shows, trust the file.
- **Every report gets a "Data Quality Notes" section**, even if it just says "nothing
  anomalous observed."

## Anomaly-check decision tree — run only the branches that apply

**Always (any dataset):**
- Missing values — pattern (random vs. systematic per condition), count per sample,
  unexpected NAs. *(automated: `audit_results.py` missingness)*
- Intensity distributions — one sample far off in median/spread ⇒ failed injection or
  loading error. *(judgment: read the signal-distribution / per-sample figures)*
- Identifier sanity — organism match, deprecated accessions, mixed naming, duplicate IDs.
- Contaminants in the top-20 by intensity (keratin, trypsin autolysis, BSA) ⇒ prep issue.
  *(automated: `audit_results.py` contamination; deeper: `sample_quality.py`)*
- File-format sanity — delimiter, encoding, column-count consistency, header alignment
  across files meant to be compared.

**If multiple files / conditions:**
- Replicate ID overlap < 60%, or between-condition overlap < 40% ⇒ investigate before DE.
- Systematic mass / RT / CCS shifts between files ⇒ different reagent or method, **not**
  biology — flag immediately.
- **Batch confounded with condition** (all treatment day 1, all control day 2) ⇒ DE is
  invalid. *(automated: `audit_results.py` confounding — a FAIL)*
- Unbalanced group sizes (e.g. 5 vs 2) ⇒ note the power implication. *(automated: group_balance)*

**If statistical testing (DE, enrichment):**
- P-value distribution shape — flat = no signal; spike at 0 **and** 1 = model problem.
- Largest fold-changes all in low-abundance proteins ⇒ noise-driven.
- log2FC > 5 in a non-knockout experiment ⇒ probably imputation or single-replicate artifact.
- No FDR correction reported ⇒ flag.
- PCA/MDS shows a sample clustering with the wrong group ⇒ possible **sample swap** — flag
  immediately. *(judgment: read the PCA figure vs. conditions)*
- ORA/GSEA without an explicit background ⇒ inflated significance.

**If XL-MS:** bridge mass consistent across files? (same pair, different mass = different
cross-linker) · decoy-implied FDR > ~10% is too hot even for XL-MS · one lysine in
disproportionately many links = hyperreactive site, not many true contacts · structure-mapped
distance violations > 30% = problem.

**If phospho / PTM:** localization probability < 0.75 ⇒ don't report as confidently
assigned · non-phosphopeptides dominant ⇒ enrichment failed · known active-pathway sites
absent ⇒ sanity-check the experiment.

**If non-model organism:** using human/mouse GO/KEGG without ortholog mapping ⇒ annotation
mismatch (map orthologs first).

## Contamination diagnostics — `sample_quality.py`

Beyond lab contaminants (keratin/trypsin), the **sample-level biological contaminations**
that repeatedly produce wrong-but-plausible DE:
```
python3 scripts/sample_quality.py --matrix de_results/Expression_Matrix.csv \
    --conditions conditions.csv --report search_out/report.parquet --out SAMPLE_QUALITY.md
```
- **HEMOLYSIS** (plasma/serum RBC lysis: hemoglobin, carbonic anhydrase, catalase,
  peroxiredoxin, spectrin) — drove artefactual plasma DE in real cases.
- **SKELETAL_MUSCLE** (muscle debris in a non-muscle biopsy: myosins, troponins, CK-M,
  myoglobin) — in one cow-liver study, uniform muscle contamination of one breed's biopsies
  produced 3,045 "DE" proteins that were contamination, not biology.
- **EPIDERMIS** (skin/hair squames).

It scores each panel per sample (z across samples) and checks whether the score is
**confounded with group**. Two hard-won lessons it encodes:
- **A contamination panel that separates the groups is the danger signal** — DE between
  those groups may be the contamination gradient. **Protein-level marker removal does NOT
  fix a confounded contrast** (dropping muscle markers once *increased* the DE count); the
  fix is at the sample/design level (drop/flag the confounded samples, or model it).
- **The limpa/DPC-Quant "complete matrix" depth trap:** a DPC-Quant expression matrix has
  ~no NAs (the detection model fills every cell), so counting non-empty cells gives a
  **constant**, not per-sample depth. Real depth = detected protein groups per run
  (`--report report.parquet`).

### ⚠ Keratin-matrix samples — keratin is the ANALYTE, not a contaminant

For **nail, hair, wool, skin, feather** (and other epidermal/keratinaceous) samples,
keratins (KRT/KRTAP) are *what you are measuring* — **never flag them as contamination.**
Pass **`--keratin-sample`** to both `sample_quality.py` and `audit_results.py`: the
EPIDERMIS/keratin panel is then reported for QC but excluded from contamination flags,
and the audit stops warning on KRT/KRTAP/keratin (trypsin/casein/BSA are still flagged).
The generic "contaminants in the top-20 by intensity (keratin, trypsin, BSA)" check means
prep contamination **only for non-keratin matrices**; on a keratin sample, keratin at the
top is expected biology. (Blood/hemolysis and muscle panels still apply — a nail sample
can still be blood-contaminated.)

## Conditional expert review — spawn only on triggers

Most analyses don't need a review panel. Spawn one when **any** trigger fires:
- the session produces a **statistical claim going to a collaborator** (DE, enrichment,
  interaction claims);
- any **log2FC > 5** or **p < 1e-10** is reported;
- a **new organism, instrument, or reagent** combination for this project;
- **cross-tool discordance** surfaced (DIA-NN vs. Spectronaut, etc.);
- the user asks for review.

When triggered, spawn **three agents in parallel** (one message, three `Agent` calls),
each given the report + input/output paths + experimental context. Each returns a bulleted
list of *issue · severity (critical/warning/note) · suggested action*:
- **Proteomics reviewer** — search-output interpretation, FDR/decoy handling, quantification
  fit, instrument-specific artifacts, mass/RT/CCS consistency, acquisition-vs-processing match.
- **Biology reviewer** — do the results make sense for this sample/organism? Expected proteins
  present? Contradictions with known biology? **Could contamination explain the finding?**
- **Statistics reviewer** — normalization fit, multiple-testing correction, confounding /
  pseudoreplication, sample size, missing-value handling, test choice, PCA/MDS for swaps.

Consolidate (dedupe, sort by severity) into a **`## Expert Review Notes`** section of the
report; surface every **critical** issue to the user before finalizing. **Skip review** for
format conversion, single-column extraction, exploratory looks with no formal report, FASTA
prep, and annotation correction against a single reference.

(These three lenses also match the release-review agents the project runs; here they are
gated to per-analysis triggers rather than run every time.)

## Collaborator deliverable — 3 reproduction options

The reproducibility bundle (step 10) is mandatory; when the deliverable is collaborator-facing,
its `REPRODUCE.md` should offer **three** explicit paths: **(A)** re-fit the statistics only,
**(B)** re-render every figure and rebuild the report, **(C)** inspect the raw data without
re-running anything — plus background, headline findings, next steps, provenance (cluster
paths/md5 for raw `.d`/`.raw`), and a contact email.
