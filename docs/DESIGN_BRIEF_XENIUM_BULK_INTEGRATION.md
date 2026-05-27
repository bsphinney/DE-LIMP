# Xenium + Bulk LC-MS Proteomics Integration — Design Brief

**Status:** design / scoping
**Audience:** DE-LIMP devs + UC Davis Proteomics Core
**Companion to:** DNA Tech Core Xenium service offering, Proteomics Core LC-MS/DIA-NN workflow

This is a design brief, not an implementation. It captures the integration plan, biological motivation, lit review, and a roadmap so we don't lose context across sessions.

---

## Why this exists

The UC Davis DNA Technologies Core offers PacBio Iso-Seq, Xenium spatial transcriptomics (10X Genomics), Visium HD, and single-cell RNA-seq (10X, Parse). The Genome Center director is enthusiastic about Xenium as a strategic platform. For the Proteomics Core to position itself as a necessary complement (rather than a competing service), we need a clear integration story that shows what bulk LC-MS proteomics adds to Xenium that Xenium fundamentally cannot answer.

This document is the technical roadmap for that integration, with the user-facing positioning bundled in for sales / collaboration discussions.

---

## What bulk LC-MS uniquely provides over Xenium

Xenium gives subcellular-resolution transcript counts for a 300–5,000-gene panel with spatial coordinates. It does NOT see:

1. **Protein abundance** — RNA-protein correlation is ~0.4–0.6; half the proteins in a cell don't track their mRNA because of translation rate, half-life, and ubiquitin-proteasome turnover
2. **Post-translational modifications** — phospho, glyco, acetyl, ubiquitin, methyl. These are the regulatory layer. Drug response, signaling, cell cycle, aggregation. Xenium is silent.
3. **Untargeted discovery** — Xenium needs a pre-designed panel; bulk LC-MS finds what's there. Drug target ID, biomarker discovery, stress response. DE-LIMP's proteogenomics workflow finds novel ORFs/isoforms NOT on any panel.
4. **Secretome and ECM** — Xenium reads intracellular transcripts; bulk reads secreted factors, ECM proteins, growth factors, complement, cytokines. Critical for tumor / inflammation / fibrosis research.
5. **Proteoforms** — splice variants, truncations, fusion proteins, neoantigens. Transcript splicing ≠ stable protein outcome.
6. **Cell-type-call validation** — Xenium classifies cells from transcript markers. Phenotypic discordance ("transcript+, protein−") is biologically real and reveals translation-level regulation.
7. **Absolute quantification** — Xenium gives per-cell relative counts; bulk LC-MS provides ng-per-mg-tissue. Required for biomarker thresholds, dosing, clinical translation.
8. **Turnover dynamics** — pulse-SILAC bulk gives half-lives. Xenium is a snapshot.
9. **Cell-of-origin for serum biomarkers** — Xenium says which cells transcribe; matched serum proteomics says which proteins left. Together: cell-of-origin attribution for circulating biomarkers.

---

## Lit review (April–May 2026)

The published proof-of-concept work is mostly **Xenium + IMC/CODEX (antibody-based imaging proteomics)**, not Xenium + bulk LC-MS. Some key papers:

- **Greenwald et al. (Cell 2024)** — *Integrative spatial analysis reveals a multi-layered organization of glioblastoma*. doi:10.1016/j.cell.2024.03.029. Combines spatial transcriptomics + spatial proteomics. Three modes of cellular organization defined by joint analysis.
- **Massoni-Badosa et al. (Immunity 2024)** — *An atlas of cells in the human tonsil*. doi:10.1016/j.immuni.2024.01.006. Five modalities integrated: scRNA + epigenome + proteome + immune repertoire + spatial transcriptomics. 556,000 cells.
- **Pan-cancer CAF atlas (Cancer Cell 2025)** — *Conserved spatial subtypes and cellular neighborhoods of cancer-associated fibroblasts revealed by single-cell spatial multi-omics*. doi:10.1016/j.ccell.2025.03.004. 14 million cells, 10 cancer types, 7 spatial platforms. Establishes that the highest-impact spatial cancer atlases are co-pairing transcriptomics with proteomics.
- **Lee et al. (Angew Chem Int Ed 2025)** — *Integrating Ambient Ionization Mass Spectrometry Imaging and Spatial Transcriptomics on the Same Cancer Tissues to Identify RNA-Metabolite Correlations*. doi:10.1002/anie.202502028. Direct MS imaging + spatial transcriptomics on the same section. The closest published model for what bulk-LC-MS-friendly integration could look like.
- **Lhumeau et al. (Cell Syst 2024)** — *PanIN and CAF transitions in pancreatic carcinogenesis revealed with spatial data integration*. doi:10.1016/j.cels.2024.07.001. Imaging + ST + scRNA-seq pipeline.
- **iCCA spatial atlas (Hepatology 2024)** — IMC + spatial proteomics + spatial transcriptomics. Cancer microenvironment.

**The gap:** bulk LC-MS proteomics is rarely formally integrated with spatial transcriptomics in published workflows. The cell-resolution mismatch (Xenium subcellular, bulk LC-MS whole-tissue) makes this an open research-software opportunity rather than a solved problem.

For scRNA-seq + bulk proteomics specifically:
- **Stewart et al. (Cell 2024)** — bone marrow niche atlas using scRNA + CODEX proteomic imaging. doi:10.1016/j.cell.2024.04.013
- **Han et al. (Cell Metab 2024)** — osteokine atlas integrating bone proteomics + scRNA-seq. doi:10.1016/j.cmet.2024.03.006
- Standard tools: MOFA+, BayesPrism, MuSiC, SCDC, CIBERSORTx (cell-type deconvolution).

---

## How the data structures fit together

| | Xenium | Bulk LC-MS |
|---|---|---|
| Unit of measurement | Single cell (x, y coordinate) | Whole tissue section |
| Quantification | Transcript counts per cell | Protein/peptide abundance (log) |
| Identity | 300–5,000 genes (chosen panel) | 3,000–8,000 proteins (unbiased) |
| Modifications | None | Phospho/glyco/acetyl |
| Statistical unit | One section | N samples per condition |
| What it answers | "Which cells transcribe what, where?" | "What proteins are present, and modified how?" |

**Join points:**
- Gene → Protein identity (most genes map 1-1 to a UniProt protein; isoforms complicate this)
- Sample identity (same tissue block, ideally adjacent serial sections cut within minutes, or LCM same-section sequential workflows)

---

## Xenium output format (what DE-LIMP would ingest)

Xenium ships a bundle directory (~5–50 GB per slide). Key files:

| File | Format | Contents | DE-LIMP needs it? |
|---|---|---|---|
| `cell_feature_matrix.h5` | HDF5 | Sparse cells × genes count matrix | **Yes** — primary input |
| `cells.parquet` | Parquet | Per-cell coords, area, total counts | Yes (for spatial patterns) |
| `transcripts.parquet` | Parquet | Per-transcript coords + cell ID + quality | No (too granular for v1) |
| `cell_boundaries.parquet` | Parquet | Cell segmentation polygons | Optional (for visualization) |
| `nucleus_boundaries.parquet` | Parquet | Nucleus polygons | No |
| `morphology.ome.tif` | OME-TIFF | DAPI + multi-channel imaging | No |
| `analysis/clustering/*.csv` | CSV | Pre-computed cell-type clusters | Optional (use 10X defaults) |
| `analysis/umap/*.csv` | CSV | UMAP coordinates | Optional |
| `gene_panel.json` | JSON | Panel gene list | Yes (for matching against bulk) |
| `metrics_summary.csv` | CSV | QC | Yes (for UI display) |

**Loaders already exist:**
- R: `Seurat::LoadXenium("path/to/bundle/")` → Seurat object
- R: `SpatialFeatureExperiment::read10xXeniumSFE()` → Bioconductor flavor
- Python: `squidpy.read.xenium(...)`, `spatialdata_io.xenium(...)`

DE-LIMP would use the Seurat loader (R-native).

---

## Five integration patterns, ranked by effort and value

### Pattern 1 — mRNA-protein concordance check ⭐ build first

**Question:** Which genes have discordant transcript vs protein levels in your samples?

**Approach:**
- Aggregate Xenium counts per sample (sum across all cells in section)
- Match to bulk LC-MS protein abundance by gene/protein identity
- Compute per-gene Pearson correlation across the sample cohort
- Rank genes by translation-discordance metric

**Output:** scatter plot of mRNA vs protein per gene; ranked table of "translationally regulated" candidates; per-sample concordance summary.

**Biological value:** lights up post-transcriptional regulation, protein degradation, secretion. The simplest, most defensible result.

**Effort:** 2–3 days. No spatial coordinates needed for this one.

---

### Pattern 2 — Cell-type signature transfer / bulk deconvolution

**Question:** What cell-type composition explains each bulk-proteomics sample?

**Approach:**
- Xenium → per-cell-type expression signatures (means by cluster)
- Bulk LC-MS abundance matrix + signatures → cell-type proportions via BayesPrism or MuSiC
- Outputs: per-sample proportion vector; cell-type-specific protein abundance estimates

**Output:** stacked bar chart of cell composition per sample; per-cell-type-protein heatmap.

**Biological value:** gives bulk samples a spatial/cell-type identity they didn't have. Useful for tumor purity, immune infiltrate quantification.

**Effort:** 1 week. Tooling exists (BayesPrism, MuSiC are R packages).

---

### Pattern 3 — Pseudo-spatial proteomics

**Question:** Where in the tissue does each detected protein localize?

**Approach:** reverse direction from Pattern 2. Bulk LC-MS protein abundance + Xenium cell-type-by-location map → back-infer spatial protein distribution.

For each protein:
1. Find Xenium-detected cell types that express the cognate gene
2. Weight by where those cells live (x, y)
3. Generate a heatmap

**Output:** spatial heatmap per protein on the tissue section.

**Biological value:** visually striking, but methodologically uncertain when transcript-protein correlation is weak. Best for highly-correlated proteins.

**Effort:** 2 weeks. Requires plotting infrastructure for tissue overlays.

---

### Pattern 4 — PTM → kinase → cell-type triangulation

**Question:** Which cells are driving a given phospho-signaling pathway?

**Approach:**
- Bulk phospho → KSEA → upstream kinases
- Xenium → cells expressing those kinases
- Output: cell-type-resolved kinase activity map

**Output:** kinase-by-cell-type matrix; spatial map of inferred kinase activity.

**Biological value:** mechanistic — answers "which cells drive this signal" that neither modality alone can answer.

**Effort:** 2 weeks. Phospho already in DE-LIMP; KSEA available.

---

### Pattern 5 — Proteogenomics × spatial isoform validation

**Question:** Are novel proteoforms detected by bulk LC-MS expressed in specific cell types?

**Approach:**
- DE-LIMP-built proteogenomics FASTA detects novel ORFs/isoforms in bulk
- For each, check if Xenium panel has probes spanning the novel splice junction
- Cross-validate at cell-type level

**Output:** novel-isoform table annotated with cell-type expression from Xenium.

**Biological value:** new science territory. Pair with the existing proteogenomics workflow (v3.11.0).

**Effort:** 3 weeks. Requires Xenium panel probe-coordinate matching.

---

## Recommended roadmap

**Phase 1 (1 week):** Build Pattern 1 (mRNA-protein concordance) as a new "Spatial × Bulk" tab in DE-LIMP. Minimum viable integration:
- File picker for Xenium bundle directory
- Load `cell_feature_matrix.h5` via Seurat
- Aggregate to per-sample pseudobulk counts
- Match to existing DE-LIMP bulk proteomics results by gene/UniProt mapping
- Concordance scatter plot, ranked discordance table

**Phase 2 (1 week):** Add Pattern 2 (cell-type deconvolution) — BayesPrism integration. This is the most published-precedent pattern.

**Phase 3 (deferred):** Patterns 3, 4, 5 once Phase 1+2 are validated on real data.

**Pilot project idea:** before any of this is built, take one of the Genome Center's existing Xenium datasets, do matched bulk proteomics on adjacent serial sections, and demonstrate Pattern 1 manually (analyst-driven, not productized). That's the "Xenium + LC-MS together tells you something neither does alone" demo for the Genome Center director. From there, productize in DE-LIMP.

---

## What this means for the Proteomics Core

The pitch to the Genome Center director:

> "Xenium tells you what cells are planning to do. Bulk LC-MS tells you what they actually did, including:
> - Modifications (the regulatory layer)
> - Stable proteoforms (not transcripts)
> - Secreted output (what the rest of the body sees)
> - Drug targets and binding partners
>
> The highest-impact Xenium papers (CAF atlas, glioblastoma, tonsil atlas) all pair spatial transcriptomics with proteomics. UC Davis Genome Center can position itself as offering the complete picture by routing Xenium projects through the Proteomics Core for paired LC-MS — and the new DE-LIMP integration tab will make the joint analysis as easy as the single-modality one."

---

## Open questions / TODO

- [ ] Confirm Xenium output format details against an actual delivered bundle (probably ask 10X / DNA Tech Core for a sample dataset)
- [ ] Decide on Seurat vs Bioconductor SpatialFeatureExperiment as the in-app data structure (Seurat is more widely known; SFE integrates better with limpa)
- [ ] BayesPrism vs MuSiC for cell-type deconvolution — pick one
- [ ] Talk to the DNA Tech Core about whether Xenium runs come with matched bulk LC-MS sample availability (serial sections, LCM)
- [ ] Identify a willing collaborator with paired Xenium + bulk LC-MS data for the pilot
- [ ] Add `xenium_bundle_dir` to the Search History schema so we remember which Xenium dataset was paired with which DIA-NN run

---

*Last updated: 2026-05-26 (v3.11.0 cycle)*
