# DE-LIMP TODO

## Phosphoproteomics — Phase 2 (Kinase Activity & Motifs)
- [ ] **KSEA integration** (`KSEAapp` CRAN package): Infer upstream kinase activity from phosphosite fold-changes using PhosphoSitePlus + NetworKIN database. Horizontal bar plot of kinase z-scores.
- [ ] **Sequence logo / Motif analysis** (`ggseqlogo` CRAN package): Extract ±7 flanking residues around significant phosphosites, display as sequence logos. Requires FASTA upload.
- [ ] **Kinase Activity tab** in phospho results navset: Run KSEA button, bar plot, results table
- [ ] **Motif Analysis tab** in phospho results navset: Logos for up/down regulated sites
- [ ] Dockerfile: Add `KSEAapp`, `ggseqlogo` to CRAN install list

## Phosphoproteomics — Phase 3 (Advanced)
- [ ] **Protein-level abundance correction**: Subtract protein logFC from phosphosite logFC
- [ ] **PhosR integration** (Bioconductor): RUVphospho normalization, kinase-substrate scoring
- [ ] **AI context for phospho**: Append phosphosite DE results and KSEA kinase activities to Gemini chat
- [ ] **Phospho-specific FASTA upload**: Map peptide-relative positions to protein-relative positions

## MOFA2 — Next Steps
- [ ] **MEFISTO integration**: Temporal/spatial MOFA for time-course experiments
- [ ] **Factor annotation**: Link factors to GO terms based on top weights
- [ ] **DIA-NN report processing**: Process raw DIA-NN .parquet as MOFA view via existing pipeline
- [ ] **Dockerfile**: Add MOFA2 + basilisk to Docker image

## Core Facility Mode — Next Steps
- [ ] **QC run ingestion**: Auto-record QC metrics when loading HeLa digest report.parquet
- [ ] **Report template polish**: Add GSEA section, MOFA variance explained, configurable logo/header
- [ ] **Report comparison**: Side-by-side QC bracket + DE summary for two reports
- [ ] **HF state upload/download**: Upload `.rds` state to HF Spaces for shareable live links
- [ ] **Template application on search submit**: Auto-apply saved search preset
- [ ] **Audit log**: Track who generated which report, when, with what parameters
- [ ] **Multi-instrument QC alerts**: Flag instruments where protein count drops below rolling mean - 2*SD
- [ ] **End-to-end testing**: Test full flow with real DIA-NN search → QC ingest → report generation

## DIA-NN Search
- [x] **Shared speclib cache**: Move `~/.delimp_speclib_cache.rds` to shared volume (`/Volumes/proteomics-grp/dia-nn/`) so all lab members benefit from cached predicted libraries. Fall back to local home dir if volume not mounted.
- [x] **NCBI proteome download**: Download FASTA from NCBI with gene symbol mapping via E-utilities (v3.7)
- [x] **SSH file browser**: Visual directory browser for remote mode (v3.7)
- [x] **Load from HPC**: One-click download and load of completed search results (v3.7)
- [x] **No-replicates mode**: Quantification completes, DE skipped gracefully (v3.7)
- [ ] **Job queue GUI accuracy after retries**: When parallel jobs are retried/resumed, the step progress display shows "0/54, 0%" because it tracks the new job entry (which has no completed tasks) instead of aggregating original + retry. Old failed jobs should be hidden or marked "superseded". The progress counter should sum completed tasks from both original and retry runs.
- [ ] **End-to-end Docker testing**: Test full Docker submit → monitor → auto-load flow with real data
- [ ] **Thermo .raw TIC extraction**: Extend chromatography QC to Thermo files
- [ ] **XIC viewer over SSH**: Currently requires local file access. Need SCP download of `_xic/*.xic.parquet` files from HPC. Large files (100+ MB/sample) — consider streaming or on-demand per-protein download.

## DPC-Quant Detection Transparency (per statistician review)
- [x] **Expression Grid tooltips**: Hover shows nObs, SE, 95% CI per cell via JS rowCallback with hidden columns (v3.7)
- [x] **Detection_Class export column**: Detected_All/Detected_Partial/Inferred_All based on nObs across samples. `compute_detection_class()` in helpers.R, used in Expression Grid CSV, session export, and AI export (v3.7)
- [ ] **Expression Grid saturation overlay**: Toggle (off by default) — opacity scales with nObs/maxNobs. NOT red/green (implies "bad"). Neutral visual: saturation/opacity only.
- [ ] **Volcano "detection-driven DE" markers**: Triangle shape for proteins with nObs=0 in all samples of one condition. These are DE calls driven by the detection probability model — scientifically interesting.
- [x] **Violin plot hollow markers**: Open circles (shape 21) for nObs=0 inferred estimates, filled circles for detected. Subtitle shows count of inferred values (v3.7)
- ~~Evidence Score 0-100~~: REJECTED — double-counts info (SE already incorporates nObs). Use SE directly if single number needed.
- ~~Filterable high-confidence subset~~: REJECTED — contradicts DPC-Quant's design. At most export-only option with warning.

## Cascadia De Novo Sequencing
- [x] **R integration (Phases 2-4)**: SSL parsing, peptide classification, DIAMOND BLAST, sbatch generation, server module, UI (feature/cascadia-denovo branch)
- [x] **Bruker native loader**: `bruker_augment.py` using timsrust_pyo3 for native `.d` file reading (~2-5 min vs 45-90 min mzML conversion)
- [x] **Own navbar dropdown**: Moved from DE Dashboard sub-tab to top-level "De Novo" dropdown
- [ ] **Per-residue amino acid coloring**: Color each amino acid in the sequence column by its confidence probability. Requires modifying Cascadia's output to export per-residue softmax probabilities from the transformer beam search (not in SSL format currently).
- [ ] **Submit Cascadia from GUI**: Wire up the "Submit Cascadia Job" tab with SSH job submission, conda env path, model checkpoint path, GPU partition selection
- [ ] **Cascadia routing patch**: Apply the 6-line routing change to Cascadia's `cascadia.py` on HIVE to auto-detect `.d` files and use `bruker_augment.py`
- [ ] **End-to-end test**: Run Cascadia on the same `.d` files as a DIA-NN search, load both results, verify cross-referencing works
- [ ] **DIAMOND BLAST integration**: Wire up the Run DIAMOND button with `module load diamond` on HIVE

## Run Comparator
- [x] **Spectronaut 20+ RunOverview format**: Key-value pair format (Parameter/Value columns) now auto-detected alongside older wide-table format in `parse_spectronaut_run_summaries()` (v3.7)

## Contaminant Tracking & Benchmarking
- [ ] **Contamination level database**: Record per-sample contaminant % in activity log on every pipeline run. Build reference distribution across all analyses (percentiles).
- [ ] **Benchmarking badge**: After pipeline, show "Your contaminant level (2.1%) is in the 35th percentile of all samples processed" — green/yellow/red badge.
- [ ] **Core Facility QC report section**: Add contaminant summary to generated reports. Flag samples above 90th percentile.
- [ ] **Instrument-specific baselines**: Track contaminant levels per instrument (from instrument_metadata). Different instruments have different typical contamination.
- [ ] **Keratin trend monitoring**: Track keratin contamination over time to detect sample prep workflow degradation.

## Data Explorer
- [x] **Abundance Profiles (Quartile Analysis)**: Heatmap of top 10 proteins per intensity quartile with per-sample consistency (v3.7)
- [x] **Sample-Sample Scatter**: Pairwise comparison with correlation, outlier labeling, contaminant overlay (v3.7)
- [ ] **History download**: Download .rds session files from History tab for sharing with collaborators

## Deployment (v3.7 — Complete)
- [x] **Docker launcher for Windows**: `Launch_DE-LIMP_Docker.bat` with shared PC support
- [x] **SSH auto-connect**: Auto-connect to HPC on startup when SSH key detected
- [x] **Environment badge**: Colored badge showing Docker/HPC/Local/HF mode
- [x] **SLURM proxy for Apptainer**: All 9 command paths proxied
- [x] **Shared HPC storage**: All files on `/quobyte/proteomics-grp/de-limp/`
- [x] **Per-user HPC directories**: Multi-user support without conflicts
- [x] **Container detection**: Skip BiocManager validation offline
- [x] **Home directory quota warning**: Startup check for HPC quota limits

## CV Analysis Tab Redesign (Complete)
- [x] Replace broken DT table with plotly scatter plot (logFC vs Avg CV, color-coded by CV category)
- [x] Add Avg CV (%) column to DE Results Table (inline computation, no reactive dependency)
- [x] Simplify CSV export (removed toggle filter, exports all significant proteins)
- [x] Update info modal for new design (scatter plot, summary stats, Results Table column)
- [x] **Fix summary stats cards**: Replaced fragile plotly annotation cards with ggplot subtitle (per-group median CV + % below 20%)
- [x] **Fix scatter plot compression**: Wrapped CV Analysis tab in scrollable div with min-height on scatter plot container

## Volcano Plot Fixes (Complete — v3.1.1)
- [x] Fix P.Value vs adj.P.Val mismatch: y-axis raw P.Value, dashed line at FDR-equivalent threshold
- [x] Color significance by adj.P.Val only (not logFC cutoff) — logFC lines are visual guides
- [x] Add DE protein count annotation ("78 DE proteins (X up, Y down)")
- [x] Default logFC cutoff changed from 1.0 (2FC) to 0.6 (~1.5FC)

## Publication Export (per biological researcher & proteomics expert review)
- [ ] **Vector figure export (SVG/PDF)**: All plots (Volcano, Heatmap, PCA, CV, GSEA) need SVG/PDF export for publication-quality figures. PNG is raster and blurry at journal scale.
- [ ] **Excel workbook export (.xlsx)**: Single workbook with multiple sheets (DE, CV, GSEA, Contaminant, Metadata) for researcher convenience
- [ ] **Customizable figure sizing**: User-specified dimensions for journal column widths

## Documentation & Education (per biological researcher review)
- [ ] **Glossary tab**: In-app definitions for logFC, adj.P.Val, CV, FDR, DPC-Quant, mass accuracy, ppm. Link to external resources.
- [ ] **DPC-Quant methodology documentation**: Explain distributional assumptions, when imputation kicks in, interaction with limma eBayes. Link to limpa vignette.
- [ ] **Volcano P.Value vs adj.P.Val explanation**: Info modal clarifying y-axis uses raw P.Value for spread, coloring uses FDR-adjusted threshold
- [ ] **GSEA mapping efficiency display**: Show "Mapped X/Y proteins (Z%)" after bitr(), warn if <80% mapped
- [ ] **Power calculation**: Post-hoc display: "With n=3, minimum detectable FC = X at 80% power"

## Biology Features (per biological researcher review)
- [ ] **Protein-protein interaction networks**: Query STRING/BioGRID for top DE proteins, visualize with igraph
- [ ] **Subcellular localization overlay**: Fetch UniProt compartment annotations, add to results table
- [ ] **Multi-contrast biomarker panel**: Find proteins consistently DE across multiple contrasts
- [ ] **Batch effect warning**: Auto-detect when all samples of one group ran on same date, flag with warning

## Core Facility Enhancements (per proteomics expert review)
- [ ] **Real-time QC dashboard**: 30-day rolling plots (proteins, signal per instrument) with outlier detection
- [ ] **Search parameters in reports**: Add enzyme, mass_acc, normalization, DIA-NN version to report metadata
- [ ] **Instrument-specific QC baselines**: Rolling 30-day median per instrument, flag outliers
- [ ] **GSEA contrast-specific visualization**: Lollipop plot of top pathways per contrast with NES comparison

## Statistical Transparency (per statistician review)
- [ ] **Uncertainty quantification in Expression Grid**: Add SE and 95% CI columns (from DPC-Quant posterior)
- [ ] **No-replicates warning banner**: Show "No statistical inference possible" on Expression Grid when in no-replicates mode
- [ ] **Comparator Rule 3 (Quant3) quantitative threshold**: Add explicit statement about t-statistic inflation factor
- [ ] **Comparator Rule 4 (Variance) threshold**: Define what SD ratio constitutes "mismatch"

## Automation
- [ ] **Nightly documentation GitHub Action**: Auto-generate daily changelog summary from git commits. Runs at 9 PM Pacific, updates CHANGELOG.md if new commits exist, commits and pushes. Replaces session-only Claude Code cron which dies on terminal close.

## General
- [ ] Grid View: Open violin plot on protein click with bar plot toggle
- [x] Sample correlation heatmap (Replicate Consistency tab)
- [x] Venn diagram of significant proteins across comparisons (→ Run Comparator protein universe)
- [ ] Sample CV distribution plots
- [ ] Protein numbers bar plot per sample
- [ ] Absence/presence table for on/off proteins
