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
- [x] **bruker_patch.py MS1 frame reading**: Correct TOF-to-m/z conversion using Sage formula from timsrust `src/converters.rs`. Reads calibration from GlobalMetadata SQLite table. (April 2026)
- [x] **DIA window matching fix**: Use `SpectrumReader.new_with_span_step()` for 174k spectra (7.7x expansion from 21k raw) with correct `precursor_mz` and `isolation_width` per DIA window. (April 2026)
- [x] **b/y ion quality filter**: 10% relative intensity threshold for timsTOF training data. Removes noise peaks while retaining real fragment ions. (April 2026)
- [x] **IM model architecture**: 4th embedding channel for ion mobility (1/K0), zero-init trick for backward compatibility with pre-trained checkpoint. 5-column ASF format. (April 2026)
- [x] **IM ASF pipeline verified**: `new_with_span_step()` produces 174k spectra from single .d file with IM values. All 5 unit tests pass (creation, forward, zero-init match, IM sensitivity, mixed batch). (April 2026)
- [x] **Native .d validation**: 738 peptides from native .d path vs 44 from mzML (16.8x improvement). bruker_patch v3 working. (April 2026)
- [x] **Training pipeline audit**: Verified all hyperparameters against Cascadia source. Key findings: (1) Cascadia does NOT use peak filtering — `max_num_peaks=200` is a depthcharge default, not what the pretrained model expects; (2) `configure_optimizers()` has hidden CosineWarmupScheduler that must be overridden for fine-tuning; (3) batch_size=1 + grad_accum=16 needed for full unfiltered spectra (median 9,558 peaks, max 113k). (April 2026)
- [ ] **Run IM-enhanced training**: Submit IM training (5-column ASF) after baseline fine-tuning validates. Compare IM vs baseline on Zhao validation set.
- [ ] **Compare baseline vs IM model on test data**: Use held-out Zhao dataset (43 .d files). Metrics: peptide count, sequence accuracy, score distribution.
- [ ] **Download and process ddaPASEF pre-training data**: PXD014777 (Prianichnikov 2020, HeLa) and PXD010012 (Meier 2018, HeLa). Clean isolated precursor spectra for timsTOF-specific pre-training.
- [ ] **Propose Noble Lab collaboration**: Working prototype with 738 peptides (16.8x mzML), IM integration, ddaPASEF training pipeline. Demonstrate value of native Bruker support.
- [ ] **Per-residue amino acid coloring**: Color each amino acid in the sequence column by its confidence probability. Requires modifying Cascadia's output to export per-residue softmax probabilities from the transformer beam search (not in SSL format currently).
- [ ] **Submit Cascadia from GUI**: Wire up the "Submit Cascadia Job" tab with SSH job submission, conda env path, model checkpoint path, GPU partition selection
- [ ] **Cascadia routing patch**: Apply the 6-line routing change to Cascadia's `cascadia.py` on HIVE to auto-detect `.d` files and use `bruker_augment.py`
- [ ] **End-to-end test**: Run Cascadia on the same `.d` files as a DIA-NN search, load both results, verify cross-referencing works
- [x] **DIAMOND BLAST integration**: Wire up the Run DIAMOND button with `module load diamond` on HIVE. Auto-BLASTs against SwissProt after Casanovo completion. (April 2026)
- [x] **Sage DDA pipeline (Phase 2)**: Full Sage search pipeline — helpers_dda.R, server_dda.R, SLURM submission, polling, result loading. (April 2026)
- [x] **Casanovo integration (Phase 3)**: Casanovo de novo GPU job, mztab parsing, confirmed/novel classification, DIAMOND BLAST. (April 2026)
- [x] **Casanovo IM architecture**: Precursor-level IM encoding for DDA, zero-init from checkpoint, all 7 tests pass. (April 2026)
- [x] **Mobility-filtered extraction (Mode B)**: Mobilogram peak detection from diaPASEF frames, 5x more spectra with real 1/K0 values. (April 2026)

## De Novo Visualization (from biologist + proteomics expert reviews)
- [x] **Interactive confidence slider**: Draggable threshold (0.5-1.0) that updates all de novo tables, charts, and counts in real time. (April 2026)
- [x] **Contaminant filtering before species charts**: Checkbox to exclude contaminant proteins from species donut/bar. (April 2026)
- [x] **Per-residue confidence heatmap**: Click peptide row → colored sequence bar (green >0.95, yellow 0.7-0.95, red <0.7). Cross-reference with BLAST substitutions. (April 2026)
- [x] **Peptide length/charge distribution QC**: Side-by-side confirmed vs novel histograms. Flag peptides <7 aa or >25 aa, charge 1+. (April 2026)
- [x] **Cross-species comparison table**: Protein × sample matrix with Venn diagram and species heatmap. (April 2026)
- [x] **Protein family grouping**: 16 families (keratins, collagens, histones, etc.). Stacked bar + treemap. (April 2026)
- [x] **BLAST alignment view for near-matches**: Color-coded alignment with per-residue confidence cross-reference. Green=variant, Red=error. (April 2026)
- [x] **Target-decoy BLAST FDR**: Reversed peptide BLAST as SLURM job. FDR curve + hit count plots. (April 2026)
- [x] **Modification tracking (deamidation)**: Parse modification masses from Casanovo sequences. N/Q deamidation ratio for paleoproteomics authenticity. (April 2026)
- [x] **Protein sequence coverage maps**: Top 20 proteins, color-coded by identity (green=100%, orange=90-99%, red=<90%). (April 2026)
- [x] **Manuscript summary statistics card (Table 1)**: Per-sample breakdown with CSV download. (April 2026)
- [x] **GO/pathway annotation**: 11 functional categories from protein name patterns. Bar chart + table. (April 2026)
- [x] **Sage vs Casanovo disagreement analysis**: Match by scan, classify I/L swap, single/multiple substitutions, completely different. (April 2026)
- [x] **Species resolution bar chart**: Identity gap (delta) between best and 2nd-best species per peptide. (April 2026)
- [x] **Taxonomic coverage dot plot**: Peptide identity across species, grouped by protein. (April 2026)
- [x] **Top diagnostic peptides card**: Highest species-resolution peptides for species identification. (April 2026)
- [x] **Info (?) buttons on all sub-tabs**: 12 detailed help modals covering every de novo visualization. (April 2026)
- [x] **SVG export on all plots**: Camera icon for publication-quality SVG download on every plotly chart. (April 2026)
- [ ] **Novel peptide clustering**: Cluster similar novel peptides by sequence similarity. 3+ clustered peptides with same BLAST protein family = candidate novel protein variant. Biologist #5.
- [ ] **FASTA/mzTab export**: Download novel peptides as FASTA (for BLAST2GO, InterProScan), mzTab format for PRIDE/ProteomeXchange deposition. Proteomics expert #9.
- [ ] **Quality flags column**: Traffic light (green/yellow/red) per peptide combining: score, length, BLAST hit, mass error. Biologist #10.
- [ ] **BLAST database selector**: Let user choose SwissProt (~2min) / TrEMBL (~30min) / Custom FASTA with estimated run times.

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

## Data Tracking
- [ ] **Organism column in activity log**: Add `organism` field to activity_log.csv (column 34). Auto-detect from FASTA filename (UniProt suffix `_HUMAN`, `_BOVIN`, `_PIG`, `_MOUSE`; NCBI accession → taxonomy via E-utilities). Store as common name + taxonomy ID. Useful for AI training dataset curation, cross-species analysis, and core facility reporting.
- [ ] **Instrument column in activity log**: Add `instrument` field from `values$instrument_metadata$instrument_model`. Track which instrument produced each dataset for instrument-specific QC baselines.
- [ ] **Species-matched contaminant handling**: When searching bovine sample against bovine FASTA + Universal Contaminants, bovine proteins in the contaminant DB (serum albumin, keratins, trypsin) get `Cont_` prefix but are real endogenous proteins. Options: (a) auto-detect species overlap and strip `Cont_` prefix, (b) add checkbox "Organism matches contaminant species — treat Cont_ as endogenous", (c) warn user in the Contaminant Analysis tab.
- [ ] **Parse search_info.md on Load from HPC**: When loading report.parquet from HPC, also download and parse search_info.md to restore instrument metadata (instrument model, LC system, gradient) and search parameters. Currently Methods text is missing LC/MS info when loading from report.parquet instead of session.rds.
- [ ] **Excluded files in search_info.md**: `values$excluded_files` is tracked in-app and exported in Claude/Comparator ZIPs as `excluded_files.csv`, but NOT written to search_info.md on the HPC. Anyone looking at the output directory can't see why 54/60 files were searched. Add an "Excluded Files" section to `generate_search_info()` listing filenames and reasons.

## Search Performance
- [ ] **Adaptive CPUs on public queue**: When auto-switching step 2/4 to publicgrp/low, increase CPUs from 16 to 64 (public nodes have 128 CPUs). Faster per-file completion reduces preemption risk. Use scontrol update NumCPUs=64 alongside the partition move.
- [ ] **DIA-NN 8-CPU mode**: Test 8 CPUs × 8 concurrent (instead of 16 × 4) on genome-center-grp/high. Could improve total throughput ~30% if DIA-NN scaling is sublinear. Need benchmarks on real data.
- [ ] **Cascadia batch_size=64**: Increase from 32 to 64 on A100 (80GB VRAM). Could halve inference time.

## DDA/DIA Pipeline Refactor
- [ ] **Extract shared UI components**: Create reusable helper functions (`fasta_selector_ui()`, `group_assignment_ui()`, `slurm_status_ui()`) that both DIA and DDA pipelines call. Currently FASTA selection, SSH file browser, SLURM polling, group assignment, and limma DE are partially duplicated between server_search.R/server_data.R (DIA) and server_dda.R (DDA).
- [ ] **Shared FASTA management module**: Single FASTA selector with UniProt/NCBI download, contaminant append, SSH file browser — used by both DIA and DDA search tabs.
- [ ] **Thermo .raw support for Casanovo**: Add ThermoRawFileParser or msconvert pre-conversion step in the Casanovo sbatch. Sage already reads .raw natively, but Casanovo needs MGF input.
- [ ] **IM-aware Casanovo**: Same 5th embedding channel as Cascadia IM model — add ion mobility to Casanovo for DDA de novo on timsTOF data.
- [ ] **Mobility-filtered extraction for Cascadia**: Mode B mobilogram peak detection from CASCADIA_MOBILITY_FILTER_ADDENDUM.md — produces cleaner pseudo-DDA spectra from diaPASEF windows. Expected 2-3x improvement in de novo calls.

## Claude Code Configuration
- [ ] **Restructure CLAUDE.md into `.claude/rules/`**: Move verbose reference material (gotchas table, UI patterns, SSH patterns, Spectronaut parsing, DIA-NN flags, comparator details) into scoped rule files under `.claude/rules/`. Keep CLAUDE.md under ~150 lines with just project overview, architecture, working preferences, key commands. Do on a branch, test by verifying gotcha knowledge in fresh conversation. Per expert consensus: short CLAUDE.md + modular rules > monolithic file.
- [ ] **Add post-edit hook**: Auto-format R files on save (e.g., `styler::style_file`)
- [ ] **Claude export education section**: Ensure the education section is always included in reports generated from Claude export ZIPs. Currently it's in the PROMPT.md template but was missed in the standalone muscle report.

## General
- [ ] Grid View: Open violin plot on protein click with bar plot toggle
- [x] Sample correlation heatmap (Replicate Consistency tab)
- [x] Venn diagram of significant proteins across comparisons (→ Run Comparator protein universe)
- [ ] Sample CV distribution plots
- [ ] Protein numbers bar plot per sample
- [ ] Absence/presence table for on/off proteins
