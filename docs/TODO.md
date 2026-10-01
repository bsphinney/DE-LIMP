# DE-LIMP TODO

## Proteogenomics (v3.11.0 shipped — these are follow-ups)
- [x] **Auto-`assemble` SLURM stage** chained to every new build (v3.11.0)
- [x] **`submit_assemble_only()` + per-row Assemble button** for legacy builds (v3.11.0)
- [x] **UniProt + NCBI download integration** with auto-submit-after-download (v3.11.0)
- [x] **FASTA library auto-registration** on assemble completion (v3.11.0)
- [x] **"Restore from Hive" + "Discover from Hive"** for multi-user catalog (v3.11.0)
- [x] **Active-builds persistence** at `~/.delimp_proteog_builds.rds` (v3.11.0)
- [x] **Explain this workflow modal** in the green header (v3.11.0)
- [x] **Source-tag column + Last-polled column** in Proteog DBs / Active Builds (v3.11.0)
- [x] **Reference genome builder script** + 5 additional species (pig + rat done; bovine, arabidopsis, maize rerunning 2026-05-27) — merge `registry_pending → registry.json` via `references/scripts/merge_registry_pending.sh` when done.
- [ ] **NCBI gene_map.tsv at result-load time**. Currently uploaded alongside the FASTA but not yet wired into DE-LIMP's `Genes` column population in `server_data.R`. Without it, NCBI-derived proteins in a proteog DB show as accessions instead of gene symbols in DE results. Touches the FASTA-resolution path that loads results from a search.
- [ ] **search_settings on proteog catalog entries**. Currently NULL — picking a proteog DB doesn't pre-configure DIA-NN enzymes/mods like Database Library entries do. Less important since proteog FASTAs work with default DIA-NN settings.
- [ ] **End-to-end test** of the auto-assemble pipeline. All 4 legacy builds were submitted pre-v3.11.0. Need to submit one fresh build (small SRA test like the smoketest) and watch it go through all 11 stages → catalog entry → main-page picker without manual intervention.

### Nanopore long-read RNA-seq track for proteogenomics (requested 2026-06-01)
Add a long-read (nanopore cDNA / direct-RNA) alternative to the current short-read STAR track for building the custom proteogenomics FASTA — long reads resolve full-length isoforms → better ORFs → better novel peptides. Parallel track, not a tweak: **minimap2** (splice-aware) instead of STAR → isoform assembly (**IsoQuant / FLAIR / StringTie2 long-read mode**) → ORF calling (TransDecoder) → existing assemble stage / FASTA-library catalog. Read length is very different, so STAR/QC adaptive thresholds (see [[feedback_adaptive_star_thresholds]] / [[feedback_gate_unique_not_total]]) don't apply — long-read QC keys off read-quality + transcript completeness. Prior art: Wilburn lab `NanoporeAnalysis` (github.com/wilburnlab/constellation). Touches helpers_rnaseq.R + the SLURM stage chain in helpers_proteog_assembly.R.

### Novel-peptide validation panes (requested 2026-06-01 — build after de novo work)
Two QC sub-tabs in the proteogenomics section to vet novel/variant peptides (standard proteogenomics validation). Data source: DIA-NN proteog output (peptide → protein → intensity; "novel" = peptide from a custom-FASTA / non-reference entry).
- **Prior art to align with before building:** Wilburn lab (OSU) — github.com/wilburnlab, esp. `constellation` (`Cartographer` MS + `NanoporeAnalysis` cDNA); ASMS 2026 talk "Transcriptome-informed peptide search spaces redirect low-confidence identifications toward biologically supported peptides in DIA" (Edgington, Belyy, Wilburn) — method not yet public as of 2026-06-01. Also Sequoia + SPIsnake (Mol Cell Proteomics 2025, S1535947625001380) for RNA-informed search-space construction.
- [ ] **Protein-level corroboration (priority).** Per protein carrying a novel peptide: # novel vs # canonical (reference) peptides; flag proteins whose ONLY evidence is a single novel peptide as low-confidence (the classic false-positive pattern). Ideally show the novel peptide's position relative to canonical peptides along the protein.
- [ ] **Novel-peptide abundance context.** Novel-peptide intensities vs the whole peptide intensity distribution, with a "highly abundant" reference band. Derive the band from the data's own top-N most-abundant proteins (NOT a hardcoded housekeeping gene list) so it generalizes to non-model organisms. Novel peptides at top-abundance = suspicious; at noise floor = weak.
- [ ] **Bioshare URL support** in Step 1 — SLIMS retired May 30, 2025. **Deferred** because Bioshare requires username/password authentication; not easily automatable.
- [ ] **PacBio Iso-Seq pipeline mode**. Brett asked about it (2026-05-26). Same 11-stage scaffold but stages 1–5 differ: `pbskera` + `isoseq cluster` + `minimap2 -ax splice:hq` + `pbisoseq collapse` instead of fastp → STAR → stringtie. Stages 6–11 (gffcompare, gffread, TransDecoder, rewrite, assemble) stay the same. ~1–2 days dev.
- [ ] **Xenium × bulk LC-MS integration** — see `docs/DESIGN_BRIEF_XENIUM_BULK_INTEGRATION.md`. Pattern 1 (mRNA-protein concordance) is the cheapest demo to build.
- [ ] **scRNA-seq cell-type deconvolution** of bulk proteomics — BayesPrism integration as a new analysis tab. See Xenium brief for context.
- [ ] **Variant calling → patient-specific neoantigen mode** — long-term proteogenomics extension. STAR's WASP-like + GATK + custom OPF format.

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

## DDA / Sage quantification (found by the skill 2.9 Sage review, 2026-09-29; app code not yet changed)
The skill fixed the same issues in `run_search.py adapt_sage` / `sage_lfq_check.py` (see `docs/SKILL_OPEN_DEFECTS.md`).
Sage 0.14.7 facts, from its source:
- `lfq.tsv` has no decoys (`output.rs write_lfq`) but keeps targets at ANY `q_value`.
- An MS1 peak counts at `q_value` ≤ 0.05 (`fdr.rs picked_precursor`, compared in f32).
- LFQ integrates only within ±`lfq_settings.ppm_tolerance`, default 5 ppm (`lfq.rs`).

Items:
- [ ] **`parse_sage_results()` (R/helpers_dda.R) sums every target peptide in `lfq.tsv` with no `q_value` filter.** Keep rows with `q_value` ≤ 0.05 (compare as float32) before the per-protein sum.
- [ ] **Its "legacy long-format" branch** (`proteins`/`filename`/`intensity`, `dcast … fun.aggregate = sum`) checks no `is_decoy`. Fed `lfq.parquet`-shaped data, it would add the decoy MS1 peaks, which carry the target's `proteins` string, to the protein. Drop `is_decoy` there too.
- [ ] **`generate_sage_config()` leaves `lfq_settings.ppm_tolerance` at Sage's 5 ppm and never checks the MS1 offset.** A cohort whose MS1 is +7 ppm off identifies normally and quantifies nothing: gabrig's HeL50 HCDOT kept 0 MS1 peaks, silently. Port the skill's check: the median `precursor_ppm` per run from `results.sage.tsv`, plus Sage's "target MS1 peaks at 5% FDR" count. Warn, and offer a wider window.

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
- [ ] **Skill: a "handling" tag for keratin samples** (keratin-review of skill 2.9, 2026-09-29). `fetch_fasta.py --keratin-sample` removes every keratin-family contaminant entry, so on a non-human appendage (wool, feather, mouse fur) the human skin keratins (KRT1/2/9/10) -- genuine handling contamination there -- become indistinguishable from the sample's keratins. Design: keep those four human entries under a separate tag (not `Cont_`), NOT excluded from quantification by `--cont-quant-exclude` or `run_de.R`, and report their share per run as a handling-contamination QC readout (next to `QC_contaminant_share.csv`). Same question for skin samples and FLG/FLG2/HRNR/DMKN (analytes there, contaminants elsewhere), and for species whose proteome lacks KRTAPs (the contaminant KRTAPs may be the only ones in the search). Documented as limits in SKILL.md step 3 for 2.9.

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
- [x] **Thermo .raw support for Casanovo** (v3.11.2): Sage sbatch now pre-converts .raw → mzML via msconvert apptainer; Casanovo v5 reads mzML directly, so no separate MGF conversion needed for Thermo data. Single conversion reused across both engines.
- [ ] **IM-aware Casanovo**: Same 5th embedding channel as Cascadia IM model — add ion mobility to Casanovo for DDA de novo on timsTOF data.
- [ ] **Mobility-filtered extraction for Cascadia**: Mode B mobilogram peak detection from CASCADIA_MOBILITY_FILTER_ADDENDUM.md — produces cleaner pseudo-DDA spectra from diaPASEF windows. Expected 2-3x improvement in de novo calls.

## DDA dropdown — Peptidomics & Immunopeptidomics modes (next after ocelot ships)
Rename current "De Novo" navbar dropdown → **"DDA"** with sub-entries for the different analysis types. All ride the same Sage + Casanovo + DIAMOND backend (both Thermo .raw and Bruker .d via the existing v3.11.2 chain — Sage's msconvert step for .raw, `bruker_to_mgf.py` for .d). Search params + result framing differ per mode; instrument-agnostic.

**Borrowed from Michael Krawitzky's Ziggy/BOWIE** (numbers + tested observations only — his repos are mostly README + UI shells, no working HLA implementation):
- HLA-I/II preset numbers (length, charge, mods) come straight from BOWIE's published preset table — match standard immunopeptidomics practice.
- `precursor_charge_range = [1, 3]` for HLA on TOF instruments — Michael's `z=1 native` observation; most MHC-I peptides come off singly-charged on TOF.

Skipping (aspirational, no source to crib): MHC-I binding corridor viz, HLA allele inference (Ziggy's "HLA Discovery" — README-only), NetMHCpan integration (license-gated), single-cell HLA. Build only if a user actually asks.

- [ ] **Rename `De Novo` dropdown → `DDA`** in `R/ui.R` navbar
- [ ] **Single shared Run Search form** with `selectInput("dda_mode")` at top — choices: `De Novo Search` (default, current behavior), `Peptidomics`, `HLA / MHC class I`, `HLA / MHC class II`. Mode drives downstream search-param defaults + which results page opens after Load.
- [ ] **Sage preset variants** in `generate_sage_config()`:
  - **Peptidomics**: `enzyme.cleave_at = ""` (nonspecific), `peptide_min_len=5, max_len=25, peptide_min_mass≈400, max_mass≈5000`. Variable mods: pyro-Glu (Q/E N-term), C-term amidation (−0.984), N-term acetylation. Walltime bump → 8h, mem → 128G (nonspecific search is ~10–50× slower than tryptic).
  - **MHC class I**: `enzyme.cleave_at = ""`, `peptide_min_len=8, max_len=12` (BOWIE-compat, slightly broader than 8–11), `peptide_min_mass≈700, max_mass≈1500`. Variable mods: oxidation + deamidation. `precursor_charge_range = [1, 3]` (per Michael's `z=1 native` observation for TOF; safe default for Orbitrap too). Tighter `precursor_tol = ±5 ppm`. Default walltime OK.
  - **MHC class II**: `enzyme.cleave_at = ""`, `peptide_min_len=13, max_len=25`, `peptide_min_mass≈1300, max_mass≈3000`. Variable mods: oxidation + deamidation. `precursor_charge_range = [1, 3]`. Default walltime OK.
- [ ] **Mode-specific results pages** (three new sub-entries in the DDA dropdown). Build from published methods, not from Ziggy:
  - **Peptidomics Results**: peptide-source-protein chart (parent contributions), N- and C-terminal cleavage motif logos (`ggseqlogo`, ~20 LOC), mass distribution histogram, PTM landscape (pyroGlu/amidation/oxidation rates).
  - **HLA/MHC Results**: length distribution histogram (expect sharp 9-mer peak for class I — the diagnostic plot for a clean HLA-I prep), P2 + PΩ anchor residue sequence logos (`ggseqlogo`), source-protein analysis (re-skin of de novo's protein-family panel).
  - **De Novo Results** (existing): species attribution, BLAST alignments, coverage maps, deamidation tracking — unchanged.
- [ ] **Auto-route Load Results → correct page** based on the mode tag stored in `~/.delimp_dda_queue.rds`. User can also manually switch lenses for cross-comparison.
- [ ] **`ggseqlogo` to Dockerfile** for sequence logo rendering (small CRAN dep).
- [ ] **timsTOF-only bonus (HLA mode)**: pull `1/K₀` column from Sage output for .d-data searches and add a "CCS vs peptide length" panel. HLA peptides have tight CCS clustering by charge — quick sanity-check plot. Easy add since the data is already in Sage's output, just an extra ggplot.

### Deferred / build only if a user asks
- NetMHCpan-4.1 binding-affinity scoring (license-gated, ~1 day work to wire on Hive). Adds neoantigen prioritization. Only worthwhile for tumor immunopeptidomics users.
- HLA allele inference (PSSM motif matching to predict donor HLA type from peptidome motifs). Published methods exist (e.g. GibbsCluster, MixMHCpred deconvolution). Not in any of Michael's code despite the marketing — would need to implement from papers.

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

## Code Audit Follow-ups (post-v3.10.8)
Moved out of CLAUDE.md's Version History during the v3.x doc refactor. These are guardrail-driven cleanups (see CLAUDE.md "Architectural rules").
- [ ] **Covariate-name single source of truth** — refactor the 22 read sites of `values$cov1_name %||% "Covariate1"` into one helper (Architectural Rule #3).
- [ ] **`<<-` in `add_covariate()` closure** — refactor away the `<<-` per the `withProgress` gotcha.
- [ ] **`DIANN_DEFAULTS` constant extraction** in `helpers_search.R` — stop scattering default mass-acc/Q-value fallbacks.
- [ ] **`detect_organism_db()` silent-fallback refactor** — currently returns `"org.Hs.eg.db"` as a default with no signal to callers; should return `NULL` or `list(db, method)` so callers (PROMPT.md builder, Run Comparator, Explorer prompt, GSEA tab) handle uncertainty explicitly per Architectural Rule #2. Caused the v3.10.6 → v3.10.7 "Organism: Human" bug on a Peromyscus dataset. Sweep for similar silent-fallback functions (`coalesce_setting`, default Q-value cutoffs, default mass-acc fallbacks).
- [ ] **Convert export bundlers to `safe_section()`** — sweep all export paths for the `if (!is.null(f)) files_to_zip <- c(...)` and `tryCatch(error = function(e) NULL)` patterns (Architectural Rule #4). The v3.10.4 → v3.10.8 hotfix train (5 patches in one day) was almost entirely silent-failure-masquerading-as-success bugs in export paths. Files to audit: `R/server_ai.R` (Claude AI export), `R/server_viz.R` (Explorer export), `R/server_phospho.R` (Phospho export), `R/server_mofa.R` (MOFA export), `R/server_comparator.R` (Comparator export). The Complete Analysis bundle in `R/server_session.R` was made clean in v3.10.8; the others are still likely broken the same way.

## Skill (ucdavis-proteomics-core-pipeline) — backlog
- [ ] **Pulldown set tests: background subtraction in place of the 4-fold minimum (`pooled_bgsub`)** — after 2.9. The 4-fold interactome minimum makes the bait-normalised tests one-sided: a set that falls in one condition keeps 4.5 / 2.9 of 20 members (Old background +0 / +0.6 log2) and reaches the 10-member minimum in 0.9% / 0.0% of replicates, so losses from a complex are close to undetectable (2.9 says so in every reading). Candidate: subtract the condition's mean control (IgG) background on the linear scale, with no minimum enrichment, and floor what falls to or below the background. In rel-stats' simulation (check6, 1,000 reps) a falling set then keeps 19.4 / 19.1 members with fry power 0.98 / 0.96; null random sets 0.065 / 0.085; the 20 least-enriched interactors 0.094 / 0.094. **That simulation uses the same additive background + signal model the method assumes, so it cannot validate it: validate on real pulldowns first** (e.g. PROT_0756, and a spike-in or known-interactor pulldown), and decide the floor there. Details: `skill/ucdavis-proteomics-core-pipeline/references/set-tests.md` "Known limitations".
