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
- [ ] **End-to-end Docker testing**: Test full Docker submit → monitor → auto-load flow with real data
- [ ] **Thermo .raw TIC extraction**: Extend chromatography QC to Thermo files
- [ ] **XIC viewer over SSH**: Currently requires local file access. Need SCP download of `_xic/*.xic.parquet` files from HPC. Large files (100+ MB/sample) — consider streaming or on-demand per-protein download.

## DPC-Quant Detection Transparency (per statistician review)
- [ ] **Expression Grid tooltips**: Hover shows nObs, SE, 95% CI per cell (zero risk, high value)
- [ ] **Detection_Class export column**: Complete/Partial/Sparse/Inferred based on nObs across samples
- [ ] **Expression Grid saturation overlay**: Toggle (off by default) — opacity scales with nObs/maxNobs. NOT red/green (implies "bad"). Neutral visual: saturation/opacity only.
- [ ] **Volcano "detection-driven DE" markers**: Triangle shape for proteins with nObs=0 in all samples of one condition. These are DE calls driven by the detection probability model — scientifically interesting.
- [ ] **Violin plot hollow markers**: Open circles for nObs=0 estimates with tooltip "inferred from detection pattern"
- ~~Evidence Score 0-100~~: REJECTED — double-counts info (SE already incorporates nObs). Use SE directly if single number needed.
- ~~Filterable high-confidence subset~~: REJECTED — contradicts DPC-Quant's design. At most export-only option with warning.

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

## General
- [ ] Grid View: Open violin plot on protein click with bar plot toggle
- [ ] Publication-quality plot exports (SVG/PNG/TIFF with size controls)
- [x] Sample correlation heatmap (Replicate Consistency tab)
- [x] Venn diagram of significant proteins across comparisons (→ Run Comparator protein universe)
- [ ] Sample CV distribution plots
- [ ] Protein numbers bar plot per sample
- [ ] Absence/presence table for on/off proteins
