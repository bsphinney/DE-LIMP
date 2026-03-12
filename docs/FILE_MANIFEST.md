# DE-LIMP File Manifest

All files created by DE-LIMP, organized by location and purpose.

## Per-Dataset Files (in raw data directory)

These files live alongside your raw mass spec data. Shared across all lab members who access the same directory.

| File | Purpose | Created When |
|------|---------|-------------|
| `.delimp_tic_cache.rds` | Cached TIC traces and QC metrics | After TIC extraction (Extract TIC button) |

## Per-Search Files (in output directory)

Created by DIA-NN searches. Each search gets its own output directory.

| File | Purpose | Created When |
|------|---------|-------------|
| `session.rds` | Full analysis state (raw data, DE fits, GSEA, MOFA2, metadata, excluded files) | After pipeline completion (auto-saved) |
| `search_info.md` | Search metadata — parameters, job IDs, file list, log paths | At job submission |
| `report.parquet` | DIA-NN quantification output | When DIA-NN search completes |
| `*.sbatch` | SLURM submission scripts (HPC only) | At job submission |
| `file_list.txt` | Raw file paths for parallel search | At parallel job submission |
| `submit_all.sh` | Launcher script chaining SLURM steps | At parallel job submission |
| `logs/` | Directory containing SLURM stdout/stderr or local execution logs | At job submission |
| `logs/diann_*.out` | SLURM job stdout | During search execution |
| `logs/diann_*.err` | SLURM job stderr | During search execution |
| `logs/diann_*.log` | Local/Docker execution log | During search execution |
| `quant_step2/` | Per-file quant output from parallel Step 2 | During parallel search |
| `quant_step2_orig/` | Backup of Step 2 quant files (before Step 3 overwrites) | During parallel Step 3 |
| `quant_step4/` | Final per-file quant output from parallel Step 4 | During parallel search |
| `*.tsv`, `*.parquet` | Other DIA-NN outputs (matrices, stats, libraries) | When DIA-NN search completes |

## User Config Files (`~/.delimp_*`)

App-level configuration and tracking. Per-user, not shared. These are NOT tied to any specific dataset.

| File | Purpose | Created When |
|------|---------|-------------|
| `~/.delimp_activity_log.csv` | Unified activity log — search submissions, completions, analysis runs, session loads (33 columns) | First search submission or analysis |
| `~/.delimp_job_queue.rds` | Persistent job queue — tracks all DIA-NN jobs across app restarts | First job submission |
| `~/.delimp_speclib_cache.rds` | Registry of predicted spectral libraries with metadata (FASTA, params, path). Migrates to shared volume when available. | First library-free search |
| `~/.delimp_cluster_usage_history.csv` | Historical cluster resource snapshots (CPU/memory utilization over time) | When SSH connected and cluster monitor active |
| `~/.delimp_job_wait_log.csv` | Per-job queue wait times (queued → running transition) for grant reporting | When a queued job starts running |
| `~/.delimp_per_user_usage.csv` | Per-user CPU/memory usage snapshots across lab members | When cluster monitor polls (every 60s) |
| `~/.delimp_lab_members.json` | List of lab member SLURM usernames for per-user resource monitoring | Manually configured |

## Shared Volume Files (optional)

When a shared volume is mounted (e.g., `/Volumes/proteomics-grp/` or `/quobyte/proteomics-grp/`), some files prefer shared locations:

| File | Purpose | Fallback |
|------|---------|----------|
| `{shared}/dia-nn/speclib_cache.rds` | Shared spectral library cache so all lab members benefit from cached predicted libraries | `~/.delimp_speclib_cache.rds` |
| `{shared}/dia-nn/fasta_library/` | Shared FASTA library directory with `catalog.rds` index | `~/.delimp_fasta_library/` |

## Claude/Gemini Export ZIP Contents

Generated on-demand via the Export button. Not persisted on disk.

| File | Purpose |
|------|---------|
| `PROMPT.md` | Full AI analysis prompt with inline data summaries |
| `DE_Results_Full.csv` | Complete DE statistics for all proteins across all contrasts |
| `Expression_Matrix.csv` | Log2 expression values for all proteins × samples |
| `QC_Metrics.csv` | Per-sample QC (precursor/protein counts, MS1 signal) |
| `Sample_Groups.csv` | Group assignments from the metadata table |
| `Search_Parameters.csv` | DIA-NN search settings |
| `GSEA_Results.csv` | Gene set enrichment results (if GSEA was run) |
| `Phospho_DE_Results.csv` | Phosphosite DE results (if phospho detected) |
| `Instrument_Metadata.csv` | Instrument model, m/z range, LC info |
| `TIC_QC_Metrics.csv` | Per-run chromatography QC metrics |
| `Excluded_Files.csv` | Files excluded from analysis with reasons, groups, and user notes |
| `Methods_and_References.txt` | Full methodology text with citations |
| `Reproducibility_Code.R` | R code log reproducing every analysis step |
| `Session_State.rds` | Complete session state for exact reproduction |

## Comparator Claude ZIP Contents

Generated from the Run Comparator tab. Not persisted on disk.

| File | Purpose |
|------|---------|
| `claude_prompt.md` | Comparator-specific AI analysis prompt |
| `settings_diff.csv` | Parameter comparison between runs |
| `protein_universe.csv` | Protein overlap (shared/A-only/B-only) |
| `de_results_combined.csv` | Side-by-side DE statistics for both runs |
| `discordant_proteins.csv` | Proteins with disagreeing DE calls + hypotheses |
| `comparison_context.md` | Tool context and Gemini analysis (if generated) |
| `excluded_files.csv` | Excluded files (if any) |
| `spectronaut_run_qc.csv` | Per-sample QC from Spectronaut (Mode B only) |
| `spectronaut_library_info.csv` | Spectronaut library stats (Mode B only) |
| `search_params_run_a.md` | DIA-NN search parameters for Run A (if available) |

## Design Principles

1. **Derived data stays with source data** — cached/computed results (TIC cache, session.rds) go in the raw data or output directory, not the user's home. This ensures any lab member scanning the same directory gets cached results without re-computation.

2. **App config in `~/.delimp_*`** — user-level settings and tracking (activity log, job queue, cluster usage) belong in home directory. These are per-user and not tied to any specific dataset.

3. **No mounted drive dependency** — all app state uses local `~/.delimp_*` paths. SMB mounts may be absent, slow, or disappear. Shared storage is preferred for spectral library and FASTA caches but always has a local fallback.

4. **Hidden dotfiles for caches** — cache files use `.delimp_` prefix (hidden on Unix) to avoid cluttering directory listings.
