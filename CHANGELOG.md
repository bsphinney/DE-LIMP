# Changelog

All notable changes to DE-LIMP will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.11.13] — 2026-05-29

### Changed
- **Renamed "Novel Peptides" → "De novo only"** (summary card, results tab, info modal). Clearer than "novel" (which sounds like novel proteins/discoveries) — these are peptides found by Casanovo de novo sequencing but NOT by the Sage database search.
- **DIAMOND BLAST now queries short peptides too** — lowered the de novo peptide extraction length gate from ≥7 aa to ≥5 aa in the sbatch generator (`R/server_dda.R`). 5–6mers (often the top-scoring de novo calls) were silently excluded from the BLAST query; per the "a de novo peptide with no BLAST hit is probably an incorrect call" heuristic, every peptide should be attempted. (DIAMOND's seed length means the very shortest may still return no hit — which is itself the signal.)

### Notes
- Diagnosed the "BLAST Hits: 0" report: ground-truth on disk shows **5,852 of the novel peptides DO have BLAST hits** in the fresh `blast_results.tsv` — the app was displaying STALE blast data because the post-restart reload skipped the SSH-gated blast download. Fix on the user side: reconnect SSH, then Load from HPC. (Follow-up: make the loader warn instead of silently skipping the blast download when SSH is down.)

## [3.11.12] — 2026-05-29

### Added
- **"Protein filter" dropdown on the DDA/de novo Results page** with three modes: *All proteins* (default), *Skin & hair only* (keep keratins/KRTAP/collagen/trichohyalin/filaggrin… — for paleoproteomics / hair / feather work), and *Exclude skin & hair* (treat keratins as contaminants — for normal tissue/cell projects). Applied consistently across the DIAMOND BLAST tab and the Cross-Species / Protein Families / Coverage / Alignment views via one shared helper. Opt-in (off by default) and **name-based only when the user chooses it** — keratins are never filtered automatically. New globals in `R/helpers_dda.R`: `is_skin_hair_protein()`, `dda_apply_protein_filter()`.

### Changed
- **Removed the automatic keratin/name-based contaminant heuristic from the DIAMOND BLAST filter.** For hair/feather/skin de novo work keratins are the signal, not contaminants — dropping them by name was wrong. Contaminant exclusion is now `Cont_`-tag only (the searched contaminant DB); skin/hair removal is user-controlled via the new dropdown.
- **Sage DB hits table + count now drop `Cont_` contaminants** consistently. `filtered_classification()` in `R/server_denovo_controls.R` removes Cont_-only confirmed peptides (and recomputes the confirmed count, confirm rate, and protein list) when the contaminant toggle is on. Uses the shared `is_contaminant_protein_group()` global (Architectural Rule #3).

## [3.11.11] — 2026-05-29

### Added
- **"Exclude contaminants" toggle on the DDA/de novo Results page.** A checkbox in the file-filter panel (default ON) drops PSMs whose entire protein group is `Cont_`-tagged (the contaminant FASTA appended at search time — keratins, trypsin, albumin, etc.). Wired into `dda_filtered_psms()` so every downstream view (Length & Motifs, Peptide × File matrix, anchor/cleavage logos) inherits it; the filter summary now reports how many contaminant PSMs are hidden.

### Changed
- **DIAMOND BLAST tab now shows NOVEL peptides only.** De novo BLAST is only meaningful for peptides Sage's DB search missed, so `blast_filtered()` now restricts the tab to novel-peptide hits (and still respects the contaminant toggle).
- **Relabeled "Confirmed (in Sage)" → "Sage DB hits"** (summary card, results tab, info modal) to make it clear these are database-search identifications.

### Fixed
- **De novo peptide cross-references silently returned 0 ("BLAST Hits: 0", novel∩BLAST empty).** Three different mod-strippers disagreed on named modifications: Casanovo's `[Acetyl]`/`[Carbamidomethyl]` were mangled differently by (a) R's `seq_stripped` (`gsub("\\[|\\]","")` left the mod *name* text), (b) the Sage `db_stripped` path, and (c) the DIAMOND sbatch awk (`gsub(/[^A-Z]/,"")` leaked the mod's leading capital). Introduced one canonical `build_dda_canonical_peptide()` in `R/helpers_dda.R` (strips `[...]`, `(UniMod:N)`, `+57.02`, keeps only A-Z) used by `seq_stripped`, the Sage `db_stripped` path, and matched by the regenerated sbatch awk. Per CLAUDE.md Architectural Rule #3 (one definition per concept). _Existing BLAST results for **modified** peptides were written with the old mangling, so a BLAST re-run is needed to recover those; unmodified novel hits resolve immediately._

## [3.11.10] — 2026-05-29

### Fixed
- **DDA/de novo result load failed with `value.var values [intensity] are not found in 'data'`.** The Sage LFQ parser in `R/helpers_dda.R` assumed `lfq.tsv` was long-format (`proteins`/`filename`/`intensity` columns) and called `dcast(..., value.var = "intensity")`. But Sage v0.14.7's actual `lfq.tsv` is **wide + peptide-level**: id columns (`peptide`, `charge`, `proteins`, `q_value`, `score`, `spectral_angle`) followed by one raw-intensity column per run file. The parser now detects the wide format, treats the per-run columns as the intensity matrix, and rolls peptides up to protein groups by summing intensities per run (legacy long-format input still supported for back-compat). Caught loading the Ocelot Sage search (9 runs).

## [3.11.9] — 2026-05-29

### Fixed
- **DIAMOND BLAST failed on de novo peptides via DIAMOND's DNA auto-detect ("only contain DNA letters").** DIAMOND samples the query FASTA to guess protein-vs-nucleotide; because the Casanovo de novo peptide set is short and front-loaded with poly-alanine runs, DIAMOND mistook the input for DNA and aborted in ~13s with `Error: The sequences are expected to be proteins but only contain DNA letters. Use the option --ignore-warnings to proceed.` — leaving a 0-byte `blast_results_swissprot.tsv`. Distinct from the v3.11.8 mod-tag strip fix (that worked; the FASTA is clean amino acids). Added `--ignore-warnings` to **all** de novo `diamond blastp` invocations: the auto-chain SwissProt→TrEMBL pair and the single-DB variant in `R/server_dda.R`, plus the forward/reversed target-decoy FDR pair in `R/server_denovo_viz.R`. Diagnosed and re-run against the Ocelot search (316,313 unique de novo peptides). _Note: the 5 call sites are near-duplicate `diamond blastp` blocks — a future refactor should centralize them per CLAUDE.md Architectural Rule #3._

## [3.11.8] — 2026-05-29

### Added
- **Per-file filter on DDA Results panel.** A multi-select dropdown at the top of the Results page lets the user pick which source mzML/.d files contribute to every downstream view (length hist, anchor logos, cleavage motifs, peptide × file matrix). Empty selection = combined view (default, same as before). Selecting a subset narrows everything in one place. Summary line shows e.g. "Filtered — 145,632 PSMs from 3 of 9 files."
- **Combined / Per-file compare mode.** A radio next to the file filter switches between aggregating across selected files (default) and rendering one trace per file (multi-line length hist, multi-bar anchor + cleavage charts). Lets users compare across replicates or experimental conditions in the same search.
- **Peptide × File matrix sub-tab on DDA Results.** New `DT::renderDT` table where rows are unique stripped peptides, columns are source files, and cells are PSM counts. Pre-sorted by total PSM count. Has search/filter row on every column + Copy/CSV/Excel export. Lets the user spot peptides unique to one condition (high count in one file, 0 in others). Fixed leftmost 4 columns (peptide, length, n_files, total) so they don't scroll off when there are many files.
- **`dda_filtered_psms()` reactive** is the new single source of truth for the Length & Motifs renderers, so wiring a future viz against it inherits filtering + compare-mode for free.

### Fixed
- **DIAMOND BLAST awk strip rejected Casanovo mod tags.** The auto-chain DIAMOND sbatch in `R/server_dda.R` (and a one-off variant I submitted manually for the ocelot run) only stripped `[+\-0-9.()]` from the Casanovo peptide column, leaving `[Oxidation]` / `[Deamidation]` brackets intact in the FASTA query. DIAMOND then errored at line 10 with `Invalid character ([) in sequence` and exited in 14 seconds. Replaced the strip with `gsub(/[^A-Z]/, "", $2)` — keep only uppercase amino acid letters. Same fix applied to the v3.11.7 standalone DIAMOND fallback. Verified against ocelot Casanovo mztabs (~441k PSMs).

## [3.11.7] — 2026-05-28

### Added
- **Length & Motifs sub-tab on the DDA Results panel** — first viz lens for the new mode-aware framing.
  - **Universal**: peptide length distribution bar chart (always shown). Bars colored by mode signature: HLA-I highlights 8–11; HLA-II highlights 13–18; peptidomics highlights 5–25; standard tryptic shows all neutral.
  - **HLA-only** (only renders for `hla_class_i` / `hla_class_ii` modes): P2 + PΩ anchor residue frequency bar chart, grouped by position. Helps fingerprint donor HLA alleles (e.g. A*02:01 → L at P2 + L/V at PΩ).
  - **Peptidomics-only** (only renders for `peptidomics` mode): N-term + C-term flanking residue frequencies, grouped. Diagnoses endogenous protease activity (e.g. high C-term K/R = trypsin contamination).
- **`current_dda_mode()` reactive** finds the loaded search's mode by matching `values$dda_output_dir` against `j$output_dir` in the queue rows. Defaults to "standard" if no match. Drives the mode-specific viz visibility flags `output$dda_mode_is_hla` and `output$dda_mode_is_peptidomics`.
- All three plot renderers use `plotly::renderPlotly` per CLAUDE.md's bslib safety rule (renderPlot crashes with `invalid quartz()` inside `navset_card_tab` on macOS). SVG export wired via `toImageButtonOptions`.

### Fixed
- **Casanovo v5 `--force_overwrite` doesn't survive glob conflicts.** When a stale partial mztab existed in the output dir from a prior crashed run, Casanovo v5 refused to overwrite it (despite the flag) and exited with `FileExistsError`. Patched `generate_casanovo_sbatch()` to `rm -f "$MZTAB_DIR/${BASENAME}_sequence."*` at task start, making each array task idempotent regardless of leftover artifacts from earlier crashes. Caught when ocelot search tasks 1+2 failed after the Blackwell-GPU run left 3.7-KB stub mztabs in the dir.

## [3.11.6] — 2026-05-28

### Added
- **Mode-aware DDA Export ZIP for HuggingFace Spaces viewing.** New `generate_dda_export_zip()` helper in `R/helpers_dda.R` bundles the current loaded DDA search into a flat .zip designed to render natively in HF Spaces (markdown + CSV + TSV + mztab — all formats HF previews inline). Uses the `safe_section()` pattern from CLAUDE.md architectural rule #4, so any missing artifact gets logged into MANIFEST.txt with a `[SKIPPED]` line instead of silently disappearing.
  - **Common to every mode**: `methods.md`, `settings.json` (the exact Sage v0.14.7 config that ran), `sage_results.tsv`, `peptide_length_distribution.csv`, `casanovo/*.mztab` (if Casanovo ran), `diamond_hits.tsv` (if BLAST ran), `MANIFEST.txt`, `PROMPT.md`.
  - **Mode-specific summaries**: HLA modes get `hla_anchor_residues.csv` (P2 + PΩ frequency table — the signature HLA fingerprint). Peptidomics mode gets `peptidomics_cleavage_residues.csv` (N-/C-terminal flank percentages — protease motif fingerprint).
  - **PROMPT.md** is mode-aware: includes interpretation hints specific to the analysis ("a clean HLA-I prep shows a sharp 9-mer peak"; "high C-term K/R = trypsin contamination").
- **`downloadButton("dda_export_zip")` on the De Novo / DDA Results panel**, visible once results have loaded (`output.denovo_has_data` gates it). Filename includes mode + timestamp so multiple exports don't collide.
- **Plot SVG/PNG bundling deferred** to v3.11.7+ once the mode-specific viz pages land — current bundle is data-only, but data-only is enough for HF viewing since CSVs render in HF's dataset viewer.

### Changed
- **Navbar "De Novo" dropdown renamed to "DDA workflows"** with three workflow-specific landing panels (De Novo Search, Peptidomics, HLA / MHC) above the existing Results panel. Each landing panel describes the workflow's enzyme/length/mods/charge profile, lists the downstream visualizations the user will get, and has a single "Configure → search" button that flips to the Run Search tab with `acquisition_mode = "dda"` and the correct `dda_preset` pre-selected. HLA panel also has a Class I / Class II radio so the right length window gets locked in. This turns the abstract "mode" dropdown into a discoverable entry point — a researcher who isn't sure where to start can browse the workflows and pick the one that matches their experiment.

## [3.11.5] — 2026-05-28

### Added
- **Peptidomics and HLA / MHC class I + II Sage presets.** Added 3 new entries to the DDA "Analysis mode" dropdown alongside the existing Standard / Phospho / TMT options. Each preset sets the right enzyme rule (nonspecific = `cleave_at = ""`), peptide length window, mass range, and variable mod set per published immunopeptidomics / peptidomics practice:
  - **Peptidomics**: nonspecific, 5–25 AA, 400–5000 Da. Mods: oxidation, pyro-Glu (N-term Q/E), C-term amidation (−0.984), N-term acetylation.
  - **HLA class I**: nonspecific, 8–12 AA, 700–1500 Da. Mods: oxidation + deamidation (N/Q). Numbers from BOWIE's preset table (Michael Krawitzky's Ziggy).
  - **HLA class II**: nonspecific, 13–25 AA, 1300–3000 Da. Mods: oxidation + deamidation (N/Q).
- **`generate_sage_config()` refactored to a preset table.** Each preset is now a single list-entry with explicit `cleave_at`, `min_len`, `max_len`, `peptide_min/max_mass`, `variable_mods`, `add_tmt_static`. Adding a new analysis mode is a 1-row change to the table — no scattered `if (preset == ...)` branches. Backward-compatible: standard/phospho/tmt produce the same JSON as before.
- **Mode hint under the preset selector**: live text describing what each mode does (enzyme, length, mass window, mods) so picking "HLA Class I" isn't magical.
- **`mode` field added to the DDA queue entry**: stamped at submit time from the preset selector, persisted to `~/.delimp_dda_queue.rds`. Will drive the future mode-specific results-page routing (peptidomics-results / HLA-results pages, to be wired next).

## [3.11.4] — 2026-05-28

### Fixed
- **DDA queue per-engine columns stayed "running" forever.** Brett noticed the Sage column never flipped to "done" after his ocelot search completed even though the pipeline status moved on to loading Casanovo. Root cause: the Sage and Casanovo polling observers in `R/server_dda.R` only updated the pipeline-level `values$dda_status` / `values$dda_casanovo_status`, never writing back to the queue entry's per-engine `sage_status` / `casanovo_status` fields that the `engine_badge()` cell reads. Added two helpers (`update_queue_sage_status`, `update_queue_casanovo_status`) that find the matching row by Sage job_id and write `"done"` / `"partial"` / `"failed"` into the per-engine column at the COMPLETED/FAILED detection point. Also covered the single-job-not-array Casanovo failure branch.

## [3.11.3] — 2026-05-28

### Fixed
- **Sage sbatch hit `syntax error near unexpected token '||'`** because the line-continuation backslashes after the apptainer-msconvert invocation got stripped during the v3.11.2 edit — bash then read each line as its own command, and the standalone `|| { ... }` errored at parse time, killing the job in 1 second. Collapsed the msconvert call to a single line and replaced the `\` continuations with an `if !  ... ; then  ...  fi` block. `bash -n` on the rendered sbatch now passes.
- **DDA queue showed `DDA_<timestamp>` instead of the experiment name you typed.** The queue-row builder in `server_dda.R` read `input$dda_analysis_name`, but the UI control is actually `input$dda_experiment_name` — the names drifted apart at some point and the queue always silently fell back to the auto timestamp. The SLURM job names and on-disk output dir were always correct because they were already wired through `exp_name`; only the queue display was wrong. Now uses `exp_name` directly. Renamed the Casanovo-prep SBATCH job from `delimp_mgf_*` to `delimp_casanovo_prep_*` (and the log filenames) since it no longer does MGF conversion for .raw inputs — Sage's mzML output is reused.
- **"Refresh GPU queue" link took ~45s or appeared to do nothing** — two bugs in `probe_gpu_queue()`: (1) it was reading `values$ssh_config` (wrong name — the DDA module uses the `dda_ssh_config()` reactive that builds the list from `input$ssh_host/user/key_path/port/modules`), so the probe never had credentials and silently bailed; (2) it didn't pass `sbatch_path` to `check_casanovo_gpu_queue()`, forcing each of the 3 SLURM calls into login-shell mode (~15s each on Hive). Now uses `dda_ssh_config()` + `values$ssh_sbatch_path` (full path to sbatch → skips login shell). Also shows a "Probing gpu-a100 queue..." placeholder while the probe is in flight so the user can see the click registered.

## [3.11.2] — 2026-05-28

### Fixed
- **Sage on Thermo `.raw` data now actually works.** The DDA submit flow was passing `.raw` paths straight into Sage v0.14.7's `mzml_paths`. Sage silently completes with 0 PSMs because it can only parse mzML (gotcha documented in CLAUDE.md but the code path didn't enforce it). Brett's first ocelot test "completed" in 17s with 0 hits because of this. `generate_sage_sbatch()` now auto-detects `.raw` files in the input directory, calls msconvert via the pwiz Apptainer container (with the required `--bind /quobyte:/quobyte`) to convert each to mzML under `$OUTPUT_DIR/mzml/`, and passes the converted files to Sage on the command line (which overrides `mzml_paths` per Sage v0.14.7 docs). Skip-if-present so reruns don't reconvert. .d files still go straight to Sage as before.
- **DDA Casanovo MGF conversion no longer silently exits 2.** The conversion sbatch used `EXISTING_MGF=$(ls -1 "$RAW_DIR"/*.mgf 2>/dev/null | wc -l)` under `set -euo pipefail` — when no .mgf files exist, `ls` exits 2 and `pipefail` propagates through `wc -l`, killing the script before it gets to .d/.raw conversion. Replaced both `ls -1 *.mgf` sites and `ls -1d *.d` and `ls -1 *.raw` with `find -maxdepth 1` (always exits 0), and added `shopt -s nullglob` for the `for RAW_FILE in *.raw` loop. Also added the missing `--bind /quobyte:/quobyte` to the MGF-convert msconvert invocation (same Apptainer-doesn't-auto-bind gotcha).
- **Experiment name no longer silently defaults to "dda_search".** The textInput had a hardcoded default of `"dda_search"` that any new submission would land in, mixing into prior runs in the same folder. Brett's first ocelot test landed in an April 2 working directory that had stale `mzml/`, `denovo/`, and `results.sage.profile.tsv` from previous experiments. UI now starts blank with a placeholder; server rejects submit if no name given.

### Changed
- **Casanovo reuses Sage's mzML output instead of double-converting.** Previously the DDA pipeline ran msconvert twice on every .raw file — once for Sage (.raw → mzML, now baked into Sage sbatch) and once for Casanovo (.raw → MGF). Casanovo v5 reads mzML directly, so the Casanovo file-list-build step now just `find`s `$OUTPUT_DIR/mzml/*.mzML` (produced by Sage) and concatenates them with any `.d → MGF` outputs into a unified Casanovo input list. The launcher now passes `--dependency=afterok:SAGE_ID` to the convert step so Casanovo waits for Sage's mzML conversion to finish. Net effect: half the msconvert work, half the disk usage, simpler dependency chain. The Casanovo array task now strips `.mgf` / `.mzML` / `.mzml` extensions to derive output basenames.

## [3.11.1] — 2026-05-28

### Added
- **Casanovo v5 support in DDA pipeline.** New routing in `R/server_dda.R` reads `input$dda_casanovo_model` and auto-resolves the conda env + checkpoint + CLI version:
  - `casanovo_v5_*` → `/quobyte/proteomics-grp/conda_envs/casanovo5` (Casanovo 5.0.0, Python 3.13, depthcharge-ms 0.4.8 with `depthcharge.tokenizers`).
  - `casanovo_v4_*` → `/quobyte/proteomics-grp/conda_envs/cassonovo_env` (Casanovo 4.3.0, Python 3.10) — unchanged behavior.
  - The v4 env CANNOT load v5 ckpts (fails at unpickle with `ModuleNotFoundError: depthcharge.tokenizers`), so the two must be paired explicitly.
- **`generate_casanovo_sbatch(casanovo_version = "v4" | "v5")` branches the CLI** since v5 broke API compatibility: v4 uses `--output FULL.mztab`, v5 splits it into `--output_dir DIR --output_root NAME` and adds `--force_overwrite`. Verified end-to-end against PXD048461 mammoth Orbitrap DDA data (smoke job 14541990 on Hive `casanovo5` env, v5 ckpt loads cleanly + runs predictions, 4,728 spectra sequenced in 6 min on 8 CPU cores).
- **GPU vs CPU compute toggle for Casanovo** with live queue probe. New `check_casanovo_gpu_queue()` helper in `R/helpers_dda.R` issues 3 cheap SLURM calls (squeue PD + squeue R + sinfo with Gres/GresUsed) and returns pending/running counts plus per-node GPU allocation, counting only GPUs on nodes that actually accept new jobs (drain/down/maint nodes excluded from "free" count). The DDA Casanovo accordion gained a `radioButtons("dda_casanovo_compute")` (GPU / CPU) that auto-flips to whichever the probe recommends — heuristic: ≤3 pending + ≥1 free usable GPU → GPU; ≥10 pending or 0 usable free → CPU. Compute mode threads through `generate_casanovo_sbatch(compute_mode = ...)` which swaps the SBATCH header to `high` partition, drops `--gres=gpu:1`, bumps CPUs to 16 and walltime to 8h, and exports `CUDA_VISIBLE_DEVICES=""` so PyTorch stays on CPU. Probe fires when the user ticks "Run Casanovo" on and on a manual Refresh link.
- **Casanovo model dropdown** (replaces freeform textInput) with explicit choices `casanovo_v5_0_0` (default, recommended) and `casanovo_v4_2_0` — eliminates typos that would route to a nonexistent ckpt.

### Fixed
- **Sage fragment tolerance unit + value.** The DDA UI control was labeled "Fragment tol. (Da)" with default 0.05, and `generate_sage_config()` was applying a fake `* 1000` conversion to fake-build a ppm range. Sage v0.14.7 takes ppm directly (per `DOCS.md`: `fragment_tol: Dictionary with either ppm or da as keys`). Renamed the parameter `fragment_tol_da` → `fragment_tol_ppm`, set the UI default to 20 ppm (typical for Orbitrap Exploris-class DDA), removed the bogus multiplier, and updated search_info.md text + queue entry defaults to match. Anyone who ran a DDA search since the helper was first written was searching ±50 ppm instead of the intended fragment window.

### Documentation
- Added Casanovo v5 env gotcha + CLI delta to `CLAUDE.md`.
- Split Casanovo v4/v5 env entries in `docs/HPC_PATHS.md` with explicit ckpt pairing.

## [3.11.0] — 2026-05-26 / 27

Reference registry expanded (2026-05-27):
- **Pig** (`susScr_Sscrofa11.1`) — Ensembl 110
- **Rat** (`rn7_mRatBN7.2`) — Ensembl 110
- **Arabidopsis** (`arabidopsis_TAIR10`) — Ensembl Plants 58
- Bovine + maize building; will be added when complete

Reference Genome dropdown on the Build Database page now offers 5 species (human, mouse + the 3 above). New reference build script (`references/scripts/build_reference_genome.sh`) downloads genome + GTF + ncRNA from Ensembl, extracts rRNA by biotype, builds bowtie2 rRNA index + STAR genome index, writes a per-species pending JSON entry. A second script (`references/scripts/merge_registry_pending.sh`) safely merges pending entries into `registry.json` with jq + backups, moves applied entries into `registry_pending/applied/` for a paper trail. Idempotent and re-runnable.

Proteogenomics workflow is now feature-complete and integrated end-to-end into the main DE-LIMP FASTA picker.

### Added
- **Proteogenomics DBs as a first-class FASTA source on the main search page.** New entry in the `FASTA Database` sidebar dropdown (between `Database Library` and `Pre-staged on server`) opens a dedicated modal listing only `content_type == "proteogenomics"` catalog entries with proteog-specific columns (Project, Organism, Samples, Reference, Sequences, Built). Detail panel surfaces the full self-describing build metadata: project dir, pipeline ID (`proteogenomics_v1.1`), reference key, read-length tier, sample names, UniProt input, FASTA path on Hive, sequence count, file size, and the methods paragraph. "Use This Database" wires the selection into the DIA-NN search the same way Database Library does.
- **Auto-`assemble` SLURM stage chained to every new proteogenomics build.** `submit_proteogenomics_build()` now generates an `assemble.sbatch` after `rewrite` that concatenates the predicted ORFs with an optional UniProt FASTA, attempts `seqkit rmdup -s` if available (degrades to plain `cat` if not), and writes the final FASTA to `/quobyte/proteomics-grp/de-limp/databases/proteogenomics/<project>_proteogenomics_<YYYY_MM>.fasta`. The orchestrator no longer leaves an `"unknown"` stage 11 behind.
- **`submit_assemble_only()` + per-row "Assemble" button** in the Active Builds table — covers legacy builds whose SLURM chain finished pre-auto-assemble. Click opens a modal with the same UniProt source dropdown; submit generates and submits only the assemble.sbatch, patches the stage's `job_id` + `status` in `status.json`, and the poller picks it up.
- **UniProt + NCBI download integration on the proteogenomics page.** Step 4 has a new "UniProt FASTA" dropdown with options `None / Download from UniProt / Download from NCBI / Enter path on Hive`. Reuses `search_uniprot_proteomes()` / `download_uniprot_fasta()` and `ncbi_search_assemblies()` / `ncbi_download_proteome()` via proteog-prefixed observers so they don't clobber the main search page's `values$fasta_info`. Downloads land in `/quobyte/proteomics-grp/de-limp/databases/uniprot/` on Hive; subsequent picks detect the cached file and skip the redownload. NCBI downloads also upload the side-car `_gene_map.tsv` alongside the FASTA.
- **Auto-submit assemble after a UniProt/NCBI download in the per-row Assemble flow.** When the user picks "Download from UniProt" from the Assemble modal, the assemble job fires automatically as soon as the upload to Hive completes — no second click needed.
- **FASTA library auto-registration on assemble completion.** `poll_proteog_build_status()` detects the transition to `current_stage == "complete"` and adds a new entry to `~/.delimp_fasta_library/catalog.rds` matching the existing schema plus proteogenomics extension fields (`proteog_pipeline_id`, `proteog_project_dir`, `proteog_methods_paragraph`, `proteog_sample_names`, `proteog_reference_key`, etc.). The `library_entry_id` is stored back in `status.json` so polling doesn't re-register.
- **"Discover from Hive" button** in the Proteogenomics Databases modal — scans `/quobyte/proteomics-grp/de-limp/rnaseq/*/status.json` over SSH, registers every `assemble == "complete"` build that isn't yet in the user's local catalog. Solves the multi-user case: catalog stays per-user (no shared-write conflicts) but every lab member can populate it from the shared FASTA storage on demand.
- **"Restore from Hive" button** on the Build Database page's Active Builds card — same scan pattern, restores in-progress build entries to `values$proteog_build_jobs` after a Shiny restart or on a fresh machine. Defensive parsing throughout.
- **Active builds list persistence** at `~/.delimp_proteog_builds.rds`. Persist observer is gated against overwriting a non-empty file with an empty list.

### Changed
- **`launch_slims_download()`, `launch_ena_download()`, and `poll_download_status()` are now SSH-aware.** Accept a `ssh_config = NULL` parameter; when non-NULL they `ssh_exec mkdir`, `scp_upload` the status JSON + shell script, and detach the background process on Hive via `nohup bash ... </dev/null &`. The `</dev/null` is mandatory — without it the SSH connection waits on the background process's stdin and never returns. Brett's workflow (DE-LIMP on Mac, SSH to Hive) now works end-to-end.
- **Deferred build submit for sra/slims source modes.** The submit handler no longer calls `submit_proteogenomics_build()` immediately for downloaded sources; instead it stores the request in `pending_build_submits` + adds a placeholder row to the Active Builds list (blue "downloading" badge). A 15-second `reactivePoll` observer watches `download_status.json` and fires the build the moment the state flips to `"complete"`. Failed downloads get a `dl-<state>` badge instead of vanishing.
- **Skip-if-present check on submit.** Before launching any download for sra/slims modes, the submit handler `find`s for `*.fastq.gz` under the target project_dir on Hive. If files are already there, it skips the download and treats the submit as local-mode.

### Fixed
- **`.empty_or_str` crashed on `status.json` fields serialized as `{}`.** `jsonlite::fromJSON()` parses empty JSON objects as `list()` (length 0); `nzchar(as.character(list()))` returns `logical(0)`, crashing the `if`. The helper now coerces every degenerate shape (`NULL`, `list()`, `character(0)`, scalar `NA`, multi-element vectors with NAs) to `""` up front. Hoisted from function-local to file scope.
- **Stages with no job_id were wrongly summarized as "complete".** Pre-auto-assemble status.json had stage 11 (`assemble`) with `status: "unknown"` and `job_id: {}`. The poll loop called `.sacct_state("")` → garbage → `any_running` stayed FALSE → `current_stage = "complete"`. Stages with no SLURM job_id are now treated as still-pending.
- **Restore/Discover handlers no longer crash on edge cases in `status.json`.** Every `if`, `%in%`, `nzchar()`, `is.na()` is now wrapped to handle the full set of degenerate shapes JSON can return. Failures fall back to a red toast.
- **Reactive value access outside a reactive consumer crashed every new Shiny session.** The initial proteog-builds restore-from-disk block read `values$proteog_build_jobs` at server function entry; Shiny rejects that with "Can't access reactive value outside of reactive consumer", killing the entire `server_proteog_builder` module for that session. Wrapped in `isolate({})`.

### Internal
- `generate_assemble_sbatch()` in `helpers_rnaseq.R` — small SLURM job that does `cat predicted_orfs.fasta [uniprot.fasta] > out` plus optional `seqkit rmdup -s`.
- File-scope `.empty_or_str()` and `is_scalar_char_safe()` helpers in `server_proteog_builder.R` for shared NA/NULL/list-coercion across orchestrator, poller, and library-register paths.
- New reactiveVals inside `server_proteog_builder()`: `pending_build_submits`, `proteog_assemble_target`, `proteog_uniprot_state`, `proteog_ncbi_state`.
- New protected nav tab value: `build_database_tab`.

## [3.10.33] — 2026-05-22

### Fixed
- **GSEA failed with "Bioconductor version cannot be validated; no internet connection" when the organism's annotation package wasn't already installed.** Organism detection worked (e.g. Bos taurus → `org.Bt.eg.db`), but the install step called `BiocManager::install()`, whose online version-validation check throws that error even when the machine has internet (the UniProt organism lookup right before it had just succeeded). On R 4.6 / Bioc 3.23 this hit any non-preinstalled organism (cow, dog, chicken, pig, zebrafish, …). Fix: install `org.*.eg.db` annotation packages directly from the Bioconductor annotation repo URL (`https://bioconductor.org/packages/<bioc>/data/annotation`), bypassing BiocManager's validation entirely — the same pattern `app.R` already uses for limpa. Falls back to BiocManager (validation suppressed), then to a clear, actionable error showing the exact one-line manual install command if all else fails. Bioc version resolved from BiocManager when known, else mapped from the R version (R 4.6 → 3.23, R 4.5 → 3.21).

## [3.10.32] — 2026-05-08

### Changed
- **Default WSL data directory moved from `~/.delimp/data` (hidden) to `~/DE-LIMP/` (visible).** The old default lived under our config dir, which was invisible in `ls`, file browsers, and our own shinyFiles picker. Brett's friend hit this — couldn't navigate to his own data files in the DE-LIMP file browser even though they existed. Now: app state stays in `~/.delimp/` (queue, cache, config — should stay hidden because it's not user-facing), user data goes in `~/DE-LIMP/` (visible because it's the user's stuff). The setup script's data-directory prompt now offers `~/DE-LIMP` (visible folder inside WSL) as the blank-default, with the same fallback chain (env var > config file > default). Existing installs with `~/.delimp/data_dir` config file already pointing somewhere are unaffected — they keep using whatever path they have.

## [3.10.31] — 2026-05-08

### Fixed
- **WSL deployments mis-labeled "Docker"**. The `deploy_env` detection at `app.R:312` checked `local_diann && nzchar(delimp_data_dir)` — both conditions are satisfied by both Docker mode AND the WSL setup script (which sets `DELIMP_DATA_DIR` and installs `diann-linux`). So the navbar badge + browser tab title said "DE-LIMP (Docker)" on a fresh WSL install. Added `is_wsl` detection ahead of the Docker check via `WSL_DISTRO_NAME` env, `WSL_INTEROP` env, or `Microsoft|WSL` in `/proc/version`. WSL deployments now show a purple **"WSL"** badge. Also added an explicit `is_docker` check (`/.dockerenv` file) so containers actually inside Docker are correctly labeled regardless of `DELIMP_DATA_DIR`.
- **shinyFiles file browser couldn't see directories starting with `.`** — affected anyone storing data in `~/.delimp/data` (the WSL setup script's default location). Added `hidden = TRUE` to all `shinyDirChoose` / `shinyFileChoose` calls so dotfile dirs are visible. Without it, the friend reported seeing the volume list but no `.delimp` subdirectory.

## [3.10.30] — 2026-05-07

### Fixed
- **PCA crashed with "cannot rescale a constant/zero column to unit variance" under MaxLFQ + limma pipeline.** `prcomp(t(E), scale. = TRUE)` divides each column by its standard deviation; a row of `E` (= column of `t(E)`) where all samples have the same value has SD = 0, triggering the error. `complete.cases()` filtered NAs but not zero-variance rows. Hit more often under MaxLFQ + limma where some proteins have identical intensities across all samples post-normalization. Two call sites fixed:
  1. `R/server_viz.R:1020` — DE Dashboard PCA tab (the live UI)
  2. `R/server_session.R:2022` — `figures/pca.svg` in Complete Analysis ZIP
  
  Both now also drop zero-variance rows via `apply(mat, 1, var, na.rm=TRUE) > 0` after `complete.cases`, before passing to `prcomp`.

## [3.10.29] — 2026-05-07

### Fixed
- **Docker mode also needed the .NET SDK fix.** `build_diann_docker.sh` had `FROM mcr.microsoft.com/dotnet/runtime:8.0-bookworm-slim` — runtime only, same bug as v3.10.18-26's WSL setup. DIA-NN 2.x running inside a Docker-mode search would have hit the same "cannot read .raw files, please install .NET Runtime .NET SDK 8.0.407 or later" error. Switched to `mcr.microsoft.com/dotnet/sdk:8.0-bookworm-slim`. SDK image is bigger (~700 MB vs ~200 MB) but that's the cost of Thermo .raw support. Updated the comment in `Dockerfile.search` (which COPYs from diann:2.0) to reflect that the upstream image now provides the SDK, not the runtime — DE-LIMP search containers automatically inherit the fix once `diann:2.0` is rebuilt with `bash build_diann_docker.sh`.

### Action required for Docker users
- Existing `diann:2.0` Docker image: rebuild via `bash build_diann_docker.sh` to pull the SDK base image. Existing search containers will keep running on the old runtime-only image until rebuilt.

## [3.10.28] — 2026-05-07

### Fixed
- **v3.10.27 detected the missing SDK but never installed it.** Two reasons:
  1. **Tier 1 of `install_dotnet8_runtime` short-circuited on runtime presence**, not SDK presence — Brett's box had the runtime from v3.10.18-26 (`Microsoft.NETCore.App 8.0.26`), so tier 1 returned early without checking whether the SDK was actually there. Fixed: tier 1 now checks `dotnet --list-sdks | grep '^8\.'` first; only returns early if SDK is present. If SDK is missing (runtime-only state), proceeds through tiers to install.
  2. **The verify-only auto-dispatch path didn't call `install_dotnet8_runtime`** — only ran `install_dotnet_system_deps` and `verify_diann_runtime`. So even with the tier 1 fix, an existing-binary install path skipped the SDK install entirely. Added `install_dotnet8_runtime || true` to the verify path so existing installs upgrade to SDK on next launch.

  Brett's box now upgrades on launch: tier 1 sees runtime-but-no-SDK, proceeds to tier 4 (`dotnet-install.sh --channel 8.0` without `--runtime`), installs the SDK, verification passes, search should work.

## [3.10.27] — 2026-05-07

### Fixed
- **DIA-NN 2.x requires the .NET 8 SDK, not just the runtime** — the actual root cause of the community user's "No MS2 spectra: aborting" / "cannot read .raw files" bug. From a real DIA-NN log line: *"ERROR: cannot read .raw files, please download and install .NET Runtime .NET SDK 8.0.407 or later"*. We've been installing only the runtime (`dotnet-install.sh --runtime dotnet`) since v3.10.18 — `dotnet --list-runtimes` reported `Microsoft.NETCore.App 8.0.26`, verification passed all 5 checks, but DIA-NN's runtime check looks for SDK presence and our runtime-only install fails it. Two fixes:
  1. Tier 4 now installs the full SDK (drops `--runtime dotnet` from the `dotnet-install.sh` command — without that flag the script installs the SDK by default, which includes the runtime).
  2. Verification block now checks `dotnet --list-sdks | grep '^8\.'` in addition to runtimes. Fails with a specific "DIA-NN 2.x requires the SDK" message if the runtime is present but no SDK is.

  This is the actual fix for the bug class. Twelve install-path versions from v3.10.15 to land at the real cause.

## [3.10.26] — 2026-05-07

### Fixed
- **Smoke test was treating `EXIT=0 + empty stdout` as failure** — but DIA-NN 2.x's `--help` apparently exits cleanly without printing anything. The "DIA-NN runtime verification FAILED" Brett saw on v3.10.25 was a false positive: the binary actually ran fine, the test was just wrong. Replaced with two cleaner checks:
  1. `ldd "${DIANN_DIR}/diann-linux" | grep 'not found'` — direct check that all dynamic libraries resolve. If anything's missing, prints the actual library names so the user knows exactly what to apt-install.
  2. `diann-linux --help` exit code only — if it exits non-zero, fail with stderr captured. Empty stdout is no longer treated as a problem.

  This finally aligns the verification with reality: Brett's binary IS healthy on v3.10.25; we were misreporting it as broken.

## [3.10.25] — 2026-05-07

### Fixed
- **v3.10.24's `.NET 8 system deps` install was scoped to tier 4 only** — meaning on every subsequent run, tier 1 short-circuited (".NET 8 already installed via dotnet-install.sh") and the deps install was skipped. Brett's box still had the missing libicu/libssl3/etc., still failed the smoke test. Hoisted to a separate `install_dotnet_system_deps()` function that runs (a) inside `install_diann()` after `install_dotnet8_runtime()` succeeds and (b) on the verify-only auto-dispatch path before calling `verify_diann_runtime`. Apt is fast (~1s) when packages are already present; running it every launch is fine.
- **Smoke test now sets runtime env (LD_LIBRARY_PATH + DOTNET_ROOT) before invoking `diann-linux --help`.** DIA-NN's bundled libs need `LD_LIBRARY_PATH=DIANN_DIR`; .NET 8 needs `DOTNET_ROOT=/usr/share/dotnet`. Both are set in `${DELIMP_BASE}/env.sh` for normal app launches but the bare smoke-test invocation didn't source any env, so the binary was crashing during library resolution before printing anything. Plus: capture exit code separately so "crashed silently with no output" is distinguishable from "exited 0 with no output", and surface the exit code + a manual-fix hint (`ldd ...| grep 'not found'`) when the binary fails so users can diagnose missing libs themselves.

## [3.10.24] — 2026-05-07

### Fixed
- **`dotnet-install.sh` doesn't install .NET 8 system dependencies, so DIA-NN's binary couldn't run.** Brett's box: tier 4 succeeded (`.NET 8.0.26 installed via dotnet-install.sh`), DIA-NN binary extracted, RawFileReader DLLs present. But `diann-linux --help` produced **no output** — meaning the binary started but couldn't load the .NET runtime libraries because system deps (libicu, libssl3, libgssapi-krb5-2, libstdc++6, libunwind8, liblttng-ust1, etc.) weren't pulled. The script's own warning `Note that the script does not resolve dependencies during installation` was the giveaway.
  - Tier 4 now also runs `apt-get install` for the standard .NET 8 runtime deps after `dotnet-install.sh` succeeds.
  - libicu version varies by Ubuntu release (libicu76 on 26.04, libicu74 on 24.04, libicu72 on 23.x, etc.) — tries each in order until one matches.
- **Smoke test now captures and displays stderr.** Previously `2>&1 | head -1` collapsed the output, so when the binary failed silently the user just got "produced no output" with no actionable info. New smoke test redirects stderr to a tempfile, surfaces the first 10 lines if the binary failed, and adds a manual-fix hint pointing at the apt-get command. So if there's a NEW dep missing on a future Ubuntu, the user sees the actual error instead of guessing.

## [3.10.23] — 2026-05-07

### Fixed
- **`dotnet-install.sh --version 8.0` was the wrong flag** — `--version` expects an exact version like `8.0.11`; the script tried to fetch `dotnet-runtime-8.0-linux-x64.tar.gz` (literal "8.0" in the filename), got 404 from both primary and secondary CDN URLs, and aborted with `Could not find .NET Core Runtime with version = 8.0`. For "latest 8.x release" the correct flag is `--channel 8.0`. Brett caught this when v3.10.22 finally got the install path to actually fire on his Ubuntu 26.04 box (after the v3.10.21/22 dispatch fixes) — first real attempt revealed the .NET install command itself was wrong.

## [3.10.22] — 2026-05-07

### Fixed
- **Auto-dispatch silently did nothing when license-accepted-but-binary-missing.** v3.10.21's gate was `if ! -x diann-linux && ! -f license_flag` for install, `elif -x diann-linux` for verify. Brett's box had the license flag (accepted earlier) but no binary (because the v3.10.16 .NET install had aborted before reaching the download step). Both gate conditions were false → no install, no verify, no diagnostic output. Removed the license-flag check from the gate — `install_diann()` itself handles the license prompt skip internally when the flag exists. New gate: missing binary → install_diann (which now installs .NET + binary + verifies); present binary → verify_diann_runtime independently. Brett can now see why his Ubuntu 26.04 DIA-NN install never completed.

## [3.10.21] — 2026-05-07

### Fixed
- **`sync_repo` was never called on subsequent runs in `auto` mode** — the gate `if [ ! -d "${REPO_DIR}/.git" ]; then sync_repo; fi` meant the WSL-side clone only got cloned on first install, and was never `git pull`'d again. Users who ran the launcher repeatedly stayed on whatever code was current when they FIRST installed, even with the host-side `git pull` working. Brett saw v3.10.16 stuck in the app banner across 4+ launcher runs while origin was at v3.10.20. Removed the gate — `sync_repo()` is now called every run; it handles both the clone-fresh and pull-existing cases internally.
- **`verify_diann_runtime()` was never called on subsequent runs** — was nested inside `install_diann()`, only running on the rare runs that re-installed DIA-NN (i.e. `! -x diann-linux && ! -f license_flag`). Hoisted both `install_dotnet8_runtime` and `verify_diann_runtime` to top-level functions. Verification now runs (a) at the end of `install_diann()` after the binary download, and (b) separately on every auto-mode launch when DIA-NN is already installed. Result: silent .NET / RawFileReader drift (e.g. dotnet upgraded out of 8.x, DIA-NN reinstalled without RawFileReader DLLs) gets caught immediately on the next launch instead of failing during a real search.
- **`verify_diann_runtime` was called BEFORE the DIA-NN binary was downloaded** — would always fail on first install. Moved the call to the end of `install_diann()` after the extract step.

## [3.10.20] — 2026-05-07

### Fixed
- **`sync_repo()` in `delimp_wsl_setup.sh` was running stale code without telling anyone**. The launcher's WSL-side DE-LIMP clone (at `~/.delimp/DE-LIMP/`) is independent of any Windows-side clone the user's `git pull`-ing on the host. The setup script `git pull --ff-only`'d the WSL clone but on failure (shallow-clone history diverged, force-pushed tags, local edits) it just printed `git pull failed — continuing with existing code.` and ran the stale version. Brett's box `git pull`'d the Windows clone successfully, the launcher reported "running" — but the running app showed `v3.10.16` while the Windows clone was at `v3.10.19` because the WSL clone had silently stayed pinned. Two-tier fix:
  1. If `git pull --ff-only` fails, try `git fetch --depth 1 origin main && git reset --hard origin/main` — covers the shallow-divergence case that's most common.
  2. If even that fails, print a loud `err` block telling the user how to force a clean re-clone (`rm -rf ~/.delimp/DE-LIMP/`).
- **Print the running version + short commit SHA after sync** so users can see at a glance which code is about to run. Was previously invisible — you'd only see the version after the app boots and prints its banner, which is too late if you wanted to verify the pull worked before sitting through a 30-min install.

## [3.10.19] — 2026-05-07

### Added
- **Post-install DIA-NN runtime verification in `delimp_wsl_setup.sh`** — catches the "No MS2 spectra: aborting" bug at install time instead of letting the user discover it 5 minutes into a search. After the .NET 8 + DIA-NN binary install, runs four checks:
  1. `dotnet` is on PATH
  2. `dotnet --list-runtimes` reports a `Microsoft.NETCore.App 8.x` runtime (NOT just any version)
  3. `diann-linux` binary is at the expected path and executable
  4. RawFileReader DLLs are bundled in the DIA-NN install dir (the Thermo .raw reader)
  5. `diann-linux --help` returns output without crashing (smoke test that the binary + .NET marriage actually works)

  Each check prints ✓ / ✗ to the install log. If any fails, the script aborts with a clear message: *"This is the bug class behind 'No MS2 spectra: aborting' errors during searches. Fix the issues above before submitting a search, or DIA-NN will fail to read Thermo .raw files."* — so users know exactly what's broken and why before they invest 1+ hour in a search that's going to fail.

## [3.10.18] — 2026-05-07

### Fixed
- **Robust .NET 8 runtime install in `delimp_wsl_setup.sh`** — likely root cause of the community-reported "No MS2 spectra: aborting" Thermo .raw read failure. DIA-NN's `RawFileReader` requires .NET 8; if it's missing or wrong-version, library prediction still works (pure C++) but raw file reading silently fails. Old install logic had two failure modes: (a) on Ubuntu 26.04 (released April 2026) the script added Microsoft's apt repo but Microsoft hadn't yet published `dotnet-runtime-8.0` for that release, so apt-get failed and the entire setup script aborted; (b) if `dotnet` command was already on PATH at any version (e.g. stale 6.x or broken install), the script's `command -v dotnet` short-circuited and skipped reinstall — DIA-NN then ran with a wrong-version dotnet and silently couldn't read .raw. Replaced with a 4-tier install:
  1. **Detect existing dotnet 8.x** specifically (via `dotnet --list-runtimes`), skip only if a real 8.x runtime is present
  2. **apt with multiple package-name candidates** (`dotnet-runtime-8.0`, `dotnet-runtime-8`, `dotnet8`) — covers Ubuntu naming shifts across 22.04 / 24.04 / 26.04
  3. **Microsoft's apt repo + same package-name sweep**
  4. **Last-resort: Microsoft's official `dotnet-install.sh`** — works on any Linux distro/version regardless of apt channel availability. Installs to `/usr/share/dotnet`, symlinks `/usr/local/bin/dotnet` so DIA-NN finds it.

  Verified to install cleanly on Ubuntu 26.04. Brett caught it stress-testing the v3.10.15-17 install path on a fresh Windows 26.04-WSL box.

## [3.10.17] — 2026-05-07

### Fixed
- **WSL launcher's auto-browser-open was port-blind on first install.** It used a hardcoded 90-second `timeout` then opened `http://localhost:3838` — but a fresh install spends 20–30 minutes compiling R + Bioconductor packages, so the browser hit the port long before Shiny was listening, got "can't connect", and confused the user into thinking something was broken (it wasn't — install just hadn't finished). Replaced with a port-aware PowerShell poll: every 2 seconds, check `Get-NetTCPConnection -LocalPort 3838 -State Listen`. Opens the browser only once the port is actually accepting connections. 60-minute timeout (silent exit, no false-positive open) covers even the slowest first install.
- **First-install banner in the launcher console now warns prominently** that the install takes 20-30 min and that the browser-open is port-aware (won't fire prematurely). Includes the manual-open URL (`http://localhost:3838`) so users who want to open it themselves can.

## [3.10.16] — 2026-05-07

### Fixed
- **`Launch_DE-LIMP_WSL.bat` falsely reported "Ubuntu installed" on a fresh Windows box with no WSL distro.** The Ubuntu existence probe used `wsl -d Ubuntu -e true >nul 2>&1` followed by `if errorlevel 1`. Windows batch's `if errorlevel N` evaluates as "errorlevel >= N", and WSL returns negative-ish exit codes for `WSL_E_DISTRO_NOT_FOUND` that the test misinterprets as success. So the launcher skipped the auto-install path and proceeded to copy + run, then failed downstream with `"There is no distribution with the supplied name. Error code: Wsl/Service/WSL_E_DISTRO_NOT_FOUND"`. Replaced with a sentinel-string probe: run `echo __DELIMP_UBUNTU_OK__` inside Ubuntu and `findstr` for the literal sentinel — exit-code-independent, works regardless of Windows / WSL version quirks. If the sentinel isn't in the output, trigger `wsl --install -d Ubuntu` as designed. Caught while doing a fresh Windows install of v3.10.15.

## [3.10.15] — 2026-05-07

### Changed (de-UCD-specific structural defaults)
- **Introduced `R/helpers_site.R`** with `delimp_site()` — single source of truth for site-specific structural defaults (storage prefixes, shared FASTA library paths, shared activity log path, SLURM primary/fallback partition pair, gene-map search dirs). UCD defaults preserved exactly; non-UCD sites override via env vars (`DELIMP_*`) or `~/.delimp_site.yaml`.
- **`translate_storage_path()`** now reads `storage_local` / `storage_hpc` from the site config instead of the hardcoded `/Volumes/proteomics-grp/` ↔ `/quobyte/proteomics-grp/` regex pair. UCD users see no change. Non-UCD users with `DELIMP_STORAGE_LOCAL=...` and `DELIMP_STORAGE_HPC=...` get path translation that actually works for their layout.
- **`speclib_cache_path()`, `speclib_cache_is_shared()`, `fasta_library_path()`, `fasta_library_is_shared()`** now consult `delimp_site()$shared_diann_*` and `shared_fasta_lib_*` instead of hardcoded UCD paths.
- **`activity_log_path()`** now uses `delimp_site()$shared_activity_log` (UCD default preserved).
- **Auto-queue partition picker (`select_best_partition()`)** now reads primary/fallback from `delimp_site()$slurm_account` + `slurm_partition` and `slurm_fallback_account` + `slurm_fallback_partition`. UCD users get the same `genome-center-grp/high` ↔ `publicgrp/low` behavior; non-UCD users with their own SLURM setups can configure their pair.
- **`server_de.R` / `server_gsea.R` gene-map search paths** now use `delimp_site()$gene_map_dirs` (colon-separated `DELIMP_GENE_MAP_DIRS` env var). UCD users get `c("/data/fasta", "/quobyte/proteomics-grp/de-limp/fasta")`; non-UCD users can configure.
- **`resolve_fasta_dir()`** new helper replaces three sites of `getOption("delimp.fasta_dir", default = "/quobyte/proteomics-grp/de-limp/fasta")`. Tries (in order): explicit programmatic override → `delimp_site()$fasta_dir_local` → `delimp_site()$fasta_dir_hpc` → `~/.delimp_fasta` (created on demand). Fixes the community-reported bug where non-UCD users had `/quobyte/proteomics-grp/de-limp/fasta/` silently created on their local WSL/Docker filesystem because `dir.create(..., recursive = TRUE)` succeeded against a path that looked-like-but-wasn't the UCD HPC mount.

### Backwards compatibility
- All defaults in `delimp_site()` preserve historic UCD behavior. Brett's lab and other UC Davis users see zero functional change.
- UI textbox defaults (SSH host `hive.hpc.ucdavis.edu`, SLURM account `genome-center-grp`, SIF path, raw-data placeholder) **were intentionally left as UCD strings** — they're user-overrideable in the UI, so they're not "structural" hardcoding. Non-UCD users type over them once and DE-LIMP persists their choice via session settings.

## [3.10.14] — 2026-05-06

### Added
- **Complete Analysis ZIP's PROMPT.md now includes "Appendix A: How This Analysis Works"** — an educational background section ported from the cascadia-denovo branch's de novo Claude prompt. Targets a PhD student or biologist with no mass-spec / bioinformatics background. Covers LC-MS/MS basics, DIA vs DDA, DIA-NN, the statistical framework, and key terms (logFC, P-value, FDR, volcano plot, CV, normalization) — all in plain language with analogies. Pipeline-aware: under MaxLFQ + limma it explains the Moschem 2025 protocol and quantile normalization; under DPC-Quant it explains limpa's detection-probability model. The downstream LLM uses this as a guide when writing manuscript Methods + biological interpretation, so every Complete Analysis export now ships with built-in educational scaffolding for the reader.

## [3.10.13] — 2026-05-06

### Fixed
- **`parse_search_info_md()` was extracting nothing — Methods text after Load-from-HPC had no instrument settings.** Three bugs in the parser, all fixed by re-reading and stress-testing against Brett's real `search_info.md`:
  1. **Greedy character class ate the bold markers.** The regex `^[\\s\\-*]*\\*\\*([^*]+)\\*\\*...` had `*` in the leading class — and `*` is greedy, so it consumed the `**` of `**Instrument**`, breaking the match. Every line failed; `kv` stayed empty; `instrument_metadata` was NULL. Replaced with a two-step extract (locate the bold-key token, then peel value off the rest of the line).
  2. **Whitespace stripped before underscore conversion.** "Acquisition mode" → my key normalization first stripped all whitespace (`gsub("[*:[:space:]]", "", ...)`) then tried to convert whitespace to underscores. By that point there was no whitespace left, so the key became `acquisitionmode` and didn't match `kv$acquisition_mode`. Now strips only `*:` first, THEN collapses whitespace runs to `_`.
  3. **FASTA paths weren't extracted.** They live under a `### FASTA Files (N)` heading as one bullet per line with backtick-quoted paths — not as a `**FASTA**:` key-value. Added a section parser that switches into "FASTA Files" mode on the heading and pulls each backtick-quoted path until the next heading.

  Verified against Gemma_set2's real `search_info.md`: instrument_model = "timsTOF HT", serial, acquisition_mode = "dia-PASEF", DIA windows = 37, m/z range = "100-1700", and both FASTA files all extract correctly.

## [3.10.12] — 2026-05-06

### Fixed
- **Queue still showed phase-substep names like `diann_Gemma_set2_s5_report` after the v3.10.11 collapse.** Two issues remained: (1) array-task entries (`#13828143_0`) from earlier broken Recovers were never filtered from existing queues, and (2) the v3.10.11 collapse keyed on `(base_name + output_dir)` — but substep entries often have different or empty output_dirs and end up in distinct groups, so dedup was incomplete; surviving entries also kept their substep name. The startup queue cleanup now: (a) drops any entry whose `job_id` matches `^\d+_\d+$` (array tasks); (b) groups by `base_name` only, ignoring output_dir; (c) rewrites the surviving entry's `name` to the clean base name so the queue UI shows "Gemma_set2" instead of "diann_Gemma_set2_s5_report". Console message reports `dropped N array-task entries, collapsed M substep entries -> K logical searches`.

## [3.10.11] — 2026-05-06

### Fixed
- **CRITICAL — v3.10.10 broke SSH connection entirely.** Changed `env = c("current", MallocStackLogging = "")` to `MallocStackLogging = NA_character_` thinking it would actually unset the var. But `processx::run()` rejects NA values in env: `is.null(env) || is_env_vector(env) is not TRUE`. Reverted to empty string. SSH works again. The MallocStackLogging console noise stays (it's harmless macOS chatter from RStudio's parent process forking children — DE-LIMP can't suppress it from R-land without breaking processx).
- **Queue + History tabs were full of useless parallel-pipeline substep entries.** v3.10.10's "collapse phase substeps into one logical search" used a lazy regex `^diann_(.+?)(_s[1-5]_[a-z]+)?$` — but R's default `sub()` is POSIX ERE which doesn't support lazy quantifiers (`.+?` is interpreted as `.+` followed by literal `?`). The regex silently failed to match anything, so dedup was a no-op. Replaced with two simple non-lazy `sub()` calls: `sub("^diann_", "")` then `sub("_s[1-5]_[a-z]+$", "")`. Verified to actually match all five phase suffixes (`_s1_libpred`, `_s2_firstpass`, `_s3_assembly`, `_s4_finalpass`, `_s5_report`).
- **Existing queue duplicates now get cleaned up at startup.** Users who ran v3.10.10 (or earlier) accumulated phase-substep entries in `~/.delimp_job_queue.rds`. Added a one-time collapse pass in the queue-load observer that runs the same dedup-by-(base_name + output_dir) algorithm. Console message reports how many entries were merged.

## [3.10.10] — 2026-05-06

### Fixed
- **"Recover" button now scopes to the SSH-connected user.** `recover_slurm_jobs()` previously ran `sacct` with no `-u` flag — depending on cluster policy, this could return other lab members' jobs. Now passes `-u <ssh_user>` from `cfg$user`. Per Brett's directive: "only recover jobs submitted by the person in the username box that tested the connection with HIVE."
- **Recover no longer creates 23 queue entries per parallel-pipeline search** (5 phase substeps + 20 array-task rows). Two filters added in `recover_slurm_jobs()`:
  1. `grep -v '^[0-9]\+_'` drops array task rows (e.g. `13828143_0`) — they're substeps of the parent.
  2. The Recover handler now collapses the remaining `diann_<NAME>_s<N>_<phase>` substep rows into one logical search per unique base name (preferring `_s5_report` as canonical for output_dir).
  Net: `Gemma_set2`'s 23 sacct rows → 1 queue entry.
- **Recover no longer freezes the queue render with a 490-row sacct dump.** The handler used to `values$diann_jobs <- c(values$diann_jobs, list(job_entry))` *inside* a tight loop. With N rows: O(N²) memory churn, N reactive invalidations, N persistence writes. Brett saw "490 jobs recovered" but the queue stayed empty (render couldn't keep up). Refactored to accumulate new entries in a local list, then assign once at the end → one invalidation, one persistence write, one render pass.
- **Dedup by `output_dir` as well as `job_id`** — re-running Recover doesn't pile up duplicate entries for the same logical search even when the original entry's `job_id` was the array parent and the recovered row has the report-step ID.
- **Recover and Load-from-HPC now read `search_info.md`** so LC + mass spec settings flow through. Added `parse_search_info_md()` helper in `R/helpers_search.R` — parses the markdown key-value lines into the same shape that submit-time code populates into `values$diann_search_settings`, plus an `instrument_metadata` block. Both the Recover handler and the Load-from-HPC handler now SCP `search_info.md` from the search's `output_dir` and merge the parsed settings. Brett's question: "Load from HPC isn't supposed to read the search LC and mass spec settings?" — yes it is, now it does.
- **macOS `MallocStackLogging` console noise suppressed.** `processx::run(env = c("current", MallocStackLogging = ""))` in three places set the var to empty-string, which macOS still complains about ("can't turn off malloc stack logging because it was not enabled"). Switched to `MallocStackLogging = NA_character_` which actually unsets the var — the warnings stop.

## [3.10.9] — 2026-05-06

### Added
- **Complete Analysis ZIP now includes a `figures/` subdirectory with 9 publication-quality SVG figures**: volcano, heatmap_top20, violin_top10_up, violin_top10_down, pca, qc_group_distribution, normalization_density, data_completeness, sample_correlation, pvalue_distribution. Ported from cascadia-denovo branch (commits `38c9b3b` + `c2329c8` — that work was on `R/server_ai.R` and never merged to main; v3.10.4 then consolidated all Claude exports into Complete Analysis, but missed bringing the figures along). Each figure uses the `svg() + print() + dev.off()` pattern (works on headless Linux/HF where `ggsave()`'s display path is unreliable). Each figure wrapped in `safe_section()` so failures land in MANIFEST instead of getting silently dropped.
- **PROMPT.md updated**: now references the figures/ subdirectory with a per-figure descriptor table, and instructs the LLM to use markdown image syntax (`![title](figures/X.svg)`) in its analysis. Adds analytical-context cues: reference violin plots when discussing DE proteins, sample_correlation for batch/replicate concerns, normalization_density for normalization, pvalue_distribution for study-power discussions.
- **`pvalue_distribution.svg` figure (new, not in cascadia)**: per-contrast raw P-value histogram with adj.P.Val < 0.05 callout. Spike at 0 = real signal; flat = no signal; spike at 1 = something off (model misspecification, batch effect not in design, etc.). Surfaces the diagnostic that lets a reviewer judge the experiment at a glance.
- **ZIP path-preservation fix**: `zip(file, basename(files_to_zip))` was stripping the `figures/` prefix off SVG paths, dumping them at the zip root. Now uses `normalizePath()` + relative-path computation so subdirectory structure is preserved in the archive.

## [3.10.8] — 2026-05-06

### Fixed (Complete Analysis export — silently skipped DIA-NN files)
- **`search_info.md`, `report.pg_matrix.tsv`, `report.stats.tsv` were silently skipped on Mac when the search ran on HPC.** The fetch helper checked `file.exists(file.path(output_dir_local, filename))` against the raw HPC path stored at submit time (`/quobyte/proteomics-grp/...`). On a local Mac that path doesn't exist (the share is mounted at `/Volumes/proteomics-grp/...`), so the local check returned FALSE. SSH wasn't connected at export time, so the SSH fallback also failed. The caller did `if (!is.null(f)) files_to_zip <- c(...)` — silently dropping the file with no MANIFEST entry. Two CLAUDE.md Architectural Rule violations in one section: Rule #4 (silent catch in export paths) and a missed `translate_storage_path()` call.
- **Fix #1 (path translation)**: `fetch_diann_file()` now also tries the translated local path via `translate_storage_path(output_dir, to = "local")` (e.g. `/quobyte/...` → `/Volumes/...`) before falling back to SSH. Three lookup tiers: original `output_dir` → translated local mount → SSH SCP.
- **Fix #2 (no more silent skips)**: the three fetch calls (`search_info.md`, `report.pg_matrix.tsv`, `report.stats.tsv`) are now wrapped in `safe_section()`. Failures get a `[SKIPPED] <filename> -- <reason>` line in MANIFEST.txt with the actual paths that were tried. Inlined (no `for`-loop or `local()` wrapper) because either would create a new env layer and break `files_to_zip <-` propagation — same trap as the v3.10.5 `<<-` issue.

## [3.10.7] — 2026-05-06

### Fixed (Complete Analysis export — three issues from real downstream review)
- **`parameters.txt` listed stale covariates from a previous analysis.** The export read `values$cov1_name` / `values$cov2_name` unconditionally, but those reactiveValues persist across pipeline runs even when the user un-checked the "include in model" boxes for the current analysis. Brett's Gemma export listed `Covariate 1: Student` and `Covariate 2: RunOrder` from an earlier session that wasn't using the current sample groups. Now we only emit a covariate line if the covariate name actually appears in `colnames(values$design)` — the canonical "what was used" source. CLAUDE.md Architectural Rule #2 (tagged %||% defaults) and Rule #3 (single source of truth for covariate name).
- **`detection_matrix.csv` skipped with `unimplemented type 'list' in 'EncodeElement'`.** The per-protein detection-count construction used `do.call(rbind, lapply(..., c(... as.list(detected))))`, which produced a row matrix where every cell was a length-1 list. `as.data.frame()` preserved the list-columns and `write.csv()` rejected them. Rewrote the section to build a proper numeric matrix via `vapply(..., numeric(ncol(raw_mat)))`, attach `Detected_*` integer columns one at a time, and write a clean atomic data.frame. CLAUDE.md gotcha: "nested lists in data.frame() = silent breakage."
- **`PROMPT.md` confidently stated organism = Human even when detection had silently fallen back to the default.** `detect_organism_db()` matches UniProt SwissProt suffixes (`_HUMAN`/`_MOUSE`/etc.); when the IDs are NCBI RefSeq (XP_/NP_/WP_) or any non-suffixed format, it falls through to `org.Hs.eg.db`. The PROMPT then asserted "Organism: Human" with no caveat — Brett's Peromyscus californicus (mouse-mapped) export read as a human dataset and the downstream LLM almost did GSEA against `org.Hs.eg.db`. Per CLAUDE.md Architectural Rule #2, fallback values that hit user-facing text must be tagged. The PROMPT.md organism line now: (a) prints the FASTA filename inline as a hint, (b) checks whether any UniProt organism suffix is present in the IDs, (c) when no suffix found, replaces the confident "Organism: X" line with a `**WARNING — auto-detection FELL BACK to human default**` block telling the downstream LLM to verify from the FASTA name before any GSEA/KEGG/GO work, with a stronger NCBI-specific message when the IDs are XP_/NP_/WP_.

## [3.10.6] — 2026-05-06

### Fixed
- **HOTFIX: "Recover" button never recovered any DIA-NN jobs**. `recover_slurm_jobs()` in `R/helpers_search.R` queries `sacct` for diann-named jobs over the last 14 days, then ran `grep -v '\\.'` to filter out `.batch` / `.extern` substep rows. But the sacct format includes the StdOut path (`/path/.../diann_*_%j.out`), and **every** line has a dot somewhere — so the filter dropped all 23+ matching rows, leaving zero jobs to recover. Verified against the live HIVE cluster: with the broken filter `wc -l` returned 0; the parent diann jobs were obviously present. Replaced with `grep -v '^[^|]*[.][^|]*|'` which only checks for a dot in the JobID field (field 1, before the first `|`) — so `13828146.batch` is filtered, but `13828146` plus its `/path/diann_*.out` StdOut field passes through. Caught by Brett — clicked Recover after the Gemma_set2 search completed and nothing appeared in the queue.

## [3.10.5] — 2026-05-06

### Fixed
- **HOTFIX: Export Complete Analysis ZIP was missing 9 files that MANIFEST.txt claimed were `[OK]`.** The new sections added in v3.10.4 (DE_Results_Full.csv, QC_Metrics.csv, Phospho_DE_Results.csv, group_assignments.csv, parameters.txt, PROMPT.md, detection_matrix.csv, quartile_profiles.csv, variable_proteins.csv) used `files_to_zip <<- c(files_to_zip, file)` inside their `safe_section({...})` bodies. CLAUDE.md gotcha: `<<-` inside `withProgress()` falls through to the global env (the parent chain is broken because `withProgress` uses `eval(substitute(expr), env)`). Each `<<-` was silently writing a phantom global named `files_to_zip` instead of mutating the inner accumulator. The bodies didn't error so `safe_section()` recorded `[OK]`, but the actual `zip()` call only saw the OG sections that used plain `<-`. Replaced all `<<-` with `<-` in those sections — promise semantics evaluate the safe_section body in the caller's env (the `withProgress` inner eval env where `files_to_zip` lives) and `<-` writes there directly. Caught by Brett's downstream Claude analysis after a v3.10.4 export was uploaded.

## [3.10.4] — 2026-05-06

### Changed
- **Export Complete Analysis is now a true superset of "Export for Claude"**: previously the Output > Export Complete Analysis ZIP advertised PROMPT.md + detection_matrix.csv + quartile_profiles.csv + variable_proteins.csv but the handler never wrote them — the description in `R/ui.R` was aspirational. The handler in `R/server_session.R` now writes:
  - All previous files (expression_matrix, sample_metadata, methods.txt, reproducibility_log.R + sessionInfo, search_info.md, report.pg_matrix.tsv, report.stats.tsv, protein_confidence.csv, data_quality_summary.csv, contaminant_summary.csv, session.rds)
  - **Plus**: detection_matrix.csv, quartile_profiles.csv, variable_proteins.csv, group_assignments.csv, parameters.txt, **PROMPT.md** (DE-aware — adapts wording for DE vs exploratory), **MANIFEST.txt** (per-section [OK]/[SKIPPED] log), **DE_Results_Full.csv** (when `values$fit` exists), **QC_Metrics.csv** (when QC stats exist), **Phospho_DE_Results.csv** (when phospho ran).
  - Pipeline-aware throughout: PROMPT.md and parameters.txt use `pipeline_label(values$y_protein)` and `is_maxlfq()` instead of hardcoded "DPC-Quant" strings (CLAUDE.md Architectural Rule #1).
  - Every section uses `safe_section()` from `R/helpers.R` so a single failure no longer silently drops files (Architectural Rule #4) — MANIFEST.txt records what was included vs skipped and why.
- **Three redundant "Export for Claude" buttons hidden** (kept handlers behind them — no orphan-removal): Data Explorer header, AI Summary tab, AI Chat tab. Output > Export Complete Analysis is now the single export entry point. Anyone landing on the AI tab who expected an LLM bundle will instead see the consolidated export linked from Output.

### Fixed
- **FASTA browse / SSH-scan was silently selecting every `.fasta` in the directory.** In shared dirs like `/quobyte/proteomics-grp/de-limp/fasta` (which holds many species-specific FASTAs side-by-side), hitting "Scan" combined all of them into one DIA-NN search — almost always wrong. Now: 1 file → use directly; ≥2 files → modal with `checkboxGroupInput` so the user explicitly picks one (or several to combine intentionally). Same fix applied to both `fasta_browse_dir` (local shinyFiles) and `ssh_scan_fasta_btn` (remote SSH) handlers in `R/server_search.R`.

### Renamed
- **On/Off Proteins panel column rename** for clarity: `n_in_group1` / `n_in_group2` → **`detected_g1` / `detected_g2`** (the count of samples in each group where the protein was detected). `total_in_group1` / `total_in_group2` → **`total_g1` / `total_g2`** (group sizes — these are properties of group assignment, not of any individual protein, so they're constant across all rows of a given contrast). Header comment in `compute_onoff_proteins()` and the column-order vector in the DT renderer (`R/server_de.R`) updated to match.

## [3.10.3] — 2026-05-06

### Fixed
- **On/Off Proteins panel crashed with `ncol(E) == length(group_factor) is not TRUE`** when the metadata table had extra rows from excluded-files tracking (which renders excluded files as additional display rows). The reactive read `values$metadata$Group` directly, picking up the excluded rows that aren't in `values$y_protein$E`. Fixed by aligning metadata to matrix columns via `match(colnames(E), values$metadata$File.Name)` before passing to `compute_onoff_proteins()`. The on/off panel now ignores excluded files (they were already excluded from the analysis upstream) and computes against the active sample set only.

## [3.10.2] — 2026-05-06

### Fixed
- **HOTFIX: `server_viz.R` was syntactically broken in v3.10.1**: the v3.10.1 edit injected R `if/else` expressions inline into a multi-line single-quoted `paste0()` template (the no-replicates Data Exploration prompt), but didn't close the string properly — leaving a stranded `|` token at line 2467 that errored at parse time. Refactored: the dynamic table rows are now computed into named character variables (`expr_matrix_note`, `dpc_only_protconf_row`, `diann_pg_note`, `dpc_only_detmat_row`) immediately before the `paste0()` call, then inserted as plain `', var, '` interpolation. Same pipeline-aware behaviour, no syntax break. Lesson learned (added to CLAUDE.md by inference): never embed `if (cond) X else Y` inside a multi-line single-quoted R string template — pre-compute and substitute.

## [3.10.1] — 2026-05-06

### Fixed
- **Run Comparator hypothesis-engine text now pipeline-aware** (audit item #10): `assign_hypothesis()` accepts a new `run_a_id` arg (default `"dpc"` for back-compat). Tool-comparison strings inside Rule 2 (normalization offset), Rule 3 (variance/missing-value), Rule 5 (peptide rollup) now read whether Run A used MaxLFQ + limma vs DPC-Quant + limma and emit the correct contrast.
- **Settings-diff body strings** now derive Run-A peptide-usage / rollup descriptions from the descriptor (`PG.MaxLFQ (DIA-NN)` vs `DPC-Quant: empirical Bayes precursor aggregation`) and the "Critical:" note text adapts. Was previously hardcoded "DPC-Quant uses all detected precursors".
- **Violin plot Inferred → Missing under MaxLFQ** (audit item #14): label, hover text, subtitle, and the inferred-points filter all flip on `is_maxlfq()`. Under DPC-Quant the existing "Inferred (nObs=0, SE=...)" label stays.
- **Data Completeness modal + warning banner** (audit item #13) branch on pipeline. Under MaxLFQ the modal explains "Cells without precursor evidence are NA in the MaxLFQ matrix; limma drops them per row at fit time" instead of claiming DPC-Quant inferred them. Title and warning-banner wording adapt.
- **DPC Fit info modal** (audit item #12) shows a "not applicable under MaxLFQ" message when DPC-Quant didn't run, pointing the user to the Pipeline Diagnostic + filter waterfall instead.
- **Methods README in Complete Export** (audit item #15) drops the `protein_confidence.csv` and `detection_matrix.csv` table rows under MaxLFQ (they aren't written), and adapts the `expression_matrix.csv` description to mention NAs honestly under MaxLFQ vs the always-complete DPC-Quant matrix.

## [3.10.0] — 2026-05-06

### Changed (Run Comparator now pipeline-aware)
- **Run Comparator AI prompts and Claude export now describe the actual pipeline that produced Run A**, instead of always claiming "DPC-Quant". `parse_delimp_session()` now derives `rollup_method`, `de_engine`, `pipeline_id`, and `pipeline_label` from `y_protein$other$descriptor` (the canonical pipeline descriptor added in v3.9.17). Added `run_a` to the `comp_results` payload so prompt builders can read those fields.
- **`build_gemini_comparator_prompt()`** rebuilt: defines `run_a_pipeline_label`, `run_a_rollup_text`, `run_a_de_engine_text`, `run_a_missing_text`, `run_a_norm_text`, `run_a_norm_short` once at function entry from `comp_results$run_a$settings`, then substitutes them into the Spectronaut and FragPipe-Analyst comparison narratives. Under MaxLFQ + limma the prompt now correctly tells Gemini that Run A used PG.MaxLFQ + quantile normalization + plain `lmFit` (NAs in place) instead of falsely claiming DPC-Quant + DPC-CN + detection-probability modelling.
- **`build_claude_comparator_prompt()`** rebuilt: `tool_label` derives `DIA-NN/MaxLFQ/limma` vs `DIA-NN/DPC-Quant/limma` from the same descriptor; the long methodology paragraph branches on the pipeline so reviewers see the right method.

This is the v3.10.0 release because it eliminates the last large class of "wrong-pipeline-name in user-facing exports" surfaces. Some lower-priority hardcoded strings still exist in the rule-engine hypothesis text (e.g. "FragPipe-Analyst uses Perseus-style imputation; DE-LIMP uses DPC-Quant") and the settings-diff body — those will be cleaned up in v3.10.1.

## [3.9.18] — 2026-05-06

### Changed (v3.10.0 prep — sweep `values$pipeline_mode_used` reads to `is_maxlfq()`)
- **11 read sites swept** to use the canonical `is_maxlfq(values$y_protein)` accessor (added v3.9.17) instead of the volatile `values$pipeline_mode_used` reactiveVal:
  - `R/server_data.R` × 2 — meta-alignment branch + lmFit-vs-dpcDE branch
  - `R/server_session.R` × 2 — methods text branch + protein_confidence export guard
  - `R/server_ai.R` × 3 — protein_confidence + detection_matrix export guards
  - `R/server_qc.R` × 4 — Norm QC plot branch + filter waterfall visibility + stacked-bar title + stacked-bar legend label
  - `R/server_de.R` × 1 — On/Off table empty-state message
- The accessor reads from `y_protein$other$descriptor$pipeline_id` (durable) with legacy fallback to `y_protein$other$pipeline` and a final fallback to "DPC-Quant" for very old session.rds files. Adding a third pipeline (e.g. DDA) now requires editing zero of these files.
- Local boolean variables named `is_maxlfq` were removed where they shadowed the helper of the same name; call sites now invoke `is_maxlfq(values$y_protein)` inline.

`values$pipeline_mode_used` is still set during run_pipeline (back-compat for in-flight sessions) but no longer read except internally by `pipeline_descriptor()` as a fallback. Final removal in v3.10.0 after `server_comparator.R` is swept.

## [3.9.17] — 2026-05-06

### Added (v3.10.0 prep — single source of truth for pipeline metadata)
- **Pipeline descriptor objects** in `R/helpers.R` (`dpc_pipeline_descriptor()`, `maxlfq_pipeline_descriptor()`, `pipeline_descriptor(y_protein)`, `is_maxlfq(y_protein)`, `pipeline_label(y_protein)`). Each pipeline now carries its own self-describing record (id, display label, rollup method, normalization, DE engine, missing-value policy, citation) attached to `y_protein$other$descriptor` at quantification time. Downstream code can read `is_maxlfq(values$y_protein)` instead of consulting the volatile `values$pipeline_mode_used` reactiveVal in 7 different files.
- **Both pipelines now populate `$other$descriptor`**: `build_maxlfq_pipeline()` attaches `maxlfq_pipeline_descriptor()`; the post-`dpcQuant` block in `server_data.R` attaches `dpc_pipeline_descriptor()`. Legacy `y_protein` objects (from older session.rds files) gracefully fall back to the DPC descriptor.

This commit only adds the helpers and the wiring; existing `values$pipeline_mode_used` checks still work (back-compat). Subsequent v3.9.18+ commits will sweep call sites file-by-file to use the canonical accessors.

## [3.9.16] — 2026-05-06

### Fixed
- **Gene-map TSV parse errors no longer silent** (audit item #7): when an NCBI gene_map.tsv file is present but malformed, both the live Expression-Grid lookup (`R/server_viz.R:293,325`) and the Claude-export gene lookup (`R/server_viz.R:1963,1986`) used to return NULL silently, leaving the entire DE table labeled with bare accessions instead of gene symbols. Now: file-absent still falls through silently (correct behaviour), but file-present-but-unparseable raises a yellow Shiny notification naming the file and the parse error, plus logs a `[Grid]` / `[Export]` console line. The user finds out within 12 seconds instead of "weeks later when a reviewer wonders why your DE table has no gene names."
- **`req(input$contrast_selector)` accepted empty string** (audit item #11): standard `shiny::req()` considers `""` truthy, so an empty selectInput leaked through to `topTable(fit, coef = "")` which errors deep in limma. Replaced 8 sites in `server_de.R` and `server_gsea.R` with `req_nzchar(input$contrast_selector)` (helper added in v3.9.14) that explicitly rejects empty strings.

## [3.9.15] — 2026-05-05

### Added
- **`MANIFEST.txt` in the Claude Export ZIP**: every export sub-step (Phospho results, GSEA, instrument metadata, TIC QC, excluded files, etc.) now records `[OK] <name>` or `[SKIPPED] <name> -- <reason>` to a manifest written into the ZIP root. So reviewers and the user can see at a glance which sections succeeded, which were skipped, and why — instead of being silently shipped a partial archive. Also surfaces the DE-LIMP version + pipeline used (`dpc` / `maxlfq`) at the top of the manifest. Final console line now reads `Claude export: N files, prompt M chars, K section(s) skipped (see MANIFEST.txt)` so failures are visible during the run too.
- **5 silent `tryCatch(error = function(e) NULL)` blocks** in `claude_export_content` (Phospho_DE_Results.csv, GSEA_Results.csv, Instrument_Metadata.csv, TIC_QC_Metrics.csv, Excluded_Files.csv) replaced with `safe_section()` calls. Pre-v3.9.15, any of these throwing dropped the section silently from the ZIP; now they're recorded with the actual exception message.

## [3.9.14] — 2026-05-05

### Added
- **`safe_section(manifest, name, expr)` helper in `R/helpers.R`**: replaces the dozens of `tryCatch(error = function(e) NULL)` blocks in export bundlers with a manifest-aware wrapper. On success records `[OK] <name>`; on failure records `[SKIPPED] <name> -- <reason>` to a `MANIFEST.txt` written into the ZIP root. Reviewers can see what's missing and why instead of being silently shipped a 9-file ZIP that should have had 12 files. (Wiring into `claude_export_content` and the Complete export comes in v3.9.15 — this commit only adds the helper + the architectural rules so future code is held to them.)
- **`req_nzchar(...)` helper**: like `shiny::req()` but treats `""` as missing, so empty selectInputs / textInputs don't slip through into `topTable(fit, coef = "")` and similar deep failures.

### Changed
- **CLAUDE.md gains an "Architectural rules" section** codifying the 4 root-cause patterns from the critic audit (no hardcoded pipeline descriptions, tag `%||%` defaults in user-facing text, one definition per concept, silent catch is banned in export paths). The rules exist because all four were violated in early DE-LIMP and produced misleading exports for real analyses.
- **Persistent memory entry** added (`feedback_no_hardcoded_pipeline_descriptions.md`) so future sessions surface the same rules at startup.

## [3.9.13] — 2026-05-05

### Added
- **Reproducibility log now includes a Provenance block** (v3.9.13) with: DE-LIMP version, pipeline used (`dpc` / `maxlfq` / `dpc_with_filter_experimental`), timestamp, input parquet absolute path, input parquet MD5 hash, input parquet size in bytes, and the full `sessionInfo()` printed as comments. So a reviewer running the script in a clean R session can confirm they have the **exact same input file** (`tools::md5sum('report.parquet')` should match) and same R/package versions as the original analysis.

### Fixed
- **DPC-Quant log now includes the upstream `readDIANN()` call** with the user's actual `q.cutoffs` value substituted in, plus `library(limpa); library(limma)`. Previously the log started at `dpcfit <- dpcCN(dat)` assuming `dat` already existed — breaking reproducibility for anyone running the script in a clean R session.
- **Contrast log now also emits the per-contrast `topTable()` + `write.csv()` loop**, so the script produces the same `DE_<contrast>.csv` files DE-LIMP exports.

## [3.9.12] — 2026-05-05

### Fixed
- **Reproducibility log (`Reproducibility_Code.R`) was hardcoded to DPC-Quant**: the `add_to_log("Run Pipeline (Main Analysis)", pipeline_code)` call always emitted `dpcCN()` / `dpcQuant()` / `dpcDE()` template strings regardless of which pipeline actually ran. So users who ran MaxLFQ + limma got a reproducibility log claiming they ran DPC-Quant — exactly the opposite of the file's purpose. Now branches on `input$pipeline_mode`: under MaxLFQ + limma, the log emits the actual code path (`arrow::open_dataset → filter → group_by/summarise → pivot_wider → log2 → limma::normalizeBetweenArrays → coverage filter → limma::lmFit`), with the user's actual Q-value / eQ / pgQ / coverage cutoffs substituted in. The DPC-Quant + experimental-override branch keeps the existing dpcDE template.

## [3.9.11] — 2026-05-05

### Fixed
- **Methods text + Claude/AI export now branch on the durable pipeline flag**: previously these checked `values$pipeline_mode_used`, a reactiveVal that gets overwritten if the user runs a different pipeline after MaxLFQ. The more reliable signal is `values$y_protein$other$pipeline == "maxlfq"`, set inside `build_maxlfq_pipeline()` and travelling with the data matrix itself. Both signals are now consulted (OR), so the methods text always describes the pipeline that produced the matrix the user is looking at.
- **`protein_confidence.csv` and `detection_matrix.csv` no longer included in MaxLFQ exports**: both rely on DPC-Quant artifacts (`$other$standard.error`, precursor-count `n.observations`) that don't apply under MaxLFQ — `n.observations` becomes a 0/1 detection mask, and `$standard.error` is NULL. Including them in MaxLFQ exports produced misleading data. Now skipped via `is_maxlfq_export` guard in both `Complete Export` (server_session.R) and `Claude Export` (server_ai.R).

## [3.9.10] — 2026-05-05

### Fixed
- **Last 4 stale unit tests**:
  - `test-helpers_search.R:61` — `ControlPersist` value was lowered from 300 → 60 (zombie-mux risk reduction); test now asserts on the option name, not the literal value.
  - `test-resume_launcher.R:74–75` — current launcher inserts `--kill-on-invalid-dep=yes` between the dependency and the script path, and uses `afterany` for the 4→5 transition. Test now matches the actual format.
  - `test-search_history.R:120` — the unified activity log doesn't have a `completed_at` column; completion is signalled via `event_type = "search_completed"`. Test now asserts on event_type instead of the missing field.

## [3.9.9] — 2026-05-05

### Fixed
- **Stale unit tests now pass** (revealed when v3.9.8 fixed the testthat-install bug):
  - `test-search_history.R`: path test now expects `activity_log.csv$` not `search_history.csv$` (the search history was unified into the activity log months ago; `search_history_path()` is just an alias for `activity_log_path()`).
  - `test-resume_launcher.R`: dependency-chain test now matches the current `afterany` / `afterok` mix (steps 2→3 and 4→5 use `afterany` so a few OOM/timeout array tasks don't collapse the pipeline; step 3→4 stays `afterok`).
- **Real bug in `update_search_status()` shim**: it accepted a `completed_at` argument but only set `event_type = "search_completed"` when non-NA — the timestamp itself was never written to the row. Now properly persists `completed_at` to the activity log row, matching what callers expect.
- **Added `search_history_headers` back-compat alias** for `activity_log_headers` so legacy code/tests that reference the old name still resolve.

## [3.9.8] — 2026-05-05

### Fixed
- **CI test workflow had been silently failing since pre-v3.9.2**: `install.packages(testthat, ...)` was writing to a per-step library that the next step's `Rscript` couldn't see — so every CI run errored with `Error in loadNamespace: there is no package called 'testthat'` and was reported as a test failure even though no actual test was run. Replaced the manual install step with `r-lib/actions/setup-r-dependencies@v2` which coordinates `.libPaths()` between steps and caches dependencies. None of v3.9.2 → v3.9.7 actually had test regressions; the failures were all this CI plumbing bug.

## [3.9.7] — 2026-05-05

### Fixed
- **QuantUMS sliders were silently filtering DPC-Quant input**: v3.8.0 wired `filter_quantums_parquet()` into the load handlers (load example, load local file, Load from HPC). v3.9.x added pipeline branching but didn't undo this — so when a user loaded a parquet with the sliders set to 0.75 (e.g. from a prior MaxLFQ run), the QuantUMS filter ran at file-load time and limpa's DPC-Quant got pre-filtered input. Output was paper-faithful for MaxLFQ but biased for DPC-Quant. The console gave it away with `[QuantUMS filter] Empirical.Quality >= 0.75 — kept N / M precursors` showing up after a fresh load. All three load handlers now read the **unfiltered** parquet unconditionally; QuantUMS filtering happens only at pipeline run-time inside `build_maxlfq_pipeline()`. DPC-Quant always sees the full parquet regardless of slider position. Switching pipeline modes no longer needs a re-load.

## [3.9.6] — 2026-05-05

### Fixed
- **PCA / MDS "Color by" dropdown still showed "Batch" after rename**: when the user renamed Batch → Year in the covariate panel, the PCA and MDS color selectors stayed on the canonical "Batch" label. Both observers now consult `values$batch_name` and use `setNames("Batch", display_name)` so the dropdown shows the user's chosen label while the underlying value (used for `meta$Batch` lookup) stays canonical.
- **Heatmap top annotation legend** also now reflects the renamed covariate label.

## [3.9.5] — 2026-05-05

### Changed
- **MaxLFQ filter waterfall: cleaner display.** Added a "% of input" column so the cumulative drop from the original parquet is visible at every stage (the existing column is "% of prior" — what survived the previous filter). Stages that don't actually drop anything (e.g. excluded-runs filter when no runs were excluded) are now omitted instead of showing a noisy "100% kept" row.

## [3.9.4] — 2026-05-05

### Added
- **QuantUMS / FDR filter waterfall** in the **Data Completeness** sub-tab when the MaxLFQ pipeline ran. Compact table showing precursor row counts at every filter stage: `Input → after FDR → after eQ → after pgQ → after excluded-runs`, with a "% kept" column at each step. So users can see exactly how aggressive their filters are. The same waterfall is now logged to the R console (one line per stage with absolute counts and percentages of the prior pool). Helper `build_maxlfq_pipeline()` materialises per-filter row counts via lazy Arrow `summarise(n=n())` scans (cheap, doesn't collect the whole filtered dataset).

## [3.9.3] — 2026-05-05

### Fixed
- **Heatmap crashed under MaxLFQ with `NA/NaN/Inf in foreign function call`**: `hclust()` doesn't tolerate NaN, but the MaxLFQ Z-score matrix (per-row centring of a matrix with NAs) produces NaN for rows that have any missing value. Now drops rows with `< 2` non-NA Z-scores and zero-fills the rest before passing to `Heatmap()` (clustering only — values are still NA in the source matrix). Applied across all four heatmap render/export sites in `R/server_de.R`.
- **On/Off Proteins table threw a Shiny renderWidget error** when results were present: stripped suspect attrs/options from the DT call (dropped `filter = "top"`, dropped `htmltools::tags$caption()` wrapper in favour of a plain string caption, force-cast every column to a simple atomic vector before passing to `DT::datatable`).

### Changed
- **Stacked-bar title in Data Completeness flips per pipeline**: under DPC-Quant the bar still says "Detected vs Inferred Proteins per Sample" (missing values are filled in by the probability model). Under MaxLFQ + limma the bar now says "Detected vs Missing Proteins per Sample" — those cells are genuinely missing, not imputed. Legend label and tooltip text track the title.

## [3.9.2] — 2026-05-05

### Fixed (BLOCKERS)
- **MaxLFQ branch ignored `values$excluded_files`**: pivoted across every Run in the parquet, so users who excluded samples in the metadata table silently got those samples back, with the unwanted runs flowing through quantile normalization and lmFit. Fixed: `build_maxlfq_pipeline()` accepts a `keep_runs =` argument; the run_pipeline observer passes `meta$File.Name` so the matrix matches the metadata exactly.
- **Sample-to-group misalignment under MaxLFQ**: `meta <- meta[colnames(dat$E), ]` ran unconditionally, but under MaxLFQ `dat` is the wrong matrix. With `lmFit` matching by column position, this could silently assign samples to wrong groups. Now branches on `values$pipeline_mode_used`: under MaxLFQ uses `colnames(values$y_protein$E)`. Hard fail-fast notification if any matrix sample lacks a metadata row.
- **`compute_onoff_proteins()` contrast-string fragility**: was reconstructed by splitting `forms` on `" - "`, which could break if any group level contained that substring. Now accepts the `combs` matrix from `combn(levels(groups), 2)` directly. Rows of `combs` are flipped before passing so the on/off Contrast string matches limma's `g2 - g1` convention.

### Added (statistician's recommendation + reviewer HIGH #4)
- **Coverage filter for the MaxLFQ pipeline**: new sidebar slider `coverage_min_frac` (default 0.5 = 50% of samples must have a non-NA value, matching the UC Davis Bioinformatics Core's limma-proteomics tutorial). Proteins below the threshold are dropped from the limma fit (so eBayes isn't moderating against rows with 1-2 finite values) but still appear in the On/Off Proteins sub-tab as presence/absence calls. Console logs the kept/dropped counts. Set to 0 to disable.

### Changed
- **Methods text branches on the pipeline that actually ran** (reviewer HIGH #6). Under MaxLFQ + limma, the methodology paragraph now describes the Moschem 2025 pipeline accurately (filter → PG.MaxLFQ pivot → log2 → quantile-norm → lmFit + eBayes; coverage filter percentage; on/off panel for fully-missing proteins) and cites the paper explicitly. Under DPC-Quant the existing limpa methodology stands.
- **`compute_onoff_proteins()` accepts both list-of-pairs and matrix-of-contrasts** for caller convenience.

## [3.9.1] — 2026-05-05

### Fixed
- **MaxLFQ pipeline was under-normalized**: v3.9.0 only median-centered the log2(PG.MaxLFQ) matrix before handing it to limma. Median-centering aligns medians but lets between-sample variance differences leak through, and that crushed eBayes power — users reported only a handful of DE proteins on a 230-sample dataset where dozens-to-hundreds were expected. Switched to **`limma::normalizeBetweenArrays(method = "quantile")`**, which is the standard cross-sample normalization used by FragPipe-Analyst, Spectronaut/MSstats, and DIA-NN's own analyzer. The pre-normalization matrix is preserved as `y_protein$other$E_log2_raw` so the Norm QC tab can show before/after.
- **Norm QC tab showed stale DPC-Quant output under MaxLFQ pipeline**: the pipeline diagnostic (`generate_norm_diagnostic_plot`) was hard-wired to `values$raw_data$E` (precursor input) vs `values$y_protein$E` (DPC-Quant output) regardless of which pipeline ran. Now branches on `values$pipeline_mode_used`: under MaxLFQ it shows pre-quantile-norm log2(PG.MaxLFQ) vs post-quantile-norm matrix; under DPC-Quant it shows the original DIA-NN→DPC-Quant view. Subtitle adapts to indicate which pipeline produced the plot.

## [3.9.0] — 2026-05-04

### Added
- **MaxLFQ + limma pipeline** (paper-faithful Moschem et al. 2025). New radio in **Pipeline Settings → Quantification method**: `DPC-Quant (limpa, default)` vs `MaxLFQ + limma (Moschem 2025)`. When MaxLFQ is selected, DE-LIMP **bypasses limpa entirely**: reads the parquet via Arrow, applies the Q.Value + QuantUMS filters, pivots to `Protein.Group × Run` using DIA-NN's already-computed `PG.MaxLFQ`, log2-transforms, median-normalises, and runs plain `limma::lmFit` with NAs left in place per the paper. New helper `build_maxlfq_pipeline()` in `R/helpers.R`.
- **On/Off Proteins sub-tab** in the DE Dashboard. Surfaces proteins detected in ≥ N samples of one group AND zero in the other — these get NA logFC under limma and silently drop from the volcano. Slider for the N threshold (default 2), DT table with filter / sort / CSV download, helpful empty-state message under DPC-Quant ("its missing-data model fills these in"). New helper `compute_onoff_proteins()` in `R/helpers.R`.
- **Experimental escape hatch**: when MaxLFQ + limma is selected, a checkbox **"Run filtered precursors through limpa anyway (experimental)"** lets users force DPC-Quant on the QuantUMS-filtered parquet. Yellow warning banner appears on the run because neither paper tested this combination. Default unchecked.
- **Pipeline-aware QuantUMS panel banner**: a dynamic note inside the QuantUMS-filters details block tells the user whether the sliders actually do anything based on the active pipeline mode (greyed-out warning under DPC-Quant, green confirmation under MaxLFQ + limma).
- **Methodology disclosure**: when MaxLFQ pipeline runs, the console prints `[DE-LIMP] MaxLFQ pipeline: N proteins x M runs, X cells missing (P%). Filters: ...` for full reproducibility.

### Changed
- **QuantUMS info modal** rewritten for v3.9 — explains the three resulting paths (DPC-Quant / MaxLFQ + limma / experimental combo), and points users at the new On/Off Proteins panel for proteins lost to all-missing-in-one-condition behaviour.
- **Methodological note in the modal**: DPC-Quant's design assumption (use low-quality precursors via probability modelling) is incompatible with pre-filtering, so under DPC-Quant the QuantUMS filters are forced to 0 internally regardless of slider value. Switching pipelines is the user-facing way to enable filtering.

## [3.8.5] — 2026-05-04

### Fixed
- **Template Import broken when covariate columns had been renamed**: Exporting the template after renaming Batch → Year (etc.) wrote the renamed column headers into the CSV, but the import handler hard-required the canonical names (`Batch`, `Covariate1`, `Covariate2`), so the round-trip failed with "Template must have columns: ID, File.Name, Group, Batch, Covariate1, Covariate2". Now: (a) export writes **canonical** column headers regardless of the user's display rename and prefixes a `# Display labels: …` hint comment for human readability; (b) import is tolerant — accepts canonical names, the user's current display-rename names, OR positional fallback (columns 4/5/6 mapped to Batch/Cov1/Cov2). The error message now also lists the actual columns it found, so when the import does fail the user can see why.
- **Covariate panel layout broken on some browsers (CSS grid alignment)**: Three checkboxes were rendering in the wrong rows or below their text inputs, leaving Batch with no visible checkbox. Replaced the CSS-grid layout with explicit per-row flexbox containers (`.cov-row`) plus a tighter Bootstrap-checkbox CSS reset (`margin: 0 !important; padding: 0 !important; min-height: 0 !important;` on `.form-group`, `.checkbox`, and the wrapping `label`). Renders consistently on Safari / Chrome / Firefox now.

## [3.8.4] — 2026-05-04

### Fixed
- **Cryptic "NA/NaN/Inf in 'y'" when a covariate caused a rank-deficient design**: A user ticked "In model" for a `Run order` covariate (per-sample numeric IDs like 707, 708, 813, 16437…) and `Student` (with one level appearing in only one sample). The old design builder factor-expanded every covariate, producing a 200+-column rank-deficient design matrix; limma then printed "Coefficients not estimable" for nearly all of them, generated NaN coefficients, and `eBayes()` died with `NA/NaN/Inf in 'y'`. The user-facing error said only "Differential expression failed" with no clue which covariate to fix.

### Added
- **Auto-detect numeric vs categorical covariates** (`coerce_covariate_column()` in `R/helpers.R`). When a covariate column parses cleanly as numeric AND has ≥ 5 distinct values, DE-LIMP now enters it into the design matrix as a single continuous coefficient instead of factor-expanding it. So `Run order = 707, 708, 813, …` becomes one coefficient (a linear drift term) rather than 230 columns. Console message: `[DE-LIMP] Covariate 'Run order' treated as numeric (228 distinct values).`
- **Singleton-level detection** for factor covariates. If any level of a factor covariate appears in only one sample (e.g. `Student = A` for a single row), DE-LIMP now skips that covariate with a named warning instead of letting it silently break the design: `Covariate 'Student' has 1 level(s) that occur in only one sample (A) — these break the model. Either drop those rows or merge them into another level.`
- **Pre-flight design-rank check** (`diagnose_design_rank()` in `R/helpers.R`) runs immediately before `limpa::dpcDE()`. If the design is rank-deficient, DE-LIMP refuses to fit and surfaces a notification that names the offending coefficients, suggests the fix (untick the covariate, drop singleton-level rows), and notes that QC / Expression Grid / PCA still work.

## [3.8.3] — 2026-05-04

### Added
- **Soft "star DE-LIMP" footer in three high-leverage exports**, complementing the v3.8.2 About-tab nudge:
  - **Methods text** (`build_methodology_text()`) — appended to the CITATION block, so anyone who exports a methods section sees the link.
  - **AI Summary HTML/Markdown report** — single italic line in the footer next to the existing "Generated by DE-LIMP" line.
  - **Claude Export ZIP `PROMPT.md`** — one-line italic footer at the bottom.
  All three use the same wording: "If DE-LIMP helped your work, a star on GitHub helps other proteomics labs find it." No popups, no nags — just a contextual link in places where the user has already extracted value.

## [3.8.2] — 2026-05-04

### Added
- **Soft "star DE-LIMP on GitHub" nudge** in the About → Community tab. One small line under the existing GitHub-stats cards, no popup, no timer, no nagging — just a contextual link sitting next to the data it relates to. Wording: "If DE-LIMP helped your work, a star on GitHub helps other proteomics labs find it. ★ Star DE-LIMP →".

## [3.8.1] — 2026-05-04

### Added
- **Navbar version badge**: A pill-shaped `vX.Y.Z` badge sits at the far right of the navbar, just before the gear icon. It's visible on every screen — useful for HF / WSL / Docker users so they can confirm at a glance which release they're running. Clicking the badge opens the GitHub CHANGELOG in a new tab. Reads directly from the `VERSION` file at UI-build time (no reactive plumbing needed).

## [3.8.0] — 2026-05-04

### Added
- **QuantUMS quality filters (eQ + pgQ) on parquet load**: Optional precursor pre-filter that drops rows below user-set cutoffs on DIA-NN's `Empirical.Quality` and `PG.MaxLFQ.Quality` columns *before* `limpa::readDIANN()` ingests the file. Implements the recommendations of Moschem et al. *J. Proteome Res.* 2025, 24:3860 ([10.1021/acs.jproteome.5c00009](https://doi.org/10.1021/acs.jproteome.5c00009)). The qQ score is intentionally not exposed — the paper demonstrates it has negligible impact.
  - New helper `filter_quantums_parquet()` in `R/helpers.R` reads the parquet via Arrow, drops sub-threshold rows, writes survivors to a temp parquet. No-ops when both cutoffs are 0 (the default), so existing behaviour is unchanged.
  - Two new sidebar inputs under the new **QuantUMS quality filters** collapsible panel: `eq_cutoff` (Empirical Quality ≥) and `pgq_cutoff` (PG.MaxLFQ Quality ≥), both defaulting to 0.
  - New `?` info modal explains the paper, recommended thresholds (0.75 / 0.75 = best ROC AUC; 0.9 = too aggressive), the difference between QuantUMS and the existing identification Q-value, and the caveat that pgQ is MaxLFQ-derived while DE-LIMP uses DPC-Quant for rollup.
  - Wired into all three load paths: example data, local upload, and Load from HPC. Console logs how many precursors survive, e.g. `[QuantUMS filter] Empirical.Quality >= 0.75 — kept 89,432 / 114,166 precursors (78.3%)`.
  - Filter description (e.g. `c("Empirical.Quality >= 0.75")`) is stored in `values$quantums_filter_applied` for downstream display in Methods text and AI export.

### Documentation / clarification (no code change in 3.8.0)
- Confirmed DE-LIMP's existing identification Q-value filtering is correctly handed to limpa: `q.cutoffs = input$q_cutoff` (default 0.01) is applied to `c("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value")` at the precursor row level before protein rollup — the same 1% FDR strategy Moschem et al. use, with the small column-choice difference that limpa filters per-run *precursor* FDR (`Q.Value`) while the paper filters per-run *protein-group* FDR (`PG.Q.Value`); both anchor 1% FDR and neither is wrong.

## [3.7.9] — 2026-05-04

### Changed
- **Covariate panel UX overhaul** on the *Assign Groups & Run* sub-tab: the previous compact strip put two unlabelled checkboxes next to a “Covariates:” heading and three rename text inputs that wrapped onto separate lines, so it was unclear which checkbox controlled which column or what either control actually did. Replaced with a labelled 2-column grid: an `In model` column header sits above three checkboxes (each tooltipped “Add this covariate to the DE design matrix”), a `Column name (click to rename)` header sits above three uniformly-sized text inputs (each tooltipped to clarify they only rename, not include). A new `?` info button next to the panel title opens a modal explaining covariate semantics: tick = added to design matrix; text box = column rename; when to include batch/biological covariates; and the R/limma equivalent (`design <- model.matrix(~ 0 + Group + Year)`). One-line tip below the grid summarises both controls so casual users don’t have to open the modal.

## [3.7.8] — 2026-05-04

### Fixed
- **Batch covariate label rename didn't update the metadata table column header**: Renaming "Covariate1" / "Covariate2" via the sidebar text inputs propagated to the rhandsontable column headers, but renaming "Batch" (e.g. to "Year") did not — `colnames(display_df)` and the `hot_col(...)` calls hardcoded the literal string `"Batch"` while cov1/cov2 used the dynamic display strings. The `input$batch_label` value is now read alongside `cov1_label`/`cov2_label`, stored in `values$batch_name`, and used as the displayed column header. Internal column name in the metadata data.frame stays `Batch` (the round-trip in `colnames(tbl) <-` already normalises it back), so downstream code referencing `meta$Batch` is unaffected.

## [3.7.7] — 2026-05-04

### Fixed
- **`limpa::readDIANN()` failed with "nanoparquet required but not installed"**: limpa lists `nanoparquet` in `Suggests` rather than `Imports`, so the direct-repo limpa install (Path 2 fallback added in v3.7.2) didn't pull it in, and the first parquet read on a fresh install errored. Added `nanoparquet` to `core_pkgs` so it's installed automatically alongside the rest of the stack.

## [3.7.6] — 2026-05-04

### Fixed
- **Load from HPC appeared frozen for large reports**: With a 1.4 GB / 231-sample `report.parquet`, the post-download phase (`get_diann_stats_r`, `limpa::readDIANN`, normalization detect, phospho detect) is genuinely synchronous and CPU-bound, so Shiny's progress bar can't tick while it runs. The console also had only a single message at the end of the SCP, making it impossible to tell whether the load was hung or just slow. Added per-phase tick lines (`[DE-LIMP] <phase> ... [HH:MM:SS]` then `↳ <phase> done in N.Ns`) for QC stats, expression-matrix read, normalization detect, and phospho detect, with `flush.console()` so each line appears immediately even while the next phase is blocking. Now you see live heartbeats in the RStudio console rather than silence.

## [3.7.5] — 2026-05-04

### Fixed
- **`Load from HPC` failed silently on large `report.parquet` files**: `scp_download()` and `scp_upload()` had a hard-coded 60-second `processx` timeout, which killed any transfer larger than what fit in 60 seconds — for example the 1.4 GB short-course `report.parquet`. The on-screen "SCP download failed:" message was also blank because the caller printed `dl_result$stderr` while the helper returns `stdout` (stderr is folded into it). Now: (a) default timeout raised to 1800 s (30 min), with an optional `timeout` argument so callers can tune; (b) timeouts surface a real message ("Transfer exceeded N s timeout — increase the timeout argument") instead of an empty string; (c) the Load-from-HPC notification reads the correct `stdout` field, so future transfer errors show the actual reason. Toast duration bumped to 15 s so the message is readable.

## [3.7.4] — 2026-05-04

### Documentation
- **Beginner-friendly Mac install walkthrough in README_GITHUB.md**: Replaced the lean two-command install with a four-step guide for first-time users — install R (with arm64 vs x86_64 guidance), install RStudio, install Git via Xcode CLT, then `git clone` + `shiny::runApp`. Added a "Troubleshooting first launch" subsection covering the BiocManager-stale failure mode (now self-healing as of v3.7.3), R-too-old, permission-denied, and port-in-use. Power-user single-command path retained below the walkthrough.

## [3.7.3] — 2026-05-04

### Fixed
- **Stale-BiocManager fallback now also covers core/optional Bioc packages**: v3.7.2 added a direct-repo install path for limpa, but the next install block (`ComplexHeatmap`, `AnnotationDbi`, `ggridges`, `clusterProfiler`, `enrichplot`, `org.Hs.eg.db`, `org.Mm.eg.db`, `MOFA2`, …) still went through `BiocManager::install()` only. On R 4.6 with stale BiocManager that block silently failed and the app died on "Missing required packages". The R↔Bioc mapping and direct-repo helper are now hoisted to module scope and reused: missing-packages install runs Path 1 (BiocManager) → Path 2 (direct Bioc repo URL) before checking what's still missing. Should self-heal for everyone after a fresh R release.

## [3.7.2] — 2026-05-04

### Fixed
- **limpa install on freshly released R versions**: When R is newer than BiocManager's hardcoded R↔Bioc map (common right after a major R release — e.g. R 4.6 / Bioc 3.23), `BiocManager::version()` returns an unresolved value and the install branch silently fell through, then printed a misleading "R upgrade needed" message. The install branch now (1) correctly probes BiocManager's resolved-vs-unresolved state, (2) falls back to a direct-repo `install.packages("limpa", repos = ".../packages/<bioc>/bioc")` install bypassing BiocManager, and (3) when it does fail, the error message names the real cause (R-too-old vs BiocManager-stale vs network) and prints a copy-pasteable fix command sized to the user's R version.

## [3.7.1] — 2026-05-04

### Added
- **Startup version banner**: `app.R` now prints `DE-LIMP vX.Y.Z | R x.y.z | timestamp` to the console before any package work, so the running version is visible at a glance in the RStudio console.

### Fixed
- **Misleading "R upgrade needed" message**: When `BiocManager::version()` couldn't reach Bioconductor's online validator, the limpa install branch fell through and printed a generic R-version-too-old notice. Network/validator hiccups are now no longer conflated with version mismatches (see follow-up: refine the failure-message branching).

## [Unreleased] — Post-3.7.1 Development

### Added
- **DPC-Quant Detection Transparency**: Expression Grid tooltips show nObs, SE, and 95% CI per cell. Violin plots mark inferred estimates (nObs=0) with hollow markers. New `Detection_Class` column (Complete/Partial/Sparse/Inferred) in exported data. `protein_confidence.csv` included in all exports (session, Claude ZIP).
- **Run Comparator Claude Export Enhancements**: All new data files (protein_confidence.csv, comparator context, library info, run QC) included in Claude ZIP export. DPC-Quant methodology note added to comparator prompt.
- **Normalization Mismatch Detection**: Run Comparator detects when two runs used different normalization strategies and flags it in settings diff.

### Fixed
- **Hugging Face subscript error**: Data Completeness visualization crashed on HF due to column subsetting with string names on a matrix. Switched to numeric index subsetting.
- **Load from HPC visible on HF**: "Load from HPC" button now hidden when running on Hugging Face (no SSH available).
- **Violin modal broken "Back to Grid" button**: Removed non-functional back button from Expression Grid violin modal.
- **Jaccard dendrogram rendering**: Fixed y-axis range clipping and invisible legend markers in Data Completeness dendrogram plot. Added `ggdendro` to Dockerfile.
- **Comparator maxLFQ reference**: Corrected incorrect maxLFQ terminology in comparator output.
- **FASTA info for prebuilt speclib**: Populated `fasta_info` from library catalog when using a prebuilt spectral library (was showing blank).
- **Comparator Claude export row mismatch**: Collapsed nested `library_info` list entries that caused extra rows in exported CSV.
- **Spectronaut 20+ RunOverview format**: Parser now handles Spectronaut 20+ key-value format in RunOverview (changed from tabular to key-value layout).
- **'Imputation' terminology**: Replaced incorrect "imputation" with "probabilistic estimation" throughout, matching limpa's DPC-Quant methodology.

### Performance
- **SSH connection test**: Uses 10-second timeouts and tries fast probes first, reducing initial connection time from 30s to ~5s.
- **SSH file browser**: Starts browsing at specific subdirectories (raw data, FASTA) instead of filesystem root.
- **Expression Grid pagination**: Shows 50 rows with vertical scroll instead of rendering all rows at once.
- **Async cluster resource check**: Cluster CPU/memory check on SSH connect no longer blocks the Shiny event loop.

### SLURM Queue Switching (March 30)
- **Auto-queue switch fixed**: Now queries `sinfo -p low` directly each monitoring cycle instead of stale snapshot from SSH test.
- **QOS set on partition move**: `slurm_move_job()` sets `{account}-{partition}-qos` alongside Account/Partition. Fixes `InvalidQOS` on moved jobs.
- **Requeue on low partition**: Moved jobs get `Requeue=1` so preempted tasks auto-restart.
- **PREEMPTED state handling**: `check_slurm_status()` maps PREEMPTED/REQUEUED → "queued", NODE_FAIL → "failed".
- **QOSMaxCpuPerUserLimit detection**: Per-user CPU limits now trigger auto-switch (not just InvalidQOS).
- **Partial move tracking**: Split-partition jobs (steps 2/4 on low, 1/3/5 on high) can now move back when capacity returns.
- **pending_reason for queued step**: Fetches reason for first queued step, not first non-completed (running) step.
- **Retry dependency chain**: After partial step 2 retry, updates Step 3's SLURM dependency via `scontrol update` to wait for retry completion.
- **Retry/queue events logged to search_info.md**: Timestamp, reason, new job IDs, failed tasks, memory changes recorded.
- **Job queue GUI accuracy**: Old failed jobs hidden when retry creates new entry. Progress counter aggregates original + retry completed tasks.

### Spaces in Paths (March 30)
- **Launcher script**: Quoted all sbatch script paths in `submit_all.sh` for directories with spaces.
- **sbatch scripts**: Quoted paths in `#SBATCH -o/-e` directives and `apptainer exec --bind` arguments.
- **SSH launcher execution**: Quoted remote path in `bash` command.

### Publication & Export (March 27-30)
- **SVG vector export**: Camera icon on Volcano, PCA, CV scatter (plotly). SVG download buttons on Heatmap and Violin plots (ggplot2/ComplexHeatmap).
- **AI Summary export**: Changed from HTML to Markdown format.
- **NCBI gene symbols in DE table**: DE Results Table and Volcano now use gene_map.tsv fallback for NCBI RefSeq accessions.
- **GSEA with NCBI data**: Converts RefSeq accessions to gene symbols via gene_map before SYMBOL→ENTREZID mapping.
- **Default Gemini model**: Changed to `gemini-2.5-flash` (production, was preview).
- **Live GitHub stats on HF**: Stars/forks fetched from GitHub API at startup (always fresh).
- **GSEA on HF**: Installed clusterProfiler/enrichplot/ggtangle/org.Hs.eg.db in Dockerfile.

### Drift Test Enhancements (March 24-30)
- **Daily runs**: Changed from weekly to daily for ASMS poster data collection.
- **Model tracking**: Captures model ID, input/output tokens, response time from API.
- **Per-protein fold changes**: Extracts FC values from 4 patterns (logFC=N, N logFC, range, parenthetical).
- **Language analysis**: Hedging vs confident language counts + example quotes.
- **Gene stability**: Core/frequent/one-off gene tracking, genes gained/lost between runs.
- **Trend assessment**: Linear regression on overlap, word count, gene count (4+ runs).
- **Overall verdict**: HEALTHY/ACCEPTABLE/ATTENTION NEEDED (4+ runs).
- **Fail loudly**: API key check, baseline creation check, no silent `|| true`.
- **Golden baselines committed**: `git add -f` bypasses .gitignore for .rds files.

### Documentation
- Added SSH XIC viewer and DPC-Quant confidence overlay plans to TODO.
- Revised confidence overlay design per statistician review (saturation/opacity, not red/green).
- **QUEUE_SWITCHING.md**: Comprehensive documentation of auto-queue logic, known issues, state mapping.
- **DRIFT_TEST_METHODOLOGY.md**: Full methodology for ASMS poster — metrics, ground truth, trend assessment.
- **TODO from expert reviews**: 20 items from biologist, proteomics expert, statistician reviews.
- **CI workflows fixed**: Added missing R packages (processx, stringr, dplyr) to drift test and unit test workflows.

## [3.7.0] - 2026-03-17

### Added
- **Docker Launcher for Windows** (`Launch_DE-LIMP_Docker.bat`): One-click batch file for Windows lab PCs. Auto-detects SSH keys, supports shared PC accounts (multiple Windows users), copies SSH keys into Docker-accessible volume, handles permissions. Skips rebuild on every launch for faster startup.
- **SSH Auto-Connect on Startup**: When an SSH key is detected (Docker volume mount or `~/.ssh/`), the app automatically connects to HPC on startup — no manual "Test Connection" click needed. Stale ControlMaster socket detection prevents hangs.
- **Environment Badge**: Colored badge in navbar showing deployment mode — Docker (red), HPC/Apptainer (green), Local (blue), Hugging Face (orange). Auto-detects Apptainer container environment.
- **SLURM Proxy for Apptainer**: All 9 SLURM command paths (`sbatch`, `squeue`, `scancel`, `sacct`, `sinfo`, `sacctmgr`, `scontrol`, `srun`, `sbatch --test-only`) proxied from inside Apptainer container to host via a relay process. Cluster monitor works inside containers.
- **Shared Storage for All HPC Files**: All DE-LIMP files (container, R packages, git repo, data, results) live on `/quobyte/proteomics-grp/de-limp/` to avoid home directory quota limits. Per-user subdirectories (`users/{username}/logs/`, `users/{username}/jobs/`) prevent multi-user conflicts.
- **Per-User Directories for Multi-User HPC**: Multiple users can run DE-LIMP simultaneously without conflicts. SLURM logs and generated scripts go to per-user dirs on shared storage.
- **Apptainer Cache Redirect**: `APPTAINER_CACHEDIR` set to shared storage, avoiding home directory quota exhaustion during container pulls.
- **Code Update Detection Banner**: App detects when local code is behind the git repo and shows an update notification banner.
- **Home Directory Quota Warning**: Startup check warns if home directory usage exceeds 80% (common HPC issue).
- **Container Detection**: Skips BiocManager package validation when running inside a container without internet access (offline-safe startup).
- **NCBI Proteome Download**: New "Download from NCBI" option in FASTA selector. Search NCBI Datasets by organism name, select a proteome, download RefSeq protein FASTA. Supports all organisms with NCBI reference proteomes.
- **NCBI Gene Symbol Mapping**: Batch E-utilities lookup maps RefSeq accessions (XP_, NP_, WP_) to proper gene symbols. Gene map TSV auto-downloaded to HPC via SSH for Docker users who lack direct E-utilities access. Mapping applied to expression grid, volcano plots, and all downstream analysis.
- **Contaminant Analysis Subtab**: New subtab in Data Overview with summary cards (contaminant count, % of total, median intensity ratio, keratin count), per-sample stacked bar chart, top contaminants table with keratin flagging, and contaminant heatmap (top 20 by median intensity).
- **Expression Grid Contaminant Highlighting**: Contaminant protein rows highlighted pink/red in the Expression Grid for visual identification. `Cont_` prefix proteins flagged automatically.
- **Signal Distribution Contaminant Overlay**: Checkbox to overlay contaminant proteins in orange on the Signal Distribution plot, showing where contaminants fall relative to endogenous proteins.
- **SSH File Browser**: Visual directory browser for Remote (SSH) mode. Browse buttons for raw data and FASTA directories open a modal with clickable breadcrumbs, Up/Home navigation, color-coded entries (folders blue, data files green, other grey). Replaces manual path entry.
- **Load from HPC Button**: "Load from HPC" button in sidebar opens SSH file browser filtered for `.parquet` files. SCP downloads the selected file and automatically loads it through the pipeline.
- **Remote Activity Log**: Activity log stored on shared HPC storage for multi-user visibility. History tab reads remote activity log via SSH when connected. Source badge shows whether history comes from local or remote storage. File locking for concurrent multi-user writes.
- **Configurable File Browser Roots** (`DELIMP_EXTRA_ROOTS` env var): Additional root directories for the SSH file browser, allowing any HPC to configure institution-specific data paths.
- **Pre-Staged FASTA Directory**: FASTA files on shared storage (`/quobyte/proteomics-grp/de-limp/fasta/`) appear as a dropdown option — fastest way to select commonly used organisms on HPC.

### Fixed
- **No-replicates mode**: When groups have fewer than 2 replicates, quantification completes normally but DE analysis is skipped gracefully instead of crashing. Users see an informational message explaining that DE requires replicates.
- **Expression Grid without DE results**: Grid now renders with quantified data even when DE has not been run (e.g., no replicates). Missing P.Value column added as fallback.
- **PCA tab visible without DE**: PCA works on quantified expression data and no longer requires `values$fit` to be non-NULL.
- **Parallel pipeline quant_verify_block**: Verification block no longer corrupts `file_list.txt` by appending to it. Uses a separate temp file for the missing files list.
- **Step 4 skips failed Step 2 files**: Array jobs in Step 4 now skip files that failed in Step 2 instead of failing the entire step.
- **Auto-adjust search CPUs**: Search CPU count automatically reduced to match available per-user SLURM limits, preventing job rejection.
- **Default to HPC backend in Docker**: When SSH key is detected inside Docker container, defaults to HPC backend instead of local.
- **Hardcoded username/paths removed**: All hardcoded personal directory references and test paths removed from UI defaults.
- **NCBI protein links**: Proteins with NCBI RefSeq accessions (XP_, NP_, WP_) link to NCBI Protein instead of UniProt. `Cont_` prefixed proteins link to their source database.
- **Gene map download via SSH**: When NCBI gene map TSV is not found locally (Docker users), it is automatically downloaded from HPC via SSH.

### Changed
- **Recommended deployment**: Docker on Windows + SSH to HPC is now the primary recommended approach. Apptainer on HPC is documented as an alternative.
- **HPC directory layout**: All files moved from `~/containers`, `~/DE-LIMP`, `~/R/delimp-lib` to shared storage at `/quobyte/proteomics-grp/de-limp/`.
- **File browser performance**: SSH file browser uses specific subdirectory roots instead of scanning entire filesystem. Optimized for large HPC directory structures.
- App version bumped to v3.7.0.

## [3.6.1] - 2026-03-11

### Fixed
- **Spectronaut Candidates.tsv parsing**: Protein column regex now matches `Group` (Spectronaut's actual column name, not `ProteinGroup`). Comparison column regex matches `Comparison (group1/group2)` format (removed `$` anchor). Both fixes restore proper DE concordance analysis.
- **DE Concordance NaN crash**: Spectronaut proteins with 0 `# of Ratios` have NaN for logFC/Pvalue/Qvalue. `classify_de()` and `assign_hypothesis()` now handle NaN/non-finite values safely, preventing "missing value where TRUE/FALSE needed" errors.
- **Spectronaut version showing "unknown"**: Added fallback chain through `spectronaut_version` and `Software Version` keys from AnalysisOverview parsing.
- **Spectronaut precursors "not available"**: Parse `AnalysisOverview.txt` (handles Spectronaut's `AnalyisOverview` typo) for "N of M" format precursor/protein counts. Enriches `library_info`.
- **DE Concordance error message referenced FragPipe in Spectronaut mode**: Mode-aware messaging now shows correct tool name.
- **History tab slow to populate**: Replaced 7 independent `activity_log_read()` calls (network CSV reads) with single `cached_activity_log()` reactive that caches per invalidation cycle.
- **LC/EvoSep info missing from Methods text**: `format_instrument_methods_text()` now includes EvoSep SPD and gradient length, with deduplication when method name already contains SPD info.
- **Job queue crash (`vapply: values must be length 1`)**: Corrupt job entries with NULL/empty `status` field crashed vapply and switch. Added `sanitize_job()` validator on load/save with null-safe guards on all vapply/switch calls.
- **Methodology text overflow**: Text didn't fit in window and couldn't scroll. Added `pre-wrap`/`word-wrap` CSS and scrollable wrapper div.
- **Instrument metadata lost on history load**: `values$instrument_metadata` was NULL after loading from history because session.rds was remote-only. Added fallback recovery from job queue entries.
- **Mounted drive dependency removed**: All app state files (activity log, cluster usage, lab members) now always use local `~/.delimp_*` paths. SMB mounts may be absent, slow, or disappear — no longer a failure point.
- **FDR mislabeling in comparator**: Split `fdr_threshold` into `de_significance` and `identification_fdr` — Spectronaut's `Protein Group FDR: 0.01` is identification FDR, not DE significance. Previously caused Gemini to incorrectly flag as main discordance driver.
- **SSH auto-connect hanging**: Stale ControlMaster socket from previous app instance caused SSH to hang. Added `ssh -O check` probe before auto-connect, removes dead sockets.
- **Cluster monitor CSV header mismatch**: `cluster_usage_headers` had 13 columns but data rows had 16 (schema evolved). Added auto-repair on read.
- **"Prepare Next Analysis" not working on HF**: Observer was in `server_search.R` which early-returns when `search_enabled=FALSE`. Moved to `server_session.R`.
- **Faceted TIC plot unreadable with many files**: With 209 files, panels were tiny and labels overlapped. Now shows only flagged runs when >40 files, uses adaptive column count, truncated labels, status-colored annotations, and dynamic plot height.

### Added
- **Rescue stats for Spectronaut 0-ratio proteins**: Detects proteins with 0 computable ratios in Spectronaut that DE-LIMP could still test. Imputation-aware messaging (None/enabled/unknown). New Rule 0 ("Untestable in Spectronaut") as highest-priority hypothesis. Summary banner shows rescue count below main stats.
- **Contrast mismatch warning**: Fuzzy matching detects when Spectronaut conditions don't match DE-LIMP contrasts. Amber warning div displayed in DE Concordance sub-tab.
- **Instrument context in comparator prompts**: Both `build_gemini_comparator_prompt()` and `build_claude_comparator_prompt()` now include brief instrument/LC context (model, LC system, method, SPD, gradient length) when `instrument_metadata` is available.
- **Debiased Gemini comparator prompt**: Rewrote prompt to be objective — balanced tool descriptions, structured 8-section output template (Factual Observations, Sources of Disagreement, Case for A/B, Settings Audit, Concordant Biology, Synthesis, Follow-ups), debiasing guidelines, neutral Quant3 framing. Adds pre-filtering context (Spectronaut LFC candidate filter, imputation strategy).
- **Remote HyStarMetadata.xml extraction**: SSH file scan now downloads `HyStarMetadata.xml` alongside `analysis.tdf` for LC method/system/runtime extraction from remote timsTOF data.
- **Unified activity log** (v3.6.0): Single CSV replacing dual search/analysis history CSVs + projects.json. 33 columns, append-only with file locking. One-time migration from old format.
- **Per-user cluster monitoring**: Tracks lab member CPU/memory usage on dual partitions with historical snapshots. Queue wait time display. Auto queue switching (genome-center-grp/high → publicgrp/low).
- **Compare from History**: Select 2 analyses in history table → auto-load into Run Comparator.
- **Session auto-save**: Deterministic RDS saved to `{output_dir}/session.rds` after pipeline completion.
- **View Prompt button**: Users can inspect the exact Gemini prompt in the comparator AI Analysis tab before/after running analysis. Includes copy-to-clipboard.
- **SLURM estimated start time**: Queued jobs display `squeue --format=%S` estimated start time in the job queue UI.
- **Per-job wait time logging**: Records queue-to-running transition time per job for grant justification. `record_job_wait()` / `read_job_wait_log()` in helpers_search.R.
- **Exclude Failed TIC button**: One-click removal of failed TIC runs from search file list.
- **TIC trace recovery**: Job queue entries store `metadata$tic_traces` and `metadata$tic_metrics`. On file scan, checks if matching TIC data exists in queue to avoid re-extraction.

### Changed
- App version bumped to v3.6.1.

## [3.5.1] - 2026-03-10

### Fixed
- **TopN Effect scatter blank**: Column name mismatch (`log2_mean_A`/`log2_mean_B` vs actual `mean_a`/`mean_b`) caused the filter to return 0 rows, silently failing the `req()`. Fixed in both the scatter renderer and interpretation block.

### Added
- **Sub-tab info modals**: Added `?` help buttons to all 5 Run Comparator sub-tabs (Settings Diff, Protein Universe, Quantification, DE Concordance, AI Analysis) with detailed explanations of each visualization, color coding, and interpretation guidance.
- **Spectronaut parsing improvements**: Tree-character stripping for ExperimentSetupOverview, `[N]` column prefix handling, suffix-based fallback sample matching, per-sample peptide count aggregation via `rowMeans()`, gene column case-insensitive lookup, contextual TopN extraction, FASTA database deduplication, normalization chain (search_settings → AnalysisLog → setup_overview).

### Changed
- Quantification and Settings Diff sub-tabs now have scrollable wrappers with `overflow-y: auto`.
- TopN Effect section and Per-Sample QC section added to Quantification and Settings Diff sub-tabs respectively.
- App version bumped to v3.5.1.

## [3.5.0] - 2026-03-09

### Added
- **Run Comparator** (`R/server_comparator.R`): New module for comparing two analyses of the same dataset across tools. Three modes: DE-LIMP vs DE-LIMP (Mode A), DE-LIMP vs Spectronaut (Mode B), DE-LIMP vs FragPipe (Mode C, with or without FragPipe-Analyst DE stats).
  - **4-layer diagnostic pipeline**: Settings Diff (parameter comparison with highlighting), Protein Universe (overlap bar chart with summary cards), Quantification (scatter plot, per-sample correlation, bias density), DE Concordance (3x3 matrix, volcano overlay, discordant protein table).
  - **7-rule hypothesis engine**: Per-protein diagnostic explaining *why* each discordant protein disagrees. Tool-aware rules with context for Spectronaut and FragPipe structural differences. Categories: Direction reversal, Normalization offset, Variance estimation, Missing values, Peptide count, FC magnitude, Borderline.
  - **3x3 concordance matrix**: Classifies proteins as Up/Down/NS in each run for nuanced concordance analysis.
  - **Optional DIA-NN log upload** (Mode A): Upload DIA-NN log files to enrich Settings Diff with search-derived parameters (pg-level, proteoforms, library precursor count, pipeline step). Amber warning for library prediction logs; blue info banner for >1.2x library size mismatch.
  - **Optional MOFA2 decomposition**: Treats Run A and Run B as two views, decomposes joint variance. Variance heatmap, factor weights scatter, top weights table.
  - **Tool-aware Gemini prompt**: Includes structural differences between compared tools; DIA-NN library size context when log-derived.
  - **Claude ZIP export**: Settings diff, protein universe, DE results combined, discordant proteins with hypotheses, DIA-NN log parameters, comparison context, claude_prompt.md.
  - **Summary banner**: One-line overview with concordance rate, bias badge, dominant cause badge.
  - **Session persistence**: All comparator state saved/loaded with session .rds.

- **Search History** (`server_session.R`, `helpers_search.R`): Track all DIA-NN searches with 26-field parameter logging. CSV-based with file locking and shared volume support. Features:
  - **Expandable detail rows**: Click to view enzyme, mass accuracy, scan window, MBR, normalization, extra CLI flags, output dir, job ID.
  - **Import Settings button**: Apply search parameters from a previous search to the current search UI.
  - **Import Results button**: Load completed search results (report.parquet) directly from history, with SSH/SCP support for remote files. Auto-runs phospho detection and records to Analysis History.
  - **View Log button**: Display search_info.md from output directory (SSH or local).
  - **Cross-reference**: Link icon navigates between Search History and Analysis History via shared `output_dir`.

- **Analysis History & Projects** (`server_session.R`, `helpers_search.R`): Track pipeline runs with expandable detail rows, project assignment, and filtering.
  - **Projects JSON**: Group analyses by project name. Summary cards when filtered. `selectizeInput(create=TRUE)` for existing/new project names.
  - **DT expandable rows**: Click row to show full metadata. Action buttons (Info/Load/Assign) with `event.stopPropagation()`.

- **Chromatography QC**: TIC extraction from timsTOF .d files before search submission.
  - **Three plot views**: Faceted (per-run with median overlay), Overlay (normalized 0-1), Metrics (AUC bar chart + diagnostics table).
  - **Per-run diagnostics**: Shape deviation, RT shift, loading anomaly, file size outlier, late elution, elevated baseline, narrow gradient. MAD-based outlier detection.
  - **SSH mode**: SCP downloads `analysis.tdf` to temp, extracts locally.

- **DIA-NN log parser enhancements**: Extended `parse_diann_log()` with pg-level, proteoforms, library precursor count, lib/out-lib paths, pipeline step detection from SLURM job name.

- **FASTA database library**: Shared catalog with auto-upload to HPC when local-only paths detected. Fragment m/z range (`min_fr_mz`/`max_fr_mz`) recorded per entry. Path validation blocks HPC submission with local-only FASTA paths.

- **Smart partition selection**: Queries SLURM QOS limits (`sacctmgr show qos`) for per-user CPU limits (e.g., 64 CPUs on genome-center-grp/high). Auto-switches to publicgrp/low when user at capacity. Shows "Your usage: X/64 CPUs" in partition selector.

- Added `glue` and `data.table` to core package dependencies.

### Fixed
- **FASTA library local path bug**: Catalog entries stored macOS-local paths (`/Users/...`) in `remote_dir`, causing Apptainer bind mount failures on HPC. Now validates remote paths and auto-uploads via SCP.
- **Partition selector "no limit info"**: `sacctmgr show assoc` returned empty limits; limits are on QOS objects, not associations. Changed to `sacctmgr show qos where name={account}-{partition}-qos`.
- **Array progress inflated counts**: `sacct` counted parent job + `.extern`/`.batch` substeps. Now filters to `JOBID_N` format entries only.
- **`parse_diann_log` fr_mz/pr_charge**: Fragment m/z and precursor charge flags routed to `extra_cli_flags` instead of `params`. Now parsed via `value_map`.
- **Docker container name with special characters**: Sanitized with `gsub("[^a-zA-Z0-9_.-]", "_", ...)`.
- **`max_pr_mz` default was 1200**: Changed to DIA-NN's actual default of 1800.
- **Parallel search OOM on timsTOF**: Default `mem_per_file` bumped from 32 GB to 64 GB.
- **TIC extraction auto-triggered**: `observeEvent(list(btn, trigger))` fired on button render. Fixed with separate `reactiveVal` trigger pattern.
- **Removed default raw data path**: SSH raw data directory input no longer pre-filled with a test path.

### Changed
- All comparator visualizations use `plotlyOutput`/`renderPlotly` (bslib sub-tab safety).
- Protein universe uses plotly stacked bar instead of ComplexUpset (simpler for 2-set comparison).
- App version bumped to v3.5.0.

## [3.3.0] - 2026-03-06

*Note: v3.3.0 and v3.4.0 were development milestones rolled into v3.5.0.*

### Added
- **Chromatography QC**: TIC extraction from timsTOF .d files, run diagnostics, instrument metadata export.
- **Run Comparator**: Initial implementation (completed in v3.5.0).
- **Default parallel memory**: Bumped to 64 GB (timsTOF OOM fix).

## [3.2.1] - 2026-03-05

### Added
- **Search history**: Track all DIA-NN searches with full parameter logging (enzyme, mass accuracy, scan window, MBR, normalization, extra CLI flags). CSV-based with file locking, shared volume support. Expandable detail rows, Import Settings button, View Log button, cross-reference to Analysis History via `output_dir`.
- **DIA-NN log parser** (`parse_diann_log`): Extract search parameters from DIA-NN log files — version, FASTA, enzyme, mass accuracy, scan window, MBR, modifications, fragment m/z range, precursor charge range. Inverse of `build_diann_flags()`. 107 unit tests.
- **Claude export enhancements**: Export for Claude now includes DIA-NN search settings (version, mass accuracy, modifications, scan window, etc.) with prompt for publication-ready Methods section; per-group missingness summary; MOFA2 variance explained per factor/view; phosphosite DE summary with top sites; covariate metadata per sample.
- **Parallel job consistency check**: Validates step dependency chain integrity before monitoring parallel search jobs.
- **Search history unit tests**: 51 tests covering `record_search()`, `update_search_status()`, `search_history_path()`, `backfill_search_history()`.

### Fixed
- **Array progress inflated counts**: `sacct` for parallel array jobs counted parent job and `.extern`/`.batch` substeps, inflating progress (e.g., 51/41 instead of 37/41). Now filters to only array task entries (`JOBID_N` format).
- **`sacct` `.extern` step false COMPLETED**: `check_slurm_status()` now uses `--format=JobID,State` and filters out `.extern`/`.batch` substep lines that report COMPLETED even when the main job is PENDING/FAILED.
- **`parse_diann_log` fr_mz/pr_charge**: `--max-fr-mz`, `--min-fr-mz`, `--min-pr-charge`, `--max-pr-charge` were incorrectly routed to `extra_cli_flags`. Now parsed via `value_map` and flow properly into `search_params`.
- **Docker container name with special characters**: `analysis_name` is now sanitized via `gsub("[^a-zA-Z0-9_.-]", "_", ...)` before building container name.
- **`max_pr_mz` default wrong**: UI and all fallbacks used 1200 instead of DIA-NN's actual default of 1800. FASTA library entries and searches recorded incorrect precursor m/z range when Advanced Options wasn't opened.

### Changed
- **AI Summary export buttons**: "Export Report" renamed to "Download as HTML"; both "Download as HTML" and "Export for Claude" are hidden until AI summary is generated (progressive reveal via `shinyjs`).
- **Analysis name field**: No longer has a default value — users must provide a name before submitting a DIA-NN search.
- **HF Space**: Search History and Analysis History tabs hidden on Hugging Face deployment (not useful in ephemeral container).
- **DIA-NN search settings in Analysis_Parameters.txt**: Expanded from 7 fields to ~20 (version, mass accuracy MS1/MS2, mode, scan window, variable mods, fragment m/z range, precursor charge range, RT profiling, normalization, extra CLI flags).

## [3.1.1] - 2026-02-26

### Fixed
- **Volcano plot significance mismatch**: Y-axis used raw P.Value but coloring used adj.P.Val, causing ~746 proteins to appear above the threshold line while colored as "Not Sig". Now the horizontal dashed line is drawn at the raw P.Value corresponding to adj.P.Val = 0.05, so the line and coloring agree visually.
- **Volcano significance coloring**: Removed logFC cutoff from significance determination — all proteins with adj.P.Val < 0.05 are now colored red. logFC vertical lines remain as visual guides.
- **Default logFC cutoff**: Changed from 1.0 (2-fold change) to 0.6 (~1.5-fold change).
- **CV Analysis cards compressed**: Plotly annotation cards above the scatter plot kept overlapping and rendering incorrectly in bslib sub-tabs. Replaced with a ggplot subtitle showing per-group median CV and % proteins below 20% CV.
- **CV Analysis tab layout**: Wrapped in scrollable div with min-height to prevent bslib from compressing the scatter plot. Content stacks vertically with overflow-y scroll.

### Added
- **Volcano DE protein count**: Info box now shows "78 DE proteins (X up, Y down)" so users can immediately see the count.
- **Export Data panel**: New "Export Data" tab in the Output dropdown with prominent buttons for downloading Results CSV and CV Analysis CSV.
- **AI Summary HTML export**: "Export Report" button on AI Summary tab generates a styled standalone HTML report with gradient header, formatted tables, and print-friendly CSS.
- **Docker update scripts**: `update_docker.sh` (bash) and `update_docker.ps1` (PowerShell) for one-command pull + rebuild on Windows/Mac/Linux.

## [3.1.0] - 2026-02-24

### Added
- **UI Overhaul**:
  - `page_navbar()` layout with dark navbar (`#2c3e50`), hover-activated dropdown menus with smooth animations
  - Dropdown section labels (Setup / Results / AI) in the Analysis menu via JS injection
  - Active tab teal underline indicator (`--flatly-success`)
  - `nav_spacer()` + native `nav_item()` for gear icon (replaces JS injection)

- **Accordion Sidebar**:
  - Three collapsible `accordion_panel()` sections: Upload Data (open by default), Pipeline Settings, AI Chat
  - Phospho and XIC sections in separate conditional accordions outside the main group
  - Session buttons and core facility controls outside accordion for always-visible access

- **DE Dashboard Sub-tabs**:
  - Replaced grid+accordion layout with `navset_card_tab(id = "de_dashboard_subtabs")`
  - Four sub-tabs: Volcano (with heatmap below), Results Table, PCA, Robust Changes
  - PCA moved from Data Overview into DE Dashboard

- **Core Facility Mode**:
  - Activated by `DELIMP_CORE_DIR` env var pointing to directory with `staff.yml`
  - SQLite database (`delimp.db`) with 4 tables: searches, qc_runs, reports, templates
  - WAL mode for concurrent SQLite access
  - Staff YAML config auto-fills SSH host, username, key path, SLURM account/partition
  - Search DB tab: full-width job history with 6 filters (text, lab, status, staff, instrument, LC method)
  - Instrument QC dashboard: protein/precursor/TIC trend plots with ±2SD control lines, instrument filter, date range
  - Quarto report generation: standalone HTML with metadata, QC bracket, volcanos, DE stats, top proteins
  - Template system for saving/loading search presets
  - New files: `R/helpers_facility.R`, `R/server_facility.R`, `report_template.qmd`, `seed_test_db.R`

### Changed
- `R/ui.R`: Rewritten outer wrapper from `page_sidebar()` + `navset_card_tab()` to `page_navbar()` with direct nav items
- `R/ui.R`: Sidebar width reduced from 320 to 300px
- `R/ui.R`: Data Overview reduced from 6 to 5 sub-tabs (PCA moved to DE Dashboard)
- `R/ui.R`: Removed `.de-dashboard-grid` CSS class and accordion heatmap wrapper
- `CLAUDE.md`: Condensed from ~500 to ~130 lines; detailed patterns moved to `docs/PATTERNS.md`, TODOs to `docs/TODO.md`
- `README_GITHUB.md`, `README_HF.md`, `USER_GUIDE.md`: Updated for v3.1 features
- App version bumped to v3.1

### Fixed
- **Navbar text invisible on dark background**: Flatly theme renders dark text on `navbar-inverse`. Fixed with CSS `!important` overrides for white text on `.navbar .nav-link` and `.navbar-brand`.
- **Hidden tab fragments visible**: `nav_hide()` on `page_navbar` leaves letter fragments. Fixed with `width: 0 !important; overflow: hidden !important` on hidden nav items.
- **bslib deprecation warning**: `page_navbar(bg=...)` deprecated in bslib 0.9.0+. Changed to `navbar_options = navbar_options(bg = "#2c3e50")`.

## [3.0.0] - 2026-02-20

### Added
- **Multi-Omics MOFA2 Integration**:
  - New **Multi-Omics MOFA2** tab for unsupervised integration of 2-6 data views using MOFA2
  - Dynamic view cards: add/remove views, smart RDS parser (DE-LIMP sessions, limma EList/MArrayLM, matrices, data frames)
  - CSV/TSV/Parquet matrix upload with auto-log2 detection
  - Phospho tab integration as data source
  - Sample matching with overlap statistics and color-coded status
  - MOFA training via `callr::r()` subprocess (isolates basilisk/Python from Shiny's event loop)
  - 5 results tabs: Variance Explained heatmap, Factor Weights browser (plotly), Sample Scores scatter, Top Features table (DT), Factor-DE Correlation
  - **Mouse Brain example dataset** (2-view): proteomics + phosphoproteomics, 16 samples
  - **TCGA Breast Cancer example dataset** (3-view): mRNA + miRNA + protein, 150 samples, 3 subtypes
  - Session save/load, methodology text, reproducibility logging
  - New packages: `MOFA2`, `basilisk`, `callr`

- **DIA-NN Docker Local Backend**:
  - `Dockerfile.search`: Multi-stage build embeds DIA-NN binaries from user's pre-built `diann:2.0` image into DE-LIMP container
  - `docker-compose.yml`: One-command deployment (`docker compose up`) with data volume mounts
  - DIA-NN runs inside the container as a background process — no Docker-in-Docker required
  - Windows Docker Desktop support: intermediate files written to container-internal `/tmp` to avoid FUSE layer issues with large speclib writes
  - `build_diann_docker.sh` / `build_diann_docker.ps1`: User-facing scripts to build DIA-NN image (license compliance)

- **Windows Docker Deployment**:
  - `WINDOWS_DOCKER_INSTALL.md`: Step-by-step guide for Windows users (Docker Desktop + WSL2)
  - Zero R installation required — everything runs in Docker
  - `MACOS_DOCKER_INSTALL.md` and `LINUX_DOCKER_INSTALL.md`: Platform-specific guides

### Changed
- `R/ui.R`: Renamed "Multi-View Integration" tab to "Multi-Omics MOFA2"
- `R/helpers_search.R`: Docker command uses `--entrypoint sh -c` wrapper to write intermediate files to `/tmp/diann_work` then copy final outputs to mounted volume
- `Dockerfile` and `Dockerfile.search`: Added MOFA2, basilisk, callr packages; pre-initialize basilisk Python env at build time
- `docs/index.html`: Comprehensive update for v3.0 features, deployment options, corrected Getting Started
- `README_GITHUB.md`, `README_HF.md`, `USER_GUIDE.md`: Added MOFA2 documentation, fixed Quick Start
- App version bumped to v3.0

### Fixed
- **R session crash during MOFA2 training**: basilisk's Python subprocess conflicts with Shiny's httpuv event loop. Fixed by running entire MOFA pipeline in isolated subprocess via `callr::r()`
- **Variance Explained heatmap empty**: `r2_per_factor` is a matrix (factors × views), not a nested named list. Fixed iteration logic.
- **Factor-DE Correlation tab blank**: Controls required `values$fit` which is NULL when using test data. Now shows info message when no DE results available.
- **Docker DIA-NN "Could not save" on Windows**: Large intermediate `.predicted.speclib` files fail on Windows Docker FUSE layer. Fixed by writing to container-internal `/tmp` then copying final outputs.
- **Variance explained download handler**: Same matrix iteration bug as the heatmap render.

## [2.5.0] - 2026-02-18

### Added
- **GSEA Expansion — Multi-Database Enrichment**:
  - **Four enrichment databases**: GO Biological Process (BP), GO Molecular Function (MF), GO Cellular Component (CC), and KEGG Pathways
  - **Ontology selector**: Dropdown to switch between BP/MF/CC/KEGG on the GSEA tab
  - **Per-ontology caching**: Results cached per database; switching back loads instantly without re-computation
  - **Contrast indicator**: Shows active contrast with stale-results warning when contrast changes
  - **UniProt API organism detection**: Queries `rest.uniprot.org` to determine organism from accession when no suffix is present — works automatically for human, mouse, rat, and 9 other species
  - **Robust ID mapping**: Handles multiple protein ID formats (pipe-separated, isoform suffixes, organism suffixes); fallback from UNIPROT to SYMBOL ID types
  - **KEGG organism mapping**: Supports 11 species with automatic organism code detection
  - **Dynamic plot titles**: All GSEA plots show which database was used
  - Updated info modals with all 4 database descriptions
  - Session save/load for GSEA cache, last contrast, and organism DB

- **AI Summary — All Comparisons Analysis**:
  - AI Summary now analyzes **all contrasts** simultaneously, not just the selected one
  - Cross-comparison biomarker detection: identifies proteins significant in ≥2 comparisons
  - Enhanced prompt with 5 structured sections: Overview, Key Findings Per Comparison, Cross-Comparison Biomarkers, High-Confidence Biomarker Insights, Biological Interpretation
  - Adaptive token budget: top 30/20/10 proteins per contrast (scales with number of contrasts)
  - `?` info modal explaining what data is/isn't sent to Gemini API

- **MDS Plot Coloring**:
  - Color MDS plot by Group, Batch, Covariate1, or Covariate2
  - Colorblind-friendly Okabe-Ito palette
  - Dynamic dropdown updates when metadata changes

- **Complete Dataset Export**:
  - Download button on Dataset Summary tab
  - Exports: protein IDs, gene symbols, DE stats for ALL contrasts (suffixed columns), per-sample expression values, metadata as header comment rows

- **Phosphoproteomics Phase 2 — Kinase Activity & Motif Analysis**:
  - **KSEA kinase-substrate enrichment analysis**: Infers upstream kinase activity from phosphosite fold-changes using `KSEAapp` CRAN package with PhosphoSitePlus + NetworKIN database
  - **KSEA bar plot**: Horizontal bar chart of kinase z-scores (top 15 activated + top 15 inhibited), colored by direction, with substrate count annotations
  - **KSEA results table**: Filterable, sortable DT datatable with kinase gene, z-score, FDR, substrate count. Downloadable as CSV.
  - **Sequence logo motif analysis**: Displays amino acid enrichment around regulated phosphosites using `ggseqlogo`. Separate logos for up-regulated and down-regulated sites.
  - **FASTA upload**: Sidebar file input for protein FASTA. Parses UniProt-format headers to extract accessions. Enables accurate flanking sequence extraction for motif analysis.
  - New "Kinase Activity" and "Motif Analysis" tabs in phospho results navset
  - New packages: `KSEAapp` (CRAN), `ggseqlogo` (CRAN)

- **Phosphoproteomics Phase 3 — Advanced Features**:
  - **Protein-level abundance correction**: Checkbox to subtract protein-level logFC from phosphosite logFC, isolating phosphorylation stoichiometry changes
  - **AI context integration**: Phosphosite DE results and KSEA kinase activities appended to Data Chat context when phospho analysis is active
  - **Session persistence**: KSEA results, FASTA sequences, and all Phase 2/3 state saved/loaded with sessions

### Changed
- `R/server_gsea.R`: Rewritten from 144 to ~400 lines (multi-DB, caching, organism detection)
- `R/server_ai.R`: AI Summary rewritten for all-contrast analysis; send_chat scaled for large datasets
- `R/helpers_phospho.R`: Extended from 210 to ~380 lines (5 new helper functions)
- `R/server_phospho.R`: Extended from 650 to ~950 lines (KSEA, motifs, protein correction)
- `Dockerfile`: Added `KSEAapp` and `ggseqlogo` CRAN package installation
- App version bumped to v2.5

### Fixed
- **Export CSV crash**: "Column name `Protein.Group` must not be duplicated" when limpa's topTable includes Protein.Group column
- **Gemini token limit**: Scale data sent to AI based on sample count; group-level Mean/SD for >100 samples
- **P-value histogram y-axis**: Cap y-axis to show distribution shape when first bin dominates; annotate clipped bin count
- **P-value dropdown clipping**: Removed card_body wrapper and added z-index stacking to prevent plot from overlapping dropdown
- **Comparison dropdown width**: Full-width dropdowns on all comparison banners; buttons moved inline

### Planned — Future
- PhosR integration (RUVphospho normalization, kinase profiling, signalome)
- FASTA-based protein-relative position mapping for Path B sites

## [2.4.0] - 2026-02-17

### Added
- **Phosphoproteomics Tab (Phase 1)**: Site-level differential phosphorylation analysis
  - **Auto-detection**: Scans `Modified.Sequence` for `UniMod:21` on file upload; shows blue banner in Data Overview with "Open Phospho Tab" button
  - **Two input paths**:
    - Path A (recommended): Upload DIA-NN 1.9+ `site_matrix_0.9.parquet` or `site_matrix_0.99.parquet`
    - Path B: Parse phosphosites directly from `report.parquet` with configurable localization confidence threshold (0.5–1.0)
  - **Site extraction algorithm** (Path B): Character-by-character Modified.Sequence parser to locate `(UniMod:21)` positions, expand multiply-phosphorylated peptides, aggregate per SiteID × Run via max intensity ("Top 1" method)
  - **Site-level limma DE**: Filter sites (≥2 non-NA per group), tail-based imputation (Perseus-style: mean − 1.8 SD, width 0.3 SD), optional normalization (none/median/quantile), standard limma pipeline
  - **Phospho Volcano**: ggplot2 volcano with ggrepel labels formatted as "Gene Residue+Position" (e.g., "MAPK1 T185"). Colors: Significant=#E63946, FDR-only=#457B9D, NS=gray70. Downloadable as PDF.
  - **Site Table**: DT datatable with SiteID, Gene, Residue, Position, logFC, adj.P.Val, localization confidence. Filterable, sortable, downloadable as CSV.
  - **Residue Distribution**: Grouped bar chart (S/T/Y) comparing "All quantified" vs "Significant". Subtitle with expected ~85% Ser / ~14% Thr / ~1% Tyr.
  - **QC: Completeness**: Histogram of per-site % samples quantified with red dashed line at 50% threshold.
  - **Sidebar controls**: Conditional on phospho detection — input mode, localization slider, normalization radio, "Run Phosphosite Analysis" button
  - **Normalization warning**: Yellow alert for phospho-enriched data explaining DIA-NN normalization assumptions
  - **Educational expandable**: Explains site-level vs protein-level analysis, localization confidence, imputation approach
  - **Session save/load**: All phospho state persisted and restored
  - **Reproducibility logging**: Pipeline steps logged with parameters
- New files: `R/helpers_phospho.R` (210 lines), `R/server_phospho.R` (650 lines)

## [2.3.0] - 2026-02-17

### Changed
- **Modularization**: Split 5,139-line monolith into `app.R` orchestrator + 12 `R/` module files
- **Upload limit**: Increased from 500 MB to 5 GB
- **Dockerfile**: Updated for directory-based `runApp()` and `COPY R/` directive

## [2.2.0] - 2026-02-17

### Added
- **Contextual Help System**: 15 info modal buttons (`?`) across all major tabs providing in-context guidance
  - QC Plots: Normalization Diagnostic, DPC Fit, MDS Plot, Group Distribution, P-value Distribution
  - Data Overview: Signal Distribution, Expression Grid
  - DE Dashboard: Volcano/table interaction guide with threshold explanation
  - Consistent DE: High-Consistency Table (%CV explained), CV Distribution
  - QC Trends: Metric definitions, sort order explanation, drift detection tips
  - Gene Set Enrichment: GSEA overview + Results Table column definitions (NES, p.adjust, etc.)
  - Reproducibility: Methodology (LIMPA pipeline, limma/eBayes, covariates)
  - Data Chat: Privacy info, API key instructions, selection integration
- **Volcano Plot → Results Table Filtering**: Selecting proteins in the volcano plot now filters the results table to show only selected proteins (bidirectional sync)
- **MDS Plot Legend**: Visible legend at bottom-right with white background box (was being clipped off-screen)
- **Heatmap Expanded by Default**: DE Dashboard heatmap accordion now opens expanded

### Changed
- **Normalization Diagnostic**: "What am I looking at?" moved from in-page expandable to modal dialog (no longer interferes with plot)
- **P-value Distribution**: "How do I interpret this?" moved from in-page expandable to modal dialog; guidance banner moved below plot (no longer overlaps comparison dropdown)

### Fixed
- **Bad Merge Recovery**: Restored 333 lines lost in merge commit `1361b62` including MS2 Intensity Alignment and XIC auto-load features
- **XIC Viewer Facet Error**: Added guard for "Faceting variables must have at least one value" warning when modal first opens
- **DE Table Row Index Mapping**: Fixed row selection observer to correctly index into filtered data when volcano selection is active

## [2.1.1] - 2026-02-16

### Added
- **XIC Chromatogram Viewer**: On-demand fragment-level chromatogram inspection for differentially expressed proteins
  - Sidebar section "5. XIC Viewer" with directory path input and load button
  - XIC buttons on DE Dashboard results table and Grid View modal
  - Three display modes: Facet by sample, Facet by fragment, Intensity alignment
  - Split-axis MS1/MS2 view with independent y-axes (MS1 top, fragments bottom)
- **MS2 Intensity Alignment**: Spectronaut-style stacked bar chart for fragment ion ratio consistency
  - Each bar = one sample, colored segments = relative fragment proportions
  - Automatic inconsistency detection: flags samples with deviation > mean + 2×SD
  - Green/amber guidance banners with sample IDs and possible causes
  - Bars ordered by experimental group with dashed separators
  - Cosine similarity and deviation scores in tooltips
  - Precursor selector, group filter, and MS1 toggle controls
  - Prev/Next protein navigation through significant DE proteins
  - Download handler for PNG export (14×10 inch, 150 DPI)
  - Info panel with protein stats, RT range, and DE statistics
- **DIA-NN 2.x Format Support**: Auto-detects and handles both DIA-NN 1.x (wide) and 2.x (long) XIC formats
- **Ion Mobility / Mobilogram Support**: Detects mobilogram files, checks for non-zero data (timsTOF/PASEF only)
  - Blue gradient toggle with bolt icon; prominent banner when IM mode is active
- **XIC Directory Auto-Population**: Auto-detects `_xic` sibling directory when data is loaded
  - Smart path resolution: accepts `.parquet` file paths or directories without `_xic` suffix
- **Precursor Map from In-Memory Data**: Builds protein→precursor mapping from `values$raw_data` (no file I/O)

- **XIC Auto-Load**: XICs automatically load when `_xic` directory is detected on data upload (no manual button click needed)

### Fixed
- **Assign Groups Layout**: Fixed Run Pipeline button pushed off-screen on MacBook (CSS Grid → Flexbox)
- **Arrow/dplyr Conflicts**: `arrow::select` masking `dplyr::select` — use explicit `dplyr::select()` in XIC code
- **Tidy Evaluation Issues**: `rlang::sym("pr")` and `rename(Precursor.Id = pr)` replaced with base R equivalents
- **Mobilogram Detection Pattern**: Fixed file pattern to match DIA-NN naming convention (`_mobilogram.parquet` not `.mobilogram.parquet`)
- **Windows Path Compatibility**: Fixed regex in XIC smart path resolution to use `basename()` instead of forward-slash pattern
- **HF XIC Visibility**: XIC Viewer sidebar, buttons, and auto-load logic hidden on Hugging Face Spaces (detected via `SPACE_ID` env var); replaced with info note linking to GitHub for local download

## [2.1.0] - 2026-02-13

### Added
- **Four-Way Comparison Selector Synchronization**: Signal Distribution, Expression Grid, P-value Distribution, and DE Dashboard selectors now sync automatically when any one changes
- **P-value Distribution Diagnostic Plot**: New QC Plots sub-tab with automated pattern detection
  - Color-coded guidance banners (healthy, inflation, low power, model issues)
  - Expandable interpretation guide
  - Actionable recommendations for troubleshooting
- **CV Distribution Histogram**: New Consistent DE sub-tab showing protein variability distribution per experimental group
- **Dataset Summary DE Counts**: Shows differential expression protein counts for all comparisons with explicit directional language ("X proteins higher in GroupA")
- **AI Summary Sub-Tab**: Dedicated tab in Data Overview (moved from DE Dashboard button)
  - Inline display instead of modal popup
  - Cleaner workflow for AI-generated summaries
- **Volcano Plot Annotations**:
  - Colored threshold lines (blue FDR, orange logFC)
  - Significance criteria legend box in upper-left corner
  - Uses plotly native annotations for reliable rendering
- **Comparison Selectors**: Purple gradient banners on Signal Distribution, Expression Grid, and P-value Distribution tabs

### Changed
- **Signal Distribution Plot**: Now always shows DE coloring when results are available (removed manual toggle buttons)
- **Signal Distribution Plot**: Uses dedicated `contrast_selector_signal` instead of main selector
- **Responsive Plot Heights**: All plots now use viewport-relative units (`vh`, `calc()`) for optimal viewing on any screen size
- **Assign Groups Interface**: Converted from modal popup to permanent sub-tab in Data Overview
- **DE Dashboard**: Removed AI Summary button (moved to Data Overview tab)

### Fixed
- **Volcano Plot Legend**: Fixed text rendering issues by switching from ggplot annotations to plotly native system
- **Reactive Loops**: Signal Distribution table no longer filters based on selection-derived reactives

### Technical
- All comparison selectors (`contrast_selector`, `contrast_selector_signal`, `contrast_selector_grid`, `contrast_selector_pvalue`) sync bidirectionally
- Plotly annotations use paper coordinates (`xref = "paper", yref = "paper"`) for absolute positioning
- No new R package dependencies (all features use existing packages)

## [2.0.0] - 2025-XX-XX

### Initial Release
- Interactive differential expression analysis for DIA-NN proteomics data
- Limpa pipeline integration (DPC-CN normalization, DPC-Quant protein quantification, limma statistics)
- Google Gemini AI chat integration
- Session save/load functionality
- GSEA integration (GO Biological Process)
- Reproducibility code logging
- Example data with Affinisep vs Evosep comparison
- Education tab with embedded resources
