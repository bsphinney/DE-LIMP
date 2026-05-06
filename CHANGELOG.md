# Changelog

All notable changes to DE-LIMP will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
