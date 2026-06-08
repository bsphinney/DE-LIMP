# DE-LIMP Gotchas Reference

Quick-reference tables for known issues and their solutions. **This is the single source of truth for gotchas** — CLAUDE.md keeps only a category index that points here. When you hit (or fix) a non-obvious bug, add a row to the relevant table below.

## R Shiny / bslib

| Problem | Solution |
|---------|----------|
| Navbar text invisible on dark bg | Flatly theme CSS override: `.navbar .nav-link { color: rgba(255,255,255,0.75) !important; }` |
| Hidden tabs show letter fragments | `.navbar .nav-item[style*='display: none'] { width: 0 !important; overflow: hidden !important; }` |
| `page_navbar(bg=...)` deprecation | Use `navbar_options = navbar_options(bg = ...)` (bslib 0.9.0+) |
| `source()` doesn't start app | Use `shiny::runApp()` instead |
| Server-code edits don't take effect after restarting the app | `runApp()` in the **same R session** can reuse cached function definitions — the VERSION banner updates (read fresh in `app.R`) while `R/server_*.R` functions stay stale, so changes look "not applied." Fully **restart the R session** (RStudio: Session → Restart R / ⌘⇧F10), then `runApp()`. Cost us ~3 confused rounds in v3.11.38–40. |
| Selections disappear after clicking | Reactive loop — table must not depend on selection-derived reactives |
| bslib `card()` doesn't render | Use plain `div()` for top-level nav_panel content |
| `uiOutput` vanishes in `navset_card_tab` | Use static HTML + `shinyjs::html("div_id", content)`. `plotlyOutput` with `req()` is safe. |
| `renderPlot` crashes in hidden sub-tab | `invalid quartz() device size` on macOS (0-width container). Use `plotlyOutput`/`renderPlotly`. |
| `return()` inside `withProgress` | Exits `withProgress` not enclosing function. Use flat `tryCatch`. |
| `<<-` inside `withProgress` fails | `withProgress` uses `eval(substitute(expr), env)` — `<<-` can't find parent vars. Use `new.env()` + `<-` instead. |
| Shiny hidden input not registered by JS | Use `div(style="display:none;", radioButtons(...))` for `conditionalPanel` |
| `observeEvent(list(btn, trigger))` auto-fires | Fires when button first renders (NULL→0). Use a separate `reactiveVal` trigger with `ignoreInit = TRUE`. |
| Quarto `output_file` path error | Pass filename only, then `file.rename()` to target dir |
| SQLite `Parameter N does not have length 1` | Use `NA_character_` instead of `NULL` |
| History tab slow with network CSV | Multiple `activity_log_read()` per render. Use `cached_activity_log()` reactive (read once per invalidation). |

## DIA-NN

| Problem | Solution |
|---------|----------|
| DIA-NN `Genes` column has accessions | Not gene symbols (e.g. `A0A075B6K5;P80748`). Validate with length/pattern. Real genes from `bitr()` UNIPROT → SYMBOL. |
| `readDIANN` data.table column error | Must pass `format="parquet"` for .parquet files |
| DIA-NN empirical lib is `.parquet` not `.speclib` | DIA-NN 2.0+ saves empirical libs in parquet. Use `empirical.parquet` in `--lib`/`--out-lib`. Predicted libs stay `.predicted.speclib`. |
| `--quant-ori-names` required on ALL steps | Per Vadim: preserves original filenames in `.quant` files. Without it, container bind-mount path differences cause naming mismatches between steps. |
| `--fasta-search`/`--predictor` Step 1 only | Including in Steps 2–5 causes full FASTA re-digest. `generate_parallel_scripts()` strips these from `step_flags`. |
| Auto mass acc + `--use-quant` | Produces different results. `generate_parallel_scripts()` forces `mass_acc_mode = "manual"`. Full flag reference in PATTERNS.md. |
| `max_pr_mz` default was 1200 not 1800 | DIA-NN default for `--max-pr-mz` is 1800. UI + all fallbacks were wrong, causing wrong range when Advanced Options not opened. |
| Parallel search OOM on timsTOF | Default `mem_per_file` was 32 GB, insufficient for DIA-PASEF. Now 64 GB. |
| Log import ignores `fr_mz`/`pr_charge` | `parse_diann_log` now parses `--max-fr-mz`/`--min-fr-mz`/`--min-pr-charge`/`--max-pr-charge` via `value_map` into `search_params` (was dumping to `extra_cli_flags`). |
| Two DIA-NN containers, only one reads `.raw` | `/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif` has .NET, reads Thermo `.raw`. `/quobyte/proteomics-grp/apptainers/diann2.3.0.sif` has NO .NET — `.raw` silently skipped. Always use the `dia-nn/` version unless only `.d`/`.mzML`. |
| DIA-NN binary path inside container | `/diann-2.3.0/diann-linux`, NOT just `diann`. `apptainer exec image.sif /diann-2.3.0/diann-linux ...`. |

## Data & Columns

| Problem | Solution |
|---------|----------|
| `nrow(raw_data$E)` counts precursors not proteins | Use `length(unique(raw_data$genes$Protein.Group))` for protein groups (~3k vs ~40k). `y_protein$E` rows ARE protein groups (post-pipeline). |
| `y_protein` `colSums` error | It's a limma `EList`. Extract `$E` for the expression matrix. |
| `arrow::select` masks `dplyr::select` | Use `dplyr::select()` explicitly |
| R regex `\\s` invalid | Use `[:space:]` in base R regex (POSIX ERE) |
| `unlist()` on nested lists causes row mismatch | Nested elements expand to multiple rows. Use `vapply(x, function(v) paste(v, collapse="; "), character(1))`. |
| Character matrix subsetting fails on Linux | Empty-string rownames break `mat[ids, ]` on Linux (works on macOS). Use numeric indices via `match(ids, rownames(mat))`. |
| Volcano P.Value vs adj.P.Val mismatch | Y-axis uses raw `P.Value` for spread; dashed line at `max(P.Value)` among `adj.P.Val < 0.05` proteins. |
| MOFA2 views need same sample names | Subset to matched pairs, assign common labels (`Sample_1`, `Sample_2`, ...). |
| Contaminant proteins have `Cont_` prefix | After `--fasta` contaminant lib, DIA-NN prefixes IDs with `Cont_`. Detect via `grepl("^Cont_", Protein.Group)`. Link to NCBI Protein (not UniProt). |
| No-replicates mode skips DE | Groups with <2 reps get quant only. `values$fit` stays NULL. Expression Grid, PCA, Signal Distribution still work; Volcano/DE/GSEA need replicates. |
| NCBI RefSeq accessions need gene mapping | `Genes` column empty/accession-only for NCBI FASTA. Use batch E-utilities (`esummary` on protein UIDs). Cache as TSV alongside FASTA. |
| FASTA library `remote_dir` stored local paths | `fasta_library_file_paths()` validates remote paths; auto-uploads via SCP if local-only. Blocks HPC submission with local-only FASTA paths. |
| Older TDF missing `SummedIntensities` | `extract_tic_timstof()` auto-detects intensity column: `SummedIntensities` → `AccumulatedIntensity` → `MaxIntensity` → any `*ntensit*`. |

## SSH / HPC / SLURM

| Problem | Solution |
|---------|----------|
| Symlinks in container bind mounts don't resolve | If you bind `selected_raw:/work/data` and it has symlinks to `../8min/file.raw`, the container can't follow them (target not mounted). Bind the parent dir (`--bind /base:/work`) or `cp` instead of `ln -sf`. |
| SSH output encoding crash | `iconv(..., sub="")` in `ssh_exec`/`scp_download`/`scp_upload` |
| SSH rapid connections rejected (255) | HPC `MaxStartups` throttling. Batch operations into fewer SSH calls; use ControlMaster multiplexing. |
| macOS SSH ControlPath too long | Unix domain sockets limited to 104 bytes on macOS. Use `/tmp/.delimp_<user>_<host>` (R's `tempdir()` is ~105 chars). |
| SSH file browser Unix socket path | File browser modal creates SSH connections — same ControlMaster socket-length constraint applies. |
| `parse_sbatch_output` returns dirty ID | SSH stdout may have trailing `\r`/whitespace. Always `trimws()` parsed job IDs (else `--dependency=afterok:12345\r` silently fails). |
| SSH auto-connect blocks event loop | Connection test on startup takes 10–30s. Run via `later::later()` or fast-fail timeout. Stale sockets detected with `ssh -O check`. |
| `nohup bash ... &` over SSH won't return | Remote bash inherits SSH stdin; even with `nohup`, SSH waits for it to close. Add `< /dev/null` between `bash` and `&`. |
| SLURM limits on QOS not associations | `sacctmgr show assoc` returns empty limits. Use `sacctmgr show qos where name={account}-{partition}-qos` for `GrpTRES`/`MaxTRESPU`. |
| Per-user CPU limit (not account) is binding | `MaxTRESPU` (e.g. 64 CPUs) constrains individual users. `GrpTRES` (e.g. 616) is shared. `select_best_partition()` uses the per-user limit. |
| `sacct` `.extern` step falsely reports COMPLETED | `.extern`/`.batch` substeps COMPLETE even when main job is PENDING/FAILED. Request `JobID,State` and filter out lines containing `.`. |
| Array progress sacct inflated counts | `sacct -j ARRAY_ID` returns parent + substeps per task. Filter to `JOBID_N` only: `grepl("_", jid) && !grepl("\\.", jid)`. |
| Partial retry dependency chain | After retrying failed step-2 tasks, `scontrol update` step 3's dependency to the retry job ID, else step 3 starts before retries finish. |
| SLURM proxy inside Apptainer | SLURM commands unavailable inside container. Proxy process outside relays via temp files. All 9 command paths covered. |
| Paths with spaces in SLURM scripts | Quote everything: `#SBATCH -o "path"`, `--bind "path":/work`. Launcher uses `q()` helper. |
| Python stdout buffered in SLURM | Output not flushed to `.out` until process ends. `sys.stdout.reconfigure(line_buffering=True)`. |
| Docker container name rejected | Spaces/special chars fail Docker naming `[a-zA-Z0-9][a-zA-Z0-9_.-]*`. Sanitize with `gsub("[^a-zA-Z0-9_.-]", "_", name)`. |
| Docker SSH key permissions | Windows bind mounts lose Unix perms. `Launch_DE-LIMP_Docker.bat` copies keys to internal volume with `chmod 600`. |
| Load from HPC needs build-time guard | `conditionalPanel` alone insufficient on HF (button flashes before JS hides it). Wrap in `if (!is_hf_space)` in `build_ui()`. |
| **NEVER use mounted drives for app state** | SMB mounts (`/Volumes/proteomics-grp/`) may vanish. App state (activity log, cluster usage, lab members) MUST use local `~/.delimp_*`. Cross-user sharing via SSH/SCP. |
| **Derived data stays with source data** | Cached/computed data (TIC cache, `session.rds`, `search_info.md`) lives in the raw-data or output dir, NOT home dir — so any lab member scanning the dir gets it, and `delete data = delete cache`. Only app-level config goes in `~/.delimp_*`. |

## Spectronaut (Run Comparator)

| Problem | Solution |
|---------|----------|
| Trailing dots in sample names | Spectronaut appends `.` to labels ending in digits ("AD12." → "AD12"). `match_samples()` strips with `gsub("\\.$", "", x)`. |
| `PG.UniProtIds` fallback | Some exports lack `PG.ProteinGroups`. Protein-column regex includes `UniProtIds`; Q-value regex includes `Q.Value` variant. |
| Quant3 inflates significance | "Use All MS-Level Quantities" doubles observation count (21v20 → 42v40 in t-test). Detected via `parse_spectronaut_search_settings()`; red "severe" row in settings diff. |
| `Group` not `ProteinGroup` | Candidates.tsv uses `Group` for accessions. Regex must include `^Group$`. |
| `Comparison (group1/group2)` format | Parenthetical suffix; `^Comparison$` fails — remove the `$` anchor. |
| 0-ratio proteins have NaN | Proteins with 0 `# of Ratios` → NaN logFC/Pvalue/Qvalue. `classify_de()` uses `is.finite()`; `assign_hypothesis()` coerces to safe defaults (0 logFC, 1 adjP). |
| `AnalyisOverview.txt` typo | Spectronaut misspells the filename. Regex handles both: `analy.?is.?overview`. |
| Spectronaut 20+ RunOverview format | Written as 2-column key-value pairs, not wide table. `parse_spectronaut_run_overview()` handles both via `ncol()` check. |

## Proteogenomics (v3.11.0)

| Problem | Solution |
|---------|----------|
| `jsonlite::fromJSON` on `NA_character_` round-trips as `list()` | `.init_status_json()` serializes `started_at = NA_character_` → `{}` → re-reads as `list()` (length 0). `nzchar(as.character(list()))` returns `logical(0)`, crashing `if`. Use the `.empty_or_str()` helper in `server_proteog_builder.R` for any value read from a parsed status.json. |
| Stage with no `job_id` was falsely marked complete | Stage 11 (`assemble`) is non-SLURM — old builds had `status:"unknown"` + `job_id:{}`. `.sacct_state("")` → garbage → `any_running` FALSE → wrongly "complete". Stages with no SLURM job_id are now treated as still-pending. |
| reactiveValues access outside a reactive consumer crashes session | `values$proteog_build_jobs` at module-entry (restore-from-disk) throws "Can't access reactive value outside of reactive consumer", killing all module observers. Wrap eager `values$*` reads at init in `isolate({...})`. |
| Persistence observer wrote empty list on startup, clobbering disk | Naive `observe({ saveRDS(values$proteog_build_jobs, ...) })` fires once on startup with the default empty list before restore. Gate: skip save if `length(jobs)==0 && file.exists(path) && length(readRDS(path))>0`. |
| Ensembl URL patterns are species-specific | `dna.primary_assembly.fa.gz` only exists for human + mouse. Most species use `dna.toplevel.fa.gz`. Default to `toplevel` for new reference genome additions. |
| Ensembl Plants moved to ebi.ac.uk mirror | `ftp.ensemblgenomes.org` returns nothing (DNS/refused). Use `ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-<N>/...`; path structure mirrors the old layout. |
| NCBI RefSeq accessions lack embedded gene symbols | `XP_*`/`NP_*`/`WP_*` headers have no gene names → empty `Genes`. `ncbi_download_proteome()` writes a side-car `<basename>_gene_map.tsv`; the proteog NCBI flow uploads BOTH files to Hive. |
| Proteog FASTA library catalog is per-user, FASTAs are shared | `~/.delimp_fasta_library/catalog.rds` is local (per "never mount drives for app state"); FASTAs live on shared Hive. Another member's catalog won't auto-populate. "Discover from Hive" scans `PROTEOG_RNASEQ_ROOT/*/status.json` and registers completed builds (idempotent). |

## Sage / Casanovo / Cascadia — Deployment (DDA pipeline)

| Problem | Solution |
|---------|----------|
| **Sage v0.14.7 `quant.lfq` schema** | We ship Sage v0.14.7 at `/quobyte/proteomics-grp/de-limp/cascadia/sage-v0.14.7.../sage`. In v0.14.7 `quant.lfq` is a **boolean** and settings live under `quant.lfq_settings`. v0.15+ removed the boolean and renamed the object `lfq`. Emitting the v0.15 shape → `Error: invalid type: map, expected a boolean`. `generate_sage_config()` in `helpers_dda.R` is correct for v0.14.7. Source: github.com/lazear/sage/blob/v0.14.7/DOCS.md. |
| **Sage can't read Thermo `.raw` directly** | Field is called `mzml_paths` but Sage v0.14.7 only parses mzML. Feeding `.raw` → `unhandled XML error: Unexpected EOF` + `0 spectra/s`, but Sage **still exits 0** (silent empty result). Always convert `.raw` → mzML via msconvert first. |
| **msconvert apptainer needs `--bind /quobyte:/quobyte`** | Apptainer doesn't auto-bind `/quobyte`. Without it, msconvert errors "no files found matching". Pattern: `apptainer exec --bind /quobyte:/quobyte $MSCONVERT_SIF wine msconvert ...`. |
| Casanovo v4 conda env typo | v4 env at `/quobyte/proteomics-grp/conda_envs/cassonovo_env` (`casso`, not `casa`) — Casanovo 4.3.0 + Python 3.10. Pair with `casanovo_v4_2_0.ckpt`. |
| **Casanovo v5 needs a SEPARATE env from v4** | `/quobyte/proteomics-grp/conda_envs/casanovo5/` = Casanovo 5.0.0 + Python 3.13 + depthcharge-ms 0.4.8 (`depthcharge.tokenizers`). v4 env CANNOT load a v5 ckpt (`ModuleNotFoundError: depthcharge.tokenizers`). CLI changed too: v4 `--model X --output FULL.mztab peak.mgf`; v5 `--model X --output_dir DIR --output_root NAME --force_overwrite peak.mgf`. `server_dda.R` routing reads `input$dda_casanovo_model`; `generate_casanovo_sbatch(casanovo_version=...)` branches the CLI. |
| `set -euo pipefail` + glob with no matches = silent exit 2 | `EXISTING_MGF=$(ls "$DIR"/*.mgf 2>/dev/null \| wc -l)` looks safe but `ls` exit 2 propagates through the pipe under `pipefail`, and `set -e` kills the script silently. Use `shopt -s nullglob`, OR `\|\| true`, OR `find`. |
| **De novo→homology FDR: calibrate by BLAST bitscore, NOT Casanovo score** | The decoy-spectra FDR must rank DIAMOND hits by **bitscore** (or E-value for standard searches), not by Casanovo peptide score. The Casanovo score barely separates real from decoy spectra — on the ocelot data only ~171 of 312,820 peptides clear 1% FDR when gated on Casanovo score, vs ~5,630 when gated on bitscore. Use `build_denovo_bitscore_fdr()` (helpers_dda.R), not the score-bin helpers. The "94–96% hit under EXTRA" is the RAW any-bitscore≥1 rate, NOT FDR-controlled. Full validation: `docs/denovo_decoy_method.html`. |
| **Two decoy types calibrate differently — decoy-spectra is the conservative/correct one** | Decoy-**spectra** (NovoBoard scrambled spectra) controls de novo *sequencing* error + chance matches; decoy-**database** (reversed proteins) controls only chance matches and over-accepts short peptides (at 1% FDR it accepts 6,184 peptides ≤20aa vs the spectra null's 101). For de novo→homology ID use the decoy-spectra null. Recommended: decoy-spectra + EXTRA sensitivity + bitscore @1% FDR + scramble FRAC≥0.8. |

## Cascadia / Casanovo — Training

| Problem | Solution |
|---------|----------|
| `preprocessing_fn` override replaces defaults | Cascadia uses `[scale_intensity("root"), scale_to_unit_norm]` only — NO peak filtering. `max_num_peaks=200` is a depthcharge default, not what the pretrained model expects. |
| Hidden LR scheduler in `configure_optimizers()` | CosineWarmupScheduler auto-activates. Override for fine-tuning with a flat LR. |
| OOM with full spectra | Median 9,558 peaks, max 113k. Use `batch_size=1` + `grad_accum=16`, or mobility-filtered data (~88 peaks). |
| PyTorch Lightning precision string | Old PL (1.x): `precision=16` (int), not `"16-mixed"` (PL 2.x). |
| spectrum_utils filter doesn't sync extra arrays | `filter_intensity`/`set_mz_range` only filter mz+intensity. Cascadia's rt/level/im/fragment arrays need manual sync. Patched in `primitives.py`. |
| Casanovo env missing deps | `pip check` and install: PyJWT, urllib3, Deprecated, rich. |
| Casanovo CLI version differences | Installed version uses `-o` (output prefix), not `-d` (directory). Check source, not docs. |
