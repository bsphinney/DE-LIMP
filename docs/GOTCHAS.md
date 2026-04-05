# DE-LIMP Gotchas Reference

Quick-reference table for known issues and their solutions. Referenced from CLAUDE.md.

## R Shiny / bslib

| Problem | Solution |
|---------|----------|
| Navbar text invisible on dark bg | Flatly theme CSS override: `.navbar .nav-link { color: rgba(255,255,255,0.75) !important; }` |
| Hidden tabs show letter fragments | `.navbar .nav-item[style*='display: none'] { width: 0 !important; overflow: hidden !important; }` |
| `page_navbar(bg=...)` deprecation | Use `navbar_options = navbar_options(bg = ...)` (bslib 0.9.0+) |
| `source()` doesn't start app | Use `shiny::runApp()` instead |
| Selections disappear after clicking | Reactive loop — table must not depend on selection-derived reactives |
| bslib `card()` doesn't render | Use plain `div()` for top-level nav_panel content |
| `uiOutput` vanishes in `navset_card_tab` | Use static HTML + `shinyjs::html("div_id", content)`. `plotlyOutput` with `req()` is safe. |
| `return()` inside `withProgress` | Exits `withProgress` not enclosing function. Use flat `tryCatch`. |
| `<<-` inside `withProgress` fails | `withProgress` uses `eval(substitute(expr), env)`. Use `new.env()` + `<-` instead. |
| Shiny hidden input not registered by JS | Use `div(style="display:none;", radioButtons(...))` for `conditionalPanel` |
| TIC extraction auto-triggered | `observeEvent(list(btn, trigger))` fires when button renders (NULL→0). Use separate `reactiveVal`. |

## DIA-NN

| Problem | Solution |
|---------|----------|
| DIA-NN `Genes` column has accessions | Not gene symbols. Real genes from `bitr()` UNIPROT → SYMBOL. |
| `readDIANN` data.table column error | Must pass `format="parquet"` for .parquet files |
| DIA-NN empirical lib is `.parquet` not `.speclib` | DIA-NN 2.0+ saves empirical libraries in parquet format |
| `--quant-ori-names` required on ALL steps | Preserves original filenames in `.quant` files across container bind mounts |
| `--fasta-search`/`--predictor` Step 1 only | Including in Steps 2-5 causes full FASTA re-digest |
| Auto mass acc + `--use-quant` | Produces different results. Force `mass_acc_mode = "manual"` |
| `max_pr_mz` default was 1200 not 1800 | DIA-NN default is 1800. Fixed in UI and fallbacks. |
| Parallel search OOM on timsTOF | Default `mem_per_file` was 32 GB, now 64 GB |
| Two DIA-NN containers with similar names | `/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif` (underscore, `dia-nn/` dir) has .NET, reads .raw. `/quobyte/proteomics-grp/apptainers/diann2.3.0.sif` (no underscore, `apptainers/` dir) does NOT read .raw. Always use the `dia-nn/` directory version. |

## Data & Columns

| Problem | Solution |
|---------|----------|
| `nrow(raw_data$E)` counts precursors not proteins | Use `length(unique(raw_data$genes$Protein.Group))` for protein groups |
| y_protein `colSums` error | It's limma `EList`. Extract `$E` for expression matrix. |
| `arrow::select` masks `dplyr::select` | Use `dplyr::select()` explicitly |
| R regex `\\s` invalid | Use `[:space:]` in base R regex (POSIX ERE) |
| `unlist()` on nested lists causes row mismatch | Use `vapply(x, function(v) paste(v, collapse="; "), character(1))` |
| Character matrix subsetting fails on Linux | Use numeric indices via `match()`, not rowname subsetting |
| Volcano P.Value vs adj.P.Val mismatch | Y-axis uses raw P.Value; dashed line at `max(P.Value)` among adj.P.Val < 0.05 |
| MOFA2 views need same sample names | Subset to matched pairs, assign common labels |
| Contaminant proteins have `Cont_` prefix | Detect via `grepl("^Cont_", Protein.Group)` |
| No-replicates mode skips DE | `values$fit` remains NULL. Expression Grid, PCA still work. |
| NCBI RefSeq accessions need gene mapping | Use batch E-utilities for gene symbols. Cache as TSV alongside FASTA. |

## SSH / HPC

| Problem | Solution |
|---------|----------|
| Symlinks in container bind mounts don't resolve | If you bind `selected_raw:/work/data` and it contains symlinks pointing to `../8min/file.raw`, the container can't follow them — the symlink target path isn't mounted. Fix: bind the parent dir (`--bind /base:/work`) so all subdirs are visible, or use `cp` instead of `ln -sf`. |
| SSH output encoding crash | `iconv(..., sub="")` in `ssh_exec`/`scp_download`/`scp_upload` |
| SSH rapid connections rejected (255) | Batch operations; use ControlMaster multiplexing |
| macOS SSH ControlPath too long | Use `/tmp/.delimp_<user>_<host>` (104 byte limit) |
| `parse_sbatch_output` returns dirty ID | Always `trimws()` parsed job IDs |
| SSH auto-connect blocks event loop | Run via `later::later()` or fast-fail timeout |
| SLURM limits on QOS not associations | Use `sacctmgr show qos` not `sacctmgr show assoc` for limits |
| Per-user CPU limit is binding | `MaxTRESPU` constrains users, not `GrpTRES` |
| `sacct` `.extern` step falsely reports COMPLETED | Filter out substep lines (those containing `.`) |
| Array progress sacct inflated counts | Filter to only `JOBID_N` format: `grepl("_", jid) && !grepl("\\.", jid)` |
| Partial retry dependency chain | `scontrol update` step 3's dependency after retrying step 2 |
| SLURM proxy inside Apptainer | Proxy process outside container relays commands via temp files |
| Paths with spaces in SLURM scripts | Quote all paths in sbatch |
| Docker container name rejected | Sanitize with `gsub("[^a-zA-Z0-9_.-]", "_", name)` |
| Docker SSH key permissions | Copy keys to container-internal volume with `chmod 600` |
| Log import ignores `fr_mz`/`pr_charge` | Parse via `value_map` in `parse_diann_log` |
| Load from HPC needs build-time guard | Wrap in `if (!is_hf_space)` in `build_ui()` |
| **NEVER use mounted drives for app state** | Use local `~/.delimp_*` paths. Cross-user sharing via SSH/SCP. |
| **Derived data stays with source data** | TIC cache, session.rds in data/output dir, not home dir |
| Python stdout buffered in SLURM | Fix: `sys.stdout.reconfigure(line_buffering=True)` |

## Spectronaut

| Problem | Solution |
|---------|----------|
| Trailing dots in sample names | Strip with `gsub("\\.$", "", x)` |
| `PG.UniProtIds` fallback | Protein column regex includes `UniProtIds` |
| Quant3 inflates significance | Doubles observation count. Red "severe" row in settings diff. |
| `Group` not `ProteinGroup` | Candidates.tsv uses `Group`. Regex must include `^Group$`. |
| `Comparison (group1/group2)` format | Don't anchor regex with `$` |
| 0-ratio proteins have NaN | `classify_de()` uses `is.finite()`, `assign_hypothesis()` coerces to safe defaults |
| `AnalyisOverview.txt` typo | Regex handles both spellings: `analy.?is.?overview` |
| Spectronaut 20+ RunOverview format | Key-value pairs, auto-detected via `ncol()` check |

## Cascadia / Casanovo Training

| Problem | Solution |
|---------|----------|
| `preprocessing_fn` override replaces defaults | Cascadia uses `[scale_intensity("root"), scale_to_unit_norm]` only — NO peak filtering |
| Hidden LR scheduler in configure_optimizers() | CosineWarmupScheduler auto-activates. Override for fine-tuning with flat LR. |
| OOM with full spectra | Median 9,558 peaks, max 113k. batch_size=1 + grad_accum=16, or use mobility-filtered data (~88 peaks) |
| PyTorch Lightning precision string | Old PL (1.x): `precision=16` (int), not `"16-mixed"` (PL 2.x) |
| spectrum_utils filter doesn't sync extra arrays | `filter_intensity`/`set_mz_range` only filter mz+intensity. Cascadia's rt/level/im/fragment arrays need manual sync. Patched in primitives.py. |
| Casanovo env missing deps | `pip check` and install: PyJWT, urllib3, Deprecated, rich |
| Casanovo CLI version differences | Installed version uses `-o` (output prefix), not `-d` (directory). Check source, not docs. |
