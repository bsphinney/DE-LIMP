# DE-LIMP Implementation Plan
## timsrust Test → Three-Mode Switcher → DDA Workflow → XL-MS Module

> **Philosophy:** Test first, build on confirmed ground. Each phase has a hard
> checkpoint before the next phase starts. Nothing gets committed to the app until
> it works standalone. Claude Code sessions are scoped to single files or tightly
> related file pairs — never more than ~400 lines per session to stay within
> context limits.
>
> **Reading order for Claude Code:** This plan → relevant spec → implement.
> Never implement from memory of a previous session — always re-read the spec.

---

## Pre-Implementation Checklist (You, ~30 min, before any coding)

Run these manually on Hive before starting. Record the results — they determine
which branches of the plan apply.

```bash
# 1. Rust toolchain
which cargo && cargo --version
# If missing: curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# 2. Java
module load java/17 && java -version
# Need 17+

# 3. MeroX
ls /share/proteomics/tools/xlms/MeroX_2025.jar
java -jar /share/proteomics/tools/xlms/MeroX_2025.jar --help 2>&1 | head -5

# 4. xiSearch + xiFDR
ls /share/proteomics/tools/xlms/xiSearch.jar
ls /share/proteomics/tools/xlms/xiFDR.jar

# 5. Sage
which sage || ls /share/proteomics/tools/sage
sage --version 2>/dev/null || echo "need to install"

# 6. A real .d file for testing
ls /path/to/your/xlms_data/*.d | head -3
ls /path/to/your/dda_data/*.d | head -3

# 7. DE-LIMP current state
cd /path/to/delimp
git status
git log --oneline -5
Rscript -e "library(limpa); library(visNetwork); library(jsonlite); cat('deps OK\n')"
```

Record:
- [ ] Rust available? Y/N — if N, install before Phase 0
- [ ] Java 17 module name on Hive
- [ ] MeroX JAR path (confirm exact path)
- [ ] xiSearch + xiFDR JAR paths
- [ ] Sage binary path (or "need to install")
- [ ] Path to a real XL-MS `.d` file
- [ ] Path to a real DDA `.d` file
- [ ] DE-LIMP git branch to work on

---

## Phase 0 — timsrust MGF Test (Tomorrow Morning, ~1 hr)

**Goal:** Confirm timsrust `MGFWriter` produces MGF files that MeroX accepts.
This is the single biggest unknown. Everything else in the XL-MS pipeline is
proven technology. Do not skip this phase.

**Do manually on Hive — not with Claude Code.**

### Step 0.1 — Write the test binary (~10 min)

```bash
mkdir ~/timsrust_test && cd ~/timsrust_test
cargo init --name timsrust_mgf_test
```

`Cargo.toml`:
```toml
[package]
name = "timsrust_mgf_test"
version = "0.1.0"
edition = "2021"

[dependencies]
timsrust = "0.4.2"
```

`src/main.rs`:
```rust
use timsrust::readers::SpectrumReader;
use timsrust::writers::MGFWriter;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: timsrust_mgf_test <input.d> <output.mgf>");
        std::process::exit(1);
    }

    let d_path  = &args[1];
    let mgf_out = &args[2];

    println!("Reading: {}", d_path);
    let reader = SpectrumReader::new(d_path.as_str())
        .expect("Failed to open .d file");

    println!("Found {} spectra", reader.len());

    let spectra: Vec<_> = reader.get_all()
        .into_iter()
        .filter_map(|r| r.ok())
        .collect();

    println!("Loaded {} spectra successfully", spectra.len());

    MGFWriter::write_spectra(mgf_out.as_str(), &spectra);
    println!("Written to {}", mgf_out);

    // Print first 30 lines so we can inspect the format
    let preview = std::fs::read_to_string(mgf_out).unwrap();
    println!("\n--- First 30 lines of MGF ---");
    preview.lines().take(30).for_each(|l| println!("{}", l));
}
```

### Step 0.2 — Build and run (~10 min)

```bash
cd ~/timsrust_test

# Build (takes ~2-3 min first time, downloading crates)
cargo build --release 2>&1 | tail -5

# Run on a real .d file
./target/release/timsrust_mgf_test \
  /path/to/real_xlms_sample.d \
  /tmp/test_output.mgf
```

### Step 0.3 — Inspect the MGF (~5 min)

Check these fields are present and populated:

```bash
# Count spectra
grep -c "BEGIN IONS" /tmp/test_output.mgf

# Check first spectrum block
head -40 /tmp/test_output.mgf

# Verify required fields
grep "PEPMASS" /tmp/test_output.mgf | head -3
grep "CHARGE"  /tmp/test_output.mgf | head -3
grep "RTINSECONDS\|TITLE" /tmp/test_output.mgf | head -3

# Check m/z and intensity pairs look real (non-zero)
awk '/BEGIN IONS/{f=1} f && /^[0-9]/{print; exit}' /tmp/test_output.mgf
```

### Step 0.4 — Decision gate

**If MGF looks good** (PEPMASS present, intensities non-zero, spectrum count reasonable):
- Copy the binary to `/share/proteomics/tools/timsrust_mgf`
- Record the exact binary path
- Proceed to Phase 1

**If MGF is malformed or missing fields:**
- Switch plan to timsconvert → mzML (update XLMS spec §1 accordingly)
- timsconvert command: `timsconvert --input sample.d --output mzml/ --ms2_only True`
- Both MeroX and xiSearch accept mzML — no other changes needed to the spec
- Proceed to Phase 1 with mzML conversion instead

**If timsrust fails to compile:**
- Check `cargo build` error — likely a missing C library
- Try: `module load gcc/11` then retry
- Fallback to timsconvert as above

---

## Phase 1 — Shared Infrastructure (~1 hr, Claude Code Session 1)

**Goal:** Three-mode switcher, config.yml, feature flags, `app.R` reactiveValues.
This is the foundation everything else builds on. Do it first so DDA and XL-MS
sessions don't step on each other.

**Spec to read:** `THREE_MODE_SWITCHER_ADDENDUM.md` (entire document)

**Claude Code prompt:**

> "Read THREE_MODE_SWITCHER_ADDENDUM.md in full. Then implement the following
> in DE-LIMP:
>
> 1. In `config.yml`: add `features.enable_dda` and `features.enable_xlms` flags,
>    and the full `tools:` block from §5 of the addendum with placeholder paths
>    (I will fill in real paths).
>
> 2. In `app.R` reactiveValues: add the `xlms` namespace block from
>    XLMS_INTEGRATION_SPEC.md §4. Add `acquisition_mode = 'dia'` to the existing
>    values block. Source `R/server_xlms.R` and `R/server_dda.R` (stubs OK for now).
>
> 3. In `ui.R`: replace the existing two-choice acquisition mode switcher (from the
>    DDA spec §2) with the three-choice version from addendum §1. Wire the
>    `mode_choices` gating logic from addendum §5 so XL-MS and DDA options only
>    appear when `is_hive` is TRUE.
>
> 4. In `server_dda.R` (or wherever the mode observer lives): extend the
>    `observeEvent(input$acquisition_mode)` handler to handle 'xlms' per addendum §3.
>
> 5. Create empty stub files: `R/server_xlms.R` and `R/ui_xlms.R` with a comment
>    header only.
>
> Do NOT implement any DDA or XL-MS pipeline logic yet — stubs only.
> Run the app after each file change to confirm it starts without errors."

**Checkpoint:** App starts, mode switcher shows three options on Hive, one option
(DIA only) on HF Spaces. No functionality yet — just the switcher.

---

## Phase 2 — DDA: Sage Search Pipeline (~3 hrs, Claude Code Sessions 2-4)

Work through the DDA spec in three focused sessions. Each session ends with a
working checkpoint before moving on.

**Spec to read:** `DDA_SAGE_CASANOVO_SPEC.md` — read the entire document before
starting Session 2.

### Session 2 — Sage Config + SLURM Submission (~1 hr)

> "Read DDA_SAGE_CASANOVO_SPEC.md §5 (Phase 1: Sage Search) in full.
>
> Create `R/helpers_dda.R` containing:
> - `generate_sage_config()` from §5.1
> - `parse_sage_results()` from §5.3
> - `build_dda_elist()` from §7
>
> Create `R/server_dda.R` containing:
> - The `dda$` reactive namespace observer from §9
> - `submit_sage_job()` SLURM submission from §11
> - The job polling loop from §9
> - The results loader that calls helpers_dda.R functions on completion
>
> Do NOT implement the UI yet. Test `generate_sage_config()` standalone:
> ```r
> source('R/helpers_dda.R')
> cfg <- generate_sage_config(
>   fasta_path = '/tmp/test.fasta',
>   raw_paths  = c('/tmp/sample1.d', '/tmp/sample2.d'),
>   output_dir = '/tmp/sage_test',
>   preset     = 'standard'
> )
> cat(readLines(cfg), sep='\n')
> ```
> Confirm JSON is valid with `jsonlite::fromJSON(cfg)`."

**Checkpoint:** `generate_sage_config()` produces valid Sage JSON. `parse_sage_results()`
runs without error on a real `results.sage.tsv` file.

### Session 3 — DDA UI Panel + QC Integration (~1 hr)

> "Read DDA_SAGE_CASANOVO_SPEC.md §10 (UI Spec) and §16 (Normalization & Imputation).
>
> In `ui.R`, add the DDA conditionalPanel inside the Search tab per §10.2.
>
> In `R/helpers_dda.R`, add:
> - `normalize_dda_matrix()` from §16.2
> - `filter_dda_valid_values()` from §16.3
> - `impute_dda_matrix()` from §16.3
> - `run_dda_pipeline()` from §16.5
>
> In `server_dda.R`, add the normalization/imputation UI controls from §16.6
> and wire `run_dda_pipeline()` to trigger after Sage results load.
>
> In `server_qc.R`, add the DDA QC summary card from §10.3 (mode-conditional).
>
> Run the app and confirm the DDA panel appears when mode switcher = DDA."

**Checkpoint:** Full DDA panel visible. Submit a real Sage job on Hive. Confirm
`results.sage.tsv` and `lfq.tsv` are produced. Load into DE-LIMP. QC card shows
correct PSM/protein counts.

### Session 4 — DDA → DE Pipeline Wiring (~1 hr)

> "Read DDA_SAGE_CASANOVO_SPEC.md §7 (R-Side Integration) and §12 (Session Save/Load).
>
> In `server_de.R`, add mode-awareness so the limma DE pipeline runs on
> `values$dda$elist` when `values$acquisition_mode == 'dda'`, and on the
> existing `values$y_protein` EList when `acquisition_mode == 'dia'`.
> The DE pipeline code itself does NOT change — only which EList it reads from.
>
> In `server_session.R`, add DDA state to save/load per §12.
>
> Test end-to-end: DDA mode → Sage search → load results → run DE → confirm
> volcano plot renders with DDA protein names."

**Checkpoint:** Full DDA workflow runs end-to-end. Volcano plot, DE table, and
QC plots all show DDA data when in DDA mode. Switch back to DIA mode — confirm
zero regression (existing DIA session loads cleanly, volcano still works).

---

## Phase 3 — XL-MS: Core Search Pipeline (~3 hrs, Claude Code Sessions 5-7)

**Spec to read:** `XLMS_INTEGRATION_SPEC.md` — read the entire document before
starting Session 5.

### Session 5 — Settings Generators + SLURM Arrays (~1 hr)

> "Read XLMS_INTEGRATION_SPEC.md §7.1 through §7.3 in full.
>
> Create `inst/xlms/crosslinkers.json` from §5.
> Create `inst/xlms/template_merox.mxf` — a valid MeroX XML template with
> {CROSSLINKER_ID}, {ENZYME}, {MISSED_CLEAVAGES}, {PRECURSOR_TOLERANCE},
> {FRAGMENT_TOLERANCE}, {FDR_THRESHOLD} placeholders. Base this on MeroX's
> documented settings format (precursor_ppm, fragment_ppm, enzyme, crosslinker_id,
> fdr fields in XML).
> Create `inst/xlms/template_xi.config` — a plain text xiSearch config template
> with the same placeholders.
>
> In `R/server_xlms.R`, implement:
> - `crosslinker_library` load from JSON at top of file
> - `generate_merox_mxf()` from §7.2
> - `generate_xi_config()` from §7.2
> - `submit_conversion_array_job()` from §7.3
> - `submit_merox_array_job()` from §7.3
> - `submit_xi_array_job()` from §7.3
> - The main submit `observeEvent(input$xlms_submit)` handler from §7.3
>
> Test the settings generators standalone:
> ```r
> source('R/server_xlms.R')
> # Check MeroX settings file generates without error
> generate_merox_mxf(
>   list(crosslinker='DSBU', enzyme='trypsin', missed_cleavages=3,
>        ms1_ppm=10, ms2_ppm=10, fdr_threshold=0.05),
>   '/tmp/test_settings.mxf'
> )
> cat(readLines('/tmp/test_settings.mxf'), sep='\n')
> ```"

**Checkpoint:** Both settings generators produce valid output files. SLURM array
submission functions are written (not yet tested on a real job — that comes next).

### Session 6 — Job Monitor + Result Parsers + Merge (~1 hr)

> "Read XLMS_INTEGRATION_SPEC.md §7.4 through §7.6 in full.
>
> In `R/server_xlms.R`, implement:
> - `check_array_job()` from §7.4
> - The `xlms_poll_timer` reactive and observe block from §7.4
> - `output$xlms_job_monitor` renderUI from §7.4
> - `parse_merox_results()` from §7.5
> - `parse_xi_results()` from §7.5
> - `make_site_key()` from §7.5
> - `merge_crosslinks()` from §7.6
> - `build_protein_pairs()` from §7.6
> - `parse_and_merge_results()` auto-trigger from §7.8
>
> Test the parsers against real output files if available, or construct minimal
> test CSVs matching the expected column structure.
> Test `merge_crosslinks()` with two small tibbles — verify GOLD assignment when
> site_key appears in both, SILVER when in one only."

**Checkpoint:** `merge_crosslinks()` correctly assigns GOLD/SILVER. `make_site_key()`
order-normalizes correctly (swap A/B — same key). Job monitor renderUI produces
the five-stage badge display without error.

### Session 7 — XL-MS UI + Setup Panel (~1 hr)

> "Read XLMS_INTEGRATION_SPEC.md §6 (UI) in full.
>
> Create `R/ui_xlms.R` implementing:
> - `xlms_setup_ui()` from §6.1 (complete, including shinyDirButton, crosslinker
>   dropdown, engine checkboxes, HPC resource inputs, submit button)
> - `xlms_results_ui()` from §6.3 (navset_tab with QC, Network, Crosslinks Table)
> - `xlms_export_ui()` from §6.4
>
> In `server_xlms.R`, implement:
> - `output$xlms_crosslinker_info` reactive info card from §6.2
> - `output$xlms_summary_banner` from §9
>
> In `ui.R`, wire the XL-MS conditionalPanel to show `xlms_setup_ui()`,
> `xlms_results_ui()`, and `xlms_export_ui()` when `acquisition_mode == 'xlms'`.
>
> Run the app — confirm XL-MS panel appears when mode switcher = XL-MS.
> Click through all sub-panels. Confirm no JS errors in browser console."

**Checkpoint:** XL-MS UI fully visible and navigable. Submit button present.
No app crashes switching between all three modes repeatedly.

---

## Phase 4 — XL-MS: Visualization + Export (~2 hrs, Claude Code Sessions 8-9)

### Session 8 — visNetwork + QC Plots (~1 hr)

> "Read XLMS_INTEGRATION_SPEC.md §8 (visNetwork), §9 (summary banner), and the
> QC plots section of §6.3.
>
> In `server_xlms.R`, implement:
> - `output$xlms_network` visNetwork render from §8
> - `output$xlms_edge_drilldown` edge click handler from §8
> - `output$xlms_xsm_bar` — bar chart: XSMs per sample, colored by engine
>   (simple ggplot, use `xlms$qc_stats` data)
> - `output$xlms_engine_venn` — engine overlap visualization using
>   `xlms$qc_stats$n_gold`, `n_silver_merox`, `n_silver_xi`
>   (use ggvenn or a simple ggplot-based Venn with three circles)
> - `output$xlms_crosslink_table` DT table from GOLD or all crosslinks
>   filtered by `input$xlms_table_filter`
>
> Test with synthetic `xlms$consensus` data (create a small fake tibble with
> 10 crosslinks, 3 GOLD and 7 SILVER) to verify the network renders and
> edge clicks return drilldown data."

**Checkpoint:** Network renders with test data. Gold edges are solid, silver dashed.
Clicking an edge shows the drill-down table. QC bar chart and Venn render.

### Session 9 — Claude Export ZIP (~1 hr)

> "Read XLMS_INTEGRATION_SPEC.md §10 (Claude Export ZIP) in full, including
> all three PROMPT.md template functions (§10.3).
>
> In `server_xlms.R`, implement:
> - `save_network_png()` from §11 using ggraph + igraph
> - `save_engine_overlap_png()` — save the same engine Venn as a PNG
>   using `ggsave()` on the same plot as `output$xlms_engine_venn`
> - `build_xlms_methods_text()` from §10.4
> - `build_xlms_prompt()` from §10.3 (all three report types: brief, full, manuscript)
> - `output$xlms_claude_export` downloadHandler from §10.2
>
> In `server_ai.R`, add the XL-MS export handler (or confirm it lives in
> server_xlms.R — keep it there for isolation).
>
> Test: click Download Claude ZIP with test data. Unzip locally.
> Verify all 8 expected files are present. Open PROMPT.md — confirm it contains
> real numbers from `xlms$qc_stats`, not placeholder strings."

**Checkpoint:** ZIP downloads. PROMPT.md contains correct experiment name,
crosslinker, GOLD count, top protein pairs. Network PNG is a real figure.
Methods_draft.txt is grammatically correct.

---

## Phase 5 — End-to-End Test on Real Data (~2 hrs, Manual + Claude Code)

**Do this yourself on Hive with real data before calling it done.**

### Real XL-MS test

```bash
# 1. Start DE-LIMP on Hive in a tmux session
tmux new -s delimp_test
cd /path/to/delimp
Rscript app.R  # or however you launch it

# 2. In browser:
#    - Switch to XL-MS mode
#    - Select raw .d files (real DSBU or DSSO XL-MS data)
#    - Set crosslinker = DSBU, enzyme = Trypsin, 3 missed cleavages
#    - Select both MeroX and xiSearch engines
#    - Click Submit

# 3. Watch job monitor — confirm all five stages show RUNNING then COMPLETED
squeue -u $USER  # should show array jobs

# 4. When complete:
#    - Check GOLD crosslink count is non-zero
#    - Network renders with real protein names
#    - Download Claude ZIP, upload to Claude.ai, verify report generates
```

### Real DDA test

```bash
# 1. Switch to DDA mode in the same DE-LIMP session
# 2. Point at real ddaPASEF .d files
# 3. Select standard preset, submit Sage job
# 4. After completion:
#    - Check PSM count is reasonable (thousands)
#    - Load metadata, run DE
#    - Confirm volcano renders with real proteins
```

### Regression test — DIA mode

```bash
# 1. Switch back to DIA mode
# 2. Load an existing DIA-NN report.parquet or session RDS
# 3. Run through full DIA pipeline
# 4. Confirm NOTHING broke — all existing tabs work identically
```

---

## Phase 6 — Casanovo (Optional, Add After Phase 5 Passes)

Casanovo is P1 priority in the DDA spec — do not block Phase 5 on it.

**Session 10 — Casanovo (if GPU available and depthcharge supports .d)**

> "Read DDA_SAGE_CASANOVO_SPEC.md §6 (Phase 2: Casanovo De Novo Overlay) in full.
> Check Phase 0 recon results for depthcharge .d support status.
> Implement whichever injection strategy (1, 2, or 3) was confirmed working."

---

## Risk Register and Mitigations

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| timsrust MGF malformed for MeroX | Medium | High | Fall back to timsconvert → mzML (15 min change) |
| MeroX headless `--export-csv` flag doesn't exist in 2025 version | Medium | Medium | Parse `.zhrm` directly via MeroX GUI on one file to understand output format, then parse programmatically |
| xiSearch per-file CSV header varies between versions | Low | Medium | Flexible `grep()` column matching already in parser |
| SLURM array `afterok` dependency syntax varies on Hive | Low | Low | Test dependency syntax with a trivial two-job test first |
| Sage doesn't read timsTOF `.d` natively on current Hive version | Low | Medium | Add mzML conversion step with timsconvert (already written) |
| visNetwork edge click event doesn't fire in Shiny context | Low | Low | Fall back to row selection in the crosslinks table for drill-down |
| DDA and DIA reactive pipelines contaminate each other | Low | High | Enforced by `values$acquisition_mode` gate in `server_de.R` — verify in Phase 5 regression test |
| xiFDR combined CSV header mismatch (column from one file missing in another) | Medium | Medium | `fill = TRUE` in `data.table::fread()` handles ragged CSVs |

---

## File Change Summary

| File | Phase | Change Type |
|---|---|---|
| `config.yml` | 1 | Add feature flags + tool paths |
| `app.R` | 1 | Add `xlms` reactiveValues, source new modules |
| `ui.R` | 1 | Three-choice mode switcher; XL-MS conditionalPanel |
| `R/server_dda.R` | 1+2 | Mode observer extension; full DDA pipeline |
| `R/helpers_dda.R` | 2 | New file: Sage config, parse, normalize, impute, EList |
| `R/server_qc.R` | 2 | DDA QC summary card (mode-conditional) |
| `R/server_de.R` | 2 | Mode-aware EList selection |
| `R/server_session.R` | 2+3 | DDA + XL-MS state save/load |
| `R/server_xlms.R` | 3+4 | New file: full XL-MS pipeline |
| `R/ui_xlms.R` | 3 | New file: XL-MS UI panels |
| `inst/xlms/crosslinkers.json` | 3 | New file: crosslinker library |
| `inst/xlms/template_merox.mxf` | 3 | New file: MeroX settings template |
| `inst/xlms/template_xi.config` | 3 | New file: xiSearch config template |
| `R/server_ai.R` | 4 | XL-MS Claude export (or keep in server_xlms.R) |

**Total new files:** 6
**Modified files:** 7
**Untouched files:** Everything else (server_data.R, server_gsea.R, server_phospho.R, etc.)

---

## Session Discipline Rules for Claude Code

These maximize the chance of each session succeeding:

1. **Always start by reading the relevant spec section** — never assume Claude Code
   remembers from a previous session. Paste the spec section or reference the file.

2. **One module per session max** — don't combine `server_xlms.R` and `ui_xlms.R`
   in the same session. Keep file scope tight.

3. **Test after every function** — don't write 400 lines then test. Write a function,
   test it standalone in R, then write the next.

4. **Commit after each checkpoint** — `git commit -m "Phase N checkpoint: [description]"`
   before starting the next session. Easy rollback if something goes wrong.

5. **Never modify `server_de.R` and `server_data.R` in the same session** — these are
   the highest-risk files for DIA regression. Touch them one at a time.

6. **Keep the DIA regression test fast** — load a saved `.rds` session and check that
   the volcano still renders. Do this after every session that touches `server_de.R`,
   `ui.R`, or `app.R`.

---

## Estimated Timeline

| Phase | Sessions | Wall time | Notes |
|---|---|---|---|
| Pre-flight checklist | — | 30 min | You, manual |
| Phase 0: timsrust test | — | 1 hr | You, manual on Hive |
| Phase 1: Mode switcher | 1 | 1 hr | Claude Code |
| Phase 2: DDA pipeline | 3 | 3 hrs | Claude Code |
| Phase 3: XL-MS search | 3 | 3 hrs | Claude Code |
| Phase 4: XL-MS viz + export | 2 | 2 hrs | Claude Code |
| Phase 5: End-to-end test | — | 2 hrs | You + real data |
| Buffer / debugging | — | 2 hrs | Always needed |
| **Total** | **9 sessions** | **~14 hrs** | Spread over 2-3 days |

Phase 0 tomorrow morning unblocks everything. Phases 1-4 can run in parallel once
Phase 0 is decided (timsrust vs timsconvert fallback). Phase 5 requires real XL-MS
data on Hive — schedule that for when you know the data is accessible.
