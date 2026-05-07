# Input Picker Design — Per-File Raw Selection + Output Directory Picker

Status: **draft for review, 2026-05-07**
Companion to: `HISTORY_DB_DESIGN.md`
Trigger: Brett's UX request after the v3.10.4–v3.10.29 install-stack stabilization

---

## 1. Two related improvements

| # | Feature | Where | Status |
|---|---|---|---|
| **A** | Per-file selection in the raw-data scan flow (modal with checkboxes when N>1, like the v3.10.4 FASTA picker) | New Search panel — both local browser and SSH scan paths | Design ready, low risk |
| **B** | Output directory picker button + "Create new subdirectory" option | SSH HPC mode (currently has a hidden input + auto-derived path); also enrich Local + Docker pickers with the new-subdir affordance | Design needs decision; medium risk |

---

## 2. Design — Per-file raw selection (A)

### 2.1 Current behavior

Two scan paths, both auto-select **every** raw file in the chosen directory:

| Path | Trigger | Handler | State mutation |
|---|---|---|---|
| **Local** (Docker / local-sbatch HPC) | `shinyFiles::shinyDirButton("raw_data_dir")` → user picks dir | `R/server_search.R:2277` | `values$diann_raw_files <- scan_raw_files(dir)` |
| **SSH** (HPC + SSH mode) | `actionButton("ssh_scan_raw_btn")` → reads `input$ssh_raw_data_dir` | `R/server_search.R:1922` | `values$diann_raw_files <- ssh_scan_raw_files(cfg, dir)` |

Both call `scan_raw_files()` / `ssh_scan_raw_files()` (in `helpers_search.R`) which return a data.frame with `filename`, `size_mb`, `mode`, `mtime` columns. The full data.frame is stored in `values$diann_raw_files`.

### 2.2 Consumers of `values$diann_raw_files`

12+ read sites across `server_search.R`. None mutate; all just iterate or count. Confirmed no consumer assumes a specific row count or order — they all just operate on whatever's there.

This is the **safe-to-change** signal: if we replace the full data.frame with a filtered subset, every consumer keeps working.

### 2.3 Proposed flow

When scan finds **N > 1** files:

1. Pop a `modalDialog` listing all files
2. `checkboxGroupInput`, **all pre-checked** (preserves current "use everything" default)
3. Helper buttons in the modal: **Select All** / **Select None** / **Invert**
4. **Confirm** → filter the data.frame to picked rows, write to `values$diann_raw_files`, `removeModal()`
5. **Cancel** → leave state unchanged
6. **Confirm with 0 selected** → notification "pick at least one", modal stays open

When scan finds **N == 1**: skip the modal entirely (current "use it" behavior, like the FASTA picker).

When scan finds **N == 0**: existing "no raw files found" notification, no change.

### 2.4 UX details

- **Many files (50+)**: wrap the `checkboxGroupInput` in a `div(style = "max-height: 60vh; overflow-y: auto;")`. Keeps modal usable.
- **File metadata in the labels**: render each row as `filename (size MB, mode)` so users have context. e.g. `Sample_03.d (4,800 MB, dia-PASEF)`.
- **Sort order**: alphabetical by filename — predictable.
- **Search/filter box** (nice-to-have, not blocking): `textInput` at top of modal that filters the visible list. Skipping for v1.

### 2.5 Files touched (estimated)

| File | Lines | Risk |
|---|---|---|
| `R/server_search.R:1922-1945` | ~30 changed (SSH scan path) | Low — same shape as FASTA picker |
| `R/server_search.R:2277-2284` | ~30 changed (local scan path) | Low |
| New observers for confirm/cancel buttons | ~50 added | Low |

Total: ~110 lines added, ~40 modified. ~1 hour of focused work + test pass.

### 2.6 Risks

- **Metadata extraction during scan**: if `scan_raw_files()` does heavy work (TIC extraction, instrument metadata) for ALL files, that's wasted on excluded ones. Inspect `scan_raw_files()` before implementation; if metadata is light (just `file.info` stat calls), no change needed. If heavy, defer metadata to post-pick.
- **Re-scan after a selection**: clicking scan again must fully replace state. The simple "write the data.frame" pattern handles this naturally. Verify in test.
- **TIC extraction trigger**: TIC trace generation may be wired off `values$diann_raw_files` change. If so, TIC fires only for picked files — that's actually the desired behavior, but worth confirming.

---

## 3. Design — Output directory picker (B)

### 3.1 Current state

Three modes, three different output-dir handlings:

| Mode | UI | Behavior |
|---|---|---|
| **Local** (`input.search_backend == 'local'` && no container) | `shinyFiles::shinyDirButton("local_output_dir_browse")` + `verbatimTextOutput("local_output_path")` | User picks a directory; output goes to `<chosen>/<analysis_name>_<timestamp>/` |
| **Local container** (`delimp_data_dir` env set) | `textInput("local_output_dir", value = file.path(delimp_data_dir, "output"))` (read-only-ish) | Fixed to `${DELIMP_DATA_DIR}/output/<analysis_name>_<timestamp>/` |
| **Docker** | `shinyFiles::shinyDirButton("docker_output_dir")` + `verbatimTextOutput("docker_output_path")` | Same as Local: pick parent, app appends `<analysis_name>_<timestamp>/` |
| **HPC + SSH** | `textInput("ssh_output_base_dir")` **wrapped in `div(style = "display: none;")`** | **Hidden by design.** Output auto-derived from `<ssh_raw_data_dir>/<analysis_name>_<timestamp>/`. User cannot override. |

The Local and Docker modes already let the user pick. **Only SSH is hidden + auto-derived.**

### 3.2 Open question for Brett

Brett's request: *"output directory select button by the output directory box that allows you to select the output directory like the scan raw data folder. You should add the option to make a new sub directory too for the above feature."*

Two interpretations:

- **Option 1 — unhide the SSH input + add a Browse button** so SSH HPC users can override the auto-derived output path, browsing remote HPC dirs via the existing SSH file browser
- **Option 2 — keep SSH hidden; the request is about Local/Docker, where the picker exists but doesn't have a "create new subdir" affordance**

The auto-derivation in SSH mode is load-bearing for the path-resolution logic v3.10.20–25 just stabilized (Load-from-HPC, activity-log linking, Recover button). Letting users override paths means the resolver has to handle arbitrary paths, not just the conventional `<raw_dir>/<analysis_name>_<timestamp>` shape.

### 3.3 If Option 1 (unhide SSH + Browse)

**UI changes** (`ui.R:1037`):
- Remove the `div(style = "display: none;")` wrapper around `ssh_output_base_dir`
- Place the input + a "Browse" button + a "Create New Subdirectory" button in a horizontal flex row, matching the SSH raw-data row pattern
- A `verbatimTextOutput("ssh_output_path_display")` showing the resolved full path (`<base>/<analysis_name>_<timestamp>/`)

**Browse handler**:
- Reuses the existing SSH file browser modal (the one used for raw / FASTA browsing)
- Restricted to directories (no file pick)
- Confirm sets `ssh_output_base_dir` value via `updateTextInput`

**Create-new-subdir**:
- Modal with `textInput` for the subdir name + parent path readout
- Validate: alphanumeric + `_-` only, no slashes, max 100 chars
- On confirm: SSH `mkdir -p <parent>/<new>` with error handling for permission denied
- On success: `updateTextInput("ssh_output_base_dir", value = <parent>/<new>)`

**Risks**:
- Path-resolution code in Load-from-HPC, activity log, Recover — already handles arbitrary `output_dir` values now (per v3.10.20+). Should be fine.
- Users may pick a path the SSH user can't write to. `mkdir -p` failure surfaces in a notification.
- The auto-derivation fallback should still work when the input is empty (preserves backward compat for users who haven't typed anything).

### 3.4 If Option 2 (Local/Docker only — add subdir affordance)

**UI changes** to existing pickers:
- Below `verbatimTextOutput("local_output_path")`, add a small inline `textInput("local_output_subdir")` + a "Create" button
- Same for Docker mode

**Behavior**:
- User picks parent dir via existing `shinyDirButton`
- Optionally types a subdir name
- "Create" button does `dir.create(file.path(parent, subdir), recursive = TRUE)` and updates the displayed path

**Risks**:
- Lower risk than Option 1 — picker already exists, just adding a small input + button
- File-system errors on create (permissions, disk full) need to surface

### 3.5 My recommendation

I'd lean **Option 1 + the subdir affordance from Option 2 in all three pickers** — uniformity across modes. Estimated scope: ~150 lines across `ui.R` and `server_search.R`. ~1.5 hours of work with careful testing.

But this depends on Brett's intent. If the request was really about Local/Docker only, Option 2 alone is half the scope and lower risk.

---

## 4. Combined rollout plan

**Phase 1 — Raw file picker** (low risk, value/effort favorable)

- Single PR touching `server_search.R` only
- Reuses FASTA-picker pattern verbatim
- Ship as **v3.11.0a** (or v3.10.30 if we're still in patch territory)
- Stress-test before moving to Phase 2

**Phase 2 — Output dir picker** (depends on Brett's answer to §3.2)

- Separate PR after Phase 1 lands and stress-tests cleanly
- Ship as **v3.11.0b** (or v3.10.31)

Two-PR split keeps each diff reviewable + isolates regressions.

---

## 5. Decision log

- **2026-05-07**: design drafted by Claude after Brett's UX request post-v3.10.29 install-stack stabilization. Brett asked to audit before coding given today's 26-hotfix run.
- Open: Brett to confirm Option 1 vs 2 vs 1+2 for the output-dir picker scope (§3.2).
- Open: Brett to confirm "shinyFiles SSH file browser" already exists for the proposed Browse button — if not, that's an additional design item.

---

## 6. Pre-implementation checklist

Before any code lands:

- [ ] Brett confirms Option 1 / 2 / 1+2 for the output-dir picker
- [ ] Inspect `scan_raw_files()` to confirm metadata-extraction cost during scan
- [ ] Confirm `values$diann_raw_files` is the single source of truth (no parallel state to keep in sync)
- [ ] Run a pre-existing search end-to-end on current code so we have a known-good baseline
- [ ] Implement Phase 1 only, ship, stress-test
- [ ] Then Phase 2
