# Input Picker + Filename UX Design

Status: **draft for review, 2026-05-07**
Companion to: `HISTORY_DB_DESIGN.md`
Trigger: Brett's UX request after the v3.10.4–v3.10.29 install-stack stabilization

---

## 1. Problems we're solving

| # | Problem | Phase |
|---|---|---|
| **1** | Users can't find their output after a search completes — the default location is buried (`${DELIMP_DATA_DIR}/output/...`) and not exposed in any obvious "open this folder" affordance | **A** |
| **2** | The raw-data scan auto-selects every file in the folder; users with mixed datasets in one directory have no way to pick a subset | **B** |
| **3** | Local-container output mode (when `DELIMP_DATA_DIR` env is set) is a textbox-only input with no Browse picker. Plus pickers can't create new subdirectories on the fly. | **C** |
| **4** | Filenames with spaces / parens / shell-special characters silently break searches (DIA-NN command quoting, sbatch script paths). No detection, no warning, no auto-fix. | **D** |

Scope: **Local + Docker modes only.** SSH HPC mode keeps its current auto-derived output path — that's load-bearing for the v3.10.20–25 path-resolution stability and we don't risk it.

---

## 2. Phase A — Find-output UX (highest payoff)

### 2.1 The actual user complaint

> *"Some people can't find the default output location."*

That's the bug to fix. Picker + subdir doesn't directly solve "I lost the output." Adding visibility does.

### 2.2 Three changes

**A1 — Default location moves somewhere discoverable.**
Today: `${DELIMP_DATA_DIR}/output/` if env set, else `~/.delimp_output/` (hidden dir). Both are non-obvious.
New default: `~/Documents/DE-LIMP/searches/<analysis_name>_<timestamp>/`. Visible in Finder/Explorer. Configurable via `DELIMP_OUTPUT_ROOT` env var or `delimp_site()$output_root` (folds into the v3.10.15 site-config refactor).

**A2 — "Open Output Folder" button.** Appears next to the resolved-path display in the New Search panel AND on every queue/history row. Clicking opens the dir in the OS file browser:
```r
open_in_file_browser <- function(path) {
  if (Sys.info()[["sysname"]] == "Darwin")        system2("open", path)
  else if (Sys.info()[["sysname"]] == "Windows")  shell.exec(path)
  else                                            system2("xdg-open", path)
}
```
For containerized deployments we open the **host-side path**, not the container path.

**A3 — Show host-path AND container-path** when running inside a container. Read the bind-mount mapping from the launcher env (set in `launch_delimp.sh` / `Launch_DE-LIMP_Docker.bat`) and translate. Display as:
```
Container: /data/output/Test_2026...
Host:      C:\Users\brett\DE-LIMP\output\Test_2026...   [ Open Folder ]
```

### 2.3 Files touched (A)

- `R/helpers_site.R` — add `output_root` to `delimp_site()`
- `R/server_search.R` — `open_in_file_browser()` helper, "Open Folder" button observer
- `R/ui.R` — add the button next to output path displays in New Search, Queue, History
- `R/server_session.R` — same button on Complete Analysis exports

### 2.4 Risks / edges

- `xdg-open` not always available on minimal Linux. Fall back to copy-path-to-clipboard with a notification.
- Host-path translation: if the launcher didn't set the bind-mount mapping in the container env, gracefully fall back to "Container: ... (host path unknown — check your launcher's bind mount)".
- `~/Documents` may not exist on minimal Linux installs. Use `xdg-user-dir DOCUMENTS` if available; fall back to `~/DE-LIMP-output/`.

---

## 3. Phase B — Per-file raw picker

### 3.1 Current state

| Path | Trigger | Handler | State |
|---|---|---|---|
| **Local** (Docker / non-SSH HPC) | `shinyDirButton("raw_data_dir")` | `R/server_search.R:2277` | `values$diann_raw_files <- scan_raw_files(dir)` (auto-all) |
| **SSH** (HPC + SSH) | `actionButton("ssh_scan_raw_btn")` | `R/server_search.R:1922` | `values$diann_raw_files <- ssh_scan_raw_files(cfg, dir)` (auto-all) |

12+ read sites consume `values$diann_raw_files`. None mutate. **Safe to filter** — same shape (data.frame), just fewer rows.

### 3.2 Proposed flow

When scan finds **N > 1** files:
1. Pop a `modalDialog` listing all files
2. `checkboxGroupInput`, **all pre-checked** (preserves current default)
3. Helper buttons: **Select All** / **Select None** / **Invert**
4. **Confirm** → filter to picked, write `values$diann_raw_files`, `removeModal()`
5. **Cancel** → leave state unchanged
6. **Confirm with 0 selected** → notification, modal stays open

When **N == 1**: skip modal (use it directly).
When **N == 0**: existing "no raw files found" notification.

### 3.3 UX details

- Many files (50+): wrap in `div(style = "max-height: 60vh; overflow-y: auto;")`
- Labels: `filename (size MB, mode)` — e.g. `Sample_03.d (4,800 MB, dia-PASEF)` for context
- Sort: alphabetical
- Search/filter box: nice-to-have, defer to v2

### 3.4 Files touched (B)

- `R/server_search.R:1922-1945` (SSH scan, ~30 lines changed)
- `R/server_search.R:2277-2284` (local scan, ~30 lines changed)
- New observers for confirm/cancel (~50 lines added)

### 3.5 Risks / edges

- Inspect `scan_raw_files()` first — if it does heavy per-file metadata (TIC extraction), defer that to post-pick to avoid wasted work on excluded files.
- TIC extraction observer wired off `values$diann_raw_files` change should naturally re-fire on the filtered list — confirm in test.

---

## 4. Phase C — Output dir Browse + Create-new-subdir

### 4.1 Two specific gaps

**C1 — Local container mode has no Browse button.** When `DELIMP_DATA_DIR` is set (containerized launchers), the output dir UI is just a `textInput` with no picker. Add a `shinyDirButton` next to it.

**C2 — Existing pickers don't let you create subdirectories on the fly.** User picks a parent → can only place output in pre-existing dirs.

### 4.2 Proposed flow

For all three modes (Local-native, Local-container, Docker):
- `shinyDirButton` to pick parent (existing or new for Local-container)
- `textInput("output_subdir")` next to the picker — optional "subdir name to create inside parent"
- "Create" button: validates subdir name (see Phase D rules), runs `dir.create(file.path(parent, subdir), recursive = TRUE)`, updates the displayed path
- On submit, if subdir typed but Create not clicked, auto-create at submit time (backstop)

### 4.3 Files touched (C)

- `R/ui.R:850-880` — Local container output (add Browse), all three modes (add subdir input + Create button)
- `R/server_search.R` — handlers for the new Create buttons + dir.create logic

### 4.4 Risks / edges

- Permissions: `dir.create()` may fail (read-only mount, disk full). Surface in notification.
- Path normalization: trailing slashes, `~` expansion, double-slashes — handle via `normalizePath()`.
- SSH mode untouched (preserves auto-derivation).

---

## 5. Phase D — Filename validation + auto-rename

### 5.1 Why this matters

Spaces and shell-special characters in filenames silently break:
1. **DIA-NN command line**: quoting layers (R → bash → SLURM → DIA-NN) break, files-not-found
2. **Output paths from `analysis_name`**: spaces → broken Docker container names (existing CLAUDE.md gotcha) and broken sbatch `mkdir -p` chains

Detection + warning catches a class of silent failures at submit-time instead of mid-search.

### 5.1b Scope: files AND folders

Validation + auto-rename applies to **both files AND directories**. Common cases:
- **`.d` directories** (Bruker timsTOF data) — these are dirs, not files. Same rules.
- **Output parent dirs** picked via shinyDirButton — if `~/My Search Results/` has a space, the path leaks into sbatch scripts and DIA-NN args. Detect at pick-time, offer rename.
- **Raw data parent dirs** — `~/Sample Set 1/` with raw files inside. Same issue.
- **FASTA dirs** — same.
- **Subdirs typed in Phase C's "Create new subdir" field** — validate before creating.

The auto-rename modal preview lists files AND dirs together, with a column noting type (`file`/`directory`/`.d directory`).

For nested cases (e.g. `~/Sample Set 1/sample (rep 1).raw`): rename the parent dir first, then the file. Modal shows the dependency in order.

### 5.2 Sanitization rules (matches FragPipe / msconvert conventions)

| Character | Replaced with |
|---|---|
| spaces | `_` |
| `()` `[]` | removed |
| single/double quotes | removed |
| `&` `\|` `;` `$` `*` `?` `<` `>` | `_` |
| non-ASCII | closest ASCII via `iconv` if possible, else `_` |
| leading `-` | `n_` prepended (avoids being parsed as a flag) |
| consecutive `_` | collapsed to single `_` |
| trailing `_` before extension | trimmed |

Examples:
- `Sample (rep 1) - 200ng.raw` → `Sample_rep_1_200ng.raw`
- `R&D test_file.raw` → `R_D_test_file.raw`
- `-control.raw` → `n_control.raw`
- `My Test Run #1` (analysis_name) → `My_Test_Run_1`

### 5.3 Validation flow

A new helper `validate_filename_for_shell(path)` returns a list:
```r
list(
  ok = logical(1),       # FALSE if any issues
  issues = character(),   # human-readable list ("contains space", "leading dash", ...)
  proposed = character(1) # the sanitized name
)
```

Called at:
- **FASTA picker confirm** (`fasta_browse_confirm`, `ssh_fasta_browse_confirm`)
- **Raw picker confirm** (Phase B's new confirm observers)
- **`analysis_name` field** — live validation on `input$analysis_name` change, surface a small inline warning + "use sanitized" link
- **`raw_data_dir` shinyDirButton** auto-load path (when files load without going through the picker)

### 5.4 Auto-rename flow

When validator detects ≥1 problem files, the warning banner gains an **"Auto-rename"** button.

Click flow:
1. Modal opens with a preview table:

   | Old name | New name |
   |---|---|
   | `Sample (rep 1).raw` | `Sample_rep_1.raw` |
   | `R&D test.raw` | `R_D_test.raw` |

2. Buttons: **Rename Files** / **Cancel**
3. On confirm:
   - For each file: `file.rename(old, new)` (local) or `ssh_exec("mv ...")` (SSH)
   - **Skip with warning if target name already exists**
   - **Roll back previous renames if any single rename fails** (best-effort)
   - Save `rename_log_<timestamp>.csv` to `~/.delimp_rename_logs/` so users can audit
4. Refresh the file list to show new names
5. Notification: "Renamed N files. Log saved to ~/.delimp_rename_logs/<file>.csv"

### 5.5 Special handling for Bruker `.d` directories

`.d` is a **directory**, not a file. Rules:
- Rename the directory itself (e.g. `Sample (1).d` → `Sample_1.d`)
- **Do NOT touch internal files** — `.d` directories have specific internal filenames (`analysis.tdf`, `analysis.tdf_bin`, etc.) that Bruker tooling expects. Renaming any of those would corrupt the data.

### 5.6 Files touched (D)

- `R/helpers_search.R` — new `validate_filename_for_shell()` and `sanitize_filename()` helpers
- `R/server_search.R` — observers for FASTA + raw confirm; analysis_name watcher; auto-rename modal + handler
- `R/ui.R` — warning-banner `uiOutput` slots; small inline warnings on textInputs

### 5.7 Risks / edges

- **Cross-platform file locking**: Windows holds locks on files in use; rename may fail with `EBUSY`. Surface clearly.
- **SMB / NFS mounts**: rename may not be atomic; `cp + rm` fallback if `mv` fails.
- **Permissions**: read-only filesystem → fail loudly.
- **User accidentally renames** then realizes they wanted originals back: the rename log CSV is the audit trail; reverse-rename is a future tool.
- **`.d` rename atomicity**: directory rename is atomic on the same filesystem; cross-filesystem requires copy. We'd error on cross-fs `.d` rename rather than half-rename.

---

## 6. Phase E — Local-search FIFO queue

### 6.1 Problem

Today, submitting a second local DIA-NN search while the first is still running produces undefined behavior:
- Both compete for CPU/RAM (slowing both, possibly OOM-killing one)
- Docker mode: the second `docker run` may fail with "container name already in use" if the analysis_name collides
- The job queue UI shows both as `running` but only one is making real progress

Users naturally want **submit-and-forget** for multiple searches: queue them, run sequentially, get notified as each completes.

### 6.2 Proposed shape

A `local_search_queue` reactive (or list in `values$diann_jobs` filtered to `backend == "local" | backend == "docker"` and `status == "queued" | "running"`).

**State machine per local search:**
| Status | Meaning | Transition |
|---|---|---|
| `queued` | Submitted, waiting for an open slot | When concurrency permits → `running` |
| `running` | Active subprocess (one slot occupied) | Subprocess exits → `completed` / `failed` / `cancelled` |
| `completed` / `failed` / `cancelled` | Terminal state | Removed from queue; next `queued` job promoted |

**Concurrency**: strict FIFO, **always 1 at a time**. Local boxes don't have the resources to run two DIA-NN searches concurrently — a single search already saturates 8-16 threads + 32-64 GB of RAM. Running two in parallel produces no speedup (resource-bound either way) and increases OOM risk. No concurrency-cap site config is exposed; if a power user with a workstation truly wants parallelism, they can run two DE-LIMP instances on different ports.

**Submit flow:**
1. User clicks Submit
2. `submit_diann` handler creates the job entry with `status = "queued"`
3. A new `observe()` reactive watches the queue: if a slot is open and ≥1 entry is queued, promote the oldest queued entry to `running` via the existing `run_local_diann()` / `run_docker_diann()` path
4. On subprocess exit (already monitored in the existing job-status observer), move to terminal state and re-check for queued entries

**UI changes:**
- Job Queue panel shows queue position for queued jobs ("Queued (#2 of 3)")
- "Reorder" buttons (up/down) to bump priority — v2; v1 is strict FIFO
- "Cancel" button on queued jobs removes them without launching
- Optional: estimated wait time based on average prior search duration

### 6.3 What about multi-step parallel searches (HPC mode)?

HPC parallel pipelines already get explicit SLURM-side queueing — the cluster handles concurrency. Phase E is **local-only**.

The existing `auto_load` flow on completion (already implemented for HPC) is the model: when a job's subprocess exits cleanly, fire the next-action observer. We're just adding the "promote next queued local job" trigger to that same exit observer.

### 6.4 Files touched (E)

- `R/server_search.R` — `submit_diann` handler (set status `queued` instead of immediately launching for local/docker), new queue-promotion observer
- `R/helpers_search.R` — concurrency cap reading from site config
- `R/ui.R` — queue position + cancel-queued button in the existing Job Queue panel

### 6.5 Speclib cache reuse (and contribute)

**Already in place**: `server_search.R:4430-4468` runs **before** the backend dispatch — every search (HPC, local, Docker) checks the predicted-library cache via `speclib_cache_lookup(fasta + params + mode)`. If a matching library exists, the search skips library prediction and reuses it. This means the local queue's first search may take 30-60 minutes (predicting library), but subsequent queued searches with the **same FASTA + same search params** reuse the library and skip directly to raw-file analysis. Single-FASTA labs running 50 batched searches save hours.

**Gap for local-only users**: cache **registration** (`speclib_cache_register`) is currently gated on the HPC parallel pipeline's `step_status[["step1"]] == "completed"` (`server_search.R:5566-5592`). Local/Docker single-job searches **never register** the libraries they generate. So a lab that only runs local searches has an empty cache forever, and Phase E's queue can't take advantage of cache reuse.

**Phase E adds**: a local-completion observer that registers the generated library on successful local/Docker exit. The library file is at `<output_dir>/report-lib.predicted.speclib` (already written by DIA-NN's `--out-lib`). Detect file existence on subprocess exit, call `speclib_cache_register()` with the same `(fasta_files, search_params, search_mode, ...)` tuple, mark the job entry's `speclib_cached = TRUE`. After the first search in the queue completes, every subsequent queued search with the same FASTA reuses the library.

**Cache key components** (already defined in `speclib_cache_key()` line 3430): FASTA file hashes (or filenames + sizes) + `enzyme` + `missed_cleavages` + peptide length range + precursor m/z range + charge range + variable mods. Different FASTA = different key = no false reuse.

**Cross-machine note**: when `speclib_cache_path()` resolves to a shared location (UCD's `/Volumes/proteomics-grp/dia-nn/`), local-search registration contributes to the lab-wide cache. Non-UCD users get a per-machine cache by default.

### 6.6 Risks / edges

- **Mid-flight cancel**: cancelling a `running` local job needs to kill the subprocess and free the slot. `processx::process` already has `$kill()` — wire it to the cancel button.
- **App restart with queued jobs**: queued jobs are persisted to `~/.delimp_job_queue.rds` (v3.10.20+). On restart, queued items should resume in order. `running` items should be retroactively reconciled — if the subprocess died when the app was killed, transition to `failed` with a notification.
- **Concurrent submits from rapid clicks**: debounce the Submit button (disable until submission is recorded) to avoid race conditions in the queue.
- **Disk-space check for the queue**: if 5 queued searches each need 50 GB, may run out of space midway. Pre-flight disk check shows a warning but doesn't block.

---

## 7. Phased rollout

| Phase | Scope | Effort | Ship as |
|---|---|---|---|
| **A** | Find-output UX (Open Folder + better default + host-path display) | ~1h | v3.11.0a |
| **B** | Per-file raw picker | ~1h | v3.11.0b |
| **C** | Output dir Browse (Local-container) + Create-new-subdir affordance | ~1.5h | v3.11.0c |
| **D** | Filename validation (files AND folders) + auto-rename + analysis_name sanitizer | ~2h | v3.11.0d |
| **E** | Local-search FIFO queue with concurrency cap | ~2h | v3.11.0e |

Total: ~7.5 hours across five small PRs, each independent and reviewable.

**Recommended order**:
1. **A** — highest payoff for the "I lost my output" complaint
2. **D** — prevents a class of silent failures (closest cousin to today's stability work)
3. **E** — quality-of-life: submit a batch, walk away
4. **B** — UX polish
5. **C** — UX polish

---

## 8. Open questions

1. **Default `output_root`**: `~/Documents/DE-LIMP/searches/`? Or `~/DE-LIMP/output/`? Or something else? Whatever it is, env var override via `DELIMP_OUTPUT_ROOT`.
2. **Auto-rename log location**: `~/.delimp_rename_logs/<timestamp>.csv` (separate dir) or written into the search output dir? I'd lean separate dir — survives if the user deletes the search output.
3. **Should auto-rename also handle `analysis_name`** (sanitizing the typed string), or just file paths? Probably yes — gives users a one-click "fix it for me" for the analysis name field too.
4. **Phase ordering** — A first, D next? Or all four concurrently?

---

## 9. Pre-implementation checklist

Before any code lands:

- [ ] Brett confirms phase order
- [ ] Brett picks default `output_root` (Q1)
- [ ] Brett confirms rename-log location (Q2)
- [ ] Inspect `scan_raw_files()` to confirm metadata-extraction cost during scan (Phase B risk)
- [ ] Confirm `values$diann_raw_files` is the single source of truth (no parallel state)
- [ ] Run a known-good search end-to-end on current code as a baseline
- [ ] Implement Phase A first, ship, stress-test
- [ ] Then D (validation only — auto-rename can be a follow-up if scope is too big)
- [ ] Then B
- [ ] Then C

---

## 10. Decision log

- **2026-05-07**: design drafted by Claude after Brett's UX request post-v3.10.29 install-stack stabilization. Brett asked to audit before coding given today's 26-hotfix run.
- **2026-05-07**: Brett confirmed scope is **Local + Docker only** — SSH HPC mode stays auto-derived. Real complaint is "users can't find their output," so Phase A added as highest-priority.
- **2026-05-07**: Brett confirmed filename validation should detect spaces/special chars; auto-rename feature added as Phase D modeled on FragPipe's rename tool — **non-destructive preview-then-confirm flow with rename log for audit**.
- **2026-05-07**: Brett added scope: validation also applies to **folder names** (not just files) — `.d` directories, parent dirs picked via shinyDirButton, subdirs typed in the Phase C "Create new subdir" field. Same sanitizer.
- **2026-05-07**: Brett added Phase E — local-search FIFO queue. Today, submitting a second local search while the first runs produces undefined behavior (concurrent processes fight for resources, Docker name collisions). Want submit-and-forget batch behavior with strict FIFO by default; concurrency cap configurable via site config for power users.
- **2026-05-07**: Brett added speclib cache reuse to Phase E — "like we do on the cluster." Cache **lookup** is already backend-agnostic (line 4430), so local searches automatically reuse cached libraries from anywhere in the lab. But cache **registration** is currently HPC-parallel-only — local searches never contribute to the cache. Phase E adds a local-completion observer that calls `speclib_cache_register()` after the search exits successfully, so a lab batching local searches gets cache reuse from the second search onward.
