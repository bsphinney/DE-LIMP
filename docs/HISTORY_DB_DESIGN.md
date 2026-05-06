# DE-LIMP History Database — Design Doc (v3.11.0)

Status: **draft for review, 2026-05-06**
Author: Claude (with Brett)
Supersedes: the activity-log CSV at `~/.delimp_activity_log.csv`

---

## 1. Why we're rebuilding history

The current activity log (`~/.delimp_activity_log.csv`, 33 columns, append-only) has been the source of repeated friction across the v3.10.x hotfix run:

- **Local-only**: Brett's laptop, the Docker container, and the HF Space each have a different copy. No way to find a search submitted from another machine.
- **Fragile**: a single corrupt row breaks parsing for everyone (CSV with no schema enforcement).
- **Substep leakage**: parallel-pipeline phase substeps (`_s1_libpred`, `_s2_firstpass`, etc.) ended up as separate "search" entries from the broken Recover handler. The CSV had no concept of "logical search" so the duplicate rows persisted.
- **Doesn't scale**: 33 columns with most NULL most of the time; can't index for fast lookup; no JOIN with samples or proteins.
- **Duplicates source-of-truth**: `search_info.md` next to each parquet already has every search parameter and instrument detail. The CSV duplicates them and can drift out of sync.

The primary user need we're solving for: **"I want to find a past search by sample name (or gene, or date, or submitter) and one-click load it back into DE-LIMP."** Bonus: a lab-wide proteome aggregate ("which proteins has our lab observed across all 5,000 of our past searches?") falls out for free once the per-search protein detection table exists.

---

## 2. Architectural principles

1. **`search_info.md` next to the parquet is the canonical record.** Anything we can re-derive from there does not get duplicated.
2. **The DB is a searchable index, not a source of truth.** If `~/.delimp_history.db` is deleted, running "Discover" against the data dirs rebuilds it lossless.
3. **No app state on mounted drives.** The DB lives at `~/.delimp_history.db` (local). The DATA being indexed lives wherever (mount, SSH, etc.).
4. **Lab-wide via shared filesystem walk, not via centralized server.** Every DE-LIMP instance discovers the same `SERVICE/` root and ends up with the same picture in its local DB. No sync server, no shared SQLite.
5. **Every silent failure becomes visible.** Discovery reports per-section status; ingest errors are logged; missing files are recorded but never crash the walk. (Hard-learned lesson from the v3.10.4–v3.10.13 silent-failure-masquerading-as-success pattern.)

---

## 3. Schema

```sql
-- One row per logical search. Keyed on output_dir.
CREATE TABLE searches (
  id              INTEGER PRIMARY KEY,
  search_name     TEXT NOT NULL,
  output_dir      TEXT UNIQUE NOT NULL,    -- canonical key — a search "is" its output dir
  submitted_at    INTEGER,                  -- epoch seconds
  completed_at    INTEGER,
  status          TEXT,                     -- queued/running/completed/failed/cancelled
  duration_min    REAL,
  job_id          TEXT,                     -- SLURM array parent ID (for parallel pipelines)
  user            TEXT,                     -- submitter (from search_info.md "Submitted by" or file owner)
  project         TEXT,
  notes           TEXT,
  n_files         INTEGER,
  fasta           TEXT,                     -- comma-joined basenames
  instrument      TEXT,                     -- "timsTOF HT" etc.
  acquisition     TEXT,                     -- "dia-PASEF" etc.
  pipeline_id     TEXT,                     -- dpc | maxlfq (set when analyzed, NULL until then)
  app_version     TEXT,                     -- DE-LIMP version that ingested this row
  parquet_path    TEXT,                     -- canonical local path (or HPC if local missing)
  last_seen       INTEGER,                  -- when discover last saw this dir
  source          TEXT                      -- submit | discover | recover | load_from_hpc
);
CREATE INDEX idx_s_user       ON searches(user);
CREATE INDEX idx_s_project    ON searches(project);
CREATE INDEX idx_s_status     ON searches(status);
CREATE INDEX idx_s_submitted  ON searches(submitted_at DESC);
CREATE INDEX idx_s_instrument ON searches(instrument);

-- One row per raw file in a search.
CREATE TABLE samples (
  id            INTEGER PRIMARY KEY,
  search_id     INTEGER REFERENCES searches(id) ON DELETE CASCADE,
  file_name     TEXT NOT NULL,             -- basename of the raw file
  group_label   TEXT,                       -- group assigned at analysis time (NULL until then)
  precursors    INTEGER,                    -- from report.stats.tsv
  proteins      INTEGER,
  ms1_signal    REAL,
  rt_range      TEXT,                       -- "0-120 min"
  mass_acc_ms2  REAL,                       -- median ppm error
  detected_pct  REAL                        -- % of total precursors detected
);
CREATE INDEX idx_samples_filename ON samples(file_name);
CREATE INDEX idx_samples_search   ON samples(search_id);

-- Pointers to canonical files (parquet, pg_matrix, log, search_info, session.rds).
-- Used by the "Load this search" path-resolution logic.
CREATE TABLE search_files (
  id          INTEGER PRIMARY KEY,
  search_id   INTEGER REFERENCES searches(id) ON DELETE CASCADE,
  kind        TEXT,                         -- parquet | pg_matrix | stats | log | search_info | empirical_lib | session_rds
  path        TEXT NOT NULL,
  size_bytes  INTEGER,
  mtime       INTEGER,
  is_remote   INTEGER DEFAULT 0             -- 0 = local resolvable, 1 = HPC-only (need SSH)
);
CREATE UNIQUE INDEX idx_files_search_kind ON search_files(search_id, kind);

-- Per-search per-protein detection summary. Keyed on (search_id, protein_group).
-- Powers the lab-wide proteome aggregate.
CREATE TABLE search_proteins (
  search_id           INTEGER REFERENCES searches(id) ON DELETE CASCADE,
  protein_group       TEXT NOT NULL,
  gene                TEXT,                 -- gene symbol when resolvable
  protein_names       TEXT,                  -- description
  n_samples_detected  INTEGER,               -- in this search, how many runs detected it
  avg_log2_intensity  REAL,
  is_contaminant      INTEGER DEFAULT 0,
  PRIMARY KEY (search_id, protein_group)
);
CREATE INDEX idx_sp_gene ON search_proteins(gene);
CREATE INDEX idx_sp_pg   ON search_proteins(protein_group);

-- Materialized lab-wide aggregate. Refreshed incrementally during ingest.
CREATE TABLE lab_proteome (
  identifier             TEXT PRIMARY KEY, -- gene symbol if available, else protein_group
  is_gene                INTEGER,           -- 1 if `identifier` is a real gene symbol
  protein_group_example  TEXT,              -- canonical example
  protein_names          TEXT,
  n_searches             INTEGER,           -- in how many searches detected
  n_samples_total        INTEGER,           -- summed across searches
  first_seen             INTEGER,           -- epoch of earliest search
  last_seen              INTEGER,           -- epoch of most recent
  median_log2_intensity  REAL,
  pct_contaminant        REAL               -- 0..1
);
CREATE INDEX idx_lp_n_searches ON lab_proteome(n_searches DESC);
CREATE INDEX idx_lp_last_seen  ON lab_proteome(last_seen DESC);

-- Free-form events log (per-search analyses, session restores, etc.).
-- Replaces the catch-all CSV; users mostly won't query this directly.
CREATE TABLE events (
  id          INTEGER PRIMARY KEY,
  ts          INTEGER NOT NULL,
  event_type  TEXT NOT NULL,                -- search_submitted | search_completed | analysis_completed | session_restored | discover_run
  search_id   INTEGER,                       -- NULL for discover_run, etc.
  details     TEXT                           -- JSON catch-all
);
CREATE INDEX idx_events_ts   ON events(ts DESC);
CREATE INDEX idx_events_type ON events(event_type);

-- Schema version tracking (for migrations).
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
INSERT INTO meta VALUES ('schema_version', '1');
INSERT INTO meta VALUES ('created_at',     CAST(strftime('%s','now') AS TEXT));
```

### Storage budget

- Per search: ~5,000 protein rows × ~80 bytes + ~50 sample rows × ~80 bytes + 1 search row + ~6 file rows ≈ **400–500 KB/search**
- 1,000 searches: ~500 MB
- 10,000 searches: ~5 GB
- 50,000 searches (50 lab members × 200/year × 5 years): ~20–25 GB

SQLite handles this without breaking a sweat. If we ever need to trim, we add an "archive searches older than N years" knob that moves rows to `~/.delimp_history_archive.db`.

---

## 4. Discovery (the "find search results" part)

### Walk strategy

```
discover(roots = c("/Volumes/proteomics-grp/SERVICE", ...)):
  for each root:
    if root is unmounted/inaccessible:
      if SSH connected: walk via `ssh ... find /quobyte/.../SERVICE -name search_info.md`
      else: skip with [SKIPPED] reason
    else: walk via list.files(recursive=TRUE, pattern="search_info.md")
  
  for each search_info.md found:
    output_dir = dirname(search_info_path)
    mtime     = file.info(search_info_path)$mtime
    if (output_dir, mtime) already in DB and unchanged: skip [INCREMENTAL]
    parse search_info.md -> (search_params, instrument, fasta, ...)
    upsert into `searches` keyed on output_dir
    parse report.stats.tsv if present -> per-sample rows in `samples`
    parse report.pg_matrix.tsv if present -> per-protein rows in `search_proteins`
    list output_dir -> insert paths into `search_files`
    update `lab_proteome` for the genes that appeared
  
  emit summary: "Found N new + M updated + K unchanged across L directories in 12s"
```

### User attribution

Three fallback tiers:
1. **`**Submitted by**: <username>` line in search_info.md** — preferred. Backwards-incompatible upgrade: `generate_search_info()` writes this going forward (one-line addition; v3.11.0a includes the change).
2. **File owner of search_info.md** via `file.info(...)$uname` — works for older searches on POSIX filesystems.
3. **Path-based heuristic** — last resort, parse `SERVICE/<group>/<submitter>/...` if the lab convention encodes submitter in the path. Stays NULL if no convention applies.

### Performance / scaling

- A walk of `SERVICE/` with 1,000 searches: ~10–20 seconds full ingest, ~1–2 seconds incremental (skipping unchanged dirs).
- Discovery runs in `withProgress(...)` so the user sees per-dir progress.
- DB writes happen inside a single transaction per search → atomic commit, no partial state on crash.

### Default roots

```r
# Configurable via Settings → History → Discovery roots.
# Or env var: DELIMP_DISCOVER_ROOTS="path1,path2,..."
default_roots <- c(
  "/Volumes/proteomics-grp/SERVICE",                          # Mac SMB mount
  "/quobyte/proteomics-grp/SERVICE",                           # HPC (when running as Apptainer)
  "C:/Users/<user>/Documents/dia-nn"                          # Windows fallback (TODO)
)
```

---

## 5. Robust "Load this search" — the primary UX

When the user clicks **Load** on a History row, the resolver tries paths in order:

```
load_search(search_id):
  files = SELECT * FROM search_files WHERE search_id = ?
  parquet_paths = files where kind == 'parquet'
  session_paths = files where kind == 'session_rds'
  
  # 1. Pick canonical resolved path for the parquet
  resolved_parquet = NULL
  for path_record in parquet_paths sorted by (is_remote ASC, mtime DESC):
    candidate = path_record$path
    
    # Try original
    if file.exists(candidate): resolved_parquet = candidate; break
    
    # Try translated local mount (HPC -> /Volumes)
    translated = translate_storage_path(candidate, to = "local")
    if file.exists(translated): resolved_parquet = translated; break
    
    # Try cached local copy
    cache_path = file.path("~/.delimp_cache", digest(output_dir), "report.parquet")
    if file.exists(cache_path): resolved_parquet = cache_path; break
    
    # Try SSH
    if path_record$is_remote && ssh_connected:
      tmp = tempfile(fileext=".parquet")
      if scp_download(cfg, candidate, tmp): resolved_parquet = tmp; break
  
  if resolved_parquet is NULL:
    showNotification("Couldn't find report.parquet anywhere. Tried: <list>. Connect SSH or mount the share.", type="error")
    return
  
  # 2. Pick session.rds with same logic if present, then prompt the user
  resolved_session = (same logic for session_paths)
  
  if resolved_session AND resolved_parquet:
    showModal(dialog(
      title = "Load mode",
      "session.rds detected — restore the full analysis state, or load just the parquet?",
      footer = c(
        actionButton("load_session", "Restore Full Session"),
        actionButton("load_parquet", "Load Parquet (Fresh Analysis)")
      )
    ))
  else if resolved_session:
    restore_session(resolved_session)
  else:
    load_parquet(resolved_parquet)
  
  # 3. Always populate values$diann_search_settings + values$instrument_metadata
  #    from the search_info.md path (if available in search_files).
  if a search_info path resolves: 
    values$diann_search_settings <- parse_search_info_md(...)
    values$instrument_metadata <- parsed$instrument_metadata
```

This is the path-resolution logic that v3.10.8 introduced for Complete Analysis exports, generalized and made the canonical Load mechanism. It also subsumes the current "Load from HPC" button — that just becomes a special case where the user clicked a row that has only HPC paths.

### Show-details modal

Click the row (not the Load button) → modal with:
- Parsed `search_info.md` rendered as HTML
- Sample list (joined `samples` table)
- File listing with sizes (joined `search_files` table)
- Tail of the DIA-NN log (live SSH read if remote)
- Status badge + computed duration
- "Recompute discovery for this dir" button (if the user added new files post-search)

---

## 6. Find UX — what the History tab gets

### Filter bar

A persistent header above the table:

| Element | Behavior |
|---|---|
| Free-text search | Single text box. Fuzzy across `search_name`, `project`, `notes`, `samples.file_name`, `search_proteins.gene`. Returns row if ANY match. |
| Quick pills | "My searches", "Lab searches", "Last week", "Last month", "Failed only", "DPC-Quant only", "MaxLFQ only", "With phospho" |
| Date range picker | Two date inputs |
| Instrument filter | Multi-select (timsTOF HT, Orbitrap Fusion, …) |
| FASTA filter | Multi-select |
| Submitter filter | Multi-select |
| Saved searches | Bookmark a filter combo with a name; reapply with one click |

### Row actions

Each row in the table has a button group:

| Button | Action |
|---|---|
| **Load** | The resolver above. Default action. |
| **Settings** | Populate the New Search panel with this search's params (re-run with the same setup) |
| **Folder** | Open `output_dir` in Finder / Explorer / `xdg-open` |
| **Log** | Modal with DIA-NN log content (tail) |
| **Notes** | Edit `searches.notes` inline |

### Lab Proteome tab (v3.11.2)

Separate tab, sourced from the `lab_proteome` aggregate. Sketch:

- **Top identified proteins** (data table): identifier, n_searches, n_samples_total, last_seen, median intensity. Sortable. Click a row → modal with all searches that detected it, intensity heatmap.
- **Search by gene** (text input): enter a gene symbol → list of searches. Clickable to Load.
- **Coverage panel**: pie / bar of `count(distinct gene) / known_proteome_size` per organism. "Our lab has observed 67% of the bovine proteome over the past 18 months."
- **Trends**: line chart of distinct proteins observed per month. Per instrument. Per submitter.
- **Compare searches** (pick two): Venn of overlapping protein sets, scatter of shared-protein log-intensities.

---

## 7. Migration from CSV

Run once at first startup of v3.11.0a:

```
migrate_csv_to_db():
  csv_path = "~/.delimp_activity_log.csv"
  if !file.exists(csv_path): return
  
  rows = read.csv(csv_path)
  
  # Filter out cruft from earlier broken Recovers
  rows = filter(rows, !grepl("_s[1-5]_[a-z]+$", rows$search_name))
  rows = filter(rows, !grepl("^[0-9]+_[0-9]+$", rows$job_id))
  
  # Group by output_dir to dedup substeps that survived the filter
  searches_grouped = group_by(rows, output_dir)
  
  for each group:
    upsert searches keyed on output_dir
    upsert events for each row in the group
  
  # Backup old CSV
  file.rename(csv_path, paste0(csv_path, ".pre_v3.11.bak"))
  
  showNotification(sprintf("Migrated %d searches from CSV to DB. Old CSV backed up.", n))
  log_event("migration", details=...)
```

**Migration is idempotent** — if the user re-runs it (e.g., re-installs DE-LIMP and finds an old CSV), it skips already-upserted entries.

---

## 8. Phased rollout

| Version | Scope | Time | UI changes? |
|---|---|---|---|
| **v3.11.0a** | DB schema + helpers + migration. Writes the DB in parallel; CSV remains primary. | ~1.5h | None |
| **v3.11.0b** | History tab queries DB. Filter bar, row actions including new Load resolver. | ~2h | History tab rewrite |
| **v3.11.0c** | Discover command exposed in UI; Settings panel for roots. | ~1h | Discover button + Settings |
| **v3.11.1** | `search_proteins` + `lab_proteome` populated during ingest. | ~1.5h | None |
| **v3.11.2** | New "Lab Proteome" tab. | ~2h | New tab |
| **v3.11.3** | Drop CSV writes; add archive-old-searches knob. | ~30m | None |

Each version is shippable on its own; the user gets value at every step. By v3.11.2 the lab GPMDB use case is real.

---

## 9. Known risks & open questions

### Risks

- **DB grows unbounded** — at 50K searches we're at ~20 GB. Mitigation: archive-after-N-years knob in v3.11.3; index on `submitted_at` to make the cutoff query fast.
- **Cold discovery is slow** — first walk of a 1,000-search SERVICE/ takes ~10–20 s. Mitigation: incremental ingest (mtime-keyed), background discovery on app startup with cached results.
- **Path translation edge cases** — what if SERVICE/ moves? The DB has stale `output_dir`s. Mitigation: discover updates `last_seen` and `parquet_path`; the resolver tries the translated path before failing.
- **Cross-machine drift** — Brett's laptop and the HF Space might run discover on different days and end up out of sync briefly. Acceptable: both eventually converge after the next discover; no centralized state to corrupt.
- **search_info.md format drift** — older versions might not have all the bullets. Mitigation: parser already handles missing keys gracefully; ingest never crashes on a missing field.

### Open questions for review

1. **User attribution lookup**: do you want path-based fallback (parse `SERVICE/<lab>/<submitter>/...`)? Some lab conventions encode it; others don't.
2. **Default discovery roots on Windows**: ?
3. **Should `lab_proteome` be lab-wide or per-user?** Earlier you said lab-wide. Confirming: `lab_proteome.n_searches` counts EVERY search in DB regardless of submitter. If you'd rather have a "my proteome" view alongside, we can add a parameterized query (filtered by user); both views from one table.
4. **Auto-discover cadence**: walk SERVICE/ on every app startup (slow if mount is slow), once per day, or only manual? Default = manual + a non-blocking background scan once per app launch with the previous results displayed immediately?
5. **What to do when the DB schema changes** (e.g., you ask for new columns in v3.11.5)? `meta.schema_version` lets us run incremental ALTER TABLE migrations on bump. Worth implementing the migration framework now in v3.11.0a.

---

## 10. Testing plan

For each PR:

- **Unit tests**: parser correctness (search_info.md → searches row); SQL upserts idempotent; load-resolver picks correct path tier.
- **Migration test**: feed a real `activity_log.csv.bak` through the migrator; confirm row counts match expectations and no substep cruft survives.
- **Discovery stress**: walk a directory with 100+ search_info.md files; measure ingest time; confirm second walk is incremental.
- **Robust-load test**: simulate (a) mount missing, (b) SSH disconnected, (c) cached copy only — Load button should pick the right path or surface a clear error.
- **Lab Proteome correctness**: ingest a known set of N searches with known protein detections; confirm `lab_proteome.n_searches` matches.

---

## 11. Decision log

- **2026-05-06**: design doc drafted by Claude after the v3.10.4–v3.10.13 hotfix run made the CSV's fragility undeniable. Brett's primary ask: "robust way to find past searches and load them into DE-LIMP." Bonus ask: store proteins for a lab-wide GPMDB-style aggregate.
- **2026-05-06**: lab-wide aggregation chosen over per-user. Walk SERVICE/ recursively as the default discovery root.
- **2026-05-06**: SQLite over filesystem-of-JSONs because we want indexed queries (find by gene, find by sample). SQLite is bundled with R, ACID, file-based, no server.
- **2026-05-06**: `search_info.md` next to data is canonical; DB is just an index. Hard rule from CLAUDE.md: "Derived data stays with source data."
