# Output packaging, re-analysis & comparison

Every run is packaged into a tidy **session directory** so people can find things.

## Where the session goes — ask the user
The orchestrator asks where results should live (SKILL.md step 3b):
- **Default (recommended): in the folder with the raw data** being analyzed. The
  session folder is created right inside that directory, so results sit next to the
  files they came from. Pass `--raw <globs>` and omit `--base`.
- **A central location** the user prefers (e.g. `~/Documents/DataAnalysis`, or any
  path they give): pass `--base <path>` → the session is created under
  `<path>/sessions/`.

`session.py init` reports `placement` (`with-raw-data` | `central` | `reanalysis`).

## Session layout (`session.py`)
```
<YYYY-MM-DD>_<DescriptiveName>/    # inside the raw-data folder, or under <base>/sessions/
  README.html               # OPEN THIS: summary, links, where it all lives on HIVE (finalize)
  README.md                 # the same content as text (one source: session_docs.py)
  AGENTS.md                 # a guide to the folder for an AI agent, from the session's records
  MANIFEST.txt              # every finalize part as [OK] / [SKIPPED] <reason> / [INFO]
  input/                    # conditions.csv, search.fasta, params.*, wf/workflow.manifest.json,
                            #   raw_files.txt (raw data is referenced, NOT copied — too large),
                            #   submission.json + samples.tsv (Core data: the CoreOmics record,
                            #   allowlisted — submission_report.py attach)
  session.json              # session metadata; `coreomics` names the submission (Core data)
  output/
    search/                 # the normalized search report.parquet (+ search_provenance.json, logs)
    tables/                 # DE_*.csv, Expression_Matrix.csv, methods.txt, sessionInfo.txt,
                            #   de_provenance.json, and reproducibility_log.R — the analysis as
                            #   plain R, runnable with just R + limpa/limma (point users here
                            #   when they ask for "the code")
    figures/                # volcano / PCA / heatmap / p-value / QC PNGs + figures.json
    reproducibility/        # the pinned bundle (reproduce.sh, env lock, sessionInfo, skill.txt, checksums)
    AI_Analysis_Report.md   # the biological interpretation, with figures (read first)
    Analysis_Report.html    # THE report of record: one self-contained page (no Word copy)
    AUDIT.md                # results audit — common-mistake checks (PASS/WARN/FAIL)
    OUTPUT_FILES.md         # catalog of every file
    comparison/             # (re-analyses) COMPARISON.md + concordance CSVs
  scripts/                  # a copy of the skill scripts that ran this analysis (self-contained)
  logs/                     # commands.log + engine logs
```

- `session.py init --name "..." --raw <globs> [--base <path>] [--reanalysis-of <prior>]`
  makes the folders and prints a `paths` map; **route every step's
  `--out`/`--outdir`/`--dest` into those paths**.
- `session.py finalize --dir <session> [--zip]` writes `README.md` + `README.html` +
  `AGENTS.md`, moves any loose tables/figures into their subdirs, and (for a re-analysis)
  writes `DIFFERENCES.md`. `session.py docs --dir <session> [--as <real location>]` writes only
  the three documents (e.g. for a session finalized before they existed, or a copy of one).

### README.html, README.md and AGENTS.md (`session_docs.py`)
- **README.html** is what collaborators open: one self-contained page (inline CSS, no external
  assets) rendered from the same text as README.md by make_deposit's Markdown renderer, with links
  to the analysis report, the Word files, `HOW_TO_SUBMIT.html` and the tables. A double-clicked
  page opened from inside a zip loses its links (Windows extracts only that file): unzip first.
- **AGENTS.md** is for an AI agent given the folder: the study (organism, groups, contrasts,
  instrument, engine + version), which file is authoritative for what, the columns of `DE_*.csv`
  / `Expression_Matrix.csv` / `QC_detected_vs_inferred.csv`, the traps (significance as
  `de_provenance.json` records it, inferred vs measured values, contaminants, pull-down controls
  by group name), the AUDIT / SAMPLE_QUALITY notes, how to reproduce, and what not to do. The
  pipeline is described in `de_provenance.json`'s own words -- never a description written here.
- **"Where this lives on HIVE"** (both): the session folder, the session the analysis ran in when
  this is a copy, the raw data folder(s) + file count, the search output, the FASTA the search
  read, the Core run-registry record and the FRAN hand-off entry -- each as its HIVE path and its
  Windows (`\\128.120.208.24\proteomics\...`) and Mac (`/Volumes/proteomics/...`) equivalent.
  Which share is which HIVE path is `scripts/hive_shares.tsv`, the one table `hive_path.sh` also
  reads (`share_map.py` is the Python side). Anything the records do not give is "not recorded".
- **input/raw_files.txt** is written at finalize when missing (a hive_remote session is initialised
  without `--raw`), from `search_provenance.json` `files`, else `output/search/file_list.txt`.
- The registry record is looked up read-only before the zip (`record_run.locate()`); when the
  run-log hook creates it, the three documents are rewritten with its path and go into the zip
  after the hook, like `MANIFEST.txt`.

## Re-analysis of the same dataset
Re-running the same raw data (different engine, version, parameters, FASTA, or
design) is common and must not clobber or be confused with the original.

- **Detection:** `session.py find-prior --raw <globs>` scans existing sessions'
  `input/raw_files.txt` for an overlapping raw-file set and reports matches
  (`same_dataset` true when the sets are identical).
- **Placement:** with `--reanalysis-of <prior>`, the new run nests under
  `<prior>/reanalysis/<date>_<name>/` — same internal layout — so all re-analyses
  live with their original.
- **`DIFFERENCES.md`** (written at finalize) states exactly what changed vs the
  original: engine + version, DE method + thresholds, contrasts, FASTA, the pinned
  workflow commit, a unified diff of the search parameters, and the
  significant-protein counts per contrast. Unchanged settings are omitted.

## Comparing analyses (`compare_analyses.R`)
`DIFFERENCES.md` says what *settings* changed; the Comparator shows how the
*results* changed. `compare_analyses.R` is a faithful port of DE-LIMP's Run
Comparator core (`normalize_protein_id`, `classify_de`, the 3×3 concordance) so the
skill stays self-contained (it can't depend on the DE-LIMP repo being present).

For each shared contrast across ≥2 analyses it reports:
- **protein-universe overlap** (proteins found + significant per analysis),
- the **3×3 Up/Down/NS concordance** matrix on shared proteins,
- **direction concordance** on co-significant proteins, and
- **logFC correlation** on shared proteins.

Outputs: `COMPARISON.md`, `concordance_summary.csv`, and per-pair 3×3 +
merged-protein CSVs. Use it whenever two analyses of the same dataset exist
(re-analysis vs original, or two engines/parameter sets side by side).
