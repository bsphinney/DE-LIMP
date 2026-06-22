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
  README.md                 # what this was + where everything is (written at finalize)
  input/                    # conditions.csv, search.fasta, params.*, wf/workflow.manifest.json,
                            #   raw_files.txt (raw data is referenced, NOT copied — too large)
  output/
    search/                 # the normalized search report.parquet (+ engine logs)
    tables/                 # DE_*.csv, methods.txt, sessionInfo.txt, de_provenance.json, QC
    figures/                # plots (reserved)
    reproducibility/        # the full bundle (reproduce.sh, env lock, checksums)
    AI_Analysis_Report.md   # the biological interpretation (read first)
    AI_Analysis_Report.docx # the same report as a Word document
    OUTPUT_FILES.md         # catalog of every file
    comparison/             # (re-analyses) COMPARISON.md + concordance CSVs
  scripts/                  # a copy of the skill scripts that ran this analysis (self-contained)
  logs/                     # commands.log + engine logs
```

- `session.py init --name "..." --raw <globs> [--base <path>] [--reanalysis-of <prior>]`
  makes the folders and prints a `paths` map; **route every step's
  `--out`/`--outdir`/`--dest` into those paths**.
- `session.py finalize --dir <session> [--zip]` writes `README.md`, moves any loose
  tables/figures into their subdirs, and (for a re-analysis) writes `DIFFERENCES.md`.

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
