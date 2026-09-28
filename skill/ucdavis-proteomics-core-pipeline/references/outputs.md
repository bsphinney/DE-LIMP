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
    tables/                 # DE_<method>_<contrast>.csv (one per comparison), Expression_Matrix.csv,
                            #   Detection_Matrix.csv (per protein x run: DPC = precursors
                            #   observed, 0 = inferred; MaxLFQ = 1/0 quantified/missing; which
                            #   applies is de_provenance.json detection_matrix), QC_detected_vs_inferred.csv,
                            #   QC_contaminant_share.csv + contaminants_removed.csv (the Cont_
                            #   filter), methods.txt, sessionInfo.txt, de_provenance.json,
                            #   DE-LIMP_session.rds (DPC runs: load it in the DE-LIMP app), and
                            #   reproducibility_log.R — the analysis as plain R, runnable with
                            #   just R + limpa/limma (point users here when they ask for "the code")
                            #   Each DE table carries Detected_<group> (k/n measured runs per
                            #   compared group) and Evidence -- "The DE tables" below
    figures/                # volcano / top-protein violins / PCA / heatmap / p-value / QC PNGs
                            #   + figures.json (captions) + sample_labels.csv (the short
                            #   sample names on the plots -> their full run names)
    reproducibility/        # the pinned bundle (reproduce.sh, env lock, sessionInfo, skill.txt, checksums)
    AI_Analysis_Report.md   # the interpretation the agent wrote (step 9); an input to the report
    Analysis_Report.html    # THE report of record: one self-contained page (no Word copy)
    Analysis_Report.md      # its plain-text twin, captions and numbers written out (NotebookLM)
    Analysis_Report.pdf     # the same page printed by a headless Chrome/Chromium/Edge, when one
                            #   is installed (make_analysis_html.py, or finalize via html_to_pdf.py)
    methods.md / .docx      # publication Methods + the instrument-grant acknowledgment (step 9d;
                            #   finalize writes them when missing)
    AUDIT.md                # results audit — common-mistake checks (PASS/WARN/FAIL)
    SAMPLE_QUALITY.md       # biological sample-quality panels (contamination that mimics biology)
    DATA_SUBMISSION/        # PRIDE / MassIVE deposit package (finalize): HOW_TO_SUBMIT.md/.html,
                            #   sdrf.tsv, protocols.txt, files_to_upload.tsv, prepare_upload.sbatch
    podcast/                # OPTIONAL audio discussion (make_podcast.py; references/podcast.md):
                            #   podcast.m4a, podcast_script.md (with its claims ledger),
                            #   transcript.html, podcast.json, check.txt; .cache/ is scratch
    OUTPUT_FILES.md         # catalog of every file
    comparison/             # (re-analyses) COMPARISON.md + concordance CSVs
  scripts/                  # a copy of the skill scripts that ran this analysis (self-contained)
  logs/                     # commands.log + engine logs
```

- `session.py init --name "..." --raw <globs> [--base <path>] [--reanalysis-of <prior>]`
  makes the folders and prints a `paths` map; **route every step's
  `--out`/`--outdir`/`--dest` into those paths**.
- `session.py finalize --dir <session> [--zip]` ensures the Methods, writes the deposit package,
  `README.md` + `README.html` + `AGENTS.md` and `MANIFEST.txt`, prints the report PDF when it is
  missing and a browser is there, moves any loose tables/figures into their subdirs, (for a
  re-analysis) writes `DIFFERENCES.md`, and with `--zip` zips the session, then logs the run and
  posts to the Core's Slack channel (Core runs). `session.py docs --dir <session> [--as <real location>]` writes only
  the three documents (e.g. for a session finalized before they existed, or a copy of one).

### README.html, README.md and AGENTS.md (`session_docs.py`)
- **README.html** is what collaborators open: one self-contained page (inline CSS, no external
  assets) rendered from the same text as README.md by make_deposit's Markdown renderer, with links
  to the analysis report (HTML, and its PDF and Markdown twin when present), the Methods
  (`methods.md` / `.docx`), `HOW_TO_SUBMIT.html` and the tables. Its "Where everything is"
  table names `Detection_Matrix.csv`, `QC_detected_vs_inferred.csv`, `contaminants_removed.csv`,
  `AUDIT.md` and `SAMPLE_QUALITY.md` one by one, each in its own record's words where it has
  them. When there is no PDF it points at the "Report PDF" line of `MANIFEST.txt`. A
  double-clicked page opened from inside a zip loses its links (Windows extracts only that
  file): unzip first.
- **AGENTS.md** is for an AI agent given the folder: the study (organism, groups, contrasts,
  instrument, engine + version), which file is authoritative for what, the columns of `DE_*.csv`
  / `Expression_Matrix.csv` / `Detection_Matrix.csv` / `QC_detected_vs_inferred.csv`, the traps (significance as
  `de_provenance.json` records it, inferred vs measured values, contaminants, pull-down controls
  by group name), the AUDIT / SAMPLE_QUALITY notes, how to reproduce, and what not to do. The
  pipeline is described in `de_provenance.json`'s own words -- never a description written here.
- **"Where this lives on HIVE"** (both): the session folder, the session the analysis ran in when
  this is a copy, the raw data folder(s) + file count, the search output, the FASTA the search
  read, the Core run-registry record and the FRAN hand-off entry -- each as its HIVE path and its
  Windows (`\\128.120.208.24\proteomics\...`) and Mac (`/Volumes/proteomics/...`) equivalent.
  Which share is which HIVE path is `scripts/hive_shares.tsv`, the one table `hive_path.sh` also
  reads (`share_map.py` is the Python side). Anything the records do not give is "not recorded".
  A row on the Core's own storage (`/quobyte/proteomics-grp`) says so: only Core members can open
  it, so its Windows cell reads "Proteomics Core storage: only Core members can open it -- ask the
  Core for a copy or for access" (the `access` column of `hive_shares.tsv`) and its Mac path is
  marked "(Core members)". The Mac Connect-to-Server route `smb://128.120.208.24/proteomics` is
  given only for paths under `/Volumes/proteomics` (Flinders).
- README.md, README.html and AGENTS.md are each written through a `.part` file and renamed into
  place, so a render that fails leaves the previous file whole (never a 0-byte README.html), and
  finalize zips only the ones it wrote this time.
- **input/raw_files.txt** is written at finalize when missing (a hive_remote session is initialised
  without `--raw`), from `search_provenance.json` `files`, else `output/search/file_list.txt`.
- The registry record is looked up read-only before the zip (`record_run.locate()`); when the
  run-log hook creates it, the three documents are rewritten with its path and go into the zip
  after the hook, like `MANIFEST.txt`.

### The optional podcast (`output/podcast/`, `make_podcast.py`)
An AI-generated audio discussion of the finished report, made only when the user asks
(`references/podcast.md`). It is a **derivative, not a record**: README and AGENTS.md say so,
and the report and tables stay authoritative.
- `podcast.m4a` (or `podcast.wav` when there is no AAC encoder): the episode.
- `podcast_script.md`: the script, with its **"Claims beyond the report"** ledger.
- `transcript.html`: the disclosure, a player, the transcript, the claims and how it was made.
- `podcast.json`: show, title, hosts and voices, TTS backend and exact model, cloud consent,
  script sha256, sources with their sha256, words, duration, `ai_generated: true`.
- `check.txt`: the fidelity check the render was gated on.
- `verify.txt`, `verify_transcript.txt`: the ASR round trip (`make_podcast.py verify`) -- word
  match ratio, spans not heard, numbers not heard with transcript context, segments to listen
  to; its result is also the `verify` block in `podcast.json`.

`make_podcast.py link <session>/output` adds a "Listen" card near the top of
`Analysis_Report.html` (hidden when printed, which prints the file name instead), a line under
the Markdown report's title, and entries in README and AGENTS.md. Each addition sits between
`<!-- podcast:start -->` / `<!-- podcast:end -->` and is replaced on a re-run.
`make_analysis_html.py` and `session_docs.py` add the same card and entries on their own when
`podcast/podcast.json` exists, so regenerating the report or re-finalizing keeps them. An
`Analysis_Report.pdf` older than the edited HTML is reprinted by `link` (html_to_pdf.py), or
flagged with an `[INFO]` line saying how to reprint it.
`podcast/.cache/` (per-chunk TTS audio, ~60 MB for 20 min) is kept on disk for resuming.
`scripts/scratch_files.py` is the one rule for scratch: every `.cache` folder, every `*.part`
file, and `podcast.wav` beside `podcast.m4a`. The session zip leaves them out (`zip_excluded`)
and `make_report.py` does not list them. Never deposit them. The run-registry record
(`record_run.py`) copies what the Listen card points at: the audio `podcast.json` names,
`transcript.html`, `podcast_script.md`, `check.txt` and `podcast.json`. A `podcast.json` that is
unreadable or has a wrong field never stops the report: it is reported as a `[WARN]` and the
report, README and AGENTS.md are made without the podcast.

## Re-analysis of the same dataset
Re-running the same raw data (different engine, version, parameters, FASTA, or
design) is common and must not clobber or be confused with the original.

- **Detection:** `session.py find-prior --raw <globs>` scans existing sessions'
  `input/raw_files.txt` for an overlapping raw-file set and reports matches
  (`same_dataset` true when the sets are identical).
- **Placement:** with `--reanalysis-of <prior>`, the new run nests under
  `<prior>/reanalysis/<date>_<name>/` — same internal layout — so all re-analyses
  live with their original.
- **`DIFFERENCES.md`** (written at finalize) states exactly what changed vs the original:
  engine + version, DE method, q / logFC value and role / adj.P, contrasts, the design with its
  covariates, the blocking factor (column, effect, scope), the DE contaminant policy, the FASTA
  and its sidecar state (`fetch_fasta.sidecar_state`), the skill and limpa versions, the raw-file
  list, a unified diff of the search parameters (the resolved cfg, else `input/wf/params.cfg`),
  and the significant-protein counts per contrast. Unchanged settings are omitted. A setting
  neither session records is listed as "not compared", never counted as unchanged.

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

`session.py init --date YYYY-MM-DD` dates the session folder (default: today) — for a session
created after the day the analysis actually ran.
