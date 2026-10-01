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
                            #   compared group) and Evidence -- "The DE tables" below.
                            #   set_test_inputs.rds: the DE model's inputs for run_sets.R.
                            #   Protein-set tests (run_sets.R, when run): Sets_<contrast>.csv,
                            #   Sets_baitnorm_<contrast>.csv (pulldowns), sets_members.csv,
                            #   sets_provenance.json
                            #   (sources, settings, counts, column meanings, how to read them),
                            #   sets_methods.txt -- references/set-tests.md
    figures/                # volcano / top-protein violins / PCA / heatmap / QC PNGs, one
                            #   sets_<contrast>.png per set-tested comparison (run_sets.R), one
                            #   p-value panel (qc_pvalue_panel.png) + figures.json ("figures":
                            #   captions; "failed": what was not drawn and why) +
                            #   sample_labels.csv (short sample names -> full run names)
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
    Analysis_Report_with_audio.html  # with a podcast: the report with the audio built in (send this one file)
    OUTPUT_FILES.md         # catalog of every file
    comparison/             # (re-analyses) COMPARISON.md + concordance CSVs
  scripts/                  # a copy of the skill scripts that ran this analysis (self-contained)
  logs/                     # commands.log + engine logs
    decisions.md            #   the decisions log (log_decision.py): what was decided, and why
    conversation/           #   CORE-INTERNAL (save_transcript.py): <session-id>.jsonl (redacted
                            #   Claude Code transcripts), index.json, conversation.md
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
- **The Core's feedback survey -- Core runs only.** A run that answers a CoreOmics submission
  (`submission_report.core_run`, the one test every place asks, of the session alone: a record
  that names a submission loads, or session.json announces one by its number or id whose record
  is missing; anything else, unreadable included, is not a Core run. A `--submission` given to
  make_analysis_html feeds only the Submission section) ends
  `Analysis_Report.html`, its `.md` twin and PDF with "**How did we do?** Tell us what you
  thought of this report (and the podcast): a 5-minute survey". It also has that line near
  the end of README and a "**Feedback:**" line in AGENTS.md. The link is
  `core_submission.feedback_url(src, prot)`, the one definition:
  `https://feedback-ucd-proteomics.azurewebsites.net/?prot=PROT_0807&src=report|readme|email|share`
  (share: `Analysis_Report_with_audio.html`, the copy made to be forwarded). `prot`
  is read with the one submission-number pattern (`normalize_submission`) and left out when it
  is not one. The print stylesheet writes the address after the link, so the PDF carries it,
  and `Analysis_Report_with_audio.html` has it because it is the same report. A run outside the
  Core carries it nowhere; when the test cannot be made, the answer is no. A report re-made
  from its own `.md` twin carries the line once (`strip_feedback`). "And the podcast" follows
  `make_podcast.has_podcast`: an invalid `podcast.json` is no podcast. The report is
  usually made before the podcast, so `make_podcast.py link` rewords an existing line to ask
  about the podcast too, in the HTML and the `.md`, before it reprints the PDF. It never adds
  a line to a report that has none.
- **The GitHub star line -- every report.** Just before the survey line (or last, outside the
  Core), `Analysis_Report.html`, its `.md` twin and PDF, and so `Analysis_Report_with_audio.html`,
  end with "**Found this report useful?** Please star the DE-LIMP repository on GitHub -- it
  helps other labs find these tools." Unlike the survey it is on every run: users outside the
  Core are who should find the tools. The words are `core_submission.star_line`, and the address
  is plugin.json's `repository` (no address there, no line). The PDF prints the address after
  the link. `strip_feedback` takes it out of a report re-made from its `.md` twin, so the
  line is there once; `link` leaves it as it is. It is not in README, AGENTS.md or the podcast.
- README.md, README.html and AGENTS.md are each written through a `.part` file and renamed into
  place, so a render that fails leaves the previous file whole (never a 0-byte README.html), and
  finalize zips only the ones it wrote this time.
- **input/raw_files.txt** is written at finalize when missing (a hive_remote session is initialised
  without `--raw`), from `search_provenance.json` `files`, else `output/search/file_list.txt`.
- The registry record is looked up read-only before the zip (`record_run.locate()`); when the
  run-log hook creates it, the three documents are rewritten with its path and go into the zip
  after the hook, like `MANIFEST.txt`.

### The DE tables (`output/tables/DE_<method>_<contrast>.csv`, `run_de.R`)
One row per protein, sorted by `adj.P.Val`. limma's columns (`logFC`, `AveExpr`, `t`,
`P.Value`, `adj.P.Val` — THE significance column —, `B`), the protein annotation
(`Protein.Group`, `Genes`, `Protein.Names`, and on dpc limpa's `NPeptides` / `PropObs`), then
two kinds of detection column, from `Detection_Matrix.csv` (the one definition of
"measured"; a run measured a protein when that matrix is > 0):
- **`Detected_<group>`**, one per group the contrast compares (`Detected_Old_JPH3`,
  `Detected_Old_IgG`): `k/n` = the protein was measured in `k` of that group's `n` runs.
  On dpc the other runs' values were inferred by the detection-probability model; on maxlfq
  they are missing.
- **`Evidence`**: `measured in both` (every run of every group in the contrast measured it) ·
  `presence call` (never measured in at least one group — the difference there is the
  model's, not a measurement: read it as present/absent) · `partly inferred` (dpc) /
  `partly missing` (maxlfq) otherwise. Empty when the protein has no detection record.

The exact wording of both is in `de_provenance.json` → `detection_matrix.de_columns`, which
AGENTS.md reads.

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
`Analysis_Report.pdf` older than the edited HTML is reprinted by `link`
(`html_to_pdf.print_report`); if it cannot be, it is renamed `Analysis_Report.stale.pdf` and
reported as `[SKIPPED]` with the reason. A `*.stale.pdf` never reaches a collaborator: it is
scratch (below), and the next finalize that has a current PDF deletes it (`[OK] … removed the
stale copy superseded by the current PDF`); until then its `[SKIPPED]` line says it is kept on
disk but not shipped.
`output/Analysis_Report_with_audio.html` (`make_podcast.py share`; `session.py finalize` makes
it when a podcast exists) is the same report with the audio and its transcript built in: the
**one file to send** when someone wants the report with the audio, since the card's links are
relative and `Analysis_Report.html` alone has no audio. README ("Start here") and AGENTS.md
name it. `deliver` copies it only while it is current (its `podcast:share` mark matches the
report, `podcast.json` and the audio). One that finalize cannot bring up to date is renamed
`Analysis_Report_with_audio.stale.html` and never sent. Above ~18 MB (~25 MB as an email
attachment, base64), send it through Bioshare.
`podcast/.cache/` (per-chunk TTS audio, ~60 MB for 20 min) is kept on disk for resuming.
`scripts/scratch_files.py` is the one rule for scratch: every `.cache` folder, every `*.part`
file, every `*.stale.pdf` and `*.stale.html`, and `podcast.wav` beside `podcast.m4a`. The session zip leaves them out (`zip_excluded`)
and `make_report.py` does not list them. Never deposit them. The run-registry record
(`record_run.py`) copies what the Listen card points at: the audio `podcast.json` names,
`transcript.html`, `podcast_script.md`, `check.txt` and `podcast.json`. A `podcast.json` that is
unreadable or has a wrong field never stops the report: it is reported as a `[WARN]` and the
report, README and AGENTS.md are made without the podcast.

### The record for a reviewer: `logs/decisions.md` and `logs/conversation/`
Two records let an AI (or a person) check later how the analysis was done; AGENTS.md's
"Reviewing this analysis" section is the checklist (groups and contrasts as confirmed, every
number traced to a command, nothing skipped or re-run unrecorded, no warning ignored).
- `logs/decisions.md` (`log_decision.py`): one entry per decision point, any agent.
- `logs/conversation/` (`save_transcript.py`, Claude Code only): each conversation's transcript,
  found by `CLAUDE_CODE_SESSION_ID` under `$CLAUDE_CONFIG_DIR` or `~/.claude/projects/` and
  **redacted** before it is written (the skill's secret patterns, the secret values this
  computer holds, anything named like a key/token/password/webhook); `index.json` (id,
  first/last timestamp, entries, bytes, sha256, when saved, Claude Code version); and
  `conversation.md`, a best-effort reading (an entry it does not understand is one
  "[unrecognised entry]" line). `session.py init` records the conversation for the plugin
  hook (`hooks/hooks.json`), which saves it before a compaction and at the end of the
  session; it is also saved after the search is submitted, after DE, and by `finalize` (a
  MANIFEST line).
  Images become "[image]". Redaction cannot catch a secret in ordinary prose ("my password is
  …") or one no pattern knows, so both records are **Core-internal**: `deliver` never ships
  them (any `conversation/`, `conversation.md`, `decisions.md` or transcript `.jsonl` is
  `[SKIPPED]`), the session zip leaves `logs/conversation/` and `logs/decisions.md` out, the
  run registry record only points to them, and they are 0640 in a 0750 folder. The collaborator's
  AGENTS.md has no "Reviewing this analysis" section.
- One conversation, several analyses (an orchestrator, staff checking another client): the
  conversation is a run of segments, one analysis each. `session.py init`, `checkpoint.py
  status --resume` and `save_transcript.py --remember` open a segment for the analysis named,
  and so close the one before (an earlier analysis too: X, Y, back to X). Recording the active
  one again, a plain `status` (a peek) and a save change nothing. Each analysis's copy is its own
  segments, joined by "[other work in this conversation omitted]". The first segment starts at
  the user's request that made Claude load the skill (when that request names no other
  submission), else at the load. It starts at `init` instead, with a warning, when the talk
  before `init` names a submission that is not this session's (path or attached record).
  `--transcript-from-now` / `--from-now` start it at `init` on purpose. A copy of the whole
  conversation (one segment, from its first line) keeps every entry. Any other copy keeps an
  ALLOWLIST: the user's and assistant's messages, tool results and compaction boundaries. Every
  other entry type (custom-title, last-prompt, queue-operation, attachment, a system recap, and
  whatever Claude Code adds next) becomes one "[N other entries omitted]" line per run. A
  compaction summary is kept only when the copy holds every entry before it; otherwise it is
  "[compaction summary omitted: it covers conversation outside this analysis's record]".
  `--transcript-from-now` starts the copy at init; compaction summaries and other
  whole-conversation entries are then left out of it, because they cover what came before. The
  user's request that loaded the skill opens the copy only when the skill call answered it
  directly (no tool ran in between).
  **Another client mentioned only by name, not by PROT number, before `init` is not detected
  and lands in this analysis's copy. When the conversation has touched another client, always
  use `--transcript-from-now`.** A peek's output (a plain `checkpoint.py status` on another
  session) and the compaction summaries of a single-analysis conversation also land in the
  active analysis's copy, and any tool output in this analysis's part that lists other clients
  (`ls` of the service tree, `checkpoint.py find`, `locate`) stays in its copy. `index.json`
  and `conversation.md` are rebuilt from the `.jsonl` files present at every save. A copy is
  replaced by a later save when it was made from a different selection (after --from-now or a
  switch), or from the same one and the new save reaches as far into the conversation; an older
  snapshot finishing late is refused (index.json keeps each copy's `selection` and `reach`).
- `MANIFEST.txt`: `[OK]` with a saved conversation; `[SKIPPED]` on Claude Code when nothing
  was saved; `[INFO]` with only the decisions log (another agent, or finalize on HIVE);
  `[SKIPPED]` with neither.

## Re-analysis of the same dataset
Re-running the same raw data (different engine, version, parameters, FASTA, or
design) is common and must not clobber or be confused with the original.

- **Detection:** `session.py find-prior --raw <globs>` scans existing sessions'
  `input/raw_files.txt` for an overlapping raw-file set and reports matches
  (`same_dataset` true when the sets are identical).
- **Fresh by default:** a match is mentioned to the user in one line, nothing more
  (`"default": "fresh"`). The new session reads nothing from the earlier one and does not nest
  under it. It becomes a re-analysis only when the user asks for one.
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
