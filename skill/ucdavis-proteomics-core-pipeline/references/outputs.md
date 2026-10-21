# Output packaging, re-analysis & comparison

Every run is packaged into a tidy **session directory** so people can find things.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 20: Where the session goes — ask the user
- 31: Session layout (`session.py`)
- 296: Re-analysis of the same dataset
- 327: Comparing analyses (`compare_analyses.R`)
- 346: 1b. A prior analysis: a new session is fresh (step 1b)
- 358: 1b. A re-analysis, only when the user asks for one (step 1b)
- 367: 3b. Creating the analysis session (step 3b)
- 457: 11. The output-files report (step 11)
- 478: 11b. A re-analysis: compare to the original (step 11b)
- 495: 12. Finalize the session, and the summary to the user (step 12)
- 571: Keeping the record of the analysis: the decisions log and the conversation
<!-- /contents -->

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
                            #   acquisition.json (step 2's detect_acquisition.py output),
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
    podcast/                # audio discussion (make_podcast.py; references/podcast.md), by default for a CoreOmics submission:
                            #   podcast.m4a, podcast_script.md (with its claims ledger),
                            #   transcript.html, podcast.json, check.txt; .cache/ is scratch
    Analysis_Report_with_audio.html  # with a podcast: the report with the audio built in (send this one file)
    OUTPUT_FILES.md         # catalog of every file
    comparison/             # (re-analyses) COMPARISON.md + concordance CSVs
  scripts/                  # a copy of the skill scripts that ran this analysis (self-contained)
  logs/                     # commands.log + engine logs; Core-internal records (never delivered,
                            #   not in the session zip): decisions.md, conversation/, and the
                            #   staff-only QC check (qc_bracket.py, step 8e): qc_bracket.json
                            #   (schema_version 1), qc_bracket.md, qc_bracket_ack.json
    decisions.md            #   the decisions log (log_decision.py): what was decided, and why
    conversation/           #   CORE-INTERNAL (save_transcript.py): <session-id>.jsonl (redacted
                            #   Claude Code transcripts), index.json, conversation.md
```

- `session.py init --name "..." --raw <globs> [--base <path>] [--reanalysis-of <prior> [--beside]]`
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
  `presence call` (never measured in at least one group but measured in another — the
  difference there is the model's, not a measurement: read it as present/absent) ·
  `not measured in either group` (measured only in other groups of the study: **not tested
  in this contrast** — see below) · `partly inferred` (dpc) / `partly missing` (maxlfq)
  otherwise. Empty when the protein has no detection record.

The exact wording of both is in `de_provenance.json` → `detection_matrix.de_columns`, which
AGENTS.md reads.

**Rows without a p-value.** A `not measured in either group` row has `logFC`, `t`,
`P.Value`, `adj.P.Val` and `B` empty (NA): it was not tested in that contrast and is not in
its Benjamini–Hochberg adjustment, which runs over the tested proteins only (2.12.0;
`references/de-analysis.md`). On maxlfq a protein with no value at all in one compared group
has no p-value either (limma cannot estimate it). Proteins *tested* in a contrast are the rows
with a `P.Value`; `de_provenance.json` → `de_tables` gives each table's `n_proteins`,
`n_tested`, `n_not_measured` and `n_significant`, and `not_measured` the rule.

### The podcast (`output/podcast/`, `make_podcast.py`)
An AI-generated audio discussion of the finished report: made by default for a CoreOmics
submission, and only on a yes for internal Core work and anyone else (`references/podcast.md`). It is a **derivative, not a record**: README and AGENTS.md say so,
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
  live with their original. With `--beside` as well, it goes **next to** the prior, in the
  prior's parent folder, as `<date>_<name>_v<N>`: N is one above the prior's version and every
  `_v<M>` of that name already there, so an existing folder is never reused, and nothing is
  written into the prior (a delivered, read-only-by-policy folder; a Core re-analysis, 2026-09-25).
- **The link:** `.reanalysis_of` (the prior's path) and `reanalysis.json` (`reanalysis_of`,
  `version`, `prior_version`, `placement`, and the prior's `search_provenance.json` /
  `de_provenance.json` / README / raw list as they were). `session.reanalysis_of()` is the one
  reader: finalize's `DIFFERENCES.md`, the README and AGENTS.md ("Re-analysis of … (this is
  version N)") and the run log's entry use it.
- **Finalize again** (to add a podcast, say): an unchanged analysis adds nothing to the Core run
  log and is not posted to Slack a second time (`references/run-registry.md`).
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

## 1b. A prior analysis: a new session is fresh (step 1b)

```
python3 scripts/session.py find-prior --raw /path/to/*.d
```

**A new session is a fresh analysis by default**, match or no match. When `matches` is not
empty, tell the user in one line that these files were analysed before, and where. Then carry
on fresh: read nothing from that session (conditions, FASTA, parameters, search output,
results) and do not nest under it. A Core staff member starting a fresh run was steered into
the old results by the earlier wording here.

## 1b. A re-analysis, only when the user asks for one (step 1b)

It is a **re-analysis** only when the user asks for one ("re-analyse", "redo that analysis
with…", "compare with the earlier run"). Then pass that session to `session.py init
--reanalysis-of` (step 3b), so the run nests under `<prior>/reanalysis/` and gets a
`DIFFERENCES.md`, and run the Comparator at the end (step 11b). When the earlier session is a
delivered folder that must not be written to, add `--beside`: the run goes **beside** it as
`<date>_<name>_v<N>` (never an existing folder), linked to it by `reanalysis.json`.

## 3b. Creating the analysis session (step 3b)

**Ask the user where the results should go** — two choices:
- **In the folder with their raw data** (default, recommended — keeps results next
  to the data): pass `--raw` and omit `--base`.
- **A central location** (e.g. their Documents folder, or a path they give): pass
  `--base <that path>` (results land under `<path>/sessions/`).

Then scaffold the session and route **everything** into it:

```
# default — results live with the raw data:
python3 scripts/session.py init --name "<short study name>" --raw /path/to/*.d \
    [--reanalysis-of <prior session dir -- only when the user asked to re-analyse it, step 1b>
     [--beside]]
# or central, if the user chose one:
python3 scripts/session.py init --name "<short study name>" --raw /path/to/*.d \
    --base ~/Documents/DataAnalysis
```

**Core staff:** add step 0a's answer to either form: `--submission PROT_0807` (once per number,
the main one first) or `--internal '<its title on the board>'`. Its JSON's `transcript` says
whether this conversation was recorded for the session; `recorded: false` with a `reason` other
than "not a Claude Code conversation here": tell the user, and do what the reason says.
**`--call-id` whenever any subagent or teammate of this conversation ran in the last 6 hours**
(a team, or step 8d's reviewers -- so `finalize` after them too): add `--call-id <token>` to
`init`, `finalize`, `checkpoint.py status` and `save_transcript.py`. It is how they tell whose
conversation is calling. Type the token literally in the command (never `$(uuidgen)` or a
variable from an earlier call), fresh for every call and unlikely to collide:
`<your agent name, or main>-<6 random letters/digits>`, e.g. `main-k7f3x9`. Run the call in the
foreground: a background call is never counted.

**Core data with a CoreOmics submission (step 1)? Attach it now -- BEFORE the experiment type below**, where the session lives (`experiment_type.py set --source submission` reads the attached record and refuses without one; set before attach, the submission's ASK never fired):

```
python3 scripts/submission_report.py attach --session <session> --record ~/core/PROT_0807   # fetch's folder
python3 scripts/submission_report.py attach --session <session> --given '{"internal_id": "PROT_0807", ...}'  # no token
```

It stores ONE allowlisted record (`input/submission.json` + `session.json`; never an email,
phone or billing field). The report's Submission section, the Methods' Sample preparation,
the analysis brief and the run log all read it. It is also what makes this a Core run, so the
report, README, AGENTS.md, delivery README and email draft end with the Core's "How did we do?"
feedback survey (`core_submission.feedback_url`); a run outside the Core never shows it.
The record's `notes` are Data Quality Notes for the report: organism vs the FASTA, blank
UniProt, sheet IDs vs raw files, conditions vs the design analysed, and pairing in the sheet (e.g. every mouse under all five IPs — samples from
one mouse are not independent). Re-run `submission_report.py notes --session <S>` once the
FASTA and `conditions.csv` exist.

**Record the experiment type now** (step 3's confirmed answer):

```
python3 scripts/experiment_type.py set --session <session> --type <whole_proteome|ip|proximity|secretome|fractions_depth|fractions_separation|other> \
    --source <submission|user> [--stated "<their words>"] [--bait <gene>] [--controls <group,...>]
```

**Ask the bait once.** With sequences added by `fetch_fasta.py --add-fasta` (EGFP, a TurboID
fusion), ask it **after step 6** (the database is built): `experiment_type.py propose`'s
`ask_bait`, recorded with `experiment_type.py set-bait --bait <accession>` (or `--bait none`), kept
as `bait_confirmed` (who, when). Otherwise ask it here, with `--bait`. `propose` asks nothing
(`ask_bait` null, `bait_recorded` says why) when set-bait was answered or the `--bait` recorded
here names one of the added sequences. Never leave the one added sequence as an unconfirmed bait.
**Record the sample type now** (step 3's answer, in the user's words — or that they did not
say), so the brief and report never fill it in from file names:

```
python3 scripts/sample_type.py set --session <session> --stated "<what the user said the samples are>"
python3 scripts/sample_type.py set --session <session> --not-stated     # asked; the user does not know
python3 scripts/sample_type.py set --session <session> --not-asked      # Core staff, any analysis: not asked (step 3)
```

A later answer (a flagged panel in 8c) replaces it with `--stated`.
**Save step 2's JSON as `<session>/input/acquisition.json`** (the whole `detect_acquisition.py`
output). Step 8e reads each file's instrument and acquisition time from it instead of reading
every `.raw` again, and the run record takes its detection facts from it.

`init` also records this conversation for the transcript hook, and from here on **log each
decision** with `log_decision.py` (see "Keeping the record of the analysis" below).
Core staff, internal work: board.md §3 row "session made" (a submission's is at step 1c.7).

The output's `placement` tells you which was used. **Use the printed `paths` map for
every later step** — put
`conditions.csv` and the FASTA in `paths.input_dir`, search output in
`paths.search_out`, DE results in `paths.de_dir` (= `output/tables`), the
reproducibility bundle in `paths.repro_dir`, the report in
`paths.analysis_report`, and append commands to `paths.commands_log`. Raw files
are recorded (in `input/raw_files.txt`), never copied. With `--reanalysis-of`, the
session nests under `<prior>/reanalysis/<date>_<name>/`; with `--beside` too, it sits next to the
prior as `<date>_<name>_v<N>` and nothing is written into the prior.

## 11. The output-files report (step 11)

`session.py finalize` (step 12) builds this catalog itself when step 11 was not run, from the
same folders, and says so in MANIFEST.txt (`[OK] OUTPUT_FILES.md -- written at finalize`); one
that exists is never rewritten, only appended to. Run step 11 by hand when the user wants the
catalog before finalize:

Catalog everything the run produced so the user knows what each file is:

```
python3 scripts/make_report.py --out OUTPUT_FILES.md \
  --search-out ./search_out --de-dir ./de_results --repro ./reproducibility \
  --extra ./conditions.csv ./search.fasta ./wf ./figures ./AUDIT.md ./SAMPLE_QUALITY.md \
          ./AI_Analysis_Report.md ./Analysis_Report.html ./Analysis_Report.md \
          ./Analysis_Report.pdf ./methods.md
```

`OUTPUT_FILES.md` lists every file (figures, audit, search/DE outputs, the bundle)
with its size and a plain-language description, grouped by purpose, and flags
anything unrecognized (never silently omitted).

## 11b. A re-analysis: compare to the original (step 11b)

When the user asked to re-analyse an earlier session (step 1b), compare the two with the
Comparator (a faithful port of DE-LIMP's Run Comparator):

```
Rscript scripts/compare_analyses.R --out <session>/output/comparison \
  --adjp 0.05 --logfc 0 \
  --analysis "Original:<prior>/output/tables" \
  --analysis "Reanalysis:<this session>/output/tables"
```

It writes `COMPARISON.md` + CSVs: protein-universe overlap, the 3×3 Up/Down/NS
concordance per shared contrast, direction concordance on co-significant proteins,
and logFC correlation. (`session.py finalize` already wrote `DIFFERENCES.md` —
what *settings* changed; the Comparator shows how the *results* changed.)

## 12. Finalize the session, and the summary to the user (step 12)

First check the notes from the Core again (step 0: `notes.py check --json --session <dir>`).

```
python3 scripts/session.py finalize --dir <session> --zip
```

(On HIVE's login node it moves itself to a compute node: run it in the background.) This tidies `output/` (tables→`tables/`, figures→`figures/`) and makes sure the session carries
**publication Methods**: `output/methods.md` + `.docx`, covering LC-MS acquisition, the search
(engine, pinned version, parameters, FDR), the sequence database and contaminants, the DE, and the
instrument-grant acknowledgment. If step 9d did not write them, `make_methods.py` does it here. An
existing `methods.md` is never overwritten; if it lacks a section, a complete draft is written
beside it as `methods_complete_draft.md`. Finalize then writes the **repository-deposit package**
`output/DATA_SUBMISSION/`: `HOW_TO_SUBMIT.md`/`.html`, a pre-filled SDRF `sdrf.tsv`,
`protocols.txt`, `files_to_upload.tsv`, and `prepare_upload.sbatch`. Never run
`prepare_upload.sbatch` yourself: it reads every raw file, so the user submits it with `sbatch`
when ready. Finalize also writes, at the session root, **`README.html`** (the page collaborators
open: summary, links, and where the folder, raw data, search output and FASTA are on HIVE with
their Windows/Mac paths), the same text as `README.md`, **`AGENTS.md`** (a guide for an AI agent
given the folder), `input/raw_files.txt` if it was missing, `DIFFERENCES.md` for a re-analysis,
and **`MANIFEST.txt`**, which lists every part as `[OK]` or `[SKIPPED] <name> -- <reason>`
(→ this file). Then it zips the session, leaving out the raw data, `upload_staging/`, DIA-NN's
`.quant` intermediates (~30 MB per run) and the predicted spectral library (`*.predicted.speclib`,
~0.7 GB; rebuilt from the FASTA + params). Both stay on disk; `zip_excluded` counts them.
`--no-deposit` skips only the package. → detail: `references/deposit.md`.
Finally it logs the run in the Core's run log (`record_run.py`) and posts **"analysis
complete"** to the Core's Slack channel. The post carries the session, instrument, engine +
version, significant proteins per contrast, and where the folder, the zip and
`HOW_TO_SUBMIT.md` are. Neither step can fail finalize, and `--no-notify` skips only the post.
Their outcomes are the last two lines of `MANIFEST.txt`. → `references/notifications.md`.

**Read `MANIFEST.txt` and relay every `[SKIPPED]` line with its reason** (e.g. "the Word copy of
the Methods was skipped: pandoc/python-docx not installed"). `[INFO]` lines are notices, not
missing parts, so do not relay them. Examples are the run log and Slack lines such as "Core
notification -- not configured for this user". Then summarize: data type (instrument
+ acquisition), engine + **pinned version**, mass accuracy **and its source**, the skill
version (`reproducibility/environment/skill.txt`), FASTA source, DE method, and per-contrast significant counts. Point them at the
**session folder** and its **`README.html`** — tell collaborators to open that (double-click);
`AGENTS.md` is for an AI assistant they hand the folder to. Then:

- **`output/Analysis_Report.html`** (the report of record: QC, figures and the interpretation in
  one page; README.html links it first), with `Analysis_Report.pdf` (when a browser printed it)
  and `Analysis_Report.md` (the plain-text twin for NotebookLM)
- `output/OUTPUT_FILES.md` (what every file is)
- **`output/methods.md`** (the publication Methods; its `[... — confirm]` tags must be resolved
  before publishing)
- `output/tables/methods.txt` verbatim (the DE record — don't paraphrase)
- **`output/tables/reproducibility_log.R`** (the analysis as plain R — say this is where the code
  is; it is what most people mean when they ask)
- `output/reproducibility/REPRODUCE.md` (the pinned recipe, for re-running the search too)
- **`output/DATA_SUBMISSION/HOW_TO_SUBMIT.md`** (how to deposit the data in PRIDE or MassIVE from
  HIVE). Say how many `TO-FILL` cells `sdrf.tsv` has (`deposit.sdrf_to_fill` in finalize's
  output), and that these are facts only the user knows — tissue, disease, sex, age, what the
  groups are — which the skill never guesses.
- for a re-analysis, `DIFFERENCES.md` + `comparison/COMPARISON.md`.
- if one was made (9e), `output/podcast/podcast.m4a` and its `transcript.html`: say it is
  AI-generated, that you could not listen to it, and that someone should before it is shared.
  To send the report with the audio, they send `output/Analysis_Report_with_audio.html` alone
  (say its size; above ~18 MB, ~25 MB attached, Bioshare rather than email).
  If none was made: a submission's was skipped, failed or left out (say which, in one line);
  internal Core work, and anyone else whose results go to a collaborator or the lab, are
  offered it once (step 9e).

If anything was recorded with `report_issue.sh` this session, say so in one line and where
it went (`report_issue.sh --where`), so the user knows the Core will see it.

For a Core HIVE run, close with the **FRAN handover** (step 7c) in one line — that the search
was staged for FRAN's ingest cron, or, if it was not eligible, the reason in plain words
("this is a collaborator account, so it stays out of the Core corpus"). Say "handed over,
FRAN ingests on its next pass" rather than implying it is already in the corpus.
If `stage` returned a `health_warning`, add it to that line: the search is handed over, but
FRAN's ingest is stuck (or its code is stale) on FRAN's side.

Core staff, internal work: board.md §3 row "finished".

## Keeping the record of the analysis: the decisions log and the conversation

Two records let an AI (or a person) review later how the analysis was done — another
reproducibility check. Both stay in the session's `logs/`, and AGENTS.md's "Reviewing this
analysis" checklist reads them.

- **The decisions log** — any agent, every analysis. At each decision point, one line:

```
python3 scripts/log_decision.py --session <S> --step '3. design' --what-file - --why-file - <<'EOF'
<what was decided>
---
<the user's words, or the evidence>
EOF
# no usable local Python: the same on HIVE (the text goes over ssh on standard input)
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/log_decision.py --session <S> --step "3. design" --what-file - --why-file -' <<'EOF'
<what was decided>
---
<the user's words, or the evidence>
EOF
```

  The text goes in on **standard input** — the what, a line `---`, the why — in a heredoc
  quoted `<<'EOF'` (`EOF` alone at the start of its line), where nothing is expanded: the
  user's own words are full of quotes, and in `--what "…"` bash would run anything in
  backticks or `$(…)`, here and again on HIVE.
  Log: submissions analysed together and which is the main one (step 1c), the conditions
  confirmed (step 3), the defaults confirmed or changed (step 4), the
  contrasts chosen (step 8), every override (`--unchecked`, `--force`, a gate accepted), and
  anything the user declined or asked for beyond the defaults. → `logs/decisions.md`. It is the
  quickest thing for a reviewer to read, so write it on Claude Code too.
- **The conversation** — Claude Code only. `save_transcript.py` copies this conversation's
  transcript into `logs/conversation/`, **redacted** (the skill's secret patterns, the tokens and
  keys this computer holds, and anything named like a key, token, password or webhook), with
  `index.json` and a readable `conversation.md`. Redaction cannot catch a secret written as
  ordinary prose ("my password is …"), so the conversation and the decisions log are
  **Core-internal**: `deliver` never ships them, the session zip leaves them out, and they are
  readable by the owner and group only. Run it after the search is submitted (step 7), after
  DE (step 8); `session.py finalize` runs it itself and records it in `MANIFEST.txt`:
  ```
  python3 scripts/save_transcript.py <S>                    # the session is on this computer
  python3 scripts/save_transcript.py --hive "$S"            # hive_remote: redacted here, then put on HIVE
  ```

  In hive_remote mode run `save_transcript.py --hive "$S" --remember [--call-id <token>]` right after
  `session.py init` (on HIVE), and the `--hive` save before `finalize`: `init` on HIVE cannot
  see this conversation. **No usable local Python** (Windows staff): skip it. It must redact
  the conversation here, before anything leaves this computer, so it has no HIVE route; the
  decisions log is then the record, as `finalize`'s MANIFEST.txt says. Between those calls, the
  plugin's hook (`hooks/hooks.json` → `transcript_hook.py`) saves it before a compaction and
  when the session ends — a no-op in any conversation not recorded for an analysis.
  No transcript (another agent), or none found: it says so and exits 0. It never fails a step.
- **Run by a subagent or teammate** (an Agent call): it shares the parent's session id, and
  nothing in its environment names it. Whenever any subagent of the conversation ran in the last
  6 hours, every call takes `--call-id <token>` (typed into the command, fresh each time, run in
  the foreground; step 3b): the transcript whose RUNNING call carries it is saved -- a subagent's
  own (`agent-<id>.jsonl`), never the parent's, which holds other work. Without a token, or with
  one found in no running call or in two, it records nothing and says so (`init` and
  `checkpoint.py status` put that in their JSON's `transcript`): re-run with a fresh
  `--call-id`, or name the transcript, as its message says (a subagent's:
  `--transcript <…/subagents/agent-<id>.jsonl> --session-id agent-<id>`).
  The hook saves only the main conversation, so a subagent saves its own before it reports back.
- **One conversation, several analyses** (an orchestrator, or staff looking at another client's
  run): each analysis keeps only its own part.
  - The first part starts at the user's request that loaded the skill (only when the skill call
    answered it directly), or at the load itself. It starts at `init` instead when the talk
    before it names another submission number.
  - `session.py init`, `checkpoint.py status --resume` and `save_transcript.py <S> --remember`
    say "the conversation is now working on <S>". A plain `status` is a peek and moves nothing.
  - A copy that is not the whole conversation keeps only the messages, tool results and
    compaction boundaries; every other entry type is one "[N other entries omitted]" line. A
    compaction summary is kept only when the copy holds everything before it.
    `--transcript-from-now` starts the copy at init; compaction summaries and other
    whole-conversation entries are then left out of it, because they cover what came before.
  - **Another client mentioned only by name, not by PROT number, before `init` is not detected
    and lands in this analysis's copy. When the conversation has touched another client, always
    use `--transcript-from-now`** (`session.py init … --transcript-from-now`, or
    `save_transcript.py <S> --remember --from-now`).
  - These also land in the active analysis's copy: a peek's output (a plain `checkpoint.py
    status` on another session), the compaction summaries of a single-analysis conversation,
    and any tool output in this analysis's part that lists other clients (`ls` of the service
    tree, `checkpoint.py find`, `locate`) stays in its copy.
