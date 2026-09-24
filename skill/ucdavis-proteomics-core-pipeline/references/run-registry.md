# Run registry: every Core search, recorded on HIVE

`scripts/record_run.py` records every search the UC Davis Proteomics Core runs with this skill
in `/quobyte/proteomics-grp/skill_runs/`. The registry follows the conventions of Brett's
DataAnalysis sessions:

- one **session folder** per search;
- an append-only **master log**, `data_analysis.md`;
- an **activity log**, `activity_log.csv`.

The record is written the moment a search ends, or fails, and updated when the analysis is
finalized. The registry sits beside `skill_issues/`, which is where `report_issue.sh` writes, and
uses the same routing.

It exists because the evidence for one search ends up in many places:

| Evidence | Where it lives |
|---|---|
| `search_provenance.json` | the search folder |
| FASTA sidecar | beside the FASTA |
| parameters | `wf/` |
| SLURM logs | wherever the job wrote them |
| report, session zip | the instrument share, or a laptop |

Answering "what did the Core search last month, for which submission, with which DIA-NN, and how
deep did it go?" meant opening each of these, if anyone still knew where they were. Now the
master log answers it for people, and `record_run.py list` answers it as a table.

## Who is recorded, and who must not be

**Proteomics Core runs only.** The gate is **write permission** on the registry, which lives
inside `/quobyte/proteomics-grp` (mode 2770, group `proteomics-grp`).

- A collaborator's HIVE account is outside the group, so it cannot write there, and their runs
  stay out.
- The filesystem enforces this. The script only reports it, as `"reason": "not_core_member"`.
- It is the same gate as `fran_deposit.py` and `report_issue.sh`.

| From | Route | What happens |
|---|---|---|
| HIVE, registry writable | `direct` | written directly (a SLURM job, a login shell, `hive_exec.sh`) |
| HIVE, registry not writable | `not_core` | **not recorded**: `"reason": "not_core_member"` |
| Off HIVE, HIVE login set up (`HIVE_ENV_FILE` / `hive.env`, or `HIVE_USER` + `HIVE_KEY`) | `ssh` | the script and anything that exists only here (a local session and its zip) are uploaded; HIVE writes the record (and checks the gate again); the upload is deleted |
| Anything else | `none` | not recorded: `"reason": "not_on_hive"` |

**Switches:**
- `RECORD_RUN=off` turns it off (`SKILL_RUNS_DIR=off` is an alias). It returns
  `{"recorded": false, "reason": "disabled"}` before anything is read, written or sent. No SSH is
  attempted.
- `SKILL_RUNS_DIR` is the destination on every route, including on HIVE. Tests point it at a
  temporary folder.
- `record_run.py --where` prints the route from here and does nothing else.

## When it runs

| Moment | Command | Called by |
|---|---|---|
| A search's final job ends, or any job fails (for an array, the first failing task) | `record_run.py search-done --out <search out dir> --status completed\|failed --exit-code N [--session <dir>]` | the job-end hook (`notify_slack.py`); the agent at step 7d |
| `session.py finalize --zip` has run | `record_run.py analysis-done --session <dir> [--out <dir>] [--zip <zip>]` | the finalize hook; the agent at step 12 |

Options on both commands:
- `--prot PROT_0807` (also `0807`, `807`, or the 12-hex CoreOmics id; repeatable).
- `--name` sets the folder name for a new record.
- `--issues-tag` gives the `report_issue.sh` tag when it is not the session name.
- `--dry-run` writes nothing and prints the SEARCH_LOG.md it would write to stderr.

### When the chain dies before its last job, the agent records the failure

The job-end hook runs inside a job, so there are two cases in which no hook fires:
- **SLURM kills the job:** OOM, TIMEOUT, node failure or `scancel`. The job's own trap may never
  run.
- **An `afterok` chain stops:** once a step fails, the jobs queued after it never start
  (`DependencyNeverSatisfied`).

In both cases no record is written unless the agent writes one. **So whenever `watch_run.sh`
reports `failed`, and before resubmitting, record it**. On HIVE, through `hive_exec.sh` in
hive_remote:

```bash
python3 scripts/record_run.py search-done --out <search out dir> --status failed \
    [--exit-code <N from sacct>] [--step <the failed step's job name, e.g. s3_assembly>]
```

- **Safe to repeat.** If the hook did record it, this updates the same folder: the
  `search_failed` activity row says `re-recorded`, and the master log does not repeat the entry.
- **Inferred values.** Without `--step`, the failing step is the first job sacct shows in a
  failed state (FAILED, TIMEOUT, OUT_OF_MEMORY, CANCELLED, NODE_FAIL, ...). Without `--exit-code`, the exit code is the last job's
  sacct ExitCode.
- **After a fix.** When the resubmitted chain completes, its last job's hook records
  `completed`. That adds a dated update line to the master log and a `search_completed` row, in
  the same folder.

**Output contract.** stdout is ONE JSON object:
- `{"recorded": true, "path": "<session folder>", "status": ..., "prot": ..., ...}` when a record
  was written;
- `{"recorded": false, "reason": "<code>", "detail": "..."}` when it was not.

The reason codes are `disabled`, `not_core_member`, `not_on_hive`, `out_not_found`,
`session_not_found`, `bad_input`, `ssh_failed`, `timeout`, `error` and `dry_run`. The exit status
is 0 whatever happened; only an argparse usage error exits non-zero. Diagnostics go to stderr.

**Bounded.** The job-end hook allows 60 s per call, so `--timeout` defaults to 45 s, with a
SIGALRM backstop at 50 s. It never walks a directory: `.quant` and XIC folders are listed and
counted, not traversed. There is a per-file cap and a per-record copy budget.

**Idempotent.** Calling again for the same search updates its folder:
- **The folder:** a later failure, a resubmission that completes, and the finalize all land in
  one folder.
- **`run_record.json`:** keeps the `history` of every call.
- **Status without `--status`:** it is inferred, with sacct supplying the job states and times.
  A failed first attempt followed by a newer report counts as completed. An inference that
  cannot tell, or that agrees with what the job said, keeps the job's own word: its exit code
  and the step that reported.
- **Failing step:** inside a job it comes from `$SLURM_JOB_NAME` (or `--step`). A failure
  recorded from step 3 of the 5-step chain is filed under the same folder the final step would
  use.

## Layout

```
/quobyte/proteomics-grp/skill_runs/
  data_analysis.md            master log (below)
  activity_log.csv            timestamp,session,action,tool,target,status,notes
  sessions/<YYYY-MM-DD_Short-Description>/
    SEARCH_LOG.md             the run, for people
    run_record.json           the same facts, machine-readable (schema_version 2)
    README.md  MANIFEST.txt   the session's own (after finalize)
    <session>.zip             the session zip, without per-run .quant files, when under the cap
    input/                    conditions.csv, <fasta>.meta.json, params.cfg + .rationale.json,
                              workflow.manifest.json, raw_files.txt -- NEVER raw data or the FASTA
    output/                   *.docx (the report of record + Methods), methods.md,
                              AI_Analysis_Report.md, AUDIT.*, SAMPLE_QUALITY.*, tables/,
                              figures/ (files up to 5 MB), DATA_SUBMISSION/HOW_TO_SUBMIT.md
    output/search/            search_provenance.json, report.stats.tsv, report.log.txt, one SLURM
                              log per job (per-file task logs only when the search failed), FRAN
                              receipt, window.json / massacc.txt / params.resolved.cfg (5-step
                              chain), jobs.txt, submit.sh -- and SYMLINKS report.parquet and
                              search_out to the real search
    scripts/                  commands.log, reproduce.sh, REPRODUCE.md, run_manifest.json
  .index/out_<16 hex>         -> ../sessions/<name>   (and session_<16 hex>)
```

- **The folder name** is the skill's session name, which already has the DataAnalysis form
  (`2026-09-24_Silva_LRS_JPH_Kv21_RyR`). With no session it is `<submission date>_<search dir
  name>`; a generic name such as `search_out` climbs to its parent (`2026-09-23_run1`).
- **No hash in the name.** The identity is the search out dir's real path. It is stored in
  `run_record.json` and indexed in `.index/`, so a re-record finds its folder whatever it is
  called. The name is kept once chosen.
- **A name collision.** A different search that wants the same name gets `<name>_2`, `_3`, and
  so on. Creating the folder claims the name, so two writers never share one.
- **Who ran it** is the owner of the search folder, so a maintainer recording someone else's
  search does not re-attribute it.

## SEARCH_LOG.md

The log's sections, in order:

1. **Header:** status, exit code and failing step, the engine and the version that ran, the file
   count, instrument and acquisition, and who ran it.
2. **CoreOmics submission.**
3. **When:** submitted, started and finished, from sacct.
4. **Data Quality Notes** (below).
5. **Errors:** an excerpt from the engine and job logs, when the search failed.
6. **SLURM jobs:** an array collapses to its task counts.
7. **Data:** the files, vendor, instrument, acquisition, and where the raw files live.
8. **Engine:** version and where that came from, the binary, and the route.
9. **Key parameters**, each with where it came from:
   - the precursor m/z range;
   - MS1/MS2 mass accuracy and its source, plus what the engine log says it optimised;
   - the Orbitrap resolution and its source;
   - the scan window;
   - FDR;
   - the enzyme;
   - peptide length and charge;
   - modifications;
   - the library and MBR.

   The values are read by `make_methods.search_record()`. That is the one reader the Methods and
   the SDRF use too, so the three cannot disagree.
10. **Sequence database:** the organism, the proteome and its entries, and the UniProt release.
    Also the contaminants and how they were handled (`--cont-quant-exclude`,
    `fetch_fasta.py --enzyme`).
11. **Results:** the median precursors and proteins per run, and a per-run table from
    `report.stats.tsv`.
12. **Where the outputs are:** report, engine log, `.quant`, XIC, library.
13. **FRAN:** the hand-over receipt.
14. **Skill issues** recorded for the session.
15. **Analysis**, after finalize:
    - the Word documents;
    - the DE method and thresholds, and significant proteins per contrast;
    - the overall audit result;
    - where `methods.md`, `DATA_SUBMISSION/` and `REPRODUCE.md` are;
    - `MANIFEST.txt` counts;
    - what happened to the zip.
16. **Expert Review Notes**, copied from the session README when it has them.
17. **Copies, not-copied items and links.**
18. **Record history.**

### Data Quality Notes: always present

DataAnalysis rule: every session gets a "Data Quality Notes" section, even if it just says
"Nothing anomalous observed". The skill's session README does not have one yet, so for now it
lives here.

Each note uses the five-part form where the source gives enough: *what · where · why it matters
· likely cause · suggested fix*. Each is marked **CRITICAL**, **WARNING** or **NOTE**.

The notes are gathered from what the skill already produced:
- a failed search, with its failing step and the error lines;
- runs with zero IDs;
- `detect_acquisition.py` output (per-file warnings, `orbitrap_resolution_unknown`,
  `resolution_mixed`, `ms2_ion_trap`, low confidence, mixed m/z range), or a rationale that says
  `orbitrap_generic`;
- FASTA sidecar warnings, contaminants dropped as identical to a target, and a **legacy
  database**: a sidecar with contaminants but no `contaminant_target_rule`, built before
  identical-to-target contaminants were removed;
- `AUDIT.json` WARN/FAIL findings;
- `SAMPLE_QUALITY.json` flags (a flag confounded with a group is CRITICAL);
- `session_zip_contains_quant`;
- every `report_issue.sh` file for the session;
- **"CoreOmics submission: not recorded"**, when it is not.

## CoreOmics submission

DataAnalysis rule: every session links to its CoreOmics submission, `PROT_####` plus a 12-hex
id. The script looks for it in this order:

1. `--prot` (`PROT_0807`, `0807`, `807`, `#807`, a 12-hex id, or both);
2. a `coreomics` / `prot` / `internal_id` / `coreomics_id` key in the session's `input/*.json` or
   top-level `*.json`;
3. a `.core_submission.json` receipt in the session or up to three folders above it;
4. the search's `fran_deposit.json` / `search_provenance.json`.

It is **never taken from a folder name**, because those are free-form. When the submission is
found, it is written into SEARCH_LOG.md, `run_record.json`, the master-log entry and the activity
rows. When it is not, the log says `not recorded` and adds a Data Quality Note. Re-recording with
`--prot` fills it in.

## Locking: `mkdir`, not `flock`

Several writers on different HIVE nodes share `data_analysis.md`, `activity_log.csv` and a
run's `run_record.json`. **`flock` does not work across nodes on /quobyte.** Measured on
2026-09-24: two SLURM jobs on different nodes wrote 400 times each to one file. With `flock`,
578 of the 800 updates were lost, about the same as with no lock at all (567), and readers saw
torn JSON.

The registry uses FRAN's lock instead (`ingest/auto_ingest_state.py`), which lost 0 of 800 in
the same test. Re-measured for this registry's own appends on 2026-09-24 (job 23991970, two
nodes, 400 appends each to one `activity_log.csv`-shaped file on /quobyte):

| Lock | Result |
|---|---|
| mkdir lock | 800/800 intact, 0 lost, 0 torn, 0 lock timeouts (about 0.17 s per append under contention) |
| none | one node's writer failed with `EIO` on `close()`, and all 400 of its rows were lost |

The lock works like this:

1. Take the lock by creating `<file>.lock.d` with `mkdir`, retrying with jitter. Write a unique
   owner token into `<file>.lock.d/owner`.
2. A lock older than 60 s belonged to a writer that died. Breaking it is **judge -> rename ->
   verify** (FRAN commit ad81863):
   - note the owner token of the lock judged dead;
   - rename the lock directory away (of several breakers, only one rename succeeds);
   - check that the renamed directory still carries that token. If it does not, the lock changed
     hands in between: A died, C broke A's lock and took a fresh one, and B, which had judged
     A's lock, renamed C's live one. B then puts it back and waits.

   Without the verify step, two writers can hold the lock at once when three or more contend.
3. Open the file only after the lock is held, so the append lands at the true end of the file.
4. Release the lock by removing the directory, **but only if its owner token is still ours**. A
   writer that held the lock past 60 s may have had it broken and taken over; removing the
   directory then would free someone else's lock. Tests cover the A/B/C interleaving and the
   ownership check.

Waiting is bounded at 10 s (`RECORD_RUN_LOCK_WAIT`). Past that, the write happens unlocked and
the JSON result lists the file under `"lock_timeout"`: a late record is better than none.

`run_record.json` and `SEARCH_LOG.md` are written through a temporary file, `fsync` and a rename,
so a reader sees the old file or the new one, never half of one. A reader retries a JSON file
that does not parse, because across nodes even a renamed file can read torn for a moment. A
`run_record.json` that stays unparseable for about 6 s is **moved aside** to
`run_record.json.unreadable-<time>`, never deleted or overwritten. The rebuilt record says so
in a Data Quality Note.

The record's lock is held only while the record is read, merged and written, never during a
long copy:
- **(A)** under the lock, merge this event into the record and write it;
- **(B)** copy, unlocked;
- **(C)** under the lock again, re-read the record, apply this call's copy results to it, append
  to the shared logs, and write.

## The master log: `data_analysis.md`

The file is append-only. Each entry is a single `write()` under the lock above. The header is
written only when the file is created: a temporary file is hard-linked into place, so no writer
can append before the header exists.

- **The first event of a run** is an entry of the form `## <date it happened>: DIA-NN 2.7.0
  search, 3 × .raw (Orbitrap Fusion Lumos, DIA), Homo sapiens -- COMPLETED`. It links the
  session folder and gives the CoreOmics submission, who ran it, the search folder and the
  median depth.
- **Every later event is one dated line**, for example `### 2026-09-25 update -- <name>: search
  completed` after a failure, or `... analysis complete`. The analysis-complete line lists the
  **Word report first**, then the Methods, the significant counts, the zip or reproducibility
  path, `REPRODUCE.md` and the submission.
- **No duplicates.** Each entry ends with a `<!-- record_run <key> <event> <status> -->`
  marker. The marker is checked under the same lock, so the same event is never logged twice.

## The activity log: `activity_log.csv`

The columns are exactly `timestamp,session,action,tool,target,status,notes`. Timestamps are ISO
8601 to the minute, with the offset (`2026-09-24T12:45-07:00`). Each row is quoted by `csv`,
kept on one line, and short. It is written like the master log: one `write()` per row, under
the lock, append-only, with the header written only at creation.

| action | when | logged |
|---|---|---|
| `search_submitted`, `search_started` | from sacct, timestamped when they happened | once per run |
| `search_completed` / `search_failed` | every `search-done` (a repeat says `re-recorded`) | every call |
| `fran_staged` / `fran_skipped` | receipt found / no receipt at finalize | once |
| `issue_recorded` | each `report_issue.sh` entry for the session | once each |
| `analysis_completed` | every `analysis-done` | every call |

## What is never copied

- **Raw data, the FASTA itself, and per-run `.quant` files.** Single-shot searches keep `.quant`
  in `quant/`; the 5-step chain keeps them in `quant_step2/`, `quant_step2_orig/` and
  `quant_step4/`. The record lists where they are and how many there are.
- **Anything over the per-file cap:** 20 MB (`--file-cap-mb`), or 5 MB for figures.
- **Anything beyond the per-record budget:** 200 MB (`RECORD_RUN_COPY_BUDGET_MB`). These files
  are listed under "Not copied" with the reason.
- **Anything that looks like a credential.**
  - By name: `*token*`, `*webhook*`, `*secret*`, `*passw*`, `*credential*`, `.env`, `hive.env`,
    `id_rsa`/`id_ed25519`, `*.pem`, `*.key`, `.netrc`, `.pgpass`.
  - By content: a private-key block, a GitHub or Hugging Face token, an `Authorization:` header,
    a Slack token or webhook URL.

  The whole registry is readable by the Core group, so a refused file is listed with the reason
  and its contents are never quoted.

## The session zip

`analysis-done` copies `<session>.zip` to the folder's top level when it is under the cap: 5 GB
by default (`--zip-cap-gb`), or 1 GB when it has to be uploaded from a laptop
(`RECORD_RUN_SSH_ZIP_CAP_GB`). Over the cap, the record points at the zip and says why.

**Session zips DO contain the search's `.quant` files, predicted library and XICs** whenever the
search ran into `output/search/`. `session.py finalize --zip` leaves out only
`output/raw_data/` and `DATA_SUBMISSION/upload_staging/`. This was verified by running finalize
on a synthetic session.

A real 15-file Lumos chain session holds 2.44 GB under `output/search`:

| Content | Size |
|---|---|
| 45 `.quant` files (step 2, its `_orig` backup, and step 4) | 1.41 GB |
| `step1.predicted.speclib` | 0.69 GB |
| XIC | 0.23 GB |

`.quant` files are DIA-NN intermediates. The analysis is reproduced from `report.parquet` and the
parameters, not from them, so `session.py finalize` is being changed to leave `*.quant` and
`quant*/` out of the zip (it will record the count and path in `zip_excluded`). Whatever a zip
holds, the registry copies it **without** its `.quant` members. It moves the kept
members' compressed bytes as they are, with no recompression, so the copy is disk I/O rather
than CPU on a login node, and it checks the copy before using it. It also records a
`session_zip_contains_quant` finding.

## The directory README

`record_run.py` writes `skill_runs/README.md` itself. It writes the file when it is missing,
or when the version marker at its foot is older than `README_VERSION`, so the directory never
describes a layout it no longer has. The text lives in one place, `README_TEXT` in
`record_run.py`, and is quoted verbatim below. `tests/test_record_run.py` fails when this copy
and the script disagree. Edit the script, raise `README_VERSION`, and paste the result here.

<!-- README:BEGIN -->
```markdown
# Skill run registry

Every search the `ucdavis-proteomics-core-pipeline` skill runs for the UC Davis Proteomics Core,
recorded by `scripts/record_run.py` when the search ends (or fails) and updated when the analysis
is finalized. Laid out like the Core's DataAnalysis sessions.

- `data_analysis.md` -- the master log: one entry per run; later events (a failure fixed, the
  analysis finalized) are appended as dated lines, never edited in.
- `activity_log.csv` -- `timestamp,session,action,tool,target,status,notes`, one row per action
  (search submitted / started / completed / failed, FRAN staged or skipped, skill issue recorded,
  analysis completed).
- `sessions/<YYYY-MM-DD_Short-Description>/` -- one folder per search, named after the skill's
  session (else `<date>_<search folder name>`; a different search wanting the same name gets
  `_2`). Start with `SEARCH_LOG.md`: the CoreOmics submission, Data Quality Notes, status, engine
  and the version that ran, key parameters and where each came from, results, and where every
  output is. Beside it: `run_record.json` (the same, machine-readable), the session `README.md`
  and zip, `input/` (conditions, FASTA sidecar, parameters, `raw_files.txt`), `output/` (the
  Word report of record, Methods, tables, `search/` logs and a link to `report.parquet`) and
  `scripts/` (commands, reproduce script).
- `.index/` -- how recording the same search again finds its folder (by the search folder's real
  path, not its name). Leave it alone.
- `*.lock.d` -- a writer's lock, held for a second or two; one older than 60 s is broken
  automatically. Do not use `flock` on these files: it does not lock across HIVE nodes.

Never copied: raw data, FASTAs, per-run `.quant` files, anything over 20 MB, or anything that
looks like a credential -- the record says where they are.

Only Proteomics Core (proteomics-grp) accounts can write here, so collaborators' runs are never
recorded. To keep a run out: `RECORD_RUN=off`. To list the registry:
`python3 <skill>/scripts/record_run.py list [--since YYYY-MM-DD] [--user U] [--status failed]`.
Skill problems go to the sibling folder `../skill_issues/`.

<!-- record_run.py README v1: written by the skill; edit README_TEXT in record_run.py, not this
file -->
```
<!-- README:END -->

## `list`

```bash
python3 scripts/record_run.py list                                   # a table
python3 scripts/record_run.py list --since 2026-09-01 --user gabrig  # filters (--user/--status repeat)
python3 scripts/record_run.py list --status failed --tsv             # every column, for a spreadsheet
python3 scripts/record_run.py list --json
```

`list` reads `sessions/*/run_record.json`. Off HIVE, with a HIVE login, it lists over SSH.

## Tests

`tests/test_record_run.py` needs no network, HIVE or SLURM:
- `SKILL_RUNS_DIR` always points into a temporary folder;
- `HIVE_EXEC` is a stub over another temporary folder;
- `sacct` is either disabled (`RECORD_RUN_SACCT=none`) or replaced by a stub;
- `SLURM_*` and `HIVE_*` are removed from the environment.

No test can reach the real registry.
