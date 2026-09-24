# The job-end hook and Slack notifications (`notify_slack.py`)

Every search job this skill writes ends by running one hook, from the compute node, however the
job ends. It does three things, in this order:

1. **Logs the run.** `record_run.py search-done` writes to the Core's run log,
   `/quobyte/proteomics-grp/skill_runs`. It runs on success and on failure.
2. **Hands a finished search to FRAN.** This happens only in the last job of a route that
   ends with a completeness guard, and only on success. Those routes are the DIA-NN ones:
   `report_guard` in the single-shot and 2-job searches, and the quant count in step 5 of the
   chain. The job runs `fran_deposit.py stage --out <out>`, with these added when they apply:
   - `--fasta-meta <session>/input/search.fasta.meta.json`, when that file exists;
   - `--name "<--fran-name>"`;
   - `--qc` / `--not-qc`.

   The argv is built in exactly the shape of `fran_deposit.stage_argv()`, and a test pins it.
   It is built here rather than imported, so the hook never has to import `fran_deposit.py`.
   `stage` makes its own eligibility check.
   Its one-line JSON result goes into the job log. The hook does not stage in these cases, and
   step 7c's check + stage does it instead:
   - **Every other route** (Radiant, FragPipe, Sage, AlphaDIA) logs
     `fran: left_to_agent (no completeness guard on this route)`.
   - **Within 120 s of the job's time limit** (`squeue -h -j $SLURM_JOB_ID -o %e`), the result
     is `near_time_limit`. Otherwise `stage` may take at most the time left minus 30 s.
3. **Posts to the Core's Slack channel.** The post says what steps 1 and 2 did: `Run log: yes`
   and `Staged for FRAN: yes`, or `skipped: …`. When FRAN's own ingest is unhealthy, it says
   `yes, but FRAN's ingest is unhealthy: …`.

The hook runs whether or not anyone is watching. A Core DIA-NN search still reaches FRAN, and
every search reaches the run log, when the user's laptop is closed and no agent session is
alive.

`session.py finalize` does the same for the analysis. After the zip, it runs
`record_run.py analysis-done --session <dir>` and then posts "analysis complete".

**Every step is non-fatal.**
- Nothing in the hook changes the job's exit status or what the report guard decided.
- A missing helper is skipped cleanly.
- A helper that crashes or hangs is reported as `error`. The caps are `record_run` 60 s,
  `stage` 300 s (or less near the limit) and the whole Slack delivery 25 s, webhook lookup
  included, all inside a `timeout 480` around the hook.

## Telling the user, and opting out

Before submitting a Core search, the agent says in one line that its end will be posted to the
Core's Slack channel and handed to FRAN. If the user says not to, the agent passes the flags at
**generation**:

| What | How | Effect |
|---|---|---|
| No Slack post | `run_search.py … --no-notify` (or `SKILL_SLACK=0` when generating); `finalize --no-notify` | The job carries `--no-slack`. The run is still logged and staged. |
| No FRAN hand-over | `run_search.py … --no-fran` (or `FRAN_DEPOSIT=off` when generating) | The job carries `--no-fran`. On success it runs `stage --out <out> [--name …] --skip`, with the same time bounds as any stage, which stages nothing but records `opted_out` in `<out>/fran_deposit.json`, so a later `fran_deposit.py backfill` leaves the search alone. The post says `skipped: off for this search`. Step 7c must then not stage it either. |
| An instrument QC / standard run | `run_search.py … --qc` | The job never stages it. On success it runs `stage --out <out> --name … --qc`, which stages nothing and records the QC decision, so a later stage or backfill honours it. The post says `skipped: instrument QC / standard run`. The decision is made at generation, because a QC session named without a "QC" token would otherwise be staged by the hook before the agent's own `stage --name "... QC"`. |
| FRAN's name for the search | `run_search.py … --fran-name "<the session's descriptive name>"` (always) | Passed to the job's `stage --name`. `--not-qc` (the user said it is not a QC run, whatever its name) is passed through to `stage` as well. |
| The run log | — | It always runs where `record_run.py` is installed, which is what makes it complete. |

**Generation-time flags, not runtime environment.** In `hive_remote` every `hive_exec.sh` call
is a fresh shell, so an environment variable set while generating never reaches the job. The
choice is baked into the hook as a flag. It is also recorded in
`search_provenance.json` → `job_end_hook`, e.g.
`{"run_log": "on", "slack": "on", "fran": "stage" | "left_to_agent" | "off", "fran_name":
"Mouse liver KO vs WT", "qc": true | false | null}`.

`diann_parallel.py` and `radiant_parallel.py` take all of these flags, and `run_search.py`
passes them on. At run time, `SKILL_SLACK=0` and `FRAN_DEPOSIT=off` in the job's own environment still work
as well.

## Setting up the webhook (once, by a Core admin)

1. In Slack, create the channel, then an app with **Incoming Webhooks** turned on, and add a
   webhook for that channel. Copy its URL (`https://hooks.slack.com/services/…`).
2. On HIVE, store it where only `proteomics-grp` can read it. `read -rs` keeps the URL out of
   the shell history and off the screen:
   ```bash
   install -d -m 2750 -g proteomics-grp /quobyte/proteomics-grp/.config
   read -rs HOOK      # paste the URL, press Enter (nothing is echoed)
   ( umask 027; printf '%s\n' "$HOOK" > /quobyte/proteomics-grp/.config/skill_slack_webhook )
   unset HOOK
   chgrp proteomics-grp /quobyte/proteomics-grp/.config/skill_slack_webhook
   chmod 640 /quobyte/proteomics-grp/.config/skill_slack_webhook
   ```
3. Send one test message. From a laptop set up for `hive_remote`, in the skill folder, run
   `bash scripts/hive_exec.sh 'python3 - --test' < scripts/notify_slack.py`. On HIVE, run
   `python3 scripts/notify_slack.py --test`. It prints `sent` or `not sent: …`. Adding
   `--dry-run` names the webhook source it would use, and says why it is off if it is, without
   sending.

**To rotate the URL,** replace the file's contents. Jobs read the file when they post.

## Where the webhook is looked for (first one set wins)

| Order | Source | For |
|---|---|---|
| 1 | `$SKILL_SLACK_WEBHOOK` | a one-off override |
| 2 | `~/.config/ucdavis-proteomics/slack_webhook` (`chmod 600`) | a personal channel, or a laptop |
| 3 | `/quobyte/proteomics-grp/.config/skill_slack_webhook` | the Core channel (on HIVE; group-readable only) |
| — | none of these | Slack is off: `off: Core notification not configured for this user` |

Only a `https://hooks.slack.com/` URL is used. Any other value switches Slack off; the skill
does not fall through to the next source or post somewhere else. Job logs, status lines and
`MANIFEST.txt` never name the files or the group, because a collaborator's copies must not
describe the Core's internals. Only `--dry-run` / `--test` show the lookup detail.

**Laptop finalize, webhook only on HIVE.** Sometimes `finalize` runs on the user's computer,
no webhook resolves there, and a HIVE login is saved (`HIVE_USER`, or
`~/.config/ucdavis-proteomics/hive.env`). Then the notifier sends the message facts to itself
on HIVE through `hive_exec.sh` and posts from there.

## What finalize writes to MANIFEST.txt

The run log and the post each add one line, run log first, as the last lines of `MANIFEST.txt`
and of the zip's copy. Neither is an export part.

| Line | When |
|---|---|
| `[OK]      Core run log -- logged` / `[OK] Slack notification (Core channel) -- sent` | done |
| `[INFO]    Core run log -- not configured for this user`, `[INFO] Core notification -- not configured for this user`, `… -- off (--no-notify)`, `… -- off (SKILL_SLACK=0)`, `… -- off (RECORD_RUN=off)` | not configured, or opted out. A notice: the agent does **not** relay `[INFO]` lines as missing parts. |
| `[SKIPPED] Slack notification (Core channel) -- not sent: …` / `[SKIPPED] Core run log -- record_run.py failed …` | attempted and failed |

**Two fallbacks keep the manifest sound:**
- If the final append of `MANIFEST.txt` into the zip fails, finalize retries once with the
  manifest as it was before these two lines. The result's `zip_manifest` says `added`,
  `added without the run-log/Slack lines (…)`, or `MISSING: …`.
- Every note is flattened to one line.

## Which job does what

| Route | Last job | Earlier jobs |
|---|---|---|
| DIA-NN, library-free (2 jobs) | `job_2_search.sh`: logs, posts; **stages** on success | `job_1_lib.sh`: logs + posts a failure |
| DIA-NN, one job | the job: logs, posts; **stages** on success | — |
| DIA-NN 5-step chain | `step5_report.sbatch`: logs, posts; **stages** on success | steps 1, 1b, 2, 3, 4: log + post a failure |
| Radiant 3-step chain | `step3_fulcrum.sbatch`: logs, posts; FRAN `left_to_agent` | steps 1, 2: log + post a failure |
| Sage / FragPipe / AlphaDIA / single Radiant | the job: logs, posts; FRAN `left_to_agent` | — |

A job that another job waits on (`afterok`) reports its failure because the chain stops there.
An array reports only its first failing task. That task claims the directory
`<out>/.slack_failed_<array job id>` with `mkdir`, and the other tasks see it. `mkdir` is the
primitive measured to be atomic across HIVE nodes on Quobyte: `flock` lost 578 of 800 updates
between 2 nodes, and `mkdir` lost 0. The hook's own marker was tested the same way on
2026-09-24: two tasks on different nodes raced over 200 rounds, and each round had exactly one
winner.

**Nothing is done:**
- For a job stopped by SIGTERM before its time limit (a `scancel`, or a preemption that is
  requeued).
- For an intermediate job that succeeded.
- Outside SLURM.

A time-limit TERM is logged and posted as failed, and never staged.

## What a message contains

| Search done / failed | Analysis complete |
|---|---|
| HIVE user; session name | HIVE user (or local login); session name |
| instrument + acquisition (workflow manifest) | instrument + acquisition |
| number of runs | runs searched / samples in the DE |
| engine + the version that ran, route | engine + version; DE method (`de_provenance.json` label) |
| status (the failing step and its exit code), time since set-up, this job's run time | significant proteins per contrast, with `run_de.R`'s own rule |
| median precursors and proteins per run (DIA-NN `report.stats.tsv`), plus how many runs identified nothing | session folder, the zip, `HOW_TO_SUBMIT.md`; how many export parts were skipped |
| `Run log: …` · `Staged for FRAN: …` | `Run log: …` |
| search folder; the job's log when it failed | — |
| skill problems `report_issue.sh` recorded for that user since the run began (count + file) | the same |

**Never posted:**
- Raw file names, conditions or other sample annotations. Contrast names are the only group
  labels.
- Log excerpts.
- The webhook.

Significant counts are the ones `run_de.R` wrote; they are never recounted.

## Rules the code keeps (from STAN's `stan/notify.py`)

1. **A notifier never breaks its caller.** Every function returns and raises nothing. All of
   `deliver()` runs under one deadline, so a stalled mount under the webhook lookup or the
   issue count cannot hold a job, finalize, or `send_alert`.
2. **The webhook is a bearer credential, and so is anything shaped like a secret.**
   - The webhook is read only inside `notify_slack.py`. It is never passed on a command line,
     written into a job script, printed or logged, and it is scrubbed from every error.
   - **Every string in every payload** (alert bodies and titles, error tails, session names),
     and every log, status and MANIFEST line, first goes through the patterns `report_issue.sh`
     refuses, then through `_scrub`. Those patterns are:
     - private-key blocks;
     - `ghp_` / `github_pat_` / `hf_` tokens;
     - `Authorization: Bearer|Token …`;
     - `password=` / `password:`;
     - `postgres(ql)://user:pass@`;
     - any `hooks.slack.com/services/…` URL.

     A match becomes `[redacted]`, and the message still goes out.
3. **The top-level `text`** is escaped like the blocks (`&`, `<`, `>`), so `<!channel>` or
   `<url|x>` in a session name cannot ping or link.

## Posting from another script (`send_alert`)

`notify_slack.send_alert(text, *, title=None, with_status=False)` posts a one- or two-line
alert, for example from `fran_deposit.py health --alert`.
- **Returns** `True` when Slack accepted the post, `False` otherwise. With `with_status=True`
  it returns `(sent, status line)` instead, and that tuple is always truthy.
- Never raises, prints nothing, and does not relay.
- Bounded at about 25 s, fact-gathering included.
- Text and title are redacted.

Guard the import, because an older install has no `notify_slack.py`.

## How the hook is wired into a job (`wrap_job_script`)

Every generator runs the job script through `notify_slack.wrap_job_script()`.
- **Header.** Every leading comment or blank line (the shebang, `#SBATCH`, a comment between
  them) stays in the header, so no directive can end up below the hook.
- **Body.** Everything after it runs, unchanged, inside a background subshell `( … )&`, which
  the batch shell `wait`s on with an `EXIT` trap and a `TERM` trap set.

The design follows two measurements on HIVE (bash 5.1.16, 2026-09-24):
- **At the time limit, SLURM sends SIGTERM to the batch shell only.** Bash runs a trap only
  after the foreground command returns. With DIA-NN in the foreground, the trap never ran, and
  SIGKILL followed 130 s later (`KillWait`). `wait` is interrupted by a trapped signal, so the
  handler runs at once and then re-raises SIGTERM: the job still dies by the signal (job
  23990802: `TIMEOUT`, batch `0:15`).
- **Without a TERM trap, bash's EXIT trap sees `$? = 0` on SIGTERM.** A timed-out search would
  have been treated as finished, and staged.

The EXIT trap exits with the status `wait` returned, so exit codes, `set -e` and the report
guards are unchanged. The body keeps its own traps.

## Running the skill's tests

Every test that executes a generated job script, or runs finalize, builds its environment
with `tests/job_env.py`. That environment:
- unsets `SLURM_*`;
- sets `SKILL_SLACK=0`, `FRAN_DEPOSIT=off`, `RECORD_RUN=off`;
- points `SKILL_RUNS_DIR` / `FRAN_DROP_DIR` into the test's temp dir;
- sets `HIVE_ENV_FILE=/nonexistent`.

Isolation therefore does not depend on `$SLURM_JOB_ID` being unset: a suite run inside a HIVE
allocation cannot log, stage or post. `tests/test_job_env_guard.py` fails any test module
that runs a generated script without it. `test_slack_notify.py` switches single steps back on,
against a loopback webhook, a fake `record_run.py`, and a temp FRAN drop directory.
