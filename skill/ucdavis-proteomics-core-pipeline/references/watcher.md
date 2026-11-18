# Run watcher — monitor searches and auto-correct errors

A search is long (often hours) and runs detached (a background process locally, a
SLURM job on HIVE). **Always watch every search — this is mandatory, not optional.**
The orchestrator must never start a search and walk away, and must **auto-correct
stalls and failures on its own, without waiting for the user** (the user has standing
authorization for this — see golden rule #1's compute confirmation is for *starting*
the search, not for recovering a run that's already approved and in flight). Surface
what you did afterwards; don't ask permission to unstick a job.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 33: Watch the WHOLE chain, not its last job
- 55: Two failure modes it catches
- 67: Loop pattern
- 84: What a poll reads, and says
- 90: Error classes → fixes (what `watch_run.sh` detects)
- 122: Warnings that are not failures: proteins the contaminant list took
- 132: Warnings that are not failures: a DIA-NN --dda search on its first pass
- 143: Warnings that are not failures: Sage LFQ
- 160: Warnings that are not failures: protein inference in step 5
- 175: SLURM specifics (HIVE)
- 194: Node faults: the node check every chain job runs first
- 209: Keeping a search off nodes known to be bad (`--exclude`)
- 234: Auto-recovery playbook (apply autonomously; log every action to `commands.log`)
- 262: Don't
- 280: Progress reporting (`--out`, `--poll`)
- 315: 0c. Is this a resume? (step 0c)
- 351: A closed laptop: `RECOVERY.md` (step 0c)
- 360: 7b. The watch loop as SKILL.md runs it (step 7b)
<!-- /contents -->

## Watch the WHOLE chain, not its last job

For any multi-step chain (`diann_parallel.py`, `radiant_parallel.py`), watch with:

```bash
bash scripts/watch_run.sh --all <search out dir> --hive
```

It reads `jobs.txt` (written by `submit.sh`) and reports every job's state at once. It takes the
first all-digit word of each line as the job id, so the Skyline check's `<label><TAB><id>` lines
(`skyline_work/jobs.txt`) work too, and it reads a failed job's log from the folder or its `logs/`
(where the Skyline jobs write `sky_<step>_<id>.log`).
**Do not** watch only the final job. If an earlier step dies, the final job sits
`PENDING` on a dependency that can never be satisfied — which is indistinguishable
from "still queued" if that is the only job you are looking at. This actually
happened: a step-4 array hit `TIMEOUT`, step 5 went `CANCELLED`, and a watcher
pointed at step 5 alone reported "pending" for hours. `--all` returns
`failed:true` with `first_failed_job` in that situation.

`watch_run.sh` does one poll + diagnosis (state, **stall**, error class, fix); loop it
until the run finishes, and on `failed` **or** `stalled` apply the fix and resubmit.

## Two failure modes it catches

1. **Hard failure** — `sacct` reports a state that ended the job without finishing it
   (FAILED / TIMEOUT / OUT_OF_MEMORY / NODE_FAIL / CANCELLED / PREEMPTED / DEADLINE /
   BOOT_FAIL / REVOKED / SPECIAL_EXIT: `scripts/slurm_states.py`, the one list chain mode and
   single-job mode both read). Matched against known error signatures → a fix.
2. **Stall** — the job stays **RUNNING** but its log stops advancing (a hung file, e.g. a
   pathological DIA-NN run whose RT-window/search phase never returns). `sacct`/`squeue`
   will happily show RUNNING forever. `watch_run.sh --log <log> [--stall-min N]` flags
   this when the log's mtime is older than N minutes (default 15) while RUNNING, emitting
   `stalled: true`. **This is the NA41-class failure the first nail cohort hit.**

## Loop pattern

```
# start the search (background locally, or sbatch on HIVE), capture the log path / job id
n=0
while true; do
  n=$((n+1))      # --poll: one more each poll; --out: real progress
  status=$(bash scripts/watch_run.sh --slurm <jobid> --log <log> --out <out> --poll $n --hive)   # or --log <log> --out <out> --poll $n locally
  echo "$status" | grep -qE '"done": ?true' && break      # no python needed (Windows)
  sleep 300       # SLURM: 5 min or more (each poll is several logins + sacct on the login
                  # node); a local run: 60. For long jobs, schedule a wake-up instead
done
# if failed: read error_class + fix, apply it, resubmit, watch again.
```
The agent itself is the watcher — `watch_run.sh` is its eyes. Surface progress and
any auto-fix you applied to the user.

## What a poll reads, and says

`watch_run.sh` reads the **Slurm state** (`sacct`; `squeue` only for PENDING/RUNNING), never the log alone. `--out <search folder>` is enough: with no `--slurm`/`--all` it takes the job ids from `jobs.txt` (the whole chain) or from the `<step>_<jobid>.log` names, so a FAILED job can no longer pass for a running one. A job that ended is reported as it ended (`FAILED`, `OUT_OF_MEMORY`, `TIMEOUT`, `CANCELLED`, `NODE_FAIL`) with `exit_code`, the log's tail (`log_tail`, from `--out` when there is no `--log`) and an `error_class`; its `progress.summary` never says "running". A normalisation check that stopped to ask is **not** one of these: its job exits 0 (`normalization_check.py run --job-exit0`) and `<S>/output/norm_check/NEEDS_DECISION.json` makes the poll say `needs_decision: true` (`--session <S>`, or a `--log` in `<S>/logs`) -- a question for the user, not a failure, never resubmitted. Only for the job that wrote the marker (its `job_id`): any other finished job of that session keeps its own summary and carries a warning that a decision is owed; a marker that is there but cannot be read is a warning too, never "no decision owed".

`search_provenance.json` `queue` is the queue that was **requested**; `queue_actual` (written by the job-end hook, `queue_actual.py`, from `sacct`) is where the jobs **ran** -- its `moved` is true when a pending job was moved (e.g. to `publicgrp`/`low`). Report the second in the Methods and to the user.

## Error classes → fixes (what `watch_run.sh` detects)
| error_class | signal | fix |
|---|---|---|
| `node_fault` | a FAILED job whose node could not reach the storage: exit **75** (the job's node check: `NODE_FAULT: node=… check=… path=… detail=…`), SLURM `NODE_FAIL`, a storage I/O error in its log (`Transport endpoint is not connected`, `Stale file handle`, `Input/output error`), or a black-hole node: `FAILED`, ExitCode `0:53` (signal 53), Elapsed seconds and **no log at all** (a signal-53 job that wrote a log ran, and is not one; when the log cannot be read — no `--log` and no `--out`, a wrong path, an unreadable file, a failed `--hive` read — the watcher cannot tell, does not call it one, and its `fix` says the log could not be read) — never over OOM/TIMEOUT/CANCELLED. Checked first; the JSON carries `node_fault` (node, evidence, `say`) | **a node problem, not DIA-NN and not the data — say so to the user.** Run `node_fault.py retry --out <search out dir> --job <id>` on HIVE (from a laptop: `bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/node_fault.py retry --out … --job …'`): it resubmits the failed step — an array's failed tasks only, once the array has finished — and every step waiting on it with `--exclude=<node>` (every node that failed this search stays excluded), scancels the dependants left `DependencyNeverSatisfied`, rewrites `jobs.txt` and the session's `.recovery.json`, and records it in `<out>/node_faults.json` and `search_provenance.json` `node_faults`. Each failed task is judged on its own log and sacct row. At most **2** retries per step; it refuses (exit 3) a third, failures spread over more than 3 nodes (the storage, not a node), a step with ANY failed task that is not a node fault — an `OUT_OF_MEMORY` or `TIMEOUT` task, or one with no node-fault evidence: nothing is resubmitted, and `not_node_faults` lists each with its fix (`--force` overrides) — and a chain it cannot map; then tell the user and check the HIVE status page. `complete` is true only when every failed task and every step after it was resubmitted |
| `sage_stack_overflow` | a Sage job that died with `fatal runtime error: stack overflow` (exit **134**), or whose `<out>/sage_failure.json` says `sage_stack_overflow` | an **intermittent** Rust crash, not memory and not the data: **re-run the same job once**, unchanged (jobs from 2.15 export `RUST_MIN_STACK`; regenerate an older one with `run_search.py --sbatch`). Report it if it crashes again with the stack raised |
| `out_of_memory` | oom-kill, `std::bad_alloc`, OUT_OF_MEMORY | more memory, **re-generated, never by editing the sbatch by hand** (the job would run at the new size while `search_provenance.json` still records the old one): a single-job search `run_search.py --mem GB` (e.g. 64→128; recorded as `job_mem`); the 5-step chain `--mem-per-file` / `--assembly-mem`. DIA-NN: fewer threads. Resubmit the job it wrote |
| `preempted` | state PREEMPTED (not requeued: a requeued job reads PENDING/RUNNING/REQUEUED) | not the search, not the data: a higher-priority job took the node on the preemptible queue. Resubmit the step on the high queue (not `--preemptible`), and the chain steps waiting on it |
| `slurm_ended` | state DEADLINE / REVOKED / SPECIAL_EXIT (BOOT_FAIL is a `node_fault`) | Slurm ended it before it finished; check `sacct -j <id> -X`, then resubmit the step |
| `timeout` | TIMEOUT / "DUE TO TIME LIMIT" | more time, re-generated the same way, never by hand: a single-job search `run_search.py --time HOURS` (recorded as `job_time`); the chain `--time-per-file`; or split the run. Resubmit the job it wrote |
| `cancelled` | state CANCELLED, not for time or memory (those are `timeout` / `out_of_memory`): someone ran `scancel`, or Slurm cancelled it because a job it waited on failed. Checked before the log signatures below, which a cancelled job's log must not set (a cancelled FragPipe job once read `empty_results` from its own counts) | not a fault the search's log shows. Find out which (`sacct -j <id> -X -o JobID,State,Reason`; "CANCELLED by <uid>" names who) before resubmitting; in a chain, fix the step that failed first |
| `diann_no_dotnet` | `dotnet: not found` | wrong DIA-NN container (no .NET → `.raw` silently skipped). Use the HIVE **native** build `build_<v>/diann-<v>/diann-linux` or a .NET image |
| `diann_zero_ids` | `Number of IDs at 0.01 FDR: 0` | Completed but identified nothing — a **silent null result**, not a crash. **Preserve `report.log.txt` before re-running**; it is the only evidence. Check: library generated precursors → mzML really are DIA with isolation windows (`detect_acquisition.py`) → FASTA matches organism → mass accuracy. DIA-NN's warning about generating the predicted library in a separate step is **benign** (per its author) — don't chase it. |
| `spark_heap` | `java.lang.OutOfMemoryError: Java heap space` / `Answer from Java side is empty` | JVM heap OOM in Fulcrum/Spark. **Raising `--mem` does not help** — Spark sizes its driver heap independently of the SLURM allocation. Set `spark.driver.memory` via the workflow TOML's `spark_config`. |
| `empty_results` | a job that FAILED and whose log says `0 proteins` (as a count, not the end of `=1340 proteins`) / `No fragment ions` / `No precursors` / `no spectra`; never a completed job, never the bare word "empty" | FASTA/organism mismatch, mass-accuracy, or wrong acquisition type — check inputs |
| `gpu` | CUDA / no kernel image | AlphaDIA needs a GPU: submit to a GPU node (`--gres=gpu:1`) or reduce batch size |
| `sage_no_mzml` | `FAILED: no mzML from …` / msconvert not found / needs mzML | the job's own conversion failed (it converts `.raw` with ThermoRawFileParser and `.d` with msconvert before Sage/Radiant). Read the converter's error above the FAILED line: for ThermoRawFileParser usually .NET (`bash scripts/ensure_dotnet8.sh`); for `.raw` through msconvert, a Linux build without vendor readers (bioconda's). Fix it, regenerate with `run_search.py --sbatch`, resubmit |
| `disk` | Disk quota / No space left | free space or point `--out` elsewhere; resubmit |
| `missing_input` | No such file / fasta/raw not found | a path is wrong — re-check Windows→WSL/HIVE path translation; resubmit |
| `stalled` | job **RUNNING** but log frozen > `--stall-min` (default 15 min) | a hung file. `scancel` the task/job, **retry it once on a fresh node**; if it stalls again, **stop and ask the user** — never leave the file out yourself (see playbook) |
| `held` | job PENDING with reason **`JobHeldUser`** / **`JobHeldAdmin`** | it will **never start by itself**, and it is **not** failed — do **not** resubmit (that duplicates it). `JobHeldUser`: `scontrol release <jobid>` once whatever it was held for is fixed; `JobHeldAdmin`: ask the HIVE admins. (A PENDING job otherwise reports "queued, not started" only when nothing has run yet; a chain step waiting on earlier steps keeps its file count.) |
| `dependency_failed` | job PENDING with reason **`DependencyNeverSatisfied`** (an upstream step failed) | the chain is dead but **sits PENDING forever** (never leaves the queue). `sacct -j <arrayjob>` to find the failed step, fix it, resubmit downstream steps reusing completed outputs. If it was **step 1b**, `scancel` the pending steps 2–5 (ids in `<out>/jobs.txt`), then resubmit `step1b_window.sbatch` and steps 2–5 chained `afterok` on its new id (reusing `step1.predicted.speclib`), or re-run `submit.sh`: steps 2–5 refuse to start without the `window.txt` / `massacc.txt` step 1b writes (without them DIA-NN would optimise per file). When only the probe's own machinery failed (its log unreadable — e.g. an NFS `Stale file handle` — a crash, or a time limit), step 1b retries a crash or an unreadable log once. **From 2.15.1 a probe that measures the window then fails** (the chain stops: resubmit, or pin a `--window` measured on an earlier search with the same LC method); 2.9–2.15.0 fell back to `window.txt` = `auto` instead. A probe that measures only the mass accuracy (`--window` pinned) still **falls back** to the documented level for every step (`<out>/probe_fallback.json`, mode `fallback_default`), and the chain runs on — `watch_run.sh --all` then reports `probe_fallback: true` with a CAUTION: say so to the user, with its `reason`. It still fails, as before, on everything else — the job log names which: a **refused** mass accuracy (implausible: check the FASTA, species and calibration first), no .NET / DIA-NN cannot start, the probe's arguments, runs DIA-NN finished without logging a radius, a signal |
| `unknown_failure` | job FAILED with no known signature | read the full log; diagnose; fix; resubmit (and add the new signature here) |
| `watcher_query_failed` | the watcher could not ask HIVE (unreachable, `HIVE_USER`/`HIVE_KEY` unset, or — with no local python3 — the skill not on HIVE: exit 3), or for `--slurm <id>` sacct and squeue both answered nothing (`state: query_failed`: SLURM tools not on PATH, or a wrong or aged-out id). A chain (`--all` / `--out` with `jobs.txt`) says the same for the ids neither answers for (`unknown_jobs`; `state: query_failed` when it is every id), and is never `done` while one is unknown. Always `failed: false` | **not a job failure: never resubmit on it.** Fix the connection (or `bash scripts/hive_exec.sh --put-skill`), or check the id with `hive_exec.sh 'sacct -j <id> -X'`, then poll again |
| `watcher_no_python` | no usable `python3` on this computer and no `--hive` | **not a job failure** (`failed: false`). Install Python 3, or for a search on HIVE add `--hive` |

**No python3 on this computer** (Windows, where `python3` is often the Microsoft Store stub
that prints "Python was not found"): with `--hive` the watcher runs the whole poll on HIVE,
through `hive_exec.sh`, and prints HIVE's JSON — the skill must be on HIVE (step 0 puts it).
Before 2.13.2 it printed `failed: true, done: true` for a healthy search. With `--hive`, Git
Bash's path conversion is off (`MSYS_NO_PATHCONV=1`), so a `/nfs/...` folder is never
rewritten into `C:/Program Files/Git/nfs/...` (which read as `no_jobs_file`).

## Warnings that are not failures: proteins the contaminant list took

Every DIA-NN search's last job runs `contaminant_competition.py` (`references/audit.md`). Once
the search has finished, the watcher reads `<out>/contaminant_competition.json` and adds
`contaminant_competition` (`status`, `competed_out`, `siphoned`, `drop_contaminant`) and, when a
gene was taken, its message to `warnings`. A record whose report has since changed (a re-search
into the same folder) reads `status: stale` and is not shown; `unchecked` when the report could
not be looked at (poll again). Not a failure, nothing to
resubmit: tell the user which genes.

## Warnings that are not failures: a DIA-NN --dda search on its first pass

A DIA-NN `--dda` search whose second pass (match-between-runs) crashes, but which wrote its
complete first-pass report (`report-first-pass.parquet`), falls back to it in the job
(`dda_first_pass.py`): the first pass is copied onto `report.parquet`, and
`search_provenance.json` says `report_used: "first_pass"`. The job then succeeds. With
`--out` (single job) or `--all` (the two-job search), the watcher reports it `done`, not
`failed`, with `dda_first_pass: true` and the caveat in `caution` (single job: also in
`warnings`). Nothing to resubmit: tell the user, quoting the caveat; the Methods and AUDIT.md
carry the same words. A crash with no first pass still fails the job, as before.

## Warnings that are not failures: Sage LFQ

A Sage search can finish cleanly and still have quantities it could not measure. With
`--out <search out dir>` the watcher reads `quantification_check.json` (2.14.x: `sage_lfq_check.json`; written by the Sage job,
and again by `run_search.py --adapt-only`); without it, it reads the log tail. It adds:

- `sage_lfq`: `status` (`ok` / `warn` / `unchecked`), `target_ms1_peaks_5pct_fdr` (Sage's own
  count, "discovered N target MS1 peaks at 5% FDR"), `ppm_tolerance`,
  `runs_outside_lfq_window`, `suggested_ppm_tolerance`.
- `warnings`: the check's message when it warns (or could not run).

This is **not** a failure: `failed` and `error_class` are unchanged, and there is nothing
to resubmit automatically. Tell the user before any DE on those quantities, quoting the
message. It names the fix (a wider `quant.lfq_settings.ppm_tolerance`, then
`--adapt-only`). Re-run Sage only once they agree: a changed LFQ setting is a changed
search.

## Warnings that are not failures: protein inference in step 5

Every 5-step chain's step 5 log (and its `report.log.txt`) prints, on DIA-NN 2.7.0:
`WARNING: protein inference is enabled while empirical library is used; make sure the library
does not contain already inferred protein groups or disable protein inference`.
**It is expected and does not matter here; do not pass `--no-prot-inf`.** Checked on HIVE on a
real chain (HIVE, 2026-10-16, read-only): step 3's `empirical.parquet` keeps each precursor's full
`Protein.Ids`, the same as the first-pass report for all 65,597 precursors both hold, and 2,986
of 3,000 sampled precursors map to exactly those proteins in the FASTA (the other 14 also match,
as a bare string, proteins the library does not list: a string search ignores where trypsin
cuts; not checked one by one). Its `Protein.Group` column is the inferred group, but the
library has not lost the mapping step 5 infers from. Step 5 run with `--no-prot-inf` instead
gave each precursor the same protein group but a different protein-level FDR (9,486 rows in
the protein-group matrix against 9,414). Tell staff it is benign, in one line, if they ask.

## SLURM specifics (HIVE)
- State comes from `sacct -j <jobid> -o State` (falls back to `squeue`, in chain mode too).
  What each state means is `scripts/slurm_states.py` (also read by `record_run.py` and
  `checkpoint.py`): COMPLETED (done); FAILED / TIMEOUT / OUT_OF_MEMORY / NODE_FAIL / CANCELLED
  / PREEMPTED / DEADLINE / BOOT_FAIL / REVOKED / SPECIAL_EXIT (done + failed); PENDING /
  RUNNING / REQUEUED / SUSPENDED and the other live states (not done). A state the list does
  not know is never read as done.
- **Preempted, then requeued.** `sacct -j` without `-D` returns only a job id's most recent
  record, and the watcher never passes `-D`. A job preempted on the low queue and requeued by
  Slurm therefore reads as its new instance (PENDING, then RUNNING; `REQUEUED` while Slurm puts
  it back: `progress.summary` says it will run again): not done, not failed, nothing to do.
  `PREEMPTED` in that newest record means Slurm did not requeue it: `error_class: preempted`,
  resubmit the step on the high queue.
- The job log is the sbatch `--output` file (the skill's `emit_sbatch` names it
  `<out>/<job>_<jobid>.log`).
- A fix that changes memory, time or nodes is **re-generated** with `run_search.py` (the
  `out_of_memory` / `timeout` rows above), never made by editing the sbatch by hand; the
  command goes into `commands.log` as usual.

## Node faults: the node check every chain job runs first

Every job of the DIA-NN chain starts with a **node check** (`node_fault.preflight_lines`, run
inside the job-end hook): the DIA-NN binary (or image), the search folder (writable), the FASTA
and the raw-data folders must answer on this node within 120 s. They all existed on the login
node when the chain was generated, so one this node cannot see — a dropped `/quobyte` FUSE mount
(`Transport endpoint is not connected`, hive-dc-7-4-50, 2026-09-28), or a mount that never answers
(16 min on hive-dc-7-4-54 the same day) — ends the job with exit **75** and a
`NODE_FAULT: node=<n> check=<what> path=<p> detail=<why>` line, **before** DIA-NN runs, so the
failure is never reported as "DIA-NN exited 0 but did not write …". The job-end post says
"the compute node could not reach the storage". It is not the probe's ESTALE handling: step 1b's
probe rereads its own log on `Stale file handle` and retries; a mass-accuracy-only probe then
falls back (`probe_fallback.py`) and that job succeeds, while a window probe fails (2.15.1). Only a
failed job is classified.

## Keeping a search off nodes known to be bad (`--exclude`)

When the user, the HIVE status page or an earlier failure says which nodes are bad, keep the
whole search off them **when it is generated**:
`run_search.py ... --exclude hive-dc-7-4-50,hive-dc-7-4-[52-54]` (names, or SLURM ranges in
brackets). Every job script it writes — both single-shot jobs, every step of the DIA-NN chain,
every step of the Radiant chain, the one Sage/AlphaDIA/FragPipe job — gets
`#SBATCH --exclude=<the nodes>`, so even a job that starts before anyone is watching cannot
land on them. The list is checked first (anything that is not a node name or range is refused
and nothing is written) and recorded in `search_provenance.json` as `queue.exclude`.
`diann_parallel.py` and `radiant_parallel.py` take the same `--exclude` when run by hand.

- **`export SBATCH_EXCLUDE=...` does nothing.** sbatch reads no such variable — it is not in the
  list of input environment variables in `man sbatch` (Slurm 26.05.4 on HIVE) — so the jobs
  show `ExcNodeList=(null)`. Checked on HIVE 2026-10-05: `SBATCH_EXCLUDE=<no such node> sbatch
  --test-only` was accepted, `--exclude=<no such node>` was refused.
- **A misspelt node name stops the submission.** sbatch refuses an unknown name ("Invalid node
  name specified"), so `submit.sh` stops at its first job and nothing is queued.
- **Jobs already queued** are not changed by a new `--exclude`: for those,
  `scontrol update JobId=<id> ExcNodeList=<nodes>` while they are pending.
- **It works with `node_fault.py retry`.** A `--exclude` on the sbatch command line *replaces*
  the script's own line (the command line overrides the script), so the retry passes both: the
  nodes the search was kept off and the node that failed, each once. `node_faults.json` records
  them as `excluded`, with `excluded_by_search` for the ones from the script.

## Auto-recovery playbook (apply autonomously; log every action to `commands.log`)
- **`node_fault`:** `node_fault.py retry --out <out> --job <id>` (above). Tell the user in one
  line that a compute node could not reach the storage and the step was resubmitted on another
  node — the `say` field. Not a retry of the same class against the 2-attempt limit below: the
  node retry is bounded by its own (2 per step).
- **`stalled` (hung file):** `scancel <arrayjob>_<taskid>` (or the job). Retry that one file
  once on a fresh node. **If it stalls again, stop and ask the user** (SKILL.md golden rule 9:
  runs left out are a person's decision). Tell them which file hung twice and what it costs to
  leave it out (one sample fewer in its group). Never leave it out on your own, and never edit
  `file_list.txt`, a step's `--f` list or the run count in `step5_report.sbatch`: the chain's run
  list is fixed when it is generated, so step 5 fails on an edited one, and the record would
  no longer match the search. What the user can choose:
  - **Leave the file out.** Record their words (`log_decision.py`), then generate the search
    again with `run_search.py` on the remaining files (a new chain; tell them it starts again
    from step 1). At step 8 the run has no row in `conditions.csv`.
  - **Check the file first** (copied completely? readable?), or ask Core staff, then retry.
- **Broken `afterok` after you killed/dropped an array task:** the downstream steps go
  `DependencyNeverSatisfied`. Cancel them and **resubmit steps 3→4→5 fresh, reusing the
  completed `.quant`** (don't re-run step 1/library-prediction — reuse `step1.predicted.speclib`).
- **`out_of_memory`:** raise `--mem` and resubmit only the failed step (reuse prior `.quant`).
- **Missing `--temp` dir** (`cannot find the temp folder`): the temp dirs must pre-exist —
  `mkdir -p <out>/quant_step2 <out>/quant_step4` and resubmit. **Every generated step now does
  this itself**, so a step resubmitted directly (as this playbook tells you to) cannot hit it;
  you should only ever see this from a hand-written sbatch or a cfg-supplied `--temp` on the
  single-shot path (which `run_search.ensure_temp_dirs` also creates).
- **Never restart the whole cohort** — resume from the earliest incomplete step,
  reusing everything already computed.

## Don't

- Don't declare a run "done" on a non-COMPLETED state, and don't proceed to DE until
  `report.parquet` exists.
- Don't silently retry forever — after 2 failed auto-fix attempts of the same class,
  stop and tell the user what's wrong. A file that stalls twice is the user's decision, not a
  third retry and not a file to leave out yourself (playbook above).
- Don't wait for the user to notice a stall/failure — you are the monitor; recover, then
  report what you did.
- **Don't use a passive "wait until the job leaves the queue" monitor for a chain.** If an
  upstream step fails, the downstream job sits **PENDING forever** (`DependencyNeverSatisfied`)
  and such a monitor hangs indefinitely and never alerts — a whole overnight chain can die
  silently. The monitor must poll **STATE** (incl. `squeue -o %r` reason) and exit/act on
  `failed` / `dependency_failed` / `stalled`, not just on queue-exit.
- **A job's SLURM State=COMPLETED does NOT mean the tool succeeded.** If the sbatch's last
  line is `echo ...` (exit 0), FragPipe/DIA-NN can crash yet the job shows COMPLETED. Verify
  the expected OUTPUT file exists (`report.parquet`, `combined_peptide.tsv`), not just the state.

## Progress reporting (`--out`, `--poll`)

`--out <search out dir>` makes the watcher report **real** progress instead of "still
running". It counts what the 5-step chain leaves on disk — one `.quant` per finished file —
and infers the stage from which artefacts exist:

| Present in `--out` | Reported stage |
|---|---|
| nothing yet | step 1/5, library prediction |
| `step1.predicted.speclib` | step 2/5, first pass — progress = `quant_step2/*.quant` ÷ `file_list.txt` |
| all first-pass `.quant` present | step 3/5, empirical assembly |
| `empirical.parquet` | step 4/5, final pass — progress = `quant_step4/*.quant` |
| `report.parquet` | step 5/5, done |

Emptiness counts as absent (`test -s`), so a zero-byte or truncated output is not mistaken
for a finished step. All checks run through the same `run()` helper, so they work over SSH
on HIVE exactly as locally. Without `--out` the watcher behaves as before.

`--poll <n>` is a counter you increment each poll. It rotates the stage explanation's
companion note (`pipeline_notes.py`) so a multi-hour search does not repeat itself; notes
relevant to the current stage come first, then the rest of the list.

**Relay `progress.summary` + `doing` to the user every poll**, and pass `note` along with
its `note_source`. Silence during a multi-hour search reads as a hang. The note is general
background — never present it as a finding about the user's data.

Every fact in `pipeline_notes.py` carries a `source`. If you add one and cannot attribute
it, do not add it.

`checkpoint.py record` flags are written by the generated `submit.sh`: `--jobs` (the SLURM ids),
`--desc`, `--next` (the command to run once the stage completes), `--watch-job` / `--watch-log`
(the job and log the watcher follows for that stage). `check_report_runs.py --files-list <file>` gives the input paths
one per line; `pipeline_notes.py --index <n>` rotates the progress notes (`--no-fact` leaves out
the proteomics fact).

## 0c. Is this a resume? (step 0c)

A cluster search runs for hours and the user **will** close their laptop. The SLURM
jobs survive that; the conversation does not. So before doing any work, check whether
an earlier session already submitted something (step 0's version check,
`skill_version.sh --check-update` and in `hive_remote` `--check-hive`, runs on a resume too,
first; then its notes check, with `--session <dir>` once the session is found):

```
python3 scripts/checkpoint.py find --base <the project folder>   # never ~ or a whole share; HIVE: in the background
python3 scripts/checkpoint.py status --session <session dir> --resume   # live sacct re-query
```

`status` returns each stage's real state, whether the expected output file exists, and
`next_commands` — the exact commands still outstanding.

**A run found here is a resume only when it is unfinished work the user is coming back to**
("pick up where we left off", "is my search done?"). A new request (new settings, or "analyse
these files") is fresh even when an unfinished run exists: mention it in one line, and still
never resubmit one that is RUNNING without asking. A finished earlier analysis is never a resume
(step 1b); nothing is reused from either unless the user asks. For a resume:

- **Still RUNNING** → say so, give the job ids, and offer to wait (`watch_run.sh`) rather
  than resubmitting. **Never resubmit a search that is already running** — it wastes hours
  of cluster time and the user's budget.
- **COMPLETED and the output exists** → skip the search entirely and continue from
  `next_commands` (usually DE onward).
- **FAILED or stalled** → apply the `watch_run.sh` playbook (step 7b) and resubmit only the
  incomplete steps; finished `.quant` files and the predicted library are reused.

**Core staff: a resume asks no submission question.** `status` prints the session's
`analysis_for` (the submissions or the internal work's board title, from `session.py init`);
`bash scripts/board.sh threads` (no usable local Python: `bash scripts/hive_exec.sh 'bash
~/proteomics-pipeline/scripts/board.sh threads'`) gives its thread for the next milestone's
record (`references/board.md` section 3).

## A closed laptop: `RECOVERY.md` (step 0c)

`submit.sh` writes `RECOVERY.md` + `.recovery.json` into the session directory
automatically (a search whose `--out` is not inside a session gets them beside `--out`, and
`submit.sh` says where; `node_fault.py retry` updates them in the same place), so this works
even if the previous session ended mid-sentence. Tell the
user, in plain words, that they can close their terminal and that the jobs will keep
going — they do not need `tmux` or `screen` for this.

## 7b. The watch loop as SKILL.md runs it (step 7b)

**The moment a search starts, watch it — never leave a multi-hour search
unmonitored.** Poll with `watch_run.sh` in a loop until it finishes:

```
# local: python3 scripts/run_search.py ... run_in_background, then loop:
bash scripts/watch_run.sh --log <search log> --out <search out dir> --poll <n>
# SLURM on HIVE:
bash scripts/watch_run.sh --slurm <jobid> --log <job log> --out <out> --poll <n> --hive
# multi-step chains (diann_parallel / radiant_parallel): watch EVERY step, never just the last
bash scripts/watch_run.sh --all <search out dir> --hive
```

With `--hive` and no usable `python3` on this computer (Windows' Store stub), the watcher
runs the poll on HIVE and prints HIVE's answer. `watcher_query_failed` and
`watcher_no_python` are never a job failure: never resubmit on them; fix the connection and
poll again (this file).

**Always pass `--out` and increment `--poll` each time.** `--out` is what turns "still
running" into real progress: the chain drops a `.quant` per finished file, so the watcher
reports the actual step and file count rather than guessing. `--poll` rotates the note.

**Report progress to the user on every poll — a search is hours long and silence reads as
a hang.** Give them one short line from `progress.summary` plus `doing`, and pass along
`note` (with `note_source`) so there's something to read while waiting. For example:

> Step 2 of 5 — searching each file against the predicted library. 34 of 66 files done
> (52%). *While you wait: decoy sequences exist to be wrong on purpose. How often a search
> prefers a decoy is what calibrates the false discovery rate. (Elias & Gygi 2007)*

Keep it to a couple of lines; don't dump the JSON. Never present the note as a finding
about **their** data — it is general background, and it must not be mistaken for a result.
It returns `{state, done, failed, stalled, error_class, fix}`. While `done` is false,
keep polling: about every minute for a local run; **on SLURM every 5 minutes or more**, or a
scheduled wake-up (each `--hive` poll is several logins and `sacct` calls on HIVE's shared login
node, and a stall is only called after 15 minutes anyway).
Pass `--log` so it can also catch **stalls** — a job stuck in `RUNNING` with its log
frozen (a hung file; `sacct` will show RUNNING forever). **If `failed` OR `stalled` is
true, apply the `fix` and resubmit AUTONOMOUSLY — do not wait for the user.** Starting
the search needed confirmation (golden rule #1); recovering a run that's already in
flight does not. **Except** (each of these is told to the user, not resubmitted):
`cancelled` (find out who cancelled it first: parallel sessions share one HIVE account),
`held`, `needs_decision`, `watcher_query_failed`, a search that failed its run check
(`check_report_runs.py`: never resubmitted, the user decides about the missing run), a file
that stalled twice, a third `node_fault` retry, and a third fix of one class. A stalled job
is RUNNING but hung: `scancel` it first, then retry; a RUNNING job that is progressing is
never resubmitted. Common auto-fixes (→ the playbook above):

- **stalled/hung file:** `scancel` that task, retry once on a fresh node; if it stalls
  again, **stop and ask the user** (the playbook above): never leave the file out yourself.
- **failed step in the 5-step chain:** resubmit from the earliest incomplete step,
  **reusing completed `.quant` and `step1.predicted.speclib`** — never restart the whole
  cohort; broken `afterok` after a killed task → resubmit steps 3→4→5 fresh.
- **`node_fault`** (a node could not reach the storage: exit 75 from the job's node check,
  `NODE_FAIL`, or `Transport endpoint is not connected` / `Stale file handle` in a failed
  job): tell the user it is a **node problem, not DIA-NN and not their data**, and run
  `node_fault.py retry --out <out> --job <id>` on HIVE — it resubmits that step and the steps
  waiting on it with `--exclude=<node>`, at most twice per step.
- **Nodes already known to be bad** (the user says so, the HIVE status page, an earlier
  failure): generate the search with `run_search.py ... --exclude <node[,node...]>` (ranges
  like `hive-dc-7-4-[50-54]` work). Every job script then carries `#SBATCH --exclude`, so no
  job can land there, and `search_provenance.json` records it. Do not `export SBATCH_EXCLUDE`:
  sbatch ignores it. A later `node_fault.py retry` keeps these nodes excluded as well
- OOM → raise `--mem`; timeout → raise `--time`; missing temp dir → `mkdir -p` first. A single-job search takes both from `run_search.py --mem GB --time HOURS` (recorded as `job_mem`/`job_time`); a non-specific or semi-specific Sage or FragPipe search already defaults to its broad-digest memory and time (`run_search.BROAD_DIGEST_SAGE_MEM_GB`, `BROAD_DIGEST_FRAGPIPE_MEM_GB`, `BROAD_DIGEST_HOURS`), recorded as that default — with the database's sequence count known (fetch_fasta's `.meta.json`) the memory is computed from it (`run_search.broad_digest_mem_gb`, from the measured points `SAGE_NONSPECIFIC_PEAKS`; recorded as `job_mem_estimate_gb` / `job_n_sequences`, and `job_mem_asked_gb` for your `--mem`), and an estimate over the queue's per-job cap is refused with the number — otherwise raise `--mem` for a database bigger than the one they were measured on (`BROAD_DIGEST_SEQUENCES`; `BROAD_DIGEST_MEASURED` says how). Every Sage job exports `RUST_MIN_STACK`; a stack overflow (exit 134) writes `sage_failure.json` with a re-run hint. A semi-specific FragPipe job's memory follows the measured peaks of semi-specific searches over a 12,685-sequence database (`FRAGPIPE_SEMISPECIFIC_PEAKS`: 23 runs 196 GiB, 36 runs 184, 58 runs 230, 183 runs 251), x1.15, rounded up to 8 (`run_search.fragpipe_broad_mem_gb`: 232G for 8 to 36 runs, 272G for 58, 296G for 183, 312G for 300; the job's memory is the larger of that and the database-sized estimate); a non-specific or tryptic one grows by `FRAGPIPE_MEM_GB_PER_RUN` a run past `FRAGPIPE_MEM_RUNS_BASE` runs (`run_search.fragpipe_runs_mem_gb`; a lower-bound estimate, the 200G default was full at 219 runs). The record says which numbers were measured. A FragPipe job's time follows its runs (`run_search.fragpipe_run_hours`: `FRAGPIPE_MINUTES_PER_RUN` a run × `FRAGPIPE_RUN_HEADROOM`, scaled up under `FRAGPIPE_RATE_CPUS` CPUs; the larger of that and the default), recorded as that default.
- **A failed search is recorded in the Core run registry. When the chain dies before its last
  job, record it yourself.** A job that SLURM kills (OOM, TIMEOUT, node failure, `scancel`)
  never runs its own job-end hook. The jobs queued after it never start at all. So whenever the
  watcher reports `failed`, run this on HIVE, then fix and resubmit. Recording again is safe:
  it updates the same record.
  `python3 scripts/record_run.py search-done --out <search out dir> --status failed
  [--exit-code <N>] [--step <failed step's job name>]`.
  → `references/run-registry.md`

**Every search — single-shot or parallel array — must run under this watch loop.** Stop
and tell the user only after 2 failed auto-fixes of the same class. Report what you
recovered. Only proceed to DE
once `COMPLETED` and `report.parquet` exists.

**5-step chain: check what the final pass kept before DE.** Step 5 compares its report with the
chain's first pass (step 3's `step3_assembly.parquet`) run by run, prints the table in its log,
and writes it to `<out>/pass_comparison.json` and `search_provenance.json` (`pass_comparison`).
A run whose final pass kept less than half its first-pass precursors is flagged with a
`WARNING`. DIA-NN's README says to compare the two and use the first pass when it performs
better. **Those counts are before protein-group FDR; compare the post-filter rows before
switching.** The table's "rows after q" columns count each report's rows that pass
`run_de.R`'s q-value filter (the same definition, from `diann_q_columns.py`; they need pyarrow,
and read "n/a" without it). On SET28 (hair DDA) 4 of 6 runs fell by more than half (126 -> 0,
165 -> 1). At 1% on every q-column those four had 0 rows in BOTH reports; with `run_de.R`'s own
filter, which holds PG.Q.Value at 5%, the first pass kept 109–160 rows for them against 0–9 --
protein groups with a run-specific q between 1% and 5%, all of the gain. Tell the user where a
gain comes from. **Show the table, and offer the first-pass report only for the runs
`switch_recommended_for` names** -- those with more post-filter rows in the first pass:
`run_de.R --input <out>/step3_assembly.parquet …`. That report is each run searched
once against the predicted library, with no match-between-runs, and it has no DDA MS1 q columns.
If they take it, `make_methods.py` says so in the Methods, from `de_provenance.json`'s `input`
(never write that by hand), and `audit_results.py` reports the comparison in AUDIT.md either way.

**Every DIA-NN search: the proteins the contaminant list took.** Its last job runs
`contaminant_competition.py` (a report, never a gate); the watcher's `contaminant_competition`
names COMPETED OUT genes (no quantity left) and SIPHONED ones. Tell the user; AUDIT.md carries
it to the report. → `references/audit.md`.

Core staff: board.md §3 row "search finished", once its last job has ended.
