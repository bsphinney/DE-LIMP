# Run watcher — monitor searches and auto-correct errors

A search is long (often hours) and runs detached (a background process locally, a
SLURM job on HIVE). **Always watch every search — this is mandatory, not optional.**
The orchestrator must never start a search and walk away, and must **auto-correct
stalls and failures on its own, without waiting for the user** (the user has standing
authorization for this — see golden rule #1's compute confirmation is for *starting*
the search, not for recovering a run that's already approved and in flight). Surface
what you did afterwards; don't ask permission to unstick a job.

## Watch the WHOLE chain, not its last job

For any multi-step chain (`diann_parallel.py`, `radiant_parallel.py`), watch with:

```bash
bash scripts/watch_run.sh --all <search out dir> --hive
```

It reads `jobs.txt` (written by `submit.sh`) and reports every job's state at once.
**Do not** watch only the final job. If an earlier step dies, the final job sits
`PENDING` on a dependency that can never be satisfied — which is indistinguishable
from "still queued" if that is the only job you are looking at. This actually
happened: a step-4 array hit `TIMEOUT`, step 5 went `CANCELLED`, and a watcher
pointed at step 5 alone reported "pending" for hours. `--all` returns
`failed:true` with `first_failed_job` in that situation.

`watch_run.sh` does one poll + diagnosis (state, **stall**, error class, fix); loop it
until the run finishes, and on `failed` **or** `stalled` apply the fix and resubmit.

## Two failure modes it catches
1. **Hard failure** — `sacct` reports FAILED / TIMEOUT / OUT_OF_MEMORY / NODE_FAIL /
   CANCELLED. Matched against known error signatures → a fix.
2. **Stall** — the job stays **RUNNING** but its log stops advancing (a hung file, e.g. a
   pathological DIA-NN run whose RT-window/search phase never returns). `sacct`/`squeue`
   will happily show RUNNING forever. `watch_run.sh --log <log> [--stall-min N]` flags
   this when the log's mtime is older than N minutes (default 15) while RUNNING, emitting
   `stalled: true`. **This is the NA41-class failure the first nail cohort hit.**

## Loop pattern
```
# start the search (background locally, or sbatch on HIVE), capture the log path / job id
while true; do
  status=$(bash scripts/watch_run.sh --slurm <jobid> --log <log> --hive)   # or --log <log> locally
  done=$(echo "$status" | python3 -c 'import sys,json;print(json.load(sys.stdin)["done"])')
  [ "$done" = "True" ] && break
  sleep 60        # for long SLURM jobs, schedule a wake-up instead of busy-waiting
done
# if failed: read error_class + fix, apply it, resubmit, watch again.
```
The agent itself is the watcher — `watch_run.sh` is its eyes. Surface progress and
any auto-fix you applied to the user.

## Error classes → fixes (what `watch_run.sh` detects)
| error_class | signal | fix |
|---|---|---|
| `out_of_memory` | oom-kill, `std::bad_alloc`, OUT_OF_MEMORY | raise sbatch `--mem` (e.g. 64G→128G); DIA-NN: fewer threads; resubmit |
| `timeout` | TIMEOUT / "DUE TO TIME LIMIT" | raise `--time`, or split the run; resubmit |
| `diann_no_dotnet` | `dotnet: not found` | wrong DIA-NN container (no .NET → `.raw` silently skipped). Use the HIVE **native** build `build_<v>/diann-<v>/diann-linux` or a .NET image |
| `diann_zero_ids` | `Number of IDs at 0.01 FDR: 0` | Completed but identified nothing — a **silent null result**, not a crash. **Preserve `report.log.txt` before re-running**; it is the only evidence. Check: library generated precursors → mzML really are DIA with isolation windows (`detect_acquisition.py`) → FASTA matches organism → mass accuracy. DIA-NN's warning about generating the predicted library in a separate step is **benign** (per its author) — don't chase it. |
| `spark_heap` | `java.lang.OutOfMemoryError: Java heap space` / `Answer from Java side is empty` | JVM heap OOM in Fulcrum/Spark. **Raising `--mem` does not help** — Spark sizes its driver heap independently of the SLURM allocation. Set `spark.driver.memory` via the workflow TOML's `spark_config`. |
| `empty_results` | 0 proteins / no fragment ions | FASTA/organism mismatch, mass-accuracy, or wrong acquisition type — check inputs |
| `gpu` | CUDA / no kernel image | AlphaDIA needs a GPU: submit to a GPU node (`--gres=gpu:1`) or reduce batch size |
| `sage_no_mzml` | msconvert not found / needs mzML | convert `.d`/`.raw` → mzML first (Linux/HIVE), then re-run Sage |
| `disk` | Disk quota / No space left | free space or point `--out` elsewhere; resubmit |
| `missing_input` | No such file / fasta/raw not found | a path is wrong — re-check Windows→WSL/HIVE path translation; resubmit |
| `stalled` | job **RUNNING** but log frozen > `--stall-min` (default 15 min) | a hung file. `scancel` the task/job, **retry it once on a fresh node**; if it stalls again, **drop that file** and continue (see playbook), note it in Data Quality Notes |
| `dependency_failed` | job PENDING with reason **`DependencyNeverSatisfied`** (an upstream step failed) | the chain is dead but **sits PENDING forever** (never leaves the queue). `sacct -j <arrayjob>` to find the failed step, fix it, resubmit downstream steps reusing completed outputs |
| `unknown_failure` | job FAILED with no known signature | read the full log; diagnose; fix; resubmit (and add the new signature here) |

## SLURM specifics (HIVE)
- State comes from `sacct -j <jobid> -o State` (falls back to `squeue`). Terminal
  states: COMPLETED (done), FAILED / TIMEOUT / OUT_OF_MEMORY / NODE_FAIL / CANCELLED
  (done + failed).
- The job log is the sbatch `--output` file (the skill's `emit_sbatch` names it
  `<out>/<job>_<jobid>.log`).
- When resubmitting after a fix, edit the sbatch (`--mem`/`--time`/`--gres`) and note
  the change in `commands.log` so the reproducibility bundle records what happened.

## Auto-recovery playbook (apply autonomously; log every action to `commands.log`)
- **`stalled` (hung file):** `scancel <arrayjob>_<taskid>` (or the job). Retry that one file
  once on a fresh node. If it stalls again, **drop it** — it's pathological (often a file
  the facility itself re-ran). In the 5-step parallel chain, step 4 already auto-skips a
  file with no step-2 `.quant`; just remove it from step 3/step 5's `--f` list and resume.
- **Broken `afterok` after you killed/dropped an array task:** the downstream steps go
  `DependencyNeverSatisfied`. Cancel them and **resubmit steps 3→4→5 fresh, reusing the
  completed `.quant`** (don't re-run step 1/library-prediction — reuse `step1.predicted.speclib`).
- **`out_of_memory`:** raise `--mem` and resubmit only the failed step (reuse prior `.quant`).
- **Missing `--temp` dir** (`cannot find the temp folder`): the temp dirs must pre-exist —
  `mkdir -p <out>/quant_step2 <out>/quant_step4` and resubmit (the generator now does this).
- **Whole-cohort restart is rarely needed** — resume from the earliest incomplete step,
  reusing everything already computed.

## Don't
- Don't declare a run "done" on a non-COMPLETED state, and don't proceed to DE until
  `report.parquet` exists.
- Don't silently retry forever — after 2 failed auto-fix attempts of the same class,
  stop and tell the user what's wrong. (Dropping 1 pathological file out of many and
  continuing is *success*, not a failed retry.)
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
