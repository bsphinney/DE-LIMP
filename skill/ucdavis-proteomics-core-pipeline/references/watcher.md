# Run watcher — monitor searches and auto-correct errors

A search is long (often hours) and runs detached (a background process locally, a
SLURM job on HIVE). **Always watch it** — the orchestrator must not start a search
and walk away. `watch_run.sh` does one poll + diagnosis; loop it until the run
finishes, and on failure apply the fix and resubmit.

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
| `empty_results` | 0 proteins / no fragment ions | FASTA/organism mismatch, mass-accuracy, or wrong acquisition type — check inputs |
| `gpu` | CUDA / no kernel image | AlphaDIA needs a GPU: submit to a GPU node (`--gres=gpu:1`) or reduce batch size |
| `sage_no_mzml` | msconvert not found / needs mzML | convert `.d`/`.raw` → mzML first (Linux/HIVE), then re-run Sage |
| `disk` | Disk quota / No space left | free space or point `--out` elsewhere; resubmit |
| `missing_input` | No such file / fasta/raw not found | a path is wrong — re-check Windows→WSL/HIVE path translation; resubmit |
| `unknown_failure` | job FAILED with no known signature | read the full log; diagnose; fix; resubmit (and add the new signature here) |

## SLURM specifics (HIVE)
- State comes from `sacct -j <jobid> -o State` (falls back to `squeue`). Terminal
  states: COMPLETED (done), FAILED / TIMEOUT / OUT_OF_MEMORY / NODE_FAIL / CANCELLED
  (done + failed).
- The job log is the sbatch `--output` file (the skill's `emit_sbatch` names it
  `<out>/<job>_<jobid>.log`).
- When resubmitting after a fix, edit the sbatch (`--mem`/`--time`/`--gres`) and note
  the change in `commands.log` so the reproducibility bundle records what happened.

## Don't
- Don't declare a run "done" on a non-COMPLETED state, and don't proceed to DE until
  `report.parquet` exists.
- Don't silently retry forever — after 2 failed auto-fix attempts of the same class,
  stop and tell the user what's wrong.
