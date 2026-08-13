# DIA-NN parallel search (5-step SLURM chain)

DIA-NN's high-throughput workflow — **poorly documented upstream**, so it is encoded
here (`diann_parallel.py`), ported faithfully from DE-LIMP's `generate_parallel_scripts()`.

## When it is used — automatic above 5 files

`run_search.py` decides this for you; you rarely call `diann_parallel.py` by hand. It
routes to the chain when **all** of these hold, and prints the reason either way:

| Condition | Why it's required |
|---|---|
| engine is **DIA-NN** | the chain is DIA-NN-specific |
| **more than 5 files** (`--parallel-threshold`, default 5) | below that, chain overhead (library prediction + two array round-trips) outweighs the win |
| **SLURM present** (`sbatch` on PATH) | steps 2 and 4 are job arrays — there is no non-cluster equivalent |
| **mass accuracy fixed** in the `--cfg` | steps 3/5 reuse the `.quant` files from 2/4; auto-calibration would differ between passes and corrupt the cross-run report |

Any condition unmet → single-shot search, reason printed and written to
`search_provenance.json` (`search_mode`, `parallel_routing_reason`). The fixable one is
almost always mass accuracy: re-run `estimate_params.py` with the **real instrument** so
it pins the DIA-NN recommended values, and parallel enables itself. Force the decision
with `--no-parallel` or `--parallel-threshold N`.

`diann_parallel.py` **refuses to run** on an auto/0 mass-accuracy cfg rather than
producing a silently inconsistent report (`--allow-auto-mass-acc` overrides, for testing
only).

## The 5 steps (chained with `afterok` dependencies)
1. **Library prediction** — single job, no raw: predict a spectral library from the
   FASTA (`--fasta-search --predictor --gen-spec-lib --out-lib step1.speclib`).
2. **First pass** — SLURM **array**, one file per task: search each raw vs the
   predicted library (`--lib step1.predicted.speclib --temp quant_step2 --gen-spec-lib
   --quant-ori-names`) → per-file `.quant`.
3. **Empirical-library assembly** — single job: `--use-quant` over the step-2 `.quant`
   to build the empirical library (`--out-lib empirical.parquet`).
4. **Final pass** — SLURM array: search each raw vs the empirical library
   (`--lib empirical.parquet --temp quant_step4 --no-ifs-removal --quant-ori-names`).
5. **Cross-run report** — single job: `--use-quant --matrices` over the step-4
   `.quant` → `report.parquet` (the DE contract).

**Why it's faster:** the per-file passes (2 & 4) run as a SLURM array across many
nodes simultaneously instead of one long single-node job; MBR is replaced by the
empirical-library round-trip.

## Critical details (don't change these)
- **Mass accuracy is FIXED (manual), never auto.** Steps 3/5 reuse `.quant` files, so
  auto-calibration would be inconsistent (per DIA-NN dev guidance). So the `--cfg` you
  pass **must** have real `--mass-acc`/`--mass-acc-ms1` values — i.e. estimate params
  from a **known instrument** (timsTOF → 15/15, Astral → 4/10, Orbitrap by resolution;
  the DIA-NN-recommended table in `estimate_params.py`). Don't run parallel with an
  `--mass-acc 0` (auto) cfg.
- **No MBR** (`--reanalyse` is dropped) — the 5-step replaces it.
- **`--quant-ori-names`** on every step so `.quant` files are `<basename>.quant`.
- Step 4 **skips** files that failed step 2 (missing `.quant`).
- **XICs require `--out` on step 4, not just `--xic`.** DIA-NN names the XIC folder after
  the `--out` report basename. Step 4 otherwise has no `--out`, so DIA-NN resolves the
  folder against the filesystem **root** and aborts with
  `cannot create directory: Permission denied [/report_xic]` — *after* writing the
  `.quant`, and **still exiting 0**. Verified on DIA-NN 2.6.0: zero `.xic.parquet`
  produced, SLURM records `COMPLETED`, and the `afterok` chain advances as if XICs
  existed. Step 4 therefore gets `--out <out>/xic/t<TASKID>.parquet` whenever `--xic` is
  in the cfg, and DIA-NN writes `<out>/xic/t<TASKID>_xic/<run>.xic.parquet`.
  *Why step 4 and not step 5:* extraction needs the raw spectra. Step 5 runs
  `--use-quant`, which never re-reads them — DIA-NN accepts `--xic` there, logs that it
  will extract, and writes nothing.
- Step 4 hands each task **its own copy** of the empirical library
  (`libpriv/t<TASKID>/lib.parquet`, removed by a `trap` on exit). This is not an
  optimisation — DIA-NN re-saves the library it is given as `<lib>.skyline.speclib`,
  written *next to* `--lib`. Pointing every array task at one shared `empirical.parquet`
  makes them all write the same file; most win the race in seconds and the losers block
  until the wall clock kills them. Observed on a 399-file run: 51 tasks `TIMEOUT` at 3 h,
  then after raising the limit to 12 h, 8 tasks stalled for over 3.5 h at
  `[0:04] Saving the library`. With private copies the same 8 finished in 2 minutes.
  **Raising `--time-per-file` does not fix this** — it only moves the stall further out.
  Step 2 is unaffected because it is handed `step1.predicted.speclib`, already in
  DIA-NN's processed form, so there is nothing to re-save. Step 5 is a single job and
  cannot contend with itself.
- **The `--temp` folders MUST pre-exist.** DIA-NN aborts immediately with
  `ERROR: cannot find the temp folder .../quant_step2. Specify an existing folder` if
  the `--temp` dir is missing — it will **not** create it. `submit.sh` therefore does
  `mkdir -p <out>/quant_step2 <out>/quant_step4` before submitting. (This bit us on the
  first DDA cohort run — every first-pass array task died in ~7 s until the dirs existed.
  If you ever hand-edit or hand-run a step script, create the temp dirs first.)
- **DDA:** put `--dda` in the `--cfg` (it is not stripped, so it flows into every step);
  `.raw` inputs get a `.NET 8` export prefix automatically (see `ensure_dotnet8.sh`).

## Usage
```
python3 scripts/diann_parallel.py \
  --diann '<DIA-NN command from tools.json>' \   # e.g. the HIVE 2.6 native binary
  --raw /path/to/*.d --fasta search.fasta --out ./diann_parallel \
  --cfg params.cfg \                              # estimate_params.py output (known instrument!)
  --threads-per-file 16 --mem-per-file 64 --time-per-file 2 \
  --assembly-cpus 64 --assembly-mem 128 --assembly-time 12 \
  --partition high --account genome-center-grp --max-simultaneous 20
# then submit the chain (on HIVE, over hive_exec.sh):
bash scripts/hive_exec.sh 'bash <out>/submit.sh'
```
It writes `file_list.txt`, `step{1..5}_*.sbatch`, and `submit.sh` (which submits all
five with dependencies and prints the job ids). **Watch the final job** with
`watch_run.sh --all <out> --hive` (reads `jobs.txt`; covers all five steps — watching
only step 5 hides an upstream failure as a permanent `PENDING`), or for one step
`watch_run.sh --slurm <jid5> --log <out>/s5_report_<jid5>.log --hive`. When step 5
completes, point `run_de.R` at `<out>/report.parquet`.

## Parallelizing semi-tryptic / non-specific / InfinDIA searches (`--seed-lib`)

A **semi-tryptic** or **non-specific** DIA-NN search can't use the normal Step-1
predicted library — the predicted semi/non-specific library is enormous (a human
semi-tryptic `.speclib` was **19.6 GB**; searching it directly is impractical and
hogs the cluster). DIA-NN's answer is **InfinDIA** (`--pre-search`): an index-based
engine that builds a *small empirical* library from the data itself. It also
processes **DDA** (`--dda`, DIA-NN ≥2.3) — confirmed on the UC Davis nail cohort.

> **InfinDIA DDA is beta and can SEGFAULT on a single file.** On the 66-file nail
> cohort a single-shot `--pre-search --dda --semi` over all 66 crashed with
> `Segmentation fault (core dumped)` at file 20 — losing **all** work. Worse, the
> sbatch's trailing `echo "...exit $?"` masked the crash's non-zero code, so SLURM
> reported `State=COMPLETED exit 0:0` with **no `report.parquet`** (verify the output
> file, never trust State — and never end an sbatch with a bare `echo`; make the tool
> the last line or `rc=$?; …; exit $rc`). **This is the core reason to run InfinDIA
> DDA via the two-phase array below:** Phase 2 processes each file as an independent
> array task, so one bad file kills only its own task (the chain skips files with no
> `.quant`) instead of nuking the whole run.

**InfinDIA does NOT fan out per-file** like the 5-step chain. Its empirical library
is built from the whole experiment, so if you split `--pre-search` per file you'd get
N *different* libraries whose `.quant` files can't be merged. So parallelize in **two
phases**:

1. **Phase 1 — build the empirical library (one InfinDIA job).** Run
   `--pre-search --dda --semi --out-lib empirical.parquet` over a **representative
   subset** (DIA-NN's own tip: 20–100 high-quality runs) to keep it quick. InfinDIA
   **forces MBR + empirical-library generation** even without `--reanalyse` (it logs
   `enabling MBR and empirical spectral library generation, as required by InfinDIA
   pre-search`), so you always get a refined empirical library. Give it a fixed
   `--mass-acc`/`--mass-acc-ms1` (well-calibrated data) or a small `--ref` calibration
   library. For a search that's stuck **PENDING on `low`** (publicgrp is preemptible /
   congested), move it to `high`: `scontrol update jobid=<id> partition=high
   qos=genome-center-grp-high-qos account=genome-center-grp`.
2. **Phase 2 — fan the per-file passes out (`diann_parallel.py --seed-lib`).** Feed
   the empirical library from Phase 1 as the seed; Step 1 (prediction) is skipped and
   the small library drives the per-file first pass → assembly → final pass →
   cross-run, exactly the tested Steps 2–5. Fully parallel; this is where the speed is.
   **Do NOT put `--semi` in the Phase-2 `--cfg`** — the semi precursors are already in
   the empirical library, and re-adding `--semi` would re-expand the giant FASTA space
   you just avoided. Phase-2 cfg is a plain library search: `--dda` + `--mass-acc` +
   `--qvalue` + the var-mods that match the library.
   ```
   python3 scripts/diann_parallel.py \
     --raw-list all_runs.txt --fasta search.fasta --out ./p2 --diann '<binary>' \
     --cfg p2.cfg --seed-lib ./empirical.parquet \
     --seed-dep <phase1_jobid> \          # optional: first pass waits afterok on Phase 1
     --partition low --account publicgrp --qos publicgrp-low-qos \
     --threads-per-file 8 --max-simultaneous 66
   bash ./p2/submit.sh
   ```
   `--seed-dep` chains Phase 2 to start when Phase 1 finishes. **Verify the empirical
   library exists and is non-empty first** (State=COMPLETED ≠ a valid file) — a small
   launcher that globs for the `.parquet`, checks `>1 MB`, then generates+submits is the
   robust pattern (a watcher can trigger it on Phase-1 completion, no intervention).
- **`--qos`** is required to target `low`/publicgrp (`publicgrp-low-qos`); `high`/
  genome-center-grp uses its default qos. `--seed-lib` skips Step 1; combine with
  `--seed-dep` to auto-chain onto the InfinDIA lib-build job.

**Fairness note for cross-engine comparisons:** a subset-derived library, then used to
search all runs, does not bias detection of anything that recurs across samples (e.g.
keratin in hair/nail — the same peptides appear in every run). It can under-sample rare
precursors unique to un-subsetted runs. When that matters, build Phase 1 over all runs
(slower) or cross-check against a full single-shot InfinDIA run.

## Notes
- The generator emits **real absolute paths** (the HIVE DIA-NN 2.6 build is a native
  binary that reads `/quobyte` directly). If you instead use an Apptainer `.sif`, the
  `--diann` command must be the full `apptainer exec --bind … <sif> /diann-*/diann-linux`
  and the paths must be inside the bound mounts.
- Like the other engine paths, **validate the first real parallel run** end-to-end on
  HIVE before trusting it for production.

## Scan window MUST be pinned (not just mass accuracy)

Steps 3/5 reuse the `.quant` files written by steps 2/4. DIA-NN warns:

> WARNING: combining reuse of .quant files with automatic optimisation of mass
> accuracies **or scan window** will lead to results that are different from those of
> the original analysis that produced the .quant files and is strongly not recommended

`estimate_params.py` pins mass accuracy but omits `--window`, which means DIA-NN
optimises the radius **per file**. On a real 18-file poplar run that gave a radius of
**7 for seventeen files and 8 for one**, and the chain combined them anyway.

`mass_acc_status()` now treats a missing or `0` `--window` as NOT parallel-safe, so the
chain refuses to generate rather than producing a quietly-inconsistent report.

**The chain does this for you.** When the cfg has no `--window`, `diann_parallel.py`
inserts **step 1b** after library prediction: it runs `probe_window.py` on the first
file, writes the radius to `<out>/window.txt`, and steps 2–5 read that file at runtime,
so every pass uses the identical value. `--no-probe-window` disables it, in which case a
pinned `--window` in the cfg is required or the chain refuses to generate.

To measure it yourself instead, never by guessing — it depends on the acquisition scheme
(cycle time vs peak width), not the instrument model:

```bash
# after step 1 has produced the predicted library
python3 scripts/probe_window.py --diann "<diann cmd>" --raw <one file> \
    --fasta <f.fasta> --lib <step1.predicted.speclib> --write-cfg <params.cfg>
```

It terminates as soon as DIA-NN logs `Scan window radius set to N` (printed during
calibration, well before the search completes), so it costs minutes, not a full pass.
One file is enough for a set acquired with the same method; re-probe for a different
gradient, cycle time, or instrument.

## DIA-NN exits 0 on fatal errors — never trust the job state alone

Verified on DIA-NN 2.6.0: a run against a nonexistent `.mzML` and a nonexistent library
prints `ERROR: ...` and **returns exit code 0**.

DIA-NN is the last command in every generated step, so the step inherits that 0. SLURM
records `COMPLETED`, `watch_run.sh` reports success, and the `afterok` dependency
releases the next step. A step-4 array task that dies this way leaves no `.quant`, and
step 5 then builds the cross-run report from whatever survived — **a silently dropped
sample, reported as a clean run.**

Every step therefore asserts its own artefact:

| step | asserts |
|---|---|
| 1 | the predicted spectral library exists and is non-empty |
| 2 | this task's `.quant` was written to `quant_step2/` |
| 3 | the empirical spectral library exists |
| 4 | this task's `.quant` was written to `quant_step4/` |
| 5 | `report.parquet` exists **and** `quant_step4/` holds exactly one `.quant` per input |

Step 5's count check is the backstop: it is the only place that can notice a sample went
missing several steps earlier. On failure it prints the loop that names the missing file.
