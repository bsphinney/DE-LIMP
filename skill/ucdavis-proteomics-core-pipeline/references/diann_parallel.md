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

`diann_parallel.py` **refuses to run** on a cfg whose mass accuracy is omitted (auto) or
invalid rather than producing a silently inconsistent report. `--allow-auto-mass-acc`
overrides the *omitted* case only, for testing; it never overrides an invalid value (`0`,
which DIA-NN reads as a literal 0 ppm tolerance, negative, non-numeric, or set twice
differently), a bad `--window`, or an unparseable cfg. The refusal and `run_search.py`'s
routing decline come from the same `parallel_safe()` verdict and name the same fix.

**`--sbatch` does not apply to the chain.** It is six jobs chained by `submit.sh`, so
`run_search.py --sbatch job.sh` writes no `job.sh`, generates the chain, then moves an
existing regular file of that name aside to `job.sh.stale-<time>`, and exits **3** — so a
chained `&& sbatch job.sh` cannot resubmit a stale script. If `--sbatch` names a directory or
anything else that is not a regular file, it refuses up front and changes nothing; if
generation fails, the existing file is left where it was. Submit with `bash <out>/submit.sh`.

## The 5 steps (chained with `afterok` dependencies)
1. **Library prediction** — single job, no raw: predict a spectral library from the
   FASTA (`--fasta-search --predictor --gen-spec-lib --out-lib step1.speclib`).
2. **First pass** — SLURM **array**, one file per task: search each raw vs the
   predicted library (`--lib step1.predicted.speclib --temp quant_step2 --gen-spec-lib
   --quant-ori-names`) → per-file `.quant`.
3. **Empirical-library assembly** — single job: `--use-quant` over the step-2 `.quant`
   to build the empirical library (`--out-lib empirical.parquet`).
4. **Final pass** — SLURM array: search each raw vs the empirical library
   (`--lib empirical.parquet --temp quant_step4 --quant-ori-names`).
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
  the DIA-NN-recommended table in `estimate_params.py`). `--mass-acc 0` is **not** auto:
  DIA-NN fixes the tolerance at a literal 0 ppm and returns 0 IDs. Auto is the flags
  omitted, and neither can run as the chain.
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
- **Queue: detected, never hardcoded — and never demoted in silence.** Every sbatch the
  skill emits derives `--partition/--account/--qos` from `run_search.slurm_queue()`, which
  asks SLURM what this account is entitled to (`genome-center-grp/high` for facility
  members, `publicgrp/low` for everyone else) and then considers *utilisation*, not just
  entitlement — the per-user CPU cap on `high` means a big array can start sooner on `low`.
  If detection ever fails, the fallback is `publicgrp/low` **with a warning on stderr**: it
  is preemptible, and quietly putting a facility member's multi-hour search there is a real
  demotion. If you see that warning, pass `--partition high --account genome-center-grp`.
- **`high` needs no explicit `--qos`.** Measured on HIVE 2026-08-25: a job submitted with
  `--partition high --account genome-center-grp` and no `--qos` is accepted, and SLURM
  assigns `genome-center-grp-high-qos` by itself. `publicgrp/low` **does** need
  `--qos publicgrp-low-qos` named, and the generator adds it. Do not "fix" this by
  hardcoding the facility QOS — an account outside the group cannot submit with it.
- **The `--temp` folders MUST pre-exist.** DIA-NN aborts immediately with
  `ERROR: cannot find the temp folder .../quant_step2. Specify an existing folder` if
  the `--temp` dir is missing — it will **not** create it, and it aborts *before* doing any
  work, so a whole submission cycle is lost to a missing directory. **Every step now
  `mkdir -p`s its own** — the watcher playbook tells you to resubmit individual steps after a
  failure (`sbatch step4_finalpass.sbatch`), and that path never runs `submit.sh`.
  `submit.sh` also does
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

`parallel_safe()` is the single rule for "may this cfg run as the chain?", and both
`diann_parallel.py` (generate?) and `run_search.py` (auto-route?) gate on its `ok` — they
used to decide separately and drifted, which silently demoted a 310-file cohort to one
sequential search. It reads the cfg through one reader (`cfg_tokens`) that the step flags,
`params.base.cfg`, and `run_search.py`'s single-shot search and XIC check also use. It splits
words by **bash's** quoting rules — `#` starts a comment only at the start of a word, so
`/data/run#1` stays whole — because bash is what finally reads them. A `--window 0` mid-line
or after a tab, or a trailing `# comment`, therefore means the same thing to the gate as to
the generated steps. On the way out every value is quoted unless bash cannot change it:
`--var-mod "Phospho(STY),79.966331,STY"` used to be re-emitted bare and killed every step
with a bash syntax error after the gate had approved it, and `x{1,2}`, `~/libs` and globs
would be rewritten. Only `$NAME`/`${NAME}` still expand; `$(...)` and backticks are
literal. A cfg path that is not a file is reported as `cfg not found`, not as unpinned mass
accuracy. Unpinned or invalid **mass accuracy** is not
parallel-safe and refuses, rather than producing a quietly-inconsistent report.

`--window` must be a **positive integer**. DIA-NN does not accept `0` — it logs
`scan window radius should be a positive integer` and optimises per file (the poplar run
above) — so a missing or `0` `--window` is treated the same: recoverable, and step 1b
measures the radius (below). Anything else that is not a positive integer (`0.5`, `7.0`,
`-1`, `nan`, `wide`, or two different values) is a typo and refuses — the chain will not
measure over a mistake.

**The chain does this for you.** When the cfg has no usable `--window`, `diann_parallel.py`
inserts **step 1b** after library prediction: it hands `probe_window.py` the whole cohort
(`file_list.txt`), which measures the radius on **representative runs** and pins their
**median** (below — not the first file). It writes the radius to `<out>/window.txt` and
every probe to `<out>/window.json`, and steps 2–5 read `window.txt` at runtime, so every
pass uses the identical value. Any `--window` in the cfg is dropped from those
steps so the measured value cannot collide with it. For Thermo `.raw` inputs step 1b
exports the same .NET 8 environment every other step's DIA-NN gets — without it DIA-NN
cannot read `.raw`, no radius is logged, and steps 2–5 wait on `afterok` for ever.

Only after a radius is measured does step 1b write `<out>/params.resolved.cfg` — the cfg
plus the measured `--window` — so a "resolved" cfg without a window can never exist; a
resubmitted step 1b first removes the previous run's `window.txt`, `window.json`, resolved
cfg and probe logs. If no radius comes back, step 1b exits non-zero with `FAILED:` (keeping
`window.json` as the evidence) and steps 2–5 never
start: they were submitted `afterok` on that job id, so they sit `DependencyNeverSatisfied`
even if a fixed step 1b is resubmitted and succeeds. Cancel them and resubmit step 1b and
steps 2–5 chained on the new ids (ids in `jobs.txt`; reuse `step1.predicted.speclib`), or
re-run `submit.sh`, which also repeats step 1 (→ `references/watcher.md`).
`search_provenance.json` records the file as `resolved_params_file` with
`resolved_params_produced: "runtime"` (and `scan_window` saying it is measured by step 1b,
with `evidence_file` pointing at `window.json` and `probe_rule` stating the rule), because
at generation time it does not exist yet. Its job id is in `jobs.txt`
along with the rest of the chain, so `watch_run.sh --all` sees it fail.
`--no-probe-window` disables step 1b, in which case a pinned `--window` in the cfg is
required or the chain refuses to generate.

To measure it yourself instead, never by guessing — it depends on the acquisition scheme
(cycle time vs peak width), not the instrument model:

```bash
# after step 1 has produced the predicted library. Give it ALL the runs -- it chooses.
# Thermo .raw: export DOTNET_ROOT first, or DIA-NN cannot open the files.
export DOTNET_ROOT="$(bash scripts/ensure_dotnet8.sh | tail -1)"; export PATH="$DOTNET_ROOT:$PATH"
python3 scripts/probe_window.py --diann "<diann cmd>" --raw <all runs> \
    --fasta <f.fasta> --lib <step1.predicted.speclib> --write-cfg <params.cfg> \
    --workdir <scratch dir> [-- <the search's DIA-NN flags>] > window.json
```

It stops each DIA-NN as soon as it logs `Scan window radius set to N` (printed during
calibration, well before the search completes), so it costs minutes per run, not a full
pass. Re-probe for a different gradient, cycle time, or instrument.

### Which runs step 1b measures

Not the first file. DIA-NN's README ("Changing default settings") says automatic
optimisation "is inherently noisy: even replicate injections may not produce identical
results, and therefore the analysis results will depend on which run is first in the
list", and recommends running "on several representative runs". The first file of a
listing is whatever sorts first — often a blank, a QC injection, or a failed acquisition.
`probe_window.select_representative()`:

1. **Never probes a Bruker `.d` whose frame index cannot be trusted** — its `analysis.tdf`
   header is in WAL mode (bytes 18–19 = 2,2), a non-empty `analysis.tdf-wal` or `-journal`
   sits beside it, or the last indexed frame block (its `TimsId` offset plus the uint32
   block size stored there) ends before 99.9% of `analysis.tdf_bin` or past its end. The
   header and side files are read as bytes; the tdf is only ever opened with
   `file:<tdf>?mode=ro&immutable=1` (`mode=ro` alone still reads a stale WAL and leaves
   `-wal`/`-shm` files beside it). These runs are named as `WARNING` in the job log and
   listed under `excluded_damaged` in `window.json`: **they are still searched** by steps 2–5,
   so look at them.

   **One last resort.** A WAL-mode header is a sign of the damage, not the damage: a run can
   carry it with a complete index and nothing stale beside it, and a whole cohort can (6 of 6
   in review). When **nothing else in the cohort is probeable**, those runs — and only those:
   complete index, no `-wal`/`-journal`, nothing else wrong — are probed, listed under
   `probed_despite_wal_header` in `window.json` and named as a `WARNING` in the job log.
   Refusing them fails step 1b and leaves steps 2–5 in `DependencyNeverSatisfied` over runs
   those steps would have **searched** regardless. A `.d` whose `Frames` table has no readable
   `Time` (no such column, or NULL) is **not** damaged either — its spectra are all there; it
   just cannot be ranked by acquisition time, so the cohort falls back to size
   (`no_acquisition_time` in `window.json`). Naming runs with `--raw` overrides none of these
   checks: it narrows the cohort they are applied to.
2. **Never probes a run under half a typical run's size, measured per kind of run.** Typical is
   the median of the **larger half** of the runs *of that kind* — a blank, wash or failed
   acquisition falls under it. Not the median of all runs: with a blank after every sample the
   blanks are the majority and nothing would be excluded (6 samples + 7 blanks probed blank3,
   blank6 and s2). **Per kind**, because `size_bytes` is not one quantity: a `.d`'s is the
   indexed bytes of its `analysis.tdf_bin`, a `.raw`'s is the file. One floor across both
   excluded all three 0.9 GB Orbitrap `.raw` of a `.raw` + 8 GB `.d` cohort as blanks, flipped
   the ranking to time because only `.d` were left, and pinned a two-instrument cohort from one
   instrument. A cohort holding more than one kind gets `mixed_cohort` in `window.json` and a
   `WARNING` instead — one radius is not valid for two acquisition schemes; run the chain once
   per instrument. All directories share one kind, `.d` or not, because a failed Bruker
   acquisition is a `.d` with neither file in it and it is the floor against the real `.d`
   beside it that has to drop it. Then it **ranks what is left by how much was acquired**:
   when every run is a `.d` with a readable `Frames.Time`, by acquisition time, with the same
   half-of-typical floor on time (a short wash the size of a real run); otherwise by size (for
   a directory that is not a TDF run, every file under it — an Agilent `.d` keeps its data in
   `AcqData/`, and a top-level-only sum tied them all at ~0 bytes and sorted by name). The size
   floor comes first so a failed acquisition with no data at all (a `.d` with neither file)
   cannot switch a timsTOF cohort to size ranking.

   The floor is measured **from the cohort**, so a cohort that is *mostly* blanks or washes
   drags it down with them and excludes none of them (12 washes and 2 samples probes wash07,
   wash03 and wash10). Nothing here can tell a wash from a sample, so when the median run is
   under half the largest run left, `window.json` says so in `median_much_smaller` and the job
   log prints a `WARNING`: hand step 1b the sample runs, not the folder.
3. **Probes the median run first**, then the lower and upper quartile runs (positions ⌊m/4⌋
   and ⌈3m/4⌉ of the ranked eligible runs, distinct for every n ≥ 3), independent of input
   order.
4. **Replaces a run that logs no radius** with the remaining run nearest the median *in rank
   position* (`|i - n//2|` in the ranking, not the nearest value — they differ wherever the
   ranked values are spaced unevenly) — the
   fall-through from PR #70, applied to representative runs rather than to the next file.
   After 3 runs without a radius step 1b gives up. DIA-NN's missing-.NET error is not
   replaced: no other `.raw` can succeed in that environment, so it stops at once and prints
   the export.
5. **Pins the median** of the measured radii (the lower middle of two). If fewer than three
   answered — the replacements ran out, or the time budget did — it still pins, as the
   first-file fall-through did, but `window.json` says `incomplete: true` and the job log
   prints a `WARNING`. `radii_agree` is `true`/`false` only when **more than one** radius was
   measured — `null` otherwise, because one measurement corroborates nothing (a cohort of two
   runs is probed once) — and radii that differ by more than 2 get their own `WARNING` in the
   job log rather than only `radii_agree: false` in the JSON. `window.json` records every probe
   (role, size, time, radius, seconds, log) and every run left out and why.

**Thermo `.raw`** is ranked by size only. Reading a `.raw`'s acquisition length or checking
its integrity needs Thermo's RawFileReader (.NET), the library DIA-NN itself loads, so the
probe does not; a `.raw` DIA-NN cannot open logs no radius and is replaced.

Why these rules, measured on HIVE:

- **The truncated index.** 342 of 39,374 Bruker `.d` on HIVE have an `analysis.tdf` whose
  frame index stops early while `analysis.tdf_bin` is complete (a stale mid-acquisition
  `-wal` was checkpointed into the finished tdf by a read-write open); every one has a WAL
  header, intact ones are 1,1. The FRAN pilot's size rule picked one:
  `8aug25_KogantiLys7uL_60spd_6_S4-B6_1_15704.d`, `tdf_bin` 2.40 GB, index 1,451 frames /
  134 s covering **0.74%** of it (read with the immutable URI), where its 11 siblings index
  13,736–13,743 frames / 21.0 min covering 99.9999–100%. DIA-NN read 121 cycles and 80
  precursors. On that cohort (srun jobs 23536639 and 23536675) the rule excludes run 6 by its WAL header,
  ranks the other 11 by time, and picks B11 / B3 / B9 in 0.07 s; two intact blanks with a
  4.2 MB stale `-wal` (`apr26/blankDia_S1-H11_1_21474.d`, `June26/blankDia_S1-H2_1_22321.d`)
  are refused without being opened, and no file in any of the 14 `.d` changed.
- **Half, and of the larger half** — sizes in `raw_data/Lumos1/noi25`: washes 120–240 MB
  and 60-min washes 527–567 MB beside 80 90-min DIA runs of 0.93–2.85 GB (median 2.12 GB;
  median of the larger half 2.42 GB). The 60-min washes are 25–27% of that DIA median, so a
  25% floor let them in. At 50% of the larger half's median they are out; 2 of the 80 DIA
  runs (38–40%) are not probed (still searched). The rule keeps blanks out while they are at
  most 60–75% of the list; the 16 60-min runs of one project plus every wash in that folder
  (86% washes) probe three washes, so hand the search the cohort, not the folder.
- **Several runs, median** — DIA-NN 2.7.0, 200-protein library / full predicted library:

  | set | runs probed | full library | 200-protein library |
  |---|---|---|---|
  | Exploris 480 `.raw`, 5 runs | TT34 / TT33 / TT32 | 7, 7, 7 → 7 | 7, 8, 8 → 8 |
  | timsTOF Pro `.d`, 8 runs + a 0-byte failed acquisition (excluded) | A3 / A7 / A4 | 14, 10, 11 → 11 | 15, 12, 12 → 12 |

  A run's radius is reproducible (A3 gave 14 in three separate probes); it is the runs that
  differ, so a one-file probe of that timsTOF cohort pins 10, 11 or 14 depending on the run.
  The old first-file step 1b handed DIA-NN the failed acquisition, never got a radius, and
  DIA-NN wrote `report.parquet` into `<out>` — the path step 5's check reads.
- **This step 1b, run as a compute node would** (DIA-NN 2.7.0, 200-protein library, 16 CPUs):
  - timsTOF, the 8 runs + failed acquisition above (srun job 23538112): the failed acquisition
    is excluded by size, the 8 runs are ranked by acquisition time, A6 / A2 / A1 give 12, 12, 11
    → **12** in 252 s, 129 DIA-NN lines streamed into the job log, no `diann-linux` left. Every
    run there is 21.0 min, so which full-length runs the quartiles land on comes down to
    sub-second differences in the last frame's time; ranked by size instead (srun job
    23537237, before the size floor came first) it probed A7 / A3 / A4 → 12, 15, 12 → **12**.
  - Exploris `.raw`, 5 runs (srun job 23537237): TT33 / TT34 / TT32 → 8, 7, 8 → **8** in 67 s.
    With the `.NET` export line removed: exit 1 after one probe in 1.8 s (`stopped_because:
    environment`), no `window.txt`, no resolved cfg, no traceback.

**Each probe gets its own `--temp` and `--out`** under `<out>/window_probe/probe<N>_<run>/`
(with its `probe.log`), placed before `--threads` so the search's flags stay the tail of
DIA-NN's argv, as in steps 2–5. Without them a probe that never logs a radius writes its
report into `<out>` and, if the search completes, its `.quant` next to the raw file.

**One at a time, live in the job log.** Measured on a 16-CPU `high` allocation with full
mouse predicted libraries, the same three runs each way: timsTOF `.d` 1083 s one at a time
vs 1314 s concurrent (loading each `.d` took 227–238 s concurrent against 39–43 s alone);
Exploris `.raw` 289 s vs 250 s. Neither wins, concurrent needs three probes' memory (~30 GB
on those timsTOF runs), and one at a time is what lets a failed run be replaced. Each probe
is announced (`[probe_window] probe 1/3: <run> (median, 2.20 GB, 21.0 min), 16 threads`) and
DIA-NN's output is copied into the step-1b log as it arrives, tagged `[probe N]`: 1083 s of
silence is past `watch_run.sh`'s 15-minute stall rule, whose playbook is to cancel the job.

**Time.** `--timeout` (3600 s) is per probe, so one hung run is cut and replaced — one
timsTOF probe once made no progress for 28 minutes while other jobs read the same storage.
`--budget` covers all probes, replacements included: step 1b passes its 4-hour wall clock
less 10 minutes (13,800 s), so the probe stops itself and writes `window.json` before SLURM
would kill the job. Three probes at the full 3600 s still fit.

**Stopping DIA-NN means all of it.** `--diann` is not always the engine: without a native
build `acquire_tools.sh` records `apptainer exec --bind /quobyte:/quobyte <sif>
/diann-*/diann-linux`, and sites wrap binaries in scripts. DIA-NN runs in its own process
group, writes to `probe.log` rather than a pipe a grandchild could hold open, and every stop
(radius found: SIGTERM then SIGKILL after 30 s; timeout, budget, or SIGTERM/SIGINT/SIGHUP to
the probe: SIGKILL) goes to the whole group. On HIVE (apptainer 1.5.3, srun job 23521026) with
`--diann "apptainer exec … bash -c 'sleep 40; echo Scan window radius set to 9'"`: a 5 s
timeout returned after 5.3 s with 0 container processes left, and SIGTERM to the probe left
0 (a probe that signalled only its child returned after 41.2 s and left 2).

## DIA-NN exits 0 on fatal errors — never trust the job state alone

Verified on DIA-NN 2.6.0: a run against a nonexistent `.mzML` and a nonexistent library
prints `ERROR: ...` and **returns exit code 0**.

DIA-NN is the last command in every generated step, so the step inherits that 0. SLURM
records `COMPLETED`, `watch_run.sh` reports success, and the `afterok` dependency
releases the next step. A step-4 array task that dies this way leaves no `.quant`, and
step 5 then builds the cross-run report from whatever survived — **a silently dropped
sample, reported as a clean run.**

Every step therefore asserts its own artefact — and first **deletes** it, because an
existence check cannot tell this run's file from the previous run's. Re-run a chain into
the same `--out`, or resubmit a step as `references/watcher.md` says to, and a DIA-NN that
exits 0 having written nothing leaves the old file for the check to find. Reproduced on
DIA-NN 2.7.0 (HIVE, review srun 23512013): a re-run with no .NET logged `ERROR: cannot read
.raw files`, exited 0, and left the previous `report.parquet` byte-identical.

| step | deletes first | asserts |
|---|---|---|
| 1 | `step1.predicted.speclib` | the predicted spectral library exists and is non-empty |
| 1b (window measured) | `window.txt`, `window.json`, `params.resolved.cfg` and its `.tmp` | a positive-integer radius was measured, and `window.txt` and `params.resolved.cfg` are non-empty regular files |
| 2 | this task's `quant_step2/<run>.quant` | this task's `.quant` was written to `quant_step2/` |
| 3 | `empirical.parquet` | the empirical spectral library exists |
| 4 | this task's `quant_step4/<run>.quant` (before the no-step-2-quant skip) | this task's `.quant` was written to `quant_step4/` |
| 5 | `report.parquet` and `report.stats.tsv` | `report.parquet` exists **and** every `.quant` the runs in `file_list.txt` map to is present and non-empty in `quant_step4/` |

Step 5's count check is the backstop: it is the only place that can notice a sample went
missing several steps earlier. It counts the `.quant` **files** this chain's inputs map to —
not every `.quant` in the folder, because a previous search's `.quant` for some other run
would make up the number, and not one per input **line**, because DIA-NN names a run (and
this generator names its `.quant`) by the file name alone: `/plate1/s1.d` and `/plate2/s1.d`
share one `s1.quant`, and counting per line counts that single surviving file twice and
passes a report in which the two samples were merged. The names are deduplicated (`sort -u`),
so a collision makes the count fall short and the job fails. It names each missing file
(`MISSING: <run>.quant -- no step-4 task wrote it`), and reports a short count of distinct
names separately (`N inputs map to only M distinct run names`).

Colliding run names are refused before any of that, at generation: `run_search.py` checks the
input list it routes, and `diann_parallel.py` checks `--raw`/`--raw-list` in its own `main()`
— the chain is run directly too, and that route used to generate a chain, spend its SLURM
hours and merge two samples into one column.

`--out` is refused if it contains `$`, a backtick, `"`, a backslash or a newline. Every guard
above quotes its path with **double** quotes on purpose, so that an array task's `$QUANT`
expands; a `$(...)` in `--out` would therefore be command substitution that runs when the
**job** runs — inside the `rm -f --` the step performs before DIA-NN starts, so it would also
choose what gets deleted.

The **single-shot** search `run_search.py` writes for ≤ 5 files has the same contract: its
library job deletes, then asserts, `diann_lib.predicted.speclib`, and its search job deletes
`report.parquet` and `report.stats.tsv`, then asserts the report and runs
`scripts/check_report_runs.py`, which fails unless the report holds one `Run` per input
(→ `references/search-engines.md`). Re-measured on DIA-NN 2.7.0:
a search over one readable `.raw`, one truncated `.raw` and one nonexistent path exited 0
with a report holding 1 of the 3 runs.

The input list that check reads is written next to the job and named after it
(`<sbatch stem>_input_files.txt`), so a second search generated into the same `--out` cannot
rewrite the list a first, still-unsubmitted job would read. Running the search **inline**
(`--allow-inline`, no `--sbatch`) renames the previous `report.parquet`, `report.stats.tsv`
and predicted library to `<name>.stale-<timestamp>` rather than deleting them — the job route
can only delete, because `clear_stale` runs on a compute node long after generation, but
inline we are the ones re-running and a re-run whose DIA-NN dies should still leave the user
the report they had. Inline also checks the report with `getsize > 0`, matching the job
route's `[ -s ]`: DIA-NN creating `report.parquet` and then dying leaves it at 0 bytes.
