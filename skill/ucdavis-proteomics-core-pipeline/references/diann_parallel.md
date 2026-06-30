# DIA-NN parallel search (5-step SLURM chain)

DIA-NN's high-throughput workflow — **poorly documented upstream**, so it is encoded
here (`diann_parallel.py`), ported faithfully from DE-LIMP's `generate_parallel_scripts()`.
Use it for **many DIA files on a cluster** (roughly ≥ ~8–10, or whenever a single
job would be too slow). For a few files, the single-shot `run_search.py` is fine.

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
`watch_run.sh --slurm <jid5> --log <out>/s5_report_<jid5>.log --hive`. When step 5
completes, point `run_de.R` at `<out>/report.parquet`.

## Notes
- The generator emits **real absolute paths** (the HIVE DIA-NN 2.6 build is a native
  binary that reads `/quobyte` directly). If you instead use an Apptainer `.sif`, the
  `--diann` command must be the full `apptainer exec --bind … <sif> /diann-*/diann-linux`
  and the paths must be inside the bound mounts.
- Like the other engine paths, **validate the first real parallel run** end-to-end on
  HIVE before trusting it for production.
