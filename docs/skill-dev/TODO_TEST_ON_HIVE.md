# First-run validation checklist (HIVE)

These parts of the skill were built/verified **without** a real run on HIVE (no GPU /
cluster from the dev box). Each was grounded in an authoritative source (DE-LIMP code,
the facility's HIVE files, DIA-NN/AlphaDIA docs) and unit-tested where possible, but
the first real run should confirm the items below. Tick them off, and fix + note any
surprise in the relevant `references/*.md`.

Run everything via `hive_exec.sh` (Claude Code local, HIVE over SSH) and submit heavy
work with `sbatch` — never the login node.

---

## 1. DIA-NN 2.6 single-shot on HIVE
- [ ] `acquire_tools.sh hpc` (with `PIN_ENGINE=diann PIN_VERSION=2.6.0`) sets
      `diann` to `/quobyte/proteomics-grp/dia-nn/build_260/diann-2.6.0/diann-linux`
      (NOT the old 2.3.0 `.sif`).
- [ ] The native binary runs on a compute node (deps OK) — compare invocation to the
      facility's `/quobyte/proteomics-grp/dia-nn/run_diann_*.sbatch`.
- [ ] A small search produces `report.parquet` that `run_de.R` reads.

## 2. DIA-NN 5-step parallel chain (the big one)
`diann_parallel.py` is ported from DE-LIMP's 2.3.0-era logic; confirm it still holds for
the **2.6 native binary**. Generate the chain, run `submit.sh`, and check each step:
- [ ] **Step 1** writes `step1.predicted.speclib` (the exact name steps 2/3 expect — a
      DIA-NN naming quirk; if 2.6 names it differently, fix the `predicted`/`empirical`
      paths in `diann_parallel.py`).
- [ ] **Step 2** array → per-file `<basename>.quant` in `quant_step2/` (`--quant-ori-names`).
- [ ] **Step 3** `--use-quant` → `empirical.parquet`.
- [ ] **Step 4** array → `quant_step4/` (skips files missing a step-2 `.quant`).
- [ ] **Step 5** `--use-quant --matrices` → `report.parquet`.
- [ ] `afterok` dependencies chain correctly; the array honors `%max_simultaneous`.
- [ ] Protein/precursor counts are in line with a single-shot run on the same files.

## 3. AlphaDIA (commercial-OK engine)
- [ ] On HIVE: `acquire_tools.sh` finds `/quobyte/proteomics-grp/apptainers/alphadia.sif`;
      off-HIVE: `pip install alphadia` works.
- [ ] A small DIA run produces `pg.matrix.parquet`; **confirm `adapt_alphadia` column
      assumptions** — id column (`pg`/`pg.name`), sample columns = run names, values =
      intensity. Fix the column detection if the real headers differ.
- [ ] Confirm the limitation: AlphaDIA does **not** do timsTOF whole-proteome directDIA
      (it should be routed away from that case).
- [ ] AlphaDIA needs a **GPU** — submit to a GPU node (`--gres=gpu:1`).

## 4. Sage / FragPipe adapters (already flagged)
- [ ] Sage `lfq.parquet` → `adapt_sage` → `report.parquet` columns correct on real data.
- [ ] `run_search.py --engine sage --sbatch` on `.raw` (2.9): no converter runs at generation;
      the job converts with the Core's ThermoRawFileParser (`-f=2`, `DOTNET_ROOT` baked in),
      every mzML ends in `</indexedmzML>`, and Sage reads them. Only fixtures and stand-ins
      have tested this so far.
- [ ] `sage_lfq_check.py` on gabrig's HeL50 Fusion Lumos Sage output (2026-09-29): the +7.1/+7.7
      ppm runs WARN (Sage logged 0 target MS1 peaks); of the +1.8 to +5.8 ppm runs, those whose
      median + 2 ppm exceeds the 5 ppm window are named. The median `precursor_ppm` it prints should match what
      she measured from `results.sage.parquet`.
- [ ] FragPipe `combined_protein.tsv` → `adapt_fragpipe` correct (if used).

## 5. HIVE remote execution + access gate
- [ ] `check_access.sh <user> <key>` correctly reports `hive_remote` + `facility_software`.
- [ ] `hive_exec.sh` runs commands + `--put`/`--get` staging with the user's key.
- [ ] **Rebuild-on-HIVE** path (non-Core): `setup.sh` + `acquire_tools.sh` +
      `fetch_fasta.py` build a working env in the user's home.
- [ ] `fetch_fasta.py fetch --hive` picks `/quobyte/proteomics-grp/MRS/UP000005640_9606.fasta`
      (20,663 seqs) and **not** `UP000005640_9606_plus_universal_contam.fasta` — the
      latter already contains contaminants and would double-append them. Verified
      against a local mock of the MRS layout; confirm on the real mount.
- [ ] `--contaminants universal --hive` resolves to MRS
      `kg_nov/fasta/Universal Protein Contaminants.fasta` (381 seqs, `Cont_`-tagged),
      and a set NOT staged there (e.g. `rat_tissue`) falls through to the Hao lab
      GitHub download rather than silently substituting the wrong set.

## 6. The watcher
- [ ] `watch_run.sh --slurm <jid> --hive` returns the right `state`/`done`/`failed`.
- [ ] A deliberately under-resourced job (`--mem` too low) is caught as `out_of_memory`
      with the correct fix; resubmit-after-fix loop works.

## 7. End-to-end
- [ ] For each engine, the adapted `report.parquet` → `run_de.R` → figures → report →
      reproducibility bundle, with the search command + version recorded in
      `search_provenance.json` and the bundle.

---
When an item passes, note any deviation (paths, flag names, column headers) and update
the matching `references/*.md` so the next person/agent doesn't re-discover it.
