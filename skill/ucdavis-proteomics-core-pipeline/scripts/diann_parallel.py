#!/usr/bin/env python3
"""
diann_parallel.py  --  Generate DIA-NN's canonical 5-step PARALLEL search as a SLURM
job chain. This is the high-throughput DIA-NN workflow that is poorly documented
upstream; ported faithfully from DE-LIMP's generate_parallel_scripts() (R/helpers_
search.R) and the facility's usage.

The 5 steps (each chained `afterok` on the previous):
  1  library prediction  single job, no raw — predict a spectral library from the FASTA
  2  first pass          SLURM array (1 file/task) — search vs the predicted lib -> .quant
  3  empirical assembly  single job — `--use-quant` over step-2 .quant -> empirical lib
  4  final pass          SLURM array — search vs the empirical lib -> .quant
  5  cross-run report    single job — `--use-quant --matrices` -> report.parquet

Why it's faster: the per-file passes (steps 2 & 4) run as a SLURM **array** across many
nodes at once instead of one long single-node job; MBR is replaced by the empirical-
library round-trip. **Mass accuracy is FIXED (manual), not auto** — steps 3/5 reuse the
.quant files and auto-calibration would be inconsistent (per DIA-NN dev guidance).

Writes into <out>: `file_list.txt`, `step{1..5}_*.sbatch`, and `submit.sh` (submits the
chain with dependencies). Run `submit.sh` on the cluster (or via `hive_exec.sh`). All
heavy work runs on compute nodes through the array — never the login node.

Usage:
  python3 diann_parallel.py --diann '<diann binary | apptainer exec ... diann-linux>' \
      --raw /data/*.d --fasta /path/search.fasta --out ./diann_parallel \
      --cfg params.cfg [--threads-per-file 16] [--mem-per-file 64] [--time-per-file 2] \
      [--assembly-cpus 64] [--assembly-mem 128] [--assembly-time 12] \
      [--partition high] [--account genome-center-grp] [--max-simultaneous 20] [--no-norm]
"""
import os, sys, glob, argparse

# flags that are step-specific or auto-determined — never carry them into every step
STRIP = ("--fasta-search", "--predictor", "--gen-spec-lib", "--matrices", "--reanalyse",
         "--rt-profiling", "--no-norm", "--xic", "--out-lib", "--lib", "--out", "--f",
         "--fasta", "--threads", "--temp")


def read_cfg_flags(cfg):
    """Read a diann.cfg into a flat flag string, dropping step-specific flags."""
    if not cfg or not os.path.exists(cfg):
        return ""
    out = []
    for raw in open(cfg):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if any(line == s or line.startswith(s + " ") for s in STRIP):
            continue
        out.append(line)
    return " ".join(out)


def header(name, cpus, mem_gb, hours, partition, account, array=None):
    h = ["#!/bin/bash -l",
         f"#SBATCH --job-name={name}",
         f"#SBATCH --cpus-per-task={cpus}",
         f"#SBATCH --mem={mem_gb}G",
         f"#SBATCH --time={hours}:00:00",
         f"#SBATCH --partition={partition}",
         f"#SBATCH --account={account}",
         f"#SBATCH -o {name}_%j.log", f"#SBATCH -e {name}_%j.log"]
    if array:
        h.insert(2, f"#SBATCH --array={array}")
        h = [x.replace("_%j.log", "_%A_%a.log") for x in h]
    return "\n".join(h)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diann", required=True, help="DIA-NN command (native binary path, or 'apptainer exec --bind … <sif> /diann-*/diann-linux')")
    ap.add_argument("--raw", nargs="+", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cfg", help="diann.cfg with the search params (estimate_params.py output)")
    ap.add_argument("--threads-per-file", type=int, default=16)
    ap.add_argument("--mem-per-file", type=int, default=64)
    ap.add_argument("--time-per-file", type=int, default=2)
    ap.add_argument("--assembly-cpus", type=int, default=64)
    ap.add_argument("--assembly-mem", type=int, default=128)
    ap.add_argument("--assembly-time", type=int, default=12)
    ap.add_argument("--libpred-cpus", type=int, default=16)
    ap.add_argument("--libpred-mem", type=int, default=64)
    ap.add_argument("--libpred-time", type=int, default=4)
    ap.add_argument("--partition", default="high")
    ap.add_argument("--account", default="genome-center-grp")
    ap.add_argument("--max-simultaneous", type=int, default=20)
    ap.add_argument("--no-norm", action="store_true")
    a = ap.parse_args()

    raws = []
    for p in a.raw:
        raws.extend(sorted(glob.glob(p)) or [p])
    raws = [os.path.abspath(r.rstrip("/")) for r in raws]
    n = len(raws)
    if n < 2:
        sys.exit("Parallel search needs >= 2 raw files; use the single-shot run_search.py for 1.")
    out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    fasta = os.path.abspath(a.fasta)
    DN = a.diann
    flags = read_cfg_flags(a.cfg)
    D = out  # all DIA-NN intermediate/output lives here (real paths; native binary reads them directly)
    report = "no_norm_report.parquet" if a.no_norm else "report.parquet"
    norm = "--no-norm" if a.no_norm else ""

    # file list (1 raw path per line) — array tasks index into it
    open(os.path.join(out, "file_list.txt"), "w").write("\n".join(raws) + "\n")
    all_f = " ".join(f"--f {r}" for r in raws)
    array = f"0-{n-1}%{a.max_simultaneous}"
    predicted = f"{D}/step1.predicted.speclib"
    empirical = f"{D}/empirical.parquet"

    def write(name, body):
        p = os.path.join(out, name)
        open(p, "w").write(body + "\n"); os.chmod(p, 0o755); return name

    # array preamble: pick this task's raw file
    pick = ('FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ' + f'{D}/file_list.txt)\n'
            'if [ -z "$FILE" ]; then echo "no file for task $SLURM_ARRAY_TASK_ID"; exit 1; fi\n'
            'echo "Processing: $FILE"\n')

    # Step 1 — library prediction (single job)
    s1 = write("step1_libpred.sbatch", "\n".join([
        header("s1_libpred", a.libpred_cpus, a.libpred_mem, a.libpred_time, a.partition, a.account), "",
        f'echo "Step 1/5 library prediction"; date',
        f'{DN} --fasta {fasta} --fasta-search --predictor --gen-spec-lib \\',
        f'  --out-lib {D}/step1.speclib --out {D}/step1_lib.parquet \\',
        f'  --threads {a.libpred_cpus} {flags}']))

    # Step 2 — first pass (array): predicted lib -> per-file .quant
    s2 = write("step2_firstpass.sbatch", "\n".join([
        header("s2_firstpass", a.threads_per_file, a.mem_per_file, a.time_per_file, a.partition, a.account, array=array), "",
        f'echo "Step 2/5 first pass, task ${{SLURM_ARRAY_TASK_ID}} of {n}"; date', pick,
        f'{DN} --f "$FILE" --fasta {fasta} --lib {predicted} \\',
        f'  --temp {D}/quant_step2 --rt-profiling --gen-spec-lib --quant-ori-names \\',
        f'  --threads {a.threads_per_file} {flags}']))

    # Step 3 — empirical library assembly (single job, --use-quant)
    s3 = write("step3_assembly.sbatch", "\n".join([
        header("s3_assembly", a.assembly_cpus, a.assembly_mem, a.assembly_time, a.partition, a.account), "",
        f'echo "Step 3/5 empirical library assembly"; date',
        f'cp -r {D}/quant_step2 {D}/quant_step2_orig 2>/dev/null || true   # backup for resume',
        f'{DN} {all_f} --fasta {fasta} --lib {predicted} --use-quant --quant-ori-names \\',
        f'  --rt-profiling --gen-spec-lib --out-lib {empirical} \\',
        f'  --temp {D}/quant_step2 --out {D}/step3_assembly.parquet \\',
        f'  --threads {a.assembly_cpus} {flags}']))

    # Step 4 — final pass (array): empirical lib -> per-file .quant
    s4 = write("step4_finalpass.sbatch", "\n".join([
        header("s4_finalpass", a.threads_per_file, a.mem_per_file, a.time_per_file, a.partition, a.account, array=array), "",
        f'echo "Step 4/5 final pass, task ${{SLURM_ARRAY_TASK_ID}} of {n}"; date', pick,
        'QUANT="${FILE##*/}"; QUANT="${QUANT%.*}.quant"',
        f'if [ ! -f "{D}/quant_step2/$QUANT" ]; then echo "SKIP: no step-2 quant for $QUANT"; exit 0; fi',
        f'{DN} --f "$FILE" --fasta {fasta} --lib {empirical} \\',
        f'  --temp {D}/quant_step4 --no-ifs-removal --quant-ori-names \\',
        f'  --threads {a.threads_per_file} {flags}']))

    # Step 5 — cross-run report (single job, --use-quant --matrices)
    s5 = write("step5_report.sbatch", "\n".join([
        header("s5_report", a.assembly_cpus, a.assembly_mem, a.assembly_time, a.partition, a.account), "",
        f'echo "Step 5/5 cross-run report"; date',
        f'{DN} {all_f} --fasta {fasta} --lib {empirical} --use-quant --quant-ori-names \\',
        f'  --temp {D}/quant_step4 --matrices --out {D}/{report} \\',
        f'  --threads {a.assembly_cpus} {norm} {flags}']))

    # submit.sh — chain the 5 steps with afterok dependencies
    submit = "\n".join([
        "#!/bin/bash", "set -euo pipefail", f'cd "{out}"',
        'jid1=$(sbatch --parsable %s)' % s1,
        'jid2=$(sbatch --parsable --dependency=afterok:$jid1 %s)' % s2,
        'jid3=$(sbatch --parsable --dependency=afterok:$jid2 %s)' % s3,
        'jid4=$(sbatch --parsable --dependency=afterok:$jid3 %s)' % s4,
        'jid5=$(sbatch --parsable --dependency=afterok:$jid4 %s)' % s5,
        'echo "submitted: lib=$jid1 firstpass=$jid2 assembly=$jid3 finalpass=$jid4 report=$jid5"',
        f'echo "final report will be {D}/{report}; watch with: watch_run.sh --slurm $jid5 --log {D}/s5_report_${{jid5}}.log"'])
    write("submit.sh", submit)

    import json
    print(json.dumps({
        "out": out, "n_files": n, "report": f"{D}/{report}",
        "scripts": [s1, s2, s3, s4, s5, "submit.sh"],
        "submit": f"bash {out}/submit.sh   (or: hive_exec.sh 'bash {out}/submit.sh')",
        "report_jobid_var": "jid5",
        "note": "5-step DIA-NN parallel chain. Mass accuracy is fixed (manual) — ensure --cfg has --mass-acc/--mass-acc-ms1 set, not 0/auto. Watch the run with watch_run.sh. Then point run_de.R at the report.",
    }, indent=2))


if __name__ == "__main__":
    main()
