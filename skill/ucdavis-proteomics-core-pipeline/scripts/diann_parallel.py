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
      [--partition <auto>] [--account <auto>] [--max-simultaneous 20] [--no-norm]
"""
import os, sys, glob, argparse, shlex, subprocess

# flags that are step-specific or auto-determined — never carry them into every step.
# NOTE: --dda is intentionally NOT stripped — for DDA data put --dda in the --cfg and it
# flows into every step (DIA-NN 2.6 searches DDA per file exactly as it does DIA).
STRIP = ("--fasta-search", "--predictor", "--gen-spec-lib", "--matrices", "--reanalyse",
         "--rt-profiling", "--no-norm", "--xic", "--out-lib", "--lib", "--out", "--f",
         # NOTE: --xic is stripped here on purpose and re-added to step 4 ONLY (see
         # xic_flag() below) -- step 2 IDs are not final, and step 5 runs --use-quant,
         # which never re-reads the raw spectra so --xic is silently a no-op there.
         "--fasta", "--threads", "--temp")


def dotnet_prefix(raws):
    """If any input is Thermo .raw, the DIA-NN 2.6 native binary needs a .NET 8 runtime
    (>= 8.0.17) on PATH to read it. Resolve/install via ensure_dotnet8.sh and return an
    'export DOTNET_ROOT=...; export PATH=...; ' prefix to put in front of every DIA-NN
    invocation (each array task/step runs it; the shared install is read on the node).
    Returns "" for mzML/.d-only inputs. Run the generator on a login node (internet)."""
    if not any(r.lower().endswith(".raw") for r in raws):
        return ""
    helper = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ensure_dotnet8.sh")
    try:
        root = subprocess.check_output(["bash", helper], text=True).strip().splitlines()[-1]
    except Exception as e:
        sys.stderr.write(f"[diann_parallel] ensure_dotnet8.sh failed ({e}); DIA-NN may not "
                         "read .raw. Provide mzML or install .NET 8 >= 8.0.17.\n")
        return ""
    return f'export DOTNET_ROOT={root}; export PATH={root}:"$PATH"; '


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


def xic_flag(cfg):
    """Return the '--xic N' flag from the cfg, or '' if XICs were not requested.

    XIC extraction must happen on the FINAL per-file pass (step 4): step 2 works
    against the predicted library so its IDs are not final, and step 5 runs with
    --use-quant, which reuses .quant files without re-reading the raw spectra --
    DIA-NN accepts --xic there and logs that it will extract, but writes nothing.
    """
    if not cfg or not os.path.exists(cfg):
        return ""
    try:
        toks = shlex.split(open(cfg).read(), comments=True)
    except ValueError:
        toks = open(cfg).read().split()
    for i, t in enumerate(toks):
        if t == "--xic":
            if i + 1 < len(toks) and not toks[i + 1].startswith("-"):
                return f"--xic {toks[i + 1]}"
            return "--xic"
    return ""

def mass_acc_status(cfg):
    """Is this cfg PARALLEL-SAFE? Steps 3/5 reuse the .quant files from steps 2/4, so
    anything DIA-NN auto-optimises per file is applied inconsistently between passes
    and then stitched together. DIA-NN says so itself:

        WARNING: combining reuse of .quant files with automatic optimisation of mass
        accuracies OR SCAN WINDOW will lead to results that are different from those
        of the original analysis that produced the .quant files and is strongly not
        recommended

    So this checks BOTH mass accuracy and --window. Checking only mass accuracy was a
    real defect: on an 18-file poplar run with mass accuracy correctly pinned, DIA-NN
    still inferred a scan-window radius of 7 for seventeen files and 8 for one, emitted
    the warning above, and the chain combined them anyway.

    This is the single definition of "is this cfg parallel-safe"; run_search.py imports
    it to decide whether it may auto-route. Returns {fixed, ms2, ms1, window, reason}."""
    vals = {"--mass-acc": None, "--mass-acc-ms1": None, "--window": None}
    if cfg and os.path.exists(cfg):
        try:
            toks = shlex.split(open(cfg).read(), comments=True)
        except ValueError:                      # unbalanced quotes -- fall back to raw split
            toks = open(cfg).read().split()
        for i, t in enumerate(toks):
            if t in vals and i + 1 < len(toks):
                try:
                    vals[t] = float(toks[i + 1])
                except ValueError:
                    pass
    ms2, ms1, win = vals["--mass-acc"], vals["--mass-acc-ms1"], vals["--window"]
    pairs = (("--mass-acc", ms2), ("--mass-acc-ms1", ms1), ("--window", win))
    missing = [k for k, v in pairs if v is None]
    auto = [k for k, v in pairs if v == 0]      # 0 means "optimise automatically"
    if missing or auto:
        why = []
        if missing:
            why.append("not set: " + ", ".join(missing))
        if auto:
            why.append("auto (0): " + ", ".join(auto))
        return {"fixed": False, "ms2": ms2, "ms1": ms1, "window": win,
                "reason": "; ".join(why)}
    return {"fixed": True, "ms2": ms2, "ms1": ms1, "window": win,
            "reason": f"fixed (MS1 {ms1} ppm / MS2 {ms2} ppm / window {int(win)})"}


# DIA-NN 2.6 EXITS 0 ON A FATAL ERROR -- verified: a run against a nonexistent .mzML and
# a nonexistent library prints "ERROR: ..." and returns exit code 0. Since DIA-NN is the
# last command in every step, SLURM records COMPLETED, watch_run.sh reports success, and
# the afterok dependency releases the next step. A step-4 task that dies this way simply
# leaves no .quant, and step 5 then builds the cross-run report from the survivors and
# calls it done -- a silently DROPPED SAMPLE. Exit status is therefore not trustworthy;
# every step must assert that the artefact it was supposed to produce actually exists.
def must_exist(path, what):
    """Bash that fails the job unless `path` exists and is non-empty."""
    return (f'if [ ! -s "{path}" ]; then '
            f'echo "FAILED: DIA-NN exited 0 but did not write {what}: {path}" >&2; '
            f'echo "(DIA-NN 2.6 returns 0 on fatal errors -- check the log above for ERROR:)" >&2; '
            f'exit 1; fi')


def header(name, cpus, mem_gb, hours, partition, account, qos=None, array=None):
    h = ["#!/bin/bash -l",
         f"#SBATCH --job-name={name}",
         f"#SBATCH --cpus-per-task={cpus}",
         f"#SBATCH --mem={mem_gb}G",
         f"#SBATCH --time={hours}:00:00",
         f"#SBATCH --partition={partition}",
         f"#SBATCH --account={account}"]
    if qos:
        h.append(f"#SBATCH --qos={qos}")
    # publicgrp/low is PREEMPTIBLE: without --requeue a preempted task is simply lost.
    if (qos or "").startswith("public") or partition == "low":
        h.append("#SBATCH --requeue")
    h += [f"#SBATCH -o {name}_%j.log", f"#SBATCH -e {name}_%j.log"]
    if array:
        h.insert(2, f"#SBATCH --array={array}")
        h = [x.replace("_%j.log", "_%A_%a.log") for x in h]
    return "\n".join(h)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diann", required=True, help="DIA-NN command (native binary path, or 'apptainer exec --bind … <sif> /diann-*/diann-linux')")
    ap.add_argument("--raw", nargs="+", default=[], help="raw paths/globs (or use --raw-list)")
    ap.add_argument("--raw-list", help="file with one raw path per line — handles spaces in paths")
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
    ap.add_argument("--seed-lib", help="Skip step-1 prediction; use this existing library "
                    "(e.g. an InfinDIA empirical .parquet/.speclib) as the first-pass seed. "
                    "This is how you PARALLELIZE a semi-tryptic / non-specific / InfinDIA search: "
                    "build the small empirical library once (InfinDIA --pre-search), then fan the "
                    "per-file passes out across the cluster against it.")
    ap.add_argument("--seed-dep", help="SLURM job id the first pass should wait for (afterok) — "
                    "e.g. the InfinDIA lib-build job that produces --seed-lib.")
    ap.add_argument("--partition")
    ap.add_argument("--account")
    ap.add_argument("--qos", default=None, help="SLURM QOS. Needed for publicgrp/low "
                    "(publicgrp-low-qos); high/genome-center-grp uses its default.")
    ap.add_argument("--max-simultaneous", type=int, default=20)
    ap.add_argument("--no-norm", action="store_true")
    ap.add_argument("--probe-window", action=argparse.BooleanOptionalAction, default=True,
                    help="measure the scan-window radius after step 1 and pin it for all "
                         "steps (default: on). --no-probe-window requires --window in the cfg.")
    ap.add_argument("--allow-auto-mass-acc", action="store_true",
                    help="proceed even though the cfg leaves mass accuracy on auto. Results "
                         "across steps become inconsistent -- only for deliberate testing.")
    a = ap.parse_args()

    raws = []
    if a.raw_list:
        with open(a.raw_list) as fh:
            raws.extend(line.strip() for line in fh if line.strip())
    for p in a.raw:
        raws.extend(sorted(glob.glob(p)) or [p])
    raws = [os.path.abspath(r.rstrip("/")) for r in raws]

    # Detect the queue from the submitting user's SLURM associations rather than
    # assuming facility membership. genome-center-grp/high for members; publicgrp/low
    # for everyone else (incl. class accounts) — where `high` caps at 8 CPUs/job, so a
    # 32-CPU request would never start.
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from run_search import slurm_queue
        a.partition, a.account, _q = slurm_queue(a.partition, a.account, None)
        if not a.qos and a.partition == "low" and a.account == "publicgrp":
            a.qos = "publicgrp-low-qos"
    except Exception:
        a.partition = a.partition or "low"
        a.account = a.account or "publicgrp"
    n = len(raws)
    if n < 2:
        sys.exit("Parallel search needs >= 2 raw files (pass --raw or --raw-list); "
                 "use the single-shot run_search.py for 1.")
    out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    fasta = os.path.abspath(a.fasta)
    DN = dotnet_prefix(raws) + a.diann     # .NET 8 export prefix when inputs are Thermo .raw
    flags = read_cfg_flags(a.cfg)

    # The chain is only valid with pinned mass accuracy AND scan window (steps 3/5
    # reuse .quant files; see mass_acc_status for DIA-NN's own warning).
    ma = mass_acc_status(a.cfg)
    # An unpinned --window is recoverable: the chain can MEASURE it after step 1 and
    # feed the same value to steps 2 and 4 (--probe-window, on by default). Only an
    # unpinned MASS ACCURACY is unrecoverable here, because DIA-NN calibrates that per
    # run against the library and there is no single value to carry forward.
    win_probe = (a.probe_window and ma.get("window") in (None, 0)
                 and ma["ms2"] not in (None, 0) and ma["ms1"] not in (None, 0)
                 and not a.seed_lib)
    if win_probe:
        ma = dict(ma, fixed=True,
                  reason=f"MS1 {ma['ms1']} ppm / MS2 {ma['ms2']} ppm; "
                         "scan window measured by step 1b and pinned for steps 2+4")
    if not ma["fixed"] and not a.allow_auto_mass_acc:
        sys.exit(
            f"Not parallel-safe: {a.cfg or '(no --cfg given)'} -- {ma['reason']}.\n"
            "The 5-step chain reuses .quant files across steps, so anything DIA-NN\n"
            "auto-optimises PER FILE (mass accuracy AND scan window) is applied\n"
            "inconsistently between passes and then stitched together.\n"
            "  mass accuracy -> re-run estimate_params.py with the real instrument\n"
            "                   (timsTOF 15/15, Astral 4/10, Orbitrap by resolution)\n"
            "  --window      -> run `probe_window.py` on ONE of these files to measure\n"
            "                   the radius DIA-NN infers, then pin that value. Do NOT\n"
            "                   guess it: it depends on the acquisition scheme.\n"
            "Or run the single-shot search instead, which calibrates safely on its own.\n"
            "To override deliberately: --allow-auto-mass-acc.")
    if not ma["fixed"]:
        sys.stderr.write(f"[diann_parallel] WARNING: proceeding with auto mass accuracy "
                         f"({ma['reason']}) -- steps will not be mutually consistent.\n")
    D = out  # all DIA-NN intermediate/output lives here (real paths; native binary reads them directly)
    report = "no_norm_report.parquet" if a.no_norm else "report.parquet"
    norm = "--no-norm" if a.no_norm else ""
    xic = xic_flag(a.cfg)          # step 4 only -- see xic_flag() docstring

    # file list (1 raw path per line) — array tasks index into it
    open(os.path.join(out, "file_list.txt"), "w").write("\n".join(raws) + "\n")
    all_f = " ".join(f"--f {shlex.quote(r)}" for r in raws)   # quote — data paths may contain spaces
    array = f"0-{n-1}%{a.max_simultaneous}"
    seed = os.path.abspath(a.seed_lib) if a.seed_lib else None
    predicted = seed if seed else f"{D}/step1.predicted.speclib"
    empirical = f"{D}/empirical.parquet"

    # QUEUE CHOICE, ported from DE-LIMP's select_best_partition()
    # (R/helpers_search.R) + docs/QUEUE_SWITCHING.md:
    #   steps 2 & 4 are ARRAYS -- embarrassingly parallel, one file per task, so a
    #     preemption costs one task. They are also exactly what the per-user CPU cap on
    #     the priority queue throttles, so prefer publicgrp/low when it has idle CPUs
    #     (it routinely has thousands).
    #   steps 1, 3 & 5 are SINGLE jobs that cannot restart mid-way (3 consumes every
    #     .quant file, 5 needs all of step 4), so a preemption throws the stage away --
    #     keep them on the priority queue.
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from run_search import slurm_queue as _sq
        _pa, _aa, _qa = _sq(a.partition, a.account, a.qos,
                            peak_cpus=a.threads_per_file, preemptible_ok=True)
        _ps, _as_, _qs = _sq(a.partition, a.account, a.qos, peak_cpus=a.assembly_cpus)
    except Exception:
        _pa, _aa, _qa = a.partition, a.account, a.qos
        _ps, _as_, _qs = a.partition, a.account, a.qos

    def write(name, body):
        p = os.path.join(out, name)
        open(p, "w").write(body + "\n"); os.chmod(p, 0o755); return name

    # array preamble: pick this task's raw file
    pick = ('FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ' + f'{D}/file_list.txt)\n'
            'if [ -z "$FILE" ]; then echo "no file for task $SLURM_ARRAY_TASK_ID"; exit 1; fi\n'
            'echo "Processing: $FILE"\n')

    # Step 1 — library prediction (single job) — SKIPPED when --seed-lib is given
    s1 = None
    if not seed:
        s1 = write("step1_libpred.sbatch", "\n".join([
            header("s1_libpred", a.libpred_cpus, a.libpred_mem, a.libpred_time, _ps, _as_, qos=_qs), "",
            f'echo "Step 1/5 library prediction"; date',
            f'{DN} --fasta {fasta} --fasta-search --predictor --gen-spec-lib \\',
            f'  --out-lib {D}/step1.speclib --out {D}/step1_lib.parquet \\',
            f'  --threads {a.libpred_cpus} {flags}',
            must_exist(predicted, "the predicted spectral library")]))

    # Step 1b — measure the scan-window radius ONCE, so steps 2 and 4 share it.
    # DIA-NN optimises the radius per file when --window is absent, and steps 3/5 then
    # combine .quant files produced under different windows -- which DIA-NN's own
    # warning calls "strongly not recommended". Measuring beats guessing: the radius
    # depends on the acquisition scheme (cycle time vs peak width), not the instrument.
    s1b = None
    if win_probe:
        probe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "probe_window.py")
        s1b = write("step1b_window.sbatch", "\n".join([
            header("s1b_window", a.threads_per_file, a.mem_per_file, 2, _ps, _as_, qos=_qs), "",
            'echo "Step 1b/5 measuring scan-window radius"; date',
            f'python3 "{probe}" --diann {shlex.quote(a.diann)} --raw "{raws[0]}" \\',
            f'  --fasta {fasta} --lib {predicted} --threads {a.threads_per_file} \\',
            f'  --extra {shlex.quote(read_cfg_flags(a.cfg))} > {D}/window.json',
            f'python3 -c "import json;print(json.load(open(\'{D}/window.json\'))[\'window_radius\'])" > {D}/window.txt',
            must_exist(f"{D}/window.txt", "the measured scan-window radius"),
            f'echo "scan window radius = $(cat {D}/window.txt) (pinned for steps 2 and 4)"']))
    # steps 2 and 4 read the measured radius at RUNTIME so both use the identical value
    wflag = f'--window $(cat {D}/window.txt) ' if win_probe else ''

    # Step 2 — first pass (array): predicted lib -> per-file .quant
    s2 = write("step2_firstpass.sbatch", "\n".join([
        header("s2_firstpass", a.threads_per_file, a.mem_per_file, a.time_per_file, _pa, _aa, qos=_qa, array=array), "",
        f'echo "Step 2/5 first pass, task ${{SLURM_ARRAY_TASK_ID}} of {n}"; date', pick,
        f'{DN} --f "$FILE" --fasta {fasta} --lib {predicted} \\',
        f'  --temp {D}/quant_step2 --rt-profiling --gen-spec-lib --quant-ori-names \\',
        f'  --threads {a.threads_per_file} {wflag}{flags}',
        'QOUT="${FILE##*/}"; QOUT="${QOUT%.*}.quant"',
        must_exist(f'{D}/quant_step2/$QOUT', "this file's .quant")]))

    # Step 3 — empirical library assembly (single job, --use-quant)
    s3 = write("step3_assembly.sbatch", "\n".join([
        header("s3_assembly", a.assembly_cpus, a.assembly_mem, a.assembly_time, _ps, _as_, qos=_qs), "",
        f'echo "Step 3/5 empirical library assembly"; date',
        f'cp -r {D}/quant_step2 {D}/quant_step2_orig 2>/dev/null || true   # backup for resume',
        f'{DN} {all_f} --fasta {fasta} --lib {predicted} --use-quant --quant-ori-names \\',
        f'  --rt-profiling --gen-spec-lib --out-lib {empirical} \\',
        f'  --temp {D}/quant_step2 --out {D}/step3_assembly.parquet \\',
        f'  --threads {a.assembly_cpus} {wflag}{flags}',
        must_exist(empirical, "the empirical spectral library")]))

    # Step 4 — final pass (array): empirical lib -> per-file .quant
    # --xic alone is not enough. DIA-NN names the XIC folder after the --out report
    # basename; with no --out it resolves that against the filesystem ROOT and dies with
    # "cannot create directory: Permission denied [/report_xic]" -- AFTER writing the
    # .quant, and STILL EXITING 0. Verified on DIA-NN 2.6.0 (single file, --xic, no
    # --out): zero .xic.parquet produced, SLURM records COMPLETED, afterok advances.
    # So step 4 needs its own --out whenever XICs are requested. DIA-NN then writes
    # <out>/xic/t<TASKID>_xic/<run>.xic.parquet, one folder per array task.
    # Both carry their own leading space so the command line has no double space (and no
    # trailing space before the line continuation) when XICs are off.
    xic_arg = f' {xic}' if xic else ''
    xic_out = f' --out {D}/xic/t${{SLURM_ARRAY_TASK_ID}}.parquet' if xic else ''

    s4 = write("step4_finalpass.sbatch", "\n".join([
        header("s4_finalpass", a.threads_per_file, a.mem_per_file, a.time_per_file, _pa, _aa, qos=_qa, array=array), "",
        f'echo "Step 4/5 final pass, task ${{SLURM_ARRAY_TASK_ID}} of {n}"; date', pick,
        'QUANT="${FILE##*/}"; QUANT="${QUANT%.*}.quant"',
        f'if [ ! -f "{D}/quant_step2/$QUANT" ]; then echo "SKIP: no step-2 quant for $QUANT"; exit 0; fi',
        # Splat an empty list, not an empty string: a conditional string leaves a stray
        # blank line in the generated sbatch when XICs are off.
        *([f'mkdir -p {D}/xic'] if xic else []),
        # DIA-NN re-saves the library it is handed as "<lib>.skyline.speclib", written
        # NEXT TO --lib. With one shared path every concurrent array task writes the same
        # file; most win the race in seconds, the losers block until the wall clock kills
        # them. Give each task its own copy so there is nothing to contend on.
        f'LIBPRIV={D}/libpriv/t${{SLURM_ARRAY_TASK_ID}}',
        'mkdir -p "$LIBPRIV"',
        f'cp -f {empirical} "$LIBPRIV/lib.parquet"',
        'trap \'rm -rf "$LIBPRIV"\' EXIT',
        f'{DN} --f "$FILE" --fasta {fasta} --lib "$LIBPRIV/lib.parquet" \\',
        f'  --temp {D}/quant_step4 --no-ifs-removal --quant-ori-names{xic_arg}{xic_out} \\',
        f'  --threads {a.threads_per_file} {wflag}{flags}',
        must_exist(f'{D}/quant_step4/$QUANT', "this file's final-pass .quant")]))

    # Step 5 — cross-run report (single job, --use-quant --matrices)
    s5 = write("step5_report.sbatch", "\n".join([
        header("s5_report", a.assembly_cpus, a.assembly_mem, a.assembly_time, _ps, _as_, qos=_qs), "",
        f'echo "Step 5/5 cross-run report"; date',
        f'{DN} {all_f} --fasta {fasta} --lib {empirical} --use-quant --quant-ori-names \\',
        f'  --temp {D}/quant_step4 --matrices --out {D}/{report} \\',
        f'  --threads {a.assembly_cpus} {norm} {wflag}{flags}',
        must_exist(f'{D}/{report}', "the cross-run report"),
        # A step-4 task that failed silently leaves no .quant, and step 5 happily
        # reports on whatever survived. Count them: fewer quants than inputs means a
        # sample was dropped, which must never pass as success.
        f'NQ=$(ls -1 {D}/quant_step4/*.quant 2>/dev/null | wc -l | tr -d " ")',
        f'if [ "$NQ" -ne {n} ]; then '
        f'echo "FAILED: report built from $NQ of {n} runs -- a step-4 task produced no .quant." >&2; '
        f'echo "Find it: for f in \\$(cat {D}/file_list.txt); do b=\\$(basename \\"\\$f\\"); '
        f'[ -f {D}/quant_step4/\\${{b%.*}}.quant ] || echo MISSING \\$b; done" >&2; '
        f'exit 1; fi',
        f'echo "OK: report built from all {n} runs"']))

    # submit.sh — chain the steps with afterok dependencies
    sub_lines = ["#!/bin/bash", "set -euo pipefail", f'cd "{out}"',
                 f'mkdir -p "{D}/quant_step2" "{D}/quant_step4"   # DIA-NN --temp dirs MUST pre-exist']
    if seed:
        # Step 1 skipped — the InfinDIA/empirical seed library IS the first-pass lib.
        # First pass optionally waits (afterok) on the lib-build job that produces it.
        dep2 = f'--dependency=afterok:{a.seed_dep} ' if a.seed_dep else ''
        sub_lines += [
            f'echo "Step 1/5 SKIPPED — seeding first pass with {predicted}"',
            'jid2=$(sbatch --parsable %s%s)' % (dep2, s2)]
    else:
        sub_lines += ['jid1=$(sbatch --parsable %s)' % s1]
        if s1b:
            sub_lines += [
                'jid1b=$(sbatch --parsable --dependency=afterok:$jid1 %s)' % s1b,
                'jid2=$(sbatch --parsable --dependency=afterok:$jid1b %s)' % s2]
        else:
            sub_lines += ['jid2=$(sbatch --parsable --dependency=afterok:$jid1 %s)' % s2]
    sub_lines += [
        'jid3=$(sbatch --parsable --dependency=afterok:$jid2 %s)' % s3,
        'jid4=$(sbatch --parsable --dependency=afterok:$jid3 %s)' % s4,
        'jid5=$(sbatch --parsable --dependency=afterok:$jid4 %s)' % s5,
        'echo "submitted: firstpass=$jid2 assembly=$jid3 finalpass=$jid4 report=$jid5"',
        f'echo "final report will be {D}/{report}; watch with: watch_run.sh --slurm $jid5 --log {D}/s5_report_${{jid5}}.log"']

    # Record the submission so a disconnected session can pick the run back up.
    # SLURM keeps the jobs alive after the user closes their terminal; without this
    # the next session has no idea what was submitted or what is still outstanding.
    # <out> is normally <session>/output/search, so the session is two levels up.
    sess = os.path.dirname(os.path.dirname(out))
    ck = os.path.join(os.path.dirname(os.path.abspath(__file__)), "checkpoint.py")
    jobs = '$jid2,$jid3,$jid4,$jid5' if seed else '$jid1,$jid2,$jid3,$jid4,$jid5'
    sub_lines += [
        f'python3 "{ck}" record --session "{sess}" --stage search \\',
        f'  --jobs "{jobs}" --desc "DIA-NN 5-step parallel chain ({n} files)" \\',
        f'  --report "{D}/{report}" --watch-job "$jid5" --watch-log "{D}/s5_report_${{jid5}}.log" \\',
        f'  --next "Rscript run_de.R --input {D}/{report} --metadata {sess}/input/conditions.csv --method dpc --outdir {sess}/output/tables" \\',
        '  >/dev/null 2>&1 || true',
        f'printf "%s\\n" {jobs.replace(",", " ")} > "{D}/jobs.txt"',
        f'echo "all chain job ids -> {D}/jobs.txt  (watch the WHOLE chain: watch_run.sh --all {D})"',
        f'echo "recovery notes written to {sess}/RECOVERY.md — you can safely close your terminal"']
    write("submit.sh", "\n".join(sub_lines))

    import json
    print(json.dumps({
        "out": out, "n_files": n, "report": f"{D}/{report}",
        "mass_acc": ma,
        "seeded": bool(seed), "seed_lib": predicted if seed else None,
        "scripts": [x for x in [s1, s2, s3, s4, s5, "submit.sh"] if x],
        "submit": f"bash {out}/submit.sh   (or: hive_exec.sh 'bash {out}/submit.sh')",
        "report_jobid_var": "jid5",
        "note": "5-step DIA-NN parallel chain. Mass accuracy is fixed (manual) — ensure --cfg has --mass-acc/--mass-acc-ms1 set, not 0/auto. Watch the run with watch_run.sh. Then point run_de.R at the report.",
    }, indent=2))


if __name__ == "__main__":
    main()
