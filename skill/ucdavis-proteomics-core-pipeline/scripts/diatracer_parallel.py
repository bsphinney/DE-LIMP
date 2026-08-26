#!/usr/bin/env python3
"""
diatracer_parallel.py -- run diaTracer conversion as a SLURM array, one task per file.

WHY
---
FragPipe runs diaTracer serially inside its own job: ~20 min per .d (paper: 34 files,
32 cores, 19.4 min each). A nine-file cohort is therefore ~3 hours of wall-clock before
MSFragger even starts. But the conversion is per-file and completely independent, so on
a cluster it should be an array: nine tasks at once finish in the time of the slowest
single file.

The pseudo-MS/MS mzML are reusable. Once they exist, diatracer_stage.py emits the reuse
form of the manifest (mzML as DDA + the original .d as DIA-Quant), and the subsequent
FragPipe run skips conversion entirely and goes straight to MSFragger.

  stage symlinks -> THIS (array, ~1 file-time) -> re-stage -> FragPipe (no conversion)

⚠ diaTracer's standalone main() parses every numeric option UNCONDITIONALLY and does
NOT apply the defaults its own --help advertises. Omit one and it dies immediately with
"NullPointerException: Cannot invoke String.trim() because in is null" inside
Double.parseDouble. So every option below is passed explicitly, at the documented
default. Do not "simplify" this by dropping them.

Usage
  python3 diatracer_parallel.py --stage <staging dir> --out <dir for scripts+logs> \\
      [--jar <diatracer.jar>] [--threads 16] [--mem 96] [--time 6] \\
      [--partition low] [--account publicgrp] [--max-simultaneous 9]
  bash <out>/submit_diatracer.sh
"""
import argparse, glob, json, os, shutil, sys

# documented defaults from `java -jar diatracer.jar --help` (diaTracer 2.2.1)
DEFAULTS = {"-dI": "0.01", "-dR": "3", "-mC": "0.3",
            "-mF": "1", "-mO": "0.1", "-r": "0", "-rM": "500"}

HIVE_JAR = "/quobyte/proteomics-grp/fragpipe24/fragpipe-24.0/tools/diatracer-2.2.1.jar"


def find_jar(explicit):
    if explicit:
        return explicit
    if os.path.exists(HIVE_JAR):
        return HIVE_JAR
    env = os.environ.get("FRAGPIPE_TOOLS_FOLDER")
    if env:
        hits = sorted(glob.glob(os.path.join(env, "**", "diatracer*.jar"), recursive=True))
        if hits:
            return hits[-1]
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True, help="staging dir of .d symlinks (diatracer_stage.py)")
    ap.add_argument("--out", required=True, help="where to write the sbatch + logs")
    ap.add_argument("--jar")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--mem", type=int, default=96, help="GB per task")
    ap.add_argument("--time", type=int, default=6, help="hours per task")
    # No hardcoded queue. Defaults of None mean "detect it" (see below) -- `low`/publicgrp was
    # wrong in both directions: it put facility members on the PREEMPTIBLE queue when they had a
    # non-preemptible one, and it emitted no --qos, which publicgrp/low requires on HIVE.
    ap.add_argument("--partition", default=None)
    ap.add_argument("--account", default=None)
    ap.add_argument("--qos", default=None,
                    help="SLURM QOS. publicgrp/low needs publicgrp-low-qos; "
                         "genome-center-grp/high uses its default. Detected when omitted.")
    ap.add_argument("--max-simultaneous", type=int, default=0,
                    help="cap concurrent tasks (0 = all at once)")
    a = ap.parse_args()

    stage = os.path.abspath(a.stage)
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)

    dfiles = sorted(glob.glob(os.path.join(stage, "*.d")))
    if not dfiles:
        sys.exit(f"no .d entries in {stage} — run diatracer_stage.py first")

    todo = [d for d in dfiles if not os.path.exists(d[:-2] + "_diatracer.mzML")]
    n = len(dfiles)

    jar = find_jar(a.jar)
    if not jar:
        sys.exit("diaTracer jar not found. Pass --jar, or set FRAGPIPE_TOOLS_FOLDER. "
                 "It is licence-gated (MSFragger/IonQuant/diaTracer) and cannot be "
                 "downloaded by script.")

    if not todo:
        print(json.dumps({"stage": stage, "n_files": n, "n_to_convert": 0,
                          "note": "all pseudo-spectra already present — skip straight to FragPipe"},
                         indent=2))
        return

    # ONE definition of the queue: reuse run_search.slurm_queue() rather than re-deriving it.
    # It asks SLURM what this account is actually entitled to -- genome-center-grp/high for
    # facility members, publicgrp/low for everyone else. A hardcoded facility queue is REJECTED
    # outright for an account outside the group, and a hardcoded `low` needlessly exposes a
    # facility member's 20-min-per-file conversion to preemption.
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from run_search import slurm_queue
        a.partition, a.account, _q = slurm_queue(a.partition, a.account, a.qos)
        a.qos = a.qos or _q
        if not a.qos and a.partition == "low" and a.account == "publicgrp":
            a.qos = "publicgrp-low-qos"
    except Exception as e:                 # no SLURM/sacctmgr: emit nothing and let SLURM decide
        # Say so. Emitting no queue is the safe choice, but it is NOT the same as having chosen
        # one, and a facility member should know they may land on the cluster default.
        sys.stderr.write(f"[diatracer_parallel] queue detection unavailable "
                         f"({type(e).__name__}: {e}); emitting no partition/account/qos and "
                         f"letting SLURM apply its default\n")

    cap = f"%{a.max_simultaneous}" if a.max_simultaneous else ""
    acct = f"#SBATCH --account={a.account}\n" if a.account else ""
    part = f"#SBATCH --partition={a.partition}\n" if a.partition else ""
    qos = f"#SBATCH --qos={a.qos}\n" if a.qos else ""
    # --requeue only where it means something: `low` is preemptible, so without it a preempted
    # conversion is simply lost. On a non-preemptible queue it is noise.
    rq = "#SBATCH --requeue\n" if a.partition == "low" else ""
    opts = " ".join(f"{k} {v}" for k, v in DEFAULTS.items())

    script = f"""#!/bin/bash -l
#SBATCH --job-name=diatracer
#SBATCH --array=0-{n - 1}{cap}
#SBATCH --cpus-per-task={a.threads}
#SBATCH --mem={a.mem}G
#SBATCH --time={a.time}:00:00
{part}{acct}{qos}{rq}#SBATCH -o {out}/diatracer_%A_%a.log
#SBATCH -e {out}/diatracer_%A_%a.log
set -uo pipefail

mapfile -t DFILES < <(find "{stage}" -maxdepth 1 -name '*.d' | sort)
D="${{DFILES[$SLURM_ARRAY_TASK_ID]}}"
[ -n "${{D:-}}" ] || {{ echo "no .d for task $SLURM_ARRAY_TASK_ID"; exit 1; }}
OUT="${{D%.d}}_diatracer.mzML"

# reusable: never redo a conversion that already succeeded
if [ -s "$OUT" ]; then echo "SKIP (already converted): $OUT"; exit 0; fi

echo "=== task $SLURM_ARRAY_TASK_ID  host $(hostname)  $(date -u +%FT%TZ) ==="
echo "input : $D"

# every numeric option passed explicitly — diaTracer does NOT apply its own defaults
java -Xmx{max(a.mem - 6, 8)}G -jar "{jar}" -d "$D" -w "{stage}" \\
  -t "${{SLURM_CPUS_PER_TASK:-{a.threads}}}" {opts}
rc=$?

echo "=== finished $(date -u +%FT%TZ) exit=$rc ==="
if [ -s "$OUT" ]; then ls -la "$OUT"; else echo "WARNING: expected output missing"; fi
exit $rc
"""
    sp = os.path.join(out, "diatracer_array.sbatch")
    with open(sp, "w") as fh:
        fh.write(script)
    os.chmod(sp, 0o755)

    sub = os.path.join(out, "submit_diatracer.sh")
    with open(sub, "w") as fh:
        fh.write(f"""#!/bin/bash
set -euo pipefail
jid=$(sbatch --parsable "{sp}")
echo "diatracer array = $jid"
echo "watch: python3 checkpoint.py progress --session <session>"
""")
    os.chmod(sub, 0o755)

    print(json.dumps({
        "stage": stage, "out": out, "jar": jar,
        "n_files": n, "n_to_convert": len(todo),
        "n_already_converted": n - len(todo),
        "sbatch": sp, "submit": f"bash {sub}",
        "expected_outputs": [d[:-2] + "_diatracer.mzML" for d in dfiles],
        "note": ("One SLURM task per file instead of FragPipe's serial loop: wall-clock "
                 "becomes the slowest single file, not the sum. Re-stage afterwards so "
                 "the manifest takes the reuse form and FragPipe skips conversion."),
    }, indent=2))


if __name__ == "__main__":
    main()
