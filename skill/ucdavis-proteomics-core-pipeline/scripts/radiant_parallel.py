#!/usr/bin/env python3
"""
radiant_parallel.py -- Run Radiant DIA per-file as a SLURM array, then rescore the
whole set in ONE Fulcrum job.

WHY THIS EXISTS
---------------
Radiant's own Fulcrum search backend is serial: `radiant_fulcrum_search.search`
raises `NotImplementedError("TODO: parallel processing")` for anything but
`mode="serial"`, so a 20-file study searches one file at a time even on a 64-core
node. But the search IS embarrassingly parallel -- each file is independent, and
only the downstream stages (rescoring, FDR, protein inference, rollup) need to see
every file at once.

So we split it exactly there:

  step 1  DIA-NN predicts the spectral library                     (1 job)
  step 2  RadiantDIA runs on ONE file per array task               (N tasks)
  step 3  Fulcrum consumes the per-file results and does the       (1 job)
          rescoring / FDR / inference / rollup / output

Step 3 does NOT re-search. `_execute_radiant()` computes its result path as
`<output_location>/<input name>.radiantDIA` and returns it untouched when
`reuse_existing` is set, so pointing Fulcrum at the directory the array wrote makes
it skip straight to the downstream stages. That is the whole trick, and it is why
step 2 must write into the same `radiant-results/` directory Fulcrum expects.

Emits scripts only -- it never submits. The orchestrator submits with dependencies
so heavy compute never touches a login node.

Usage:
  python3 radiant_parallel.py --runtime docker --image seerbio/radiant-fulcrum:2.3.3 \
      --raw-list files.txt --fasta search.fasta --config default.radiantConfig \
      --out ./search_out [--library lib.tsv] [--diann '<diann cmd>'] [--threads-per-file 16]
"""
import sys, os, json, argparse, shlex

HERE = os.path.dirname(os.path.abspath(__file__))


def _mounts_for(runtime, pairs):
    flag = "-v" if runtime == "docker" else "--bind"
    return " ".join(f"{flag} {shlex.quote(h)}:{c}" for h, c in pairs)


def container_prefix(runtime, image, pairs):
    """docker/apptainer invocation with the given bind mounts."""
    if runtime == "docker":
        return f"docker run --rm {_mounts_for(runtime, pairs)} {shlex.quote(image)}"
    return f"apptainer exec {_mounts_for(runtime, pairs)} {shlex.quote(image)}"


def detect_queue(partition, account, qos):
    """Resolve the queue from the submitting user's own SLURM associations.

    ONE definition: reuse run_search.slurm_queue() rather than re-deriving it. Emitting
    no partition at all is not neutral -- SLURM then applies the cluster default, which
    on HIVE is `high`, precisely the queue a non-facility account cannot submit to (the
    job is REJECTED, not merely slowed). That is the bug emit_sbatch was fixed for; this
    generator must not reintroduce it."""
    try:
        sys.path.insert(0, HERE)
        from run_search import slurm_queue
        return slurm_queue(partition, account, qos)
    except Exception as e:
        sys.stderr.write(f"[radiant_parallel] queue detection unavailable ({e}); "
                         "emitting the caller's values verbatim\n")
        return partition, account, qos


def header(job, cpus, mem, hours, partition, account, qos, array=None):
    L = ["#!/bin/bash -l", f"#SBATCH --job-name={job}",
         f"#SBATCH --cpus-per-task={cpus}", f"#SBATCH --mem={mem}G",
         f"#SBATCH --time={hours}:00:00",
         f"#SBATCH --output=%x_%j.log"]
    if partition:
        L.append(f"#SBATCH --partition={partition}")
    if account:
        L.append(f"#SBATCH --account={account}")
    if qos:
        L.append(f"#SBATCH --qos={qos}")
        if qos.startswith("public") or qos == "low":
            # publicgrp/low is preemptible -- requeue rather than losing the task.
            L.append("#SBATCH --requeue")
    if array:
        L.append(f"#SBATCH --array={array}")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runtime", required=True, choices=["docker", "apptainer"])
    ap.add_argument("--image", required=True, help="image ref or .sif path")
    ap.add_argument("--raw-list", required=True, help="file with one mzML path per line")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--config", required=True, help="the .radiantConfig to search with")
    ap.add_argument("--out", required=True)
    ap.add_argument("--library", help="DIA-NN .tsv library; generated in step 1 if omitted")
    ap.add_argument("--diann", help="DIA-NN command (needed only when --library is omitted)")
    # Memory defaults are MEASURED, not guessed (3-run HeLa, 9.9 GB TSV library,
    # 2026-08-03): each RadiantDIA task peaked at ~19 GB, Fulcrum at ~14 GB — roughly
    # 2x the library size plus overhead. Over-requesting is not free: asking 200 GB
    # made SLURM schedule the array tasks one at a time, turning a parallel array back
    # into a near-serial queue (14 min wall vs 21 min serial, instead of ~7). Scale
    # --mem-per-file with the library if you use a much larger one.
    ap.add_argument("--threads-per-file", type=int, default=16)
    ap.add_argument("--mem-per-file", type=int, default=32)
    ap.add_argument("--time-per-file", type=int, default=4)
    ap.add_argument("--fulcrum-cpus", type=int, default=32)
    ap.add_argument("--fulcrum-mem", type=int, default=48)
    ap.add_argument("--fulcrum-time", type=int, default=8)
    ap.add_argument("--libpred-cpus", type=int, default=16)
    ap.add_argument("--libpred-mem", type=int, default=64)
    ap.add_argument("--libpred-time", type=int, default=4)
    ap.add_argument("--partition"); ap.add_argument("--account"); ap.add_argument("--qos")
    ap.add_argument("--max-simultaneous", type=int, default=20)
    ap.add_argument("--fdr-thresh", type=float, default=0.01)
    ap.add_argument("--mbr", action="store_true")
    a = ap.parse_args()

    files = [ln.strip() for ln in open(a.raw_list) if ln.strip()]
    if not files:
        sys.exit(f"{a.raw_list} listed no input files.")
    bad = [f for f in files if f.rstrip("/").lower().endswith(".d")]
    if bad:
        sys.exit("Radiant cannot read Bruker .d — convert to mzML or use DIA-NN/diaTracer.\n"
                 f"  offending: {bad}")

    out = os.path.abspath(a.out)
    D = os.path.join(out, "radiant_parallel")
    results = os.path.join(out, "radiant_results")
    # This exact directory is the contract between step 2 and step 3: Fulcrum's
    # radiant workflow looks for reusable results in <results_dir>/radiant-results.
    rdir = os.path.join(results, "radiant-results")
    os.makedirs(D, exist_ok=True); os.makedirs(rdir, exist_ok=True)

    # DIA-NN IGNORES the extension given to --out-lib for a PREDICTED library and always
    # writes <stem>.predicted.speclib (its compact binary format). Radiant reads .speclib
    # directly (FragLibReader accepts .fragLibFF/.tsv/.csv/.speclib), so no conversion is
    # needed -- but the filename must match what DIA-NN really produces or step 2 looks
    # for a file that was never written.
    lib_stem = os.path.join(D, "diann_predicted_lib")
    lib = a.library or (lib_stem + ".predicted.speclib")
    listing = os.path.join(D, "file_list.txt")
    with open(listing, "w") as fh:
        fh.write("\n".join(os.path.abspath(f) for f in files) + "\n")

    # Bind every distinct host dir the container must see. Mount parents rather than
    # copying so large mzML stay where they are.
    dirs = {os.path.dirname(os.path.abspath(f)) for f in files}
    dirs |= {os.path.dirname(os.path.abspath(x)) for x in (a.fasta, a.config, lib)}
    pairs, cmap = [], {}
    for i, d in enumerate(sorted(dirs)):
        cmap[d] = f"/mnt/in{i}"
        pairs.append((d, cmap[d]))
    pairs.append((results, "/mnt/results"))
    cres = "/mnt/results"

    def cp(host):
        h = os.path.abspath(host)
        return f"{cmap[os.path.dirname(h)]}/{os.path.basename(h)}"

    prefix = container_prefix(a.runtime, a.image, pairs)
    part, acct, qos = detect_queue(a.partition, a.account, a.qos)
    q = dict(partition=part, account=acct, qos=qos)

    def write(name, text):
        p = os.path.join(D, name)
        with open(p, "w") as fh:
            fh.write(text + "\n")
        os.chmod(p, 0o755)
        return p

    scripts = {}

    # ---- step 1: DIA-NN predicted library (skipped when --library is given) -----
    if not a.library:
        if not a.diann:
            sys.exit("--library not given, so step 1 must generate one, but --diann was "
                     "not provided. Radiant needs a DIA-NN library (it reads DIA-NN's TSV "
                     "schema natively).")
        scripts["step1_libpred"] = write("step1_libpred.sbatch", "\n".join([
            header("r1_libpred", a.libpred_cpus, a.libpred_mem, a.libpred_time, **q), "",
            'echo "Step 1/3 DIA-NN library prediction"; date',
            f'{a.diann} --fasta {shlex.quote(os.path.abspath(a.fasta))} --fasta-search '
            f'--predictor --gen-spec-lib --out-lib {shlex.quote(lib_stem)} '
            f'--threads {a.libpred_cpus}',
            # Check for what DIA-NN ACTUALLY writes, and say what it did write if absent.
            f'if [ ! -s {shlex.quote(lib)} ]; then',
            f'  echo "expected {os.path.basename(lib)}; DIA-NN wrote:"; ls -la {shlex.quote(D)}',
            f'  exit 1',
            f'fi']))

    # ---- step 2: one RadiantDIA per file (array) --------------------------------
    n = len(files)
    array = f"0-{n-1}%{a.max_simultaneous}"
    # Two parallel lists: the array indexes the CONTAINER-path list directly. Rewriting
    # host->container paths in bash would need string substitution whose delimiter is
    # `/` -- the same character the paths are full of -- so we resolve it here instead.
    clisting = os.path.join(D, "file_list_container.txt")
    with open(clisting, "w") as fh:
        fh.write("\n".join(cp(f) for f in files) + "\n")

    pick = (f'FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" {shlex.quote(clisting)})\n'
            'if [ -z "$FILE" ]; then echo "no file for task $SLURM_ARRAY_TASK_ID"; exit 1; fi\n'
            'BASE=$(basename "$FILE")\n'
            'echo "Searching: $FILE"\n'
            # Idempotent: a requeued preemptible task must not redo finished work.
            f'if [ -s "{rdir}/$BASE.radiantDIA" ]; then echo "already searched; skipping"; exit 0; fi\n')
    # The container-path list lives on the host but is read by the SUBMITTING shell
    # (sed runs outside the container), so mount it nowhere -- it is a host file.
    scripts["step2_search"] = write("step2_search.sbatch", "\n".join([
        header("r2_search", a.threads_per_file, a.mem_per_file, a.time_per_file,
               array=array, **q), "",
        f'echo "Step 2/3 Radiant search, task ${{SLURM_ARRAY_TASK_ID}} of {n}"; date',
        pick,
        f'{prefix} RadiantDIA {cp(lib)} {cp(a.fasta)} {cp(a.config)} "$FILE" '
        f'--output-folder {cres}/radiant-results']))

    # ---- step 3: ONE Fulcrum job over all per-file results ----------------------
    # reuse_existing=true is what makes this a rescoring pass, not a re-search.
    toml = "\n".join([
        '# Generated by radiant_parallel.py — step 3 rescores the per-file results',
        '# produced by the step-2 array. reuse_existing skips the search itself.',
        'workflow = "radiant"',
        f'library = {json.dumps(cp(lib))}',
        f'fasta = {json.dumps(cp(a.fasta))}',
        f'config = {json.dumps(cp(a.config))}',
        'mzml = [' + ", ".join(json.dumps(cp(f)) for f in files) + ']',
        f'mbr = {"true" if a.mbr else "false"}',
        f'fdr_thresh = {a.fdr_thresh}',
        f'threads = {a.fulcrum_cpus}',
        f'results_dir = {json.dumps(cres)}',
        '',
        '[overrides.search]',
        'reuse_existing = true',
    ])
    toml_path = os.path.join(results, "fulcrum_rescore.toml")
    with open(toml_path, "w") as fh:
        fh.write(toml + "\n")

    scripts["step3_fulcrum"] = write("step3_fulcrum.sbatch", "\n".join([
        header("r3_fulcrum", a.fulcrum_cpus, a.fulcrum_mem, a.fulcrum_time, **q), "",
        'echo "Step 3/3 Fulcrum rescoring + FDR + rollup over all files"; date',
        f'N=$(ls -1 {shlex.quote(rdir)}/*.radiantDIA 2>/dev/null | wc -l)',
        f'echo "found $N per-file result(s), expected {n}"',
        f'if [ "$N" -lt {n} ]; then echo "MISSING per-file results — not rescoring a '
        f'partial set"; exit 1; fi',
        f'{prefix} fulcrum -v --toml-file {cres}/fulcrum_rescore.toml']))

    print(json.dumps({
        "mode": "radiant_parallel_3step",
        "n_files": n,
        "scripts": scripts,
        "results_dir": results,
        "per_file_results": rdir,
        "toml": toml_path,
        "library": lib,
        "queue": {"partition": part, "account": acct, "qos": qos},
        "submit": ("Submit with dependencies, e.g.:\n"
                   + ("  j1=$(sbatch --parsable step1_libpred.sbatch)\n"
                      "  j2=$(sbatch --parsable --dependency=afterok:$j1 step2_search.sbatch)\n"
                      if not a.library else
                      "  j2=$(sbatch --parsable step2_search.sbatch)\n")
                   + "  j3=$(sbatch --parsable --dependency=afterok:$j2 step3_fulcrum.sbatch)"),
        "note": ("Step 2 parallelises the part Radiant itself runs serially; step 3 is the "
                 "only stage that needs all files at once. Adapt the output with "
                 "run_search.py --adapt-only --engine radiant --out <out>."),
    }, indent=2))


if __name__ == "__main__":
    main()
