#!/usr/bin/env python3
"""
run_search.py  --  Run the selected search engine and normalize its output to
the DE-input contract (a DIA-NN-shaped report.parquet / matrix the DE step
consumes; see references/de-analysis.md §8.3).

Routing (PLAN.md §7b):
  default by acquisition: DIA -> diann, DDA -> sage; --engine overrides;
  FragPipe only when the bundle names it or the user asks.

Per engine:
  diann    <cmd> --cfg <bundle .cfg> --f <files> --fasta <fasta>
           --out report.parquet --threads N        (native contract, no adapter)
  sage     convert .d/.raw -> mzML if needed (msconvert), then
           <cmd> <bundle sage_config.json> -f <fasta> -o <out> --parquet
           --disable-telemetry-i-dont-want-to-improve-sage
           then adapt lfq.parquet -> DIA-NN-shaped report for --method maxlfq
  fragpipe <cmd> --headless --workflow <.workflow> --manifest <m> --workdir <out>
           then adapt combined_protein.tsv -> DIA-NN-shaped report

On HIVE the diann/sage command from tools.json is already Apptainer-wrapped.
Pass --sbatch to EMIT an sbatch script instead of running inline (so heavy
compute never lands on a login node) — the orchestrator submits it.

Usage:
  python3 run_search.py --tools tools.json --bundle wf/workflow.manifest.json \
      --params wf/diann.cfg --fasta search.fasta --out search_out \
      --files /data/*.raw --threads 16 [--engine diann|sage|fragpipe] [--sbatch job.sh]
"""
import sys, os, json, glob, shlex, argparse, subprocess, shutil

DIANN_CONTRACT = ["Run", "Protein.Group", "PG.MaxLFQ",
                  "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value"]


def sh(cmd, **kw):
    print(f"  $ {cmd}", flush=True)
    return subprocess.run(cmd, shell=True, check=True, **kw)


def expand_files(patterns):
    files = []
    for p in patterns:
        hits = sorted(glob.glob(p))
        files.extend(hits or [p])
    if not files:
        sys.exit("No input files matched.")
    return files


def pick_engine(args, bundle):
    if args.engine:
        return args.engine
    name = (bundle.get("engine", {}) or {}).get("name")
    if name:
        return name
    return "diann" if bundle.get("acquisition", "").upper() == "DIA" else "sage"


# ----------------------------------------------------------------- DIA-NN -----
def run_diann(cmd, params, files, fasta, out, threads, sbatch):
    os.makedirs(out, exist_ok=True)
    report = os.path.join(out, "report.parquet")
    f_args = " ".join(f"--f {shlex.quote(f)}" for f in files)
    full = (f"{cmd} --cfg {shlex.quote(params)} {f_args} "
            f"--fasta {shlex.quote(fasta)} --out {shlex.quote(report)} "
            f"--threads {threads}")
    if sbatch:
        emit_sbatch(sbatch, full, out, threads, job="diann_search")
        return {"engine": "diann", "report": report, "submitted": sbatch, "ran": False}
    sh(full)
    if not os.path.exists(report):
        sys.exit(f"DIA-NN finished but {report} is missing.")
    return {"engine": "diann", "report": report, "ran": True}


# ------------------------------------------------------------------- Sage -----
def ensure_mzml(files, out):
    """Sage is mzML-first. Convert .d/.raw via msconvert if present."""
    msconvert = shutil.which("msconvert")
    converted, need = [], []
    for f in files:
        low = f.lower()
        if low.endswith((".mzml", ".mzml.gz")):
            converted.append(f)
        else:
            need.append(f)
    if need and not msconvert:
        sys.exit("Sage needs mzML. Found non-mzML inputs but no msconvert on PATH.\n"
                 "  Convert .d/.raw to mzML first (ProteoWizard), or use a Bruker-reader Sage build.\n"
                 f"  Inputs needing conversion: {need}")
    mzdir = os.path.join(out, "mzml")
    if need:
        os.makedirs(mzdir, exist_ok=True)
        for f in need:
            sh(f"{shlex.quote(msconvert)} {shlex.quote(f)} --mzML --zlib -o {shlex.quote(mzdir)}")
            base = os.path.splitext(os.path.basename(f.rstrip('/')))[0]
            converted.append(os.path.join(mzdir, base + ".mzML"))
    return converted


def run_sage(cmd, params, files, fasta, out, threads, sbatch):
    os.makedirs(out, exist_ok=True)
    mzml = ensure_mzml(files, out)
    files_args = " ".join(shlex.quote(m) for m in mzml)
    full = (f"{cmd} {shlex.quote(params)} -f {shlex.quote(fasta)} -o {shlex.quote(out)} "
            f"--parquet --disable-telemetry-i-dont-want-to-improve-sage {files_args}")
    if sbatch:
        emit_sbatch(sbatch, full, out, threads, job="sage_search")
        return {"engine": "sage", "out": out, "submitted": sbatch, "ran": False,
                "note": "After the job runs, re-run with --adapt-only to build report.parquet."}
    sh(full)
    report = adapt_sage(out)
    return {"engine": "sage", "report": report, "ran": True}


def adapt_sage(out):
    """Map Sage lfq.parquet -> a DIA-NN-shaped protein x run report.parquet.

    This adapter is the part flagged for real-data testing (Sage VALIDATION.md).
    Sage's lfq.parquet has, per (protein, filename), an LFQ intensity. We emit
    the minimal DIA-NN contract columns the MaxLFQ DE path needs.
    """
    try:
        import pyarrow.parquet as pq
        import pyarrow as pa
    except ImportError:
        sys.exit("pyarrow required to adapt Sage output. pip install pyarrow.")

    lfq = _find(out, ["lfq.parquet"])
    if not lfq:
        sys.exit(f"No lfq.parquet under {out}; was Sage run with quant.lfq=true?")
    t = pq.read_table(lfq)
    cols = {c.lower(): c for c in t.column_names}

    def col(*cands):
        for c in cands:
            if c.lower() in cols:
                return cols[c.lower()]
        return None

    c_prot = col("proteins", "protein", "protein_group")
    c_run = col("filename", "run", "file")
    c_int = col("intensity", "lfq", "abundance")
    if not all([c_prot, c_run, c_int]):
        sys.exit(f"Sage lfq.parquet missing expected columns; saw {t.column_names}")

    prot = t.column(c_prot).to_pylist()
    run = [os.path.splitext(os.path.basename(str(r)))[0] for r in t.column(c_run).to_pylist()]
    inten = t.column(c_int).to_pylist()

    n = len(prot)
    out_tbl = pa.table({
        "Run": run,
        "Protein.Group": [str(p) for p in prot],
        "PG.MaxLFQ": [float(x) if x is not None else float("nan") for x in inten],
        "Q.Value": [0.0] * n,            # Sage already FDR-filtered at write time
        "Lib.Q.Value": [0.0] * n,
        "Lib.PG.Q.Value": [0.0] * n,
    })
    report = os.path.join(out, "report.parquet")
    pq.write_table(out_tbl, report)
    print(f"  [adapt] Sage -> {report}  ({n} protein×run rows)")
    return report


# --------------------------------------------------------------- FragPipe -----
def run_fragpipe(cmd, bundle, params, files, fasta, out, threads, sbatch):
    os.makedirs(out, exist_ok=True)
    acq = bundle.get("acquisition", "DDA").upper()
    manifest = os.path.join(out, "fragpipe.fp-manifest")
    dtype = "DIA" if acq == "DIA" else "DDA"
    with open(manifest, "w") as fh:
        for f in files:
            fh.write(f"{os.path.abspath(f)}\t\t\t{dtype}\n")
    tools = os.environ.get("FRAGPIPE_TOOLS_FOLDER", "")
    tools_arg = f"--config-tools-folder {shlex.quote(tools)}" if tools else ""
    full = (f"{cmd} --headless --workflow {shlex.quote(params)} "
            f"--manifest {shlex.quote(manifest)} --workdir {shlex.quote(out)} {tools_arg}")
    if sbatch:
        emit_sbatch(sbatch, full, out, threads, job="fragpipe_search")
        return {"engine": "fragpipe", "out": out, "submitted": sbatch, "ran": False}
    sh(full)
    report = adapt_fragpipe(out)
    return {"engine": "fragpipe", "report": report, "ran": True}


def adapt_fragpipe(out):
    """combined_protein.tsv (IonQuant MaxLFQ) -> DIA-NN-shaped report.parquet."""
    try:
        import pyarrow as pa, pyarrow.parquet as pq
    except ImportError:
        sys.exit("pyarrow required to adapt FragPipe output.")
    import csv
    cp = _find(out, ["combined_protein.tsv"])
    if not cp:
        sys.exit(f"No combined_protein.tsv under {out}.")
    with open(cp, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit("combined_protein.tsv is empty.")
    # MaxLFQ intensity columns look like "<sample> MaxLFQ Intensity"
    lfq_cols = [c for c in rows[0] if c.endswith("MaxLFQ Intensity")]
    if not lfq_cols:
        lfq_cols = [c for c in rows[0] if c.endswith("Intensity") and c != "Intensity"]
    if not lfq_cols:
        sys.exit("No per-sample MaxLFQ Intensity columns in combined_protein.tsv.")
    pid_col = "Protein" if "Protein" in rows[0] else "Protein ID"
    runs, prots, ints = [], [], []
    for r in rows:
        pg = r.get(pid_col, "").strip()
        if not pg:
            continue
        for c in lfq_cols:
            sample = c.replace(" MaxLFQ Intensity", "").replace(" Intensity", "")
            val = r.get(c, "")
            try:
                v = float(val)
            except ValueError:
                v = float("nan")
            runs.append(sample); prots.append(pg); ints.append(v if v > 0 else float("nan"))
    n = len(prots)
    tbl = pa.table({"Run": runs, "Protein.Group": prots, "PG.MaxLFQ": ints,
                    "Q.Value": [0.0]*n, "Lib.Q.Value": [0.0]*n, "Lib.PG.Q.Value": [0.0]*n})
    report = os.path.join(out, "report.parquet")
    pq.write_table(tbl, report)
    print(f"  [adapt] FragPipe -> {report}  ({n} protein×run rows)")
    return report


# ------------------------------------------------------------------ helpers ---
def _find(root, names):
    for dp, _, fns in os.walk(root):
        for fn in fns:
            if fn in names:
                return os.path.join(dp, fn)
    return None


def emit_sbatch(path, command, out, threads, job):
    """Emit a minimal SLURM script (login-node-safe). Orchestrator submits it.
    Mirrors DE-LIMP queue policy: genome-center-grp/high, publicgrp/low fallback."""
    script = f"""#!/bin/bash
#SBATCH --job-name={job}
#SBATCH --output={os.path.join(out, job)}_%j.log
#SBATCH --cpus-per-task={threads}
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=high
#SBATCH --qos=genome-center-grp-high-qos
set -euo pipefail
cd {shlex.quote(os.path.abspath(out))}
{command}
"""
    with open(path, "w") as fh:
        fh.write(script)
    print(f"  [sbatch] wrote {path} — submit with: sbatch {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tools", required=True, help="tools.json from acquire_tools.sh")
    ap.add_argument("--bundle", required=True, help="workflow.manifest.json from fetch_workflows pull")
    ap.add_argument("--params", required=True, help="engine params file (diann.cfg / sage_config.json / .workflow)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", default="search_out")
    ap.add_argument("--files", nargs="+", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--engine", choices=["diann", "sage", "fragpipe"])
    ap.add_argument("--sbatch", help="emit an sbatch script at this path instead of running inline")
    ap.add_argument("--adapt-only", action="store_true",
                    help="skip the search; just build report.parquet from an existing engine output dir")
    a = ap.parse_args()

    tools = json.load(open(a.tools))
    bundle = json.load(open(a.bundle))
    engine = pick_engine(a, bundle)
    files = expand_files(a.files)

    if a.adapt_only:
        report = {"sage": adapt_sage, "fragpipe": adapt_fragpipe}.get(engine, lambda o: None)(a.out)
        print(json.dumps({"engine": engine, "report": report, "ran": False, "adapt_only": True}, indent=2))
        return

    cmd = tools.get(engine)
    if not cmd:
        sys.exit(f"tools.json has no command for engine '{engine}'. "
                 f"Re-run acquire_tools.sh, or check its notes:\n  "
                 + "\n  ".join(tools.get("notes", [])))

    print(f"[run_search] engine={engine}  files={len(files)}  threads={a.threads}  "
          f"{'(emit sbatch)' if a.sbatch else '(inline)'}")
    if engine == "diann":
        res = run_diann(cmd, a.params, files, a.fasta, a.out, a.threads, a.sbatch)
    elif engine == "sage":
        res = run_sage(cmd, a.params, files, a.fasta, a.out, a.threads, a.sbatch)
    elif engine == "fragpipe":
        res = run_fragpipe(cmd, bundle, a.params, files, a.fasta, a.out, a.threads, a.sbatch)
    else:
        sys.exit(f"unknown engine {engine}")

    # always record what was run (engine + version + exact command) for reproducibility
    try:
        os.makedirs(a.out, exist_ok=True)
        version = (tools.get("versions", {}) or {}).get(engine) \
            or (bundle.get("engine", {}) or {}).get("version")
        with open(os.path.join(a.out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": engine, "version": version, "resolved_command": cmd,
                       "params_file": a.params, "fasta": a.fasta, "threads": a.threads,
                       "n_files": len(files), "files": files,
                       "submitted_sbatch": a.sbatch or None, "result": res}, fh, indent=2)
    except Exception as e:
        sys.stderr.write(f"[run_search] could not write search_provenance.json: {e}\n")

    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
