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
           adds --dda for DDA acquisition (DIA-NN 2.6+); if inputs are Thermo .raw,
           auto-provisions a .NET 8 runtime (ensure_dotnet8.sh) so 2.6 can read them
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
      --files /data/*.raw --threads 16 [--sbatch job.sh]
      [--engine diann|alphadia|sage|fragpipe|radiant]
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
def dotnet_env_for(files):
    """DIA-NN 2.6's NATIVE binary needs a .NET 8 runtime (>= 8.0.17) to read Thermo
    .raw. If any input is .raw, resolve/install one via ensure_dotnet8.sh and return
    a shell snippet ('export DOTNET_ROOT=...; export PATH=...;') to prepend to the
    DIA-NN command (inline) or the sbatch (compute node reads the shared install).
    Returns "" for mzML/.d-only inputs (no .NET needed). Runs the helper where
    run_search.py runs -- i.e. on the HIVE login node in the hive_remote model."""
    if not any(f.lower().endswith(".raw") for f in files):
        return ""
    helper = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ensure_dotnet8.sh")
    try:
        root = subprocess.check_output(["bash", helper], text=True).strip().splitlines()[-1]
    except Exception as e:
        sys.stderr.write(
            f"[run_diann] ensure_dotnet8.sh could not provide a .NET 8 runtime ({e}).\n"
            "  DIA-NN 2.6 will not read .raw. Options: run ensure_dotnet8.sh on a login\n"
            "  node (needs internet), or feed mzML instead of .raw.\n")
        return ""
    return (f"export DOTNET_ROOT={shlex.quote(root)}; "
            f'export PATH={shlex.quote(root)}:"$PATH";')


# Per-user CPU cap on HIVE's genome-center-grp/high (docs/QUEUE_SWITCHING.md).
HIVE_USER_CPU_CAP = int(os.environ.get('HIVE_USER_CPU_CAP', '64'))


def slurm_available():
    return shutil.which("sbatch") is not None


def _diann_parallel_mod():
    """Import the sibling generator so the parallel-safety rule has ONE definition."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import diann_parallel
    return diann_parallel


def parallel_decision(engine, files, params, a):
    """Should this DIA-NN run use the 5-step SLURM chain instead of one job?

    Automatic above --parallel-threshold files (default 5, matching the facility's
    DE-LIMP practice), but only where the chain is actually valid: it needs a cluster
    and it needs pinned mass accuracy. When a precondition fails we fall back to the
    single-shot search rather than erroring -- routing was our choice, not the user's.
    Returns (use_parallel, reason)."""
    n = len(files)
    if engine != "diann":
        return False, f"engine is {engine}; the chain is DIA-NN only"
    if a.no_parallel:
        return False, "--no-parallel was given"
    if n <= a.parallel_threshold:
        return False, f"{n} file(s), at or below the threshold of {a.parallel_threshold}"
    if not slurm_available():
        return False, (f"{n} files, but no SLURM here (sbatch not on PATH) -- the 5-step "
                       "chain runs as job arrays, so a single search is the only option")
    ma = _diann_parallel_mod().mass_acc_status(params)
    if not ma["fixed"]:
        return False, (f"{n} files, but mass accuracy is not pinned in {params} "
                       f"({ma['reason']}); steps 3/5 reuse .quant files, so the chain needs "
                       "fixed --mass-acc/--mass-acc-ms1. Re-run estimate_params.py with the "
                       "real instrument to enable parallel")
    return True, (f"{n} files > {a.parallel_threshold}, SLURM present, mass accuracy "
                  f"{ma['reason']}")


def run_diann_parallel(cmd, params, files, fasta, out, threads, a):
    """Generate DIA-NN's 5-step SLURM chain (does not submit -- the orchestrator does,
    then watches the step-5 job with watch_run.sh)."""
    os.makedirs(out, exist_ok=True)
    listing = os.path.join(out, "parallel_input_files.txt")
    with open(listing, "w") as fh:                    # a list file survives spaces in paths
        fh.write("\n".join(files) + "\n")
    argv = [sys.executable,
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "diann_parallel.py"),
            "--diann", cmd, "--raw-list", listing, "--fasta", fasta,
            "--out", out, "--cfg", params, "--threads-per-file", str(threads)]
    for flag, val in (("--partition", a.partition), ("--account", a.account),
                      ("--qos", a.qos), ("--max-simultaneous", a.max_simultaneous)):
        if val:
            argv += [flag, str(val)]
    res = subprocess.run(argv, capture_output=True, text=True)
    if res.stderr:
        sys.stderr.write(res.stderr)
    if res.returncode != 0:
        sys.exit(f"diann_parallel.py failed (exit {res.returncode}). "
                 "Re-run with --no-parallel to fall back to a single search.")
    info = json.loads(res.stdout)
    info.update({"engine": "diann", "mode": "parallel_5step", "ran": False})
    return info


a_globals = None   # set in main(); carries --one-step


def run_diann(cmd, params, files, fasta, out, threads, sbatch, acquisition=""):
    os.makedirs(out, exist_ok=True)
    report = os.path.join(out, "report.parquet")
    f_args = " ".join(f"--f {shlex.quote(f)}" for f in files)
    # DIA-NN 2.6 supports DDA via --dda (must NOT be used on DIA data). QuantUMS is
    # auto-disabled on DDA; for DDA quant DIA-NN recommends extra MS1 filtering on
    # Ms1.Global.Q.Value / Ms1.Global.Quality (see references/search-engines.md).
    dda = " --dda" if (acquisition or "").upper() == "DDA" else ""

    # Library-free runs are split into TWO SLURM JOBS: predict the library, then search
    # against it as an afterok dependency. Not because DIA-NN's one-step warning is
    # fatal -- its author says that warning is benign -- but because it is the right
    # shape regardless: the prediction is a single-threaded-ish CPU job with different
    # resource needs from the search, it is expensive to redo, and a separate job means
    # a failed search can be requeued against the SAME library instead of rebuilding it.
    # The 5-step parallel chain already worked this way; this makes the single-shot path
    # match. `--one-step` collapses them back if you ever want the old behaviour.
    cfg_txt = ""
    try:
        cfg_txt = open(params).read()
    except OSError:
        pass
    libfree = all(f in cfg_txt for f in ("--fasta-search", "--gen-spec-lib")) \
        and not getattr(a_globals, "one_step", False)
    dnet = dotnet_env_for(files)          # .NET 8 for reading Thermo .raw, if needed

    lib = os.path.join(out, "diann_lib")
    strip = ("--fasta-search", "--gen-spec-lib", "--predictor", "--reanalyse",
             "--matrices", "--rt-profiling")
    search_cfg = " ".join(t for t in cfg_txt.split() if t not in strip)
    lib_cmd = (f"{cmd} --cfg {shlex.quote(params)} --fasta {shlex.quote(fasta)} "
               f"--out-lib {shlex.quote(lib)} --threads {threads}")
    search_cmd = (f"{cmd} {search_cfg} {f_args} --fasta {shlex.quote(fasta)} "
                  f"--lib {shlex.quote(lib)}.predicted.speclib --reanalyse --matrices "
                  f"--out {shlex.quote(report)} --threads {threads}{dda}")
    onecmd = (f"{cmd} --cfg {shlex.quote(params)} {f_args} "
              f"--fasta {shlex.quote(fasta)} --out {shlex.quote(report)} "
              f"--threads {threads}{dda}")

    # TWO JOBS + dependency when emitting sbatch: the library is expensive and
    # reusable, so a failed search requeues against it instead of rebuilding.
    if libfree and sbatch:
        lib_sh = sbatch.replace(".sh", "") + "_1_lib.sh"
        srch_sh = sbatch.replace(".sh", "") + "_2_search.sh"
        emit_sbatch(lib_sh, lib_cmd, out, threads, job="diann_libpred", preamble=dnet)
        emit_sbatch(srch_sh, search_cmd, out, threads, job="diann_search", preamble=dnet)
        submit = os.path.join(out, "submit.sh")
        with open(submit, "w") as fh:
            fh.write("#!/bin/bash -l\nset -euo pipefail\n"
                     f"j1=$(sbatch --parsable {shlex.quote(lib_sh)})\n"
                     f"j2=$(sbatch --parsable --dependency=afterok:$j1 {shlex.quote(srch_sh)})\n"
                     'echo "submitted: libpred=$j1 search=$j2"\n'
                     f'echo "report will be {report}"\n')
        os.chmod(submit, 0o755)
        print(f"  [sbatch] two-job chain: {lib_sh} -> {srch_sh}; submit with: bash {submit}")
        return {"engine": "diann", "report": report, "submitted": submit,
                "mode": "two_job_libfree", "library": lib + ".predicted.speclib",
                "ran": False, "dda": bool(dda), "raw_dotnet": bool(dnet)}

    full = (lib_cmd + " && " + search_cmd) if libfree else onecmd
    if sbatch:
        emit_sbatch(sbatch, full, out, threads, job="diann_search", preamble=dnet)
        return {"engine": "diann", "report": report, "submitted": sbatch,
                "ran": False, "dda": bool(dda), "raw_dotnet": bool(dnet)}
    sh((dnet + " " if dnet else "") + full)
    if not os.path.exists(report):
        sys.exit(f"DIA-NN finished but {report} is missing.")
    return {"engine": "diann", "report": report, "ran": True, "dda": bool(dda)}


# --------------------------------------------------------------- AlphaDIA -----
# Apache-2.0 (commercial use OK) — the open-source DIA alternative to DIA-NN,
# whose free "Academia" build is academic/non-profit only. Library-free:
#   alphadia -o <out> -f <raw> [-f ...] --fasta <fasta> [-c <config.yaml>]
def run_alphadia(cmd, config, files, fasta, out, threads, sbatch):
    os.makedirs(out, exist_ok=True)
    f_args = " ".join(f"-f {shlex.quote(f)}" for f in files)
    cfg = f"-c {shlex.quote(config)} " if config and os.path.exists(config) else ""
    full = (f"{cmd} -o {shlex.quote(out)} {f_args} --fasta {shlex.quote(fasta)} {cfg}").strip()
    if sbatch:
        emit_sbatch(sbatch, full, out, threads, job="alphadia_search")
        return {"engine": "alphadia", "out": out, "submitted": sbatch, "ran": False,
                "note": "After the job runs, re-run with --adapt-only to build report.parquet."}
    sh(full)
    return {"engine": "alphadia", "report": adapt_alphadia(out), "ran": True}


def adapt_alphadia(out):
    """AlphaDIA pg.matrix.parquet (protein-group × run) -> DIA-NN-shaped report.parquet.
    Falls back to precursors.parquet (raw.name, pg.name, pg.intensity). Like the Sage
    adapter, this is the part to confirm on real data the first time."""
    try:
        import pyarrow.parquet as pq, pyarrow as pa
    except ImportError:
        sys.exit("pyarrow required to adapt AlphaDIA output. pip install pyarrow.")

    runs, prots, ints = [], [], []
    pgm = _find(out, ["pg.matrix.parquet"])
    if pgm:
        t = pq.read_table(pgm); cols = t.column_names
        id_col = next((c for c in cols if c.lower() in
                       ("pg", "pg.name", "protein", "proteins", "protein.group", "proteingroup")), cols[0])
        sample_cols = [c for c in cols if c != id_col]
        ids = [str(x) for x in t.column(id_col).to_pylist()]
        for sc in sample_cols:
            rn = os.path.splitext(os.path.basename(str(sc)))[0]
            for pid, v in zip(ids, t.column(sc).to_pylist()):
                runs.append(rn); prots.append(pid)
                ints.append(float(v) if v not in (None, 0) else float("nan"))
    else:
        pr = _find(out, ["precursors.parquet"])
        if not pr:
            sys.exit(f"No pg.matrix.parquet or precursors.parquet under {out}.")
        t = pq.read_table(pr); cols = {c.lower(): c for c in t.column_names}
        def col(*c):
            for x in c:
                if x.lower() in cols: return cols[x.lower()]
            return None
        c_run, c_pg, c_int = col("raw.name", "run"), col("pg.name", "pg", "protein.group"), col("pg.intensity", "intensity")
        if not all([c_run, c_pg, c_int]):
            sys.exit(f"AlphaDIA precursors.parquet missing expected columns; saw {t.column_names}")
        best = {}
        for r, p, v in zip(t.column(c_run).to_pylist(), t.column(c_pg).to_pylist(), t.column(c_int).to_pylist()):
            if v is None: continue
            rn = os.path.splitext(os.path.basename(str(r)))[0]
            best[(rn, str(p))] = max(best.get((rn, str(p)), 0.0), float(v))
        for (rn, p), v in best.items():
            runs.append(rn); prots.append(p); ints.append(v)

    n = len(prots)
    report = os.path.join(out, "report.parquet")
    pq.write_table(pa.table({
        "Run": runs, "Protein.Group": prots, "PG.MaxLFQ": ints,
        "Q.Value": [0.0]*n, "Lib.Q.Value": [0.0]*n, "Lib.PG.Q.Value": [0.0]*n,
    }), report)
    print(f"  [adapt] AlphaDIA -> {report}  ({n} protein×run rows)")
    return report


# -------------------------------------------------- Radiant DIA + Fulcrum -----
# Seer ships Radiant only as a container, and the container CLI reads mzML or
# Parquet -- NOT Bruker .d and NOT Thermo .raw (verified against
# radiant_fulcrum_search/search.py in the 2.3.3 image). So this route is scoped to
# Thermo Orbitrap DIA with a .raw -> mzML conversion in front of it.
#
# Radiant also always needs a spectral library: the `--libfree` flag is a MISNOMER
# -- in the click definition it is the same switch as `--no-mbr`
# ("--mbr/--no-mbr", "--no-libfree/--libfree"), so it selects single-pass vs
# match-between-runs and has nothing to do with running without a library.
# `--library` is `required=True` either way. We generate that library with DIA-NN's
# predictor, because Radiant's TSV reader takes DIA-NN's library schema directly
# (FragLibTsvReader test fixture header == DIA-NN report-lib.tsv columns).

def _container_argv(tools, mounts, inner):
    """Build a docker/apptainer invocation, binding each host dir given in `mounts`.

    mounts: list of (host_dir, container_dir). Docker and Apptainer spell bind
    mounts differently, which is why acquire_tools.sh records the runtime.
    """
    prefix = tools.get("radiant")
    runtime = tools.get("radiant_runtime")
    image = tools.get("radiant_image")
    if not (prefix and runtime and image):
        sys.exit("tools.json has no usable Radiant runtime/image. Re-run:\n"
                 "  ACQUIRE_RADIANT=1 PIN_ENGINE=radiant PIN_VERSION=<ver> "
                 "bash scripts/acquire_tools.sh <platform_class>")
    flag = "-v" if runtime == "docker" else "--bind"
    argv = prefix.split()
    for host, cont in mounts:
        argv += [flag, f"{os.path.abspath(host)}:{cont}"]
    argv += [image] + inner
    return argv


# Radiant's library loader accepts exactly these four, dispatching on the suffix
# (FragLibReader.cpp: FRAG_LIB_FF / TSV / CSV / SPEC_LIB suffixes). DIA-NN's own
# .predicted.speclib therefore needs NO conversion -- it ends in .speclib.
RADIANT_LIB_SUFFIXES = (".fraglibff", ".tsv", ".csv", ".speclib")


def _find_diann_library(stem_dir, stem):
    """Locate whatever DIA-NN actually wrote, in Radiant-preference order.

    DIA-NN IGNORES the extension you give --out-lib for a PREDICTED library: it always
    writes <stem>.predicted.speclib, its compact binary format (DIA-NN docs, "Output
    library": "For predicted library generation, however, the output file takes the
    .predicted.speclib extension"). Asking for .tsv and then checking for .tsv fails
    even though the run succeeded -- so look for what it really produces.
    """
    for cand in (f"{stem}.predicted.speclib", f"{stem}.speclib",
                 f"{stem}.parquet", f"{stem}.tsv"):
        p = os.path.join(stem_dir, cand)
        if os.path.exists(p) and os.path.getsize(p) > 1000:
            return p
    return None


def radiant_library(tools, fasta, out, threads, params=None):
    """Build the spectral library Radiant searches against, from DIA-NN's predictor.

    NOT a straight hand-off: DIA-NN 2.x writes .predicted.speclib v-10/-11 and Radiant's
    reader only supports v>=-3, so its binary is rejected outright. The conversion lives
    in make_radiant_library.py.
    """
    libdir = os.path.join(out, "radiant_lib")
    ready = os.path.join(libdir, "radiant_library.tsv")
    if os.path.exists(ready) and os.path.getsize(ready) > 1_000_000:
        print(f"  [radiant] reusing existing library {ready}")
        return ready
    dn = tools.get("diann")
    if not dn:
        sys.exit("Radiant needs a spectral library, but tools.json has no DIA-NN command. "
                 "Acquire DIA-NN first (it is the library generator for this route), or "
                 "pass --library with an existing .tsv/.csv/.fragLibFF library.")
    os.makedirs(libdir, exist_ok=True)
    # DIA-NN 2.x and Radiant share NO library format directly: DIA-NN writes
    # .predicted.speclib v-10/-11, Radiant's reader supports v>=-3 ("ERROR: version is
    # not supported11"), and DIA-NN cannot emit TSV. make_radiant_library.py does the
    # required predict -> parquet -> renamed TSV hop. See its docstring for the evidence.
    helper = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "make_radiant_library.py")
    res = subprocess.run([sys.executable, helper, "--diann", dn, "--fasta", fasta,
                          "--out-dir", libdir, "--threads", str(threads)],
                         capture_output=True, text=True)
    sys.stderr.write(res.stderr)
    if res.returncode != 0:
        sys.exit(f"make_radiant_library.py failed (exit {res.returncode}).")
    info = json.loads(res.stdout[res.stdout.index("{"):])
    print(f"  [radiant] library -> {info['library']} ({info.get('rows', '?')} rows)")
    return info["library"]


def run_radiant_parallel(tools, params, files, fasta, out, threads, a, library=None):
    """On a cluster, search each file as its own array task, then rescore once.

    Radiant's Fulcrum backend is SERIAL (`NotImplementedError` for parallel mode), so
    an N-file study otherwise costs N x one-file wall-clock even on a big node. The
    search is per-file independent though; only the downstream rescoring/FDR/rollup
    needs the whole set. radiant_parallel.py splits exactly there and emits a 3-step
    chain. Emits scripts; the orchestrator submits them with dependencies."""
    os.makedirs(out, exist_ok=True)
    listing = os.path.join(out, "radiant_input_files.txt")
    with open(listing, "w") as fh:
        fh.write("\n".join(files) + "\n")
    argv = [sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "radiant_parallel.py"),
            "--runtime", tools.get("radiant_runtime") or "docker",
            "--image", tools.get("radiant_image") or "",
            "--raw-list", listing, "--fasta", fasta, "--config", params,
            "--out", out, "--threads-per-file", str(threads)]
    if library:
        argv += ["--library", library]
    elif tools.get("diann"):
        argv += ["--diann", tools["diann"]]
    for flag, val in (("--partition", a.partition), ("--account", a.account),
                      ("--qos", a.qos), ("--max-simultaneous", a.max_simultaneous)):
        if val:
            argv += [flag, str(val)]
    if getattr(a, "mbr", False):
        argv.append("--mbr")
    res = subprocess.run(argv, capture_output=True, text=True)
    if res.stderr:
        sys.stderr.write(res.stderr)
    if res.returncode != 0:
        sys.exit(f"radiant_parallel.py failed (exit {res.returncode}). "
                 "Re-run with --no-parallel for a single serial search.")
    info = json.loads(res.stdout)
    info.update({"engine": "radiant", "ran": False})
    return info


def run_radiant(tools, params, files, fasta, out, threads, sbatch, library=None, mbr=True):
    os.makedirs(out, exist_ok=True)
    bad = [f for f in files if f.rstrip("/").lower().endswith(".d")]
    if bad:
        sys.exit("Radiant/Fulcrum does not read Bruker .d in this container — it takes "
                 "mzML or Parquet. This route is for Thermo Orbitrap DIA.\n"
                 f"  Bruker inputs: {bad}\n"
                 "  Use the DIA-NN or FragPipe/diaTracer route for timsTOF data.")
    mzml = ensure_mzml(files, out)          # .raw -> mzML (msconvert), same path Sage uses
    if library and not library.lower().endswith(RADIANT_LIB_SUFFIXES):
        sys.exit(f"Radiant cannot read {os.path.basename(library)} — its loader accepts "
                 f"only {', '.join(RADIANT_LIB_SUFFIXES)}.")
    if library and library.lower().endswith(".speclib"):
        sys.stderr.write(
            "[run_search] WARNING: Radiant only reads .speclib format v>=-3, and every\n"
            "  DIA-NN 2.x library is v-10/-11 — it will abort with 'version is not\n"
            "  supported'. If this is a DIA-NN library, convert it first:\n"
            "    python3 scripts/make_radiant_library.py --from-speclib <lib> "
            "--diann '<cmd>' --out-dir <dir>\n")
    lib = library or radiant_library(tools, fasta, out, threads, params)

    results = os.path.join(out, "radiant_results")
    os.makedirs(results, exist_ok=True)

    # Bind each distinct host directory the container must see. Mounting parents
    # (rather than copying) keeps large mzML in place.
    mounts, cmap = [], {}
    def cpath(host, tag):
        d = os.path.dirname(os.path.abspath(host))
        if d not in cmap:
            cmap[d] = f"/mnt/{tag}{len(cmap)}"
            mounts.append((d, cmap[d]))
        return f"{cmap[d]}/{os.path.basename(host)}"

    c_lib, c_fa = cpath(lib, "in"), cpath(fasta, "in")
    c_files = [cpath(f, "in") for f in mzml]
    c_res = "/mnt/results"
    mounts.append((results, c_res))

    inner = ["radiant_fulcrum", "-v",
             "--mbr" if mbr else "--libfree",
             "--library", c_lib, "--fasta", c_fa,
             "--results-dir", c_res, "--threads", str(threads)]
    if params and params.lower().endswith((".radiantconfig", ".toml", ".pythiaconfig")):
        inner += ["--config", cpath(params, "in")]
    inner += c_files

    argv = _container_argv(tools, mounts, inner)
    full = " ".join(shlex.quote(x) for x in argv)
    if sbatch:
        emit_sbatch(sbatch, full, out, threads, job="radiant_search")
        return {"engine": "radiant", "out": out, "submitted": sbatch, "ran": False}
    sh(full)
    report = adapt_radiant(out)
    return {"engine": "radiant", "report": report, "ran": True, "library": lib}


def _radiant_run_name(raw):
    """Fulcrum reports Run as a full URI of the per-file result, e.g.
    `file:///mnt/results/radiant-results/Sample_1.mzML.radiantDIA`.

    A single splitext leaves `Sample_1.mzML`, which does NOT match the bare run names
    every other engine emits — so a conditions.csv keyed on sample names would silently
    fail to assign groups. Strip the URI, the container path, and BOTH extensions.
    """
    s = str(raw)
    if "://" in s:
        s = s.split("://", 1)[1]
    s = os.path.basename(s.rstrip("/"))
    for _ in range(3):                       # .mzML.radiantDIA, .d.radiantDIA, ...
        stem, ext = os.path.splitext(s)
        if ext.lower() in (".radiantdia", ".mzml", ".raw", ".d", ".gz", ".parquet"):
            s = stem
        else:
            break
    return s


def adapt_radiant(out):
    """Fulcrum `combined` output -> the report.parquet contract.

    Fulcrum's combined backend already emits DIA-NN-style column names (Run,
    Protein.Group, PG.Quantity/PG.Normalised, Q.Value, Global.PG.Q.Value ...), but
    writes them as a SPARK PARQUET DIRECTORY of part-files, not a single file --
    so read it as a dataset.
    """
    try:
        import pyarrow.parquet as pq, pyarrow as pa, pyarrow.dataset as ds
    except ImportError:
        sys.exit("pyarrow required to adapt Radiant output. pip install pyarrow.")

    # Fulcrum writes a DIRECTORY of spark part-files, so _find() (files only) can't
    # locate it — walk for the directory name instead. A single-file parquet is
    # accepted too, in case a future backend writes one.
    root = None
    for cand in ("fulcrum-results", "fulcrum-proteins"):
        for dp, dns, fns in os.walk(out):
            if cand in dns:
                root = os.path.join(dp, cand)
                break
            if cand in fns:
                root = os.path.join(dp, cand)
                break
        if root:
            break
    if not root:
        sys.exit(f"No fulcrum-results/ under {out}. Did the Fulcrum workflow finish? "
                 f"Check the Radiant logs in {out}.")

    t = ds.dataset(root, format="parquet").to_table()
    cols = {c.lower(): c for c in t.column_names}

    def col(*cands):
        for c in cands:
            if c.lower() in cols:
                return cols[c.lower()]
        return None

    c_run = col("Run", "filename", "File.Name")
    c_pg = col("Protein.Group", "ProteinGroup", "PG")
    # Prefer the normalised protein-group quantity; that is the MaxLFQ analogue here.
    c_int = col("PG.Normalised", "PG.Quantity", "PG.MaxLFQ", "Precursor.Normalised",
                "Precursor.Quantity")
    if not all([c_run, c_pg, c_int]):
        sys.exit(f"Radiant output missing expected columns; saw {t.column_names}")
    c_q = col("Global.PG.Q.Value", "Q.Value")

    runs = [_radiant_run_name(r) for r in t.column(c_run).to_pylist()]
    pgs = [str(p) for p in t.column(c_pg).to_pylist()]
    vals = t.column(c_int).to_pylist()
    qs = t.column(c_q).to_pylist() if c_q else [0.0] * len(pgs)

    # The combined report is PSM-level (many precursors per protein x run); collapse
    # to one protein x run row, keeping the best q-value seen.
    best = {}
    for r, p, v, q in zip(runs, pgs, vals, qs):
        if v is None:
            continue
        k = (r, p)
        prev = best.get(k)
        qq = float(q) if q is not None else 0.0
        if prev is None or float(v) > prev[0]:
            best[k] = (float(v), min(qq, prev[1]) if prev else qq)
        else:
            best[k] = (prev[0], min(prev[1], qq))

    runs2 = [k[0] for k in best]
    pgs2 = [k[1] for k in best]
    ints = [v[0] for v in best.values()]
    qv = [v[1] for v in best.values()]
    n = len(pgs2)
    report = os.path.join(out, "report.parquet")
    pq.write_table(pa.table({
        "Run": runs2, "Protein.Group": pgs2, "PG.MaxLFQ": ints,
        "Q.Value": qv, "Lib.Q.Value": qv, "Lib.PG.Q.Value": qv,
    }), report)
    print(f"  [adapt] Radiant/Fulcrum -> {report}  ({n} protein×run rows)")
    return report


# ------------------------------------------------------------------- Sage -----
def ensure_mzml(files, out):
    """mzML-first engines (Sage, Radiant). Convert .d/.raw via msconvert if present."""
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

    if acq == "DIA" and any(f.rstrip("/").endswith(".d") for f in files):
        # diaTracer writes its pseudo-MS/MS mzML NEXT TO THE INPUT, so pointing it at the
        # shared raw files would have two users racing to write the same output (and would
        # fail outright on a read-only share). Stage per-user symlinks instead: FragPipe
        # normalizes but does not resolve them, so the output lands in our own directory.
        # The stager also reuses any conversion that already exists.
        stage = os.path.join(out, "diatracer_stage")
        res = subprocess.run(
            [sys.executable,
             os.path.join(os.path.dirname(os.path.abspath(__file__)), "diatracer_stage.py"),
             "--raw", *files, "--stage", stage, "--manifest", manifest],
            capture_output=True, text=True)
        if res.returncode != 0:
            sys.stderr.write(res.stderr)
            sys.exit("diatracer_stage.py failed — cannot build a safe FragPipe manifest.")
        info = json.loads(res.stdout)
        print(f"  [diatracer] staged {info['n_files']} file(s): "
              f"{info['n_to_convert']} to convert, {info['n_reused']} reused -> {stage}")
        for n in info.get("notes", []):
            print(f"  [diatracer] {n}")
    else:
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


def adapt_fragpipe_dia(out):
    """FragPipe's DIA route (diaTracer -> MSFragger -> DIA-NN) writes DIA-NN's own
    output to <workdir>/dia-quant-output/ (report.parquet + report.tsv +
    report.stats.tsv -- verified in FragPipe 24.0 CmdDiann.java). report.parquet is
    ALREADY the DE contract, so prefer it and only fall back to converting the TSV.
    Returns None if this doesn't look like a DIA run, so the caller can try DDA."""
    # Look for the file INSIDE dia-quant-output specifically. A plain _find() would also
    # match the report.parquet this function itself writes into <out> on an earlier run,
    # whose parent is <out> rather than dia-quant-output -- so the DIA branch would miss
    # and needlessly re-convert the TSV every time.
    def in_dia_out(name):
        for root, dirs, files in os.walk(out):
            if os.path.basename(root) == "dia-quant-output" and name in files:
                return os.path.join(root, name)
        return None

    pq_path = in_dia_out("report.parquet")
    if pq_path:
        print(f"  [adapt] FragPipe DIA: DIA-NN output already meets the contract -> {pq_path}")
        return pq_path
    tsv = in_dia_out("report.tsv")
    if not tsv:
        return None
    try:
        import pyarrow as pa, pyarrow.parquet as pq
    except ImportError:
        sys.exit("pyarrow required to adapt FragPipe DIA output.")
    import csv
    with open(tsv, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"{tsv} is empty — DIA-NN produced no quantification.")
    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return float("nan")
    cols = {c.lower(): c for c in rows[0]}
    def col(*names):
        for n in names:
            if n.lower() in cols:
                return cols[n.lower()]
        return None
    c_run, c_pg, c_q = col("Run"), col("Protein.Group"), col("PG.MaxLFQ", "PG.Quantity")
    if not all([c_run, c_pg, c_q]):
        sys.exit(f"{tsv} lacks Run / Protein.Group / PG.MaxLFQ — cannot build the DE contract.")
    # Carry the q-value columns through when present; the DE step filters on them.
    cq, clq, clpq = col("Q.Value"), col("Lib.Q.Value"), col("Lib.PG.Q.Value")
    n = len(rows)
    tbl = pa.table({
        "Run": [r[c_run] for r in rows],
        "Protein.Group": [r[c_pg] for r in rows],
        "PG.MaxLFQ": [num(r[c_q]) for r in rows],
        "Q.Value": [num(r[cq]) if cq else 0.0 for r in rows],
        "Lib.Q.Value": [num(r[clq]) if clq else 0.0 for r in rows],
        "Lib.PG.Q.Value": [num(r[clpq]) if clpq else 0.0 for r in rows],
    })
    report = os.path.join(out, "report.parquet")
    pq.write_table(tbl, report)
    print(f"  [adapt] FragPipe DIA: {os.path.basename(tsv)} -> {report}  ({n} rows)")
    return report


def adapt_fragpipe(out):
    """FragPipe -> DIA-NN-shaped report.parquet. Tries the DIA route first (diaTracer
    leaves DIA-NN output in dia-quant-output/), then the DDA route (IonQuant's
    combined_protein.tsv). The two produce different files, so which one exists is
    what tells us which route ran."""
    dia = adapt_fragpipe_dia(out)
    if dia:
        return dia
    return adapt_fragpipe_dda(out)


def adapt_fragpipe_dda(out):
    """combined_protein.tsv (IonQuant MaxLFQ) -> DIA-NN-shaped report.parquet."""
    try:
        import pyarrow as pa, pyarrow.parquet as pq
    except ImportError:
        sys.exit("pyarrow required to adapt FragPipe output.")
    import csv
    cp = _find(out, ["combined_protein.tsv"])
    if not cp:
        sys.exit(f"No FragPipe output found under {out}: neither dia-quant-output/report.* "
                 "(the diaTracer DIA route) nor combined_protein.tsv (the IonQuant DDA "
                 "route). Check the FragPipe log — a headless run can exit 0 on a crash.")
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



def _partition_idle_cpus(part):
    """Idle CPUs in a partition. `sinfo -h -o %C` gives Alloc/Idle/Other/Total."""
    try:
        out = subprocess.run(["sinfo", "-p", part, "-h", "-o", "%C"],
                             capture_output=True, text=True, timeout=20).stdout.strip()
        return int(out.split("/")[1]) if "/" in out else 0
    except Exception:
        return 0


def _lab_cpus_available():
    """CPUs left under MY per-user cap on the priority queue, or None if unknown.

    The per-user limit is the binding constraint on genome-center-grp/high (not the
    much larger account limit, which is shared), so this counts what I am already
    running there rather than what the group is."""
    user = os.environ.get("USER", "")
    try:
        used = 0
        out = subprocess.run(["squeue", "-h", "-u", user, "-t", "RUNNING",
                              "-p", "high", "-o", "%C"],
                             capture_output=True, text=True, timeout=20).stdout
        for ln in out.split():
            try:
                used += int(ln)
            except ValueError:
                pass
        return max(HIVE_USER_CPU_CAP - used, 0)
    except Exception:
        return None


def slurm_queue(partition=None, account=None, qos=None,
                peak_cpus=None, preemptible_ok=False):
    """Pick a SLURM partition/account/qos the CURRENT USER can actually submit to.

    Never hardcode a queue. The old behaviour emitted
    `--partition=high --qos=genome-center-grp-high-qos` unconditionally, which a user
    outside genome-center-grp cannot submit to at all — the job is REJECTED, not merely
    slowed. And on HIVE `high` caps publicgrp at 8 CPUs / 128 GB per job, so a 32-CPU
    request there would never start (QOSMaxCpuPerJobLimit).

    Ask SLURM what this account is entitled to, then prefer, in order:
      1. an explicit override the caller passed
      2. genome-center-grp on `high`   (facility members: no per-job cap, not preemptible)
      3. publicgrp on `low`            (everyone else, incl. class accounts: no per-job
                                        cap either, preemptible — add --requeue)
    Returns (partition, account, qos); any may be None, and a None is simply omitted
    from the script so SLURM applies its own default."""
    if partition and account:
        return partition, account, qos
    assoc = []
    # sacctmgr is frequently absent from PATH in a non-login shell, so look for it
    # explicitly. Failing to find it must NOT silently emit an empty queue: SLURM would
    # then use the cluster default partition, which on HIVE is `high` — precisely the
    # queue a non-facility account cannot use.
    sacctmgr = shutil.which("sacctmgr")
    if not sacctmgr:
        for c in ("/usr/bin/sacctmgr", "/usr/local/bin/sacctmgr",
                  "/cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin/sacctmgr"):
            if os.path.exists(c):
                sacctmgr = c
                break
    try:
        if not sacctmgr:
            raise FileNotFoundError("sacctmgr not found")
        out = subprocess.run(
            [sacctmgr, "-nP", "show", "assoc",
             f"user={os.environ.get('USER', '')}", "format=account,partition,qos"],
            capture_output=True, text=True, timeout=30).stdout
        for line in out.splitlines():
            f = line.split("|")
            if len(f) >= 3 and f[0]:
                assoc.append((f[0].strip(), f[1].strip(), f[2].strip()))
    except Exception:
        pass                                  # no SLURM, or sacctmgr unavailable

    def find(acct, part):
        for a, p, q in assoc:
            if a == acct and p == part:
                return a, p, (q or None)
        return None

    lab, pub = find("genome-center-grp", "high"), find("publicgrp", "low")

    # Port of DE-LIMP's select_best_partition() (R/helpers_search.R). Entitlement is
    # not the question -- UTILISATION is. The priority queue has a PER-USER CPU cap
    # (64 on HIVE), and once you are at it your own jobs queue behind each other:
    # an 18-task array on `high` starves everything else you submit (QOSGrpCpuLimit,
    # observed). publicgrp/low is preemptible but has thousands of idle CPUs, so for
    # work that is safe to preempt it starts sooner and finishes sooner.
    if lab and pub and not partition:
        need = min(peak_cpus or 16, 16)          # at least one array task's worth
        avail = _lab_cpus_available()
        idle = _partition_idle_cpus("low")
        if avail is not None and avail < need and idle >= need:
            a, p, q = pub
            print(f"[slurm_queue] priority queue at capacity ({avail} CPUs free, need "
                  f"{need}); publicgrp/low has {idle} idle -> using low (preemptible, "
                  f"--requeue is added)", file=sys.stderr)
            return p, a, q
        if preemptible_ok and idle >= need and (avail is None or avail < need * 2):
            a, p, q = pub
            print(f"[slurm_queue] preemption-safe step and low has {idle} idle CPUs "
                  f"-> using publicgrp/low for throughput", file=sys.stderr)
            return p, a, q

    for acct, part in (("genome-center-grp", "high"), ("publicgrp", "low")):
        hit = find(acct, part)
        if hit:
            a, p, q = hit
            return partition or p, account or a, qos or q
    if assoc:                                  # entitled to something unanticipated
        a, p, q = assoc[0]
        return partition or (p or None), account or a, qos or (q or None)
    # Could not detect. Do NOT fall through to the cluster default — on HIVE that is
    # `high`, which rejects non-facility accounts. publicgrp/low is submittable by
    # everyone who has any allocation at all, so it is the safe floor.
    return partition or "low", account or "publicgrp", qos


def emit_sbatch(path, command, out, threads, job, preamble="",
                partition=None, account=None, qos=None, mem="64G", hours=12):
    """Emit a minimal SLURM script (login-node-safe). Orchestrator submits it.
    The queue is DETECTED from the submitting user's own SLURM associations — see
    slurm_queue(). `preamble` runs before the command (e.g. the DOTNET_ROOT exports
    that let DIA-NN 2.6 read Thermo .raw)."""
    part, acct, q = slurm_queue(partition, account, qos)
    pre = (preamble + "\n") if preamble else ""
    lines = [
        "#!/bin/bash -l",
        f"#SBATCH --job-name={job}",
        f"#SBATCH --output={os.path.join(out, job)}_%j.log",
        f"#SBATCH --cpus-per-task={threads}",
        f"#SBATCH --mem={mem}",
        f"#SBATCH --time={hours}:00:00",
    ]
    if part:  lines.append(f"#SBATCH --partition={part}")
    if acct:  lines.append(f"#SBATCH --account={acct}")
    if q:     lines.append(f"#SBATCH --qos={q}")
    # `low` is preemptible; without --requeue a preempted search is simply lost.
    if part == "low":
        lines.append("#SBATCH --requeue")
    lines += ["set -euo pipefail", f"cd {shlex.quote(os.path.abspath(out))}", f"{pre}{command}", ""]
    script = "\n".join(lines)
    with open(path, "w") as fh:
        fh.write(script)
    print(f"  [sbatch] wrote {path} (partition={part or 'default'}, "
          f"account={acct or 'default'}) — submit with: sbatch {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tools", required=True, help="tools.json from acquire_tools.sh")
    ap.add_argument("--bundle", required=True, help="workflow.manifest.json from fetch_workflows pull")
    ap.add_argument("--params", required=True, help="engine params file (diann.cfg / sage_config.json / .workflow)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", default="search_out")
    ap.add_argument("--files", nargs="+", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--engine", choices=["diann", "alphadia", "sage", "fragpipe", "radiant"])
    ap.add_argument("--library", help="Radiant: an existing DIA-NN .tsv spectral library. "
                                      "Omit and one is generated with DIA-NN's predictor.")
    ap.add_argument("--one-step", action="store_true",
                    help="DIA-NN library-free in ONE command instead of a "
                         "library job + dependent search job")
    ap.add_argument("--allow-inline", action="store_true",
                    help="permit an inline search on a host where sbatch exists "
                         "(only inside an salloc/srun allocation)")
    # Default ON: the DIA-NN chain shares information across runs (step 3 builds an
    # empirical library from ALL files, step 4 re-searches against it). Radiant at
    # mbr=false would be the only engine searching each file in isolation, which
    # understates it in any comparison. --no-mbr restores Seer's shipped default.
    ap.add_argument("--mbr", action=argparse.BooleanOptionalAction, default=True,
                    help="Radiant: match-between-runs / two-pass (default: on, for "
                         "parity with the DIA-NN two-pass chain). --no-mbr = single-pass.")
    ap.add_argument("--sbatch", help="emit an sbatch script at this path instead of running inline")
    ap.add_argument("--parallel-threshold", type=int, default=5,
                    help="DIA-NN: use the 5-step SLURM chain above this many files (default 5)")
    ap.add_argument("--no-parallel", action="store_true",
                    help="force a single-shot DIA-NN search regardless of file count")
    ap.add_argument("--partition", help="SLURM partition for the parallel chain")
    ap.add_argument("--account", help="SLURM account for the parallel chain")
    ap.add_argument("--qos", help="SLURM QOS for the parallel chain")
    ap.add_argument("--max-simultaneous", type=int,
                    help="cap concurrent array tasks in the parallel chain")
    ap.add_argument("--adapt-only", action="store_true",
                    help="skip the search; just build report.parquet from an existing engine output dir")
    a = ap.parse_args()

    global a_globals
    a_globals = a
    tools = json.load(open(a.tools))
    bundle = json.load(open(a.bundle))
    engine = pick_engine(a, bundle)
    files = expand_files(a.files)

    # Absolutize EVERY path before it is baked into a command. emit_sbatch() writes
    # `cd <abspath(out)>` and then the command verbatim, so a caller-relative
    # --params/--fasta/--out (exactly what SKILL.md documents: `--params ./wf/x.cfg`)
    # resolves against the WRONG directory once the job runs — the search dies on a
    # missing params file, or worse silently looks at ./out/out/. Same hazard for the
    # container routes, whose bind mounts are derived from these paths.
    a.out = os.path.abspath(a.out)
    for attr in ("params", "fasta", "library"):
        v = getattr(a, attr, None)
        if v:
            setattr(a, attr, os.path.abspath(v))
    files = [os.path.abspath(f) for f in files]

    if a.adapt_only:
        report = {"sage": adapt_sage, "fragpipe": adapt_fragpipe,
                  "alphadia": adapt_alphadia,
                  "radiant": adapt_radiant}.get(engine, lambda o: None)(a.out)
        print(json.dumps({"engine": engine, "report": report, "ran": False, "adapt_only": True}, indent=2))
        return

    cmd = tools.get(engine)
    if not cmd:
        sys.exit(f"tools.json has no command for engine '{engine}'. "
                 f"Re-run acquire_tools.sh, or check its notes:\n  "
                 + "\n  ".join(tools.get("notes", [])))

    use_parallel, why = parallel_decision(engine, files, a.params, a)
    if engine == "diann":
        print(f"[run_search] parallel routing: {'YES' if use_parallel else 'no'} -- {why}")

    # HARD STOP: never start a multi-hour search inline on a cluster LOGIN NODE.
    # Golden rule #3 says every heavy step goes through the scheduler, but nothing
    # enforced it -- and the failure is silent and expensive. It bit for real: a
    # DIA-NN run whose mass accuracy was unpinned fell out of the parallel chain into
    # the single-shot path and, with no --sbatch, launched diann-linux on the HIVE
    # login node at 16 threads (twice, because an earlier invocation was still alive).
    # sbatch present + no SLURM_JOB_ID == we are on a submit host, not a compute node.
    if not use_parallel and not a.sbatch and slurm_available() \
            and not os.environ.get("SLURM_JOB_ID") and not a.allow_inline:
        sys.exit(
            "REFUSING to run the search inline: this looks like a cluster login/submit "
            "node (sbatch is on PATH and SLURM_JOB_ID is unset).\n"
            f"  engine={engine}  files={len(files)}  threads={a.threads}\n"
            f"  parallel routing declined because: {why}\n"
            "  Re-run with --sbatch <script> and submit it, e.g.:\n"
            f"    ... --sbatch ./{engine}_job.sh && sbatch ./{engine}_job.sh\n"
            "  (--allow-inline overrides this, e.g. inside an salloc/srun session.)")

    print(f"[run_search] engine={engine}  files={len(files)}  threads={a.threads}  "
          f"{'(5-step chain)' if use_parallel else '(emit sbatch)' if a.sbatch else '(inline)'}")
    if use_parallel:
        res = run_diann_parallel(cmd, a.params, files, a.fasta, a.out, a.threads, a)
    elif engine == "diann":
        res = run_diann(cmd, a.params, files, a.fasta, a.out, a.threads, a.sbatch,
                        acquisition=bundle.get("acquisition", ""))
    elif engine == "alphadia":
        res = run_alphadia(cmd, a.params, files, a.fasta, a.out, a.threads, a.sbatch)
    elif engine == "sage":
        res = run_sage(cmd, a.params, files, a.fasta, a.out, a.threads, a.sbatch)
    elif engine == "fragpipe":
        res = run_fragpipe(cmd, bundle, a.params, files, a.fasta, a.out, a.threads, a.sbatch)
    elif engine == "radiant":
        # Takes the whole tools dict: it needs the container runtime + image to build
        # bind mounts, and DIA-NN to generate the spectral library.
        # On a cluster with >1 file, split the SERIAL search into a per-file array and
        # rescore once — Radiant's own backend cannot parallelise, but the search is
        # per-file independent, so this is the difference between N x t and ~t.
        if (len(files) > 1 and slurm_available() and not a.no_parallel):
            print(f"[run_search] radiant: {len(files)} files on SLURM -> per-file array "
                  f"+ one Fulcrum rescoring job (Radiant's own search is serial)")
            res = run_radiant_parallel(tools, a.params, files, a.fasta, a.out,
                                       a.threads, a, library=a.library)
        else:
            res = run_radiant(tools, a.params, files, a.fasta, a.out, a.threads, a.sbatch,
                              library=a.library, mbr=a.mbr)
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
                       "search_mode": "parallel_5step" if use_parallel else "single_shot",
                       "parallel_routing_reason": why,
                       "submitted_sbatch": a.sbatch or None, "result": res}, fh, indent=2)
    except Exception as e:
        sys.stderr.write(f"[run_search] could not write search_provenance.json: {e}\n")

    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
