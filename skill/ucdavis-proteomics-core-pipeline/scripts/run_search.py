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

Every input Bruker .d is checked with bruker_tdf.tdf_integrity() before anything is
provisioned or submitted, and a .d that is not `ok` is REFUSED (--allow-damaged-tdf
searches it anyway). detect_acquisition.py checks the same thing in step 2, but a caller
who already has the paths -- "re-run this search", "the paths are already known" -- comes
straight here.
"""
import sys, os, json, glob, re, shlex, argparse, subprocess, shutil, stat, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# ONE definition of the q-value columns -- see diann_q_columns.py.
from diann_q_columns import FDR_REQUIRED, PROTEIN_Q_PREFERENCE
# No Bruker .d is searched before its analysis.tdf has been checked -- see bruker_tdf.py.
from bruker_tdf import STATUSES, integrity_warning, tdf_integrity

ALLOW_DAMAGED_TDF = "--allow-damaged-tdf"

# The minimum an adapted report must carry for the DE step to run. The q-columns
# are exactly the REQUIRED filter set: an adapter that emits fewer leaves limpa
# silently applying no filter for the missing ones.
DIANN_CONTRACT = ["Run", "Protein.Group", "PG.MaxLFQ"] + FDR_REQUIRED


def sh(cmd, **kw):
    print(f"  $ {cmd}", flush=True)
    return subprocess.run(cmd, shell=True, check=True, **kw)


def tdf_problems(files):
    """Every Bruker .d in `files` whose analysis.tdf is not `ok`, worst first.

    detect_acquisition.py runs the same check, but step 2 is exactly what the entry points
    SKILL.md advertises for a bare search skip: "re-run this search", "re-search with
    different parameters", "the paths are already known" all arrive here with the files in
    hand and no detection pass behind them. A .d whose frame index covers part of its
    tdf_bin is then searched with no error anywhere -- DIA-NN reads what the index points
    at, and one HIVE run's index covered 0.7% of a 2.4 GB tdf_bin. So the gate has to be on
    THIS path too, not only on the one the orchestrator is told to walk.

    Per file this is a 100-byte header read, one indexed query and a 4-byte read -- nothing
    that needs a compute node, and nothing next to a multi-hour search."""
    out = []
    for f in files:
        p = f.rstrip("/")
        if p.lower().endswith(".d") and os.path.isfile(os.path.join(p, "analysis.tdf")):
            r = tdf_integrity(p)
            if r["status"] != "ok":
                out.append((p, r))
    out.sort(key=lambda t: (STATUSES.index(t[1]["status"]), t[0]))
    return out


def refuse_damaged_tdf(files, allow):
    """Stop before the search when any input .d is not `ok`. Prints and returns the list."""
    bad = tdf_problems(files)
    if not bad:
        return bad
    lines = [f"  {p}\n      {integrity_warning(r)}" for p, r in bad]
    if allow:
        print(f"[run_search] {ALLOW_DAMAGED_TDF}: searching {len(bad)} .d whose analysis.tdf "
              f"is not `ok` anyway:\n" + "\n".join(lines), flush=True)
        return bad
    sys.exit(f"REFUSING to search: {len(bad)} of {len(files)} input file(s) are a Bruker .d "
             f"whose analysis.tdf is not `ok`. A search of one of these reports no error and "
             f"no missing data -- it silently covers only the part of the run the frame index "
             f"still points at.\n" + "\n".join(lines) +
             f"\n  Fix or drop those files, or re-run with {ALLOW_DAMAGED_TDF} to search them "
             f"as they are.")


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


# What a version number looks like: digits and dots, optionally led by `v` -- 2.6.1, v0.14.7,
# 24.0. NOTHING else is accepted as one. Naming the strings that are NOT versions cannot work:
# tools.json `versions` is written by acquire_tools.sh, and every branch of it that cannot
# determine a build is one edit away from writing some new word there. A list of known words
# ("latest", "env") passed "nightly", "dev", "n/a", "TBD", "<unknown>" and "garbage" straight
# through to search_provenance.json `version`, which fran_deposit.py hands to FRAN as the
# engine version of a deposited search. The check has to be the shape of an answer, not a list
# of known wrong answers.
_VERSION_RE = re.compile(r"v?\d+(?:\.\d+)+")

# How a build folder / image file spells each engine's name. Everything else uses the
# engine key itself (fragpipe-24.0/ is how the FragPipe zip unpacks).
_ENGINE_NAME_RE = {"diann": r"dia-?nn", "radiant": r"radiant(?:-fulcrum)?"}


def _concrete_version(v):
    """'v0.14.7' / '0.14.7' -> '0.14.7'; anything that is not shaped like a version -> None.

    The strings this has to refuse are open-ended -- "latest" and "env" are only the two
    acquire_tools.sh writes today -- so this asks whether `v` IS a version, never whether it
    is one of the known non-versions."""
    if not isinstance(v, str) or not _VERSION_RE.fullmatch(v.strip()):
        return None
    return v.strip().lstrip("vV")


def command_engine_versions(engine, tools):
    """The versions the command tools.json runs for `engine` names IN ITSELF, sorted, distinct.

      * a build folder or image file: <name>[-_]<version>, EVERY one in the string -- the rule
        acquire_tools.sh's diann_version_of() uses, e.g. build_260/diann-2.6.0/diann-linux,
        diann_2.3.0.sif, radiant-fulcrum-2.3.3.sif, fragpipe-24.0/bin/fragpipe,
        sage-v0.14.7-x86_64-unknown-linux-gnu/sage;
      * an image tag: a word ending :<version>, e.g. proteomics-pipeline/diann:2.7.0.
    tools.json keeps Radiant's image apart from its runtime prefix (`radiant_image`), so that
    is read too. A bare version DIRECTORY (sage/0.14.6/sage, diann/2.7.0/) is deliberately not
    matched: acquire_tools.sh names cache folders after the REQUEST, and a request-keyed
    folder can hold another release (it writes a note when a pinned Sage cache does)."""
    name = _ENGINE_NAME_RE.get(engine, re.escape(engine))
    found = set()
    for s in ((tools or {}).get(engine), (tools or {}).get(f"{engine}_image")):
        if not isinstance(s, str):
            continue
        # EVERY match, not just the last one: a command that names two different builds
        # (`apptainer exec .../diann_2.3.0.sif /diann-2.6.1/diann-linux`) is exactly the case
        # the `len(cvs) > 1` refusal below exists for, and keeping only the last hid it.
        found.update(re.findall(rf"(?i){name}[-_]v?(\d+(?:\.\d+)+)", s))
        for word in s.split():
            m = re.search(r":v?(\d+(?:\.\d+)+)$", word)
            if m:
                found.add(m.group(1))
    return sorted(found)


def _reacquire_hint(engine, pin, tools):
    """How to get the manifest's pin into THIS tools.json: a lead-in ending in ':' and then
    each command on a line of its own, indented, so the line a user or agent copies is
    exactly a command. (The first version put the explanation after the command on the same
    line, and pasting it failed with `syntax error near unexpected token '('`.)

    The root matters: the pilot ran --tools .../pilot_tools/tools.json, and acquire_tools.sh
    without a root rewrites ~/.proteomics-pipeline/tools/tools.json instead, leaving the
    mismatch in place."""
    t = tools or {}
    cls = t.get("platform_class") or "<platform_class>"
    root = f" {shlex.quote(t['tools_root'])}" if t.get("tools_root") else ""
    # acquire_tools.sh re-resolves EVERY engine into that tools.json, and an engine it is not
    # pinning goes to "latest": on HIVE the pilot's pinned 2.7.0 DIA-NN entry would become
    # the Core's 2.6.1 after a Sage-only re-acquire. Say so rather than let the fix move it.
    caveat = (f"this rewrites every engine's entry in that tools.json; engines other than "
              f"{engine} are re-resolved unpinned")
    if engine == "diann" and cls == "mac":
        # acquire_tools.sh only wraps $DIANN_DOCKER_IMAGE on macOS; re-running it cannot
        # change the version inside the image. The image has to be built for the pin, and
        # build_diann_docker.sh tags it proteomics-pipeline/diann:<pin>. Passing that tag on
        # the acquire line works without re-sourcing activate.sh first.
        image = f"proteomics-pipeline/diann:{pin}"
        cmds = [f"bash scripts/build_diann_docker.sh {pin}",
                f"DIANN_DOCKER_IMAGE={image} bash scripts/acquire_tools.sh mac{root}"]
        caveat += f"; afterwards re-source activate.sh so new shells export {image}"
    else:
        cmds = [f"PIN_ENGINE={engine} PIN_VERSION={pin} bash scripts/acquire_tools.sh {cls}{root}"]
    return f"({caveat}):" + "".join(f"\n    {c}" for c in cmds)


def engine_version_record(engine, tools, bundle):
    """Which engine build this search ran, and whether it is the one the manifest asked for --
    the `engine_version` record of search_provenance.json.

    Three things name a version, and until this nothing compared them:
      * workflow.manifest.json `engine.version` -- resolve_defaults.py's PIN. A request.
      * tools.json `versions.<engine>`          -- what acquire_tools.sh resolved, beside the
                                                   command run_search.py executes VERBATIM. It
                                                   is read only if it is SHAPED like a version
                                                   (_concrete_version): "latest", "env",
                                                   "nightly" and the rest name no build.
      * that command itself                     -- a build folder, .sif name or image tag.
    The FRAN pilot (2026-09-16) is why they must be compared: manifest 2.6.1, tools.json 2.7.0,
    and the compute-node banner said "DIA-NN 2.7.0 Academia".

    Shaped like `scan_window`, its neighbour in that file: one object whose `value` is the
    answer and whose `source` says in words where it came from, with every input kept as
    written beside it:

      value             the build that runs: tools.json, else the ONE version the command
                        names, else None
      source            "tools.json versions.<engine> ..." | "named by the command ..." |
                        "unknown -- <why>"
      tools_json        tools.json `versions.<engine>` as written, even "latest" / "env"
      named_by_command  every version the command names
      manifest_pin      the manifest's pin as written; None unless `engine.name` IS this engine
      mismatch          True/False only when `value` and the pin are both known, else None

    `value` is only ever a build something CONFIRMS runs. The manifest's pin is NEVER promoted
    to it, because main() also writes it as top-level `version`, which
    fran_deposit.detect_engine() forwards -- and nothing else from this record -- to FRAN as
    the engine version. Promoting the pin whenever tools.json said "latest" would stamp an old
    tools.json pointing at build_260/diann-2.6.0 as 2.6.1, and every conda-env Sage as whatever
    the manifest pins. An unknown version is null; a wrong one is a false claim about published
    results. `mismatch` is null when the two cannot be compared -- "false" would claim they were
    checked and agree. A tools.json that contradicts its own command records `value` null with
    a WARNING: one of them is wrong and nothing here can tell which.
    """
    t = tools or {}
    tools_raw = (t.get("versions") or {}).get(engine)
    beng = (bundle or {}).get("engine") or {}
    # The manifest pins ONE engine. acquire_tools.sh writes `versions` keys only for diann,
    # sage and radiant, so `tools.versions or manifest.version` stamped `--engine fragpipe` or
    # `--engine alphadia` under a DIA-NN manifest with DIA-NN's version.
    # `== engine`, not `in (None, engine)`: resolve_defaults.py always writes `name`, and a
    # block that does not name one pins NOTHING. A nameless {"engine": {"version": "2.6.1"}}
    # used to be read as a pin for every engine at once, which stamped a false
    # `mismatch: true` and a "the manifest pins sage 2.6.1" WARNING onto a Sage search.
    manifest_raw = beng.get("version") if beng.get("name") == engine else None
    tv, mv = _concrete_version(tools_raw), _concrete_version(manifest_raw)
    cvs = command_engine_versions(engine, t)
    # What runs, as printed: Radiant's runtime prefix alone ("apptainer exec") says nothing.
    cmd = " ".join(str(x) for x in (t.get(engine), t.get(f"{engine}_image")) if x)
    out = sys.stderr.write

    if tv and cvs and tv not in cvs:
        version = None
        source = (f"unknown -- tools.json says {tv}, but the command it runs names "
                  f"{', '.join(cvs)} ({cmd})")
        out(f"[run_search] WARNING: tools.json says {engine} {tv}, but the command it runs names "
            f"{', '.join(cvs)} ({cmd}). One of them is wrong, so search_provenance.json records "
            f"`version: null` (both values are kept in `engine_version`). Re-run "
            f"acquire_tools.sh to rewrite tools.json from what it finds.\n")
    elif tv:
        version = tv
        source = (f"tools.json versions.{engine}, as acquire_tools.sh recorded it"
                  + ("; the command names it too" if cvs else
                     "; the command names no version to check it against"))
    elif len(cvs) == 1:
        version = cvs[0]
        source = f"named by the command tools.json runs ({cmd})"
        out(f"[run_search] NOTE: tools.json names no {engine} build ({tools_raw!r}); the "
            f"command it runs names {version} ({cmd}), recorded as `version`.\n")
    else:
        version = None
        # Two different failures, said differently: tools.json named nothing at all, or it
        # named something that is not a version ("latest", "env", and whatever a later branch
        # of acquire_tools.sh writes next). The second is the one a reader has to be able to
        # see, because that string is what would otherwise have been deposited as the version.
        said = (f"tools.json records {engine} {tools_raw!r}, which is not a version"
                if isinstance(tools_raw, str) and tools_raw.strip()
                else f"tools.json names no {engine} build ({tools_raw!r})")
        why = (f"the command it runs names several ({', '.join(cvs)})" if cvs
               else "the command it runs names none either")
        source = f"unknown -- {said}; {why}"
        pin = (f" The manifest's pin {mv} is kept as `engine_version.manifest_pin` but is not "
               f"recorded as `version`: it is what was asked for, not evidence of what ran."
               if mv else "")
        if tools_raw == "env":
            fix = (" (A sage on PATH has no release record, and `sage --version` cannot supply "
                   "one: the v0.14.7 release binary prints 0.14.6.)")
        elif mv:
            fix = f" To record the build, re-acquire it {_reacquire_hint(engine, mv, t)}"
        else:
            fix = ""
        out(f"[run_search] NOTE: {said}; {why}; "
            f"search_provenance.json records `version: null`.{pin}{fix}\n")

    mismatch = (version != mv) if (version and mv) else None
    if mismatch:
        out(f"[run_search] WARNING: engine version mismatch -- workflow.manifest.json pins "
            f"{engine} {mv}, but {version} is what runs ({cmd}). {version} is recorded as "
            f"`version`; the pin is kept beside it in search_provenance.json "
            f"(`engine_version.mismatch: true`). If {mv} was intended, re-acquire it "
            f"{_reacquire_hint(engine, mv, t)}\n")
    return {"value": version, "source": source,
            "tools_json": tools_raw or None,               # as written, even "latest"
            "named_by_command": cvs,
            "manifest_pin": manifest_raw,
            "mismatch": mismatch}


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
    # THE routing rule lives in diann_parallel.parallel_safe(), which the generator consumes
    # too, so the router can no longer decline a cfg the chain would happily run. That drift
    # is the whole bug: an unpinned --window -- which estimate_params.py omits BY DESIGN,
    # because the scan-window radius has to be measured against a real file, not guessed --
    # read as "not parallel-safe" here while diann_parallel was ready to measure it in step
    # 1b. The router won, and the fallback is not a slower path but a different order of
    # magnitude: --threads parallelises WITHIN a run, not across runs, so at ~30 min/file a
    # 310-file cohort is ~155 h sequential against a few hours for the chain. Nothing errors;
    # SLURM reports success and the user waits a week.
    #
    # The defaults below match what run_diann_parallel() actually passes today (probing on,
    # no seed library). If this ever grows --seed-lib or --no-probe-window passthroughs, they
    # must be forwarded here as well, or the two answers diverge again.
    #
    # The fix comes from the same verdict (`remedy`, keyed on WHY it declined). This used to
    # say "re-run estimate_params.py with the instrument table" for every decline -- including
    # a typo'd --window, which estimate_params.py cannot fix because it never writes one.
    safe = _diann_parallel_mod().parallel_safe(params)
    if safe["code"] in ("cfg_missing", "cfg_unparseable"):
        return False, f"{n} files, but {safe['reason']}. To fix: {safe['remedy']}"
    if not safe["ok"]:
        return False, (f"{n} files, but {params} is not parallel-safe for the 5-step chain "
                       f"(steps 3/5 reuse .quant files): {safe['reason']}. To enable it: "
                       f"{safe['remedy']}")
    return True, (f"{n} files > {a.parallel_threshold}, SLURM present, {safe['reason']}")


# Exit status when the search routed to the 5-step chain but --sbatch asked for ONE job
# script. The chain is generated; the file the caller is about to `sbatch` is not.
SBATCH_NOT_WRITTEN = 3


def _not_a_regular_file(path):
    """None if `path` is absent or a regular file; otherwise what it is (for the refusal)."""
    try:
        mode = os.lstat(path).st_mode
    except FileNotFoundError:
        return None
    if stat.S_ISREG(mode):
        return None
    return ("a directory" if stat.S_ISDIR(mode) else "a symlink" if stat.S_ISLNK(mode)
            else "not a regular file")


def set_aside(path):
    """Rename an existing REGULAR FILE out of the way -- never delete it. Returns the new path,
    or None when there was nothing there. The name records that it is stale and when it was
    moved. Anything else (a directory, a symlink, a device) raises ValueError: `--sbatch proj`
    once renamed a whole project folder, including the cfg the search was about to read."""
    kind = _not_a_regular_file(path)
    if kind:
        raise ValueError(f"{path} is {kind}, not a job script -- refusing to move it")
    if not os.path.lexists(path):
        return None
    base = f"{path}.stale-{time.strftime('%Y%m%dT%H%M%S')}"
    new, k = base, 1
    while os.path.lexists(new):
        new, k = f"{base}.{k}", k + 1
    os.rename(path, new)
    return new


def scan_window_record(engine, params, res):
    """What set the DIA-NN scan window, for search_provenance.json.

    The chain describes itself (res["scan_window"]: measured at run time by step 1b, or pinned
    in the cfg). The single-shot path is described from the cfg it ran. estimate_params.py
    cannot say which happened -- it runs before routing -- so its rationale points here."""
    if isinstance(res, dict) and res.get("scan_window"):
        return res["scan_window"]
    if engine != "diann":
        return None
    dp = _diann_parallel_mod()
    try:
        # the same description the chain uses when it does not probe -- exactly what was
        # passed, "unverified" where DIA-NN's handling has not been measured
        return dp.window_record(dp.mass_acc_status(params))
    except dp.CfgError as e:
        return {"source": f"unknown -- {e}", "value": None, "passed": None}


def ensure_temp_dirs(params, out):
    """Create any `--temp` directory the DIA-NN cfg names. Returns the list created.

    DIA-NN aborts with "cannot find the temp folder" when `--temp` points at a directory that
    does not exist -- it will NOT create it -- and it aborts BEFORE doing any work, so on a
    cluster the whole submission cycle is lost to a missing directory: queue, start, die, read
    the log, resubmit. This cost a real user a cycle. `mkdir -p` is idempotent and free.

    The 5-step parallel chain makes its own (each step, plus submit.sh). This covers the
    single-shot path, where `--temp` can only arrive from the cfg -- a workflow bundle's or the
    user's own. Relative paths resolve against the output directory, which is where DIA-NN runs.
    """
    made = []
    dp = _diann_parallel_mod()
    try:
        groups = dp.cfg_groups(dp.cfg_tokens(params))    # the one cfg reader
    except dp.CfgError:
        return made                         # ensure_xic, which runs first, reports it
    for flag, vals in groups:
        if flag == "--temp" and vals and not vals[0].startswith("-"):
            d = vals[0]
            d = d if os.path.isabs(d) else os.path.join(out, d)
            try:
                os.makedirs(d, exist_ok=True)
                made.append(d)
            except OSError as e:
                # Report it; do not raise. DIA-NN's own error names the folder, and failing here
                # would hide a permissions problem behind an unrelated traceback.
                sys.stderr.write(f"[run_search] could not create --temp {d}: {e}\n")
    if made:
        print(f"[run_search] --temp ready: {', '.join(made)} "
              f"(DIA-NN aborts rather than creating it)")
    return made


XIC_WINDOW_DEFAULT = 10          # seconds; DIA-NN's own default window


def ensure_xic(params, out):
    """Guarantee the DIA-NN cfg extracts XICs. Returns the cfg path to actually use.

    Chromatograms are what let a person *look* at an identification instead of trusting a
    q-value, and they are the input to FRAN's XIC lane — so every DIA-NN search this skill runs
    produces them. `estimate_params.py` already writes `--xic 10 --mobilograms`, but a cfg can
    also arrive from a workflow bundle or straight from the user, and those have no reason to
    carry it. Without this, whether a search has chromatograms depends on where its cfg came
    from — which is invisible until someone goes looking for a trace that was never extracted.

    The user's file is never edited: an augmented copy is written into the output directory, so
    the cfg that ran is recorded next to the results and the original stays exactly as given.

    `--mobilograms` rides along because `--xic` alone allocates the mobilogram parquets and
    leaves them full of zeros — silently, at plausible file size. DIA-NN ignores it on
    instruments without ion mobility.
    """
    # Read through the one cfg reader. `txt.split()` counted a commented-out `# --xic 10` as
    # present and added nothing, while the chain's reader saw no --xic -- so step 4 extracted
    # no XICs and nothing said so.
    dp = _diann_parallel_mod()
    try:
        flags = {f for f, _ in dp.cfg_groups(dp.cfg_tokens(params))}
    except dp.CfgError as e:
        if e.code == "cfg_missing":
            return params                   # nothing to augment; the caller will fail louder
        sys.exit(f"[run_search] {e}. Every route splices this cfg's flags into a command line "
                 f"or reads them to decide routing, so it must parse -- close the quote.")
    if "--xic" in flags:
        return params
    txt = open(params).read()
    os.makedirs(out, exist_ok=True)
    aug = os.path.join(out, "params_with_xic.cfg")
    extra = f"--xic {XIC_WINDOW_DEFAULT}" + ("" if "--mobilograms" in flags else " --mobilograms")
    with open(aug, "w") as fh:
        fh.write(txt.rstrip("\n") + f"\n{extra}\n")
    # The copy is the same parameters plus XICs, so it keeps the cfg's estimate_params.py
    # rationale: without it the mass-accuracy plan (measure_with_diann) would be lost here and the
    # search would fall back to DIA-NN's first-run auto mode (diann_parallel.mass_acc_measure_plan).
    # A cfg with NO sidecar must not inherit one an earlier search left in the same --out: that
    # stale plan would have a hand-written cfg's missing mass accuracy measured instead of reported.
    if os.path.exists(params + ".rationale.json"):
        shutil.copyfile(params + ".rationale.json", aug + ".rationale.json")
    elif os.path.lexists(aug + ".rationale.json"):
        os.remove(aug + ".rationale.json")
    print(f"[run_search] cfg had no --xic; using {aug} (added: {extra}). "
          f"Every DIA-NN search extracts chromatograms.")
    return aug


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
    # Chain resources are forwarded only when given, so diann_parallel.py's own defaults stay
    # the single definition of them. They are exposed because steps 3 and 5 ask for 64 CPUs
    # by default (--assembly-cpus), a whole node's worth that can sit pending for a long time
    # on a busy preemptible queue, and until now the only way to ask for less was to bypass
    # run_search.py and hand-run diann_parallel.py (FRAN pilot, publicgrp/low, 2026-09-16).
    for flag, val in (("--partition", a.partition), ("--account", a.account),
                      ("--qos", a.qos), ("--max-simultaneous", a.max_simultaneous),
                      ("--libpred-cpus", getattr(a, "libpred_cpus", None)),
                      ("--assembly-cpus", getattr(a, "assembly_cpus", None)),
                      ("--assembly-mem", getattr(a, "assembly_mem", None)),
                      ("--time-per-file", getattr(a, "time_per_file", None))):
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

# Dropped from the cfg for the single-shot SEARCH job, because that job supplies its own:
# --fasta-search/--predictor/--gen-spec-lib belong to the library job (in the search they would
# re-digest the FASTA instead of using the predicted library), and --reanalyse/--matrices are
# re-added explicitly. --rt-profiling is deliberately NOT here. It sets "the empirical library
# generation mode to IDs, RT and IM profiling" (DIA-NN README, command-line reference), which the
# README calls "strongly recommended for almost all workflows", and with --reanalyse this job is
# where the empirical library is generated -- MBR's first pass "creates an empirical spectral
# library from the data". Stripping it is not a no-op: measured on DIA-NN 2.7.0 (HIVE,
# 2026-09-16), the search log prints "The spectral library (if generated) will retain the
# original spectra but will include empirically-aligned RTs" with the flag and nothing without
# it, so the library mode silently changed. The 5-step chain passes it to exactly its
# library-building steps (2 and 3), and DE-LIMP's own single search keeps it. It was stripped
# here since e20ae63 with no stated reason. The rest of the cfg goes through
# diann_parallel.bash_flags(), the one emitter the chain uses too.
#
# DERIVED from the chain's STRIP rather than listed again. A hand-written list had drifted:
# --out, --f, --fasta, --threads, --lib and --out-lib were missing from it, so a cfg carrying
# any of them emitted the flag TWICE on the search command line -- and a cfg --out pointed
# DIA-NN's report somewhere the job neither clears beforehand nor checks afterwards, while
# clear_stale() deleted the report at the path the job believes in. The chain strips exactly
# these because the STEP supplies them, and so does this job.
#
# KEEP is the difference, and every entry is a flag this job does NOT supply itself, so the
# cfg is its only source:
#   --rt-profiling        the paragraph above
#   --temp                ensure_temp_dirs() creates the cfg's temp directory for this very
#                         command; the chain strips it only because each step passes its own
#   --xic / --mobilograms ensure_xic() puts them IN the cfg so that every search this skill
#                         runs extracts chromatograms; the chain strips --xic only to re-add
#                         it to step 4 alone (xic_flag())
#   --no-norm             step 5 re-adds it from diann_parallel's own --no-norm flag; the
#                         single-shot path has no such flag, so stripping it would silently
#                         discard what the user's cfg asked for
SINGLE_SHOT_SEARCH_KEEP = ("--rt-profiling", "--temp", "--xic", "--mobilograms", "--no-norm")
SINGLE_SHOT_SEARCH_STRIP = tuple(f for f in _diann_parallel_mod().STRIP
                                 if f not in SINGLE_SHOT_SEARCH_KEEP)


def report_guard(report, listing):
    """Bash that fails the job unless DIA-NN wrote the report AND it holds every input run.

    DIA-NN exits 0 on fatal errors, so the job state alone says nothing. must_exist() is the
    chain's own assertion; check_report_runs.py is the single-shot stand-in for the chain's
    step-5 .quant count, which has no per-file directory to count here. It runs under the
    interpreter that generated the job (on a cluster that path is on the shared filesystem,
    exactly like the scripts/ directory this references), falling back to python3."""
    checker = os.path.join(os.path.dirname(os.path.abspath(__file__)), "check_report_runs.py")
    return "\n".join([
        _diann_parallel_mod().must_exist(report, "the report"),
        f'PY={shlex.quote(sys.executable)}; [ -x "$PY" ] || PY=python3',
        f'"$PY" {shlex.quote(checker)} --report {shlex.quote(report)} '
        f'--files-list {shlex.quote(listing)}',
    ])


def single_shot_mass_acc(params, listing, fasta, predicted, out, threads, cmd, libfree):
    """Measure an undocumented Orbitrap mass accuracy before a single-shot search.

    For a cfg estimate_params.py planned `measure_with_diann` (diann_parallel.
    mass_acc_measure_plan: an Orbitrap level outside DIA-NN's table), the 5-step chain measures
    it in step 1b. A single-shot search used to get DIA-NN's first-run auto mode instead -- and
    every machine WITHOUT SLURM searches single-shot at any cohort size (parallel_decision
    returns "no SLURM here" first). DIA-NN 2.7.0 labels that mode "use this mode for preliminary
    analyses only", the result depends on which file sorts first, and nothing records the value
    it chose. SKILL.md golden rule 7: HIVE is a fast path, never a requirement. probe_window.py
    needs no SLURM, so the same probe runs here, between the library and the search: the same
    representative runs (and replacements), the same flags step 1b gives it -- as bash words
    after `--` -- and the same budget; the median of each measured level and the documented
    value of the other, pinned for the search via massacc.txt. params.resolved.cfg is built in a
    .tmp and moved into place only once the value is measured, as in step 1b.

    `listing` is the caller's per-job file list -- the probe reads it at RUN time, so it must
    be the job's own, never a shared name (see the comment on it below).

    Returns (bash lines to run before the search, the search's extra flags, mass-accuracy
    record, resolved-params record). All empty/None when nothing is planned. When this search
    has no library of its own to measure against -- --one-step, an external --lib, or a cfg that
    neither predicts nor supplies one -- nothing is measured, DIA-NN optimises on the first run,
    and it is SAID, naming which of the three it is, in the job output and in the record."""
    dp = _diann_parallel_mod()
    plan = dp.mass_acc_measure_plan(params)
    if not plan:
        return [], "", None, None
    if not libfree:
        # Three different reasons land here and they have three different fixes. Saying "no
        # library to measure against" and "drop --one-step" for all of them is wrong for the
        # commonest of them: a cfg with an external --lib HAS a library, and there is no
        # --one-step to drop.
        try:
            present = {f for f, _ in dp.cfg_groups(dp.cfg_tokens(params))}
        except dp.CfgError:
            present = set()
        if getattr(a_globals, "one_step", False):
            cause = ("--one-step was given, so the library is predicted inside the search itself "
                     "and does not exist before it")
            fix = "Drop --one-step to measure it."
        elif "--lib" in present:
            cause = (f"{params} searches against a library it did not build (--lib), and this "
                     "route measures only against the library the search predicts, so the probe "
                     "has nothing wired up to run against")
            fix = ("Pin --mass-acc and --mass-acc-ms1 in the cfg, or pass "
                   "--ms1-resolution/--ms2-resolution to estimate_params.py so DIA-NN's own "
                   "Orbitrap table can supply them.")
        else:
            cause = (f"{params} neither predicts a library (--fasta-search --gen-spec-lib) nor "
                     "supplies one (--lib), so there is nothing to measure against")
            fix = ("Add --fasta-search --gen-spec-lib to the cfg, or pin --mass-acc and "
                   "--mass-acc-ms1 in it.")
        why = (f"{params} plans to measure mass accuracy with DIA-NN, but {cause}. DIA-NN will "
               "optimise it on the first run of the search instead -- 'use this mode for "
               f"preliminary analyses only'. {fix}")
        sys.stderr.write(f"[run_diann] NOTE: {why}\n")
        return [], "", {"fixed": False, "measured": False, "reason": why}, None
    q = shlex.quote
    # The CALLER's list file, named after the job (run_diann: `{stem}_input_files.txt`), not a
    # fixed search_input_files.txt of our own. This probe is generated now and read by the job
    # later, so it has exactly the hazard report_guard()'s listing has: generating a second
    # search into the same --out would rewrite the list the FIRST, still unsubmitted, job
    # probes from, and that job would measure mass accuracy on the other cohort and pin it for
    # its own runs. A job, its report guard and its probe now rise and fall together.
    probe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "probe_window.py")
    massacc = os.path.join(out, "massacc.txt")
    evidence = os.path.join(out, "mass_acc.json")
    resolved = os.path.join(out, "params.resolved.cfg")
    tmp = resolved + ".tmp"
    workdir = os.path.join(out, "mass_acc_probe")
    documented = plan["documented"]
    doc_args = (" " + dp.probe_mass_acc_args(documented)) if documented else ""
    # the same flags the chain's step 1b probes with: the cfg minus the step-specific ones
    flags = dp.read_cfg_flags(params, drop=dp.MASS_ACC_FLAGS)
    pattern = dp.MEASURED_FILE_RE["massacc.txt"]
    failed = (f'{{ echo "FAILED: no mass accuracy measured -- see the messages above; per-run '
              f'evidence in {evidence}. Nothing was searched." >&2; rm -f {q(tmp)}; exit 1; }}')
    lines = [
        'echo "measuring mass accuracy with DIA-NN on representative runs before the search"',
        # an earlier run's value, evidence or resolved cfg must not survive
        f"rm -f {q(evidence)} {q(massacc)} {q(resolved)} {q(tmp)}",
        f"rm -rf {q(workdir)}",
        f"cp {q(params)} {q(tmp)}",
        f"python3 {q(probe)} --diann {q(cmd)} --raw-list {q(listing)} --fasta {q(fasta)} "
        f"--lib {q(predicted)} --threads {threads} --measure mass-acc{doc_args} "
        f"--max-probes {dp.PROBE_CANDIDATES} --max-failures {dp.PROBE_MAX_FAILURES} "
        f"--timeout {dp.PROBE_TIMEOUT_S} --budget {dp.PROBE_BUDGET_S} "
        f"--workdir {q(workdir)} --write-cfg {q(tmp)} -- {flags} > {q(evidence)} || {failed}",
        # the two flags, in exactly the shape the guard below accepts
        "python3 -c \"import json,re,sys; m=json.load(open(sys.argv[1]))['mass_acc']['pin_as']; "
        f"assert re.fullmatch(sys.argv[2], m); print(m)\" {q(evidence)} {q(pattern)} "
        f"> {q(massacc)} || {failed}",
        dp.needs_measured(massacc, "the two mass-accuracy flags",
                          producer="the probe before this search"),
        f"mv -f {q(tmp)} {q(resolved)}",
        f'echo "mass accuracy = $(cat {q(massacc)}) (pinned for the search; evidence {evidence})"',
    ]
    record = {"fixed": True, "measured": True, "documented": documented,
              "ms1": documented.get("--mass-acc-ms1"), "ms2": documented.get("--mass-acc"),
              "source": "measured with DIA-NN on representative runs by the single-shot search "
                        "job's probe before the search (probe_window.py --measure mass-acc); "
                        "planned by estimate_params.py (measure_with_diann) for an Orbitrap level "
                        "with no documented DIA-NN value",
              "value_file": massacc, "evidence_file": evidence,
              "sop_floor": dict(dp.SOP_MASS_ACC_FLAGS),
              "floor_note": dp.floor_note(evidence, massacc),
              "reason": "not in the cfg; measured before the search and passed to it"}
    resolved_params = {"file": resolved, "produced": "runtime",
                       "by": "the single-shot search's pre-search probe",
                       "note": "written only after mass accuracy is measured, just before the "
                               "DIA-NN search; absent until then, so a missing file after the "
                               "search job means the measurement failed"}
    return lines, f" $(cat {q(massacc)})", record, resolved_params


def run_diann(cmd, params, files, fasta, out, threads, sbatch, acquisition="", queue=None):
    """Single-shot DIA-NN. `queue` is {partition, account, qos} from the command line, handed
    to every emit_sbatch() call -- without it slurm_queue() re-detects a queue and overrides
    the one the user chose."""
    queue = queue or {}
    os.makedirs(out, exist_ok=True)
    report = os.path.join(out, "report.parquet")
    f_args = " ".join(f"--f {shlex.quote(f)}" for f in files)
    # Named after the JOB, not fixed, because report_guard() bakes this path into the script
    # and the script is what decides whether the search was complete. One shared
    # search_input_files.txt meant generating a second search into the same --out rewrote the
    # list the FIRST, still unsubmitted, job would read -- so that job would check its report
    # against a different cohort than the one it searches. A job and its list now rise and
    # fall together: regenerating the same --sbatch rewrites both, and a different --sbatch
    # gets its own. The inline path keeps the old name; nothing outlives the process there.
    stem = os.path.splitext(os.path.basename(sbatch))[0] if sbatch else "search"
    listing = os.path.join(out, f"{stem}_input_files.txt")
    with open(listing, "w") as fh:                    # a list file survives spaces in paths
        fh.write("\n".join(files) + "\n")
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
    # The cfg is read by the same reader, and its flags emitted by the same quoting, as the
    # parallel chain's. `cfg_txt.split()` spliced raw words into bash: a `# comment` commented
    # out --f/--fasta/--out, and `--cut K*,R*` went in bare, where one matching file in the
    # output directory rewrites the digest rule. Comment-free cfgs emit the same command as
    # before apart from that quoting and the --rt-profiling this job now keeps
    # (SINGLE_SHOT_SEARCH_STRIP); tests/test_cfg_reader_quoting.py pins both.
    dp = _diann_parallel_mod()
    try:
        groups = dp.cfg_groups(dp.cfg_tokens(params))
    except dp.CfgError as e:
        if e.code != "cfg_missing":
            sys.exit(f"[run_diann] {e}")
        groups = []                         # onecmd below hands DIA-NN the path; it reports it
    present = {f for f, _ in groups}
    libfree = {"--fasta-search", "--gen-spec-lib"} <= present \
        and not getattr(a_globals, "one_step", False)
    dnet = dotnet_env_for(files)          # .NET 8 for reading Thermo .raw, if needed

    lib = os.path.join(out, "diann_lib")
    # DIA-NN ignores the extension asked of --out-lib for a PREDICTED library and always writes
    # <name>.predicted.speclib (README, "Output library") -- the file the search job reads.
    predicted = lib + ".predicted.speclib"
    search_cfg = dp.bash_flags(groups, drop=SINGLE_SHOT_SEARCH_STRIP)
    lib_cmd = (f"{cmd} --cfg {shlex.quote(params)} --fasta {shlex.quote(fasta)} "
               f"--out-lib {shlex.quote(lib)} --threads {threads}")
    # An undocumented Orbitrap mass accuracy is measured between the library and the search
    # (single_shot_mass_acc); `mflag` then carries the pinned flags into the search.
    measure_lines, mflag, mass_acc, resolved_params = single_shot_mass_acc(
        params, listing, fasta, lib + ".predicted.speclib", out, threads, cmd, libfree)
    # recorded only when there is something to say, so a pinned cfg's result is unchanged
    ma_rec = {} if mass_acc is None else {"mass_acc": mass_acc}
    if resolved_params:
        ma_rec["resolved_params"] = resolved_params
    # --temp, always. Without it DIA-NN writes every run's .quant NEXT TO THE RAW FILE -- on
    # HIVE that is the instrument archive (/nfs/lssc0/flinders/.../raw_data), shared by every
    # user who searches those runs, so two searches race on the same .quant files and the
    # archive fills with search byproducts (HIVE e2e test 2026-09-23: DIA-NN 2.7.0 job held
    # before it ran; the FL*.raw.quant beside older HeLa raws show it had happened). The 5-step
    # chain passes its own per step; a --temp in the cfg still wins (ensure_temp_dirs).
    tmp_arg = ""
    if "--temp" not in present:
        quant_dir = os.path.join(os.path.abspath(out), "quant")
        os.makedirs(quant_dir, exist_ok=True)           # DIA-NN aborts rather than create it
        tmp_arg = f" --temp {shlex.quote(quant_dir)}"
    search_cmd = (f"{cmd} {search_cfg}{mflag} {f_args} --fasta {shlex.quote(fasta)} "
                  f"--lib {shlex.quote(lib)}.predicted.speclib --reanalyse --matrices "
                  f"--out {shlex.quote(report)} --threads {threads}{dda}{tmp_arg}")
    onecmd = (f"{cmd} --cfg {shlex.quote(params)} {f_args} "
              f"--fasta {shlex.quote(fasta)} --out {shlex.quote(report)} "
              f"--threads {threads}{dda}{tmp_arg}")
    # DIA-NN exits 0 on fatal errors, so each job asserts the artefact it exists to make --
    # the same contract every step of the 5-step chain has (references/diann_parallel.md).
    # And each job first DELETES that artefact, because an existence check cannot tell this
    # run's file from the previous search's in the same --out. Reproduced on HIVE, DIA-NN
    # 2.7.0 (review srun 23512013): a re-run with no DOTNET_ROOT logged "ERROR: cannot read
    # .raw files", exited 0 and wrote nothing; the old report.parquet and report.stats.tsv
    # were byte-identical afterwards, and both guards passed on them. The stats file goes too:
    # it is check_report_runs.py's fallback evidence when there is no parquet reader.
    import check_report_runs
    report_files = (report, check_report_runs.stats_path(report))
    lib_job = "\n".join([dp.clear_stale(predicted), lib_cmd,
                         dp.must_exist(predicted, "the predicted spectral library")])
    guard = report_guard(report, listing)
    search_job = "\n".join([dp.clear_stale(*report_files), *measure_lines, search_cmd, guard])
    one_job = "\n".join([dp.clear_stale(*report_files), onecmd, guard])

    # TWO JOBS + dependency when emitting sbatch: the library is expensive and
    # reusable, so a failed search requeues against it instead of rebuilding.
    if libfree and sbatch:
        # ABSOLUTE paths: submit.sh lives in <out>, the job scripts beside --sbatch, and
        # run_search.py prints `bash <out>/submit.sh` -- run from anywhere but the --sbatch
        # folder, relative names made sbatch fail "Unable to open file job_1_lib.sh" (HIVE e2e
        # test 2026-09-23). splitext, not replace(".sh", ""): that also ate a ".sh" inside a
        # folder name.
        stem = os.path.splitext(os.path.abspath(sbatch))[0]
        lib_sh, srch_sh = stem + "_1_lib.sh", stem + "_2_search.sh"
        emit_sbatch(lib_sh, lib_job, out, threads, job="diann_libpred", preamble=dnet,
                    submit_hint=False, **queue)
        # The measurement belongs to the SEARCH job: a search requeued against the same library
        # measures again rather than trusting a massacc.txt from a run it cannot vouch for.
        emit_sbatch(srch_sh, search_job, out, threads, job="diann_search", preamble=dnet,
                    hours=search_job_hours(measure_lines), submit_hint=False, **queue)
        submit = os.path.join(os.path.abspath(out), "submit.sh")
        jobs_txt = os.path.join(os.path.abspath(out), "jobs.txt")
        # jobs.txt, as the 5-step chain's submit.sh writes it: `watch_run.sh --all <out>`
        # reads it, and without one it answered failed/no_jobs_file for a healthy running
        # search -- which step 7b says to resubmit (HIVE e2e test 2026-09-23).
        with open(submit, "w") as fh:
            fh.write("#!/bin/bash -l\nset -euo pipefail\n"
                     f"j1=$(sbatch --parsable {shlex.quote(lib_sh)})\n"
                     f"j2=$(sbatch --parsable --dependency=afterok:$j1 {shlex.quote(srch_sh)})\n"
                     f'printf "%s\\n" "$j1" "$j2" > {shlex.quote(jobs_txt)}\n'
                     'echo "submitted: libpred=$j1 search=$j2"\n'
                     f'echo "both job ids -> {jobs_txt}  (watch both: watch_run.sh --all '
                     f'{os.path.abspath(out)})"\n'
                     f'echo "report will be {report}"\n')
        os.chmod(submit, 0o755)
        print(f"  [sbatch] two-job chain: {lib_sh} -> {srch_sh}; submit with: bash {submit}")
        print(f"  [sbatch] --sbatch {sbatch} itself is NOT written for a library-free search: "
              f"submit.sh submits both jobs.")
        return {"engine": "diann", "report": report, "submitted": submit,
                "mode": "two_job_libfree", "library": predicted,
                "ran": False, "dda": bool(dda), "raw_dotnet": bool(dnet), **ma_rec}

    if sbatch:                          # libfree + sbatch returned above, so this is onecmd
        emit_sbatch(sbatch, one_job, out, threads, job="diann_search", preamble=dnet,
                    hours=search_job_hours(measure_lines), **queue)
        return {"engine": "diann", "report": report, "submitted": sbatch,
                "ran": False, "dda": bool(dda), "raw_dotnet": bool(dnet), **ma_rec}
    pre = (dnet + " ") if dnet else ""
    # The sbatch route CAN only delete (clear_stale runs on a compute node, long after this
    # process is gone). Here we are the ones re-running the search, so the previous results are
    # RENAMED out of the way instead -- same effect on the guards below, and a re-run whose
    # DIA-NN dies still leaves the user the report they had. set_aside() is the idiom --sbatch
    # already uses (`existing_file_moved_to`), including its refusal to move anything that is
    # not a regular file: `--out` pointing at a directory called report.parquet must not be
    # renamed away.
    moved = []
    for p in ((predicted,) if libfree else ()) + report_files:   # see clear_stale above
        try:
            old = set_aside(p)
        except ValueError as e:
            sys.exit(f"[run_diann] {e}. Move it yourself, or search into a different --out.")
        if old:
            moved.append(old)
    if moved:
        print(f"  [run_diann] previous artefacts set aside: {', '.join(moved)}")
    if libfree:
        sh(pre + lib_cmd)
        if not (os.path.exists(predicted) and os.path.getsize(predicted) > 0):
            sys.exit(f"DIA-NN exited 0 but did not write the predicted library {predicted} "
                     "-- check its output above for ERROR:")
        if measure_lines:
            # Between the library and the search, in ONE shell that stops at the first
            # failure, so nothing is searched with an unmeasured mass accuracy.
            try:
                sh((dnet + "\n" if dnet else "") + "\n".join(["set -e", *measure_lines]))
            except subprocess.CalledProcessError as e:
                sys.exit(f"single-shot DIA-NN search stopped (exit {e.returncode}) -- see "
                         "the messages above; nothing was searched with an unmeasured mass "
                         "accuracy")
        sh(pre + search_cmd)
    else:
        sh(pre + onecmd)
    # `-s`, not `-e`: the sbatch route's must_exist() tests `[ -s ]`, and a 0-byte report is
    # exactly what DIA-NN leaves when it creates the file and then dies -- which reached
    # verify() as an unhandled parquet-reader traceback instead of this sentence.
    if not (os.path.exists(report) and os.path.getsize(report) > 0):
        sys.exit(f"DIA-NN finished but {report} is missing or empty "
                 "-- check its output above for ERROR: (DIA-NN exits 0 on a fatal error).")
    ok, msg = check_report_runs.verify(report, files)
    if not ok:
        sys.exit(msg)
    print(f"  [run_diann] {msg}")
    return {"engine": "diann", "report": report, "ran": True, "dda": bool(dda),
            "previous_artefacts_moved_to": moved, **ma_rec}


# --------------------------------------------------------------- AlphaDIA -----
# Apache-2.0 (commercial use OK) — the open-source DIA alternative to DIA-NN,
# whose free "Academia" build is academic/non-profit only. Library-free:
#   alphadia -o <out> -f <raw> [-f ...] --fasta <fasta> [-c <config.yaml>]
def run_alphadia(cmd, config, files, fasta, out, threads, sbatch, queue=None):
    os.makedirs(out, exist_ok=True)
    f_args = " ".join(f"-f {shlex.quote(f)}" for f in files)
    cfg = f"-c {shlex.quote(config)} " if config and os.path.exists(config) else ""
    full = (f"{cmd} -o {shlex.quote(out)} {f_args} --fasta {shlex.quote(fasta)} {cfg}").strip()
    if sbatch:
        emit_sbatch(sbatch, full, out, threads, job="alphadia_search", **(queue or {}))
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
    # AlphaDIA has already applied its own FDR, so the contract's q-columns are
    # emitted as 0.0 placeholders -- the downstream filter is deliberately a
    # no-op here rather than a second, different FDR. Names come from the shared
    # definition so this adapter cannot fall behind the contract it satisfies.
    pq.write_table(pa.table({
        "Run": runs, "Protein.Group": prots, "PG.MaxLFQ": ints,
        **{c: [0.0] * n for c in FDR_REQUIRED},
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


def run_radiant(tools, params, files, fasta, out, threads, sbatch, library=None, mbr=True,
                queue=None):
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
        emit_sbatch(sbatch, full, out, threads, job="radiant_search", **(queue or {}))
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
    # Preference order, NOT a filter set: the first available column wins. Shared
    # with compare_searches.py so the two cannot disagree about which q-value a
    # given report is being judged on (SKILL_OPEN_DEFECTS #2). Widens the old
    # two-column chain -- Global.Q.Value / Lib.PG.Q.Value / PG.Q.Value are now
    # tried before falling all the way back to the run-level Q.Value.
    c_q = col(*PROTEIN_Q_PREFERENCE)

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
    # Fulcrum reports ONE q-value per protein x run, so it is broadcast to every
    # column of the contract rather than invented per column. Names derived from
    # the shared definition -- see the AlphaDIA adapter above.
    pq.write_table(pa.table({
        "Run": runs2, "Protein.Group": pgs2, "PG.MaxLFQ": ints,
        **{c: qv for c in FDR_REQUIRED},
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


def run_sage(cmd, params, files, fasta, out, threads, sbatch, queue=None):
    os.makedirs(out, exist_ok=True)
    mzml = ensure_mzml(files, out)
    files_args = " ".join(shlex.quote(m) for m in mzml)
    full = (f"{cmd} {shlex.quote(params)} -f {shlex.quote(fasta)} -o {shlex.quote(out)} "
            f"--parquet --disable-telemetry-i-dont-want-to-improve-sage {files_args}")
    if sbatch:
        emit_sbatch(sbatch, full, out, threads, job="sage_search", **(queue or {}))
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
        # Sage already FDR-filtered at write time, so these are 0.0 placeholders.
        **{c: [0.0] * n for c in FDR_REQUIRED},
    })
    report = os.path.join(out, "report.parquet")
    pq.write_table(out_tbl, report)
    print(f"  [adapt] Sage -> {report}  ({n} protein×run rows)")
    return report


# --------------------------------------------------------------- FragPipe -----
def run_fragpipe(cmd, bundle, params, files, fasta, out, threads, sbatch, queue=None):
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
        emit_sbatch(sbatch, full, out, threads, job="fragpipe_search", **(queue or {}))
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
    # Carry each contract q-column through when the source has it, else 0.0.
    # Driven off FDR_REQUIRED so adding a column to the contract cannot leave
    # this adapter silently emitting one fewer.
    src_q = {c: col(c) for c in FDR_REQUIRED}
    n = len(rows)
    tbl = pa.table({
        "Run": [r[c_run] for r in rows],
        "Protein.Group": [r[c_pg] for r in rows],
        "PG.MaxLFQ": [num(r[c_q]) for r in rows],
        **{c: [num(r[src_q[c]]) if src_q[c] else 0.0 for r in rows]
           for c in FDR_REQUIRED},
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
                    # already FDR-filtered upstream -> 0.0 placeholders
                    **{c: [0.0] * n for c in FDR_REQUIRED}})
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
    A PARTIAL override (say `--partition low` alone) is completed from an association that
    matches every field given, never from the preferred one: filling `--partition low` with
    genome-center-grp's account and QOS wrote a header HIVE rejects ("Invalid account or
    account/partition combination specified", srun --test-only, 2026-09-16), and `--qos
    publicgrp-low-qos` alone on `high` gave "Invalid qos specification". When no association
    matches, this EXITS naming the ones that exist -- there is nothing valid to fill in.
    A COMPLETE `--partition X --account Y` pair is still honoured whatever the associations
    say (nothing has to be invented, and a reservation may exist that sacctmgr does not show),
    but it is checked against them and a mismatch is said out loud -- it used to skip the check
    entirely, so the only overrides that got validated were the incomplete ones.
    Returns (partition, account, qos); any may be None, and a None is simply omitted
    from the script so SLURM applies its own default."""
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

    # An empty partition or QOS field in an association means "no restriction", which is a
    # fine DEFAULT but not a licence to confirm a value the user typed: a blank QOS field used
    # to match anything, so `--qos totally-bogus-qos` was approved here and written into the
    # header, where SLURM rejects it hours later ("Invalid qos specification"). A blank field
    # now matches only when nothing was given for it -- which is every detection call, so a
    # cluster laid out differently from HIVE is still not refused on a guess.
    def fits(a, p, q):
        return ((not account or a == account)
                and (not partition or p == partition)
                and (not qos or qos in (q or "").split(",")))

    given = " ".join(f"--{k} {v}" for k, v in (("partition", partition),
                     ("account", account), ("qos", qos)) if v)
    shown = ", ".join("|".join(x) for x in assoc)

    if partition and account:
        # Complete: nothing is guessed, so it is honoured either way. But it used to return
        # here BEFORE the associations were even read, which is why `--partition low --account
        # genome-center-grp` sailed through generation and died at submit time with "Invalid
        # account or account/partition combination specified". Checked, and said out loud.
        if assoc and not [x for x in assoc if fits(*x)]:
            sys.stderr.write(
                f"[slurm_queue] WARNING: no SLURM association of user "
                f"{os.environ.get('USER', '?')} has {given}; SLURM will most likely reject "
                f"this job. Associations (account|partition|qos): {shown}. Using it anyway "
                f"because you named a complete queue -- check the #SBATCH header.\n")
        return partition, account, _public_low_qos(partition, account, qos)

    def find(acct, part):
        for a, p, q in assoc:
            if a == acct and p == part:
                return a, p, (q or None)
        return None

    lab, pub = find("genome-center-grp", "high"), find("publicgrp", "low")

    if (partition or account or qos) and assoc:
        hits = [x for x in assoc if fits(*x)]
        if not hits:
            sys.exit(f"[slurm_queue] no SLURM association of user "
                     f"{os.environ.get('USER', '?')} has {given}, so SLURM would reject the "
                     f"job. Associations (account|partition|qos): "
                     f"{shown}. Pass --partition/--account/"
                     f"--qos from ONE of them, or none to have one chosen.")
        preferred = [h for key in (("genome-center-grp", "high"), ("publicgrp", "low"))
                     for h in hits if h[:2] == key]
        cands = preferred or hits
        # Several associations match what was given, and `[0]` picked whichever sacctmgr
        # happened to list first -- silently. Two lab associations on `low` resolved to labA
        # with no output at all, and that decides WHO IS BILLED for a multi-hour search. A
        # differing ACCOUNT therefore stops here; a differing partition/QOS under one account
        # is a scheduling choice, so it is announced and the preferred order stands.
        uniq = sorted(set(cands))
        if len({x[0] for x in uniq}) > 1:
            sys.exit(f"[slurm_queue] {given} matches more than one SLURM ACCOUNT "
                     f"({', '.join(sorted({x[0] for x in uniq}))}), and choosing one decides "
                     f"which is billed for this run. Matching associations: "
                     f"{', '.join('|'.join(x) for x in uniq)}. Add --account to say which.")
        if len(uniq) > 1:
            sys.stderr.write(f"[slurm_queue] WARNING: {given} matches {len(uniq)} associations "
                             f"({', '.join('|'.join(x) for x in uniq)}); using "
                             f"{'|'.join(cands[0])}. Pass --partition/--qos to choose.\n")
        a, p, q = cands[0]
        p = partition or p or None
        # a comma-separated QOS list names no single QOS to write; SLURM picks the default
        q = qos or (q if q and "," not in q else None)
        return p, a, _public_low_qos(p, a, q)

    # Port of DE-LIMP's select_best_partition() (R/helpers_search.R). Entitlement is
    # not the question -- UTILISATION is. The priority queue has a PER-USER CPU cap
    # (64 on HIVE), and once you are at it your own jobs queue behind each other:
    # an 18-task array on `high` starves everything else you submit (QOSGrpCpuLimit,
    # observed). publicgrp/low is preemptible but has thousands of idle CPUs, so for
    # work that is safe to preempt it starts sooner and finishes sooner.
    if lab and pub and not (partition or account or qos):
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
    if partition or account or qos:
        print("[slurm_queue] WARNING: cannot read SLURM associations here, so the partial "
              f"queue (partition={partition}, account={account}, qos={qos}) is completed "
              "with publicgrp/low unchecked; read the #SBATCH header before submitting",
              file=sys.stderr)
    p, a = partition or "low", account or "publicgrp"
    return p, a, _public_low_qos(p, a, qos)


def _public_low_qos(partition, account, qos):
    """The QOS for publicgrp on `low` when none was given.

    diann_parallel.py and diatracer_parallel.py already add `publicgrp-low-qos` there ("low
    DOES need its qos named", tests/test_hive_submission_guards.py); emit_sbatch() did not, so
    `--partition low --account publicgrp` wrote a --qos line in the 5-step chain and none in
    the single-shot jobs. Here, every caller gets it. `high` is left alone on purpose: a
    facility job with no --qos is accepted and SLURM assigns genome-center-grp-high-qos
    (measured 2026-08-25, test_high_needs_no_explicit_qos)."""
    if not qos and partition == "low" and account == "publicgrp":
        return "publicgrp-low-qos"
    return qos


def _positive_int(v):
    """argparse type: a CPU count or an hour count of 0 is a job SLURM rejects outright."""
    try:
        n = int(v)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a whole number, got {v!r}")
    if n < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1, got {n}")
    return n


# Wall clock for a search job with nothing in front of it. A job that ALSO measures mass accuracy
# gets search_job_hours() instead -- see there.
SEARCH_WALL_HOURS = 12


def search_job_hours(measure_lines):
    """Hours for a --sbatch search job, given the lines that run before the search.

    The pre-search mass-accuracy probe is handed diann_parallel.PROBE_BUDGET_S, which is sized
    against step 1b's OWN wall clock (PROBE_WALL_HOURS, 3 full probes plus an hour) -- step 1b is
    a job of its own and nothing else runs in it. Here the same probe runs INSIDE the search job,
    so at the default 12 h a probe that used its whole budget would have eaten 3h50 of the
    search's wall before DIA-NN started, and a long search would be cut by SLURM with the
    measurement done and nothing to show for it. The budget stays what step 1b proved; the job
    gets that much more wall, so the search still has its full 12 hours."""
    if not measure_lines:
        return SEARCH_WALL_HOURS
    return SEARCH_WALL_HOURS + _diann_parallel_mod().PROBE_WALL_HOURS


def emit_sbatch(path, command, out, threads, job, preamble="",
                partition=None, account=None, qos=None, mem="64G", hours=SEARCH_WALL_HOURS,
                submit_hint=True):
    """Emit a minimal SLURM script (login-node-safe). Orchestrator submits it.
    The queue is DETECTED from the submitting user's own SLURM associations — see
    slurm_queue() — unless the caller passes one, which then wins. Every caller must forward
    the command line's --partition/--account/--qos: on the 2026-09-16 FRAN pilot run_diann()
    did not, so `--partition low --account publicgrp --qos publicgrp-low-qos` came out as
    genome-center-grp/high in both job headers and had to be hand-edited before submission.
    `preamble` runs before the command (e.g. the DOTNET_ROOT exports that let DIA-NN 2.6
    read Thermo .raw)."""
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
    # Preemptible queue: requeue a preempted search rather than lose it (HIVE's JobRequeue=1
    # does this by default; another cluster's may not). The rule is
    # diann_parallel.needs_requeue(), shared with the chain's step headers.
    requeue = _diann_parallel_mod().needs_requeue(part, q)
    if requeue:
        lines.append("#SBATCH --requeue")
    lines += ["set -euo pipefail", f"cd {shlex.quote(os.path.abspath(out))}", f"{pre}{command}", ""]
    script = "\n".join(lines)
    with open(path, "w") as fh:
        fh.write(script)
    print(f"  [sbatch] wrote {path} (partition={part or 'default'}, "
          f"account={acct or 'default'}, qos={q or 'default'}"
          f"{', requeue' if requeue else ''})"
          # A script that is one link of a chain must not be advertised as submittable on its
          # own: `sbatch job_2_search.sh` skips the afterok on the library job.
          + (f" — submit with: sbatch {path}" if submit_hint else " — submitted by submit.sh"))


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
    ap.add_argument("--partition", help="SLURM partition for every job this writes (the "
                    "single-shot --sbatch job(s) and the parallel chain); detected if omitted")
    ap.add_argument("--account", help="SLURM account for every job this writes")
    ap.add_argument("--qos", help="SLURM QOS for every job this writes")
    ap.add_argument("--max-simultaneous", type=int,
                    help="cap concurrent array tasks in the parallel chain")
    ap.add_argument("--libpred-cpus", type=_positive_int,
                    help="parallel chain: CPUs for step 1, library prediction "
                         "(diann_parallel.py default: 16)")
    ap.add_argument("--assembly-cpus", type=_positive_int,
                    help="parallel chain: CPUs for steps 3 and 5, assembly and report "
                         "(diann_parallel.py default: 64, a whole node's worth that can "
                         "wait a long time on a busy preemptible queue)")
    ap.add_argument("--assembly-mem", type=_positive_int,
                    help="parallel chain: memory in GB for steps 3 and 5 (diann_parallel.py "
                         "default: 128); lower it with --assembly-cpus, or a smaller job "
                         "still waits for 128 GB")
    ap.add_argument("--time-per-file", type=_positive_int,
                    help="parallel chain: wall-clock hours per array task in steps 2 and 4 "
                         "(diann_parallel.py default: 2)")
    ap.add_argument("--adapt-only", action="store_true",
                    help="skip the search; just build report.parquet from an existing engine output dir")
    ap.add_argument(ALLOW_DAMAGED_TDF, action="store_true",
                    help="search a Bruker .d whose analysis.tdf is not `ok` (truncated, "
                         "stale -wal/-journal beside it, WAL-mode header, unreadable) "
                         "anyway; every one is printed first")
    a = ap.parse_args()

    global a_globals
    a_globals = a
    # The queue the user asked for, handed to every emit_sbatch() -- see emit_sbatch().
    queue = {"partition": a.partition, "account": a.account, "qos": a.qos}
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
    # Every DIA-NN job this writes puts its --out inside DOUBLE quotes (clear_stale(),
    # must_exist()), on purpose, so that an array task's $QUANT expands. A `$(...)` in --out
    # therefore runs when the JOB runs -- and it lands inside the `rm -f --` the search does
    # before DIA-NN starts, so it also chooses what gets deleted. One rule, in the generator
    # that emits those quotes.
    _diann_parallel_mod().refuse_unsafe_path(a.out, prog="run_search")
    for attr in ("params", "fasta", "library"):
        v = getattr(a, attr, None)
        if v:
            setattr(a, attr, os.path.abspath(v))
    files = [os.path.abspath(f) for f in files]

    # DIA-NN names a run by its file name without the folder, so /plate1/s1.raw and
    # /plate2/s1.raw become ONE Run: two samples merged in the report (and, in the chain, two
    # array tasks writing the same .quant). That is knowable now, from the input list, so it
    # stops here instead of failing check_report_runs.py after the search has run.
    if engine == "diann" and not a.adapt_only:
        import check_report_runs
        dupes = check_report_runs.duplicate_run_names(files)
        if dupes:
            sys.exit(f"[run_search] inputs share a run name: {', '.join(dupes)}. DIA-NN names a "
                     "run by its file name without the folder, so they would be merged into "
                     "one Run. Rename them, or search them separately.")

    if a.adapt_only:
        report = {"sage": adapt_sage, "fragpipe": adapt_fragpipe,
                  "alphadia": adapt_alphadia,
                  "radiant": adapt_radiant}.get(engine, lambda o: None)(a.out)
        print(json.dumps({"engine": engine, "report": report, "ran": False, "adapt_only": True}, indent=2))
        return

    # Before anything is provisioned, submitted or run: no damaged .d gets searched here.
    # --adapt-only is above this on purpose -- it reads an engine's finished output and
    # never touches a raw file.
    refuse_damaged_tdf(files, a.allow_damaged_tdf)

    cmd = tools.get(engine)
    if not cmd:
        sys.exit(f"tools.json has no command for engine '{engine}'. "
                 f"Re-run acquire_tools.sh, or check its notes:\n  "
                 + "\n  ".join(tools.get("notes", [])))
    # Before anything is submitted, so a version the user did not confirm is on screen
    # while it can still be stopped -- not only in a file read after the search.
    ver_rec = engine_version_record(engine, tools, bundle)

    # A DIA-NN cfg that is not there is reported as exactly that, first -- not as whatever a
    # later reader makes of an empty flag list ("mass accuracy is not pinned").
    if engine == "diann" and not os.path.isfile(a.params):
        sys.exit(f"[run_search] cfg not found: {a.params}")

    # Do this BEFORE parallel_decision: that reads the cfg to check mass accuracy, and the
    # augmented copy is the cfg the run will actually use.
    if engine == "diann":
        a.params = ensure_xic(a.params, a.out)
        ensure_temp_dirs(a.params, a.out)
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
            "  Re-run with --sbatch <script>, then submit what it prints under "
            "\"submit with:\", e.g.:\n"
            f"    ... --sbatch ./{engine}_job.sh   # then: bash <out>/submit.sh, or sbatch "
            f"./{engine}_job.sh if that one script was written\n"
            "  (--allow-inline overrides this, e.g. inside an salloc/srun session.)")

    # These size the 5-step chain only. Said out loud when the route is single-shot, where
    # they change nothing: the flag being accepted reads as the jobs having got smaller.
    chain_only = [f for f, v in (("--libpred-cpus", a.libpred_cpus),
                                 ("--assembly-cpus", a.assembly_cpus),
                                 ("--assembly-mem", a.assembly_mem),
                                 ("--time-per-file", a.time_per_file),
                                 ("--max-simultaneous", a.max_simultaneous)) if v]
    if chain_only and not use_parallel:
        sys.stderr.write(f"[run_search] NOTE: {', '.join(chain_only)} size the 5-step chain "
                         f"only; this search is single-shot ({why}), so they are ignored. "
                         f"Its job(s) request --threads {a.threads} CPUs.\n")
    print(f"[run_search] engine={engine}  files={len(files)}  threads={a.threads}  "
          f"{'(5-step chain)' if use_parallel else '(emit sbatch)' if a.sbatch else '(inline)'}")
    sbatch_refused = None
    if use_parallel and a.sbatch:
        # --sbatch asks for ONE script to submit; the chain is six jobs chained by its own
        # submit.sh, so there is nothing to write there. A NOTE and exit 0 were not enough:
        # SKILL.md documents `run_search.py ... --sbatch job.sh && sbatch job.sh`, and when a
        # job.sh is left over from an earlier run -- say a sequential 310-file search -- the
        # `&& sbatch` resubmits THAT, silently. So: generate the chain, THEN move an existing
        # job script aside (renamed, never deleted), then exit non-zero so the `&&` stops.
        #
        # Anything but a regular file is refused before anything happens. `--sbatch proj`
        # renamed a whole project folder -- the cfg inside it with it -- and generation then
        # failed blaming mass accuracy.
        kind = _not_a_regular_file(a.sbatch)
        if kind:
            sys.exit(f"[run_search] REFUSING --sbatch {a.sbatch}: it is {kind}, not a job "
                     f"script, and nothing was changed. This search routed to the 5-step "
                     f"chain, which never writes --sbatch: drop the flag and submit "
                     f"{os.path.join(a.out, 'submit.sh')} once generated.")
        sbatch_refused = {"requested": os.path.abspath(a.sbatch), "written": False,
                          "existing_file_moved_to": None,
                          "exit_status": SBATCH_NOT_WRITTEN,
                          "why": "routed to the 5-step chain, which submits itself via submit.sh"}
    if use_parallel:
        res = run_diann_parallel(cmd, a.params, files, a.fasta, a.out, a.threads, a)
        if sbatch_refused:
            # Only now that the chain exists: a failed generation must leave the user's file
            # exactly where it was. run_diann_parallel exits on failure, so reaching here is
            # success.
            moved = set_aside(a.sbatch)
            sbatch_refused["existing_file_moved_to"] = moved and os.path.abspath(moved)
    elif engine == "diann":
        res = run_diann(cmd, a.params, files, a.fasta, a.out, a.threads, a.sbatch,
                        acquisition=bundle.get("acquisition", ""), queue=queue)
        generated = {os.path.abspath(p) for p in (res.get("submitted"),) if p} \
            if isinstance(res, dict) else set()
        generated |= {os.path.splitext(os.path.abspath(a.sbatch or "x"))[0] + s
                      for s in ("_1_lib.sh", "_2_search.sh")}
        if (a.sbatch and isinstance(res, dict) and res.get("mode") == "two_job_libfree"
                and os.path.abspath(a.sbatch) not in generated):
            # (--sbatch naming submit.sh or one of the two job scripts would otherwise have
            # its fresh file renamed to .stale-* -- review 2026-09-23.)
            # --sbatch names ONE script, and a library-free search is two jobs chained by
            # submit.sh, so that script is never written. One left over from an earlier search
            # would be what `sbatch job.sh` resubmits -- the hazard the chain route guards
            # against the same way: rename it, never delete it.
            try:
                moved = set_aside(a.sbatch)
            except ValueError as e:
                moved = None
                sys.stderr.write(f"[run_search] note: {e}\n")
            if moved:
                res["existing_sbatch_moved_to"] = os.path.abspath(moved)
                sys.stderr.write(f"[run_search] {a.sbatch} was left from an earlier search and "
                                 f"does not describe this one; moved to {moved}. Submit this "
                                 f"search with: bash {res['submitted']}\n")
    elif engine == "alphadia":
        res = run_alphadia(cmd, a.params, files, a.fasta, a.out, a.threads, a.sbatch,
                           queue=queue)
    elif engine == "sage":
        res = run_sage(cmd, a.params, files, a.fasta, a.out, a.threads, a.sbatch, queue=queue)
    elif engine == "fragpipe":
        res = run_fragpipe(cmd, bundle, a.params, files, a.fasta, a.out, a.threads, a.sbatch,
                           queue=queue)
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
                              library=a.library, mbr=a.mbr, queue=queue)
    else:
        sys.exit(f"unknown engine {engine}")

    # always record what was run (engine + version + exact command) for reproducibility
    try:
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, "search_provenance.json"), "w") as fh:
            # Where the FULLY-resolved parameters are, and WHEN they exist -- as the generator
            # reports it, not assumed here. On the probe path the chain measures the scan
            # window at run time (step 1b) and only then writes the file, so recording it as
            # already resolved would be false until step 1b succeeds.
            rp = res.get("resolved_params") if isinstance(res, dict) else None
            # `version` is what fran_deposit.py forwards to FRAN as the engine version, inside
            # fran_manifest.json -- THIS file is deliberately not staged into a drop entry,
            # because FRAN's scanner reads the name `search_provenance.json` as a **Radiant**
            # marker and would relabel every DIA-NN search (fran_deposit.py LINK_ITEMS).
            # `engine_version` is the record `version` comes from (engine_version_record:
            # only ever a version number, never "latest", never the manifest's pin).
            json.dump({"engine": engine, "version": ver_rec["value"],
                       "engine_version": ver_rec, "resolved_command": cmd,
                       "params_file": a.params,
                       "resolved_params_file": (rp or {}).get("file") or a.params,
                       "resolved_params_produced": ((rp or {}).get("produced")
                                                    or "before the search (the params file as given)"),
                       "resolved_params_note": (rp or {}).get("note"),
                       "scan_window": scan_window_record(engine, a.params, res),
                       "fasta": a.fasta, "threads": a.threads,
                       "n_files": len(files), "files": files,
                       "search_mode": "parallel_5step" if use_parallel else "single_shot",
                       "parallel_routing_reason": why,
                       # what was actually written to submit -- for a library-free search
                       # that is <out>/submit.sh, never the --sbatch name it did not write
                       "submitted_sbatch": None if sbatch_refused else (
                           (res.get("submitted") if isinstance(res, dict) else None)
                           or a.sbatch or None),
                       "sbatch_refused": sbatch_refused, "result": res}, fh, indent=2)
    except Exception as e:
        sys.stderr.write(f"[run_search] could not write search_provenance.json: {e}\n")

    print(json.dumps(res, indent=2))

    if sbatch_refused:
        moved = sbatch_refused["existing_file_moved_to"]
        sys.stderr.write(
            f"[run_search] --sbatch {a.sbatch} was NOT written: this search routed to the "
            f"5-step parallel chain, which is six SLURM jobs chained by its own submit.sh.\n"
            + (f"[run_search] The existing {a.sbatch} does not describe this search and was "
               f"moved to {moved}, so `sbatch {a.sbatch}` cannot resubmit it.\n" if moved else "")
            + f"[run_search] The chain IS generated. Submit it with: bash "
              f"{os.path.join(a.out, 'submit.sh')}  (or hive_exec.sh 'bash .../submit.sh')\n"
            f"[run_search] Exiting {SBATCH_NOT_WRITTEN} so `... --sbatch {a.sbatch} && sbatch "
            f"{a.sbatch}` stops here. For ONE job script instead, re-run with --no-parallel.\n")
        sys.exit(SBATCH_NOT_WRITTEN)


if __name__ == "__main__":
    main()
