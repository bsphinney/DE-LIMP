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
           DDA is searched with --dda, which estimate_params.py writes into the cfg;
           a cfg whose --dda disagrees with the bundle's acquisition is refused on both
           routes (diann_parallel.dda_mismatch). If inputs are Thermo .raw,
           auto-provisions a .NET 8 runtime (ensure_dotnet8.sh) so 2.6 can read them
  sage     convert .raw -> mzML (ThermoRawFileParser) / .d -> mzML (msconvert under Wine in the
           Core's pinned ProteoWizard container, --pwiz-sif) if needed -- inside the job with
           --sbatch, never before submission; a cfg with quant.lfq false (identifications
           only) hands a timsTOF .d to Sage as it is -- then
           <cmd> <bundle sage_config.json> -f <fasta> -o <out> --parquet
           --disable-telemetry-i-dont-want-to-improve-sage
           then sage_lfq_check.py (LFQ window vs the runs' MS1 mass error: warns; runs outside
           the window are a gate -- --adapt-only refuses report.parquet until the recorded
           corrected-window re-run is used or --accept-lfq-window; never re-runs by itself)
           then adapt lfq.parquet -> DIA-NN-shaped report for --method maxlfq, keeping only
           Sage's valid MS1 quantities (not decoy, q_value <= 0.05; counts in sage_adapt.json).
           Identifications only: no LFQ check, and results.sage.parquet -> identifications.tsv
           (targets at peptide_q <= 0.01), no report.parquet
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
import sys, os, json, glob, re, shlex, argparse, subprocess, shutil, stat, time, math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# ONE definition of the q-value columns -- see diann_q_columns.py.
from diann_q_columns import FDR_REQUIRED, PROTEIN_Q_PREFERENCE
# the one reading of a protein identifier
from protein_ids import fragpipe_group, group_accessions, group_entry_names
# No Bruker .d is searched before its analysis.tdf has been checked -- see bruker_tdf.py.
from bruker_tdf import (STATUSES, INTEGRITY_RULES, fingerprint, integrity_warning,
                        tdf_integrity)
# Every file a job runs or reads is written without following a planted link: see safe_write.py.
from safe_write import open_text, write_file
# An identifications-only Sage search's output file: named where its LFQ record says so.
from sage_lfq_check import IDENTIFICATIONS as SAGE_IDS

ALLOW_DAMAGED_TDF = "--allow-damaged-tdf"

# The minimum an adapted report must carry for the DE step to run. The q-columns
# are exactly the REQUIRED filter set: an adapter that emits fewer leaves limpa
# silently applying no filter for the missing ones.
DIANN_CONTRACT = ["Run", "Protein.Group", "PG.MaxLFQ"] + FDR_REQUIRED


def sh(cmd, **kw):
    print(f"  $ {cmd}", flush=True)
    return subprocess.run(cmd, shell=True, check=True, **kw)


# More .d than this still to read (no matching record from step 2) on a login node, and the
# whole command moves to a compute node (on_compute.ensure_compute_node) -- HPC@UCD's rules: no
# I/O-heavy work on a login node. Five is a few seconds of reads.
TDF_LOGIN_MAX = 5
# Set in the environment of the run_search.py that ensure_compute_node() re-runs under srun for
# the tdf check: that srun is a short, small allocation, so a search never runs inline in it.
TDF_REEXEC_ENV = "RUN_SEARCH_MOVED_FOR_TDF_CHECK"


def tdf_files(files):
    """The inputs that are a Bruker .d with an analysis.tdf: the ones the gate reads."""
    return [f for f in files
            if f.rstrip("/").lower().endswith(".d")
            and os.path.isfile(os.path.join(f.rstrip("/"), "analysis.tdf"))]


def find_acquisition_json(out, given=None):
    """step 2's detect_acquisition.py output: `given` (--acquisition-json), else the first
    input/acquisition.json walking up from --out, at most 3 folders (a session's
    output/search -> <session>/input/acquisition.json). None when there is none."""
    if given:
        return given
    d = os.path.abspath(out)
    for _ in range(4):
        cand = os.path.join(d, "input", "acquisition.json")
        if os.path.isfile(cand):
            return cand
        up = os.path.dirname(d)
        if up == d:
            break
        d = up
    return None


def recorded_integrity(files, acq_json):
    """{path: tdf_integrity} for each input .d whose record in `acq_json` (step 2) still
    describes it: made under today's rules (bruker_tdf.INTEGRITY_RULES), its fingerprint
    (bruker_tdf.fingerprint: sizes, mtimes, ctimes, inodes) equals the files' now, and so do the
    sizes of any -wal/-journal beside the tdf -- a side file that appeared since is exactly what
    the check is for. A record without a fingerprint or rules (before 2.14.2, or other rules), an
    unreadable or missing acquisition.json: not recorded, so that .d is read again (never a
    crash)."""
    if not acq_json:
        return {}
    try:
        with open(acq_json) as fh:
            doc = json.load(fh)
        recs = doc.get("files") or []
    except (OSError, ValueError, AttributeError) as e:
        sys.stderr.write(f"[run_search] NOTE: {acq_json} could not be read ({type(e).__name__}: "
                         f"{e}); every .d's analysis.tdf is checked here instead.\n")
        return {}
    by_path = {}
    for r in recs if isinstance(recs, list) else []:
        if isinstance(r, dict) and isinstance(r.get("file"), str) \
                and isinstance(r.get("tdf_integrity"), dict):
            by_path[os.path.realpath(r["file"].rstrip("/"))] = r["tdf_integrity"]
    out = {}
    for f in tdf_files(files):
        p = f.rstrip("/")
        integ = by_path.get(os.path.realpath(p))
        if (not integ or integ.get("rules") != INTEGRITY_RULES or not integ.get("fingerprint")
                or integ["fingerprint"] != fingerprint(p)):
            continue
        tdf = os.path.join(p, "analysis.tdf")
        side = {"wal_bytes": tdf + "-wal", "journal_bytes": tdf + "-journal"}
        if any((os.path.getsize(v) if os.path.isfile(v) else 0) != (integ.get(k) or 0)
               for k, v in side.items()):
            continue
        out[f] = integ
    return out


def tdf_check_minutes(n):
    """The srun's time for reading n analysis.tdf: about 4 s each, plus 5 minutes, at least 10."""
    return max(10, 5 + math.ceil(n * 4 / 60))


def check_tdf(files, out, acquisition_json, allow, queue=None):
    """The tdf gate of main(): step 2's records reused where they still describe the file, the
    rest read here -- or, more than TDF_LOGIN_MAX of them on a login node, this whole command
    re-run on a compute node first (never returns then), on the --partition / --account / --qos
    it was given (`queue`), else on_compute's choice; the check only reads, so a preemptible
    queue may take it. Refuses as refuse_damaged_tdf()."""
    acq = find_acquisition_json(out, acquisition_json)
    recorded = recorded_integrity(files, acq)
    unchecked = [f for f in tdf_files(files) if f not in recorded]
    if recorded:
        print(f"[run_search] analysis.tdf: {len(recorded)} .d checked by step 2 and unchanged "
              f"since ({acq}); {len(unchecked)} read here", flush=True)
    import on_compute
    if len(unchecked) > TDF_LOGIN_MAX and on_compute.on_login_node():
        os.environ[TDF_REEXEC_ENV] = "1"
        q = queue or {}
        on_compute.ensure_compute_node(
            f"run_search.py reads the analysis.tdf of {len(unchecked)} .d",
            minutes=tdf_check_minutes(len(unchecked)), mem="2G", script=__file__,
            preemptible_ok=True, partition=q.get("partition"), account=q.get("account"),
            qos=q.get("qos"))
        os.environ.pop(TDF_REEXEC_ENV, None)    # not moved (main() called in-process): read here
    return refuse_damaged_tdf(files, allow, recorded)


def tdf_problems(files, recorded=None):
    """Every Bruker .d in `files` whose analysis.tdf is not `ok`, worst first. A .d in
    `recorded` (recorded_integrity: step 2's check, still describing the file) is not read
    again; its recorded status counts exactly as a fresh one.

    detect_acquisition.py runs the same check, but step 2 is exactly what the entry points
    SKILL.md advertises for a bare search skip: "re-run this search", "re-search with
    different parameters", "the paths are already known" all arrive here with the files in
    hand and no detection pass behind them. A .d whose frame index covers part of its
    tdf_bin is then searched with no error anywhere -- DIA-NN reads what the index points
    at, and one HIVE run's index covered 0.7% of a 2.4 GB tdf_bin. So the gate has to be on
    THIS path too, not only on the one the orchestrator is told to walk.

    Per file this is a 100-byte header read, one indexed query and a 4-byte read -- nothing
    that needs a compute node, and nothing next to a multi-hour search."""
    recorded = recorded or {}
    out = []
    for f in tdf_files(files):
        p = f.rstrip("/")
        r = recorded.get(f) or tdf_integrity(p)
        if r["status"] != "ok":
            out.append((p, r))
    out.sort(key=lambda t: (STATUSES.index(t[1]["status"]), t[0]))
    return out


def refuse_damaged_tdf(files, allow, recorded=None):
    """Stop before the search when any input .d is not `ok`. Prints and returns the list."""
    bad = tdf_problems(files, recorded)
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
            "  DIA-NN 2.6 will not read .raw. Options: re-run ensure_dotnet8.sh (needs\n"
            "  internet; on a cluster login node it moves to a compute node), or feed mzML.\n")
        return ""
    return (f"export DOTNET_ROOT={shlex.quote(root)}; "
            f'export PATH={shlex.quote(root)}:"$PATH";')


# Per-user CPU cap on HIVE's genome-center-grp/high (docs/QUEUE_SWITCHING.md). The fallback
# only: user_limits() reads the QOS's own MaxTRESPU from sacctmgr first (cpu=64,mem=1T on
# genome-center-grp-high-qos; publicgrp-low-qos has none -- `sacctmgr show qos`, 2026-10-01).
HIVE_USER_CPU_CAP = int(os.environ.get('HIVE_USER_CPU_CAP', '64'))
# The fewest CPUs array_task_cpus() gives an array task it sizes itself. MEASURED (HIVE,
# 2026-10-01, jobs 24276579-82): the chain's step-2 first pass of one 50 ng HeLa QC run
# (timsTOF HT diaPASEF, 60 SPD, 1.9 GB .d; the Core's human 2026-09 predicted library, 3.5 M
# precursors; DIA-NN 2.7.0), each thread count alone on an AMD EPYC 7532 node (zen2, the
# commonest CPU on `high`), identical IDs at every count (43,123 precursors, 5,287 groups):
#   threads   4      8      16     32
#   wall      39.3   21.4   10.9   6.7 min    (doubling: 1.84x, 1.96x, 1.62x)
#   CPU x s   9,437  10,283 10,469 12,957     (what one file costs a per-user cap)
#   peak RSS  7.4    7.9    8.7    9.7 GB
# Under the 64-CPU cap a 29-file array then takes 79 / 86 / 87 / 101 min at 4 / 8 / 16 / 32
# per file. 8 is the floor: within 9% of the cheapest per file (4) at half its wall time, and
# 32 per file -- what --threads 32 used to give every task -- costs 26% more per file.
# references/diann_parallel.md ("CPUs per array task") has the table.
ARRAY_TASK_MIN_CPUS = 8


_SBATCH_OFF_PATH_WARNED = False


def sbatch_path():
    """sbatch, or None when there is no SLURM here: on PATH, else at HPC@UCD's fixed path
    (on_compute: HPCCF_SLURM_BIN, or SLURM_BIN_DEFAULT) when it is executable there. A
    non-login shell on HIVE -- an agent's `ssh hive 'cmd'`, or a shell after `module purge` --
    has no sbatch on PATH, and `which` alone then routed a 6-file cohort to one single-shot
    search and let an inline search start on the login node (I13, I44). Said once, on stderr,
    when only the fixed path has it: submit.sh calls `sbatch` by name."""
    global _SBATCH_OFF_PATH_WARNED
    found = shutil.which("sbatch")
    if found:
        return found
    import on_compute
    d = on_compute._slurm_bin(os.environ)
    p = os.path.join(d, "sbatch")
    if not os.access(p, os.X_OK):
        return None
    if not _SBATCH_OFF_PATH_WARNED:
        _SBATCH_OFF_PATH_WARNED = True
        out = getattr(a_globals, "out", None)
        submit = shlex.quote(os.path.join(out, "submit.sh")) if out else "<out>/submit.sh"
        sys.stderr.write(f"[run_search] WARNING: sbatch is not on PATH in this shell (a non-login "
                         f"shell); found {p}. Submit with: PATH={shlex.quote(d)}:$PATH bash "
                         f"{submit}\n")
    return p


def slurm_available():
    return sbatch_path() is not None


def cpu_hours_mod():
    """cpu_hours.py: THE reserved CPU-hours of a search (the most it can use)."""
    import cpu_hours
    return cpu_hours


def search_cpu_hours(res):
    """`res` with cpu_hours_reserved for a search that wrote ONE job script (every --sbatch route
    but the two-job library-free search, the 5-step chain and the Radiant chain, which say it
    themselves). Unchanged for an inline search: nothing was reserved."""
    if (isinstance(res, dict) and "cpu_hours_reserved" not in res
            and isinstance(res.get("submitted"), str) and os.path.isfile(res["submitted"])):
        cpu = cpu_hours_mod().reservation([res["submitted"]])
        res.update(cpu_hours_reserved=cpu["cpu_hours_reserved"],
                   gpu_hours_reserved=cpu["gpu_hours_reserved"], cpu_hours_reservation=cpu)
    return res


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
        import on_compute
        return False, (f"{n} files, but no SLURM here (no sbatch on PATH or in "
                       f"{on_compute._slurm_bin(os.environ)}) -- the 5-step chain runs as job "
                       "arrays, so a single search is the only option")
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
    """diann_parallel._not_a_regular_file -- the one definition."""
    return _diann_parallel_mod()._not_a_regular_file(path)


def set_aside(path):
    """diann_parallel.set_aside -- the one definition: rename a regular file `.stale-<time>`,
    never delete it; anything else raises ValueError."""
    return _diann_parallel_mod().set_aside(path)


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
        return {"mode": "unknown", "source": f"unknown -- {e}", "value": None, "passed": None}


def mass_acc_provenance(engine, params, res):
    """What set DIA-NN's mass accuracy, for search_provenance.json `mass_acc` -- present for every
    DIA-NN search, with its stable `mode` (diann_parallel.MASS_ACC_MODES), as scan_window is: the
    search's own record (result.mass_acc) when it has one, else read from the cfg it ran."""
    if isinstance(res, dict) and isinstance(res.get("mass_acc"), dict):
        return res["mass_acc"]
    if engine != "diann":
        return None
    dp = _diann_parallel_mod()
    try:
        return dp.mass_acc_record(dp.mass_acc_status(params), params)
    except dp.CfgError as e:
        return {"mode": "unknown", "fixed": None, "reason": f"unknown -- {e}"}


def protein_inference_provenance(engine, params):
    """DIA-NN's protein inference as the cfg set it, for search_provenance.json
    `protein_inference` (diann_parallel.protein_inference_record: the one reading, which the
    Methods use too). None for another engine."""
    if engine != "diann":
        return None
    return _diann_parallel_mod().protein_inference_record(params)


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
    write_file(aug, txt.rstrip("\n") + f"\n{extra}\n")
    # The copy is the same parameters plus XICs, so it keeps the cfg's estimate_params.py
    # rationale: without it the mass-accuracy plan (measure_with_diann) would be lost here and the
    # search would fall back to DIA-NN's first-run auto mode (diann_parallel.mass_acc_measure_plan).
    # A cfg with NO sidecar must not inherit one an earlier search left in the same --out: that
    # stale plan would have a hand-written cfg's missing mass accuracy measured instead of reported.
    if os.path.exists(params + ".rationale.json"):
        with open(params + ".rationale.json", encoding="utf-8") as fh:
            write_file(aug + ".rationale.json", fh.read())
    elif os.path.lexists(aug + ".rationale.json"):
        os.remove(aug + ".rationale.json")
    print(f"[run_search] cfg had no --xic; using {aug} (added: {extra}). "
          f"Every DIA-NN search extracts chromatograms.")
    return aug


def refuse_dda_mismatch(params, acquisition):
    """Stop before anything is generated when the cfg's --dda disagrees with the bundle's
    acquisition -- the same check on both DIA-NN routes (diann_parallel.dda_mismatch)."""
    dp = _diann_parallel_mod()
    try:
        why = dp.dda_mismatch(params, acquisition)
    except dp.CfgError as e:
        if e.code == "cfg_missing":
            return          # each route reports a missing cfg its own way, as before
        sys.exit(f"[run_search] {e}")
    if why:
        sys.exit(f"[run_search] REFUSING the DIA-NN search: {why} Nothing was generated.")


def run_diann_parallel(cmd, params, files, fasta, out, threads, a, acquisition=""):
    """Generate DIA-NN's 5-step SLURM chain (does not submit -- the orchestrator does,
    then watches the step-5 job with watch_run.sh)."""
    # Before the generator: its own failure message says "re-run with --no-parallel", and the
    # single-shot route refuses the same cfg.
    refuse_dda_mismatch(params, acquisition)
    os.makedirs(out, exist_ok=True)
    listing = os.path.join(out, "parallel_input_files.txt")
    write_file(listing, "\n".join(files) + "\n")     # a list file survives spaces in paths
    argv = [sys.executable,
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "diann_parallel.py"),
            "--diann", cmd, "--raw-list", listing, "--fasta", fasta,
            "--out", out, "--cfg", params, "--threads-max", str(threads)]
    # Chain resources are forwarded only when given, so diann_parallel.py's own defaults stay
    # the single definition of them. They are exposed because steps 3 and 5 ask for 64 CPUs
    # by default (--assembly-cpus), a whole node's worth that can sit pending for a long time
    # on a busy preemptible queue, and until now the only way to ask for less was to bypass
    # run_search.py and hand-run diann_parallel.py (FRAN pilot, publicgrp/low, 2026-09-16).
    for flag, val in (("--partition", a.partition), ("--account", a.account),
                      ("--qos", a.qos), ("--max-simultaneous", a.max_simultaneous),
                      ("--exclude", _exclude_arg(a)),
                      # pins the array tasks' CPUs; without it diann_parallel.py sizes them to
                      # the queue's per-user cap, at most --threads (array_task_cpus)
                      ("--threads-per-file", getattr(a, "threads_per_file", None)),
                      ("--libpred-cpus", getattr(a, "libpred_cpus", None)),
                      ("--assembly-cpus", getattr(a, "assembly_cpus", None)),
                      ("--assembly-mem", getattr(a, "assembly_mem", None)),
                      ("--time-per-file", getattr(a, "time_per_file", None)),
                      ("--mem-per-file", getattr(a, "mem_per_file", None))):
        if val:
            argv += [flag, str(val)]
    argv += _job_end_argv(a)
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


def single_shot_mass_acc(params, listing, fasta, predicted, out, threads, cmd, libfree,
                         provenance=None):
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
    .tmp and moved into place only once the value is measured -- or, as in step 1b, once a
    probe that measured nothing twice has fallen back (diann_parallel.probe_attempts(),
    probe_fallback.py): the documented level as given, the other at the SOP, tagged DEFAULT,
    recorded in <out>/probe_fallback.json and, when `provenance` is given (the --sbatch routes,
    whose search_provenance.json exists before the job runs), in it.

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
    # A DDA cfg that still plans a measurement: refused now, with the chain's own words, rather
    # than carried into the job for probe_window.py to refuse at run time.
    refusal = dp.mass_acc_dda_refusal(params)
    if refusal:
        sys.exit(f"[run_diann] REFUSING the search: {refusal['reason']}. Fix: "
                 f"{refusal['remedy']}. Nothing was generated.")
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
        return [], "", {"mode": "auto", "fixed": False, "measured": False, "reason": why}, None
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
    fallback_out = os.path.join(out, dp.FALLBACK_RECORD)
    # A probe whose own machinery failed no longer stops the search: probe_attempts() retries
    # it (a crash, an unreadable log) and falls back (probe_fallback.py: the documented level as
    # given, the other at the SOP, tagged DEFAULT). Anything else stops it before DIA-NN runs.
    failed = ("{ echo \"FAILED: no mass accuracy measured, and no fallback (only a failure of "
              "the probe's own machinery may fall back):\" >&2; "
              + " ".join(ln.strip() for ln in dp.probe_failure_lines()) + ";"
              + f' echo "Per-run evidence in {evidence}. Nothing was searched." >&2; '
              f'rm -f {q(tmp)} {q(massacc)}; exit 1; }}')
    lines = [
        'echo "measuring mass accuracy with DIA-NN on representative runs before the search"',
        # an earlier run's value, evidence, fallback or resolved cfg must not survive
        f"rm -f {q(evidence)} {q(evidence + '.attempt1')} {q(massacc)} {q(fallback_out)} "
        f"{q(resolved)} {q(tmp)}",
        f"rm -rf {q(workdir)} {q(workdir + '.attempt1')}",
        f"cp {q(params)} {q(tmp)}",
        *dp.probe_attempts(
            f"python3 {q(probe)} --diann {q(cmd)} --raw-list {q(listing)} --fasta {q(fasta)} "
            f"--lib {q(predicted)} --threads {threads} --measure mass-acc{doc_args} "
            f"--max-probes {dp.PROBE_CANDIDATES} --max-failures {dp.PROBE_MAX_FAILURES} "
            f"--timeout {dp.PROBE_TIMEOUT_S} --workdir {q(workdir)} --write-cfg {q(tmp)}",
            flags, evidence, ["mass-acc"], documented,
            reset=dp.keep_attempt(workdir) + [f"cp {q(params)} {q(tmp)}"],
            massacc_file=massacc, write_cfg=tmp, provenance=provenance,
            fallback_out=fallback_out),
        'if [ "$PROBE_RC" -eq 0 ]; then',
        # the two flags, in exactly the shape the guard below accepts
        "  python3 -c \"import json,re,sys; m=json.load(open(sys.argv[1]))['mass_acc']['pin_as']; "
        f"assert re.fullmatch(sys.argv[2], m); print(m)\" {q(evidence)} {q(pattern)} "
        f"> {q(massacc)} || {failed}",
        # probe_fallback.py wrote massacc.txt and the cfg's two flags
        f'elif [ "$PROBE_FALLBACK" -ne 1 ]; then {failed}; fi',
        dp.needs_measured(massacc, "the two mass-accuracy flags",
                          producer="the probe before this search"),
        f"mv -f {q(tmp)} {q(resolved)}",
        'if [ "$PROBE_FALLBACK" -eq 1 ]; then',
        f'  echo "mass accuracy = $(cat {q(massacc)}) -- FALLBACK: the probe measured nothing, '
        'so the documented level as given and the other at the facility SOP (DEFAULT, not '
        f'measured); record in {fallback_out}"',
        "else",
        f'  echo "mass accuracy = $(cat {q(massacc)}) (pinned for the search; evidence '
        f'{evidence})"',
        "fi",
    ]
    record = {"mode": "measured", "fixed": True, "measured": True, "documented": documented,
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
                       "note": "written only after mass accuracy is measured -- or set by the "
                               "fallback when the probe measured nothing (probe_fallback.json) -- "
                               "just before the DIA-NN search; absent until then, so a missing "
                               "file after the search job means the probe was refused or "
                               "stopped"}
    return lines, f" $(cat {q(massacc)})", record, resolved_params


def run_diann(cmd, params, files, fasta, out, threads, sbatch, acquisition="", queue=None,
              raw_files=None, keratin_sample=None):
    """Single-shot DIA-NN. `queue` is {partition, account, qos} from the command line, handed
    to every emit_sbatch() call -- without it slurm_queue() re-detects a queue and overrides
    the one the user chose. `raw_files`: the inputs as given, when some of `files` are links to
    them (diann_parallel.dashed_input_plan) -- the job's node check looks at their folders."""
    queue = queue or {}
    # `--cfg <params>` goes on the job's command line verbatim, and the job runs after
    # `cd <out>`: a relative path would be read there -- another file, or none -- and DIA-NN logs
    # it as given, so the run's own record could not say which file it read. Absolute, and it
    # must exist now, before anything is written (main() checks too; this writes the command).
    params = os.path.abspath(params)
    if not os.path.isfile(params):
        sys.exit(f"[run_diann] cfg not found: {params}")
    refuse_dda_mismatch(params, acquisition)
    # SILAC (--channels): refused when it would label twice or stop on channel normalisation
    import silac
    why = silac.refusal(params, prog="run_diann")
    if why:
        sys.exit(why)
    os.makedirs(out, exist_ok=True)
    # The job's first lines check that its node reaches DIA-NN, --out, the FASTA and the raw
    # data (node_fault.py): a dropped mount is a node fault (exit 75), never "DIA-NN failed".
    if sbatch:
        import node_fault
        queue = dict(queue, preflight=node_fault.preflight_lines(
            node_fault.chain_checks(cmd, os.path.abspath(out), fasta, raw_files or files)))
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
    write_file(listing, "\n".join(files) + "\n")     # a list file survives spaces in paths
    # DIA-NN 2.6 supports DDA via --dda (must NOT be used on DIA data). QuantUMS is
    # auto-disabled on DDA; for DDA quant DIA-NN recommends extra MS1 filtering on
    # Ms1.Global.Q.Value / Ms1.Global.Quality (see references/search-engines.md).
    # The flag comes from the cfg (estimate_params.py writes it), as it does for the chain. This
    # route used to append it to the command instead -- so the chain never had it, and a cfg that
    # already carried it got it twice. A cfg that disagrees with the bundle is refused (at the
    # top of this function, before anything is written).

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
        sys.exit(f"[run_diann] {e}")
    present = {f for f, _ in groups}
    dda = dp.DIANN_DDA_FLAG in present          # recorded in the result; the cfg carries it
    libfree = {"--fasta-search", "--gen-spec-lib"} <= present \
        and not getattr(a_globals, "one_step", False)
    dnet = dotnet_env_for(files)          # .NET 8 for reading Thermo .raw, if needed

    # With an extension, as the chain's step 1 names its own: DIA-NN cuts the name it derives at
    # the last '.' of the whole path, so a bare `diann_lib` in a folder like `P2_X.01dup2` put the
    # library and its log one folder up (dp.diann_output_path). It ignores the extension asked of
    # --out-lib for a PREDICTED library and always writes <name>.predicted.speclib (README,
    # "Output library") -- the file the search job reads.
    lib = os.path.join(out, "diann_lib.speclib")
    predicted = dp.predicted_speclib(lib)
    search_cfg = dp.bash_flags(groups, drop=SINGLE_SHOT_SEARCH_STRIP)
    lib_cmd = (f"{cmd} --cfg {shlex.quote(params)} --fasta {shlex.quote(fasta)} "
               f"--out-lib {shlex.quote(lib)} --threads {threads}")
    # An undocumented Orbitrap mass accuracy is measured between the library and the search
    # (single_shot_mass_acc); `mflag` then carries the pinned flags into the search.
    # the --sbatch routes write search_provenance.json now, before the job runs, so a fallback
    # in the job records itself there; inline, it is written after the search, from ma_rec
    measure_lines, mflag, mass_acc, resolved_params = single_shot_mass_acc(
        params, listing, fasta, predicted, out, threads, cmd, libfree,
        provenance=os.path.join(out, "search_provenance.json") if sbatch else None)
    # recorded only when there is something to say, so a pinned cfg's result is unchanged
    ma_rec = {} if mass_acc is None else {"mass_acc": mass_acc}
    if mass_acc is None and groups:
        # ...and when the cfg pins a DDA level at the SOP, say it is a DEFAULT (rule 2)
        st = dp.mass_acc_status(params)
        if dp.mass_acc_defaults(params, st):
            ma_rec["mass_acc"] = dp.mass_acc_record(st, params)
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
                  f"--lib {shlex.quote(predicted)} --reanalyse --matrices "
                  f"--out {shlex.quote(report)} --threads {threads}{tmp_arg}")
    onecmd = (f"{cmd} --cfg {shlex.quote(params)} {f_args} "
              f"--fasta {shlex.quote(fasta)} --out {shlex.quote(report)} "
              f"--threads {threads}{tmp_arg}")
    # DIA-NN exits 0 on fatal errors, so each job asserts the artefact it exists to make --
    # the same contract every step of the 5-step chain has (references/diann_parallel.md).
    # And each job first DELETES that artefact, because an existence check cannot tell this
    # run's file from the previous search's in the same --out. Reproduced on HIVE, DIA-NN
    # 2.7.0 (review srun 23512013): a re-run with no DOTNET_ROOT logged "ERROR: cannot read
    # .raw files", exited 0 and wrote nothing; the old report.parquet and report.stats.tsv
    # were byte-identical afterwards, and both guards passed on them. The stats file goes too:
    # it is check_report_runs.py's fallback evidence when there is no parquet reader.
    import check_report_runs
    import contaminant_competition
    import dda_first_pass
    # ...and the contaminant record, which must never describe a previous search's report
    report_files = (report, check_report_runs.stats_path(report),
                    os.path.join(out, contaminant_competition.RECORD))
    if dda:
        # A --dda search whose second pass crashes falls back to the first-pass report DIA-NN
        # wrote (dda_first_pass.py), so an earlier search's first pass and record must not
        # survive into this one: with them, a crash before this first pass would use THAT one.
        report_files += tuple(dda_first_pass.stale_files(report))
    provenance = os.path.join(out, "search_provenance.json")

    def searched(command):
        """The DIA-NN search line; for --dda, wrapped in the first-pass fallback (a DIA search's
        job is unchanged)."""
        return dda_first_pass.job_lines(command, report, provenance) if dda else [command]
    lib_job = "\n".join([dp.clear_stale(predicted), lib_cmd,
                         dp.must_exist(predicted, "the predicted spectral library")])
    guard = report_guard(report, listing)
    # then which sample proteins the contaminant list took away: a report, never a gate
    competition = contaminant_competition.job_lines(report, fasta, out)
    # ...and, for --dda, an earlier run's fallback leaves the record (dda_first_pass.forget):
    # a resubmitted job whose second pass succeeds must not read as a first pass (L8R-F1)
    cleared = [dp.clear_stale(*report_files)] + (dda_first_pass.forget_lines(provenance)
                                                 if dda else [])
    search_job = "\n".join([*cleared, *measure_lines, *searched(search_cmd), guard, *competition])
    one_job = "\n".join([*cleared, *searched(onecmd), guard, *competition])

    # Every refusal is behind us (single_shot_mass_acc's is the last): an earlier search's probe
    # outputs must not describe this one -- set aside, never deleted (dda-review N1, R1). Inline,
    # they go with the report files in the loop below, after its own refusal.
    if sbatch:
        dp.set_aside_probe_outputs(out)
        contaminant_competition.set_aside_record(out)    # and an earlier search's record
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
                    submit_hint=False, notify="fail",
                    stage="library prediction (job 1 of 2)", **queue)
        # The measurement belongs to the SEARCH job: a search requeued against the same library
        # measures again rather than trusting a massacc.txt from a run it cannot vouch for.
        emit_sbatch(srch_sh, search_job, out, threads, job="diann_search", preamble=dnet,
                    hours=search_job_hours(measure_lines), submit_hint=False,
                    stage="search (job 2 of 2)", fran_guarded=True, **queue)
        submit = os.path.join(os.path.abspath(out), "submit.sh")
        jobs_txt = os.path.join(os.path.abspath(out), "jobs.txt")
        # jobs.txt, as the 5-step chain's submit.sh writes it: `watch_run.sh --all <out>`
        # reads it, and without one it answered failed/no_jobs_file for a healthy running
        # search -- which step 7b says to resubmit (HIVE e2e test 2026-09-23).
        cpu = cpu_hours_mod().reservation([lib_sh, srch_sh])
        with open_text(submit, 0o755) as fh:
            fh.write("#!/bin/bash -l\n" + cpu_hours_mod().submit_comment(cpu)
                     + "\nset -euo pipefail\n"
                     f"j1=$(sbatch --parsable {shlex.quote(lib_sh)})\n"
                     f"j2=$(sbatch --parsable --dependency=afterok:$j1 {shlex.quote(srch_sh)})\n"
                     f'printf "%s\\n" "$j1" "$j2" > {shlex.quote(jobs_txt)}\n'
                     'echo "submitted: libpred=$j1 search=$j2"\n'
                     f'echo "both job ids -> {jobs_txt}  (watch both: watch_run.sh --all '
                     f'{os.path.abspath(out)})"\n'
                     f'echo "report will be {report}"\n')
        print(f"  [sbatch] two-job chain: {lib_sh} -> {srch_sh}; submit with: bash {submit}")
        print(f"  [sbatch] --sbatch {sbatch} itself is NOT written for a library-free search: "
              f"submit.sh submits both jobs.")
        return {"engine": "diann", "report": report, "submitted": submit,
                "mode": "two_job_libfree", "library": predicted,
                "ran": False, "dda": bool(dda), "raw_dotnet": bool(dnet),
                "cpu_hours_reserved": cpu["cpu_hours_reserved"],
                "gpu_hours_reserved": cpu["gpu_hours_reserved"], "cpu_hours_reservation": cpu,
                **ma_rec}

    if sbatch:                          # libfree + sbatch returned above, so this is onecmd
        emit_sbatch(sbatch, one_job, out, threads, job="diann_search", preamble=dnet,
                    hours=search_job_hours(measure_lines), fran_guarded=True, **queue)
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
    # ...and their probe outputs with them: kept as that search's evidence, never this one's
    dp.set_aside_probe_outputs(out)
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
            # the probe measured nothing and fell back: the record says so, not "measured"
            fb = os.path.join(out, dp.FALLBACK_RECORD)
            if os.path.isfile(fb):
                with open(fb) as fh:
                    rec = json.load(fh)
                ma_rec["mass_acc"] = dict(rec["mass_acc"], value_file=mass_acc["value_file"],
                                          evidence_file=mass_acc["evidence_file"])
                ma_rec["probe_fallback"] = rec
    command = pre + (search_cmd if libfree else onecmd)
    fell_back = {}
    if dda:
        fell_back = dda_search_inline(command, report)
    else:
        sh(command)
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
            "previous_artefacts_moved_to": moved, **ma_rec, **fell_back,
            "contaminant_competition": competition_inline(report, fasta, out, keratin_sample)}


def dda_search_inline(command, report):
    """An inline --dda search, as the job runs it (dda_first_pass.job_lines): a missing or empty
    report beside a first pass makes the first pass the report, recorded in the result, which
    main() writes into search_provenance.json (`report_used`, `dda_first_pass`). Otherwise
    DIA-NN's non-zero status stops the search. -> {} or {report_used, dda_first_pass}."""
    import dda_first_pass
    print(f"  $ {command}", flush=True)
    rc = subprocess.run(command, shell=True).returncode
    fp = dda_first_pass.first_pass_path(report)
    nonempty = dda_first_pass.nonempty
    if nonempty(report) or not nonempty(fp):
        if rc != 0:
            sys.exit(f"DIA-NN search failed (exit {rc}) -- check its output above for ERROR:"
                     + ("" if nonempty(report) else f"; there is no first-pass report ({fp}) "
                                                    "to fall back to"))
        return {}
    try:
        rec = dda_first_pass.record(report, fp, None, rc)
    except (OSError, ValueError) as e:
        sys.exit(f"DIA-NN exited {rc} without {report}, and its first-pass report could not be "
                 f"used: {e}")
    rec.pop("provenance_updated", None)
    print(f"  [run_diann] {dda_first_pass.MARKER}: DIA-NN exited {rc} without {report}; its "
          f"first-pass report {fp} (no match-between-runs) is used instead.\n"
          f"  [run_diann] CAUTION: {rec['caveat']}", flush=True)
    return {"report_used": rec["report_used"], "dda_first_pass": rec}


def search_records(res):
    """What a search's result adds to search_provenance.json beside `result`: the records some
    routes make while the search runs (sage_lfq_check, sage_adapt, probe_fallback, cpu_sizing), and
    an inline DIA-NN search's contaminant record -- only a record with a status, never the
    chain's plan of one (diann_parallel's {file, produced: runtime}), so the key has one shape (a
    job merges its own in here)."""
    if not isinstance(res, dict):
        return {}
    out = {k: res[k] for k in ("sage_lfq_check", "sage_adapt", "probe_fallback", "cpu_sizing",
                               "sage_input", "sage_identifications", "pwiz_sif_arg",
                               # FragPipe's Java heap in its job (fragpipe_ram_gb), its DDA
                               # stage in the search folder and the mono its job reads .raw with
                               "fragpipe_ram_gb", "fragpipe_stage", "fragpipe_mono",
                               # the database it searched and the workflow it ran (FU3)
                               "fragpipe_database", "fragpipe_workflow",
                               # what its workflow quantifies (Ruling 73a)
                               "quantification", "identifications",
                               # an inline --dda search that fell back to its first pass (a
                               # job's dda_first_pass.py writes these here itself)
                               "report_used", "dda_first_pass")
           if res.get(k)}
    if isinstance(res.get("quantified"), bool):     # a Sage search: false = identifications only
        out["quantified"] = res["quantified"]
    cc = res.get("contaminant_competition")
    if isinstance(cc, dict) and cc.get("status"):
        out["contaminant_competition"] = cc
    return out


def competition_inline(report, fasta, out, keratin_sample=None):
    """An inline search's contaminant competition check, as the jobs run it
    (contaminant_competition.job_lines): <out>/contaminant_competition.json, and the record
    returned for search_provenance.json, which main() writes after the search -- so this
    search's keratin answer is handed in, never read from an earlier search's record there.
    Never stops the search: a check that fails is said, and recorded as not run."""
    import contaminant_competition as cc
    # check() and write() never raise: a check that could not run is a not_run record
    rec = cc.check(report, fasta, cc.resolve_inputs(None, report, fasta)[2],
                   keratin_sample=keratin_sample)
    _written, problems = cc.write(rec, os.path.join(out, cc.RECORD))
    cc.say(rec)
    for p in problems:
        sys.stderr.write(f"{cc.TAG} NOTE: could not write {p}; the record is in "
                         "search_provenance.json only.\n")
    return rec


# ------------------------------------------------- quantity declarations -----
# Every adapter says what its report's numbers ARE (architectural rule 1), as parquet schema
# metadata under QUANTITY_META; build_maxlfq.R declared_quantity() reads it and
# maxlfq_descriptor() -- the one description the methods, the report, the AI brief and the
# reproducibility log all read -- is built from it. No downstream file names an engine. Keys:
#   level           "protein" (one value per protein and run, the engine's own) | "peptide"
#   label           the short name the DE's display label uses
#   value           what the number is, with engine version and source file
#   identification  the FDR the engine applied (the report's q-columns may be placeholders)
#   citation        the engine / quantification papers, Crossref-verified
#   caveat          optional: a limitation the report and AUDIT.md must carry
# DIA-NN's own report (native, or FragPipe's DIA route) declares nothing: its PG.MaxLFQ IS
# DIA-NN's MaxLFQ.
QUANTITY_META = "delimp.quantity."
# Crossref-verified 2026-09-29 (doi in brackets)
CITE = {
    "sage": "Sage (Lazear 2023, J Proteome Res 22:3652-3659) [10.1021/acs.jproteome.3c00486]",
    "msfragger": "MSFragger (Kong et al. 2017, Nat Methods 14:513-520) [10.1038/nmeth.4256]",
    "ionquant": "IonQuant (Yu et al. 2021, Mol Cell Proteomics 20:100077) "
                "[10.1016/j.mcpro.2021.100077]",
    "maxlfq": "MaxLFQ (Cox et al. 2014, Mol Cell Proteomics 13:2513-2526) "
              "[10.1074/mcp.M113.031591]",
    "alphadia": "AlphaDIA (Wallmann et al. 2025, Nat Biotechnol) [10.1038/s41587-025-02791-w]",
    "directlfq": "directLFQ (Ammar et al. 2023, Mol Cell Proteomics 22:100581) "
                 "[10.1016/j.mcpro.2023.100581]",
    "radiant": "Radiant DIA (Just et al. 2026, bioRxiv) [10.64898/2026.04.29.721743]",
}
PLACEHOLDER_Q = ("The report's q-value columns are 0.0 placeholders, so the DE's q-value filter "
                 "keeps every row")
# Set by main() before a search runs: the engine version it recorded (engine_version_record),
# for an adapter whose engine output does not carry its own version.
ENGINE_VERSION = None


def declare_quantity(tbl, declared):
    """`tbl` with `declared` as its QUANTITY_META schema metadata (None values left out)."""
    return tbl.replace_schema_metadata(
        {QUANTITY_META + k: str(v) for k, v in declared.items() if v is not None})


def recorded_engine_version(out):
    """The engine version for an adapter's declaration: main()'s, else the one
    search_provenance.json recorded when the job was generated (--adapt-only), else None."""
    if ENGINE_VERSION:
        return ENGINE_VERSION
    try:
        with open(os.path.join(out, "search_provenance.json")) as fh:
            return json.load(fh).get("version")
    except (OSError, ValueError, AttributeError):
        return None


def _vtag(name, ver):
    return f"{name} {ver}" if ver else f"{name} (version not recorded)"


def fragpipe_quantity_declaration(out, maxlfq_columns):
    """combined_protein.tsv's quantity, read from the fragpipe.workflow FragPipe saves beside it:
    its '# FragPipe version' / '# IonQuant version' header, the Philosopher filter it ran, and
    whether IonQuant matched between runs (ionquant.mbr)."""
    import fragpipe_params                  # the one reader of a FragPipe workflow (make_methods too)
    wf = _find(out, [fragpipe_params.SAVED_WORKFLOW])
    try:
        props, vers = fragpipe_params.read_workflow(wf) if wf else ({}, {})
    except OSError:
        props, vers = {}, {}                # said below: "not recorded"
    filt = props.get(fragpipe_params.PHILOSOPHER_FILTER_KEY)
    # the filter AS RUN, from FragPipe's log, with over which runs its rates hold -- the one
    # statement the Methods' search paragraph makes too (fragpipe_params.fdr_statement; C2
    # proteomics M1: this quoted the workflow's line, which FragPipe changes when it runs it)
    st = fragpipe_params.fdr_statement(out, filt)
    # IonQuant's own cross-run normalisation (ionquant.normalization, 1 in FragPipe's shipped
    # workflows), for the Methods' words before any between-run step of the DE (C2 stat M2)
    norm = props.get("ionquant.normalization")
    eng_norm = ("IonQuant's cross-run normalisation (ionquant.normalization=1 in "
                "fragpipe.workflow)" if norm == "1" else
                "no cross-run normalisation by IonQuant (ionquant.normalization=0 in "
                "fragpipe.workflow)" if norm == "0" else
                "whether IonQuant normalised across runs is NOT RECORDED ("
                + ("no ionquant.normalization in fragpipe.workflow" if wf else
                   "no fragpipe.workflow") + ") [confirm]")
    mbr = props.get("ionquant.mbr")
    what = "MaxLFQ intensity" if maxlfq_columns else "protein intensity (not MaxLFQ)"
    cols = "'<sample> MaxLFQ Intensity'" if maxlfq_columns else "'<sample> Intensity'"
    mbr_txt = ("match-between-runs on" if mbr == "1" else "match-between-runs off" if mbr == "0"
               else "match-between-runs not recorded" + ("" if wf else " (no fragpipe.workflow)"))
    return {
        "level": "protein",
        "engine_normalisation": eng_norm,
        "label": "FragPipe IonQuant" + (" MaxLFQ" if maxlfq_columns else ""),
        "value": f"{_vtag('FragPipe', vers.get('FragPipe'))} / "
                 f"{_vtag('IonQuant', vers.get('IonQuant'))} {what}, {mbr_txt} "
                 f"(combined_protein.tsv {cols})",
        "identification": (f"FragPipe's own FDR: Philosopher filter {st['filter']} ("
                           + ("FragPipe's log, the filter as run" if st["from_log"] else
                              "fragpipe.workflow") + "). " + st["sentence"] + " "
                           if st else "FragPipe's own FDR (its Philosopher filter is not "
                           "recorded: no fragpipe.workflow in the output). ") + PLACEHOLDER_Q,
        "q_columns": "placeholder",
        "citation": "; ".join([CITE["msfragger"], CITE["ionquant"]]
                              + ([CITE["maxlfq"]] if maxlfq_columns else [])),
    }


def _version_tuple(v):
    m = re.match(r"v?(\d+)\.(\d+)", str(v or ""))
    return (int(m.group(1)), int(m.group(2))) if m else None


def _frozen_config(out):
    """AlphaDIA's <output>/frozen_config.yaml (search_step.py writes it for every search) ->
    {"found", "version", "normalization_method", "fdr"}; a key the file lacks is None.

    `version` is the package version only from AlphaDIA 2.0.3 on (search_step.py sets it); up
    to 2.0.0 the file carries the config SCHEMA version `1` (constants/default.yaml), so only a
    number shaped like a release (\\d+.\\d+) is taken as AlphaDIA's version. `schema_v1` records
    that the file was the old shape."""
    rec = {"found": False, "version": None, "normalization_method": None, "fdr": None,
           "normalize_directlfq": None, "schema_v1": False}
    p = _find(out, ["frozen_config.yaml"])
    if not p:
        return rec
    try:
        with open(p) as fh:
            text = fh.read()
    except OSError:
        return rec
    rec["found"] = True
    try:
        import yaml
        cfg = yaml.safe_load(text) or {}
        raw = (cfg.get("version"), (cfg.get("search_output") or {}).get("normalization_method"),
               (cfg.get("fdr") or {}).get("fdr"),
               (cfg.get("search_output") or {}).get("normalize_directlfq"))
    except Exception:
        # no PyYAML: the keys are plain scalars (alphadia/constants/default.yaml)
        def key(k, indented):
            m = re.search(rf"^{'[ ]+' if indented else ''}{k}:\s*['\"]?([^'\"#\n]+)", text, re.M)
            return m.group(1).strip() if m else None
        raw = (key("version", False), key("normalization_method", True), key("fdr", True),
               key("normalize_directlfq", True))
    ver, rec["normalization_method"], rec["fdr"], nd = raw
    # AlphaDIA 2.1's search_output.normalize_directlfq (a real 2.1.2 run records `true`)
    rec["normalize_directlfq"] = (nd if isinstance(nd, bool) else
                                  {"true": True, "false": False}.get(str(nd).strip().lower())
                                  if nd is not None else None)
    if _version_tuple(ver) and "." in str(ver):
        rec["version"] = str(ver)
    elif str(ver).strip() == "1":
        rec["schema_v1"] = True
    return rec


def alphadia_quantity_declaration(out, source):
    fc = _frozen_config(out)
    ver = fc["version"] or recorded_engine_version(out)
    m = (fc["normalization_method"] or "").lower()
    vt = _version_tuple(ver)
    if m in ("directlfq", "quantselect"):
        algo = "directLFQ" if m == "directlfq" else "QuantSelect"
    elif fc["found"] and (fc["schema_v1"] or (vt and vt < (2, 1))):
        # AlphaDIA before 2.1 has no normalization_method: directLFQ is its only protein LFQ
        algo, m = "directLFQ", "directlfq"
    elif fc["found"]:
        algo = "its protein LFQ (method not recorded in frozen_config.yaml)"
    else:
        algo = "its protein LFQ (method not recorded: no frozen_config.yaml)"
    fdr = fc["fdr"]
    no_fdr = ("no fdr.fdr in frozen_config.yaml" if fc["found"] else "no frozen_config.yaml")
    # the engine's own cross-run normalisation, as its config records it, for the Methods' words
    # before any between-run step of the DE (C2 statistician M2) -- never assumed
    nd = fc["normalize_directlfq"]
    eng_norm = ("directLFQ's between-run normalisation (search_output.normalize_directlfq: true "
                "in frozen_config.yaml)" if m == "directlfq" and nd is True else
                "no directLFQ normalisation (search_output.normalize_directlfq: false in "
                "frozen_config.yaml)" if m == "directlfq" and nd is False else
                "whether AlphaDIA normalised its protein quantities across runs is NOT RECORDED "
                + ("(no search_output.normalize_directlfq in frozen_config.yaml)" if fc["found"]
                   else "(no frozen_config.yaml)") + " [confirm]")
    return {
        "level": "protein",
        "engine_normalisation": eng_norm,
        "label": "AlphaDIA" + (f" {algo}" if m in ("directlfq", "quantselect") else ""),
        "value": f"{_vtag('AlphaDIA', ver)} protein-group LFQ intensity ({algo}; {source})",
        "identification": (f"AlphaDIA's own FDR: precursors and protein groups at {fdr} "
                           f"(frozen_config.yaml fdr.fdr)" if fdr is not None else
                           f"AlphaDIA's own FDR (threshold not recorded: {no_fdr})")
                          + f". {PLACEHOLDER_Q}",
        "q_columns": "placeholder",
        "citation": "; ".join([CITE["alphadia"]] + ([CITE["directlfq"]] if m == "directlfq"
                                                     else [])),
    }


def radiant_quantity_declaration(out, intensity_col, q_col):
    protein = intensity_col.lower().startswith("pg.")
    ver = recorded_engine_version(out)
    ident = (f"Fulcrum's {q_col}, one q-value per protein and run, copied into each of the "
             f"report's q-value columns and filtered at the DE's q-value cutoff" if q_col else
             "none applied here: the Fulcrum output has no q-value column, so the report's "
             "q-value columns are 0.0 placeholders and the DE's q-value filter keeps every row")
    if protein:
        return {"level": "protein", "label": "Radiant Fulcrum",
                "value": f"{_vtag('Radiant', ver)} Fulcrum {intensity_col} per protein group",
                "identification": ident, "q_columns": "real" if q_col else "placeholder",
                "citation": CITE["radiant"]}
    # Fulcrum reported only a precursor-level quantity and the adapter kept the highest one per
    # protein and run: that is the peptide-level rollup, and it is declared as one
    return {"level": "peptide", "label": "Radiant Fulcrum",
            "value": f"{_vtag('Radiant', ver)} Fulcrum {intensity_col} precursor intensity",
            "identification": ident, "q_columns": "real" if q_col else "placeholder",
            "citation": CITE["radiant"]}


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
    """AlphaDIA pg.matrix (protein-group × run) -> DIA-NN-shaped report.parquet.
    Falls back to precursors.parquet (raw.name, pg.name, pg.intensity). The matrix is parquet
    or TSV, as the run's search_output.file_format chose: HIVE has both (AlphaDIA 2.0.1 service
    runs wrote pg.matrix.parquet; 2.1.2 runs with file_format tsv wrote pg.matrix.tsv, id column
    pg.name). Like the Sage adapter, this is the part to confirm on real data the first time."""
    try:
        import pyarrow.parquet as pq, pyarrow as pa
    except ImportError:
        sys.exit("pyarrow required to adapt AlphaDIA output. pip install pyarrow.")

    runs, prots, ints = [], [], []
    pgm = _find(out, ["pg.matrix.parquet", "pg.matrix.tsv"])
    if pgm:
        if pgm.endswith(".tsv"):
            import pyarrow.csv as pacsv
            t = pacsv.read_csv(pgm, parse_options=pacsv.ParseOptions(delimiter="\t"))
        else:
            t = pq.read_table(pgm)
        cols = t.column_names
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
            sys.exit(f"No pg.matrix.parquet / pg.matrix.tsv or precursors.parquet under {out}.")
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
    # The declaration says what the numbers are (rule 1): AlphaDIA's, not DIA-NN's.
    pq.write_table(declare_quantity(pa.table({
        "Run": runs, "Protein.Group": prots, "PG.MaxLFQ": ints,
        **{c: [0.0] * n for c in FDR_REQUIRED},
    }), alphadia_quantity_declaration(
        out, os.path.basename(pgm) if pgm else "precursors.parquet pg.intensity")), report)
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


def radiant_cys(params):
    """The --cys-alkylation choice for Radiant's DIA-NN library: Radiant takes every modification
    from the library it searches (its config has none), so cysteine is fixed there, as
    make_presets.py recorded it in <params>.provenance.json. A config without that record gets the
    default, cam, said on stderr. Before 2.12 no Radiant library fixed anything on C: on HIVE
    (DIA-NN 2.7.0, job 24807967) 0 of 6,044 C-peptides carried UniMod:4."""
    try:
        with open(f"{params}.provenance.json") as fh:
            value = (json.load(fh).get("cys_alkylation") or {}).get("value")
    except (OSError, ValueError, TypeError, AttributeError):
        value = None
    if value is None:
        sys.stderr.write(f"[run_search] NOTE: {params} has no make_presets.py record of the "
                         "cysteine alkylation, so Radiant's library fixes carbamidomethyl on C "
                         "(the default, not user-confirmed); make the config with make_presets.py "
                         "--cys-alkylation to choose\n")
    return value or "cam"


def radiant_library(tools, fasta, out, threads, params=None):
    """Build the spectral library Radiant searches against, from DIA-NN's predictor.

    NOT a straight hand-off: DIA-NN 2.x writes .predicted.speclib v-10/-11 and Radiant's
    reader only supports v>=-3, so its binary is rejected outright. The conversion lives
    in make_radiant_library.py, which also decides whether a library already there can be reused
    (only one built with the same cysteine modification).
    """
    libdir = radiant_lib_dir(out)
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
    res = subprocess.run([sys.executable, *radiant_library_argv(dn, fasta, libdir, threads,
                                                                radiant_cys(params))],
                         capture_output=True, text=True)
    sys.stderr.write(res.stderr)
    if res.returncode != 0:
        sys.exit(f"make_radiant_library.py failed (exit {res.returncode}).")
    info = json.loads(res.stdout[res.stdout.index("{"):])
    print(f"  [radiant] library -> {info['library']} ({info.get('rows', '?')} rows)")
    return info["library"]


def radiant_library_job(tools, fasta, out, threads, params=None):
    """radiant_library() for a --sbatch search: (the library path, the job lines that build it).
    The same helper writes the same file (make_radiant_library.py keeps its output name, and
    reuses one already there), but on the compute node: predicting a library is a DIA-NN run,
    and it used to run here, on the login node, before the job was written."""
    dn = tools.get("diann")
    if not dn:
        sys.exit("Radiant needs a spectral library, but tools.json has no DIA-NN command. "
                 "Acquire DIA-NN first (it is the library generator for this route), or "
                 "pass --library with an existing .tsv/.csv/.fragLibFF library.")
    return radiant_library_lines(dn, fasta, radiant_lib_dir(out), threads, radiant_cys(params))


def radiant_lib_dir(out):
    """Where a search's Radiant library is built: <out>/radiant_lib, on every route."""
    return os.path.join(out, "radiant_lib")


def radiant_library_argv(diann, fasta, libdir, threads, cys):
    """make_radiant_library.py and its arguments: THE argv that builds a Radiant library (DIA-NN's
    predictor -> parquet -> the TSV Radiant reads) into `libdir`, for the inline route
    (radiant_library) and the job lines (radiant_library_lines)."""
    helper = os.path.join(os.path.dirname(os.path.abspath(__file__)), "make_radiant_library.py")
    return [helper, "--diann", diann, "--fasta", fasta, "--out-dir", libdir,
            "--threads", str(threads), "--cys-alkylation", cys]


def radiant_library_lines(diann, fasta, libdir, threads, cys):
    """The job lines that build a Radiant library on a compute node: radiant_library_argv() run
    with the job's python, then a guard that it wrote the TSV. -> (the TSV, the lines). The
    single-job route (radiant_library_job) and the parallel route's step 1 (radiant_parallel.py)
    both use it: that step handed Radiant DIA-NN 2.x's .predicted.speclib, a format Radiant
    rejects (S8)."""
    from make_radiant_library import LIBRARY_TSV
    q = shlex.quote
    tsv = os.path.join(libdir, LIBRARY_TSV)
    return tsv, [
        'echo "building the Radiant spectral library (DIA-NN predictor) on this compute node"',
        _job_python(),
        '"$_py" ' + " ".join(q(w) for w in radiant_library_argv(diann, fasta, libdir, threads,
                                                                  cys)),
        f"[ -s {q(tsv)} ] || {{ echo {q(f'FAILED: make_radiant_library.py wrote no {tsv}')} "
        ">&2; exit 1; }"]


def run_radiant_parallel(tools, params, files, fasta, out, threads, a, library=None):
    """On a cluster, search each file as its own array task, then rescore once.

    Radiant's Fulcrum backend is SERIAL (`NotImplementedError` for parallel mode), so
    an N-file study otherwise costs N x one-file wall-clock even on a big node. The
    search is per-file independent though; only the downstream rescoring/FDR/rollup
    needs the whole set. radiant_parallel.py splits exactly there and emits a 3-step
    chain. Emits scripts; the orchestrator submits them with dependencies."""
    os.makedirs(out, exist_ok=True)
    listing = os.path.join(out, "radiant_input_files.txt")
    write_file(listing, "\n".join(files) + "\n")
    argv = [sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "radiant_parallel.py"),
            "--runtime", tools.get("radiant_runtime") or "docker",
            "--image", tools.get("radiant_image") or "",
            "--raw-list", listing, "--fasta", fasta, "--config", params,
            "--out", out, "--threads-per-file", str(threads)]
    if library:
        argv += ["--library", library]
    elif tools.get("diann"):
        argv += ["--diann", tools["diann"], "--cys-alkylation", radiant_cys(params)]
    for flag, val in (("--partition", a.partition), ("--account", a.account),
                      ("--qos", a.qos), ("--max-simultaneous", a.max_simultaneous),
                      ("--exclude", _exclude_arg(a))):
        if val:
            argv += [flag, str(val)]
    if getattr(a, "mbr", False):
        argv.append("--mbr")
    argv += _job_end_argv(a)
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
    # .raw -> mzML, planned as for Sage (mzml_plan): converted in the job with --sbatch
    mzml, steps = mzml_plan(files, out)
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
    if not sbatch:
        convert_inline(steps)               # an inline search converts here, now
        steps = []
    # The predicted library is a DIA-NN run: in the job too with --sbatch, not here.
    lib_lines = []
    if library:
        lib = library
    elif sbatch:
        lib, lib_lines = radiant_library_job(tools, fasta, out, threads, params)
    else:
        lib = radiant_library(tools, fasta, out, threads, params)

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
        job = "\n".join(conversion_lines(steps) + lib_lines + [full])
        emit_sbatch(sbatch, job, out, threads, job="radiant_search", **(queue or {}))
        return {"engine": "radiant", "out": out, "submitted": sbatch, "ran": False,
                "library": lib, "library_built": "in the job" if lib_lines else "given",
                "mzml_conversion": conversion_record(steps, "in the job, before Radiant")}
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
    # as DIA-NN writes Protein.Group (bare accessions; group_accessions leaves a bare one as it
    # is), so a FASTA-header-shaped id cannot hide a Cont_ entry from the contaminant filter
    pgs = [group_accessions(p) for p in t.column(c_pg).to_pylist()]
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
    # The declaration says what the numbers are (rule 1): Fulcrum's, not DIA-NN's.
    pq.write_table(declare_quantity(pa.table({
        "Run": runs2, "Protein.Group": pgs2, "PG.MaxLFQ": ints,
        **{c: qv for c in FDR_REQUIRED},
    }), radiant_quantity_declaration(out, c_int, c_q)), report)
    print(f"  [adapt] Radiant/Fulcrum -> {report}  ({n} protein×run rows)")
    return report


# ------------------------------------------------------------ mzML inputs -----
# Sage and Radiant read mzML. The conversion is PLANNED here (mzml_plan -- nothing is run) and
# RUN wherever the search runs: in this process for an inline search (convert_inline), inside
# the job for --sbatch (conversion_lines), on the compute node. It used to run right here,
# before the job script was even written: gabrig 2026-09-29, 4 Fusion Lumos .raw (~2.4 GB) on
# HIVE login2 -- `run_search.py --engine sage --sbatch` started msconvert on the LOGIN NODE for
# every file (golden rule 3). It failed at once only because bioconda's Linux msconvert is built
# from ProteoWizard's source WITHOUT the vendor readers (recipe: pwiz-src-without-v), so it
# cannot open a .raw at all.
#
# .raw -> ThermoRawFileParser (-f=2, indexed mzML): the parser step 2 already reads .raw with,
# found the same way (detect_acquisition.locate_trfp: $THERMORAWFILEPARSER, PATH, the pipeline
# env, the Core's shared copy) and given the same checked DOTNET_ROOT (trfp_launch -- directory
# checks only, nothing is started here). msconvert only when there is no usable parser.
# .d -> searched as it is for a Sage search for identifications only (dotd="direct": Sage 0.14.7
# reads a timsTOF .d itself, but has no MS1 quantification on one); otherwise msconvert in the
# Core's pinned vendor-enabled ProteoWizard container (pwiz_step), never msconvert from PATH:
# bioconda's has no Bruker reader (I45).
# Anything else not mzML -> msconvert.
MZML_SUFFIXES = (".mzml", ".mzml.gz")
TRFP_UNSET_ENV = ("DOTNET_ROOT_X64", "DOTNET_ROOT_ARM64", "DOTNET_ROOT_X86")  # see _child_env


def _trfp_converter():
    """(argv prefix, {env}, label) for ThermoRawFileParser, or (None, None, why not)."""
    import detect_acquisition as da
    cmd, where = da.locate_trfp()
    if not cmd:
        return None, None, da._trfp_not_found()
    launch = da.trfp_launch(cmd)
    if launch["problem"]:
        return None, None, launch["problem"]
    env = {"DOTNET_ROOT": launch["dotnet_root"]} if launch["dotnet_root"] else {}
    return launch["cmd"], env, f"ThermoRawFileParser ({where})"


# The container that converts a timsTOF .d for a quantifying Sage search: by default the Core's
# pinned vendor-enabled ProteoWizard image, the one skyline_peaks.py runs SkylineCmd from
# (resolve_defaults.SKYLINE_PIN, read here, never changed). msconvert is in it under Wine.
PWIZ_SIF_FLAG = "--pwiz-sif"

# msconvert's options for a timsTOF .d: the Core's Windows msconvert settings (ProteoWizard
# 3.0.26058), as recorded in the <dataProcessing> of the mzML it made for a Core submission (arm C), and
# written as MSConvertGUI writes them (sumMs1=0 is "false"). Without the combine and the summing,
# msconvert writes one spectrum per mobility scan: on HIVE one run gave 11,139,896 spectra ("frame=1
# scan=1") and passed 9.9 GB, against 82,905 combined ones ("merged=0 frame=1 scanStart=1
# scanEnd=811") -- H2-F1, 2026-11-01. --combineIonMobilitySpectra is an option, not a filter;
# vendor peakPicking must be the FIRST filter (MSConvertGUI: "it only works on raw, untransformed
# data"). The step's argv and every record of it read this one tuple.
PWIZ_TIMS_MSCONVERT = ("--mzML", "--zlib", "--combineIonMobilitySpectra",
                       "--filter", "peakPicking vendor msLevel=1-",
                       "--filter", "scanSumming precursorTol=0.05 scanTimeTol=10 "
                                   "ionMobilityTol=0.1 sumMs1=0")
# sumMs1=0 is kept because it is the settings' source and the HIVE acceptance ran with it -- but
# the pinned 3.0.26266's scanSumming usage (msconvert --help, lines 410-416) has no sumMs1, and
# its mzML records sumMS1 value="true" (Windows 3.0.26058: "false"). MS1 is not summed in effect:
# spectrum counts, IDs, MS1 peaks and LFQ all equal the Windows output (H2-7, 2026-11-01). The
# records say so beside the options (PWIZ_SUMMS1_NOTE); the Methods sentence does not need it.
PWIZ_SUMMS1_NOTE = ("the pinned ProteoWizard 3.0.26266 does not take sumMs1 and records "
                    "sumMS1=true in its mzML; MS1 spectrum counts equal the Windows 3.0.26058 "
                    "output, so MS1 is not summed in effect (HIVE check, 2026-11-01)")
# (in records, which are delivered: no submission number -- L7R-F10)
PWIZ_TIMS_MSCONVERT_SOURCE = ("the Core's Windows msconvert settings (ProteoWizard 3.0.26058) for "
                              "timsTOF HT DDA-PASEF, as recorded in the <dataProcessing> of an "
                              "mzML it made (2026-10)")


def pwiz_msconvert_options(windows=False):
    """PWIZ_TIMS_MSCONVERT as one string: shell-quoted for the records, or with double quotes for
    a command a person types into a Windows prompt (the refusal's way out)."""
    if windows:
        return " ".join(f'"{x}"' if " " in x else x for x in PWIZ_TIMS_MSCONVERT)
    return " ".join(shlex.quote(x) for x in PWIZ_TIMS_MSCONVERT)


def pwiz_rerun_args(sif):
    """What a re-run of this search must add to get the same container: nothing for the pin's
    own path (the default), else --pwiz-sif <sif> (search_provenance.json `pwiz_sif_arg`, read
    by sage_lfq_check.rerun_plan and provenance.py's reproduce.sh)."""
    if os.path.abspath(sif) == os.path.abspath(_skyline_pin()["sif"]):
        return []
    return [PWIZ_SIF_FLAG, sif]


def _skyline_pin():
    from resolve_defaults import SKYLINE_PIN
    return SKYLINE_PIN


def _pwiz_pinned(sif):
    """The Core's pinned image (by its file name, as skyline_peaks.verify_container decides)."""
    return os.path.basename(sif) == os.path.basename(_skyline_pin()["sif"])


def _pwiz_sha_failed(sif):
    """The one message for a pinned-name container whose bytes are not the pin's."""
    return f"FAILED: the ProteoWizard container {sif} does not match its pinned sha256"


def _pwiz_sif_or_refuse(sif):
    """The container's absolute path; refused, before anything is written, when it is not there.
    Off HIVE it usually is not: the refusal names the two ways out."""
    sif = os.path.abspath(sif or _skyline_pin()["sif"])
    if not os.path.isfile(sif):
        which = (f"the Core's pinned image on HIVE; another can be given with {PWIZ_SIF_FLAG}"
                 if sif == _skyline_pin()["sif"] else
                 f"{PWIZ_SIF_FLAG}; the Core's on HIVE is {_skyline_pin()['sif']}")
        sys.exit(
            f"[run_search] REFUSED -- Nothing was converted, written or submitted. Sage 0.14.7 has "
            f"no MS1 quantification on a timsTOF .d, so a quantifying search converts each .d to "
            f"mzML in the job with a vendor-enabled ProteoWizard container, and there is none at "
            f"{sif} ({which}). Two ways out:\n  1. convert each .d to mzML first with a "
            f"vendor-enabled msconvert (e.g. ProteoWizard on Windows: msconvert <run>.d "
            f"{pwiz_msconvert_options(windows=True)}) and search the mzML; or\n"
            f"  2. search for identifications "
            f"only: estimate_params.py --engine sage --ids-only writes a cfg (quant.lfq false) "
            f"with which Sage reads the .d itself -- no quantities, so no DE.")
    return sif


def pwiz_step(d_path, partial, final, sif):
    """The conversion step for one timsTOF .d in the ProteoWizard container (see mzml_plan for
    the step's keys; `sif` added). apptainer as skyline_peaks.apptainer_argv runs it: the .d's
    folder read-only (as given and as its real path, for a .d reached through a link), the
    partial mzML folder read-write, parents bound first. msconvert is given plain Linux paths, as
    ProteoWizard's own container README runs it (`wine msconvert /data/file.raw`)."""
    d = os.path.abspath(d_path.rstrip("/"))
    partial_dir = os.path.dirname(partial)
    binds = [(p, ":ro") for p in dict.fromkeys([os.path.dirname(d),
                                                 os.path.dirname(os.path.realpath(d))])]
    binds.append((partial_dir, ""))
    for b, _mode in binds:
        if ":" in b or "," in b:
            sys.exit(f"[run_search] REFUSED -- Nothing was converted, written or submitted. The "
                     f"folder {b} holds ':' or ',', which an apptainer --bind cannot carry, so "
                     f"{d} cannot be converted in the ProteoWizard container. Search it from a "
                     f"folder whose path has neither (a link to it is enough), or with "
                     f"estimate_params.py --ids-only, where Sage reads the .d itself.")
    binds.sort(key=lambda b: (b[0].rstrip("/").count("/"), b[0], b[1]))
    # 2.15: build this with skyline_peaks.apptainer_argv, the one definition (L7R-F12); kept as is
    # in 2.14.2 because the HIVE acceptance run used this argv.
    argv = ["apptainer", "exec", "--cleanenv", "--env", "WINEDEBUG=-all"]
    for b, mode in binds:
        argv += ["--bind", f"{b}:{b}{mode}"]
    tool = (f"ProteoWizard msconvert (pinned Skyline container {_skyline_pin()['version']})"
            if _pwiz_pinned(sif) else
            f"ProteoWizard msconvert (container {os.path.basename(sif)}, not the Core's pinned "
            f"image)")
    # centroided MS1 and MS2 (Sage wants peaks, not profile), mobility scans combined
    return {"tool": tool, "sif": sif,
            "argv": argv + [sif, "wine", "msconvert", d, *PWIZ_TIMS_MSCONVERT,
                            "-o", partial_dir],
            "env": {}, "complete": "</indexedmzML>", "input": d_path, "mzml": final,
            "partial": partial}


def pwiz_sif_lines(sif):
    """The job's check of the container, before any conversion: the Core's pinned image by its
    sha256; any other sif is used as given (skyline_peaks.verify_container's rule)."""
    if not _pwiz_pinned(sif):
        return []
    q = shlex.quote
    return [f"echo {q(_skyline_pin()['sif_sha256'] + '  ' + sif)} | sha256sum -c --quiet - "
            f"|| {{ echo {q(_pwiz_sha_failed(sif))} >&2; exit 1; }}"]


def pwiz_container_record(sif):
    """What search_provenance.json says of the container. Only the pinned image's version is
    known; another sif's is not checked, and none is written for it."""
    pin = _skyline_pin()
    how = {"msconvert_options": pwiz_msconvert_options(),
           "msconvert_options_source": PWIZ_TIMS_MSCONVERT_SOURCE}
    if _pwiz_pinned(sif):
        how["msconvert_options_note"] = PWIZ_SUMMS1_NOTE
    if _pwiz_pinned(sif):
        ver = pin.get("msconvert_version")
        return {"sif": sif, "image": pin["image"], "tag": pin["tag"], "version": pin["version"],
                "sif_sha256": pin["sif_sha256"], "docker_digest": pin["docker_digest"],
                "msconvert_version": ver,
                "msconvert_version_source": (
                    "resolve_defaults.SKYLINE_PIN msconvert_version, read from the pinned image "
                    "(`wine msconvert --help`, its ProteoWizard release line)" if ver
                    else "not recorded [DEFAULT -- the pinned image's own msconvert version has "
                         "not been read from it yet]"),
                # written at generation: what will be checked, not that it was (L7R-F14)
                "verified": ("sha256 checked by the job before the first conversion (it stops on "
                             "a mismatch); an inline search checks it before converting"), **how}
    return {"sif": sif, "version": None, "sif_sha256": None, "msconvert_version": None,
            "msconvert_version_source": "not recorded: not the Core's pinned image",
            "verified": "not checked: not the Core's pinned image, so its ProteoWizard version "
                        "and sha256 are not known", **how}


def sage_run_name(path):
    """The run an input is named by in Sage's output, as the adapters write it: its file name
    without the folder and the last extension (S1.d, S1.raw -> S1.mzML, S1.mzML -> S1). The ONE
    rule for mzml_plan's duplicate check, adapt_sage and adapt_sage_ids (L7R-F8)."""
    return os.path.splitext(os.path.basename(str(path).rstrip("/\\")))[0]


def mzml_plan(files, out, dotd="convert", pwiz_sif=None):
    """(what each input is searched as, in order; the conversion steps). Nothing is run or
    created. A step is {"input", "mzml", "partial", "tool", "argv", "env", "complete"}: `argv`
    (with `env` added) writes `partial`, which is moved to `mzml` only once `complete` -- the
    closing tag an indexed mzML ends with, or None for "non-empty" -- is found in its tail.

    `dotd`: "convert" (every .d becomes an mzML, in the ProteoWizard container `pwiz_sif`,
    default the Core's pinned one -- refused when it is not there) or "direct" (a .d is searched
    as given, with no step: run_sage() for a search with quant.lfq false). .raw and mzML are the
    same either way."""
    if dotd not in ("convert", "direct"):
        raise ValueError(f"mzml_plan: dotd must be 'convert' or 'direct', not {dotd!r}")
    mzdir = os.path.join(out, "mzml")
    partial_dir = os.path.join(mzdir, ".partial")
    mzml, steps, trfp, msconvert, sif = [], [], None, None, None
    for f in files:
        if f.lower().endswith(MZML_SUFFIXES):
            mzml.append(f)
            continue
        if dotd == "direct" and _is_dotd(f):
            mzml.append(f)
            continue
        base = sage_run_name(f)
        final, partial = os.path.join(mzdir, base + ".mzML"), os.path.join(partial_dir,
                                                                            base + ".mzML")
        if _is_dotd(f):
            if sif is None:
                sif = _pwiz_sif_or_refuse(pwiz_sif)
            steps.append(pwiz_step(f, partial, final, sif))
            mzml.append(final)
            continue
        tool = None
        if f.lower().endswith(".raw"):
            if trfp is None:
                trfp = _trfp_converter()
            argv0, env, label = trfp
            if argv0:
                # every value inlined: TRFP takes optional parameters only as -option=value
                tool = {"tool": label, "argv": argv0 + [f"-i={f}", f"-b={partial}", "-f=2"],
                        "env": env, "complete": "</indexedmzML>"}
        if tool is None:
            if msconvert is None:
                msconvert = shutil.which("msconvert") or ""
            if not msconvert:
                why = (f"ThermoRawFileParser: {trfp[2]}\n  " if trfp and not trfp[0] else "")
                sys.exit("[run_search] Sage/Radiant read mzML, and this input is not mzML: "
                         f"{f}\n  Nothing was converted, written or submitted.\n  {why}"
                         "msconvert: not on PATH (other formats need it).\n  Convert to mzML "
                         "first, or install ThermoRawFileParser for .raw (bash scripts/setup.sh "
                         "puts bioconda's in the pipeline env).")
            if sys.platform.startswith("linux"):
                why = (f" because ThermoRawFileParser is not usable here ({trfp[2]})"
                       if f.lower().endswith(".raw") and trfp else "")
                sys.stderr.write(
                    f"[run_search] WARNING: converting {os.path.basename(f.rstrip('/'))} with "
                    f"msconvert ({msconvert}){why}. bioconda's Linux msconvert has no vendor "
                    "readers and cannot open .raw or .d; this works only with a vendor-enabled "
                    "build. (A .raw normally goes through ThermoRawFileParser. A timsTOF .d never "
                    "comes here: Sage reads it itself for identifications only -- "
                    "estimate_params.py --ids-only -- and to quantify, the job converts it in the "
                    "Core's pinned ProteoWizard container.)\n")
            # centroided MS1 and MS2 (Sage and Radiant want peaks, not profile), and the same
            # completeness check as the parser: an indexed mzML ends with its closing tag
            tool = {"tool": f"msconvert ({msconvert})",
                    "argv": [msconvert, f, "--mzML", "--zlib",
                             "--filter", "peakPicking vendor msLevel=1-", "-o", partial_dir],
                    "env": {}, "complete": "</indexedmzML>"}
        steps.append(dict(tool, input=f, mzml=final, partial=partial))
        mzml.append(final)
    # One input per run name: two inputs called s1 (from two folders) would overwrite each
    # other's conversion, and a run1.d searched as it is beside a run1.raw converted to run1.mzML
    # would be one run, run1, in the adapters' output (sage_run_name).
    seen = {}
    for m in mzml:
        seen.setdefault(sage_run_name(m), []).append(m)
    dupes = {k: v for k, v in seen.items() if len(v) > 1}
    if dupes:
        sys.exit("[run_search] inputs share a Sage run name, so they would overwrite each other "
                 "or merge into one run: " + "; ".join(f"{k}: {v}" for k, v in dupes.items())
                 + ". Rename them, or search them separately. Nothing was converted, written or "
                   "submitted.")
    return mzml, steps


def _step_argv(step):
    """The step as one argv, its environment applied with env(1) -- in the job and inline alike."""
    if not step["env"]:
        return list(step["argv"])
    unset = [x for k in TRFP_UNSET_ENV for x in ("-u", k)]
    return (["env"] + unset + [f"{k}={v}" for k, v in step["env"].items()]
            + list(step["argv"]))


def conversion_lines(steps):
    """Bash lines that do the conversions on the compute node, first thing in the job. A
    converter that exits 0 without a complete file stops the job here, before the search."""
    if not steps:
        return []
    q = shlex.quote
    tools = sorted({s["tool"] for s in steps})
    # the ProteoWizard container first: a pinned image that is not the pin stops the job here
    lines = [ln for sif in dict.fromkeys(s["sif"] for s in steps if s.get("sif"))
             for ln in pwiz_sif_lines(sif)]
    lines += [f'echo "converting {len(steps)} input file(s) to mzML on this compute node: '
              f'{"; ".join(tools)}"',
              f"mkdir -p {q(os.path.dirname(steps[0]['partial']))}"]
    for s in steps:
        # quoted as one word: an input path is data, never shell. `|| { ... }` on the converter
        # itself: under set -e a converter that exits non-zero would end the job before any
        # FAILED line was printed (and the watcher would have nothing to classify).
        fail = (f"{{ echo {q(_no_mzml(s))} >&2; exit 1; }}")
        # `case` on the file's tail, not `tail | grep -q`: under pipefail a grep that stops
        # early can fail the pipeline with tail's SIGPIPE
        check = ([f"case \"$(tail -c 256 {q(s['partial'])} 2>/dev/null)\" in "
                  f"*{q(s['complete'])}*) ;; *) {fail} ;; esac"] if s["complete"] else
                 [f"[ -s {q(s['partial'])} ] || {fail}"])
        lines += [f"rm -f {q(s['mzml'])} {q(s['partial'])}",
                  " ".join(q(x) for x in _step_argv(s)) + f" || {fail}",
                  *check,
                  f"mv -f {q(s['partial'])} {q(s['mzml'])}"]
    return lines


def _no_mzml(step):
    """The one failure message (job and inline): 'no mzML' is what watch_run.sh classifies as
    sage_no_mzml."""
    return (f"FAILED: no mzML from {step['tool']} for {step['input']} (expected a complete "
            f"{step['partial']}). Nothing was searched.")


def convert_inline(steps):
    """Run the conversions here (an inline search: never on a login node -- main() refuses)."""
    for sif in dict.fromkeys(s["sif"] for s in steps if s.get("sif")):
        if _pwiz_pinned(sif):               # pwiz_sif_lines' check, in Python
            import skyline_peaks
            if skyline_peaks.sha256_file(sif) != _skyline_pin()["sif_sha256"]:
                sys.exit(_pwiz_sha_failed(sif))
    for s in steps:
        os.makedirs(os.path.dirname(s["partial"]), exist_ok=True)
        for p in (s["mzml"], s["partial"]):
            if os.path.lexists(p):
                os.remove(p)
        try:
            sh(" ".join(shlex.quote(x) for x in _step_argv(s)))
        except subprocess.CalledProcessError as e:
            sys.exit(f"{_no_mzml(s)} (it exited {e.returncode})")
        ok = os.path.isfile(s["partial"]) and os.path.getsize(s["partial"]) > 0
        if ok and s["complete"]:
            with open(s["partial"], "rb") as fh:
                fh.seek(max(0, os.path.getsize(s["partial"]) - 256))
                ok = s["complete"].encode() in fh.read()
        if not ok:
            sys.exit(_no_mzml(s))
        os.replace(s["partial"], s["mzml"])


def conversion_record(steps, where):
    """What converted what, for search_provenance.json."""
    if not steps:
        return None
    rec = {"where": where, "tools": sorted({s["tool"] for s in steps}),
           "files": [{"input": s["input"], "mzml": s["mzml"],
                      "command": " ".join(shlex.quote(x) for x in _step_argv(s))}
                     for s in steps]}
    sifs = list(dict.fromkeys(s["sif"] for s in steps if s.get("sif")))
    if sifs:
        rec["containers"] = [pwiz_container_record(sif) for sif in sifs]
    return rec


# ------------------------------------------------------------------- Sage -----
def refuse_sage_fragment_mismatch(bundle, params):
    """Stop before anything is generated when the Sage cfg's fragment window does not fit the MS2
    analyzer workflow.manifest.json records (estimate_params.sage_fragment_mismatch() decides).
    An ion-trap manifest with a ppm cfg -- --ms2-analyzer given in step 4 but not in 6b -- would
    otherwise search every spectrum at +/-10 ppm and match almost no ion-trap fragments."""
    import estimate_params
    analyzer = (bundle.get("resolution") or {}).get("ms2_analyzer")
    try:
        with open(params) as fh:
            cfg = json.load(fh)
    except (OSError, ValueError):
        return                      # not a readable JSON: Sage itself rejects it, with its reason
    hit = estimate_params.sage_fragment_mismatch(
        estimate_params.ms2_analyzer_arg(analyzer) if analyzer else None, cfg)
    if hit and hit[0] == "refuse":
        sys.exit(f"[run_search] REFUSED -- nothing was converted, written or submitted. "
                 f"{params}: {hit[1]}")
    if hit:
        print(f"[run_search] WARNING: {params}: {hit[1]}", file=sys.stderr)


def search_digest(engine, params):
    """This search's in-silico digest in estimate_params.DIANN_DIGEST's shape -- only the keys
    the params file states -- or None when it cannot be read. DIA-NN: the cfg's --cut /
    --missed-cleavages / --min-/--max-pep-len / --met-excision. Sage: its database.enzyme block,
    with cleave_at/restrict written in DIA-NN's --cut syntax (KR + restrict P = K*,R*,!*P)."""
    try:
        if engine == "diann":
            groups = dict(_diann_parallel_mod().cfg_groups(_diann_parallel_mod().cfg_tokens(params)))
            got = {"met_excision": "--met-excision" in groups}
            for flag, key, conv in (("--cut", "cut", str), ("--missed-cleavages", "missed_cleavages", int),
                                    ("--min-pep-len", "min_pep_len", int),
                                    ("--max-pep-len", "max_pep_len", int)):
                if groups.get(flag):
                    got[key] = conv(groups[flag][0])
            return got
        if engine == "sage":
            with open(params) as fh:
                return sage_digest((json.load(fh).get("database") or {}).get("enzyme") or {})
    except (OSError, ValueError, TypeError):
        return None
    return None


def sage_digest(enz):
    """A Sage database.enzyme block in DIANN_DIGEST's shape (keys it states only)."""
    got = {k2: enz[k1] for k1, k2 in (("missed_cleavages", "missed_cleavages"),
                                       ("min_len", "min_pep_len"),
                                       ("max_len", "max_pep_len")) if k1 in enz}
    if enz.get("cleave_at") == "":
        got["cut"] = "(none: a non-specific digest)"   # compared, never dropped (FU3, H7 bug 5)
    elif enz.get("cleave_at") == "$":
        got["cut"] = "(none: no digestion)"
    elif enz.get("cleave_at"):
        got["cut"] = ",".join(f"{c}*" for c in enz["cleave_at"]) + (
            f",!*{enz['restrict']}" if enz.get("restrict") else "")
    return got


def warn_contaminant_digest(engine, params, fasta, sidecar=None):
    """Warn when this search digests differently from the digest fetch_fasta.py judged the
    contaminants on (the sidecar's contaminant_digest, estimate_params.DIANN_DIGEST): its peptide
    rule dropped the contaminant entries with fewer than min_unique_peptides peptides of their
    own IN THAT DIGEST, so under another one -- an --overrides change, or Sage's own 2 missed
    cleavages and proline rule -- an entry near that cut may be judged differently. Sage's known
    default difference is one INFO line; anything else is a WARNING. Never a stop: the database
    is still correct for the identity rule, and the list is on record. Returns the message (or
    None). `sidecar`: fasta_sidecar_or_refuse's (default <fasta>.meta.json)."""
    sidecar = sidecar or fasta + ".meta.json"
    try:
        with open(sidecar) as fh:
            meta = json.load(fh)
    except (OSError, ValueError):
        return None
    want = meta.get("contaminant_digest")
    if not isinstance(want, dict) or not meta.get("min_unique_peptides"):
        return None                 # no peptide rule was applied, so no digest to disagree with
    got = search_digest(engine, params)
    if not got:
        return None
    differ = [k for k in sorted(got) if k in want and got[k] != want[k]]
    if not differ:
        return None
    # Sage's own default digest (estimate_params.SAGE_ENZYME) always differs from DIA-NN's: that
    # is KNOWN, and a WARNING on every Sage search would teach people to ignore the warning. It
    # gets one INFO line; a difference the default does not explain (an --overrides change of
    # cut, missed cleavages or length) stays a WARNING.
    import estimate_params
    known = sage_digest(estimate_params.SAGE_ENZYME) if engine == "sage" else {}
    expected = [k for k in differ if k in known and got[k] == known[k]]
    unexpected = [k for k in differ if k not in expected]
    fmt = lambda ks: "; ".join(f"{k} {got[k]!r} (FASTA: {want[k]!r})" for k in ks)   # noqa: E731
    rule = (f"its peptide rule dropped contaminant entries with fewer than "
            f"{meta['min_unique_peptides']} peptides of their own in THAT digest, so an entry near "
            f"the cut may be judged differently here (contaminants_dropped_as_target in {sidecar} "
            f"has each entry's peptide counts)")
    if not unexpected:
        msg = (f"[run_search] INFO: Sage's default digest ({fmt(expected)}) differs from the "
               f"DIA-NN digest fetch_fasta.py judged the contaminants on -- expected for every Sage "
               f"search; {rule}.")
    else:
        msg = (f"[run_search] WARNING: this search's digest differs from the one fetch_fasta.py "
               f"judged the contaminants on ({sidecar}, contaminant_digest): {fmt(unexpected)}"
               + (f" (besides Sage's known default: {fmt(expected)})" if expected else "")
               + f". {rule[0].upper()}{rule[1:]}.")
    print(msg, file=sys.stderr)
    return msg


def fasta_sidecar_or_refuse(a, engine):
    """The sidecar describing --fasta (fetch_fasta.sidecar_for, the one reading), looked for also
    in the session's input/ -- or REFUSED, before anything is converted, written or submitted,
    when --fasta is a Philosopher target+decoy .fas whose target FASTA's sidecar is not found
    (H7b-1: the record then said keratin_sample_source null, and nothing that reads the sidecar
    from the record's `fasta` found one). Never null facts recorded silently. -> the record
    (search_provenance.json `fasta_meta`)."""
    import fetch_fasta as ff
    look_in = []
    try:
        import session as session_mod
        sess = session_mod.session_of_search(a.out)
        if sess:
            look_in.append(os.path.join(sess, "input"))
    except ImportError:
        pass                                 # a partial scripts/: the .fas's own folder only
    tag = (_workflow_value(a.params, "database.decoy-tag") if engine == "fragpipe" and a.params
           else None) or "rev_"
    try:
        rec = ff.sidecar_for(a.fasta, getattr(a, "fasta_meta", None), look_in, tag)
    except OSError as e:
        sys.exit(f"[run_search] REFUSED -- nothing was converted, written or submitted. --fasta "
                 f"{a.fasta} could not be read to find its sidecar ({e}).")
    if rec["problem"]:
        sys.exit(f"[run_search] REFUSED -- nothing was converted, written or submitted. "
                 f"{rec['problem']}. Its facts (keratin sample, contaminant set, organism) would "
                 f"be recorded as unknown. Give the target FASTA's sidecar: --fasta-meta "
                 f"<target FASTA>.meta.json (fetch_fasta.py wrote it beside the FASTA the .fas "
                 f"was built from), or --fasta-meta {ff.NO_SIDECAR} if the database was not "
                 f"built by fetch_fasta.py.")
    if rec["path"] and not os.path.isfile(rec["path"]):
        sys.exit(f"[run_search] REFUSED -- nothing was converted, written or submitted. "
                 f"--fasta-meta {rec['path']}: no such file.")
    return rec


def keratin_sample_check(fasta, flag, engine="diann", sidecar=None):
    """A keratin sample (hair, wool, feather, skin, nail ...) is never searched against a database
    that still holds keratin-family Cont_ entries. Keratin is the analyte there, and those entries
    take its peptides away: DIA-NN's --cont-quant-exclude Cont_ keeps every peptide they share with
    the sample's keratins out of quantification, and run_de.R removes every precursor that names a
    Cont_ entry. msalemi's hair benchmark (SET28, skill 2.8.0) lost them exactly this way, with
    --keratin-sample reaching only the auditors after the search.

    The sample is a keratin sample when --keratin-sample says so OR the FASTA's sidecar was built
    with fetch --keratin-sample (fetch_fasta.keratin_sample_recorded). The FASTA's own headers are
    checked with fetch_fasta's one keratin rule, sidecar or not. -> the record search_provenance.json
    keeps; exits (REFUSED) before anything is converted, written or submitted. `sidecar`: the one
    fasta_sidecar_or_refuse found (a Philosopher .fas: its target FASTA's), else <fasta>.meta.json."""
    import fetch_fasta as ff
    sidecar = sidecar or fasta + ".meta.json"
    try:
        with open(sidecar) as fh:
            meta = json.load(fh)
    except (OSError, ValueError):
        meta = None
    recorded = ff.keratin_sample_recorded(meta)
    rec = {"value": bool(flag) or recorded is True,
           "source": "--keratin-sample" if flag else (sidecar if recorded is True else None),
           # was the question asked? "user" (this flag, or the sidecar's answer), "default"
           # (the sidecar says nobody answered), None (nothing recorded) -- rule 2
           "sample_source": ("user" if flag else
                             (meta.get("keratin_sample_source") or "default")
                             if recorded is not None else None),
           "sidecar_keratin_sample": recorded, "fasta": fasta,
           "keratin_contaminants_in_database": None, "contaminant_tags_in_database": None}
    if not rec["value"]:
        return rec
    stop = "[run_search] REFUSED -- nothing was converted, written or submitted."
    try:
        found, rec["contaminant_tags_in_database"] = ff.fasta_contaminant_summary(fasta)
    except OSError as e:
        sys.exit(f"{stop} This is a keratin sample ({rec['source']}), and {fasta} could not be "
                 f"read to check it for keratin-family contaminant entries ({e}).")
    rec["keratin_contaminants_in_database"] = len(found)
    if found:
        m = meta if isinstance(meta, dict) else {}
        src = (m.get("source") or "")
        kind, _, where = src.partition(":")
        base = (f"--path {where}" if kind == "override"
                else f"--ncbi-accession {where}" if kind == "ncbi_refseq"
                else f"--proteome {m['proteome']}" if m.get("proteome")
                else "--proteome <UPID>")
        if kind in ("override", "ncbi_refseq"):
            # the organism the user named for it (fetch_fasta.user_organism; --path needs one)
            base += (f" --organism {shlex.quote(m['organism'])}" if m.get("organism")
                     else " --organism none" if str(m.get("organism_source") or "").startswith(
                         "user: no single organism")
                     else " --organism '<scientific name>'")
            base += f" --taxid {int(m['taxid'])}" if m.get("taxid") else ""
        cont = m.get("contaminant_set")
        cont = cont if cont and cont not in ("none", "already_in_supplied_database") else "universal"
        enz = ff.enzyme_arg(m)          # the build's --enzyme, "none" too (fetch_fasta's reader)
        content = m.get("content_requested") or "one_per_gene"
        new = os.path.splitext(fasta)[0] + "_keratin.fasta"
        sys.exit(
            f"{stop} This is a keratin sample ({rec['source']}), but {fasta} still holds "
            f"{len(found)} keratin-family contaminant entr"
            f"{'y' if len(found) == 1 else 'ies'}: {ff._entry_list(found)}. Keratin is the analyte "
            f"here: "
            # what THIS engine would do with them -- only DIA-NN is given --cont-quant-exclude
            + ("the search engine (DIA-NN --cont-quant-exclude) would keep every peptide they "
               "share with the sample's keratins out of quantification, and "
               if engine == "diann" else "")
            + f"run_de.R's contaminant filter would remove those precursors from the DE. Build the "
            f"database with them removed and search that one:\n"
            f"  python3 scripts/fetch_fasta.py fetch {base} --content {content} --contaminants "
            f"{cont} --enzyme {enz} --keratin-sample --out {new}"
            + ("  [--hive as before]" if m.get("staged_file") else "")
            + (f"\n(values from {sidecar}; check them)" if m else
               f"\n({fasta} has no fetch_fasta.py sidecar: fill in the proteome and settings it "
               f"was built with)"))
    if not flag:
        print(f"[run_search] keratin sample: {sidecar} says the database was built with "
              f"--keratin-sample, so this search is recorded as one (keratin quantified as the "
              f"analyte).", file=sys.stderr)
    return rec


def _workflow_value(params, key):
    """A FragPipe workflow's value for `key` as FragPipe reads it (Java .properties: a backslash
    escapes the next character), or None when it is absent or empty."""
    try:
        with open(params, encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                if ln.startswith(key + "="):
                    v = re.sub(r"\\(.)", r"\1", ln.split("=", 1)[1].strip())
                    return v or None
    except OSError:
        pass
    return None


def fragpipe_databases(params, fasta):
    """The database(s) a FragPipe search reads: its workflow's database.db-path (make_presets.py
    writes --fasta there) and --fasta, existing files only, each once."""
    paths = [fasta, _workflow_value(params, "database.db-path")]
    return [p for p in dict.fromkeys(os.path.abspath(p) for p in paths if p) if os.path.isfile(p)]


# The workflow a FragPipe search RUNS: the one given, with database.db-path set to --fasta
# (fragpipe_workflow_for_run), in the search folder. Not FragPipe's own "fragpipe.workflow",
# which it saves there itself.
FRAGPIPE_RUN_WORKFLOW = "fragpipe_input.workflow"


def refuse_fragpipe_database(params, fasta):
    """Before anything is written (H7 bug 1): FragPipe's shipped workflows name no
    database.db-path, and headless FragPipe then takes the one cached in the account's
    ~/.config/FragPipe/.../fragpipe-ui.cache -- the LAST search's database (H5 finding 11). So
    the route sets it to --fasta itself (fragpipe_workflow_for_run) and refuses a workflow that
    names a DIFFERENT database. -> the record of the database searched (search_provenance.json
    fragpipe_database): path, sha256, entries, decoys under the workflow's decoy tag, and the
    db-path the workflow gave. A FASTA without decoys is said (MSFragger refuses one)."""
    import hashlib
    import fetch_fasta as ff
    fas = os.path.abspath(fasta)
    given = _workflow_value(params, "database.db-path")
    if given and os.path.realpath(os.path.abspath(given)) != os.path.realpath(fas):
        sys.exit(f"[run_search] REFUSED -- nothing was written or submitted. The workflow "
                 f"{params} names database.db-path={given}, but --fasta is {fas}: FragPipe "
                 f"would search the workflow's database, not the one recorded for this search. "
                 f"Give --fasta {given}, or a workflow without database.db-path (run_search.py "
                 f"sets it from --fasta).")
    tag = _workflow_value(params, "database.decoy-tag") or "rev_"
    try:
        heads = ff._fasta_headers(fas)
        h = hashlib.sha256()
        with open(fas, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                h.update(block)
        sha = h.hexdigest()
    except OSError as e:
        sys.exit(f"[run_search] REFUSED -- nothing was written or submitted. --fasta {fas} could "
                 f"not be read ({e}).")
    n_dec = sum(1 for ln, _ in heads if ln[1:].startswith(tag))
    if not n_dec:
        sys.stderr.write(f"[run_search] NOTE: {fas} has no decoys (no entry starts with {tag}, "
                         f"the workflow's database.decoy-tag): MSFragger refuses a target-only "
                         f"FASTA (\"No decoys found in the FASTA file\"). Add them with FragPipe's "
                         f"Philosopher, without --contam: philosopher database --custom {fas} "
                         f"--prefix {tag}, then give --fasta the .fas it writes.\n")
    return {"path": fas, "sha256": sha, "n_entries": len(heads), "n_decoys": n_dec,
            "decoy_tag": tag, "workflow_db_path": given,
            "source": "--fasta, written as the workflow's database.db-path by run_search.py"}


def fragpipe_workflow_for_run(params, db, out):
    """<out>/FRAGPIPE_RUN_WORKFLOW: the workflow given, line for line, with database.db-path set
    to the database recorded (refuse_fragpipe_database) -- so FragPipe never fills it from its
    cache. -> its path."""
    with open(params, encoding="utf-8", errors="replace", newline="") as fh:
        lines = fh.read().splitlines(keepends=True)
    nl = "\r\n" if lines and lines[0].endswith("\r\n") else "\n"
    value = db["path"].replace("\\", "\\\\").replace(":", "\\:").replace("=", "\\=")
    new = f"database.db-path={value}{nl}"
    out_lines, done = [], False
    for ln in lines:
        if ln.startswith("database.db-path="):
            if not done:
                out_lines.append(new)
                done = True
            continue
        out_lines.append(ln)
    if not done:
        if out_lines and not out_lines[-1].endswith(("\n", "\r")):
            out_lines[-1] += nl
        out_lines.append(new)
    path = os.path.join(out, FRAGPIPE_RUN_WORKFLOW)
    write_file(path, "".join(out_lines))
    return path


def fragpipe_contam_entries(databases):
    """[(database, number of FragPipe `contam_`-tagged entries it holds)] for the readable ones.
    The ONE count of them: the refusal before a FragPipe search and the warning when an older
    FragPipe DIA output is adapted both read it."""
    import fetch_fasta as ff
    out = []
    for path in databases:
        try:
            heads = ff._fasta_headers(path)
        except OSError:
            continue
        out.append((path, sum(1 for h, _l in heads
                              if (ff._contaminant_ids(h) or ("",))[0] == ff.FRAGPIPE_CONT_TAG)))
    return out


def refuse_fragpipe_contam_database(params, fasta):
    """A FragPipe search on a database holding `contam_` entries (built with Philosopher's --contam,
    not by fetch_fasta.py) is REFUSED before anything is written. keratin-review, 2026-09-29,
    verified on HIVE (haggerman_brains_dia_NN): FragPipe's DIA route writes library.tsv without the
    tag, so its DIA-NN report.tsv carries those contaminants as bare accessions -- nothing
    downstream (DIA-NN --cont-quant-exclude, run_de.R's filter, a keratin sample's exemption) can
    see them as contaminants, and they go into the DE. The skill's own set is tagged Cont_ inside
    the accession field, which survives both routes. FragPipe itself adds no contaminants when it
    searches: they come only from such a database."""
    import fetch_fasta as ff
    for path, n in fragpipe_contam_entries(fragpipe_databases(params, fasta)):
        if n:
            sys.exit(
                f"[run_search] REFUSED -- nothing was converted, written or submitted. {path} holds "
                f"{n} contaminant entries tagged {ff.FRAGPIPE_CONT_TAG} (a database built with "
                f"FragPipe/Philosopher's --contam, not by the skill). FragPipe's DIA route drops that "
                f"tag when it builds its library, so those contaminants reach report.tsv as plain "
                f"accessions that neither DIA-NN nor run_de.R's contaminant filter can recognise, and "
                f"they would be tested as sample proteins. Build the database with the skill's "
                f"contaminant set (tagged {ff.CONT_TAG}) instead:\n"
                f"  python3 scripts/fetch_fasta.py fetch --proteome <UPID> --contaminants universal "
                f"--out ./search.fasta [--keratin-sample | --no-keratin-sample]\n"
                # S9: that FASTA is targets only, and MSFragger refuses one (references/
                # fragpipe-diatracer.md); the decoys come from FragPipe's own Philosopher
                f"then add the decoys MSFragger needs (it refuses a target-only FASTA: \"No decoys "
                f"found in the FASTA file\") with FragPipe's Philosopher, without --contam, which "
                f"would add the {ff.FRAGPIPE_CONT_TAG} entries back:\n"
                f"  <fragpipe>/tools/Philosopher/philosopher-<v> workspace --init\n"
                f"  <fragpipe>/tools/Philosopher/philosopher-<v> database --custom ./search.fasta "
                f"--prefix rev_    # -> <date>-decoys-search.fasta.fas (targets + decoys)\n"
                f"and point the workflow's database.db-path at that target+decoy .fas:\n"
                f"  python3 scripts/make_presets.py --engine fragpipe ... --fasta "
                f"./<date>-decoys-search.fasta.fas")


# The FragPipe workflow key that runs DIA-NN as its quantification step (FragPipe 24.0's shipped
# workflows: true in DIA_SpecLib_Quant*.workflow, false in LFQ-MBR and Stellar-DDA).
FRAGPIPE_RUNS_DIANN = "diann.run-dia-nn"


def fragpipe_runs_diann(params):
    """True when the FragPipe workflow at `params` runs DIA-NN (FRAGPIPE_RUNS_DIANN=true)."""
    try:
        with open(params, encoding="utf-8", errors="replace") as fh:
            return any(ln.split("=", 1)[0].strip() == FRAGPIPE_RUNS_DIANN
                       and ln.split("=", 1)[1].strip().lower() == "true"
                       for ln in fh if "=" in ln)
    except OSError:
        return False


def refuse_fragpipe_dashed(params, files, fasta, out, acquisition):
    """A FragPipe search whose DIA-NN step would be handed a path with '--' is REFUSED before
    anything is written. FragPipe 24.0 (CmdDiann.java) writes each run as a `--f <absolute path>`
    line of filelist_diann.txt and hands DIA-NN `--cfg <filelist>--`, and DIA-NN cuts a cfg's text
    at every '--' exactly as it does its command line: on HIVE (srun 24807938) DIA-NN 2.7.0 and
    FragPipe 24.0's bundled 1.8.2 beta 8 both logged "unrecognised option [--sampA_S3.mzML]" for a
    `x_DIA--sampA_S3.mzML` in that list, skipped it, and searched the other run. The path DIA-NN
    gets is the run's link in <out>/diatracer_stage/ on the diaTracer route (dia-PASEF .d,
    run_fragpipe), so only its FILE NAME matters there; any other input goes in as given. The
    --workdir (DIA-NN's --out and the filelist live under it) and the database it is given count
    too. There is no renaming on this route: FragPipe names its own files after the run."""
    if not fragpipe_runs_diann(params):
        return
    import check_report_runs as crr
    staged = (acquisition or "").upper() == "DIA" and any(
        f.rstrip("/").endswith(".d") for f in files)
    seen = {f: (os.path.join(out, "diatracer_stage", os.path.basename(f.rstrip("/"))) if staged
                else f) for f in files}
    runs = [f for f in files if crr.dashed(seen[f])]
    paths = [(flag, p) for flag, p in [("--out", out)] + [("database", d) for d in
                                                         fragpipe_databases(params, fasta)]
             if crr.dashed(p)]
    if not (runs or paths):
        return
    lines = [f"  {flag} {p}" for flag, p in paths]
    named = crr.dashed_file_names(runs)
    lines += [ln for ln in crr.dashed_names_text(named).splitlines()]
    lines += [f"  {f}  (a folder holds '--')" for f in runs if f not in {r["file"] for r in named}]
    sys.exit(
        f"[run_search] REFUSING: this FragPipe workflow runs DIA-NN ({FRAGPIPE_RUNS_DIANN}=true in "
        f"{params}), and DIA-NN would be handed a path with '--' in it. FragPipe lists the runs "
        f"for DIA-NN in a file DIA-NN reads like a command line, cutting it at every '--', so it "
        f"would skip the run and search the others (DIA-NN 2.7.0 and FragPipe's bundled 1.8.2 "
        f"both did, HIVE srun 24807938):\n" + "\n".join(lines) + "\n"
        f"Ways out: search with DIA-NN directly (--engine diann, which links such inputs and, with "
        f"{crr.RENAME_DASHED_FLAG}, renames a run whose file name holds '--' as shown); or link "
        f"each file under a name and a folder without '--' and give FragPipe the links (use the "
        f"same names in conditions.csv); a --workdir or database with '--' needs another path. "
        f"Nothing was written or submitted.")


# Sage's own output (its log lines are on stderr), teed here by the job and by an inline run:
# sage_lfq_check.py reads its "discovered N target MS1 peaks at 5% FDR" from it.
SAGE_LOG = "sage.log"

# A Sage search for IDENTIFICATIONS ONLY: its cfg says "quant": {"lfq": false} (written by
# estimate_params.py --ids-only). Sage 0.14.7 reads a timsTOF .d itself for identifications but has
# no MS1 quantification on one, so this is the route that searches a .d without converting it:
# A Core submission on HIVE, 2026-10-30, 896 peptides at 1% from the two .d in 22 s against 753 from an
# 11.6 GB Windows mzML conversion (sage015_test/RESULTS.md, arm B). Decided by the cfg key alone.
# The search gives identifications.tsv (SAGE_IDS, imported from sage_lfq_check; adapt_sage_ids)
# and no report.parquet, so no DE.
SAGE_IDS_RECORD = "sage_identifications.json"
# A run's identification (L7R-F2): a target PSM of rank 1 at Sage's spectrum_q AND peptide_q at 1%.
# peptide_q alone is one value per peptide over all the runs, so a peptide confidently seen in
# one run would be listed as identified in another where its only spectrum is poor; and a rank-2
# (chimeric) PSM is not the spectrum's identification. This is how sage015_test/RESULTS.md counts.
SAGE_IDS_Q = 0.01
SAGE_IDS_RULE = (f"a target PSM (is_decoy false) of rank 1 with Sage's spectrum_q <= "
                 f"{SAGE_IDS_Q:g} and peptide_q <= {SAGE_IDS_Q:g}")
# ...and the same rule in the words of a Methods section (sage_identifications.methods_sentence,
# which make_methods.py quotes)
# (the q-values are global: a run's subset of a globally filtered list carries no 1% FDR of its
# own, so the sentence says "over all runs together, listed per run" -- C2 statistician nit)
SAGE_IDS_METHODS = (f"Identifications were Sage's own at {SAGE_IDS_Q * 100:g}% FDR over all runs "
                    f"together, listed per run where a peptide-spectrum match passed: a target "
                    f"peptide-spectrum match of rank 1 with spectrum- and peptide-level q-values "
                    f"<= {SAGE_IDS_Q:g}.")
SAGE_IDS_COLUMNS = ("run", "peptide", "proteins", "charge", "spectrum_q", "peptide_q",
                    "protein_q")
# No protein-level filter is applied to the table (C2 statistician H1 + proteomics L1: on the
# two-run DDA test the kept PSMs map to 198 protein groups, of which 143 pass Sage's own 1%
# picked-protein FDR), so the Methods say so, with that count from the table itself
# (sage_ids_protein_groups; the record's totals protein_groups / protein_groups_protein_q).
SAGE_IDS_PROTEIN_METHODS = (
    "No protein-level FDR filter was applied to the identification list: each peptide is listed "
    "with every database protein that contains it, and with Sage's picked-protein q-value; "
    "{passing} of the {groups} protein groups listed have a protein-level q-value <= {q:g}.")


def sage_ids_protein_groups(rows):
    """(protein groups listed, of them with a Sage protein_q <= SAGE_IDS_Q) over
    identifications.tsv rows (SAGE_IDS_COLUMNS order): a group is Sage's `proteins` string, at
    its best protein_q."""
    best = {}
    for r in rows:
        q = r[6]
        if r[2] not in best or (q is not None and (best[r[2]] is None or q < best[r[2]])):
            best[r[2]] = q
    return len(best), sum(1 for q in best.values() if q is not None and q <= SAGE_IDS_Q)


def sage_ids_methods_sentence(groups, passing):
    """The Methods' sentences for a Sage identifications-only search: the PSM/peptide rule
    (SAGE_IDS_METHODS) and that no protein-level filter was applied, with the count."""
    return (SAGE_IDS_METHODS + " " + SAGE_IDS_PROTEIN_METHODS.format(
        passing=passing, groups=groups, q=SAGE_IDS_Q))


def _is_dotd(path):
    """A Bruker timsTOF .d (a folder named <run>.d)."""
    return str(path).rstrip("/").lower().endswith(".d")


def _sage_cfg_lfq(params):
    """(the cfg's quant.lfq when it is a bool, else None; what the cfg says, in words)."""
    try:
        with open(params) as fh:
            cfg = json.load(fh)
    except (OSError, ValueError, TypeError) as e:
        return None, f"the cfg {params} could not be read ({type(e).__name__})"
    quant = cfg.get("quant") if isinstance(cfg, dict) else None
    lfq = quant.get("lfq") if isinstance(quant, dict) else None
    if isinstance(lfq, bool):
        return lfq, f"quant.lfq is {'true' if lfq else 'false'} in the cfg"
    return None, "quant.lfq is not set (true or false) in the cfg"


def sage_lfq_setting(params):
    """(whether this Sage search quantifies, how that was decided). Only a cfg whose quant.lfq is
    false is identifications only. Not set, or a cfg that cannot be read, counts as quantifying
    -- the DEFAULT, said so -- for a search with no .d; with a .d, refuse_sage_dotd_route stops
    it first (L7R-F6: Sage's own default is no LFQ, so the route would be guessed)."""
    lfq, said = _sage_cfg_lfq(params)
    if lfq is None:
        return True, f"DEFAULT -- {said}, so the search is treated as quantifying"
    return lfq, said


def refuse_sage_dotd_route(bundle, params, files):
    """Stop, before anything is converted, written or submitted, a Sage search of a timsTOF .d
    whose route is not decided, or does not fit the data:
      * the cfg's quant.lfq is not true or false (L7R-F6): it decides whether each .d is
        converted to mzML (hours, tens of GB) or read as it is, and Sage's own default (no LFQ)
        is not what the skill would otherwise assume;
      * a quantifying search of a diaPASEF .d (L7R-F7): the conversion's settings
        (PWIZ_TIMS_MSCONVERT) are DDA-PASEF's -- its scan summing would add up each isolation
        window's MS2 over about ten cycles, which DIA scoring cannot use."""
    dots = [f for f in files if _is_dotd(f)]
    if not dots:
        return
    lfq, said = _sage_cfg_lfq(params)
    if lfq is None:
        sys.exit(f"[run_search] REFUSED -- Nothing was converted, written or submitted. {said}, "
                 f"and it decides how a timsTOF .d goes into Sage ({len(dots)} input(s) are .d). "
                 f"Set quant.lfq true (each .d converted to mzML and quantified) or false "
                 f"(identifications only: Sage reads the .d itself, no DE -- estimate_params.py "
                 f"--engine sage --ids-only writes it) in {params}.")
    acq = str(bundle.get("acquisition") or "").upper()
    if lfq and acq == "DIA":
        sys.exit(f"[run_search] REFUSED -- Nothing was converted, written or submitted. This is a "
                 f"quantifying Sage search of timsTOF .d whose acquisition is DIA (the workflow "
                 f"manifest). The .d conversion's settings are DDA-PASEF's (msconvert "
                 f"{pwiz_msconvert_options()}): their scan summing would add up each diaPASEF "
                 f"window's MS2 across cycles. Search diaPASEF with DIA-NN (--engine diann).")


def sage_lfq_enabled(params):
    """The cfg's quant.lfq, default True (an unreadable cfg too): see sage_lfq_setting."""
    return sage_lfq_setting(params)[0]


def pwiz_methods_sentence(sif):
    """The Methods sentence for a .d converted in this container, from PWIZ_TIMS_MSCONVERT and the
    container record -- a version that is not recorded is said to be so (L7R-F3)."""
    c = pwiz_container_record(sif)
    ver = c.get("msconvert_version")
    prog = (f"ProteoWizard msconvert {ver}" if ver else
            "ProteoWizard msconvert (version not recorded -- confirm)")
    where = (f"the ProteoWizard container {c['image']}:{c['tag']}" if c.get("tag") else
             f"the container {os.path.basename(sif)} (not the Core's pinned image)")
    return (f"Each timsTOF .d file was converted to mzML with {prog} in {where}, under Wine, with "
            f"vendor peak picking at MS1 and MS2, ion-mobility spectra combined, and MS2 spectra "
            f"of similar precursor m/z, retention time and ion mobility summed (msconvert options: "
            f"{pwiz_msconvert_options()}).")


def sage_input_record(files, params, steps):
    """search_provenance.json `sage_input`: how each timsTOF .d went into Sage, and why. `dotd` is
    "direct" (searched as it is), "converted" (to mzML first) or null when no input is a .d.
    `methods_sentence` says the route in a Methods section's words (make_methods.py quotes it;
    None when there is nothing to add)."""
    dots = [f for f in files if _is_dotd(f)]
    quantified, how = sage_lfq_setting(params)
    ids_only = "for identifications only, with no MS1 quantification"
    if not dots:
        return {"dotd": None, "files": [], "why": "no timsTOF .d among the inputs",
                "methods_sentence": None if quantified else f"The search was {ids_only}."}
    if not quantified:
        return {"dotd": "direct", "files": dots,
                "why": (f"identifications only ({how}): Sage 0.14.7 reads a timsTOF .d itself "
                        f"for identifications, so each .d was searched as it is, with no mzML "
                        f"conversion"),
                "methods_sentence": (f"Sage read each timsTOF .d file directly, without "
                                     f"conversion to mzML; the search was {ids_only}.")}
    tools = sorted({s["tool"] for s in steps if s["input"] in dots})
    sifs = list(dict.fromkeys(s["sif"] for s in steps if s.get("sif")))
    return {"dotd": "converted", "files": dots, "tools": tools,
            "why": (f"quantification ({how}): Sage 0.14.7 has no MS1 quantification on a .d, "
                    f"so each .d was converted to mzML first, with {'; '.join(tools)}: "
                    f"msconvert {pwiz_msconvert_options()} ({PWIZ_TIMS_MSCONVERT_SOURCE})"),
            "methods_sentence": " ".join(pwiz_methods_sentence(sif) for sif in sifs) or None}


# The pipeline env's interpreter, as acquire_tools.sh finds it -- left for the JOB's shell to
# expand (the cluster's $HOME, not the generator's). ${HOME:-}: a job without $HOME under set -u
# falls through to the next choice instead of ending.
PIPELINE_ENV_PYTHON = ("${PROTEOMICS_PIPELINE_HOME:-${HOME:-}/.proteomics-pipeline}"
                       "/micromamba/envs/proteomics-pipeline/bin/python")


def _job_python():
    """A job's interpreter for the skill's own scripts (the LFQ check, the Radiant library): the
    first executable of the pipeline env's python (which has pyarrow), this generator's, then
    python3. Not this one alone: a search generated under /usr/bin/python3 gave the job no
    pyarrow, and its Sage LFQ check could only say "unchecked" (a Core submission, 2026-10-29)."""
    return (f'_py="{PIPELINE_ENV_PYTHON}"; [ -x "$_py" ] || '
            f'_py={shlex.quote(sys.executable or "python3")}; [ -x "$_py" ] || _py=python3')


def sage_check_line(out, params):
    """The job's LFQ check, right after Sage. Never fails the job (Sage succeeded); if the check
    itself cannot run it says so, and --adapt-only runs it again."""
    q = shlex.quote
    check = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sage_lfq_check.py")
    return (f'"$_py" {q(check)} --out {q(out)} --params {q(params)} '
            f'--log {q(os.path.join(out, SAGE_LOG))} > /dev/null '
            '|| echo "[sage_lfq_check] NOTE: the LFQ check itself failed (exit $?); '
            'run_search.py --adapt-only runs it again" >&2')


def sh_logged(cmd, log):
    """sh(), with everything the command prints also written to `log`."""
    print(f"  $ {cmd}", flush=True)
    with open(log, "w") as fh:
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, errors="replace")
        for line in p.stdout:
            sys.stdout.write(line)
            fh.write(line)
        rc = p.wait()
    if rc:
        raise subprocess.CalledProcessError(rc, cmd)


def sage_lfq(out, params, log=None, accept=None):
    """Run the LFQ mass-window check (sage_lfq_check.py), record it, say it. The record. A record
    the user accepted before is kept accepted while nothing it accepted has changed; `accept`
    (--accept-lfq-window's reason) accepts a gated one now."""
    import sage_lfq_check
    rec = sage_lfq_check.carry_acceptance(sage_lfq_check.load(out),
                                          sage_lfq_check.check(out, params, log))
    if accept and sage_lfq_check.gated(rec):
        sage_lfq_check.accept(out, rec, accept)
    try:
        rec["written_to"] = sage_lfq_check.record(out, rec)
    except OSError as e:
        sys.stderr.write(f"{sage_lfq_check.TAG} NOTE: could not write "
                         f"{sage_lfq_check.RECORD} in {out}: {e}\n")
    sage_lfq_check.say(rec)
    return rec


def sage_ids_lfq_record(out):
    """An identifications-only search's LFQ check: not applicable, written where a quantifying
    search's record would be (sage_lfq_check.json + search_provenance.json), so the watcher and
    the audit never read an earlier search's record as this one's. Called only once
    adapt_sage_ids has accepted the output as identifications only."""
    rec = sage_lfq_check_mod().not_applicable()
    try:
        rec["written_to"] = sage_lfq_check_mod().record(out, rec)
    except OSError as e:
        sys.stderr.write(f"{sage_lfq_check_mod().TAG} NOTE: could not write "
                         f"{sage_lfq_check_mod().RECORD} in {out}: {e}\n")
    return rec


def sage_set_aside_other_kind(out, quantified):
    """Set aside (renamed .stale-<time>, never deleted) what an earlier Sage search of the OTHER
    kind left in `out`, when this one is generated -- as run_diann's previous_artefacts_moved_to
    (L7R-F9). Until --adapt-only they described the folder: the watcher showed an earlier LFQ
    record, and run_de.R could gate on it and run a DE on an earlier report.parquet beside a
    provenance saying quantified: false. Identifications only: the quantities and their records
    go, and the LFQ record is written now (sage_ids_lfq_record). Quantifying: identifications.tsv
    and its record go, and a not-applicable LFQ record; a quantifying search's own LFQ record
    stays (the job's check carries an earlier acceptance over). Called only once every refusal
    has passed. -> the new paths."""
    slc = sage_lfq_check_mod()
    if quantified:
        names = [SAGE_IDS, SAGE_IDS_RECORD]
        if (slc.load(out) or {}).get("status") == "not_applicable":
            names.append(slc.RECORD)
    else:
        names = ["report.parquet", "lfq.parquet", SAGE_ADAPT, slc.RECORD]
    moved = []
    for name in names:
        try:
            new = set_aside(os.path.join(out, name))
        except ValueError as e:
            sys.exit(f"[run_sage] {e}. Move it yourself, or search into a different --out. "
                     f"Nothing was submitted.")
        if new:
            moved.append(new)
    if moved:
        print(f"  [run_sage] an earlier search's outputs set aside: {', '.join(moved)}")
    return moved


def _load_sage_record(out, name):
    try:
        with open(os.path.join(out, name)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None                     # the step that writes it said why it could not


def run_sage(cmd, params, files, fasta, out, threads, sbatch, queue=None, pwiz_sif=None):
    quantified = sage_lfq_enabled(params)
    # mzml_plan creates nothing and its refusals say "Nothing was converted, written or
    # submitted": the folder is made only once it has passed, so a refusal leaves no empty
    # search folder behind (H2c N1)
    mzml, steps = mzml_plan(files, out, dotd="convert" if quantified else "direct",
                            pwiz_sif=pwiz_sif)
    os.makedirs(out, exist_ok=True)
    files_args = " ".join(shlex.quote(m) for m in mzml)
    full = (f"{cmd} {shlex.quote(params)} -f {shlex.quote(fasta)} -o {shlex.quote(out)} "
            f"--parquet --disable-telemetry-i-dont-want-to-improve-sage {files_args}")
    log = os.path.join(out, SAGE_LOG)
    said = {"quantified": quantified, "sage_input": sage_input_record(files, params, steps),
            "previous_artefacts_moved_to": sage_set_aside_other_kind(out, quantified)}
    if not quantified:
        said["sage_lfq_check"] = sage_ids_lfq_record(out)      # this search's, from now on
    sifs = [s["sif"] for s in steps if s.get("sif")]
    if sifs and pwiz_rerun_args(sifs[0]):
        said["pwiz_sif_arg"] = pwiz_rerun_args(sifs[0])      # a re-run keeps this container
    if sbatch:
        # Conversion, search and LFQ check all run in the job, on the compute node: nothing
        # heavy happens before submission (golden rule 3). pipefail keeps Sage's exit status.
        # Identifications only: no LFQ, so no LFQ check.
        job = "\n".join(conversion_lines(steps)
                        + [f"{full} 2>&1 | tee {shlex.quote(log)}"]
                        + ([_job_python(), sage_check_line(out, params)] if quantified else []))
        emit_sbatch(sbatch, job, out, threads, job="sage_search", **(queue or {}))
        res = {"engine": "sage", "out": out, "submitted": sbatch, "ran": False,
               "mzml_conversion": conversion_record(steps, "in the job, before Sage"), **said}
        if quantified:
            res["note"] = ("After the job runs, re-run with --adapt-only to build report.parquet "
                           "(it also re-runs the LFQ mass-window check, sage_lfq_check.json).")
        else:
            res["note"] = (f"Identifications only (quant.lfq is false): after the job runs, "
                           f"re-run with --adapt-only to write {SAGE_IDS}. There is no "
                           f"report.parquet, so no differential expression.")
        return res
    convert_inline(steps)
    sh_logged(full, log)
    if not quantified:
        ids = adapt_sage_ids(out)
        return {"engine": "sage", "report": None, "identifications": ids, "ran": True,
                "mzml_conversion": conversion_record(steps, "inline, before Sage"), **said,
                "sage_identifications": _load_sage_record(out, SAGE_IDS_RECORD)}
    rec = sage_lfq(out, params, log, accept=getattr(a_globals, "accept_lfq_window", None))
    if sage_lfq_check_mod().gated(rec):
        # the gate (sage_lfq_check.refusal): no report.parquet from quantities known to be wrong;
        # main() records the search, then says how to re-run it and exits non-zero
        return {"engine": "sage", "report": None, "ran": True, "lfq_window_refused": True,
                "mzml_conversion": conversion_record(steps, "inline, before Sage"), **said,
                "sage_lfq_check": rec}
    report = adapt_sage(out)
    rec = sage_lfq_check_mod().load(out) or rec          # with the quantified coverage (M2)
    return {"engine": "sage", "report": report, "ran": True,
            "mzml_conversion": conversion_record(steps, "inline, before Sage"), **said,
            "sage_lfq_check": rec, "sage_adapt": _load_sage_record(out, SAGE_ADAPT)}


SAGE_ADAPT = "sage_adapt.json"


def sage_lfq_check_mod():
    import sage_lfq_check
    return sage_lfq_check


def sage_version(out):
    """The Sage version for the methods: the one run_search.py recorded for the engine it ran,
    else Sage's own results.json -- whose number is the CRATE version: the v0.14.7 release
    writes "0.14.6" (crates/sage-cli/Cargo.toml at tag v0.14.7), so that one is said as such."""
    ver = recorded_engine_version(out)
    if ver:
        return ver
    try:
        with open(_find(out, ["results.json"]) or "") as fh:
            rj = json.load(fh).get("version")
    except (OSError, ValueError, AttributeError):
        return None
    if rj == "0.14.6":
        return "0.14.6 per its results.json (the v0.14.7 release also reports 0.14.6)"
    return rj or None


def sage_quantity_declaration(out, qmax):
    """What adapt_sage's rows are, in words the DE's methods can quote."""
    return {
        "level": "peptide",
        "label": "Sage",
        "value": f"{_vtag('Sage', sage_version(out))} label-free MS1 peptide intensity "
                 f"(lfq.parquet)",
        "identification": (f"Sage's own FDR -- peptides passing its 1% peptide-level q-value "
                           f"enter its LFQ, and MS1 peaks are kept at a stored LFQ q_value <= "
                           f"{qmax:g} with its decoy peaks dropped (run_search.py adapt_sage; the "
                           f"stored q is shared with the precursor's decoy, so this keeps a "
                           f"conservative subset of Sage's own 5% count). {PLACEHOLDER_Q}"),
        "q_columns": "placeholder",
        "citation": CITE["sage"],
    }


# An empty Sage result is refused, never written as one (H2-2, rule 4): on HIVE an un-combined
# timsTOF mzML gave Sage 0 peptides at 1%, and --adapt-only exited 0 with a 0-row report.parquet
# and no warning. Exit 2, so `... --adapt-only && run_de.R ...` stops there.
SAGE_EMPTY_EXIT = 2


def _sage_refuse_empty(out, product, record, key, rec, reason):
    """Refuse an empty adapter result: no <product> is written, and one already in `out` is set
    aside (it is not this search's, and a DE or a delivery would read it as one); the reason is
    recorded in <record> + search_provenance.json `key` (status "refused"), said, and the
    process exits SAGE_EMPTY_EXIT."""
    path = os.path.join(out, product)
    try:
        moved = set_aside(path)
    except ValueError as e:
        moved = None
        reason += f" ({path} is not a regular file, so it was left where it is: {e})"
    rec.update(status="refused", refused=reason, earlier_output_moved_to=moved)
    try:
        rec["written_to"] = sage_lfq_check_mod().write_record(out, record, key, rec)
    except OSError as e:
        sys.stderr.write(f"[adapt] NOTE: could not write {record} in {out}: {e}\n")
    sys.stderr.write(f"[adapt] REFUSED: {reason}. No {product} was written"
                     + (f"; the earlier one was moved to {moved}" if moved else "")
                     + f" (recorded in {record}).\n")
    sys.exit(SAGE_EMPTY_EXIT)


def adapt_sage(out):
    """Sage lfq.parquet -> the DE-input report.parquet, one row per (peptide, file) that Sage
    counts as a valid MS1 quantity. Sage quantifies PEPTIDES; build_maxlfq.R makes the protein x
    run cell from these rows downstream.

    lfq.parquet is NOT filtered by Sage. Its writer (sage-cloudpath parquet.rs serialize_lfq,
    v0.14.7) writes every traced MS1 peak: each target AND its decoy -- the same peptide measured
    in a +11.06 Da-shifted window (lfq.rs build_feature_map), with the SAME `proteins` string,
    so an unfiltered decoy row reads as that protein's quantity -- at any q_value. Sage's own
    lfq.tsv writer (sage-cli output.rs write_lfq) drops the decoys and filters on no q_value;
    Sage counts an MS1 peak as discovered only at q_value <= 0.05 (fdr.rs picked_precursor,
    logged "target MS1 peaks at 5% FDR"). So a row is kept when is_decoy is false (as write_lfq
    keeps it) AND its stored q_value <= sage_lfq_check.LFQ_Q_MAX. The q_value Sage stores is
    shared with the precursor's decoy (the lower-scoring member's q -- see LFQ_Q_MAX), so this
    keeps a conservative subset of Sage's own count; Sage's logged count is recorded beside what
    was kept. Until 2.9 every row went through, decoys included (gabrig 2026-09-29: a cohort whose
    targets all sat at q 0.128 still produced a full matrix). What was kept and dropped, per file,
    goes to <out>/sage_adapt.json and search_provenance.json `sage_adapt`.
    """
    try:
        import pyarrow.parquet as pq
        import pyarrow as pa
    except ImportError:
        sys.exit("pyarrow required to adapt Sage output. pip install pyarrow.")
    import sage_lfq_check

    lfq = _find(out, ["lfq.parquet"])
    if not lfq:
        if sage_ran_lfq(out) is False:
            sys.exit(f"No lfq.parquet under {out}: this Sage search was for identifications only "
                     f"(its results.json says quant.lfq false), so there is nothing to build "
                     f"report.parquet from. Its identifications go to {SAGE_IDS}: run_search.py "
                     f"--adapt-only with the cfg the search ran with (quant.lfq false) writes it.")
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
    c_dec, c_q = col("is_decoy"), col("q_value")
    if not all([c_prot, c_run, c_int]):
        sys.exit(f"Sage lfq.parquet missing expected columns; saw {t.column_names}")
    if not (c_dec and c_q):
        # every Sage lfq.parquet since v0.14.0 has both: without them decoys cannot be told
        # from targets, and adapting anyway is how they reached the DE before
        sys.exit(f"{lfq} has no is_decoy / q_value column ({t.column_names}), so its decoy MS1 "
                 "peaks cannot be told from the targets. Nothing was adapted. Re-run the search "
                 "with Sage >= 0.14.")

    # Nothing identified at 1% (sage_lfq_check's one count): nothing to quantify, whatever
    # lfq.parquet holds -- refused, not adapted (H2-2)
    res = _find(out, ["results.sage.parquet"])
    peptides, why_not = None, (f"no results.sage.parquet under {out}" if not res else None)
    if res:
        try:
            peptides = sage_lfq_check.per_run_mass_error(res)[1]
        except Exception as e:      # noqa: BLE001 -- said in the record, never dropped
            why_not = f"{res} unreadable ({type(e).__name__}: {e})"
        if peptides == 0:
            _sage_refuse_empty(out, "report.parquet", SAGE_ADAPT, "sage_adapt",
                               {"source": lfq, "report": None, "sage_results": res,
                                "target_peptides_1pct": 0},
                               f"Sage identified no target peptide at "
                               f"{sage_lfq_check.PSM_Q * 100:g}% FDR (peptide_q <= "
                               f"{sage_lfq_check.PSM_Q}) in {res}, so this search has nothing to "
                               f"quantify")
        if peptides is None:
            sys.stderr.write(f"[adapt] NOTE: the target peptides at 1% could not be counted: "
                             f"{why_not}\n")

    qmax = sage_lfq_check.LFQ_Q_MAX
    per_file, runs, prots, ints, names, feats = {}, [], [], [], [], []
    # the precursors kept, to set beside Sage's own logged count (one q per precursor, not per row)
    n_rows = t.num_rows
    peps = t.column(col("peptide")).to_pylist() if col("peptide") else [None] * n_rows
    chgs = t.column(col("charge")).to_pylist() if col("charge") else [None] * n_rows
    kept_precursors = set()
    for p, r, v, dec, ok, pep, chg in zip(
            t.column(c_prot).to_pylist(), t.column(c_run).to_pylist(),
            t.column(c_int).to_pylist(), t.column(c_dec).to_pylist(),
            sage_lfq_check.q_passes(t.column(c_q)).to_pylist(), peps, chgs):
        rn = sage_run_name(r)
        n = per_file.setdefault(rn, {"rows": 0, "decoy_dropped": 0, "q_value_dropped": 0,
                                     "kept": 0})
        n["rows"] += 1
        if dec:
            n["decoy_dropped"] += 1
            continue
        if not ok:
            n["q_value_dropped"] += 1
            continue
        n["kept"] += 1
        kept_precursors.add((pep, chg))
        runs.append(rn)
        # Sage's `proteins` is a ';' list of full FASTA IDs (sp|Cont_P02769|ALBU_BOVIN): written
        # as DIA-NN writes Protein.Group -- bare accessions -- so the contaminant filter, whose
        # rule is "an accession starts with Cont_/contam_", sees them. Written as Sage gave it,
        # no Sage DE ever had a contaminant removed (gabrig HeL50 UnvPe: 171 Cont_ groups tested,
        # methods "Contaminants : none"). The entry names go where DIA-NN puts them.
        prots.append(group_accessions(p))
        names.append(group_entry_names(p))
        # the item the DE counts (contaminants.R CONTAMINANT_FEATURE_COLUMNS): Sage's quantified
        # peptide -- its charge only when Sage kept charge states apart (combine_charge_states
        # false). Without it every peptide x run row was counted as a "precursor".
        feats.append(str(pep) if chg is None else f"{pep}/{chg}")
        ints.append(float(v) if v is not None else float("nan"))

    k = len(prots)
    out_tbl = pa.table({
        "Run": runs, "Protein.Group": prots, "Protein.Names": names, "Peptide.Id": feats,
        "PG.MaxLFQ": ints,
        # The rows are already Sage's valid MS1 quantities (decoys out, q_value <= LFQ_Q_MAX), so
        # the contract's q-columns are 0.0: the DE's own q filter must not apply a second,
        # different threshold to an MS1-integration q-value that none of them describes.
        **{c: [0.0] * k for c in FDR_REQUIRED},
    }, schema=pa.schema([("Run", pa.string()), ("Protein.Group", pa.string()),
                         ("Protein.Names", pa.string()), ("Peptide.Id", pa.string()),
                         ("PG.MaxLFQ", pa.float64())]
                        + [(c, pa.float64()) for c in FDR_REQUIRED]))
    # PG.MaxLFQ here is the DE contract's column NAME, not what the numbers are: the report says
    # what they are, and build_maxlfq.R's methods text reads it from there.
    declared = sage_quantity_declaration(out, qmax)
    out_tbl = declare_quantity(out_tbl, declared)
    report = os.path.join(out, "report.parquet")
    totals = {key: sum(f[key] for f in per_file.values())
              for key in ("rows", "decoy_dropped", "q_value_dropped", "kept")}
    if not k:
        _sage_refuse_empty(out, "report.parquet", SAGE_ADAPT, "sage_adapt",
                           {"source": lfq, "report": None, "q_value_max": qmax,
                            "per_file": per_file, "totals": totals},
                           f"no Sage MS1 peak passed q_value <= {qmax:g} ({totals['rows']} "
                           f"lfq.parquet rows: {totals['decoy_dropped']} decoy, "
                           f"{totals['q_value_dropped']} above it), so this search has nothing "
                           f"to quantify -- see sage_lfq_check.json")
    pq.write_table(out_tbl, report)
    logp = sage_lfq_check.find_log(out)
    rec = {"source": lfq, "report": report,
           "kept_rule": (f"is_decoy false and stored q_value <= {qmax:g} (the stored q is shared "
                         f"with the precursor's decoy, so this is a conservative subset of "
                         f"Sage's own count)"),
           "decoy_rule_source": "Sage 0.14.7 sage-cli output.rs write_lfq (lfq.tsv) drops "
                                "decoys; lfq.parquet keeps them",
           "q_value_max": qmax,
           "q_value_source": ("Sage 0.14.7 fdr.rs picked_precursor / main.rs: an MS1 peak is "
                              "discovered at its own q_value <= 0.05 ('target MS1 peaks at 5% "
                              "FDR'); the q it then STORES is keyed by precursor, shared by the "
                              "target and its decoy (the lower-scoring member's), so filtering "
                              "it keeps a conservative subset of that count"),
           "target_precursors_kept": len(kept_precursors),
           # the protein groups report.parquet has amounts for (sage_lfq_check's coverage)
           "protein_groups_kept": len(set(prots)),
           "sage_logged_target_ms1_peaks_5pct": (sage_lfq_check.ms1_peaks_from_log(logp)
                                                 if logp else None),
           "target_peptides_1pct": peptides,
           "target_peptides_1pct_not_counted": why_not,
           "sage_log": logp,
           "row_unit": "one row per (peptide, file); protein rollup downstream",
           # the rule above in a Methods section's words (make_methods.py quotes it)
           "methods_sentence": sage_lfq_check.methods_sentence(qmax),
           "declared_quantity": declared,
           "per_file": per_file, "totals": totals}
    try:
        rec["written_to"] = sage_lfq_check.write_record(out, SAGE_ADAPT, "sage_adapt", rec)
    except OSError as e:
        sys.stderr.write(f"[adapt] NOTE: could not write {SAGE_ADAPT} in {out}: {e}\n")
    # how much of what Sage identified this report quantifies: into the LFQ check's record, a
    # WARNING when low (C2 proteomics M2; a warning, never a gate)
    try:
        sage_lfq_check.add_coverage(out, rec)
    except OSError as e:
        sys.stderr.write(f"[adapt] NOTE: could not record the quantified coverage in "
                         f"{sage_lfq_check.RECORD}: {e}\n")
    cov = sage_lfq_check.quantified_coverage(rec)
    logged = rec["sage_logged_target_ms1_peaks_5pct"]
    print(f"  [adapt] Sage -> {report}  ({k} peptide x run rows kept of {totals['rows']}: "
          f"{totals['decoy_dropped']} decoy, {totals['q_value_dropped']} above stored q_value "
          f"{qmax:g} dropped; {len(kept_precursors)} target precursors kept"
          + (f" of the {logged} Sage logged at 5% FDR" if logged is not None else "")
          + (f", for {cov['target_peptides_identified_1pct']} target peptides identified at 1% "
             f"FDR ({cov['fraction'] * 100:.1f}%)" if cov else "")
          + f"; per file in {SAGE_ADAPT})")
    return report


def sage_ran_lfq(out):
    """quant.lfq as Sage ran this search (its results.json holds every parameter as run), or
    None when there is no readable results.json."""
    try:
        with open(_find(out, ["results.json"]) or "") as fh:
            quant = json.load(fh).get("quant")
    except (OSError, ValueError, AttributeError):
        return None
    lfq = quant.get("lfq") if isinstance(quant, dict) else None
    return lfq if isinstance(lfq, bool) else None


def adapt_sage_ids(out):
    """Sage results.sage.parquet -> <out>/identifications.tsv, for a search for identifications
    only (quant.lfq false). No report.parquet: there are no quantities, so there is no DE.

    Kept: SAGE_IDS_RULE -- a target PSM of rank 1 at spectrum_q and peptide_q <= SAGE_IDS_Q (Sage
    computes its q-values over all the runs together). Written: one row per (run, modified
    peptide, charge), at the best spectrum_q of its PSMs -- `run` is the file name without its
    extension, as adapt_sage names it; `proteins` is Sage's own ';' list of FASTA IDs. What was
    kept and dropped, per run, goes to <out>/sage_identifications.json and search_provenance.json
    `sage_identifications`. A report.parquet already in <out> is from an earlier search (this one
    wrote none), so it is set aside: the DE would read it as this search's."""
    try:
        import pyarrow.parquet as pq
    except ImportError:
        sys.exit("pyarrow required to read Sage's results.sage.parquet. pip install pyarrow.")
    import csv
    import sage_lfq_check

    if sage_ran_lfq(out) is True:
        sys.exit(f"[adapt] {out} holds a QUANTIFIED Sage search (its results.json says quant.lfq "
                 f"true), but the cfg given says identifications only (quant.lfq false). Nothing "
                 f"was written or moved. Re-run --adapt-only with the cfg the search ran with.")
    res = _find(out, ["results.sage.parquet"])
    if not res:
        sys.exit(f"No results.sage.parquet under {out}; did the Sage search finish?")
    need = ("filename", "peptide", "proteins", "charge", "rank", "is_decoy", "spectrum_q",
            "peptide_q", "protein_q")
    have = pq.read_schema(res).names
    missing = [c for c in need if c not in have]
    if missing:
        sys.exit(f"{res} has no {', '.join(missing)} column ({have}); it is not Sage 0.14 output. "
                 f"Nothing was written.")
    cols = list(need) + (["stripped_peptide"] if "stripped_peptide" in have else [])
    t = pq.read_table(res, columns=cols).to_pydict()
    stripped = t.get("stripped_peptide") or t["peptide"]
    qmax = SAGE_IDS_Q
    best, per_run, peptides = {}, {}, set()

    def passes(v):
        return v is not None and v <= qmax

    for i, fn in enumerate(t["filename"]):
        run = sage_run_name(fn)
        n = per_run.setdefault(run, {"psms": 0, "decoy_dropped": 0, "rank_dropped": 0,
                                     "spectrum_q_dropped": 0, "peptide_q_dropped": 0,
                                     "psms_kept": 0, "precursors": 0})
        n["psms"] += 1
        # each PSM is counted under the first rule it fails, in SAGE_IDS_RULE's order
        if t["is_decoy"][i]:
            n["decoy_dropped"] += 1
            continue
        if t["rank"][i] != 1:
            n["rank_dropped"] += 1
            continue
        if not passes(t["spectrum_q"][i]):
            n["spectrum_q_dropped"] += 1
            continue
        q = t["peptide_q"][i]
        if not passes(q):
            n["peptide_q_dropped"] += 1
            continue
        n["psms_kept"] += 1
        peptides.add(stripped[i])
        key = (run, t["peptide"][i], t["charge"][i])
        row = (run, t["peptide"][i], t["proteins"][i], t["charge"][i], t["spectrum_q"][i], q,
               t["protein_q"][i])
        if key not in best:
            n["precursors"] += 1
        if key not in best or row[4] < best[key][4]:
            best[key] = row

    if not best:
        _sage_refuse_empty(out, SAGE_IDS, SAGE_IDS_RECORD, "sage_identifications",
                           {"source": res, "identifications": None, "kept_rule": SAGE_IDS_RULE,
                            "per_run": per_run, "quantified": False},
                           f"no Sage PSM in {res} is {SAGE_IDS_RULE} (of "
                           f"{sum(r['psms'] for r in per_run.values())} PSMs read), so this "
                           f"search has no identifications to list")
    report = os.path.join(out, "report.parquet")
    moved = None
    if os.path.lexists(report):
        try:
            moved = set_aside(report)
        except ValueError as e:
            sys.exit(f"[adapt] {report} is from an earlier search (this one wrote none) and "
                     f"could not be set aside: {e}. Move it away, then re-run --adapt-only.")
        sys.stderr.write(f"[adapt] NOTE: {report} was left by an earlier search -- this one is "
                         f"for identifications only and writes none -- so it was moved to "
                         f"{moved}, where no DE reads it.\n")

    def cell(v):
        return "" if v is None else (f"{v:.6g}" if isinstance(v, float) else str(v))

    path = os.path.join(out, SAGE_IDS)
    with open_text(path) as fh:                     # written whole, or not at all
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(SAGE_IDS_COLUMNS)
        for key in sorted(best, key=lambda k: (k[0], str(k[1]), k[2] or 0)):
            w.writerow([cell(v) for v in best[key]])
    totals = {k: sum(r[k] for r in per_run.values())
              for k in ("psms", "decoy_dropped", "rank_dropped", "spectrum_q_dropped",
                        "peptide_q_dropped", "psms_kept")}
    totals.update(rows=len(best), precursors=len({k[1:] for k in best}), peptides=len(peptides))
    groups, passing = sage_ids_protein_groups(best.values())
    totals.update(protein_groups=groups, protein_groups_protein_q=passing)
    rec = {"source": res, "identifications": path, "columns": list(SAGE_IDS_COLUMNS),
           "kept_rule": SAGE_IDS_RULE,
           "methods_sentence": sage_ids_methods_sentence(groups, passing),
           "protein_filter": (f"none: `proteins` lists every database protein containing the "
                              f"peptide; Sage's picked-protein q-value is in `protein_q` "
                              f"({passing} of {groups} protein groups have protein_q <= "
                              f"{SAGE_IDS_Q:g})"),
           "q_value_scope": "Sage's own q-values, computed over all the runs together",
           "row_unit": "one row per (run, modified peptide, charge)",
           "quantified": False, "per_run": per_run, "totals": totals,
           "stale_report_moved_to": moved}
    try:
        rec["written_to"] = sage_lfq_check.write_record(out, SAGE_IDS_RECORD,
                                                        "sage_identifications", rec)
    except OSError as e:
        sys.stderr.write(f"[adapt] NOTE: could not write {SAGE_IDS_RECORD} in {out}: {e}\n")
    print(f"  [adapt] Sage identifications only -> {path}  ({totals['rows']} precursor x run rows "
          f"from {totals['psms_kept']} PSMs kept as {SAGE_IDS_RULE}: {totals['precursors']} "
          f"precursors, {totals['peptides']} peptides, {len(per_run)} runs; {groups} protein "
          f"groups, {passing} at protein_q <= {SAGE_IDS_Q:g}, none filtered on it; no "
          f"report.parquet)")
    return path


# --------------------------------------------------------------- FragPipe -----
# FragPipe reads a Thermo .raw through tools/ext/thermo/BatmassIoThermoServer.exe, a .NET
# Framework program: on Linux it needs mono (references/search-engines.md "Reading Thermo .raw";
# the .NET 8 runtime does not run it). HIVE has no mono module (`module avail mono` and
# `module search mono` find nothing, 2026-11-03); the Core keeps one at FRAGPIPE_MONO_SHARED
# (conda-forge mono 6.12.0.199), which the 38-run hydrolysate FragPipe search ran with.
FRAGPIPE_MONO_ENV = "FRAGPIPE_MONO_BIN"           # a folder holding `mono`; set = used or refused
FRAGPIPE_MONO_SHARED = "/quobyte/proteomics-grp/conda_envs/mono/bin"


def fragpipe_mono():
    """(the folder holding mono, where it was found) or (None, why there is none): $FRAGPIPE_MONO_BIN
    when set (it is then the only answer), else PATH, else the Core's shared copy."""
    env = os.environ.get(FRAGPIPE_MONO_ENV)
    if env is not None and env.strip():
        d = os.path.abspath(env.strip())
        if _is_exe(os.path.join(d, "mono")):
            return d, f"${FRAGPIPE_MONO_ENV}"
        return None, f"${FRAGPIPE_MONO_ENV}={d} holds no executable mono"
    hit = shutil.which("mono")
    if hit:
        return os.path.dirname(os.path.abspath(hit)), "PATH"
    if _is_exe(os.path.join(FRAGPIPE_MONO_SHARED, "mono")):
        return FRAGPIPE_MONO_SHARED, "the Core's shared mono (FRAGPIPE_MONO_SHARED)"
    return None, (f"no mono on PATH, ${FRAGPIPE_MONO_ENV} not set, and none at "
                  f"{FRAGPIPE_MONO_SHARED}")


def _is_exe(path):
    return os.path.isfile(path) and os.access(path, os.X_OK)


def refuse_fragpipe_raw_without_mono(bundle, files, sbatch):
    """Before anything is written: a FragPipe DDA search of Thermo .raw that will run on Linux (a
    SLURM job, or here on Linux) needs mono; -> {"bin", "source"} for the job's PATH, None when no
    input needs it, or exit (REFUSED) with the way out -- never a job that fails mid-search."""
    raws = [f for f in files if str(f).lower().endswith(".raw")]
    if not raws or str(bundle.get("acquisition", "DDA")).upper() == "DIA":
        return None
    if not sbatch and not sys.platform.startswith("linux"):
        return None                       # Windows / macOS FragPipe reads .raw itself
    where, why = fragpipe_mono()
    if where:
        return {"bin": where, "source": why}
    sys.exit(f"[run_search] REFUSED -- nothing was written or submitted. FragPipe reads Thermo "
             f".raw through tools/ext/thermo/BatmassIoThermoServer.exe, a .NET Framework program "
             f"that needs mono on Linux, and none was found ({why}); {len(raws)} input(s) are "
             f".raw. Either set {FRAGPIPE_MONO_ENV} to a folder holding mono (on HIVE, for the "
             f"Core group: {FRAGPIPE_MONO_SHARED}; outside the Core group, or on another "
             f"cluster: install conda-forge's mono, e.g. `micromamba create -n mono -c "
             f"conda-forge mono`, and set {FRAGPIPE_MONO_ENV}=<that env>/bin), or convert each "
             f".raw to mzML first and search those -- ThermoRawFileParser, as run_search.py's "
             f"Sage route converts them: ThermoRawFileParser -i=<run>.raw -b=<folder>/<run>.mzML "
             f"-f=2  -- then run_search.py --engine fragpipe ... --files <folder>/*.mzML.")


def refuse_fragpipe_dia_unstaged(bundle, files, out):
    """Before anything is written (review F18): FragPipe's DIA route stages only a timsTOF .d
    (diatracer_stage.py). On .raw or mzML it would hand MSFragger the raw paths, and MSFragger
    writes its uncalibrated mzML/mzBIN beside them -- the harm the DDA stage removes (N1) -- with
    no mono check either. Refused, with the ways out; staged properly in 2.15."""
    if str(bundle.get("acquisition", "")).upper() != "DIA":
        return
    other = [f for f in files if not str(f).rstrip("/").lower().endswith(".d")]
    if not other:
        return
    q = shlex.quote
    stage = os.path.abspath(out).rstrip(os.sep) + "_inputs"
    sys.exit(f"[run_search] REFUSED -- nothing was written or submitted. FragPipe's DIA route "
             f"stages only timsTOF .d (diatracer_stage.py); {len(other)} input(s) are .raw or "
             f"mzML, and MSFragger would write its uncalibrated mzML/mzBIN (~1 GB a run) beside "
             f"them, in the raw folder. Either search them with DIA-NN (--engine diann; it never "
             f"writes beside its inputs), or stage them yourself in a folder of your own and give "
             f"FragPipe those:\n  mkdir -p {q(stage)} && cp -p <each .raw> {q(stage)}/ && "
             f"ln -s <each .mzML> {q(stage)}/\n  run_search.py ... --engine fragpipe --files "
             f"{q(stage)}/*\n(a .raw needs mono on Linux: {FRAGPIPE_MONO_ENV}).")


# What MSFragger writes beside each staged input (<run>_uncalibrated.mzML / .mzBIN, <run>.mzBIN):
# ~1.1 GB a run, measured (a 38-run Orbitrap hydrolysate search's raw folder gained ~1.1 GB a run;
# H5: 2.2 GB for two timsTOF .d). The free-space check counts 1.2 GB a run on top of the copies.
STAGE_SIDE_FILES_BYTES_PER_RUN = int(1.2 * 1024 ** 3)


def _gb(n):
    return f"{n / 1024 ** 3:.1f} GB"


def _existing_parent(path):
    path = os.path.abspath(path)
    while not os.path.exists(path):
        up = os.path.dirname(path)
        if up == path:
            break
        path = up
    return path


def fragpipe_stage_problem(stage, out):
    """Why `stage` is not a safe FragPipe staging folder, or None (C2 security H2). The stage is
    emptied after a finished search, so it must be a REAL directory -- lstat, never a link -- at
    <realpath(--out)>/raw_stage, or not exist yet. A raw_stage linked to the raw folder (`ln -s`)
    would have every regular file in that folder deleted. The one test, used before anything is
    written, before the copies and again before the clean-up (stage_dir_is_safe_sh: the job's)."""
    import scratch_files
    want = os.path.join(os.path.realpath(out), scratch_files.FRAGPIPE_STAGE_DIR)
    try:
        st = os.lstat(stage)
    except FileNotFoundError:
        st = None
    except OSError as e:
        return f"{stage} could not be checked ({e})"
    if st is not None and stat.S_ISLNK(st.st_mode):
        return (f"{stage} is a link (to {os.path.realpath(stage)}), not a folder of its own: the "
                f"search empties its stage after it finishes, and through a link that would "
                f"delete the files the link points at")
    if st is not None and not stat.S_ISDIR(st.st_mode):
        return f"{stage} exists and is not a folder"
    if os.path.realpath(stage) != want:
        return (f"{stage} resolves to {os.path.realpath(stage)}, not to {want} inside the search "
                f"folder")
    return None


def refuse_fragpipe_stage_link(stage, out, now="nothing was written or submitted"):
    """sys.exit with fragpipe_stage_problem()'s reason and the way out, when there is one."""
    why = fragpipe_stage_problem(stage, out)
    if why:
        sys.exit(f"[run_search] REFUSED -- {now}. {why}. Move it out of the way (a link: `rm` "
                 f"on the link itself, never `rm -r` through it) and re-run: the search makes its "
                 f"own raw_stage folder and copies the .raw into it.")


def stage_dir_is_safe_sh(d):
    """The job's form of fragpipe_stage_problem() for the shell-quoted stage `d`: true when it is
    a real folder or does not exist yet, never a link (`[ -L ]` does not follow one)."""
    return f"{{ ! [ -L {d} ] && {{ ! [ -e {d} ] || [ -d {d} ]; }}; }}"


def refuse_fragpipe_stage(bundle, files, out):
    """Before anything is written, for a FragPipe DDA search (review F2, F10): refuse a raw_stage
    that is a link or not a folder of its own (fragpipe_stage_problem, C2 security H2); refuse an
    input given from inside <out>/raw_stage -- compared by path and by realpath, so an alias of
    the stage is caught -- the job empties that folder after a finished search; and refuse when
    the folder the stage goes in has not the room for the .raw copies plus MSFragger's files
    (STAGE_SIDE_FILES_BYTES_PER_RUN a run), with the numbers and the ways out.
    -> {"to_copy_bytes", "need_bytes", "free_bytes", "where"}, or None for a DIA search."""
    import scratch_files
    # FragPipe's manifest is one tab-separated line a run (both routes): a tab or a line break in
    # an input's path would split its row into a wrong one (C2 security n3)
    split = [repr(str(f)) for f in files if any(c in str(f) for c in "\t\n\r")]
    if split:
        sys.exit(f"[run_search] REFUSED -- nothing was written or submitted. {len(split)} input "
                 f"path(s) hold a tab or a line break, which would split FragPipe's manifest row "
                 f"(one tab-separated line a run): {', '.join(split[:3])}"
                 f"{' ...' if len(split) > 3 else ''}. Rename them, or link them under plain "
                 f"names, and give those.")
    if str(bundle.get("acquisition", "DDA")).upper() == "DIA":
        return None
    stage = os.path.join(os.path.abspath(out), scratch_files.FRAGPIPE_STAGE_DIR)
    refuse_fragpipe_stage_link(stage, out)
    real_stage = os.path.realpath(stage).rstrip(os.sep) + os.sep
    inside = [f for f in files
              if os.path.abspath(str(f)).startswith(stage.rstrip(os.sep) + os.sep)
              or os.path.realpath(str(f)).startswith(real_stage)]
    if inside:
        sys.exit(f"[run_search] REFUSED -- nothing was written or submitted. {len(inside)} input(s) "
                 f"lie inside {stage}, the search's own staging folder, which the job empties "
                 f"after a finished search: {', '.join(inside[:3])}{' ...' if len(inside) > 3 else ''}. "
                 f"Give the original raw files (the job copies them in itself).")
    to_copy = 0
    for f in files:
        src = os.path.abspath(str(f).rstrip("/"))
        if not src.lower().endswith(".raw"):
            continue
        dest = os.path.join(stage, os.path.basename(src))
        try:
            st = os.stat(src)
            d = os.lstat(dest)
            if not os.path.islink(dest) and d.st_size == st.st_size and d.st_mtime >= st.st_mtime:
                continue                               # already staged: the job reuses it
        except OSError:
            pass
        try:
            to_copy += os.path.getsize(src)
        except OSError:
            pass
    need = to_copy + len(files) * STAGE_SIDE_FILES_BYTES_PER_RUN
    where = _existing_parent(stage)
    free = shutil.disk_usage(where)[2]
    rec = {"to_copy_bytes": to_copy, "need_bytes": need, "free_bytes": free, "where": where}
    if free < need:
        sys.exit(f"[run_search] REFUSED -- nothing was written or submitted. The FragPipe stage in "
                 f"{stage} needs {_gb(need)} ({_gb(to_copy)} of .raw copies and ~1.2 GB a run of "
                 f"MSFragger's files for {len(files)} run(s)); {where} has {_gb(free)} free. Free "
                 f"space there, search into another --out on a filesystem with room, or convert "
                 f"the .raw to mzML first (an mzML is linked into the stage, not copied).")
    return rec


def fragpipe_dda_stage(files, out, manifest, per_run=True):
    """Stage a FragPipe DDA search's inputs inside the search folder (<out>/raw_stage,
    scratch_files.FRAGPIPE_STAGE_DIR) and write the manifest for the staged paths, one
    experiment per run (L11 N1, N4). MSFragger writes <run>_uncalibrated.mzML/.mzBIN and
    <run>.mzBIN beside each input (~1 GB a run), so the raw folder is never handed to it:
      * a Thermo .raw is COPIED -- in the job (golden rule 3: no heavy I/O here) -- because
        MSFragger cannot read a .raw through a link ("Instrument index not available for
        requested device");
      * a .d or an mzML is LINKED, as diatracer_stage.py links a .d (its link_into).
    `per_run` (the workflow quantifies): one experiment per run, named as the adapters and
    conditions.csv name it (sage_run_name) -- with the column empty, FragPipe refuses IonQuant
    MBR and writes no combined_protein.tsv (H5 gap 1). Otherwise (an identifications-only
    workflow) the column stays empty: ONE pooled experiment, so Philosopher filters all runs
    together, as FragPipe did before 2.14.2, and identifications.tsv takes each run from the
    PSMs' Spectrum (follow-up 3). -> the record (search_provenance.json `fragpipe_stage`)."""
    import diatracer_stage
    import scratch_files
    stage = os.path.join(os.path.abspath(out), scratch_files.FRAGPIPE_STAGE_DIR)
    # never a link: the inline route reaches here without refuse_fragpipe_stage (C2 security H2)
    refuse_fragpipe_stage_link(stage, out, now="nothing was staged or searched")
    os.makedirs(stage, exist_ok=True)
    refuse_fragpipe_stage_link(stage, out, now="nothing was staged or searched")
    rows, copied, linked, scratch = [], [], [], []
    for f in files:
        src = os.path.abspath(str(f).rstrip("/"))
        if src.lower().endswith(".raw"):
            dest = os.path.join(stage, os.path.basename(src))
            copied.append({"src": src, "dest": dest})
            scratch.append(dest)
        else:
            dest = diatracer_stage.link_into(stage, src)
            linked.append({"src": src, "link": dest})
        stem = os.path.splitext(dest)[0]
        scratch += [stem + "_uncalibrated.mzML", stem + "_uncalibrated.mzBIN", stem + ".mzBIN"]
        rows.append((dest, sage_run_name(src) if per_run else ""))
    write_file(manifest, "".join(f"{d}\t{e}\t\tDDA\n" for d, e in rows))
    return {"dir": stage, "copied": copied, "linked": linked,
            "experiments": ({sage_run_name(c["src"]): c["dest"] for c in copied + [
                {"src": x["src"], "dest": x["link"]} for x in linked]} if per_run else {}),
            "experiment_rule": ("one experiment per run: the workflow quantifies (IonQuant needs "
                                "experiments)" if per_run else
                                "none: one pooled experiment, because the workflow runs no "
                                "quantification (FragPipe's filter then pools all runs, as before "
                                "2.14.2; the run of each identification comes from its PSMs)"),
            # every regular file in the stage goes after a finished search (review F17): the
            # copies, MSFragger's and IonQuant's side files, the console log; links stay
            "removed_after_a_finished_search": ("every regular file in the stage (the .raw "
                                                "copies, MSFragger's and IonQuant's side files, "
                                                "the console log); links stay"),
            "msfragger_side_files": scratch,
            "why": ("MSFragger writes <run>_uncalibrated.mzML/.mzBIN beside each input, so the "
                    "inputs are staged here and the raw folder is never written to; a .raw is "
                    "copied (MSFragger cannot read one through a link), a .d / mzML linked. The "
                    "job removes the copies and MSFragger's files once FragPipe reports ALL JOBS "
                    "DONE; after a search that did not finish they are kept for a re-run. Scratch: "
                    "never zipped, catalogued or delivered (scratch_files.py)")}


FRAGPIPE_IDS_RECORD = "fragpipe_identifications.json"


def fragpipe_quantification(params):
    """What the workflow this FragPipe search runs quantifies (fragpipe_params.
    workflow_file_quantification, the one reading) -- {"quantified": True/False/None, ...}."""
    import fragpipe_params
    return fragpipe_params.workflow_file_quantification(params)


def fragpipe_not_quantified_record(quant):
    """The quantification record of a FragPipe search whose workflow runs none: identifications
    only, in sage_lfq_check.not_applicable()'s one shape -- the file run_de.R asks the gate about,
    so a DE is refused with this record's own words (Ruling 73a)."""
    import fragpipe_params
    return sage_lfq_check_mod().not_applicable(
        engine="fragpipe", why=quant["why"], ran="FragPipe ran no quantification",
        de_refusal=(f"REFUSED: this FragPipe search ran no quantification ({quant['why']}), so "
                    f"there is no report.parquet and no DE. Its identifications are in {SAGE_IDS} "
                    f"(run_search.py --adapt-only writes it). To quantify, search with a workflow "
                    f"that has {fragpipe_params.LFQ_KEY}=true (FragPipe's LFQ-MBR), or with Sage "
                    f"and quant.lfq true."))


def fragpipe_record_quantification(out, quant):
    """At generation: an identifications-only search writes its record now (the gate's, and
    search_provenance.json `quantification`), and an earlier quantifying search's report.parquet
    in `out` is set aside; a quantifying one sets aside an earlier identifications-only search's
    record and tables -- renamed .stale-<time>, never deleted (as sage_set_aside_other_kind).
    -> (the record or None, the paths set aside)."""
    slc = sage_lfq_check_mod()
    if quant.get("quantified") is False:
        names, rec = ["report.parquet"], fragpipe_not_quantified_record(quant)
    elif quant.get("quantified"):
        names = [SAGE_IDS, FRAGPIPE_IDS_RECORD]
        rec = {"status": "quantified", "why": quant["why"], "switches": quant.get("switches")}
        if (slc.load(out) or {}).get("status") == "not_applicable":
            names.append(slc.RECORD)
    else:
        return None, []
    moved = []
    for name in names:
        try:
            new = set_aside(os.path.join(out, name))
        except ValueError as e:
            sys.exit(f"[run_fragpipe] {e}. Move it yourself, or search into a different --out. "
                     f"Nothing was submitted.")
        if new:
            moved.append(new)
    if rec and rec["status"] == "not_applicable":
        try:
            rec["written_to"] = slc.write_record(out, slc.RECORD, "quantification", rec)
        except OSError as e:
            sys.stderr.write(f"[run_fragpipe] NOTE: could not write {slc.RECORD} in {out}: {e}\n")
    elif rec:
        _provenance_set(out, "quantification", rec)     # the gate's file is for IDs-only only
    if moved:
        print(f"  [run_fragpipe] an earlier search's outputs set aside: {', '.join(moved)}")
    return rec, moved


def _provenance_set(out, key, rec):
    """search_provenance.json[key] = rec in `out`, when that file exists (merged, never
    replaced); a write that fails is said."""
    prov = os.path.join(out, "search_provenance.json")
    try:
        with open(prov) as fh:
            p = json.load(fh)
    except (OSError, ValueError):
        return
    if isinstance(p, dict):
        p[key] = rec
        try:
            with open(prov + ".tmp", "w") as fh:
                json.dump(p, fh, indent=2)
            os.replace(prov + ".tmp", prov)
        except OSError as e:
            sys.stderr.write(f"[run_fragpipe] NOTE: could not record {key} in {prov}: {e}\n")


def fragpipe_run_quantification(out, params=None):
    """What a FragPipe output's search quantified: the workflow FragPipe saved in it
    (fragpipe.workflow), else `params`, as fragpipe_params reads it -- and when that does not
    say (no label-free switch; review F7), what the search wrote: combined_protein.tsv is
    quantified, identification tables without one are identifications only. None when neither
    tells."""
    import fragpipe_params
    wf = fragpipe_params.find(out, fragpipe_params.SAVED_WORKFLOW) or (
        params if params and str(params).endswith(".workflow") else None)
    quant = (fragpipe_params.workflow_file_quantification(wf) if wf else
             {"quantified": None, "label_free": False, "runs": [], "switches": {},
              "reports": {}, "why": "no workflow was found for this output"})
    if quant.get("quantified") is not None:
        return quant
    has_table = _written(out, [fragpipe_params.LFQ_TABLE]) == []
    intensities = _combined_protein_intensities(out) if has_table else None
    if has_table and intensities:
        return dict(quant, quantified=True, label_free=True,
                    why=f"the search wrote combined_protein.tsv ({quant['why']})")
    if _written(out, ["psm.tsv", "peptide.tsv"]) == []:
        return dict(quant, quantified=False,
                    why=("the search wrote identification tables and "
                         + ("a combined_protein.tsv with no intensities (spectral counts only)"
                            if has_table else "no combined_protein.tsv")
                         + f" ({quant['why']})"))
    return None


def _combined_protein_intensities(out):
    """Whether the combined_protein.tsv at the top of `out` holds any intensity > 0 in a
    '<sample> Intensity' / '<sample> MaxLFQ Intensity' column. With several experiments,
    Philosopher's abacus writes one for an identifications-only search too, with spectral counts
    and every intensity 0 (H7 A) -- not quantities (follow-up 3). False when unreadable."""
    import csv
    try:
        with open(os.path.join(out, "combined_protein.tsv"), newline="", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                for k, v in r.items():
                    if k and k.endswith("Intensity") and k != "Intensity":
                        try:
                            if float(v) > 0:
                                return True
                        except (TypeError, ValueError):
                            continue
    except OSError:
        return False
    return False


def adapt_fragpipe_ids(out, quant):
    """A FragPipe search whose workflow ran no quantification -> <out>/identifications.tsv (the
    name the Sage identifications-only route writes, SAGE_IDS) and <out>/fragpipe_identifications
    .json, and NO report.parquet: there are no quantities, so no DE (the gate refuses one).
    identifications.tsv is every experiment's peptide.tsv as Philosopher filtered it, one row per
    (run, peptide): `run` the run name (FragPipe's experiment name mapped back,
    fragpipe_experiment_names), `peptide`, `modifications` (Assigned Modifications), `charges`,
    `probability`, `proteins` (Protein ID, then Mapped Proteins). -> its path."""
    import csv
    import fragpipe_params
    back = fragpipe_experiment_names(out)
    tables = []                                           # each folder's psm.tsv
    for dp, dns, fns in os.walk(out):
        depth = os.path.relpath(dp, out).count(os.sep) + (dp != out)
        dns[:] = [d for d in dns if depth < 1 and not os.path.islink(os.path.join(dp, d))]
        if "psm.tsv" in fns:
            tables.append(os.path.join(dp, "psm.tsv"))
    if not tables:
        sys.exit(f"No psm.tsv under {out}: this FragPipe search ran no quantification "
                 f"({quant['why']}), and its identification tables are not there either. Check "
                 f"the FragPipe log -- a headless run can exit 0 on a crash.")
    # One row per (run, peptide, modifications): the run from each PSM's Spectrum (<run>.<scan>.
    # <scan>.<charge>, the input file's name), so it holds whether Philosopher filtered each
    # experiment or all runs together (the tables then at the workdir root, one set for all
    # runs); a peptide only when it is in the same folder's peptide.tsv (Philosopher's
    # peptide-level filter).
    best, psm_only = {}, []
    for t in sorted(tables):
        folder = os.path.dirname(t)
        exp = os.path.basename(folder) if folder != out else ""
        ok = None
        try:
            with open(os.path.join(folder, "peptide.tsv"), newline="", errors="replace") as fh:
                ok = {r.get("Peptide", "") for r in csv.DictReader(fh, delimiter="\t")}
        except OSError:
            ok = None                                     # said in the record: PSM level only
            psm_only.append(os.path.relpath(folder, out) if folder != out else ".")
        with open(t, newline="", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                pep = r.get("Peptide", "")
                if ok is not None and pep not in ok:
                    continue
                spec = r.get("Spectrum") or ""
                run = (spec.rsplit(".", 3)[0] if spec.count(".") >= 3 else back.get(exp, exp))
                key = (run, pep, r.get("Assigned Modifications", ""))
                mapped = (r.get("Mapped Proteins") or "").strip()
                row = best.setdefault(key, {
                    "run": run, "peptide": pep, "modifications": key[2], "charges": set(),
                    "probability": 0.0,
                    "proteins": ";".join(x for x in [r.get("Protein ID", "")] +
                                         [m.strip() for m in mapped.split(",")] if x)})
                row["charges"].add(r.get("Charge", ""))
                try:
                    row["probability"] = max(row["probability"], float(r.get("Probability")))
                except (TypeError, ValueError):
                    pass
    rows = sorted(({**v, "charges": ",".join(sorted(c for c in v["charges"] if c))}
                   for v in best.values()), key=lambda x: (x["run"], x["peptide"]))
    per_run = {}
    for x in rows:
        per_run[x["run"]] = per_run.get(x["run"], 0) + 1
    path = os.path.join(out, SAGE_IDS)
    cols = ["run", "peptide", "modifications", "charges", "probability", "proteins"]
    with open(path + ".part", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    os.replace(path + ".part", path)
    wf = fragpipe_params.find(out, fragpipe_params.SAVED_WORKFLOW)
    wf_filt = None
    if wf:
        try:
            wf_filt = fragpipe_params.read_workflow(wf)[0].get(
                fragpipe_params.PHILOSOPHER_FILTER_KEY)
        except OSError:
            wf_filt = None
    # the filter AS RUN, from FragPipe's log, else the workflow's line, said which: the one
    # reading the Methods use (fragpipe_params.filter_as_recorded; follow-up 3, H7)
    fr = fragpipe_params.filter_as_recorded(out, wf_filt)
    filt = fr["filter"]
    if not filt:
        filt_src = None
    elif fr["from_log"]:
        filt_src = "FragPipe's log, the filter as run"
    else:
        filt_src = ("the workflow's phi-report.filter line; FragPipe may add options when it "
                    "runs it, and its log line for it "
                    + (f"could not be read ({fr['log_error']})" if fr["log_error"]
                       else "was not found"))
    rec = {"source": sorted(tables), "identifications": path, "columns": cols,
           # what the inputs had, never more (follow-up 2, verify R2): the peptide-level filter
           # only where a folder had a peptide.tsv
           "kept_rule": ("the PSMs of psm.tsv, as FragPipe's Philosopher filter wrote it"
                         + (f" (filter {filt}, from {filt_src})" if filt
                            else " (its filter is not recorded here)")
                         + (", whose peptide is in the same folder's peptide.tsv"
                            if not psm_only else
                            f"; PSM level only for {', '.join(psm_only)} (no peptide.tsv there, "
                            f"so no peptide-level filter)"
                            if len(psm_only) == len(tables) else
                            f", whose peptide is in the same folder's peptide.tsv -- PSM level "
                            f"only for {', '.join(psm_only)} (no peptide.tsv there)")
                         + "; one row per run (from each PSM's Spectrum) and peptide, at its "
                           "best probability"),
           "peptide_filter_applied": (True if not psm_only else
                                      False if len(psm_only) == len(tables) else "partial"),
           "psm_level_only": psm_only,
           "quantification": quant, "per_run": per_run, "rows": len(rows),
           "filter": filt, "filter_source": filt_src,
           "methods_sentence": (fragpipe_params.philosopher_methods_sentence(filt) if filt
                                else None)}
    try:
        rec["written_to"] = sage_lfq_check_mod().write_record(out, FRAGPIPE_IDS_RECORD,
                                                              "fragpipe_identifications", rec)
    except OSError as e:
        sys.stderr.write(f"[adapt] NOTE: could not write {FRAGPIPE_IDS_RECORD} in {out}: {e}\n")
    fragpipe_record_quantification(out, quant)
    print(f"  [adapt] FragPipe identifications only -> {path} ({len(rows)} peptide rows, "
          f"{len(per_run)} run(s)); no report.parquet: {quant['why']}")
    return path


def fragpipe_stage_job(stage, command, log, expect=None, out=None):
    """The job's lines around FragPipe for a staged DDA search: copy each .raw in (into .part,
    then renamed; a copy already there of the same size and no older is reused), run FragPipe
    with its output also in `log`, then -- only when FragPipe exited 0, printed ALL JOBS DONE
    (headless FragPipe can exit 0 after a crash) AND wrote what its workflow makes (`expect`,
    fragpipe_params.expected_outputs: the identification tables, plus combined_protein.tsv when it
    ran label-free quantification -- an identifications-only workflow is a finished search
    without one, Ruling 73a) -- remove the staged copies and MSFragger's files; otherwise keep
    them, say so, and fail the job."""
    q = shlex.quote
    lines = ["# [raw_stage] the inputs are staged in the search folder, never written beside "
             "the raw files: a .raw is copied (MSFragger cannot read one through a link)",
             "_stage_copy() {",
             # a link left in the stage (a hand attempt) is not a copy: MSFragger cannot read a
             # .raw through one (review F10)
             '  if [ -L "$2" ]; then rm -f -- "$2"; fi',
             '  if [ -f "$2" ] && [ "$(wc -c < "$1")" -eq "$(wc -c < "$2")" ] && ! [ "$1" -nt "$2" ]; then',
             '    echo "[raw_stage] $2 already staged"; return 0',
             "  fi",
             '  cp -p -- "$1" "$2.part" && mv -f -- "$2.part" "$2"',
             "}",
             # a copy that fails (a full share): its partial file goes, and the job says what is
             # kept and how to free it (review F2)
             "_stage_fail() {",
             '  rm -f -- "$2.part"',
             f'  echo "[raw_stage] $1 could not be copied (cp exit $3): its partial copy was '
             f'removed; the copies already staged are KEPT for a re-run in" {q(stage["dir"])} '
             '"-- delete them to free the space" >&2',
             "  exit 1",
             "}"]
    d = q(stage["dir"])
    # a raw_stage made a link since generation would have the copies "found" in the raw folder,
    # MSFragger's files written there and every file in it removed after the search (C2
    # security H2): refused before anything is copied, and checked again before the clean-up
    lines += [f"if ! {stage_dir_is_safe_sh(d)}; then",
              f'  echo "[raw_stage] REFUSED: the stage" {d} "is a link or not a folder of its '
              f'own -- nothing was copied or searched. Remove the link itself (rm, never rm -r '
              f'through it) and resubmit." >&2',
              "  exit 1",
              "fi",
              f"mkdir -p {d}"]
    if stage["copied"]:
        # the room checked at generation, again here: time has passed (review F2)
        need_kb = (stage.get("space") or {}).get("need_bytes", 0) // 1024
        lines += [f"_free_kb=$(df -Pk {d} | awk 'NR==2 {{print $4}}')",
                  f'if [ -n "$_free_kb" ] && [ "$_free_kb" -lt {need_kb} ]; then',
                  f'  echo "[raw_stage] REFUSED: the stage needs {need_kb} KB, $_free_kb KB are '
                  f'free in" {d} "-- nothing was copied" >&2',
                  "  exit 1",
                  "fi"]
    lines += [f"_stage_copy {q(c['src'])} {q(c['dest'])} || _stage_fail {q(c['src'])} "
              f"{q(c['dest'])} $?" for c in stage["copied"]]
    if stage["copied"]:
        lines.append(f'echo "[raw_stage] {len(stage["copied"])} .raw staged in" {d}')
    lines += ["set +e", f"{command} 2>&1 | tee {q(log)}", "_fp_rc=${PIPESTATUS[0]}", "set -e",
              # what this workflow must have written (fragpipe_params.expected_outputs): an
              # identifications-only workflow is a finished search without combined_protein.tsv
              '_fp_missing=""']
    # `-print -quit` in a test, never `find | grep -q .`: under pipefail, grep -q closing the
    # pipe early gives find SIGPIPE (exit 141) once it prints more than a pipe block -- ~30 runs'
    # tables -- and a finished search was called unfinished (C2 errors_ux M3; as the mzML check)
    for t in (expect or ()):
        lines.append(f'[ -n "$(find {q(out)} -maxdepth 2 -name {q(t)} -size +0 -print -quit '
                     f'2>/dev/null)" ] || _fp_missing="$_fp_missing {t}"')
    lines += [f'if [ "$_fp_rc" -eq 0 ] && grep -q "ALL JOBS DONE" {q(log)} && '
              f'[ -z "$_fp_missing" ]; then',
              # every regular file in the stage -- copies, MSFragger's and IonQuant's side
              # files, the console copy (the job log has it all); links stay (review F17).
              # Only in a real folder: find -P never follows the stage or a link in it (C2 H2)
              f"  if {stage_dir_is_safe_sh(d)} && [ -d {d} ]; then",
              f"    find -P {d} -maxdepth 1 -type f -exec rm -f -- {{}} +",
              f'    echo "[raw_stage] the search finished: removed the staged copies and '
              f'MSFragger\'s uncalibrated mzML/mzBIN from" {d}',
              "  else",
              f'    echo "[raw_stage] WARNING: the search finished, but the stage" {d} "is no '
              f'longer a folder of its own (a link?): NOTHING in it was removed -- check it by '
              f'hand" >&2',
              "  fi",
              "else",
              '  _why="FragPipe exit $_fp_rc"',
              f'  grep -q "ALL JOBS DONE" {q(log)} || _why="$_why, no ALL JOBS DONE"',
              '  [ -z "$_fp_missing" ] || _why="$_why, not written:$_fp_missing"',
              f'  echo "[raw_stage] the search did not finish ($_why): the staged copies are KEPT '
              f'for a re-run in" {d} >&2',
              '  [ "$_fp_rc" -ne 0 ] || _fp_rc=1',
              '  exit "$_fp_rc"',
              "fi"]
    return "\n".join(lines)


def run_fragpipe(cmd, bundle, params, files, fasta, out, threads, sbatch, queue=None, mono=None,
                 space=None, database=None):
    """`mono`: refuse_fragpipe_raw_without_mono()'s {"bin", "source"} -- put on the job's PATH for
    a Thermo .raw input."""
    os.makedirs(out, exist_ok=True)
    acq = bundle.get("acquisition", "DDA").upper()
    manifest = os.path.join(out, "fragpipe.fp-manifest")
    dtype = "DIA" if acq == "DIA" else "DDA"
    dda_stage = None
    # What the workflow quantifies decides what a finished search has written and whether it is
    # identifications only (recorded now: the gate's record, search_provenance quantification)
    import fragpipe_params
    quant = fragpipe_quantification(params)
    quant_rec, _ = fragpipe_record_quantification(out, quant) if dtype == "DDA" else (None, [])

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
    elif dtype == "DDA":
        # never into the raw folder: the inputs go into <out>/raw_stage (fragpipe_dda_stage)
        # one experiment per run only when the workflow quantifies (follow-up 3)
        dda_stage = fragpipe_dda_stage(files, out, manifest,
                                       per_run=quant.get("quantified") is True)
        dda_stage["space"] = space                     # refuse_fragpipe_stage's numbers
        print(f"  [raw_stage] {len(dda_stage['copied'])} .raw to copy in the job, "
              f"{len(dda_stage['linked'])} linked -> {dda_stage['dir']}")
    else:
        write_file(manifest, "".join(f"{os.path.abspath(f)}\t\t\t{dtype}\n" for f in files))
    tools = os.environ.get("FRAGPIPE_TOOLS_FOLDER", "")
    tools_arg = f"--config-tools-folder {shlex.quote(tools)}" if tools else ""
    # The job's size reaches FragPipe: its Java heap (--ram, GB) inside the job's memory -- with
    # no --ram FragPipe sizes the heap from the NODE, not the job's cgroup (H5 gap 4) -- and its
    # threads from the CPUs SLURM gave the job (L11 N1)
    ram = (fragpipe_ram_gb(int(str((queue or {}).get("mem") or JOB_MEM_DEFAULT).rstrip("Gg")))
           if sbatch else None)
    size = (f' --ram {ram} --threads "${{SLURM_CPUS_PER_TASK:-{int(threads)}}}"' if sbatch
            else f" --threads {int(threads)}")
    # the workflow it RUNS: the given one with database.db-path = --fasta (H7 bug 1)
    database = database or refuse_fragpipe_database(params, fasta)
    run_wf = fragpipe_workflow_for_run(params, database, out)
    full = (f"{cmd} --headless --workflow {shlex.quote(run_wf)} "
            f"--manifest {shlex.quote(manifest)} --workdir {shlex.quote(out)} {tools_arg}{size}")
    stage = dda_stage
    said = {"fragpipe_stage": stage, "fragpipe_mono": mono, "fragpipe_database": database,
            "fragpipe_workflow": {"given": os.path.abspath(params), "run": run_wf}}
    if dtype == "DDA" and isinstance(quant.get("quantified"), bool):
        said["quantified"] = quant["quantified"]          # false = identifications only
        said["quantification"] = quant_rec or {"status": "quantified", "why": quant["why"],
                                               "switches": quant["switches"]}
    expect = fragpipe_params.expected_outputs(quant)
    # mono on PATH for a .raw (refuse_fragpipe_raw_without_mono found it), as DIA-NN's job gets
    # its .NET 8 (dotnet_env_for)
    pre = (f'export PATH={shlex.quote(mono["bin"])}:"$PATH"; '
           f'echo "[fragpipe] mono: $(command -v mono)"' if mono else "")
    log = os.path.join(stage["dir"], "fragpipe_console.log") if stage else None
    if sbatch:
        job = fragpipe_stage_job(stage, full, log, expect, out) if stage else full
        emit_sbatch(sbatch, job, out, threads, job="fragpipe_search", preamble=pre,
                    **(queue or {}))
        res = {"engine": "fragpipe", "out": out, "submitted": sbatch, "ran": False,
               "fragpipe_ram_gb": ram, **said}
        if quant.get("quantified") is False and dtype == "DDA":
            res["note"] = (f"Identifications only ({quant['why']}): after the job runs, re-run "
                           f"with --adapt-only to write {SAGE_IDS}. There is no report.parquet, "
                           f"so no differential expression.")
        return res
    if stage:
        # `pre` ends in an echo: joined with a space it would print the command instead of
        # running it (review F1), so it is its own command
        fragpipe_stage_inline(stage, f"{pre}; {full}" if pre else full, log, expect, out)
    else:
        sh(full)
    if quant.get("quantified") is False and dtype == "DDA":
        ids = adapt_fragpipe_ids(out, quant)
        return {"engine": "fragpipe", "report": None, "identifications": ids, "ran": True,
                **said}
    report = adapt_fragpipe(out)
    return {"engine": "fragpipe", "report": report, "ran": True, **said}


def _written(out, names):
    """The names in `names` that no non-empty file of that name holds at depth <= 2 of `out`."""
    found = set()
    for dp, dns, fns in os.walk(out):
        depth = os.path.relpath(dp, out).count(os.sep) + (dp != out)
        dns[:] = [d for d in dns if depth < 1]
        found |= {f for f in fns if f in names and os.path.getsize(os.path.join(dp, f)) > 0}
    return [n for n in names if n not in found]


def empty_fragpipe_stage(stage_dir, out):
    """Remove every regular file directly in the stage, never through a link (C2 security H2):
    the stage must pass fragpipe_stage_problem(), is then opened O_NOFOLLOW (a stage swapped for a
    link since the check fails to open), and each entry is unlinked relative to that open folder
    only when it is a regular file and not a link. -> the names removed, or a str: why NOTHING
    was removed."""
    why = fragpipe_stage_problem(stage_dir, out)
    if why:
        return why
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(stage_dir, flags)
    except OSError as e:
        return f"{stage_dir} could not be opened as a folder of its own ({e})"
    removed = []
    try:
        for e in os.scandir(fd):
            if e.is_file(follow_symlinks=False) and not e.is_symlink():
                os.unlink(e.name, dir_fd=fd)
                removed.append(e.name)
    finally:
        os.close(fd)
    return removed


def fragpipe_stage_inline(stage, command, log, expect=(), out=None):
    """fragpipe_stage_job() for a search run here (inside an allocation): the same copy, run,
    and remove-only-after-a-finished-search (ALL JOBS DONE, and `expect` written). A stage that is
    a link or not a folder of its own is refused before anything is copied, and never emptied
    (C2 security H2)."""
    stage_out = out or os.path.dirname(os.path.abspath(stage["dir"]))
    refuse_fragpipe_stage_link(stage["dir"], stage_out, now="nothing was copied or searched")
    for c in stage["copied"]:
        if os.path.islink(c["dest"]):
            os.remove(c["dest"])                       # a link is not a copy (review F10)
        st_src = os.stat(c["src"])
        try:
            st = os.stat(c["dest"])
            if st.st_size == st_src.st_size and st.st_mtime >= st_src.st_mtime:
                print(f"  [raw_stage] {c['dest']} already staged")
                continue
        except OSError:
            pass
        try:
            shutil.copy2(c["src"], c["dest"] + ".part")
        except OSError as e:
            try:
                os.remove(c["dest"] + ".part")
            except OSError:
                pass
            sys.exit(f"[raw_stage] {c['src']} could not be copied ({e}): its partial copy was "
                     f"removed; the copies already staged are KEPT for a re-run in "
                     f"{stage['dir']} -- delete them to free the space")
        os.replace(c["dest"] + ".part", c["dest"])
    try:
        sh_logged(command, log)
        with open(log, errors="replace") as fh:
            done = "ALL JOBS DONE" in fh.read()
    except subprocess.CalledProcessError:
        done = False
    missing = _written(out, list(expect)) if done and out else []
    if not done or missing:
        sys.exit(f"[raw_stage] the search did not finish"
                 + (f" (not written: {', '.join(missing)})" if missing else "")
                 + f": the staged copies in {stage['dir']} are KEPT for a re-run (FragPipe's "
                 f"output: {log})")
    # every regular file; links stay (F17); never through a link (C2 security H2)
    removed = empty_fragpipe_stage(stage["dir"], stage_out)
    if isinstance(removed, str):
        print(f"  [raw_stage] WARNING: the search finished, but NOTHING was removed from the "
              f"stage: {removed} -- check it by hand", file=sys.stderr)
        return
    print(f"  [raw_stage] the search finished: removed the staged copies and MSFragger's "
          f"uncalibrated mzML/mzBIN from {stage['dir']} ({len(removed)} file(s))")


def warn_fragpipe_dia_contam(out):
    """A FragPipe DIA output searched on a database with `contam_` entries (Philosopher's --contam)
    is WARNED about when it is adapted: its report carries those contaminants as bare accessions
    (FragPipe's DIA route drops the tag in library.tsv), so run_de.R's contaminant filter cannot
    see them and they are tested as sample proteins. run_search.py refuses such a database before
    a new FragPipe search (refuse_fragpipe_contam_database); an output made before 2.9, or
    outside the skill, only gets this. The database is the one the output's own fragpipe.workflow
    names. -> the message, or None."""
    wf = _find(out, ["fragpipe.workflow"])
    if not wf:
        return None
    hits = [(p, n) for p, n in fragpipe_contam_entries(fragpipe_databases(wf, "")) if n]
    if not hits:
        return None
    msg = ("[adapt] WARNING: this FragPipe DIA search used a database with FragPipe `contam_` "
           "contaminant entries (" + "; ".join(f"{p}: {n}" for p, n in hits) + "). Its DIA route "
           "drops that tag, so report.parquet carries them as plain accessions: run_de.R's "
           "contaminant filter cannot recognise them, and they will be tested as sample proteins "
           "(a keratin sample's exemption cannot see them either). Re-search on a database built "
           "by fetch_fasta.py (contaminants tagged Cont_), or name them from the database and "
           "drop them before the DE.")
    print(msg, file=sys.stderr)
    return msg


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
    tsv = None if pq_path else in_dia_out("report.tsv")
    if pq_path or tsv:
        warn_fragpipe_dia_contam(out)
    if pq_path:
        fragpipe_dia_runs_check(out, pq_path)
        print(f"  [adapt] FragPipe DIA: DIA-NN output already meets the contract -> {pq_path}")
        return pq_path
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
    fragpipe_dia_runs_check(out, report)
    print(f"  [adapt] FragPipe DIA: {os.path.basename(tsv)} -> {report}  ({n} rows)")
    return report


# The manifest data types FragPipe hands its DIA-NN step: CmdDiann.java (FragPipe 24.0) drops DDA,
# GPF-DIA, DIA-Lib and DDA+ rows -- on the diaTracer reuse form the pseudo-spectra mzML are DDA
# rows and the .d the DIA-Quant ones.
FRAGPIPE_DIANN_TYPES = ("DIA", "DIA-Quant")


def fragpipe_dia_runs_check(out, report):
    """Stop, naming them, when FragPipe's DIA-NN report leaves out runs of the manifest this
    search wrote -- the check the DIA-NN routes run in their jobs (check_report_runs.verify, the
    one rule). DIA-NN exits 0 when it skips a run, inside FragPipe too. It counts runs, so
    FragPipe's renaming of a Thermo .raw (it hands DIA-NN <name>_uncalibrated.mzML) does not
    matter. A FragPipe output without the skill's manifest is said, not checked."""
    import check_report_runs
    manifest = os.path.join(out, "fragpipe.fp-manifest")
    try:
        with open(manifest, encoding="utf-8") as fh:
            rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    except OSError:
        sys.stderr.write(f"[run_search] NOTE: no {manifest}, so the FragPipe DIA report's runs "
                         "were not checked against the inputs\n")
        return
    inputs = [r[0] for r in rows if len(r) > 3 and r[3].strip() in FRAGPIPE_DIANN_TYPES]
    if not inputs:
        return
    try:
        ok, msg = check_report_runs.verify(report, inputs)
    except Exception as e:                       # an unreadable report is a failed search
        ok, msg = False, f"FAILED: could not read {report}: {type(e).__name__}: {e}"
    if not ok:
        sys.exit(f"[run_search] FragPipe DIA: {msg}")
    print(f"  [adapt] FragPipe DIA: {msg}")


def adapt_fragpipe(out):
    """FragPipe -> DIA-NN-shaped report.parquet. Tries the DIA route first (diaTracer
    leaves DIA-NN output in dia-quant-output/), then the DDA route (IonQuant's
    combined_protein.tsv). The two produce different files, so which one exists is
    what tells us which route ran."""
    dia = adapt_fragpipe_dia(out)
    if dia:
        return dia
    return adapt_fragpipe_dda(out)


def _fp_manifest(path):
    """An .fp-manifest -> {file path: experiment}, {} when it cannot be read."""
    try:
        with open(path, errors="replace") as fh:
            rows = [ln.rstrip("\r\n").split("\t") for ln in fh if ln.strip()]
    except OSError:
        return {}
    return {r[0]: r[1] for r in rows if len(r) >= 2}


def fragpipe_experiment_names(out):
    """{experiment as FragPipe named it -> the run name our manifest gave that file}. FragPipe
    turns '-' into '_' in an experiment name (its own copy of the manifest, fragpipe-files.
    fp-manifest, in the workdir: an experiment named S-1 became S_1; H5 gap 8), so
    combined_protein.tsv's sample columns are not the run names conditions.csv uses. Matched by
    file path; {} when either manifest is missing (an output from before 2.14.2 or from outside
    the skill keeps FragPipe's names)."""
    ours = _fp_manifest(os.path.join(out, "fragpipe.fp-manifest"))
    theirs = _fp_manifest(os.path.join(out, "fragpipe-files.fp-manifest"))
    return {e: ours[p] for p, e in theirs.items() if e and ours.get(p) and ours[p] != e}


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
    # Protein.Group is the accession, as DIA-NN reports it: "Protein ID" (P12345), not
    # "Protein" (sp|P12345|ALBU_HUMAN), which compare_searches.py could never match to DIA-NN
    # (0 shared proteins). The entry name and gene go where DIA-NN puts them. The group is every
    # member, as DIA-NN writes one: the leading protein, then FragPipe's Indistinguishable
    # Proteins (protein_ids.fragpipe_group), so the contaminant rule sees a Cont_/contam_ entry
    # listed only there -- it reads every accession of a group, as on the DIA-NN path. The
    # leading protein stays first, which is what compare_searches.py compares on. Genes and
    # Protein.Names stay the leading protein's (FragPipe gives no per-member gene).
    names = {"Genes": "Gene", "Protein.Names": "Entry Name"}
    names = {k: v for k, v in names.items() if v in rows[0]}
    back_to_runs = fragpipe_experiment_names(os.path.dirname(os.path.abspath(cp)))
    runs, prots, ints = [], [], []
    ann = {k: [] for k in names}
    for r in rows:
        pg = fragpipe_group(r.get("Protein", ""), r.get("Protein ID", ""),
                            r.get("Indistinguishable Proteins", ""))
        if not pg:
            continue
        for c in lfq_cols:
            sample = c.replace(" MaxLFQ Intensity", "").replace(" Intensity", "")
            sample = back_to_runs.get(sample, sample)        # FragPipe's name -> the run's
            val = r.get(c, "")
            try:
                v = float(val)
            except ValueError:
                v = float("nan")
            runs.append(sample); prots.append(pg); ints.append(v if v > 0 else float("nan"))
            for k, src in names.items():
                ann[k].append((r.get(src) or "").strip())
    n = len(prots)
    tbl = pa.table({"Run": runs, "Protein.Group": prots, "PG.MaxLFQ": ints, **ann,
                    # already FDR-filtered upstream -> 0.0 placeholders
                    **{c: [0.0] * n for c in FDR_REQUIRED}})
    # what these numbers are, for the DE's methods (rule 1): IonQuant's, not DIA-NN's
    tbl = declare_quantity(tbl, fragpipe_quantity_declaration(
        out, maxlfq_columns=lfq_cols[0].endswith("MaxLFQ Intensity")))
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
    running there rather than what the group is. The cap is the QOS's own (user_limits)."""
    user = os.environ.get("USER", "")
    cap = user_limits("high", "genome-center-grp").get("cpu") or HIVE_USER_CPU_CAP
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
        return max(cap - used, 0)
    except Exception:
        return None


def _sacctmgr_path():
    """sacctmgr, or None. It is frequently absent from PATH in a non-login shell, so the usual
    install locations are looked in too."""
    found = shutil.which("sacctmgr")
    if found:
        return found
    for c in ("/usr/bin/sacctmgr", "/usr/local/bin/sacctmgr",
              "/cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin/sacctmgr"):
        if os.path.exists(c):
            return c
    return None


def _associations():
    """The current user's SLURM associations, [(account, partition, qos)]; [] when sacctmgr is
    missing or cannot answer (no SLURM here)."""
    assoc = []
    sacctmgr = _sacctmgr_path()
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
    return assoc


def _tres(text, key):
    """One TRES out of a sacctmgr TRES string: 'cpu=64,gres/gpu=0,mem=1T' -> 64 for "cpu", and
    for "mem" 1024 (GB; SLURM's unit is MB when none is given). None when it is not there."""
    m = re.search(rf"(?:^|,){re.escape(key)}=(\d+(?:\.\d+)?)([KMGTP]?)(?:,|$)", text or "")
    if not m:
        return None
    if key != "mem":
        return int(float(m.group(1)))
    mb = float(m.group(1)) * {"K": 1 / 1024, "": 1, "M": 1, "G": 1024, "T": 1024 ** 2,
                              "P": 1024 ** 3}[m.group(2)]
    return int(mb // 1024)


def user_limits(partition, account, qos=None):
    """The limits a job on this queue runs under: the PER-USER caps of its QOS (MaxTRESPU: every
    job of yours on it counts) and the PER-JOB cap (MaxTRESPerJob: a bigger job never starts).

    Read from SLURM, not assumed: `sacctmgr show qos` gives genome-center-grp-high-qos
    MaxTRESPU cpu=64,mem=1T, publicgrp-high-qos MaxTRESPerJob cpu=8,mem=128G and
    publicgrp-low-qos neither (HIVE, 2026-10-01). A queue with no QOS named is looked up in your
    associations (a facility job on `high` writes no --qos, and SLURM assigns the association's).
    When sacctmgr cannot be read, genome-center-grp/high falls back to HIVE_USER_CPU_CAP and any
    other queue is UNKNOWN -- said in `source`, never guessed.
    Returns {"qos", "cpu", "mem_gb", "per_job_cpu", "per_job_mem_gb", "source"}; a limit that
    does not exist is None."""
    lim = {"qos": qos, "cpu": None, "mem_gb": None, "per_job_cpu": None, "per_job_mem_gb": None,
           "source": None}
    if not qos:
        hits = {q for a, p, q in _associations() if a == account and p == partition and q}
        if len(hits) == 1 and "," not in next(iter(hits)):
            lim["qos"] = qos = next(iter(hits))
    sacctmgr = _sacctmgr_path()
    out = None
    if qos and sacctmgr:
        try:
            r = subprocess.run([sacctmgr, "-nP", "show", "qos", qos,
                                "format=name,MaxTRESPU,MaxTRESPerJob"],
                               capture_output=True, text=True, timeout=30)
            out = r.stdout if r.returncode == 0 else None
        except (OSError, subprocess.SubprocessError):
            out = None
    row = next((ln.split("|") for ln in (out or "").splitlines()
                if ln.split("|")[0].strip() == qos and len(ln.split("|")) >= 3), None)
    if row:
        lim.update(cpu=_tres(row[1], "cpu"), mem_gb=_tres(row[1], "mem"),
                   per_job_cpu=_tres(row[2], "cpu"), per_job_mem_gb=_tres(row[2], "mem"),
                   source=f"sacctmgr show qos {qos} (MaxTRESPU / MaxTRESPerJob)")
    elif (account, partition) == ("genome-center-grp", "high"):
        lim.update(cpu=HIVE_USER_CPU_CAP,
                   source=(f"HIVE_USER_CPU_CAP={HIVE_USER_CPU_CAP} (sacctmgr could not be read "
                           f"for {qos or 'this queue'}; docs/QUEUE_SWITCHING.md)"))
    else:
        lim["source"] = (f"unknown: sacctmgr could not be read for {qos or 'this queue'} "
                         f"({account}/{partition}), so no per-user cap is assumed")
    return lim


def fit_to_queue(cpus, limits):
    """`cpus` lowered to what one job on this queue can ever get: the per-job cap, and the
    per-user cap (a job asking for more than your whole cap waits for ever, QOSMaxCpuPerUser)."""
    caps = [c for c in (limits.get("cpu"), limits.get("per_job_cpu")) if c]
    return min([cpus] + caps)


def array_task_cpus(n_tasks, ceiling, *, limits, mem_per_task_gb, max_simultaneous,
                    pinned=None):
    """CPUs per task of an array step (the DIA-NN chain's steps 2 and 4), and how many tasks run
    at once. The ONE sizing rule; diann_parallel.py calls it and records what it returns
    (search_provenance.json `cpu_sizing`).

    Under a per-user CPU cap C, k tasks of t CPUs run at once with k = C // t, so throughput is
    C / (t x T(t)) for a per-file time T(t). t x T(t) -- the CPU-seconds one file costs -- does
    not fall as t rises (no program scales better than linearly), so fewer CPUs per task never
    lowers throughput while there are files waiting. 32 per task under a 64 cap ran 2 of 29
    files at a time (2026-09-25: the rest PENDING with QOSMaxCpuPerUserLimit). So, on a
    capped queue, t = C // (the tasks that could run at once), never below ARRAY_TASK_MIN_CPUS
    and never above `ceiling` (what the caller asked for: run_search.py's --threads). Few files
    get more CPUs each, as many as the cap allows them all at once. Without a cap (publicgrp/
    low) nothing competes for one, so t = `ceiling`. `pinned` (--threads-per-file) is honoured
    as given, lowered only to fit the queue -- and said when it costs concurrency."""
    cap, mem_cap = limits.get("cpu"), limits.get("mem_gb")
    want = max(1, min(n_tasks, max_simultaneous or n_tasks))

    def at_once(t):
        k = [want]
        if cap:
            k.append(cap // t)
        if mem_cap and mem_per_task_gb:
            k.append(mem_cap // mem_per_task_gb)
        return max(1, min(k))

    asked = pinned or ceiling
    if pinned:
        t = fit_to_queue(pinned, limits)
        why = (f"pinned with --threads-per-file {pinned}"
               + (f", lowered to {t} to fit the queue" if t < pinned else ""))
    elif cap:
        t = fit_to_queue(min(ceiling, max(cap // want, ARRAY_TASK_MIN_CPUS)), limits)
        why = (f"sized to the per-user cap: {cap} CPUs shared by up to {want} task(s), "
               f"at least {ARRAY_TASK_MIN_CPUS} and at most {ceiling} (--threads) each")
    else:
        t = fit_to_queue(ceiling, limits)
        why = f"no per-user CPU cap on this queue: {t} CPUs each, as asked (--threads)"
    k, k_asked = at_once(t), at_once(fit_to_queue(asked, limits))
    rec = {"cpus": t, "concurrent": k, "tasks": n_tasks, "requested": asked,
           "pinned": bool(pinned), "concurrent_at_requested": k_asked,
           "per_user_cpu_cap": cap, "per_user_mem_cap_gb": mem_cap,
           "per_job_cpu_cap": limits.get("per_job_cpu"), "mem_per_task_gb": mem_per_task_gb,
           "max_simultaneous": max_simultaneous, "qos": limits.get("qos"),
           "cap_source": limits.get("source"), "reason": why}
    rec["summary"] = (f"{t} CPUs per file, up to {k} of {n_tasks} files at once"
                      + (f" ({limits.get('qos') or 'this queue'} caps you at {cap} CPUs"
                         + (f" / {mem_cap} GB" if mem_cap else "") + ")" if cap else "")
                      + (f"; {asked} per file would run {k_asked} at once" if asked != t else ""))
    return rec


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
    # Failing to read the associations must NOT silently emit an empty queue: SLURM would
    # then use the cluster default partition, which on HIVE is `high` — precisely the
    # queue a non-facility account cannot use.
    assoc = _associations()

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
                  f"{need}); publicgrp/low has {idle} idle -> using low (preemptible: a job "
                  f"script gets --requeue)", file=sys.stderr)
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


def _qc_choice(a):
    """--qc -> True, --not-qc -> False, neither -> None (FRAN's stage decides from the names)."""
    return True if getattr(a, "qc", False) else False if getattr(a, "not_qc", False) else None


def _exclude_arg(a):
    """--exclude for the chain generators (diann/radiant_parallel.py): the checked node list,
    comma-joined (they check it again; it is their own option too), or None."""
    nodes = getattr(a, "exclude_nodes", None)
    return ",".join(nodes) if nodes else None


def _job_end_argv(a):
    """The job-end choices as flags for the chain generators (diann/radiant_parallel.py), through
    notify_slack.job_end_flags -- the one mapping the re-run commands use too."""
    import notify_slack
    return notify_slack.job_end_flags(slack=not getattr(a, "no_notify", False),
                                      fran=not getattr(a, "no_fran", False),
                                      fran_name=getattr(a, "fran_name", None), qc=_qc_choice(a))


def _job_end_plan(a, engine):
    """The job-end hook's plan for this search (notify_slack.job_end_plan), for the record."""
    try:
        import notify_slack
        return notify_slack.job_end_plan(slack=not getattr(a, "no_notify", False),
                                         fran=not getattr(a, "no_fran", False),
                                         fran_guarded=(engine == "diann"),
                                         fran_name=getattr(a, "fran_name", None),
                                         qc=_qc_choice(a))
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


# Memory of a single-job search's job(s) when --mem is not given; recorded "64G (default)".
JOB_MEM_DEFAULT = "64G"

# The job of a Sage or FragPipe search with a BROAD digest (non-specific, or semi-specific), when
# --mem / --time are not given (L11 N1/N2). MEASURED on HIVE (sacct, 2026-10-30 and 11-02), a
# 38-run Orbitrap DDA hydrolysate peptidomics search against a 6,300-sequence database, 32 CPUs:
#   * Sage 0.14.7, cleave_at "" (7-30 residues: 165 M peptides, 3.87 G fragments), .raw converted
#     in the job: 0:45:11 and 1:01:52 elapsed, MaxRSS 126.2 and 121.9 GiB (132,333,884 and
#     127,832,764 KiB) at --mem 192G. 192G is 1.5x the larger.
#   * FragPipe 24.0 Nonspecific-peptidome (7-65 residues): 1:55:15 elapsed (FragPipe: ALL JOBS
#     DONE IN 112.6 MINUTES, MSFragger 97.4 of them), MaxRSS 200.0 GiB at --mem 200G --ram 190:
#     the Java heap grows to what it is given, so this is the request, not a need below it; 200G
#     is what that search finished with.
#   * Time: both are far under 12 h; it scales with the runs (FragPipe 2.96 min a run here), so
#     24 h keeps >= 10x headroom at 38 runs and covers ~480 runs of FragPipe at that rate.
# A semi-specific digest gets the same (smaller than non-specific; not measured on its own).
# Memory follows the DATABASE (peptides x fragments), not the runs: a bigger database needs more
# -- the generation NOTE says what this was measured on.
BROAD_DIGEST_SAGE_MEM_GB = 192
BROAD_DIGEST_FRAGPIPE_MEM_GB = 200
BROAD_DIGEST_HOURS = 24
BROAD_DIGEST_MEASURED = ("measured on a 38-run hydrolysate search against a 6,300-sequence "
                         "database: Sage non-specific 7-30 residues peaked at 126.2 GiB (two "
                         "runs of it took 0:45 and 1:02 on 32 CPUs), FragPipe "
                         "Nonspecific-peptidome finished in 1:55 at 200G")
BROAD_DIGEST_SEQUENCES = 6300


def broad_digest_source(digest):
    """The words a digest-sized default is recorded with (rule 2: a default, said so)."""
    if digest == "non-specific":
        return "default for a non-specific digest (measured on a 38-run hydrolysate search)"
    return (f"default for a {digest} digest (the non-specific one's, measured on a 38-run "
            f"hydrolysate search)")


def digest_breadth(engine, params):
    """(how broad this Sage / FragPipe search's digest is -- "non-specific", "semi-specific" or
    None for an enzyme rule -- and where that was read), from its parameters file as
    cleavage_rules reads it. A file that cannot be read is None, said."""
    import cleavage_rules
    try:
        if engine == "sage":
            with open(params) as fh:
                enz = ((json.load(fh) or {}).get("database") or {}).get("enzyme")
            if not isinstance(enz, dict) or not enz:      # no block: Sage's trypsin default
                return None, f"{os.path.basename(params)} (no enzyme block: Sage's trypsin)"
            spec = cleavage_rules.sage_cleavage(enz, params)["specificity"]
        elif engine == "fragpipe":
            import fragpipe_params as fpp
            if os.path.basename(params) == fpp.FRAGGER_PARAMS:
                kv = fpp.read_fragger_params(params)
            else:
                props, _ = fpp.read_workflow(params)
                kv = {k[len("msfragger."):]: v for k, v in props.items()
                      if k.startswith("msfragger.")}
            spec = fpp.digest(kv)["specificity"]
        else:
            return None, None
    except (OSError, ValueError, AttributeError) as e:
        return None, f"{params} could not be read ({type(e).__name__}: {e})"
    return cleavage_rules.BROAD.get(spec), os.path.basename(params)


def single_job_sizing(engine, params, mem_gb=None, hours=None):
    """The memory and wall time of a single-job search's job(s), and where each came from: the
    user's --mem / --time, else a broad digest's default (Sage, FragPipe), else the route's own
    default (64G; 12 h -- the DIA-NN single shot adds its probe's hours, search_job_hours).
    -> {"digest", "digest_source", "mem" ("<GB>G" or None), "mem_source", "hours" (int or
    None), "hours_source"}. None = the route's own default."""
    digest, where = (digest_breadth(engine, params) if engine in ("sage", "fragpipe")
                     else (None, None))
    broad = broad_digest_source(digest) if digest else None
    default_mem = (BROAD_DIGEST_SAGE_MEM_GB if engine == "sage" else BROAD_DIGEST_FRAGPIPE_MEM_GB)
    mem = (f"{mem_gb}G", "user (--mem)") if mem_gb else \
        (f"{default_mem}G", broad) if broad else (None, "default")
    hrs = (int(hours), "user (--time)") if hours else \
        (BROAD_DIGEST_HOURS, broad) if broad else (None, "default")
    return {"digest": digest, "digest_source": where, "mem": mem[0], "mem_source": mem[1],
            "hours": hrs[0], "hours_source": hrs[1]}


def sizing_note(sizing):
    """The generation NOTE for a digest that could not be read (review F15): the default then
    applies without anyone being told why. "" when there is nothing to say."""
    src = (sizing or {}).get("digest_source") or ""
    if (sizing or {}).get("digest") is None and "could not be read" in src:
        return (f"[run_search] NOTE: the digest could not be read for the job's size ({src}); "
                f"the job gets the route's default memory and time -- a non-specific or "
                f"semi-specific search needs more (--mem GB / --time HOURS).")
    return ""


def refuse_job_over_queue_cap(mem, limits, mem_source):
    """Before anything is written (review F11): a single job asking for more memory than its
    QOS lets one job have (user_limits per_job_mem_gb: publicgrp-high-qos MaxTRESPerJob
    mem=128G on HIVE) would never start -- refused, with the options. A cap that is not known
    checks nothing. -> None when it fits."""
    cap = (limits or {}).get("per_job_mem_gb")
    try:
        gb = int(str(mem).rstrip("Gg"))
    except (TypeError, ValueError):
        return None
    if not cap or gb <= cap:
        return None
    sys.exit(f"[run_search] REFUSED -- nothing was written or submitted. The job would ask for "
             f"{gb}G ({mem_source}), but {limits.get('qos') or 'this queue'} lets one job have at "
             f"most {cap}G ({limits.get('source')}), so it would never start. Either search on "
             f"a queue without that cap (--partition/--account/--qos: genome-center-grp high for "
             f"the Core; publicgrp low has no per-job memory cap but is preemptible), or pass "
             f"--mem {cap} and accept that the search may run out of memory (a non-specific "
             f"search was {BROAD_DIGEST_MEASURED}).")


def job_size_record(sizing, emitted):
    """search_provenance.json job_mem / job_time (the last job script's header: the search job),
    each with its source: "128G" / "36h" for the user's flag, "64G (default)" / "12h (default)",
    "192G (default for a non-specific digest)". All None without a single-job script."""
    if not sizing or not emitted:
        return {"job_mem": None, "job_mem_source": None, "job_time": None,
                "job_time_source": None, "job_digest": None, "job_digest_source": None}
    last = emitted[-1]

    def tag(source):
        return ("" if source.startswith("user") else " (default)" if source == "default"
                else f" (default for a {sizing['digest']} digest)")
    return {"job_mem": f"{last['mem']}{tag(sizing['mem_source'])}",
            "job_mem_source": sizing["mem_source"],
            "job_time": f"{last['hours']}h{tag(sizing['hours_source'])}",
            "job_time_source": sizing["hours_source"],
            # what the size was decided from: the digest read, and where (review F15)
            "job_digest": sizing.get("digest"), "job_digest_source": sizing.get("digest_source")}


def _fasta_sequences(fasta):
    """The database's sequence count as fetch_fasta.py's sidecar records it (n_sequences), or None
    -- the FASTA itself is not read here (a login node)."""
    try:
        with open(str(fasta) + ".meta.json") as fh:
            n = json.load(fh).get("n_sequences")
        return int(n) if isinstance(n, int) and not isinstance(n, bool) else None
    except (OSError, ValueError, AttributeError):
        return None


def fragpipe_ram_gb(mem_gb):
    """FragPipe's --ram (its Java heap, GB) for a job of `mem_gb`: what is left after a margin for
    the rest of the job (Python helpers, Percolator, the JVM's own memory) -- 10 GB or 5%, the
    margin the 200G / --ram 190 hydrolysate search finished with -- never under half the job."""
    margin = max(10, -(-mem_gb * 5 // 100))
    return max(mem_gb - margin, mem_gb // 2, 1)


# What emit_sbatch() wrote in this run: [{"path", "mem", "hours"}] -- search_provenance.json
# job_mem / job_time are read from it, so the record says what the job scripts say.
EMITTED_JOBS = []


def emit_sbatch(path, command, out, threads, job, preamble="",
                partition=None, account=None, qos=None, mem=JOB_MEM_DEFAULT,
                hours=SEARCH_WALL_HOURS,
                submit_hint=True, notify="final", stage="search", fran_guarded=False,
                preflight=None, exclude=None, time=None):
    """Emit a minimal SLURM script (login-node-safe). Orchestrator submits it.
    The queue is DETECTED from the submitting user's own SLURM associations — see
    slurm_queue() — unless the caller passes one, which then wins. Every caller must forward
    the command line's --partition/--account/--qos: on the 2026-09-16 FRAN pilot run_diann()
    did not, so `--partition low --account publicgrp --qos publicgrp-low-qos` came out as
    genome-center-grp/high in both job headers and had to be hand-edited before submission.
    `exclude`: the nodes --exclude keeps the search off (node_fault.exclude_nodes), carried in
    the same queue dict. `preamble` runs before the command (e.g. the DOTNET_ROOT exports that
    let DIA-NN 2.6 read Thermo .raw).

    `notify` is the job-end hook (notify_slack.wrap_job_script: run log -> FRAN -> Slack):
    "final" for the job that ends the search (success and failure; FRAN on success), "fail" for
    a job another waits on with afterok (a failure only -- the search stops there), None for no
    hook. --no-notify / SKILL_SLACK=0 at generation switch off only the Slack post; --no-fran /
    FRAN_DEPOSIT=off at generation, the FRAN hand-over. `fran_guarded`: this job ends with a
    completeness guard (report_guard), so a finished search may be staged for FRAN from the job;
    without one it is left to step 7c. `preflight`: node_fault.preflight_lines(), run first in the
    job (through the hook). `time`: run_search.py --time HOURS (carried in the queue dict, as
    --mem is), which wins over the route's `hours`."""
    if time:
        hours = int(time)
    EMITTED_JOBS.append({"path": path, "mem": mem, "hours": hours})
    part, acct, q = slurm_queue(partition, account, qos)
    pre = (preamble + "\n") if preamble else ""
    lines = [
        "#!/bin/bash -l",
        f"#SBATCH --job-name={job}",
        # in double quotes: a session folder with a space made sbatch refuse the whole script
        # ("Invalid directive found in batch script"); quoted, HIVE's SLURM accepts it and writes
        # the log there (job 24808508, 2026-10-15). --out cannot hold a '"' (refuse_unsafe_path).
        f'#SBATCH --output="{os.path.join(out, job)}_%j.log"',
        f"#SBATCH --cpus-per-task={threads}",
        f"#SBATCH --mem={mem}",
        f"#SBATCH --time={hours}:00:00",
    ]
    if part:  lines.append(f"#SBATCH --partition={part}")
    if acct:  lines.append(f"#SBATCH --account={acct}")
    if q:     lines.append(f"#SBATCH --qos={q}")
    if exclude:
        import node_fault
        lines.append(node_fault.exclude_directive(exclude))
    # Preemptible queue: requeue a preempted search rather than lose it (HIVE's JobRequeue=1
    # does this by default; another cluster's may not). The rule is
    # diann_parallel.needs_requeue(), shared with the chain's step headers.
    requeue = _diann_parallel_mod().needs_requeue(part, q)
    if requeue:
        lines.append("#SBATCH --requeue")
    lines += ["set -euo pipefail", f"cd {shlex.quote(os.path.abspath(out))}", f"{pre}{command}", ""]
    script = "\n".join(lines)
    if notify:
        import notify_slack
        script = notify_slack.wrap_job_script(
            script, os.path.abspath(out), final=(notify == "final"), time_limit_h=hours,
            stage=stage, slack=not getattr(a_globals, "no_notify", False),
            fran=not getattr(a_globals, "no_fran", False), fran_guarded=fran_guarded,
            fran_name=getattr(a_globals, "fran_name", None), qc=_qc_choice(a_globals),
            preflight=preflight)
    write_file(path, script, 0o755)
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
    ap.add_argument("--threads", type=int, default=8,
                    help="CPUs for the search job(s); on the 5-step chain, the MOST an array task "
                         "gets -- each is sized to the queue's per-user cap (--threads-per-file)")
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
    ap.add_argument("--no-notify", action="store_true",
                    help="no Slack post from the generated job(s) when the search ends or "
                         "fails (same as SKILL_SLACK=0). The jobs still log the run and hand a "
                         "finished Core search to FRAN; see references/notifications.md")
    ap.add_argument("--no-fran", action="store_true",
                    help="never hand this search to FRAN from its jobs (same as FRAN_DEPOSIT=off "
                         "when generating; baked into the job, because a hive_remote shell's "
                         "environment does not reach it)")
    ap.add_argument("--fran-name",
                    help="the session's descriptive name: FRAN's corpus name for this search, "
                         "handed to the job's `fran_deposit.py stage --name`")
    qc = ap.add_mutually_exclusive_group()
    qc.add_argument("--qc", action="store_true",
                    help="an instrument QC / standard run: the job never stages it for FRAN")
    qc.add_argument("--not-qc", action="store_true",
                    help="not a QC run, whatever its name (passed on to stage as --not-qc)")
    ap.add_argument("--parallel-threshold", type=int, default=5,
                    help="DIA-NN: use the 5-step SLURM chain above this many files (default 5)")
    ap.add_argument("--no-parallel", action="store_true",
                    help="force a single-shot DIA-NN search regardless of file count")
    ap.add_argument("--partition", help="SLURM partition for every job this writes (the "
                    "single-shot --sbatch job(s) and the parallel chain); detected if omitted")
    ap.add_argument("--account", help="SLURM account for every job this writes")
    ap.add_argument("--qos", help="SLURM QOS for every job this writes")
    ap.add_argument("--exclude", metavar="NODE[,NODE...]",
                    help="SLURM nodes to keep every job of this search off, e.g. "
                         "hive-dc-7-4-50,hive-dc-7-4-[52-54]: written as #SBATCH --exclude into "
                         "every job script (both single-shot jobs, every step of the chain) and "
                         "recorded in search_provenance.json (queue.exclude). sbatch reads no "
                         "SBATCH_EXCLUDE variable, so exporting one does nothing")
    ap.add_argument("--mem", type=_positive_int, metavar="GB",
                    help="memory in GB for the job(s) of a single-job search -- Sage, DIA-NN "
                         "single shot (both of its jobs), FragPipe, AlphaDIA, Radiant -- written "
                         f"as #SBATCH --mem=<GB>G (default: {JOB_MEM_DEFAULT}); recorded as "
                         "search_provenance.json job_mem. The 5-step chain sizes its own jobs: "
                         "--mem-per-file, --assembly-mem. A non-specific or semi-specific Sage "
                         f"or FragPipe search defaults to {BROAD_DIGEST_SAGE_MEM_GB}G (Sage) / "
                         f"{BROAD_DIGEST_FRAGPIPE_MEM_GB}G (FragPipe), recorded as that default")
    ap.add_argument("--time", type=_positive_int, metavar="HOURS",
                    help="wall-clock hours for the job(s) of a single-job search -- the same "
                         "routes as --mem -- written as #SBATCH --time=<HOURS>:00:00 (default: "
                         f"{SEARCH_WALL_HOURS}; the DIA-NN single shot adds its mass-accuracy "
                         f"probe's hours; a non-specific or semi-specific Sage or FragPipe search "
                         f"{BROAD_DIGEST_HOURS}); recorded as search_provenance.json job_time. The "
                         "5-step chain sizes its own jobs: --time-per-file")
    ap.add_argument("--max-simultaneous", type=int,
                    help="cap concurrent array tasks in the parallel chain")
    ap.add_argument("--threads-per-file", type=_positive_int,
                    help="parallel chain: PIN the CPUs of each array task (steps 2 and 4) and "
                         "step 1b. Default: sized to the queue's per-user CPU cap so more files "
                         "run at once -- on genome-center-grp/high (64 CPUs per user) 8 each, "
                         "8 files at a time, never more than --threads (array_task_cpus; "
                         "recorded as cpu_sizing in search_provenance.json)")
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
                         "(diann_parallel.py default: 4 h at 16 CPUs per task -- real runs' "
                         "p99 is 160 min, p99.5 248 min -- scaled up for fewer: 8 h at 8; "
                         "array_task_hours)")
    ap.add_argument("--mem-per-file", type=_positive_int,
                    help="parallel chain: memory in GB per array task in steps 2 and 4, and step "
                         "1b (diann_parallel.py default: 64 GB for steps 1b and 2, 48 GB for step "
                         "4 -- from real Core runs' sacct memory); raise it for a much larger "
                         "library")
    ap.add_argument("--adapt-only", action="store_true",
                    help="skip the search; just build report.parquet from an existing engine output dir")
    ap.add_argument(PWIZ_SIF_FLAG, metavar="SIF",
                    help="Sage, quantifying a timsTOF .d: the vendor-enabled ProteoWizard "
                         "container the job converts it to mzML with (default: the Core's "
                         "pinned Skyline image on HIVE, resolve_defaults.SKYLINE_PIN, whose "
                         "sha256 the job checks first). Not used for identifications only "
                         "(quant.lfq false), where Sage reads the .d itself")
    ap.add_argument("--accept-lfq-window", metavar="REASON",
                    help="Sage: build report.parquet although the LFQ window did not fit the "
                         "runs' MS1 mass error (sage_lfq_check.py's gate), because the user "
                         "accepts those quantities as they are -- give WHY (no name: who is the "
                         "login that runs this, named to clients by role); recorded in "
                         "sage_lfq_check.json and search_provenance.json, and AUDIT.md says so")
    ap.add_argument("--fasta-meta", metavar="SIDECAR",
                    help="fetch_fasta.py's sidecar describing --fasta's targets, when --fasta has "
                         "none of its own (a FragPipe search's Philosopher target+decoy .fas: its "
                         "target FASTA's <name>.meta.json is found beside it or in the session's "
                         "input/ when it matches; else give it here). 'none': the database has "
                         "no sidecar (not built by fetch_fasta.py), recorded so")
    ap.add_argument("--keratin-sample", action="store_true",
                    help="the samples ARE keratin (hair, wool, feather, skin, nail ...): refuse a "
                         "--fasta that still holds keratin-family Cont_ entries (build it with "
                         "fetch_fasta.py fetch --keratin-sample). Also implied by a sidecar built "
                         "that way; recorded in search_provenance.json as keratin_sample")
    ap.add_argument(_diann_parallel_mod().RENAME_DASHED_FLAG, action="store_true",
                    help="DIA-NN: search an input whose FILE NAME holds '--' through a link "
                         "named with each run of dashes made one, so its run is renamed (DIA-NN "
                         "cannot be given a path with '--'); without it such an input is refused. "
                         "A '--' in a folder only is linked under the same file name either way")
    ap.add_argument("--acquisition-json", metavar="JSON",
                    help="step 2's detect_acquisition.py output (input/acquisition.json): a .d "
                         "whose analysis.tdf check it records is not read again while its files "
                         "are unchanged. Default: the first input/acquisition.json up to 3 "
                         "folders above --out")
    ap.add_argument(ALLOW_DAMAGED_TDF, action="store_true",
                    help="search a Bruker .d whose analysis.tdf is not `ok` (truncated, "
                         "stale -wal/-journal beside it, WAL-mode header, unreadable) "
                         "anyway; every one is printed first")
    a = ap.parse_args()

    global a_globals
    a_globals = a
    # Nodes to keep every job off: checked before anything is written, because the value lands
    # on a #SBATCH line of every job (node_fault.exclude_nodes, the one rule for it).
    exclude = None
    if a.exclude is not None:
        import node_fault
        try:
            exclude = node_fault.exclude_nodes(a.exclude)
        except ValueError as e:
            sys.exit(f"[run_search] REFUSING --exclude {a.exclude!r}: {e}. Nothing was written.")
    a.exclude_nodes = exclude                  # for the chain generators (_exclude_arg)
    # The queue the user asked for, handed to every emit_sbatch() -- see emit_sbatch() -- and
    # recorded as search_provenance.json `queue`; `exclude` is in every job script's header.
    queue = {"partition": a.partition, "account": a.account, "qos": a.qos, "exclude": exclude}
    # --mem: every single-job route's emit_sbatch() takes it from here; without it, its default
    if a.mem:
        queue["mem"] = f"{a.mem}G"
    # --time, the same way (emit_sbatch's `time`); the queue holds only what the user asked for --
    # a digest-sized default goes to the route alone (route_queue, below), so a re-run built from
    # the record (sage_lfq_check.rerun_plan) recomputes it instead of replaying it as a choice
    if a.time:
        queue["time"] = a.time
    del EMITTED_JOBS[:]
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

    # Repeats in the input list, for every engine (check_report_runs: the one rule and the one
    # wording). One file listed more than once is searched ONCE and flagged -- on stderr and in
    # search_provenance.json `repeated_paths`, which the analysis report shows. Different files
    # that share a run name are flagged too (`repeated_names`) -- and stop here, before any
    # engine: DIA-NN names a run by its file name without the folder, so /plate1/s1.raw and
    # /plate2/s1.raw would become ONE Run (two samples merged; in the chain two array tasks
    # writing one .quant), and Sage converts both to one mzML name. Knowable from the list, so
    # it stops now rather than in check_report_runs.py after the search has run.
    import check_report_runs
    files, repeated_paths = check_report_runs.distinct_inputs(files)
    repeated_names = check_report_runs.repeated_names(files)
    if repeated_paths or repeated_names:
        sys.stderr.write(check_report_runs.repeats_note(repeated_paths, "run_search", repeated_names))
    if repeated_names and not a.adapt_only:
        sys.exit(check_report_runs.names_stop(repeated_names, "run_search"))

    if a.adapt_only:
        if engine == "sage" and not sage_lfq_enabled(a.params):
            # identifications only (quant.lfq false): no LFQ check and no report.parquet. The
            # record only after adapt_sage_ids: it refuses a quantified output, whose LFQ record
            # (a gate run_de.R reads) must then stay as it is.
            ids = adapt_sage_ids(a.out)
            sage_ids_lfq_record(a.out)
            print(json.dumps({"engine": "sage", "report": None, "identifications": ids,
                              "quantified": False, "ran": False, "adapt_only": True}, indent=2))
            return
        fq = fragpipe_run_quantification(a.out, a.params) if engine == "fragpipe" else None
        if fq and fq.get("quantified") is False:
            # a workflow that ran no quantification: its identifications, no report.parquet
            ids = adapt_fragpipe_ids(a.out, fq)
            print(json.dumps({"engine": "fragpipe", "report": None, "identifications": ids,
                              "quantified": False, "ran": False, "adapt_only": True}, indent=2))
            return
        if fq and fq.get("quantified"):
            # the record follows what ran: an IDs-only record from generation is set aside
            fragpipe_record_quantification(a.out, fq)
        # After a Sage job: the LFQ mass-window check again, here, where the orchestrator sees
        # it (the job ran it too, into its log) -- and its record into search_provenance.json.
        lfq = sage_lfq(a.out, a.params, accept=a.accept_lfq_window) if engine == "sage" else None
        refused = sage_lfq_check_mod().refusal(a.out, lfq) if lfq else None
        if refused:
            sys.exit(f"[run_search] {refused}")
        report = {"sage": adapt_sage, "fragpipe": adapt_fragpipe,
                  "alphadia": adapt_alphadia,
                  "radiant": adapt_radiant}.get(engine, lambda o: None)(a.out)
        if lfq:
            lfq = sage_lfq_check_mod().load(a.out) or lfq    # with the quantified coverage (M2)
        print(json.dumps({"engine": engine, "report": report, "ran": False, "adapt_only": True,
                          **({"sage_lfq_check": lfq} if lfq else {})}, indent=2))
        return

    if engine == "sage":
        refuse_sage_fragment_mismatch(bundle, a.params)
        refuse_sage_dotd_route(bundle, a.params, files)
    if engine == "fragpipe":
        refuse_fragpipe_contam_database(a.params, a.fasta)
        refuse_fragpipe_dashed(a.params, files, a.fasta, a.out, bundle.get("acquisition", ""))
        # Thermo .raw needs mono where the search runs (Linux): found now, or refused now
        refuse_fragpipe_dia_unstaged(bundle, files, a.out)
        # the database it searches: --fasta, never one FragPipe cached from an earlier search
        a.fragpipe_database = refuse_fragpipe_database(a.params, a.fasta)
        a.fragpipe_mono = refuse_fragpipe_raw_without_mono(bundle, files, a.sbatch)
        # the stage: no input from inside it, and room for it -- before anything is written
        a.fragpipe_stage_space = refuse_fragpipe_stage(bundle, files, a.out)
    # DIA-NN cuts its command line at every '--', quoted or not (diann_parallel.DIANN_OPTION_MARK):
    # an output, database or cfg path holding one is refused here, and an input run is searched
    # through a link in <out>/input_links/ -- made just before the search is generated, once
    # nothing is left to refuse -- or refused when the link would have to rename the run.
    links_plan = {"links": [], "refused": None}
    if engine == "diann":
        dp = _diann_parallel_mod()
        for flag, path in (("--out", a.out), ("--fasta", a.fasta), ("--params", a.params),
                           ("--library", a.library)):
            if path:
                dp.refuse_dashed_path(path, flag, prog="run_search")
        if os.path.isfile(a.params):
            dp.refuse_dashed_cfg(a.params, prog="run_search")
        links_plan = dp.dashed_input_plan(files, a.out, rename=a.rename_dashed_runs)
        if links_plan["refused"]:
            sys.exit(f"[run_search] REFUSING: {links_plan['refused']} Nothing was written or "
                     "submitted.")

    # Before anything is provisioned, submitted or run: no damaged .d gets searched here.
    # --adapt-only is above this on purpose -- it reads an engine's finished output and
    # never touches a raw file. Step 2's check is reused where it still describes the file;
    # many .d left to read on a login node move this command to a compute node (check_tdf).
    check_tdf(files, a.out, a.acquisition_json, a.allow_damaged_tdf, queue=queue)
    # Same place, same reason: a keratin sample is never searched on a database whose keratin-
    # family Cont_ entries would take its analyte out of quantification.
    fasta_meta_rec = fasta_sidecar_or_refuse(a, engine)
    keratin_rec = keratin_sample_check(a.fasta, a.keratin_sample, engine,
                                       sidecar=fasta_meta_rec["path"])

    cmd = tools.get(engine)
    if not cmd:
        sys.exit(f"tools.json has no command for engine '{engine}'. "
                 f"Re-run acquire_tools.sh, or check its notes:\n  "
                 + "\n  ".join(tools.get("notes", [])))
    # Before anything is submitted, so a version the user did not confirm is on screen
    # while it can still be stopped -- not only in a file read after the search.
    ver_rec = engine_version_record(engine, tools, bundle)
    global ENGINE_VERSION
    ENGINE_VERSION = ver_rec.get("value")     # for an adapter's quantity declaration

    # A DIA-NN cfg that is not there is reported as exactly that, first -- not as whatever a
    # later reader makes of an empty flag list ("mass accuracy is not pinned").
    if engine == "diann" and not os.path.isfile(a.params):
        sys.exit(f"[run_search] cfg not found: {a.params}")

    # Do this BEFORE parallel_decision: that reads the cfg to check mass accuracy, and the
    # augmented copy is the cfg the run will actually use.
    if engine == "diann":
        a.params = ensure_xic(a.params, a.out)
        ensure_temp_dirs(a.params, a.out)
    if engine in ("diann", "sage"):
        warn_contaminant_digest(engine, a.params, a.fasta, sidecar=fasta_meta_rec["path"])
    use_parallel, why = parallel_decision(engine, files, a.params, a)
    if engine == "diann":
        print(f"[run_search] parallel routing: {'YES' if use_parallel else 'no'} -- {why}")

    # HARD STOP: never start a multi-hour search inline on a cluster LOGIN NODE.
    # Golden rule #3 says every heavy step goes through the scheduler, but nothing
    # enforced it -- and the failure is silent and expensive. It bit for real: a
    # DIA-NN run whose mass accuracy was unpinned fell out of the parallel chain into
    # the single-shot path and, with no --sbatch, launched diann-linux on the HIVE
    # login node at 16 threads (twice, because an earlier invocation was still alive).
    # A login node is on_compute.on_login_node(), THE definition: no SLURM job here, and SLURM
    # (on PATH or at the fixed path) or HPCCF's cluster file is -- "sbatch on PATH" alone missed
    # a non-login shell, which has none (I44).
    import on_compute
    # Radiant on SLURM with >1 file: a per-file array + one Fulcrum job (radiant_parallel.py).
    # Decided HERE: like the DIA-NN chain it only writes job scripts, so the refusal below must
    # know it -- decided after it, a multi-file Radiant search was refused inline (L8R-F10).
    radiant_chain = (engine == "radiant" and len(files) > 1 and slurm_available()
                     and not a.no_parallel)
    # ...nor in the short, small srun check_tdf() moved this command to for the tdf check
    moved = bool(os.environ.get(TDF_REEXEC_ENV))
    if not (use_parallel or radiant_chain) and not a.sbatch \
            and (on_compute.on_login_node() or moved) and not a.allow_inline:
        where = ("this command was moved to a compute node only to read the .d files' "
                 "analysis.tdf, in a short, small allocation (check_tdf)" if moved else
                 "this is a cluster login/submit node (on_compute.on_login_node: no SLURM job "
                 "here, and SLURM or HPCCF's cluster file is)")
        sys.exit(
            f"REFUSING to run the search inline: {where}.\n"
            f"  engine={engine}  files={len(files)}  threads={a.threads}\n"
            f"  parallel routing declined because: {why}\n"
            "  Re-run with --sbatch <script>, then submit what it prints under "
            "\"submit with:\", e.g.:\n"
            f"    ... --sbatch ./{engine}_job.sh   # then: bash <out>/submit.sh, or sbatch "
            f"./{engine}_job.sh if that one script was written\n"
            "  (--allow-inline overrides this, e.g. inside an salloc/srun session.)")

    # These size the 5-step chain only. Said out loud when the route is single-shot, where
    # they change nothing: the flag being accepted reads as the jobs having got smaller.
    chain_only = [f for f, v in (("--threads-per-file", a.threads_per_file),
                                 ("--libpred-cpus", a.libpred_cpus),
                                 ("--assembly-cpus", a.assembly_cpus),
                                 ("--assembly-mem", a.assembly_mem),
                                 ("--time-per-file", a.time_per_file),
                                 ("--mem-per-file", a.mem_per_file),
                                 ("--max-simultaneous", a.max_simultaneous)) if v]
    if chain_only and not use_parallel:
        sys.stderr.write(f"[run_search] NOTE: {', '.join(chain_only)} size the 5-step chain "
                         f"only; this search is single-shot ({why}), so they are ignored. "
                         f"Its job(s) request --threads {a.threads} CPUs.\n")
    if exclude and not (use_parallel or radiant_chain or a.sbatch):
        sys.stderr.write(f"[run_search] NOTE: --exclude keeps SLURM jobs off nodes; this search "
                         f"runs here, inline, with no job script, so it changes nothing.\n")
    # --mem sizes the job of a single-job search; said when it changes nothing
    if a.mem and (use_parallel or radiant_chain):
        sizing = ("the 5-step chain's are --mem-per-file (steps 1b, 2, 4) and --assembly-mem "
                  "(steps 3, 5)" if use_parallel else
                  "Radiant's per-file array is sized by radiant_parallel.py itself (--mem-per-file "
                  "per file, --fulcrum-mem for Fulcrum, --libpred-mem for the library), which "
                  "run_search does not pass on: run radiant_parallel.py with them to change it")
        sys.stderr.write("[run_search] NOTE: --mem sets the memory of a single-job search; this "
                         "search is a chain whose generator sizes its own jobs, so it is "
                         f"ignored -- {sizing}.\n")
    elif a.mem and not a.sbatch:
        sys.stderr.write("[run_search] NOTE: --mem sets the memory of a SLURM job; this search "
                         "runs here, inline, with no job script, so it changes nothing.\n")
    # --time likewise
    if a.time and (use_parallel or radiant_chain):
        sys.stderr.write("[run_search] NOTE: --time sets the wall time of a single-job search; "
                         "this search is a chain whose generator sizes its own jobs, so it is "
                         "ignored -- " + ("the 5-step chain's array tasks take --time-per-file"
                                          if use_parallel else
                                          "run radiant_parallel.py with its own options") + ".\n")
    elif a.time and not a.sbatch:
        sys.stderr.write("[run_search] NOTE: --time sets the wall time of a SLURM job; this search "
                         "runs here, inline, with no job script, so it changes nothing.\n")
    # The job size of a single-job search, and where it came from (single_job_sizing). A broad
    # digest's default reaches the route's job script through route_queue, never `queue`.
    sizing = (single_job_sizing(engine, a.params, a.mem, a.time)
              if a.sbatch and not (use_parallel or radiant_chain) else None)
    route_queue = dict(queue)
    if sizing_note(sizing):
        print(sizing_note(sizing), file=sys.stderr, flush=True)
    if sizing:
        # the job's memory against the queue's per-job cap, where SLURM says it (review F11)
        _p, _a, _q = slurm_queue(a.partition, a.account, a.qos)
        _mem = (a.mem and f"{a.mem}G") or sizing["mem"] or JOB_MEM_DEFAULT
        refuse_job_over_queue_cap(_mem, user_limits(_p, _a, _q),
                                  sizing["mem_source"] if (a.mem or sizing["mem"]) else "default")
    if sizing and sizing["digest"]:
        if not a.mem:
            route_queue["mem"] = sizing["mem"]
        if not a.time:
            route_queue["time"] = sizing["hours"]
        n_seq = _fasta_sequences(a.fasta)
        bigger = (f" This database has {n_seq:,} sequences, {n_seq / BROAD_DIGEST_SEQUENCES:.1f}x "
                  f"that one: memory follows the database, so raise --mem with it."
                  if n_seq and n_seq > BROAD_DIGEST_SEQUENCES else "")
        print(f"[run_search] a {sizing['digest']} digest ({sizing['digest_source']}): the job "
              f"asks for {route_queue['mem']} ({'your --mem' if a.mem else 'the default'}) and "
              f"{route_queue['time']} h ({'your --time' if a.time else 'the default'})"
              + ("" if a.mem and a.time else
                 f"; the default for one is {BROAD_DIGEST_MEASURED}") + "."
              + bigger + " --mem GB / --time HOURS set your own.", flush=True)
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
    # The '--' links (planned above): DIA-NN is handed them, every record keeps the inputs as given
    search_files = files
    if links_plan["links"]:
        try:
            linked = _diann_parallel_mod().make_input_links(links_plan)
        except (OSError, ValueError) as e:
            sys.exit(f"[run_search] could not link the inputs whose path holds '--': {e}. "
                     "Nothing was submitted.")
        sys.stderr.write(_diann_parallel_mod().input_links_note(links_plan, "run_search"))
        search_files = [linked.get(f, f) for f in files]
    if use_parallel:
        res = run_diann_parallel(cmd, a.params, search_files, a.fasta, a.out, a.threads, a,
                                 acquisition=bundle.get("acquisition", ""))
        if sbatch_refused:
            # Only now that the chain exists: a failed generation must leave the user's file
            # exactly where it was. run_diann_parallel exits on failure, so reaching here is
            # success.
            moved = set_aside(a.sbatch)
            sbatch_refused["existing_file_moved_to"] = moved and os.path.abspath(moved)
    elif engine == "diann":
        res = run_diann(cmd, a.params, search_files, a.fasta, a.out, a.threads, a.sbatch,
                        acquisition=bundle.get("acquisition", ""), queue=queue, raw_files=files,
                        keratin_sample=keratin_rec["value"])
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
        res = run_sage(cmd, a.params, files, a.fasta, a.out, a.threads, a.sbatch,
                       queue=route_queue, pwiz_sif=a.pwiz_sif)
    elif engine == "fragpipe":
        res = run_fragpipe(cmd, bundle, a.params, files, a.fasta, a.out, a.threads, a.sbatch,
                           queue=route_queue, mono=getattr(a, "fragpipe_mono", None),
                           space=getattr(a, "fragpipe_stage_space", None),
                           database=getattr(a, "fragpipe_database", None))
    elif engine == "radiant":
        # Takes the whole tools dict: it needs the container runtime + image to build
        # bind mounts, and DIA-NN to generate the spectral library.
        # On a cluster with >1 file, split the SERIAL search into a per-file array and
        # rescore once — Radiant's own backend cannot parallelise, but the search is
        # per-file independent, so this is the difference between N x t and ~t.
        if radiant_chain:
            print(f"[run_search] radiant: {len(files)} files on SLURM -> per-file array "
                  f"+ one Fulcrum rescoring job (Radiant's own search is serial)")
            res = run_radiant_parallel(tools, a.params, files, a.fasta, a.out,
                                       a.threads, a, library=a.library)
        else:
            res = run_radiant(tools, a.params, files, a.fasta, a.out, a.threads, a.sbatch,
                              library=a.library, mbr=a.mbr, queue=queue)
    else:
        sys.exit(f"unknown engine {engine}")

    # the most this search can use (cpu_hours.py), for every route that wrote job scripts
    res = search_cpu_hours(res)
    if isinstance(res, dict) and res.get("cpu_hours_reserved") is not None:
        print(f"[run_search] cpu_hours_reserved: {res['cpu_hours_reserved']}  "
              f"gpu_hours_reserved: {res.get('gpu_hours_reserved')} -- {cpu_hours_mod().WHAT}")
    # always record what was run (engine + version + exact command) for reproducibility
    try:
        os.makedirs(a.out, exist_ok=True)
        with open_text(os.path.join(a.out, "search_provenance.json")) as fh:
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
                       # read back by notify_slack.py (instrument + acquisition in the post)
                       "bundle": os.path.abspath(a.bundle),
                       # with the queue, what a re-run of this search needs (sage_lfq_check's
                       # corrected-window command)
                       "tools": os.path.abspath(a.tools), "queue": queue,
                       # the memory and wall time of the single-job search's job(s), as the
                       # job scripts say (EMITTED_JOBS), and where each came from
                       # (single_job_sizing: the user's flag, a broad digest's default -- rule
                       # 2, said so -- or the route's default); null for a chain (its generator
                       # sizes its jobs) and for an inline search (no job)
                       **job_size_record(sizing, EMITTED_JOBS),
                       # what the jobs' end hook will do, as baked into them at generation
                       # (only DIA-NN routes have a completeness guard to stage from)
                       "job_end_hook": _job_end_plan(a, engine),
                       "resolved_params_file": (rp or {}).get("file") or a.params,
                       "resolved_params_produced": ((rp or {}).get("produced")
                                                    or "before the search (the params file as given)"),
                       "resolved_params_note": (rp or {}).get("note"),
                       # scan_window.mode / mass_acc.mode: STABLE values FRAN ingests
                       # (diann_parallel.SCAN_WINDOW_MODES / MASS_ACC_MODES)
                       "scan_window": scan_window_record(engine, a.params, res),
                       "mass_acc": mass_acc_provenance(engine, a.params, res),
                       # --pg-level / --no-prot-inf as the cfg set them, or DIA-NN's default
                       "protein_inference": protein_inference_provenance(engine, a.params),
                       "fasta": a.fasta, "threads": a.threads,
                       # the sidecar describing it (fetch_fasta.sidecar_for: a Philosopher .fas
                       # has none of its own; its target FASTA's) -- what to read for its facts
                       "fasta_meta": fasta_meta_rec,
                       # read by run_de.R (contaminants.R): keratin kept as the analyte
                       "keratin_sample": keratin_rec["value"],
                       "keratin_sample_source": keratin_rec["sample_source"],
                       "keratin_sample_check": keratin_rec,
                       "n_files": len(files), "files": files,
                       # inputs whose path holds '--', searched through links (DIA-NN cuts its
                       # command line at every '--'): {original, link, run_name, renamed}
                       "input_links": links_plan["links"],
                       # inputs listed more than once, each searched once, with the count
                       # (check_report_runs.distinct_inputs); the analysis report shows them
                       "repeated_paths": repeated_paths,
                       "search_mode": "parallel_5step" if use_parallel else "single_shot",
                       "parallel_routing_reason": why,
                       # what was actually written to submit -- for a library-free search
                       # that is <out>/submit.sh, never the --sbatch name it did not write
                       "submitted_sbatch": None if sbatch_refused else (
                           (res.get("submitted") if isinstance(res, dict) else None)
                           or a.sbatch or None),
                       "sbatch_refused": sbatch_refused, "result": res,
                       # the most this search can use: CPUs x time limit x tasks, every job
                       # (cpu_hours.py); null for a search run inline
                       "cpu_hours_reserved": (res.get("cpu_hours_reserved")
                                              if isinstance(res, dict) else None),
                       "gpu_hours_reserved": (res.get("gpu_hours_reserved")
                                              if isinstance(res, dict) else None),
                       # an inline Sage search's LFQ mass-window check (sage_lfq_check.py) and
                       # which lfq.parquet rows the adapter kept; a Sage job's are merged in here
                       # by the job and by --adapt-only. probe_fallback: an inline search's
                       # pre-search probe measured nothing (a job's probe_fallback.py writes it
                       # here itself)
                       # cpu_sizing: the chain's CPUs per array task and how many run at once,
                       # and why (array_task_cpus) -- `threads` above is what was ASKED
                       **search_records(res)},
                      fh, indent=2)
    except Exception as e:
        sys.stderr.write(f"[run_search] could not write search_provenance.json: {e}\n")

    print(json.dumps(res, indent=2))

    if isinstance(res, dict) and res.get("lfq_window_refused"):
        # recorded above, so the re-run command can be built from the record (rerun_plan)
        sys.stderr.write(f"[run_search] {sage_lfq_check_mod().refusal(a.out, res['sage_lfq_check'])}\n")
        sys.exit(1)

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
