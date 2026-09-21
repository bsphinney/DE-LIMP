#!/usr/bin/env python3
"""
detect_acquisition.py  --  Decide DDA vs DIA for each input file, so the
pipeline can auto-route to the right engine (DIA -> DIA-NN, DDA -> Sage by
default; FragPipe selectable for either).

Detection per format (best-effort, with a confidence score):
  Bruker .d   read analysis.tdf (SQLite). Presence of a DiaFrameMsMsInfo /
              DiaFrameMsMsWindowGroups table, or Frames.MsMsType == 9,
              => dia-PASEF. PasefFrameMsMsInfo / MsMsType == 8 => ddaPASEF.
              (Mirrors DE-LIMP's helpers_instrument.R logic.)
  mzML(.gz)   stream isolation windows of MS2 scans. Wide windows (median
              >= 3 Da) over a small repeating set of centers => DIA; narrow
              windows (<= 2 Da) with many distinct centers => DDA.
  Thermo .raw not natively readable here; read through ThermoRawFileParser
              (https://github.com/compomics/ThermoRawFileParser/releases, v1.4.0+).
              Metadata JSON -> instrument model and scan counts; a `query` of a
              mid-run slice of scans -> isolation windows + the filter string's
              data-dependent flag -> DIA/DDA and the acquired bounds. Parser missing
              or failing => 'unknown'/low, a stated reason, and a stderr WARNING.

Every Bruker .d is also checked for a truncated or at-risk analysis.tdf
(bruker_tdf.tdf_integrity): WAL-mode header, a non-empty -wal/-journal beside it,
or a frame index that ends short of analysis.tdf_bin. Anything but `ok` becomes a
line in that file's `warnings` and sets `needs_confirmation`.

Output: JSON to stdout. The orchestrator MUST confirm with the user whenever
confidence != "high" or any file has `warnings` (a Thermo .raw the parser could not read,
a Bruker analysis.tdf that is not `ok`), before launching a multi-hour search. Each file's
`warnings` and per-.raw progress also go to stderr.

On a cluster login node (sbatch on PATH, SLURM_JOB_ID unset) more than LOGIN_NODE_MAX_RAW
Thermo .raw files are refused before any is read: run it under srun, or pass
--allow-login-node.

Usage: python3 detect_acquisition.py [--allow-login-node] FILE [FILE ...]
       THERMORAWFILEPARSER="dotnet /opt/trfp/ThermoRawFileParser.dll" \
       python3 detect_acquisition.py run.raw        # parser not on PATH as one executable
"""
import sys, os, json, glob, gzip, math, sqlite3, statistics, shutil, subprocess, shlex, tempfile, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# Every analysis.tdf is opened through bruker_tdf.connect_tdf (read-only AND immutable):
# a read-write open truncates a tdf with a stale -wal beside it (the state of 342 tdfs on
# HIVE), and mode=ro alone reads through the stale -wal.
from bruker_tdf import connect_tdf, tdf_integrity, integrity_warning  # noqa: E402

# ---------------------------------------------------------------------------
# Instrument extraction (PLAN.md §7d) — mirrors DE-LIMP R/helpers_instrument.R.
# The matcher in fetch_workflows.py scores on instrument; null is allowed and
# falls back to score 0 (user confirms). Best-effort, never fatal.
# ---------------------------------------------------------------------------
def instrument_bruker_d(path):
    tdf = os.path.join(path, "analysis.tdf")
    if not os.path.exists(tdf):
        return None
    try:
        con = connect_tdf(tdf)
        cur = con.cursor()
        # GlobalMetadata is a key/value table; InstrumentName holds e.g. "timsTOF Pro"
        rows = dict(cur.execute("SELECT Key, Value FROM GlobalMetadata"))
        con.close()
        for k in ("InstrumentName", "InstrumentSourceType", "Instrument"):
            if rows.get(k):
                return rows[k]
    except sqlite3.Error:
        return None
    return None

def instrument_thermo_raw(path):
    """Instrument model of a Thermo .raw, e.g. "Orbitrap Exploris 480", or None."""
    meta, _err = trfp_metadata(path)
    return _cv(meta, "InstrumentProperties", CV_MODEL) if meta else None

def detect_instrument(path):
    low = path.lower().rstrip("/")
    if low.endswith(".d"):
        return instrument_bruker_d(path)
    if low.endswith(".raw") and not os.path.isdir(path):     # a .raw FOLDER is Waters
        return instrument_thermo_raw(path)
    return None

def detect_bruker_d(path):
    tdf = os.path.join(path, "analysis.tdf")
    if not os.path.exists(tdf):
        return ("unknown", "low", "no analysis.tdf in .d folder", None)
    con = None
    try:
        con = connect_tdf(tdf)
        cur = con.cursor()
        tables = {r[0] for r in cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        if {"DiaFrameMsMsInfo", "DiaFrameMsMsWindowGroups"} & tables:
            # Also read the ACQUIRED precursor m/z range while we are in here.
            # The isolation windows live in DiaFrameMsMsWindows (IsolationMz,
            # IsolationWidth); DiaFrameMsMsWindowGroups holds only an Id, so
            # detecting on the latter and reading from the former is deliberate
            # -- schema verified against a real analysis.tdf on 2026-08-17.
            # Without this the caller falls back to a hardcoded 380-980, which
            # on a 299.5-1200.5 dia-PASEF method silently discards both tails.
            rng = None
            if "DiaFrameMsMsWindows" in tables:
                try:
                    lo, hi = cur.execute(
                        "SELECT MIN(IsolationMz - IsolationWidth/2.0), "
                        "       MAX(IsolationMz + IsolationWidth/2.0) "
                        "FROM DiaFrameMsMsWindows").fetchone()
                    if lo is not None and hi is not None and hi > lo:
                        rng = (float(lo), float(hi))
                except sqlite3.Error:
                    rng = None
            return ("DIA", "high", "dia-PASEF window tables present", rng)
        if "PasefFrameMsMsInfo" in tables:
            return ("DDA", "high", "PasefFrameMsMsInfo present (ddaPASEF)", None)
        # fall back to MsMsType histogram
        try:
            rows = dict(cur.execute(
                "SELECT MsMsType, COUNT(*) FROM Frames GROUP BY MsMsType"))
            if rows.get(9, 0) > 0:  return ("DIA", "high", "Frames.MsMsType==9", None)
            if rows.get(8, 0) > 0:  return ("DDA", "medium", "Frames.MsMsType==8", None)
        except sqlite3.Error:
            pass
        return ("unknown", "low", "no DIA/DDA markers in tdf", None)
    except sqlite3.Error as e:
        return ("unknown", "low", f"tdf read error: {e}", None)
    finally:
        try: con.close()
        except Exception: pass

def _iter_mzml(path, cap=3000):
    """Yield (ms_level, iso_target, iso_width, lo, hi) for spectra, streaming.

    lo/hi are the ACQUIRED isolation bounds (target -/+ the CV offsets). They are
    what the precursor m/z search range should be derived from; previously only
    the width was kept and the bounds were discarded."""
    import xml.etree.ElementTree as ET
    opn = gzip.open if path.endswith(".gz") else open
    ms_level = iso_target = lower = upper = None
    n = 0
    with opn(path, "rb") as fh:
        for ev, el in ET.iterparse(fh, events=("end",)):
            tag = el.tag.rsplit("}", 1)[-1]
            if tag == "cvParam":
                acc = el.get("accession"); val = el.get("value")
                if acc == "MS:1000511": ms_level = int(float(val))
                elif acc == "MS:1000827" and val: iso_target = float(val)
                elif acc == "MS:1000828" and val: lower = float(val)
                elif acc == "MS:1000829" and val: upper = float(val)
            elif tag == "spectrum":
                if ms_level == 2:
                    w = (lower or 0) + (upper or 0)
                    lo = hi = None
                    if iso_target is not None:
                        lo = iso_target - (lower or 0)
                        hi = iso_target + (upper or 0)
                    yield (2, iso_target, w if w else None, lo, hi)
                    n += 1
                ms_level = iso_target = lower = upper = None
                el.clear()
                if n >= cap: break

def _cap_medium(conf):
    """`high` is a claim about the READ as well as the verdict: downgrade it."""
    return "medium" if conf == "high" else conf


def classify_isolation_windows(widths, targets, los, his):
    """DIA vs DDA from MS2 isolation windows -> (kind, confidence, reason, range).

    The one definition of the width/centre rule, shared by the mzML reader and the Thermo
    .raw reader so the two formats cannot drift into classifying the same run differently.
    `targets` are window centres, `los`/`his` the ACQUIRED window edges.

    A width needs only the two offsets, while an EDGE PAIR also needs the window centre
    (MS:1000827), so `los` can cover fewer windows than `widths`: a valid mzML whose upper
    windows carry no isolation target, or a parser that answers with the first part of the
    scans asked for. The range is then built from a SUBSET of the method and comes out
    clipped -- 350.0-793.0 for an Exploris method that acquired 350.0-1201.0 -- and
    estimate_params.py would emit that as `--max-pr-mz 793` tagged "measured from the
    acquired isolation windows". A wrong number wearing the measured label is worse than the
    FALLBACK tag it replaces, so a partial read never carries `high`, and the reason says how
    much of the method it saw.
    """
    if not widths:
        return ("unknown", "low", "no MS2 isolation windows found", None)
    targets = [round(t * 2) / 2 for t in targets]
    med = statistics.median(widths)
    distinct = len(set(targets))
    n = len(widths)
    rng = (min(los), max(his)) if los else None
    partial = len(los) != n
    gap = (f"; isolation window edges for only {len(los)} of {n} windows, so any range read "
           f"from them may be clipped" if partial else "")
    if med >= 3.0 and distinct <= max(80, n // 50):
        conf = "high" if med >= 4 else "medium"
        return ("DIA", _cap_medium(conf) if partial else conf,
                f"median isolation width {med:.1f} Da over {distinct} distinct centers{gap}",
                rng)
    if med <= 2.0:
        conf = "high" if distinct > n // 10 else "medium"
        return ("DDA", _cap_medium(conf) if partial else conf,
                f"median isolation width {med:.1f} Da, {distinct} distinct precursors{gap}",
                None)
    return ("unknown", "low",
            f"ambiguous: median width {med:.1f} Da, {distinct} centers{gap}", None)

def detect_mzml(path):
    widths, targets, los, his = [], [], [], []
    try:
        for lvl, tgt, w, lo, hi in _iter_mzml(path):
            if w: widths.append(w)
            if tgt is not None: targets.append(tgt)
            if lo is not None and hi is not None and hi > lo:
                los.append(lo); his.append(hi)
    except Exception as e:
        return ("unknown", "low", f"mzML parse error: {e}", None)
    return classify_isolation_windows(widths, targets, los, his)


# ---------------------------------------------------------------------------
# Thermo .raw, read through ThermoRawFileParser (TRFP).
#
# Public source (golden rule 7): https://github.com/compomics/ThermoRawFileParser/releases
# -- self-contained Linux/macOS/Windows builds, a .NET 8 framework-dependent DLL, and
# `conda install -c bioconda thermorawfileparser`. Nothing here assumes HIVE.
#
# The command lines below are the only ones this parser accepts. The previous code called
# `ThermoRawFileParser metadata -i <raw>` and `ThermoRawFileParser query -i <raw>`; on
# TRFP 2.0.0.0 (HIVE srun job 23510571) both exit 255 -- "Unexpected extra arguments"
# (there is no `metadata` subcommand) and "specify a valid scan range" (query needs -n).
# The exit codes were ignored and the empty stdout grepped for "dia", so EVERY .raw came
# back unknown / instrument null / range null and estimate_params.py searched 380-980 on
# methods that acquired 350.0-1201.0 (Exploris 480) and 350.05-1200.95 (Fusion Lumos).
#
# Grammar read from MainClass.cs for every release v1.4.0 .. v.2.0.0-dev, and run for real on
# 2.0.0.0 (the only build exercised on data -- 1.4.x needs Mono). The README states that
# optional parameters "only work in the -option=value format", so every value is inlined.
#   metadata:  -i=<raw> -m=0 -o=<dir>            -> <dir>/<stem>-metadata.json
#              (-m alone writes no spectra: the indexed-mzML default applies only when
#               neither -m nor -f is given. `-f=4` means the same but is documented only
#               from 1.4.4, so it is left out rather than relied on.)
#   spectra:   query -i=<raw> -n=<a>-<b> -b=<file>  -> JSON PROXI spectra with isolation
#              target + lower/upper offsets and the filter string (see CV_FILTER_* below)
# ---------------------------------------------------------------------------
TRFP_RELEASES = "https://github.com/compomics/ThermoRawFileParser/releases"
TRFP_ENV = "THERMORAWFILEPARSER"     # a full command, e.g. "dotnet /opt/trfp/ThermoRawFileParser.dll"
TRFP_NAMES = ("ThermoRawFileParser", "thermorawfileparser")   # release binary / bioconda link
# Measured on HIVE for 3.5 GB raws: metadata 2.6-4.8 s; a query of 200-2000 scans 1.1-5.8 s
# (10.8 s once, for the first 1000 scans of a DDA run). The metadata call walks every scan
# header, so a slow network mount can take far longer -- but a parser that has not answered
# in 5 min is stuck, and stuck must not look like done.
TRFP_TIMEOUT_S = 300
# JSON metadata terms, keyed by ACCESSION: TRFP's own `name` strings are not stable
# (the source labels the MS2 count "Number of MS1 spectra").
CV_MODEL = "MS:1000494"          # InstrumentProperties: Thermo Scientific instrument model
CV_SCAN_RANGE = "PRIDE:0000479"  # ScanSettings: "first:last"
CV_N_MS1 = "PRIDE:0000481"       # MsData
CV_N_MS2 = "PRIDE:0000482"       # MsData
# `query` spectrum attributes. The filter string carries the instrument's data-dependent
# flag, and its accession depends on the build: every release v1.3.0-v1.4.4 writes it as
# "MS:10000512" -- one zero too many -- and v1.4.5 corrected it to MS:1000512
# (Query/ProxiSpectrumReader.cs at each tag). Bioconda still serves 1.3.2-1.4.4. Reading only
# the correct spelling silently loses the flag on those builds, and with it the check that
# stops narrow-window DIA being called DDA (300 x 2 m/z DIA came back DDA/high, no range, no
# confirmation). ms level and the isolation offsets are spelled correctly in every release.
CV_FILTER = "MS:1000512"
CV_FILTER_PRE_1_4_5 = "MS:10000512"
# NOT used, on purpose: MsData "MS min MZ"/"MS max MZ" (PRIDE:0000476/7) are the lowest and
# highest isolation window CENTRES (367.5/1183.5 on the Exploris run above, whose windows
# span 350.0-1201.0). Taking them as the range clips half a window off each end.

# How many scans to query. A DIA cycle is one MS1 plus its MS2 windows, so the metadata's
# MS2/MS1 count ratio is the cycle length (78352/3135 -> 26 scans on the Exploris). Four
# cycles guarantees every window is seen, including staggered schemes that alternate
# window sets between cycles; the floor covers methods with extra MS1 scans per cycle, the
# ceiling bounds the output (TRFP writes every peak: 1000 mid-run Exploris DIA scans were
# 43 MB of JSON, 200 are ~9 MB).
QUERY_CYCLES, QUERY_MIN_SCANS, QUERY_MAX_SCANS = 4, 200, 1500
QUERY_SCANS_WITHOUT_METADATA = 1000

FALLBACK_CONSEQUENCE = ("precursor m/z range NOT measured, so estimate_params.py will search "
                        "its 380-980 FALLBACK -- pass --precursor-mz-range if the method "
                        "acquired wider")


def find_trfp():
    """The ThermoRawFileParser command as an argv prefix, or None.

    $THERMORAWFILEPARSER wins, because two public distributions are not one executable on
    PATH: the framework-dependent release is `dotnet ThermoRawFileParser.dll`, and 1.4.x on
    Linux/macOS is `mono ThermoRawFileParser.exe`.
    """
    env = os.environ.get(TRFP_ENV, "").strip()
    if env:
        toks = shlex.split(env, posix=(os.name != "nt"))
        return [t[1:-1] if len(t) > 1 and t[0] == t[-1] == '"' else t for t in toks]
    for name in TRFP_NAMES:
        hit = shutil.which(name)
        if hit:
            return [hit]
    return None


def _trfp_not_found():
    return (f"ThermoRawFileParser not found (${TRFP_ENV} is unset and neither "
            f"{' nor '.join(TRFP_NAMES)} is on PATH), so this .raw was not read: acquisition "
            f"and instrument unknown, and {FALLBACK_CONSEQUENCE}. Install it from "
            f"{TRFP_RELEASES} (self-contained Linux/macOS/Windows builds, or "
            f"`conda install -c bioconda thermorawfileparser`), or convert the .raw to mzML")


def _trfp_error_lines(stdout):
    """log4net `ERROR` lines from the parser's console appender (stdout)."""
    return [ln.strip() for ln in (stdout or "").splitlines() if " ERROR " in f" {ln} "]


def _trfp_message(stdout, stderr):
    """The parser's own words. Usage errors go to stderr; processing errors are log4net
    `ERROR` lines on stdout (its console appender), so look in both."""
    usage = [ln.strip() for ln in (stderr or "").splitlines() if ln.strip()]
    msg = " | ".join(_trfp_error_lines(stdout)[:2] + usage[:1]) or "no message"
    return msg if len(msg) <= 400 else msg[:400] + "..."


def _run_trfp(cmd, args):
    """(True, None, stdout) or (False, "exit N: <parser message>", stdout). Never raises.

    stdout comes back on SUCCESS too, because exit 0 is not proof the parser did the work:
    a processing error is a log4net `ERROR` line on the console appender and need not move
    the exit code (see _trfp_message). A caller that reads only the exit code accepts a
    partial answer as a complete one.
    """
    try:
        p = subprocess.run(cmd + args, capture_output=True, text=True, errors="replace",
                           timeout=TRFP_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return False, f"no answer after {TRFP_TIMEOUT_S} s", ""
    except OSError as e:                       # includes FileNotFoundError / PermissionError
        return False, f"could not run {' '.join(cmd)!r}: {e}", ""
    if p.returncode != 0:
        return False, f"exit {p.returncode}: {_trfp_message(p.stdout, p.stderr)}", p.stdout
    return True, None, p.stdout


_TRFP_VERSIONS = {}

def trfp_version(cmd):
    """`--version` output (e.g. "2.0.0.0"), cached per command; None if it will not say."""
    key = tuple(cmd)
    if key not in _TRFP_VERSIONS:
        try:
            p = subprocess.run(cmd + ["--version"], capture_output=True, text=True,
                               errors="replace", timeout=60)
            lines = [ln.strip() for ln in p.stdout.splitlines() if ln.strip()]
            _TRFP_VERSIONS[key] = lines[-1] if p.returncode == 0 and lines else None
        except (OSError, subprocess.TimeoutExpired):
            _TRFP_VERSIONS[key] = None
    return _TRFP_VERSIONS[key]


def _cv(meta, section, accession):
    for term in (meta or {}).get(section) or []:
        if isinstance(term, dict) and term.get("accession") == accession:
            v = term.get("value")
            return v.strip() if isinstance(v, str) and v.strip() else None
    return None


def trfp_metadata(path, cmd=None):
    """(metadata dict, None) or (None, why)."""
    cmd = cmd or find_trfp()
    if not cmd:
        return None, _trfp_not_found()
    with tempfile.TemporaryDirectory(prefix="trfp_meta_") as tmp:
        ok, why, stdout = _run_trfp(cmd, [f"-i={path}", "-m=0", f"-o={tmp}"])
        if not ok:
            return None, f"ThermoRawFileParser metadata call failed ({why})"
        if _trfp_error_lines(stdout):
            return None, ("ThermoRawFileParser metadata call reported an error and still "
                          f"exited 0 ({_trfp_message(stdout, None)})")
        hits = [f for f in os.listdir(tmp) if f.endswith("-metadata.json")]
        if not hits:
            return None, "ThermoRawFileParser metadata call wrote no *-metadata.json"
        try:
            with open(os.path.join(tmp, hits[0]), encoding="utf-8-sig") as fh:
                meta = json.load(fh)
        except (OSError, ValueError) as e:
            return None, f"ThermoRawFileParser metadata JSON unreadable ({e})"
    if not isinstance(meta, dict):
        return None, "ThermoRawFileParser metadata JSON is not an object"
    return meta, None


def trfp_query(path, first, last, cmd):
    """(list of PROXI spectra, None, [notes]) or (None, why, [notes]).

    `notes` are what did not stop the read but must reach the user, because they change what
    the answer is worth: chiefly an answer SHORTER than the scans asked for. The scan range
    is clamped to the run's own (query_scan_range), so a short answer means the parser
    stopped early -- the first 13 of an Exploris method's 25 windows read as a complete
    350.0-793.0 method. The caller refuses `high` on any note.
    """
    notes = []
    with tempfile.TemporaryDirectory(prefix="trfp_query_") as tmp:
        dest = os.path.join(tmp, "query.json")
        ok, why, stdout = _run_trfp(cmd, ["query", f"-i={path}", f"-n={first}-{last}",
                                          f"-b={dest}"])
        if not ok:
            return None, f"ThermoRawFileParser query of scans {first}-{last} failed ({why})", notes
        if _trfp_error_lines(stdout):
            # exit 0 with an ERROR line: the parser hit a problem, wrote what it had, and
            # said so only on stdout. Whatever is in that file is not the slice we asked for.
            return None, (f"ThermoRawFileParser query of scans {first}-{last} reported an "
                          f"error and still exited 0 ({_trfp_message(stdout, None)})"), notes
        try:
            with open(dest, encoding="utf-8-sig") as fh:
                spectra = json.load(fh)
        except (OSError, ValueError) as e:
            return None, f"ThermoRawFileParser query output for scans {first}-{last} unreadable ({e})", notes
    if not isinstance(spectra, list):
        return None, "ThermoRawFileParser query output is not a JSON list of spectra", notes
    asked = last - first + 1
    if len(spectra) < asked:
        notes.append(
            f"ThermoRawFileParser answered the query of scans {first}-{last} with "
            f"{len(spectra)} of the {asked} scans asked for, so only part of the requested "
            f"slice was read: the acquisition and any precursor m/z range below come from "
            f"that part and may miss windows the method also acquired")
    return spectra, None, notes


def query_scan_range(meta):
    """(first, last, note) scans to query: several acquisition cycles from the middle of the
    run, where the method is at steady state -- not the loading/wash start.

    `note` is None, or says that the slice is NOT mid-run after all. Without the run's own
    scan range there is nothing to take a middle of, so the fallback reads from scan 1 --
    the void volume -- and no other message says so: "instrument unknown" is about the
    instrument, not about which part of the run was classified.
    """
    not_mid_run = (
        f"the run's scan range was not readable from the ThermoRawFileParser metadata, so "
        f"scans 1-{QUERY_SCANS_WITHOUT_METADATA} were read -- the START of the run (void "
        f"volume and equilibration), not a mid-run slice at steady state. A method segment "
        f"that differs there, or a start with no MS2 at all, is classified as the whole run")
    try:
        first, last = (int(v) for v in _cv(meta, "ScanSettings", CV_SCAN_RANGE).split(":"))
    except (AttributeError, ValueError):
        return 1, QUERY_SCANS_WITHOUT_METADATA, not_mid_run
    if last < first:
        return 1, QUERY_SCANS_WITHOUT_METADATA, not_mid_run
    try:
        cycle = int(_cv(meta, "MsData", CV_N_MS2)) / int(_cv(meta, "MsData", CV_N_MS1)) + 1
        want = min(QUERY_MAX_SCANS, max(QUERY_MIN_SCANS, math.ceil(QUERY_CYCLES * cycle)))
    except (TypeError, ValueError, ZeroDivisionError):
        want = QUERY_SCANS_WITHOUT_METADATA      # still mid-run, just a fixed-width slice
    a = max(first, (first + last) // 2 - want // 2)
    return a, min(last, a + want - 1), None


def thermo_isolation_windows(spectra):
    """MS2 isolation windows + data-dependent evidence from TRFP `query` JSON.

    Edges are target -/+ the lower/upper offsets TRFP reports, which already fold in the
    instrument's isolation-width offset -- the same values it writes to mzML, so the two
    readers agree: detect_mzml() on the full TRFP 2.0.0.0 mzML of the runs above gave
    350.0-1201.0 and 350.04999389648435-1200.9499755859374 (HIVE srun jobs 23508203/23508208),
    this reader 350.0-1201.0 and 350.05-1200.95 after rounding (job 23511567).
    """
    w = {"widths": [], "centres": [], "los": [], "his": [],
         "n_ms2": 0, "n_filter": 0, "n_dependent": 0}
    for s in spectra:
        attrs = {}
        for a in (s.get("attributes") or []) if isinstance(s, dict) else []:
            if isinstance(a, dict) and a.get("accession"):
                attrs.setdefault(a["accession"], a.get("value"))
        try:
            if int(float(attrs.get("MS:1000511"))) != 2:
                continue
        except (TypeError, ValueError):
            continue
        w["n_ms2"] += 1
        filt = attrs.get(CV_FILTER) or attrs.get(CV_FILTER_PRE_1_4_5)
        if isinstance(filt, str) and filt.strip():
            w["n_filter"] += 1
            # Thermo filter grammar: a standalone `d` token before the mass list means
            # "data-dependent scan" ("FTMS + c NSI d Full ms2 572.3181@hcd30.00 [...]");
            # DIA windows carry none ("FTMS + p NSI Full ms2 367.5000@hcd30.00 [...]").
            if "d" in filt.split("[", 1)[0].split():
                w["n_dependent"] += 1
        try:
            tgt = float(attrs["MS:1000827"])
            lo_off = float(attrs.get("MS:1000828") or 0)
            hi_off = float(attrs.get("MS:1000829") or 0)
        except (KeyError, TypeError, ValueError):
            continue
        if lo_off + hi_off <= 0:
            continue
        w["widths"].append(lo_off + hi_off)
        w["centres"].append(tgt)
        w["los"].append(tgt - lo_off)
        w["his"].append(tgt + hi_off)
    return w


def classify_thermo_windows(w):
    """(kind, confidence, reason, range) from thermo_isolation_windows() output.

    The width/centre rule first (the same one mzML gets), then the instrument's own
    data-dependent flag, which outranks window shape whenever the parser reports it.
    """
    kind, conf, why, rng = classify_isolation_windows(w["widths"], w["centres"],
                                                      w["los"], w["his"])
    # Every MS2 scan the parser answered with, against the ones a window could be read from.
    # A shortfall is the same defect classify_isolation_windows guards against for mzML --
    # a method judged, and a range measured, from part of its windows.
    unread = w["n_ms2"] - len(w["widths"])
    capped = unread > 0          # things about the READ that forbid `high`, applied last so
    if unread > 0:               # no later branch can promote a partial read back to high
        why += (f"; {unread} of {w['n_ms2']} MS2 scans in the slice had no readable "
                f"isolation window")
    if w["n_filter"] == 0:
        # No filter string, so the instrument's own flag is unavailable and
        # references/search-engines.md's "none flagged -> never DDA" cannot be applied: the
        # width rule decided alone, and a 300 x 2 m/z DIA method looks exactly like DDA to
        # it. Every shipping release writes the string (v1.3.0-v1.4.4 under the misspelled
        # MS:10000512, v1.4.5+ correctly), so a parser that does not is one this code has
        # never been run against. That is not a `high` answer.
        dep = ("data-dependent flag not available (the parser reported no filter string, "
               f"{CV_FILTER} or {CV_FILTER_PRE_1_4_5}), so the window widths decided alone")
        return kind, _cap_medium(conf), f"{why}; {dep}", rng
    dep = f"{w['n_dependent']}/{w['n_filter']} MS2 scans flagged data-dependent"
    if w["n_dependent"] == 0 and kind != "DIA" and w["widths"]:
        # Nothing was data-dependent, so this is not DDA whatever the widths look like:
        # narrow-window DIA (2-3 m/z isolation, Astral-style) or targeted PRM. Medium, so
        # the user confirms which.
        kind, conf = "DIA", "medium"
        rng = (min(w["los"]), max(w["his"])) if w["los"] else None
        dep += " -- so not DDA despite the window widths (narrow-window DIA or PRM?)"
    elif w["n_dependent"] == w["n_filter"]:
        if kind == "DDA":
            conf = "high"                      # instrument flag and window shape agree
        else:
            kind, conf, rng = "DDA", "medium", None
            # There is no "despite the widths" to talk about when there were no widths.
            dep += (" -- so DDA despite the window widths" if w["widths"] else
                    " -- so DDA on the flag alone: no isolation windows were readable")
    elif w["n_dependent"]:
        conf = "medium" if conf == "high" else conf
        dep += " -- mixed dependent and independent MS2 scans (hybrid method?)"
    return kind, (_cap_medium(conf) if capped else conf), f"{why}; {dep}", rng


def read_thermo_raw(path):
    """Everything detect_acquisition reports for one Thermo .raw, from at most two TRFP calls.

    Returns a dict: acquisition, confidence, reason, precursor_mz_range, instrument,
    warnings (every read failure, in words), reader.
    """
    out = {"acquisition": "unknown", "confidence": "low", "reason": "",
           "precursor_mz_range": None, "instrument": None, "warnings": [], "reader": None}
    cmd = find_trfp()
    if not cmd:
        out["reason"] = _trfp_not_found()
        out["warnings"].append(out["reason"])
        return out
    version = trfp_version(cmd)
    out["reader"] = f"ThermoRawFileParser {version or '(version unknown)'} [{' '.join(cmd)}]"

    meta, meta_err = trfp_metadata(path, cmd)
    if meta is None:
        # Not fatal for DIA/DDA, but the instrument decides mass accuracy downstream.
        out["warnings"].append(f"instrument unknown: {meta_err}")
    else:
        out["instrument"] = _cv(meta, "InstrumentProperties", CV_MODEL)
        if not out["instrument"]:
            out["warnings"].append("instrument unknown: ThermoRawFileParser metadata has no "
                                   f"instrument model ({CV_MODEL})")

    first, last, slice_note = query_scan_range(meta)
    if slice_note:
        out["warnings"].append(slice_note)
    spectra, q_err, q_notes = trfp_query(path, first, last, cmd)
    if spectra is None:
        out["reason"] = f"{q_err}: acquisition unknown and {FALLBACK_CONSEQUENCE}"
        out["warnings"].append(out["reason"])
        return out
    out["warnings"].extend(q_notes)

    w = thermo_isolation_windows(spectra)
    kind, conf, why, rng = classify_thermo_windows(w)
    if q_notes:
        conf = _cap_medium(conf)     # a partial answer is not a high-confidence read
    if w["n_filter"] == 0:
        out["warnings"].append(
            f"ThermoRawFileParser reported no filter string ({CV_FILTER} or "
            f"{CV_FILTER_PRE_1_4_5}) for any MS2 scan, so the instrument's own "
            f"data-dependent flag could not be read and the isolation widths decided DIA vs "
            f"DDA alone -- confidence is capped at medium: confirm the acquisition before "
            f"searching, and note that a narrow-window (2-3 m/z) DIA method is "
            f"indistinguishable from DDA on widths alone")
    if rng:
        # TRFP serialises the isolation target as a float32 -- 372.8999938964844 for the
        # Lumos method's 372.9 -- so the raw edge comes out 350.0499938964844. Near m/z 1200
        # a float32 is only good to ~1e-4, so three decimals is all the precision there is,
        # and it returns the method's own 350.05 / 1200.95.
        rng = (round(rng[0], 3), round(rng[1], 3))
    out.update(acquisition=kind, confidence=conf, precursor_mz_range=rng,
               reason=f"ThermoRawFileParser query of scans {first}-{last}: {why}")
    if kind == "DIA" and not rng:
        out["warnings"].append(f"DIA, but no isolation window edges were readable: "
                               f"{FALLBACK_CONSEQUENCE}")
    elif kind == "unknown":
        out["warnings"].append(f"{out['reason']}: {FALLBACK_CONSEQUENCE}")
    return out


def detect_thermo_raw(path):
    """(kind, confidence, reason, precursor m/z range) -- the shape every detector returns."""
    t = read_thermo_raw(path)
    return (t["acquisition"], t["confidence"], t["reason"], t["precursor_mz_range"])


def classify(path):
    p = path.rstrip("/")
    low = p.lower()
    mz_range = None
    warnings, reader, instrument = [], None, None
    if low.endswith(".d"):
        kind, conf, why, mz_range = detect_bruker_d(p); vendor = "Bruker"
    elif low.endswith((".mzml", ".mzml.gz")):
        kind, conf, why, mz_range = detect_mzml(p); vendor = "mzML"
    elif low.endswith(".raw") and os.path.isdir(p):
        # Waters .raw is a FOLDER; handing it to a Thermo reader only produces a confusing
        # "not a valid RAW file" from the wrong vendor's tool.
        kind, conf, why, vendor = ("unknown", "low",
                                   "Waters .raw folder: convert to mzML to classify", "Waters")
    elif low.endswith(".raw"):
        t = read_thermo_raw(p); vendor = "Thermo"
        kind, conf, why, mz_range = (t["acquisition"], t["confidence"], t["reason"],
                                     t["precursor_mz_range"])
        # instrument comes from the same metadata call -- do not run the parser twice
        instrument, warnings, reader = t["instrument"], t["warnings"], t["reader"]
    elif low.endswith((".wiff", ".wiff2")):
        kind, conf, why, vendor = "unknown", "low", "SCIEX .wiff: convert to mzML to classify", "SCIEX"
    else:
        kind, conf, why, vendor = "unknown", "low", "unrecognized extension", "?"
    if vendor != "Thermo":
        instrument = detect_instrument(p)
    # A .d whose analysis.tdf indexes only part of its tdf_bin is searched without an
    # error -- DIA-NN reads what the index points at. Nothing downstream can tell, so it
    # is checked here, before the engine is chosen. None for anything but a Bruker .d.
    integrity = (tdf_integrity(p) if vendor == "Bruker"
                 and os.path.exists(os.path.join(p, "analysis.tdf")) else None)
    if integrity and integrity["status"] != "ok":
        # the same per-file list the Thermo reader fills: one place to look, one gate
        warnings = warnings + [integrity_warning(integrity)]
    return {"file": p, "vendor": vendor, "acquisition": kind,
            "confidence": conf, "reason": why, "instrument": instrument,
            # ACQUIRED precursor m/z bounds, or null when the format cannot tell
            # us. estimate_params.py searches this range instead of a hardcoded
            # 380-980 -- see its --precursor-mz-range flag.
            "precursor_mz_range": (list(mz_range) if mz_range else None),
            # Bruker .d only (None otherwise): bruker_tdf.tdf_integrity() of its analysis.tdf
            "tdf_integrity": integrity,
            # every problem with this file the user must hear before a search starts, in
            # words: a Thermo .raw read failure or a Bruker analysis.tdf that is not ok
            # (also printed to stderr; any one sets needs_confirmation)
            "warnings": warnings,
            # the external reader and its version, when one was needed (Thermo .raw)
            "reader": reader}

# ---------------------------------------------------------------------------
# A cohort of .raw on a cluster login node.
#
# Reading a .raw is two parser processes that read the raw. Measured on HIVE compute nodes
# (srun jobs 23511567, 23512661, 23513356, 23515380; TRFP 2.0.0.0, 3.5 GB raws): 3.1-7.0 s
# per file, ~1 core (68-98% CPU), 286-354 MB peak RSS, a 1.2-8.8 MB query JSON in $TMPDIR
# (deleted per file); 6 files in one call took 27-33 s. With the old command line the calls died at
# argument parsing, so step 2 cost nothing -- and SKILL.md 0a drives step 2 over SSH, i.e.
# on the login node. A 200-file cohort there is ~20 min of parser I/O on a shared host,
# which golden rule #3 forbids. The limit is a choice, not a measurement: 5 files keeps the
# login-node work to about half a minute. The login-node test is run_search.py's own
# (sbatch present, not inside a SLURM job), so the two refuse on the same hosts. No
# scheduler => no refusal: a laptop is not a shared login node, and gets progress lines
# instead. No parser => no refusal either: each .raw fails at once with the not-found
# warning, and a detour through srun would only delay that message.
# ---------------------------------------------------------------------------
LOGIN_NODE_MAX_RAW = 5
ALLOW_LOGIN_NODE = "--allow-login-node"
RAW_SECONDS_EACH = 7            # upper end of the per-file measurement above


def is_thermo_raw(path):
    p = path.rstrip("/")
    return p.lower().endswith(".raw") and not os.path.isdir(p)     # a .raw FOLDER is Waters


def on_cluster_login_node():
    """sbatch present + no SLURM_JOB_ID == a submit host, not a compute node."""
    return shutil.which("sbatch") is not None and not os.environ.get("SLURM_JOB_ID")


def login_node_refusal(n_raw):
    secs = n_raw * RAW_SECONDS_EACH
    est = f"{secs} s" if secs < 120 else f"{math.ceil(secs / 60)} min"
    minutes = max(10, math.ceil(2 * secs / 60))
    return (
        f"REFUSING to read {n_raw} Thermo .raw files on what looks like a cluster login/submit "
        f"node (sbatch is on PATH and SLURM_JOB_ID is unset).\n"
        f"  Each .raw is two ThermoRawFileParser processes reading the raw -- measured ~1 core, "
        f"~350 MB and 3-7 s per file -- so up to ~{est} of I/O here, on a shared host "
        f"(SKILL.md golden rule #3).\n"
        f"  Run it on a compute node instead, e.g.:\n"
        f"    srun --cpus-per-task=1 --mem=2G --time={minutes} python3 "
        f"{os.path.abspath(__file__)} <the same files>\n"
        f"  adding the --account/--partition/--qos your cluster requires "
        f"(`sacctmgr show assoc user=$USER` lists yours),\n"
        f"  or re-run with {ALLOW_LOGIN_NODE} to read them here anyway "
        f"({LOGIN_NODE_MAX_RAW} or fewer need no flag).")


def main(argv):
    allow_login_node = ALLOW_LOGIN_NODE in argv
    argv = [a for a in argv if a != ALLOW_LOGIN_NODE]
    usage = f"usage: detect_acquisition.py [{ALLOW_LOGIN_NODE}] FILE [FILE ...]"
    if not argv:
        print(usage, file=sys.stderr)
        sys.exit(2)
    # Anything dash-led that is not that one flag used to fall through to the glob, match
    # nothing, and enter files[] as itself: `--allow_login_node` and `--allow-login-node=1`
    # came back as an "unrecognized extension" file, which turned `overall` into "mixed" and
    # asked the user to confirm a phantom input -- while the login-node guard they were
    # trying to lift stayed on. A typo must not look like another input file.
    unknown = [a for a in argv if a.startswith("-")]
    if unknown:
        print(f"detect_acquisition.py: unknown option(s): {' '.join(unknown)} "
              f"(the only option is {ALLOW_LOGIN_NODE}; a file whose name starts with '-' "
              f"can be passed as ./{unknown[0]})", file=sys.stderr)
        print(usage, file=sys.stderr)
        sys.exit(2)
    files = []
    for a in argv:
        files.extend(sorted(glob.glob(a)) or [a])
    n_raw = sum(1 for f in files if is_thermo_raw(f))
    # Refuse BEFORE reading anything: a refusal after the first 100 files is no refusal.
    if (n_raw > LOGIN_NODE_MAX_RAW and on_cluster_login_node() and not allow_login_node
            and find_trfp()):
        sys.exit(login_node_refusal(n_raw))
    results, k = [], 0
    for f in files:
        if not is_thermo_raw(f):
            results.append(classify(f))
            continue
        k += 1
        t0 = time.monotonic()
        results.append(classify(f))
        # Seconds per file, as it happens: minutes of parser work in silence looks like a hang.
        print(f"[detect_acquisition] Thermo .raw {k}/{n_raw} read in "
              f"{time.monotonic() - t0:.1f} s: {f}", file=sys.stderr, flush=True)
    # stdout is JSON for the caller; problems also go to stderr, where a person sees them
    # even when a script only keeps the JSON.
    for r in results:
        for w in r.get("warnings") or []:
            print(f"[detect_acquisition] WARNING: {r['file']}: {w}", file=sys.stderr)
    kinds = {r["acquisition"] for r in results}
    overall = (next(iter(kinds)) if len(kinds) == 1 else "mixed")
    low_conf = [r["file"] for r in results if r["confidence"] != "high"]
    instruments = sorted({r["instrument"] for r in results if r.get("instrument")})
    # one instrument string for the matcher: unambiguous only if all files agree
    instrument = instruments[0] if len(instruments) == 1 else ""
    # Union the per-file acquired ranges: search everything any run acquired.
    # A narrower intersection would silently drop data from the wider methods,
    # which is the failure this whole field exists to prevent.
    ranges = [r["precursor_mz_range"] for r in results if r.get("precursor_mz_range")]
    mz_range = [min(x[0] for x in ranges), max(x[1] for x in ranges)] if ranges else None
    mixed_ranges = len({tuple(x) for x in ranges}) > 1
    print(json.dumps({
        "overall": overall,
        "instrument": instrument,
        "instruments_seen": instruments,
        # Feed straight to estimate_params.py --precursor-mz-range LO HI.
        "precursor_mz_range": mz_range,
        "precursor_mz_range_mixed": mixed_ranges,
        "precursor_mz_range_files_without": [
            r["file"] for r in results if not r.get("precursor_mz_range")],
        "needs_confirmation": (bool(low_conf) or overall in ("mixed", "unknown")
                               or len(instruments) > 1
                               or any(r.get("warnings") for r in results)),
        "low_confidence_files": low_conf,
        # No top-level list of problem files: files[].warnings is the one place a problem is
        # reported, whichever reader found it (a vendor-specific list reads "nothing wrong"
        # whenever the problem is another vendor's). tdf status: files[].tdf_integrity.
        "files": results,
    }, indent=2))

if __name__ == "__main__":
    main(sys.argv[1:])
