#!/usr/bin/env python3
"""probe_window.py -- measure the scan-window radius DIA-NN infers for an acquisition.

WHY THIS EXISTS
---------------
The 5-step parallel chain reuses .quant files across steps. DIA-NN warns:

    WARNING: combining reuse of .quant files with automatic optimisation of mass
    accuracies or scan window will lead to results that are different from those of
    the original analysis that produced the .quant files and is strongly not
    recommended

Mass accuracy is pinned by estimate_params.py. The scan window was NOT: the flag was
omitted, so DIA-NN optimised it PER FILE. On a real 18-file poplar run that produced a
radius of 7 for seventeen files and 8 for one, which the chain then stitched together.

The radius is a property of the ACQUISITION SCHEME (cycle time vs chromatographic peak
width), not of the individual sample -- but it must be MEASURED, not guessed, because it
depends on the gradient and method.

WHICH RUNS
----------
Not simply the first one. DIA-NN's README ("Changing default settings") says its automatic
optimisation "is inherently noisy: even replicate injections may not produce identical
results, and therefore the analysis results will depend on which run is first in the list",
and recommends running "on several representative runs". Measured on HIVE (DIA-NN 2.7.0,
2026-09-16): one-file probes of one timsTOF cohort gave 10, 11 or 14 depending on the run,
while each run's own radius was reproducible. So, given the whole cohort (--raw ... /
--raw-list):

  * never probe a run whose data DIA-NN cannot fully read (BRUKER INDEX, below);
  * never probe a run under HALF the size of a typical run (a .d: the indexed bytes of its
    analysis.tdf_bin), where "typical" is the median of the LARGER half of the runs (see
    MIN_FRACTION) -- blanks, washes, failed acquisitions;
  * rank what is left by how much was acquired: when every run is a .d with a frame index, by
    acquisition time (max Frames.Time), with the same half-of-typical floor on it (a short wash
    the size of a real run); otherwise by size;
  * probe the median run first, then the lower and upper quartile runs (when >= 3 remain);
  * a run that logs no radius is REPLACED by the remaining run nearest the median, and the
    radius pinned is the MEDIAN of the radii measured. Every probe is recorded.

BRUKER INDEX
------------
A .d is only as long as its analysis.tdf frame index says. 342 .d on HIVE have a truncated
index beside a complete analysis.tdf_bin: the copier left a stale mid-acquisition -wal beside
the finished tdf, and a READ-WRITE sqlite open checkpointed those pages into it. DIA-NN then
reads the run silently short -- in the pilot, 121 cycles and 80 precursors -- and a bytes rule
picked exactly such a run (tdf_bin full size, index 0.7%). So a .d is never probed when:

  * its analysis.tdf header is in WAL mode (bytes 18-19 = 2,2; every truncated file on HIVE is,
    intact ones are 1,1),
  * a non-empty analysis.tdf-wal or analysis.tdf-journal sits beside it (the next read-write
    open by anything would rewrite it), or
  * the last frame block (its TimsId offset + the uint32 block size stored there) ends before
    99.9% of analysis.tdf_bin, or past its end.

The header and side files are checked with plain byte reads, and the tdf is only opened with
`file:<tdf>?mode=ro&immutable=1` (tdf_uri). `mode=ro` alone still reads a stale WAL and leaves
-wal/-shm files beside it. Such runs are named in the job log: they are still SEARCHED by the
chain, and someone should look at them.

THERMO .RAW AND OTHER FILES
---------------------------
Ranked by file size only. Reading a .raw's acquisition length or checking its integrity needs
Thermo's RawFileReader (.NET), the same library DIA-NN loads, so it is not done here. A .raw
that DIA-NN cannot open logs no radius and is replaced by the next representative run, with its
log in window.json. The one failure that is NOT replaced is DIA-NN's missing-.NET error: no
other .raw can succeed in that environment, so probing stops at once and says how to fix it.

PROBES RUN ONE AT A TIME
------------------------
Measured on HIVE (16-CPU allocation, full mouse predicted library, the same three runs each
way): timsTOF .d took 1083 s one after another at 16 threads against 1314 s with three
concurrent probes at 5 threads each; Orbitrap .raw 289 s against 250 s. Neither order wins, but
concurrent probes need three probes' memory at once (~10 GB each on timsTOF) and triple the
simultaneous reads of the raw data, and one at a time is what lets a failed run be replaced.
Each probe is announced on stderr -- the job log -- and DIA-NN's output is copied there as it
arrives (and to the probe's probe.log): a silent 18-minute job log is past watch_run.sh's
15-minute stall rule, whose playbook is to cancel the job.

TIME
----
--timeout is per probe (default 3600 s): a run that hangs is cut there, so the next run still
gets its turn. --budget is for ALL probes together: no probe starts after it, and a running one
is cut at it, so the probes cannot outlive the job's wall clock and the evidence is always
written. The chain's step 1b passes its own wall clock less 10 minutes.

STOPPING DIA-NN
---------------
--diann is not always the engine itself: acquire_tools.sh records `apptainer exec ... <sif>
/diann-*/diann-linux` where there is no native build, and sites wrap binaries in scripts. So
DIA-NN runs in its own session (process group), writes to probe.log rather than a pipe a
grandchild could hold open, and every stop -- radius found (SIGTERM, then SIGKILL after 30 s),
timeout or budget (SIGKILL), SIGTERM/SIGINT/SIGHUP to this script -- goes to the whole group.
Measured on HIVE under apptainer 1.5.3: a 5 s timeout returned after 5.3 s with 0 container
processes left (a probe that signalled only its child returned after 41.2 s and left 2).

USAGE
-----
    probe_window.py --diann <diann cmd> --raw <run> [<run> ...] | --raw-list <file> \\
                    --fasta <f.fasta> --lib <predicted.speclib> [--threads 16] \\
                    [--timeout 3600] [--budget S] [--max-probes 3] [--max-failures 3] \\
                    [--workdir DIR] [--write-cfg CFG] [-- <DIA-NN flags of the real search>]

The library must already exist -- run this after step 1 of the chain, not before. Everything
after a bare `--` goes to DIA-NN verbatim (the chain passes the cfg flags this way, so the probe
runs under exactly the flags steps 2-5 run). --diann is exec'd without a shell, so the .NET 8
exports Thermo .raw needs must be in THIS script's environment (ensure_dotnet8.sh, next to
this script, prints DOTNET_ROOT). Prints JSON: the pinned radius and every probe's evidence --
also on failure, with window_radius null and a non-zero exit status.
"""
import argparse, json, os, re, shlex, signal, sqlite3, statistics, struct, subprocess, sys
import tempfile, time
from urllib.parse import quote

HERE = os.path.dirname(os.path.abspath(__file__))

# DIA-NN prints e.g. "Scan window radius set to 7". Match loosely (case-insensitive,
# tolerant of the leading '[m:ss]' timestamp) but require the integer.
WINDOW_RE = re.compile(r"window\s+radius\s+set\s+to\s+(\d+)", re.I)

# Printed by DIA-NN when it has no usable .NET 8 runtime -- and it still exits 0. Measured on
# HIVE with DIA-NN 2.7.0 (srun job 23509324). Not a property of the run: nothing else will read.
DOTNET_RE = re.compile(r"cannot read \.raw files", re.I)

# A run under this fraction of a typical run is treated as a blank, wash or failed injection and
# never probed. It only decides which runs are PROBED; every run is still searched.
#
# "Typical" is the median of the LARGER half of the runs, not of all of them: with a blank
# injected after every sample the blanks are the majority, the overall median is a blank's
# size, and a floor measured from it excludes nothing (review reproduction: 6 samples and 7
# blanks -> probes blank3, blank6 and s2). Measured sizes on HIVE, raw_data/Lumos1/noi25:
# washes 120-240 MB and 60-min washes 527-567 MB, beside 80 90-min DIA runs of 0.93-2.85 GB
# (median 2.12 GB; median of the larger half 2.42 GB). The 60-min washes are 25-27% of that DIA
# median, so a 25% floor let them in. At 50% of the larger half's median they are out; 2 of the
# 80 DIA runs (38-40%) are not probed. They are still searched. The rule keeps blanks out while
# they are at most 60-75% of the list; on a folder that is 86% washes it fails, so hand the
# search the cohort, not the folder.
MIN_FRACTION = 0.5

# A .d whose last indexed frame block ends before this fraction of analysis.tdf_bin is treated
# as a truncated index (the pilot's was 0.7%).
INDEX_COVERAGE_MIN = 0.999

SELECTION_RULE = (
    "never probe a Bruker .d whose analysis.tdf is in WAL mode, has a non-empty -wal/-journal "
    "beside it, or whose frame index ends before 99.9% of analysis.tdf_bin; drop runs under "
    f"{MIN_FRACTION:.0%} of the median size of the larger half of the runs (.d: indexed bytes of "
    "analysis.tdf_bin); when every run left is a .d with a frame index, rank by acquisition time "
    f"(max Frames.Time) and drop runs under {MIN_FRACTION:.0%} of its larger-half median too, "
    "otherwise rank by size; probe the median run, then the lower and upper quartile runs; "
    "replace a run that logs no radius with the remaining run nearest the median; pin the "
    "median of the measured radii")

POLL_S = 0.2          # how often DIA-NN's log is read; nothing against minutes of calibration
STOP_GRACE_S = 30     # SIGTERM -> SIGKILL, once a radius is in hand
_HARD = getattr(signal, "SIGKILL", signal.SIGTERM)

_LIVE = set()    # running DIA-NN processes, so a signal to the probe can stop them
_STOP = []       # signals received; once set, no probe starts and running ones are stopped


# ------------------------------------------------------------------------------------------
# Bruker .d integrity
# ------------------------------------------------------------------------------------------
def tdf_uri(tdf):
    """The ONLY way this script opens an analysis.tdf. immutable=1: SQLite neither reads a -wal
    nor creates -wal/-shm, and cannot write. `mode=ro` alone still reads (and depends on) a
    stale WAL -- the state that truncated 342 tdfs on HIVE when something opened them
    read-write. Percent-encoded, so a path with a space, `#` or `?` is still that path."""
    return "file:" + quote(os.path.abspath(tdf)) + "?mode=ro&immutable=1"


def bruker_tdf_status(dotd):
    """What a Bruker TDF .d's frame index says, and whether it can be trusted.

    Returns None when `dotd` has neither analysis.tdf nor analysis.tdf_bin (not a TDF run).
    Otherwise {time_s, n_frames, indexed_bytes, tdf_bin_bytes, index_coverage, problems}:
    time_s is max(Frames.Time) in seconds, indexed_bytes where the last indexed frame block
    ends in analysis.tdf_bin, and `problems` is empty only for a run safe to probe. The tdf is
    not opened at all when its header or side files already disqualify it."""
    d = dotd.rstrip("/") or dotd
    tdf, tdf_bin = os.path.join(d, "analysis.tdf"), os.path.join(d, "analysis.tdf_bin")
    has_tdf, has_bin = os.path.isfile(tdf), os.path.isfile(tdf_bin)
    if not (has_tdf or has_bin):
        return None
    st = {"time_s": None, "n_frames": None, "indexed_bytes": None,
          "tdf_bin_bytes": os.path.getsize(tdf_bin) if has_bin else None,
          "index_coverage": None, "problems": []}
    problems = st["problems"]
    if not has_tdf:
        problems.append("analysis.tdf_bin without analysis.tdf (no frame index to read it by)")
        return st
    if not has_bin:
        problems.append("analysis.tdf without analysis.tdf_bin (no spectra)")

    # Plain bytes, never SQLite: nothing here may touch the file's journal state.
    with open(tdf, "rb") as fh:
        head = fh.read(100)
    if len(head) < 100 or not head.startswith(b"SQLite format 3\x00"):
        problems.append("analysis.tdf has no SQLite header")
        return st
    if head[18] == 2 or head[19] == 2:
        problems.append(f"analysis.tdf is in WAL mode (header bytes 18-19 = {head[18]},{head[19]}; "
                        "every truncated tdf found on HIVE is, intact ones are 1,1)")
    for side in ("-wal", "-journal"):
        path = tdf + side
        if os.path.isfile(path) and os.path.getsize(path) > 0:
            problems.append(f"non-empty analysis.tdf{side} ({os.path.getsize(path)} bytes) beside "
                            "it -- a read-write open would rewrite (and may truncate) the index")
    if problems:
        return st

    try:
        con = sqlite3.connect(tdf_uri(tdf), uri=True)
        try:
            t_max, tims_max, n = con.execute(
                "SELECT MAX(Time), MAX(TimsId), COUNT(*) FROM Frames").fetchone()
        finally:
            con.close()
    except sqlite3.Error as e:
        problems.append(f"analysis.tdf frame index cannot be read ({e})")
        return st
    if not n or t_max is None or tims_max is None:
        problems.append("analysis.tdf indexes no frames")
        return st
    st["time_s"], st["n_frames"] = float(t_max), int(n)

    size = st["tdf_bin_bytes"]
    with open(tdf_bin, "rb") as fh:
        fh.seek(int(tims_max))
        raw = fh.read(4)
    if len(raw) < 4:
        problems.append(f"the last indexed frame starts at byte {tims_max}, past the end of "
                        f"analysis.tdf_bin ({size} bytes) -- tdf_bin is truncated")
        return st
    end = int(tims_max) + struct.unpack("<I", raw)[0]
    st["indexed_bytes"] = end
    st["index_coverage"] = end / size if size else 0.0
    if end > size:
        problems.append(f"the last indexed frame block ends at byte {end}, past the end of "
                        f"analysis.tdf_bin ({size} bytes) -- tdf_bin is truncated")
    elif end < INDEX_COVERAGE_MIN * size:
        problems.append(f"the frame index covers only {st['index_coverage']:.1%} of "
                        f"analysis.tdf_bin ({n} frames, {t_max / 60:.1f} min) -- a truncated "
                        "index; DIA-NN would read only that part of the run")
    return st


# ------------------------------------------------------------------------------------------
# Which runs
# ------------------------------------------------------------------------------------------
def measure_run(path):
    """{file, readable, size_bytes, time_s, problems} for one input run.

    A Bruker TDF .d: time_s from its frame index, size_bytes the indexed part of
    analysis.tdf_bin, problems from bruker_tdf_status(). Any other directory: the bytes of its
    top-level files. A file (.raw, .mzML, ...): its size, time_s None."""
    p = path.rstrip("/") or path
    m = {"file": path, "readable": False, "size_bytes": None, "time_s": None, "problems": []}
    try:
        if os.path.isdir(p):
            st = bruker_tdf_status(p)
            if st is None:
                m["size_bytes"] = sum(os.path.getsize(os.path.join(p, x)) for x in os.listdir(p)
                                      if os.path.isfile(os.path.join(p, x)))
            else:
                m["time_s"], m["problems"] = st["time_s"], list(st["problems"])
                m["size_bytes"] = (st["indexed_bytes"] if st["indexed_bytes"] is not None
                                   else st["tdf_bin_bytes"] or 0)
                m["index_coverage"] = st["index_coverage"]
            m["readable"] = True
        elif os.path.isfile(p):
            m["size_bytes"], m["readable"] = os.path.getsize(p), True
    except OSError as e:
        m["readable"], m["error"] = False, str(e)
    return m


class NoProbeableRun(ValueError):
    """No input run can be probed. `selection` holds what was found, for the evidence JSON."""

    def __init__(self, msg, selection):
        super().__init__(msg)
        self.selection = selection


def _typical(values):
    """Median of the larger half -- see MIN_FRACTION."""
    v = sorted(values)
    return statistics.median(v[len(v) // 2:])


def select_representative(paths, max_probes=3, min_fraction=MIN_FRACTION):
    """Choose the runs to probe. Deterministic and independent of input order (ties are broken
    by path), because "which run is first" is exactly the dependence being removed.

    Returns {chosen, reserves, rank_by, reference_size_bytes, min_size_bytes,
    reference_time_s, min_time_s, excluded_small, excluded_damaged, unreadable, n_inputs,
    n_eligible, rule}. `chosen` is in probing order: median, lower quartile, upper quartile
    (positions n//2, floor(m/4), ceil(3m/4) of the m+1 ranked eligible runs -- distinct for every
    n >= 3; on an even count the upper-middle run, the less likely to be a weak injection).
    `reserves` are the other eligible runs, nearest the median first (the larger on a tie).
    Raises NoProbeableRun (a ValueError) when nothing is eligible."""
    measured = [measure_run(p) for p in paths]
    usable = [m for m in measured if m["readable"] and not m["problems"]]
    sel = {"chosen": [], "reserves": [], "rank_by": None,
           "reference_size_bytes": None, "min_size_bytes": None,
           "reference_time_s": None, "min_time_s": None,
           "excluded_small": [],
           "excluded_damaged": [{"file": m["file"], "problems": m["problems"]}
                                for m in measured if m["readable"] and m["problems"]],
           "unreadable": [m["file"] for m in measured if not m["readable"]],
           "n_inputs": len(paths), "n_eligible": 0, "rule": SELECTION_RULE}
    if not usable:
        raise NoProbeableRun(
            f"none of the {len(paths)} input runs can be probed ({len(sel['unreadable'])} "
            f"unreadable, {len(sel['excluded_damaged'])} with a damaged Bruker index), so there "
            "is nothing to measure the scan window on: " + ", ".join(paths[:5]), sel)

    def small(m, why):
        sel["excluded_small"].append({"file": m["file"], "size_bytes": m["size_bytes"],
                                      "time_s": m["time_s"], "why": why})

    # The size floor first, over every usable run: it is what removes a failed acquisition with
    # no data at all -- a .d with neither analysis.tdf nor tdf_bin has no acquisition time, and
    # deciding the ranking before it was gone made a whole timsTOF cohort rank by bytes (HIVE,
    # garg cohort). Then, when every run left has a frame index, the time floor and time ranking.
    ref_size = _typical([m["size_bytes"] for m in usable])
    sel["reference_size_bytes"], sel["min_size_bytes"] = ref_size, min_fraction * ref_size
    eligible = []
    for m in usable:
        if m["size_bytes"] < sel["min_size_bytes"]:
            small(m, f"size {m['size_bytes'] / 1e9:.2f} GB is under {min_fraction:.0%} of a "
                     f"typical run's {ref_size / 1e9:.2f} GB")
        else:
            eligible.append(m)
    by_time = bool(eligible) and all(m["time_s"] is not None for m in eligible)
    sel["rank_by"] = "time" if by_time else "size"
    if by_time:
        ref_time = _typical([m["time_s"] for m in eligible])
        sel["reference_time_s"], sel["min_time_s"] = ref_time, min_fraction * ref_time
        timed, eligible = eligible, []
        for m in timed:
            if m["time_s"] < sel["min_time_s"]:
                small(m, f"acquisition time {m['time_s'] / 60:.1f} min is under "
                         f"{min_fraction:.0%} of a typical run's {ref_time / 60:.1f} min")
            else:
                eligible.append(m)
    if not eligible:
        raise NoProbeableRun("no run clears the size and acquisition-time floors", sel)

    if by_time:
        eligible.sort(key=lambda m: (m["time_s"], m["size_bytes"], m["file"]))
    else:
        eligible.sort(key=lambda m: (m["size_bytes"], m["file"]))
    n, last = len(eligible), len(eligible) - 1
    mid = n // 2
    if max_probes >= 3 and n >= 3:
        picks = [("median", mid), ("lower_quartile", last // 4),
                 ("upper_quartile", -(-3 * last // 4))]
    else:
        picks = [("median", mid)]
    taken = {i for _, i in picks}
    rest = sorted((i for i in range(n) if i not in taken), key=lambda i: (abs(i - mid), -i))

    def entry(m, role):
        return {"file": m["file"], "size_bytes": m["size_bytes"], "time_s": m["time_s"],
                "role": role}

    sel["chosen"] = [entry(eligible[i], role) for role, i in picks]
    sel["reserves"] = [entry(eligible[i], "reserve") for i in rest]
    sel["n_eligible"] = n
    return sel


# ------------------------------------------------------------------------------------------
# Running DIA-NN
# ------------------------------------------------------------------------------------------
def _signal_group(p, sig):
    """Send `sig` to DIA-NN and everything it started (its process group). Where there are no
    process groups (Windows) only the direct child can be reached."""
    try:
        if hasattr(os, "killpg"):
            os.killpg(p.pid, sig)
        elif p.poll() is None:
            p.kill()
    except OSError:          # the whole group is already gone
        pass


def _group_alive(pgid):
    try:
        os.killpg(pgid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def _end_group(p, grace=STOP_GRACE_S):
    """SIGTERM DIA-NN's whole group, then SIGKILL whatever is left after `grace` s. A wrapper
    that forks (a bash script without `exec`; `apptainer exec`) leaves the engine as a
    grandchild, which signalling only the direct child never reaches. A process that moves
    itself into a NEW session escapes this."""
    _signal_group(p, signal.SIGTERM)
    if not hasattr(os, "killpg"):
        try:
            p.wait(timeout=grace)
        except subprocess.TimeoutExpired:
            p.kill()
        return
    deadline = time.time() + grace
    while time.time() < deadline:
        p.poll()                           # reap our child so it stops counting as a member
        if not _group_alive(p.pid):
            return
        time.sleep(0.1)
    _signal_group(p, _HARD)
    p.wait()


def _kill_now(p):
    """Timeout, budget or a signal: the whole group, at once."""
    _signal_group(p, _HARD)
    try:
        p.wait(timeout=STOP_GRACE_S)
    except subprocess.TimeoutExpired:
        pass


def _on_signal(signum, _frame):
    """SIGTERM (scancel, a SLURM time limit), SIGINT, SIGHUP. DIA-NN is in its own process group,
    so a signal to the job's group does not reach it; stop it here. The probes then return and
    main() exits with 128 + signum."""
    _STOP.append(signum)
    for p in list(_LIVE):
        _signal_group(p, _HARD)


def _say(text):
    """One line into the job log (stderr), now -- not when the job ends."""
    sys.stderr.write(text + "\n")
    sys.stderr.flush()


def run_probe(diann, raw, fasta, lib, threads, timeout, extra="", extra_args=(), workdir=None,
              tag=""):
    """Run DIA-NN on one run until it logs the radius.

    Returns {radius, lines, timed_out, log, environmental}; `environmental` means the failure
    is not about this run (DIA-NN cannot read .raw here, or could not be started).

    `workdir` gives DIA-NN its own --temp and --out and holds probe.log. Without them a probe
    that never logs a radius runs on, writes its report into the working directory (in the
    chain, <out> -- the path step 5's check reads) and, if the search completes, its .quant
    NEXT TO THE RAW FILE. They go before --threads, so the flags of the real search stay the
    tail of DIA-NN's argv, as in steps 2-5."""
    workdir = workdir or tempfile.mkdtemp(prefix="probe_window_")
    os.makedirs(workdir, exist_ok=True)
    flags = shlex.split(extra) + list(extra_args)
    cmd = shlex.split(diann) + ["--f", raw, "--fasta", fasta, "--lib", lib]
    if "--temp" not in flags:
        cmd += ["--temp", workdir]
    if "--out" not in flags:
        cmd += ["--out", os.path.join(workdir, "report.parquet")]
    cmd += ["--threads", str(threads)] + flags
    log_path = os.path.join(workdir, "probe.log")
    res = {"radius": None, "lines": [], "timed_out": False, "log": log_path,
           "environmental": False}
    deadline = time.time() + max(0.0, timeout)

    def note(msg):
        res["lines"].append(msg)
        _say(tag + msg)
        with open(log_path, "a") as fh:
            fh.write(msg + "\n")

    if _STOP or time.time() >= deadline:
        open(log_path, "w").close()
        res["timed_out"] = not _STOP
        note("[probe_window] stopped by a signal before this probe started" if _STOP else
             "[probe_window] timeout: no time left to start this probe")
        return res

    try:
        with open(log_path, "wb") as sink:
            p = subprocess.Popen(cmd, stdout=sink, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL, start_new_session=True)
    except OSError as e:
        res["environmental"] = True
        note(f"[probe_window] could not start DIA-NN ({cmd[0]}): {e}")
        return res
    _LIVE.add(p)
    if _STOP:                       # a signal arrived while it was starting
        _signal_group(p, _HARD)
    try:
        with open(log_path, "rb") as src:
            pending = b""
            while res["radius"] is None:
                exited = p.poll() is not None   # before reading: the read then has it all
                pending += src.read()
                if exited and pending and not pending.endswith(b"\n"):
                    pending += b"\n"              # DIA-NN is gone: its last line is complete
                *complete, pending = pending.split(b"\n")
                for raw_line in complete:
                    line = raw_line.decode("utf-8", "replace").rstrip()
                    res["lines"].append(line)
                    _say(tag + line)
                    m = WINDOW_RE.search(line)
                    if m and int(m.group(1)) > 0:
                        res["radius"] = int(m.group(1))
                        break              # got it -- no need to finish the search
                if res["radius"] is not None or exited or _STOP:
                    break
                if time.time() >= deadline:
                    res["timed_out"] = True
                    break
                time.sleep(POLL_S)
            if res["radius"] is None and pending.strip():
                res["lines"].append(pending.decode("utf-8", "replace").rstrip())
                _say(tag + res["lines"][-1])
    finally:
        if p.poll() is None and not (res["timed_out"] or _STOP):
            _end_group(p)                   # radius in hand: let it stop gracefully
        else:
            _kill_now(p)                    # timeout, signal, or the leader already gone
        _LIVE.discard(p)
    if res["timed_out"]:
        note(f"[probe_window] timeout: DIA-NN stopped after {timeout:.0f} s without logging a "
             "scan-window radius")
    res["environmental"] = res["environmental"] or any(DOTNET_RE.search(ln) for ln in res["lines"])
    return res


def probe(diann, raw, fasta, lib, threads, timeout, extra="", extra_args=(), workdir=None,
          tag=""):
    """(radius | None, log lines) for one run -- see run_probe()."""
    r = run_probe(diann, raw, fasta, lib, threads, timeout, extra, extra_args, workdir, tag)
    return r["radius"], r["lines"]


# ------------------------------------------------------------------------------------------
def _gb(b):
    return "?" if b is None else f"{b / 1e9:.2f} GB"


def _min(t):
    return "" if t is None else f", {t / 60:.1f} min"


def _base(path):
    return os.path.basename(path.rstrip("/"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diann", required=True,
                    help="DIA-NN command (binary path or 'apptainer exec ... diann-linux')")
    ap.add_argument("--raw", nargs="+", default=[],
                    help="the cohort's runs (all of them -- representative runs are chosen here)")
    ap.add_argument("--raw-list", help="file with one run path per line (handles spaces in paths)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--lib", required=True, help="predicted spectral library (step 1 output)")
    ap.add_argument("--threads", type=int, default=16, help="DIA-NN threads per probe")
    ap.add_argument("--timeout", type=int, default=3600,
                    help="seconds for ONE probe (default 3600); a run that hangs is cut here "
                         "and replaced")
    ap.add_argument("--budget", type=int, default=None,
                    help="seconds for ALL probes together (default: no limit). The chain passes "
                         "its step's wall clock less a margin, so the evidence is always written")
    ap.add_argument("--max-probes", type=int, choices=(1, 3), default=3,
                    help="radii to measure: 3 (default) = median + quartile runs; 1 = the median run")
    ap.add_argument("--max-failures", type=int, default=3,
                    help="stop after this many runs logged no radius (default 3)")
    ap.add_argument("--workdir", help="per-probe DIA-NN --temp/--out and logs "
                    "(default: a fresh temporary directory)")
    ap.add_argument("--extra", default="", help="extra DIA-NN flags to match the real search, "
                    "as ONE shlex-quoted string (for hand use; the chain passes them after --)")
    ap.add_argument("--write-cfg", help="append '--window N' to this cfg file on success")
    # Everything after a bare `--` goes to DIA-NN verbatim, as separate arguments, so bash quotes
    # and expands them exactly as it does for steps 2-5.
    argv = sys.argv[1:]
    extra_args = []
    if "--" in argv:
        k = argv.index("--")
        argv, extra_args = argv[:k], argv[k + 1:]
    a = ap.parse_args(argv)

    started = time.time()
    deadline = started + a.budget if a.budget is not None else None
    raws = list(a.raw)
    if a.raw_list:
        with open(a.raw_list) as fh:
            raws += [ln.strip() for ln in fh if ln.strip()]
    if not raws:
        sys.exit("no runs given (--raw or --raw-list)")
    for f, what in ((a.fasta, "fasta"), (a.lib, "library")):
        if not os.path.exists(f):
            sys.exit(f"{what} not found: {f}")

    def report(result):
        print(json.dumps(result, indent=2))
        sys.stdout.flush()

    try:
        sel = select_representative(raws, max_probes=a.max_probes)
    except NoProbeableRun as e:
        report({"window_radius": None, "pin_as": None, "radii": [], "radii_agree": None,
                "incomplete": None, "failed": [], "stopped_because": "no_probeable_run",
                "probes": [], "selection": e.selection})
        sys.exit(str(e))

    for f in sel["unreadable"]:
        _say(f"[probe_window] WARNING: cannot read {f}; not considered")
    for x in sel["excluded_damaged"]:
        _say(f"[probe_window] WARNING: not probing {_base(x['file'])}: "
             + "; ".join(x["problems"]) + " -- it is still SEARCHED by the chain; check it")
    for x in sel["excluded_small"]:
        _say(f"[probe_window] not probing {_base(x['file'])}: {x['why']} -- likely a blank, "
             "wash or failed injection")
    _say(f"[probe_window] {sel['n_eligible']} of {sel['n_inputs']} runs eligible, ranked by "
         + ("acquisition time" if sel["rank_by"] == "time" else "size")
         + f"; measuring {len(sel['chosen'])}, {len(sel['reserves'])} in reserve")

    for name in ("SIGTERM", "SIGINT", "SIGHUP"):
        if hasattr(signal, name):
            signal.signal(getattr(signal, name), _on_signal)

    workdir = a.workdir or tempfile.mkdtemp(prefix="probe_window_")
    targets, reserves = sel["chosen"], list(sel["reserves"])
    queue = [dict(c) for c in targets]
    probes, tails, failures, stopped, k = [], {}, 0, None, 0
    while queue:
        if _STOP:
            stopped = "signal"
            break
        if deadline is not None and time.time() >= deadline:
            stopped = "budget"
            _say(f"[probe_window] budget: {a.budget} s spent; no further probe started")
            break
        c = queue.pop(0)
        k += 1
        base = _base(c["file"])
        if c["role"] == "reserve":
            label = f"probe {k}: {base} (reserve for {_base(c['replaces'])}"
        else:
            label = f"probe {k}/{len(targets)}: {base} ({c['role']}"
        _say(f"[probe_window] {label}, {_gb(c['size_bytes'])}{_min(c['time_s'])}), "
             f"{a.threads} threads")
        limit = a.timeout
        budget_cut = deadline is not None and deadline - time.time() < a.timeout
        if budget_cut:
            limit = max(0.0, deadline - time.time())
        t0 = time.time()
        wd = os.path.join(workdir, f"probe{k}_{base}")
        r = run_probe(a.diann, c["file"], a.fasta, a.lib, a.threads, limit, a.extra,
                      extra_args, workdir=wd, tag=f"[probe {k}] ")
        secs = round(time.time() - t0, 1)
        probes.append(dict(c, radius=r["radius"], seconds=secs, threads=a.threads,
                           timed_out=r["timed_out"], log=r["log"]))
        if r["radius"] is not None:
            _say(f"[probe_window] probe {k}: {base} -> radius {r['radius']} in {secs} s")
            continue
        tails[k] = r["lines"]
        failures += 1
        _say(f"[probe_window] probe {k}: {base} -> NO radius"
             + (" (timed out)" if r["timed_out"] else "") + f" in {secs} s")
        if _STOP:
            stopped = "signal"
            break
        if r["environmental"]:
            stopped = "environment"
            break
        if r["timed_out"] and budget_cut:
            stopped = "budget"
            _say(f"[probe_window] budget: {a.budget} s spent during this probe; no further "
                 "probe started")
            break
        if failures >= a.max_failures:
            stopped = "max_failures"
            break
        if reserves:
            nxt = dict(reserves.pop(0), replaces=c["file"])
            queue.append(nxt)
            _say(f"[probe_window] {base} gave no radius; {_base(nxt['file'])}, the next run "
                 "nearest the median, takes its place")

    radii = [p["radius"] for p in probes if p["radius"] is not None]
    if stopped is None:
        stopped = "measured" if len(radii) >= len(targets) else "no_more_runs"
    pin = bool(radii) and stopped not in ("environment", "signal")
    radius = statistics.median_low(radii) if pin else None
    result = {
        "window_radius": radius,
        "pin_as": f"--window {radius}" if radius is not None else None,
        "radii": radii,
        "radii_agree": (len(set(radii)) == 1) if radii else None,
        "incomplete": (len(radii) < len(targets)) if radius is not None else None,
        "failed": [_base(p["file"]) for p in probes if p["radius"] is None],
        "stopped_because": stopped,
        "seconds": round(time.time() - started, 1),
        "probes": probes,
        "selection": dict({x: v for x, v in sel.items() if x not in ("chosen", "reserves")},
                          planned=[c["file"] for c in targets],
                          reserves=[c["file"] for c in sel["reserves"]]),
        "note": ("Pin this in the cfg for EVERY step of the parallel chain. It is a property of "
                 "the acquisition scheme, so it is valid for all files acquired with the same "
                 "method -- but re-probe for a different gradient, cycle time, or instrument."),
    }

    _say("[probe_window] per-run radii:")
    for p in probes:
        _say(f"[probe_window]   {p['role']:<15} {_base(p['file'])}  {_gb(p['size_bytes'])}"
             f"{_min(p['time_s'])}  "
             + (f"radius {p['radius']}" if p["radius"] is not None else
                "NO radius" + (" (timed out)" if p["timed_out"] else "")))

    if stopped == "signal":
        report(result)
        sys.exit(128 + _STOP[0])

    if radius is None:
        report(result)                     # the evidence first: window.json is what is read
        for i, p in enumerate(probes, 1):
            if p["radius"] is None:
                _say(f"\n--- {p['file']}: no scan-window radius; DIA-NN log tail "
                     f"(full log {p['log']}) ---\n" + "\n".join(tails.get(i, [])[-25:]))
        if any(DOTNET_RE.search(ln) for t in tails.values() for ln in t):
            _say("\nDIA-NN could not open Thermo .raw: there is no .NET 8 runtime in this "
                 "environment. Export DOTNET_ROOT before running this probe:\n"
                 f'    export DOTNET_ROOT="$(bash {shlex.quote(os.path.join(HERE, "ensure_dotnet8.sh"))} '
                 '| tail -1)"; export PATH="$DOTNET_ROOT:$PATH"\n'
                 "(diann_parallel.py puts that export into every generated step.)")
        why = {"environment": "a failure no other run can fix -- DIA-NN could not start, or "
                              "cannot read .raw in this environment (above) -- so no other run "
                              "was tried",
               "max_failures": f"{failures} runs logged no radius",
               "budget": f"the {a.budget} s budget ran out",
               "no_more_runs": "every eligible run was tried"}.get(stopped, stopped)
        sys.exit("Could not read the scan-window radius: " + why + " ("
                 + ", ".join(result["failed"]) + "). Do NOT guess a value -- an inconsistent "
                 "window across files is exactly the defect this probe exists to prevent. Fix the "
                 "cause above, or re-run with --raw naming runs to use.")

    if a.write_cfg:
        with open(a.write_cfg, "a") as fh:
            fh.write(f"\n--window {radius}\n")
    report(result)
    if result["incomplete"]:
        _say(f"[probe_window] WARNING: pinned --window {radius} from {len(radii)} of the "
             f"{len(targets)} runs planned ({stopped}; no radius from "
             + ", ".join(result["failed"]) + "). It is the median of the runs that answered; "
             "the per-run evidence is in the JSON.")


if __name__ == "__main__":
    main()
