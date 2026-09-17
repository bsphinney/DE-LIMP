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
width), not of the individual sample, so one file is enough to measure it -- but it
must be MEASURED, not guessed, because it depends on the gradient and method.

WHAT IT DOES
------------
Runs DIA-NN on one file, watches the log for

    Scan window radius set to N

and terminates as soon as that line appears (it is printed during calibration, long
before the search finishes, so this costs a couple of minutes rather than a full pass).
Prints JSON with the radius; pin it as `--window N` in the cfg for every step.

USAGE
-----
    probe_window.py --diann <diann cmd> --raw <one file> --fasta <f.fasta> \
                    --lib <predicted.speclib> [--threads 16] [--timeout 3600]

The library must already exist -- run this after step 1 of the chain, not before.

--diann is exec'd directly (no shell), so an environment it needs -- e.g. the .NET 8 exports
DIA-NN 2.6 needs to read Thermo .raw -- must be exported BEFORE running this script; it cannot
be spliced into --diann. The chain's step 1b does exactly that.
"""
import argparse, json, os, re, shlex, signal, subprocess, sys, threading, time

# DIA-NN prints e.g. "Scan window radius set to 7". Match loosely (case-insensitive,
# tolerant of the leading '[m:ss]' timestamp) but require the integer.
WINDOW_RE = re.compile(r"window\s+radius\s+set\s+to\s+(\d+)", re.I)


def _signal_group(pgid, sig):
    try:
        os.killpg(pgid, sig)
    except (ProcessLookupError, PermissionError):
        pass                               # the whole group is already gone


def _group_alive(pgid):
    try:
        os.killpg(pgid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def _end_group(p, grace=30):
    """SIGTERM DIA-NN's whole process group, then SIGKILL whatever is still there after `grace` s.

    Signalling only the direct child is not enough: when --diann is a wrapper that forks (a bash
    script without `exec`; `apptainer exec` very likely too), the grandchild that is DIA-NN keeps
    running -- holding the stdout pipe open, so the read never ends and the next file is never
    tried, and burning the job's CPUs as an orphan. Verified with a forking bash wrapper. A
    process that moves itself into a NEW session escapes this; whether apptainer's starter does
    is unverified."""
    _signal_group(p.pid, signal.SIGTERM)
    deadline = time.time() + grace
    while time.time() < deadline:
        p.poll()                           # reap our child so it stops counting as a member
        if not _group_alive(p.pid):
            return
        time.sleep(0.1)
    _signal_group(p.pid, signal.SIGKILL)
    p.wait()


def probe(diann, raw, fasta, lib, threads, timeout, extra="", extra_args=()):
    cmd = shlex.split(diann) + [
        "--f", raw, "--fasta", fasta, "--lib", lib,
        "--threads", str(threads),
    ] + shlex.split(extra) + list(extra_args)
    # Its own session, so DIA-NN and anything it forks share one process group we can end.
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1, start_new_session=True)
    radius, lines, deadline = None, [], time.time() + timeout
    # The deadline check below only runs when DIA-NN prints a line, so a DIA-NN that goes
    # silent would block this loop until SLURM's wall clock killed the whole job -- and step 1b
    # would never get to try its next file. The timer ends the read from outside, for the whole
    # group (see _end_group).
    watchdog = threading.Timer(timeout, _signal_group, (p.pid, signal.SIGKILL))
    watchdog.daemon = True
    watchdog.start()
    try:
        for line in p.stdout:
            lines.append(line.rstrip())
            m = WINDOW_RE.search(line)
            if m and int(m.group(1)) > 0:
                radius = int(m.group(1))
                break                      # got it -- no need to finish the search
            if time.time() > deadline:
                break
    finally:
        watchdog.cancel()
        _end_group(p)
    return radius, lines


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diann", required=True, help="DIA-NN command (binary path or 'apptainer exec ... diann-linux')")
    ap.add_argument("--raw", required=True, help="ONE raw/mzML file representative of the acquisition")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--lib", required=True, help="predicted spectral library (step 1 output)")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--timeout", type=int, default=3600,
                    help="give up after N seconds (default 3600)")
    ap.add_argument("--extra", default="", help="extra DIA-NN flags to match the real search, "
                    "as ONE shlex-quoted string (for hand use)")
    ap.add_argument("--write-cfg", help="append '--window N' to this cfg file on success")
    # Everything after a bare `--` goes to DIA-NN verbatim, as separate arguments. Step 1b
    # passes the cfg flags this way, so bash quotes and expands them exactly as it does for
    # steps 2-5 -- a second parse through shlex (--extra) differs from bash on `$VAR` and
    # backslashes, and the probe would measure under different flags than the search runs.
    argv = sys.argv[1:]
    extra_args = []
    if "--" in argv:
        k = argv.index("--")
        argv, extra_args = argv[:k], argv[k + 1:]
    a = ap.parse_args(argv)

    for f, what in ((a.raw, "raw"), (a.fasta, "fasta"), (a.lib, "library")):
        if not os.path.exists(f):
            sys.exit(f"{what} not found: {f}")

    radius, lines = probe(a.diann, a.raw, a.fasta, a.lib, a.threads, a.timeout, a.extra,
                          extra_args)
    if radius is None:
        sys.stderr.write("\n".join(lines[-25:]) + "\n")
        sys.exit("Could not read the scan-window radius from DIA-NN's output (see log tail "
                 "above). Do NOT guess a value -- an inconsistent window across files is "
                 "exactly the defect this probe exists to prevent.")

    if a.write_cfg:
        with open(a.write_cfg, "a") as fh:
            fh.write(f"\n--window {radius}\n")

    print(json.dumps({
        "window_radius": radius,
        "probed_file": os.path.basename(a.raw),
        "pin_as": f"--window {radius}",
        "note": ("Pin this in the cfg for EVERY step of the parallel chain. It is a "
                 "property of the acquisition scheme, so it is valid for all files "
                 "acquired with the same method -- but re-probe for a different "
                 "gradient, cycle time, or instrument."),
    }, indent=2))


if __name__ == "__main__":
    main()
