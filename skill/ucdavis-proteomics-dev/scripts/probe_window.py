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
"""
import argparse, json, os, re, shlex, subprocess, sys, threading, time

# DIA-NN prints e.g. "Scan window radius set to 7". Match loosely (case-insensitive,
# tolerant of the leading '[m:ss]' timestamp) but require the integer.
WINDOW_RE = re.compile(r"window\s+radius\s+set\s+to\s+(\d+)", re.I)


def probe(diann, raw, fasta, lib, threads, timeout, extra=""):
    cmd = shlex.split(diann) + [
        "--f", raw, "--fasta", fasta, "--lib", lib,
        "--threads", str(threads),
    ] + shlex.split(extra)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    radius, lines, deadline = None, [], time.time() + timeout
    try:
        for line in p.stdout:
            lines.append(line.rstrip())
            m = WINDOW_RE.search(line)
            if m:
                radius = int(m.group(1))
                break                      # got it -- no need to finish the search
            if time.time() > deadline:
                break
    finally:
        if p.poll() is None:
            p.terminate()
            try:
                p.wait(timeout=30)
            except subprocess.TimeoutExpired:
                p.kill()
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
    ap.add_argument("--extra", default="", help="extra DIA-NN flags to match the real search")
    ap.add_argument("--write-cfg", help="append '--window N' to this cfg file on success")
    a = ap.parse_args()

    for f, what in ((a.raw, "raw"), (a.fasta, "fasta"), (a.lib, "library")):
        if not os.path.exists(f):
            sys.exit(f"{what} not found: {f}")

    radius, lines = probe(a.diann, a.raw, a.fasta, a.lib, a.threads, a.timeout, a.extra)
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
