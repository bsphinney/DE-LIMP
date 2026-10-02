#!/usr/bin/env python3
"""
probe_fallback.py  --  What the search runs with when its pre-search probe measured nothing.

The chain's step 1b (and the single-shot search's probe before the search) runs probe_window.py
to measure the scan-window radius and, where planned, the mass accuracy. A probe that measured
nothing used to fail step 1b, and afterok then left steps 2-5 DependencyNeverSatisfied: the whole
search died over one measurement. Between 2026-09-29 17:00 and 09-30 08:00, 62 of fran-5b's 261
step-1b jobs (24%) did exactly that -- the probe's log tail hit ESTALE on Flinders NFS while
DIA-NN itself succeeded -- leaving 183 downstream jobs DependencyNeverSatisfied.

So a probe whose OWN MACHINERY failed no longer blocks the search: a log it could not read, a
crash, a time limit (probe_window.RETRY_ON / FALLBACK_ON; diann_parallel.probe_attempts). The job
retries a crash or an unreadable log once, and then this script sets what the search runs with
instead:

  scan window    left to DIA-NN: window.txt says `auto`, and steps 2-5 pass no --window, so
                 DIA-NN chooses the radius itself, per run
  mass accuracy  a level with a DIA-NN table value as documented (--ms1-ppm / --ms2-ppm, exactly
                 what the probe was given), the other level at the facility SOP
                 (estimate_params.SOP_MASS_ACC) -- the value a measured level is floored to
                 anyway -- tagged DEFAULT, not user-confirmed

Neither is presented as measured (DE-LIMP rule 2). window.txt says `auto`; <out>/probe_fallback.json
holds the record; search_provenance.json gets `probe_fallback`, and its scan_window (and, when mass
accuracy was to be measured, mass_acc and result.mass_acc) say fallback with the reason, with the
stable `mode` values fallback_auto / fallback_default (diann_parallel.SCAN_WINDOW_MODES /
MASS_ACC_MODES); the Methods say DIA-NN set the window automatically because the measurement
failed, and AUDIT.md and the run record carry caution() as a WARNING.

Nothing else is ever passed here -- the job fails as before: a mass accuracy REFUSED as
implausible (probe_window.EXIT_REFUSED; any run's logged value is checked, complete or not), the
environment (no .NET, DIA-NN cannot start), the probe's own arguments, runs DIA-NN finished
without logging it (the data answering), a signal. No fallback can fix those, and one would hide
them.

Usage:
  probe_fallback.py --measure window [mass-acc] [--ms1-ppm X | --ms2-ppm Y] --exit-code N
      [--attempts N] [--evidence window.json] [--window-file window.txt]
      [--massacc-file massacc.txt] [--write-cfg params.resolved.cfg.tmp]
      [--provenance search_provenance.json] [--out probe_fallback.json]
Exit 0 once everything is written; 1 when a file the search needs could not be.
"""
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from estimate_params import SOP_MASS_ACC, SOP_MASS_ACC_SOURCE  # noqa: E402  the one SOP
from probe_window import MEASURES, ppm_text, EXIT_CRASH  # noqa: E402
# massacc.txt in exactly the shape steps 2-5's guard (needs_measured) accepts; the stable modes
from diann_parallel import (MEASURED_FILE_RE, SCAN_WINDOW_MODES, MASS_ACC_MODES,  # noqa: E402
                            FALLBACK_RECORD)

# The flags the search passes per level, and the level each documented value belongs to.
LEVELS = (("--mass-acc", "ms2_ppm", "MS2"), ("--mass-acc-ms1", "ms1_ppm", "MS1"))

WINDOW_SOURCE = ("fallback (probe failed: {reason}) -- the scan window was NOT measured: steps "
                 "2-5 pass no --window and DIA-NN chooses the radius itself, per run")
MASS_ACC_SOURCE = ("fallback (probe failed: {reason}) -- the mass accuracy was NOT measured: "
                   "{levels}")

WINDOW_FALLBACK, MASS_ACC_FALLBACK = "fallback_auto", "fallback_default"
assert WINDOW_FALLBACK in SCAN_WINDOW_MODES and MASS_ACC_FALLBACK in MASS_ACC_MODES

FAILED_BECAUSE = {"io_error": "the probe could not read its own DIA-NN log (an I/O error)",
                  "timed_out": "the runs hit the probe's time limit before logging it"}

STOPPED = {"environment": "a failure no other run could fix (DIA-NN could not start or read the "
                          "runs, or would not optimise under these flags)",
           "max_failures": "the runs tried did not log it",
           "budget": "the time budget ran out",
           "no_more_runs": "every eligible run was tried",
           "no_probeable_run": "no run could be probed",
           "signal": "stopped by a signal",
           # never fallen back past (probe_window.EXIT_OOM fails the job); said if ever read
           "oom": "DIA-NN ran out of memory"}


def reason_from(evidence, exit_code, attempts):
    """One line saying why the probe measured nothing, from its evidence JSON when there is
    one -- a probe that crashed may have left none."""
    try:
        with open(evidence) as fh:
            ev = json.load(fh)
    except (OSError, ValueError, TypeError):
        ev = None
    tried = f"{attempts} attempt{'s' if attempts != 1 else ''}"
    if not isinstance(ev, dict):
        return (f"the probe crashed (exit {exit_code}) without writing its evidence ({tried})")
    stopped = ev.get("stopped_because")
    why = FAILED_BECAUSE.get(ev.get("failure")) or STOPPED.get(stopped, stopped or "unknown")
    if exit_code == EXIT_CRASH:
        why = f"the probe crashed after writing its evidence ({why})"
    probes = [p for p in ev.get("probes") or [] if isinstance(p, dict)]
    io = sorted({p["io_error"] for p in probes if p.get("io_error")})
    failed = ev.get("failed") or []
    # What the runs DID log before the probe failed -- a lone measured mass accuracy included --
    # is not pinned, but it is said: the fallback must not discard it unseen.
    logged = []
    for p in probes:
        got = ", ".join(x for x in (
            f"radius {p['radius']}" if p.get("radius") is not None else None,
            f"MS2 {ppm_text(p['ms2_ppm'])} ppm" if p.get("ms2_ppm") is not None else None,
            f"MS1 {ppm_text(p['ms1_ppm'])} ppm" if p.get("ms1_ppm") is not None else None) if x)
        if got:
            logged.append(f"{os.path.basename(str(p.get('file', '?')).rstrip('/'))} {got}")
    return (f"{why}" + (f" ({', '.join(failed)})" if failed else "")
            + (f"; I/O error on the probe log: {'; '.join(io)}" if io else "")
            + (f"; logged before it failed, not pinned: {'; '.join(logged)}" if logged else "")
            + f"; {tried}, last exit {exit_code}")


def fallback(measure, documented, reason):
    """The record of what the search runs with. `documented` is {"ms1_ppm": X} / {"ms2_ppm": Y}."""
    rec = {"fallback": True, "reason": reason, "measured": list(measure)}
    if "window" in measure:
        rec["window"] = {"mode": WINDOW_FALLBACK, "value": "auto", "fallback": True,
                         "reason": reason,
                         "source": WINDOW_SOURCE.format(reason=reason)}
    if "mass-acc" in measure:
        vals, default, doc = {}, {}, {}
        for flag, key, _ in LEVELS:
            if key in documented:
                vals[flag] = documented[key]
                doc[flag] = documented[key]
            else:
                vals[flag] = SOP_MASS_ACC[key]
                default[flag] = SOP_MASS_ACC[key]
        pin = (f"--mass-acc {ppm_text(vals['--mass-acc'])} "
               f"--mass-acc-ms1 {ppm_text(vals['--mass-acc-ms1'])}")
        assert re.fullmatch(MEASURED_FILE_RE["massacc.txt"], pin), pin
        levels = "; ".join(
            f"{lvl} {ppm_text(vals[flag])} ppm "
            + ("as documented from DIA-NN's Orbitrap resolution table" if flag in doc else
               f"-- DEFAULT, not user-confirmed: {SOP_MASS_ACC_SOURCE}")
            for flag, _, lvl in LEVELS)
        rec["mass_acc"] = {"mode": MASS_ACC_FALLBACK, "fixed": True, "measured": False,
                           "fallback": True, "reason": reason,
                           "ms2": vals["--mass-acc"], "ms1": vals["--mass-acc-ms1"],
                           "documented": doc, "default": default, "pin_as": pin,
                           "levels": levels,
                           "source": MASS_ACC_SOURCE.format(reason=reason, levels=levels)}
    return rec


def _dict(x):
    return x if isinstance(x, dict) else {}


def fallback_modes(prov):
    """(scan_window mode, mass_acc mode) of a search_provenance.json -- THE reader of whether a
    search fell back, for everything that says so: the Methods (make_methods.
    probe_fallback_record, and so the run record's rows) and the CAUTION (provenance_caution():
    AUDIT.md, the report, the run record's notes, the status views). The stable `mode`
    (diann_parallel.SCAN_WINDOW_MODES / MASS_ACC_MODES), top level first, then under `result`.
    A provenance written before the modes existed (7f907ed / 8290ab0: `probe_fallback`, no
    `mode`) gets the fallback modes derived from `probe_fallback` -- here, and nowhere else.
    Never raises on a malformed record: every piece is type-checked."""
    prov = _dict(prov)
    res = _dict(prov.get("result"))
    sw = _dict(prov.get("scan_window")) or _dict(res.get("scan_window"))
    ma = _dict(prov.get("mass_acc")) or _dict(res.get("mass_acc"))
    fb = _dict(prov.get("probe_fallback"))
    window = sw.get("mode") if isinstance(sw.get("mode"), str) else (
        WINDOW_FALLBACK if isinstance(fb.get("window"), dict) else None)
    mass_acc = ma.get("mode") if isinstance(ma.get("mode"), str) else (
        MASS_ACC_FALLBACK if isinstance(fb.get("mass_acc"), dict) else None)
    return window, mass_acc


def caution(window_mode, mass_acc_mode, reason=None, pin_as=None):
    """The CAUTION for a search that ran on a fallback, from its stable modes (fallback_modes()),
    or None. The one wording, for AUDIT.md -- and so the report's Audit & caveats -- the report's
    Results at a glance, the run record's Data Quality Notes and the status views."""
    parts = []
    if window_mode == WINDOW_FALLBACK:
        parts.append("the DIA-NN scan window was NOT measured: the measurement before the search "
                     "failed, so DIA-NN chose the window itself for each run, and DIA-NN advises "
                     "against combining runs searched that way (its own warning: \"strongly not "
                     "recommended\")")
    if mass_acc_mode == MASS_ACC_FALLBACK:
        parts.append("the mass accuracy was NOT measured: the measurement before the search "
                     "failed, so the search used the documented value and the facility SOP "
                     "default (" + (pin_as or "see the Methods") + ")")
    if not parts:
        return None
    return ("CAUTION: " + "; and ".join(parts) + ". Re-run the search if the cause was transient "
            "(e.g. an NFS `Stale file handle`). Why the measurement failed: "
            + (reason or "not recorded") + ".")


def provenance_caution(prov):
    """caution() for a search_provenance.json (or anything shaped like one), or None."""
    window, mass_acc = fallback_modes(prov)
    prov = _dict(prov)
    fb = _dict(prov.get("probe_fallback"))
    ma = _dict(prov.get("mass_acc")) or _dict(_dict(prov.get("result")).get("mass_acc"))
    pin = ma.get("pin_as") or _dict(fb.get("mass_acc")).get("pin_as")
    return caution(window, mass_acc, fb.get("reason") if isinstance(fb.get("reason"), str)
                   else None, pin if isinstance(pin, str) else None)


def caution_for(out_dir):
    """caution() for a search output directory, or None when nothing fell back there -- for the
    status views (checkpoint.py status, notify_slack.py) that have the directory. From its
    search_provenance.json, through fallback_modes() like every other reader (dda-review N1: a
    stale probe_fallback.json must not outvote the provenance); from probe_fallback.json
    (FALLBACK_RECORD) only when there is no provenance there."""
    prov_path = os.path.join(out_dir, "search_provenance.json")
    if os.path.isfile(prov_path):
        try:
            with open(prov_path) as fh:
                return provenance_caution(json.load(fh))
        except (OSError, ValueError) as e:
            return (f"CAUTION: {prov_path} could not be read ({type(e).__name__}: {e}), so "
                    "whether the search fell back instead of measuring is not known")
    try:
        with open(os.path.join(out_dir, FALLBACK_RECORD)) as fh:
            rec = json.load(fh)
    except FileNotFoundError:
        return None
    except (OSError, ValueError) as e:
        return (f"CAUTION: {FALLBACK_RECORD} in {out_dir} could not be read ({type(e).__name__}: "
                f"{e}) -- the search may have fallen back instead of measuring")
    rec = _dict(rec)
    return provenance_caution({"probe_fallback": rec, "scan_window": rec.get("window"),
                               "mass_acc": rec.get("mass_acc")})


def record_in_provenance(path, rec, files):
    """search_provenance.json: `probe_fallback`, and the scan_window / mass_acc it replaces (top
    level, and under result).
    Written through a .tmp and a rename. Returns False when there is no such file."""
    if not path or not os.path.isfile(path):
        return False
    with open(path) as fh:
        prov = json.load(fh)
    prov["probe_fallback"] = rec
    res = prov.get("result") if isinstance(prov.get("result"), dict) else None
    if "window" in rec:
        sw = dict(rec["window"], value=None, value_file=files.get("window"),
                  evidence_file=files.get("evidence"))
        prov["scan_window"] = sw
        if res is not None:
            res["scan_window"] = sw
    if "mass_acc" in rec:
        ma = dict(rec["mass_acc"], value_file=files.get("massacc"),
                  evidence_file=files.get("evidence"))
        prov["mass_acc"] = ma
        if res is not None:
            res["mass_acc"] = ma
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(prov, fh, indent=2)
    os.replace(tmp, path)
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--measure", nargs="+", choices=tuple(MEASURES), required=True)
    ap.add_argument("--ms1-ppm", type=float, help="the documented MS1 level the probe was given")
    ap.add_argument("--ms2-ppm", type=float, help="the documented MS2 level the probe was given")
    ap.add_argument("--exit-code", type=int, required=True, help="the probe's last exit status")
    ap.add_argument("--attempts", type=int, default=1)
    ap.add_argument("--evidence", help="the probe's evidence JSON (may be missing or partial)")
    ap.add_argument("--window-file", help="window.txt: gets `auto`")
    ap.add_argument("--massacc-file", help="massacc.txt: gets the fallback pair")
    ap.add_argument("--write-cfg", help="the resolved cfg being built: gets the mass-acc flags")
    ap.add_argument("--provenance", help="search_provenance.json to record the fallback in")
    ap.add_argument("--out", help="write the record here as JSON")
    a = ap.parse_args(argv)
    measure = [m for m in MEASURES if m in a.measure]
    documented = {k: v for k, v in (("ms1_ppm", a.ms1_ppm), ("ms2_ppm", a.ms2_ppm))
                  if v is not None}
    reason = reason_from(a.evidence, a.exit_code, a.attempts)
    rec = fallback(measure, documented, reason)
    files = {"window": a.window_file, "massacc": a.massacc_file, "evidence": a.evidence}
    try:
        if "window" in measure and a.window_file:
            with open(a.window_file, "w") as fh:
                fh.write("auto\n")
        if "mass-acc" in measure:
            if a.massacc_file:
                with open(a.massacc_file, "w") as fh:
                    fh.write(rec["mass_acc"]["pin_as"] + "\n")
            if a.write_cfg:
                with open(a.write_cfg, "a") as fh:
                    fh.write(f"\n--mass-acc {ppm_text(rec['mass_acc']['ms2'])}\n"
                             f"--mass-acc-ms1 {ppm_text(rec['mass_acc']['ms1'])}\n")
        if a.out:
            # indent=2, as steps 2-5's guard greps it: `"mode": "fallback_auto"`
            # (diann_parallel.needs_measured) accepts window.txt `auto` only beside this record
            with open(a.out, "w") as fh:
                json.dump(rec, fh, indent=2)
    except OSError as e:
        sys.exit(f"[probe_fallback] FAILED: the fallback could not be written ({e}); the search "
                 "cannot run without it")
    recorded = record_in_provenance(a.provenance, rec, files)
    sys.stderr.write(
        f"WARNING: the probe measured nothing ({reason}). FALLBACK: "
        + "; ".join(x for x in (
            "the scan window is left to DIA-NN, per run (window.txt: auto)"
            if "window" in measure else None,
            f"mass accuracy {rec['mass_acc']['pin_as']} ({rec['mass_acc']['levels']})"
            if "mass-acc" in measure else None) if x)
        + (". Recorded in search_provenance.json (probe_fallback)." if recorded else
           ". No search_provenance.json to record it in: the record is "
           + (a.out or "only in this log") + ".") + "\n")
    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
