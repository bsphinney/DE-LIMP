#!/usr/bin/env python3
"""
sage_lfq_check.py -- does Sage's label-free quantification window fit the data's MS1 mass error?

Sage integrates each identified peptide's MS1 isotope envelope only within
+/- quant.lfq_settings.ppm_tolerance of the THEORETICAL mass (Sage 0.14.7, crates/sage/src/lfq.rs
build_feature_map: Tolerance::Ppm(-ppm_tolerance, ppm_tolerance) around `calcmass`; default 5.0,
LfqSettings::default). The search itself uses precursor_tol (+/-10 ppm here), so a run whose MS1
sits a few ppm off identifies normally and quantifies badly -- and nothing says so.

gabrig 2026-09-29 (HeLa, Fusion Lumos, Sage 0.14.7, skill 2.6.0): runs with a systematic MS1
offset of +7.1/+7.7 ppm searched fine (29,846 PSMs) but LFQ logged "discovered 0 target MS1 peaks
at 5% FDR" -- target and decoy MS1 intensities were the same (1.85e5 vs 1.89e5). Runs at +1.8 to
+5.8 ppm quantified, but with run-dependent loss: protein CV 50% against MSFragger's 14% on the
same files. No warning anywhere.

So after a Sage search this reads, per run, the median precursor mass error of its confident
target PSMs (results.sage.parquet) and the number of target MS1 peaks Sage kept at 5% FDR (its
log, else lfq.parquet), and WARNS when

  * a run's median |precursor_ppm| + MARGIN_PPM exceeds the LFQ ppm_tolerance, or
  * Sage kept 0 target MS1 peaks, or fewer than LOW_PEAK_FRACTION of the target peptides it
    identified at 1% FDR (the peptides LFQ tries to quantify).

It never changes or re-runs anything: the warning names the fix (a wider
quant.lfq_settings.ppm_tolerance, or recalibrated mzML). The record goes to
<out>/sage_lfq_check.json and into <out>/search_provenance.json (`sage_lfq_check`); watch_run.sh
--out, checkpoint.py status and audit_results.py --search-out (-> AUDIT.md -> the report's
"Audit & caveats") read it from there.

Columns, from Sage 0.14.7 crates/sage-cloudpath/src/parquet.rs:
  results.sage.parquet: filename, is_decoy, rank, peptide, peptide_q, precursor_ppm (ABSOLUTE:
    scoring.rs takes .abs()), expmass, calcmass, isotope_error (Da). The signed offset is
    (expmass - isotope_error - calcmass) / calcmass * 1e6 and is reported for direction only.
  lfq.parquet: peptide, charge, is_decoy, q_value, filename, intensity.
  The log line is sage-cli/src/main.rs: "discovered {} target MS1 peaks at 5% FDR".

Usage:
  python3 sage_lfq_check.py --out <sage output dir> [--params sage_config.json] [--log sage.log]
Exit 0 whatever it finds (a warning is not a failed search); 2 on a usage error.
"""
import argparse
import datetime
import glob
import json
import math
import os
import re
import sys

# Sage 0.14.7 LfqSettings::default().ppm_tolerance (crates/sage/src/lfq.rs). Used only when neither
# Sage's results.json nor the config states one, and then tagged as the default in the record.
SAGE_DEFAULT_LFQ_PPM = 5.0
# How much room the median needs inside the window: the envelope is spread around the median,
# so a median right at the edge leaves about half of it outside. 2 ppm is the margin gabrig's
# diagnosis used (runs at +1.8..+5.8 ppm already quantified worse than MSFragger against 5 ppm).
MARGIN_PPM = 2.0
# Heuristic, not measured: a healthy search quantifies most of what it identified, so under 10%
# is a collapse. Zero is always a warning.
LOW_PEAK_FRACTION = 0.10
# Fewer confident PSMs than this in a run and its median is not used to judge it.
MIN_PSMS = 50
# Sage warns above this itself (sage-cli/src/input.rs: "lfq_settings.ppm_tolerance is higher
# than expected").
SAGE_WARNS_ABOVE_PPM = 20.0
PSM_Q = 0.01            # the peptide_q Sage's LFQ feature map takes (lfq.rs build_feature_map)
# Sage's own line for a valid MS1 peak: picked_precursor counts a target as discovered at
# q_value <= 0.05 (crates/sage/src/fdr.rs), logged as "discovered N target MS1 peaks at 5% FDR"
# (sage-cli main.rs). Sage documents no other LFQ threshold (sage-docs results/lfq.mdx: q_value is
# the "MS1-integration specific q-value"), and its lfq.tsv writer filters on none. ONE definition:
# this check counts peaks with it, and run_search.adapt_sage keeps LFQ rows with it.
#
# The q_value Sage STORES is not the one it counted with: picked_precursor keys its q map by
# PrecursorId alone, so a target and its decoy end with ONE shared q -- the lower-scoring member's
# (the map is filled in score order, and a later insert wins). Stored q >= the target's own q, so
# filtering the stored q keeps a conservative SUBSET of Sage's logged count: gabrig's HeL50 UnvPe
# (2026-09-29) had 1,489 of 1,489 target/decoy pairs with identical q, and 16,649 targets kept
# against Sage's logged 17,288. No FDR harm; Sage's count is recorded beside what was kept.
LFQ_Q_MAX = 0.05


def q_passes(q):
    """A pyarrow boolean array: q_value <= LFQ_Q_MAX, compared in float32 as Sage compares it
    (fdr.rs: `q_min <= 0.05` on an f32; lfq.parquet stores q_value as float). In float64 a q of
    exactly 0.05 -- 0.0500000007 as an f32 -- would fail the line Sage says it passes. Nulls fail."""
    import pyarrow as pa
    import pyarrow.compute as pc
    return pc.fill_null(pc.less_equal(pc.cast(q, pa.float32()),
                                      pa.scalar(LFQ_Q_MAX, pa.float32())), False)
RECORD = "sage_lfq_check.json"
PROVENANCE = "search_provenance.json"
TAG = "[sage_lfq_check]"
LOG_RE = re.compile(r"discovered (\d+) target MS1 peaks at 5% FDR")


def _find(out, name):
    """<out>/<name>, else the first one under <out> (Sage writes into -o, but a hand-run
    search may have used the config's output_directory below it)."""
    p = os.path.join(out, name)
    if os.path.isfile(p):
        return p
    for dp, _dns, fns in os.walk(out):
        if name in fns:
            return os.path.join(dp, name)
    return None


def lfq_settings(out, params=None):
    """(lfq enabled or None, ppm_tolerance, where the tolerance came from).

    Sage's own results.json first -- it is every parameter as Sage actually ran it, defaults
    included -- then the config handed to Sage, then Sage's built-in default, said so."""
    for path, label in ((_find(out, "results.json"), "Sage's results.json (as run)"),
                        (params, "the Sage config")):
        if not path:
            continue
        try:
            with open(path) as fh:
                cfg = json.load(fh)
        except (OSError, ValueError):
            continue
        quant = cfg.get("quant") if isinstance(cfg, dict) else None
        if not isinstance(quant, dict):
            continue
        lfq = quant.get("lfq")
        tol = (quant.get("lfq_settings") or {}).get("ppm_tolerance")
        if isinstance(tol, (int, float)):
            return lfq, abs(float(tol)), f"{label}: quant.lfq_settings.ppm_tolerance"
        return lfq, SAGE_DEFAULT_LFQ_PPM, (
            f"DEFAULT -- not set in {label}; Sage 0.14.7's built-in LfqSettings default")
    return None, SAGE_DEFAULT_LFQ_PPM, (
        "DEFAULT -- no Sage results.json or config was readable; Sage 0.14.7's built-in "
        "LfqSettings default")


def _median(xs):
    xs = sorted(xs)
    n = len(xs)
    if not n:
        return None
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


def per_run_mass_error(results):
    """({run: {"psms", "median_abs_ppm", "median_offset_ppm"}}, target peptides at 1%).

    Confident target PSMs only: not decoy, peptide_q <= PSM_Q (what LFQ quantifies), rank 1 for
    the mass error (a chimeric second match has its own precursor)."""
    import pyarrow.parquet as pq
    import pyarrow.compute as pc
    cols = ["filename", "is_decoy", "rank", "peptide", "peptide_q", "precursor_ppm",
            "expmass", "calcmass", "isotope_error"]
    t = pq.read_table(results, columns=cols)
    t = t.filter(pc.and_(pc.invert(t["is_decoy"]), pc.less_equal(t["peptide_q"], PSM_Q)))
    peptides = len(set(t["peptide"].to_pylist()))
    t = t.filter(pc.equal(t["rank"], 1))
    runs = {}
    for f, ppm, em, cm, iso in zip(t["filename"].to_pylist(), t["precursor_ppm"].to_pylist(),
                                   t["expmass"].to_pylist(), t["calcmass"].to_pylist(),
                                   t["isotope_error"].to_pylist()):
        r = runs.setdefault(os.path.basename(str(f)), ([], []))
        if ppm is not None:
            r[0].append(abs(float(ppm)))
        if em is not None and cm:
            r[1].append((float(em) - float(iso or 0.0) - float(cm)) / float(cm) * 1e6)
    return {f: {"psms": len(a),
                "median_abs_ppm": _round(_median(a)),
                "median_offset_ppm": _round(_median(s))}
            for f, (a, s) in sorted(runs.items())}, peptides


def _round(x):
    return None if x is None else round(x, 2)


def ms1_peaks_from_log(log):
    """The LAST 'discovered N target MS1 peaks at 5% FDR' in Sage's log, or None."""
    try:
        with open(log, errors="replace") as fh:
            hits = LOG_RE.findall(fh.read())
    except OSError:
        return None
    return int(hits[-1]) if hits else None


def ms1_peaks_from_lfq(lfq):
    """Target precursors with stored q_value <= LFQ_Q_MAX in lfq.parquet (one row per precursor x
    run). A LOWER BOUND on Sage's own count: the stored q is shared with the decoy (above)."""
    import pyarrow.parquet as pq
    import pyarrow.compute as pc
    t = pq.read_table(lfq, columns=["peptide", "charge", "is_decoy", "q_value"])
    t = t.filter(pc.and_(pc.invert(t["is_decoy"]), q_passes(t["q_value"])))
    return len(set(zip(t["peptide"].to_pylist(), t["charge"].to_pylist())))


def find_log(out, log=None):
    """The Sage log: --log, else <out>/sage.log (run_search.py tees Sage into it), else the
    newest SLURM log of a Sage job run_search.py wrote (<out>/sage_search_<id>.log)."""
    if log:
        return log if os.path.isfile(log) else None
    p = os.path.join(out, "sage.log")
    if os.path.isfile(p):
        return p
    jobs = sorted(glob.glob(os.path.join(out, "sage_search_*.log")), key=os.path.getmtime)
    return jobs[-1] if jobs else None


def check(out, params=None, log=None):
    """The record (see the module docstring). Never raises for a missing or unreadable input:
    that is status "unchecked", with the reason, so it is said rather than skipped."""
    rec = {"check": "sage_lfq_window", "out": os.path.abspath(out),
           "checked_at": datetime.datetime.now().isoformat(timespec="seconds"),
           "margin_ppm": MARGIN_PPM, "low_peak_fraction": LOW_PEAK_FRACTION,
           "min_psms_per_run": MIN_PSMS}
    lfq_on, tol, tol_src = lfq_settings(out, params)
    rec.update({"lfq_enabled": lfq_on, "ppm_tolerance": tol, "ppm_tolerance_source": tol_src})
    if lfq_on is False:
        rec.update(status="not_applicable", message="Sage ran without LFQ (quant.lfq false).")
        return rec
    results = _find(out, "results.sage.parquet")
    rec["sage_results"] = results
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        return _unchecked(rec, "pyarrow is not installed in this Python")
    if not results:
        return _unchecked(rec, f"no results.sage.parquet under {out} (did Sage finish, and "
                               "with --parquet?)")
    try:
        runs, peptides = per_run_mass_error(results)
    except Exception as e:          # an unreadable parquet is reported, never a crash
        return _unchecked(rec, f"could not read {results} ({type(e).__name__}: {e})")
    rec["runs"] = runs
    rec["target_peptides_1pct"] = peptides

    logp = find_log(out, log)
    peaks, src = (ms1_peaks_from_log(logp), f"Sage's log {logp}") if logp else (None, None)
    if peaks is None:
        lfq = _find(out, "lfq.parquet")
        if lfq:
            try:
                peaks, src = ms1_peaks_from_lfq(lfq), (
                    f"{lfq}: target precursors with stored q_value <= {LFQ_Q_MAX} -- a lower "
                    f"bound on Sage's own count (the stored q is shared with the decoy)")
            except Exception as e:
                src = f"{lfq} unreadable ({type(e).__name__}: {e})"
        else:
            src = "neither Sage's log nor lfq.parquet was found"
    rec["target_ms1_peaks_5pct_fdr"] = peaks
    rec["target_ms1_peaks_source"] = src

    flagged = sorted(f for f, r in runs.items()
                     if r["psms"] >= MIN_PSMS and r["median_abs_ppm"] is not None
                     and r["median_abs_ppm"] + MARGIN_PPM > tol)
    too_few = sorted(f for f, r in runs.items() if r["psms"] < MIN_PSMS)
    for f, r in runs.items():
        r["outside_lfq_window"] = f in flagged
    low_peaks = peaks is not None and (peaks == 0 or (
        peptides and peaks < LOW_PEAK_FRACTION * peptides))
    rec["runs_outside_lfq_window"] = flagged
    rec["runs_too_few_psms"] = too_few
    rec["low_ms1_peaks"] = bool(low_peaks)
    unjudged = (f" ({len(too_few)} run(s) with fewer than {MIN_PSMS} confident target PSMs "
                f"not judged: {', '.join(too_few[:6])})" if too_few else "")
    if not (flagged or low_peaks):
        if len(too_few) == len(runs):
            return _unchecked(rec, f"no run has {MIN_PSMS} confident target PSMs (peptide_q <= "
                                   f"{PSM_Q}) to measure its precursor mass error from")
        worst = max((r["median_abs_ppm"] for f, r in runs.items()
                     if f not in too_few and r["median_abs_ppm"] is not None), default=None)
        rec.update(status="ok", message=(
            f"Sage LFQ window fits the data: every run's median precursor mass error "
            f"({_ppm(worst)} at most) is at least {MARGIN_PPM:g} ppm inside the +/-{tol:g} ppm "
            f"LFQ window{unjudged}; {_peaks_text(peaks, peptides)}."))
        return rec
    rec["status"] = "warn"
    rec["suggested_ppm_tolerance"] = suggest(runs, flagged, tol)
    rec["message"] = warn_message(rec, runs, flagged, peaks, peptides, tol, tol_src) + (
        f" Not judged{unjudged.replace(' not judged', '')}." if too_few else "")
    return rec


def suggest(runs, flagged, tol):
    """The smallest whole-ppm window this check accepts for the worst run (never below the
    current one + 1). A wider window is Sage's own knob for this; above 20 ppm Sage warns."""
    worst = max((runs[f]["median_abs_ppm"] for f in flagged), default=None)
    if worst is None:
        return None
    return max(int(math.ceil(worst + MARGIN_PPM)), int(math.floor(tol)) + 1)


def _ppm(x):
    return "unknown" if x is None else f"{x:.1f} ppm"


def _peaks_text(peaks, peptides):
    if peaks is None:
        return "the number of MS1 peaks Sage kept is unknown"
    return (f"Sage kept {peaks:,} target MS1 peaks at 5% FDR"
            + (f" for {peptides:,} target peptides identified at 1% FDR" if peptides else ""))


def warn_message(rec, runs, flagged, peaks, peptides, tol, tol_src):
    parts = []
    if flagged:
        show = "; ".join(f"{f} {runs[f]['median_abs_ppm']:.1f} ppm"
                         + (f" ({runs[f]['median_offset_ppm']:+.1f} signed)"
                            if runs[f]["median_offset_ppm"] is not None else "")
                         for f in flagged[:6])
        more = f" and {len(flagged) - 6} more" if len(flagged) > 6 else ""
        parts.append(f"{len(flagged)} of {len(runs)} run(s) have a median precursor mass error "
                     f"within {MARGIN_PPM:g} ppm of, or beyond, the +/-{tol:g} ppm window Sage "
                     f"integrates MS1 in ({show}{more}), so part or all of their isotope "
                     f"envelopes fall outside it")
    parts.append(_peaks_text(peaks, peptides)
                 + (" -- LFQ collapsed" if rec.get("low_ms1_peaks") else ""))
    s = rec.get("suggested_ppm_tolerance")
    if not flagged:
        # the mass error fits, so a wider window is not the known fix -- say so, not guess
        fix = ("the runs' median precursor mass errors fit the window, so this check cannot "
               "name the cause; inspect the MS1 data (centroided MS1 present in the mzML?) and "
               "Sage's RT alignment before trusting the quantities")
    else:
        fix = ("re-run Sage with a wider window: add \"lfq_settings\": {\"ppm_tolerance\": "
               f"{s}}} to the \"quant\" block of the Sage config" if s else
               "re-run Sage with a wider quant.lfq_settings.ppm_tolerance")
        if s and s > SAGE_WARNS_ABOVE_PPM:
            fix += (f" (Sage warns above {SAGE_WARNS_ABOVE_PPM:g} ppm; recalibrating the mzML "
                    "is the better fix at that size)")
        else:
            fix += ", or recalibrate the mzML"
    return (f"Sage label-free quantification is unreliable: {'; '.join(parts)}. The window is "
            f"{tol:g} ppm ({tol_src}). The identifications are not affected, only the MS1 "
            f"quantities. "
            + (f"Fix: {fix}, then run_search.py --adapt-only again." if flagged else
               f"Next: {fix}.")
            + " Nothing was re-run automatically.")


def _unchecked(rec, why):
    rec.update(status="unchecked", reason=why,
               message=f"The Sage LFQ mass-window check could not run: {why}. Look for "
                       f"'target MS1 peaks at 5% FDR' in Sage's log by hand; 0 there means the "
                       f"search has no usable quantities.")
    return rec


def load(out):
    """The record written for this Sage output dir, or None (no Sage search, or not checked)."""
    try:
        with open(os.path.join(out, RECORD)) as fh:
            rec = json.load(fh)
    except (OSError, ValueError):
        return None
    return rec if isinstance(rec, dict) else None


def audit_finding(rec):
    """(status, message, detail) for audit_results.py, or None. ONE mapping: warn -> WARN,
    ok -> PASS, unchecked -> INFO (a check that could not run never changes the verdict)."""
    if not rec or rec.get("status") not in ("warn", "ok", "unchecked"):
        return None
    status = {"warn": "WARN", "ok": "PASS", "unchecked": "INFO"}[rec["status"]]
    detail = {k: rec.get(k) for k in ("ppm_tolerance", "ppm_tolerance_source",
                                      "target_ms1_peaks_5pct_fdr", "target_peptides_1pct",
                                      "runs_outside_lfq_window", "suggested_ppm_tolerance")
              if rec.get(k) is not None}
    return status, rec.get("message") or "", detail


def say(rec, stream=None):
    """The one line (or two) a person sees: WARNING for warn, NOTE for unchecked, INFO for ok."""
    stream = stream or sys.stderr
    level = {"warn": "WARNING", "unchecked": "NOTE", "ok": "INFO"}.get(rec.get("status"))
    if level and rec.get("message"):
        print(f"{TAG} {level}: {rec['message']}", file=stream, flush=True)


def write_record(out, name, key, rec):
    """Write <out>/<name> and put the same record into search_provenance.json under `key` when
    that exists (merged, never replaced). The ONE way a Sage step records itself after the
    search: the job and --adapt-only run long after run_search.py wrote the provenance.
    Returns the list of files written."""
    written = []
    path = os.path.join(out, name)
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh, indent=2)
    os.replace(path + ".tmp", path)
    written.append(path)
    prov = os.path.join(out, PROVENANCE)
    try:
        with open(prov) as fh:
            p = json.load(fh)
    except (OSError, ValueError):
        return written
    if isinstance(p, dict):
        p[key] = rec
        with open(prov + ".tmp", "w") as fh:
            json.dump(p, fh, indent=2)
        os.replace(prov + ".tmp", prov)
        written.append(prov)
    return written


def record(out, rec):
    """<out>/sage_lfq_check.json + search_provenance.json `sage_lfq_check`."""
    return write_record(out, RECORD, "sage_lfq_check", rec)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="the Sage output directory (-o)")
    ap.add_argument("--params", help="the Sage config the search ran with")
    ap.add_argument("--log", help="Sage's log (default: <out>/sage.log, else the newest "
                                  "<out>/sage_search_*.log)")
    a = ap.parse_args(argv)
    rec = check(a.out, a.params, a.log)
    try:
        rec["written_to"] = record(a.out, rec)
    except OSError as e:
        print(f"{TAG} NOTE: could not write {RECORD} in {a.out}: {e}", file=sys.stderr)
    say(rec)
    print(json.dumps(rec, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
