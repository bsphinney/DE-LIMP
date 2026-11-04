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

It never changes or re-runs the search. The record goes to <out>/sage_lfq_check.json and into
<out>/search_provenance.json (`sage_lfq_check`); watch_run.sh --out, checkpoint.py status and
audit_results.py --search-out (-> AUDIT.md -> the report's "Audit & caveats") read it from there.

When it warns BECAUSE OF THE MASS ERROR (runs outside the window), the quantities are known to
be wrong and the fix is known, so it is a GATE, not only a warning (2.10; it only warned in 2.9
and the quantities went on into the DE):
  * the corrected window is computed from the measured offsets (suggest(): the worst run's
    median + MARGIN_PPM) and written into a copy of the config the search ran with,
    <out>/sage_config.lfq_ppm<N>.json (`corrected_params`);
  * refusal() records the exact run_search.py command that repeats the search with it into a
    new folder, <out>_lfq_ppm<N> (`rerun`), from what search_provenance.json recorded;
  * run_search.py --adapt-only (and an inline search) will not build report.parquet, and
    run_de.R will not read one beside this record, until either that re-run is used or the
    user accepts these quantities as they are: run_search.py --adapt-only --accept-lfq-window
    "<why>" (accept(), recorded as `accepted`; who is the login that ran it, named to clients by
    role -- never typed into the reason, which reaches AUDIT.md and commands.log verbatim).
Nothing is re-run automatically: a wider window also admits more co-eluting interference, so
the re-run's target/decoy separation has to be compared, and a changed search setting is a new
search to confirm. A low peak count with a fitting window stays a warning: its cause is unknown,
so there is no corrected run to point to.

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


def methods_sentence(q_value_max=LFQ_Q_MAX):
    """A quantifying Sage search's identification FDR in a Methods section's words: the rule
    run_search.adapt_sage keeps report rows with (PSM_Q, the line Sage's feature map takes -- lfq.rs
    build_feature_map, `peptide_q <= 0.01 && label == 1` at v0.14.7 -- and the MS1-peak q it kept,
    `q_value_max`, LFQ_Q_MAX unless a record says otherwise). adapt_sage writes it into
    sage_adapt.json `methods_sentence`; make_methods.py quotes it, or builds it with this from an
    older record's own q_value_max. The identifications-only route has its own
    (run_search.SAGE_IDS_METHODS)."""
    return (f"Peptides were quantified from Sage's label-free MS1 quantification: Sage traced "
            f"the MS1 signal of the peptides it identified at {PSM_Q * 100:g}% FDR (peptide-level "
            f"q-value <= {PSM_Q:g}, computed over all runs together), and an MS1 peak was kept "
            f"when it was a target peak with Sage's MS1-integration q-value <= {q_value_max:g}; "
            f"no protein-level FDR filter was applied.")


def q_passes(q):
    """A pyarrow boolean array: q_value <= LFQ_Q_MAX, compared in float32 as Sage compares it
    (fdr.rs: `q_min <= 0.05` on an f32; lfq.parquet stores q_value as float). In float64 a q of
    exactly 0.05 -- 0.0500000007 as an f32 -- would fail the line Sage says it passes. Nulls fail."""
    import pyarrow as pa
    import pyarrow.compute as pc
    return pc.fill_null(pc.less_equal(pc.cast(q, pa.float32()),
                                      pa.scalar(LFQ_Q_MAX, pa.float32())), False)
RECORD = "sage_lfq_check.json"
# What a Sage search for identifications only (quant.lfq false) writes instead of report.parquet
# (run_search.py adapt_sage_ids; run_search.SAGE_IDS is this name).
IDENTIFICATIONS = "identifications.tsv"
NOT_APPLICABLE_WHY = "identifications only (quant.lfq is false)"
PROVENANCE = "search_provenance.json"
# The config with the corrected window, beside the search, and the folder its re-run goes to.
CORRECTED_CONFIG = "sage_config.lfq_ppm{ppm}.json"
RERUN_SUFFIX = "_lfq_ppm{ppm}"
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
        return not_applicable(rec)
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
    # Nothing identified: there is no run whose mass error could be measured, and no quantity to
    # judge -- said as such, never as a window that "fits" or identifications "not affected"
    # (H2-3: an un-combined timsTOF mzML gave 0 peptides, and the record said both).
    if not peptides:
        return _unchecked(rec, f"Sage identified no target peptide at {PSM_Q * 100:g}% FDR "
                               f"(peptide_q <= {PSM_Q}) in any run, so no run's precursor mass "
                               f"error could be measured and the LFQ window was not checked; "
                               f"with no identifications there is nothing to quantify either")
    if not runs:
        return _unchecked(rec, f"no run has a rank-1 target PSM at peptide_q <= {PSM_Q}, so no "
                               f"run's precursor mass error could be measured")

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
    if flagged and rec["suggested_ppm_tolerance"]:
        rec["corrected_params"], why = corrected_config(out, params,
                                                        rec["suggested_ppm_tolerance"])
        if why:
            rec["corrected_params_why"] = why
    rec["message"] = warn_message(rec, runs, flagged, peaks, peptides, tol, tol_src) + (
        f" Not judged{unjudged.replace(' not judged', '')}." if too_few else "")
    return rec


def _load_json(path):
    try:
        with open(path) as fh:
            v = json.load(fh)
    except (OSError, ValueError, TypeError):
        return None
    return v if isinstance(v, dict) else None


def corrected_config(out, params, ppm):
    """Write <out>/sage_config.lfq_ppm<N>.json: the config the search ran with (`params`, else
    the params_file search_provenance.json names) with quant.lfq_settings.ppm_tolerance = N, and
    nothing else changed -> (path, None), or (None, why) when no config could be read."""
    src = params or (_load_json(os.path.join(out, PROVENANCE)) or {}).get("params_file")
    cfg = _load_json(src) if src else None
    if cfg is None:
        return None, (f"the Sage config the search ran with ({src}) could not be read"
                      if src else "the Sage config the search ran with is not known (no "
                                  "--params, and search_provenance.json names none)")
    quant = cfg.setdefault("quant", {})
    if not isinstance(quant, dict):
        return None, f"{src} has a quant entry that is not a block"
    settings = quant.setdefault("lfq_settings", {})
    if not isinstance(settings, dict):
        return None, f"{src} has a quant.lfq_settings entry that is not a block"
    settings["ppm_tolerance"] = float(ppm)
    path = os.path.join(out, CORRECTED_CONFIG.format(ppm=ppm))
    with open(path + ".tmp", "w") as fh:
        json.dump(cfg, fh, indent=2)
        fh.write("\n")
    os.replace(path + ".tmp", path)
    # ...and the params side record (estimate_params.py's <params>.rationale.json: cysteine
    # alkylation, SILAC labels, defaults) beside it, so the re-run's Methods keep what it says
    # (follow-up 3, H7 bug 3: "Cysteine was searched unmodified" was lost). Only one that exists
    # is carried, with what changed; none is invented.
    side = _load_json(f"{src}.rationale.json")
    if isinstance(side, dict):
        side = dict(side, derived_from={"params": src, "changed": {
            "quant.lfq_settings.ppm_tolerance": float(ppm)}})
        with open(path + ".rationale.json.tmp", "w") as fh:
            json.dump(side, fh, indent=2)
            fh.write("\n")
        os.replace(path + ".rationale.json.tmp", path + ".rationale.json")
    return path, None


def rerun_plan(out, corrected, ppm):
    """The run_search.py command that repeats this search with the corrected window, into a new
    folder beside it, from what search_provenance.json recorded. -> {"out", "command"}, with
    command None and `why` when the record lacks what the command needs (a search from before
    2.10 records no tools.json)."""
    import shlex
    out = os.path.abspath(out)
    new_out = out.rstrip(os.sep) + RERUN_SUFFIX.format(ppm=ppm)
    prov = _load_json(os.path.join(out, PROVENANCE)) or {}
    need = {"tools": prov.get("tools"), "bundle": prov.get("bundle"),
            "fasta": prov.get("fasta"), "files": prov.get("files")}
    missing = [k for k, v in need.items() if not v]
    if not corrected:
        missing.insert(0, "the corrected config")
    if missing:
        return {"out": new_out, "command": None,
                "why": f"not recorded: {', '.join(missing)} (search_provenance.json of a search "
                       f"from before 2.10, or none). Run the search's own run_search.py command "
                       f"again with --params {corrected or '<the config, with quant.lfq_settings.ppm_tolerance ' + str(ppm) + '>'} "
                       f"and --out {new_out}"}
    argv = ["python3", os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_search.py"),
            "--tools", need["tools"], "--bundle", need["bundle"], "--engine", "sage",
            "--params", corrected, "--fasta", need["fasta"], "--out", new_out,
            "--threads", str(prov.get("threads") or 8)]
    if prov.get("keratin_sample") is True:
        argv.append("--keratin-sample")
    for k in ("partition", "account", "qos"):
        v = (prov.get("queue") or {}).get(k)
        if v:
            argv += [f"--{k}", str(v)]
    # the nodes the search was kept off (run_search.py --exclude): the re-run stays off them too
    if (prov.get("queue") or {}).get("exclude"):
        argv += ["--exclude", ",".join(prov["queue"]["exclude"])]
    # the memory it asked for (run_search.py --mem <GB>, recorded as queue.mem "<GB>G"; nothing is
    # recorded for the default, and nothing is passed then)
    mem = str((prov.get("queue") or {}).get("mem") or "")
    if mem:
        argv += ["--mem", mem[:-1] if mem.upper().endswith("G") else mem]
    # ...and the wall time (run_search.py --time HOURS, queue.time); a digest-sized default is not
    # in `queue`, so the re-run recomputes it from its own config
    hours = (prov.get("queue") or {}).get("time")
    if hours:
        argv += ["--time", str(hours)]
    # a .d converted in a container other than the pin's: the re-run converts with it again
    # (run_search.py records the flag and its value, pwiz_rerun_args)
    argv += [str(x) for x in (prov.get("pwiz_sif_arg") or [])]
    # ...and its job-end choices (Slack, FRAN, QC), as recorded, through the one mapping
    # (notify_slack.job_end_argv; follow-up 3, H7 bug 4)
    je_note = None
    try:
        import notify_slack
        argv += notify_slack.job_end_argv(prov.get("job_end_hook"))
    except ImportError:
        # a partial copy of the skill: the command is still printed, and says what it lacks
        je_note = ("notify_slack.py is not beside sage_lfq_check.py, so the search's job-end "
                   "choices (search_provenance.json job_end_hook) are not in this command: add "
                   "--no-notify / --no-fran / --qc or --not-qc as the search had them")
    job = prov.get("submitted_sbatch")
    if job:
        stem, ext = os.path.splitext(job)
        argv += ["--sbatch", stem + RERUN_SUFFIX.format(ppm=ppm) + (ext or ".sh")]
    argv += ["--files", *need["files"]]
    plan = {"out": new_out, "command": " ".join(shlex.quote(str(x)) for x in argv)}
    if je_note:
        plan["note"] = je_note
    return plan


def gated(rec):
    """True when report.parquet must not be built from this search: the window did not fit the
    runs' MS1 mass error (runs outside it, so the quantities are known to be wrong and the fix
    is known) and nobody has accepted the quantities as they are."""
    return bool(rec and rec.get("status") == "warn" and rec.get("runs_outside_lfq_window")
                and not rec.get("accepted"))


def refusal(out, rec):
    """The gate's message for a gated() record, with the exact re-run command (rerun_plan),
    recorded into the record and search_provenance.json (`rerun`, `refusal`). None when the gate
    is open."""
    if not gated(rec):
        return None
    s = rec.get("suggested_ppm_tolerance")
    plan = rec.get("rerun") or rerun_plan(out, rec.get("corrected_params"), s)
    rec["rerun"] = plan
    worst = max((rec["runs"][f]["median_abs_ppm"] for f in rec["runs_outside_lfq_window"]
                 if f in (rec.get("runs") or {})), default=None)
    how = (f"\n  {plan['command']}\nthen run_search.py --adapt-only on {plan['out']}" if
           plan.get("command") else f" {plan.get('why')}")
    rec["refusal"] = (
        f"REFUSED: report.parquet is not built from this Sage search -- {len(rec['runs_outside_lfq_window'])} "
        f"run(s) sit outside its +/-{rec.get('ppm_tolerance'):g} ppm LFQ window, so their MS1 "
        f"quantities are known to be wrong. The corrected window, from the measured offsets: "
        f"+/-{s} ppm (the worst run's median {_ppm(worst)} + {MARGIN_PPM:g} ppm, rounded up)"
        + (f", in {rec['corrected_params']}" if rec.get("corrected_params") else
           f" ({rec.get('corrected_params_why')})")
        + f". Re-run the search with it:{how}. Before using the re-run, compare its target MS1 "
        "peaks and target/decoy separation with this one: a wider window also admits more "
        "interference. Or, if the user accepts these quantities as they are, run this "
        "--adapt-only again with --accept-lfq-window \"<why, in the user's words>\" (recorded; "
        "who accepted is the login that runs it, named to clients by role -- put no name in "
        "the reason).")
    try:
        record(out, rec)
    except OSError as e:
        print(f"{TAG} NOTE: could not record the re-run plan in {out}: {e}", file=sys.stderr)
    return rec["refusal"]


def accept(out, rec, reason, by=None):
    """The user accepts a gated record's quantities as they are: `accepted` = who, when and why,
    recorded; the record stays a warning (AUDIT.md says it was accepted, and by whom)."""
    if not (reason and reason.strip()):
        raise ValueError("--accept-lfq-window needs the reason, in the user's words")
    if by is None:
        import getpass
        try:
            by = getpass.getuser()
        except (KeyError, OSError):
            by = None
    rec["accepted"] = {"reason": reason.strip(), "by": by,
                       "by_source": "the login that ran run_search.py --accept-lfq-window",
                       "at": datetime.datetime.now().isoformat(timespec="seconds")}
    record(out, rec)
    return rec


def carry_acceptance(prev, rec):
    """A re-check of the same search keeps an earlier acceptance while what was accepted -- the
    window and the runs outside it -- is unchanged."""
    if (prev and prev.get("accepted") and rec.get("status") == "warn"
            and prev.get("ppm_tolerance") == rec.get("ppm_tolerance")
            and prev.get("runs_outside_lfq_window") == rec.get("runs_outside_lfq_window")):
        rec["accepted"] = prev["accepted"]
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
        # the mass error fits, so a wider window is not the known fix -- say so, not guess; and
        # when no run could be measured, say THAT, not that it fits (H2-3)
        measured = [f for f in runs if f not in (rec.get("runs_too_few_psms") or [])]
        fix = (("the runs' median precursor mass errors fit the window, so this check cannot "
                "name the cause" if measured else
                f"no run's precursor mass error could be measured (each has fewer than "
                f"{MIN_PSMS} confident target PSMs), so this check cannot tell whether the "
                f"window fits, nor name the cause")
               + "; inspect the MS1 data (centroided MS1 present in the mzML?) and Sage's RT "
                 "alignment before trusting the quantities")
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


def not_applicable(rec=None, engine=None, why=None, ran=None, de_refusal=None):
    """The record of a search without quantification: identifications only, so there are no
    quantities to check, and no DE (gate_main refuses one). The ONE definition, for every engine:
    check() and run_search.py's Sage identifications-only route write it as Sage's (the defaults,
    unchanged); run_search.py's FragPipe route writes it for a workflow that runs no
    quantification (L11, Ruling 73a) with its own `engine`, `why`, `ran` and the DE refusal
    (`de_refusal`, which gate_main says). Its file keeps the name RECORD: run_de.R asks the gate
    about that file."""
    rec = dict(rec or {})
    why = why or NOT_APPLICABLE_WHY
    rec.update(status="not_applicable", why=why,
               message=(f"{ran or 'Sage ran without LFQ (quant.lfq false)'}: {why}, so there "
                        f"are no MS1 quantities to check and no DE. Its identifications are in "
                        f"{IDENTIFICATIONS}."))
    if engine:
        rec["engine"] = engine
    if de_refusal:
        rec["de_refusal"] = de_refusal
    return rec


def identifications_only(rec):
    """THE predicate "this search did not quantify", on a search's quantification record (the
    not_applicable() shape, any engine): the DE gate (gate_main) refuses on it, and every client
    text that says whether the proteins were quantified asks it through search_quantified()."""
    return isinstance(rec, dict) and rec.get("status") == "not_applicable"


def search_quantified(out, record=None):
    """True when the search in folder `out` quantified, False when it was for identifications
    only: identifications_only() of its record -- <out>/RECORD, the file the DE gate reads, else
    the copies search_provenance.json carries (write_record: `sage_lfq_check` for Sage,
    `quantification` for FragPipe, and Sage's `quantified: false`), as a delivered or laptop
    copy of the folder holds only search_provenance.json. A folder with neither (a DIA-NN
    search, which always quantifies) is True. The delivery README and email, the session
    README/AGENTS.md, the report page, the analysis brief and the reproducibility bundle
    (provenance.identified_only) all ask this one question (C2 B1). `record`: the search record
    the caller has already read (provenance.py --search-record), asked in place of
    <out>/search_provenance.json; None reads that file."""
    if identifications_only(load(out)):
        return False
    prov = record if isinstance(record, dict) else _load_json(os.path.join(out, PROVENANCE))
    if isinstance(prov, dict) and (prov.get("quantified") is False or any(
            identifications_only(prov.get(k)) for k in ("sage_lfq_check", "quantification"))):
        return False
    return True


def _unchecked(rec, why):
    rec.update(status="unchecked", reason=why,
               message=f"The Sage LFQ mass-window check could not run: {why}. Look for "
                       f"'target MS1 peaks at 5% FDR' in Sage's log by hand; 0 there means the "
                       f"search has no usable quantities.")
    return rec


def quantified_coverage(adapt):
    """How much of what Sage identified the adapter's report.parquet quantifies (C2 proteomics
    M2), from sage_adapt.json -- the one place both counts are made (run_search.adapt_sage):
    `target_precursors_kept` (the precursors in report.parquet) against `target_peptides_1pct`
    (the target peptides Sage identified at 1% FDR). Precursors against peptides is generous to
    the search (one peptide can be several precursors). -> {"target_precursors_quantified",
    "target_peptides_identified_1pct", "protein_groups_quantified", "fraction", "low",
    "low_below"}, or None when either count is not recorded. `protein_groups_quantified`: the
    adapter's `protein_groups_kept` (the groups report.parquet has amounts for), None in an older
    sage_adapt.json. `low`: under LOW_PEAK_FRACTION, the check's own "LFQ collapsed" line."""
    adapt = adapt if isinstance(adapt, dict) else {}
    kept, ids = adapt.get("target_precursors_kept"), adapt.get("target_peptides_1pct")
    if not isinstance(kept, int) or not isinstance(ids, int) or ids <= 0:
        return None
    groups = adapt.get("protein_groups_kept")
    frac = kept / ids
    return {"target_precursors_quantified": kept, "target_peptides_identified_1pct": ids,
            "protein_groups_quantified": groups if isinstance(groups, int) else None,
            "fraction": round(frac, 4), "low": frac < LOW_PEAK_FRACTION,
            "low_below": LOW_PEAK_FRACTION,
            "source": "sage_adapt.json: target_precursors_kept / target_peptides_1pct "
                      "(protein_groups_kept: the groups with amounts)"}


def coverage_note(cov):
    """The sentence AUDIT.md and the record read for a low coverage: the numbers, no verdict on
    the cause, and no task -- AUDIT.md goes into the client's bundle (COVERAGE_STAFF_TASK is said
    to staff only, by add_coverage)."""
    g = cov.get("protein_groups_quantified")
    return (f"Only {cov['target_precursors_quantified']:,} target precursors"
            + (f" in {g:,} protein groups" if isinstance(g, int) else "")
            + f" were quantified (report.parquet) for the "
            f"{cov['target_peptides_identified_1pct']:,} target peptides Sage identified at 1% FDR "
            f"({cov['fraction'] * 100:.1f}%, under {cov['low_below'] * 100:g}%): a differential "
            f"expression on this report sees a small part of what was identified. This is a "
            f"warning, not a refusal.")


def coverage_plain(cov):
    """What a low coverage means for the results, in the client's words (the audit table): most
    identified proteins get no amount at all, and those that do rest on few peptides."""
    g = cov.get("protein_groups_quantified")
    k, m = cov.get("target_precursors_quantified"), cov.get("target_peptides_identified_1pct")
    if isinstance(g, int):
        return (f"The quantification settings fit these runs, but only a small share of what was "
                f"identified could be quantified: the comparison has amounts for {g} protein "
                f"groups, from {k} quantified peptide signals for the {m} peptides identified. "
                f"Many identified proteins therefore have no amount in the comparison, and the "
                f"{g} that do rest on few peptides; the identifications are not affected.")
    return (f"The quantification settings fit these runs, but only a small share of the "
            f"identified peptides could be quantified ({k} quantified peptide signals for {m} "
            f"identified peptides), so many identified proteins have no amount in the "
            f"comparison, and those that do rest on few peptides; the identifications are not "
            f"affected.")


# The staff's next step for a low coverage: said where staff read it (add_coverage's stderr, the
# --adapt-only print), never in AUDIT.md or the client's table.
COVERAGE_STAFF_TASK = ("before relying on this report, compare its protein count with another "
                       "engine's search of the same runs (DIA-NN --dda, FragPipe)")


def add_coverage(out, adapt):
    """After the adapter wrote sage_adapt.json: put quantified_coverage() into this search's
    LFQ check record (and search_provenance.json), and WARN -- status warn, said -- when it is
    low. Never a gate: gated() still asks only about runs outside the window. -> the record, or
    None when there is no checked record or no counts to compare."""
    rec = load(out)
    cov = quantified_coverage(adapt)
    if not rec or cov is None or rec.get("status") not in ("ok", "warn"):
        return None
    rec["quantified_coverage"] = cov
    if cov["low"] and not rec.get("coverage_warning"):
        rec["coverage_warning"] = coverage_note(cov)
        rec["message"] = ((rec.get("message") or "") + " " + rec["coverage_warning"]).strip()
        rec["status"] = "warn"
    record(out, rec)
    if cov["low"]:
        say(rec)
        print(f"{TAG} STAFF: {COVERAGE_STAFF_TASK}.", file=sys.stderr, flush=True)
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
    """(status, message, detail, plain) for audit_results.py (its add(): `plain` is the client's
    words when they differ from the check's usual ones, else None), or None. ONE mapping: warn ->
    WARN, ok -> PASS, unchecked -> INFO (a check that could not run never changes the verdict).
    A warning ONLY for a low quantified coverage (C2 proteomics M2) carries its own plain words:
    the usual ones say the settings did not fit the mass error, which is then not so."""
    if not rec or rec.get("status") not in ("warn", "ok", "unchecked"):
        return None
    status = {"warn": "WARN", "ok": "PASS", "unchecked": "INFO"}[rec["status"]]
    detail = {k: rec.get(k) for k in ("ppm_tolerance", "ppm_tolerance_source",
                                      "target_ms1_peaks_5pct_fdr", "target_peptides_1pct",
                                      "runs_outside_lfq_window", "suggested_ppm_tolerance",
                                      "corrected_params", "rerun", "accepted",
                                      "quantified_coverage")
              if rec.get(k) is not None}
    acc = dict(rec.get("accepted") or {})
    if acc.get("by"):
        # AUDIT.md is delivered: the role, never the login (C2 security M1/L5; staff.py's rule)
        import staff
        if acc["by"] not in (staff.ROLE, staff.USER_ROLE):
            acc["by"] = staff.role_for(acc["by"])
            acc["by_source"] = staff.ACCEPTED_BY_SOURCE
        detail["accepted"] = acc
    who = acc.get("by") or "who: not recorded"
    cov = rec.get("quantified_coverage") if isinstance(rec.get("quantified_coverage"), dict) else {}
    plain = None
    if (rec.get("status") == "warn" and rec.get("coverage_warning")
            and not rec.get("runs_outside_lfq_window") and not rec.get("low_ms1_peaks")):
        plain = coverage_plain(cov)
    return status, (rec.get("message") or "") + (
        f" The user accepted these quantities as they are ({who}, {acc.get('at')}): "
        f"{acc.get('reason')}." if acc else ""), detail, plain


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
    rec = carry_acceptance(load(a.out), check(a.out, a.params, a.log))
    try:
        rec["written_to"] = record(a.out, rec)
    except OSError as e:
        print(f"{TAG} NOTE: could not write {RECORD} in {a.out}: {e}", file=sys.stderr)
    say(rec)
    msg = refusal(a.out, rec)
    if msg:
        # the job does not fail -- Sage did its part -- but the next step will refuse, so say how
        print(f"{TAG} {msg}", file=sys.stderr, flush=True)
    print(json.dumps(rec, indent=2))
    return 0


def gate_main(argv):
    """`sage_lfq_check.py gate --out <dir>` -> {"checked", "refused", "message"} on stdout, exit 0:
    run_de.R's question (the gate's one rule is gated()). Reads the record; re-checks nothing."""
    ap = argparse.ArgumentParser(prog="sage_lfq_check.py gate")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    rec = load(a.out)
    msg = (rec.get("refusal") or refusal(a.out, rec)) if gated(rec) else None
    if identifications_only(rec):
        # identifications only: there is nothing to run a DE on (run_de.R stops with this) --
        # the engine's own words when its record carries them (a FragPipe workflow without
        # quantification), else Sage's
        msg = rec.get("de_refusal") or (
            f"REFUSED: this Sage search was for {rec.get('why') or NOT_APPLICABLE_WHY}: there "
            f"are no quantities, so no DE. Its identifications are in {IDENTIFICATIONS} "
            f"(run_search.py --adapt-only, with the cfg the search ran with, writes it).")
    print(json.dumps({"checked": rec is not None, "refused": bool(msg), "message": msg}))
    return 0


if __name__ == "__main__":
    sys.exit(gate_main(sys.argv[2:]) if sys.argv[1:2] == ["gate"] else main())
