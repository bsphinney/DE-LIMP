#!/usr/bin/env python3
"""
normalization_check.py -- is cross-run normalisation right for THIS experiment? The safety net
under the experiment type's default (experiment_type.py), whatever the user said.

Global normalisation -- DIA-NN's included -- assumes most proteins do not change between
samples. An IP / pull-down, separation fractions, and sometimes proximity labelling or a
secretome, violate
that: the controls carry little protein, normalisation scales them up several-fold, and the
enrichment is hidden or turned into depletion. A TurboID cohort was being analysed when this was
written. The type gives a DEFAULT (experiment_type.DEFAULTS); this check, run for EVERY type,
decides whether to stop and ask:

  B1 factor_confound  the per-run normalisation factor (log2 Precursor.Normalised /
                      Precursor.Quantity, DIA-NN: normalised = Normalisation.Factor x
                      non-normalised) against the groups: the gap between group medians, whether
                      the groups separate completely, and an exact (or seeded) permutation p.
                      The deciding signal: LARGE (>= FACTOR_GAP_LOG2) AND confounded with group
                      (the groups separate completely) -- e.g. the controls scaled up several-fold.
  B2 raw_shift        before normalisation (Precursor.Quantity), each simple contrast's
                      per-protein log2 shift: centred on zero, or most proteins moved one way?
  B3 id_gap           identified precursors per run, group medians: very different depth.
  B4 instability      DIA-NN's Normalisation.Instability (check_report_runs.py, the one reader).
  V  volcano          from the DE tables behind the volcano, for each contrast, on BOTH the
                      normalised and the non-normalised DE (run_de.R --quantities): the centre
                      (median log2FC of the non-significant proteins), the up:down ratio, the
                      fraction significant, the bait's rank and log2FC, and the p-value
                      histogram's shape. Reported always, so staff learn what normal looks like.
  carboxylases        proximity labelling only: how constant the endogenously biotinylated
                      carboxylases are under each option -- a reference, never a normaliser.

WHAT TRIPS depends on the default (one rule, trips()):
  default normalised  B1 large+confounded, B2 lopsided, B3 gap, B4, or a volcano anomaly in the
                      normalised DE -- the data say the normalisation assumption fails.
  default raw         run factors spanning >= FACTOR_GAP_LOG2 that the groups do NOT explain
                      (loading variation the raw quantities would keep), or a volcano anomaly in
                      the non-normalised DE.
A trip STOPS the analysis before the final DE (exit 3): NORMALIZATION_CHECK.md holds the
side-by-side summary and a recommendation, the user chooses (`decide`), and the decision is
recorded. Nothing switches silently in either direction: run_de.R --normalization-check refuses
any other quantities (`gate`). The record is written whether or not anything trips, with fixed
keys and a schema_version, so projects can be compared over months (does TurboID really behave
differently from a classic IP?); run_de.R copies it into de_provenance.json
"normalization_check".

  # the two DEs the check compares (scratch folders), then the check
  Rscript run_de.R --input report.parquet --metadata conditions.csv --quantities normalised --outdir <S>/output/norm_check/normalised
  Rscript run_de.R --input report.parquet --metadata conditions.csv --quantities raw        --outdir <S>/output/norm_check/raw
  python3 normalization_check.py run --report report.parquet --conditions conditions.csv \\
      --de-normalised <S>/output/norm_check/normalised --de-raw <S>/output/norm_check/raw \\
      --session <S> [--bait GENE] [--out <S>/output/norm_check]
      # the bait: --bait, else the recorded one, else the ONE --add-fasta sequence (unconfirmed)
  # tripped (exit 3): show NORMALIZATION_CHECK.md, ask, then record THEIR answer -- the decision
  # is a Core staff member's, never Claude's: --by is their HIVE login (staff.py checks it)
  python3 normalization_check.py decide --check <S>/output/norm_check/normalization_check.json \\
      --quantities raw|normalised --by <their HIVE login> --reason "<why, in their words>"
  # the final DE, held to the decision
  Rscript run_de.R ... --quantities <decided> --normalization-check <the json> --outdir <S>/output/tables
  # a DE from before 2.10 (no check record) that staff say stands: lets the delivery through
  python3 normalization_check.py ack-legacy --session <S> --by <their HIVE login> --note "<why>"

MaxLFQ has no non-normalised protein quantity in a normal report (PG.MaxLFQ is normalised unless
the search ran with --no-norm): its raw DE reads a report searched with --no-norm, given to the
check as --report-raw; the B1-B4 checks still read --report. dpc needs no second search.

The client documents name who decided by role (staff.ROLE, "Core staff"); the name and login are
in the staff-only STAFF_RECORD (logs/normalization_check.staff.json).

Stdlib, plus pyarrow to read a parquet report (without it B1-B3 are "not assessed", said so).
"""
import argparse
import csv
import datetime
import hashlib
import itertools
import json
import math
import os
import random
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import experiment_type as et  # noqa: E402  the type, its default, the pull-down design
import staff  # noqa: E402  who may decide (--by), and the role client documents name instead
from check_report_runs import run_name, normalisation_instability  # noqa: E402

SCHEMA, SCHEMA_VERSION = "normalization_check", 1
RECORD, SUMMARY = "normalization_check.json", "NORMALIZATION_CHECK.md"
# Who decided, by name and login (decide, ack-legacy): staff-only, in the session's logs/ (or
# beside the check outside a session). The check record -- copied into de_provenance.json, the
# Methods and the zip -- names the role (staff.ROLE) and points here.
STAFF_RECORD = "normalization_check" + staff.STAFF_SUFFIX
NOTE_MAX = 300                   # ack-legacy --note: one line, as qc_bracket.py ack's
TRIPPED_EXIT = 3

# ---- thresholds: each one's basis, set on the synthetic cases in
# tests/test_normalization_check.py and on what whole-proteome loading looks like ----------------
# B1: a run's normalisation factor corrects its loading. In a whole-proteome cohort loading
# varies by tens of percent and does not follow the groups; the IP signature is the controls
# scaled up several-fold, all of them. 2-fold (1 log2) between group medians, with the groups
# separating completely, is far outside loading noise and well inside that signature. Its false-
# trip rate grows as the groups shrink: the 2.10 review simulated whole-proteome cohorts of 2 vs 2
# runs and found 0.1% / 4% / 13% tripping at a per-run loading SD of 0.3 / 0.5 / 0.7 log2 (two
# runs per group separate by chance a third of the time). A trip only asks -- it never switches --
# so B1 stays; with n = 2 per group expect to be asked more often.
FACTOR_GAP_LOG2 = 1.0
# B2: the same 2-fold, as the median per-protein shift before normalisation, with at least 80%
# of the proteins moved that way -- a loading difference or an enrichment, not biology of a few.
RAW_SHIFT_LOG2 = 1.0
RAW_SHIFT_FRACTION = 0.80
# B3: groups of one proteome identify nearly the same precursors; a 1.5-fold gap in the median
# means the groups do not share most of their proteome -- the assumption normalisation needs.
ID_GAP_RATIO = 1.5
# V: after a sound normalisation the non-significant majority sits at log2FC ~0; 0.5 (1.41-fold)
# off is a shifted volcano. More than 5:1 up:down (or down:up) among >= 20 significant proteins is
# lopsided -- for an enrichment-vs-control contrast only the DOWN side counts (enrichment is up by
# design; many proteins DOWN there is background made to look depleted by scaled-up controls).
# Over half the proteins significant is implausible (audit_results.py's de_signal says the same).
CENTRE_LOG2 = 0.5
RATIO = 5.0
MIN_SIG_FOR_RATIO = 20
FRAC_SIG = 0.5
BAIT_TOP = 10          # the bait should be among the 10 most enriched, and significant
# The recommendation keeps the experiment type's default unless that default's own volcano has at
# least this many MORE anomalies than the alternative's: one anomaly can sit on a threshold, two
# more is a clear difference.
CLEARLY_WORSE = 2
# p-value histogram, 20 bins: anti-conservative = the first bin >= 2x the median density of the
# upper half; odd = the last bin >= 2x it (excess near 1: a conservative or broken test).
P_BINS, P_SPIKE = 20, 2.0
PERMUTATIONS = 20000   # B1's p: exact up to this many labellings, else this many, seeded
Q_MAX = 0.01           # the identification filter for B1-B3's precursors (Q.Value)


def _median(xs):
    xs = [x for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x))]
    return statistics.median(xs) if xs else None


def _r(x, d=3):
    return None if x is None else round(float(x), d)


def _fmt(x, spec):
    """x formatted, or n/a when there is no value (e.g. no non-significant protein left)."""
    return "n/a" if x is None else format(x, spec)


def read_conditions(path):
    """{run name: group} from conditions.csv (File.Name, Group)."""
    with open(path, newline="", encoding="utf-8") as fh:
        return {run_name(r["File.Name"]): r["Group"].strip() for r in csv.DictReader(fh)
                if (r.get("File.Name") or "").strip() and (r.get("Group") or "").strip()}


# ---------------------------------------------------------------- the report (B1-B3) --
def read_report(report, runs):
    """The report's per-run and per-protein numbers B1-B3 need, or (None, why)."""
    try:
        import pyarrow as pa
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
    except ImportError:
        return None, "pyarrow is not installed in this Python, so the report could not be read"
    want = ["Run", "Precursor.Id", "Protein.Group", "Precursor.Quantity", "Precursor.Normalised",
            "Q.Value"]
    try:
        if report.endswith(".parquet"):
            have = pq.read_schema(report).names
            t = pq.read_table(report, columns=[c for c in want if c in have])
        else:
            import pyarrow.csv as pcsv
            t = pcsv.read_csv(report, parse_options=pcsv.ParseOptions(delimiter="\t"),
                              convert_options=pcsv.ConvertOptions(include_columns=want,
                                                                  include_missing_columns=True))
            have = [c for c in want if t.column(c).null_count < t.num_rows]
    except (OSError, ValueError) as e:
        return None, f"{report} could not be read ({type(e).__name__}: {e})"
    missing = [c for c in want if c not in have]
    if missing:
        return None, (f"{os.path.basename(report)} has no {', '.join(missing)}: not a DIA-NN "
                      "precursor report, so its normalisation cannot be read")
    t = t.filter(pc.and_(pc.less_equal(t["Q.Value"], Q_MAX),
                         pc.and_(pc.greater(t["Precursor.Quantity"], 0),
                                 pc.greater(t["Precursor.Normalised"], 0))))
    names = [run_name(r) for r in t["Run"].to_pylist()]
    keep = pa.array([n in runs for n in names])
    t = t.filter(keep).append_column("run", pa.array([n for n in names if n in runs]))
    lq = pc.log2(pc.cast(t["Precursor.Quantity"], pa.float64()))
    lf = pc.subtract(pc.log2(pc.cast(t["Precursor.Normalised"], pa.float64())), lq)
    t = t.append_column("lq", lq).append_column("lf", lf)
    per_run = t.group_by("run").aggregate([("lf", "approximate_median"),
                                          ("Precursor.Id", "count_distinct")]).to_pydict()
    factors = dict(zip(per_run["run"], per_run["lf_approximate_median"]))
    ids = dict(zip(per_run["run"], per_run["Precursor.Id_count_distinct"]))
    pr = t.group_by(["Protein.Group", "run"]).aggregate([("lq", "approximate_median")]).to_pydict()
    protein_run = {}
    for p, r, v in zip(pr["Protein.Group"], pr["run"], pr["lq_approximate_median"]):
        protein_run.setdefault(p, {})[r] = v
    return {"factors": factors, "ids": ids, "protein_run": protein_run,
            "factor_basis": "log2(Precursor.Normalised / Precursor.Quantity) per precursor (DIA-NN: "
                            "normalised = Normalisation.Factor x non-normalised), median per run "
                            f"(t-digest), precursors at Q.Value <= {Q_MAX}"}, None


def permutation_p(values, labels):
    """P(R^2 >= observed) over relabellings of the groups: exact when there are at most
    PERMUTATIONS distinct labellings, else PERMUTATIONS seeded draws."""
    def r2(lab):
        grand = statistics.mean(values)
        tot = sum((v - grand) ** 2 for v in values)
        if tot == 0:
            return 0.0
        groups = {}
        for v, g in zip(values, lab):
            groups.setdefault(g, []).append(v)
        between = sum(len(vs) * (statistics.mean(vs) - grand) ** 2 for vs in groups.values())
        return between / tot
    obs = r2(labels)
    n = len(values)
    sizes = [labels.count(g) for g in dict.fromkeys(labels)]
    total = math.factorial(n)
    for k in sizes:
        total //= math.factorial(k)
    if total <= PERMUTATIONS:
        hits = cnt = 0
        uniq = list(dict.fromkeys(labels))

        def assign(idx, remaining, cur):
            nonlocal hits, cnt
            if idx == len(uniq) - 1:
                lab = list(cur)
                for i in remaining:
                    lab[i] = uniq[idx]
                cnt += 1
                hits += r2(lab) >= obs - 1e-12
                return
            for comb in itertools.combinations(remaining, sizes[idx]):
                lab = list(cur)
                for i in comb:
                    lab[i] = uniq[idx]
                assign(idx + 1, [i for i in remaining if i not in comb], lab)
        assign(0, list(range(n)), [None] * n)
        return obs, hits / cnt, "exact"
    rng = random.Random(20261001)
    hits = sum(r2(rng.sample(labels, n)) >= obs - 1e-12 for _ in range(PERMUTATIONS))
    return obs, (hits + 1) / (PERMUTATIONS + 1), f"{PERMUTATIONS} seeded permutations"


def factor_confound(data, groups):
    if data is None:
        return {"status": "not_assessed"}
    by = {}
    for r, f in data["factors"].items():
        if f is not None and r in groups:
            by.setdefault(groups[r], []).append(f)
    by = {g: v for g, v in by.items() if v}
    if len(by) < 2:
        return {"status": "not_assessed", "why": "fewer than two groups with factors"}
    per = {g: {"median_log2_factor": _r(_median(v)), "min": _r(min(v)), "max": _r(max(v)),
               "n_runs": len(v)} for g, v in sorted(by.items())}
    hi = max(per, key=lambda g: per[g]["median_log2_factor"])
    lo = min(per, key=lambda g: per[g]["median_log2_factor"])
    gap = per[hi]["median_log2_factor"] - per[lo]["median_log2_factor"]
    separated = min(by[hi]) > max(by[lo])
    vals = [f for g in sorted(by) for f in by[g]]
    labs = [g for g in sorted(by) for _ in by[g]]
    r2, p, how = permutation_p(vals, labs)
    spread = max(vals) - min(vals)
    return {"per_group": per, "highest": hi, "lowest": lo, "gap_log2": _r(gap),
            "gap_fold": _r(2 ** gap, 2), "separated": separated, "r2": _r(r2),
            "permutation_p": _r(p, 4), "permutation": how, "run_spread_log2": _r(spread),
            "large_and_confounded": gap >= FACTOR_GAP_LOG2 and separated,
            "large_unexplained": spread >= FACTOR_GAP_LOG2 and not (gap >= FACTOR_GAP_LOG2
                                                                    and separated),
            "basis": data["factor_basis"],
            "threshold": f"gap >= {FACTOR_GAP_LOG2:g} log2 ({2 ** FACTOR_GAP_LOG2:g}-fold) between "
                         "group medians with the two groups separating completely"}


def simple_contrasts(contrasts, groups):
    """[(contrast, positive group, negative group)] for the contrasts of the form A-B."""
    gs = set(groups.values())
    out = []
    for c in contrasts:
        m = re.fullmatch(r"\s*([^+\-*/()]+?)\s*-\s*([^+\-*/()]+?)\s*", c)
        if m and m.group(1) in gs and m.group(2) in gs:
            out.append((c, m.group(1), m.group(2)))
    return out


def raw_shift(data, groups, contrasts):
    if data is None:
        return {"status": "not_assessed"}
    out = {}
    for c, a, b in simple_contrasts(contrasts, groups):
        shifts = []
        for p, rv in data["protein_run"].items():
            va = [v for r, v in rv.items() if groups.get(r) == a]
            vb = [v for r, v in rv.items() if groups.get(r) == b]
            na = sum(1 for g in groups.values() if g == a)
            nb = sum(1 for g in groups.values() if g == b)
            if len(va) * 2 >= na and len(vb) * 2 >= nb and va and vb:
                shifts.append(statistics.mean(va) - statistics.mean(vb))
        if not shifts:
            continue
        med = _median(shifts)
        frac_up = sum(1 for x in shifts if x > 0) / len(shifts)
        out[c] = {"median_log2_shift": _r(med), "fraction_up": _r(frac_up),
                  "n_proteins": len(shifts),
                  "lopsided": abs(med) >= RAW_SHIFT_LOG2 and
                              max(frac_up, 1 - frac_up) >= RAW_SHIFT_FRACTION}
    if not out:
        return {"status": "not_assessed", "why": "no contrast of the form A-B between two groups"}
    return {"per_contrast": out,
            "basis": "per protein: mean over each group's runs of the median log2 "
                     "Precursor.Quantity (non-normalised), proteins seen in at least half the "
                     "runs of both groups",
            "threshold": f"|median shift| >= {RAW_SHIFT_LOG2:g} log2 with >= "
                         f"{RAW_SHIFT_FRACTION:.0%} of proteins shifted the same way"}


def id_gap(data, groups):
    if data is None:
        return {"status": "not_assessed"}
    by = {}
    for r, n in data["ids"].items():
        if r in groups:
            by.setdefault(groups[r], []).append(n)
    per = {g: {"median_precursors": _median(v), "n_runs": len(v)} for g, v in sorted(by.items())}
    meds = [v["median_precursors"] for v in per.values() if v["median_precursors"]]
    if len(meds) < 2:
        return {"status": "not_assessed", "why": "fewer than two groups"}
    ratio = max(meds) / min(meds)
    return {"per_group": per, "ratio": _r(ratio, 2), "gap": ratio >= ID_GAP_RATIO,
            "basis": f"distinct precursors per run at Q.Value <= {Q_MAX}",
            "threshold": f"highest / lowest group median >= {ID_GAP_RATIO:g}"}


def instability(report, runs):
    ni = normalisation_instability(report, list(runs))
    if not ni or not ni.get("column") or not ni.get("per_run"):
        return {"status": "not_assessed", "why": "no Normalisation.Instability for these runs"}
    return {"median": _r(ni["median"]), "max": _r(ni["max"]), "flagged": ni["flagged"],
            "threshold": ni["threshold"], "basis": ni["basis"]}


# ---------------------------------------------------------------------- the volcano --
def read_de(de_dir):
    """(de_provenance.json, {contrast: rows}) of one run_de.R output, or (None, why)."""
    try:
        with open(os.path.join(de_dir, "de_provenance.json"), encoding="utf-8") as fh:
            prov = json.load(fh)
    except (OSError, ValueError) as e:
        return None, f"no readable de_provenance.json in {de_dir} ({e})"
    tables = {}
    for c, t in (prov.get("de_tables") or {}).items():
        path = os.path.join(de_dir, t.get("file") or "")
        try:
            with open(path, newline="", encoding="utf-8") as fh:
                tables[c] = list(csv.DictReader(fh))
        except OSError as e:
            return None, f"{path} could not be read ({e})"
    return {"prov": prov, "tables": tables}, None


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def _tokens(cell):
    return {t.strip().upper() for t in re.split(r"[;,]", cell or "") if t.strip()}


def pvalue_shape(ps):
    if len(ps) < 50:
        return {"shape": "not_assessed", "why": "fewer than 50 p-values"}
    bins = [0] * P_BINS
    for p in ps:
        bins[min(int(p * P_BINS), P_BINS - 1)] += 1
    upper = statistics.median(bins[P_BINS // 2:]) or 1
    first, last = bins[0] / upper, bins[-1] / upper
    shape = ("odd: excess near 1" if last >= P_SPIKE else
             "anti-conservative" if first >= P_SPIKE else "uniform")
    return {"shape": shape, "first_bin_ratio": _r(first, 2), "last_bin_ratio": _r(last, 2),
            "bins": bins}


def volcano(rows, adjp, enrichment, bait):
    lfc, sig, ps, kept = [], [], [], []
    for r in rows:
        l, a, p = _f(r.get("logFC")), _f(r.get("adj.P.Val")), _f(r.get("P.Value"))
        if l is None or a is None:
            continue
        kept.append(r)
        lfc.append(l)
        sig.append(a < adjp)
        if p is not None:
            ps.append(p)
    if not lfc:
        return {"status": "not_assessed", "why": "no tested proteins"}
    up = sum(1 for l, s in zip(lfc, sig) if s and l > 0)
    down = sum(1 for l, s in zip(lfc, sig) if s and l < 0)
    nonsig = [l for l, s in zip(lfc, sig) if not s]
    m = {"n": len(lfc), "n_sig": up + down, "up": up, "down": down,
         "up_down_ratio": _r(up / down, 2) if down else None,
         "frac_sig": _r((up + down) / len(lfc)),
         "centre_all": _r(_median(lfc)), "centre_nonsig": _r(_median(nonsig)),
         "enrichment_contrast": enrichment, "pvalue": pvalue_shape(ps)}
    if bait:
        order = sorted(((l, i) for i, l in enumerate(lfc)), reverse=True)
        rank = {i: k + 1 for k, (_, i) in enumerate(order)}
        b = bait.upper()
        hit = [i for i, r in enumerate(kept)
               if b in _tokens(r.get("Genes")) | _tokens(r.get("Protein.Group"))]
        m["bait"] = ({"name": bait, "rank": rank[hit[0]], "log2fc": _r(lfc[hit[0]]),
                      "significant": sig[hit[0]]} if hit else {"name": bait, "rank": None,
                                                               "found": False})
    a = []
    if m["centre_nonsig"] is not None and abs(m["centre_nonsig"]) >= CENTRE_LOG2:
        a.append(f"centre shift {m['centre_nonsig']:+.2f} log2 (non-significant proteins)")
    if m["n_sig"] >= MIN_SIG_FOR_RATIO:
        if enrichment:
            if down >= MIN_SIG_FOR_RATIO and down > up:
                a.append(f"{down} significantly DOWN vs {up} up in an enrichment-vs-control "
                         "contrast: background made to look depleted")
        elif max(up, down) >= RATIO * max(min(up, down), 1):
            a.append(f"lopsided: {up} up vs {down} down")
    if m["frac_sig"] > FRAC_SIG:
        a.append(f"{m['frac_sig']:.0%} of proteins significant")
    if bait and enrichment:
        bt = m["bait"]
        if not bt.get("rank") or bt["rank"] > BAIT_TOP or not bt.get("significant") \
                or (bt.get("log2fc") or 0) <= 0:
            a.append(f"the bait {bait} is not among the {BAIT_TOP} most enriched and significant"
                     + (f" (rank {bt['rank']}, log2FC {bt['log2fc']:+.2f})" if bt.get("rank")
                        else " (not found)"))
    if m["pvalue"]["shape"].startswith("odd"):
        a.append("p-value histogram: excess near 1")
    m["anomalies"] = a
    return m


def volcano_version(de, design, bait, raw=False, etype=None):
    if de is None:
        return None
    adjp = _f(de["prov"].get("adjp")) or 0.05
    out = {}
    for c, rows in de["tables"].items():
        enrich = c in design["vs_control"]
        v = volcano(rows, adjp, enrich, bait)
        if raw and etype in RAW_VOLCANO_EXPECTED:
            # separation / profiling fractions: each fraction holds a different proteome, so
            # without normalisation a shifted, lopsided, mostly-significant volcano is the design
            # (DIA-NN's README: no normalisation for such fractions) -- reported, not counted
            v["expected"], v["anomalies"] = v.get("anomalies", []), []
        elif raw and enrich:
            # without normalisation an enrichment contrast's background sits above zero and much
            # of it tests up (the controls carry less protein): expected there, not an anomaly
            v["anomalies"] = [x for x in v.get("anomalies", [])
                              if not (x.startswith("centre shift +") or "of proteins significant" in x)]
        out[c] = v
    return {"adjp": adjp, "contrasts": out}


# Types whose non-normalised volcano is SUPPOSED to look anomalous (see volcano_version).
RAW_VOLCANO_EXPECTED = ("fractions_separation",)


def carboxylase_stability(de_dir):
    """Per carboxylase found, the SD of its log2 values across samples (Expression_Matrix.csv)."""
    try:
        with open(os.path.join(de_dir, "Expression_Matrix.csv"), newline="",
                  encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
    except OSError:
        return None
    if not rows:
        return None
    ids = [c for c in ("Protein.Group", "Genes", "Protein.Names") if c in rows[0]]
    samples = [c for c in rows[0] if c not in ids]
    found = {}
    for r in rows:
        genes = _tokens(r.get("Genes"))
        for g in et.CARBOXYLASES:
            if g in genes and g not in found:
                vals = [v for v in (_f(r.get(s)) for s in samples) if v is not None]
                if len(vals) >= 2:
                    found[g] = _r(statistics.stdev(vals))
    return {"sd_log2": found, "median_sd_log2": _r(_median(list(found.values()))),
            "n_found": len(found)} if found else {"n_found": 0}


def sha256(path):
    """sha256 of a file (a report can be GBs: read in chunks), or None when it cannot be read."""
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return None
    return h.hexdigest()


def check_dir(session):
    return os.path.join(os.path.abspath(session), "output", "norm_check")


def step8_commands(report, conditions, session, method="dpc", scripts=None, report_raw=None):
    """Step 8 for a session, as commands: the two DEs the check compares (--check-input), the
    check, and the final DE held to its decision. THE one wording -- the chain's recovery step
    (diann_parallel.py --next) and run_de.R's refusal both print it. -> {"check", "final",
    "next"}: `next` runs the first three (the final DE waits for the check's answer).
    `report_raw`: MaxLFQ's non-normalised DE reads a report searched with --no-norm (the module
    docstring); the check is told so with --report-raw."""
    import shlex
    q = shlex.quote
    sc = scripts or HERE
    nd = check_dir(session)

    def de(rep):
        return (f"Rscript {q(os.path.join(sc, 'run_de.R'))} --input {q(rep)} --metadata "
                f"{q(conditions)} --method {method}")
    check = [f"{de(report_raw if k == et.RAW and report_raw else report)} --quantities {k} "
             f"--check-input --outdir {q(os.path.join(nd, k))}" for k in (et.NORMALISED, et.RAW)]
    check.append(f"python3 {q(os.path.join(sc, 'normalization_check.py'))} run --report "
                 f"{q(report)} --conditions {q(conditions)} --de-normalised "
                 f"{q(os.path.join(nd, et.NORMALISED))} --de-raw {q(os.path.join(nd, et.RAW))} "
                 f"--session {q(session)}" + (f" --report-raw {q(report_raw)}" if report_raw else ""))
    final = (f"{de(report)} --quantities <the decided quantities> --normalization-check "
             f"{q(os.path.join(nd, RECORD))} --outdir "
             f"{q(os.path.join(os.path.abspath(session), 'output', 'tables'))}"
             + (f"  # raw decided: --input {q(report_raw)}" if report_raw else ""))
    nxt = (" && ".join(check) + "  # exit 3: stop, show NORMALIZATION_CHECK.md, ask (SKILL.md "
           "step 8); then the final DE: " + final)
    return {"check": check, "final": final, "next": nxt}


MAXLFQ_RAW_NOTE = (
    "\n(maxlfq: a normal report has no non-normalised protein quantity -- PG.MaxLFQ is normalised "
    "unless the search ran with --no-norm -- so the raw DE above is refused on it. Either run all "
    "three DEs with --method dpc, which reads both quantities from this report; or, to keep "
    "MaxLFQ, search again with --no-norm (diann_parallel.py --no-norm writes "
    "no_norm_report.parquet), run the raw DE on that report, and pass it to the check as "
    "--report-raw. SKILL.md step 8.)")


def required(report, conditions, method="dpc"):
    """Must a DE of these files be held to a normalisation check? Yes when they belong to a session
    that records an experiment type (experiment_type.session_of): that type has a default the DE
    must not override silently. -> {"required", "session", "type", "default", "message"}."""
    session = et.session_of(conditions, report)
    rec = et.load(session) if session else None
    if not rec:
        return {"required": False, "session": session, "type": None, "default": None,
                "message": None}
    d = et.default_for(rec["type"])
    cmds = step8_commands(report, conditions, session, method)
    return {"required": True, "session": session, "type": rec["type"],
            "default": d["quantities"],
            "message": (f"REFUSED: this session records the experiment type {rec['label']} "
                        f"(default: {d['quantities']} quantities), and this DE was given no "
                        f"--normalization-check -- it would read whichever quantities it was "
                        f"told, unchecked. Run step 8's check first:\n  "
                        + "\n  ".join(cmds["check"])
                        + f"\nthen, once decided:\n  {cmds['final']}\n(The two DEs the check "
                          "compares pass --check-input.)" + (MAXLFQ_RAW_NOTE if method == "maxlfq"
                                                             else ""))}


def searched_no_norm(report, search_dir=None):
    """Did the DIA-NN search behind `report` run with --no-norm? From its provenance, never by
    reading the report (run_de.R reads the data only when this says it is not recorded): the
    chain's --no-norm report name, else the search's search_provenance.json -- the 5-step chain
    run_search.py runs never has --no-norm (it passes none and the chain strips the cfg's), a
    single-shot search has it when its cfg does. `search_dir`: run_de.R's search folder (the
    report's, or its parent). -> {"no_norm": True | False | None, "source"}."""
    import diann_parallel as dp
    if os.path.basename(report) == dp.NO_NORM_REPORT:
        return {"no_norm": True, "source": f"the report's name ({dp.NO_NORM_REPORT}, what "
                                           "diann_parallel.py --no-norm writes)"}
    sp = os.path.join(search_dir or os.path.dirname(os.path.abspath(report)),
                      "search_provenance.json")
    try:
        with open(sp, encoding="utf-8") as fh:
            prov = json.load(fh)
    except (OSError, ValueError):
        return {"no_norm": None, "source": f"not recorded (no readable {sp})"}
    if not isinstance(prov, dict) or prov.get("engine") != "diann":
        return {"no_norm": None, "source": f"not recorded ({sp} is not a DIA-NN search's)"}
    mode = prov.get("search_mode")
    if mode == "parallel_5step":
        return {"no_norm": False, "source": f"{sp}: run_search.py's 5-step chain, which runs "
                                            "without --no-norm"}
    if mode == "single_shot":
        cfg = prov.get("resolved_params_file") or prov.get("params_file")
        if not cfg:
            return {"no_norm": None, "source": f"not recorded ({sp} names no cfg)"}
        try:
            has = "--no-norm" in dp.cfg_tokens(cfg)
        except (dp.CfgError, OSError) as e:
            return {"no_norm": None, "source": f"not recorded (the search's cfg could not be "
                                               f"read: {e})"}
        return {"no_norm": has, "source": f"{sp}: a single-shot search whose cfg ({cfg}) "
                                          + ("has --no-norm" if has else "has no --no-norm")}
    return {"no_norm": None, "source": f"not recorded ({sp} has search_mode {mode!r})"}


def de_provenance_path(session):
    return os.path.join(os.path.abspath(session), "output", "tables", "de_provenance.json")


def staff_record_path(session):
    return os.path.join(os.path.abspath(session), "logs", STAFF_RECORD)


def _read_prov(path):
    with open(path, encoding="utf-8") as fh:
        prov = json.load(fh)
    if not isinstance(prov, dict):
        raise ValueError("not a JSON object")
    return prov


def _legacy_acks(session, prov_sha):
    """The ack-legacy entries in the session's staff record that cover THIS DE record."""
    return [e for e in staff.read(staff_record_path(session)) if e.get("event") == "ack-legacy"
            and prov_sha and e.get("de_provenance_sha256") == prov_sha]


def delivery_gate(session):
    """Hold an analysis delivery whose DE did not run under a decided normalisation check -- the
    second gate core_submission.py deliver applies, beside qc_bracket.delivery_gate, in the same
    form: {"proceed", "reason", "status", "hint"}. A DE from before 2.10 (no check record) passes
    once Core staff acknowledge that it stands (ack-legacy), for that DE record as it is now."""
    hint = ("run step 8's normalisation check and the final DE with --normalization-check "
            "(SKILL.md step 8); python3 scripts/normalization_check.py required --input <report> "
            "--metadata <conditions> prints the commands")
    path = de_provenance_path(session)
    try:
        prov = _read_prov(path)
    except (OSError, ValueError) as e:
        return {"proceed": False, "status": "unreadable", "hint": hint,
                "reason": f"the DE record {path} could not be read ({e})"}
    nc = prov.get("normalization_check")
    status = nc.get("status") if isinstance(nc, dict) else None
    if status == "decided" and (nc.get("decision") or {}).get("quantities"):
        return {"proceed": True, "status": status, "hint": None,
                "reason": f"{nc['decision']['quantities']} quantities, decided by "
                          f"{nc['decision'].get('by')}"}
    if "normalization_check" not in prov:
        try:
            acks = _legacy_acks(session, sha256(path))
        except ValueError as e:
            return {"proceed": False, "status": "unreadable", "hint": hint,
                    "reason": f"the staff record {staff_record_path(session)} could not be read ({e})"}
        if acks:
            return {"proceed": True, "status": "legacy_acknowledged", "hint": None,
                    "reason": f"the DE predates the normalisation check (before skill 2.10); "
                              f"{staff.ROLE} acknowledged that it stands (logs/{STAFF_RECORD})"}
        hint = ("this DE predates the check (before skill 2.10): run step 8 on it (" + hint +
                "), or, if Core staff decide it stands as it is, record that: python3 "
                "scripts/normalization_check.py ack-legacy --session <S> --by <their HIVE login> "
                "--note \"<why this DE stands>\"")
    why = {None: "the DE record predates the normalisation check (no normalization_check)",
           "not_run": "the DE ran without the normalisation check",
           "check_input": "the delivered DE is one of the check's two inputs, not the final DE"}
    return {"proceed": False, "status": status, "hint": hint,
            "reason": why.get(status, f"the normalisation check's status is {status!r}") +
                      ": whether its quantities suit the experiment is not checked"}


def ack_legacy(session, by, note):
    """Core staff say a DE from before 2.10 -- one with no normalization_check record at all --
    stands as it is. Valid only for such a DE; recorded, by name and login, in the staff-only
    record, tied to the DE record's sha256 (a new DE needs a new acknowledgement). -> the gate."""
    note = (note or "").strip()
    if not note or "\n" in note or "\r" in note or len(note) > NOTE_MAX:
        raise ValueError(f"--note: one line, 1-{NOTE_MAX} characters, saying why this DE stands")
    path = de_provenance_path(session)
    try:
        prov = _read_prov(path)
    except (OSError, ValueError) as e:
        raise ValueError(f"no DE record to acknowledge: {path} could not be read ({e})")
    if "normalization_check" in prov:
        st = (prov["normalization_check"] or {}).get("status") \
            if isinstance(prov["normalization_check"], dict) else None
        raise ValueError(f"this DE ran under the normalisation check (status {st!r}): ack-legacy "
                         "is only for a DE from before skill 2.10, which has no check record. "
                         "Run step 8's check and the final DE (SKILL.md step 8).")
    c = staff.require(by)
    staff.append(staff_record_path(session), {
        "event": "ack-legacy", "by": c["by"], "staff_check": c["checked"],
        "account": staff.login(), "at": staff.now(), "note": note,
        "de_provenance": path, "de_provenance_sha256": sha256(path)})
    return delivery_gate(session)


# ------------------------------------------------------------------------ the verdict --
def trips(rec):
    """The checks that trip, given the default (see the module docstring)."""
    c, out = rec["checks"], []
    default = rec["default"]["quantities"]
    vol = rec["volcano"].get(default)
    if default == et.NORMALISED:
        if c["factor_confound"].get("large_and_confounded"):
            fc = c["factor_confound"]
            out.append(f"B1 normalisation factors {fc['gap_fold']:g}-fold apart between groups "
                       f"({fc['highest']} vs {fc['lowest']}), separating completely")
        for k, v in (c["raw_shift"].get("per_contrast") or {}).items():
            if v["lopsided"]:
                out.append(f"B2 {k}: before normalisation the median protein shifted "
                           f"{v['median_log2_shift']:+.2f} log2, {max(v['fraction_up'], 1 - v['fraction_up']):.0%} one way")
        if c["id_gap"].get("gap"):
            out.append(f"B3 identifications differ {c['id_gap']['ratio']:g}-fold between groups")
        if c["instability"].get("flagged"):
            out.append(f"B4 Normalisation.Instability above {c['instability']['threshold']:g} in "
                       f"{len(c['instability']['flagged'])} run(s)")
    else:
        fc = c["factor_confound"]
        if fc.get("large_unexplained"):
            out.append(f"B1 normalisation factors span {fc['run_spread_log2']:.2f} log2 across runs "
                       "without following the groups: loading differences the non-normalised "
                       "quantities would keep")
    for cn, v in ((vol or {}).get("contrasts") or {}).items():
        for a in v.get("anomalies") or []:
            out.append(f"volcano ({default}) {cn}: {a}")
    return out


def run(report, conditions, de_norm, de_raw, session=None, bait=None, out=None, report_raw=None):
    """The check. `report_raw`: the report the non-normalised DE read when it is not `report`
    (MaxLFQ: a search with --no-norm); B1-B4 read `report`, and a raw decision holds the final DE
    to `report_raw`."""
    groups = read_conditions(conditions)
    rec_type = et.load(session) if session else None
    default = et.default_for((rec_type or {}).get("type"))
    norm, why_n = read_de(de_norm) if de_norm else (None, "no --de-normalised given")
    raw, why_r = read_de(de_raw) if de_raw else (None, "no --de-raw given")
    contrasts = (norm or raw or {}).get("prov", {}).get("contrasts") or []
    design = et.pulldown_design(rec_type, sorted(set(groups.values())), contrasts)
    bait, bait_source, bait_cands = et.bait_for_check(session, bait, rec_type)
    data, why_data = read_report(report, set(groups))
    checks = {"factor_confound": factor_confound(data, groups),
              "raw_shift": raw_shift(data, groups, contrasts),
              "id_gap": id_gap(data, groups),
              "instability": instability(report, groups)}
    for k in ("factor_confound", "raw_shift", "id_gap"):
        if data is None:
            checks[k]["why"] = why_data
        checks[k].setdefault("status", "assessed")
    checks["instability"].setdefault("status", "assessed")
    rec = {
        "schema": SCHEMA, "schema_version": SCHEMA_VERSION,
        "checked_at": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "report": os.path.abspath(report), "conditions": os.path.abspath(conditions),
        # what the decision is about: gate() refuses a DE of any other report or design
        "report_sha256": sha256(report), "conditions_sha256": sha256(conditions),
        # the non-normalised DE's own report, when it is another one (MaxLFQ: --no-norm)
        "report_raw": os.path.abspath(report_raw) if report_raw else None,
        "report_raw_sha256": sha256(report_raw) if report_raw else None,
        # where decide records who decided (STAFF_RECORD in its logs/)
        "session": os.path.abspath(session) if session else None,
        "experiment_type": {"type": (rec_type or {}).get("type"),
                            "label": (rec_type or {}).get("label") or default["label"],
                            "source": (rec_type or {}).get("source") or "NOT RECORDED -- ask",
                            # the submitter's loading answer and the question it raises, if any
                            "submission_normalisation": (rec_type or {}).get("submission_normalisation"),
                            "reconcile": (rec_type or {}).get("reconcile")},
        "default": {"quantities": default["quantities"], "why": default["why"]},
        "design": design, "bait": bait, "bait_source": bait_source,
        # the sequences the user added to the database (--add-fasta), the bait's candidates
        "bait_candidates": [c["accession"] for c in bait_cands],
        "groups": {g: sum(1 for x in groups.values() if x == g) for g in sorted(set(groups.values()))},
        "checks": checks,
        "volcano": {"normalised": volcano_version(norm, design, bait,
                                                  etype=(rec_type or {}).get("type")),
                    "raw": volcano_version(raw, design, bait, raw=True,
                                           etype=(rec_type or {}).get("type"))},
        "volcano_not_assessed": {k: w for k, w in (("normalised", why_n), ("raw", why_r)) if w},
        "carboxylases": ({"normalised": carboxylase_stability(de_norm) if norm else None,
                          "raw": carboxylase_stability(de_raw) if raw else None,
                          "genes": list(et.CARBOXYLASES),
                          "note": "reference only: endogenously biotinylated carboxylases should "
                                  "be roughly constant across samples; never used to normalise"}
                         if (rec_type or {}).get("type") == "proximity" else None),
        "thresholds": {"factor_gap_log2": FACTOR_GAP_LOG2, "raw_shift_log2": RAW_SHIFT_LOG2,
                       "raw_shift_fraction": RAW_SHIFT_FRACTION, "id_gap_ratio": ID_GAP_RATIO,
                       "centre_log2": CENTRE_LOG2, "ratio": RATIO,
                       "min_sig_for_ratio": MIN_SIG_FOR_RATIO, "frac_sig": FRAC_SIG,
                       "bait_top": BAIT_TOP, "p_spike": P_SPIKE},
    }
    rec["trips"] = trips(rec)
    rec["tripped"] = bool(rec["trips"])
    rec["data_check_agrees"] = not rec["tripped"]
    # No experiment type recorded: the default above is the skill's, not the user's, so nothing
    # is decided for them -- the check asks, whatever the data say.
    rec["needs_decision"] = rec["tripped"] or rec_type is None
    if rec_type is None:
        rec["needs_decision_why"] = ("no experiment type is recorded for this session, so there is "
                                     "no default to apply: ask the user what kind of experiment "
                                     "this is (experiment_type.py set) and which quantities to use")
    if design.get("note"):
        rec["design_note"] = design["note"]
    rec["decision"] = None if rec["needs_decision"] else {
        "quantities": default["quantities"], "source": "default",
        "by": "the skill: the default for this experiment type, and the data check agreed",
        "reason": default["why"], "at": rec["checked_at"]}
    rec["recommendation"] = recommend(rec)
    out = out or os.path.dirname(os.path.abspath(de_norm or de_raw or report))
    os.makedirs(out, exist_ok=True)
    write(rec, out)
    return rec, out


def recommend(rec):
    """One recommendation from the record, said as a recommendation and never applied: the
    experiment type's DEFAULT, unless that default's own volcano is clearly the worse one
    (CLEARLY_WORSE more anomalies than the alternative's). The data checks that tripped are named
    as the questions to ask, not turned into a verdict against the type -- a whole proteome whose
    loading happens to follow the groups is not told "raw", an IP whose factors wander is not told
    to normalise. With no type recorded, the volcanoes decide. For proximity labelling the
    biotinylated carboxylases' stability under each option is added, as a reference."""
    pick, why = _recommend(rec)
    return f"{LABEL[pick]} quantities: {why}{_carboxylase_note(rec, pick)}"


LABEL = {et.RAW: "non-normalised (raw)", et.NORMALISED: "normalised"}


def _anomaly_counts(rec):
    """Counted anomalies per version (volcano_version already left out the expected ones)."""
    etype = (rec.get("experiment_type") or {}).get("type")
    out = {}
    for k in (et.NORMALISED, et.RAW):
        if k == et.RAW and etype in RAW_VOLCANO_EXPECTED:
            out[k] = 0
            continue
        out[k] = sum(len(v.get("anomalies") or []) for v in
                     ((rec.get("volcano") or {}).get(k) or {}).get("contrasts", {}).values())
    return out


def _recommend(rec):
    default = (rec.get("default") or {}).get("quantities") or et.NORMALISED
    other = et.RAW if default == et.NORMALISED else et.NORMALISED
    etype = (rec.get("experiment_type") or {}).get("type")
    if not rec.get("tripped") and etype is not None:
        return default, "the default for this experiment type; the data check agreed"
    n = _anomaly_counts(rec)
    fc = (rec.get("checks") or {}).get("factor_confound") or {}
    if etype is None:
        pick = other if n[default] > n[other] else default
        return pick, (f"no experiment type is recorded, so this follows the volcanoes "
                      f"({n[pick]} anomalies vs {n[other if pick == default else default]}); "
                      "ask what kind of experiment this is before choosing")
    if n[default] >= n[other] + CLEARLY_WORSE:
        return other, (f"the default for this type is {default}, but its own volcano is clearly "
                       f"the worse one ({n[default]} anomalies vs {n[other]})")
    why = (f"the default for this experiment type ({n[default]} volcano anomalies vs "
           f"{n[other]} for {other})")
    if default == et.NORMALISED and fc.get("large_and_confounded"):
        why += (f"; but the normalisation factors are {fc.get('gap_fold')}-fold apart and follow "
                "the groups -- ask whether one group carries less material by design (an "
                "enrichment): if so, the experiment type is wrong and raw quantities fit")
    elif default == et.RAW and fc.get("large_unexplained"):
        why += ("; but the normalisation factors vary across runs without following the groups -- "
                "ask whether the samples were loaded unevenly, which non-normalised quantities keep")
    return default, why


def _carboxylase_note(rec, pick):
    cx = rec.get("carboxylases") or {}
    n, r = ((cx.get(k) or {}).get("median_sd_log2") for k in (et.NORMALISED, et.RAW))
    if n is None or r is None:
        return ""
    steadier = et.NORMALISED if n < r else et.RAW
    note = (f". Reference only: the biotinylated carboxylases are steadier with {steadier} "
            f"quantities (median SD {min(n, r):.2f} vs {max(n, r):.2f} log2 across samples)")
    if steadier != pick:
        return note + f" -- which argues for {steadier} quantities: weigh it before choosing"
    if steadier == et.NORMALISED and rec.get("tripped"):
        return note + (" -- consistent with the controls carrying less material, which "
                       "normalisation corrects")
    return note


def write(rec, out):
    path = os.path.join(out, RECORD)
    with open(path + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=2)
        fh.write("\n")
    os.replace(path + ".tmp", path)
    with open(os.path.join(out, SUMMARY), "w", encoding="utf-8") as fh:
        fh.write(summary_md(rec))
    return path


def summary_md(rec):
    """The side-by-side summary the user chooses from (tolerant of a record missing a section:
    what is not there is said to be not recorded)."""
    rec = dict(rec)
    rec.setdefault("experiment_type", {"label": "NOT RECORDED", "source": "NOT RECORDED"})
    rec["default"] = dict({"why": "not recorded"}, **(rec.get("default") or {}))
    for k, v in (("trips", []), ("tripped", bool(rec.get("trips"))), ("volcano", {}),
                 ("volcano_not_assessed", {}), ("recommendation", "not recorded"),
                 ("decision", None)):
        rec.setdefault(k, v)
    rec["checks"] = {k: (rec.get("checks") or {}).get(k) or {"why": "not recorded"}
                     for k in ("factor_confound", "raw_shift", "id_gap", "instability")}
    L = ["# Normalisation check", "",
         f"- **Experiment type:** {rec['experiment_type']['label']} ({rec['experiment_type']['source']})",
         f"- **Default for it:** {rec['default']['quantities']} quantities -- {rec['default']['why']}",
         f"- **Bait:** " + (f"{rec['bait']} ({rec.get('bait_source') or 'source not recorded'})"
                            if rec.get("bait") else
                            f"none ({rec.get('bait_source') or 'not recorded'})"),
         f"- **Data check:** " + ("agreed with the default" if not rec["tripped"] else
                                  "**TRIPPED -- stop and ask before the final DE**"), ""]
    if rec["experiment_type"].get("reconcile"):
        L[-1:-1] = [f"- **Ask the user (the submission):** {rec['experiment_type']['reconcile']}"]
    if rec.get("needs_decision_why"):
        L[-1:-1] = [f"- **Ask the user:** {rec['needs_decision_why']}"]
    for t in rec["trips"]:
        L.append(f"  - {t}")
    if rec.get("design_note"):
        L += ["", f"- **Design:** {rec['design_note']}"]
    c = rec["checks"]
    fc = c["factor_confound"]
    L += ["", "## Checks", "",
          "| check | value | threshold |", "|---|---|---|",
          f"| B1 normalisation factor by group | "
          + (", ".join(f"{g} {v['median_log2_factor']:+.2f}" for g, v in fc["per_group"].items())
             + f" log2; gap {fc['gap_fold']:g}-fold, separated {fc['separated']}, "
               f"permutation p {fc['permutation_p']}" if fc.get("per_group")
             else f"not assessed ({fc.get('why')})")
          + f" | {fc.get('threshold', '')} |"]
    rs = c["raw_shift"]
    for k, v in (rs.get("per_contrast") or {}).items():
        L.append(f"| B2 shift before normalisation, {k} | median {v['median_log2_shift']:+.2f} "
                 f"log2, {v['fraction_up']:.0%} up ({v['n_proteins']} proteins) | "
                 f"{rs['threshold']} |")
    if not rs.get("per_contrast"):
        L.append(f"| B2 shift before normalisation | not assessed ({rs.get('why')}) | |")
    ig = c["id_gap"]
    L.append(f"| B3 precursors per group (median) | "
             + (", ".join(f"{g} {v['median_precursors']:g}" for g, v in ig["per_group"].items())
                + f"; ratio {ig['ratio']:g}" if ig.get("per_group") else
                f"not assessed ({ig.get('why')})") + f" | {ig.get('threshold', '')} |")
    ni = c["instability"]
    L.append(f"| B4 Normalisation.Instability | "
             + (f"median {ni['median']}, max {ni['max']}, {len(ni['flagged'])} above "
                f"{ni['threshold']:g}" if "median" in ni else f"not assessed ({ni.get('why')})")
             + " | |")
    L += ["", "## Side by side: normalised vs non-normalised", "",
          "| contrast | version | significant (up / down) | % significant | centre (non-sig) "
          "| up:down | bait rank (log2FC) | p-values | anomalies |",
          "|---|---|---|---|---|---|---|---|---|"]
    names = sorted(set((((rec["volcano"].get("normalised") or {}).get("contrasts") or {}).keys())
                       | set(((rec["volcano"].get("raw") or {}).get("contrasts") or {}).keys())))
    for cn in names:
        for k in (et.NORMALISED, et.RAW):
            v = ((rec["volcano"].get(k) or {}).get("contrasts") or {}).get(cn)
            if not v or "n" not in v:
                L.append(f"| {cn} | {k} | not assessed ({rec['volcano_not_assessed'].get(k, '')}) "
                         "| | | | | | |")
                continue
            b = v.get("bait") or {}
            L.append(f"| {cn} | {k} | {v['n_sig']} ({v['up']} / {v['down']}) | "
                     f"{_fmt(v['frac_sig'], '.0%')} | {_fmt(v['centre_nonsig'], '+.2f')} | "
                     f"{v['up_down_ratio'] if v['up_down_ratio'] is not None else 'n/a'} | "
                     + (f"{b.get('rank')} ({_fmt(b.get('log2fc'), '+.2f')})" if b.get("rank") else
                        ("not found" if b else "--"))
                     + f" | {v['pvalue']['shape']} | {'; '.join(v['anomalies']) or 'none'} |")
    cx = rec.get("carboxylases")
    if cx:
        L += ["", "## Biotinylated carboxylases (reference only, never a normaliser)", ""]
        for k in (et.NORMALISED, et.RAW):
            v = cx.get(k)
            L.append(f"- {k}: " + (("median SD " + str(v["median_sd_log2"]) + " log2 across samples ("
                                    + ", ".join(f"{g} {s}" for g, s in v["sd_log2"].items()) + ")")
                                   if v and v.get("n_found") else "none found"))
    L += ["", f"**Recommendation:** {rec['recommendation']}", ""]
    if rec["decision"]:
        d = rec["decision"]
        said = (d.get("statement")
                or f"{d['quantities']} quantities -- {d.get('by')}: {d.get('reason')}")
        L.append(f"**Decision:** {said}")
    else:
        L.append("**Decision: not made.** Show this to the Core staff member, ask which quantities "
                 "to use -- the decision is theirs, never Claude's -- and record their answer: "
                 "`normalization_check.py decide --check <this json> --quantities raw|normalised "
                 "--by <their HIVE login> --reason \"<why, in their words>\"`.")
    return "\n".join(L) + "\n"


def load(path):
    with open(path, encoding="utf-8") as fh:
        rec = json.load(fh)
    if not isinstance(rec, dict) or rec.get("schema") != SCHEMA:
        raise ValueError(f"{path} is not a {SCHEMA} record")
    return rec


def statement(quantities, reason, who=staff.ROLE):
    """THE sentence client documents give a person's decision: by role, never by name."""
    r = re.sub(r"^because\s+", "", (reason or "").strip(), flags=re.I).rstrip(". ")
    return f"{who} chose {LABEL.get(quantities, quantities)} quantities because {r}."


def decision_staff_record(path, rec):
    """Where `decide` records who decided -> (the file, as the check record names it): the
    session's logs/ (the session the check ran with, else the one its folder is in), else beside
    the check."""
    here = os.path.dirname(os.path.abspath(path))
    s = rec.get("session") or (os.path.dirname(os.path.dirname(here))
                               if here.endswith(os.path.join("output", "norm_check")) else None)
    if s:
        return staff_record_path(s), os.path.join("logs", STAFF_RECORD)
    return os.path.join(here, STAFF_RECORD), STAFF_RECORD


def decide(path, quantities, by, reason):
    """Record a Core staff member's decision -- theirs, never Claude's: --by is their HIVE login
    (staff.require). The check record (and so de_provenance.json, the Methods, the zip) names the
    role; their name and login go to the staff-only record."""
    if quantities not in (et.RAW, et.NORMALISED):
        raise ValueError("--quantities must be raw or normalised")
    if not (reason and reason.strip()):
        raise ValueError("--reason is required: why they chose it, in their words")
    rec = load(path)
    c = staff.require(by)
    at = staff.now()
    where, named = decision_staff_record(path, rec)
    staff.append(where, {"event": "decide", "by": c["by"], "staff_check": c["checked"],
                         "account": staff.login(), "at": at, "quantities": quantities,
                         "reason": reason.strip(), "check": os.path.abspath(path)})
    rec["decision"] = {"quantities": quantities, "source": "user", "by": staff.ROLE,
                       "reason": reason.strip(), "statement": statement(quantities, reason),
                       "at": at, "default_was": rec["default"]["quantities"],
                       "data_check_tripped": rec["tripped"],
                       # the person, by name and login: staff-only, never delivered
                       "staff_record": named}
    write(rec, os.path.dirname(os.path.abspath(path)))
    return rec


def gate(path, quantities, report=None, conditions=None):
    """(ok, message, record) for run_de.R: the DE may read only the decided quantities, of the
    report and design the decision was made on (their sha256, recorded by run())."""
    rec = load(path)
    # MaxLFQ's non-normalised DE reads its own (--no-norm) report: a raw DE is held to that one
    rkey = "report_raw" if quantities == et.RAW and rec.get("report_raw_sha256") else "report"
    for what, given, key in (("report", report, rkey), ("conditions", conditions, "conditions")):
        if given is None:
            continue
        want, got = rec.get(key + "_sha256"), sha256(given)
        if not want or want != got:
            return False, (f"REFUSED: the normalisation check {path} was made on "
                           f"{rec.get(key) or 'a ' + what + ' it did not record'}"
                           + (" (its non-normalised DE's report, --report-raw)"
                              if key == "report_raw" else "")
                           + ("" if want else " (a check from before 2.10 records no sha256)")
                           + f", not on this DE's {what} {given} (sha256 differs). Run the check "
                           "again on this DE's report and conditions -- a decision about one "
                           "report never authorises another."), rec
    d = rec.get("decision")
    if not d:
        why = "; ".join(rec.get("trips") or []) or rec.get("needs_decision_why") or "undecided"
        return False, (f"REFUSED: the normalisation check needs a decision ({why}) and none is "
                       f"recorded. Show {SUMMARY} to the user, ask, and record the answer with "
                       "normalization_check.py decide."), rec
    if d["quantities"] != quantities:
        return False, (f"REFUSED: the recorded decision is {d['quantities']} quantities "
                       f"({d['by']}: {d['reason']}), and this DE asked for {quantities}. Run it "
                       f"with --quantities {d['quantities']}, or record a new decision -- never a "
                       "silent switch."), rec
    return True, (f"{quantities} quantities, as decided ({d['by']}); data check "
                  + ("agreed" if not rec["tripped"] else "tripped: " + "; ".join(rec["trips"]))), rec


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--report", required=True)
    r.add_argument("--conditions", required=True)
    r.add_argument("--de-normalised")
    r.add_argument("--de-raw")
    r.add_argument("--session")
    r.add_argument("--bait")
    r.add_argument("--out")
    r.add_argument("--report-raw", help="the report the non-normalised DE read, when it is not "
                                        "--report (MaxLFQ: a search with --no-norm)")
    d = sub.add_parser("decide", help="record a Core staff member's decision (theirs, never "
                                      "Claude's)")
    d.add_argument("--check", required=True)
    d.add_argument("--quantities", required=True, choices=(et.RAW, et.NORMALISED))
    d.add_argument("--by", required=True, help="their HIVE login, listed in the Core's staff file "
                                               "(staff.py)")
    d.add_argument("--reason", required=True, help="why, in their words -- shown to the client, "
                                                   "so no names")
    k = sub.add_parser("ack-legacy", help="a DE from before 2.10 (no check record) stands: let "
                                          "its delivery through")
    k.add_argument("--session", required=True)
    k.add_argument("--by", required=True, help="the Core staff member's HIVE login (staff.py)")
    k.add_argument("--note", required=True, help="why this DE stands, one line")
    g = sub.add_parser("gate")
    g.add_argument("--check", required=True)
    g.add_argument("--quantities", required=True)
    g.add_argument("--report", help="the DE's --input: must be the report the check was made on")
    g.add_argument("--conditions", help="the DE's --metadata: must be the check's conditions")
    nn = sub.add_parser("no-norm", help="did the search behind a report run with --no-norm, by "
                                        "its provenance? (run_de.R's MaxLFQ)")
    nn.add_argument("--report", required=True)
    nn.add_argument("--search-dir")
    rq = sub.add_parser("required", help="must a DE of these files run under a check? (run_de.R)")
    rq.add_argument("--input", required=True)
    rq.add_argument("--metadata", required=True)
    rq.add_argument("--method", default="dpc")
    a = ap.parse_args()
    try:
        if a.cmd == "run":
            if not (a.de_normalised or a.de_raw):
                raise ValueError("give the DE folders to compare: --de-normalised and --de-raw")
            rec, out = run(a.report, a.conditions, a.de_normalised, a.de_raw, a.session, a.bait,
                           a.out, a.report_raw)
            print(json.dumps({"check": os.path.join(out, RECORD),
                              "summary": os.path.join(out, SUMMARY),
                              "tripped": rec["tripped"], "trips": rec["trips"],
                              "needs_decision": rec["needs_decision"],
                              "default": rec["default"], "decision": rec["decision"],
                              "recommendation": rec["recommendation"]}, indent=2))
            if rec["needs_decision"]:
                print(f"[normalization_check] {'TRIPPED' if rec['tripped'] else 'NO EXPERIMENT TYPE'}"
                      f" -- stop before the final DE: show {os.path.join(out, SUMMARY)} to the "
                      "user and ask.", file=sys.stderr)
                return TRIPPED_EXIT
        elif a.cmd == "decide":
            rec = decide(a.check, a.quantities, a.by, a.reason)
            print(json.dumps({"decision": rec["decision"]}, indent=2))
        elif a.cmd == "ack-legacy":
            g = ack_legacy(a.session, a.by, a.note)
            print(json.dumps({"gate": g}, indent=2))
            if not g["proceed"]:
                return 2
        elif a.cmd == "no-norm":
            print(json.dumps(searched_no_norm(a.report, a.search_dir)))
        elif a.cmd == "required":
            print(json.dumps(required(a.input, a.metadata, a.method)))
        else:
            ok, msg, rec = gate(a.check, a.quantities, a.report, a.conditions)
            print(json.dumps({"ok": ok, "message": msg, "record": rec}))
    except (ValueError, OSError) as e:
        sys.exit(f"[normalization_check] {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
