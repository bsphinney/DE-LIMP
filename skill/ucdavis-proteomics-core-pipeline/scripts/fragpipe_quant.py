#!/usr/bin/env python3
"""
fragpipe_quant.py -- the protein quantity of a FragPipe DDA search, and directLFQ (with its guard)
for every search with more than one experiment that IonQuant gave no MaxLFQ.

FragPipe's IonQuant computes MaxLFQ only from 3 experiments on: with 2 it logs "There are only 2
experiments. Will not calculate MaxLFQ intensity." and combined_protein.tsv carries only each
experiment's protein `Intensity` (FragPipe issues #1389, #1770). That Intensity sums the ions IonQuant
saw in each run, including ions seen in one run only, so a deeper run looks more abundant: on a
two-run DDA case it put the second run +1.67 log2 above the first, where directLFQ (Ammar et al.
2023), on the ions measured in both runs, gave +0.14 (Ruling 64; 38 protein groups). So:

  * IonQuant's own `<exp> MaxLFQ Intensity` when it wrote one for every experiment (3 or more
    experiments, FragPipe's default) -- unchanged;
  * otherwise, with 2 or more experiments (2 experiments, or MaxLFQ turned off or missing for some):
    directLFQ's protein intensities from combined_ion.tsv are the report's quantity, and IonQuant's
    `<exp> Intensity` goes beside them (COMPANION_COLUMN) for every protein group -- the Intensity
    is inflated by depth whatever the run count, and directLFQ is built for many runs (Ruling 102);
  * 1 experiment: IonQuant's `<exp> Intensity`; there is no between-run step, and the declaration
    says so.

The quantity columns are EXACTLY "<experiment> MaxLFQ Intensity" for every experiment, else exactly
"<experiment> Intensity" -- never a mix, and never a Unique / Total / Razor variant. An experiment of
the run map with neither column stops the adapter, naming it. plan() returns the one self-describing
quantity object (architectural rule 1): what the quantity is, how many experiments, why it is not
IonQuant's MaxLFQ; run_search.fragpipe_quantity_declaration() is the one path that turns it into the
report's declaration.

directLFQ runs in its OWN pinned venv (DIRECTLFQ_PINS; built once, in a job, by
build_directlfq_venv.sh), never in the pipeline env: directLFQ 0.3.3 under pandas 3 silently skips
its between-sample normalisation (copy-on-write: it normalises a copy and hands back the original;
H6). The guard (Ruling 90, defined once: guard()) refuses the result unless ALL of these hold:
  1. pandas < 3;
  2. directLFQ's COPY_NUMPY_ARRAYS_DERIVED_FROM_PANDAS is False after its own check (the code path
     that discards the normalisation is the True one);
  3. the shift APPLIED equals the shift directLFQ FOUND, per sample: the normalised protein table
     against the one from deactivate_normalization=True, columns matched by name, relative to a
     reference sample, compared with the shifts its normaliser computed in its own frame. A real
     shift can be 0.021 log2 (C0), so there is no size threshold: only GUARD_TOL for float noise.

This module imports only the standard library at the top: the adapter imports it in the pipeline
env, and the venv runs it as a script (`directlfq`, `check-env`).

  python3 fragpipe_quant.py directlfq --out <FragPipe output dir> [--cores N]   (inside the venv)
  python3 fragpipe_quant.py check-env                  guard checks 1-2 of this python, as JSON
  python3 fragpipe_quant.py pins                       the venv's requirements, one per line
  python3 fragpipe_quant.py venv [--default-prefix]    where the venv is, or why there is none

Exit: 0 ok; 3 refused by the guard (the record says why); 2 usage / missing input.
"""
import argparse
import datetime
import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from safe_write import sha256_file  # noqa: E402  the ONE chunked sha256

# IonQuant writes MaxLFQ columns only from this many experiments (fcyu, FragPipe #1389 and #1770;
# H5: "There are only 2 experiments. Will not calculate MaxLFQ intensity.").
MIN_MAXLFQ_EXPERIMENTS = 3
# directLFQ's own venv (H6/C0 recipe; Rulings 86, 90): 0.3.3 imports yaml without declaring it;
# pandas 2.2.3 is the version directLFQ's own "stable" extra pins, and pandas 3 skips its
# normalisation.
DIRECTLFQ_VERSION = "0.3.3"
DIRECTLFQ_PINS = (f"directlfq=={DIRECTLFQ_VERSION}", "pandas==2.2.3", "pyyaml")
VENV_ENV = "DIRECTLFQ_VENV"
# THE Core's venv on HIVE (Ruling 102: built by build_directlfq_venv.sh in a job, versioned by its
# pins, group-readable, not world-writable; the 20 GB home is no place for a venv: ucd-hpc). Off
# HIVE, this account's own tools folder (a laptop, another cluster).
GROUP_VENV = "/quobyte/proteomics-grp/tools/directlfq/0.3.3-pandas2.2.3"
HOME_VENV = os.path.join("~", ".proteomics-pipeline", "venvs", "directlfq-0.3.3-pandas2.2.3")
VENV_MARKER = "directlfq_venv.json"          # written last by build_directlfq_venv.sh
BUILDER = "build_directlfq_venv.sh"

ION_TABLE = "combined_ion.tsv"               # IonQuant's ion table: directLFQ's input
DIRECTLFQ_DIR = "directlfq"                  # in the search folder: directLFQ's tables + record
DIRECTLFQ_RECORD = "directlfq_record.json"
RECORD = "fragpipe_quantity.json"            # the adapter's record (+ search_provenance.json key)
PROVENANCE_KEY = "fragpipe_quantity"
COMPANION_COLUMN = "PG.IonQuant.Intensity"   # report.parquet, beside PG.MaxLFQ = directLFQ
ANNOTATION_COLUMNS = ("Genes", "Protein.Names")   # report.parquet, the leading protein's
GUARD_TOL = 1e-4                             # log2; float noise only, never a size threshold

QUANTITY_DIRECTLFQ = "directLFQ"
QUANTITY_MAXLFQ = "IonQuant MaxLFQ Intensity"
QUANTITY_INTENSITY = "IonQuant Intensity"
# Variants IonQuant/Philosopher may write per experiment that are never the protein quantity.
QUALIFIERS = ("Unique", "Total", "Razor")
IONQUANT_NO_MAXLFQ = re.compile(
    r"There (?:are|is) only \d+ experiments?\.\s*Will not calculate MaxLFQ intensit(?:y|ies)\.?")
NO_MAP_CAVEAT = ("No run map (the search's fragpipe.fp-manifest and FragPipe's own "
                 "fragpipe-files.fp-manifest): the Run names are FragPipe's experiment names, "
                 "which can differ from the file names (FragPipe turns '-' into '_'), and the runs "
                 "were not checked against the search's inputs.")
COMPANION_TEXT = ("IonQuant's protein Intensity (combined_protein.tsv '<experiment> Intensity') "
                  "for every protein group, beside directLFQ: it sums ions seen in one run only, "
                  "so it is not comparable between runs whose depth differs.")


class QuantError(Exception):
    """A FragPipe quantity that cannot be delivered as declared: the adapter stops with this."""


# ------------------------------------------------------------- which column -----
def quantity_columns(header, experiments):
    """{"kind", "columns" {experiment: column}, "maxlfq_missing" [experiments]}: exactly
    "<exp> MaxLFQ Intensity" for EVERY experiment when all are there, else exactly
    "<exp> Intensity" for every experiment -- never a mix, never another column. Raises QuantError
    naming each experiment that has neither (a run missing from the table)."""
    if not experiments:
        raise QuantError("no experiments to quantify (the run map and combined_protein.tsv's "
                         "columns name none)")
    cols = set(header)
    maxlfq = {e: f"{e} MaxLFQ Intensity" for e in experiments}
    plain = {e: f"{e} Intensity" for e in experiments}
    if all(c in cols for c in maxlfq.values()):
        return {"kind": "MaxLFQ Intensity", "columns": maxlfq, "maxlfq_missing": []}
    missing = [e for e in experiments if plain[e] not in cols]
    if missing:
        raise QuantError(
            f"combined_protein.tsv has no quantity column for {len(missing)} of the "
            f"{len(experiments)} runs: {', '.join(missing)} (neither '<experiment> MaxLFQ "
            f"Intensity' nor '<experiment> Intensity'). FragPipe left these runs out of its "
            f"quantification; nothing was built. Check FragPipe's log for them.")
    return {"kind": "Intensity", "columns": plain,
            "maxlfq_missing": [e for e in experiments if maxlfq[e] not in cols]}


def experiments_from_header(header):
    """The experiments of a combined_protein.tsv WITHOUT a run map (an output from before 2.14.2,
    or made outside the skill), in table order: each "<exp> Intensity" / "<exp> MaxLFQ Intensity"
    column's <exp>, never a Unique / Total / Razor variant read as an experiment."""
    found = []
    for c in header:
        if not c.endswith(" Intensity"):
            continue
        p = c[:-len(" Intensity")]
        if p.endswith(" MaxLFQ"):
            p = p[:-len(" MaxLFQ")]
        words = p.rsplit(" ", 1)
        if len(words) == 2 and words[1] in QUALIFIERS:
            continue
        if p and p != "Combined" and p not in found:
            found.append(p)
    return found


def quantity_column_names(header):
    """Every column of `header` that can be a protein quantity: exactly "<exp> Intensity" and
    "<exp> MaxLFQ Intensity" of experiments_from_header() -- for a yes/no "did it quantify"
    (run_search._combined_protein_intensities), so a Unique / Total variant never decides it."""
    cols = set(header)
    return [c for e in experiments_from_header(header)
            for c in (f"{e} Intensity", f"{e} MaxLFQ Intensity") if c in cols]


def ionquant_maxlfq_line(folder):
    """IonQuant's own words for not computing MaxLFQ, from FragPipe's newest log in `folder`
    (log_*.txt, then fragpipe_stdout.txt), or None."""
    import fragpipe_params                      # the one list of FragPipe's logs, newest first
    logs = fragpipe_params.logs(folder) + [os.path.join(folder or "", "fragpipe_stdout.txt")]
    for log in logs:
        try:
            with open(log, errors="replace") as fh:
                for line in fh:
                    m = IONQUANT_NO_MAXLFQ.search(line)
                    if m:
                        return m.group(0).strip()
        except OSError:
            continue
    return None


def _workflow_maxlfq(folder):
    """ionquant.maxlfq of the workflow FragPipe saved in `folder` ("1"/"0"), or None."""
    try:
        import fragpipe_params
        wf = fragpipe_params.find(folder, fragpipe_params.SAVED_WORKFLOW)
        return fragpipe_params.read_workflow(wf)[0].get("ionquant.maxlfq") if wf else None
    except Exception:                       # the reason is then said without the workflow
        return None


# One experiment (FragPipe's grouping of runs: the skill gives each run its own, so one means one
# quantified column) has no between-run step. Said in plain words wherever that quantity is
# described -- the Methods (make_methods.fragpipe_quantity_sentence) and the report
# (make_analysis_html.search_quantity_clause) -- never the bare word "experiment" (C2 bio M1).
ONE_EXPERIMENT_WHY = ("one experiment (FragPipe's grouping of runs) was quantified, so there is "
                      "no between-run step (no MaxLFQ, no directLFQ): the quantity is IonQuant's "
                      "protein Intensity as it reported it")
NOT_NORMALISED = "not normalised between runs; compare runs with care"


# The one sentence for "directLFQ because it was asked for": make_presets.py --directlfq (or an
# --overrides ionquant.maxlfq=0) writes ionquant.maxlfq=0 into the workflow, no vendor workflow
# does. Said once here; the Methods, the report and the record read it from the plan.
REQUESTED_WORDS = ("directLFQ was requested instead of IonQuant's MaxLFQ (the workflow sets "
                   "ionquant.maxlfq=0).")


def reason_phrase(q):
    """A few words for a plan/declaration `q`'s reason for directLFQ, for a one-line clause (the
    report): "requested instead of IonQuant's MaxLFQ", "2 experiments: IonQuant computes MaxLFQ
    only from 3 on", or None when there is no reason to give."""
    why = (q or {}).get("why_not_maxlfq") or ""
    n = (q or {}).get("n_experiments")
    if why.startswith(REQUESTED_WORDS.split(" (")[0]):
        return "requested instead of IonQuant's MaxLFQ"
    if isinstance(n, int) and n < MIN_MAXLFQ_EXPERIMENTS and why:
        return f"{n} experiments: IonQuant computes MaxLFQ only from {MIN_MAXLFQ_EXPERIMENTS} on"
    return None


def why_not_maxlfq(n, ionquant, folder):
    """Why the report's quantity is not IonQuant's MaxLFQ (None when it is)."""
    if ionquant["kind"] == "MaxLFQ Intensity" and n >= MIN_MAXLFQ_EXPERIMENTS:
        return None
    line = ionquant_maxlfq_line(folder)
    quoted = f' IonQuant\'s log: "{line}"' if line else ""
    requested = _workflow_maxlfq(folder) == "0"      # no vendor workflow ships 0: someone chose it
    if n < MIN_MAXLFQ_EXPERIMENTS:
        return ((REQUESTED_WORDS + " " if requested else "")
                + f"IonQuant computes MaxLFQ only from {MIN_MAXLFQ_EXPERIMENTS} experiments on "
                f"(FragPipe issues #1389, #1770); this search has {n}.{quoted}")
    if ionquant["maxlfq_missing"] and len(ionquant["maxlfq_missing"]) < n:
        return (f"IonQuant wrote MaxLFQ columns for only some experiments (none for "
                f"{', '.join(ionquant['maxlfq_missing'])}), so the plain Intensity is used for "
                f"all of them, never a mix.{quoted}")
    if requested:
        return f"{REQUESTED_WORDS}{quoted}"
    return ("IonQuant wrote no MaxLFQ columns" + (f".{quoted}" if line else
            "; its log does not say why [confirm]."))


def plan(header, run_map=None, folder=None):
    """THE quantity of a FragPipe DDA search (rule 1: what it is, how many experiments, why not
    MaxLFQ), from combined_protein.tsv's `header` and the run map (run_search.
    fragpipe_experiment_map: [{"experiment", "run", "input"}], None without one). `folder`: where
    FragPipe's log and workflow are. Raises QuantError (a run missing from the table, a map that
    names one experiment or one run twice)."""
    caveats = []
    if run_map:
        unsearched = [m["run"] for m in run_map if not m.get("experiment")]
        if unsearched:
            raise QuantError(f"FragPipe's own manifest (fragpipe-files.fp-manifest) does not list "
                             f"{len(unsearched)} run(s) of this search: {', '.join(unsearched)}; "
                             f"FragPipe did not search them. Nothing was built.")
        exps = [m["experiment"] for m in run_map]
        runs = [m["run"] for m in run_map]
        for what, xs in (("FragPipe experiment", exps), ("run", runs)):
            twice = sorted({x for x in xs if xs.count(x) > 1})
            if twice:
                raise QuantError(f"the run map names a {what} more than once ({', '.join(twice)}): "
                                 f"FragPipe would have pooled those runs into one experiment. "
                                 f"Nothing was built; search them with distinct names.")
        experiments = [{"experiment": m["experiment"], "run": m["run"],
                        "input": os.path.basename(m.get("input") or "") or None} for m in run_map]
    else:
        experiments = [{"experiment": e, "run": e, "input": None}
                       for e in experiments_from_header(header)]
        caveats.append(NO_MAP_CAVEAT)
    names = [e["experiment"] for e in experiments]
    ionquant = quantity_columns(header, names)          # raises: a run missing from the table
    n = len(names)
    why = None
    if n >= MIN_MAXLFQ_EXPERIMENTS and ionquant["kind"] == "MaxLFQ Intensity":
        quantity = QUANTITY_MAXLFQ
    elif n >= 2:
        # no MaxLFQ for every experiment (2 experiments, or MaxLFQ off / missing for some):
        # directLFQ is the quantity at any run count, IonQuant's plain Intensity beside it
        quantity = QUANTITY_DIRECTLFQ
        missing = [e for e in names if f"{e} Intensity" not in set(header)]
        if missing:
            raise QuantError(f"combined_protein.tsv has no '<experiment> Intensity' column for "
                             f"{', '.join(missing)}, which goes beside directLFQ; nothing was "
                             f"built.")
        ionquant = {"kind": "Intensity", "columns": {e: f"{e} Intensity" for e in names},
                    "maxlfq_missing": [e for e in names
                                       if f"{e} MaxLFQ Intensity" not in set(header)]}
        why = why_not_maxlfq(n, ionquant, folder)
    else:
        quantity = QUANTITY_INTENSITY
        why = ONE_EXPERIMENT_WHY
    return {"quantity": quantity, "n_experiments": n, "experiments": experiments,
            "ionquant": {"kind": ionquant["kind"], "columns": ionquant["columns"]},
            "why_not_maxlfq": why, "run_map": bool(run_map),
            "companion_column": COMPANION_COLUMN if quantity == QUANTITY_DIRECTLFQ else None,
            "caveat": " ".join(caveats) or None}


def dda_report_columns(companion=None):
    """THE column list of the protein x run report.parquet adapt_fragpipe_dda writes (one row per
    protein group and run; the q-value columns are 0.0 placeholders: the identification rule is
    FragPipe's Philosopher filter, whose protein level the record's `protein_fdr` states -- with
    --prot 1 the protein groups carry no FDR at all). The adapter's table and the analysis brief
    both read this."""
    from diann_q_columns import FDR_REQUIRED
    return (["Run", "Protein.Group", "PG.MaxLFQ", *ANNOTATION_COLUMNS, *FDR_REQUIRED]
            + ([companion] if companion else []))


# ------------------------------------------------------------------- the venv --
def venv_candidates():
    """[(prefix, where)] in the order looked at: $DIRECTLFQ_VENV alone when set (never a silent
    fall-back past it), else the Core's group location, then this account's tools folder."""
    env = os.environ.get(VENV_ENV)
    if env:
        return [(os.path.expanduser(env), f"${VENV_ENV}")]
    return [(GROUP_VENV, "the Core's group location"),
            (os.path.expanduser(HOME_VENV), "this account's tools folder")]


def default_prefix():
    """Where build_directlfq_venv.sh builds when given no --prefix: $DIRECTLFQ_VENV, else the
    group location when its share exists (HIVE), else this account's tools folder."""
    env = os.environ.get(VENV_ENV)
    if env:
        return os.path.expanduser(env)
    share = os.path.dirname(os.path.dirname(os.path.dirname(GROUP_VENV)))   # the group share
    if os.path.isdir(share):
        return GROUP_VENV
    return os.path.expanduser(HOME_VENV)


def read_marker(prefix):
    try:
        with open(os.path.join(prefix, VENV_MARKER)) as fh:
            m = json.load(fh)
    except (OSError, ValueError):
        return None
    return m if isinstance(m, dict) else None


def venv_python():
    """(python, where) of directLFQ's pinned venv, or (None, why): ready only with the builder's
    marker, whose pins equal DIRECTLFQ_PINS (a venv built to other pins is not this one)."""
    looked = []
    for prefix, where in venv_candidates():
        py = os.path.join(prefix, "bin", "python")
        m = read_marker(prefix)
        if m and list(m.get("pins") or []) == list(DIRECTLFQ_PINS) and os.access(py, os.X_OK):
            return py, f"{where} ({prefix})"
        looked.append(f"{prefix} ({where}): " + (
            "not there" if not os.path.isdir(prefix) else
            f"no {VENV_MARKER} (not built by {BUILDER}, or the build did not finish)" if not m
            else f"built to other pins {m.get('pins')}" if list(m.get('pins') or [])
            != list(DIRECTLFQ_PINS) else "no runnable bin/python"))
    return None, "; ".join(looked)


def builder_command():
    """The one line that builds the venv, for a refusal's fix."""
    return (f"python3 {os.path.join(HERE, 'on_compute.py')} --time 30 --mem 8G --cpus 2 "
            f"--why 'build directLFQ venv' -- bash {os.path.join(HERE, BUILDER)}")


# ------------------------------------------------------- reading the tables -----
def read_protein_table(path):
    """directLFQ's <input>.protein_intensities.tsv -> {protein: {sample: intensity > 0}} with the
    samples BY COLUMN NAME (pandas 2 and 3 write them in different orders; C0). 0 and blanks =
    not quantified, left out."""
    import csv
    out, samples = {}, []
    with open(path, newline="", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        if not hdr or hdr[0] != "protein":
            raise QuantError(f"{os.path.basename(path)}: not a directLFQ protein table (first "
                             f"column {hdr[:1]}, expected 'protein')")
        samples = hdr[1:]
        for row in rd:
            if not row:
                continue
            vals = {}
            for s, v in zip(samples, row[1:]):
                try:
                    x = float(v)
                except ValueError:
                    continue
                if x > 0 and math.isfinite(x):
                    vals[s] = x
            out[row[0]] = vals
    return out, samples


def _median(xs):
    xs = sorted(xs)
    n = len(xs)
    if not n:
        return None
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


# ---------------------------------------------------------------- the guard -----
def pandas_major(version):
    m = re.match(r"\s*(\d+)\.", str(version or ""))
    return int(m.group(1)) if m else None


def env_checks(env):
    """Guard checks 1 and 2 on an environment record (environment()): pandas < 3; directLFQ's
    COPY flag False after its own check. Also the pinned directLFQ version (the guard reads 0.3.3's
    normaliser). -> {name: {"passed", "detail"}}."""
    major = pandas_major(env.get("pandas"))
    return {
        "pandas_below_3": {
            "passed": major is not None and major < 3,
            "detail": f"pandas {env.get('pandas')}" + (
                "" if major is not None and major < 3 else
                " -- directLFQ 0.3.3 silently skips its normalisation under pandas 3")},
        "copy_flag_false": {
            "passed": env.get("copy_flag") is False,
            "detail": f"directlfq.config.COPY_NUMPY_ARRAYS_DERIVED_FROM_PANDAS = "
                      f"{env.get('copy_flag')}" + (
                          "" if env.get("copy_flag") is False else
                          " -- the path on which 0.3.3 hands back the un-normalised table "
                          "('Some numpy arrays derived from pandas will be copied.')")},
        "directlfq_pinned": {
            "passed": env.get("directlfq") == DIRECTLFQ_VERSION,
            "detail": f"directLFQ {env.get('directlfq')} (pinned {DIRECTLFQ_VERSION})"},
    }


def shift_check(found, normalised, unnormalised, tol=GUARD_TOL):
    """Guard check 3: per sample, the shift APPLIED (log2 normalised - log2 unnormalised protein
    intensity, relative to a reference sample, median over the proteins quantified in both tables
    in both samples; samples matched BY NAME) equals the shift FOUND (directLFQ's normaliser, its
    own frame, relative to the same reference). `found`: {sample: log2 shift or None};
    `normalised` / `unnormalised`: read_protein_table()[0]. No size threshold: `tol` is float
    noise only (a real shift can be 0.021 log2)."""
    samples_n = sorted({s for v in normalised.values() for s in v})
    samples_u = sorted({s for v in unnormalised.values() for s in v})
    res = {"passed": False, "tolerance_log2": tol, "found": {}, "applied": {}, "reference": None,
           "max_abs_difference": None, "detail": ""}
    if set(samples_n) != set(samples_u):
        res["detail"] = (f"the normalised and unnormalised tables name different samples "
                         f"({samples_n} vs {samples_u})")
        return res
    if not found:
        res["detail"] = ("directLFQ's normaliser reported no shifts (it did not run on the "
                         "samples)")
        return res
    if set(found) != set(samples_n) and samples_n:
        extra = sorted(set(samples_n) - set(found))
        if extra:
            res["detail"] = f"no shift found for samples {extra}: the normaliser did not see them"
            return res
    lg = math.log2

    def pairs(s, ref):
        out = []
        for p, nv in normalised.items():
            uv = unnormalised.get(p) or {}
            if s in nv and s in uv and ref in nv and ref in uv:
                out.append((lg(nv[s]) - lg(uv[s])) - (lg(nv[ref]) - lg(uv[ref])))
        return out

    # the reference: the sample sharing the most proteins with the others (ties: by name)
    def shared(s):
        return sum(1 for p, nv in normalised.items() if s in nv and s in (unnormalised.get(p) or {}))
    ref = max(samples_n, key=lambda s: (shared(s), [-ord(c) for c in s])) if samples_n else None
    res["reference"] = ref
    worst, bad = 0.0, []
    fref = found.get(ref) or 0.0
    for s in samples_n:
        f = None if found.get(s) is None else found[s] - fref
        a = 0.0 if s == ref else _median(pairs(s, ref))
        res["found"][s] = None if f is None else round(f, 6)
        res["applied"][s] = None if a is None else round(a, 6)
        if f is None and a is None:
            continue
        if f is None or a is None:
            bad.append(f"{s}: found {f}, applied {a} (no proteins in common with {ref})")
            continue
        d = abs(a - f)
        worst = max(worst, d)
        if d > tol:
            bad.append(f"{s}: shift found {f:+.4f} log2, applied {a:+.4f}")
    res["max_abs_difference"] = round(worst, 6)
    res["passed"] = not bad
    res["detail"] = ("shift applied = shift found for every sample (max |difference| "
                     f"{worst:.2g} log2, reference {ref})" if not bad else
                     "the normalisation directLFQ found was not applied: " + "; ".join(bad))
    return res


def guard(env, found, normalised, unnormalised, tol=GUARD_TOL):
    """THE directLFQ guard (Ruling 90): passed only when every check holds. -> {"passed",
    "checks", "refusal"}; `refusal` is the message a refused run stops with."""
    checks = env_checks(env)
    checks["shift_applied_equals_found"] = shift_check(found, normalised, unnormalised, tol)
    failed = [k for k, v in checks.items() if not v["passed"]]
    refusal = None
    if failed:
        refusal = ("REFUSED: directLFQ's result is not used -- its between-sample normalisation "
                   "cannot be shown to have run: " + "; ".join(
                       checks[k]["detail"] for k in failed)
                   + f". Use the pinned venv ({', '.join(DIRECTLFQ_PINS)}), built by {BUILDER}.")
    return {"passed": not failed, "checks": checks, "refusal": refusal}


# ------------------------------------------------- inside the venv (directLFQ) ---
def environment():
    """This python's directLFQ environment, after directLFQ's own copy check (what run_lfq()
    calls first)."""
    import platform
    import numpy
    import pandas
    import directlfq
    import directlfq.config as config
    try:
        import numba
        numba_v = numba.__version__
    except ImportError:
        numba_v = None
    config.check_wether_to_copy_numpy_arrays_derived_from_pandas()
    return {"directlfq": getattr(directlfq, "__version__", None), "pandas": pandas.__version__,
            "numpy": numpy.__version__, "numba": numba_v, "python": platform.python_version(),
            "copy_flag": config.COPY_NUMPY_ARRAYS_DERIVED_FROM_PANDAS}


def _install_spy(normalization):
    """Record, per sample, the shift directLFQ's between-sample normaliser computes, in its own
    frame (before 0.3.3's last line can hand back the original under pandas 3). -> (found list,
    undo)."""
    cls = normalization.NormalizationManager
    orig = cls._run_normalization
    target = normalization.NormalizationManagerSamplesOnSelectedProteins
    found = []

    def spy(self):
        if not isinstance(self, target):
            return orig(self)
        before = self.complete_dataframe.copy(deep=True)       # samples x ions, log2
        orig(self)
        med = (self.complete_dataframe - before).median(axis=1, skipna=True)
        found.append({str(k): (None if v != v else float(v)) for k, v in med.items()})
        return None

    cls._run_normalization = spy
    return found, (lambda: setattr(cls, "_run_normalization", orig))


def run_in_venv(out, cores=1):
    """`directlfq` subcommand: directLFQ on <out>/combined_ion.tsv twice (its normalisation on,
    then off), the guard, and the record <out>/directlfq/directlfq_record.json. -> exit code."""
    ion = os.path.join(out, ION_TABLE)
    if not os.path.isfile(ion):
        print(f"[fragpipe_quant] no {ION_TABLE} in {out}", file=sys.stderr)
        return 2
    d = os.path.join(out, DIRECTLFQ_DIR)
    os.makedirs(d, exist_ok=True)
    rec = {"status": "refused", "tool": "directLFQ", "input": ION_TABLE,
           "input_sha256": sha256_file(ion), "pins": list(DIRECTLFQ_PINS),
           "at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")}
    env = environment()
    rec["versions"] = {k: env[k] for k in ("directlfq", "pandas", "numpy", "numba", "python")}
    link = os.path.join(d, ION_TABLE)
    for p in glob.glob(link + "*"):                 # an earlier run's tables, never the input
        if os.path.isfile(p) or os.path.islink(p):
            os.remove(p)
    try:
        os.symlink(os.path.join("..", ION_TABLE), link)
    except OSError:
        shutil.copyfile(ion, link)
    import directlfq.config as config
    import directlfq.lfq_manager as lm
    import directlfq.normalization as normalization
    found, undo = _install_spy(normalization)
    try:
        lm.run_lfq(link, num_cores=cores, deactivate_normalization=False,
                   log_processed_proteins=False, compile_normalized_ion_table=False)
        copy_norm = config.COPY_NUMPY_ARRAYS_DERIVED_FROM_PANDAS
        lm.run_lfq(link, num_cores=cores, deactivate_normalization=True,
                   log_processed_proteins=False, compile_normalized_ion_table=False)
        copy_raw = config.COPY_NUMPY_ARRAYS_DERIVED_FROM_PANDAS
    finally:
        undo()
    env["copy_flag"] = copy_norm if copy_norm == copy_raw else True
    reformat = glob.glob(link + ".*.aq_reformat.tsv")
    rec["input_type"] = (re.match(re.escape(ION_TABLE) + r"\.(.+)\.aq_reformat\.tsv$",
                                  os.path.basename(reformat[0])).group(1) if reformat else None)
    norm_f, raw_f = link + ".protein_intensities.tsv", link + ".no_norm.protein_intensities.tsv"
    normalised, samples = read_protein_table(norm_f)
    unnormalised, _ = read_protein_table(raw_f)
    g = guard(env, found[0] if len(found) == 1 else {}, normalised, unnormalised)
    if len(found) != 1:
        g["checks"]["shift_applied_equals_found"]["detail"] += (
            f" (the normaliser ran {len(found)} times on the samples; expected once)")
    for p in reformat + [link]:                     # directLFQ's reformatted copy of the input
        try:
            os.remove(p)
        except OSError:
            pass
    rec.update(status="ok" if g["passed"] else "refused", samples=samples, guard=g,
               outputs={"normalised": f"{DIRECTLFQ_DIR}/{os.path.basename(norm_f)}",
                        "unnormalised": f"{DIRECTLFQ_DIR}/{os.path.basename(raw_f)}"},
               proteins_quantified=sum(1 for v in normalised.values() if v),
               proteins_in_every_sample=sum(1 for v in normalised.values()
                                            if len(v) == len(samples)))
    with open(os.path.join(d, DIRECTLFQ_RECORD + ".tmp"), "w") as fh:
        json.dump(rec, fh, indent=2)
    os.replace(os.path.join(d, DIRECTLFQ_RECORD + ".tmp"), os.path.join(d, DIRECTLFQ_RECORD))
    print(json.dumps({"status": rec["status"], "samples": samples,
                      "proteins_quantified": rec["proteins_quantified"],
                      "guard": {"passed": g["passed"], "refusal": g["refusal"]}}, indent=2))
    if not g["passed"]:
        print(f"[fragpipe_quant] {g['refusal']}", file=sys.stderr)
        return 3
    return 0


# ------------------------------------------------ from the adapter (pipeline env) ---
def load_directlfq_record(folder):
    try:
        with open(os.path.join(folder, DIRECTLFQ_DIR, DIRECTLFQ_RECORD)) as fh:
            rec = json.load(fh)
    except (OSError, ValueError):
        return None
    return rec if isinstance(rec, dict) else None


def _record_fits(rec, sha):
    g = (rec or {}).get("guard") or {}
    return bool(rec and rec.get("status") == "ok" and g.get("passed") is True
                and rec.get("input_sha256") == sha and rec.get("pins") == list(DIRECTLFQ_PINS))


def directlfq_resources(ion_bytes, n_experiments):
    """The compute node for one directLFQ run (two run_lfq calls: normalised and not). DEFAULTS, not
    measured beyond the two-run case (330 KB ion table, 2 runs: under 1 min, 330 MB): 1 CPU and
    8G for up to 3 experiments; else 4 CPUs, 8G or 12x the ion table, whichever is larger; 30
    minutes plus 1 per experiment, at most 12 h."""
    gb = ion_bytes / 1e9
    big = n_experiments > MIN_MAXLFQ_EXPERIMENTS
    return {"cpus": 4 if big else 1, "mem": f"{max(8, math.ceil(12 * gb))}G",
            "minutes": min(30 + n_experiments, 720), "source": "default"}


def run_directlfq(folder, experiments, runner=None):
    """directLFQ for this search: its record when one for THIS combined_ion.tsv passed the guard,
    else a run in the pinned venv (on a compute node: on_compute.py), then its record. Raises
    QuantError when the venv is missing, the run fails, the guard refuses, or an experiment is
    missing from directLFQ's table. `runner`: subprocess.run (tests inject one)."""
    ion = os.path.join(folder, ION_TABLE)
    if not os.path.isfile(ion):
        raise QuantError(f"no {ION_TABLE} in the FragPipe output: directLFQ (the quantity of a "
                         f"search whose IonQuant wrote no MaxLFQ for every experiment) reads "
                         f"IonQuant's ion table. Nothing was built.")
    sha = sha256_file(ion)
    rec = load_directlfq_record(folder)
    if not _record_fits(rec, sha):
        py, where = venv_python()
        if not py and rec and rec.get("input_sha256") == sha and rec.get("status") == "refused":
            # directLFQ's last run on THIS table was refused, and there is no venv to run again
            raise QuantError(((rec.get("guard") or {}).get("refusal") or "REFUSED: directLFQ's "
                              "last run on this combined_ion.tsv was refused")
                             + f" (its record, {DIRECTLFQ_DIR}/{DIRECTLFQ_RECORD}). No pinned venv "
                             f"to run it again ({where}): {builder_command()}. Nothing was built.")
        if not py:
            raise QuantError(
                f"directLFQ's pinned venv is not installed ({where}). Without IonQuant's "
                f"MaxLFQ for every experiment, directLFQ is this search's protein quantity. "
                f"Build the venv once, in a job: {builder_command()} -- then run "
                f"run_search.py --adapt-only again. Nothing was built.")
        res = directlfq_resources(os.path.getsize(ion), len(experiments))
        cmd = [sys.executable or "python3", os.path.join(HERE, "on_compute.py"),
               "--time", str(res["minutes"]), "--mem", res["mem"], "--cpus", str(res["cpus"]),
               "--name", "directlfq", "--why", "directLFQ on combined_ion.tsv", "--",
               py, "-I", os.path.abspath(__file__), "directlfq", "--out",
               os.path.abspath(folder), "--cores", str(res["cpus"])]
        p = (runner or subprocess.run)(cmd, capture_output=True, text=True)
        rec = load_directlfq_record(folder)
        if p.returncode not in (0, 3) or not rec or rec.get("input_sha256") != sha:
            raise QuantError(f"directLFQ did not run (exit {p.returncode}, venv {where}): "
                             f"{(p.stderr or p.stdout or '').strip()[-1500:]}. Nothing was built.")
        if not _record_fits(rec, sha):
            raise QuantError(((rec.get("guard") or {}).get("refusal")
                              or f"directLFQ's record says {rec.get('status')}")
                             + " Nothing was built.")
    missing = [e for e in experiments if e not in (rec.get("samples") or [])]
    if missing:
        raise QuantError(f"directLFQ's table has no column for {', '.join(missing)} (its samples: "
                         f"{', '.join(rec.get('samples') or [])}): IonQuant's ion table left "
                         f"those runs out. Nothing was built.")
    return rec


def directlfq_values(folder, rec):
    """{(protein, sample): directLFQ intensity} from the normalised table the record names."""
    table, _ = read_protein_table(os.path.join(folder, rec["outputs"]["normalised"]))
    return {(p, s): v for p, vals in table.items() for s, v in vals.items()}


def guard_sentence(rec):
    """One line for the declaration: the directLFQ guard's verdict and its evidence."""
    if not rec:
        return None
    g = rec.get("guard") or {}
    checks = g.get("checks") or {}
    sc = checks.get("shift_applied_equals_found") or {}
    v = rec.get("versions") or {}

    def held(name):
        """each claim from its own check (review L3): True, False, or None (not recorded)"""
        return (checks.get(name) or {}).get("passed")

    pandas = {True: "(< 3)", False: "(NOT < 3)"}.get(held("pandas_below_3"), "(not checked)")
    copy = {True: "directLFQ's copy flag False",
            False: "directLFQ's copy flag NOT False"}.get(held("copy_flag_false"),
                                                          "directLFQ's copy flag not checked")
    shift = {True: "applied as found", False: "NOT applied as found"}.get(
        held("shift_applied_equals_found"), "not checked")
    return (f"{'passed' if g.get('passed') else 'REFUSED'}: pandas {v.get('pandas')} {pandas}; "
            f"{copy}; its between-sample normalisation {shift} "
            f"(max |applied - found| {sc.get('max_abs_difference')} log2 over "
            f"{len(sc.get('found') or {})} samples)")


_ABS_PATH = re.compile(r"(?<![\w.])/[^\s'\"():;,]+")


def path_free(text):
    """`text` with each absolute path cut to its last name: a record that can be delivered
    (fragpipe_quantity.json is in SEARCH_FILES) never carries the Core's paths; the full text
    goes to the person on stderr."""
    return _ABS_PATH.sub(lambda m: os.path.basename(m.group(0).rstrip("/")) or "/",
                         str(text or ""))


def record(folder, quantity, directlfq=None):
    """<folder>/fragpipe_quantity.json + search_provenance.json `fragpipe_quantity` (merged):
    the quantity object, plus the directLFQ record's verdict and versions when it ran. No paths
    beyond the search folder's own file names. -> the record."""
    rec = dict(quantity)
    if rec.get("refusal"):
        rec["refusal"] = path_free(rec["refusal"])
    if directlfq:
        rec["directlfq"] = {k: directlfq.get(k) for k in (
            "status", "versions", "input", "input_sha256", "input_type", "samples", "outputs",
            "proteins_quantified", "proteins_in_every_sample", "at")}
        rec["directlfq"]["guard"] = {"passed": (directlfq.get("guard") or {}).get("passed"),
                                     "sentence": guard_sentence(directlfq)}
    import sage_lfq_check                        # the one writer of a search record + provenance
    sage_lfq_check.write_record(folder, RECORD, PROVENANCE_KEY, rec)
    return rec


# --------------------------------------------------------------------- the CLI ---
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("directlfq", help="run directLFQ + the guard (inside the pinned venv)")
    d.add_argument("--out", required=True, help="the FragPipe output folder (combined_ion.tsv)")
    d.add_argument("--cores", type=int,
                   default=int(os.environ.get("SLURM_CPUS_PER_TASK") or 1))
    sub.add_parser("check-env", help="guard checks 1-2 of this python, as JSON")
    sub.add_parser("pins", help="the venv's requirements, one per line")
    v = sub.add_parser("venv", help="where the pinned venv is, or why there is none")
    v.add_argument("--default-prefix", action="store_true",
                   help="print where build_directlfq_venv.sh builds without --prefix")
    wm = sub.add_parser("write-marker", help=f"(the builder) write {VENV_MARKER}")
    wm.add_argument("--prefix", required=True)
    wm.add_argument("--base-python", required=True)
    wm.add_argument("--freeze", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "pins":
        print("\n".join(DIRECTLFQ_PINS))
        return 0
    if a.cmd == "venv":
        if a.default_prefix:
            print(default_prefix())
            return 0
        py, where = venv_python()
        print(json.dumps({"python": py, "where" if py else "why": where}, indent=2))
        return 0 if py else 1
    if a.cmd == "check-env":
        env = environment()
        checks = env_checks(env)
        ok = all(c["passed"] for c in checks.values())
        print(json.dumps({"passed": ok, "env": env, "checks": checks}, indent=2))
        return 0 if ok else 3
    if a.cmd == "write-marker":
        env = environment()
        checks = env_checks(env)
        if not all(c["passed"] for c in checks.values()):
            print(json.dumps(checks, indent=2), file=sys.stderr)
            return 3
        marker = {"pins": list(DIRECTLFQ_PINS), "versions": env, "checks": checks,
                  "base_python": a.base_python, "freeze": os.path.basename(a.freeze),
                  "built": datetime.datetime.now(datetime.timezone.utc).isoformat(
                      timespec="seconds"),
                  "builder": BUILDER}
        path = os.path.join(a.prefix, VENV_MARKER)
        with open(path + ".tmp", "w") as fh:
            json.dump(marker, fh, indent=2)
        os.replace(path + ".tmp", path)
        print(path)
        return 0
    return run_in_venv(a.out, a.cores)


if __name__ == "__main__":
    sys.exit(main())
