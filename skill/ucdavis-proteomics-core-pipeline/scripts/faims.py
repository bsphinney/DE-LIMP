#!/usr/bin/env python3
"""
faims.py -- which runs were acquired with FAIMS, and what a mix of them does to a FragPipe search.

The ONE definition of "this run has FAIMS" (project rule 3): detect_acquisition.py records it per
run, run_search.py decides on it, make_methods.py / session.py say it. Nothing else parses a
compensation voltage.

Where a run says so (checked 2026-11-11 against real Orbitrap mzML on HIVE, converted from Thermo
.raw by msconvert/ThermoRawFileParser):
  * mzML: every scan carries  <cvParam accession="MS:1001581" name="FAIMS compensation voltage"
    value="-55.0"/>  and its filter string (MS:1000512) reads "FTMS + p NSI cv=-55.00 Full ms
    [400.0000-1600.0000]". A run without FAIMS has neither: "FTMS + p NSI Full ms [375.0000-1500.0000]".
  * Thermo .raw: the same "cv=-55.00" token in the scan filter string (ThermoRawFileParser's
    PROXI/mzML output of the filter), read by detect_acquisition.thermo_isolation_windows.
  * Bruker .d (timsTOF): no FAIMS (ion mobility is TIMS); not applicable, None.

Why it matters: FragPipe 24.0's IonQuant indexes each run by its compensation voltage and logs
"<run> does not have FAIMS. If there are any other runs having FAIMS, IonQuant will crash." then
dies with "Cannot find CV -50 from all indexed CVs: [0]" -- after MSFragger, MSBooster, Percolator,
ProteinProphet and Philosopher have all finished (a 2.6 h, 16-CPU search lost, 2026-11-11). So a
search that quantifies with IonQuant must never be handed a mix. Stdlib only.

A run's status is one of:
  []                  spectra were read and none carries a compensation voltage -- no FAIMS
  [-50.0, -70.0]      the voltages seen (sorted), FAIMS
  None                not readable here (no spectrum in the slice, a .raw nobody read, a vendor
                      folder, a missing file) -- never guessed as "no FAIMS"
"""
import gzip
import json
import os
import re

# MS:1001581 is the PSI-MS term "FAIMS compensation voltage"
CV_ACCESSION = "MS:1001581"
HEAD_BYTES = (1_000_000, 8_000_000)   # read the first 1 MB; only if it holds no spectrum, 8 MB
# The most statuses() may read on a cluster LOGIN node (HPC@UCD's ucd-hpc rules: no I/O-heavy work
# there), counted at the worst case of HEAD_BYTES[-1] a run: 40 MB, i.e. five runs -- the five .d
# run_search.py's TDF_LOGIN_MAX reads there. More, and run_search.py moves the whole command to a
# compute node to read them (fragpipe_faims_gate), never marks the rest "unknown": an unknown run
# only gets a NOTE from the plain gate, so a mix hidden in it would crash IonQuant hours later.
LOGIN_READ_BYTES = 5 * HEAD_BYTES[-1]
MZML_SUFFIXES = (".mzml", ".mzml.gz")
GROUP_NONE = "no FAIMS"
GROUP_FAIMS = "FAIMS"
SPLIT_FLAG = "--faims-split"
GROUP_FLAG = "--faims-group"          # set by --faims-split on each child search; not for users
SPLIT_RECORD = "faims_split.json"

# an mzML cvParam, in either attribute order, and Thermo's "cv=-55.00" filter token
_CVPARAM_RE = re.compile(
    rb'<cvParam\b[^>]*?accession="' + CV_ACCESSION.encode() + rb'"[^>]*?>', re.S)
_VALUE_RE = re.compile(rb'\bvalue="\s*(-?[0-9]+(?:\.[0-9]*)?)\s*"')
_FILTER_CV_RE = re.compile(r"(?<![A-Za-z0-9])cv=(-?[0-9]+(?:\.[0-9]+)?)(?![0-9])")
_SPECTRUM_RE = re.compile(rb"<spectrum\b")


def cvs_from_filter(filt):
    """The compensation voltages in a Thermo scan filter string ("FTMS + p NSI cv=-55.00 Full ms
    ..."), as floats, sorted and distinct; [] for a filter without one (or not a string)."""
    if not isinstance(filt, str):
        return []
    return sorted({float(v) for v in _FILTER_CV_RE.findall(filt)})


def cvs_from_mzml_bytes(blob):
    """(cvs, saw_spectrum) from the head of an mzML as bytes: the MS:1001581 values and the
    voltages in the filter strings, sorted and distinct, and whether any <spectrum> was in it
    (a head with none says nothing about FAIMS). Works on a slice cut mid-element."""
    found = set()
    for tag in _CVPARAM_RE.findall(blob):
        m = _VALUE_RE.search(tag)
        if m:
            found.add(float(m.group(1)))
    for m in re.finditer(rb'accession="MS:1000512"[^>]*?value="([^"]*)"', blob):
        found.update(cvs_from_filter(m.group(1).decode("utf-8", "replace")))
    return sorted(found), bool(_SPECTRUM_RE.search(blob))


def run_cvs(path, head_bytes=HEAD_BYTES):
    """The FAIMS status of ONE run (see the module doc): [] / [cv, ...] / None. Only an mzML
    (.mzML, .mzML.gz) is read here -- a Thermo .raw needs ThermoRawFileParser (detect_acquisition
    reads it and records `faims_cvs`; use from_records); a Bruker .d has no FAIMS (None). Reads
    the run's head only (head_bytes, smallest first): a login node reads a MB a run, not the run."""
    p = str(path).rstrip("/")
    if not p.lower().endswith(MZML_SUFFIXES):
        return None
    opn = gzip.open if p.lower().endswith(".gz") else open
    for n in ((head_bytes,) if isinstance(head_bytes, int) else head_bytes):
        try:
            with opn(p, "rb") as fh:
                blob = fh.read(n)
        except (OSError, EOFError):
            return None
        cvs, saw = cvs_from_mzml_bytes(blob)
        if cvs:
            return cvs
        if saw:
            return []
    return None


def from_records(records):
    """{file: faims_cvs} from detect_acquisition's per-file records (acquisition.json `files`),
    for files whose record carries the field (a record from before 2.15 does not)."""
    out = {}
    for r in records if isinstance(records, list) else []:
        if isinstance(r, dict) and isinstance(r.get("file"), str) and "faims_cvs" in r:
            v = r["faims_cvs"]
            out[os.path.realpath(r["file"].rstrip("/"))] = v if isinstance(v, list) else None
    return out


def load_records(acquisition_json):
    """from_records of an acquisition.json path; {} when absent or unreadable."""
    if not acquisition_json:
        return {}
    try:
        with open(acquisition_json) as fh:
            return from_records((json.load(fh) or {}).get("files"))
    except (OSError, ValueError, AttributeError):
        return {}


def to_read(files, recorded=None):
    """The runs statuses() would read here: an mzML with no record from step 2 (a .raw or .d is
    never read by this module)."""
    recorded = recorded or {}
    return [f for f in files if str(f).rstrip("/").lower().endswith(MZML_SUFFIXES)
            and os.path.realpath(str(f).rstrip("/")) not in recorded]


def read_bytes_max(files, recorded=None):
    """The most statuses() reads for `files`: HEAD_BYTES[-1] for each run it would read."""
    return len(to_read(files, recorded)) * HEAD_BYTES[-1]


def over_login_budget(files, recorded=None):
    """True when reading `files` here could pass LOGIN_READ_BYTES -- the one budget."""
    return read_bytes_max(files, recorded) > LOGIN_READ_BYTES


def statuses(files, recorded=None, reader=run_cvs):
    """{file: status} -- the record step 2 left for the file when there is one, else the run read
    here (`reader`, run_cvs). One rule for every caller."""
    recorded = recorded or {}
    out = {}
    for f in files:
        key = os.path.realpath(str(f).rstrip("/"))
        out[f] = recorded[key] if key in recorded else reader(f)
    return out


def group(status_by_file):
    """{"none": [files], "faims": [files], "unknown": [files]} from statuses()."""
    g = {"none": [], "faims": [], "unknown": []}
    for f, s in status_by_file.items():
        g["unknown" if s is None else "faims" if s else "none"].append(f)
    return g


def is_mixed(grouped):
    """True when FAIMS and non-FAIMS runs are both present (IonQuant crashes on that)."""
    return bool(grouped["faims"]) and bool(grouped["none"])


def all_cvs(status_by_file):
    """Every compensation voltage seen across the runs, sorted."""
    return sorted({c for s in status_by_file.values() if s for c in s})


def mix_message(grouped, cvs, quantifier="IonQuant"):
    """The named refusal for a quantifying FragPipe search handed a mix: the counts and the ways
    out. Words only; run_search.py exits with it."""
    n_f, n_n, n_u = (len(grouped[k]) for k in ("faims", "none", "unknown"))
    unk = (f" ({n_u} more could not be read and are not counted)" if n_u else "")
    return (f"FAIMS and non-FAIMS runs are mixed in one {quantifier} quantification: "
            f"{n_f} run(s) were acquired with FAIMS (compensation voltage "
            f"{', '.join(f'{c:g}' for c in cvs) or 'unknown'} V) and {n_n} without{unk}. "
            f"FragPipe's IonQuant indexes runs by compensation voltage and crashes on a mix "
            f"(\"Cannot find CV -50 from all indexed CVs: [0]\") after the whole search has run "
            f"-- hours of compute lost. Choose one: (1) {SPLIT_FLAG} -- one search and one "
            f"quantification per group ({GROUP_FAIMS} {n_f} runs, {GROUP_NONE} {n_n} runs), "
            f"each in its own folder; quantities are NOT comparable across the groups, and the "
            f"Methods, the report and MANIFEST.txt say so; (2) search only one group's runs; "
            f"(3) an identifications-only workflow (quantify off), which has no IonQuant step "
            f"to crash.")


def split_record(group_name, cvs, n_runs, other_counts):
    """The self-describing record one group's search carries (search_provenance.json `faims`
    and faims_split.json): rule 1 -- consumers read this, they do not re-derive it."""
    return {"group": group_name, "cvs": list(cvs), "n_runs": n_runs,
            "other_groups": dict(other_counts), "quantified_separately": True,
            "comparable_across_groups": False,
            "sentence": methods_sentence(group_name, n_runs, other_counts)}


def methods_sentence(group_name, n_runs, other_counts):
    """The one sentence the Methods, the report and MANIFEST.txt use for a split search."""
    others = "; ".join(f"{n} {name} run(s)" for name, n in other_counts.items())
    return (f"This dataset mixed FAIMS and non-FAIMS acquisitions, so FragPipe was run "
            f"separately for each group; this search holds the {n_runs} {group_name} run(s) "
            f"(the other group(s): {others}). Peptide and protein quantities were computed "
            f"within each group only and are not comparable across the groups.")


def split_plan(grouped, cvs_by_file):
    """[(group_name, key, files, cvs)] for --faims-split; empty groups left out."""
    plan = []
    for key, name in (("none", GROUP_NONE), ("faims", GROUP_FAIMS)):
        files = grouped[key]
        if files:
            cvs = sorted({c for f in files for c in (cvs_by_file.get(f) or [])})
            plan.append((name, key, files, cvs))
    return plan


def design_check(plan, conditions_path):
    """Is a condition of the session's design wholly inside one FAIMS group (C2 review stat M4)?
    Then its comparison with a condition measured in the other group would be a FAIMS difference,
    not a biological one -- each group is searched and quantified on its own. A WARNING, never a
    refusal: the split is still the only way to quantify a mix.
    plan: split_plan()'s [(group_name, key, files, cvs)]. conditions_path: the session's
    conditions.csv (File.Name, Group), or None when the session has none yet.
    -> {"checked", "why_not", "conditions", "table" {condition: {group: n_runs}},
        "confounded" [{"condition", "group"}], "runs_not_in_design", "warning"}.
    Never raises: anything that stops the check is said in why_not, never read as "no confound"."""
    res = {"checked": False, "why_not": None, "conditions": conditions_path, "table": {},
           "confounded": [], "runs_not_in_design": [], "warning": None}
    if not conditions_path or not os.path.isfile(conditions_path):
        res["why_not"] = ("the session's design (input/conditions.csv) did not exist yet when the "
                          "search was split, so whether a condition lies wholly in one FAIMS group "
                          "was NOT checked; check it before any comparison")
        return res
    try:
        import csv
        import collect_conditions                     # run_name: THE run name of a raw file
        with open(conditions_path, encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.DictReader(fh))
        group_of = {collect_conditions.run_name(f): name
                    for name, _key, files, _cvs in plan for f in files}
        for r in rows:
            run = collect_conditions.run_name((r.get("File.Name") or "").strip())
            cond = (r.get("Group") or "").strip()
            if not run or not cond:
                continue
            g = group_of.get(run)
            if g is None:
                continue
            res["table"].setdefault(cond, {}).setdefault(g, 0)
            res["table"][cond][g] += 1
        named = {collect_conditions.run_name((r.get("File.Name") or "").strip()) for r in rows}
        res["runs_not_in_design"] = sorted(x for x in group_of if x not in named)
    except Exception as e:                            # noqa: BLE001 -- said, never "no confound"
        res["why_not"] = (f"{conditions_path} could not be read ({type(e).__name__}: {e}), so "
                          "whether a condition lies wholly in one FAIMS group was NOT checked")
        return res
    if not res["table"]:
        res["why_not"] = (f"no run of this search is named in {conditions_path}, so whether a "
                          "condition lies wholly in one FAIMS group was NOT checked")
        return res
    res["checked"] = True
    groups = [name for name, *_ in plan]
    res["confounded"] = [{"condition": c, "group": next(iter(t))}
                         for c, t in sorted(res["table"].items()) if len(t) == 1 and len(groups) > 1]
    if res["confounded"]:
        res["warning"] = (
            "WARNING: " + "; ".join(f"condition {x['condition']!r} lies wholly in the "
                                    f"{x['group']} group" for x in res["confounded"])
            + ". Each FAIMS group is searched and quantified on its own, so comparing such a "
              "condition with one measured in the other group would be a FAIMS difference, not "
              "a biological one: compare conditions only within one group, and never combine "
              "the two groups' quantities.")
    return res


def design_table_text(check):
    """The condition-by-group table design_check() found, as lines for the terminal."""
    groups = sorted({g for t in check.get("table", {}).values() for g in t})
    lines = ["condition | " + " | ".join(groups)]
    for c, t in sorted(check.get("table", {}).items()):
        lines.append(f"{c} | " + " | ".join(str(t.get(g, 0)) for g in groups))
    return lines

