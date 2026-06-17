#!/usr/bin/env python3
"""
collect_conditions.py  --  Build, MAP, and validate the experimental-design
metadata the DE step needs.

The user provides conditions however is easiest — they tell the agent in words
("the first three are control, the rest treated") or upload a file (any CSV/TSV
with a sample column + a group column, however named). The agent turns that into
an intent, and `--map` does the deterministic part: fuzzy-match each named sample
to the REAL raw filenames, then report exactly what's ambiguous so the agent can
confirm only those with the user. Filename matching is grounded in the actual
runs — never guessed.

metadata CSV schema:  File.Name,Group[,Batch,Covariate1,Covariate2]
  File.Name must match the Run names in the DIA-NN report (or raw basenames).

Modes:
  # list the real run names the agent must map to
  python3 collect_conditions.py --list-runs --from-dir /data --glob '*.d'

  # MAP an uploaded conditions file onto the real runs
  python3 collect_conditions.py --map conditions.csv --from-dir /data --glob '*.d' \
          --from-file user_conditions.csv

  # MAP an agent-built mapping (from the user's free-text description) onto runs
  python3 collect_conditions.py --map conditions.csv --from-report report.parquet \
          --mapping-json '{"groups": {"control": ["A1","A2"], "treated": ["B1","B2"]}}'

  # emit a blank template (fallback when there's nothing to map)
  python3 collect_conditions.py --emit-template conditions.csv --from-dir /data --glob '*.d'

  # validate a finished design against the search output
  python3 collect_conditions.py --validate conditions.csv --against report.parquet
"""
import sys, os, csv, glob, json, re, argparse
from collections import Counter, defaultdict

COV_COLS = ["Batch", "Covariate1", "Covariate2"]
SAMPLE_HEADERS = {"file.name", "filename", "file", "run", "sample", "sample name",
                  "samplename", "name", "raw", "raw file", "rawfile", "id"}
GROUP_HEADERS = {"group", "condition", "treatment", "class", "type", "category",
                 "cohort", "phenotype"}
BATCH_HEADERS = {"batch", "block", "plate", "run order", "runorder"}


# ----------------------------------------------------------------- run lists --
def runs_from_report(path):
    if path.endswith(".parquet"):
        try:
            import pyarrow.parquet as pq
            t = pq.read_table(path, columns=["Run"])
            return sorted(set(t.column("Run").to_pylist()))
        except Exception as e:
            sys.exit(f"Could not read Run column from parquet: {e}\nInstall pyarrow, or export a TSV report.")
    seen = set()
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if "Run" not in (rd.fieldnames or []):
            sys.exit("No 'Run' column in report.")
        for row in rd:
            seen.add(row["Run"])
    return sorted(seen)


def runs_from_dir(d, pattern):
    files = sorted(glob.glob(os.path.join(d, pattern)))
    if not files:
        sys.exit(f"No files matched {pattern} in {d}")
    return [os.path.splitext(os.path.basename(f.rstrip("/")))[0] for f in files]


def get_runs(a):
    if a.from_report:  return runs_from_report(a.from_report)
    if a.from_dir:     return runs_from_dir(a.from_dir, a.glob)
    if a.runs:         return [r.strip() for r in a.runs.split(",") if r.strip()]
    sys.exit("need --from-report, --from-dir, or --runs to know the real run names")


# ------------------------------------------------------------------ matching --
def norm(s):
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())


def match_to_runs(identifier, runs, run_norms):
    """Return runs matching an identifier: exact (ci) > substring either way."""
    idn = norm(identifier)
    if not idn:
        return []
    exact = [r for r, rn in zip(runs, run_norms) if rn == idn]
    if exact:
        return exact
    return [r for r, rn in zip(runs, run_norms) if idn and (idn in rn or rn in idn)]


# ---------------------------------------------------------- conditions source --
def parse_conditions_file(path):
    """Parse an uploaded CSV/TSV: detect sample + group (+batch/covariate) columns
    however they're named. Returns list of dicts: {sample, group, extras{}}."""
    delim = "\t" if path.lower().endswith((".tsv", ".txt")) else None
    with open(path, newline="") as fh:
        sample = fh.read(4096); fh.seek(0)
        if delim is None:
            try: delim = csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
            except csv.Error: delim = ","
        rd = csv.DictReader(fh, delimiter=delim)
        headers = rd.fieldnames or []
        hmap = {h: norm(h) for h in headers}
        scol = next((h for h in headers if hmap[h] in {norm(x) for x in SAMPLE_HEADERS}), None)
        gcol = next((h for h in headers if hmap[h] in {norm(x) for x in GROUP_HEADERS}), None)
        bcol = next((h for h in headers if hmap[h] in {norm(x) for x in BATCH_HEADERS}), None)
        if scol is None or gcol is None:
            sys.exit(json.dumps({"error": "could not find sample and group columns",
                                 "headers": headers,
                                 "hint": "rename a column to 'sample'/'file' and 'group'/'condition', "
                                         "or have the agent pass --mapping-json instead"}))
        other = [h for h in headers if h not in (scol, gcol, bcol) and h]
        out = []
        for row in rd:
            s = (row.get(scol) or "").strip()
            g = (row.get(gcol) or "").strip()
            if not s:
                continue
            extras = {}
            if bcol and (row.get(bcol) or "").strip():
                extras["Batch"] = row[bcol].strip()
            for i, h in enumerate(other[:2]):
                if (row.get(h) or "").strip():
                    extras[COV_COLS[i + 1]] = row[h].strip()
            out.append({"sample": s, "group": g, "extras": extras})
        return out


def intent_from_json(blob):
    """Accept {'mapping': {sample: group}} or {'groups': {group: [samples]}}."""
    d = json.loads(blob)
    out = []
    if "mapping" in d:
        for s, g in d["mapping"].items():
            out.append({"sample": s, "group": g, "extras": {}})
    elif "groups" in d:
        for g, samples in d["groups"].items():
            for s in samples:
                out.append({"sample": s, "group": g, "extras": {}})
    else:
        sys.exit("--mapping-json must have a 'mapping' or 'groups' key")
    return out


# ----------------------------------------------------------------------- map --
def do_map(out_path, runs, intent):
    run_norms = [norm(r) for r in runs]
    run_groups = defaultdict(set)      # run -> {groups}
    run_extras = {}                    # run -> extras
    multi_match = {}                   # identifier -> [runs] (one label, several files)
    unmatched = []                     # identifiers matching no run

    for item in intent:
        hits = match_to_runs(item["sample"], runs, run_norms)
        if not hits:
            unmatched.append(item["sample"]); continue
        if len(hits) > 1:
            multi_match[item["sample"]] = hits
        for r in hits:
            run_groups[r].add(item["group"])
            if item["extras"]:
                run_extras.setdefault(r, {}).update(item["extras"])

    assigned, conflicting = {}, {}
    for r in runs:
        gs = run_groups.get(r, set())
        if len(gs) == 1:
            assigned[r] = next(iter(gs))
        elif len(gs) > 1:
            conflicting[r] = sorted(gs)
    unassigned = [r for r in runs if r not in run_groups]

    # write the proposed CSV for the unambiguous part (every run gets a row;
    # ambiguous ones get a blank Group so the agent fills it after confirming)
    cov_used = sorted({k for e in run_extras.values() for k in e}, key=lambda c: COV_COLS.index(c) if c in COV_COLS else 99)
    cols = ["File.Name", "Group"] + cov_used
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in runs:
            row = [r, assigned.get(r, "")]
            for c in cov_used:
                row.append(run_extras.get(r, {}).get(c, ""))
            w.writerow(row)

    sizes = Counter(assigned.values())
    singletons = [g for g, c in sizes.items() if c < 2]
    needs_conf = bool(unassigned or conflicting or unmatched or singletons)
    print(json.dumps({
        "proposed_csv": os.path.abspath(out_path),
        "n_runs": len(runs),
        "assigned": assigned,
        "groups": dict(sizes),
        "ambiguities": {
            "unassigned_runs": unassigned,
            "conflicting_runs": conflicting,
            "unmatched_identifiers": unmatched,
            "multi_match_identifiers": multi_match,
            "singleton_groups": singletons,
        },
        "needs_confirmation": needs_conf,
        "guidance": "Confirm every item under 'ambiguities' with the user, then finalize "
                    "the CSV and run --validate. Do NOT proceed to a search while runs are "
                    "unassigned or conflicting.",
    }, indent=2))


# ------------------------------------------------------------------ template --
def emit_template(out, names, covariates):
    cols = ["File.Name", "Group"] + covariates
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for n in names:
            w.writerow([n] + [""] * (len(cols) - 1))
    print(json.dumps({"template": out, "n_samples": len(names), "columns": cols, "samples": names}, indent=2))


# ------------------------------------------------------------------ validate --
def validate(meta_path, report_path):
    with open(meta_path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    problems = []
    if not rows: problems.append("metadata is empty")
    cols = rows[0].keys() if rows else []
    for req in ("File.Name", "Group"):
        if req not in cols: problems.append(f"missing required column '{req}'")
    blanks = [r["File.Name"] for r in rows if not r.get("Group", "").strip()]
    if blanks: problems.append(f"{len(blanks)} rows have no Group: {blanks[:5]}...")
    sizes = Counter(r["Group"].strip() for r in rows if r.get("Group", "").strip())
    singletons = [g for g, c in sizes.items() if c < 2]
    if singletons: problems.append(f"groups with <2 replicates (no within-group variance): {singletons}")
    if report_path:
        report_runs = set(runs_from_report(report_path))
        meta_runs = {r["File.Name"] for r in rows}
        missing = report_runs - meta_runs
        extra = meta_runs - report_runs
        if missing: problems.append(f"{len(missing)} report runs missing from metadata: {sorted(missing)[:5]}...")
        if extra:   problems.append(f"{len(extra)} metadata rows not in report: {sorted(extra)[:5]}...")
    ok = not problems
    print(json.dumps({"valid": ok, "n_samples": len(rows), "groups": dict(sizes), "problems": problems}, indent=2))
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--list-runs", action="store_true", help="print the real run names and exit")
    ap.add_argument("--emit-template")
    ap.add_argument("--map", dest="map_out", help="write a proposed conditions.csv by mapping intent onto runs")
    ap.add_argument("--from-file", help="user-uploaded conditions file (CSV/TSV, any column names)")
    ap.add_argument("--mapping-json", help="agent-built intent: {'groups':{g:[samples]}} or {'mapping':{sample:g}}")
    ap.add_argument("--from-report")
    ap.add_argument("--from-dir")
    ap.add_argument("--glob", default="*.raw")
    ap.add_argument("--runs", help="comma-separated run names (alternative to --from-*)")
    ap.add_argument("--covariates", default="", help="comma list for --emit-template, e.g. Batch,Covariate1")
    ap.add_argument("--validate")
    ap.add_argument("--against", help="report to validate File.Name against")
    a = ap.parse_args()

    if a.list_runs:
        runs = get_runs(a)
        print(json.dumps({"n_runs": len(runs), "runs": runs}, indent=2))
    elif a.map_out:
        runs = get_runs(a)
        if a.from_file:      intent = parse_conditions_file(a.from_file)
        elif a.mapping_json: intent = intent_from_json(a.mapping_json)
        else: sys.exit("--map needs --from-file (uploaded file) or --mapping-json (agent intent)")
        do_map(a.map_out, runs, intent)
    elif a.emit_template:
        names = get_runs(a)
        covs = [c.strip() for c in a.covariates.split(",") if c.strip()]
        emit_template(a.emit_template, names, covs)
    elif a.validate:
        validate(a.validate, a.against)
    else:
        ap.print_help(); sys.exit(2)


if __name__ == "__main__":
    main()
