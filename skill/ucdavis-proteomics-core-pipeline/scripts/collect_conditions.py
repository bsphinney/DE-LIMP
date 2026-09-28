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

metadata CSV schema:  File.Name,Group[,Batch,Covariate1,Covariate2][,<subject column>]
  File.Name must match the Run names in the DIA-NN report (or raw basenames).
  A subject / animal column (Mouse, Animal, Subject, Patient, Donor ...) is kept UNDER ITS
  OWN NAME when its values look like one (>= 3 subjects, not just the groups relabelled, and
  -- where subjects span groups -- genuinely repeating across them); run DE with
  --block <that name>. A header that says "subject" over values that do not (Subject = M/F,
  Patient = Yes/No, Donor identical to Group) is NOT trusted: it stays a covariate and
  --map asks for confirmation (subject_ambiguous). --subject-column <header> confirms one.

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
# "block" is NOT a batch header: filed as Batch, a sample sheet's Block column became the
# fixed covariate Batch, and `--block Batch` is then refused (a column cannot be both). It
# names a blocking unit, so it goes through the subject detection below instead.
BATCH_HEADERS = {"batch", "plate", "run order", "runorder"}
# The unit samples come from. Kept under its own name so run_de.R can fit it as a block
# (--block); as Covariate1 it was a fixed effect, and nested in the groups (mice within
# age) that makes the design rank-deficient. Matched after dropping an id/number suffix,
# so "Mouse ID", "animal_no" and "Patient #" count too. Its VALUES still have to look like
# subjects (subject_assessment) or it stays a covariate and --map asks.
SUBJECT_HEADERS = {"mouse", "mice", "animal", "rat", "subject", "patient", "donor",
                   "individual", "participant", "pair", "block"}


def subject_header(h):
    """True when a column header names the animal / subject / patient samples came from."""
    n = norm(h)
    return n in SUBJECT_HEADERS or re.sub(r"(id|no|num|number)$", "", n) in SUBJECT_HEADERS


def subject_assessment(subject_of, group_of):
    """Why a subject column's VALUES do not look like the animal / patient runs came from
    (empty list = they do). The header alone is not evidence: 'Subject' can hold sex,
    'Patient' yes/no, 'Donor' the group itself -- and a wrong block silently drops that
    factor from the model."""
    runs = [r for r in subject_of if r in group_of]
    subjects = {subject_of[r] for r in runs}
    reasons = []
    if len(subjects) < 3:
        reasons.append(f"only {len(subjects)} distinct value(s) ({', '.join(sorted(subjects))}): "
                       f"a category such as sex or yes/no, not an animal / patient identifier")
    groups_of = defaultdict(set)
    subjects_of = defaultdict(set)
    for r in runs:
        groups_of[subject_of[r]].add(group_of[r])
        subjects_of[group_of[r]].add(subject_of[r])
    if runs and all(len(g) == 1 for g in groups_of.values()) and \
            all(len(v) == 1 for v in subjects_of.values()):
        reasons.append("its values are the groups relabelled (one value per group): as a block "
                       "it would be the grouping itself")
    spanning = [k for k, g in groups_of.items() if len(g) > 1]
    if spanning and len(spanning) * 2 < len(groups_of):
        reasons.append(f"only {len(spanning)} of {len(groups_of)} values recur across groups: "
                       f"not a subject measured under several conditions")
    return reasons


def column_name(h):
    """A header as a metadata column name run_de.R's --block can take: 'Mouse ID' -> 'Mouse_ID'."""
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9_.]", "_", h.strip())).strip("_") or "Subject"


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


def run_name(path):
    """The run name a raw file gets in File.Name: its basename without the extension.
    Defined once -- core_submission.py imports it rather than re-deriving it, so a
    submission's conditions.csv can never disagree with one built from --from-dir."""
    return os.path.splitext(os.path.basename(str(path).rstrip("/")))[0]


def runs_from_dir(d, pattern):
    files = sorted(glob.glob(os.path.join(d, pattern)))
    if not files:
        sys.exit(f"No files matched {pattern} in {d}")
    return [run_name(f) for f in files]


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
def parse_conditions_file(path, subject_column=None):
    """Parse an uploaded CSV/TSV: detect sample + group (+batch/covariate) columns
    however they're named. Returns (list of {sample, group, extras{}, subject},
    subject column name, subject confirmed?, columns{}). subject_column: the header the user
    confirmed as the animal / patient ("none" = detect nothing). columns: which header fills
    each covariate slot, and the extra headers that did NOT fit (conditions.csv has two
    covariate slots) -- reported, never dropped silently."""
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
        confirmed = bool(subject_column) and subject_column.lower() != "none"
        if confirmed:
            subj = next((h for h in headers if h == subject_column or norm(h) == norm(subject_column)), None)
            if subj is None:
                sys.exit(json.dumps({"error": f"--subject-column {subject_column!r}: no such column",
                                     "headers": headers}))
        elif subject_column:        # "none": the user said there is no subject column
            subj = None
        else:
            subj = next((h for h in headers if h not in (scol, gcol, bcol) and h and subject_header(h)), None)
        if scol is None or gcol is None:
            sys.exit(json.dumps({"error": "could not find sample and group columns",
                                 "headers": headers,
                                 "hint": "rename a column to 'sample'/'file' and 'group'/'condition', "
                                         "or have the agent pass --mapping-json instead"}))
        other = [h for h in headers if h not in (scol, gcol, bcol, subj) and h]
        columns = {"sources": {**({"Batch": bcol} if bcol else {}),
                               **{COV_COLS[i + 1]: h for i, h in enumerate(other[:2])}},
                   "not_written": other[2:]}
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
            out.append({"sample": s, "group": g, "extras": extras,
                        "subject": (row.get(subj) or "").strip() if subj else ""})
        return out, (column_name(subj) if subj else None), confirmed, columns


def intent_from_json(blob):
    """Accept {'mapping': {sample: group}} or {'groups': {group: [samples]}}, optionally with
    {'subjects': {sample: subject}, 'subject_column': 'Mouse'} when the user said which
    animal / patient each sample came from ("five IPs per mouse")."""
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
    subjects = d.get("subjects") or {}
    for item in out:
        item["subject"] = str(subjects.get(item["sample"], "")).strip()
    # subjects the agent took from the user's own words count as confirmed
    return (out, (column_name(d.get("subject_column") or "Subject") if subjects else None),
            bool(subjects), {"sources": {}, "not_written": []})


# ----------------------------------------------------------------------- map --
def do_map(out_path, runs, intent, subject_col=None, subject_confirmed=False, columns=None):
    run_norms = [norm(r) for r in runs]
    run_groups = defaultdict(set)      # run -> {groups}
    run_extras = {}                    # run -> extras
    run_subject = defaultdict(set)     # run -> {subjects}
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
            if item.get("subject"):
                run_subject[r].add(item["subject"])

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
    # The subject column, under its own name, last. A run matched to two subjects is a
    # conflict to confirm, like a run matched to two groups.
    subject_of = {r: next(iter(v)) for r, v in run_subject.items() if len(v) == 1}
    subject_conflicts = {r: sorted(v) for r, v in run_subject.items() if len(v) > 1}
    if not run_subject:
        subject_col = None
    # Do the VALUES look like the animal / patient each run came from? Unconfirmed and they
    # do not: keep the column as a covariate (the old behaviour) and ask, rather than
    # suggest a --block that would silently drop the factor from the model.
    subj_reasons = subject_assessment(subject_of, assigned) if subject_col else []
    ambiguous = bool(subject_col and subj_reasons and not subject_confirmed)
    columns = columns or {"sources": {}, "not_written": []}
    sources = dict(columns["sources"])
    subject_ambiguous = None
    if ambiguous:
        slot = next((c for c in COV_COLS[1:] if c not in cov_used), None)
        taken = ", ".join(f"{c} = {sources[c]}" for c in COV_COLS[1:] if c in sources)
        subject_ambiguous = {
            "column": subject_col, "reasons": subj_reasons, "written_as": slot,
            "written": bool(slot),
            "to_confirm": (f"if {subject_col} really is the animal / patient each run came from, "
                           f"re-run --map with --subject-column '{subject_col}' to keep it as the "
                           f"block; otherwise it stays " +
                           (f"a covariate ({slot})." if slot else
                            f"OUT of conditions.csv: it was NOT written, because both covariate "
                            f"slots are taken ({taken}) and conditions.csv has only two. To keep "
                            f"it as a covariate instead, remove one of those columns from the "
                            f"sample sheet and re-run --map."))}
        if slot:
            sources[slot] = subject_col
            for r, v in subject_of.items():
                run_extras.setdefault(r, {})[slot] = v
            cov_used = sorted(set(cov_used) | {slot}, key=COV_COLS.index)
        subject_col = None
    cols = ["File.Name", "Group"] + cov_used + ([subject_col] if subject_col else [])
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in runs:
            row = [r, assigned.get(r, "")]
            for c in cov_used:
                row.append(run_extras.get(r, {}).get(c, ""))
            if subject_col:
                row.append(subject_of.get(r, ""))
            w.writerow(row)

    sizes = Counter(assigned.values())
    singletons = [g for g, c in sizes.items() if c < 2]
    # --block needs every subject to hold >= 2 runs (run_de.R stops otherwise). All
    # singletons = the column is a sample id, not a blocking unit: no --block. The groups
    # relabelled can never be a block (run_de.R stops on it), confirmed or not.
    subj_sizes = Counter(subject_of.values()) if subject_col else Counter()
    subj_single = sorted(k for k, c in subj_sizes.items() if c < 2)
    relabelled = any("groups relabelled" in x for x in subj_reasons)
    block_ok = bool(subject_col and subj_sizes and not subj_single and not relabelled
                    and len(subject_of) == len(runs) and not subject_conflicts)
    subj_partial = bool(subject_col and subj_single and len(subj_single) < len(subj_sizes))
    subj_missing = [r for r in runs if subject_col and r not in subject_of and r not in subject_conflicts]
    needs_conf = bool(unassigned or conflicting or unmatched or singletons or ambiguous
                      or columns["not_written"]
                      or subject_conflicts or subj_partial or (subject_col and subj_missing))
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
            **({"subject_conflicting_runs": subject_conflicts,
                "runs_without_subject": subj_missing,
                "single_run_subjects": subj_single if subj_partial else []}
               if subject_col else {}),
            **({"subject_ambiguous": subject_ambiguous} if subject_ambiguous else {}),
            # extra sample-sheet columns beyond the two covariate slots: named, never dropped silently
            **({"columns_not_written": {
                "columns": columns["not_written"],
                "why": "conditions.csv carries Batch plus two covariates (Covariate1/2); these "
                       "extra columns did not fit and were NOT written",
                "fix": "if one matters to the model, remove a less important column from the "
                       "sample sheet (or pass the animal / patient column with "
                       "--subject-column so it is kept under its own name) and re-run --map"}}
               if columns["not_written"] else {}),
        },
        # which sample-sheet header each covariate slot holds
        "covariate_columns": sources,
        # The column naming the animal / subject each run came from, kept under its own
        # name. block_suggested: every run has one, every subject holds >= 2 runs, and the
        # values look like subjects (or the user confirmed them).
        "block_column": subject_col,
        "block_suggested": block_ok,
        "subject_confirmed": bool(subject_col and subject_confirmed),
        "subjects": dict(subj_sizes) if subject_col else {},
        "needs_confirmation": needs_conf,
        "guidance": "Confirm every item under 'ambiguities' with the user, then finalize "
                    "the CSV and run --validate. Do NOT proceed to a search while runs are "
                    "unassigned or conflicting."
                    + (f" Runs sharing a {subject_col} come from one source: run DE with "
                       f"--block {subject_col} so they are not treated as independent."
                       if block_ok else ""),
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
    ap.add_argument("--subject-column", help="--map --from-file: the header the user CONFIRMED names the "
                    "animal / patient each run came from (kept under its own name for --block); "
                    "'none' = there is no such column")
    ap.add_argument("--validate")
    ap.add_argument("--against", help="report to validate File.Name against")
    a = ap.parse_args()

    if a.list_runs:
        runs = get_runs(a)
        print(json.dumps({"n_runs": len(runs), "runs": runs}, indent=2))
    elif a.map_out:
        runs = get_runs(a)
        if a.from_file:
            intent, subject_col, confirmed, columns = parse_conditions_file(a.from_file, a.subject_column)
        elif a.mapping_json:
            intent, subject_col, confirmed, columns = intent_from_json(a.mapping_json)
        else: sys.exit("--map needs --from-file (uploaded file) or --mapping-json (agent intent)")
        do_map(a.map_out, runs, intent, subject_col, confirmed, columns)
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
