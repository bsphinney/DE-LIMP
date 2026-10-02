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

  # MAP an uploaded conditions file onto the real runs (CSV, TSV, or an Excel .xlsx such as
  # the Core LIMS export PROT_####.samples.<date>.xlsx -- read with the stdlib, no openpyxl)
  python3 collect_conditions.py --map conditions.csv --from-dir /data --glob '*.d' \
          --from-file user_conditions.csv [--sheet <name>] [--replicate-labels auto|collapse|keep]

  # MAP an agent-built mapping (from the user's free-text description) onto runs
  python3 collect_conditions.py --map conditions.csv --from-report report.parquet \
          --mapping-json '{"groups": {"control": ["A1","A2"], "treated": ["B1","B2"]}}'

  # emit a blank template (fallback when there's nothing to map)
  python3 collect_conditions.py --emit-template conditions.csv --from-dir /data --glob '*.d'

  # validate a finished design against the search output
  python3 collect_conditions.py --validate conditions.csv --against report.parquet
"""
import sys, os, csv, glob, json, re, argparse, io, zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

COV_COLS = ["Batch", "Covariate1", "Covariate2"]
# "unique id" / "sample id": the Core LIMS export (PROT_####.samples.<date>.xlsx: internal_id,
# internal_notes, sample_name, unique_id, condition_name, amt_2_inject)
SAMPLE_HEADERS = {"file.name", "filename", "file", "run", "sample", "sample name",
                  "samplename", "name", "raw", "raw file", "rawfile", "id", "unique id",
                  "sample id"}
GROUP_HEADERS = {"group", "condition", "treatment", "class", "type", "category",
                 "cohort", "phenotype", "condition name", "group name"}
# Core LIMS bookkeeping, never a covariate of the design: the submission's PROT number (the same
# on every row), staff notes, the amount injected. Named in columns_skipped, never dropped silently.
LIMS_BOOKKEEPING = {"internal id": "the CoreOmics submission number, the same for every sample",
                    "internal notes": "Core staff notes",
                    "amt 2 inject": "the amount injected -- an instrument setting, not a "
                                    "condition"}
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


# ------------------------------------------------------------------ encoding --
# conditions.csv is WRITTEN as UTF-8 always: open() without an encoding writes the platform's
# (cp1252 on Windows), and one accented sample name then cost the Methods and the SDRF
# downstream (deposit-builder, 2.8.0). What we READ may come from anywhere: UTF-8 (Excel's
# "CSV UTF-8" adds a BOM, dropped here) or Windows-1252 (Excel's plain "CSV" on Windows).
def _read_text(path):
    """A small text file's contents whatever wrote it: UTF-8 (BOM dropped), else
    Windows-1252, else Latin-1 (which cannot fail). Returns (text, encoding)."""
    with open(path, "rb") as fh:
        raw = fh.read()
    for enc in ("utf-8-sig", "cp1252"):
        try:
            return raw.decode(enc), enc
        except UnicodeDecodeError:
            pass
    return raw.decode("latin-1"), "latin-1"


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
    # a DIA-NN report can be gigabytes: streamed, UTF-8 (BOM dropped), never decoded whole
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as fh:
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


# A plate position in a run name: the autosampler slot and well, `_S3-A3_` (timsTOF / Evosep
# names). It says where a tube sat on the CORE's plate, never which sample it is: a submitter who
# labelled tubes A1..A4 matched the runs in wells A1..A4 by substring, and conditions followed the
# Core's wells with no question (review of 2.10, 2026-10-01).
PLATE_POSITION_RE = re.compile(r"(?<![A-Za-z0-9])S\d{1,2}-[A-P]\d{1,2}(?![A-Za-z0-9])", re.I)


def tokens(s):
    """Lower-case alphanumeric words: 'DIA-KG1_S3' -> ['dia', 'kg1', 's3']."""
    return re.findall(r"[a-z0-9]+", str(s).lower())


def run_tokens(run):
    """A run name's words, with the plate position's words (slot and well) masked as None so no
    identifier can match them."""
    low = str(run).lower()
    masked = [m.span() for m in PLATE_POSITION_RE.finditer(low)]
    return [None if any(a <= m.start() < b for a, b in masked) else m.group(0)
            for m in re.finditer(r"[a-z0-9]+", low)]


def _token_eq(want, have):
    """An identifier word against a run word: equal, or -- for a word of letters only -- the run
    word is it plus digits ('ctrl' names 'ctrl1'). A word with digits must match whole: 'a1' never
    names 'a10', 'kg1' never 'kg12'."""
    if have is None:
        return False
    return want == have or (want.isalpha() and have.startswith(want)
                            and have[len(want):].isdigit())


def _contains(seq, sub, eq):
    n = len(sub)
    return n > 0 and any(all(eq(sub[j], seq[i + j]) for j in range(n))
                         for i in range(len(seq) - n + 1))


def match_to_runs(identifier, runs):
    """Runs an identifier names: the whole name exactly (case and punctuation aside), else its
    words as a run of whole words of the run name (or the run name's inside the identifier -- a
    full file name names its run). Never a bare substring, never a plate position: 'A1' does not
    name '..._S3-A1_...', and 'KG1' does not name 'KG12'."""
    idn = norm(identifier)
    if not idn:
        return []
    exact = [r for r in runs if norm(r) == idn]
    if exact:
        return exact
    want = tokens(identifier)
    return [r for r in runs
            if _contains(run_tokens(r), want, _token_eq)
            or _contains(want, tokens(r), lambda a, b: a == b)]


# ---------------------------------------------------------------------- xlsx --
# A sample sheet saved from Excel, read with the standard library: an .xlsx is a zip of XML
# (ECMA-376). The pipeline env has no openpyxl (staff, 2026-09-25 and 09-28: the Core's LIMS
# export had to be parsed by hand). Members bigger than this are refused, not decompressed: a
# sample sheet is kilobytes.
XLSX_MEMBER_LIMIT = 50 * 1024 * 1024


def _local(tag):
    """An XML tag or attribute name without its namespace (transitional and strict OOXML alike)."""
    return tag.rsplit("}", 1)[-1]


def _xml(z, name):
    info = z.getinfo(name)
    if info.file_size > XLSX_MEMBER_LIMIT:
        raise ValueError(f"{name} is {info.file_size} bytes unpacked: not a sample sheet")
    return ET.fromstring(z.read(name))


def _col_index(ref):
    """'B12' -> 1 (0-based column)."""
    n = 0
    for ch in re.match(r"[A-Za-z]*", ref or "").group(0).upper():
        n = n * 26 + ord(ch) - 64
    return n - 1


def _cell_text(c, shared):
    t = c.get("t")
    if t == "inlineStr":
        return "".join(x.text or "" for x in c.iter() if _local(x.tag) == "t")
    v = next((x.text for x in c if _local(x.tag) == "v"), None)
    if v is None:
        return ""
    if t == "s":
        # a hand-edited or broken workbook: said as such, never an IndexError traceback
        try:
            return shared[int(v)]
        except (ValueError, IndexError):
            raise ValueError(f"cell {c.get('r') or '?'} points at shared string {v!r}, but the "
                             f"workbook holds {len(shared)}")
    if t == "b":
        return "TRUE" if v.strip() == "1" else "FALSE"
    if t in ("str", "e"):
        return v
    try:                                           # a number: 55 stays 55, never 55.0
        f = float(v)
        return str(int(f)) if f.is_integer() else v
    except ValueError:
        return v


def _cell_ref(ref):
    """'B12' -> (row 12, column 1); None for anything else."""
    m = re.fullmatch(r"([A-Za-z]{1,3})(\d+)", ref or "")
    return (int(m.group(2)), _col_index(m.group(1))) if m else None


def read_xlsx(path, sheet=None):
    """(headers, rows as dicts, info) of one worksheet of an .xlsx: `sheet` by name, else the
    first one that is not hidden. The header is the first row with anything in it; empty rows
    are skipped. A MERGED range reads as Excel shows it -- its top-left value in every cell of
    it -- not as one value and blanks (a group typed once over five merged rows gave four runs
    the group "", review of 2.10). info says which sheet was read, which others there are, and
    how many merged ranges were filled."""
    try:
        z = zipfile.ZipFile(path)
    except zipfile.BadZipFile:
        raise ValueError(f"{path} is not an .xlsx (not a zip archive)")
    with z:
        names = set(z.namelist())
        shared = []
        if "xl/sharedStrings.xml" in names:
            for si in _xml(z, "xl/sharedStrings.xml"):
                # rich text: every run's <t>, never the phonetic guide (<rPh>)
                shared.append("".join(t.text or "" for r in si.iter()
                                      if _local(r.tag) in ("si", "r") for t in r
                                      if _local(t.tag) == "t"))
        wb = _xml(z, "xl/workbook.xml")
        rels = {}
        if "xl/_rels/workbook.xml.rels" in names:
            for r in _xml(z, "xl/_rels/workbook.xml.rels"):
                tgt = r.get("Target") or ""
                rels[r.get("Id")] = tgt.lstrip("/") if tgt.startswith("/") else "xl/" + tgt
        sheets = []
        for sh in wb.iter():
            if _local(sh.tag) != "sheet":
                continue
            rid = next((v for k, v in sh.attrib.items() if _local(k) == "id" and "}" in k), None)
            sheets.append({"name": sh.get("name"), "path": rels.get(rid),
                           "hidden": sh.get("state") in ("hidden", "veryHidden")})
        if not sheets:
            raise ValueError(f"{path} has no worksheet")
        if sheet:
            pick = next((x for x in sheets if x["name"] == sheet), None)
            if pick is None:
                raise ValueError(f"no sheet {sheet!r} in {path} (sheets: "
                                 f"{', '.join(x['name'] for x in sheets)})")
        else:
            pick = next((x for x in sheets if not x["hidden"]), sheets[0])
        if not pick["path"] or pick["path"] not in names:
            raise ValueError(f"sheet {pick['name']!r} of {path} has no worksheet part")
        ws = _xml(z, pick["path"])
        cells, rnum = {}, 0            # (row, column) -> text, rows numbered as in the sheet
        for row in ws.iter():
            if _local(row.tag) != "row":
                continue
            rnum = int(row.get("r")) if (row.get("r") or "").isdigit() else rnum + 1
            nxt = 0
            for c in row:
                if _local(c.tag) != "c":
                    continue
                i = _col_index(c.get("r")) if c.get("r") else nxt
                cells[(rnum, i)] = _cell_text(c, shared).strip()
                nxt = i + 1
        merged = 0
        for m in ws.iter():
            if _local(m.tag) != "mergeCell":
                continue
            a, _, b = (m.get("ref") or "").partition(":")
            top, end = _cell_ref(a), _cell_ref(b or a)
            if not (top and end) or not cells.get(top):
                continue
            merged += 1
            for r in range(top[0], end[0] + 1):
                for col in range(top[1], end[1] + 1):
                    if not cells.get((r, col)):
                        cells[(r, col)] = cells[top]
        grid = []
        for r in sorted({k[0] for k in cells}):
            vals = {col: v for (row_, col), v in cells.items() if row_ == r}
            if any(vals.values()):
                grid.append([vals.get(i, "") for i in range(max(vals) + 1)])
    if not grid:
        raise ValueError(f"sheet {pick['name']!r} of {path} is empty")
    header = grid[0]
    headers = [h or f"column_{i + 1}" for i, h in enumerate(header)]
    rows = [{h: (r[i] if i < len(r) else "") for i, h in enumerate(headers)} for r in grid[1:]]
    return headers, rows, {"format": "xlsx", "sheet": pick["name"],
                           "other_sheets": [x["name"] for x in sheets if x is not pick],
                           **({"merged_ranges_filled": merged} if merged else {})}


def read_table(path, sheet=None):
    """(headers, rows as dicts, info) of a sample sheet: CSV / TSV (any encoding _read_text
    takes) or an Excel .xlsx / .xlsm. An old binary .xls is refused with what to do."""
    low = path.lower()
    if low.endswith((".xlsx", ".xlsm")):
        try:
            return read_xlsx(path, sheet)
        except (ValueError, KeyError, ET.ParseError, OSError) as e:
            sys.exit(json.dumps({"error": f"could not read {path} as an Excel workbook: {e}",
                                 "hint": "save it from Excel as .xlsx or CSV and try again"}))
    if low.endswith(".xls"):
        sys.exit(json.dumps({"error": f"{path} is an old binary Excel file (.xls), which this "
                                      "cannot read", "hint": "save it from Excel as .xlsx "
                                      "(Excel Workbook) or CSV, then pass that file"}))
    delim = "\t" if low.endswith((".tsv", ".txt")) else None
    text, enc = _read_text(path)
    with io.StringIO(text, newline="") as fh:
        sample = fh.read(4096); fh.seek(0)
        if delim is None:
            try: delim = csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
            except csv.Error: delim = ","
        rd = csv.DictReader(fh, delimiter=delim)
        headers = rd.fieldnames or []
        rows = [dict(r) for r in rd]
    return headers, rows, {"format": "tsv" if delim == "\t" else "csv", "encoding": enc}


# ---------------------------------------------------------- conditions source --
def parse_conditions_file(path, subject_column=None, runs=None, sheet=None, sample_column=None):
    """Parse an uploaded sample sheet -- CSV, TSV or Excel .xlsx (read_table) -- and detect the
    sample + group (+batch/covariate) columns however they're named. Returns (list of {sample,
    group, extras{}, subject}, subject column name, subject confirmed?, columns{}).
    subject_column: the header the user confirmed as the animal / patient ("none" = detect
    nothing). runs: the real run names -- with several sample-identifier columns, the one naming
    the most runs is used, and a TIE is asked (sample_column_tie), never settled by the sheet's
    order; sample_column: the header the user named. columns: which header fills each covariate
    slot, the extra headers
    that did NOT fit (conditions.csv has two covariate slots), and the ones skipped on purpose
    (LIMS bookkeeping, a second identifier column) with why -- reported, never dropped
    silently."""
    headers, table, source = read_table(path, sheet)
    hmap = {h: norm(h) for h in headers}
    lims = {norm(k): v for k, v in LIMS_BOOKKEEPING.items()}
    skipped = {h: lims[hmap[h]] for h in headers if hmap[h] in lims}
    snames = [h for h in headers if hmap[h] in {norm(x) for x in SAMPLE_HEADERS}]
    scol, scores, tie = (snames[0] if snames else None), {}, None
    if sample_column:
        scol = next((h for h in headers if h == sample_column
                     or norm(h) == norm(sample_column)), None)
        if scol is None:
            sys.exit(json.dumps({"error": f"--sample-column {sample_column!r}: no such column",
                                 "headers": headers}))
        snames = [scol] + [h for h in snames if h != scol]
    elif len(snames) > 1 and runs:
        # Several sample-identifier columns (the LIMS sheet has sample_name AND unique_id): the
        # one whose values name the most runs, grounded in the real run names. A TIE is the
        # user's to settle: the sheet's order once put a submitter's well-like names (A1..A4)
        # ahead of the Core's ids, and every run got another sample's condition.
        for h in snames:
            scores[h] = sum(1 for row in table if (row.get(h) or "").strip()
                            and match_to_runs(row[h], runs))
        best = max(scores.values())
        top = [h for h in snames if scores[h] == best]
        scol = top[0]
        if best and len(top) > 1:
            tie = {"columns": {h: scores[h] for h in top}, "used_for_now": scol,
                   "to_confirm": (f"ASK which column names the samples: {', '.join(top)} each "
                                  f"name {best} run(s). Re-run --map with --sample-column "
                                  "<header>; until then no run is given a group.")}
    for h in snames:
        if h != scol:
            skipped[h] = (f"another sample-identifier column; {scol} was used"
                          + (f" (it names {scores[scol]} run(s), this one {scores[h]})"
                             if scores else ""))
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
        subj = next((h for h in headers if h not in (scol, gcol, bcol) and h
                     and h not in skipped and subject_header(h)), None)
    if scol is None or gcol is None:
        others = source.get("other_sheets") or []
        sys.exit(json.dumps({"error": "could not find sample and group columns"
                                      + (f" on sheet {source['sheet']!r}" if source.get("sheet")
                                         else ""),
                             "headers": headers,
                             **({"sheet_read": source["sheet"], "other_sheets": others}
                                if source.get("sheet") else {}),
                             "hint": ((f"the workbook has other sheets ({', '.join(others)}): "
                                       "re-run with --sheet <name> for the one that holds the "
                                       "samples; or ") if others else "")
                                     + "rename a column to 'sample'/'file' and 'group'/"
                                       "'condition', or have the agent pass --mapping-json "
                                       "instead"}))
    other = [h for h in headers if h not in (scol, gcol, bcol, subj) and h and h not in skipped]
    columns = {"sources": {**({"Batch": bcol} if bcol else {}),
                           **{COV_COLS[i + 1]: h for i, h in enumerate(other[:2])}},
               "not_written": other[2:], "skipped": skipped,
               "sample_column": scol, "group_column": gcol,
               "sample_column_matches": scores or None, "sample_column_tie": tie,
               "source": source}
    out = []
    for row in table:
        s = (row.get(scol) or "").strip()
        # a tie gives no run a group: the CSV must not carry a column-order guess
        g = "" if tie else (row.get(gcol) or "").strip()
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


# Runs of one biological sample (core_submission.py locate --reinjections all keeps every
# injection) are technical replicates. Their sample is written under THIS column, whose name is
# the declaration: run_de.R blocks on it (random effect), so the injections are never counted as
# independent samples (2.10 review, MED 7). Mirrored by run_de.R's TECH_REPLICATE_COLUMN.
TECH_REPLICATE_COLUMN = "Sample"


def intent_from_json(blob):
    """Accept {'mapping': {sample: group}} or {'groups': {group: [samples]}}, optionally with
    {'subjects': {sample: subject}, 'subject_column': 'Mouse'} when the user said which
    animal / patient each sample came from ("five IPs per mouse"), and
    {'technical_replicates': {run name: biological sample}} (TECH_REPLICATE_COLUMN)."""
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
            bool(subjects), {"sources": {}, "not_written": [],
                             "technical_replicates": d.get("technical_replicates") or {}})


# ------------------------------------------------------------- replicate labels --
# A condition column that names each REPLICATE, not each condition: the Core LIMS export's
# condition_name is <sample>_mix_1 .. _5 (Core staff, 2026-09-25 and 09-28). Taken
# literally it is ten singleton groups and no DE at all; the condition is the label without its
# trailing replicate number. Read that way only when it is unambiguous -- otherwise ASKED.
_REP_TAIL = re.compile(r"^(?P<stem>.*?)(?P<n>\d+)$")
_REP_SEP = re.compile(r"[\s_.\-]+$")
_REP_WORD = re.compile(r"^(?P<stem>.+?)[\s_.\-]+(?:replicate|rep)$", re.I)


def split_replicate(label):
    """(condition, replicate number) for a label ending in a replicate number:
    'CtrlA_mix_3' -> ('CtrlA_mix', 3), 'Control1' -> ('Control', 1), 'WT_rep2' -> ('WT', 2).
    The number must follow a separator, or a letter. None when the label has no trailing number,
    or nothing would be left of it."""
    m = _REP_TAIL.match((label or "").strip())
    if not m:
        return None
    stem, n = m.group("stem"), int(m.group("n"))
    sep = _REP_SEP.search(stem)
    if sep:
        stem = stem[:sep.start()]
    elif not stem or not stem[-1].isalpha():
        return None
    w = _REP_WORD.match(stem)
    if w:
        stem = w.group("stem")
    return (stem, n) if stem else None


def collapse_replicates(intent, mode="auto"):
    """Read a per-replicate condition column as conditions + replicate numbers.

    `auto` collapses, in the PROPOSED csv, only when every sample has its own label (a
    per-replicate column) and the reading is unambiguous: every label ends in a replicate number,
    they make >= 2 conditions of >= 2 samples, no two claim the same replicate of a condition, no
    two conditions differ only in case or punctuation, and each condition's numbers run 1..n (or
    0..n-1) -- and even then it is the user's to confirm (replicate_labels_to_confirm): numbered
    labels can be real levels, and Day1..3 beside Ctl1..3 passes every one of those tests as a
    time course (review of 2.10). Anything doubtful is left as given and returned with the
    reasons (replicate_labels_ambiguous). `collapse` applies the proposal (the user said yes),
    `keep` uses the labels as they are. Returns (intent, record or None); a collapsed item keeps
    its label as `label` and its number as `replicate`."""
    labels = [it["group"] for it in intent if it.get("group")]
    counts = Counter(labels)
    per_sample = len(labels) >= 4 and len(counts) == len(labels)
    if mode == "keep" or not labels or (mode == "auto" and not per_sample):
        return intent, ({"detected": per_sample, "mode": mode, "collapsed": False}
                        if mode == "keep" and per_sample else None)
    parsed = {lab: split_replicate(lab) for lab in counts}
    groups, dup = defaultdict(dict), []
    for lab, v in sorted(parsed.items()):
        if v is None:
            continue
        if v[1] in groups[v[0]]:
            dup.append(f"{groups[v[0]][v[1]]} / {lab} (both {v[0]} replicate {v[1]})")
        else:
            groups[v[0]][v[1]] = lab
    unparsed = sorted(lab for lab, v in parsed.items() if v is None)
    reasons = []
    if unparsed:
        reasons.append(f"{len(unparsed)} label(s) end in no replicate number: "
                       + ", ".join(unparsed[:8]))
    if dup:
        reasons.append("two labels give the same condition and replicate: " + "; ".join(dup[:5]))
    if len(groups) < 2:
        reasons.append(f"they make {len(groups)} condition(s) ({', '.join(sorted(groups)) or '-'})"
                       ": nothing to compare -- the numbers may be the conditions themselves "
                       "(time points, doses, sample IDs)")
    small = sorted(g for g, m in groups.items() if len(m) < 2)
    if small:
        reasons.append("condition(s) with one sample: " + ", ".join(small)
                       + " -- the number may be part of the condition's name")
    variants = defaultdict(list)
    for g in groups:
        variants[norm(g)].append(g)
    clash = [sorted(v) for v in variants.values() if len(v) > 1]
    if clash:
        reasons.append("conditions that differ only in case or punctuation: "
                       + "; ".join(" / ".join(v) for v in clash))
    for g, m in sorted(groups.items()):
        nums = sorted(m)
        if len(nums) > 1 and (nums[0] not in (0, 1)
                              or nums != list(range(nums[0], nums[0] + len(nums)))):
            reasons.append(f"{g}: replicate numbers {', '.join(map(str, nums))} do not run "
                           f"{nums[0]}..{nums[0] + len(nums) - 1} -- a missing replicate, or not "
                           "replicate numbers")
    # forced (the user confirmed) or unambiguous -- and only if any label has a number at all
    collapse = (mode == "collapse" or not reasons) and bool(groups)
    rec = {"detected": True, "per_sample_labels": per_sample, "mode": mode,
           "collapsed": collapse,
           "conditions": {g: {str(n): lab for n, lab in sorted(m.items())}
                          for g, m in sorted(groups.items())},
           "not_numbered": unparsed, "reasons": reasons,
           "rule": "condition = the label without its trailing replicate number (and a "
                   "separator, or a 'rep' word, before it); conditions.csv keeps the label as "
                   "given in its Label column"}
    if collapse:
        if mode == "auto":
            rec["to_confirm"] = ("ASK the user: the proposed csv reads these labels as the "
                                 "conditions under `conditions` plus replicate numbers -- but "
                                 "the numbers may be the conditions themselves (time points, "
                                 "doses). Yes: re-run --map with --replicate-labels collapse. No: "
                                 "--replicate-labels keep, and the groups are the labels as given.")
        out = []
        for it in intent:
            v = parsed.get(it.get("group"))
            out.append(dict(it, group=v[0], label=it["group"], replicate=v[1]) if v else dict(it))
        return out, rec
    rec["to_confirm"] = ("ASK the user whether these labels are replicates of the conditions "
                         "proposed under `conditions` (reasons above). Yes: re-run --map with "
                         "--replicate-labels collapse. No: --replicate-labels keep, and the "
                         "groups are the labels as given.")
    return intent, rec


# ----------------------------------------------------------------------- map --
def decisions_path(out_path):
    """The record of the sample-identity answers the user gave (--confirm-multi), beside the CSV;
    provenance.py copies it into the bundle with conditions.csv."""
    return out_path + ".decisions.json"


def do_map(out_path, runs, intent, subject_col=None, subject_confirmed=False, columns=None,
           confirm_multi=(), replicate_mode="auto"):
    intent, replicates = collapse_replicates(intent, replicate_mode)
    run_label, run_rep = {}, {}        # run -> its per-replicate label / number, when collapsed
    run_groups = defaultdict(set)      # run -> {groups}
    run_extras = {}                    # run -> extras
    run_subject = defaultdict(set)     # run -> {subjects}
    multi_match = {}                   # identifier -> [runs] (one label, several files)
    unmatched = []                     # identifiers matching no run
    # A label that names several runs, or that two samples share, is a question about which run
    # IS which sample -- the user's to answer, never the agent's (a staff search, 2026-09-30: a
    # duplicate label was resolved by the agent with no question, before a multi-hour search,
    # and sample identity really was in doubt). Until the answer comes, those runs get NO group,
    # so --validate refuses the CSV: --confirm-multi '<label>' says the label names every one of
    # its runs (a group prefix); otherwise the mapping names the one run.
    confirmed = {norm(x) for x in confirm_multi if norm(x)}
    seen = Counter(norm(item["sample"]) for item in intent)
    duplicates = defaultdict(list)     # identifier -> [groups] (one label, several samples)
    pending = {}                       # identifier -> [runs], awaiting the user's answer
    confirmed_multi = {}               # identifier -> [runs], answered with --confirm-multi
    run_pending = defaultdict(set)     # run -> {identifiers awaiting an answer}

    for item in intent:
        hits = match_to_runs(item["sample"], runs)
        if not hits:
            unmatched.append(item["sample"]); continue
        if len(hits) > 1:
            multi_match[item["sample"]] = hits
        if seen[norm(item["sample"])] > 1:
            duplicates[item["sample"]].append(item["group"])
        if seen[norm(item["sample"])] > 1 or (len(hits) > 1 and norm(item["sample"]) not in confirmed):
            pending[item["sample"]] = hits
            for r in hits:
                run_pending[r].add(item["sample"])
            continue
        if len(hits) > 1:
            confirmed_multi[item["sample"]] = hits
        for r in hits:
            run_groups[r].add(item["group"])
            if item.get("label"):
                run_label.setdefault(r, item["label"])
                run_rep.setdefault(r, item.get("replicate"))
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
    awaiting = [r for r in runs if r in run_pending and r not in run_groups]
    unassigned = [r for r in runs if r not in run_groups and r not in run_pending]
    decisions = []
    for ident, hits in pending.items():
        if ident in duplicates:
            decisions.append({
                "identifier": ident, "kind": "duplicate_identifier", "runs": hits,
                "question": f"'{ident}' is given for {len(duplicates[ident])} samples (groups: "
                            f"{', '.join(map(str, duplicates[ident]))}), so it cannot say which run "
                            f"is which sample. Which run is each? Map each sample by something only "
                            f"its run carries."})
        else:
            decisions.append({
                "identifier": ident, "kind": "label_names_several_runs", "runs": hits,
                "question": f"'{ident}' matches {len(hits)} runs: {', '.join(hits)}. Does it name "
                            f"every one of them (a group label such as a replicate prefix), or ONE "
                            f"sample -- and then which run is it? Every one: re-run --map with "
                            f"--confirm-multi '{ident}'. One: name that run in the mapping."})
    unused = sorted(x for x in confirm_multi if norm(x) and norm(x) not in
                    {norm(k) for k in confirmed_multi})
    # Technical replicates: each run's biological sample (a run not named is its own sample).
    tech_map = (columns or {}).get("technical_replicates") or {}
    tech_of = {r: str(tech_map.get(r) or r) for r in runs} if tech_map else {}
    tech_runs = defaultdict(list)
    for r, smp in tech_of.items():
        tech_runs[smp].append(r)
    tech = ({"column": TECH_REPLICATE_COLUMN,
             "samples": {k: v for k, v in tech_runs.items() if len(v) > 1},
             "handling": "run_de.R blocks on the Sample column: the sample is a random blocking "
                         "factor, every contrast comes from that fit, and the injections are "
                         "never counted as independent samples"}
            if any(len(v) > 1 for v in tech_runs.values()) else None)
    dpath = decisions_path(out_path)
    record = {}
    if confirmed_multi:
        record.update(confirmed_multi_match=confirmed_multi,
                      how="--confirm-multi: the user said each label names every one of these runs "
                          "(collect_conditions.py --map)")
    if tech:
        record["technical_replicates"] = tech
    if record:
        with open(dpath, "w", encoding="utf-8") as fh:
            json.dump(dict(record, conditions_csv=os.path.abspath(out_path)), fh, indent=2)
    elif os.path.exists(dpath):
        os.remove(dpath)               # an older answer does not describe this mapping

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
    # Label: the per-replicate label as given, when the column was read as conditions +
    # replicates (make_figures.R names samples by it; unique per run, so never a block hint)
    cols = (["File.Name", "Group"] + cov_used + (["Label"] if run_label else [])
            + ([subject_col] if subject_col else []) + ([TECH_REPLICATE_COLUMN] if tech_of else []))
    if tech_of and cols.count(TECH_REPLICATE_COLUMN) > 1:
        sys.exit(json.dumps({"error": f"the sample sheet already has a column named "
                                      f"{TECH_REPLICATE_COLUMN!r}, which conditions.csv keeps for "
                                      f"technical replicates; rename that column"}))
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in runs:
            row = [r, assigned.get(r, "")]
            for c in cov_used:
                row.append(run_extras.get(r, {}).get(c, ""))
            if run_label:
                row.append(run_label.get(r, ""))
            if subject_col:
                row.append(subject_of.get(r, ""))
            if tech_of:
                row.append(tech_of[r])
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
    # numbered labels read as conditions + replicates: always the user's to confirm, and a
    # doubtful reading is not applied at all
    rep_unconfirmed = bool(replicates and replicates.get("collapsed")
                           and replicates.get("mode") == "auto")
    rep_ambiguous = bool(replicates and replicates.get("reasons")
                         and not replicates.get("collapsed"))
    tie = columns.get("sample_column_tie")
    needs_conf = bool(unassigned or conflicting or unmatched or singletons or ambiguous
                      or awaiting or pending or rep_unconfirmed or rep_ambiguous or tie
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
            # labels given for two or more samples: their runs await an answer too
            "duplicate_identifiers": dict(duplicates),
            # runs written with NO group until the user answers decisions_required
            "awaiting_decision_runs": awaiting,
            "singleton_groups": singletons,
            **({"subject_conflicting_runs": subject_conflicts,
                "runs_without_subject": subj_missing,
                "single_run_subjects": subj_single if subj_partial else []}
               if subject_col else {}),
            **({"subject_ambiguous": subject_ambiguous} if subject_ambiguous else {}),
            # a per-replicate condition column read as conditions + replicates: the user confirms
            # it; one that could not be read unambiguously is not applied: ASK either way
            **({"replicate_labels_to_confirm": replicates} if rep_unconfirmed else {}),
            **({"replicate_labels_ambiguous": replicates} if rep_ambiguous else {}),
            # two identifier columns name as many runs as each other: ASK which is the sample
            **({"sample_column_tie": tie} if tie else {}),
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
        # The sample-identity questions only the user can answer -- ask each, listing its runs;
        # never resolve one by another key (a plate well, a run number) without asking.
        "decisions_required": decisions,
        # answered: labels the user said name every one of their runs (--confirm-multi), also
        # recorded in <csv>.decisions.json
        "confirmed_multi_match": confirmed_multi,
        **({"decisions_file": os.path.abspath(dpath)} if record else {}),
        # several runs of one biological sample: their Sample column, which run_de.R blocks on
        **({"technical_replicates": tech} if tech else {}),
        **({"confirm_multi_unused": unused} if unused else {}),
        # which sample-sheet header each covariate slot holds
        "covariate_columns": sources,
        # what the sheet was and which of its columns were used / left out on purpose (LIMS
        # bookkeeping, a second identifier column), with why
        **({"source": columns["source"]} if columns.get("source") else {}),
        **({"sample_column": columns["sample_column"],
            "group_column": columns.get("group_column")} if columns.get("sample_column") else {}),
        **({"sample_column_matches": columns["sample_column_matches"]}
           if columns.get("sample_column_matches") else {}),
        **({"columns_skipped": columns["skipped"]} if columns.get("skipped") else {}),
        # a per-replicate condition column read as conditions + replicate numbers (or not)
        **({"replicate_labels": replicates} if replicates else {}),
        **({"replicates": run_rep} if run_rep else {}),       # run -> replicate number
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
                    + (" Ask every question under decisions_required, listing its runs -- "
                       "which run is which sample is the user's answer, never yours." if decisions
                       else "")
                    + (f" Runs sharing a {subject_col} come from one source: run DE with "
                       f"--block {subject_col} so they are not treated as independent."
                       if block_ok else ""),
    }, indent=2))


# ------------------------------------------------------------------ template --
def emit_template(out, names, covariates):
    cols = ["File.Name", "Group"] + covariates
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for n in names:
            w.writerow([n] + [""] * (len(cols) - 1))
    print(json.dumps({"template": out, "n_samples": len(names), "columns": cols, "samples": names}, indent=2))


# ------------------------------------------------------------------ validate --
def validate(meta_path, report_path):
    with io.StringIO(_read_text(meta_path)[0], newline="") as fh:
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
    ap.add_argument("--from-file", help="user-uploaded conditions file (CSV/TSV, or Excel .xlsx "
                    "such as the Core LIMS export; any column names)")
    ap.add_argument("--sheet", help="--from-file .xlsx: the worksheet to read (default: the first "
                    "one that is not hidden)")
    ap.add_argument("--replicate-labels", choices=("auto", "collapse", "keep"), default="auto",
                    help="a condition column naming each REPLICATE (X_1 .. X_5): auto reads it as "
                         "conditions + replicate numbers in the proposed csv when that is "
                         "unambiguous and ASKS the user to confirm it (replicate_labels_to_confirm), "
                         "and otherwise applies nothing and asks with the reasons "
                         "(replicate_labels_ambiguous); collapse applies the reading (the user "
                         "said yes); keep uses the labels as given")
    ap.add_argument("--mapping-json", help="agent-built intent: {'groups':{g:[samples]}} or {'mapping':{sample:g}}")
    ap.add_argument("--from-report")
    ap.add_argument("--from-dir")
    ap.add_argument("--glob", default="*.raw")
    ap.add_argument("--runs", help="comma-separated run names (alternative to --from-*)")
    ap.add_argument("--covariates", default="", help="comma list for --emit-template, e.g. Batch,Covariate1")
    ap.add_argument("--sample-column", help="--map --from-file: the header the user said names the "
                    "samples, when two identifier columns name as many runs (sample_column_tie)")
    ap.add_argument("--subject-column", help="--map --from-file: the header the user CONFIRMED names the "
                    "animal / patient each run came from (kept under its own name for --block); "
                    "'none' = there is no such column")
    ap.add_argument("--confirm-multi", action="append", default=[], metavar="LABEL",
                    help="--map: the USER said this label names every one of the runs it matches "
                         "(a group label such as a replicate prefix), not one sample. Repeatable; "
                         "recorded in <csv>.decisions.json")
    ap.add_argument("--validate")
    ap.add_argument("--against", help="report to validate File.Name against")
    a = ap.parse_args()

    if a.list_runs:
        runs = get_runs(a)
        print(json.dumps({"n_runs": len(runs), "runs": runs}, indent=2))
    elif a.map_out:
        runs = get_runs(a)
        if a.from_file:
            intent, subject_col, confirmed, columns = parse_conditions_file(
                a.from_file, a.subject_column, runs=runs, sheet=a.sheet,
                sample_column=a.sample_column)
        elif a.mapping_json:
            intent, subject_col, confirmed, columns = intent_from_json(a.mapping_json)
        else: sys.exit("--map needs --from-file (uploaded file) or --mapping-json (agent intent)")
        do_map(a.map_out, runs, intent, subject_col, confirmed, columns,
               confirm_multi=a.confirm_multi, replicate_mode=a.replicate_labels)
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
