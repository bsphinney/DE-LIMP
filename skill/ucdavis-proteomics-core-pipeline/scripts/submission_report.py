#!/usr/bin/env python3
"""
submission_report.py  --  The CoreOmics submission a session answers: stored once in the
session, shown in every report.

A Core analysis answers a request. A lab filled in a CoreOmics form (PROT_0756) saying what
the samples are, who prepared them, which organism, and what they want done. That record is
the only thing that ties a report back to what was asked for, so the report shows it by
default, and the Methods, the analysis brief and the run log read the same record.

ONE RECORD (DE-LIMP rule 3). `attach` writes it into the session once:
    <session>/input/submission.json   the allowlisted record (schema submission_record/1)
    <session>/input/samples.tsv       the sample sheet as a table, for people to read; code
                                      reads the json
    <session>/session.json            {"coreomics": {internal_id, id, url, source, record}}
make_analysis_html.py, make_methods.py, analysis_prompt.py and session.py read it through
load(); record_run.py finds the PROT number in session.json.

ALLOWLIST, NOT BLOCKLIST. A CoreOmics record also carries emails, phone numbers, payment and
PPMS order fields, contacts and internal notes. Reports go to collaborators and get
forwarded, so sanitize() copies only the fields named below. A field CoreOmics adds later is
dropped until someone adds it here on purpose. Free text is also scrubbed of anything shaped
like an email address or a phone number, because submitters type them into descriptions.

WHO PREPARED THE SAMPLES decides the Sample preparation Methods. When the lab sent peptides,
every step before LC-MS/MS was theirs, and the Methods must say so rather than carry
placeholders the Core can never fill. prepared_by() is the one reading of the form.

SOURCE. A record read from CoreOmics says `source: CoreOmics`. With no CoreOmics token the
agent asks the key facts, and `attach --given` records them as `source: given by the user`.
They are never presented as the CoreOmics record.

    python3 submission_report.py attach --session <S> --record ~/core/PROT_0756
    python3 submission_report.py attach --session <S> --given '{"internal_id": "PROT_0756",
        "organism": "mouse", "prot_or_pep": "peptides", "sample_prep": "lab"}'
    python3 submission_report.py show  --session <S> [--format md|html|json]
    python3 submission_report.py notes --session <S>

Each prints one JSON object (`show --format md|html` prints the section itself). Exit codes:
0 ok, 2 refused or nothing to show (the JSON says why), 1 a bug.
"""
import argparse
import csv
import html
import json
import os
import re
import sys
from collections import Counter, OrderedDict, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import core_submission as cs  # noqa: E402  ids, dates and sample-id matching live there
from session import paths_for, read_raw_list  # noqa: E402  the session layout, one place

SCHEMA = "submission_record/1"
SOURCE_COREOMICS = "CoreOmics"
SOURCE_USER = "given by the user"

# The CoreOmics form fields that may leave CoreOmics: (key in the record, key in the form's
# submission_data). Verified against PROT_0756's form schema (2026-09-24).
FORM_FIELDS = (
    ("organism", "organism"),               # "Organism and Organism Expressed In"
    ("uniprot", "uniprot"),                 # free text: a link, an accession, or blank
    ("description", "description"),         # "Please describe what you want us to do"
    ("prot_or_pep", "prot_or_pep"),         # Intact Proteins | peptides | other/discuss
    ("sample_prep", "sample_prep"),         # who does extraction / digestion (enum)
    ("buffer", "buffer"),                   # "buffers/solvents"
    ("beads", "magbead_yes_no"),            # "beads"
    ("normalisation", "volume_or_mass"),    # "Normalization"
    ("data_analysis", "data_analysis"),     # the Core analyses, or raw data only
    ("instrument_wanted", "mass_spec_wanted"),
)
TEXT_KEYS = tuple(k for k, _ in FORM_FIELDS)
SAMPLE_KEYS = ("unique_id", "sample_name", "condition_name")

# sample_prep's two answers on the form (PROT_0756's schema enum). A record given by the user
# says "lab" or "core" instead.
LAB_PREP = "i have prepped my samples"
CORE_PREP = "i want the proteomics core to prepare"

EMAIL = re.compile(r"[\w.+-]+@[\w-]+(?:\.[\w-]+)+")
# A number introduced as one: "tel: 555-0100", "phone 752 1234".
LABELLED_PHONE = re.compile(r"\b(?:tel|phone|ph|cell|mobile|fax)\b\.?\s*[:#]?\s*\+?[\d\s().-]{6,}\d", re.I)
# International: "+", then 8-15 digits with the usual separators ("+44 20 7946 0958").
INTL_PHONE = re.compile(r"(?<![\w+])\+\d[\d\s().-]{6,22}\d(?!\d)")
# North American, with or without separators or a leading 1: "(530)555-0100", "5305550100",
# "530-5550100", "1 530 555 0100". Not a quantity: "100-200-3000 ug" is a dilution series.
PHONE = re.compile(r"(?<![\w+])(?:1[\s.-]?)?\(?\d{3}\)?[\s.-]?\d{3}[\s.-]?\d{4}(?!\d)"
                   r"(?!\s?(?:[unpmµ]?[gLlM]\b|mM\b|%))")
# The only links shown: a CoreOmics submission page, and a UniProt page.
COREOMICS_URL = re.compile(r"^https://[a-z0-9.-]+/submissions/[0-9a-f]{12}$")
UNIPROT_URL = re.compile(r"^https://(?:www\.|rest\.)?uniprot\.org/[\w./?=&%:-]*$")
MISSING, BLANK = "not in the record", "left blank on the form"
FROM_SUMMARY = "not in submission_summary.json (attach fetch's submission.json)"


class RecordError(Exception):
    """A submission record exists but cannot be read -- shown, never silently skipped."""


# ------------------------------------------------------------------------- the record --
def scrub(s):
    """Text with anything shaped like an email address or a phone number removed."""
    s = EMAIL.sub("[email removed]", str(s))
    s = LABELLED_PHONE.sub("[phone removed]", s)
    s = INTL_PHONE.sub(lambda m: "[phone removed]" if 8 <= sum(c.isdigit() for c in m.group(0)) <= 15
                       else m.group(0), s)
    return PHONE.sub("[phone removed]", s)


def clean(v):
    """Free text as it may appear in a report: contact details removed, whitespace tidied.
    None stays None (the record does not have the field); "" means it was left blank."""
    if v is None:
        return None
    if isinstance(v, (dict, list, tuple)):
        return None
    s = scrub(v)
    lines = [" ".join(ln.split()) for ln in s.replace("\r\n", "\n").split("\n")]
    return "\n".join(ln for ln in lines if ln).strip()


def _name(first, last):
    return " ".join(x for x in (clean(first), clean(last)) if x) or None


def _internal_id(v):
    try:
        kind, key = cs.normalize_submission(v)
    except ValueError:
        return None
    return key if kind == "internal_id" else None


def _hex_id(v):
    s = cs._s(v).lower()
    return s if cs.HEX_ID.match(s) else None


def _url(v):
    """A CoreOmics submission page, or None: nothing else -- no query string, no other site."""
    s = cs._s(v)
    return s if COREOMICS_URL.match(s) else None


def _types(v):
    vals = v if isinstance(v, list) else ([v] if cs._s(v) else [])
    return [t for t in (clean(x) for x in vals) if t]


def _samples(rows):
    out = []
    for r in rows if isinstance(rows, list) else []:
        if isinstance(r, dict):
            out.append({k: clean(r.get(k)) or "" for k in SAMPLE_KEYS})
    return out


def _base(source):
    return OrderedDict([("schema", SCHEMA), ("source", source), ("fetched_at", None),
                        ("internal_id", None), ("id", None), ("url", None),
                        ("submitted_date", None), ("pi", {"name": None, "department": None,
                                                          "institution": None}),
                        ("submitter", {"name": None}), ("basis", None)]
                       + [(k, None) for k in TEXT_KEYS]
                       + [("experiment_types", []), ("samples", [])])


def _from_coreomics(rec):
    sd = rec["submission_data"]
    pi = rec.get("pi") if isinstance(rec.get("pi"), dict) else {}
    dept = pi.get("department")
    dept = dept.get("name") if isinstance(dept, dict) else dept
    _campus, institution, _basis = cs.classify_campus(rec)
    day = cs.parse_date(rec.get("submitted"))
    out = _base(SOURCE_COREOMICS)
    out.update(basis="record", internal_id=_internal_id(rec.get("internal_id")), id=_hex_id(rec.get("id")),
               url=_url(rec.get("url")), submitted_date=day.isoformat() if day else None,
               pi={"name": _name(rec.get("pi_first_name") or pi.get("first_name"),
                                 rec.get("pi_last_name") or pi.get("last_name")),
                   "department": clean(dept) or None, "institution": clean(institution) or None},
               submitter={"name": _name(rec.get("first_name"), rec.get("last_name"))},
               experiment_types=_types(sd.get("proteomics_type")),
               samples=_samples(sd.get("samples")))
    for key, form_key in FORM_FIELDS:
        # Unanswered (null) is blank; absent from the form, or not a plain answer, is missing.
        v = sd.get(form_key)
        out[key] = (None if form_key not in sd or isinstance(v, (dict, list)) else
                    "" if v is None else clean(v))
    return out


def _from_summary(s):
    """core_submission.py fetch's submission_summary.json. It carries no buffer, beads,
    UniProt, sample type or normalisation -- attach the fetch folder's submission.json."""
    pi = s.get("pi") if isinstance(s.get("pi"), dict) else {}
    sub = s.get("submitter") if isinstance(s.get("submitter"), dict) else {}
    org = s.get("organism_as_submitted")
    out = _base(SOURCE_COREOMICS)
    out.update(basis="summary", fetched_at=clean(s.get("fetched_at")),
               internal_id=_internal_id(s.get("internal_id")),
               id=_hex_id(s.get("id")), url=_url(s.get("url")),
               submitted_date=clean(s.get("submitted_date")),
               pi={"name": clean(pi.get("name")) or None,
                   "department": clean(pi.get("department")) or None,
                   "institution": clean(pi.get("institution")) or None},
               submitter={"name": clean(sub.get("name")) or None},
               organism=clean(org.get("value") if isinstance(org, dict) else org),
               description=clean(s.get("description")), sample_prep=clean(s.get("sample_prep")),
               data_analysis=clean(s.get("data_analysis")),
               instrument_wanted=clean(s.get("instrument_wanted")),
               experiment_types=_types(s.get("experiment_types")),
               samples=_samples(s.get("samples")))
    return out


def _from_record(r):
    """Our own schema (re-projected, so sanitize is idempotent), or facts given by the user."""
    source = r.get("source") if r.get("schema") == SCHEMA else None
    out = _base(source if source in (SOURCE_COREOMICS, SOURCE_USER) else SOURCE_USER)
    pi = r.get("pi") if isinstance(r.get("pi"), dict) else {"name": r.get("pi")}
    sub = r.get("submitter") if isinstance(r.get("submitter"), dict) else {"name": r.get("submitter")}
    out.update(basis=r.get("basis") if r.get("basis") in ("record", "summary") and source else None,
               fetched_at=clean(r.get("fetched_at")), internal_id=_internal_id(r.get("internal_id")),
               id=_hex_id(r.get("id")), url=_url(r.get("url")),
               submitted_date=clean(r.get("submitted_date")),
               pi={k: clean(pi.get(k)) or None for k in ("name", "department", "institution")},
               submitter={"name": clean(sub.get("name")) or None},
               experiment_types=_types(r.get("experiment_types")),
               samples=_samples(r.get("samples")))
    for key in TEXT_KEYS:
        out[key] = clean(r.get(key))
    return out


def sanitize(obj):
    """Any submission record -> the allowlisted record. Accepts the raw CoreOmics record,
    core_submission.py's summary, this module's own record, or facts given by the user."""
    if not isinstance(obj, dict):
        raise RecordError("a submission record must be a JSON object")
    if isinstance(obj.get("submission_data"), dict):
        return _from_coreomics(obj)
    if cs._s(obj.get("schema")).startswith("core_submission"):
        return _from_summary(obj)
    return _from_record(obj)


def label(rec):
    return (rec or {}).get("internal_id") or (rec or {}).get("id") or "submission"


def prepared_by(rec):
    """("lab" | "core" | None, why). The form's sample_prep answer decides; the proteins /
    peptides answer is used when sample_prep says nothing. When the two disagree the answer
    is None -- the form contradicts itself and nobody should guess which half is right."""
    sp = cs._s(rec.get("sample_prep")).casefold()
    pp = cs._s(rec.get("prot_or_pep")).casefold()
    by_prep = ("lab" if sp.startswith(LAB_PREP) or sp == "lab" else
               "core" if sp.startswith(CORE_PREP) or sp == "core" else None)
    # Peptides can only have been made by the lab; intact proteins say nothing about who
    # extracted them, so they only ever count against a "the lab made peptides" answer.
    by_type = "lab" if pp == "peptides" else ("proteins" if pp == "intact proteins" else None)
    if by_prep and by_type and (by_prep, by_type) in (("lab", "proteins"), ("core", "lab")):
        return None, (f"the form contradicts itself: sample prep says "
                      f"\u201c{rec.get('sample_prep')}\u201d but it was sent as "
                      f"\u201c{rec.get('prot_or_pep')}\u201d")
    if by_prep:
        return by_prep, "sample prep answer"
    if by_type == "lab":
        return "lab", "sent as peptides"
    return None, "the form does not say who prepared the samples"


def sent_as_peptides(rec):
    """Only when the proteins/peptides answer says so: nothing else on the form does."""
    return cs._s(rec.get("prot_or_pep")).casefold() == "peptides"


def raw_data_only(rec):
    return cs.RAW_ONLY_MARKER in cs._s(rec.get("data_analysis")).lower()


# ------------------------------------------------------------------- session storage --
def _read_json(path):
    """The JSON at path, or None when there is no file. Errors name the file, not its full path:
    they can end up in a report that goes to a collaborator."""
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except FileNotFoundError:
        return None
    except OSError as e:
        raise RecordError(f"cannot read {os.path.basename(path)}: {e.strerror or type(e).__name__}")
    except ValueError as e:
        raise RecordError(f"{os.path.basename(path)} is not valid JSON ({e})")


def _write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, path)


def read_source(path):
    """A record file, or a folder holding one (fetch's: submission.json is preferred over
    submission_summary.json because only the full record has buffer, beads and UniProt)."""
    path = os.path.abspath(os.path.expanduser(path))
    if os.path.isdir(path):
        for name in ("submission.json", "submission_summary.json"):
            if os.path.isfile(os.path.join(path, name)):
                path = os.path.join(path, name)
                break
        else:
            raise RecordError(f"no submission.json or submission_summary.json in the folder "
                              f"{os.path.basename(path)}")
    obj = _read_json(path)
    if obj is None:
        raise RecordError(f"not found: {os.path.basename(path)}")
    return sanitize(obj), path


def _session_block(p):
    sj = _read_json(p["session_json"]) or {}
    if not isinstance(sj, dict):
        raise RecordError("session.json is not a JSON object")
    return sj, (sj.get("coreomics") if isinstance(sj.get("coreomics"), dict) else None)


def is_session(path):
    p = paths_for(path)
    return os.path.isfile(p["session_json"]) or os.path.isfile(p["submission_record"])


def load(session):
    """The session's submission record, or None when it has none. A record that session.json
    announces but that cannot be read raises RecordError: the report must say so, not drop it."""
    if not session:
        return None
    p = paths_for(os.path.expanduser(session))
    _sj, block = _session_block(p)
    obj = _read_json(p["submission_record"])
    if obj is None:
        if block:
            raise RecordError(f"session.json names a CoreOmics submission, but "
                              f"{os.path.relpath(p['submission_record'], p['session_dir'])} is missing")
        return None
    return sanitize(obj)


def resolve(source=None, session=None):
    """(record or None, session dir or None) from a --submission argument (a session dir, a
    record file, or fetch's folder) and/or a session dir."""
    if source:
        s = os.path.abspath(os.path.expanduser(source))
        if os.path.isdir(s) and is_session(s):
            return load(s), s
        return read_source(s)[0], session
    session = os.path.abspath(os.path.expanduser(session)) if session else None
    return load(session), session


def attach(session, rec, replace=False):
    """Write the record into the session (see the module docstring). Refuses to replace a
    different submission's record unless replace=True. The record is sanitized here too, so no
    caller can write anything but the allowlisted fields into a session."""
    rec = sanitize(rec)
    p = paths_for(os.path.expanduser(session))
    if not os.path.isdir(p["input_dir"]):
        raise RecordError(f"not a session directory (no input/): {p['session_dir']}")
    sj, old = _session_block(p)
    if old and not replace:
        same = (old.get("internal_id") and old.get("internal_id") == rec.get("internal_id")) or \
               (old.get("id") and old.get("id") == rec.get("id"))
        if not same:
            raise RecordError(f"this session is already attached to "
                              f"{old.get('internal_id') or old.get('id')}; pass --replace only if "
                              f"that was wrong")
    _write_json(p["submission_record"], rec)
    with open(p["submission_samples"], "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(SAMPLE_KEYS)
        for r in rec["samples"]:
            w.writerow([re.sub(r"[\t\r\n]+", " ", r[k]) for k in SAMPLE_KEYS])
    rel = lambda k: os.path.relpath(p[k], p["session_dir"]).replace(os.sep, "/")  # noqa: E731
    sj["coreomics"] = OrderedDict([("internal_id", rec.get("internal_id")), ("id", rec.get("id")),
                                   ("url", rec.get("url")), ("source", rec.get("source")),
                                   ("record", rel("submission_record")),
                                   ("samples", rel("submission_samples"))])
    _write_json(p["session_json"], sj)
    return {k: p[k] for k in ("session_json", "submission_record", "submission_samples")}


# --------------------------------------------------------------------- the design --
# A dash with a space on at least one side separates fields ("Old - IgG- Mouse 1"); a dash
# inside a word ("Wild-type", "Kv2.1") does not.
FIELD_SPLIT = re.compile(r"\s+-\s*|\s*-\s+")
# A field naming an individual: the blocking unit. "Rep 2" is NOT one -- a replicate index
# says nothing about which animal a sample came from.
UNIT = re.compile(r"^(?P<noun>mouse|mice|rat|animal|pig|dog|monkey|fish|patient|donor|subject|"
                  r"individual|participant|plant|tree)\s*#?\s*(?P<n>\d+)$", re.I)
REPLICATE_FIELD = re.compile(r"^(?:(?:bio(?:logical)?|tech(?:nical)?)\s*)?(?:rep(?:licate)?|r)?\s*#?\s*\d+$",
                             re.I)
REPLICATE_SUFFIX = re.compile(r"^(?P<base>.*?[A-Za-z)])[\s_-]*(?:rep(?:licate)?\s*)?#?\d+$", re.I)


def _key(s):
    return " ".join(cs._s(s).split()).casefold()


def _fields(name):
    return [f.strip() for f in FIELD_SPLIT.split(name)]


def replicate_base(name):
    """The condition without a trailing replicate index ("Old - A - Rep 2" -> "Old - A",
    "Control 3" -> "Control"), or None when it has none."""
    f = _fields(name)
    if len(f) > 1 and REPLICATE_FIELD.match(f[-1]):
        return " - ".join(f[:-1])
    m = REPLICATE_SUFFIX.match(cs._s(name))
    return m.group("base").strip() if m else None


def _ranges(nums):
    nums, out = sorted(set(nums)), []
    for n in nums:
        if out and n == out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return ", ".join(f"{a}–{b}" if b > a + 1 else (f"{a}, {b}" if b > a else str(a))
                     for a, b in out)


def _unit_design(rows, col, noun):
    """rows = [(uid, fields)], col = the unit field. A unit NUMBER alone names an individual
    only when the sheet numbers them across the whole study: when some other field keeps each
    number to one of its levels ("Mouse 1-3 Old, 4-6 Young"). Without that, "Mouse 1" under Old
    and under Young may be two animals whose numbering restarted, so a repeat is `unclear`,
    never read as pairing."""
    k = len(rows[0][1])
    nums = [int(UNIT.match(f[col]).group("n")) for _, f in rows]
    others = [c for c in range(k) if c != col]
    levels, by_level = {}, {}
    for c in others:
        levels[c] = list(OrderedDict((_key(f[c]), f[c]) for _, f in rows).values())
        by_level[c] = OrderedDict((lv, set()) for lv in levels[c])
        for (_uid, f), n in zip(rows, nums):
            by_level[c][next(lv for lv in levels[c] if _key(lv) == _key(f[c]))].add(n)
    multi = [c for c in others if len(levels[c]) > 1]
    between = [c for c in multi if sum(len(s) for s in by_level[c].values())
               == len(set().union(*by_level[c].values()))]
    # Who is who: with study-wide numbering, a number within its between-level; without it, every
    # combination of levels is its own individual (the cautious reading).
    ident = between if between else multi
    unit_of = [tuple(_key(f[c]) for c in ident) + (n,) for (_uid, f), n in zip(rows, nums)]
    factors = []
    for c in multi:
        complete = False
        if c in between:
            kind = "between"
        elif not between:
            kind = "unclear"
        else:
            spans = defaultdict(set)
            for (_uid, f), u in zip(rows, unit_of):
                spans[u].add(_key(f[c]))
            sizes = {len(v) for v in spans.values()}
            kind = ("within" if min(sizes) > 1 else "mixed")
            complete = sizes == {len(levels[c])}
        factors.append({"levels": levels[c], "kind": kind,
                        "complete": kind == "within" and complete,
                        "units_by_level": {lv: sorted(s) for lv, s in by_level[c].items()}})
    groups = OrderedDict()
    for uid, f in rows:
        g = " - ".join(x for i, x in enumerate(f) if i != col)
        groups.setdefault(_key(g), {"name": g, "samples": []})["samples"].append(uid)
    counts = defaultdict(int)
    for u in unit_of:
        counts[u] += 1
    return {"noun": noun, "field": col, "numbers": sorted(set(nums)),
            "by_sample": {uid: u for (uid, _f), u in zip(rows, unit_of)},
            "repeated": any(v > 1 for v in counts.values()),
            "factors": factors, "groups": list(groups.values())}


def design(samples):
    """What the sheet's condition names say about the design. Most names share one pattern of
    fields; when one field names an individual ("Mouse 3") it is the blocking unit and the rest
    are the conditions. Names that do not follow the pattern ("Pool") are listed, not dropped."""
    named = [(s["unique_id"], s["condition_name"]) for s in samples if s.get("condition_name")]
    out = {"n_samples": len(samples), "n_named": len(named),
           "n_distinct": len({_key(c) for _, c in named}), "unit": None, "outliers": []}
    if not named:
        return out
    split = [(uid, _fields(c)) for uid, c in named]
    k, n_k = Counter(len(f) for _, f in split).most_common(1)[0]
    if k < 2 or n_k * 2 <= len(split):
        return out
    rows = [(uid, f) for uid, f in split if len(f) == k]
    for col in range(k):
        ms = [UNIT.match(f[col]) for _, f in rows]
        nouns = {m.group("noun").lower().replace("mice", "mouse") for m in ms if m}
        if all(ms) and len(nouns) == 1:
            out["unit"] = _unit_design(rows, col, nouns.pop())
            out["outliers"] = [uid for uid, f in split if len(f) != k]
            return out
    return out


def sheet_groups(rec):
    """unique_id -> the group the sheet puts it in: the condition with the blocking unit set
    aside when the names carry one, else without a trailing replicate index when that leaves
    real groups, else the condition as written."""
    named = {s["unique_id"]: s["condition_name"] for s in rec["samples"] if s["condition_name"]}
    d = design(rec["samples"])
    if d["unit"]:
        out = dict(named)
        for g in d["unit"]["groups"]:
            out.update({uid: g["name"] for uid in g["samples"]})
        return out
    bases = {uid: replicate_base(c) or c for uid, c in named.items()}
    return bases if 1 < len({_key(b) for b in bases.values()}) < len(named) else named


# ------------------------------------------------------------ matching sheet to runs --
def run_owners(samples, names):
    """{name: {unique_id, ...}} by core_submission's matching: delimited tokens, the timsTOF
    sample field only, the longest id owns a run. Two ids on one run is an ambiguous run."""
    universe, cache = cs.id_universe(samples, [], "", None), {}
    return {n: {o["uid"] for o in cs.file_owners(
                cs.prepare_entry({"path": n, "name": os.path.basename(n.rstrip("/\\"))}),
                universe, cache)} for n in names}


def _read_conditions(p):
    try:
        with open(p["conditions"], newline="", encoding="utf-8-sig") as fh:
            rows = list(csv.DictReader(fh))
    except FileNotFoundError:
        return None
    return [r for r in rows if cs._s(r.get("File.Name"))]


def _de_design(p):
    """(columns the DE models besides Group, a block column or None, de_provenance or None).
    From de_provenance.json when the DE has run; before that, the columns run_de.R reads from
    conditions.csv (collect_conditions.COV_COLS). A `block` key in de_provenance names a
    column the DE blocked on (random effect / duplicateCorrelation)."""
    prov = _read_json(os.path.join(p["de_dir"], "de_provenance.json"))
    if isinstance(prov, dict) and cs._s(prov.get("design")).startswith("~"):
        terms = [x.strip() for x in prov["design"].split("~", 1)[1].split("+")]
        block = prov.get("block") if isinstance(prov.get("block"), str) else None
        return [x for x in terms if x not in ("0", "1", "groups", "")], block, prov
    import collect_conditions
    return list(collect_conditions.COV_COLS), None, None


# ------------------------------------------------------------------ quality notes --
def _examples(items, n=6):
    items = list(items)
    return ", ".join(items[:n]) + (f" and {len(items) - n} more" if len(items) > n else "")


def _binomial(name):
    return " ".join(cs._s(name).casefold().replace("(", " ").split()[:2])


def _organism_note(rec, p):
    org = rec.get("organism")
    if org is None:
        return None
    if not org:
        return ("organism_missing", "Organism: left blank on the form, so the submitter never "
                                    "stated which organism the samples are from.")
    fm = _read_json(p["fasta_meta"]) if p else None
    if not isinstance(fm, dict) or not (fm.get("organism") or fm.get("taxid")):
        return None
    f_org, f_tax = cs._s(fm.get("organism")), fm.get("taxid")
    try:
        import fetch_fasta
    except ImportError as e:
        return ("organism_unchecked", f"Organism: the form says “{org}”; it could not be "
                                      f"compared with the search database ({e}).")
    parts = [org] + [x for x in re.split(r"[,;/()]|\band\b|\bin\b", org) if x.strip()]
    taxa = {t for t in (fetch_fasta.alias_taxid(x) for x in parts) if t}
    names = {_binomial(fetch_fasta.ORGANISM_TAXIDS[t][0]) for t in taxa if t in fetch_fasta.ORGANISM_TAXIDS}
    if (f_tax and f_tax in taxa) or (f_org and (_binomial(f_org) in names
                                                or _binomial(f_org) in _binomial(org))):
        return None
    db = f"{f_org or '?'}" + (f" (taxid {f_tax})" if f_tax else "")
    if taxa:
        return ("organism_mismatch", f"Organism: the form says “{org}” but the search "
                                     f"database is {db}. Confirm which is right before using "
                                     f"these results.")
    return ("organism_unchecked", f"Organism: the form says “{org}”, which could not be "
                                  f"matched automatically to the search database ({db}). Confirm "
                                  f"they agree.")


def _levels_phrase(levels):
    return (f"both {levels[0]} and {levels[1]}" if len(levels) == 2
            else f"all {len(levels)} of {', '.join(levels)}")


def _model_sentence(rec, u, p):
    """What the DE does with the unit -- read from de_provenance.json / conditions.csv, never
    assumed. None when there is nothing to read yet."""
    rows = _read_conditions(p)
    cols, block, prov = _de_design(p)
    if not rows and not prov:
        return None
    noun = u["noun"]
    unit_of_run = {}
    for run, uids in run_owners(rec["samples"], [r["File.Name"] for r in rows or []]).items():
        if len(uids) == 1 and next(iter(uids)) in u["by_sample"]:
            unit_of_run[run] = u["by_sample"][next(iter(uids))]

    def carries(col):
        pairs = {(unit_of_run[r["File.Name"]], cs._s(r.get(col))) for r in rows or []
                 if r["File.Name"] in unit_of_run}
        return len(pairs) == len({x[0] for x in pairs}) == len({x[1] for x in pairs}) > 1

    if block and carries(block):
        return f"The analysis blocked on the {noun} (“{block}”)."
    fixed = [c for c in cols if carries(c)]
    if fixed:
        return f"The design analysed includes the {noun} as “{fixed[0]}”, a fixed effect."
    if prov:
        return (f"The design analysed ({prov['design']}) has no term for the {noun}, so it treated "
                f"the {prov.get('n_samples') or len(rows or [])} samples as independent.")
    return (f"No column the DE reads from input/conditions.csv (Group, {', '.join(cols)}) "
            f"identifies the {noun}, so as set up it will treat the samples as independent.")


def _pairing_note(rec, p):
    d = design(rec["samples"])
    u = d["unit"]
    if not u:
        return None
    noun, cap = u["noun"], u["noun"].capitalize()
    plural = {"mouse": "mice", "fish": "fish"}.get(noun, noun + "s")
    order = ("within", "unclear", "between", "mixed")
    facs = sorted(u["factors"], key=lambda f: order.index(f["kind"]))
    paired = any(f["kind"] == "within" for f in facs)
    if not (paired or u["repeated"] or any(f["kind"] == "unclear" for f in facs)):
        return None                       # every sample its own individual: nothing to pair
    bits = [f"Pairing in the sample sheet: the condition names carry the {noun} each sample came "
            f"from ({cap} {_ranges(u['numbers'])})."]
    for f in facs:
        if f["kind"] == "within":
            bits.append(f"{'Each' if f['complete'] else 'Some'} {noun} gave samples under "
                        f"{_levels_phrase(f['levels']) if f['complete'] else 'several of ' + ', '.join(f['levels'])}"
                        f", so comparisons among those are within-{noun} (paired).")
        elif f["kind"] == "unclear":
            lv = f["levels"]
            bits.append(f"The sheet uses the same {noun} numbers under {_levels_phrase(lv)}, so it "
                        f"does not say whether {cap} {u['numbers'][0]} under {lv[0]} and under "
                        f"{lv[1]} is the same {noun}: if so, comparisons among them are paired; "
                        f"if the numbering restarts per group, they are not. Ask the submitter.")
        elif f["kind"] == "between":
            by = "; ".join(f"{lv}: {cap} {_ranges(us)}" for lv, us in f["units_by_level"].items())
            bits.append(f"{' vs '.join(f['levels'])} differs between {plural} ({by}).")
    if not paired and u["repeated"] and not any(f["kind"] == "unclear" for f in facs):
        bits.append(f"Some {plural} have more than one sample under the same condition.")
    sizes = sorted({len(g["samples"]) for g in u["groups"]})
    bits.append(f"Setting the {noun} aside, the sheet defines {len(u['groups'])} groups of "
                f"{'/'.join(map(str, sizes))}.")
    if paired or u["repeated"]:
        bits.append(f"Samples from one {noun} are not independent.")
    if d["outliers"]:
        bits.append(f"{len(d['outliers'])} sample(s) do not follow this naming and are not "
                    f"part of this reading ({_examples(d['outliers'])}).")
    if p and (paired or u["repeated"]):
        m = _model_sentence(rec, u, p)
        if m:
            bits.append(m)
    return ("pairing", " ".join(bits))


def _sheet_vs_raw(rec, p):
    names = read_raw_list(p["session_dir"])
    if not names or not rec["samples"]:
        return []
    owners = run_owners(rec["samples"], names)
    found = {uid for o in owners.values() for uid in o}
    loose = [n for n, o in owners.items() if not o]
    uids = [s["unique_id"] for s in rec["samples"] if s["unique_id"]]
    none = [u for u in uids if u not in found]
    out = []
    if none:
        out.append(("sheet_ids_without_raw", f"Sample sheet vs raw files: {len(none)} of "
                                             f"{len(uids)} sample IDs on the sheet match no raw "
                                             f"file analysed here ({_examples(none)})."))
    if loose:
        out.append(("raw_without_sheet_id", f"Sample sheet vs raw files: {len(loose)} of "
                                            f"{len(names)} raw files carry no sample ID from the "
                                            f"sheet ({_examples(os.path.basename(n.rstrip('/')) for n in loose)})."))
    weak = [u for u in uids if cs.weak_reason(u)]
    if weak and (none or loose or len(weak) == len(uids)):
        out.append(("weak_sample_ids", f"{len(weak)} sample IDs are too short or numeric to "
                                       f"match file names reliably ({_examples(weak)})."))
    return out


def _conditions_vs_analysed(rec, p):
    rows = _read_conditions(p)
    d = design(rec["samples"])
    sheet = sheet_groups(rec)
    if not rows or not sheet or (not d["unit"] and len({_key(g) for g in sheet.values()}) == len(sheet)):
        return []
    owners = run_owners(rec["samples"], [r["File.Name"] for r in rows])
    groups_of, in_design = defaultdict(set), set()
    for r in rows:
        o = owners[r["File.Name"]]
        in_design |= o
        if len(o) == 1:                   # a run carrying two sheet ids says nothing here
            groups_of[next(iter(o))].add(cs._s(r.get("Group")))
    a2s, s2a = defaultdict(set), defaultdict(set)
    for uid, groups in groups_of.items():
        if uid in sheet:
            for g in groups:
                a2s[g].add(sheet[uid])
                s2a[sheet[uid]].add(g)
    bits = []
    for a, ss in sorted(a2s.items()):
        if len(ss) > 1:
            bits.append(f"group “{a}” pools the sheet's {_examples(sorted(ss), 4)}")
    for s, aa in sorted(s2a.items()):
        if len(aa) > 1:
            bits.append(f"the sheet's “{s}” is split across groups {_examples(sorted(aa), 4)}")
    left_out = [u for u in sheet if u not in in_design]
    if left_out and in_design:
        n = len(left_out)
        bits.append(f"{n} sheet sample{'s are' if n > 1 else ' is'} not in the analysed design "
                    f"({_examples(left_out)})")
    if not bits:
        return []
    return [("conditions_differ", "Conditions analysed vs the sample sheet: " + "; ".join(bits) + ".")]


def _no_groups(rec, d):
    """True when the sheet's condition names give no replicate groups at all: every name is
    its own, even with a trailing replicate index set aside ("sample 1..30" counts as none)."""
    if d["unit"] or d["n_named"] < 2 or d["n_distinct"] < d["n_named"]:
        return False
    names = [s["condition_name"] for s in rec["samples"] if s["condition_name"]]
    bases = {_key(replicate_base(c) or c) for c in names}
    return len(bases) in (1, len(names))


def quality_notes(rec, session=None):
    """[{"id", "text"}]: gaps in the submission, and places where it disagrees with what was
    analysed. Each is checked against the session's own files; a check with nothing to read
    (no FASTA record yet, no conditions.csv) says nothing rather than guessing."""
    p = paths_for(os.path.expanduser(session)) if session else None
    notes = []
    if rec["source"] == SOURCE_USER:
        notes.append(("source_user", "These submission details were given by the user during "
                                     "the analysis, not read from the CoreOmics record."))
    n = _organism_note(rec, p)
    if n:
        notes.append(n)
    if rec.get("uniprot") == "":
        notes.append(("uniprot_blank", "UniProt: left blank on the form, so the sequence database "
                                       "was chosen by the Core rather than specified by the "
                                       "submitter."))
    who, why = prepared_by(rec)
    if who is None and why.startswith("the form contradicts"):
        notes.append(("prep_unclear", f"Sample preparation: {why}."))
    if rec.get("beads") == "" and any(re.search(r"affinity|bead|enrich|bioid|turboid|immuno|pull", t, re.I)
                                      for t in rec["experiment_types"]):
        notes.append(("beads_blank", "Beads: left blank on the form although the experiment type is "
                                     f"{_examples(rec['experiment_types'], 3)}, so the bead or "
                                     f"enrichment step is recorded nowhere."))
    blank = [s["unique_id"] or "?" for s in rec["samples"] if not s["condition_name"]]
    if blank:
        notes.append(("sheet_conditions_blank", f"Sample sheet: {len(blank)} of {len(rec['samples'])} "
                                                f"samples have no condition ({_examples(blank)})."))
    d = design(rec["samples"])
    if _no_groups(rec, d):
        notes.append(("sheet_conditions_unique", "Sample sheet: every sample has its own condition "
                                                 "name, so the sheet gives no replicate groups. The "
                                                 "groups analysed were not taken from it and are "
                                                 "not corroborated by it."))
    if p:
        notes += _sheet_vs_raw(rec, p)
        notes += _conditions_vs_analysed(rec, p)
    pn = _pairing_note(rec, p)
    if pn:
        notes.append(pn)
    return [{"id": i, "text": t} for i, t in notes]


# ---------------------------------------------------------------------- rendering --
def _shown(v, rec=None):
    if v is None:
        return FROM_SUMMARY if (rec or {}).get("basis") == "summary" else MISSING
    return v or BLANK


def _prep_text(rec):
    who, _why = prepared_by(rec)
    if who == "lab":
        return "the lab prepared peptides" if sent_as_peptides(rec) else "the lab prepared the samples"
    return "the Core prepared the samples" if who == "core" else "not clear from the form"


def fields(rec):
    """(label, value) rows of the Submission section, in order. Values are plain text."""
    pi = rec["pi"]
    who = ", ".join(x for x in (pi.get("name"), pi.get("department"), pi.get("institution")) if x)
    return [
        ("PI", who or MISSING),
        ("Submitted by", _shown(rec["submitter"].get("name"))),
        ("Submitted", _shown(rec.get("submitted_date"))),
        ("Organism (as specified)", _shown(rec.get("organism"), rec)),
        ("UniProt", _shown(rec.get("uniprot"), rec)),
        ("Experiment type", "; ".join(rec["experiment_types"]) or MISSING),
        ("Sent as", _shown(rec.get("prot_or_pep"), rec)),
        ("Sample preparation", _prep_text(rec) + (f" (\u201c{rec['sample_prep']}\u201d)"
                                                  if rec.get("sample_prep") not in (None, "", "lab", "core") else "")),
        ("Buffer", _shown(rec.get("buffer"), rec)),
        ("Beads", _shown(rec.get("beads"), rec)),
        ("Normalisation", _shown(rec.get("normalisation"), rec)),
        ("Data analysis requested", _shown(rec.get("data_analysis"), rec)),
    ]


def _source_line(rec):
    if rec["source"] == SOURCE_USER:
        return (f"Submission {label(rec)}, as given by the user during the analysis \u2014 not "
                f"read from CoreOmics.")
    return (f"From CoreOmics submission {label(rec)}, as the submitter wrote it. Contact and "
            f"billing details are left out.")


def _h(v):
    return html.escape(cs._s(v), quote=True)


def _link(rec, text=None):
    t = _h(text or label(rec))
    return f"<a href='{_h(rec['url'])}'>{t}</a>" if rec.get("url") else t


def _uniprot_html(v, rec):
    if v is None or v == "":
        return f"<span class='blank'>{_h(_shown(v, rec))}</span>"
    return f"<a href='{_h(v)}'>{_h(v)}</a>" if UNIPROT_URL.match(v) else _h(v)


def render_html(rec, notes=()):
    """The Submission section as HTML (no heading: the caller adds it to its contents). Its look
    (.subm) is report_style.py's, with the rest of the report's."""
    rec = sanitize(rec)
    rows = [("Submission", _link(rec))]
    for k, v in fields(rec):
        if k == "UniProt":
            rows.append((k, _uniprot_html(rec.get("uniprot"), rec)))
        else:
            cls = " class='blank'" if v in (MISSING, BLANK, FROM_SUMMARY) else ""
            rows.append((k, f"<span{cls}>{_h(v)}</span>" if cls else _h(v)))
    out = ["<div class='subm'>", f"<p class='sub'>{_h(_source_line(rec))}</p>", "<dl>"]
    out += [f"<dt>{_h(k)}</dt><dd>{v}</dd>" for k, v in rows]
    out.append("</dl>")
    if rec.get("description"):
        out.append(f"<h3>Description, as written</h3><blockquote>{_h(rec['description'])}</blockquote>")
    if rec["samples"]:
        named = any(s["sample_name"] and s["sample_name"] != s["unique_id"] for s in rec["samples"])
        out.append(f"<h3>Sample sheet ({len(rec['samples'])} samples)</h3><ol class='sheet'>")
        for s in rec["samples"]:
            extra = f" <span class='blank'>({_h(s['sample_name'])})</span>" if named and s["sample_name"] else ""
            cond = _h(s["condition_name"]) if s["condition_name"] else "<span class='blank'>no condition</span>"
            out.append(f"<li><code>{_h(s['unique_id'] or '?')}</code>{cond}{extra}</li>")
        out.append("</ol>")
    if notes:
        out.append("<h3>Data quality notes</h3><ul>")
        out += [f"<li>{_h(n['text'])}</li>" for n in notes]
        out.append("</ul>")
    out.append("</div>")
    return "".join(out)


def _md(v):
    return cs._s(v).replace("|", "\\|").replace("\n", " ")


def render_markdown(rec, notes=(), heading="## Submission"):
    rec = sanitize(rec)
    ident = f"[{label(rec)}]({rec['url']})" if rec.get("url") else label(rec)
    L = [heading, "", f"*{_source_line(rec)}*", "", "| | |", "|---|---|",
         f"| Submission | {_md(ident)} |"]
    L += [f"| {k} | {_md(v)} |" for k, v in fields(rec)]
    if rec.get("description"):
        L += ["", "**Description, as written:**", ""] + [f"> {ln}" for ln in rec["description"].split("\n")]
    if rec["samples"]:
        L += ["", f"**Sample sheet ({len(rec['samples'])} samples):**", "", "| Sample ID | Condition |",
              "|---|---|"]
        L += [f"| {_md(s['unique_id'] or '?')} | {_md(s['condition_name'] or 'no condition')} |"
              for s in rec["samples"]]
    if notes:
        L += ["", "**Data quality notes:**", ""] + [f"- {n['text']}" for n in notes]
    return "\n".join(L) + "\n"


def one_line(rec):
    """For a README: the submission and where it came from."""
    ident = f"[{label(rec)}]({rec['url']})" if rec.get("url") else label(rec)
    return f"CoreOmics submission: {ident} (source: {rec['source']})"


def html_section(source=None, session=None):
    """The report's Submission section, or "" when the run has no submission. A record that
    exists but cannot be read is shown as such rather than left out."""
    try:
        rec, sess = resolve(source, session)
    except RecordError as e:
        import report_style
        return report_style.callout("warning", f"<p>The submission record could not be read: "
                                               f"{_h(e)}</p>")
    if rec is None:
        return ""
    return render_html(rec, quality_notes(rec, sess))


# --------------------------------------------------------------------------- CLI --
def _emit(obj, code=0):
    print(json.dumps(obj, indent=2))
    return code


def cmd_attach(a):
    if bool(a.record) == bool(a.given):
        return _emit({"error": "give exactly one of --record or --given"}, 2)
    warnings = []
    if a.record:
        rec, path = read_source(a.record)
        if os.path.basename(path) == "submission_summary.json":
            warnings.append("attached from submission_summary.json, which has no buffer, beads, "
                            "UniProt, sample type or normalisation: attach fetch's submission.json "
                            "to show them")
        if rec["source"] != SOURCE_COREOMICS:
            warnings.append(f"{path} is not a CoreOmics record; recorded as {rec['source']}")
    else:
        try:
            given = json.loads(a.given)
        except ValueError as e:
            return _emit({"error": f"--given is not JSON: {e}"}, 2)
        if not isinstance(given, dict):
            return _emit({"error": "--given must be a JSON object"}, 2)
        allowed = {"internal_id", "id", "url", "submitted_date", "pi", "submitter",
                   "experiment_types", "samples"} | set(TEXT_KEYS)
        ignored = sorted(set(given) - allowed)
        if ignored:
            warnings.append(f"ignored (not recorded): {', '.join(ignored)}")
        rec = sanitize({k: v for k, v in given.items() if k in allowed})
        rec["source"] = SOURCE_USER
    if not (rec.get("internal_id") or rec.get("id")):
        return _emit({"error": "the record has no PROT number or CoreOmics id"}, 2)
    session = os.path.abspath(os.path.expanduser(a.session))
    paths = attach(session, rec, a.replace)
    who, why = prepared_by(rec)
    return _emit({"attached": label(rec), "source": rec["source"], "prepared_by": who,
                  "prepared_by_basis": why, "n_samples": len(rec["samples"]),
                  "notes": quality_notes(rec, session),
                  "warnings": warnings, "outputs": paths})


def cmd_show(a):
    rec, sess = resolve(a.submission, a.session)
    if rec is None:
        return _emit({"error": "no submission is attached to this session",
                      "hint": "submission_report.py attach --session <S> --record <fetch dir>"}, 2)
    notes = quality_notes(rec, sess)
    if a.format == "html":
        print(render_html(rec, notes))
        return 0
    if a.format == "md":
        print(render_markdown(rec, notes), end="")
        return 0
    return _emit({"record": rec, "prepared_by": prepared_by(rec)[0], "notes": notes})


def cmd_notes(a):
    rec, sess = resolve(a.submission, a.session)
    if rec is None:
        return _emit({"error": "no submission is attached to this session"}, 2)
    return _emit({"submission": label(rec), "notes": quality_notes(rec, sess)})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    at = sub.add_parser("attach", help="store the submission record in a session")
    at.add_argument("--session", required=True)
    at.add_argument("--record", help="fetch's folder, its submission.json, or a summary")
    at.add_argument("--given", help="JSON of facts the user gave (no CoreOmics token)")
    at.add_argument("--replace", action="store_true",
                    help="replace a DIFFERENT submission already attached (only if it was wrong)")
    at.set_defaults(func=cmd_attach)
    for name, fn, hlp in (("show", cmd_show, "print the Submission section"),
                          ("notes", cmd_notes, "print the data-quality notes")):
        p = sub.add_parser(name, help=hlp)
        p.add_argument("--session")
        p.add_argument("--submission", help="a record file or fetch's folder, instead of the session's")
        if name == "show":
            p.add_argument("--format", choices=("json", "md", "html"), default="json")
        p.set_defaults(func=fn)
    a = ap.parse_args(argv)
    try:
        return a.func(a)
    except RecordError as e:
        return _emit({"error": str(e)}, 2)


if __name__ == "__main__":
    sys.exit(main())
