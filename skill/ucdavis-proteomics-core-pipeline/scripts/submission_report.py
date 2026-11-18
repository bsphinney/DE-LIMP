#!/usr/bin/env python3
"""
submission_report.py  --  The CoreOmics submission a session answers: stored once in the
session, shown in every report.

A Core analysis answers a request. A lab filled in a CoreOmics form (PROT_0807) saying what
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

    python3 submission_report.py attach --session <S> --record ~/core/PROT_0807
    python3 submission_report.py attach --session <S> --given '{"internal_id": "PROT_0807",
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
import shlex
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
# submission_data). Verified against a real submission's form schema (2026-09-24).
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

# sample_prep's two answers on the form (a real submission's schema enum). A record given by the user
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
# session.json: Core staff marked a session with no CoreOmics submission as Core client work
# (session.py core-client), so core_run answers yes and its report carries the Core's survey
CORE_CLIENT = "core_client"

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


# The staff-only record of who marked a session as Core client work (session.py core-client): the
# login, the time, the reason. A mark in session.json counts only with its entry here.
CORE_CLIENT_STAFF_RECORD = os.path.join("logs", "core_client.staff.json")


def core_client_mark(sj, session_dir=None):
    """session.json's Core client work mark (session.py core-client), when it is one: an object
    with a `marked` time, whose staff-only record (CORE_CLIENT_STAFF_RECORD, in `session_dir`)
    ends with the matching entry -- event "core-client", by a login, at that time. session.py
    core-client writes both, only for Core staff (staff.core_staff); a mark typed into
    session.json by hand has no such entry and is ignored (survey_off_note says so). Anything
    else is None: fail closed."""
    b = sj.get(CORE_CLIENT) if isinstance(sj, dict) else None
    if not (isinstance(b, dict) and isinstance(b.get("marked"), str) and b["marked"].strip()
            and session_dir):
        return None
    import staff
    try:
        entries = staff.read(os.path.join(session_dir, CORE_CLIENT_STAFF_RECORD))
    except (OSError, ValueError):
        return None
    ours = [e for e in entries if cs._s(e.get("event")).startswith("core-client")]
    last = ours[-1] if ours else {}
    if (last.get("event") != "core-client" or last.get("at") != b["marked"]
            or not cs._s(last.get("by"))):
        return None
    return b


def core_run(session=None, record=None):
    """Is this a UC Davis Core run -- does it answer a CoreOmics submission, or did Core staff mark
    it as Core client work? THE one test for what only the Core's own clients see (its feedback
    survey), asked of the SESSION alone by the report and its .md twin, README, AGENTS.md and
    make_podcast.py (check, link, share, and the podcast's standing consent); a --submission given
    to make_analysis_html feeds only the Submission section. `record` (a delivery's summary,
    already read) is for the delivery README and its email. Yes when a record that names a
    submission (a PROT number or a CoreOmics id) loads -- `record`, or the one attached to
    `session` -- or when session.json parses and announces one, by that number or id, whose
    record file is missing (the report says so; still the Core's). Yes too, with no number, when
    session.json carries Core staff's mark that the session is Core client work with no
    submission (core_client_mark: with its staff record; session.py core-client, 2.15): a pilot
    or a rush job before its submission exists. Anything else is NO, fail closed: a damaged
    session.json or record, a block or record naming no submission, a mark without its staff
    record. The skill is used outside the Core too, and those users must never be sent the
    Core's survey. -> {"prot": "PROT_0807" | None}, {"prot": None, "core_client": True} when the
    yes rests on the mark alone (no CoreOmics submission: nothing to deliver through CoreOmics),
    or None."""
    def named(rec):
        return ({"prot": rec.get("internal_id")}
                if rec and (rec.get("internal_id") or rec.get("id")) else None)
    if record is not None:
        try:
            return named(sanitize(record))
        except RecordError:
            return None
    if not session:
        return None
    p = paths_for(os.path.expanduser(session))
    try:
        sj = _read_json(p["session_json"])
    except RecordError:
        return None                                  # damaged session.json
    if sj is not None and not isinstance(sj, dict):
        return None
    mark = ({"prot": None, "core_client": True} if core_client_mark(sj, p["session_dir"])
            else None)
    block = (sj or {}).get("coreomics")
    if not os.path.lexists(p["submission_record"]):
        prot = _internal_id(block.get("internal_id")) if isinstance(block, dict) else None
        if prot or (isinstance(block, dict) and _hex_id(block.get("id"))):
            return {"prot": prot}                    # announced by number or id, record missing
        return mark
    try:
        return named(load(session)) or mark
    except RecordError:
        return mark


def _core_work(p):
    """Do the session's own records say it is Core work? session.json with what init records for
    Core staff only (analysis_for: --submission / --internal, SKILL.md step 0a), a CoreOmics block
    or a Core client mark, or an attached record. Only then is the staff list asked
    (survey_off_note): anyone else's session costs no check, and no SSH from a laptop."""
    try:
        sj = _read_json(p["session_json"])
    except RecordError:
        sj = None
    if isinstance(sj, dict) and (isinstance(sj.get("analysis_for"), dict) or CORE_CLIENT in sj
                                 or isinstance(sj.get("coreomics"), dict)):
        return True
    return os.path.lexists(p["submission_record"])


def core_acquisition(session):
    """Did the UC Davis Proteomics Core acquire this session's data, by the session's own records?
    THE test for printing the Core's instrument acknowledgment (make_methods.py; 2.15.1 review
    finding 4). A reason (str) when yes, None when not recorded -- fail closed:
      * core_run(session): a CoreOmics submission attached (or announced with its record
        missing), or Core staff's client-work mark;
      * session.json "analysis_for" (session.py init, asked of Core staff only, SKILL.md 0a):
        {"submissions": [...]} -- the analysis answers a CoreOmics submission whose record is
        not attached yet -- or {"internal": "<title>"} -- the Core's own work.
    A damaged or missing session.json, or one that says none of these, is None: data analysed
    outside the Core never get the Core's sentence. A standalone Methods (no session) asks the
    user and passes make_methods.py --acquired-at-core (SKILL.md 9d)."""
    if not session:
        return None
    c = core_run(session)
    if c:
        return (f"CoreOmics submission {c['prot']}" if c.get("prot")
                else "Core staff's client-work mark (session.py core-client)")
    p = paths_for(os.path.expanduser(session))
    try:
        sj = _read_json(p["session_json"])
    except RecordError:
        return None
    af = sj.get("analysis_for") if isinstance(sj, dict) else None
    if not isinstance(af, dict):
        return None
    subs = [str(x) for x in (af.get("submissions") or []) if str(x).strip()]
    if subs:
        return (f"session.json analysis_for: submission {', '.join(subs)} (session.py init "
                f"--submission; its CoreOmics record is not attached)")
    title = af.get("internal")
    if isinstance(title, str) and title.strip():
        return f"session.json analysis_for: the Core's internal work ({title.strip()})"
    return None


def core_work(session):
    """Is this session Core work, by its own records or by core_run? The test that keeps a
    decision staff-only (normalization_check.own_analysis_signer): yes for a session init marked
    for Core staff (--submission / --internal), with a CoreOmics block, a Core client mark (backed
    or not) or an attached submission record -- and, failing closed, for a session.json that
    cannot be read OR does not exist (2.15.1 review minor 8: `--session <an empty folder>` opened
    the analyst's own signer). No only for a session whose session.json exists, reads, and says
    nothing of the Core."""
    p = paths_for(os.path.expanduser(session))
    try:
        sj = _read_json(p["session_json"])
    except RecordError:
        return True                                  # damaged: never loosen a staff-only rule
    if sj is None or not isinstance(sj, dict):
        return True                                  # missing / not an object: the same
    return bool(_core_work(p) or core_run(session))


def _mark_unbacked(sj, p):
    """A core_client block in session.json with no matching staff record (core_client_mark
    says no): typed by hand, or a laptop copy pulled without logs/core_client.staff.json."""
    return (isinstance(sj, dict) and CORE_CLIENT in sj
            and not core_client_mark(sj, p["session_dir"]))


def _survey_off_reason(p):
    """Why core_run says no, from the session's records: the clause survey_off_note opens with."""
    try:
        sj = _read_json(p["session_json"])
    except RecordError:
        return "its session.json cannot be read"
    if sj is not None and not isinstance(sj, dict):
        return "its session.json is not a JSON object"
    if _mark_unbacked(sj, p):
        return (f"session.json's Core client work mark is ignored (it has no matching entry in "
                f"{CORE_CLIENT_STAFF_RECORD}, and session.py core-client writes both)")
    rel = os.path.relpath(p["submission_record"], p["session_dir"])
    if os.path.lexists(p["submission_record"]):
        try:
            load(p["session_dir"])
        except RecordError:
            return f"its CoreOmics record ({rel}) cannot be read"
        return f"its CoreOmics record ({rel}) names no submission"
    from session import read_analysis_for
    named = (read_analysis_for(p["session_dir"])[0] or {}).get("submissions") or []
    if named:
        return (f"session.json names {', '.join(str(x) for x in named)} (session.py init "
                f"--submission), but no CoreOmics record is attached")
    return "this session answers no CoreOmics submission"


SURVEY_STEPS = ("report", "finalize", "status")


def survey_off_note(session, step="report"):
    """THE line a Core staff member sees when a session's report, README and AGENTS.md will NOT
    carry the Core's feedback survey (core_run says no), with the two ways to add it -- printed by
    make_analysis_html.py at the report step (`step` "report") and by session.py finalize, and
    listed by session.py status as a note (there it also says to run finalize again: finalize has
    written README, AGENTS.md and the zip without it). A warning, never a gate. None when core_run
    says yes, without a session, when the session's records say nothing of Core work (_core_work),
    or when this account is not Core staff (staff.core_staff, asked once): a user outside the Core
    never sees it."""
    if step not in SURVEY_STEPS:
        raise ValueError(f"survey_off_note step {step!r}: one of {', '.join(SURVEY_STEPS)}")
    if not session:
        return None
    p = paths_for(os.path.abspath(os.path.expanduser(session)))
    if core_run(p["session_dir"]) or not _core_work(p):
        return None
    if step != "report" and _report_carries_survey(p)[0]:
        return None                      # survey_missing_note says the report is out of date
    import staff
    if not staff.core_staff()[0]:
        return None
    q = shlex.quote(p["session_dir"])
    try:
        unbacked = _mark_unbacked(_read_json(p["session_json"]), p)
    except RecordError:
        unbacked = False
    lead = (f"If it was marked on HIVE, pull {CORE_CLIENT_STAFF_RECORD} with session.json and "
            "rebuild the report. Otherwise, for" if unbacked else "For")
    then = ("then rebuild the report." if step == "report" else
            "then rebuild the report (make_analysis_html.py; a report built on a laptop needs "
            f"session.json and {CORE_CLIENT_STAFF_RECORD} pulled again first) and run finalize "
            "again.")
    return (f"The Core's feedback survey will NOT be included: {_survey_off_reason(p)}, so its "
            f"report, README and AGENTS.md go out without it. {lead} work for a Core client, "
            f"attach the submission (python3 scripts/submission_report.py attach --session {q} "
            f"--record <its fetch folder>) or mark the session as Core client work (python3 "
            f"scripts/session.py core-client --session {q} --why-file - <<'EOF', then the "
            f"reason, then EOF), " + then)


def _report_carries_survey(p):
    """(carries, cannot): does output/Analysis_Report.html carry the Core's survey line
    (core_submission.FEEDBACK_URL)? carries None when there is no report; cannot = why it could
    not be read (carries None then too)."""
    page = os.path.join(p["output_dir"], "Analysis_Report.html")
    try:
        with open(page, encoding="utf-8", errors="replace") as fh:
            return cs.FEEDBACK_URL in fh.read(), None
    except FileNotFoundError:
        return None, None
    except OSError as e:
        return None, f"{type(e).__name__}: {e.strerror or e}"


def survey_missing_note(session):
    """THE line when the report and the session disagree about the Core's survey line, said by
    session.py finalize and status. Either way round:
      * the session IS a Core run (core_run says yes) but output/Analysis_Report.html does not
        carry the survey: built before the submission was attached or the session marked -- or on
        a laptop from a session.json pulled before (2.15 review H4c); refresh_feedback only
        rewords a line that is there, so nothing else would ever add it;
      * the session is NOT a Core run but the report carries it: built while the session was
        marked, before session.py core-client --clear (2.15 review2 L3) -- a survey sent to
        someone it is not for.
    None when there is no report yet or the two agree."""
    if not session:
        return None
    p = paths_for(os.path.abspath(os.path.expanduser(session)))
    core = core_run(p["session_dir"])
    carries, cannot = _report_carries_survey(p)
    pull = (f"a report built on a laptop needs session.json, input/submission.json and "
            f"{CORE_CLIENT_STAFF_RECORD} pulled again first")
    if not core:
        if not carries:
            return None
        return ("The report (output/Analysis_Report.html) carries the Core's feedback survey, "
                "but this session is no longer a Core run (it answers no CoreOmics submission and "
                "carries no Core client work mark -- the mark was cleared, or the record taken "
                f"away, after the report was built). Rebuild it (make_analysis_html.py; {pull}), "
                "then run finalize again.")
    if carries or (carries is None and not cannot):
        return None
    what = (f"CoreOmics submission {core['prot']}" if core.get("prot") else
            "Core client work (session.py core-client)" if core.get("core_client") else
            "a CoreOmics submission")
    head = (f"The report (output/Analysis_Report.html) cannot be read ({cannot}), so whether it "
            f"carries the Core's feedback survey is not known" if cannot else
            f"The report (output/Analysis_Report.html) does not carry the Core's feedback survey, "
            f"although this session is {what}: it was built before that was recorded")
    return head + f". Rebuild it (make_analysis_html.py; {pull}), then run finalize again."


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
# run_owners' "matched on" for a run whose name stage_map cannot place
UNKNOWN_SIDE = "unknown side"


def run_owners(samples, names, internal_id=None, smap=None, maybe=None):
    """({name: {unique_id, ...}}, {name: what it was matched on}) by core_submission's rule for
    one run (run_matches -- the rule locate and the HT manifest use): delimited tokens, the
    timsTOF sample field only, the longest id owns a run; a run of this submission's HT plate
    (core_submission.ht_plate_run) in its sample field, by unique_id or else by sample_name. Two
    ids on one run is an ambiguous run. `smap` (stage_map): a run Core staff renamed when staging
    is matched by its delivered (corrected) name -- staff's correction of which sample it is --
    whichever side of the rename the analysis used; a run whose name stage_map cannot place
    (UNKNOWN_SIDE) is matched to no sample, never by a guess. `maybe` (a dict): filled with each
    such run's possible samples (by either of its names), so the notes call them undetermined,
    not missing (2.15 review2 L4)."""
    universe = cs.id_universe(samples, [], "", None)
    plate = cs.ht_pattern(internal_id)
    snames = cs.name_universe(samples) if plate else []
    owners, on = {}, {}
    for n in names:
        r = _stage_run(n, smap)
        if r and r["side"] == "unknown":
            owners[n], on[n] = set(), UNKNOWN_SIDE
            if maybe is not None:
                maybe[n] = {o["uid"] for x in (r.get("either") or []) if x
                            for o in cs.run_matches(x, universe, snames, plate)[0]}
            continue
        found, how = cs.run_matches(r["delivered"] if r and r["side"] else n, universe, snames,
                                    plate)
        owners[n], on[n] = {o["uid"] for o in found}, how
    return owners, on


def _stage_run(name, smap):
    from collect_conditions import run_name
    return ((smap or {}).get("runs") or {}).get(run_name(cs._s(name)))


def _place(names, pairs):
    """{run name: {"side", "acquired", "delivered", "renamed"}} for each name analysed here, against
    the stage record's pairs [(name as acquired, delivered name, target path)]. The side is which
    of the pair the analysis read: "acquired" (the instrument's file -- its path is the pair's
    target, or its name is a name as acquired and no delivered one), "delivered" (a corrected
    link: its path leads to the target under the delivered name, or its name is only a delivered
    one), "unchanged" (both names are one), None (in no pair). A name that is one run's acquired
    name and ANOTHER's delivered name (a swap) takes the side the session's other runs used; when
    they say nothing it is "unknown": never matched by a guess (2.15 review H2)."""
    from collect_conditions import run_name

    def placed(side, pair):
        return {"side": side, "acquired": pair[0], "delivered": pair[1],
                "renamed": run_name(pair[0]) != run_name(pair[1])}
    acq, dlv, tgt = {}, {}, {}
    for a, d, t in pairs:
        acq.setdefault(run_name(a), (a, d))
        dlv.setdefault(run_name(d), (a, d))
        if t:                              # input/raw_names.tsv has no paths
            tgt.setdefault(os.path.normpath(t), (a, d))
    runs, undecided = {}, []
    for n in names:
        s, rn = cs._s(n).rstrip("/\\"), run_name(cs._s(n))
        if rn in runs or rn in undecided:
            continue
        base = os.path.basename(s)
        exact = tgt.get(os.path.normpath(s)) if os.path.isabs(s) else None
        if exact is None and os.path.isabs(s) and os.path.lexists(s):
            real = tgt.get(os.path.realpath(s))       # a link to the run: which name it carries
            if real and base == real[1] and base != real[0]:
                runs[rn] = placed("delivered", real)
                continue
            exact = real if real and base == real[0] else None
        a_hit, d_hit = acq.get(rn), dlv.get(rn)
        if exact is not None:
            runs[rn] = placed("acquired" if run_name(exact[0]) != run_name(exact[1]) else
                              "unchanged", exact)
        elif a_hit and d_hit and a_hit != d_hit:
            undecided.append(rn)
        elif a_hit:
            runs[rn] = placed("unchanged" if a_hit == d_hit else "acquired", a_hit)
        elif d_hit:
            runs[rn] = placed("delivered", d_hit)
        else:
            runs[rn] = {"side": None, "acquired": None, "delivered": None, "renamed": False}
    used = {r["side"] for r in runs.values() if r["side"] in ("acquired", "delivered")}
    for rn in undecided:
        if used == {"acquired"}:
            runs[rn] = placed("acquired", acq[rn])
        elif used == {"delivered"}:
            runs[rn] = placed("delivered", dlv[rn])
        else:                             # renamed either way: both of its pairs are renames
            # `either`: the name it is matched by under each reading (its pair as acquired, as
            # delivered) -- for the samples it may be, never for a match (review2 L4)
            runs[rn] = {"side": "unknown", "acquired": None, "delivered": None, "renamed": True,
                        "either": [acq[rn][1], dlv[rn][1]]}
    return runs


def stage_map(rec, p, names):
    """How the runs analysed here (`names`: input/raw_files.txt entries, or run names) stand in
    THIS submission's stage record: the service project's .core_submission.json, whose `links`
    name each raw/ link (a corrected name when staff fixed the instrument's) and its target (the
    run as acquired) -- core_submission.find_service_record, the one lookup (deliver's; without
    its slow search of the service tree), and raw_names. Without one above the session -- a copy
    pulled to a laptop for the report (step 9) -- the session's own input/raw_names.tsv
    (write_name_map), which travels with it. -> {"source": "stage record" | "input/raw_names.tsv"
    | None, "problem": None | "not found" | "unreadable" (said in the notes for a CoreOmics
    record, which the Core stages), "foreign": another submission's record above the session,
    "runs": _place, "renamed": [run names analysed here that staff renamed]}."""
    out = {"source": None, "problem": None, "foreign": None, "runs": {}, "renamed": []}
    if not p or not names or not rec.get("id"):
        return out
    stage, where, out["foreign"] = cs.find_service_record(
        {"id": rec["id"], "internal_id": rec.get("internal_id")}, p["session_dir"], None,
        search_tree=False)
    if stage is not None:
        pairs = [(os.path.basename(cs._s(t).rstrip("/")), d, cs._s(t))
                 for d, t in cs.raw_names(stage) if d and t]
        out["source"] = "stage record"
    else:
        try:
            with open(p["raw_names"], encoding="utf-8") as fh:
                pairs = [(a, d, "") for a, d, _pos in cs.parse_raw_name_map(fh.read())]
            out["source"] = "input/raw_names.tsv"
        except FileNotFoundError:
            pairs = None
        except (OSError, ValueError) as e:          # said, never taken for "no rename"
            pairs, out["map_problem"] = None, f"{type(e).__name__}: {e}"
        if pairs is None:
            if rec.get("source") == SOURCE_COREOMICS:
                out["problem"] = ("unreadable" if cs.STAGE_RECORD_UNREADABLE in cs._s(where)
                                  else "not found")
            return out
    out["runs"] = _place(names, pairs)
    out["renamed"] = [k for k, r in out["runs"].items() if r["renamed"]]
    return out


def _run_label(name, smap):
    """A raw file's name as the notes list it: its base name -- and, when staff renamed runs
    analysed here (stage_map), which side of the rename it is, with the other name."""
    base = os.path.basename(cs._s(name).rstrip("/\\"))
    if not (smap or {}).get("renamed"):
        return base
    r = _stage_run(base, smap)
    if r is None or r["side"] is None:
        return f"{base} (not in the stage record)"
    if r["side"] == "unknown":
        return (f"{base} (one run's name as acquired and another's delivered name: which this "
                "is, the records do not say)")
    if not r["renamed"]:
        return f"{base} (as acquired and as delivered)"
    if r["side"] == "delivered":
        return f"{base} (as delivered; acquired as {r['acquired']})"
    return f"{base} (as acquired; delivered as {r['delivered']})"


def session_locate(p):
    """Staff's locate decisions for the session: input/locate_decisions.json (attach writes it
    from locate.json), or that locate.json itself when it is newer -- locate re-run after
    attach. A staff report: attach at 12:15, the `--choose` locate at 12:18, and the session's
    copy still said no run was matched. None when the session has neither."""
    d = _read_json(p["locate_decisions"]) if os.path.isfile(p["locate_decisions"]) else None
    if not isinstance(d, dict):
        return None
    src, at = d.get("source"), d.get("source_mtime")
    if isinstance(src, str) and isinstance(at, (int, float)):
        try:
            newer = os.path.getmtime(src) > at
        except OSError:                   # not here (a laptop's copy of the session): the copy
            newer = False
        if newer:
            fresh = locate_decisions(src)
            # the fetch folder's locate.json may since have been run for another submission
            # (review of 2.13.2): only the same submission's replaces the session's copy
            rec = _read_json(p["submission_record"]) if os.path.isfile(
                p["submission_record"]) else None
            want = d.get("internal_id") or (rec or {}).get("internal_id")
            got = (fresh or {}).get("internal_id")
            if fresh is not None and not (got and want and got != want):
                return fresh
    return d


def run_assignment(rec, p, names, smap=None):
    """({name: {unique_id, ...}}, basis, how): which sheet sample each run is. A run is matched by
    collect_conditions.run_name, so raw paths and conditions.csv's File.Name both work. Staff's
    locate result settles every run it assigned (matched by itself or by `--choose`; a record
    from skill 2.13.1 or older lists only the staff choices): the one record of sample identity,
    and the only one when the names do not carry the ids (sheet S01, runs DIA-S1). Every other
    run is matched by name (run_owners: by unique_id, or an HT plate run by sample_name; a run
    staff renamed when staging by its delivered name): an empty or partial locate record never
    erases a name match (review of 2.13.2). basis: "locate" when locate assigned every run here,
    "mixed" when some, "names" when none. how: {"on": {name: what each run matched by name was
    matched on (core_submission.MATCHED_ON_*, or UNKNOWN_SIDE)}, "smap": stage_map}."""
    from collect_conditions import run_name
    loc = session_locate(p) if p else None
    uids = {s["unique_id"] for s in rec["samples"] if s["unique_id"]}

    def uid_of(key):                      # a sample key KG1#2 is the unique_id KG1's
        key = cs._s(key)
        return key if key in uids else re.sub(r"#\d+$", "", key)

    by_run = defaultdict(set)
    loc = loc if isinstance(loc, dict) else {}
    given = loc.get("assignments") if isinstance(loc.get("assignments"), list) else \
        loc.get("staff_choices") or []
    for a in given:
        if isinstance(a, dict) and a.get("sample") and a.get("file"):
            by_run[run_name(a["file"])].add(uid_of(a["sample"]))
    if smap is None:
        smap = stage_map(rec, p, names)
    maybe = {}
    by_name, on = run_owners(rec["samples"], [n for n in names if run_name(n) not in by_run],
                             rec.get("internal_id"), smap, maybe)
    owners = {n: set(by_run[run_name(n)]) if run_name(n) in by_run else by_name[n] for n in names}
    covered = sum(run_name(n) in by_run for n in names)
    return owners, ("locate" if names and covered == len(names) else
                    "mixed" if covered else "names"), {"on": on, "smap": smap, "maybe": maybe}


def _read_conditions(p):
    try:
        with open(p["conditions"], newline="", encoding="utf-8-sig") as fh:
            rows = list(csv.DictReader(fh))
    except FileNotFoundError:
        return None
    return [r for r in rows if cs._s(r.get("File.Name"))]


def _de_design(p):
    """(columns the design models besides Group, the column the DE blocked on or None,
    de_provenance or None). From de_provenance.json when the DE has run; before that, the
    columns run_de.R reads from conditions.csv without --block (collect_conditions.COV_COLS).
    The block is run_de.R's record: `block_column` exists exactly when the design is blocked
    (block.applied); `block` is an object, and {"applied": false} is an UNblocked run -- never
    read its truthiness."""
    prov = _read_json(os.path.join(p["de_dir"], "de_provenance.json"))
    if isinstance(prov, dict) and cs._s(prov.get("design")).startswith("~"):
        terms = [x.strip() for x in prov["design"].split("~", 1)[1].split("+")]
        b = prov.get("block")
        block = cs._s(prov.get("block_column")) or (
            cs._s(b.get("column")) if isinstance(b, dict) and b.get("applied") is True else "")
        return [x for x in terms if x not in ("0", "1", "groups", "")], block or None, prov
    import collect_conditions
    return list(collect_conditions.COV_COLS), None, None


def _plural(noun):
    return {"mouse": "mice", "fish": "fish"}.get(noun, noun + "s")


# ------------------------------------------------------------------ quality notes --
def _examples(items, n=6):
    items = list(items)
    return ", ".join(items[:n]) + (f" and {len(items) - n} more" if len(items) > n else "")


def _binomial(name):
    return " ".join(cs._s(name).casefold().replace("(", " ").split()[:2])


def same_species(a, b):
    """Two organism names for one species: the same first two words, an abbreviated genus
    ("S. enterica serovar Typhimurium LT2 (attenuated)") matching a spelled-out one
    ("Salmonella enterica ... LT2") -- a form and its database once read as different organisms
    for that."""
    x, y = _binomial(a).split(), _binomial(b).split()
    if len(x) < 2 or len(y) < 2 or x[1] != y[1]:
        return False
    abbr = lambda g: len(g) <= 2 and g.endswith(".")              # noqa: E731
    return x[0] == y[0] or (abbr(x[0]) and y[0].startswith(x[0][0])) \
        or (abbr(y[0]) and x[0].startswith(y[0][0]))


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
                                                or _binomial(f_org) in _binomial(org)
                                                or same_species(f_org, org))):
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
    for run, uids in run_assignment(rec, p, [r["File.Name"] for r in rows or []])[0].items():
        if len(uids) == 1 and next(iter(uids)) in u["by_sample"]:
            unit_of_run[run] = u["by_sample"][next(iter(uids))]

    def carries(col):
        pairs = {(unit_of_run[r["File.Name"]], cs._s(r.get(col))) for r in rows or []
                 if r["File.Name"] in unit_of_run}
        return len(pairs) == len({x[0] for x in pairs}) == len({x[1] for x in pairs}) > 1

    b = (prov or {}).get("block")
    b = b if isinstance(b, dict) else {}
    other = ""
    if block and carries(block):
        if b.get("effect") == "fixed":
            return (f"The analysis modelled the {noun} (“{block}”) as a fixed effect (--block), so "
                    f"every comparison is made within each {noun}.")
        rho = b.get("consensus_correlation")
        between = [c for c, k in (b.get("contrast_structure") or {}).items() if k == "between"]
        return (f"The analysis blocked on the {noun} (“{block}”, a random effect"
                + (f"; consensus within-{noun} correlation {rho:.2f}"
                   if isinstance(rho, (int, float)) else "") + ")."
                + (f" {len(between)} comparison(s) ({_examples(between, 3)}) set different "
                   f"{_plural(noun)} against each other: their evidence is the number of "
                   f"{_plural(noun)}, not of samples." if between else ""))
    if block:
        other = f"The analysis blocked on “{block}”, which does not identify the {noun}. "
    fixed = [c for c in cols if carries(c)]
    if fixed:
        return other + f"The design analysed includes the {noun} as “{fixed[0]}”, a fixed effect."
    if prov:
        if "block" not in prov:
            return (f"Whether the analysis modelled the {noun} is not recorded: its "
                    f"de_provenance.json predates run_de.R's --block, and the design "
                    f"({prov['design']}) has no term for the {noun}.")
        given = f" ({b['note']})" if b.get("column") and b.get("note") else ""
        return (other + f"The design analysed ({prov['design']}) has no term for the {noun}, so "
                f"it treated the {prov.get('n_samples') or len(rows or [])} samples as "
                f"independent{given}.")
    import collect_conditions
    subj = [h for h in (rows[0] if rows else {}) if collect_conditions.subject_header(h)
            and carries(h)]
    if subj:
        return (f"input/conditions.csv carries the {noun} as “{subj[0]}”: the DE models it only "
                f"when run with --block {subj[0]}; without it, it treats the samples as "
                f"independent.")
    return (f"No column the DE reads from input/conditions.csv (Group, {', '.join(cols)}) "
            f"identifies the {noun}, so as set up it will treat the samples as independent.")


def _pairing_note(rec, p):
    d = design(rec["samples"])
    u = d["unit"]
    if not u:
        return None
    noun, cap = u["noun"], u["noun"].capitalize()
    plural = _plural(noun)
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


def _sheet_vs_raw(rec, p, smap=None):
    """The sheet against the runs analysed here (input/raw_files.txt), each note with the basis
    it was matched on: "locate" / "mixed" (staff's assignment), "unique_id" (an ordinary
    submission: the run names carry the sample IDs), "sample_name" (an HT plate: its runs carry
    the sample names, never the sheet's IDs), "unique_id or sample_name", or "none" (an HT plate
    whose sheet gives no sample names: nothing can match by design, said in one neutral line).
    A run matched by name to more than one sample (a sample name several samples share, two IDs
    of one length on one run) is AMBIGUOUS, as locate's curated_rows has it: it counts for none of
    them, and no note claims a full match then (2.15 review H1). Raw names are listed as
    _run_label gives them; `names` says which names the runs were matched by ("delivered" when
    the stage record renamed runs analysed here, else "as analysed": the names in
    input/raw_files.txt, with no record that any was changed). A run stage_map cannot place is
    neither matched nor "loose": _stage_notes says so."""
    names = read_raw_list(p["session_dir"])
    if not names or not rec["samples"]:
        return []
    owners, basis, how = run_assignment(rec, p, names, smap)
    smap = how["smap"]
    tag = {"names": "delivered" if smap["renamed"] else "as analysed"}
    on = how["on"]                        # the runs matched by name, and on what
    HT = (cs.MATCHED_ON_HT_ID, cs.MATCHED_ON_HT_NAME)
    shared = [n for n in names if n in on and len(owners[n]) > 1]
    sole = {u for n in names if n not in on or len(owners[n]) == 1 for u in owners[n]}
    in_shared = {u for n in shared for u in owners[n]}
    loose = [n for n in names if not owners[n] and on.get(n) != UNKNOWN_SIDE]
    uids = [s["unique_id"] for s in rec["samples"] if s["unique_id"]]
    # samples whose only possible runs are renamed runs stage_map could not place: undetermined,
    # not missing (review2 L4; _stage_notes says which runs)
    undecided = {u for v in (how.get("maybe") or {}).values() for u in v} - sole - in_shared
    none = [u for u in uids if u not in sole and u not in in_shared and u not in undecided]
    sheet_name = {s["unique_id"]: s["sample_name"] for s in rec["samples"] if s["unique_id"]}
    such = f", such as {uids[0]}" if uids else ""          # never a made-up example (rule 2)
    out = []
    files = lambda ns: _examples(_run_label(n, smap) for n in ns)   # noqa: E731

    def note(nid, text, on_key):
        out.append((nid, text, dict(tag, basis=on_key)))

    def ambiguous(key, on_key):
        if not shared:
            return
        if key == "sample name":
            counts = Counter(sheet_name.get(u) or u for u in in_shared)
            who = ", ".join(f"{k} on {v}" for k, v in sorted(counts.items()))
        else:
            who = _examples(sorted(in_shared))
        note("sheet_ambiguous",
             f"Sample sheet vs raw files: {len(shared)} of {len(names)} raw files carry a {key} "
             f"that more than one sample on the sheet has ({len(in_shared)} samples: {who}), so "
             f"the names cannot tell which of those samples each of these runs is "
             f"({files(shared)}). The Core's assignment of runs to samples (core_submission.py "
             f"locate --choose) settles it.", on_key)

    if undecided:
        note("sheet_undetermined",
             f"Sample sheet vs raw files: {len(undecided)} of {len(uids)} samples on the sheet "
             f"cannot be told from the runs analysed here: their runs are ones Core staff renamed "
             f"when staging, and which side of the rename this analysis used could not be placed "
             f"(the raw file names note says which runs) ({_examples(sorted(undecided))}).",
             "stage record")
    if basis == "mixed":                  # some runs by staff's assignment, the rest by name
        key = ("sample ID or sample name" if any(on.get(n) in HT for n in names)
               else "sample ID")
        ambiguous(key, "mixed")
        if none:
            note("sheet_ids_without_raw",
                 f"Sample sheet vs raw files: {len(none)} of {len(uids)} sample IDs on the "
                 f"sheet have no raw file analysed here, by the Core's assignment of runs "
                 f"to samples or, for the runs it did not assign, by name "
                 f"({_examples(none)}).", "mixed")
        if loose:
            note("raw_without_sheet_id",
                 f"Sample sheet vs raw files: {len(loose)} of {len(names)} raw files "
                 f"analysed here were assigned to no sample and carry no {key} from the "
                 f"sheet ({files(loose)}).", "mixed")
        return out
    if basis == "locate":                 # staff assigned the runs: say so, not "match"
        if none:
            note("sheet_ids_without_raw",
                 f"Sample sheet vs raw files: {len(none)} of {len(uids)} sample IDs on the "
                 f"sheet have no raw file analysed here in the Core's assignment of runs to "
                 f"samples ({_examples(none)}).", "locate")
        if loose:
            note("raw_without_sheet_id",
                 f"Sample sheet vs raw files: {len(loose)} of {len(names)} raw files "
                 f"analysed here are assigned to no sample on the sheet ({files(loose)}).",
                 "locate")
        return out                        # how well ids match names no longer matters
    ht = any(on.get(n) in HT for n in names)
    by_id = any(owners[n] and on.get(n) in (cs.MATCHED_ON_ID, cs.MATCHED_ON_HT_ID)
                for n in names)
    by_name = any(owners[n] and on.get(n) == cs.MATCHED_ON_HT_NAME for n in names)
    has_names = any(s["sample_name"] for s in rec["samples"])
    if ht and not has_names and not by_id:
        note("sheet_not_compared",
             f"Sample sheet vs raw files: not compared by name. HT plate runs carry the plate "
             f"position and a sample name, never the sheet's sample IDs (unique_id{such}), and "
             f"the sheet gives no sample names.", "none")
        return out
    if ht and has_names and (by_name or not by_id):      # an HT plate: by sample name
        if by_name and by_id:
            key, on_key = "sample ID or sample name", "unique_id or sample_name"
            lead = ("matched by sample ID, and the HT plate runs that carry none by sample name "
                    "(sample_name)")
        elif by_name:
            key, on_key = "sample name", "sample_name"
            lead = (f"matched by sample name. HT plate runs carry the sheet's sample names "
                    f"(sample_name), not its sample IDs (unique_id{such})")
        else:
            key, on_key = "sample name", "sample_name"
            lead = (f"compared by sample name (HT plate runs carry sample names, not the sheet's "
                    f"sample IDs{such}), and none of the runs analysed here carries a sample "
                    f"name from the sheet")
        # no full match while the stage record is missing or unreadable: the names analysed
        # may be the instrument's wrong ones (2.15 review M2)
        full = (by_name and not none and not loose and not shared and not undecided
                and not smap["problem"])
        note("sheet_match_basis",
             "Sample sheet vs raw files: " + lead
             + ("; every sample on the sheet has a run analysed here, and every run a sample."
                if full else ".")
             + (" Runs Core staff renamed when staging were matched by the name they are "
                "delivered under." if tag["names"] == "delivered" else ""), on_key)
        ambiguous(key, on_key)
        if none:
            note("sheet_ids_without_raw",
                 f"Sample sheet vs raw files: {len(none)} of {len(uids)} samples on the sheet "
                 f"match no raw file analysed here by {key} ("
                 + _examples(f"{u} ({sheet_name[u]})" if sheet_name.get(u) else u for u in none)
                 + ").", on_key)
        if loose:
            note("raw_without_sheet_id",
                 f"Sample sheet vs raw files: {len(loose)} of {len(names)} raw files carry no "
                 f"{key} from the sheet ({files(loose)}).", on_key)
        return out
    ambiguous("sample ID", "unique_id")
    if none:
        note("sheet_ids_without_raw", f"Sample sheet vs raw files: {len(none)} of "
                                      f"{len(uids)} sample IDs on the sheet match no raw "
                                      f"file analysed here ({_examples(none)}).", "unique_id")
    if loose:
        note("raw_without_sheet_id", f"Sample sheet vs raw files: {len(loose)} of "
                                     f"{len(names)} raw files carry no sample ID from the "
                                     f"sheet ({files(loose)}).", "unique_id")
    weak = [u for u in uids if cs.weak_reason(u)]
    if weak and (none or loose or len(weak) == len(uids)):
        note("weak_sample_ids", f"{len(weak)} sample IDs are too short or numeric to "
                                f"match file names reliably ({_examples(weak)}).", "unique_id")
    return out


def write_name_map(session, service_project=None):
    """input/raw_names.tsv from this submission's stage record (core_submission.raw_name_map and
    format_raw_name_map: the one format), when Core staff renamed runs when staging: the pairs
    that travel with the session -- the step-9 pull, the zip -- and that deliver delivers as
    RAW_NAMES.tsv. Written by attach and finalize, where the stage record is (HIVE). -> (level,
    note): ("OK", what) written; ("SKIPPED", why) the stage record or the submission record
    cannot be read; (None, why) nothing to write -- no submission, no stage record above the
    session (a laptop copy keeps the map it was pulled with), or nothing renamed (an old map is
    then removed: it would describe renames the stage record no longer has). `service_project`
    (attach --service-project, as deliver's): the staged project of a session that is not inside
    it -- its map is then written here and the notes read it (2.15 review2 L7)."""
    p = paths_for(os.path.abspath(os.path.expanduser(session)))
    try:
        rec = load(p["session_dir"])
    except RecordError as e:
        return "SKIPPED", f"the submission record cannot be read ({e})"
    if not rec or not rec.get("id"):
        return None, "no CoreOmics submission attached"
    stage, where, _foreign = cs.find_service_record(
        {"id": rec["id"], "internal_id": rec.get("internal_id")}, p["session_dir"],
        os.path.abspath(os.path.expanduser(service_project)) if service_project else None,
        search_tree=False)
    if stage is None:
        if service_project:
            return "SKIPPED", f"--service-project: {where}"
        if cs.STAGE_RECORD_UNREADABLE in cs._s(where):
            return "SKIPPED", f"the stage record above the session ({cs.MARKER}) cannot be read"
        return None, "no stage record above the session"
    rows = cs.raw_name_map(stage)
    renamed = cs.raw_name_map_renamed(rows)
    if not renamed:
        if os.path.exists(p["raw_names"]):
            os.remove(p["raw_names"])
        return None, "the stage record renames no run"
    tmp = p["raw_names"] + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(cs.format_raw_name_map(rows))
    os.replace(tmp, p["raw_names"])
    return "OK", f"{len(rows)} staged runs, {len(renamed)} renamed by Core staff when staging"


def _stage_notes(rec, p, smap):
    """What the stage record says about the runs' names (stage_map): that it was not found above
    the session or cannot be read (a CoreOmics submission's -- the Core stages those -- so the
    names in the notes are only as analysed); or, when staff corrected names when staging, how
    many runs, which names the analysis used, and how the notes label them."""
    names = read_raw_list(p["session_dir"])
    if not names:
        return []
    if smap["problem"]:
        other = (f" The stage record above it is {smap['foreign']}'s." if smap["foreign"] else "")
        if smap.get("map_problem"):
            other += (f" input/raw_names.tsv cannot be read either ({smap['map_problem']}).")
        if smap["problem"] == "unreadable":
            return [("stage_record_unreadable",
                      f"Raw file names: the stage record of {label(rec)} above this session "
                      f"({cs.MARKER}) cannot be read, so whether Core staff corrected run names "
                      f"when staging is not known here; the names in these notes are as "
                      f"analysed. Repair it, or stage again (core_submission.py stage --apply).")]
        return [("stage_record_not_found",
                  f"Raw file names: no stage record of {label(rec)} ({cs.MARKER}, written by "
                  f"core_submission.py stage) was found in the folders above this session, so "
                  f"whether Core staff corrected run names when staging is not known here; the "
                  f"names in these notes are as analysed." + other)]
    if not smap["renamed"]:
        return []
    from collect_conditions import run_name
    runs = [smap["runs"].get(run_name(n)) or {} for n in names]
    sides = Counter(r.get("side") for r in runs if r.get("renamed") and r.get("side") != "unknown")
    if not sides:
        used = "names the stage record cannot place (below)"
    elif set(sides) == {"acquired"}:
        used = "the names as acquired"
    elif set(sides) == {"delivered"}:
        used = "the corrected (delivered) names"
    else:
        used = (f"the names as acquired for {sides['acquired']} of them and the delivered names "
                f"for {sides['delivered']}")
    unknown = sum(r.get("side") == "unknown" for r in runs)
    absent = sum(r.get("side") is None for r in runs)
    return [("raw_names_corrected",
             f"Raw file names: Core staff corrected the names of {len(smap['renamed'])} of the "
             f"{len(names)} runs analysed here when staging them, and the delivered raw/ folder "
             f"holds those runs under the corrected names. input/raw_files.txt, the search and its "
             f"tables use {used}. Raw-file names in these notes are labelled as acquired or as "
             f"delivered; join on the delivered name or the plate position, never on a name as "
             f"acquired. "
             + (f"input/raw_names.tsv (delivered as {cs.RAW_NAME_MAP}, beside the report) lists "
                f"every pair, with the plate position: join on it."
                if os.path.isfile(p["raw_names"]) else
                f"Every pair is in the stage record ({cs.MARKER}); finalize writes them to "
                f"input/raw_names.tsv, delivered as {cs.RAW_NAME_MAP}.")
             + (f" {unknown} run(s) carry a name that is one run's name as acquired and "
                f"another's delivered name, and nothing here says which: they are matched to no "
                f"sample." if unknown else "")
             + (f" {absent} run(s) are not in the stage record." if absent else ""))]


def _conditions_vs_analysed(rec, p):
    rows = _read_conditions(p)
    d = design(rec["samples"])
    sheet = sheet_groups(rec)
    if not rows or not sheet or (not d["unit"] and len({_key(g) for g in sheet.values()}) == len(sheet)):
        return []
    owners = run_assignment(rec, p, [r["File.Name"] for r in rows])[0]
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


def locate_decisions(locate_json):
    """locate.json -> what staff decided about the runs: the re-injection policy and its result
    per sample, each run chosen by --choose, and the whole assignment of runs to samples (every
    matched row: `assignments`, which run_assignment reads). Runs by file NAME only -- the notes
    a collaborator reads are made from it; `source` (and its time) is the path of the
    locate.json it came from, kept in the session's input/ (never delivered), so a locate re-run
    after attach for the same submission is read from there (session_locate). None when there
    is no locate.json."""
    loc = _read_json(locate_json)
    if not isinstance(loc, dict):
        return None
    base = lambda p: os.path.basename(cs._s(p).rstrip("/"))   # noqa: E731
    acc = loc.get("accepted") if isinstance(loc.get("accepted"), dict) else {}
    reinj, chosen, assigned = {}, [], []
    for r in loc.get("samples") or []:
        if not isinstance(r, dict):
            continue
        key = cs._s(r.get("sample_key")) or cs._s(r.get("unique_id"))
        rec = r.get("reinjections") if isinstance(r.get("reinjections"), dict) else None
        if rec:
            reinj.setdefault(key, {"policy": rec.get("policy"),
                                   "kept": [base(p) for p in rec.get("kept") or []],
                                   "set_aside": [base(p) for p in rec.get("set_aside") or []]})
        if cs._s(r.get("chosen_by")) and r.get("file"):
            chosen.append({"sample": key, "file": base(r["file"]), "by": cs._s(r.get("chosen_by"))})
        if r.get("status") == "matched" and r.get("file") and cs._s(r.get("unique_id")):
            assigned.append({"sample": cs._s(r.get("unique_id")), "file": base(r["file"])})
    return {"reinjections_policy": acc.get("reinjections"), "reinjections": reinj,
            "staff_choices": chosen, "accept_ambiguous": bool(acc.get("accept_ambiguous")),
            "allow_partial": bool(acc.get("allow_partial")), "assignments": assigned,
            "internal_id": loc.get("internal_id"),
            "source": os.path.abspath(locate_json), "source_mtime": os.path.getmtime(locate_json)}


def _locate_notes(p, smap=None):
    d = session_locate(p)
    if not isinstance(d, dict):
        return []
    out = []
    lab = lambda ns: ", ".join(_run_label(n, smap) for n in ns)   # noqa: E731
    reinj = d.get("reinjections") or {}
    if reinj:
        pol = d.get("reinjections_policy")
        each = "; ".join(f"{k}: kept {lab(v.get('kept') or [])}"
                         + (f", set aside {lab(v['set_aside'])}" if v.get("set_aside") else "")
                         for k, v in reinj.items())
        out.append(("reinjections",
                    f"Re-injections: {len(reinj)} sample(s) had several runs of their own, and staff "
                    + ("chose the newest injection of each (--reinjections latest)" if pol == "latest"
                       else "kept every injection as technical replicates (--reinjections all), "
                            "which are not independent samples" if pol == "all"
                       else f"decided --reinjections {pol}") + f" -- {each}."))
    if d.get("staff_choices"):
        out.append(("staff_run_choices",
                    f"Run choices: Core staff assigned the run by hand for "
                    f"{len(d['staff_choices'])} sample(s), where the file names alone did not "
                    f"settle which run is the sample ("
                    + "; ".join(f"{c['sample']}: {_run_label(c['file'], smap)}"
                                for c in d["staff_choices"]) + ")."))
    return out


def quality_notes(rec, session=None):
    """[{"id", "text"}]: gaps in the submission, and places where it disagrees with what was
    analysed. Each is checked against the session's own files; a check with nothing to read
    (no FASTA record yet, no conditions.csv) says nothing rather than guessing. The sheet-vs-raw
    notes also say what they matched on ("basis") and by which names ("names": "as analysed",
    or "delivered" when staff renamed runs when staging -- raw_names_corrected then says so)."""
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
        smap = stage_map(rec, p, read_raw_list(p["session_dir"]))
        notes += _stage_notes(rec, p, smap)
        notes += _sheet_vs_raw(rec, p, smap)
        notes += _conditions_vs_analysed(rec, p)
        notes += _locate_notes(p, smap)
    pn = _pairing_note(rec, p)
    if pn:
        notes.append(pn)
    # a note may carry more than its text: the sheet check's basis and which names it matched by
    return [dict({"id": n[0], "text": n[1]}, **(n[2] if len(n) > 2 else {})) for n in notes]


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


def report_section(source=None, session=None):
    """The report's Submission section -- {"html", "md", "label"} -- or None when the run has
    no submission. html and md say the same thing (the HTML report and its Markdown twin). A
    record that exists but cannot be read is shown as such rather than left out; its label is
    then None, so no header names a submission nobody could read."""
    try:
        rec, sess = resolve(source, session)
    except RecordError as e:
        import report_style
        msg = f"The submission record could not be read: {e}"
        return {"html": report_style.callout("warning", f"<p>{_h(msg)}</p>"),
                "md": f"> **Warning.** {msg}", "label": None}
    if rec is None:
        return None
    notes = quality_notes(rec, sess)
    return {"html": render_html(rec, notes),
            "md": render_markdown(rec, notes, heading="").lstrip("\n"), "label": label(rec)}


def html_section(source=None, session=None):
    """The report's Submission section as HTML, or "" when the run has no submission."""
    sec = report_section(source, session)
    return sec["html"] if sec else ""


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
    if a.record:
        # locate writes locate.json beside the record (`locate --out` = the fetch folder): the
        # staff decisions behind the file list come with it
        dec = locate_decisions(os.path.join(os.path.dirname(path), "locate.json"))
        if dec is not None:
            lp = paths_for(session)["locate_decisions"]
            _write_json(lp, dec)
            paths["locate_decisions"] = lp
    nm_level, nm_note = write_name_map(session, a.service_project)
    if nm_level:
        paths["raw_names"] = {"level": nm_level, "note": nm_note}
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
    at.add_argument("--service-project", default=None,
                    help="the staged service project, for a session that is not inside it (as "
                         "deliver's): its raw-name map is then copied to input/raw_names.tsv")
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
