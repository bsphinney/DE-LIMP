#!/usr/bin/env python3
"""
qc_bracket.py  --  Was the instrument healthy while this project's samples were acquired? The
QC runs around the project, from STAN.

The Core injects a QC standard on each instrument between projects, and STAN
(https://ucd.stan-proteomics.org) searches and scores every one of them. A project's results
are only as good as the instrument was while its samples ran, so this finds, for each
instrument the project used:

  * the nearest QC run BEFORE the first sample, and how long before;
  * every QC run DURING the project (between the first and the last sample);
  * the nearest QC run AFTER the last sample, or "none yet";

and grades each one against STAN's Instrument Performance Score (IPS) and against that
instrument's own QC runs over the 30 days before the project. Then it gives one verdict:
good / check / concern / no_qc_on_record.

STAFF-ONLY (Brett, 2026-10-01). The verdict and everything behind it are for Core staff: the
record and its staff page live in the session's logs/ (Core-internal: never delivered, kept out
of the session zip) and in the run registry (record_run.py). No client deliverable carries a
verdict word -- the client report has no QC section, and the Methods get one fixed sentence
(METHODS_SENTENCE, a statement of the Core's practice), the same whatever the verdict. A verdict
of check or concern, or any file that could not be placed in time, HOLDS the client delivery
(core_submission.py deliver asks delivery_gate()) until a staff member records an
acknowledgement -- who, when, a one-line note -- for this record. A record for a different
file list from the session's (files added or removed since) is refused outright: re-run.
  python3 qc_bracket.py ack --session <S> --by <HIVE login> --note "<what was reviewed, decided>"
      [--no-record-reason "<why the check could not run>"]   # only when there is no record
  python3 qc_bracket.py gate --session <S>     # exit 0 = may go out, 2 = needs an acknowledgement
`--by` is checked by staff.py, the one check behind every Core decision: a login on the Core
staff list ($CORE_STAFF_FILE, default /quobyte/proteomics-grp/.config/core_staff.txt, trusted
only when a Core admin owns it and its folder); with no list, agent-like names are refused and
the acknowledgement is marked NOT VERIFIED. This stops an agent acknowledging by accident; that a
PERSON did it comes later, with the board's Duo approval.
Every staff-only output carries STAFF_ONLY_MARKER, and the session zip and deliver leave out any
file carrying it, whatever its name; deliver writes the full gate to logs/qc_gate.json and only
proceed / record_sha256 / ack_at to delivery.json.
Run it only on a Core STAFF route (SKILL.md 8e): never for a client, never shown to one.

WHAT IT READS
  The project's files: their instrument and acquisition START time (acq_time.py) -- from
  --acquisition-json (detect_acquisition.py's output, step 2: nothing is re-read), else from the
  files (a .d's analysis.tdf; a .raw through ThermoRawFileParser, ~3-7 s each).
  STAN: GET <url>/api/runs?instrument=<model>&limit=500&offset=<k>&qc_only=true, anonymous,
  newest first, paged until the rows are older than the baseline (stan/dashboard/server.py
  api_runs -> stan/db.py get_runs). qc_only is STAN's own QC-name rule (stan/watcher/
  qc_filter.py); runs an operator hid in STAN are left out, as STAN's dashboard leaves them out.
  qc_only filters each page after fetching 3x its limit, so an EMPTY page can come mid-history;
  paging stops only after EMPTY_PAGES_STOP empty pages in a row.

THE MAINTENANCE LOG
  STAN's GET <url>/api/instruments/<model>/events (anonymous, newest first): column changes and
  clogs, source cleaning, calibration, LC service, downtime (with an end date). Long QC gaps are
  often the instrument being down, not QC missed. For each instrument:
  * every event within 30 days either side of the project is listed ("column change
    2026-07-30 (date only)", "instrument downtime 2026-06-25 to 2026-07-21");
  * every gap of more than 7 days between the QC runs around the project is explained by the
    events that overlap it -- quoted from the log, never a guessed cause -- or said plainly:
    "no QC for N days; STAN has no maintenance record for this period" (and, for an instrument
    with an empty log, that nothing is logged for it at all: on 2026-10-01 the Exploris 480
    had no events and the Fusion Lumos none from May to September, so no record is not no
    downtime);
  * the first QC run after an event is marked as the post-maintenance check. When project
    files ran after the event with no other QC run between them and it, it speaks for the
    project however far away it is, and a flag on it is never discounted as isolated. Files
    acquired after an event and before any QC run, or inside a logged downtime, are a check.
  A date with no time (STAN stores one as noon UTC) is a whole Pacific day: files on that day
  are "order unknown", and its first QC "on or after the day" may predate the work. An event's
  first_run, when it names one of this project's files or one of STAN's QC runs, pins the
  change to that run. Notes, operator and creator are never copied: they can name people.
  A log that cannot be read leaves the verdict to the QC runs and says the gaps are unexplained.

WHAT IT DOES NOT USE, ON PURPOSE
  `gate_result`: STAN records "pass" for every run, including ones that identified nothing
  (STAN docs/qc_gating_and_slack_summary.md: no thresholds file exists, so the gate never
  fails). IPS and the run's own numbers are the QC signal.

HOW A QC RUN IS GRADED (each threshold is a constant below, with its basis)
  concern  no identifications; IDs under 50% of the 30-day median
  check    IDs under 80% of the 30-day median; MS1 or MS2 mass error over 1.5x the median (and
           at least 1 ppm above it); chromatographic peaks (FWHM) over 1.5x the median width;
           too few baseline runs to compare it with (fewer than 3 same-mode QC runs in the 30
           days before -- a sparse schedule, or the month after a downtime): nothing to compare
           with is not "normal"
  IPS      shown with STAN's dashboard colour (green >= 80, amber >= 60, red), per run and as
           the instrument's 30-day median -- NOT graded. IPS compares a run with STAN's April
           2026 reference cohort, not with the instrument's recent self, and on 2026-10-01 62%
           of STAN's 4,683 QC runs scored red; the Exploris 480's QC at 38 samples/day had a
           median IPS of 28 over the past year. Graded, it would call most projects on that
           instrument a concern while its QC ran exactly as usual.
  The 30-day median is taken over the same instrument's QC runs of the same acquisition mode
  (DIA/DDA) and throughput (samples per day): a 35-min and a 120-min gradient differ by 1.7x
  in IDs on the same instrument (STAN docs). With fewer than 3 such runs it widens to every
  run of that mode, and says so. Runs named as blanks are left out even when STAN's QC-name
  rule matched them ("blank...HELA": a blank before a QC, zero IDs by design).

THE VERDICT, per instrument and overall (the worst)
  The QC runs that speak for the project are the ones DURING it and the nearest before and
  after, if within 7 days. A flagged run whose neighbours in STAN's sequence -- the QC runs
  just before and just after it -- are both normal AND within 7 days of it (ISOLATION_DAYS) is an
  ISOLATED failure (a bad injection; the instrument worked on both sides of it): it is shown,
  and counts one level lower. A normal run weeks away clears nothing.
  concern          a QC run that speaks for the project graded concern
  check            one graded check (or an isolated concern); or no QC run within 7 days of
                   the project at all, so the instrument's state while it ran was not measured
  no_qc_on_record  STAN has no QC run for that instrument
  good             otherwise
  Overall order: concern > no_qc_on_record > check > good.
  Replayed over STAN's last year of QC and its maintenance log (2026-10-01; 1,000 random
  projects of 2 h to 3 days per instrument; re-run after the 2.10 review's fixes), these rules
  gave good / check / concern 43 / 29 / 28% on the timsTOF HT (13% of its QC runs identified under half its usual,
  often several in a row), 40 / 35 / 24% on the Fusion Lumos and 41 / 46 / 13% on the Exploris
  480 (16% of Exploris projects had no QC run within 7 days).

WHAT IT CANNOT KNOW, AND SAYS
  * STAN records the instrument MODEL, not its serial number, so files are matched to QC runs
    by model ("timsTOF HT"). This is a check for runs acquired on the UC Davis Core's
    instruments; the same model elsewhere is not the same instrument.
  * STAN does not record whether a run's time came from the file's header or, when it could
    not read that, from the file's modification time (stan/pipeline/hive_process.py). For the
    Core's Orbitrap QC runs it is often the modification time -- the END of the run, up to a
    gradient length after its start (two checked on 2026-10-01 were). The project's own files
    are timed from the START of acquisition.
  * A Thermo .raw records no time zone (acq_time.py: read as the Core's Pacific time).
  * Files with no readable time or instrument are listed, not guessed.

EXIT CODES -- the orchestrator branches on these
  0  checked: the record says the verdict
  2  usage error (no raw file list, no files), or refused (too many .raw to read on a cluster
     login node) -- NO record is written
  3  STAN could not be reached, or did not answer with its run list: NO verdict. The record
     says so; staff must hear the check could not be made -- never "good"
  4  nothing to check: no file gave both an instrument and an acquisition time
  (`gate` / `ack`: 0 = the client report may go out, 2 = it needs an acknowledgement / bad input)

The record (JSON, schema_version 1 -- the Core's board reads it) goes to stdout and to --out
(with --session: <session>/logs/qc_bracket.json), with the staff page beside it
(logs/qc_bracket.md); record_run.py copies both into the run registry, make_methods.py takes only
its methods_sentence, provenance.py only a pointer and checksum. The plain-English summary also
goes to stderr.

Usage:
  python3 qc_bracket.py --session <session dir>          # input/raw_files.txt (else the
                                                         # search's record); input/acquisition.json
  python3 qc_bracket.py --files /data/*.d /data/*.raw [--out qc_bracket.json]
      [--acquisition-json step2.json] [--stan-url https://ucd.stan-proteomics.org]
      [--timeout 30] [--allow-login-node]
"""
import argparse
import glob
import http.client
import json
import os
import re
import statistics
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import acq_time  # noqa: E402  -- when a run was acquired, and its time zone: one reader
import staff  # noqa: E402  -- who may sign a Core decision: the one check

SCHEMA = "qc_bracket/1"
# The board reads these records: a field is never renamed or removed without bumping this.
SCHEMA_VERSION = 1
STAN_URL = "https://ucd.stan-proteomics.org"
STAN_URL_ENV = "STAN_URL"
RUNS_PATH = "/api/runs"
PAGE = 500                  # rows per request; verified to page on 2026-10-01
MAX_PAGES = 40              # per instrument: 20,000 QC rows (STAN held 1,379-1,726 each)
TIMEOUT_S = 30
BASELINE_DAYS = 30          # the instrument's "usual", the month before the project
MIN_BASELINE = 3            # fewer same-method runs than this: widen to the whole mode; fewer
                            # same-mode runs than this: no comparison, and the run is a CHECK
                            # (2.10 review: it used to grade "ok", so a run at 10% of the
                            # usual IDs read GOOD after a sparse month or a downtime)
NEAR_DAYS = 7               # a QC run this close to the project speaks for it
# A flagged run is ISOLATED only when the QC runs on both sides of it are normal AND within this
# many days of it (2.10 review: neighbours 20 days away cleared a check on the only QC run near
# the project, and the verdict read GOOD).
ISOLATION_DAYS = 7
# With qc_only, STAN fetches 3 x limit rows at `offset`, keeps the QC-named ones and returns the
# first `limit` (stan/db.py get_runs): a page is EMPTY wherever 3 x limit consecutive rows are not
# QC runs, not only at the end of the table (2.10 review). Paging stops after this many empty
# pages in a row -- with PAGE 500, 3,000+ consecutive rows holding no QC run.
EMPTY_PAGES_STOP = 4
# STAN rows dated before this, or more than a day in the future, are not acquisition times
# (STAN holds a Fusion Lumos run dated 1980-01-02): listed, never bracketed.
OLDEST_PLAUSIBLE = datetime(2000, 1, 1, tzinfo=timezone.utc)
# IPS bands, exactly as STAN's dashboard colours them (stan/dashboard/public/index.html
# IpsBadge: >= 80 green, >= 60 amber, else red). 60 is a median run of STAN's April 2026
# reference cohort (STAN docs/ips_metric.md).
IPS_GREEN, IPS_AMBER = 80, 60
# docs/ips_metric.md: "<30 = you underperformed the bottom 10% -- something is wrong". Named
# in a run's notes; not graded (see the docstring).
IPS_FAILING = 30
# A run name with one of these is not a QC standard, whatever STAN's QC-name rule matched.
NOT_QC_WORDS = ("blank",)
# Fractions of the 30-day median. Measured on STAN's 4,683 QC runs (2026-10-01), each against
# its own instrument's same-throughput runs in the 30 days before it: 5-10% of runs fall under
# 0.5 (failed or near-failed injections) and 10-25% under 0.8.
IDS_CONCERN, IDS_CHECK = 0.5, 0.8
# Mass error (STAN's median |ppm|) and peak width (FWHM): about the 95th percentile of the
# same ratio in that data. The ppm floor keeps a 0.5 -> 0.8 ppm wobble on an Orbitrap quiet.
PPM_RATIO, PPM_FLOOR = 1.5, 1.0
FWHM_RATIO = 1.5
# STAN's maintenance log: GET /api/instruments/{instrument}/events (stan/dashboard/server.py
# api_events -> db.get_events), anonymous, newest first, no date filter. On 2026-10-01 it held 11
# events for the timsTOF HT, 4 for the Fusion Lumos (none May-September) and 0 for the Exploris
# 480: Orbitrap downtime is not logged there yet, so "no record" is not "no downtime".
EVENTS_PATH = "/api/instruments/{}/events"
EVENTS_LIMIT = 500
# QC runs around the project further apart than this: a gap the report explains from the log.
GAP_DAYS = 7
# A bare date, or the noon-UTC time STAN's older Trends form stores a date as (server.py
# api_log_event: "2026-08-10T12:00:00Z" vs a date input's "2026-08-10"), is a calendar day in
# the lab's time zone, not an instant.
DATE_ONLY = re.compile(r"^(\d{4}-\d\d-\d\d)(?:T12:00:00(?:\.0+)?(?:Z|\+00:00))?$")
# stan/db.py EVENT_TYPES, as a reader says them. A type not here is shown as it is spelled.
EVENT_LABEL = {"column_change": "column change", "column_clog": "column clog",
               "capillary_change": "capillary change", "emitter_change": "emitter change",
               "source_clean": "ion source cleaning", "calibration": "mass calibration",
               "pm": "preventive maintenance", "lc_service": "LC service",
               "downtime": "instrument downtime", "other": "maintenance (other)"}
DOWNTIME_TYPES = ("downtime",)          # stan/db.py DOWNTIME_EVENT_TYPES
# Never copied out of an event: notes, operator, created_by -- they can name people and
# customers (STAN shares none of them without an opt-in) -- and first_run, which is MATCHED to
# this project's files and STAN's QC runs and shown only as one of those.
EXIT_OK, EXIT_USAGE, EXIT_STAN, EXIT_NOTHING = 0, 2, 3, 4
RANK = {"good": 0, "check": 1, "no_qc_on_record": 2, "concern": 3}
LABEL = {"good": "GOOD", "check": "CHECK", "concern": "CONCERN",
         "no_qc_on_record": "NO QC ON RECORD"}
GRADE_RANK = {"ok": 0, "check": 1, "concern": 2}
ONE_LOWER = {"concern": "check", "check": "ok", "ok": "ok"}


# ----------------------------------------------------------------------------- the files --
def _norm(p):
    return str(p).rstrip("/\\")


def load_acquisition_json(path):
    """detect_acquisition.py's output -> ({file: entry}, problem). Entries without the
    acquired_at fields (an older skill's step 2) are dropped, so those files are read."""
    if not path:
        return {}, None
    try:
        with open(path, encoding="utf-8-sig") as fh:
            d = json.load(fh)
    except (OSError, ValueError) as e:
        return {}, f"{path} could not be read ({type(e).__name__}: {e}); the files are read instead"
    files = d.get("files") if isinstance(d, dict) else None
    if not isinstance(files, list):
        return {}, f"{path} has no `files` list (not detect_acquisition.py output)"
    return ({_norm(f["file"]): f for f in files if isinstance(f, dict) and f.get("file")
             and "acquired_at" in f}, None)


def identify(path, known=None):
    """{file, instrument, acquired_at, acquired_at_source, acquired_at_note, read_from} for one
    file: from step 2's entry when it has one, else read here."""
    out = {"file": path, "instrument": None, **acq_time.result(), "read_from": None}
    if known:
        out.update({k: known.get(k) for k in ("instrument", "acquired_at", "acquired_at_source",
                                              "acquired_at_note")})
        out["read_from"] = "detect_acquisition.py output (--acquisition-json)"
        return out
    import detect_acquisition as da      # only when a file must be read: it is heavy
    if not os.path.exists(path):
        out["acquired_at_note"] = "not found on this computer"
        return out
    low = path.lower()
    if low.endswith(".d") and os.path.isdir(path):
        out["instrument"] = da.detect_instrument(path)
        out.update(acq_time.read_bruker(path))
        out["read_from"] = "the .d (analysis.tdf GlobalMetadata)"
    elif da.is_thermo_raw(path):
        meta, why = da.trfp_metadata(path)
        if meta is None:
            out["acquired_at_note"] = f"not read: {why}"
        else:
            out["instrument"] = da._cv(meta, "InstrumentProperties", da.CV_MODEL)
            out.update(da.thermo_acquired_at(meta, path))
        out["read_from"] = "the .raw (ThermoRawFileParser metadata)"
    else:
        out["acquired_at_note"] = "not read: only a Bruker .d and a Thermo .raw carry a time here"
    return out


# ------------------------------------------------------------------------------- STAN --
class StanUnreachable(Exception):
    pass


def get_json(url, timeout):
    """GET url -> parsed JSON. The one network call; the tests replace it."""
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "ucdavis-proteomics-skill qc_bracket.py"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def stan_runs(base, instrument, oldest_needed, timeout, log):
    """Every QC run STAN has for `instrument` from now back past `oldest_needed`, newest first,
    de-duplicated (with qc_only, STAN filters each page after fetching it, so consecutive
    offsets can repeat a row: offsets step by PAGE, never by the rows returned). Raises
    StanUnreachable with the reason; `log` gets one line per page. An empty page is not the end
    (see EMPTY_PAGES_STOP)."""
    rows, seen, pages, empty = [], set(), 0, 0
    for k in range(MAX_PAGES):
        qs = urllib.parse.urlencode({"instrument": instrument, "limit": PAGE,
                                     "offset": k * PAGE, "qc_only": "true"})
        url = f"{base}{RUNS_PATH}?{qs}"
        try:
            page = get_json(url, timeout)
        except urllib.error.HTTPError as e:
            raise StanUnreachable(f"{base} answered HTTP {e.code} for {RUNS_PATH} "
                                  f"(instrument {instrument!r})")
        except urllib.error.URLError as e:
            raise StanUnreachable(f"cannot reach {base}: {e.reason}")
        except (OSError, ValueError, http.client.HTTPException) as e:  # timeout, reset, cut short, not JSON
            raise StanUnreachable(f"{base}{RUNS_PATH} did not answer with JSON "
                                  f"({type(e).__name__}: {e})")
        if not isinstance(page, list) or not all(isinstance(r, dict) for r in page):
            raise StanUnreachable(f"{base}{RUNS_PATH} answered something that is not a list of "
                                  f"runs ({type(page).__name__})")
        pages += 1
        if not page:
            empty += 1
            if empty >= EMPTY_PAGES_STOP:
                break               # the end of the table, or 3,000+ rows with no QC run
            continue                # an empty page mid-history: older QC runs may follow
        empty = 0
        for r in page:
            if r.get("id") not in seen:
                seen.add(r.get("id"))
                rows.append(r)
        times = [_aware(t) for t in (acq_time.parse_iso(r.get("run_date")) for r in page) if t]
        log(f"  {instrument}: page {pages}, {len(page)} rows, oldest "
            f"{min(times).date() if times else '?'}")
        if times and min(times) < oldest_needed:
            break
    else:
        log(f"  {instrument}: stopped after {MAX_PAGES} pages")
    return rows, pages


def stan_events(base, instrument, timeout):
    """(STAN's maintenance events for `instrument`, None) or (None, why). A log that cannot be
    read does not stop the check -- the QC runs still speak -- but every gap says so."""
    url = (f"{base}{EVENTS_PATH.format(urllib.parse.quote(instrument, safe=''))}?"
           f"{urllib.parse.urlencode({'limit': EVENTS_LIMIT})}")
    try:
        got = get_json(url, timeout)
    except urllib.error.HTTPError as e:
        return None, f"{base} answered HTTP {e.code} for the maintenance log"
    except urllib.error.URLError as e:
        return None, f"cannot reach {base}: {e.reason}"
    except (OSError, ValueError, http.client.HTTPException) as e:
        return None, f"the maintenance log did not answer with JSON ({type(e).__name__}: {e})"
    if not isinstance(got, list) or not all(isinstance(e, dict) for e in got):
        return None, "the maintenance log answered something that is not a list of events"
    return got, None


def _aware(t):
    return t if t.tzinfo else t.replace(tzinfo=timezone.utc)


# -------------------------------------------------------------------------- grading --
def _mode(r):
    m = str(r.get("mode") or "").lower()
    return "dda" if "dda" in m else "dia" if "dia" in m else m or None


def ids_of(r):
    """(identifications, what they count): precursors for DIA, PSMs for DDA."""
    if _mode(r) == "dda":
        return (r.get("n_psms") if r.get("n_psms") is not None else r.get("n_precursors")), "PSMs"
    return (r.get("n_precursors") if r.get("n_precursors") is not None else r.get("n_psms")), \
        "precursors"


def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if v == v else None


def ips_band(ips):
    if ips is None:
        return None
    return "green" if ips >= IPS_GREEN else "amber" if ips >= IPS_AMBER else "red"


def baseline_for(r, base):
    """(runs, label): the baseline runs this QC run is compared with -- same mode and
    throughput, else same mode -- and how they were chosen."""
    same_mode = [b for b in base if _mode(b) == _mode(r)]
    if r.get("spd") is not None:
        same = [b for b in same_mode if b.get("spd") == r.get("spd")]
        if len(same) >= MIN_BASELINE:
            return same, f"{len(same)} {(_mode(r) or '').upper()} QC runs at {r['spd']} samples/day"
    if len(same_mode) >= MIN_BASELINE:
        return same_mode, (f"{len(same_mode)} {(_mode(r) or '').upper()} QC runs, all throughputs"
                           + (f" (fewer than {MIN_BASELINE} at {r['spd']} samples/day)"
                              if r.get("spd") is not None else ""))
    return [], (f"fewer than {MIN_BASELINE} {(_mode(r) or 'same-mode').upper()} QC runs in the "
                f"{BASELINE_DAYS} days before the project: no comparison")


def _median(vals):
    vals = [v for v in vals if v is not None and v > 0]
    return statistics.median(vals) if vals else None


def grade_run(r, base):
    """One QC run, graded: its numbers, each against the baseline median, and its flags."""
    ids, unit = ids_of(r)
    ips = _num(r.get("ips_score"))
    ips = int(ips) if ips is not None else None
    runs, cohort = baseline_for(r, base)
    med_ids = _median([_num(ids_of(b)[0]) for b in runs])
    ms1, ms2 = _num(r.get("median_mass_acc_ms1_ppm")), _num(r.get("median_mass_acc_ms2_ppm"))
    fwhm = _num(r.get("fwhm_rt_min"))
    med = {"ids": med_ids,
           "ms1_ppm": _median([_num(b.get("median_mass_acc_ms1_ppm")) for b in runs]),
           "ms2_ppm": _median([_num(b.get("median_mass_acc_ms2_ppm")) for b in runs]),
           "fwhm_s": _median([(_num(b.get("fwhm_rt_min")) or 0) * 60 for b in runs])}
    flags = []

    def flag(code, severity, text):
        flags.append({"code": code, "severity": severity, "text": text})

    if not ids:
        flag("zero_ids", "concern", f"no {unit} identified: the QC injection failed, or its "
                                    f"search found nothing")
    if ips is None:
        flag("ips_missing", "note", "STAN recorded no IPS for this run")
    else:
        flag(f"ips_{ips_band(ips)}", "note", f"IPS {ips} ({ips_band(ips)} on STAN's dashboard"
             + (f"; under {IPS_FAILING}, STAN's 'something is wrong' level" if ips < IPS_FAILING
                else "") + ")")
    ratio = (ids / med_ids) if ids and med_ids else None
    if ids and ratio is None:
        # nothing to compare with is not "normal": a staff member looks
        flag("no_baseline", "check", f"too few baseline runs to compare it with ({cohort}), so "
                                     f"it cannot be called normal")
    if ratio is not None:
        if ratio < IDS_CONCERN:
            flag("ids_low", "concern", f"{ids:,} {unit}: {ratio:.0%} of the instrument's 30-day "
                                       f"median ({med_ids:,.0f})")
        elif ratio < IDS_CHECK:
            flag("ids_low", "check", f"{ids:,} {unit}: {ratio:.0%} of the instrument's 30-day "
                                     f"median ({med_ids:,.0f})")
    for key, val, name in (("ms1_ppm", ms1, "MS1"), ("ms2_ppm", ms2, "MS2")):
        m = med[key]
        if ids and val and m and val > max(PPM_RATIO * m, m + PPM_FLOOR):
            flag(f"{key}_high", "check", f"{name} mass error {val:.2f} ppm vs {m:.2f} ppm typical "
                                         f"over the 30 days before")
    fwhm_s = fwhm * 60 if fwhm else None
    if ids and fwhm_s and med["fwhm_s"] and fwhm_s > FWHM_RATIO * med["fwhm_s"]:
        flag("fwhm_wide", "check", f"chromatographic peaks {fwhm_s:.1f} s wide (FWHM) vs "
                                   f"{med['fwhm_s']:.1f} s typical")
    worst = max((f["severity"] for f in flags if f["severity"] in GRADE_RANK),
                key=GRADE_RANK.get, default="ok")
    return {"stan_id": r.get("id"), "run_name": r.get("run_name"),
            "run_date": acq_time.iso_utc(_aware(acq_time.parse_iso(r["run_date"]))),
            "local_time": acq_time.local_text(_aware(acq_time.parse_iso(r["run_date"]))),
            "mode": r.get("mode"), "samples_per_day": r.get("spd"),
            "gradient_min": r.get("gradient_length_min"),
            "ips": ips, "ips_band": ips_band(ips), "ids": ids, "ids_unit": unit,
            "ids_vs_median": round(ratio, 3) if ratio is not None else None,
            "ms1_ppm": ms1, "ms2_ppm": ms2, "fwhm_s": round(fwhm_s, 1) if fwhm_s else None,
            "baseline": {"runs": cohort, "median_ids": med_ids, "median_ms1_ppm": med["ms1_ppm"],
                         "median_ms2_ppm": med["ms2_ppm"],
                         "median_fwhm_s": round(med["fwhm_s"], 1) if med["fwhm_s"] else None},
            "grade": worst, "flags": flags}


def samples_around(x, order, graded_run, files):
    """The project's files acquired between the last NORMAL QC run before a flagged one and the
    first normal QC run after it: the samples a staff member should look at."""
    i = next(k for k, r in enumerate(order) if r.get("id") == x["stan_id"])
    prev = next((r for r in reversed(order[:i]) if graded_run(r)["grade"] == "ok"), None)
    nxt = next((r for r in order[i + 1:] if graded_run(r)["grade"] == "ok"), None)
    lo = _aware(acq_time.parse_iso(prev["run_date"])) if prev else None
    hi = _aware(acq_time.parse_iso(nxt["run_date"])) if nxt else None
    hit = sorted((acq_time.parse_iso(f["acquired_at"]), f["file"]) for f in files
                 if (lo is None or acq_time.parse_iso(f["acquired_at"]) > lo)
                 and (hi is None or acq_time.parse_iso(f["acquired_at"]) < hi))
    return {"from": acq_time.local_text(lo) if lo else "no normal QC run before it",
            "to": acq_time.local_text(hi) if hi else "no normal QC run after it yet",
            "n_files": len(hit), "files": [os.path.basename(f) for _t, f in hit]}


def _days(a, b):
    return round(abs((a - b).total_seconds()) / 86400, 1)


def _gap_text(days):
    if days < 1:
        h = round(days * 24)
        return "under an hour" if h < 1 else f"{h} hour{'s' if h != 1 else ''}"
    d = round(days)
    return f"{d} day{'s' if d != 1 else ''}"


# ---------------------------------------------------------------- the maintenance log --
def _event_time(text, end=False):
    """(aware datetime, date_only) for a STAN event_date / end_date, or (None, None). A calendar
    day (DATE_ONLY) gives its start in the lab's zone, or with end=True its end; a time with no
    offset is the lab's local time (make_methods.py reads event_date the same way)."""
    if not isinstance(text, str) or not text.strip():
        return None, None
    m = DATE_ONLY.match(text.strip())
    if m:
        day = datetime.strptime(m.group(1), "%Y-%m-%d") + timedelta(days=1 if end else 0)
        t, _how = acq_time.localize(day)
        return (t, True) if t is not None else (None, None)
    t = acq_time.parse_iso(text)
    if t is not None and t.tzinfo is None:
        t, _how = acq_time.localize(t)
    return (t, False) if t is not None else (None, None)


def _stem(name):
    b = os.path.basename(str(name or "").strip().rstrip("/\\")).lower()
    return os.path.splitext(b)[0] if b.endswith((".d", ".raw", ".mzml")) else b


def _day(t):
    return str(acq_time.local_date(t))


def read_events(events, files, dated):
    """STAN's maintenance events, as spans: [{type, label, start, end, ready, date_only,
    anchor}], and the ones that could not be read. `ready` is when the instrument is in its
    after-the-event state: the first run on the new column when `first_run` names one of this
    project's files or one of STAN's QC runs, else the end of a logged downtime, else the
    event (for a calendar day, its start: the first QC run on or after that day)."""
    by_file = {_stem(f["file"]): f for f in files}
    by_qc = {_stem(r.get("run_name")): (t, r) for t, r in dated}
    out, bad = [], []
    for e in events:
        start, date_only = _event_time(e.get("event_date"))
        if start is None:
            bad.append({"event_id": e.get("id"), "why": "event_date is not a date"})
            continue
        end_given = bool(e.get("end_date"))
        end, end_day = (_event_time(e.get("end_date"), end=True) if end_given else (None, None))
        if end is None:
            end_given = False
            end, end_day = (_event_time(e.get("event_date"), end=True) if date_only
                            else (start, False))
        ready = end - timedelta(days=1) if end_day else end
        if end_day and not end_given:
            ready = start
        anchor, fr = None, _stem(e.get("first_run"))
        if fr and fr in by_file:
            f = by_file[fr]
            anchor = {"what": f"this project's file {os.path.basename(f['file'])}",
                      "at": acq_time.parse_iso(f["acquired_at"])}
        elif fr and fr in by_qc:
            anchor = {"what": "a QC run in STAN", "at": by_qc[fr][0]}
        if anchor:
            ready = anchor["at"]
        etype = str(e.get("event_type") or "other")
        label = EVENT_LABEL.get(etype, etype.replace("_", " "))
        if etype == "column_change" and e.get("column_model"):
            label += f" ({e['column_model']})"          # equipment, never a person
        a = _day(start) if date_only else acq_time.local_text(start)
        b = (_day(end - timedelta(days=1)) if end_day else acq_time.local_text(end)) \
            if end_given else a
        when = (f"{a} to {b}" if b != a else
                a + (" (date only)" if date_only and not anchor else ""))
        out.append({"type": etype, "label": label, "when": when, "start": start, "end": end,
                    "ready": ready, "date_only": bool(date_only) and not anchor,
                    "anchor": anchor, "downtime": etype in DOWNTIME_TYPES})
    return out, bad


def maintenance(events, error, files, dated, graded_run, first, last, now):
    """The maintenance log around the project: its events within BASELINE_DAYS either side, each
    with how this project's files fall around it and the first QC run after it (the
    post-maintenance check), and every events-related verdict note."""
    rec = {"log": "unreadable" if error else "read", "error": error,
           "n_logged": None if error else len(events or []), "events": [],
           "unreadable_events": [], "unchecked_files": 0, "files_during_downtime": 0}
    if error:
        return rec, []
    spans, rec["unreadable_events"] = read_events(events or [], files, dated)
    lo, hi = first - timedelta(days=BASELINE_DAYS), last + timedelta(days=BASELINE_DAYS)
    ftimes = sorted((acq_time.parse_iso(f["acquired_at"]), f["file"]) for f in files)
    for ev in sorted((x for x in spans if x["start"] <= hi and x["end"] >= lo),
                     key=lambda x: x["start"]):
        cut = ev["anchor"]["at"] if ev["anchor"] else ev["end"]
        before_ev = ev["anchor"]["at"] if ev["anchor"] else ev["start"]
        after = [f for t, f in ftimes if t >= cut]
        inside = [f for t, f in ftimes if before_ev <= t < cut]
        k = next((i for i, (t, _r) in enumerate(dated) if t >= ev["ready"]), None)
        post = dated[k] if k is not None else None
        nxt = dated[k + 1][0] if k is not None and k + 1 < len(dated) else None
        # the files this QC run is the check for: acquired after the event, with no other QC
        # run between them and it (a change a month before, with QC every day since, is not)
        gated = [f for t, f in ftimes if t >= cut and (nxt is None or t < nxt)]
        first_qc = None
        if post:
            first_qc = dict(graded_run(post[1]), after_maintenance=ev["label"],
                            days_from_event=_days(post[0], ev["ready"]))
        unchecked = [f for t, f in ftimes if t >= cut and (post is None or t < post[0])]
        rec["events"].append({
            "type": ev["type"], "label": ev["label"], "when": ev["when"],
            "start": acq_time.iso_utc(ev["start"]), "end": acq_time.iso_utc(ev["end"]),
            "date_only": ev["date_only"],
            "anchored_by": ev["anchor"]["what"] if ev["anchor"] else None,
            "where": ("before the project" if ev["end"] <= first else
                      "after the project" if ev["start"] > last else "during the project"),
            "files_after": len(after),
            # project files whose QC check this event's first QC run is (see `gated`)
            "files_checked_by_first_qc": len(gated) if post else 0,
            "files_same_day": ([os.path.basename(f) for f in inside]
                               if ev["date_only"] and not ev["downtime"] else []),
            "files_during_downtime": ([os.path.basename(f) for f in inside]
                                      if ev["downtime"] else []),
            "files_unchecked_after": [os.path.basename(f) for f in unchecked],
            "first_qc_after": first_qc})
        rec["unchecked_files"] += len(unchecked)
        rec["files_during_downtime"] += len(inside) if ev["downtime"] else 0
    return rec, spans


def qc_gaps(o, spans, mrec, now):
    """Every stretch of more than GAP_DAYS without a QC run around the project, each explained
    from the maintenance log when an event overlaps it -- or said plainly when none does. The
    log is quoted, never a cause guessed."""
    seq = sorted((acq_time.parse_iso(x["run_date"]) for x in
                  [o["before"]] + list(o["during"]) + [o["after"]] if x))
    gaps = [(a, b, False) for a, b in zip(seq, seq[1:]) if b - a > timedelta(days=GAP_DAYS)]
    if seq and o["after"] is None and now - seq[-1] > timedelta(days=GAP_DAYS):
        gaps.append((seq[-1], now, True))
    out = []
    for a, b, open_ended in gaps:
        n = round((b - a).total_seconds() / 86400)
        head = (f"No QC since {_day(a)} ({n} days)" if open_ended else
                f"No QC for {n} days ({_day(a)} to {_day(b)})")
        hits = [x for x in spans if x["start"] <= b and x["end"] >= a]
        if mrec["log"] != "read":
            text = (f"{head}; STAN's maintenance log could not be read ({mrec['error']}), so "
                    f"this gap is not explained.")
        elif hits:
            text = (f"{head}. STAN's maintenance log has, in this period: "
                    + "; ".join(f"{x['label']} {x['when']}" for x in hits) + ".")
        elif not mrec["n_logged"]:
            text = (f"{head}; STAN has no maintenance record for this period (none is logged "
                    f"for this instrument at all).")
        else:
            text = f"{head}; STAN has no maintenance record for this period."
        out.append({"from": acq_time.iso_utc(a), "to": None if open_ended else acq_time.iso_utc(b),
                    "days": n, "events": [f"{x['label']} {x['when']}" for x in hits],
                    "text": text})
    return out


def bracket(instrument, files, rows, now, events=None, events_error=None):
    """The verdict for one instrument: its files' window, the QC runs around it, graded, and
    STAN's maintenance log around it."""
    times = sorted(acq_time.parse_iso(f["acquired_at"]) for f in files)
    first, last = times[0], times[-1]
    win = {"first": acq_time.iso_utc(first), "last": acq_time.iso_utc(last),
           "first_local": acq_time.local_text(first), "last_local": acq_time.local_text(last),
           "days": _days(first, last), "n_files": len(files)}
    ignored, dated = [], []
    for r in rows:
        t = acq_time.parse_iso(r.get("run_date"))
        if t is None:
            ignored.append({"stan_id": r.get("id"), "run_date": r.get("run_date"),
                            "why": "run_date is not a time"})
        elif not (OLDEST_PLAUSIBLE <= _aware(t) <= now + timedelta(days=1)):
            ignored.append({"stan_id": r.get("id"), "run_date": r.get("run_date"),
                            "why": "run_date is not a plausible acquisition time"})
        elif any(w in str(r.get("run_name") or "").lower() for w in NOT_QC_WORDS):
            ignored.append({"stan_id": r.get("id"), "run_date": r.get("run_date"),
                            "why": "named as a blank: not a QC standard"})
        else:
            dated.append((_aware(t), r))
    dated.sort(key=lambda x: x[0])
    base = [r for t, r in dated if first - timedelta(days=BASELINE_DAYS) <= t < first]
    before = [(t, r) for t, r in dated if t < first]
    during = [(t, r) for t, r in dated if first <= t <= last]
    after = [(t, r) for t, r in dated if t > last]
    out = {"instrument": instrument, "window": win, "stan_rows": len(rows),
           "ignored_rows": ignored,
           "baseline_window": {"from": acq_time.iso_utc(first - timedelta(days=BASELINE_DAYS)),
                               "to": acq_time.iso_utc(first), "n_runs": len(base)},
           "before": None, "during": [], "after": None}
    if not dated:
        mrec, _spans = maintenance(events, events_error, files, dated, None, first, last, now)
        out["maintenance"] = mrec
        out["verdict"] = "no_qc_on_record"
        out["summary_lines"] = [
            f"{instrument}: {LABEL['no_qc_on_record']}. STAN has no QC run on record for this "
            f"instrument" + (f" ({len(ignored)} row(s) with no usable date or named as blanks "
                             f"left out)" if ignored else "")
            + ", so its performance during this project cannot be checked from QC."]
        out["summary_lines"] += maintenance_lines(mrec)
        out["summary"] = " ".join(out["summary_lines"])
        return out
    order = [r for _t, r in dated]
    when = [t for t, _r in dated]
    pos = {id(r): i for i, r in enumerate(order)}
    graded = {}

    def graded_run(r):
        """grade_run(), once per run, plus whether a flagged run is isolated: the QC runs on
        either side of it in STAN's sequence are both normal and both within ISOLATION_DAYS of
        it -- a normal run weeks away says nothing about the instrument at the time."""
        if id(r) not in graded:
            g = grade_run(r, base)
            graded[id(r)] = g
            if g["grade"] != "ok":
                i = pos[id(r)]
                sides = [j for j in (i - 1, i + 1) if 0 <= j < len(order)]
                g["isolated"] = len(sides) == 2 and all(
                    abs(when[j] - when[i]) <= timedelta(days=ISOLATION_DAYS)
                    and grade_run(order[j], base)["grade"] == "ok" for j in sides)
                g["counts_as"] = ONE_LOWER[g["grade"]] if g["isolated"] else g["grade"]
            else:
                g["isolated"], g["counts_as"] = False, "ok"
        return graded[id(r)]

    if before:
        t, r = before[-1]
        out["before"] = dict(graded_run(r), days_from_project=_days(t, first))
    out["during"] = [graded_run(r) for _t, r in during]
    if after:
        t, r = after[0]
        out["after"] = dict(graded_run(r), days_from_project=_days(t, last))
    else:
        out["after_note"] = (f"none yet: STAN has no QC run on this instrument after the last "
                             f"sample ({_gap_text(_days(now, last))} ago)")
    ips = [x for x in (_num(r.get("ips_score")) for r in base) if x is not None]
    out["baseline_window"]["median_ips"] = statistics.median(ips) if ips else None
    out["baseline_window"]["median_ips_band"] = ips_band(statistics.median(ips)) if ips else None
    mrec, spans = maintenance(events, events_error, files, dated, graded_run, first, last, now)
    out["maintenance"] = mrec
    near = list(out["during"]) + [x for x in (out["before"], out["after"])
                                  if x and x["days_from_project"] <= NEAR_DAYS]
    out["desert"] = not near
    # The first QC run after a maintenance event is the check for every file acquired after the
    # event: it speaks for the project however far away it is, it is marked where it is shown,
    # and a flag on it is never discounted as an isolated injection (its neighbour before it is
    # the instrument before the work).
    shown = [x for x in [out["before"], out["after"]] + list(out["during"]) if x]
    for ev in mrec["events"]:
        q = ev["first_qc_after"]
        if not q:
            continue
        for x in shown:                         # marked wherever it is shown
            if x["stan_id"] == q["stan_id"]:
                x["after_maintenance"] = ev["label"]
        if not ev["files_checked_by_first_qc"]:
            continue                            # other QC runs came between it and the files
        q["counts_as"] = q["grade"]
        for x in shown:
            if x["stan_id"] == q["stan_id"]:
                x["counts_as"] = x["grade"]
        if not any(x["stan_id"] == q["stan_id"] for x in near):
            near.append(q)
    for x in near:
        if x["grade"] != "ok":
            x["samples_between_normal_qc"] = samples_around(x, order, graded_run, files)
    out["gaps"] = qc_gaps(out, spans, mrec, now)
    if near:
        worst = max((x["counts_as"] for x in near), key=GRADE_RANK.get)
        out["verdict"] = {"ok": "good", "check": "check", "concern": "concern"}[worst]
    else:
        out["verdict"] = "check"
    if (mrec["unchecked_files"] or mrec["files_during_downtime"]) and out["verdict"] == "good":
        out["verdict"] = "check"
    out["summary_lines"] = instrument_summary(out)
    out["summary"] = " ".join(out["summary_lines"])
    return out


def _run_line(x, where):
    bits = [f"IPS {x['ips']} ({x['ips_band']})" if x["ips"] is not None else "no IPS",
            f"{x['ids']:,} {x['ids_unit']}" if x["ids"] else f"0 {x['ids_unit']}"]
    if x["ids_vs_median"] is not None:
        bits[-1] += f" ({x['ids_vs_median']:.0%} of the 30-day median)"
    problems = [f["text"] for f in x["flags"] if f["severity"] in ("check", "concern")
                and f["code"] not in ("ids_low",)]
    graded = {"ok": "", "check": " CHECK.", "concern": " CONCERN."}[x["grade"]]
    if x.get("isolated") and x["counts_as"] != x["grade"]:
        graded += (" Isolated: the QC runs just before and after it were normal, so it counts "
                   "as " + {"check": "a check", "ok": "normal"}[x["counts_as"]] + ".")
    elif x.get("isolated"):
        graded += (" The QC runs just before and after it were normal, but as the first QC run "
                   "after maintenance it counts in full.")
    sa = x.get("samples_between_normal_qc")
    if sa:
        n = sa["n_files"]
        graded += (f" {n} of this project's files {'was' if n == 1 else 'were'} acquired between "
                   f"the normal QC runs around it ({sa['from']} to {sa['to']})."
                   if n else " None of this project's files were acquired between the normal "
                   "QC runs around it.")
    if x.get("after_maintenance") and "after the" not in where:
        graded += f" It is the first QC run after the {x['after_maintenance']}."
    return (f"{where} {x['local_time']}: {', '.join(bits)}"
            + (f"; {'; '.join(problems)}" if problems else "") + "." + graded)


def _files(n):
    return f"{n} of this project's file{'s' if n != 1 else ''} {'was' if n == 1 else 'were'}"


def maintenance_lines(m):
    """The maintenance log's lines of the summary: what it holds around the project, the first
    QC run after each event that project files ran after, and files nothing checked."""
    if not m:
        return []
    if m["log"] != "read":
        return [f"STAN's maintenance log could not be read ({m['error']}): maintenance around "
                f"this project is not shown."]
    if not m["n_logged"]:
        return ["STAN has no maintenance events logged for this instrument."]
    if not m["events"]:
        return [f"STAN's maintenance log has nothing for this instrument within {BASELINE_DAYS} "
                f"days of the project."]
    L = [f"Maintenance in STAN's log within {BASELINE_DAYS} days of the project: "
         + "; ".join(f"{e['label']} {e['when']}" for e in m["events"]) + "."]
    for e in m["events"]:
        q, n = e["first_qc_after"], e["files_checked_by_first_qc"]
        what = (f"First QC run on or after the day of the {e['label']} of {e['when']} (it may "
                f"have run before the work: STAN records only the date)" if e["date_only"] else
                f"First QC run after the {e['label']} of {e['when']}")
        if q and n:     # the post-maintenance check for these files (no other QC in between)
            L.append(_run_line(q, f"{what}, the QC check for {n} file{'s' if n != 1 else ''} "
                                  f"acquired after it with no other QC run in between:"))
        elif not q and e["files_after"]:
            L.append(f"No QC run after the {e['label']} of {e['when']} yet, and "
                     f"{_files(e['files_after'])} acquired after it.")
        if e["files_unchecked_after"] and q:
            L.append(f"{_files(len(e['files_unchecked_after']))} acquired after the {e['label']} "
                     f"and before the first QC run after it. CHECK.")
        if e["files_same_day"]:
            L.append(f"{_files(len(e['files_same_day']))} acquired on the day of the "
                     f"{e['label']} ({e['when'].replace(' (date only)', '')}); STAN records only "
                     f"the date, so whether they ran before or after it is not known.")
        if e["files_during_downtime"]:
            L.append(f"{_files(len(e['files_during_downtime']))} acquired while STAN's log has "
                     f"the instrument down ({e['when']}). CHECK.")
    return L


def instrument_summary(o):
    """The instrument's lines of the plain-English summary, one fact per line."""
    w = o["window"]
    span = (f"{w['first_local'][:10]}" if w["first_local"][:10] == w["last_local"][:10]
            else f"{w['first_local'][:10]} to {w['last_local'][:10]}")
    L = [f"{o['instrument']} ({w['n_files']} file{'s' if w['n_files'] != 1 else ''}, acquired "
         f"{span}): {LABEL[o['verdict']]}."]
    b, a = o["before"], o["after"]
    if b:
        L.append(_run_line(b, f"Nearest QC before ({_gap_text(b['days_from_project'])} before "
                              f"the first sample):"))
    else:
        L.append("No QC run on record before the first sample.")
    if o["during"]:
        n = len(o["during"])
        bad = [x for x in o["during"] if x["grade"] != "ok"]
        L.append(f"{n} QC run{'s' if n != 1 else ''} during the project"
                 + (", all within normal range." if not bad else f"; {len(bad)} flagged:"))
        L += [_run_line(x, "During the project,") for x in bad]
    else:
        L.append("No QC run during the project.")
    if a:
        L.append(_run_line(a, f"Nearest QC after ({_gap_text(a['days_from_project'])} after "
                              f"the last sample):"))
    else:
        L.append(f"Nearest QC after: {o['after_note']}.")
    if o["desert"]:
        L.append(f"No QC run within {NEAR_DAYS} days of the project, so the instrument's state "
                 f"while these samples ran was not measured.")
    L += [g["text"] for g in o.get("gaps") or []]
    L += maintenance_lines(o.get("maintenance"))
    bw = o["baseline_window"]
    if bw.get("median_ips") is not None:
        L.append(f"Over the {BASELINE_DAYS} days before, this instrument's {bw['n_runs']} QC runs "
                 f"had a median IPS of {bw['median_ips']:g} ({bw['median_ips_band']}).")
    return L


# The ONE sentence a client deliverable says about this check (make_methods.py prints it). It is
# the same whatever the verdict -- the verdict is staff-only (Brett, 2026-10-01) -- and it states
# the Core's PRACTICE, which is true whatever this project's QC looked like: "acquired regularly
# on each instrument used" was false for an instrument with a 27-day QC gap, or none on record
# (2.10 review). It is in EVERY record (a check that could not run, an instrument with no QC on
# record included): a sentence that came and went with the check's outcome would itself tell
# the client something (2.10 safety review).
METHODS_SENTENCE = ("The Core monitors instrument performance with routine HeLa digest "
                    "quality-control standards; QC records are kept by the UC Davis Proteomics "
                    "Core.")


def methods_sentence(_insts=None):
    """METHODS_SENTENCE, whatever the instruments and their verdicts."""
    return METHODS_SENTENCE


# --------------------------------------------------------- the staff record and gate --
# STAFF-ONLY (Brett, 2026-10-01): the verdict, the section and the record never reach a client
# deliverable. They live in the session's logs/ -- Core-internal already: deliver never ships it,
# the session zip leaves these files out -- and in the run registry (record_run.py), which only
# the Core group reads. The client's Methods get METHODS_SENTENCE and nothing else.
REPORT_TITLE = "QC runs around this project"
RECORD_NAME, STAFF_MD_NAME, ACK_NAME = "qc_bracket.json", "qc_bracket.md", "qc_bracket_ack.json"
GATE_NAME = "qc_gate.json"        # deliver's snapshot of the gate it applied (write_gate_snapshot)
# EVERY staff-only output starts with this marker -- the record's and the acknowledgements' first
# JSON key, the staff page's first line, deliver's gate snapshot -- and the session zip and
# deliver leave out ANY file carrying it, whatever it is called: `--out` can name a record
# anything, anywhere (2.10 safety review: `--out qc.json` escaped the exact-name rule).
STAFF_ONLY_MARKER = "qc_bracket:staff-only"
MARKER_PROBE_BYTES = 4096
# ... and by NAME, for records with no marker (2.10.0 wrote none): the record, its page, the
# acknowledgements and deliver's snapshot -- never scripts/qc_bracket.py itself.
STAFF_ONLY_NAMES = re.compile(r"^(?:qc_bracket(?:_ack)?\.(?:json|md)|qc_gate\.json)$")
FILE_MODE = 0o640                 # as logs/decisions.md: the Core's group reads it on HIVE
MEANING = {
    "good": "The QC runs nearest to and during this project were within this instrument's "
            "normal range.",
    "check": "Something about the QC around this project needs a look (below). It does not by "
             "itself mean the results are wrong.",
    "concern": "A QC run near or during this project performed well below this instrument's "
               "normal. The samples acquired around it should be reviewed before the results "
               "are relied on.",
    "no_qc_on_record": "STAN has no QC run for this instrument, so its performance during the "
                       "project could not be checked from QC.",
}
MAX_TABLE_DURING = 20
# Verdicts that hold the client delivery until a staff member records an acknowledgement.
NEEDS_ACK = ("check", "concern")
ACK_SCHEMA = "qc_bracket_ack"
NOTE_MAX = 300
# Who may acknowledge: staff.require() -- a HIVE login on the Core's staff list, trusted only when
# a Core admin owns it and its folder; with no such list, agent-like names are refused and the
# acknowledgement says the list is not set up. This stops an ACCIDENTAL self-acknowledgement by an
# agent; it is not proof a person did it -- that comes with the board's Duo approval.


def is_staff_only_file(path):
    """True for a staff-only QC file: one of its names (STAFF_ONLY_NAMES), or any file carrying
    STAFF_ONLY_MARKER near its start. The one test session.py's zip and core_submission.py deliver
    use to keep the QC verdict out of what a client receives."""
    if STAFF_ONLY_NAMES.match(os.path.basename(str(path))):
        return True
    try:
        with open(path, "rb") as fh:
            return STAFF_ONLY_MARKER.encode("ascii") in fh.read(MARKER_PROBE_BYTES)
    except OSError:
        return False


def public_gate(gate):
    """What deliver may write where others can read it (<session>/delivery.json): whether the
    report went out, under which record, and when it was acknowledged -- never the verdict, the
    reason, who acknowledged or their note."""
    g = gate or {}
    return {"proceed": g.get("proceed"), "record_sha256": g.get("record_sha256"),
            "ack_at": (g.get("ack") or {}).get("at")}


def write_gate_snapshot(session_dir, gate, now=None):
    """The gate deliver applied, in full, staff-only: <session>/logs/qc_gate.json (0640)."""
    path = os.path.join(record_paths(session_dir)["dir"], GATE_NAME)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    _write_internal(path, json.dumps({"staff_only_marker": STAFF_ONLY_MARKER,
                                      "written_at": acq_time.iso_utc(now or datetime.now(timezone.utc)),
                                      "gate": gate}, indent=2) + "\n")
    return path


def record_paths(session_dir):
    """Where a session's QC record, its staff Markdown and the acknowledgements live."""
    logs = os.path.join(os.path.abspath(session_dir), "logs")
    return {"dir": logs, "json": os.path.join(logs, RECORD_NAME),
            "md": os.path.join(logs, STAFF_MD_NAME), "ack": os.path.join(logs, ACK_NAME),
            "gate": os.path.join(logs, GATE_NAME)}


def _cell_run(x, where):
    if x.get("after_maintenance"):
        where += f"; first QC after the {x['after_maintenance']}"
    ppm = " / ".join(f"{v:.2f}" if v is not None else "-" for v in (x["ms1_ppm"], x["ms2_ppm"]))
    ids = (f"{x['ids']:,} {x['ids_unit']}" if x["ids"] else f"0 {x['ids_unit']}") + (
        f" ({x['ids_vs_median']:.0%})" if x["ids_vs_median"] is not None else "")
    grade = {"ok": "normal", "check": "**check**", "concern": "**concern**"}[x["grade"]]
    if x.get("isolated"):
        grade += " (isolated)"
    return (f"| {where} | {x['local_time']} | "
            f"{x['ips'] if x['ips'] is not None else '-'}"
            f"{' (' + x['ips_band'] + ')' if x['ips_band'] else ''} | {ids} | {ppm} | "
            f"{str(x['fwhm_s']) + ' s' if x['fwhm_s'] else '-'} | {grade} |")


def staff_markdown(rec, record_path=None):
    """The staff page (logs/qc_bracket.md, and the run registry's section): the verdict, the QC
    runs around the project, the maintenance log, and what could not be checked. Never part of
    anything delivered to a client."""
    where = f" (`{os.path.basename(record_path)}`)" if record_path else ""
    L = [f"<!-- {STAFF_ONLY_MARKER} -->", f"# {REPORT_TITLE} (staff only -- never delivered)", ""]
    if not isinstance(rec, dict):
        return "\n".join(L + [f"**Not checked.** The QC record{where} could not be read. Re-run "
                               f"`scripts/qc_bracket.py`."])
    if rec.get("status") != "ok":
        return "\n".join(L + [f"**Not checked.** {rec.get('summary') or 'The check did not complete.'}",
                               "", "This is not a pass: the QC runs around this project were not "
                               "looked at. Re-run `scripts/qc_bracket.py` once the reason above "
                               "is fixed."])
    L += [f"**Verdict: {LABEL[rec['verdict']]}.** {MEANING[rec['verdict']]}", "",
          "The Core runs a QC standard on each instrument between projects, and STAN scores "
          "every one. These are the QC runs nearest to this project's samples, each compared "
          f"with the same instrument's QC runs over the {BASELINE_DAYS} days before the project. "
          f"IPS is STAN's Instrument Performance Score (green 80+, amber 60-79, red under 60); "
          f"it is shown but not graded. IDs: % of that 30-day median. Times are Pacific "
          f"(the instruments' clock).", ""]
    for o in rec.get("instruments") or []:
        L += [f"## {o['instrument']}: {LABEL[o['verdict']]}", ""]
        if o["verdict"] == "no_qc_on_record":
            L += [" ".join(o["summary_lines"]), ""]
            continue
        w = o["window"]
        L += [f"{w['n_files']} file{'s' if w['n_files'] != 1 else ''} acquired "
              f"{w['first_local']} to {w['last_local']}.", "",
              "| QC run | When | IPS | IDs | MS1 / MS2 ppm | FWHM | Grade |",
              "|---|---|---|---|---|---|---|"]
        if o["before"]:
            L.append(_cell_run(o["before"], f"Nearest before ({_gap_text(o['before']['days_from_project'])})"))
        else:
            L.append("| Nearest before | none on record | | | | | |")
        shown = sorted(o["during"], key=lambda x: (x["grade"] == "ok", x["run_date"]))
        for x in sorted(shown[:MAX_TABLE_DURING], key=lambda x: x["run_date"]):
            L.append(_cell_run(x, "During"))
        if len(o["during"]) > MAX_TABLE_DURING:
            L.append(f"| During | {len(o['during']) - MAX_TABLE_DURING} more, all normal | | | | | |")
        if not o["during"]:
            L.append("| During | **none** | | | | | |")
        if o["after"]:
            L.append(_cell_run(o["after"], f"Nearest after ({_gap_text(o['after']['days_from_project'])})"))
        else:
            L.append(f"| Nearest after | {o.get('after_note', 'none yet')} | | | | | |")
        shown_ids = {x["stan_id"] for x in [o["before"], o["after"]] + list(o["during"]) if x}
        for ev in (o.get("maintenance") or {}).get("events") or []:
            q = ev["first_qc_after"]
            if q and ev["files_checked_by_first_qc"] and q["stan_id"] not in shown_ids:
                L.append(_cell_run(dict(q, after_maintenance=None),
                                   f"First QC after the {ev['label']}"))
                shown_ids.add(q["stan_id"])
        L.append("")
        notes = [ln for ln in o["summary_lines"][1:]
                 if not ln.startswith(("Nearest QC before", "Nearest QC after"))
                 or "CHECK" in ln or "CONCERN" in ln]
        L += [f"- {ln}" for ln in notes] + [""]
    if rec.get("not_checked"):
        L += ["**Not checked** (no acquisition time or instrument could be read):", ""]
        L += [f"- `{os.path.basename(f['file'])}`: {f['why']}" for f in rec["not_checked"][:20]]
        if len(rec["not_checked"]) > 20:
            L.append(f"- and {len(rec['not_checked']) - 20} more (see `{RECORD_NAME}`)")
        L.append("")
    L += ["**How this was checked:** " + " ".join(
        a if a.endswith(".") else a + "." for a in rec.get("assumptions") or [])
        + f" QC runs from {rec['stan']['url']} (checked {rec['checked_at']}); the full record is "
          f"`{RECORD_NAME}`."]
    return "\n".join(L)


def _sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_acks(path):
    """The acknowledgements on record ([] when none), or raise ValueError for a damaged file --
    a damaged acknowledgement never counts as one."""
    if not os.path.isfile(path):
        return []
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    if not (isinstance(d, dict) and d.get("schema") == ACK_SCHEMA and isinstance(d.get("acks"), list)):
        raise ValueError(f"{path} is not a {ACK_SCHEMA} record")
    return [x for x in d["acks"] if isinstance(x, dict)]


def covers_other_files(rec, session_dir):
    """Why `rec` does not cover the session's current raw-file list (None when it does): files
    added, removed or re-pointed since the check ran. Re-run step 8e."""
    have = (rec.get("files") or {}).get("list_sha256")
    try:
        now, _src = session_files(session_dir)
    except Exception as e:                    # never a pass by accident
        return f"the session's raw-file list could not be read ({type(e).__name__}: {e})"
    if not have:
        return ("the QC record does not say which files it checked (an earlier version of this "
                "check made it): re-run step 8e (qc_bracket.py --session)")
    if have != list_sha256(now):
        return (f"the QC record covers a different file list from the session's ("
                f"{(rec.get('files') or {}).get('n')} file(s) checked, {len(now)} in the "
                f"session now): re-run step 8e (qc_bracket.py --session)")
    return None


def delivery_gate(session_dir):
    """May this session's client report go out? THE one rule (core_submission.py deliver and
    record_run.py both ask it here):
      * a record for another file list than the session's (files added or removed since the
        check) is refused outright: re-run step 8e -- no acknowledgement releases it;
      * a verdict of check or concern on any instrument, any file that could not be placed in
        time (not_checked) -- or no record, an unreadable one, or a check that could not run
        (STAN unreachable, nothing placeable) -- needs a staff acknowledgement recorded for THIS
        record (its sha256: a re-run needs a new one);
      * good and no_qc_on_record proceed; instruments with no QC on record are listed.
    -> {proceed, needs_ack, needs_rerun, reason, status, verdict, record, record_sha256, ack,
        no_qc_on_record, ack_file}."""
    p = record_paths(session_dir)
    out = {"proceed": True, "needs_ack": False, "needs_rerun": False, "reason": None,
           "status": None, "verdict": None, "record": p["json"], "record_sha256": None,
           "ack": None, "no_qc_on_record": [], "ack_file": p["ack"]}
    if not os.path.isfile(p["json"]):
        reason = (f"the QC runs around this project were not checked (no logs/{RECORD_NAME}): run "
                  f"qc_bracket.py --session first")
    else:
        out["record_sha256"] = _sha256(p["json"])
        try:
            with open(p["json"], encoding="utf-8") as fh:
                rec = json.load(fh)
        except (OSError, ValueError):
            rec = None
        readable = isinstance(rec, dict) and str(rec.get("schema", "")).startswith("qc_bracket/")
        other = covers_other_files(rec, session_dir) if readable else None
        if not readable:
            reason = f"the QC record logs/{RECORD_NAME} could not be read"
        elif other:
            # a record for another file list: re-run, never acknowledged away
            out.update(status=rec.get("status"), verdict=rec.get("verdict"), needs_rerun=True,
                       proceed=False, reason=other)
            return out
        elif rec.get("status") != "ok":
            out["status"] = rec.get("status")
            reason = (f"the QC runs around this project could not be checked "
                      f"({rec.get('status')}: {(rec.get('stan') or {}).get('error') or 'see the record'})")
        else:
            insts = rec.get("instruments") or []
            out.update(status="ok", verdict=rec.get("verdict"),
                       no_qc_on_record=[o.get("instrument") for o in insts
                                        if o.get("verdict") == "no_qc_on_record"])
            flagged = [o for o in insts if o.get("verdict") in NEEDS_ACK]
            unplaced = rec.get("not_checked") or []
            if not flagged and not unplaced:
                return out
            reason = "; ".join(
                (["QC verdict " + "; ".join(f"{o.get('instrument')}: {LABEL[o['verdict']]}"
                                            for o in flagged)] if flagged else [])
                + ([f"{len(unplaced)} file{'s' if len(unplaced) != 1 else ''} could not be "
                    f"placed in time, so the QC around {'them' if len(unplaced) != 1 else 'it'} "
                    f"was not checked"] if unplaced else [])) + (
                f" -- a staff member must review logs/{STAFF_MD_NAME} and record an "
                f"acknowledgement")
    out.update(needs_ack=True, reason=reason, proceed=False)
    try:
        acks = _read_acks(p["ack"])
    except (OSError, ValueError) as e:
        out["reason"] += f" (the acknowledgement file is unreadable: {e})"
        return out
    match = [x for x in acks if x.get("record_sha256") == out["record_sha256"]]
    if match:
        out.update(ack=match[-1], proceed=True)
    return out


def acknowledge(session_dir, by, note, now=None, no_record_reason=None):
    """Record a staff acknowledgement of this session's QC record as it is now. -> the gate.
    Refused: an acknowledger who is not staff (staff.require), a record for another file
    list (re-run), and no record at all unless `no_record_reason` says why the check could not
    run. An acknowledgement applies to the record it was made for and no other: one made with
    no record never releases a record written later."""
    note = (note or "").strip()
    c = staff.require(by)               # raises ValueError; prints its own NOT VERIFIED warning
    by, staff_list = c["by"], c["checked"]
    if not note or "\n" in note or "\r" in note or len(note) > NOTE_MAX:
        raise ValueError(f"--note: one line, 1-{NOTE_MAX} characters, saying what was reviewed "
                         f"and decided")
    gate = delivery_gate(session_dir)
    if gate["needs_rerun"]:
        raise ValueError(f"nothing to acknowledge yet: {gate['reason']}")
    no_record_reason = (no_record_reason or "").strip() or None
    if gate["record_sha256"] is None:
        if not no_record_reason or "\n" in no_record_reason:
            raise ValueError("there is no QC record: run step 8e (qc_bracket.py --session); only "
                             "when it cannot run, say why in one line with --no-record-reason")
    elif no_record_reason:
        raise ValueError("--no-record-reason is for a session with no QC record; this one has one")
    p = record_paths(session_dir)
    acks = _read_acks(p["ack"])
    import getpass
    try:
        account = getpass.getuser()
    except Exception:
        account = None
    acks.append({"by": by, "account": account, "staff_list": staff_list,
                 "no_record_reason": no_record_reason,
                 "at": acq_time.iso_utc(now or datetime.now(timezone.utc)), "note": note,
                 "record_sha256": gate["record_sha256"], "status": gate["status"],
                 "verdict": gate["verdict"], "reason_acknowledged": gate["reason"]})
    os.makedirs(os.path.dirname(p["ack"]), exist_ok=True)
    _write_internal(p["ack"], json.dumps({"staff_only_marker": STAFF_ONLY_MARKER,
                                          "schema": ACK_SCHEMA, "schema_version": 1,
                                          "acks": acks}, indent=2) + "\n")
    return delivery_gate(session_dir)


def _write_internal(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    try:
        os.chmod(path, FILE_MODE)
    except OSError:
        pass


# --------------------------------------------------------------------------- the run --
def overall_summary(rec):
    if rec["status"] == "stan_unreachable":
        return (f"QC around this project could NOT be checked: {rec['stan']['error']}. This is "
                f"not a pass -- re-run qc_bracket.py when STAN is reachable.")
    if rec["status"] == "nothing_checked":
        return ("QC around this project could NOT be checked: no file gave both an instrument "
                "and an acquisition time (see not_checked).")
    L = [f"QC runs around this project: {LABEL[rec['verdict']]}."]
    for o in rec["instruments"]:
        L += o["summary_lines"]
    if rec["not_checked"]:
        n = len(rec["not_checked"])
        L.append(f"{n} file{'s' if n != 1 else ''} could not be placed in time and "
                 f"{'were' if n != 1 else 'was'} not checked (see not_checked).")
    return "\n".join(L)


def session_files(session_dir):
    """(the session's raw file paths, normalised, and where the list came from). The one reader
    the check and the delivery gate share: session_docs.raw_record."""
    from session import paths_for
    from session_docs import raw_record           # the one reader of a session's raw list
    raws, src = raw_record(paths_for(session_dir))
    return [_norm(r) for r in raws], src


def list_sha256(files):
    """Which files a record covers: the sha256 of the sorted, normalised paths, one per line."""
    import hashlib
    return hashlib.sha256("\n".join(sorted(_norm(f) for f in files)).encode("utf-8")).hexdigest()


def collect_files(a):
    """(file paths, where the list came from, the step-2 JSON path or None)."""
    if a.session:
        from session import paths_for
        raws, src = session_files(a.session)
        acq = a.acquisition_json or next(
            (c for c in (os.path.join(paths_for(a.session)["input_dir"], "acquisition.json"),)
             if os.path.isfile(c)), None)
        return raws, src, acq
    files = []
    for pat in a.files or []:
        files.extend(sorted(glob.glob(pat)) or [pat])
    return [_norm(f) for f in files], "--files", a.acquisition_json


def _say(msg):
    print(f"[qc_bracket] {msg}", file=sys.stderr, flush=True)


def write(rec, out):
    """The record to stdout and `out`, with the staff page beside it (<out stem>.md), both
    readable by the owner and the Core's group only."""
    text = json.dumps(rec, indent=2)
    if out:
        os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
        _write_internal(out, text + "\n")
        _write_internal(os.path.splitext(out)[0] + ".md", staff_markdown(rec, out) + "\n")
    print(text)


def gate_main(argv, now=None):
    """`qc_bracket.py gate --session S`: the delivery gate as JSON; exit 0 = may proceed, 2 = a
    staff acknowledgement is needed. `qc_bracket.py ack --session S --by NAME --note "..."`:
    record one, for the record as it is now."""
    cmd = argv[0]
    ap = argparse.ArgumentParser(prog=f"qc_bracket.py {cmd}")
    ap.add_argument("--session", required=True)
    if cmd == "ack":
        ap.add_argument("--by", required=True, help="the HIVE login of the staff member "
                                                    "acknowledging (on the Core staff list)")
        ap.add_argument("--note", required=True, help="one line: what was reviewed and decided")
        ap.add_argument("--no-record-reason", help="only when there is no QC record: one line "
                                                   "saying why the check could not run")
    a = ap.parse_args(argv[1:])
    if cmd == "ack":
        try:
            gate = acknowledge(a.session, a.by, a.note, now, a.no_record_reason)
        except (ValueError, OSError) as e:
            _say(f"ERROR: {e}")
            return EXIT_USAGE
        _say(f"acknowledged by {a.by}: {gate['reason'] or 'nothing needed it'}")
    else:
        gate = delivery_gate(a.session)
        _say("the client report may go out" if gate["proceed"] else
             f"HOLD: {gate['reason']} (qc_bracket.py ack --session ... --by ... --note ...)")
    print(json.dumps(gate, indent=2))
    return EXIT_OK if gate["proceed"] else EXIT_USAGE


def main(argv=None, now=None):
    """`now` is for the tests: "none yet" and the open gap are measured to it."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in ("ack", "gate"):
        return gate_main(argv, now)
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--session", help="session dir: its raw file list, input/acquisition.json, "
                                       "and logs/qc_bracket.json (+ .md) as --out")
    src.add_argument("--files", nargs="+", help="raw files (.d / .raw) or globs")
    ap.add_argument("--acquisition-json", help="detect_acquisition.py's output for these files "
                                               "(step 2): their instrument and time, not re-read")
    ap.add_argument("--out", help="where to write the record, and the staff page beside it "
                                  "(default with --session: <session>/logs/qc_bracket.json; "
                                  "else stdout only)")
    ap.add_argument("--stan-url", default=None,
                    help=f"STAN dashboard (default ${STAN_URL_ENV}, else {STAN_URL})")
    ap.add_argument("--timeout", type=float, default=TIMEOUT_S, help="seconds per STAN request")
    ap.add_argument("--allow-login-node", action="store_true",
                    help="read more than a handful of .raw here even on a cluster login node")
    a = ap.parse_args(argv)
    now = now or datetime.now(timezone.utc)
    base = (a.stan_url or os.environ.get(STAN_URL_ENV) or STAN_URL).rstrip("/")
    out_path = a.out or (record_paths(a.session)["json"] if a.session else None)

    files, files_from, acq_path = collect_files(a)
    if not files:
        _say(f"no raw files to check ({'no raw file list in the session' if a.session else 'none given'})")
        return EXIT_USAGE
    known, acq_problem = load_acquisition_json(acq_path)
    to_read = [f for f in files if f not in known]
    n_raw = sum(1 for f in to_read if f.lower().endswith(".raw") and os.path.isfile(f))
    if n_raw:
        import detect_acquisition as da
        if n_raw > da.LOGIN_NODE_MAX_RAW and da.on_cluster_login_node() and not a.allow_login_node:
            _say(f"REFUSING to read {n_raw} Thermo .raw on a cluster login node (each is a "
                 f"ThermoRawFileParser process, ~3-7 s): pass --acquisition-json with step 2's "
                 f"output, or run this under srun, or add --allow-login-node.")
            return EXIT_USAGE
    if acq_problem:
        _say(f"WARNING: {acq_problem}")
    ident = [identify(f, known.get(f)) for f in files]
    for f in ident:
        # a time this cannot place on the clock is no time: never guessed into one
        t = acq_time.parse_iso(f["acquired_at"]) if f["acquired_at"] else None
        if f["acquired_at"] and (t is None or t.tzinfo is None):
            f["acquired_at_note"] = (f"acquired_at {f['acquired_at']!r} is not an ISO 8601 time "
                                     f"with a time zone")
            f["acquired_at"] = None

    not_checked = [{"file": f["file"], "why": (f["acquired_at_note"] or "no acquisition time")
                    if not f["acquired_at"] else "no instrument name could be read"}
                   for f in ident if not (f["acquired_at"] and f["instrument"])]
    by_inst = {}
    for f in ident:
        if f["acquired_at"] and f["instrument"]:
            by_inst.setdefault(f["instrument"].strip(), []).append(f)
    notes = sorted({f["acquired_at_note"] for f in ident
                    if f["acquired_at"] and f["acquired_at_note"]})
    rec = {"staff_only_marker": STAFF_ONLY_MARKER,         # first: is_staff_only_file()
           "schema": SCHEMA, "schema_version": SCHEMA_VERSION, "tool": "qc_bracket.py",
           "staff_only": True, "checked_at": acq_time.iso_utc(now),
           "status": None, "verdict": None, "verdict_label": None, "summary": None,
           "stan": {"url": base, "endpoint": RUNS_PATH, "qc_only": True,
                    "hidden_runs": "left out, as STAN's dashboard leaves them out",
                    "pages": 0, "error": None},
           "files": {"n": len(files), "from": files_from, "acquisition_json": acq_path,
                     # which files this record covers: the delivery gate compares it with the
                     # session's raw list (a GOOD for 3 files must not release 4)
                     "list_sha256": list_sha256(files),
                     "acquisition_json_problem": acq_problem,
                     "n_checked": sum(len(v) for v in by_inst.values())},
           "rules": {"baseline_days": BASELINE_DAYS, "near_days": NEAR_DAYS,
                     "min_baseline_runs": MIN_BASELINE,
                     "ips_bands": {"green": f">= {IPS_GREEN}", "amber": f">= {IPS_AMBER}",
                                   "red": f"< {IPS_AMBER}"},
                     "concern": ["no identifications", f"IDs < {IDS_CONCERN:.0%} of the 30-day "
                                 f"median"],
                     "check": [f"IDs < {IDS_CHECK:.0%} of the 30-day median",
                               f"MS1/MS2 mass error > {PPM_RATIO}x the median and > median + "
                               f"{PPM_FLOOR} ppm", f"FWHM > {FWHM_RATIO}x the median",
                               f"no QC run within {NEAR_DAYS} days of the project"],
                     "isolated": f"a flagged run whose neighbouring QC runs in STAN are both "
                                 f"normal and within {ISOLATION_DAYS} days of it counts one "
                                 f"level lower",
                     "no_baseline": f"fewer than {MIN_BASELINE} same-mode QC runs in the "
                                    f"{BASELINE_DAYS} days before: check",
                     "ips": "shown per run and as the 30-day median; not graded",
                     "maintenance": f"STAN's maintenance log ({EVENTS_PATH}) within "
                                    f"{BASELINE_DAYS} days either side: the first QC run after "
                                    f"an event that project files ran after speaks for them "
                                    f"and is never discounted as isolated; files acquired after "
                                    f"an event and before any QC, or during a logged downtime, "
                                    f"are a check; a gap of more than {GAP_DAYS} days without "
                                    f"QC is explained from the log or said to have no record",
                     "not_used": "gate_result (STAN records 'pass' for every run)"},
           "assumptions": [
               "Files are matched to STAN's QC runs by instrument MODEL: STAN records no serial "
               "number. Valid for runs acquired on the UC Davis Proteomics Core's instruments.",
               "STAN does not record whether a QC run's time came from the file's header or its "
               "modification time (the end of the run); the project's files are timed from the "
               "start of acquisition."] + notes,
           "instruments": [], "not_checked": not_checked, "methods_sentence": methods_sentence()}

    if not by_inst:
        rec.update(status="nothing_checked")
        rec["summary"] = overall_summary(rec)
        for ln in rec["summary"].splitlines():
            _say(ln)
        write(rec, out_path)
        return EXIT_NOTHING

    _say(f"asking {base}{RUNS_PATH} for {', '.join(sorted(by_inst))}")
    fetched = {}
    try:
        for inst in sorted(by_inst):
            oldest = min(acq_time.parse_iso(f["acquired_at"]) for f in by_inst[inst])
            rows, pages = stan_runs(base, inst, oldest - timedelta(days=BASELINE_DAYS),
                                    a.timeout, _say)
            fetched[inst] = rows
            rec["stan"]["pages"] += pages
    except StanUnreachable as e:
        rec.update(status="stan_unreachable")
        rec["stan"]["error"] = str(e)
        rec["summary"] = overall_summary(rec)
        for ln in rec["summary"].splitlines():
            _say(ln)
        write(rec, out_path)
        return EXIT_STAN

    logs = {}
    for inst in sorted(by_inst):
        logs[inst] = stan_events(base, inst, a.timeout)
        if logs[inst][1]:
            _say(f"WARNING: {inst}: {logs[inst][1]} -- the QC runs still speak; the gaps are "
                 f"not explained")
    rec["stan"]["events_endpoint"] = EVENTS_PATH
    rec["instruments"] = [bracket(inst, by_inst[inst], fetched[inst], now, *logs[inst])
                          for inst in sorted(by_inst)]
    verdict = max((o["verdict"] for o in rec["instruments"]), key=RANK.get)
    rec.update(status="ok", verdict=verdict, verdict_label=LABEL[verdict],
               methods_sentence=methods_sentence(rec["instruments"]))
    rec["summary"] = overall_summary(rec)
    for ln in rec["summary"].splitlines():
        _say(ln)
    write(rec, out_path)
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
