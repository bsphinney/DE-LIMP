#!/usr/bin/env python3
"""
core_submission.py  --  Take a CoreOmics submission (PROT_0807) from "search the data" to a
delivered, shareable result.

A Core service run does not start as a folder of raw files. It starts as a CoreOmics record
-- a PI, a sample sheet with conditions, an organism typed into a form -- and it ends as a
Bioshare folder the collaborator opens. Every step in between is a place where the work can
silently attach to the wrong thing:

  * The sample sheet's `unique_id`s are the ONLY link to the raw files, and they are short,
    reused by other submissions (BN1-6 in both PROT_0794 and PROT_0776), and collide with
    plate wells (A3, H10) that match thousands of files.
  * The service directory is organized by human-named folders (`McDonald karen`,
    `UCSF/Feeley_lab`), so "where does this project live" is a lookup, not a formula.
  * A share built from links depends on server settings. PROT_0793's links into /quobyte
    served nothing over https because Bioshare's Apache file-streaming allowed /quobyte on
    port 80 but not 443 (fixed 2026-09-16); absolute /nfs links are invisible over SMB. So deliverables are
    copied as real files, and the only links this script makes are relative raw/ links.

So this script does the bookkeeping deterministically and hands every judgment call back to
a human as an exit code plus the exact question. The search and DE in between are the
ordinary skill flow -- nothing here runs an engine.

WHERE EACH SUBCOMMAND RUNS
    identify     LOCAL    which submission is this data? PROT/hex ids in the names or the
                          user's message, else the sample ids in the file names (token)
    fetch        LOCAL    the CoreOmics token lives on the staff member's computer
                          (~/.coreomics_token); HIVE has none
    locate       ON HIVE  reads the Flinders raw_data tree
    stage        ON HIVE  symlinks the raw files into the Core service directory
    conditions   either   wraps collect_conditions.py
    deliver      ON HIVE  copies the deliverables into the submission's Bioshare folder
                          (--mode raw-only for "I only require raw data" submissions)
    bioshare     LOCAL    status | ensure | send, through the CoreOmics Bioshare plugin
    email-draft  LOCAL    writes a draft to the submitter; never sends

    python3 core_submission.py fetch 807 --out ~/core/PROT_0807
    bash hive_exec.sh --put ~/core/PROT_0807/submission_summary.json '~/core/PROT_0807/'
    bash hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/core_submission.py locate \\
        --summary ~/core/PROT_0807/submission_summary.json --out ~/core/PROT_0807'

Every subcommand prints ONE JSON object to stdout; human notes go to stderr. `stage`,
`deliver` and `bioshare ensure|send` are DRY RUNS unless --apply is given.

PATHS. Bioshare and CoreOmics know a share only by its HIVE path, so every server-side path
(`share_dir` in the summary, Bioshare's `link_to_path`) is built with posixpath from the HIVE
Flinders root and never re-derived from CORE_FLINDERS_ROOT, which may be an SMB mount, a
Windows drive, or a test directory. CORE_FLINDERS_ROOT only says where to do file work.

CONFIGURATION (environment -- which is also how the tests point it at temp dirs)
    COREOMICS_BASE_URL   default https://ucdavis.coreomics.com/server/api
    COREOMICS_TOKEN      else the contents of ~/.coreomics_token
    CORE_FLINDERS_ROOT   default: the Flinders share's HIVE path in hive_shares.tsv
                         (/nfs/lssc0/flinders/proteomics)  (local file work only)
    CORE_WORK_ROOT       default: the proteomics-grp share's SERVICE tree (hive_shares.tsv)

EXIT CODES -- the orchestrator branches on these:
    0  ok. From `identify`: exactly one candidate -- confirm it with the user in one line
    2  needs a human decision, or a hard gate failed. Proposal files are still written.
       From `deliver --apply`: the delivery is NOT verified -- do not share it.
       From `identify`: none or several candidates -- ask the user for the number
    3  CoreOmics or the filesystem is unreachable, auth failed (server detail printed), or
       the scripts directory on this machine is incomplete
    4  (locate) the files look like an HT plate -- use ht_manifest.py instead
    5  (deliver) refused by the size guard -- sbatch the deliver_job.sh it wrote
    1  an unexpected error, i.e. a bug -- read stderr
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import errno
import glob
import hashlib
import importlib
import json
import os
import posixpath
import re
import shlex
import shutil
import stat
import subprocess
import sys
import traceback
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
try:                        # imported now, while os.path is the real one (tests swap in ntpath)
    import share_map
except ImportError:         # a partial copy of scripts/: server_flinders_root() says so
    share_map = None

DEFAULT_BASE_URL = "https://ucdavis.coreomics.com/server/api"
TOKEN_FILE = "~/.coreomics_token"
LAB = "PROTEOMICS"
HTTP_TIMEOUT = 60
SCHEMA = "core_submission/1"
MARKER = ".core_submission.json"
DELIVERY_MARKER = ".core_delivery.json"
GENERATED_MARK = "<!-- generated by core_submission.py stage -- do not hand-edit; notes go in SUBMISSION.notes.md -->"
SYNC_HINT = ("the scripts directory is incomplete: sync the whole scripts/ directory, and "
             ".claude-plugin/ beside it, to ~/proteomics-pipeline/ (bash scripts/hive_exec.sh "
             "--put-skill)")
RAW_WHITELIST_WARNING = (
    "Check once: Brett reported (2026-09-16) that PROT_0793's raw/ links download from Bioshare, "
    "and relative links resolve to the same files, but nobody has checked a relative raw/ link "
    "in Bioshare itself. Before `bioshare send`, open the link and download one raw file.")

EXIT_OK, EXIT_DECIDE, EXIT_UNREACHABLE, EXIT_HT, EXIT_GUARD = 0, 2, 3, 4, 5

# Label reuse is only a risk between submissions close in time; a year-old BN1 is not a
# plausible owner of this month's file. fetch records neighbours inside this window, and
# locate refuses a --max-days wider than it.
NEIGHBOR_WINDOW_DAYS = 240
NEIGHBOR_MAX_PAGES = 200            # 20,000 submissions -- a runaway guard, not a limit
# An id whose token names more raw files than this, across every instrument and year, is
# too common to identify one submission's runs (a plate position or a generic label).
WEAK_MAX_FILES = 10

RAW_ONLY_MARKER = "only require raw data"   # data_analysis: "I only require raw data and ..."

# CoreOmics `mass_spec_wanted` -> the raw_data instrument folder. Scanned first, never
# exclusively: staff sometimes run a sample on a different instrument than requested.
INSTRUMENT_FOLDERS = {"timstof ht": "tTOF_HT", "exploris 480": "Exploris480",
                      "fusion lumos": "Lumos1"}
RAW_SUFFIXES = (".d", ".raw", ".mzml", ".mzml.gz")


class Stop(Exception):
    """End a subcommand with a specific exit code and the JSON that explains it."""

    def __init__(self, code: int, payload: dict):
        super().__init__(payload.get("error", ""))
        self.code, self.payload = code, payload


class ApiError(Exception):
    def __init__(self, status, detail: str, url: str | None):
        super().__init__(f"HTTP {status}: {detail}" if status else detail)
        self.status, self.detail, self.url = status, detail, url


class UnsafePath(Exception):
    """A write would pass through a symlink or leave the Flinders root."""


def emit(obj) -> None:
    print(json.dumps(obj, indent=2, default=str))


def note(msg: str) -> None:
    print(f"[core_submission] {msg}", file=sys.stderr)


def _s(v) -> str:
    return "" if v is None else str(v).strip()


def now_utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _fold(s) -> str:
    return "".join(c for c in unicodedata.normalize("NFKD", _s(s)) if not unicodedata.combining(c))


def normkey(s) -> str:
    """Letters and digits only, case-folded: `EB_001`, `eb-001` and `EB 001` are one key."""
    return re.sub(r"[^a-z0-9]", "", _fold(s).casefold())


def sibling(module: str):
    """Import a script that must sit next to this one. On HIVE the scripts are a copy; a
    partial copy must say so, not end in a traceback."""
    try:
        return importlib.import_module(module)
    except ImportError as e:
        raise Stop(EXIT_UNREACHABLE, {"error": f"{module}.py could not be imported from {HERE}: {e}",
                                      "hint": SYNC_HINT})


# ---------------------------------------------------------------- configuration --
def base_url() -> str:
    return (os.environ.get("COREOMICS_BASE_URL") or DEFAULT_BASE_URL).rstrip("/")


def server_flinders_root() -> str:
    """The HIVE path of the Flinders share -- the only spelling Bioshare accepts in link_to_path
    (existing shares: /nfs/lssc0/flinders/proteomics/coreomics/projects/2026/08/<id>/share).
    From hive_shares.tsv, never a second copy of it here."""
    sm = _share_map()
    try:
        hive = sm.hive_root(sm.FLINDERS_SHARE)
    except OSError as e:
        raise Stop(EXIT_UNREACHABLE, {"error": f"hive_shares.tsv could not be read: {e}",
                                      "hint": SYNC_HINT})
    if not hive:
        raise Stop(EXIT_UNREACHABLE, {"error": f"hive_shares.tsv has no HIVE path for the "
                                               f"'{sm.FLINDERS_SHARE}' (Flinders) share",
                                      "hint": SYNC_HINT})
    return hive


def _share_map():
    """share_map.py: which share is which HIVE path, and the Flinders trees (one definition)."""
    return share_map or sibling("share_map")


def flinders_root() -> str:
    """Where THIS machine does file work. Never used to build a path sent to a server."""
    return os.path.abspath(os.environ.get("CORE_FLINDERS_ROOT") or server_flinders_root())


def raw_root() -> str:
    return os.path.join(flinders_root(), *_share_map().FLINDERS_RAW)


def service_root() -> str:
    return os.path.join(flinders_root(), *_share_map().FLINDERS_SERVICE)


def server_work_root() -> str:
    """The Core's SERVICE tree on Quobyte: the proteomics-grp share's HIVE path from
    hive_shares.tsv plus share_map.QUOBYTE_SERVICE -- the same tree fran_deposit's backfill walks.
    Never a second copy of the path here."""
    sm = _share_map()
    try:
        hive = sm.hive_root(sm.QUOBYTE_SHARE)
    except OSError as e:
        raise Stop(EXIT_UNREACHABLE, {"error": f"hive_shares.tsv could not be read: {e}",
                                      "hint": SYNC_HINT})
    if not hive:
        raise Stop(EXIT_UNREACHABLE, {"error": f"hive_shares.tsv has no HIVE path for the "
                                               f"'{sm.QUOBYTE_SHARE}' (Quobyte) share",
                                      "hint": SYNC_HINT})
    return "/".join((hive.rstrip("/"),) + sm.QUOBYTE_SERVICE)


def work_root() -> str:
    return os.path.abspath(os.environ.get("CORE_WORK_ROOT") or server_work_root())


def is_under(path: str, root: str) -> bool:
    """True when `path` is `root` or inside it. Compares REAL paths: a link's visibility over
    SMB/Bioshare is decided by where it actually lands, not by how it was spelled."""
    p, r = os.path.realpath(path), os.path.realpath(root)
    return p == r or p.startswith(r.rstrip(os.sep) + os.sep)


def skill_version():
    """skill_version.py's reading, the one reader of .claude-plugin/plugin.json: the version, or
    its "(unknown -- plugin.json not found)" tag when scripts/ went up to HIVE without
    .claude-plugin/ (`hive_exec.sh --put-skill` puts both). Never a guessed version; None only
    when skill_version.py itself is missing from a partial copy of scripts/ (SYNC_HINT)."""
    try:
        import skill_version as sv
    except ImportError:
        return None
    return sv.skill_version()


# ------------------------------------------------------------ ids, dates, campus --
HEX_ID = re.compile(r"^[0-9a-f]{12}$")
SUBMISSION_NUMBER = re.compile(r"^(?:prot)?[\s_\-]*0*(\d{1,5})$", re.I)


def normalize_submission(text) -> tuple:
    """User input -> ("internal_id", "PROT_0807") or ("id", "99922f5337f8").

    Accepts PROT_0807, prot-807, 0807, 807, #807, or a 12-hex CoreOmics id. Look-up is by
    the exact `internal_id` filter, never `search=` -- `search=0807` returns two records."""
    s = _s(text).lstrip("#").strip()
    if HEX_ID.match(s.lower()):
        return "id", s.lower()
    m = SUBMISSION_NUMBER.match(s)
    if m and int(m.group(1)) > 0:
        return "internal_id", "PROT_%04d" % int(m.group(1))
    raise ValueError(f"not a submission number or CoreOmics id: {text!r} "
                     f"(try 807, PROT_0807, or the 12-character id)")


def parse_date(value):
    """The literal YYYY-MM-DD at the start of a timestamp, or None. No timezone conversion:
    `2026-08-31T23:30:00-07:00` is 2026-08-31 -- the same rule the CoreOmics builder uses."""
    m = re.match(r"^\s*(\d{4})-(\d{2})-(\d{2})", _s(value))
    if not m:
        return None
    try:
        return dt.date(int(m.group(1)), int(m.group(2)), int(m.group(3)))
    except ValueError:
        return None


def _year_month(sub_id: str, submitted: str) -> tuple:
    """YYYY and MM as the LITERAL first two `-` fields of `submitted` -- exactly coreomics_fs
    build_views.ensure_canonical: `y, m = submitted.split()[0].split("-")[:2]`. Converting to
    UTC would file a late-evening submission under the wrong month (verified: all 6
    last-evening-of-month submissions sit under their LOCAL month)."""
    parts = _s(submitted).split()
    fields = parts[0].split("-")[:2] if parts else []
    if len(fields) != 2 or not re.match(r"^\d{4}$", fields[0]) or not re.match(r"^\d{2}$", fields[1]):
        raise ValueError(f"cannot take YYYY-MM from submitted={submitted!r}")
    if not HEX_ID.match(_s(sub_id)):
        raise ValueError(f"not a CoreOmics id: {sub_id!r}")
    return fields[0], fields[1]


def server_project_dir(sub_id: str, submitted: str) -> str:
    """The HIVE path of <projects>/<YYYY>/<MM>/<id>, always with forward slashes."""
    y, m = _year_month(sub_id, submitted)
    return posixpath.join(server_flinders_root(), "coreomics", "projects", y, m, sub_id)


def local_project_dir(sub_id: str, submitted: str) -> str:
    """The same directory as seen by this machine's CORE_FLINDERS_ROOT."""
    y, m = _year_month(sub_id, submitted)
    return os.path.join(flinders_root(), "coreomics", "projects", y, m, sub_id)


def classify_campus(rec: dict) -> tuple:
    """(campus, institution, basis). On campus only when the PI's institution is exactly
    "UC Davis" (153 of the 300 newest records); everything else is off campus."""
    inst, basis = None, "missing"
    pi = rec.get("pi")
    if isinstance(pi, dict):
        i = pi.get("institution")
        name = i.get("name") if isinstance(i, dict) else (i if isinstance(i, str) else None)
        if _s(name):
            inst, basis = _s(name), "pi.institution.name"
    if inst is None and _s(rec.get("institute")):
        inst, basis = _s(rec.get("institute")), "institute"
    return ("on_campus" if inst == "UC Davis" else "off_campus"), inst, basis


def submission_label(summary: dict) -> str:
    """Folder-name identity: internal_id, or the hex id when a record has none."""
    return _s(summary.get("internal_id")) or _s(summary.get("id"))


def server_share_dir(summary: dict) -> str:
    """The share dir exactly as the summary recorded it (server-side). Never re-derived from
    a local root. A summary without one is refused rather than guessed."""
    d = _s(summary.get("share_dir"))
    if not d.startswith("/"):
        raise Stop(EXIT_DECIDE, {"error": f"summary has no server-side share_dir ({d!r})",
                                 "hint": "re-run `core_submission.py fetch`"})
    return d.rstrip("/")


def local_share_dir(summary: dict, warnings: list) -> str:
    """Where this machine writes the share. Checked against the summary's server path, so a
    summary for one submission can never steer files into another's share."""
    try:
        server = posixpath.join(server_project_dir(_s(summary.get("id")), _s(summary.get("submitted"))), "share")
        local = os.path.join(local_project_dir(_s(summary.get("id")), _s(summary.get("submitted"))), "share")
    except ValueError as e:
        raise Stop(EXIT_DECIDE, {"error": f"cannot derive the Bioshare share dir: {e}"})
    if server != server_share_dir(summary):
        raise Stop(EXIT_DECIDE, {"error": "summary share_dir does not match its id and submitted date",
                                 "summary_share_dir": summary.get("share_dir"), "derived": server})
    return local


# ------------------------------------------------------------------- CoreOmics API --
class _NoRedirectForWrites(urllib.request.HTTPRedirectHandler):
    """urllib silently turns a redirected POST into a GET, so `send --apply` once reported
    success with the share LISTING as its response. A write that redirects is an error."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if req.get_method() != "GET":
            raise urllib.error.HTTPError(req.full_url, code,
                                         f"{req.get_method()} was redirected to {newurl}; refusing to follow",
                                         headers, fp)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


_OPENER = urllib.request.build_opener(_NoRedirectForWrites)


def api_token() -> str:
    tok = _s(os.environ.get("COREOMICS_TOKEN"))
    if tok:
        return tok
    path = os.path.expanduser(TOKEN_FILE)
    try:
        with open(path) as fh:
            tok = fh.read().strip()
    except OSError as e:
        raise ApiError(None, f"no CoreOmics token: set COREOMICS_TOKEN or save your token in "
                             f"{TOKEN_FILE} ({e.strerror}). This runs on YOUR computer -- "
                             f"HIVE has no CoreOmics token.", None)
    if not tok:
        raise ApiError(None, f"{TOKEN_FILE} is empty", None)
    return tok


def api_url(path: str, params: dict | None = None) -> str:
    url = path if path.startswith(("http://", "https://")) else base_url() + "/" + path.lstrip("/")
    if params:
        url += ("&" if "?" in url else "?") + urllib.parse.urlencode(params)
    return url


def _error_detail(e: urllib.error.HTTPError) -> str:
    if 300 <= e.code < 400:
        return _s(e.msg)
    try:
        body = e.read().decode("utf-8", "replace")
    except Exception:
        body = ""
    try:
        j = json.loads(body)
        if isinstance(j, dict) and set(j) == {"detail"}:
            return _s(j["detail"])
        return json.dumps(j)
    except ValueError:
        return (body.strip()[:500] or _s(e.reason))


def api_call(method: str, path: str, params: dict | None = None, payload=None):
    url = api_url(path, params)
    headers = {"Authorization": f"Token {api_token()}", "Accept": "application/json"}
    data = None
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, method=method, headers=headers)
    try:
        with _OPENER.open(req, timeout=HTTP_TIMEOUT) as r:
            body = r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        raise ApiError(e.code, _error_detail(e), url) from None
    except urllib.error.URLError as e:
        raise ApiError(None, f"cannot reach {base_url()}: {e.reason}", url) from None
    except OSError as e:                                 # socket timeout, reset
        raise ApiError(None, f"cannot reach {base_url()}: {e}", url) from None
    if not body.strip():
        return {}
    try:
        return json.loads(body)
    except ValueError:
        raise ApiError(None, f"non-JSON response: {body[:200]}", url) from None


def api_get_all(path: str, params: dict | None = None, max_pages: int = 50) -> list:
    """A DRF endpoint may answer with a bare list or `{results, next}`. Handle both."""
    out, url, prm = [], path, params
    for _ in range(max_pages):
        data = api_call("GET", url, prm)
        if isinstance(data, list):
            return out + data
        if not isinstance(data, dict) or not isinstance(data.get("results"), list):
            raise ApiError(None, f"unexpected payload shape from {path}", api_url(path, params))
        out += data["results"]
        if not data.get("next"):
            return out
        url, prm = data["next"], None
    raise ApiError(None, f"more than {max_pages} pages from {path}", api_url(path, params))


def get_submission(kind: str, key: str) -> dict:
    if kind == "id":
        try:
            return api_call("GET", f"submissions/{key}/")
        except ApiError as e:
            if e.status == 404:
                raise Stop(EXIT_DECIDE, {"error": f"no CoreOmics submission with id {key}"})
            raise
    data = api_call("GET", "submissions/", {"lab": LAB, "internal_id": key})
    results = data.get("results") if isinstance(data, dict) else data
    # Filter again client-side: if the server ever ignored the filter we would otherwise
    # take the first submission of a whole page and call it PROT_0807.
    hits = [r for r in (results or []) if isinstance(r, dict) and _s(r.get("internal_id")) == key]
    if not hits:
        raise Stop(EXIT_DECIDE, {"error": f"no {LAB} submission with internal_id {key}",
                                 "hint": "check the number with the staff member"})
    if len(hits) > 1:
        raise Stop(EXIT_DECIDE, {"error": f"{len(hits)} submissions share internal_id {key}",
                                 "candidates": [{"id": h.get("id"), "submitted": h.get("submitted")}
                                                for h in hits],
                                 "hint": "re-run fetch with the 12-character id"})
    return hits[0]


def sample_rows(rec: dict) -> list:
    sd = rec.get("submission_data") or {}
    rows = []
    for smp in sd.get("samples") or []:
        if isinstance(smp, dict):
            rows.append({"unique_id": _s(smp.get("unique_id")),
                         "sample_name": _s(smp.get("sample_name")),
                         "condition_name": _s(smp.get("condition_name"))})
    return rows


def contact_rows(rec: dict) -> list:
    """The submission's extra contacts -- people `bioshare send` also grants access to."""
    raw = rec.get("contacts")
    if isinstance(raw, dict):
        raw = [raw]
    out = []
    for c in raw if isinstance(raw, list) else []:
        if isinstance(c, dict):
            name = _s(c.get("name")) or " ".join(x for x in (_s(c.get("first_name")), _s(c.get("last_name"))) if x)
            email = _s(c.get("email"))
        else:
            text = _s(c)
            name, email = ("", text) if "@" in text else (text, "")
        if name or email:
            out.append({"name": name or None, "email": email or None})
    return out


def submissions_between(lo: dt.date, hi: dt.date) -> tuple:
    """([record, ...] submitted in [lo, hi], truncated). The list is newest-first, so paging
    stops once a page reaches past the old edge. List results already carry the samples."""
    out, page = [], 1
    while True:
        data = api_call("GET", "submissions/", {"lab": LAB, "page": page, "page_size": 100,
                                                "ordering": "-submitted"})
        results = data.get("results") if isinstance(data, dict) else data
        if not isinstance(results, list):
            raise ApiError(None, "unexpected submissions list payload", api_url("submissions/"))
        oldest = None
        for r in results:
            d = parse_date(r.get("submitted")) if isinstance(r, dict) else None
            if d is None:
                continue
            oldest = d if oldest is None or d < oldest else oldest
            if lo <= d <= hi:
                out.append(r)
        if not results or (oldest is not None and oldest < lo):
            return out, False
        if not isinstance(data, dict) or not data.get("next"):
            return out, False
        page += 1
        if page > NEIGHBOR_MAX_PAGES:
            return out, True


def neighbor_row(r: dict) -> dict:
    """A submission as `locate` and `identify` weigh it: its ids, date and sample ids."""
    return {"internal_id": r.get("internal_id"), "id": r.get("id"), "submitted": r.get("submitted"),
            "unique_ids": [x["unique_id"] for x in sample_rows(r) if x["unique_id"]]}


def fetch_neighbors(ours: dict, ours_date: dt.date, window_days: int) -> tuple:
    """Submissions within +/-window of ours, with their unique_ids -- what `locate` needs to
    tell "our BN1" from another lab's BN1."""
    recs, truncated = submissions_between(ours_date - dt.timedelta(days=window_days),
                                          ours_date + dt.timedelta(days=window_days))
    return [neighbor_row(r) for r in recs if _s(r.get("id")) != _s(ours.get("id"))], truncated


def compact_share(s: dict) -> dict:
    keys = ("id", "submission", "bioshare_id", "name", "notes", "sub_folder", "link_to_path", "url")
    return {k: s.get(k) for k in keys}


def list_shares(sub_id: str) -> list:
    return [compact_share(s) for s in
            api_get_all(f"plugins/bioshare/submissions/{sub_id}/submission_shares/")
            if isinstance(s, dict)]


def build_summary(rec: dict, neighbors: list, shares, shares_error, window_days: int,
                  truncated: bool, warnings: list) -> dict:
    sd = rec.get("submission_data") or {}
    samples = sample_rows(rec)
    campus, institution, basis = classify_campus(rec)
    pi = rec.get("pi") if isinstance(rec.get("pi"), dict) else {}
    dept = pi.get("department")
    dept = dept.get("name") if isinstance(dept, dict) else dept
    pi_first, pi_last = _s(rec.get("pi_first_name")), _s(rec.get("pi_last_name"))
    sub_first, sub_last = _s(rec.get("first_name")), _s(rec.get("last_name"))
    types = sd.get("proteomics_type")
    types = [_s(t) for t in types if _s(t)] if isinstance(types, list) else ([_s(types)] if _s(types) else [])
    data_analysis = _s(sd.get("data_analysis"))
    groups = Counter(x["condition_name"] for x in samples if x["condition_name"])
    sub_date = parse_date(rec.get("submitted"))
    try:
        project = server_project_dir(_s(rec.get("id")), _s(rec.get("submitted")))
    except ValueError as e:
        project = None
        warnings.append(f"no canonical project dir: {e}")
    if basis == "missing":
        warnings.append("record has no PI institution and no institute -- campus defaulted to "
                        "off_campus; confirm with staff before staging")
    if not _s(rec.get("internal_id")):
        warnings.append("record has no internal_id; folders will use the hex id")
    organism = _s(sd.get("organism"))
    summary = {
        "schema": SCHEMA,
        "fetched_at": now_utc(),
        "internal_id": rec.get("internal_id"),
        "id": rec.get("id"),
        "label": _s(rec.get("internal_id")) or _s(rec.get("id")),
        "url": rec.get("url"),
        "status": rec.get("status"),
        "submitted": rec.get("submitted"),
        "submitted_date": sub_date.isoformat() if sub_date else None,
        "submitter": {"name": " ".join(x for x in (sub_first, sub_last) if x) or None,
                      "first": sub_first or None, "last": sub_last or None,
                      "email": _s(rec.get("email")) or None},
        "pi": {"name": " ".join(x for x in (pi_first, pi_last) if x) or None,
               "first": pi_first or None, "last": pi_last or None,
               "email": _s(rec.get("pi_email")) or None,
               "institution": institution, "department": _s(dept) or None},
        "contacts": contact_rows(rec),
        "campus": campus,
        "campus_basis": basis,
        # Golden rule #4: the organism is ASKED and CONFIRMED, never assumed. What the
        # submitter typed is a good first guess to put in front of the user -- nothing more.
        "organism_as_submitted": {
            "value": organism or None,
            "source": "submitter's free text on the CoreOmics form",
            "confirmed": False,
            "note": "a suggestion to CONFIRM with the user and resolve with "
                    "fetch_fasta.py resolve -- never a resolved organism"},
        "experiment_types": types,
        "instrument_wanted": _s(sd.get("mass_spec_wanted")) or None,
        "data_analysis": data_analysis or None,
        "raw_data_only": RAW_ONLY_MARKER in data_analysis.lower(),
        "description": _s(sd.get("description")) or None,
        "sample_prep": _s(sd.get("sample_prep")) or None,
        "n_samples": len(samples),
        "samples": samples,
        "conditions": {"present": bool(groups), "groups": dict(groups),
                       "missing": [x["unique_id"] for x in samples if not x["condition_name"]]},
        "server_root": server_flinders_root(),
        "canonical_project_dir": project,
        "share_dir": posixpath.join(project, "share") if project else None,
        "neighbor_window_days": window_days,
        "neighbors": neighbors,
        "neighbors_truncated": truncated,
        "existing_shares": shares,
    }
    if shares_error:
        summary["existing_shares_error"] = shares_error
    return summary


def write_tsv(path: str, header: list, rows: list) -> None:
    def clean(v):
        return re.sub(r"[\t\r\n]+", " ", _s(v))
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([clean(r.get(h)) for h in header])


def read_tsv(path: str) -> list:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_json_file(path: str, what: str) -> dict:
    try:
        with open(path) as fh:
            d = json.load(fh)
    except (OSError, ValueError) as e:
        raise Stop(EXIT_DECIDE, {"error": f"cannot read {what} {path}: {e}"})
    if not isinstance(d, dict):
        raise Stop(EXIT_DECIDE, {"error": f"{path} is not a {what}"})
    return d


def hive_summary(summary: dict) -> dict:
    """The summary as it may leave the staff member's computer: no email address, no contacts,
    and free text scrubbed of anything shaped like an email or phone number. locate, stage and
    deliver need none of those; `bioshare` and `email-draft` do, and they run locally on the
    full summary. (The record's allowlist is submission_report.py's; this is the summary's.)"""
    scrub = sibling("submission_report").scrub

    def walk(v):
        if isinstance(v, dict):
            return {k: walk(x) for k, x in v.items() if k not in ("email", "contacts")}
        if isinstance(v, list):
            return [walk(x) for x in v]
        return scrub(v) if isinstance(v, str) else v
    out = walk(summary)
    out["redacted_for_hive"] = True
    return out


def load_summary(path: str) -> dict:
    try:
        with open(path) as fh:
            s = json.load(fh)
    except (OSError, ValueError) as e:
        raise Stop(EXIT_DECIDE, {"error": f"cannot read summary {path}: {e}",
                                 "hint": "run `core_submission.py fetch` first (locally), then "
                                         "--put the summary to HIVE"})
    if not isinstance(s, dict) or not s.get("id"):
        raise Stop(EXIT_DECIDE, {"error": f"{path} is not a submission_summary.json"})
    return s


def cmd_fetch(a) -> int:
    try:
        kind, key = normalize_submission(a.submission)
    except ValueError as e:
        raise Stop(EXIT_DECIDE, {"error": str(e)})
    rec = get_submission(kind, key)
    warnings = []
    ours_date = parse_date(rec.get("submitted"))
    if ours_date is None:
        raise Stop(EXIT_DECIDE, {"error": f"submission {key} has no usable `submitted` date "
                                          f"({rec.get('submitted')!r})"})
    neighbors, truncated = fetch_neighbors(rec, ours_date, a.neighbor_days)
    if truncated:
        warnings.append(f"neighbour paging stopped at {NEIGHBOR_MAX_PAGES} pages; label-reuse "
                        f"checks in `locate` may be incomplete")
    shares, shares_error = [], None
    try:
        shares = list_shares(_s(rec.get("id")))
    except ApiError as e:
        # Reading shares is a courtesy for `bioshare status` later; it must not block the run.
        shares_error = f"HTTP {e.status}: {e.detail}" if e.status else e.detail
        warnings.append(f"could not list Bioshare shares: {shares_error}")

    summary = build_summary(rec, neighbors, shares, shares_error, a.neighbor_days,
                            truncated, warnings)
    os.makedirs(a.out, exist_ok=True)
    paths = {"submission_json": os.path.join(a.out, "submission.json"),
             "samples_tsv": os.path.join(a.out, "samples.tsv"),
             "summary": os.path.join(a.out, "submission_summary.json")}
    with open(paths["submission_json"], "w") as fh:
        json.dump(rec, fh, indent=2)
    write_tsv(paths["samples_tsv"], ["unique_id", "sample_name", "condition_name"], summary["samples"])
    with open(paths["summary"], "w") as fh:
        json.dump(summary, fh, indent=2)
    # What goes to HIVE: only these two. The raw record and the full summary hold emails and
    # billing (PPMS) fields and stay on this computer.
    hive = os.path.join(a.out, "hive")
    os.makedirs(hive, exist_ok=True)
    paths["hive_summary"] = os.path.join(hive, "submission_summary.json")
    paths["hive_record"] = os.path.join(hive, "submission.json")
    with open(paths["hive_summary"], "w") as fh:
        json.dump(hive_summary(summary), fh, indent=2)
    with open(paths["hive_record"], "w") as fh:
        json.dump(sibling("submission_report").sanitize(rec), fh, indent=2)
    for w in warnings:
        note(f"WARN {w}")
    emit({"internal_id": summary["internal_id"], "id": summary["id"],
          "submitted_date": summary["submitted_date"], "campus": summary["campus"],
          "pi": summary["pi"]["name"], "institution": summary["pi"]["institution"],
          "n_samples": summary["n_samples"], "conditions": summary["conditions"],
          "organism_as_submitted": summary["organism_as_submitted"],
          "instrument_wanted": summary["instrument_wanted"],
          "raw_data_only": summary["raw_data_only"], "share_dir": summary["share_dir"],
          "n_neighbors": len(neighbors), "n_existing_shares": len(shares),
          "warnings": warnings, "outputs": paths})
    return EXIT_OK if summary["share_dir"] else EXIT_DECIDE


# ------------------------------------------------------------------------- locate --
DATE_RULES = (
    # (pattern, digit layout). Order matters: the instrument prefixes are tried first so
    # `Ex08312026` is never read as a bare leading date.
    (re.compile(r"^ex(\d{8})(?!\d)", re.I), ("MMDDYYYY",)),
    (re.compile(r"^ex(\d{6})(?!\d)", re.I), ("DDMMYY",)),
    (re.compile(r"^fl(\d{6})(?!\d)", re.I), ("DDMMYY",)),
    (re.compile(r"^(\d{8})(?!\d)"), ("YYYYMMDD", "MMDDYYYY")),
)
EARLIEST_ACQUISITION = dt.date(2015, 1, 1)


def _layout_date(s: str, layout: str):
    try:
        if layout == "YYYYMMDD":
            return dt.date(int(s[:4]), int(s[4:6]), int(s[6:8]))
        if layout == "MMDDYYYY":
            return dt.date(int(s[4:8]), int(s[:2]), int(s[2:4]))
        if layout == "DDMMYY":
            return dt.date(2000 + int(s[4:6]), int(s[2:4]), int(s[:2]))
    except ValueError:
        return None
    return None


def date_from_name(name: str, today: dt.date | None = None):
    """Acquisition date encoded in a raw filename, or None.

    timsTOF: leading YYYYMMDD or MMDDYYYY (never both valid inside the accepted range).
    Exploris: `Ex` + 8 digits = MMDDYYYY, `Ex` + 6 = DDMMYY. Lumos: `FL` + 6 = DDMMYY.
    Only a real calendar date in [2015-01-01, today+1] counts -- the 9-digit typos seen
    on the share (`070622026__`) and undated Lumos names fall through to mtime."""
    latest = (today or dt.date.today()) + dt.timedelta(days=1)
    for rx, layouts in DATE_RULES:
        m = rx.match(name)
        if not m:
            continue
        for layout in layouts:
            d = _layout_date(m.group(1), layout)
            if d and EARLIEST_ACQUISITION <= d <= latest:
                return d
        return None
    return None


def acquisition_date(path: str, name: str | None = None) -> tuple:
    d = date_from_name(name or os.path.basename(path.rstrip("/")))
    if d:
        return d, "filename"
    try:
        return dt.date.fromtimestamp(os.stat(path).st_mtime), "mtime"
    except OSError:
        return None, "unknown"


def _kind(c: str) -> str:
    return "digit" if c.isdigit() else ("letter" if c.isalpha() else "other")


def token_pattern(uid: str):
    """A unique_id as a DELIMITED, case-insensitive token. `-`, `_` and space are the same
    separator on both sides, so EB_001 finds EB-001; letters/digits may not touch either
    end, so KG1 does not find KG13.

    Where letters meet digits a separator may also be added or dropped: PROT_0756's sheet
    says LRS96 and every one of its 30 runs says DIA-LRS-96 (measured on the share), so LRS96
    finds LRS-96 and EB_001 finds EB001. Between two runs of digits, or two of letters, the
    separator stays required -- SG-001-2 must not find SG0012."""
    parts = [p for p in re.split(r"[-_\s]+", _s(uid)) if p]
    if not parts:
        return None
    optional, required = r"[-_\s]?", r"[-_\s]+"

    def boundary(a: str, b: str) -> bool:
        return {_kind(a), _kind(b)} == {"letter", "digit"}

    def part(p: str) -> str:
        return "".join((optional if i and boundary(p[i - 1], c) else "") + re.escape(c)
                       for i, c in enumerate(p))

    body = part(parts[0])
    for prev, p in zip(parts, parts[1:]):
        body += (r"[-_\s]*" if boundary(prev[-1], p[0]) else required) + part(p)
    return re.compile(r"(?<![A-Za-z0-9])" + body + r"(?![A-Za-z0-9])", re.I)


# timsTOF runs carry the sample id in ONE fixed field:
#   MMDDYYYY__60SPD_DIA-<unique_id>_S3-<well>_1_<acq#>.d
# Matching anywhere else in such a name lets an id like A1-1 hit the `_S3-A1_1_` plate
# position that EVERY run has -- and succeed with another lab's files.
TIMS_NAME = re.compile(r"^\d{8,9}_+[^_]*_(?:DIA|DDA)-(?P<field>.+?)_S\d+-[A-Za-z]\d{1,2}_\d+_\d+(?:\.d)?$",
                       re.I)


def match_space(name: str) -> str:
    """The part of a raw filename a sample id may match: the DIA-<uid> field for names that
    follow the timsTOF convention, the whole name otherwise."""
    m = TIMS_NAME.match(name)
    return m.group("field") if m else name


WELL = re.compile(r"^[A-Ha-h](?:[1-9]|1[0-2])$")


def weak_reason(uid: str):
    """Why an id cannot be trusted to find its own files, or None. Measured: A3, H10 and 001
    each match hundreds to thousands of raw files. (Frequency across the whole listing is
    checked separately in match_samples.)"""
    u = _s(uid)
    if not u:
        return "blank unique_id"
    if WELL.match(u):
        return "looks like a plate well (A1-H12)"
    alnum = re.sub(r"[^A-Za-z0-9]", "", u)
    if len(alnum) < 3:
        return "fewer than 3 letters/digits"
    if alnum.isdigit():
        return "purely numeric"
    return None


def is_raw_name(name: str) -> bool:
    return name.lower().endswith(RAW_SUFFIXES)


def instrument_folder(wanted):
    return INSTRUMENT_FOLDERS.get(_s(wanted).lower())


def list_raw(root: str, first=None) -> tuple:
    """Every raw entry under <root>/<instrument>/<month folder>/. Month folder names are
    free-form (sep26, JUL26, jan25AndPM, Std_He_...), so ALL of them are scanned."""
    entries, unreadable = [], []
    try:
        insts = sorted(os.listdir(root))
    except OSError as e:
        raise Stop(EXIT_UNREACHABLE, {"error": f"cannot list {root}: {e.strerror}"})
    if first in insts:
        insts.remove(first)
        insts.insert(0, first)
    for inst in insts:
        ip = os.path.join(root, inst)
        if inst.startswith(".") or not os.path.isdir(ip):
            continue
        try:
            children = sorted(os.listdir(ip))
        except OSError as e:
            unreadable.append(f"{ip}: {e.strerror}")
            continue
        for c in children:
            cp = os.path.join(ip, c)
            if is_raw_name(c):
                entries.append({"path": cp, "name": c, "instrument_folder": inst})
                continue
            if c.startswith(".") or not os.path.isdir(cp):
                continue
            try:
                grand = sorted(os.listdir(cp))
            except OSError as e:
                unreadable.append(f"{cp}: {e.strerror}")
                continue
            for g in grand:
                if is_raw_name(g):
                    entries.append({"path": os.path.join(cp, g), "name": g, "instrument_folder": inst})
    return entries, unreadable


def prepare_entry(e: dict) -> dict:
    e["space"] = match_space(e["name"])
    e["norm"] = normkey(e["space"])
    return e


def ht_pattern(internal_id):
    """HT plate filenames put the submission number straight after the leading date:
    `20260827_793_100spd_...`, `20260901_0793_rerun_...`. Anchoring there (rather than
    anywhere in the name) keeps run counters such as Exploris `Ex08312026_380_JE21` from
    impersonating submission 0380."""
    m = re.search(r"(\d+)$", _s(internal_id))
    if not m or int(m.group(1)) == 0:
        return None
    n = int(m.group(1))
    nums = sorted({str(n), "%04d" % n}, key=len, reverse=True)
    return re.compile(r"^\d{8}_(?:" + "|".join(nums) + r")(?![A-Za-z0-9])", re.I)


def _mtime(path: str) -> float:
    try:
        return os.stat(path).st_mtime
    except OSError:
        return 0.0


def _row(smp: dict, **kw) -> dict:
    r = {"unique_id": smp.get("unique_id", ""), "sample_name": smp.get("sample_name", ""),
         "condition_name": smp.get("condition_name", ""), "file": "", "acquired": "",
         "date_source": "", "status": "", "alternates": [], "note": "",
         "instrument_folder": None, "ambiguous": [], "accepted_ambiguous": False,
         "out_of_window": 0, "shadowed": []}
    r.update(kw)
    return r


def id_universe(samples: list, neighbors: list, ours_label: str, ours_date) -> list:
    """Every sample id in play -- ours and every neighbour's -- as a compiled pattern."""
    universe = []

    def add(uid, sub, submitted, ours, smp=None):
        pat, key = token_pattern(uid), normkey(uid)
        if pat and key:
            universe.append({"uid": _s(uid), "norm": key, "pattern": pat, "sub": sub,
                             "submitted": submitted, "ours": ours, "sample": smp})

    for smp in samples:
        add(smp.get("unique_id"), ours_label, ours_date, True, smp)
    for n in neighbors or []:
        nid, nd = _s(n.get("internal_id")) or _s(n.get("id")), parse_date(n.get("submitted"))
        for u in n.get("unique_ids") or []:
            add(u, nid, nd, False)
    return universe


def file_owners(e: dict, universe: list, cache: dict, keep=None) -> list:
    """The ids that own a file: of every id whose pattern matches it (and that `keep` accepts),
    the LONGEST. `DH1` matches `DH1-1`'s run too, but the run belongs to DH1-1."""
    if keep is not None or e["path"] not in cache:
        cands = [u for u in universe if u["norm"] in e["norm"] and u["pattern"].search(e["space"])
                 and (keep is None or keep(u))]
        top = max((len(u["norm"]) for u in cands), default=0)
        owners = [u for u in cands if len(u["norm"]) == top]
        if keep is not None:
            return owners
        cache[e["path"]] = owners
    return cache[e["path"]]


def match_samples(samples: list, entries: list, neighbors: list, lo: dt.date, hi: dt.date,
                  ours_label: str = "", accept_ambiguous: bool = False) -> list:
    """The deterministic core of `locate`: one row per sample.

    A file in the window is AMBIGUOUS when another submission -- older, same day or newer --
    uses the same label and was submitted on or before the file was acquired: either lab
    could have sent that sample, and the label alone cannot say which."""
    for e in entries:
        if "space" not in e:
            prepare_entry(e)
    universe = id_universe(samples, neighbors, ours_label, lo)
    by_norm = {}
    for u in universe:
        by_norm.setdefault(u["norm"], set()).add((u["submitted"], u["sub"]))
    owner_cache, date_cache = {}, {}

    def dated(e):
        if e["path"] not in date_cache:
            date_cache[e["path"]] = acquisition_date(e["path"], e["name"])
        return date_cache[e["path"]]

    def earlier_unambiguous(label_norm, sub, before, owned):
        """Files with this label that could only be `sub`'s: acquired after it was submitted
        and before any other user of the label was."""
        users = sorted((d, s) for d, s in by_norm.get(label_norm, ()) if d is not None)
        mine = [d for d, s in users if s == sub]
        if not mine:
            return []
        start = min(mine)
        others = [d for d, s in users if s != sub]
        if any(d <= start for d in others):
            return []                                 # someone else used it first or same day
        end = min([d for d in others if d > start] + [before])
        return [e["path"] for e in owned if dated(e)[0] is not None and start <= dated(e)[0] < end]

    rows = []
    for smp in samples:
        uid = smp.get("unique_id", "")
        pat, key = token_pattern(uid), normkey(uid)
        reason = weak_reason(uid)
        if reason or not pat:
            n = sum(1 for e in entries if key and key in e["norm"] and pat.search(e["space"])) if pat else 0
            rows.append(_row(smp, status="weak",
                             note=f"{reason}; {n} raw file(s) carry this token -- not auto-assigned"))
            continue
        hits = [e for e in entries if key in e["norm"] and pat.search(e["space"])]
        owned, shadowed = [], []
        for e in hits:
            owners = file_owners(e, universe, owner_cache)
            if any(o["norm"] == key for o in owners):
                owned.append(e)
            else:
                shadowed.append({"file": e["path"], "owned_by": sorted({f"{o['uid']} ({o['sub']})" for o in owners})})
        if len(owned) > WEAK_MAX_FILES:
            rows.append(_row(smp, status="weak", shadowed=shadowed[:5],
                             note=f"this id names {len(owned)} raw files across all instruments and "
                                  f"dates (more than {WEAK_MAX_FILES}) -- too common to pick out this "
                                  f"submission's runs; not auto-assigned"))
            continue
        in_win, n_out = [], 0
        for e in owned:
            d, src = dated(e)
            if d is not None and lo <= d <= hi:
                in_win.append((e, d, src))
            else:
                n_out += 1
        free, amb = [], []
        for e, d, src in in_win:
            contenders = {}
            for o in file_owners(e, universe, owner_cache):
                if not o["ours"] and o["submitted"] is not None and o["submitted"] <= d:
                    contenders.setdefault((o["sub"], o["norm"]), o)
            (amb if contenders else free).append((e, d, src, list(contenders.values())))
        amb_info = []
        for e, d, src, cons in amb:
            detail = []
            for o in cons:
                earlier = earlier_unambiguous(o["norm"], o["sub"], d, owned)
                detail.append({"submission": o["sub"], "submitted": o["submitted"].isoformat(),
                               "label": o["uid"], "has_earlier_unambiguous_files": bool(earlier),
                               "earlier_files": len(earlier), "earlier_examples": earlier[:3]})
            amb_info.append({"file": e["path"], "acquired": d.isoformat(), "contenders": detail})
        pool = free or (amb if accept_ambiguous else [])
        if pool:
            pool.sort(key=lambda t: (t[1], _mtime(t[0]["path"]), t[0]["name"]), reverse=True)
            e, d, src, _ = pool[0]
            accepted = not free
            rows.append(_row(smp, status="matched", file=e["path"], acquired=d.isoformat(),
                             date_source=src, instrument_folder=e["instrument_folder"],
                             alternates=[x[0]["path"] for x in pool[1:]], ambiguous=amb_info,
                             accepted_ambiguous=accepted, out_of_window=n_out, shadowed=shadowed[:5],
                             note="ambiguous label accepted with --accept-ambiguous" if accepted else ""))
        elif amb:
            amb.sort(key=lambda t: (t[1], t[0]["name"]), reverse=True)
            e, d, src, _ = amb[0]
            subs = sorted({c["submission"] for i in amb_info for c in i["contenders"]})
            rows.append(_row(smp, status="ambiguous_label", file=e["path"], acquired=d.isoformat(),
                             date_source=src, instrument_folder=e["instrument_folder"],
                             ambiguous=amb_info, out_of_window=n_out, shadowed=shadowed[:5],
                             note=f"{', '.join(subs)} also use this label and were submitted before "
                                  f"these runs -- the label alone cannot say whose runs they are"))
        else:
            bits = []
            if n_out:
                bits.append(f"{n_out} file(s) carry this id but were acquired outside the window")
            if shadowed:
                bits.append(f"{len(shadowed)} file(s) carry it inside a longer id")
            rows.append(_row(smp, status="unmatched", out_of_window=n_out, shadowed=shadowed[:5],
                             note="; ".join(bits) or "no raw file carries this id"))
    return rows


def curated_rows(samples: list, paths: list, raw: str) -> list:
    """--files-from: staff chose the files. Map each back to the sample whose id is the
    longest one it carries; a curated list makes even a weak id usable."""
    universe = id_universe(samples, [], "", None)
    by_sample, loose = {}, []
    for p in paths:
        e = prepare_entry({"path": p, "name": os.path.basename(p.rstrip("/"))})
        owners = file_owners(e, universe, {})
        who = {id(o["sample"]): o["sample"] for o in owners}
        if len(who) == 1:
            by_sample.setdefault(next(iter(who)), []).append(p)
        else:
            loose.append((p, sorted({o["uid"] for o in owners})))
    rows = []
    for smp in samples:
        files = by_sample.get(id(smp), [])
        if not files:
            rows.append(_row(smp, status="unmatched", note="no file in the curated list carries this id"))
        for p in files:
            d, src = acquisition_date(p)
            folder = os.path.relpath(p, raw).split(os.sep)[0] if is_under(p, raw) else None
            rows.append(_row(smp, status="matched", file=p, acquired=d.isoformat() if d else "",
                             date_source=src, instrument_folder=folder, note="curated"))
    for p, cands in loose:
        d, src = acquisition_date(p)
        rows.append(_row({}, status="unassigned", file=p, acquired=d.isoformat() if d else "",
                         date_source=src,
                         note=("matches several samples: " + ", ".join(cands)) if cands
                         else "matches no sample id"))
    return rows


def _gate(name, status, detail="", **kw):
    g = {"gate": name, "status": status}
    if detail:
        g["detail"] = detail
    g.update(kw)
    return g


def locate_gates(rows: list, files: list, allow_partial: bool, wanted_folder, unreadable: list,
                 curated: bool, missing_paths: list, accept_ambiguous: bool = False) -> list:
    """PASS/WARN/INFO/FAIL per check. Any FAIL is a hard gate: exit 2, nothing staged."""
    gates = []
    by = {}
    for r in rows:
        by.setdefault(r["status"], []).append(r)
    if unreadable:
        gates.append(_gate("listing", "WARN", f"{len(unreadable)} folder(s) could not be read",
                           examples=unreadable[:5]))
    if curated:
        gates.append(_gate("paths_exist", "FAIL" if missing_paths else "PASS",
                           f"{len(missing_paths)} curated path(s) do not exist" if missing_paths else "",
                           examples=missing_paths[:10]))
    gates.append(_gate("no_files", "FAIL" if not files else "PASS",
                       "no raw files chosen" if not files else "", n=len(files)))
    unmatched = by.get("unmatched", [])
    gates.append(_gate("unmatched_samples",
                       ("WARN" if allow_partial else "FAIL") if unmatched else "PASS",
                       (f"{len(unmatched)} sample(s) have no raw file"
                        + (" -- accepted with --allow-partial" if allow_partial else
                           ". If staff confirm they were not run, re-run with --allow-partial"))
                       if unmatched else "",
                       samples=[{"unique_id": r["unique_id"], "note": r["note"]} for r in unmatched]))
    if not curated:
        weak = by.get("weak", [])
        gates.append(_gate("weak_ids", "FAIL" if weak else "PASS",
                           (f"{len(weak)} id(s) too short, numeric, well-like or too common to find "
                            f"their own files. Staff pick the files, then re-run with --files-from")
                           if weak else "",
                           samples=[{"unique_id": r["unique_id"], "note": r["note"]} for r in weak]))
        amb = by.get("ambiguous_label", [])
        accepted = [r for r in by.get("matched", []) if r["accepted_ambiguous"]]
        if amb:
            gates.append(_gate("ambiguous_label", "FAIL",
                               f"{len(amb)} sample(s) whose only candidate runs carry a label another "
                               f"submission (submitted on or before the run) also uses. Staff decide: "
                               f"--files-from with the right files, or --accept-ambiguous if these runs "
                               f"are this submission's",
                               samples=[{"unique_id": r["unique_id"], "files": r["ambiguous"]} for r in amb]))
        elif accepted:
            gates.append(_gate("ambiguous_label", "WARN",
                               f"{len(accepted)} sample(s) use runs with a shared label -- accepted with "
                               f"--accept-ambiguous",
                               samples=[{"unique_id": r["unique_id"], "chosen": r["file"],
                                         "files": r["ambiguous"]} for r in accepted]))
        else:
            gates.append(_gate("ambiguous_label", "PASS"))
        partial = [r for r in by.get("matched", []) if r["ambiguous"] and not r["accepted_ambiguous"]]
        if partial:
            gates.append(_gate("ambiguous_files_excluded", "WARN",
                               f"{len(partial)} sample(s) also have runs with a shared label; those were "
                               f"set aside and an unambiguous run was chosen",
                               samples=[{"unique_id": r["unique_id"], "chosen": r["file"],
                                         "set_aside": r["ambiguous"]} for r in partial]))
        shadowed = [r for r in rows if r["shadowed"]]
        if shadowed:
            gates.append(_gate("shadowed_by_longer_id", "INFO",
                               f"{len(shadowed)} id(s) also appear inside a longer id's run names "
                               f"(DH1 inside DH1-1); those runs were left to the longer id",
                               samples=[{"unique_id": r["unique_id"], "files": r["shadowed"]} for r in shadowed]))
        alts = [r for r in by.get("matched", []) if r["alternates"]]
        if alts:
            gates.append(_gate("alternates", "WARN",
                               f"{len(alts)} sample(s) have several files (re-injections?); the most "
                               f"recent was chosen", samples=[{"unique_id": r["unique_id"],
                                                               "chosen": r["file"],
                                                               "alternates": r["alternates"]} for r in alts]))
        n_out = sum(r["out_of_window"] for r in rows)
        if n_out:
            gates.append(_gate("out_of_window", "INFO",
                               f"{n_out} file(s) carry a sample id but were acquired before the "
                               f"submission or after the window (older runs of the same label)"))
    else:
        loose = by.get("unassigned", [])
        if loose:
            gates.append(_gate("unassigned_files", "WARN",
                               f"{len(loose)} curated file(s) map to no single sample; they are "
                               f"searched, but `conditions` will ask which group they belong to",
                               files=[{"file": r["file"], "note": r["note"]} for r in loose]))
    counts = Counter(r["file"] for r in rows if r["status"] == "matched")
    dupes = sorted(f for f, c in counts.items() if c > 1)
    gates.append(_gate("duplicate_assignment", "FAIL" if dupes else "PASS",
                       f"{len(dupes)} file(s) were assigned to more than one sample" if dupes else "",
                       files=dupes[:10]))
    mt = [r for r in rows if r["file"] and r["status"] == "matched" and r["date_source"] != "filename"]
    if mt:
        gates.append(_gate("mtime_dates", "WARN",
                           f"{len(mt)} chosen file(s) have no parseable date in the name; dated by "
                           f"modification time", files=[r["file"] for r in mt][:10]))
    if wanted_folder:
        elsewhere = sorted({r["instrument_folder"] for r in rows if r["status"] == "matched"
                            and r["instrument_folder"] and r["instrument_folder"] != wanted_folder})
        if elsewhere:
            gates.append(_gate("instrument", "WARN",
                               f"submitter asked for {wanted_folder}; some files were acquired on "
                               f"{', '.join(elsewhere)}"))
    return gates


def cmd_locate(a) -> int:
    s = load_summary(a.summary)
    ours = parse_date(s.get("submitted_date") or s.get("submitted"))
    if ours is None:
        raise Stop(EXIT_DECIDE, {"error": "summary has no usable submitted date"})
    window = s.get("neighbor_window_days")
    if not isinstance(window, int):
        raise Stop(EXIT_DECIDE, {"error": "summary does not record its neighbour window",
                                 "hint": "re-run `core_submission.py fetch`"})
    if a.max_days > window:
        raise Stop(EXIT_DECIDE, {"error": f"--max-days {a.max_days} is wider than the {window}-day "
                                          f"neighbour window fetch recorded; label reuse could not be "
                                          f"checked for the extra days",
                                 "hint": f"re-run fetch with --neighbor-days {a.max_days}, or keep "
                                         f"--max-days <= {window}"})
    lo, hi = ours, ours + dt.timedelta(days=a.max_days)
    samples = s.get("samples") or []
    raw = os.path.abspath(a.raw_root or raw_root())
    wanted = instrument_folder(s.get("instrument_wanted"))
    os.makedirs(a.out, exist_ok=True)
    outputs = {"files_txt": os.path.join(a.out, "files.txt"),
               "sample_files_tsv": os.path.join(a.out, "sample_files.tsv"),
               "locate_json": os.path.join(a.out, "locate.json")}
    result = {"internal_id": s.get("internal_id"), "id": s.get("id"),
              "window": {"from": lo.isoformat(), "to": hi.isoformat(), "max_days": a.max_days,
                         "neighbor_window_days": window},
              "mode": "files_from" if a.files_from else "auto", "raw_root": raw,
              "accepted": {"allow_partial": bool(a.allow_partial),
                           "accept_ambiguous": bool(a.accept_ambiguous)}}

    if a.files_from:
        try:
            with open(a.files_from) as fh:
                paths = [ln.strip() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
        except OSError as e:
            raise Stop(EXIT_DECIDE, {"error": f"cannot read --files-from {a.files_from}: {e.strerror}"})
        paths = [os.path.abspath(os.path.expanduser(p.rstrip("/"))) for p in paths]
        missing = [p for p in paths if not os.path.exists(p)]
        existing = list(dict.fromkeys(p for p in paths if os.path.exists(p)))
        rows = curated_rows(samples, existing, raw)
        files = existing
        gates = locate_gates(rows, files, a.allow_partial, wanted, [], True, missing)
    else:
        if not os.path.isdir(raw):
            raise Stop(EXIT_UNREACHABLE, {
                "error": f"raw data root not reachable: {raw}",
                "hint": "locate runs ON HIVE (hive_exec.sh); elsewhere set CORE_FLINDERS_ROOT or --raw-root"})
        entries, unreadable = list_raw(raw, wanted)
        result["n_raw_entries"] = len(entries)
        pat = ht_pattern(s.get("internal_id"))
        ht = []
        if pat:
            for e in entries:
                if pat.match(e["name"]):
                    d, _src = acquisition_date(e["path"], e["name"])
                    if d is not None and lo <= d <= hi:
                        ht.append(e["path"])
        if ht:
            result.update(ht_plate=True, n_ht_files=len(ht), ht_examples=ht[:10],
                          next_step=("this is a high-throughput plate: get the file list from STAN -- "
                                     "python3 ~/proteomics-pipeline/scripts/ht_manifest.py fetch "
                                     + re.sub(r"\D", "", _s(s.get("internal_id"))) + " --out <dir>"),
                          outputs={"locate_json": outputs["locate_json"]})
            with open(outputs["locate_json"], "w") as fh:
                json.dump(result, fh, indent=2)
            note(f"{len(ht)} file(s) carry the plate token -- use ht_manifest.py (references/ht-submissions.md)")
            emit(result)
            return EXIT_HT
        rows = match_samples(samples, entries, s.get("neighbors") or [], lo, hi,
                             submission_label(s), a.accept_ambiguous)
        files = list(dict.fromkeys(r["file"] for r in rows if r["status"] == "matched"))
        gates = locate_gates(rows, files, a.allow_partial, wanted, unreadable, False, [],
                             a.accept_ambiguous)
        if s.get("neighbors_truncated"):
            gates.append(_gate("neighbors", "WARN", "fetch stopped paging neighbours early; label "
                                                    "reuse may be under-reported"))

    hard_fail = any(g["status"] == "FAIL" for g in gates)
    with open(outputs["files_txt"], "w") as fh:
        fh.write("".join(f + "\n" for f in files))
    tsv_rows = [dict(r, alternates=";".join(r["alternates"])) for r in rows]
    write_tsv(outputs["sample_files_tsv"],
              ["unique_id", "sample_name", "condition_name", "file", "acquired", "date_source",
               "status", "alternates", "note"], tsv_rows)
    result.update(counts=dict(Counter(r["status"] for r in rows), samples=len(samples), files=len(files)),
                  gates=gates, hard_fail=hard_fail, files=files, samples=rows, outputs=outputs)
    with open(outputs["locate_json"], "w") as fh:
        json.dump(result, fh, indent=2)
    for g in gates:
        if g["status"] not in ("PASS",):
            note(f"[{g['status']}] {g['gate']}: {g.get('detail', '')}")
    emit({k: v for k, v in result.items() if k not in ("samples", "files")})
    if hard_fail:
        note("HARD GATE FAILED -- show the staff member the proposal; nothing was staged")
        return EXIT_DECIDE
    return EXIT_OK


# ----------------------------------------------------------------------- identify --
# Which submission is this data? Asked at the START of a run (SKILL.md step 1), so every report
# can carry the submission. Named ids are read only where they are written on purpose: a
# PROT_#### token in a path or message, a 12-hex CoreOmics id, or "submission 807" typed by the
# user. A bare number in a FILE name is never read as one -- Exploris runs carry counters
# (Ex08312026_380_JE21). Everything else goes through the same sample-id matching as `locate`.
# 3-4 digits after PROT: "Total_prot_10ug" and "WT_prot1" are protein amounts and names, not
# submissions. A hex id preceded by "xxxx-" is the tail of a UUID/GUID, not a CoreOmics id.
PROT_TOKEN = re.compile(r"(?<![A-Za-z0-9])prot[_\-# ]?(\d{3,4})(?![A-Za-z0-9])", re.I)
SUBMISSION_WORD = re.compile(r"\bsubmission\s*(?:number|no\.?)?\s*#?\s*(\d{3,4})(?![0-9])", re.I)
HEX_TOKEN = re.compile(r"(?<![0-9a-z])(?<![0-9a-f]{4}-)([0-9a-f]{12})(?![0-9a-z])")
# identify's exit 0 needs more than one lucky id: a generic id ("HeLa", "Blank") can name one
# file of thirty. The candidate must own at least half the files and two distinct sheet ids.
MIN_FILE_SHARE, MIN_SHEET_IDS = 0.5, 2
STALE_SNAPSHOT = "/quobyte/proteomics-grp/coreomics/.submissions_db"
TOKEN_HELP = (
    "Looking a submission up needs a CoreOmics API token with staff access, saved on THIS "
    f"computer as {TOKEN_FILE} (chmod 600) or exported as COREOMICS_TOKEN; HIVE has none. How "
    "to obtain one is not documented in this skill -- ask the Core's CoreOmics administrator.")


def ids_in(text: str, typed: bool = False) -> list:
    """[("internal_id", "PROT_0756") | ("id", "0022066cd85f"), ...] written in a path or name,
    or -- typed=True -- in what the user wrote. A 12-hex token counts only with a letter AND a
    digit in it: a 12-digit timestamp is not a CoreOmics id."""
    s, hits = _s(text), []
    for m in PROT_TOKEN.finditer(s):
        if int(m.group(1)):
            hits.append(("internal_id", "PROT_%04d" % int(m.group(1))))
    if typed:
        for m in SUBMISSION_WORD.finditer(s):
            if int(m.group(1)):
                hits.append(("internal_id", "PROT_%04d" % int(m.group(1))))
    for m in HEX_TOKEN.finditer(s):
        h = m.group(1)
        if re.search(r"[a-f]", h) and re.search(r"\d", h):
            hits.append(("id", h))
    return list(dict.fromkeys(hits))


def attribute_files(names: list, subs: list, max_days: int) -> dict:
    """Which submission's sample ids each raw file name carries -- `locate`'s rules run in
    reverse, on its own pieces (id_universe, file_owners). Weak ids (A3, 001) are no evidence.
    A file dated in its NAME counts only for a submission made on or before that date and
    within `max_days` of it; an undated name (or a modification time, which a copy resets)
    counts for every submission using the id. A file two submissions could own is ambiguous,
    exactly like `locate`'s ambiguous_label. `subs` are neighbor_row()s."""
    universe = [u for u in id_universe([], subs, "", None) if not weak_reason(u["uid"])]
    by_sub, ambiguous, loose = {}, [], []
    for n in names:
        e = prepare_entry({"path": n, "name": os.path.basename(_s(n).rstrip("/\\"))})
        d = date_from_name(e["name"])
        owners = file_owners(e, universe, {}, keep=lambda u: d is None or (
            u["submitted"] is not None and u["submitted"] <= d <= u["submitted"] + dt.timedelta(days=max_days)))
        who = sorted({u["sub"] for u in owners})
        if len(who) == 1:
            b = by_sub.setdefault(who[0], {"files": [], "sample_ids": set()})
            b["files"].append(e["name"])
            b["sample_ids"].update(u["uid"] for u in owners)
        elif who:
            ambiguous.append({"file": e["name"], "submissions": who, "dated": bool(d)})
        else:
            loose.append(e["name"])
    return {"by_submission": by_sub, "ambiguous": ambiguous, "unattributed": loose}


def _names_from(paths: list, files_from) -> list:
    """The paths as given, plus the raw entries one level inside any folder given (a folder
    alone has no file names to match)."""
    names = list(paths or [])
    if files_from:
        try:
            with open(files_from) as fh:
                names += [ln.strip() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
        except OSError as e:
            raise Stop(EXIT_DECIDE, {"error": f"cannot read --files-from {files_from}: {e.strerror}"})
    for p in list(names):
        if os.path.isdir(p) and not is_raw_name(os.path.basename(p.rstrip("/\\"))):
            try:
                names += [os.path.join(p, c) for c in sorted(os.listdir(p)) if is_raw_name(c)]
            except OSError:
                pass                        # an unreadable folder adds no names; the path still counts
    return list(dict.fromkeys(names))


def cmd_identify(a) -> int:
    names = _names_from(a.paths, a.files_from)
    raws = [n for n in names if is_raw_name(os.path.basename(_s(n).rstrip("/\\")))]
    named = {}
    for n in names:
        for hit in ids_in(n):
            named.setdefault(hit, []).append(f"in the path {n}")
    for hit in ids_in(a.text or "", typed=True):
        named.setdefault(hit, []).append("in the user's message")
    base = {"never_search": f"{STALE_SNAPSHOT} is a stale snapshot -- never search it for a "
                            f"submission (a search for 0756 there matched an unrelated 2019 record)"}
    ask_none = ("Which CoreOmics submission is this data from? Give the PROT number (e.g. "
                "PROT_0756) or the 12-character CoreOmics id.")

    def label_of(kind, key):
        return key if kind == "internal_id" else f"CoreOmics id {key}"

    def lookup_failed(e):
        emit(dict(base, status="lookup_failed", detail=e.detail, token_help=TOKEN_HELP,
                  ask=ask_none + " (CoreOmics could not be asked: " + e.detail + ")"))
        return EXIT_UNREACHABLE

    have_token = False
    if not a.no_lookup:
        try:
            api_token()
            have_token = True
        except ApiError:
            pass

    if named and have_token:
        # Look every named id up: a PROT number and a hex id are often one submission, and a
        # named submission whose sample ids appear in NONE of the file names is asked about.
        merged, found = {}, {}
        try:
            for (kind, key), where in named.items():
                try:
                    rec = get_submission(kind, key)
                except Stop as e:
                    emit(dict(base, status="not_found", submission=key, detail=e.payload.get("error"),
                              ask=f"CoreOmics has no submission {label_of(kind, key)} ({where[0]}). "
                                  f"Which submission is this data from?"))
                    return EXIT_DECIDE
                ident = _s(rec.get("internal_id")) or _s(rec.get("id"))
                k = ("internal_id", ident) if ident.upper().startswith("PROT_") else ("id", ident)
                merged.setdefault(k, []).extend(where)
                found[k] = rec
        except ApiError as e:
            return lookup_failed(e)
        named = merged
        if len(named) == 1:
            (kind, key), where = next(iter(named.items()))
            rec = found[(kind, key)]
            res = attribute_files(raws, [neighbor_row(rec)], a.max_days) if raws else None
            hit = (res or {}).get("by_submission", {}).get(key)
            n_sheet = len(sample_rows(rec))
            if raws and not hit:
                emit(dict(base, status="named_unconfirmed", submission=key, kind=kind,
                          evidence=where[:5],
                          ask=f"{label_of(kind, key)} is named {where[0]}, but none of its "
                              f"{n_sheet} sample IDs appear in these {len(raws)} file names. Is "
                              f"this data from {label_of(kind, key)}?"))
                return EXIT_DECIDE
            cover = (f": {len(hit['files'])} of {len(raws)} file names carry its sample IDs"
                     if hit else "")
            emit(dict(base, status="named", submission=key, kind=kind, evidence=where[:5],
                      verified=bool(hit), files_matched=len(hit["files"]) if hit else None,
                      ask=f"This looks like CoreOmics submission {label_of(kind, key)} "
                          f"({where[0]}{cover}). Is that right?"))
            return EXIT_OK
    if len(named) == 1:
        (kind, key), where = next(iter(named.items()))
        emit(dict(base, status="named", submission=key, kind=kind, evidence=where[:5], verified=False,
                  ask=f"This looks like CoreOmics submission {label_of(kind, key)} "
                      f"({where[0]}). Is that right?"))
        return EXIT_OK
    if len(named) > 1:
        cands = [label_of(k, v) for k, v in named]
        emit(dict(base, status="ambiguous", candidates=[{"kind": k, "submission": v, "evidence": w[:5]}
                                                        for (k, v), w in named.items()],
                  ask=f"The names and message mention {', '.join(cands)} -- which submission is "
                      f"this data from?"))
        return EXIT_DECIDE

    if a.no_lookup or not raws:
        emit(dict(base, status="none", looked_up=False, ask=ask_none))
        return EXIT_DECIDE
    if not have_token:
        emit(dict(base, status="needs_token", looked_up=False,
                  ask=ask_none + " (No PROT number is in the names, and without a CoreOmics token "
                                 "I cannot look it up by sample ID.)", token_help=TOKEN_HELP))
        return EXIT_UNREACHABLE
    dates = [d for d in (date_from_name(os.path.basename(_s(n).rstrip("/\\"))) for n in raws) if d]
    today = dt.date.today()
    lo = (min(dates) if dates else today - dt.timedelta(days=a.max_days)) - dt.timedelta(days=a.max_days)
    hi = max(dates) if dates else today
    try:
        recs, truncated = submissions_between(lo, hi)
    except ApiError as e:
        return lookup_failed(e)
    subs = [neighbor_row(r) for r in recs]
    res = attribute_files(raws, subs, a.max_days)
    window = {"from": lo.isoformat(), "to": hi.isoformat(), "submissions_checked": len(subs),
              "truncated": truncated, "files": len(raws), "files_dated_by_name": len(dates)}
    by = res["by_submission"]
    sizes = {(_s(s.get("internal_id")) or _s(s.get("id"))): len(s["unique_ids"]) for s in subs}
    cands = [{"submission": k, "files_matched": len(v["files"]), "files": len(raws),
              "sample_ids_matched": len(v["sample_ids"]), "samples_on_sheet": sizes.get(k),
              "examples": v["files"][:3]} for k, v in sorted(by.items())]
    out = dict(base, window=window, candidates=cands, ambiguous_files=res["ambiguous"][:20],
               n_unattributed=len(res["unattributed"]), unattributed_examples=res["unattributed"][:5])
    one = len(by) == 1 and not res["ambiguous"]
    if one and not truncated and cands[0]["files_matched"] >= MIN_FILE_SHARE * len(raws) \
            and cands[0]["sample_ids_matched"] >= min(MIN_SHEET_IDS, len(raws)):
        c = cands[0]
        emit(dict(out, status="matched", submission=c["submission"],
                  ask=f"These files look like CoreOmics submission {c['submission']}: "
                      f"{c['files_matched']} of {c['files']} file names carry its sample IDs "
                      f"({c['sample_ids_matched']} of {c['samples_on_sheet']} on its sheet). Is that right?"))
        return EXIT_OK
    if one:
        c = cands[0]
        why = ("the submission list was cut off before every submission in the window was checked"
               if truncated else f"only {c['files_matched']} of {c['files']} file names carry its "
                                 f"sample IDs")
        emit(dict(out, status="weak", submission=c["submission"],
                  ask=f"These files may be CoreOmics submission {c['submission']}, but {why}. "
                      f"Which submission is this data from?"))
        return EXIT_DECIDE
    if by or res["ambiguous"]:
        who = sorted(set(by) | {s for f in res["ambiguous"] for s in f["submissions"]})
        emit(dict(out, status="ambiguous",
                  ask=f"These file names carry sample IDs from {', '.join(who)} -- which submission "
                      f"is this data from?"))
        return EXIT_DECIDE
    emit(dict(out, status="none", looked_up=True, ask=ask_none))
    return EXIT_DECIDE


# -------------------------------------------------------------------------- stage --
INSTITUTION_STOPWORDS = {
    "university", "univ", "uni", "of", "the", "college", "inc", "llc", "ltd", "lab", "labs",
    "laboratory", "laboratories", "uc", "and", "at", "for", "school", "center", "centre",
    "institute", "department", "dept", "corp", "corporation", "company", "research",
    "medical", "medicine", "health", "sciences", "science",
}
INSTITUTION_ALIASES = (
    ("university of california san francisco", "ucsf"), ("uc san francisco", "ucsf"),
    ("university of california santa barbara", "ucsb"), ("uc santa barbara", "ucsb"),
    ("university of california berkeley", "berkeley"), ("uc berkeley", "berkeley"),
    ("university of california santa cruz", "ucsc"), ("uc santa cruz", "ucsc"),
    ("university of california san diego", "ucsd"), ("uc san diego", "ucsd"),
    ("university of california los angeles", "ucla"), ("uc los angeles", "ucla"),
    ("university of california irvine", "uci"), ("uc irvine", "uci"),
    ("university of california riverside", "ucr"), ("uc riverside", "ucr"),
    ("university of california merced", "ucm"), ("uc merced", "ucm"),
)
# Never enough on their own to identify a surname: "de la Cruz" must not match UC_Santa_Cruz.
SURNAME_PARTICLES = {"de", "la", "le", "van", "von", "der", "den", "da", "di", "du", "del",
                     "della", "dos", "das", "st", "ter", "ten"}
# Words a group folder carries besides a person's name.
GENERIC_FOLDER_WORDS = {"lab", "labs", "laboratory", "group", "grp", "the", "and", "core", "service",
                        "services", "proteomics", "project", "projects", "data", "ucd", "ucdavis",
                        "davis", "uc", "prot", "samples"}


def name_tokens(s, min_len: int = 3) -> list:
    return [t for t in re.findall(r"[a-z0-9]+", _fold(s).lower()) if len(t) >= min_len]


def institution_tokens(s) -> set:
    text = " ".join(re.findall(r"[a-z0-9]+", _fold(s).lower()))
    for long, short in INSTITUTION_ALIASES:
        text = re.sub(r"\b" + re.escape(long) + r"\b", short, text)
    return {t for t in text.split() if len(t) >= 3 and t not in INSTITUTION_STOPWORDS}


def safe_component(s, fallback: str = "") -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", _fold(s)).strip("_") or fallback


def _subdirs(d: str) -> list:
    try:
        return sorted(os.path.join(d, n) for n in os.listdir(d)
                      if not n.startswith(".") and os.path.isdir(os.path.join(d, n)))
    except OSError:
        return []


def resolve_group(summary: dict, sroot: str) -> dict:
    """Find the PI's existing folder in the service directory, or propose a new one.
    Exactly one clean candidate is reused, none proposes a folder, and anything else -- several
    candidates, or one that names somebody else -- asks a human."""
    campus = summary.get("campus")
    cdir = os.path.join(sroot, campus)
    pi = summary.get("pi") or {}
    last, first, inst = _s(pi.get("last")), _s(pi.get("first")), _s(pi.get("institution"))
    parts = re.findall(r"[a-z0-9]+", _fold(last).lower())
    core = [t for t in parts if t not in SURNAME_PARTICLES]
    if not core:
        return {"mode": "unresolvable", "reason": "the submission has no usable PI last name",
                "candidates": []}
    joined = "".join(parts)
    particles = len(core) < len(parts)
    min_len = min(3, min(len(t) for t in core))       # Li, Wu: match the whole short name
    first_tokens = set(re.findall(r"[a-z0-9]+", _fold(first).lower()))
    itoks = institution_tokens(inst)

    def pi_match(path):
        base = os.path.basename(path)
        if particles:                                  # the whole surname, particles included
            return joined in normkey(base)
        return set(core) <= set(name_tokens(base, min_len))

    def contradiction(path, parent_inst_check):
        base = os.path.basename(path)
        extra = {t for t in name_tokens(base, 3)
                 if not t.isdigit() and t not in joined and t not in GENERIC_FOLDER_WORDS
                 and t not in itoks and t not in SURNAME_PARTICLES}
        if extra and not (extra & first_tokens):
            return (f"'{base}' also names {', '.join(sorted(extra))}, which is not the PI's first "
                    f"name ({first or 'unknown'}) -- possibly a different person with the same surname")
        if parent_inst_check:
            parent = os.path.basename(os.path.dirname(path))
            if not (institution_tokens(parent) & itoks):
                return f"'{parent}' does not look like the PI's institution ({inst or 'unknown'})"
        return None

    depth1 = _subdirs(cdir)
    pi_last_safe = safe_component(last, "PI")
    if campus == "on_campus":
        found = [(d, contradiction(d, False)) for d in depth1 if pi_match(d)]
    else:
        depth2 = [c for d in depth1 for c in _subdirs(d)]
        found = [(d, contradiction(d, False)) for d in depth1 if pi_match(d)]
        found += [(d, contradiction(d, True)) for d in depth2 if pi_match(d)]
    cands = [d for d, _ in found]
    reasons = {d: r for d, r in found if r}
    if len(found) == 1 and not reasons:
        return {"mode": "existing", "group_dir": cands[0], "basis": "PI last name", "candidates": cands}
    if found:
        return {"mode": "ambiguous", "basis": "PI last name", "candidates": cands,
                "reason": "; ".join(reasons.values()) or "several folders carry the PI's surname",
                "contradictions": reasons}
    if campus == "on_campus":
        return {"mode": "new", "group_dir": os.path.join(cdir, pi_last_safe), "basis": "PI last name",
                "candidates": []}
    # An institution folder must carry EXACTLY the institution's distinctive words:
    # "University of Washington" is not "Washington_State_Univ".
    icands = [d for d in depth1 if itoks and institution_tokens(os.path.basename(d)) == itoks]
    if len(icands) == 1:
        return {"mode": "new", "group_dir": os.path.join(icands[0], pi_last_safe),
                "basis": "institution folder + PI last name", "candidates": icands}
    if len(icands) > 1:
        return {"mode": "ambiguous", "basis": "institution", "candidates": icands,
                "reason": "several folders match the institution"}
    inst_safe = safe_component(inst)
    if not inst_safe:
        return {"mode": "unresolvable", "reason": "off-campus submission with no institution",
                "candidates": []}
    return {"mode": "new", "group_dir": os.path.join(cdir, inst_safe, pi_last_safe),
            "basis": "new institution folder + PI last name", "candidates": []}


def resolve_explicit_service_dir(arg: str, sroot: str, campus: str) -> str:
    p = os.path.expanduser(arg)
    if not os.path.isabs(p):
        first = p.replace("\\", "/").split("/")[0]
        p = os.path.join(sroot if first in ("on_campus", "off_campus") else os.path.join(sroot, campus), p)
    p = os.path.abspath(p)
    rel = os.path.relpath(os.path.realpath(p), os.path.realpath(sroot))
    if rel.startswith("..") or os.path.isabs(rel) or len(rel.split(os.sep)) < 2:
        raise Stop(EXIT_DECIDE, {"error": f"--service-dir must be a group folder under "
                                          f"{sroot}/<on_campus|off_campus>/, got {p}"})
    return p


def _chmod_add(path: str, bits: int, warnings: list) -> None:
    try:
        os.chmod(path, stat.S_IMODE(os.lstat(path).st_mode) | bits)
    except OSError as e:
        warnings.append(f"could not set permissions on {path}: {e.strerror}")


def make_shared_dirs(path: str, warnings: list) -> None:
    """mkdir -p, then g+rwxs,o+rx on the directories THIS call created: the service tree is
    shared by the whole Core group (existing dirs are 2775), and a 0700 folder -- what umask
    077 produces -- locks the next staff member out."""
    created, p = [], os.path.abspath(path)
    while not os.path.lexists(p):
        created.append(p)
        parent = os.path.dirname(p)
        if parent == p:
            break
        p = parent
    os.makedirs(path, exist_ok=True)
    for d in reversed(created):
        _chmod_add(d, 0o2775, warnings)


def read_file_list(path: str) -> list:
    try:
        with open(path) as fh:
            return [ln.strip() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    except OSError as e:
        raise Stop(EXIT_DECIDE, {"error": f"cannot read {path}: {e.strerror}"})


def _md_cell(v) -> str:
    return _s(v).replace("|", "\\|").replace("\n", " ") or "—"


def submission_md(s: dict, label: str, links: list, sample_rows_: list, work_dir: str,
                  share_dir: str) -> str:
    pi, sub = s.get("pi") or {}, s.get("submitter") or {}
    url = _s(s.get("url"))
    org = (s.get("organism_as_submitted") or {}).get("value")
    L = [GENERATED_MARK, "",
         f"# {label} — {pi.get('name') or 'PI unknown'} ({pi.get('institution') or 'institution unknown'})",
         "",
         f"_Staged by core_submission.py on {dt.date.today().isoformat()}. Staff-facing: this file "
         f"stays in the service directory and is never delivered._", "",
         "| | |", "|---|---|",
         f"| CoreOmics | {('[' + label + '](' + url + ')') if url else _md_cell(label)} |",
         f"| Status | {_md_cell(s.get('status'))} |",
         f"| PI | {_md_cell(pi.get('name'))} — {_md_cell(pi.get('institution'))}"
         f"{(', ' + _s(pi.get('department'))) if pi.get('department') else ''} |",
         f"| Submitter | {_md_cell(sub.get('name'))} |",
         f"| Submitted | {_md_cell(s.get('submitted_date'))} |",
         f"| Organism (as submitted — UNCONFIRMED) | {_md_cell(org)} |",
         f"| Experiment type | {_md_cell(', '.join(s.get('experiment_types') or []))} |",
         f"| Instrument requested | {_md_cell(s.get('instrument_wanted'))} |",
         f"| Analysis requested | {_md_cell(s.get('data_analysis'))} |",
         f"| Samples | {s.get('n_samples', '—')} |",
         f"| Raw files linked | {len(links)} (in `raw/`) |",
         f"| Compute work dir | `{work_dir}` |",
         f"| Bioshare share dir | `{share_dir}` |", ""]
    for title, key in (("Description", "description"), ("Sample prep", "sample_prep")):
        if _s(s.get(key)):
            L += [f"## {title}", "", _s(s.get(key)), ""]
    L += ["## Samples", "", "| unique_id | sample name | condition | raw file | acquired | status |",
          "|---|---|---|---|---|---|"]
    if sample_rows_:
        for r in sample_rows_:
            L.append(f"| {_md_cell(r.get('unique_id'))} | {_md_cell(r.get('sample_name'))} | "
                     f"{_md_cell(r.get('condition_name'))} | "
                     f"{_md_cell(os.path.basename(_s(r.get('file'))))} | {_md_cell(r.get('acquired'))} | "
                     f"{_md_cell(r.get('status'))} |")
    else:
        for x in s.get("samples") or []:
            L.append(f"| {_md_cell(x.get('unique_id'))} | {_md_cell(x.get('sample_name'))} | "
                     f"{_md_cell(x.get('condition_name'))} | — | — | — |")
    return "\n".join(L) + "\n"


def locate_gate_check(files_path: str, warnings: list) -> dict:
    """stage --apply refuses a file list whose locate run hard-failed: files.txt is written
    even then (it is a proposal), and nothing else stops it being staged."""
    lj = os.path.join(os.path.dirname(os.path.abspath(files_path)), "locate.json")
    if not os.path.exists(lj):
        warnings.append(f"no locate.json next to {files_path}; its gates could not be checked")
        return {"checked": False}
    loc = load_json_file(lj, "locate.json")
    info = {"checked": True, "locate_json": lj, "mode": loc.get("mode"), "hard_fail": bool(loc.get("hard_fail")),
            "accepted": loc.get("accepted")}
    if loc.get("hard_fail") and loc.get("mode") != "files_from":
        raise Stop(EXIT_DECIDE, {
            "error": "locate hard-failed for this file list; staging it would search the wrong or partial files",
            "failing_gates": [g for g in loc.get("gates") or [] if g.get("status") == "FAIL"],
            "hint": "resolve with staff: locate --files-from <curated list>, or re-run locate with "
                    "--allow-partial / --accept-ambiguous as they decide"})
    return info


def cmd_stage(a) -> int:
    s = load_summary(a.summary)
    label = submission_label(s)
    warnings = []
    files = [os.path.abspath(os.path.expanduser(f.rstrip("/"))) for f in read_file_list(a.files)]
    if not files:
        raise Stop(EXIT_DECIDE, {"error": f"{a.files} lists no files"})
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        raise Stop(EXIT_DECIDE, {"error": f"{len(missing)} listed file(s) do not exist", "missing": missing[:20]})
    names = Counter(os.path.basename(f) for f in files)
    clash = sorted(n for n, c in names.items() if c > 1)
    if clash:
        raise Stop(EXIT_DECIDE, {"error": "two listed files share a basename; raw/ cannot hold both",
                                 "basenames": clash})
    campus = s.get("campus")
    if campus not in ("on_campus", "off_campus"):
        raise Stop(EXIT_DECIDE, {"error": f"summary campus is {campus!r}"})
    sroot = service_root()
    if not os.path.isdir(os.path.join(sroot, campus)):
        raise Stop(EXIT_UNREACHABLE, {"error": f"service directory not reachable: {os.path.join(sroot, campus)}",
                                      "hint": "stage runs ON HIVE (hive_exec.sh)"})
    gate_info = locate_gate_check(a.files, warnings) if a.apply else {"checked": False, "note": "checked on --apply"}

    if a.service_dir:
        group = resolve_explicit_service_dir(a.service_dir, sroot, campus)
        resolution = {"mode": "explicit", "group_dir": group,
                      "exists": os.path.isdir(group), "candidates": []}
        if not is_under(group, os.path.join(sroot, campus)):
            warnings.append(f"--service-dir is outside {campus}/ although CoreOmics says {campus}")
    else:
        resolution = resolve_group(s, sroot)
        if resolution["mode"] in ("ambiguous", "unresolvable"):
            raise Stop(EXIT_DECIDE, {
                "error": f"cannot pick the service folder for {label}: {resolution.get('reason') or 'several candidates'}",
                "group_resolution": resolution,
                "hint": "ask staff which folder, then re-run with --service-dir <folder>"})
        group = resolution["group_dir"]

    project = os.path.join(group, label)
    marker = os.path.join(project, MARKER)
    prior = None
    if os.path.exists(marker):
        prior = load_json_file(marker, MARKER)
        if _s(prior.get("id")) != _s(s.get("id")):
            raise Stop(EXIT_DECIDE, {"error": f"{project} already belongs to submission "
                                              f"{prior.get('internal_id') or prior.get('id')}",
                                     "hint": "pick another folder with --service-dir"})
    elif os.path.isdir(project) and os.listdir(project):
        warnings.append(f"{project} already exists without {MARKER} (staged by hand?); adopting it -- "
                        f"existing files are never overwritten")

    froot = flinders_root()
    raw_dir = os.path.join(project, "raw")
    rel_root = os.path.relpath(os.path.realpath(group), os.path.realpath(sroot))
    work_dir = os.path.join(work_root(), rel_root, label)
    share = server_share_dir(s)

    plan = []
    for f in files:
        target = os.path.realpath(f)
        dst = os.path.join(raw_dir, os.path.basename(f))
        if is_under(target, froot):
            link = os.path.relpath(target, os.path.realpath(raw_dir))
        else:
            link = target
            warnings.append(f"{os.path.basename(f)}: target is outside the Flinders root, so the link "
                            f"is absolute -- not visible over SMB/Bioshare")
        if os.path.islink(dst):
            action = "existing" if os.readlink(dst) == link else "replace"
        elif os.path.lexists(dst):
            action = "conflict"
        else:
            action = "create"
        plan.append({"name": os.path.basename(f), "target": target, "link": link,
                     "relative": not os.path.isabs(link), "action": action})
    conflicts = [p["name"] for p in plan if p["action"] == "conflict"]
    if conflicts:
        warnings.append(f"{len(conflicts)} name(s) in raw/ are real files, not links; left untouched")

    md_path = os.path.join(project, "SUBMISSION.md")
    if os.path.lexists(md_path):
        try:
            with open(md_path, errors="replace") as fh:
                ours = fh.readline().strip() == GENERATED_MARK
        except OSError:
            ours = False
        if not ours or os.path.islink(md_path):
            md_path = os.path.join(project, "SUBMISSION.core.md")
            warnings.append("SUBMISSION.md was written by hand; it is left alone and the generated "
                            "record goes to SUBMISSION.core.md")

    out = {"applied": bool(a.apply), "internal_id": s.get("internal_id"), "service_dir": project,
           "group_dir": group, "group_resolution": resolution, "work_dir": work_dir, "share_dir": share,
           "links_created": 0, "links_existing": sum(p["action"] == "existing" for p in plan),
           "links_replaced": 0, "conflicts": conflicts, "locate_gates": gate_info, "warnings": warnings}
    if not a.apply:
        out.update(would_create=sum(p["action"] == "create" for p in plan),
                   would_replace=sum(p["action"] == "replace" for p in plan), plan=plan,
                   submission_md=md_path,
                   next_step="show staff this folder together with the file list in ONE confirmation; "
                             "then re-run with --apply")
        emit(out)
        return EXIT_OK

    make_shared_dirs(raw_dir, warnings)
    make_shared_dirs(work_dir, warnings)
    broken = []
    for p in plan:
        dst = os.path.join(raw_dir, p["name"])
        if p["action"] == "conflict":
            continue
        if p["action"] == "replace":
            os.unlink(dst)
        if p["action"] in ("create", "replace"):
            os.symlink(p["link"], dst)
            out["links_created" if p["action"] == "create" else "links_replaced"] += 1
        if not (os.path.exists(dst) and os.path.realpath(dst) == p["target"]):
            broken.append(p["name"])

    sample_files = a.sample_files or os.path.join(os.path.dirname(os.path.abspath(a.files)), "sample_files.tsv")
    rows = []
    if os.path.exists(sample_files):
        rows = read_tsv(sample_files)
    else:
        warnings.append(f"no {sample_files}; SUBMISSION.md lists samples without their raw files")
    with open(md_path, "w") as fh:
        fh.write(submission_md(s, label, plan, rows, work_dir, share))
    record = {"schema": SCHEMA, "id": s.get("id"), "internal_id": s.get("internal_id"), "label": label,
              "coreomics_url": s.get("url"), "campus": campus, "service_dir": project, "group_dir": group,
              "work_dir": work_dir, "share_dir": share, "files": [p["target"] for p in plan],
              "links": plan, "created": (prior or {}).get("created") or now_utc(), "updated": now_utc(),
              "skill_version": skill_version()}
    with open(marker, "w") as fh:
        json.dump(record, fh, indent=2)
    # The same record in the work dir lets `deliver` find the service project by walking up
    # from the session, which lives under work_dir.
    with open(os.path.join(work_dir, MARKER), "w") as fh:
        json.dump(record, fh, indent=2)
    out.update(submission_md=md_path, marker=marker, broken_links=broken)
    emit(out)
    if broken:
        note(f"{len(broken)} link(s) do not resolve to their target -- see broken_links")
        return EXIT_DECIDE
    return EXIT_OK


# --------------------------------------------------------------------- conditions --
def cmd_conditions(a) -> int:
    run_name = sibling("collect_conditions").run_name     # File.Name is defined in ONE place
    cc_script = os.path.join(HERE, "collect_conditions.py")
    s = load_summary(a.summary)
    try:
        rows = read_tsv(a.sample_files)
    except OSError as e:
        raise Stop(EXIT_DECIDE, {"error": f"cannot read {a.sample_files}: {e.strerror}"})
    matched = [r for r in rows if r.get("status") == "matched" and _s(r.get("file"))]
    loose = [r for r in rows if r.get("status") == "unassigned" and _s(r.get("file"))]
    if not matched and not loose:
        raise Stop(EXIT_DECIDE, {"error": "no located files in sample_files.tsv -- run `locate` first",
                                 "needs_user_input": True})
    runs = list(dict.fromkeys(run_name(r["file"]) for r in matched + loose))
    bad = [r for r in runs if "," in r]
    if bad:
        raise Stop(EXIT_DECIDE, {"error": "run names containing commas cannot be passed to "
                                          "collect_conditions.py --runs", "runs": bad})
    mapping = {run_name(r["file"]): _s(r.get("condition_name")) for r in matched if _s(r.get("condition_name"))}
    cmd = [sys.executable, cc_script, "--map", a.out,
           "--runs", ",".join(runs), "--mapping-json", json.dumps({"mapping": mapping})]
    p = subprocess.run(cmd, capture_output=True, text=True)
    try:
        cc = json.loads(p.stdout)
    except ValueError:
        cc = None
    if p.returncode != 0 or not isinstance(cc, dict):
        emit({"error": "collect_conditions.py failed", "returncode": p.returncode,
              "stderr": p.stderr[-2000:], "stdout": p.stdout[-2000:]})
        return 1

    per_sample = {}
    for r in matched:
        per_sample.setdefault(_s(r.get("unique_id")), _s(r.get("condition_name")))
    with_cond = {u: c for u, c in per_sample.items() if c}
    groups = Counter(with_cond.values())
    findings, questions = {}, []
    blank = sorted(u for u, c in per_sample.items() if not c)
    if blank:
        findings["blank_conditions"] = blank
        questions.append(f"{len(blank)} sample(s) have no condition in CoreOmics ({', '.join(blank[:12])}"
                         f"{'...' if len(blank) > 12 else ''}). Which group does each belong to?")
    variants = {}
    for g in groups:
        variants.setdefault(normkey(g), []).append(g)
    clashes = {k: sorted(v) for k, v in variants.items() if len(v) > 1}
    if clashes:
        findings["case_variants"] = clashes
        questions.append("These conditions differ only in capitals or punctuation: "
                         + "; ".join(" / ".join(f"'{x}'" for x in v) for v in clashes.values())
                         + ". Is each set one group?")
    if len(groups) == 1:
        only = next(iter(groups))
        findings["single_condition"] = only
        questions.append(f"Every sample's condition is '{only}', so there is nothing to compare. "
                         f"What are the groups?")
    if len(with_cond) >= 2 and len(groups) == len(with_cond):
        findings["all_unique"] = True
        questions.append("Every sample has its own condition -- that usually means sample names were "
                         "typed into the condition column. Which samples are replicates of each other?")
    singles = sorted(g for g, n in groups.items() if n < 2)
    if singles and not findings.get("all_unique"):
        findings["singleton_groups"] = singles
        questions.append(f"These groups have a single sample, so no within-group variance: "
                         f"{', '.join(singles)}. Are there replicates missing, or should groups be merged?")
    elif singles:
        findings["singleton_groups"] = singles
    if cc.get("needs_confirmation"):
        amb = {k: v for k, v in (cc.get("ambiguities") or {}).items() if v}
        findings["collect_conditions_ambiguities"] = amb
        if amb.get("unassigned_runs"):
            questions.append(f"These raw files have no condition yet: {', '.join(amb['unassigned_runs'][:12])}"
                             f"{'...' if len(amb['unassigned_runs']) > 12 else ''}. Which group is each?")
        if amb.get("conflicting_runs"):
            questions.append(f"These raw files matched two groups: {json.dumps(amb['conflicting_runs'])}. "
                             f"Which is right?")
    needs = bool(findings)
    emit({"proposed_csv": cc.get("proposed_csv"), "needs_user_input": needs, "questions": questions,
          "findings": findings, "groups": dict(groups), "n_samples": len(per_sample), "n_runs": len(runs),
          "internal_id": s.get("internal_id"), "collect_conditions": cc,
          "next_step": ("ask the staff member exactly the questions above, then fix the CSV and run "
                        "collect_conditions.py --validate") if needs else
                       "conditions come straight from CoreOmics; confirm the groups in the one compute confirmation"})
    return EXIT_DECIDE if needs else EXIT_OK


# ------------------------------------------------------------------------ deliver --
# No AI_Analysis_Report.docx: Analysis_Report.html is the report of record, and an older
# session's Word copy (Word mangled its figures) is kept there but never handed out.
DELIVER_FILES = (("Analysis_Report.html", True), ("Analysis_Report.pdf", False),
                 ("Analysis_Report.md", False), ("AI_Analysis_Report.md", False), ("methods.md", False), ("methods.docx", False),
                 ("OUTPUT_FILES.md", False), ("AUDIT.md", False), ("SAMPLE_QUALITY.md", False))
DELIVER_DIRS = ("tables", "figures", "reproducibility")
# The optional audio discussion (make_podcast.py): the audio and its transcript only. Its
# chunk cache, script, check/verify logs and podcast.json (consent notes) stay in the session.
# The report's "Listen" card links podcast/<file> relative to Analysis_Report.html.
PODCAST_FILES = ("podcast.m4a", "podcast.wav", "transcript.html")
SEARCH_FILES = ("report.parquet", "report.pg_matrix.tsv", "report.pr_matrix.tsv", "report.gg_matrix.tsv",
                "report.unique_genes_matrix.tsv", "report.stats.tsv", "report.log.txt",
                "search_provenance.json")
EXCLUDED_DIRS = {"xic", "report_xic", "temp", "tmp", ".tmp", "__pycache__", "raw_data"}
EXCLUDED_SUFFIXES = (".quant", ".speclib")


def _excluded_dir(name: str) -> bool:
    n = name.lower()
    return n in EXCLUDED_DIRS or n.endswith("_xic") or n.startswith(("temp_", "tmp_"))


def _plan_file(src: str, name: str, required: bool = False) -> dict:
    item = {"name": name, "src": src, "required": required}
    try:
        st = os.stat(src)
    except OSError as e:
        reason = ("broken symlink" if os.path.islink(src) else
                  "not found in the session" if e.errno == errno.ENOENT else e.strerror)
        return dict(item, status="SKIPPED", reason=reason)
    if not stat.S_ISREG(st.st_mode):
        return dict(item, status="SKIPPED", reason="not a regular file")
    if not os.access(src, os.R_OK):
        return dict(item, status="SKIPPED", reason="unreadable (permission denied)")
    return dict(item, status="OK", bytes=st.st_size)


def plan_delivery(output_dir: str) -> list:
    items = [_plan_file(os.path.join(output_dir, n), n, req) for n, req in DELIVER_FILES]
    for d in DELIVER_DIRS:
        root = os.path.join(output_dir, d)
        if not os.path.isdir(root):
            items.append({"name": d + "/", "status": "SKIPPED", "reason": "not found in the session",
                          "required": False})
            continue

        def unreadable(e, _items=items):
            # os.walk skips a directory it cannot list unless told; that must not vanish silently.
            rel = os.path.relpath(e.filename, output_dir) if e.filename else "?"
            _items.append({"name": rel.rstrip("/") + "/", "status": "SKIPPED",
                           "reason": f"unreadable directory ({e.strerror})", "required": False})

        seen = set()
        for dirpath, dirnames, filenames in os.walk(root, followlinks=True, onerror=unreadable):
            real = os.path.realpath(dirpath)
            if real in seen:                           # a symlink loop
                dirnames[:] = []
                continue
            seen.add(real)
            rel_dir = os.path.relpath(dirpath, output_dir)
            keep = []
            for n in sorted(dirnames):
                if _excluded_dir(n):
                    items.append({"name": os.path.join(rel_dir, n) + "/", "status": "SKIPPED",
                                  "reason": "excluded: chromatograms/temp data are not a deliverable",
                                  "required": False})
                else:
                    keep.append(n)
            dirnames[:] = keep
            for fn in sorted(filenames):
                rel = os.path.join(rel_dir, fn)
                if fn.lower().endswith(EXCLUDED_SUFFIXES):
                    items.append({"name": rel, "status": "SKIPPED",
                                  "reason": "excluded: engine intermediate", "required": False})
                    continue
                items.append(_plan_file(os.path.join(dirpath, fn), rel))
    pdir = os.path.join(output_dir, "podcast")
    if not os.path.isdir(pdir):
        items.append({"name": "podcast/", "status": "SKIPPED", "required": False,
                      "reason": "none was made (optional)"})
    else:
        audio = [fn for fn in PODCAST_FILES if os.path.exists(os.path.join(pdir, fn))]
        if not any(fn != "transcript.html" for fn in audio):
            items.append({"name": "podcast/", "status": "SKIPPED", "required": False,
                          "reason": "no rendered audio in output/podcast"})
        items += [_plan_file(os.path.join(pdir, fn), "podcast/" + fn) for fn in audio]
    for fn in SEARCH_FILES:
        items.append(_plan_file(os.path.join(output_dir, "search", fn), "search/" + fn))
    return items


def find_service_record(summary: dict, session_dir, explicit) -> tuple:
    """(record, where, foreign) for the staged service project. `foreign` names a DIFFERENT
    submission whose stage record sits above the session -- a session that belongs to it."""
    sid = _s(summary.get("id"))

    def read(path):
        try:
            with open(path) as fh:
                r = json.load(fh)
            return r if isinstance(r, dict) else None
        except (OSError, ValueError):
            return None

    foreign = None
    if session_dir:
        d = os.path.abspath(session_dir)
        while True:
            r = read(os.path.join(d, MARKER))
            if r and _s(r.get("id")) != sid:
                foreign = foreign or (r.get("internal_id") or r.get("id"))
            elif r and r.get("service_dir") and os.path.exists(os.path.join(r["service_dir"], MARKER)):
                return r, r["service_dir"], foreign
            parent = os.path.dirname(d)
            if parent == d:
                break
            d = parent
    if explicit:
        r = read(os.path.join(explicit, MARKER))
        if r and _s(r.get("id")) == sid:
            return r, explicit, foreign
        return None, f"no {MARKER} for this submission in {explicit}", foreign
    label = submission_label(summary)
    hits = []
    for campus in ("on_campus", "off_campus"):
        for depth in ("*", os.path.join("*", "*")):
            for m in glob.glob(os.path.join(service_root(), campus, depth, label, MARKER)):
                r = read(m)
                if r and _s(r.get("id")) == sid:
                    hits.append(os.path.dirname(m))
    if len(hits) == 1:
        return read(os.path.join(hits[0], MARKER)), hits[0], foreign
    if len(hits) > 1:
        return None, "several staged service projects: " + ", ".join(hits) + " (pass --service-project)", foreign
    return None, "no staged service project found (run `stage --apply` first, or pass --service-project)", foreign


def session_ownership(session_dir: str, record: dict) -> tuple:
    """A session belongs to a submission when it lives under the staged work dir, or when it
    searched exactly the staged raw files. Anything else is another submission's results."""
    wd = record.get("work_dir")
    if wd and is_under(session_dir, wd):
        return True, "the session is under the staged work dir"
    # session.read_raw_list: the one reader, which also reads a list skill 2.7 wrote in cp1252
    try:
        raws = {os.path.realpath(r) for r in sibling("session").read_raw_list(session_dir)}
    except OSError:
        raws = set()
    staged = {os.path.realpath(f) for f in record.get("files") or []}
    if raws and raws == staged:
        return True, "the session's input/raw_files.txt matches the staged files"
    return False, (f"the session is not under the staged work dir ({wd}) and its raw_files.txt "
                   f"({len(raws)} file(s)) does not match the {len(staged)} staged file(s)")


def symlinked_components(path: str, root: str) -> list:
    """Why writing at `path` would escape: a component below `root` that is a symlink, or a
    real location outside `root`. Empty when the path is safe to create or write."""
    ap, ar = os.path.abspath(path), os.path.abspath(root)
    rel = os.path.relpath(ap, ar)
    if rel == os.pardir or rel.startswith(os.pardir + os.sep):
        return [f"{ap} is not under the Flinders root {ar}"]
    problems, cur = [], ar
    for part in ([] if rel == os.curdir else rel.split(os.sep)):
        cur = os.path.join(cur, part)
        if os.path.islink(cur):
            problems.append(f"{cur} is a symlink (-> {os.readlink(cur)})")
            break
        if not os.path.exists(cur):
            break
    if not is_under(cur if os.path.lexists(cur) else os.path.dirname(cur), ar):
        problems.append(f"{cur} resolves outside the Flinders root")
    return problems


def mkdirs_nofollow(path: str, root: str, warnings: list) -> None:
    """mkdir -p that refuses to pass through a symlink: makedirs would follow one planted in
    the share and create the delivery somewhere else entirely."""
    problems = symlinked_components(path, root)
    if problems:
        raise UnsafePath("; ".join(problems))
    ar = os.path.abspath(root)
    rel = os.path.relpath(os.path.abspath(path), ar)
    cur = ar
    for part in ([] if rel == os.curdir else rel.split(os.sep)):
        cur = os.path.join(cur, part)
        if os.path.islink(cur):
            raise UnsafePath(f"{cur} is a symlink")
        if not os.path.lexists(cur):
            os.mkdir(cur)
            _chmod_add(cur, 0o2775, warnings)
        elif not os.path.isdir(cur):
            raise UnsafePath(f"{cur} is not a directory")


def _open_nofollow(dst: str) -> int:
    if os.path.islink(dst):
        os.unlink(dst)                                  # never write THROUGH a planted link
    return os.open(dst, os.O_WRONLY | os.O_CREAT | os.O_TRUNC | getattr(os, "O_NOFOLLOW", 0), 0o664)


def copy_nofollow(src: str, dst: str, root: str, warnings: list) -> None:
    """Copy the BYTES of src (following a symlinked source) to dst without following any
    symlink on the destination side, then make it group- and world-readable: Bioshare reads
    as another user, and copy2 would have kept a 0600 source mode."""
    mkdirs_nofollow(os.path.dirname(dst), root, warnings)
    st = os.stat(src)
    with open(src, "rb") as inp, os.fdopen(_open_nofollow(dst), "wb") as out:
        shutil.copyfileobj(inp, out, 1 << 20)
    os.utime(dst, (st.st_atime, st.st_mtime))
    _chmod_add(dst, 0o664 | (0o111 if st.st_mode & 0o111 else 0), warnings)


def write_text_nofollow(dst: str, text: str, root: str, warnings: list) -> None:
    mkdirs_nofollow(os.path.dirname(dst), root, warnings)
    with os.fdopen(_open_nofollow(dst), "w", encoding="utf-8") as out:
        out.write(text)
    _chmod_add(dst, 0o664, warnings)


def tables_description(delivered: set) -> str:
    """Name only what tables/ really holds -- a README promising an expression matrix that is
    not there sends the collaborator looking for it."""
    parts = []
    if any(x.startswith("tables/DE_") for x in delivered):
        parts.append("one DE_*.csv per comparison")
    for name, text in (("tables/Expression_Matrix.csv", "the expression matrix"),
                       ("tables/methods.txt", "methods.txt"),
                       ("tables/reproducibility_log.R", "reproducibility_log.R (the analysis as R code)")):
        if name in delivered:
            parts.append(text)
    return "result tables" + (": " + ", ".join(parts) if parts else "")


def human_size(n: int) -> str:
    for unit in ("bytes", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:,} bytes" if unit == "bytes" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def _submission_details(summary: dict) -> str:
    extra = []
    if summary.get("submitted_date"):
        extra.append(f"submitted {summary['submitted_date']}")
    if summary.get("n_samples"):
        extra.append(f"{summary['n_samples']} samples")
    if summary.get("experiment_types"):
        extra.append("experiment type as submitted: " + ", ".join(summary["experiment_types"]))
    return f" ({'; '.join(extra)})" if extra else ""


def build_readme(summary: dict, delivered: set, raw_names: list, mode: str = "analysis") -> str:
    """Collaborator-facing. Every claim comes from the submission record or from a file that
    was actually delivered -- no internal paths, no email addresses, no numbers that are not
    in the inputs, and no "we compared your groups" without a comparison table beside it."""
    label = submission_label(summary)
    pi_name = (summary.get("pi") or {}).get("name")
    fixed = ("MANIFEST.txt", "checksums.sha256")
    if mode == "raw-only":
        L = [f"# {label} — raw data from the UC Davis Proteomics Core", ""]
        if pi_name:
            L += [f"Prepared for {pi_name}.", ""]
        why = (" The submission asked for raw data only, so no database search or statistical "
               "analysis was run." if summary.get("raw_data_only") else
               " This delivery holds the raw files only; it contains no search or statistical results.")
        L += [f"This share holds the raw instrument files for CoreOmics submission {label}"
              f"{_submission_details(summary)}.{why}", "", "## What is here", "",
              "| item | what it is |", "|---|---|"]
        if raw_names:
            L.append(f"| `../raw/` | {len(raw_names)} raw instrument file(s), in the `raw` folder next to this one |")
        for name, text in (("methods.md", "a Methods section describing the acquisition, with the instrument acknowledgment"),
                           ("methods.docx", "the Methods section as a Word document"),
                           ("MANIFEST.txt", "what was included and, for anything not included, why"),
                           ("checksums.sha256", "SHA-256 checksums of the files in this folder "
                                                "(`shasum -a 256 -c checksums.sha256`)")):
            if name in delivered or name in fixed:
                L.append(f"| `{name}` | {text} |")
        exts = {os.path.splitext(n.lower())[1] for n in raw_names}
        if exts & {".d", ".raw"}:
            L += ["", "## Opening the files", ""]
            if ".d" in exts:
                L.append("- `.d` entries are Bruker timsTOF data **folders** — download each one as a whole folder.")
            if ".raw" in exts:
                L.append("- `.raw` files are Thermo instrument files.")
    else:
        L = [f"# {label} — results from the UC Davis Proteomics Core", ""]
        if pi_name:
            L += [f"Prepared for {pi_name}.", ""]
        claims = []
        if any(x.startswith("search/") for x in delivered):
            claims.append("The mass-spectrometry data were searched to identify and quantify proteins; "
                          "the search output is in `search/`.")
        if any(x.startswith("tables/DE_") for x in delivered):
            claims.append("Protein abundances were compared between sample groups; there is one "
                          "table per comparison in `tables/`.")
        L += [f"This folder holds what the UC Davis Proteomics Core prepared for CoreOmics submission "
              f"{label}{_submission_details(summary)}." + ("" if not claims else " " + " ".join(claims)), ""]
        if "Analysis_Report.html" in delivered:
            L += ["## Start here", "",
                  "**Open `Analysis_Report.html` first.** It is the whole report in one file — every figure "
                  "is embedded, so it opens by double-clicking in any web browser, with no internet "
                  "connection and nothing to install. The quality-control panels come before the results "
                  "on purpose: they show how much weight the results can carry.", ""]
        L += ["## What is in this folder", "", "| item | what it is |", "|---|---|"]
        desc = (
            ("Analysis_Report.html", "the full report: quality control, figures and interpretation (open first)"),
            ("Analysis_Report.pdf", "the same report with its figures as a PDF, for printing or NotebookLM"),
            ("Analysis_Report.md", "the same report as plain text, with each figure's caption and numbers "
                                   "written out, for NotebookLM or other AI notebooks"),
            ("AI_Analysis_Report.md", "the interpretation alone, as plain text"),
            ("methods.md", "a Methods section ready to adapt for a manuscript, with the instrument acknowledgment"),
            ("methods.docx", "the Methods section as a Word document"),
            ("tables/", tables_description(delivered)),
            ("figures/", "every figure as a separate image file"),
            ("podcast/", "an optional AI-generated audio discussion of these results (synthetic "
                         "voices) and its transcript; the report links it near the top"),
            ("reproducibility/", "the pinned record of the run: software versions, parameters, checksums, "
                                 "and reproduce.sh"),
            ("search/", "the search engine's own output (report.parquet and the matrices delivered)"),
            ("AUDIT.md", "automated checks of the experimental design and results, with caveats"),
            ("SAMPLE_QUALITY.md", "sample-quality checks (e.g. contamination that can mimic biology)"),
            ("OUTPUT_FILES.md", "a catalogue of every file and what it is"),
            ("MANIFEST.txt", "what was included in this folder and, for anything not included, why"),
            ("checksums.sha256", "SHA-256 checksums to verify the download (`shasum -a 256 -c checksums.sha256`)"),
        )
        for name, text in desc:
            present = (name in delivered or name in fixed or
                       (name.endswith("/") and any(x.startswith(name) for x in delivered)))
            if present:
                L.append(f"| `{name}` | {text} |")
        if raw_names:
            L.append("| `../raw/` | your raw instrument files, in the `raw` folder next to this one |")
        methods = []
        if "methods.md" in delivered:
            methods.append("- **Methods text:** `methods.md`" + (" (and `methods.docx`)" if "methods.docx" in delivered
                                                                 else "")
                           + ", generated from the instrument files and the settings that actually ran.")
        if "tables/methods.txt" in delivered:
            methods.append("- **Statistical methods paragraph:** `tables/methods.txt`.")
        if "tables/reproducibility_log.R" in delivered:
            methods.append("- **The analysis as R code:** `tables/reproducibility_log.R` — every step with every "
                           "value written out; it re-runs with R and the limpa/limma packages.")
        if "reproducibility/REPRODUCE.md" in delivered:
            methods.append("- **Re-running everything, search included:** `reproducibility/REPRODUCE.md`.")
        L += ["", "## Methods and code", ""] + (methods or ["- Ask the Core for the methods text for this analysis."])
    L += ["", "## Acknowledging the Core", "",
          "If these data appear in a publication, poster or talk, please acknowledge the UC Davis "
          "Proteomics Core" + (" and the instrument grant named at the end of `methods.md`."
                               if "methods.md" in delivered else
                               ", and ask us for the instrument grant acknowledgment that applies."),
          "", "## Questions", "", "Contact the UC Davis Proteomics Core.", ""]
    return "\n".join(L)


def sha256_of(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_delivery(delivery_dir: str, share_dir: str, froot: str, expected=None, own_raw=None):
    """What a collaborator would actually get. Returns (violations, notes).

    Strict about what THIS run made: nothing inside the delivery folder may be a symlink, every
    file there must be one this run delivered (no stale leftovers) and readable by group and
    other (Bioshare reads as another user), and every raw/ link this run created (`own_raw`,
    by name) must be relative and resolve under the Flinders root.

    Links already in the share that this run did not create are reported as NOTES and left
    alone. PROT_0793's share/search -> /quobyte links are the model case: they are valid, and
    stopped serving only because of a Bioshare server setting (fixed 2026-09-16), so treating
    them as a failure would block every later delivery into that share for something this run
    neither made nor can fix."""
    own_raw = set(own_raw or ())
    bad = [f"unsafe delivery path: {p}" for p in symlinked_components(delivery_dir, froot)]
    notes = []
    share_raw = os.path.abspath(os.path.join(share_dir, "raw"))
    delivery_abs = os.path.abspath(delivery_dir)

    def unreadable(e):
        bad.append(f"cannot read {e.filename}: {e.strerror}")

    if os.path.isdir(share_dir) and not os.path.islink(share_dir):
        for dirpath, dirnames, filenames in os.walk(share_dir, onerror=unreadable):
            for n in dirnames + filenames:
                p = os.path.join(dirpath, n)
                rel = os.path.relpath(p, share_dir)
                if os.path.islink(p):
                    t = os.readlink(p)
                    inside_delivery = os.path.abspath(p) == delivery_abs or \
                        os.path.abspath(p).startswith(delivery_abs + os.sep)
                    if inside_delivery:
                        bad.append(f"symlink inside this delivery (deliverables must be real files): {rel} -> {t}")
                    elif os.path.abspath(dirpath) == share_raw and n in own_raw:
                        if os.path.isabs(t):
                            bad.append(f"absolute link in raw/ (invisible over SMB): {rel} -> {t}")
                        elif not os.path.exists(p):
                            bad.append(f"broken link in raw/: {rel} -> {t}")
                        elif not is_under(p, froot):
                            bad.append(f"raw/ link leaves the Flinders root: {rel}")
                    elif not os.path.exists(p):
                        notes.append(f"existing link left untouched, does not resolve on this machine "
                                     f"(its target may not be mounted on this machine, e.g. /quobyte): {rel} -> {t}")
                    elif os.path.isabs(t) and not is_under(p, froot):
                        notes.append(f"existing link left untouched, points outside the Flinders root -- it "
                                     f"works only where the file server and Bioshare can reach its target: {rel} -> {t}")
                    elif os.path.isabs(t):
                        notes.append(f"existing absolute link left untouched (invisible over SMB): {rel} -> {t}")
                elif not (os.path.isfile(p) or os.path.isdir(p)):
                    bad.append(f"not a regular file or directory: {rel}")
    if expected is not None and os.path.isdir(delivery_dir) and not os.path.islink(delivery_dir):
        if stat.S_IMODE(os.lstat(delivery_dir).st_mode) & 0o055 != 0o055:
            bad.append("the delivery folder is not group- and other-readable")
        for dirpath, dirnames, filenames in os.walk(delivery_dir):
            for n in filenames:
                p = os.path.join(dirpath, n)
                rel = os.path.relpath(p, delivery_dir)
                if os.path.islink(p) or rel == DELIVERY_MARKER:
                    continue
                if rel not in expected:
                    bad.append(f"file not in this delivery's plan (stale or foreign): {rel}")
                if stat.S_IMODE(os.lstat(p).st_mode) & 0o044 != 0o044:
                    bad.append(f"not group- and other-readable: {rel}")
            for n in dirnames:
                p = os.path.join(dirpath, n)
                if not os.path.islink(p) and stat.S_IMODE(os.lstat(p).st_mode) & 0o055 != 0o055:
                    bad.append(f"folder not group- and other-readable: {os.path.relpath(p, delivery_dir)}")
    return bad, notes


def write_deliver_job(a, session_dir, where_dir: str, label: str, total: int, warnings: list) -> str:
    path = os.path.join(where_dir, "deliver_job.sh")
    logs = os.path.join(where_dir, "logs")
    part, acct, qos = sibling("run_search").slurm_queue(peak_cpus=1, preemptible_ok=True)
    argv = [sys.executable if os.path.isabs(sys.executable) else "python3", os.path.abspath(__file__),
            "deliver", "--summary", os.path.abspath(a.summary), "--mode", a.mode,
            "--include-raw", a.include_raw, "--apply", "--no-size-guard"]
    if session_dir:
        argv += ["--session", session_dir]
    if a.label:
        argv += ["--label", a.label]
    if a.service_project:
        argv += ["--service-project", os.path.abspath(a.service_project)]
    if a.force:
        argv += ["--force"]
    L = ["#!/bin/bash", f"#SBATCH --job-name=deliver_{label}",
         f"#SBATCH --output={os.path.join(logs if os.path.isdir(logs) else where_dir, 'deliver_%j.log')}",
         "#SBATCH --cpus-per-task=1", "#SBATCH --mem=2G", "#SBATCH --time=24:00:00"]
    if part:
        L.append(f"#SBATCH --partition={part}")
    if acct:
        L.append(f"#SBATCH --account={acct}")
    if qos:
        L.append(f"#SBATCH --qos={qos}")
    if part == "low":
        L.append("#SBATCH --requeue")
    L += [f"# Written by core_submission.py deliver: {total / 1024 ** 3:.1f} GB is too much to copy on",
          "# the login node. Re-runs the SAME delivery with --apply. A preempted, requeued job resumes",
          "# the folder it started (it carries .core_delivery.json); verification still runs at the end.",
          "set -euo pipefail"]
    for var in ("CORE_FLINDERS_ROOT", "CORE_WORK_ROOT"):
        if os.environ.get(var):
            L.append(f"export {var}={shlex.quote(os.environ[var])}")
    L.append(" ".join(shlex.quote(x) for x in argv))
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    os.chmod(path, 0o755)
    return path


def cmd_deliver(a) -> int:
    s = load_summary(a.summary)
    label = submission_label(s)
    warnings = []
    froot = flinders_root()
    mode = a.mode if a.mode != "auto" else ("raw-only" if s.get("raw_data_only") is True else "analysis")
    session_dir = os.path.abspath(a.session) if a.session else None
    if mode == "analysis" and not session_dir:
        raise Stop(EXIT_DECIDE, {"error": "--session is required for an analysis delivery",
                                 "hint": "a raw-data-only submission uses --mode raw-only"})
    if session_dir and not os.path.isdir(session_dir):
        raise Stop(EXIT_DECIDE, {"error": f"session dir not found: {session_dir}"})
    output_dir = sibling("session").paths_for(session_dir)["output_dir"] if session_dir else None
    if mode == "raw-only" and a.include_raw == "no":
        raise Stop(EXIT_DECIDE, {"error": "--mode raw-only delivers the raw links; --include-raw no contradicts it"})
    server_share = server_share_dir(s)
    share = local_share_dir(s, warnings)

    record, where, foreign = find_service_record(s, session_dir, a.service_project)
    if foreign:
        raise Stop(EXIT_DECIDE, {"error": f"this session belongs to {foreign}, not {label}",
                                 "session": session_dir})
    ownership = None
    if session_dir:
        if record:
            owned, why = session_ownership(session_dir, record)
            if not owned:
                raise Stop(EXIT_DECIDE, {"error": f"the session does not belong to {label}: {why}",
                                         "hint": "deliver the session made for this submission"})
            ownership = why
        elif not a.force:
            raise Stop(EXIT_DECIDE, {"error": f"no stage record for {label} ({where}), so nothing proves "
                                              f"this session is {label}'s",
                                     "hint": "run `stage --apply`, or pass --force if staff confirm the session"})
        else:
            ownership = "NOT CHECKED (--force)"
            warnings.append("--force: no stage record, so the session was not checked against the submission")
    want_raw = mode == "raw-only" or a.include_raw == "yes" or (a.include_raw == "auto" and s.get("raw_data_only") is True)
    if want_raw and not record:
        raise Stop(EXIT_DECIDE, {"error": f"raw data was requested but the staged service project cannot be found: {where}"})

    today = dt.date.today().isoformat()
    suffix = safe_component(a.label) if a.label else (f"raw_data_{today}" if mode == "raw-only" else f"analysis_{today}")
    folder = f"{label}_{suffix}"
    delivery = os.path.join(share, folder)
    share_raw = os.path.join(share, "raw")
    unsafe = symlinked_components(delivery, froot) + (symlinked_components(share_raw, froot) if want_raw else [])
    if unsafe:
        raise Stop(EXIT_DECIDE, {"error": "refusing to write: the delivery path passes through a symlink "
                                          "or leaves the Flinders root", "problems": unsafe})
    resume = False
    if os.path.isdir(delivery) and os.listdir(delivery):
        try:
            with open(os.path.join(delivery, DELIVERY_MARKER)) as fh:
                m = json.load(fh)
        except (OSError, ValueError):
            m = {}
        if m.get("session") == session_dir and m.get("folder") == folder and m.get("mode") == mode:
            resume = True
            warnings.append("resuming an interrupted delivery into its own folder")
        else:
            raise Stop(EXIT_DECIDE, {"error": f"{delivery} already has files; mixing two deliveries in one "
                                              f"folder leaves stale results beside new ones",
                                     "hint": "deliver again under a new name: --label <suffix>"})

    if mode == "raw-only":
        items = ([_plan_file(os.path.join(output_dir, n), n) for n in ("methods.md", "methods.docx")]
                 if output_dir else [{"name": "methods.md", "status": "SKIPPED", "required": False,
                                      "reason": "no --session given (run make_methods.py on the raw files "
                                                "and pass its session to include it)"}])
    else:
        items = plan_delivery(output_dir)
    ok = [i for i in items if i["status"] == "OK"]
    total = sum(i["bytes"] for i in ok)
    missing_required = [i["name"] for i in items if i.get("required") and i["status"] != "OK"]

    raw_plan = []
    if want_raw:
        for f in record.get("files") or []:
            name = os.path.basename(f.rstrip("/"))
            target = os.path.realpath(f)
            if not os.path.exists(target):
                raw_plan.append({"name": name, "status": "SKIPPED", "reason": "staged raw file no longer exists"})
            elif not is_under(target, froot):
                raw_plan.append({"name": name, "status": "SKIPPED",
                                 "reason": "outside the Flinders root -- the file server and Bioshare may not see it"})
            else:
                raw_plan.append({"name": name, "status": "OK", "target": target,
                                 "link": os.path.relpath(target, os.path.realpath(share_raw))})
        warnings.append(RAW_WHITELIST_WARNING)
    if not record and not want_raw:
        warnings.append(f"service project: {where}")
    plan = {"applied": bool(a.apply), "mode": mode, "internal_id": s.get("internal_id"), "id": s.get("id"),
            "delivery_dir": delivery, "share_dir": server_share, "share_dir_local": share,
            "session": session_dir, "session_ownership": ownership,
            "n_files": len(ok), "bytes": total, "gb": round(total / 1024 ** 3, 3), "max_gb": a.max_gb,
            "items": [{k: v for k, v in i.items() if k != "src"} for i in items],
            "raw": {"wanted": want_raw, "include_raw": a.include_raw, "links": raw_plan,
                    "service_project": where if record else None},
            "warnings": warnings}
    if missing_required:
        plan["error"] = "required deliverable missing: " + ", ".join(missing_required)
        plan["hint"] = "run step 9 (make_analysis_html.py) before delivering"
        emit(plan)
        return EXIT_DECIDE
    record_dir = session_dir or (record or {}).get("work_dir")
    if total > a.max_gb * 1024 ** 3 and not a.no_size_guard:
        plan["deliver_job"] = write_deliver_job(a, session_dir, record_dir, label, total, warnings)
        plan["error"] = (f"{total / 1024 ** 3:.1f} GB exceeds --max-gb {a.max_gb}: too big for the login "
                         f"node. Submit the job instead: sbatch {plan['deliver_job']}")
        emit(plan)
        return EXIT_GUARD
    if not a.apply:
        plan["next_step"] = "re-run with --apply to copy"
        emit(plan)
        return EXIT_OK

    # ---- apply. Everything below records what happened; nothing may leave a folder with no
    # MANIFEST, and any failure ends as verified:false with a non-zero exit.
    expected, fatal = set(), None
    raw_result = {"created": 0, "existing": 0, "replaced": 0, "skipped": 0}
    results_link = None
    readme_ok = False
    try:
        mkdirs_nofollow(delivery, froot, warnings)
        write_text_nofollow(os.path.join(delivery, DELIVERY_MARKER),
                            json.dumps({"session": session_dir, "folder": folder, "mode": mode,
                                        "started": now_utc()}, indent=2), froot, warnings)
        for i in ok:
            try:
                copy_nofollow(i["src"], os.path.join(delivery, i["name"]), froot, warnings)
                expected.add(i["name"])
            except (OSError, UnsafePath) as e:
                i.update(status="SKIPPED", reason=f"copy failed: {getattr(e, 'strerror', None) or e}")
        if raw_plan:
            mkdirs_nofollow(share_raw, froot, warnings)
        for r in raw_plan:
            if r["status"] != "OK":
                raw_result["skipped"] += 1
                continue
            dst = os.path.join(share_raw, r["name"])
            try:
                if os.path.islink(dst):
                    if os.readlink(dst) == r["link"]:
                        raw_result["existing"] += 1
                        continue
                    os.unlink(dst)
                    os.symlink(r["link"], dst)
                    raw_result["replaced"] += 1
                elif os.path.lexists(dst):
                    r.update(status="SKIPPED", reason="a real file already has this name in raw/")
                    raw_result["skipped"] += 1
                else:
                    os.symlink(r["link"], dst)
                    raw_result["created"] += 1
            except OSError as e:
                r.update(status="SKIPPED", reason=f"link failed: {e.strerror or e}")
                raw_result["skipped"] += 1
        if record:
            sp = record.get("service_dir") or where
            try:
                problems = symlinked_components(sp, froot) if is_under(sp, froot) else [f"{sp} is outside the Flinders root"]
                if problems:
                    raise UnsafePath("; ".join(problems))
                rel = os.path.relpath(os.path.realpath(delivery), os.path.realpath(sp))
                base, n = f"results_{today}", 1
                while True:
                    dst = os.path.join(sp, base if n == 1 else f"{base}_{n}")
                    if os.path.islink(dst) and os.readlink(dst) == rel:
                        results_link = {"link": dst, "state": "existing"}
                        break
                    if not os.path.lexists(dst):
                        os.symlink(rel, dst)
                        results_link = {"link": dst, "state": "created"}
                        break
                    n += 1
            except (OSError, UnsafePath) as e:
                results_link = {"state": "SKIPPED", "reason": f"{getattr(e, 'strerror', None) or e}"}
                warnings.append(f"results_ link in the service project not made: {results_link['reason']}")
        raw_linked = [r["name"] for r in raw_plan if r["status"] == "OK"]
        write_text_nofollow(os.path.join(delivery, "README.md"),
                            build_readme(s, expected, raw_linked, mode), froot, warnings)
        expected.add("README.md")
        readme_ok = True
    except Exception as e:                               # recorded in MANIFEST + delivery.json, exit 2
        fatal = f"{type(e).__name__}: {e}"
        note("delivery stopped by an error:\n" + traceback.format_exc())

    manifest_lines = [f"# MANIFEST -- {label} delivery ({now_utc()})",
                      "# [OK] = included.  [SKIPPED] <name> -- <reason> = not included, and why."]
    for i in items:
        manifest_lines.append(f"[OK] {i['name']}" if i["status"] == "OK" and i["name"] in expected else
                              f"[SKIPPED] {i['name']} -- {i.get('reason') or 'not copied (delivery stopped early)'}")
    for r in raw_plan:
        manifest_lines.append(f"[OK] ../raw/{r['name']} -- linked to the raw instrument file"
                              if r["status"] == "OK" else f"[SKIPPED] ../raw/{r['name']} -- {r['reason']}")
    if not want_raw:
        manifest_lines.append(f"[SKIPPED] ../raw/ -- not requested (--include-raw {a.include_raw}"
                              + (", CoreOmics says the Core does the analysis)" if a.include_raw == "auto" else ")"))
    if results_link and results_link.get("state") == "SKIPPED":
        manifest_lines.append(f"[SKIPPED] results link in the Core service project (staff bookkeeping) -- "
                              f"{results_link['reason']}")
    manifest_lines.append("[OK] README.md" if readme_ok else "[SKIPPED] README.md -- not written (delivery stopped early)")
    if fatal:
        manifest_lines.append(f"[SKIPPED] remainder of this delivery -- stopped by an error: {fatal}")
    manifest_ok, errors = False, []

    def write_manifest(checksum_line):
        write_text_nofollow(os.path.join(delivery, "MANIFEST.txt"),
                            "\n".join(manifest_lines + [checksum_line]) + "\n", froot, warnings)

    try:
        write_manifest("[OK] checksums.sha256")
        manifest_ok = True
        expected.add("MANIFEST.txt")
        sums = []
        for dirpath, dirnames, filenames in os.walk(delivery):
            dirnames.sort()
            for fn in sorted(filenames):
                p = os.path.join(dirpath, fn)
                rel = os.path.relpath(p, delivery)
                if rel not in ("checksums.sha256", DELIVERY_MARKER) and os.path.isfile(p) and not os.path.islink(p):
                    sums.append(f"{sha256_of(p)}  {rel}")
        write_text_nofollow(os.path.join(delivery, "checksums.sha256"), "\n".join(sums) + "\n", froot, warnings)
        expected.add("checksums.sha256")
    except (OSError, UnsafePath) as e:
        errors.append(f"MANIFEST/checksums: {getattr(e, 'strerror', None) or e}")
        if manifest_ok:
            try:
                write_manifest(f"[SKIPPED] checksums.sha256 -- {getattr(e, 'strerror', None) or e}")
            except (OSError, UnsafePath) as e2:
                errors.append(f"MANIFEST rewrite: {getattr(e2, 'strerror', None) or e2}")

    failed_required = [i["name"] for i in items if i.get("required") and i["name"] not in expected]
    violations, share_notes = verify_delivery(
        delivery, share, froot, expected, own_raw={r["name"] for r in raw_plan if r["status"] == "OK"})
    warnings.extend(share_notes)
    verified = bool(manifest_ok and not violations and not failed_required and not fatal and not errors)
    if verified:
        try:
            os.unlink(os.path.join(delivery, DELIVERY_MARKER))
        except OSError as e:
            warnings.append(f"could not remove {DELIVERY_MARKER}: {e.strerror}")
    files_now = [os.path.join(dp, f) for dp, _, fs in os.walk(delivery) for f in fs] if os.path.isdir(delivery) else []
    result = {"applied": True, "mode": mode, "internal_id": s.get("internal_id"), "id": s.get("id"),
              "label": label, "delivered_at": now_utc(), "delivery_dir": delivery,
              "share_dir": server_share, "share_dir_local": share, "session": session_dir,
              "session_ownership": ownership, "resumed": resume,
              "n_files": len(files_now),
              "bytes": sum(os.lstat(p).st_size for p in files_now if not os.path.islink(p)),
              "contents": sorted(x for x in os.listdir(delivery) if x != DELIVERY_MARKER) if os.path.isdir(delivery) else [],
              "skipped": [f"{i['name']} -- {i.get('reason', '')}" for i in items if i["name"] not in expected]
                         + [f"../raw/{r['name']} -- {r['reason']}" for r in raw_plan if r["status"] != "OK"]
                         + ([f"results link -- {results_link['reason']}"] if results_link and results_link.get("state") == "SKIPPED" else []),
              "raw_links": raw_result, "raw_whitelist_unverified": bool(raw_plan), "results_link": results_link,
              "verified": verified, "violations": violations, "failed_required": failed_required,
              "error": fatal, "errors": errors, "warnings": warnings}
    djson = os.path.join(record_dir, "delivery.json") if record_dir else None
    try:
        if not djson:
            raise OSError("no session or work dir to write delivery.json into")
        with open(djson, "w") as fh:
            json.dump(result, fh, indent=2)
        result["delivery_json"] = djson
    except OSError as e:
        result["delivery_json"] = None
        result["errors"].append(f"delivery.json not written: {e}")
        result["verified"] = verified = False
    emit(result)
    if not verified:
        note("delivery is NOT verified -- do not share it")
        return EXIT_DECIDE
    return EXIT_OK


# ----------------------------------------------------------------------- bioshare --
SHARE_NAME_OK = re.compile(r"^[\w\d\s'\"\.!\?\-:,]+$")


def share_name(summary: dict) -> str:
    """`<PI last>, <PI first>: <internal_id>` (the plugin model's own default), reduced to
    what the serializer accepts: ^[\\w\\d\\s'".!?\\-:,]+$."""
    pi = summary.get("pi") or {}
    last, first, label = _s(pi.get("last")), _s(pi.get("first")), submission_label(summary)
    if last and first:
        raw = f"{last}, {first}: {label}"
    elif last:
        raw = f"{last}: {label}"
    else:
        raw = f"UC Davis Proteomics Core: {label}"
    clean = re.sub(r"\s+", " ", re.sub(r"[^\w\s'\"\.!\?\-:,]", " ", raw)).strip()
    return clean if SHARE_NAME_OK.match(clean) else f"UC Davis Proteomics Core: {safe_component(label, 'submission')}"


def linked_share(shares: list, share_dir: str):
    want = share_dir.rstrip("/")
    for sh in shares:
        if _s(sh.get("link_to_path")).rstrip("/") == want:
            return sh
    return None


def submission_contacts(summary: dict, summary_path: str):
    if isinstance(summary.get("contacts"), list):
        return summary["contacts"], "submission_summary.json"
    sj = os.path.join(os.path.dirname(os.path.abspath(summary_path)), "submission.json")
    try:
        with open(sj) as fh:
            return contact_rows(json.load(fh)), sj
    except (OSError, ValueError):
        return None, f"not known: no contacts in the summary and no readable {sj}"


def delivery_gate(summary: dict, delivery_path, share: str) -> dict:
    """send --apply shares only a delivery that verified, for THIS submission, in THIS share."""
    if not delivery_path:
        raise Stop(EXIT_DECIDE, {"error": "send --apply needs --delivery <delivery.json>, pulled from HIVE, so a "
                                          "delivery that failed verification is never shared",
                                 "hint": "bash scripts/hive_exec.sh --get '<session>/delivery.json' ./delivery.json"})
    d = load_json_file(delivery_path, "delivery.json")
    problems = []
    if d.get("verified") is not True:
        problems.append("the delivery did not verify (verified is not true)")
    if _s(d.get("internal_id")) != _s(summary.get("internal_id")):
        problems.append(f"delivery.json is for {d.get('internal_id')!r}, not {summary.get('internal_id')!r}")
    if _s(d.get("share_dir")).rstrip("/") != share:
        problems.append(f"delivery.json share_dir {d.get('share_dir')!r} is not this submission's share {share!r}")
    if problems:
        raise Stop(EXIT_DECIDE, {"error": "refusing to share: " + "; ".join(problems), "delivery": delivery_path})
    return d


def cmd_bioshare(a) -> int:
    s = load_summary(a.summary)
    sid = _s(s.get("id"))
    share = server_share_dir(s)                         # never re-derived from a local root
    shares = list_shares(sid)
    linked = linked_share(shares, share)
    base = f"plugins/bioshare/submissions/{sid}/submission_shares/"
    listing = [dict(sh, linked=sh is linked) for sh in shares]

    if a.action == "status":
        emit({"internal_id": s.get("internal_id"), "share_dir": share, "shares": listing,
              "linked": linked, "url": (linked or {}).get("url")})
        return EXIT_OK

    if a.action == "ensure":
        if linked:
            emit({"already_exists": True, "share": linked, "url": linked.get("url"), "applied": False})
            return EXIT_OK
        payload = {"submission": sid, "name": share_name(s),
                   "notes": f"UC Davis Proteomics Core results for {submission_label(s)}",
                   "link_to_path": share}
        if not a.apply:
            emit({"applied": False, "would_post": {"method": "POST", "url": api_url(base), "json": payload},
                  "next_step": "re-run with --apply to create the share (the share dir must already exist "
                               "on HIVE -- `deliver --apply` makes it)"})
            return EXIT_OK
        created = api_call("POST", base, payload=payload)
        if not isinstance(created, dict) or created.get("id") in (None, ""):
            raise ApiError(None, f"creating the share returned no share record: {json.dumps(created)[:300]}",
                           api_url(base))
        emit({"applied": True, "share": created, "url": created.get("url")})
        return EXIT_OK

    # send
    if not linked:
        raise Stop(EXIT_DECIDE, {"error": f"no Bioshare share is linked to {share}",
                                 "hint": "run `bioshare ensure --apply` first"})
    pi, sub = s.get("pi") or {}, s.get("submitter") or {}
    contacts, contacts_source = submission_contacts(s, a.summary)
    recipients = {"submitter": {"name": sub.get("name"), "email": sub.get("email")},
                  "pi": {"name": pi.get("name"), "email": pi.get("email")},
                  "contacts": contacts, "contacts_source": contacts_source}
    url = f"{base}{linked.get('id')}/share/"
    payload = {"email": bool(a.email)}
    if not a.apply:
        out = {"applied": False, "recipients": recipients, "share": linked,
               "would_post": {"method": "POST", "url": api_url(url), "json": payload},
               "notification_email": payload["email"],
               "next_step": "OUTWARD-FACING: only with the staff member's explicit yes, re-run with --apply "
                            "--delivery <delivery.json>"}
        if a.delivery:
            try:
                d = delivery_gate(s, a.delivery, share)
                out["delivery_gate"] = "PASS"
                if d.get("raw_whitelist_unverified"):
                    out["check_first"] = RAW_WHITELIST_WARNING
            except Stop as e:
                out["delivery_gate"] = e.payload
        emit(out)
        return EXIT_OK
    delivery_gate(s, a.delivery, share)
    resp = api_call("POST", url, payload=payload)
    if not isinstance(resp, dict) or "results" in resp:
        raise ApiError(None, f"unexpected response to the share request: {json.dumps(resp)[:300]}", api_url(url))
    emit({"applied": True, "recipients": recipients, "notification_email": payload["email"],
          "response": resp, "url": linked.get("url")})
    return EXIT_OK


# -------------------------------------------------------------------- email-draft --
def cmd_email_draft(a) -> int:
    s = load_summary(a.summary)
    delivery = load_json_file(a.delivery, "delivery.json") if a.delivery else {}
    url = _s(a.share_url)
    if not url:
        for sh in s.get("existing_shares") or []:
            if s.get("share_dir") and _s(sh.get("link_to_path")).rstrip("/") == _s(s["share_dir"]).rstrip("/"):
                url = _s(sh.get("url"))
    label = submission_label(s)
    sub, pi = s.get("submitter") or {}, s.get("pi") or {}
    placeholders = []
    to = _s(sub.get("email"))
    if not to:
        placeholders.append("submitter email")
    cc = _s(pi.get("email")) if _s(pi.get("email")).lower() != to.lower() else ""
    link = url or "[BIOSHARE LINK -- run `core_submission.py bioshare ensure` and paste the share URL]"
    if not url:
        placeholders.append("Bioshare link")
    contents = set(delivery.get("contents") or [])
    raw = delivery.get("raw_links") or {}
    n_raw = (raw.get("created") or 0) + (raw.get("existing") or 0) + (raw.get("replaced") or 0)
    raw_only = delivery.get("mode") == "raw-only"
    L = ["<!-- DRAFT -- not sent. Review it, then send it from your own mail client. -->", "",
         f"To: {to or '[SUBMITTER EMAIL]'}"]
    if cc:
        L.append(f"Cc: {cc}")
    samples = f" ({s['n_samples']} samples)" if s.get("n_samples") else ""
    if raw_only:
        L += [f"Subject: Raw data for {label} is ready", "", f"Hi {sub.get('first') or 'there'},", "",
              f"The raw instrument files for your submission {label}{samples} are ready on Bioshare:", "",
              link, "", "They are in the raw folder of the share. As requested, no search or analysis was run."]
    else:
        L += [f"Subject: Proteomics results for {label} are ready", "", f"Hi {sub.get('first') or 'there'},", "",
              f"The UC Davis Proteomics Core has finished the analysis for your submission {label}{samples}. "
              f"Your results are on Bioshare:", "", link, ""]
        if not contents or "Analysis_Report.html" in contents:
            L.append("Please start with Analysis_Report.html. It is the whole report in one file and opens in "
                     "any web browser -- the quality-control panels come first, then the results.")
        if contents:
            L += ["", "The folder also contains:"]
            for name, text in (("tables", "the result tables"),
                               ("figures", "every figure as an image file"),
                               ("Analysis_Report.pdf", "the report as a PDF, for printing"),
                               ("methods.md", "a Methods section you can adapt for a manuscript"),
                               ("reproducibility", "a full record of the software and settings used"),
                               ("search", "the search engine output")):
                if name in contents:
                    L.append(f"  - {text}")
        if n_raw:
            L += ["", "Your raw instrument files are in the raw folder of the same share."]
    if delivery.get("n_files") and delivery.get("bytes"):
        L.append(f"  ({delivery['n_files']} files, {human_size(delivery['bytes'])} in the delivery folder)")
    ack = ("the wording, including the instrument grant, is at the end of methods.md."
           if "methods.md" in contents or (not contents and not raw_only) else
           "reply to this email and we will send the wording, including the instrument grant.")
    L += ["", "If you use these data in a publication, please acknowledge the UC Davis Proteomics Core; " + ack, "",
          "Let us know if anything is unclear.", "", "Best regards,", "UC Davis Proteomics Core", ""]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\n".join(L))
    emit({"wrote": os.path.abspath(a.out), "to": to or None, "cc": cc or None, "share_url": url or None,
          "placeholders": placeholders, "sent": False})
    return EXIT_DECIDE if placeholders else EXIT_OK


# --------------------------------------------------------------------------- main --
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="LOCAL: fetch a submission + neighbours + Bioshare shares")
    f.add_argument("submission", help="807, 0807, PROT_0807, prot-807, #807, or the 12-hex CoreOmics id")
    f.add_argument("--out", required=True, help="directory for submission.json, samples.tsv, submission_summary.json")
    f.add_argument("--neighbor-days", type=int, default=NEIGHBOR_WINDOW_DAYS,
                   help=f"record submissions within +/- this many days (default {NEIGHBOR_WINDOW_DAYS})")
    f.set_defaults(func=cmd_fetch)

    idn = sub.add_parser("identify", help="LOCAL: which submission is this data? (names first, "
                                          "then sample ids; never guesses)")
    idn.add_argument("paths", nargs="*", help="raw files and/or their folder, as the user gave them")
    idn.add_argument("--files-from", default=None, help="a file listing raw paths, one per line")
    idn.add_argument("--text", default="", help="the user's message, for 'PROT_0756' / 'submission 756'")
    idn.add_argument("--max-days", type=int, default=240,
                     help="a run counts for a submission made up to this many days before it (default 240)")
    idn.add_argument("--no-lookup", action="store_true",
                     help="names only: do not ask CoreOmics which submission uses these sample ids")
    idn.set_defaults(func=cmd_identify)

    lo = sub.add_parser("locate", help="ON HIVE: find each sample's raw file (writes a proposal)")
    lo.add_argument("--summary", required=True)
    lo.add_argument("--out", required=True, help="directory for files.txt, sample_files.tsv, locate.json")
    lo.add_argument("--raw-root", default=None, help="default <CORE_FLINDERS_ROOT>/Data/raw_data")
    lo.add_argument("--max-days", type=int, default=240,
                    help="accept files acquired up to this many days after submission (default 240; "
                         "never more than the neighbour window fetch recorded)")
    lo.add_argument("--allow-partial", action="store_true",
                    help="staff confirmed the unmatched samples were not run: downgrade that gate to WARN")
    lo.add_argument("--accept-ambiguous", action="store_true",
                    help="staff confirmed runs whose label another submission also uses are THIS "
                         "submission's (recorded in locate.json)")
    lo.add_argument("--files-from", default=None,
                    help="staff-curated file list: skip matching, still validate and write the TSV")
    lo.set_defaults(func=cmd_locate)

    st = sub.add_parser("stage", help="ON HIVE: link the raw files into the service directory (dry run)")
    st.add_argument("--summary", required=True)
    st.add_argument("--files", required=True, help="files.txt from locate")
    st.add_argument("--sample-files", default=None, help="sample_files.tsv (default: next to --files)")
    st.add_argument("--service-dir", default=None,
                    help="the group folder to use when resolution is ambiguous (under the service root)")
    st.add_argument("--apply", action="store_true", help="create the links and folders")
    st.set_defaults(func=cmd_stage)

    co = sub.add_parser("conditions", help="build conditions.csv from CoreOmics conditions")
    co.add_argument("--summary", required=True)
    co.add_argument("--sample-files", required=True, help="sample_files.tsv from locate")
    co.add_argument("--out", required=True, help="conditions.csv to write")
    co.set_defaults(func=cmd_conditions)

    de = sub.add_parser("deliver", help="ON HIVE: copy deliverables into the Bioshare share dir (dry run)")
    de.add_argument("--summary", required=True)
    de.add_argument("--session", default=None, help="the analysis session directory (optional for raw-only)")
    de.add_argument("--mode", choices=["auto", "analysis", "raw-only"], default="auto",
                    help="auto = raw-only when CoreOmics says 'I only require raw data', else analysis")
    de.add_argument("--label", default=None, help="folder suffix instead of analysis_<date> / raw_data_<date>")
    de.add_argument("--include-raw", choices=["auto", "yes", "no"], default="auto",
                    help="auto = only when CoreOmics says the submitter wants raw data")
    de.add_argument("--service-project", default=None, help="staged service project (default: discovered)")
    de.add_argument("--force", action="store_true",
                    help="deliver a session with NO stage record (staff confirmed it is this submission's)")
    de.add_argument("--apply", action="store_true", help="copy, link, verify")
    de.add_argument("--max-gb", type=float, default=5.0,
                    help="above this, write deliver_job.sh and exit 5 instead of copying (default 5)")
    de.add_argument("--no-size-guard", action="store_true", help="used by deliver_job.sh itself")
    de.set_defaults(func=cmd_deliver)

    bs = sub.add_parser("bioshare", help="LOCAL: status | ensure | send (ensure/send are dry runs)")
    bs.add_argument("action", choices=["status", "ensure", "send"])
    bs.add_argument("--summary", required=True)
    bs.add_argument("--delivery", default=None,
                    help="send: delivery.json from deliver (required with --apply; must be verified)")
    bs.add_argument("--apply", action="store_true", help="perform the POST")
    g = bs.add_mutually_exclusive_group()
    g.add_argument("--email", dest="email", action="store_true",
                   help="send: have Bioshare email the recipients")
    g.add_argument("--no-email", dest="email", action="store_false", help="send: grant access silently (default)")
    bs.set_defaults(func=cmd_bioshare, email=False)

    em = sub.add_parser("email-draft", help="LOCAL: draft the results email (never sends)")
    em.add_argument("--summary", required=True)
    em.add_argument("--delivery", default=None, help="delivery.json from deliver")
    em.add_argument("--share-url", default=None)
    em.add_argument("--out", required=True)
    em.set_defaults(func=cmd_email_draft)

    a = ap.parse_args(argv)
    try:
        return a.func(a)
    except Stop as e:
        emit(e.payload)
        note(e.payload.get("error", "stopped"))
        return e.code
    except ApiError as e:
        emit({"error": "CoreOmics request failed", "status": e.status, "detail": e.detail, "url": e.url})
        note(f"CoreOmics: {e}")
        return EXIT_UNREACHABLE


if __name__ == "__main__":
    sys.exit(main())
