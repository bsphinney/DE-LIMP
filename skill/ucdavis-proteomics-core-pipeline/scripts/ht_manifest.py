"""
ht_manifest.py  --  Turn a Core high-throughput submission number into a search file list.

A high-throughput submission is not a folder. Submission 0793 can occupy two trays, and
the second tray's filenames may never mention 0793 at all -- the extent of a submission is
inferred from the acquisition counter, not stored as a directory. It also contains runs
that must NOT be searched as customer samples (well blanks, HeLa standards), and runs the
operator already knows came out badly. Globbing a directory gets all of this wrong,
silently, and you find out six hours later.

STAN knows which files a submission is. This skill knows how to search them. Neither should
learn the other's job, so this script does exactly one thing: it asks STAN, checks the
answer is safe to search, and writes a file list the rest of the flow consumes normally.
After that, nothing about the run is HT-specific -- step 2 onward is the ordinary flow.

**Runs ON HIVE**, like fran_deposit.py -- the orchestrator invokes it through
`hive_exec.sh`, because that is where STAN, its database credential and the raw files are.

    bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py fetch 0793 --out ~/ht0793'
    bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py link 0793'

The share token and the Entra cookie are read from a FILE (--share-token-file, --cookie-file).
On the command line they are in the process list and in the session's commands.log, which the
run registry copies into a folder the whole Core group can read (release review, 2.8.0). Every
line this script prints goes through _say(), and ht_manifest.json through _scrubbed(): the
token, cookie and PG credential are masked, whatever the server echoes back -- in an error or
in a 200 reply's plates, counts or example paths (release verification, 2.8.0).

`fetch` writes `files.txt` (one absolute .d path per line, for --raw/--files) and
`ht_manifest.json` (the full STAN payload plus the gate results).

STAN's list is CHECKED, not trusted (2026-10-05, 2.10.0): STAN assigns runs to a submission by
tray and injection window, whoever's they are, and a 96-sample plate came back as 136 files --
95 of its runs, one real run missing, and 41 earlier runs of other projects on the same tray, at
other methods, 5 of them DDA -- with every gate PASS. So files.txt holds only runs NAMED for the
submission (`named_for_submission`), the runs must declare one acquisition method
(`one_method`), and with the CoreOmics record (`--summary`) every sample must have a run and
every run a sample (`sample_list`). Each FAILS closed, naming every file or sample it is about.

    bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py fetch <n> \\
        --summary ~/core/PROT_<n>/submission_summary.json --out ~/ht<n>'

`link` only READS. There is deliberately no write-back: STAN v1.0.42 added an
`ht_searches` table and an `ht-record-search` command, and v1.0.43 reverted both, because
a table in STAN is a second copy of what FRAN already knows and the copy is the one that
goes stale. The submission number is in the raw filenames, so FRAN's
`/api/internal/submission/{id}` already answers "what was searched for 0793". `link` just
confirms that loop closed.

EXIT CODES -- the orchestrator must honour these:
    0  safe to proceed
    2  HARD GATE FAILED. Do not search. The run would silently cover a subset.
    3  STAN could not be reached, or is too old to have `ht-manifest`.
    4  (`link` only) not visible in FRAN yet. NOT a search failure.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_report_runs import (distinct_inputs, repeated_names, repeats_note,  # noqa: E402
                               repeats_text)   # repeated inputs: the one rule

# Which runs are this submission's (core_submission.ht_pattern, 2.9.1: the number right after the
# leading date, bare or PROT_) and which run is which sample (curated_rows, the rule `locate
# --files-from` applies) have ONE definition each, as does the method token in a run name
# (fran_deposit.RUN_METHOD_RE, mirrored in FRAN). A partial copy of scripts/ without them cannot
# check STAN's list, and those gates then FAIL, saying why -- they never pass for want of a rule.
try:
    import core_submission as _cs
    CORE_UNAVAILABLE = None
except ImportError as _e:
    _cs = None
    CORE_UNAVAILABLE = f"core_submission.py could not be imported ({type(_e).__name__}: {_e})"
try:
    from fran_deposit import RUN_METHOD_RE
    METHOD_UNAVAILABLE = None
except ImportError as _e:
    RUN_METHOD_RE = None
    METHOD_UNAVAILABLE = f"fran_deposit.py could not be imported ({type(_e).__name__}: {_e})"

# DIA or DDA, when a run name writes it as its own word right after the method token
# (`60SPD_DIA-`, `60-spd-dda-`). A plate run's sample name follows the speed directly
# (`100spd_Diabetic1`), so a word that merely begins with the letters does not count.
MODE_AFTER_METHOD = re.compile(r"[-_]+(dia|dda)(?![a-z0-9])", re.I)

# STAN's own venv on HIVE. Overridable because nothing should hardcode one person's path.
DEFAULT_STAN = "/quobyte/proteomics-grp/brett/stan_venv/bin/stan"

# Where to look for the PG Farm credential, in order. OWNER_TOKEN is NOT a short-lived token,
# whatever its name says: it holds the long-lived 512-character secret of the
# `genome-proteomics-service-account` (STAN CLAUDE.md, "PG Farm auth"; 512 bytes, mode 0600 on
# HIVE, 2026-10-01), which STAN exchanges for a fresh JWT on every use. Whoever can read it can
# act as the service account -- write access to STAN's tables, no per-person trail -- for as
# long as the secret lives. It must NEVER be copied or made group-readable. Skill 2.9 and
# earlier called it a "7-day token" and suggested a group-readable copy at
# /quobyte/proteomics-grp/etc/pgfarm_token, searched first; that was wrong, and that path is no
# longer searched. Everyone but the owner uses the HTTP path with a share token (--http,
# --share-token-file), which needs no database credential at all.
OWNER_TOKEN = "/quobyte/proteomics-grp/brett/.pgfarm_token"
TOKEN_CANDIDATES = ("~/.pgfarm_token", OWNER_TOKEN)
DEFAULT_TOKEN = OWNER_TOKEN  # kept for the error message's "tried" list

# A tray is 96 wells. A submission reporting a handful of samples has probably been
# mistyped -- warn, don't block, because a genuine small submission is legal.
MIN_PLAUSIBLE_SAMPLES = 12
MAX_PLAUSIBLE_PLATES = 2


def _stan_bin(explicit: str | None) -> str:
    for cand in (explicit, os.environ.get("STAN_BIN"), DEFAULT_STAN, "stan"):
        if not cand:
            continue
        if cand == "stan" or os.path.exists(cand):
            return cand
    return DEFAULT_STAN


def _env(token_path: str | None) -> dict:
    """STAN's Postgres backend needs a credential. Fail loudly and specifically if the
    caller cannot read it -- the default token is mode 0600 and owned by one person, so
    'permission denied' is the EXPECTED failure for every other Core member, and it must
    not be mistaken for 'the submission does not exist'."""
    env = dict(os.environ)
    env["STAN_DB_BACKEND"] = env.get("STAN_DB_BACKEND", "pg")
    if env.get("PGPASSWORD"):
        _remember(env["PGPASSWORD"])
        return env

    explicit = token_path or env.get("STAN_PG_TOKEN")
    candidates = [explicit] if explicit else list(TOKEN_CANDIDATES)

    tried = []
    for cand in candidates:
        p = os.path.expanduser(cand)
        try:
            tok = open(p).read().strip()
        except PermissionError:
            tried.append(f"{p} — permission denied (owned by another user, mode 0600)")
            continue
        except (FileNotFoundError, NotADirectoryError):
            tried.append(f"{p} — not found")
            continue
        if tok:
            env["PGPASSWORD"] = tok
            _remember(tok)
            return env
        tried.append(f"{p} — empty")

    detail = "\n".join(f"    {t}" for t in tried)
    sys.exit(
        f"[ht_manifest] no usable PG Farm credential. Tried:\n{detail}\n"
        f"\n"
        f"  The CLI path reads STAN's database credential, and only its owner can: the file\n"
        f"  holds the long-lived service-account SECRET (not a short-lived token), mode 0600 on\n"
        f"  purpose. NEVER copy it or make it group-readable -- whoever can read it can act as\n"
        f"  the service account, which can write to STAN's tables.\n"
        f"\n"
        f"  Use the hosted dashboard instead (works for every Core member, headless):\n"
        f"    --http https://ucd.stan-proteomics.org --share-token-file <file>\n"
        f"  (a file holding the submission's share token from its HT tab, mode 600 -- never the\n"
        f"  token on the command line). Your OWN credential, if you have one, can be given as\n"
        f"  STAN_PG_TOKEN=<a file you own> or PGPASSWORD.")


def _read_secret_file(path: str, what: str) -> str:
    """A token or cookie from a file: stripped, never echoed. A file other users can read is
    used, with a warning -- the point of the file is that nobody else sees the value."""
    p = os.path.expanduser(path)
    try:
        with open(p) as fh:
            val = fh.read().strip()
        mode = os.stat(p).st_mode
    except OSError as e:
        sys.exit(f"[ht_manifest] cannot read the {what} file {p}: {e.strerror}")
    if not val:
        sys.exit(f"[ht_manifest] the {what} file {p} is empty")
    if mode & 0o077:
        _say(f"[ht_manifest] WARNING: {p} can be read by other users (mode {mode & 0o777:o}); "
              f"chmod 600 it.", file=sys.stderr)
    return val


def _credentials(a) -> tuple:
    """(share token, cookie), files first. The argv forms still work, with a warning."""
    tok = cookie = None
    if getattr(a, "share_token_file", None):
        tok = _read_secret_file(a.share_token_file, "share token")
    elif getattr(a, "share_token", None):
        tok = a.share_token
        _say("[ht_manifest] WARNING: --share-token puts the token in the process list and in "
              "commands.log; use --share-token-file.", file=sys.stderr)
    else:
        tok = os.environ.get("STAN_HT_SHARE_TOKEN")
    if getattr(a, "cookie_file", None):
        cookie = _read_secret_file(a.cookie_file, "cookie")
    elif getattr(a, "cookie", None):
        cookie = a.cookie
        _say("[ht_manifest] WARNING: --cookie puts the session cookie in the process list and in "
              "commands.log; use --cookie-file.", file=sys.stderr)
    _remember(tok, cookie)
    return tok, cookie


def _masked(text, *secrets) -> str:
    """`text` with each secret value replaced, then notify_slack's patterns applied (the skill's
    one list; it catches `token=...` even when the value is not one this run knows)."""
    out = str(text)
    for v in secrets:
        if v:
            out = out.replace(v, "[redacted]")
    try:
        from notify_slack import redact
    except ImportError:          # copied here alone: the exact-value mask above still applied
        return out
    return redact(out)


# Every secret value this run has read -- share token, cookie, PG credential, and each one
# URL-encoded -- so _say() and _scrubbed() can mask them wherever they turn up.
_SECRETS: list = []


def _remember(*values):
    for v in values:
        if v and v not in _SECRETS:
            _SECRETS.extend([v, urllib.parse.quote_plus(v)])


def _say(text="", file=None):
    """print(), masked: the only way this script writes a line."""
    print(_masked(text, *_SECRETS), file=file or sys.stdout)


def _scrubbed(obj):
    """`obj` with every string in it -- keys too -- masked like a printed line."""
    if isinstance(obj, str):
        return _masked(obj, *_SECRETS)
    if isinstance(obj, dict):
        return {_scrubbed(k): _scrubbed(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_scrubbed(v) for v in obj]
    return obj


def _run_stan(argv: list, env: dict) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(argv, capture_output=True, text=True, env=env, timeout=180)
    except FileNotFoundError:
        sys.exit(f"[ht_manifest] STAN not found at {argv[0]}. Pass --stan or set STAN_BIN.")
    except subprocess.TimeoutExpired:
        sys.exit(f"[ht_manifest] STAN timed out after 180s running: {' '.join(argv[:3])}")


def _manifest_over_http(a) -> dict | int:
    """Ask the hosted dashboard instead of the local database.

    This is the path that makes the workflow usable by the whole Core rather than by
    whoever owns the Postgres token. `/api/ht/manifest` returns the same payload under the
    same access rules, and it accepts a per-submission **share token** — which is the only
    one of its three auth paths a headless HIVE job can actually use. Entra sign-in needs a
    browser; a share link does not.
    """
    import urllib.error
    import urllib.request

    base = (a.http or os.environ.get("STAN_HT_URL") or "").rstrip("/")
    tok, cookie = _credentials(a)
    qs = {"q": a.submission, "include": a.include}
    if tok:
        qs["token"] = tok
    url = f"{base}/api/ht/manifest?" + urllib.parse.urlencode(qs)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    if cookie:
        req.add_header("Cookie", cookie)

    def err(msg):
        # every line, not only the one that used to print `url`: a server may echo the request
        _say(msg, file=sys.stderr)

    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        body = ""
        try:
            body = e.read().decode()[:400]
        except Exception:
            pass
        if e.code in (401, 403):
            err(f"[ht_manifest] {base} refused the request ({e.code}).\n"
                f"  HT data is not public — it carries customer sample names and paths.\n"
                f"  Use ONE of:\n"
                f"    --share-token-file <f>  a file holding the per-submission share token from\n"
                f"                            the HT tab (the only option that works headless)\n"
                f"    --cookie-file <f>       a file holding an Entra session cookie from a\n"
                f"                            signed-in browser\n"
                f"  Sign-in is Microsoft Entra ({base}/.auth/login/aad), not CAS.\n"
                f"  server said: {body}")
        else:
            err(f"[ht_manifest] {url} returned HTTP {e.code}: {body}")
        return 3
    except urllib.error.URLError as e:
        err(f"[ht_manifest] cannot reach {base}: {e.reason}")
        return 3
    except json.JSONDecodeError:
        err(f"[ht_manifest] {base} returned non-JSON")
        return 3


def _manifest_over_cli(a) -> dict | int:
    env = _env(a.token)
    stan = _stan_bin(a.stan)
    p = _run_stan([stan, "ht-manifest", a.submission, "--include", a.include], env)
    if p.returncode != 0:
        err = (p.stderr or p.stdout or "").strip()
        if "No such command" in err:
            _say(f"[ht_manifest] the STAN on this machine has no 'ht-manifest' command.\n"
                  f"  {stan}\n"
                  f"  HT support landed in STAN v1.0.43; this venv predates it. Update the\n"
                  f"  venv (pip install -U from the STAN repo) and re-run. Nothing was searched.",
                  file=sys.stderr)
            return 3
        _say(f"[ht_manifest] stan ht-manifest failed (rc={p.returncode}):\n{err}", file=sys.stderr)
        return 3
    try:
        return json.loads(p.stdout)
    except json.JSONDecodeError:
        _say(f"[ht_manifest] STAN returned non-JSON:\n{p.stdout[:500]}", file=sys.stderr)
        return 3


def _base(path: str) -> str:
    return os.path.basename(str(path).rstrip("/"))


def declared_method(name: str):
    """What a run name says about how it was acquired -- '100spd', '60spd DIA', '90m' -- or None
    when it names no method (a re-injection: `<date>_<n>_rerun_<sample>_...`)."""
    m = RUN_METHOD_RE.search(name)
    if not m:
        return None
    mode = MODE_AFTER_METHOD.match(name, m.end())
    return re.sub(r"[-_]", "", m.group(0)).lower() + (f" {mode.group(1).upper()}" if mode else "")


def submission_plate(submission):
    """The test a run name passes when it is named for `submission` (core_submission.ht_pattern),
    or None when there is none: the rule could not be imported, or `submission` is not a
    submission number (normalize_submission: 0793, 793 and PROT_0793 are; a CoreOmics hex id is
    not -- plate runs never carry it)."""
    if _cs is None:
        return None
    try:
        kind, key = _cs.normalize_submission(submission)
    except ValueError:
        return None
    return _cs.ht_pattern(key) if kind == "internal_id" else None


def named_gate(a, files: list, own: list, foreign: list, entries: dict, plate) -> dict:
    """`named_for_submission`: the files STAN listed whose names do not carry the submission
    number, every one of them, with STAN's class, plate, well and injection for each."""
    if plate is None:
        why = CORE_UNAVAILABLE or f"{a.submission!r} is not a submission number"
        return {"gate": "named_for_submission", "status": "FAIL", "n": len(files),
                "detail": f"{why}, so no file name could be checked against the submission; "
                          f"files.txt is left empty"
                          + (". Put the whole skill on this machine again (hive_exec.sh "
                             "--put-skill) and re-run" if CORE_UNAVAILABLE else "")}
    if not foreign:
        return {"gate": "named_for_submission", "status": "PASS", "n": 0}
    # the instrument's acquisition counter orders one instrument's runs: a run before the first
    # one named for the submission cannot be a later tray of it
    first = min((c for c in map(_cs.acquisition_counter, own) if c is not None), default=None)
    before = [f for f in foreign if first is not None
              and (_cs.acquisition_counter(f) or first) < first]
    listed = []
    for f in foreign:
        e = entries.get(os.path.realpath(f)) or {}
        listed.append({"file": f, **{k: e.get(k) for k in ("class", "plate", "well", "injection")}})
    return {"gate": "named_for_submission", "status": "FAIL", "n": len(foreign),
            "detail": f"{len(foreign)} of the {len(files)} file(s) STAN listed are not named for "
                      f"submission {a.submission} (its runs are named <date>_[PROT_]<n>_...). STAN "
                      f"assigns runs to a submission by tray -- every run on a tray the submission "
                      f"used, inside its injection window, whoever's it is -- so these may be "
                      f"another project's runs"
                      + (f"; {len(before)} of them ran before the first run named for it, so they "
                         f"are no continuation of this plate" if before else "")
                      + ". They are NOT in files.txt. Ask the operator whose they are; runs of "
                        "this submission named without its number (a continuation tray) are "
                        "searched only from a curated list (core_submission.py locate --files-from)",
            "files": listed,
            "listed": [f"{x['file']}  (STAN class: {x['class'] or 'not given'})" for x in listed]}


def method_gate(files: list, own: list) -> dict:
    """`one_method`: the acquisition methods the file names declare, grouped. More than one is a
    stop; a name declaring none (a re-injection) is flagged."""
    if RUN_METHOD_RE is None:
        return {"gate": "one_method", "status": "FAIL",
                "detail": f"{METHOD_UNAVAILABLE}, so the files' acquisition methods could not be "
                          f"compared. Put the whole skill on this machine again "
                          f"(hive_exec.sh --put-skill) and re-run"}
    mine = set(own)
    groups, undeclared = {}, []
    for f in files:
        meth = declared_method(_base(f))
        if meth is None:
            undeclared.append(f)
        else:
            groups.setdefault(meth, []).append(f)
    grouped = [{"method": k, "n": len(v), "n_named_for_submission": sum(f in mine for f in v),
                "files": v} for k, v in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    lines = [f"{g['method']}: {g['n']} file(s), {g['n_named_for_submission']} named for the "
             f"submission, e.g. {_base(g['files'][0])}" for g in grouped]
    if len(grouped) > 1:
        return {"gate": "one_method", "status": "FAIL", "n": len(grouped), "groups": grouped,
                "undeclared": undeclared, "listed": lines,
                "detail": f"the files declare {len(grouped)} acquisition methods "
                          f"({', '.join(g['method'] for g in grouped)}). A plate is acquired one "
                          f"way; runs at another throughput or in another mode (DIA vs DDA) are "
                          f"another experiment and cannot be searched as one"}
    if undeclared:
        return {"gate": "one_method", "status": "WARN", "n": len(undeclared), "groups": grouped,
                "undeclared": undeclared, "listed": undeclared,
                "detail": f"{len(undeclared)} file name(s) declare no acquisition method (a "
                          f"re-injection is named <date>_<n>_rerun_...), so they could not be "
                          f"compared with the rest"
                          + (f" ({grouped[0]['method']})" if grouped else "")
                          + ". Step 2 (detect_acquisition.py) still reads DIA vs DDA from the "
                            "data; confirm the throughput with the operator"}
    return {"gate": "one_method", "status": "PASS", "n": len(grouped), "groups": grouped}


def sample_list_gate(a, own: list, entries: dict) -> dict:
    """`sample_list`: the CoreOmics record's samples (--summary, core_submission.py fetch's
    submission_summary.json) against the runs named for the submission -- every sample without a
    run and every run without a sample, named. Matched by core_submission.curated_rows, the rule
    `locate --files-from` applies to the same list, so the two can never disagree."""
    g = {"gate": "sample_list"}
    if not a.summary:
        return dict(g, status="WARN",
                    detail="not checked: no CoreOmics record given. STAN's list cannot show a run "
                           "STAN never recorded, so a sample with no run is invisible here; pass "
                           "--summary <submission_summary.json> (core_submission.py fetch), or "
                           "check it with core_submission.py locate --files-from")
    if a.include in ("rerun", "standards"):
        return dict(g, status="INFO", detail=f"not checked for --include {a.include}: that list "
                                             f"is not meant to hold every sample")
    if _cs is None:
        return dict(g, status="FAIL", detail=f"{CORE_UNAVAILABLE}, so the runs could not be "
                                             f"matched to the record's samples")
    try:
        s = _cs.load_json_file(a.summary, "--summary")
        ours = _cs.normalize_submission(a.submission)
    except _cs.Stop as e:
        return dict(g, status="FAIL", detail=e.payload["error"])
    except ValueError as e:
        return dict(g, status="FAIL", detail=str(e))
    try:
        theirs = _cs.normalize_submission(s.get("internal_id"))
    except ValueError:
        theirs = None
    if theirs != ours:
        return dict(g, status="FAIL",
                    detail=f"--summary {a.summary} is the record of "
                           f"{s.get('internal_id') or 'no submission number'}, not of {ours[1]}")
    samples = s.get("samples") or []
    if not samples:
        return dict(g, status="WARN", detail=f"not checked: the record of {ours[1]} lists no "
                                             f"samples")
    try:
        rows = _cs.curated_rows(samples, own, _cs.raw_root(), None, ours[1], "ask", entries)
    except _cs.Stop as e:
        return dict(g, status="FAIL", detail=e.payload["error"])
    missing = [{"unique_id": r["unique_id"], "sample_name": r["sample_name"]}
               for r in rows if r["status"] == "unmatched"]
    # STAN's blanks and standards (--include all) are no sample's, by design
    extra = [{"file": r["file"], "note": r["note"]} for r in rows if r["status"] == "unassigned"
             and ((entries.get(os.path.realpath(r["file"])) or {}).get("class") or "sample")
             == "sample"]
    several = sorted({r["unique_id"] for r in rows if r["status"] == "needs_choice"})
    g.update(n_samples=len(samples), n_runs=len(own), missing=missing, extra=extra,
             several_runs=several)
    if missing or extra:
        return dict(g, status="FAIL",
                    listed=[f"no run: {x['unique_id']} ({x['sample_name']})" for x in missing]
                    + [f"no sample: {x['file']}" for x in extra],
                    detail=f"the record of {ours[1]} lists {len(samples)} sample(s): "
                           f"{len(missing)} sample(s) have no run named for the submission, and "
                           f"{len(extra)} run(s) named for it belong to no sample on the record. "
                           f"A sample with no run: find its "
                           f"run (STAN may never have recorded it) and stage a curated list with "
                           f"core_submission.py locate --files-from <list> (--allow-partial there "
                           f"only if staff confirm it was not run). A run of no sample: ask the "
                           f"operator whose it is")
    if several:
        return dict(g, status="INFO",
                    detail=f"every sample has a run; {len(several)} sample(s) have several runs "
                           f"or share a name, and core_submission.py locate --files-from asks "
                           f"staff which run is which")
    return dict(g, status="PASS")


def fetch(a) -> int:
    use_http = bool(a.http or os.environ.get("STAN_HT_URL"))
    m = _manifest_over_http(a) if use_http else _manifest_over_cli(a)
    if isinstance(m, int):
        return m
    m.setdefault("source", "http" if use_http else "cli")

    files = m.get("files") or []
    gates, hard_fail = [], False

    # --- Repeats (check_report_runs: the one rule and the one wording). Brett (2026-10-01):
    # repeats are okay, but they are flagged. STAN listed one run 4 times and another twice (100
    # lines for 96 files) with every gate PASS (a staff plate, 2026-10-01). The same file is in
    # files.txt ONCE -- no engine can take it twice -- and flagged (`repeated_paths`, with the
    # count). Different files that share a run name (a rerun in another folder, one name on two
    # plates) are ALL kept and flagged (`repeated_names`); run_search.py stops before searching
    # them as given, because an engine would merge them into one run. Both are WARN, never a stop
    # here. STAN entries that name one file but DISAGREE (another well, run name, class or
    # verdict) are a hard fail: which sample that file is cannot be read from the manifest.
    listed = len(files)
    files, repeats = distinct_inputs(files)
    shared = repeated_names(files)
    by_path = {}
    for e in m.get("entries") or []:
        if isinstance(e, dict) and e.get("raw_path"):
            by_path.setdefault(os.path.realpath(str(e["raw_path"]).rstrip("/")), []).append(e)
    clash = []
    for v in by_path.values():
        if len(v) > 1 and any(x != v[0] for x in v[1:]):
            differ = sorted({k for x in v[1:] for k in set(x) | set(v[0]) if x.get(k) != v[0].get(k)})
            clash.append({"path": v[0]["raw_path"], "differ_in": differ,
                          "entries": [{k: x.get(k) for k in differ} for x in v]})
    if clash:
        # Not a repeat: a CONTRADICTION about which sample the file is. Brett's "repeats are
        # okay, but flagged" covers plain repeats only, so this stays a stop (2026-10-01).
        hard_fail = True
        gates.append({"gate": "repeated_paths", "status": "FAIL", "n": len(clash),
                      "detail": f"{len(clash)} file(s) are listed by STAN more than once with "
                                f"CONFLICTING entries. That is not a repeat but a contradiction "
                                f"about which sample the file is, so it stops here (repeats alone "
                                f"are only flagged). Ask the operator which entry is right",
                      "conflicts": clash,
                      "examples": [f"{c['path']}: " + " vs ".join(
                          ", ".join(f"{k}={e.get(k)}" for k in c["differ_in"]) for e in c["entries"])
                          for c in clash[:5]]})
    elif repeats:
        gates.append({"gate": "repeated_paths", "status": "WARN", "n": listed - len(files),
                      "detail": f"STAN listed {listed} runs for {len(files)} files. "
                                + repeats_text(repeats)[0] + ". STAN's counts include the repeats",
                      "examples": [f"{r['path']} x{r['times']}" for r in repeats[:5]]})
    else:
        gates.append({"gate": "repeated_paths", "status": "PASS", "n": 0})
    if shared:
        gates.append({"gate": "repeated_names", "status": "WARN", "n": len(shared),
                      "detail": repeats_text((), shared)[0] + ". run_search.py stops before "
                                "searching them as given (an engine would merge them into one "
                                "run): give one another name, or search them separately",
                      "examples": [f"{r['run_name']}: {', '.join(r['paths'])}" for r in shared[:5]]})
    else:
        gates.append({"gate": "repeated_names", "status": "PASS", "n": 0})

    # --- HARD: every file is NAMED for this submission. STAN assigns runs by tray: every run on
    # a tray the submission used, inside its injection window, joins, whoever's it is (STAN's
    # expand_submission_runs opens the window at the tray's earliest run nearby, not at the
    # submission's first). 2026-10-05: 136 files for a 96-sample plate, 41 of them another
    # project's earlier runs on the same tray, every gate PASS. Only the runs named for the
    # submission go into files.txt; the rest are listed, every one.
    plate = submission_plate(a.submission)
    own = [f for f in files if plate.match(_base(f))] if plate else []
    mine = set(own)
    others = [f for f in files if f not in mine]
    stan_entries = {k: v[0] for k, v in by_path.items()}
    named = named_gate(a, files, own, others, stan_entries, plate)
    gates.append(named)
    hard_fail = hard_fail or named["status"] == "FAIL"

    # --- HARD: one acquisition method across the list (100spd beside 60SPD DIA and 30SPD DDA,
    # 2026-10-05). Read from the names; DIA vs DDA is read again from the data at step 2.
    method = method_gate(files, own)
    gates.append(method)
    hard_fail = hard_fail or method["status"] == "FAIL"

    m["files"] = own
    m["n_files_listed"], m["n_files"] = listed, len(own)
    m["not_named_for_submission"] = others if plate else None     # None: names not checked
    m["repeated_paths"], m["repeated_names"] = repeats, shared

    # --- HARD: runs STAN knows about but has no raw path for. They are EXCLUDED from
    # `files`, so proceeding means searching a subset and then reporting success.
    missing = m.get("missing_paths") or []
    if missing:
        hard_fail = True
        gates.append({"gate": "missing_paths", "status": "FAIL", "n": len(missing),
                      "detail": f"{len(missing)} run(s) have no raw path and are excluded "
                                f"from the file list", "examples": missing[:5]})
    else:
        gates.append({"gate": "missing_paths", "status": "PASS", "n": 0})

    # --- HARD: nothing to search.
    if not files:
        hard_fail = True
        gates.append({"gate": "n_files", "status": "FAIL", "n": 0,
                      "detail": f"no files for submission {a.submission} with --include {a.include}"})
    elif not own:
        hard_fail = True
        gates.append({"gate": "n_files", "status": "FAIL", "n": 0,
                      "detail": (f"none of the {len(files)} file(s) STAN listed is named for "
                                 f"submission {a.submission}, so files.txt is empty: STAN's list "
                                 f"cannot be checked by name. Staff choose the files: write them "
                                 f"to a list and stage it with core_submission.py locate "
                                 f"--files-from <list> (references/core-submissions.md)")
                      if plate else "no file could be checked (see named_for_submission)"})
    else:
        gates.append({"gate": "n_files", "status": "PASS", "n": len(own)})

    # --- HARD: a path STAN resolved but the filesystem does not have. Cheap to check here;
    # otherwise a 120-file SLURM array dies partway through, hours in.
    unreadable = [f for f in own if not os.path.exists(f)]
    if unreadable:
        hard_fail = True
        gates.append({"gate": "paths_exist", "status": "FAIL", "n": len(unreadable),
                      "detail": f"{len(unreadable)} listed path(s) do not exist on this "
                                f"filesystem", "examples": unreadable[:5]})
    else:
        gates.append({"gate": "paths_exist", "status": "PASS", "n": len(own)})

    # --- WARN: shape checks. A submission fills one or two trays of 96.
    plates = m.get("plates") or []
    if len(plates) > MAX_PLAUSIBLE_PLATES:
        gates.append({"gate": "plates", "status": "WARN", "n": len(plates),
                      "detail": f"{len(plates)} plates ({', '.join(map(str, plates))}) — more "
                                f"than a submission usually spans. Confirm the extent with the "
                                f"operator; it is inferred from the acquisition counter."})
    else:
        gates.append({"gate": "plates", "status": "PASS", "n": len(plates), "plates": plates})

    counts = m.get("counts") or {}
    n_sample = counts.get("sample", len(files))
    if n_sample < MIN_PLAUSIBLE_SAMPLES:
        gates.append({"gate": "counts", "status": "WARN", "n": n_sample,
                      "detail": f"only {n_sample} customer samples — a mistyped submission "
                                f"number looks exactly like this. Confirm before searching."})
    else:
        gates.append({"gate": "counts", "status": "PASS", "n": n_sample, "counts": counts})

    # --- HARD with the CoreOmics record: every sample has a run, every run is a sample's. A run
    # STAN never recorded is in no STAN gate (2026-10-05: one of the 96 was missing).
    samples_gate = sample_list_gate(a, own, stan_entries)
    gates.append(samples_gate)
    hard_fail = hard_fail or samples_gate["status"] == "FAIL"

    n_rerun = m.get("n_needs_rerun") or 0
    if n_rerun:
        gates.append({"gate": "needs_rerun", "status": "INFO", "n": n_rerun,
                      "detail": f"{n_rerun} sample(s) flagged by STAN as wanting another "
                                f"injection. They ARE included here (--include samples). Use "
                                f"--include rerun to search only those."})

    m["gates"] = gates
    m["hard_fail"] = hard_fail

    outdir = a.out or "."
    os.makedirs(outdir, exist_ok=True)
    fl = os.path.join(outdir, "files.txt")
    mf = os.path.join(outdir, "ht_manifest.json")
    with open(fl, "w") as fh:
        fh.write("".join(f + "\n" for f in own))
    with open(mf, "w") as fh:
        json.dump(_scrubbed(m), fh, indent=2)

    why = (["repeats removed"] if listed != len(files) else []) + (
        [f"{len(others)} not named for submission {a.submission}" if plate else "names not checked"]
        if others else [])
    _say(f"submission {m.get('submission', a.submission)}  include={m.get('include', a.include)}")
    _say(f"  files      : {len(own)}" + (f" ({listed} listed by STAN; {'; '.join(why)})"
                                         if why else ""))
    _say(f"  plates     : {', '.join(map(str, plates)) or '(none reported)'}")
    _say(f"  counts     : {counts}")
    _say(f"  files.txt  : {fl}")
    _say(f"  manifest   : {mf}")
    for g in gates:
        if g["status"] != "PASS":
            _say(f"  [{g['status']}] {g['gate']}: {g.get('detail', '')}")
            # `listed`: every file or sample the gate is about; `examples`: the first five
            for ex in g.get("listed") or g.get("examples", [])[:5]:
                _say(f"        {ex}")
    if repeats or shared:
        _say(repeats_note(repeats, "ht_manifest", shared).rstrip("\n"), file=sys.stderr)
    if hard_fail:
        _say("\nHARD GATE FAILED — do not search. Resolve the above with the operator first.",
              file=sys.stderr)
        if plate and own and (others or samples_gate["status"] == "FAIL"):
            _say(f"  files.txt holds only the {len(own)} run(s) named for submission "
                 f"{a.submission}. Once the operator has said which runs are this submission's, "
                 f"stage that list with core_submission.py locate --files-from <list> (staff "
                 f"choices are recorded in locate.json). Never search STAN's list as given.",
                 file=sys.stderr)
        return 2
    return 0


def link(a) -> int:
    """Confirm the finished search is visible under its submission in FRAN.

    **There is no write-back, by design.** STAN v1.0.42 briefly had an `ht_searches` table
    and an `ht-record-search` command; v1.0.43 reverted both. The reasoning is worth
    keeping, because it is the reason this function only reads: a table in STAN would be a
    second copy of what FRAN already knows, and the copy is the one that goes stale --- it
    depends on the skill remembering to call it, and says nothing when that step is
    skipped. It would also mean re-implementing FRAN's authorization.

    So the link is implicit: the submission number is in the raw filenames, so FRAN's
    `/api/internal/submission/{id}` already answers "what was searched for 0793" with the
    search id, engine, organism and counts. Nothing needs recording. This just *checks*
    the loop closed, so a plate that silently failed to ingest is visible.
    """
    import urllib.error
    import urllib.request

    base = (a.fran or os.environ.get("FRAN_URL") or "https://fran.stan-proteomics.org").rstrip("/")
    url = f"{base}/api/internal/submission/{urllib_quote(a.submission)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    _, cookie = _credentials(a)
    if cookie:
        req.add_header("Cookie", cookie)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            payload = json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            _say(f"[ht_manifest] FRAN returned 404 for submission {a.submission}.\n"
                  f"  That endpoint is INTERNAL — 404 is also what it returns to a caller\n"
                  f"  who is not signed in, or to a lab user asking about another lab's\n"
                  f"  submission. So this is 'not visible to you', not proof of 'not ingested'.\n"
                  f"  Sign in at {base} (Microsoft Entra) or pass --cookie.", file=sys.stderr)
            return 4
        _say(f"[ht_manifest] FRAN returned HTTP {e.code} for {url}", file=sys.stderr)
        return 4
    except urllib.error.URLError as e:
        _say(f"[ht_manifest] cannot reach {base}: {e.reason}", file=sys.stderr)
        return 4

    body = payload.get("data", payload) if isinstance(payload, dict) else {}
    searches = (body or {}).get("searches") or []
    _say(f"submission {a.submission} — FRAN knows {len(searches)} search(es)")
    for s in searches:
        _say(f"  {str(s.get('search_name') or s.get('id'))[:44]:<46}"
              f"{str(s.get('search_engine') or '-'):<13}"
              f"prec={s.get('n_precursors_total')} pg={s.get('n_protein_groups_total')}")
    if not searches:
        _say("\n  No search under this submission yet. If fran_deposit.py reported\n"
              "  staged_pending_cron, that is expected — FRAN ingests on its next scan.\n"
              "  Re-run this later; it is not a failure.", file=sys.stderr)
        return 4
    return 0


def urllib_quote(s: str) -> str:
    import urllib.parse
    return urllib.parse.quote(str(s), safe="")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--stan", default=None, help=f"path to the stan binary (default {DEFAULT_STAN})")
    ap.add_argument("--token", default=None,
                    help=f"file holding the STAN Postgres credential (default {DEFAULT_TOKEN}; "
                         f"or set STAN_PG_TOKEN / PGPASSWORD)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="get + validate the file list for a submission")
    f.add_argument("submission", help="submission number, e.g. 0793 (793 also works)")
    f.add_argument("--include", default="samples",
                   choices=["samples", "rerun", "standards", "all"],
                   help="samples (default) excludes blanks and HeLa standards")
    f.add_argument("--out", default=None, help="directory for files.txt + ht_manifest.json")
    f.add_argument("--summary", default=None, metavar="FILE",
                   help="the submission's CoreOmics record (core_submission.py fetch's "
                        "submission_summary.json): every sample on it must have a run, and every "
                        "run a sample (gate sample_list). Without it that check is not made")
    f.add_argument("--http", default=None, metavar="BASE_URL",
                   help="ask the hosted dashboard instead of the local DB, e.g. "
                        "https://ucd.stan-proteomics.org (or set STAN_HT_URL). Use this when "
                        "you cannot read the Postgres token — which is every Core member "
                        "except its owner.")
    f.add_argument("--share-token-file", default=None, metavar="FILE",
                   help="file holding the per-submission share token from the HT tab (mode 600). "
                        "The only auth path that works headless.")
    f.add_argument("--share-token", default=None,
                   help="the token itself (or STAN_HT_SHARE_TOKEN). Avoid: it lands in the process "
                        "list and commands.log -- use --share-token-file")
    f.add_argument("--cookie-file", default=None, metavar="FILE",
                   help="file holding an Entra session cookie from a signed-in browser")
    f.add_argument("--cookie", default=None,
                   help="the cookie itself. Avoid, as for --share-token: use --cookie-file")
    f.set_defaults(func=fetch)

    r = sub.add_parser("link", help="check the finished search is visible under its "
                                    "submission in FRAN (read-only; there is no write-back)")
    r.add_argument("submission")
    r.add_argument("--fran", default=None, metavar="BASE_URL",
                   help="FRAN base URL (or FRAN_URL; default https://fran.stan-proteomics.org)")
    r.add_argument("--cookie-file", default=None, metavar="FILE",
                   help="file holding an Entra session cookie — the submission endpoint is "
                        "internal-only")
    r.add_argument("--cookie", default=None,
                   help="the cookie itself. Avoid: it lands in commands.log -- use --cookie-file")
    r.set_defaults(func=link)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
