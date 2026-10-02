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
import subprocess
import sys
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_report_runs import (distinct_inputs, repeated_names, repeats_note,  # noqa: E402
                               repeats_text)   # repeated inputs: the one rule

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
    m["files"] = files
    m["n_files_listed"], m["n_files"] = listed, len(files)
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
    else:
        gates.append({"gate": "n_files", "status": "PASS", "n": len(files)})

    # --- HARD: a path STAN resolved but the filesystem does not have. Cheap to check here;
    # otherwise a 120-file SLURM array dies partway through, hours in.
    unreadable = [f for f in files if not os.path.exists(f)]
    if unreadable:
        hard_fail = True
        gates.append({"gate": "paths_exist", "status": "FAIL", "n": len(unreadable),
                      "detail": f"{len(unreadable)} listed path(s) do not exist on this "
                                f"filesystem", "examples": unreadable[:5]})
    else:
        gates.append({"gate": "paths_exist", "status": "PASS", "n": len(files)})

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
        fh.write("".join(f + "\n" for f in files))
    with open(mf, "w") as fh:
        json.dump(_scrubbed(m), fh, indent=2)

    _say(f"submission {m.get('submission', a.submission)}  include={m.get('include', a.include)}")
    _say(f"  files      : {len(files)}" + (f" ({listed} listed by STAN; repeats removed)"
                                           if listed != len(files) else ""))
    _say(f"  plates     : {', '.join(map(str, plates)) or '(none reported)'}")
    _say(f"  counts     : {counts}")
    _say(f"  files.txt  : {fl}")
    _say(f"  manifest   : {mf}")
    for g in gates:
        if g["status"] != "PASS":
            _say(f"  [{g['status']}] {g['gate']}: {g.get('detail', '')}")
            for ex in g.get("examples", [])[:5]:
                _say(f"        {ex}")
    if repeats or shared:
        _say(repeats_note(repeats, "ht_manifest", shared).rstrip("\n"), file=sys.stderr)
    if hard_fail:
        _say("\nHARD GATE FAILED — do not search. Resolve the above with the operator first.",
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
